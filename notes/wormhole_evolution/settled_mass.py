#!/usr/bin/env python3
"""Settled black-hole mass at the Ellis threshold (#41): the Misner-Sharp mass on the apparent horizon, M_AH(t), followed past the
first trapped cell. Uses dn.py, checks.py and threshold.py of #39 and #40 unchanged.

A uniform grid cannot follow the horizon: along it e^{-f} grows exponentially, the horizon is a layer of width ~ r e^f in x-, and the
cell update loses it (#39's F_cut = 25). The grid here is a chain of levels. Level 0 is #39's uniform square grid (h, h). Level k >= 1 is a
box that starts on the column x+ = v_k, covers Nc cells of width h/2^k in x- centred on the horizon, and runs to the end of x+. Its data
are the row of the previous level at the box's lower edge and, on x+ = v_k, the previous level's column interpolated (cubic) to the
finer nodes. Characteristic evolution has no upstream dependence, so a level is exact in its domain of dependence. The next level
starts at the first column where e^{-f} on the horizon reaches F_{k+1} = 2 F_k, so that e^{-f} h_- is the same at the end of every level.

    python3 settled_mass.py planted   [--mutate=NAME]       planted exact solutions through the same chain (criteria PL)
    python3 settled_mass.py job --h=H --delta=D --out=DIR   one run: collapse at p*_h - delta, family S, D = 4
    python3 settled_mass.py analyze --raw=DIR [--mutate=NAME]

Units: a = 1, 8 pi G = 1, D = 4 (n = 2). M = r^(n-1)/2 (1 + 2 e^f r_+ r_-); on r_+ = 0 it is r^(n-1)/2.
"""
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dn           # noqa: E402

SQ2 = math.sqrt(2.0)
MUT = None


def mut(name):
    return MUT == name


# --------------------------------------------------------------------------------------------------------------------
# kernel: dn.row_step with separate spacings hp (x+) and hm (x-)
# --------------------------------------------------------------------------------------------------------------------
if dn.HAVE_NUMBA:
    import numba as nb
    from dn import _F

    @nb.njit(cache=True)
    def row_step2(old, new, dead_old, dead_new, hp, hm, par, eps, rcut, fcut):
        """new[:, 0] and dead_new[0] are set by the caller; fills new[:, 1:] from old and new[:, :-1]. Same arithmetic as dn.row_step."""
        nf, n1 = old.shape
        hh = hp * hm
        ncorr = int(par[4])
        c = np.empty(nf)
        up = np.empty(nf)
        um = np.empty(nf)
        nh = np.empty(nf)
        g = np.empty(nf)
        for i in range(n1 - 1):
            if dead_old[i] or dead_old[i + 1] or dead_new[i]:
                dead_new[i + 1] = 1
                for q in range(nf):
                    new[q, i + 1] = old[q, i + 1]
                continue
            for q in range(nf):
                S = old[q, i]
                E = old[q, i + 1]
                W = new[q, i]
                c[q] = 0.5 * (E + W)
                up[q] = (E - S) / hp
                um[q] = (W - S) / hm
            _F(c, up, um, par, eps, g)
            for q in range(nf):
                nh[q] = old[q, i + 1] + new[q, i] - old[q, i] + hh * g[q]
            for _ in range(ncorr):
                for q in range(nf):
                    S = old[q, i]
                    E = old[q, i + 1]
                    W = new[q, i]
                    c[q] = 0.25 * (S + E + W + nh[q])
                    up[q] = (E + nh[q] - S - W) / (2 * hp)
                    um[q] = (W + nh[q] - S - E) / (2 * hm)
                _F(c, up, um, par, eps, g)
                for q in range(nf):
                    nh[q] = old[q, i + 1] + new[q, i] - old[q, i] + hh * g[q]
            bad = False
            for q in range(nf):
                if not (nh[q] == nh[q]) or abs(nh[q]) > 1e300:
                    bad = True
            if bad or nh[0] < rcut or math.exp(-nh[1]) > fcut:
                dead_new[i + 1] = 1
                for q in range(nf):
                    new[q, i + 1] = old[q, i + 1]
            else:
                dead_new[i + 1] = 0
                for q in range(nf):
                    new[q, i + 1] = nh[q]
else:
    raise SystemExit("settled_mass.py needs numba (the chain has no numpy kernel)")


class Level:
    """A box of nodes G[q, i, j], x+ = xp0 + i hp, x- = xm0 + j hm, with dead flags."""

    def __init__(self, k, G, dead, hp, hm, xp0, xm0, i0=0):
        self.k, self.G, self.dead, self.hp, self.hm, self.xp0, self.xm0, self.i0 = k, G, dead, hp, hm, xp0, xm0, i0
        self.nf, self.Np, self.Nm = G.shape[0], G.shape[1] - 1, G.shape[2] - 1


def evolve_box(Sp, Sm, dSp, dSm, hp, hm, n, eps, rcut=0.05, fcut=1e300, ncorr=1, detect=None):
    """Characteristic evolution on the rectangle from the lower ray Sp[q, i] (x- = x0m) and the left ray Sm[q, j] (x+ = x0p).
    detect = tol: also return the first trapped cell (r_+ and r_- both below -tol; smallest row, then smallest column) as (j, i), as
    dn.events_from_grid, or None."""
    eps = np.asarray(eps, float)
    par = dn.kernel_params(n, eps, ncorr)
    nf, np1 = Sp.shape
    nm1 = Sm.shape[1]
    G = np.empty((nf, np1, nm1))
    dead = np.zeros((np1, nm1), np.uint8)
    G[:, :, 0] = Sp
    dead[:, 0] = dSp
    old = np.ascontiguousarray(Sp.copy())
    dold = np.ascontiguousarray(dSp.copy())
    new = np.empty_like(old)
    dnew = np.zeros_like(dold)
    hm_k = hp if mut("hm_swap") else hm
    event = None
    for j in range(nm1 - 1):
        new[:, 0] = Sm[:, j + 1]
        dnew[0] = dSm[j + 1]
        row_step2(old, new, dold, dnew, hp, hm_k, par, eps, rcut, fcut)
        if detect is not None and event is None:
            S_, E_, W_, N_ = old[0, :-1], old[0, 1:], new[0, :-1], new[0, 1:]
            rp_ = (E_ + N_ - S_ - W_) / (2 * hp)
            rm_ = (W_ + N_ - S_ - E_) / (2 * hm_k)
            ok_ = ~(dold[:-1].astype(bool) | dold[1:].astype(bool) | dnew[:-1].astype(bool) | dnew[1:].astype(bool))
            tr_ = ok_ & (rp_ < -detect) & (rm_ < -detect)
            if tr_.any():
                event = (j, int(np.argmax(tr_)))
        G[:, :, j + 1] = new
        dead[:, j + 1] = dnew
        old, new = new, old
        dold, dnew = dnew, dold
    return (G, dead) if detect is None else (G, dead, event)


def zoom(prev, i0, jc, Nc):
    """Level k+1 from level k: starts at node column i0 of prev, Nc cells of width prev.hm/2, lower edge at the prev row nearest
    jc - Nc/4 (jc fractional prev row of the horizon). Lower ray: prev's row. Left ray: prev's column, cubic midpoints at odd nodes."""
    j0 = int(round(jc - Nc / 4.0))
    if j0 < 1 or j0 + Nc // 2 + 2 > prev.Nm:
        raise RuntimeError(f"horizon too close to the edge of level {prev.k}: row {jc:.1f} of {prev.Nm}, box needs [{j0 - 1}, {j0 + Nc // 2 + 2}]")
    Sp = prev.G[:, i0:, j0].copy()
    dSp = prev.dead[i0:, j0].copy()
    col = prev.G[:, i0, :]
    dcol = prev.dead[i0, :]
    nm = Nc + 1
    Sm = np.empty((prev.nf, nm))
    dSm = np.zeros(nm, np.uint8)
    ev = np.arange(0, nm, 2)
    od = np.arange(1, nm, 2)
    Sm[:, ev] = col[:, j0 + ev // 2]
    dSm[ev] = dcol[j0 + ev // 2]
    js = j0 + (od - 1) // 2
    t = 0.25 if mut("interp_shift") else 0.5                      # position between nodes js and js+1 (mutation: wrong)
    w = np.array([-t * (1 - t) * (2 - t) / 6, (1 + t) * (1 - t) * (2 - t) / 2, (1 + t) * t * (2 - t) / 2, -(1 + t) * t * (1 - t) / 6])
    Sm[:, od] = w[0] * col[:, js - 1] + w[1] * col[:, js] + w[2] * col[:, js + 1] + w[3] * col[:, js + 2]
    dSm[od] = dcol[js - 1] | dcol[js] | dcol[js + 1] | dcol[js + 2]
    return Sp, Sm, dSp, dSm, j0, prev.hm / 2


def cell_window(L, c, jlo, jhi):
    """Cell-centre r, f, r_+, r_-, scalar derivatives and the alive mask of cell column c, rows jlo..jhi (cell rows)."""
    G, dead = L.G, L.dead
    S, E = G[:, c, jlo:jhi + 1], G[:, c + 1, jlo:jhi + 1]
    W, Nn = G[:, c, jlo + 1:jhi + 2], G[:, c + 1, jlo + 1:jhi + 2]
    cc = 0.25 * (S + E + W + Nn)
    up = (E + Nn - S - W) / (2 * L.hp)
    um = (W + Nn - S - E) / (2 * L.hm)
    al = ~(dead[c, jlo:jhi + 1].astype(bool) | dead[c + 1, jlo:jhi + 1].astype(bool) | dead[c, jlo + 1:jhi + 2].astype(bool) | dead[c + 1, jlo + 1:jhi + 2].astype(bool))
    return cc, up, um, al


def track(L, n, c_start, jguess, stop_F=None, tol=1e-2, win=60, require_trapped=True, c_end=None, stop_xp=None):
    """Follow the zero of r_+ (+ to - along x-, with r_- < -tol) down the cell columns of L from c_start. jguess: fractional row of the
    crossing at c_start (cell centres are at j + 1/2). Stops at the first column whose e^{-f} on the horizon is >= stop_F; that column is
    the last one returned and its index is `stop`. Per column: x+ (cell centre), x-, r, F = e^{-f}, r_-, the scalar derivatives phi_+,
    phi_-, the interpolated cell Misner-Sharp mass (identity check), the row and its distance to the box edges."""
    keys = ("c", "xp", "xm", "r", "F", "rm", "phip", "phim", "Mcell", "row", "edge")
    out = {k: [] for k in keys}
    stop = None
    jg = jguess
    c_end = L.Np if c_end is None else c_end
    ex = n if mut("mass_exponent") else n - 1
    sgn = -1.0 if mut("finder_sign") else 1.0
    exhausted = True
    for c in range(c_start, c_end):
        jlo = max(0, int(jg) - win)
        jhi = min(L.Nm - 2, int(jg) + win)
        cc, up, um, al = cell_window(L, c, jlo, jhi)
        rp = up[0]
        idx = np.where(al[:-1] & al[1:] & (sgn * rp[:-1] > 0) & (sgn * rp[1:] <= 0))[0]
        cand = []
        for jj in idx:
            fr = rp[jj] / (rp[jj] - rp[jj + 1])
            rmi = um[0][jj] + fr * (um[0][jj + 1] - um[0][jj])
            if require_trapped and not rmi < -tol:
                continue
            cand.append((abs(jlo + jj + 0.5 + fr - jg), jj, fr))
        if not cand:
            exhausted = False
            break
        _, jj, fr = min(cand)

        def lin(a):
            return a[jj] + fr * (a[jj + 1] - a[jj])
        r, f = lin(cc[0]), lin(cc[1])
        row = jlo + jj + 0.5 + fr
        out["c"].append(c)
        out["xp"].append(L.xp0 + (c if mut("col_shift") else c + 0.5) * L.hp)
        out["xm"].append(L.xm0 + row * L.hm)
        out["r"].append(r)
        out["F"].append(math.exp(-f))
        out["rm"].append(lin(um[0]))
        out["phip"].append(lin(up[2]) if L.nf > 2 else 0.0)
        out["phim"].append(lin(um[2]) if L.nf > 2 else 0.0)
        out["Mcell"].append(lin(cc[0] ** ex * (1 + 2 * np.exp(cc[1]) * up[0] * um[0]) / 2))
        out["row"].append(row)
        out["edge"].append(min(row, L.Nm - row))
        jg = row
        if (stop_F is not None and math.exp(-f) >= stop_F) or (stop_xp is not None and out["xp"][-1] >= stop_xp):
            stop = c
            exhausted = False
            break
    out = {k: np.array(v) for k, v in out.items()}
    out["M"] = out["r"] ** ex / 2
    return out, stop, exhausted


# --------------------------------------------------------------------------------------------------------------------
# the chain
# --------------------------------------------------------------------------------------------------------------------
def chain(n, eps, hp, Sp, Sm, xp0, xm0, nm0, c0, jg0, F1, Fmax, Nc, forced=None, tol=1e-2, require_trapped=True, rcut=0.05, win=60, log=print):
    """Level 0 over len(Sp)-1 columns and nm0 rows from the rays Sp, Sm (corner at (xp0, xm0)), horizon followed from column c0 near row jg0.
    Level k >= 1 starts where e^{-f} on the horizon reaches F1 2^(k-1) (or at forced[k-1], an x+ value), has Nc cells of width hp/2^k, and
    the chain ends when e^{-f} reaches Fmax (or after len(forced) levels). c0 = None: start at the first trapped cell. Returns levels, per-level records, the joined series, status, event."""
    t0 = time.time()
    rcut = rcut * 1.5 if mut("rcut_off") else rcut
    G, dead, event = evolve_box(Sp, Sm[:, :nm0 + 1], np.zeros(Sp.shape[1], np.uint8), np.zeros(nm0 + 1, np.uint8), hp, hp, n, eps, rcut, detect=tol)
    L = Level(0, G, dead, hp, hp, xp0, xm0)
    levels, recs = [L], []
    if c0 is None:                                               # start at the first trapped cell (#40's event)
        if event is None:
            return dict(levels=levels, recs=[], ser={}, status="no trapped cell in level 0", secs=time.time() - t0, event=None)
        c0, jg0 = event[1], event[0] + 0.5
    status = "ok"
    c_start, jg = c0, jg0
    k = 0
    while True:
        if forced is not None:
            stop_F, stop_xp = None, (forced[k] if k < len(forced) else None)
        else:
            stop_F, stop_xp = F1 * 2.0 ** k, None
        last = (forced is not None and k >= len(forced)) or (forced is None and stop_F >= Fmax * (1 - 1e-12))
        if last:
            stop_F = Fmax if forced is None else None
            stop_xp = None
        rec, stop, exhausted = track(L, n, c_start, jg, stop_F=stop_F, tol=tol, win=win, require_trapped=require_trapped, stop_xp=stop_xp)
        rec["level"] = np.full(len(rec["r"]), k)
        if len(rec["r"]) == 0:
            status = f"level {k}: no horizon at the first column"
            recs.append(rec)
            break
        if stop is None:
            status = ("end" if (exhausted and last and forced is not None) else
                      f"level {k}: {'columns exhausted' if exhausted else 'horizon lost'} at x+ = {rec['xp'][-1]:.4f}, F = {rec['F'][-1]:.4g}")
            recs.append(rec)
            break
        cut = {kk: v[:-1] for kk, v in rec.items()}                  # the stop column belongs to the next level
        recs.append(cut if not last else cut)
        if last:
            status = "limit"
            break
        Spz, Smz, dSp, dSm, j0, hmz = zoom(L, stop, rec["row"][-1], Nc)
        Gz, deadz = evolve_box(Spz, Smz, dSp, dSm, hp, hmz, n, eps, rcut)
        L2 = Level(k + 1, Gz, deadz, hp, hmz, L.xp0 + (L.i0 + stop) * hp, L.xm0 + j0 * L.hm)
        L2.i0 = 0
        log(f"  level {k + 1}: x+ from {L2.xp0:.4f}, x- [{L2.xm0:.5f}, {L2.xm0 + Nc * hmz:.5f}], hm {hmz:.3e}, {Gz.shape[1] - 1} x {Nc} cells, "
            f"F at start {rec['F'][-1]:.4g} ({time.time() - t0:.1f} s)")
        if L.k >= 1:
            L.G, L.dead = None, None                              # a finished level is not needed again (level 0 is kept for PL3)
        levels.append(L2)
        jg = (rec["xm"][-1] - L2.xm0) / hmz
        L, c_start = L2, 0
        k += 1
    ser = {kk: np.concatenate([r[kk] for r in recs if len(r["r"])]) if any(len(r["r"]) for r in recs) else np.array([]) for kk in recs[0]}
    return dict(levels=levels, recs=recs, ser=ser, status=status, secs=time.time() - t0, event=event)


# --------------------------------------------------------------------------------------------------------------------
# planted data
# --------------------------------------------------------------------------------------------------------------------
def lambertw0(x):
    from scipy.special import lambertw
    return np.real(lambertw(x, 0))


def schw_exact(lam, kp, xh, xp, xm, xs=3.0):
    """Exact Schwarzschild, D = 4, horizon radius lam (M = lam/2), in the gauge x- = xh + U (Kruskal U, scaled), x+ with V = lam V0' e^{kp x+}:
    the horizon is x- = xh and e^{-f} on it is e^{kp (x+ - xs)}. Returns r, f at (xp, xm) arrays. Unit solution: (r'-1) e^{r'-1} = -U' V',
    e^{-f'} = (2/r') e^{1-r'}; scaling x -> lam x' leaves f unchanged and scales r."""
    xp, xm = np.broadcast_arrays(np.asarray(xp, float), np.asarray(xm, float))
    Up = (xm - xh) / lam
    lnVp = -math.log(2.0 * kp * lam) + kp * (xp - xs)             # ln V'; e^{-f} on the horizon is e^{kp (x+ - xs)}
    Zs = -Up * np.exp(lnVp)                                        # (r'-1) e^{r'-1} = Zs (exterior Zs > 0)
    with np.errstate(over="ignore", invalid="ignore"):
        w = np.where(Zs < 1e300, lambertw0(np.minimum(Zs, 1e300)), 0.0)
    # large Zs: W(z) by the asymptotic series is not needed below 1e300; r' = 1 + W
    rp = 1.0 + w
    lnemf = math.log(2.0) - np.log(rp) + 1.0 - rp + math.log(kp * lam) + lnVp   # ln e^{-f~} = ln(2/r') + 1 - r' + ln(kp lam V')
    return lam * rp, -lnemf


def schw_rays(lam, kp, xh, hp, Np, nm0, xs=3.0):
    xs_ = np.arange(Np + 1) * hp
    ms = np.arange(nm0 + 1) * hp
    Sp = np.zeros((2, Np + 1))
    Sm = np.zeros((2, nm0 + 1))
    Sp[0], Sp[1] = schw_exact(lam, kp, xh, xs_, 0.0 * xs_, xs)
    Sm[0], Sm[1] = schw_exact(lam, kp, xh, 0.0 * ms, ms, xs)
    return Sp, Sm


def roberts_planted(h, Nc, forced=(3.0, 3.08), p_end=3.12, a=0.75, p0=1.0, m0=1.0):
    """The exact ghost solution r^2 = a (p^2 + m^2) - p m (D = 4, f = 0) through the chain with forced zoom levels. Its contour r_+ = 0 is
    m = 2 a p, a moving line with r_- > 0 there (not trapped), so the finder runs with require_trapped = False. Returns the series and errors."""
    Np = int(round((p_end - p0) / h))
    nm0 = int(round((2 * a * p_end + 0.3 - m0 + Nc * h / 4) / h)) + 8
    d = dn.roberts_rays(max(Np, nm0), h, a, p0, m0)
    out = chain(2, d.eps, h, d.Sp[:, :Np + 1], d.Sm, p0, m0, nm0, 0, (2 * a * p0 - m0) / h, 1.0, 1.0, Nc, forced=list(forced),
                require_trapped=False, win=60, log=lambda *x: None)
    s = out["ser"]
    if len(s) == 0 or len(s["r"]) == 0:
        return out, dict(e_r=np.array([]), e_m=np.array([]), e_M=np.array([]), r_ex=np.array([]))
    r_ex = np.sqrt(a * (s["xp"] ** 2 + (2 * a * s["xp"]) ** 2) - s["xp"] * 2 * a * s["xp"])
    return out, dict(e_r=np.abs(s["r"] - r_ex), e_m=np.abs(s["xm"] - 2 * a * s["xp"]), e_M=np.abs(s["M"] - r_ex / 2), r_ex=r_ex)


def schw_planted(h, F1, K, Nc, lam=0.6, kp=2.0, xh=5.0, xs=3.0):
    """Exact Schwarzschild in a gauge where e^{-f} on the horizon grows as e^{kp (x+ - xs)} (the stretching of the real run), through the
    same chain. Returns the series and the errors of M and of the horizon position."""
    Fmax = F1 * 2.0 ** K
    Xf = xs + math.log(Fmax) / kp + 1.0
    Np = int(round(Xf / h))
    nm0 = int(round(xh / h)) + Nc // 4 + 6
    Sp, Sm = schw_rays(lam, kp, xh, h, Np, nm0, xs)
    out = chain(2, [], h, Sp, Sm, 0.0, 0.0, nm0, 0, xh / h, F1, Fmax, Nc, tol=0.0, log=lambda *x: None)
    s = out["ser"]
    return out, dict(e_M=np.abs(s["M"] - lam / 2) if len(s["M"]) else np.array([np.inf]), e_x=np.abs(s["xm"] - xh) if len(s["M"]) else np.array([np.inf]))


# --------------------------------------------------------------------------------------------------------------------
# criteria and constants (declared in negative-mass-raw I41-rpi5/README.md before the runs of record)
# --------------------------------------------------------------------------------------------------------------------
import threshold as th    # noqa: E402

CRIT = dict(
    n=2, hs=(0.02, 0.01, 0.005, 0.0025), deltas=(1e-6, 1e-8), X=28.0,  # D = 4, family S, the three h of #40 and a fourth; X as #40's rectangle
    F1=1.5, K=20, W1=2.4, rcut=0.05, tol=1e-2, dxp=12.0, dm=3.0,      # chain: F_max = F1 2^K = 1572864; level 1 is W1 wide: Nc = 2 W1 / h cells
    bisect_tol=1e-14,
)
I40 = os.path.expanduser("~/Source/negative-mass-raw/I40-dh610/runs")


def nc_of(h):
    """Cells per box: level 1 is W1 = 2.4 wide (the horizon moves 0.24 in x- inside it) at width h/2, so Nc = 2 W1 / h (240, 480, 960, 1920)."""
    return int(round(2.0 * CRIT["W1"] / h / 4.0)) * 4


def i40_ref(h):
    """#40's p*_h and its formation values for D = 4, S (dh610-fedora), or None."""
    p = f"{I40}/D4_S_h{h:g}"
    if not os.path.exists(p + "/summary.json"):
        return None
    s = json.load(open(p + "/summary.json"))
    st = json.load(open(p + "/stall.json"))["-1"]
    return dict(ps=0.5 * (s["lo"] + s["hi"]) + (1e-10 if mut("i40_shift") else 0.0), lo=s["lo"], hi=s["hi"], event=st["event"], T=st["T"], M_hor=st["M_hor"], M_cell=st["M_cell"], r_cell=st["r_cell"])


def pstar_cmd(h, out):
    """Bisection of the fate of family S to hi - lo <= 1e-14 on this machine, around #40's p*_h (or the h^2 extrapolation of it)."""
    os.makedirs(out, exist_ok=True)
    ref = i40_ref(h)
    c = ref["ps"] if ref else 2.9711800e-06 * (h / 0.005) ** 2
    sink_f = open(os.path.join(out, f"pstar_h{h:g}.jsonl"), "w")

    def run(p):
        r = th.point(CRIT["n"], "S", h, CRIT["X"], p)
        sink_f.write(json.dumps(r) + "\n")
        sink_f.flush()
        return r
    t0 = time.time()
    for w in ((2e-10, 2e-9, 2e-8, 2e-7) if ref else (1e-8, 1e-7)):
        try:
            lo, hi, stopped = th.bisect(run, c - w, c + w, CRIT["bisect_tol"], lambda r: None)
            break
        except RuntimeError as e:
            print("bracket", w, "failed:", e, flush=True)
    else:
        raise SystemExit("no bracket")
    res = dict(h=h, lo=lo, hi=hi, ps=0.5 * (lo + hi), stopped=stopped, i40_ps=ref["ps"] if ref else None,
               diff_i40=(0.5 * (lo + hi) - ref["ps"]) if ref else None, bracket=w, secs=round(time.time() - t0, 1), runs=sum(1 for _ in open(sink_f.name)))
    json.dump(res, open(os.path.join(out, f"pstar_h{h:g}.json"), "w"), indent=1)
    print(json.dumps(res), flush=True)



# --------------------------------------------------------------------------------------------------------------------
# the real runs
# --------------------------------------------------------------------------------------------------------------------
def run_real(h, delta, ps, F1=None, K=None, Nc=None, log=lambda *a: None, need_up=True):
    """Family S at p*_h - delta (collapse side): #40's own point() for the fate and the first trapped cell, then the chain from there."""
    F1 = CRIT["F1"] if F1 is None else F1
    K = CRIT["K"] if K is None else K
    Nc = nc_of(h) if Nc is None else Nc
    A = ps - delta
    pulses, sym, eps = th.family("S", A)
    pt = th.point(CRIT["n"], "S", h, CRIT["X"], A)
    up = th.point(CRIT["n"], "S", h, CRIT["X"], ps + delta) if need_up else dict(fate=1)
    if pt["fate"] != -1 or up["fate"] != 1:
        raise RuntimeError(f"fates at p*_h -/+ delta: {pt['fate']}, {up['fate']} (wanted -1, +1)")
    j_ev, i_ev = pt["j"], pt["i"]
    Np = i_ev + int(round(CRIT["dxp"] / h))
    nm0 = j_ev + int(round((CRIT["dm"] + Nc * h / 4) / h))
    d = dn.ellis_rays(CRIT["n"], max(Np, nm0), h, pulses=pulses, symmetric=sym, eps=eps)
    out = chain(CRIT["n"], eps, h, d.Sp[:, :Np + 1], d.Sm, 0.0, 0.0, nm0, None, None, F1, F1 * 2.0 ** K, Nc, tol=CRIT["tol"], rcut=CRIT["rcut"], log=log)
    out.update(point=pt, up=up, A=A, Np=Np, nm0=nm0, F1=F1, K=K, Nc=Nc, h=h, delta=delta, ps=ps, data=d)
    return out


def t_of(xp, xm):
    return (xp + xm) if mut("time_unit") else (xp + xm) / SQ2


def job_cmd(h, delta, out, pstar_dir):
    import resource
    ps = json.load(open(os.path.join(pstar_dir, f"pstar_h{h:g}.json")))["ps"]
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    r = run_real(h, delta, ps, log=lambda m: print(m, flush=True))
    s = r["ser"]
    np.savez(os.path.join(out, "series.npz"), **s)
    summ = dict(h=h, delta=delta, ps=ps, A=r["A"], status=r["status"], event=r["event"], point=r["point"], up=r["up"], Np=r["Np"], nm0=r["nm0"],
                F1=r["F1"], K=r["K"], Nc=r["Nc"], levels=len(r["levels"]), cols=int(len(s["r"])), xp0=float(s["xp"][0]), xp1=float(s["xp"][-1]),
                F_end=float(s["F"][-1]), M_start=float(s["M"][0]), M_end=float(s["M"][-1]), min_edge=float(s["edge"].min()), secs=round(time.time() - t0, 1),
                maxrss_MB=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 1))
    json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print(json.dumps(summ), flush=True)


def load_run(raw, h, delta):
    p = os.path.join(raw, "runs", f"h{h:g}_d{delta:g}")
    s = dict(np.load(os.path.join(p, "series.npz")))
    return s, json.load(open(os.path.join(p, "summary.json")))


def get_run(raw, h, delta):
    """The stored series, or, under a mutation, a fresh run from the stored p*_h."""
    if MUT is None:
        return load_run(raw, h, delta)
    ps = json.load(open(os.path.join(raw, "runs", f"pstar_h{h:g}.json")))["ps"]
    K = 0 if mut("no_zoom") else (12 if mut("short_limit") else None)
    ps = ps + 2e-8 if mut("pstar_off") else ps
    r = run_real(h, 3e-2 if (mut("delta_swap") and delta == 1e-6) else delta, ps, K=K)
    s = r["ser"]
    return s, dict(status=r["status"], event=r["event"], point=r["point"], up=r["up"], Nc=r["Nc"], F1=r["F1"], K=r["K"], h=h, delta=delta, ps=ps)


# --------------------------------------------------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------------------------------------------------
results = []
BAND_P = (1.75, 2.25)           # planted: #39's band for the observed order
BAND_R = (1.5, 2.5)             # real series: wider, set from the development runs (series noise of order 1e-5)


def check(cid, text, ok, detail=""):
    ok = None if ok is None else bool(ok)
    results.append((cid, ok))
    print(f"[{'SKIP' if ok is None else ('pass' if ok else 'FAIL')}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def mx(a):
    return float(np.max(a)) if len(a) else float("inf")


def order(e12, e23):
    try:
        return math.log2(e12 / e23) if e12 > 0 and e23 > 0 else float("inf")
    except (ValueError, ZeroDivisionError, OverflowError):
        return float("inf")


def flux_integral(s, i0, i1):
    """Integral of dM/dx+ = (r^n/n) e^f r_- Sum(eps phi_+^2) = -(r^2/2)(1/F)|r_-| phi_+^2 (one ghost, n = 2) along the horizon, columns i0..i1."""
    sg = -1.0 if mut("flux_sign") else 1.0
    f = sg * 0.5 * s["r"] ** 2 / s["F"] * s["rm"] * s["phip"] ** 2
    x = s["xp"][i0:i1 + 1]
    return float(np.sum(0.5 * (f[i0:i1] + f[i0 + 1:i1 + 1]) * np.diff(x)))


def planted_checks(hs=None, log=print):
    hs = CRIT["hs"] if hs is None else hs
    F1, K = CRIT["F1"], CRIT["K"]
    if mut("no_zoom"):
        K = 0
    eM, eX, st, Fend = [], [], [], []
    for h in hs:
        out, e = schw_planted(h, F1, K, nc_of(h))
        s = out["ser"]
        eM.append(mx(e["e_M"]))
        eX.append(mx(e["e_x"]))
        st.append(out["status"])
        Fend.append(float(s["F"][-1]) if len(s["F"]) else 0.0)
        log(f"    P1 h = {h:g}: status {out['status']}, columns {len(s['M'])}, F to {Fend[-1]:.4g}, max |M - 0.3| {eM[-1]:.3e}, max |x- - x_h| {eX[-1]:.3e}")
    o = [order(eM[i], eM[i + 1]) for i in range(len(hs) - 1)]
    check("PL1", "exact Schwarzschild (M = 0.3) in the stretching gauge, through the chain: M recovered, second order, limit reached",
          eM[2] <= 1e-4 and all(BAND_P[0] <= x <= BAND_P[1] for x in o) and all(x == "limit" for x in st) and all(f >= 0.5 * F1 * 2.0 ** CRIT["K"] for f in Fend),
          f"max |M - 0.3| at h = {hs[2]:g}: {eM[2]:.2e} (bar 1e-4); orders {[round(x, 3) for x in o]} in {BAND_P}; status {st}; F reached {[round(f) for f in Fend]} of {F1 * 2.0 ** CRIT['K']:g}")
    eR = []
    for h in hs:
        out, e = roberts_planted(h, nc_of(h))
        s = out["ser"]
        eR.append(mx(e["e_r"]))
        log(f"    P2 h = {h:g}: columns {len(s['r'])} (levels {[int((s['level'] == k).sum()) for k in range(3)] if len(s['r']) else []}), max |r - r_exact| on r_+ = 0: {eR[-1]:.3e}, max |x- - 2 a x+| {mx(e['e_m']):.3e}, min edge {(s['edge'].min() if len(s['r']) else float('nan')):.0f} cells")
    o = [order(eR[i], eR[i + 1]) for i in range(len(hs) - 1)]
    check("PL2", "exact Roberts ghost solution, moving contour r_+ = 0, two forced zoom levels: r on the contour, second order",
          eR[2] <= 5e-5 and all(BAND_P[0] <= x <= BAND_P[1] for x in o), f"max error at h = {hs[2]:g}: {eR[2]:.2e} (bar 5e-5); orders {[round(x, 3) for x in o]} in {BAND_P}")
    return eM, eR


def kernel_identity(raw, h=0.02, delta=1e-6):
    """Level 0 of the chain against dn.evolve (numba row sweep) on the same data: bit-identical nodes and dead flags, and the same first trapped cell."""
    ps = json.load(open(os.path.join(raw, "runs", f"pstar_h{h:g}.json")))["ps"]
    r = run_real(h, delta, ps)
    d, L0 = r["data"], r["levels"][0]
    res = dn.evolve(d, store=True, stop=False, rcut=CRIT["rcut"], fcut=1e300, tol=CRIT["tol"])
    same_G = np.array_equal(res.G[:, :L0.Np + 1, :L0.Nm + 1], L0.G)
    same_d = np.array_equal(res.dead[:L0.Np + 1, :L0.Nm + 1], L0.dead)
    ev = res.events.get("collapse")
    check("PL3", "level 0 of the chain equals dn.evolve bit for bit (nodes, dead flags) and finds the same first trapped cell",
          same_G and same_d and ev == r["event"], f"h = {h:g}, delta = {delta:g}: nodes equal {same_G}, dead equal {same_d}, events {ev} and {r['event']}")


def band_orders(S, hs, windows, delta_tag):
    """Max-norm differences of M(x+) between consecutive h (cubic splines on common x+), per window of x+ - x+_f."""
    from scipy.interpolate import CubicSpline
    cs = {h: CubicSpline(S[h]["xp"], S[h]["M"]) for h in hs}
    vf = S[hs[-1]]["xp"][0]
    lo_all = max(S[h]["xp"][0] for h in hs) + max(hs)
    hi_all = min(S[h]["xp"][-1] for h in hs) - max(hs)
    rows = []
    for a, b in windows:
        v = np.linspace(max(vf + a, lo_all), min(vf + b, hi_all), 2000)
        E = [float(np.max(np.abs(cs[hs[i]](v) - cs[hs[i + 1]](v)))) for i in range(len(hs) - 1)]
        rows.append((a, b, E, [order(E[i], E[i + 1]) for i in range(len(E) - 1)]))
    return rows



def fate_aniso(h, hm, A, Np, nrows, tol=1e-2, rexp=2.0):
    """Fate of family S at amplitude A on the uniform anisotropic grid (h, hm), streaming (no stored grid): collapse = first trapped cell,
    expansion = the cell on the diagonal x+ = x- anti-trapped with r >= rexp (as #40's events). Returns dict(fate, j, i)."""
    pulses, sym, eps = th.family("S", A)
    d1 = dn.ellis_rays(CRIT["n"], Np, h, pulses=pulses, symmetric=sym, eps=eps)
    d2 = dn.ellis_rays(CRIT["n"], nrows, hm, pulses=pulses, symmetric=sym, eps=eps)
    par = dn.kernel_params(CRIT["n"], eps, 1)
    old, dold = d1.Sp.copy(), np.zeros(Np + 1, np.uint8)
    new, dnew = np.empty_like(old), np.zeros(Np + 1, np.uint8)
    for j in range(nrows):
        new[:, 0] = d2.Sm[:, j + 1]
        dnew[0] = 0
        row_step2(old, new, dold, dnew, h, hm, par, eps, CRIT["rcut"], 1e300)
        S_, E_, W_, N_ = old[0, :-1], old[0, 1:], new[0, :-1], new[0, 1:]
        rp = (E_ + N_ - S_ - W_) / (2 * h)
        rm = (W_ + N_ - S_ - E_) / (2 * hm)
        ok = ~(dold[:-1].astype(bool) | dold[1:].astype(bool) | dnew[:-1].astype(bool) | dnew[1:].astype(bool))
        tr = ok & (rp < -tol) & (rm < -tol)
        if tr.any():
            return dict(fate=-1, j=j, i=int(np.argmax(tr)))
        i = int(((j + 0.5) * hm) / h)
        if i < Np and ok[i] and 0.25 * (S_[i] + E_[i] + W_[i] + N_[i]) >= rexp and rp[i] > tol and rm[i] > tol:
            return dict(fate=1, j=j, i=i)
        old, new = new, old
        dold, dnew = dnew, dold
    return dict(fate=0, j=None, i=None)


AN_H, AN_Q, AN_ROWS_X, AN_COLS_X, AN_TOL = 0.02, 32, 16.5, 20.5, 1e-10


def pstar_an_cmd(out, h=AN_H, Q=AN_Q):
    """p*_AN of family S on the anisotropic grid (h, h/Q), bisected to 1e-14 (the check AN compares at p*_AN - delta)."""
    hm = h / Q
    Np, nrows = int(round(AN_COLS_X / h)), int(round(AN_ROWS_X / hm))
    t0 = time.time()
    c = i40_ref(h)["ps"]
    for lo, hi in ((-6 * c, 6 * c), (-20 * c, 20 * c)):
        try:
            lo_, hi_, stopped = th.bisect(lambda p: fate_aniso(h, hm, p, Np, nrows), lo, hi, AN_TOL, lambda r: None)
            break
        except RuntimeError as e:
            print("bracket", lo, hi, "failed:", e, flush=True)
    else:
        raise SystemExit("no bracket")
    res = dict(h=h, Q=Q, lo=lo_, hi=hi_, ps=0.5 * (lo_ + hi_), stopped=stopped, tol=AN_TOL, secs=round(time.time() - t0, 1), square_ps=c)
    os.makedirs(out, exist_ok=True)
    json.dump(res, open(os.path.join(out, f"pstar_an_h{h:g}_Q{Q}.json"), "w"), indent=1)
    print(json.dumps(res), flush=True)


def an_grid_check(raw, h=AN_H, Q=AN_Q, delta=1e-6, Fstop=40.0):
    """Zoom-free control: one uniform anisotropic grid (h, h/Q) over the whole rectangle (same kernel, no boxes, no interpolation), at its own
    p*_AN - delta; the horizon series of the chain (at p*_h - delta) against it where e^{-f} on the horizon is below Fstop."""
    pa = json.load(open(os.path.join(raw, "runs", f"pstar_an_h{h:g}_Q{Q}.json")))
    ps = json.load(open(os.path.join(raw, "runs", f"pstar_h{h:g}.json")))["ps"]
    r = run_real(h, delta, ps)
    s = r["ser"]
    if len(s["F"]) == 0 or not (s["F"] < Fstop).any():
        return dict(h=h, Q=Q, delta=delta, ev=None, ev_chain=r["event"], n=0, dM_max=float("inf"), dX_max=float("inf"), F_to=0.0)
    A = pa["ps"] - delta
    pulses, sym, eps = th.family("S", A)
    hm = h / Q
    Np = r["Np"]
    Nm = int(round((float(s["xm"][s["F"] < Fstop].max()) + 0.8) / hm))
    d1 = dn.ellis_rays(CRIT["n"], Np, h, pulses=pulses, symmetric=sym, eps=eps)
    d2 = dn.ellis_rays(CRIT["n"], Nm, hm * (1.01 if mut("an_nodes") else 1.0), pulses=pulses, symmetric=sym, eps=eps)
    G, dead, ev = evolve_box(d1.Sp, d2.Sm, np.zeros(Np + 1, np.uint8), np.zeros(Nm + 1, np.uint8), h, hm, CRIT["n"], eps, CRIT["rcut"], detect=CRIT["tol"])
    if ev is None:
        return dict(h=h, Q=Q, delta=delta, ev=None, ev_chain=r["event"], n=0, dM_max=float("inf"), dX_max=float("inf"), F_to=0.0)
    L = Level(0, G, dead, h, hm, 0.0, 0.0)
    rec, stop, exh = track(L, CRIT["n"], ev[1], ev[0] + 0.5, stop_F=Fstop, tol=CRIT["tol"], win=4 * Q)
    if len(rec["r"]) < 10:
        return dict(h=h, Q=Q, delta=delta, ev=ev, ev_chain=r["event"], n=0, dM_max=float("inf"), dX_max=float("inf"), F_to=0.0)
    from scipy.interpolate import CubicSpline
    m = s["F"] < Fstop
    lo_x = max(rec["xp"][0], s["xp"][0]) + h
    hi_x = min(rec["xp"][-1], s["xp"][m][-1]) - h
    v = np.linspace(lo_x, hi_x, 4000)
    dM = CubicSpline(rec["xp"], rec["M"])(v) - CubicSpline(s["xp"], s["M"])(v)
    dX = CubicSpline(rec["xp"], rec["xm"])(v) - CubicSpline(s["xp"], s["xm"])(v)
    return dict(h=h, Q=Q, delta=delta, ev=ev, ev_chain=r["event"], n=len(v), dM_max=float(np.max(np.abs(dM))), dX_max=float(np.max(np.abs(dX))), F_to=float(rec["F"][-1]),
                M_AN_end=float(rec["M"][-1]), x_end=float(rec["xp"][-1]))


MUTATIONS = {           # name: the checks that must fire; computed in the development runs (raw I41-rpi5 dev/mut/), declared in its README before the runs of record
    "mass_exponent": ["PL1", "RC3", "RC4"],
    "finder_sign": ["AN", "PL1", "PL2", "RC2", "RC3", "RC4", "RC5", "RC6", "RC6b", "RC7", "RC8"],
    "hm_swap": ["AN", "PL1", "PL2", "RC4", "RC5", "RC6", "RC6b", "RC7", "RC8"],
    "interp_shift": ["AN", "PL1", "PL2", "RC6", "RC6b", "RC7", "RC8"],
    "flux_sign": ["RC4"],
    "time_unit": ["RC2"],
    "no_zoom": ["PL1", "RC5", "RC6", "RC6b", "RC7", "RC8"],
    "short_limit": ["RC5", "RC8"],
    "pstar_off": ["RC2", "RC3", "RC4", "RC5", "RC6", "RC6b", "RC7", "RC8"],
    "theta_coef": ["AN", "PL1", "PL2", "PL3", "RC2", "RC3", "RC4", "RC5", "RC6", "RC6b", "RC7", "RC8"],
    "rcut_off": ["PL3"],
    "an_nodes": ["AN"],
    "i40_shift": ["RC1"],
    "col_shift": ["PL2", "RC6", "RC6b"],
    "delta_swap": ["RC7"],
}


def settled_summary(s, h):
    """Late-time values of one series: M at the limit, plateau means over the last 1.0 and 0.5 in x+ and their range, the drift over the last 1.0,
    and the mean horizon flux |dM/dx+| over the last 0.5 (the first law's integrand)."""
    v, M = s["xp"], s["M"]
    out = dict(M_end=float(M[-1]), F_end=float(s["F"][-1]))
    for w in (1.0, 0.5):
        m = v > v[-1] - w
        out[f"plat{w:g}"] = float(M[m].mean())
        out[f"rng{w:g}"] = float(M[m].max() - M[m].min())
    m = v > v[-1] - 1.0
    out["drift1"] = float(M[m][0] - M[-1])
    mp = (v > v[-1] - 2.0) & (v <= v[-1] - 1.0)
    out["plat_prev"] = float(M[mp].mean())
    out["dplat"] = out["plat_prev"] - out["plat1"]            # change of the plateau mean between the last two windows of 1.0
    flux = 0.5 * s["r"] ** 2 / s["F"] * np.abs(s["rm"]) * s["phip"] ** 2
    out["flux05"] = float(flux[v > v[-1] - 0.5].mean())
    return out


def analyze_cmd(raw, an=True):
    t0 = time.time()
    hs, ds = CRIT["hs"], CRIT["deltas"]
    h3 = [h for h in hs if i40_ref(h)]                 # the three h of #40
    Fmax = CRIT["F1"] * 2.0 ** CRIT["K"]
    print("\n== Planted exact solutions through the chain")
    done = {c for c, _ in results}
    try:
        planted_checks()
    except Exception as e:
        for cid in ("PL1", "PL2"):
            if cid not in {c for c, _ in results}:
                check(cid, "planted case", False, "failed: " + repr(e))
    try:
        kernel_identity(raw)
    except Exception as e:
        check("PL3", "level 0 against dn.evolve", False, "failed: " + repr(e))
    if an:
        try:
            r = an_grid_check(raw)
        except Exception as e:
            r = dict(h=AN_H, delta=1e-6, ev=None, ev_chain=None, dM_max=float("inf"), dX_max=float("inf"), n=0, F_to=0.0, err=repr(e))
        check("AN", "chain against a zoom-free uniform anisotropic grid (h, h/32): horizon M and x- where e^{-f} < 40",
              r["dM_max"] <= 3e-4, f"h = {r['h']:g}, delta = {r['delta']:g}: events {r['ev']} / {r['ev_chain']}, max |dM| {r['dM_max']:.2e} (bar 3e-4), max |dx-| {r['dX_max']:.2e}, "
              f"{r['n']} columns, AN series to F = {r['F_to']:.4g}")
    print("\n== p*_h on this machine against #40's (dh610-fedora)")
    ok = True
    for h in hs:
        pj = json.load(open(os.path.join(raw, "runs", f"pstar_h{h:g}.json")))
        ref = i40_ref(h)
        d = pj["ps"] - ref["ps"] if ref else float("nan")
        print(f"    h = {h:g}: p*_h = {pj['ps']:.16e}, hi - lo = {pj['hi'] - pj['lo']:.2e}, minus #40's {d:+.2e}, {pj['runs']} runs, {pj['secs']:.0f} s")
        ok &= (pj["hi"] - pj["lo"] <= 1e-14) and (ref is None or abs(d) <= 1e-11)
    check("RC1", "p*_h bisected to hi - lo <= 1e-14 and, where #40 has it, within 1e-11 of #40's", ok)
    S, I = {}, {}
    fail = {}
    for dl in ds:
        for h in hs:
            try:
                S[(h, dl)], I[(h, dl)] = get_run(raw, h, dl)
                if len(S[(h, dl)]["r"]) == 0:
                    raise RuntimeError("empty series: " + str(I[(h, dl)]["status"]))
            except Exception as e:
                fail[(h, dl)] = repr(e)
                print(f"    run h = {h:g}, delta = {dl:g} failed: {e!r}")
    complete = not fail
    if complete:
        ev_ok = all(list(I[(h, 1e-8)]["event"]) == list(i40_ref(h)["event"]) for h in h3)
        t_ok = all(abs(t_of((I[(h, 1e-8)]["event"][1] + 0.5) * h, (I[(h, 1e-8)]["event"][0] + 0.5) * h) - i40_ref(h)["T"]) <= SQ2 * h for h in h3)
        fate_ok = all(I[(h, dl)]["point"]["fate"] == -1 and I[(h, dl)]["up"]["fate"] == 1 for h in hs for dl in ds)
        check("RC2", "fates (collapse at p*_h - delta, expansion at p*_h + delta; all h, both delta); the first trapped cell and its time against #40's at delta = 1e-8",
              fate_ok and ev_ok and t_ok, "events " + str({h: list(I[(h, 1e-8)]["event"]) for h in hs}) + " against #40's " + str({h: i40_ref(h)["event"] for h in h3}))
    else:
        check("RC2", "fates, first trapped cell and its time", False, "a run failed")
        for cid, txt in (("RC3", "formation value"), ("RC4", "first law"), ("RC5", "limit"), ("RC6", "convergence"), ("RC6b", "convergence, finer triple"), ("RC8", "settled"), ("RC7", "delta independence")):
            check(cid, txt, False, "a run failed")
        return finish(t0)
    dm = [abs(S[(h, 1e-8)]["M"][0] - i40_ref(h)["M_hor"]) for h in h3]
    check("RC3", "M_AH at the first trapped cell against #40's formation value at p*_h - 1e-8", max(dm) <= 2e-3,
          "|dM| " + ", ".join(f"{x:.1e}" for x in dm) + f" (bar 2e-3); M_AH {[round(float(S[(h, 1e-8)]['M'][0]), 5) for h in h3]}, #40 {[round(i40_ref(h)['M_hor'], 5) for h in h3]}")
    worst = 0.0
    for key, s in S.items():
        i1 = int(np.argmax(s["F"] >= 20.0))
        dMs = s["M"][i1] - s["M"][0]
        worst = max(worst, abs(dMs - flux_integral(s, 0, i1)) / abs(dMs))
    check("RC4", "first law on the horizon: M(x+_b) - M(x+_a) = integral of (r^2/2) e^f r_- phi_+^2 along it, from the first trapped cell to F = 20", worst <= 0.02,
          f"max relative difference {worst:.2e} over all runs (bar 2e-2)")
    ok5 = True
    det = []
    for key, s in S.items():
        h = key[0]
        fine = np.allclose(np.diff(s["xp"]), h, atol=1e-9)
        ok5 &= I[key]["status"] == "limit" and s["F"][-1] >= 0.5 * Fmax and s["F"][-1] < Fmax and s["edge"].min() >= I[key]["Nc"] / 8 and fine
        det.append(f"h={key[0]:g},d={key[1]:g}: F_end {s['F'][-1]:.0f}, edge {s['edge'].min():.0f}/{I[key]['Nc']}")
    check("RC5", f"every run ends at the declared limit (F = e^-f on the horizon, {0.5 * Fmax:g} <= F_end < {Fmax:g}), horizon at least Nc/8 cells from the box edges, no gaps", ok5, "; ".join(det))
    print("\n== Second-order convergence of M_AH(x+): max |M_h - M_h/2| over windows of x+ - x+_f (cubic splines)")
    for cid, trip in (("RC6", hs[:3]), ("RC6b", hs[1:4])):
        ok6 = True
        for dl in ds:
            rows = band_orders({h: S[(h, dl)] for h in trip}, trip, [(0.0, 1.5), (1.5, 3.0), (3.0, 99.0)], dl)
            for a_, b_, E, o in rows:
                good = all(BAND_R[0] <= x <= BAND_R[1] for x in o)
                ok6 &= good
                print(f"    {cid} h = {trip}, delta = {dl:g}, x+ - x+_f in [{a_:g}, {b_:g}]: E12 {E[0]:.2e}, E23 {E[1]:.2e}, order {o[0]:.2f} {'ok' if good else 'OUT'}")
        check(cid, f"M_AH(x+) converges at second order at the three h {trip}, three windows, both deltas (observed order in {BAND_R})", ok6)
    dr = {k: settled_summary(s, k[0]) for k, s in S.items()}
    check("RC8", "M_AH has settled at the limit: the plateau mean over the last 1.0 of x+ differs from the previous 1.0 by <= 1e-4 and the mean horizon flux over the last 0.5 is <= 1e-5, every run",
          all(abs(q["dplat"]) <= 1e-4 and q["flux05"] <= 1e-5 for q in dr.values()),
          "plateau change " + ", ".join(f"{q['dplat']:+.1e}" for q in dr.values()) + "; flux " + ", ".join(f"{q['flux05']:.1e}" for q in dr.values()))
    dd = [abs(dr[(h, 1e-6)]["plat1"] - dr[(h, 1e-8)]["plat1"]) for h in hs]
    check("RC7", "the settled value does not depend on delta: plateau means over the last 1.0 at 1e-6 and 1e-8 agree", max(dd) <= 1e-4, "|dM| " + ", ".join(f"{x:.1e}" for x in dd) + " (bar 1e-4)")
    print("\n== M_AH(x+) against x+ - x+_f (a = 1, 8 pi G = 1; M = r_AH / 2); t = (x+ + x-)/sqrt2 at the horizon; finest h")
    for dl in ds:
        s, hf = S[(hs[-1], dl)], hs[-1]
        print(f"  delta = {dl:g}, h = {hf:g}: x+_f = {s['xp'][0]:.4f}, t_f = {t_of(s['xp'][0], s['xm'][0]):.4f}")
        print("    dx+     x+       x-         t       r_AH     M_AH      e^-f       r_-")
        for dv in (0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.4):
            i = int(np.argmin(np.abs(s["xp"] - (s["xp"][0] + dv))))
            if abs(s["xp"][i] - s["xp"][0] - dv) > 2 * hf:
                break
            print(f"    {s['xp'][i] - s['xp'][0]:5.2f} {s['xp'][i]:8.4f} {s['xm'][i]:9.5f} {t_of(s['xp'][i], s['xm'][i]):8.4f} {s['r'][i]:8.5f} {s['M'][i]:8.5f} {s['F'][i]:10.4g} {s['rm'][i]:11.4g}")
    print("\n== Late-time values (reports)")
    print("  delta   h        M_end     plat(1.0)  range(1.0) plat(0.5)  drift(1.0)  plat change  flux(0.5)   2 M_end")
    SS = {}
    for dl in ds:
        for h in hs:
            q = settled_summary(S[(h, dl)], h)
            SS[(h, dl)] = q
            print(f"  {dl:6g} {h:7g}  {q['M_end']:.6f}  {q['plat1']:.6f}  {q['rng1']:.1e}  {q['plat0.5']:.6f}  {q['drift1']:+.2e}  {q['dplat']:+.2e}  {q['flux05']:.2e}  {2 * q['M_end']:.5f}")
    for dl in ds:
        for trip in (hs[:3], hs[1:4]):
            a = [SS[(h, dl)]["M_end"] for h in trip]
            p = [SS[(h, dl)]["plat1"] for h in trip]
            print(f"  delta = {dl:g}, h = {trip}: M_end {[round(x, 6) for x in a]}, differences {a[0] - a[1]:+.2e}, {a[1] - a[2]:+.2e}, Richardson (second order) M_end(0) = {a[2] - (a[1] - a[2]) / 3:.6f}; "
                  f"plateau(1.0) {[round(x, 6) for x in p]}, differences {p[0] - p[1]:+.2e}, {p[1] - p[2]:+.2e}, Richardson {p[2] - (p[1] - p[2]) / 3:.6f}")
    return finish(t0)


def finish(t0):
    nfail = sum(1 for _, ok in results if ok is False)
    print(f"\n{len(results)} checks, {nfail} failed; {time.time() - t0:.0f} s")
    if MUT is not None or dn.MUT is not None:
        fired = sorted({c for c, ok in results if ok is False})
        want = sorted(MUTATIONS.get(MUT, []))
        print(f"mutation {MUT}: fired {fired}, declared {want}: {'MATCH' if fired == want else 'DIFFERENT'}")
        return 0 if fired == want else 1
    return 1 if nfail else 0


# --------------------------------------------------------------------------------------------------------------------
# supplementary reports (no acceptance bar): the final mass against the distance from threshold, and [SH02]'s own data
# --------------------------------------------------------------------------------------------------------------------
def initial_energy(Sp, h):
    """[SH02] (3.3): E(x+ = X, x- = 0) from the ray data (D = 4, f = 0): E = r/2 + r r_+ r_-, with r_- integrated along the ray from
    r_{+-} = -(r_+ r_- + 1/2)/r, r_-(0) = 0 (Heun). Returns the array E(x+) at the nodes."""
    r = Sp[0]
    rp = np.gradient(r, h, edge_order=2)
    rm = np.zeros_like(r)
    for i in range(len(r) - 1):
        k1 = -(rp[i] * rm[i] + 0.5) / r[i]
        pred = rm[i] + h * k1
        k2 = -(rp[i + 1] * pred + 0.5) / r[i + 1]
        rm[i + 1] = rm[i] + 0.5 * h * (k1 + k2)
    return 0.5 * r + r * rp * rm


def sh02_rays(c_a, c_b, c_c, hp, Np, hm, Nm):
    """Characteristic data of [SH02] Sec. III.B in this code's variables (D = 4, a = 1, f = 0): on the '+' ray
    wp(x) = 1/sqrt(2 + x^2) + c_a exp(-c_b (l - c_c)^2), l = x/sqrt2, phi' = sqrt2 wp / r, r'' = wp^2 / r (RK4 from r = 1, r' = 0, phi = 0); the '-' ray
    is the static wormhole. The Gaussian is in l as printed in [SH02] (3.1); with it E0 below reproduces their Table I to a constant factor 0.963 (a Gaussian
    in x+ at c_c does not: 0.83 to 0.85)."""
    def wp(x):
        return 1.0 / np.sqrt(2.0 + x ** 2) + c_a * np.exp(-c_b * (x / SQ2 - c_c) ** 2)

    def rhs(x, y):
        return np.array([y[1], wp(x) ** 2 / y[0], SQ2 * wp(x) / y[0]])
    Sp = np.zeros((3, Np + 1))
    y = np.array([1.0, 0.0, 0.0])
    Sp[0, 0], Sp[2, 0] = y[0], y[2]
    for i in range(Np):
        x = i * hp
        k1 = rhs(x, y)
        k2 = rhs(x + 0.5 * hp, y + 0.5 * hp * k1)
        k3 = rhs(x + 0.5 * hp, y + 0.5 * hp * k2)
        k4 = rhs(x + hp, y + hp * k3)
        y = y + hp / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        Sp[0, i + 1], Sp[2, i + 1] = y[0], y[2]
    d2 = dn.ellis_rays(2, Nm, hm, pulses=[], eps=(-1.0,))
    return Sp, d2.Sm


def e0_of_family(h, A):
    """[SH02] (3.3) for family S at amplitude A: E(20, 0) minus the static value, from the ray data."""
    pulses, sym, eps = th.family("S", A)
    N = int(round(20.0 / h))
    d = dn.ellis_rays(CRIT["n"], N, h, pulses=pulses, symmetric=sym, eps=eps)
    d0 = dn.ellis_rays(CRIT["n"], N, h, pulses=[], eps=eps)
    return float(initial_energy(d.Sp, h)[-1] - initial_energy(d0.Sp, h)[-1])


def scan_cmd(h, out, pstar_dir):
    """Supplementary: M at the limit against the distance from threshold, family S (no acceptance bar)."""
    ps = json.load(open(os.path.join(pstar_dir, f"pstar_h{h:g}.json")))["ps"]
    os.makedirs(out, exist_ok=True)
    rows = []
    print("  delta        A = p*_h - delta     E0 (SH02 3.3)    x+_f     M_start   M_end     r_AH(end)   F_end   status")
    for dl in (1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-2, 1e-1):
        try:
            r = run_real(h, dl, ps, need_up=False)
            s = r["ser"]
            E0 = e0_of_family(h, r["A"])
            row = dict(delta=dl, A=r["A"], E0=E0, xpf=float(s["xp"][0]), M_start=float(s["M"][0]), M_end=float(s["M"][-1]), r_end=float(s["r"][-1]), F_end=float(s["F"][-1]), status=r["status"])
        except Exception as e:
            row = dict(delta=dl, status="failed: " + repr(e))
        rows.append(row)
        print("  " + "  ".join(f"{row[k]:.4g}" if isinstance(row.get(k), float) else str(row.get(k, "")) for k in ("delta", "A", "E0", "xpf", "M_start", "M_end", "r_end", "F_end", "status")), flush=True)
    json.dump(rows, open(os.path.join(out, f"scan_h{h:g}.json"), "w"), indent=1)


SH02_ROWS = [(-1e-4, 3.0, 4.39e-5, 8.42, 0.296), (-1e-3, 3.0, 4.38e-4, 6.98, 0.310), (-1e-2, 3.0, 4.32e-3, 5.64, 0.325), (-1e-1, 3.0, 3.70e-2, 4.46, 0.420),
             (-1e-1, 6.0, 2.60e-2, 4.82, 0.38), (-1e-1, 9.0, 2.12e-2, 4.90, 0.35)]


def sh02_cmd(h, out):
    """Supplementary: [SH02]'s Table I data (Sec. III.B) through the chain: E0, the final horizon x-_H and M, against their printed values."""
    os.makedirs(out, exist_ok=True)
    print("  c_a      c_b   E0 here    E0 SH02    x-_H here  x-_H SH02   M_end here  M SH02   x+_f   F_end   status   F* (our e^-f where M_AH = M SH02)  dx+ since formation")
    rows = []
    for c_a, c_b, E0p, xHp, Mp in SH02_ROWS:
        Np = int(round((xHp + 13.0) / h))
        nm0 = int(round((xHp + 3.5 + nc_of(h) * h / 4) / h))
        n20 = int(round(20.0 / h))
        Nray = max(Np, nm0, n20)
        Sp, Sm = sh02_rays(c_a, c_b, 3.0, h, Nray, h, Nray)
        E0 = float(initial_energy(Sp[:, :n20 + 1], h)[-1] - initial_energy(sh02_rays(0.0, c_b, 3.0, h, n20, h, 4)[0], h)[-1])
        o = chain(CRIT["n"], np.array([-1.0]), h, Sp[:, :Np + 1], Sm, 0.0, 0.0, nm0, None, None, CRIT["F1"], CRIT["F1"] * 2.0 ** CRIT["K"], nc_of(h), tol=CRIT["tol"], rcut=CRIT["rcut"], log=lambda *a: None)
        s = o["ser"]
        Mser = s["M"]
        i_ = int(np.argmax(Mser <= Mp)) if len(Mser) and (Mser <= Mp).any() else 0
        if i_ > 0:
            Fs = float(np.interp(Mp, Mser[i_ - 1:i_ + 1][::-1], s["F"][i_ - 1:i_ + 1][::-1]))
            dxs = float(np.interp(Mp, Mser[i_ - 1:i_ + 1][::-1], s["xp"][i_ - 1:i_ + 1][::-1]) - s["xp"][0])
        else:
            Fs, dxs = float("nan"), float("nan")
        row = dict(c_a=c_a, c_b=c_b, E0=E0, E0_SH02=E0p, F_star=Fs, dx_star=dxs, xH=float(s["xm"][-1]) if len(s["xm"]) else None, xH_SH02=xHp, M_end=float(s["M"][-1]) if len(s["M"]) else None, M_SH02=Mp,
                   xpf=float(s["xp"][0]) if len(s["xp"]) else None, F_end=float(s["F"][-1]) if len(s["F"]) else None, status=o["status"])
        rows.append(row)
        print(f"  {c_a:7g} {c_b:4g}  {E0:.3e}  {E0p:.3e}  {row['xH']:9.4f}  {xHp:9.2f}   {row['M_end']:.5f}   {Mp:.3f}  {row['xpf']:.3f}  {row['F_end']:.0f}  {o['status']}   {Fs:.2f}   {dxs:.2f}", flush=True)
    json.dump(rows, open(os.path.join(out, f"sh02_h{h:g}.json"), "w"), indent=1)


def sh02_uniform_cmd(out, N=9600, X=20.0):
    """Supplementary: [SH02]'s Table I data on one uniform grid of their size (N = 9600 over x+ in [0, 20], h = 2.08e-3), this code's kernel, no zoom and no
    excision of large e^{-f}: the horizon mass the uniform grid reports at x+ = 12, 16, 20 (where the finder still finds it), against the chain's and theirs."""
    os.makedirs(out, exist_ok=True)
    h = X / N
    print(f"  uniform grid h = {h:.5f}; c_a c_b: M_AH on the uniform grid at x+ = 12, 16, 20 (e^-f there), the chain's M_end, SH02's M")
    rows = []
    for c_a, c_b, E0p, xHp, Mp in SH02_ROWS:
        Nm = int(round((xHp + 1.5) / h))
        Sp, Sm = sh02_rays(c_a, c_b, 3.0, h, N, h, Nm)
        G, dead, ev = evolve_box(Sp, Sm, np.zeros(N + 1, np.uint8), np.zeros(Nm + 1, np.uint8), h, h, CRIT["n"], np.array([-1.0]), CRIT["rcut"], detect=CRIT["tol"])
        row = dict(c_a=c_a, c_b=c_b, M_SH02=Mp, event=ev)
        if ev is not None:
            L = Level(0, G, dead, h, h, 0.0, 0.0)
            rec, stop, exh = track(L, CRIT["n"], ev[1], ev[0] + 0.5, win=400)
            for xq in (12.0, 16.0, 20.0):
                if len(rec["xp"]) and rec["xp"][-1] >= xq - h:
                    i = int(np.argmin(np.abs(rec["xp"] - xq)))
                    row[f"M{int(xq)}"], row[f"F{int(xq)}"] = float(rec["M"][i]), float(rec["F"][i])
                else:
                    row[f"M{int(xq)}"], row[f"F{int(xq)}"] = None, None
            row["x_last"], row["F_last"], row["M_last"] = (float(rec["xp"][-1]), float(rec["F"][-1]), float(rec["M"][-1])) if len(rec["xp"]) else (None, None, None)
        rows.append(row)
        print("  " + json.dumps(row), flush=True)
        del G, dead
    json.dump(rows, open(os.path.join(out, f"sh02_uniform_N{N}.json"), "w"), indent=1)


def main():
    global MUT
    argv = sys.argv[1:]
    opts = dict(a.lstrip("-").split("=", 1) if "=" in a else (a.lstrip("-"), True) for a in argv[1:])
    cmd = argv[0] if argv else "help"
    if "mutate" in opts:
        MUT = opts["mutate"]
        dn.MUT = MUT if MUT in ("theta_coef", "ghost_sign", "scalar_coef") else None
    import platform
    la = open("/proc/loadavg").read().split()[:3] if platform.system() == "Linux" else ["?"]
    print(f"settled_mass.py {cmd} {opts}: numba {nb.__version__}, numpy {np.__version__}, python {platform.python_version()}, {platform.node()}, "
          f"load average at start {' '.join(la)}, mutation {MUT}, {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    if cmd == "pstar":
        pstar_cmd(float(opts["h"]), opts["out"])
    elif cmd == "job":
        job_cmd(float(opts["h"]), float(opts["delta"]), opts["out"], opts["pstar"])
    elif cmd == "planted":
        planted_checks()
        sys.exit(1 if any(ok is False for _, ok in results) else 0)
    elif cmd == "analyze":
        sys.exit(analyze_cmd(opts["raw"], an=not opts.get("noan")))
    elif cmd == "sh02_uniform":
        sh02_uniform_cmd(opts["out"])
    elif cmd == "pstar_an":
        pstar_an_cmd(opts["out"])
    elif cmd == "scan":
        scan_cmd(float(opts["h"]), opts["out"], opts["pstar"])
    elif cmd == "sh02":
        sh02_cmd(float(opts["h"]), opts["out"])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()

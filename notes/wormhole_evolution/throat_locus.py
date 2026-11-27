#!/usr/bin/env python3
"""Lifetime of the one-sided threshold family G along the throat locus (#225), and the expanding side on the locus's proper time.
Uses dn.py, checks.py, threshold.py and expansion.py of #39, #40 and #213 unchanged (imported, never modified).

    python3 throat_locus.py ...   (commands: see main(); criteria are CRIT and MUTATIONS below, as declared in negative-mass-raw I225-dh610/README.md)

Conventions (a = 1, 8 pi G = 1, t = (x+ + x-)/sqrt2, l = (x+ - x-)/sqrt2, metric -2 e^{-f} dx+ dx- + r^2 dOmega_n^2):
    the throat locus  = on each slice t = const (the anti-diagonal i + j = s of the grid nodes), the minimum of r with a parabolic refinement
                        (dn.throat_min's definition), its position l(t), its drift v = dl/dt;
    proper time on it = d tau = e^{-f/2} sqrt(dt^2 - dl^2) = e^{-f/2} sqrt(1 - v^2) dt   (on the diagonal this is dn.finish_diag's dtau).
The sweep is dn.row_step unchanged, row by row, on the triangle i + j <= s_max (or the whole rectangle), with the locus of each slice tracked while the
rows are produced: no grid is stored. #40's own events (first trapped cell in a row, expansion on the diagonal) are evaluated in the same loop, so
that #40's diagonal lifetime is reproduced by the same runs (the control).
"""
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import checks       # noqa: E402
import dn           # noqa: E402

SQ2 = math.sqrt(2.0)
W_EXACT = checks.W_EXACT
LMUT = None         # this file's mutation (dn.MUT holds dn.py's)


def lmut(name):
    return LMUT == name


# --------------------------------------------------------------------------------------------------------------------
# data: the static Ellis wormhole with the pair rescaling x+ -> cp x+, x- -> cm x- (cp cm = 1 is a gauge map; cp = cm = 1 is dn.ellis_rays)
# --------------------------------------------------------------------------------------------------------------------
def scaled_rays(n, N, h, lam=1.0, pulses=(), eps=(-1.0,), symmetric=False):
    """dn.ellis_rays with the whole solution mapped by x+ = x~+/lam, x- = lam x~-, that is data(s) = data~(lam s) on the ray x- = 0 and data~(s/lam) on the
    ray x+ = 0. The metric -2 e^{-f} dx+ dx- is invariant (f unchanged), so this is a coordinate change of the same spacetime: the static throat, at
    x~+ = x~-, moves in the grid with the velocity v = (1 - lam^2)/(1 + lam^2) (l = (x+ - x-)/sqrt2, t = (x+ + x-)/sqrt2), the lapse f stays 0.
    pulses: as dn.ellis_rays, given in the unscaled variable (a pulse centred at xc is centred at xc/lam on the ray x- = 0 and at xc lam on x+ = 0).
    The constraint r'' = -(S/n) r is solved by RK4 in s, S = sum eps phi_s^2, as dn.ellis_rays does. lam = 1 reproduces dn.ellis_rays bit for bit."""
    eps = list(eps)
    k = len(eps)
    pulses = [dict(p) for p in pulses]
    if symmetric:
        for p in list(pulses):
            if p["ray"] == "+":
                q = dict(p)
                q["ray"] = "-"
                q["A"] = p["A"] * (-1.0 if p["field"] == 0 else 1.0)
                pulses.append(q)
    sc = {"+": lam, "-": 1.0 / lam}
    sol, q = dn.ellis_dense(n, N * h * max(sc.values()) / SQ2 + 1.0)
    q = q * (1.1 if dn.mut("ellis_charge") else 1.0)
    cs = -1.0 if dn.mut("data_constraint_sign") else 1.0

    def phi_s(s, ray, field):
        c = sc[ray]
        out = np.zeros_like(s)
        if field == 0:
            l = c * s / SQ2
            r = sol(l)[0]
            out = out + (1.0 if ray == "+" else -1.0) * c * q / r ** n / SQ2
        for p in pulses:
            if p["ray"] == ray and p["field"] == field:
                out = out + c * dn._pulse(c * s, p["A"], p["xc"], p["w"], p.get("kind", "phi"))[1]
        return out

    def phi_val(s, ray, field):
        c = sc[ray]
        out = np.zeros_like(s)
        if field == 0:
            out = out + (1.0 if ray == "+" else -1.0) * sol(c * s / SQ2)[2]
        for p in pulses:
            if p["ray"] == ray and p["field"] == field:
                out = out + dn._pulse(c * s, p["A"], p["xc"], p["w"], p.get("kind", "phi"))[0]
        return out

    rays = []
    nodes = np.arange(N + 1) * h
    for ray in ("+", "-"):
        def S_of(s):
            return sum(eps[f] * phi_s(s, ray, f) ** 2 for f in range(k))
        S0, Sh, S1 = S_of(nodes[:-1]), S_of(nodes[:-1] + 0.5 * h), S_of(nodes[1:])
        r = np.empty(N + 1)
        v = np.empty(N + 1)
        r[0], v[0] = 1.0, 0.0
        for i in range(N):
            y = (r[i], v[i])

            def d_(y, S):
                return (y[1], -cs * (S / n) * y[0])
            k1 = d_(y, S0[i])
            y2 = (y[0] + 0.5 * h * k1[0], y[1] + 0.5 * h * k1[1])
            k2 = d_(y2, Sh[i])
            y3 = (y[0] + 0.5 * h * k2[0], y[1] + 0.5 * h * k2[1])
            k3 = d_(y3, Sh[i])
            y4 = (y[0] + h * k3[0], y[1] + h * k3[1])
            k4 = d_(y4, S1[i])
            r[i + 1] = y[0] + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            v[i + 1] = y[1] + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        U = np.zeros((2 + k, N + 1))
        U[0] = r
        for f in range(k):
            U[2 + f] = phi_val(nodes, ray, f)
        rays.append(U)
    return dn.Data(n, eps, N, h, rays[0], rays[1], name="ellis_scaled")


def scaled_exact(n, d, lam):
    """Exact static fields at the nodes of d's rectangle, mapped by the pair rescaling: r(x+, x-) = r_s(l~), l~ = (x+ lam - x-/lam)/sqrt2, f = 0, phi = +-phi_s."""
    sol, q = dn.ellis_dense(n, d.N * d.h * max(lam, 1.0 / lam) / SQ2 + 1.0)
    i = np.arange(d.N + 1)[:, None]
    j = np.arange(d.N + 1)[None, :]
    l = ((i * lam) - (j / lam)) * d.h / SQ2
    y = sol(np.abs(l).ravel())
    return y[0].reshape(l.shape), np.sign(l) * y[2].reshape(l.shape)


# --------------------------------------------------------------------------------------------------------------------
# the sweep: dn.row_step row by row, the locus of every slice t = const tracked on the way
# --------------------------------------------------------------------------------------------------------------------
def sweep(d, smax=None, rect=False, tol=1e-2, rexp=2.0, rcut=0.05, fcut=25.0, col_levels=(0.75, 0.5), exp_levels=(2.0, 3.0), stop=True, s_stop=None,
          prof_s=(), log=None, need_loc=True):
    """Evolve d over the triangle i + j <= smax (rect=True: the whole rectangle) with dn.row_step (unchanged) and return a dict:
    the locus series (one entry per slice s = i + j with a finite minimum away from the ends and no dead node on it), #40's events, the diagonal series, the
    profiles on the slices prof_s (all nodes of the slice), the level crossings of the locus's r_min.
    stop: end the sweep once #40's event and the locus events (a collapse level, or both expansion levels; need_loc=False: #40's event alone) are found,
    and s >= s_stop if given."""
    N, h, nf, eps, n = d.N, d.h, d.nf, d.eps, d.n
    kk = nf - 2
    par = dn.kernel_params(n, eps, 1)
    smax = N if (rect or smax is None) else min(int(smax), N)
    nrows = (N if rect else smax) + 1
    nS = 2 * N + 1 if rect else smax + 1
    nan = np.nan
    cr = np.full(nS, np.inf)
    crow = np.full(nS, -1, np.int64)
    ci = np.zeros(nS, np.int64)
    cf = np.full(nS, nan)
    crp = np.full(nS, nan)
    crm = np.full(nS, nan)
    cphi = np.full((max(kk, 1), nS), nan)
    cphp = np.full((max(kk, 1), nS), nan)
    cphm = np.full((max(kk, 1), nS), nan)
    pr = np.full(nS, nan)
    pf = np.full(nS, nan)
    nr = np.full(nS, nan)
    nfv = np.full(nS, nan)
    sdead = np.zeros(nS, bool)
    ser = {k: [] for k in ("s", "t", "l", "r", "f", "ic", "dk", "rp", "rm")}
    ser["phi"], ser["php"], ser["phm"] = [], [], []
    cross = {}
    ev40 = {}
    dg = {k: np.full(nrows, nan) for k in ("r", "f", "rp", "rm")}
    dg["cls"] = np.zeros(nrows, np.int8)
    prof = {int(s): dict(r=np.full(s + 1, nan), f=np.full(s + 1, nan), phi=np.full((max(kk, 1), s + 1), nan), rp=np.full(s + 1, nan), rm=np.full(s + 1, nan),
                         php=np.full((max(kk, 1), s + 1), nan), phm=np.full((max(kk, 1), s + 1), nan), dead=np.zeros(s + 1, bool)) for s in prof_s}
    lev_c, lev_e = sorted(col_levels, reverse=True), sorted(exp_levels)
    state = dict(done=False, rows=0)

    def finalize(s):
        lo, hi = (max(0, s - N), min(N, s)) if rect else (0, s)
        if sdead[s] or crow[s] < 0:
            return
        ic = int(ci[s])
        if ic <= lo or ic >= hi:
            return
        y1, y0, y2 = cr[s], nr[s], pr[s]
        if not (np.isfinite(y0) and np.isfinite(y2) and np.isfinite(y1)):
            return
        den = y0 - 2 * y1 + y2
        dk = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        rmin = y1 - 0.25 * (y0 - y2) * dk
        f0, f1, f2 = nfv[s], cf[s], pf[s]
        fl = f1 + dk * 0.5 * (f2 - f0) + 0.5 * dk * dk * (f0 - 2 * f1 + f2)
        ip = ic + dk
        ser["s"].append(s)
        ser["t"].append(s * h / SQ2)
        ser["l"].append((ip - (s - ip)) * h / SQ2)
        ser["r"].append(rmin)
        ser["f"].append(fl)
        ser["ic"].append(ic)
        ser["dk"].append(dk)
        ser["rp"].append(crp[s])
        ser["rm"].append(crm[s])
        ser["phi"].append(cphi[:kk, s].copy())
        ser["php"].append(cphp[:kk, s].copy())
        ser["phm"].append(cphm[:kk, s].copy())
        k = len(ser["s"]) - 1
        for lv in lev_c:
            if ("c", lv) not in cross and rmin <= lv:
                cross[("c", lv)] = k
        for lv in lev_e:
            if ("e", lv) not in cross and rmin >= lv:
                cross[("e", lv)] = k

    def loc_done():
        if any(("c", lv) in cross for lv in lev_c[-1:]):        # the deepest collapse level
            return True
        return all(("e", lv) in cross for lv in lev_e)

    old = np.ascontiguousarray(d.Sp[:, :(N if rect else smax) + 1].copy())
    dold = np.zeros(old.shape[1], np.uint8)
    older = None

    def node_row(jn, row, drow, oldr, olderr):
        """bookkeeping of the nodes of row jn: row (nf, K+1), oldr = row jn-1 (nf, >= K+2 or K+1 in rect mode), olderr = row jn-2."""
        K = row.shape[1] - 1
        R, F, P = row[0], row[1], row[2:]
        idx = np.arange(K + 1)
        s_idx = idx + jn
        if jn >= 1:
            sel = crow[s_idx] == jn - 1
            if sel.any():
                nr[s_idx[sel]] = R[idx[sel]]
                nfv[s_idx[sel]] = F[idx[sel]]
            s_c = s_idx - 1
            sel2 = crow[s_c] == jn - 1
            if sel2.any() and olderr is not None:
                ii = idx[sel2]
                ss = s_c[sel2]
                crm[ss] = (R[ii] - olderr[0][ii]) / (2 * h)
                if kk:
                    cphm[:kk, ss] = (P[:, ii] - olderr[2:][:, ii]) / (2 * h)
        ok = drow == 0
        sdead[s_idx[~ok]] = True
        better = ok & (R < cr[s_idx])
        if better.any():
            ii = idx[better]
            ss = s_idx[better]
            cr[ss] = R[ii]
            crow[ss] = jn
            ci[ss] = ii
            cf[ss] = F[ii]
            Rp = np.concatenate([[nan], R, [nan]])
            crp[ss] = (Rp[ii + 2] - Rp[ii]) / (2 * h)
            crm[ss] = nan
            if kk:
                Pp = np.concatenate([np.full((kk, 1), nan), P, np.full((kk, 1), nan)], axis=1)
                cphi[:kk, ss] = P[:, ii]
                cphp[:kk, ss] = (Pp[:, ii + 2] - Pp[:, ii]) / (2 * h)
                cphm[:kk, ss] = nan
            nr[ss] = nan
            nfv[ss] = nan
            if oldr is not None:
                w = oldr.shape[1]
                jj = ii + 1
                pr[ss] = np.where(jj < w, oldr[0][np.minimum(jj, w - 1)], nan)
                pf[ss] = np.where(jj < w, oldr[1][np.minimum(jj, w - 1)], nan)
            else:
                pr[ss] = nan
                pf[ss] = nan
        # profiles
        for s_k, pk in prof.items():
            if jn <= s_k:
                i = s_k - jn
                if i <= K:
                    pk["r"][i], pk["f"][i] = R[i], F[i]
                    pk["dead"][i] = bool(drow[i])
                    if kk:
                        pk["phi"][:kk, i] = P[:, i]
                    if 1 <= i <= K - 1:
                        pk["rp"][i] = (R[i + 1] - R[i - 1]) / (2 * h)
                        if kk:
                            pk["php"][:kk, i] = (P[:, i + 1] - P[:, i - 1]) / (2 * h)
            if jn >= 2 and jn - 1 <= s_k:
                i = s_k - (jn - 1)
                if i <= K and olderr is not None:
                    pk["rm"][i] = (R[i] - olderr[0][i]) / (2 * h)
                    if kk:
                        pk["phm"][:kk, i] = (P[:, i] - olderr[2:][:, i]) / (2 * h)

    node_row(0, old, dold, None, None)
    if not rect:
        pass
    jn = 0
    for jn in range(1, nrows):
        K = N if rect else smax - jn
        old_c = old[:, :K + 1] if rect else np.ascontiguousarray(old[:, :K + 1])
        dold_c = dold[:K + 1] if rect else np.ascontiguousarray(dold[:K + 1])
        new = np.empty((nf, K + 1))
        dnew = np.zeros(K + 1, np.uint8)
        new[:, 0] = d.Sm[:, jn]
        dn.row_step(old_c, new, dold_c, dnew, h, par, eps, rcut, fcut, 1)
        c, up, um = dn.strip_cells(old_c, new, h)
        alive = ~(dold_c[:-1].astype(bool) | dold_c[1:].astype(bool) | dnew[:-1].astype(bool) | dnew[1:].astype(bool))
        j = jn - 1
        if j < K:
            if alive[j]:
                dg["r"][j], dg["f"][j], dg["rp"][j], dg["rm"][j] = c[0, j], c[1, j], up[0, j], um[0, j]
            dg["cls"][j] = dn.cell_class(up[0, j], um[0, j], tol) if alive[j] else 0
        if not ev40:
            cc = dn.cell_class(up[0], um[0], tol)
            trapped = alive & (cc == -1)
            if trapped.any():
                ev40["collapse"] = (j, int(np.argmax(trapped)))
            elif j < K and alive[j] and cc[j] == 1 and c[0, j] >= rexp:
                ev40["expand"] = (j, j)
        node_row(jn, new, dnew, old, older)
        if not rect:
            finalize(jn)
        elif jn <= N:
            finalize(jn)
        older, old, dold = old, new, dnew
        state["rows"] = jn + 1
        if log and jn % max(1, nrows // 10) == 0:
            log(f"sweep row {jn} of {nrows - 1}")
        if stop and ev40 and (loc_done() or not need_loc) and (s_stop is None or jn >= s_stop):
            break
    if rect:
        for s in range(N + 1, 2 * N):
            finalize(s)
    out = {k: np.array(v) for k, v in ser.items() if k not in ("phi", "php", "phm")}
    for k in ("phi", "php", "phm"):
        out[k] = np.array(ser[k]).reshape(len(ser["s"]), kk)
    out.update(N=N, h=h, n=n, eps=np.array(eps, float), smax=smax, rect=rect, ev40=ev40, cross={f"{a}:{b}": v for (a, b), v in cross.items()},
               diag=dg, prof=prof, rows=state["rows"])
    return out


def ev40_time(L, key):
    """#40's T for the event (dn.event_time with x0 = 0)."""
    j, i = L["ev40"][key]
    h = L["h"]
    return ((i + 0.5) * h + (j + 0.5) * h) / SQ2


# --------------------------------------------------------------------------------------------------------------------
# the locus series: proper time, drift, H, events
# --------------------------------------------------------------------------------------------------------------------
CRIT = dict(
    planted_tol=1e-13,            # bisection tolerance of the planted series
    vwin=0.5,                     # time window (each side) of the local linear fit of l(t) that gives the drift v = dl/dt
)


def locus_series(L, vwin=None):
    """From the sweep's raw locus: drift v(t) by a local linear fit of l(t) over +- vwin in t, proper time tau along the locus
    (d tau = e^{-f/2} sqrt(1 - v^2) dt, f at the locus, tau starting from t at the first valid slice), and H = d ln r_min / d tau."""
    vwin = CRIT["vwin"] if vwin is None else vwin
    t, l, r, f = L["t"], L["l"], L["r"], L["f"]
    m = len(t)
    v = np.full(m, np.nan)
    if m >= 3:
        dts = np.median(np.diff(t))
        k = max(1, int(round(vwin / dts)))
        for i in range(m):
            a, b = max(0, i - k), min(m, i + k + 1)
            if b - a >= 3:
                v[i] = np.polyfit(t[a:b], l[a:b], 1)[0]
    # tau_diag: the diagonal's proper time (v = 0) instead of the locus's; tau_nolapse: f dropped from d tau
    vv = np.zeros_like(v) if lmut("tau_diag") else v
    lapse = np.zeros_like(f) if lmut("tau_nolapse") else f
    g = np.exp(-lapse / 2) * np.sqrt(np.clip(1 - vv ** 2, 0, None))       # d tau / d t along the locus
    tau = (t[0] * g[0] + np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(t))])) if m else g
    sc = SQ2 if lmut("time_unit") else 1.0
    L.update(v=v, tau=tau * sc, dtaudt=g)
    return L


def crossing(L, key):
    """Level crossing of r_min: (t, tau) by linear interpolation between the two slices that bracket the level, or None."""
    k = L["cross"].get(key)
    if k is None:
        return None
    lv = float(key.split(":")[1])
    if k == 0:
        return L["t"][0], L["tau"][0]
    r0, r1 = L["r"][k - 1], L["r"][k]
    w = (lv - r0) / (r1 - r0) if r1 != r0 else 1.0
    return (L["t"][k - 1] + w * (L["t"][k] - L["t"][k - 1]), L["tau"][k - 1] + w * (L["tau"][k] - L["tau"][k - 1]))


# --------------------------------------------------------------------------------------------------------------------
# the static ghost wormholes with a lapse (two-parameter family) and the asymmetry from a slice
# --------------------------------------------------------------------------------------------------------------------
def static_member(n, rth, sth, rmax=40.0):
    """ds^2 = -N^2 dt^2 + dl^2 + r^2 dOmega_n^2 with a massless ghost scalar (8 pi G = 1): N' = K/r^n, r r'' = -(n-1)(r'^2 - 1) - r r' N'/N (the angular
    and the lapse equations; the l l equation fixes the charge q^2 = n (n-1) N_th^2 r_th^{2(n-1)}). Parameters: the throat radius rth and
    sth = r_th N'/N at the throat (= gamma/sqrt(1 + gamma^2) in GGS09a's family, D = 4). N_th = 1. Returns f(sign) -> (l, r, r', N) arrays on the side
    sign = +1 (l > 0, K as given) or -1 (l < 0, by l -> -l and K -> -K), integrated until r = rmax, and the helper M(r) = r^{n-1} (1 - r'^2)/2."""
    from scipy.integrate import solve_ivp
    K = sth * rth ** (n - 1)

    def side(sg):
        Ks = sg * K

        def rhs(l, y):
            r, p, N = y
            return [p, (-(n - 1.0) * (p * p - 1.0) - r * p * Ks / (r ** n * N)) / r, Ks / r ** n]
        ev = lambda l, y: y[0] - rmax      # noqa: E731
        ev.terminal = True
        sol = solve_ivp(rhs, [0.0, 400.0], [rth, 0.0, 1.0], method="DOP853", rtol=1e-12, atol=1e-14, events=ev, dense_output=True)
        l = np.linspace(0.0, sol.t[-1], 4001)
        y = sol.sol(l)
        return l, y[0], y[1], y[2]
    return side


def member_M(n, rth, sth, rs, rmax=None):
    """M at areal radii rs (arrays) on the + side and on the - side of the member: M = r^{n-1}(1 - r'^2)/2 (the Misner-Sharp mass of a static state)."""
    rs = np.asarray(rs, float)
    side = static_member(n, rth, sth, rmax=rmax or float(rs.max()) * 1.05 + 1.0)
    out = []
    for sg in (1, -1):
        l, r, p, N = side(sg)
        M = r ** (n - 1) * (1.0 - p * p) / 2
        out.append(np.interp(rs, r, M))
    return out[0], out[1]


def asym_masses(n, rth, sth):
    """The two asymptotic masses m_+, m_- of the member (the limit of M(r) as r -> infinity, Richardson in 1/r: M = m + c/r^{n-1}... fitted on two large radii)."""
    out = []
    for sg in (0, 1):        # M(r) = m + c1/r^{n-1} + c2/r^{2(n-1)} on three large radii
        rs = np.array([40.0, 80.0, 160.0])
        Ms = member_M(n, rth, sth, rs, rmax=170.0)[sg]
        A = np.array([[1.0, rr ** -(n - 1.0), rr ** -(2 * n - 2.0)] for rr in rs])
        out.append(float(np.linalg.solve(A, Ms)[0]))
    return out[0], out[1]


def ggs_closed(g):
    """GGS09a's D = 4 family in b = 1: r_th, s_th = r_th N'/N, m_+, m_- (B = r_th e^{g atan g}/sqrt(1 + g^2))."""
    rth = math.sqrt(1 + g * g) * math.exp(-g * math.atan(g))
    B = rth * math.exp(g * math.atan(g)) / math.sqrt(1 + g * g)
    return rth, g / math.sqrt(1 + g * g), B * g * math.exp(-g * math.pi / 2), -B * g * math.exp(g * math.pi / 2)


def fit_asym(n, rth, r_pts, Mp, Mm):
    """Least-squares s_th of the member (given its throat radius) from M on the two sides at the areal radii r_pts: minimizes the sum of (DM_pred - DM_meas)^2,
    DM = M_+ - M_-, over s_th in (-0.95, 0.95). Returns (s_th, rms residual of DM)."""
    from scipy.optimize import minimize_scalar
    r_pts = np.asarray(r_pts, float)
    dm = np.asarray(Mp) - np.asarray(Mm)
    sgn = -1.0 if lmut("asym_sign") else 1.0

    def cost(s):
        a, b = member_M(n, rth, s, r_pts)
        return float(np.sum((sgn * (a - b) - dm) ** 2))
    o = minimize_scalar(cost, bounds=(-0.95, 0.95), method="bounded", options=dict(xatol=1e-9))
    return float(o.x), math.sqrt(o.fun / len(r_pts))


def slice_masses(L, key, n, r_pts):
    """From a profile slice of the sweep (L['prof'][s]) and its locus: M = r^{n-1}(1 + 2 e^f r_+ r_-)/2 at the nodes, and M(r) on each side of the locus
    (the node of the minimum of r), interpolated at the areal radii r_pts. Returns (Mp, Mm, r_th, node of the minimum, per-node arrays)."""
    pk = L["prof"][key]
    r, f, rp, rm = pk["r"], pk["f"], pk["rp"], pk["rm"]
    ok = np.isfinite(r) & np.isfinite(rp) & np.isfinite(rm) & ~pk["dead"]
    M = np.full_like(r, np.nan)
    M[ok] = r[ok] ** (n - 1) * (1 + 2 * np.exp(f[ok]) * rp[ok] * rm[ok]) / 2
    ic = int(np.nanargmin(np.where(ok, r, np.inf)))
    out = []
    for sg in (+1, -1):
        idx = np.arange(ic, len(r)) if sg > 0 else np.arange(ic, -1, -1)
        idx = idx[ok[idx]]
        rr, MM = r[idx], M[idx]
        mono = np.maximum.accumulate(rr)
        _, u = np.unique(mono, return_index=True)
        out.append(np.interp(r_pts, rr[u], MM[u], left=np.nan, right=np.nan))
    return out[0], out[1], float(r[ic]), ic, M


# --------------------------------------------------------------------------------------------------------------------
# planted data: the exact static wormhole in another double-null gauge
# --------------------------------------------------------------------------------------------------------------------
GAUGE = dict(c=0.15, b=3.0)       # F(x) = x (1 - c (1 - e^{-x/b})): F(0) = 0, F'(0) = 1, F' -> 1 - c


def gauge_F(x, c=None, b=None):
    c, b = GAUGE["c"] if c is None else c, GAUGE["b"] if b is None else b
    x = np.asarray(x, float)
    return x * (1 - c * (1 - np.exp(-x / b)))


def gauge_dF(x, c=None, b=None):
    c, b = GAUGE["c"] if c is None else c, GAUGE["b"] if b is None else b
    x = np.asarray(x, float)
    e = np.exp(-x / b)
    return 1 - c * (1 - e) - c * (x / b) * e


def gauge_rays(n, N, h):
    """Exact data of the static Ellis wormhole in the gauge x~+ = F(x+), x~- = x- of its own null coordinates (x~ = x of dn.ellis_rays): on x- = 0
    r = r_s(F(x+)/sqrt2), f = -ln F'(x+), phi = phi_s(F(x+)/sqrt2); on x+ = 0 the unchanged static data (F'(0) = 1). The whole solution is
    r = r_s(l~), l~ = (F(x+) - x-)/sqrt2, f = -ln F'(x+), phi = phi_s(l~): a coordinate change, so the throat is the curve x- = F(x+)."""
    sol, q = dn.ellis_dense(n, float(gauge_F(N * h)) / SQ2 + 2.0)
    x = np.arange(N + 1) * h
    xt = gauge_F(x)
    Sp = np.zeros((3, N + 1))
    Sm = np.zeros((3, N + 1))
    y = sol(xt / SQ2)
    Sp[0], Sp[1], Sp[2] = y[0], -np.log(gauge_dF(x)), y[2]
    y = sol(x / SQ2)
    Sm[0], Sm[1], Sm[2] = y[0], 0.0, -y[2]
    return dn.Data(n, [-1.0], N, h, Sp, Sm, name="ellis_gauge")


def gauge_locus_exact(n, t):
    """The exact locus of gauge_rays at times t (arrays): x+_* from F(x+) + x+ = sqrt2 t, x-_* = F(x+_*); l, r_min = 1, f, proper time tau = sqrt2 x-_*, v = dl/dt."""
    from scipy.optimize import brentq
    t = np.asarray(t, float)
    xp = np.array([brentq(lambda x: float(gauge_F(x)) + x - SQ2 * tt, 0.0, 2 * SQ2 * tt + 1.0, xtol=1e-14, rtol=1e-14) for tt in t])
    xm = gauge_F(xp)
    l = (xp - xm) / SQ2
    dl = np.gradient(l, t)
    return dict(l=l, f=-np.log(gauge_dF(xp)), tau=SQ2 * xm, v=dl)


def boost_locus_exact(lam, t):
    v = (1 - lam ** 2) / (1 + lam ** 2)
    return dict(l=v * t, f=np.zeros_like(t), tau=t * math.sqrt(1 - v * v), v=np.full_like(t, v))


def lam_of_v(v):
    return math.sqrt((1 - v) / (1 + v))


# --------------------------------------------------------------------------------------------------------------------
# one run of a family member: #40's events and the locus events
# --------------------------------------------------------------------------------------------------------------------
import threshold    # noqa: E402  (read-only: family(), fit_set(), bisect(), CRIT)

LEVELS = dict(col=(0.75, 0.5), exp=(2.0, 3.0))       # r_min levels of the locus events; primary: 0.5 (collapse) and 2.0 (expansion)
TCAP = 34.0                                            # the triangle covers t <= TCAP (s_max = sqrt2 TCAP / h); the two bracket ends (p*_h -/+ 6e-15) need TCAP_BR
TCAP_BR = 44.0
RAW40 = os.path.expanduser("~/Source/negative-mass-raw/I40-dh610/runs")


def make_data(n, fam, h, N, A, lam=1.0, bare=False):
    """Family member of #40 (threshold.family) on the static Ellis wormhole; lam != 1 applies the pair rescaling (scaled_rays). bare: dn.ellis_rays itself."""
    pulses, sym, eps = threshold.family(fam, A)
    if bare or lam == 1.0:
        return dn.ellis_rays(n, N, h, pulses=pulses, symmetric=sym, eps=eps)
    return scaled_rays(n, N, h, lam, pulses=pulses, symmetric=sym, eps=eps)


def record(L, with_series=False):
    """The point record: #40's event (fate, T, cell), the locus's level crossings (t, tau), the drift and the lapse at the primary crossing."""
    key = next(iter(L["ev40"]), None)
    rec = dict(fate40=(-1 if key == "collapse" else 1 if key else 0), T40=ev40_time(L, key) if key else None,
               cell40=list(L["ev40"][key]) if key else None, rows=L["rows"], nslices=len(L["s"]))
    cr = {}
    for k in L["cross"]:
        c = crossing(L, k)
        cr[k] = [float(c[0]), float(c[1])]
    rec["loc"] = cr
    kc, ke = L["cross"].get("c:0.5"), L["cross"].get("e:2.0")
    rec["fate_loc"] = -1 if kc is not None and (ke is None or kc < ke) else (1 if ke is not None else 0)
    kk = kc if rec["fate_loc"] == -1 else ke
    if kk is not None:
        rec["v_ev"], rec["f_ev"], rec["l_ev"] = float(L["v"][kk]), float(L["f"][kk]), float(L["l"][kk])
    return rec


def run_point(n, fam, h, A, smax, lam=1.0, bare=False, **kw):
    t0 = time.time()
    d = make_data(n, fam, h, smax, A, lam, bare)
    L = sweep(d, smax=smax, stop=True, col_levels=LEVELS["col"], exp_levels=LEVELS["exp"], **kw)
    locus_series(L)
    rec = record(L)
    rec.update(A=A, secs=round(time.time() - t0, 2))
    return rec, L


def tcap_smax(h, tcap=None):
    return int((TCAP if tcap is None else tcap) * SQ2 / h)


def pstar_raw40(D, fam, hk):
    """p*_h of #40's job (D, family, h index) from its raw summary: the same bisection, not repeated here."""
    n = D - 2
    h = threshold.CRIT["H"][n][hk]
    s = json.load(open(os.path.join(RAW40, f"D{D}_{fam}_h{h:g}", "summary.json")))
    assert s["h"] == h, (s["h"], h)
    return 0.5 * (s["lo"] + s["hi"]), s


def points40(D, fam, hk):
    """#40's recorded points of that job keyed by the exact amplitude."""
    h = threshold.CRIT["H"][D - 2][hk]
    pts = [json.loads(l) for l in open(os.path.join(RAW40, f"D{D}_{fam}_h{h:g}", "points.jsonl"))]
    return {p["A"]: p for p in pts}


class Sink:
    def __init__(self, path):
        self.f = open(path, "a")

    def __call__(self, rec):
        self.f.write(json.dumps(rec) + "\n")
        self.f.flush()


def job_series(D, hk, out, part=(0, 1), fam="G", es=None):
    """G (or the family fam) of #40 at D, h index hk: the series p*_h -/+ 10^-e (e = threshold.CRIT['es']) with p*_h from #40's raw (bit for bit the same
    amplitudes), each run to #40's event and to the locus events; the two bracket ends as well. Records go to out/points_<part>.jsonl."""
    n = D - 2
    h = threshold.CRIT["H"][n][hk]
    ps, s40 = pstar_raw40(D, fam, hk)
    old = points40(D, fam, hk)
    os.makedirs(out, exist_ok=True)
    sink = Sink(os.path.join(out, f"points_{part[0]}of{part[1]}.jsonl"))
    smax = tcap_smax(h)
    jobs = [("bracket", None, -1, s40["lo"]), ("bracket", None, 1, s40["hi"])]
    for e in (es or threshold.CRIT["es"]):
        for sd in (-1, 1):
            jobs.append(("series", float(e), sd, ps + sd * 10.0 ** (-e)))
    t0 = time.time()
    for k, (phase, e, sd, A) in enumerate(jobs):
        if k % part[1] != part[0]:
            continue
        sm = tcap_smax(h, TCAP_BR) if phase == "bracket" else smax
        rec, L = run_point(n, fam, h, A, sm, bare=True)
        o = old.get(A)
        rec.update(phase=phase, e=e, side=sd, h=h, D=D, family=fam, smax=sm, old_T=(o or {}).get("T"), old_fate=(o or {}).get("fate"),
                   old_cell=[(o or {}).get("j"), (o or {}).get("i")] if o else None)
        sink(rec)
        print(f"{phase} e={e} side={sd} fate40={rec['fate40']} T40={rec['T40']} (#40: {rec['old_T']}) fate_loc={rec['fate_loc']} "
              f"tau_loc={rec['loc'].get('c:0.5' if sd < 0 else 'e:2.0', [None, None])[1]} {rec['secs']} s ({time.time() - t0:.0f} s)", flush=True)
    return True


# --------------------------------------------------------------------------------------------------------------------
# the member at the stall: two runs p*_h -/+ 10^-e kept whole (locus series, profiles every time unit)
# --------------------------------------------------------------------------------------------------------------------
STALL = dict(e=9.0, t_lo=2.0, dt=1.0, r_pts=(1.25, 1.5, 2.0, 3.0, 4.0), rate_win=0.5)      # r_pts in units of the throat radius


def save_run(path, L, meta):
    arrs = {k: L[k] for k in ("s", "t", "l", "r", "f", "ic", "dk", "rp", "rm", "phi", "php", "phm", "v", "tau")}
    for k, v in L["diag"].items():
        arrs["diag_" + k] = v
    for s_k, pk in L["prof"].items():
        for kk, vv in pk.items():
            arrs[f"prof{s_k}_{kk}"] = vv
    meta = dict(meta, ev40={k: list(v) for k, v in L["ev40"].items()}, cross=L["cross"], rows=int(L["rows"]), N=int(L["N"]), h=L["h"], n=int(L["n"]),
                eps=list(map(float, L["eps"])), prof=sorted(int(s) for s in L["prof"]))
    np.savez_compressed(path, meta=json.dumps(meta), **arrs)


def load_run(path):
    z = np.load(path, allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    L = {k: z[k] for k in ("s", "t", "l", "r", "f", "ic", "dk", "rp", "rm", "phi", "php", "phm", "v", "tau")}
    L["diag"] = {k[5:]: z[k] for k in z.files if k.startswith("diag_")}
    L["prof"] = {s: {kk: z[f"prof{s}_{kk}"] for kk in ("r", "f", "phi", "rp", "rm", "php", "phm", "dead")} for s in meta["prof"]}
    L["ev40"] = {k: tuple(v) for k, v in meta["ev40"].items()}
    L["cross"], L["rows"], L["N"], L["h"], L["n"], L["eps"] = meta["cross"], meta["rows"], meta["N"], meta["h"], meta["n"], np.array(meta["eps"])
    L["meta"] = meta
    return L


def job_stall(D, hk, out, fam="G", e=None):
    n = D - 2
    h = threshold.CRIT["H"][n][hk]
    ps, s40 = pstar_raw40(D, fam, hk)
    e = STALL["e"] if e is None else e
    smax = tcap_smax(h)
    os.makedirs(out, exist_ok=True)
    tk = np.arange(STALL["t_lo"], TCAP - 3.0, STALL["dt"])
    prof_s = sorted({int(round(t * SQ2 / h)) for t in tk})
    for sd in (-1, 1):
        A = ps + sd * 10.0 ** (-e)
        t0 = time.time()
        d = make_data(n, fam, h, smax, A, bare=True)
        L = sweep(d, smax=smax, stop=True, col_levels=LEVELS["col"], exp_levels=LEVELS["exp"], prof_s=prof_s, log=lambda m: print(m, flush=True))
        locus_series(L)
        meta = dict(D=D, family=fam, h=h, A=A, side=sd, e=e, pstar=ps, smax=smax, secs=time.time() - t0, rec=record(L))
        save_run(os.path.join(out, f"stall_D{D}_h{h:g}_side{sd}.npz"), L, meta)
        print(f"stall D={D} h={h:g} side={sd}: {time.time() - t0:.0f} s, slices {len(L['s'])}, profiles {len(prof_s)}, cross {L['cross']}", flush=True)
    return True


def pair_rate(La, Lb, win=None):
    """Growth rate of the unstable mode on the locus from the two runs that bracket p*_h: lambda(tau) = d ln|r_min(b) - r_min(a)| / d tau, the local linear
    fit over +- win in tau on the slices both runs have (the same slice index s, so the same t)."""
    win = STALL["rate_win"] if win is None else win
    sa, sb = La["s"], Lb["s"]
    common = np.intersect1d(sa, sb)
    ia, ib = np.searchsorted(sa, common), np.searchsorted(sb, common)
    dr = Lb["r"][ib] - La["r"][ia]
    tau = La["tau"][ia]
    ok = np.abs(dr) > 0
    ln = np.where(ok, np.log(np.where(ok, np.abs(dr), 1.0)), np.nan)
    lam = np.full(len(common), np.nan)
    for k in range(len(common)):
        m = (np.abs(tau - tau[k]) <= win) & ok
        if m.sum() >= 5:
            lam[k] = np.polyfit(tau[m], ln[m], 1)[0]
    return common, tau, dr, lam


def charge_ratio(L, key, n, r_th, lmax=2.0):
    """(d phi_k)^2 = -2 e^f phi_+ phi_- of each field on a profile slice, at the nodes within lmax throat radii of the locus on either side: median of the
    ratio of the normal field's to the ghost's. In a static state it is q_n^2/q_g^2 and constant in l."""
    pk = L["prof"][key]
    ok = np.isfinite(pk["php"][0]) & np.isfinite(pk["phm"][0]) & np.isfinite(pk["r"]) & ~pk["dead"]
    near = ok & (pk["r"] <= (1 + lmax) * r_th)
    if L["php"].shape[1] < 2 or near.sum() < 5:
        return float("nan"), float("nan")
    g = -2 * np.exp(pk["f"][near]) * pk["php"][:, near] * pk["phm"][:, near]
    ratio = g[1] / g[0]
    return float(np.median(ratio)), float(np.max(ratio) - np.min(ratio))


def member_table(La, Lb, D, r_pts=None):
    """Per profile slice of the run La (collapse side): t, tau, r_th, drift v, lapse f, the two masses at the areal radii r_pts r_th, the fitted asymmetry s_th and
    the member's e-fold (D = 4: ellis_efold.T_of x r_th; D = 5: r_th/W, the symmetric member), the normal field's charge ratio, the pair rate lambda."""
    n = D - 2
    r_pts = np.array(STALL["r_pts"] if r_pts is None else r_pts, float)
    common, tau_p, dr, lam = pair_rate(La, Lb)
    rows = []
    for s_k in sorted(La["prof"]):
        k = np.searchsorted(La["s"], s_k)
        if k >= len(La["s"]) or La["s"][k] != s_k:
            continue
        rth = float(La["r"][k])
        if lmut("locus_diag"):          # locus_diag: the diagonal's radius instead of the locus's
            rth = float(La["diag"]["r"][s_k // 2])
        Mp, Mm, r_node, ic, M = slice_masses(La, s_k, n, r_pts * rth)
        good = np.isfinite(Mp) & np.isfinite(Mm)
        if good.sum() >= 3:
            sth, rms = fit_asym(n, rth, (r_pts * rth)[good], Mp[good], Mm[good])
        else:
            sth, rms = float("nan"), float("nan")
        kq, kq_sp = charge_ratio(La, s_k, n, rth)
        kp = np.searchsorted(common, s_k)
        rows.append(dict(s=int(s_k), t=float(La["t"][k]), tau=float(La["tau"][k]), r_th=rth, v=float(La["v"][k]), f=float(La["f"][k]), l=float(La["l"][k]),
                         r_pts=(r_pts * rth)[good].tolist(), Mp=Mp[good].tolist(), Mm=Mm[good].tolist(), sth=sth, sth_rms=rms, qratio=kq, qratio_spread=kq_sp,
                         lam=float(lam[kp]) if kp < len(common) and common[kp] == s_k else float("nan")))
    return rows


_T_CACHE = {}


def efold_ggs(g):
    """T(gamma_1) = tau_throat / r_throat of GGS09a's D = 4 member (ellis_efold.py, method A, printed W; read-only import) at gamma_1 = g."""
    if g not in _T_CACHE:
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))      # notes/: ellis_efold.py (read-only)
        import ellis_efold as ee
        lam = ee.lam_of(g)
        beta = ee.method_A(g, ee.W_A, L=30.0, lam=lam)[0]
        _T_CACHE[g] = ee.T_of(beta, g)
    return _T_CACHE[g]


# --------------------------------------------------------------------------------------------------------------------
# the expanding side: #213's observer quantities on the locus (rho, p, w in the frame of the locus's tangent)
# --------------------------------------------------------------------------------------------------------------------
import expansion as xp    # noqa: E402  (read-only: classify, slope_test, final_window, CRIT)


def locus_observer(L):
    """H = d ln r_min / d tau, rho = T(u, u), p_par = T(e', e'), p_perp = T^theta_theta and w = (p_par + n p_perp)/((n + 1) rho) on the locus, u = tangent of the locus
    (velocity v = dl/dt relative to the grid's t direction), e' the unit spatial vector orthogonal to u in the (t, l) plane. With T_+- = 0 and the node values of
    phi_+ (centred in x+) and phi_- (centred in x-). At v = 0 these are expansion.throat_stress's."""
    n, eps = int(L["n"]), L["eps"]
    r, f, v, tau = L["r"], L["f"], L["v"], L["tau"]
    ph_p, ph_m = L["php"], L["phm"]
    Tpp = np.sum(eps * ph_p ** 2, axis=1)
    Tmm = np.sum(eps * ph_m ** 2, axis=1)
    ef = np.exp(f)
    Ttt = ef * (Tpp + Tmm) / 2
    Ttl = ef * (Tpp - Tmm) / 2
    g2 = 1.0 / np.clip(1 - v ** 2, 1e-12, None)
    vv = np.zeros_like(v) if lmut("rho_diag") else v
    rho = g2 * ((1 + vv ** 2) * Ttt + 2 * vv * Ttl)
    p_par = rho
    sg = -1.0 if lmut("w_sign") else 1.0
    p_perp = sg * ef * np.sum(eps * ph_p * ph_m, axis=1)
    with np.errstate(all="ignore"):
        w = (p_par + n * p_perp) / ((n + 1) * rho)
        H = np.gradient(np.log(r), tau)
        Hf = np.exp(f / 2) * (L["rp"] + L["rm"]) / (SQ2 * r) / np.sqrt(np.clip(1 - v ** 2, 1e-12, None))
    return dict(H=H, Hf=Hf, rho=rho, ppar=p_par, pperp=p_perp, w=w)


def expansion_summary(L, W):
    """#213's declared quantities on the locus: the class on the last e-fold of r, k and tau_s of 1/H linear in tau, H/W at r/a = 10, 25, 50 (the slope test)."""
    o = locus_observer(L)
    ok = np.isfinite(L["r"]) & np.isfinite(o["H"]) & np.isfinite(o["rho"]) & np.isfinite(o["w"]) & np.isfinite(L["tau"]) & (o["H"] != 0)
    ix = np.where(ok)[0]
    cl = xp.classify(L["tau"][ix], L["r"][ix], o["H"][ix], o["rho"][ix], o["w"][ix])
    met, vals = xp.slope_test(L["r"][ix], o["H"][ix], W)
    return dict(cl=cl, slope_met=met, slope_vals=vals, obs=o, ix=ix)


def job_expand(D, hk, member, X, out, fam="G"):
    """The expanding branch of family fam (G: #40's, or S) on the locus: the rectangle run to X, member = 'tuned' (p*_h + 1e-6, p*_h from #40's raw) or an amplitude."""
    n = D - 2
    h = threshold.CRIT["H"][n][hk] if fam == "G" else xp.CRIT["H"][n][hk]
    if member == "tuned":
        ps = pstar_raw40(D, fam, hk)[0] if fam == "G" else None
        if ps is None:
            ps, _k = xp.pstar(n, h)
        A = ps + xp.CRIT["delta"]
    else:
        A = float(member)
    N = int(round(X / h))
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    d = make_data(n, fam, h, N, A, bare=True)
    L = sweep(d, rect=True, stop=False, log=lambda m: print(m, flush=True))
    locus_series(L)
    name = f"{fam}_D{D}_{member}_h{h:g}"
    meta = dict(D=D, family=fam, h=h, A=A, member=str(member), X=X, secs=time.time() - t0, rec=record(L))
    save_run(os.path.join(out, name + ".npz"), L, meta)
    print(f"expand {name}: A = {A:.12e}, {time.time() - t0:.0f} s, N = {N}, slices {len(L['s'])}, ev40 {L['ev40']}", flush=True)
    return True


# --------------------------------------------------------------------------------------------------------------------
# planted lifetime series: the boosted symmetric family S through bisection, series and fit
# --------------------------------------------------------------------------------------------------------------------
def planted_series(D, hk, v, es=None, log=None, tol=None, tcap_bis=44.0, tcap=30.0):
    """Family S in the pair-rescaled coordinates (drift v of the static throat): p*_h by bisection of the locus's fate (the first of the levels r_min = 0.5 and 2 reached;
    #40's own events are not trusted at a large drift, the diagonal is far from the throat) to `tol`, then the series p*_h -/+ 10^-e read at the locus (and, beside it,
    on #40's events). Returns the rows (T = locus tau) and the diagonal rows (T = #40's t), p*_h and the bracket."""
    n = D - 2
    h = threshold.CRIT["H"][n][hk] if isinstance(hk, int) else float(hk)       # an index of #40's h, or an explicit h (the reduced suite)
    lam = lam_of_v(v)
    smax = tcap_smax(h, tcap)
    smax_b = tcap_smax(h, tcap_bis)

    def run(p):
        rec, _L = run_point(n, "S", h, p, smax_b, lam=lam, need_loc=True)
        return dict(A=p, fate=rec["fate_loc"], T=None)
    lo, hi, stopped = threshold.bisect(run, *threshold.CRIT["bracket"]["S"], CRIT["planted_tol"] if tol is None else tol)
    ps = 0.5 * (lo + hi)
    rows, drows = [], []
    for e in (es or threshold.CRIT["es"]):
        for sd in (-1, 1):
            A = ps + sd * 10.0 ** (-e)
            rec, L = run_point(n, "S", h, A, smax, lam=lam)
            key = "c:0.5" if sd < 0 else "e:2.0"
            c = rec["loc"].get(key)
            rows.append(dict(A=A, side=sd, e=float(e), fate=rec["fate_loc"], T=c[1] if c else None))
            drows.append(dict(A=A, side=sd, e=float(e), fate=rec["fate40"], T=rec["T40"]))
            if log:
                log(f"PS e={e} side={sd} fate_loc {rec['fate_loc']} fate40 {rec['fate40']} T40 {rec['T40']} tau_loc {c and c[1]} v {rec.get('v_ev')}")
    return dict(D=D, hk=hk, v=v, lam=lam, lo=lo, hi=hi, pstar=ps, stopped=stopped, rows=rows, drows=drows)


# --------------------------------------------------------------------------------------------------------------------
# declared criteria (negative-mass-raw I225-dh610/README.md, before the runs of record) and mutations
# --------------------------------------------------------------------------------------------------------------------
CRIT.update(
    # runs of record: G of #40 at D = 4, 5 and its three h (index 0, 1, 2 of threshold.CRIT['H']); the two finest are the pair of RA
    series_hk=(0, 1, 2), RA=0.005, fm_e=11.0,
    ctl={(4, 0): 0.898291, (4, 1): 0.898043, (4, 2): 0.899089, (5, 0): 0.798948, (5, 1): 0.800044, (5, 2): 0.800277}, ctl_tol=1e-6,   # threshold.md, 6 digits
    # planted static throats (PB: pair-rescaled with drift v; PG: another double-null gauge), errors on t in [t_lo, t_c] at the finest of three h
    pb=dict(cases=((2, (0.04, 0.02, 0.01), 4.0), (3, (0.02, 0.01, 0.005), 2.7)), vs=(0.1, 0.5, -0.3), t_lo=0.8,
            bar=dict(l=2e-4, r=5e-4, f=1e-3, tau=5e-4, v=1e-3), order=(1.75, 2.25)),
    # planted lifetime series: S in the pair-rescaled coordinates; |slope W - 1| bar for the locus tau
    ps=dict(cases=((4, 0, 0.0), (4, 0, 0.1), (4, 0, 0.5), (5, 0, 0.13)), bar=1e-3),
    # the member machinery
    pm=dict(g=(0.3, 0.5), m_bar=1e-5, s_bar=1e-6, evolved=dict(D=4, h=0.01, v=0.1, t=5.0, s_bar=0.02, M_bar=2e-3)),
    po=dict(vs=(0.0, 0.3, -0.5), bar=1e-12),
    # the member at the stall (D = 4: T(gamma_1) r_th of GGS09a's family; D = 5: r_th/W of the symmetric member), compared at #40's stall time
    cmp=dict(bar=0.02, t_stall={4: 12.0, 5: 6.0}, q_bar=1e-3, win={4: (7.0, 15.0), 5: (7.0, 13.0)}),
    dr=dict(win={4: (10.0, 15.0), 5: (7.0, 13.0)}, d4_max=0.01, d5_min=0.05),     # the throat radius of the locus: a plateau at D = 4, a drifting scale at D = 5
    # the expanding side: as #213, on the locus
    X={4: 200.0, 5: 100.0}, members={4: ("tuned", 0.05, 0.2), 5: ("tuned", 0.05)}, exp_hk={4: (0, 1), 5: (0, 1)},
    kx=0.01, k_S={4: 2.0000, 5: 3.9875}, res=dict(H_rel=1e-2, k_rel=1e-2, ts_rel=1e-2),
    xc=dict(D=4, hk=0, A=0.05, X=200.0, k_bar=1e-3, ts_bar=1e-4, H_bar=1e-3),
)

MUTATIONS = {}      # name: the checks that must fire; filled from the development runs (dev_NN in the raw folder) before the runs of record

results = []


def check(cid, text, ok, detail=""):
    ok = None if ok is None else bool(ok)
    results.append((cid, ok))
    tag = "SKIP" if ok is None else ("pass" if ok else "FAIL")
    print(f"[{tag}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def fmt(xs, f="{:.2e}"):
    return "[" + ", ".join(f.format(x) for x in xs) + "]"


def order_of(errs):
    return [math.log2(errs[i] / errs[i + 1]) if errs[i + 1] > 0 and errs[i] > 0 else float("inf") for i in range(len(errs) - 1)]


# --------------------------------------------------------------------------------------------------------------------
# planted checks
# --------------------------------------------------------------------------------------------------------------------
def pb_errors(kind, n, hs, tc, par):
    """Errors of the locus (l, r_min - 1, f, tau, v) against the exact one on t in [t_lo, tc], for each h."""
    out = []
    for h in hs:
        smax = int((tc + 2) * SQ2 / h)
        d = scaled_rays(n, smax, h, lam_of_v(par)) if kind == "boost" else gauge_rays(n, smax, h)
        L = sweep(d, smax=smax, stop=False)
        locus_series(L)
        t = L["t"]
        m = (t >= CRIT["pb"]["t_lo"]) & (t <= tc)
        ex = boost_locus_exact(lam_of_v(par), t) if kind == "boost" else gauge_locus_exact(n, t)
        out.append({k: float(np.max(np.abs(a[m] - b[m]))) for k, a, b in (("l", L["l"], ex["l"]), ("r", L["r"], np.ones_like(t)), ("f", L["f"], ex["f"]),
                                                                       ("tau", L["tau"], ex["tau"]), ("v", L["v"], ex["v"]))})
    return out


def chk_PB_PG(cases=None):
    c = CRIT["pb"]
    bars, lo, hi = c["bar"], c["order"][0], c["order"][1]
    for cid, kind, pars in (("PB", "boost", c["vs"]), ("PG", "gauge", (None,))):
        ok, rows = True, []
        for n, hs, tc in (cases or c["cases"]):
            for par in pars:
                errs = pb_errors(kind, n, hs, tc, par)
                fin = errs[-1]
                good = all(fin[k] <= bars[k] for k in bars)
                ords = {k: order_of([e[k] for e in errs]) for k in (("l", "r", "f", "tau", "v") if kind == "boost" else ("l", "r", "f"))}
                good &= all(lo <= x <= hi for v in ords.values() for x in v)
                ok &= good
                rows.append(f"D = {n + 2}, {'v = ' + format(par, 'g') if kind == 'boost' else 'F(x) = x (1 - c (1 - e^(-x/b))), c = %g, b = %g' % (GAUGE['c'], GAUGE['b'])}, "
                            f"h = {hs[-1]:g}, t in [{c['t_lo']}, {tc}]: " + ", ".join(f"{k} {fin[k]:.1e}" for k in ("l", "r", "f", "tau", "v")) + "; orders " +
                            ", ".join(f"{k} {fmt(v, '{:.2f}')}" for k, v in ords.items()))
        for x in rows:
            print(f"    {cid} {x}")
        check(cid, "planted static wormhole " + ("in the pair-rescaled coordinates (drift v)" if kind == "boost" else "in another double-null gauge (f != 0, drifting throat)") +
              ": the locus's l, r_min, f, tau" + (", v" if kind == "boost" else "") + " against the exact ones, second order in h (PG: orders of l, r, f; tau and v carry the drift window's smoothing floor)", ok,
              f"bars {bars}, orders in [{lo}, {hi}]")


def chk_PM():
    c = CRIT["pm"]
    worst_m, worst_s = 0.0, 0.0
    for g in c["g"]:
        rth, sth, mp, mm = ggs_closed(g)
        a, b = asym_masses(2, rth, sth)
        worst_m = max(worst_m, abs(a - mp), abs(b - mm) / max(1.0, abs(mm)))
        rp = np.array([1.2, 1.5, 2.0, 3.0, 5.0]) * rth
        for n in (2, 3):
            Mp, Mm = member_M(n, rth, sth, rp)
            s, rms = fit_asym(n, rth, rp, Mp, Mm)
            worst_s = max(worst_s, abs(s - sth))
    print(f"    PM the ODE members against GGS09a's closed form (gamma_1 = {c['g']}): masses {worst_m:.1e}; the fit recovers s_th (n = 2, 3) to {worst_s:.1e}")
    e = c["evolved"]
    n, h = e["D"] - 2, e["h"]
    smax = int((e["t"] + 2) * SQ2 / h)
    s_k = int(round(e["t"] * SQ2 / h))
    L = sweep(scaled_rays(n, smax, h, lam_of_v(e["v"])), smax=smax, stop=False, prof_s=[s_k])
    locus_series(L)
    rth = float(L["r"][np.searchsorted(L["s"], s_k)])
    r_pts = np.array(STALL["r_pts"]) * rth
    Mp, Mm, rn, ic, M = slice_masses(L, s_k, n, r_pts)
    ex = 1.0 / (2 * r_pts)
    s, rms = fit_asym(n, rth, r_pts, Mp, Mm)
    dM = max(np.max(np.abs(Mp - ex)), np.max(np.abs(Mm - ex)))
    print(f"    PM the evolved static wormhole (D = {e['D']}, h = {h}, v = {e['v']}, t = {e['t']}): r_th {rth:.6f}; M on the two sides against 1/(2r): {dM:.1e}; fitted s_th {s:+.2e}")
    check("PM", "member machinery: the ODE members reproduce GGS09a's closed form, the fit recovers s_th from M(r) on both sides, and an evolved symmetric static wormhole gives "
          "M = 1/(2r) and s_th = 0", worst_m <= c["m_bar"] and worst_s <= c["s_bar"] and abs(s) <= e["s_bar"] and dM <= e["M_bar"],
          f"bars: masses {c['m_bar']:g}, s_th {c['s_bar']:g}, evolved |s_th| {e['s_bar']}, M {e['M_bar']:g}")


def chk_PO():
    """The locus-frame rho and w on planted fields: a single null stream (phi_- = 0) and a two-sided state, exact algebra for a drift v."""
    c = CRIT["po"]
    worst = 0.0
    for v in c["vs"]:
        for n in (2, 3):
            for (a, b) in ((0.7, 0.0), (0.4, -0.9)):
                f = 0.3
                L = dict(n=n, eps=np.array([-1.0]), r=np.full(3, 2.0), f=np.full(3, f), v=np.full(3, v), tau=np.array([1.0, 1.1, 1.2]), php=np.full((3, 1), a), phm=np.full((3, 1), b),
                         rp=np.zeros(3), rm=np.zeros(3))
                o = locus_observer(L)
                o = {k: x[1:2] for k, x in o.items()}
                g2 = 1 / (1 - v * v)
                Tpp, Tmm = -a * a, -b * b
                rho_ex = g2 * math.exp(f) * ((1 + v * v) * (Tpp + Tmm) / 2 + v * (Tpp - Tmm))
                pperp = -math.exp(f) * a * b
                w_ex = (rho_ex + n * pperp) / ((n + 1) * rho_ex)
                worst = max(worst, abs(o["rho"][0] - rho_ex) / abs(rho_ex), abs(o["w"][0] - w_ex))
    check("PO", "locus frame: rho, p_par = rho and w on planted fields (one null stream and a two-sided state) at drifts v = 0, 0.3, -0.5, n = 2, 3", worst <= c["bar"],
          f"max deviation {worst:.1e}, bar {c['bar']:g}")


def chk_PL():
    ok, rows = True, []
    for name, cls_exp, sl_exp, (tau, r, H, rho, w) in xp.planted_series():
        out = xp.classify(tau, r, H, rho, w)
        met, vals = xp.slope_test(r, H, W_EXACT[2])
        ok &= out["class"] == cls_exp and met == sl_exp
        rows.append(f"{name}: class {out['class']} (expected {cls_exp}), slope test {'met' if met else 'not met'} (expected {'met' if sl_exp else 'not met'})")
    for x in rows:
        print("    PL " + x)
    check("PL", "planted series with known classes (#213's five) through the classification used on the locus", ok)


def chk_PS(P):
    """P: list of planted_series results (stored or fresh)."""
    c = CRIT["ps"]
    ok, rows = True, []
    for r in P:
        D, W = r["D"], W_EXACT[r["D"] - 2]
        rr = [x for x in r["rows"] if x["T"] is not None and x["fate"] == x["side"]]
        dd = [x for x in r["drows"] if x["T"] is not None and x["fate"] == x["side"]]
        fits, prim, env = threshold.fit_set(rr, pstar_used(r["pstar"]))
        fd, pd_, ed = threshold.fit_set(dd, pstar_used(r["pstar"]))
        good = abs(prim["slope"] * W - 1) <= c["bar"] and r["stopped"] is None
        ok &= good
        rows.append(f"D = {D}, h index {r['hk']}, drift v = {r['v']:g} (lam = {r['lam']:.6f}): p*_h = {r['pstar']:.12e}, bracket {r['hi'] - r['lo']:.1e}; locus tau slope {prim['slope']:.6f} +- {prim['se']:.6f} "
                    f"(OLS), rel to 1/W {prim['slope'] * W - 1:+.2e}, spread [{env[0]:.6f}, {env[1]:.6f}]; diagonal t (#40's events) {pd_['slope']:.6f} +- {pd_['se']:.6f}, rel {pd_['slope'] * W - 1:+.2e}, "
                    f"spread [{ed[0]:.6f}, {ed[1]:.6f}]")
    for x in rows:
        print("    PS " + x)
    check("PS", "planted lifetime series (family S in the pair-rescaled coordinates, throat drifting at v): the slope of the locus's tau against ln(1/|p - p*_h|) is 1/W",
          ok, f"|slope W - 1| <= {c['bar']:g} in every case")


# --------------------------------------------------------------------------------------------------------------------
# the series of G at the locus: control, fates, slopes, resolution
# --------------------------------------------------------------------------------------------------------------------
def pstar_used(ps):
    return ps + 1e-9 if lmut("pstar_off") else ps


def series_rows(recs, key_fn):
    """Rows (side, e, fate, A, T) of the series records for threshold.fit_set; T from key_fn(rec) (None: no event)."""
    rows = []
    for r in recs:
        if r["phase"] != "series":
            continue
        T = key_fn(r)
        rows.append(dict(A=r["A"], side=r["side"], e=r["e"], fate=r["fate40"], T=T))
    return [x for x in rows if x["T"] is not None and x["fate"] == x["side"]]


def fit_summary(rows, ps, W):
    fits, prim, env = threshold.fit_set(rows, pstar_used(ps))
    w0 = threshold.CRIT["windows"][0]
    sc = [f for f in fits if f["window"] == w0 and f["sides"] == (-1,)][0]
    ex = [f for f in fits if f["window"] == w0 and f["sides"] == (1,)][0]
    return dict(prim=prim, env=env, col=sc["slope"], exp=ex["slope"], rel=prim["slope"] * W - 1, inside=threshold.verdict_inside(env, 1.0 / W))


def analyze_series(S, full=True):
    """S: {(D, hk): records}. Prints each job's table; checks JOB, CTL, FA, ED per job and RA per D. Returns {(D, hk): summary}."""
    out = {}
    for (D, hk), recs in sorted(S.items()):
        n = D - 2
        h = threshold.CRIT["H"][n][hk]
        W = W_EXACT[n]
        ps, _s40 = pstar_raw40(D, "G", hk)
        if lmut("drop_bracket"):         # drop_bracket: the bracket ends are not read
            recs = [r for r in recs if r["phase"] != "bracket"]
        ser = [r for r in recs if r["phase"] == "series"]
        brk = [r for r in recs if r["phase"] == "bracket"]
        nexp = (len(threshold.CRIT["es"]) * 2, 2)
        tag = f"D={D} G h={h:g}"
        if full:
            check(f"JOB[{tag}]", "all series points and both bracket ends are present and each ran to #40's event", len(ser) == nexp[0] and len(brk) == nexp[1]
                  and all(r["fate40"] != 0 for r in recs), f"{len(ser)} of {nexp[0]} series, {len(brk)} of {nexp[1]} bracket")
        # control: the same runs give #40's events bit for bit
        both = [r for r in recs if r.get("old_T") is not None]
        bad = [r for r in both if r["T40"] != r["old_T"] or r["fate40"] != r["old_fate"] or (r["cell40"] != r["old_cell"] and r["old_cell"] is not None)]
        rows40 = series_rows(recs, lambda r: r["T40"])
        s40 = fit_summary(rows40, ps, W)
        ref = CRIT["ctl"].get((D, hk))
        okc = not bad and (not full or ref is None or abs(s40["prim"]["slope"] - ref) <= CRIT["ctl_tol"])
        check(f"CTL[{tag}]", "control: the same amplitudes run through the locus code give #40's event, cell and time bit for bit, and #40's diagonal slope" + (f" {ref}" if ref else ""), okc,
              f"{len(both) - len(bad)} of {len(both)} points equal; diagonal slope {s40['prim']['slope']:.6f} +- {s40['prim']['se']:.6f} (rel to 1/W {s40['rel']:+.4f})")
        fm = [r for r in ser if r["e"] <= CRIT["fm_e"]]
        offs = [r for r in fm if not (r["fate_loc"] == r["fate40"] == r["side"])]
        check(f"FA[{tag}]", f"every series point with e <= {CRIT['fm_e']:g} has the same fate on #40's events and on the locus's levels, on its own side of p*_h", not offs,
              f"{len(offs)} of {len(fm)} points off" + (": " + ", ".join(f"e={r['e']} side {r['side']}: {r['fate40']}/{r['fate_loc']}" for r in offs[:4]) if offs else ""))
        key_p = lambda r: (r["loc"].get("c:0.5" if r["side"] < 0 else "e:2.0") or [None, None])[1]       # noqa: E731
        key_t = lambda r: (r["loc"].get("c:0.5" if r["side"] < 0 else "e:2.0") or [None, None])[0]       # noqa: E731
        key_a = lambda r: (r["loc"].get("c:0.75" if r["side"] < 0 else "e:3.0") or [None, None])[1]      # noqa: E731
        sp, st, sa = (fit_summary(series_rows(recs, k), ps, W) for k in (key_p, key_t, key_t if lmut("ed_t") else key_a))     # ed_t: the alternative read in t
        vs = [r["v_ev"] for r in ser if r.get("v_ev") is not None]
        print(f"    SL {tag}: p*_h = {ps:.12e}; locus tau slope {sp['prim']['slope']:.6f} +- {sp['prim']['se']:.6f} (OLS SE), n = {sp['prim']['n']}; collapse {sp['col']:.6f}, expansion {sp['exp']:.6f}; "
              f"spread [{sp['env'][0]:.6f}, {sp['env'][1]:.6f}]; 1/W = {1 / W:.6f}, rel {sp['rel']:+.2e}")
        print(f"       locus t slope {st['prim']['slope']:.6f} +- {st['prim']['se']:.6f} (rel {st['rel']:+.2e}); diagonal t of #40 {s40['prim']['slope']:.6f} +- {s40['prim']['se']:.6f} "
              f"(rel {s40['rel']:+.2e}, spread [{s40['env'][0]:.6f}, {s40['env'][1]:.6f}]); alternative levels (0.75, 3) tau {sa['prim']['slope']:.6f}; drift at the events {fmt([min(vs), max(vs)], '{:.4f}') if vs else 'n/a'}")
        check(f"ED[{tag}]", "the locus tau slope does not depend on the event level: the alternative levels (r_min = 0.75, 3) give a slope inside the primary fit spread",
              sp["env"][0] <= sa["prim"]["slope"] <= sp["env"][1], f"alternative {sa['prim']['slope']:.6f}, primary spread [{sp['env'][0]:.6f}, {sp['env'][1]:.6f}]")
        out[(D, hk)] = dict(sp=sp, st=st, s40=s40, sa=sa, ps=ps, h=h, W=W)
    for D in sorted({k[0] for k in out}):
        ks = sorted(k for k in out if k[0] == D)
        if len(ks) >= 2:
            a, b = out[ks[-2]]["sp"]["prim"]["slope"], out[ks[-1]]["sp"]["prim"]["slope"] * (1.01 if lmut("ra_scale") else 1.0)     # ra_scale: the finer h's slope scaled by 1.01
            check(f"RA[D={D} G]", f"the locus tau slopes at the two finest h agree within {CRIT['RA']:.1%}", abs(a / b - 1) <= CRIT["RA"],
                  f"{a:.6f} (h = {out[ks[-2]]['h']:g}), {b:.6f} (h = {out[ks[-1]]['h']:g}); rel {a / b - 1:+.2e}")
    return out


# --------------------------------------------------------------------------------------------------------------------
# the member at the stall
# --------------------------------------------------------------------------------------------------------------------
def member_report(D, hk, La, Lb, sp=None):
    """Prints the member table of a stall pair and checks CMP: at #40's stall time the pair rate against the member's e-fold (D = 4: T(gamma_1) r_th; D = 5: r_th/W of the symmetric
    member) and, if the series' slope is given, that slope against the member's e-fold."""
    n = D - 2
    h = threshold.CRIT["H"][n][hk]
    W = W_EXACT[n]
    rows = member_table(La, Lb, D)
    tag = f"D={D} G h={h:g}"
    print(f"    MB {tag}: p*_h -/+ 10^-{STALL['e']:g}; per profile slice: t, tau, r_th, drift v, lapse f, M(r) on the + and - side at r = {STALL['r_pts']} r_th, fitted s_th = r_th N'/N "
          f"(gamma_1 = s_th/sqrt(1 - s_th^2) in D = 4), normal/ghost charge ratio q_n^2/q_g^2, pair rate lambda (1/lambda), member e-fold")
    for r in rows:
        if r["t"] < 3.5 or not np.isfinite(r["sth"]):
            continue
        s = r["sth"]
        if D == 4:
            g1 = s / math.sqrt(max(1e-12, 1 - s * s))
            mem = efold_ggs(round(g1, 3)) * r["r_th"]
        else:
            g1, mem = float("nan"), r["r_th"] / W
        r["gamma1"], r["member"] = g1, mem
        print(f"       t {r['t']:5.1f} tau {r['tau']:6.2f} r_th {r['r_th']:.4f} v {r['v']:+.4f} f {r['f']:+.4f} s_th {s:+.4f} (rms {r['sth_rms']:.1e}) gamma_1 {g1:.4f} "
              f"q ratio {r['qratio']:+.1e} lambda {r['lam']:.4f} (1/lambda {1 / r['lam'] if r['lam'] > 0 else float('nan'):.4f}) member {mem:.4f} (1/lambda / member {1 / r['lam'] / mem:.4f}); "
              f"M(1.5, 2, 3 r_th) + {fmt(r['Mp'][1:4], '{:.4f}')} - {fmt(r['Mm'][1:4], '{:.4f}')}")
    ts = CRIT["cmp"]["t_stall"][D]
    cand = [r for r in rows if "member" in r]
    if not cand:
        check(f"CMP1[{tag}]", "the stall member's e-fold against the pair rate", False, "no profile slice")
        return rows
    r0 = min(cand, key=lambda r: abs(r["t"] - ts))
    mem_txt = "T(gamma_1) r_th, GGS09a" if D == 4 else "r_th/W of the symmetric member; T(gamma_1) is not available at D = 5"
    d1 = 1 / r0["lam"] / r0["member"] - 1
    check(f"CMP1[{tag}]", f"the member at #40's stall time (t = {ts}): its e-fold ({mem_txt}) against the growth rate of the unstable mode from the two runs that bracket p*_h, within {CRIT['cmp']['bar']:.0%}",
          abs(d1) <= CRIT["cmp"]["bar"], f"t = {r0['t']:.1f}: 1/lambda {1 / r0['lam']:.4f}, member {r0['member']:.4f} (r_th {r0['r_th']:.4f}, gamma_1 {r0['gamma1']:.3f}), rel {d1:+.3f}")
    if sp is not None:
        d2 = sp["prim"]["slope"] / r0["member"] - 1
        check(f"CMP2[{tag}]", f"the member at #40's stall time (t = {ts}): its e-fold ({mem_txt}) against the locus tau slope of the series, within {CRIT['cmp']['bar']:.0%}",
              abs(d2) <= CRIT["cmp"]["bar"], f"series slope {sp['prim']['slope']:.4f}, member {r0['member']:.4f}, rel {d2:+.3f}; the member's e-fold over t in [{ts - 4:g}, {ts + 8:g}]: "
              + fmt([r["member"] for r in cand if ts - 4 <= r["t"] <= ts + 8], "{:.3f}"))
    wa, wb = CRIT["cmp"]["win"][D]
    inw = [r for r in cand if wa <= r["t"] <= wb]
    dev = [1 / r["lam"] / r["member"] - 1 for r in inw]
    check(f"CMP3[{tag}]", f"at every profile slice with t in [{wa:g}, {wb:g}] the growth rate of the unstable mode is the member's ({mem_txt}) within {CRIT['cmp']['bar']:.0%}: the member follows the drifting throat",
          bool(inw) and max(abs(x) for x in dev) <= CRIT["cmp"]["bar"], f"{len(inw)} slices; 1/lambda / member - 1 in " + (f"[{min(dev):+.4f}, {max(dev):+.4f}]" if dev else "n/a"))
    wa, wb = CRIT["dr"]["win"][D]
    ra = [r for r in rows if abs(r["t"] - wa) < 0.5]
    rb = [r for r in rows if abs(r["t"] - wb) < 0.5]
    if ra and rb:
        dr_ = rb[0]["r_th"] / ra[0]["r_th"] - 1
        ok = abs(dr_) <= CRIT["dr"]["d4_max"] if D == 4 else dr_ >= CRIT["dr"]["d5_min"]
        check(f"DR[{tag}]", f"the throat radius of the locus from t = {wa:g} to {wb:g}: " + (f"a plateau, |change| <= {CRIT['dr']['d4_max']:.0%}" if D == 4 else f"a drifting scale, growth >= {CRIT['dr']['d5_min']:.0%}"),
              ok, f"r_th {ra[0]['r_th']:.4f} -> {rb[0]['r_th']:.4f}, change {dr_:+.4f}")
    qs = [r["qratio"] for r in rows if ts - 1 <= r["t"] <= ts + 3 and np.isfinite(r["qratio"])]
    print(f"       the normal field's charge ratio q_n^2/q_g^2 over t in [{ts - 1:g}, {ts + 3:g}]: " + (fmt(qs, '{:+.1e}') if qs else "n/a") + "; reported, not gated")
    return rows
    r0 = min(cand, key=lambda r: abs(r["t"] - ts))
    d1 = 1 / r0["lam"] / r0["member"] - 1
    ok = abs(d1) <= CRIT["cmp"]["bar"]
    det = f"t = {r0['t']:.1f}: 1/lambda {1 / r0['lam']:.4f}, member {r0['member']:.4f} (r_th {r0['r_th']:.4f}, gamma_1 {r0['gamma1']:.3f}), rel {d1:+.3f}"
    if sp is not None:
        d2 = sp["prim"]["slope"] / r0["member"] - 1
        ok &= abs(d2) <= CRIT["cmp"]["bar"]
        det += f"; series slope {sp['prim']['slope']:.4f} against the member: rel {d2:+.3f}"
    check(f"CMP[{tag}]", f"the member at #40's stall time (t = {ts}): its e-fold ({'T(gamma_1) r_th, GGS09a' if D == 4 else 'r_th/W of the symmetric member, T(gamma_1) not available at D = 5'}) "
          f"against the pair rate" + (" and against the series slope" if sp is not None else "") + f" within {CRIT['cmp']['bar']:.0%}", ok, det)
    qs = [abs(r["qratio"]) for r in rows if 5.0 <= r["t"] <= ts + 2 and np.isfinite(r["qratio"])]
    print(f"       the normal field's charge ratio |q_n^2/q_g^2| over t in [5, {ts + 2:g}]: max {max(qs) if qs else float('nan'):.1e} "
          f"({'below' if qs and max(qs) < CRIT['cmp']['q_bar'] else 'not below'} {CRIT['cmp']['q_bar']:g}; reported, not gated)")
    return rows


# --------------------------------------------------------------------------------------------------------------------
# the expanding side on the locus
# --------------------------------------------------------------------------------------------------------------------
def diag_H(L):
    """The diagonal observer of #213 on the same run: tau, r, H of the diagonal cells (dn.finish_diag's definitions)."""
    h = L["h"]
    d = L["diag"]
    ok = np.isfinite(d["r"])
    f, r, rp, rm = d["f"][ok], d["r"][ok], d["rp"][ok], d["rm"][ok]
    j = np.where(ok)[0]
    dtau = SQ2 * h * np.exp(-d["f"] / 2)
    tau = np.nancumsum(dtau) - 0.5 * np.nan_to_num(dtau)
    return tau[ok], r, np.exp(f / 2) * (rp + rm) / (SQ2 * r), j


def expand_lines(E, W, n):
    """The printed table for one expanding run: the locus at r/a = 1.5, 3, 10, 25, 50, 100."""
    L, S = E["L"], E["S"]
    o, ix = S["obs"], S["ix"]
    t, tau, r, l, v = L["t"][ix], L["tau"][ix], L["r"][ix], L["l"][ix], L["v"][ix]
    H, rho, w = o["H"][ix], o["rho"][ix], o["w"][ix]
    lines = [f"    {'r/a':>6s} {'t':>8s} {'tau':>9s} {'l':>8s} {'v':>8s} {'H':>9s} {'H/W':>7s} {'rho':>11s} {'|rho|/H^2':>10s} {'w_bar':>8s}"]
    for r0 in (1.5, 3.0, 10.0, 25.0, 50.0, 100.0):
        jj = np.where(r >= r0)[0]
        if len(jj):
            j = jj[0]
            lines.append(f"    {r[j]:6.2f} {t[j]:8.2f} {tau[j]:9.4f} {l[j]:+8.3f} {v[j]:+8.4f} {H[j]:9.4f} {H[j] / W:7.3f} {rho[j]:+11.3e} {abs(rho[j]) / H[j] ** 2:10.4f} {w[j]:8.5f}")
    return lines


def analyze_expand(E, ctrl213=None, full=True):
    """E: list of dict(D, fam, member, h, L, meta). Prints each run's class, k, tau_s, H/W; checks KX (G's k against S's of #213), RES (two h), XC (S through the locus
    code against #213's diagonal run)."""
    summ = {}
    for e in E:
        n = e["D"] - 2
        W = W_EXACT[n]
        e["S"] = expansion_summary(e["L"], W)
    for e in sorted(E, key=lambda e: (e["fam"], e["D"], str(e["member"]), -e["h"])):
        n, W = e["D"] - 2, W_EXACT[e["D"] - 2]
        S, L = e["S"], e["L"]
        cl = S["cl"]
        tag = f"{e['fam']} D={e['D']} {e['member']} h={e['h']:g}"
        print(f"  --- {tag} (A = {e['meta']['A']:.10e}, X = {e['meta']['X']:g}, {e['meta'].get('secs', float('nan')):.0f} s) ---")
        for x in expand_lines(e, W, n):
            print(x)
        ix = S["ix"]
        print(f"    final window (last e-fold of r, {cl['npts']} points, r/a = {cl['r0']:.1f} to {cl['r1']:.1f}): H max/min {cl['Hratio']:.3f}, H end/start {cl['Hgrow']:.3f}, |rho| max/min "
              f"{cl['rhoratio']:.3f}, w_bar in [{cl['wmin']:.4f}, {cl['wmax']:.4f}]; 1/H linear in tau: k = {cl['k']:.4f}, tau_s = {cl['tau_s']:.4f} (last point {cl['tau_end']:.4f}), max residual {cl['ft_resid']:.1e} -> class {cl['class']}")
        met, vals = S["slope_met"], S["slope_vals"]
        print(f"    slope test (#211): H/W at r/a = " + ", ".join(f"{k:g}: {v:.3f}" for k, v in vals.items()) + f"; bar {xp.CRIT['SL']['bar']:.0%}: {'met' if met else 'not met'}")
        v = L["v"][ix]
        print(f"    locus drift v = dl/dt: {np.nanmin(v[len(v) // 4:]):+.4f} to {np.nanmax(v[len(v) // 4:]):+.4f} over the last three quarters; l at the end {L['l'][ix][-1]:+.3f}; "
              f"max |H / H_fields - 1| over the last half {np.nanmax(np.abs(S['obs']['H'][ix][len(ix) // 2:] / S['obs']['Hf'][ix][len(ix) // 2:] - 1)):.1e}")
        try:
            tauD, rD, HD, jD = diag_H(L)
            jj = np.where(np.isfinite(HD) & (HD > 0))[0]
            vd = {r0: float(HD[jj[np.where(rD[jj] >= r0)[0][0]]] / W) for r0 in (10.0, 25.0, 50.0) if (rD[jj] >= r0).any()}
            fw = xp.final_window(tauD[jj], rD[jj], 1.0)
            b, a = np.polyfit(tauD[jj][fw], 1 / HD[jj][fw], 1)
            print(f"    the diagonal of the same run (#213's observer): H/W at r/a = " + ", ".join(f"{k:g}: {v_:.3f}" for k, v_ in vd.items()) + f"; k = {-1 / b:.4f}, tau_s = {-a / b:.4f}")
        except Exception as exc:       # a report
            print(f"    the diagonal's report is not available: {exc!r}")
        summ[(e["fam"], e["D"], str(e["member"]), e["h"])] = cl
    if E:
        okk, rows = True, []
        for e in E:
            if e["fam"] != "G":
                continue
            cl, kS = e["S"]["cl"], CRIT["k_S"][e["D"]]
            good = abs(cl["k"] / kS - 1) <= CRIT["kx"] and cl["class"] in ("FT",)
            okk &= good
            rows.append(f"D{e['D']} {e['member']} h={e['h']:g}: class {cl['class']}, k {cl['k']:.4f} against S's {kS} ({cl['k'] / kS - 1:+.1e})")
        for x in rows:
            print("    KX " + x)
        if rows:
            check("KX", f"G on the locus: class FT and k within {CRIT['kx']:.0%} of the symmetric throat's (#213: D = 4: {CRIT['k_S'][4]}, D = 5: {CRIT['k_S'][5]})", okk)
        groups = {}
        for e in E:
            groups.setdefault((e["fam"], e["D"], str(e["member"])), []).append(e)
        allok, lines = True, []
        for key, es in sorted(groups.items()):
            es = sorted(es, key=lambda e: -e["h"])
            for a, b in zip(es[:-1], es[1:]):
                ca, cb = a["S"]["cl"], b["S"]["cl"]
                dH = []
                for r0 in (10.0, 25.0, 50.0):
                    ga, gb = a["S"], b["S"]
                    ra, rb = a["L"]["r"][ga["ix"]], b["L"]["r"][gb["ix"]]
                    if (ra >= r0).any() and (rb >= r0).any():
                        dH.append(abs(np.interp(r0, ra, ga["obs"]["H"][ga["ix"]]) / np.interp(r0, rb, gb["obs"]["H"][gb["ix"]]) - 1))
                dk, dts = abs(ca["k"] / cb["k"] - 1), abs(ca["tau_s"] / cb["tau_s"] - 1)
                ok = ca["class"] == cb["class"] and (max(dH) if dH else 0) <= CRIT["res"]["H_rel"] and dk <= CRIT["res"]["k_rel"] and dts <= CRIT["res"]["ts_rel"]
                allok &= ok
                lines.append(f"{key[0]} D{key[1]} {key[2]}: h {a['h']:g} vs {b['h']:g}: class {ca['class']} / {cb['class']}; tau_s {ca['tau_s']:.4f} / {cb['tau_s']:.4f}; k {ca['k']:.4f} / {cb['k']:.4f}; H at fixed r/a {fmt(dH)}")
        for x in lines:
            print("    RES " + x)
        if lines:
            check("RES", f"two resolutions agree: same class, H at fixed r/a within {CRIT['res']['H_rel']:.0%}, k and tau_s within {CRIT['res']['k_rel']:.0%}", allok)
    if ctrl213 is not None:
        c = CRIT["xc"]
        e = [x for x in E if x["fam"] == "S" and x["D"] == c["D"] and str(x["member"]) == str(c["A"])]
        if e:
            cl = e[0]["S"]["cl"]
            dk, dts = abs(cl["k"] - ctrl213["k"]), abs(cl["tau_s"] * (1.01 if lmut("xc_tau") else 1.0) / ctrl213["tau_s"] - 1)     # xc_tau: tau_s scaled by 1.01
            dH = max(abs(e[0]["S"]["slope_vals"][r0] / ctrl213["vals"][r0] - 1) for r0 in (10.0, 25.0, 50.0))
            check("XC", "control: family S (A = 0.05) through the locus code gives #213's diagonal run: class, k, tau_s and H/W at r/a = 10, 25, 50",
                  cl["class"] == ctrl213["class"] and dk <= c["k_bar"] and dts <= c["ts_bar"] and dH <= c["H_bar"],
                  f"class {cl['class']} / {ctrl213['class']}; k {cl['k']:.5f} / {ctrl213['k']:.5f}; tau_s {cl['tau_s']:.5f} / {ctrl213['tau_s']:.5f}; H/W rel dev {dH:.1e}")
        else:
            check("XC", "control: S through the locus code against #213", False, "the S run is missing")
    return summ


def ctrl213_values():
    """#213's classification of its stored S run (D = 4, A = 0.05, h = 0.02), recomputed with expansion.py's own functions from its raw file."""
    c = CRIT["xc"]
    path = os.path.expanduser(f"~/Source/negative-mass-raw/I213-dh610/runs/D{c['D']}_A{c['A']:g}_h{xp.CRIT['H'][c['D'] - 2][c['hk']]:g}.npz")
    R = xp.load_run(path)
    g = xp.diag_arrays(R)
    cl = xp.classify(g["tau"], g["r"], g["H"], g["rho"], g["w"])
    met, vals = xp.slope_test(g["r"], g["H"], W_EXACT[c["D"] - 2])
    return dict(class_=cl["class"], k=cl["k"], tau_s=cl["tau_s"], vals=vals, **{"class": cl["class"]})


# --------------------------------------------------------------------------------------------------------------------
# suites: the runs of record (from the raw folder) and the reduced suite for the mutation runs
# --------------------------------------------------------------------------------------------------------------------
def job_planted(D, hk, v, out):
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    P = planted_series(D, hk, v, log=lambda m: print(m, flush=True))
    P["secs"] = time.time() - t0
    json.dump(P, open(os.path.join(out, f"ps_D{D}_hk{hk}_v{v:g}.json"), "w"), indent=1)
    print(f"planted series D = {D} hk = {hk} v = {v}: {P['secs']:.0f} s, p*_h = {P['pstar']:.12e}, bracket {P['hi'] - P['lo']:.1e}", flush=True)
    return True


def load_record(raw):
    import glob
    S = {}
    for D in (4, 5):
        for hk in CRIT["series_hk"]:
            h = threshold.CRIT["H"][D - 2][hk]
            recs = []
            for f in sorted(glob.glob(os.path.join(raw, "series", f"D{D}_G_h{h:g}", "points_*.jsonl"))):
                recs += [json.loads(l) for l in open(f)]
            if recs:
                S[(D, hk)] = recs
    stalls = {}
    for D in (4, 5):
        for hk in CRIT["series_hk"]:
            h = threshold.CRIT["H"][D - 2][hk]
            a, b = (os.path.join(raw, "stall", f"stall_D{D}_h{h:g}_side{sd}.npz") for sd in (-1, 1))
            if os.path.exists(a) and os.path.exists(b):
                stalls[(D, hk)] = (load_run(a), load_run(b))
    E = []
    for f in sorted(glob.glob(os.path.join(raw, "expand", "*.npz"))):
        L = load_run(f)
        m = L["meta"]
        E.append(dict(D=m["D"], fam=m["family"], member=m["member"], h=m["h"], L=L, meta=m))
    P = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(raw, "planted", "ps_*.json")))]
    return dict(series=S, stall=stalls, expand=E, ps=P, mode="record")


def mini_suite():
    """The suite at reduced size for the mutation runs (D = 4): the planted checks at one case each, the G series at h = 0.02 for four e (4 x 2 points and the bracket), the
    stall pair at h = 0.02, the expanding branch of S and G (A = 0.05) at X = 40, h = 0.02 and 0.01, and the planted lifetime series at h = 0.02 (v = 0.1, nine e). A failure to build a part
    is recorded and reported as the check BUILD."""
    S, stalls, E, P, errs = {}, {}, [], [], []
    try:
        ps, s40 = pstar_raw40(4, "G", 0)
        h = threshold.CRIT["H"][2][0]
        old = points40(4, "G", 0)
        recs = []
        smax = tcap_smax(h)
        for phase, e, sd, A in [("bracket", None, -1, s40["lo"]), ("bracket", None, 1, s40["hi"])] + [("series", float(e), sd, ps + sd * 10.0 ** (-e)) for e in (4.0, 6.5, 9.0, 11.0) for sd in (-1, 1)]:
            rec, _L = run_point(2, "G", h, A, tcap_smax(h, TCAP_BR) if phase == "bracket" else smax, bare=True)
            o = old.get(A)
            rec.update(phase=phase, e=e, side=sd, h=h, D=4, family="G", smax=smax, old_T=(o or {}).get("T"), old_fate=(o or {}).get("fate"),
                       old_cell=[(o or {}).get("j"), (o or {}).get("i")] if o else None)
            recs.append(rec)
        S[(4, 0)] = recs
    except Exception as exc:
        errs.append(f"series: {exc!r}")
    try:
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            job_stall(4, 0, td)
            stalls[(4, 0)] = (load_run(os.path.join(td, f"stall_D4_h{threshold.CRIT['H'][2][0]:g}_side-1.npz")), load_run(os.path.join(td, f"stall_D4_h{threshold.CRIT['H'][2][0]:g}_side1.npz")))
    except Exception as exc:
        errs.append(f"stall: {exc!r}")
    try:
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            for fam, member in (("S", "0.05"), ("G", "0.05"), ("G", "tuned")):
                for hk in (0, 1):
                    job_expand(4, hk, member, 40.0, td, fam=fam)
            import glob
            for f in sorted(glob.glob(os.path.join(td, "*.npz"))):
                L = load_run(f)
                m = L["meta"]
                E.append(dict(D=m["D"], fam=m["family"], member=m["member"], h=m["h"], L=L, meta=m))
    except Exception as exc:
        errs.append(f"expand: {exc!r}")
    try:
        P.append(planted_series(4, 0, 0.1, es=[3.0 + k for k in range(9)], tol=1e-11))
    except Exception as exc:
        errs.append(f"planted series: {exc!r}")
    return dict(series=S, stall=stalls, expand=E, ps=P, mode="mini", build_error="; ".join(errs) or None)


def analyze_suite(suite):
    rec = suite["mode"] == "record"
    print(f"throat_locus.py analyze: mode {suite['mode']}, mutation {LMUT or dn.MUT}, numba {dn.nb.__version__ if dn.HAVE_NUMBA else 'absent'}, numpy {np.__version__}, {time.strftime('%Y-%m-%d %H:%M:%S')}")
    if suite.get("build_error"):
        check("BUILD", "the runs of the suite complete", False, suite["build_error"])
    for fn, kw in ((chk_PB_PG, dict(cases=None if rec else [(2, (0.04, 0.02, 0.01), 4.0)])), (chk_PM, {}), (chk_PO, {}), (chk_PL, {})):
        try:
            fn(**kw)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            check(fn.__name__[4:], "check raised an exception", False, repr(exc))
    if suite["ps"]:
        chk_PS(suite["ps"])
    elif rec:
        check("PS", "planted lifetime series", False, "no stored planted series")
    sumr = analyze_series(suite["series"], full=rec)
    for (D, hk), (La, Lb) in sorted(suite["stall"].items()):
        member_report(D, hk, La, Lb, sp=sumr.get((D, hk), {}).get("sp"))
    if suite["expand"]:
        analyze_expand(suite["expand"], ctrl213=ctrl213_values() if rec else None)
    bad = [c for c, ok in results if ok is False]
    print(f"\n{sum(1 for _, ok in results if ok)} of {len(results)} checks pass; {len(bad)} fail" + (f": {', '.join(bad)}" if bad else ""))
    return bad


def main():
    global LMUT
    argv = sys.argv[1:]
    if not argv:
        sys.exit(__doc__)
    cmd, rest = argv[0], argv[1:]
    a, mut = {}, None
    for x in rest:
        if x.startswith("--mutate="):
            mut = x.split("=", 1)[1]
        elif x.startswith("--") and "=" in x:
            k, v = x[2:].split("=", 1)
            a[k] = v
    known = set(MUTATIONS) | {"time_unit", "tau_diag", "tau_nolapse", "rho_diag", "w_sign", "asym_sign", "pstar_off", "class_loose", "ellis_charge", "theta_coef", "scalar_coef",
                              "ghost_sign", "no_corrector", "data_constraint_sign", "finder_sign", "constraint_sign", "drop_bracket", "ed_t", "ra_scale", "xc_tau", "locus_diag"}
    if mut:
        if mut not in known:
            sys.exit(f"unknown mutation {mut}; known: {', '.join(sorted(known))}")
        if mut in ("ellis_charge", "theta_coef", "scalar_coef", "ghost_sign", "no_corrector", "data_constraint_sign", "finder_sign", "constraint_sign"):
            dn.MUT = mut
        elif mut == "class_loose":
            xp.EMUT = mut
        else:
            LMUT = mut
    import platform
    la = open("/proc/loadavg").read().split()[:3] if platform.system() == "Linux" else ["?"]
    print(f"throat_locus.py {cmd} {a}: numba {dn.nb.__version__ if dn.HAVE_NUMBA else 'absent'}, numpy {np.__version__}, python {platform.python_version()}, {platform.node()}, "
          f"load average at start {' '.join(la)}, mutation {mut}, {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    t0 = time.time()
    if cmd == "series":
        job_series(int(a["D"]), int(a["hk"]), a["out"], part=tuple(int(z) for z in a.get("part", "0/1").split("/")))
    elif cmd == "stall":
        job_stall(int(a["D"]), int(a["hk"]), a["out"])
    elif cmd == "expand":
        job_expand(int(a["D"]), int(a["hk"]), a["member"], float(a["X"]), a["out"], fam=a.get("family", "G"))
    elif cmd == "planted":
        job_planted(int(a["D"]), int(a["hk"]), float(a["v"]), a["out"])
    elif cmd == "analyze":
        bad = analyze_suite(load_record(a["raw"]))
        print(f"{time.time() - t0:.1f} s")
        sys.exit(1 if bad else 0)
    elif cmd in ("mutate", "dev"):
        bad = analyze_suite(mini_suite())
        if cmd == "mutate":
            want = set(MUTATIONS.get(mut, []))
            fired = {c.split("[")[0] for c in bad}
            print(f"mutation {mut}: declared target {sorted(want)}; fired {sorted(fired)}; " + ("targets all fired" if want <= fired else f"SILENT: {sorted(want - fired)}")
                  + ("" if fired <= want else f"; extra: {sorted(fired - want)}"))
        print(f"{time.time() - t0:.1f} s")
    else:
        sys.exit(f"unknown command {cmd}")


if __name__ == "__main__":
    main()

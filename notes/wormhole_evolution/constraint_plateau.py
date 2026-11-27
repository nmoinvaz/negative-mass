#!/usr/bin/env python3
"""Why the null-constraint maximum of #39's double-null code stops falling at the finest h (#226).

Max |C_++|, |C_--| (dn.constraints' formulas) over a fixed spacetime region, the slab {nodes of the rectangle [0, X]^2 with
t = (x+ + x-)/sqrt2 <= Tc}, streamed (the grid is never stored), on the static Ellis wormhole (no pulse), on #40's family S at
p*_h -/+ 1e-8 and at the two amplitudes that bracket p*_h. The candidate causes are tested by mutations of the run (noise, data,
spike, precision) and by a closed-form model of the constraint carried by roundoff. dn.py, checks.py and threshold.py are used unchanged
(read only). Criteria, regions, h and mutations: negative-mass-raw I226-rpi5/README.md, committed before the runs of record.

    python3 constraint_plateau.py controls [--mutate=NAME]                      K0-K3: the new code against dn.py's
    python3 constraint_plateau.py pstar --D=6 --h=H --X=X --guess=P --out=BASE  bisect p*_h of family S (extension h)
    python3 constraint_plateau.py run --job=NAME --table=FILE --out=DIR         one job of jobs() (json, npz of the slice)
    python3 constraint_plateau.py analyze --raw=DIR [--mutate=NAME]             every check and table from the stored runs
    python3 constraint_plateau.py jobs                                          the job names and their cost class

Engines: `row` is dn.row_step (numba, float64); `diag` is the numpy anti-diagonal sweep of dn.cell_np in float64 (bit-identical to the row
engine, K1) or in np.longdouble, which is IEEE binary128 on aarch64 (the higher-precision run); noise and spikes are added only there, to
a node as it is computed, as roundoff is.
Status words as in AGENTS.md: the cause is *numerical evidence*; the noise model is compared with the runs, not proved.
"""
import argparse
import glob
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dn           # noqa: E402
import threshold    # noqa: E402

SQ2 = math.sqrt(2.0)
SQ2H = 1.0 / SQ2
CRIT = threshold.CRIT
TW = 0.5            # width of a time bin of the profile
MUT = None          # --mutate=NAME of controls and analyze; names in MUTATIONS

# declared in negative-mass-raw I226-rpi5/README.md before the runs of record
CP = dict(
    H={5: [0.01 * SQ2H, 0.005 * SQ2H, 0.0025 * SQ2H, 0.00125 * SQ2H], 6: [0.005, 0.0025, 0.00125, 0.000625]},   # the fourth is the extension
    NDECL=3,                                                  # the first three are the issue's range
    Tc=dict(static=4.0, S={5: 10.0, 6: 8.0}, B={5: 16.0, 6: 12.0}),
    X=dict(static=6.0, S={5: 15.0, 6: 12.0}, B={5: 18.0, 6: 15.0}),
    off=1e-8,                                                 # the offset of family S from p*_h
    band=(2.0 ** 1.75, 2.0 ** 2.25),                          # ratio of a halving, as CRIT order_lo, order_hi
    slice_back=2,                                             # the slice of L1, M1 is at diagonal Sc - 2 (L1) or t = model_t (M1)
    model_t=6.0, model_band=(0.5, 2.0), model_min_nodes=20,
    c2_move=3.0, c2_coarse=1.10,                              # a = 256 raises the finest-h maximum by >= 3; a = 16 moves the coarsest by <= 1.10
    d1=(0.7, 1.4), l1_count=5, l2_dist=1.0,
    spike=(0.35, 0.12, 1e-10), corner=(2, 2, 1e-10), raypert=1e-7, noise_mut=256.0,
)
BANDS = [(0.0, 0.25), (0.25, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 3.0), (3.0, 5.0)]    # |l| bands of M1


# --------------------------------------------------------------------------------------------------------------------
# data: dn.ellis_rays with its RK4 constraint solve exposed
# --------------------------------------------------------------------------------------------------------------------
def pulses_of(family, A):
    """(pulses, symmetric, eps) of the member; 'static' is the unperturbed Ellis wormhole (no pulse)."""
    if family == "static":
        return [], False, (-1.0,)
    return threshold.family(family, A)


def ellis_rays_x(n, N, h, pulses=(), symmetric=False, eps=(-1.0,), dtype=np.float64, sub=1, dense=False, ray_pert=0.0):
    """dn.ellis_rays with: dtype (np.longdouble is IEEE binary128 on aarch64), sub (RK4 steps per grid step, the solve's own step h/sub),
    dense (r(l) from the dense static solution instead of the RK4; no pulses), ray_pert (a planted relative bump in r on both rays: the
    data-constraint mutation). dtype float64, sub 1, dense False, ray_pert 0 reproduces dn.ellis_rays bit for bit (K2)."""
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
    sol, q = dn.ellis_dense(n, N * h / SQ2 + 1.0)

    def phi_s(s, ray, field):
        out = np.zeros_like(s)
        if field == 0:
            r = sol(s / SQ2)[0]
            out = out + (1.0 if ray == "+" else -1.0) * q / r ** n / SQ2
        for p in pulses:
            if p["ray"] == ray and p["field"] == field:
                out = out + dn._pulse(s, p["A"], p["xc"], p["w"], p.get("kind", "phi"))[1]
        return out

    def phi_val(s, ray, field):
        out = np.zeros_like(s)
        if field == 0:
            out = out + (1.0 if ray == "+" else -1.0) * sol(s / SQ2)[2]
        for p in pulses:
            if p["ray"] == ray and p["field"] == field:
                out = out + dn._pulse(s, p["A"], p["xc"], p["w"], p.get("kind", "phi"))[0]
        return out

    nodes = np.arange(N + 1) * h
    rays = []
    for ray in ("+", "-"):
        def S_of(s):
            return sum(eps[f] * phi_s(s, ray, f) ** 2 for f in range(k))
        if dense:
            if pulses:
                raise ValueError("dense data needs the static case")
            r = np.asarray(sol(nodes / SQ2)[0], dtype)
        else:
            hs = h / sub
            fine = np.arange(N * sub + 1) * hs
            S0, Sh, S1 = S_of(fine[:-1]), S_of(fine[:-1] + 0.5 * hs), S_of(fine[1:])
            r = np.empty(N + 1, dtype)
            y0, v0 = dtype(1.0), dtype(0.0)
            r[0] = y0
            hh = dtype(hs)
            for i in range(N * sub):
                a0, b0 = y0, v0
                k1 = (b0, -(dtype(S0[i]) / n) * a0)
                a2, b2 = a0 + dtype(0.5) * hh * k1[0], b0 + dtype(0.5) * hh * k1[1]
                k2 = (b2, -(dtype(Sh[i]) / n) * a2)
                a3, b3 = a0 + dtype(0.5) * hh * k2[0], b0 + dtype(0.5) * hh * k2[1]
                k3 = (b3, -(dtype(Sh[i]) / n) * a3)
                a4, b4 = a0 + hh * k3[0], b0 + hh * k3[1]
                k4 = (b4, -(dtype(S1[i]) / n) * a4)
                y0 = a0 + hh / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
                v0 = b0 + hh / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
                if (i + 1) % sub == 0:
                    r[(i + 1) // sub] = y0
        if ray_pert:
            r = r * (1.0 + ray_pert * np.exp(-((nodes - 5.0) / 0.5) ** 2))
        U = np.zeros((2 + k, N + 1), dtype)
        U[0] = r
        for f in range(k):
            U[2 + f] = phi_val(nodes, ray, f)
        rays.append(U)
    return dn.Data(n, eps, N, h, rays[0], rays[1], name="ellis")


# --------------------------------------------------------------------------------------------------------------------
# constraints, as dn.constraints writes them, on one row or one diagonal
# --------------------------------------------------------------------------------------------------------------------
def sgn():
    return -1.0 if MUT == "cons_sign" else 1.0


def terms_pp(R, h, n, eps):
    """C_++ at nodes 1 .. N-1 of one row R (nf, N+1), and the sum of the absolute values of its terms."""
    r, f = R[0], R[1]
    rpp = (r[2:] - 2 * r[1:-1] + r[:-2]) / h ** 2
    fp = (f[2:] - f[:-2]) / (2 * h)
    rp = (r[2:] - r[:-2]) / (2 * h)
    sp = sum(eps[k] * ((R[2 + k][2:] - R[2 + k][:-2]) / (2 * h)) ** 2 for k in range(len(eps))) if len(eps) else 0.0
    t3 = sgn() * (r[1:-1] / n) * sp
    return rpp + fp * rp + t3, np.abs(rpp) + np.abs(fp * rp) + np.abs(t3)


def terms_mm(Rm, R0, Rp, h, n, eps):
    """C_-- at nodes 0 .. N of the row R0 from the rows j-1, j, j+1, and the sum of the absolute values of its terms."""
    r = R0[0]
    rmm = (Rp[0] - 2 * r + Rm[0]) / h ** 2
    fm = (Rp[1] - Rm[1]) / (2 * h)
    rm = (Rp[0] - Rm[0]) / (2 * h)
    sm = sum(eps[k] * ((Rp[2 + k] - Rm[2 + k]) / (2 * h)) ** 2 for k in range(len(eps))) if len(eps) else 0.0
    t3 = sgn() * (r / n) * sm
    return rmm + fm * rm + t3, np.abs(rmm) + np.abs(fm * rm) + np.abs(t3)


class Stats:
    """Row engine: per row j and per prefix edge If (nodes i <= If): max |C|, its node and the max sum of terms; per time bin: max and node."""

    def __init__(self, N, h, X, extra=(), xstep=0.5):
        self.N, self.h = N, h
        self.xf = np.arange(xstep, X + 1e-9, xstep)
        self.nx = len(self.xf)
        self.xstep = xstep
        self.If = np.concatenate([np.floor(self.xf / h + 1e-9).astype(int), np.array(extra, int)])
        K = len(self.If)
        z = lambda dt: np.zeros((N + 1, K), dt)
        self.pv, self.pi, self.ps = z(float), z(np.int32), z(float)
        self.mv, self.mi, self.ms = z(float), z(np.int32), z(float)
        self.nT = int(math.ceil(2 * N * h / SQ2 / TW)) + 2
        self.tv = {w: np.zeros(self.nT) for w in ("pp", "mm")}
        self.ti = {w: np.zeros(self.nT, np.int64) for w in ("pp", "mm")}
        self.tj = {w: np.zeros(self.nT, np.int64) for w in ("pp", "mm")}

    def add(self, which, j, C, S, dead_mask, off):
        """C, S over the nodes i = off + arange(L) of row j; dead_mask True where undefined."""
        a = np.where(dead_mask, 0.0, np.abs(C))
        s = np.where(dead_mask, 0.0, S)
        L = len(a)
        pm = np.maximum.accumulate(a)
        rec = np.maximum.accumulate(np.where(a == pm, np.arange(L), -1))
        sm = np.maximum.accumulate(s)
        idx = np.minimum(self.If - off, L - 1)
        ok = (self.If - off) >= 0
        v, ps_, pi_ = (self.pv, self.ps, self.pi) if which == "pp" else (self.mv, self.ms, self.mi)
        v[j, ok] = pm[idx[ok]]
        ps_[j, ok] = sm[idx[ok]]
        pi_[j, ok] = rec[idx[ok]] + off
        tb = np.floor(((off + np.arange(L)) + j) * self.h / SQ2 / TW).astype(int)
        starts = np.concatenate([[0], np.flatnonzero(np.diff(tb)) + 1])
        seg = np.maximum.reduceat(a, starts)
        bins = tb[starts]
        for q in np.flatnonzero(seg > self.tv[which][bins]):
            lo = starts[q]
            hi = starts[q + 1] if q + 1 < len(starts) else L
            k = lo + int(np.argmax(a[lo:hi]))
            b = bins[q]
            self.tv[which][b], self.ti[which][b], self.tj[which][b] = a[k], off + k, j

    def region(self, which, X_f=None, k=None, Jf=None):
        """(max |C|, node i, row j, max sum of terms) over the square [0, X_f]^2, or over column k of the prefix edges and rows j <= Jf."""
        if k is None:
            k = int(round(X_f / self.xstep)) - 1
        Jf = self.If[k] if Jf is None else Jf
        v, ii, ss = (self.pv, self.pi, self.ps) if which == "pp" else (self.mv, self.mi, self.ms)
        col = v[: Jf + 1, k]
        j = int(np.argmax(col))
        return float(col[j]), int(ii[j, k]), j, float(ss[: Jf + 1, k].max())


def stream(d, stop=False, jmax=None, extra=(), slices=(), tol=None, rexp=None, rcut=None, fcut=None, excise=True):
    """Evolve d with dn.row_step (the production kernel), three rows in memory. stop: end at the first event as dn.evolve(stop=True).
    jmax: end after row jmax. Returns (Stats, info); info['slices'][w][s] is C along the diagonal i + j = s."""
    N, h, nf, eps, n = d.N, d.h, d.nf, d.eps, d.n
    tol = tol or CRIT["trap_tol"]
    rexp = rexp or CRIT["rexp"]
    rcut = CRIT["rcut"] if rcut is None else rcut
    fcut = CRIT["fcut"] if fcut is None else fcut
    par = dn.kernel_params(n, eps, 1)
    st = Stats(N, h, N * h, extra)
    old = np.ascontiguousarray(d.Sp.copy())
    new = np.empty_like(old)
    dold = np.zeros(N + 1, np.uint8)
    dnew = np.zeros(N + 1, np.uint8)
    prev = prevd = None
    info = dict(events={}, rows=0, slices={w: {int(s0): np.full(N + 1, np.nan) for s0 in slices} for w in ("pp", "mm")})

    def take(w, j, C, off):
        for s0 in slices:
            i = int(s0) - j
            if 0 <= i - off < len(C) and 0 <= i <= N:
                info["slices"][w][int(s0)][i] = C[i - off]

    def dmask3(a, b, c):
        return (a | b | c).astype(bool)

    C, S = terms_pp(old, h, n, eps)
    st.add("pp", 0, C, S, dmask3(dold[2:], dold[1:-1], dold[:-2]), 1)
    take("pp", 0, C, 1)
    for j in range(N):
        new[:, 0] = d.Sm[:, j + 1]
        dnew[0] = 0
        dn.row_step(old, new, dold, dnew, h, par, eps, rcut, fcut, int(excise))
        C, S = terms_pp(new, h, n, eps)
        st.add("pp", j + 1, C, S, dmask3(dnew[2:], dnew[1:-1], dnew[:-2]), 1)
        take("pp", j + 1, C, 1)
        if prev is not None:                                   # C_-- of row j needs the rows j-1, j, j+1
            C, S = terms_mm(prev, old, new, h, n, eps)
            st.add("mm", j, C, S, dmask3(prevd, dold, dnew), 0)
            take("mm", j, C, 0)
        if stop and not info["events"]:
            c, up, um = dn.strip_cells(old, new, h)
            alive = ~(dold[:-1].astype(bool) | dold[1:].astype(bool) | dnew[:-1].astype(bool) | dnew[1:].astype(bool))
            cc = dn.cell_class(up[0], um[0], tol)
            trapped = alive & (cc == -1)
            if trapped.any():
                info["events"]["collapse"] = (j, int(np.argmax(trapped)))
            elif alive[j] and cc[j] == 1 and c[0, j] >= rexp:
                info["events"]["expand"] = (j, j)
            if info["events"]:
                info["rows"] = j + 1
                break
        prev, prevd = old.copy(), dold.copy()
        old, new = new, old
        dold, dnew = dnew, dold
        info["rows"] = j + 1
        if jmax is not None and j + 1 >= jmax:
            break
    return st, info


def stream_diag(d, dtype=np.float64, Sc=None, slices=(), noise=0.0, seed=0, spike=None):
    """Anti-diagonal sweep of the reference kernel dn.cell_np in `dtype`, three diagonals in memory; C_++, C_-- (in dtype) over the nodes
    with i + j <= Sc. float64 is the same arithmetic as the numba row sweep (K1); np.longdouble is binary128 here. noise: relative uniform
    noise of that amplitude added to every node as it is computed (the roundoff mutation); spike (i, j, rel): one node of r, once.
    Returns time-bin maxima (as Stats.tv, ti, tj) and the C fields on the diagonals `slices`."""
    N, h, nf, eps, n = d.N, d.h, d.nf, d.eps, d.n
    par = dn.kernel_params(n, eps, 1)
    Sc = 2 * N if Sc is None else min(int(Sc), 2 * N)
    rng = np.random.default_rng(seed) if noise else None
    Sp, Sm = d.Sp.astype(dtype), d.Sm.astype(dtype)
    D = {0: np.zeros((nf, N + 1), dtype)}
    D[0][:, 0] = Sp[:, 0]
    dead = {0: np.zeros(N + 1, bool)}
    nT = int(math.ceil((Sc + 1) * h / SQ2 / TW)) + 2
    tv = {w: np.zeros(nT) for w in ("pp", "mm")}
    ti = {w: np.zeros(nT, np.int64) for w in ("pp", "mm")}
    tj = {w: np.zeros(nT, np.int64) for w in ("pp", "mm")}
    sl = {w: {int(s0): np.full(N + 1, np.nan) for s0 in slices} for w in ("pp", "mm")}
    hd = dtype(h * (1.01 if (MUT == "quad_h" and dtype is np.longdouble) else 1.0))
    sg = sgn()

    def cons(s):
        out = {}
        Dm, D0, Dp = D[s - 1], D[s], D[s + 1]
        dm, d0, dp = dead[s - 1], dead[s], dead[s + 1]
        I = np.arange(max(1, s - N), min(N - 1, s) + 1)
        if len(I):
            r = D0[0, I]
            rpp = (Dp[0, I + 1] - 2 * r + Dm[0, I - 1]) / hd ** 2
            fp = (Dp[1, I + 1] - Dm[1, I - 1]) / (2 * hd)
            rp = (Dp[0, I + 1] - Dm[0, I - 1]) / (2 * hd)
            sp = sum(eps[k] * ((Dp[2 + k, I + 1] - Dm[2 + k, I - 1]) / (2 * hd)) ** 2 for k in range(len(eps)))
            C = rpp + fp * rp + sg * (r / n) * sp
            out["pp"] = (I, np.where(dm[I - 1] | d0[I] | dp[I + 1], np.nan, C))
        I = np.arange(max(0, s - N + 1), min(N, s - 1) + 1)
        if len(I):
            r = D0[0, I]
            rmm = (Dp[0, I] - 2 * r + Dm[0, I]) / hd ** 2
            fm = (Dp[1, I] - Dm[1, I]) / (2 * hd)
            rm = (Dp[0, I] - Dm[0, I]) / (2 * hd)
            sm = sum(eps[k] * ((Dp[2 + k, I] - Dm[2 + k, I]) / (2 * hd)) ** 2 for k in range(len(eps)))
            C = rmm + fm * rm + sg * (r / n) * sm
            out["mm"] = (I, np.where(dm[I] | d0[I] | dp[I], np.nan, C))
        return out

    for s in range(1, Sc + 2):
        nd = np.zeros((nf, N + 1), dtype)
        dd = np.zeros(N + 1, bool)
        I = np.arange(max(1, s - N), min(N, s - 1) + 1)
        if len(I):
            Sv = D[s - 1][:, I - 1] if MUT == "diag_S" else D[s - 2][:, I - 1]
            Ev, Wv = D[s - 1][:, I], D[s - 1][:, I - 1]
            ds = dead[s - 2][I - 1] | dead[s - 1][I] | dead[s - 1][I - 1]
            with np.errstate(all="ignore"):
                out = dn.cell_np(Sv, Ev, Wv, h, par, eps)
                if noise:
                    out = out * (1.0 + noise * rng.uniform(-1.0, 1.0, out.shape))
                bad = ~np.all(np.isfinite(out), axis=0)
                ds = ds | bad | (out[0] < CRIT["rcut"]) | (np.exp(-out[1]) > CRIT["fcut"])
            nd[:, I] = np.where(ds[None, :], Ev, out)
            dd[I] = ds
        if s <= N:
            nd[:, 0] = Sm[:, s]
            nd[:, s] = Sp[:, s]
        if spike is not None and spike[0] + spike[1] == s:
            nd[0, spike[0]] *= 1.0 + spike[2]
        D[s], dead[s] = nd, dd
        sc = s - 1                                                # the diagonal whose neighbours are now all known
        if sc >= 1:
            for w, (Iw, Cw) in cons(sc).items():
                a = np.abs(np.where(np.isnan(Cw), 0.0, Cw)).astype(np.float64)
                b = int(math.floor(sc * h / SQ2 / TW))
                k = int(np.argmax(a))
                if a[k] > tv[w][b]:
                    tv[w][b], ti[w][b], tj[w][b] = a[k], int(Iw[k]), sc - int(Iw[k])
                if sc in sl[w]:
                    sl[w][sc][Iw] = Cw.astype(np.float64)
        D.pop(s - 2, None)
        dead.pop(s - 2, None)
    return dict(tv=tv, ti=ti, tj=tj, slices=sl, nT=nT)


# --------------------------------------------------------------------------------------------------------------------
# one run over a slab
# --------------------------------------------------------------------------------------------------------------------
def slab_run(n, h, family, A, Tc, X, engine="row", dtype="f8", noise=0.0, seed=1, sub=1, dense=False, ray_pert=0.0, spike=None, slices=()):
    """Max |C_++|, |C_--| over the slab {nodes (i, j) <= X/h, t = (i + j) h / sqrt2 <= Tc}, per time bin and overall, with the node of each
    maximum. Returns a json-able record and the slices {'pp'|'mm': {s: array}}."""
    t0 = time.time()
    N = int(round(X / h))
    Sc = int(math.floor(Tc * SQ2 / h + 1e-9))
    pulses, sym, eps = pulses_of(family, A)
    ld = dtype == "f16"
    d = ellis_rays_x(n, N, h, pulses=pulses, symmetric=sym, eps=eps, dtype=np.longdouble if ld else np.float64, sub=sub, dense=dense,
                     ray_pert=ray_pert)
    if engine == "row":
        if ld or noise or spike:
            raise ValueError("the row engine is float64 without noise or spike")
        st, info = stream(d, jmax=Sc, slices=slices)
        tv, ti, tj, sl = st.tv, st.ti, st.tj, info["slices"]
    else:
        r = stream_diag(d, np.longdouble if ld else np.float64, Sc=Sc, slices=slices, noise=noise, seed=seed, spike=spike)
        tv, ti, tj, sl = r["tv"], r["ti"], r["tj"], r["slices"]
    nb = int(round(Tc / TW))
    bins = []
    for b in range(nb):
        rec = dict(t=(b + 0.5) * TW)
        for w in ("pp", "mm"):
            rec[w] = [float(tv[w][b]), int(ti[w][b]), int(tj[w][b])]
        bins.append(rec)
    v, w, b = max((bins[b][w][0], w, b) for b in range(nb) for w in ("pp", "mm"))
    i, j = bins[b][w][1], bins[b][w][2]
    rec = dict(D=n + 2, h=h, N=N, X=X, family=family, A=A, Tc=Tc, Sc=Sc, engine=engine, dtype=dtype, noise=noise, seed=seed, sub=sub, dense=dense,
               ray_pert=ray_pert, spike=spike, secs=round(time.time() - t0, 1), max=v, which=w, i=i, j=j, xp=i * h, xm=j * h,
               t=(i + j) * h / SQ2, bins=bins)
    return rec, sl


# --------------------------------------------------------------------------------------------------------------------
# the table of p*_h and bracket amplitudes, and the jobs
# --------------------------------------------------------------------------------------------------------------------
def pstar_bisect(n, h, X, guess, tol, out=None, rel=0.02):
    """Bisection of the fate of family S to hi - lo <= tol (threshold.bisect, threshold.point: #40's events and fates), from guess (1 -/+ rel)."""
    log = open(out, "a") if out else None

    def run(p):
        r = threshold.point(n, "S", h, X, p)
        if log:
            log.write(json.dumps(r) + "\n")
            log.flush()
        return r
    lo, hi, stopped = threshold.bisect(run, guess * (1 - rel), guess * (1 + rel), tol)
    return dict(n=n, h=h, X=X, lo=lo, hi=hi, stopped=stopped, pstar=0.5 * (lo + hi))


def jobs():
    """name -> spec. k indexes CP['H'][D]; the fourth h is the extension. Cost class: s (seconds), m (minutes), l (10 min or more)."""
    J = {}
    for D in (5, 6):
        for k in range(4):
            J[f"st_D{D}_k{k}"] = dict(D=D, k=k, set="static", engine="row")
        for k in range(3):
            J[f"stq_D{D}_k{k}"] = dict(D=D, k=k, set="static", engine="diag", dtype="f16")
            for a in (1, 16, 256):
                J[f"stn{a}_D{D}_k{k}"] = dict(D=D, k=k, set="static", engine="diag", noise=a)
        for side, tag in ((-1, "m"), (1, "p")):
            for k in range(4):
                J[f"s{tag}_D{D}_k{k}"] = dict(D=D, k=k, set="S", side=side, engine="row", slice=(k == 2 and side == -1))
            for k in range(3 if tag == "p" else 4):
                J[f"sq{tag}_D{D}_k{k}"] = dict(D=D, k=k, set="S", side=side, engine="diag", dtype="f16")
        for k in (0, 2):
            for a in (1, 16, 256):
                J[f"sn{a}m_D{D}_k{k}"] = dict(D=D, k=k, set="S", side=-1, engine="diag", noise=a)
        J[f"ssub4m_D{D}_k2"] = dict(D=D, k=2, set="S", side=-1, engine="row", sub=4)
        J[f"ssub16m_D{D}_k2"] = dict(D=D, k=2, set="S", side=-1, engine="row", sub=16)
        J[f"smutpert_D{D}_k2"] = dict(D=D, k=2, set="S", side=-1, engine="row", ray_pert=CP["raypert"])
        J[f"smutspike_D{D}_k2"] = dict(D=D, k=2, set="S", side=-1, engine="diag", spike=CP["spike"], slice=True)
        J[f"smutcorner_D{D}_k2"] = dict(D=D, k=2, set="S", side=-1, engine="diag", spike=CP["corner"])
        for side in ("lo", "hi"):
            for k in range(3):
                J[f"b{side}_D{D}_k{k}"] = dict(D=D, k=k, set="B", side=side, engine="row", slice=(k == 2 and side == "lo"))
        for k in range(3):
            J[f"bqlo_D{D}_k{k}"] = dict(D=D, k=k, set="B", side="lo", engine="diag", dtype="f16")
        for k in (1, 2):
            J[f"model_D{D}_k{k}"] = dict(D=D, k=k, set="S", side=-1, kind="model")
    return J


def est_secs(spec, table):
    """Cost estimate of a job on one rpi5 core (0.4 us per cell for the row engine, 4 us float64 and 6 us binary128 for the diagonal engine)."""
    n, h, fam, A, Tc, X = spec_args(spec, table)
    N = int(round(X / h))
    Sc = int(math.floor(Tc * SQ2 / h))
    if spec.get("kind") == "model":
        s0 = int(round(CP["model_t"] * SQ2 / h))
        return 6e-6 * s0 * s0 / 2 + 8 * 4e-6 * s0 * s0 / 2 * 0.5 + 60.0
    if spec.get("engine", "row") == "row":
        return 0.4e-6 * N * min(N, Sc) * (spec.get("sub", 1) ** 0.2) + 5.0
    cells = min(Sc, N) ** 2 / 2 + max(0, Sc - N) * N
    return (6e-6 if spec.get("dtype") == "f16" else 1.5e-6) * cells + 5.0


def spec_args(spec, table):
    """(n, h, family, A, Tc, X) of a job."""
    D, k, st = spec["D"], spec["k"], spec["set"]
    T = table[str(D)]
    h = T["h"][k]
    if st == "static":
        return D - 2, h, "static", 0.0, CP["Tc"]["static"], CP["X"]["static"]
    if st == "S":
        return D - 2, h, "S", T["pstar"][k] + spec["side"] * CP["off"], CP["Tc"]["S"][D], CP["X"]["S"][D]
    return D - 2, h, "S", T[spec["side"]][k], CP["Tc"]["B"][D], CP["X"]["B"][D]


def run_job(name, table, out):
    spec = jobs()[name]
    if spec.get("kind") == "model":
        rec = model_job(spec, table)
        rec["spec"] = spec
        json.dump(rec, open(os.path.join(out, name + ".json"), "w"), indent=1)
        return rec
    n, h, fam, A, Tc, X = spec_args(spec, table)
    Sc = int(math.floor(Tc * SQ2 / h + 1e-9))
    sl = (Sc - CP["slice_back"],) if spec.get("slice") else ()
    sp = spec.get("spike")
    if sp and isinstance(sp[0], float):                         # fractions of the slab's last diagonal -> nodes
        sp = (int(sp[0] * Sc), int(sp[1] * Sc), sp[2])
    rec, slices = slab_run(n, h, fam, A, Tc, X, engine=spec.get("engine", "row"), dtype=spec.get("dtype", "f8"), noise=spec.get("noise", 0.0) * 2.0 ** -52,
                           sub=spec.get("sub", 1), ray_pert=spec.get("ray_pert", 0.0), spike=sp, slices=sl)
    rec["spec"] = spec
    rec["spike"] = sp
    json.dump(rec, open(os.path.join(out, name + ".json"), "w"), indent=1)
    if sl:
        np.savez_compressed(os.path.join(out, name + ".npz"), s0=sl[0], pp=slices["pp"][sl[0]], mm=slices["mm"][sl[0]])
    return rec


# --------------------------------------------------------------------------------------------------------------------
# the model: the constraint carried by roundoff
# --------------------------------------------------------------------------------------------------------------------
def model_pred(n, h, N, s0, expo):
    """rms of C_++ on the diagonal i + j = s0 from uniform relative noise of amplitude 2^-52 on r in every node, along the characteristics:
    C r^(n-1) is conserved along x- (the evolution equations), so the noise injected in row j' reaches (i, j) with the factor (r(i,j')/r(i,j))^(n-1);
    a node's noise d gives the pair +d/h^2, -d/h^2 in C_++ at i-1, i on every later row. rms = sqrt(2/3) eps/h^2 r(i,j)^-expo (sum_j' r(i,j')^(2n))^(1/2)
    with expo = n - 1 (r(i, j') is the static r(l), l = (i - j') h / sqrt2)."""
    sol, q = dn.ellis_dense(n, N * h / SQ2 + 1.0)
    rs = lambda l: sol(np.abs(l))[0]
    I = np.arange(max(1, s0 - N), min(N - 1, s0 - 1) + 1)
    J = s0 - I
    pred = np.zeros(len(I))
    for k, (i, j) in enumerate(zip(I, J)):
        rr = rs((i - np.arange(0, j)) * h / SQ2)
        pred[k] = math.sqrt(2.0 / 3.0) * 2.0 ** -52 / h ** 2 * rs((i - j) * h / SQ2) ** (-expo) * math.sqrt(np.sum(rr ** (2 * n)))
    return I, (I - J) * h / SQ2, pred


def model_job(spec, table):
    """Natural roundoff (float64 minus binary128), and injected noise a = 16, 256 (noisy minus float64), on the diagonal t = model_t of
    family S at p*_h - 1e-8: rms of the C_++ and C_-- difference in |l| bands against model_pred."""
    n, h, fam, A, Tc, X = spec_args(spec, table)
    N = int(round(X / h))
    s0 = int(round(CP["model_t"] * SQ2 / h))
    pulses, sym, eps = pulses_of(fam, A)
    d8 = ellis_rays_x(n, N, h, pulses=pulses, symmetric=sym, eps=eps)
    d16 = ellis_rays_x(n, N, h, pulses=pulses, symmetric=sym, eps=eps, dtype=np.longdouble)
    base = stream_diag(d8, np.float64, Sc=s0 + 2, slices=(s0,))["slices"]
    quad = stream_diag(d16, np.longdouble, Sc=s0 + 2, slices=(s0,))["slices"]
    fields = {"natural": {w: base[w][s0] - quad[w][s0] for w in ("pp", "mm")}}
    amp = {"natural": 1.0}
    for a in (16, 256):
        nz = stream_diag(d8, np.float64, Sc=s0 + 2, slices=(s0,), noise=a * 2.0 ** -52, seed=1)["slices"]
        fields[f"a{a}"] = {w: nz[w][s0] - base[w][s0] for w in ("pp", "mm")}
        amp[f"a{a}"] = float(a)
    I, l, pred = model_pred(n, h, N, s0, n - 1)
    _, _, pred_wrong = model_pred(n, h, N, s0, n)
    out = dict(D=n + 2, h=h, s0=s0, t=s0 * h / SQ2, bands=[])
    for lo, hi in BANDS:
        m = (np.abs(l) >= lo) & (np.abs(l) < hi)
        row = dict(lo=lo, hi=hi, nodes=int(m.sum()))
        for nm, F in fields.items():
            for w in ("pp", "mm"):
                f = F[w][I][m]
                f = f[np.isfinite(f)]
                row[f"{nm}_{w}"] = float(np.sqrt(np.mean(f ** 2))) if len(f) else None
        row["model"] = float(np.sqrt(np.mean(pred[m] ** 2))) if m.any() else None
        row["model_wrong"] = float(np.sqrt(np.mean(pred_wrong[m] ** 2))) if m.any() else None
        out["bands"].append(row)
    out["amp"] = amp
    return out


# --------------------------------------------------------------------------------------------------------------------
# controls
# --------------------------------------------------------------------------------------------------------------------
RESULTS = []


def check(cid, text, ok, detail=""):
    ok = None if ok is None else bool(ok)
    RESULTS.append((cid, ok))
    print(f"[{'SKIP' if ok is None else ('pass' if ok else 'FAIL')}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def controls():
    """K0 streaming constraints = dn.constraints; K1 float64 diagonal sweep = numba row sweep; K2 ellis_rays_x = dn.ellis_rays;
    K3 binary128 and float64 sweeps agree where truncation dominates."""
    worst = 0.0
    for fam, A, n, h, X in (("static", 0.0, 3, 0.05, 6.0), ("S", 0.02, 3, 0.05, 8.0), ("G", 0.02, 4, 0.05, 6.0)):
        pulses, sym, eps = pulses_of(fam, A) if fam != "G" else threshold.family("G", A)
        N = int(round(X / h))
        d = dn.ellis_rays(n, N, h, pulses=pulses, symmetric=sym, eps=eps)
        res = dn.evolve(d, store=True, stop=False)
        Cpp, Cmm = dn.constraints(res)
        st, info = stream(d, stop=False)
        for k in range(st.nx):
            If = st.If[k]
            a = np.abs(np.nan_to_num(Cpp, nan=0.0))[: max(If, 0), : If + 1]
            worst = max(worst, abs((a.max() if a.size else 0.0) - st.pv[: If + 1, k].max()))
            a = np.abs(np.nan_to_num(Cmm, nan=0.0))[: If + 1, : If]
            worst = max(worst, abs((a.max() if a.size else 0.0) - st.mv[: If + 1, k].max()))
    check("K0", "streamed prefix maxima of |C_++|, |C_--| equal dn.constraints' (static, family S, family G; numba kernel)", worst == 0.0, f"largest difference {worst:.1e}")
    ok1, bad1 = True, 0.0
    for fam, A, n, h, X in (("static", 0.0, 3, 0.05, 6.0), ("S", 0.02, 3, 0.05, 8.0)):
        pulses, sym, eps = pulses_of(fam, A)
        N = int(round(X / h))
        d = dn.ellis_rays(n, N, h, pulses=pulses, symmetric=sym, eps=eps)
        st, info = stream(d, stop=False, slices=(40, 100))
        dg = stream_diag(d, np.float64, slices=(40, 100))
        nb = int(N * h / SQ2 / TW) - 1
        for w in ("pp", "mm"):
            ok1 = ok1 and np.array_equal(st.tv[w][:nb], dg["tv"][w][:nb]) and np.array_equal(st.ti[w][:nb], dg["ti"][w][:nb]) and np.array_equal(st.tj[w][:nb], dg["tj"][w][:nb])
            for s0 in (40, 100):
                ok1 = ok1 and np.array_equal(np.nan_to_num(info["slices"][w][s0]), np.nan_to_num(dg["slices"][w][s0]))
                bad1 = max(bad1, float(np.nanmax(np.abs(info["slices"][w][s0] - dg["slices"][w][s0]))))
    check("K1", "the float64 anti-diagonal sweep reproduces the numba row sweep: time-bin maxima, their nodes and two diagonals, bit for bit", ok1, f"largest slice difference {bad1:.1e}")
    ok2, bad2 = True, 0.0
    for fam, A, n, h, X in (("static", 0.0, 3, 0.02, 6.0), ("S", 0.02, 3, 0.02, 6.0), ("G", 0.02, 4, 0.02, 6.0)):
        pulses, sym, eps = pulses_of(fam, A) if fam != "G" else threshold.family("G", A)
        N = int(round(X / h))
        d0 = dn.ellis_rays(n, N, h, pulses=pulses, symmetric=sym, eps=eps)
        d1 = ellis_rays_x(n, N, h, pulses=pulses, symmetric=sym, eps=eps, sub=2 if MUT == "rk4_step" else 1)
        bad2 = max(bad2, float(np.abs(d0.Sp - d1.Sp).max()), float(np.abs(d0.Sm - d1.Sm).max()))
    check("K2", "ellis_rays_x (float64, own RK4 step = the grid's) reproduces dn.ellis_rays bit for bit (static, S, G)", bad2 == 0.0, f"largest difference {bad2:.1e}")
    n, h, X, Tc = 4, 0.005, 6.0, 4.0
    N, Sc = int(round(X / h)), int(Tc * SQ2 / h)
    d8 = ellis_rays_x(n, N, h)
    d16 = ellis_rays_x(n, N, h, dtype=np.longdouble)
    a8 = stream_diag(d8, np.float64, Sc=Sc)
    a16 = stream_diag(d16, np.longdouble, Sc=Sc)
    nb = int(Tc / TW)
    m8 = max(a8["tv"][w][:nb].max() for w in ("pp", "mm"))
    m16 = max(a16["tv"][w][:nb].max() for w in ("pp", "mm"))
    check("K3", "binary128 and float64 sweeps give the same slab maximum where truncation dominates (D = 6, h = 0.005, static, t <= 4)", abs(m16 / m8 - 1) < 2e-3,
          f"float64 {m8:.6e}, binary128 {m16:.6e}, ratio - 1 = {m16 / m8 - 1:.1e}, tolerance 2e-3")
    return RESULTS


# --------------------------------------------------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------------------------------------------------
MUTATIONS = {   # name: the checks that must fire; the effect of each is computed in the development runs (negative-mass-raw I226-rpi5/dev)
    "cons_sign": ["K0"], "diag_S": ["K1"], "rk4_step": ["K2"], "quad_h": ["K3"],
    "static_noise256": ["R1"], "quad_as_double": ["C1"], "noise_inert": ["C2"], "ray_pert": ["D1"], "spike": ["L1"], "corner": ["L2"], "model_exp": ["M1"],
}
CONTROL_MUTATIONS = ("cons_sign", "diag_S", "rk4_step", "quad_h")


def load(raw):
    R = {}
    for fn in glob.glob(os.path.join(raw, "*.json")):
        nm = os.path.basename(fn)[:-5]
        try:
            R[nm] = json.load(open(fn))
        except Exception:
            pass
    return R


def ratios(maxes):
    return [maxes[k] / maxes[k + 1] for k in range(len(maxes) - 1)]


def in_band(x):
    return CP["band"][0] <= x <= CP["band"][1]


def where(r):
    return f"{r['which']} at (i {r['i']}, j {r['j']}), x+ {r['xp']:.3f}, x- {r['xm']:.3f}, t {r['t']:.3f}"


def table(R, title, fmt, D, label, ndecl=3):
    """Print the slab maxima of the jobs fmt.format(D=D, k=k), their locations and the ratios; return the declared-h ratios."""
    recs = [R.get(fmt.format(D=D, k=k)) for k in range(4)]
    print(f"  {title}, D = {D}, region {label}")
    mx = []
    for k, r in enumerate(recs):
        if r is None:
            if k < ndecl:
                print(f"    h = {CP['H'][D][k]:.8f}: not run")
            continue
        mx.append(r["max"])
        tag = "  (extension)" if k >= ndecl else ""
        print(f"    h = {r['h']:.8f}: max |C| {r['max']:.4e}  {where(r)}{tag}")
    rs = ratios(mx)
    print("    ratios per halving: " + ", ".join(f"{x:.3f}{'' if in_band(x) else ' (outside the band)'}" for x in rs) + f"; band [{CP['band'][0]:.3f}, {CP['band'][1]:.3f}]")
    return rs


def analyze(raw):
    R = load(raw)
    mut = MUT
    fired = set()
    print(f"constraint_plateau.py analyze: {len([k for k in R if k in jobs()])} of {len(jobs())} jobs found in {raw}" + (f"; mutation {mut}" if mut else ""))
    print("slab = nodes of [0, X]^2 with t = (x+ + x-)/sqrt2 <= Tc; ratios are max|C|(h)/max|C|(h/2) over the same slab; second order is 4, the band is CRIT's order 1.75 to 2.25.\n")

    def verdict(cid, text, ok, detail=""):
        if not ok:
            fired.add(cid)
        print(f"[{'pass' if ok else 'FAIL'}] {cid}: {text}" + (f"  ({detail})" if detail else ""))

    # R1: the static case
    print("== R1 the unperturbed static wormhole (no pulse)")
    r1_ok, r1_det = True, []
    for D in (5, 6):
        fmt = "stn256_D{D}_k{k}" if mut == "static_noise256" else "st_D{D}_k{k}"
        rs = table(R, "static" + (" with 256 x 2^-52 noise (mutation)" if mut == "static_noise256" else ""), fmt, D, f"t <= {CP['Tc']['static']}, X = {CP['X']['static']}")
        okD = all(in_band(x) for x in rs[: CP["NDECL"] - 1])
        r1_ok = r1_ok and okD
        r1_det.append(f"D = {D}: " + ", ".join(f"{x:.3f}" for x in rs[: CP["NDECL"] - 1]))
    verdict("R1", "static case: max|C| over the fixed slab falls by a second-order ratio at each halving of the three declared h, D = 5 and 6", r1_ok, "; ".join(r1_det))

    # R2: family S at p* -/+ 1e-8 (a measurement; its band verdict is reported, not gated)
    print("\n== R2 family S at p*_h -/+ 1e-8")
    for D in (5, 6):
        for tag, nm in (("m", "p* - 1e-8"), ("p", "p* + 1e-8")):
            rs = table(R, f"S {nm}", f"s{tag}_D{{D}}_k{{k}}", D, f"t <= {CP['Tc']['S'][D]}, X = {CP['X']['S'][D]}")
            print(f"    band verdict at the declared h: {'met' if all(in_band(x) for x in rs[:2]) else 'not met'}")
    print("\n== R3 the two amplitudes that bracket p*_h (#40's lo and hi; their fates depend on roundoff and differ between machines)")
    for D in (5, 6):
        for side in ("lo", "hi"):
            rs = table(R, f"bracket {side}", f"b{side}_D{{D}}_k{{k}}", D, f"t <= {CP['Tc']['B'][D]}, X = {CP['X']['B'][D]}")
            print(f"    band verdict at the declared h: {'met' if all(in_band(x) for x in rs[:2]) else 'not met'}")

    # C1: binary128
    print("\n== C1 the same slabs in binary128 (the evolution and the constraint in 113-bit precision, data included)")
    c1_ok, c1_det = True, []
    for D in (5, 6):
        for tag, label, dfmt, qfmt in (("static", "static", "st_D{D}_k{k}", "stq_D{D}_k{k}"), ("S -", "S, p* - 1e-8", "sm_D{D}_k{k}", "sqm_D{D}_k{k}"),
                                         ("S +", "S, p* + 1e-8", "sp_D{D}_k{k}", "sqp_D{D}_k{k}"), ("B lo", "bracket lo", "blo_D{D}_k{k}", "bqlo_D{D}_k{k}")):
            if mut == "quad_as_double":
                qfmt = dfmt
            q = [R.get(qfmt.format(D=D, k=k)) for k in range(4)]
            dd = [R.get(dfmt.format(D=D, k=k)) for k in range(4)]
            idx = [k for k in range(4) if q[k] is not None and dd[k] is not None]
            if not idx:
                continue
            qs = [q[k]["max"] for k in idx]
            ds = [dd[k]["max"] for k in idx]
            print(f"  {label}, D = {D}: " + "; ".join(f"h = {q[k]['h']:.6f}: float64 {dd[k]['max']:.3e}, binary128 {q[k]['max']:.3e} (x{dd[k]['max'] / q[k]['max']:.3f}) {where(q[k])}" for k in idx))
            qr = ratios(qs[:3])
            ok = all(in_band(x) for x in qr)
            c1_ok = c1_ok and ok
            c1_det.append(f"{label} D={D} " + ", ".join(f"{x:.3f}" for x in qr))
            if len(idx) > 3 or (len(idx) == 4):
                print(f"    extension halving: binary128 ratio {qs[2] / qs[3]:.3f}, float64 ratio {ds[2] / ds[3]:.3f}")
    verdict("C1", "binary128: second-order ratios at every halving of the declared h, in every set (static, S -, S +, bracket lo), D = 5 and 6", c1_ok, "; ".join(c1_det))

    # C2: the roundoff mutation
    print("\n== C2 added noise, uniform relative amplitude a x 2^-52 on every node as it is computed (a = 1 doubles the roundoff), family S at p* - 1e-8;")
    print("   the numbers for a = 1 and 16 are reported; the realization-to-realization scatter of the maximum is about 13 % (D1)")
    c2_ok, c2_det = True, []
    for D in (5, 6):
        base = {k: R.get(f"sm_D{D}_k{k}") for k in (0, 2)}
        nz = {(a, k): R.get(f"sn{a}m_D{D}_k{k}") for a in (1, 16, 256) for k in (0, 2)}
        if mut == "noise_inert":
            nz = {(a, k): base[k] for a in (1, 16, 256) for k in (0, 2)}
        if any(v is None for v in list(base.values()) + list(nz.values())):
            continue
        line = []
        for k in (0, 2):
            line.append(f"h = {base[k]['h']:.6f}: a = 0 {base[k]['max']:.3e}, " + ", ".join(f"a = {a} {nz[(a, k)]['max']:.3e} (x{nz[(a, k)]['max'] / base[k]['max']:.3f})" for a in (1, 16, 256)))
        print(f"  D = {D}: " + "; ".join(line))
        mov = nz[(256, 2)]["max"] / base[2]["max"]
        crs = nz[(16, 0)]["max"] / base[0]["max"]
        ok = mov >= CP["c2_move"] and crs <= CP["c2_coarse"]
        c2_ok = c2_ok and ok
        c2_det.append(f"D={D}: a=256 at the finest declared h x{mov:.3f} (>= {CP['c2_move']}), a=16 at the coarsest x{crs:.3f} (<= {CP['c2_coarse']})")
    verdict("C2", "added roundoff moves the finest declared h's maximum (a = 256) and leaves the coarsest (a = 16)", c2_ok, "; ".join(c2_det))
    print("  static slab, same noise:")
    for D in (5, 6):
        for a in (1, 16, 256):
            mx = [R.get(f"stn{a}_D{D}_k{k}") for k in range(3)]
            if all(mx):
                print(f"    D = {D}, a = {a}: max |C| " + ", ".join(f"{x['max']:.3e}" for x in mx) + "; ratios " + ", ".join(f"{x:.3f}" for x in ratios([x['max'] for x in mx])))

    # D1: the ray data
    print("\n== D1 the ray data's RK4 constraint solve at its own step (h/sub), family S at p* - 1e-8, finest declared h")
    d1_ok, d1_det = True, []
    for D in (5, 6):
        base = R.get(f"sm_D{D}_k2")
        for nm in (("smutpert" if mut == "ray_pert" else None), "ssub4m", "ssub16m"):
            if nm is None:
                continue
            r = R.get(f"{nm}_D{D}_k2")
            if base is None or r is None:
                continue
            x = r["max"] / base["max"]
            ok = CP["d1"][0] <= x <= CP["d1"][1]
            d1_ok = d1_ok and ok
            d1_det.append(f"D={D} {nm}: x{x:.3f}")
            print(f"  D = {D}, {nm}: max |C| {r['max']:.4e} against {base['max']:.4e} (x{x:.3f}); {where(r)}")
    verdict("D1", f"the maximum is unchanged by the data's solve step (h/4, h/16) within {CP['d1']}", d1_ok, "; ".join(d1_det))

    # L1, L2: a single cell, the corner
    print("\n== L1, L2 where the maximum is")
    l1_ok, l1_det = True, []
    srcs = {5: ("smutspike_D5_k2" if mut == "spike" else "sm_D5_k2", "blo_D5_k2"), 6: ("smutspike_D6_k2" if mut == "spike" else "sm_D6_k2", "blo_D6_k2")}
    for D in (5, 6):
        for nm in srcs[D]:
            r = R.get(nm)
            fn = os.path.join(raw, nm + ".npz")
            if r is None or not os.path.exists(fn):
                continue
            z = np.load(fn)
            for w in ("pp", "mm"):
                a = np.abs(np.nan_to_num(z[w]))
                cnt = int((a >= 0.5 * a.max()).sum())
                ok = cnt >= CP["l1_count"]
                l1_ok = l1_ok and ok
                l1_det.append(f"{nm} {w}: {cnt}")
                print(f"  {nm}, diagonal i + j = {int(z['s0'])} (t = {int(z['s0']) * r['h'] / SQ2:.3f}), C_{w}: max {a.max():.3e}, nodes with |C| >= half the maximum: {cnt}")
    verdict("L1", f"the maximum on the last diagonal of the slab is not a single cell: at least {CP['l1_count']} nodes within a factor 2 of it", l1_ok, "; ".join(l1_det))
    l2_ok, l2_det = True, []
    for D in (5, 6):
        names = [f"st_D{D}_k{k}" for k in range(3)] + [f"sm_D{D}_k{k}" for k in range(3)] + [f"sp_D{D}_k{k}" for k in range(3)] + [f"blo_D{D}_k{k}" for k in range(3)] + [f"bhi_D{D}_k{k}" for k in range(3)]
        if mut == "corner":
            names = [f"smutcorner_D{D}_k2"]
        for nm in names:
            r = R.get(nm)
            if r is None:
                continue
            dist = min(r["xp"], r["xm"])
            ok = dist >= CP["l2_dist"]
            l2_ok = l2_ok and ok
            if not ok or mut == "corner":
                l2_det.append(f"{nm}: min(x+, x-) = {dist:.3f}")
                print(f"  {nm}: the maximum is at x+ {r['xp']:.3f}, x- {r['xm']:.3f}")
    verdict("L2", f"the slab maximum is at least {CP['l2_dist']} from the corner (x+ and x- both >= {CP['l2_dist']}), every declared-h run of every set", l2_ok, "; ".join(l2_det))

    # M1: the model
    print("\n== M1 the noise model against the runs: rms of the extra C_++, C_-- on the diagonal t = 6, by |l| band")
    m1_ok, m1_det = True, []
    for D in (5, 6):
        for k in (1, 2):
            r = R.get(f"model_D{D}_k{k}")
            if r is None:
                continue
            print(f"  D = {D}, h = {r['h']:.6f}, diagonal t = {r['t']:.3f}: |l| band, nodes, model (a = 1), then observed/model for natural roundoff, a = 16, a = 256, each C_++ / C_--")
            for b in r["bands"]:
                if b["nodes"] < CP["model_min_nodes"]:
                    continue
                mdl = b["model_wrong"] if mut == "model_exp" else b["model"]
                cells = []
                for nm, a in (("natural", 1.0), ("a16", 16.0), ("a256", 256.0)):
                    for w in ("pp", "mm"):
                        ratio = b[f"{nm}_{w}"] / (a * mdl)
                        cells.append(f"{ratio:.3f}")
                        ok = CP["model_band"][0] <= ratio <= CP["model_band"][1]
                        m1_ok = m1_ok and ok
                        if not ok:
                            m1_det.append(f"D={D} k={k} band [{b['lo']},{b['hi']}) {nm} {w}: {ratio:.3f}")
                print(f"    [{b['lo']:4.2f},{b['hi']:4.2f}) {b['nodes']:5d}  {mdl:.3e}   " + "  ".join(cells))
    verdict("M1", f"the extra constraint (natural roundoff, a = 16, a = 256) has the rms the model predicts, within {CP['model_band']} in every |l| band", m1_ok, "; ".join(m1_det[:6]))

    print("\n== summary")
    gated = {"R1", "C1", "C2", "D1", "L1", "L2", "M1"}
    if mut:
        want = set(MUTATIONS[mut])
        print(f"mutation {mut}: expected to fire {sorted(want)}; fired {sorted(fired & gated)}")
        return 0 if (fired & gated) == want else 1
    print("fired (validation checks): " + (", ".join(sorted(fired & gated)) or "none"))
    return 1 if fired & gated else 0


def main():
    global MUT
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["controls", "pstar", "run", "analyze", "jobs", "table"])
    ap.add_argument("--kmax", type=int, default=None)
    ap.add_argument("--kmin", type=int, default=None)
    ap.add_argument("--emit", default=None, help="jobs: print queue lines whose output goes to this directory")
    ap.add_argument("--raw40", default=os.path.expanduser("~/Source/negative-mass-raw/I40-dh610/runs"))
    ap.add_argument("--ext", default=None, help="table: the pstar json of the extension h, D = 5 then D = 6, comma separated")
    ap.add_argument("--mutate", default=None)
    ap.add_argument("--raw", default=None)
    ap.add_argument("--job", default=None)
    ap.add_argument("--table", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--D", type=int)
    ap.add_argument("--h", type=float)
    ap.add_argument("--X", type=float)
    ap.add_argument("--guess", type=float)
    ap.add_argument("--tol", type=float, default=2e-10)
    o = ap.parse_args()
    if o.mutate and o.mutate not in MUTATIONS:
        raise SystemExit(f"unknown mutation {o.mutate}; known: {sorted(MUTATIONS)}")
    MUT = o.mutate
    if o.cmd == "jobs":
        table = json.load(open(o.table)) if o.table else None
        rows = []
        for nm, sp in jobs().items():
            if (o.kmax is not None and sp.get("k", 0) > o.kmax) or (o.kmin is not None and sp.get("k", 0) < o.kmin):
                continue
            est = est_secs(sp, table) if table else 0.0
            rows.append((est, nm, sp))
        for est, nm, sp in sorted(rows, key=lambda r: -r[0]):          # longest first
            if o.emit:
                print(f"{nm}|timeout 7200 {sys.executable} -u constraint_plateau.py run --job={nm} --table={o.table} --out={o.emit}")
            else:
                print(f"{nm}  est {est:7.0f} s  {json.dumps(sp)}")
        return 0
    if o.cmd == "controls":
        if MUT and MUT not in CONTROL_MUTATIONS:
            raise SystemExit("--mutate for controls: " + ", ".join(CONTROL_MUTATIONS))
        res = controls()
        bad = {c for c, ok in res if ok is False}
        if MUT:
            want = set(MUTATIONS[MUT])
            print(f"mutation {MUT}: expected to fire {sorted(want)}; fired {sorted(bad)}")
            return 0 if bad == want else 1
        return 1 if bad else 0
    if o.cmd == "table":
        ext = dict(zip((5, 6), o.ext.split(","))) if o.ext else {}
        T = {}
        for D in (5, 6):
            T[str(D)] = dict(h=CP["H"][D], pstar=[], lo=[], hi=[], pstar_src=[])
            for k in range(3):
                nm = f"D{D}_S_h{CP['H'][D][k]:.8g}"
                cand = [d for d in glob.glob(os.path.join(o.raw40, f"D{D}_S_h*")) if abs(float(os.path.basename(d).split("_h")[1]) / CP["H"][D][k] - 1) < 1e-5]
                s = json.load(open(os.path.join(cand[0], "summary.json")))
                assert abs(s["h"] / CP["H"][D][k] - 1) < 1e-12 and s["family"] == "S" and s["D"] == D
                T[str(D)]["lo"].append(s["lo"])
                T[str(D)]["hi"].append(s["hi"])
                T[str(D)]["pstar"].append(0.5 * (s["lo"] + s["hi"]))
                T[str(D)]["pstar_src"].append(f"I40-dh610 {os.path.basename(cand[0])}")
            if D in ext:
                e = json.load(open(ext[D]))
                assert abs(e["h"] / CP["H"][D][3] - 1) < 1e-12
                T[str(D)]["pstar"].append(e["pstar"])
                T[str(D)]["pstar_src"].append(f"this issue, bisection to {e['hi'] - e['lo']:.1e}")
        print(json.dumps(T, indent=1))
        return 0
    if o.cmd == "pstar":
        rec = pstar_bisect(o.D - 2, o.h, o.X, o.guess, o.tol, out=o.out + ".points.jsonl" if o.out else None)
        if o.out:
            json.dump(rec, open(o.out + ".json", "w"), indent=1)
        print(json.dumps(rec))
        return 0
    if o.cmd == "run":
        table = json.load(open(o.table))
        os.makedirs(o.out, exist_ok=True)
        rec = run_job(o.job, table, o.out)
        print(json.dumps({k: v for k, v in rec.items() if k not in ("bins", "bands")}))
        return 0
    return analyze(o.raw)


if __name__ == "__main__":
    sys.exit(main())

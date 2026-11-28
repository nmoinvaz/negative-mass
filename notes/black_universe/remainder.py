#!/usr/bin/env python3
"""The remainder Y_rem of smooth c = -1 data at fixed rho against the advanced time (#322).

Smooth data of the c = -1 black universe have Y = c1 R + c2 V R + c2' U R + Y_rem (#302), with R = R(rho(Pi)) the areal radius, Pi = -V U, and Y_rem vanishing on both horizons.
This script evolves the classes of #248 and #302 in the Kruskal chart of #248 to V = VMAX (a streaming march: only the horizon line, the columns around S and a window of nodes
around each fixed-rho curve are kept), forms Y_rem / R = Y / R - (c1 + c2 V + c2' U) with (c1, c2, c2') read at S, and fits its rate against ln V = kappa v at three fixed rho.

Reads horizon_sector.py, horizon_basis.py, full_slice.py, static_c1.py, linear.py (read-only imports). Python: ~/Source/negative-mass/.venv/bin/python.
"""
import math
import os
import sys
import time

import numba as nb
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
_ARGV = sys.argv
sys.argv = sys.argv[:1]                     # the imported modules parse --mutate / --only themselves
import horizon_sector as hs  # noqa: E402
import horizon_basis as hb  # noqa: E402
sys.argv = _ARGV

SQ2 = math.sqrt(2.0)
C = -1.0
KAPPA = 2.0 / math.pi
NW = 4                                      # half width of the interpolation window: 2 NW = 8 nodes (7th degree Lagrange in U)


# ----------------------------------------------------------------------------------------------------------------------
# the streaming march: hs.march_K's cells, two columns in memory
# ----------------------------------------------------------------------------------------------------------------------
@nb.njit(cache=True)
def march_stream(UL, UR, V0, U0, par, tab, h, Pcap, i_h, Pis, nw, snap_js):
    """March the rectangle V in [V0, V0 + Nv h], U in [U0, U0 + Nu h] (data UL on V = V0, UR on U = U0) as hs.march_K does, restricted to the past set Pi = -V U >= -Pcap (V, U > 0 part),
    and keep only
      hp[:, m, j]    the column j at the nodes i_h - 2 + m (U = (m - 2) h, m = 2 is the horizon U = 0), all j;
      win[j, k, m]   Y / R at the nodes ilo[j, k] + m (m = 0 .. 2 nw - 1) of the column j around U* = -Pis[k] / V_j (ilo = -1: the window leaves the marched region or V_j <= 0);
      snaps[n]       the whole column snap_js[n]."""
    Nu = UL.shape[1] - 1
    Nv = UR.shape[1] - 1
    K = len(Pis)
    U1 = U0 + Nu * h
    hp = np.zeros((3, 5, Nv + 1))
    win = np.full((Nv + 1, K, 2 * nw), np.nan)
    ilo = np.full((Nv + 1, K), -1, dtype=np.int64)
    snaps = np.zeros((len(snap_js), 3, Nu + 1))
    prev = np.zeros((3, Nu + 1))
    cur = np.zeros((3, Nu + 1))
    B0 = np.zeros((3, 3))
    Bp = np.zeros((3, 3))
    Bm = np.zeros((3, 3))
    M = np.zeros((3, 3))
    rhs = np.zeros(3)
    out = np.zeros(3)
    Sv = np.zeros(3)
    Ev = np.zeros(3)
    Wv = np.zeros(3)
    for q in range(3):
        for i in range(Nu + 1):
            prev[q, i] = UL[q, i]
    itop_prev = Nu
    for j in range(Nv + 1):
        V = V0 + j * h
        if j > 0:
            Ucap = U1
            if V > 0.0 and Pcap / V < Ucap:
                Ucap = Pcap / V
            itop = min(Nu, int(math.floor((Ucap - U0) / h + 1e-9)))
            for q in range(3):
                cur[q, 0] = UR[q, j]
            guess = 0.0
            for i in range(itop):
                for q in range(3):
                    Sv[q] = prev[q, i]
                    Ev[q] = cur[q, i]
                    Wv[q] = prev[q, i + 1]
                # hs._cellK(Sv, Ev, Wv, V, U, ...): its (Sv, Ev, Wv) are (j, i), (j + 1, i), (j, i + 1) at the cell with corner (V_{j-1}, U_i)
                guess = hs._cellK(Sv, Ev, Wv, V - h, U0 + i * h, h, h, par, tab, guess, B0, Bp, Bm, M, rhs, out)
                for q in range(3):
                    cur[q, i + 1] = out[q]
            # the next column becomes the previous one
            for q in range(3):
                for i in range(Nu + 1):
                    prev[q, i] = cur[q, i]
            itop_prev = itop
        # prev is now the column j
        if i_h >= 2 and i_h + 2 <= itop_prev:
            for q in range(3):
                for m in range(5):
                    hp[q, m, j] = prev[q, i_h - 2 + m]
        if V > 0.0:
            for k in range(K):
                Us = -Pis[k] / V
                i0 = int(math.floor((Us - U0) / h)) - nw + 1
                if i0 >= 0 and i0 + 2 * nw - 1 <= itop_prev:
                    ilo[j, k] = i0
                    for m in range(2 * nw):
                        U = U0 + (i0 + m) * h
                        x, d, r, rV, rU, pV, pU, V0_, V1_, V2_, E, FV, FU, rhoV, rhoU = hs.jetsK(V, U, par, tab, 0.0)
                        win[j, k, m] = (prev[0, i0 + m] - x * r / 1.4142135623730951 * prev[2, i0 + m]) / r
        for n in range(len(snap_js)):
            if snap_js[n] == j:
                for q in range(3):
                    for i in range(Nu + 1):
                        snaps[n, q, i] = prev[q, i]
    return hp, win, ilo, snaps


# ----------------------------------------------------------------------------------------------------------------------
# settings, mutations, bookkeeping
# ----------------------------------------------------------------------------------------------------------------------
MUT = next((a.split("=", 1)[1] for a in _ARGV[1:] if a.startswith("--mutate=")), None)
ONLY = next((a.split("=", 1)[1].split(",") for a in _ARGV[1:] if a.startswith("--only=")), None)
CACHE = next((a.split("=", 1)[1] for a in _ARGV[1:] if a.startswith("--cache=")), None)
KNOWN_D = ("k_W", "k_J")
KNOWN_N = ("no_c2", "c2_sign", "no_c1", "c2p_V", "pi_off", "interp_near", "fit_semilog", "J_sign", "plant_alpha", "hp_sign", "g_nophi", "f0_W")
KNOWN = KNOWN_D + KNOWN_N
assert MUT is None or MUT in KNOWN, MUT

RECT0 = (-0.1, 3.0, -2.0, 0.06)             # #302's rectangle (V0, V1, U0, U1); V1 is replaced by VMAX
VMAX = 24.0
PCAP = 0.30                                 # the past set Pi >= -0.30 (the interior ends at Pi = -0.37, rho -> -oo)
PIS = (-0.1, 0.5, 1.0)                      # the three fixed-rho curves: Pi = -V U; rho = -0.3114, 0.6819, 1.0
PIS_ALL = PIS + (2.0,)                      # Pi = 2 (rho = 1.3698) is also recorded
HL = (0.00125, 0.000625, 0.0003125)
VA, VB = 8.0, 24.0                          # the fit window in V = e^{kappa v}
BETA = 2e-3                                 # amplitude of the exact modes: c2 = beta/2 = 1e-3, c2' = -1e-3
RECT_G = (0.05, VMAX, -2.0, 0.06)           # the exact growing mode: V >= 0.05
RECT_D = (-0.06, VMAX, -2.0, -0.05)         # its time reverse: U <= -0.05
UNIT = 1e-3
_KC = {}
RESULTS = []


def kchart():
    if "kc" not in _KC:
        _KC["kc"] = hs.KChart(C, pi_max=400.0)
    return _KC["kc"]


def rho_R(Pi):
    """rho and R at Pi (the chart at V = 1, U = -Pi)."""
    kc = kchart()
    x, d, r = hs.jetsK(1.0, -Pi, kc.par, kc.tab, 0.0)[:3]
    return float(x), float(r)


def sci(x):
    return f"{x:.3e}"


def want(cid):
    return ONLY is None or cid in ONLY


def check(cid, text, ok, detail=""):
    import re
    if re.search(r"(?<![a-z])(nan|-?inf)(?![a-z])", str(detail).lower()):
        ok = False
    RESULTS.append((cid, bool(ok)))
    print(f"[{'pass' if ok else 'FAIL'}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def rich(a, b):
    return (4.0 * b - a) / 3.0


# ----------------------------------------------------------------------------------------------------------------------
# readouts of a streamed run
# ----------------------------------------------------------------------------------------------------------------------
def stream_run(UL, UR, rect, h, Pis, centers, pcap=PCAP, nw=NW):
    """UL, UR: ray data of rect = (V0, V1, U0, U1). centers: [(name, V, U)]: 5 x 5 patches of dr around those nodes. Returns the readout dict (arrays only)."""
    V0, V1, U0, U1 = rect
    kc = kchart()
    ih = int(round(-U0 / h)) if U0 < -1e-9 and U1 > 1e-9 else -1
    cols = {}
    for nm, Vc, Uc in centers:
        jc, ic = int(round((Vc - V0) / h)), int(round((Uc - U0) / h))
        assert abs(V0 + jc * h - Vc) < 1e-9 and abs(U0 + ic * h - Uc) < 1e-9, (nm, Vc, Uc)
        cols[nm] = (jc, ic)
    js = sorted({j + a for jc, _ in cols.values() for j in (jc,) for a in range(-2, 3)})
    t0 = time.time()
    hp, win, ilo, snaps = march_stream(UL, UR, V0, U0, kc.par, kc.tab, h, pcap, ih, np.asarray(Pis, float), nw, np.array(js, dtype=np.int64))
    Nv = UR.shape[1] - 1
    ro = dict(h=np.array(h), V0=np.array(V0), U0=np.array(U0), hp=hp, win=win, ilo=ilo, wall=np.array(time.time() - t0), Pis=np.asarray(Pis, float), nw=np.array(nw))
    for nm, (jc, ic) in cols.items():
        ro["patch_" + nm] = np.array([[snaps[js.index(jc + a), 0, ic + b] for b in range(-2, 3)] for a in range(-2, 3)])     # f5[a + 2, b + 2] = dr(V + a h, U + b h)
    ro["V"] = V0 + h * np.arange(Nv + 1)
    return ro


def cached(key, fn):
    """readouts of evolved rectangles are stored in CACHE (mutations of the readings reuse them); no mutation changes an evolution."""
    import hashlib
    if CACHE:
        os.makedirs(CACHE, exist_ok=True)
        path = os.path.join(CACHE, hashlib.sha1(key.encode()).hexdigest()[:16] + ".npz")
        if os.path.exists(path):
            z = np.load(path, allow_pickle=False)
            return {k: z[k] for k in z.files}
        out = fn()
        np.savez(path + ".tmp.npz", **out)
        os.replace(path + ".tmp.npz", path)
        return out
    return fn()


_MEM = {}


def run_data(spec, h, tag):
    """evolved class / random datum of #248, #302 on RECT0 with V1 = VMAX, readouts at S and around the fixed-Pi curves."""
    key = f"data|{tag}|{sorted(spec.items())}|{h}|{VMAX}|{PCAP}|{PIS_ALL}|{NW}"
    if key not in _MEM:
        def go():
            rect = (RECT0[0], VMAX, RECT0[2], RECT0[3])
            su = hs.KSetup(C, h, *rect, spec, kc=kchart())
            return stream_run(su.UL, su.UR, rect, h, PIS_ALL, [("S", 0.0, 0.0)])
        _MEM[key] = cached(key, go)
    return _MEM[key]


def run_exact(kind, h):
    """the exact growing mode (kind 'g', V >= 0.05) or its time reverse (kind 'd', U <= -0.05) as ray data of the Kruskal solver (#302), streamed to V = VMAX."""
    key = f"exact|{kind}|{h}|{VMAX}|{PCAP}|{PIS_ALL}|{NW}"
    if key not in _MEM:
        def go():
            if kind == "g":
                rect = RECT_G
                UL, UR = hb.growing_ray_data(rect, h, BETA)
                cen = [("Hp05", 0.5, 0.0), ("Hp1", 1.0, 0.0)]
            else:
                rect = RECT_D
                UL, UR = hb.decaying_ray_data(rect, h, BETA)
                cen = [("Hm05", 0.0, -0.5), ("Hm1", 0.0, -1.0)]
            return stream_run(UL, UR, rect, h, PIS_ALL, cen)
        _MEM[key] = cached(key, go)
    return _MEM[key]


def c_read(ro, kind):
    """(c1, c2, c2') by the slice functionals: kind 'S' at S; 'g' on H+ at V0 = 0.5, 1 (c2' = 0, c1 = dr - V0 c2); 'd' on H- at U0 = -0.5, -1 (c2 = 0, c1 = dr - U0 c2')."""
    h = float(ro["h"])
    if kind == "S":
        f = ro["patch_S"]
        c2, c2p = hb.d_patch(f, h)
        return float(f[2, 2]), float(c2), float(c2p)
    if kind == "g":
        out = []
        for nm, V0 in (("Hp05", 0.5), ("Hp1", 1.0)):
            f = ro["patch_" + nm]
            c2 = hb.d_patch(f, h)[0]
            out.append((float(f[2, 2] - V0 * c2), float(c2), 0.0))
        return tuple(np.mean(out, axis=0))
    out = []
    for nm, U0 in (("Hm05", -0.5), ("Hm1", -1.0)):
        f = ro["patch_" + nm]
        c2p = hb.d_patch(f, h)[1]
        out.append((float(f[2, 2] - U0 * c2p), 0.0, float(c2p)))
    return tuple(np.mean(out, axis=0))


def apply_c_mutation(cs):
    c1, c2, c2p = cs
    if MUT == "no_c2":
        c2 = 0.0
    if MUT == "c2_sign":
        c2 = -c2
    if MUT == "no_c1":
        c1 = 0.0
    return c1, c2, c2p


def lagr(t, n):
    """weights of the degree n - 1 Lagrange interpolation on the nodes 0 .. n - 1 at t (array)."""
    t = np.asarray(t, float)
    w = np.ones((len(t), n))
    for m in range(n):
        for l in range(n):
            if l != m:
                w[:, m] *= (t - l) / (m - l)
    return w


def series(ro, k, cs, plant=None):
    """(j, V_j, w_rem(V_j)) with w_rem = Y_rem / R at the fixed Pi = Pis[k]: interpolate (8-point Lagrange in U) the node residuals Y/R - (c1 + c2 V + c2' U) to U* = -Pi / V.
    plant = (eps, alpha): add the planted Y/R = eps Pi V^alpha (Pi = -V U at the node; vanishes on both horizons) to every node before the residual."""
    c1, c2, c2p = apply_c_mutation(cs)
    h, V0, U0, nw = float(ro["h"]), float(ro["V0"]), float(ro["U0"]), int(ro["nw"])
    ilo, win, V = ro["ilo"][:, k], ro["win"][:, k, :], ro["V"]
    ok = (ilo >= 0) & (V > 0.0)
    j = np.where(ok)[0]
    Vj, i0 = V[j], ilo[j]
    Pi = float(ro["Pis"][k]) * (1.05 if MUT == "pi_off" else 1.0)
    Us = -Pi / Vj
    nodesU = U0 + (i0[:, None] + np.arange(2 * nw)[None, :]) * h
    u = win[j].copy()
    if plant is not None:
        eps, alpha = plant
        if MUT == "plant_alpha":
            alpha = alpha + 0.1
        u = u + eps * (-Vj[:, None] * nodesU) * Vj[:, None] ** alpha
    lin = c1 + c2 * Vj[:, None] + c2p * (Vj[:, None] if MUT == "c2p_V" else nodesU)
    resid = u - lin
    t = (Us - (U0 + i0 * h)) / h
    if MUT == "interp_near":
        near = np.clip(np.rint(t).astype(int), 0, 2 * nw - 1)
        return j, Vj, resid[np.arange(len(j)), near]
    w = lagr(t, 2 * nw)
    return j, Vj, np.sum(w * resid, axis=1)


def rich_series(sc_, sf_):
    """Richardson of the series at h (sc_) and h/2 (sf_) on their common nodes (fine index = 2 coarse index)."""
    jc, Vc, wc = sc_
    jf, Vf, wf = sf_
    common, ic, i_f = np.intersect1d(2 * jc, jf, return_indices=True)
    return common // 2, Vc[ic], (4.0 * wf[i_f] - wc[ic]) / 3.0


# ----------------------------------------------------------------------------------------------------------------------
# the fits
# ----------------------------------------------------------------------------------------------------------------------
def fit_power(V, y, Va=VA, Vb=VB):
    """least squares ln|y| = alpha ln V + b over Va <= V <= Vb (mutation fit_semilog: against V). Returns (alpha, std error, amplitude, sign constant)."""
    m = (V >= Va) & (V <= Vb)
    X, Y = V[m], y[m]
    sgn = bool(np.all(Y > 0) or np.all(Y < 0))
    x = X if MUT == "fit_semilog" else np.log(X)
    ly = np.log(np.maximum(np.abs(Y), 1e-300))
    A_ = np.vstack([x, np.ones_like(x)]).T
    co, res, _, _ = np.linalg.lstsq(A_, ly, rcond=None)
    n = len(x)
    s2 = float(np.sum((A_ @ co - ly) ** 2) / max(n - 2, 1))
    se = math.sqrt(s2 / max(float(np.sum((x - x.mean()) ** 2)), 1e-300))
    return float(co[0]), se, float(math.exp(co[1])), sgn


def fit_poly3(V, y, Va=VA, Vb=VB):
    """y = a V + b + d / V over Va <= V <= Vb: (a, b, d, rms residual / rms y)."""
    m = (V >= Va) & (V <= Vb)
    X, Y = V[m], y[m]
    A_ = np.vstack([X, np.ones_like(X), 1.0 / X]).T
    co = np.linalg.lstsq(A_, Y, rcond=None)[0]
    return float(co[0]), float(co[1]), float(co[2]), float(np.sqrt(np.mean((A_ @ co - Y) ** 2)) / max(np.sqrt(np.mean(Y ** 2)), 1e-300))


# ----------------------------------------------------------------------------------------------------------------------
# the closed form of #262: the profile of the growing solution that is outgoing at infinity, and the regular static solution
# ----------------------------------------------------------------------------------------------------------------------
I0_262 = 1.0149036099284089                 # #262: I0 = int_0^oo 2 K rho^2 e^{-2 kappa x} / (A R^2) drho
_CF = {}


def _A(r):
    return hb.sc.A_(r)


def _Ap(r):
    return hb.sc.Ap_(r)


def closed_form():
    """H(rho) = int_0^rho (K/A - 1/(2r)) dr and J(rho) = int_0^rho 2 K r e^{-2H} / (A R^2) dr (e^{2 kappa x} = rho e^{2H}, #262's normalization), by DOP853 from |rho| = 1e-4 with the series
    H = rho/(4K) + rho^2 (1/(4K)) (1/(2K) - 2K/3)/2, J = rho + O(rho^3); dense output on [-2.9, 400]. Mutation k_J: e^{+2H}."""
    if "J" in _CF:
        return _CF
    from scipy.integrate import solve_ivp
    K, r0 = KAPPA, 1e-4
    sg = -1.0 if MUT == "k_J" else 1.0

    def rhs(r, y):
        A = _A(r)
        return [(2 * K * r - A) / (2 * r * A), 2 * K * r * math.exp(-sg * 2 * y[0]) / (A * (1 + r * r))]

    sols = []
    for s, rend in ((1.0, 400.0), (-1.0, -2.9)):
        rr = s * r0
        H0 = rr / (4 * K) + 0.5 * rr * rr * (1 / (4 * K)) * (1 / (2 * K) - 2 * K / 3)
        sols.append(solve_ivp(rhs, [rr, rend], [H0, rr], method="DOP853", rtol=1e-13, atol=1e-18, dense_output=True))
    _CF["pos"], _CF["neg"] = sols
    _CF["I0"] = float(sols[0].y[1, -1])

    def J(r):
        r = float(r)
        if abs(r) < r0:
            return r
        return float((_CF["pos"] if r > 0 else _CF["neg"]).sol(r)[1])

    def H(r):
        r = float(r)
        if abs(r) < r0:
            return r / (4 * K)
        return float((_CF["pos"] if r > 0 else _CF["neg"]).sol(r)[0])

    _CF["J"], _CF["H"] = J, H
    return _CF


def Q_pred(rho):
    """the late-time profile of Y_rem / (c2 V R): Y_inf / Y_ref - 1 = -J(rho) / I0 (#262: Y_inf = Y_ref (1 - J / I0))."""
    cf = closed_form()
    return (+1.0 if MUT == "J_sign" else -1.0) * cf["J"](rho) / cf["I0"]


def ode_resid(rho, y, yp, ypp, W):
    """residual of Y'' = P Y' + Q Y (#262) for Y = E y, E = e^{kappa x} (E'/E = K/A): returns the residual over E and a scale."""
    K = KAPPA
    A, Ap = _A(rho), _Ap(rho)
    R2 = 1 + rho * rho
    P = 2 / rho - Ap / A
    Q = 1 / R2 ** 2 - Ap / (A * rho * R2) - 2 / (A * R2) + W * W / A ** 2
    e1 = K / A
    e1p = -K * Ap / A ** 2
    Y, Yp, Ypp = y, yp + e1 * y, ypp + 2 * e1 * yp + (e1p + e1 * e1) * y
    return Ypp - P * Yp - Q * Y, abs(Ypp) + abs(P * Yp) + abs(Q * Y)


def static_profile(rhos):
    """the solution of #262's equation at W = 0 that vanishes linearly at the horizon (Y = rho + ..., the free coefficient y_11 of a VU term), by DOP853 from |rho| = 1e-4."""
    from scipy.integrate import solve_ivp
    K = KAPPA

    def rhs(r, y):
        A, Ap = _A(r), _Ap(r)
        R2 = 1 + r * r
        P = 2 / r - Ap / A
        Q = 1 / R2 ** 2 - Ap / (A * r * R2) - 2 / (A * R2) + (0.05 ** 2 / A ** 2 if MUT == "f0_W" else 0.0)
        return [y[1], P * y[1] + Q * y[0]]

    out = {}
    for s, rend in ((1.0, 6.0), (-1.0, -0.95)):
        sol = solve_ivp(rhs, [s * 1e-4, rend], [s * 1e-4, 1.0], method="DOP853", rtol=1e-12, atol=1e-18, dense_output=True)
        for r in rhos:
            if r * s > 0:
                out[r] = float(sol.sol(r)[0])
    return out


def G0():
    """G0 = d rho / d Pi at S = 2 F_Pi(0) / pi (full_slice's series at s = 1e-6)."""
    return 2.0 * hb.fs.bgc(C).at(1e-6)["Fp"] / math.pi


def horizon_series(ro):
    """(j, V_j, d_U Y(V_j, 0)): the transverse slope of Y on H+ from the five nodes U = -2h .. 2h of the column (4th-order central d_U dr) and
    d_U Y = d_U dr + (G0 V / sqrt2) dphi at U = 0 (g = rho R / sqrt2 vanishes there, g_U = -G0 V / sqrt2). Independent of the fixed-Pi windows."""
    h, V, hp = float(ro["h"]), ro["V"], ro["hp"]
    drU = (-hp[0, 4] + 8 * hp[0, 3] - 8 * hp[0, 1] + hp[0, 0]) / (12.0 * h)
    if MUT == "hp_sign":
        drU = -drU
    YU = drU + (0.0 if MUT == "g_nophi" else 1.0) * G0() * V / SQ2 * hp[2, 2]
    j = np.where(V >= 1.0)[0]
    return j, V[j], YU[j]


# ----------------------------------------------------------------------------------------------------------------------
# data: classes of #248, #302 and random data
# ----------------------------------------------------------------------------------------------------------------------
CLASS_NAMES = ("A0", "A1R", "A1L", "L", "G", "Lm", "Rb", "Gb")
NRAND = 4


def data_items():
    rng = np.random.default_rng(3221)
    items = [(n, hb.CLASSES[n]) for n in CLASS_NAMES]
    items += [(f"rand{k}", hb._rand_spec(rng)) for k in range(NRAND)]
    return items


def datum_series(spec, tag, hl=HL):
    """per h: the readout, (c1, c2, c2'), and the series at the fixed Pi of PIS_ALL."""
    out = []
    for h in hl:
        ro = run_data(spec, h, tag)
        cs = c_read(ro, "S")
        out.append((ro, cs, [series(ro, k, cs) for k in range(len(PIS_ALL))]))
    return out


def analyse(hs_runs, k, va=VA, vb=VB):
    """for the fixed-Pi curve k: Richardson series of the two finest pairs, the fits, c2 (Richardson of the two finest)."""
    s = [r[2][k] for r in hs_runs]
    R1, R2 = rich_series(s[0], s[1]), rich_series(s[1], s[2])
    c2 = rich(hs_runs[1][1][1], hs_runs[2][1][1])
    out = dict(c2=c2)
    for nm, ser_ in (("h2", s[1]), ("h3", s[2]), ("R1", R1), ("R2", R2)):
        a, se, amp, sg = fit_power(ser_[1], ser_[2], va, vb)
        p3 = fit_poly3(ser_[1], ser_[2], va, vb)
        out[nm] = dict(alpha=a, se=se, amp=amp, sign=sg, a=p3[0], b=p3[1], d=p3[2], res=p3[3])
    out["V"], out["w"] = R2[1], R2[2]
    return out


# ----------------------------------------------------------------------------------------------------------------------
# derive: K1 the closed form of #262
# ----------------------------------------------------------------------------------------------------------------------
def derive_K1():
    """K1a: I0 of the quadrature against #262's value. K1b: Y_ref = R e^{kappa x} and Y_b = Y_ref J solve #262's equation at W = kappa (residual over the size of its terms). K1c: Y_inf = Y_ref (1 - J / I0)
    is outgoing: Y_inf / Y_ref <= 1e-6 at rho = 20."""
    cf = closed_form()
    K = KAPPA
    ok_a = abs(cf["I0"] - I0_262) <= 1e-8 * I0_262
    check("K1a", "I0 = int_0^oo 2K r e^{-2H} / (A R^2) dr (H from rho = 1e-4, DOP853 1e-13) equals #262's 1.0149036099284089 to 1e-8 relative", ok_a, f"I0 = {cf['I0']:.12g}, relative {cf['I0'] / I0_262 - 1:.1e}")
    W = K * (1.01 if MUT == "k_W" else 1.0)
    sg = -1.0 if MUT == "k_J" else 1.0
    worst_ref, worst_b = 0.0, 0.0
    for rho in (-0.8, -0.3, 0.4, 1.0, 2.0, 5.0):
        R = math.sqrt(1 + rho * rho)
        Rp, Rpp = rho / R, 1 / R ** 3
        r1, s1 = ode_resid(rho, R, Rp, Rpp, W)
        worst_ref = max(worst_ref, abs(r1) / s1)
        A, Ap = _A(rho), _Ap(rho)
        J = cf["J"](rho)
        w = 2 * K * rho * math.exp(-2 * sg * cf["H"](rho)) / (A * (1 + rho * rho))
        # w'/w = 1/r - 2 sg H' - A'/A - 2r/R^2 with H' = K/A - 1/(2r)
        wp = w * (1 / rho - 2 * sg * (K / A - 0.5 / rho) - Ap / A - 2 * rho / (1 + rho * rho))
        y, yp, ypp = R * J, Rp * J + R * w, Rpp * J + 2 * Rp * w + R * wp
        r2, s2 = ode_resid(rho, y, yp, ypp, W)
        worst_b = max(worst_b, abs(r2) / s2)
    check("K1b", "Y_ref = R e^{kappa x} and Y_b = Y_ref J solve #262's equation Y'' = P Y' + Q Y at W = kappa (rho = -0.8 ... 5): residual over the size of its terms <= 1e-9", worst_ref <= 1e-9 and worst_b <= 1e-9,
          f"Y_ref {worst_ref:.1e}, Y_b {worst_b:.1e}")
    tail = 1.0 - cf["J"](20.0) / cf["I0"]
    check("K1c", "Y_inf / Y_ref = 1 - J(rho) / I0 is the outgoing solution: <= 1e-6 at rho = 20", abs(tail) <= 1e-6, f"{tail:.2e}")


def run_derive():
    t0 = time.time()
    derive_K1()
    bad = [c for c, ok in RESULTS if not ok]
    print(f"derive: {len(RESULTS) - len(bad)} of {len(RESULTS)} pass ({time.time() - t0:.0f} s)")
    return 1 if bad else 0


# ----------------------------------------------------------------------------------------------------------------------
# checks on the planted exact modes (known remainders)
# ----------------------------------------------------------------------------------------------------------------------
BAR_P1 = 1e-4          # planted zero: max |Y_rem / R| / (|c2| V) over the window at the finest h
EPS_P2 = 0.01          # P2a: the growing mode with 1 % of its c2 left in the remainder
BAR_RATIO = (3.4, 4.6) # the h^2 convergence of the truncation floor (last ratio)


def exact_series(kind, hl=HL, mod=None):
    """per h: series of the exact mode for the three Pi of PIS_ALL, with (c1, c2, c2') read on the horizon crossings; mod(cs) changes what is subtracted."""
    out = []
    for h in hl:
        ro = run_exact(kind, h)
        cs = c_read(ro, kind)
        if mod:
            cs = mod(cs)
        out.append((ro, cs, [series(ro, k, cs) for k in range(len(PIS_ALL))]))
    return out


def window_for(kind, Pi):
    return (VA, VB) if kind == "g" else (4.0, min(VB, 18.0 * Pi))


def check_P1():
    """planted zero: the exact growing mode has Y_rem = 0; what is left is the truncation floor."""
    runs = exact_series("g")
    worst, rt, out = [0.0, 0.0, 0.0], [], []
    for ih, r in enumerate(runs):
        for k in range(len(PIS)):
            j, V, w = r[2][k]
            m = (V >= VA) & (V <= VB)
            worst[ih] = max(worst[ih], float(np.max(np.abs(w[m]) / (abs(r[1][1]) * V[m]))))
    rt = [worst[0] / worst[1], worst[1] / worst[2]]
    c2s = [r[1][1] for r in runs]
    out = f"max |Y_rem/R| / (c2 V) over V in [8, 24], three Pi, at h = {', '.join(str(h) for h in HL)}: {', '.join(sci(x) for x in worst)}; ratios {rt[0]:.2f}, {rt[1]:.2f}; c2 {c2s[-1]:.9e} (exact {BETA / 2:.1e})"
    check("P1", f"planted zero: the exact growing mode (beta = 2e-3, c2 = 1e-3) through the stream marcher, c1, c2 by the slice functional at V0 = 0.5, 1: Y_rem / R <= {BAR_P1:.0e} c2 V at the finest h, "
          f"falling by {BAR_RATIO[0]} to {BAR_RATIO[1]} per halving", worst[2] <= BAR_P1 and BAR_RATIO[0] <= rt[1] <= BAR_RATIO[1], out)


def check_P2():
    """planted rates +kappa and -kappa: the same fit on remainders of known rate. (a) the growing mode with (1 - eps) c2 subtracted: Y_rem/R = eps c2 V, alpha = 1, a/c2 = eps. (b) the decaying mode with
    c2' U R left in the remainder: Y_rem/R = c2' U* = -c2' Pi / V, alpha = -1, amplitude |c2'| Pi."""
    outs, ok = [], True
    ra = exact_series("g", hl=HL[1:], mod=lambda cs: (cs[0], (1.0 - EPS_P2) * cs[1], cs[2]))
    for k in range(len(PIS)):
        R2 = rich_series(ra[0][2][k], ra[1][2][k])
        al, se, amp, sg = fit_power(R2[1], R2[2], VA, VB)
        a, b, d, res = fit_poly3(R2[1], R2[2], VA, VB)
        c2 = rich(ra[0][1][1] / (1 - EPS_P2), ra[1][1][1] / (1 - EPS_P2))
        good = sg and abs(al - 1.0) <= 0.01 and abs(a / c2 - EPS_P2) <= 0.02 * EPS_P2
        ok = ok and good
        outs.append(f"(a) Pi = {PIS[k]}: alpha {al:.5f} (se {se:.1e}), a/c2 {a / c2:.6f} (exact {EPS_P2})")
    rd = exact_series("d", hl=HL[1:], mod=lambda cs: (cs[0], 0.0, 0.0))
    for k in range(len(PIS)):
        Pi = PIS[k]
        if Pi <= 0:
            continue
        va, vb = window_for("d", Pi)
        R2 = rich_series(rd[0][2][k], rd[1][2][k])
        al, se, amp, sg = fit_power(R2[1], R2[2], va, vb)
        c2p = rich(rd[0][1][2], rd[1][1][2])
        good = sg and abs(al + 1.0) <= 0.01 and abs(amp / (abs(c2p) * Pi) - 1.0) <= 0.01
        ok = ok and good
        outs.append(f"(b) Pi = {Pi}: alpha {al:.5f} (se {se:.1e}), amplitude / (|c2'| Pi) {amp / (abs(c2p) * Pi):.6f}")
    check("P2", "planted rates through the same fit (Richardson of h = 0.000625, 0.0003125): (a) growing mode with 99 % of c2 subtracted: alpha = 1 +- 0.01, a/c2 = 0.01 +- 2 %; (b) decaying mode with c2' U R left in: "
          "alpha = -1 +- 0.01, amplitude |c2'| Pi +- 1 %", ok, "; ".join(outs))


def check_P3():
    """planted array-level remainders eps Pi V^alpha (vanishing on both horizons) added to the exact growing mode's windows: the same series and fit recover alpha and eps Pi."""
    ra = exact_series("g", hl=HL[1:])
    outs, ok = [], True
    for alpha in (2.0, 1.0, 0.0, -1.0):
        eps = 2e-4 / 16.0 ** alpha
        for k in range(len(PIS)):
            if PIS[k] <= 0:
                continue
            sers = [series(r[0], k, r[1], plant=(eps, alpha)) for r in ra]
            R2 = rich_series(sers[0], sers[1])
            al, se, amp, sg = fit_power(R2[1], R2[2], VA, VB)
            good = sg and abs(al - alpha) <= 0.01 and abs(amp / (eps * PIS[k]) - 1.0) <= 0.01
            ok = ok and good
            outs.append(f"alpha {alpha:+.0f}, Pi {PIS[k]}: {al:+.4f}, amplitude ratio {amp / (eps * PIS[k]):.4f}")
    check("P3", "planted remainders eps Pi V^alpha (alpha = 2, 1, 0, -1; 2e-4 at V = 16, Pi = 1) added at the nodes of the exact growing mode: the fit returns alpha +- 0.01 and the amplitude +- 1 %", ok, "; ".join(outs))


# ----------------------------------------------------------------------------------------------------------------------
# checks on the evolved data
# ----------------------------------------------------------------------------------------------------------------------
C2MIN = 1e-4            # data with |c2| below this (class Rb: the truncation seed only) have no c2 V R to compare with
BAR_Q = 1e-2            # |a / c2 - Q_pred| of the three-term fit
BAR_QC = 5e-3           # |Q(R1) - Q(R2)|: the two Richardson pairs
ALPHA_LO, ALPHA_HI = 0.8, 1.02
BAR_S = 0.05            # the static piece: spread of y1 over the three Pi


def datum_summary(name, spec):
    runs = datum_series(spec, name)
    c2 = rich(runs[1][1][1], runs[2][1][1])
    c2_1 = rich(runs[0][1][1], runs[1][1][1])
    c1 = rich(runs[1][1][0], runs[2][1][0])
    c2p = rich(runs[1][1][2], runs[2][1][2])
    out = dict(name=name, c1=c1, c2=c2, c2p=c2p, k=[])
    gs = [horizon_series(r[0]) for r in runs]
    gR = rich_series(gs[1], gs[2])
    m = (gR[1] >= VA) & (gR[1] <= VB)
    Vg, yg = gR[1][m], gR[2][m]
    Ag = np.vstack([Vg ** 2, Vg, np.ones_like(Vg)]).T
    cg = np.linalg.lstsq(Ag, yg, rcond=None)[0]
    out["gA"], out["gB"], out["gC"] = float(cg[0]), float(cg[1]), float(cg[2])
    out["gres"] = float(np.sqrt(np.mean((Ag @ cg - yg) ** 2)) / max(np.sqrt(np.mean(yg ** 2)), 1e-300))
    out["gA_ratio"] = out["gA"] / (c2 * G0() / I0_262) if abs(c2) >= C2MIN else float("nan")
    fs0 = static_profile([rho_R(Pi)[0] for Pi in PIS])
    for k in range(len(PIS)):
        an = analyse(runs, k)
        rho, R = rho_R(PIS[k])
        qp = Q_pred(rho)
        an["rho"], an["R"], an["Qpred"] = rho, R, qp
        an["Q2"], an["Q1"] = an["R2"]["a"] / c2, an["R1"]["a"] / c2_1
        an["y1"] = an["R2"]["b"] * R / fs0[rho]
        an["aerr"] = max(abs(an["R2"]["alpha"] - an["R1"]["alpha"]), abs(an["h2"]["alpha"] - an["h3"]["alpha"]), an["R2"]["se"])
        out["k"].append(an)
    return out


_SUM = {}


def summaries():
    if "s" not in _SUM:
        _SUM["s"] = [datum_summary(n, sp) for n, sp in data_items()]
    return _SUM["s"]


def check_D1():
    """the late growth at fixed rho: Y_rem / R = a V + b + d / V over V in [8, 24] (Richardson of the two finest pairs); a / c2 against #262's -J(rho) / I0, the window exponent of |Y_rem|."""
    ok, outs = True, []
    for s in summaries():
        if abs(s["c2"]) < C2MIN:
            continue
        dev = max(abs(a["Q2"] - a["Qpred"]) for a in s["k"])
        conv = max(abs(a["Q1"] - a["Q2"]) for a in s["k"])
        al = [a["R2"]["alpha"] for a in s["k"]]
        sg = all(a["R2"]["sign"] for a in s["k"])
        good = dev <= BAR_Q and conv <= BAR_QC and all(ALPHA_LO <= x <= ALPHA_HI for x in al) and sg
        ok = ok and good
        outs.append(f"{s['name']}: |a/c2 - Qpred| {dev:.1e}, pair change {conv:.1e}, alpha {', '.join(f'{x:.3f}' for x in al)}{'' if good else ' <-- FAIL'}")
    check("D1", f"smooth data, evolved to V = {VMAX:g}: at Pi = -0.1, 0.5, 1 (rho = -0.311, 0.682, 1) Y_rem / R = a V + b + d / V on V in [8, 24] (Richardson, h = 0.000625, 0.0003125): a / c2 equals #262's -J(rho) / I0 to {BAR_Q:g}, "
          f"changes by <= {BAR_QC:g} from the pair h = 0.00125, 0.000625, the window exponent alpha of |Y_rem| in [{ALPHA_LO}, {ALPHA_HI}] with one sign; data with |c2| >= {C2MIN:g}", ok, "; ".join(outs[:6]) + " ...")


def check_D2():
    """the V-independent part: b R / f0(rho) with f0 the solution of #262's equation at W = 0 that vanishes at the horizon: the same number at the three Pi."""
    ok, outs = True, []
    for s in summaries():
        ys = [a["y1"] for a in s["k"]]
        sp = (max(ys) - min(ys)) / max(abs(np.mean(ys)), 1e-300)
        if abs(s["c2"]) < C2MIN and max(abs(a["R2"]["b"]) for a in s["k"]) < 1e-9:
            continue
        good = sp <= BAR_S
        ok = ok and good
        outs.append(f"{s['name']}: y1 {', '.join(sci(x) for x in ys)}, spread {sp:.1e}")
    check("D2", f"the V-independent part b of the fit has the profile of the static solution f0 that vanishes at the horizon: y1 = b R / f0(rho) agrees at the three Pi to {BAR_S:g} (spread over mean)", ok, "; ".join(outs[:6]) + " ...")


BAR_A = 1e-2            # the horizon readout: |A / (c2 G0 / I0) - 1|


def check_D3():
    """independent readout on H+: d_U Y(V, 0) = A V^2 + B V + C over V in [8, 24] (Richardson of the two finest pairs): A = c2 G0 / I0 (Y_inf's profile: Y / R = c2 V (1 + Q), Q = -rho / I0 + ..., rho = -G0 V U)."""
    ok, outs = True, []
    for s in summaries():
        if abs(s["c2"]) < C2MIN:
            continue
        good = abs(s["gA_ratio"] - 1.0) <= BAR_A
        ok = ok and good
        outs.append(f"{s['name']}: A/(c2 G0/I0) {s['gA_ratio']:.5f} (fit residual {s['gres']:.1e})")
    check("D3", f"the transverse slope d_U Y on the future horizon (five nodes across U = 0, no fixed-Pi window) grows as A V^2 + B V + C with A = c2 G0 / I0 to {BAR_A:g}; data with |c2| >= {C2MIN:g}", ok, "; ".join(outs[:6]) + " ...")


CHECKS = [("K1", None), ("P1", check_P1), ("P2", check_P2), ("P3", check_P3), ("D1", check_D1), ("D2", check_D2), ("D3", check_D3)]
SCOPE = ("K1a", "K1b", "K1c", "P1", "P2", "P3", "D1", "D2", "D3")
DECLARED = {}


def run_checks():
    t0 = time.time()
    for cid, f in CHECKS:
        if f is None or not want(cid):
            continue
        t1 = time.time()
        f()
        print(f"    ({cid}: {time.time() - t1:.0f} s)", flush=True)
    bad = [c for c, ok in RESULTS if not ok]
    print(f"check: {len(RESULTS) - len(bad)} of {len(RESULTS)} pass ({time.time() - t0:.0f} s)")
    return 1 if bad else 0


def report_mutations(raw, prefix):
    bad = 0
    for name in KNOWN:
        path = os.path.join(raw, f"{prefix}{name}.log")
        if not os.path.exists(path):
            print(f"{name}: no log {path}")
            bad += 1
            continue
        fired = set()
        for line in open(path):
            if line.startswith("[FAIL] "):
                fired.add(line[7:].split(":")[0])
        inscope = {c for c in fired if c in SCOPE}
        ok = inscope == DECLARED.get(name, None)
        bad += 0 if ok else 1
        print(f"{name}: fired {sorted(inscope)}; declared {sorted(DECLARED.get(name, []))}; {'equal' if ok else 'DIFFERENT'}")
    print(f"{len(KNOWN) - bad} of {len(KNOWN)} fired sets equal the declared ones")
    return 1 if bad else 0


def report_table():
    return 0


if __name__ == "__main__":
    cmd = next((a for a in _ARGV[1:] if not a.startswith("--")), "check")
    if cmd == "derive":
        sys.exit(run_derive())
    if cmd == "check":
        sys.exit(run_checks())
    if cmd == "table":
        sys.exit(report_table())
    if cmd == "mutations":
        i = _ARGV.index("mutations")
        sys.exit(report_mutations(_ARGV[i + 1], _ARGV[i + 2]))
    sys.exit(f"unknown command {cmd}")

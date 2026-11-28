#!/usr/bin/env python3
"""The horizon basis (c1, c2, c2') of the c = -1 black universe: the decaying partner of #262's mode (#302).

At the bifurcation surface S the regular numbers of Cauchy data are c1 = Y(S) (static mass mode, #259), c2 = d_V dr(S) (growing mode Y = c2 V R(rho(Pi)), #262) and
c2' = d_U dr(S) (its time reverse Y = c2' U R(rho(Pi))). This script
  * derives that the exact modes are not smooth solutions through S (d_N dphi ~ -(pi/sqrt2) c2 ln|s|), and that this is a property of the solution, not of the gauge (S1, L1, L2);
  * puts the three exact modes through the Kruskal solver of #248 (P2 static, P3 growing, P4 decaying: exact double-null ray data built from #269's family) and reads
    c1, c2, c2' by the slice functionals and by the horizon fits;
  * tests the law Y_h(U) = c1 + c2' U on the past horizon V = 0 and the three-mode reproduction of Y_h on both horizons for evolved classes and random data (H1 - H3).

    horizon_basis.py derive [--mutate=NAME]                 sympy / mpmath checks (S-series)
    horizon_basis.py check [--only=P2,..] [--mutate=NAME] [--cache=DIR]   numerical checks; exit 1 when one fails
    horizon_basis.py table                                  the numbers of the note
    horizon_basis.py mutations DIR PREFIX                   fired sets of the stored logs DIR/PREFIX<name>.log against the declared ones

Reads horizon_sector.py, full_slice.py, static_c1.py, linear.py (read-only imports). Python: ~/Source/negative-mass/.venv/bin/python.
"""
import hashlib
import math
import os
import re
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
_ARGV = sys.argv
sys.argv = sys.argv[:1]                     # the imported modules parse --mutate / --only themselves
import horizon_sector as hs  # noqa: E402
import full_slice as fs  # noqa: E402
import static_c1 as sc  # noqa: E402
sys.argv = _ARGV

MUT = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mutate=")), None)
ONLY = next((a.split("=", 1)[1].split(",") for a in sys.argv[1:] if a.startswith("--only=")), None)
CACHE = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--cache=")), None)
KNOWN_D = ("s_zeta", "s_kappa", "s_noflip", "s_swapfn")
KNOWN_N = ("fn_swap", "fn_aT", "law_c1", "law_x", "law_off", "gauge_dF", "gauge_dphi", "alpha_off", "sc_off", "zeta_F")
KNOWN = KNOWN_D + KNOWN_N
assert MUT is None or MUT in KNOWN, MUT

C = -1.0
KAPPA = 2.0 / math.pi
SQ2 = math.sqrt(2.0)
RECT = hs.RECT                                          # (V0, V1, U0, U1) = (-0.05, 3, -2, 0.06)
RECT_H = (-0.1, 3.0, -2.0, 0.06)                        # data classes: the past strip V in [-0.1, 0] is wide enough for a pulse on the ray U = -2 to reach V = 0 (Pi = -0.2 at the corner)
RECT_G = (0.05, 3.05, -2.0, 0.06)                       # growing mode: V > 0 (the mode is smooth there, not at V = 0)
RECT_D = (-0.06, 2.0, -3.0, -0.05)                      # decaying mode: U < 0
HL = (0.005, 0.0025, 0.00125)
HF = (0.00125, 0.000625)                                 # the finest pair: Richardson of the readings of H2, H3, L1
BETA = 2e-3                                             # the mode amplitude of the planted cases: c2 = beta/2 = 1e-3 (growing), c2' = -1e-3 (decaying)
CLASSES = dict(hs.CLASSES)
CLASSES["Lm"] = dict(L=(1e-3, -1.0, 0.3))               # pulse on the ray V = V0 centred on U = -1 (in the past of the past horizon V = 0, away from S)
CLASSES["Rb"] = dict(R=(1e-3, -0.05, 0.03))             # pulse on the ray U = U0 centred in the strip V in [-0.1, 0] (in the past of the horizon point (0, U))
CLASSES["Gb"] = dict(L=(1e-3, -1.0, 0.3), R=(1e-3, -0.05, 0.03))
UNIT = 1e-3                                             # the data amplitude: every datum is specified at amplitude <= 1e-3
PTS = {"S": (0.0, 0.0), "Hp05": (0.5, 0.0), "Hp1": (1.0, 0.0), "Hp2": (2.0, 0.0), "Hm05": (0.0, -0.5), "Hm1": (0.0, -1.0), "Hm2": (0.0, -2.0)}


def sci(x):
    return f"{x:.3e}"


# ----------------------------------------------------------------------------------------------------------------------
# bookkeeping
# ----------------------------------------------------------------------------------------------------------------------
RESULTS = []


def want(cid):
    return ONLY is None or cid in ONLY


def check(cid, text, ok, detail=""):
    if re.search(r"(?<![a-z])(nan|-?inf)(?![a-z])", str(detail).lower()):
        ok = False
    RESULTS.append((cid, bool(ok)))
    print(f"[{'pass' if ok else 'FAIL'}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def rich(a, b):
    return (4.0 * b - a) / 3.0


def ratios(vals):
    return [abs(vals[k]) / max(abs(vals[k + 1]), 1e-300) for k in range(len(vals) - 1)]


# ----------------------------------------------------------------------------------------------------------------------
# readings: slice functionals on 5 x 5 patches, horizon fits
# ----------------------------------------------------------------------------------------------------------------------
def d_patch(f5, h):
    """(d_V, d_U) at the centre node of a 5 x 5 patch f5[a + 2, b + 2] = f(j + a, i + b): 4th-order T = d_V - d_U along the slice nodes (1, -1) and
    N = d_V + d_U along (1, 1) (hs.fd4); d_V = (T + N)/2, d_U = (N - T)/2. Mutations fn_swap (d_U reads d_V) and fn_aT (tangential half only)."""
    T = (-f5[4, 0] + 8 * f5[3, 1] - 8 * f5[1, 3] + f5[0, 4]) / (12.0 * h)
    N = (-f5[4, 4] + 8 * f5[3, 3] - 8 * f5[1, 1] + f5[0, 0]) / (12.0 * h)
    dV = 0.5 * (T + N)
    dU = 0.5 * (N - T)
    if MUT == "fn_swap":
        dU = dV
    if MUT == "fn_aT":
        dU = -T
    return dV, dU


def fit_h(vals, x):
    """horizon_fit of a horizon line: vals = c1 + c2 x. Mutation law_x fits against -U (V and U swapped through the time reversal (V, U) -> (-U, -V) without the sign of c2')."""
    if MUT == "law_x":
        x = -x
    return hs.horizon_fit(vals, x)


def readout(U, V0, U0, h, V1=None):
    """small arrays of an evolved rectangle: the horizon lines, 5 x 5 patches of dr at PTS, and the band around the slice V + U = 0 (s from 0 to 0.6)."""
    Nv, Nu = U.shape[1] - 1, U.shape[2] - 1
    V = V0 + h * np.arange(Nv + 1)
    Ua = U0 + h * np.arange(Nu + 1)
    out = dict(V=V, Ua=Ua, h=np.array(h))
    ih = int(round(-U0 / h)) if (U0 < -1e-9 and Ua[-1] > 1e-9) else None            # the node U = 0 (future horizon line)
    jS = int(round(-V0 / h)) if (V0 < -1e-9 and V[-1] > 1e-9) else None              # the node V = 0 (past horizon line)
    assert ih is None or abs(U0 + ih * h) < 1e-9
    assert jS is None or abs(V0 + jS * h) < 1e-9
    if ih is not None:
        out["Hp"] = U[0, :, ih].copy()
    if jS is not None:
        out["Hm"] = U[0, jS, :].copy()
        out["Hm_off"] = U[0, jS + int(round(0.025 / h)), :].copy()                 # the line V = 0.025 (not a horizon): mutation law_off
    for name, (Vp, Up) in PTS.items():
        j, i = int(round((Vp - V0) / h)), int(round((Up - U0) / h))
        if 2 <= j <= Nv - 2 and 2 <= i <= Nu - 2 and abs(V0 + j * h - Vp) < 1e-9 and abs(U0 + i * h - Up) < 1e-9:
            out["patch_" + name] = U[0, j - 2: j + 3, i - 2: i + 3].copy()
    if V0 < 0 < V[-1] and U0 < 0 < Ua[-1]:
        ks = [k for k in range(int(round((-0.05 - V0) / h)), int(round((0.6 - V0) / h)) + 1)]
        band = np.zeros((3, len(ks), 5))
        sv = np.zeros(len(ks))
        for a, k in enumerate(ks):
            j = k
            i = int(round((-(V0 + j * h) - U0) / h))
            sv[a] = V0 + j * h
            for d in range(-2, 3):
                band[:, a, d + 2] = U[:, j + d, i + d]
        out["band"], out["band_s"] = band, sv
    return out


def cached(key, fn):
    """store the readout of a run in CACHE (mutations of the readings reuse the evolutions); the evolution itself has no mutation."""
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


def run_spec(spec, h, rect=RECT_H, tag=""):
    key = f"spec|{sorted(spec.items())}|{h}|{rect}|{tag}"
    if key not in _MEM:
        def go():
            kc = hs.kchart(C)
            su = hs.KSetup(C, h, *rect, spec, kc=kc)
            U = su.run()
            return readout(U, rect[0], rect[2], h)
        _MEM[key] = cached(key, go)
    return _MEM[key]


def run_arrays(UL, UR, rect, h, tag):
    key = f"arr|{tag}|{h}|{rect}"
    if key not in _MEM:
        def go():
            kc = hs.kchart(C)
            U = hs.march_K(UL, UR, rect[0], rect[2], kc.par, kc.tab, h, h)
            ro = readout(U, rect[0], rect[2], h)
            return ro, U
        if CACHE:
            def go2():
                ro, U = go()
                ro["Yerr_pending"] = np.array(0.0)
                return ro
            _MEM[key] = cached(key, go2)
        else:
            _MEM[key] = go()[0]
    return _MEM[key]


def readings(ro, h):
    """c1 at S, slice functionals at S, fits on both horizons."""
    out = {}
    if "patch_S" in ro:
        f = ro["patch_S"]
        out["c1S"] = float(f[2, 2])
        dV, dU = d_patch(f, h)
        out["c2S"], out["c2pS"] = float(dV), float(dU)
    for nm in ("Hp05", "Hp1", "Hp2"):
        if "patch_" + nm in ro:
            out["c2_" + nm] = float(d_patch(ro["patch_" + nm], h)[0])
    for nm in ("Hm05", "Hm1", "Hm2"):
        if "patch_" + nm in ro:
            out["c2p_" + nm] = float(d_patch(ro["patch_" + nm], h)[1])
    if "Hp" in ro:
        c1, c2, res = hs.horizon_fit(ro["Hp"], ro["V"])
        out["fitp"] = (float(c1), float(c2), float(res))
    if "Hm" in ro:
        c1, c2, res = fit_h(ro["Hm_off"] if MUT == "law_off" else ro["Hm"], ro["Ua"])
        out["fitm"] = (float(c1), float(c2), float(res))
    return out


# ----------------------------------------------------------------------------------------------------------------------
# the exact growing mode in a double-null gauge (from #269's ingoing-null family) and its time reverse
# ----------------------------------------------------------------------------------------------------------------------
def _chart():
    return hs.kchart(C)


def _jets(V, U):
    kc = _chart()
    j = hs.jetsK(V, U, kc.par, kc.tab, 0.0)
    return j[0], j[10]                                   # rho, e^{-F}


def _ode_rhs(Vv, y, U):
    """y = (xi^u, d_U xi^u) at fixed U along V (beta = 1): xi^u_V = W/kappa, (d_U xi^u)_V = -V W_rho G/kappa, W = (1 + kappa rho)/A, G = kappa e^{-F}."""
    rho, E = _jets(Vv, U)
    A, Ap = sc.A_(rho), sc.Ap_(rho)
    W = (1.0 + KAPPA * rho) / A
    Wr = (KAPPA * A - (1.0 + KAPPA * rho) * Ap) / (A * A)
    G = KAPPA * E
    return [W / KAPPA, -Vv * Wr * G / KAPPA]


def xi_line(Vs, U, Vb):
    """(xi^u, d_U xi^u) at the points Vs (sorted or not, all > 0) on the line U = const, with xi^u = 0 on V = Vb."""
    Vs = np.asarray(Vs, float)
    out = np.zeros((2, len(Vs)))
    for sgn in (+1, -1):
        sel = np.where(sgn * (Vs - Vb) > 1e-14)[0]
        if len(sel) == 0:
            continue
        order = sel[np.argsort(sgn * Vs[sel])]
        sol = solve_ivp(_ode_rhs, [Vb, Vs[order[-1]]], [0.0, 0.0], args=(U,), method="DOP853", rtol=1e-13, atol=1e-18, t_eval=Vs[order])
        out[:, order] = sol.y
    return out


def g_fields(V, U, xi, dxi):
    """(dr, dF, dphi) of the exact growing mode (beta = 1) in the double-null gauge from xi^u = xi, d_U xi^u = dxi at (V, U)."""
    rho, E = _jets(V, U)
    R = math.sqrt(1.0 + rho * rho)
    A, Ap = sc.A_(rho), sc.Ap_(rho)
    dr = 0.5 * R * V - xi * A * rho / (2.0 * R)
    dphi = 0.0 if MUT == "gauge_dphi" else -xi * A / (SQ2 * R * R)
    dF = 0.5 * Ap * xi + (0.0 if MUT == "gauge_dF" else KAPPA * U * dxi)
    return dr, dF, dphi


def growing_ray_data(rect, h, beta):
    """data on V = V0 (UL) and U = U0 (UR) of the growing mode, gauge xi^u = 0 on V = V0."""
    V0, V1, U0, U1 = rect
    Nv, Nu = int(round((V1 - V0) / h)), int(round((U1 - U0) / h))
    UL, UR = np.zeros((3, Nu + 1)), np.zeros((3, Nv + 1))
    for i in range(Nu + 1):
        UL[:, i] = beta * np.array(g_fields(V0, U0 + i * h, 0.0, 0.0))
    Vs = V0 + h * np.arange(Nv + 1)
    xi = xi_line(Vs, U0, V0)
    for j in range(Nv + 1):
        UR[:, j] = beta * np.array(g_fields(Vs[j], U0, xi[0, j], xi[1, j]))
    return UL, UR


def decaying_ray_data(rect, h, beta):
    """time reverse (V, U) -> (-U, -V) of the growing mode, whose gauge is xi^u = 0 on V' = Vb = -U0 (the data ray U = U0 is the image of that ray)."""
    V0, V1, U0, U1 = rect
    Nv, Nu = int(round((V1 - V0) / h)), int(round((U1 - U0) / h))
    UL, UR = np.zeros((3, Nu + 1)), np.zeros((3, Nv + 1))
    Vb = -U0
    for j in range(Nv + 1):                              # ray U = U0: V' = -U0, U' = -V
        UR[:, j] = beta * np.array(g_fields(Vb, -(V0 + j * h), 0.0, 0.0))
    Vp = -(U0 + h * np.arange(Nu + 1))                   # ray V = V0: U' = -V0, V' = -U
    xi = xi_line(Vp, -V0, Vb)
    for i in range(Nu + 1):
        UL[:, i] = beta * np.array(g_fields(Vp[i], -V0, xi[0, i], xi[1, i]))
    return UL, UR


def y_exact(V, U, amp):
    rho, _ = _jets(V, U)
    return amp * V * math.sqrt(1.0 + rho * rho)


def y_num_err(U, rect, h, amp, kind):
    """max |Y_num - Y_exact| / max |Y_exact| on a sub-grid, Y = dr - (rho R/sqrt2) dphi; kind 'g': Y = amp V R, 'd': Y = amp (-U) R."""
    V0, U0 = rect[0], rect[2]
    st = max(1, int(round(0.05 / h)))
    err, scale = 0.0, 0.0
    kc = _chart()
    for j in range(0, U.shape[1], st):
        for i in range(0, U.shape[2], st):
            V, Uu = V0 + j * h, U0 + i * h
            x, d, r, rV, rU, pV, pU, V0_, V1_, V2_, E, FV, FU, rhoV, rhoU = hs.jetsK(V, Uu, kc.par, kc.tab, 0.0)
            Yn = U[0, j, i] - x * r / SQ2 * U[2, j, i]
            Ye = amp * (V if kind == "g" else -Uu) * r
            err = max(err, abs(Yn - Ye))
            scale = max(scale, abs(Ye))
    return err / scale


_PL = {}


def planted(kind, h):
    """kind 'g', 'd': evolved exact mode, readout and the error of Y."""
    key = (kind, h, MUT if MUT in ("gauge_dF", "gauge_dphi") else None)
    if key not in _PL:
        rect = RECT_G if kind == "g" else RECT_D
        amp = BETA / 2.0
        kc = _chart()
        UL, UR = (growing_ray_data if kind == "g" else decaying_ray_data)(rect, h, BETA)
        U = hs.march_K(UL, UR, rect[0], rect[2], kc.par, kc.tab, h, h)
        ro = readout(U, rect[0], rect[2], h)
        ro["Yerr"] = np.array(y_num_err(U, rect, h, amp, kind))
        _PL[key] = ro
        del U
    return _PL[key]


# ----------------------------------------------------------------------------------------------------------------------
# derive: S1 the obstruction at S (sympy), S2 the master equation for both modes (sympy), S3 time reversal of the slice code
# ----------------------------------------------------------------------------------------------------------------------
def derive_S1():
    """Frobenius structure of the slice constraint M at S. With F_Pi, a = g/s^2 from the series of rho(Pi): the coefficient of ln|s| in d_N dphi of the exact mode."""
    import sympy as sp
    rho, Pi, G0, s = sp.symbols("rho Pi G0 s", positive=True)
    g1, g2 = sp.symbols("g1 g2")
    kap = 2 / sp.pi
    A = 1 - (rho ** 2 + 1) + kap * (rho + (rho ** 2 + 1) * sp.atan(rho))
    A1, A2 = sp.diff(A, rho).subs(rho, 0), sp.diff(A, rho, 2).subs(rho, 0)
    ok_a = sp.simplify(A1 - 2 * kap) == 0 and sp.simplify(A2 + 2) == 0
    rs = G0 * Pi + g1 * Pi ** 2 / 2 + g2 * Pi ** 3 / 6                            # rho(Pi), d rho/d Pi = A/(2 kappa Pi)
    Aser = sp.series(A.subs(rho, rs), Pi, 0, 4).removeO()
    eq = sp.expand(sp.diff(rs, Pi) * 2 * kap * Pi - Aser)
    sol = sp.solve([eq.coeff(Pi, 2), eq.coeff(Pi, 3)], [g1, g2], dict=True)[0]
    ok_g = sp.simplify(sol[g1] - A2 * G0 ** 2 / (2 * kap)) == 0
    G = sp.diff(rs, Pi)
    FPi0 = sp.simplify((-sp.diff(sp.log(G / kap), Pi)).subs(sol).subs(Pi, 0))      # F = -ln(G/kappa)
    if MUT == "s_zeta":
        FPi0 = 2 * FPi0
    ok_F = sp.simplify(FPi0 - G0 / kap) == 0
    a = sp.simplify(sp.limit((rs * sp.sqrt(1 + rs ** 2) / sp.sqrt(2)).subs(Pi, s ** 2) / s ** 2, s, 0))
    ok_a2 = sp.simplify(a - G0 / sp.sqrt(2)) == 0
    zeta = sp.simplify(-FPi0 / a)                                                # d z/d ln s per unit c2
    ok_z = sp.simplify(zeta + sp.pi / sp.sqrt(2)) == 0
    # exact identity of the M equation: g' + eps r s phi_Pi = 2 sqrt2 s G rho^2 / R (eps = -1; g = rho R/sqrt2, rho' = 2 s G, phi_Pi = sqrt2 G/(1 + rho^2))
    r_, G_ = sp.symbols("r_ G_", positive=True)
    R = sp.sqrt(1 + r_ ** 2)
    gprime = sp.diff(r_ * R / sp.sqrt(2), r_) * 2 * s * G_
    ident = sp.simplify(gprime - R * s * sp.sqrt(2) * G_ / (1 + r_ ** 2) - 2 * sp.sqrt(2) * s * G_ * r_ ** 2 / R)
    ok_i = ident == 0
    ok = bool(ok_a and ok_g and ok_F and ok_a2 and ok_z and ok_i)
    check("S1", "series at S: A'(0) = 2 kappa, A''(0) = -2, F_Pi(0) = G0/kappa, a = g/s^2 = G0/sqrt2, g' + eps r s phi_Pi = 2 sqrt2 s G rho^2/R exactly; so M forces "
          "g z' = -(Y_N' + s F_Pi Y_N) + O(s^3) and the exact mode (Y_N = c2 R) has z = -(pi/sqrt2) c2 ln|s| + O(1)", ok, f"F_Pi(0) = {FPi0}, a = {a}, zeta = {zeta}, identity residual {ident}")


def derive_S2():
    """Y = V R(rho(Pi)) and Y = U R(rho(Pi)) solve the Kruskal master equation 2 kappa^2 Pi u_VU + (1/2) Vhat u = 0, u = Y/rho: both reduce to the same ODE for f = R/rho."""
    import sympy as sp
    rho = sp.symbols("rho", positive=True)
    kap = 2 / sp.pi * (sp.Rational(101, 100) if MUT == "s_kappa" else 1)
    K = 2 / sp.pi
    A = 1 - (rho ** 2 + 1) + K * (rho + (rho ** 2 + 1) * sp.atan(rho))
    R = sp.sqrt(1 + rho ** 2)
    Ap = sp.diff(A, rho)
    Vh = 2 * A ** 2 / rho ** 2 - (1 + 1 / R ** 2) * A * Ap / rho + A ** 2 / R ** 4 - 2 * A / R ** 2
    f = R / rho
    fr = sp.diff(f, rho)
    # u = V f(Pi): u_VU = -V (Pi f)''; u = U f: u_VU = -U (Pi f)''; with d/dPi = (A/(2 kappa Pi)) d/drho the condition is kappa A f' + (A/2)(A f')' - Vhat f/2 = 0
    res = kap * A * fr + A * sp.diff(A * fr, rho) / 2 - Vh * f / 2
    vals = []
    import mpmath as mp
    mp.mp.dps = 50
    fn = sp.lambdify(rho, res, "mpmath")
    fs_ = sp.lambdify(rho, A * fr / 2, "mpmath")
    for x in ("0.3", "0.9", "2.5", "7.0", "-0.4", "-1.2"):
        vals.append(abs(fn(mp.mpf(x))) / max(abs(fs_(mp.mpf(x))), mp.mpf("1e-30")))
    simp = sp.simplify(res.subs(kap, 2 / sp.pi)) if MUT is None else None
    ok = max(vals) < 1e-40 or (simp is not None and simp == 0)
    check("S2", "Y = V R(rho(Pi)) and Y = U R(rho(Pi)) both solve the master equation 2 kappa^2 Pi u_VU + Vhat u/2 = 0 (u = Y/rho): the same ODE kappa A f' + A (A f')'/2 = Vhat f/2 for f = R/rho "
          "(its Pi-dependence drops out); evaluated at 50 digits at 6 points", ok, f"max relative residual {float(max(vals)):.1e}; symbolic simplify: {simp}")


def _rand_compact(rng):
    def parts(k):
        return fs.Data([(rng.uniform(-1, 1) * 1e-3, rng.uniform(-0.6, 0.6), rng.uniform(0.2, 0.45)) for _ in range(k)])
    return parts(2), parts(2), fs.Data([(1.0, 0.0, 0.6)]), fs.Data([(1.0, 0.3, 0.3), (-1.0, -0.3, 0.3)])


def derive_S3():
    """The slice code under time reversal T: (V, U) -> (-U, -V): (dr, dphi, dr', dphi') kept, (n, dphi_N) -> -(n, dphi_N). Then c1(TD) = c1(D), c2(TD) = -c2'(D), c2'(TD) = -c2(D)."""
    sl = fs.Slice(fs.bgc(C))
    if MUT == "s_swapfn":
        fs.MUT = "fn_sign"
    rng = np.random.default_rng(302)
    worst, lines = 0.0, []
    sets = [fs.data_G()] + [_rand_compact(rng) for _ in range(3)]
    for k, (sp_, spN, cp, cN) in enumerate(sets):
        d = fs.Datum(sl, sp_, spN, cp, cN, name=f"S3-{k}")
        flip = 1.0 if MUT == "s_noflip" else -1.0
        t = fs.Datum(sl, sp_, spN.scaled(flip), cp, cN.scaled(flip), name=f"S3T-{k}")
        sc_ = max(abs(d.c2), abs(d.c2p), 1e-300)
        dev = max(abs(t.c2 + d.c2p), abs(t.c2p + d.c2), abs(t.c1 - d.c1)) / sc_
        worst = max(worst, dev)
        lines.append(f"{dev:.1e} (c2 {d.c2:.4e}, c2' {d.c2p:.4e})")
    fs.MUT = None
    check("S3", "time reversal of the slice code (n and dphi_N flipped): c2(TD) = -c2'(D), c2'(TD) = -c2(D), c1(TD) = c1(D) on #258's datum G and 3 random compact data, relative to max(|c2|, |c2'|)",
          worst < 1e-10, "; ".join(lines))


def run_derive():
    t0 = time.time()
    for cid, f in (("S1", derive_S1), ("S2", derive_S2), ("S3", derive_S3)):
        if want(cid):
            f()
    bad = [c for c, ok in RESULTS if not ok]
    print(f"derive: {len(RESULTS) - len(bad)} of {len(RESULTS)} pass ({time.time() - t0:.0f} s)")
    return 1 if bad else 0


# ----------------------------------------------------------------------------------------------------------------------
# numerical checks
# ----------------------------------------------------------------------------------------------------------------------
def check_P1():
    """#258's exact growing mode data on [0.3, 1.6] (M = 0): it carries no mass (Theta_gen = 0); d_N dphi of the mode has the log slope of S1."""
    bg = fs.bgc(C)
    sl = fs.Slice(bg)
    if MUT == "alpha_off":
        fs.MUT = "mode_M"
    md = fs.ModeData(sl, c2=1e-3)
    fs.MUT = None
    worst = 0.0
    for s in (0.4, 0.8, 1.2, 1.5):
        U, UT, UN, UTT = md.fn(s)
        th = sc.theta_gen(sl, s, U[0], UT[0], U[2], 0.0)
        worst = max(worst, abs(th) / max(sc.theta_scale(sl, s, U[0], UT[0], U[2], 0.0), 1e-300))
    mscale = md.Mmax / max(md.Mscale, 1e-300)
    # d z / d ln s of the M equation for the exact Y_N = R (c2 = 1), integrated inward from s = 0.3
    K = fs.KAPPA

    def rhs(s, y):
        b = bg.at(s)
        rho, R = b["rho"], b["r"]
        rs_ = 2 * s * K * b["E"]
        gr = (1 + 2 * rho * rho) / (SQ2 * R)
        g, gs, Rs = rho * R / SQ2, gr * rs_, rho / R * rs_
        Fp = b["Fp"] * (2.0 if MUT == "zeta_F" else 1.0)
        return [(-(Rs + s * Fp * R) - (gs + fs.EPS * b["r"] * s * b["php"] + s * Fp * g) * y[0]) / g]
    sol = solve_ivp(rhs, [0.3, 1e-4], [0.0], dense_output=True, method="DOP853", rtol=1e-12, atol=1e-16)
    zeta = -math.pi / SQ2
    e = 1e-3
    slopes = []
    for s in (1e-2, 1e-3, 1e-4):
        slopes.append((sol.sol(s * (1 + e))[0] - sol.sol(s * (1 - e))[0]) / (math.log(1 + e) - math.log(1 - e)) / zeta)
    ok = worst <= 1e-8 and mscale <= 1e-9 and abs(slopes[-1] - 1.0) <= 1e-5 and abs(slopes[-2] - 1.0) <= 1e-3
    check("P1", "exact growing mode data of #258 (M = 0, [0.3, 1.6]): |Theta_gen| over its largest term <= 1e-8 (the mode carries no mass), M <= 1e-9 of its scale before the fix; "
          "dz/d ln s of the M equation with Y_N = c2 R over -(pi/sqrt2) c2 equals 1 to 1e-5 at s = 1e-4 and 1e-3 at s = 1e-3", ok,
          f"Theta {worst:.1e}; M {mscale:.1e}; slope ratio {slopes[0]:.7f}, {slopes[1]:.7f}, {slopes[2]:.7f}")


def _sc_ratio(ro, h):
    """Y_N''(0) / (-F_Pi(0) Y_N(0)) from the band of an evolved solution (Y = dr - g dphi, N = d_V + d_U, derivatives along the slice by 5-point stencils)."""
    band, sv = ro["band"], ro["band_s"]
    kc = _chart()
    k0 = int(np.argmin(np.abs(sv)))
    YN = []
    for k in range(k0 - 2, k0 + 3):
        Y = []
        for d in range(5):
            V = sv[k] + (d - 2) * h
            U = -sv[k] + (d - 2) * h
            x, dd, r, rV, rU, pV, pU, V0_, V1_, V2_, E, FV, FU, rhoV, rhoU = hs.jetsK(V, U, kc.par, kc.tab, 0.0)
            Y.append(band[0, k, d] - x * r / SQ2 * band[2, k, d])
        YN.append((-Y[4] + 8 * Y[3] - 8 * Y[1] + Y[0]) / (12 * h))
    YN0 = YN[2]
    YN2 = (-YN[4] + 16 * YN[3] - 30 * YN[2] + 16 * YN[1] - YN[0]) / (12 * h * h)
    FPi = fs.bgc(C).at(1e-6)["Fp"] * (0.0 if MUT == "sc_off" else 1.0)
    return YN2 / (-FPi * YN0), YN0


def check_L1():
    """The smoothness condition of S1 on evolved smooth solutions: Y_N''(0) = -F_Pi(0) Y_N(0) (any C^2 representation); the exact mode has Y_N''(0) = 0."""
    out, ok = [], True
    for name in ("A0", "A1R", "A1L", "L", "G", "Lm", "Rb", "Gb"):
        pr = [_sc_ratio(run_spec(CLASSES[name], h), h) for h in HF]
        vals = [x[0] for x in pr]
        r = rich(vals[-2], vals[-1])
        yn0 = rich(pr[-2][1], pr[-1][1])
        scored = abs(yn0) >= 1e-4 * UNIT * 1e3                         # |Y_N(0)| >= 1e-4: the ratio is conditioned only when Y_N(0) is not at the truncation seed
        good = (abs(r - 1.0) <= 1e-3) or not scored
        ok = ok and good
        out.append(f"{name}: {vals[0]:.6f}, {vals[1]:.6f} -> {r:.6f} (Y_N(0) {yn0:.2e}{'' if scored else ', not scored'})")
    # the exact mode: Y_N = c2 R(rho(Pi)): second derivative along the slice at S from the closed form
    bg = fs.bgc(C)
    e = 2e-3
    RR = [bg.at(k * e)["r"] for k in (-2, -1, 0, 1, 2)]
    d2 = (-RR[4] + 16 * RR[3] - 30 * RR[2] + 16 * RR[1] - RR[0]) / (12 * e * e)
    ok = ok and abs(d2) <= 1e-6
    check("L1", "smoothness condition on gauge-invariant Y: Y_N''(0)/(-F_Pi(0) Y_N(0)) = 1 to 1e-3 (Richardson, h = 0.00125, 0.000625) for the evolved smooth classes A0, A1R, A1L, L, G, Lm, Rb, Gb with |Y_N(0)| >= 1e-4; "
          "the exact mode Y_N = c2 R has Y_N''(0) = 0 (|.| <= 1e-6 c2)", ok, "; ".join(out) + f"; exact mode R'' (s = 0) = {d2:.1e}, against -F_Pi = {-fs.bgc(C).at(1e-6)['Fp']:.4f}")


def check_L2():
    """The exact growing mode in the double-null gauge of #269's family: dphi has a V ln V term at V -> 0 whose coefficient is that of S1; changing the base ray (the residual gauge xi^u -> xi^u + eta(U)) leaves it."""
    zeta = -math.pi / SQ2 * 0.5                           # beta = 1: c2 = 1/2
    out, ok = [], True
    for U in (-0.5, -1.0):
        for Vb in (0.05, 1.0):
            sl = []
            for V in (1e-3, 1e-4, 1e-5, 1e-6):
                e = 5e-3
                Vs = V * (1 + e * np.array([-4.0, -2.0, 0.0, 2.0, 4.0]))
                outs = []
                for fac in (math.exp(-0.1), math.exp(0.1)):
                    Vq = Vs * fac
                    xi = xi_line(Vq, U, Vb)
                    ph = np.array([g_fields(Vq[k], U, xi[0, k], xi[1, k])[2] for k in range(5)])
                    outs.append((-ph[4] + 8 * ph[3] - 8 * ph[1] + ph[0]) / (6.0 * (Vq[2] - Vq[0])))     # 4th-order d_V: nodes at +-1, +-2 steps of (Vq[2] - Vq[0])/2
                sl.append((outs[1] - outs[0]) / 0.2 / zeta)
            good = abs(sl[-1] - 1.0) <= 1e-3 and abs(sl[-2] - 1.0) <= 1e-2
            ok = ok and good
            out.append(f"U {U}, Vb {Vb}: {sl[0]:.5f}, {sl[1]:.5f}, {sl[2]:.5f}, {sl[3]:.5f}")
    check("L2", "exact double-null fields (#269's EF family + gauge vector): d_V dphi = zeta ln V + const with zeta = -(pi/sqrt2) c2 (S1's coefficient) to 1e-3 at V = 1e-6 and 1e-2 at V = 1e-5, "
          "for two U and two base rays (the residual gauge only shifts the constant)", ok, "; ".join(out))


def _rel(a, b, floor):
    return abs(a - b) / max(abs(b), floor)


def check_P2():
    """the static mass mode through the Kruskal solver (#259 T8): (c1, c2, c2') = (amp/2, 0, 0) on both horizons and by the slice functionals."""
    amp = 2e-3
    rds = []
    for h in HL[1:]:
        UL, UR = sc.mode_arrays(h, amp)
        ro = run_arrays(UL, UR, RECT, h, tag=f"static{amp}")
        rds.append(readings(ro, h))
    c1 = amp / 2.0
    R = lambda key, idx=None: rich(*(((r[key][idx] if idx is not None else r[key]) for r in rds)))
    vals = dict(c1S=R("c1S"), c2S=R("c2S"), c2pS=R("c2pS"), c1p=R("fitp", 0), c1m=R("fitm", 0), c2p=R("fitp", 1), c2m=R("fitm", 1), c2pH1=R("c2p_Hm1"))
    ok = (abs(vals["c1S"] - c1) <= 1e-6 * c1 and abs(vals["c1p"] - c1) <= 1e-6 * c1 and abs(vals["c1m"] - c1) <= 1e-6 * c1 and max(abs(vals[k]) for k in ("c2S", "c2pS", "c2p", "c2m", "c2pH1")) <= 1e-6 * c1)
    check("P2", "static mass mode (amplitude dL = 2e-3, exact ray data of #259) through the solver, Richardson of h = 0.0025, 0.00125: c1 = 1e-3 at S and on both horizons to 1e-6 relative; c2, c2' by the slice "
          "functionals and both fits <= 1e-6 c1", ok, ", ".join(f"{k} = {v:.4e}" for k, v in vals.items()))


def check_P3():
    """the growing mode (exact double-null ray data, V > 0.05) through the solver: Y error, c2 by the horizon fit and by the slice functional at V0 = 0.5, 1."""
    c2 = BETA / 2.0
    rds = [planted("g", h) for h in HL]
    errs = [float(r["Yerr"]) for r in rds]
    rt = ratios(errs)
    rr = [readings(r, h) for r, h in zip(rds, HL)]
    fit = rich(rr[-2]["fitp"][1], rr[-1]["fitp"][1])
    sls = [rich(rr[-2][k], rr[-1][k]) for k in ("c2_Hp05", "c2_Hp1")]
    dev = max([abs(fit - c2) / c2] + [abs(x - c2) / c2 for x in sls])
    c1fit = rich(rr[-2]["fitp"][0], rr[-1]["fitp"][0])
    ok = errs[-1] <= 1e-3 and abs(rt[-1] - 4.0) <= 0.3 and dev <= 1e-6 and abs(c1fit) <= 1e-6 * c2
    check("P3", "growing mode Y = c2 V R from exact ray data on V >= 0.05: error of Y over its maximum <= 1e-3 at h = 0.00125, last ratio 4.0 +- 0.3; c2 = beta/2 by the horizon fit and by the slice functional "
          "at V0 = 0.5, 1 to 1e-6 relative (Richardson h = 0.0025, 0.00125); c1 of the fit <= 1e-6 c2", ok,
          f"Y errors {', '.join(f'{e:.2e}' for e in errs)}; ratio {rt[-1]:.2f}; fit {fit:.9e}, slices {sls[0]:.9e}, {sls[1]:.9e} (exact {c2:.1e}); c1 {c1fit:.1e}")


def check_P4():
    """the decaying mode (time reverse of the growing mode's data) through the solver: c2' by the fit on V = 0 and by the slice functional at U0 = -0.5, -1."""
    c2p = -BETA / 2.0
    rds = [planted("d", h) for h in HL]
    errs = [float(r["Yerr"]) for r in rds]
    rt = ratios(errs)
    rr = [readings(r, h) for r, h in zip(rds, HL)]
    fit = rich(rr[-2]["fitm"][1], rr[-1]["fitm"][1])
    sls = [rich(rr[-2][k], rr[-1][k]) for k in ("c2p_Hm05", "c2p_Hm1")]
    dev = max([abs(fit - c2p) / abs(c2p)] + [abs(x - c2p) / abs(c2p) for x in sls])
    c1fit = rich(rr[-2]["fitm"][0], rr[-1]["fitm"][0])
    ok = errs[-1] <= 1e-3 and abs(rt[-1] - 4.0) <= 0.3 and dev <= 1e-6 and abs(c1fit) <= 1e-6 * abs(c2p)
    check("P4", "decaying mode Y = c2' U R from the time reversed exact ray data (U <= -0.05): error of Y over its maximum <= 1e-3 at h = 0.00125, last ratio 4.0 +- 0.3; c2' = -beta/2 by the fit on V = 0 "
          "and by the slice functional d_U at U0 = -0.5, -1 to 1e-6 relative; c1 of the fit <= 1e-6 |c2'|", ok,
          f"Y errors {', '.join(f'{e:.2e}' for e in errs)}; ratio {rt[-1]:.2f}; fit {fit:.9e}, slices {sls[0]:.9e}, {sls[1]:.9e} (exact {c2p:.1e}); c1 {c1fit:.1e}")


def _rand_spec(rng):
    return dict(a0=rng.uniform(-1e-3, 1e-3), a1R=rng.uniform(-1e-3, 1e-3), a1L=rng.uniform(-1e-3, 1e-3),
                L=(rng.uniform(-1e-3, 1e-3), rng.uniform(-1.2, 0.4), rng.uniform(0.2, 0.5)), R=(rng.uniform(-1e-3, 1e-3), rng.uniform(-0.08, 0.0), rng.uniform(0.02, 0.05)))


def _horizon_runs(spec, hl):
    return [run_spec(spec, h) for h in hl]


def check_H1():
    out, ok = [], True
    for name in ("L", "G", "A0", "A1R", "A1L", "Lm", "Rb", "Gb"):
        ros = _horizon_runs(CLASSES[name], HL)
        rr = [readings(r, h) for r, h in zip(ros, HL)]
        res = [x["fitm"][2] for x in rr]
        rt = ratios(res)
        scale = float(np.max(np.abs(ros[-1]["Hm"])))
        big = scale >= 1e-2 * UNIT                                         # a class with max|dr| below 1% of the amplitude on V = 0 has only the truncation seed there (as #248's class R)
        good = res[-1] * scale <= 2e-5 * max(scale, UNIT) and (all(3.5 <= q <= 4.5 for q in rt) or not big)
        ok = ok and good
        out.append(f"{name}: res {res[0]:.1e}, {res[1]:.1e}, {res[2]:.1e} (ratios {rt[0]:.2f}, {rt[1]:.2f}), max|dr| on V = 0 {scale:.2e}{'' if big else ' (seed only, ratios not scored)'}")
    check("H1", "horizon law on the past horizon V = 0: dr(0, U) is linear in U on [-2, 0.06]; least-squares residual <= 2e-5 of max(max|dr|, 1e-3) at h = 0.00125, and falling by 4 +- 0.5 per halving where max|dr| >= 1e-5; classes L, G, A0, A1R, A1L, Lm, Rb, Gb",
          ok, "; ".join(out))


def check_H2():
    out, ok = [], True
    rng = np.random.default_rng(3021)
    items = [(n, CLASSES[n]) for n in ("L", "G", "A0", "A1R", "A1L", "Lm", "Rb", "Gb")] + [(f"rand{k}", _rand_spec(rng)) for k in range(8)]
    c2ps = []
    worst = 0.0
    for name, spec in items:
        ros = _horizon_runs(spec, HF)
        rr = [readings(r, h) for r, h in zip(ros, HF)]
        fit = rich(rr[0]["fitm"][1], rr[1]["fitm"][1])
        scale = max(abs(fit), 1e-2 * UNIT)
        devs = [abs(rich(rr[0][k], rr[1][k]) - fit) / scale for k in ("c2pS", "c2p_Hm05", "c2p_Hm1")]
        worst = max(worst, max(devs))
        c2ps.append(abs(fit))
        out.append(f"{name}: c2' {fit:.6e}, slices U0 = 0, -0.5, -1: {devs[0]:.1e}, {devs[1]:.1e}, {devs[2]:.1e}")
    ok = worst <= 1e-6 and min(c2ps[8:]) >= 1e-6
    check("H2", "c2' = d_U dr(S): the slice functional (d_U from the slice values and the normal derivative) at the slices crossing V = 0 at U0 = 0 (through S), -0.5, -1 equals the fit on V = 0 to 1e-6 of max(|c2'|, 1e-5) "
          "(Richardson, h = 0.00125, 0.000625) for 8 classes and 8 random data; no random datum has |c2'| < 1e-6", ok, f"worst {worst:.1e}; |c2'| of the random data from {min(c2ps[8:]):.2e} to {max(c2ps[8:]):.2e}; " + "; ".join(out[:8]))


def check_H3():
    """the three modes reproduce Y_h on both horizons: Y_h(V) = c1 + c2 V on U = 0 and Y_h(U) = c1 + c2' U on V = 0, with (c1, c2, c2') read at S by the slice."""
    rng = np.random.default_rng(3021)
    items = [(n, CLASSES[n]) for n in ("L", "G", "A0", "A1R", "A1L", "Lm", "Rb", "Gb")] + [(f"rand{k}", _rand_spec(rng)) for k in range(8)]
    worst, out = 0.0, []
    for name, spec in items:
        ros = _horizon_runs(spec, HF)
        rr = [readings(r, h) for r, h in zip(ros, HF)]
        c1 = rich(rr[0]["c1S"], rr[1]["c1S"])
        c2 = rich(rr[0]["c2S"], rr[1]["c2S"])
        c2p = rich(rr[0]["c2pS"], rr[1]["c2pS"])
        Hp = rich_arr(ros[0]["Hp"], ros[1]["Hp"])
        Hm = rich_arr(ros[0]["Hm_off" if MUT == "law_off" else "Hm"], ros[1]["Hm_off" if MUT == "law_off" else "Hm"])
        V, Ua = ros[0]["V"], ros[0]["Ua"]
        devp = np.max(np.abs(Hp - (c1 + c2 * V)))
        c1m = 0.0 if MUT == "law_c1" else c1
        devm = np.max(np.abs(Hm - (c1m + c2p * Ua)))
        scale = max(float(np.max(np.abs(Hp))), float(np.max(np.abs(Hm))), UNIT)
        d = max(devp, devm) / scale
        worst = max(worst, d)
        out.append(f"{name}: c1 {c1:.3e}, c2 {c2:.3e}, c2' {c2p:.3e}; dev {devp / scale:.1e} (U = 0), {devm / scale:.1e} (V = 0)")
    check("H3", "the three modes reproduce the horizon values of evolved data: dr(V, 0) - (c1 + c2 V) and dr(0, U) - (c1 + c2' U), with c1, c2, c2' from the slice at S (Richardson, h = 0.00125, 0.000625), "
          "over max(max|dr_h|, 1e-3) <= 1e-6 at every node of both horizons, for 8 classes and 8 random data", worst <= 1e-6, f"worst {worst:.1e}; " + "; ".join(out[:5]) + " ...")


def rich_arr(a, b):
    """Richardson of two arrays on grids h and h/2 (the coarse nodes are every second fine node)."""
    return (4.0 * b[::2] - a) / 3.0


CHECKS = [("P1", check_P1), ("L1", check_L1), ("L2", check_L2), ("P2", check_P2), ("P3", check_P3), ("P4", check_P4), ("H1", check_H1), ("H2", check_H2), ("H3", check_H3)]
SCOPE = ("S1", "S2", "S3", "P1", "L1", "L2", "P2", "P3", "P4", "H1", "H2", "H3")
DECLARED = {
    "s_zeta": {"S1"}, "s_kappa": {"S2"}, "s_noflip": {"S3"}, "s_swapfn": {"S3"},
    "fn_swap": {"P2", "P4", "H2", "H3"}, "fn_aT": {"P2", "P4", "H2", "H3"},
    "law_c1": {"H3"}, "law_x": {"P4", "H2"}, "law_off": {"P2", "P4", "H1", "H2", "H3"},
    "gauge_dF": {"P3", "P4"}, "gauge_dphi": {"L2", "P3", "P4"},
    "alpha_off": {"P1"}, "sc_off": {"L1"}, "zeta_F": {"P1"},
}


def run_checks():
    t0 = time.time()
    for cid, f in CHECKS:
        if want(cid):
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
    print("== the three exact modes through the Kruskal solver (Richardson of h = 0.0025, 0.00125); exact values in the last column ==")
    amp = 2e-3
    rds = []
    for h in HL[1:]:
        UL, UR = sc.mode_arrays(h, amp)
        rds.append(readings(run_arrays(UL, UR, RECT, h, tag=f"static{amp}"), h))
    R = lambda key, idx=None: rich(*(((r[key][idx] if idx is not None else r[key]) for r in rds)))
    print(f"static (dL = 2e-3):  c1 = {R('c1S'):.9e} (S), {R('fitp', 0):.9e} (H+ fit), {R('fitm', 0):.9e} (H- fit); c2 = {R('c2S'):.2e}, {R('fitp', 1):.2e}; c2' = {R('c2pS'):.2e}, {R('fitm', 1):.2e};  exact (1e-3, 0, 0)")
    for kind, keys, exact in (("g", ("c2_Hp05", "c2_Hp1", "c2_Hp2"), BETA / 2), ("d", ("c2p_Hm05", "c2p_Hm1", "c2p_Hm2"), -BETA / 2)):
        rr = [readings(planted(kind, h), h) for h in HL]
        fi = "fitp" if kind == "g" else "fitm"
        errs = [float(planted(kind, h)["Yerr"]) for h in HL]
        sl = [rich(rr[-2][k], rr[-1][k]) for k in keys]
        nm = "growing  " if kind == "g" else "decaying "
        print(f"{nm}(beta = 2e-3): Y error {errs[0]:.2e}, {errs[1]:.2e}, {errs[2]:.2e}; {'c2' if kind == 'g' else 'c2' + chr(39)}: fit {rich(rr[-2][fi][1], rr[-1][fi][1]):.9e}, slices at 3 crossings " + ", ".join(f"{x:.9e}" for x in sl) + f"; c1 fit {rich(rr[-2][fi][0], rr[-1][fi][0]):.1e};  exact {exact:.3e}")
    print("== classes and random data: c1, c2, c2' by the slice at S, by the fits (Richardson of h = 0.00125, 0.000625) ==")
    print("| datum | c1 (S) | c1 (fit U = 0) | c1 (fit V = 0) | c2 (S) | c2 (fit) | c2' (S) | c2' (fit) |")
    print("|---|---|---|---|---|---|---|---|")
    rng = np.random.default_rng(3021)
    items = [(n, CLASSES[n]) for n in ("L", "G", "A0", "A1R", "A1L", "Lm", "Rb", "Gb")] + [(f"rand{k}", _rand_spec(rng)) for k in range(8)]
    for name, spec in items:
        rr = [readings(run_spec(spec, h), h) for h in HF]
        g = lambda key, idx=None: rich(*(((r[key][idx] if idx is not None else r[key]) for r in rr)))
        print(f"| {name} | {g('c1S'):.6e} | {g('fitp', 0):.6e} | {g('fitm', 0):.6e} | {g('c2S'):.6e} | {g('fitp', 1):.6e} | {g('c2pS'):.6e} | {g('fitm', 1):.6e} |")
    print("== the smoothness condition Y_N''(0)/(-F_Pi(0) Y_N(0)), Richardson of h = 0.00125, 0.000625 ==")
    for name in ("A0", "A1R", "A1L", "L", "G", "Lm", "Rb", "Gb"):
        pr = [_sc_ratio(run_spec(CLASSES[name], h), h) for h in HF]
        print(f"{name}: {pr[0][0]:.6f}, {pr[1][0]:.6f} -> {rich(pr[0][0], pr[1][0]):.6f}; Y_N(0) = {rich(pr[0][1], pr[1][1]):.4e}")
    print(f"F_Pi(0) = {fs.bgc(C).at(1e-6)['Fp']:.8f} (pi G0/2 with G0 = {2.0 * fs.bgc(C).at(1e-6)['Fp'] / math.pi:.8f}), -pi/sqrt2 = {-math.pi / SQ2:.8f}")
    return 0


if __name__ == "__main__":
    cmd = next((a for a in sys.argv[1:] if not a.startswith("--")), "check")
    if cmd == "derive":
        sys.exit(run_derive())
    if cmd == "check":
        sys.exit(run_checks())
    if cmd == "table":
        sys.exit(report_table())
    if cmd == "mutations":
        i = sys.argv.index("mutations")
        sys.exit(report_mutations(sys.argv[i + 1], sys.argv[i + 2]))
    sys.exit(f"unknown command {cmd}")

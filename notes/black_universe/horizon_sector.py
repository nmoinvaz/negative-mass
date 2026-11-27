#!/usr/bin/env python3
"""The horizon sector of the c = -1 black universe in Cauchy data across the horizon (#248).

#232 found that at c = -1 the linear double-null response grows at kappa for data with Y_h = c1 + c2 e^{b p} on the horizon, c2 != 0. This module asks whether c2
is free in regular data on a spacelike slice through the horizon, in the Kruskal chart of the member, which contains the bifurcation surface S and the interior (rho < 0).

    horizon_sector.py derive [--mutate=NAME]        sympy / mpmath checks of the formulation (D-series)
    horizon_sector.py check  [--only=..] [--mutate=NAME]   numerical checks, planted case and control (N-series); exit 1 when one fails
    horizon_sector.py table                         the c2 table of the classes through S and across the horizon
    horizon_sector.py mutations DIR PREFIX          fired sets of the stored logs DIR/PREFIX<name>.log against the declared ones

Chart. Kruskal null coordinates x+ = V = e^{b p}, x- = U = -Q(sigma) (b = sqrt2 kappa), Pi = -V U = e^{2 kappa x} (static region V > 0 > U), the metric
-2 e^{-F} dV dU + r^2 dOmega^2 with e^{-F} = G(rho)/kappa, G = d rho/d Pi = A/(2 kappa Pi), r = R(rho), phi = sqrt2 arctan rho, rho = rho(Pi) analytic through
Pi = 0 (both horizons and S) into Pi < 0 (the interior). The linear system U_VU = B0 U + Bp U_V + Bm U_U, U = (dr, dF, dphi), is that of linear.py (derived in
linear_derive.py) with the background jets of this chart; the cell scheme is linear.py's.
"""
import math
import os
import sys
import time

import numba as nb
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import linear as lin  # noqa: E402

SQ2 = math.sqrt(2.0)
MUT = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mutate=")), None)
ONLY = next((a.split("=", 1)[1].split(",") for a in sys.argv[1:] if a.startswith("--only=")), None)
# extra parameter slots after linear.py's NPAR = 29 (mutation multipliers of this module)
P_MEF, P_MFV, P_KIND2, P_MX = 29, 30, 31, 32          # e^{-F} scale, F_V sign, 1 = Schwarzschild (planted case), spare
NPAR2 = 33


# ----------------------------------------------------------------------------------------------------------------------
# the background in the Kruskal chart: the member (rho(Pi) through the horizon) or Schwarzschild (the planted case, r_s = 1)
# ----------------------------------------------------------------------------------------------------------------------
@nb.njit(cache=True)
def jetsK(V, U, par, tab, guess):
    """(x, d, r, r_V, r_U, phi_V, phi_U, V0, V1, V2, e^{-F}, F_V, F_U, rho_V, rho_U) at (V, U), Pi = -V U.
    Member: x = rho, delta = Pi E(rho), G = d rho/d Pi = E (A/delta)/(2 kappa), e^{-F} = G/kappa, F = ln kappa - ln G.
    Schwarzschild (r_s = 2M = 1, kappa = 1/2): (r - 1) e^r = Pi, G = e^{-r}/r, e^{-F} = G/kappa, phi = 0, V = 0."""
    Pi = -V * U
    if par[P_KIND2] == 1.0:
        r = 1.0 + Pi * math.exp(-1.0) if guess <= 0.0 else guess
        if r < 1e-3:
            r = 1e-3
        for _ in range(100):
            g = (r - 1.0) * math.exp(r) - Pi
            gp = r * math.exp(r)
            rn = r - g / gp
            if rn <= 0.0:
                rn = 0.5 * r
            if abs(rn - r) < 1e-15 * rn:
                r = rn
                break
            r = rn
        G = math.exp(-r) / r
        dl = -1.0 - 1.0 / r
        rV = -U * G
        rU = -V * G
        return r, r, r, rV, rU, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0 * G * par[P_MEF], -dl * rV * par[P_MFV], -dl * rU * par[P_MFV], rV, rU
    d = lin._solve_delta(Pi, par, tab, guess)
    rho = par[lin.P_RH] + d
    E = math.exp(lin._lnE(rho, par, tab))
    aod = lin._Aod(d, par)
    G = E * aod / (2.0 * par[lin.P_KAP])
    rV = -U * G
    rU = -V * G
    dl = lin._gfun(d, par) + lin._dAod(d, par) / aod            # d ln G / d rho
    R = math.sqrt(rho * rho + 1.0)
    Rp = rho / R
    Phip = SQ2 / (1.0 + rho * rho)
    V0, V1, V2 = lin.pot3(rho, par)
    return rho, d, R, Rp * rV, Rp * rU, Phip * rV, Phip * rU, V0, V1, V2, G / par[lin.P_KAP] * par[P_MEF], -dl * rV * par[P_MFV], -dl * rU * par[P_MFV], rV, rU


@nb.njit(cache=True)
def coefK(V, U, par, tab, guess, B0, Bp, Bm):
    """B0, Bp (= B_V), Bm (= B_U) of the linear system at (V, U); same entries as linear.coef, background jets of the Kruskal chart."""
    x, d, r, rp, rm, php, phm, V0, V1, V2, E, FV, FU, rho_V, rho_U = jetsK(V, U, par, tab, guess)
    mv, mvp, mvpp, meps, mcross, mde, mrm = par[lin.P_MV], par[lin.P_MVP], par[lin.P_MVPP], par[lin.P_MEPS], par[lin.P_MCROSS], par[lin.P_MDE], par[lin.P_MRM]
    V0 *= mv
    V1 *= mvp
    V2 *= mvpp
    for i in range(3):
        for j in range(3):
            B0[i, j] = 0.0
            Bp[i, j] = 0.0
            Bm[i, j] = 0.0
    r2 = r * r
    B0[0, 0] = E * V0 / 2 + E / (2 * r2) + rm * rp / r2
    B0[0, 1] = mde * (-E * V0 * r / 2 + E / (2 * r))
    B0[0, 2] = E * V1 * r / 2
    Bp[0, 0] = -mrm * rm / r
    Bm[0, 0] = -mrm * rp / r
    B0[1, 0] = 2 * E / (r2 * r) + 4 * rm * rp / (r2 * r)
    B0[1, 1] = mde * E / r2
    Bp[1, 0] = -2 * rm / r2
    Bm[1, 0] = -2 * rp / r2
    Bp[1, 2] = -meps * phm
    Bm[1, 2] = -meps * php
    B0[2, 0] = (phm * rp + php * rm) / r2
    B0[2, 1] = -mde * meps * E * V1 / 2
    B0[2, 2] = meps * E * V2 / 2
    Bp[2, 0] = -phm / r
    Bm[2, 0] = -php / r
    Bp[2, 2] = -mcross * rm / r
    Bm[2, 2] = -mcross * rp / r
    return d


@nb.njit(cache=True)
def _cellK(Sv, Ev, Wv, V, U, hV, hU, par, tab, guess, B0, Bp, Bm, M, rhs, out):
    """linear._cell with the Kruskal background: N = E + W - S + hV hU (B0 Ubar + Bp Ubar_V + Bm Ubar_U), cell-centre averages."""
    d = coefK(V + 0.5 * hV, U + 0.5 * hU, par, tab, guess, B0, Bp, Bm)
    H = hV * hU
    for i in range(3):
        s = Ev[i] + Wv[i] - Sv[i]
        for j in range(3):
            M[i, j] = (1.0 if i == j else 0.0) - H * B0[i, j] / 4 - 0.5 * hU * Bp[i, j] - 0.5 * hV * Bm[i, j]
            s += H * (B0[i, j] * (Sv[j] + Ev[j] + Wv[j]) / 4 + Bp[i, j] * (Ev[j] - Sv[j] - Wv[j]) / (2 * hV)
                      + Bm[i, j] * (Wv[j] - Sv[j] - Ev[j]) / (2 * hU))
        rhs[i] = s
    lin._solve3(M, rhs, out)
    return d


@nb.njit(cache=True)
def march_K(UL, UR, V0, U0, par, tab, hV, hU):
    """The whole rectangle V in [V0, V0 + Nv hV], U in [U0, U0 + Nu hU]. UL[:, i] = data on V = V0 (U_i = U0 + i hU), UR[:, j] = data on U = U0 (V_j = V0 + j hV).
    Returns U[:, j, i]."""
    Nu = UL.shape[1] - 1
    Nv = UR.shape[1] - 1
    Ua = np.zeros((3, Nv + 1, Nu + 1))
    for q in range(3):
        for i in range(Nu + 1):
            Ua[q, 0, i] = UL[q, i]
        for j in range(Nv + 1):
            Ua[q, j, 0] = UR[q, j]
    B0 = np.zeros((3, 3))
    Bp = np.zeros((3, 3))
    Bm = np.zeros((3, 3))
    M = np.zeros((3, 3))
    rhs = np.zeros(3)
    out = np.zeros(3)
    Sv = np.zeros(3)
    Ev = np.zeros(3)
    Wv = np.zeros(3)
    for j in range(Nv):
        guess = 0.0
        for i in range(Nu):
            for q in range(3):
                Sv[q] = Ua[q, j, i]
                Ev[q] = Ua[q, j + 1, i]
                Wv[q] = Ua[q, j, i + 1]
            guess = _cellK(Sv, Ev, Wv, V0 + j * hV, U0 + i * hU, hV, hU, par, tab, guess, B0, Bp, Bm, M, rhs, out)
            for q in range(3):
                Ua[q, j + 1, i + 1] = out[q]
    return Ua


# ----------------------------------------------------------------------------------------------------------------------
# the chart object
# ----------------------------------------------------------------------------------------------------------------------
class KChart:
    """Member c in the Kruskal chart. The lnE table of linear.Chart is reused (Pi = 1 at rho_c = 1, the corner of #232); it is extended to cover Pi in [-pi_min, pi_max]."""

    def __init__(self, c, pi_max=60.0, pi_min=2.0, mut=None):
        self.c = float(c)
        # linear.Chart sizes its table by (p_max, q_int): rho_hi = rho_c + 1.25 p_max/sqrt2 + 3, rho_lo = rho_h - max(4 q_int, 0.5); choose both generously
        p_max = 2.0 * math.log(max(pi_max, 2.0)) / (SQ2 * 1.0) + 8.0
        mut = dict(mut or {})
        self.ch = lin.Chart(c, rho_c=1.0, sig_star=16.0, p_max=p_max, q_int=max(pi_min, 1.0) * 2.0, mut={k: v for k, v in mut.items() if k < lin.NPAR})
        self.par = np.concatenate([self.ch.par, np.array([1.0, 1.0, 0.0, 0.0])])      # P_MEF, P_MFV, P_KIND2, P_MX
        self.tab = self.ch.tab
        for k, v in mut.items():
            self.par[k] = v
        self.mem = self.ch.mem
        self.kappa = self.mem.kappa
        self.b = SQ2 * self.kappa


# ----------------------------------------------------------------------------------------------------------------------
# exact pure-gauge solutions L_xi(member), xi = xi(V) d_V + eta(U) d_U
# ----------------------------------------------------------------------------------------------------------------------
@nb.njit(cache=True)
def gaugeK(par, tab, V, U, xi, xip, eta, etap):
    x, d, r, rV, rU, pV, pU, V0, V1, V2, E, FV, FU, rho_V, rho_U = jetsK(V, U, par, tab, 0.0)
    out = np.zeros(3)
    out[0] = xi * rV + eta * rU
    out[2] = xi * pV + eta * pU
    out[1] = xi * FV + eta * FU - xip - etap
    return out


def gauge_arrays_K(kc, V0, U0, hV, hU, Nv, Nu, xi_f, eta_f):
    U = np.zeros((3, Nv + 1, Nu + 1))
    for j in range(Nv + 1):
        xi, xip = xi_f(V0 + j * hV)
        for i in range(Nu + 1):
            eta, etap = eta_f(U0 + i * hU)
            U[:, j, i] = gaugeK(kc.par, kc.tab, V0 + j * hV, U0 + i * hU, xi, xip, eta, etap)
    return U


# ----------------------------------------------------------------------------------------------------------------------
# characteristic data on the two rays: dphi a Gaussian pulse, dF = 0, dr from the linearized constraint along the ray with corner constants
# ----------------------------------------------------------------------------------------------------------------------
@nb.njit(cache=True)
def _ray_rhsK(s, axis, fixed, par, tab, guess, y0, y1, amp, s0, w):
    """a'' = r phi' e' - F' a' + phi'^2 a/2 along V (axis 0, U = fixed) or along U (axis 1, V = fixed); e is a Gaussian pulse in the ray coordinate."""
    if axis == 0:
        x, d, r, rV, rU, pV, pU, V0, V1, V2, E, FV, FU, rhoV, rhoU = jetsK(s, fixed, par, tab, guess)
        phd, fd = pV, FV
    else:
        x, d, r, rV, rU, pV, pU, V0, V1, V2, E, FV, FU, rhoV, rhoU = jetsK(fixed, s, par, tab, guess)
        phd, fd = pU, FU
    xx = (s - s0) / w
    ep_ = amp * (-2.0 * xx / w) * math.exp(-xx * xx)
    return y1, r * phd * ep_ - par[lin.P_MCON] * fd * y1 + 0.5 * phd * phd * y0, d


@nb.njit(cache=True)
def ray_dataK(par, tab, axis, fixed, s_start, n, h, amp, s0, w, a0, a1, nsub):
    """(3, n+1): dr, dF = 0, dphi along a ray with nodes s = s_start + i h; RK4 with nsub sub-steps per node. dphi vanishes at the first node."""
    out = np.zeros((3, n + 1))
    y0, y1 = a0, a1
    guess = 0.0
    out[0, 0] = y0
    xx0 = (s_start - s0) / w
    g0 = amp * math.exp(-xx0 * xx0)
    hh = h / nsub
    for i in range(n):
        for k in range(nsub):
            s = s_start + i * h + k * hh
            k1a, k1b, guess = _ray_rhsK(s, axis, fixed, par, tab, guess, y0, y1, amp, s0, w)
            k2a, k2b, guess = _ray_rhsK(s + 0.5 * hh, axis, fixed, par, tab, guess, y0 + 0.5 * hh * k1a, y1 + 0.5 * hh * k1b, amp, s0, w)
            k3a, k3b, guess = _ray_rhsK(s + 0.5 * hh, axis, fixed, par, tab, guess, y0 + 0.5 * hh * k2a, y1 + 0.5 * hh * k2b, amp, s0, w)
            k4a, k4b, guess = _ray_rhsK(s + hh, axis, fixed, par, tab, guess, y0 + hh * k3a, y1 + hh * k3b, amp, s0, w)
            y0 += hh / 6 * (k1a + 2 * k2a + 2 * k3a + k4a)
            y1 += hh / 6 * (k1b + 2 * k2b + 2 * k3b + k4b)
        out[0, i + 1] = y0
        xx = (s_start + (i + 1) * h - s0) / w
        out[2, i + 1] = amp * math.exp(-xx * xx) - g0
    return out


class KSetup:
    """A characteristic rectangle V in [V0, V1], U in [U0, U1] of the Kruskal chart (V0 = 0: the left ray V = 0 passes through S when U0 < 0 < U1).
    spec = dict(R=(amp, s0, w) pulse on the ray U = U0 in V, L=(amp, s0, w) pulse on the ray V = V0 in U, a0 = dr at the corner, a1R = dr_V at the corner, a1L = dr_U at the corner)."""

    def __init__(self, c, h, V0, V1, U0, U1, spec, kc=None, nsub=8):
        self.c, self.h, self.V0, self.V1, self.U0, self.U1 = float(c), h, V0, V1, U0, U1
        self.kc = kc or KChart(c)
        self.Nv = int(round((V1 - V0) / h))
        self.Nu = int(round((U1 - U0) / h))
        par, tab = self.kc.par, self.kc.tab
        a0 = spec.get("a0", 0.0)
        Rr = spec.get("R", (0.0, 3.0, 0.5))
        Ll = spec.get("L", (0.0, 3.0, 0.5))
        self.UR = ray_dataK(par, tab, 0, U0, V0, self.Nv, h, Rr[0], Rr[1], Rr[2], a0, spec.get("a1R", 0.0), nsub)
        self.UL = ray_dataK(par, tab, 1, V0, U0, self.Nu, h, Ll[0], Ll[1], Ll[2], a0, spec.get("a1L", 0.0), nsub)
        self.ih = int(round(-U0 / h)) if U0 < 0 < U1 else None           # the horizon node U = 0
        assert self.ih is None or abs(U0 + self.ih * h) < 1e-9

    def run(self):
        t0 = time.time()
        self.U = march_K(self.UL, self.UR, self.V0, self.U0, self.kc.par, self.kc.tab, self.h, self.h)
        self.wall = time.time() - t0
        return self.U


# ----------------------------------------------------------------------------------------------------------------------
# c1, c2 on the horizon: the fit of the horizon values, and the functional of the data on a spacelike slice
# ----------------------------------------------------------------------------------------------------------------------
def horizon_fit(ah, V):
    """Least-squares a_h = c1 + c2 V over all nodes. Returns (c1, c2, max |a_h - fit| / max |a_h|)."""
    M = np.vstack([np.ones_like(V), V]).T
    co = np.linalg.lstsq(M, ah, rcond=None)[0]
    return co[0], co[1], float(np.max(np.abs(M @ co - ah)) / max(np.max(np.abs(ah)), 1e-300))


def fd4(f, j, i, dj, di, h):
    """4th-order central derivative of the grid function f along the node direction (dj, di), in units of the step h per unit of (dj, di)."""
    return (-f[j + 2 * dj, i + 2 * di] + 8 * f[j + dj, i + di] - 8 * f[j - dj, i - di] + f[j - 2 * dj, i - 2 * di]) / (12.0 * h)


def slice_V_derivative(f, j0, i0, h):
    """d_V f at the node (j0, i0) from the data of the spacelike slice Sigma = {V + U = const} through it and of its normal derivative:
    d_V = (T + N)/2 with T = d_V - d_U along the slice (values on Sigma) and N = d_V + d_U normal to it (values on Sigma and on the slices at +-1, +-2 nodes).
    The slice nodes are (j0 + k, i0 - k); the normal derivative uses (j0 + k, i0 + k)."""
    fT = fd4(f, j0, i0, 1, -1, h)
    fN = fd4(f, j0, i0, 1, 1, h)
    mut = MUT
    if mut == "fn_aT":
        return fT                                  # the tangential half only
    if mut == "fn_dirU":
        return 0.5 * (fN - fT)                     # d_U instead of d_V
    return 0.5 * (fT + fN)


def c2_from_slice(U, kc, j0, ih, h, V0, which="a"):
    """The functional: c2 = d_V Y on the future horizon at the point where the slice crosses it, from the slice data. which = 'a': d_V dr; 'Y': d_V (dr - g dphi),
    g = rho R/sqrt2 constant on the horizon (g_V = 0 there)."""
    x, d, r, rV, rU, pV, pU, V0_, V1_, V2_, E, FV, FU, rhoV, rhoU = jetsK(V0 + j0 * h, 0.0, kc.par, kc.tab, 0.0)
    aV = slice_V_derivative(U[0], j0, ih, h)
    if which == "a":
        return aV
    g = x * r / SQ2
    return aV - g * slice_V_derivative(U[2], j0, ih, h)


def richardson2(c_h, c_h2):
    return (4.0 * c_h2 - c_h) / 3.0


RESULTS = []


def want(cid):
    return ONLY is None or cid in ONLY or any(cid.startswith(o) for o in ONLY)


def check(cid, text, ok, detail=""):
    RESULTS.append((cid, bool(ok)))
    print(f"[{'pass' if ok else 'FAIL'}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


# ----------------------------------------------------------------------------------------------------------------------
# data classes, runs (cached), and the numerical checks N1 - N10
# ----------------------------------------------------------------------------------------------------------------------
CLASSES = {                                         # amplitudes 1e-3, as #232; corner = (V0, U0) of the rectangle
    "L": dict(L=(1e-3, 0.0, 0.3)),                  # pulse on the left ray V = V0, centred on the horizon U = 0
    "R": dict(R=(1e-3, 1.0, 0.3)),                  # pulse on the ray U = U0 (exterior), left ray zero
    "G": dict(L=(1e-3, 0.0, 0.3), R=(1e-3, 1.0, 0.3)),
    "A0": dict(a0=1e-3), "A1R": dict(a1R=1e-3), "A1L": dict(a1L=1e-3),
}
RECT = (-0.05, 3.0, -2.0, 0.06)                     # V0, V1, U0, U1: contains S, the exterior (U < 0) and the interior (U > 0, -Pi = V U <= 0.18)
HS = (0.01, 0.005, 0.0025, 0.00125)
SLICES = (0.0, 0.5, 1.0)                            # V where the slice V + U = const crosses the future horizon U = 0 (0 = the bifurcation surface S)
SAMPLES = {"S": (0.0, 0.0), "interior1": (1.0, 0.04), "interior2": (2.0, 0.06), "exterior": (1.0, -0.5), "left": (-0.03, 0.05)}
_CACHE = {}
_KC = {}


def kchart(c, **kw):
    key = (c, MUT in ("k_emF", "k_FV", "k_cross"), tuple(sorted(kw.items())))
    if key not in _KC:
        mut = {}
        if MUT == "k_cross":
            mut[lin.P_MCROSS] = -1.0
        if MUT == "k_emF":
            mut[P_MEF] = 1.01
        if MUT == "k_FV":
            mut[P_MFV] = -1.0
        if MUT == "k_con" and kw.get("data", True):
            pass
        _KC[key] = KChart(c, mut=mut)
    return _KC[key]


def run_class(c, name, h, spec=None, rect=RECT):
    key = (c, name, h, rect, MUT if MUT in ("k_emF", "k_FV", "k_cross", "k_con", "fn_aT", "fn_dirU") else None)
    if key in _CACHE:
        return _CACHE[key]
    kc = kchart(c)
    spec = spec if spec is not None else CLASSES[name]
    if MUT == "k_con":
        kc.par[lin.P_MCON] = -1.0                  # the ray constraint with the wrong sign: the data violate C_VV, C_UU
    su = KSetup(c, h, *rect, spec, kc=kc)
    if MUT == "k_con":
        kc.par[lin.P_MCON] = 1.0                   # the march itself is unmutated
    U = su.run()
    V = su.V0 + h * np.arange(su.Nv + 1)
    c1, c2, res = horizon_fit(U[0, :, su.ih], V)
    out = dict(c1=c1, c2=c2, res=res, slice={}, samples={}, wall=su.wall)
    for V0 in SLICES:
        j0 = int(round((V0 - su.V0) / h))
        out["slice"][V0] = c2_from_slice(U, kc, j0, su.ih, h, V0)
    for k, (Vp, Up) in SAMPLES.items():
        j, i = int(round((Vp - su.V0) / h)), int(round((Up - su.U0) / h))
        x, d, r, rV, rU, pV, pU, V0_, V1_, V2_, E, FV, FU, rhoV, rhoU = jetsK(Vp, Up, kc.par, kc.tab, 0.0)
        out["samples"][k] = (U[0, j, i], U[1, j, i], U[2, j, i], U[0, j, i] - x * r / SQ2 * U[2, j, i])
    _CACHE[key] = out
    del U
    return out


def conv_ratio(vals):
    """Cauchy ratios |v_h - v_h/2| / |v_h/2 - v_h/4| of a sequence at halving steps."""
    return [abs(vals[k] - vals[k + 1]) / max(abs(vals[k + 1] - vals[k + 2]), 1e-300) for k in range(len(vals) - 2)]


def check_N1():
    worst = 0.0
    for c in (-1.0, -1.5):
        kc = kchart(c)
        ch, b = kc.ch, kc.ch.b
        ss = ch.par[lin.P_SS]
        for p, sg in ((0.3, 1.0), (1.0, 3.0), (2.0, 10.0), (3.0, 15.9), (4.0, 16.0), (3.0, 16.3), (2.0, 16.6), (5.0, 12.0)):
            rho, d, rp, rs, emf, fp, fs = lin.bg(p, sg, ch.par, ch.tab, 0.0)
            Q = math.sinh(b * (ss - sg)) / math.sinh(b * ss)
            Qd = -b * math.cosh(b * (ss - sg)) / math.sinh(b * ss)
            V, U = math.exp(b * p), -Q
            x, dd, r, rV, rU, pV, pU, V0_, V1_, V2_, E, FV, FU, rhoV, rhoU = jetsK(V, U, kc.par, kc.tab, 0.0)
            # (p, sigma) jets from the Kruskal ones: d_p = b V d_V, d_sigma = -Qd d_U (U = -Q), f = F - ln(b V |Q'|)
            ref = [rho, emf, rp, rs, fp, fs]
            Qdd = b * b * math.sinh(b * (ss - sg)) / math.sinh(b * ss)           # f_sigma = F_U dU/dsigma - Q''/Q'
            got = [x, E * b * V * (-Qd), rhoV * b * V, rhoU * (-Qd), FV * b * V - b, FU * (-Qd) - Qdd / Qd]
            for a_, b_ in zip(got, ref):
                worst = max(worst, abs(a_ - b_) / (1e-3 + abs(b_)))
    check("N1", "the Kruskal-chart jets (rho, e^{-F}, r_V, r_U, F_V, F_U) reproduce linear.bg through V = e^{b p}, U = -Q(sigma) at 8 points each of c = -1 and -1.5 (static region, horizon, interior)",
          worst < 1e-12, f"max relative difference {worst:.1e}")


def check_N2():
    out, ok = [], True
    for tag, c, schw in (("c=-1", -1.0, False), ("c=-1.5", -1.5, False), ("Schwarzschild", -1.0, True)):
        kc = kchart(c)
        if schw:
            kc = KChart(c)
            kc.par[P_KIND2] = 1.0
            if MUT == "k_FV":
                kc.par[P_MFV] = -1.0
        errs, rels = [], []
        for h in (0.02, 0.01, 0.005):
            V0, V1, U0, U1 = -0.05, 2.0, -1.5, 0.05
            Nv, Nu = int(round((V1 - V0) / h)), int(round((U1 - U0) / h))
            xi_f = lambda V: (0.01 * math.sin(0.7 * V + 0.3), 0.007 * math.cos(0.7 * V + 0.3))
            eta_f = lambda U: (0.02 * math.cos(0.5 * U) + 0.01, -0.01 * math.sin(0.5 * U))
            Uex = gauge_arrays_K(kc, V0, U0, h, h, Nv, Nu, xi_f, eta_f)
            Un = march_K(np.ascontiguousarray(Uex[:, 0, :]), np.ascontiguousarray(Uex[:, :, 0]), V0, U0, kc.par, kc.tab, h, h)
            errs.append(float(np.abs(Un - Uex).max()))
            rels.append(errs[-1] / float(np.abs(Uex).max()))
        ratios = [errs[0] / errs[1], errs[1] / errs[2]]
        good = 3.7 < ratios[1] < 4.3 and rels[2] < 5e-4
        ok = ok and good
        out.append(f"{tag}: {errs[0]:.2e}, {errs[1]:.2e}, {errs[2]:.2e} (ratios {ratios[0]:.2f}, {ratios[1]:.2f}; relative to the largest field {rels[2]:.1e})")
    check("N2", "exact pure-gauge solutions L_xi(member), xi(V), eta(U) smooth, reproduced by the Kruskal-chart scheme at second order through S and the interior (h = 0.02, 0.01, 0.005; the last ratio in [3.7, 4.3], error below 5e-4 of the largest field at 0.005): c = -1, -1.5, Schwarzschild",
          ok, "; ".join(out))


def check_N3():
    out, ok = [], True
    for name in ("L", "G", "A0", "A1R", "A1L"):
        res = [run_class(-1.0, name, h)["res"] for h in HS[1:]]
        ratios = [res[k] / res[k + 1] for k in range(2)]
        good = all(3.5 < q < 4.5 for q in ratios) and res[-1] < 2e-5
        ok = ok and good
        out.append(f"{name}: {res[0]:.1e}, {res[1]:.1e}, {res[2]:.1e} ({ratios[0]:.2f}, {ratios[1]:.2f})")
    check("N3", "c = -1, data through S: the horizon values dr_h(V) = Y_h(V) are linear in V on the whole horizon V in [-0.05, 3] (least-squares residual relative to max |dr_h| at h = 0.005, 0.0025, 0.00125 below 2e-5 at the finest, falling by 4 +- 0.5 per halving): classes L, G, A0, A1R, A1L",
          ok, "; ".join(out))


def check_N4():
    out, ok = [], True
    for name in ("L", "G", "A0", "A1R", "A1L"):
        runs = [run_class(-1.0, name, h) for h in HS]
        rf = [richardson2(runs[k]["c2"], runs[k + 1]["c2"]) for k in range(3)]
        worst = 0.0
        txt = []
        for V0 in SLICES:
            rs = [richardson2(runs[k]["slice"][V0], runs[k + 1]["slice"][V0]) for k in range(3)]
            dev = abs(rs[-1] - rf[-1]) / abs(rf[-1])
            worst = max(worst, dev)
            txt.append(f"{dev:.1e}")
        good = worst < 1e-6
        ok = ok and good
        out.append(f"{name}: c2 = {rf[-1]:.8e}, slices V = 0, 0.5, 1: " + ", ".join(txt))
    check("N4", "c = -1: c2 from the data of the spacelike slice V + U = const crossing the future horizon at V = 0 (through S), 0.5, 1.0 (4th-order d_V of dr from the slice values and its normal derivative) equals c2 of the horizon fit to 1e-6, Richardson in h (0.005, 0.0025, 0.00125)",
          ok, "; ".join(out))


def lin_L_runs(s0, hs_, P=10.0, c=-1.0):
    """#232's own linear solution at c = -1 (linear.py, sigma* = 16, left-ray pulse (1e-3, s0, 0.5), the setup of its check L3a), on p in [0, P]:
    the full array and c2 of the horizon fit (coefficient of e^{b p})."""
    b = lin.ep.Member(c).kappa * SQ2
    out = {}
    for h in hs_:
        key = ("232", c, s0, h, P, MUT if MUT in ("fn_aT", "fn_dirU", "fn_b") else None)
        if key not in _CACHE:
            su = lin.mksetup(c, h, P, dict(L=(1e-3, s0, 0.5)))
            U = lin.march_full(su.UL, su.UR, su.chart.par, su.chart.tab, su.hp, su.hs)
            ih = su.jh
            p = h * np.arange(su.Np + 1)
            c1, c2, res = horizon_fit(U[0, :, ih], np.exp(b * p))
            two = lin.c12_from(dict(Yh=U[0, :, ih]), b, 0.0, 2.0, h)[1]
            bb = b * (1.01 if MUT == "fn_b" else 1.0)
            sl = {pj: slice_V_derivative(U[0], int(round(pj / h)), ih, h) / (bb * math.exp(bb * pj)) for pj in (1.0, 2.0, 4.0, 6.0)}
            _CACHE[key] = dict(c2=c2, two=two, res=res, slice=sl, mem=U.nbytes / 1e6)
            del U
        out[h] = _CACHE[key]
    return out


def check_N5():
    out, ok = [], True
    hs_ = (0.01, 0.005, 0.0025)
    for s0 in (2.5, 3.5, 4.5):
        r = lin_L_runs(s0, hs_)
        rf = richardson2(r[0.005]["c2"], r[0.0025]["c2"])
        r2 = richardson2(r[0.005]["two"], r[0.0025]["two"])
        devs = [abs(richardson2(r[0.005]["slice"][pj], r[0.0025]["slice"][pj]) - rf) / abs(rf) for pj in (1.0, 2.0, 4.0, 6.0)]
        raw = [abs(r[0.005]["slice"][pj] - rf) / abs(rf) for pj in (1.0, 2.0, 4.0, 6.0)]
        dev2 = abs(rf - r2) / abs(rf)
        good = max(devs) < 1e-6 and dev2 < 1e-6
        ok = ok and good
        out.append(f"s0 = {s0}: c2 = {rf:.9e}; slices at p = 1, 2, 4, 6 deviate by " + ", ".join(f"{d:.1e}" for d in devs) + f" (raw at h = 0.005: {max(raw):.1e}); #232's two-point c2 {dev2:.1e}; {r[0.0025]['mem']:.0f} MB")
    check("N5", "#232's own linear solutions at c = -1 (linear.py, L3a's three left-ray pulses, h = 0.005, 0.0025 Richardson): c2 from the data of the spacelike slices p + sigma = const crossing the horizon at p = 1, 2, 4, 6 equals c2 of the horizon fit and #232's two-point c2 to 1e-6",
          ok, "; ".join(out))


def mp_G_horizon(c):
    """G(rho_h) = E(rho_h) A'(rho_h)/(2 kappa) = E(rho_h) at c < -1 (kappa = A'/2), E(rho_c) = rho_c - rho_h: mpmath quadrature of ln E' = 1/delta - 2 kappa/A
    (60 digits, from delta = 1e-12 with the limit g(rho_h) = A''/(2 A') on the first sliver). Returns (rho_h, kappa, G(rho_h))."""
    import mpmath as mp
    mp.mp.dps = 60
    cc = mp.mpf(c)
    K = -2 * cc / mp.pi
    A = lambda r: 1 + cc * (r * r + 1) + K * (r + (r * r + 1) * mp.atan(r))
    Ap = lambda r: 2 * cc * r + 2 * K * (1 + r * mp.atan(r))
    App = lambda r: 2 * cc + 2 * K * (mp.atan(r) + r / (1 + r * r))
    rh = mp.findroot(A, mp.mpf("0.34"))
    kap = Ap(rh) / 2
    g = lambda r: 1 / (r - rh) - 2 * kap / A(r)
    eps = mp.mpf("1e-12")
    integral = mp.quad(g, [rh + eps, (rh + 1) / 2, 1]) + eps * App(rh) / (2 * Ap(rh))
    return float(rh), float(kap), float(mp.exp(mp.log(1 - rh) - integral))


def check_N6():
    """Planted: pure gauge L_xi(member) has a known dr_h(V) = -eta(0) r_U-coefficient V on the horizon. c = -1.5 (rho_h != 0) and Schwarzschild; the Y-version is zero."""
    out, ok = [], True
    eta0 = 0.03
    xi_f = lambda V: (0.01 * math.sin(0.7 * V + 0.3), 0.007 * math.cos(0.7 * V + 0.3))
    eta_f = lambda U: (0.02 * math.cos(0.5 * U) + 0.01, -0.01 * math.sin(0.5 * U))
    rh, kap, Gh = mp_G_horizon(-1.5)
    Rp = rh / math.sqrt(1 + rh * rh)
    sgn = -1.0 if MUT == "plant_sign" else 1.0
    cases = [("c=-1.5", -1.5, False, sgn * (-eta0 * Rp * Gh), 0.0), ("Schwarzschild", -1.0, True, sgn * (-eta0 * math.exp(-1.0)), None)]
    for tag, c, schw, want_a, want_Y in cases:
        kc = KChart(c) if schw else kchart(c)
        if schw:
            kc.par[P_KIND2] = 1.0
        V0, V1, U0, U1 = -0.05, 2.0, -1.5, 0.05
        ex, mar = [], []
        for h in (0.01, 0.005, 0.0025):
            Nv, Nu = int(round((V1 - V0) / h)), int(round((U1 - U0) / h))
            Uex = gauge_arrays_K(kc, V0, U0, h, h, Nv, Nu, xi_f, eta_f)
            Un = march_K(np.ascontiguousarray(Uex[:, 0, :]), np.ascontiguousarray(Uex[:, :, 0]), V0, U0, kc.par, kc.tab, h, h)
            ih = int(round(-U0 / h))
            jj = int(round((0.5 - V0) / h))
            if schw:
                x, d, r, rV, rU, pV, pU, V0_, V1_, V2_, E, FV, FU, rhoV, rhoU = jetsK(0.5, 0.0, kc.par, kc.tab, 0.0)
                g_h = 0.0
            else:
                g_h = rh * math.sqrt(1 + rh * rh) / SQ2
            fa = slice_V_derivative(Uex[0], jj, ih, h)
            fY = fa - g_h * slice_V_derivative(Uex[2], jj, ih, h)
            ma = slice_V_derivative(Un[0], jj, ih, h)
            mY = ma - g_h * slice_V_derivative(Un[2], jj, ih, h)
            ex.append((fa, fY))
            mar.append((ma, mY))
        da_ex = abs(ex[-1][0] - want_a) / abs(want_a)
        rich_a = richardson2(mar[1][0], mar[2][0])
        da_mar = abs(rich_a - want_a) / abs(want_a)
        rich_Y = richardson2(mar[1][1], mar[2][1])
        if want_Y is not None:
            dY_ex, dY_mar = abs(ex[-1][1]) / abs(want_a), abs(rich_Y) / abs(want_a)
        else:
            dY_ex = dY_mar = 0.0
        good = da_ex < 1e-8 and da_mar < 1e-6 and dY_ex < 1e-8 and dY_mar < 1e-6
        ok = ok and good
        out.append(f"{tag}: known dr_h coefficient {want_a:.8e}; exact data {da_ex:.1e}, marched (Richardson) {da_mar:.1e}; Y-version {dY_ex:.1e}, {dY_mar:.1e} of it")
    # Schwarzschild test field: dr = dF = 0 stays zero, dphi evolves; its horizon values are not linear in V (no h^2 convergence to a linear law)
    kc = KChart(-1.0)
    kc.par[P_KIND2] = 1.0
    amax, rese = 0.0, []
    for h in (0.01, 0.005):
        su = KSetup(-1.0, h, -0.05, 2.0, -1.5, 0.05, dict(R=(1e-3, 1.0, 0.3)), kc=kc)
        U = su.run()
        Vv = su.V0 + h * np.arange(su.Nv + 1)
        rese.append(horizon_fit(U[2, :, su.ih], Vv)[2])
        amax = max(amax, float(np.abs(U[0]).max()))
    out.append(f"Schwarzschild test field (right-ray pulse): max |dr| = {amax:.1e} (decoupled); the horizon values of dphi have a linear-fit residual {rese[0]:.2f}, {rese[1]:.2f} of their maximum at h = 0.01, 0.005 (no linear law)")
    ok = ok and amax < 1e-14 and min(rese) > 0.1 and rese[0] / rese[1] < 1.5
    check("N6", "planted: pure-gauge data have a known horizon coefficient: dr_h = -eta(0) R'(rho_h) G(rho_h) V (c = -1.5; G from an mpmath quadrature) and -eta(0) e^{-1} V (Schwarzschild); the functional gives it to 1e-8 on exact data and 1e-6 on the marched solution, and the Y-version gives zero. A Schwarzschild test field leaves dr = 0 and has no linear horizon law (fit residual above 0.1 of the maximum at h = 0.01 and 0.005, not falling)",
          ok, "; ".join(out))


def check_N7():
    """Control: data with c2 = 0 give zero."""
    out, ok = [], True
    # (i) exact pure gauge at c = -1: dr vanishes to first order on the horizon, c2 = 0 exactly; the functional's 4th-order stencil leaves h^4
    kc = kchart(-1.0)
    V0, V1, U0, U1 = -0.05, 2.0, -1.5, 0.05
    xi_f = lambda V: (0.01 * math.sin(0.7 * V + 0.3), 0.007 * math.cos(0.7 * V + 0.3))
    eta_f = lambda U: (0.02 * math.cos(0.5 * U) + 0.01, -0.01 * math.sin(0.5 * U))
    rel_g = []
    for h in (0.01, 0.005):
        Nv, Nu = int(round((V1 - V0) / h)), int(round((U1 - U0) / h))
        Uex = gauge_arrays_K(kc, V0, U0, h, h, Nv, Nu, xi_f, eta_f)
        ih = int(round(-U0 / h))
        vals = [abs(slice_V_derivative(Uex[0], int(round((Vc - V0) / h)), ih, h)) for Vc in (0.0, 0.5, 1.0)]
        rel_g.append(max(vals) / float(np.abs(Uex[0]).max()))
    out.append(f"exact pure gauge at c = -1: |c2| / max |dr| = {rel_g[0]:.1e}, {rel_g[1]:.1e} at h = 0.01, 0.005 (ratio {rel_g[0] / rel_g[1]:.1f})")
    ok = ok and rel_g[1] < 1e-5 and 10.0 < rel_g[0] / rel_g[1] < 22.0
    # (ii) class R (right-ray pulse, left ray zero): zero in the continuum, the truncation seed falls by 4 per halving
    rR = [run_class(-1.0, "R", hh) for hh in HS]
    rL = [run_class(-1.0, "L", hh) for hh in HS]
    cL = richardson2(rL[2]["slice"][0.5], rL[3]["slice"][0.5])
    seeds = [abs(r["slice"][0.5]) for r in rR]
    ratios = [seeds[k] / seeds[k + 1] for k in range(3)]
    rel_R = abs(rR[3]["slice"][0.5]) / abs(cL)
    out.append(f"class R: |c2| at h = 0.01 ... 0.00125 = " + ", ".join(f"{x:.1e}" for x in seeds) + f" (ratios {ratios[0]:.2f}, {ratios[1]:.2f}, {ratios[2]:.2f}); the finest is {rel_R:.1e} of class L's c2")
    ok = ok and rel_R < 1e-5 and all(3.0 < q < 5.0 for q in ratios)
    # (iii) a superposition with c2 = 0 by construction: lambda from the horizon fits (Richardson) of L and A1R, then the functional on the slice data of the sum
    cfit = {nm: richardson2(run_class(-1.0, nm, HS[2])["c2"], run_class(-1.0, nm, HS[3])["c2"]) for nm in ("L", "A1R")}
    lam = -cfit["L"] / cfit["A1R"] * (1.1 if MUT == "ctl_lambda" else 1.0)
    spec = dict(L=(1e-3, 0.0, 0.3), a1R=lam * 1e-3)
    sums = [run_class(-1.0, "tuned", hh, spec=spec) for hh in HS[2:]]
    res = {}
    for V0c in SLICES:
        res[V0c] = abs(richardson2(sums[0]["slice"][V0c], sums[1]["slice"][V0c])) / abs(cfit["L"])
    fitres = abs(richardson2(sums[0]["c2"], sums[1]["c2"])) / abs(cfit["L"])
    out.append(f"tuned L + {lam:.6f} A1R: |c2| / |c2 of L| = {fitres:.1e} (fit), " + ", ".join(f"{res[v]:.1e}" for v in SLICES) + " (slices V = 0, 0.5, 1)")
    ok = ok and max(res.values()) < 1e-5 and fitres < 1e-5
    check("N7", "control, c2 = 0 by construction: exact pure gauge at c = -1 gives |c2| / max |dr| below 1e-5 at h = 0.005, falling by 16 (10 to 22) from h = 0.01 (the stencil's h^4); the right-ray class R gives |c2| below 1e-5 of class L's at h = 0.00125 with the seed falling by 4 (3 to 5) per halving; the superposition L + lambda A1R with lambda from the horizon fits gives |c2| below 1e-5 of class L's on the slices V = 0, 0.5, 1",
          ok, "; ".join(out))


def check_N8():
    """Regularity: the fields at S, in the interior and on the other side converge at second order."""
    out, ok = [], True
    for name in ("G", "A0"):
        runs = [run_class(-1.0, name, hh) for hh in HS]
        for k in SAMPLES:
            for q, qn in enumerate(("dr", "dF", "dphi", "Y")):
                vals = [r["samples"][k][q] for r in runs]
                if max(abs(v) for v in vals) < 1e-9:
                    continue
                rat = conv_ratio(vals)
                good = 3.5 < rat[1] < 4.5 and 3.0 < rat[0] < 5.5
                ok = ok and good
                if not good or (name == "G" and qn == "Y"):
                    out.append(f"{name} {k} {qn}: ratios {rat[0]:.2f}, {rat[1]:.2f}")
    check("N8", "c = -1, classes G and A0: dr, dF, dphi and Y at S, two points in the interior (-Pi = V U = 0.04, 0.12), the static region and the left side converge at second order (Cauchy ratios over h = 0.01 ... 0.00125: 3.0 to 5.5 for the coarse pair, 3.5 to 4.5 for the fine pair)",
          ok, "; ".join(out))


def check_N10():
    """The data are physical: the linearized constraints C_VV, C_UU of the evolved solution (second-order differences, every 0.05 of V and U, through S and the interior) converge at second order."""
    kc = kchart(-1.0)
    out, ok = [], True
    for name in ("G", "A0"):
        vals = []
        for h in (0.01, 0.005):
            su = KSetup(-1.0, h, *RECT, CLASSES[name], kc=kc)
            U = su.run()
            st = int(round(0.05 / h))
            cvv, cuu, sc = [], [], 1e-300
            for j in range(2 * st, su.Nv - 2, st):
                for i in range(2 * st, su.Nu - 2, st):
                    x, d, r, rV, rU, pV, pU, V0_, V1_, V2_, E, FV, FU, rhoV, rhoU = jetsK(su.V0 + j * h, su.U0 + i * h, kc.par, kc.tab, 0.0)
                    a, dd, e = U[0], U[1], U[2]
                    aV = (a[j + 1, i] - a[j - 1, i]) / (2 * h)
                    aVV = (a[j + 1, i] - 2 * a[j, i] + a[j - 1, i]) / h ** 2
                    dV_ = (dd[j + 1, i] - dd[j - 1, i]) / (2 * h)
                    eV = (e[j + 1, i] - e[j - 1, i]) / (2 * h)
                    aU = (a[j, i + 1] - a[j, i - 1]) / (2 * h)
                    aUU = (a[j, i + 1] - 2 * a[j, i] + a[j, i - 1]) / h ** 2
                    dU_ = (dd[j, i + 1] - dd[j, i - 1]) / (2 * h)
                    eU = (e[j, i + 1] - e[j, i - 1]) / (2 * h)
                    cvv.append(aVV + FV * aV + rV * dV_ - 0.5 * pV ** 2 * a[j, i] - r * pV * eV)
                    cuu.append(aUU + FU * aU + rU * dU_ - 0.5 * pU ** 2 * a[j, i] - r * pU * eU)
                    sc = max(sc, abs(aVV), abs(FV * aV), abs(aUU), abs(FU * aU))
            vals.append((float(np.max(np.abs(cvv))), float(np.max(np.abs(cuu))), sc))
        rat = [vals[0][0] / vals[1][0], vals[0][1] / vals[1][1]]
        rel = max(vals[1][0], vals[1][1]) / vals[1][2]
        good = all(3.3 < q < 4.7 for q in rat) and rel < 5e-3
        ok = ok and good
        out.append(f"{name}: C_VV {vals[0][0]:.1e} -> {vals[1][0]:.1e}, C_UU {vals[0][1]:.1e} -> {vals[1][1]:.1e} (ratios {rat[0]:.2f}, {rat[1]:.2f}; relative to the largest term {rel:.1e})")
    check("N10", "the data are physical: the linearized constraints C_VV, C_UU of the evolved Kruskal-chart solution (classes G and A0, through S and the interior) shrink by a factor in [3.3, 4.7] from h = 0.01 to 0.005 and are below 5e-3 of the largest term at 0.005",
          ok, "; ".join(out))


def check_N11():
    """Generic data: random smooth data on the two rays (random corner constants and pulses) have c2 != 0, and the slice functional equals the horizon fit on each."""
    rng = np.random.default_rng(248)
    kc = kchart(-1.0)
    n, devs, c2s = 20, [], []
    for k in range(n):
        spec = dict(a0=rng.uniform(-1e-3, 1e-3), a1R=rng.uniform(-1e-3, 1e-3), a1L=rng.uniform(-1e-3, 1e-3),
                    L=(rng.uniform(-1e-3, 1e-3), rng.uniform(-1.2, 0.4), rng.uniform(0.2, 0.5)), R=(rng.uniform(-1e-3, 1e-3), rng.uniform(0.3, 2.5), rng.uniform(0.2, 0.5)))
        runs = [run_class(-1.0, f"rand{k}", h, spec=spec) for h in HS[2:]]
        cf = richardson2(runs[0]["c2"], runs[1]["c2"])
        dev = max(abs(richardson2(runs[0]["slice"][v], runs[1]["slice"][v]) - cf) for v in (0.0, 0.5))
        devs.append(dev / max(abs(cf), 1e-3))
        c2s.append(cf)
    c2s = np.array(c2s)
    ok = max(devs) < 1e-6 and float(np.min(np.abs(c2s))) > 1e-5
    check("N11", "generic data: for 20 random smooth data (random corner constants dr, dr_V, dr_U and random Gaussian pulses of dphi on both rays, amplitudes up to 1e-3) c2 from the slices V = 0 and 0.5 equals the horizon fit to 1e-6 of max(|c2|, 1e-3), and no draw has |c2| below 1e-5",
          ok, f"worst slice-fit deviation {max(devs):.1e}; |c2| from {float(np.min(np.abs(c2s))):.2e} to {float(np.max(np.abs(c2s))):.2e}; median {float(np.median(np.abs(c2s))):.2e}")


def check_N9():
    """The Kruskal-chart solver against linear.py's solver (a different chart, grid and background code) on the same data: corner constants at (V = 1, U = -1) = (p = 0, sigma = 0)."""
    out, ok = [], True
    c = -1.0
    mem = lin.ep.Member(c)
    b = mem.kappa * SQ2
    ss = 16.0
    coth = math.cosh(b * ss) / math.sinh(b * ss)
    cases = {"A0": (dict(a0=1e-3), dict(a0=1e-3)), "A1R": (dict(a1R=1e-3), dict(a1R=1e-3 / b)), "A1L": (dict(a1L=1e-3), dict(a1L=1e-3 / (b * coth)))}
    P = 3.0
    for name, (spec_lin, spec_K) in cases.items():
        cl = {}
        for h in (0.01, 0.005, 0.0025):
            su = lin.mksetup(c, h, P, spec_lin)
            U = lin.march_full(su.UL, su.UR, su.chart.par, su.chart.tab, su.hp, su.hs)
            p = h * np.arange(su.Np + 1)
            cl[h] = horizon_fit(U[0, :, su.jh], np.exp(b * p))[1]
            del U
        c_lin = richardson2(cl[0.005], cl[0.0025])
        ck = {}
        for h in (0.005, 0.0025, 0.00125):
            kc = kchart(c)
            su = KSetup(c, h, 1.0, math.exp(b * P), -1.0, 0.02, spec_K, kc=kc)
            U = su.run()
            V = su.V0 + h * np.arange(su.Nv + 1)
            ck[h] = horizon_fit(U[0, :, su.ih], V)[1]
            del U
        c_K = richardson2(ck[0.0025], ck[0.00125])
        dev = abs(c_K - c_lin) / abs(c_lin)
        ok = ok and dev < 1e-5
        out.append(f"{name}: c2 = {c_lin:.8e} (linear.py), {c_K:.8e} (Kruskal), relative {dev:.1e}")
    check("N9", "the Kruskal-chart solver and linear.py's solver, on the same corner-constant data at (V, U) = (1, -1) = (p, sigma) = (0, 0), give the same c2 (Richardson of three resolutions each) to 1e-5, p in [0, 3]",
          ok, "; ".join(out))


# ----------------------------------------------------------------------------------------------------------------------
# the derivation (sympy / mpmath): D1 - D8
# ----------------------------------------------------------------------------------------------------------------------
def _member_sym(c):
    import sympy as sp
    rho = sp.Symbol("rho", real=True)
    c = sp.nsimplify(c)
    K = -2 * c / sp.pi
    A = 1 + c * (rho ** 2 + 1) + K * (rho + (rho ** 2 + 1) * sp.atan(rho))
    R = sp.sqrt(rho ** 2 + 1)
    Ph = SQ2 * sp.atan(rho) if False else sp.sqrt(2) * sp.atan(rho)
    psi = sp.atan(rho)
    Vp = -c * (2 - sp.cos(2 * psi)) - K * (sp.Rational(3, 2) * sp.sin(2 * psi) + psi * (2 - sp.cos(2 * psi)))
    Vp1 = sp.diff(Vp, rho) / sp.diff(Ph, rho)
    return rho, A, R, Ph, Vp, Vp1


def derive_D1():
    """The member in the Kruskal chart solves the nonlinear equations: e^{-F} = G/kappa, d rho/d Pi = G, Pi = A/(2 kappa G); residual identically zero for symbolic kappa."""
    import mpmath as mp_
    import sympy as sp
    V, U = sp.symbols("V U", real=True)
    kap = sp.Symbol("kappa", positive=True)
    worst = []
    for c in (-1, sp.Rational(-3, 2), sp.Rational(-3, 10)):
        rho, A, R, Ph, Vp, Vp1 = _member_sym(c)
        G = sp.Function("G")(rho)
        L = (sp.diff(A, rho) - 2 * kap) / A * (1.01 if MUT == "d_Gode" else 1)          # d ln G / d rho = (A' - 2 kappa)/A  (G = A/(2 kappa Pi))

        def sub(e):
            for _ in range(3):
                e = e.subs(sp.Derivative(G, (rho, 3)), sp.diff(G * L, rho, 2)).subs(sp.Derivative(G, (rho, 2)), sp.diff(G * L, rho)).subs(sp.Derivative(G, rho), G * L)
            return e
        DV = lambda f: sub(sp.diff(f, rho) * (-U * G) + sp.diff(f, V))
        DU = lambda f: sub(sp.diff(f, rho) * (-V * G) + sp.diff(f, U))
        emF = G / kap * (sp.Rational(101, 100) if MUT == "d_emF" else 1)
        F_ = sp.log(kap) - sp.log(G)
        rV, rU, pV, pU = DV(R), DU(R), DV(Ph), DU(Ph)
        res = [DU(rV) - (-(rV * rU + emF / 2) / R + Vp * R * emF / 2),
               DU(DV(F_)) - (-2 * (rV * rU + emF / 2) / R ** 2 - pV * pU),
               DU(pV) - (-(rV * pU + rU * pV) / R + sp.Rational(1, 2) * emF * Vp1),
               DV(rV) + DV(F_) * rV - R / 2 * pV ** 2,
               DU(rU) + DU(F_) * rU - R / 2 * pU ** 2]
        for e in res:
            e = sp.expand(e).subs(U, -A / (2 * kap * G * V))            # U V = -Pi = -A/(2 kappa G)
            # G(rho) and rho are independent at a point (G carries the free normalization of Pi): a numeric test first, the exact simplification when it is zero
            mp_.mp.dps = 50
            gs, kv = sp.Symbol("Gv", positive=True), sp.Rational(3, 5)
            en = e.subs(G, gs).subs(kap, kv)
            fn = sp.lambdify((rho, V, gs), en, "mpmath")
            numeric = max(abs(fn(mp_.mpf(r_), mp_.mpf(v_), mp_.mpf(g_))) for r_, v_, g_ in (('0.37', '1.3', '0.8'), ('-0.6', '0.7', '1.9'), ('1.4', '-0.4', '0.3')))
            worst.append(numeric > 1e-25 or sp.simplify(e) != 0)
    check("D1", "the Kruskal-chart member (r = R(rho), phi = sqrt2 arctan rho, e^{-F} = G/kappa, rho_V = -U G, rho_U = -V G, G = A/(2 kappa Pi)) solves the three evolution equations and both constraints, residual identically zero (sympy, symbolic kappa): c = -1, -3/2, -3/10",
          not any(worst), f"{sum(worst)} of {len(worst)} residuals nonzero")


def derive_D2():
    """Horizon law: on U = 0, linearized C_VV is a_VV = 0 (r_V = phi_V = F_V = 0 there); at c = -1 also g = rho R/sqrt2 = 0, so Y_h = dr_h; the generic-jet propagation identity."""
    import sympy as sp
    V, U = sp.symbols("V U", real=True)
    kap = sp.Symbol("kappa", positive=True)
    rho, A, R, Ph, Vp, Vp1 = _member_sym(-1)
    G = sp.Function("G")(rho)
    L = (sp.diff(A, rho) - 2 * kap) / A
    sub = lambda e: e.subs(sp.Derivative(G, rho), G * L)
    DV = lambda f: sub(sp.diff(f, rho) * (-U * G) + sp.diff(f, V))
    rV, pV, FV = DV(R), DV(Ph), DV(sp.log(kap) - sp.log(G))
    a, d, e = (sp.Function(n)(V, U) for n in "ade")
    eps = -1
    CVV = sp.diff(a, V, 2) + FV * sp.diff(a, V) + rV * sp.diff(d, V) + sp.Rational(eps, 2) * pV ** 2 * a + eps * R * pV * sp.diff(e, V)
    on_h = sp.simplify(CVV.subs(U, sp.Rational(1, 10) if MUT == "d_hz" else 0).doit())
    ok1 = sp.simplify(on_h - sp.diff(a, V, 2).subs(U, 0)) == 0
    g = rho * R / sp.sqrt(2)
    ok2 = sp.simplify(g.subs(rho, 0)) == 0
    check("D2", "on the future horizon U = 0 of the c = -1 member (rho = 0) the linearized constraint C_VV is a_VV = 0: dr_h = c1 + c2 V; and g = rho R/sqrt2 = 0 there, so Y_h = dr_h (sympy)", ok1 and ok2, "")
    # propagation of the constraints (generic jets, exact rational check)
    ok3 = _propagation(MUT == "d_prop")
    check("D4", "the linearized constraints propagate: d_U C_VV = -(r_U/r) C_VV and d_V C_UU = -(r_V/r) C_UU on solutions of the evolution equations, for any background solving the field equations (generic jets, sympy, zero residual): the first-order jet at a point is not restricted", ok3[0], ok3[1])


def _propagation(mutate):
    """Generic-jet identity d_U C_VV = -(r_U/r) C_VV, d_V C_UU = -(r_V/r) C_UU on solutions (exact rational check at random jets)."""
    import sympy as sp
    import random
    import sys
    argv, sys.argv = sys.argv, [sys.argv[0]]
    try:
        import linear_derive as LD
    finally:
        sys.argv = argv
    EPS = -1
    r, rV, rU, E, FV, FU, pV, pU, V0, V1, V2, V3 = sp.symbols("r rV rU E FV FU pV pU V0 V1 V2 V3")
    pVV, FVV, pUU, FUU = sp.symbols("pVV FVV pUU FUU")
    rUV = -(rV*rU + E/2)/r + V0*r*E/2
    FUVx = -2*(rV*rU + E/2)/r**2 + EPS*pV*pU
    pUVx = -(rV*pU + rU*pV)/r - sp.Rational(EPS, 2)*E*V1
    rVV = -FV*rV - sp.Rational(EPS, 2)*r*pV**2
    rUU = -FU*rU - sp.Rational(EPS, 2)*r*pU**2
    dV = {r: rV, rV: rVV, rU: rUV, E: -FV*E, FV: FVV, FU: FUVx, pV: pVV, pU: pUVx, V0: V1*pV, V1: V2*pV, V2: V3*pV}
    dU = {r: rU, rV: rUV, rU: rUU, E: -FU*E, FV: FUVx, FU: FUU, pV: pUVx, pU: pUU, V0: V1*pU, V1: V2*pU, V2: V3*pU}
    a, aV, aU, aVV, aUU, d, dV_, dU_, dVV, dUU, e, eV, eU, eVV, eUU = sp.symbols("a aV aU aVV aUU d dV_ dU_ dVV dUU e eV eU eVV eUU")
    back = dict(zip((LD.r, LD.rp, LD.rm, LD.E, LD.php, LD.phm, LD.V0, LD.V1, LD.V2), (r, rV, rU, E, pV, pU, V0, V1, V2)))
    B0, Bp, Bm = [Mx.subs(back, simultaneous=True) for Mx in (LD.B0, LD.Bp, LD.Bm)]
    Uv, UV_, UU_ = sp.Matrix([a, d, e]), sp.Matrix([aV, dV_, eV]), sp.Matrix([aU, dU_, eU])
    rhs = B0*Uv + Bp*UV_ + Bm*UU_
    fV = {a: aV, aV: aVV, aU: rhs[0], d: dV_, dV_: dVV, dU_: rhs[1], e: eV, eV: eVV, eU: rhs[2]}
    fU = {a: aU, aV: rhs[0], aU: aUU, d: dU_, dV_: rhs[1], dU_: dUU, e: eU, eV: rhs[2], eU: eUU}
    D = lambda expr, bgt, ft: sum(sp.diff(expr, k)*v for k, v in list(bgt.items()) + list(ft.items()))
    fU[aVV] = D(rhs[0], dV, fV)
    fV[aUU] = D(rhs[0], dU, fU)
    eps_c = -EPS if mutate else EPS
    CVV = aVV + FV*aV + rV*dV_ + sp.Rational(EPS, 2)*pV**2*a + eps_c*r*pV*eV
    CUU = aUU + FU*aU + rU*dU_ + sp.Rational(EPS, 2)*pU**2*a + eps_c*r*pU*eU
    out = []
    for C, bg_, ft, key, lam in ((CVV, dU, fU, aVV, -rU/r), (CUU, dV, fV, aUU, -rV/r)):
        res = sp.expand(D(C, bg_, ft) - lam*C)
        syms = sorted(res.free_symbols, key=str)
        random.seed(3)
        worst = 0
        for _ in range(4):
            vals = {s_: sp.Rational(random.randint(1, 9), random.randint(1, 9))*random.choice([1, -1]) for s_ in syms}
            worst = max(worst, abs(res.subs(vals)))
        out.append(worst)
    return (max(out) == 0, f"residuals {out[0]}, {out[1]} at 4 random rational points each (exact)")



def derive_D3():
    """Y is invariant under every linear diffeomorphism; at c = -1 the horizon is r-stationary, so c2 is not a gauge artifact; at c != -1 it is (planted)."""
    import sympy as sp
    rho = sp.Symbol("rho", real=True)
    R = sp.sqrt(rho ** 2 + 1)
    Phi = sp.sqrt(2) * sp.atan(rho)
    g = rho * R / (2 if MUT == "d_gauge" else sp.sqrt(2))
    inv = sp.simplify(sp.diff(R, rho) - g * sp.diff(Phi, rho))
    ok1 = inv == 0
    check("D3a", "Y = dr - (rho R/sqrt2) dphi is invariant under every linear diffeomorphism: L_xi r = xi^mu d_mu(rho) R', L_xi phi = xi^mu d_mu(rho) Phi', and R' - (rho R/sqrt2) Phi' = 0 identically (sympy)", ok1, f"R' - g Phi' = {inv}")
    # the horizon: r_U = R'(rho) rho_U = -R'(rho_h) G V on U = 0
    Rp = sp.diff(R, rho)
    c15 = sp.Rational(-3, 2)
    _, A15, _, _, _, _ = _member_sym(c15)
    rh = sp.nsolve(A15, rho, 0.34)
    ok2 = sp.simplify(Rp.subs(rho, 0)) == 0 and abs(float(Rp.subs(rho, rh))) > 0.3
    check("D3b", "r_U on the future horizon is -R'(rho_h) G(rho_h) V: zero at c = -1 (R'(0) = 0, the horizon is the area minimum, so no residual gauge shifts dr_h) and nonzero at c = -1.5 (R'(rho_h) = %.4f), where dr_h = c2 V is the horizon-position gauge" % float(Rp.subs(rho, rh)), ok2)


def derive_D6():
    """Formal series of Y near S: the equation b^2 Pi (rho^2 Y_VU - rho (rho_U Y_V + rho_V Y_U) + Y (2 rho_U rho_V - rho rho_VU)) + Vhat rho^2 Y/2 = 0 (u_pm = -Vhat u/2, u = Y/rho, d_p d_m = b^2 Pi d_V d_U),
    Y = sum y_mn V^m U^n with rho(Pi) and Vhat(rho) series at c = -1: the free coefficients are exactly y_mn with m = 1 or n = 1; y_00 and all y_0n, y_m0 (n, m != 1) vanish."""
    import sympy as sp
    rho, Pi = sp.symbols("rho Pi")
    Vs, Us = sp.symbols("V U")
    N = 7
    c = -1
    K = -2 * sp.Integer(c) / sp.pi
    A = 1 + c * (rho ** 2 + 1) + K * (rho + (rho ** 2 + 1) * sp.atan(rho))
    As = sp.series(A, rho, 0, N + 2).removeO()
    A1 = As.coeff(rho, 1)
    kap = A1 / 2
    R2 = rho ** 2 + 1
    Ap = sp.diff(As, rho)
    Vh = sp.series(sp.simplify(2 * As ** 2 / rho ** 2 - (1 + 1 / R2) * As * Ap / rho + As ** 2 / R2 ** 2 - 2 * As / R2), rho, 0, N).removeO()
    integrand = sp.series(2 * kap / As - 1 / rho, rho, 0, N).removeO()
    Pi_of_rho = sp.series(rho * sp.exp(sp.integrate(integrand, rho)), rho, 0, N).removeO()           # normalization: Pi = rho + ...
    rs = sp.symbols("r1:%d" % N)
    rho_series = sum(rs[k] * Pi ** (k + 1) for k in range(N - 1))
    eq = sp.expand(sp.series(Pi_of_rho.subs(rho, rho_series) - Pi, Pi, 0, N).removeO())
    sol = {}
    for k in range(1, N):
        sol[rs[k - 1]] = sp.solve(eq.coeff(Pi, k).subs(sol), rs[k - 1])[0]
    rho_Pi = sp.expand(rho_series.subs(sol))
    D = 7
    b2 = 2 * kap ** 2 * (1.1 if MUT == "d_res" else 1)

    def trunc(expr, deg):
        p = sp.Poly(sp.expand(expr), Vs, Us)
        return sum(cf * Vs ** m * Us ** n for (m, n), cf in p.terms() if m + n <= deg)
    rho_VU = trunc(rho_Pi.subs(Pi, -Vs * Us), D + 6)
    Vhat_VU = trunc(Vh.subs(rho, rho_VU), D + 6)
    rV, rU, rVU = sp.diff(rho_VU, Vs), sp.diff(rho_VU, Us), sp.diff(rho_VU, Vs, Us)
    ys = {(m, n): sp.Symbol(f"y_{m}_{n}") for m in range(D + 1) for n in range(D + 1 - m)}
    Y = sum(ys[k] * Vs ** k[0] * Us ** k[1] for k in ys)
    expr = b2 * (-Vs * Us) * (rho_VU ** 2 * sp.diff(Y, Vs, Us) - rho_VU * (rU * sp.diff(Y, Vs) + rV * sp.diff(Y, Us)) + Y * (2 * rU * rV - rho_VU * rVU)) + Vhat_VU * rho_VU ** 2 * Y / 2
    poly = sp.Poly(sp.expand(expr), Vs, Us)
    eqs = [cf for (m, n), cf in poly.terms() if m + n <= D + 4]            # unknowns of degree <= D enter at monomial degree <= D + 4 only through these
    res = [k for k in ys if k[0] == 1 or k[1] == 1]
    nonres = [k for k in ys if k not in res]
    unk_nr = [ys[k] for k in nonres]
    solved = sp.solve(eqs, unk_nr, dict=True)
    ok = len(solved) == 1 and all(u in solved[0] for u in unk_nr)
    zero_set, bad = [], []
    if ok:
        s0 = solved[0]
        # (i) consistency: with the resonant coefficients as free parameters no equation is left over
        left = [sp.simplify(e.subs(s0)) for e in eqs]
        ok = all(e == 0 for e in left)
        # (ii) the non-resonant coefficients that vanish identically: y_00, y_0n (n != 1), y_m0 (m != 1)
        want_zero = [k for k in nonres if k[0] == 0 or k[1] == 0]
        ok = ok and all(sp.simplify(s0[ys[k]]) == 0 for k in want_zero)
        zero_set = want_zero
    check("D6", "formal series of the Y equation at the bifurcation surface (c = -1, M = 0, degree <= 7): with the %d coefficients y_mn, m = 1 or n = 1, as free parameters every other coefficient is determined and the system is consistent; y_00 and every y_0n, y_m0 (n, m != 1) vanish, so Y|_(U=0) = y_10 V and Y|_(V=0) = y_01 U exactly (sympy)" % len(res),
          ok, f"{len(eqs)} equations, {len(unk_nr)} determined unknowns, {len(zero_set)} forced to zero")


def derive_D7():
    """The growing solution at W = kappa: Frobenius series of #212's Y equation at c = -1 with exponent 1/2 has no logarithm; Y_ref = rho^{1/2} sum g_k rho^k, so Y_ref e^{-kappa x} is analytic in rho (Pi)."""
    import sympy as sp
    rho = sp.Symbol("rho")
    c = -1
    K = -2 * sp.Integer(c) / sp.pi
    N = 7
    A = sp.series(1 + c * (rho ** 2 + 1) + K * (rho + (rho ** 2 + 1) * sp.atan(rho)), rho, 0, N + 4).removeO()
    kap = A.coeff(rho, 1) / 2
    W = kap * (sp.Rational(11, 10) if MUT == "d_log" else 1)
    R2 = rho ** 2 + 1
    Ap = sp.diff(A, rho)
    gs = sp.symbols("g0:%d" % (N + 1))
    s = sp.Rational(1, 2)
    Yx = sum(gs[k] * rho ** (s + k) for k in range(N + 1))
    resid = sp.diff(Yx, rho, 2) - (2 / rho - Ap / A) * sp.diff(Yx, rho) - (1 / R2 ** 2 - Ap / (A * rho * R2) - 2 / (A * R2) + W ** 2 / A ** 2) * Yx
    ser = sp.expand(sp.series(sp.expand(resid * rho ** (sp.Rational(3, 2)) * 1), rho, 0, N - 1).removeO())      # leading power rho^{-3/2} * rho^{3/2}
    sol = {gs[0]: 1, gs[1]: 0}                               # the resonance at k = 1 (exponent 3/2): g1 is the free multiple of the other solution
    consistent = True
    for k in range(N - 1):
        eqk = sp.simplify(ser.coeff(rho, k).subs(sol))
        unknown = [g for g in gs if g in eqk.free_symbols and g not in sol]
        if not unknown:
            consistent = consistent and sp.simplify(eqk) == 0
            continue
        sl = sp.solve(eqk, unknown[-1])
        if not sl:
            consistent = False
            break
        sol[unknown[-1]] = sl[0]
    check("D7", "c = -1, W = kappa: the Frobenius series of Y with exponent 1/2 has no logarithm (every order solvable, the order-1 equation consistent), so Y_ref = rho^{1/2}(1 + g2 rho^2 + ...) and Y_ref e^{-kappa x}, e^{kappa x} = rho^{1/2}(1 + O(rho)), is analytic in rho and in Pi = -VU (sympy, orders <= %d)" % (N - 2),
          consistent, f"g2 = {sp.nsimplify(sol.get(gs[2], 'n/a'))}" if consistent else "inconsistent at some order")


def derive_D8():
    """Y_ref e^{-kappa x}: the ODE has no singular point on rho < 0 (A < 0, R >= 1): integrate the W = kappa solution from the series into the interior and the static region."""
    from scipy.integrate import solve_ivp
    mem = lin.ep.Member(-1.0)
    kap = mem.kappa
    K = mem.K
    A = lambda r: 1 - (r * r + 1) + K * (r + (r * r + 1) * math.atan(r))
    Ap = lambda r: 2 * K * (1 + r * math.atan(r)) - 2 * r
    # A < 0 on a dense grid of rho < 0 (no second horizon), A > 0 for rho > 0
    gr = np.linspace(-200.0, -1e-3, 200000)
    neg_ok = bool(np.all(np.array([A(x) for x in gr[::50]]) < 0)) and all(A(x) > 0 for x in np.linspace(1e-3, 200.0, 4000))
    # g0 = 1, g1 = 0 and g2 from the ODE at W = kappa (from derive_D7: use the series value computed numerically by Richardson-free shooting from rho = 0.05)

    def f(r, y, W=kap):
        P = 2 / r - Ap(r) / A(r)
        Q = 1 / (r * r + 1) ** 2 - Ap(r) / (A(r) * r * (r * r + 1)) - 2 / (A(r) * (r * r + 1)) + W * W / A(r) ** 2
        return [y[1], P * y[1] + Q * y[0]]
    out = {}
    for sgn in (1.0, -1.0):
        r0 = sgn * 0.05
        # start vector from rho^{1/2}(1 + g2 rho^2): g2 unknown to first order; shoot with g2 = 0 and assess smoothness (g2 rho^2 ~ 1e-3 relative): the exponent-3/2 admixture is irrelevant for finiteness
        y0 = [abs(r0) ** 0.5, 0.5 * sgn * abs(r0) ** -0.5]
        sol = solve_ivp(f, [r0, sgn * 40.0], y0, method="DOP853", rtol=1e-11, atol=1e-14)
        out[sgn] = (sol.status, float(np.max(np.abs(sol.y[0]))), sol.t[-1])
    ok = neg_ok and out[1.0][0] == 0 and out[-1.0][0] == 0 and np.isfinite(out[1.0][1]) and np.isfinite(out[-1.0][1])
    check("D8", "the Y equation at c = -1, W = kappa has no singular point on rho in [-200, 0) or (0, 200] (A < 0 for rho < 0, A > 0 for rho > 0, R >= 1): the solution with Y ~ rho^{1/2} integrates from rho = +-0.05 to +-40 with finite values (the interior is reached by the same analytic function)",
          ok, f"A<0 on rho<0 and A>0 on rho>0: {neg_ok}; integration status {out[1.0][0]}, {out[-1.0][0]}; max |Y| {out[1.0][1]:.3g} (static), {out[-1.0][1]:.3g} (interior)")


# ----------------------------------------------------------------------------------------------------------------------
# mutations: each breaks one thing; the checks it must fire were computed in development runs (negative-mass-raw I248-dh610/dev) before they were declared
# ----------------------------------------------------------------------------------------------------------------------
MUTATIONS = {
    # name: checks that must fire (development runs dev07 .. dev09 and dev04/05 of negative-mass-raw I248-dh610, before they were declared)
    "d_emF": ["D1"],                 # e^{-F} = 1.01 G/kappa in the symbolic Kruskal background
    "d_Gode": ["D1"],                # d ln G/d rho scaled by 1.01 (the ODE of G)
    "d_hz": ["D2"],                  # the constraint C_VV evaluated at U = 0.1, off the horizon
    "d_prop": ["D4"],                # the sign of the eps r phi_V e_V term of the constraints
    "d_gauge": ["D3a"],              # g = rho R/2 instead of rho R/sqrt2 in Y
    "d_res": ["D6"],                 # b^2 scaled by 1.1 in the formal series (the resonances move)
    "d_log": ["D7"],                 # W = 1.1 kappa in the Frobenius series (the logarithm appears)
    "k_emF": ["N1", "N2", "N3", "N4", "N6", "N7", "N9", "N10", "N11"],     # e^{-F} 1% off in the Kruskal-chart solver
    "k_FV": ["N1", "N2", "N3", "N4", "N6", "N7", "N9", "N10", "N11"],      # sign of F_V, F_U in the Kruskal-chart background
    "k_cross": ["N2", "N3", "N4", "N6", "N7", "N9", "N10", "N11"],         # sign of the (r_V/r) dphi_U + (r_U/r) dphi_V terms of the solver
    "k_con": ["N3", "N4", "N7", "N11"],                                     # the ray constraint with the wrong sign: the data violate C_VV, C_UU
    "fn_aT": ["N4", "N5", "N6", "N7", "N11"],                               # the functional uses the tangential half d_T only
    "fn_dirU": ["N4", "N5", "N6", "N7", "N11"],                             # the functional takes d_U instead of d_V
    "fn_b": ["N5"],                                                         # b = sqrt2 kappa 1% off in the conversion e^{bp} -> V of the (p, sigma) solutions
    "plant_sign": ["N6"],                                                   # the sign of the planted known coefficient
    "ctl_lambda": ["N7"],                                                   # lambda of the tuned superposition 10% off
}


DERIVE = [("D1", derive_D1), ("D2", derive_D2), ("D3", derive_D3), ("D6", derive_D6), ("D7", derive_D7), ("D8", derive_D8)]
NUMERIC = [("N1", check_N1), ("N2", check_N2), ("N3", check_N3), ("N4", check_N4), ("N5", check_N5), ("N6", check_N6), ("N7", check_N7), ("N8", check_N8), ("N9", check_N9), ("N10", check_N10), ("N11", check_N11)]


def run_all(items):
    t0 = time.time()
    for cid, fn in items:
        if want(cid) or any(want(x) for x in (cid + "a", cid + "b", "D4" if cid == "D2" else cid)):
            t1 = time.time()
            fn()
            print(f"    ({cid}: {time.time() - t1:.0f} s)", flush=True)
    bad = [cid for cid, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS)} checks, {len(RESULTS) - len(bad)} pass; failed: {bad}; mutation {MUT}; {time.time() - t0:.0f} s")
    if MUT is not None and MUT in MUTATIONS:
        want_ = MUTATIONS[MUT]
        print(f"mutation {MUT}: fired {bad}, declared {want_}: {'equal' if sorted(bad) == sorted(want_) else 'DIFFERENT'}")
    return bad


def table():
    """The numbers of the note: c1, c2 of the classes, the slice values, the interior, and the exterior comparison."""
    print("c = -1 black universe, Kruskal chart, rectangle V in [-0.05, 3], U in [-2, 0.06]; Richardson of h = 0.0025, 0.00125 (c2 of the horizon fit; of the slices crossing at V = 0, 0.5, 1)")
    print(f"{'class':5s} {'c1':>13s} {'c2 (fit)':>14s} {'c2 slice V=0':>14s} {'V=0.5':>14s} {'V=1':>14s}   fit residual h=0.00125")
    for name in ("L", "G", "A0", "A1R", "A1L", "R"):
        runs = [run_class(-1.0, name, h) for h in HS]
        rf = richardson2(runs[2]["c2"], runs[3]["c2"])
        r1 = richardson2(runs[2]["c1"], runs[3]["c1"])
        sl = [richardson2(runs[2]["slice"][v], runs[3]["slice"][v]) for v in SLICES]
        print(f"{name:5s} {r1:13.6e} {rf:14.6e} {sl[0]:14.6e} {sl[1]:14.6e} {sl[2]:14.6e}   {runs[3]['res']:.1e}")
    print("\nsamples of the converged fields (class G, h = 0.00125): dr, dF, dphi, Y")
    run = run_class(-1.0, "G", HS[3])
    for k, (Vp, Up) in SAMPLES.items():
        print(f"  {k:10s} (V, U) = ({Vp:6.2f}, {Up:6.2f})  " + "  ".join(f"{v:12.5e}" for v in run["samples"][k]))


def report_mutations(raw, prefix):
    """Fired sets of stored mutation logs (the 'failed: [...]' line of each) against the declared ones."""
    import ast
    import re
    print("== Mutations: fired sets in %s<name>.log against the declared ones ==" % prefix)
    n_eq = 0
    for name, declared in MUTATIONS.items():
        fn = os.path.join(raw, f"{prefix}{name}.log")
        if not os.path.exists(fn):
            print(f"  {name:11s} NO LOG")
            continue
        m = re.findall(r"failed: (\[.*?\]); mutation", open(fn).read())
        fired = ast.literal_eval(m[-1]) if m else None
        eq = fired is not None and sorted(fired) == sorted(declared)
        n_eq += eq
        print(f"  {name:11s} fired {fired}  declared {declared}  {'equal' if eq else 'DIFFERENT'}")
    print(f"{n_eq} of {len(MUTATIONS)} equal")
    return n_eq == len(MUTATIONS)


def main():
    sub = next((a for a in sys.argv[1:] if not a.startswith("--")), "check")
    if sub == "derive":
        bad = run_all(DERIVE)
    elif sub == "check":
        bad = run_all(NUMERIC)
    elif sub == "table":
        table()
        return
    elif sub == "mutations":
        args = [a for a in sys.argv[2:] if not a.startswith("--")]
        sys.exit(0 if report_mutations(args[0], args[1] if len(args) > 1 else "") else 1)
    else:
        sys.exit(f"unknown command {sub}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

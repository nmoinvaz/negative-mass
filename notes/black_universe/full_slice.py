#!/usr/bin/env python3
"""A compactly supported datum on a full slice of the c = -1 black universe, constraints solved, evolved (#258).

The slice T = 0 is V + U = 0 of the Kruskal chart of #248 (V = s, U = -s), through the bifurcation surface S (s = 0) from the left to the right
exterior. Free data: dphi(s) and dphi_N(s) = d_N dphi, C^inf bumps. Gauge: dF = dF_N = 0 on the slice. The two linearized constraints are reduced to
two ODEs for dr(s) and n(s) = d_N dr (derivation in full_slice.md, checked here on exact pure-gauge solutions). The datum is evolved with the cell scheme of
horizon_sector.py (linear.py) from two diagonal layers, and c2 = d_V dr(S) read from the slice is compared with the fit of the evolved horizon values.

    full_slice.py check [--only=Q1,..] [--mutate=NAME]    checks; exit 1 when one fails
    full_slice.py table                                   the data tables (c1, c2, Y at S, parity, the first integral)
    full_slice.py mutations DIR PREFIX                    fired sets of the stored logs DIR/PREFIX<name>.log against the declared ones

Mutations: H_P H_eps M_phi M_F fn_aT fn_sign l1_sign l1_noquad dm_noT mass_off odd_off asym_S left_mass kg_eps.
Reads horizon_sector.py and linear.py (read-only imports). Python: ~/Source/negative-mass/.venv/bin/python.
"""
import math
import os
import sys
import time

import numba as nb
import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
_ARGV = sys.argv
sys.argv = sys.argv[:1]                     # the imported modules parse --mutate / --only themselves
import horizon_sector as hs  # noqa: E402
import linear  # noqa: E402,F401
sys.argv = _ARGV

MUT = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mutate=")), None)
ONLY = next((a.split("=", 1)[1].split(",") for a in sys.argv[1:] if a.startswith("--only=")), None)
EPS = -1.0
SL, SR, SFAR = -1.05, 1.05, 3.0e4            # integration interval of the constraint ODEs, far end of the tails
HLIST = (0.01, 0.005, 0.0025, 0.00125)
VMAX = 3.0                                    # the right-half march reaches V = 3
JD = 0.5                                      # the diamond over |s| <= 0.5
KAPPA = 2.0 / math.pi
KNOWN = ("H_P", "H_eps", "M_phi", "M_F", "fn_aT", "fn_sign", "l1_sign", "l1_noquad", "dm_noT", "mass_off", "odd_off", "asym_S", "left_mass", "kg_eps")
assert MUT is None or MUT in KNOWN, MUT


# ----------------------------------------------------------------------------------------------------------------------
# background along the slice
# ----------------------------------------------------------------------------------------------------------------------
class Bg:
    """The member c (or Schwarzschild, r_s = 1) in the Kruskal chart of #248, evaluated on the slice V = s, U = -s (functions of Pi = s^2)."""

    def __init__(self, c=-1.0, schw=False, pi_max=1e12):
        self.kc = hs.KChart(c, pi_max=pi_max)
        self.par = self.kc.par.copy()
        if schw:
            self.par[hs.P_KIND2] = 1.0
        self.tab = self.kc.tab
        self.kappa = self.kc.kappa

    def at(self, s):
        a = max(abs(s), 1e-8)
        x, d, r, rV, rU, pV, pU, V0, V1, V2, E, FV, FU, rhoV, rhoU = hs.jetsK(a, -a, self.par, self.tab, 0.0)
        return dict(rho=x, r=r, rp=rV / a, php=pV / a, Fp=FV / a, E=E, V0=V0, V1=V1)


_BG = {}


def bgc(c=-1.0, schw=False):
    key = (c, schw)
    if key not in _BG:
        _BG[key] = Bg(c, schw, pi_max=1e12 if c == -1.0 and not schw else 100.0)
    return _BG[key]


def bump(x):
    return math.e * math.exp(-1.0 / (1.0 - x * x)) if abs(x) < 1.0 else 0.0


def dbump(x):
    return bump(x) * (-2.0 * x / (1.0 - x * x) ** 2) if abs(x) < 1.0 else 0.0


def ddbump(x):
    if abs(x) >= 1.0:
        return 0.0
    u = 1.0 - x * x
    return bump(x) * ((2.0 * x / u ** 2) ** 2 - 2.0 / u ** 2 - 8.0 * x * x / u ** 3)


class Data:
    """A sum of amp * bump((s - c)/w); parts = [(amp, c, w), ...]."""

    def __init__(self, parts):
        self.parts = list(parts)

    def f(self, s):
        return sum(a * bump((s - c) / w) for a, c, w in self.parts)

    def fT(self, s):
        return sum(a * dbump((s - c) / w) / w for a, c, w in self.parts)

    def fTT(self, s):
        return sum(a * ddbump((s - c) / w) / w ** 2 for a, c, w in self.parts)

    def __add__(self, o):
        return Data(self.parts + o.parts)

    def scaled(self, k):
        return Data([(a * k, c, w) for a, c, w in self.parts])

    def amp(self):
        return max(abs(self.f(s)) for s in np.linspace(-1, 1, 801))


# ----------------------------------------------------------------------------------------------------------------------
# the constraint ODEs on the slice (gauge dF = dF_N = 0 in production; dF sources for the planted exact solutions)
# ----------------------------------------------------------------------------------------------------------------------
class Slice:
    """H: dr'' + P dr' + Q dr = -[E V' r dphi + eps r s phi' dphi'] - [E (1/r - V r) dF + s r' dF'];   M: n' + s F' n = -s r' dF_N - eps r s phi' dphi_N."""

    def __init__(self, bg, mut=None):
        self.bg = bg
        self.mut = mut

    def PQ(self, s):
        b = self.bg.at(s)
        m = self.mut
        P = s * b["Fp"] + (0.0 if m == "H_P" else 2 * s * b["rp"] / b["r"])
        Q = b["E"] * b["V0"] + b["E"] / b["r"] ** 2 - 2 * s * s * b["rp"] ** 2 / b["r"] ** 2 + (-1.0 if m == "H_eps" else 1.0) * EPS * s * s * b["php"] ** 2
        return P, Q, b

    def rhs_even(self, s, y, dphi, dphiT, dF=0.0, dFT=0.0):
        P, Q, b = self.PQ(s)
        S = -(b["E"] * b["V1"] * b["r"] * dphi + EPS * b["r"] * s * b["php"] * dphiT) - (b["E"] * (1 / b["r"] - b["V0"] * b["r"]) * dF + s * b["rp"] * dFT)
        return [y[1], S - P * y[1] - Q * y[0]]

    def rhs_odd(self, s, n, dphiN, dFN=0.0):
        b = self.bg.at(s)
        m = self.mut
        return -(-1.0 if m == "M_F" else 1.0) * s * b["Fp"] * n - s * b["rp"] * dFN - (-1.0 if m == "M_phi" else 1.0) * EPS * b["r"] * s * b["php"] * dphiN

    def even(self, a, b_, D, y0, **kw):
        return solve_ivp(lambda s, y: self.rhs_even(s, y, D.f(s), D.fT(s)), [a, b_], list(y0), dense_output=True, method="DOP853", rtol=1e-12, atol=1e-18, **kw)

    def odd(self, a, b_, D, n0, **kw):
        return solve_ivp(lambda s, y: [self.rhs_odd(s, y[0], D.f(s))], [a, b_], [n0], dense_output=True, method="DOP853", rtol=1e-12, atol=1e-18, **kw)

    def hom_even(self, a, b_, y0):
        zero = Data([])
        return self.even(a, b_, zero, y0)

    def hom_odd(self, a, b_, n0):
        return self.odd(a, b_, Data([]), n0)

    # first integral Theta = a dr + mu dr' + c dphi,  mu = s r,  a = mu P - mu',  c = eps s^2 r^2 phi'  (Theta' = mu H)
    def theta_coef(self, s):
        P, Q, b = self.PQ(s)
        mu = s * b["r"]
        a = mu * P - (b["r"] + 2 * s * s * b["rp"])
        c = (-1.0 if self.mut == "kg_eps" else 1.0) * EPS * s * s * b["r"] ** 2 * b["php"]
        return mu, a, c, Q, b

    def theta(self, s, dr, drs, dphi=0.0):
        mu, a, c, Q, b = self.theta_coef(s)
        return a * dr + mu * drs + c * dphi


def dm_of(bg, s, dr, drs, dF=0.0):
    """Misner-Sharp mass perturbation on the slice: dm = dr (1/2 - e^F s^2 r'^2) - r e^F s r' (s r' dF + dr')."""
    b = bg.at(s)
    eF = 1.0 / b["E"]
    grp = 0.0 if MUT == "dm_noT" else 1.0
    return dr * (0.5 - eF * s * s * b["rp"] ** 2) - grp * b["r"] * eF * s * b["rp"] * (s * b["rp"] * dF + drs)


def c_functionals(drs, n):
    """c2 = d_V dr(S) and c2' = d_U dr(S) from the slice data at S."""
    if MUT == "fn_aT":
        return 0.5 * drs, -0.5 * drs
    if MUT == "fn_sign":
        return 0.5 * (drs - n), 0.5 * (n + drs)
    return 0.5 * (drs + n), 0.5 * (n - drs)


# ----------------------------------------------------------------------------------------------------------------------
# constructing a datum: corrections, end conditions, and the mutated asymptotic conditions
# ----------------------------------------------------------------------------------------------------------------------
class Datum:
    """Seed plus two correction bumps; the solution of H and M with the left start and, as the cross-check, with the right start."""

    def __init__(self, sl, seed_phi, seed_phiN, corr_phi, corr_phiN, mode=None, name=""):
        self.sl, self.name, self.mode = sl, name, mode
        self.bg = sl.bg
        mass, odd = mode != "mass_off", mode != "odd_off"
        A0 = sl.even(SL, SR, seed_phi, (0, 0), max_step=0.01).y[0, -1]
        Ac = sl.even(SL, SR, corr_phi, (0, 0), max_step=0.01).y[0, -1] if corr_phi.parts else 1.0
        N0 = sl.odd(SL, SR, seed_phiN, 0.0, max_step=0.01).y[0, -1]
        Nc = sl.odd(SL, SR, corr_phiN, 0.0, max_step=0.01).y[0, -1] if corr_phiN.parts else 1.0
        self.km = -A0 / Ac if (mass and corr_phi.parts) else 0.0
        self.ko = -N0 / Nc if (odd and corr_phiN.parts) else 0.0
        self.dphi = seed_phi + corr_phi.scaled(self.km)
        self.dphiN = seed_phiN + corr_phiN.scaled(self.ko)
        D, DN = self.dphi, self.dphiN
        if mode == "asym_S":
            self.pe = (sl.even(0.0, SR, D, (0, 0), max_step=0.01), sl.even(0.0, SL, D, (0, 0), max_step=0.01))
            self.po = (sl.odd(0.0, SR, DN, 0.0, max_step=0.01), sl.odd(0.0, SL, DN, 0.0, max_step=0.01))
            self.se, self.so = self.pe, self.po
            self.ev = lambda s: self._pw(self.pe, self.po, s)
            self.ev2 = self.ev
        else:
            y0 = (0.0, 0.0)
            if mode == "left_mass":
                y0 = tuple(1e-4 * self.mass_start(sl))
            pe = sl.even(SL, SR, D, y0, max_step=0.01)
            po = sl.odd(SL, SR, DN, 0.0, max_step=0.01)
            se = sl.even(SR, SL, D, (0, 0), max_step=0.01)
            so = sl.odd(SR, SL, DN, 0.0, max_step=0.01)
            self.pe, self.po, self.se, self.so = pe, po, se, so
            self.ev = lambda s: self._one(pe, po, s)
            self.ev2 = lambda s: self._one(se, so, s)
        self.c1, drs, n0 = self.ev(0.0)[0], self.ev(0.0)[1], self.ev(0.0)[2]
        self.c2, self.c2p = c_functionals(drs, n0)
        r2 = self.ev2(0.0)
        self.c2_2, self.c2p_2 = c_functionals(r2[1], r2[2])
        self.c1_2 = r2[0]

    @staticmethod
    def _one(sE, sO, s):
        if abs(s) > SR + 1e-12:
            return np.zeros(3)
        y = sE.sol(s)
        return np.array([y[0], y[1], sO.sol(s)[0]])

    @staticmethod
    def _pw(pe, po, s):
        if abs(s) > SR + 1e-12:
            return np.zeros(3)
        i = 0 if s >= 0 else 1
        y = pe[i].sol(s)
        return np.array([y[0], y[1], po[i].sol(s)[0]])

    @staticmethod
    def mass_start(sl):
        """start vector (dr, dr') at SL of the homogeneous solution with unit delta m and no dr tail at SFAR."""
        h = [sl.hom_even(SL, SFAR, y0) for y0 in ((1.0, 0.0), (0.0, 1.0))]
        M = np.array([[dm_of(sl.bg, SFAR, *h[k].sol(SFAR)), h[k].sol(SFAR)[0]] for k in (0, 1)]).T
        return np.linalg.solve(M, np.array([1.0, 0.0]))

    def slice_fn(self):
        """s -> (U, U_T, U_N, U_TT) for U = (dr, dF, dphi), the data of the layers."""
        sl, D, DN = self.sl, self.dphi, self.dphiN

        def fn(s):
            if abs(s) > SR:
                z = np.zeros(3)
                return z, z, z, z
            dr, drT, n = self.ev(s)
            P, Q, b = sl.PQ(s)
            S = -(b["E"] * b["V1"] * b["r"] * D.f(s) + EPS * b["r"] * s * b["php"] * D.fT(s))
            drTT = S - P * drT - Q * dr
            return (np.array([dr, 0.0, D.f(s)]), np.array([drT, 0.0, D.fT(s)]), np.array([n, 0.0, DN.f(s)]), np.array([drTT, 0.0, D.fTT(s)]))
        return fn

    def Y(self, s):
        """gauge-invariant Y = dr - g dphi, g = rho R/sqrt2 (rho >= 0 on both exteriors)."""
        dr = self.ev(s)[0]
        b = self.bg.at(s)
        return dr - b["rho"] * b["r"] / math.sqrt(2.0) * self.dphi.f(s)


def tails(sl, end, direction, hom_s=None):
    """homogeneous continuation of (dr, dr', n) from SR (direction +1) or SL (-1) to +-SFAR; returns [(s, dr, dr', n)] at the sample points."""
    s0 = SR if direction > 0 else SL
    far = SFAR * direction
    e = sl.hom_even(s0, far, (end[0], end[1]))
    o = sl.hom_odd(s0, far, end[2])
    out = []
    for s in (s0, 10.0 * direction, 100.0 * direction, 1e3 * direction, 1e4 * direction, far):
        y = e.sol(s)
        out.append((s, y[0], y[1], o.sol(s)[0]))
    return out


# ----------------------------------------------------------------------------------------------------------------------
# the evolution: two diagonal layers, the cell scheme of horizon_sector.py
# ----------------------------------------------------------------------------------------------------------------------
@nb.njit(cache=True)
def march_half(L0, L1, J, h, par, tab):
    """Right half-slice, nodes (j, i) = (V/h, U/h), i <= 0. L0[q, j] at (j, -j), j = 0..J; L1[q, j] at (j, 1 - j), j = 1..J. Ua[q, j, i + J]."""
    Ua = np.zeros((3, J + 1, J + 2))
    for q in range(3):
        for j in range(J + 1):
            Ua[q, j, -j + J] = L0[q, j]
        for j in range(1, J + 1):
            Ua[q, j, 1 - j + J] = L1[q, j]
    B0 = np.zeros((3, 3))
    Bp = np.zeros((3, 3))
    Bm = np.zeros((3, 3))
    M = np.zeros((3, 3))
    rhs = np.zeros(3)
    out = np.zeros(3)
    Sv = np.zeros(3)
    Ev = np.zeros(3)
    Wv = np.zeros(3)
    for ll in range(2, J + 1):
        guess = 0.0
        for j in range(ll, J + 1):
            i = ll - j
            for q in range(3):
                Sv[q] = Ua[q, j - 1, i - 1 + J]
                Ev[q] = Ua[q, j, i - 1 + J]
                Wv[q] = Ua[q, j - 1, i + J]
            guess = hs._cellK(Sv, Ev, Wv, (j - 1) * h, (i - 1) * h, h, h, par, tab, guess, B0, Bp, Bm, M, rhs, out)
            for q in range(3):
                Ua[q, j, i + J] = out[q]
    return Ua


@nb.njit(cache=True)
def march_diamond(L0, L1, J, h, par, tab):
    """The domain of dependence of the slice nodes j = -J..J (V <= Jh, U <= Jh, V + U >= 0). Ua[q, j + J, i + J]."""
    n = 2 * J + 1
    Ua = np.zeros((3, n, n))
    for q in range(3):
        for j in range(-J, J + 1):
            Ua[q, j + J, -j + J] = L0[q, j + J]
        for j in range(1 - J, J + 1):
            Ua[q, j + J, 1 - j + J] = L1[q, j + J]
    B0 = np.zeros((3, 3))
    Bp = np.zeros((3, 3))
    Bm = np.zeros((3, 3))
    M = np.zeros((3, 3))
    rhs = np.zeros(3)
    out = np.zeros(3)
    Sv = np.zeros(3)
    Ev = np.zeros(3)
    Wv = np.zeros(3)
    for ll in range(2, 2 * J + 1):
        guess = 0.0
        for j in range(ll - J, J + 1):
            i = ll - j
            for q in range(3):
                Sv[q] = Ua[q, j - 1 + J, i - 1 + J]
                Ev[q] = Ua[q, j + J, i - 1 + J]
                Wv[q] = Ua[q, j - 1 + J, i + J]
            guess = hs._cellK(Sv, Ev, Wv, (j - 1) * h, (i - 1) * h, h, h, par, tab, guess, B0, Bp, Bm, M, rhs, out)
            for q in range(3):
                Ua[q, j + J, i + J] = out[q]
    return Ua


def layer_pair(bg, h, js, fn):
    """L0 at s = j h (j in js), L1 at the points (h/2)(1, 1) off the slice midpoints s = (j - 1/2) h (j in js[1:] when js starts at the first node)."""
    B0 = np.zeros((3, 3))
    Bp = np.zeros((3, 3))
    Bm = np.zeros((3, 3))
    L0 = np.zeros((3, len(js)))
    L1 = np.zeros((3, len(js)))
    sgn = -1.0 if MUT == "l1_sign" else 1.0
    quad = 0.0 if MUT == "l1_noquad" else 1.0
    for k, j in enumerate(js):
        L0[:, k] = fn(j * h)[0]
    for k, j in enumerate(js):
        if k == 0:
            continue
        s = (j - 0.5) * h
        U, UT, UN, UTT = fn(s)
        hs.coefK(s, -s, bg.par, bg.tab, 0.0, B0, Bp, Bm)
        UV = 0.5 * (UT + UN)
        UU = 0.5 * (UN - UT)
        UNN = UTT + 4.0 * (B0 @ U + Bp @ UV + Bm @ UU)
        L1[:, k] = U + sgn * 0.5 * h * UN + quad * h * h / 8.0 * UNN
    return L0, L1


def run_half(bg, fn, h, vmax=VMAX):
    J = int(round(vmax / h))
    L0, L1 = layer_pair(bg, h, list(range(0, J + 1)), fn)
    return march_half(L0, L1, J, h, bg.par, bg.tab), J


def run_diamond(bg, fn, h, half=JD):
    J = int(round(half / h))
    L0, L1 = layer_pair(bg, h, list(range(-J, J + 1)), fn)
    return march_diamond(L0, L1, J, h, bg.par, bg.tab), J


def fit_half(Ua, J, h):
    V = h * np.arange(J + 1)
    return hs.horizon_fit(Ua[0, :, J], V)


def richardson(vals, hl):
    return [(4.0 * vals[hl[k + 1]] - vals[hl[k]]) / 3.0 for k in range(len(hl) - 1)]


CPTS = [(0.6, -0.2), (1.0, -0.5), (1.5, -0.3), (2.0, -1.2), (0.8, -0.4), (0.4, -0.1), (2.5, -0.2), (1.2, -1.0)]


def constraints_at(bg, Ua, J, h, pts=CPTS):
    """C_VV, C_UU of the evolved fields (lemma of #261: dr_VV + F_V dr_V + r_V dF_V + (eps/2) phi_V^2 dr + eps r phi_V dphi_V) by 4th-order differences;
    returns max |C| and the largest term at the points."""
    mc, mt = 0.0, 0.0
    for V, U in pts:
        j, i = int(round(V / h)), int(round(U / h))
        x, d, r, rV, rU, pV, pU, V0_, V1_, V2_, E, FV, FU, rhoV, rhoU = hs.jetsK(j * h, i * h, bg.par, bg.tab, 0.0)

        def f(q, dj, di):
            return Ua[q, j + dj, i + di + J]

        d1V = lambda q: (-f(q, 2, 0) + 8 * f(q, 1, 0) - 8 * f(q, -1, 0) + f(q, -2, 0)) / (12 * h)
        d1U = lambda q: (-f(q, 0, 2) + 8 * f(q, 0, 1) - 8 * f(q, 0, -1) + f(q, 0, -2)) / (12 * h)
        d2V = lambda q: (-f(q, 2, 0) + 16 * f(q, 1, 0) - 30 * f(q, 0, 0) + 16 * f(q, -1, 0) - f(q, -2, 0)) / (12 * h * h)
        d2U = lambda q: (-f(q, 0, 2) + 16 * f(q, 0, 1) - 30 * f(q, 0, 0) + 16 * f(q, 0, -1) - f(q, 0, -2)) / (12 * h * h)
        for d2, d1, Fx, px, rx in ((d2V, d1V, FV, pV, rV), (d2U, d1U, FU, pU, rU)):
            terms = [d2(0), Fx * d1(0), rx * d1(1), 0.5 * EPS * px ** 2 * f(0, 0, 0), EPS * r * px * d1(2)]
            mc = max(mc, abs(sum(terms)))
            mt = max(mt, max(abs(t) for t in terms))
    return mc, mt


@nb.njit(cache=True)
def march_segment(L0, L1, j0, j1, h, par, tab):
    """Domain of dependence of the slice nodes j = j0..j1 (j0 > 0): nodes (j, i), j0 <= j <= j1, -j1 <= i <= -j0, j + i >= 0. Ua[q, j - j0, i + j1]."""
    n = j1 - j0 + 1
    Ua = np.zeros((3, n, n))
    for q in range(3):
        for j in range(j0, j1 + 1):
            Ua[q, j - j0, -j + j1] = L0[q, j - j0]
        for j in range(j0 + 1, j1 + 1):
            Ua[q, j - j0, 1 - j + j1] = L1[q, j - j0]
    B0 = np.zeros((3, 3))
    Bp = np.zeros((3, 3))
    Bm = np.zeros((3, 3))
    M = np.zeros((3, 3))
    rhs = np.zeros(3)
    out = np.zeros(3)
    Sv = np.zeros(3)
    Ev = np.zeros(3)
    Wv = np.zeros(3)
    for ll in range(2, j1 - j0 + 1):
        guess = 0.0
        for j in range(ll + j0, j1 + 1):
            i = ll - j
            for q in range(3):
                Sv[q] = Ua[q, j - 1 - j0, i - 1 + j1]
                Ev[q] = Ua[q, j - j0, i - 1 + j1]
                Wv[q] = Ua[q, j - 1 - j0, i + j1]
            guess = hs._cellK(Sv, Ev, Wv, (j - 1) * h, (i - 1) * h, h, h, par, tab, guess, B0, Bp, Bm, M, rhs, out)
            for q in range(3):
                Ua[q, j - j0, i + j1] = out[q]
    return Ua


S0M, S1M = 0.3, 1.6          # the slice segment of the planted growing mode (exterior only: g = rho R/sqrt2 vanishes at S)


class ModeData:
    """The exact growing mode Y = c2 V R(rho(Pi)) of #262 (W = kappa, exact for c = -1) on the slice segment [0.3, 1.6], in the gauge dF = dF_N = 0.
    Y = dr - g dphi, Y_N = n - g dphi_N, g = rho R/sqrt2. The two constraint ODEs, solved for dphi and dphi_N with dr = Y + g dphi and n = Y_N + g dphi_N, give the data.
    Their homogeneous solutions (Y = 0) are not all pure gauge: the static term M of #232 (u_pm + Vhat u/2, u = Y/rho, zero on the exact mode) is carried by one
    combination of them. The initial value dphi(s0) is fixed so that M = 0 (dphi'(s0) = dphi_N(s0) = 0); then the data are those of the exact mode up to gauge."""

    def __init__(self, sl, c2=1e-3, ic=None):
        self.sl, self.c2 = sl, c2
        bg = sl.bg
        self.s0, self.s1 = S0M, S1M

        def aux(s):
            b = bg.at(s)
            rho, R = b["rho"], b["r"]
            rs = 2 * s * KAPPA * b["E"]                              # d rho / d s = 2 s G, G = kappa E
            rss = 2 * KAPPA * b["E"] * (1 - 2 * s * s * b["Fp"])
            gr = (1 + 2 * rho * rho) / (math.sqrt(2) * R)
            grr = rho * (3 + 2 * rho * rho) / (math.sqrt(2) * R ** 3)
            g, gs, gss = rho * R / math.sqrt(2), gr * rs, grr * rs * rs + gr * rss
            Rs = rho / R * rs
            Rss = rs * rs / R ** 3 + rho / R * rss
            return dict(g=g, gs=gs, gss=gss, Y=self.c2 * s * R, Ys=self.c2 * (R + s * Rs), Yss=self.c2 * (2 * Rs + s * Rss), YN=self.c2 * R, YNs=self.c2 * Rs,
                        rho=rho, R=R, E=b["E"], Fp=b["Fp"], rs=rs, gr=gr, grr=grr)
        self.aux = aux

        def rhs(s, y):
            a = aux(s)
            dphi, dphiT, z = y
            dr, drs = a["Y"] + a["g"] * dphi, a["Ys"] + a["gs"] * dphi + a["g"] * dphiT
            dr2 = sl.rhs_even(s, [dr, drs], dphi, dphiT)[1]
            dphiTT = (dr2 - a["Yss"] - 2 * a["gs"] * dphiT - a["gss"] * dphi) / a["g"]
            n = a["YN"] + a["g"] * z
            nT = sl.rhs_odd(s, n, z)
            zT = (nT - a["YNs"] - a["gs"] * z) / a["g"]
            return [dphiT, dphiTT, zT]
        self._rhs = rhs
        if ic is not None:
            self.sol = self._integrate(ic)
            return
        m0 = ModeData(sl, c2, ic=(0.0, 0.0, 0.0)).M_at(1.0)
        m1 = ModeData(sl, 0.0, ic=(1e-3, 0.0, 0.0)).M_at(1.0)
        self.alpha = -m0 / m1                                           # the data with M = 0 are mode(0) + alpha * e1
        self.sol = self._integrate((self.alpha * 1e-3, 0.0, 0.0))
        self.Mmax = max(abs(self.M_at(s)) for s in np.linspace(self.s0, self.s1, 14))
        self.Mscale = max(abs(m0), 1e-300)

    @staticmethod
    def aux_static(bg, s):
        b = bg.at(s)
        rho, R = b["rho"], b["r"]
        rs = 2 * s * KAPPA * b["E"]
        rss = 2 * KAPPA * b["E"] * (1 - 2 * s * s * b["Fp"])
        gr = (1 + 2 * rho * rho) / (math.sqrt(2) * R)
        grr = rho * (3 + 2 * rho * rho) / (math.sqrt(2) * R ** 3)
        return dict(g=rho * R / math.sqrt(2), gs=gr * rs, gss=grr * rs * rs + gr * rss, rho=rho, R=R, E=b["E"], Fp=b["Fp"], rs=rs, gr=gr, grr=grr)

    def _integrate(self, ic):
        return solve_ivp(self._rhs, [self.s0, self.s1], list(ic), dense_output=True, method="DOP853", rtol=1e-12, atol=1e-18)

    def fn(self, s):
        a = self.aux(s)
        dphi, dphiT, z = self.sol.sol(s)
        dr, drs = a["Y"] + a["g"] * dphi, a["Ys"] + a["gs"] * dphi + a["g"] * dphiT
        _, dr2 = self.sl.rhs_even(s, [dr, drs], dphi, dphiT)
        dphiTT = (dr2 - a["Yss"] - 2 * a["gs"] * dphiT - a["gss"] * dphi) / a["g"]
        n = a["YN"] + a["g"] * z
        return (np.array([dr, 0.0, dphi]), np.array([drs, 0.0, dphiT]), np.array([n, 0.0, z]), np.array([dr2, 0.0, dphiTT]))

    def M_at(self, s):
        """M = (u_tt - u_xx + Vhat u)/2 = b^2 Pi u_VU + Vhat u/2 with u = Y/rho (#232 derive C1, C2), from the slice data at s."""
        bg = self.sl.bg
        U, UT, UN, UTT = self.fn(s)
        a = self.aux(s)
        rho, R, E, Fp = a["rho"], a["R"], a["E"], a["Fp"]
        G = KAPPA * E
        rhoV, rhoU = s * G, -s * G
        rhoVU = -G * (1 - s * s * Fp)
        g, gr, grr = a["g"], a["gr"], a["grr"]
        gV, gU = gr * rhoV, gr * rhoU
        gVU = grr * rhoV * rhoU + gr * rhoVU
        Vd, Ud = 0.5 * (UT + UN), 0.5 * (UN - UT)
        B0, Bp, Bm = np.zeros((3, 3)), np.zeros((3, 3)), np.zeros((3, 3))
        hs.coefK(s, -s, bg.par, bg.tab, 0.0, B0, Bp, Bm)
        UVU = B0 @ U + Bp @ Vd + Bm @ Ud
        Y = U[0] - g * U[2]
        YV = Vd[0] - g * Vd[2] - gV * U[2]
        YU = Ud[0] - g * Ud[2] - gU * U[2]
        YVU = UVU[0] - g * UVU[2] - gV * Ud[2] - gU * Vd[2] - gVU * U[2]
        uVU = YVU / rho - (YV * rhoU + YU * rhoV) / rho ** 2 + Y * (2 * rhoU * rhoV / rho ** 3 - rhoVU / rho ** 2)
        A, Ap = hs.lin._A(rho, bg.par), hs.lin._Ap(rho, bg.par)
        Vhat = 2 * A * A / rho ** 2 - (1 + 1 / R ** 2) * A * Ap / rho + A * A / R ** 4 - 2 * A / R ** 2
        return 2 * KAPPA ** 2 * s * s * uVU + 0.5 * Vhat * Y / rho

    def exact_Y(self, bg, V, U):
        x, d, r, rV, rU, pV, pU, V0, V1, V2, E, FV, FU, rhoV, rhoU = hs.jetsK(V, U, bg.par, bg.tab, 0.0)
        return self.c2 * V * r

    def y_num(self, bg, Ua, j0, j1, j, i):
        x, d, r, rV, rU, pV, pU, V0, V1, V2, E, FV, FU, rhoV, rhoU = hs.jetsK(j * self._h, i * self._h, bg.par, bg.tab, 0.0)
        return Ua[0, j - j0, i + j1] - x * r / math.sqrt(2.0) * Ua[2, j - j0, i + j1]


def run_mode(bg, md, h):
    j0, j1 = int(round(S0M / h)), int(round(S1M / h))
    js = list(range(j0, j1 + 1))
    L0, L1 = layer_pair(bg, h, js, md.fn)
    Ua = march_segment(L0, L1, j0, j1, h, bg.par, bg.tab)
    md._h = h
    err, scale = 0.0, 0.0
    for j in range(j0, j1 + 1, max(1, (j1 - j0) // 60)):
        for i in range(-j1, -j0 + 1, max(1, (j1 - j0) // 60)):
            if j + i < 0 or i < -j1 or j + i > j1 - j0 + 0 and False:
                continue
            if i > -j0 or i < -j1 or j + i < 0:
                continue
            ex = md.exact_Y(bg, j * h, i * h)
            err = max(err, abs(md.y_num(bg, Ua, j0, j1, j, i) - ex))
            scale = max(scale, abs(ex))
    return err / scale


# ----------------------------------------------------------------------------------------------------------------------
# exact pure-gauge solutions L_xi(member) for the planted cases
# ----------------------------------------------------------------------------------------------------------------------
XI = lambda V: (0.3 + 0.7 * V + 0.2 * V * V + 0.1 * V ** 3, 0.7 + 0.4 * V + 0.3 * V * V)
ETA = lambda U: (0.03 + 0.5 * U - 0.3 * U * U + 0.15 * U ** 3, 0.5 - 0.6 * U + 0.45 * U * U)


def gauge_fun(bg):
    def gf(V, U):
        return hs.gaugeK(bg.par, bg.tab, V, U, *XI(V), *ETA(U))
    return gf


def fd_slice_fn(gf, hfd=1e-3):
    """slice data (U, U_T, U_N, U_TT) of an exact solution gf(V, U) by 4th-order differences with step hfd."""
    def f(s, a, b):
        return gf(s + a, -s + b)

    def fn(s):
        U0 = f(s, 0, 0)

        def d1(dV, dU):
            v = [f(s, m * hfd * dV, m * hfd * dU) for m in (-2, -1, 1, 2)]
            return (v[0] - 8 * v[1] + 8 * v[2] - v[3]) / (12 * hfd)
        UT, UN = d1(1, -1), d1(1, 1)
        v = [f(s, m * hfd, -m * hfd) for m in (-2, -1, 0, 1, 2)]
        UTT = (-v[0] + 16 * v[1] - 30 * v[2] + 16 * v[3] - v[4]) / (12 * hfd ** 2)
        return U0, UT, UN, UTT
    return fn


def known_c2_gauge(c, eta0=0.03):
    """c2 of dr for the pure gauge of ETA at c: -eta(0) R'(rho_h) G(rho_h), G = dRho/dPi at Pi = 0 = E(rho_h) with Pi(rho = 1) = 1 (mpmath, 60 digits; #248's N6 quotes -1.12925016e-2 for c = -1.5)."""
    import mpmath as mp
    mp.mp.dps = 60
    cc = mp.mpf(c)
    K = -2 * cc / mp.pi
    A = lambda r: 1 + cc * (r * r + 1) + K * (r + (r * r + 1) * mp.atan(r))
    Ap = lambda r: 2 * cc * r + 2 * K * (1 + r * mp.atan(r))
    rh = mp.findroot(A, (mp.mpf(0), mp.mpf(1)), solver="illinois")
    kap = Ap(rh) / 2
    lnE_h = mp.log(1 - rh) + mp.quad(lambda r: 1 / (r - rh) - 2 * kap / A(r), [1, (1 + rh) / 2, rh], method="gauss-legendre")   # interior nodes only: the integrand is a difference of poles
    G = mp.e ** lnE_h
    return float(-eta0 * (rh / mp.sqrt(1 + rh * rh)) * G)


# ----------------------------------------------------------------------------------------------------------------------
# the data of the checks
# ----------------------------------------------------------------------------------------------------------------------
def data_G():
    return (Data([(1e-3, 0.30, 0.45), (-0.5e-3, -0.25, 0.45)]), Data([(1e-3, 0.30, 0.45), (0.5e-3, -0.25, 0.45)]),
            Data([(1.0, 0.0, 0.6)]), Data([(1.0, 0.3, 0.3), (-1.0, -0.3, 0.3)]))


def data_R():
    return (Data([(1e-3, 0.775, 0.225), (-0.6e-3, 0.70, 0.15)]), Data([(1e-3, 0.775, 0.225), (0.4e-3, 0.70, 0.15)]),
            Data([(1.0, 0.775, 0.225)]), Data([(1.0, 0.70, 0.15), (-1.0, 0.85, 0.15)]))


def data_O1():
    return (Data([(1e-3, 0.8, 0.2), (-1e-3, -0.8, 0.2)]), Data([(1e-3, 0.8, 0.2), (1e-3, -0.8, 0.2)]), Data([]), Data([]))


def data_X20():
    rng = np.random.default_rng(258)
    corr = (Data([(1.0, 0.0, 0.6)]), Data([(1.0, 0.3, 0.3), (-1.0, -0.3, 0.3)]))
    out = []
    for k in range(20):
        sp = Data([(rng.normal() * 1e-3, rng.uniform(-0.5, 0.5), rng.uniform(0.2, 0.45)) for _ in range(3)])
        spN = Data([(rng.normal() * 1e-3, rng.uniform(-0.5, 0.5), rng.uniform(0.2, 0.45)) for _ in range(3)])
        out.append((sp, spN) + corr)
    return out


_DAT = {}


def datum(name, mode=None):
    """cached Datum; the mutations asym_S, mass_off, odd_off, left_mass act on G only, the mutated ODEs act on all."""
    key = (name, mode)
    if key not in _DAT:
        sl = Slice(bgc(), MUT if MUT in ("H_P", "H_eps", "M_phi", "M_F", "kg_eps") else None)
        if name == "G":
            _DAT[key] = Datum(sl, *data_G(), mode=mode, name="G")
        elif name == "R":
            _DAT[key] = Datum(sl, *data_R(), name="R")
        elif name == "O1":
            _DAT[key] = Datum(sl, *data_O1(), name="O1")
    return _DAT[key]


def G_datum():
    return datum("G", MUT if MUT in ("mass_off", "odd_off", "asym_S", "left_mass") else None)


_EVO = {}


def evolve_G(h):
    """right-half evolution of G at step h: (c1, c2, linearity residual, max |C|, largest term)."""
    if h not in _EVO:
        d = G_datum()
        bg = bgc()
        Ua, J = run_half(bg, d.slice_fn(), h)
        c1, c2, res = fit_half(Ua, J, h)
        mc, mt = constraints_at(bg, Ua, J, h) if h <= 0.005 else (float("nan"), float("nan"))
        _EVO[h] = (c1, c2, res, mc, mt)
        del Ua
    return _EVO[h]


# ----------------------------------------------------------------------------------------------------------------------
# checks
# ----------------------------------------------------------------------------------------------------------------------
RESULTS = []
SCOPE = ("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "O1")


def want(cid):
    return ONLY is None or cid in ONLY


def check(cid, text, ok, detail=""):
    RESULTS.append((cid, bool(ok)))
    print(f"[{'pass' if ok else 'FAIL'}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def ratios(vals):
    return [vals[k] / vals[k + 1] for k in range(len(vals) - 1)]


def sci(x):
    return f"{x:.3e}"


def check_Q1():
    worst = 0.0
    det = []
    for c in (-1.0, -1.5):
        bg = bgc(c)
        sl = Slice(bg, MUT if MUT in ("H_P", "H_eps", "M_phi", "M_F") else None)
        fn = fd_slice_fn(gauge_fun(bg))
        cache = {}

        def F(s):
            k = round(s, 12)
            if k not in cache:
                cache[k] = fn(s)
            return cache[k]
        s0, s1 = 0.3, 1.6
        U0, UT0, UN0, _ = F(s0)
        tol = dict(dense_output=True, method="DOP853", rtol=1e-11, atol=1e-16)
        sE = solve_ivp(lambda s, y: sl.rhs_even(s, y, F(s)[0][2], F(s)[1][2], F(s)[0][1], F(s)[1][1]), [s0, s1], [U0[0], UT0[0]], **tol)
        sO = solve_ivp(lambda s, y: [sl.rhs_odd(s, y[0], F(s)[2][2], F(s)[2][1])], [s0, s1], [UN0[0]], **tol)
        ss = np.linspace(s0, s1, 14)
        ex = np.array([F(s)[0][0] for s in ss])
        exn = np.array([F(s)[2][0] for s in ss])
        eA = np.max(np.abs(sE.sol(ss)[0] - ex)) / np.max(np.abs(ex))
        eB = np.max(np.abs(sO.sol(ss)[0] - exn)) / np.max(np.abs(exn))
        worst = max(worst, eA, eB)
        det.append(f"c={c}: dr {eA:.1e}, n {eB:.1e}")
    check("Q1", "the slice ODEs reproduce the exact pure-gauge dr and n (with dF sources), relative", worst <= 1e-9, "; ".join(det))


def check_Q2():
    bg = bgc(-1.0, schw=True)
    eps_s = 1e-3

    def gf(V, U):
        x, d, r, rV, rU, pV, pU, V0, V1, V2, E, FV, FU, rhoV, rhoU = hs.jetsK(V, U, bg.par, bg.tab, 0.0)
        g = hs.gaugeK(bg.par, bg.tab, V, U, *XI(V), *ETA(U))
        return g + np.array([eps_s * r, -2 * eps_s, 0.0])
    fn = fd_slice_fn(gf)
    errs = []
    for s in (0.2, 1.0, 5.0):
        U0, UT, UN, _ = fn(s)
        errs.append(abs(dm_of(bg, s, U0[0], UT[0], U0[1]) - eps_s / 2))
    check("Q2", "delta m = eps/2 on the Schwarzschild scaling solution plus pure gauge, s = 0.2, 1, 5", max(errs) <= 1e-10, "errors " + ", ".join(f"{e:.1e}" for e in errs))


def check_Q3Q4Q5Q7():
    d = G_datum()
    sl = d.sl
    endR = d.ev(SR)
    endL2 = d.ev2(SL)
    tR = tails(sl, endR, +1)
    tL = tails(sl, endL2, -1)
    v = lambda t: max(max(abs(x[1]), abs(x[2]), abs(KAPPA * x[0] * x[3])) for x in t)
    tail = max(v(tR), v(tL))
    dmax = max(abs(dm_of(sl.bg, x[0], x[1], x[2])) for x in tR[1:] + tL[1:])
    dmR = dm_of(sl.bg, SFAR, tR[-1][1], tR[-1][2])
    if want("Q3"):
        check("Q3", "the solution vanishes outside the support: tails of |dr|, |dr'|, |kappa s n| on 1.05 <= |s| <= 3e4", tail <= 1e-12,
              f"right {v(tR):.2e}, left {v(tL):.2e}")
    if want("Q4"):
        check("Q4", "delta M_ADM: |dm(s)| for 10 <= |s| <= 3e4 at both ends", dmax <= 1e-12, f"max {dmax:.2e}, dm_R(3e4) = {dmR:.3e}")
    if want("Q5"):
        dc2 = abs(d.c2 - d.c2_2)
        check("Q5", "left start against right start: |c2 - c2'|, |c1| (both starts)", max(dc2, abs(d.c1), abs(d.c1_2)) <= 1e-12,
              f"c2 {d.c2:.9e} against {d.c2_2:.9e}, c1 {d.c1:.1e}, {d.c1_2:.1e}")
    if want("Q7"):
        thR = abs(sl.theta(SR, endR[0], endR[1]))
        thL = abs(sl.theta(SL, endL2[0], endL2[1]))
        h0 = sl.hom_even(SL, SFAR, (1.0, 0.0))
        h1 = sl.hom_even(SL, SFAR, (0.0, 1.0))
        var = 0.0
        for hh in (h0, h1):
            t0 = sl.theta(SR, *hh.sol(SR))
            for s in (10.0, 1e2, 1e3, 1e4, SFAR):
                var = max(var, abs(sl.theta(s, *hh.sol(s)) - t0) / abs(t0))
        check("Q7", "first integral Theta: conserved on the homogeneous solutions (relative), and zero past the support", var <= 1e-9 and max(thR, thL) <= 1e-12,
              f"variation {var:.1e}, Theta(1.05) = {thR:.1e}, Theta(-1.05, right start) = {thL:.1e}")


def check_Q6():
    sl = Slice(bgc(), MUT if MUT in ("H_P", "H_eps", "kg_eps") else None)
    worst = 0.0

    def fd(fun, s, h=1e-3):
        return (-fun(s + 2 * h) + 8 * fun(s + h) - 8 * fun(s - h) + fun(s - 2 * h)) / (12 * h)
    for s in (0.1, 0.5, 1.0, 2.0, 4.0, 8.0):
        mu, a, c, Q, b = sl.theta_coef(s)
        a_s = fd(lambda t: sl.theta_coef(t)[1], s)
        c_s = fd(lambda t: sl.theta_coef(t)[2], s)
        worst = max(worst, abs(a_s - mu * Q) / abs(mu * Q), abs(c_s - mu * b["E"] * b["V1"] * b["r"]) / abs(mu * b["E"] * b["V1"] * b["r"]))
    check("Q6", "Theta' = mu H: a' = mu Q and c' = mu E V' r (the Einstein and the static scalar equation), relative", worst <= 1e-7, f"worst {worst:.1e}")


def check_E1E2():
    # E1: planted evolution on exact pure gauge; E2: planted c2
    det1, ok1 = [], True
    det2, ok2 = [], True
    for c in (-1.0, -1.5):
        bg = bgc(c)
        gf = gauge_fun(bg)
        fn = fd_slice_fn(gf)
        errs, c2f = [], {}
        for h in HLIST:
            Ua, J = run_half(bg, fn, h)
            scale = max(abs(fn(j * h)[0][0]) for j in range(0, J + 1, max(1, J // 150)))
            e = 0.0
            for (V, U) in CPTS[:7] + [(3.0, 0.0)]:
                j, i = int(round(V / h)), int(round(U / h))
                e = max(e, float(np.max(np.abs(Ua[:, j, i + J] - gf(j * h, i * h)))))
            errs.append(e / scale)
            c1, c2, rr = fit_half(Ua, J, h)
            c2f[h] = (c2, scale)
            del Ua
        rt = ratios(errs)
        ok = errs[-1] <= 3e-5 and abs(rt[-1] - 4.0) <= 0.3
        ok1 &= ok
        det1.append(f"c={c}: err/scale {errs[-1]:.1e}, last ratio {rt[-1]:.2f}")
        U0, UT0, UN0, _ = fn(0.0)
        c2s = c_functionals(UT0[0], UN0[0])[0]
        if c == -1.5:
            kn = known_c2_gauge(-1.5)
            rich = richardson({h: c2f[h][0] for h in HLIST}, [0.005, 0.0025])[0]
            e_sl, e_fit = abs(c2s - kn) / abs(kn), abs(rich - kn) / abs(kn)
            ok2 &= e_sl <= 1e-5 and e_fit <= 1e-5
            det2.append(f"c=-1.5 known {kn:.8e}: slice {e_sl:.1e}, Rich(0.005,0.0025) fit {e_fit:.1e}")
        else:
            rel = {h: abs(c2f[h][0]) / c2f[h][1] for h in HLIST}
            rt2 = ratios([rel[h] for h in HLIST])
            ok2 &= rel[0.0025] <= 1e-5 and abs(rt2[-1] - 4.0) <= 0.3
            det2.append(f"c=-1: |c2 fit|/max|dr| at 0.0025 {rel[0.0025]:.1e}, last ratio {rt2[-1]:.2f}")
    if want("E1"):
        check("E1", "planted evolution of exact pure gauge: error/scale at h = 0.00125 and last Cauchy ratio", ok1, "; ".join(det1))
    if want("E2"):
        check("E2", "planted c2: functional and evolved fit against the known coefficient (c = -1.5), seed falling by 4 (c = -1)", ok2, "; ".join(det2))


def check_E345():
    res = {h: evolve_G(h) for h in HLIST}
    c2s = G_datum().c2
    rich = richardson({h: res[h][1] for h in HLIST}, list(HLIST))
    rel = [abs(r - c2s) / max(abs(c2s), 1e-300) for r in rich]
    if want("E3"):
        lin = [res[h][2] for h in HLIST]
        rt = ratios(lin)
        check("E3", "the evolved horizon values are linear in V: residual/max at h = 0.00125 and last ratio", lin[-1] <= 1e-5 and abs(rt[-1] - 4.0) <= 0.3,
              f"residuals {', '.join(f'{x:.1e}' for x in lin)}; last ratio {rt[-1]:.2f}")
    if want("E4"):
        check("E4", "c2 from the slice against c2 of the evolved horizon fit, relative, Rich(0.005, 0.0025)", rel[1] <= 1e-6,
              f"slice {c2s:.9e}; raw fits {', '.join(f'{res[h][1]:.9e}' for h in HLIST)}; Rich relative errors {', '.join(f'{x:.1e}' for x in rel)}; c1 fit {res[0.00125][0]:.1e}")
    if want("E5"):
        mcs = [res[h][3] for h in (0.005, 0.0025, 0.00125)]
        mts = [res[h][4] for h in (0.005, 0.0025, 0.00125)]
        rt = ratios(mcs)
        check("E5", "constraints of the evolved fields at 8 fixed points: |C|/largest term at h = 0.0025 and last ratio", mcs[1] / mts[1] <= 1e-3 and abs(rt[-1] - 4.0) <= 0.3,
              f"|C| {', '.join(f'{x:.2e}' for x in mcs)}; term {mts[1]:.2e}; ratios {', '.join(f'{x:.2f}' for x in rt)}")


def check_E6():
    d = G_datum()
    bg = bgc()
    fn = d.slice_fn()
    ca, cb = {}, {}
    for h in (0.01, 0.005, 0.0025):
        Ua, J = run_diamond(bg, fn, h)
        V = h * np.arange(0, J + 1)
        Hh = np.array([Ua[0, j + J, J] for j in range(0, J + 1)])
        Hl = np.array([Ua[0, J, i + J] for i in range(0, J + 1)])
        ca[h] = hs.horizon_fit(Hh, V)[1]
        cb[h] = hs.horizon_fit(Hl, V)[1]
        del Ua
    ra = richardson(ca, [0.01, 0.005, 0.0025])
    rb = richardson(cb, [0.01, 0.005, 0.0025])
    ea, eb = abs(ra[1] - d.c2) / max(abs(d.c2), 1e-300), abs(rb[1] - d.c2p) / max(abs(d.c2p), 1e-300)
    check("E6", "diamond |s| <= 0.5 through S and the interior: c2 and c2' fits on both future horizons against the slice, relative, Rich(0.005, 0.0025)", max(ea, eb) <= 1e-6,
          f"c2 slice {d.c2:.9e} error {ea:.1e}; c2' slice {d.c2p:.9e} error {eb:.1e}")


def Yvalues(d):
    return dict(Y0=d.Y(0.0), dVY=d.c2, dUY=d.c2p, Ym=d.Y(-0.3), Yp=d.Y(0.3), Ymax=max(abs(d.Y(s)) for s in np.linspace(-1.05, 1.05, 211)),
                Ygap=max(abs(d.Y(s)) for s in np.linspace(-0.5, 0.5, 101)))


def check_E7():
    d = datum("R")
    bg = bgc()
    fn = d.slice_fn()
    fit = {}
    for h in HLIST:
        Ua, J = run_half(bg, fn, h)
        fit[h] = fit_half(Ua, J, h)[1]
        del Ua
    rt = ratios([abs(fit[h]) for h in HLIST])
    y = Yvalues(d)
    ok = abs(d.c2) <= 1e-12 and abs(fit[0.00125]) <= 1e-7 and abs(rt[-1] - 4.0) <= 0.3 and y["Ygap"] <= 1e-12
    check("E7", "one-sided datum R (support [0.55, 1.0]): c2 = 0 on the slice, evolved fit is a seed falling by 4, Y vanishes on |s| <= 0.5", ok,
          f"c2 slice {d.c2:.1e}; fit {', '.join(f'{abs(fit[h]):.2e}' for h in HLIST)}, last ratio {rt[-1]:.2f}; Y(S) {y['Y0']:.1e}, Y(+-0.3) {y['Yp']:.1e}, {y['Ym']:.1e}, max|Y| on |s|<=0.5 {y['Ygap']:.1e}")


def check_O1():
    d = datum("O1")
    bg = bgc()
    fn = d.slice_fn()
    fit = {}
    for h in HLIST[:3]:
        Ua, J = run_half(bg, fn, h)
        fit[h] = fit_half(Ua, J, h)[1]
        del Ua
    rich = richardson(fit, list(HLIST[:3]))
    rel = abs(rich[1] - d.c2) / max(abs(d.c2), 1e-300)
    y = Yvalues(d)
    ends = d.ev(SR)
    ok = abs(d.c2) >= 1e-6 and rel <= 1e-6 and y["Ygap"] >= 1e-6 and max(abs(ends[0]), abs(ends[1]), abs(ends[2])) <= 1e-12
    check("O1", "two-sided datum supported at 0.6 <= |s| <= 1.0 (no condition needed): c2 != 0, fit equal to the slice, Y nonzero in the gap", ok,
          f"c2 slice {d.c2:.6e} (amplitude 1e-3), fit error {rel:.1e}; Y(S) {y['Y0']:.1e}, dV Y(S) {y['dVY']:.4e}, dU Y(S) {y['dUY']:.4e}, Y(+-0.3) {y['Yp']:.3e}, {y['Ym']:.3e}, max|Y| on |s|<=0.5 {y['Ygap']:.3e}; end vector {max(abs(ends)):.1e}")


def check_E8():
    bg = bgc()
    sl = Slice(bg, MUT if MUT in ("H_P", "H_eps", "M_phi", "M_F", "kg_eps") else None)
    worst, cmin, cs = 0.0, 1.0, []
    for k, parts in enumerate(data_X20()):
        d = Datum(sl, *parts, name=f"X{k}")
        fn = d.slice_fn()
        fit = {}
        for h in (0.005, 0.0025):
            Ua, J = run_half(bg, fn, h)
            fit[h] = fit_half(Ua, J, h)[1]
            del Ua
        r = richardson(fit, [0.005, 0.0025])[0]
        worst = max(worst, abs(r - d.c2))
        cmin = min(cmin, abs(d.c2))
        cs.append(d.c2)
    check("E8", "20 random compact data: |Rich(0.005, 0.0025) - slice| absolute, min |c2|", worst <= 1e-6 and cmin >= 1e-5,
          f"max error {worst:.2e}; |c2| min {cmin:.2e}, median {np.median(np.abs(cs)):.2e}, max {np.max(np.abs(cs)):.2e}")


def check_E9():
    bg = bgc()
    sl = Slice(bg, MUT if MUT in ("H_P", "H_eps", "M_phi", "M_F") else None)
    md = ModeData(sl)
    errs = [run_mode(bg, md, h) for h in HLIST]
    rt = ratios(errs)
    check("E9", "planted exact growing mode Y = c2 V R (#262), data from Y by the two constraint ODEs on [0.3, 1.6] with M = 0, evolved: error of Y over its maximum at h = 0.00125, last ratio",
          errs[-1] <= 1e-5 and abs(rt[-1] - 4.0) <= 0.3 and md.Mmax <= 1e-9 * md.Mscale,
          f"errors {', '.join(f'{e:.2e}' for e in errs)}; last ratio {rt[-1]:.2f}; M of the data: max {md.Mmax:.1e} against {md.Mscale:.1e} before the fix (alpha = {md.alpha:.6f})")


CHECKS = [("Q1", check_Q1), ("Q2", check_Q2), ("Q3", check_Q3Q4Q5Q7), ("Q6", check_Q6), ("E1", check_E1E2), ("E3", check_E345), ("E6", check_E6),
          ("E7", check_E7), ("O1", check_O1), ("E8", check_E8), ("E9", check_E9)]
GROUPS = {"Q1": ["Q1"], "Q2": ["Q2"], "Q3": ["Q3", "Q4", "Q5", "Q7"], "Q6": ["Q6"], "E1": ["E1", "E2"], "E3": ["E3", "E4", "E5"], "E6": ["E6"], "E7": ["E7"], "O1": ["O1"], "E8": ["E8"], "E9": ["E9"]}


def run_checks():
    t0 = time.time()
    for key, fnc in CHECKS:
        if ONLY is not None and not any(c in ONLY for c in GROUPS[key]):
            continue
        t1 = time.time()
        fnc()
        print(f"        ({key} group {time.time() - t1:.0f} s)", flush=True)
    bad = [c for c, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS) - len(bad)} of {len(RESULTS)} checks pass" + (f"; fired: {', '.join(bad)}" if bad else "") + f"  [{time.time() - t0:.0f} s]  mutation: {MUT}")
    return 1 if bad else 0


# ----------------------------------------------------------------------------------------------------------------------
# declared fired sets and the comparison with stored logs
# ----------------------------------------------------------------------------------------------------------------------
DECLARED = {                      # computed in the development run (dev1) before the runs of record; every check is in scope
    "H_P": {"Q1", "Q3", "Q4", "Q5", "Q6", "Q7", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "O1"},
    "H_eps": {"Q1", "Q3", "Q4", "Q5", "Q6", "Q7", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "O1"},
    "M_phi": {"Q1", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "O1"},
    "M_F": {"Q1", "Q3", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "O1"},
    "fn_aT": {"E2", "E4", "E6", "E8", "O1"},
    "fn_sign": {"E2", "E4", "E6", "E8", "O1"},
    "l1_sign": {"E1", "E2", "E4", "E6", "E8", "E9", "O1"},
    "l1_noquad": {"E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "O1"},
    "dm_noT": {"Q2"},
    "mass_off": {"Q3", "Q4", "Q5", "E3", "E4"},
    "odd_off": {"Q3", "Q5", "E3", "E4"},
    "asym_S": {"Q3", "Q4", "E3", "E4", "E6"},
    "left_mass": {"Q3", "Q4", "Q5", "Q7", "E3", "E4"},
    "kg_eps": {"Q6"},
}


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
        extra = sorted(fired - set(SCOPE))
        ok = inscope == DECLARED[name]
        bad += 0 if ok else 1
        print(f"{name}: fired {sorted(inscope)}; declared {sorted(DECLARED[name])}; {'equal' if ok else 'DIFFERENT'}" + (f"; also fired outside the declared scope: {extra}" if extra else ""))
    print(f"{len(KNOWN) - bad} of {len(KNOWN)} fired sets equal the declared ones")
    return 1 if bad else 0


# ----------------------------------------------------------------------------------------------------------------------
# the table
# ----------------------------------------------------------------------------------------------------------------------
def report_table():
    bg = bgc()
    sl = Slice(bg)
    print("== data: c1 = dr(S), c2 = d_V dr(S), c2' = d_U dr(S), Y at S and in the gap (amplitude 1e-3) ==")
    for nm in ("G", "R", "O1"):
        d = datum(nm)
        y = Yvalues(d)
        print(f"{nm}: support amp {d.dphi.amp():.2e}, k_m {d.km:.4e}, k_o {d.ko:.4e}; c1 {d.c1:.2e}; c2 {d.c2:.9e}; c2' {d.c2p:.9e}; Y(S) {y['Y0']:.2e}; Y(+0.3) {y['Yp']:.4e}; Y(-0.3) {y['Ym']:.4e}; max|Y| (|s|<=0.5) {y['Ygap']:.3e}")
    print("== the static term M of #232 on the slice (zero on the exact mode) ==")
    for nm, mode in (("G", None), ("G, C1 not imposed (mass_off)", "mass_off"), ("G, delta M_L = 1e-4 (left_mass)", "left_mass")):
        d = datum("G") if mode is None else Datum(sl, *data_G(), mode=mode, name=nm)
        md = ModeData.__new__(ModeData)
        md.sl, md.c2 = sl, 0.0
        md.fn = d.slice_fn()
        md.aux = lambda s_, _b=bg: ModeData.aux_static(_b, s_)
        ms = [(s_, md.M_at(s_)) for s_ in (-0.9, -0.5, -0.2, 0.2, 0.5, 0.9)]
        print(f"{nm}: M(s) at s = -0.9, -0.5, -0.2, 0.2, 0.5, 0.9: " + ", ".join(f"{m:.3e}" for _, m in ms) + f"; at s = 1.0 (outside the support): {md.M_at(1.0):.3e}")
    print("== parity: P-odd dphi + P-even dphi_N need no condition and carry c2; P-even dphi + P-odd dphi_N need C1, C2 and carry none ==")
    sp, spN, cp, cN = data_G()

    def proj(D, sgn):
        return Data([(a / 2, c, w) for a, c, w in D.parts] + [(sgn * a / 2, -c, w) for a, c, w in D.parts])
    for lab, A, B, sA, sB in (("P-odd dphi + P-even dphi_N", proj(sp, -1), proj(spN, +1), -1, +1), ("P-even dphi + P-odd dphi_N", proj(sp, +1), proj(spN, -1), +1, -1)):
        d = Datum(sl, A, B, cp, cN, name=lab)
        endv = d.ev(SR)
        print(f"{lab}: corrections k_m {d.km:.2e}, k_o {d.ko:.2e}; c2 {d.c2:.6e}; dr'(0) {d.ev(0.0)[1]:.2e}, n(0) {d.ev(0.0)[2]:.2e}; end vector {max(abs(endv)):.1e}")
    print("== first integral: Theta of the homogeneous solutions (conserved), delta m, -1/kappa ratio at rho -> infinity ==")
    for y0 in ((1.0, 0.0), (0.0, 1.0)):
        h = sl.hom_even(SL, SFAR, y0)
        th = [sl.theta(s, *h.sol(s)) for s in (SR, 1e2, 1e4, SFAR)]
        dm = dm_of(bg, SFAR, *h.sol(SFAR))
        print(f"hom start {y0}: Theta at s = 1.05, 1e2, 1e4, 3e4: {', '.join(f'{t:.8e}' for t in th)}; delta m(3e4) = {dm:.5e}; Theta/dm = {th[-1] / dm:.4f} (-1/kappa = {-1 / KAPPA:.4f})")
    print("== c1 and the left start: left_mass with delta M = 1e-4 and 2e-4 ==")
    for k in (1.0, 2.0):
        w = k * 1e-4 * Datum.mass_start(sl)
        sol = sl.even(SL, SR, data_G()[0] + data_G()[2].scaled(datum("G").km), tuple(w), max_step=0.01)
        print(f"delta M = {k * 1e-4:.0e}: c1 = {sol.sol(0.0)[0]:.6e}, c1/delta M = {sol.sol(0.0)[0] / (k * 1e-4):.6f}")
    return 0


if __name__ == "__main__":
    cmd = next((a for a in sys.argv[1:] if not a.startswith("--")), "check")
    if cmd == "check":
        sys.exit(run_checks())
    if cmd == "table":
        sys.exit(report_table())
    if cmd == "mutations":
        sys.exit(report_mutations(sys.argv[sys.argv.index("mutations") + 1], sys.argv[sys.argv.index("mutations") + 2]))
    sys.exit(f"unknown command {cmd}")

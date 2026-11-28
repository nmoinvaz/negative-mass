#!/usr/bin/env python3
"""The static coefficient c1 of the c = -1 black universe's horizon data (#259).

c1 = Y(S) = dr at the bifurcation surface S. #258 found a first integral Theta of the sum constraint on the slice V + U = 0 (gauge dF = 0) with Theta = -c1 and
Theta(infinity) = -dM_ADM/kappa. Here Theta is extended to dF != 0 (Theta_gen), put through the evolved data of #248's classes and the exact static solution r = sqrt(L) R
(#269, N6), and the exterior solution is written in closed form in rho so that dM_ADM is read at rho = 1e8.

    static_c1.py derive [--mutate=NAME]                    sympy checks (S-series)
    static_c1.py check [--only=T1,..] [--mutate=NAME]      numerical checks (T-series); exit 1 when one fails
    static_c1.py table                                     the numbers of the note
    static_c1.py mutations DIR PREFIX                      fired sets of the stored logs DIR/PREFIX<name>.log against the declared ones

Reads horizon_sector.py, full_slice.py and linear.py (read-only imports). Python: ~/Source/negative-mass/.venv/bin/python.
"""
import math
import os
import re
import sys
import time

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
_ARGV = sys.argv
sys.argv = sys.argv[:1]                     # the imported modules parse --mutate / --only themselves
import horizon_sector as hs  # noqa: E402
import full_slice as fs  # noqa: E402
import linear as lin  # noqa: E402
sys.argv = _ARGV

MUT = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mutate=")), None)
ONLY = next((a.split("=", 1)[1].split(",") for a in sys.argv[1:] if a.startswith("--only=")), None)
KNOWN_D = ("s_q", "s_pot")
KNOWN_N = ("law_norm", "th_noF", "th_nophi", "th_asign", "dm_noT", "kappa_off", "jump_off", "xform_p", "xform_q", "mode_khat", "mode_dF0", "fn_c1_noV", "fn_c1_aT", "fn_c1_sign")
KNOWN = KNOWN_D + KNOWN_N
assert MUT is None or MUT in KNOWN, MUT

C = -1.0
KC = -2.0 * C / math.pi                      # coefficient of Z = rho + R^2 arctan rho in A
KAPPA = 2.0 / math.pi
KAPPA_C = KAPPA * (1.01 if MUT == "kappa_off" else 1.0)   # the kappa that multiplies c1 when it is compared with dM_ADM
A1 = 2.0 * KC                                # A'(0) = 2 kappa
KHAT = -math.pi ** 2 / 8.0 * (0.0 if MUT == "mode_khat" else 1.0)   # kappa'/kappa of the family F_L = A + 1/L - 1 at L = 1: A''(0)/A'(0)^2
RHO_MAX = 1e8
S_CUTS = (1.0, 1.9)
S_CUTL = (-0.03,)                              # the left strip of the evolved rectangle is 0.05 wide
S_LIST = (-0.03, 0.0, 0.1, 0.25, 0.5, 1.0, 1.5, 1.9)
HL = (0.01, 0.005, 0.0025, 0.00125)
SQ2 = math.sqrt(2.0)
RECT = hs.RECT
V0R, U0R = RECT[0], RECT[2]


# ----------------------------------------------------------------------------------------------------------------------
# A(rho), A', A'' of the c = -1 member, stable at large rho (series in t = 1/rho)
# ----------------------------------------------------------------------------------------------------------------------
def Am1(rho):
    """A - 1 = (2/pi)[rho - (1 + rho^2) arccot rho]."""
    if rho >= 6.0:
        t = 1.0 / rho
        t2 = t * t
        s, pw = 0.0, t
        for j in range(45):
            s += (-1) ** j * pw * 2.0 / ((2 * j + 1) * (2 * j + 3))
            pw *= t2
        return -(2.0 / math.pi) * s
    return C * (rho * rho + 1.0) + KC * (rho + (rho * rho + 1.0) * math.atan(rho))


def A_(rho):
    return 1.0 + Am1(rho)


def Ap_(rho):
    if rho >= 6.0:
        t = 1.0 / rho
        t2 = t * t
        s, pw = 0.0, t2
        for k in range(1, 45):
            s += (-1) ** (k + 1) * pw / (2 * k + 1)
            pw *= t2
        return (4.0 / math.pi) * s
    return 2.0 * C * rho + 2.0 * KC * (1.0 + rho * math.atan(rho))


def App_(rho):
    if rho >= 6.0:
        t = 1.0 / rho
        t2 = t * t
        s, pw = 0.0, t2 * t
        for k in range(1, 45):
            s += (-1) ** k * 2.0 * k * pw / (2 * k + 1)
            pw *= t2
        return (4.0 / math.pi) * s
    return 2.0 * C + 2.0 * KC * (math.atan(rho) + rho / (1.0 + rho * rho))


# ----------------------------------------------------------------------------------------------------------------------
# Theta and the Misner-Sharp mass perturbation on the slice (s-form, chart of #248), with the mutations
# ----------------------------------------------------------------------------------------------------------------------
_SL = {}


def slice_obj(c=C):
    if c not in _SL:
        _SL[c] = fs.Slice(fs.bgc(c))
    return _SL[c]


def theta_gen(sl, s, dr, drs, dphi, dF):
    """Theta_gen = a dr + mu dr' + c dphi + mu s r_Pi dF (dF = 0 gives #258's Theta). Conserved on solutions of H with dF sources."""
    mu, a, c, Q, b = sl.theta_coef(s)
    t = (-1.0 if MUT == "th_asign" else 1.0) * a * dr + mu * drs
    if MUT != "th_nophi":
        t += c * dphi
    if MUT != "th_noF":
        t += mu * s * b["rp"] * dF
    return t


def theta_scale(sl, s, dr, drs, dphi, dF):
    """the largest of the four terms of Theta_gen (absolute value), the scale of its roundoff."""
    mu, a, c, Q, b = sl.theta_coef(s)
    return max(abs(a * dr), abs(mu * drs), abs(c * dphi), abs(mu * s * b["rp"] * dF))


def dm_s(bg, s, dr, drs, dF):
    """delta m = dr (1/2 - e^F s^2 r'^2) - r e^F s r' (s r' dF + dr') (#258)."""
    b = bg.at(s)
    eF = 1.0 / b["E"]
    grp = 0.0 if MUT == "dm_noT" else 1.0
    return dr * (0.5 - eF * s * s * b["rp"] ** 2) - grp * b["r"] * eF * s * b["rp"] * (s * b["rp"] * dF + drs)


# ----------------------------------------------------------------------------------------------------------------------
# the exterior in rho (dphi = dF = 0): dr_xx + p dr_x + q dr = 0, d/dx = A d/drho; Theta = (R/kappa)(dr_x - A' dr/2) is its first integral
# ----------------------------------------------------------------------------------------------------------------------
def xform_pq(rho):
    A, A1_, A2 = A_(rho), Ap_(rho), App_(rho)
    R2 = 1.0 + rho * rho
    p = -A1_ / 2.0 + (0.0 if MUT == "xform_p" else A * rho / R2)
    q = (1.0 if MUT == "xform_q" else -1.0) * (A / 2.0) * (A2 + rho * A1_ / R2)
    return p, q


def theta_x(rho, dr, drx):
    return math.sqrt(1.0 + rho * rho) / KAPPA * (drx - Ap_(rho) * dr / 2.0)


def dm_x(rho, dr, drx):
    """delta m at dphi = dF = 0: (1/2)(1 - A rho^2/R^2) dr - rho dr_x, written with 1 - A = -Am1 (no cancellation)."""
    R2 = 1.0 + rho * rho
    grp = 0.0 if MUT == "dm_noT" else 1.0
    return 0.5 * (-Am1(rho) + A_(rho) / R2) * dr - grp * rho * drx


def tail_ode(rho0, dr0, drx0, rho1=RHO_MAX):
    """(dr, dr_x) on [rho0, rho1] from the ODE in u = ln rho."""
    def rhs(u, y):
        rho = math.exp(u)
        A = A_(rho)
        p, q = xform_pq(rho)
        return [rho * y[1] / A, rho * (-p * y[1] - q * y[0]) / A]
    return solve_ivp(rhs, [math.log(rho0), math.log(rho1)], [dr0, drx0], method="DOP853", rtol=1e-13, atol=1e-20, dense_output=True)


def tail_closed(rho0, dr0, drx0, rho):
    """dr(rho) = sqrt(A)[C0 + kappa Theta I(rho)], I = int dr'/(R A^{3/2}) ((dr/sqrt A)' = kappa Theta/(R A^{3/2})), dr_x = A' dr/2 + kappa Theta/R (the first integral solved for dr_x)."""
    th = theta_x(rho0, dr0, drx0)
    C0 = dr0 / math.sqrt(A_(rho0))
    f = lambda u: math.exp(u) / (math.sqrt(1.0 + math.exp(2 * u)) * A_(math.exp(u)) ** 1.5)
    I = quad(f, math.log(rho0), math.log(rho), epsabs=1e-15, epsrel=1e-13, limit=800)[0] if rho > rho0 else 0.0
    dr = math.sqrt(A_(rho)) * (C0 + KAPPA * th * I)
    return dr, Ap_(rho) * dr / 2.0 + KAPPA * th / math.sqrt(1.0 + rho * rho)


def read_mass(bg, s_cut, dr, drs, dphi, dF):
    """Hard cut of the sources at s_cut: dphi = dF = 0 beyond, the jump of dr' from H, the exterior solution to rho = RHO_MAX.
    Returns dict(theta_before, theta_after, dm(rho_k), rho_k, dM = dm(RHO_MAX), sol)."""
    b = bg.at(s_cut)
    sl = slice_obj()
    th_before = theta_gen(sl, s_cut, dr, drs, dphi, dF)
    jump = 0.0 if MUT == "jump_off" else fs.EPS * b["r"] * s_cut * b["php"] * dphi + s_cut * b["rp"] * dF
    drs_after = drs + jump
    rho0 = b["rho"]
    drx0 = fs.KAPPA * s_cut * drs_after
    sol = tail_ode(rho0, dr, drx0)
    rhos = [rho0, 1e2, 1e4, 1e6, RHO_MAX]
    dms = [dm_x(r_, *sol.sol(math.log(r_))) for r_ in rhos]
    th_after = theta_x(rho0, dr, drx0)
    return dict(th_before=th_before, th_after=th_after, rhos=rhos, dms=dms, dM=dms[-1], rho0=rho0, dr0=dr, drx0=drx0, sol=sol)


# ----------------------------------------------------------------------------------------------------------------------
# the exact static mode r = sqrt(L) R, F_L = A + 1/L - 1 (#269, N6), d/dL at L = 1, in the Kruskal chart at fixed Pi
# ----------------------------------------------------------------------------------------------------------------------
def _a_series(n):
    a = np.zeros(n + 1)
    a[1] = A1
    a[2] = C * 1.0 if n >= 2 else 0.0                    # -rho^2 (c = -1): A = c rho^2 + K Z
    for k in range(1, n):
        if 2 * k + 1 <= n:
            a[2 * k + 1] = (-1) ** (k + 1) * 2.0 * KC / (4.0 * k * k - 1.0)
    return a


class StaticMode:
    """q = d rho/dL at fixed Pi solves 2 kappa Pi q' = A'(rho) q - 1 - khat A(rho), regular at Pi = 0 (q(0) = 1/A'(0)); normalization q'(0) = 0 (the homogeneous
    solution is the scaling gauge xi = V d_V + U d_U). dr = R/2 + R' q, dphi = sqrt2 q/R^2, dF = khat - dq/drho, c1 = dr(S) = 1/2, dM_ADM = 1/pi per unit dL."""
    NS = 20
    PS = 0.03

    def __init__(self):
        kc = hs.kchart(C)
        n = self.NS
        a = _a_series(n + 1)
        ap = np.array([m * a[m] for m in range(1, n + 2)])           # A'(rho) coefficients (index m-1)

        def compose(coefs, rho):
            out = np.zeros(n + 1)
            pw = np.zeros(n + 1)
            pw[0] = 1.0
            for m in range(len(coefs)):
                if m > 0:
                    pw = np.convolve(pw, rho)[: n + 1]
                out += coefs[m] * pw
            return out
        rho = np.zeros(n + 1)
        rho[1] = 1.0                                               # rho_1 = 1; the chart's rho_1 is a scale z = rho_1 Pi
        for k in range(2, n + 1):
            Ak = compose(np.concatenate([[0.0], a[1: n + 2]]), rho)    # coefficients a_0 = 0, a_1, ...
            # [A(rho)]_k with rho_k = 0: the a_1 rho_k term is excluded by rho[k] = 0
            rho[k] = Ak[k] / (A1 * (k - 1))
        Aser = compose(np.concatenate([[0.0], a[1: n + 2]]), rho)
        Apser = compose(ap, rho)
        d = np.zeros(n + 1)
        d[0] = 1.0 / A1
        for k in range(2, n + 1):
            rhs = sum((Apser[j] - (A1 if j == 0 else 0.0)) * d[k - j] for j in range(1, k + 1)) - KHAT * Aser[k]
            d[k] = rhs / (A1 * (k - 1))
        self.rt, self.dt = rho, d
        self.drt = np.array([k * rho[k] for k in range(1, n + 1)])
        self.ddt = np.array([k * d[k] for k in range(1, n + 1)])
        P0 = 1e-3
        chart_rho = hs.jetsK(1.0, -P0, kc.par, kc.tab, 0.0)[0]
        self.r1 = brentq(lambda r1: np.polynomial.polynomial.polyval(r1 * P0, rho) - chart_rho, 1e-3, 20.0, xtol=1e-15, rtol=1e-15)
        self.kc = kc
        w0 = math.log(self.PS)
        rhs_w = lambda w, y: [A_(y[0]) / (2 * KAPPA), (Ap_(y[0]) * y[1] - 1.0 - KHAT * A_(y[0])) / (2 * KAPPA)]
        self.sol = {}
        for sgn, w1 in ((+1, math.log(1.6e12)), (-1, math.log(0.25))):
            r0, q0 = self._series(sgn * self.PS)[0], self._series(sgn * self.PS)[2]
            self.sol[sgn] = solve_ivp(rhs_w, [w0, w1], [r0, q0], method="DOP853", rtol=1e-13, atol=1e-18, dense_output=True)
        # large rho: q(rho) from the ODE in u = ln rho, started from the end of the Pi ODE
        rho_e, q_e = self.sol[+1].y[0, -1], self.sol[+1].y[1, -1]
        rhs_u = lambda u, y: [math.exp(u) * (Ap_(math.exp(u)) * y[0] - 1.0 - KHAT * A_(math.exp(u))) / A_(math.exp(u))]
        self.big = solve_ivp(rhs_u, [math.log(rho_e), math.log(RHO_MAX)], [q_e], method="DOP853", rtol=1e-13, atol=1e-18, dense_output=True)
        self.rho_e = rho_e

    def _series(self, Pi):
        z = self.r1 * Pi
        P = np.polynomial.polynomial.polyval
        return P(z, self.rt), self.r1 * P(z, self.drt), P(z, self.dt), self.r1 * P(z, self.ddt)

    def rq(self, Pi):
        """(rho, rho_Pi, q, q_Pi) at Pi in [-0.25, 1.6e12]."""
        if abs(Pi) <= self.PS:
            return self._series(Pi)
        rho, q = self.sol[1 if Pi > 0 else -1].sol(math.log(abs(Pi)))
        A, Ap = A_(rho), Ap_(rho)
        return rho, A / (2 * KAPPA * Pi), q, (Ap * q - 1.0 - KHAT * A) / (2 * KAPPA * Pi)

    def fields(self, Pi):
        """(dr, dr_Pi, dphi, dF) of the mode for unit dL."""
        rho, rho_P, q, q_P = self.rq(Pi)
        R = math.sqrt(1.0 + rho * rho)
        Rp, Rpp = rho / R, 1.0 / R ** 3
        dr = R / 2.0 + Rp * q
        dr_P = Rp * rho_P / 2.0 + Rpp * rho_P * q + Rp * q_P
        dphi = SQ2 * q / R ** 2
        dF = 0.0 if MUT == "mode_dF0" else KHAT - q_P / rho_P
        return dr, dr_P, dphi, dF

    def dm_exact(self, rho, q):
        """d_L m at fixed Pi: m_L(rho) = sqrt(L)(1 + L rho^2 (1 - A))/(2R), dm = d_L m|rho + m_rho q."""
        R = math.sqrt(1.0 + rho * rho)
        A, Ap = A_(rho), Ap_(rho)
        omA = -Am1(rho)
        dmL = 1.0 / (4 * R) + 3.0 * rho ** 2 * omA / (4 * R)
        m_rho = (2 * rho * omA - rho ** 2 * Ap) / (2 * R) - rho * (1.0 + rho ** 2 * omA) / (2 * R ** 3)
        return dmL + m_rho * q

    def q_of_rho(self, rho):
        return self.big.sol(math.log(rho))[0] if rho > self.rho_e else None


def exact_dM_mp():
    """d/dL of lim_{rho->oo} m_L(rho), m_L = sqrt(L)(1 + L rho^2 (1 - A))/(2R), by mpmath at rho = 1e14 (central difference in L; the O(1/rho) remainder is 1e-14)."""
    import mpmath as mp
    mp.mp.dps = 100
    rho = mp.mpf(10) ** 14
    c, K = mp.mpf(C), 2 / mp.pi
    A = 1 + c * (rho ** 2 + 1) + K * (rho + (rho ** 2 + 1) * mp.atan(rho))
    R = mp.sqrt(1 + rho ** 2)
    m = lambda L: mp.sqrt(L) * (1 + L * rho ** 2 * (1 - A)) / (2 * R)
    e = mp.mpf(10) ** -15
    return float((m(1 + e) - m(1 - e)) / (2 * e))


# ----------------------------------------------------------------------------------------------------------------------
# evolved slices of the Kruskal solver: readings at the slice V + U = 0
# ----------------------------------------------------------------------------------------------------------------------
def slice_vals(U, h, s):
    j = int(round((s - V0R) / h))
    i = int(round((-s - U0R) / h))
    return U[0, j, i], hs.fd4(U[0], j, i, 1, -1, h), U[2, j, i], U[1, j, i]


def kruskal_readings(U, h, ih, amp=1.0):
    """c1 and c2 of the horizon fit, dr at S, the slice functional of c2 and c1 at V0 = 0.5, Theta_gen on the slice, and the mass of the cut extensions."""
    bg = fs.bgc(C)
    sl = slice_obj()
    V = V0R + h * np.arange(U.shape[1])
    c1, c2, res = hs.horizon_fit(U[0, :, ih], V)
    jS = int(round(-V0R / h))
    out = dict(c1fit=c1, c2fit=c2, res=res, c1S=U[0, jS, ih], theta={}, cut={})
    for s in S_LIST:
        dr, drs, dphi, dF = slice_vals(U, h, s)
        out["theta"][s] = theta_gen(sl, s, dr, drs, dphi, dF)
    for sc in S_CUTS + S_CUTL:
        dr, drs, dphi, dF = slice_vals(U, h, sc)
        rm = read_mass(bg, sc, dr, drs, dphi, dF)
        out["cut"][sc] = dict(dM=rm["dM"], th_before=rm["th_before"], th_after=rm["th_after"], dms=rm["dms"])
    for V0 in (0.5, 1.0):
        j0 = int(round((V0 - V0R) / h))
        c2s = c2_slice(U[0], j0, ih, h)
        out[f"c2slice{V0}"] = c2s
        out[f"c1slice{V0}"] = c1_functional(U[0][j0, ih], c2s, V0)
    return out


_RUNS = {}


def class_run(name, h):
    key = (name, h)
    if key not in _RUNS:
        su = hs.KSetup(C, h, *RECT, hs.CLASSES[name], kc=hs.kchart(C))
        U = su.run()
        _RUNS[key] = kruskal_readings(U, h, su.ih)
        del U
    return _RUNS[key]


_MODE = {}


def mode():
    if "m" not in _MODE:
        _MODE["m"] = StaticMode()
    return _MODE["m"]


def mode_arrays(h, amp):
    """exact ray data of the static mode (3, n+1) on V = V0R (UL) and U = U0R (UR); returns also the evolved rectangle."""
    m = mode()
    Nv = int(round((RECT[1] - RECT[0]) / h))
    Nu = int(round((RECT[3] - RECT[2]) / h))
    UL = np.zeros((3, Nu + 1))
    UR = np.zeros((3, Nv + 1))
    for i in range(Nu + 1):
        dr, _, dphi, dF = m.fields(-V0R * (U0R + i * h))
        UL[:, i] = amp * np.array([dr, dF, dphi])
    for j in range(Nv + 1):
        dr, _, dphi, dF = m.fields(-(V0R + j * h) * U0R)
        UR[:, j] = amp * np.array([dr, dF, dphi])
    return UL, UR


def mode_run(h, amp=2e-3):
    key = ("mode", h, amp)
    if key not in _RUNS:
        kc = hs.kchart(C)
        UL, UR = mode_arrays(h, amp)
        U = hs.march_K(UL, UR, V0R, U0R, kc.par, kc.tab, h, h)
        ih = int(round(-U0R / h))
        rd = kruskal_readings(U, h, ih)
        m = mode()
        st = max(1, int(round(0.05 / h)))
        err, scale = 0.0, 0.0
        for j in range(0, U.shape[1], st):
            for i in range(0, U.shape[2], st):
                Pi = -(V0R + j * h) * (U0R + i * h)
                if Pi < -0.24 or Pi > 1e9:
                    continue
                ex = amp * np.array([m.fields(Pi)[0], m.fields(Pi)[3], m.fields(Pi)[2]])
                err = max(err, float(np.max(np.abs(U[:, j, i] - ex))))
                scale = max(scale, float(np.max(np.abs(ex))))
        rd["err"], rd["scale"] = err, scale
        _RUNS[key] = rd
        del U
    return _RUNS[key]


def rich(v1, v2):
    return (4.0 * v2 - v1) / 3.0


def ratios(vals):
    return [abs(vals[k]) / max(abs(vals[k + 1]), 1e-300) for k in range(len(vals) - 1)]


# ----------------------------------------------------------------------------------------------------------------------
# #232's chart: c1 as a functional of the slice data
# ----------------------------------------------------------------------------------------------------------------------
P232 = 7.0
H232 = (0.01, 0.005, 0.0025)
SPEC232 = {"A0": dict(a0=1e-3), "A1R": dict(a1R=1e-3), "A1L": dict(a1L=1e-3), "L": dict(L=(1e-3, 3.5, 0.5)), "R": dict(R=(1e-3, 3.5, 0.5))}
PJ232 = (0.5, 1.0, 2.0, 4.0, 6.0)


def c2_slice(f, j, ih, h):
    """d_V of f on the horizon node (j, ih) from the data of the slice through it and its normal derivative (hs.slice_V_derivative); the mutation keeps the tangential half only."""
    if MUT == "fn_c1_aT":
        return hs.fd4(f, j, ih, 1, -1, h)
    return hs.slice_V_derivative(f, j, ih, h)


def c1_functional(dr_h, c2s, V0):
    """c1 = dr(V0, 0) - V0 c2 with c2 the slice functional (the horizon law dr_V V = c1 + c2 V extended to S)."""
    if MUT == "fn_c1_noV":
        return dr_h - c2s
    if MUT == "fn_c1_sign":
        return dr_h + V0 * c2s
    return dr_h - V0 * c2s


def run232(name, h):
    key = ("232", name, h, MUT if MUT in ("fn_c1_noV", "fn_c1_aT", "fn_c1_sign") else None)
    if key not in _RUNS:
        b = lin.ep.Member(C).kappa * SQ2
        su = lin.mksetup(C, h, P232, SPEC232[name])
        U = lin.march_full(su.UL, su.UR, su.chart.par, su.chart.tab, su.hp, su.hs)
        ih = su.jh
        p = h * np.arange(su.Np + 1)
        c1, c2, res = hs.horizon_fit(U[0, :, ih], np.exp(b * p))
        sl = {}
        for pj in PJ232:
            j = int(round(pj / h))
            V0 = math.exp(b * pj)
            c2s = c2_slice(U[0], j, ih, h) / (b * V0)
            sl[pj] = (c1_functional(U[0][j, ih], c2s, V0), c2s)
        _RUNS[key] = dict(c1=c1, c2=c2, res=res, sl=sl)
        del U
    return _RUNS[key]


# ----------------------------------------------------------------------------------------------------------------------
# bookkeeping
# ----------------------------------------------------------------------------------------------------------------------
RESULTS = []


def want(cid):
    return ONLY is None or cid in ONLY


def check(cid, text, ok, detail=""):
    if re.search(r"(?<![a-z])(nan|-?inf)(?![a-z])", str(detail).lower()):   # max() and <= hide a NaN
        ok = False
    RESULTS.append((cid, bool(ok)))
    print(f"[{'pass' if ok else 'FAIL'}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def sci(x):
    return f"{x:.2e}"


# ----------------------------------------------------------------------------------------------------------------------
# derive: sympy
# ----------------------------------------------------------------------------------------------------------------------
def derive_S1():
    import sympy as sp
    r = sp.symbols("rho", real=True)
    Af = sp.Function("A")(r)
    u = sp.Function("u")(r)
    R = sp.sqrt(1 + r ** 2)
    App = 2 * (Af - 1) / R ** 2                                              # the static equation R^2 A'' = 2 (A - 1)
    V = -App / 2 - r * sp.diff(Af, r) / R ** 2                              # (v rho): V = -A''/2 - rho A'/R^2
    if MUT == "s_pot":
        V = V + 1 / R ** 2
    q = (Af / 2) * (V + (1 - Af) / R ** 2) * (-1 if MUT == "s_q" else 1)     # q = (A/2)(V + (1 - A)/R^2)
    p = -sp.diff(Af, r) / 2 + Af * r / R ** 2
    # d/dx = A d/drho: A^2 u'' + A A' u' + p A u' + q u = 0
    upp = -((Af * sp.diff(Af, r) + p * Af) * sp.diff(u, r) + q * u) / Af ** 2
    Th = R * (Af * sp.diff(u, r) - sp.diff(Af, r) * u / 2)                  # kappa Theta
    dTh = sp.diff(Th, r).subs(sp.diff(u, r, 2), upp).subs(sp.diff(Af, r, 2), App)
    res = sp.simplify(dTh)
    check("S1", "x-form of H at dphi = dF = 0 with q = (A/2)(V + (1 - A)/R^2), V from (v rho), R^2 A'' = 2 (A - 1): d/drho[R (A dr' - A' dr/2)] = 0 (symbolic)",
          res == 0, f"residual {res}")


def derive_S2():
    import sympy as sp
    t = sp.symbols("t", positive=True)
    rho = 1 / t
    atn = sp.pi / 2 - sp.atan(t)                                             # arctan(1/t)
    K = 2 / sp.pi
    A = 1 - (1 + rho ** 2) + K * (rho + (1 + rho ** 2) * atn)
    r = sp.symbols("r", positive=True)
    Ar = 1 - (1 + r ** 2) + K * (r + (1 + r ** 2) * sp.atan(r))
    Apr = sp.diff(Ar, r).subs(r, rho).subs(sp.atan(1 / t), atn)
    cb = (1 - A * rho ** 2 / (1 + rho ** 2)) / 2 - rho * Apr / 2
    ser = sp.expand(sp.series(sp.simplify(cb), t, 0, 4).removeO())
    ser2 = sp.expand(sp.series(sp.sqrt(1 + rho ** 2) - rho, t, 0, 4).removeO())
    co = [sp.simplify(ser.coeff(t, k)) for k in range(3)]
    co2 = [sp.simplify(ser2.coeff(t, k)) for k in range(2)]
    ok = co[0] == 0 and co[1] == 0 and co[2] == sp.Rational(1, 2) and co2[0] == 0 and co2[1] == sp.Rational(1, 2)
    check("S2", "c_b = (1 - A rho^2/R^2)/2 - rho A'/2 = 1/(2 rho^2) + O(rho^-3) and R - rho = 1/(2 rho) + O(rho^-3) at rho -> infinity (c = -1): dm + kappa Theta = c_b dr + (R - rho) dr_x -> 0 for dr = O(ln rho)",
          bool(ok), f"c_b = {ser}; R - rho = {ser2}")


def derive_S3():
    import sympy as sp
    r, c = sp.symbols("rho c", real=True)
    R2 = 1 + r ** 2
    K = -2 * c / sp.pi if MUT != "s_pot" else -2 * c / sp.pi + c ** 2
    A = 1 + c * R2 + K * (r + R2 * sp.atan(r))
    V = -sp.diff(A, r, 2) / 2 - r * sp.diff(A, r) / R2
    dc = sp.simplify(sp.diff(sp.simplify(V / c), c))
    check("S3", "the potential V = -A''/2 - rho A'/R^2 of the member c on the profile phi = sqrt2 arctan rho is proportional to c (dV/c has no c): a neighbouring member solves a different theory",
          dc == 0, f"d(V/c)/dc = {dc}")


def derive_S4():
    import sympy as sp
    s, kap = sp.symbols("s kappa", positive=True)
    x = sp.symbols("rho", real=True)
    Af = sp.Function("A")(x)
    R = sp.sqrt(1 + x ** 2)
    ds = lambda f: (Af / (kap * s)) * sp.diff(f, x)              # d/ds = rho_s d/drho, rho_s = A/(kappa s) (Pi = s^2, d rho/d Pi = A/(2 kappa Pi))
    rho_Pi = Af / (2 * kap * s ** 2)
    mu = s * R
    mu_s_rPi = sp.simplify(mu * s * (x / R) * rho_Pi)             # mu s r_Pi = rho A/(2 kappa)
    E = Af / (2 * kap ** 2 * s ** 2)
    V = (1 - sp.diff(x * Af, x)) / R ** 2 + (1 / R ** 2 if MUT == "s_pot" else 0)   # (theta theta): V R^2 = 1 - (rho A)'
    lhs = ds(mu_s_rPi)
    rhs = mu * E * (1 / R - V * R)
    res = sp.simplify(lhs - rhs)
    check("S4", "the background identity of T1: mu s r_Pi = rho A/(2 kappa) and (mu s r_Pi)' = mu E (1/r - V r) with V R^2 = 1 - (rho A)' ((theta theta), #269), symbolic in A(rho)",
          res == 0 and sp.simplify(mu_s_rPi - x * Af / (2 * kap)) == 0, f"residual {res}")


def run_derive():
    t0 = time.time()
    for fn in (derive_S1, derive_S2, derive_S3, derive_S4):
        fn()
    bad = [c for c, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS) - len(bad)} of {len(RESULTS)} checks pass" + (f"; fired: {', '.join(bad)}" if bad else "") + f"  [{time.time() - t0:.0f} s]  mutation: {MUT}")
    return 1 if bad else 0


# ----------------------------------------------------------------------------------------------------------------------
# numerical checks
# ----------------------------------------------------------------------------------------------------------------------
def check_T1():
    """background identity (mu s r_Pi)' = mu E (1/r - V r) (what makes dF's coefficient in Theta_gen)."""
    worst = 0.0
    det = []
    for c in (-1.0, -1.5):
        bg = fs.bgc(c)
        w = 0.0
        for s in (0.05, 0.3, 0.7, 1.0, 1.5, 3.0, 10.0, 100.0):
            def B(x):
                b = bg.at(x)
                return x * x * b["r"] * b["rp"]
            hh = 1e-4 * s
            dB = (-B(s + 2 * hh) + 8 * B(s + hh) - 8 * B(s - hh) + B(s - 2 * hh)) / (12 * hh)
            b = bg.at(s)
            rhs = s * b["r"] * b["E"] * (1.0 / b["r"] - b["V0"] * b["r"])
            w = max(w, abs(dB - rhs) / abs(rhs))
        worst = max(worst, w)
        det.append(f"c = {c}: {w:.1e}")
    check("T1", "background identity (mu s r_Pi)' = mu E (1/r - V r), relative, 8 points of 0.05 <= s <= 100", worst <= 1e-8, "; ".join(det))


def check_T2():
    """the rho-form of H against the chart's s-form, the closed-form tail against the ODE, Theta conserved in rho."""
    bg = fs.bgc(C)
    sl = slice_obj()
    wp = wq = 0.0
    for s in (0.1, 0.5, 1.0, 2.0, 5.0, 30.0, 1e3, 1e5):
        P, Q, b = sl.PQ(s)
        rho = b["rho"]
        p, q = xform_pq(rho)
        wp = max(wp, abs((KAPPA * s * P - KAPPA) - p) / abs(p))
        wq = max(wq, abs(KAPPA ** 2 * s * s * Q - q) / abs(q))
    h = sl.hom_even(fs.SL, 1e5, (1.0, 0.0))
    wth = wdm = 0.0
    for s in (1.05, 3.0, 30.0, 1e3, 1e4, 1e5):
        dr, drs = h.sol(s)
        rho = bg.at(s)["rho"]
        thx = theta_x(rho, dr, fs.KAPPA * s * drs)
        ths = sl.theta(s, dr, drs)
        wth = max(wth, abs(thx - ths) / abs(ths))
        wdm = max(wdm, abs(dm_x(rho, dr, fs.KAPPA * s * drs) - fs.dm_of(bg, s, dr, drs)) / abs(fs.dm_of(bg, s, dr, drs)))
    s0 = 1.05
    dr0, drs0 = h.sol(s0)
    rho0 = bg.at(s0)["rho"]
    drx0 = fs.KAPPA * s0 * drs0
    sol = tail_ode(rho0, dr0, drx0)
    th0 = theta_x(rho0, dr0, drx0)
    wcl = wcon = 0.0
    for rho in (3.0, 16.0, 1e2, 1e4, 1e6, RHO_MAX):
        a_ = sol.sol(math.log(rho))
        b_ = tail_closed(rho0, dr0, drx0, rho)
        wcl = max(wcl, abs(a_[0] - b_[0]) / abs(b_[0]), abs(a_[1] - b_[1]) / abs(b_[1]))
        wcon = max(wcon, abs(theta_x(rho, *a_) - th0) / abs(th0))
    ok = max(wp, wq, wth, wdm) <= 1e-9 and wcl <= 1e-9 and wcon <= 1e-9
    check("T2", "rho-form p, q against the chart's P, Q (8 points); Theta and dm in rho against the s-form on a homogeneous solution (6 points); the closed-form tail against the ODE; Theta conserved on [rho(1.05), 1e8]",
          ok, f"p {wp:.1e}, q {wq:.1e}, Theta {wth:.1e}, dm {wdm:.1e}; closed form {wcl:.1e}; conservation {wcon:.1e}")


def check_T3():
    """Theta_gen = 0 on exact pure gauge with dF != 0 (c1 = 0 there), both signs of s."""
    worst = 0.0
    wold = 0.0
    det = []
    for c in (-1.0, -1.5):
        bg = fs.bgc(c)
        sl = fs.Slice(bg)
        fn = fs.fd_slice_fn(fs.gauge_fun(bg))
        w = 0.0
        for s in (-3.0, -1.5, -0.3, 0.1, 0.5, 1.0, 3.0, 10.0):
            U0, UT, UN, _ = fn(s)
            w = max(w, abs(theta_gen(sl, s, U0[0], UT[0], U0[2], U0[1])) / theta_scale(sl, s, U0[0], UT[0], U0[2], U0[1]))
            mu, a, cc, Q, b = sl.theta_coef(s)
            wold = max(wold, abs(a * U0[0] + mu * UT[0] + cc * U0[2]) / abs(U0[0]))
        worst = max(worst, w)
        det.append(f"c = {c}: {w:.1e}")
    check("T3", "Theta_gen vanishes on exact pure gauge L_xi(member), dF != 0, s = -3 ... 10 (finite-difference jets, step 1e-3), |Theta| over the largest of its four terms", worst <= 1e-8, "; ".join(det) + f"; Theta without the dF term: {wold:.1e}")


S_MODE = (0.05, 0.1, 0.3, 0.7, 1.0, 2.0, 5.0, 20.0, 1e2, 1e3, 1e4, 3e4, 1e5, 1e6)


def check_T4():
    """the exact static mode: c1, Theta_gen, dm, dM_ADM, the cut extensions."""
    m = mode()
    bg = fs.bgc(C)
    sl = slice_obj()
    c1 = m.fields(0.0)[0]
    wth = wdm = 0.0
    for s in S_MODE:
        dr, dr_P, dphi, dF = m.fields(s * s)
        drs = 2 * s * dr_P
        wth = max(wth, abs(theta_gen(sl, s, dr, drs, dphi, dF) + c1) / abs(c1))
        rho, rho_P, q, q_P = m.rq(s * s)
        dx = m.dm_exact(rho, q)
        wdm = max(wdm, abs(dm_s(bg, s, dr, drs, dF) - dx) / abs(dx))
    rho_e = m.rho_e
    q_big = m.big.sol(math.log(RHO_MAX))[0]
    dM_big = m.dm_exact(RHO_MAX, q_big)
    dM_mp = exact_dM_mp()
    cut_err = []
    for sc in (0.5, 2.0, 20.0, 1e3, -0.5, -2.0, -20.0):
        dr, dr_P, dphi, dF = m.fields(sc * sc)
        rm = read_mass(bg, sc, dr, 2 * sc * dr_P, dphi, dF)
        cut_err.append(abs(rm["dM"] - KAPPA_C * c1) / (KAPPA_C * c1))
    e_dM = abs(dM_big - dM_mp) / dM_mp
    e_mp = abs(dM_mp - 1 / math.pi) * math.pi
    e_kc = abs(KAPPA_C * c1 - dM_mp) / dM_mp
    ok = c1 == 0.5 or abs(c1 - 0.5) <= 1e-12
    ok = ok and wth <= 1e-9 and wdm <= 1e-9 and e_dM <= 1e-6 and e_mp <= 1e-9 and e_kc <= 1e-6 and max(cut_err) <= 1e-6
    check("T4", "exact static mode r = sqrt(L) R (d/dL at L = 1): c1 = 1/2; Theta_gen = -c1 at 14 values 0.05 <= s <= 1e6 (relative); dm of the slice code against d_L m at fixed Pi; dM_ADM from the family at rho = 1e8 against mpmath's d/dL of lim m_L = 1/pi; kappa c1 = dM_ADM; seven cut extensions (four right, three left) read at rho = 1e8",
          ok, f"c1 {c1:.12f}; Theta {wth:.1e}; dm {wdm:.1e}; dM(1e8) {dM_big:.12f} vs mp {dM_mp:.12f} ({e_dM:.1e}), mp vs 1/pi {e_mp:.1e}; kappa c1 vs dM {e_kc:.1e}; cut extensions {max(cut_err):.1e}")


def check_T5():
    """#258's compact datum G through Theta (Theta = 0), and a homogeneous start with mass through the rho-form (the 0.2% of #258 closed)."""
    sl = fs.Slice(fs.bgc(C))
    d = fs.datum("G")
    w = 0.0
    for s in np.linspace(-1.05, 1.05, 22):
        dr, drs, n = d.ev(s)
        w = max(w, abs(theta_gen(sl, s, dr, drs, d.dphi.f(s), 0.0)))
    # delta M = 1e-4 and 2e-4 left starts: c1 against dm at rho = 1e8
    bg = fs.bgc(C)
    errs = []
    cs = []
    for k in (1.0, 2.0):
        w0 = k * 1e-4 * fs.Datum.mass_start(sl)
        sol = sl.even(fs.SL, fs.SR, fs.Data([]), tuple(w0), max_step=0.01)
        c1 = sol.sol(0.0)[0]
        dr, drs = sol.sol(fs.SR)
        rm = read_mass(bg, fs.SR, dr, drs, 0.0, 0.0)
        errs.append(abs(KAPPA_C * c1 - rm["dM"]) / abs(rm["dM"]))
        cs.append((c1, rm["dM"], rm["dms"][1], rm["dms"][2]))
    ok = w <= 1e-12 and max(errs) <= 1e-6
    check("T5", "compact datum G (#258, ODE path): |Theta_gen| <= 1e-12 at 22 points of |s| <= 1.05; homogeneous left starts of dM = 1e-4 and 2e-4: kappa c1 against dm at rho = 1e8 (relative)",
          ok, f"max |Theta| {w:.1e}; kappa c1 vs dM: {errs[0]:.1e}, {errs[1]:.1e}; c1 {cs[0][0]:.6e}, {cs[1][0]:.6e}; dm(1e8) {cs[0][1]:.6e}, {cs[1][1]:.6e}; dm(1e2), dm(1e4) of the first: {cs[0][2]:.6e}, {cs[0][3]:.6e}")


def class_summary(name, hl):
    """Richardson of the two finest steps: c1 fit, c1 at S, Theta_gen(s), dM of the cuts; also the sequence for the ratios."""
    r = {h: class_run(name, h) for h in hl}
    h1, h2 = hl[-2], hl[-1]
    out = dict(c1=rich(r[h1]["c1fit"], r[h2]["c1fit"]), c1S=rich(r[h1]["c1S"], r[h2]["c1S"]), c2=rich(r[h1]["c2fit"], r[h2]["c2fit"]))
    out["theta"] = {s: rich(r[h1]["theta"][s], r[h2]["theta"][s]) for s in S_LIST}
    out["dM"] = {sc: rich(r[h1]["cut"][sc]["dM"], r[h2]["cut"][sc]["dM"]) for sc in S_CUTS + S_CUTL}
    out["thseq"] = {s: [r[h]["theta"][s] for h in hl] for s in S_LIST}
    out["c1seq"] = [r[h]["c1fit"] for h in hl]
    out["thraw"] = {s: [r[h]["theta"][s] for h in hl] for s in S_LIST}
    return out


QUOTED = {"A0": "1.658047e-03", "A1L": "2.094707e-03", "A1R": "5.236768e-05"}


def check_T6():
    """#248's classes A0, A1R, A1L (c1 != 0): Theta_gen on the evolved slice = -c1, and the cut extensions give dM_ADM = kappa c1."""
    md = open(os.path.join(HERE, "horizon_sector.md")).read()
    worst_th, worst_m, worst_q = 0.0, 0.0, 0.0
    det = []
    for name in ("A0", "A1R", "A1L"):
        sm = class_summary(name, HL)
        c1 = sm["c1"]
        e_th = max(abs(sm["theta"][s] + c1) / abs(c1) for s in S_LIST)
        e_S = abs(sm["c1S"] - c1) / abs(c1)
        e_m = max(abs(sm["dM"][sc] - KAPPA_C * c1) / (KAPPA_C * abs(c1)) for sc in S_CUTS + S_CUTL)
        e_q = abs(abs(c1) - abs(float(QUOTED[name]))) / abs(c1)
        r1, r2 = class_run(name, HL[-2]), class_run(name, HL[-1])
        e_f = max(abs(rich(r1[f"c1slice{V0}"], r2[f"c1slice{V0}"]) - c1) / abs(c1) for V0 in (0.5, 1.0))
        quoted_ok = QUOTED[name] in md
        worst_th, worst_m, worst_q = max(worst_th, e_th, e_S), max(worst_m, e_m), max(worst_q, e_q) if quoted_ok else 1.0
        det.append(f"{name}: c1 {c1:.7e}, Theta {e_th:.1e}, dr(S) {e_S:.1e}, (info) functional at V0 = 0.5, 1 {e_f:.1e}, dM {e_m:.1e}, vs #248's table {e_q:.1e}")
    check("T6", "evolved classes A0, A1R, A1L (Richardson of h = 0.0025, 0.00125): Theta_gen(s) = -c1 at s = -0.03 ... 1.9 and dr(S) = c1 to 1e-6 (relative; the slice functional at V0 = 0.5, 1 is printed, not a criterion: its error is V0 times that of c2, 20 V0 times c1's for A1R); dM_ADM of the cut extensions (s_cut = 1.0, 1.9 on the right and -0.03 on the left, read at rho = 1e8) = kappa c1 to 1e-6; c1 against #248's table to 1e-6",
          worst_th <= 1e-6 and worst_m <= 1e-6 and worst_q <= 1e-6, "; ".join(det))


def check_T7():
    """controls with c1 = 0: evolved L, G, R (Theta_gen -> 0 at the truncation rate) ."""
    det, ok = [], True
    for name in ("L", "G", "R"):
        sm = class_summary(name, HL)
        w = max(abs(sm["theta"][s]) for s in S_LIST)
        raw = [max(abs(sm["thraw"][s][k]) for s in S_LIST) for k in range(len(HL))]
        rt = ratios(raw)
        good = w <= 1e-12 and all(3.0 <= x <= 5.0 for x in rt[-2:])
        ok = ok and good
        det.append(f"{name}: |Theta|max {w:.1e} (Richardson); raw {', '.join(f'{x:.1e}' for x in raw)}; ratios {', '.join(f'{x:.2f}' for x in rt)}")
    check("T7", "evolved controls L, G, R (c1 = 0 in the continuum): max |Theta_gen(s)| <= 1e-12 after Richardson (the issue's bar) and the raw values fall by 4 +- 1 per halving at the last two steps (T5 gives the compact datum)",
          ok, "; ".join(det))


def check_T8():
    """the exact mode through the evolution: pointwise error, c1, c2, Theta_gen, the mass reading."""
    amp = 2e-3
    rs = {h: mode_run(h, amp) for h in HL}
    errs = [rs[h]["err"] / rs[h]["scale"] for h in HL]
    rt = ratios(errs)
    h1, h2 = HL[-2], HL[-1]
    c1_exact = amp * 0.5
    c1 = rich(rs[h1]["c1fit"], rs[h2]["c1fit"])
    c2 = rich(rs[h1]["c2fit"], rs[h2]["c2fit"])
    e_c1 = abs(c1 - c1_exact) / c1_exact
    e_S = abs(rich(rs[h1]["c1S"], rs[h2]["c1S"]) - c1_exact) / c1_exact
    e_th = max(abs(rich(rs[h1]["theta"][s], rs[h2]["theta"][s]) + c1_exact) / c1_exact for s in S_LIST)
    dM_exact = amp / math.pi
    e_m = max(abs(rich(rs[h1]["cut"][sc]["dM"], rs[h2]["cut"][sc]["dM"]) - dM_exact) / dM_exact for sc in S_CUTS + S_CUTL)
    e_f1 = max(abs(rich(rs[h1][f"c1slice{V0}"], rs[h2][f"c1slice{V0}"]) - c1_exact) / c1_exact for V0 in (0.5, 1.0))
    ok = errs[-1] <= 1e-3 and abs(rt[-1] - 4.0) <= 0.3 and e_c1 <= 1e-6 and e_S <= 1e-6 and e_th <= 1e-6 and e_m <= 1e-6 and e_f1 <= 1e-6 and abs(c2) <= 1e-6 * c1_exact
    check("T8", "exact static mode (amplitude dL = 2e-3, c1 = 1e-3, dM_ADM = dL/pi) as ray data through the Kruskal solver: error/scale at h = 0.00125, last ratio; c1 of the fit, dr(S), the slice functional c1 = dr - V0 c2 at V0 = 0.5, 1; Theta_gen(s); dM_ADM of the cut extensions; |c2|",
          ok, f"errors {', '.join(sci(e) for e in errs)}; last ratio {rt[-1]:.2f}; c1 fit {e_c1:.1e}, dr(S) {e_S:.1e}, functional {e_f1:.1e}, Theta {e_th:.1e}, dM {e_m:.1e}; c2 {c2:.1e}")


def check_T9():
    """#232's own linear solutions: c1 as a functional of the slice data."""
    det, ok = [], True
    for name in ("A0", "A1R", "A1L", "L", "R"):
        r = {h: run232(name, h) for h in H232}
        h1, h2 = H232[-2], H232[-1]
        c1f = rich(r[h1]["c1"], r[h2]["c1"])
        c2f = rich(r[h1]["c2"], r[h2]["c2"])
        devs = []
        for pj in PJ232:
            c1s = rich(r[h1]["sl"][pj][0], r[h2]["sl"][pj][0])
            devs.append(abs(c1s - c1f))
        if name in ("A0", "A1R", "A1L"):
            rel = [d / abs(c1f) for d in devs]
            good = max(rel) <= 1e-6
            det.append(f"{name}: c1 {c1f:.7e}, c2 {c2f:.5e}, functional vs fit at p = 0.5 ... 6: " + ", ".join(sci(x) for x in rel))
        else:
            good = max(abs(rich(r[h1]["sl"][pj][0], r[h2]["sl"][pj][0])) for pj in PJ232) <= 1e-9 and abs(c1f) <= 1e-9
            det.append(f"{name}: c1 fit {c1f:.1e}, functional " + ", ".join(f"{rich(r[h1]['sl'][pj][0], r[h2]['sl'][pj][0]):.1e}" for pj in PJ232))
        ok = ok and good
    check("T9", "#232's solver (linear.py, c = -1, h = 0.005, 0.0025 Richardson, P = 7): c1 = dr(p0) - e^{b p0} c2_slice from the slices p + sigma = const crossing the horizon at p0 = 0.5, 1, 2, 4, 6 equals the horizon fit to 1e-6 (relative to |c1|; classes A0, A1R, A1L), and |c1| <= 1e-9 for the pulse classes L, R",
          ok, "; ".join(det))


def check_T10():
    """the first law dM = (kappa_xi/8 pi) dA along the exact static family at finite L (mpmath, 40 digits): dA = 8 pi r_h dr_h is the area change that c1 measures."""
    import mpmath as mp
    mp.mp.dps = 40
    c, K = mp.mpf(C), 2 / mp.pi
    A = lambda r: 1 + c * (r * r + 1) + K * (r + (r * r + 1) * mp.atan(r))
    Ap = lambda r: 2 * c * r + 2 * K * (1 + r * mp.atan(r))

    def rh(L):
        return mp.findroot(lambda r: A(r) - (1 - 1 / L), (mp.pi / 4) * (1 - 1 / L))

    area = lambda L: 4 * mp.pi * L * (1 + rh(L) ** 2)
    kap = lambda L: (1 if MUT == "law_norm" else mp.sqrt(L)) * Ap(rh(L)) / 2           # xi = sqrt(L) d_t is unit at infinity (F_L -> 1/L)
    M = lambda L: 2 * L ** mp.mpf(1.5) / (3 * mp.pi)
    worst, rows = 0.0, []
    for L in ("0.5", "0.9", "1", "1.5", "2", "4"):
        L = mp.mpf(L)
        ratio = kap(L) * mp.diff(area, L) / (8 * mp.pi) / mp.diff(M, L)
        worst = max(worst, abs(float(ratio) - 1.0))
        rows.append(f"{float(L)}: {float(ratio):.12f}")
    check("T10", "first law along the exact static family r = sqrt(L) R (finite L, mpmath): dM/dL = kappa_xi (dA/dL)/(8 pi), A = 4 pi L R(rho_h)^2, kappa_xi = sqrt(L) A'(rho_h)/2, at L = 0.5 ... 4 (ratio of the two sides)",
          worst <= 1e-12, "ratios " + ", ".join(rows))


CHECKS = [("T1", check_T1), ("T2", check_T2), ("T3", check_T3), ("T4", check_T4), ("T5", check_T5), ("T6", check_T6), ("T7", check_T7), ("T8", check_T8), ("T9", check_T9), ("T10", check_T10)]
SCOPE = tuple(c for c, _ in CHECKS)
SCOPE_D = ("S1", "S2", "S3", "S4")


def run_checks():
    t0 = time.time()
    for key, fnc in CHECKS:
        if ONLY is not None and key not in ONLY:
            continue
        t1 = time.time()
        fnc()
        print(f"        ({key} {time.time() - t1:.0f} s)", flush=True)
    bad = [c for c, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS) - len(bad)} of {len(RESULTS)} checks pass" + (f"; fired: {', '.join(bad)}" if bad else "") + f"  [{time.time() - t0:.0f} s]  mutation: {MUT}")
    return 1 if bad else 0


# ----------------------------------------------------------------------------------------------------------------------
# declared fired sets and the comparison with stored logs
# ----------------------------------------------------------------------------------------------------------------------
DECLARED = {                      # computed in the development runs d4 and d5 (raw I259-mac/dev/logs) before the runs of record; every check is in scope
    "s_q": {"S1"},
    "s_pot": {"S1", "S3", "S4"},
    "th_noF": {"T3", "T4", "T6", "T7", "T8"},
    "th_nophi": {"T3", "T4", "T5", "T6", "T7", "T8"},
    "th_asign": {"T3", "T4", "T5", "T6", "T7", "T8"},
    "dm_noT": {"T2", "T4", "T5", "T6", "T8"},
    "kappa_off": {"T4", "T5", "T6"},
    "jump_off": {"T4", "T6", "T8"},
    "xform_p": {"T2", "T4", "T5", "T6", "T8"},
    "xform_q": {"T2", "T4", "T5", "T6", "T8"},
    "mode_khat": {"T4", "T8"},
    "mode_dF0": {"T4", "T8"},
    "fn_c1_noV": {"T9"},
    "fn_c1_aT": {"T8", "T9"},
    "fn_c1_sign": {"T9"},
    "law_norm": {"T10"},
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
        scope = set(SCOPE_D if name in KNOWN_D else SCOPE)
        inscope = fired & scope
        ok = inscope == DECLARED.get(name, set())
        bad += 0 if ok else 1
        print(f"{name}: fired {sorted(inscope)}; declared {sorted(DECLARED.get(name, set()))}; {'equal' if ok else 'DIFFERENT'}")
    print(f"{len(KNOWN) - bad} of {len(KNOWN)} fired sets equal the declared ones")
    return 1 if bad else 0


def report_table():
    bg = fs.bgc(C)
    sl = slice_obj()
    m = mode()
    c1 = m.fields(0.0)[0]
    print("== exact static mode r = sqrt(L) R, F_L = A + 1/L - 1 (#269 N6), d/dL at L = 1, fixed Pi ==")
    print(f"c1 = dr(S) = {c1:.12f}; c2 = 0; dM_ADM = d_L [2 L^(3/2)/(3 pi)] = 1/pi = {1 / math.pi:.12f}; kappa c1 = {KAPPA * c1:.12f}; mpmath d_L lim m_L = {exact_dM_mp():.12f}")
    print(f"{'s':>9} {'rho':>9} {'Theta_gen + c1':>15} {'Theta (no dF)':>15} {'dm slice':>15} {'dm exact':>15}")
    for s in S_MODE:
        dr, dr_P, dphi, dF = m.fields(s * s)
        rho, rho_P, q, q_P = m.rq(s * s)
        drs = 2 * s * dr_P
        mu, a, cc, Q, b = sl.theta_coef(s)
        print(f"{s:9.3g} {rho:9.4f} {theta_gen(sl, s, dr, drs, dphi, dF) + c1:15.3e} {a * dr + mu * drs + cc * dphi:15.6e} {dm_s(bg, s, dr, drs, dF):15.10f} {m.dm_exact(rho, q):15.10f}")
    print("dm(rho) at fixed Pi, rho > 21.6, minus 1/pi, times rho (the 1/rho convergence of the mass at fixed Pi):")
    for rho in (30.0, 1e2, 1e3, 1e4, 1e6, 1e8):
        q = m.big.sol(math.log(rho))[0]
        d = m.dm_exact(rho, q)
        print(f"   rho = {rho:8.1e}: dm = {d:.12f}, (dm - 1/pi) rho = {(d - 1 / math.pi) * rho:.6f}")
    print("== #248's classes (corner (-0.05, -2) of the rectangle, amplitude 1e-3), Richardson of h = 0.0025, 0.00125 ==")
    print(f"{'class':>4} {'c1 fit':>14} {'dr(S)':>14} {'c2':>14} {'kappa c1':>14} {'dM cut 1.0':>14} {'dM cut 1.9':>14} {'dM cut -0.03':>14} {'Theta(0.5)+c1':>14} {'Theta(1.9)+c1':>14}")
    for name in ("A0", "A1R", "A1L", "G", "L", "R"):
        sm = class_summary(name, HL)
        print(f"{name:>4} {sm['c1']:14.7e} {sm['c1S']:14.7e} {sm['c2']:14.7e} {KAPPA * sm['c1']:14.7e} {sm['dM'][1.0]:14.7e} {sm['dM'][1.9]:14.7e} {sm['dM'][-0.03]:14.7e} {sm['theta'][0.5] + sm['c1']:14.2e} {sm['theta'][1.9] + sm['c1']:14.2e}")
    print("c1 of the horizon fit and Theta_gen(s = 0, 0.5, 1.9) against h = 0.01, 0.005, 0.0025, 0.00125:")
    for name in ("A0", "A1L"):
        r = [class_run(name, h) for h in HL]
        print(f"   {name}: c1 fit " + ", ".join(f"{x['c1fit']:.7e}" for x in r) + "; Theta(0) " + ", ".join(f"{x['theta'][0.0]:.7e}" for x in r) + "; Theta(1.9) " + ", ".join(f"{x['theta'][1.9]:.7e}" for x in r))
    print("== the formal neighbouring member d/dc at c = -1 (a solution of the potential V_c, not of V_-1): mpmath, 40 digits ==")
    import mpmath as mp
    mp.mp.dps = 40

    def Ac_(cc, r):
        return 1 + cc * (r * r + 1) - 2 * cc / mp.pi * (r + (r * r + 1) * mp.atan(r))

    def M_of(cc):
        r = mp.mpf(10) ** 14
        return (1 + r * r * (1 - Ac_(cc, r))) / (2 * mp.sqrt(1 + r * r))

    def area_S(cc):
        rr = mp.findroot(lambda x: Ac_(cc, x), (mp.pi / 4) * (-(cc + 1)) if cc != -1 else mp.mpf(0))
        return mp.sqrt(1 + rr * rr)                    # r at S = R(rho_h(c))
    dMdc = mp.diff(M_of, mp.mpf(-1))
    c1c = mp.diff(area_S, mp.mpf(-1))
    print(f"dM_ADM/dc = {float(dMdc):.12f} (-2/(3 pi) = {-2 / (3 * math.pi):.12f}); c1 = d_c R(rho_h(c)) = {float(c1c):.2e} (R'(0) = 0); kappa c1 = 0 against dM_ADM/dc: the identity does not hold for it, and it is not a perturbation of the fixed theory (S3)")
    print("== the evolved exact mode (amplitude 2e-3), Richardson of h = 0.0025, 0.00125 ==")
    rs = {h: mode_run(h, 2e-3) for h in HL}
    h1, h2 = HL[-2], HL[-1]
    print(f"c1 fit {rich(rs[h1]['c1fit'], rs[h2]['c1fit']):.9e} (exact 1e-3); dr(S) {rich(rs[h1]['c1S'], rs[h2]['c1S']):.9e}; c2 {rich(rs[h1]['c2fit'], rs[h2]['c2fit']):.2e}; dM cuts "
          + ", ".join(f"{rich(rs[h1]['cut'][sc]['dM'], rs[h2]['cut'][sc]['dM']):.9e}" for sc in S_CUTS + S_CUTL) + f" (exact {2e-3 / math.pi:.9e}); error/scale " + ", ".join(sci(rs[h]["err"] / rs[h]["scale"]) for h in HL))
    print("== #258's ODE path: compact datum G and homogeneous left starts with mass ==")
    d = fs.datum("G")
    th = [theta_gen(sl, s, *d.ev(s)[:2], d.dphi.f(s), 0.0) for s in np.linspace(-1.05, 1.05, 22)]
    print(f"G: c1 = {d.c1:.2e}, c2 = {d.c2:.9e}; max |Theta_gen| over 22 points = {max(abs(t) for t in th):.2e}")
    for k in (1.0, 2.0):
        w0 = k * 1e-4 * fs.Datum.mass_start(sl)
        sol = sl.even(fs.SL, fs.SR, fs.Data([]), tuple(w0), max_step=0.01)
        c1h = sol.sol(0.0)[0]
        dr, drs = sol.sol(fs.SR)
        rm = read_mass(bg, fs.SR, dr, drs, 0.0, 0.0)
        print(f"homogeneous start (dM(3e4) = {k * 1e-4:.0e} of #258): c1 = {c1h:.9e}, kappa c1 = {KAPPA * c1h:.9e}; dm at rho = " + ", ".join(f"{r_:.0e}: {v:.9e}" for r_, v in zip(rm["rhos"][1:], rm["dms"][1:])) + f"; -kappa Theta = {-KAPPA * rm['th_after']:.9e}; dm(3e4) of the s-form = {fs.dm_of(bg, fs.SFAR, *sl.hom_even(fs.SR, fs.SFAR, (dr, drs)).sol(fs.SFAR)):.9e}")
    print("== #232's solver (linear.py, c = -1, P = 7), Richardson of h = 0.005, 0.0025: c1 as a functional of the slice data ==")
    for name in ("A0", "A1R", "A1L", "L", "R"):
        r = {h: run232(name, h) for h in H232}
        h1, h2 = H232[-2], H232[-1]
        c1f, c2f = rich(r[h1]["c1"], r[h2]["c1"]), rich(r[h1]["c2"], r[h2]["c2"])
        print(f"{name:>3}: c1 fit {c1f:.9e} (h = 0.01, 0.005, 0.0025: " + ", ".join(f"{r[h]['c1']:.7e}" for h in H232) + f"), c2 fit {c2f:.7e}; functional at p0 = 0.5, 1, 2, 4, 6: "
              + ", ".join(f"{rich(r[h1]['sl'][pj][0], r[h2]['sl'][pj][0]):.9e}" for pj in PJ232))
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
        sys.exit(report_mutations(sys.argv[sys.argv.index("mutations") + 1], sys.argv[sys.argv.index("mutations") + 2]))
    sys.exit(f"unknown command {cmd}")

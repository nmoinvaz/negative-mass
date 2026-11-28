#!/usr/bin/env python3
"""kappa r_S c1 = dM_ADM at the black universes c = -1.5, c = -0.3 and on Schwarzschild (#301), with c = -1 as the control of #259.

Members: the Ellis profile R = sqrt(1 + rho^2), phi = sqrt2 arctan rho, A = 1 + c R^2 + K (rho + R^2 arctan rho), K = -2c/pi (#212), and Schwarzschild (R = rho, A = 1 - 1/rho,
r_S = 2M = 1, no scalar). The perturbation is the static mass mode: for the Ellis members the exact family r = sqrt(L) R, F = A + 1/L - 1 (#269, N6; same potential), for
Schwarzschild the homothety r_S -> r_S (1 + eps), both written in the Kruskal chart of #248 at fixed Pi and differentiated at the member. #259's first integral Theta_gen
(a dr + mu dr' + c dphi + mu s r_Pi dF, conserved on the slice V + U = 0) is read at one point and gives c1 = dr(S) = -Theta_gen / r_S, because a(0) = -r_S; beyond a cut the
exterior in rho gives dM_ADM = -kappa Theta, so kappa r_S c1 = dM_ADM (Misner-Sharp mass at infinity, G = 1, kappa = A'(rho_h)/2 for the Killing field with A -> 1).

    static_c1_members.py derive [--mutate=NAME]               sympy checks (S-series)
    static_c1_members.py check [--only=T1,..] [--mutate=NAME] numerical checks (T-series); exit 1 when one fails
    static_c1_members.py table                                the numbers of the note
    static_c1_members.py mutations DIR PREFIX                 fired sets of the stored logs DIR/PREFIX<name>.log against the declared ones

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
sys.argv = _ARGV

MUT = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mutate=")), None)
ONLY = next((a.split("=", 1)[1].split(",") for a in sys.argv[1:] if a.startswith("--only=")), None)
KNOWN_D = ("s_q", "s_pot", "s_k")
KNOWN_N = ("rs_one", "ctl_true", "c1_noshift", "th_noF", "th_nophi", "th_asign", "dm_noT", "kappa_off", "jump_off", "xform_p", "xform_q", "rho_ellis",
           "mode_khat", "mode_dF0", "law_norm", "law_mass", "gauge_noRp")
KNOWN = KNOWN_D + KNOWN_N
assert MUT is None or MUT in KNOWN, MUT

SQ2 = math.sqrt(2.0)
EPS = fs.EPS
RHO_MAX = 1e8
S_CUTS = (1.0, 1.9)
S_CUTL = (-0.03,)
S_LIST = (-0.03, 0.0, 0.1, 0.25, 0.5, 1.0, 1.5, 1.9)
HL = (0.005, 0.0025, 0.00125, 0.000625)
RECT = hs.RECT
V0R, U0R = RECT[0], RECT[2]
BAR = 1e-6                                      # the issue's bar for kappa r_S c1 = dM_ADM


# ----------------------------------------------------------------------------------------------------------------------
# the members: A, A', A'', R, R' in double precision (series in 1/rho for rho >= 6), horizon, kappa, r_S
# ----------------------------------------------------------------------------------------------------------------------
class Ellis:
    """The member c of #212: A - 1 = c R^2 + K (rho + R^2 arctan rho) = K (rho - R^2 arccot rho) for rho > 0."""
    kind = "ellis"

    def __init__(self, c):
        self.c = float(c)
        self.K = -2.0 * self.c / math.pi
        self.tag = f"c={c:g}"
        self.rh = 0.0 if self.c == -1.0 else self._root()
        self.rS = math.sqrt(1.0 + self.rh ** 2)                       # r_S = R(rho_h), the radius of the bifurcation surface
        self.kappa = 0.5 * self.Ap(self.rh)                           # A'(rho_h)/2, the Killing field d_t with A -> 1
        self.khat = self.App(self.rh) / self.Ap(self.rh) ** 2         # kappa'/kappa along the family F_L = A + 1/L - 1
        self.Rph = self.rh / self.rS

    def _root(self):
        xs = np.linspace(-6.0, 6.0, 1201)
        vs = [self.A(x) for x in xs]
        for k in range(len(xs) - 1):
            if vs[k] < 0.0 <= vs[k + 1]:
                return brentq(self.A, xs[k], xs[k + 1], xtol=1e-15, rtol=1e-15)
        raise ValueError("no horizon")

    def Am1(self, rho):
        if rho >= 6.0:
            t = 1.0 / rho
            t2 = t * t
            s, pw = 0.0, t
            for j in range(45):
                s += (-1) ** j * pw * 2.0 / ((2 * j + 1) * (2 * j + 3))
                pw *= t2
            return -self.K * s
        return self.c * (rho * rho + 1.0) + self.K * (rho + (rho * rho + 1.0) * math.atan(rho))

    def A(self, rho):
        return 1.0 + self.Am1(rho)

    def Ap(self, rho):
        if rho >= 6.0:
            t = 1.0 / rho
            t2 = t * t
            s, pw = 0.0, t2
            for k in range(1, 45):
                s += (-1) ** (k + 1) * pw / (2 * k + 1)
                pw *= t2
            return 2.0 * self.K * s
        return 2.0 * self.c * rho + 2.0 * self.K * (1.0 + rho * math.atan(rho))

    def App(self, rho):
        if rho >= 6.0:
            t = 1.0 / rho
            t2 = t * t
            s, pw = 0.0, t2 * t
            for k in range(1, 45):
                s += (-1) ** k * 2.0 * k * pw / (2 * k + 1)
                pw *= t2
            return 2.0 * self.K * s
        return 2.0 * self.c + 2.0 * self.K * (math.atan(rho) + rho / (1.0 + rho * rho))

    def R(self, rho):
        return math.sqrt(1.0 + rho * rho)

    def Rp(self, rho):
        return rho / math.sqrt(1.0 + rho * rho)

    def omAR2(self, rho):
        """1 - A R'^2 = (1 - (A - 1) rho^2)/R^2, no cancellation."""
        return (1.0 - self.Am1(rho) * rho * rho) / (1.0 + rho * rho)


class Schw:
    """Schwarzschild with r_S = 2M = 1: R = rho, A = 1 - 1/rho, no scalar (the planted case of #248). rho_ellis puts the Ellis R into the exterior form."""
    kind = "schw"
    c = None
    K = None
    tag = "Schwarzschild"
    rh = 1.0
    rS = 1.0
    kappa = 0.5
    khat = None
    Rph = 1.0

    def Am1(self, rho):
        return -1.0 / rho

    def A(self, rho):
        return 1.0 - 1.0 / rho

    def Ap(self, rho):
        return 1.0 / rho ** 2

    def App(self, rho):
        return -2.0 / rho ** 3

    def R(self, rho):
        return math.sqrt(1.0 + rho * rho) if MUT == "rho_ellis" else rho

    def Rp(self, rho):
        return rho / math.sqrt(1.0 + rho * rho) if MUT == "rho_ellis" else 1.0

    def omAR2(self, rho):
        return 1.0 / rho


_MEM = {}
_BGS = {}
_SLS = {}


def member(tag):
    if tag not in _MEM:
        _MEM[tag] = Schw() if tag == "schw" else Ellis(float(tag))
    return _MEM[tag]


def bg_of(m):
    """The member in the Kruskal chart of #248 (Pi = 1 at rho = 1), tabulated to Pi = 1e12 (Schwarzschild: closed form, Pi <= 1e3)."""
    if m.tag not in _BGS:
        _BGS[m.tag] = fs.Bg(-1.0, True) if m.kind == "schw" else fs.Bg(m.c, False, pi_max=1e12)
    return _BGS[m.tag]


def slice_of(m):
    if m.tag not in _SLS:
        _SLS[m.tag] = fs.Slice(bg_of(m))
    return _SLS[m.tag]


def solver_chart(m):
    if m.kind == "schw":
        kc = hs.KChart(-1.0)
        kc.par[hs.P_KIND2] = 1.0
        return kc
    return hs.kchart(m.c)


def kap_c(m):
    """the kappa that multiplies c1 when it is compared with dM_ADM (the mutation scales it)."""
    return m.kappa * (1.01 if MUT == "kappa_off" else 1.0)


def rs_c(m):
    """the r_S of the identity kappa r_S c1 = dM_ADM (the mutation sets it to 1)."""
    return 1.0 if MUT == "rs_one" else m.rS


# ----------------------------------------------------------------------------------------------------------------------
# Theta_gen and the Misner-Sharp mass perturbation on the slice (s-form, chart of #248), with the mutations
# ----------------------------------------------------------------------------------------------------------------------
def theta_gen(sl, s, dr, drs, dphi, dF):
    """Theta_gen = a dr + mu dr' + c dphi + mu s r_Pi dF (#259). a(0) = -r_S, so c1 = dr(S) = -Theta_gen / r_S."""
    mu, a, c, Q, b = sl.theta_coef(s)
    t = (-1.0 if MUT == "th_asign" else 1.0) * a * dr + mu * drs
    if MUT != "th_nophi":
        t += c * dphi
    if MUT != "th_noF":
        t += mu * s * b["rp"] * dF
    return t


def theta_scale(sl, s, dr, drs, dphi, dF):
    mu, a, c, Q, b = sl.theta_coef(s)
    return max(abs(a * dr), abs(mu * drs), abs(c * dphi), abs(mu * s * b["rp"] * dF))


def dm_s(bg, s, dr, drs, dF):
    """delta m = dr (1/2 - e^F s^2 r'^2) - r e^F s r' (s r' dF + dr') (#258)."""
    b = bg.at(s)
    eF = 1.0 / b["E"]
    grp = 0.0 if MUT == "dm_noT" else 1.0
    return dr * (0.5 - eF * s * s * b["rp"] ** 2) - grp * b["r"] * eF * s * b["rp"] * (s * b["rp"] * dF + drs)


# ----------------------------------------------------------------------------------------------------------------------
# the exterior in rho (dphi = dF = 0): y_xx + p y_x + q y = 0, y = dr, d/dx = A d/drho; kappa Theta = R (y_x - A' y/2) is its first integral
# ----------------------------------------------------------------------------------------------------------------------
def xform_pq(m, rho):
    A, A1, A2, R, Rp = m.A(rho), m.Ap(rho), m.App(rho), m.R(rho), m.Rp(rho)
    p = -A1 / 2.0 + (0.0 if MUT == "xform_p" else A * Rp / R)
    q = (1.0 if MUT == "xform_q" else -1.0) * (A / 2.0) * (A2 + A1 * Rp / R)
    return p, q


def ktheta_x(m, rho, y, yx):
    """kappa Theta = R (y_x - A' y / 2)."""
    return m.R(rho) * (yx - m.Ap(rho) * y / 2.0)


def dm_x(m, rho, y, yx):
    """delta m at dphi = dF = 0: (1/2)(1 - A R'^2) y - R R' y_x."""
    grp = 0.0 if MUT == "dm_noT" else 1.0
    return 0.5 * m.omAR2(rho) * y - grp * m.R(rho) * m.Rp(rho) * yx


class Tail:
    """(y, y_x) on [rho0, rho1] from the ODE: in rho up to 4, in u = ln rho beyond."""
    RT = 4.0

    def __init__(self, m, rho0, y0, yx0, rho1=RHO_MAX):
        def f_rho(rho, Y):
            A = m.A(rho)
            p, q = xform_pq(m, rho)
            return [Y[1] / A, -(p * Y[1] + q * Y[0]) / A]

        def f_u(u, Y):
            rho = math.exp(u)
            A = m.A(rho)
            p, q = xform_pq(m, rho)
            return [rho * Y[1] / A, -rho * (p * Y[1] + q * Y[0]) / A]
        kw = dict(method="DOP853", rtol=1e-13, atol=1e-20, dense_output=True)
        self.s1 = None
        Y1, start = [y0, yx0], rho0
        if rho0 < self.RT:
            self.s1 = solve_ivp(f_rho, [rho0, self.RT], Y1, **kw)
            Y1, start = list(self.s1.y[:, -1]), self.RT
        self.start = start
        self.s2 = solve_ivp(f_u, [math.log(start), math.log(rho1)], Y1, **kw)

    def at(self, rho):
        if self.s1 is not None and rho <= self.RT:
            return self.s1.sol(rho)
        return self.s2.sol(math.log(rho))


def tail_closed(m, rho0, y0, yx0, rho):
    """y = sqrt(A) [y0/sqrt(A0) + kappa Theta int drho/(R A^{3/2})] and y_x = A' y/2 + kappa Theta/R (the first integral solved for y and y_x)."""
    kth = ktheta_x(m, rho0, y0, yx0)
    f = lambda r: 1.0 / (m.R(r) * m.A(r) ** 1.5)
    g = lambda u: math.exp(u) * f(math.exp(u))
    I = 0.0
    lo = rho0
    if rho0 < Tail.RT and rho > rho0:
        hi = min(rho, Tail.RT)
        I += quad(f, lo, hi, epsabs=1e-15, epsrel=1e-13, limit=400)[0]
        lo = hi
    if rho > lo:
        I += quad(g, math.log(lo), math.log(rho), epsabs=1e-15, epsrel=1e-13, limit=800)[0]
    y = math.sqrt(m.A(rho)) * (y0 / math.sqrt(m.A(rho0)) + kth * I)
    return y, m.Ap(rho) * y / 2.0 + kth / m.R(rho)


def read_mass(m, s_cut, dr, drs, dphi, dF):
    """Hard cut of the sources at s_cut: dphi = dF = 0 beyond, the jump of dr' from H (#259), the exterior ODE to rho = RHO_MAX. dM = dm(RHO_MAX)."""
    bg, sl = bg_of(m), slice_of(m)
    b = bg.at(s_cut)
    th_before = theta_gen(sl, s_cut, dr, drs, dphi, dF)
    jump = 0.0 if MUT == "jump_off" else EPS * b["r"] * s_cut * b["php"] * dphi + s_cut * b["rp"] * dF
    drs_after = drs + jump
    rho0 = b["rho"]
    yx0 = m.kappa * s_cut * drs_after
    tl = Tail(m, rho0, dr, yx0)
    rhos = [rho0, 1e2, 1e4, 1e6, RHO_MAX]
    dms = [dm_x(m, r_, *tl.at(r_)) for r_ in rhos]
    kth_after = ktheta_x(m, rho0, dr, yx0)
    return dict(th_before=th_before, kth_after=kth_after, rhos=rhos, dms=dms, dM=dms[-1], rho0=rho0, dr0=dr, yx0=yx0, tail=tl)


# ----------------------------------------------------------------------------------------------------------------------
# the exact static mass mode in the Kruskal chart at fixed Pi, d/dL at L = 1 (Ellis) or d/d r_S at r_S = 1 (Schwarzschild)
# ----------------------------------------------------------------------------------------------------------------------
def a_series(m, n):
    """Taylor coefficients of A about rho_h from R^2 A'' = 2 (A - 1): beta_0 = A(rho_h) - 1 = -1, beta_1 = A'(rho_h) = 2 kappa, a two-term recursion."""
    s0, rh = 1.0 + m.rh ** 2, m.rh
    beta = np.zeros(n + 3)
    beta[0], beta[1] = -1.0, 2.0 * m.kappa
    for k in range(0, n + 1):
        beta[k + 2] = (2.0 * beta[k] - k * (k - 1) * beta[k] - 2.0 * rh * k * (k + 1) * beta[k + 1]) / (s0 * (k + 2) * (k + 1))
    al = beta.copy()
    al[0] = 0.0
    return al


class EllisMode:
    """q = d rho/dL at fixed Pi solves A q' = A' q - 1 - khat A (d/drho), regular at rho_h with q(rho_h) = 1/A'(rho_h); the free constant is the scaling gauge xi = V d_V + U d_U (q_1 = 0).
    dr = R/2 + R' q, dphi = sqrt2 q/R^2, dF = khat - q', c1 = dr(S) = R_h/2 + R'_h/(2 kappa) per unit dL, dM_ADM = K/2."""
    NT = 40
    D0 = 0.2

    def __init__(self, m):
        self.m = m
        self.khat = m.khat * (0.0 if MUT == "mode_khat" else 1.0)
        al = a_series(m, self.NT + 1)
        q = np.zeros(self.NT + 1)
        q[0] = 1.0 / al[1]
        for n in range(2, self.NT + 1):
            s1 = sum(al[k] * (n - k + 1) * q[n - k + 1] for k in range(2, n + 1))
            s2 = sum(k * al[k] * q[n - k + 1] for k in range(2, n + 2))
            q[n] = -(s1 - s2 + self.khat * al[n]) / (al[1] * (n - 1))
        self.qs = q
        self.dqs = np.array([k * q[k] for k in range(1, self.NT + 1)])
        rhs = lambda rho, y: [(m.Ap(rho) * y[0] - 1.0 - self.khat * m.A(rho)) / m.A(rho)]
        kw = dict(method="DOP853", rtol=1e-13, atol=1e-18, dense_output=True)
        P = np.polynomial.polynomial.polyval
        self.out = solve_ivp(rhs, [m.rh + self.D0, 40.0], [P(self.D0, q)], **kw)
        self.inn = solve_ivp(rhs, [m.rh - self.D0, m.rh - 2.0], [P(-self.D0, q)], **kw)
        rhs_u = lambda u, y: [math.exp(u) * (m.Ap(math.exp(u)) * y[0] - 1.0 - self.khat * m.A(math.exp(u))) / m.A(math.exp(u))]
        self.big = solve_ivp(rhs_u, [math.log(40.0), math.log(RHO_MAX)], [self.out.y[0, -1]], **kw)

    def q_qp(self, rho):
        m = self.m
        d = rho - m.rh
        P = np.polynomial.polynomial.polyval
        if abs(d) <= self.D0:
            return P(d, self.qs), P(d, self.dqs)
        if d > 0:
            q = self.out.sol(rho)[0] if rho <= 40.0 else self.big.sol(math.log(rho))[0]
        else:
            q = self.inn.sol(rho)[0]
        return q, (m.Ap(rho) * q - 1.0 - self.khat * m.A(rho)) / m.A(rho)

    def fields(self, Pi, par, tab):
        """(dr, dr_Pi, dphi, dF) of the mode for unit dL at Pi of the chart (par, tab)."""
        m = self.m
        x, d_, r, rV, rU, pV, pU, V0, V1, V2, E, FV, FU, rhoV, rhoU = hs.jetsK(1.0, -Pi, par, tab, 0.0)
        rho, rho_P = x, m.kappa * E
        q, qp = self.q_qp(rho)
        R, Rp = m.R(rho), m.Rp(rho)
        dr = R / 2.0 + Rp * q
        dr_P = rho_P * (Rp / 2.0 + q / R ** 3 + Rp * qp)
        dphi = SQ2 * q / R ** 2
        dF = 0.0 if MUT == "mode_dF0" else self.khat - qp
        return dr, dr_P, dphi, dF

    def rho_q(self, Pi, par, tab):
        x = hs.jetsK(1.0, -Pi, par, tab, 0.0)[0]
        return x, self.q_qp(x)[0]

    def dm_exact(self, rho, q):
        """d_L m at fixed Pi: m_L(rho) = sqrt(L)(1 + L rho^2 (1 - A))/(2R), dm = d_L m|rho + m_rho q."""
        m = self.m
        R = math.sqrt(1.0 + rho ** 2)
        omA = -m.Am1(rho)
        dmL = 1.0 / (4 * R) + 3.0 * rho ** 2 * omA / (4 * R)
        m_rho = (2 * rho * omA - rho ** 2 * m.Ap(rho)) / (2 * R) - rho * (1.0 + rho ** 2 * omA) / (2 * R ** 3)
        return dmL + m_rho * q

    def q_big(self, rho):
        return self.big.sol(math.log(rho))[0]


class SchwMode:
    """r_S -> r_S (1 + eps) at fixed Pi, (r/r_S - 1) e^{r/r_S} = Pi: dr = r, dF = -2, dphi = 0 per unit d r_S (r_S = 1); c1 = r_S, dM_ADM = 1/2."""

    def __init__(self, m):
        self.m = m

    def fields(self, Pi, par, tab):
        x, d_, r, rV, rU, pV, pU, V0, V1, V2, E, FV, FU, rhoV, rhoU = hs.jetsK(1.0, -Pi, par, tab, 0.0)
        return r, 0.5 * E, 0.0, -2.0

    def rho_q(self, Pi, par, tab):
        return hs.jetsK(1.0, -Pi, par, tab, 0.0)[0], 0.0

    def dm_exact(self, rho, q):
        return 0.5


_MODE = {}


def mode_of(m):
    if m.tag not in _MODE:
        _MODE[m.tag] = SchwMode(m) if m.kind == "schw" else EllisMode(m)
    return _MODE[m.tag]


def mode_fields(m, Pi, chart=None):
    """(dr, dr_Pi, dphi, dF) of the mode at Pi in the member's chart `chart` (default: the tabulated slice chart)."""
    bg = chart or bg_of(m)
    return mode_of(m).fields(Pi, bg.par, bg.tab)


# ----------------------------------------------------------------------------------------------------------------------
# references: mpmath (independent of the chart and of the double-precision series)
# ----------------------------------------------------------------------------------------------------------------------
def mp_family(c, dps=120):
    """closures of the Ellis member c at dps digits: A, A', R, horizon rho_h(L) of F_L = A + 1/L - 1, area, M_MS(L) at rho = 1e20."""
    import mpmath as mp
    mp.mp.dps = dps
    c = mp.mpf(c)
    K = -2 * c / mp.pi
    A = lambda r: 1 + c * (r * r + 1) + K * (r + (r * r + 1) * mp.atan(r))
    Ap = lambda r: 2 * c * r + 2 * K * (1 + r * mp.atan(r))

    def rhL(L):
        return mp.findroot(lambda r: A(r) - (1 - 1 / L), (mp.mpf(-60), mp.mpf(60)), solver="illinois", tol=mp.mpf(10) ** (-dps + 10), maxsteps=400)

    def M_MS(L):
        rho = mp.mpf(10) ** 20
        return mp.sqrt(L) * (1 + L * rho ** 2 * (1 - A(rho))) / (2 * mp.sqrt(1 + rho ** 2))

    return mp, A, Ap, rhL, M_MS, K


def exact_ref(m):
    """(c1, dM) of the exact mode per unit dL (Ellis) or per unit d r_S (Schwarzschild): c1 = d_L [sqrt(L) R(rho_h(L))] and dM = d_L lim m_L, central differences at 120 digits."""
    if m.kind == "schw":
        return 1.0, 0.5
    mp, A, Ap, rhL, M_MS, K = mp_family(m.c)
    e = mp.mpf(10) ** -20
    c1 = (mp.sqrt(1 + e) * mp.sqrt(1 + rhL(1 + e) ** 2) - mp.sqrt(1 - e) * mp.sqrt(1 + rhL(1 - e) ** 2)) / (2 * e)
    dM = (M_MS(1 + e) - M_MS(1 - e)) / (2 * e)
    return float(c1), float(dM)


def c1_closed(m):
    """c1 per unit dL: R_h/2 + R'_h rho_h'(L), rho_h' = 1/A'(rho_h) = 1/(2 kappa) (implicit differentiation of A(rho_h) = 1 - 1/L)."""
    if m.kind == "schw":
        return 1.0
    shift = 0.0 if MUT == "c1_noshift" else m.Rph / (2.0 * m.kappa)
    return m.rS / 2.0 + shift


# ----------------------------------------------------------------------------------------------------------------------
# evolved readings of the Kruskal solver: the exact mode as ray data, #248's classes
# ----------------------------------------------------------------------------------------------------------------------
def slice_vals(U, h, s):
    j = int(round((s - V0R) / h))
    i = int(round((-s - U0R) / h))
    return U[0, j, i], hs.fd4(U[0], j, i, 1, -1, h), U[2, j, i], U[1, j, i]


def kruskal_readings(m, U, h, ih):
    """c1 and c2 of the horizon fit, dr at S, Theta_gen on the slice, and dM_ADM of the cut extensions."""
    sl = slice_of(m)
    V = V0R + h * np.arange(U.shape[1])
    c1, c2, res = hs.horizon_fit(U[0, :, ih], V)
    jS = int(round(-V0R / h))
    out = dict(c1fit=c1, c2fit=c2, res=res, c1S=U[0, jS, ih], theta={}, cut={})
    for s in S_LIST:
        dr, drs, dphi, dF = slice_vals(U, h, s)
        out["theta"][s] = theta_gen(sl, s, dr, drs, dphi, dF)
    for sc in S_CUTS + S_CUTL:
        dr, drs, dphi, dF = slice_vals(U, h, sc)
        rm = read_mass(m, sc, dr, drs, dphi, dF)
        out["cut"][sc] = dict(dM=rm["dM"], th_before=rm["th_before"], kth_after=rm["kth_after"])
    return out


_RUNS = {}


def class_run(m, name, h):
    key = (m.tag, name, h)
    if key not in _RUNS:
        su = hs.KSetup(m.c if m.kind == "ellis" else -1.0, h, *RECT, hs.CLASSES[name], kc=solver_chart(m))
        U = su.run()
        _RUNS[key] = kruskal_readings(m, U, h, su.ih)
        del U
    return _RUNS[key]


def mode_arrays(m, h, amp, ch):
    md = mode_of(m)
    Nv = int(round((RECT[1] - RECT[0]) / h))
    Nu = int(round((RECT[3] - RECT[2]) / h))
    UL = np.zeros((3, Nu + 1))
    UR = np.zeros((3, Nv + 1))
    for i in range(Nu + 1):
        dr, _, dphi, dF = md.fields(-V0R * (U0R + i * h), ch.par, ch.tab)
        UL[:, i] = amp * np.array([dr, dF, dphi])
    for j in range(Nv + 1):
        dr, _, dphi, dF = md.fields(-(V0R + j * h) * U0R, ch.par, ch.tab)
        UR[:, j] = amp * np.array([dr, dF, dphi])
    return UL, UR


def mode_run(m, h, amp=2e-3):
    key = (m.tag, "mode", h, amp)
    if key not in _RUNS:
        ch = solver_chart(m)
        UL, UR = mode_arrays(m, h, amp, ch)
        U = hs.march_K(UL, UR, V0R, U0R, ch.par, ch.tab, h, h)
        ih = int(round(-U0R / h))
        rd = kruskal_readings(m, U, h, ih)
        md = mode_of(m)
        st = max(1, int(round(0.05 / h)))
        err, scale = 0.0, 0.0
        for j in range(0, U.shape[1], st):
            for i in range(0, U.shape[2], st):
                Pi = -(V0R + j * h) * (U0R + i * h)
                if Pi < -0.24 or Pi > 1e9:
                    continue
                f = md.fields(Pi, ch.par, ch.tab)
                ex = amp * np.array([f[0], f[3], f[2]])
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


def class_summary(m, name, hl=HL):
    r = {h: class_run(m, name, h) for h in hl}
    h1, h2 = hl[-2], hl[-1]
    out = dict(c1=rich(r[h1]["c1fit"], r[h2]["c1fit"]), c1S=rich(r[h1]["c1S"], r[h2]["c1S"]), c2=rich(r[h1]["c2fit"], r[h2]["c2fit"]))
    out["theta"] = {s: rich(r[h1]["theta"][s], r[h2]["theta"][s]) for s in S_LIST}
    out["dM"] = {sc: rich(r[h1]["cut"][sc]["dM"], r[h2]["cut"][sc]["dM"]) for sc in S_CUTS + S_CUTL}
    out["res"] = r[h2]["res"]
    return out


# ----------------------------------------------------------------------------------------------------------------------
# bookkeeping
# ----------------------------------------------------------------------------------------------------------------------
RESULTS = []


def check(cid, text, ok, detail=""):
    if re.search(r"(?<![a-z])(nan|-?inf)(?![a-z])", str(detail).lower()):   # max() and <= hide a NaN
        ok = False
    RESULTS.append((cid, bool(ok)))
    print(f"[{'pass' if ok else 'FAIL'}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def sci(x):
    return f"{x:.2e}"


TAGS = ("-1", "-1.5", "-0.3", "schw")                # c = -1 is #259's member, the control; the other three are the issue's


# ----------------------------------------------------------------------------------------------------------------------
# derive: sympy
# ----------------------------------------------------------------------------------------------------------------------
def derive_S1():
    """The ρ-form of the sum constraint from the chart's Q, and its first integral, for the Ellis members (generic A with R^2 A'' = 2(A - 1)) and Schwarzschild."""
    import sympy as sp
    r = sp.symbols("rho", positive=True)
    u = sp.Function("u")(r)
    out = []
    Af = sp.Function("A")(r)
    R = sp.sqrt(1 + r ** 2)
    cases = []
    App = 2 * (Af - 1) / R ** 2
    Ve = -App / 2 - sp.diff(R, r) * sp.diff(Af, r) / R                    # (v rho): V = -A''/2 - A' R'/R
    cases.append(("Ellis", R, Af, Ve, 2 / R ** 4, App))
    As = 1 - 1 / r
    cases.append(("Schwarzschild", r, As, sp.Integer(0), sp.Integer(0), sp.diff(As, r, 2)))
    for name, Rr, Aa, V, php2, App_ in cases:
        Rp = sp.diff(Rr, r)
        q = (Aa / 2) * (V + 1 / Rr ** 2 - Aa * Rp ** 2 / Rr ** 2 - Aa * php2 / 2) * (-1 if MUT == "s_q" else 1)   # kappa^2 Pi Q of the chart
        p = -sp.diff(Aa, r) / 2 + Aa * Rp / Rr
        upp = -((Aa * sp.diff(Aa, r) + p * Aa) * sp.diff(u, r) + q * u) / Aa ** 2        # A^2 u'' + A A' u' + p A u' + q u = 0
        Th = Rr * (Aa * sp.diff(u, r) - sp.diff(Aa, r) * u / 2)
        dTh = sp.diff(Th, r).subs(sp.diff(u, r, 2), upp)
        if name == "Ellis":
            dTh = dTh.subs(sp.diff(Af, r, 2), App_)
        out.append((name, sp.simplify(dTh)))
    ok = all(res == 0 for _, res in out)
    check("S1", "chart Q -> rho-form q = -(A/2)(A'' + A' R'/R), p = -A'/2 + A R'/R; d/drho[R (A y' - A' y/2)] = 0 on its solutions (symbolic), Ellis members (generic A, R^2 A'' = 2(A - 1), V from (v rho)) and Schwarzschild",
          ok, "; ".join(f"{n}: residual {res}" for n, res in out))


def derive_S2():
    """a = -R A'/(2 kappa) (so a(0) = -r_S) and (mu s r_Pi)' = mu E (1/r - V r) for generic R(rho), A(rho), V R^2 = 1 - (A (R^2)')'/2 ((theta theta))."""
    import sympy as sp
    s, kap = sp.symbols("s kappa", positive=True)
    x = sp.symbols("rho", real=True)
    Af = sp.Function("A")(x)
    Rf = sp.Function("R")(x)
    ds = lambda f: (Af / (kap * s)) * sp.diff(f, x)                       # d/ds = rho_s d/drho, rho_s = A/(kappa s), Pi = s^2
    rho_Pi = Af / (2 * kap * s ** 2)
    mu = s * Rf
    E = Af / (2 * kap ** 2 * s ** 2)
    mu_s_rPi = sp.simplify(mu * s * sp.diff(Rf, x) * rho_Pi)              # mu s r_Pi = R R' A/(2 kappa)
    V = (1 - sp.diff(Af * sp.diff(Rf ** 2, x), x) / 2) / Rf ** 2 + (1 / Rf ** 2 if MUT == "s_pot" else 0)
    res1 = sp.simplify(ds(mu_s_rPi) - mu * E * (1 / Rf - V * Rf))
    dF_Pi = -sp.diff(Af, x) * rho_Pi / Af + 1 / s ** 2                    # F = -ln E = -ln A + ln(2 kappa^2 Pi), E = e^{-F} = A/(2 kappa^2 Pi)
    a = Rf * (s ** 2 * dF_Pi - 1)                                         # a = mu P - mu_s = r (s^2 F_Pi - 1)
    res2 = sp.simplify(a + Rf * sp.diff(Af, x) / (2 * kap))
    check("S2", "background identity (mu s r_Pi)' = mu E (1/r - V r) with V from (theta theta), and a = r (s^2 F_Pi - 1) = -R A'/(2 kappa) (so a(0) = -R(rho_h) since A'(rho_h) = 2 kappa), for generic R(rho), A(rho)",
          res1 == 0 and res2 == 0, f"residuals {res1}, {res2}")


def derive_S3():
    """kappa r_S c1 = K/2 = dM_ADM/dL for every Ellis member; K is constant on the static solutions; M = K/3 at L = 1."""
    import sympy as sp
    r, c = sp.symbols("rho c", real=True)
    Af = sp.Function("A")(r)
    R2 = 1 + r ** 2
    Kexpr = (R2 * sp.diff(Af, r) - 2 * r * Af) / 2 + (0 if MUT == "s_k" else r)
    dK = sp.simplify(sp.diff(Kexpr, r).subs(sp.diff(Af, r, 2), 2 * (Af - 1) / R2))
    Rh, rh, A1 = sp.symbols("R_h rho_h A1", positive=True)
    c1 = Rh / 2 + (rh / Rh) / A1                                           # d_L [sqrt(L) R(rho_h(L))], rho_h' = 1/A'(rho_h)
    kap = A1 / 2
    Kh = sp.Rational(1, 2) * Rh ** 2 * A1 + rh                             # K at A = 0
    ident = sp.simplify(kap * Rh * c1 - Kh / 2)
    Kc = sp.symbols("K", positive=True)
    t = sp.symbols("t", positive=True)
    rho = 1 / t
    Ac = 1 + (-Kc * sp.pi / 2) * (1 + rho ** 2) + Kc * (rho + (1 + rho ** 2) * (sp.pi / 2 - sp.atan(t)))
    Mser = sp.series(sp.simplify(rho * (1 - Ac) / 2), t, 0, 2).removeO()
    rs, M = sp.symbols("r_S M", positive=True)
    schw = sp.simplify((1 / (2 * rs)) * rs * 1 * 2 - 1)                    # kappa r_S c1 / dM at r_S = 2M with c1 = d r_S = 2 dM, kappa = 1/(4M): (1/(2 rs)) rs (d rs) = d rs / 2 = dM
    check("S3", "K = (R^2 A' - 2 rho A)/2 + rho is constant on R^2 A'' = 2 (A - 1); kappa r_S c1 = K/2 with c1 = R_h/2 + R'_h/A'(rho_h); lim rho (1 - A)/2 = K/3 so dM/dL = K/2 (M ~ L^{3/2}); Schwarzschild kappa r_S c1 = dM",
          dK == 0 and ident == 0 and sp.simplify(Mser - Kc / 3) == 0 and schw == 0, f"dK/drho = {dK}; identity {ident}; M(L=1) = {Mser}; Schwarzschild {schw}")


def run_derive():
    t0 = time.time()
    for fn in (derive_S1, derive_S2, derive_S3):
        fn()
    bad = [c for c, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS) - len(bad)} of {len(RESULTS)} checks pass" + (f"; fired: {', '.join(bad)}" if bad else "") + f"  [{time.time() - t0:.0f} s]  mutation: {MUT}")
    return 1 if bad else 0


# ----------------------------------------------------------------------------------------------------------------------
# numerical checks
# ----------------------------------------------------------------------------------------------------------------------
def check_T1():
    """members: constants against mpmath, c1 closed form against d_L r_S(L), kappa r_S c1 against -c/pi."""
    import mpmath as mp
    worst_c, worst_id, worst_cl = 0.0, 0.0, 0.0
    det = []
    for tag in TAGS:
        m = member(tag)
        if m.kind == "schw":
            c1x, dMx = exact_ref(m)
            ident = abs(kap_c(m) * rs_c(m) * c1_closed(m) - dMx) / dMx
            worst_id = max(worst_id, ident)
            det.append(f"Schwarzschild: kappa {m.kappa}, r_S {m.rS}, c1 {c1_closed(m):.12f}, dM {dMx}, identity {ident:.1e}")
            continue
        mp.mp.dps = 40
        c = mp.mpf(m.c)
        K = -2 * c / mp.pi
        A = lambda r: 1 + c * (r * r + 1) + K * (r + (r * r + 1) * mp.atan(r))
        Ap = lambda r: 2 * c * r + 2 * K * (1 + r * mp.atan(r))
        rh = mp.findroot(A, (mp.mpf(m.rh) - mp.mpf("0.01"), mp.mpf(m.rh) + mp.mpf("0.01")), solver="illinois", tol=mp.mpf(10) ** -35)
        kap = Ap(rh) / 2
        rS = mp.sqrt(1 + rh ** 2)
        e_c = max(abs(m.rh - float(rh)), abs(m.kappa - float(kap)) / float(kap), abs(m.rS - float(rS)) / float(rS))
        c1x, dMx = exact_ref(m)
        e_cl = abs(c1_closed(m) - c1x) / c1x
        ident = abs(kap_c(m) * rs_c(m) * c1_closed(m) - dMx) / dMx
        e_K = abs(dMx - (-c / mp.pi)) / abs(float(-c / mp.pi))
        worst_c, worst_cl, worst_id = max(worst_c, e_c), max(worst_cl, e_cl, e_K), max(worst_id, ident)
        det.append(f"{m.tag}: rho_h {m.rh:.9f}, kappa {m.kappa:.9f}, r_S {m.rS:.9f}, R'_h {m.Rph:.9f}, c1 {c1_closed(m):.12f} (d_L r_S {c1x:.12f}), dM {dMx:.12f} (-c/pi {float(-c / mp.pi):.12f}), kappa r_S c1 vs dM {ident:.1e}")
    check("T1", "members: rho_h, kappa, r_S against mpmath (40 digits); c1 = R_h/2 + R'_h/(2 kappa) against the derivative of r_S(L) = sqrt(L) R(rho_h(L)) (120 digits); dM_ADM against -c/pi; kappa r_S c1 = dM_ADM from the constants",
          worst_c <= 1e-12 and worst_cl <= 1e-12 and worst_id <= 1e-12, "; ".join(det))


def check_T2():
    """first law dM = kappa_xi dA/(8 pi) along each family at finite L, mpmath at 120 digits."""
    import mpmath as mp
    worst, det = 0.0, []
    for tag in ("-1.5", "-0.3", "-1"):
        c = float(tag)
        mp_, A, Ap, rhL, M_MS, K = mp_family(c)
        mp_.mp.dps = 120
        e = mp_.mpf(10) ** -15
        area = lambda L: 4 * mp_.pi * L * (1 + rhL(L) ** 2)
        kap = lambda L: (1 if MUT == "law_norm" else mp_.sqrt(L)) * Ap(rhL(L)) / 2          # xi = sqrt(L) d_t is unit at infinity (F_L -> 1/L)
        w, rows, dm_route = 0.0, [], 0.0
        for Ls in ("0.5", "0.9", "1", "1.5", "2", "4"):
            L = mp_.mpf(Ls)
            Mc = lambda L_: (K * L_ ** (2 if MUT == "law_mass" else mp_.mpf(1.5)) / 3)
            dM = (M_MS(L + e) - M_MS(L - e)) / (2 * e)
            dMc = (Mc(L + e) - Mc(L - e)) / (2 * e)
            dA = (area(L + e) - area(L - e)) / (2 * e)
            ratio = kap(L) * dA / (8 * mp_.pi) / dM
            w = max(w, abs(float(ratio) - 1.0))
            dm_route = max(dm_route, abs(float(dM / dMc) - 1.0))
            rows.append(f"{Ls}: {float(ratio):.13f}")
        worst = max(worst, w, dm_route)
        det.append(f"{tag}: {', '.join(rows)}; M_MS vs K L^(3/2)/3 {dm_route:.1e}")
    for rs in ("0.5", "1", "2"):
        mp.mp.dps = 60
        r = mp.mpf(rs)
        Mf = lambda x: x / 2                                                                 # m = (r/2)(1 - (1 - r_S/r)) = r_S/2
        kapf = lambda x: 1 / (2 * x)
        e = mp.mpf(10) ** -20
        ratio = kapf(r) * (4 * mp.pi * ((r + e) ** 2 - (r - e) ** 2)) / (2 * e) / (8 * mp.pi) / ((Mf(r + e) - Mf(r - e)) / (2 * e))
        worst = max(worst, abs(float(ratio) - 1.0))
        det.append(f"Schwarzschild r_S = {rs}: {float(ratio):.13f}")
    check("T2", "first law dM = kappa_xi dA/(8 pi) along the exact families (F_L = A + 1/L - 1, r = sqrt(L) R, kappa_xi = sqrt(L) A'(rho_h(L))/2, A = 4 pi L R(rho_h)^2; M = Misner-Sharp limit at rho = 1e20, 120 digits) at L = 0.5 ... 4; Schwarzschild at r_S = 0.5, 1, 2; |ratio - 1| <= 1e-12 (bar of the issue)",
          worst <= 1e-12, "ratios " + "; ".join(det))


def check_T3():
    """rho-form of H against the chart's P, Q; a = -R A'/(2 kappa); Theta and dm in rho against the s-form; closed-form tail against the ODE; Theta conserved in rho."""
    det, worst = [], 0.0
    for tag in TAGS:
        m = member(tag)
        bg, sl = bg_of(m), slice_of(m)
        smax = 12.0 if m.kind == "schw" else 1e4
        wp = wq = wa = 0.0
        for s in (0.1, 0.5, 1.0, 2.0, 5.0, 10.0) + ((30.0, 1e3, 1e5) if m.kind != "schw" else ()):
            P, Q, b = sl.PQ(s)
            rho = b["rho"]
            p, q = xform_pq(m, rho)
            wp = max(wp, abs((m.kappa * s * P - m.kappa) - p) / abs(p))
            wq = max(wq, abs(m.kappa ** 2 * s * s * Q - q) / abs(q))
            mu, a, cc, Q_, b_ = sl.theta_coef(s)
            wa = max(wa, abs(a + m.R(rho) * m.Ap(rho) / (2.0 * m.kappa)) / abs(a))
        # a(0): at the smallest s the coefficient is -r_S to O(s)
        mu0, a0, _, _, b0 = sl.theta_coef(1e-6)
        e_a0 = abs(a0 + m.rS) / m.rS
        h = sl.hom_even(fs.SL, smax, (1.0, 0.0))
        wth = wdm = 0.0
        for s in (1.05, 3.0, 10.0) + ((30.0, 1e3, 1e4) if m.kind != "schw" else ()):
            dr, drs = h.sol(s)
            rho = bg.at(s)["rho"]
            yx = m.kappa * s * drs
            kth_x = ktheta_x(m, rho, dr, yx)
            kth_s = m.kappa * theta_gen(sl, s, dr, drs, 0.0, 0.0)
            wth = max(wth, abs(kth_x - kth_s) / abs(kth_s))
            wdm = max(wdm, abs(dm_x(m, rho, dr, yx) - dm_s(bg, s, dr, drs, 0.0)) / abs(dm_s(bg, s, dr, drs, 0.0)))
        s0 = 1.05
        dr0, drs0 = h.sol(s0)
        rho0 = bg.at(s0)["rho"]
        yx0 = m.kappa * s0 * drs0
        tl = Tail(m, rho0, dr0, yx0)
        kth0 = ktheta_x(m, rho0, dr0, yx0)
        wcl = wcon = 0.0
        for rho in (3.0, 16.0, 1e2, 1e4, 1e6, RHO_MAX):
            a_ = tl.at(rho)
            b_ = tail_closed(m, rho0, dr0, yx0, rho)
            wcl = max(wcl, abs(a_[0] - b_[0]) / abs(b_[0]), abs(a_[1] - b_[1]) / abs(b_[1]))
            wcon = max(wcon, abs(ktheta_x(m, rho, *a_) - kth0) / abs(kth0))
        worst = max(worst, wp, wq, wa, e_a0, wth, wdm, wcl, wcon)
        det.append(f"{m.tag}: p {wp:.1e}, q {wq:.1e}, a {wa:.1e}, a(0)+r_S {e_a0:.1e}, Theta {wth:.1e}, dm {wdm:.1e}; closed form {wcl:.1e}; conservation {wcon:.1e}")
    check("T3", "rho-form p, q against the chart's P, Q (6-9 points); a = -R A'/(2 kappa) and a(s = 1e-6) = -r_S; kappa Theta and dm in rho against the s-form on a homogeneous solution of H (6 points); closed-form tail against the ODE; kappa Theta conserved on [rho(1.05), 1e8]; bar 1e-8",
          worst <= 1e-8, "; ".join(det))


S_MODE = (0.05, 0.1, 0.3, 0.7, 1.0, 2.0, 5.0, 10.0, 20.0, 1e2, 1e3, 1e4, 3e4, 1e5)
S_MODE_SCHW = (0.05, 0.1, 0.3, 0.7, 1.0, 2.0, 5.0, 10.0)
S_CUT_MODE = (0.5, 2.0, 20.0, 1e3, -0.5, -2.0, -20.0)
S_CUT_MODE_SCHW = (0.5, 2.0, 10.0, -0.5, -2.0, -10.0)


def check_T4():
    """the exact static mode through the slice code, per member: c1 from Theta_gen at one point, the mass of the cut extensions, kappa r_S c1 = dM_ADM."""
    det, ok = [], True
    for tag in TAGS:
        m = member(tag)
        md, bg, sl = mode_of(m), bg_of(m), slice_of(m)
        c1x, dMx = exact_ref(m)
        c1 = md.fields(0.0, bg.par, bg.tab)[0]
        wth = wdm = wc1 = 0.0
        for s in (S_MODE_SCHW if m.kind == "schw" else S_MODE):
            dr, dr_P, dphi, dF = md.fields(s * s, bg.par, bg.tab)
            drs = 2 * s * dr_P
            th = theta_gen(sl, s, dr, drs, dphi, dF)
            wth = max(wth, abs(th + m.rS * c1x) / (m.rS * c1x))
            wc1 = max(wc1, abs(-th / m.rS - c1x) / c1x)
            rho, q = md.rho_q(s * s, bg.par, bg.tab)
            dx = md.dm_exact(rho, q)
            wdm = max(wdm, abs(dm_s(bg, s, dr, drs, dF) - dx) / abs(dx))
        e_S = abs(c1 - c1x) / c1x
        errs, errs_ex = [], []
        for sc in (S_CUT_MODE_SCHW if m.kind == "schw" else S_CUT_MODE):
            dr, dr_P, dphi, dF = md.fields(sc * sc, bg.par, bg.tab)
            rm = read_mass(m, sc, dr, 2 * sc * dr_P, dphi, dF)
            errs.append(abs(rm["dM"] - kap_c(m) * rs_c(m) * c1x) / dMx)
            errs_ex.append(abs(rm["dM"] - dMx) / dMx)
        e_kc = abs(kap_c(m) * rs_c(m) * c1x - dMx) / dMx
        good = e_S <= 1e-9 and wth <= 1e-9 and wc1 <= 1e-9 and wdm <= 1e-9 and e_kc <= BAR and max(errs) <= BAR and max(errs_ex) <= BAR
        ok = ok and good
        det.append(f"{m.tag}: c1(S) {c1:.12f} (exact {c1x:.12f}, {e_S:.1e}); c1 from Theta_gen at 9-14 s {wc1:.1e}; Theta + r_S c1 {wth:.1e}; dm slice vs exact {wdm:.1e}; dM exact {dMx:.12f}; kappa r_S c1 vs dM {e_kc:.1e}; "
                   f"cut extensions vs kappa r_S c1 {max(errs):.1e}, vs dM {max(errs_ex):.1e}")
    check("T4", "exact static mass mode through the slice code (planted truth), per member: c1 = dr(S) against the mpmath d_L r_S; Theta_gen = -r_S c1 and c1 = -Theta_gen/r_S at 9-14 values 0.05 <= s <= 1e5 (relative, bar 1e-9); dm of the slice code against d_L m at fixed Pi; "
          "kappa r_S c1 = dM_ADM (mpmath d_L lim m_L) and the 6-7 cut extensions read at rho = 1e8 (bar 1e-6, the issue's)", ok, "; ".join(det))


def check_T5():
    """negative control: with r_S set to 1 the identity fails at c != -1 (and at Schwarzschild r_S != 1), and holds at the members where r_S = 1."""
    det, ok = [], True
    for tag in TAGS:
        m = member(tag)
        c1x, dMx = exact_ref(m)
        r_used = m.rS if MUT == "ctl_true" else 1.0
        e = abs(kap_c(m) * r_used * c1x - dMx) / dMx
        if tag in ("-1.5", "-0.3"):
            good = e > 1e-3
        else:
            good = e <= 1e-12                                              # r_S = 1 there: the "wrong" value is the right one
        ok = ok and good
        det.append(f"{m.tag}: r_S {m.rS:.6f}, kappa c1 vs dM {e:.2e} ({'fails' if e > 1e-3 else 'holds'})")
    rs = 2.0                                                               # Schwarzschild 2M = 2 per unit d r_S: kappa = 1/(2 r_S) = 1/4, c1 = 1, dM = 1/2 (closed forms: the solver chart is fixed at r_S = 1)
    e2 = abs(0.25 * (rs if MUT == "ctl_true" else 1.0) * 1.0 - 0.5) / 0.5
    ok = ok and e2 > 1e-3
    det.append(f"Schwarzschild 2M = 2: kappa c1 vs dM {e2:.2e} ({'fails' if e2 > 1e-3 else 'holds'})")
    check("T5", "wrong r_S (set to 1) in kappa r_S c1 = dM_ADM fails by more than 1e-3 at c = -1.5, -0.3 and at Schwarzschild 2M = 2 (an added case: at 2M = 1, r_S = 1 and nothing can fail), and holds to 1e-12 at c = -1 and Schwarzschild 2M = 1 where r_S = 1",
          ok, "; ".join(det))


def gauge_c2_mp(c, eta0):
    """c2 = d_V dr on U = 0 of the pure gauge xi = eta(U) d_U with eta(0) = eta0: -eta0 R'(rho_h) G(rho_h), G = d rho/d Pi at Pi = 0 = E(rho_h) with Pi(rho = 1) = 1 (mpmath, 60 digits, from A alone)."""
    import mpmath as mp
    mp.mp.dps = 60
    cc = mp.mpf(c)
    K = -2 * cc / mp.pi
    A = lambda r: 1 + cc * (r * r + 1) + K * (r + (r * r + 1) * mp.atan(r))
    Ap = lambda r: 2 * cc * r + 2 * K * (1 + r * mp.atan(r))
    rh = mp.mpf(0) if c == -1.0 else mp.findroot(A, (mp.mpf(-3), mp.mpf(3)), solver="illinois", tol=mp.mpf(10) ** -50)
    kap = Ap(rh) / 2
    lnE_h = mp.log(1 - rh) + mp.quad(lambda r: 1 / (r - rh) - 2 * kap / A(r), [1, (1 + rh) / 2, rh], method="gauss-legendre")   # interior nodes only: a difference of poles
    G = mp.e ** lnE_h
    Rp = 1.0 if MUT == "gauge_noRp" else rh / mp.sqrt(1 + rh ** 2)
    return float(-eta0 * Rp * G)


def check_T6():
    """gauge: Theta_gen = 0 on the exact pure gauge L_xi(member) (dF != 0), dr(S) = 0, and c2 shifts by -eta(0) R'(rho_h) G(0) while c1 does not."""
    det, ok = [], True
    for tag in TAGS:
        m = member(tag)
        bg, sl = bg_of(m), slice_of(m)
        gf = fs.gauge_fun(bg)
        fn = fs.fd_slice_fn(gf)
        w = 0.0
        for s in (-3.0, -1.5, -0.3, 0.1, 0.5, 1.0, 3.0, 8.0):
            U0, UT, UN, _ = fn(s)
            w = max(w, abs(theta_gen(sl, s, U0[0], UT[0], U0[2], U0[1])) / theta_scale(sl, s, U0[0], UT[0], U0[2], U0[1]))
        eta0 = fs.ETA(0.0)[0]
        d_S = abs(gf(0.0, 0.0)[0])
        Vs = np.linspace(0.0, 2.0, 9)
        dh = np.array([gf(V, 0.0)[0] for V in Vs])
        c1h, c2h, res = hs.horizon_fit(dh, Vs)
        if m.kind == "schw":
            good = w <= 1e-8 and d_S <= 1e-14 and abs(c1h) <= 1e-12 and abs(c2h) > 1e-4
            det.append(f"{m.tag}: Theta {w:.1e}, dr(S) {d_S:.1e}, c1 of the gauge solution {c1h:.1e}, c2 {c2h:.6e} (R'_h = 1; pure gauge)")
        else:
            c2x = gauge_c2_mp(m.c, eta0)
            e2 = abs(c2h - c2x) / max(abs(c2x), 1e-300) if m.Rph != 0.0 else abs(c2h)
            good = w <= 1e-8 and d_S <= 1e-14 and abs(c1h) <= 1e-12 and (e2 <= 1e-8 if m.Rph != 0.0 else e2 <= 1e-12)
            det.append(f"{m.tag}: Theta {w:.1e}, dr(S) {d_S:.1e}, c1 of the gauge solution {c1h:.1e}, c2 {c2h:.6e} vs -eta(0) R'_h G {c2x:.6e} ({e2:.1e})")
        ok = ok and good
    check("T6", "exact pure gauge xi = xi(V) d_V + eta(U) d_U (eta(0) = 0.03, delta F != 0): Theta_gen = 0 at s = -3 ... 8 (over its largest term, bar 1e-8); dr(S) = 0 and c1 = 0 of the horizon fit (bar 1e-12); c2 = -eta(0) R'(rho_h) G(0) (mpmath, bar 1e-8; 0 at c = -1 where R'(rho_h) = 0; Schwarzschild: c2 != 0)",
          ok, "; ".join(det))


def check_T7():
    """the exact mode as ray data through the Kruskal solver, per member."""
    amp = 2e-3
    det, ok = [], True
    for tag in TAGS:
        m = member(tag)
        rs_ = {h: mode_run(m, h, amp) for h in HL}
        errs = [rs_[h]["err"] / rs_[h]["scale"] for h in HL]
        rt = ratios(errs)
        h1, h2 = HL[-2], HL[-1]
        c1x, dMx = exact_ref(m)
        c1_exact, dM_exact = amp * c1x, amp * dMx
        c1 = rich(rs_[h1]["c1fit"], rs_[h2]["c1fit"])
        c2 = rich(rs_[h1]["c2fit"], rs_[h2]["c2fit"])
        e_c1 = abs(c1 - c1_exact) / c1_exact
        e_S = abs(rich(rs_[h1]["c1S"], rs_[h2]["c1S"]) - c1_exact) / c1_exact
        e_th = max(abs(rich(rs_[h1]["theta"][s], rs_[h2]["theta"][s]) + m.rS * c1_exact) / (m.rS * c1_exact) for s in S_LIST)
        cuts = {sc: rich(rs_[h1]["cut"][sc]["dM"], rs_[h2]["cut"][sc]["dM"]) for sc in S_CUTS + S_CUTL}
        e_m = max(abs(v - dM_exact) / dM_exact for v in cuts.values())
        e_k = max(abs(v - kap_c(m) * rs_c(m) * c1) / (kap_c(m) * rs_c(m) * abs(c1)) for v in cuts.values())
        c2_ok = abs(c2) <= 1e-6 * c1_exact if m.kind == "ellis" else True
        good = errs[-1] <= 5e-2 and abs(rt[-1] - 4.0) <= 0.3 and abs(rt[-2] - 4.0) <= 0.3 and e_c1 <= BAR and e_S <= BAR and e_th <= BAR and e_m <= BAR and e_k <= BAR and c2_ok
        ok = ok and good
        det.append(f"{m.tag}: errors {', '.join(sci(e) for e in errs)}, ratio {rt[-1]:.2f}; c1 fit {e_c1:.1e}, dr(S) {e_S:.1e}, Theta_gen {e_th:.1e}, dM cuts vs exact {e_m:.1e}, kappa r_S c1 vs dM cuts {e_k:.1e}; c2 {c2:.1e}")
    check("T7", "exact static mode (amplitude 2e-3) as ray data through the Kruskal solver at each member, h = 0.005 ... 0.000625 (Richardson of the last two): error/scale at the finest step <= 5e-2 with the last two ratios 4.0 +- 0.3; c1 of the horizon fit, dr(S), Theta_gen(s) = -r_S c1 at 8 values of s, dM_ADM of the three cut extensions against the exact value and against kappa r_S c1, all <= 1e-6 (relative); c2 = 0 for Ellis",
          ok, "; ".join(det))


CLASSES_T8 = ("A0", "A1R", "A1L", "G")


def check_T8():
    """#248's classes through the slice code at each member: Theta_gen(s) = -r_S c1, kappa r_S c1 = dM_ADM of the cut extensions."""
    det, ok = [], True
    for tag in TAGS:
        m = member(tag)
        for name in CLASSES_T8:
            if m.kind == "schw" and name == "G":
                continue                                                    # no scalar: G is the zero datum for dr
            sm = class_summary(m, name)
            c1 = sm["c1"]
            scale = max(abs(c1), 1e-5)                                      # relative to |c1|; a class with |c1| < 1e-5 (1% of the amplitude) is held to the absolute bar 1e-6 * 1e-5
            e_th = max(abs(sm["theta"][s] + m.rS * c1) for s in S_LIST) / (m.rS * scale)
            e_S = abs(sm["c1S"] - c1) / scale
            e_m = max(abs(v - kap_c(m) * rs_c(m) * c1) for v in sm["dM"].values()) / (kap_c(m) * rs_c(m) * scale)
            good = e_th <= BAR and e_S <= BAR and e_m <= BAR
            det.append(f"{m.tag} {name}: c1 {c1:.7e}, kappa r_S c1 {kap_c(m) * rs_c(m) * c1:.7e}, Theta {e_th:.1e}, dr(S) {e_S:.1e}, dM cuts {e_m:.1e}")
            ok = ok and good
    check("T8", "#248's classes A0, A1R, A1L (and G with the scalar) at each member, Richardson of h = 0.00125, 0.000625: Theta_gen(s) = -r_S c1 at s = -0.03 ... 1.9, dr(S) = c1 of the horizon fit, kappa r_S c1 = dM_ADM of the cut extensions (s_cut = 1.0, 1.9, -0.03, read at rho = 1e8), all <= 1e-6 relative to max(|c1|, 1e-5)",
          ok, "; ".join(det))


CHECKS = [("T1", check_T1), ("T2", check_T2), ("T3", check_T3), ("T4", check_T4), ("T5", check_T5), ("T6", check_T6), ("T7", check_T7), ("T8", check_T8)]
SCOPE = tuple(c for c, _ in CHECKS)
SCOPE_D = ("S1", "S2", "S3")


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
DECLARED = {                      # computed in the development run d1 (raw I301-mac/dev) before the runs of record; every check is in scope
    "s_q": {"S1"},
    "s_pot": {"S2"},
    "s_k": {"S3"},
    "rs_one": {"T1", "T4", "T7", "T8"},
    "ctl_true": {"T5"},
    "c1_noshift": {"T1"},
    "th_noF": {"T4", "T6", "T7", "T8"},
    "th_nophi": {"T4", "T6", "T7", "T8"},
    "th_asign": {"T3", "T4", "T6", "T7", "T8"},
    "dm_noT": {"T4", "T7", "T8"},
    "kappa_off": {"T1", "T4", "T5", "T7", "T8"},
    "jump_off": {"T4", "T7", "T8"},
    "xform_p": {"T3", "T4", "T7", "T8"},
    "xform_q": {"T3", "T4", "T7", "T8"},
    "rho_ellis": {"T3", "T4", "T7", "T8"},
    "mode_khat": {"T4", "T7"},
    "mode_dF0": {"T4", "T7"},
    "law_norm": {"T2"},
    "law_mass": {"T2"},
    "gauge_noRp": {"T6"},
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
    print("== members (mpmath references; r_S = R(rho_h), kappa = A'(rho_h)/2 for the Killing field with A -> 1, G = 1, M = Misner-Sharp mass at infinity) ==")
    print(f"{'member':>14} {'rho_h':>12} {'kappa':>12} {'r_S':>12} {'Rp(rho_h)':>10} {'c1 per unit':>13} {'dM per unit':>13} {'kappa r_S c1':>13} {'kappa c1':>13} {'-c/pi':>13}")
    for tag in TAGS:
        m = member(tag)
        c1x, dMx = exact_ref(m)
        print(f"{m.tag:>14} {m.rh:12.9f} {m.kappa:12.9f} {m.rS:12.9f} {m.Rph:10.7f} {c1x:13.10f} {dMx:13.10f} {m.kappa * m.rS * c1x:13.10f} {m.kappa * c1x:13.10f} {(f'{-m.c / math.pi:13.10f}' if m.c is not None else '            -')}")
    for tag in TAGS:
        m = member(tag)
        md, bg, sl = mode_of(m), bg_of(m), slice_of(m)
        c1x, dMx = exact_ref(m)
        print(f"== exact mode at {m.tag}: Theta_gen + r_S c1 (per unit dL or d r_S), dm of the slice code against the exact d m at fixed Pi ==")
        print(f"{'s':>9} {'rho':>10} {'Theta_gen+r_S c1':>17} {'Theta (no dF)':>15} {'dm slice':>15} {'dm exact':>15}")
        for s in (S_MODE_SCHW if m.kind == "schw" else S_MODE):
            dr, dr_P, dphi, dF = md.fields(s * s, bg.par, bg.tab)
            drs = 2 * s * dr_P
            mu, a, cc, Q, b = sl.theta_coef(s)
            rho, q = md.rho_q(s * s, bg.par, bg.tab)
            print(f"{s:9.3g} {rho:10.4f} {theta_gen(sl, s, dr, drs, dphi, dF) + m.rS * c1x:17.3e} {a * dr + mu * drs + cc * dphi:15.6e} {dm_s(bg, s, dr, drs, dF):15.10f} {md.dm_exact(rho, q):15.10f}")
        cuts = []
        for sc in (S_CUT_MODE_SCHW if m.kind == "schw" else S_CUT_MODE):
            dr, dr_P, dphi, dF = md.fields(sc * sc, bg.par, bg.tab)
            cuts.append((sc, read_mass(m, sc, dr, 2 * sc * dr_P, dphi, dF)["dM"]))
        print("   cut extensions, dM at rho = 1e8 (kappa r_S c1 = %.12f): " % (m.kappa * m.rS * c1x) + ", ".join(f"{sc:g}: {v:.12f}" for sc, v in cuts))
        if m.kind == "ellis":
            for rho in (30.0, 1e2, 1e4, 1e8):
                q = md.q_big(rho) if rho > 40 else md.q_qp(rho)[0]
                d = md.dm_exact(rho, q)
                print(f"   rho = {rho:8.1e}: dm = {d:.12f}, (dm - dM) rho = {(d - dMx) * rho:.6f}")
    print("== the evolved exact mode (amplitude 2e-3), Richardson of h = 0.00125, 0.000625 ==")
    for tag in TAGS:
        m = member(tag)
        rs_ = {h: mode_run(m, h, 2e-3) for h in HL}
        h1, h2 = HL[-2], HL[-1]
        c1x, dMx = exact_ref(m)
        print(f"{m.tag:>14}: c1 fit {rich(rs_[h1]['c1fit'], rs_[h2]['c1fit']):.9e} (exact {2e-3 * c1x:.9e}); dr(S) {rich(rs_[h1]['c1S'], rs_[h2]['c1S']):.9e}; c2 {rich(rs_[h1]['c2fit'], rs_[h2]['c2fit']):.2e}; dM cuts "
              + ", ".join(f"{rich(rs_[h1]['cut'][sc]['dM'], rs_[h2]['cut'][sc]['dM']):.9e}" for sc in S_CUTS + S_CUTL) + f" (exact {2e-3 * dMx:.9e}); error/scale " + ", ".join(sci(rs_[h]["err"] / rs_[h]["scale"]) for h in HL))
    print("== #248's classes (corner (-0.05, -2) of the rectangle, amplitude 1e-3), Richardson of h = 0.00125, 0.000625 ==")
    print(f"{'member':>14} {'class':>4} {'c1 fit':>14} {'dr(S)':>14} {'c2':>14} {'kappa r_S c1':>14} {'dM cut 1.0':>14} {'dM cut 1.9':>14} {'dM cut -0.03':>14} {'Theta(0)+r_S c1':>15} {'Theta(1.9)+r_S c1':>16} {'kappa c1':>14} {'fit res':>9}")
    for tag in TAGS:
        m = member(tag)
        for name in CLASSES_T8:
            if m.kind == "schw" and name == "G":
                continue
            sm = class_summary(m, name)
            c1 = sm["c1"]
            print(f"{m.tag:>14} {name:>4} {c1:14.7e} {sm['c1S']:14.7e} {sm['c2']:14.7e} {m.kappa * m.rS * c1:14.7e} {sm['dM'][1.0]:14.7e} {sm['dM'][1.9]:14.7e} {sm['dM'][-0.03]:14.7e} {sm['theta'][0.0] + m.rS * c1:15.2e} {sm['theta'][1.9] + m.rS * c1:16.2e} {m.kappa * c1:14.7e} {sm['res']:9.1e}")
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

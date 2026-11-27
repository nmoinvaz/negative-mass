#!/usr/bin/env python3
"""Radial (l = 0) modes of the sourced black bounce (Simpson-Visser metric with the phantom scalar and magnetic NED field of
[BW22], arXiv:2112.13198) and of the magnetic black universes of [BBS12] (arXiv:1208.4619). Method of #212 (stability.py).
See README.md in this directory, section "Sourced black bounce and magnetic black universes (#223)".

    bounce.py check    [--mutate=M]     background, linearization with the NED/Maxwell term, family coefficients
    bounce.py planted  [--res=R] [--mutate=M]   P1: Ellis wormhole through the S family; W^2 a^2 against w2_output.txt
    bounce.py planted_low [--res=R] [--mutate=M]   P3: a weakly bound state (W = 3e-5) of a known potential through the low-W stage's machinery
    bounce.py count ID [--res=R1|R2] [--stage=main|low] [--out=FILE.json] [--mutate=M]    growing modes of one member (ID as listed by `members`);
                                       main: W in [1e-3, 8]; low: W in [1e-5, 1e-3]; the Sturm count (x_h >= 0) is included
    bounce.py members                      the declared members
    bounce.py report DIR                   table, P2 (zero charge against #212), R1 v R2, low stage, Sturm, from DIR/*.json
    bounce.py spectral DIR                 independent Chebyshev check of the counts in DIR/*_R2.json
    bounce.py controls [--base=DIR] [--only=M1,M2]   runs every --mutate=M and asserts that its target fires (DIR: finished R2 results to reuse)

Units: b = min R = 1 (a = 1 for S).  Exit status 1 when a check fails (a mutated run is meant to exit 1).
"""
import glob
import json
import math
import os
import re
import subprocess
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stability as S  # noqa: E402  (run_chunks, make_ev, cheb: the generic integrator helpers of #212)

# ------------------------------------------------------------------------------------------------ declared settings
PLANTED_TOL = {"R1": 1e-7, "R2": 1e-9}                     # absolute on W^2 a^2, as #212
TABLE_TOL = {"R1": 1e-7, "R2": 2e-9}                       # relative on W, zero charge against #212's table (10 printed decimals)
TABLE_FILE = HERE + "/stability_output.txt"
W_LO, W_HI = 1e-3, 8.0                                     # range searched, units 1/b
RES = {
    "R1": dict(rtol=1e-9, z0=1e-9, rmax0=150.0, rmaxW=30.0, delta=0.3, nW=80),
    "R2": dict(rtol=1e-12, z0=1e-13, rmax0=300.0, rmaxW=45.0, delta=0.2, nW=200),
}
ROOT_AGREE = 1e-6                                          # relative, R1 against R2
# low-W stage (amendment A2, ruling of the coordinator): shooting on [1e-5, 1e-3] with its own outer radius r0 = max(rmax0, rmaxL/W)
W_LOW_LO, W_LOW_HI = 1e-5, 1e-3
LOW = {"R1": dict(nW=40, rmaxL=20.0), "R2": dict(nW=100, rmaxL=20.0)}
ROOT_AGREE_LOW = 1e-4                                      # relative, low stage R1 against R2
PT_TOL = {"R1": 1e-3, "R2": 1e-4}                          # relative on W, planted Poeschl-Teller weak bound state
PT_KAPPA, PT_NU = 1.0, 3e-5                                # V = -nu (nu + 1) kappa^2 sech^2(kappa x): one bound state, W = nu kappa
LAM_MIN = W_LO**2                                          # counted real eigenvalues of the spectral check: lambda > W_LO^2
SPEC = {"S-A": dict(N=401, X0=45.0, kX=12.0, tol_nothroat=1e-2, tol_throat=1e-2),
        "S-B": dict(N=801, X0=90.0, kX=24.0, tol_nothroat=1e-4, tol_throat=2e-2)}
SER_X = 30.0                                                # amendment A7: low stage, B family, A and A' from the series in 1/x for x >= SER_X
MAP_X, MAP_A = 150.0, 5.0                                  # amendment A4: spectral grid mapped (sinh) where X > MAP_X
MUT = set()


def mutating(name):
    return name in MUT


# ------------------------------------------------------------------------------------------------ members
S_M = [0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]
B1_M = 2 / math.pi
B1_Q = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.74, 0.76, 0.77, 0.775]
B2_Q = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.78, 0.8, 0.81]
B3_M = 0.2 / math.pi
B3_Q = [0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.2, 0.22, 0.23, 0.24]
P2_C = [-10.0, -1.1, -0.9, -0.01]                          # the zero-charge members that are not the q = 0 members of B1 (c = -3), B2 (-1), B3 (-0.3)
P2_ALIAS = {"B1_q0": -3.0, "B2_q0": -1.0, "B3_q0": -0.3}   # these three are compared with #212's table too


def m_throat(q):
    """mass at which the horizon of B sits at the throat, A(0) = 0:  m = (2/3pi)(1 + q^2 + pi^2 q^2 / 4)"""
    return (2 / (3 * math.pi)) * (1 + q * q + math.pi**2 * q * q / 4)


def member_ids():
    ids = ["S_m%g" % m for m in S_M]
    ids += ["B1_q%g" % q for q in B1_Q] + ["B2_q%g" % q for q in B2_Q] + ["B3_q%g" % q for q in B3_Q]
    ids += ["P2_c%g" % c for c in P2_C]
    return ids


def member_spec(mid):
    s, v = mid.split("_")
    if s == "S":
        return dict(fam="sv", m=float(v[1:]), q=0.0, at_throat=False)
    if s == "B1":
        return dict(fam="bbs", m=B1_M, q=float(v[1:]), at_throat=False)
    if s == "B2":
        q = float(v[1:])
        return dict(fam="bbs", m=m_throat(q), q=q, at_throat=True)
    if s == "B3":
        return dict(fam="bbs", m=B3_M, q=float(v[1:]), at_throat=False)
    if s == "P2":
        c = float(v[1:])
        return dict(fam="bbs", m=-2 * c / (3 * math.pi), q=0.0, at_throat=(c == -1.0))
    raise ValueError(mid)


def planted_value():
    return S.planted_value()


# ------------------------------------------------------------------------------------------------ background (numerics)
def mp_A(fam, m, q):
    """closed-form A(x) in mpmath (b = 1)"""
    import mpmath as mp
    m = mp.mpf(m); q = mp.mpf(q)
    if fam == "sv":
        return lambda x: 1 - 2 * m / mp.sqrt(x * x + 1)
    p = 3 * m - mp.pi * q * q
    B0 = -mp.pi * p / 2 - mp.pi**2 * q * q / 4

    def A(x):
        B = B0 + (1 + q * q + p * x) / (1 + x * x) + (p + 2 * q * q * x / (1 + x * x)) * mp.atan(x) + q * q * mp.atan(x)**2
        return (1 + x * x) * B
    return A


class Bg:
    """static background, b = 1: A, A', R L_R (R times dLambda/dR of the NED energy density), the outermost horizon x_h.
    kind: 'twoend' (no horizon), 'simple' (A'(x_h) != 0), 'ext0' (S, a = 2m: double zero of A at x = 0)."""

    def __init__(self, fam, m, q=0.0, res="R2", at_throat=False, stage="main"):
        self.fam, self.m, self.q = fam, float(m), float(q)
        self.res = dict(RES[res])
        if stage == "low":
            self.res.update(rmaxW=LOW[res]["rmaxL"], nW=LOW[res]["nW"])
            if mutating("rmax_low"):
                self.res.update(rmax0=4.0, rmaxW=0.0)
        if mutating("rmax"):
            self.res.update(rmax0=4.0, rmaxW=0.0)
        if mutating("hz_start"):
            self.res.update(z0=0.05)
        self.ned = 0.0 if mutating("ned_off") else 1.0
        self.xh = self.A1 = None
        self.taylor = None
        self.series = None
        if fam == "sv":
            self.xh2 = 4 * self.m**2 - 1
            if self.m > 0.5:
                self.xh = math.sqrt(self.xh2); self.kind = "simple"
            elif self.m == 0.5:
                self.xh = 0.0; self.kind = "ext0"
            else:
                self.kind = "twoend"
            if self.xh is not None:
                self.A1 = float(self.Ap(self.xh))
        else:
            self.p = 3 * self.m - math.pi * self.q**2
            self.kind = "simple"
            if stage == "low":
                self._series_setup()
            self.xh = 0.0 if at_throat else self._horizon()
            self.A1 = 2 * self.p if at_throat else float(self.Ap(self.xh))
            self._taylor_setup(at_throat)

    # -- A, A', R L_R, for real and complex x
    def _B(self, x):
        q2, p = self.q**2, self.p
        x = np.asarray(x)
        if np.iscomplexobj(x):
            at = np.arctan(x)
            return -math.pi * p / 2 - math.pi**2 * q2 / 4 + (1 + q2 + p * x) / (1 + x * x) + (p + 2 * q2 * x / (1 + x * x)) * at + q2 * at**2
        xs = np.where(x >= 1.0, x, 2.0)                                    # t = arctan(1/x) form, no cancellation for x >= 1
        t = np.arctan(1 / xs)
        Bt = ((1 + q2 + p * xs) / (1 + xs * xs) + math.pi * q2 * xs / (1 + xs * xs) - p * t - 2 * q2 * xs * t / (1 + xs * xs)
              - math.pi * q2 * t + q2 * t * t)
        at = np.arctan(x)
        Bo = -math.pi * p / 2 - math.pi**2 * q2 / 4 + (1 + q2 + p * x) / (1 + x * x) + (p + 2 * q2 * x / (1 + x * x)) * at + q2 * at**2
        return np.where(x >= 1.0, Bt, Bo)

    def _series_setup(self, K=40):
        """B = sum_k b_k y^k, y = 1/x (amendment A7): B(x) = (1+q^2+px)/(1+x^2) + pi q^2 x/(1+x^2) - (p + pi q^2) t - 2 q^2 x t/(1+x^2) + q^2 t^2,
        t = atan(y); no cancellation at large x (the closed form loses 1e-16 x relative accuracy in A - 1 through B ~ 1/x^2 from O(p/x) terms)."""
        q2, p = self.q**2, self.p
        S = np.zeros(K); S[0::2] = (-1.0) ** np.arange(len(S[0::2]))                       # 1/(1+y^2)
        t = np.zeros(K); t[1::2] = (-1.0) ** np.arange(len(t[1::2])) / np.arange(1, K, 2)   # atan(y)
        mul = lambda a, b: np.convolve(a, b)[:K]
        yS = np.concatenate([[0.0], S[:-1]])                                              # y/(1+y^2)
        y2S = np.concatenate([[0.0, 0.0], S[:-2]])                                        # y^2/(1+y^2)
        b = (1 + q2) * y2S + p * yS + math.pi * q2 * yS - (p + math.pi * q2) * t - 2 * q2 * mul(yS, t) + q2 * mul(t, t)
        if mutating("ser_coef"):
            b[3] *= 1.001                                                                 # control: wrong coefficient of the series
        b[0] = b[1] = 0.0                                                                 # vanish analytically (B = O(y^2)); exact zeros, not roundoff
        c = np.zeros(K)                                                                   # A = (1 + y^-2) B = sum_n c_n y^n, c_n = b_{n+2} + b_n
        c[:K - 2] = b[2:]
        c += b
        self.ser_c = c
        self.series = True

    def _A_series(self, x):
        y = 1.0 / np.asarray(x, dtype=float)
        A = np.polynomial.polynomial.polyval(y, self.ser_c)
        n = np.arange(1, len(self.ser_c))
        Ap = -y**2 * np.polynomial.polynomial.polyval(y, self.ser_c[1:] * n)               # dA/dx = -y^2 dA/dy
        return A, Ap

    def A(self, x):
        if self.series and np.isrealobj(x) and np.all(np.asarray(x) >= SER_X):
            return self._A_series(x)[0]
        if self.fam == "sv":
            R = np.sqrt(x * x + 1)
            return (x * x - self.xh2) / (R * (R + 2 * self.m))
        return (x * x + 1) * self._B(x)

    def Ap(self, x):
        if self.series and np.isrealobj(x) and np.all(np.asarray(x) >= SER_X):
            return self._A_series(x)[1]
        if self.fam == "sv":
            return 2 * self.m * x / np.sqrt(x * x + 1)**3
        q2, p = self.q**2, self.p
        Bp = 2 / (1 + x * x)**2 * (p - x + 2 * q2 * np.arctan(x))
        return 2 * x * self._B(x) + (x * x + 1) * Bp

    def RL(self, x):
        """R * dLambda/dR; S: Lambda = 6m/(5 R^5); B: Lambda = q^2/R^4"""
        R = np.sqrt(x * x + 1)
        return -6 * self.m / R**5 if self.fam == "sv" else -4 * self.q**2 / R**4

    def Vhat(self, x):
        """u_xx = (Vhat + W^2) u, u = Y/x, dx_* = dx/A"""
        A = self.A(x); Ap = self.Ap(x); R2 = x * x + 1
        return 2 * A * A / x**2 - (1 + 1 / R2) * A * Ap / x + A * A / R2**2 - 2 * A / R2 - self.ned * A * self.RL(x)

    # -- horizon
    def _horizon(self):
        xs = np.linspace(60.0, -40.0, 400001)
        ys = self.A(xs)
        for i in range(len(xs) - 1):
            if ys[i] > 0 and ys[i + 1] <= 0:
                r = brentq(lambda t: float(self.A(t)), xs[i + 1], xs[i], xtol=1e-15, rtol=1e-15)
                return r
        raise RuntimeError("no horizon found")

    def _taylor_setup(self, at_throat):
        import mpmath as mp
        mp.mp.dps = 60
        A = mp_A("bbs", mp.mpf(self.m), mp.mpf(self.q))
        xh = mp.mpf(0) if at_throat else mp.findroot(A, mp.mpf(self.xh))
        if not at_throat:
            self.xh = float(xh)
        self.taylor = [float(c) for c in mp.taylor(A, xh, 14)]
        if mutating("az_series"):
            self.taylor[2] *= 1.01                                      # control: wrong second Taylor coefficient of A
        if at_throat:
            self.taylor[0] = 0.0
        self.A1 = self.taylor[1]

    def a_over_z(self, z):
        """A(x_h + z)/z without cancellation"""
        if self.fam == "sv":
            R = np.sqrt((self.xh + z)**2 + 1)
            return (2 * self.xh + z) / (R * (R + 2 * self.m))
        if abs(z) < 0.05:
            return sum(c * z**(k - 1) for k, c in enumerate(self.taylor) if k >= 1)
        return self.A(self.xh + z) / z


def coef_x(bg, r, A, Ar, Ap):
    """Y_xx = Cx Y_x + Dx Y + W^2 Y (dx = d x_*), Ar = A / x. eps = -1 (phantom); the control `eps` flips the kinetic sign in the
    generic form (Cx = (-2 eps A/x + 2 A x)/R^2, Dx = -eps A^2/R^4 + eps Ar A'/R^2 - 2A/R^2 - A R Lam_R)"""
    R2 = r * r + 1
    if mutating("eps"):
        e = 1.0
        return (-2 * e * Ar + 2 * A * r) / R2, -e * A * A / R2**2 + e * Ar * Ap / R2 - 2 * A / R2 - bg.ned * A * bg.RL(r)
    return 2 * Ar, A * A / R2**2 - Ar * Ap / R2 - 2 * A / R2 - bg.ned * A * bg.RL(r)


def rhs_x(bg, W):
    def f(x, s):
        r, Y, Yx = s
        A = bg.A(r); Ap = bg.Ap(r)
        Cx, Dx = coef_x(bg, r, A, A / r, Ap)
        return [A, Yx, Cx * Yx + Dx * Y + W * W * Y]
    return f


def rhs_zeta(bg, W):
    def f(x, s):
        zeta, Y, Yx = s
        z = np.exp(zeta); r = bg.xh + z
        at = bg.a_over_z(z); A = z * at
        Ar = at if bg.xh == 0.0 else A / r
        Cx, Dx = coef_x(bg, r, A, Ar, bg.Ap(r))
        return [at, Yx, Cx * Yx + Dx * Y + W * W * Y]
    return f


def rhs_arc(bg, W, d):
    """x-form along x = d exp(i theta): Y'' = P1 Y' + P0 Y + W^2 Y / A^2"""
    e = 1.0 if mutating("eps") else -1.0

    def f(th, s):
        Y, Yr = s
        r = d * np.exp(1j * th)
        A = bg.A(r); Ap = bg.Ap(r); R2 = r * r + 1
        P1 = -2 * e / (r * R2) + 2 * r / R2 - Ap / A
        P0 = -e / R2**2 + e * Ap / (A * r * R2) - 2 / (A * R2) - bg.ned * bg.RL(r) / A
        return [1j * r * Yr, 1j * r * (P1 * Yr + P0 * Y + W * W * Y / A**2)]
    return f


def match_point(bg):
    return 1.0 if (bg.xh is None or bg.xh < 0.5) else bg.xh + 1.0


def has_throat(bg):
    return bg.kind == "twoend" or (bg.xh is not None and bg.xh < 0.0)


def left_solution(bg, W, rm, nodes=False):
    """(x, Y, Y_x, nodes, flag) at x = rm of the solution ingoing at the horizon (or decaying at x -> -inf), continued through the throat"""
    p = bg.res; rtol = p["rtol"]
    node_ev = S.make_ev(lambda x, s: s[1], 0, False) if nodes else None
    n_nodes = 0
    atol = [1e-12, 1e-30, 1e-30]
    if mutating("throat_wall") and has_throat(bg):
        # control: independent sides, Dirichlet wall at the throat (Y ~ x^3), right side only
        r0 = 1e-4
        s = np.array([r0, r0**3, bg.A(r0) * 3 * r0**2]); s[1:] /= np.hypot(s[1], s[2])
        stop = S.make_ev(lambda x, ss: ss[0] - rm, 1, True)
        s, ln, nn = S.run_chunks(rhs_x(bg, W), s, +1, W, rtol, atol, stop, node_ev)
        return s[0], s[1], s[2], nn, None
    if bg.kind == "twoend":
        r0 = -max(p["rmax0"], p["rmaxW"] / W)
        lg = W + bg.Vhat(r0) / (2 * W)
        s = np.array([r0, 1.0, bg.A(r0) / r0 + lg]); s[1:] /= np.hypot(s[1], s[2])
        stop = S.make_ev(lambda x, ss: ss[0] + p["delta"], 1, True)
        s, ln, nn = S.run_chunks(rhs_x(bg, W), s, +1, W, rtol, atol, stop, node_ev)
        n_nodes += nn
        rstart = s[0]
    elif bg.kind == "ext0":
        x0 = float(np.clip(0.04 * W, 1e-4, 0.1))
        lg = W + bg.Vhat(x0) / (2 * W)
        s = np.array([x0, 1.0, bg.A(x0) / x0 + lg]); s[1:] /= np.hypot(s[1], s[2])
        stop = S.make_ev(lambda x, ss: ss[0] - rm, 1, True)
        s, ln, nn = S.run_chunks(rhs_x(bg, W), s, +1, W, rtol, atol, stop, node_ev)
        return s[0], s[1], s[2], nn, None
    else:
        if bg.xh == 0.0:
            s = np.array([math.log(p["z0"]), 1.0, bg.A1 + W])          # Y ~ x^(1 + W/A1)
        else:
            s = np.array([math.log(p["z0"]), 1.0, W])                  # Y ~ z^(W/A1)
        s[1:] /= np.hypot(s[1], s[2])
        rstop = rm if bg.xh >= 0 else -min(p["delta"], 0.5 * abs(bg.xh))
        stop = S.make_ev(lambda x, ss: bg.xh + np.exp(ss[0]) - rstop, 1, True)
        s, ln, nn = S.run_chunks(rhs_zeta(bg, W), s, +1, W, rtol, atol, stop, node_ev)
        n_nodes += nn
        if bg.xh >= 0:
            return bg.xh + np.exp(s[0]), s[1], s[2], n_nodes, None
        s = np.array([bg.xh + np.exp(s[0]), s[1], s[2]])
        rstart = s[0]
    d = -rstart
    Yl, Yxl = s[1], s[2]
    sol = solve_ivp(lambda th, y: rhs_arc(bg, W, d)(th, y), [np.pi, 0.0], [Yl + 0j, Yxl / bg.A(rstart) + 0j], method="DOP853",
                    rtol=rtol, atol=1e-30)
    Yc, Yrc = sol.y[:, -1]
    flag = (np.sign(Yl) == np.sign(Yc.real), abs(Yc.imag) / abs(Yc))
    s = np.array([d, Yc.real, (Yrc * bg.A(d)).real]); s[1:] /= np.hypot(s[1], s[2])
    stop = S.make_ev(lambda x, ss: ss[0] - rm, 1, True)
    s, ln, nn = S.run_chunks(rhs_x(bg, W), s, +1, W, rtol, atol, stop, node_ev)
    return s[0], s[1], s[2], n_nodes + nn, flag


def right_solution(bg, W, rm, nodes=False):
    p = bg.res; rtol = p["rtol"]
    r0 = max(p["rmax0"], p["rmaxW"] / W)
    lg = -W - bg.Vhat(r0) / (2 * W)
    s = np.array([r0, 1.0, bg.A(r0) / r0 + lg]); s[1:] /= np.hypot(s[1], s[2])
    node_ev = S.make_ev(lambda x, ss: ss[1], 0, False) if nodes else None
    stop = S.make_ev(lambda x, ss: ss[0] - rm, -1, True)
    s, ln, nn = S.run_chunks(rhs_x(bg, W), s, -1, W, rtol, [1e-12, 1e-30, 1e-30], stop, node_ev)
    return s[0], s[1], s[2], nn


def mismatch(bg, W, nodes=False):
    rm = match_point(bg)
    _, Yl, Ylx, nl, flag = left_solution(bg, W, rm, nodes)
    _, Yr, Yrx, nr = right_solution(bg, W, rm, nodes)
    D = Yl * Yrx - Ylx * Yr
    return (D, nl + nr, flag) if nodes else D


def sturm_count(bg):
    """negative eigenvalues of -d^2/dx_*^2 + Vhat on the static region (x_h >= 0): nodes of the zero-energy solution that is regular
    (bounded) at the horizon, plus one if it will still cross zero in the linear tail u ~ a + b x_*. Integrated in the Y form at W = 0
    (u = Y/x, x > 0, same nodes), which has no cancellation at x_h = 0 where Vhat -> 0 by cancelling O(1) terms."""
    if bg.xh is None or bg.xh < 0:
        return None
    z0 = 1e-8 if bg.kind != "ext0" else 1e-4
    r_end = 2e4
    f0 = rhs_zeta(bg, 0.0)
    flip = -1.0 if mutating("sturm_sign") else 1.0                     # control: sign of the potential term

    def f(x, s):
        d = f0(x, s)
        zeta, Y, Yx = s
        z = np.exp(zeta); r = bg.xh + z
        A = z * bg.a_over_z(z)
        Cx, Dx = coef_x(bg, r, A, bg.a_over_z(z) if bg.xh == 0.0 else A / r, bg.Ap(r))
        return [d[0], Yx, Cx * Yx + flip * Dx * Y]
    if bg.xh == 0.0:
        s = np.array([math.log(z0), 1.0, float(bg.a_over_z(z0))])         # Y ~ x: u = const
    else:
        s = np.array([math.log(z0), 1.0, 0.0])                            # Y ~ const at a simple horizon
    s[1:] /= np.hypot(s[1], s[2])
    nodes = 0
    stop = S.make_ev(lambda x, ss: bg.xh + np.exp(ss[0]) - r_end, 1, True)
    zero = S.make_ev(lambda x, ss: ss[1], 0, False)
    x = 0.0
    while True:
        sol = solve_ivp(f, [x, x + 50.0], s, method="DOP853", rtol=1e-12, atol=[1e-12, 1e-30, 1e-30], events=[stop, zero])
        s = sol.y[:, -1].copy(); x = sol.t[-1]
        nodes += len(sol.t_events[1])
        n = np.hypot(s[1], s[2]); s[1:] /= n
        if sol.status == 1:
            break
    z = np.exp(s[0]); r = bg.xh + z
    A = z * bg.a_over_z(z)
    u = s[1] / r
    ux = s[2] / r - s[1] * A / r**2
    if u * ux < 0:
        nodes += 1
    return nodes


def count_member(spec, res="R2", W_lo=W_LO, W_hi=W_HI, nW=None, stage="main"):
    if stage == "low":
        W_lo, W_hi = W_LOW_LO, W_LOW_HI
    bg = Bg(spec["fam"], spec["m"], spec["q"], res, spec["at_throat"], stage)
    if mutating("scan_range"):
        W_lo = 3.0
    p = bg.res
    nW = nW or p["nW"]
    Ws = np.exp(np.linspace(math.log(W_lo), math.log(W_hi), nW))
    t0 = time.time()
    D = np.array([mismatch(bg, W) for W in Ws])
    roots = []
    for i in range(len(Ws) - 1):
        if D[i] == 0 or np.sign(D[i]) != np.sign(D[i + 1]):
            roots.append(brentq(lambda W: mismatch(bg, W), Ws[i], Ws[i + 1], xtol=1e-14, rtol=1e-13))
    info = []
    for W in roots:
        Dn, nn, flag = mismatch(bg, W, nodes=True)
        info.append(dict(W=W, W2=W * W, nodes=int(nn), arc_same_sign=None if flag is None else bool(flag[0]),
                         arc_imag=None if flag is None else float(flag[1])))
    st = sturm_count(bg)
    return dict(fam=bg.fam, m=bg.m, q=bg.q, at_throat=spec["at_throat"], res=res, kind=bg.kind, xh=bg.xh, A1=bg.A1,
                kappa=None if bg.A1 is None else bg.A1 / 2, p=getattr(bg, "p", None), stage=stage, rmaxW=p["rmaxW"], n_roots=len(roots), roots=info, sturm=st,
                D_min_abs=float(np.min(np.abs(D))), D_first=float(D[0]), D_last=float(D[-1]),
                sign_changes=int(np.sum(np.sign(D[:-1]) != np.sign(D[1:]))), W_lo=W_lo, W_hi=W_hi, nW=nW,
                seconds=time.time() - t0, mutate=sorted(MUT))


def planted_command(res_list=("R1", "R2")):
    ok_all = True
    ref = planted_value()
    import mpmath as mp
    mp.mp.dps = 30
    for res in res_list:
        o = count_member(member_spec("S_m0"), res)
        n = o["n_roots"]
        d = abs(mp.mpf(o["roots"][0]["W2"]) - mp.mpf(ref)) if n else None
        ok = (n == 1) and d is not None and d <= PLANTED_TOL[res]
        ok_all &= ok
        if n:
            print(("PASS  " if ok else "FAIL  ") + f"P1 planted Ellis (S, m = 0) {res}: growing modes {n} (declared 1); W^2 a^2 = "
                  f"{o['roots'][0]['W2']:.15f}; reference {ref}; |difference| = {float(d):.3e} (declared <= {PLANTED_TOL[res]:g}); "
                  f"nodes {o['roots'][0]['nodes']}", flush=True)
        else:
            print(f"FAIL  P1 planted Ellis (S, m = 0) {res}: growing modes 0 (declared 1); no growing mode found", flush=True)
    return ok_all


def planted_low_command(res_list=("R1", "R2")):
    """P3: a weakly bound state of a known potential through the low-W stage's machinery (outer radius rule r0 = max(rmax0, rmaxL/W), chunked
    integration with renormalization, log grid, sign change of the Wronskian mismatch, Brent): u'' = (V + W^2) u,
    V = -nu (nu + 1) kappa^2 sech^2(kappa x); the single bound state has W = nu kappa = 3e-5. Shares run_chunks and the scan with the
    runs; not the background, the horizon or the throat code."""
    ok_all = True
    nu = PT_NU * (2.0 if mutating("lowW_V") else 1.0)                    # control: wrong potential depth
    kap = PT_KAPPA
    ref = PT_NU * PT_KAPPA
    V = lambda x: -nu * (nu + 1) * kap**2 / np.cosh(np.clip(kap * x, -300, 300))**2
    xw = 40.0 + 1000.0                                                   # the well is |x| < 40 (V < 1e-34 beyond): free flow with run_chunks up to |x| = xw, max_step inside
    for res in res_list:
        p = dict(RES[res]); rtol = p["rtol"]; nW = LOW[res]["nW"]; rmaxL = LOW[res]["rmaxL"]
        lo = 1e-4 if mutating("lowW_range") else W_LOW_LO                   # control: a range that excludes the state

        def mism(W):
            r0 = max(p["rmax0"], rmaxL / W)
            f = lambda x, s: [1.0, s[2], (V(s[0]) + W * W) * s[1]]
            sl = np.array([-r0, 1.0, W]); sl[1:] /= np.hypot(sl[1], sl[2])
            sl, _, _ = S.run_chunks(f, sl, +1, W, rtol, [1e-12, 1e-30, 1e-30], S.make_ev(lambda x, ss: ss[0] + xw, 1, True))
            sol = solve_ivp(f, [sl[0], 0.0], sl, method="DOP853", rtol=rtol, atol=[1e-12, 1e-30, 1e-30], max_step=0.5)
            sl = sol.y[:, -1]
            sr = np.array([r0, 1.0, -W]); sr[1:] /= np.hypot(sr[1], sr[2])
            sr, _, _ = S.run_chunks(f, sr, -1, W, rtol, [1e-12, 1e-30, 1e-30], S.make_ev(lambda x, ss: ss[0] - xw, -1, True))
            sol = solve_ivp(f, [sr[0], 0.0], sr, method="DOP853", rtol=rtol, atol=[1e-12, 1e-30, 1e-30], max_step=0.5)
            sr = sol.y[:, -1]
            return sl[1] * sr[2] - sl[2] * sr[1]
        Ws = np.exp(np.linspace(math.log(lo), math.log(W_LOW_HI), nW))
        D = np.array([mism(W) for W in Ws])
        roots = [brentq(mism, Ws[i], Ws[i + 1], xtol=1e-16, rtol=1e-13) for i in range(nW - 1) if np.sign(D[i]) != np.sign(D[i + 1])]
        ok = len(roots) == 1 and abs(roots[0] - ref) / ref <= PT_TOL[res]
        ok_all &= ok
        print(("PASS  " if ok else "FAIL  ") + f"P3 planted weak bound state {res}: roots {len(roots)} (declared 1), W = "
              f"{roots[0] if roots else float('nan'):.10e}, reference {ref:.10e}, relative difference "
              f"{abs(roots[0] - ref) / ref if roots else float('nan'):.3e} (declared <= {PT_TOL[res]:g}); range [{lo:g}, {W_LOW_HI:g}], {nW} points", flush=True)
    return ok_all


# ------------------------------------------------------------------------------------------------ symbolic checks
def sym_funcs(fam, mut=()):
    """closed forms in sympy (b = 1): returns dict with x, m, q, A, V, Lam (Lambda as a function of x), R, P0"""
    import sympy as sp
    x = sp.Symbol("x", real=True)
    m, q = sp.symbols("m q", positive=True)
    R = sp.sqrt(x**2 + 1)
    P0 = sp.sqrt(2) * sp.atan(x)
    if fam == "sv":
        cV = sp.Rational(3, 5) if "sv_V" in mut else sp.Rational(4, 5)             # [BW22] (22), Q = 0: V = 4 m cos^5(phi)/5
        cL = sp.Rational(4, 5) if "sv_L" in mut else sp.Rational(6, 5)             # [BW22] (21), Q = 0: Lambda = L/2 = 6 m/(5 R^5)
        A = 1 - 2 * m / R
        V = cV * m / R**5
        Lam = cL * m / R**5
    else:
        p = 3 * m - sp.pi * q**2
        B0 = -sp.pi * p / 2 - sp.pi**2 * q**2 / 4
        k = 2 if "bbs_A" in mut else 1                                              # control: wrong coefficient of the arctan^2 term
        B = B0 + (1 + q**2 + p * x) / (1 + x**2) + (p + 2 * q**2 * x / (1 + x**2)) * sp.atan(x) + k * q**2 * sp.atan(x)**2
        A = (1 + x**2) * B
        V = -sp.diff(sp.diff(A, x) * (1 + x**2), x) / (2 * (1 + x**2)) + q**2 / (1 + x**2)**2        # [BBS12] (9)
        Lam = q**2 / R**4
    return dict(sp=sp, x=x, m=m, q=q, R=R, P0=P0, A=A, V=V, Lam=Lam)


def check_background_symbolic(mut=()):
    """the four independent background equations of the derivation (R'' from xx-tt, A'' from tt-theta-theta, V from tt, scalar
    equation), on the closed forms of S and of B, residuals 0; for B also [BBS12] (8)-(12), (15) as printed and V of (20)"""
    out = {}
    for fam in ("sv", "bbs"):
        F = sym_funcs(fam, mut)
        sp, x, R, P0, A, V, Lam = (F[k] for k in ("sp", "x", "R", "P0", "A", "V", "Lam"))
        eps = -1
        Rp, Rpp = sp.diff(R, x), sp.diff(R, x, 2)
        Ap, App = sp.diff(A, x), sp.diff(A, x, 2)
        Pp = sp.diff(P0, x)
        LamR = sp.diff(Lam, x) / Rp
        Vp = sp.diff(V, x) / Pp
        eqs = {
            "R'' + eps phi'^2 R/2": Rpp + eps * Pp**2 * R / 2,
            "A'' - [(A (R^2)'' - 2)/R^2 - R Lam_R]": App - ((A * (2 * Rp**2 + 2 * R * Rpp) - 2) / R**2 - R * LamR),
            "V + Lam + (R/2) Lam_R + (A' R^2)'/(2 R^2)": V + Lam + R * LamR / 2 + sp.diff(Ap * R**2, x) / (2 * R**2),
            "V_phi - eps (A phi'' + 2 A phi' R'/R + A' phi')": Vp - eps * (A * sp.diff(P0, x, 2) + 2 * A * Pp * Rp / R + Ap * Pp),
        }
        for k, e in eqs.items():
            out[f"{fam.upper()} {k}"] = sp.simplify(e)
        if fam == "bbs":
            # [BBS12] as printed: phi = arctan x (their normalization), eps = -1, r = R
            phi = sp.atan(x); ph1 = sp.diff(phi, x)
            q, m = F["q"], F["m"]
            r = R; r1 = Rp
            Vx = V
            printed = {
                "(8) 2 (A r^2 phi')' - eps r^2 dV/dphi": 2 * sp.diff(A * r**2 * ph1, x) - eps * r**2 * sp.diff(Vx, x) / ph1,
                "(9) (A' r^2)' + 2 r^2 V - 2 q^2/r^2": sp.diff(Ap * r**2, x) + 2 * r**2 * Vx - 2 * q**2 / r**2,
                "(10) r''/r + eps phi'^2": Rpp / r + eps * ph1**2,
                "(11) A (r^2)'' - r^2 A'' - 2 + 4 q^2/r^2": A * sp.diff(r**2, x, 2) - r**2 * App - 2 + 4 * q**2 / r**2,
                "(12) -1 + A' r r' + A r'^2 - r^2 (eps A phi'^2 - V) + q^2/r^2": -1 + Ap * r * r1 + A * r1**2 - r**2 * (eps * A * ph1**2 - Vx) + q**2 / r**2,
            }
            for k, e in printed.items():
                out["BBS printed " + k] = sp.simplify(e)
            p = 3 * m - sp.pi * q**2
            out["BBS printed (15) B' - 2 (p - x + 2 q^2 atan x)/(1+x^2)^2"] = sp.simplify(sp.diff(A / (1 + x**2), x) - 2 / (1 + x**2)**2 * (p - x + 2 * q**2 * sp.atan(x)))
    return out


def check_bbs_printed_V20(mut=()):
    """[BBS12] (20) as printed (transcribed from the PDF) against V of (9) on the closed-form A, mpmath 30 digits, members (m, q) and points"""
    import mpmath as mp
    mp.mp.dps = 30
    worst = mp.mpf(0)
    for m, q in (("0.7", "0.4"), ("0.05", "0.1"), ("0.2122", "0.0")):
        m = mp.mpf(m); q = mp.mpf(q)
        A = mp_A("bbs", m, q)
        V9 = lambda x: -mp.diff(lambda t: mp.diff(A, t) * (1 + t * t), x) / (2 * (1 + x * x)) + q * q / (1 + x * x)**2
        den = (lambda x: 2 * (1 + x * x)**2) if "bbs_V20" in mut else (lambda x: 2 * (1 + x * x))     # control: wrong reading of the denominator

        def V20(x):
            at = mp.atan(x)
            br = (2 * q * q * (3 * x * x + 1) * at**2
                  + (18 * x * x * m - 6 * mp.pi * x * x * q * q + 12 * x * q * q - 2 * mp.pi * q * q + 6 * m) * at
                  + q * q * (mp.mpf(3) / 2 * mp.pi**2 * x * x - 6 * mp.pi * x + mp.pi**2 / 2 + 6) - m * (9 * mp.pi * x * x - 18 * x + 3 * mp.pi))
            return q * q / (1 + x * x)**2 - br / den(x)
        for x in ("-2", "-0.3", "0.3", "1.7"):
            x = mp.mpf(x)
            worst = max(worst, abs(V20(x) - V9(x)))
    return worst


def check_bw_printed_vs_einstein(mut=()):
    """[BW22] (11)-(13), (16), (17), (21), (22) at Q = 0, a = 1, in their conventions (signature +---, G^nu_mu = -T^nu_mu, which is minus the
    standard Ricci sign): the printed V(phi) and L(F) with the closed-form metric, all four diagonal components, residual 0 (sympy).
    Independent of the derivation of the linearized equations (own Einstein tensor routine, mixed components)."""
    import sympy as sp
    t, x, th, ph = sp.symbols("t x theta varphi", real=True)
    m, q = sp.symbols("m q", positive=True)
    X = [t, x, th, ph]
    R = sp.sqrt(x**2 + 1); A = 1 - 2 * m / R
    g = sp.diag(A, -1 / A, -R**2, -R**2 * sp.sin(th)**2)
    gi = g.inv()
    Gm = [[[sum(gi[i, d] * (sp.diff(g[d, j], X[k]) + sp.diff(g[d, k], X[j]) - sp.diff(g[j, k], X[d])) for d in range(4)) / 2
            for k in range(4)] for j in range(4)] for i in range(4)]

    def Ric(j, k):
        s = 0
        for i in range(4):
            s += sp.diff(Gm[i][j][k], X[i]) - sp.diff(Gm[i][j][i], X[k])
            for d in range(4):
                s += Gm[i][i][d] * Gm[d][j][k] - Gm[i][k][d] * Gm[d][j][i]
        return s
    Rl = sp.Matrix(4, 4, lambda j, k: Ric(j, k))
    Rs = sum(gi[j, k] * Rl[j, k] for j in range(4) for k in range(4))
    Gstd = gi * Rl - sp.eye(4) * Rs / 2
    eps = -1
    phi = sp.atan(x)                                                           # [BW22] (19), a = 1
    cV = sp.Rational(3, 5) if "sv_V" in mut else sp.Rational(4, 5)
    cL = sp.Rational(8, 5) if "sv_L" in mut else sp.Rational(12, 5)
    V = cV * m * sp.cos(phi)**5                                                # (22) at Q = 0, a = 1
    Fs = sp.Symbol("Fs", positive=True)
    Lf = cL * m / (2 * q**2 / Fs)**sp.Rational(5, 4)                           # (21) at Q = 0, a = 1
    LF = sp.diff(Lf, Fs)
    Fv = 2 * q**2 / R**4
    dphi = [sp.diff(phi, xx) for xx in X]
    dphi2 = sum(gi[i, i] * dphi[i]**2 for i in range(4))
    Tphi = lambda mu: 2 * eps * dphi[mu] * gi[mu, mu] * dphi[mu] - (eps * dphi2 - V)             # (12)
    Fl = sp.zeros(4, 4); Fl[2, 3] = q * sp.sin(th); Fl[3, 2] = -q * sp.sin(th)
    Fu = gi * Fl * gi
    TF = lambda mu: -2 * LF.subs(Fs, Fv) * sum(Fl[mu, s] * Fu[mu, s] for s in range(4)) + Lf.subs(Fs, Fv) / 2             # (13)
    return [sp.simplify(-Gstd[mu, mu] + Tphi(mu) + TF(mu)) for mu in range(4)]


def derive_master(mut=()):
    """Linearize R_mn - eps d_m phi d_n phi - V g_mn - g_mn S_m = 0 (S_t = S_x = Lam + r Lam'/2, S_th = S_ph = Lam) and the scalar equation
    about the static background (generic A, R, phi0, Lam(R)); perturbation ~ exp(W t). The background relations used for the elimination
    are written out (not solved from the same equations). Returns residuals and the derived equations."""
    import sympy as sp
    t, rho, th, ph, eta, W = sp.symbols("t rho theta varphi eta W")
    eps = sp.Symbol("epsilon")
    A = sp.Function("A")(rho); R = sp.Function("R")(rho); P0 = sp.Function("P0")(rho)
    Vs = sp.Function("Vs")(rho); Vp = sp.Function("Vp")(rho); Vpp = sp.Function("Vpp")(rho)
    Lam = sp.Function("Lam")(rho); LamR = sp.Function("LamR")(rho); LamRR = sp.Function("LamRR")(rho)
    LamRRR = sp.Symbol("LamRRR")
    n1, l1, y1, f1 = [sp.Function(s)(t, rho) for s in ("n1", "l1", "y1", "f1")]
    X = [t, rho, th, ph]
    nu = sp.log(A) / 2 + eta * n1; lam = -sp.log(A) / 2 + eta * l1; r = R + eta * y1; phi = P0 + eta * f1
    g = sp.diag(-sp.exp(2 * nu), sp.exp(2 * lam), r**2, r**2 * sp.sin(th)**2)
    gi = g.inv()
    G = [[[sum(gi[a, d] * (sp.diff(g[d, b_], X[c_]) + sp.diff(g[d, c_], X[b_]) - sp.diff(g[b_, c_], X[d])) for d in range(4)) / 2
           for c_ in range(4)] for b_ in range(4)] for a in range(4)]

    def Ric(b_, c_):
        s = 0
        for a in range(4):
            s += sp.diff(G[a][b_][c_], X[a]) - sp.diff(G[a][b_][a], X[c_])
            for d in range(4):
                s += G[a][a][d] * G[d][b_][c_] - G[a][c_][d] * G[d][b_][a]
        return s
    Vexp = Vs + eta * Vp * f1
    Lr = Lam + eta * LamR * y1                                              # Lambda(r)
    rLr = (R + eta * y1) * (LamR + eta * LamRR * y1)                        # r Lambda'(r)
    sgnT = -1 if "ned_T" in mut else 1                                      # control: wrong sign of the transverse NED source
    S0 = Lr + rLr / 2
    S2 = sgnT * Lr
    dphi = [sp.diff(phi, x) for x in X]
    Sd = {0: S0, 1: S0, 2: S2, 3: S2}
    E = {(a, b_): Ric(a, b_) - eps * dphi[a] * dphi[b_] - Vexp * g[a, b_] - (Sd[a] * g[a, b_] if a == b_ else 0)
         for (a, b_) in [(0, 0), (0, 1), (1, 1), (2, 2)]}
    sq = sp.sqrt(-g.det())
    E["S"] = eps * sum(sp.diff(sq * gi[a, a] * dphi[a], X[a]) for a in range(4)) / sq - (Vp + eta * Vpp * f1)
    Nn, Ll, Yy, Ff = [sp.Function(s)(rho) for s in "NLYF"]
    Rp, Ap, Pp = sp.diff(R, rho), sp.diff(A, rho), sp.diff(P0, rho)
    # background relations, written from [BBS12] (9)-(11) generalized to Lambda(R) (for Maxwell Lambda = q^2/R^4: (9), (11) exactly)
    Rpp_bg = -eps * Pp**2 * R / 2                                                                     # (10)
    App_bg = (A * (2 * Rp**2 + 2 * R * Rpp_bg) - 2) / R**2 - R * LamR                                  # (11): A (R^2)'' - R^2 A'' = 2 + R^3 Lam_R
    Vs_bg = eps * A * Pp**2 / 2 - A * Rp**2 / R**2 - Lam - Ap * Rp / R + 1 / R**2                       # (9) + (12)
    Vp_bg = eps * (A * sp.diff(P0, rho, 2) + 2 * A * Pp * Rp / R + Ap * Pp)                             # scalar equation (8)

    def bg(e):
        e = e.subs(Vpp, sp.diff(Vp, rho) / Pp) if e.has(Vpp) else e
        e = e.subs(Vp, Vp_bg).doit()
        e = e.subs(Vs, Vs_bg)
        for _ in range(5):
            e = e.subs(sp.diff(Lam, rho), LamR * Rp).subs(sp.diff(LamR, rho), LamRR * Rp).subs(sp.diff(LamRR, rho), LamRRR * Rp)
            for k in (4, 3, 2):
                e = e.subs(sp.diff(R, rho, k), sp.diff(Rpp_bg, rho, k - 2))
            for k in (4, 3, 2):
                e = e.subs(sp.diff(A, rho, k), sp.diff(App_bg, rho, k - 2))
            e = e.doit()
        return sp.simplify(e)
    res = {}
    for k in [(0, 0), (0, 1), (1, 1), (2, 2), "S"]:
        res["background " + str(k)] = bg(E[k].subs(eta, 0).doit())

    def lin(k, **kw):
        e1 = sp.diff(E[k], eta).subs(eta, 0).doit()
        e1 = e1.subs({n1: kw.get("N", Nn) * sp.exp(W * t), l1: kw.get("L", Ll) * sp.exp(W * t), y1: kw.get("Y", Yy) * sp.exp(W * t),
                      f1: kw.get("F", Ff) * sp.exp(W * t)}).doit()
        return bg(sp.simplify(e1 * sp.exp(-W * t)))
    Np = sp.Symbol("Np")
    zero = sp.Integer(0)
    # --- gauge dphi = 0
    E1 = {k: lin(k, F=zero) for k in [(0, 0), (0, 1), (1, 1), (2, 2), "S"]}
    Ls = (sp.diff(Yy, rho) - Ap * Yy / (2 * A)) / Rp                              # from the t-rho equation
    res["G1 t-rho equation solved for delta lambda"] = sp.simplify(sp.solve(E1[(0, 1)], Ll)[0] - Ls)
    Npsol = sp.solve(E1["S"].subs(sp.diff(Nn, rho), Np), Np)[0]
    e22 = bg(E1[(2, 2)].subs(sp.diff(Nn, rho), Np).subs(Np, Npsol).subs(Ll, Ls).doit())
    Ypp = sp.simplify(sp.solve(e22, sp.diff(Yy, rho, 2))[0])
    if "master_p" in mut:
        Ypp = Ypp + Rp / R * sp.diff(Yy, rho)                                      # control: wrong coefficient of Y'
    if "master_ned" in mut:
        Ypp = Ypp + LamR * R * Yy / A                                              # control: NED term doubled
    for key in [(0, 0), (1, 1)]:
        e = E1[key].subs(sp.diff(Nn, rho, 2), sp.diff(Npsol, rho)).subs(sp.diff(Nn, rho), Npsol).subs(Ll, Ls).doit()
        for _ in range(3):
            e = e.subs(sp.diff(Yy, rho, 3), sp.diff(Ypp, rho)).subs(sp.diff(Yy, rho, 2), Ypp).subs(Ll, Ls).doit()
        res["G1 %s residual after the master equation" % (key,)] = bg(e)
    # --- gauge dR = 0
    E2 = {k: lin(k, Y=zero) for k in [(0, 0), (0, 1), (1, 1), (2, 2), "S"]}
    Lb = eps * Pp * R * Ff / (2 * Rp)
    res["BK t-rho equation solved for delta lambda"] = sp.simplify(sp.solve(E2[(0, 1)], Ll)[0] - Lb)
    e22 = E2[(2, 2)].subs(sp.diff(Nn, rho), Np).subs(Ll, Lb).doit()
    Npb = sp.solve(e22, Np)[0]
    S2_ = bg(E2["S"].subs(sp.diff(Nn, rho), Np).subs(Np, Npb).subs(Ll, Lb).doit())
    Fpp = sp.simplify(sp.solve(S2_, sp.diff(Ff, rho, 2))[0])
    cF1 = sp.simplify(Fpp.coeff(sp.diff(Ff, rho)))
    cF0 = sp.simplify(sp.expand(Fpp - cF1 * sp.diff(Ff, rho)).coeff(Ff))
    Ut = cF0 - W**2 / A**2
    for key in [(0, 0), (1, 1)]:
        e = E2[key].subs(sp.diff(Nn, rho, 2), sp.diff(Npb, rho)).subs(sp.diff(Nn, rho), Npb).subs(Ll, Lb).doit()
        for _ in range(3):
            e = e.subs(sp.diff(Ff, rho, 3), sp.diff(Fpp, rho)).subs(sp.diff(Ff, rho, 2), Fpp).subs(Ll, Lb).doit()
        res["BK %s residual after the master equation" % (key,)] = bg(e)
    # --- [BKZ12] (41) with V -> V + Lambda in the first term (the NED extension, found here): A Ut = eps phi'^2/R'^2 (R^2 (V + Lam) - 1) + ...
    sgn = -1 if "bkz41_sign" in mut else 1
    rhs41 = eps * Pp**2 / Rp**2 * (R**2 * (Vs + Lam) - 1) + sgn * 2 * Pp * R * Vp / Rp + eps * Vpp
    res["(41) with V -> V + Lambda: A*Ut - rhs, eps = -1"] = bg((A * Ut - rhs41).subs(eps**2, 1)).subs(eps, -1)
    res["  ... eps = +1"] = bg((A * Ut - rhs41).subs(eps**2, 1)).subs(eps, 1)
    rhs41_0 = eps * Pp**2 / Rp**2 * (R**2 * Vs - 1) + 2 * Pp * R * Vp / Rp + eps * Vpp
    ned_term = sp.simplify(bg((A * Ut - rhs41_0).subs(eps**2, 1)).subs(eps, -1))
    # --- gauge map
    Ym = sp.Function("Y")(rho)
    Fm = -Pp * Ym / Rp
    ex = (sp.diff(Fm, rho, 2) - Fpp.subs(Ff, Fm).subs(sp.diff(Ff, rho), sp.diff(Fm, rho))).doit()
    Yp2 = Ypp.subs(Yy, Ym)
    for _ in range(3):
        ex = ex.subs(sp.diff(Ym, rho, 3), sp.diff(Yp2, rho)).subs(sp.diff(Ym, rho, 2), Yp2).doit()
    res["gauge map F = -phi0' Y/R' residual"] = bg(ex)
    return res, dict(Ypp=Ypp, Fpp=Fpp, Ut=Ut, cF1=cF1, A=A, R=R, P0=P0, rho=rho, W=W, eps=eps, Yy=Yy, Vs=Vs, Vp=Vp, Vpp=Vpp, Lam=Lam,
                     LamR=LamR, ned_term=ned_term)


def check_ned_source_and_bianchi(mut=()):
    """(i) T_mn[F] = 2 L_F F_ms F_n^s - g_mn L/2 for F_thph = q sin(th) on the perturbed metric, generic L(F): T^t_t = T^x_x = -L/2,
    T^th_th = F L_F - L/2, its trace-reversed form equals g_mn S_m of the derivation, and F = 2 q^2/r^4 exactly: the NED source depends on
    the metric only through r. (ii) dF = 0 with F_thph = f(t,x) sin th forces f' = 0 = f_t. (iii) the NED equation
    d_m(sqrt(-g) L_F F^{mn}) = 0 holds for the perturbed metric with F_thph = q sin th. All residuals 0 (sympy)."""
    import sympy as sp
    t, x, th, ph, eta = sp.symbols("t x theta varphi eta")
    q = sp.Symbol("q", positive=True)
    X = [t, x, th, ph]
    nu = sp.Function("nu")(t, x); lam = sp.Function("lam")(t, x); r = sp.Function("r")(t, x)
    g = sp.diag(-sp.exp(2 * nu), sp.exp(2 * lam), r**2, r**2 * sp.sin(th)**2)
    gi = g.inv()
    Fl = sp.zeros(4, 4); Fl[2, 3] = q * sp.sin(th); Fl[3, 2] = -q * sp.sin(th)
    Fu = gi * Fl * gi
    Fin = sp.simplify(sum(Fl[a, b] * Fu[a, b] for a in range(4) for b in range(4)))
    out = {}
    out["F = F_mn F^mn - 2 q^2/r^4"] = sp.simplify(Fin - 2 * q**2 / r**4)
    Fs = sp.Symbol("Fs", positive=True)
    Lf = sp.Function("Lf")
    L_, LF_ = Lf(Fin), sp.diff(Lf(Fs), Fs).subs(Fs, Fin)
    Tmix = sp.Matrix(4, 4, lambda mu, nu_: 2 * LF_ * sum(Fl[mu, s] * Fu[nu_, s] for s in range(4)) - sp.KroneckerDelta(mu, nu_) * L_ / 2)
    # (-+++): T^t_t = -rho; Fl[mu, s] F^{nu s}: lower first index, upper second
    out["T^t_t + L/2"] = sp.simplify(Tmix[0, 0] + L_ / 2)
    out["T^x_x + L/2"] = sp.simplify(Tmix[1, 1] + L_ / 2)
    out["T^th_th - (F L_F - L/2)"] = sp.simplify(Tmix[2, 2] - (Fin * LF_ - L_ / 2))
    # trace-reversed: R_mn^(F) = T_mn - g_mn T/2 = g_mn S_m, S_t = S_x = L/2 - F L_F, S_th = L/2; derivation: S_t = Lam + r Lam'/2, S_th = Lam with Lam = L/2
    Ttr = sum(Tmix[i, i] for i in range(4))
    S_t = sp.simplify(Tmix[0, 0] - Ttr / 2); S_th = sp.simplify(Tmix[2, 2] - Ttr / 2)
    rr = sp.Symbol("rr", positive=True)
    Lam_r = Lf(2 * q**2 / rr**4) / 2
    sgnT = -1 if "ned_T" in mut else 1
    out["S_t - (Lam + r Lam'/2)"] = sp.simplify((S_t - (Lam_r + rr * sp.diff(Lam_r, rr) / 2).subs(rr, r)).subs(Fin, 2 * q**2 / r**4).doit())
    out["S_th - Lam"] = sp.simplify((S_th - sgnT * Lam_r.subs(rr, r)).subs(Fin, 2 * q**2 / r**4).doit())
    # (ii) Bianchi
    f = sp.Function("f")(t, x)
    Fg = sp.zeros(4, 4); Fg[2, 3] = f * sp.sin(th); Fg[3, 2] = -f * sp.sin(th)
    bianchi = {}
    for (a, b, c) in [(0, 2, 3), (1, 2, 3), (0, 1, 2), (0, 1, 3), (1, 2, 3)]:
        bianchi[(a, b, c)] = sp.simplify(sp.diff(Fg[b, c], X[a]) + sp.diff(Fg[c, a], X[b]) + sp.diff(Fg[a, b], X[c]))
    out["dF = 0: (t,th,ph) component is f_t sin th"] = sp.simplify(bianchi[(0, 2, 3)] - sp.diff(f, t) * sp.sin(th))
    out["dF = 0: (x,th,ph) component is f_x sin th"] = sp.simplify(bianchi[(1, 2, 3)] - sp.diff(f, x) * sp.sin(th))
    # (iii) NED field equation
    sq = sp.sqrt(-g.det())
    for nu_ in range(4):
        e = sum(sp.diff(sq * LF_ * Fu[mu, nu_], X[mu]) for mu in range(4))
        out["NED equation, nu = %d" % nu_] = sp.simplify(e)
    return out


def family_coeffs_generic(D, fam, mval, qval, mut=()):
    """returns (generic - family) of the Y'' equation, mpmath 30 digits, at several x"""
    import sympy as sp
    import mpmath as mp
    mp.mp.dps = 30
    rho = D["rho"]; A = D["A"]; R = D["R"]; P0 = D["P0"]; W = D["W"]; Yy = D["Yy"]
    Rf = sp.sqrt(rho**2 + 1); Pf = sp.sqrt(2) * sp.atan(rho)
    if fam == "sv":
        mm = sp.Rational(mval)
        Af = 1 - 2 * mm / Rf
        LamRf = -6 * mm / Rf**6                                           # Lambda = 6m/(5 R^5): dLambda/dR
    else:
        mm = sp.Rational(mval[0], mval[1]) if isinstance(mval, tuple) else sp.nsimplify(mval)
        qq = sp.Rational(qval)
        p = 3 * mm - sp.pi * qq**2
        B0 = -sp.pi * p / 2 - sp.pi**2 * qq**2 / 4
        Af = (1 + rho**2) * (B0 + (1 + qq**2 + p * rho) / (1 + rho**2) + (p + 2 * qq**2 * rho / (1 + rho**2)) * sp.atan(rho) + qq**2 * sp.atan(rho)**2)
        LamRf = -4 * qq**2 / Rf**5
    Ypp = D["Ypp"].subs(D["eps"], -1)
    Y0, Y1 = sp.symbols("Y0 Y1")
    ex = Ypp.subs(sp.diff(Yy, rho), Y1).subs(Yy, Y0)
    ex = ex.subs({D["LamR"]: LamRf}).subs({A: Af, R: Rf, P0: Pf}).doit()
    ned = 0 if "family_ned" in mut else 1
    # family: Y'' = (2/x - A'/A) Y' + (1/R^4 - A'/(A x R^2) - 2/(A R^2) - R Lam_R/A + W^2/A^2) Y
    fam_e = (2 / rho - sp.diff(Af, rho) / Af) * Y1 + (1 / Rf**4 - sp.diff(Af, rho) / (Af * rho * Rf**2) - 2 / (Af * Rf**2)
                                                        - ned * Rf * LamRf / Af + W**2 / Af**2) * Y0
    d = sp.lambdify((rho, Y0, Y1, W), ex - fam_e, "mpmath")
    return max(abs(d(mp.mpf(r), mp.mpf("0.37"), mp.mpf("-0.81"), mp.mpf("0.7"))) for r in ("-2.3", "-0.6", "0.35", "1.4", "5.0"))


def check_vhat_vs_bk(D, fam, mval, qval, mut=()):
    """V-hat of the u-form against A*Ut of the delta R = 0 gauge plus R_{x*x*}/R (mpmath 30 digits). Ut from the derived generic
    expression, V, V_phi, V_phiphi, Lambda of the family."""
    import sympy as sp
    import mpmath as mp
    mp.mp.dps = 30
    F = sym_funcs(fam)
    x, R, P0, A, V, Lam = F["x"], F["R"], F["P0"], F["A"], F["V"], F["Lam"]
    sub = {F["m"]: sp.nsimplify(mval) if fam == "bbs" else sp.Rational(mval)}
    if fam == "bbs":
        sub[F["q"]] = sp.Rational(qval)
    Pp = sp.diff(P0, x)
    Vp = sp.diff(V, x) / Pp
    Vpp = sp.diff(Vp, x) / Pp
    Rp = sp.diff(R, x)
    ned = 0 if "veff_ned" in mut else 1
    AUt = -Pp**2 / Rp**2 * (R**2 * (V + ned * Lam) - 1) + 2 * Pp * R * Vp / Rp - Vpp                    # eps = -1, derived (41) with V -> V + Lambda
    Rxx = (A * sp.diff(A, x) * Rp + A**2 * sp.diff(Rp, x)) / R
    Veff = A * AUt + Rxx                                                                           # U = A^2 Ut = A * (A Ut)
    Rr = R
    Vhat = 2 * A**2 / x**2 - (1 + 1 / Rr**2) * A * sp.diff(A, x) / x + A**2 / Rr**4 - 2 * A / Rr**2 - A * Rr * sp.diff(Lam, x) / Rp
    f = sp.lambdify(x, (Veff - Vhat).subs(sub), "mpmath")
    return max(abs(f(mp.mpf(v))) for v in ("-1.3", "-0.4", "0.7", "2.4"))


def check_numerics_vs_mp(mut=()):
    """float64 implementation (A, A', series of A/z at the horizon, large-x form) against mpmath at 40 digits: |dA|/|A| <= 1e-11 where
    |A| > 1e-3, |dA'| <= 1e-11 absolute (A' enters the potential multiplied by A/x), |d(A/z)|/|A/z| <= 1e-11 with the exact horizon"""
    import mpmath as mp
    mp.mp.dps = 40
    worst = 0.0
    cases = [("sv", 0.3, 0.0, False), ("sv", 0.75, 0.0, False), ("sv", 10.0, 0.0, False), ("sv", 0.5, 0.0, False),
             ("bbs", B1_M, 0.5, False), ("bbs", m_throat(0.6), 0.6, True), ("bbs", B3_M, 0.2, False)]
    for fam, m, q, thr in cases:
        bg = Bg(fam, m, q, "R2", thr)
        A = mp_A(fam, mp.mpf(m), mp.mpf(q))
        for x in (-3.0, -0.4, 0.3, 1.0, 2.5, 40.0, 3000.0):
            a = float(A(mp.mpf(x)))
            worst = max(worst, abs(float(bg.A(x)) - a) / abs(a) if abs(a) > 1e-3 else abs(float(bg.A(x)) - a))
            ap = float(mp.diff(A, mp.mpf(x)))
            worst = max(worst, abs(float(bg.Ap(x)) - ap))
        if bg.xh is not None and bg.kind != "ext0":
            if fam == "sv":
                xh = mp.sqrt(4 * mp.mpf(m)**2 - 1)
            elif thr:                                                   # the horizon is exactly at x = 0 at the exact throat mass
                qm = mp.mpf(q)
                A = mp_A("bbs", 2 / (3 * mp.pi) * (1 + qm**2 + mp.pi**2 * qm**2 / 4), qm)
                xh = mp.mpf(0)
            else:
                xh = mp.findroot(A, mp.mpf(bg.xh))
            for z in (1e-9, 1e-5, 1e-3, 0.02, 0.049, 0.06, 0.3):
                a = float(A(xh + mp.mpf(z)) / mp.mpf(z))
                worst = max(worst, abs(bg.a_over_z(z) - a) / abs(a))
    # low stage, B: A and A' from the series in 1/x for x >= 30 (amendment A7): |dA| <= 1e-14, |dA'| <= 1e-12 |A'|
    for m, q, thr in ((B1_M, 0.5, False), (m_throat(0.7), 0.7, True), (B3_M, 0.2, False), (-2 * (-10.0) / (3 * math.pi), 0.0, False)):
        bg = Bg("bbs", m, q, "R2", thr, "low")
        A = mp_A("bbs", mp.mpf(m), mp.mpf(q))
        for x in (30.0, 100.0, 1e3, 1e4, 1e5, 2e6):
            a = float(A(mp.mpf(x)))
            ap = float(mp.diff(A, mp.mpf(x), h=mp.mpf(x) * mp.mpf("1e-20")))
            worst = max(worst, abs(float(bg.A(x)) - a) / 1e-3, abs(float(bg.Ap(x)) - ap) / abs(ap) / 1e-1)
    return worst


def check_horizons(mut=()):
    """S: A'(x_h) = 2 m x_h / R_h^3 at x_h = sqrt(4 m^2 - 1) against mpmath; B2: A(0) = 0 and A'(0) = 2p at the throat mass m_c(q)"""
    import mpmath as mp
    mp.mp.dps = 40
    worst = mp.mpf(0)
    for m in (0.55, 1.0, 3.0, 10.0):
        A = mp_A("sv", mp.mpf(m), 0)
        xh = mp.sqrt(4 * mp.mpf(m)**2 - 1)
        k = 2 * mp.mpf(m) * xh / (xh * xh + 1)**mp.mpf(1.5)
        if "bg_kappa" in mut:
            k = k / 2
        worst = max(worst, abs(A(xh)), abs(mp.diff(A, xh) - k))
    for q in (0.0, 0.3, 0.6, 0.81):
        m = 2 / (3 * mp.pi) * (1 + mp.mpf(q)**2 + mp.pi**2 * mp.mpf(q)**2 / 4)
        A = mp_A("bbs", m, mp.mpf(q))
        p = 3 * m - mp.pi * mp.mpf(q)**2
        worst = max(worst, abs(A(mp.mpf(0))), abs(mp.diff(A, 0) - 2 * p))
    return worst


def check_command():
    ok_all = True

    def line(tag, ok, msg):
        nonlocal ok_all
        ok_all &= ok
        print(("PASS  " if ok else "FAIL  ") + tag + ": " + msg, flush=True)
    t0 = time.time()
    for k, e in check_background_symbolic(MUT).items():
        line("S1 background (sympy, closed form)", e == 0, f"{k}: residual {e}")
    for mu, e in enumerate(check_bw_printed_vs_einstein(MUT)):
        line("S1 [BW22] printed source against the Einstein tensor (sympy, Q = 0)", e == 0, f"component {'txth'[mu] if mu < 2 else 'th'} ({mu}): residual {e}")
    w20 = check_bbs_printed_V20(MUT)
    line("S1 [BBS12] printed (20) for V against (9) (mpmath 30 digits, 3 members, 4 points)", w20 < 1e-25, f"max abs difference {float(w20):.3e} (declared <= 1e-25)")
    mpx = check_numerics_vs_mp(MUT)
    line("S2 float64 background (A, A', A/z series, large x; low-stage series of B) against mpmath 40 digits, 7 members", mpx < 1e-11, f"max scaled difference {mpx:.3e} (declared <= 1e-11)")
    import mpmath as mp
    worst = check_horizons(MUT)
    line("S3 horizons: S A'(x_h) = 2 m x_h/R_h^3; B2 A(0) = 0, A'(0) = 2p (mpmath 40 digits)", worst < 1e-30, f"max abs difference {float(worst):.3e} (declared <= 1e-30)")
    print(f"      elapsed {time.time() - t0:.0f} s", flush=True)
    for k, e in check_ned_source_and_bianchi(MUT).items():
        line("S4 NED source, Bianchi and NED equation (sympy)", e == 0, f"{k}: residual {e}")
    res, D = derive_master(MUT)
    for k, e in res.items():
        line("S4 linearized equations (sympy, generic A, R, phi0, Lambda(R))", e == 0, f"{k}: {e if e != 0 else 0}")
    print(f"      elapsed {time.time() - t0:.0f} s", flush=True)
    print("      NED extension of [BKZ12] (41): A*Ut - (41 as printed) = " + str(D["ned_term"]) + "   (eps = -1)", flush=True)
    for fam, mv, qv in (("sv", "3/10", None), ("sv", "7/10", None), ("bbs", "7/10", "2/5"), ("bbs", "1/20", "1/10")):
        w = family_coeffs_generic(D, fam, mv, qv, MUT) if fam == "sv" else family_coeffs_generic(D, fam, (7, 10) if mv == "7/10" else (1, 20), qv, MUT)
        line(f"S5 family coefficients vs generic derived equation (mpmath 30 digits, {fam} m = {mv}" + (f", q = {qv})" if qv else ")"), w < 1e-25, f"max abs difference {float(w):.3e} (declared <= 1e-25)")
    for fam, mv, qv in (("sv", "3/10", None), ("sv", "13/10", None), ("bbs", "7/10", "2/5"), ("bbs", "1/20", "1/10")):
        w = check_vhat_vs_bk(D, fam, mv if fam == "sv" else (0.7 if mv == "7/10" else 0.05), qv, MUT)
        line(f"S8 Vhat of u = Y/x against A Ut + R_x*x*/R of the delta R = 0 gauge (mpmath 30 digits, {fam} m = {mv}" + (f", q = {qv})" if qv else ")"), w < 1e-25, f"max abs difference {float(w):.3e} (declared <= 1e-25)")
    print(f"      elapsed {time.time() - t0:.0f} s", flush=True)
    return ok_all


# ------------------------------------------------------------------------------------------------ independent spectral check
def xstar_map(bg, xs, X):
    """radius x at each tortoise value xs (array in [-X, X]); x_* = 0 at x = 0 (throat cases) or at x_h + 1 (no throat)"""
    xs = np.asarray(xs)
    if bg.kind == "twoend":
        out = np.empty_like(xs)
        for sgn_, mask in ((1, xs >= 0), (-1, xs < 0)):
            pts = xs[mask]
            if len(pts) == 0:
                continue
            order = np.argsort(np.abs(pts))
            sol = solve_ivp(lambda s, y: [bg.A(y[0])], [0, sgn_ * X * 1.0001], [0.0], method="DOP853", rtol=1e-13, atol=1e-14,
                            t_eval=pts[order])
            tmp = np.empty(len(pts)); tmp[order] = sol.y[0]
            out[mask] = tmp
        return out
    xref = 0.0 if bg.xh < 0 else bg.xh + 1.0
    zref = xref - bg.xh
    z = np.empty_like(xs)
    neg = xs < 0
    if neg.any():
        pts = xs[neg]
        sol = solve_ivp(lambda s, q: [bg.a_over_z(np.exp(q[0]))], [0, -X * 1.0001], [math.log(zref)], method="DOP853", rtol=1e-13, atol=1e-13,
                        t_eval=np.sort(pts)[::-1])
        z[neg] = np.exp(sol.y[0])[np.argsort(np.argsort(-pts))]
    pos = ~neg
    sol = solve_ivp(lambda s, q: [bg.a_over_z(np.exp(q[0]))], [0, X * 1.0001], [math.log(zref)], method="DOP853", rtol=1e-13, atol=1e-13,
                    t_eval=np.sort(xs[pos]))
    z[pos] = np.exp(sol.y[0])[np.argsort(np.argsort(xs[pos]))]
    return bg.xh + z


def spectral_member(spec, which, W_ref):
    """eigenvalues lambda = W^2 of  u_xx - Vhat u = lambda u  (no throat in the static region: x_h >= 0) or of the Y-form
    Y_xx = Cx Y_x + Dx Y + lambda Y (throat in the region), Chebyshev collocation on x_* in [-X, X], Dirichlet ends. Shares with the shooting
    code the background and the coefficients; nothing of the boundary conditions, the throat crossing, the matching or the root finding."""
    bg = Bg(spec["fam"], spec["m"], spec["q"], "R2", spec["at_throat"])
    cfg = SPEC[which]
    N = cfg["N"]
    X = max(cfg["X0"], cfg["kX"] / W_ref) if W_ref else cfg["X0"]
    D, xc = S.cheb(N)
    if X > MAP_X:                                                       # amendment A4: mapped grid x_* = X sinh(a xi)/sinh(a) where the plain grid cannot
        xs = X * np.sinh(MAP_A * xc) / np.sinh(MAP_A)                   # resolve the structure at the centre (long tails, 1/W up to 3600)
        jac = X * MAP_A * np.cosh(MAP_A * xc) / np.sinh(MAP_A)
    else:
        xs = X * xc
        jac = X * np.ones_like(xc)
    rho = xstar_map(bg, xs, X)
    A = bg.A(rho); Ap = bg.Ap(rho)
    if bg.xh is not None and bg.xh == 0.0:
        Ar = np.array([bg.a_over_z(r) for r in rho])
    else:
        Ar = A / rho
    D1 = D / jac[:, None]
    D2 = D1 @ D1
    sgn = -1.0 if mutating("spec_sign") else 1.0
    if not has_throat(bg):
        R2 = rho * rho + 1
        Vh = 2 * Ar**2 - (1 + 1 / R2) * Ar * Ap + A * A / R2**2 - 2 * A / R2 - bg.ned * A * bg.RL(rho)
        Lop = D2 - sgn * np.diag(Vh)
        form = "u"
    else:
        Cx, Dx = coef_x(bg, rho, A, Ar, Ap)
        Lop = D2 - Cx[:, None] * D1 - sgn * np.diag(Dx)
        form = "Y"
    ev = np.linalg.eigvals(Lop[1:-1, 1:-1])
    good = ev[(ev.real > LAM_MIN) & (ev.real < 100)]
    return form, np.sort_complex(good), X, bg


def load_results(d):
    rows = {}
    for fn in sorted(glob.glob(d + "/*.json")):
        o = json.load(open(fn))
        key = o["res"] if o.get("stage", "main") == "main" else "low_" + o["res"]
        rows.setdefault(o["id"], {})[key] = o
    return rows


def spectral_command(d):
    ok_all = True
    rows = load_results(d)
    for mid in sorted(rows, key=lambda s: member_ids().index(s) if s in member_ids() else 0):
        o = rows[mid].get("R2") or rows[mid].get("R1")
        spec = member_spec(mid)
        ref = [r["W"] for r in o["roots"]]
        for which in ("S-A", "S-B"):
            form, ev, X, bg = spectral_member(spec, which, min(ref) if ref else None)               # box from the smallest rate
            realpos = [e.real for e in ev if abs(e.imag) < 1e-6 * abs(e)]
            cplx = [e for e in ev if abs(e.imag) >= 1e-6 * abs(e) and e.real > 1e-3]
            tol = SPEC[which]["tol_throat" if has_throat(bg) else "tol_nothroat"]
            ok = len(realpos) == len(ref) and not cplx and all(abs(math.sqrt(a) - b) / b <= tol for a, b in zip(sorted(realpos), sorted(ref)))
            ok_all &= ok
            rel = [abs(math.sqrt(a) - b) / b for a, b in zip(sorted(realpos), sorted(ref))]
            print(("PASS  " if ok else "FAIL  ") + f"S6 spectral {which} ({form}-form, X = {X:.0f}) {mid}: positive real eigenvalues {len(realpos)} "
                  f"(shooting {len(ref)}), complex with Re > 1e-3: {len(cplx)}, W = {[round(math.sqrt(a), 9) for a in sorted(realpos)]}, "
                  f"|dW|/W = {['%.2e' % r for r in rel]} (declared <= {tol:g})", flush=True)
    return ok_all


# ------------------------------------------------------------------------------------------------ report
def table_212():
    """c -> (n_modes, W) from stability_output.txt (the report table of #212)"""
    tab = {}
    for line in open(TABLE_FILE):
        m = re.match(r"\s+(-?\d+(?:\.\d+)?)\s+(\S+)\s+(\S+)\s*\|\s*(\d+)\s+(\d+)\s*\|\s*(none|\d\.\d+)", line)
        if m:
            tab[float(m.group(1))] = (int(m.group(5)), None if m.group(6) == "none" else float(m.group(6)))
    return tab


def agree_resolutions(o1, o2):
    same = o1["n_roots"] == o2["n_roots"]
    rel = max([abs(a["W"] - b["W"]) / b["W"] for a, b in zip(o1["roots"], o2["roots"])] or [0.0])
    return same and rel <= ROOT_AGREE, rel


def report_command(d):
    rows = load_results(d)
    ok_all = True
    tab = table_212()
    print("Radial modes, b = 1 (W in units 1/b). Main stage: W in [%g, %g], R1 and R2; low stage: W in [%g, %g], R1 and R2; Sturm: nodes of the zero-energy "
          "solution (x_h >= 0, no range)" % (W_LO, W_HI, W_LOW_LO, W_LOW_HI))
    print("%-10s %7s %8s %9s %9s %7s | %-3s %-3s | %-34s | %-3s %-3s | %-6s | %-5s | %s" % (
        "id", "m", "q", "x_h", "kappa", "p", "nR1", "nR2", "W main (R2)", "lR1", "lR2", "Sturm", "count", "verdicts"))
    total = {}
    for mid in member_ids():
        if mid not in rows:
            continue
        r1, r2, l1, l2 = (rows[mid].get(k) for k in ("R1", "R2", "low_R1", "low_R2"))
        ref = r2 or r1
        notes = []
        if r1 and r2:
            ok, rel = agree_resolutions(r1, r2)
            ok_all &= ok
            notes.append(f"R1 v R2 {rel:.1e} " + ("ok" if ok else "FAIL"))
        n_main = ref["n_roots"]
        n_low = None
        if l1 and l2:
            same = l1["n_roots"] == l2["n_roots"]
            rel = max([abs(a["W"] - b["W"]) / b["W"] for a, b in zip(l1["roots"], l2["roots"])] or [0.0])
            good = same and rel <= ROOT_AGREE_LOW
            ok_all &= good
            n_low = l2["n_roots"]
            notes.append(f"low R1 v R2 {rel:.1e} " + ("ok" if good else "FAIL"))
        st = ref["sturm"]
        count = None if n_low is None else n_main + n_low
        if st is not None:
            s9 = (st == n_main)
            ok_all &= s9
            notes.append("S9 " + ("ok" if s9 else f"NOT MET (Sturm {st}, shooting in range {n_main})"))
            if n_low is not None:
                if st == n_main + n_low:
                    count = st
                elif st > n_main + n_low:
                    count = st
                    notes.append(f"{st - n_main - n_low} mode(s) with W < {W_LOW_LO:g}: Sturm count is the count, rate unresolved")
                else:
                    ok_all = False
                    notes.append(f"FAIL: shooting finds {n_main + n_low} modes, Sturm {st}")
        label = mid if mid not in P2_ALIAS else mid
        c212 = P2_ALIAS.get(mid) if mid in P2_ALIAS else (float(mid.split("_c")[1]) if mid.startswith("P2_") else None)
        if c212 is not None:
            n212, W212 = tab[c212]
            for rr, o in (("R1", r1), ("R2", r2)):
                if not o:
                    continue
                good = o["n_roots"] == n212 and (W212 is None or abs(o["roots"][0]["W"] - W212) / W212 <= TABLE_TOL[rr])
                ok_all &= good
                dW = "-" if W212 is None else "%.1e" % (abs(o["roots"][0]["W"] - W212) / W212)
                notes.append(f"P2 {rr} vs #212 c = {c212:g} (n {n212}, W {W212}): dW/W {dW} " + ("ok" if good else "FAIL"))
        Ws = "; ".join(f"{x['W']:.10f}" for x in ref["roots"]) or "none"
        Wl = "; ".join(f"{x['W']:.6e}" for x in (l2 or l1 or {"roots": []})["roots"])
        if Wl:
            Ws += "  | low " + Wl
        total[mid] = count
        print("%-10s %7.4f %8.4f %9s %9s %7s | %-3s %-3s | %-34s | %-3s %-3s | %-6s | %-5s | %s" % (
            label, ref["m"], ref["q"], "-" if ref["xh"] is None else "%.5f" % ref["xh"], "-" if ref["kappa"] is None else "%.5f" % ref["kappa"],
            "-" if ref["p"] is None else "%.4f" % ref["p"], r1["n_roots"] if r1 else "-", r2["n_roots"] if r2 else "-", Ws[:34] if len(Ws) <= 34 else Ws,
            l1["n_roots"] if l1 else "-", l2["n_roots"] if l2 else "-", "-" if st is None else st, "-" if count is None else count, "; ".join(notes)))
    print("members with at least one growing mode (count, else main stage): %d of %d" % (
        sum(1 for m in rows if (total[m] if total.get(m) is not None else (rows[m].get("R2") or rows[m].get("R1"))["n_roots"]) > 0), len(rows)))
    return ok_all


# ------------------------------------------------------------------------------------------------ controls
CONTROLS = {
    # mutation: (kind, target); check/planted: the command must exit 1 with a FAIL line containing the target;
    # agree: R1 with the mutation against an unmutated R2 must fail agree_resolutions; sturm: the shooting count and the Sturm count
    # (one of them mutated) must differ; spectral: `spectral` on one member must give a FAIL line containing the target
    "sv_V": ("check", "S1"), "sv_L": ("check", "S1"), "bbs_A": ("check", "S1 background"), "bbs_V20": ("check", "S1 [BBS12] printed (20)"),
    "ned_T": ("check", "S4 NED source"), "master_ned": ("check", "S4 linearized"), "master_p": ("check", "S4 linearized"),
    "bkz41_sign": ("check", "S4 linearized"), "family_ned": ("check", "S5 family"), "veff_ned": ("check", "S8 Vhat"),
    "bg_kappa": ("check", "S3 horizons"), "az_series": ("check", "S2 float64"), "ser_coef": ("check", "S2 float64"),
    "eps": ("planted", "P1 planted"), "throat_wall": ("planted", "P1 planted"), "rmax": ("planted", "P1 planted"),
    "hz_start": ("agree", "R1 v R2"), "ned_off": ("agree", "R1 v R2"),
    "lowW_V": ("plantedlow", "P3 planted"), "lowW_range": ("plantedlow", "P3 planted"), "rmax_low": ("agreelow", "low R1 v R2"),
    "scan_range": ("sturm", "S9 Sturm"), "sturm_sign": ("sturm", "S9 Sturm"), "spec_sign": ("spectral", "S6 spectral"),
}
CONTROL_MEMBER = "B3_q0.15"
CONTROL_MEMBER_STURM = "B1_q0.3"
CONTROL_MEMBER_LOW = "S_m1"


def run_count(mid, res, mut=None, out=None):
    args = [sys.executable, __file__, "count", mid, "--res=" + res] + (["--mutate=" + mut] if mut else []) + (["--out=" + out] if out else [])
    r = subprocess.run(args, capture_output=True, text=True)
    return json.loads(r.stdout[r.stdout.index("{"):])


def controls_command(basedir=None, only=None):
    ok_all = True
    base = None
    lowbase = None
    for m, (kind, target) in CONTROLS.items():
        if only and m not in only:
            continue
        ret = 0
        if kind in ("check", "planted"):
            args = [sys.executable, __file__, kind] + (["--res=R1"] if kind == "planted" else []) + ["--mutate=" + m]
            r = subprocess.run(args, capture_output=True, text=True)
            txt = r.stdout + r.stderr
            fired = r.returncode == 1 and any(target in l and l.startswith("FAIL") for l in txt.splitlines())
            detail = "; ".join(l[:140] for l in txt.splitlines() if l.startswith("FAIL"))[:300]
            ret = r.returncode
        elif kind == "agree":
            o = run_count(CONTROL_MEMBER, "R1", m)
            if base is None:
                fn = None if not basedir else basedir + "/%s_R2.json" % CONTROL_MEMBER
                base = json.load(open(fn)) if fn and os.path.exists(fn) else run_count(CONTROL_MEMBER, "R2")
            good, rel = agree_resolutions(o, base)
            fired = not good
            detail = f"{CONTROL_MEMBER}: n_roots {o['n_roots']} against {base['n_roots']}; max |dW|/W = {rel:.3e} (declared <= {ROOT_AGREE:g})"
        elif kind == "plantedlow":
            r = subprocess.run([sys.executable, __file__, "planted_low", "--res=R1", "--mutate=" + m], capture_output=True, text=True)
            txt = r.stdout + r.stderr
            fired = r.returncode == 1 and any(target in l and l.startswith("FAIL") for l in txt.splitlines())
            detail = "; ".join(l[:160] for l in txt.splitlines() if l.startswith("FAIL"))[:300]
            ret = r.returncode
        elif kind == "agreelow":
            o = json.loads(subprocess.run([sys.executable, __file__, "count", CONTROL_MEMBER_LOW, "--res=R1", "--stage=low", "--mutate=" + m],
                                          capture_output=True, text=True).stdout)
            fn = None if not basedir else basedir + "/%s_R2_low.json" % CONTROL_MEMBER_LOW
            if lowbase is None:
                lowbase = json.load(open(fn)) if fn and os.path.exists(fn) else json.loads(subprocess.run(
                    [sys.executable, __file__, "count", CONTROL_MEMBER_LOW, "--res=R2", "--stage=low"], capture_output=True, text=True).stdout)
            same = o["n_roots"] == lowbase["n_roots"]
            rel = max([abs(a["W"] - b["W"]) / b["W"] for a, b in zip(o["roots"], lowbase["roots"])] or [0.0])
            fired = not (same and rel <= ROOT_AGREE_LOW)
            detail = f"{CONTROL_MEMBER_LOW} low stage: n_roots {o['n_roots']} against {lowbase['n_roots']}; max |dW|/W = {rel:.3e} (declared <= {ROOT_AGREE_LOW:g})"
        elif kind == "sturm":
            o = run_count(CONTROL_MEMBER_STURM, "R1", m)
            fired = o["n_roots"] != o["sturm"]
            detail = f"{CONTROL_MEMBER_STURM}: shooting n_roots {o['n_roots']} in [{W_LO:g}, {W_HI:g}], Sturm {o['sturm']}"
        elif kind == "spectral":
            import tempfile
            tmp = tempfile.mkdtemp(prefix="nm_i223_ctl_")
            run_count(CONTROL_MEMBER_STURM, "R1", None, out=tmp + "/%s_R1.json" % CONTROL_MEMBER_STURM)
            r = subprocess.run([sys.executable, __file__, "spectral", tmp, "--mutate=" + m], capture_output=True, text=True)
            txt = r.stdout + r.stderr
            fired = r.returncode == 1 and any(target in l and l.startswith("FAIL") for l in txt.splitlines())
            detail = "; ".join(l[:140] for l in txt.splitlines() if l.startswith("FAIL"))[:300]
            ret = r.returncode
        print(("fires   " if fired else "SILENT  ") + f"--mutate={m} (target {target}): exit {ret}; {detail}", flush=True)
        ok_all &= fired
    return ok_all


def main():
    global MUT
    argv = [a for a in sys.argv[1:] if not a.startswith("--mutate=")]
    for a in sys.argv[1:]:
        if a.startswith("--mutate="):
            MUT |= set(a.split("=", 1)[1].split(","))
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if a.startswith("--") and "=" in a}
    pos = [a for a in argv if not a.startswith("--")]
    cmd = pos[0]
    if cmd == "check":
        sys.exit(0 if check_command() else 1)
    if cmd == "planted":
        sys.exit(0 if planted_command((opts["--res"],) if "--res" in opts else ("R1", "R2")) else 1)
    if cmd == "planted_low":
        sys.exit(0 if planted_low_command((opts["--res"],) if "--res" in opts else ("R1", "R2")) else 1)
    if cmd == "count":
        o = count_member(member_spec(pos[1]), opts.get("--res", "R2"), W_LO, W_HI, None, opts.get("--stage", "main"))
        o["id"] = pos[1]
        s = json.dumps(o, indent=1)
        if "--out" in opts:
            open(opts["--out"], "w").write(s)
        print(s)
        sys.exit(0)
    if cmd == "members":
        print("\n".join(member_ids()))
        sys.exit(0)
    if cmd == "report":
        sys.exit(0 if report_command(pos[1]) else 1)
    if cmd == "spectral":
        sys.exit(0 if spectral_command(pos[1]) else 1)
    if cmd == "controls":
        sys.exit(0 if controls_command(opts.get("--base"), opts["--only"].split(",") if "--only" in opts else None) else 1)
    raise SystemExit(__doc__)


if __name__ == "__main__":
    main()

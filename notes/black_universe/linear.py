#!/usr/bin/env python3
"""Linear double-null response of the phantom black universes about the static member, across the horizon (#232).

The three evolution equations of #217 (evolve_derive.py, `R_mn - eps d_m phi d_n phi - V g_mn = 0`, eps = -1, 8 pi G = 1, D = 4)

    r_pm   = -(r_p r_m + e^{-f}/2)/r + V r e^{-f}/2
    f_pm   = -2 (r_p r_m + e^{-f}/2)/r^2 + eps phi_p phi_m
    phi_pm = -(r_p phi_m + r_m phi_p)/r - (eps/2) e^{-f} V'(phi)

are linearized about the static member (derivation and checks: linear_derive.py). The linear system U_pm = B0 U + Bp U_p + Bm U_m,
U = (dr, df, dphi), is marched on a characteristic rectangle with an implicit second-order cell scheme. It shares with evolve.py only the
static member and the potential (evolve_potential.py).

Chart. The null coordinates are p = x+ (the tortoise advanced time of the nonlinear code, t = (p + m)/sqrt2, x = (p - m)/sqrt2, rho(x = 0) = rho_c)
and sigma, a null coordinate of the second family with Q(sigma) = sinh(b (s* - sigma))/sinh(b s*), b = sqrt2 kappa. For sigma well below s* it is
the tortoise coordinate m (Q = e^{-b m}); Q vanishes linearly at sigma = s* (the future horizon) and is negative beyond it, in the interior.
The member in this chart is analytic: Pi = e^{b p} Q = e^{2 kappa x} (Kruskal product), rho(Pi) smooth through the horizon, e^{-f} = -(Q'/b) e^{bp} E (A/delta),
with delta = rho - rho_h and E = delta/Pi (lnE' = 1/delta - 2 kappa/A, a regular equation). At c = 0 (no horizon) the chart is the tortoise one.

    linear.py check [--only=C1,..] [--mutate=NAME]       code checks, planted cases, acceptance and location checks (exit 1 when one fails)
    linear.py table                                      the comparison of growth rates at c = -1.5, -1.01, -1, -0.99, -0.3
    linear.py c1                                         c = -1: data classes, horizon sector, projection, profile, seeds
"""
import math
import os
import sys
import time

import numba as nb
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import evolve_potential as ep  # noqa: E402

SQ2 = math.sqrt(2.0)
MUT = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--mutate=")), None)

# parameter vector indices
P_C, P_K, P_RH, P_B, P_KAP, P_A1, P_A2, P_A3, P_A4, P_A5, P_KIND, P_SS, P_SINH, P_RLO, P_DR, P_NT, P_RHOC = range(17)
# mutation multipliers (all 1.0 = the equations as derived)
P_MV, P_MVP, P_MVPP, P_MEPS, P_MCROSS, P_MDE, P_MRM, P_MBK, P_MIRR, P_MCON, P_MVH, P_MISC = range(17, 29)
NPAR = 29


# ----------------------------------------------------------------------------------------------------------------------
# background
# ----------------------------------------------------------------------------------------------------------------------
@nb.njit(cache=True)
def _A(rho, par):
    return 1.0 + par[P_C] * (rho * rho + 1.0) + par[P_K] * (rho + (rho * rho + 1.0) * math.atan(rho))


@nb.njit(cache=True)
def _Ap(rho, par):
    return 2.0 * par[P_C] * rho + 2.0 * par[P_K] * (1.0 + rho * math.atan(rho))


@nb.njit(cache=True)
def _Aod(delta, par):
    """A(rho_h + delta)/delta, series for |delta| < 1e-3."""
    if abs(delta) < 1e-3:
        return par[P_A1] + delta * (par[P_A2] / 2 + delta * (par[P_A3] / 6 + delta * (par[P_A4] / 24 + delta * par[P_A5] / 120)))
    return _A(par[P_RH] + delta, par) / delta


@nb.njit(cache=True)
def _dAod(delta, par):
    """d/d rho of A/delta."""
    if abs(delta) < 1e-3:
        return par[P_A2] / 2 + delta * (par[P_A3] / 3 + delta * (par[P_A4] / 8 + delta * par[P_A5] / 30))
    return (_Ap(par[P_RH] + delta, par) * delta - _A(par[P_RH] + delta, par)) / (delta * delta)


@nb.njit(cache=True)
def _gfun(delta, par):
    """d lnE/d rho = 1/delta - 2 kappa/A, regular at the horizon."""
    if abs(delta) < 1e-3:
        S = par[P_A2] / 2 + delta * (par[P_A3] / 6 + delta * (par[P_A4] / 24 + delta * par[P_A5] / 120))
        return S / (par[P_A1] + delta * S)
    return 1.0 / delta - 2.0 * par[P_KAP] / _A(par[P_RH] + delta, par)


@nb.njit(cache=True)
def _lnE(rho, par, tab):
    """cubic Hermite of lnE on the uniform rho grid; tab[0] = lnE nodes, tab[1] = slopes."""
    x = (rho - par[P_RLO]) / par[P_DR]
    n = int(par[P_NT])
    i = int(math.floor(x))
    if i < 0:
        i = 0
    if i > n - 2:
        i = n - 2
    t = x - i
    h = par[P_DR]
    y0, y1, m0, m1 = tab[0, i], tab[0, i + 1], tab[1, i] * h, tab[1, i + 1] * h
    t2 = t * t
    t3 = t2 * t
    return (2 * t3 - 3 * t2 + 1) * y0 + (t3 - 2 * t2 + t) * m0 + (-2 * t3 + 3 * t2) * y1 + (t3 - t2) * m1


@nb.njit(cache=True)
def _solve_delta(Pi, par, tab, guess):
    """delta with delta = Pi E(rho_h + delta); G(delta) = delta - Pi E is increasing (G' = 2 kappa delta/A > 0): safeguarded Newton."""
    if Pi == 0.0:
        return 0.0
    lo = par[P_RLO] - par[P_RH]
    hi = par[P_RLO] + par[P_DR] * (par[P_NT] - 1) - par[P_RH]
    d = guess
    if d < lo or d > hi:
        d = 0.5 * (lo + hi)
    for _ in range(200):
        rho = par[P_RH] + d
        E = math.exp(_lnE(rho, par, tab))
        G = d - Pi * E
        Gp = 1.0 - Pi * E * _gfun(d, par)
        if G > 0:
            hi = d
        else:
            lo = d
        dn = d - G / Gp if Gp != 0.0 else 0.5 * (lo + hi)
        if not (lo < dn < hi):
            dn = 0.5 * (lo + hi)
        if abs(dn - d) <= 1e-15 * abs(dn) + 1e-300:
            d = dn
            break
        d = dn
    return d


@nb.njit(cache=True)
def bg(p, sig, par, tab, guess):
    """The member at (p, sigma): (rho, delta, rho_p, rho_s, emf = e^{-f}, f_p, f_s)."""
    if par[P_KIND] == 0.0:
        return par[P_RHOC] + (p - sig) / SQ2_NB, 0.0, 1.0 / SQ2_NB, -1.0 / SQ2_NB, 1.0, 0.0, 0.0
    b = par[P_B]
    u = b * (par[P_SS] - sig)
    Q = math.sinh(u) / par[P_SINH]
    Qd = -b * math.cosh(u) / par[P_SINH]
    ebp = math.exp(b * p)
    Pi = ebp * (abs(Q) if par[P_MIRR] != 0.0 else Q) * (par[P_MISC] if Q < 0.0 else 1.0)
    d = _solve_delta(Pi, par, tab, guess)
    rho = par[P_RH] + d
    E = math.exp(_lnE(rho, par, tab))
    aod = _Aod(d, par)
    rhoPi = E * aod / (2.0 * par[P_KAP])
    rp = rhoPi * b * Pi
    rs = rhoPi * ebp * Qd
    emf = -(Qd / b) * ebp * E * aod
    dl = _gfun(d, par) + _dAod(d, par) / aod
    fp = -(b + dl * rp)
    fs = -(-b * math.tanh(u) + dl * rs)
    return rho, d, rp, rs, emf, fp, fs


SQ2_NB = 1.4142135623730951


@nb.njit(cache=True)
def pot3(rho, par):
    """V, V', V'' (derivatives with respect to phi) on the member at rho."""
    psi = math.atan(rho)
    s2 = math.sin(2.0 * psi)
    c2 = math.cos(2.0 * psi)
    c, K = par[P_C], par[P_K]
    V = -c * (2.0 - c2) - K * (1.5 * s2 + psi * (2.0 - c2))
    Vpsi = -2.0 * c * s2 - 2.0 * K * (1.0 + c2 + psi * s2)
    Vpsipsi = -4.0 * c * c2 + 2.0 * K * s2 - 4.0 * K * psi * c2
    return V, Vpsi / SQ2_NB, 0.5 * Vpsipsi


class Chart:
    """The member c in the (p, sigma) chart; c = 0 is the tortoise chart."""

    def __init__(self, c, rho_c=1.0, sig_star=16.0, p_max=14.0, q_int=0.05, mut=None):
        self.mem = ep.Member(c)
        self.c = float(c)
        self.rho_c = float(rho_c)
        par = np.zeros(NPAR)
        par[P_MV:P_MBK + 1] = 1.0
        par[P_MIRR] = 0.0
        par[P_MCON] = 1.0
        par[P_MVH] = 1.0
        par[P_MISC] = 1.0
        par[P_C], par[P_K] = self.mem.c, self.mem.K
        par[P_RHOC] = rho_c
        self.mut = mut or {}
        for k, v in self.mut.items():
            par[k] = v
        if self.c == 0.0:
            par[P_KIND] = 0.0
            self.tab = np.zeros((2, 2))
            self.par = par
            self.b = 0.0
            return
        mem = self.mem
        rh, kap = mem.rh, mem.kappa
        K = mem.K
        s = 1.0 + rh * rh
        par[P_RH], par[P_KAP] = rh, kap
        par[P_B] = kap * SQ2 * par[P_MBK]
        par[P_A1] = float(mem.Ap(rh))
        par[P_A2] = 2 * self.c + 2 * K * (math.atan(rh) + rh / s)
        par[P_A3] = 4 * K / s ** 2
        par[P_A4] = -16 * K * rh / s ** 3
        par[P_A5] = -16 * K * (1 - 5 * rh * rh) / s ** 4
        par[P_KIND] = 1.0
        par[P_SS] = sig_star
        b = par[P_B]
        par[P_SINH] = math.sinh(b * sig_star)
        self.b = b
        # lnE on a uniform rho grid, from the corner (E(rho_c) = rho_c - rho_h since Pi = 1 there)
        from scipy.integrate import solve_ivp
        a1_est = 1.0                                     # only used for the extent of the interior
        rho_hi = rho_c + 1.25 * p_max / SQ2 + 3.0
        # interior depth: Pi_min = e^{b p_max} Q(sigma_max) = -q_int; delta ~ -E(rho_h) q_int, E(rho_h) <= (rho_c - rho_h) e^{..}; take a generous margin
        rho_lo = rh - max(4.0 * q_int, 0.5) * 1.0
        dr = 0.0025
        n = int(math.ceil((rho_hi - rho_lo) / dr)) + 1
        par[P_RLO], par[P_DR], par[P_NT] = rho_lo, dr, n
        self.par = par
        grid = rho_lo + dr * np.arange(n)

        def rhs(r, y):
            return [float(_gfun(r - rh, par))]
        y0 = math.log(rho_c - rh) if rho_c > rh else None
        if y0 is None:
            raise ValueError("the corner must lie in the static region")
        up = solve_ivp(rhs, [rho_c, grid[-1]], [y0], method="DOP853", rtol=1e-13, atol=1e-15, t_eval=grid[grid >= rho_c])
        dn = solve_ivp(rhs, [rho_c, grid[0]], [y0], method="DOP853", rtol=1e-13, atol=1e-15, t_eval=grid[grid < rho_c][::-1])
        lnE = np.concatenate([dn.y[0][::-1], up.y[0]])
        assert len(lnE) == n, (len(lnE), n)
        g = np.array([float(_gfun(r - rh, par)) for r in grid])
        self.tab = np.vstack([lnE, g])

    def sigma_max(self, p_max, q_int=0.05):
        """sigma beyond the horizon at which Pi(p_max) = -q_int."""
        if self.c == 0.0:
            return None
        b, ss = self.b, self.par[P_SS]
        return ss + math.asinh(q_int * math.exp(-b * p_max) * math.sinh(b * ss)) / b

    def sigma_of_m(self, m):
        """sigma with Q(sigma) = e^{-b m}."""
        if self.c == 0.0:
            return np.asarray(m, float)
        b, ss = self.b, self.par[P_SS]
        return ss - np.arcsinh(np.exp(-b * np.asarray(m, float)) * math.sinh(b * ss)) / b

    def m_of_sigma(self, sig):
        if self.c == 0.0:
            return np.asarray(sig, float)
        b, ss = self.b, self.par[P_SS]
        Q = np.sinh(b * (ss - np.asarray(sig, float))) / math.sinh(b * ss)
        return -np.log(Q) / b


# ----------------------------------------------------------------------------------------------------------------------
# the linear system U_pm = B0 U + Bp U_p + Bm U_m, U = (dr, df, dphi); m stands for sigma
# ----------------------------------------------------------------------------------------------------------------------
@nb.njit(cache=True)
def coef(p, sig, par, tab, guess, B0, Bp, Bm):
    rho, d, rho_p, rho_s, E, fp, fs = bg(p, sig, par, tab, guess)
    R = math.sqrt(rho * rho + 1.0)
    Rp = rho / R
    r = R
    rp = Rp * rho_p
    rm = Rp * rho_s
    Phip = SQ2_NB / (1.0 + rho * rho)
    php = Phip * rho_p
    phm = Phip * rho_s
    V0, V1, V2 = pot3(rho, par)
    mv, mvp, mvpp, meps, mcross, mde, mrm = par[P_MV], par[P_MVP], par[P_MVPP], par[P_MEPS], par[P_MCROSS], par[P_MDE], par[P_MRM]
    V0 *= mv
    V1 *= mvp
    V2 *= mvpp
    for i in range(3):
        for j in range(3):
            B0[i, j] = 0.0
            Bp[i, j] = 0.0
            Bm[i, j] = 0.0
    r2 = r * r
    # row dr
    B0[0, 0] = E * V0 / 2 + E / (2 * r2) + rm * rp / r2
    B0[0, 1] = mde * (-E * V0 * r / 2 + E / (2 * r))
    B0[0, 2] = E * V1 * r / 2
    Bp[0, 0] = -mrm * rm / r
    Bm[0, 0] = -mrm * rp / r
    # row df
    B0[1, 0] = 2 * E / (r2 * r) + 4 * rm * rp / (r2 * r)
    B0[1, 1] = mde * E / r2
    Bp[1, 0] = -2 * rm / r2
    Bm[1, 0] = -2 * rp / r2
    Bp[1, 2] = -meps * phm
    Bm[1, 2] = -meps * php
    # row dphi
    B0[2, 0] = (phm * rp + php * rm) / r2
    B0[2, 1] = -mde * meps * E * V1 / 2
    B0[2, 2] = meps * E * V2 / 2
    Bp[2, 0] = -phm / r
    Bm[2, 0] = -php / r
    Bp[2, 2] = -mcross * rm / r
    Bm[2, 2] = -mcross * rp / r
    # the V1 term of the dr row carries no eps; the V', V'' terms of the dphi row carry eps (meps) through the equation
    return d


@nb.njit(cache=True)
def _solve3(M, rhs, out):
    """Gaussian elimination with partial pivoting, 3 x 3."""
    a = np.empty((3, 4))
    for i in range(3):
        for j in range(3):
            a[i, j] = M[i, j]
        a[i, 3] = rhs[i]
    for k in range(3):
        piv = k
        for i in range(k + 1, 3):
            if abs(a[i, k]) > abs(a[piv, k]):
                piv = i
        if piv != k:
            for j in range(4):
                t = a[k, j]
                a[k, j] = a[piv, j]
                a[piv, j] = t
        for i in range(k + 1, 3):
            f = a[i, k] / a[k, k]
            for j in range(k, 4):
                a[i, j] -= f * a[k, j]
    for i in range(2, -1, -1):
        s = a[i, 3]
        for j in range(i + 1, 3):
            s -= a[i, j] * out[j]
        out[i] = s / a[i, i]


@nb.njit(cache=True)
def _cell(Sv, Ev, Wv, p, sig, hp, hs, par, tab, guess, B0, Bp, Bm, M, rhs, out):
    """N = E + W - S + hp hs (B0 U + Bp U_p + Bm U_s) with the cell-centre averages, U_p = (E + N - S - W)/(2 hp), U_s = (W + N - S - E)/(2 hs)."""
    d = coef(p + 0.5 * hp, sig + 0.5 * hs, par, tab, guess, B0, Bp, Bm)
    H = hp * hs
    for i in range(3):
        s = Ev[i] + Wv[i] - Sv[i]
        for j in range(3):
            M[i, j] = (1.0 if i == j else 0.0) - H * B0[i, j] / 4 - 0.5 * hs * Bp[i, j] - 0.5 * hp * Bm[i, j]
            s += H * (B0[i, j] * (Sv[j] + Ev[j] + Wv[j]) / 4 + Bp[i, j] * (Ev[j] - Sv[j] - Wv[j]) / (2 * hp)
                      + Bm[i, j] * (Wv[j] - Sv[j] - Ev[j]) / (2 * hs))
        rhs[i] = s
    _solve3(M, rhs, out)
    return d


@nb.njit(cache=True)
def march_full(UL, UR, par, tab, hp, hs):
    """The whole rectangle. UL[:, i] = data on p = 0 (sigma_i = i hs), UR[:, j] = data on sigma = 0 (p_j = j hp). Returns U[:, j, i]."""
    Ns = UL.shape[1] - 1
    Np = UR.shape[1] - 1
    U = np.zeros((3, Np + 1, Ns + 1))
    for q in range(3):
        for i in range(Ns + 1):
            U[q, 0, i] = UL[q, i]
        for j in range(Np + 1):
            U[q, j, 0] = UR[q, j]
    B0 = np.zeros((3, 3))
    Bp = np.zeros((3, 3))
    Bm = np.zeros((3, 3))
    M = np.zeros((3, 3))
    rhs = np.zeros(3)
    out = np.zeros(3)
    Sv = np.zeros(3)
    Ev = np.zeros(3)
    Wv = np.zeros(3)
    for j in range(Np):
        guess = 0.0
        for i in range(Ns):
            for q in range(3):
                Sv[q] = U[q, j, i]
                Ev[q] = U[q, j + 1, i]
                Wv[q] = U[q, j, i + 1]
            guess = _cell(Sv, Ev, Wv, j * hp, i * hs, hp, hs, par, tab, guess, B0, Bp, Bm, M, rhs, out)
            for q in range(3):
                U[q, j + 1, i + 1] = out[q]
    return U


@nb.njit(cache=True)
def Ynode(U3, p, sig, par, tab, guess):
    """Y = dr - (rho R/sqrt2) dphi at a node."""
    rho, d, a, b_, c_, e_, f_ = bg(p, sig, par, tab, guess)
    R = math.sqrt(rho * rho + 1.0)
    return U3[0] - rho * R / SQ2_NB * U3[2], d


@nb.njit(cache=True)
def march_obs(UL, UR, par, tab, hp, hs, obs_sig, jh, krow, keep_rows):
    """As march_full but only two rows are held. Records, per row j, Y at the observers (cubic Lagrange in sigma through four nodes around obs_sig[k, j]),
    Y at the node jh (the horizon node), and U on every krow-th row (keep_rows[(j // krow)])."""
    Ns = UL.shape[1] - 1
    Np = UR.shape[1] - 1
    nobs = obs_sig.shape[0]
    Yobs = np.full((nobs, Np + 1), np.nan)
    Yh = np.full(Np + 1, np.nan)
    old = UL.copy()
    new = np.zeros_like(old)
    B0 = np.zeros((3, 3))
    Bp = np.zeros((3, 3))
    Bm = np.zeros((3, 3))
    M = np.zeros((3, 3))
    rhs = np.zeros(3)
    out = np.zeros(3)
    Sv = np.zeros(3)
    Ev = np.zeros(3)
    Wv = np.zeros(3)
    yy = np.zeros(4)
    for j in range(Np + 1):
        p = j * hp
        if j > 0:
            for q in range(3):
                new[q, 0] = UR[q, j]
            guess = 0.0
            for i in range(Ns):
                for q in range(3):
                    Sv[q] = old[q, i]
                    Ev[q] = new[q, i]
                    Wv[q] = old[q, i + 1]
                guess = _cell(Sv, Ev, Wv, (j - 1) * hp, i * hs, hp, hs, par, tab, guess, B0, Bp, Bm, M, rhs, out)
                for q in range(3):
                    new[q, i + 1] = out[q]
            old, new = new, old
        # now `old` holds row j
        if jh >= 0:
            yh, dd = Ynode(old[:, jh], p, jh * hs, par, tab, 0.0)
            Yh[j] = yh
        for k in range(nobs):
            s = obs_sig[k, j]
            if np.isnan(s):
                continue
            i0 = int(math.floor(s / hs)) - 1
            if i0 < 0 or i0 + 3 > Ns:
                continue
            g = 0.0
            for a_ in range(4):
                yy[a_], g = Ynode(old[:, i0 + a_], p, (i0 + a_) * hs, par, tab, g)
            x = s / hs - i0           # in [1, 2]
            Yobs[k, j] = (yy[0] * (-(x - 1) * (x - 2) * (x - 3) / 6) + yy[1] * (x * (x - 2) * (x - 3) / 2)
                          + yy[2] * (-x * (x - 1) * (x - 3) / 2) + yy[3] * (x * (x - 1) * (x - 2) / 6))
        if krow > 0 and j % krow == 0:
            idx = j // krow
            for q in range(3):
                for i in range(Ns + 1):
                    keep_rows[idx, q, i] = old[q, i]
    return Yobs, Yh


# ----------------------------------------------------------------------------------------------------------------------
# data on the two rays: dphi given (a Gaussian pulse), df = 0, dr from the constraint ODE with the corner constants
# ----------------------------------------------------------------------------------------------------------------------
@nb.njit(cache=True)
def _ray_rhs(s, axis, par, tab, guess, y0, y1, amp, s0, w):
    """a'' = r phi' e' - f' a' + phi'^2 a / 2 along p (axis 0, sigma = 0) or along sigma (axis 1, p = 0); e is a Gaussian pulse in the ray coordinate."""
    if axis == 0:
        rho, d, rho_p, rho_s, E, fp, fs = bg(s, 0.0, par, tab, guess)
        rd, fd = rho_p, fp
    else:
        rho, d, rho_p, rho_s, E, fp, fs = bg(0.0, s, par, tab, guess)
        rd, fd = rho_s, fs
    r = math.sqrt(rho * rho + 1.0)
    phd = SQ2_NB / (1.0 + rho * rho) * rd
    x = (s - s0) / w
    ep_ = amp * (-2.0 * x / w) * math.exp(-x * x)
    return y1, r * phd * ep_ - par[P_MCON] * fd * y1 + 0.5 * phd * phd * y0, d


@nb.njit(cache=True)
def ray_data(par, tab, axis, n, h, amp, s0, w, a0, a1, nsub):
    """(3, n+1): dr, df = 0, dphi along a ray with the nodes s = i h; RK4 with nsub sub-steps per node."""
    out = np.zeros((3, n + 1))
    y0, y1 = a0, a1
    guess = 0.0
    out[0, 0] = y0
    xx0 = -s0 / w
    g0 = amp * math.exp(-xx0 * xx0)
    hh = h / nsub
    for i in range(n):
        for k in range(nsub):
            s = i * h + k * hh
            k1a, k1b, guess = _ray_rhs(s, axis, par, tab, guess, y0, y1, amp, s0, w)
            k2a, k2b, guess = _ray_rhs(s + 0.5 * hh, axis, par, tab, guess, y0 + 0.5 * hh * k1a, y1 + 0.5 * hh * k1b, amp, s0, w)
            k3a, k3b, guess = _ray_rhs(s + 0.5 * hh, axis, par, tab, guess, y0 + 0.5 * hh * k2a, y1 + 0.5 * hh * k2b, amp, s0, w)
            k4a, k4b, guess = _ray_rhs(s + hh, axis, par, tab, guess, y0 + hh * k3a, y1 + hh * k3b, amp, s0, w)
            y0 += hh / 6 * (k1a + 2 * k2a + 2 * k3a + k4a)
            y1 += hh / 6 * (k1b + 2 * k2b + 2 * k3b + k4b)
        out[0, i + 1] = y0
        xx = ((i + 1) * h - s0) / w
        out[2, i + 1] = amp * math.exp(-xx * xx) - g0           # dphi vanishes at the corner on both rays
    return out


@nb.njit(cache=True)
def gauge_solution(par, tab, p, sig, xi, xip, eta, etap):
    """The pure-gauge perturbation (L_xi of the member) at (p, sigma): dr = xi r_p + eta r_s, dphi likewise, df = xi f_p + eta f_s - xi' - eta' (xi, eta, and
    their derivatives evaluated at p and sigma by the caller)."""
    rho, d, rho_p, rho_s, E, fp, fs = bg(p, sig, par, tab, 0.0)
    R = math.sqrt(rho * rho + 1.0)
    Rp = rho / R
    Phip = SQ2_NB / (1.0 + rho * rho)
    out = np.zeros(3)
    out[0] = xi * Rp * rho_p + eta * Rp * rho_s
    out[2] = xi * Phip * rho_p + eta * Phip * rho_s
    out[1] = xi * fp + eta * fs - xip - etap
    return out


def gauge_arrays(chart, hp, hs, Np, Ns, xi_f, eta_f):
    """Exact pure-gauge solution on all nodes (3, Np+1, Ns+1) with xi(p), eta(sigma) given as functions returning (value, derivative)."""
    par, tab = chart.par, chart.tab
    U = np.zeros((3, Np + 1, Ns + 1))
    for j in range(Np + 1):
        xi, xip = xi_f(j * hp)
        for i in range(Ns + 1):
            eta, etap = eta_f(i * hs)
            U[:, j, i] = gauge_solution(par, tab, j * hp, i * hs, xi, xip, eta, etap)
    return U


# ----------------------------------------------------------------------------------------------------------------------
# runs
# ----------------------------------------------------------------------------------------------------------------------
class Setup:
    """Grid and data of one run. spec = dict(R=(amp, s0, w), L=(amp, s0, w), a0=dr at the corner, a1R=dr_p at the corner, a1L=dr_sigma at the corner)."""

    def __init__(self, c, h, P, spec, rho_c=1.0, mut=None, sig_star=16.0, q_int=0.05, hp=None, hs=None, jh_shift=0.0):
        self.c, self.h, self.P = float(c), h, P
        self.hp = hp or h
        self.hs = hs or h
        self.chart = Chart(c, rho_c=rho_c, sig_star=sig_star, p_max=P, q_int=q_int, mut=mut)
        ch = self.chart
        self.Np = int(round(P / self.hp))
        if c == 0.0:
            self.sig_max = P + 2.0
            self.jh = -1
        else:
            self.sig_max = ch.sigma_max(P, q_int)
            self.jh = int(round(sig_star / self.hs))
            assert abs(self.jh * self.hs - sig_star) < 1e-9, "sigma* must be a node"
            self.jh -= int(round(jh_shift / self.hs))          # (mutation proj_node: the 'horizon' values are taken jh_shift outside the horizon)
        self.Ns = int(math.ceil(self.sig_max / self.hs - 1e-9))
        a0 = spec.get("a0", 0.0)
        R = spec.get("R", (0.0, 3.0, 0.5))
        Lf = spec.get("L", (0.0, 3.0, 0.5))
        self.UR = ray_data(ch.par, ch.tab, 0, self.Np, self.hp, R[0], R[1], R[2], a0, spec.get("a1R", 0.0), 8)
        self.UL = ray_data(ch.par, ch.tab, 1, self.Ns, self.hs, Lf[0], Lf[1], Lf[2], a0, spec.get("a1L", 0.0), 8)

    def observers(self, xs):
        """sigma of the observers x = const at each row p_j (nan outside the chart)."""
        p = np.arange(self.Np + 1) * self.hp
        out = np.full((len(xs), self.Np + 1), np.nan)
        for k, x in enumerate(xs):
            m = p - SQ2 * x
            ok = m > 0
            sg = np.full_like(m, np.nan)
            sg[ok] = self.chart.sigma_of_m(m[ok])
            sg[(sg > self.sig_star_cut()) | (sg < 0)] = np.nan
            out[k] = sg
        return out

    def sig_star_cut(self):
        return (self.Ns - 3) * self.hs

    def run(self, xs=(0.0,), krow=0):
        obs = self.observers(xs)
        nkeep = (self.Np // krow + 1) if krow > 0 else 1
        keep = np.zeros((nkeep, 3, self.Ns + 1))
        t0 = time.time()
        Yobs, Yh = march_obs(self.UL, self.UR, self.chart.par, self.chart.tab, self.hp, self.hs, obs, self.jh, krow, keep)
        self.wall = time.time() - t0
        p = np.arange(self.Np + 1) * self.hp
        t = np.array([SQ2 * p - x for x in xs])
        return dict(t=t, Y=Yobs, Yh=Yh, p=p, keep=keep)


def window_rate(t, y, t1, t2):
    """Slope of ln|y| against t over [t1, t2] (least squares)."""
    m = (t >= t1) & (t <= t2) & np.isfinite(y) & (y != 0)
    if m.sum() < 4:
        return float("nan")
    return float(np.polyfit(t[m], np.log(np.abs(y[m])), 1)[0])


# ----------------------------------------------------------------------------------------------------------------------
# the 1D master equation of #212 in characteristic form, on the tortoise grid (c <= -1: no throat in the static region)
#     u_pm = -Vhat u / 2,   u = Y/rho,   Vhat = 2A^2/rho^2 - (1 + 1/R^2) A A'/rho + A^2/R^4 - 2A/R^2   (README.md Sec. 2)
# ----------------------------------------------------------------------------------------------------------------------
@nb.njit(cache=True)
def rho_of_x(x, par, tab, guess):
    """rho on the member with rho(0) = rho_c at tortoise x (Pi = e^{2 kappa x}); c = 0: rho_c + x."""
    if par[P_KIND] == 0.0:
        return par[P_RHOC] + x, 0.0
    Pi = math.exp(par[P_B] * SQ2_NB * x)       # 2 kappa x = b sqrt2 x
    d = _solve_delta(Pi, par, tab, guess)
    return par[P_RH] + d, d


@nb.njit(cache=True)
def vhat(rho, par):
    A = _A(rho, par)
    Ap = _Ap(rho, par)
    R2 = rho * rho + 1.0
    return par[P_MVH] * (2 * A * A / (rho * rho) - (1.0 + 1.0 / R2) * A * Ap / rho + A * A / (R2 * R2) - 2 * A / R2)


@nb.njit(cache=True)
def master_march(uL, uR, h, par, tab):
    """Characteristic march of u_pm = -vhat u/2 on a uniform (p, m) grid; uL[k] = u(0, m_k), uR[j] = u(p_j, 0). Returns the diagonal u(p_j, m_j) and the last row."""
    N = len(uR) - 1
    old = uL.copy()
    new = np.zeros_like(old)
    diag = np.zeros(N + 1)
    diag[0] = uL[0]
    for j in range(1, N + 1):
        new[0] = uR[j]
        guess = 0.0
        for k in range(N):
            xc = ((j - 0.5) - (k + 0.5)) * h / SQ2_NB
            rho, guess = rho_of_x(xc, par, tab, guess)
            V = vhat(rho, par)
            S, E, W = old[k], new[k], old[k + 1]
            new[k + 1] = (E + W - S - h * h * V * (S + E + W) / 8.0) / (1.0 + h * h * V / 8.0)
        old, new = new, old
        diag[j] = old[j]
    return diag


def master_from_setup(su):
    """u data from the double-null ray data of a Setup (Y/rho on the two rays, cubic interpolation of the left ray from the sigma nodes to the tortoise nodes),
    marched on the tortoise grid h = hp. Returns t, Y(x = 0) of the 1D equation."""
    ch, par, tab = su.chart, su.chart.par, su.chart.tab
    h, N = su.hp, su.Np
    Y = {}

    def Yray(U, ps, ss):
        out = np.zeros(U.shape[1])
        g = 0.0
        for i in range(U.shape[1]):
            out[i], g = Ynode(U[:, i], ps[i], ss[i], par, tab, g)
        return out
    Np1 = N + 1
    pR = np.arange(Np1) * h
    YR = Yray(su.UR[:, :Np1], pR, np.zeros(Np1))
    sig_nodes = np.arange(su.UL.shape[1]) * su.hs
    YL_nodes = Yray(su.UL, np.zeros(su.UL.shape[1]), sig_nodes)
    m = np.arange(Np1) * h
    sg = ch.sigma_of_m(m)
    from scipy.interpolate import CubicSpline
    YL = CubicSpline(sig_nodes, YL_nodes)(sg)
    uR = np.array([YR[j] / rho_of_x(pR[j] / SQ2, par, tab, 0.0)[0] for j in range(Np1)])
    uL = np.array([YL[k] / rho_of_x(-m[k] / SQ2, par, tab, 0.0)[0] for k in range(Np1)])
    uL[0] = uR[0] = YR[0] / ch.rho_c
    diag = master_march(uL, uR, h, par, tab)
    t = SQ2 * pR
    return t, ch.rho_c * diag


# ----------------------------------------------------------------------------------------------------------------------
# mutations (each breaks one thing; the checks in MUTATIONS must fire; the targets were computed in the development runs dev_80 to dev_99 before they were declared)
# ----------------------------------------------------------------------------------------------------------------------
_COEF = ["C1", "C3", "C5", "R1", "R3", "L1", "L2", "L3a", "L3b", "L5"]
_PHYS = ["C1", "C3", "C5", "P1a", "P1c", "R1", "R3", "L1", "L2", "L3a", "L3b", "L5"]
MUTATIONS = {
    "V_sign": _COEF + ["L3c"],                       # sign of the V terms of the dr row
    "Vp_sign": _COEF + ["L3c"],                      # sign of the V' terms (dr and dphi rows)
    "Vpp_zero": list(_COEF),                         # V'' term of the dphi row dropped
    "eps_sign": _PHYS + ["L3c"],                     # kinetic sign eps in the df cross terms and the V', V'' terms of the dphi row
    "cross_sign": list(_PHYS),                       # sign of the (r_p/r) dphi_m + (r_m/r) dphi_p terms of the dphi row
    "df_zero": list(_PHYS),                          # the df couplings of the dr, df, dphi rows dropped
    "rm_sign": list(_PHYS),                          # sign of the first-derivative terms of the dr row
    "bg_b": ["C2", "C2h", "C3", "C5", "R1", "R3", "L1", "L2", "L3a", "L3b", "L5"],      # b = sqrt2 kappa of the chart 1% off (the chart is not the member)
    "mirror": ["C2", "C3", "C5"],                    # Pi = e^{bp}|Q|: the interior reflected, the chart not smooth at the horizon
    "interior_scale": ["C2", "C3", "C5"],            # Pi scaled by 1.5 inside the horizon
    "con_sign": ["C5", "R1", "R3", "L2", "L3a", "L4", "L5"],                              # sign of the f' dr' term in the ray constraint (the data violate it)
    "vhat_sign": ["L2"],                             # the 1D master equation with V-hat -> -V-hat
    "data_shift": ["L1"],                            # the nonlinear kernel's left ray one node off
    "seed_h": ["L5"],                                # h_p held at 0.01 while h_sigma is halved
    "proj_node": ["L3a"],                            # the 'horizon' values of Y taken one unit of sigma outside the horizon
    "freq_V": ["P1b"],                               # -1/R^2 instead of -2/R^2 in the frequency-domain equation at c = 0
    "quote": ["Q1"],                                 # a quoted slope with an extra digit
    "prof_rate": ["L4"],                             # the profile of L4 divided by e^{1.1 kappa t}
}        # the targets were computed in the development runs dev_92 and dev_93 (negative-mass-raw I232-mac) before they were declared


def chart_mut():
    return {"V_sign": {P_MV: -1.0}, "Vp_sign": {P_MVP: -1.0}, "Vpp_zero": {P_MVPP: 0.0}, "eps_sign": {P_MEPS: -1.0}, "cross_sign": {P_MCROSS: -1.0},
            "df_zero": {P_MDE: 0.0}, "rm_sign": {P_MRM: -1.0}, "bg_b": {P_MBK: 1.01}, "mirror": {P_MIRR: 1.0}, "con_sign": {P_MCON: -1.0}, "vhat_sign": {P_MVH: -1.0}, "interior_scale": {P_MISC: 1.5}}.get(MUT, {})


def mkchart(c, **kw):
    return Chart(c, mut=chart_mut(), **kw)


def mksetup(c, h, P, spec, **kw):
    return Setup(c, h, P, spec, mut=chart_mut(), **kw)


# ----------------------------------------------------------------------------------------------------------------------
# #212's table and #217's slopes, quoted from files of the repository
# ----------------------------------------------------------------------------------------------------------------------
W212 = {-0.3: 0.9700800261, -1.5: 0.4941711675, -1.01: 0.6331476235, -0.99: 0.6401216241}      # README.md Sec. 4 (stability_output.txt)
W2A2 = 1.39705243433350684                                                                      # ../wormhole_largeD/w2_output.txt


def quoted(fname, strings):
    """True when every string occurs in the file of the repository (quoted numbers come from a file)."""
    txt = open(os.path.join(HERE, fname)).read()
    return all(x in txt for x in strings)


SLOPE217 = {-0.3: ("1.03101", 0.9700800261), -1.5: ("2.02338", 0.4941711675), -1.0: ("1.57042", None)}        # evolve_output.txt: primary slope at h = 0.0025


# ----------------------------------------------------------------------------------------------------------------------
# the master equation at c = 0 in the frequency domain (the planted case), written here from README.md Sec. 2
# ----------------------------------------------------------------------------------------------------------------------
def master_freq_c0(wlo=0.3, whi=3.0, rho_max=14.0):
    """Even solution Y = 1 + y2 rho^2 + ... of Y'' = (2/rho) Y' + (1/R^4 - 2/R^2 + W^2) Y at c = 0 (A = 1), shot to rho_max; roots of the sign of Y(rho_max) e^{-W rho_max}
    in W: the growing modes of the Ellis wormhole. Returns the list of W."""
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq

    def f(W):
        r0 = 1e-3
        q0 = W * W - 1.0
        y2 = -q0 / 2.0
        y4 = q0 * y2 / 4.0          # n (n - 3) y_n = q0 y_{n-2}; Q = q0 + O(rho^4) at c = 0
        y0 = [1.0 + y2 * r0 ** 2 + y4 * r0 ** 4, 2.0 * y2 * r0 + 4.0 * y4 * r0 ** 3]
        sol = solve_ivp(lambda r, y: [y[1], 2.0 / r * y[1] + (1.0 / (1 + r * r) ** 2 - (1.0 if MUT == "freq_V" else 2.0) / (1 + r * r) + W * W) * y[0]], [r0, rho_max], y0,
                        method="DOP853", rtol=1e-13, atol=1e-30)
        return sol.y[0, -1] * math.exp(-W * rho_max)
    ws = np.linspace(wlo, whi, 60)
    fv = [f(w) for w in ws]
    roots = []
    for i in range(len(ws) - 1):
        if fv[i] * fv[i + 1] < 0:
            roots.append(brentq(f, ws[i], ws[i + 1], xtol=1e-14, rtol=1e-14))
    return roots


# ----------------------------------------------------------------------------------------------------------------------
# checks
# ----------------------------------------------------------------------------------------------------------------------
RESULTS = []
ONLY = next((a.split("=", 1)[1].split(",") for a in sys.argv if a.startswith("--only=")), None)


def want(cid):
    return ONLY is None or cid in ONLY or any(cid.startswith(o) for o in ONLY)


def check(cid, text, ok, detail=""):
    RESULTS.append((cid, bool(ok)))
    print(f"[{'pass' if ok else 'FAIL'}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def jets(chart, p, sig):
    par, tab = chart.par, chart.tab
    rho, d, rho_p, rho_s, E, fp, fs = bg(p, sig, par, tab, 0.0)
    R = math.sqrt(rho * rho + 1.0)
    Rp = rho / R
    Phip = SQ2 / (1.0 + rho * rho)
    V0, V1, V2 = pot3(rho, par)
    return dict(rho=rho, r=R, rp=Rp * rho_p, rm=Rp * rho_s, E=E, php=Phip * rho_p, phm=Phip * rho_s, V0=V0, V1=V1, V2=V2, fp=fp, fs=fs, rho_p=rho_p, rho_s=rho_s)


def check_C1():
    """coef() against the sympy derivation of linear_derive.py at random points of the chart."""
    argv = sys.argv
    sys.argv = [argv[0]]
    try:
        import linear_derive as LD
    finally:
        sys.argv = argv
    rng = np.random.default_rng(232)
    worst = 0.0
    for c in (-1.5, -1.0, -0.3, 0.0):
        ch = mkchart(c, sig_star=16.0, p_max=6.0)
        par, tab = ch.par, ch.tab
        smax = ch.sigma_max(6.0) if c != 0 else 8.0
        for _ in range(10):
            p = rng.uniform(0.0, 6.0)
            sg = rng.uniform(0.0, smax - 0.05) if c != 0 else rng.uniform(0.0, 8.0)
            if c != 0 and rng.random() < 0.3:
                sg = ch.par[P_SS] + rng.uniform(-0.2, 0.0)
            B0 = np.zeros((3, 3)); Bp = np.zeros((3, 3)); Bm = np.zeros((3, 3))
            coef(p, sg, par, tab, 0.0, B0, Bp, Bm)
            jt = jets(ch, p, sg)
            # the unmutated parameter vector gives the V's; the derivation's matrices use the same V0, V1, V2
            Bd = LD.derived_matrices(jt["r"], jt["rp"], jt["rm"], jt["E"], jt["php"], jt["phm"], jt["V0"], jt["V1"], jt["V2"])
            for X, Y_ in zip((B0, Bp, Bm), Bd):
                worst = max(worst, float(np.max(np.abs(X - Y_) / (1.0 + np.abs(Y_)))))
    check("C1", "coef() equals the sympy-derived matrices B0, Bp, Bm at 40 random points of the chart (c = -1.5, -1, -0.3, 0; interior and horizon included)",
          worst < 1e-12, f"max relative difference {worst:.1e}")


def _fields(ch, p, sig):
    rho, d, rho_p, rho_s, E, fp, fs = bg(p, sig, ch.par, ch.tab, 0.0)
    return math.sqrt(rho * rho + 1.0), -math.log(E), SQ2 * math.atan(rho), rho, E


def check_C2():
    """The member in the chart satisfies the nonlinear equations (finite differences of the chart's closed forms, Richardson in the step), through the horizon."""
    worst = 0.0
    worst_h = 0.0
    for c in (-1.5, -1.0, -0.3):
        ch = mkchart(c, sig_star=16.0, p_max=8.0)
        par, tab = ch.par, ch.tab
        for p in (0.5, 3.0, 7.0):
            for sg in (1.0, 9.0, 15.5, 16.0, 16.06):
                res = []
                for h in (0.01, 0.005):
                    F = lambda dp_, ds: _fields(ch, p + dp_ * h, sg + ds * h)
                    f0 = F(0, 0)
                    fpp, fmm, fpm_, fmp_ = F(1, 0), F(-1, 0), F(0, 1), F(0, -1)
                    a = [F(1, 1), F(1, -1), F(-1, 1), F(-1, -1)]
                    out = []
                    for q in range(3):
                        Dp_ = (fpp[q] - fmm[q]) / (2 * h)
                        Ds = (fpm_[q] - fmp_[q]) / (2 * h)
                        Dpp = (fpp[q] - 2 * f0[q] + fmm[q]) / h ** 2
                        Dss = (fpm_[q] - 2 * f0[q] + fmp_[q]) / h ** 2
                        Dps = (a[0][q] - a[1][q] - a[2][q] + a[3][q]) / (4 * h * h)
                        out.append((Dp_, Ds, Dpp, Dss, Dps))
                    r_, f_, ph_, rho_, E_ = f0
                    V0, V1, V2 = pot3(rho_, ch.par)
                    (rp_, rs_, rpp_, rss_, rps_), (fp_, fs_, fpp_, fss_, fps_), (pp_, ps_, ppp_, pss_, pps_) = out
                    res.append([rps_ - (-(rp_ * rs_ + E_ / 2) / r_ + 0.5 * V0 * r_ * E_),
                                fps_ - (-2 * (rp_ * rs_ + E_ / 2) / r_ ** 2 + EPS_ * pp_ * ps_),
                                pps_ - (-(rp_ * ps_ + rs_ * pp_) / r_ - 0.5 * EPS_ * E_ * V1),
                                rpp_ + fp_ * rp_ + 0.5 * r_ * EPS_ * pp_ ** 2,
                                rss_ + fs_ * rs_ + 0.5 * r_ * EPS_ * ps_ ** 2])
                res = np.array(res)
                rich = (4 * res[1] - res[0]) / 3
                worst = max(worst, float(np.max(np.abs(rich))))
        # the horizon node: r_p = phi_p = 0, f_p = -b
        for p in (1.0, 5.0):
            jt = jets(ch, p, ch.par[P_SS])
            worst_h = max(worst_h, abs(jt["rp"]), abs(jt["php"]), abs(jt["fp"] + SQ2 * ch.mem.kappa))
    check("C2", "the member in the chart solves the three evolution equations and both constraints, across the horizon (interior nodes included): c = -1.5, -1, -0.3, 27 points each",
          worst < 1e-6, f"max residual {worst:.1e}")
    check("C2h", "on the horizon node r_p = phi_p = 0 and f_p = -sqrt2 kappa of the member (the Raychaudhuri equation there is dr_pp - b dr_p = 0)", worst_h < 1e-9, f"max {worst_h:.1e}")


EPS_ = -1.0


def check_C3():
    """The linear scheme reproduces exact pure-gauge solutions at second order (through the horizon)."""
    out = []
    ok = True
    for c in (-1.5, -1.0, -0.3, 0.0):
        errs = []
        for h in (0.04, 0.02, 0.01):
            ch = mkchart(c, rho_c=1.0, sig_star=16.0, p_max=3.0)
            Np = int(round(3.0 / h))
            smax = ch.sigma_max(3.0) if c != 0 else 6.0
            Ns = int(math.ceil(smax / h - 1e-9)) if c != 0 else int(round(6.0 / h))
            xi_f = lambda p: (0.01 * math.sin(0.7 * p + 0.3), 0.007 * math.cos(0.7 * p + 0.3))
            eta_f = lambda s_: (0.02 * math.cos(0.5 * s_) + 0.01, -0.01 * math.sin(0.5 * s_))
            Uex = gauge_arrays(ch, h, h, Np, Ns, xi_f, eta_f)
            U = march_full(np.ascontiguousarray(Uex[:, 0, :]), np.ascontiguousarray(Uex[:, :, 0]), ch.par, ch.tab, h, h)
            errs.append(float(np.abs(U - Uex).max()))
        ratios = [errs[0] / errs[1], errs[1] / errs[2]]
        good = all(3.5 < q < 4.5 for q in ratios) and errs[2] < 3e-4
        ok = ok and good
        out.append(f"c = {c}: {errs[0]:.2e}, {errs[1]:.2e}, {errs[2]:.2e} (ratios {ratios[0]:.2f}, {ratios[1]:.2f})")
    check("C3", "the scheme on exact pure-gauge solutions (xi(p), eta(sigma) smooth) is second order: errors at h = 0.04, 0.02, 0.01, ratios in [3.5, 4.5], error at 0.01 below 3e-4",
          ok, "; ".join(out))


def check_C5():
    """Both constraints stay at truncation level and converge at second order along the evolution."""
    out = []
    ok = True
    for c in (-1.5, -1.0):
        vals = []
        for h in (0.02, 0.01):
            su = mksetup(c, h, 6.0, dict(R=(1e-3, 3.0, 0.5), L=(1e-3, 3.0, 0.5)), q_int=0.05)
            U = march_full(su.UL, su.UR, su.chart.par, su.chart.tab, su.hp, su.hs)
            par, tab = su.chart.par, su.chart.tab
            Cpp = []
            Css = []
            sc = 1e-300
            stp = int(round(0.28 / h))          # the same physical nodes at both resolutions
            for j in range(stp, su.Np - 3, stp):
                for i in range(stp, su.Ns - 3, stp):
                    jt = jets(su.chart, j * su.hp, i * su.hs)
                    a, dd, e = U[0], U[1], U[2]
                    ap = (a[j + 1, i] - a[j - 1, i]) / (2 * su.hp)
                    app = (a[j + 1, i] - 2 * a[j, i] + a[j - 1, i]) / su.hp ** 2
                    dp = (dd[j + 1, i] - dd[j - 1, i]) / (2 * su.hp)
                    ep = (e[j + 1, i] - e[j - 1, i]) / (2 * su.hp)
                    asg = (a[j, i + 1] - a[j, i - 1]) / (2 * su.hs)
                    ass = (a[j, i + 1] - 2 * a[j, i] + a[j, i - 1]) / su.hs ** 2
                    ds = (dd[j, i + 1] - dd[j, i - 1]) / (2 * su.hs)
                    es = (e[j, i + 1] - e[j, i - 1]) / (2 * su.hs)
                    Cpp.append(app + jt["rp"] * dp + jt["fp"] * ap + 0.5 * EPS_ * jt["php"] ** 2 * a[j, i] + EPS_ * jt["r"] * jt["php"] * ep)
                    Css.append(ass + jt["rm"] * ds + jt["fs"] * asg + 0.5 * EPS_ * jt["phm"] ** 2 * a[j, i] + EPS_ * jt["r"] * jt["phm"] * es)
                    sc = max(sc, abs(app), abs(jt["fp"] * ap), abs(ass))
            vals.append((float(np.max(np.abs(Cpp))), float(np.max(np.abs(Css))), sc))
        rat = [vals[0][0] / vals[1][0], vals[0][1] / vals[1][1]]
        rel = [vals[1][0] / vals[1][2], vals[1][1] / vals[1][2]]
        good = all(3.3 < q < 4.7 for q in rat) and all(q < 5e-3 for q in rel)
        ok = ok and good
        out.append(f"c = {c}: C_pp {vals[0][0]:.1e} -> {vals[1][0]:.1e}, C_ss {vals[0][1]:.1e} -> {vals[1][1]:.1e} (ratios {rat[0]:.2f}, {rat[1]:.2f}; relative to the largest term {rel[0]:.1e}, {rel[1]:.1e})")
    check("C5", "the linearized constraints C_pp, C_ss on the solution (pulses on both rays) shrink by a factor in [3.3, 4.7] from h = 0.02 to 0.01, and are below 5e-3 of the largest term at 0.01",
          ok, "; ".join(out))


# ---- rates ------------------------------------------------------------------------------------------------------------
def late_rate(c, h, spec, t1, t2, P=16.0, x=0.0):
    su = mksetup(c, h, P, spec)
    r = su.run(xs=(x,))
    return L_rate(r["t"][0], r["Y"][0], t1, t2), su, r


def L_rate(t, y, t1, t2):
    return window_rate(t, y, t1, t2)


def richardson(r1, r2):
    """second order in h, r1 at h and r2 at h/2"""
    return (4.0 * r2 - r1) / 3.0


SPEC_G = dict(R=(1e-3, 3.0, 0.5), L=(1e-3, 3.0, 0.5))
SPEC_R = dict(R=(1e-3, 3.0, 0.5))
SPEC_L = dict(L=(1e-3, 3.0, 0.5))
SPEC_A0 = dict(a0=1e-3)
WIN = (16.0, 20.0)


def rate_ext(c, spec, win=WIN, hs=(0.01, 0.005), P=16.0):
    rs = [late_rate(c, h, spec, win[0], win[1], P=P)[0] for h in hs]
    return richardson(rs[0], rs[1]), rs


def check_P1():
    """Planted: c = 0 (Ellis, V = 0): W^2 a^2 = 1.39705243433350684 from the linear solver, the frequency-domain master equation, and #39's nonlinear kernel."""
    # (a) linear solver, late window, Richardson
    rext, rs = rate_ext(0.0, dict(R=(1e-3, 3.0, 0.5)), win=(10.0, 14.0), P=14.0)
    e1 = abs(rext ** 2 / W2A2 - 1.0)
    check("P1a", "c = 0: the linear double-null solver, rate of Y(x = 0) over t in [10, 14], Richardson of h = 0.01, 0.005: W^2 a^2 within 1e-4 of 1.39705243433350684",
          e1 < 1e-4, f"W^2 = {rext ** 2:.8f}, relative difference {e1:.1e}; raw rates {rs[0]:.6f}, {rs[1]:.6f}")
    # (b) frequency-domain master equation
    roots = master_freq_c0()
    e2 = abs(roots[0] ** 2 / W2A2 - 1.0) if roots else 1.0
    check("P1b", "c = 0: the master equation of #212 (Y'' = (2/rho) Y' + (1/R^4 - 2/R^2 + W^2) Y, even solution shot to rho = 14): one growing mode in W in [0.3, 3], W^2 a^2 within 1e-8",
          len(roots) == 1 and e2 < 1e-8, f"roots {['%.10f' % x for x in roots]}, W^2 {(roots[0] ** 2 if roots else float('nan')):.10f}, relative difference {e2:.1e}")
    # (c) #39's kernel (via evolve.py): the odd response of member +- eps data against the linear solver in the window [6, 10], Richardson over h
    diffs = []
    for h in (0.01, 0.005):
        su = mksetup(0.0, h, 14.0, dict(R=(1e-6, 3.0, 0.5)))
        r = su.run(xs=(0.0,))
        tn, yn, ye = nl_odd(su)
        diffs.append(window_rate(tn, yn, 6.0, 10.0) - window_rate(r["t"][0], r["Y"][0], 6.0, 10.0))
    d0 = richardson(diffs[0], diffs[1])
    W0 = math.sqrt(W2A2)
    err = abs(2.0 * W0 * d0) / W2A2
    check("P1c", "c = 0: #39's kernel (evolve.py) on member +- the same data: the odd response's rate over [6, 10] minus the linear solver's, Richardson of h = 0.01, 0.005, moves W^2 a^2 by less than 1e-4",
          err < 1e-4, f"differences {diffs[0]:+.2e}, {diffs[1]:+.2e}, extrapolated {d0:+.2e}, W^2 shift {err:.1e}")


def nl_odd(su, amp_scale=1.0, tend=None, shift=0):
    """Odd part (+eps data minus -eps data)/2 of the gauge-invariant Y on the diagonal of #217's nonlinear kernel (evolve.py), member +- the data of the linear solver,
    and the even part (the background's own drift)."""
    import evolve as ev
    import dn
    from scipy.interpolate import CubicSpline
    c, h, N = su.c, su.hp, su.Np
    mem = ev.ep.Member(c)
    rc = su.chart.rho_c
    d0 = ev.member_rays(mem, rc, N, h)
    st = ev.Static(mem, rc, N, h)
    par = ev.params(mem)
    m = (np.arange(N + 1) + (1.0 if MUT == "data_shift" else 0.0)) * h
    sg = su.chart.sigma_of_m(m)
    sig_nodes = np.arange(su.UL.shape[1]) * su.hs
    UL_m = np.array([CubicSpline(sig_nodes, su.UL[q])(sg) for q in range(3)])
    UR = su.UR[:, :N + 1]
    out = {}
    for sgn in (+1.0, -1.0):
        d = dn.Data(2, [-1.0], N, h, d0.Sp + sgn * amp_scale * UR, d0.Sm + sgn * amp_scale * UL_m, name="lin")
        res = ev.evolve_dep(d, par, st, stop=False, record=True, tend=tend)
        out[sgn] = res["rec"]
    return out[1.0]["t"], (out[1.0]["Y"] - out[-1.0]["Y"]) / 2.0, (out[1.0]["Y"] + out[-1.0]["Y"]) / 2.0


# ---- acceptance: rates against #212's W ------------------------------------------------------------------------------
def check_R1():
    """Acceptance: at c = -1.5 and -0.3 the linear rate equals #212's W to 1e-3 (also -1.01, -0.99). Data G: Gaussian pulses of amplitude 1e-3 on both rays."""
    rows = []
    ok_req, ok_all = True, True
    for c in (-1.5, -0.3, -1.01, -0.99):
        rext, rs = rate_ext(c, SPEC_G)
        rel = abs(rext / W212[c] - 1.0)
        rows.append(f"c = {c}: {rext:.6f} against W = {W212[c]:.6f} (relative {rel:.1e}; h = 0.01: {rs[0]:.6f}, h = 0.005: {rs[1]:.6f})")
        ok_all = ok_all and rel < 1e-3
        if c in (-1.5, -0.3):
            ok_req = ok_req and rel < 1e-3
    check("R1", "linear double-null solver, Y(x = 0) rate over t in [16, 20], Richardson of h = 0.01, 0.005, equals #212's W to 1e-3 at c = -1.5 and -0.3 (the issue's criterion) and at -1.01, -0.99",
          ok_req and ok_all, "; ".join(rows))


def check_R3():
    """c = -1: the linear solver grows at kappa (data G)."""
    mem = ep.Member(-1.0)
    rext, rs = rate_ext(-1.0, SPEC_G)
    rel = abs(rext / mem.kappa - 1.0)
    check("R3", "c = -1, data G (pulses on both rays): the rate of Y(x = 0) over [16, 20] (Richardson) equals kappa = 2/pi to 1e-3 (it grows; #212 has no mode)",
          rel < 1e-3, f"{rext:.6f} against kappa = {mem.kappa:.6f} (relative {rel:.1e}; raw {rs[0]:.6f}, {rs[1]:.6f})")


def check_Q217():
    ok = quoted("evolve_output.txt", [v[0] + ("9" if MUT == "quote" else "") for v in SLOPE217.values()])
    check("Q1", "the slopes of #217's lifetime law quoted below (1.03101 at c = -0.3, 2.02338 at c = -1.5, 1.57042 at c = -1; h = 0.0025) occur in evolve_output.txt", ok)


# ---- location: kernel ------------------------------------------------------------------------------------------------
def check_L1():
    """The nonlinear kernel's response to the same data equals the linear solver's until the background's own drift matters."""
    out = []
    ok = True
    for c in (-1.5, -1.0):
        for name, spec in (("G", dict(R=(1e-6, 3.0, 0.5), L=(1e-6, 3.0, 0.5))), ("A0", dict(a0=1e-6))):
            su = mksetup(c, 0.01, 14.0, spec)
            r = su.run(xs=(0.0,))
            tn, yn, ye = nl_odd(su)
            tl, yl = r["t"][0], r["Y"][0]
            okm = np.isfinite(yl)
            from scipy.interpolate import CubicSpline
            sp = CubicSpline(tl[okm], yl[okm])
            ratio = yn / sp(tn)
            Tk = tn[np.where(np.abs(ye) <= 5e-4)[0][-1]]
            sel = (tn >= 3.0) & (tn <= Tk)
            dev = float(np.max(np.abs(ratio[sel] - 1.0)))
            ok = ok and dev < 3e-3
            out.append(f"c = {c} {name}: |ratio - 1| <= {dev:.1e} for t in [3, {Tk:.1f}] (drift {abs(ye[np.argmin(abs(tn - Tk))]):.1e} there)")
    check("L1", "#217's nonlinear kernel (evolve.py) on member +- eps data (eps = 1e-6): the odd part of Y on the diagonal equals the linear solver's to 3e-3 from t = 3 until the even part (the untuned member's drift) reaches 5e-4; c = -1.5, -1; data G and a0",
          ok, "; ".join(out))


# ---- location: the master equation is the same equation ---------------------------------------------------------------
def check_L2():
    out = []
    ok = True
    mem1 = ep.Member(-1.0)
    for c in (-1.5, -1.01, -1.0):
        for name, spec, tmax in (("R", SPEC_R, 10.0 if c == -1.0 else 19.0), ("L", SPEC_L, 19.0), ("G", SPEC_G, 19.0)):
            su = mksetup(c, 0.01, 14.0, spec)
            r = su.run(xs=(0.0,))
            t1, y1 = master_from_setup(su)
            t, y = r["t"][0], r["Y"][0]
            sel = (t >= 2.0) & (t <= tmax) & np.isfinite(y)
            dev = float(np.max(np.abs(y1[sel] - y[sel]) / np.maximum.accumulate(np.abs(y[sel]))))        # against the running maximum: Y crosses zero
            ok = ok and dev < 2e-3
            out.append(f"c = {c} {name}: {dev:.1e}")
    check("L2", "the 1D master equation of #212 in characteristic form (u = Y/rho on the two rays of the same data, marched on the tortoise grid) reproduces the double-null Y(x = 0) to 2e-3 of its running maximum for pulse data R, L, G at c = -1.5, -1.01, -1 (R at c = -1 to t = 10; its later values are the seed)",
          ok, "max |difference| / running maximum: " + "; ".join(out))


# ---- location: the horizon sector at c = -1 ---------------------------------------------------------------------------
def horizon_sets(c, h=0.01, P=16.0):
    sets = {}
    for k, s0 in (("A", 2.5), ("B", 3.5), ("C", 4.5)):
        su = mksetup(c, h, P, dict(L=(1e-3, s0, 0.5)), jh_shift=(1.0 if MUT == "proj_node" else 0.0))
        r = su.run(xs=(0.0,))
        sets[k] = (su, r)
    return sets


def c12_from(r, b, p0=0.0, p1=2.0, h=0.01):
    """c1, c2 of Y_h = c1 + c2 exp(b p) from the horizon values at p0 and p1."""
    yh = r["Yh"]
    i0, i1 = int(round(p0 / h)), int(round(p1 / h))
    M = np.array([[1.0, math.exp(b * p0)], [1.0, math.exp(b * p1)]])
    return np.linalg.solve(M, np.array([yh[i0], yh[i1]]))


def projection(c, h=0.01, P=16.0, bscale=1.0):
    mem = ep.Member(c)
    b = mem.kappa * SQ2 * bscale
    sets = horizon_sets(c, h, P)
    cA, cB, cC = (c12_from(sets[k][1], b, h=h) for k in "ABC")
    lam = np.linalg.solve(np.array([[cB[0], cC[0]], [cB[1], cC[1]]]), -cA)
    rA, rB, rC = (sets[k][1] for k in "ABC")
    t = rA["t"][0]
    Yc = rA["Y"][0] + lam[0] * rB["Y"][0] + lam[1] * rC["Y"][0]
    Yhc = rA["Yh"] + lam[0] * rB["Yh"] + lam[1] * rC["Yh"]
    i14 = np.argmin(abs(t - 14.0))
    return dict(t=t, YA=rA["Y"][0], Yc=Yc, Yh=rA["Yh"], Yhc=Yhc, red=abs(rA["Y"][0][i14]) / abs(Yc[i14]), rate_c=window_rate(t, Yc, 10.0, 14.0),
                rate_A=window_rate(t, rA["Y"][0], 10.0, 14.0), sets=sets, b=b, fits=[cA, cB, cC])


def check_L3():
    bscale = 1.0
    P1 = projection(-1.0, bscale=bscale)
    b = P1["b"] / bscale
    # (a) Y_h = c1 + c2 e^{bp} at c = -1: least-squares residual and the late exponent, at two resolutions
    resid, slopes = {}, {}
    for h in (0.01, 0.005):
        rr, ss = [], []
        for s0 in (2.5, 3.5, 4.5):
            su = mksetup(-1.0, h, 16.0, dict(L=(1e-3, s0, 0.5)), jh_shift=(1.0 if MUT == "proj_node" else 0.0))
            r = su.run(xs=(0.0,))
            p, yh = r["p"], r["Yh"]
            m = p >= 1.0
            Mx = np.vstack([np.ones(m.sum()), np.exp(b * p[m])]).T
            co = np.linalg.lstsq(Mx / np.abs(yh[m])[:, None], yh[m] / np.abs(yh[m]), rcond=None)[0]
            rr.append(float(np.max(np.abs(Mx @ co - yh[m]) / np.abs(yh[m]))))
            ss.append(float(np.polyfit(p[p > 8], np.log(np.abs(yh[p > 8])), 1)[0]))
        resid[h], slopes[h] = max(rr), ss
    ratio = resid[0.01] / resid[0.005]
    sl = [richardson(slopes[0.01][k], slopes[0.005][k]) for k in range(3)]
    relb = max(abs(x / b - 1.0) for x in sl)
    check("L3a", "c = -1: the horizon values of Y for three left-ray pulses (sigma0 = 2.5, 3.5, 4.5) obey Y_h(p) = c1 + c2 exp(b p), b = sqrt2 kappa: the least-squares residual over p in [1, 16] is 5e-4 or less at h = 0.005 and falls by 4 (+-0.5) from h = 0.01, and the late exponent of |Y_h| (Richardson) equals b to 1e-5",
          resid[0.005] < 5e-4 and 3.5 < ratio < 4.5 and relb < 1e-5, f"residual {resid[0.01]:.1e} (h = 0.01), {resid[0.005]:.1e} (h = 0.005), ratio {ratio:.2f}; late exponent {sl[0]:.6f}, {sl[1]:.6f}, {sl[2]:.6f} against b = {b:.6f} (relative {relb:.1e})")
    # (b) projection kills the growth at c = -1, does not at c = -1.5
    P2 = projection(-1.5, bscale=bscale)
    kap = ep.Member(-1.0).kappa
    check("L3b", "c = -1: the combination of the three pulses with c1 = c2 = 0 (from Y_h at p = 0 and p = 2) has |Y(x = 0, t = 14)| smaller by a factor of at least 1e3 than one pulse, and does not grow at kappa (rate over [10, 14] below 0.25 kappa)",
          P1["red"] > 1e3 and P1["rate_c"] < 0.25 * kap, f"reduction {P1['red']:.2e}, rate of the combination {P1['rate_c']:.4f}, of one pulse {P1['rate_A']:.4f}, kappa {kap:.4f}; max |Y_h| of the combination {np.max(np.abs(P1['Yhc'])):.1e} against {np.max(np.abs(P1['Yh'])):.1e} for one pulse")
    check("L3c", "control, c = -1.5: the same operation (the form of Y_h at c = -1 with b of c = -1.5) reduces |Y(x = 0, t = 14)| by less than 100 and the rate stays within 15% of W (the horizon sector is specific to c = -1)",
          P2["red"] < 100 and abs(P2["rate_c"] - W212[-1.5]) < 0.15 * W212[-1.5], f"reduction {P2['red']:.2e}, rate of the combination {P2['rate_c']:.4f}, W = {W212[-1.5]:.4f}")


def profile_slope(c, rate, p_row=8.0, h=0.01, P=16.0, lo=2e-4, hi=2e-3):
    su = mksetup(c, h, P, SPEC_G)
    krow = int(round(p_row / h))
    obs = su.observers((0.0,))
    keep = np.zeros((su.Np // krow + 1, 3, su.Ns + 1))
    march_obs(su.UL, su.UR, su.chart.par, su.chart.tab, su.hp, su.hs, obs, su.jh, krow, keep)
    row = keep[1]
    sig = np.arange(su.jh) * su.hs
    Ys = np.zeros(su.jh)
    rho = np.zeros(su.jh)
    g = 0.0
    for i in range(su.jh):
        Ys[i], g = Ynode(row[:, i], p_row, sig[i], su.chart.par, su.chart.tab, g)
        rho[i] = bg(p_row, sig[i], su.chart.par, su.chart.tab, g)[0]
    t = (p_row + su.chart.m_of_sigma(sig)) / SQ2
    delta = rho - su.chart.mem.rh
    Psi = Ys * np.exp(-rate * t)
    m = (delta > lo) & (delta < hi)
    return float(np.polyfit(np.log(delta[m]), np.log(np.abs(Psi[m])), 1)[0]), int(m.sum())


def check_L4():
    mem1, mem15 = ep.Member(-1.0), ep.Member(-1.5)
    s1, n1 = profile_slope(-1.0, mem1.kappa * (1.1 if MUT == "prof_rate" else 1.0))
    s15, n15 = profile_slope(-1.5, W212[-1.5])
    e15 = W212[-1.5] / (2.0 * mem15.kappa)
    check("L4", "late-time profile near the horizon, Y e^{-W t} against delta = rho - rho_h at the row p = 8, delta in [2e-4, 2e-3]: c = -1 has exponent 1/2 (the branch 1 - W/(2 kappa), W = kappa) and not 3/2 (the branch u -> 0, #212's); control c = -1.5 has W/A'(rho_h)",
          abs(s1 - 0.5) < 0.02 and abs(s15 - e15) < 0.02, f"c = -1: {s1:.4f} ({n1} nodes), expected 0.5, ingoing would be 1.5; c = -1.5: {s15:.4f} ({n15} nodes), expected {e15:.4f}")


def check_L5():
    c = -1.0
    yh = []
    for h in (0.02, 0.01, 0.005):
        hp = h if MUT != "seed_h" else 0.01
        su = mksetup(c, h, 14.0, SPEC_R, hp=hp)
        r = su.run(xs=(0.0,))
        yh.append(abs(r["Yh"][-1]))
    q = [yh[0] / yh[1], yh[1] / yh[2]]
    dr = []
    for h in (0.02, 0.01, 0.005):
        su = mksetup(c, h, 14.0, dict(R=(1e-9, 3.0, 0.5)))
        tn, yn, ye = nl_odd(su)
        dr.append(abs(ye[np.argmin(abs(tn - 10.0))]))
    qn = [dr[0] / dr[1], dr[1] / dr[2]]
    check("L5", "c = -1, data R (right-ray pulse, left ray zero, c2 = 0 in the continuum): the horizon value Y_h(p = 14) of the linear solver and the nonlinear kernel's drift of the untuned member at t = 10 both fall by 4 (+-0.3) per halving of h: the growth from truncation seeds is the h^2 content of c2",
          all(3.7 < x < 4.3 for x in q + qn), f"linear Y_h: {yh[0]:.2e}, {yh[1]:.2e}, {yh[2]:.2e} (ratios {q[0]:.2f}, {q[1]:.2f}); nonlinear drift: {dr[0]:.2e}, {dr[1]:.2e}, {dr[2]:.2e} (ratios {qn[0]:.2f}, {qn[1]:.2f})")


ALL_CHECKS = [("C1", check_C1), ("C2", check_C2), ("C3", check_C3), ("C5", check_C5), ("P1", check_P1), ("R1", check_R1), ("R3", check_R3),
              ("Q1", check_Q217), ("L1", check_L1), ("L2", check_L2), ("L3", check_L3), ("L4", check_L4), ("L5", check_L5)]


def run_checks():
    t0 = time.time()
    for cid, fn in ALL_CHECKS:
        if want(cid):
            t1 = time.time()
            fn()
            print(f"    ({cid}: {time.time() - t1:.0f} s)", flush=True)
    bad = [cid for cid, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS)} checks, {len(RESULTS) - len(bad)} pass; failed: {bad}; mutation {MUT}; {time.time() - t0:.0f} s")
    if MUT is not None and MUT in MUTATIONS:
        want_ = MUTATIONS[MUT]
        print(f"mutation {MUT}: fired {bad}, declared {want_}: {'equal' if sorted(bad) == sorted(want_) else 'DIFFERENT'}")
    return 1 if bad else 0


def master_time_domain(c):
    """#217's evolve_master1d.py (compact Gaussian u data on the t = 0 slice, leapfrog, Dirichlet ends at +-90): rate of ln|u(x = 0)| over t in [45, 60]."""
    import evolve_master1d as em1
    ts, ys, mem = em1.run(c)
    m = (ts >= 45) & (ts <= 60) & (ys > 0)
    return float(np.polyfit(ts[m], np.log(ys[m]), 1)[0])


def nl_untuned_rate(c, h=0.01, win=(8.0, 14.0)):
    su = Setup(c, h, 16.0, dict(R=(1e-9, 3.0, 0.5)))
    tn, yn, ye = nl_odd(su)
    return window_rate(tn, ye, *win)


def report_table():
    print("== Comparison of the growth rates of Y on the diagonal (x = 0), units b = 1 ==")
    print("   linear: the linear double-null solver, data G (Gaussian pulses of amplitude 1e-3, centre 3, width 0.5, on both rays), rate over t in [16, 20], Richardson of h = 0.01 and 0.005")
    print("   1D char: #212's master equation in characteristic form, u = Y/rho of the same data (c <= -1 only: no throat in the static region)")
    print("   1D time: evolve_master1d.py (compact data on a t = 0 slice), rate over [45, 60]")
    print("   nl slope: #217's lifetime law, 1/slope at h = 0.0025 (evolve_output.txt); nl seed: the untuned member's own growth of Y, rate over t in [8, 14], h = 0.01")
    print("   nl same: #217's kernel on member +- the same data (eps = 1e-6): odd part, rate over [6, 10] at h = 0.01 minus the linear solver's over the same window (not asymptotic)")
    hdr = f"{'c':>7} {'kappa':>9} {'#212 W':>10} {'1D time':>9} {'1D char':>9} {'linear':>9} {'lin/W-1':>9} {'nl slope':>9} {'nl seed':>9} {'nl same - lin':>14}"
    print(hdr)
    for c in (-1.5, -1.01, -1.0, -0.99, -0.3):
        mem = ep.Member(c)
        w = W212.get(c)
        rl, _ = rate_ext(c, SPEC_G)
        if c <= -1.0:
            rcs = []
            for h in (0.01, 0.005):
                su = Setup(c, h, 16.0, SPEC_G)
                t1, y1 = master_from_setup(su)
                rcs.append(window_rate(t1, y1, *WIN))
            rc = richardson(*rcs)
            r1d = master_time_domain(c)
        else:
            rc = r1d = float("nan")
        sl = SLOPE217.get(c)
        rslope = 1.0 / float(sl[0]) if sl else float("nan")
        rseed = nl_untuned_rate(c)
        su = Setup(c, 0.01, 14.0, dict(R=(1e-6, 3.0, 0.5), L=(1e-6, 3.0, 0.5)))
        r = su.run(xs=(0.0,))
        tn, yn, ye = nl_odd(su)
        dsame = window_rate(tn, yn, 6.0, 10.0) - window_rate(r["t"][0], r["Y"][0], 6.0, 10.0)
        wtxt = f"{w:10.6f}" if w else f"{'none':>10}"
        rel = f"{rl / w - 1:9.1e}" if w else f"{rl / mem.kappa - 1:8.1e}k"
        print(f"{c:7.2f} {mem.kappa:9.6f} {wtxt} {r1d:9.5f} {rc:9.6f} {rl:9.6f} {rel} {rslope:9.5f} {rseed:9.5f} {dsame:14.2e}", flush=True)
    print("   (c = -1: 'lin/W-1' is the relative difference from kappa, marked k; c = -0.99 and -0.3 have a throat in the static region, where u = Y/rho is singular: no 1D rate is computed there)")


def report_c1():
    kap = ep.Member(-1.0).kappa
    print("== c = -1: data classes (linear solver, h = 0.01 and 0.005, rate of Y(x = 0) over [16, 20], Richardson) ==")
    print(f"   kappa = {kap:.6f}; Y_h(0) is the value of Y at the horizon node on the left ray, 'late exponent' that of |Y_h(p)| over p in [8, 16]; 1D/DN: the 1D master equation's Y over the double-null Y at t = 14")
    classes = [("R  right pulse only", SPEC_R), ("L  left pulse only", SPEC_L), ("G  both pulses", SPEC_G), ("A0 corner dr = 1e-3", SPEC_A0),
               ("A1R dr_p(corner) = 1e-3", dict(a1R=1e-3)), ("A1L dr_s(corner) = 1e-3", dict(a1L=1e-3))]
    print(f"{'class':>26} {'Y_h(0)':>11} {'rate Y(x=0)':>12} {'late exponent Y_h':>18} {'rate 1D char':>13} {'1D/DN at t=14':>14}")
    for name, spec in classes:
        rates, exps, rcs = [], [], []
        for h in (0.01, 0.005):
            su = Setup(-1.0, h, 16.0, spec)
            r = su.run(xs=(0.0,))
            rates.append(window_rate(r["t"][0], r["Y"][0], *WIN))
            p, yh = r["p"], r["Yh"]
            exps.append(float(np.polyfit(p[p > 8], np.log(np.abs(yh[p > 8])), 1)[0]) if np.all(yh[p > 8] != 0) else float("nan"))
            t1, y1 = master_from_setup(su)
            rcs.append(window_rate(t1, y1, *WIN))
            if h == 0.01:
                yh0 = yh[0]
                i = np.argmin(abs(r["t"][0] - 14.0))
                ratio = y1[i] / r["Y"][0][i]
        print(f"{name:>26} {yh0:11.3e} {richardson(*rates):12.5f} {richardson(*exps):18.5f} {richardson(*rcs):13.5f} {ratio:14.4f}", flush=True)
    print("   (the pulse classes R, L, G have the static term M = 0 and the 1D equation reproduces them; the corner constants carry M != 0 (linear_derive.py C2))")
    print(f"   b = sqrt2 kappa = {kap * SQ2:.6f}: a growing Y_h(p) = c2 exp(b p) is the horizon sector")
    P1 = projection(-1.0)
    print("== c = -1: projection of the horizon sector (three left pulses at sigma0 = 2.5, 3.5, 4.5; c1, c2 of Y_h = c1 + c2 exp(b p) cancelled from Y_h at p = 0, 2) ==")
    print(f"   c1, c2 of the three: " + "; ".join(f"({x[0]:+.3e}, {x[1]:+.3e})" for x in P1["fits"]))
    for tt in (6, 10, 14, 18):
        i = np.argmin(abs(P1["t"] - tt))
        print(f"   t = {P1['t'][i]:5.2f}: |Y(x = 0)| one pulse {abs(P1['YA'][i]):.3e}, combination {abs(P1['Yc'][i]):.3e}; Y_h one pulse {P1['Yh'][i]:+.3e}, combination {P1['Yhc'][i]:+.3e}")
    print(f"   rate over [10, 14]: one pulse {P1['rate_A']:.4f}, combination {P1['rate_c']:.4f}; the combination's Y is also the 1D master equation's (u -> 0 at the horizon): it decays")
    s1, n1 = profile_slope(-1.0, kap)
    print(f"== late-time profile near the horizon (row p = 8, delta = rho - rho_h in [2e-4, 2e-3]): d ln|Y e^(-W t)| / d ln delta = {s1:.4f} at c = -1 (W = kappa) ==")
    print("   the branches of the Y equation there: 1 - W/(2 kappa) = 0.5 (regular, Y_h != 0) and 1 + W/(2 kappa) = 1.5 (u -> 0, #212)")


def report_mutations(raw, prefix):
    """The fired sets of the mutated runs stored as RAW/<prefix>_<name>.log against MUTATIONS."""
    import re
    print(f"== Mutations: fired sets in {prefix}_<name>.log against the declared ones ==")
    bad = 0
    for name, want_ in MUTATIONS.items():
        f = os.path.join(raw, f"{prefix}_{name}.log")
        if not os.path.exists(f):
            print(f"  {name:15s} no log"); bad += 1
            continue
        m = re.search(r"failed: \[(.*?)\]", open(f).read())
        got = [x.strip().strip("'") for x in m.group(1).split(",") if x.strip()] if m else None
        ok = got is not None and sorted(got) == sorted(want_)
        bad += 0 if ok else 1
        print(f"  {name:15s} fired {got}  declared {want_}  {'equal' if ok else 'DIFFERENT'}")
    print(f"{len(MUTATIONS) - bad} of {len(MUTATIONS)} equal")
    return 1 if bad else 0


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "check"
    if cmd == "check":
        sys.exit(run_checks())
    if cmd == "table":
        report_table()
        sys.exit(0)
    if cmd == "c1":
        report_c1()
        sys.exit(0)
    if cmd == "mutations":
        sys.exit(report_mutations(sys.argv[2], sys.argv[3]))
    sys.exit(f"unknown command {cmd}")

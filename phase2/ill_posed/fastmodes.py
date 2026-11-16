"""
Vectorised float64 versions of the objects in modes.py (the mp functions are the reference; check_fastmodes.py compares them).
Used to build potentials and data on grids of 1e4-1e5 points.  Near r_s, where Vhat = g^2 - g_x subtracts two 1/(x-x_s)^2
poles, the mp value is used (|x - x_s| omega < 0.02).
"""
import functools

import numpy as np
import mpmath as mp

import modes as M


@functools.lru_cache(None)
def _coef(l):
    ch = np.array([float(a) for a in M.chi_coeffs(l)])
    pi = np.array([float(a) for a in M.pi_coeffs(l)])
    return ch, pi


def _horner(cs, r):
    v = np.zeros_like(r)
    for a in cs[::-1]:
        v = v * r + a
    return v


def _der(cs):
    return cs[1:] * np.arange(1, len(cs))


def x_of_r(r):
    u = np.asarray(r, float) / 2
    small = u < 0.05
    ph = np.empty_like(u)
    us = u[small]
    s = np.zeros_like(us)
    for k in range(24, 1, -1):          # u - log1p(u) = sum_{k>=2} (-1)^k u^k / k
        s = s * us + (-1.0)**k / k
    ph[small] = s * us * us
    ph[~small] = u[~small] - np.log1p(u[~small])
    return 2 * ph


def r_of_x(x):
    x = np.asarray(x, float)
    r = 2 * np.sqrt(x) + x
    for _ in range(60):
        dr = (x_of_r(r) - x) * (r + 2) / r
        r = r - dr
        if np.all(np.abs(dr) <= 1e-15 * r):
            break
    return r


def vhat(l, r, omega_x=None):
    """Vhat(r) for an array r (float), mp fallback within 0.02/omega (in x) of r_s."""
    r = np.asarray(r, float)
    ch, _ = _coef(l)
    c = M.cc(l)
    X = _horner(ch, r)
    X1 = _horner(_der(ch), r)
    X2 = _horner(_der(_der(ch)), r)
    f = 1 + 2 / r
    a = X1 / X - c / (c * r - 6)
    da = X2 / X - (X1 / X)**2 + c**2 / (c * r - 6)**2
    g = f * a
    gr = (-2 / r**2) * a + f * da
    V = g * g - f * gr
    w = float(M.omega(l))
    xs = float(M.x_of(M.rs(l).numerator / mp.mpf(M.rs(l).denominator)))
    near = np.abs(x_of_r(r) - xs) * w < 0.02
    for i in np.nonzero(near)[0]:
        V[i] = float(M.Vhat(l, mp.mpf(float(r[i]))))
    return V


def psihat(l, r):
    """e^{-omega x} r^4 Pi/chi (unnormalised), float."""
    r = np.asarray(r, float)
    ch, pi = _coef(l)
    w = float(M.omega(l))
    return np.exp(-w * x_of_r(r)) * r**4 * _horner(pi, r) / _horner(ch, r)


def kappa_c(l, r):
    """d ln phi0/dx = -g, float (valid away from r_s)."""
    r = np.asarray(r, float)
    ch, _ = _coef(l)
    c = M.cc(l)
    a = _horner(_der(ch), r) / _horner(ch, r) - c / (c * r - 6)
    return -(1 + 2 / r) * a

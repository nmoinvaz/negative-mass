"""The potential V(phi) and the static member of the phantom "black universe" family of [BF06], for the double-null evolution (#217).

Importable by other evolution codes (#222): no dependence on dn.py, numpy and numba only (numba optional).

    ds^2 = -A dt^2 + d rho^2 / A + R^2 dOmega^2,  R = sqrt(rho^2 + 1),  phi = sqrt2 arctan(rho),  units b = 1,
    A = 1 + c (rho^2 + 1) + K (rho + (rho^2 + 1) arctan rho),   K = 3m = -2c/pi,
    R_mn - eps d_m phi d_n phi - V(phi) g_mn = 0,  eps = -1,  8 pi G = 1  (notes/black_universe/README.md, evolve_derive.py).

With psi = phi / sqrt2 (psi in (-pi/2, pi/2) on the member):
    V(psi)   = -c (3 - 2 cos^2 psi) - K (3 sin psi cos psi + psi (3 - 2 cos^2 psi))    [BF06] (14)
    dV/dpsi  = -2c sin 2psi - 2K (1 + cos 2psi + psi sin 2psi)
The formulas are used as written for every phi, including |psi| > pi/2 (the grid should not get there; evolve.py reports max |psi|).
V = 0 for c = 0 (the Ellis wormhole). V(psi = pi/2) = 0 and V(-pi/2) = -6c.

    pot(phi, c, K)  ->  (V, V'(phi))      numba-compatible scalar function (pure Python when numba is absent)
"""
import math

import numpy as np

try:
    import numba as nb
    _jit = nb.njit(cache=True)
except Exception:      # numba is optional
    nb = None

    def _jit(f):
        return f

SQ2 = math.sqrt(2.0)


@_jit
def pot(phi, c, K):
    """(V(phi), dV/dphi) for the member (c, K), b = 1."""
    psi = phi / 1.4142135623730951
    s2 = math.sin(2.0 * psi)
    c2 = math.cos(2.0 * psi)
    V = -c * (2.0 - c2) - K * (1.5 * s2 + psi * (2.0 - c2))
    Vpsi = -2.0 * c * s2 - 2.0 * K * (1.0 + c2 + psi * s2)
    return V, Vpsi / 1.4142135623730951


class Member:
    """The static member c of the family: background functions in double precision, horizon and surface gravity.

    rho_h = the root of A on the side where it exists (c < 0), kappa = A'(rho_h)/2, None for c >= 0 (A > 0 everywhere at c = 0)."""

    def __init__(self, c):
        self.c = float(c)
        self.K = 0.0 if self.c == 0.0 else -2.0 * self.c / math.pi
        self.rh = self.kappa = None
        if self.c < 0.0:
            self.rh = self._horizon()
            self.kappa = 0.5 * float(self.Ap(self.rh))

    def A(self, r):
        return 1.0 + self.c * (r * r + 1.0) + self.K * (r + (r * r + 1.0) * np.arctan(r))

    def Ap(self, r):
        return 2.0 * self.c * r + 2.0 * self.K * (1.0 + r * np.arctan(r))

    def R(self, r):
        return np.sqrt(r * r + 1.0)

    def phi(self, r):
        return SQ2 * np.arctan(r)

    def V(self, r):
        """V along the member, as a function of rho (not phi)."""
        return np.array([pot(SQ2 * math.atan(float(x)), self.c, self.K)[0] for x in np.atleast_1d(r)])

    def _horizon(self):
        from scipy.optimize import brentq
        if self.c == -1.0:
            return 0.0
        lo = -1.0
        while self.A(lo) > 0:
            lo *= 2.0
        return brentq(self.A, lo, self.K + 1.0, xtol=1e-15, rtol=1e-15)

    def a_over_z(self, z):
        """A(rh + z) / z, by a Taylor series for small |z| (no cancellation at the horizon)."""
        z = np.asarray(z, float)
        out = np.empty_like(z)
        sm = np.abs(z) < 1e-3
        rh, K, c, A1 = self.rh, self.K, self.c, 2.0 * self.kappa
        s = 1.0 + rh * rh
        A2 = 2.0 * c + 2.0 * K * (math.atan(rh) + rh / s)
        A3 = 4.0 * K / s ** 2
        A4 = -16.0 * K * rh / s ** 3
        A5 = -16.0 * K * (1.0 - 5.0 * rh * rh) / s ** 4
        zs = z[sm]
        out[sm] = A1 + A2 * zs / 2 + A3 * zs ** 2 / 6 + A4 * zs ** 3 / 24 + A5 * zs ** 4 / 120
        zl = z[~sm]
        out[~sm] = self.A(rh + zl) / zl
        return out

"""
Shared definitions for the Dotti-Gleiser (arXiv:0809.3615) transformed even-parity operator
    H^ = -d^2/dx^2 + V^(x),   V^ = V - 2 g' = g^2 - g',   g = psi0'/psi0   (' = d/dx),
on negative-mass Schwarzschild, units |M| = 1 (M = -1), so f = 1 + 2/r, dx/dr = 1/f = r/(r+2),
x = r - 2 ln(1 + r/2)  (x = 0 at r = 0).

psi0 = chi(r)/(6M + c r), c = (l+2)(l-1), chi = DG eq. (55) (same code as ../vhat_negativity.py).
Since psi0 is a zero-energy solution of the Zerilli operator, V = g' + g^2 and V^ = g^2 - g', i.e.
formally H^ = (-d/dx + g)(d/dx + g).  Hence
    phi0 = 1/psi0 = (c r - 6)/chi
is an exact, smooth (chi > 0 on r > 0) zero-energy solution of H^.  It decays like r^-(l+1) at
infinity (recessive there), behaves like -4/r ~ x^(-1/2) at r -> 0 (the non-admissible branch), and
has exactly one zero, at r_s = 6/c.

In r the equation H^ Psi = E Psi reads
    -f (f Psi_r)_r + V^ Psi = E Psi,
and, with Psi = r^3 u (r^3 ~ 8 x^(3/2): the admissible branch at x = 0),
    -r u'' + (-6 + 2/(r+2)) u' + C(r) u = E r^3/(r+2)^2 u,
    C(r) = r^3 V^/(r+2)^2 - 6(r+1)/(r(r+2))           (regular at r = 0; built exactly in sympy).
"""
import functools
import numpy as np
import sympy as sp

rs_ = sp.Symbol('r', positive=True)


def chi_poly(l):
    M = sp.Integer(-1)
    s = 0
    for n in range(1, l + 3):
        c = sp.Integer(n - 2) * ((n - 4) * l * (l + 1) + n - 1) * sp.gamma(l + n - 1) \
            / (2**n * sp.gamma(n)**2 * sp.gamma(l - n + 3) * (-M)**n)
        s += c * rs_**n
    return sp.expand(s)


def kval(l):
    return (l - 1) * l * (l + 1) * (l + 2) // 6


def r_s(l):
    return 6.0 / ((l + 2) * (l - 1))


@functools.lru_cache(None)
def symbolic(l):
    r = rs_
    c = (l + 2) * (l - 1)
    f = 1 + 2 / r
    chi = chi_poly(l)
    psi0 = chi / (-6 + c * r)
    g = sp.cancel(f * sp.diff(psi0, r) / psi0)
    Vh = sp.cancel(g**2 - f * sp.diff(g, r))
    C = sp.cancel(r**3 * Vh / (r + 2)**2 - 6 * (r + 1) / (r * (r + 2)))
    phi0 = (c * r - 6) / chi
    dlogphi0_dx = sp.cancel(f * sp.diff(phi0, r) / phi0)
    return dict(chi=chi, Vh=Vh, C=C, phi0=phi0, dlogphi0_dx=dlogphi0_dx)


def _ratfun(expr):
    """Return float evaluator for a rational function via polynomial coefficients (Horner)."""
    num, den = sp.fraction(sp.cancel(expr))
    pn = [float(a) for a in sp.Poly(num, rs_).all_coeffs()]
    pd = [float(a) for a in sp.Poly(den, rs_).all_coeffs()]
    return lambda r: np.polyval(pn, r) / np.polyval(pd, r)


@functools.lru_cache(None)
def numeric(l):
    s = symbolic(l)
    return dict(Vh=_ratfun(s['Vh']), C=_ratfun(s['C']), phi0=_ratfun(s['phi0']),
                kappa_c=_ratfun(s['dlogphi0_dx']), chi=_ratfun(s['chi']))


def x_of_r(r):
    r = np.asarray(r, dtype=float)
    return r - 2 * np.log1p(r / 2)


def r_of_x(x):
    """Invert x = r - 2 ln(1+r/2) by Newton (vectorised)."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    r = np.where(x < 1, 2 * np.sqrt(x), x + 2 * np.log1p(x / 2) + 1)
    for _ in range(60):
        F = x_of_r(r) - x
        r = np.maximum(r - F * (r + 2) / r, 1e-300)
    return r


def exact_mode(l, r):
    """DG unstable mode mapped by the intertwiner, Psi^ = A^ Psi_DG, for l = 2:
    Psi^ ~ r^3 (r+2)^4 e^{-2r}/(r^3 - 3 r^2 + 6)  (|M| = 1)."""
    assert l == 2
    r = np.asarray(r, float)
    return r**3 * (r + 2)**4 * np.exp(-2 * r) / (r**3 - 3 * r**2 + 6)


def cheb(N):
    """Chebyshev-Gauss-Lobatto differentiation matrix (Trefethen), points t_j = cos(pi j/N)."""
    j = np.arange(N + 1)
    t = np.cos(np.pi * j / N)
    cc = np.ones(N + 1); cc[0] = cc[-1] = 2
    cc = cc * (-1.0)**j
    T = np.tile(t, (N + 1, 1)).T
    dT = T - T.T
    D = np.outer(cc, 1 / cc) / (dT + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    return D, t

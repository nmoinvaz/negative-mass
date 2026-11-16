"""
Exact modes of the Dotti-Gleiser operator, in mpmath.  Units |M| = 1, f = 1 + 2/r, x = r - 2 ln(1 + r/2), d/dx = f d/dr.

This is a re-implementation from the formulas of phase2/dg_proof/RESULTS.md (Lemma A for chi, L4 for the mode).  It shares
with phase2/dg_spectrum/common.py only the formula for chi (DG eq. 55 = Lemma A); nothing else.

  c = (l+2)(l-1), k = c(c+2)/6, omega = k/2 = (l-1)l(l+1)(l+2)/12,  r_s = 6/c
  chi(r) = sum_{n=1}^{l+2} c_n r^n,    psi0 = chi/(c r - 6),    g = psi0_x/psi0,    Vhat = g^2 - g_x = 2 g^2 - V_Zerilli
  Psi_DG = r (r+2)^k e^{-k r/2}/(c r - 6),    Psihat = (-d/dx + g) Psi_DG.

Closed form used for evaluation (no pole at r_s; derived and checked in check_modes.py, not assumed):
  Psihat = 2^k e^{-omega x} r^4 Pi(r)/chi(r),   Pi = Q/(c r - 6),   Q = [Nm]/r^5,
  Nm = sum_m (omega c_{m-2} + (m-2) c_{m-1} + 2 (m-1) c_m) r^m  (exact Fractions; Nm = O(r^5) and (c r - 6) | Q are asserted).
"""
import functools
import math
from fractions import Fraction as Fr

import mpmath as mp

mp.mp.dps = 40

MUT = set()   # mutations (check_*.py --mutate=...): 'chi' perturbs c_3, 'omega' shifts omega, 'jac' drops dx/dr in the norm, 'vhat' scales Vhat by 1 + 1e-6


def cc(l):
    return (l + 2) * (l - 1)


def kk(l):
    return (l - 1) * l * (l + 1) * (l + 2) // 6


def omega(l):
    """growth rate k/2 of Theorem 1 (units 1/|M|), exact."""
    w = Fr((l - 1) * l * (l + 1) * (l + 2), 12)
    return w + 1 if 'omega' in MUT else w


def rs(l):
    return Fr(6, cc(l))


def poch(a, m):
    p = 1
    for j in range(m):
        p *= a + j
    return p


@functools.lru_cache(None)
def chi_coeffs(l):
    """c_n, n = 0..l+2, exact (Lemma A).  c_1 = 3/2 is the n = 1 limit of the formula."""
    out = [Fr(0)] * (l + 3)
    out[1] = Fr(3, 2)
    for n in range(2, l + 3):
        out[n] = Fr((n - 2) * ((n - 4) * l * (l + 1) + n - 1) * poch(l - n + 3, 2 * n - 4),
                    2**n * math.factorial(n - 1)**2)
    if 'chi' in MUT and l >= 3:
        out[3] *= Fr(101, 100)
    return tuple(out)


@functools.lru_cache(None)
def pi_coeffs(l):
    """Coefficients of Pi (exact), with the two divisibility facts asserted."""
    c = chi_coeffs(l)
    w = omega(l)
    cn = lambda i: c[i] if 0 <= i <= l + 2 else Fr(0)
    nm = [w * cn(m - 2) + (m - 2) * cn(m - 1) + 2 * (m - 1) * cn(m) for m in range(l + 5)]
    assert all(v == 0 for v in nm[:5]), f"Nm is not O(r^5) at l={l}: {nm[:5]}"
    q = nm[5:]
    cl = cc(l)
    pi = []
    prev = Fr(0)
    for i in range(len(q) - 1):
        prev = (cl * prev - q[i]) / 6 if i else -q[0] / 6
        pi.append(prev)
    assert cl * pi[-1] == q[-1], f"(c r - 6) does not divide Q at l={l}"
    return tuple(pi)


def _mpl(cs):
    return [mp.mpf(a.numerator) / a.denominator for a in cs]


@functools.lru_cache(None)
def _mp_chi(l):
    return _mpl(chi_coeffs(l))


@functools.lru_cache(None)
def _mp_pi(l):
    return _mpl(pi_coeffs(l))


def horner(cs, r, nder=0):
    """value (and first nder derivatives) of sum cs[n] r^n at r (Horner, mp)."""
    out = []
    d = list(cs)
    for _ in range(nder + 1):
        v = mp.mpf(0)
        for a in reversed(d):
            v = v * r + a
        out.append(v)
        d = [i * d[i] for i in range(1, len(d))]
    return out


def chi(l, r, nder=2):
    return horner(_mp_chi(l), mp.mpf(r), nder)


def f_of(r):
    return 1 + 2 / r


def x_of(r):
    r = mp.mpf(r)
    return r - 2 * mp.log1p(r / 2)


def r_of_x(x):
    """Newton from the right (x(r) is convex and x <= r^2/4, so the root is >= 2 sqrt(x))."""
    x = mp.mpf(x)
    if x == 0:
        return mp.mpf(0)
    r = 2 * mp.sqrt(x) + x
    for _ in range(80):
        dr = (x_of(r) - x) * (r + 2) / r
        r -= dr
        if abs(dr) < mp.mpf(10)**(-mp.mp.dps + 4) * r:
            break
    return r


def g_and_gr(l, r):
    """g = psi0_x/psi0 and dg/dr, with psi0 = chi/(c r - 6)."""
    r = mp.mpf(r)
    c = cc(l)
    X, X1, X2 = chi(l, r, 2)
    a = X1 / X - c / (c * r - 6)
    da = X2 / X - (X1 / X)**2 + c**2 / (c * r - 6)**2
    f = f_of(r)
    fr = -2 / r**2
    return f * a, fr * a + f * da


def Vhat(l, r):
    """g^2 - g_x."""
    g, gr = g_and_gr(l, r)
    v = g * g - f_of(r) * gr
    return v * (1 + mp.mpf('1e-6')) if 'vhat' in MUT else v


def Vzer(l, r):
    """Static Zerilli potential, M = -1, lambda = c/2 (Gleiser-Dotti conventions; dg_proof Theorem 3)."""
    r = mp.mpf(r)
    lam = mp.mpf(cc(l)) / 2
    num = lam**2 * (lam + 1) * r**3 - 3 * lam**2 * r**2 + 9 * lam * r - 9
    return 2 * f_of(r) * num / (r**3 * (lam * r - 3)**2)


def Vhat2(l, r):
    """2 g^2 - V_Zerilli, an independent route to Vhat (needs psi0 to solve the static Zerilli equation)."""
    g, _ = g_and_gr(l, r)
    return 2 * g * g - Vzer(l, r)


def kappa_c(l, R):
    """d ln phi0/dx at R, phi0 = 1/psi0 = -g."""
    g, _ = g_and_gr(l, R)
    return -g


def psihat_unnorm(l, r):
    """e^{-omega x} r^4 Pi/chi (the constant 2^k dropped)."""
    r = mp.mpf(r)
    w = mp.mpf(omega(l).numerator) / omega(l).denominator
    pi = horner(_mp_pi(l), r, 0)[0]
    return mp.exp(-w * x_of(r)) * r**4 * pi / chi(l, r, 0)[0]


def psihat_direct(l, r):
    """(-d/dx + g) Psi_DG evaluated directly (pole of g and of Psi_DG cancel); valid away from r_s."""
    r = mp.mpf(r)
    k = kk(l)
    c = cc(l)
    psi = r * (r + 2)**k * mp.exp(-k * r / 2) / (c * r - 6)
    dlog = 1 / r + mp.mpf(k) / (r + 2) - mp.mpf(k) / 2 - c / (c * r - 6)
    g, _ = g_and_gr(l, r)
    return psi * (-f_of(r) * dlog + g)


def log_psihat(l, r):
    return mp.log(abs(psihat_unnorm(l, r)))


def w_of(l, r):
    """psihat_x/psihat from the closed form."""
    r = mp.mpf(r)
    w = mp.mpf(omega(l).numerator) / omega(l).denominator
    pi, pi1 = horner(_mp_pi(l), r, 1)
    X, X1 = chi(l, r, 1)
    return f_of(r) * (-w / f_of(r) + 4 / r + pi1 / pi - X1 / X)


def residual(l, r):
    """(H psihat + omega^2 psihat)/psihat = -psihat_xx/psihat + Vhat + omega^2 at r (should be 0)."""
    r = mp.mpf(r)
    w2 = (mp.mpf(omega(l).numerator) / omega(l).denominator)**2
    wx = w_of(l, r)
    wr = mp.diff(lambda s: w_of(l, s), r)
    pxx = f_of(r) * wr + wx**2
    return -pxx + Vhat(l, r) + w2


@functools.lru_cache(None)
def norm_psihat(l):
    """||psihat_unnorm||_{L^2(dx)}: integral of psihat^2 dx = psihat^2 r/(r+2) dr."""
    w = float(omega(l))
    rsc = 2 / math.sqrt(w)
    pts = [0] + [rsc * s for s in (0.5, 1, 2, 4, 8, 16, 32, 64)] + [mp.inf]
    F = (lambda r: psihat_unnorm(l, r)**2) if 'jac' in MUT else (lambda r: psihat_unnorm(l, r)**2 * r / (r + 2))
    return mp.sqrt(mp.quad(F, pts))


def psi_unit(l, r):
    return psihat_unnorm(l, r) / norm_psihat(l)

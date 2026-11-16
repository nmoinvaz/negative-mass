"""
Exact-rational machinery shared by threshold_isolation.py, large_l_proof.py and vhat_structure.py (#2).

Scaling: rho = c r, eps = 1/c, c = (l+2)(l-1) = L - 2 with L = l(l+1).
  X(rho) = c chi(rho/c) = sum_{n=1}^{l+2} a_n(eps) rho^n,   a_1 = 3/2, a_2 = 0 and, for n >= 3,
  a_n(eps) = (n-2) [(n-4) + 3 (n-3) eps] / (2^n (n-1)!^2) * prod_{a=0}^{n-3} (1 - (a-1)(a+2) eps).
The limit eps = 0 is X_inf = sum (n-2)(n-4)(rho/2)^n/(n-1)!^2.  The product vanishes at a = l, so the series stops at l+2.
(The formula is checked against common.chi_poly in every script that imports this module.)

  N(rho) = X - (rho-6) X'                      = c^2 chi^2 (d phi0/d r)/c  up to the positive factor c
  S(rho) = rho X [eps N + (2 + eps rho) N'] - (2 + eps rho) N (X + 2 rho X'),   N' = -(rho-6) X''
sign V^ = sign S on rho > 6 (derivation in RESULTS.md section 4).

All arithmetic is exact (fractions.Fraction).  A series is truncated at n1 = NTAIL terms and the omitted terms are
enclosed: for n >= 5, 0 <= a_n <= 3 n^2 / (2^n (n-1)!^2) when eps <= 1/4 (every product factor except a = 0 lies in
[0, 1], and the a = 0 factor is 1 + 2 eps <= 3/2), and the bound's consecutive-term ratio is checked <= 1/2.
"""
from fractions import Fraction as Fr
import math

NTAIL = 90          # terms kept exactly
MUT = set()         # mutation names switched on by a script's --mutate=...: coef, ntail, sterm, tail


def c_of(l):
    return (l + 2) * (l - 1)


def q(a):
    return (a - 1) * (a + 2)


def a_dict(l, nmax=None):
    """a_n for n = 1..min(l+2, nmax) as Fractions.  l = None gives the limit eps = 0 (needs nmax)."""
    if l is None:
        eps = Fr(0)
        top = nmax
    else:
        eps = Fr(1, c_of(l))
        top = l + 2 if nmax is None else min(l + 2, nmax)
    a = {1: Fr(3, 2), 2: Fr(0)}
    P = Fr(1)
    fact = Fr(1)       # (n-1)!
    for n in range(2, top + 1):
        fact *= (n - 1)
        if n >= 3:
            P *= 1 - q(n - 3) * eps
            a[n] = (n - 2) * ((n - 4) + 3 * (n - 3) * eps) / (2 ** n * fact ** 2) * P
    if 'coef' in MUT and 5 in a:
        a[5] *= 1 + Fr(1, 10 ** 8)          # mutation: perturb one coefficient
    return a, eps


def tail_bound(n1, j, rho_max):
    """sum_{n > n1} 3 n^(2+j) rho^(n-j) / (2^n (n-1)!^2) <= 2 * (first omitted term), ratio <= 1/2 asserted."""
    ratio = Fr(n1 + 1, n1) ** (2 + j) * rho_max / (2 * n1 ** 2)
    assert ratio <= Fr(1, 2), ("tail ratio", n1, j, rho_max, float(ratio))
    n = n1 + 1
    first = Fr(3 * n ** (2 + j)) * Fr(rho_max) ** (n - j) / (2 ** n * Fr(math.factorial(n - 1)) ** 2)
    return 2 * first


class Iv:
    """Closed interval with exact Fraction endpoints."""
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        self.lo = Fr(lo)
        self.hi = Fr(lo if hi is None else hi)
        assert self.lo <= self.hi

    @staticmethod
    def _c(o):
        return o if isinstance(o, Iv) else Iv(o)

    def __add__(self, o):
        o = Iv._c(o)
        return Iv(self.lo + o.lo, self.hi + o.hi)
    __radd__ = __add__

    def __neg__(self):
        return Iv(-self.hi, -self.lo)

    def __sub__(self, o):
        o = Iv._c(o)
        return Iv(self.lo - o.hi, self.hi - o.lo)

    def __rsub__(self, o):
        return Iv._c(o) - self

    def __mul__(self, o):
        o = Iv._c(o)
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return Iv(min(p), max(p))
    __rmul__ = __mul__

    def __pow__(self, k):
        out = Iv(1)
        for _ in range(k):
            out = out * self
        return out

    def __truediv__(self, k):
        k = Fr(k)
        return Iv(self.lo / k, self.hi / k) if k > 0 else Iv(self.hi / k, self.lo / k)

    def sign(self):
        if self.lo > 0:
            return 1
        if self.hi < 0:
            return -1
        return 0

    def width(self):
        return self.hi - self.lo

    def __repr__(self):
        return f"[{float(self.lo):.6g}, {float(self.hi):.6g}]"


def make_series(l, n1=None):
    """Coefficient data for X: dict a (n <= min(l+2, n1)) and whether a tail beyond n1 exists."""
    n1 = (3 if 'ntail' in MUT else NTAIL) if n1 is None else n1
    a, eps = a_dict(l, n1)
    has_tail = (l is None) or (l + 2 > n1)
    return dict(a=a, eps=eps, n1=n1, has_tail=has_tail, l=l)


def xder(ser, rho, jmax, rho_max):
    """Iv enclosures of X^(j)(rho), j = 0..jmax, for rho in (0, rho_max] given as Fraction."""
    rho = Fr(rho)
    out = []
    for j in range(jmax + 1):
        s = Fr(0)
        for n, an in ser['a'].items():
            if n < j or an == 0:
                continue
            s += an * math.perm(n, j) * rho ** (n - j)
        t = tail_bound(ser['n1'], j, Fr(rho_max)) if ser['has_tail'] else Fr(0)
        out.append(Iv(s, s + t))
    return out


def N_iv(ser, rho, rho_max=20):
    X = xder(ser, rho, 1, rho_max)
    return X[0] - (Fr(rho) - 6) * X[1]


def S_iv(ser, rho, rho_max=20):
    rho = Fr(rho)
    eps = ser['eps']
    X0, X1, X2 = xder(ser, rho, 2, rho_max)
    N = X0 - (rho - 6) * X1
    Np = -(rho - 6) * X2
    if 'sterm' in MUT:
        eps = Fr(0)                       # mutation: drop the eps-dependence of the weights in S
    return rho * X0 * (eps * N + (2 + eps * rho) * Np) - (2 + eps * rho) * N * (X0 + 2 * rho * X1)


def float_funcs(ser, nfloat=60):
    """double-precision N and S (not certified), used only to guess brackets."""
    a = {n: float(v) for n, v in ser['a'].items() if n <= nfloat}
    eps = float(ser['eps'])

    def X(rho, j):
        return sum(v * math.perm(n, j) * rho ** (n - j) for n, v in a.items() if n >= j)

    def N(rho):
        return X(rho, 0) - (rho - 6) * X(rho, 1)

    def S(rho):
        x0, x1, x2 = X(rho, 0), X(rho, 1), X(rho, 2)
        n = x0 - (rho - 6) * x1
        n1 = -(rho - 6) * x2
        return rho * x0 * (eps * n + (2 + eps * rho) * n1) - (2 + eps * rho) * n * (x0 + 2 * rho * x1)
    return N, S


def certify_root(ser, which, guess, lo_sign, rho_max=20, width0=Fr(1, 10 ** 13), maxdbl=40):
    """Return (lo, hi) Fractions with sign(F(lo)) = lo_sign and sign(F(hi)) = -lo_sign certified exactly, F = N or S.
    Starts at guess -/+ width0 and doubles the half-width until both signs certify."""
    F = N_iv if which == 'N' else S_iv
    guess = Fr(guess)
    w = width0
    for _ in range(maxdbl):
        lo, hi = guess - w, guess + w
        slo, shi = F(ser, lo, rho_max).sign(), F(ser, hi, rho_max).sign()
        if slo == lo_sign and shi == -lo_sign:
            return lo, hi
        w *= 2
    raise RuntimeError(f"could not certify a sign change near {float(guess)}")


def locate(ser):
    """Float guesses of c r_* and c r_neg (brentq on [7.9, 10.5] and [r_*+0.05, 12]); not certified."""
    from scipy.optimize import brentq
    Nf, Sf = float_funcs(ser)
    rn = brentq(Nf, 7.9, 10.5, xtol=1e-15, rtol=1e-15)
    rsg = brentq(Sf, rn + 0.05, 12.0, xtol=1e-15, rtol=1e-15)
    return rn, rsg


def frac_from_float(x):
    return Fr(x)


def refine(ser, which, lo, hi, lo_sign, steps, rho_max=20):
    """Exact bisection of a certified bracket (sign(F(lo)) = lo_sign, sign(F(hi)) = -lo_sign)."""
    F = N_iv if which == 'N' else S_iv
    for _ in range(steps):
        mid = (lo + hi) / 2
        sg = F(ser, mid, rho_max).sign()
        if sg == lo_sign:
            lo = mid
        elif sg == -lo_sign:
            hi = mid
        else:
            break            # the enclosure of F(mid) contains 0: stop refining
    return lo, hi

#!/usr/bin/env python3
"""Lemma A(i) of RESULTS.md section 1: the signs of c_n, for every l >= 2.

    c_n = (n-2)((n-4) l(l+1) + n-1) (l-n+3)_{2n-4} / (2^n ((n-1)!)^2).

Hand proof (stated here as lemmas for Lean, #135):
  A1-L1  c_1 = 3/2, c_2 = 0, c_3 = -K/32, c_4 = K/96, K = (l-1) l (l+1)(l+2) = L(L-2), L = l(l+1).
  A1-L2  For n >= 4 and l >= 2: (n-2) > 0 and (n-4) L + n - 1 > 0 (a sum of non-negative terms and n-1 > 0).
  A1-L3  (l-n+3)_{2n-4} is the product of the 2n-4 consecutive integers l-n+3, ..., l+n-2. Its largest factor l+n-2 > 0.
         If n <= l+2 the smallest factor l-n+3 >= 1, so the product is > 0. If n >= l+3 the factors run from <= 0 to > 0,
         so one factor is 0 and the product is 0.
  A1-L4  Hence c_n > 0 for 4 <= n <= l+2, c_n = 0 for n >= l+3, and c_{l+2} > 0 (since l+2 >= 4). chi has degree l+2.
  A1-L5  c_3 < 0 since K > 0 (l >= 2).
Checks below are the symbolic instances of L1 and L2 and exact data for L3, L4.
"""
import math
import sys

import sympy as sp

import review_lemma_a_common as C

l = sp.Symbol("l", positive=True)
a, b = sp.symbols("a b", nonnegative=True)

print(f"review_lemma_a_signs (mutation: {C.MUT})")

# A1.1 low orders, symbolic in l
K = (l - 1) * l * (l + 1) * (l + 2)
want = {1: sp.Rational(3, 2), 2: sp.Integer(0), 3: -K / 32, 4: K / 96}
for n, w in want.items():
    got = C.cn_sym(n, l)
    C.check(f"A1.1 c_{n} = {w}", sp.simplify(got - w) == 0, f"got {sp.factor(got)}")

# A1.2 factor positivity, symbolic: n = 4 + a, l = 2 + b
Lsym = (2 + b) * (3 + b)
for name, expr in (("(n-2)", 2 + a), ("(n-4)L+n-1", a * Lsym + (3 + a))):
    co = sp.Poly(sp.expand(expr), a, b).coeffs()
    C.check(f"A1.2 {name} > 0 for n >= 4, l >= 2: all {len(co)} coefficients positive", all(x > 0 for x in co) and len(co) > 0,
            f"coefficients {co}")
# the closed form of c_{l+2} (positive): (l)((l-2)L + l+1)(1)_{2l}... only the factors are asserted above and in A1.3

# A1.3 Pochhammer sign, exact data (literal shift 3): positive iff n <= l+2, zero iff n >= l+3
LMAX = 300


def poch_data(l_):
    for n in range(3, l_ + 7):
        v = math.prod(range(l_ - n + 3, l_ + n - 2 + 1))     # l-n+3 .. l+n-2, 2n-4 factors
        assert len(range(l_ - n + 3, l_ + n - 1)) == 2 * n - 4
        if n <= l_ + 2 and not v > 0:
            return False, f"n={n}: product {v} not > 0"
        if n >= l_ + 3 and v != 0:
            return False, f"n={n}: product {v} not 0"
    return True, ""


C.check_range("A1.3 (l-n+3)_{2n-4} > 0 for 3 <= n <= l+2, = 0 for n >= l+3", list(range(2, LMAX + 1)), poch_data)


# A1.4 signs of c_n, exact
def sign_data(l_):
    Kq = (l_ - 1) * l_ * (l_ + 1) * (l_ + 2)
    from fractions import Fraction
    for n in range(1, l_ + 7):
        v = C.cn_q(n, l_)
        if n == 1 and v != Fraction(3, 2):
            return False, f"c_1 = {v}"
        if n == 2 and v != 0:
            return False, f"c_2 = {v}"
        if n == 3 and v != Fraction(-Kq, 32):
            return False, f"c_3 = {v}, want {Fraction(-Kq, 32)}"
        if 4 <= n <= l_ + 2 and not v > 0:
            return False, f"c_{n} = {float(v):.3e} not > 0"
        if n >= l_ + 3 and v != 0:
            return False, f"c_{n} = {v} not 0"
    if not C.cn_q(l_ + 2, l_) > 0:
        return False, "c_{l+2} not > 0"
    return True, ""


C.check_range("A1.4 signs of c_n (c_1 = 3/2, c_2 = 0, c_3 = -K/32 < 0, c_n > 0 for 4..l+2, 0 beyond, c_{l+2} > 0)",
              list(range(2, LMAX + 1)), sign_data)

C.finish()

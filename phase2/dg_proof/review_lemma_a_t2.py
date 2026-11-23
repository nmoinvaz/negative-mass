#!/usr/bin/env python3
"""T2 of RESULTS.md section 3 (Corollary): N has exactly one root r_* on (r_s, inf), for every l >= 2.

Hand proof (stated here as lemmas for Lean, #135). t = r - r_s, c r - 6 = c t, d_m = chi^{(m)}(r_s)/m!.
  T2-L1  chi(r_s + t) = sum_m d_m t^m and t chi' = sum_m m d_m t^m, so N(r_s + t) = c chi - c t chi' = c sum_m (1-m) d_m t^m (finite, degree <= l+2).
  T2-L2  d_m = sum_{n >= m} c_n C(n, m) r_s^{n-m} (Taylor shift of the polynomial chi).
  T2-L3  For 4 <= m <= l+2 every term has n >= 4, so c_n >= 0, C(n,m) > 0, r_s > 0; the term n = l+2 is present with c_{l+2} > 0. So d_m > 0.
         d_m = 0 for m > l+2.
  T2-L4  d_0 = chi(r_s) = N(r_s)/c > 0 by T1 (N(r_s) >= 9/4), d_1 term has factor 1-m = 0, d_2 = 0 by L1, the m = 3 term is -2c d_3 of unknown sign.
  T2-L5  The coefficient list is (+, 0, 0, -2c d_3, then strictly negative terms up to m = l+2): one sign change whatever the sign of d_3
         (the list has at least one negative term since l+2 >= 4). Descartes: the number of roots in t > 0 (with multiplicity) is <= 1 and
         has the parity of the number of sign changes; N(r_s) > 0 and N -> -inf (leading coefficient -c(l+1)c_{l+2} < 0) give at least one. So exactly one,
         simple. N has no root on [0, r_s] by T1. So exactly one root on (r_s, inf).
"""
import sympy as sp
from fractions import Fraction
from math import comb

import review_lemma_a_common as C

print(f"review_lemma_a_t2 (mutation: {C.MUT})")
C.run_planted({"PL5"})
C.gate_planted()

r = C.R
c = sp.Symbol("c")
t = sp.Symbol("t")


def shiftfac(m):
    return (2 if C.MUT == "t2_shift" else 1) - m


# ---- T2.1 generic Taylor identity (degree 8, symbolic coefficients)
a = {k: sp.Symbol(f"a{k}") for k in range(1, 9)}
chi = sum(a[k] * r**k for k in a)
rs = 6 / c
Nr = c * chi - (c * r - 6) * sp.diff(chi, r)
lhs = sp.expand(Nr.subs(r, rs + t))
dm = {m: sp.diff(chi, r, m).subs(r, rs) / sp.factorial(m) for m in range(0, 9)}
rhs = sp.expand(c * sum(shiftfac(m) * dm[m] * t**m for m in range(0, 9)))
C.check("T2.1 N(r_s + t) = c sum_m (1-m) d_m t^m, d_m = chi^(m)(r_s)/m!, generic chi of degree 8",
        sp.simplify(lhs - rhs) == 0, f"difference {sp.factor(sp.simplify(lhs - rhs))}")

# ---- T2.2 d_m = sum_{n>=m} c_n C(n,m) r_s^{n-m}
okgen = all(sp.simplify(dm[m] - sum(a[k] * sp.binomial(k, m) * rs ** (k - m) for k in range(max(m, 1), 9))) == 0 for m in range(0, 9))


def taylor(l):
    cl = C.c_of_l(l)
    rsq = sp.Rational(6, cl)
    chip = C.chi_poly(l)
    sh = sp.Poly(chip.as_expr().subs(r, rsq + t), t, domain="QQ")
    sf = [sum(sp.Rational(C.cn_q(nn, l).numerator, C.cn_q(nn, l).denominator) * sp.binomial(nn, m) * rsq ** (nn - m) for nn in range(max(m, 1), l + 3))
          for m in range(0, l + 3)]
    got = [sh.coeff_monomial(t**m) for m in range(0, l + 3)]
    return got == sf, "Taylor shift differs from the binomial sum"


ok2 = C.check("T2.2a d_m = sum_{n>=m} c_n C(n,m) r_s^{n-m} for generic chi (degree 8)", okgen, "mismatch")
C.check_range("T2.2b the same for the actual chi (sympy Taylor shift against the sum), l = 2..40", list(range(2, 41)), taylor)

# ---- exact d_m for l = 2..120 by the binomial sum
LMAX = 120
D = {}


def dm_exact(l):
    cl = C.c_of_l(l)
    rsq = Fraction(6, cl)
    cs = {nn: C.cn_q(nn, l) for nn in range(1, l + 3)}
    pw = [rsq**k for k in range(0, l + 3)]
    return [sum(cs[nn] * comb(nn, m) * pw[nn - m] for nn in range(max(m, 1), l + 3)) for m in range(0, l + 3)]


for l in range(2, LMAX + 1):
    D[l] = dm_exact(l)


# ---- T2.3
def t23(l):
    d = D[l]
    rsq = Fraction(6, C.c_of_l(l))
    if not d[0] > 0:
        return False, f"d_0 = {d[0]}"
    if d[1] != -d[0] / (rsq * (rsq + 2)):
        return False, f"d_1 = {d[1]}, want {-d[0] / (rsq * (rsq + 2))}"
    if d[2] != 0:
        return False, f"d_2 = {d[2]}"
    return True, ""


C.check_range("T2.3 d_0 = chi(r_s) > 0, d_1 = -d_0/(r_s(r_s+2)), d_2 = 0 (exact)", list(range(2, LMAX + 1)), t23)


# ---- T2.4
def t24(l):
    d = D[l]
    for m in range(4, l + 3):
        if not d[m] > 0:
            return False, f"d_{m} = {float(d[m]):.3e} not > 0"
    return True, ""


C.check_range("T2.4 d_m > 0 for 4 <= m <= l+2 (exact; d_m = 0 beyond since chi has degree l+2)", list(range(2, LMAX + 1)), t24)
sg = {l: (D[l][3] > 0) - (D[l][3] < 0) for l in range(2, LMAX + 1)}
pos = [l for l in sg if sg[l] > 0]
neg = [l for l in sg if sg[l] < 0]
zer = [l for l in sg if sg[l] == 0]
print(f"info: sign of d_3 for l = 2..{LMAX}: positive at {len(pos)} values of l ({pos[:6]}...{pos[-2:]}), negative at {len(neg)} ({neg[:6]}...{neg[-2:]}), zero at {zer}")


# ---- T2.5
def t25(l):
    cl = C.c_of_l(l)
    d = list(D[l])
    if C.MUT == "t2_flip_d0":
        d[0] = -d[0]
    e = [cl * shiftfac(m) * d[m] for m in range(0, l + 3)]
    sc = C.sign_changes(e)
    lead = e[l + 2]
    if C.MUT != "t2_shift" and lead != -cl * (l + 1) * C.cn_q(l + 2, l):
        return False, f"leading coefficient {lead}, want {-cl * (l + 1) * C.cn_q(l + 2, l)}"
    return sc == 1 and lead < 0, f"sign changes {sc}, leading coefficient sign {(lead > 0) - (lead < 0)}"


C.check_range("T2.5 coefficients of N(r_s + t) have exactly one sign change; leading coefficient -c(l+1)c_{l+2} < 0",
              list(range(2, LMAX + 1)), t25)


# ---- T2.6 exact Sturm count
def t26(l):
    cl = C.c_of_l(l)
    chip = C.chi_poly(l).as_expr()
    Np = sp.Poly(cl * chip - (cl * r - 6) * sp.diff(chip, r), r, domain="QQ")
    rsq = sp.Rational(6, cl)
    upper = rsq * sp.Rational(6, 5) if C.MUT == "t2_window" else sp.oo
    out = Np.count_roots(rsq, upper)
    inside = Np.count_roots(0, rsq)
    return out == 1 and inside == 0 and Np.is_sqf, f"roots in (r_s, upper] {out}, in [0, r_s] {inside}, squarefree {Np.is_sqf}"


C.check_range("T2.6 exact Sturm: one root of N in (r_s, inf), none in [0, r_s], N squarefree (root is simple)", list(range(2, 41)), t26)
C.finish()

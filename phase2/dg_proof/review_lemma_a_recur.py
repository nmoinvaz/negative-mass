#!/usr/bin/env python3
"""Lemma A(ii) of RESULTS.md section 1: psi_0 = chi/(cr-6) solves the static Zerilli equation for every l.

Hand derivation (stated here as lemmas for Lean, #135). d/dx = f d/dr, f = (r+2)/r, D = cr-6, psi = chi/D.
  A2-L1  psi_xx = f (f psi')' = f^2 psi'' + f f' psi'. With psi' = chi'/D - c chi/D^2, psi'' = chi''/D - 2c chi'/D^2 + 2c^2 chi/D^3
         and V of RESULTS 4.2: (psi_xx - V psi) * (-r^4 D^2/(r+2)) = P2 chi'' + P1 chi' + P0 chi with
         P2 = -r^2(r+2)D, P1 = 2cr^3 + 6cr^2 - 12r, P0 = c^2r^2 - 6cr + 12 (uses 8[...] = c^2(c+2)r^3 - 6c^2r^2 + 36cr - 72 and
         -c^3r^3 + 12c^2r^2 - 48cr + 72 = -(cr-6)(c^2r^2 - 6cr + 12)).
  A2-L2  On r^n the operator gives (c(c - n(n-3))) r^{n+2} + (n-1)(6n - 2c(n-3)) r^{n+1} + 12(n-1)^2 r^n, hence for chi = sum c_m r^m
         the coefficient of r^m is 12(m-1)^2 c_m + (m-2)(6(m-1) - 2c(m-4)) c_{m-1} + c(c-(m-2)(m-5)) c_{m-2}, and
         c - (m-2)(m-5) = (l-m+4)(l+m-3) for c = (l+2)(l-1).
  A2-L3  For m >= 4: (l-m+3)_{2m-4} = (l-m+3)(l-m+4)(l+m-3)(l+m-2) Pi and (l-m+4)_{2m-6} = (l-m+4)(l+m-3) Pi, with Pi = (l-m+5)_{2m-8};
         Pi > 0 for m <= l+4. Dividing the recurrence by Pi (and by 2^m (m-3)!^2) leaves a polynomial in (l, m): it is 0.
  A2-L4  m = 1, 2, 3, 4 are checked directly from c_1..c_4 (c_0 = c_{-1} = 0). For m >= l+5, c_m = c_{m-1} = c_{m-2} = 0.
  A2-L5  So the coefficient of every power of r in P2 chi'' + P1 chi' + P0 chi vanishes: chi solves the ODE, psi_0 solves psi_xx = V psi.
"""
import sys

import sympy as sp

import review_lemma_a_common as C

print(f"review_lemma_a_recur (mutation: {C.MUT})")
C.run_planted({"PL1", "PL2"})
C.gate_planted()

r = C.R
c = sp.Symbol("c")
n = sp.Symbol("n")
l, m = sp.symbols("l m")

# ---- A2.1: P2, P1, P0 from V, symbolic in c, chi a generic function
chi = sp.Function("chi")
D = c * r - 6
psi = chi(r) / D
V = C.V_zerilli(r, c)
resid = C.zerilli_residual(psi, V, r)
x0_, x1_, x2_ = sp.symbols("x0 x1 x2")
expr = resid.subs(sp.Derivative(chi(r), (r, 2)), x2_).subs(sp.Derivative(chi(r), r), x1_).subs(chi(r), x0_)
pref = -(r**4) * D**2 / (r + 2)
derived = [sp.cancel(pref * sp.diff(expr, xk)) for xk in (x2_, x1_, x0_)]
hand = C.P_hand(r, c)
ok1 = all(sp.cancel(d - h) == 0 for d, h in zip(derived, hand))
C.check("A2.1 P2, P1, P0 derived from V equal the hand forms (and psi_xx - V psi = -(r+2)/(r^4 D^2) (P2 chi'' + P1 chi' + P0 chi))",
        ok1, f"derived {[sp.factor(d) for d in derived]}")
ispoly = all(sp.cancel(d).is_polynomial(r) for d in derived)
Pd = derived if ispoly else list(hand)
if not ispoly:
    print("note: the derived P are not polynomials; the following checks use the hand forms")

# ---- A2.2: P2 as RESULTS states it
C.check("A2.2 P2 = -r^2 (r+2)(cr-6) as RESULTS states", sp.expand(C.P2_results(r, c) - Pd[0]) == 0,
        f"RESULTS' P2 minus derived: {sp.factor(sp.expand(C.P2_results(r, c) - Pd[0]))}")

# ---- A2.3: recurrence on r^n extracted from the ODE, symbolic in n
B = C.three_term(Pd[0], Pd[1], Pd[2], r)
shifts_ok = sorted(B) == [0, 1, 2]
Bm = lambda s, idx: sp.sympify(B.get(s, 0)).subs(n, idx)  # noqa: E731
b0, b1, b2 = Bm(0, m), Bm(1, m - 1), Bm(2, m - 2)
h0, h1, h2 = C.rec_hand(m, c)
ok3a = shifts_ok and all(sp.expand(x - y) == 0 for x, y in ((b0, h0), (b1, h1), (b2, h2)))
C.check("A2.3a extracted recurrence equals my typed formula", ok3a, f"shifts {sorted(B)}; diffs "
        f"{[sp.factor(sp.expand(x - y)) for x, y in ((b0, h0), (b1, h1), (b2, h2))]}")
# RESULTS' form: 12(m-1)^2 c_m - 2(m-2)(...) c_{m-1} + c (l-m+4)(l+m-3) c_{m-2} = 0; the bracket (...) is elided in RESULTS.
# Mine is (...) = c(m-4) - 3(m-1). The c_{m-2} coefficient is compared after c = (l+2)(l-1).
cl = (l + 2) * (l - 1)
res0 = 12 * (m - 1) ** 2
res1 = -2 * (m - 2) * (c * (m - 4) - 3 * (m - 1))
res2 = c * (l - m + 4) * (l + m - 3)
ok3b = (sp.expand(b0 - res0) == 0 and sp.expand(b1 - res1) == 0 and sp.expand(b2.subs(c, cl) - res2.subs(c, cl)) == 0)
C.check("A2.3b extracted recurrence equals RESULTS' form (c_m, c_{m-1} bracket as derived, c_{m-2} with c = (l+2)(l-1))", ok3b,
        f"diffs {[sp.factor(sp.expand(x.subs(c, cl) - y.subs(c, cl))) for x, y in ((b0, res0), (b1, res1), (b2, res2))]}")

# ---- A2.4: low orders m = 1..4, symbolic in l (c_0 = c_{-1} = 0)
cs = {k: C.cn_sym(k, l) for k in range(1, 5)}
cs[0] = sp.Integer(0)
cs[-1] = sp.Integer(0)
for mm in range(1, 5):
    res = (b0.subs(m, mm) * cs[mm] + b1.subs(m, mm) * cs[mm - 1] + b2.subs(m, mm) * cs[mm - 2]).subs(c, cl)
    C.check(f"A2.4 m={mm} low-order residual symbolic in l", sp.simplify(res) == 0, f"residual {sp.factor(sp.simplify(res))}")

# ---- A2.5: polynomial identity in (l, m) after dividing by Pi = (l-m+5)_{2m-8}
Lsym = l * (l + 1)
F = sp.Symbol("F")          # (m-3)!
fact = {0: (m - 1) * (m - 2) * F, 1: (m - 2) * F, 2: F}       # (m-1)!, (m-2)!, (m-3)!
Rm = (l - m + 3) * (l - m + 4) * (l + m - 3) * (1 if C.MUT == "poch_ratio" else (l + m - 2))
Rm1 = (l - m + 4) * (l + m - 3)
R = {0: Rm, 1: Rm1, 2: sp.Integer(1)}
terms = 0
for d, bk in ((0, b0), (1, b1), (2, b2)):
    nn = m - d
    A = (nn - 2) * ((nn - 4) * Lsym + nn - 1) * R[d] * sp.Integer(2) ** d / (fact[d] ** 2)       # 2^m c_{nn}/Pi, 2^d for 2^m/2^nn
    terms += bk * A
poly_id = sp.cancel(sp.together((terms * F**2).subs(c, cl)))
C.check("A2.5 polynomial identity in (l, m) after dividing by (l-m+5)_{2m-8}", sp.expand(poly_id) == 0,
        f"nonzero polynomial with {len(sp.Poly(sp.expand(sp.numer(poly_id)), l, m).terms())} terms")

# ---- A2.6: the Pochhammer reductions of A2.5 on a grid, Pi > 0, vanishing for m >= l+5
def grid(l_):
    for mm in range(4, l_ + 5):
        Pi = C.poch_q(l_ - mm + 5, 2 * mm - 8)
        if not Pi > 0:
            return False, f"Pi <= 0 at m={mm}"
        r0 = (l_ - mm + 3) * (l_ - mm + 4) * (l_ + mm - 3) * (1 if C.MUT == "poch_ratio" else (l_ + mm - 2)) * Pi
        r1 = (l_ - mm + 4) * (l_ + mm - 3) * Pi
        if C.poch_q(l_ - mm + 3, 2 * mm - 4) != r0:
            return False, f"(l-m+3)_(2m-4) reduction fails at m={mm}"
        if C.poch_q(l_ - mm + 4, 2 * mm - 6) != r1:
            return False, f"(l-m+4)_(2m-6) reduction fails at m={mm}"
    for mm in range(l_ + 5, l_ + 9):
        if any(C.cn_q(k, l_) != 0 for k in (mm, mm - 1, mm - 2)):
            return False, f"c_k nonzero for k near m={mm}"
    return True, ""


C.check_range("A2.6 Pochhammer reductions exact, Pi > 0 for 4 <= m <= l+4, c_m = c_{m-1} = c_{m-2} = 0 for m >= l+5", list(range(2, 41)), grid)

# ---- A2.7: the recurrence with the actual c_n, exact
fb = [sp.lambdify((n, c), sp.sympify(B.get(s, 0)), modules="math") for s in (0, 1, 2)]


def rec_data(l_):
    cl_ = C.c_of_l(l_)
    cc = lambda k: C.cn_q(k, l_) if k >= 1 else 0      # noqa: E731
    for mm in range(0, l_ + 7):
        v = fb[0](mm, cl_) * cc(mm) + (fb[1](mm - 1, cl_) * cc(mm - 1) if mm - 1 >= 1 else 0) + \
            (fb[2](mm - 2, cl_) * cc(mm - 2) if mm - 2 >= 1 else 0)
        if v != 0:
            return False, f"residual {v} at m={mm}"
    return True, ""


C.check_range("A2.7 recurrence residual exactly 0 for m = 0..l+6", list(range(2, 61)), rec_data)

# ---- A2.8: P2 chi'' + P1 chi' + P0 chi is the zero polynomial
def poly_data(l_):
    cl_ = C.c_of_l(l_)
    chip = C.chi_poly(l_)
    P = [sp.Poly(sp.expand(p.subs(c, cl_)), C.R, domain="QQ") for p in Pd]
    tot = P[0] * chip.diff(C.R).diff(C.R) + P[1] * chip.diff(C.R) + P[2] * chip
    return tot.is_zero, f"nonzero polynomial of degree {tot.degree()}"


C.check_range("A2.8 P2 chi'' + P1 chi' + P0 chi = 0 as a polynomial", list(range(2, 41)), poly_data)

C.finish()

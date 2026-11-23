#!/usr/bin/env python3
"""T1 of RESULTS.md section 3 (Corollary): N = c chi - (cr-6) chi' > 0 on [0, r_s], r_s = 6/c, for every l >= 2.

Hand proof (stated here as lemmas for Lean, #135):
  T1-L1  N = sum_n c_n [c(1-n) r^n + 6n r^{n-1}] (term by term: c c_n r^n - (cr-6) n c_n r^{n-1}).
         n = 1: 6 c_1 = 9. n = 2: 0. n = 3: 18 c_3 r^2 - 2c c_3 r^3. n = 4: 24 c_4 r^3 - 3c c_4 r^4.
         n >= 5: (n-1) c_n r^{n-1}(6 - cr) + 6 c_n r^{n-1} (from -(n-1) c r^n + 6n r^{n-1}).
         n = 4 splits as 6 c_4 r^3 + 3 c_4 r^3 (6 - cr).
         Hence N - [9 + 18 c_3 r^2 + (6 c_4 - 2c c_3) r^3] = sum_{m>=4} (m-1) c_m r^{m-1}(6 - cr) + sum_{j>=5} 6 c_j r^{j-1}.
  T1-L2  On [0, r_s]: 6 - cr >= 0 and c_m >= 0 for m >= 4 (Lemma A(i)), so the right side is >= 0.
  T1-L3  With c_3 = -K/32, c_4 = K/96, K = c(c+2), rho = cr, e = 1/c: the bracket is B(rho) = 9 - (9/16)(1+2e) rho^2 + (1/16)(1+2e)(1+e) rho^3.
         B'(rho) = (3/16)(1+2e) rho ((1+e) rho - 6): negative on (0, 6/(1+e)), positive after. rho* = 6/(1+e) < 6, so on [0, 6] the minimum is
         B(rho*) = 9 - (27/4)(1+2e)/(1+e)^2 = 9/4 + (27/4) e^2/(1+e)^2 >= 9 - 27/4 = 9/4 > 0.
  T1-L4  So N >= B >= 9/4 on [0, r_s] (including r_s, where N = c chi(r_s)); in particular chi(r_s) > 0.
"""
import sympy as sp

import review_lemma_a_common as C

print(f"review_lemma_a_t1 (mutation: {C.MUT})")
C.run_planted({"PL6"})
C.gate_planted()

r = C.R
c = sp.Symbol("c")
n = sp.Symbol("n")
lsym = sp.Symbol("l")
rho, eps = sp.symbols("rho epsilon", positive=True)


def bracket(C3, C4, cc, rr):
    return 9 + (17 if C.MUT == "t1_coef18" else 18) * C3 * rr**2 + (6 * C4 - 2 * cc * C3) * rr**3


def rhs_terms(Cm, m, cc, rr):
    """(m-1) C_m r^{m-1}(6 - cr) [m >= 4] and 6 C_j r^{j-1} [j >= 5]."""
    t = (m - 1) * Cm * rr ** (m - 1) * (6 - cc * rr) if m >= 4 else 0
    if m >= 5 and C.MUT != "t1_drop_sum2":
        t += 6 * Cm * rr ** (m - 1)
    return t


# ---- T1.1 identity for generic coefficients C_1 = 3/2, C_2 = 0, C_3..C_9 symbols, and termwise for symbolic n
Cs = {k: sp.Symbol(f"C{k}") for k in range(3, 10)}
coef = {1: sp.Rational(3, 2), 2: sp.Integer(0), **Cs}
chi = sum(coef[k] * r**k for k in range(1, 10))
N = c * chi - (c * r - 6) * sp.diff(chi, r)
diff = sp.expand(N - bracket(coef[3], coef[4], c, r) - sum(rhs_terms(coef[m], m, c, r) for m in range(4, 10)))
Cn = sp.Symbol("Cn")
lhs_n = Cn * (c * (1 - n) * r + 6 * n)                       # term n >= 5 divided by r^{n-1}
rhs_n = (n - 1) * Cn * (6 - c * r) + (0 if C.MUT == "t1_drop_sum2" else 6 * Cn)
diff_n = sp.expand(lhs_n - rhs_n)
C.check("T1.1 identity for generic coefficients (degree 9) and termwise for symbolic n >= 5",
        diff == 0 and diff_n == 0, f"degree-9 difference {sp.factor(diff)}; termwise {sp.factor(diff_n)}")


# ---- T1.2 identity with the actual c_n, exact polynomial identity in r, l = 2..40
def ident(l):
    cl = C.c_of_l(l)
    cq = {k: C.cn_q(k, l) for k in range(1, l + 5)}
    cs = {k: sp.Rational(v.numerator, v.denominator) for k, v in cq.items()}
    chip = sum(cs[k] * r**k for k in cs)
    Np = sp.expand(cl * chip - (cl * r - 6) * sp.diff(chip, r))
    br = bracket(cs[3], cs[4], cl, r)
    rhs = sum(rhs_terms(cs[m], m, cl, r) for m in range(4, l + 5))
    d = sp.expand(Np - br - rhs)
    return d == 0, f"N - bracket - RHS = {sp.factor(d)}"


C.check_range("T1.2 N - bracket = sum_{m>=4}(m-1)c_m r^{m-1}(6-cr) + sum_{j>=5} 6 c_j r^{j-1} with the actual c_n (exact)", list(range(2, 41)), ident)

# ---- T1.3 the cubic: closed form in (rho, e), derivative, minimum
cl_sym = (lsym + 2) * (lsym - 1)
br_l = bracket(C.cn_sym(3, lsym), C.cn_sym(4, lsym), cl_sym, rho / cl_sym)
Bform = 9 - sp.Rational(9, 16) * (1 + 2 * eps) * rho**2 + sp.Rational(1, 16) * (1 + 2 * eps) * (1 + eps) * rho**3
okform = sp.simplify(br_l.subs(lsym, lsym) - Bform.subs(eps, 1 / cl_sym)) == 0
C.check("T1.3a bracket with closed-form c_3, c_4 equals 9 - (9/16)(1+2e)rho^2 + (1/16)(1+2e)(1+e)rho^3", okform,
        f"difference {sp.factor(sp.simplify(br_l - Bform.subs(eps, 1 / cl_sym)))}")
dB = sp.diff(Bform, rho)
C.check("T1.3b B'(rho) = (3/16)(1+2e) rho ((1+e) rho - 6)",
        sp.simplify(dB - sp.Rational(3, 16) * (1 + 2 * eps) * rho * ((1 + eps) * rho - 6)) == 0, f"{sp.factor(dB)}")
rstar = 6 / (1 + eps)
minval = sp.simplify(Bform.subs(rho, rstar))
C.check("T1.3c B(6/(1+e)) = 9 - (27/4)(1+2e)/(1+e)^2 = 9/4 + (27/4) e^2/(1+e)^2",
        sp.simplify(minval - (9 - sp.Rational(27, 4) * (1 + 2 * eps) / (1 + eps) ** 2)) == 0
        and sp.simplify(minval - (sp.Rational(9, 4) + sp.Rational(27, 4) * eps**2 / (1 + eps) ** 2)) == 0, f"{sp.factor(minval)}")
claimed = 3 if C.MUT == "t1_bound3" else 9 - sp.Rational(27, 4)
numer = sp.numer(sp.together(minval - claimed))
co = sp.Poly(sp.expand(numer), eps).coeffs()
C.check("T1.3d B(rho*) - (9 - 27/4) has a numerator with non-negative coefficients in e (so >= 0 for every e > 0)",
        all(x >= 0 for x in co), f"minimum minus claimed bound = {sp.factor(sp.simplify(minval - claimed))}; numerator coefficients {co}")
C.check("T1.3e rho* = 6/(1+e) lies in (0, 6) for e > 0 (so the minimum on [0, 6] is at rho*)",
        sp.simplify(6 - rstar - 6 * eps / (1 + eps)) == 0 and sp.simplify(6 * eps / (1 + eps)).is_positive is True, "")


def cubic_min(l):
    cl = C.c_of_l(l)
    e = sp.Rational(1, cl)
    p = Bform.subs(eps, e)
    v, arg = C.poly_min(p, rho, 0, 6)
    ok = v == minval.subs(eps, e) and v >= sp.Rational(9, 4)
    return ok, f"min {v} at {arg}, closed form {minval.subs(eps, e)}"


C.check_range("T1.3f exact minimum on [0, 6] of B for e = 1/c equals the closed form and is >= 9/4", list(range(2, 41)), cubic_min)


# ---- T1.4 N has no root on [0, r_s], N(0) = 9, N(r_s) = c chi(r_s) >= 9/4
def roots_in(l):
    cl = C.c_of_l(l)
    chip = C.chi_poly(l)
    Np = sp.Poly(cl * chip.as_expr() - (cl * r - 6) * sp.diff(chip.as_expr(), r), r, domain="QQ")
    rsq = sp.Rational(6, cl)
    cnt = Np.count_roots(0, rsq)
    n0, ns = Np.eval(0), Np.eval(rsq)
    ok = cnt == 0 and n0 == 9 and ns == cl * chip.eval(rsq) and ns >= sp.Rational(9, 4)
    return ok, f"roots in [0,r_s]: {cnt}, N(0) = {n0}, N(r_s) = {ns}"


C.check_range("T1.4 N has no root in [0, r_s] (exact Sturm), N(0) = 9, N(r_s) = c chi(r_s) >= 9/4", list(range(2, 41)), roots_in)
C.finish()

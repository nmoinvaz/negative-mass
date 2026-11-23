"""#142: Lemma A(iii) of RESULTS.md (chi > 0 on r > 0), written from RESULTS.md alone.

Run: python review_lemma_a3_a3.py [--mutate=NAME] > review_lemma_a3_a3_output.txt
Planted PLA (critical point of a quartic with known answer) and PLS (Sturm counter) run first.
Mutations: pl_crit_root, pl_sturm_sign, a3_c5, a3_s, a3_qclosed, a3_ident, a3_qsign, a3_sturm_neg, a3_l2_min,
a3_K120, cn_n4to3.
"""
import sys
import time
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from review_lemma_a3_common import *

run = Run('review_lemma_a3_a3: Lemma A(iii), chi > 0')
T0 = time.time()
s, K, rho = symbols('s K rho', positive=True)
a_, b_, d_ = symbols('a b d', positive=True)


# ------------------------------------------------------------------ planted problems
def crit_min(c0, a, b, d):
    """q = c0 - a r^2 + b r^3 + d r^4. q' = r (4 d r^2 + 3 b r - 2 a): returns (r0, q(r0)) for the positive root."""
    disc = 9 * b ** 2 + 32 * a * d
    sgn = -1 if mut('pl_crit_root') else 1   # mutation pl_crit_root: take the other (negative) root
    r0 = (-3 * b + sgn * sp.sqrt(disc)) / (8 * d)
    return r0, c0 - a * r0 ** 2 + b * r0 ** 3 + d * r0 ** 4


r0, q0 = crit_min(Q(3, 2), Q(3), Q(1), Q(3, 4))
run.planted('PLA.1', r0 == 1 and q0 == Q(1, 4), 'got r0=%s q(r0)=%s, want 1, 1/4' % (r0, q0))
r0, q0 = crit_min(Q(3, 2), Q(1), Q(0), Q(1, 8))
run.planted('PLA.2', r0 == 2 and q0 == Q(-1, 2) and not (q0 > 0), 'got r0=%s q(r0)=%s, want 2, -1/2' % (r0, q0))
r0s, _ = crit_min(1, a_, b_, d_)
res = sp.simplify(4 * d_ * r0s ** 2 + 3 * b_ * r0s - 2 * a_)
run.planted('PLA.3', res == 0, 'quadratic residual %s' % res)

planted_sturm(run)

# ------------------------------------------------------------------ A3.1 c3, c4, c5 closed forms (symbolic in l)
Kl = K_of(l)
c5_typed = Kl * (Kl - 24) / (6000 if mut('a3_c5') else 6144)   # mutation a3_c5: wrong denominator
ok1 = cancel(cn(1, l) - Q(3, 2)) == 0 and cn(2, l) == 0
ok3 = cancel(cn(3, l) + Kl / 32) == 0
ok4 = cancel(cn(4, l) - Kl / 96) == 0
ok5 = cancel(cn(5, l) - c5_typed) == 0
run.check('A3.1', ok1 and ok3 and ok4 and ok5,
          'c1,c2 %s; c3 %s; c4 %s; c5 %s' % (ok1, ok3, ok4, ok5))

# ------------------------------------------------------------------ A3.2 K facts
kmin = 121 if mut('a3_K120') else 120   # mutation a3_K120: claim K >= 121 for l >= 3
facts = []
facts.append(('K(2)=24', K_of(2) == 24))
facts.append(('K(l)>=%d for 3<=l<=50' % kmin, all(K_of(i) >= kmin for i in range(3, 51))))
facts.append(('K-24=(l-2)(l+3)(l^2+l+4)', sp.expand(Kl - 24 - (l - 2) * (l + 3) * (l ** 2 + l + 4)) == 0))
facts.append(('c5(2)=0', cn(5, 2) == 0))
facts.append(('c5(l)>0, 3<=l<=50', all(cn(5, i) > 0 for i in range(3, 51))))
facts.append(('K strictly increasing, K>24 for l>=3', all(K_of(i + 1) > K_of(i) for i in range(2, 50))))
run.check('A3.2', all(v for _, v in facts), '; '.join(n for n, v in facts if not v))

# ------------------------------------------------------------------ A3.3 chi/r - q coefficients
bad = []
for lv in range(3, 51):
    P = chi_over_r(lv)
    coeffs = [P.coeff_monomial(r ** p) for p in range(lv + 2)]
    Kv = K_of(lv)
    sg = -1 if mut('a3_qsign') else 1   # mutation a3_qsign: flip the sign of the cubic term of q
    qc = [Q(3, 2), 0, -Q(Kv, 32), sg * Q(Kv, 96), cn(5, lv)]
    for p in range(5):
        if coeffs[p] != qc[p]:
            bad.append((lv, p))
    for p in range(5, lv + 2):
        if coeffs[p] < 0:
            bad.append((lv, p))
    if bad:
        break
    run.check('A3.3[l=%d]' % lv, True)
if bad:
    run.check('A3.3[l=%d]' % bad[0][0], False, 'power %d of chi/r - q wrong' % bad[0][1])

# ------------------------------------------------------------------ A3.4 critical point of q (symbolic in s, K = s^2/6 + 18)
Ks = s ** 2 / 6 + (17 if mut('a3_s') else 18)   # mutation a3_s: s^2 = 6(K-17) instead of 6(K-18)
a = Ks / 32
b = Ks / 96
d = Ks * (Ks - 24) / 6144
r0 = 24 / (s + 6)
r1 = -24 / (s - 6)
quad = 4 * d * r ** 2 + 3 * b * r - 2 * a
res = cancel(together(quad - 4 * d * (r - r0) * (r - r1)))
res_root = cancel(quad.subs(r, r0))
prod_neg = cancel(-a / (2 * d) + 96 / (Ks - 24)) == 0   # product of the roots = -96/(K-24)
run.check('A3.4', res == 0 and res_root == 0 and prod_neg,
          'quadratic - 4d(r-r0)(r-r1): %s; quadratic(r0): %s; product of roots: %s' % (res, res_root, prod_neg))

# ------------------------------------------------------------------ A3.5 q(r0) closed form
qr0 = Q(3, 2) - a * r0 ** 2 + b * r0 ** 3 + d * r0 ** 4
off = 15 if mut('a3_qclosed') else 16   # mutation a3_qclosed: 18(K-16) -> 18(K-15)
closed = 4 * (Ks - 18) * (Ks * s - 18 * (Ks - off)) / (Ks - 24) ** 3
res_closed = cancel(together(qr0 - closed))
# DG09 (59), typed from the paper: chi(r0) >= 16(rho-6)((l^3-l)(l+2)-18)((l^3-l)(l+2)(rho-18)+288) / ((l^2+l+4)^4 (l-2)^4 (l+3)^4)
# with rho = sqrt(6 l(l^2-1)(l+2) - 108); chi(r0) = r0 q(r0). Reduce modulo rho^2 = 6K(l) - 108.
Kp = (l ** 3 - l) * (l + 2)
dg59 = 16 * (rho - 6) * (Kp - 18) * (Kp * (rho - 18) + 288) / ((l ** 2 + l + 4) ** 4 * (l - 2) ** 4 * (l + 3) ** 4)
r0_l = 4 * (rho - 6) / ((l - 2) * (l + 3) * (l ** 2 + l + 4))   # DG09 (58)
Kl2 = K_of(l)
qr0_l = 4 * (Kl2 - 18) * (Kl2 * rho - 18 * (Kl2 - 16)) / ((l - 2) * (l + 3) * (l ** 2 + l + 4)) ** 3   # RESULTS closed form
diff59 = together(r0_l * qr0_l - dg59)
num59 = sp.numer(diff59)
rem59 = sp.rem(Poly(sp.expand(num59), rho), Poly(rho ** 2 - (6 * Kl2 - 108), rho))
rem_ok = sp.expand(rem59.as_expr()) == 0
run.check('A3.5', res_closed == 0 and rem_ok,
          'closed form residual %s; DG09 (59) vs r0*q(r0) remainder %s' % (res_closed, rem59.as_expr() if not rem_ok else 0))

# ------------------------------------------------------------------ A3.6 the identity (Ks)^2 - (18(K-16))^2 = 6(K-24)^3
shift = 23 if mut('a3_ident') else 24   # mutation a3_ident: 6(K-24)^3 -> 6(K-23)^3
ident = (K ** 2 * 6 * (K - 18)) - (18 * (K - 16)) ** 2 - 6 * (K - shift) ** 3   # s^2 = 6(K-18)
run.check('A3.6', sp.expand(ident) == 0, 'residual %s' % sp.expand(ident))

# ------------------------------------------------------------------ A3.7 sign logic, exact, l = 3..50
bad = None
for lv in range(3, 51):
    Kv = K_of(lv)
    # K s > 18(K-16) > 0 iff K^2 * 6(K-18) > (18(K-16))^2 (both sides positive); then q(r0) = 4(K-18)(...)/(K-24)^3 > 0
    c1 = Kv > 24 and Kv > 18 and 18 * (Kv - 16) > 0
    c2 = Kv ** 2 * 6 * (Kv - 18) > (18 * (Kv - 16)) ** 2
    c3 = (Kv ** 2 * 6 * (Kv - 18) - (18 * (Kv - 16)) ** 2) == 6 * (Kv - 24) ** 3
    if not (c1 and c2 and c3):
        bad = lv
        break
run.check('A3.7', bad is None, 'sign logic fails at l=%s' % bad)
# numeric cross-check of q(r0) > 0 and equality with the closed form at 40 digits (not the decision)
worst = None
for lv in (3, 4, 10, 50):
    Kv = K_of(lv)
    sv = sp.sqrt(6 * (Kv - 18))
    r0v = 24 / (sv + 6)
    qv = Q(3, 2) - Q(Kv, 32) * r0v ** 2 + Q(Kv, 96) * r0v ** 3 + cn(5, lv) * r0v ** 4
    cl = 4 * (Kv - 18) * (Kv * sv - 18 * (Kv - 16)) / (Kv - 24) ** 3
    dif = abs(sp.N(qv - cl, 50))
    worst = dif if worst is None else max(worst, dif)
    run.note('l=%d: q(r0) = %s' % (lv, sp.N(qv, 12)))
run.check('A3.7x', worst < sp.Float('1e-40'), 'q(r0) vs closed form differ by %s' % worst)

# ------------------------------------------------------------------ A3.8 l = 2
P2 = chi_over_r(2).as_expr()
dP2 = sp.diff(P2, r)
crit = sp.solve(dP2, r)
rmin = 3 if mut('a3_l2_min') else 2   # mutation a3_l2_min: claim the minimum is at r = 3
ok_form = sp.expand(P2 - (Q(3, 2) - Q(3, 4) * r ** 2 + Q(1, 4) * r ** 3)) == 0
ok_min = (P2.subs(r, rmin) == Q(1, 2)) and (2 in crit)
lim_r0 = sp.limit(24 / (s + 6), s, 6)
Kq = s ** 2 / 6 + 18
lim_q = sp.limit(cancel(together(4 * (Kq - 18) * (Kq * s - 18 * (Kq - 16)) / (Kq - 24) ** 3)), s, 6)
run.check('A3.8', ok_form and ok_min and lim_r0 == 2 and lim_q == Q(1, 2),
          'form %s; min at %d value %s; limits r0=%s q=%s' % (ok_form, rmin, P2.subs(r, rmin), lim_r0, lim_q))

# ------------------------------------------------------------------ A3.9 exact Sturm: no positive root of chi/r, l = 2..50
for lv in range(2, 51):
    P = chi_over_r(lv)
    if mut('a3_sturm_neg'):   # mutation a3_sturm_neg: count the roots of P(-r), i.e. the negative roots of chi/r
        P = Poly(P.as_expr().subs(r, -r), r, domain='QQ')
    t = time.time()
    n1 = positive_root_count(P)
    n2 = P.count_roots(0, sp.oo)
    run.check('A3.9[l=%d]' % lv, n1 == 0 and n2 == 0, 'Sturm count %s, count_roots %s' % (n1, n2))
run.note('elapsed %.1f s' % (time.time() - T0))
run.finish()

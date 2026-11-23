"""#142: V-hat in the tortoise coordinate x near x = 0 against Dotti-Gleiser's eq. (61) (DG09, arXiv:0809.3615).

DG09 (61): V-hat = 12 M^2 r^-4 - 2M(l^2+l+3) r^-3 + O(r^-2) = 3/(4x^2) + (l^2+l-1)/(4|M|^(1/2) x^(3/2))
           - l(l+1)/(4|M| x) + O(x^(-1/2)), printed at M < 0; here M = -1.
E2.6  compares each printed coefficient with the series of RESULTS' V-hat (published(): a difference is reported as DIFFERS).
E2.6b states what the series is (checked exactly): the coefficient of 1/x is 0, and where DG's -l(l+1)/4 comes from.
E2.6c numeric cross-check (mpmath, 50 digits) of E2.6b.
Run: python review_lemma_a3_dg61.py [--mutate=NAME] > review_lemma_a3_dg61_output.txt
Mutations: pl_taylor, e_dg61b, cn_n4to3.
"""
import sys
import time
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from review_lemma_a3_common import *
import mpmath as mp

run = Run('review_lemma_a3_dg61: V-hat in x near 0 against DG09 (61)')
T0 = time.time()
planted_taylor(run)
s = Symbol('s', positive=True)
c, h2 = symbols('c h2')

# r^4 V-hat = 12 + A r + B r^2 + O(r^3) for psi = r (3/2 + h2 r^2 + h3 r^3)/(cr-6) (E2.1 of the ends script)
h3 = Symbol('h3')
psi_gen = r * (Q(3, 2) + h2 * r ** 2 + h3 * r ** 3) / (c * r - 6)
Vg = cancel(together(r ** 4 * Vhat_of_psi(psi_gen)))
ng, dg = sp.fraction(Vg)
tg = taylor(Poly(ng, r), Poly(dg, r), 3)
A = cancel(tg[1].subs(c, c_of(l)))
B = cancel(tg[2].subs({c: c_of(l), h2: -K_of(l) / 32}))
run.note('r^4 V-hat = 12 + (%s) r + (%s) r^2 + ...' % (sp.factor(A), sp.factor(B)))

# invert x(r) = r - 2 ln(1 + r/2) = s^2: r = a1 s + a2 s^2 + a3 s^3 + a4 s^4
a1, a2, a3, a4 = symbols('a1 a2 a3 a4')
rs_ser = a1 * s + a2 * s ** 2 + a3 * s ** 3 + a4 * s ** 4
xr = sp.series(r - 2 * sp.log(1 + r / 2), r, 0, 7).removeO()
eqs = sp.expand(sp.series(xr.subs(r, rs_ser), s, 0, 6).removeO() - s ** 2)
sol = sp.solve([eqs.coeff(s, j) for j in range(2, 6)], [a1, a2, a3, a4], dict=True)
sol = [d_ for d_ in sol if d_[a1] > 0][0]
r_of_s = rs_ser.subs(sol)
run.note('r(s) = %s, s = sqrt(x)' % sp.series(r_of_s, s, 0, 5))


def coeffs_of(Vr):
    """Coefficients of s^-4, s^-3, s^-2 of Vr(r(s))."""
    ser = sp.expand(sp.series(cancel(together(s ** 4 * Vr.subs(r, r_of_s))), s, 0, 3).removeO())
    return {j: ser.coeff(s, j + 4) for j in (-4, -3, -2)}


full = coeffs_of((12 + A * r + B * r ** 2) / r ** 4)
two = coeffs_of((12 + A * r) / r ** 4)          # the first two terms of the r-series only
only_B = coeffs_of(B / r ** 2)                 # the O(r^-2) term alone
printed = {-4: Q(3, 4), -3: (l ** 2 + l - 1) / 4, -2: -l * (l + 1) / 4}
LAB = {-4: 'x^-2', -3: 'x^-3/2', -2: 'x^-1'}
for j in (-4, -3, -2):
    run.published('E2.6[%s]' % LAB[j], sp.expand(full[j] - printed[j]) == 0,
                  'coefficient of %s: series gives %s, DG09 (61) prints %s' % (LAB[j], sp.factor(full[j]), sp.factor(printed[j])))

true_m2 = -l * (l + 1) / 4 if mut('e_dg61b') else 0   # mutation e_dg61b: claim the 1/x coefficient is DG's printed -l(l+1)/4
ok = (sp.expand(full[-2] - true_m2) == 0
      and sp.expand(two[-2] - printed[-2]) == 0      # DG's printed third term = the first two terms' contribution
      and sp.expand(only_B[-2] - l * (l + 1) / 4) == 0   # and the O(r^-2) term adds +l(l+1)/4 at the same order
      and sp.expand(full[-4] - Q(3, 4)) == 0 and sp.expand(full[-3] - (l ** 2 + l - 1) / 4) == 0)
run.check('E2.6b', ok, 'x^-1 coefficient %s (claimed %s); first two terms give %s, the r^-2 term gives %s' % (
    sp.factor(full[-2]), true_m2, sp.factor(two[-2]), sp.factor(only_B[-2])))

# E2.6c numeric cross-check
mp.mp.dps = 50
worst = 0
for lv in (2, 3, 5):
    cv = c_of(lv)
    V = cancel(together(Vhat_of_psi(chi_poly(lv).as_expr() / (cv * r - 6))))
    Vf = sp.lambdify(r, V, 'mpmath')
    rr = mp.mpf('1e-8')
    x = rr - 2 * mp.log(1 + rr / 2)
    rest = x * (Vf(rr) - 3 / (4 * x ** 2) - (lv ** 2 + lv - 1) / (4 * x ** mp.mpf(1.5)))
    worst = max(worst, abs(rest))
    run.note('l=%d: x (V-hat - 3/(4x^2) - (l^2+l-1)/(4 x^1.5)) at r = 1e-8: %s (DG09 prints %s)' % (lv, mp.nstr(rest, 6), -lv * (lv + 1) / 4))
claimed = 0 if not mut('e_dg61b') else 1
run.check('E2.6c', (worst < mp.mpf('1e-6')) == (claimed == 0), 'largest |x(V-hat - ...)| = %s, below 1e-6' % mp.nstr(worst, 6))
run.note('elapsed %.1f s' % (time.time() - T0))
run.finish()

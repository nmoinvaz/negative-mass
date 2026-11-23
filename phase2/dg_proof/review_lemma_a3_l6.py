"""#142: L6 of RESULTS.md (Psi-hat = A* Psi_DG is O(r^3) at 0 and regular at r_s), from RESULTS.md alone.

Psi_DG = r (r+2)^k e^{-kr/2}/(cr-6), psi_0 = chi/(cr-6), g = psi_0,x/psi_0, A* = -d/dx + g, f = 1 + 2/r.
With W = r (r+2)^k e^{-kr/2} = r E: A* Psi_DG = -f psi_0 (Psi_DG/psi_0)_r = -f (W' chi - W chi')/((cr-6) chi).
Run: python review_lemma_a3_l6.py [--mutate=NAME] > review_lemma_a3_l6_output.txt
Mutations: pl_wronsk, pl_order, pl_sturm_sign, l6_gen, l6_w, l6_c4, cn_n4to3, l6_k, l6_rs, l6_chip.
"""
import sys
import time
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from review_lemma_a3_common import *

run = Run('review_lemma_a3_l6: L6, Psi-hat = A* Psi_DG')
T0 = time.time()
LMAX = 30
c, k = symbols('c k', positive=True)
NSER = 8


def wr_order(W, chi, n=NSER):
    """(order, coefficient) at r = 0 of W' chi - W chi' (series)."""
    sgn = 1 if mut('pl_wronsk') else -1   # mutation pl_wronsk: W' chi + W chi'
    num = sp.diff(W, r) * chi + sgn * W * sp.diff(chi, r)
    ser = sp.series(num, r, 0, n).removeO()
    ser = sp.expand(ser)
    for j in range(n):
        co = ser.coeff(r, j)
        if sp.simplify(co) != 0:
            return j, sp.simplify(co)
    return None, 0


# ------------------------------------------------------------------ planted problems
run.planted('PLW.1', wr_order(r, r + r ** 3) == (3, -2), 'got %s' % (wr_order(r, r + r ** 3),))
run.planted('PLW.2', wr_order(r + 2 * r ** 2, r ** 2) == (2, -1), 'got %s' % (wr_order(r + 2 * r ** 2, r ** 2),))
run.planted('PLW.3', wr_order(sp.exp(r), sp.exp(2 * r)) == (0, -1), 'got %s' % (wr_order(sp.exp(r), sp.exp(2 * r)),))

planted_lead(run)
planted_sturm(run)

# ------------------------------------------------------------------ L6.1 generic algebra
Wf, chif = sp.Function('W')(r), sp.Function('chi')(r)
psif, Psif = chif / (c * r - 6), Wf / (c * r - 6)
g = dx(psif) / psif
Astar = -dx(Psif) + g * Psif
sg = 1 if mut('l6_gen') else -1   # mutation l6_gen: +f instead of -f
formula = sg * f_of_r * (sp.diff(Wf, r) * chif - Wf * sp.diff(chif, r)) / ((c * r - 6) * chif)
res = sp.simplify(Astar - formula)
run.check('L6.1', res == 0, 'A* Psi - formula = %s' % res)

# ------------------------------------------------------------------ L6.2 series of W / 2^k
wser = sp.series(r * sp.exp(k * (sp.log(1 + r / 2) - r / 2)), r, 0, 6).removeO()
want = {1: 1, 2: 0, 3: -k / 8, 4: (k / 12 if mut('l6_w') else k / 24), 5: k ** 2 / 128 - k / 64}   # mutation l6_w: k/12
okw = all(sp.expand(wser.coeff(r, j) - want[j]) == 0 for j in want)
dW = cancel(sp.diff(r * (r + 2) ** k * sp.exp(-k * r / 2), r) / ((r + 2) ** k * sp.exp(-k * r / 2)) - (r * (k / (r + 2) - k / 2) + 1))
run.check('L6.2', okw and sp.simplify(dW) == 0, 'series %s; W\'/E - (r e + 1) = %s' % (wser, sp.simplify(dW)))

# ------------------------------------------------------------------ L6.3 symbolic in k: Wronskian orders and the leading coefficient
Kk = 6 * k
c4 = Kk / (90 if mut('l6_c4') else 96)   # mutation l6_c4: K/90 instead of K/96
chi_ser = Q(3, 2) * r + 0 * r ** 2 - Kk / 32 * r ** 3 + c4 * r ** 4 + Kk * (Kk - 24) / 6144 * r ** 5
w_ser = sum(want_ * r ** j for j, want_ in {1: 1, 2: 0, 3: -k / 8, 4: k / 24, 5: k ** 2 / 128 - k / 64}.items())
num = sp.expand(sp.diff(w_ser, r) * chi_ser - w_ser * sp.diff(chi_ser, r))
coef = {p: sp.expand(num.coeff(r, p)) for p in range(1, 6)}
rho_hat = -f_of_r * (sp.diff(w_ser, r) * chi_ser - w_ser * sp.diff(chi_ser, r)) / ((c * r - 6) * chi_ser)
pw, cf = lead(sp.expand(rho_hat), 0) if False else (None, None)
# leading term: numerator lowest term over denominator lowest term (truncation errors are higher order)
lead_num = [p for p in range(1, 6) if coef[p] != 0]
if lead_num:
    p5 = lead_num[0]
    pw = p5 - 1 - 1 - 1 + 0        # r^p / (r * r): f ~ 2/r, chi ~ (3/2) r
    cf = cancel(-2 * coef[p5] / ((-6) * Q(3, 2)))
    pw = p5 - 2
else:
    pw, cf = None, None
ok = all(coef[p] == 0 for p in (1, 2, 3, 4)) and sp.expand(coef[5] - 3 * k ** 2 / 128) == 0 and pw == 3 and sp.expand(cf - k ** 2 / 192) == 0
run.check('L6.3', ok, 'Wronskian coefficients p=1..5: %s; Psi-hat/2^k ~ %s r^%s (want k^2/192 r^3)' % (
    [sp.factor(coef[p]) for p in range(1, 6)], sp.factor(cf) if cf is not None else None, pw))

# ------------------------------------------------------------------ L6.4 exact: l = 2..LMAX with the actual chi and k
RH = {}
for lv in range(2, LMAX + 1):
    cv = c_of(lv)
    kv = k_of(lv) + (1 if mut('l6_k') else 0)   # mutation l6_k: k + 1
    chi = chi_poly(lv).as_expr()
    e = kv / (r + 2) - kv / 2
    RH[lv] = cancel(together(-f_of_r * ((r * e + 1) * chi - r * sp.diff(chi, r)) / ((cv * r - 6) * chi)))
    pw, cf = lead(RH[lv], 0)
    run.check('L6.4[l=%d]' % lv, (pw, cf) == (3, kv ** 2 / 192) and kv == sp.Rational(K_of(lv) + (6 if mut('l6_k') else 0), 6),
              'rho-hat ~ %s r^%s at 0; want (k^2/192) r^3, k = %s' % (cf, pw, kv))

# ------------------------------------------------------------------ L6.5 regularity at r_s, symbolic in c
rs5 = Q(5) / c if mut('l6_rs') else 6 / c   # mutation l6_rs: r_s = 5/c
k5 = c * (c + 2) / 6
Wlog = 1 / rs5 + k5 / (rs5 + 2) - k5 / 2
chilog = -1 / (rs5 * (rs5 + 2))
res_ok = cancel(Wlog - chilog)
k_other = k5 + 1
res_bad = cancel(1 / rs5 + k_other / (rs5 + 2) - k_other / 2 - chilog)
run.check('L6.5', res_ok == 0 and res_bad != 0, "W'/W(r_s) - chi'/chi(r_s) = %s (k = c(c+2)/6), %s (k+1)" % (res_ok, res_bad))

# ------------------------------------------------------------------ L6.6 exact: no pole at r_s, none on (0, oo), degree
for lv in range(2, LMAX + 1):
    cv = c_of(lv)
    chi = chi_poly(lv)
    rsv = Q(6, cv)
    sh = 3 if mut('l6_chip') else 2   # mutation l6_chip: -1/(r_s (r_s+3))
    ok1 = cancel(chi.diff(r).eval(rsv) / chi.eval(rsv) + 1 / (rsv * (rsv + sh))) == 0
    num_, den_ = sp.fraction(RH[lv])
    pd = Poly(den_, r)
    ok2 = pd.eval(rsv) != 0 and pd.eval(0) != 0
    ok3 = positive_root_count(pd) == 0
    pinf, cinf = lead(RH[lv], 'inf')
    kv = k_of(lv)
    ok4 = kv > 0 and sp.fraction(RH[lv])[1] != 0
    run.check('L6.6[l=%d]' % lv, ok1 and ok2 and ok3 and ok4,
              "chi'/chi(r_s) relation %s; den(r_s) != 0 and den(0) != 0: %s; positive roots of den: %s; rho-hat ~ r^%s at inf" % (
                  ok1, ok2, positive_root_count(pd), pinf))
run.note('elapsed %.1f s' % (time.time() - T0))
run.finish()

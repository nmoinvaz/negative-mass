"""#142: the Zerilli potential V of RESULTS.md 4.2 against the published forms, at M = -1, and related printed formulas.

Sources (data; typed from the page images / text, see I142-rpi5/README.md): DG06 = Dotti-Gleiser CQG 23 (2006) 5063,
arXiv:gr-qc/0604021, eqs. (6), (8); DG09 = Dotti-Gleiser CQG 26 (2009) 215002, arXiv:0809.3615, eqs. (11), (13), (52)-(57) and
the Fig. 1 caption; KS = Kokkotas-Schmidt, Living Rev. Relativ. 2 (1999) 2, arXiv:gr-qc/9909058, eq. (24); BCS = Berti-Cardoso-
Starinets, CQG 26 (2009) 163001, arXiv:0905.2975, eqs. (14), (15) (the AdS radius of (15) taken to infinity).
Zerilli's papers (PRL 24, 737; PRD 2, 2141) and Chandrasekhar's book were not read: (unverified).
lam = (l-1)(l+2)/2 = c/2. M symbolic unless said.
Run: python review_lemma_a3_v.py [--mutate=NAME] > review_lemma_a3_v_output.txt
Mutations: pl_cmp, v_typo, v_dg06, v_ks, v_x, v_chand, v_dg54, v_dg53, v_cn, v_m2.
"""
import sys
import time
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from review_lemma_a3_common import *
import mpmath as mp

run = Run('review_lemma_a3_v: the Zerilli V of RESULTS 4.2 against the published forms')
T0 = time.time()
M, lam, c = symbols('M lam c')
chi_f = sp.Function('chi')(r)


def same(e1, e2):
    """Exactly equal as rational functions. Mutation pl_cmp: compare at r = 1 only."""
    d = e1 - e2
    if mut('pl_cmp'):
        d = d.subs(r, 1)
    return cancel(together(d)) == 0


# ------------------------------------------------------------------ planted problems
fM = 1 - 2 * M / r
VZ_poly = 2 * fM * (lam ** 2 * (lam + 1) * r ** 3 + 3 * lam ** 2 * M * r ** 2 + 9 * lam * M ** 2 * r + 9 * M ** 3) / (r ** 3 * (lam * r + 3 * M) ** 2)
VZ_fact = 2 * fM * (lam ** 2 * r ** 2 * ((lam + 1) * r + 3 * M) + 9 * M ** 2 * (lam * r + M)) / (r ** 3 * (lam * r + 3 * M) ** 2)
run.planted('PLV.1a', same(VZ_poly, VZ_fact), 'planted equal pair reported different')
run.planted('PLV.1b', not same(VZ_poly, VZ_poly + (r - 1) * M), 'planted unequal pair (differ away from r = 1) reported equal')
run.planted('PLV.2', cancel(fM.subs(M, -1) - (1 + 2 / r)) == 0, 'M -> -1 on 1 - 2M/r')

# ------------------------------------------------------------------ published forms (M, lam symbolic)
nine = 9 * M ** 2 * (lam * r + (2 * M if mut('v_dg06') else M))   # mutation v_dg06: 9M^2(lam r + 2M)
P_DG06 = 2 * (1 - 2 * M / r) * (lam ** 2 * r ** 2 * ((lam + 1) * r + 3 * M) + nine) / (r ** 3 * (lam * r + 3 * M) ** 2)         # DG06 (6)
P_DG09 = 2 * (1 - 2 * M / r) * (lam ** 2 * r ** 2 * ((lam + 1) * r + 3 * M) + 9 * M ** 2 * (lam * r + M)) / (r ** 3 * (lam * r + 3 * M) ** 2)   # DG09 (11)
nks = lam
ks18 = 18 * nks * M ** 2 * r if not mut('v_ks') else 18 * nks * M * r   # mutation v_ks: 18 n M r
P_KS = (1 - 2 * M / r) * (2 * nks ** 2 * (nks + 1) * r ** 3 + 6 * nks ** 2 * M * r ** 2 + ks18 + 18 * M ** 3) / (r ** 3 * (nks * r + 3 * M) ** 2)   # KS (24), n = lam
P_BCS = 2 * fM / r ** 3 * (9 * M ** 3 + 3 * lam ** 2 * M * r ** 2 + lam ** 2 * (1 + lam) * r ** 3 + 9 * M ** 2 * lam * r) / (3 * M + lam * r) ** 2   # BCS (15), L -> oo

# ------------------------------------------------------------------ V1.1 RESULTS' V at M = -1 against each published form
VR = V_results(lam)
for name, P in (('DG06', P_DG06), ('DG09', P_DG09), ('KS', P_KS), ('BCS', P_BCS)):
    ok_lam = same(P.subs(M, -1), VR)
    ok_l = same(P.subs(M, -1).subs(lam, (l - 1) * (l + 2) / 2), VR.subs(lam, (l - 1) * (l + 2) / 2))
    run.check('V1.1[%s]' % name, ok_lam and ok_l, 'published form at M = -1 differs from RESULTS V (symbolic in lam: %s, in l: %s)' % (ok_lam, ok_l))

# ------------------------------------------------------------------ V1.2 the published forms agree for generic M
pairs = [('DG06', P_DG06, 'DG09', P_DG09), ('DG09', P_DG09, 'KS', P_KS), ('DG09', P_DG09, 'BCS', P_BCS)]
bad = [(a, b) for a, Pa, b, Pb in pairs if not same(Pa, Pb)]
run.check('V1.2', not bad, 'forms differ: %s' % bad)

# ------------------------------------------------------------------ V1.3 tortoise coordinate
xR = r - 2 * sp.log(1 + r / 2)                                       # RESULTS
# For r > 0 and M < 0: |M| = -M and |(r-2M)/(2M)| = (r-2M)/(-2M) (r - 2M > 0 > 2M). The printed logs are read with these signs.
x06 = r + 2 * M * (sp.log(r - 2 * M) - sp.log(-2 * M))               # DG06 (8): r + 2M[ln(r-2M) - ln(2|M|)]
x09 = r + 2 * M * sp.log((r - 2 * M) / (-2 * M))                     # DG09 (13): r + 2M ln((r-2M)/(2M)), argument read as its modulus
x_rhs = r - 2 * sp.log(1 + r) if mut('v_x') else xR   # mutation v_x: ln(1 + r) instead of ln(1 + r/2)
inv_f = 1 / (1 - 2 * M / r)
dxr_06 = cancel(sp.diff(x06, r) - inv_f)
dxr_09 = cancel(sp.diff(x09, r) - inv_f)
d06 = (x06 - x_rhs).subs(M, -1)
d09 = (x09 - x_rhs).subs(M, -1)
ok06 = cancel(sp.diff(d06, r)) == 0 and sp.simplify(d06.subs(r, 0)) == 0
ok09 = cancel(sp.diff(d09, r)) == 0 and sp.simplify(d09.subs(r, 0)) == 0
okR = cancel(sp.diff(x_rhs, r) - 1 / f_of_r) == 0 and x_rhs.subs(r, 0) == 0
run.check('V1.3', dxr_06 == 0 and dxr_09 == 0 and ok06 and ok09 and okR,
          'dx/dr = 1/(1-2M/r): DG06 %s, DG09 %s; equal to RESULTS x at M = -1: %s, %s; RESULTS: %s' % (dxr_06, dxr_09, ok06, ok09, okR))

# ------------------------------------------------------------------ V1.4 Chandrasekhar's V^{+-} = +-beta F_x + beta^2 F^2 + kappa F (unverified source)
Delta = r ** 2 - 2 * M * r
mu2 = c                                                              # (l-1)(l+2) = 2 lam = c
F = Delta / (r ** 3 * (mu2 * r + 6 * M))
beta = (5 if mut('v_chand') else 6) * M    # mutation v_chand: beta = 5M
kappa = mu2 * (mu2 + 2)
Fx = fM * sp.diff(F, r)
Vplus = beta * Fx + beta ** 2 * F ** 2 + kappa * F
Vminus = -beta * Fx + beta ** 2 * F ** 2 + kappa * F
VRW = fM * ((c + 2) / r ** 2 - 6 * M / r ** 3)                      # BCS (14): f [l(l+1)/r^2 - 6M/r^3], l(l+1) = c + 2
okp = same(Vplus, P_DG09.subs(lam, c / 2))
okm = same(Vminus, VRW)
okp1 = same(Vplus.subs(M, -1), VR.subs(lam, c / 2))
run.check('V1.4', okp and okm and okp1, 'V+ = V_Z: %s; V- = V_RW: %s; V+ at M = -1 = RESULTS V: %s' % (okp, okm, okp1))

# ------------------------------------------------------------------ V1.5 V-hat from the definition against DG09 (54), chi'' eliminated with DG09 (53), M = -1
cc = Symbol('c', positive=True)
lamc = cc / 2
Mn = -1
p53 = (6 * Mn ** 2 + 2 * r * lamc * (3 * Mn - r)) / (r * (r - 2 * Mn) * (3 * Mn + lamc * r))
q53 = -(6 * Mn ** 2 + 2 * r * lamc * (3 * Mn + lamc * r)) / (r ** 2 * (r - 2 * Mn) * (3 * Mn + lamc * r))
chi2 = -p53 * sp.diff(chi_f, r) - q53 * chi_f                        # chi'' from (53)
psi_c = chi_f / (6 * Mn + r * (cc + 0) )                             # DG09 (52): chi/(6M + r (l+2)(l-1)), (l+2)(l-1) = c
Vh_def = Vhat_of_psi(psi_c).subs(sp.diff(chi_f, r, 2), chi2)
four = 3 if mut('v_dg54') else 4   # mutation v_dg54: 3 instead of 4 in the middle term
Vh_54 = (2 * (r - 2 * Mn) ** 2 / (r ** 2 * chi_f ** 2) * sp.diff(chi_f, r) ** 2
         - four * (r - 2 * Mn) ** 2 * lamc / (r ** 2 * (3 * Mn + lamc * r) * chi_f) * sp.diff(chi_f, r)
         - (r - 2 * Mn) * (6 * Mn ** 2 + lamc * r * (2 * lamc * r + 4 * Mn)) / (r ** 4 * (3 * Mn + lamc * r)))
xi, xi1 = symbols('xi xi1')
d54 = cancel(together((Vh_def - Vh_54).subs(sp.diff(chi_f, r), xi1).subs(chi_f, xi)))
run.check('V1.5', d54 == 0, 'V-hat (definition, chi\'\' from DG09 (53)) - DG09 (54) = %s' % d54)

# ------------------------------------------------------------------ V1.6 the chi-ODE derived from V and psi_0 equals DG09 (53), M = -1
p53m = (6 * Mn ** 2 * (5 if mut('v_dg53') else 6) / 6 + 2 * r * lamc * (3 * Mn - r)) / (r * (r - 2 * Mn) * (3 * Mn + lamc * r))   # mutation v_dg53: 5M^2
chi_a, chi_b, chi_c = symbols('chi_a chi_b chi_c')                   # chi, chi', chi''
Vc = V_results(lamc)
expr = dx(dx(psi_c)) - Vc * psi_c
expr = expr.subs(sp.diff(chi_f, r, 2), chi_c).subs(sp.diff(chi_f, r), chi_b).subs(chi_f, chi_a)
expr = cancel(together(expr))
co2 = sp.diff(expr, chi_c)
co1 = sp.diff(expr, chi_b)
co0 = sp.diff(expr, chi_a)
p_mine = cancel(co1 / co2)
q_mine = cancel(co0 / co2)
okp = cancel(p_mine - p53m) == 0
okq = cancel(q_mine - q53) == 0
run.check('V1.6', okp and okq, "chi'' + p chi' + q chi from V and psi_0 vs DG09 (53): p %s, q %s" % (okp, okq))

# ------------------------------------------------------------------ V1.7 M -> -1 is an ordinary substitution; the pole sits at r_s > 0 because M < 0
den = sp.fraction(cancel(together(P_DG09)))[1]
rat_in_M = P_DG09.is_rational_function(M)
root_M = sp.solve(lam * r + 3 * M, r)[0]
ok = rat_in_M and sp.degree(sp.Poly(den, M), M) >= 0 and cancel(root_M.subs(M, -1) - 3 / lam) == 0 and cancel(root_M.subs(M, -1).subs(lam, c_of(l) / 2) - 6 / c_of(l)) == 0 \
    and cancel(root_M.subs(M, 1) + 3 / lam) == 0
run.check('V1.7', ok, 'V rational in M: %s; pole at r = -3M/lam: at M = -1 equals 3/lam = 6/c; at M = +1 equals -3/lam (outside r > 2M)' % rat_in_M)

# ------------------------------------------------------------------ V1.8 c_n of RESULTS against DG09 (55), (56), (57) at M = -1
bad = None
for lv in range(2, 31):
    for n in range(1, lv + 3):
        nn4 = (n - 3) if mut('v_cn') else (n - 4)   # mutation v_cn: (n-3) in DG's bracket
        dg = sp.Integer(n - 2) * (nn4 * lv * (lv + 1) + n - 1) * sp.gamma(lv + n - 1) / (2 ** n * sp.gamma(n) ** 2 * sp.gamma(lv - n + 3))
        # RESULTS' c_n for the same index; (-M)^n = 1 at M = -1
        if sp.simplify(dg - chi_poly(lv).coeff_monomial(r ** n)) != 0:
            bad = (lv, n)
            break
    if bad:
        break
Kl = K_of(l)
dg57 = {1: Q(3, 2), 3: -l * (l + 2) * (l ** 2 - 1) / 32, 4: l * (l + 2) * (l ** 2 - 1) / 96,
        5: l * (l ** 2 + l + 4) * (l + 3) * (l ** 2 - 4) * (l ** 2 - 1) / 6144}
ok57 = all(cancel(cn(n, l) - dg57[n]) == 0 for n in dg57)
dg56 = Q(3, 2) * r - Q(3, 4) * r ** 3 + Q(1, 4) * r ** 4                  # DG09 (56) at M = -1
ok56 = sp.expand(chi_poly(2).as_expr() - dg56) == 0
run.check('V1.8', bad is None and ok57 and ok56, 'c_n of RESULTS vs DG09 (55) differ at (l, n) = %s; (57) coefficients %s; (56) %s' % (bad, ok57, ok56))

# ------------------------------------------------------------------ V1.9 DG09 Fig. 1 caption: the l = 2, M = -2 double pole at x ~ 0.761
Mn2 = -2
lam2 = Q((2 - 1) * (2 + 2), 2)
rs2 = -3 * Mn2 / lam2
x2 = rs2 + 2 * Mn2 * sp.log(sp.Abs((rs2 - 2 * Mn2) / (2 * Mn2)))
caption = Q(861, 1000) if mut('v_m2') else Q(761, 1000)   # mutation v_m2: 0.861
val = sp.N(x2, 30)
run.check('V1.9', rs2 == 3 and abs(val - caption) < Q(1, 1000), 'x(r_s) = %s at r_s = %s, caption %s' % (val, rs2, sp.N(caption)))
run.note('elapsed %.1f s' % (time.time() - T0))
run.finish()

"""
Theorem 3 (#2): for every l >= 2, V^ has exactly one zero r_neg beyond r_*, it is simple, V^ < 0 on [r_*, r_neg) and
V^ > 0 on (r_neg, oo); in particular r_* < r_neg.

Proof ingredients, all checked here.  With g = psi0_x/psi0 (x-derivative), V the static Zerilli potential (M = -1):
  (Z)  psi0_xx = V psi0, so g_x = V - g^2 and V^ = g^2 - g_x = 2 g^2 - V.
  (P1) V > 0 on (r_s, oo).                                   [coefficient positivity, symbolic in l]
  (P2) V_r < 0 on (r_s, oo), l >= 2 real.                    [coefficient positivity in (tau, m), r = (6+tau)/c, l = 2+m]
  (B)  barrier lemma (by hand, in RESULTS.md): h = sqrt(V/2); at a zero of u = g - h on (r_*, oo),
       u_x = V/2 - V_x/(4h) > 0, so u has at most one zero there; u(r_*) = -h < 0 and V^ ~ (l+1)(l+2)/r^2 > 0 at infinity.
  (S)  V^ = c^4 (2 + eps rho) S(rho) / (rho^3 (rho-6) X^2) (identity for a generic X), so sign V^ = sign S on rho > 6.
Controls: a planted toy (psi = x^2 + b, known zero sqrt(b/3)), exact zero counts of S for l = 2..16, the independent
numerator of phi0_xx from chi_poly for l = 2..12, the identity (S) against chi_poly for l = 2..9.
Mutations: --mutate=vcoef | qsign | ssign | toy.
"""
import argparse
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'dg_spectrum'))
sys.path.insert(0, HERE)
import sympy as sp
import rescaled as R
from common import chi_poly, rs_

ap = argparse.ArgumentParser()
ap.add_argument('--mutate', default='')
args = ap.parse_args()
MUT = set(m for m in args.mutate.split(',') if m)
failures = []


def check(name, ok):
    print(("PASS  " if ok else "FAIL  ") + name, flush=True)
    if not ok:
        failures.append(name)


l, r, tau, m = sp.symbols('l r tau m', positive=True)
c = (l + 2) * (l - 1)
lam = c / 2
f = 1 + 2 / r
nine = 8 if 'vcoef' in MUT else 9
V = 2 * f * (lam ** 2 * r ** 2 * ((lam + 1) * r - 3) + nine * (lam * r - 1)) / (r ** 3 * (lam * r - 3) ** 2)

print("== (Z) psi0 = chi/(c r - 6) solves psi_xx = V psi with this V (exact, l = 2..14) ==")
okZ = True
for ll in range(2, 15):
    cc = (ll + 2) * (ll - 1)
    chi = chi_poly(ll)
    psi = chi / (cc * rs_ - 6)
    Vl = V.subs(l, ll).subs(r, rs_)
    resid = sp.cancel(sp.together((1 + 2 / rs_) * sp.diff((1 + 2 / rs_) * sp.diff(psi, rs_), rs_) - Vl * psi))
    if resid != 0:
        okZ = False
check("Z: f (f psi_r)_r - V psi = 0 identically in r for l = 2..14 (so V^ = 2 g^2 - V)", okZ)

print("\n== (P1) V > 0 for r > r_s ==")
sig, mu, lamS = sp.symbols('sigma mu lambda', positive=True)       # lam r = 3 + sigma, mu = 1/lam
nz = 8 if 'vcoef' in MUT else 9
rr_ = (3 + sig) / lamS
NZs = sp.expand(sp.simplify((lamS ** 2 * rr_ ** 2 * ((lamS + 1) * rr_ - 3) + nz * (lamS * rr_ - 1)).subs(lamS, 1 / mu)))
pol = sp.Poly(NZs, sig, mu)
cof = pol.coeffs()
print("   numerator of V at r = (3+sigma)/lam, mu = 1/lam:", sp.factor(NZs))
check("P1: the numerator of V is a polynomial in (sigma, mu) with all coefficients positive, so it is > 0 for lam r >= 3", all(co > 0 for co in cof))
print("   f > 0 and (lam r - 3)^2 > 0 off r_s, so V > 0 on (r_s, oo).")

print("\n== (P2) V_r < 0 for r > r_s ==")
num, den = sp.fraction(sp.factor(sp.together(sp.diff(V, r))))
print("   denominator of V_r:", sp.factor(den))
Q = sp.expand(-num / 2)
den_ok = sp.simplify(den - r ** 5 * (c * r - 6) ** 3) == 0
Qt = sp.expand(sp.cancel(Q.subs(r, (6 + tau) / c) * c ** 5).subs(l, 2 + m))
PQ = sp.Poly(Qt, tau, m)
coef = PQ.coeffs()
if 'qsign' in MUT:
    coef = coef[:5] + [-coef[5]] + coef[6:]
print(f"   Q(tau, m) = c^5 Q((6+tau)/c) with l = 2+m: {len(coef)} coefficients, min = {min(coef)}")
check("P2: denominator of V_r is r^5 (c r - 6)^3 > 0 on r > r_s", den_ok)
check("P2: every coefficient of Q(tau, m) is positive, so V_r = -2Q/(r^5 (c r - 6)^3) < 0 for all real l >= 2, r > r_s",
      all(co > 0 for co in coef) and den_ok)

print("\n== (S) sign V^ = sign S: identity for a generic X(rho) ==")
rho_s, cc_s = sp.symbols('rho c', positive=True)
Xf = sp.Function('X')
eps_s = 1 / cc_s
X0 = Xf(rho_s)
phi = cc_s * (rho_s - 6) / X0
Dx = lambda F: (1 + 2 * cc_s / rho_s) * cc_s * sp.diff(F, rho_s)          # d/dx = f d/dr = f c d/drho
Vh_gen = sp.simplify(Dx(Dx(phi)) / phi)
x1, x2 = sp.diff(X0, rho_s), sp.diff(X0, rho_s, 2)
Nn = X0 - (rho_s - 6) * x1
Np = -(rho_s - 6) * x2
sgn = -1 if 'ssign' in MUT else 1
Sg = rho_s * X0 * (eps_s * Nn + sgn * (2 + eps_s * rho_s) * Np) - (2 + eps_s * rho_s) * Nn * (X0 + 2 * rho_s * x1)
rhs = cc_s ** 4 * (2 + eps_s * rho_s) * Sg / (rho_s ** 3 * (rho_s - 6) * X0 ** 2)
check("S: V^ = c^4 (2 + rho/c) S / (rho^3 (rho-6) X^2) for generic X (symbolic in c)", sp.simplify(Vh_gen - rhs) == 0)

print("\n== (S) against chi_poly: V^ = g^2 - g_x computed from chi directly, l = 2..9 ==")
okS = True
for ll in range(2, 10):
    cc = (ll + 2) * (ll - 1)
    a, eps = R.a_dict(ll)
    rho = sp.Symbol('rho')
    Xp = sum(sp.Rational(v.numerator, v.denominator) * rho ** n for n, v in a.items())
    e = sp.Rational(eps.numerator, eps.denominator)
    X1p, X2p = sp.diff(Xp, rho), sp.diff(Xp, rho, 2)
    Np_ = Xp - (rho - 6) * X1p
    Npp = -(rho - 6) * X2p
    S_ = rho * Xp * (e * Np_ + sgn * (2 + e * rho) * Npp) - (2 + e * rho) * Np_ * (Xp + 2 * rho * X1p)
    chi = chi_poly(ll)
    psi = chi / (cc * rs_ - 6)
    g = (1 + 2 / rs_) * sp.diff(psi, rs_) / psi
    Vh = g ** 2 - (1 + 2 / rs_) * sp.diff(g, rs_)
    rhs_l = cc ** 4 * (2 + e * rho) * S_ / (rho ** 3 * (rho - 6) * Xp ** 2)
    if sp.cancel(sp.together(Vh - rhs_l.subs(rho, cc * rs_))) != 0:
        okS = False
check("S: the polynomial S built from rescaled.a_dict reproduces V^ from chi_poly exactly, l = 2..9", okS)

print("\n== (Psi) V^ = 2 g^2 - V with g = -c f N/((rho-6) X): the quadratic form Psi used in large_l_proof.py ==")
cS = sp.Symbol('c', positive=True)
rhoP = sp.Symbol('rho', positive=True)
epsP = 1 / cS
lamP = cS / 2
rP = rhoP / cS
fP = 1 + 2 / rP
nine_P = 8 if 'vcoef' in MUT else 9
V_c = 2 * fP * (lamP ** 2 * rP ** 2 * ((lamP + 1) * rP - 3) + nine_P * (lamP * rP - 1)) / (rP ** 3 * (lamP * rP - 3) ** 2)
NZ = (1 + 2 * epsP) * rhoP ** 3 / 8 - 3 * rhoP ** 2 / 4 + 9 * rhoP / 2 - 9
V_rho = 8 * cS ** 4 * (2 + epsP * rhoP) * NZ / (rhoP ** 4 * (rhoP - 6) ** 2)
check("Psi1: V(r = rho/c) = 8 c^4 (2 + rho/c) N_Z / (rho^4 (rho-6)^2) with N_Z = (1+2/c) rho^3/8 - 3 rho^2/4 + 9 rho/2 - 9",
      sp.simplify(V_c - V_rho) == 0)
e_s = sp.Symbol('e', positive=True)
NZ_tau = sp.expand(((1 + 2 * e_s) * (6 + tau) ** 3 / 8 - 3 * (6 + tau) ** 2 / 4 + 9 * (6 + tau) / 2 - 9))
check("Psi1: N_Z(rho = 6 + tau) is a polynomial in (tau, 1/c) with positive coefficients, so N_Z > 0 for rho >= 6",
      all(co > 0 for co in sp.Poly(NZ_tau, tau, e_s).coeffs()))
XP = sp.Function('X')(rhoP)
NP = XP - (rhoP - 6) * sp.diff(XP, rhoP)
phiP = cS * (rhoP - 6) / XP
DxP = lambda F: fP * cS * sp.diff(F, rhoP)
gP = -DxP(phiP) / phiP
check("Psi2: g = psi0_x/psi0 = -phi0_x/phi0 = -c f N / ((rho-6) X) for a generic X", sp.simplify(gP + cS * fP * NP / ((rhoP - 6) * XP)) == 0)
PsiP = rhoP ** 2 * (2 + epsP * rhoP) * NP ** 2 - 4 * NZ * XP ** 2
Psign = -1 if 'ssign' in MUT else 1
check("Psi3: 2 g^2 - V = 2 c^4 (2 + rho/c) Psi / (rho^4 (rho-6)^2 X^2), Psi = rho^2 (2+rho/c) N^2 - 4 N_Z X^2 (generic X)",
      sp.simplify(2 * gP ** 2 - V_rho - Psign * 2 * cS ** 4 * (2 + epsP * rhoP) * PsiP / (rhoP ** 4 * (rhoP - 6) ** 2 * XP ** 2)) == 0)
okPs = True
for ll in range(2, 10):
    cc = (ll + 2) * (ll - 1)
    a, eps = R.a_dict(ll)
    rho = sp.Symbol('rho')
    e = sp.Rational(eps.numerator, eps.denominator)
    Xp = sum(sp.Rational(v.numerator, v.denominator) * rho ** n for n, v in a.items())
    Np_ = Xp - (rho - 6) * sp.diff(Xp, rho)
    NZl = (1 + 2 * e) * rho ** 3 / 8 - 3 * rho ** 2 / 4 + 9 * rho / 2 - 9
    Psi_l = rho ** 2 * (2 + e * rho) * Np_ ** 2 - 4 * NZl * Xp ** 2
    chi = chi_poly(ll)
    psi = chi / (cc * rs_ - 6)
    g = (1 + 2 / rs_) * sp.diff(psi, rs_) / psi
    Vh = g ** 2 - (1 + 2 / rs_) * sp.diff(g, rs_)
    rhs_l = 2 * cc ** 4 * (2 + e * rho) * Psi_l / (rho ** 4 * (rho - 6) ** 2 * Xp ** 2)
    if sp.cancel(sp.together(Vh - rhs_l.subs(rho, cc * rs_))) != 0:
        okPs = False
check("Psi4: V^ computed from chi_poly equals 2 c^4 (2+rho/c) Psi / (rho^4 (rho-6)^2 X^2) exactly, l = 2..9", okPs)
print("   Hence sign V^ = sign Psi on rho > 6, and at a zero of N (r_*) Psi = -4 N_Z X^2 < 0, so V^(r_*) < 0 without any root counting.")

print("\n== planted truth for the barrier lemma: psi = x^2 + b (V = 2/psi, V^ = 2 (3x^2 - b)/(x^2+b)^2) ==")
x = sp.Symbol('x', positive=True)
b = sp.Integer(3)
psi_t = sp.exp(x ** 2 / 2) if 'toy' in MUT else x ** 2 + b
V_t = sp.simplify(sp.diff(psi_t, x, 2) / psi_t)
g_t = sp.diff(psi_t, x) / psi_t
Vh_t = sp.simplify(g_t ** 2 - sp.diff(g_t, x))
h_t = sp.sqrt(V_t / 2)
zeros = sp.solve(sp.Eq(Vh_t, 0), x)
print("   V =", V_t, " V^ =", sp.factor(Vh_t), " zeros of V^ on x > 0:", zeros)
check("toy: V_t is decreasing on x > 0 (hypothesis of the lemma)", bool(sp.simplify(sp.diff(V_t, x)).is_negative))
check("toy: V^ has exactly one zero on x > 0, at sqrt(b/3) = 1, as the lemma predicts",
      zeros == [1] or zeros == [sp.sqrt(b / 3)])
du = sp.simplify((sp.diff(g_t, x) - sp.diff(h_t, x)).subs(x, 1))
check("toy: u_x > 0 at the zero (transversal, upward)", bool(du > 0))

print("\n== exact zero counts of S for l = 2..16 (S built from rescaled.a_dict, sympy Sturm) ==")
rho = sp.Symbol('rho')
okcnt = True
okLC = True
okneg = True
for ll in range(2, 17):
    a, eps = R.a_dict(ll)
    e = sp.Rational(eps.numerator, eps.denominator)
    Xp = sum(sp.Rational(v.numerator, v.denominator) * rho ** n for n, v in a.items())
    X1p, X2p = sp.diff(Xp, rho), sp.diff(Xp, rho, 2)
    Np_ = sp.expand(Xp - (rho - 6) * X1p)
    Npp = -(rho - 6) * X2p
    S_ = sp.Poly(sp.expand(rho * Xp * (e * Np_ + sgn * (2 + e * rho) * Npp) - (2 + e * rho) * Np_ * (Xp + 2 * rho * X1p)), rho)
    Npoly = sp.Poly(Np_, rho)
    iso = [iv for iv in Npoly.intervals(eps=sp.Rational(1, 10 ** 12)) if iv[0][0] >= 6]
    hi = iso[0][0][1]
    n_beyond = S_.count_roots(hi, None)
    n_all = S_.count_roots(sp.Rational(1, 10 ** 9), None)
    if not (len(iso) == 1 and n_beyond == 1 and n_all == 3):
        okcnt = False
        print(f"   l = {ll}: unique N root {len(iso) == 1}, S roots beyond r_* {n_beyond}, positive S roots {n_all}")
    if not S_.eval(hi) < 0:
        okneg = False
    if not S_.LC() > 0:
        okLC = False
check("count: for l = 2..16, N has one root beyond rho = 6 and S has exactly one root in [rho_*, oo) (Theorem 3)", okcnt)
check("count: S(rho_*) < 0 and the leading coefficient of S is positive, l = 2..16 (V^ < 0 at r_*, V^ > 0 at infinity)", okneg and okLC)
print("   (S has exactly 3 positive roots for l = 2..16: rho = 6, the lower edge of {V^ < 0}, and r_neg c.  Only the root beyond r_*")
print("    is proved for all l; the count 3 is exact for the l listed and not claimed beyond them.)")

print("\n== independent check: largest zero of V^ from the numerator of phi0_xx (chi_poly), l = 2..12 ==")
okind = True
for ll in range(2, 13):
    cc = (ll + 2) * (ll - 1)
    chi = chi_poly(ll)
    phi0 = (cc * rs_ - 6) / chi
    ff = 1 + 2 / rs_
    numr = sp.numer(sp.together(sp.diff(ff * sp.diff(phi0, rs_), rs_)))
    Pn = sp.Poly(sp.expand(numr), rs_)
    big = [iv for iv in Pn.intervals(eps=sp.Rational(1, 10 ** 12)) if iv[0][0] > 0][-1][0]
    a, eps = R.a_dict(ll)
    e = sp.Rational(eps.numerator, eps.denominator)
    Xp = sum(sp.Rational(v.numerator, v.denominator) * rho ** n for n, v in a.items())
    X1p, X2p = sp.diff(Xp, rho), sp.diff(Xp, rho, 2)
    Np_ = Xp - (rho - 6) * X1p
    Npp = -(rho - 6) * X2p
    S_ = sp.Poly(sp.expand(rho * Xp * (e * Np_ + sgn * (2 + e * rho) * Npp) - (2 + e * rho) * Np_ * (Xp + 2 * rho * X1p)), rho)
    isoS = [iv for iv in S_.intervals(eps=sp.Rational(1, 10 ** 12)) if iv[0][0] > 0][-1][0]
    lo_s, hi_s = isoS[0] / cc, isoS[1] / cc
    if not (lo_s <= big[1] and big[0] <= hi_s):
        okind = False
        print(f"   l = {ll}: S root in [{float(lo_s)}, {float(hi_s)}], numerator root in [{float(big[0])}, {float(big[1])}]")
check("indep: largest root of S/c equals the largest root of the numerator of phi0_xx (chi_poly route), l = 2..12", okind)

sys.exit(1 if failures else 0)

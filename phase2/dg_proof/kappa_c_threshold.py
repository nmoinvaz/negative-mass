"""
kappa_c(R) = (d ln phi0/dx)(R) = (1 + 2/R) N(R)/((cR-6) chi(R)),  N = c chi - (c r - 6) chi'.
(a) closed form vs dg_spectrum/common.py; (b) R -> oo and R -> 0 expansions (l symbolic where possible);
(c) exact root isolation: N has no root on (0, r_s) and exactly one on (r_s, oo) [=> r_* unique], l = 2..50,
    via Descartes' rule on N(r_s + t) and Sturm counts.
"""
import sys
import sympy as sp
sys.path.insert(0, '/home/claude/negative-mass/phase2/dg_spectrum')
from common import chi_poly, rs_ as r, symbolic

l, R, t = sp.symbols('l R t', positive=True)

print("(a) closed form = common.symbolic(l)['dlogphi0_dx']:")
for ll in range(2, 11):
    c = (ll + 2)*(ll - 1); chi = chi_poly(ll)
    N = sp.expand(c*chi - (c*r - 6)*sp.diff(chi, r))
    kc = (1 + 2/r)*N/((c*r - 6)*chi)
    assert sp.cancel(kc - symbolic(ll)['dlogphi0_dx']) == 0
print("   identical for l = 2..10")

print("(b) large R: kappa_c = -(l+1)/R + a2/R^2 + O(R^-3)")
for ll in range(2, 11):
    c = (ll + 2)*(ll - 1); chi = chi_poly(ll)
    kc = sp.cancel((1 + 2/r)*(c/(c*r - 6) - sp.diff(chi, r)/chi))
    s = sp.series(kc.subs(r, 1/t), t, 0, 3).removeO()
    print(f"   l={ll}: {sp.expand(s.subs(t, 1/R))}")
# general l: chi = c_{l+2} R^{l+2} (1 + rho/R + ...), rho = c_{l+1}/c_{l+2}
def cn(nn, ll):
    LL = ll*(ll + 1); return sp.Integer(nn - 2)*((nn - 4)*LL + nn - 1)*sp.rf(ll - nn + 3, 2*nn - 4)/(2**nn*sp.factorial(nn - 1)**2)
rho = sp.Symbol('rho')
cc = (l + 2)*(l - 1)
kc_gen = (1 + 2/R)*(cc/(cc*R - 6) - ((l + 2)/R + (l + 1)*rho/R**2)/(1 + rho/R))
a2 = sp.factor(sp.series(kc_gen.subs(R, 1/t), t, 0, 3).removeO().coeff(t, 2))
print("   general l: a2 =", a2, " with rho = c_{l+1}/c_{l+2}")
rho_l = [sp.nsimplify(cn(ll + 1, ll)/cn(ll + 2, ll)) for ll in range(2, 12)]
print("   rho(l) =", rho_l)
# guess closed form of rho(l)
lsym = sp.Symbol('l')
cand = sp.interpolate(list(zip(range(2, 12), rho_l)), lsym)
print("   rho(l) interpolating poly/rational:", sp.factor(cand))
print("(b') small R: x_R = R^2/4 - R^3/12 + ..., kappa_c * x_R ->")
for ll in [2, 3, 5]:
    c = (ll + 2)*(ll - 1); chi = chi_poly(ll)
    kc = sp.cancel((1 + 2/r)*(c/(c*r - 6) - sp.diff(chi, r)/chi))
    x = r - 2*sp.log(1 + r/2)
    print(f"   l={ll}: kappa_c = {sp.series(kc, r, 0, 2)};  kappa_c*x = {sp.series(kc*x, r, 0, 3)}")

print("(c) roots of N:")
bad = []
for ll in range(2, 51):
    c = (ll + 2)*(ll - 1); chi = chi_poly(ll); rs = sp.Rational(6, c)
    N = sp.Poly(sp.expand(c*chi - (c*r - 6)*sp.diff(chi, r)), r)
    n_left = N.count_roots(0, rs); n_right = N.count_roots(rs, None)
    co = sp.Poly(sp.expand(N.as_expr().subs(r, rs + t)), t).all_coeffs()
    sg = [sp.sign(x) for x in co if x != 0]
    desc = sum(1 for a, b in zip(sg, sg[1:]) if a != b)
    # N(0) = 6 chi'(0) = 9 > 0
    if not (n_left == 0 and n_right == 1 and desc == 1):
        bad.append(ll)
    if ll <= 12 or ll % 10 == 0:
        iv = [x for x in N.intervals(eps=1e-14) if x[0][0] >= rs]
        print(f"   l={ll:2d}: roots in (0,r_s): {n_left}; in (r_s,oo): {n_right}; Descartes sign changes of N(r_s+t): {desc};"
              f"  r_* = {float(iv[0][0][0]):.12f}, c r_* = {float(c*iv[0][0][0]):.8f}", flush=True)
print("   exceptions for l = 2..50:", bad)

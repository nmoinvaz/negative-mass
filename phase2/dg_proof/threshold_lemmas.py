"""
All-l lemmas for the thresholds (l symbolic):
 T1  N(r) = c chi - (c r - 6) chi' > 0 on [0, r_s]          (lower bound by telescoping, then a cubic in u = r/r_s)
 T2  N has exactly one root r_* on (r_s, oo)                 (Descartes on N(r_s + t); needs chi''(r_s), chi'''(r_s) >= 0,
                                                               obtained in closed form from the chi-ODE at its singular point r_s)
 T3  large-R coefficient a2(l) of kappa_c = -(l+1)/R + a2/R^2 + ...
"""
import sys
import sympy as sp
sys.path.insert(0, '/home/claude/negative-mass/phase2/dg_spectrum')
from common import chi_poly, rs_
l, r, u = sp.symbols('l r u', positive=True)
c = (l + 2)*(l - 1); K = c*(c + 2); rs = 6/c
L = l*(l + 1)

print("== T1 ==")
print("K - (l-1)l(l+1)(l+2) =", sp.expand(K - (l - 1)*l*(l + 1)*(l + 2)))
c3 = -K/32; c4 = K/96
low = 9 + 18*c3*r**2 + (24*c4 - 2*c*c3 - 18*c4)*r**3
p = sp.expand(sp.simplify(low.subs(r, rs*u)))
print("lower bound N >= 9 + 18 c3 r^2 + (6 c4 - 2 c c3) r^3, in u = r/r_s:", sp.factor(p))
u0 = c/(c + 1)
print("critical point u0 = c/(c+1): p'(u0) =", sp.simplify(sp.diff(p, u).subs(u, u0)))
pm = sp.factor(sp.simplify(p.subs(u, u0)))
print("p(u0) =", pm, "=  9 - (27/4) c(c+2)/(c+1)^2 :", sp.simplify(pm - (9 - sp.Rational(27, 4)*c*(c + 2)/(c + 1)**2)))
print("p(0) = 9, p(1) =", sp.factor(p.subs(u, 1)), " -> min over [0,1] >= 9 - 27/4 > 0")
# exact check of the telescoping inequality for l = 2..12: N - low = sum_{j>=5} 6 c_j r^{j-1} + sum_{m>=4} c (m-1) c_m r^{m-1} (6/c - r)... verify N - low >= 0 on [0, r_s] by root counting
for ll in range(2, 13):
    cc = (ll + 2)*(ll - 1); chi = chi_poly(ll)
    N = sp.expand(cc*chi - (cc*rs_ - 6)*sp.diff(chi, rs_))
    lw = sp.expand(low.subs(l, ll).subs(r, rs_))
    d = sp.Poly(N - lw, rs_)
    # N - low = sum_{m>=4} c_m (m-1) r^(m-1) (6 - c r) + sum_{j>=5} 6 c_j r^(j-1)
    chk = sp.expand(N - lw - sum(sp.Poly(chi, rs_).coeff_monomial(rs_**m)*((m - 1)*rs_**(m - 1)*(6 - cc*rs_)) for m in range(4, ll + 3))
                    - sum(6*sp.Poly(chi, rs_).coeff_monomial(rs_**j)*rs_**(j - 1) for j in range(5, ll + 3)))
    assert chk == 0, (ll, chk)
print("identity N - low = sum_{m>=4} (m-1) c_m r^(m-1) (6 - c r) + sum_{j>=5} 6 c_j r^(j-1): verified exactly l = 2..12")

print("\n== T2: derivatives of chi at r_s from the chi-ODE  P2 chi'' + P1 chi' + P0 chi = 0,  P2(r_s) = 0 ==")
P2 = -r**2*(r + 2)*(c*r - 6)
P1 = 2*r*(l**2*r**2 + 3*l**2*r + l*r**2 + 3*l*r - 2*r**2 - 6*r - 6)
P0 = l**4*r**2 + 2*l**3*r**2 - 3*l**2*r**2 - 6*l**2*r - 4*l*r**2 - 6*l*r + 4*r**2 + 12*r + 12
# unknowns y_j = chi^(j)(r_s)/chi(r_s)
# indicial exponents of the chi-ODE at r_s: 0 and 1 - P1(r_s)/P2'(r_s)
expo = sp.simplify(1 - P1.subs(r, rs)/sp.diff(P2, r).subs(r, rs))
print("  exponents of the chi-ODE at its regular singular point r_s: 0 and", expo,
      " (resonant: chi3(r_s) is a free parameter, fixed only globally)")
y = sp.symbols('y0:4'); sol = {y[0]: 1}
for j in range(0, 2):
    expr = 0
    for (P, o) in ((P2, 2), (P1, 1), (P0, 0)):
        for i in range(0, j + 1):
            expr += sp.binomial(j, i)*sp.diff(P, r, i).subs(r, rs)*y[o + j - i]
    sol[y[j + 1]] = sp.factor(sp.simplify(sp.solve(sp.Eq(expr.subs(sol), 0), y[j + 1])[0]))
    print(f"  chi^({j+1})(r_s)/chi(r_s) =", sol[y[j + 1]])
for ll in range(2, 13):
    chi = chi_poly(ll); rr = sp.Rational(6, (ll + 2)*(ll - 1))
    d1 = sp.simplify(sp.diff(chi, rs_).subs(rs_, rr)/chi.subs(rs_, rr) - sol[y[1]].subs(l, ll))
    d2 = sp.diff(chi, rs_, 2).subs(rs_, rr)
    d3 = sp.diff(chi, rs_, 3).subs(rs_, rr)
    assert d1 == 0 and d2 == 0
    print(f"   l={ll:2d}: chi'(r_s)/chi(r_s) formula ok, chi''(r_s) = {d2}, chi3(r_s)/chi(r_s) = {float(d3/chi.subs(rs_, rr)):+.5g}")
print("  N(r_s + t) = c sum_m (1-m) d_m t^m with d_m = chi^(m)(r_s)/m!:  t^0: c chi(r_s) > 0;  t^1: 0;  t^2: -c d_2 = 0;")
print("  t^3: -2 c d_3 (either sign);  t^m, 4 <= m <= l+2: -(m-1) c d_m < 0 because chi^(m) has nonnegative coefficients")
print("  (c_n >= 0 for n >= 4) and c_{l+2} > 0.  Sign sequence (+, +/-, -, ..., -): exactly one change -> exactly one root.")
print("\n== T3: large-R expansion of kappa_c ==")
rho = (l - 1)*(l + 1)**2*((l - 3)*L + l)/(l**2*((l - 2)*L + l + 1))
for ll in range(2, 8):
    def cn(nn):
        LL = ll*(ll + 1); return sp.Integer(nn - 2)*((nn - 4)*LL + nn - 1)*sp.rf(ll - nn + 3, 2*nn - 4)/(2**nn*sp.factorial(nn - 1)**2)
    assert sp.simplify(cn(ll + 1)/cn(ll + 2) - rho.subs(l, ll)) == 0
a2 = -(2*l**3 - l**2*rho + 4*l**2 - l*rho - 2*l + 2*rho - 10)/((l - 1)*(l + 2))
print("  rho = c_{l+1}/c_{l+2} =", sp.factor(rho))
print("  a2(l) =", sp.factor(sp.simplify(a2)), ";  a2 values l=2..6:", [sp.nsimplify(sp.simplify(a2.subs(l, v))) for v in range(2, 7)])
print("  a2 ~", sp.limit(a2/l, l, sp.oo), "* l  as l -> oo")

print("\n== T4: small-R expansion of kappa_c, general l ==")
e5 = sp.Symbol('e5')
chi_s = sp.Rational(3, 2)*r - (K/32)*r**3 + (K/96)*r**4 + e5*r**5
kc = (1 + 2/r)*(c/(c*r - 6) - sp.diff(chi_s, r)/chi_s)
ser = sp.series(kc, r, 0, 1).removeO()
print("  kappa_c =", sp.factor(ser.coeff(r, -2)), "/R^2 +", sp.factor(ser.coeff(r, -1)), "/R + O(1)")
x = r - 2*sp.log(1 + r/2)
sx = sp.series(kc*x, r, 0, 2).removeO()
print("  kappa_c * x_R =", sp.factor(sx.coeff(r, 0)), "+", sp.factor(sx.coeff(r, 1)), "R + O(R^2)")

print("\n== T3': rho = c_{l+1}/c_{l+2} for symbolic l (Gamma form) ==")
def cg(nn):
    return (nn - 2)*((nn - 4)*L + nn - 1)*sp.gamma(l + nn - 1)/(2**nn*sp.gamma(nn)**2*sp.gamma(l - nn + 3))
print("  c_{l+1}/c_{l+2} - rho =", sp.simplify(sp.gammasimp(cg(l + 1)/cg(l + 2)) - rho))

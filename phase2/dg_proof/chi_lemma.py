"""
Lemma A (all l >= 2): chi(r) of DG eq. (55) [M = -1] solves the chi-ODE obtained from the static Zerilli
equation with psi0 = chi/(c r - 6), c = (l+2)(l-1); and chi > 0 on r > 0.

Part 1: closed form of the coefficients c_n(l) and their sign.
Part 2: psi0 static Zerilli for ALL l: derive the ODE  P2 chi'' + P1 chi' + P0 chi = 0, turn it into a
        coefficient recurrence, and check it identically in (n, l) using the hypergeometric ratio c_{n+1}/c_n.
Part 3: positivity proof inequality (symbolic in l) + exact root isolation for small l.
"""
import sys
import sympy as sp
sys.path.insert(0, '/home/claude/negative-mass/phase2/dg_spectrum')
from common import chi_poly, rs_

l, n, r = sp.symbols('l n r', positive=True)
L = l*(l + 1)
K = (l - 1)*l*(l + 1)*(l + 2)            # = 6k


def cn(nn, ll=l):
    """c_n = (n-2)((n-4)l(l+1)+n-1) Gamma(l+n-1)/(2^n Gamma(n)^2 Gamma(l-n+3)), M=-1."""
    LL = ll*(ll + 1)
    return sp.Integer(nn - 2)*((nn - 4)*LL + nn - 1)*sp.rf(ll - nn + 3, 2*nn - 4) / (2**nn*sp.factorial(nn - 1)**2)


print("== Part 1: coefficients ==")
for k_ in range(1, 8):
    print(f"c_{k_} =", sp.factor(cn(k_)))
for ll in range(2, 13):
    d = sp.expand(sum(cn(k_, sp.Integer(ll))*rs_**k_ for k_ in range(1, ll + 3)) - chi_poly(ll))
    assert d == 0
print("closed form agrees with common.chi_poly for l = 2..12")
print("Sign: for n >= 4, (n-2) > 0, (n-4)l(l+1)+n-1 > 0, and rf(l-n+3, 2n-4) = (l-n+3)...(l+n-2) is a product")
print("of consecutive integers: > 0 if n <= l+2, = 0 if n >= l+3 (contains the factor 0). So c_n >= 0 for n >= 4;")
print("c_1 = 3/2 > 0, c_2 = 0, c_3 = -K/32 < 0 is the ONLY negative coefficient.  deg chi = l+2, c_{l+2} > 0.")

print("\n== Part 2: chi-ODE and recurrence, identically in l ==")
M = -1
lam = (l - 1)*(l + 2)/2
f = 1 - 2*M/r
V = 2*f*(lam**2*r**2*((lam + 1)*r + 3*M) + 9*M**2*(lam*r + M))/(r**3*(lam*r + 3*M)**2)
chi = sp.Function('chi')(r)
psi = chi/(2*lam*r + 6*M)
E = sp.together(-f*sp.diff(f*sp.diff(psi, r), r) + V*psi)
num = sp.numer(E)
num = sp.expand(num)
P2 = sp.factor(num.coeff(sp.diff(chi, r, 2)))
P1 = sp.factor(num.coeff(sp.diff(chi, r)).subs(sp.diff(chi, r, 2), 0))
P0 = sp.factor(sp.expand(num - P2*sp.diff(chi, r, 2) - P1*sp.diff(chi, r)).coeff(chi))
chk = sp.simplify(num - P2*sp.diff(chi, r, 2) - P1*sp.diff(chi, r) - P0*chi)
assert chk == 0
# remove common factor
G = sp.gcd(sp.gcd(P2, P1), P0)
P2, P1, P0 = [sp.factor(sp.cancel(p/G)) for p in (P2, P1, P0)]
print("P2 =", P2); print("P1 =", P1); print("P0 =", P0)
# action on r^n: sum_j a_j(n) r^{n+j}
act = sp.expand((P2*n*(n - 1)*r**(n - 2) + P1*n*r**(n - 1) + P0*r**n)/r**n, power_exp=False)
act = sp.expand(sp.powsimp(act))
shifts = {}
for term in sp.Add.make_args(act):
    pw = sp.degree(term, r) if term.has(r) else 0
    # handle negative powers
    e = sp.Poly(term*r**10, r).degree() - 10
    shifts.setdefault(e, 0)
    shifts[e] += sp.simplify(term/r**e)
shifts = {e: sp.factor(v) for e, v in shifts.items() if sp.simplify(v) != 0}
print("operator on r^n -> sum_j a_j(n) r^(n+j):")
for e in sorted(shifts):
    print(f"  j={e}: a_j(n) =", shifts[e])
# recurrence: coefficient of r^m in L[sum c_n r^n] = sum_j a_j(m-j) c_{m-j} = 0
js = sorted(shifts)
# generic-n check using ratio rho(n) = c_{n+1}/c_n (valid where c_n != 0)
def cgen(nn):
    LL = l*(l + 1)
    return (nn - 2)*((nn - 4)*LL + nn - 1)*sp.gamma(l + nn - 1)/(2**nn*sp.gamma(nn)**2*sp.gamma(l - nn + 3))
m = sp.Symbol('m')
expr = sum(shifts[j].subs(n, m - j)*cgen(m - j) for j in js)
ratio = sp.simplify(sp.gammasimp(expr/cgen(m - js[-1])))
print("recurrence residual / c_(m-jmax), generic m, l:", ratio)
# explicit integer check (covers small n where c_n = 0 and the generic ratio is invalid)
ok = True
for ll in range(2, 41):
    cs = {k_: cn(k_, sp.Integer(ll)) for k_ in range(-5, ll + 10)}
    cs = {k_: (v if 1 <= k_ <= ll + 2 else 0) for k_, v in cs.items()}
    for mm in range(-3, ll + 8):
        s = sum(shifts[j].subs({n: mm - j, l: ll})*cs.get(mm - j, 0) for j in js)
        if sp.simplify(s) != 0:
            ok = False; print("FAIL", ll, mm, s)
print("recurrence satisfied exactly for every coefficient, l = 2..40:", ok)

print("\n== Part 2b: recurrence as a polynomial identity in (l, m) (all l, all m >= 5) ==")
# coefficient of r^m: a0(m) c_m + a1(m-1) c_{m-1} + a2(m-2) c_{m-2} = 0.
# c_n = p(n) rf(l-n+3, 2n-4)/(2^n (n-1)!^2), p(n) = (n-2)((n-4)l(l+1)+n-1).
# Factor out rf(l-m+5, 2m-8) [the c_{m-2} factorial]:
#   rf(l-m+4, 2m-6) = rf(l-m+5,2m-8) (l-m+4)(l+m-3),   rf(l-m+3, 2m-4) = rf(l-m+4,2m-6)(l-m+3)(l+m-2).
mm = sp.Symbol('m')
a0 = shifts[0].subs(n, mm); a1 = shifts[1].subs(n, mm - 1); a2 = shifts[2].subs(n, mm - 2)
p = lambda q: (q - 2)*((q - 4)*L + q - 1)
# divide everything by rf(l-m+5,2m-8)/(2^(m-2) (m-3)!^2)
t2 = p(mm - 2)
t1 = p(mm - 1)*(l - mm + 4)*(l + mm - 3)/(2*(mm - 2)**2)
t0 = p(mm)*(l - mm + 4)*(l + mm - 3)*(l - mm + 3)*(l + mm - 2)/(4*(mm - 2)**2*(mm - 1)**2)
ident = sp.simplify(sp.factor(a0*t0 + a1*t1 + a2*t2))
print("polynomial identity residual (should be 0):", ident)
print("(valid for m >= 4 where (m-3)! and the rf factorisation make sense; m = 1..4 checked below symbolically in l)")
cs = {k_: cn(k_) for k_ in range(1, 7)}
for m_ in range(-1, 5):
    s = sum(shifts[j].subs(n, m_ - j)*(cs.get(m_ - j, 0) if m_ - j >= 1 else 0) for j in js)
    print(f"  m={m_}: residual =", sp.simplify(s))
print("Also: c_n = 0 for n >= l+3 automatically (rf contains the factor 0), so the polynomial chi of degree")
print("l+2 satisfies the full recurrence -> psi0 = chi/(c r - 6) solves the static Zerilli eq. for every l >= 2.")

print("\n== Part 3: chi > 0 on r > 0 for all l ==")
a = -cn(3); b = cn(4); c = cn(5)
print("chi/r = 3/2 - a r^2 + b r^3 + c r^4 + (non-negative terms), a =", sp.factor(a), " b =", sp.factor(b), " c =", sp.factor(c))
print("l = 2: c = 0, chi/r = 3/2 - (3/4) r^2 + (1/4) r^3 exactly; min at r = 2, value 1/2 > 0.")
ll2 = sp.expand(chi_poly(2)/rs_); print("   check:", ll2, " at r=2:", ll2.subs(rs_, 2))
# l >= 3: q(r) = 3/2 - a r^2 + b r^3 + c r^4, q' = r(-2a + 3 b r + 4 c r^2): unique positive critical point r0.
rho = sp.sqrt(9*b**2 + 32*a*c)
r0 = (-3*b + rho)/(8*c)
rhoDG = sp.sqrt(6*l*(l**2 - 1)*(l + 2) - 108)
r0DG = 4*(rhoDG - 6)/((l - 2)*(l + 3)*(l**2 + l + 4))
for lv in [3, 4, 5, 7, 10, 30]:
    print(f"  l={lv}: r0 (ours) = {float(r0.subs(l, lv)):.10f}, r0 (DG eq.58) = {float(r0DG.subs(l, lv)):.10f}")
# rho^2 = 9b^2+32ac; check rho = (K/32) * rhoDG * const ?
print("  rho^2 / rhoDG^2 =", sp.factor(sp.simplify(rho**2/rhoDG**2)))
# q(r0): use q'(r0)=0 -> 4c r0^2 = 2a - 3b r0  to reduce:
# q(r0) = 3/2 - a r0^2 + b r0^3 + c r0^4.  Exact positivity for l>=3: write l = 3 + t, t >= 0.
# Rigorous route: q(r0) > 0  <=>  (resultant-free) evaluate with s = rhoDG symbolic, s^2 = D(l).
s = sp.Symbol('s', positive=True)
Dl = 6*l*(l**2 - 1)*(l + 2) - 108
fac = sp.sqrt(sp.factor(rho**2/rhoDG**2))
r0s = (-3*b + fac*s)/(8*c)
qs = sp.expand(3*sp.Rational(1, 2) - a*r0s**2 + b*r0s**3 + c*r0s**4)
# reduce modulo s^2 = Dl
qs = sp.Poly(qs, s)
red = 0
for (e,), co in qs.terms():
    red += co*Dl**(e//2)*s**(e % 2)
red = sp.expand(red)
A0 = sp.factor(red.coeff(s, 0)); A1 = sp.factor(red.coeff(s, 1))
print("  q(r0) = A0(l) + A1(l)*s,  s = sqrt(6l(l^2-1)(l+2)-108)")
print("  A0 =", A0); print("  A1 =", A1)
Kl = (l - 1)*l*(l + 1)*(l + 2)
den = ((l - 2)*(l + 3)*(l**2 + l + 4))**3
target = 4*(Kl - 18)*(Kl*s - 18*(Kl - 16))/den
print("  q(r0) - 4(K-18)(K s - 18(K-16))/((l-2)(l+3)(l^2+l+4))^3 =", sp.simplify((A0 + A1*s - target).subs(sp.Abs(l - 1), l - 1)))
print("  (K s)^2 - (18(K-16))^2 = 6(K-24)^3 > 0 for K > 24 (l >= 3)  =>  q(r0) > 0  =>  chi/r >= q > 0.")
for lv in [3, 4, 6, 10, 25]:
    qd = (sp.Rational(3, 2) - a*r**2 + b*r**3 + c*r**4).subs(l, lv)
    r0v = float(r0.subs(l, lv))
    print(f"  l={lv}: q(r0) direct = {float(qd.subs(r, r0v)):.12g}, formula = {float(target.subs({l: lv, s: sp.sqrt(Dl.subs(l, lv))})):.12g}")

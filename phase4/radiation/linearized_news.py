"""Independent linearized-gravity computation of the news of a +/- dipole in hyperbolic motion.
Conserved pole-dipole source with zero 4-momentum (MPD: p=0, S const):
  T^{mu nu} = -d_alpha int dtau S^{alpha(mu} u^{nu)} delta^4(x - z(tau)),  S^{ab} = D (n^a u^b - u^a n^b)
(dipole Sum mu x = +D n, positive mass in front). Conservation holds for ANY worldline in the t-z plane since
u^t n^z - u^z n^t = 1 => S = D dt^dz is constant.  Worldline z = a cosh(eta), t = a sinh(eta).
Far field: r hbar^{mu nu} = -4 d_u [ S^{alpha(mu} u^{nu)} k_alpha / (k.u) ]_ret, k = d(t - n.x).
c = r h_+/2 = (1/4) sin^2(th) r hbar^{zz}.  Compare with Bonnor-Swaminarayan news (both pairs).
"""
import sympy as sp
u, a, D = sp.symbols('u a D', positive=True)
U, c0 = sp.symbols('U c0', real=True)
s2 = 1 - c0**2
# retarded rapidity: sinh(eta) - c0 cosh(eta) = U (U = u/a); X = e^eta
X = (U + sp.sqrt(U**2 + s2))/(1 - c0)
sh = (X - 1/X)/2; ch = (X + 1/X)/2
ku = ch - c0*sh      # k.u
F = D*sh*(sh*U/ku - ch)   # S^{alpha z} k_alpha u^z / (k.u)
# c1 = -s2 * dF/du = -(s2/a) dF/dU ; c1_u = -(s2/a^2) d2F/dU2
c1u = -(s2/a**2)*sp.diff(F, U, 2)
c1u_mirror = c1u.subs(c0, -c0)
BS = 3*D*U*s2/(a**2*(U**2 + s2)**sp.Rational(5, 2))   # BS news written with u = aU
import random
print("check sinh-c0 cosh = U:", sp.simplify((sh - c0*ch - U).subs({U: sp.Rational(3,7), c0: sp.Rational(1,3)})))
for _ in range(6):
    vals = {U: sp.Rational(random.randint(-30, 30), 10), c0: sp.Rational(random.randint(-9, 9), 10), a: 1, D: 1}
    one = sp.N(c1u.subs(vals), 15); two = sp.N((c1u + c1u_mirror).subs(vals), 15); bs = sp.N(BS.subs(vals), 15)
    print(f"U={vals[U]}, cos={vals[c0]}: single pair {one}, both pairs {two}, BS {bs}, ratio {sp.N(two/bs,12) if bs!=0 else None}")
# (slow symbolic simplification removed; closed form c_u = -3 D a^2 u s^2/(u^2+a^2 s^2)^{5/2} for one pair
#  is checked to 1e-22 at 200 random points in radiated_momentum.py)

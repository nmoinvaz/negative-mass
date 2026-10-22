"""Check the Bondi mass / news of the Bonnor-Swaminarayan (BS) +/- pair (Pravda & Pravdova, gr-qc/9806114).
Units: G=c=1. Pair: rear/front masses -/+ m, proper separation d, proper acceleration A = 1/a, a^2 = 2h.
Identification (small eps): m = eps^2/(2h)^{3/2}, d = eps/sqrt(2h)  =>  A = m/d^2 exactly (Bondi dipole),
D := m d = eps^3/(2h)^2.  News: c_u = 3 eps^3 u s^2 / (2h (u^2 + 2h s^2)^{5/2}) = 3 D a^2 u s^2/(u^2+a^2 s^2)^{5/2}.
Bondi: m_B = (1/2) int M sin, dm_B/du = -(1/2) int_0^pi c_u^2 sin(th) dth  =: -L(u).
"""
import sympy as sp
u, x, h, eps, D, a = sp.symbols('u x h epsilon D a', positive=True)
th = sp.symbols('theta')

# 1. identification
m1 = eps**2/(2*h*sp.sqrt(2*h+2*eps)); m2 = -eps**2/((2*h+2*eps)*sp.sqrt(2*h))
dprop = sp.sqrt(2*h+2*eps)-sp.sqrt(2*h)
print("m1 series:", sp.series(m1, eps, 0, 4).removeO())
print("m2 series:", sp.series(m2, eps, 0, 4).removeO())
print("(m1+m2) leading:", sp.simplify(sp.series(m1+m2, eps, 0, 4).removeO()))
A = 1/sp.sqrt(2*h)
mm = eps**2/(2*h)**sp.Rational(3,2); dd = eps/sp.sqrt(2*h)
print("A - m/d^2 =", sp.simplify(A - mm/dd**2), "; proper sep series:", sp.series(dprop, eps, 0, 3).removeO())
print("m d - eps^3/(2h)^2 =", sp.simplify(mm*dd - eps**3/(2*h)**2))

# 2. luminosity from the news (a=1, D=1 units; L scales as D^2 A^4 * f(A u))
s2 = 1-x**2
cu = 3*u*s2/(u**2+s2)**sp.Rational(5,2)          # D=1,a=1
Lint = sp.Rational(1,2)*sp.integrate(sp.expand(cu**2), (x, -1, 1))   # sin(th)dth = dx
L = sp.simplify(Lint)
print("L(u) [D=a=1] =", L)
# turnaround check: emission at t_e=0 seen at u=-cos(th): c_u = -3 cos sin^2 exactly
cu_turn = sp.simplify((3*u*s2/(u**2+s2)**sp.Rational(5,2)).subs(u, -x))
print("news at u=-cos th (turnaround emission):", cu_turn)
L0 = sp.Rational(1,2)*sp.integrate(cu_turn**2, (x,-1,1))
print("rest-frame luminosity L0 =", L0, "* D^2 A^4   (Phase 3 octupole: 24/35)")

# 3. Pravda-Pravdova closed form m(u) (h = 1/2 <-> a=1; o = 3 eps^3/(2h) = 3 D a^2 = 3)
hh = sp.Rational(1,2); o = 3
def mPP(U):
    t1 = o**2*U**3*(5*U**4+32*U**2*hh+64*hh**2)*sp.sqrt(2)*(-sp.log(U**2)+sp.log(sp.Abs(-2*sp.sqrt(2*hh)*sp.sqrt(U**2+2*hh)+U**2+4*hh)))/(4096*hh**sp.Rational(7,2)*(U**2+2*hh)**sp.Rational(7,2))
    t2 = o**2*(U**2+4*hh)*(15*U**4+16*U**2*hh+32*hh**2)/(3072*hh**3*U*(U**2+2*hh)**3)
    return t1+t2
mu = mPP(u)
for U0 in [sp.Rational(1,10), 1, 3, 10]:
    lhs = sp.N(sp.diff(mu,u).subs(u,U0), 20); rhs = sp.N(-L.subs(u,U0), 20)
    print(f"u={U0}: dm/du={lhs}, -L={rhs}")
# m(inf) and small-u behaviour
import mpmath as mp
mp.mp.dps = 30
Lf = sp.lambdify(u, L, 'mpmath'); mf = sp.lambdify(u, mu, 'mpmath')
for U0 in [0.01, 0.1, 1, 10, 100]:
    tail = mp.quad(Lf, [U0, 10*U0, mp.inf])
    print(f"u={U0}: m_PP={mp.nstr(mf(U0),12)}  int_u^inf L={mp.nstr(tail,12)}  u*m={mp.nstr(U0*mf(U0),10)}  u^3 m={mp.nstr(U0**3*mf(U0),8)}")
print("small-u coefficient: lim u->0 u*m =", sp.limit(u*mu, u, 0))
print("large-u: series of m in 1/u:", sp.series(mu.subs(u, 1/x), x, 0, 6))
print("m_PP(-u)+m_PP(u) at u=0.7:", sp.N(mPP(-sp.Rational(7,10))+mPP(sp.Rational(7,10))))
print("L even in u:", sp.simplify(L.subs(u,-u)-L) if True else None)

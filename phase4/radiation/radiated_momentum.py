"""Linearized single +/- dipole (D = m d, proper acceleration A = 1/a) in hyperbolic motion z=a cosh eta, t=a sinh eta.
(1) check numerically that its news equals minus the Bonnor-Swaminarayan news of Pravda-Pravdova (gr-qc/9806114 eq.22)
(2) radiated 4-momentum per unit rapidity of the emitter: expect dP^mu/deta = (L0/A) u^mu, L0 = 24/35 D^2 A^4.
Units D = a = 1. Bondi flux: dE/(du dOmega) = c_u^2/(4 pi); dP_z = cos(theta) dE.
"""
import sympy as sp, mpmath as mp, random
U, c0 = sp.symbols('U c0', real=True)
s2 = 1 - c0**2
X = (U + sp.sqrt(U**2 + s2))/(1 - c0)
sh = (X - 1/X)/2; ch = (X + 1/X)/2
F = sh*(sh*U/(ch - c0*sh) - ch)
c1u = -s2*sp.diff(F, U, 2)
BS = 3*U*s2/(U**2 + s2)**sp.Rational(5, 2)
f1 = sp.lambdify((U, c0), c1u, 'mpmath'); fBS = sp.lambdify((U, c0), BS, 'mpmath')
mp.mp.dps = 30
random.seed(1); worst = 0
for _ in range(200):
    Uv = mp.mpf(random.uniform(-5, 5)); cv = mp.mpf(random.uniform(-0.999, 0.999))
    worst = max(worst, abs(f1(Uv, cv) + fBS(Uv, cv))/(abs(fBS(Uv, cv)) + mp.mpf('1e-30')))
print("max rel |single + BS| over 200 random points:", mp.nstr(worst, 5))
# closed form for single pair: c_u = -3 U s^2 /(U^2+s^2)^{5/2}; use it below
def cu(Uv, x): return -3*Uv*(1-x**2)/(Uv**2 + 1 - x**2)**2.5
def uret(eta, x): return mp.sinh(eta) - x*mp.cosh(eta)
def dE(eta):
    return 0.5*mp.quad(lambda x: cu(uret(eta, x), x)**2*(mp.cosh(eta) - x*mp.sinh(eta)), [-1, 0, 1])
def dPz(eta):
    return 0.5*mp.quad(lambda x: x*cu(uret(eta, x), x)**2*(mp.cosh(eta) - x*mp.sinh(eta)), [-1, 0, 1])
L0 = mp.mpf(24)/35
for eta in [0, 0.5, 1, 2, 4]:
    e, p = dE(eta), dPz(eta)
    print(f"eta={eta}: dE/deta={mp.nstr(e,12)} (L0 cosh={mp.nstr(L0*mp.cosh(eta),12)}),  dPz/deta={mp.nstr(p,12)} (L0 sinh={mp.nstr(L0*mp.sinh(eta),12)})")

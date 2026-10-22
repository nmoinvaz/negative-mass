"""
AXIAL (odd-parity) l >= 2 perturbations of the signed-source boson star, Regge-Wheeler gauge:
    h_tph = h0(r) sin(th) Y'(th) e^{st},  h_rph = h1(r) sin(th) Y'(th) e^{st},  delta Phi = 0.
(The scalar has no axial part; its stress tensor enters only through the background.)

Reduction: (th ph) eq gives h0 algebraically in terms of h1, h1'; (r ph) is then a 2nd-order ODE in h1;
(t ph) must hold identically -- we instead build the master variable Q = e^{(nu-lam)/2} h1 / r and derive
    d^2Q/dr*^2 = [ s^2 + V_ax(r) ] Q,   dr*/dr = e^{(lam-nu)/2},
checking all three equations, and compare V_ax with the generic anisotropic-matter form
    V_ax = e^nu [ l(l+1)/r^2 - 6m/r^3 + 4 pi eps (rho - p_r) ]      (vacuum: Regge-Wheeler).
"""
import sympy as sp, os, pickle
from lincommon import *
HERE = os.path.dirname(os.path.abspath(__file__))

h0, h1 = sp.Function('h0')(r), sp.Function('h1')(r)
Yp = sp.diff(Yf, th)
h = sp.zeros(4, 4)
h[0, 3] = h[3, 0] = h0*sp.sin(th)*Yp*ex
h[1, 3] = h[3, 1] = h1*sp.sin(th)*Yp*ex
eqs, bg = linearize(h, sp.Integer(0), sp.Integer(0))

def clean(z):
    z = sp.expand(sp.simplify(z/ex))
    return angsub(bgsub(z))

E = {k: clean(eqs[k]) for k in [(0, 3), (1, 3), (2, 3)]}
# (t ph), (r ph) ~ sin(th) * (Y'' - cot Y') structure -> after angsub: proportional to sin * Y' etc.
for k, v in E.items():
    print(k, "free symbols:", [x for x in v.free_symbols if x in (th, Ysym, Ypsym)])
# extract radial factors
def radial(z):
    z = sp.expand(z)
    # remove the angular factor by dividing by the coefficient structure
    zz = sp.simplify(z)
    fac = sp.factor_terms(zz)
    return fac
R = {}
# Expected structure: (t ph),(r ph): f(r) * sin(th) * Y' ; (th ph): g(r) * sin(th)*(Y'' - cot Y') -> after Legendre:
# Y'' - cot Y' = -2 cot Y' - LL Y.
for k, v in E.items():
    a = sp.simplify(v/(sp.sin(th)*Ypsym))
    if a.has(th) or a.has(Ysym):
        # theta-phi: v = g * sin * (-2 cot Y' - LL Y)  =>  coefficient of Ysym is -LL g sin
        a = sp.simplify(sp.expand(v).coeff(Ysym)/(-LL*sp.sin(th)))
        chk = sp.simplify(v - a*sp.sin(th)*(-2*sp.cos(th)/sp.sin(th)*Ypsym - LL*Ysym))
        print(k, "angular check:", chk)
    R[k] = sp.simplify(a)
    print(k, "radial:", R[k])

# (th ph) is algebraic in h0 (contains h1, h1'); (r ph) contains h0' -> 2nd-order eq. in h1
th_ph = R[(2, 3)]
r_ph = R[(1, 3)]
h0sol = sp.simplify(sp.solve(th_ph, h0)[0])
print("\n(th ph) gives h0 =", h0sol)
eq_h1 = sp.simplify(bgsub((r_ph*2*r**2*sp.exp(nu0)).subs(sp.diff(h0, r), sp.diff(h0sol, r)).subs(h0, h0sol).doit()))
eq_h1 = sp.simplify(bgsub(eq_h1))
# master variable Q = e^{(nu-lam)/2} h1 / r
Q = sp.Function('Q')(r)
h1Q = r*sp.exp((lam0 - nu0)/2)*Q
eqQ = sp.simplify(bgsub(eq_h1.subs(h1, h1Q).doit()))
# write as Q'' + c1 Q' + c0 Q = 0
eqQ = sp.expand(eqQ)
c2 = sp.simplify(eqQ.coeff(sp.diff(Q, r, 2)))
c1 = sp.simplify(eqQ.coeff(sp.diff(Q, r)) / c2)
c0 = sp.simplify((eqQ - c2*sp.diff(Q, r, 2) - eqQ.coeff(sp.diff(Q, r))*sp.diff(Q, r)).coeff(Q) / c2)
# tortoise: d^2/dr*^2 = e^{nu-lam} [ d^2/dr^2 + (nu'-lam')/2 d/dr ]
print("\nQ-equation first-derivative coefficient minus tortoise form (0?):", sp.simplify(bgsub(c1 - (sp.diff(nu0, r) - sp.diff(lam0, r))/2)))
Vtot = sp.simplify(-c0*sp.exp(nu0 - lam0))       # Q_{r*r*} = (s^2 + V) Q  -> Q'' + c1 Q' - e^{lam-nu}(s^2+V) Q = 0
Vax = sp.simplify(bgsub(Vtot - s**2))
print("V_ax has s?", Vax.has(s))
# compare with generic form
a2 = sp.exp(lam0); al2 = sp.exp(nu0)
m = r*(1 - 1/a2)/2
rho = (w**2*p0**2/al2 + sp.diff(p0, r)**2/a2)/2 + Vf(p0**2)
pr = (w**2*p0**2/al2 + sp.diff(p0, r)**2/a2)/2 - Vf(p0**2)
Vgen = al2*(LL/r**2 - 6*m/r**3 + 4*sp.pi*eps*(rho - pr))
print("V_ax - e^nu[LL/r^2 - 6m/r^3 + 4 pi eps (rho - p_r)] =", sp.simplify(bgsub(Vax - Vgen)))
# remaining equation (t ph) on solutions
t_ph = R[(0, 3)]
res = t_ph.subs(sp.diff(h0, r, 2), sp.diff(h0sol, r, 2)).subs(sp.diff(h0, r), sp.diff(h0sol, r)).subs(h0, h0sol).doit()
res = bgsub(res)
h1pp = sp.solve(eq_h1, sp.diff(h1, r, 2))[0]
for _ in range(2):
    res = bgsub(res.subs(sp.diff(h1, r, 3), sp.diff(h1pp, r)).subs(sp.diff(h1, r, 2), h1pp).doit())

print("(t ph) residual on solutions:", sp.simplify(res))
print("\nh0 =", h0sol)
pickle.dump({'Vax': Vax, 'h0': h0sol}, open(os.path.join(HERE, 'axial_eqs.pkl'), 'wb'))

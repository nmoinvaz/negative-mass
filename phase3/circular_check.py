"""
Literature check of the 1PN acceleration used everywhere (Blanchet eq. 203 form), including the m1 m2 term
that the test-body check (test_body.py) cannot see: circular orbits in harmonic coordinates must satisfy
    omega^2 = (M/r^3) [1 + (nu - 3) M/r],     M = m1 + m2, nu = m1 m2 / M^2
and the EIH energy must give the gauge-invariant  E(x) = -(M nu x / 2) [1 - (3/4 + nu/12) x],  x = (M omega)^(2/3).
(Blanchet, Living Rev. Relativ. 17, 2 (2014), eqs. 228 and 233 at 1PN; standard results.)
"""
import sympy as sp

eps = sp.Symbol('epsilon', positive=True)
m1, m2 = sp.symbols('m1 m2', real=True, nonzero=True)
r, w = sp.symbols('r omega', positive=True)
M = m1 + m2; nu = m1 * m2 / M**2

def acc(ma, mb, nab, va, vb):
    return (-mb * nab / r**2
            + eps * ((5 * ma * mb / r**3 + 4 * mb**2 / r**3
                      + mb / r**2 * (sp.Rational(3, 2) * nab.dot(vb)**2 - va.dot(va) + 4 * va.dot(vb) - 2 * vb.dot(vb))) * nab
                     + mb / r**2 * (4 * nab.dot(va) - 3 * nab.dot(vb)) * (va - vb)))

# CM frame circular orbit at t = 0: x1 = (m2/M) r ex, x2 = -(m1/M) r ex (CM correction is O(eps (m1-m2)) and
# cancels in the relative acceleration for circular orbits at 1PN), v_A = omega ey x x_A
n = sp.Matrix([1, 0, 0])
v1 = sp.Matrix([0, w * r * m2 / M, 0]); v2 = sp.Matrix([0, -w * r * m1 / M, 0])
arel = acc(m1, m2, n, v1, v2) - acc(m2, m1, -n, v2, v1)
cond = sp.expand(arel[0] + w**2 * r)            # centripetal: a_rel = -omega^2 r n
w2 = sp.symbols('W2', positive=True)
cond = cond.subs(w, sp.sqrt(w2))
sol = sp.solve(sp.Eq(cond, 0), w2)[0]
sol1 = sp.series(sol, eps, 0, 2).removeO()
print("omega^2 r^3/M - [1 + (nu - 3) M/r] =", sp.simplify(sol1 * r**3 / M - (1 + eps * (nu - 3) * M / r)))
print("tangential relative acceleration (must vanish):", sp.simplify(arel[1]))

# energy
E = (m1 * v1.dot(v1) / 2 + m2 * v2.dot(v2) / 2 - m1 * m2 / r
     + eps * (sp.Rational(3, 8) * (m1 * v1.dot(v1)**2 + m2 * v2.dot(v2)**2)
              + m1 * m2 / (2 * r) * (3 * (v1.dot(v1) + v2.dot(v2)) - 7 * v1.dot(v2) - n.dot(v1) * n.dot(v2))
              + m1 * m2 * M / (2 * r**2)))
# E = sum v.p - L for the LL Lagrangian (derived symbolically in eih.py); written out here.
x = sp.Symbol('x', positive=True)
Ew = sp.expand(E.subs(w, sp.sqrt(sol1)))
# invert x = (M omega)^(2/3) -> r(x) at 1PN:  M/r = x [1 + eps (1 - nu/3) x]  (from omega^2 relation)
rx = M / (x * (1 + eps * (1 - nu / 3) * x))
chk = sp.series(sp.expand(sol1.subs(r, rx)) * M**2 - x**3, eps, 0, 2).removeO()
print("check r(x) inversion (should be 0):", sp.simplify(chk))
target = -(M * nu * x / 2) * (1 - eps * (sp.Rational(3, 4) + nu / 12) * x)
for mm1, mm2 in ((1, 1), (sp.Rational(3, 7), sp.Rational(5, 2)), (2, -sp.Rational(1, 3))):
    sub = {m1: mm1, m2: mm2}
    Ex = sp.series(sp.together(Ew.subs(sub).subs(r, rx.subs(sub))), eps, 0, 2).removeO()
    print(f"  m1={mm1}, m2={mm2}: E(x) - [-(M nu x/2)(1 - (3/4 + nu/12) x)] =", sp.simplify(Ex - target.subs(sub)))

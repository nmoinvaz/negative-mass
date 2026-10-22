"""
Collinear signed-mass pair at 1PN (fast; uses Blanchet eq. 203, which eih.py verifies equals the
Euler-Lagrange equations of the EIH Lagrangian with m_A -> mu_A).

Body 1 (mu1 = m1 > 0) at x = d (front), body 2 (mu2 = -m2 < 0) at x = 0 (rear). G = c = 1, eps = 1/c^2.
Results:
  * accelerations at rest and with velocities V1, V2 along the axis;
  * the proper-separation acceleration relative to Born-rigid (common-centre hyperbolic) motion,
        delta'' := a1 - a2 + a_mid^2 d   (at rest),
    compared with the total 1PN mass-energy E_tot = mu1 + mu2 - mu1 mu2/d (rest masses + binding);
  * transverse common velocity.
"""
import sympy as sp

eps = sp.Symbol('epsilon', positive=True)
d, m, m1, m2 = sp.symbols('d m m1 m2', positive=True)
V, V1, V2 = sp.symbols('V V1 V2', real=True)

def acc(ma, mb, nab, va, vb, r):
    """Blanchet (2014) eq. (203), 1PN harmonic gauge, 1D (n, v scalars along the axis)."""
    return (-mb * nab / r**2
            + eps * ((5 * ma * mb / r**3 + 4 * mb**2 / r**3
                      + mb / r**2 * (sp.Rational(3, 2) * (nab * vb)**2 - va**2 + 4 * va * vb - 2 * vb**2)) * nab
                     + mb / r**2 * (4 * nab * va - 3 * nab * vb) * (va - vb)))

mu1, mu2 = m1, -m2
a1 = sp.expand(acc(mu1, mu2, 1, V1, V2, d))
a2 = sp.expand(acc(mu2, mu1, -1, V2, V1, d))
print("a1 (front, mu1=+m1) =", sp.collect(a1, eps))
print("a2 (rear,  mu2=-m2) =", sp.collect(a2, eps))

rest = {V1: 0, V2: 0}
a1r, a2r = a1.subs(rest), a2.subs(rest)
print("\nat rest:  a1 =", sp.factor(a1r), "   a2 =", sp.factor(a2r))
dip = {m1: m, m2: m}
print("dipole at rest: a1 =", sp.factor(a1r.subs(dip)), ", a2 =", sp.factor(a2r.subs(dip)),
      ", a_mid =", sp.factor(((a1r + a2r) / 2).subs(dip)), ", d'' =", sp.factor((a1r - a2r).subs(dip)))

amid = (a1r + a2r) / 2
delta_dd = sp.expand(a1r - a2r + eps * (amid**2 * d).subs(eps, 0))   # Born-rigid term is O(eps)
Etot = mu1 + mu2 + eps * (-mu1 * mu2 / d)                               # rest-mass + binding (c^2 = 1/eps)
# Rest-mass difference m1 - m2 counts as O(eps) (mismatch ~ m^2/d); compare to -E_tot/d^2:
print("\nproper-separation acceleration at rest  a1 - a2 + a_mid^2 d =", sp.factor(delta_dd))
print("  -E_tot/d^2                                                 =", sp.factor(-Etot / d**2))
diff = sp.expand(delta_dd + Etot / d**2)
print("  difference =", sp.factor(diff), " -> zero at O(eps) when m2 - m1 = O(eps):",
      sp.simplify(diff.subs(m2, m1 * (1 + sp.Symbol('k') * eps * m1 / d)).series(eps, 0, 2).removeO()))
ks = sp.solve(sp.Eq(m1 - m2 + m1 * m2 / d, 0), m2)
print("  rigid (E_tot = 0) mass condition: m2 =", ks, " = m1 (1 + m1/d) + O(m1^3/d^2)")

# common velocity V along the axis (dipole)
cv = {V1: V, V2: V}
a1V, a2V = sp.factor(a1.subs(cv).subs(dip)), sp.factor(a2.subs(cv).subs(dip))
print("\ndipole, common axial velocity V: a1 =", a1V, ", a2 =", a2V, ", a_mid =", sp.factor((a1V + a2V) / 2))
# relative velocity u = V1 - V2 (small) around common V
u = sp.Symbol('u')
rel = sp.expand((a1 - a2).subs({V1: V + u / 2, V2: V - u / 2}).subs(dip))
print("dipole, d'' with V1 = V + u/2, V2 = V - u/2:", sp.collect(sp.factor_terms(rel), eps))

# transverse common velocity: need 3D; use Blanchet vector form along x with v = (0, V, 0)
n = sp.Matrix([1, 0, 0]); vv = sp.Matrix([0, V, 0])
def acc3(ma, mb, nab, va, vb, r):
    return (-mb * nab / r**2
            + eps * ((5 * ma * mb / r**3 + 4 * mb**2 / r**3
                      + mb / r**2 * (sp.Rational(3, 2) * nab.dot(vb)**2 - va.dot(va) + 4 * va.dot(vb) - 2 * vb.dot(vb))) * nab
                     + mb / r**2 * (4 * nab.dot(va) - 3 * nab.dot(vb)) * (va - vb)))
t1 = sp.simplify(acc3(m, -m, n, vv, vv, d)); t2 = sp.simplify(acc3(-m, m, -n, vv, vv, d))
print("dipole, common transverse velocity V: a1 =", list(t1), ", a2 =", list(t2))

# Shirokov campaign numbers
num = {m1: 0.014350, m2: 0.014295, d: 10, eps: 1}
print("\nShirokov's quoted masses (M+ = 0.014350, M- = -0.014295), d = 10, at rest:")
print(f"  a1 d^2/Mbar = {float((a1r * d**2).subs(num) / 0.0143225):.5f},  a2 d^2/Mbar = {float((a2r * d**2).subs(num) / 0.0143225):.5f}")
print(f"  E_tot = {float(Etot.subs(num)):.3e};  proper-separation accel = {float(delta_dd.subs(num)):.3e}"
      f"  -> delta(t) = {float(delta_dd.subs(num)) / 2:.3e} t^2 :  t=200: {float(delta_dd.subs(num)) / 2 * 4e4:+.4f},"
      f" t=400: {float(delta_dd.subs(num)) / 2 * 16e4:+.4f}")
num2 = {m1: 0.014626, m2: 0.014580, d: 10, eps: 1}
print(f"  (Phase-2 masses +0.014626/-0.014580: E_tot = {float(Etot.subs(num2)):.3e},"
      f" delta(400) = {float(delta_dd.subs(num2)) / 2 * 16e4:+.4f})")

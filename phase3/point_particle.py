"""
Task 1 - which "masses" are negative for the signed point-particle action

    S = sum_A s_A ( -m_A int dtau_A )  +  sum_A q_A int A_mu dx^mu  +  S_EH  +  S_Maxwell,
    m_A > 0, s_A = +-1,   mu_A := s_A m_A.

Checks (sympy):
  (a) Euler-Lagrange equations of the worldline action in a GENERAL metric with an EM potential:
          mu (u^nu nabla_nu u^mu) = q F^mu_nu u^nu.
      For q = 0 mu drops out (geodesic, any sign). We verify this in a general 1+1 metric
      (t, x) with arbitrary g_ab(t, x) and A_a(t, x) using the proper-time parametrisation.
  (b) Source of the metric: T^{mu nu}(x) = mu int u^mu u^nu delta^4(x - z(tau)) / sqrt(-g) dtau,
      i.e. the action is linear in mu, so the metric source scales with mu (shown analytically in
      RESULTS.md; here we check the Newtonian limit of the field equation sources rho = mu delta).
  (c) Flat-space kinetic term and canonical momentum: L = -mu sqrt(1 - v^2) -> p = mu gamma v,
      E = mu gamma.
  (d) Response to a non-gravitational force (flat space, standard-sign q int A):
          mu d(gamma v)/dt = q (E + v x B)  ->  for mu < 0 the body accelerates AGAINST the force.
      If instead the charge coupling is inside the flipped sector (s_A q_A int A), then
          m d(gamma v)/dt = q (E + v x B)  but it SOURCES Maxwell with charge s_A q_A:
      identical to a Bondi negative-mass body (mass -m) with charge -q.
"""
import sympy as sp

tau = sp.Symbol('tau')
mu, q = sp.symbols('mu q', nonzero=True)
T, X = sp.Function('T')(tau), sp.Function('X')(tau)
gtt, gtx, gxx, At, Ax = [sp.Function(n) for n in ('g_tt', 'g_tx', 'g_xx', 'A_t', 'A_x')]
t, x = sp.symbols('t x')

def ev(f):  # field evaluated on the worldline
    return f(T, X)

# (a) Worldline Lagrangian with affine parameter; use L = -mu sqrt(-g u u) + q A u
Td, Xd = T.diff(tau), X.diff(tau)
guu = ev(gtt) * Td**2 + 2 * ev(gtx) * Td * Xd + ev(gxx) * Xd**2
Lw = -mu * sp.sqrt(-guu) + q * (ev(At) * Td + ev(Ax) * Xd)

def EL(L, Q):
    return sp.diff(L.diff(Q.diff(tau)), tau) - L.diff(Q)

eqs = [EL(Lw, T), EL(Lw, X)]
# impose proper-time parametrisation AFTER differentiating: -g u u = 1
# Write each EL as  mu * G_a + q * F_a  and verify G_a = -g_ab (u.nabla u)^b, F_a = -F_ab u^b
# (signs: d/dtau(dL/du^a) - dL/dx^a = 0  =>  -mu (g_ab a^b) ... = q F_ab u^b)
gm = sp.Matrix([[gtt(t, x), gtx(t, x)], [gtx(t, x), gxx(t, x)]])
coords = (t, x)
ginv = gm.inv()
Gam = [[[sum(ginv[a, d] * (gm[d, b].diff(coords[c]) + gm[d, c].diff(coords[b]) - gm[b, c].diff(coords[d]))
             for d in range(2)) / 2 for c in range(2)] for b in range(2)] for a in range(2)]
u = [Td, Xd]
acc = [T.diff(tau, 2), X.diff(tau, 2)]
on_wl = {t: T, x: X}
geo = [acc[a] + sum(Gam[a][b][c].subs(on_wl) * u[b] * u[c] for b in range(2) for c in range(2)) for a in range(2)]
Fab = sp.Matrix([[0, Ax(t, x).diff(t) - At(t, x).diff(x)], [At(t, x).diff(x) - Ax(t, x).diff(t), 0]])
# F_ab = d_a A_b - d_b A_a :  F_tx = d_t A_x - d_x A_t
lower_geo = [sum(gm[a, b].subs(on_wl) * geo[b] for b in range(2)) for a in range(2)]
lorentz = [sum(Fab[a, b].subs(on_wl) * u[b] for b in range(2)) for a in range(2)]

norm_sub = {}
ok = True
for a in range(2):
    # expected: EL_a = mu * g_ab a^b (with sqrt(-guu) = 1 after differentiation)  - q F_ab u^b
    expr = eqs[a]
    expr = expr.subs(sp.sqrt(-guu), 1)
    # remaining derivatives of sqrt(-guu) are proportional to d/dtau(-guu) = 0 in proper time
    dnorm = sp.diff(-guu, tau)
    expr = sp.expand(expr)
    target = sp.expand(mu * lower_geo[a] - q * lorentz[a])
    diff = sp.simplify(sp.expand((expr - target).doit()))
    # eliminate dnorm terms: substitute the second derivative of T from d(-guu)/dtau = 0
    solT = sp.solve(sp.Eq(dnorm, 0), T.diff(tau, 2))
    diff2 = sp.simplify(diff.subs(T.diff(tau, 2), solT[0]))
    print(f"(a) EL_{coords[a]} - [mu g_ab a^b - q F_ab u^b]  (proper time) = {diff2}")
    ok &= diff2 == 0
print("    => mu * (covariant acceleration) = q F u :  mu cancels when q = 0 (geodesic for either sign)." if ok
      else "    MISMATCH")

# (c) flat kinetic term, canonical momentum and energy
v = sp.Symbol('v', real=True)
L0 = -mu * sp.sqrt(1 - v**2)
p = sp.diff(L0, v)
E = sp.simplify(p * v - L0)
print("(c) L =", sp.series(L0, v, 0, 5).removeO(), ";  p =", sp.simplify(p), ";  E =", E)

# (d) flat space Lorentz force with standard-sign coupling: mu d(gamma v)/dt = q E
Ef = sp.Symbol('E_field')
vt = sp.Function('v')(sp.Symbol('t'))
eom = sp.Eq(sp.diff(mu * vt / sp.sqrt(1 - vt**2), sp.Symbol('t')), q * Ef)
dv = sp.solve(eom, sp.diff(vt, sp.Symbol('t')))[0]
print("(d) dv/dt =", sp.simplify(dv), " -> sign(dv/dt) = sign(q E / mu): opposite to the force when mu < 0")
# flipped-sector coupling: s*(-m dtau + q A) with s = -1  ==  (-m) dtau  + (-q) A
s, m = sp.symbols('s m', nonzero=True)
print("    flipped-sector coupling s(-m dtau + q A): EOM coefficient ratio q_eff/mu_eff =",
      sp.simplify((s * q) / (s * m)), "(unchanged response), sources Maxwell with charge", s * q)

# (b) Newtonian limit of the source: the action's only mass parameter is mu, the metric equation
#     G_00 = 8 pi T_00 gives nabla^2 U = -4 pi mu delta  (U = -Phi). Potential of a mu<0 body is repulsive.
r = sp.Symbol('r', positive=True)
U = mu / r
lap = sp.simplify(sp.diff(r**2 * sp.diff(U, r), r) / r**2)
print("(b) U = mu/r is harmonic away from 0:", lap, "; flux -oint grad U . dS / 4pi =",
      sp.simplify(-4 * sp.pi * r**2 * sp.diff(U, r) / (4 * sp.pi)), "= active mass")

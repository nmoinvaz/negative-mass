"""
Independent check of the 1PN acceleration (Blanchet eq. 203) in the test-body limit mu1 -> 0, derived from
first principles: geodesic of the Schwarzschild metric of mass mu (either sign) in harmonic coordinates,
expanded to 1PN,
    g_00 = -(1 - 2 mu/r + 2 mu^2/r^2),  g_ij = (1 + 2 mu/r) delta_ij,  g_0i = 0.
Coordinate-time geodesic:  a^i = -Gamma^i_00 - 2 Gamma^i_0j v^j - Gamma^i_jk v^j v^k
                                 + v^i (Gamma^0_00 + 2 Gamma^0_0j v^j + Gamma^0_jk v^j v^k).
Expected (eq. 203 with m1 = 0, v2 = 0):  a = -mu n/r^2 (1 - 4 mu/r + v^2) + 4 mu (n.v) v / r^2.
"""
import sympy as sp

eps = sp.Symbol('epsilon')
mu = sp.Symbol('mu', real=True, nonzero=True)
t = sp.Symbol('t'); xs = sp.symbols('x y z', real=True); vs = sp.symbols('vx vy vz', real=True)
r = sp.sqrt(sum(c**2 for c in xs))
coords = (t,) + xs
g = sp.diag(-(1 - 2 * eps * mu / r + 2 * eps**2 * mu**2 / r**2), *(3 * [1 + 2 * eps * mu / r]))
# eps counts 1/c^2 with v ~ eps^(1/2): write g_00 with c^2 dt^2, so time derivative terms carry eps
ginv = g.inv()
def Gam(a, b, c):
    return sum(ginv[a, d] * (g[d, b].diff(coords[c]) + g[d, c].diff(coords[b]) - g[b, c].diff(coords[d]))
               for d in range(4)) / 2
U = [1] + list(vs)
acc = []
for i in range(1, 4):
    ai = -sum(Gam(i, b, c) * U[b] * U[c] for b in range(4) for c in range(4)) \
         + vs[i - 1] * sum(Gam(0, b, c) * U[b] * U[c] for b in range(4) for c in range(4))
    acc.append(ai)
# PN counting: Newtonian acceleration is O(eps) here (mu always comes with eps); keep O(eps^2) with v^2 ~ eps
lam = sp.Symbol('lam')
scal = {vs[k]: lam * vs[k] for k in range(3)}
res = []
nv = sum(x * v for x, v in zip(xs, vs)) / r
v2 = sum(v * v for v in vs)
for i in range(3):
    e = sp.series(acc[i].subs(scal).subs(eps, lam**2), lam, 0, 5).removeO().subs(lam, 1)
    expected = -mu * xs[i] / r**3 * (1 - 4 * mu / r + v2) + 4 * mu * nv * vs[i] / r**2
    res.append(sp.simplify(e - expected))
print("geodesic (harmonic Schwarzschild, 1PN) minus Blanchet test-body limit:", res)

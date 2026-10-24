"""Third order of the matched expansion (see matched.py): Omega1 in W^2 = s - 4 ln 2 + Omega1/s + O(1/s^2),
s = 2(D - 3). Closed forms (u = sqrt(R - 1), l = ln R, K(u) = int_0^u ln(1+t^2)/(1+t^2) dt):
  Y0 = 1, L0 = -1, Y1 = -l, L1 = -2 Omega0 - 2 l, Y2 = -2 Omega0 l - (3/2) l^2,
  L2 = -2 Omega1 + 2 Omega0 l + 4 l^2 - 8 Omega0 u arctan u - 16 u K(u),   Omega0 = -4 ln 2.
Then Y3 = int_1^R (L2 + L1 l + L0 l^2/2)/R' dR', h3 from the series, and G3(Omega1) = 0. Everything is done on one
grid in v = sqrt(ln R) (smooth at the throat), with cumulative Simpson integrals, at three grid sizes.

    python3 order3.py
"""
import numpy as np
import sympy as sp
from scipy.integrate import cumulative_simpson, simpson

R = sp.Symbol("R", positive=True)
eps, Om0, Om1, Om2 = sp.symbols("epsilon Omega0 Omega1 Omega2")
Ys = sp.symbols("Y0:4")
Ls = sp.symbols("L0:4")
lg = sp.log(R)
r = sp.exp(eps * lg)
Y = sum(Ys[k] * eps**k for k in range(4))
L = sum(Ls[k] * eps**k for k in range(4))
Om = Om0 + eps * Om1 + eps**2 * Om2
rhs = ((1 + eps * Om) * r * Y + L / (2 * R) + (sp.Rational(1, 2) + eps) * L * (R - 1) / R
       - (sp.Rational(1, 2) + eps) * Y * (R - 1) / (R * r) - Y / (2 * r * R))
h3 = sp.expand(sp.series(rhs, eps, 0, 4).removeO().coeff(eps, 3) - Ls[3] / 2)
assert sp.diff(h3, Ls[3]) == 0
h3f = sp.lambdify((R, Om0, Om1, Om2, *Ys[:4], *Ls[:3]), h3, "numpy")
O0 = -4 * np.log(2)


def G3(o1, V=8.0, N=400001):
    v = np.linspace(0, V, N)
    t = v * v                                   # t = ln R
    Rv = np.exp(t)
    em1 = np.expm1(t)
    u = np.sqrt(em1)
    dK = np.zeros_like(v)
    dK[1:] = v[1:] ** 3 / u[1:]                 # dK/dv = v^3 / sqrt(e^{v^2} - 1)
    K = cumulative_simpson(dK, x=v, initial=0)
    l = t
    L2 = -2 * o1 + 2 * O0 * l + 4 * l * l - 8 * O0 * u * np.arctan(u) - 16 * u * K
    L1 = -2 * O0 - 2 * l
    f = L2 + L1 * l - l * l / 2                 # R Y3'
    Y3 = cumulative_simpson(f * 2 * v, x=v, initial=0)            # dt = 2 v dv
    Y2 = -2 * O0 * l - 1.5 * l * l
    h = h3f(Rv, O0, o1, 0.0, 1.0, -l, Y2, Y3, -1.0, L1, L2)
    g = np.zeros_like(v)
    g[1:] = (h[1:] - h[0]) / em1[1:] ** 1.5 * Rv[1:] * 2 * v[1:]
    g[0] = 2 * g[1] - g[2]                      # g is smooth in v; its v = 0 value by extrapolation
    return simpson(g, x=v)


if __name__ == "__main__":
    for V, N in ((7.0, 200001), (8.0, 400001), (8.5, 800001)):
        g0, g1 = G3(0.0, V, N), G3(1.0, V, N)
        print(f"V = {V}, N = {N}: G3(0) = {g0:.12f}, dG3/dOmega1 = {g1 - g0:.12f}, Omega1 = {-g0 / (g1 - g0):.10f}",
              flush=True)
    print("expected dG3/dOmega1 = -4 pi =", -4 * np.pi)

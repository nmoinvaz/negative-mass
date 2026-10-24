"""Radial perturbations of the massless symmetric ghost-scalar (Ellis) wormhole in D = n + 2 dimensions.

Background (throat radius a = 1): ds^2 = -dt^2 + dl^2 + r(l)^2 dOmega_n^2, r'^2 = 1 - r^{-2(n-1)}, a ghost scalar
with kappa phi'^2 = n r''/r (Einstein equations in Ricci form, R_ab = -kappa d_a phi d_b phi for a ghost).

Warped-product Einstein equations for g = h_AB dx^A dx^B + r^2 dOmega_n (h two-dimensional):
  (E1) (1/2) R[h] h_AB - (n/r) nabla_A nabla_B r = -kappa d_A phi d_B phi
  (E2) (n - 1)(1 - (nabla r)^2) - r box_h r = 0
  (S)  d_A( sqrt(-h) r^n h^{AB} d_B phi ) = 0
Perturbation: h = diag(-e^{2 nu}, e^{2 lam}), nu = eps N(l) e^{W t}, lam = eps L(l) e^{W t}, r = r0 + eps Y(l) e^{W t},
in the gauge delta phi = 0 (possible everywhere since phi' never vanishes). Prints the first-order equations and the
master equation for Y.
"""
import sympy as sp

t, l, eps, W = sp.symbols("t l epsilon W")
n = sp.Symbol("n", positive=True)
r0 = sp.Function("r0")(l)
N, L, Y = (sp.Function(s)(l) for s in ("N", "L", "Y"))
E = sp.exp(W * t)
nu, lam, r = eps * N * E, eps * L * E, r0 + eps * Y * E
x = [t, l]
h = sp.diag(-sp.exp(2 * nu), sp.exp(2 * lam))
hi = h.inv()


def christoffel(a, b, c):
    return sum(hi[a, d] * (sp.diff(h[d, b], x[c]) + sp.diff(h[d, c], x[b]) - sp.diff(h[b, c], x[d])) for d in range(2)) / 2


Gm = [[[christoffel(a, b, c) for c in range(2)] for b in range(2)] for a in range(2)]


def ricci(b, c):
    return sum(sp.diff(Gm[a][b][c], x[a]) - sp.diff(Gm[a][b][a], x[c])
               + sum(Gm[a][a][d] * Gm[d][b][c] - Gm[a][c][d] * Gm[d][b][a] for d in range(2)) for a in range(2))


Rh = sum(hi[a, b] * ricci(a, b) for a in range(2) for b in range(2))


def hess(A, B, f):
    return sp.diff(f, x[A], x[B]) - sum(Gm[C][A][B] * sp.diff(f, x[C]) for C in range(2))


grad2 = sum(hi[a, b] * sp.diff(r, x[a]) * sp.diff(r, x[b]) for a in range(2) for b in range(2))
box = sum(hi[a, b] * hess(a, b, r) for a in range(2) for b in range(2))
kphi2 = n * sp.diff(r0, l, 2) / r0          # kappa phi0'^2, unperturbed (delta phi = 0)
E1 = {(A, B): Rh * h[A, B] / 2 - (n / r) * hess(A, B, r) + (kphi2 if (A, B) == (1, 1) else 0)
      for A in range(2) for B in range(A, 2)}
E2 = (n - 1) * (1 - grad2) - r * box
sqrtmh = sp.exp(nu + lam)
S = sp.diff(sqrtmh * r**n * hi[1, 1] * sp.sqrt(kphi2), l)   # phi' = sqrt(kappa phi0'^2) up to kappa, delta phi = 0


def first(expr):
    return sp.simplify(sp.diff(expr, eps).subs(eps, 0) / E)


if __name__ == "__main__":
    out = {k: first(v) for k, v in E1.items()}
    out["E2"] = first(E2)
    out["S"] = first(S)
    for k, v in out.items():
        print(k, "=", v)

"""Radial perturbations of static wormholes supported by a ghost scalar, with an optional radial electric field, in
D = n + 2 dimensions. Generalizes ../derive.py (which had g_tt = -1 and no charge).

Static background: ds^2 = -e^{2 nu0(x)} dt^2 + e^{2 lam0(x)} dx^2 + r0(x)^2 dOmega_n^2, a ghost scalar phi0(x) and a
radial electric field with charge q, F_tx = q e^{nu + lam} / r^n (Maxwell's equations; q is conserved by radial
perturbations). Units: G_ab = kappa T_ab with kappa = 1 absorbed into phi and q:
  ghost:   T_ab = -(d_a phi d_b phi - 1/2 g_ab (d phi)^2)
  Maxwell: T_AB = -(q^2 / 2 r^{2n}) h_AB on the (t, x) plane, T_ij = +(q^2 / 2 r^{2n}) r^2 gamma_ij on the sphere.
Warped-product Einstein tensor (h two-dimensional, so its own Einstein tensor vanishes):
  G_AB = -(n/r) nabla_A nabla_B r + h_AB [ (n/r) box r + n(n-1)/(2 r^2) ((grad r)^2 - 1) ]
  G_ij = r^2 gamma_ij [ -R[h]/2 + ((n-1)/r) box r + (n-1)(n-2)/(2 r^2) ((grad r)^2 - 1) ]
(checked against a direct computation for n = 2 and 3 by `python3 derive_general.py check`).

Perturbations ~ e^{W t} in the gauge delta phi = 0: nu = nu0 + eps N, lam = lam0 + eps L, r = r0 + eps Y. Then
  (tx)     L = (Y' - nu0' Y) / r0'
  (scalar) N = L - n Y / r0            (a function of t is removed by reparametrizing t)
and the angular equation is the master equation, printed by `python3 derive_general.py master`.
"""
import sys

import sympy as sp

t, x, eps, W = sp.symbols("t x epsilon W")
n = sp.Symbol("n", positive=True)
q = sp.Symbol("q", real=True)
nu0, lam0, r0 = (sp.Function(s)(x) for s in ("nu0", "lam0", "r0"))
N, L, Y = (sp.Function(s)(x) for s in ("N", "L", "Y"))
P = sp.Function("P")(x)          # phi0'


def setup(perturbed=True):
    E = sp.exp(W * t)
    e = eps if perturbed else 0
    nu, lam, r = nu0 + e * N * E, lam0 + e * L * E, r0 + e * Y * E
    X = [t, x]
    h = sp.diag(-sp.exp(2 * nu), sp.exp(2 * lam))
    hi = h.inv()
    Gm = [[[sum(hi[a, d] * (sp.diff(h[d, b], X[c]) + sp.diff(h[d, c], X[b]) - sp.diff(h[b, c], X[d]))
                for d in range(2)) / 2 for c in range(2)] for b in range(2)] for a in range(2)]

    def ric(b, c):
        return sum(sp.diff(Gm[a][b][c], X[a]) - sp.diff(Gm[a][b][a], X[c])
                   + sum(Gm[a][a][d] * Gm[d][b][c] - Gm[a][c][d] * Gm[d][b][a] for d in range(2)) for a in range(2))

    Rh = sum(hi[a, b] * ric(a, b) for a in range(2) for b in range(2))

    def hess(A, B, f):
        return sp.diff(f, X[A], X[B]) - sum(Gm[C][A][B] * sp.diff(f, X[C]) for C in range(2))

    grad2 = sum(hi[a, b] * sp.diff(r, X[a]) * sp.diff(r, X[b]) for a in range(2) for b in range(2))
    box = sum(hi[a, b] * hess(a, b, r) for a in range(2) for b in range(2))
    # ghost scalar, delta phi = 0: phi = phi0(x), phi0' = P
    dphi = [0, P]
    gradphi2 = hi[1, 1] * P**2
    rhoE = q**2 / (2 * r**(2 * n))
    eqs = {}
    for A in range(2):
        for B in range(A, 2):
            G = -(n / r) * hess(A, B, r) + h[A, B] * ((n / r) * box + n * (n - 1) / (2 * r**2) * (grad2 - 1))
            T = -(dphi[A] * dphi[B] - h[A, B] * gradphi2 / 2) - rhoE * h[A, B]
            eqs[(A, B)] = G - T
    Gang = -Rh / 2 + ((n - 1) / r) * box + (n - 1) * (n - 2) / (2 * r**2) * (grad2 - 1)
    Tang = -(-gradphi2 / 2) + rhoE          # T_ij / (r^2 gamma_ij): ghost gives +(1/2) g_ij (d phi)^2 / ... sign below
    # ghost: T_ij = -(0 - 1/2 r^2 gamma_ij (d phi)^2) = +1/2 r^2 gamma_ij (d phi)^2
    Tang = gradphi2 / 2 + rhoE
    eqs["ang"] = Gang - Tang
    # angular equation in Ricci form (no R[h], so no third derivatives after elimination):
    # R_ij = r^2 gamma_ij [(n-1)(1 - (grad r)^2)/r^2 - box r / r] = kappa (T_ij - g_ij T/(D-2)); the massless scalar
    # gives 0, the electric field 2 rho_E/n r^2 gamma_ij
    eqs["angR"] = (n - 1) * (1 - grad2) - r * box - r**2 * 2 * rhoE / n
    sqrtmh = sp.exp(nu + lam)
    eqs["scalar"] = sp.diff(sqrtmh * r**n * hi[1, 1] * P, x)
    return eqs, E


def first(expr, E):
    return sp.simplify(sp.diff(expr, eps).subs(eps, 0) / E)


def check_warped():
    """compare the warped-product G_AB, G_ij with a direct computation for n = 2, 3 (static metric, generic functions)"""
    bad = 0
    for nn in (2, 3):
        th = sp.symbols(f"th0:{nn}")
        X = [t, x, *th]
        D = nn + 2
        g = sp.zeros(D)
        nuf, lamf, rf = (sp.Function(s)(t, x) for s in ("nuf", "lamf", "rf"))
        g[0, 0] = -sp.exp(2 * nuf)
        g[1, 1] = sp.exp(2 * lamf)
        w = 1
        for i in range(nn):
            g[2 + i, 2 + i] = rf**2 * w
            w = w * sp.sin(th[i]) ** 2
        gi = g.inv()
        Gm = [[[sp.simplify(sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                                for d in range(D)) / 2) for c in range(D)] for b in range(D)] for a in range(D)]

        def ric(b, c):
            return sp.simplify(sum(sp.diff(Gm[a][b][c], X[a]) - sp.diff(Gm[a][b][a], X[c])
                                   + sum(Gm[a][a][d] * Gm[d][b][c] - Gm[a][c][d] * Gm[d][b][a] for d in range(D))
                                   for a in range(D)))
        Ric = sp.Matrix(D, D, ric)
        Rs = sp.simplify(sum(gi[a, b] * Ric[a, b] for a in range(D) for b in range(D)))
        Gd = Ric - Rs * g / 2
        # warped formulas
        h = sp.diag(-sp.exp(2 * nuf), sp.exp(2 * lamf))
        hi = h.inv()
        X2 = [t, x]
        G2 = [[[sum(hi[a, d] * (sp.diff(h[d, b], X2[c]) + sp.diff(h[d, c], X2[b]) - sp.diff(h[b, c], X2[d]))
                    for d in range(2)) / 2 for c in range(2)] for b in range(2)] for a in range(2)]

        def ric2(b, c):
            return sum(sp.diff(G2[a][b][c], X2[a]) - sp.diff(G2[a][b][a], X2[c])
                       + sum(G2[a][a][d] * G2[d][b][c] - G2[a][c][d] * G2[d][b][a] for d in range(2)) for a in range(2))
        Rh = sum(hi[a, b] * ric2(a, b) for a in range(2) for b in range(2))

        def hess(A, B):
            return sp.diff(rf, X2[A], X2[B]) - sum(G2[C][A][B] * sp.diff(rf, X2[C]) for C in range(2))
        grad2 = sum(hi[a, b] * sp.diff(rf, X2[a]) * sp.diff(rf, X2[b]) for a in range(2) for b in range(2))
        box = sum(hi[a, b] * hess(a, b) for a in range(2) for b in range(2))
        for A in range(2):
            for B in range(A, 2):
                Gw = -(nn / rf) * hess(A, B) + h[A, B] * ((nn / rf) * box + nn * (nn - 1) / (2 * rf**2) * (grad2 - 1))
                d = sp.simplify(Gd[A, B] - Gw)
                print(f"n = {nn}: G_{A}{B} direct - warped = {d}")
                bad += d != 0
        Gw = rf**2 * (-Rh / 2 + ((nn - 1) / rf) * box + (nn - 1) * (nn - 2) / (2 * rf**2) * (grad2 - 1))
        d = sp.simplify(Gd[2, 2] - Gw)
        print(f"n = {nn}: G_th0th0 direct - warped = {d}")
        bad += d != 0
    return bad


if __name__ == "__main__":
    if sys.argv[1] == "check":
        sys.exit(1 if check_warped() else 0)
    eqs, E = setup()
    out = {k: first(v, E) for k, v in eqs.items()}
    for k, v in out.items():
        print(k, "=", v)

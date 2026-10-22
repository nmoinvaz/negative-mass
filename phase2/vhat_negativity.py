"""
Where can the vacuum exterior of a negative-mass star host an unstable even-parity mode?

Dotti & Gleiser (arXiv:0809.3615) map the singular Zerilli equation to
    d^2 Psi^/dt^2 - d^2 Psi^/dx^2 + V^ Psi^ = 0,   V^ = V - 2 g',  g = psi0'/psi0  (x = tortoise),
with psi0 = chi(r) / (6M + r(l+2)(l-1)) and chi the polynomial of their eq. (55).

If a star of radius R sits outside the region where V^ < 0, then for any Robin boundary
condition Psi^'(x_R) = kappa Psi^(x_R) with kappa >= 0, the exterior quadratic form
    int (|Psi^'|^2 + V^ |Psi^|^2) dx + kappa |Psi^(x_R)|^2  >= 0,
so the exterior alone cannot support a negative eigenvalue: any instability must be sourced by
the matter (through the interface condition). This script finds r_neg(l) = sup{r : V^(r) < 0}.
Units |M| = 1 (M = -1).
"""
import sympy as sp
import numpy as np

r = sp.Symbol('r', positive=True)
M = sp.Integer(-1)

def chi_poly(l):
    n_terms = []
    for n in range(1, l + 3):
        c = sp.Integer(n - 2) * ((n - 4) * l * (l + 1) + n - 1) * sp.gamma(l + n - 1) \
            / (2**n * sp.gamma(n)**2 * sp.gamma(l - n + 3) * (-M)**n)
        n_terms.append(c * r**n)
    return sp.expand(sum(n_terms))

def zerilli_V(lam):
    f = 1 - 2*M/r
    return 2*f*(lam**2*r**2*((lam+1)*r + 3*M) + 9*M**2*(lam*r + M)) / (r**3*(lam*r + 3*M)**2)

def static_residual(psi, lam):
    f = 1 - 2*M/r
    return sp.simplify(-f*sp.diff(f*sp.diff(psi, r), r) + zerilli_V(lam)*psi)

rows = []
for l in range(2, 9):
    lam = sp.Rational((l - 1)*(l + 2), 2)
    chi = chi_poly(l)
    psi0 = chi / (6*M + r*(l + 2)*(l - 1))
    res = static_residual(psi0, lam) if l <= 5 else None   # symbolic check (slow for big l)
    f = 1 - 2*M/r
    g = f*sp.diff(psi0, r)/psi0                 # d/dx = f d/dr
    Vhat = sp.simplify(zerilli_V(lam) - 2*f*sp.diff(g, r))
    Vf = sp.lambdify(r, Vhat, 'numpy')
    chif = sp.lambdify(r, chi, 'numpy')
    rr = np.geomspace(1e-4, 200, 400001)
    vals = Vf(rr)
    neg = rr[vals < 0]
    rneg = neg.max() if neg.size else 0.0
    rmin_neg = neg.min() if neg.size else 0.0
    vmin = vals.min()
    chi_pos = bool(np.all(chif(rr) > 0))
    k = (l - 1)*l*(l + 1)*(l + 2)/6
    rows.append((l, k, rneg, rmin_neg, vmin, chi_pos, res))
    print(f"l={l}: static residual={res}  chi>0 on r>0: {chi_pos}  "
          f"V^<0 on r in [{rmin_neg:.4f}, {rneg:.4f}]|M|  min V^={vmin:.4g}")

print()
print("Compactness bound: exterior cannot host an instability (kappa>=0) when R > r_neg(l), i.e. |M|/R < 1/r_neg(l)")
print(f"{'l':>2} {'k':>5} {'r_neg/|M|':>10} {'|M|/R bound':>12} {'2sqrt(2/k)':>11}")
for l, k, rneg, _, _, _, _ in rows:
    print(f"{l:>2} {k:>5.0f} {rneg:>10.4f} {1/rneg:>12.4f} {2*np.sqrt(2/k):>11.4f}")
print(f"\nAll l in 2..8: exterior positive whenever |M|/R < {1/max(x[2] for x in rows):.4f}")

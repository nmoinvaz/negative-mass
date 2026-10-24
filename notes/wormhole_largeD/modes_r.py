"""The Ellis wormhole's growing mode in D = n + 2 dimensions, with the areal radius r as variable (throat radius 1).

With F(r) = 1 - r^{-2(n-1)} (= r0'^2) the master equation of master.py becomes a first-order system on one side of the
throat,
  dY/dr = L,      F r dL/dr = W^2 Y r + L (r F'/2 + n F) - n Y F / r - Y F'/2,
which has a regular singular point at the throat r = 1 (F = 0). The even mode (the one that moves the throat) is the
solution analytic in rho = r - 1, fixed by Y(1) = 1; it starts from a Frobenius series. Beyond r_m = 1 + X/(2(n-1))
(F = 1 - e^{-X}) the equation is Y'' - (n/r) Y' + (n/r^2) Y = W^2 Y up to e^{-X}, whose decaying solution is
Y = r^{(n+1)/2} K_{(n-1)/2}(W r). The eigenvalue condition matches Y'/Y there.

    python3 modes_r.py D [D ...]          W^2 at each D
    python3 modes_r.py --conv D           the same with changed numerical parameters
"""
import sys

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import kve

# Frobenius series about rho = 0, to order 3, from sympy
_rho, _n, _W2 = sp.symbols("rho n W2", positive=True)
_K = 4


def _series():
    r = 1 + _rho
    F = sp.series(1 - r**(-2 * (_n - 1)), _rho, 0, _K + 2).removeO()
    Fp = sp.diff(1 - r**(-2 * (_n - 1)), _rho)
    Fp = sp.series(Fp, _rho, 0, _K + 2).removeO()
    y = sp.symbols(f"y0:{_K + 1}")
    lc = sp.symbols(f"l0:{_K + 1}")
    Y = sum(y[k] * _rho**k for k in range(_K + 1))
    L = sum(lc[k] * _rho**k for k in range(_K + 1))
    eq1 = sp.expand(sp.diff(Y, _rho) - L)
    eq2 = sp.expand(F * r * sp.diff(L, _rho) - (_W2 * Y * r + L * (r * Fp / 2 + _n * F)
                                                 - _n * Y * F * sp.series(1 / r, _rho, 0, _K + 2).removeO() - Y * Fp / 2))
    sol = {y[0]: 1}
    eqs = []
    for k in range(_K):
        eqs.append(eq1.coeff(_rho, k))
        eqs.append(eq2.coeff(_rho, k))
    unknowns = [v for v in (*y[1:], *lc) if v not in sol]
    s = sp.solve([e.subs(sol) for e in eqs], unknowns[:2 * _K], dict=True)[0]
    sol.update(s)
    Ys = sp.simplify(Y.subs(sol).subs({v: 0 for v in (*y, *lc)}))
    Ls = sp.simplify(L.subs(sol).subs({v: 0 for v in (*y, *lc)}))
    return sp.lambdify((_rho, _n, _W2), Ys, "numpy"), sp.lambdify((_rho, _n, _W2), Ls, "numpy")


_YS, _LS = _series()


def F_of(rho, n):
    return -np.expm1(-2 * (n - 1) * np.log1p(rho))


def Fp_of(rho, n):
    return 2 * (n - 1) * np.exp(-(2 * n - 1) * np.log1p(rho))


def mismatch(n, W2, X=40.0, rho0_fac=1e-4, rtol=1e-12):
    """L - q_far Y at r_m for the regular solution with Y(1) = 1 (zero at the eigenvalue). Y can have a node, so the
    linear variables are integrated, not log Y."""
    s = 2 * (n - 1)
    rho0 = rho0_fac / s
    rho_m = np.expm1(X / s)          # r_m^{2(n-1)} = e^X

    def f(rho, z):
        Y, L = z
        r = 1 + rho
        F, Fp = F_of(rho, n), Fp_of(rho, n)
        rhs = W2 * Y * r + L * (r * Fp / 2 + n * F) - n * Y * F / r - Y * Fp / 2
        return [L, rhs / (F * r)]

    sol = solve_ivp(f, (rho0, rho_m), [_YS(rho0, n, W2), _LS(rho0, n, W2)], method="DOP853", rtol=rtol, atol=1e-30)
    Y_m, L_m = sol.y[0, -1], sol.y[1, -1]
    r_m = 1 + rho_m
    W = np.sqrt(W2)
    nu = (n - 1) / 2
    z = W * r_m
    # Y = r^{(n+1)/2} K_nu(W r):  Y'/Y = (n+1)/(2r) + W K_nu'/K_nu,  K_nu' = -K_{nu-1} - (nu/z) K_nu
    ratio = kve(nu - 1, z) / kve(nu, z)
    if not np.isfinite(ratio) or ratio == 0:
        import mpmath as mp
        mp.mp.dps = 30
        ratio = float(mp.besselk(nu - 1, z) / mp.besselk(nu, z))
    q_far = (n + 1) / (2 * r_m) + W * (-ratio - nu / z)
    return (L_m - q_far * Y_m) * np.exp(-X / 2)


def W2_of(D, **kw):
    n = D - 2
    g = 2 * n - 4.7 if D > 6 else (1.4 if D == 4 else 2 * n - 1)
    lo, hi = max(0.05, g - 0.4 * max(g, 1)), g + 0.4 * max(g, 1) + 1
    grid = np.linspace(lo, hi, 41)
    vals = [mismatch(n, w, **kw) for w in grid]
    for (a, fa), (b, fb) in zip(zip(grid, vals), zip(grid[1:], vals[1:])):
        if np.isfinite(fa) and np.isfinite(fb) and np.sign(fa) != np.sign(fb) and abs(fa - fb) < 1e3:
            return brentq(lambda w: mismatch(n, w, **kw), a, b, xtol=1e-13, rtol=1e-15)
    return np.nan


if __name__ == "__main__":
    args = sys.argv[1:]
    conv = "--conv" in args
    Ds = [int(a) for a in args if not a.startswith("--")]
    for D in Ds:
        w = W2_of(D)
        line = f"D = {D:6d}: W^2 a^2 = {w:.10f}   2(D-2) - W^2 = {2 * (D - 2) - w:.10f}"
        if conv:
            w2 = W2_of(D, X=50.0, rho0_fac=1e-5, rtol=1e-13)
            line += f"   (X 50, rho0/10, rtol 1e-13: {w2:.10f}, diff {w2 - w:.1e})"
        print(line, flush=True)

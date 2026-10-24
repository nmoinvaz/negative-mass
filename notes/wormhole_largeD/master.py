"""The master equation for radial perturbations of the D-dimensional Ellis wormhole, and its unstable mode.

From derive.py (gauge delta phi = 0, perturbations ~ e^{W t}, throat radius 1, n = D - 2):
  (0,1):  L = Y' / r0'
  (S):    N = L - n Y / r0            (the residual time reparametrization sets the integration function to 0)
  (E2):   Y'' - (2 r0''/r0' + n r0'/r0) Y' + (n r0'^2/r0^2 + r0''/r0) Y = W^2 Y
with the background r0'^2 = 1 - r0^{-2(n-1)}, r0'' = (n - 1) r0^{1-2n}. At the throat (l = 0) the equation has a
regular singular point with exponents 0 and 3; the radial (throat-moving) mode is the even one, Y(0) = 1.

    python3 master.py check          symbolic check that the remaining Einstein equations (tt, ll) follow
    python3 master.py modes D...     W^2 of the growing mode at each D (compare Torii & Shinkai 2013)
"""
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def check(coef_shift=0):
    """coef_shift != 0 is the control: a wrong master equation, which must leave nonzero residuals"""
    import sympy as sp
    import derive as d
    l, n = d.l, d.n
    r0, N, L, Y, W = d.r0, d.N, d.L, d.Y, d.W
    first = {k: d.first(v) for k, v in d.E1.items()}
    first["E2"] = d.first(d.E2)
    # substitute L, N in terms of Y, then Y'' from the master equation, then the background relations
    Ls = sp.diff(Y, l) / sp.diff(r0, l)
    Ns = Ls - n * Y / r0
    rp, rpp = sp.diff(r0, l), sp.diff(r0, l, 2)
    Ypp = (2 * rpp / rp + (n + coef_shift) * rp / r0) * sp.diff(Y, l) - (n * rp**2 / r0**2 + rpp / r0) * Y + W**2 * Y
    rr = sp.Symbol("rr", positive=True)
    q = sp.Symbol("q", positive=True)   # r0'
    bg_rpp = (n - 1) * rr**(1 - 2 * n)
    bg_rppp = (n - 1) * (1 - 2 * n) * rr**(-2 * n) * q
    bad = 0
    for k, e in first.items():
        e = e.subs({L: Ls, N: Ns}).doit()
        # eliminate Y''' and Y'' using the master equation and its derivative
        Yppp = sp.diff(Ypp, l)
        e = e.subs(sp.Derivative(Y, (l, 3)), Yppp).subs(sp.Derivative(Y, (l, 2)), Ypp)
        e = e.subs(sp.Derivative(Y, (l, 2)), Ypp)
        e = e.subs(sp.Derivative(r0, (l, 3)), bg_rppp).subs(sp.Derivative(r0, (l, 2)), bg_rpp)
        e = e.subs(sp.Derivative(r0, l), q).subs(r0, rr)
        e = e.subs(q, sp.sqrt(1 - rr**(-2 * (n - 1))))
        if coef_shift:
            # the control only needs to show a nonzero residual: evaluate at a generic point
            v = e.subs({sp.Derivative(Y, l): sp.Rational(37, 100)}).subs({Y: sp.Rational(11, 10)})
            v = sp.N(v.subs({n: 3, rr: sp.Rational(13, 10), W: sp.Rational(7, 10)}))
            ok = abs(v) < 1e-12
            print(k, "-> residual at a generic point", v)
        else:
            e = sp.simplify(e)
            ok = e == 0
            print(k, "->", e, "(ok)" if ok else "(NOT ZERO)")
        bad += not ok
    # control: a wrong master equation (n -> n + 1 in the Y' coefficient) must leave a nonzero residual
    return bad


def background_and_mode(n, W2, lmax, l0=1e-4, rtol=1e-11):
    """integrate background (r, r') and Y, Y' from l0 to lmax for growth rate^2 W2; returns Y(lmax) e^{-W lmax}"""
    W = np.sqrt(W2)
    c = (n - 1 - W2) / 2                       # Y = 1 + c l^2 + ... (with r0 = 1 + (n-1) l^2/2 + ...)
    r_init = 1 + (n - 1) * l0**2 / 2
    rp_init = (n - 1) * l0

    def f(l, s):
        r, rp, y, yp = s
        rpp = (n - 1) * r**(1 - 2 * n)
        ypp = (2 * rpp / rp + n * rp / r) * yp - (n * rp**2 / r**2 + rpp / r) * y + W2 * y
        return [rp, rpp, yp, ypp]

    sol = solve_ivp(f, (l0, lmax), [r_init, rp_init, 1 + c * l0**2, 2 * c * l0], method="DOP853", rtol=rtol,
                    atol=1e-14, dense_output=False)
    y = sol.y[2, -1]
    return y * np.exp(-W * lmax)


def growing_mode(D, lo=None, hi=None, lmax=None):
    n = D - 2
    lmax = lmax or 12.0
    lo = lo if lo is not None else 0.2
    hi = hi if hi is not None else 4.0 * D
    # scan for a sign change
    grid = np.linspace(lo, hi, 80)
    vals = [background_and_mode(n, w, lmax) for w in grid]
    roots = []
    for (w1, f1), (w2, f2) in zip(zip(grid, vals), zip(grid[1:], vals[1:])):
        if np.sign(f1) != np.sign(f2):
            roots.append(brentq(lambda w: background_and_mode(n, w, lmax), w1, w2, xtol=1e-12))
    return roots


if __name__ == "__main__":
    if sys.argv[1] == "check":
        bad = check()
        print("control (wrong Y' coefficient):")
        fired = check(coef_shift=1) > 0
        print("control", "fires" if fired else "DOES NOT FIRE")
        sys.exit(1 if bad or not fired else 0)
    for D in [int(a) for a in sys.argv[2:]]:
        r = growing_mode(D)
        print(f"D = {D:3d}: W^2 a^2 = {', '.join(f'{x:.6f}' for x in r)}", flush=True)

"""The master equation for radial perturbations of the general static ghost wormhole (asymmetric, charged), in
D = n + 2 dimensions, in a harmonic radial coordinate y in (-pi/2, pi/2).

Background (from derive_general.py's equations; bg_derive.py prints them): lam = nu + n beta + c, r = e^beta,
  eta := nu + (n-1) beta = -ln((n-1) e^c cos y),
  nu'' = eps^2 e^{2 nu},  eps^2 = ((n-1)/n) q^2 e^{2c}:   nu = -ln[cosh(Lam y) - gamma1 sinh(Lam y)/Lam],
  Lam^2 = gamma1^2 - eps^2 (nu(0) = 0, nu'(0) = gamma1; Lam imaginary or zero handled by analytic continuation),
  P^2 = phi0'^2 = n (1 + Lam^2) / (n - 1)   (the xx constraint).
q = 0 gives the uncharged asymmetric family (nu = gamma1 y); gamma1 = q = 0 the symmetric one.

Perturbation (gauge delta phi = 0): L = (Y' - nu' Y)/r', N = L - n Y/r, and the angular equation (Ricci form)
  A2 Y'' + A1 Y' + A0 Y = 0  (coefficients printed by `python3 master_general.py coeffs`).
`python3 master_general.py verify` checks the background equations and that the tt, xx equations hold on solutions of
the master equation (symbolic in n, with sample numeric parameters).
"""
import sys

import sympy as sp

import derive_general as d

y = d.x
n, q, W = d.n, d.q, d.W
c, g1, e2 = sp.symbols("c gamma1 eps2", real=True)
Lam = sp.sqrt(g1**2 - e2)


def background_exprs():
    if True:
        nu = -sp.log(sp.cosh(Lam * y) - g1 * sp.sinh(Lam * y) / Lam)
    eta = -sp.log((n - 1) * sp.exp(c) * sp.cos(y))
    beta = (eta - nu) / (n - 1)
    lam = nu + n * beta + c
    P = sp.sqrt(n * (1 + Lam**2) / (n - 1))
    qq = sp.sqrt(n * e2 * sp.exp(-2 * c) / (n - 1))
    return nu, lam, sp.exp(beta), P, qq


def first_order():
    eqs, E = d.setup()
    return {k: d.first(v, E) for k, v in eqs.items()}


def master_coeffs(first=None):
    first = first or first_order()
    nu0, lam0, r0, Y = d.nu0, d.lam0, d.r0, d.Y
    Ls = (sp.diff(Y, y) - sp.diff(nu0, y) * Y) / sp.diff(r0, y)
    Ns = Ls - n * Y / r0
    ang = sp.expand(first["angR"].subs({d.L: Ls, d.N: Ns}).doit())
    Ypp, Yp = sp.Derivative(Y, (y, 2)), sp.Derivative(Y, y)
    a2 = ang.coeff(Ypp)
    rest = sp.expand(ang - a2 * Ypp)
    a1 = rest.coeff(Yp)
    a0 = sp.expand(rest - a1 * Yp).coeff(Y)
    assert sp.simplify(ang - (a2 * Ypp + a1 * Yp + a0 * Y)) == 0
    return a2, a1, a0, first, Ls, Ns


if __name__ == "__main__":
    a2, a1, a0, first, Ls, Ns = master_coeffs()
    if sys.argv[1] == "coeffs":
        for name, v in (("A2", a2), ("A1", a1), ("A0", a0)):
            print(name, "=", sp.simplify(v))
    elif sys.argv[1] == "verify":
        nu, lam, r, P, qq = background_exprs()
        vals = {n: 5, c: sp.Rational(1, 3), g1: sp.Rational(7, 10), e2: sp.Rational(3, 10)}
        bg = {d.nu0: nu, d.lam0: lam, d.r0: r, d.P: P, d.q: qq}
        eqs0, _ = d.setup(perturbed=False)
        bad = 0
        for k, v in eqs0.items():
            val = v.subs(bg).doit().subs(vals)
            num = [sp.N(val.subs(y, yy)) for yy in (sp.Rational(-1, 2), sp.Rational(1, 5), sp.Rational(9, 10))]
            print("background", k, [f"{complex(z):.2e}" for z in num])
            bad += any(abs(complex(z)) > 1e-12 for z in num)
        # the tt and xx equations on a solution of the master equation: eliminate Y'' (and Y''') with it
        Y = d.Y
        Ypp_sol = -(a1 * sp.Derivative(Y, y) + a0 * Y) / a2
        for k in ((0, 0), (1, 1)):
            e = first[k].subs({d.L: Ls, d.N: Ns}).doit()
            e = e.subs(sp.Derivative(Y, (y, 3)), sp.diff(Ypp_sol, y)).subs(sp.Derivative(Y, (y, 2)), Ypp_sol)
            e = e.subs(sp.Derivative(Y, (y, 2)), Ypp_sol)
            e = e.subs(bg).doit().subs(vals)
            Yd, Yv = sp.Symbol("Yd"), sp.Symbol("Yv")
            e = e.subs(sp.Derivative(Y, y), Yd).subs(Y, Yv)
            e = e.subs({Yd: sp.Rational(-3, 7), Yv: sp.Rational(11, 10), W: sp.Rational(13, 10)})
            if e.free_symbols - {y}:
                print("leftover symbols:", e.free_symbols)
            num = [sp.N(e.subs(y, yy)) for yy in (sp.Rational(-1, 2), sp.Rational(1, 5), sp.Rational(9, 10))]
            print("on-shell", k, [f"{complex(z):.2e}" for z in num])
            bad += any(abs(complex(z)) > 1e-10 for z in num)
        sys.exit(1 if bad else 0)

#!/usr/bin/env python3
"""The ghost (test) field in choptuik-large-d's D = infinity inner problem (its README, Proposition 5), by sympy.

The inner equations are Y_t = c S, S_t = c (Y_s + Y^3/4), c_s = c Y^2/4 (s = r^2/(2n), t the time). They come from
Rozali and Way's leading-order system phi_s = c^-1 (c^-1 phi_t)_t, c_s = phi_t^2/c with Y = 2 phi_t/c, S = Y_t/c.
A second massless field psi = 0 in the background obeys, at linear order, the same wave equation with c frozen:
psi_s = c^-1 (c^-1 psi_t)_t. Checked here: with Psi = 2 psi_t/c and Sigma = Psi_t/c this is

    Psi_t = c Sigma,   Sigma_t = c (Psi_s + Y^2 Psi/4) = (c Psi)_s ,

the inner system linearized with dc = 0 and the cubic's coefficient 3/4 replaced by 1/4; psi = phi (the boost) solves
it. Control: the coefficient 3/4 (the full linearization's, with dc dropped) must leave a nonzero residual.
Exit 1 if the check fails or the control passes.
"""
import sys

import sympy as sp

s, t = sp.symbols("s t")
c = sp.Function("c")(s, t)
psi = sp.Function("psi")(s, t)
phi = sp.Function("phi")(s, t)
Y = 2 * sp.diff(phi, t) / c
Psi = 2 * sp.diff(psi, t) / c
Sig = sp.diff(Psi, t) / c
wave = lambda f: sp.diff(f, s) - sp.diff(sp.diff(f, t) / c, t) / c          # = 0 for a field on the background
rules = {sp.Derivative(c, s): c * Y ** 2 / 4}                               # c_s = phi_t^2 / c
bad = 0
for coef, want_zero in ((sp.Rational(1, 4), True), (sp.Rational(3, 4), False)):
    res = sp.diff(Sig, t) - c * (sp.diff(Psi, s) + coef * Y ** 2 * Psi)
    # eliminate psi_s through the wave equation: psi_s = c^-1 (c^-1 psi_t)_t, differentiated in t
    psi_s = sp.diff(sp.diff(psi, t) / c, t) / c
    res = res.subs(sp.Derivative(psi, s, t), sp.diff(psi_s, t)).subs(sp.Derivative(psi, t, s), sp.diff(psi_s, t))
    res = sp.simplify(res.subs(rules).doit().subs(rules))
    ok = (res == 0) == want_zero
    print(f"coefficient {coef}: residual {'0' if res == 0 else 'nonzero'} -> {'PASS' if ok else 'FAIL'}"
          f"{'' if want_zero else ' (control)'}")
    bad += not ok
# conservation form: c (Psi_s + Y^2 Psi / 4) = (c Psi)_s given c_s = c Y^2/4
cons = sp.simplify((sp.diff(c * Psi, s) - c * (sp.diff(Psi, s) + Y ** 2 * Psi / 4)).subs(rules))
print(f"(c Psi)_s - c (Psi_s + Y^2 Psi/4) = {cons}")
bad += cons != 0
# the boost: psi = phi solves the test equation exactly
print("psi = phi solves psi_s = c^-1 (c^-1 psi_t)_t: identical to the background's own equation (by construction)")
sys.exit(1 if bad else 0)

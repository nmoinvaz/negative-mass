#!/usr/bin/env python3
"""psi_0,xx = V psi_0 exactly, for l = 2..40 (RESULTS section 1 Lemma A(ii) and section 4.2 (Z)).

Direct from the definitions, with no P's and no recurrence: psi_0 = chi/(cr-6) with chi from the closed form of Lemma A,
d/dx = f d/dr with f = 1 + 2/r (M = -1), V = 2f[lambda^2 r^2((lambda+1) r - 3) + 9(lambda r - 1)]/(r^3 (lambda r - 3)^2), lambda = c/2.
The residual f (f psi')' - V psi is a rational function of r over Q, reduced exactly (sympy cancel); it must be 0.
"""
import time

import sympy as sp

import review_lemma_a_common as C

print(f"review_lemma_a_zerilli (mutation: {C.MUT})")
C.run_planted({"PL7"})
C.gate_planted()

r = C.R


def z1(l):
    cl = C.c_of_l(l)
    chi = C.chi_poly(l).as_expr()
    psi = chi / (cl * r - 6)
    V = C.V_zerilli(r, sp.Integer(cl))
    t0 = time.time()
    res = C.zerilli_residual(psi, V, r)
    dt = time.time() - t0
    print(f"  l={l}: degree of chi {sp.degree(chi, r)}, residual {'0' if res == 0 else 'NONZERO'} ({dt:.1f} s)")
    if res == 0:
        return True, ""
    num, den = sp.fraction(sp.together(res))
    return False, f"residual numerator of degree {sp.degree(num, r)}"


C.check_range("Z1 psi_0,xx = V psi_0 exactly (rational functions over Q)", list(range(2, 41)), z1)
C.finish()

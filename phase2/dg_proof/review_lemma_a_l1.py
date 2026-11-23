#!/usr/bin/env python3
"""L1 of RESULTS.md section 1: the chi-ODE at r_s = 6/c has exponents 0 and 3, hence chi'/chi(r_s) = -1/(r_s(r_s+2)) and chi''(r_s) = 0.

Hand derivation (stated here as lemmas for Lean, #135). ODE: P2 chi'' + P1 chi' + P0 chi = 0 with P2 = -r^2(r+2)(cr-6),
P1 = 2cr^3 + 6cr^2 - 12r, P0 = c^2r^2 - 6cr + 12 (derived from V in review_lemma_a_recur.py, A2.1). At r_s (c r_s = 6):
  L1-L1  P2(r_s) = 0, p := P2'(r_s) = -c r_s^2 (r_s+2) = -6 r_s (r_s+2) != 0 (simple zero), q := P1(r_s) = 12 r_s (r_s+2), P0(r_s) = 12.
  L1-L2  Regular singular point (P2 has a simple zero, so P1/P2 has a simple pole and P0/P2 a simple pole). Indicial equation
         rho(rho-1) + (q/p) rho = 0 (the P0/P2 pole is of order 1, not 2), roots 0 and 1 - q/p = 3.
  L1-L3  The ODE at r_s: q chi' + 12 chi = 0, so chi'(r_s) = -chi(r_s)/(r_s(r_s+2)) (no assumption on chi(r_s)).
  L1-L4  Differentiate once, evaluate at r_s (P2 = 0): (p + q) chi'' + (P1' + P0) chi' + P0' chi = 0, with p + q = 6 r_s (r_s+2) != 0,
         P1'(r_s) + P0(r_s) = 36 (r_s+2), P0'(r_s) = 6c = 36/r_s. With L1-L3: chi'' = [36 (r_s+2) / (r_s (r_s+2)) - 36/r_s] chi / (6 r_s (r_s+2)) = 0.
  L1-L5  Differentiate twice: the coefficient of chi''' is 2p + q = 0 (the resonance of exponents 0 and 3), and the rest vanishes with L1-L3, L1-L4;
         so the exponent-0 solution has no logarithm and chi analytic at r_s is consistent.
  L1-L6  Residue of psi_0^2 dx/dr at r_s: psi_0^2 dx/dr = chi^2 r / (c^2 (r-r_s)^2 (r+2)), residue g'(r_s) with g = chi^2 r/(c^2 (r+2)): 0 by L1-L3.
         phi_0 = (cr-6)/chi: phi_0,xx(r_s) = f^2 phi'' + f f' phi' = 0 by L1-L3 (f = (r+2)/r).
"""
import sympy as sp

import review_lemma_a_common as C

print(f"review_lemma_a_l1 (mutation: {C.MUT})")
C.run_planted({"PL3", "PL4"})
C.gate_planted()

r = C.R
c = sp.Symbol("c")
rs = 6 / c
P2, P1, P0 = C.P_hand(r, c)

# ---- L1.1a exponents, symbolic in c
exps = C.exponents(P2, P1, P0, r, rs)
C.check("L1.1a exponents of the chi-ODE at r_s are {0, 3}", dict(exps) == {0: 1, 3: 1}, f"got {dict(exps)}")

# ---- L1.1b local values
vals = {
    "P2(r_s) = 0": (P2.subs(r, rs), 0),
    "P2'(r_s) = -6 r_s (r_s+2)": (sp.diff(P2, r).subs(r, rs), -6 * rs * (rs + 2)),
    "P1(r_s) = 12 r_s (r_s+2)": (P1.subs(r, rs), 12 * rs * (rs + 2)),
    "P0(r_s) = 12": (P0.subs(r, rs), 12),
}
okb = True
msgs = []
for name, (got, want) in vals.items():
    d = sp.simplify(got - want)
    if d != 0:
        okb = False
        msgs.append(f"{name}: difference {sp.factor(d)}")
C.check("L1.1b local values p, q, P0(r_s), P2(r_s)", okb, "; ".join(msgs))

# ---- L1.2, L1.3, L1.4 from the ODE and its derivatives (generic y, built by the routine of PL4)
rel = C.local_relations(P2, P1, P0, r, rs)
claim = -1 / (rs * (rs + (1 if C.MUT == "ratio_claim" else 2)))
C.check("L1.2 chi'(r_s)/chi(r_s) = -1/(r_s (r_s+2))", sp.simplify(rel["yp"] - claim) == 0, f"ODE gives {sp.simplify(rel['yp'])}")
C.check("L1.3 chi''(r_s) = 0", sp.simplify(rel["ypp"]) == 0, f"ODE gives chi''/chi = {sp.simplify(rel['ypp'])}")
C.check("L1.4 resonance: coefficient of chi''' is 2p+q = 0 and the remainder vanishes",
        sp.simplify(rel["res_coef"]) == 0 and sp.simplify(rel["res_rest"]) == 0,
        f"coefficient {sp.simplify(rel['res_coef'])}, remainder {sp.simplify(rel['res_rest'])}")

# ---- L1.5 consequences (symbolic in c), using the relations of the ODE
chi = sp.Function("chi")
Y0, Y1, Y2 = sp.symbols("Y0 Y1 Y2")


def at_rs(expr):
    expr = expr.doit()
    expr = expr.subs(sp.Derivative(chi(r), (r, 2)), rel["ypp"] * Y0).subs(sp.Derivative(chi(r), r), rel["yp"] * Y0)
    return sp.simplify(expr.subs(chi(r), Y0).subs(r, rs))


# residue: psi_0^2 dx/dr = g(r)/(r-r_s)^2 with g = chi^2 r/(c^2 (r+2)); residue = g'(r_s)
g = chi(r) ** 2 * r / (c**2 * (r + 2))
res = at_rs(sp.diff(g, r))
C.check("L1.5a residue of psi_0^2 dx/dr at r_s is 0", res == 0, f"residue {res}")
# phi_0,xx at r_s
phi = (c * r - 6) / chi(r)
f = (r + 2) / r
phixx = f * sp.diff(f * sp.diff(phi, r), r)
val = at_rs(phixx)
C.check("L1.5b phi_0,xx(r_s) = 0", val == 0, f"value {val}")


# ---- L1.6 the actual chi, exact
def data(l):
    cl = C.c_of_l(l)
    chip = C.chi_poly(l)
    rsq = sp.Rational(6, cl)
    v0, v1, v2 = chip.eval(rsq), chip.diff(r).eval(rsq), chip.diff(r).diff(r).eval(rsq)
    if not v0 > 0:
        return False, f"chi(r_s) = {v0}"
    if v1 / v0 != -1 / (rsq * (rsq + 2)):
        return False, f"chi'/chi = {v1 / v0}, want {-1 / (rsq * (rsq + 2))}"
    if v2 != 0:
        return False, f"chi''(r_s) = {v2}"
    return True, ""


C.check_range("L1.6 actual chi: chi(r_s) > 0, chi'/chi = -1/(r_s(r_s+2)), chi''(r_s) = 0 (exact)", list(range(2, 41)), data)
C.finish()

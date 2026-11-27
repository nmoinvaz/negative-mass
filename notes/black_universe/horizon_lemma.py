#!/usr/bin/env python3
"""Symbolic checks of the horizon lemma (#261, horizon_lemma.md): a non-extremal Killing horizon at an extremum of the areal
radius carries a gauge-invariant coefficient c2 = d_V dr, constant along the horizon.

    python3 horizon_lemma.py [--mutate=nonstatic|rp-nonzero|a-not-fr|wrong-ricci|l2-short]

L0  R_VV = -(n/r)(r_VV + F_V r_V) for -2 e^{-F} dV dU + r^2 dOmega_n^2, by a brute-force Ricci tensor at n = 2 and 3.
L1  For a background that depends on (V, U) through Pi = -V U only, the linearized C_VV on U = 0 is dr_VV (any n, eps, potential).
L2  d_V (L_xi r) on U = 0 equals -R'(rho_h) rho'(0) d_V(V xi^U)(V, 0) for every vector field xi: zero for all xi iff R'(rho_h) = 0.
L3  A(rho) = f(R(rho)) with R'(rho_h) = 0 and f(R(rho_h)) = 0 gives A'(rho_h) = 0 (kappa = 0, extremal).
L4  The members: the c = -1 black universe (R'(0) = 0, kappa = 2/pi), c = -1.5 (R'(rho_h) != 0), Schwarzschild (R' = 1),
    Simpson-Visser at a = 2m (kappa = 0).
Exit 1 when a check fails.
"""
import sys

import mpmath as mp
import sympy as sp

MUT = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mutate=")), None)
RESULTS = []


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


V, U = sp.symbols("V U", real=True)
eps, n = sp.symbols("epsilon n")

# L0: R_VV of the warped double-null metric, brute force at n = 2, 3
def ricci_VV(nn):
    th = sp.symbols(f"t0:{nn}")
    X = [V, U, *th]
    r = sp.Function("r")(V, U)
    F = sp.Function("F")(V, U)
    dim = len(X)
    g = sp.zeros(dim, dim)
    g[0, 1] = g[1, 0] = -sp.exp(-F)
    w = r ** 2
    for k in range(nn):
        g[2 + k, 2 + k] = w
        w = w * sp.sin(th[k]) ** 2
    gi = g.inv()
    Gam = [[[sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(dim)) / 2
             for c in range(dim)] for b in range(dim)] for a in range(dim)]
    b = d = 0
    e = sum(sp.diff(Gam[a][b][d], X[a]) - sp.diff(Gam[a][b][a], X[d]) for a in range(dim))
    e += sum(Gam[a][a][k] * Gam[k][b][d] - Gam[a][d][k] * Gam[k][b][a] for a in range(dim) for k in range(dim))
    want = -(nn / r) * (sp.diff(r, V, 2) + sp.diff(F, V) * sp.diff(r, V))
    if MUT == "wrong-ricci":
        want = -((nn + 1) / r) * (sp.diff(r, V, 2) + sp.diff(F, V) * sp.diff(r, V))
    return sp.simplify(e - want)


res0 = [ricci_VV(k) for k in (2, 3)]
check("L0", all(x == 0 for x in res0), f"R_VV + (n/r)(r_VV + F_V r_V) at n = 2, 3: {res0}")

# L1: the linearized constraint on U = 0
Pi = -V * U
rb, Fb, pb = sp.Function("rb"), sp.Function("Fb"), sp.Function("pb")
r0, F0, p0 = rb(Pi), Fb(Pi), pb(Pi)
if MUT == "nonstatic":
    r0 = r0 + V ** 2 * sp.Function("w")(U)        # a background that depends on V other than through Pi
dr, dF, dp = (sp.Function(s)(V, U) for s in ("dr", "dF", "dp"))
a = sp.Symbol("a")
def C_VV(r, F, p):  # noqa: E302
    return sp.diff(r, V, 2) + sp.diff(F, V) * sp.diff(r, V) + (r / n) * eps * sp.diff(p, V) ** 2
lin = sp.diff(C_VV(r0 + a * dr, F0 + a * dF, p0 + a * dp), a).subs(a, 0)
lin_h = sp.simplify(lin.doit().subs(U, 0))
target = sp.diff(dr, V, 2).subs(U, 0)
check("L1", sp.simplify(lin_h - target) == 0, f"linearized C_VV at U = 0 minus dr_VV: {sp.simplify(lin_h - target)}")

# L2: gauge, d_V(L_xi r) on U = 0 for an arbitrary xi = xi^V d_V + xi^U d_U
Rf, rho = sp.Function("R"), sp.Function("rho")
rr = Rf(rho(Pi))
xV, xU = sp.Function("xiV")(V, U), sp.Function("xiU")(V, U)
Lr = xV * sp.diff(rr, V) + xU * sp.diff(rr, U)
dLr = sp.simplify(sp.diff(Lr, V).subs(U, 0).doit())
s0 = sp.Symbol("s0")
expect = -sp.diff(V * xU, V).subs(U, 0) * sp.Subs(sp.Derivative(Rf(s0), s0), s0, rho(0)).doit() * sp.Subs(sp.Derivative(rho(s0), s0), s0, 0).doit()
if MUT == "l2-short":                                # the form without the r_U d_V xi^U term (the first draft of this check)
    expect = -xU.subs(U, 0) * sp.Subs(sp.Derivative(Rf(s0), s0), s0, rho(0)).doit() * sp.Subs(sp.Derivative(rho(s0), s0), s0, 0).doit()
diff2 = sp.simplify(dLr - expect)
check("L2", diff2 == 0, f"d_V(L_xi r)|_(U=0) + R'(rho(0)) rho'(0) d_V(V xi^U)|_(U=0) = {diff2}; so c2 is invariant for every xi iff R'(rho_h) = 0")

# L3: A = f(R(rho)) at an extremum of R is extremal
x, xh = sp.symbols("x x_h")
f = sp.Function("f")
Rg = sp.Function("Rg")
A = f(Rg(x))
if MUT == "a-not-fr":
    A = A + (x - xh)
Ap = sp.diff(A, x).subs(x, xh)
Ap = Ap.subs(sp.Subs(sp.Derivative(Rg(x), x), x, xh).doit(), 0)
Ap = sp.simplify(Ap.doit().subs(sp.Derivative(Rg(xh), xh), 0))
check("L3", sp.simplify(Ap) == 0, f"A'(rho_h) for A = f(R(rho)) with R'(rho_h) = 0: {Ap}")

# L4: the members
mp.mp.dps = 30
def bu(c):  # noqa: E302
    c = mp.mpf(c)
    K = -2 * c / mp.pi
    A_ = lambda y: 1 + c * (y * y + 1) + K * (y + (y * y + 1) * mp.atan(y))  # noqa: E731
    rh = mp.findroot(A_, mp.mpf(0) if c == -1 else mp.mpf(0.5))
    return rh, rh / mp.sqrt(rh * rh + 1), mp.diff(A_, rh) / 2
rows = []
rh1, rp1, k1 = bu(-1)
rh15, rp15, k15 = bu("-1.5")
if MUT == "rp-nonzero":
    rp1 = rp1 + mp.mpf("0.01")
rows.append(("c = -1 black universe", rp1, k1, abs(rp1) < mp.mpf("1e-25") and k1 > 0, "c2 physical"))
rows.append(("c = -1.5 black universe", rp15, k15, abs(rp15) > mp.mpf("0.1") and k15 > 0, "c2 gauge"))
rows.append(("Schwarzschild (r = rho, 2M = 1)", mp.mpf(1), mp.mpf("0.5"), True, "c2 gauge"))
a_sv = mp.mpf(1)                                       # Simpson-Visser, a = 2m = 1: A = 1 - 2m/sqrt(rho^2 + a^2)
A_sv = lambda y: 1 - a_sv / mp.sqrt(y * y + a_sv * a_sv)  # noqa: E731
rows.append(("Simpson-Visser a = 2m", mp.mpf(0), mp.diff(A_sv, 0) / 2, abs(mp.diff(A_sv, 0)) < mp.mpf("1e-25"), "extremal, outside the lemma"))
for name, rp, k, ok, verdict in rows:
    print(f"  {name}: R'(rho_h) = {mp.nstr(rp, 6)}, kappa = {mp.nstr(k, 10)}: {verdict}")
check("L4", all(r[3] for r in rows) and abs(k1 - 2 / mp.pi) < mp.mpf("1e-25"),
      f"c = -1: rho_h = {mp.nstr(rh1, 3)}, R'(rho_h) = {mp.nstr(rp1, 3)}, kappa - 2/pi = {mp.nstr(k1 - 2 / mp.pi, 3)}; "
      f"c = -1.5: R'(rho_h) = {mp.nstr(rp15, 6)} (#248: 0.3216)")

bad = [nm for nm, ok in RESULTS if not ok]
print(f"\n{len(RESULTS)} checks, {len(RESULTS) - len(bad)} pass; failed: {bad}; mutation {MUT}", flush=True)
sys.exit(1 if bad else 0)

#!/usr/bin/env python3
"""Multiplicative separable shifts u = X(x) R(rho) + v with det R_3 = 0 (issue #160).

D = det R_3 of the static shift u (det_classify.md Setup) is retyped here:
  D = h (h^2 - P^2 - Q^2 + Y^2 + Z^2) - P (Y^2 - Z^2) + 2 Q Y Z,
  P = -(u_x^2 + u u_xx), Q = -(u_rho_rho + u_rho/rho)/2, Z = u_x_rho/2, Y = u_x u_rho + u u_x_rho/2, h = u_rho^2/2.
For u = X0 R0 + v: E = 8 D / R1^2 is a polynomial in the x-jets (X0, X1, X2) and the rho-jets (R0, R1, S0), S0 = R'' + R'/rho.
E is a bilinear form, E = sum_{i<=r} F_i(x) G_i(rho). If E = 0 on a box then span{F(x)} and span{G(rho)} are orthogonal in R^r,
so d_F + d_G <= r (Lemma L). The script
  C1  builds E from D and compares it with det_classify.E_mult_poly, with the sum-of-squares form and with numerical derivatives;
  C2  v = 0: the six-term separation, rank 6;
  C3  Lemma L on planted random data;
  C4  Lemma X (x-side, d_F >= 4 unless m = (X'/X)^2 is constant), C5 Lemma R (rho-side, d_G >= 3), C6 m constant;
  C7  planted subcases b(rho) and planar; C8-C12 planted additive u = X + R and its perturbations;
  C13 v != 0: ten-term separation, rank 10; C14 the rank bound alone does not decide v != 0.
Imported from det_classify: E_mult_poly, E_add_poly (comparison targets only). Shared with it: the formula for D.
Usage: mult_classify.py [--no-controls] [--mutate=NAME] [--only=c1,c2]. Exit 1 if a check fails or a control does not fire."""
import os
import random
import subprocess
import sys
import time
from fractions import Fraction

import mpmath as mp
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import det_classify as dc  # noqa: E402  (E_mult_poly, E_add_poly only)

# name: (target check, what the mutation does)
MUTATIONS = {
    "m-d-sign": ("e_from_d", "Y := u_x u_rho - u u_x_rho/2 in D"),
    "m-theta": ("decomp_static", "Theta with 19 m^2 in place of 20 m^2 in f_3"),
    "m-rank-off": ("rank_lemma", "the planted pass uses N' = N with one column not orthogonal to col(B)"),
    "m-fdrop": ("lemma_x", "f_3 and f_5 (rho Theta, Psi) set to 0 in the rank test"),
    "m-gdrop": ("lemma_rho", "g_4 := c and g_6 := 1"),
    "m-const-sign": ("const_m", "the exact product test over c = -n/10 in place of c = +n/10"),
    "m-sub-drop": ("mult_subcases", "E_sub built with S0 in place of S0^2"),
    "m-add-drop": ("add_decomp", "u = X0 in the additive jets, i.e. (X0 + R0)^2 replaced by X0^2"),
    "m-stratum": ("add_stratum", "f_4 and f_5 dropped from the elimination"),
    "m-add-freq": ("add_planted", "J0(c rho), Y0(c rho) in place of J0(sqrt2 c rho), Y0(sqrt2 c rho)"),
    "m-perturb-zero": ("add_perturbed", "delta = 0"),
    "m-xperturb-zero": ("add_xperturbed", "X = c x (no x^2/10)"),
    "m-trav-rank": ("trav_decomp", "rank expected 9 at v = 1/3"),
    "m-gap-d": ("trav_gap", "R = rho^2 replaced by R = exp(rho)"),
}

RESULTS = []
OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate")
ONLY = OPTS["only"].split(",") if "only" in OPTS else None
mp.mp.dps = 130
TOL_ZERO = mp.mpf(10) ** -40  # singular values below TOL_ZERO * largest count as zero

X0, X1, X2, R0, R1, S0, v = sp.symbols("X0 X1 X2 R0 R1 S0 v")
p1, be, c1, s2 = sp.symbols("p1 be c1 s2")
XV = (X0, X1, X2)
RV = (R0, R1, S0)


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def wanted(name):
    return ONLY is None or name in ONLY


# ---------------------------------------------------------------------------------------------------------------
# D, E and the bilinear form
# ---------------------------------------------------------------------------------------------------------------
def D_from_jets(u, ux, ur, uxx, uxr, lap, ysign=1):
    P = -(ux**2 + u * uxx)
    Q = -lap / 2
    Z = uxr / 2
    Y = ux * ur + ysign * u * uxr / 2
    h = ur**2 / 2
    return h * (h**2 - P**2 - Q**2 + Y**2 + Z**2) - P * (Y**2 - Z**2) + 2 * Q * Y * Z


def E_mult(vv, ysign=1):
    """E = 8 D / R1^2 for u = X0 R0 + vv, as a polynomial in the jets."""
    u = X0 * R0 + vv
    D = D_from_jets(u, X1 * R0, X0 * R1, X2 * R0, X1 * R1, X0 * S0, ysign)
    return sp.expand(sp.cancel(8 * D / R1**2))


def bilinear_form(E, xvars=XV, rvars=RV):
    """Coefficient matrix of the polynomial E in (x-monomials) x (rho-monomials), and the monomial lists."""
    P = sp.Poly(sp.expand(E), *xvars, *rvars)
    nx = len(xvars)
    rows, cols, ent = {}, {}, {}
    for mon, co in P.terms():
        ri = rows.setdefault(mon[:nx], len(rows))
        ci = cols.setdefault(mon[nx:], len(cols))
        ent[(ri, ci)] = co
    M = sp.zeros(len(rows), len(cols))
    for (i, j), co in ent.items():
        M[i, j] = co
    return M, list(rows), list(cols)


def mat_rank(E):
    return bilinear_form(E)[0].rank()


def frac_rank(rows):
    """Exact rank of a matrix given as a list of rows of Fractions."""
    M = [list(r) for r in rows]
    rk, ncol = 0, len(M[0]) if M else 0
    for col in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][col] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][col] != 0:
                f = M[i][col] / M[rk][col]
                M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


def num_rank(rows):
    """Numerical rank (mpmath, 130 digits): singular values <= 1e-40 * largest count as zero."""
    s = mp.svd_r(mp.matrix(rows), compute_uv=False)
    smax = max(s)
    return sum(1 for t in s if t > TOL_ZERO * smax)


def zero(e):
    return sp.expand(sp.cancel(e)) == 0


# log variables: p1 = X'/X, m = p1^2, be = (ln X)'', c = (R'/R)^2 = c1^2, s2 = S/R, U = X R
LOGSUB = {X1: X0 * p1, X2: X0 * (be + p1**2), R1: R0 * c1, S0: R0 * s2}


def sos_form(vv):
    """s^2 - delta^2 - tau^2 (the form of E/X0^2) in the log variables."""
    U = X0 * R0
    m = p1**2
    c = c1**2
    u = U + vv
    W = 3 * U + vv
    s = U**2 * c + m * (W**2 + 1) / 2
    delta = 2 * u * U * (be + m) + m * (4 * U**2 + 1 - W**2) / 2
    tau = U * s2 + m * W
    return s**2 - delta**2 - tau**2


def static_parts():
    """v = 0: Ehat = a T_- T_+ + m T_- - (s2 + 3m)^2, and the six-term lists f, g."""
    m = p1**2
    c = c1**2
    a = X0**2 * R0**2
    Tm = c + 5 * m - 2 * be
    Tp = c + 4 * m + 2 * be
    Ehat = a * Tm * Tp + m * Tm - (s2 + 3 * m) ** 2
    th = (19 if MUT == "m-theta" else 20) * m**2 + 2 * m * be - 4 * be**2
    ps = -4 * m**2 - 2 * m * be
    f = [X0**2, X0**2 * m, X0**2 * th, m, ps, sp.Integer(1)]
    g = [R0**2 * c**2, 9 * R0**2 * c, R0**2, c - 6 * s2, sp.Integer(1), -(s2**2)]
    return Ehat, f, g


# ---------------------------------------------------------------------------------------------------------------
# C1 e_from_d
# ---------------------------------------------------------------------------------------------------------------
def c1_e_from_d():
    ysign = -1 if MUT == "m-d-sign" else 1
    E = E_mult(v, ysign)
    i_ok = zero(E - dc.E_mult_poly(X0, X1, X2, R0, R1, S0, v))
    sub = E.subs(LOGSUB, simultaneous=True)
    ii_ok = zero(sp.cancel(sub / X0**2) - sos_form(v))
    Ehat = static_parts()[0]
    iii_ok = zero(sub.subs(v, 0) - R0**2 * X0**4 * Ehat)
    # (iv) numerical derivatives of u(x, rho) = X(x) R(rho) + v at five points
    xx, rr = sp.symbols("xx rr")
    Xf = sp.exp(xx) + xx**2
    Rf = rr**3 + rr / 2 + 1
    Ej = sp.lambdify((X0, X1, X2, R0, R1, S0, v), E, "mpmath")
    jx = [sp.lambdify(xx, sp.diff(Xf, xx, k), "mpmath") for k in range(3)]
    jr = [sp.lambdify(rr, sp.diff(Rf, rr, k), "mpmath") for k in range(3)]
    pts = [(mp.mpf(1) / 3, mp.mpf(1)), (mp.mpf(1) / 2, mp.mpf(3) / 2), (mp.mpf(-1) / 5, mp.mpf(2)), (mp.mpf(1), mp.mpf(5) / 4),
           (mp.mpf(-1) / 2, mp.mpf(9) / 7)]
    worst = mp.mpf(0)
    for vv in (mp.mpf(0), mp.mpf(1) / 3):
        for (x0, r0) in pts:
            def uf(x, r):
                return jx[0](x) * jr[0](r) + vv
            d = lambda i, j: mp.diff(uf, (x0, r0), (i, j))
            u0, ux, ur, uxx, uxr, urr = uf(x0, r0), d(1, 0), d(0, 1), d(2, 0), d(1, 1), d(0, 2)
            Dnum = D_from_jets(u0, ux, ur, uxx, uxr, urr + ur / r0, ysign)
            Ejet = Ej(jx[0](x0), jx[1](x0), jx[2](x0), jr[0](r0), jr[1](r0), jr[2](r0) + jr[1](r0) / r0, vv)
            Dj = jr[1](r0) ** 2 * Ejet / 8
            worst = max(worst, abs(Dnum - Dj) / abs(Dj))
    iv_ok = worst <= mp.mpf(10) ** -40
    check("e_from_d", i_ok and ii_ok and iii_ok and iv_ok,
          f"E = 8D/R'^2 equals det_classify.E_mult_poly: {i_ok}; E/X^2 = s^2 - delta^2 - tau^2 (symbolic v): {ii_ok}; v = 0: "
          f"E = R^2 X^4 Ehat: {iii_ok}; numerical-derivative D vs E at 10 points, worst relative error {mp.nstr(worst, 3)}: {iv_ok}")


# ---------------------------------------------------------------------------------------------------------------
# C2 decomp_static
# ---------------------------------------------------------------------------------------------------------------
def c2_decomp_static():
    Ehat, f, g = static_parts()
    i_ok = zero(Ehat - sum(a * b for a, b in zip(f, g)))
    E0 = E_mult(0)
    M, rows, cols = bilinear_form(E0)
    rk = M.rank()
    ii_ok = rk == 6
    check("decomp_static", i_ok and ii_ok,
          f"v = 0: Ehat = sum_(i<=6) f_i g_i with f = (X^2, X^2 m, X^2 Theta, m, Psi, 1), g = (R^2 c^2, 9 R^2 c, R^2, c - 6 s2, 1, -s2^2): "
          f"{i_ok}; exact rank of E in x-monomials ({len(rows)}) x rho-monomials ({len(cols)}) = {rk}: {ii_ok}")


# ---------------------------------------------------------------------------------------------------------------
# C3 rank_lemma (Lemma L on planted data)
# ---------------------------------------------------------------------------------------------------------------
def c3_rank_lemma():
    rng = random.Random(160)
    pts = [sp.Rational(j, 7) + 1 for j in range(12)]
    pass_ok, fail_ok, detail = True, True, []
    for draw in range(5):
        while True:
            B = sp.Matrix(6, 4, lambda i, j: rng.randint(-3, 3))
            if B.rank() == 4:
                break
        N = sp.Matrix.hstack(*(B.T).nullspace())  # 6 x 2, columns orthogonal to col(B)
        if MUT == "m-rank-off":
            while True:
                w = sp.Matrix(6, 1, lambda i, j: rng.randint(-3, 3))
                if (B.T * w).norm() != 0:
                    break
            N = sp.Matrix.hstack(N[:, 0], w)
        fx = [B * sp.Matrix([1, x, x**2, x**3]) for x in pts]
        gr = [N * sp.Matrix([1, r]) for r in pts]
        ident = all((fa.T * gb)[0] == 0 for fa in fx for gb in gr)
        dF = sp.Matrix.hstack(*fx).rank()
        dG = sp.Matrix.hstack(*gr).rank()
        ok = ident and dF == 4 and dG == 2 and dF + dG <= 6
        pass_ok &= ok
        while True:
            N3 = sp.Matrix(6, 3, lambda i, j: rng.randint(-3, 3))
            if N3.rank() == 3:
                break
        gr3 = [N3 * sp.Matrix([1, r, r**2]) for r in pts]
        nonzero = any((fa.T * gb)[0] != 0 for fa in fx for gb in gr3)
        fail_ok &= nonzero and sp.Matrix.hstack(*gr3).rank() == 3
        detail.append((ident, dF, dG, nonzero))
    check("rank_lemma", pass_ok and fail_ok,
          f"5 planted draws, r = 6: sum f_i g_i = 0 with d_F = 4, d_G = 2 (d_F + d_G <= 6): {pass_ok}; with d_G = 3 (sum 7 > 6) the sum is "
          f"nonzero at a sample point in every draw: {fail_ok}; per draw (identity, d_F, d_G, 7-case nonzero) = {detail}")


# ---------------------------------------------------------------------------------------------------------------
# C4 lemma_x
# ---------------------------------------------------------------------------------------------------------------
def moebius_f(a, b, c, d, m):
    """The six x-functions (E_x, m E_x, E_x Theta, m, Psi, 1) on the Moebius family E_x = (a m + b)/(c m + d)."""
    den = c * m + d
    Dl = a * d - b * c
    rho = (a * m + b) / den
    beta = (a * m + b) * den / Dl
    th = 20 * m**2 + 2 * m * beta - 4 * beta**2
    ps = -4 * m**2 - 2 * m * beta
    f = [rho, m * rho, rho * th, m, ps, Fraction(1)]
    if MUT == "m-fdrop":
        f[2] = Fraction(0)
        f[4] = Fraction(0)
    return f


def c4_lemma_x():
    a, b, c, d, m = sp.symbols("a b c d m")
    Dl = a * d - b * c
    rho = (a * m + b) / (c * m + d)
    beta = sp.simplify(rho / sp.diff(rho, m))
    s1_ok = zero(beta - (a * m + b) * (c * m + d) / Dl)
    th = 20 * m**2 + 2 * m * beta - 4 * beta**2
    ps = -4 * m**2 - 2 * m * beta
    t = 1 / (c * m + d)
    # c != 0: span{1, m, rho, m rho} = span{1, m, t}
    r1 = zero(t - (a - c * rho) / Dl)
    r2 = zero(m * rho - (a / c * m - Dl / c**2 + Dl * d / c**2 * t))
    r3 = zero(rho - (a / c - Dl / c * t))
    ps_poly = sp.Poly(sp.cancel(ps), m)
    coef3 = sp.simplify(ps_poly.coeff_monomial(m**3))
    cf3_ok = zero(coef3 + 2 * a * c / Dl)
    coef2_a0 = sp.simplify(ps_poly.coeff_monomial(m**2).subs(a, 0))
    cf2_ok = zero(coef2_a0 + 2)
    # c = 0
    rho0 = ((a * m + b) / d)
    th0 = sp.simplify(th.subs(c, 0))
    cub = sp.Poly(sp.cancel(rho0 * th0), m).coeff_monomial(m**3)
    cub_ok = zero(cub - 18 * a / d)
    sym_ok = s1_ok and r1 and r2 and r3 and cf3_ok and cf2_ok and cub_ok
    mv = [Fraction(k, 3) + Fraction(1, 7) for k in range(40)]
    worst, nsamp, bad = 99, 0, None
    params = [(A, B, C, Dd) for A in range(-2, 4) for B in range(-2, 4) for C in range(-2, 4) for Dd in range(-2, 4)
              if A * Dd - B * C != 0]
    rng = random.Random(160)
    for _ in range(60):
        while True:
            A, B, C, Dd = [Fraction(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(4)]
            if A * Dd - B * C != 0:
                break
        params.append((A, B, C, Dd))
    for (A, B, C, Dd) in params:
        A, B, C, Dd = Fraction(A), Fraction(B), Fraction(C), Fraction(Dd)
        pts = [mm for mm in mv if C * mm + Dd != 0][:14]
        rk = frac_rank([[fi for fi in moebius_f(A, B, C, Dd, mm)] for mm in pts])
        nsamp += 1
        if rk < worst:
            worst, bad = rk, (A, B, C, Dd)
    rk_ok = worst >= 4
    check("lemma_x", sym_ok and rk_ok,
          f"Moebius step beta = (am+b)(cm+d)/(ad-bc): {s1_ok}; c != 0: t, m rho, rho in span: {r1 and r2 and r3}, Psi cubic coefficient "
          f"-2ac/(ad-bc): {cf3_ok}, at a = 0 quadratic coefficient -2: {cf2_ok}; c = 0: cubic coefficient of rho Theta is 18a/d: {cub_ok}; "
          f"smallest exact rank of the six functions over {nsamp} parameter sets = {worst} (needs >= 4; at {bad})")


# ---------------------------------------------------------------------------------------------------------------
# C5 lemma_rho
# ---------------------------------------------------------------------------------------------------------------
rho_s = sp.Symbol("rho", positive=True)


def g_funcs(R):
    """The six rho-functions g_i for a sympy expression R(rho_s)."""
    Rp = sp.diff(R, rho_s)
    S = sp.diff(R, rho_s, 2) + Rp / rho_s
    cc = Rp**2 / R**2
    s2e = S / R
    g4 = cc - 6 * s2e
    g6 = -(s2e**2)
    if MUT == "m-gdrop":
        g4, g6 = cc, sp.Integer(1)
    return [R**2 * cc**2, 9 * R**2 * cc, R**2, g4, sp.Integer(1), g6]


def c5_lemma_rho():
    y, al, ga, p, q, kap, rr = sp.symbols("y alpha gamma p q kappa rr")
    co = sp.Poly(sp.expand((al + ga * y) ** 2 - (p * y + q * y**2)), y).all_coeffs()
    sol = sp.solve(co, [al, p, q], dict=True)
    a_ok = len(sol) > 0 and all(s_.get(al) == 0 for s_ in sol)
    Rexp = sp.exp(kap * rho_s)
    s2e = (sp.diff(Rexp, rho_s, 2) + sp.diff(Rexp, rho_s) / rho_s) / Rexp
    w = sp.Symbol("w")
    pol = sp.Poly(sp.expand(sp.simplify(s2e**2).subs(rho_s, 1 / w)), w)
    k_ok = zero(pol.coeff_monomial(w**2) - kap**2)
    fam = {"exp(rho/2)": sp.exp(rho_s / 2), "exp(-rho)": sp.exp(-rho_s), "exp(2 rho)": sp.exp(2 * rho_s), "rho^(1/2)": sp.sqrt(rho_s),
           "rho": rho_s, "rho^2": rho_s**2, "rho^3": rho_s**3, "1/rho": 1 / rho_s, "rho^(5/2)": rho_s ** sp.Rational(5, 2),
           "rho^4": rho_s**4, "rho^2+1": rho_s**2 + 1, "rho+rho^3": rho_s + rho_s**3, "cosh(rho)": sp.cosh(rho_s)}
    pts = [mp.mpf(1) + 2 * mp.mpf(j) / 13 for j in range(14)]
    ranks = {}
    for name, R in fam.items():
        gl = [sp.lambdify(rho_s, gi, "mpmath") for gi in g_funcs(R)]
        ranks[name] = num_rank([[gi(pt) for pt in pts] for gi in gl])
    rk_ok = all(r >= 3 for r in ranks.values())
    check("lemma_rho", a_ok and k_ok and rk_ok,
          f"(alpha + gamma y)^2 = p y + q y^2 forces alpha = 0: {a_ok}; coefficient of 1/rho^2 in s2^2 for exp(kappa rho) is kappa^2: {k_ok}; "
          f"numerical ranks of the six g_i (need >= 3): {ranks}")


# ---------------------------------------------------------------------------------------------------------------
# C6 const_m
# ---------------------------------------------------------------------------------------------------------------
def c6_const_m():
    lam, c = sp.symbols("lam c")
    Ehat, f, g = static_parts()
    sub = {p1: lam, be: 0, c1: sp.sqrt(c)}
    coefX2 = sp.expand(sp.expand(Ehat.subs(sub)).coeff(X0, 2))
    sym_ok = zero(coefX2 - R0**2 * (c + 4 * lam**2) * (c + 5 * lam**2))
    sgn = -1 if MUT == "m-const-sign" else 1
    prod_ok = all((sgn * Fraction(n, 10) + 4 * Fraction(k, 4) ** 2) * (sgn * Fraction(n, 10) + 5 * Fraction(k, 4) ** 2) > 0
                  for n in range(61) for k in range(1, 9))
    xs = [mp.mpf(j) / 7 for j in range(10)]
    dfs = {}
    for lv in (mp.mpf(1) / 2, mp.mpf(1), mp.mpf(3) / 2):
        X = [mp.exp(lv * x) for x in xs]
        rowsf = [[Xi**2 for Xi in X], [lv**2 * Xi**2 for Xi in X], [20 * lv**4 * Xi**2 for Xi in X], [lv**2] * 10, [-4 * lv**4] * 10,
                 [mp.mpf(1)] * 10]
        dfs[str(lv)] = num_rank(rowsf)
    d_ok = all(dd == 2 for dd in dfs.values())
    check("const_m", sym_ok and prod_ok and d_ok,
          f"coefficient of X^2 in Ehat (m = lam^2, beta = 0) is R^2 (c + 4 lam^2)(c + 5 lam^2): {sym_ok}; (c + 4 lam^2)(c + 5 lam^2) > 0 for "
          f"c = {'-' if sgn < 0 else '+'}n/10 (n <= 60), lam = k/4 (k <= 8): {prod_ok}; numerical d_f at lam in {{1/2, 1, 3/2}} = {dfs}: {d_ok}")


# ---------------------------------------------------------------------------------------------------------------
# C7 mult_subcases
# ---------------------------------------------------------------------------------------------------------------
def c7_mult_subcases():
    E = E_mult(v)
    Esub = sp.expand(E.subs({X1: 0, X2: 0}))
    sq = 1 if MUT == "m-sub-drop" else 2
    Euse = X0**6 * R1**4 - X0**4 * S0**sq
    form_ok = zero(Esub - Euse)
    M, rows, cols = bilinear_form(Euse)
    rk_ok = M.rank() == 2
    C1_, C2_, x0 = sp.symbols("C1 C2 x0", positive=True)
    sols = []
    for eps in (1, -1):
        R = -(eps / x0) * sp.log(C1_ + C2_ * sp.log(rho_s))
        Rp = sp.diff(R, rho_s)
        S = sp.diff(R, rho_s, 2) + Rp / rho_s
        Es = x0**6 * Rp**4 - x0**4 * S**sq
        sols.append(sp.simplify(sp.together(Es)) == 0)
    planar_ok = zero(D_from_jets(X0 * R0 + v, X1 * R0, 0, X2 * R0, 0, 0))
    check("mult_subcases", form_ok and rk_ok and all(sols) and planar_ok,
          f"X1 = X2 = 0: E = X0^6 R1^4 - X0^4 S0^2 (rank 2 = {rk_ok}): {form_ok}; E = 0 on R = -(eps/X0) ln(C1 + C2 ln rho), eps = +1, -1: "
          f"{sols}; planar R1 = S0 = 0: D = 0: {planar_ok}")


# ---------------------------------------------------------------------------------------------------------------
# C8 add_decomp, C9 add_stratum
# ---------------------------------------------------------------------------------------------------------------
def add_fg():
    A = R1**2 / 2
    B = S0**2 / 4
    f = [sp.Integer(1), X1**2, X2**2, X0 * X2**2, X1**4 - X0**2 * X2**2]
    g = [A**2 - B, 2 * A, -(R0**2), -2 * R0, sp.Integer(1)]
    return f, g


def c8_add_decomp():
    u = X0 if MUT == "m-add-drop" else X0 + R0
    D = D_from_jets(u, X1, R1, X2, 0, S0)
    Eadd = sp.expand(sp.cancel(D / (R1**2 / 2)))
    i_ok = zero(Eadd - dc.E_add_poly(X0, X1, X2, R0, R1, S0))
    f, g = add_fg()
    ii_ok = zero(Eadd - sum(a * b for a, b in zip(f, g)))
    rk = mat_rank(Eadd)
    check("add_decomp", i_ok and ii_ok and rk == 5,
          f"E_add = D/h equals det_classify.E_add_poly: {i_ok}; E_add = sum_(i<=5) f_i g_i: {ii_ok}; exact rank {rk} (needs 5)")


def c9_add_stratum():
    f, g = add_fg()
    ks = sp.symbols("k2 k3 k4 k5")
    use = [1, 2, 3, 4] if MUT != "m-stratum" else [1, 2]
    eqs = [f[i] - ks[i - 1] for i in use]
    gens = (X1, X2, X0) + ks
    G = sp.groebner(eqs, *gens, order="lex")
    elim = [e for e in G.exprs if not (e.has(X1) or e.has(X2))]
    cond = []
    for e in elim:
        cond += sp.Poly(e, X0).all_coeffs() if e.has(X0) else [e]
    cond = [cc for cc in cond if cc != 0]
    G2 = sp.groebner(eqs + cond, *gens, order="lex")
    x2_ok = G2.contains(X2**2)
    k5_ok = G2.contains(ks[3] - ks[0] ** 2)
    free_ok = not any(e.free_symbols <= {ks[0]} and e != 0 for e in G2.exprs)
    check("add_stratum", x2_ok and k5_ok and free_ok,
          f"elimination of {{f_i(a) - k_i}} with X0 free: conditions on k {cond}; X2^2 in the ideal (so X2 = 0): {x2_ok}; k5 - k2^2 in the ideal: {k5_ok}; "
          f"k2 (= c^2) unconstrained: {free_ok}")


# ---------------------------------------------------------------------------------------------------------------
# C10-C12 planted Bessel family
# ---------------------------------------------------------------------------------------------------------------
BESSEL_CASES = [(sp.Rational(n, d), e) for (n, d) in ((3, 10), (1, 2), (1, 1), (3, 2), (2, 1)) for e in (1, -1)]
ZPTS = [mp.mpf(2) / 5 * (1 + j) for j in range(5)]  # 0.4 .. 2.0
ZRANK = [mp.mpf(2) / 5 + mp.mpf(8) / 5 * j / 11 for j in range(12)]


def bessel_R(c, eps, pert=0, wrong_freq=False):
    fac = 1 if wrong_freq else sp.sqrt(2)
    z = fac * c * rho_s
    psi = sp.besselj(0, z) + sp.Rational(3, 10) * sp.bessely(0, z)
    return -eps * sp.log(psi) + pert


def g_add_numeric(R, c, rhos):
    Rp = sp.diff(R, rho_s)
    S = sp.diff(R, rho_s, 2) + Rp / rho_s
    A, B = Rp**2 / 2, S**2 / 4
    gl = [sp.lambdify(rho_s, gi, "mpmath") for gi in (A**2 - B, 2 * A, -(R**2), -2 * R, sp.Integer(1))]
    return [[gi(r) for r in rhos] for gi in gl]


def rho_from_z(c, zs):
    return [z / (mp.sqrt(2) * mp.mpf(c.p) / c.q) for z in zs]


def c10_add_planted():
    worst, dgs, amp_ok, dFs = mp.mpf(0), [], True, []
    for (c, eps) in BESSEL_CASES:
        cm = mp.mpf(c.p) / c.q
        cf = Fraction(c.p, c.q)
        dFs.append(frac_rank([[Fraction(1)] * 5, [cf**2] * 5, [Fraction(0)] * 5, [Fraction(0)] * 5, [cf**4] * 5]))
        rhos = rho_from_z(c, ZPTS)
        gm = g_add_numeric(bessel_R(c, eps, 0, MUT == "m-add-freq"), c, rhos)
        f0 = [1, cm**2, 0, 0, cm**4]
        for j in range(5):
            terms = [f0[i] * gm[i][j] for i in range(5)]
            worst = max(worst, abs(sum(terms)) / sum(abs(t) for t in terms))
        rr = rho_from_z(c, ZRANK)
        dgs.append(num_rank(g_add_numeric(bessel_R(c, eps, 0, MUT == "m-add-freq"), c, rr)))
        amp_ok &= c >= sp.Rational(3, 10)
    dF = max(dFs)
    ok = worst <= mp.mpf(10) ** -50 and max(dgs) <= 4 and dF == 1 and amp_ok
    check("add_planted", ok,
          f"u = c x - eps ln(J0(z) + 0.3 Y0(z)) + d, c in {{3/10, 1/2, 1, 3/2, 2}}, eps = +-1, 5 z-points each: largest relative residual of "
          f"f0.g {mp.nstr(worst, 3)} (needs <= 1e-50); d_F = {dF}; numerical d_G per case {dgs} (needs <= 4); amplitude c >= 3/10: {amp_ok}")


def c11_add_perturbed():
    amp = 0 if MUT == "m-perturb-zero" else 1
    smallest = mp.mpf(10) ** 9
    for (c, eps) in BESSEL_CASES:
        cm = mp.mpf(c.p) / c.q
        rhos = rho_from_z(c, ZPTS)
        for dl in (rho_s**2 / 10, rho_s**3 / 10):
            gm = g_add_numeric(bessel_R(c, eps, amp * dl), c, rhos)
            f0 = [1, cm**2, 0, 0, cm**4]
            for j in range(5):
                terms = [f0[i] * gm[i][j] for i in range(5)]
                smallest = min(smallest, abs(sum(terms)) / sum(abs(t) for t in terms))
    check("add_perturbed", smallest >= mp.mpf(10) ** -3,
          f"R + delta, delta in {{rho^2/10, rho^3/10}}: smallest relative residual over 100 points {mp.nstr(smallest, 3)} (needs >= 1e-3)")


def c12_add_xperturbed():
    amp = 0 if MUT == "m-xperturb-zero" else 1
    dFs, dGs = [], []
    for (c, eps) in BESSEL_CASES:
        xs = [Fraction(j, 4) for j in range(10)]
        cf = Fraction(c.p, c.q)
        rows = []
        for x in xs:
            X0v = cf * x + amp * x**2 / 10
            X1v = cf + amp * x / 5
            X2v = amp * Fraction(1, 5)
            rows.append([Fraction(1), X1v**2, X2v**2, X0v * X2v**2, X1v**4 - X0v**2 * X2v**2])
        dFs.append(frac_rank(rows))
        rr = rho_from_z(c, ZRANK)
        dGs.append(num_rank(g_add_numeric(bessel_R(c, eps), c, rr)))
    ok = all(a >= 2 for a in dFs) and all(b == 4 for b in dGs)
    check("add_xperturbed", ok,
          f"X = c x + x^2/10: exact d_F per case {dFs} (needs >= 2); numerical d_G {dGs} (needs 4); d_F + d_G > 5 so Lemma L rejects it: {ok}")


# ---------------------------------------------------------------------------------------------------------------
# C13 trav_decomp, C14 trav_gap
# ---------------------------------------------------------------------------------------------------------------
def trav_FG():
    m = p1**2
    c = c1**2
    Ph2 = -4 * v**2 * be**2 + (6 * v**2 - 2) * be * m + (12 * v**2 - 4) * m**2
    Ph3 = -4 * be**2 + 3 * be * m + 13 * m**2
    Th = 20 * m**2 + 2 * m * be - 4 * be**2
    F = [-2 * v * (1 - v**2) * X0 * m * (be + m), -2 * v * X0 * m, X0**2 * Ph2, X0**2 * m, X0**2, 2 * v * X0**3 * Ph3, 6 * v * X0**3 * m,
         X0**4, X0**4 * m, X0**4 * Th]
    G = [R0, R0 * s2, R0**2, R0**2 * ((v**2 + 1) * c - 6 * s2), -(R0**2) * s2**2, R0**3, R0**3 * c, R0**4 * c**2, 9 * R0**4 * c, R0**4]
    return F, G


def c13_trav_decomp():
    U = sp.Symbol("U")
    m = p1**2
    c = c1**2
    E1 = -2 * m * v * ((1 - v**2) * (be + m) + s2)
    E2 = (-4 * be**2 * v**2 + 6 * be * m * v**2 - 2 * be * m + c * m * (v**2 + 1) + 12 * m**2 * v**2 - 4 * m**2 - 6 * m * s2 - s2**2)
    E3 = 2 * v * (-4 * be**2 + 3 * be * m + 3 * c * m + 13 * m**2)
    E4 = (c + 5 * m - 2 * be) * (c + 4 * m + 2 * be)
    sos = sos_form(v)
    i_ok = zero(sos - sum((X0 * R0) ** k * Ek for k, Ek in enumerate((E1, E2, E3, E4), start=1)))
    F, G = trav_FG()
    ii_ok = zero(sos - sum(a * b for a, b in zip(F, G)))
    E = E_mult(v)
    ranks = {}
    want = 9 if MUT == "m-trav-rank" else 10
    for vv in (sp.Rational(1, 3), 2, sp.Rational(5, 7), sp.Rational(3, 2), -2):
        ranks[str(vv)] = mat_rank(E.subs(v, vv))
    iii_ok = all(r == want for r in ranks.values())
    obs = {str(vv): mat_rank(E.subs(v, vv)) for vv in (0, 1, -1)}
    check("trav_decomp", i_ok and ii_ok and iii_ok,
          f"E/X0^2 = sum_(k=1..4) U^k E_k: {i_ok}; ten-term separation F.G = E/X0^2: {ii_ok}; exact ranks of E at v = 1/3, 2, 5/7, 3/2, -2: "
          f"{ranks} (needs {want}). Observation, rank at v = 0, 1, -1: {obs}")


def c14_trav_gap():
    R = sp.exp(rho_s) if MUT == "m-gap-d" else rho_s**2
    F, G = trav_FG()
    vv = sp.Rational(1, 3)
    Rp = sp.diff(R, rho_s)
    S = sp.diff(R, rho_s, 2) + Rp / rho_s
    sub = {R0: R, R1: Rp, S0: S, c1: Rp / R, s2: S / R, v: vv}
    gl = [sp.lambdify(rho_s, sp.simplify(Gi.subs(sub)), "mpmath") for Gi in G]
    pts = [mp.mpf(1) + mp.mpf(j) / 4 for j in range(12)]
    rk = num_rank([[gi(pt) for pt in pts] for gi in gl])
    check("trav_gap", rk == 5,
          f"v = 1/3, R = {'exp(rho)' if MUT == 'm-gap-d' else 'rho^2'}: numerical d_G of the ten rho-functions = {rk} (needs 5): the rank bound "
          f"d_F + d_G <= 10 alone does not exclude the traveling case")


CHECKS = [("e_from_d", c1_e_from_d), ("decomp_static", c2_decomp_static), ("rank_lemma", c3_rank_lemma), ("lemma_x", c4_lemma_x),
          ("lemma_rho", c5_lemma_rho), ("const_m", c6_const_m), ("mult_subcases", c7_mult_subcases), ("add_decomp", c8_add_decomp),
          ("add_stratum", c9_add_stratum), ("add_planted", c10_add_planted), ("add_perturbed", c11_add_perturbed),
          ("add_xperturbed", c12_add_xperturbed), ("trav_decomp", c13_trav_decomp), ("trav_gap", c14_trav_gap)]


def verdict():
    ok = {n: s_ for n, s_ in RESULTS}
    stat = all(ok.get(n, False) for n in ("e_from_d", "decomp_static", "rank_lemma", "lemma_x", "lemma_rho", "const_m"))
    plant = all(ok.get(n, False) for n in ("add_decomp", "add_stratum", "add_planted", "add_perturbed", "add_xperturbed", "mult_subcases"))
    s_yes = ("no solution with X X' R R' != 0 (d_F >= 4 unless m is constant, d_G >= 3, 4 + 3 > 6; m constant: c >= 0 against a positive product)")
    print("\nverdicts:")
    print("  static v = 0: " + (s_yes if stat else "NOT established (a prerequisite check failed)"))
    print("  planted additive truth through the same code: " + ("found, perturbations rejected" if plant else "NOT established"))
    print("  traveling v != 0: not decided here (rank 10, d_F + d_G <= 10, the bound alone does not exclude it)")


def main():
    t0 = time.time()
    if MUT is not None and MUT not in MUTATIONS:
        print(f"unknown mutation {MUT}")
        sys.exit(2)
    for name, fn in CHECKS:
        if wanted(name):
            t1 = time.time()
            fn()
            print(f"  ({time.time() - t1:.1f} s)", flush=True)
    failed = [n for n, ok in RESULTS if not ok]
    verdict()
    print(f"\nchecks: {len(RESULTS) - len(failed)} passed, {len(failed)} failed {failed}. {time.time() - t0:.0f} s")
    if MUT is not None:
        print(f"control {MUT}: target check {MUTATIONS[MUT][0]}, fired {MUTATIONS[MUT][0] in failed}")
    return failed


if __name__ == "__main__":
    failed = main()
    controls_ok = True
    if MUT is None and "no-controls" not in OPTS and ONLY is None:
        print("\n== Controls: each mutation must make exactly its target check fail ==", flush=True)
        for name, (target, what) in MUTATIONS.items():
            p = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True)
            if "logdir" in OPTS:
                with open(os.path.join(OPTS["logdir"], f"ctl_{name}.log"), "w") as fh:
                    fh.write(p.stdout + p.stderr)
            fails = [ln.split(":")[0][7:] for ln in p.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = fails == [target] and p.returncode == 1
            controls_ok &= fired
            print(f"[{'FIRED' if fired else 'SILENT' if target not in fails else 'EXTRA'}] --mutate={name} ({what}): exit {p.returncode}, failed {fails}",
                  flush=True)
            if p.returncode not in (0, 1):
                print(p.stderr[-2000:])
        print(f"controls: {'all fired on their target only' if controls_ok else 'SOME DID NOT FIRE ON TARGET ONLY'}")
    sys.exit(0 if not failed and controls_ok else 1)

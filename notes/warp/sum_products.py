#!/usr/bin/env python3
"""Sums of two products u = X1(x) R1(rho) + X2(x) R2(rho) + v with det R_3 = 0 (issue #169).

D = det R_3 of the static shift u (det_classify.md Setup) is retyped here (shared with #80 and #160):
  D = h (h^2 - P^2 - Q^2 + Y^2 + Z^2) - P (Y^2 - Z^2) + 2 Q Y Z,
  P = -(u_x^2 + u u_xx), Q = -(u_rho_rho + u_rho/rho)/2, Z = u_x_rho/2, Y = u_x u_rho + u u_x_rho/2, h = u_rho^2/2.
For two products the jets are A_k = X1^(k), B_k = X2^(k) (k = 0, 1, 2) and a_0, a_1, a_2 = R1, R1', R1'' + R1'/rho, b_k likewise
for R2, and E = 8 D is a polynomial in these 12 jets and v. The script
  e_two      builds E, checks it against derivatives of undefined functions, against det_classify.E_mult_poly (X2 = 0) and
             det_classify.E_add_poly (X2 = 1, R1 = 1), and proves the compact form 4 h D = C^2 - A^2 - B^2;
  ode_form   with W_x, W_rho != 0: a'' = alpha a + beta a', s = Gamma r' + delta r, D in the 2 x 2 matrix N and these four functions;
  case_a     the excluded forms (1 in V_X and 1 in V_R) are one product plus a constant or additive;
  sep_rank   exact separation ranks of E (Lemma L bound d_F + d_G <= r) and of its sub-forms;
  slack      Lemma L does not decide: an exponential pair on both sides satisfies the bound and is not a solution;
  rho_const  rho-side with constant ODE coefficients has no solution (coefficient of iota^2, iota = 1/rho, is -u_rho^4);
  add_*      planted truth: the additive family u = c x + g(rho) of det_classify.md section 3 found through the same E.
Imported from det_classify: E_mult_poly, E_add_poly (comparison targets only).
Usage: sum_products.py [--no-controls] [--mutate=NAME] [--only=a,b] [--logdir=DIR] [--irreducible].
Exit 1 if a check fails or a control does not fire on its target only."""
import os
import subprocess
import sys
import time
from fractions import Fraction

import mpmath as mp
import sympy as sp
from sympy import QQ
from sympy.polys.matrices import DomainMatrix

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import det_classify as dc  # noqa: E402  (E_mult_poly, E_add_poly only)

# name: (target check, what the mutation does)
MUTATIONS = {
    "s-d-sign": ("e_two", "Y := u_x u_rho - u u_x_rho/2 in the jets-built E"),
    "s-compact-b": ("e_two", "B := 2hQ + 2YZ in the compact form"),
    "s-ode-alpha": ("ode_form", "alpha with a'' ^ a' replaced by a' ^ a''"),
    "s-casea-const": ("case_a", "the constant a - bc/d replaced by a"),
    "s-rank-qyz": ("sep_rank", "the term 2QYZ removed from D before the ranks"),
    "s-slack-kappa": ("slack", "the second exponent kappa = 2 (resonant: z2 = z1^2) in place of a transcendental kappa"),
    "s-rhoconst-iota": ("rho_const", "S_i = R_i'' without the R_i'/rho term"),
    "s-rhoconst-double": ("rho_const", "R2' = mu rho e^(mu rho) in place of (1 + mu rho) e^(mu rho)"),
    "s-add-freq": ("add_planted", "J0(c rho), Y0(c rho) in place of J0(sqrt2 c rho), Y0(sqrt2 c rho)"),
    "s-add-exact": ("add_planted_exact", "psi'' = -psi'/rho + 2 c^2 psi in place of the Bessel equation"),
    "s-stratum-drop": ("add_stratum", "two of the five f_k dropped from the elimination"),
    "s-pexact-zero": ("add_perturbed_exact", "delta = 0"),
    "s-xperturb-zero": ("add_xperturbed", "X = c x (no x^2/10)"),
}
# declared silent in README A1 after dev run 3: run once, reported, not counted among the fired controls
DECLARED_SILENT = {"s-rank-drop": ("sep_rank", "the term 8 h Z^2 removed from E before the ranks (changes the column counts, not the ranks)")}
RESULTS = []
OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate")
ONLY = OPTS["only"].split(",") if "only" in OPTS else None
mp.mp.dps = 130

A0, A1, A2, B0, B1, B2, a0, a1, a2, b0, b1, b2, v, io = sp.symbols("A0 A1 A2 B0 B1 B2 a0 a1 a2 b0 b1 b2 v io")
XV = (A0, A1, A2, B0, B1, B2)
RV = (a0, a1, a2, b0, b1, b2)


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def zero(e):
    return sp.expand(sp.cancel(sp.together(e))) == 0


# ---------------------------------------------------------------------------------------------------------------
# D, E and the bilinear form
# ---------------------------------------------------------------------------------------------------------------
def blocks(u, ux, ur, uxx, uxr, lap, ysign=1):
    """P, Q, Z, Y, h of #65's block."""
    return (-(ux**2 + u * uxx), -lap / 2, uxr / 2, ux * ur + ysign * u * uxr / 2, ur**2 / 2)


def D_of(P, Q, Z, Y, h):
    return h * (h**2 - P**2 - Q**2 + Y**2 + Z**2) - P * (Y**2 - Z**2) + 2 * Q * Y * Z


def two_jets(vv):
    """(u, u_x, u_rho, u_xx, u_x_rho, u_rho_rho + u_rho/rho) for u = X1 R1 + X2 R2 + vv in the jets."""
    return (A0 * a0 + B0 * b0 + vv, A1 * a0 + B1 * b0, A0 * a1 + B0 * b1, A2 * a0 + B2 * b0, A1 * a1 + B1 * b1, A0 * a2 + B0 * b2)


def E_two(vv=v, ysign=1):
    return sp.expand(8 * D_of(*blocks(*two_jets(vv), ysign=ysign)))


E_SYM = E_two()


def bil(E, xv, rv):
    """Coefficient matrix of E in (x-monomials) x (rho-monomials): (M, row monomials, column monomials)."""
    P = sp.Poly(sp.expand(E), *xv, *rv)
    nx = len(xv)
    rows, cols, ent = {}, {}, {}
    for mon, co in P.terms():
        ri = rows.setdefault(mon[:nx], len(rows))
        ci = cols.setdefault(mon[nx:], len(cols))
        ent[(ri, ci)] = co
    M = sp.zeros(len(rows), len(cols))
    for (i, j), co in ent.items():
        M[i, j] = co
    return M, list(rows), list(cols)


def qq_rank(M):
    return DomainMatrix.from_Matrix(M).convert_to(QQ).rank()


def bil_rank(E, xv, rv):
    M, rows, cols = bil(E, xv, rv)
    return len(rows), len(cols), qq_rank(M)


def mon_value(mon, vals):
    out = sp.Integer(1)
    for e, s in zip(mon, vals):
        out *= s**e
    return out


# ---------------------------------------------------------------------------------------------------------------
# e_two
# ---------------------------------------------------------------------------------------------------------------
def c_e_two():
    ysign = -1 if MUT == "s-d-sign" else 1
    E = E_two(v, ysign)
    # (i) against derivatives of undefined functions
    x_, r_ = sp.symbols("x rho", positive=True)
    X1f, X2f = sp.Function("X1")(x_), sp.Function("X2")(x_)
    R1f, R2f = sp.Function("R1")(r_), sp.Function("R2")(r_)
    uf = X1f * R1f + X2f * R2f + v
    ur = sp.diff(uf, r_)
    Dd = D_of(*blocks(uf, sp.diff(uf, x_), ur, sp.diff(uf, x_, 2), sp.diff(uf, x_, r_), sp.diff(uf, r_, 2) + ur / r_))
    s1 = {R1f.diff(r_, 2): a2 - a1 / r_, R2f.diff(r_, 2): b2 - b1 / r_, X1f.diff(x_, 2): A2, X2f.diff(x_, 2): B2}
    s2 = {X1f.diff(x_): A1, X2f.diff(x_): B1, R1f.diff(r_): a1, R2f.diff(r_): b1}
    s3 = {X1f: A0, X2f: B0, R1f: a0, R2f: b0}
    i_ok = zero(8 * Dd.subs(s1).subs(s2).subs(s3) - E)
    # (ii) compact form 4 h D = C^2 - A^2 - B^2, abstract blocks, then the two-product jets
    P_, Q_, Z_, Y_, h_ = sp.symbols("P_ Q_ Z_ Y_ h_")
    Aq = 2 * h_ * P_ + Y_**2 - Z_**2
    Bq = 2 * h_ * Q_ + (2 if MUT == "s-compact-b" else -2) * Y_ * Z_
    Cq = Y_**2 + Z_**2 + 2 * h_**2
    ii_ok = zero(4 * h_ * D_of(P_, Q_, Z_, Y_, h_) - (Cq**2 - Aq**2 - Bq**2))
    # (iii) X2 = 0: E = R1'^2 E_mult;  (iv) X2 = 1, R1 = 1: E = 4 R2'^2 E_add(X + v, ...)
    iii_ok = zero(E.subs({B0: 0, B1: 0, B2: 0}) - a1**2 * dc.E_mult_poly(A0, A1, A2, a0, a1, a2, v))
    iv_ok = zero(E.subs({B0: 1, B1: 0, B2: 0, a0: 1, a1: 0, a2: 0}) - 4 * b1**2 * dc.E_add_poly(A0 + v, A1, A2, b0, b1, b2))
    nterms = len(sp.Poly(E, *XV, *RV, v).terms())
    check("e_two", i_ok and ii_ok and iii_ok and iv_ok,
          f"E = 8D has {nterms} terms; E equals 8D from derivatives of undefined X1, X2, R1, R2: {i_ok}; 4hD = C^2 - A^2 - B^2 with A = 2hP + Y^2 - Z^2, "
          f"B = 2hQ - 2YZ, C = Y^2 + Z^2 + 2h^2 (abstract blocks): {ii_ok}; X2 = 0: E = R1'^2 det_classify.E_mult_poly: {iii_ok}; "
          f"X2 = 1, R1 = 1: E = 4 R2'^2 det_classify.E_add_poly(X + v, ...): {iv_ok}")


# ---------------------------------------------------------------------------------------------------------------
# ode_form
# ---------------------------------------------------------------------------------------------------------------
def c_ode_form():
    W = A0 * B1 - B0 * A1
    Wr = a0 * b1 - b0 * a1
    sg = -1 if MUT == "s-ode-alpha" else 1
    alpha = sg * (A2 * B1 - B2 * A1) / W
    beta = (A0 * B2 - B0 * A2) / W
    Ga = (a0 * b2 - b0 * a2) / Wr
    de = (a2 * b1 - b2 * a1) / Wr
    i_ok = (zero(A2 - (alpha * A0 + beta * A1)) and zero(B2 - (alpha * B0 + beta * B1)) and zero(a2 - (Ga * a1 + de * a0))
            and zero(b2 - (Ga * b1 + de * b0)))
    # D in N = [[u0, u_rho], [u_x, u_x_rho]] and four coefficient functions
    al, be, gm, dl = sp.symbols("al be gm dl")
    Es = E_SYM.subs({A2: al * A0 + be * A1, B2: al * B0 + be * B1, a2: gm * a1 + dl * a0, b2: gm * b1 + dl * b0}, simultaneous=True)
    n0, n1, n2, n3 = A0 * a0 + B0 * b0, A1 * a0 + B1 * b0, A0 * a1 + B0 * b1, A1 * a1 + B1 * b1
    Pn = -(n1**2 + (n0 + v) * (al * n0 + be * n1))
    Q2 = -(gm * n2 + dl * n0)
    Z2 = n3
    Y2 = 2 * n1 * n2 + (n0 + v) * n3
    h2 = n2**2
    Dt = h2**3 - 4 * h2 * Pn**2 - h2 * Q2**2 + h2 * Y2**2 + h2 * Z2**2 - 2 * Pn * (Y2**2 - Z2**2) + 2 * Q2 * Y2 * Z2
    ii_ok = zero(Es - Dt)
    # det N = W_x W_rho
    det_ok = zero(n0 * n3 - n1 * n2 - W * Wr)
    # X1' X2'' - X2' X1'' = -alpha W_x, and a constant in V_X forces alpha = 0
    p, q = sp.symbols("p q")
    iii_ok = zero(A1 * B2 - B1 * A2 + alpha * W) and zero(p * A2 + q * B2 - alpha * (p * A0 + q * B0) - beta * (p * A1 + q * B1))
    check("ode_form", i_ok and ii_ok and det_ok and iii_ok,
          f"a'' = alpha a + beta a', s = Gamma r' + delta r (Cramer, wedge formulas): {i_ok}; E(jets with a'', s eliminated) = D~(N; alpha, beta, Gamma, delta): "
          f"{ii_ok}; det N = W_x W_rho: {det_ok}; X1'X2'' - X2'X1'' = -alpha W_x and p X1'' + q X2'' = alpha (p X1 + q X2) + beta (p X1' + q X2'): {iii_ok}")


# ---------------------------------------------------------------------------------------------------------------
# case_a
# ---------------------------------------------------------------------------------------------------------------
def c_case_a():
    X, R, a_, b_, c_, d_ = sp.symbols("X R a_ b_ c_ d_")
    u = a_ + b_ * X + c_ * R + d_ * X * R
    const = a_ if MUT == "s-casea-const" else a_ - b_ * c_ / d_
    i_ok = zero(u - (d_ * (X + c_ / d_) * (R + b_ / d_) + const))
    ii_ok = zero(u.subs(d_, 0) - (a_ + b_ * X + c_ * R))
    # span{1, X} x span{1, R} as two products: X1 = X, X2 = 1, R1 = b + d R, R2 = a + c R; alpha = 0 for X2 = 1, delta = 0 for R2 = 1
    iii_ok = zero((A2 * B1 - B2 * A1).subs({B0: 1, B1: 0, B2: 0})) and zero((a2 * b1 - b2 * a1).subs({b0: 1, b1: 0, b2: 0}))
    iv_ok = zero(X * (b_ + d_ * R) + 1 * (a_ + c_ * R) - u)
    check("case_a", i_ok and ii_ok and iii_ok and iv_ok,
          f"u = a + bX + cR + dXR = d (X + c/d)(R + b/d) + (a - bc/d): {i_ok}; d = 0 is additive: {ii_ok}; X2 = 1 gives alpha = 0 and R2 = 1 gives "
          f"delta = 0: {iii_ok}; u = X (b + dR) + 1 (a + cR): {iv_ok}")


# ---------------------------------------------------------------------------------------------------------------
# sep_rank
# ---------------------------------------------------------------------------------------------------------------
RANKS = {  # declared after dev run 2 (README)
    "full": {"0": 59, "1/3": 93, "1": 87, "v": 102},
    "one product (X2 = 0)": {"0": 6, "1/3": 10, "1": 9, "v": 11},
    "additive (X2 = 1, R1 = 1)": {"0": 5, "1/3": 5, "1": 5, "v": 5},
    "(b) X2 = 1": {"0": 21, "1/3": 21, "1": 21, "v": 21},
    "(c) R2 = 1": {"0": 12, "1/3": 12, "1": 12, "v": 24},
}


def sub_forms(E):
    return {
        "full": (E, XV, RV),
        "one product (X2 = 0)": (E.subs({B0: 0, B1: 0, B2: 0}), (A0, A1, A2), (a0, a1, a2)),
        "additive (X2 = 1, R1 = 1)": (E.subs({B0: 1, B1: 0, B2: 0, a0: 1, a1: 0, a2: 0}), (A0, A1, A2), (b0, b1, b2)),
        "(b) X2 = 1": (E.subs({B0: 1, B1: 0, B2: 0}), (A0, A1, A2), (a0, a1, a2, b0, b1, b2)),
        "(c) R2 = 1": (E.subs({b0: 1, b1: 0, b2: 0}), XV, (a0, a1, a2)),
    }


def c_sep_rank():
    E0 = E_SYM
    if MUT in ("s-rank-drop", "s-rank-qyz"):
        u, ux, ur, uxx, uxr, lap = two_jets(v)
        P_, Q_, Z_, Y_, h_ = blocks(u, ux, ur, uxx, uxr, lap)
        E0 = sp.expand(E_SYM - (8 * h_ * Z_**2 if MUT == "s-rank-drop" else 16 * Q_ * Y_ * Z_))
    ok, lines = True, []
    for key, expect in (("0", sp.Integer(0)), ("1/3", sp.Rational(1, 3)), ("1", sp.Integer(1)), ("v", None)):
        E = E0 if expect is None else sp.expand(E0.subs(v, expect))
        for name, (Ef, xv, rv) in sub_forms(E).items():
            rvv = rv + ((v,) if expect is None else ())
            nr, nc, rk = bil_rank(Ef, xv, rvv)
            want = RANKS[name][key]
            ok &= rk == want
            lines.append(f"{name} v={key}: {nr}x{nc} rank {rk}" + ("" if rk == want else f" (needs {want})"))
    check("sep_rank", ok, "exact ranks r of E as sum_(i<=r) F_i(x) G_i(rho), x-monomials x rho-monomials, so d_F + d_G <= r (Lemma L): " + "; ".join(lines))


# ---------------------------------------------------------------------------------------------------------------
# slack: Lemma L does not decide
# ---------------------------------------------------------------------------------------------------------------
SLACK = {"0": (12, 22, 59), "1/3": (22, 36, 93)}  # d_F, d_G of the exponential pair, r (declared after dev run 2)


def exp_pair_ranks(vv):
    kap, mu, z1, z2, w1, w2 = sp.symbols("kap mu z1 z2 w1 w2")
    E = sp.expand(E_SYM.subs(v, vv))
    M, rows, cols = bil(E, XV, RV)
    r2 = qq_rank(M)
    resonant = MUT == "s-slack-kappa"
    if resonant:
        xsub = [z1, z1, z1, z1**2, 2 * z1**2, 4 * z1**2]
        xg = (z1,)
    else:
        xsub = [z1, z1, z1, z2, kap * z2, kap**2 * z2]
        xg = (z1, z2)
    rsub = [w1, w1, (1 + io) * w1, w2, mu * w2, (mu**2 + io * mu) * w2]
    Fm = [mon_value(m, xsub) for m in rows]
    Gm = [mon_value(m, rsub) for m in cols]
    Phi, Psi = {}, {}
    for (i, j), co in [((i, j), M[i, j]) for i in range(M.rows) for j in range(M.cols) if M[i, j] != 0]:
        Phi[j] = Phi.get(j, 0) + co * Fm[i]
        Psi[i] = Psi.get(i, 0) + co * Gm[j]

    def span_dim(funcs, gens, par):
        polys = [sp.Poly(sp.expand(f), *gens) for f in funcs.values()]
        mons = {}
        for P in polys:
            for m in P.monoms():
                mons.setdefault(m, len(mons))
        T = sp.zeros(len(polys), len(mons))
        for r_, P in enumerate(polys):
            for m, c in P.terms():
                T[r_, mons[m]] = c
        dom = QQ.frac_field(par) if T.has(par) else QQ
        return DomainMatrix.from_Matrix(T).convert_to(dom).rank()

    dF = span_dim(Phi, xg, kap)
    dG = span_dim(Psi, (w1, w2, io), mu)
    # the pair is not a solution: E on the pair is a nonzero exponential polynomial
    jets = {A0: xsub[0], A1: xsub[1], A2: xsub[2], B0: xsub[3], B1: xsub[4], B2: xsub[5], a0: rsub[0], a1: rsub[1], a2: rsub[2],
            b0: rsub[3], b1: rsub[4], b2: rsub[5]}
    nonzero = sp.expand(E.subs(jets, simultaneous=True)) != 0
    return dF, dG, r2, nonzero


def c_slack():
    ok, lines = True, []
    for key, vv in (("0", 0), ("1/3", sp.Rational(1, 3))):
        dF, dG, r2, nz = exp_pair_ranks(vv)
        want = SLACK[key]
        good = (dF, dG, r2) == want and dF + dG <= r2 and nz
        ok &= good
        lines.append(f"v={key}: d_F = {dF}, d_G = {dG}, sum {dF + dG} <= r = {r2}, E nonzero on the pair: {nz}" + ("" if good else f" (needs {want})"))
    check("slack", ok, "(X1, X2) = (e^x, e^(kappa x)), (R1, R2) = (e^rho, e^(mu rho)), kappa and mu transcendental (exact ranks over Q(kappa), Q(mu)): "
          + "; ".join(lines) + "; the bound d_F + d_G <= r holds with room although E != 0")


# ---------------------------------------------------------------------------------------------------------------
# rho_const: Theorem K
# ---------------------------------------------------------------------------------------------------------------
def c_rho_const():
    mu1, mu2, mu, w1, w2, w, rho = sp.symbols("mu1 mu2 mu w1 w2 w rho")
    ion = 0 if MUT == "s-rhoconst-iota" else io
    # (i) simple roots: r = (w1, w2) = (e^(mu1 rho), e^(mu2 rho)), S_i = (mu_i^2 + iota mu_i) w_i
    rs = {a0: w1, a1: mu1 * w1, a2: (mu1**2 + ion * mu1) * w1, b0: w2, b1: mu2 * w2, b2: (mu2**2 + ion * mu2) * w2}
    Es = sp.expand(E_SYM.subs(rs, simultaneous=True))
    Pi = sp.Poly(Es, io)
    n2 = A0 * mu1 * w1 + B0 * mu2 * w2
    c2 = Pi.coeff_monomial(io**2)
    i_ok = Pi.degree() == 2 and zero(c2 + n2**4)
    w14 = sp.Poly(c2, w1, w2).coeff_monomial(w1**4)
    w24 = sp.Poly(c2, w1, w2).coeff_monomial(w2**4)
    ii_ok = zero(w14 + (A0 * mu1) ** 4) and zero(w24 + (B0 * mu2) ** 4)
    # (ii) double root: r = (w, rho w), w = e^(mu rho); R2' = (1 + mu rho) w, S2 = R2'' + iota R2' = (3 mu + mu^2 rho + iota) w with iota rho = 1
    r2p = (mu * rho if MUT == "s-rhoconst-double" else 1 + mu * rho) * w
    rd = {a0: w, a1: mu * w, a2: (mu**2 + io * mu) * w, b0: rho * w, b1: r2p, b2: (3 * mu + mu**2 * rho + io) * w}
    Ed = sp.expand(E_SYM.subs(rd, simultaneous=True).subs(io, 1 / rho))
    low = sp.Poly(sp.expand(Ed * rho**2), rho).coeff_monomial(1)
    iii_ok = zero(low + w**4 * (mu * A0 + B0) ** 4)
    lowest = min(m[0] for m in sp.Poly(sp.expand(Ed * rho**2), rho).monoms())
    check("rho_const", i_ok and ii_ok and iii_ok and lowest == 0,
          f"simple roots: E is a polynomial of degree 2 in iota = 1/rho with iota^2 coefficient -(A0 mu1 w1 + B0 mu2 w2)^4 = -u_rho^4: {i_ok}, so the "
          f"w1^4 iota^2 and w2^4 iota^2 coefficients are -(A0 mu1)^4 and -(B0 mu2)^4: {ii_ok}; double root: the rho^-2 coefficient of E(rho) is "
          f"-w^4 (mu A0 + B0)^4 and no lower power occurs: {iii_ok and lowest == 0}")


# ---------------------------------------------------------------------------------------------------------------
# planted additive family through the same E
# ---------------------------------------------------------------------------------------------------------------
def add_parts():
    """E restricted to X2 = 1, R1 = 1, v = 0 (u = X + R), its bilinear form and a rank factorization from pivot columns."""
    Esp = sp.expand(E_SYM.subs({B0: 1, B1: 0, B2: 0, a0: 1, a1: 0, a2: 0}).subs(v, 0))
    M, rows, cols = bil(Esp, (A0, A1, A2), (b0, b1, b2))
    Rm, piv = M.rref()
    r = len(piv)
    xmon = [mon_value(m, (A0, A1, A2)) for m in rows]
    rmon = [mon_value(m, (b0, b1, b2)) for m in cols]
    f = [sp.expand(sum(M[i, piv[k]] * xmon[i] for i in range(M.rows))) for k in range(r)]
    g = [sp.expand(sum(Rm[k, j] * rmon[j] for j in range(M.cols))) for k in range(r)]
    return Esp, f, g, r


def add_check_identity(Esp, f, g):
    return zero(Esp - sum(a * b for a, b in zip(f, g)))


BESSEL_CASES = [(sp.Rational(n, d), e) for (n, d) in ((3, 10), (1, 2), (1, 1), (3, 2), (2, 1)) for e in (1, -1)]
ZPTS = [mp.mpf(2) / 5 * (1 + j) for j in range(5)]
ZRANK = [mp.mpf(2) / 5 + mp.mpf(8) / 5 * j / 11 for j in range(12)]
rho_s = sp.Symbol("rho_s", positive=True)


def bessel_R(c, eps, pert=0, wrong_freq=False):
    z = (1 if wrong_freq else sp.sqrt(2)) * c * rho_s
    return -eps * sp.log(sp.besselj(0, z) + sp.Rational(3, 10) * sp.bessely(0, z)) + pert


def g_numeric(g, R, rhos):
    """Values of the g_k (polynomials in b0, b1, b2) along u = ... + R(rho) at the points rhos."""
    R1 = sp.diff(R, rho_s)
    S = sp.diff(R, rho_s, 2) + R1 / rho_s
    gl = [sp.lambdify(rho_s, gk.subs({b0: R, b1: R1, b2: S}), "mpmath") for gk in g]
    return [[gk(r) for r in rhos] for gk in gl]


def rho_from_z(c, zs):
    return [z / (mp.sqrt(2) * mp.mpf(c.p) / c.q) for z in zs]


def num_rank(rows, tol=mp.mpf(10) ** -40):
    s = mp.svd_r(mp.matrix(rows), compute_uv=False)
    smax = max(s)
    return sum(1 for t in s if t > tol * smax)


def c_add_decomp():
    Esp, f, g, r = add_parts()
    ident = add_check_identity(Esp, f, g)
    check("add_decomp", ident and r == 5,
          f"E on X2 = 1, R1 = 1, v = 0 (u = X + R) is 4 R2'^2 E_add; E = sum_(k<=r) f_k(x) g_k(rho) from pivot columns: {ident}; exact rank r = {r} (needs 5); "
          f"f = {[str(t) for t in f]}")


def c_add_stratum():
    Esp, f, g, r = add_parts()
    ks = sp.symbols("k1:%d" % (r + 1))
    use = list(range(r)) if MUT != "s-stratum-drop" else [0, 1, 2]
    eqs = [f[i] - ks[i] for i in use]
    gens = (A1, A2, A0) + tuple(ks)
    G = sp.groebner(eqs, *gens, order="lex")
    elim = [e for e in G.exprs if not (e.has(A1) or e.has(A2))]
    cond = []
    for e in elim:
        cond += sp.Poly(e, A0).all_coeffs() if e.has(A0) else [e]
    cond = [cc for cc in cond if cc != 0]
    G2 = sp.groebner(eqs + cond, *gens, order="lex")
    x2_ok = G2.contains(A2**2)
    c = sp.Symbol("c")
    free_ok = all(sp.simplify(sp.diff(fk.subs({A1: c, A2: 0}), A0)) == 0 for fk in f)
    check("add_stratum", x2_ok and free_ok,
          f"x-curves with all f_k(x) constant (d_F = 1): elimination of {{f_k - k_k}} with A0 free gives conditions {cond} on k; X'' = A2: A2^2 in the ideal: "
          f"{x2_ok}; for every constant c the jets (A0 free, A1 = c, A2 = 0) have f independent of A0: {free_ok}")


def c_add_planted_exact():
    c, eps, rho = sp.symbols("c eps rho", nonzero=True)
    psi = sp.Function("psi")(rho)
    ok = True
    for e in (1, -1):
        R = -e * sp.log(psi)
        R1 = sp.diff(R, rho)
        S = sp.diff(R, rho, 2) + R1 / rho
        # Bessel: psi'' + psi'/rho + 2 c^2 psi = 0
        bess = (-sp.Derivative(psi, rho) / rho - 2 * c**2 * psi) if MUT != "s-add-exact" else (-sp.Derivative(psi, rho) / rho + 2 * c**2 * psi)
        S = S.subs(sp.Derivative(psi, (rho, 2)), bess)
        R1 = R1
        val = dc.E_add_poly(A0, c, 0, 0, R1, S)
        ok &= zero(sp.simplify(val))
    check("add_planted_exact", ok,
          f"R = -eps ln psi with psi'' + psi'/rho + 2 c^2 psi = 0 gives E_add(X = c x) = (R'^2/2 + c^2)^2 - S^2/4 = 0 identically, eps = +-1: {ok}")


def c_add_planted():
    Esp, f, g, r = add_parts()
    cs = sp.Symbol("c")
    worst, dgs, dFs, amp_ok = mp.mpf(0), [], [], True
    for (c, eps) in BESSEL_CASES:
        cm = mp.mpf(c.p) / c.q
        fv = [fk.subs({A0: c * x0, A1: c, A2: 0}) for fk in f for x0 in (1,)]
        f2 = [fk.subs({A0: 2 * c, A1: c, A2: 0}) for fk in f]
        dFs.append(0 if all(a == b for a, b in zip(fv, f2)) else 1)
        rhos = rho_from_z(c, ZPTS)
        R = bessel_R(c, eps, 0, MUT == "s-add-freq")
        gm = g_numeric(g, R, rhos)
        fvn = [mp.mpf(sp.Rational(t).p) / sp.Rational(t).q for t in fv]
        for j in range(len(rhos)):
            terms = [fvn[k] * gm[k][j] for k in range(r)]
            worst = max(worst, abs(sum(terms)) / sum(abs(t) for t in terms))
        dgs.append(num_rank(g_numeric(g, R, rho_from_z(c, ZRANK))))
        amp_ok &= c >= sp.Rational(3, 10)
    ok = worst <= mp.mpf(10) ** -50 and max(dgs) <= r - 1 and max(dFs) == 0 and amp_ok
    check("add_planted", ok,
          f"u = c x - eps ln(J0(z) + 0.3 Y0(z)) + d (z = sqrt2 c rho), c in {{3/10, 1/2, 1, 3/2, 2}}, eps = +-1: f(x) constant in x (d_F = 1): {max(dFs) == 0}; "
          f"largest relative residual of f.g over 50 points {mp.nstr(worst, 3)} (needs <= 1e-50); numerical d_G {dgs} (needs <= {r - 1}); "
          f"d_F + d_G = {1 + max(dgs)} = r = {r}: saturated; c >= 3/10: {amp_ok}")


def c_add_perturbed_exact():
    amp = 0 if MUT == "s-pexact-zero" else 1
    gam, r, c = sp.symbols("gamma r c")
    poly_ok = True
    for e in (1, -1):
        for dl in (r**2 / 10, r**3 / 10):
            dl = amp * dl
            d1 = sp.diff(dl, r)
            d2 = sp.diff(dl, r, 2) + d1 / r
            Ep = sp.expand(sp.cancel(dc.E_add_poly(A0, c, 0, 0, gam + d1, e * (gam**2 + 2 * c**2) + d2)))
            poly_ok &= Ep != 0
    pts = [(sp.Rational(1, 2), 1, 1, sp.Rational(1, 3)), (1, -1, sp.Rational(3, 2), 2), (2, 1, sp.Rational(1, 2), -1)]
    vals, pt_ok = [], True
    for (cv, e, rv, gv) in pts:
        for dl in (r**2 / 10, r**3 / 10):
            dl = amp * dl
            d1 = sp.diff(dl, r)
            d2 = sp.diff(dl, r, 2) + d1 / r
            val = sp.nsimplify(dc.E_add_poly(A0, cv, 0, 0, gv + d1.subs(r, rv), e * (gv**2 + 2 * cv**2) + d2.subs(r, rv)))
            vals.append(val)
            pt_ok &= val != 0
    check("add_perturbed_exact", poly_ok and pt_ok,
          f"E_add(X = c x, R + delta) with the Bessel equation for R, delta in {{rho^2/10, rho^3/10}}, eps = +-1: nonzero polynomial in (gamma, rho, c): {poly_ok}; "
          f"nonzero at the exact points (c, eps, rho, gamma) = (1/2, 1, 1, 1/3), (1, -1, 3/2, 2), (2, 1, 1/2, -1), values {[str(t) for t in vals]}: {pt_ok}")


def c_add_xperturbed():
    Esp, f, g, r = add_parts()
    amp = 0 if MUT == "s-xperturb-zero" else 1
    dFs, dGs = [], []
    for (c, eps) in BESSEL_CASES:
        rows = []
        for j in range(10):
            x = sp.Rational(j, 4)
            jets = {A0: c * x + amp * x**2 / 10, A1: c + amp * x / 5, A2: amp * sp.Rational(1, 5)}
            rows.append([fk.subs(jets) for fk in f])
        dFs.append(sp.Matrix(rows).rank())
        dGs.append(num_rank(g_numeric(g, bessel_R(c, eps), rho_from_z(c, ZRANK))))
    ok = all(a >= 2 for a in dFs) and all(b == r - 1 for b in dGs)
    check("add_xperturbed", ok,
          f"X = c x + x^2/10 with the same R: exact d_F per case {dFs} (needs >= 2), numerical d_G {dGs} (needs {r - 1}), so d_F + d_G > r = {r} and Lemma L rejects it: {ok}")


def c_irreducible():
    fl = sp.factor_list(E_SYM)
    check("irreducible", len(fl[1]) == 1 and fl[1][0][1] == 1,
          f"factor_list(E): constant {fl[0]}, {len(fl[1])} factor(s) of {[len(sp.Poly(f).terms()) for f, m in fl[1]]} terms, multiplicities {[m for f, m in fl[1]]}")


CHECKS = [("e_two", c_e_two), ("ode_form", c_ode_form), ("case_a", c_case_a), ("sep_rank", c_sep_rank), ("slack", c_slack),
          ("rho_const", c_rho_const), ("add_decomp", c_add_decomp), ("add_stratum", c_add_stratum),
          ("add_planted_exact", c_add_planted_exact), ("add_planted", c_add_planted), ("add_perturbed_exact", c_add_perturbed_exact),
          ("add_xperturbed", c_add_xperturbed)]
if "irreducible" in OPTS:
    CHECKS.append(("irreducible", c_irreducible))


def wanted(name):
    return ONLY is None or name in ONLY


def verdict():
    ok = {n: s_ for n, s_ in RESULTS}
    red = all(ok.get(n, False) for n in ("e_two", "ode_form", "case_a"))
    rk = ok.get("sep_rank", False) and ok.get("slack", False)
    plant = all(ok.get(n, False) for n in ("add_decomp", "add_stratum", "add_planted_exact", "add_planted", "add_perturbed_exact", "add_xperturbed"))
    print("\nverdicts:")
    print("  reductions (ODE normal form, excluded forms are one product plus a constant or additive): "
          + ("established" if red else "NOT established (a check failed)"))
    print("  rho-side with constant ODE coefficients (exponentials, exponential polynomials): " + ("no solution under W_x W_rho != 0"
          if ok.get("rho_const", False) else "NOT established"))
    print("  Lemma L on the two-product E: " + ("ranks exact; the bound does not decide (exponential pair has room)" if rk else "NOT established"))
    print("  planted additive truth through the same E: " + ("found, perturbed members rejected exactly" if plant else "NOT established"))
    print("  the question (no solution under (N), or an explicit solution): NOT DECIDED by this script (strata (b), (c), (d) outside rho-side constant coefficients)")


def main():
    t0 = time.time()
    if MUT is not None and MUT not in MUTATIONS and MUT not in DECLARED_SILENT:
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
        tg = MUTATIONS[MUT][0] if MUT in MUTATIONS else DECLARED_SILENT[MUT][0]
        print(f"control {MUT}: target check {tg}, fired {tg in failed}")
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
            extra = [f for f in fails if f != target and f not in failed]
            fired = target in fails and not extra and p.returncode == 1
            controls_ok &= fired
            print(f"[{'FIRED' if fired else 'SILENT' if target not in fails else 'EXTRA'}] --mutate={name} ({what}): exit {p.returncode}, failed {fails}",
                  flush=True)
            if p.returncode not in (0, 1):
                print(p.stderr[-2000:])
        for name, (target, what) in DECLARED_SILENT.items():
            p = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True)
            if "logdir" in OPTS:
                with open(os.path.join(OPTS["logdir"], f"ctl_{name}.log"), "w") as fh:
                    fh.write(p.stdout + p.stderr)
            fails = [ln.split(":")[0][7:] for ln in p.stdout.splitlines() if ln.startswith("[FAIL]")]
            print(f"[{'SILENT (declared in README A1, not counted)' if target not in fails else 'FIRED (unexpected)'}] --mutate={name} ({what}): exit {p.returncode}, "
                  f"failed {fails}", flush=True)
        print(f"controls: {'all fired on their target only' if controls_ok else 'SOME DID NOT FIRE ON TARGET ONLY'}")
    sys.exit(0 if not failed and controls_ok else 1)

#!/usr/bin/env python3
"""Time-dependent additive shifts b = X(t, x) + R(t, rho) at the one-sign stage (issue #166)

With Z = b_x_rho / 2 = 0 on a box, b = X(t, x) + R(t, rho). By the Dichotomy of #158 either R_rho = 0 (planar, decided in
#158) or h = R_rho^2 / 2 > 0. Then Lemma R (#158) with Z = 0 forces v0 v2 = 0, two branches: A (v0 = 0) and B (v2 = 0).

    1. the block of X + R from the coordinate Ricci tensor; 2. the gauge x -> x + D(t); 3. Lemma Z (the two branches);
    4. branch A, Q = 0: R = C(t) ln rho + d(t); 5. branch A, the rest: X = c(t) x + e(t), c' = 2 g c - g^2, g = C'/(2C);
    6. the solution c = |C| (k0 - K); 7. branch B, Y = 0: R_rho = m(t) f(rho); 8. branch B: no solution with v0 != 0;
    9. the planted static members, their gauge copies and the new time-dependent members through one function
       (`run_family`: 4D Ricci, the algebraic stage at exact points, the branch residuals);
   10. the general identity R_ab = lam (g_ab - rho^2 dphi_a dphi_b) - alpha_a alpha_b of the whole family;
   11. closure, the general member: [n, K2] and alpha ^ d alpha ^ dphi; 12. closure at exact members (Lemma K, `closure2`);
   13. b_xx = 0: nothing is localized in x.

    python3 notes/warp/one_sign_time.py                 all checks, then every control (exit 1 on any failure)
    python3 notes/warp/one_sign_time.py --mutate=NAME   one control (its target checks must fail, exit 1)
    python3 notes/warp/one_sign_time.py --no-controls   the checks only
    python3 notes/warp/one_sign_time.py --stop=N        sections up to N

Imported read-only. From one_sign_classify.py (#158): `blocks` (the block from jets), `block3`, `v_block` (Lemma R's
block). From same_sign_scalars.py (#78): `same_sign_form`, `rank1_lams`, `closure2`, `simp_r`. From phantom_source.py:
`ricci`. From two_scalars.py: `rq`. The branch equations, Lemma Z, the solve of each branch, the families and the
controls are new here. The README of negative-mass-raw/I166-rpi5 holds the criteria, which this code tests as worded.
"""
import random
import subprocess
import sys
import time

import sympy as sp

import one_sign_classify as osc  # noqa: E402
import phantom_source as ps  # noqa: E402
import same_sign_scalars as sss  # noqa: E402
import two_scalars as ts  # noqa: E402

# name: (checks that must fail, what the mutation does)
MUTATIONS = {
    "block-P": (["block"], "builds P without the -X_x^2 term"),
    "gauge-drop": (["gauge"], "the gauge-transformed X without - D'(t)"),
    "gauge-ode": (["gauge"], "solves D' - c D = d in place of D' + c D = d"),
    "lemmaZ-A": (["lemmaZ"], "claims the elimination ideal <Q, 2hP - Y^2>"),
    "lemmaZ-B": (["lemmaZ"], "claims the elimination ideal <Y, Q^2 - P(2h + P)>"),
    "lemmaZ-blocks": (["lemmaZ"], "builds the branch-B blocks with P = +v0^2"),
    "rho-dsolve": (["rho_ode"], "solves R'' + 2R'/rho = 0 and compares with C1 + C2 ln rho"),
    "A-P": (["branchA"], "builds P with X in place of X + R"),
    "A-ode": (["branchA"], "claims C^2 c' - c C C' - C'^2/4"),
    "c-ode": (["c_ode"], "takes K' = -g^2/A"),
    "t-dsolve": (["t_ode"], "claims G = C1 exp(int kappa) in place of C1 exp(2 int kappa)"),
    "B-claim": (["branchB"], "claims m1 - m2 in the ideal in place of m1^2 - m2^2"),
    "B-const": (["branchB"], "takes g = f0^2/4 (the 1/rho^2 dropped) in the f-constant case"),
    "plant-pert": (["planted"], "puts F1 = P1 + x^2/10 in the list of members that must be found"),
    "fam-ode": (["family"], "puts F2 in place of T2 in the list of members that must be found"),
    "ident-ode": (["identity"], "substitutes c' = 2 g c + g^2"),
    "ident-alpha": (["identity"], "takes alpha with (c - g) d rho"),
    "closure-bracket": (["closure_general"], "computes [n, K2] without the term K2(n^mu)"),
    "closure-form": (["closure_general"], "takes alpha without the b dt term"),
    "closure-C0": (["closure_points"], "replaces two of the seven members by the integrable b = x"),
    "loc-claim": (["localized"], "claims the ODE as c' = 2 g c - g^2 c"),
}
SECTION = {"block": 1, "gauge": 2, "lemmaZ": 3, "rho_ode": 4, "branchA": 5, "c_ode": 6, "t_ode": 7, "branchB": 8,
           "planted": 9, "family": 9, "identity": 10, "closure_general": 11, "closure_points": 12, "localized": 13}
STOP_AT = {k: SECTION[v[0][0]] for k, v in MUTATIONS.items()}

RESULTS = []
MUT = None


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def info(text):
    print(f"[info] {text}", flush=True)


def fin(t0):
    bad = [n for n, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS)} checks, {len(RESULTS) - len(bad)} pass" + (f", FAILED: {bad}" if bad else "") +
          f" ({time.time() - t0:.0f} s)")
    return bool(bad)


rq = ts.rq
simp_r = sss.simp_r
rho = sp.Symbol("rho", positive=True)


def zero(e):
    if isinstance(e, sp.MatrixBase):
        return all(sp.simplify(x) == 0 for x in e)
    return sp.simplify(e) == 0


def to_jets(expr, F):
    """Replaces every derivative of the undefined function F(t, .), plain or inside a Subs, and F itself, by the symbol
    F_ab (a: order in t, b: order in the second slot). Subs(symbol, ...) collapses under doit()."""
    def jet(e):
        f = e.expr if isinstance(e, sp.Derivative) else e
        a = b = 0
        if isinstance(e, sp.Derivative):
            for v, n in e.variable_count:
                if v == f.args[0]:
                    a += n
                else:
                    b += n
        return sp.Symbol(f"{F.__name__}_{a}{b}")

    expr = expr.replace(lambda e: isinstance(e, sp.Derivative) and e.expr.func == F, jet)  # derivatives first
    return expr.replace(lambda e: e.func == F, jet).doit()


def entries_X_R(Xe, Re, t, x):
    """(P, Q, Y, h) of the Reduction for b = Xe(t, x) + Re(t, rho); Z = 0."""
    b = Xe + Re
    P = Xe.diff(t, x) - b * Xe.diff(x, 2) - (0 if MUT == "block-P" else Xe.diff(x)**2)
    S = Re.diff(rho, 2) + Re.diff(rho) / rho
    Y = -Re.diff(t, rho) / 2 + Xe.diff(x) * Re.diff(rho)
    return [P, -S / 2, Y, Re.diff(rho)**2 / 2]


# 4D shift metrics, chain-rule symbols ------------------------------------------------------------------------------
tt, xx, ph = sp.symbols("t x varphi", real=True)
Lg = sp.Symbol("Lg", real=True)  # ln(rho), d Lg/d rho = 1/rho
En, Th, Ti = sp.symbols("En Th Ti", real=True)  # e^(2t), tanh t, 1/t
X4 = [tt, xx, rho, ph]


def Dfun(E, k):
    """d/dx^k of an expression in (t, x, rho, Lg, En, Th, Ti) with En' = 2 En, Th' = 1 - Th^2, Ti' = -Ti^2."""
    if k == 0:
        return sp.diff(E, tt) + sp.diff(E, En) * 2 * En + sp.diff(E, Th) * (1 - Th**2) - sp.diff(E, Ti) * Ti**2
    if k == 2:
        return sp.diff(E, rho) + sp.diff(E, Lg) / rho
    return sp.diff(E, X4[k])


def shift_g(b):
    return sp.Matrix([[-1 + b**2, b, 0, 0], [b, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, rho**2]])


FAMILIES = {
    "P1": sp.Rational(3, 2) * xx + 2 * Lg + sp.Rational(1, 5),
    "P2": xx + 3 * Lg,
    "P3": sp.Rational(3, 2) * xx + 2 * Lg + Th,
    "P4": sp.Rational(3, 2) * (xx - Th) + 2 * Lg + sp.Rational(1, 5) - (1 - Th**2),
    "T1": xx / 2 + En * Lg,
    "T2": (sp.Rational(1, 2) + En / 3) * xx + En * Lg,
    "T3": -Th * xx + En * (1 - Th**2) * Lg,
    "T4": -Ti * xx + Ti**2 * Lg,
    "F1": sp.Rational(3, 2) * xx + 2 * Lg + sp.Rational(1, 5) + xx**2 / 10,
    "F2": (sp.Rational(1, 2) + En / 3 + sp.Rational(1, 5)) * xx + En * Lg,
    "A0": xx,
}
CACHE = {}


def lam_of(b):
    """C = rho b_rho, gamma = C'/(2C), c = b_x and lam = b_rho^2/2 + (gamma - c)^2, by the chain rule (C != 0)."""
    br = Dfun(b, 2)
    Cc = sp.simplify(rho * br)
    cc = sp.simplify(Dfun(b, 1))
    if Cc == 0:
        return Cc, None, cc, br**2 / 2 + cc**2  # planar: lam = K = c^2 (b = x)
    gam = sp.simplify(Dfun(Cc, 0) / (2 * Cc))
    return Cc, gam, cc, sp.simplify(br**2 / 2 + (gam - cc)**2)


def draw_points(seed, Cc, n=3):
    random.seed(seed)
    pts = []
    while len(pts) < n:
        p = {rho: sp.Rational(random.randint(5, 30), 10), Lg: rq(), tt: rq(), xx: rq(),
             En: sp.Rational(random.randint(1, 9), random.randint(1, 5)), Th: sp.Rational(random.randint(-8, 8), 9),
             Ti: sp.Rational(random.randint(1, 9), random.randint(1, 5))}
        if sp.Rational(sp.cancel(sp.sympify(Cc).subs(p))) != 0 or Cc == 0:
            pts.append(p)
    return pts


def at(p, e):
    return sp.Rational(sp.cancel(sp.sympify(e).subs(p)))


def run_family(name, b, seed):
    """The code every planted and every new member goes through. Returns the 4D Ricci tensor, the metric, lam, the
    branch residuals (Q, 2hP + Y^2) from the jets of b by #158's `blocks`, and the algebraic stage at 3 exact points:
    lam_ok (lam is a root with inertia (0, 2) and pattern MM at every point), none_pat (no root has a pattern)."""
    g = shift_g(b)
    Ric = ps.ricci(g, X4, Dfun, simp_r)
    Cc, gam, cc, lam = lam_of(b)
    bx, br, bt = Dfun(b, 1), Dfun(b, 2), Dfun(b, 0)
    P, Q, Z, Y, h = osc.blocks(b, bx, br, bt, Dfun(bt, 1), Dfun(bx, 1), Dfun(bx, 2), Dfun(bt, 2), Dfun(br, 2), rho)
    resQ, resA = sp.simplify(Q), sp.simplify(2 * h * P + Y**2)
    pts = draw_points(seed, Cc)
    lam_ok, none_pat, det = True, True, []
    for p in pts:
        cands = sss.same_sign_form(Ric.applyfunc(lambda e: at(p, e)), g.applyfunc(lambda e: at(p, e)))
        lv = at(p, lam)
        lam_ok &= any(c["lam"] == lv and c["sig"] == (0, 2) and "MM" in c["pats"] for c in cands)
        none_pat &= not any(c["pats"] for c in cands)
        det.append((str(lv), [(float(sp.N(c["lam"], 5)), c["sig"], c["pats"]) for c in cands]))
    out = dict(name=name, b=b, g=g, Ric=Ric, Cc=Cc, gam=gam, cc=cc, lam=lam, resQ=resQ, resA=resA, lam_ok=lam_ok,
               none_pat=none_pat, det=det, seed=seed)
    CACHE[name] = out
    return out


def main():
    global MUT
    opts = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
    mut = MUT = opts.get("mutate")
    if mut is not None and mut not in MUTATIONS:
        sys.exit(f"unknown mutation {mut!r}. Known: {', '.join(MUTATIONS)}")
    t0 = time.time()
    stop = int(opts.get("stop", STOP_AT.get(mut, 99)))
    print(f"one_sign_time.py, issue #166. sympy {sp.__version__}, python {sys.version.split()[0]}. mutation: "
          f"{mut or 'none'}" + (f" ({MUTATIONS[mut][1]})" if mut else ""))
    t, x, ph_ = sp.symbols("t x phi", positive=True)
    v0, v1, v2, P_, Q_, Y_, h_ = sp.symbols("v0 v1 v2 P Q Y h")

    # 1. the block of X(t, x) + R(t, rho) ------------------------------------------------------------------------------
    print("\n== 1. The block R_3 of b = X(t, x) + R(t, rho) from the coordinate Ricci tensor ==")
    Xf, Rf = sp.Function("X")(t, x), sp.Function("R")(t, rho)
    b = Xf + Rf
    g = sp.Matrix([[-1 + b**2, b, 0, 0], [b, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, rho**2]])
    Xg = [t, x, rho, ph_]
    Rc = ps.ricci(g, Xg, lambda E, k: sp.diff(E, Xg[k]), sp.simplify)
    E = sp.Matrix([[1, 0, 0, 0], [-b, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1 / rho]])
    Rfr = (E.T * Rc * E).applyfunc(sp.simplify)
    P, Q, Y, h = entries_X_R(Xf, Rf, t, x)
    claimed = sp.Matrix([[P - h, Q, 0, 0], [Q, -P - h, Y, 0], [0, Y, h, 0], [0, 0, 0, 0]])
    ok_blk = zero(Rfr - claimed)
    ok_i2 = zero(Rfr[0, 0] + Rfr[1, 1] + 2 * Rfr[2, 2])
    bj = osc.blocks(b, Xf.diff(x), Rf.diff(rho), Xf.diff(t) + Rf.diff(t), Xf.diff(t, x), Xf.diff(x, 2), 0,
                    Rf.diff(t, rho), Rf.diff(rho, 2), rho)
    ok_158 = zero(bj[0] - (Xf.diff(t, x) - b * Xf.diff(x, 2) - Xf.diff(x)**2)) and zero(bj[1] - Q) and bj[2] == 0 and \
        zero(bj[3] - Y) and zero(bj[4] - h)
    rt = sp.Function("r")(t)
    Pp = Xf.diff(t, x) - (Xf + rt) * Xf.diff(x, 2) - Xf.diff(x)**2
    ok_pl = zero(osc.block3(*osc.blocks(Xf + rt, Xf.diff(x), 0, Xf.diff(t) + rt.diff(t), Xf.diff(t, x), Xf.diff(x, 2), 0,
                                        0, 0, rho)) - sp.diag(Pp, -Pp, 0))
    check("block", ok_blk and ok_i2 and ok_158 and ok_pl,
          f"coordinate Ricci of b = X(t, x) + R(t, rho) on (n, d_x, d_rho, d_phi/rho) minus [[P-h, Q, 0], [Q, -P-h, Y], [0, Y, "
          f"h]] + 0 with P = X_tx - (X + R) X_xx - X_x^2, Q = -(R_rr + R_r/rho)/2, Y = -R_tr/2 + X_x R_r, h = R_r^2/2 is zero: "
          f"{ok_blk}; (I2): {ok_i2}; #158's `blocks` on the jets gives (P, Q, 0, Y, h): {ok_158}; R = r(t) gives "
          f"diag(P, -P, 0): {ok_pl} ({time.time() - t0:.0f} s)")
    if stop < 2:
        return fin(t0)

    # 2. the gauge ----------------------------------------------------------------------------------------------------
    print("\n== 2. The gauge x -> x + D(t) ==")
    Xc = sp.Function("X")
    D = sp.Function("D")(t)
    sg = 0 if MUT == "gauge-drop" else 1
    Xt = Xc(t, x - D) - sg * D.diff(t)
    PQYh_new = [to_jets(e, Xc) for e in entries_X_R(Xt, Rf, t, x)]
    PQYh_old = [to_jets(e, Xc) for e in entries_X_R(Xc(t, x), Rf, t, x)]
    assert PQYh_old[0].has(sp.Symbol('X_11')) and PQYh_new[0].has(sp.Symbol('X_11')), 'vacuous gauge comparison'
    ok_ga = all(sp.simplify(a_ - b_) == 0 for a_, b_ in zip(PQYh_new, PQYh_old))
    cf, Cf, df = sp.Function("c")(t), sp.Function("C")(t), sp.Function("d")(t)
    bb = cf * x + Cf * sp.log(rho) + df
    bt_ = bb.subs(x, x - D) - D.diff(t)
    ok_gb1 = sp.simplify(sp.expand(bt_ - (cf * x + Cf * sp.log(rho) + (df - cf * D - D.diff(t))))) == 0
    Ic = sp.Integral(cf, t)
    Dsol = (sp.exp(Ic) * sp.Integral(df * sp.exp(-Ic), t)) if MUT == "gauge-ode" else \
        (sp.exp(-Ic) * sp.Integral(df * sp.exp(Ic), t))
    ok_gb2 = sp.simplify(Dsol.diff(t) + cf * Dsol - df) == 0
    check("gauge", ok_ga and ok_gb1 and ok_gb2,
          f"(a) (P, Q, Y, h) of X~ = X(t, x - D) - D', R~ = R at (t, x, rho) equal those of b at (t, x - D, rho) (undefined X, R, "
          f"D, derivatives of X replaced by jet symbols): {ok_ga}; (b) c x + C ln rho + d gives b~ = c x + C ln rho + (d - c D "
          f"- D'): {ok_gb1}, and D = e^(-int c) int d e^(int c) satisfies D' + c D - d = 0: {ok_gb2}")
    if stop < 3:
        return fin(t0)

    # 3. Lemma Z ------------------------------------------------------------------------------------------------------
    print("\n== 3. Lemma Z: the two branches of Lemma R at Z = 0 ==")
    GA = sp.groebner([P_ + v2**2, Q_, Y_ + v1 * v2, 2 * h_ - v1**2], v1, v2, P_, Q_, Y_, h_, order="lex")
    elA = [e for e in GA.exprs if not e.has(v1) and not e.has(v2)]
    claimA = [Q_, 2 * h_ * P_ + (-1 if MUT == "lemmaZ-A" else 1) * Y_**2]
    GAc = sp.groebner(claimA, P_, Q_, Y_, h_, order="lex")
    ok_za = all(GA.reduce(e)[1] == 0 for e in GAc.exprs) and all(GAc.reduce(e)[1] == 0 for e in elA)
    GB = sp.groebner([P_ + v0**2, Q_ + v0 * v1, Y_, 2 * h_ - v0**2 - v1**2], v0, v1, P_, Q_, Y_, h_, order="lex")
    elB = [e for e in GB.exprs if not e.has(v0) and not e.has(v1)]
    claimB = [Y_, Q_**2 + (-1 if MUT == "lemmaZ-B" else 1) * P_ * (2 * h_ + P_)]
    GBc = sp.groebner(claimB, P_, Q_, Y_, h_, order="lex")
    ok_zb = all(GB.reduce(e)[1] == 0 for e in GBc.exprs) and all(GBc.reduce(e)[1] == 0 for e in elB)
    random.seed(1661)
    rnz = lambda: sp.Rational(random.choice([-3, -2, -1, 1, 2, 3]), random.randint(1, 3))
    badZ = []
    for _ in range(12):
        a1, a2 = rnz(), rnz()
        Pb = osc.v_block(0, a1, a2)
        if sss.rank1_lams(osc.block3(*Pb)) != [(Pb[4] + a2**2, "MM")]:
            badZ.append(("A", a1, a2))
    for _ in range(12):
        a0, a1 = rnz(), rnz()
        Pb = osc.v_block(a0, a1, 0)
        if MUT == "lemmaZ-blocks":
            Pb = (a0**2,) + tuple(Pb[1:])
        if sss.rank1_lams(osc.block3(*Pb)) != [(Pb[4], "MM")]:
            badZ.append(("B", a0, a1))
    ngen = 0
    while ngen < 12:
        Pg, Qg, Yg, hg = rq(), rq(), rq(), abs(rq()) + sp.Rational(1, 7)
        if (Qg == 0 and 2 * hg * Pg + Yg**2 == 0) or (Yg == 0 and Qg**2 + Pg * (2 * hg + Pg) == 0):
            continue
        ngen += 1
        if sss.rank1_lams(osc.block3(Pg, Qg, 0, Yg, hg)) != []:
            badZ.append(("generic", Pg, Qg, Yg, hg))
    check("lemmaZ", ok_za and ok_zb and not badZ,
          f"(a) elimination ideal of Lemma R at v0 = 0 equals <Q, 2hP + Y^2> (basis {[str(e) for e in elA]}, claim "
          f"{[str(e) for e in claimA]}): {ok_za}; (b) at v2 = 0 equals <Y, Q^2 + P(2h + P)> (basis {[str(e) for e in elB]}): "
          f"{ok_zb}; (c) `rank1_lams` on 12 A-blocks gives [(h + v2^2, MM)], on 12 B-blocks [(h, MM)], on 12 generic Z = 0 "
          f"blocks []: wrong in {len(badZ)} {badZ[:2]}")
    if stop < 4:
        return fin(t0)

    # 4. branch A: Q = 0 -----------------------------------------------------------------------------------------------
    print("\n== 4. Branch A: Q = 0 gives R = C ln rho + d ==")
    Rr = sp.Function("R")(rho)
    C1, C2 = sp.symbols("C1 C2")
    if MUT == "rho-dsolve":
        sol = sp.dsolve(sp.Eq(Rr.diff(rho, 2) + 2 * Rr.diff(rho) / rho, 0), Rr)
    else:
        sol = sp.dsolve(sp.Eq(Rr.diff(rho, 2) + Rr.diff(rho) / rho, 0), Rr)
    ok_s1 = sp.simplify(sol.rhs - (C1 + C2 * sp.log(rho))) == 0
    ok_s2 = sp.simplify(sp.diff(sol.rhs, rho, 2) + sp.diff(sol.rhs, rho) / rho) == 0
    check("rho_ode", ok_s1 and ok_s2, f"dsolve(R'' + R'/rho = 0) = {sol.rhs}; equals C1 + C2 ln rho: {ok_s1}; satisfies "
          f"R'' + R'/rho = 0: {ok_s2}")
    if stop < 5:
        return fin(t0)

    # 5. branch A: the rest ---------------------------------------------------------------------------------------------
    print("\n== 5. Branch A: X_xx = 0 and the ODE for c(t) ==")
    Cf, df = sp.Function("C")(t), sp.Function("d")(t)
    ell = sp.Symbol("ell")
    RA = Cf * sp.log(rho) + df
    Xa = sp.Function("X")(t, x)
    PA, QA, YA, hA = entries_X_R(Xa, RA, t, x)
    if MUT == "A-P":
        PA = Xa.diff(t, x) - Xa * Xa.diff(x, 2) - Xa.diff(x)**2
    PA_true = Xa.diff(t, x) - (Xa + RA) * Xa.diff(x, 2) - Xa.diff(x)**2
    EA = rho**2 * (2 * hA * PA + YA**2)
    claim_i = Cf**2 * PA_true + (Xa.diff(x) * Cf - Cf.diff(t) / 2)**2
    ok_i = zero(QA) and zero(EA - claim_i)
    Epoly = sp.expand(sp.expand(EA).subs(sp.log(rho), ell))
    coef1 = sp.simplify(sp.Poly(Epoly, ell).coeff_monomial(ell) + Cf**3 * Xa.diff(x, 2))
    coef0 = sp.simplify(sp.Poly(Epoly, ell).coeff_monomial(1) - (Cf**2 * (Xa.diff(t, x) - (Xa + df) * Xa.diff(x, 2) - Xa.diff(x)**2)
                                                                  + (Xa.diff(x) * Cf - Cf.diff(t) / 2)**2))
    ok_ii = coef1 == 0 and coef0 == 0
    cA, eA = sp.Function("c")(t), sp.Function("e")(t)
    gam = Cf.diff(t) / (2 * Cf)
    Xlin = cA * x + eA
    PA2, QA2, YA2, hA2 = entries_X_R(Xlin, RA, t, x)
    EA2 = sp.simplify(rho**2 * (2 * hA2 * PA2 + YA2**2))
    sgn_ode = -1 if MUT == "A-ode" else 1
    ok_iii = sp.simplify(EA2 - (Cf**2 * cA.diff(t) - cA * Cf * Cf.diff(t) + sgn_ode * Cf.diff(t)**2 / 4)) == 0 and \
        sp.simplify(EA2 - Cf**2 * (cA.diff(t) - 2 * gam * cA + gam**2)) == 0
    ode = {cA.diff(t): 2 * gam * cA - gam**2}
    v1A = Cf / rho
    v2A = sp.simplify(-YA2 / v1A)
    ok_iv = sp.simplify(v2A - (gam - cA)) == 0 and sp.simplify((hA2 - PA2).subs(ode) - (Cf**2 / (2 * rho**2) + (gam - cA)**2)) == 0 \
        and sp.simplify((hA2 + v2A**2).subs(ode) - (hA2 - PA2).subs(ode)) == 0
    check("branchA", ok_i and ok_ii and ok_iii and ok_iv,
          f"R = C ln rho + d: (i) Q = 0 and rho^2 (2hP + Y^2) = C^2 P + (X_x C - C'/2)^2: {ok_i}; (ii) as a polynomial in l = ln rho "
          f"the l-coefficient is -C^3 X_xx and the l^0 part is as stated: {ok_ii}; (iii) for X = c x + e it equals C^2 c' - c C C' "
          f"+ C'^2/4 = C^2 (c' - 2 g c + g^2): {ok_iii}; (iv) v2 = -Y/v1 = g - c and, on the ODE, h - P = h + v2^2 = C^2/(2 rho^2) "
          f"+ (g - c)^2: {ok_iv}")
    if stop < 6:
        return fin(t0)

    # 6. the solution of the ODE ----------------------------------------------------------------------------------------
    print("\n== 6. c' = 2 g c - g^2: c = |C| (k0 - K), K' = g^2/|C| ==")
    A = sp.Function("A")(t)
    gm = sp.Function("gm")(t)
    K = sp.Function("K")(t)
    k0 = sp.Symbol("k0")
    sK = -1 if MUT == "c-ode" else 1
    csol = A * (k0 - K)
    resid = (csol.diff(t) - 2 * gm * csol + gm**2).subs({K.diff(t): sK * gm**2 / A}).subs({A.diff(t): 2 * gm * A})
    ok_c1 = sp.simplify(resid) == 0
    resid0 = csol.diff(t).subs({K.diff(t): sK * gm**2 / A}).subs({A.diff(t): 2 * gm * A}).subs(gm, 0)
    ok_c2 = sp.simplify(resid0) == 0
    tq = sp.Symbol("tq", real=True)
    A2, g2, K2 = sp.exp(2 * tq), sp.Integer(1), -sp.exp(-2 * tq) / 2
    c2 = A2 * (sp.Rational(1, 3) - K2)
    ok_c3 = sp.simplify(sp.diff(K2, tq) - sK * g2**2 / A2) == 0 and sp.simplify(c2 - (sp.Rational(1, 2) + A2 / 3)) == 0
    check("c_ode", ok_c1 and ok_c2 and ok_c3,
          f"with A' = 2 g A, K' = {sK} g^2/A, c = A (k0 - K): c' - 2 g c + g^2 = {sp.simplify(resid)}: {ok_c1}; for g = 0 "
          f"c' = {sp.simplify(resid0)}: {ok_c2}; T2: A = e^(2t), g = 1, K = -e^(-2t)/2 gives c = 1/2 + e^(2t)/3: {ok_c3}")
    if stop < 7:
        return fin(t0)

    # 7. branch B: Y = 0 ------------------------------------------------------------------------------------------------
    print("\n== 7. Branch B: Y = 0 gives R_rho = m(t) f(rho), m' = 2 kappa m ==")
    kap, Gf = sp.Function("kappa")(t), sp.Function("G")(t)
    solG = sp.dsolve(sp.Eq(Gf.diff(t), 2 * kap * Gf), Gf)
    Cq = sp.Symbol("C1")
    claimG = Cq * sp.exp((1 if MUT == "t-dsolve" else 2) * sp.Integral(kap, t))
    ok_t1 = sp.simplify(solG.rhs - claimG) == 0
    ok_t2 = sp.simplify(sp.diff(solG.rhs, t) - 2 * kap * solG.rhs) == 0
    check("t_ode", ok_t1 and ok_t2, f"dsolve(G' = 2 kappa G) = {solG.rhs}; equals {claimG}: {ok_t1}; satisfies the ODE: {ok_t2}")
    if stop < 8:
        return fin(t0)

    # 8. branch B: no solution with v0 != 0 ---------------------------------------------------------------------------------
    print("\n== 8. Branch B: the equation Q^2 + P(2h + P) = 0 has no solution with P < 0 ==")
    kp, mm, ee = sp.Function("kappa")(t), sp.Function("m")(t), sp.Function("e")(t)
    Phi = sp.Function("Phi")(rho)
    fF = Phi.diff(rho)
    RB = mm * Phi
    PB, QB, YB, hB = entries_X_R(kp * x + ee, RB, t, x)
    pp = kp**2 - kp.diff(t)
    gg = (fF.diff(rho) + fF / rho)**2 / 4
    ok_b1 = zero(PB + pp) and zero(QB + mm * (fF.diff(rho) + fF / rho) / 2) and zero(hB - mm**2 * fF**2 / 2) and \
        zero(QB**2 + PB * (2 * hB + PB) - (mm**2 * gg - pp * mm**2 * fF**2 + pp**2))
    p1, p2, m1, m2, Fa, Fb, ga, gb, s = sp.symbols("p1 p2 m1 m2 Fa Fb ga gb s")
    EB = lambda p_, m_, F_, g_: m_**2 * g_ - p_ * m_**2 * F_ + p_**2
    eqs = [EB(p1, m1, Fa, ga), EB(p1, m1, Fb, gb), EB(p2, m2, Fa, ga), EB(p2, m2, Fb, gb)]
    Gs = sp.groebner(eqs + [1 - s * p1 * p2 * m1 * m2 * (Fa - Fb)], s, p1, p2, m1, m2, Fa, Fb, ga, gb, order="grevlex")
    tgt_m = (m1 - m2) if MUT == "B-claim" else (m1**2 - m2**2)
    ok_b2 = Gs.reduce(p1 - p2)[1] == 0 and Gs.reduce(tgt_m)[1] == 0
    pq, mq, f0, ra, rb, ia, ib = sp.symbols("pq mq f0 ra rb ia ib")
    gconst = (lambda i_: f0**2 / 4) if MUT == "B-const" else (lambda i_: f0**2 * i_**2 / 4)
    eqs2 = [EB(pq, mq, f0**2, gconst(ia)), EB(pq, mq, f0**2, gconst(ib)), ra * ia - 1, rb * ib - 1,
            1 - s * pq * mq * f0 * (ra - rb) * (ia + ib)]
    Gu = sp.groebner(eqs2, s, pq, mq, f0, ra, rb, ia, ib, order="grevlex")
    ok_b3 = Gu.exprs == [1]
    check("branchB", ok_b1 and ok_b2 and ok_b3,
          f"(i) with X = kappa x + e, R_rho = m f: P = kappa' - kappa^2, Q = -m (f' + f/rho)/2, h = m^2 f^2/2 and Q^2 + P (2h + "
          f"P) = m^2 g - p m^2 f^2 + p^2 (p = kappa^2 - kappa', g = (f' + f/rho)^2/4): {ok_b1}; (ii) the four equations at two "
          f"times and two radii with 1 - s p1 p2 m1 m2 (Fa - Fb): p1 - p2 reduces to 0 and {tgt_m} reduces to 0: {ok_b2}; (iii) "
          f"f = f0 constant, with 1 - s p m f0 (ra - rb)(ia + ib): Groebner basis {[str(e) for e in Gu.exprs][:2]} is [1]: {ok_b3}. "
          f"By hand: (ii) for every pair of times gives p, m constant, so kappa = 0 and p = kappa^2 - kappa' = 0, contradicting "
          f"p > 0; (iii) f^2 constant; so v0 = 0 and branch B is inside branch A")
    if stop < 9:
        return fin(t0)

    # 9. the families through run_family -----------------------------------------------------------------------------------
    print("\n== 9. Planted and new members through `run_family` ==")
    seeds = {"P1": 1662, "P2": 1663, "P3": 1664, "P4": 1665, "T1": 1666, "T2": 1667, "T3": 1668, "T4": 1669, "F1": 1670,
             "F2": 1671, "A0": 1672}
    for nm in ("P1", "P2", "P3", "P4", "F1"):
        run_family(nm, FAMILIES[nm], seeds[nm])
    must = ["P1", "P2", "P3", "P4"]
    if MUT == "plant-pert":
        must[0] = "F1"
    R_ = CACHE
    ok_found = all(R_[n]["lam_ok"] and R_[n]["resQ"] == 0 and R_[n]["resA"] == 0 for n in must)
    ok_f1 = R_["F1"]["none_pat"] and R_["F1"]["resA"] != 0
    check("planted", ok_found and ok_f1,
          f"P1..P4 (static 3x/2 + 2 ln rho + 1/5, x + 3 ln rho, d(t) = tanh t, and the gauge copy b(t, x - tanh t) - (1 - tanh^2 t)) "
          f"through run_family: lam = C^2/(2 rho^2) + (g - c)^2 is a root with inertia (0, 2) and pattern MM at 3 exact points each "
          f"and Q = 2hP + Y^2 = 0 identically: {ok_found} {[(n, R_[n]['lam_ok'], str(R_[n]['resA'])[:20]) for n in must]}; planted "
          f"failure F1 = P1 + x^2/10: no root with a pattern at 3 of 3 points ({R_['F1']['none_pat']}), 2hP + Y^2 = "
          f"{str(R_['F1']['resA'])[:30]} != 0 ({time.time() - t0:.0f} s)")
    for nm in ("T1", "T2", "T3", "T4", "F2"):
        run_family(nm, FAMILIES[nm], seeds[nm])
    must = ["T1", "T2", "T3", "T4"]
    if MUT == "fam-ode":
        must[1] = "F2"
    ok_found = all(R_[n]["lam_ok"] and R_[n]["resQ"] == 0 and R_[n]["resA"] == 0 for n in must)
    ok_ode = all(sp.simplify(Dfun(R_[n]["cc"], 0) - (2 * R_[n]["gam"] * R_[n]["cc"] - R_[n]["gam"]**2)) == 0 for n in must) \
        and all(sp.simplify(Dfun(R_[n]["Cc"], 0)) != 0 for n in must)
    ok_t4 = sp.simplify(R_["T4"]["gam"] - R_["T4"]["cc"]) == 0 and sp.simplify(R_["T4"]["lam"] - Dfun(FAMILIES["T4"], 2)**2 / 2) == 0
    ok_f2 = R_["F2"]["none_pat"] and R_["F2"]["resA"] != 0
    check("family", ok_found and ok_ode and ok_t4 and ok_f2,
          f"T1..T4 through run_family: lam is a root with inertia (0, 2), MM at 3 exact points each and Q = 2hP + Y^2 = 0 identically: "
          f"{ok_found} {[(n, R_[n]['lam_ok']) for n in must]}; c' = 2 g c - g^2 from the chain rule and C' != 0 (not gauge-equivalent "
          f"to a static member): {ok_ode} (C = {[str(R_[n]['Cc']) for n in must]}, g = {[str(R_[n]['gam']) for n in must]}); T4: "
          f"g = c and lam = h: {ok_t4}; F2 (c ODE violated): no root with a pattern at 3 of 3 points ({R_['F2']['none_pat']}), "
          f"2hP + Y^2 != 0 ({R_['F2']['resA'] != 0}) ({time.time() - t0:.0f} s)")
    if stop < 10:
        return fin(t0)

    # 10. the identity of the whole family -------------------------------------------------------------------------------
    print("\n== 10. R_ab = lam (g_ab - rho^2 dphi_a dphi_b) - alpha_a alpha_b for b = c(t) x + C(t) ln rho + d(t) ==")
    cf, Cf, df = sp.Function("c")(t), sp.Function("C")(t), sp.Function("d")(t)
    bI = cf * x + Cf * sp.log(rho) + df
    gI = sp.Matrix([[-1 + bI**2, bI, 0, 0], [bI, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, rho**2]])
    RicI = ps.ricci(gI, Xg, lambda E_, k: sp.diff(E_, Xg[k]), sp.simplify)
    gmI = Cf.diff(t) / (2 * Cf)
    wI = gmI - cf
    lamI = Cf**2 / (2 * rho**2) + wI**2
    alI = sp.Matrix([Cf / rho * bI, Cf / rho, (-wI if MUT == "ident-alpha" else wI), 0])
    dphi = sp.Matrix([0, 0, 0, rho])
    MI = -alI * alI.T - lamI * dphi * dphi.T
    cprime = 2 * gmI * cf + (1 if MUT == "ident-ode" else -1) * gmI**2
    sub = {cf.diff(t): cprime}
    diffI = (RicI - (lamI * (gI - dphi * dphi.T) - alI * alI.T)).applyfunc(lambda e: sp.simplify(sp.simplify(e).subs(sub)))
    ok_I1 = all(e == 0 for e in diffI)
    nI = sp.Matrix([1, -bI, 0, 0])
    K2I = sp.Matrix([0, wI, -Cf / rho, 0])
    ok_I2 = zero(MI * nI) and zero(MI * K2I)
    check("identity", ok_I1 and ok_I2,
          f"Ricci of the shift metric minus [lam (g - rho^2 dphi^2) - alpha alpha] is the zero matrix after c' -> {cprime}: {ok_I1}; "
          f"with M = -alpha alpha - lam rho^2 dphi dphi: M n = 0 and M K2 = 0 (n = d_t - b d_x, K2 = (g - c) d_x - (C/rho) d_rho): "
          f"{ok_I2}. So R = lam g - alpha^2 - (sqrt(lam) rho dphi)^2: two phantoms for the whole family ({time.time() - t0:.0f} s)")
    if stop < 11:
        return fin(t0)

    # 11. closure, the general member ---------------------------------------------------------------------------------------
    print("\n== 11. Closure of the general member: [n, K2] and alpha ^ d alpha ^ dphi ==")
    Xc4 = [t, x, rho, ph_]
    nv = [1, -bI, 0, 0]
    K2v = [0, wI, -Cf / rho, 0]
    br_ = []
    for m_ in range(4):
        t1 = sum(nv[a_] * sp.diff(K2v[m_], Xc4[a_]) for a_ in range(4))
        t2 = 0 if MUT == "closure-bracket" else sum(K2v[a_] * sp.diff(nv[m_], Xc4[a_]) for a_ in range(4))
        br_.append(sp.simplify(t1 - t2))
    res_ = [sp.simplify(br_[m_] - (Cf.diff(t) / Cf) * K2v[m_]) for m_ in range(4)]
    r_claim = wI.diff(t) + wI * cf - 2 * gmI * wI - Cf**2 / rho**2
    ok_ca = res_[0] == 0 and res_[2] == 0 and res_[3] == 0 and sp.simplify(res_[1] - r_claim) == 0 and \
        sp.simplify(sp.diff(res_[1], rho) - 2 * Cf**2 / rho**3) == 0
    al = [0 if MUT == "closure-form" else Cf / rho * bI, Cf / rho, wI, 0]
    dal = lambda i_, j_: sp.diff(al[j_], Xc4[i_]) - sp.diff(al[i_], Xc4[j_])
    Fform = sp.simplify(al[0] * dal(1, 2) + al[1] * dal(2, 0) + al[2] * dal(0, 1))
    F_claim = Cf**3 / rho**3 + Cf / rho * (-wI.diff(t) + (2 * gmI - cf) * wI)
    ok_cb = sp.simplify(Fform - F_claim) == 0 and sp.simplify(sp.diff(rho * Fform / Cf, rho) + 2 * Cf**2 / rho**3) == 0
    check("closure_general", ok_ca and ok_cb,
          f"(a) [n, K2] - (C'/C) K2 has components (0, r, 0, 0) with r = w' + w c - 2 g w - C^2/rho^2 (w = g - c) and d_rho r = "
          f"2 C^2/rho^3: {ok_ca}; (b) the dt dx drho dphi coefficient of alpha ^ d alpha ^ dphi equals C^3/rho^3 + (C/rho)(-w' + "
          f"(2g - c) w) and d_rho(rho F/C) = -2 C^2/rho^3: {ok_cb}. Both are nonzero for C != 0: the kernel plane is not "
          f"integrable on any open set, so closure fails for every member ({time.time() - t0:.0f} s)")
    if stop < 12:
        return fin(t0)

    # 12. closure at exact members -----------------------------------------------------------------------------------------
    print("\n== 12. Closure (Lemma K) at the planted and the new members ==")
    names = ["P1", "P3", "P4", "T1", "T2", "T3", "T4"]
    if MUT == "closure-C0":
        names[3] = names[4] = "A0"
    if "A0" not in CACHE:
        run_family("A0", FAMILIES["A0"], seeds["A0"])
    out_c, okc = [], True
    for nm in names:
        d_ = CACHE[nm]
        M = (d_["Ric"] - d_["lam"] * d_["g"]).applyfunc(simp_r)
        cl = sss.closure2(M, Dfun, simp_r, d_["lam"])
        okc &= cl["dim"] == 2 and not cl["integrable"]
        out_c.append((nm, cl["dim"], cl["integrable"], [str(e_) for e_ in cl.get("bracket", []) if e_ != 0][:1]))
    d_ = CACHE["A0"]
    cl0 = sss.closure2((d_["Ric"] - d_["lam"] * d_["g"]).applyfunc(simp_r), Dfun, simp_r, d_["lam"])
    ok0 = cl0["dim"] == 2 and cl0["integrable"] and cl0["invariant"] and cl0["flat"]
    check("closure_points", okc and ok0,
          f"closure2 of M = R - lam g (name, kernel dim, integrable, nonzero bracket component): {out_c}: all have dimension 2 and "
          f"are not integrable: {okc}; positive control A0 (b = x, C = 0): dim {cl0['dim']}, integrable {cl0['integrable']}, L_X M = 0 "
          f"{cl0['invariant']}, quotient flat {cl0['flat']}: {ok0} ({time.time() - t0:.0f} s)")
    if stop < 13:
        return fin(t0)

    # 13. localization -----------------------------------------------------------------------------------------------------
    print("\n== 13. Nothing is localized in x ==")
    ok_l1 = sp.simplify(sp.diff(bI, x, 2)) == 0 and sp.simplify(sp.diff(bI, x) - cf) == 0
    gs = sp.Symbol("gs")
    ode_c = 2 * gs * cf - gs**2 * (cf if MUT == "loc-claim" else 1)
    sols = sp.solve(ode_c.subs(cf, 0), gs)
    ok_l2 = sols == [0]
    ok_l3 = sp.simplify((df - cf * D - D.diff(t)).subs(cf, 0) - (df - D.diff(t))) == 0  # C2(b) with c = 0: d~ = d - D'
    check("localized", ok_l1 and ok_l2 and ok_l3,
          f"b_xx = 0 and b_x = c(t) for b = c x + C ln rho + d: {ok_l1}; c = 0 in c' = {ode_c} gives g = {sols}: {ok_l2}; the gauge of (b) with c = 0 gives d~ = d - D' (zero iff D' = d): {ok_l3} (so c = 0 on "
          f"an interval forces C' = 0, b = C ln rho + d(t), and D' = d removes d: static b(rho), x-independent; every other "
          f"member is affine in x with slope c(t))")
    info("Answer: on a box with Z = 0 and R_rho != 0 the stage holds iff (branch A) b = c(t) x + C(t) ln rho + d(t) with C != 0 and "
         "c' = 2 g c - g^2, g = C'/(2C); branch B (v2 = 0, v0 != 0) is empty; d(t) is gauge; closure fails for every member; "
         "none is localized in x")
    return fin(t0)


if __name__ == "__main__":
    failed = main()
    controls_ok = True
    if "--mutate" not in " ".join(sys.argv) and "--no-controls" not in sys.argv and "--stop" not in " ".join(sys.argv):
        print("\n== Controls: each mutation must make its target checks fail ==", flush=True)
        for name, (want, what) in MUTATIONS.items():
            p = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True)
            fails = [ln.split(":")[0][7:] for ln in p.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = all(w in fails for w in want) and p.returncode == 1
            controls_ok &= fired
            print(f"[{'FIRED' if fired else 'SILENT'}] --mutate={name} ({what}): exit {p.returncode}, failed {fails}",
                  flush=True)
            if p.returncode not in (0, 1):
                print(p.stderr[-2000:])
        print(f"controls: {'all fired' if controls_ok else 'SOME DID NOT FIRE'}")
    sys.exit(0 if not failed and controls_ok else 1)

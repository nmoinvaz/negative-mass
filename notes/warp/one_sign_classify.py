#!/usr/bin/env python3
"""Which shifts b(t, x, rho) pass the one-sign scalar stage on an open set? (issue #158)

Two minimally coupled scalars of one sign have R_ab = lam g_ab + s (a_a a_b + b_b b_b)-form (#78). With (I1), (I2) the
algebraic stage holds iff R = 0, or lam != 0 and rank(R_3 - lam eta_3) = 1 on the frame (n, d_x, e_rho), n = d_t - b d_x
(Lemma S of #78). This script classifies the solutions inside the reductions of #158 and runs closure (Lemma K of #78) on
each family that passes:

    1. the block R_3 of a generic b, traveling = static, the gauge x -> x + D(t)
    2. Lemma R: rank-1 <=> N = -v v^T, planted forms through `rank1_lams`
    3. b(rho): b = C ln rho + b_0 (planted), fails closure; Lemma E (the ends of a localized shift)
    4. planar b(t, x): lam = K = -P, closure iff K constant, Theorem P (no localized one)
    5. additive u = X(x) + R(rho) (+ v): u = c x + C ln rho + d; closure
    6. multiplicative u = X(x) R(rho) + v: none
    7. conical u = F(x/rho) + v: none
    8. spherical profile u = w(r) (the bubble at constant speed): no open set (planted, #78)
    9. general static Z != 0: the resolved finite-type system, its compatibility condition (not classified)
   10. the linearization at b = C ln rho + b_0

    python3 notes/warp/one_sign_classify.py                 all checks, then every control (exit 1 on any failure)
    python3 notes/warp/one_sign_classify.py --mutate=NAME   one control (its target checks must fail, exit 1)
    python3 notes/warp/one_sign_classify.py --no-controls   the checks only
    python3 notes/warp/one_sign_classify.py --stop=N        sections up to N

Imported read-only. From same_sign_scalars.py (#78): `rank1_lams`, `same_sign_form`, `closure2`, `simp_r` (Lemma S's
algebraic test and Lemma K's closure test). From phantom_source.py: `ricci`. From two_scalars.py: `rq` (random
rationals). The block formulas, Lemma R, every reduction, the planted metrics and the controls are new here.
"""
import itertools
import random
import subprocess
import sys
import time

import sympy as sp

import phantom_source as ps  # noqa: E402
import same_sign_scalars as sss  # noqa: E402
import two_scalars as ts  # noqa: E402

# name: (checks that must fail, what the mutation does)
MUTATIONS = {
    "block-Y-sign": (["block"], "builds Y as +n(b_rho)/2 + b_x b_rho"),
    "travel-offset": (["travel"], "takes u = B without + v"),
    "gauge-drop": (["travel"], "the gauge-transformed shift without - D'(t)"),
    "lemmaR-sign": (["lemmaR"], "builds the planted blocks with P = +(v0^2 + v2^2)"),
    "lemmaR-F": (["lemmaR"], "F with 2 h Z Y replaced by h Z Y"),
    "brho-Q": (["brho_solve"], "claims the rank-1 solution Q = h in place of Q = 0"),
    "brho-pert": (["brho_points"], "replaces the planted solution by b = 2 ln rho + rho^2/10"),
    "brho-bracket": (["brho_closure"], "computes the kernel commutator without the term K_2(K_1^mu)"),
    "planar-sign": (["planar_stage"], "expects lam = +P in place of -P"),
    "planar-constK": (["planar_closure"], "replaces the two planar shifts that must fail closure by b = x and b = 2 x"),
    "planar-ode": (["planar_thm"], "takes the comoving ODE x'' = -K0 x + alpha"),
    "add-quad": (["additive_solve"], "claims X'' in the radical in place of (X + R) X''"),
    "add-pert": (["additive_points"], "replaces the planted solution by u = 3/2 x + 2 ln rho + x^2/10"),
    "add-bracket": (["additive_closure"], "computes the kernel commutator without the term K_2(K_1^mu)"),
    "mult-W": (["mult"], "takes W = 2 X R + v in the claimed F"),
    "conical-Q": (["conical_entries"], "takes Q = -eta w in place of -eta w/2"),
    "conical-Eb": (["conical_solve"], "E_b with (eta^2 - 1) replaced by (eta^2 + 1)"),
    "bubble-sign": (["bubble_F"], "the claimed bracket with +z in place of -z"),
    "bubble-planted": (["bubble_points"], "replaces the block at the exact points by a planted rank-1 block"),
    "compat-Phi1": (["static_system"], "builds Phi1 with the opposite sign of a^2"),
    "lin-claim": (["linearized"], "claims the first-order part of P as +u0 w_xx"),
}

# the last section a control needs: its target checks are all in sections up to this one
STOP_AT = {"block-Y-sign": 1, "travel-offset": 1, "gauge-drop": 1, "lemmaR-sign": 2, "lemmaR-F": 2, "brho-Q": 3,
           "brho-pert": 3, "brho-bracket": 3, "planar-sign": 4, "planar-constK": 4, "planar-ode": 4, "add-quad": 5,
           "add-pert": 5, "add-bracket": 5, "mult-W": 6, "conical-Q": 7, "conical-Eb": 7, "bubble-sign": 8,
           "bubble-planted": 8, "compat-Phi1": 9, "lin-claim": 10}

RESULTS = []
MUT = None
INFO = []


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def info(text):
    INFO.append(text)
    print(f"[info] {text}", flush=True)


rq = ts.rq
simp_r = sss.simp_r
rho = sp.Symbol("rho", positive=True)
eta3 = sp.diag(-1, 1, 1)


def zero(e):
    if isinstance(e, sp.MatrixBase):
        return all(sp.simplify(x) == 0 for x in e)
    return sp.simplify(e) == 0


def blocks(b, bx, br, bt, btx, bxx, bxr, btr, brr, r):
    """(P, Q, Z, Y, h) of the block R_3 from the jets of b(t, x, rho), n = d_t - b d_x."""
    b, bx, br, bt, btx, bxx, bxr, btr, brr = map(sp.sympify, (b, bx, br, bt, btx, bxx, bxr, btr, brr))
    P = btx - b * bxx - bx**2
    Q = -(brr + br / r) / 2
    Z = bxr / 2
    Y = -(btr - b * bxr) / 2 + bx * br
    return P, Q, Z, Y, br**2 / 2


def static_blocks(u, ux, ur, uxx, uxr, urr, r):
    return blocks(u, ux, ur, 0, 0, uxx, uxr, 0, urr, r)


def block3(P, Q, Z, Y, h):
    return sp.Matrix([[P - h, Q, Z], [Q, -P - h, Y], [Z, Y, h]])


def Fof(Q, Z, Y, h):
    """F = R_01 (R_0rho^2 + R_1rho^2) + 2 R_rhorho R_0rho R_1rho in the block's entries."""
    return Q * (Z**2 + Y**2) + (1 if MUT == "lemmaR-F" else 2) * h * Z * Y


def v_block(v0, v1, v2):
    """The block of Lemma R: R_3 = lam eta_3 - v v^T with lam = h + v2^2; returns (P, Q, Z, Y, h)."""
    sgn = 1 if MUT == "lemmaR-sign" else -1
    return sgn * (v0**2 + v2**2), -v0 * v1, -v0 * v2, -v1 * v2, (v0**2 + v1**2) / 2


def rad_member(gens, f, sat, vars_):
    """f lies in the radical of <gens> : sat^infinity (Rabinowitsch, exact). sat is a list of polynomials."""
    t, s = sp.symbols("t_sat s_rab")
    extra = [1 - t * sp.Mul(*sat)] if sat else []
    Gb = sp.groebner(list(gens) + extra + [1 - s * f], *([t, s] if sat else [s]), *vars_, order="grevlex")
    return Gb.exprs == [1]


def fin(t0, mut):
    bad = [n for n, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS)} checks, {len(RESULTS) - len(bad)} pass" + (f", FAILED: {bad}" if bad else "") +
          f" ({time.time() - t0:.0f} s)")
    return bool(bad)


# 4D planar and shift metrics -------------------------------------------------------------------------------------
tt, xx, ph = sp.symbols("t x varphi", real=True)
Lg = sp.Symbol("Lg", real=True)  # Lg = ln(rho) + const, d Lg/d rho = 1/rho
Tn, Th = sp.symbols("Tn Th", real=True)  # Tn = tan t, Th = tanh t
X4 = [tt, xx, rho, ph]


def Dfun4(E, k):
    """d/dx^k with the symbols Lg = ln rho, Tn = tan t, Th = tanh t as functions of the coordinates."""
    if k == 0:
        return sp.diff(E, tt) + sp.diff(E, Tn) * (1 + Tn**2) + sp.diff(E, Th) * (1 - Th**2)
    if k == 2:
        return sp.diff(E, rho) + sp.diff(E, Lg) / rho
    return sp.diff(E, X4[k])


def shift_g(b):
    return sp.Matrix([[-1 + b**2, b, 0, 0], [b, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, rho**2]])


def K_planar(b):
    """K = -P = -(b_tx - b b_xx - b_x^2) for b(t, x), with the coordinate symbols' chain rules."""
    bx = Dfun4(b, 1)
    return sp.simplify(-(Dfun4(Dfun4(b, 0), 1) - b * Dfun4(bx, 1) - bx**2))


def exact_pt(extra=()):
    p = {rho: sp.Rational(random.randint(5, 30), 10), Lg: rq(), tt: rq(), xx: rq(), Tn: rq(), Th: rq()}
    for s_ in extra:
        p[s_] = rq()
    return p


def at(p, e):
    return sp.Rational(sp.cancel(sp.sympify(e).subs(p)))


def stage_at_points(Ric, g, lam_expr, n, seed, want_sig=None):
    """The algebraic stage at n exact points. Returns (lam_ok, none_admissible, detail): lam_ok when lam_expr is a root
    with the expected inertia (want_sig, or by the sign of lam: (2,0) for lam < 0, (0,2) for lam > 0) and pattern; none_admissible
    when no candidate root has an admissible pattern."""
    random.seed(seed)
    ok, none, det = True, True, []
    for _ in range(n):
        p = exact_pt()
        cands = sss.same_sign_form(Ric.applyfunc(lambda e: at(p, e)), g.applyfunc(lambda e: at(p, e)))
        lv = at(p, lam_expr)
        sig = want_sig or ((2, 0) if lv < 0 else (0, 2))
        pat = "PP" if sig == (2, 0) else "MM"
        ok &= any(c["lam"] == lv and c["sig"] == sig and pat in c["pats"] for c in cands)
        none &= not any(c["pats"] for c in cands)
        det.append((str(lv), [(float(sp.N(c["lam"], 5)), c["sig"], c["pats"]) for c in cands]))
    return ok, none, det


def main():
    global MUT
    opts = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
    mut = MUT = opts.get("mutate")
    if mut is not None and mut not in MUTATIONS:
        sys.exit(f"unknown mutation {mut!r}. Known: {', '.join(MUTATIONS)}")
    t0 = time.time()
    stop = int(opts.get("stop", STOP_AT.get(mut, 99)))
    print(f"one_sign_classify.py, issue #158. sympy {sp.__version__}, python {sys.version.split()[0]}. mutation: "
          f"{mut or 'none'}" + (f" ({MUTATIONS[mut][1]})" if mut else ""))
    SEC = {}

    # 1. The block, traveling = static, the gauge x -> x + D(t) ---------------------------------------------------
    print("\n== 1. The block R_3 of a generic b(t, x, rho); traveling = static; the gauge ==")
    t, x, ph_ = sp.symbols("t x phi", positive=True)
    b = sp.Function("b")(t, x, rho)
    g = sp.Matrix([[-1 + b**2, b, 0, 0], [b, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, rho**2]])
    Xg = [t, x, rho, ph_]
    Rc = ps.ricci(g, Xg, lambda E, k: sp.diff(E, Xg[k]), sp.simplify)
    E = sp.Matrix([[1, 0, 0, 0], [-b, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1 / rho]])
    Rf = (E.T * Rc * E).applyfunc(sp.simplify)
    a_ = sp.diff(b, t) - b * sp.diff(b, x)
    sgn = 1 if MUT == "block-Y-sign" else -1
    P_, Q_, Z_, Y_, h_ = blocks(b, sp.diff(b, x), sp.diff(b, rho), sp.diff(b, t), sp.diff(b, t, x), sp.diff(b, x, 2),
                                sp.diff(b, x, rho), sp.diff(b, t, rho), sp.diff(b, rho, 2), rho)
    if MUT == "block-Y-sign":
        Y_ = (sp.diff(b, t, rho) - b * sp.diff(b, x, rho)) / 2 + sp.diff(b, x) * sp.diff(b, rho)
    diffs = {"R00": Rf[0, 0] - (P_ - h_), "R11": Rf[1, 1] - (-P_ - h_), "R22": Rf[2, 2] - h_, "R01": Rf[0, 1] - Q_,
             "R02": Rf[0, 2] - Z_, "R12": Rf[1, 2] - Y_, "Rphi": sum(abs(Rf[3, j]) for j in range(4)),
             "I2": Rf[0, 0] + Rf[1, 1] + 2 * Rf[2, 2], "P=a_x": P_ - sp.diff(a_, x)}
    bad = [k for k, e in diffs.items() if not zero(e)]
    check("block", not bad, f"coordinate Ricci of an undefined b(t, x, rho) on the frame (n, d_x, d_rho, d_phi/rho) equals "
          f"[[P-h, Q, Z], [Q, -P-h, Y], [Z, Y, h]] with P = n(b_x) - b_x^2 = a_x, a = b_t - b b_x; R(e_phi, .) = 0 and "
          f"(I2) hold; nonzero differences: {bad or 'none'} ({time.time() - t0:.0f} s)")
    # traveling = static, and the gauge
    B0, Bx, Br, Bxx, Bxr, Brr, vv, D1, D2 = sp.symbols("B0 Bx Br Bxx Bxr Brr v D1 D2")
    lab = blocks(B0, Bx, Br, -vv * Bx, -vv * Bxx, Bxx, Bxr, -vv * Bxr, Brr, rho)
    ustat = B0 if MUT == "travel-offset" else B0 + vv
    stat = static_blocks(ustat, Bx, Br, Bxx, Bxr, Brr, rho)
    d_trav = [sp.expand(p - q) for p, q in zip(lab, stat)]
    b0, bx, br, bt, btx, bxx, bxr, btr, brr = sp.symbols("b0 bx br bt btx bxx bxr btr brr")
    orig = blocks(b0, bx, br, bt, btx, bxx, bxr, btr, brr, rho)
    dd = 0 if MUT == "gauge-drop" else D1
    gau = blocks(b0 - dd, bx, br, bt - D1 * bx - D2, btx - D1 * bxx, bxx, bxr, btr - D1 * bxr, brr, rho)
    d_gauge = [sp.expand(p - q) for p, q in zip(gau, orig)]
    check("travel", all(e == 0 for e in d_trav + d_gauge), f"(P, Q, Z, Y, h) of b = B(x - v t, rho) minus those of the static "
          f"u = B + v: {d_trav}; of b~(t, x) = b(t, x - D(t)) - D'(t) minus those of b at (t, x - D): {d_gauge}. So a "
          f"traveling shift is a static one, and the classification is up to the gauge x -> x + D(t)")
    if stop < 2:
        return fin(t0, mut)

    # 2. Lemma R ----------------------------------------------------------------------------------------------------
    print("\n== 2. Lemma R: rank(R_3 - lam eta_3) = 1 <=> N = -v v^T (h > 0), planted forms ==")
    v0, v1, v2 = sp.symbols("v0 v1 v2")
    lam_s = sp.Symbol("lam_s")
    P, Q, Z, Y, h = v_block(v0, v1, v2)
    N = block3(P, Q, Z, Y, h) - (h + v2**2) * eta3
    ra = sp.expand(N + sp.Matrix([v0, v1, v2]) * sp.Matrix([v0, v1, v2]).T)
    ok_a = all(e == 0 for e in ra)
    ok_b = sp.expand(Fof(Q, Z, Y, h)) == 0
    R00, R01, R0r, R1r, Rrr, lm = sp.symbols("R00 R01 R0r R1r Rrr lm")
    R11 = -R00 - 2 * Rrr  # (I2)
    Rg = sp.Matrix([[R00, R01, R0r], [R01, R11, R1r], [R0r, R1r, Rrr]]) - lm * eta3
    m_d = Rg[0, 0] * Rg[1, 2] - Rg[0, 1] * Rg[0, 2]
    m_e = Rg[1, 1] * Rg[0, 2] - Rg[0, 1] * Rg[1, 2]
    cert = sp.expand(R0r * m_d + R1r * m_e + Fof(R01, R0r, R1r, Rrr))
    ok_c = cert == 0
    random.seed(1581)
    bad_d = []
    for _ in range(16):
        w3 = [sp.Rational(random.choice([-3, -2, -1, 1, 2, 3]), random.randint(1, 3)) for _ in range(3)]
        Pb = v_block(*w3)
        res = sss.rank1_lams(block3(*Pb))
        want = [(Pb[4] + w3[2]**2, "MM")]
        if res != want:
            bad_d.append(("v", w3, res))
    for k in range(9):
        Pp = 0 if k == 8 else sp.Rational(random.choice([-3, -2, -1, 1, 2, 3]), random.randint(1, 3))
        res = sss.rank1_lams(block3(Pp, 0, 0, 0, 0))
        want = [] if Pp == 0 else [(-Pp, "PP" if Pp > 0 else "MM")]
        if res != want:
            bad_d.append(("planar", Pp, res))
    for _ in range(16):
        Pg, Qg, Zg, Yg, hg = [rq() for _ in range(4)] + [abs(rq()) + sp.Rational(1, 7)]
        if Fof(Qg, Zg, Yg, hg) == 0:
            continue
        res = sss.rank1_lams(block3(Pg, Qg, Zg, Yg, hg))
        if res != []:
            bad_d.append(("generic", (Pg, Qg, Zg, Yg, hg), res))
    check("lemmaR", ok_a and ok_b and ok_c and not bad_d,
          f"(a) R_3(v) - (h + v2^2) eta_3 + v v^T = 0 identically: {ok_a}; (b) F = 0 on it: {ok_b}; (c) the certificate "
          f"R_0rho m_d + R_1rho m_e + F = 0 for a generic R_3 with (I2): {ok_c}; (d) planted: 16 random v-blocks give "
          f"exactly [(h + v2^2, MM)], 8 planar blocks [(-P, PP if P > 0 else MM)], P = 0 gives [], 16 generic blocks "
          f"give []: wrong in {len(bad_d)} {bad_d[:2]}")
    if stop < 3:
        return fin(t0, mut)

    # 3. b(rho) -----------------------------------------------------------------------------------------------------
    print("\n== 3. b(rho): the planted solution b = C ln rho + b_0 and Lemma E ==")
    bf = sp.Function("bf")(rho)
    jb = static_blocks(bf, 0, sp.diff(bf, rho), 0, 0, sp.diff(bf, rho, 2), rho)
    hh = sp.diff(bf, rho)**2 / 2
    Qb = -(sp.diff(bf, rho, 2) + sp.diff(bf, rho) / rho) / 2
    blk_ok = zero(block3(*jb) - sp.Matrix([[-hh, Qb, 0], [Qb, -hh, 0], [0, 0, hh]]))
    Qs, hs = sp.symbols("Qs hs")
    gens = [v0**2 + v2**2, Qs + v0 * v1, v0 * v2, v1 * v2, 2 * hs - v0**2 - v1**2]  # Lemma R with P = Z = Y = 0
    claim = Qs - hs if MUT == "brho-Q" else Qs
    vars_ = (v0, v1, v2, Qs, hs)
    m_v0, m_v2 = rad_member(gens, v0, [], vars_), rad_member(gens, v2, [], vars_)
    m_Q = rad_member(gens, claim, [], vars_)
    sol = sp.dsolve(sp.Eq(sp.diff(bf, rho, 2) + sp.diff(bf, rho) / rho, 0), bf)
    C1, C2 = sp.symbols("C1 C2")
    sol_ok = sp.simplify(sol.rhs - (C1 + C2 * sp.log(rho))) == 0
    check("brho_solve", blk_ok and m_v0 and m_v2 and m_Q and sol_ok,
          f"R_3 of b(rho) = [[-h, Q, 0], [Q, -h, 0], [0, 0, h]]: {blk_ok}; with P = Z = Y = 0 in Lemma R's equations v0 and v2 "
          f"lie in the radical ({m_v0}, {m_v2}) and {claim} does ({m_Q}), so lam = h and Q = 0 (Lemma E: the same for a "
          f"limit of a localized shift); dsolve(b'' + b'/rho = 0) = {sol.rhs}")
    D3b = sp.Rational(1, 3) + 2 * Lg
    Rb = ps.ricci(shift_g(D3b), X4, Dfun4, simp_r)
    b_pass = D3b + (rho**2 / 10 if MUT == "brho-pert" else 0)
    Rpass = ps.ricci(shift_g(b_pass), X4, Dfun4, simp_r) if MUT == "brho-pert" else Rb
    lam_b = sp.diff(b_pass, rho) + sp.diff(b_pass, Lg) / rho
    lam_b = lam_b**2 / 2
    ok1, _, d1 = stage_at_points(Rpass, shift_g(b_pass), lam_b, 3, 1582, (0, 2))
    bpert = D3b + rho**2 / 10
    Rp = ps.ricci(shift_g(bpert), X4, Dfun4, simp_r)
    _, none2, d2 = stage_at_points(Rp, shift_g(bpert), 0, 3, 1583)
    check("brho_points", ok1 and none2, f"b = 2 ln rho + 1/3, 4D, 3 exact points: lam = h is a root with inertia (0,2) and "
          f"pattern MM: {ok1} {d1[:1]}; planted failure b = 2 ln rho + rho^2/10: no admissible candidate at 3 points: "
          f"{none2} ({time.time() - t0:.0f} s)")
    if MUT == "brho-bracket":
        sss.MUT = "d3-bracket"
    lam3 = (sp.diff(D3b, rho) + sp.diff(D3b, Lg) / rho)**2 / 2
    M3 = (Rb - lam3 * shift_g(D3b)).applyfunc(simp_r)
    cl3 = sss.closure2(M3, Dfun4, simp_r, lam3)
    sss.MUT = None
    check("brho_closure", cl3["dim"] == 2 and not cl3["integrable"],
          f"closure2 of M = R - h g for b = 2 ln rho + 1/3: kernel dimension {cl3['dim']}, [K_1, K_2] = "
          f"{[str(x_) for x_ in cl3.get('bracket', [])]}, plane integrable: {cl3['integrable']} (Frobenius: no closed pair "
          f"exists) ({time.time() - t0:.0f} s)")
    if stop < 4:
        return fin(t0, mut)

    # 4. planar b(t, x) ---------------------------------------------------------------------------------------------
    print("\n== 4. Planar shifts b(t, x): lam = K = -P, closure iff K constant, Theorem P ==")
    bp = sp.Function("bp")(t, x)
    jp = blocks(bp, sp.diff(bp, x), 0, sp.diff(bp, t), sp.diff(bp, t, x), sp.diff(bp, x, 2), 0, 0, 0, rho)
    Pp_ = sp.diff(sp.diff(bp, t) - bp * sp.diff(bp, x), x)
    blk_p = zero(block3(*jp) - sp.diag(Pp_, -Pp_, 0)) and zero(jp[0] - Pp_)
    Pl, lm_ = sp.symbols("Pl lm_")
    Nn = sp.diag(Pl + lm_, -Pl - lm_, -lm_)
    mins = [Nn[i, k] * Nn[j, l] - Nn[i, l] * Nn[j, k] for (i, j) in itertools.combinations(range(3), 2)
            for (k, l) in itertools.combinations(range(3), 2)]
    sols = sp.solve([m for m in mins if m != 0], lm_, dict=True)
    claimed_lam = Pl if MUT == "planar-sign" else -Pl
    lam_ok = any(sp.simplify(s_[lm_] - claimed_lam) == 0 for s_ in sols) and len(sols) == 1
    planar_metrics = [("3/2 x", sp.Rational(3, 2) * xx), ("tan(t) x", Tn * xx), ("-tanh(t) x", -Th * xx)]
    allok, det = True, []
    for k, (nm, bpl) in enumerate(planar_metrics):
        Kp = K_planar(bpl)
        Rp_ = ps.ricci(shift_g(bpl), X4, Dfun4, simp_r)
        lam_e = -Kp if MUT == "planar-sign" else Kp
        ok_k, _, d_k = stage_at_points(Rp_, shift_g(bpl), lam_e, 3, 1590 + k)
        allok &= ok_k
        det.append((nm, str(Kp), ok_k))
    check("planar_stage", blk_p and lam_ok and allok,
          f"R_3 of a generic planar b = diag(P, -P, 0), P = a_x: {blk_p}; the 2 x 2 minors of diag(P, -P, 0) - lam eta_3 "
          f"vanish only at lam = {claimed_lam} (solve: {[s_[lm_] for s_ in sols]}): {lam_ok}; 4D, 3 exact points each, "
          f"lam = K with inertia (2,0) for K < 0 and (0,2) for K > 0: {det} ({time.time() - t0:.0f} s)")
    cls = []
    good_ok = True
    for nm, bpl in planar_metrics:
        Kp = K_planar(bpl)
        Rp_ = ps.ricci(shift_g(bpl), X4, Dfun4, simp_r)
        Mp = (Rp_ - Kp * shift_g(bpl)).applyfunc(simp_r)
        cl = sss.closure2(Mp, Dfun4, simp_r, Kp)
        okc = cl["dim"] == 2 and cl["integrable"] and cl["invariant"] and cl["flat"]
        good_ok &= okc
        cls.append((nm, str(Kp), cl["dim"], cl["integrable"], cl["invariant"], cl["flat"]))
    fails = [("x^2/10", sp.Rational(1, 10) * xx**2), ("1/(1 + x^2)", 1 / (1 + xx**2))]
    if MUT == "planar-constK":
        fails = [("x", xx), ("2 x", 2 * xx)]
    bad_ok = True
    clf = []
    for nm, bpl in fails:
        Kp = K_planar(bpl)
        Rp_ = ps.ricci(shift_g(bpl), X4, Dfun4, simp_r)
        Mp = (Rp_ - Kp * shift_g(bpl)).applyfunc(simp_r)
        cl = sss.closure2(Mp, Dfun4, simp_r, Kp)
        okf = cl["dim"] == 2 and cl["integrable"] and not cl["invariant"]
        bad_ok &= okf
        clf.append((nm, str(Kp), cl["dim"], cl["integrable"], cl["invariant"]))
    check("planar_closure", good_ok and bad_ok,
          f"closure2 of M = R - K g, (name, K, dim, integrable, L_X M = 0, flat): must pass for K constant: {cls}; must fail "
          f"on L_X M = 0 (K not constant) with the kernel integrable: {clf} ({time.time() - t0:.0f} s)")
    # Theorem P: comoving curves
    bt_, bx_, bb, K0, al, xs = sp.symbols("bt_ bx_ bb K0 al xs")
    a_val = -K0 * xs + al
    xpp = -(bt_ + bx_ * (-bb))  # x'' along dx/dt = -b, with x' = -b
    ode_ok = sp.expand(xpp + (bt_ - bb * bx_)) == 0 and sp.expand(-a_val - ((-K0 * xs + al) if MUT == "planar-ode" else (K0 * xs - al))) == 0
    kap, tq, t1, Xv, bet, xp1, xp1d = sp.symbols("kappa tq t1 X beta xp1 xp1d", real=True)
    kap = sp.Symbol("kappa", positive=True)
    xpf, alf = sp.Function("xp")(t), sp.Function("al")(t)
    res = []
    for sg, cs, sn, rat in ((1, sp.cosh, sp.sinh, kap * sp.tanh(kap * (tq - t1))), (-1, sp.cos, sp.sin, -kap * sp.tan(kap * (tq - t1)))):
        xX = xpf + (Xv - xp1) * cs(kap * (t - t1)) + ((bet - xp1d) / kap) * sn(kap * (t - t1))
        ex = sp.diff(xX, t, 2) - (sg * kap**2 * xX - alf)
        ex = ex.subs(sp.Derivative(xpf, (t, 2)), sg * kap**2 * xpf - alf)
        slope = sp.diff(xX, Xv, t).subs(t, tq) / sp.diff(xX, Xv).subs(t, tq)
        res.append(sp.simplify(ex) == 0 and sp.simplify(slope - rat) == 0)
    check("planar_thm", ode_ok and all(res),
          f"along dx/dt = -b the chain rule gives x'' = -a, so with a = -K0 x + alpha(t) (i.e. K = -a_x = K0 constant) "
          f"x'' = K0 x - alpha: {ode_ok}; the family x_X(t) through (t1, X) solves it for K0 = +kappa^2 and -kappa^2, and "
          f"d_X x'/d_X x = kappa tanh(kappa (t* - t1)) resp. -kappa tan(kappa (t* - t1)): {res}. So -b(t*, .) = x' grows "
          f"linearly in the comoving position: no bounded b, so no localized one (Theorem P, by hand)")
    if stop < 5:
        return fin(t0, mut)

    # 5. additive u = X(x) + R(rho) (+ v) ---------------------------------------------------------------------------
    print("\n== 5. Additive shifts u = X(x) + R(rho) + v: u = c x + C ln rho + d, closure ==")
    Xf, Rf = sp.Function("Xf")(x), sp.Function("Rf")(rho)
    uA = Xf + Rf + vv
    jA = static_blocks(uA, sp.diff(Xf, x), sp.diff(Rf, rho), sp.diff(Xf, x, 2), 0, sp.diff(Rf, rho, 2), rho)
    SA = sp.diff(Rf, rho, 2) + sp.diff(Rf, rho) / rho
    claimA = (-uA * sp.diff(Xf, x, 2) - sp.diff(Xf, x)**2, -SA / 2, 0, sp.diff(Xf, x) * sp.diff(Rf, rho),
              sp.diff(Rf, rho)**2 / 2)
    blkA = all(zero(p - q) for p, q in zip(jA, claimA))
    Uu, X1, X2, R1, Sx = sp.symbols("Uu X1 X2 R1 Sx")
    gensA = [-Uu * X2 - X1**2 + v0**2 + v2**2, -Sx / 2 + v0 * v1, v0 * v2, X1 * R1 + v1 * v2, R1**2 - v0**2 - v1**2]
    varsA = (v0, v1, v2, Uu, X1, X2, R1, Sx)
    mS = rad_member(gensA, Sx, [R1, X1], varsA)
    mv0 = rad_member(gensA, v0, [R1, X1], varsA)
    claimQ = X2 if MUT == "add-quad" else Uu * X2
    mU = rad_member(gensA, claimQ, [R1, X1], varsA)
    mX2 = rad_member(gensA, X2, [R1, X1], varsA)
    check("additive_block", blkA, "R_3 of u = X(x) + R(rho) + v equals Z = 0, Y = X' R', Q = -S/2, P = -u X'' - X'^2, "
          f"h = R'^2/2, S = R'' + R'/rho: {blkA}")
    check("additive_solve", mS and mv0 and mU,
          f"with X' R' != 0: S = R'' + R'/rho, v0 and {claimQ} lie in the radical of Lemma R's ideal saturated by R' X' "
          f"({mS}, {mv0}, {mU}); X'' alone does not ({mX2}, information). By hand: S = 0 gives R = C ln rho + R_0, u X'' = 0 "
          f"with u = X + R + v varying in rho gives X'' = 0, so X = c x + d_0 and u = c x + C ln rho + d; X' = 0 on an open set "
          f"is b(rho) (Section 3); R' = 0 is planar (Section 4)")
    uA1 = sp.Rational(3, 2) * xx + 2 * Lg + sp.Rational(1, 5)
    uA2 = xx + 3 * Lg
    uA3 = sp.Rational(3, 2) * xx + 2 * Lg + xx**2 / 10
    uA0 = xx
    lamA = lambda u_: Dfun4(u_, 2)**2 / 2 + Dfun4(u_, 1)**2
    RA = {}
    for nm, u_ in (("A1", uA1), ("A2", uA2), ("A3", uA3), ("A0", uA0)):
        RA[nm] = ps.ricci(shift_g(u_), X4, Dfun4, simp_r)
    okp = True
    for k, (nm, u_) in enumerate((("A1", uA3 if MUT == "add-pert" else uA1), ("A2", uA3 if MUT == "add-pert" else uA2))):
        okk, _, _ = stage_at_points(RA["A3" if MUT == "add-pert" else nm], shift_g(u_), lamA(u_), 3, 1600 + k, (0, 2))
        okp &= okk
    _, none3, d3_ = stage_at_points(RA["A3"], shift_g(uA3), 0, 3, 1603)
    check("additive_points", okp and none3,
          f"u = 3/2 x + 2 ln rho + 1/5 and u = x + 3 ln rho, 4D, 3 exact points each: lam = C^2/(2 rho^2) + c^2 is a root "
          f"with inertia (0,2) and pattern MM: {okp}; planted failure u = 3/2 x + 2 ln rho + x^2/10: no admissible candidate "
          f"at 3 points: {none3} ({time.time() - t0:.0f} s)")
    clA = []
    okc = True
    for nm, u_, want_int in (("A1", uA1, False), ("A2", uA2, False)):
        lam_u = lamA(u_)
        if MUT == "add-bracket":
            sss.MUT = "d3-bracket"
        cl = sss.closure2((RA[nm] - lam_u * shift_g(u_)).applyfunc(simp_r), Dfun4, simp_r, lam_u)
        sss.MUT = None
        okc &= cl["dim"] == 2 and not cl["integrable"]
        clA.append((nm, cl["dim"], cl["integrable"], [str(e_) for e_ in cl.get("bracket", [])]))
    cl0 = sss.closure2((RA["A0"] - lamA(uA0) * shift_g(uA0)).applyfunc(simp_r), Dfun4, simp_r, lamA(uA0))
    ok0 = cl0["dim"] == 2 and cl0["integrable"] and cl0["invariant"] and cl0["flat"]
    check("additive_closure", okc and ok0,
          f"closure2 of M = R - lam g: u = 3/2 x + 2 ln rho + 1/5 and u = x + 3 ln rho have a 2D kernel that is not integrable "
          f"(name, dim, integrable, [K_1, K_2]): {clA}; for C = 0, u = x: integrable {cl0['integrable']}, L_X M = 0 "
          f"{cl0['invariant']}, quotient flat {cl0['flat']} (the planar K = c^2, dS_2 x R^2) ({time.time() - t0:.0f} s)")
    if stop < 6:
        return fin(t0, mut)

    # 6. multiplicative u = X(x) R(rho) + v -------------------------------------------------------------------------
    print("\n== 6. Multiplicative shifts u = X(x) R(rho) + v: F = 0 alone leaves only X' R' = 0 ==")
    Xm, Xm1, Xm2, Rm, Rm1, Rm2, Sm_, tm = sp.symbols("Xm Xm1 Xm2 Rm Rm1 Rm2 Sm tm")
    um = Xm * Rm + vv
    Pm, Qm, Zm, Ym, hm = static_blocks(um, Xm1 * Rm, Xm * Rm1, Xm2 * Rm, Xm1 * Rm1, Xm * Rm2, rho)
    Fm = Fof(Qm, Zm, Ym, hm)
    Sm = Rm2 + Rm1 / rho
    Wm = (2 if MUT == "mult-W" else 3) * Xm * Rm + vv
    claimFm = -Xm * Xm1**2 * Rm1**2 * (Sm * (1 + Wm**2) - 2 * Xm * Rm1**2 * Wm) / 8
    ident_m = sp.expand(Fm - claimFm) == 0
    Wt = 3 * Xm * Rm + vv
    Et = Sm_ * (1 + Wt**2) - 2 * Xm * Rm1**2 * Wt
    cf = sp.Poly(sp.expand(Et), Xm).all_coeffs()
    Gm = sp.groebner(cf + [1 - tm * Rm * Rm1], tm, Rm, Rm1, Sm_, vv, order="grevlex")
    check("mult", ident_m and Gm.exprs == [1],
          f"F of u = X R + v equals -X X'^2 R'^2 (S (1 + W^2) - 2 X R'^2 W)/8, W = 3 X R + v, S = R'' + R'/rho: {ident_m}; the "
          f"three coefficients of X^2, X, 1 of S (1 + W^2) - 2 X R'^2 W with 1 - t R R' generate the unit ideal: "
          f"{Gm.exprs == [1]}, so with X' R' X != 0 and X ranging over an interval F = 0 has no solution; the remaining "
          f"cases X' = 0 (b(rho)) and R' = 0 (planar) are Sections 3 and 4")
    if stop < 7:
        return fin(t0, mut)

    # 7. conical u = F(x/rho) + v -----------------------------------------------------------------------------------
    print("\n== 7. Conical shifts u = F(x/rho) + v (static and traveling): no solution with F' != 0 ==")
    et, Uv, Gv, Hv, x_, sat = sp.symbols("eta U G H x_ sat")
    et0 = sp.Symbol("eta0")
    xq = sp.Symbol("xq", positive=True)
    uq = Uv + Gv * (xq / rho - et0) + Hv * (xq / rho - et0)**2 / 2  # jets at x = eta0 rho of F(x/rho) + v, U = F + v
    jq = static_blocks(uq, sp.diff(uq, xq), sp.diff(uq, rho), sp.diff(uq, xq, 2), sp.diff(uq, xq, rho),
                       sp.diff(uq, rho, 2), rho)
    jq = [sp.simplify(e_.subs(xq, et0 * rho) * rho**2) for e_ in jq]
    wq = Gv + et0 * Hv
    Qfac = 1 if MUT == "conical-Q" else sp.Rational(1, 2)
    claimq = [-(Gv**2 + Uv * Hv), -et0 * wq * Qfac, -wq / 2, -(et0 * Gv**2 + Uv * wq / 2), et0**2 * Gv**2 / 2]
    ent_ok = all(sp.simplify(p - q) == 0 for p, q in zip(jq, claimq))
    check("conical_entries", ent_ok and sp.simplify(jq[1] - et0 * jq[2]) == 0,
          f"rho^2 (P, Q, Z, Y, h) of F(x/rho) + v at x = eta rho equals (-(G^2 + U H), -eta w Q', -w/2, -(eta G^2 + U w/2), "
          f"eta^2 G^2/2), w = G + eta H, Q' = {Qfac}, from the jets of a quadratic F: {ent_ok}; Q = eta Z (hand claim) holds "
          f"for the computed entries: {sp.simplify(jq[1] - et0 * jq[2]) == 0}")
    Ph_ = -(Gv**2 + Uv * Hv)
    wv = Gv + et * Hv
    Qh_, Zh_, Yh_, hh_ = -et * wv / 2, -wv / 2, -(et * Gv**2 + Uv * wv / 2), et**2 * Gv**2 / 2
    gensC = [Ph_ + v0**2 + v2**2, Qh_ + v0 * v1, Zh_ + v0 * v2, Yh_ + v1 * v2, 2 * hh_ - v0**2 - v1**2]
    Ea = wv * (1 + Uv**2) + 2 * et * Uv * Gv**2
    Eb = (et**2 + 1) * et * Hv + ((et**2 + 1) if MUT == "conical-Eb" else (et**2 - 1)) * Gv
    GbC = sp.groebner(gensC + [1 - sat * et * Gv * Uv * wv], v0, v1, v2, sat, Hv, Uv, Gv, et, order="grevlex")
    redA, redB = GbC.reduce(sp.expand(Ea))[1], GbC.reduce(sp.expand(Eb))[1]
    Gb0 = sp.groebner(gensC + [wv, 1 - sat * et * Gv * Uv], v0, v1, v2, sat, Hv, Uv, Gv, et, order="grevlex")
    Gf = sp.Function("Gf")(et)
    solG = sp.dsolve(sp.Eq(Gf.diff(et), -(et**2 - 1) / (et * (et**2 + 1)) * Gf), Gf)
    Cc = sp.Symbol("Cc")
    G_sol = Cc * et / (et**2 + 1)
    ode_ok = sp.simplify(solG.rhs.subs(sp.Symbol("C1"), Cc) - G_sol) == 0
    Fs = sp.Symbol("Fs")
    Hsol = sp.diff(G_sol, et)
    Ea_on = Ea.subs({Hv: Hsol, Gv: G_sol, Uv: Fs})
    fac_ok = sp.simplify(Ea_on - (2 * G_sol / (et**2 + 1)) * (1 + Fs**2 + Cc * et**2 * Fs)) == 0
    dPhi = 2 * Fs * G_sol + 2 * Cc * et * Fs + Cc * et**2 * G_sol
    Fsol = sp.solve(dPhi, Fs)
    Fsol = Fsol[0] if len(Fsol) == 1 else None
    mism = sp.factor(sp.diff(Fsol, et) - G_sol) if Fsol is not None else None
    want_mism = -Cc * et * ((et**2 + 2)**2 + 2 * (et**2 + 1)) / ((et**2 + 1) * (et**2 + 2)**2)
    yq = sp.Symbol("yq")  # y = eta^2 >= 0: the numerator (y + 2)^2 + 2 (y + 1) has positive coefficients
    mism_ok = Fsol is not None and sp.simplify(mism - want_mism) == 0 and all(
        c_ > 0 for c_ in sp.Poly((yq + 2)**2 + 2 * (yq + 1), yq).all_coeffs())
    check("conical_solve", redA == 0 and redB == 0 and Gb0.exprs == [1] and ode_ok and fac_ok and mism_ok,
          f"Lemma R's equations in these jets, saturated by eta G U w: E_a reduces to {redA}, E_b to {redB} modulo the Groebner "
          f"basis ({len(GbC.exprs)} elements); the branch w = 0 with 1 - t eta G U is the unit ideal: {Gb0.exprs == [1]}; "
          f"E_b is G' = -(eta^2 - 1) G/(eta (eta^2 + 1)) with dsolve {solG.rhs}: {ode_ok}; then E_a = (2G/(eta^2 + 1)) (1 + "
          f"F^2 + C eta^2 F): {fac_ok}; d/d eta of that bracket vanishes only at F = {Fsol}, whose F' - G = {mism} = "
          f"-C eta (eta^4 + 6 eta^2 + 6)/(...) != 0 for C eta != 0: {mism_ok}. So no solution with F' != 0 ({time.time() - t0:.0f} s)")
    if stop < 8:
        return fin(t0, mut)

    # 8. the spherical profile u = w(r): the bubble at constant speed (planted, #78) ---------------------------------
    print("\n== 8. Spherical profile u = w(r), r = sqrt(x^2 + rho^2): the bubble at constant speed passes on no open set ==")
    r_, s_, c_, W0, W1, W2, tb = sp.symbols("r s c W0 W1 W2 tb")

    def sph_jets(W0_, W1_, W2_, r0, s0, c0):
        return (W0_, W1_ * c0, W1_ * s0, W2_ * c0**2 + W1_ * s0**2 / r0, (W2_ - W1_ / r0) * s0 * c0,
                W2_ * s0**2 + W1_ * c0**2 / r0)

    xs_, ys_ = sp.symbols("xs_ ys_", positive=True)
    wex = lambda R_: R_**3 + 2 * R_**2 + R_
    Rr_ = sp.sqrt(xs_**2 + ys_**2)
    ue = wex(Rr_)
    pt = {xs_: 3, ys_: 4}
    Wd = (wex(5), sp.diff(wex(sp.Symbol("q")), sp.Symbol("q")).subs(sp.Symbol("q"), 5),
          sp.diff(wex(sp.Symbol("q")), sp.Symbol("q"), 2).subs(sp.Symbol("q"), 5))
    jex = [ue, sp.diff(ue, xs_), sp.diff(ue, ys_), sp.diff(ue, xs_, 2), sp.diff(ue, xs_, ys_), sp.diff(ue, ys_, 2)]
    jex = [sp.simplify(e_.subs(pt)) for e_ in jex]
    jfor = sph_jets(*Wd, 5, sp.Rational(4, 5), sp.Rational(3, 5))
    jets_ok = all(sp.simplify(p - q) == 0 for p, q in zip(jex, jfor))
    uu_, ux_, ur_, uxx_, uxr_, urr_ = sph_jets(W0, W1, W2, r_, s_, c_)
    Pb_, Qb_, Zb_, Yb_, hb_ = static_blocks(uu_, ux_, ur_, uxx_, uxr_, urr_, r_ * s_)
    Fb = Fof(Qb_, Zb_, Yb_, hb_)
    a1 = r_ * W2 - W1
    a2 = W0 * a1 + 2 * r_ * W1**2
    zz = 2 * W1 + a1 * s_**2
    sz = 1 if MUT == "bubble-sign" else -1
    Bb = sz * zz * (a1**2 + a2**2) + 2 * r_ * a1 * a2 * W1**2 * s_**2
    diffb = sp.numer(sp.together(Fb - s_**2 * c_**2 * Bb / (8 * r_**3)))
    ident_b = sp.rem(sp.expand(diffb), c_**2 + s_**2 - 1, c_) == 0
    Bt = -(2 * W1 + a1 * s_**2) * (a1**2 + a2**2) + 2 * r_ * a1 * a2 * W1**2 * s_**2
    c0b, c1b = sp.Poly(sp.expand(Bt), s_).coeff_monomial(1), sp.Poly(sp.expand(Bt), s_).coeff_monomial(s_**2)
    Gbb = sp.groebner([c0b, c1b, 1 - tb * W1 * r_], tb, r_, W0, W1, W2, order="grevlex")
    check("bubble_F", jets_ok and ident_b and Gbb.exprs == [1],
          f"chain-rule jets of w(r) against a sample w at (x, rho) = (3, 4): {jets_ok}; 8 r^3 F = s^2 c^2 (-z (a1^2 + a2^2) + 2 r "
          f"a1 a2 w'^2 s^2) with a1 = r w'' - w', a2 = w a1 + 2 r w'^2, z = 2 w' + a1 s^2 (#78's B): {ident_b}; the s^0 and s^2 "
          f"coefficients of the bracket with 1 - t w' r generate the unit ideal: {Gbb.exprs == [1]} ({time.time() - t0:.0f} s)")
    random.seed(1581)
    found = []
    for _ in range(20):
        m = sp.Rational(random.choice([-5, -3, -2, 2, 3, 4, 5]), random.choice([1, 2, 3]))
        s0, c0 = 2 * m / (1 + m**2), (1 - m**2) / (1 + m**2)
        r0 = sp.Rational(random.randint(1, 9), random.randint(1, 3))
        w0_, w1_, w2_ = rq(), rq() + 2, rq()
        if MUT == "bubble-planted":
            vv3 = [sp.Rational(random.choice([-3, -2, -1, 1, 2, 3]), random.randint(1, 3)) for _ in range(3)]
            Bk = v_block(*vv3)
            res = sss.rank1_lams(block3(*Bk))
        else:
            j6 = sph_jets(w0_, w1_, w2_, r0, s0, c0)
            Bk = static_blocks(j6[0], j6[1], j6[2], j6[3], j6[4], j6[5], r0 * s0)
            res = sss.rank1_lams(block3(*Bk))
        found.append(res)
    check("bubble_points", all(res == [] for res in found),
          f"20 exact points (random w, w' != 0, w'', r and a rational angle with s c != 0): rank1_lams of the block is empty at "
          f"{sum(res == [] for res in found)} of 20")
    if stop < 9:
        return fin(t0, mut)

    # 9. general static shifts with Z != 0 ---------------------------------------------------------------------------
    print("\n== 9. General static shifts u(x, rho) with Z != 0: the resolved system and its compatibility condition ==")
    Pq, Qq, Zq, Yq, hq = sp.symbols("Pq Qq Zq Yq hq")
    gensR = [Pq + v0**2 + v2**2, Qq + v0 * v1, Zq + v0 * v2, Yq + v1 * v2, 2 * hq - v0**2 - v1**2]
    varsR = (v0, v1, v2, Pq, Qq, Zq, Yq, hq)
    A_ = Qq * (Zq**2 + Yq**2) + 2 * hq * Zq * Yq
    B_ = Pq * Qq * Yq - Zq * (Qq**2 + Yq**2)
    mA_ = rad_member(gensR, A_, [Zq, Yq, Qq], varsR)
    mB_ = rad_member(gensR, B_, [Zq, Yq, Qq], varsR)
    us, as_, bs, ss, sx, sr = sp.symbols("us as_ bs ss sx sr")
    cth, sth = (1 - ss**2) / (1 + ss**2), 2 * ss / (1 + ss**2)  # ss = tan(theta/2), v0 = beta cos, v1 = beta sin
    w0, w1, w2 = bs * cth, bs * sth, as_ / (us * cth - sth)
    Phi1 = (w0**2 + w2**2 + (1 if MUT == "compat-Phi1" else -1) * as_**2) / us
    Psi = -2 * w0 * w2
    Phi3 = -bs / rho + 2 * w0 * w1
    Pent = static_blocks(us, as_, bs, Phi1, Psi, Phi3, rho)
    ident_p = all(sp.simplify(p - q) == 0 for p, q in zip(Pent, v_block(w0, w1, w2)))

    def Dx(f, S):
        return sp.diff(f, us) * as_ + sp.diff(f, as_) * Phi1 + sp.diff(f, bs) * Psi + sp.diff(f, ss) * S

    def Dr(f, S):
        return sp.diff(f, rho) + sp.diff(f, us) * bs + sp.diff(f, as_) * Psi + sp.diff(f, bs) * Phi3 + sp.diff(f, ss) * S

    E1 = Dr(Phi1, sr) - Dx(Psi, sx)
    E2 = Dr(Psi, sr) - Dx(Phi3, sx)
    E1z, E2z = E1.subs({sx: 0, sr: 0}), E2.subs({sx: 0, sr: 0})
    Am = sp.Matrix([[sp.diff(E1, sx), sp.diff(E1, sr)], [sp.diff(E2, sx), sp.diff(E2, sr)]]).applyfunc(sp.cancel)
    Bm = sp.Matrix([E1z, E2z]).applyfunc(sp.cancel)
    sig = (-Am.inv() * Bm).applyfunc(sp.cancel)
    Ccomp = Dr(sig[0], sig[1]) - Dx(sig[1], sig[0])  # [X_x, X_rho] = C d_s
    random.seed(1591)
    pl = []
    for _ in range(6):
        p = {rho: sp.Rational(random.randint(1, 9), random.randint(1, 4)), us: rq(), as_: rq(), bs: rq()}
        p = {k: (v_ if v_ != 0 else sp.Integer(1)) for k, v_ in p.items()}
        pl.append(p)
    e12 = all(sp.cancel(e_.subs(ss, 1)) == 0 for e_ in (E1z, E2z))
    cvals = [sp.nsimplify(sp.cancel(Ccomp.subs({**p, ss: 1}))) for p in pl]
    check("static_system", mA_ and mB_ and ident_p and e12 and all(c_ == 0 for c_ in cvals),
          f"with h > 0 and Z Y Q != 0, A = Q (Z^2 + Y^2) + 2 h Z Y and B = P Q Y - Z (Q^2 + Y^2) lie in the radical of Lemma R's "
          f"ideal ({mA_}, {mB_}); with v0 = beta cos th, v1 = beta sin th, v2 = a/(u cos th - sin th) the block built from "
          f"u_xx = Phi1, u_x rho = Psi, u_rho rho = Phi3 equals Lemma R's: {ident_p}; planted, on the jets of u = c x + C ln rho + d "
          f"(s = tan(th/2) = 1): E1 = E2 = 0 with sigma = 0: {e12}; C = {cvals} at 6 exact points ({time.time() - t0:.0f} s)")
    random.seed(1592)
    gen_pts, dets, cv = [], [], []
    redraws = 0
    for s_val in (sp.Rational(1, 2), 2, 3):
        while True:  # positive rationals; a draw on a pole of the expressions (non-finite value) is redrawn
            p = {k_: sp.Rational(random.randint(1, 9), random.randint(1, 4)) for k_ in (rho, us, as_, bs)}
            p[ss] = s_val
            d_, c_ = sp.cancel(Am.det().subs(p)), sp.cancel(Ccomp.subs(p))
            if d_.is_Rational and c_.is_Rational:
                break
            redraws += 1
        dets.append(d_)
        cv.append(c_)
        gen_pts.append(p)
    check("static_compat", all(d_ != 0 for d_ in dets) and all(c_ != 0 for c_ in cv),
          f"at 3 exact points (rho, u, a, beta, s) with s not in {{0, 1, -1}}: det of the 2 x 2 matrix of (E1, E2) in (sigma_x, "
          f"sigma_rho) = {[float(d_) for d_ in dets]} != 0 and the compatibility C = {[float(c_) for c_ in cv]} != 0 "
          f"(exact rationals in the raw log). So the finite-type system is not Frobenius closed: a solution lies in C = 0 and the "
          f"general static class is NOT classified")
    info(f"points {[{str(k_): str(v_) for k_, v_ in p_.items()} for p_ in gen_pts]}; exact C = " + ", ".join(str(c_) for c_ in cv)
         + f"; exact det = " + ", ".join(str(d_) for d_ in dets) + f"; redraws {redraws}")
    if stop < 10:
        return fin(t0, mut)

    # 10. the linearization at b = C ln rho + b_0 ---------------------------------------------------------------------
    print("\n== 10. The linearization at the log solution u0 = C ln rho + b_0 ==")
    C0, b0s, eps, s0, s1 = sp.symbols("C0 b0s eps s0 s1")
    u0 = C0 * sp.log(rho) + b0s
    wf = sp.Function("wf")(x, rho)
    ue = u0 + eps * wf
    jl = static_blocks(ue, sp.diff(ue, x), sp.diff(ue, rho), sp.diff(ue, x, 2), sp.diff(ue, x, rho), sp.diff(ue, rho, 2), rho)
    zero_order = [sp.simplify(e_.subs(eps, 0)) for e_ in jl]
    ok0 = zero_order[:4] == [0, 0, 0, 0] and sp.simplify(zero_order[4] - C0**2 / (2 * rho**2)) == 0
    v1sq = s0 + eps * s1
    R1e = jl[2] * v1sq + jl[1] * jl[3]
    R2e = jl[0] * v1sq + jl[1]**2 + jl[3]**2
    d1 = sp.simplify(sp.diff(R1e, eps).subs(eps, 0) - sp.diff(wf, x, rho) / 2 * s0)
    sg = 1 if MUT == "lin-claim" else -1
    d2 = sp.simplify(sp.diff(R2e, eps).subs(eps, 0) - sg * u0 * sp.diff(wf, x, 2) * s0)
    gens1 = gensR
    r1m = rad_member(gens1, Zq * v1**2 + Qq * Yq, [v1], varsR)
    r2m = rad_member(gens1, Pq * v1**2 + Qq**2 + Yq**2, [v1], varsR)
    r3m = rad_member(gens1, v1**4 - 2 * hq * v1**2 + Qq**2, [v1], varsR)
    cc_, dd_ = sp.symbols("cc_ dd_")
    gfun = sp.Function("gfun")(rho)
    wlin = cc_ * x + gfun
    pl1 = sp.simplify(sp.diff(wlin, x, rho) / 2) == 0 and sp.simplify(u0 * sp.diff(wlin, x, 2)) == 0
    check("linearized", ok0 and d1 == 0 and d2 == 0 and r1m and r2m and r3m and pl1,
          f"at eps = 0: (P, Q, Z, Y) = 0 and h = C^2/(2 rho^2): {ok0}; Lemma R with v1 != 0 gives Z v1^2 + Q Y = 0, P v1^2 + Q^2 + Y^2 "
          f"= 0, v1^4 - 2 h v1^2 + Q^2 = 0 ({r1m}, {r2m}, {r3m}); their O(eps) parts are (w_x rho/2) v1^2 and "
          f"({sg}) u0 w_xx v1^2 with v1^2 = C^2/rho^2 != 0 ({d1 == 0}, {d2 == 0}), so the tangent space at u0 is w_x rho = 0 "
          f"and u0 w_xx = 0, i.e. w = c x + g(rho) (planted: this family satisfies both: {pl1}); no decaying w")
    return fin(t0, mut)


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

#!/usr/bin/env python3
"""Can two scalars of one sign source the Alcubierre wall? (issue #78)

Two minimally coupled scalars of one sign, (+,+) or (-,-), with any potential V, have
R_ab = lam g_ab + s (a_a a_b + b_a b_b), s = +1 or -1, lam = 8 pi V, a, b = sqrt(8 pi) d phi_I. The test has four stages,
as in two_scalars.py (#65):

    algebraic   some real lam makes M = R - lam g semidefinite of rank <= 2 (positive for s = +1, negative for s = -1)
    closure     some O(2) rotation Lambda(x) of (a, b) makes both closed. For a rank-2 form M with kernel plane field K
                this holds iff K is integrable, L_X M = 0 for X in K, and the quotient 2-metric of M by K is flat
    potential   K(lam) = 0 on the kernel, so that lam is a function of the fields
    residual    G_ab - 8 pi T_ab = 0 with the fields substituted back

With R_a3 = 0 and R_00 + R_11 + 2 R_22 = 0 ((I1), (I2) of phantom_source.py), the algebraic stage reduces to "lam != 0
and rank(R_3 - lam eta_3) = 1" (section 2). For the wall, the quantity F = R_01 (R_0rho^2 + R_1rho^2) + 2 R_rho,rho
R_0rho R_1rho must vanish, and it vanishes at finitely many angles at each radius (sections 4 and 6), so no open set of
the wall passes. Section 7 treats the shifts b(rho), where the stage passes on an open set and closure decides.

    python3 notes/warp/same_sign_scalars.py                 all checks, then every control (exit 1 on any failure)
    python3 notes/warp/same_sign_scalars.py --mutate=NAME   one control (its target checks must fail, exit 1)
    python3 notes/warp/same_sign_scalars.py --no-controls   the checks only
    python3 notes/warp/same_sign_scalars.py --stop=N        sections up to N

Imported read-only. From phantom_source.py: `ricci`, the frame Ricci tensors `sph_tetrad_ricci`, `lab_tetrad_ricci` and
`shift_tetrad_ricci`, the jet rules and the jet symbols. From two_scalars.py: `rq`, `vec4`, `lorentz`, `eta`, `lam`. From
three_scalars.py: `inertia`, `lie_K`, `wall_frame`, `tanh_jets_n`, `bump_jets_n` and the wall symbols. The algebraic test
`same_sign_form`, the rank-2 closure test `closure2`, the planted metrics, the F analysis and the controls are new here.
"""
import itertools
import os
import random
import subprocess
import sys
import time

import mpmath as mp
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phantom_source as ps  # noqa: E402
import three_scalars as t3  # noqa: E402
import two_scalars as ts  # noqa: E402

# name: (checks that must fail, what the mutation does)
MUTATIONS = {
    "planted-sign": (["planted_pass"], "builds the planted forms as lam eta + a a - b b"),
    "planted-fail-sign": (["planted_fail"], "builds the planted failures as two squares of one sign"),
    "lemmaS-weight": (["lemmaS_samples"], "builds the lemma's samples with R_00 - R_11 + 2 R_22 = 0 in place of (I2)"),
    "p1-mixed": (["p1_algebraic"], "gives P1's plane the signature (-,+): -dR^2 + R^2 dTheta^2"),
    "p2-mixed": (["p2_algebraic"], "gives P2's plane the signature (-,+)"),
    "p1-fiber": (["p1_closed"], "curves P1's plane, R^2 (1 + R^2/10)"),
    "p1-lambda": (["p1_potential"], "adds T/10 to lam in P1's potential stage"),
    "p2-lambda": (["p2_potential"], "adds T/10 to lam in P2's potential stage"),
    "p1-sign": (["p1_residual"], "flips one kinetic sign in T of P1"),
    "p2-sign": (["p2_residual"], "flips one kinetic sign in T of P2"),
    "d1-flat": (["d1_closure_fails"], "sets D1's plane curvature to 0"),
    "d2-gamma0": (["d2_algebraic", "d2_closure_fails"], "sets gamma = 0 in D2, so that M = 0"),
    "d3-bracket": (["d3_closure_fails"], "computes the kernel commutator without the term K_2(K_1^mu)"),
    "alc-drop-w2": (["wall_identities", "wall_entries", "wall_F"],
                    "drops the shift's square from g_tt in the spherical route (as in phantom_source.py)"),
    "wall-F-sign": (["wall_F"], "takes the claimed F formula with +z in place of -z"),
    "ideal-nosat": (["wall_ideal"], "drops the saturation 1 - zz w' r from the ideal of the s^2 coefficients"),
    "wall-planted": (["wall_points"], "replaces R at the off-axis exact points by a planted lam eta + a a + b b"),
    "wall-psd": (["wall_lambda0"], "replaces R at the off-axis exact points by a positive definite tensor"),
    "acc-static": (["acc_entries"], "drops every t-derivative of v_s, x_s and f from route D's jets"),
    "acc-drop-A": (["acc_entries"], "drops the A terms from the expected entries"),
    "acc-drop-odd": (["acc_ideal"], "keeps only the even powers of c in the accelerating ideal"),
    "brho-q": (["brho_rank1"], "takes the expected rank-1 solution with Q = h in place of Q = 0"),
}

# the last section a control needs: its target checks are all in sections up to this one
STOP_AT = {"planted-sign": 1, "planted-fail-sign": 1, "lemmaS-weight": 2, "p1-mixed": 3, "p2-mixed": 3, "p1-fiber": 3,
           "p1-lambda": 3, "p2-lambda": 3, "p1-sign": 3, "p2-sign": 3, "d1-flat": 3, "d2-gamma0": 3,
           "alc-drop-w2": 4, "wall-F-sign": 4, "ideal-nosat": 4, "wall-planted": 5, "wall-psd": 5,
           "acc-static": 6, "acc-drop-A": 6, "acc-drop-odd": 6, "brho-q": 7, "d3-bracket": 7}

RESULTS = []
MUT = None


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


lam = ts.lam
eta = ts.eta
rq, vec4, lorentz = ts.rq, ts.vec4, ts.lorentz
eta3 = sp.diag(-1, 1, 1)
lm = sp.Symbol("lm", real=True)


# ---------------------------------------------------------------------------------------------------------------
# The algebraic stage
# ---------------------------------------------------------------------------------------------------------------
def minors_gcd_roots(Rm, gm):
    """Distinct real roots of the gcd over lam of the 3 x 3 minors of Rm - lam gm (the lam where the rank is <= 2),
    sorted. Rm, gm have exact numeric entries."""
    A = Rm - lam * gm
    rows = list(itertools.combinations(range(A.shape[0]), 3))
    g = None
    for r_ in rows:
        for c_ in rows:
            num = sp.expand(sp.fraction(sp.together(A.extract(list(r_), list(c_)).det(method="berkowitz")))[0])
            if num == 0:
                continue
            P = sp.Poly(num, lam)
            g = P if g is None else sp.gcd(g, P)
            if g.degree() < 1:
                return []
    if g is None:
        return []
    rts = set()
    for fac, _ in sp.factor_list(g.as_expr())[1]:
        pf = sp.Poly(fac, lam)
        if pf.degree() >= 1:
            rts |= set(sp.real_roots(pf))
    return sorted(rts, key=lambda x: float(sp.N(x, 30)))


def same_sign_form(Rm, gm):
    """The algebraic stage at one point. Rm, gm are 4 x 4 with exact numeric entries. For every real root lam of the
    gcd of the 3 x 3 minors of R - lam g (rank <= 2) returns lam, the inertia (n+, n-) of M = R - lam g from 60-digit
    eigenvalues (numerical evidence when M has irrational entries), the rank and the admissible patterns: 'PP' when M is
    positive semidefinite (two normal scalars), 'MM' when it is negative semidefinite (two phantoms)."""
    out = []
    for lr in minors_gcd_roots(Rm, gm):
        M = (Rm - lr * gm) if lr.is_Rational else (Rm - lr * gm).applyfunc(lambda e: sp.N(e, 60))
        npos, nneg, _ = t3.inertia(M)
        pats = []
        if nneg == 0 and npos <= 2:
            pats.append("PP")
        if npos == 0 and nneg <= 2:
            pats.append("MM")
        out.append({"lam": lr, "sig": (npos, nneg), "rank": npos + nneg, "pats": pats})
    return out


def rank1_lams(R3):
    """The lemma's prediction from the 3 x 3 block alone: the real lam != 0 with rank(R3 - lam eta3) = 1, found as the
    real roots of the gcd over lam of the nine entries' 2 x 2 minors (exact). Returns [(lam, 'PP' or 'MM')]."""
    N = R3 - lam * eta3
    g = None
    for (i, j) in itertools.combinations(range(3), 2):
        for (k, l) in itertools.combinations(range(3), 2):
            m = sp.expand(N[i, k] * N[j, l] - N[i, l] * N[j, k])
            if m == 0:
                continue
            P = sp.Poly(m, lam)
            g = P if g is None else sp.gcd(g, P)
            if g.degree() < 1:
                return []
    if g is None:
        return []
    out = set()
    for fac, _ in sp.factor_list(g.as_expr())[1]:
        pf = sp.Poly(fac, lam)
        if pf.degree() >= 1:
            out |= {r_ for r_ in sp.real_roots(pf) if r_ != 0}
    return [(r_, "PP" if r_ < 0 else "MM") for r_ in sorted(out, key=lambda x: float(sp.N(x, 30)))]


# ---------------------------------------------------------------------------------------------------------------
# The closure stage for a rank-2 form
# ---------------------------------------------------------------------------------------------------------------
def closure2(M, Dfun, simp, lam_expr=None):
    """Closure and potential stages for a rank-2 symmetric 4 x 4 M with a two-dimensional kernel plane.

    The kernel basis is normalized to the identity on two coordinate rows (j1, j2), so that the basis vectors K_1, K_2
    commute iff the plane is integrable (their (j1, j2) components are constant). Returns a dict: 'dim' (the kernel
    dimension), 'rows' (j1, j2), 'K', 'integrable' ([K_1, K_2] = 0), 'invariant' (L_{K_i} M = 0, i = 1, 2), 'ric' (the
    Ricci tensor of the 2 x 2 block of M on the slice x^j1, x^j2 = const, the pullback of the quotient metric when M is
    invariant), 'flat', and with lam_expr 'klam' = [K_1(lam), K_2(lam)]."""
    n = 4
    ns = M.nullspace(iszerofunc=lambda e: simp(e) == 0)
    if len(ns) != 2:
        return {"dim": len(ns), "integrable": False, "invariant": False, "flat": False, "klam": [None, None]}
    Kb = sp.Matrix.hstack(*ns).applyfunc(simp)
    for (j1, j2) in itertools.combinations(range(n), 2):
        sub = Kb.extract([j1, j2], [0, 1])
        if simp(sub.det()) != 0:
            break
    Kb = (Kb * sub.inv()).applyfunc(simp)
    K1, K2 = Kb[:, 0], Kb[:, 1]
    br = []
    for m in range(n):
        t1 = sum(K1[a] * Dfun(K2[m], a) for a in range(n))
        t2 = 0 if MUT == "d3-bracket" else sum(K2[a] * Dfun(K1[m], a) for a in range(n))
        br.append(simp(t1 - t2))
    integrable = all(x == 0 for x in br)
    LK = [t3.lie_K(M, K, Dfun, simp) for K in (K1, K2)]
    invariant = all(simp(x) == 0 for L_ in LK for x in L_)
    idx = [i for i in range(n) if i not in (j1, j2)]
    h = M.extract(idx, idx)
    ric = ps.ricci(h, [0, 1], lambda E, a: Dfun(E, idx[a]), simp)
    out = {"dim": 2, "rows": (j1, j2), "K": (K1, K2), "bracket": br, "integrable": integrable, "invariant": invariant,
           "ric": ric, "flat": all(simp(x) == 0 for x in ric), "LK": LK}
    out["klam"] = [simp(sum(K[a] * Dfun(lam_expr, a) for a in range(n))) for K in (K1, K2)] \
        if lam_expr is not None else [None, None]
    return out


# ---------------------------------------------------------------------------------------------------------------
# Planted metrics: (a 2D Lorentzian space of constant curvature) x (a plane), the plane's Cartesian coordinates being
# the two scalars. P1 = AdS_2 x R^2 with two normal scalars, P2 = dS_2 x R^2 with two phantoms, in the polar
# coordinates (Ro, Tho) of the plane and in sheared coordinates (T, Z, Rn, Th) with
#   to = T + Z/3, zo = Z, Ro = Rn, Tho = Th + T/7 + Z/5,
# so that the kernel vectors K_1 = d_T - d_Th/7, K_2 = d_Z - d_Th/5 are not coordinate vectors and the rotation Lambda(x)
# between the eigenvector split of M and the Cartesian gradients depends on the point.
# ---------------------------------------------------------------------------------------------------------------
Tn, Zn, Rn, Thn = sp.symbols("T Z Rn Th", real=True)
XP = [Tn, Zn, Rn, Thn]
to_, zo_, Ro_, Tho_ = sp.symbols("to zo Ro Tho", real=True)
XO = [to_, zo_, Ro_, Tho_]
Rq = Rn
Thq = Thn + Tn / 7 + Zn / 5
OLD = {to_: Tn + Zn / 3, zo_: Zn, Ro_: Rq, Tho_: Thq}
JAC = sp.Matrix([[1, sp.Rational(1, 3), 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [sp.Rational(1, 7), sp.Rational(1, 5), 0, 1]])
simp_r = lambda e: sp.cancel(sp.together(e))
Dold = lambda E, k: sp.diff(E, XO[k])
Dnew = lambda E, k: sp.diff(E, XP[k])
Cq, Sq = sp.symbols("Cq Sq")


def geometry(go, lam_old):
    """g, Ricci and Einstein tensors and lam in the sheared coordinates, from the old-coordinate metric go."""
    Ric = ps.ricci(go, XO, Dold, simp_r)
    gi = go.inv().applyfunc(simp_r)
    Rs = simp_r(sum(gi[i, j] * Ric[i, j] for i in range(4) for j in range(4)))
    G = (Ric - go * Rs / 2).applyfunc(simp_r)
    tr = lambda m: (JAC.T * m.subs(OLD, simultaneous=True) * JAC).applyfunc(simp_r)
    return tr(go), tr(Ric), tr(G), simp_r(sp.sympify(lam_old).subs(OLD, simultaneous=True))


def red_trig(e):
    """Numerator of e with cos(Thq), sin(Thq) as symbols, reduced by Sq^2 = 1 - Cq^2: zero iff e vanishes."""
    e = sp.cancel(sp.expand(e).subs({sp.cos(Thq): Cq, sp.sin(Thq): Sq}))
    num, _ = sp.fraction(sp.together(e))
    return sp.expand(sp.rem(sp.Poly(sp.expand(num), Sq).as_expr(), Sq**2 + Cq**2 - 1, Sq))


def samp_new(kind):
    """A random exact point in the new coordinates, in the region where the metric and the semidefiniteness hold."""
    r_ = lambda lo, hi: sp.Rational(random.randint(lo, hi), random.randint(3, 7))
    if kind == "ads":  # zo = Z > 0
        return {Tn: r_(-3, 3), Zn: sp.Rational(random.randint(10, 30), 10), Rn: sp.Rational(random.randint(10, 30), 10),
                Thn: r_(-9, 9)}
    if kind == "ds":  # to = T + Z/3 > 0
        return {Tn: sp.Rational(random.randint(20, 40), 10), Zn: r_(-9, 9), Rn: sp.Rational(random.randint(10, 30), 10),
                Thn: r_(-9, 9)}
    # kind 'lor': to = T + Z/3 > |zo| = |Z|, so that M = (to^2 - zo^2)/25 (dX^2 + dY^2) is positive
    return {Tn: sp.Rational(random.randint(20, 30), 10), Zn: sp.Rational(random.randint(-10, 10), 10),
            Rn: sp.Rational(random.randint(10, 30), 10), Thn: r_(-9, 9)}


def lam_numeric_check(M, J, s, samp, tag):
    """Lambda = J A^T (A A^T)^-1 from a numerical factorization M = s A^T A (eigenvalues of s M, 60 digits) at three exact
    points: Lambda A = J, Lambda^T Lambda = 1, J^T J = s M, and Lambda differs between the points. Numerical evidence."""
    mp.mp.dps = 60
    lams, errs = [], []
    for _ in range(3):
        p = samp()
        Mp = mp.matrix([[mp.mpf(sp.N(M[i, j].subs(p), 60)) * s for j in range(4)] for i in range(4)])
        Jp = mp.matrix([[mp.mpf(sp.N(J[i, j].subs(p), 60)) for j in range(4)] for i in range(2)])
        ev, V = mp.eigsy(Mp)
        order = sorted(range(4), key=lambda i: -ev[i])[:2]
        A = mp.matrix(2, 4)
        for r_, i in enumerate(order):
            for j in range(4):
                A[r_, j] = mp.sqrt(ev[i]) * V[j, i]
        Lm = Jp * A.T * mp.inverse(A * A.T)
        mx = lambda B: max(abs(B[i, j]) for i in range(B.rows) for j in range(B.cols))
        errs.append(max(mx(Lm * A - Jp), mx(Lm.T * Lm - mp.eye(2)), mx(Jp.T * Jp - Mp)))
        lams.append(Lm)
    moves = max(abs(lams[0][i, j] - lams[1][i, j]) for i in range(2) for j in range(2))
    return max(errs), moves


def planted_stages(tag, geo, s, fields, kind, mut, deformed=None, sig_s=None):
    """The stages on a planted or deformed metric given by geometry(). s = +1 (normal scalars) or -1 (phantoms).
    fields are the Cartesian coordinates X^I of the plane in the old coordinates, a^I = dX^I, or None. deformed is None
    for a planted control (all four stages must pass), 'flat' for D1 (kernel integrable, L_X M = 0, quotient not flat),
    'lie' for D2 (L_X M != 0 and K(lam) != 0)."""
    g, Ric, G, lam_expr = geo
    t0 = time.time()
    want_sig = (2, 0) if s == 1 else (0, 2)
    res, lam_at = [], []
    for _ in range(3):  # the algebraic stage at exact points
        p = samp_new(kind)
        pr = lambda e, p=p: sp.Rational(sp.cancel(sp.sympify(e).subs(p)))
        res.append(same_sign_form(Ric.applyfunc(pr), g.applyfunc(pr)))
        lam_at.append(pr(lam_expr))
    ok_alg = all(any(c["lam"] == l_ and c["sig"] == want_sig for c in cands) for cands, l_ in zip(res, lam_at))
    check(f"{tag}_algebraic", ok_alg, f"lam = {[float(x) for x in lam_at]} is a real root of the rank-<=2 locus with "
          f"inertia {want_sig} at 3 exact points; found {[[(float(sp.N(c['lam'], 5)), c['sig']) for c in cs] for cs in res]}"
          f" ({time.time() - t0:.0f} s)")
    if not ok_alg and deformed is None:
        return None
    M = (Ric - lam_expr * g).applyfunc(simp_r)
    cl = closure2(M, Dnew, simp_r, lam_expr)
    if deformed is not None:
        if deformed == "flat":
            ok = ok_alg and cl["dim"] == 2 and cl["integrable"] and cl["invariant"] and not cl["flat"]
            why = (f"kernel dimension {cl['dim']}, integrable {cl['integrable']}, L_X M = 0: {cl['invariant']}, "
                   f"quotient flat: {cl['flat']}")
        else:
            ok = ok_alg and cl["dim"] == 2 and cl["integrable"] and not cl["invariant"] and \
                any(x != 0 for x in cl["klam"] if x is not None)
            why = (f"kernel dimension {cl['dim']}, integrable {cl['integrable']}, L_X M = 0: {cl['invariant']}, "
                   f"K(lam) = {cl['klam']}")
        check(f"{tag}_closure_fails", ok, why + f" ({time.time() - t0:.0f} s)")
        return cl
    ok_closed = ok_alg and cl["dim"] == 2 and cl["integrable"] and cl["invariant"] and cl["flat"]
    detail = (f"kernel dimension {cl['dim']}; [K_1, K_2] = 0: {cl['integrable']}; L_X M = 0 for both: {cl['invariant']}; "
              f"Ricci of the quotient block zero: {cl['flat']}")
    X1, X2 = fields
    J = sp.Matrix([[sp.diff(f.subs(OLD, simultaneous=True), x) for x in XP] for f in (X1, X2)])
    if ok_closed:
        err, moves = lam_numeric_check(M, J, s, lambda: samp_new(kind), tag)
        ok_lam = err < mp.mpf(10)**-45 and moves > mp.mpf(10)**-10
        ok_closed = ok_closed and ok_lam
        detail += (f"; Lambda = J A^T (A A^T)^-1 from a 60-digit factorization at 3 points: max error {mp.nstr(err, 3)} "
                   f"in Lambda A = J, Lambda^T Lambda = 1, J^T J = s M; Lambda differs between points by "
                   f"{mp.nstr(moves, 3)}")
    check(f"{tag}_closed", ok_closed, detail + f" ({time.time() - t0:.0f} s)")
    if not ok_closed:
        return cl
    pot = lam_expr + (Tn / 10 if mut == f"{tag}-lambda" else 0)
    klam = [simp_r(sum(K[a] * Dnew(pot, a) for a in range(4))) for K in cl["K"]]
    check(f"{tag}_potential", all(x == 0 for x in klam), f"K_1(lam) = {klam[0]}, K_2(lam) = {klam[1]}, lam = {lam_expr}")
    sl = [s, s]
    if mut == f"{tag}-sign":
        sl[0] = -sl[0]
    u = [sp.Matrix([sp.diff(f.subs(OLD, simultaneous=True), x) for x in XP]) for f in (X1, X2)]
    gi = g.inv().applyfunc(simp_r)
    kin = [simp_r((uu.T * gi * uu)[0]) for uu in u]
    Tab = sum((sl[i] * u[i] * u[i].T for i in range(2)), sp.zeros(4, 4)) - g * (sum(sl[i] * kin[i] for i in range(2)) / 2
                                                                                    + lam_expr)
    resid = G - Tab
    ok_res = all(red_trig(x) == 0 for x in resid)
    check(f"{tag}_residual", ok_res, f"G_ab - 8 pi T_ab with a^I = dX^I, signs {sl}, V = lam/(8 pi), substituted back: "
          f"{'0' if ok_res else [red_trig(x) for x in list(resid)[:3]]} ({time.time() - t0:.0f} s)")
    return cl


# ---------------------------------------------------------------------------------------------------------------
# The wall: symbols and the explicit entries of the frame block
# ---------------------------------------------------------------------------------------------------------------
RW, WJ, CC, SS = t3.RW, t3.WJ, t3.CC, t3.SS
W0_, W1_, W2_ = WJ[0], WJ[1], WJ[2]


def redsc(e):
    """e reduced modulo S^2 = 1 - c^2, as a polynomial in c and S (zero iff e vanishes, for polynomial e)."""
    num = sp.fraction(sp.together(e))[0]
    return sp.expand(sp.rem(sp.Poly(sp.expand(num), SS).as_expr(), SS**2 + CC**2 - 1, SS))


def uvz(r_, w0, w1, w2, S2):
    """u = r w'' - w', v = w u + 2 r w'^2, z = 2 w' + u S^2."""
    u = r_ * w2 - w1
    return u, w0 * u + 2 * r_ * w1**2, 2 * w1 + u * S2


def wall_entries_expected(mut):
    """The six entries of R_3 on (E_0, E_1, e_rho) at lam-independent level, in terms of u, v, z (derived by hand)."""
    u, v, z = uvz(RW, W0_, W1_, W2_, SS**2)
    R00 = -(2 * W0_ * W1_ + 2 * W0_ * u * CC**2 + RW * W1_**2 * (1 + CC**2)) / (2 * RW)
    R01 = -z / (2 * RW)
    R0p = u * SS * CC / (2 * RW)
    R1p = v * SS * CC / (2 * RW)
    Rpp = W1_**2 * SS**2 / 2
    R11 = -R00 - 2 * Rpp
    return sp.Matrix([[R00, R01, R0p], [R01, R11, R1p], [R0p, R1p, Rpp]]), (u, v, z)


def F_of(R3):
    return R3[0, 1] * (R3[0, 2]**2 + R3[1, 2]**2) + 2 * R3[2, 2] * R3[0, 2] * R3[1, 2]


def fin(t0, mut):
    failed = [n for n, ok in RESULTS if not ok]
    print(f"\nchecks: {len(RESULTS) - len(failed)} passed, {len(failed)} failed {failed}. {time.time() - t0:.0f} s")
    if mut is not None:
        want = MUTATIONS[mut][0]
        print(f"control {mut}: target checks {want}, fired {[n for n in want if n in failed]}")
    return failed


def main():
    global MUT
    opts = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
    mut = MUT = opts.get("mutate")
    if mut is not None and mut not in MUTATIONS:
        sys.exit(f"unknown mutation {mut!r}. Known: {', '.join(MUTATIONS)}")
    t0 = time.time()
    stop = int(opts.get("stop", STOP_AT.get(mut, 99)))
    print(f"same_sign_scalars.py, issue #78. sympy {sp.__version__}, mpmath {mp.__version__}, python "
          f"{sys.version.split()[0]}. mutation: {mut or 'none'}" + (f" ({MUTATIONS[mut][1]})" if mut else ""))

    # 1. Planted truth for the algebraic stage, in random Lorentz frames.
    print("\n== 1. Planted truth: R = lam eta + s (a a + b b) passes for s = +1 and -1, other forms fail ==")
    random.seed(781)
    bad = []
    for k in range(12):
        L_, lv, a, b = lorentz(), rq(), vec4(), vec4()
        s_ = 1 if k % 2 == 0 else -1
        sb = -s_ if mut == "planted-sign" else s_
        Rm = L_.T * (lv * eta + s_ * a * a.T + sb * b * b.T) * L_
        res = same_sign_form(Rm, eta)
        pat = "PP" if s_ == 1 else "MM"
        ok = any(c["lam"] == lv and pat in c["pats"] and c["rank"] == 2 for c in res)
        bad += [] if ok else [(lv, pat, [(c["lam"], c["sig"]) for c in res])]
    check("planted_pass", not bad, f"12 random lam eta + s (a a + b b), six with s = +1 and six with s = -1, transformed: "
          f"lam recovered with rank 2 and the planted pattern in all but {len(bad)} {bad[:2]}")
    bad = []
    kinds = [(1, -1, 0, 0), (-1, 1, 0, 0), (1, 1, 1, 0), (-1, -1, -1, 0), (1, 1, -1, 0), (1, -1, -1, 0)]
    for kind in kinds * 2:
        L_, lv = lorentz(), rq()
        sg = (1, 1, 0, 0) if mut == "planted-fail-sign" else kind
        vs = [vec4() for _ in range(4)]
        Rm = L_.T * (lv * eta + sum((sg[i] * vs[i] * vs[i].T for i in range(4)), sp.zeros(4, 4))) * L_
        res = same_sign_form(Rm, eta)
        bad += [] if not any(c["pats"] for c in res) else [(kind, lv, [(c["lam"], c["sig"]) for c in res])]
    check("planted_fail", not bad, f"12 random forms with squares (+,-), (-,+), (+,+,+), (-,-,-), (+,+,-), (+,-,-): "
          f"no admissible root in all but {len(bad)} {bad[:2]}")
    if stop < 2:
        return fin(t0, mut)

    # 2. The lemma: with (I1), (I2) the stage passes iff some real lam != 0 has rank(R_3 - lam eta_3) = 1.
    print("\n== 2. Lemma: with R_a3 = 0 and (I2), the form holds iff R_3 - lam eta_3 has rank 1 for a real lam != 0 ==")
    random.seed(782)
    wt = (1, -1, 2) if mut == "lemmaS-weight" else (1, 1, 2)
    wn = lambda u_: wt[0] * u_[0]**2 + wt[1] * u_[1]**2 + wt[2] * u_[2]**2
    mism, counts = [], [0, 0, 0]
    for trial in range(24):
        kind = trial % 3
        if kind == 0:  # R_3 = lam eta_3 + eps n n^T, lam fixed by the weighted trace
            n_ = sp.Matrix([rq() for _ in range(3)])
            if wn(n_) == 0:
                continue
            ep = 1 if (trial // 3) % 2 == 0 else -1
            R3 = (-ep * wn(n_) / 2) * eta3 + ep * n_ * n_.T
        elif kind == 1:  # a mixed rank-2 R_3 = a a^T - t b b^T (#65's det R_3 = 0 class)
            a_, b_ = sp.Matrix([rq() for _ in range(3)]), sp.Matrix([rq() for _ in range(3)])
            if wn(b_) == 0:
                continue
            R3 = a_ * a_.T - (wn(a_) / wn(b_)) * b_ * b_.T
        else:  # a generic R_3 with R_22 fixed by the weighted trace
            R3 = sp.zeros(3, 3)
            for i in range(3):
                for j in range(i, 3):
                    R3[i, j] = R3[j, i] = rq()
            R3[2, 2] = -(wt[0] * R3[0, 0] + wt[1] * R3[1, 1]) / wt[2]
        counts[kind] += 1
        pred = rank1_lams(R3)
        Rm = sp.zeros(4, 4)
        Rm[:3, :3] = R3
        L_ = lorentz()
        res = same_sign_form(L_.T * Rm * L_, eta)
        got = sorted([(c["lam"], tuple(c["pats"])) for c in res if c["pats"] and c["lam"] != 0],
                     key=lambda x: float(sp.N(x[0], 30)))
        want = [(l_, (p_,)) for l_, p_ in pred]
        same = len(got) == len(want) and all(abs(sp.N(g_[0] - w_[0], 40)) < 1e-30 and g_[1] == w_[1]
                                            for g_, w_ in zip(got, want))
        z0 = [c for c in res if c["lam"] == 0 and c["pats"]]
        if not same or z0:
            mism.append((trial, kind, want, got, bool(z0)))
    check("lemmaS_samples", not mism, f"24 random tensors with (I1), (I2) in random Lorentz frames ({counts[0]} of the form "
          f"lam eta_3 + eps n n^T, {counts[1]} mixed rank 2, {counts[2]} generic): the 4 x 4 test and the 3 x 3 rank-1 "
          f"criterion give the same roots and patterns (lam < 0 is (+,+), lam > 0 is (-,-)), none at lam = 0, in all but "
          f"{len(mism)} {mism[:2]}")
    if stop < 3:
        return fin(t0, mut)

    # 3. Planted and deformed metrics, through the same stages.
    print("\n== 3. Planted metrics: (2D Lorentzian space of constant curvature) x (plane), sheared polar coordinates ==")
    print("  P1: AdS_2 x R^2, ds^2 = (-dto^2 + dzo^2)/zo^2 + dRo^2 + Ro^2 dTho^2, two normal scalars, lam = -1")
    print("  P2: dS_2 x R^2, ds^2 = (-dto^2 + dzo^2)/to^2 + dRo^2 + Ro^2 dTho^2, two phantoms, lam = +1")
    random.seed(783)
    X1, X2 = Ro_ * sp.cos(Tho_), Ro_ * sp.sin(Tho_)
    def fiber(tag, eps0=0):
        """The plane's metric components (g_RR, g_ThTh): signature and curvature of the plane, for the controls."""
        sg = -1 if mut == f"{tag}-mixed" else 1
        eps = sp.Rational(1, 10) if mut == f"{tag}-fiber" else eps0
        return sg, Ro_**2 * (1 + eps * Ro_**2)
    ads = lambda sg, f2: sp.diag(-1 / zo_**2, 1 / zo_**2, sg, f2)
    dsp = lambda sg, f2: sp.diag(-1 / to_**2, 1 / to_**2, sg, f2)
    geo1 = geometry(ads(*fiber("p1")), -1)
    planted_stages("p1", geo1, 1, (X1, X2), "ads", mut)
    geo2 = geometry(dsp(*fiber("p2")), 1)
    planted_stages("p2", geo2, -1, (X1, X2), "ds", mut)
    print(f"  ({time.time() - t0:.0f} s)")
    print("\n== 3b. Deformed metrics: D1 curves the plane, D2 puts the Lorentzian factor off shell ==")
    print("  D1: AdS_2 x (plane with Ro^2 -> Ro^2 (1 + Ro^2/10)): kernel integrable, L_X M = 0, quotient not flat")
    print("  D2: ds^2 = -dto^2 + dzo^2 + E^2 (dX^2 + dY^2), E = 1 + gamma (zo^2 - to^2), lam = -4 gamma/E: "
          "R - lam g = gamma^2 4 (to^2 - zo^2) (dX^2 + dY^2)")
    geoD1 = geometry(ads(1, Ro_**2 * (1 + (0 if mut == "d1-flat" else sp.Rational(1, 10)) * Ro_**2)), -1)
    planted_stages("d1", geoD1, 1, None, "ads", mut, deformed="flat")
    gam = 0 if mut == "d2-gamma0" else sp.Rational(1, 10)
    E_ = 1 + gam * (zo_**2 - to_**2)
    geoD2 = geometry(sp.diag(-1, 1, E_**2, E_**2), -4 * gam / E_)
    planted_stages("d2", geoD2, 1, None, "lor", mut, deformed="lie")
    print(f"  ({time.time() - t0:.0f} s)")
    if stop < 4:
        return fin(t0, mut)

    # 4. The wall, constant v, general profile: the entries, F and its coefficients.
    print("\n== 4. Alcubierre wall, constant v: R_3 on (E_0 = n, E_1 = d_x, e_rho), general profile w = v (1 - f) ==")
    Rs = t3.wall_frame(mut)
    R3 = Rs[:3, :3]
    ids = [redsc(sp.simplify(Rs[i, 3])) for i in range(4)] + [redsc(sp.simplify(R3[0, 0] + R3[1, 1] + 2 * R3[2, 2]))]
    check("wall_identities", all(x == 0 for x in ids), f"R(e_phi, .) and R_00 + R_11 + 2 R_rho,rho: {ids} "
          f"({time.time() - t0:.0f} s)")
    Rexp, (u, v, z) = wall_entries_expected(mut)
    diffs = [redsc(sp.simplify(R3[i, j] - Rexp[i, j])) for i in range(3) for j in range(i, 3)]
    check("wall_entries", all(x == 0 for x in diffs), f"with u = r w'' - w', v = w u + 2 r w'^2, z = 2 w' + u s^2: "
          f"R_00 = -(2 w w' + 2 w u c^2 + r w'^2 (1 + c^2))/(2r), 2r R_01 = -z, 2r R_0rho = u s c, 2r R_1rho = v s c, "
          f"R_rho,rho = w'^2 s^2/2: differences {diffs}")
    # F, the claimed factorization and the certificate R_0rho m_d + R_1rho m_e + F = 0 for the generic N(lam)
    Fv = F_of(R3)
    sgn = 1 if mut == "wall-F-sign" else -1
    Bs = lambda sig: sgn * z.subs(SS**2, sig) * (u**2 + v**2) + 2 * RW * u * v * W1_**2 * sig
    Fexp = SS**2 * CC**2 * Bs(SS**2) / (8 * RW**3)
    dF = redsc(sp.simplify(Fv - Fexp))
    N = R3 - lam * eta3
    md = N[0, 0] * N[1, 2] - N[0, 1] * N[0, 2]
    me = N[1, 1] * N[0, 2] - N[0, 1] * N[1, 2]
    cert = redsc(sp.simplify(R3[0, 2] * md + R3[1, 2] * me + Fv))
    check("wall_F", dF == 0 and cert == 0, f"F = R_01 (R_0rho^2 + R_1rho^2) + 2 R_rho,rho R_0rho R_1rho equals "
          f"s^2 c^2 ({'-' if sgn < 0 else '+'}z (u^2 + v^2) + 2 r u v w'^2 s^2)/(8 r^3): difference {dF}; the certificate "
          f"R_0rho (N_00 N_1rho - N_01 N_0rho) + R_1rho (N_11 N_0rho - N_01 N_1rho) + F = {cert}")
    # the coefficients of the bracket in sigma = s^2: -2 w' (u^2 + v^2) and -u (u^2 + v^2) + 2 r u v w'^2
    sig = sp.Symbol("sigma")
    Bpoly = sp.Poly(sp.expand(-(2 * W1_ + u * sig) * (u**2 + v**2) + 2 * RW * u * v * W1_**2 * sig), sig)
    coefs = [sp.factor(Bpoly.coeff_monomial(sig**k)) for k in range(2)]
    zz = sp.Symbol("zz")
    gens = (zz, W0_, W2_, W1_, RW)
    sat = [] if mut == "ideal-nosat" else [1 - zz * W1_ * RW]
    Gb = sp.groebner([c_ for c_ in coefs] + sat, *gens, order="grevlex")
    check("wall_ideal", list(Gb.exprs) == [1], f"coefficients of the bracket in s^2: {coefs[0]} and {coefs[1]}; with "
          f"w' r != 0 they generate the unit ideal (Groebner basis {list(Gb.exprs)[:3]}), so at each r with w' != 0 "
          f"the bracket, linear in s^2 with the constant term 2 w'(u^2 + v^2) != 0, vanishes at one s^2 at most "
          f"({time.time() - t0:.0f} s)")
    print(f"  ({time.time() - t0:.0f} s)")
    if stop < 5:
        return fin(t0, mut)

    # 5. Exact points of the wall: the lab frame, the same same_sign_form as the planted controls.
    print("\n== 5. Exact points (lab frame route of phantom_source.py): the algebraic stage through same_sign_form ==")
    Rt = ps.lab_tetrad_ricci(None)
    random.seed(57)
    samples = []
    wall = [(12, 3, 4, 13), (-12, 4, 3, 13), (13, 0, 0, 13), (-13, 0, 0, 13), (0, 5, 12, 13), (5, 12, 0, 13),
            (-5, 0, 12, 13), (9, 12, 0, 15)]
    for (x0, y0, z0, r0) in wall:
        pnt = [sp.Rational(q_, r0) for q_ in (x0, y0, z0, r0)]
        for vv in (sp.Integer(1), sp.Rational(1, 2), sp.Integer(3)):
            samples.append(("tanh sig=8 Rb=1", pnt, vv, t3.tanh_jets_n(pnt[3], 8, 1, 7)))
    for (x0, y0, z0, r0) in [(4, 3, 0, 5), (-3, 0, 4, 5), (5, 0, 0, 5), (0, 3, 4, 5)]:
        pnt = [sp.Rational(q_, 5) for q_ in (x0, y0, z0, r0)]
        samples.append(("bump Rb=2", pnt, sp.Rational(7, 4), t3.bump_jets_n(pnt[3], 2, 7)))
    for (x0, y0, z0, r0) in [(12, 3, 4, 13), (-12, 4, 3, 13), (4, 3, 0, 5), (-3, 0, 4, 5)]:
        pnt = [sp.Rational(q_, 7) for q_ in (x0, y0, z0, r0)]
        samples.append(("random jets", pnt, sp.Rational(random.randint(1, 9), 3),
                        [sp.Rational(random.randint(-9, 9), random.randint(1, 9)) for _ in range(7)]))
    random.seed(785)
    pts, tens = [], []
    for (lab_, pnt, vv, jets) in samples:
        sub = {ps.X: pnt[0], ps.y: pnt[1], ps.z: pnt[2], ps.R: pnt[3], ps.v: vv, ps.F[0]: jets[0], ps.F[1]: jets[1],
               ps.F[2]: jets[2]}
        Rm = Rt.subs(sub)
        axis = pnt[1] == 0 and pnt[2] == 0
        if not axis and mut == "wall-planted":
            L_, a, b = lorentz(), vec4(), vec4()
            Rm = L_.T * (rq() * eta + a * a.T + b * b.T) * L_
        if not axis and mut == "wall-psd":
            A_ = sp.Matrix(4, 4, lambda i, j: rq())
            Rm = A_.T * A_ + sp.eye(4)
        pts.append((lab_, axis, pnt, vv))
        tens.append(Rm)
    ok_pts, n_off, n_ax, n_roots_off, n_adm = True, 0, 0, 0, 0
    rows = []
    for (lab_, axis, pnt, vv), Rm in zip(pts, tens):
        res = same_sign_form(Rm, eta)
        adm = [c for c in res if c["pats"]]
        ok_pts &= not adm
        n_adm += len(adm)
        n_ax += int(axis)
        n_off += int(not axis)
        n_roots_off += 0 if axis else len(res)
        rows.append((lab_, axis, [(float(sp.N(c["lam"], 5)), c["sig"]) for c in res]))
    for key in sorted(set((r_[0], r_[1]) for r_ in rows)):
        sub_ = [r_ for r_ in rows if (r_[0], r_[1]) == key]
        print(f"  {key[0]:>16}, {'on the axis' if key[1] else 'off the axis'}: {len(sub_)} points, rank-<=2 roots (lam, "
              f"inertia) at the points: {sorted(set(str(r_[2]) for r_ in sub_))}")
    check("wall_points", ok_pts and n_off > 0, f"{len(pts)} exact points ({n_off} off the axis, {n_ax} on it): the "
          f"same function as the planted forms finds no real lam with R - lam g semidefinite of rank <= 2, at lam = 0 or "
          f"lam != 0; real roots of the rank-<=2 locus at off-axis points: {n_roots_off}; admissible roots at any "
          f"point: {n_adm}")
    bad0 = []
    for (lab_, axis, pnt, vv), Rm in zip(pts, tens):
        if axis:
            continue
        npos, nneg, _ = t3.inertia(Rm)
        if npos < 1 or nneg < 1:
            bad0.append((lab_, pnt, (npos, nneg)))
    check("wall_lambda0", not bad0, f"at lam = 0, R has both signs (inertia with n+ >= 1 and n- >= 1, 60 digits) at "
          f"{n_off - len(bad0)} of {n_off} off-axis points, so M = R is not semidefinite; failures {bad0[:2]}")
    print(f"  ({time.time() - t0:.0f} s)")
    if stop < 6:
        return fin(t0, mut)

    # 6. The accelerating bubble of phantom_source.py section 6 (route D), at v_s(t) != 0.
    print("\n== 6. Accelerating bubble b = -v_s(t) f(r_s): the frame block with A = v_s'/v_s ==")
    Ra = ps.shift_tetrad_ricci(-ps.VS[0] * ps.F[0], ps.bubble_rule(mut == "acc-static"), False)
    A, U, Vs, Rr = sp.Symbol("A", real=True), sp.Symbol("U", real=True), sp.Symbol("V", real=True), ps.R
    cs, sn = sp.Symbol("c", real=True), sp.Symbol("S", real=True)
    jsub = {ps.F[0]: (U - W0_) / Vs, ps.F[1]: -W1_ / Vs, ps.F[2]: -W2_ / Vs, ps.VS[1]: A * Vs, ps.XS[1]: U, ps.VS[0]: Vs}
    R3a = sp.zeros(3, 3)
    for i in range(3):
        for j in range(i, 3):
            e = Ra[i, j].subs(ps.z, 0).subs(jsub).subs({ps.X: Rr * cs, ps.y: Rr * sn})
            R3a[i, j] = R3a[j, i] = sp.factor(sp.expand(e))
    free = all(not (e.has(U) or e.has(Vs)) for e in R3a)
    redsa = lambda e: sp.expand(sp.rem(sp.Poly(sp.expand(sp.fraction(sp.together(e))[0]), sn).as_expr(),
                                       sn**2 + cs**2 - 1, sn))
    ua = Rr * W2_ - W1_
    va = W0_ * ua + 2 * Rr * W1_**2
    za = 2 * W1_ + ua * sn**2
    Aexp = 0 if mut == "acc-drop-A" else A
    R00a = -(2 * W0_ * W1_ + 2 * W0_ * ua * cs**2 + Rr * W1_**2 * (1 + cs**2)) / (2 * Rr) + Aexp * W1_ * cs
    Rppa = W1_**2 * sn**2 / 2
    Rexa = sp.Matrix([[R00a, -za / (2 * Rr), ua * sn * cs / (2 * Rr)],
                      [-za / (2 * Rr), -R00a - 2 * Rppa, sn * (va * cs - Aexp * Rr * W1_) / (2 * Rr)],
                      [ua * sn * cs / (2 * Rr), sn * (va * cs - Aexp * Rr * W1_) / (2 * Rr), Rppa]])
    dacc = [redsa(sp.simplify(R3a[i, j] - Rexa[i, j])) for i in range(3) for j in range(i, 3)] if free else None
    check("acc_jets", free, f"route D's frame block depends on x_s', v_s, v_s', f, f', f'' only through w~, w~', w~'' and "
          f"A = v_s'/v_s: {free}")
    check("acc_entries", free and all(x == 0 for x in dacc), f"the entries equal the constant-speed ones with R_00 += A w' c, "
          f"R_11 -= A w' c and R_1rho = s (v c - A r w')/(2 r): differences {dacc}")
    # lab-frame 4 x 4 tensor of route D at a constant speed equals route C's lab-frame tensor
    ok_lim, nlim = True, 0
    random.seed(786)
    for (Xv, yv, Rv) in ((sp.Rational(12, 13), sp.Rational(5, 13), sp.Integer(1)),
                         (sp.Rational(-3, 5), sp.Rational(4, 5), sp.Integer(1)),
                         (sp.Rational(4, 5), sp.Rational(3, 5), sp.Integer(1))):
        vv = sp.Rational(random.randint(1, 9), random.randint(2, 5))
        fj = [sp.Rational(random.randint(-9, 9), random.randint(1, 9)) for _ in range(5)]
        pa = {ps.X: Xv, ps.y: yv, ps.z: 0, ps.R: Rv}
        pa.update({ps.F[k]: fj[k] for k in range(5)})
        pa.update({ps.VS[0]: vv, ps.VS[1]: 0, ps.VS[2]: 0, ps.VS[3]: 0, ps.XS[1]: vv, ps.XS[2]: 0, ps.XS[3]: 0})
        pc = {ps.X: Xv, ps.y: yv, ps.z: 0, ps.R: Rv, ps.v: vv, ps.F[0]: fj[0], ps.F[1]: fj[1], ps.F[2]: fj[2]}
        ok_lim &= (Ra.subs(pa) - Rt.subs(pc)).applyfunc(sp.simplify) == sp.zeros(4, 4)
        nlim += 1
    check("acc_const_limit", ok_lim, f"at v_s' = 0, x_s' = v_s the 4 x 4 frame Ricci tensor of route D equals route C's "
          f"lab-frame tensor at {nlim} exact points")
    # F for the bubble and its coefficients in c
    Fa = F_of(R3a) if free else None
    if free:
        B0 = A * Rr * W1_
        Bc = -za * (ua**2 * cs**2 + (va * cs - B0)**2) + 2 * Rr * W1_**2 * ua * sn**2 * cs * (va * cs - B0)
        dFa = redsa(sp.simplify(Fa - sn**2 * Bc / (8 * Rr**3)))
        Bp = sp.Poly(sp.expand(sp.expand(Bc).subs(sn**2, 1 - cs**2)), cs)
        cf = [sp.factor(Bp.coeff_monomial(cs**k)) for k in range(Bp.degree() + 1)]
        c0 = sp.simplify(cf[0] + (2 * W1_ + ua) * B0**2)
        c1 = sp.simplify(cf[1] - 2 * B0 * (va * (2 * W1_ + ua) - Rr * W1_**2 * ua))
        check("acc_F", dFa == 0 and c0 == 0 and c1 == 0, f"F = s^2 B(c)/(8 r^3) with B = -z (u^2 c^2 + (v c - A r w')^2) + "
              f"2 r w'^2 u s^2 c (v c - A r w'): difference {dFa}; c^0 coefficient {cf[0]} = -(2 w' + u)(A r w')^2: "
              f"{c0 == 0}, c^1 coefficient = 2 A r w' (v (2 w' + u) - r w'^2 u): {c1 == 0}")
        used = cf[0::2] if mut == "acc-drop-odd" else cf
        strip = lambda p: sp.Mul(*[f**m for f, m in sp.factor_list(p)[1] if f not in (A, W1_, Rr)])
        stripped = [strip(x) for x in used if x != 0]
        Ga = sp.groebner(stripped + [1 - zz * A * W1_ * Rr], zz, W0_, W2_, W1_, A, Rr, order="grevlex")
        check("acc_ideal", list(Ga.exprs) == [1], f"with A w' r != 0 the {len(used)} coefficients of B generate the unit "
              f"ideal (Groebner basis {list(Ga.exprs)[:3]}), so B vanishes on no interval of c at fixed r_s and t "
              f"({time.time() - t0:.0f} s)")
    # instants with v_s = 0: w' = 0 and R_3 = [[-k c, 0, 0], [0, k c, k s/2], [0, k s/2, 0]], k = v_s' f'
    kk = sp.Symbol("k", real=True)
    R3v = Ra[:3, :3].subs(ps.z, 0).subs({ps.VS[0]: 0}).subs({ps.X: Rr * cs, ps.y: Rr * sn})
    kexp = sp.Matrix([[-kk * cs, 0, 0], [0, kk * cs, kk * sn / 2], [0, kk * sn / 2, 0]])
    dv0 = [redsa(sp.simplify(x)) for x in (R3v - kexp.subs(kk, ps.VS[1] * ps.F[1]))]
    Nv = kexp - lm * eta3
    mins = [sp.expand(Nv[i, k] * Nv[j, l] - Nv[i, l] * Nv[j, k]) for (i, j) in itertools.combinations(range(3), 2)
            for (k, l) in itertools.combinations(range(3), 2)]
    Gv = sp.groebner(mins + [sn**2 + cs**2 - 1, 1 - zz * kk * sn], zz, lm, kk, cs, sn, order="grevlex")
    check("acc_vs0", all(x == 0 for x in dv0) and list(Gv.exprs) == [1], f"at v_s = 0 the block of route D is "
          f"[[-k c, 0, 0], [0, k c, k s/2], [0, k s/2, 0]] with k = v_s' f' (differences {dv0}), and the 2 x 2 minors of "
          f"R_3 - lam eta_3 with k s != 0 generate the unit ideal (Groebner basis {list(Gv.exprs)[:3]}), so N has rank "
          f">= 2 off the axis; on the axis (s = 0) lam = k c gives rank 1")
    random.seed(787)
    nacc, bad_acc = 0, []
    for (Xv, yv, Rv) in ((sp.Rational(12, 13), sp.Rational(5, 13), sp.Integer(1)),
                         (sp.Rational(-3, 5), sp.Rational(4, 5), sp.Integer(1)),
                         (sp.Rational(4, 5), sp.Rational(3, 5), sp.Integer(1)),
                         (sp.Rational(5, 13), sp.Rational(12, 13), sp.Integer(1)),
                         (sp.Rational(-5, 13), sp.Rational(12, 13), sp.Integer(1)),
                         (sp.Integer(0), sp.Integer(1), sp.Integer(1))):
        rn = lambda: sp.Rational(random.randint(-9, 9), random.randint(1, 9))
        vs = [sp.Rational(random.randint(1, 9), random.randint(2, 5)), rn(), rn(), rn()]
        xs = [0, rn(), rn(), rn()]
        fj = [rn() for _ in range(5)]
        pt = {ps.X: Xv, ps.y: yv, ps.z: 0, ps.R: Rv}
        pt.update({ps.F[k]: fj[k] for k in range(5)})
        pt.update({ps.VS[k]: vs[k] for k in range(4)})
        pt.update({ps.XS[k]: xs[k] for k in range(1, 4)})
        res = same_sign_form(Ra.subs(pt), eta)
        nacc += 1
        if any(c["pats"] for c in res):
            bad_acc.append((Xv, yv, [(c["lam"], c["sig"]) for c in res]))
    check("acc_points", not bad_acc, f"random jets of v_s (v_s' != 0), x_s and f at {nacc} exact points of the half-plane: "
          f"no admissible root at all but {len(bad_acc)} {bad_acc[:2]} ({time.time() - t0:.0f} s)")
    if stop < 7:
        return fin(t0, mut)

    # 7. A shift b(rho): the same-sign stage passes on an open set only for b = C ln rho + b_0, and closure fails.
    print("\n== 7. A shift b(rho) alone: R_3 = [[-h, Q, 0], [Q, -h, 0], [0, 0, h]], h = b'^2/2, Q = -(b'' + b'/rho)/2 ==")
    Rg = ps.shift_tetrad_ricci(ps.bjet(0, 0, 0), ps.general_rule(False), False)
    kill = {s_: 0 for s_ in list(ps.BJ) if isinstance(s_, sp.Symbol) and (ps.BJ[s_][0] or ps.BJ[s_][1])}
    rho = ps.RHO
    R3g = Rg[:3, :3].subs(kill).subs({ps.z: 0, ps.y: rho})
    b1, b2 = ps.bjet(0, 0, 1), ps.bjet(0, 0, 2)
    hh, QQ = sp.symbols("h Q", real=True)
    h_e, Q_e = b1**2 / 2, -(b2 + b1 / rho) / 2
    R3e = sp.Matrix([[-h_e, Q_e, 0], [Q_e, -h_e, 0], [0, 0, h_e]])
    dR3 = [sp.simplify(R3g[i, j] - R3e[i, j]) for i in range(3) for j in range(i, 3)]
    # rank-1 locus of N = R_3 - lam eta_3 = [[lam - h, Q, 0], [Q, -lam - h, 0], [0, 0, h - lam]]
    Nb = sp.Matrix([[lm - hh, QQ, 0], [QQ, -lm - hh, 0], [0, 0, hh - lm]])
    minb = [sp.expand(Nb[i, k] * Nb[j, l] - Nb[i, l] * Nb[j, k]) for (i, j) in itertools.combinations(range(3), 2)
            for (k, l) in itertools.combinations(range(3), 2)]
    Gb1 = sp.groebner(minb + [1 - zz * hh], zz, lm, QQ, hh, order="grevlex")
    want = [lm - hh, QQ - hh if mut == "brho-q" else QQ]
    Gw = sp.groebner(want + [1 - zz * hh], zz, lm, QQ, hh, order="grevlex")
    subs_ = [Gw.contains(p_) for p_ in minb]
    sups_ = [Gb1.contains(p_) for p_ in want]
    bsol = sp.Symbol("C0") * sp.log(rho) + sp.Symbol("b0")
    qsol = sp.simplify(-(sp.diff(bsol, rho, 2) + sp.diff(bsol, rho) / rho) / 2)
    check("brho_rank1", all(x == 0 for x in dR3) and all(subs_) and all(sups_) and qsol == 0,
          f"R_3 of a shift b(rho) equals [[-h, Q, 0], [Q, -h, 0], [0, 0, h]] (differences {dR3}); with h != 0 the 2 x 2 "
          f"minors of R_3 - lam eta_3 generate the same ideal as lam - h and {want[1]} (inclusions {all(subs_)}, "
          f"{all(sups_)}); Q = 0 for b = C ln rho + b_0: {qsol == 0}. So lam = h > 0 and Q = 0, pattern (-,-) "
          f"({time.time() - t0:.0f} s)")
    print("\n== 7b. D3: b = C ln rho + b_0, ds^2 = -dt^2 + (dx + b dt)^2 + drho^2 + rho^2 dphi^2, passes the algebraic "
          "stage on an open set ==")
    tt, xx, ph = sp.symbols("t x varphi", real=True)
    Lg = sp.Symbol("Lg", real=True)  # Lg = ln(rho) + const, d Lg/d rho = 1/rho
    Cc = sp.Integer(2)
    bb = Cc * Lg
    gc = sp.Matrix([[-1 + bb**2, bb, 0, 0], [bb, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, rho**2]])
    Xc = [tt, xx, rho, ph]

    def Dc(E_, k):
        return sp.diff(E_, rho) + sp.diff(E_, Lg) / rho if k == 2 else sp.diff(E_, Xc[k])
    Ricc = ps.ricci(gc, Xc, Dc, simp_r)
    lam_c = Cc**2 / (2 * rho**2)
    random.seed(788)
    res3, lam3 = [], []
    for _ in range(3):
        p = {rho: sp.Rational(random.randint(5, 30), 10), Lg: rq(), tt: rq(), xx: rq()}
        pr = lambda e, p=p: sp.Rational(sp.cancel(sp.sympify(e).subs(p)))
        res3.append(same_sign_form(Ricc.applyfunc(pr), gc.applyfunc(pr)))
        lam3.append(pr(lam_c))
    ok3 = all(any(c["lam"] == l_ and c["sig"] == (0, 2) and "MM" in c["pats"] for c in cands) for cands, l_ in zip(res3, lam3))
    check("d3_algebraic", ok3, f"lam = h = C^2/(2 rho^2) = {[float(x) for x in lam3]} is a root with inertia (0,2), two "
          f"phantoms, at 3 exact points; found {[[(float(sp.N(c['lam'], 5)), c['sig']) for c in cs] for cs in res3]}")
    M3 = (Ricc - lam_c * gc).applyfunc(simp_r)
    cl3 = closure2(M3, Dc, simp_r, lam_c)
    check("d3_closure_fails", ok3 and cl3["dim"] == 2 and not cl3["integrable"],
          f"kernel dimension {cl3['dim']}, [K_1, K_2] = {[str(x) for x in cl3.get('bracket', [])]}: the plane is not "
          f"integrable (Frobenius), so no d phi_I with the kernel as common null space exists ({time.time() - t0:.0f} s)")
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

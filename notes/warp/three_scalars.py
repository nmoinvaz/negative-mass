#!/usr/bin/env python3
"""Can three scalars of mixed signs source the Alcubierre wall? (issue #77)

Three minimally coupled scalars with kinetic signs s_I = +-1, (+,+,-) or (+,-,-), and any potential V(phi), have
R_ab = lam g_ab + sum_I s_I a^I_a a^I_b with lam = 8 pi V and a^I = sqrt(8 pi) d phi^I. The test has four stages, as in
two_scalars.py (#65).

    algebraic   some real lam makes M = R - lam g of rank <= 3 with at most 2 positive and 1 negative squares (or the
                mirror), that is, lam is a real root of det(R - lam g) and the inertia of M fits
    closure     some field-space map Lambda(x) in O(2,1) (or O(1,2)) makes all a^I closed. For rank M = 3, with kernel
                K, this holds iff L_K M = 0 and the 3 x 3 block of M on a slice transverse to K is flat: with
                da^I = -Omega^I_J a^J, Omega is the Levi-Civita connection of that block and Lambda exists iff its
                curvature vanishes (the 3D Ricci tensor is zero)
    potential   K(lam) = 0, that is d lam in span(a^I), so that V is a function of the phi^I
    residual    G_ab - 8 pi T_ab = 0 with the fields substituted back

Section 1 plants forms with known answers. Section 2 is the lemma behind the wall's algebraic stage. Section 3 runs the
test on planted three-scalar metrics (a cone and a cosh-warped product over a flat fiber, in sheared polar
coordinates) and on deformed ones. Sections 4 to 6 run it on the wall.

    python3 notes/warp/three_scalars.py                 all checks, then every control (exit 1 if any check fails)
    python3 notes/warp/three_scalars.py --mutate=NAME   one control (its target checks must fail, exit 1)
    python3 notes/warp/three_scalars.py --no-controls   the checks only

Imported read-only. From two_scalars.py: `rq`, `vec4`, `lorentz`, `eta`, `sign_of`, `ssqrt`, `lam`. From
phantom_source.py: `ricci`, the comoving spherical frame Ricci tensor `sph_tetrad_ricci`, the lab frame Ricci tensor
`lab_tetrad_ricci`, `bump_jets` and the jet symbols. The three-scalar test, the closure test and the controls are new
here.
"""
import os
import random
import subprocess
import sys
import time

import mpmath as mp
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phantom_source as ps  # noqa: E402
import two_scalars as ts  # noqa: E402

# name: (checks that must fail, what the mutation does)
MUTATIONS = {
    "planted3-sign": (["planted3_pass"], "builds the planted forms as lam eta + a a + b b + c c"),
    "planted3-fail-sign": (["planted3_fail"], "builds the planted failures with an admissible pattern"),
    "lemma3-weight": (["lemma3_samples"], "builds the lemma's samples with R_00 - R_11 + 2 R_22 = 0 in place of (I2)"),
    "p1-sign": (["p1_residual"], "flips one scalar's kinetic sign in T of the cone"),
    "p2-lambda": (["p2_potential"], "adds L/10 to lam in the potential stage of the cosh control"),
    "p2-fiber": (["p2_closed"], "curves the fiber of the cosh control"),
    "d1-flat": (["d1_closure_fails"], "sets the curvature of the deformed fiber to 0"),
    "d2-onshell": (["d2_closure_fails"], "puts the deformed warp factor on shell"),
    "alc-drop-w2": (["wall_identities", "wall_D", "wall_points"],
                    "drops the shift's square from g_tt in the wall (as in phantom_source.py)"),
    "wall-flat-sub": (["wall_flat_points"], "replaces the wall's orbit-space block by a flat one"),
    "lam-const": (["wall_lam_potential"], "takes d lam = 0 in the wall's lam != 0 potential test"),
}

RESULTS = []
MUT = None


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


lam = ts.lam
eta = ts.eta
rq, vec4, lorentz = ts.rq, ts.vec4, ts.lorentz


# ---------------------------------------------------------------------------------------------------------------
# The algebraic stage
# ---------------------------------------------------------------------------------------------------------------
def inertia(M, dps=60):
    """(n+, n-, n0) of a real symmetric matrix with exact entries, from its eigenvalues at dps digits. An eigenvalue
    counts as zero below 10^(-dps/2) times the largest one (or 1). Numerical evidence when M has irrational entries."""
    mp.mp.dps = dps
    n = M.shape[0]
    A = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            A[i, j] = mp.mpf(0) if M[i, j] == 0 else mp.mpf(sp.N(M[i, j], dps))
    ev = mp.eigsy(A, eigvals_only=True)
    tol = mp.mpf(10) ** (-(dps // 2)) * max([abs(x) for x in ev] + [mp.mpf(1)])
    return (sum(1 for x in ev if x > tol), sum(1 for x in ev if x < -tol), sum(1 for x in ev if abs(x) <= tol))


def real_lams(Rm, gm):
    """Distinct real roots of det(Rm - lam gm), sorted. None if the determinant vanishes identically."""
    P = sp.Poly(sp.expand((Rm - lam * gm).det(method="berkowitz")), lam)
    if P.is_zero:
        return None
    rts = set()
    for fac, _ in sp.factor_list(P.as_expr())[1]:
        if sp.Poly(fac, lam).degree() >= 1:
            rts |= set(sp.real_roots(sp.Poly(fac, lam)))
    return sorted(rts, key=lambda x: float(sp.N(x, 30)))


def form3(Rm, gm):
    """The algebraic stage at one point. Rm, gm are 4 x 4 with exact numeric entries. Returns one dict per real root
    lam of det(R - lam g): lam, the inertia (n+, n-), the rank and the admissible patterns, 'A' for (+,+,-) (at most
    two positive and one negative square) and 'B' for (+,-,-)."""
    out = []
    for lr in real_lams(Rm, gm) or []:
        M = (Rm - lr * gm) if lr.is_Rational else (Rm - lr * gm).applyfunc(lambda e: sp.N(e, 60))
        npos, nneg, nzero = inertia(M)
        pats = [p for p, (mp_, mn_) in (("A", (2, 1)), ("B", (1, 2))) if npos <= mp_ and nneg <= mn_]
        out.append({"lam": lr, "sig": (npos, nneg), "rank": npos + nneg, "pats": pats})
    return out


def split3(M, simp=sp.simplify, sampler=None):
    """Write a symmetric M of rank <= 3 as sum_i s_i a_i a_i^T by Lagrange's reduction. Returns [(s_i, a_i)], or None if
    the sign of a pivot is not decided. With X a vector where c = M(X, X) != 0 and q = M X, M - q q^T / c has rank one
    less than M."""
    n = M.shape[0]
    out = []
    M = M.applyfunc(simp)
    while not all(x == 0 for x in M):
        i = next((k for k in range(n) if M[k, k] != 0), None)
        X = sp.zeros(n, 1)
        if i is None:
            i, j = next((p, q) for p in range(n) for q in range(p + 1, n) if M[p, q] != 0)
            X[i] = X[j] = 1
        else:
            X[i] = 1
        q = M * X
        c = simp((X.T * M * X)[0])
        s = ts.sign_of(c, sampler)
        if s is None or s == 0:
            return None
        a = (q / ts.ssqrt(s * c)).applyfunc(simp)
        M = (M - s * a * a.T).applyfunc(simp)
        out.append((s, a))
        if len(out) > 3:
            return None
    return out


# ---------------------------------------------------------------------------------------------------------------
# The closure, potential and residual stages
# ---------------------------------------------------------------------------------------------------------------
def kernel(M, simp):
    """The kernel K of a rank-3 4 x 4 symmetric M, normalized by K^j = 1 at the first j with adj(M)_jj != 0. The
    adjugate of a rank-3 matrix is kappa K K^T. Returns (K, j) or None."""
    adj = M.adjugate().applyfunc(simp)
    j = next((k for k in range(M.shape[0]) if adj[k, k] != 0), None)
    if j is None:
        return None
    K = (adj[:, j] / adj[j, j]).applyfunc(simp)
    return K, j


def lie_K(M, K, Dfun, simp):
    """(L_K M)_{mu nu} = K^a d_a M_{mu nu} + M_{a nu} d_mu K^a + M_{mu a} d_nu K^a."""
    n = M.shape[0]
    return sp.Matrix(n, n, lambda m, v: simp(sum(K[a] * Dfun(M[m, v], a) + M[a, v] * Dfun(K[a], m)
                                                + M[m, a] * Dfun(K[a], v) for a in range(n))))


def ricci_pt(h, Dfun, pt):
    """Ricci tensor of the n x n metric h at the point pt in exact rationals. Dfun(E, a) is d_a E, the first and second
    derivatives of h are evaluated at pt, and the Christoffel formula is applied in exact arithmetic (route E of
    phantom_source.py, for any n)."""
    n = h.shape[0]

    def val(e):
        r = sp.cancel(sp.sympify(e).subs(pt))
        assert r.is_Rational, r
        return r
    g = h.applyfunc(val)
    dg = [[[val(Dfun(h[b, c], a)) for c in range(n)] for b in range(n)] for a in range(n)]
    ddg = [[[[val(Dfun(Dfun(h[b, c], e), a)) for c in range(n)] for b in range(n)] for e in range(n)]
           for a in range(n)]
    gi = g.inv()
    dgi = [-gi * sp.Matrix(n, n, lambda b, c: dg[a][b][c]) * gi for a in range(n)]
    Gam = [[[sum(gi[k, d] * (dg[b][d][c] + dg[c][d][b] - dg[d][b][c]) for d in range(n)) / 2
             for c in range(n)] for b in range(n)] for k in range(n)]
    dGam = [[[[sum(dgi[a][k, d] * (dg[b][d][c] + dg[c][d][b] - dg[d][b][c])
                   + gi[k, d] * (ddg[a][b][d][c] + ddg[a][c][d][b] - ddg[a][d][b][c]) for d in range(n)) / 2
               for c in range(n)] for b in range(n)] for k in range(n)] for a in range(n)]
    return sp.Matrix(n, n, lambda b, c: sum(dGam[a][a][b][c] - dGam[c][a][a][b] for a in range(n))
                     + sum(Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][a][b]
                           for a in range(n) for d in range(n))).applyfunc(sp.cancel)


def closure3(M, Dfun, simp=sp.simplify, pt=None, lam_expr=None):
    """Closure and potential stages for a rank-3 M with a one-dimensional kernel K.

    Returns a dict: 'K', 'slice' (the coordinate index j with K^j = 1), 'lie' (is L_K M = 0), 'ric' (the 3D Ricci tensor
    of the block of M with row and column j removed, at fixed x^j), 'flat'. With pt the derivatives are evaluated at that
    exact point, otherwise symbolically. lam_expr, if given, adds 'klam' = K(lam)."""
    n = M.shape[0]
    kj = kernel(M, simp)
    if kj is None:
        return {"status": "no kernel", "lie": False, "flat": False}
    K, j = kj
    LK = lie_K(M, K, Dfun, simp)
    idx = [i for i in range(n) if i != j]
    h = M.extract(idx, idx)

    def Dsl(E, a):
        return Dfun(E, idx[a])
    if pt is None:
        ric = ps.ricci(h, list(range(n - 1)), Dsl, simp)
        lie_zero = all(simp(x) == 0 for x in LK)
        flat = all(simp(x) == 0 for x in ric)
    else:
        ric = ricci_pt(h, Dsl, pt)
        LKv = LK.applyfunc(lambda e: sp.cancel(sp.sympify(e).subs(pt)))
        lie_zero = all(x == 0 for x in LKv)
        flat = all(x == 0 for x in ric)
    out = {"K": K, "slice": j, "lie": lie_zero, "ric": ric, "flat": flat, "LK": LK, "h": h}
    if lam_expr is not None:
        klam = simp(sum(K[a] * Dfun(lam_expr, a) for a in range(n)))
        out["klam"] = klam if pt is None else sp.cancel(sp.sympify(klam).subs(pt))
    return out


def einstein(g, X, simp):
    Ric = ps.ricci(g, X, lambda E, k: sp.diff(E, X[k]), simp)
    gi = g.inv().applyfunc(simp)
    Rs = simp(sum(gi[i, j] * Ric[i, j] for i in range(4) for j in range(4)))
    return Ric, (Ric - g * Rs / 2).applyfunc(simp)


# ---------------------------------------------------------------------------------------------------------------
# Planted three-scalar metrics: ds^2 = dl^2 + a(l)^2 (-dt^2 + dR^2 + f(R)^2 dTheta^2), a line warped over a fiber. The
# Ricci tensor is computed in these coordinates (To, R, Tho, Lo) and moved by the Jacobian to the sheared coordinates
# (T, R, Th, L): t = T + L/3, Theta = Th + L/5, l = L, so that the kernel direction of M is not a coordinate vector.
# For the flat fiber (f = R) the fields phi^I = kappa X^I / sqrt(8 pi), X = (t, R cos Theta, R sin Theta), have signs
# s_I. E = exp(L) is a symbol with dE/dL = E, so that every tensor is a rational function.
# ---------------------------------------------------------------------------------------------------------------
T, Th, L = sp.symbols("T Th L", real=True)
To, Tho, Lo = sp.symbols("To Tho Lo", real=True)
Rr = sp.Symbol("R", positive=True)
Ee = sp.Symbol("E", positive=True)
XP = [T, Rr, Th, L]
XO = [To, Rr, Tho, Lo]
JAC = sp.Matrix([[1, 0, 0, sp.Rational(1, 3)], [0, 1, 0, 0], [0, 0, 1, sp.Rational(1, 5)], [0, 0, 0, 1]])
t_old, Theta_old = T + L / 3, Th + L / 5
simp_r = lambda e: sp.cancel(sp.together(e))


def Dnew(E_, k):
    """d_k in the coordinates (T, R, Th, L), E = exp(L) a symbol."""
    if k == 3:
        return sp.diff(E_, L) + Ee * sp.diff(E_, Ee)
    return sp.diff(E_, XP[k])


def Dold(E_, k):
    if k == 3:
        return sp.diff(E_, Lo) + Ee * sp.diff(E_, Ee)
    return sp.diff(E_, XO[k])


def warp(kind, onshell=True):
    """Warp factor a, signs s, planted lam, kappa. 'P1' a = l (cone: lam = 0, signs (+,-,-), kappa^2 = 2). 'P2'
    a = cosh(l)/sqrt(2) (lam = -3, signs (-,+,+), kappa = 1). 'D2' adds e^(2 l)/5 inside the bracket, off shell."""
    if kind == "P1":
        return Lo, (1, -1, -1), 0, sp.sqrt(2)
    a = (Ee + 1 / Ee + (Ee**2 / 5 if kind == "D2" and not onshell else 0)) / (2 * sp.sqrt(2))
    return a, (-1, 1, 1), -3, sp.Integer(1)


def planted_geometry(a, eps=0):
    """g, Ricci and Einstein tensors in the sheared coordinates, from the diagonal old-coordinates metric."""
    f2 = Rr**2 * (1 + eps * Rr**2)
    go = sp.diag(-a**2, a**2, a**2 * f2, 1).applyfunc(sp.cancel)
    Ric = ps.ricci(go, XO, Dold, simp_r)
    gi = go.inv().applyfunc(simp_r)
    Rs = simp_r(sum(gi[i, j] * Ric[i, j] for i in range(4) for j in range(4)))
    G = (Ric - go * Rs / 2).applyfunc(simp_r)
    tr = lambda m: (JAC.T * m.subs(Lo, L) * JAC).applyfunc(simp_r)
    return tr(go), tr(Ric), tr(G)


def planted_fields(kap):
    c8 = kap / sp.sqrt(8 * sp.pi)
    return [c8 * t_old, c8 * Rr * sp.cos(Theta_old), c8 * Rr * sp.sin(Theta_old)]


def samp_planted():
    """A random exact point. E = exp(L) is rational, and L = log(E) where a formula needs it."""
    q = sp.Rational(random.randint(11, 29), random.randint(7, 10))
    return {Rr: sp.Rational(random.randint(5, 30), 10), Ee: q, L: sp.log(q), T: rq(), Th: rq()}


def samp_poly():
    return {Rr: sp.Rational(random.randint(5, 30), 10), L: sp.Rational(random.randint(5, 30), 10), T: rq(), Th: rq()}


def red_trig(e):
    """e with cos(Theta), sin(Theta) as symbols reduced by S^2 = 1 - C^2: its numerator after cancel."""
    C_, S_ = sp.symbols("Cth Sth")
    e = sp.cancel(sp.expand(e).subs({sp.cos(Theta_old): C_, sp.sin(Theta_old): S_}))
    num, _ = sp.fraction(sp.together(e))
    return sp.expand(sp.rem(sp.Poly(sp.expand(num), S_).as_expr(), S_**2 + C_**2 - 1, S_))


def planted_stages(tag, kind, eps, mut, onshell=True, fields=True):
    """The stages on the planted metric of the given kind, with a fiber curvature eps. Emits the checks `<tag>_*`."""
    a, s, lam0, kap = warp(kind, onshell)
    t0 = time.time()
    g, Ric, G = planted_geometry(a, eps)
    Dfun = Dnew
    samp = samp_poly if kind == "P1" else samp_planted
    res = []
    for _ in range(3):  # the algebraic stage at exact points
        p = samp()
        pr = lambda e: sp.Rational(sp.cancel(sp.sympify(e).subs(p)))
        res.append(form3(Ric.applyfunc(pr), g.applyfunc(pr)))
    admiss = [[c for c in cands if c["rank"] == 3 and c["pats"]] for cands in res]
    if not fields:
        ok_alg = all(len(cc) >= 1 for cc in admiss)
        detail = (f"admissible rank-3 roots at 3 exact points: {[[sp.N(c['lam'], 5) for c in cc] for cc in admiss]}, "
                  f"inertia {[[c['sig'] for c in cc] for cc in admiss]}")
    else:
        ok_alg = all(any(c["lam"] == lam0 and c["sig"] in ((2, 1), (1, 2)) for c in cc) for cc in admiss)
        detail = (f"lam = {lam0} is a real root with a rank-3 admissible form at 3 exact points, inertia "
                  f"{sorted(set(c['sig'] for cc in admiss for c in cc if c['lam'] == lam0))}; admissible rank-3 roots "
                  f"{sorted(set(float(sp.N(c['lam'], 4)) for cc in admiss for c in cc))}")
    check(f"{tag}_algebraic", ok_alg, detail + f" ({time.time() - t0:.0f} s)")
    if not ok_alg:
        return None
    # the root as a function of position, verified as a root: the kernel of R - lam g is d_l
    lam_expr = simp_r(sp.Integer(lam0) if kind != "D2" else -3 * Dnew(Dnew(a.subs(Lo, L), 3), 3) / a.subs(Lo, L))
    M = (Ric - lam_expr * g).applyfunc(simp_r)
    ok_root = simp_r(M.det()) == 0
    cl = closure3(M, Dfun, simp_r, lam_expr=lam_expr)
    if tag in ("d1", "d2"):
        if tag == "d1":
            ok = ok_root and cl["lie"] and not cl["flat"]
            why = f"lam a root: {ok_root}; L_K M = 0: {cl['lie']}; quotient flat: {cl['flat']}"
        else:
            ok = ok_root and not cl["lie"] and cl["klam"] != 0
            why = f"lam a root: {ok_root}; L_K M = 0: {cl['lie']}; K(lam) = {cl['klam']}"
        check(f"{tag}_closure_fails", ok, why + f" ({time.time() - t0:.0f} s)")
        return cl
    ok_closed = ok_root and cl["lie"] and cl["flat"]
    phis = planted_fields(kap)
    lam_detail = ""
    if ok_closed:
        spl = split3(M, simp_r, lambda: samp())
        if spl is None or len(spl) != 3:
            ok_closed, lam_detail = False, "; the split into three squares failed"
        else:
            A = sp.Matrix([list(a_.T) for _, a_ in spl])
            sg = [s_ for s_, _ in spl]
            J = sp.Matrix([[sp.diff(p_, x) * sp.sqrt(8 * sp.pi) for x in XP] for p_ in phis]).applyfunc(
                lambda e: e.subs(Ee, sp.exp(L)))
            Lam = (J * A.T * (A * A.T).inv().applyfunc(simp_r)).applyfunc(sp.expand)
            err1 = (Lam * A - J)
            err2 = (Lam.T * sp.diag(*s) * Lam - sp.diag(*sg))
            err3 = (J.T * sp.diag(*s) * J - M)
            ok_lam = all(red_trig(x) == 0 for x in list(err1) + list(err2) + list(err3))
            p1_, p2_ = samp(), samp()
            moves = any(sp.N((Lam.subs(p1_) - Lam.subs(p2_))[i, j], 12) != 0 for i in range(3) for j in range(3))
            ok_closed = ok_closed and ok_lam
            lam_detail = (f"; split signs {sg}; Lambda a = d phi, Lambda^T S Lambda = S' and d phi^T S d phi = M: "
                          f"{ok_lam}; Lambda varies with the point: {moves}")
    check(f"{tag}_closed", ok_closed, f"K = {list(cl['K'])}, L_K M = 0: {cl['lie']}, 3D Ricci of the quotient block "
          f"zero: {cl['flat']}{lam_detail} ({time.time() - t0:.0f} s)")
    if not ok_closed:
        return cl
    pot_lam = lam_expr + (L / 10 if mut == "p2-lambda" and kind == "P2" else 0)
    klam = simp_r(sum(cl["K"][k] * Dfun(pot_lam, k) for k in range(4)))
    check(f"{tag}_potential", klam == 0, f"K(lam) = {klam}, lam = {lam_expr}")
    sres = list(s)
    if mut == "p1-sign" and kind == "P1":
        sres[0] = -sres[0]
    u = [sp.Matrix([sp.diff(p_, x) for x in XP]).applyfunc(lambda e: e.subs(Ee, sp.exp(L))) * sp.sqrt(8 * sp.pi)
         for p_ in phis]
    gi = g.inv().applyfunc(simp_r)
    kin = [simp_r((uu.T * gi * uu)[0]) for uu in u]
    Tab = sum((sres[i] * u[i] * u[i].T for i in range(3)), sp.zeros(4, 4)) \
        - g * (sum(sres[i] * kin[i] for i in range(3)) / 2 + lam_expr)
    resid = (G - Tab)
    ok_res = all(red_trig(x) == 0 for x in resid)
    check(f"{tag}_residual", ok_res, f"G_ab - 8 pi T_ab with the fields substituted back: "
          f"{'0' if ok_res else [red_trig(x) for x in list(resid)[:3]]} ({time.time() - t0:.0f} s)")
    return cl


# ---------------------------------------------------------------------------------------------------------------
# The wall
# ---------------------------------------------------------------------------------------------------------------
RW = sp.Symbol("r", positive=True)
WJ = sp.symbols("W0:8", real=True)  # WJ[k] = w^(k)(r), w = v (1 - f)
CC, SS = sp.symbols("c S", real=True)


def wall_frame(mut):
    """The frame Ricci tensor of the comoving spherical route of phantom_source.py in (r, c, S, W0..W2)."""
    Rs = ps.sph_tetrad_ricci("alc-drop-w2" if mut == "alc-drop-w2" else None)
    sub = {ps.W[k]: WJ[k] for k in range(3)}
    return Rs.applyfunc(lambda e: sp.expand(sp.expand_trig(e).subs({sp.sin(ps.th): SS, sp.cos(ps.th): CC})
                                            .subs(ps.r_, RW).subs(sub)))


def orbit_metric(R3, w):
    """h = R_ab theta^a theta^b on the orbit space (t, r, theta) of the axisymmetric stationary comoving frame, with
    theta^0 = dt, theta^1 = d xi + w dt = w dt + c dr - r S d theta, theta^rho = d rho = S dr + r c d theta."""
    E = sp.Matrix([[1, 0, 0], [w, CC, -RW * SS], [0, SS, RW * CC]])
    return (E.T * R3 * E).applyfunc(sp.expand)


def D_jets(E, k):
    """d_k in the coordinates (t, r, theta, varphi), on expressions in r, c, S and the jets W_k of w."""
    if k in (0, 3):
        return sp.Integer(0)
    if k == 1:
        return sp.expand(sp.diff(E, RW) + sum(WJ[i + 1] * sp.diff(E, WJ[i]) for i in range(len(WJ) - 1)))
    return sp.expand(-SS * sp.diff(E, CC) + CC * sp.diff(E, SS))


def D_expl(E, k):
    """The same for an explicit profile: no jet symbols."""
    if k in (0, 3):
        return sp.Integer(0)
    if k == 1:
        return sp.expand(sp.diff(E, RW))
    return sp.expand(-SS * sp.diff(E, CC) + CC * sp.diff(E, SS))


def block4(h):
    M = sp.zeros(4, 4)
    M[:3, :3] = h
    return M


def jets_of(vfun, r0, n=7):
    """w = v (1 - f): W0 = v (1 - f), W_k = -v f^(k)."""
    return [vfun[0] * (1 - vfun[1][0])] + [-vfun[0] * vfun[1][k] for k in range(1, n)]


def tanh_jets_n(r0, sig, Rb, n):
    mp.mp.dps = 60
    fn = lambda r: (mp.tanh(sig * (r + Rb)) - mp.tanh(sig * (r - Rb))) / (2 * mp.tanh(sig * Rb))
    return [sp.Rational(mp.nstr(mp.diff(fn, mp.mpf(r0.p) / r0.q, k), 55)) for k in range(n)]


def bump_jets_n(r0, Rb, n):
    rr = sp.Symbol("rr")
    fb = (1 - (rr / Rb)**2)**3
    return [sp.diff(fb, rr, k).subs(rr, r0) if r0 < Rb else 0 for k in range(n)]


# ---------------------------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------------------------
def main():
    global MUT
    opts = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
    mut = MUT = opts.get("mutate")
    if mut is not None and mut not in MUTATIONS:
        sys.exit(f"unknown mutation {mut!r}. Known: {', '.join(MUTATIONS)}")
    random.seed(77)
    t0 = time.time()
    stop = int(opts.get("stop", 99))
    print(f"three_scalars.py, issue #77. sympy {sp.__version__}, mpmath {mp.__version__}, python "
          f"{sys.version.split()[0]}. mutation: {mut or 'none'}" + (f" ({MUTATIONS[mut][1]})" if mut else ""))

    # 1. Planted truth for the algebraic stage, in random Lorentz frames.
    print("\n== 1. Planted truth: R = lam eta + sum s_I a_I a_I passes for (+,+,-) and (+,-,-), other forms fail ==")
    bad = []
    for k in range(12):
        L_, lv = lorentz(), rq()
        pat = "A" if k % 2 == 0 else "B"
        sg = (1, 1, -1) if pat == "A" else (1, -1, -1)
        if mut == "planted3-sign":
            sg = (1, 1, 1)
        vs = [vec4() for _ in range(3)]
        Rm = L_.T * (lv * eta + sum((sg[i] * vs[i] * vs[i].T for i in range(3)), sp.zeros(4, 4))) * L_
        res = form3(Rm, eta)
        ok = any(c["lam"] == lv and pat in c["pats"] and c["rank"] == 3 for c in res)
        bad += [] if ok else [(lv, pat, [(c["lam"], c["sig"]) for c in res if c["lam"] == lv])]
    check("planted3_pass", not bad, f"12 random lam eta + sum s_I a_I a_I, six of each pattern, transformed: lam "
          f"recovered with rank 3 and the planted pattern in all but {len(bad)} {bad[:2]}")
    bad = []
    kinds = [(1, 1, 1, 0), (-1, -1, -1, 0), (1, 1, -1, -1), (1, 1, 1, -1), (-1, -1, -1, 1), (1, -1, -1, -1)]
    for kind in kinds * 2:
        L_, lv = lorentz(), rq()
        sg = (1, 1, -1, 0) if mut == "planted3-fail-sign" else kind
        vs = [vec4() for _ in range(4)]
        Rm = L_.T * (lv * eta + sum((sg[i] * vs[i] * vs[i].T for i in range(4)), sp.zeros(4, 4))) * L_
        res = form3(Rm, eta)
        ok = not any(c["lam"] == lv and c["pats"] for c in res)
        bad += [] if ok else [(kind, lv)]
    check("planted3_fail", not bad, f"12 random forms with squares (+,+,+), (-,-,-) and four squares: the planted lam "
          f"is admissible in all but {len(bad)} {bad[:2]}")

    if stop < 2:
        return fin(t0, mut)
    # 2. The lemma behind the wall's algebraic stage.
    print("\n== 2. Lemma: with R_a3 = 0 and (I2), lam = 0 gives rank 3 iff det R_3 != 0, pattern (+,+,-) iff det R_3 < 0 ==")
    wt = (1, -1, 2) if mut == "lemma3-weight" else (1, 1, 2)
    mism, n0, nA, nB, nlam = [], 0, 0, 0, 0
    for trial in range(24):
        wn = lambda u_: wt[0] * u_[0]**2 + wt[1] * u_[1]**2 + wt[2] * u_[2]**2
        if trial % 4 == 0:  # a rank-2 R_3 = a a^T - t b b^T, with t fixed by the weighted trace
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
        d3 = R3.det()
        Rm = sp.zeros(4, 4)
        Rm[:3, :3] = R3
        L_ = lorentz()
        res = form3(L_.T * Rm * L_, eta)
        z = [c for c in res if c["lam"] == 0]
        ok = bool(z)
        if ok:
            c0 = z[0]
            if d3 == 0:
                n0 += 1
                ok = c0["rank"] <= 2
            else:
                ok = c0["rank"] == 3 and c0["pats"] == (["A"] if d3 < 0 else ["B"])
                nA += int(bool(d3 < 0))
                nB += int(bool(d3 > 0))
        for c in res:  # a nonzero admissible root has the sign of its pattern
            if c["lam"] != 0 and c["pats"] and c["rank"] == 3:
                nlam += 1
                ok = ok and ((c["lam"] < 0 and c["pats"] == ["A"]) or (c["lam"] > 0 and c["pats"] == ["B"]))
        if not ok:
            mism.append((trial, d3, [(c["lam"], c["sig"], c["pats"]) for c in res][:3]))
    check("lemma3_samples", not mism, f"{nA} pattern (+,+,-) with det R_3 < 0, {nB} pattern (+,-,-) with det R_3 > 0, "
          f"{n0} with det R_3 = 0 of rank <= 2, {nlam} admissible roots lam != 0, each of the sign of its pattern; "
          f"{len(mism)} disagree with the lemma {mism[:2]}")

    if stop < 3:
        return fin(t0, mut)
    # 3. Planted three-scalar controls and deformed ones, through the same stages.
    print("\n== 3. Planted controls: warped products over a flat fiber, sheared polar coordinates (T, R, Th, L) ==")
    print("  P1: ds^2 = dl^2 + l^2 (-dt^2 + dR^2 + R^2 dTheta^2), phi = sqrt(2) (t, R cos Theta, R sin Theta) / sqrt(8 pi), "
          "signs (+,-,-), lam = 0")
    print("  P2: a = cosh(l)/sqrt(2), phi = (t, R cos Theta, R sin Theta) / sqrt(8 pi), signs (-,+,+), lam = -3")
    planted_stages("p1", "P1", 0, mut)
    planted_stages("p2", "P2", sp.Rational(1, 10) if mut == "p2-fiber" else 0, mut)
    print(f"  ({time.time() - t0:.0f} s)")
    print("\n== 3b. Deformed controls: D1 curves the fiber (R^2 -> R^2 (1 + R^2 / 10)), D2 puts a(l) off shell ==")
    planted_stages("d1", "P1", 0 if mut == "d1-flat" else sp.Rational(1, 10), mut, fields=False)
    planted_stages("d2", "D2", 0, mut, onshell=(mut == "d2-onshell"), fields=False)
    print(f"  ({time.time() - t0:.0f} s)")

    if stop < 4:
        return fin(t0, mut)

    # 4. The wall, algebraic stage.
    print("\n== 4. Alcubierre wall, constant v: the frame block R_3 on (E_0 = n, E_1 = d_x, e_rho), general profile ==")
    Rs = wall_frame(mut)
    r_, W0, W1, W2 = RW, WJ[0], WJ[1], WJ[2]
    R3 = Rs[:3, :3]
    redsc = lambda e: sp.expand(sp.rem(sp.Poly(sp.expand(e), SS).as_expr(), SS**2 + CC**2 - 1, SS))
    ids = [redsc(sp.simplify(Rs[i, 3])) for i in range(4)] + [redsc(sp.simplify(R3[0, 0] + R3[1, 1] + 2 * R3[2, 2]))]
    check("wall_identities", all(x == 0 for x in ids), f"R(e_phi, .) and R_00 + R_11 + 2 R_rho,rho: {ids} "
          f"({time.time() - t0:.0f} s)")
    # 16 r^3 D, D = det R_3, as a polynomial in c and S reduced by S^2 = 1 - c^2
    Dn, Dd = sp.fraction(sp.cancel(16 * r_**3 * R3.det(method="berkowitz")))
    Dn = sp.expand(sp.rem(sp.Poly(sp.expand(Dn), SS).as_expr(), SS**2 + CC**2 - 1, SS))
    q, rem_ = sp.div(sp.Poly(Dn.coeff(SS, 0), CC), sp.Poly(1 - CC**2, CC))
    only_even = all(k % 2 == 0 for (k,) in sp.Poly(Dn.coeff(SS, 0), CC).monoms()) and Dn.coeff(SS, 1) == 0
    check("wall_D", rem_.is_zero and only_even and Dd == 1,
          f"16 r^3 det R_3 = (1 - c^2) q(c^2) with q = {sp.factor(q.as_expr())}; no odd terms: {only_even}. So R_3 "
          f"(+) 0 has rank 3 exactly where q != 0 off the axis, with the pattern fixed by the sign of q")
    qpoly = sp.lambdify((r_, W0, W1, W2, CC), q.as_expr(), "sympy")
    print(f"  ({time.time() - t0:.0f} s)")

    # 4b. Exact points of the wall: the lab frame of phantom_source.py, the same form3 as the planted controls.
    print("\n== 4b. Exact points (lab frame route of phantom_source.py): the algebraic stage through form3 ==")
    Rt = ps.lab_tetrad_ricci(None)
    random.seed(57)
    samples = []
    wall = [(12, 3, 4, 13), (-12, 4, 3, 13), (13, 0, 0, 13), (-13, 0, 0, 13), (0, 5, 12, 13), (5, 12, 0, 13),
            (-5, 0, 12, 13), (9, 12, 0, 15)]
    for (x0, y0, z0, r0) in wall:
        pnt = [sp.Rational(q_, r0) for q_ in (x0, y0, z0, r0)]
        for vv in (sp.Integer(1), sp.Rational(1, 2), sp.Integer(3)):
            samples.append(("tanh sig=8 Rb=1", pnt, vv, tanh_jets_n(pnt[3], 8, 1, 7)))
    for (x0, y0, z0, r0) in [(4, 3, 0, 5), (-3, 0, 4, 5), (5, 0, 0, 5), (0, 3, 4, 5)]:
        pnt = [sp.Rational(q_, 5) for q_ in (x0, y0, z0, r0)]
        samples.append(("bump Rb=2", pnt, sp.Rational(7, 4), bump_jets_n(pnt[3], 2, 7)))
    for (x0, y0, z0, r0) in [(12, 3, 4, 13), (-12, 4, 3, 13), (4, 3, 0, 5), (-3, 0, 4, 5)]:
        pnt = [sp.Rational(q_, 7) for q_ in (x0, y0, z0, r0)]
        samples.append(("random jets", pnt, sp.Rational(random.randint(1, 9), 3),
                        [sp.Rational(random.randint(-9, 9), random.randint(1, 9)) for _ in range(7)]))
    random.seed(770)
    rows, nA, nB, nax, nlam = [], 0, 0, 0, 0
    lam_cases = []
    allok = True
    for (lab_, pnt, vv, jets) in samples:
        sub = {ps.X: pnt[0], ps.y: pnt[1], ps.z: pnt[2], ps.R: pnt[3], ps.v: vv, ps.F[0]: jets[0], ps.F[1]: jets[1],
               ps.F[2]: jets[2]}
        res = form3(Rt.subs(sub), eta)
        axis = pnt[1] == 0 and pnt[2] == 0
        cen = sp.sqrt(pnt[1]**2 + pnt[2]**2) / pnt[3]
        Wv = [vv * (1 - jets[0]), -vv * jets[1], -vv * jets[2]]
        dv = (1 - (pnt[0] / pnt[3])**2) * qpoly(pnt[3], Wv[0], Wv[1], Wv[2], pnt[0] / pnt[3]) / (16 * pnt[3]**3)
        z = [c for c in res if c["lam"] == 0]
        if axis:
            ok = bool(z) and z[0]["rank"] <= 2 and dv == 0
            nax += 1
        else:
            ok = bool(z) and z[0]["rank"] == 3 and z[0]["pats"] == (["A"] if dv < 0 else ["B"])
            nA += int(bool(dv < 0))
            nB += int(bool(dv > 0))
        extra = [c for c in res if c["lam"] != 0 and c["rank"] == 3 and c["pats"]]
        nlam += len(extra)
        if not axis:
            lam_cases += [(pnt, vv, jets, c) for c in extra]
        allok &= ok
        rows.append((lab_, axis, "A" if dv < 0 else ("B" if dv > 0 else "axis"), len(extra)))
    for key in sorted(set((r[0], r[1]) for r in rows)):
        sub_ = [r for r in rows if (r[0], r[1]) == key]
        print(f"  {key[0]:>16}, {'on the axis' if key[1] else 'off the axis'}: {len(sub_)} points, lam = 0 patterns "
              f"{sorted(set(r[2] for r in sub_))}, admissible rank-3 roots lam != 0 per point {[r[3] for r in sub_]}")
    check("wall_points", allok and nA >= 1 and nB >= 1,
          f"{len(samples)} exact points: lam = 0 is a root at all; off the axis rank 3 with pattern (+,+,-) at {nA} "
          f"(det R_3 < 0) and (+,-,-) at {nB} (det R_3 > 0), predicted by the sign of D at every point; on the axis "
          f"rank <= 2 at {nax}. Admissible rank-3 roots lam != 0 at the points: {nlam} ({time.time() - t0:.0f} s)")
    if stop < 5:
        return fin(t0, mut)

    # 5. The wall, closure at lam = 0.
    print("\n== 5. Closure at lam = 0: K = d_varphi, the block is h = R_ab theta^a theta^b on (t, r, theta) ==")
    hsym = orbit_metric(R3, W0)
    nflat, ntot, nlie = 0, 0, 0
    for (lab_, pnt, vv, jets) in samples:
        if pnt[1] == 0 and pnt[2] == 0:
            continue
        cc, ss = pnt[0] / pnt[3], sp.sqrt(pnt[1]**2 + pnt[2]**2) / pnt[3]
        pt = {RW: pnt[3], CC: cc, SS: ss}
        pt.update({WJ[k]: (vv * (1 - jets[0]) if k == 0 else -vv * jets[k]) for k in range(7)})
        hh = sp.diag(-1, 1, 1) if mut == "wall-flat-sub" else hsym
        cl = closure3(block4(hh), D_jets, sp.expand, pt=pt)
        ntot += 1
        nflat += int(cl["flat"])
        nlie += int(cl["lie"])
    check("wall_flat_points", ntot > 0 and nflat == 0 and nlie == ntot,
          f"K = d_varphi, L_K M = 0 at {nlie} of {ntot} off-axis exact points, and the Ricci tensor of the 3D block h is "
          f"nonzero at {ntot - nflat} of {ntot} ({time.time() - t0:.0f} s)")
    # explicit profiles, the same function: the polynomial bump f = (1 - r^2/4)^3 with several speeds
    rr = RW
    nfe, nte = 0, 0
    for vv in (sp.Integer(1), sp.Rational(1, 2), sp.Rational(7, 4)):
        wexp = vv * (1 - (1 - rr**2 / 4)**3)
        sub = {WJ[k]: sp.diff(wexp, rr, k) for k in range(3)}
        R3e = R3.subs(sub)
        he = orbit_metric(R3e, wexp)
        for (cc, ss) in ((sp.Rational(3, 5), sp.Rational(4, 5)), (sp.Rational(5, 13), sp.Rational(12, 13)),
                         (sp.Rational(-4, 5), sp.Rational(3, 5))):
            for r0 in (sp.Rational(1), sp.Rational(3, 2)):
                clx = closure3(block4(he), D_expl, sp.expand, pt={RW: r0, CC: cc, SS: ss})
                nte += 1
                nfe += int(clx["flat"])
    check("wall_flat_explicit", nte > 0 and nfe == 0, f"the explicit profile w = v (1 - (1 - r^2/4)^3), v in {{1, 1/2, "
          f"7/4}}: the 3D block has nonzero Ricci tensor at {nte - nfe} of {nte} exact points ({time.time() - t0:.0f} s)")
    if stop < 6:
        return fin(t0, mut)

    # 6. The wall, lam != 0: the potential stage, K(lam) = 0, at the admissible rank-3 roots of 4b.
    print("\n== 6. lam != 0 branches: at each admissible root, K(lam) with lam from the implicit function theorem ==")
    eta3 = sp.diag(-1, 1, 1)
    Pl = sp.expand((R3 - lam * eta3).det(method="berkowitz"))
    Pr, Pth, Plam = D_jets(Pl, 1), D_jets(Pl, 2), sp.diff(Pl, lam)
    rels, nfail = [], 0
    for (pnt, vv, jets, c) in lam_cases:
        mp.mp.dps = 80
        cc, ss = pnt[0] / pnt[3], sp.sqrt(pnt[1]**2 + pnt[2]**2) / pnt[3]
        pt = {RW: pnt[3], CC: cc, SS: ss}
        pt.update({WJ[k]: (vv * (1 - jets[0]) if k == 0 else -vv * jets[k]) for k in range(4)})
        lamv = sp.N(c["lam"], 80)
        Nn = (R3 - lamv * eta3).subs(pt)
        A = mp.matrix(3, 3)
        for i in range(3):
            for j in range(3):
                A[i, j] = mp.mpf(sp.N(Nn[i, j], 80))
        ev, V = mp.eigsy(A)
        k = min(range(3), key=lambda i: abs(ev[i]))
        Kv = [V[i, k] for i in range(3)]
        ev_ = lambda e: mp.mpf(sp.N(sp.sympify(e).subs(pt).subs(lam, lamv), 80))
        lam_r, lam_t = -ev_(Pr) / ev_(Plam), -ev_(Pth) / ev_(Plam)
        if mut == "lam-const":
            lam_r = lam_t = mp.mpf(0)
        w_, c_, s_, r0 = [mp.mpf(sp.N(x, 80)) for x in (pt[WJ[0]], cc, ss, pnt[3])]
        klam = (Kv[1] - w_ * Kv[0]) * (c_ * lam_r - s_ / r0 * lam_t) + Kv[2] * (s_ * lam_r + c_ / r0 * lam_t)
        rel = abs(klam) * r0 / (abs(mp.mpf(sp.N(lamv, 80))) * mp.sqrt(sum(x * x for x in Kv)))
        rels.append(float(rel))
        nfail += int(rel > mp.mpf(10) ** -8)
    check("wall_lam_potential", len(lam_cases) >= 1 and nfail == len(lam_cases),
          f"{len(lam_cases)} admissible rank-3 roots lam != 0 at the off-axis exact points; K(lam) != 0 at {nfail} "
          f"(|K(lam)| r / (|lam| |K|) between {min(rels or [0]):.2e} and {max(rels or [0]):.2e}), so the potential "
          f"stage fails at each ({time.time() - t0:.0f} s). Numerical evidence, 80 digits")
    if stop < 7:
        return fin(t0, mut)

    # 7. The accelerating bubble of phantom_source.py section 6 (route D), at lam = 0.
    print("\n== 7. Accelerating bubble b = -v_s(t) f(r_s): block h on (t, x, y) at z = 0, same closure function ==")
    bsym = -ps.VS[0] * ps.F[0]
    Ra = ps.shift_tetrad_ricci(bsym, ps.bubble_rule(mut == "acc-static"), False)
    R3a = Ra[:3, :3].subs(ps.z, 0)
    Eb = sp.Matrix([[1, 0, 0], [bsym, 1, 0], [0, 0, 1]])
    hacc = (Eb.T * R3a * Eb).applyfunc(sp.expand)
    jd = ps.JetD(ps.bubble_rule(False))
    Dacc = lambda E_, k: sp.Integer(0) if k == 3 else jd(E_, k)

    def acc_pt(Xv, yv, Rv, vs, xs, fj):
        pt = {ps.X: Xv, ps.y: yv, ps.R: Rv}
        pt.update({ps.F[k]: fj[k] for k in range(5)})
        pt.update({ps.VS[k]: vs[k] for k in range(4)})
        pt.update({ps.XS[k]: xs[k] for k in range(1, 4)})
        return pt
    # planted truth: at constant speed (v_s' = 0, x_s' = v) the lab-route Ricci tensor of h is the polar one, moved by
    # the Jacobian of (t, r, theta) -> (t, x, y)
    ok_lim, nlim = True, 0
    hp = orbit_metric(R3, W0)
    for (Xv, yv, Rv) in ((sp.Rational(12, 13), sp.Rational(5, 13), sp.Integer(1)),
                         (sp.Rational(-3, 5), sp.Rational(4, 5), sp.Integer(1)),
                         (sp.Rational(4, 5), sp.Rational(3, 5), sp.Integer(1))):
        vv = sp.Rational(random.randint(1, 9), random.randint(2, 5))
        fj = [sp.Rational(random.randint(-9, 9), random.randint(1, 9)) for _ in range(5)]
        ptl = acc_pt(Xv, yv, Rv, [vv, 0, 0, 0], [0, vv, 0, 0], fj)
        ric_lab = ricci_pt(hacc, Dacc, ptl)
        ptp = {RW: Rv, CC: Xv / Rv, SS: yv / Rv}
        ptp.update({WJ[k]: (vv * (1 - fj[0]) if k == 0 else -vv * fj[k]) for k in range(5)})
        ric_pol = ricci_pt(hp, D_jets, ptp)
        Jp = sp.Matrix([[1, 0, 0], [vv, Xv / Rv, -yv], [0, yv / Rv, Xv]])
        ok_lim &= all(x == 0 for x in (ric_pol - Jp.T * ric_lab * Jp))
        nlim += 1
    check("acc_const_limit", ok_lim, f"at v_s' = 0, x_s' = v_s the lab-route Ricci tensor of the block, moved by the "
          f"Jacobian of (t, r, theta) -> (t, x, y), equals the comoving polar one at {nlim} exact points "
          f"({time.time() - t0:.0f} s)")
    nz, ntot = 0, 0
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
        cl = closure3(block4(hacc), Dacc, sp.expand, pt=acc_pt(Xv, yv, Rv, vs, xs, fj))
        ntot += 1
        nz += int(not cl["flat"])
    check("acc_flat_points", nz == ntot, f"random jets of v_s (v_s' != 0), x_s and f at {ntot} exact points of the "
          f"half-plane: the Ricci tensor of the 3D block is nonzero at {nz} ({time.time() - t0:.0f} s)")
    return fin(t0, mut)


def fin(t0, mut):
    failed = [n for n, ok in RESULTS if not ok]
    print(f"\nchecks: {len(RESULTS) - len(failed)} passed, {len(failed)} failed {failed}. {time.time() - t0:.0f} s")
    if mut is not None:
        want = MUTATIONS[mut][0]
        print(f"control {mut}: target checks {want}, fired {[n for n in want if n in failed]}")
    return failed


if __name__ == "__main__":
    failed = main()
    controls_ok = True
    if "--mutate" not in " ".join(sys.argv) and "--no-controls" not in sys.argv:
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

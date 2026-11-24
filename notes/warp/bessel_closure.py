#!/usr/bin/env python3
"""Closure fails on every open set for the Bessel shift family u = c x + g(rho), c != 0 (issue #162).

Metric ds^2 = -dt^2 + (dx + u dt)^2 + d rho^2 + rho^2 d phi^2, u = c x + g(rho), on the ODE
g'' + g'/rho = eps (g'^2 + 2 c^2), eps = +-1, whose solutions are g = -eps ln(alpha J0(sqrt2 c rho) + beta Y0(sqrt2 c rho)) + d
(det_classify.md section 3). There R = a a - b b with a = (g' d rho + 2 c th1)/sqrt2, b = W nu/sqrt2, W = sqrt(g'^2 + 2 c^2),
th0 = dt, th1 = dx + u dt, th2 = d rho, nu = th0 + eps th1. With k = a + b and m = a - b, R(xi, xi) = k(xi) m(xi), so k and m
are fixed by R up to k -> lam k, m -> m/lam and a swap. Closure of two_scalars.md (some chi with e^chi k and e^-chi m closed)
needs k ^ dk = 0 and m ^ dm = 0 (wedge omega ^ k = -dk with k, and omega ^ m = dm with m). For k~ = sqrt2 k and m~ = sqrt2 m
the only component of k ^ dk and m ^ dm is O dt ^ dx ^ d rho, with, on the ODE,

    O_k = g' (g'^2 + 4 c^2 + 5 c eps W - 2 c g'/(rho W)),   O_m = g' (g'^2 + 4 c^2 - 5 c eps W + 2 c g'/(rho W)),
    O_k + O_m = 2 g' (g'^2 + 4 c^2)   (no ODE needed).

This script has these checks:
  frame      R = a a - b b on the ODE for eps = +-1 (own Ricci tensor), and a ^ b != 0
  frobenius  the coordinate k ^ dk, m ^ dm against the structure-equation formula and the closed forms
  sum        O_k + O_m = 2 g'(g'^2 + 4 c^2), whose only real zero is g' = 0 for c != 0
  isolated   g'' = 2 eps c^2 at g' = 0, so the zeros of g' are isolated; at c = 0 g = const is a member
  each       O_k and O_m separately vanish on no interval: the elimination polynomial N(G, c) of the Lie derivative
  points     exact jets through two_scalars.closure: it fails exactly where (O_k, O_m) != (0, 0)
  bessel     actual Bessel members at 80 digits, direct k ^ dk against the closed forms
  planted    closed pairs pass, a non-integrable pair fails, an integrable pair with incompatible factors fails closure
  c0         c = 0 is #65's b(rho): O = g'^3 = the coefficient of nu ^ d nu
  planar     g' = 0 (u = c x + d, outside the ODE): O = 0 and closure fails at the second stage
Imports: two_scalars (closure, dform), read-only. Nothing else from the repository.
Usage: bessel_closure.py [--no-controls] [--mutate=NAME] [--only=check1,check2] [--logdir=DIR]. Exit 1 if a check fails or
a control does not fire on its target only."""
import os
import random
import subprocess
import sys
import time

import mpmath as mp
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import two_scalars as ts  # noqa: E402

# name: (checks that must fail, what the mutation does)
MUTATIONS = {
    "frame-nu": (["frame"], "takes nu = th0 - eps th1 in b"),
    "frob-struct": (["frobenius"], "drops the c th1 ^ th0 term of d th1 (the term c beta gamma) in the coframe formula"),
    "frob-ode": (["frobenius"], "reduces on the ODE with -eps in place of eps"),
    "sum-factor": (["sum"], "expects 2 g'(g'^2 + 2 c^2) for O_k + O_m"),
    "isolated-c0": (["isolated"], "samples the step g'' = 2 eps c^2 at c = 0"),
    "each-w2": (["each"], "reduces the elimination with W^2 = G^2 + c^2"),
    "each-friction": (["each"], "drops the term -G r of G' in the Lie derivative"),
    "points-closed": (["points"], "replaces (a, b) at every jet by the closed pair (d rho, dt)"),
    "bessel-freq": (["bessel"], "uses J0(c rho), Y0(c rho) in place of J0(sqrt2 c rho), Y0(sqrt2 c rho)"),
    "plant-pass": (["planted"], "multiplies the boosted planted m by the non-integrable (dt - dx + x d rho)"),
    "plant-fail": (["planted"], "replaces the planted non-integrable pair (dt, x d rho) by the closed (dt, rho d rho)"),
    "c0-sign": (["c0"], "expects -g'^3 in place of g'^3"),
    "planar-closed": (["planar"], "replaces the planar pair (c th1, c th0) by the closed (c dx, c dt)"),
}

RESULTS = []
OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate")
ONLY = OPTS["only"].split(",") if "only" in OPTS else None
R = sp.Rational


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def wanted(name):
    return ONLY is None or name in ONLY


# ---------------------------------------------------------------------------------------------------------------
# Setup: coordinates, the coframe, the ODE, the covectors k, m
# ---------------------------------------------------------------------------------------------------------------
t, x, rho, phi = sp.symbols("t x rho phi", positive=True)
X4 = [t, x, rho, phi]
c = sp.Symbol("c", real=True)
g = sp.Function("g")(rho)
g1, g2 = sp.diff(g, rho), sp.diff(g, rho, 2)
Gs, Ws, rs = sp.symbols("G W r", real=True)  # jets: G = g', W = sqrt(G^2 + 2 c^2), r = 1/rho


def zero(e):
    return sp.simplify(e) == 0


def coframe(u):
    """th0 = dt, th1 = dx + u dt, th2 = d rho, as covector components on (t, x, rho, phi)."""
    return sp.Matrix([1, 0, 0, 0]), sp.Matrix([u, 1, 0, 0]), sp.Matrix([0, 0, 1, 0])


def g2rule(eps, flip=False):
    """g'' on the ODE g'' + g'/rho = eps (g'^2 + 2 c^2)."""
    e = -eps if flip else eps
    return e * (g1**2 + 2 * c**2) - g1 / rho


def reduce_ode(e, eps, flip=False):
    r2 = g2rule(eps, flip)
    r3 = sp.diff(r2, rho).subs(sp.Derivative(g, (rho, 2)), r2)
    return e.subs(sp.Derivative(g, (rho, 3)), r3).subs(sp.Derivative(g, (rho, 2)), r2)


def frob(k):
    """All components (k ^ dk)_{ijl}, i < j < l, in coordinates."""
    dk = ts.dform(k, X4)
    return {(i, j, l): k[i] * dk[j][l] + k[j] * dk[l][i] + k[l] * dk[i][j]
            for i in range(4) for j in range(i + 1, 4) for l in range(j + 1, 4)}


def O_of(k):
    """The coefficient of dt ^ dx ^ d rho in k ^ dk; the three components with d phi must vanish."""
    cm = frob(k)
    assert all(zero(cm[key]) for key in cm if 3 in key), "a d phi component of k ^ dk is nonzero"
    return cm[(0, 1, 2)]


def kcov(s, eps, gfun=None, cc=None, W=None):
    """k~ = sqrt2 (a + b) for s = +1, m~ = sqrt2 (a - b) for s = -1, with a = (G th2 + 2 c th1)/sqrt2, b = W nu/sqrt2,
    G = g', W = sqrt(G^2 + 2 c^2) unless given. Defaults to the symbolic g and c."""
    gfun = g if gfun is None else gfun
    cc = c if cc is None else cc
    G = sp.diff(gfun, rho)
    W = sp.sqrt(G**2 + 2 * cc**2) if W is None else W
    th0, th1, th2 = coframe(cc * x + gfun)
    return s * W * th0 + (2 * cc + s * eps * W) * th1 + G * th2


def O_closed(s, eps, G, W, r, cc):
    """The closed forms of the module docstring, r = 1/rho; s = +1 is O_k, s = -1 is O_m."""
    return G * (G**2 + 4 * cc**2 + s * (5 * cc * eps * W - 2 * cc * G * r / W))


def ricci(gm, X):
    """Ricci tensor R_bc of a metric in coordinates (own code: Christoffel symbols, then the contraction)."""
    n = len(X)
    gi = gm.inv().applyfunc(sp.simplify)
    dg = [[[sp.diff(gm[b, cc], X[a]) for cc in range(n)] for b in range(n)] for a in range(n)]
    Gam = [[[sp.simplify(sum(gi[a, d] * (dg[b][d][cc] + dg[cc][d][b] - dg[d][b][cc]) for d in range(n)) / 2)
             for cc in range(n)] for b in range(n)] for a in range(n)]
    Ric = sp.zeros(n, n)
    for b in range(n):
        for cc in range(b, n):
            s = 0
            for a in range(n):
                s += sp.diff(Gam[a][b][cc], X[a]) - sp.diff(Gam[a][a][b], X[cc])
                for d in range(n):
                    s += Gam[a][a][d] * Gam[d][b][cc] - Gam[a][cc][d] * Gam[d][a][b]
            Ric[b, cc] = Ric[cc, b] = sp.simplify(s)
    return Ric


def shift_metric(u):
    return sp.Matrix([[-1 + u**2, u, 0, 0], [u, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, rho**2]])


# ---------------------------------------------------------------------------------------------------------------
# C1 frame
# ---------------------------------------------------------------------------------------------------------------
def c1_frame():
    u = c * x + g
    th0, th1, th2 = coframe(u)
    Rc = ricci(shift_metric(u), X4)
    ok_R, ok_w = True, True
    for eps in (1, -1):
        nu = th0 - eps * th1 if MUT == "frame-nu" else th0 + eps * th1
        a2 = g1 * th2 + 2 * c * th1  # sqrt2 a
        RR = Rc - (a2 * a2.T / 2 - (g1**2 + 2 * c**2) / 2 * nu * nu.T)
        ok_R &= all(zero(reduce_ode(e, eps)) for e in RR)
        # minors of (sqrt2 a, nu) in the coframe: a ^ b = (W/2)(minors) has squares summing to W^2 * 2 W^2 / 4 > 0
        M = sp.Matrix([[0, 2 * c, g1], [1, eps, 0]])
        sq = sum(M[:, [i, j]].det() ** 2 for i, j in ((0, 1), (0, 2), (1, 2)))
        ok_w &= zero(sq - 2 * (g1**2 + 2 * c**2))
    check("frame", ok_R and ok_w, f"R = a a - b b on the ODE for eps = +-1 (own Ricci tensor, exact): {ok_R}; the minors of "
          f"(sqrt2 a, nu) in the coframe have squares summing to 2 (g'^2 + 2 c^2) > 0, so a ^ b != 0: {ok_w}")


# ---------------------------------------------------------------------------------------------------------------
# C2 frobenius
# ---------------------------------------------------------------------------------------------------------------
def c2_frobenius():
    ok_a, ok_b = True, True
    for eps in (1, -1):
        W = sp.sqrt(g1**2 + 2 * c**2)
        for s in (1, -1):
            O = O_of(kcov(s, eps))
            al, be, ga = s * W, 2 * c + s * eps * W, g1  # k~ = al th0 + be th1 + ga th2
            cterm = 0 if MUT == "frob-struct" else c * be * ga
            O_formula = be * sp.diff(al, rho) - al * sp.diff(be, rho) + g1 * be**2 - cterm
            ok_a &= zero(O - O_formula)
            ok_b &= zero(reduce_ode(O, eps, flip=(MUT == "frob-ode")) - O_closed(s, eps, g1, W, 1 / rho, c))
    check("frobenius", ok_a and ok_b, f"the coordinate coefficient of k~ ^ dk~ and m~ ^ dm~ equals beta alpha' - alpha beta' + "
          f"g' beta^2 - c beta gamma for k~ = alpha th0 + beta th1 + gamma th2, with no ODE: {ok_a}; on the ODE it equals "
          f"the closed forms O_k, O_m for eps = +-1: {ok_b}")


# ---------------------------------------------------------------------------------------------------------------
# C3 sum
# ---------------------------------------------------------------------------------------------------------------
def c3_sum():
    ok_s = True
    for eps in (1, -1):
        Ok, Om = O_of(kcov(1, eps)), O_of(kcov(-1, eps))
        expect = 2 * g1 * (g1**2 + (2 if MUT == "sum-factor" else 4) * c**2)
        ok_s &= zero(Ok + Om - expect)
    disc = sp.discriminant(Gs**2 + 4 * c**2, Gs)
    ok_d = zero(disc + 16 * c**2)
    roots = {cc: sp.real_roots(sp.Poly(Gs * (Gs**2 + 4 * cc**2), Gs)) for cc in (R(1, 2), -3, 7)}
    ok_r = all(rl == [0] for rl in roots.values())
    check("sum", ok_s and ok_d and ok_r, f"O_k + O_m = 2 g'(g'^2 + 4 c^2) with no ODE, for eps = +-1: {ok_s}; the discriminant of "
          f"G^2 + 4 c^2 in G is -16 c^2 < 0 for c != 0: {ok_d}; the only real root of 2 G (G^2 + 4 c^2) is G = 0 at "
          f"c = 1/2, -3, 7: {ok_r}")


# ---------------------------------------------------------------------------------------------------------------
# C4 isolated
# ---------------------------------------------------------------------------------------------------------------
def c4_isolated():
    ok_a = True
    samples = (0,) if MUT == "isolated-c0" else (R(1, 2), -3, 7, 1)
    for eps in (1, -1):
        Gp0 = sp.simplify(g2rule(eps).subs(g1, 0))  # G' = g'' at G = 0, from the ODE
        for cc in samples:
            val = Gp0.subs(c, cc)
            ok_a &= (val == 2 * eps * cc**2) and val != 0
    # planted: g = 1 (constant) in g'' + g'/rho - eps (g'^2 + 2 c^2): a member iff c = 0
    one = sp.Integer(1)
    resid = lambda cc, eps: sp.diff(one, rho, 2) + sp.diff(one, rho) / rho - eps * (sp.diff(one, rho) ** 2 + 2 * cc**2)  # noqa: E731
    ok_p = all(resid(0, e) == 0 for e in (1, -1)) and all(resid(cc, e) != 0 for cc in (R(1, 2), -3) for e in (1, -1))
    check("isolated", ok_a and ok_p, f"g'' = 2 eps c^2 != 0 at g' = 0 for c in {samples}, eps = +-1: {ok_a}; planted: g = const "
          f"solves the ODE at c = 0 and at no c in (1/2, -3): {ok_p}")


# ---------------------------------------------------------------------------------------------------------------
# C5 each
# ---------------------------------------------------------------------------------------------------------------
P_PRINTED = (Gs**12 - 21 * c**2 * Gs**10 - 78 * c**4 * Gs**8 + 720 * c**6 * Gs**6 + 5936 * c**8 * Gs**4
             + 13056 * c**10 * Gs**2 + 8704 * c**12)


def lie_N(eps, s):
    """H = (W/G) O_{k,m} on W^2 = G^2 + 2 c^2, the flow (G, r)' = (eps W^2 - G r, -r^2), W' = G G'/W. Returns the
    elimination of r through H = 0: the numerator of L_V H = p0 + p1 W and N = p0^2 - p1^2 W^2."""
    W2 = Gs**2 + (1 if MUT == "each-w2" else 2) * c**2
    H = (Gs**2 + 4 * c**2) * Ws + s * (5 * c * eps * W2 - 2 * c * Gs * rs)
    Gp = eps * W2 - (0 if MUT == "each-friction" else Gs * rs)
    Hp = sp.diff(H, Gs) * Gp + sp.diff(H, Ws) * (Gs * Gp / Ws) + sp.diff(H, rs) * (-(rs**2))
    F = sp.solve(H, rs)[0]
    num = sp.expand(sp.numer(sp.together(Hp.subs(rs, F))))
    p0, p1 = 0, 0
    for (e,), co in sp.Poly(num, Ws).terms():
        if e % 2 == 0:
            p0 += co * W2 ** (e // 2)
        else:
            p1 += co * W2 ** (e // 2)
    p0, p1 = sp.expand(p0), sp.expand(p1)
    return H, num, p0, p1, sp.expand(p0**2 - p1**2 * W2), W2


def c5_each():
    ok_N, ok_link, ok_num = True, True, True
    for eps in (1, -1):
        for s in (1, -1):
            H, num, p0, p1, N, W2 = lie_N(eps, s)
            ok_N &= sp.expand(N + (Gs**2 + 2 * c**2) * P_PRINTED) == 0
            link = sp.expand(O_closed(s, eps, Gs, Ws, rs, c) * Ws / Gs - H)
            ok_link &= sp.rem(sp.Poly(sp.expand(link), Ws), sp.Poly(Ws**2 - (Gs**2 + 2 * c**2), Ws)).is_zero
            pt = {Gs: R(3, 2), c: 1}
            Wv = sp.sqrt(R(9, 4) + 2)
            ok_num &= abs(sp.N(num.subs(pt).subs(Ws, Wv) - (p0 + p1 * Ws).subs(pt).subs(Ws, Wv), 60)) < sp.Float("1e-40")
    lc = sp.Poly(P_PRINTED, Gs).LC()
    rts = sorted(float(q) for q in sp.Poly(P_PRINTED.subs(c, 1), Gs).nroots(n=15) if abs(sp.im(q)) < 1e-12)
    print(f"  observation (not a criterion): real roots of the degree-12 factor at c = 1: {[round(q, 4) for q in rts]}")
    check("each", ok_N and ok_link and lc == 1 and ok_num,
          f"H = (W/G) O is the closed form on W^2 = G^2 + 2 c^2: {ok_link}; eliminating r through H = 0 gives "
          f"N = -(G^2 + 2 c^2) P(G, c), P monic of degree 12 (printed), the same for k, m and eps = +-1: {ok_N} (leading "
          f"coefficient {lc}); the reduction p0 + p1 W equals the numerator at (G, c) = (3/2, 1) to 1e-40: {ok_num}")


# ---------------------------------------------------------------------------------------------------------------
# C6 points
# ---------------------------------------------------------------------------------------------------------------
G0s, G1s = sp.symbols("G0 G1s")
JETS = [(2, 1, R(5, 7)), (6, 17, R(5, 7)), (0, R(1, 2), R(5, 7)),  # #80's (c, g'(rho0)) at rho0 = 5/7
        (-2, 1, R(4, 39)), (2, 1, R(4, 39)),  # B_k = 0 (eps = +1, c = -2) and B_m = 0 (eps = +1, c = 2)
        (1, 1, R(5, 7)), (3, R(7, 2), R(2, 3)), (-6, 7, R(3, 11)),  # irrational W, W = 11/2, W = 11
        (1, 0, R(5, 7)), (-3, 0, R(1, 2))]  # g' = 0


def c6_points():
    rows, ok_O, ok_cl = [], True, True
    for cc, gv, rho0 in JETS:
        for eps in (1, -1):
            rule2 = eps * (g1**2 + 2 * cc**2) - g1 / rho
            rule3 = sp.diff(rule2, rho).subs(sp.Derivative(g, (rho, 2)), rule2)
            pt = {rho: rho0, x: R(1, 3), t: R(1, 10)}

            def simp(e, rule2=rule2, rule3=rule3, gv=gv, pt=pt):
                e = e.subs(sp.Derivative(g, (rho, 3)), rule3).subs(sp.Derivative(g, (rho, 2)), rule2)
                e = e.subs(sp.Derivative(g, rho), G1s).subs(g, G0s)
                return sp.simplify(e.subs({G1s: gv, G0s: R(1, 5)}).subs(pt))
            kk, mm = kcov(1, eps, cc=cc), kcov(-1, eps, cc=cc)
            if MUT == "points-closed":
                a, b = sp.Matrix([0, 0, 1, 0]), sp.Matrix([1, 0, 0, 0])
                kk, mm = (a + b) * 2, (a - b) * 2
            else:
                a, b = (kk + mm) / (2 * sp.sqrt(2)), (kk - mm) / (2 * sp.sqrt(2))
            Ok, Om = simp(O_of(kk)), simp(O_of(mm))
            Wv = sp.sqrt(gv**2 + 2 * cc**2)
            pOk, pOm = (sp.simplify(O_closed(s_, eps, gv, Wv, 1 / rho0, cc)) for s_ in (1, -1))
            ok_O &= zero(Ok - pOk) and zero(Om - pOm)
            ok, om, why = ts.closure(a, b, X4, simp)
            fails = (not ok) and why.startswith("no chi")
            predicted = not (zero(pOk) and zero(pOm))
            ok_cl &= fails == predicted
            rows.append((cc, gv, str(rho0), eps, sp.nsimplify(sp.simplify(pOk)), sp.nsimplify(sp.simplify(pOm)),
                         "fails" if fails else f"pointwise: {why}"))
    print("  (c, g', rho0, eps, O_k, O_m, two_scalars.closure):")
    for rw in rows:
        print("   ", rw)
    check("points", ok_O and ok_cl, f"at {len(rows)} (jet, eps) points the coordinate O_k, O_m equal the closed forms exactly: "
          f"{ok_O}; two_scalars.closure fails with 'no chi' exactly where (O_k, O_m) != (0, 0): {ok_cl}")


# ---------------------------------------------------------------------------------------------------------------
# C7 bessel
# ---------------------------------------------------------------------------------------------------------------
def c7_bessel():
    mp.mp.dps = 80
    random.seed(162)
    worst_b, worst_s, npt = mp.mpf(0), mp.mpf(0), 0
    obs = []
    for cc, eps in ((R(1, 2), 1), (R(3, 2), -1), (R(-2), 1), (R(1), -1), (R(5, 2), 1)):
        cv = mp.mpf(cc.p) / cc.q
        al, be = R(random.randint(1, 9), 3), R(random.randint(-9, -1), 5)
        kap = abs(cc) if MUT == "bessel-freq" else sp.sqrt(2) * abs(cc)
        psi = al * sp.besselj(0, kap * rho) + be * sp.bessely(0, kap * rho)
        gg = -eps * sp.log(psi)
        G = sp.diff(gg, rho)
        W = sp.sqrt(G**2 + 2 * cc**2)
        Ok, Om = (O_of(kcov(s_, eps, gfun=gg, cc=cc)) for s_ in (1, -1))
        fGWp = sp.lambdify(rho, (G, W, psi), modules="mpmath")
        for x0 in (R(1, 3), R(-5, 7)):
            fOk = sp.lambdify(rho, Ok.subs(x, x0), modules="mpmath")
            fOm = sp.lambdify(rho, Om.subs(x, x0), modules="mpmath")
            for rv in (R(3, 7), R(5, 7), R(9, 7), R(13, 7), R(17, 7), R(23, 7)):
                rr = mp.mpf(rv.p) / rv.q
                Gv, Wv, pv = fGWp(rr)
                if abs(pv) < mp.mpf(1) / 1000:
                    continue
                ok_, om_ = fOk(rr), fOm(rr)
                worst_b = max(worst_b, abs(ok_ - O_closed(1, eps, Gv, Wv, 1 / rr, cv)),
                              abs(om_ - O_closed(-1, eps, Gv, Wv, 1 / rr, cv)))
                worst_s = max(worst_s, abs(ok_ + om_ - 2 * Gv * (Gv**2 + 4 * cv**2)))
                npt += 1
        # observation: sign changes of G, O_k, O_m on a grid where |psi| > 1e-2
        vals = []
        for i in range(260):
            rr = mp.mpf(3) / 10 + mp.mpf(i) / 30
            Gv, Wv, pv = fGWp(rr)
            vals.append(None if abs(pv) < mp.mpf(1) / 100 else
                        (Gv, O_closed(1, eps, Gv, Wv, 1 / rr, cv), O_closed(-1, eps, Gv, Wv, 1 / rr, cv)))
        changes = [sum(1 for u_, v_ in zip(vals, vals[1:]) if u_ is not None and v_ is not None and u_[i] * v_[i] < 0)
                   for i in range(3)]
        obs.append((str(cc), eps, str(al), str(be), changes))
    print("  observation (not a criterion), (c, eps, alpha, beta, sign changes of g', O_k, O_m on 260 grid points of "
          "0.3 < rho < 8.95 where |psi| > 1e-2):")
    for o in obs:
        print("   ", o)
    ok = worst_b < mp.mpf(10) ** -50 and worst_s < mp.mpf(10) ** -50 and npt >= 40
    check("bessel", ok, f"{npt} (member, x, rho) points at 80 digits: max |direct k ^ dk - closed form| = {mp.nstr(worst_b, 3)} "
          f"(threshold 1e-50), max |O_k + O_m - 2 g'(g'^2 + 4 c^2)| = {mp.nstr(worst_s, 3)}")


# ---------------------------------------------------------------------------------------------------------------
# C8 planted
# ---------------------------------------------------------------------------------------------------------------
def planted_row(a, b):
    ok, om, why = ts.closure(a, b, X4, sp.simplify)
    return sp.simplify(O_of(a + b)), sp.simplify(O_of(a - b)), ok, why, om


def c8_planted():
    chi = c * x + g
    e1 = sp.Matrix([1, 0, 0, 0])
    # P1: closed pair (d rho, dt)
    a1, b1 = sp.Matrix([0, 0, 1, 0]), e1
    # P2: boosted closed pair k = e^-chi d(t + x), m = e^chi d(t - x + rho^2/2)
    k2 = sp.exp(-chi) * sp.Matrix([1, 1, 0, 0])
    m2 = sp.exp(chi) * (sp.Matrix([1, -1, x, 0]) if MUT == "plant-pass" else sp.Matrix([1, -1, rho, 0]))
    a2, b2 = (k2 + m2) / 2, (k2 - m2) / 2
    # P3: non-integrable pair (dt, x d rho), or the closed (dt, rho d rho) under the control
    a3, b3 = e1, sp.Matrix([0, 0, rho if MUT == "plant-fail" else x, 0])
    # P4: k = dt and m = (1 + t rho) d rho are each integrable but e^chi k, e^-chi m cannot both be closed
    k4, m4 = e1, sp.Matrix([0, 0, 1 + t * rho, 0])
    a4, b4 = (k4 + m4) / 2, (k4 - m4) / 2
    r1, r2, r3, r4 = (planted_row(a, b) for a, b in ((a1, b1), (a2, b2), (a3, b3), (a4, b4)))
    ok1 = r1[0] == 0 and r1[1] == 0 and r1[2] and list(r1[4]) == [0, 0, 0, 0]
    ok2 = r2[0] == 0 and r2[1] == 0 and r2[2] and all(zero(u_ - v_) for u_, v_ in zip(r2[4], [0, c, g1, 0]))
    ok3 = r3[0] == 1 and r3[1] == -1 and (not r3[2]) and r3[3].startswith("no chi")
    ok4 = r4[0] == 0 and r4[1] == 0 and (not r4[2]) and r4[3].startswith("omega = d chi needed")
    check("planted", ok1 and ok2 and ok3 and ok4,
          f"P1 closed (d rho, dt): O = (0, 0), closure ok, omega = 0: {ok1}; P2 boosted by chi = c x + g: O = (0, 0), closure ok, "
          f"omega = d chi = (0, c, g', 0): {ok2}; P3 (dt, x d rho): O = (1, -1), closure fails 'no chi': {ok3}; P4 (dt, (1 + t rho) "
          f"d rho): O = (0, 0) and closure fails at d omega != 0, so the Frobenius test alone does not decide closure: {ok4}")


# ---------------------------------------------------------------------------------------------------------------
# C9 c = 0
# ---------------------------------------------------------------------------------------------------------------
def c9_c0():
    ok_O, ok_nu = True, True
    for eps in (1, -1):
        for sgn in (1, -1):  # sign of g'; W = |g'| = sgn g'
            for s in (1, -1):
                O = O_of(kcov(s, eps, cc=0, W=sgn * g1))
                ok_O &= zero(O - (-g1**3 if MUT == "c0-sign" else g1**3))
        nu = sp.Matrix([1 + eps * g, eps, 0, 0])
        dr = sp.Matrix([0, 0, 1, 0])
        ok_nu &= zero(O_of(dr + nu) - g1) and zero(O_of(dr - nu) - g1)
    check("c0", ok_O and ok_nu, f"c = 0 (u = g(rho), #65's b(rho)): O_k = O_m = g'^3 for both signs of g' and of eps: {ok_O}; "
          f"nu_+- ^ d nu_+- = g' dt ^ dx ^ d rho = b' dx ^ d rho ^ dt for nu_+- = d rho +- ((1 + eps b) dt + eps dx): {ok_nu}")


# ---------------------------------------------------------------------------------------------------------------
# C10 planar
# ---------------------------------------------------------------------------------------------------------------
def c10_planar():
    d = sp.Symbol("d", real=True)
    u = c * x + d
    th0, th1, th2 = coframe(u)
    Rc = ricci(shift_metric(u), X4)
    ok_R = all(zero(e) for e in Rc - (c**2 * th1 * th1.T - c**2 * th0 * th0.T))
    if MUT == "planar-closed":
        a, b = c * sp.Matrix([0, 1, 0, 0]), c * sp.Matrix([1, 0, 0, 0])
    else:
        a, b = c * th1, c * th0
    Ok, Om = sp.simplify(O_of(a + b)), sp.simplify(O_of(a - b))
    ok, om, why = ts.closure(a, b, X4, sp.simplify)
    dom = ts.dform(om, X4)
    ok_c = Ok == 0 and Om == 0 and (not ok) and why.startswith("omega = d chi needed") and zero(dom[0][1] - c**2)
    check("planar", ok_R and ok_c, f"u = c x + d (g' = 0, outside the ODE): R = c^2 (th1 th1 - th0 th0) exactly: {ok_R}; "
          f"O_k = O_m = 0, two_scalars.closure fails at d omega != 0 with (d omega)_tx = c^2: {ok_c}")


CHECKS = [("frame", c1_frame), ("frobenius", c2_frobenius), ("sum", c3_sum), ("isolated", c4_isolated), ("each", c5_each),
          ("points", c6_points), ("bessel", c7_bessel), ("planted", c8_planted), ("c0", c9_c0), ("planar", c10_planar)]


def main():
    t0 = time.time()
    if MUT is not None and MUT not in MUTATIONS:
        print(f"unknown mutation {MUT}")
        sys.exit(2)
    for name, fn in CHECKS:
        if wanted(name):
            t1 = time.time()
            fn()
            print(f"  ({time.time() - t1:.0f} s)", flush=True)
    failed = [n for n, ok in RESULTS if not ok]
    print(f"\nchecks: {len(RESULTS) - len(failed)} passed, {len(failed)} failed {failed}. {time.time() - t0:.0f} s")
    if MUT is not None:
        want = MUTATIONS[MUT][0]
        print(f"control {MUT}: target checks {want}, fired {[n for n in want if n in failed]}")
    return failed


if __name__ == "__main__":
    failed = main()
    controls_ok = True
    if MUT is None and "no-controls" not in OPTS and ONLY is None:
        print("\n== Controls: each mutation must make its target checks fail, and no others ==", flush=True)
        for name, (want, what) in MUTATIONS.items():
            p = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True)
            if "logdir" in OPTS:
                with open(os.path.join(OPTS["logdir"], f"mut_{name}.log"), "w") as fh:
                    fh.write(p.stdout + p.stderr)
            fails = [ln.split(":")[0][7:] for ln in p.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = all(w in fails for w in want) and p.returncode == 1
            extra = [f for f in fails if f not in want and f not in failed]
            controls_ok &= fired and not extra
            print(f"[{'FIRED' if fired and not extra else 'SILENT' if not fired else 'EXTRA'}] --mutate={name} ({what}): "
                  f"exit {p.returncode}, failed {fails}", flush=True)
            if p.returncode not in (0, 1):
                print(p.stderr[-2000:])
        print(f"controls: {'all fired on target only' if controls_ok else 'SOME DID NOT FIRE ON TARGET ONLY'}")
    sys.exit(0 if not failed and controls_ok else 1)

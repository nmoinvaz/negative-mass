#!/usr/bin/env python3
"""Attraction to the conical family's line of equilibria (issue #168, continues conical_limit.md / #159).

N(F, F', F''; eta) = 0 (det_classify.md section 6) in the variables G = F, p = eta F', q = eta^2 F'' + eta F', nu = 1/eta, sigma = nu^2,
t = ln eta. The note conical_attraction.md proves, by a Gronwall bootstrap and not by a cited theorem,
  A  every solution with data near the line {(G, 4G, 0)} converges to a point of it, with explicit constants,
  B  the formal series of conical_limit.md section 3 is its asymptotic expansion to every order,
  O  a real C^2 solution of N = 0 on a half line with F -> F_inf is eventually constant or in this family (no oscillation).
This script checks every constant of those proofs. Planted cases run first, through the same functions.
  P1 planted_hyp       hyp_line on a linear field with eigenvalues (0,-2,-2) and a Jordan block, on (0,-2,2), on (0,0,-2)
  P2 planted_gronwall  the bootstrap chain on zeta' = -2 zeta + ubar exp(-2t), exact solution known
  P3 planted_shape     the root classifier of the cubic P on a cubic with known roots, and on one with M < 0
  P4 planted_local     the local lemma on N_osc = nu^2 (q - p) + G (F'' + F = 0): not applicable
  P5 planted_series    the stage solver on q + p = 0 and q + 4 p = 8 nu^4
  M1 n_identity        N from D = theta-form = c3 q^3 + ... = p^6 P(q/p^2); symmetry; R = sigma^2 r^6 P
  M2 reduction         y' equation, e-hat expansion, f(zeta, sigma, g2), Jacobian eigenvalues, Jordan block
  M3 hyp_line_real     hyp_line on the real reduced system at five values of G
  M4 shape             interval constants of the shape lemma L1
  M5 ehat              the bound |e-hat| <= u (alpha y^2 + beta y)
  M6 theorem_a         rational chain of the bootstrap constants (exact exp bound)
  M7 numeric_a         reduced system and direct N = 0 integrated; every inequality of Theorem A tested
  M8 series            exact solver: uniqueness, table of conical_limit.md, residual order
  M9 stage_symbol      indicial polynomial I(k) = -128 F^3 k (k - 4)
  M10 expansion_link   zeta~ from the series: sigma t coefficient -4 (F^2 + 7), slope -4/F in C
  M11 numeric_b        mpmath integration, fitted F_inf and C, normalized errors of the truncations
  M12 local_lemma      divisibility, root bound, N(G, 0, q) = c3 q^3
  M13 regime_b         |s| <= |r|/10 for |y| <= 1/50
  M14 small_branch     (ln|omega|)' <= -1 on small roots, sign cases
  M15 huge_branch      sign and size of y' on the huge root
  M16 classify_numeric branch trichotomy of the cubic N = 0, integrated
Usage: conical_attraction.py [--no-controls] [--mutate=NAME] [--only=check1,check2] [--skip-num] [--logdir=DIR].
Exit 1 if a check fails or a control does not fire on its target only."""
import itertools
import math
import os
import random
import subprocess
import sys
import time
from fractions import Fraction as Fr

import mpmath as mp
import sympy as sp

# name: (checks that must fail, what the mutation does). Targets are declared after the development run (README amendments).
MUTATIONS = {
    "d-sign": (["n_identity"], "Y built as u u_xrho/2 - u_x u_rho"),
    "P-drop": (["n_identity"], "P without the (b - sigma) omega term of M"),
    "ehat-sign": (["reduction"], "y' with +sigma y^2"),
    "f-sigma": (["reduction"], "expects f_sigma = -4 (g2 + 6)"),
    "eig-plus": (["reduction", "hyp_line_real"], "expects the eigenvalue +2 in place of -2"),
    "shape-3": (["shape"], "tests P(-1/(2a)) < 0 in place of P(-1/(4a)) < 0"),
    "shape-om": (["shape"], "omega range [-300, 300] in place of [-100, 100]"),
    "ehat-alpha": (["ehat"], "alpha = 2 in place of 51/25"),
    "A-u0": (["theorem_a"], "u_0 = 1/100 in place of 1/1000"),
    "A-cf": (["theorem_a"], "c_f = 40 in place of 50"),
    "A-kappa": (["planted_gronwall"], "stage-1 rate kappa = 2 in place of 2 - a/2"),
    "num-rate": (["numeric_a"], "rate bound with exp(-3t) in place of exp(-2t)"),
    "series-a2": (["series"], "a_2 = -(11/5) F in place of -2 F"),
    "series-8": (["series"], "expects a_{4,1} = F (F^2 + 8)"),
    "stage-I": (["stage_symbol"], "expects I(k) = k (k - 5)"),
    "link-slope": (["expansion_link"], "expects the sigma t coefficient -4 (F^2 + 6)"),
    "numB-coef": (["numeric_b"], "series with a_{6,0} times 1.01"),
    "local-c2": (["local_lemma"], "adds q^2 to N (c2 no longer divisible by p)"),
    "regimeb-y": (["regime_b"], "|y| <= 1/5 in place of 1/50"),
    "small-sigma": (["small_branch"], "claims (ln|omega|)' <= -199/100"),
    "huge-thr": (["huge_branch"], "claims y' >= (3/2) y on [1/50, 1]"),
    "class-sign": (["classify_numeric"], "predicts convergence for the huge root with p < 0"),
    "plant-hyp": (["planted_hyp"], "planted attracting case uses the eigenvalue +2"),
    "plant-gron": (["planted_gronwall"], "c_f halved in the chain"),
    "plant-shape": (["planted_shape"], "planted cubic with the huge root at -1/(8a)"),
    "plant-local": (["planted_local"], "N_osc replaced by q^3 nu^2 (-2G) (the lemma applies)"),
    "plant-series": (["planted_series"], "forcing 0 in place of 8 nu^4"),
}
NUMERIC = ("numeric_a", "numeric_b", "classify_numeric")

RESULTS = []
OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate")
ONLY = OPTS["only"].split(",") if "only" in OPTS else None

G, p, q, nu = sp.symbols("G p q nu")
sig, r_, s_, y_, e_, g2_, w_, om_, x_ = sp.symbols("sigma r s y e g2 w omega x")
F_, C_ = sp.symbols("F C")


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def wanted(name):
    if "skip-num" in OPTS and name in NUMERIC and not (MUT and name in MUTATIONS[MUT][0]):
        return False
    return ONLY is None or name in ONLY


# ---------------------------------------------------------------------------------------------------------------
# rational interval arithmetic (exact, Fractions): the trusted code of the computer-assisted constants
# ---------------------------------------------------------------------------------------------------------------
class Iv:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = Fr(lo)
        self.hi = Fr(lo if hi is None else hi)
        assert self.lo <= self.hi

    @staticmethod
    def c(x):
        return x if isinstance(x, Iv) else Iv(x)

    def __add__(a, b):
        b = Iv.c(b)
        return Iv(a.lo + b.lo, a.hi + b.hi)

    __radd__ = __add__

    def __neg__(a):
        return Iv(-a.hi, -a.lo)

    def __sub__(a, b):
        b = Iv.c(b)
        return Iv(a.lo - b.hi, a.hi - b.lo)

    def __rsub__(a, b):
        return Iv.c(b) - a

    def __mul__(a, b):
        b = Iv.c(b)
        v = [a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi]
        return Iv(min(v), max(v))

    __rmul__ = __mul__

    def __pow__(a, n):
        assert isinstance(n, int) and n >= 0
        if n == 0:
            return Iv(1)
        if n % 2 == 1 or a.lo >= 0:
            return Iv(a.lo ** n, a.hi ** n)
        if a.hi <= 0:
            return Iv(a.hi ** n, a.lo ** n)
        return Iv(Fr(0), max(a.lo ** n, a.hi ** n))

    def mag(a):
        return max(abs(a.lo), abs(a.hi))

    def __repr__(a):
        return f"[{float(a.lo):.6g}, {float(a.hi):.6g}]"


def exp_upper(x, n=14):
    """rigorous upper bound of e^x for rational 0 <= x < 1."""
    x = Fr(x)
    tot, term = Fr(0), Fr(1)
    for k in range(n + 1):
        tot += term
        term = term * x / (k + 1)
    return tot + term / (1 - x / (n + 2))


def param_boxes(U):
    """sigma, b in [0, U] and a in [0, U/2]: a box that contains every (sigma, b, a) with sigma + b <= U, a^2 = sigma b."""
    return Iv(0, U), Iv(0, U), Iv(0, U / 2)


def cubic_params(sg, b, om):
    """kappa, mu, c, M of the cubic P; works for Fractions, Iv and sympy."""
    kap = 1 - b + sg
    mu = 1 + 3 * sg - b + 2 * sg ** 2 - 6 * b * sg
    c = 1 + 2 * sg
    M = mu + (b - sg) * om
    return kap, mu, c, M


def Pcubic(w, a, kap, M, c, om):
    return -2 * a * kap * w ** 3 - M * w ** 2 + 4 * a * c * w + c ** 2 - om ** 2


# ---------------------------------------------------------------------------------------------------------------
# N: typed theta-form, and N rebuilt from D (retyped from det_classify.md Setup, as in conical_limit.py n_build)
# ---------------------------------------------------------------------------------------------------------------
A0 = p ** 2 * (p ** 2 - q) * (p ** 2 + q)
A1 = G ** 2 * p ** 2 * q ** 2 + 4 * G * p ** 4 * q - 2 * G * q ** 3 + 4 * p ** 6 - 3 * p ** 2 * q ** 2
A2 = 2 * (G * (q - p) + p ** 2) * (G ** 2 * q ** 2 + 2 * G * p ** 3 + 2 * G * p ** 2 * q + 2 * p ** 4 - q ** 2)
N_TH = sp.expand(A0 + nu ** 2 * A1 + nu ** 4 * A2)


def D_of(P, Q, Z, Y, h):
    return h * (h ** 2 - P ** 2 - Q ** 2 + Y ** 2 + Z ** 2) - P * (Y ** 2 - Z ** 2) + 2 * Q * Y * Z


def build_N(ysign=1):
    """8 rho^6 D for u = F(x/rho), as a polynomial in (F0, F1, F2, eta)."""
    x, rho = sp.symbols("x rho", positive=True)
    eta = sp.Symbol("eta", positive=True)
    F0, F1, F2 = sp.symbols("F0 F1 F2")
    F = sp.Function("F")
    u = F(x / rho)
    ux, ur = sp.diff(u, x), sp.diff(u, rho)
    uxx, uxr, urr = sp.diff(u, x, 2), sp.diff(u, x, rho), sp.diff(u, rho, 2)
    P = -(u * uxx + ux ** 2)
    Q = -(urr + ur / rho) / 2
    Z = uxr / 2
    Y = u * uxr / 2 + ysign * ux * ur
    h = ur ** 2 / 2
    D = D_of(P, Q, Z, Y, h).doit().subs(x, eta * rho).doit()
    D = D.replace(lambda e: isinstance(e, sp.Subs), lambda e: e.doit())
    d2, d1 = sp.Derivative(F(eta), (eta, 2)), sp.Derivative(F(eta), eta)
    D = D.subs(d2, F2).subs(d1, F1).subs(F(eta), F0)
    Neta = sp.expand(sp.simplify(8 * rho ** 6 * D))
    return sp.expand(Neta.subs({F0: G, F1: p / eta, F2: (q - p) / eta ** 2}).subs(eta, 1 / nu))


# ---------------------------------------------------------------------------------------------------------------
# truncated series in (nu^n T^j), exact over Fractions, and the stage solver
# ---------------------------------------------------------------------------------------------------------------
NM = 17


class Ser:
    __slots__ = ("d",)

    def __init__(self, d=None):
        self.d = {k: v for k, v in (d or {}).items() if v != 0}

    @staticmethod
    def const(c):
        return Ser({(0, 0): c})

    def __add__(a, b):
        if not isinstance(b, Ser):
            b = Ser.const(b)
        d = dict(a.d)
        for k, v in b.d.items():
            d[k] = d.get(k, 0) + v
        return Ser(d)

    __radd__ = __add__

    def __neg__(a):
        return Ser({k: -v for k, v in a.d.items()})

    def __sub__(a, b):
        if not isinstance(b, Ser):
            b = Ser.const(b)
        return a + (-b)

    def __mul__(a, b):
        if not isinstance(b, Ser):
            return Ser({k: v * b for k, v in a.d.items()})
        d = {}
        for (n1, j1), v1 in a.d.items():
            for (n2, j2), v2 in b.d.items():
                n = n1 + n2
                if n > NM:
                    continue
                k = (n, j1 + j2)
                d[k] = d.get(k, 0) + v1 * v2
        return Ser(d)

    __rmul__ = __mul__

    def theta(a):
        d = {}
        for (n, j), v in a.d.items():
            d[(n, j)] = d.get((n, j), 0) - n * v
            if j > 0:
                d[(n, j - 1)] = d.get((n, j - 1), 0) + j * v
        return Ser(d)

    def shift(a, m):
        return Ser({(n + m, j): v for (n, j), v in a.d.items() if n + m <= NM})


def N_series(Gs):
    """N = A0 + nu^2 A1 + nu^4 A2 on a series G = F + g (the typed theta-form)."""
    ps = Gs.theta()
    qs = ps.theta()
    p2 = ps * ps
    a0 = p2 * (p2 - qs) * (p2 + qs)
    a1 = Gs * Gs * p2 * qs * qs + 4 * Gs * p2 * p2 * qs - 2 * Gs * qs * qs * qs + 4 * p2 * p2 * p2 - 3 * p2 * qs * qs
    a2 = 2 * (Gs * (qs - ps) + p2) * (Gs * Gs * qs * qs + 2 * Gs * p2 * ps + 2 * Gs * p2 * qs + 2 * p2 * p2 - qs * qs)
    return a0 + a1.shift(2) + a2.shift(4)


def rank_and_solve(Amat, bvec, ncols):
    """exact Gauss-Jordan over Fractions: (rank, particular solution or None, null space basis)."""
    rows = [list(r) + [b] for r, b in zip(Amat, bvec)]
    piv, r = [], 0
    for c in range(ncols):
        pr = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        pv = rows[r][c]
        rows[r] = [x / pv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
    if any(all(x == 0 for x in row[:ncols]) and row[ncols] != 0 for row in rows):
        return len(piv), None, None
    sol = [Fr(0)] * ncols
    for i, c in enumerate(piv):
        sol[c] = rows[i][ncols]
    null = []
    for fcol in [c for c in range(ncols) if c not in piv]:
        v = [Fr(0)] * ncols
        v[fcol] = Fr(1)
        for i, c in enumerate(piv):
            v[c] = -rows[i][fcol]
        null.append(v)
    return len(piv), sol, null


def stage_solve(Nfun, g0, ks, shift, Jof, free=None):
    """stage-by-stage formal solution g = sum a_{n,j} nu^n T^j of Nfun(g) = 0. At stage k the unknowns a_{k,0..J} enter linearly
    in the coefficient of nu^(k+shift); lower orders must vanish. free = (k, Cval): the free parameter a_{k,0} is set to Cval.
    returns g and a log of (k, rank, ncols, null space, status)."""
    g = g0
    log = []
    for k in ks:
        J = Jof(k)
        ncols = J + 1

        def resid(vec):
            gk = g + Ser({(k, j): vec[j] for j in range(ncols) if vec[j] != 0})
            Nv = Nfun(gk)
            low = [(n, j) for (n, j), v in Nv.d.items() if n < k + shift and v != 0]
            return [Nv.d.get((k + shift, i), Fr(0)) for i in range(ncols + 6)], low

        r0, low = resid([Fr(0)] * ncols)
        cols = []
        for j in range(ncols):
            e = [Fr(0)] * ncols
            e[j] = Fr(1)
            rj, _ = resid(e)
            cols.append([x - y for x, y in zip(rj, r0)])
        Amat = [[cols[j][i] for j in range(ncols)] for i in range(ncols + 6)]
        rk, sol, null = rank_and_solve(Amat, [-x for x in r0], ncols)
        if sol is None:
            log.append((k, rk, ncols, [], "obstruction", low))
            break
        if free is not None and k == free[0] and len(null) == 1 and null[0][0] != 0:
            sh = (Fr(free[1]) - sol[0]) / null[0][0]
            sol = [s + sh * n for s, n in zip(sol, null[0])]
        g = g + Ser({(k, j): sol[j] for j in range(ncols) if sol[j] != 0})
        log.append((k, rk, ncols, null, "ok", low))
    return g, log


def solve_series(Fi, Cval, kmax=9, a2=None):
    Fi = Fr(Fi)
    g0 = Ser({(2, 0): -2 * Fi if a2 is None else Fr(a2)})
    return stage_solve(lambda gk: N_series(Ser.const(Fi) + gk), g0, range(3, kmax + 1), 6, lambda k: k // 2 + 2, free=(4, Cval))


# ---------------------------------------------------------------------------------------------------------------
# P1, the hypotheses of the normally attracting line of equilibria
# ---------------------------------------------------------------------------------------------------------------
def hyp_line(J, tangent):
    """Hypotheses of a normally attracting line of equilibria at a point of x' = f(x) with Jacobian J (analyticity is checked
    separately): the tangent of the line is in ker J, 0 is a simple eigenvalue (so ker J = the tangent line), every other
    eigenvalue has Re < 0. Returns (verdict, jordan block present, eigenvalues)."""
    n = J.shape[0]
    ev = J.eigenvals()
    ker = len(J.nullspace())
    tan_ok = (J * tangent).is_zero_matrix
    if ker > 1:
        return "kernel larger than the line", False, ev
    zero_simple = ev.get(0, 0) == 1
    others = [lam for lam in ev if lam != 0]
    neg = all(sp.re(lam) < 0 for lam in others)
    jordan = any(ev[lam] > 1 and (J - lam * sp.eye(n)).rank() > n - ev[lam] for lam in others)
    if tan_ok and zero_simple and neg:
        return "attracting", jordan, ev
    return "not attracting", jordan, ev


def c_planted_hyp():
    lam = 2 if MUT == "plant-hyp" else -2
    A = sp.Matrix([[0, 0, 0], [0, lam, 1], [0, 0, lam]])
    e1 = sp.Matrix([1, 0, 0])
    v1 = hyp_line(A, e1)
    v2 = hyp_line(sp.diag(0, -2, 2), e1)
    v3 = hyp_line(sp.diag(0, 0, -2), e1)
    t = sp.Symbol("t", positive=True)
    E = sp.simplify(sp.exp(sp.Matrix([[0, 0, 0], [0, -2, 1], [0, 0, -2]]) * t))
    closed = (E - sp.Matrix([[1, 0, 0], [0, sp.exp(-2 * t), t * sp.exp(-2 * t)], [0, 0, sp.exp(-2 * t)]])).applyfunc(sp.simplify) == sp.zeros(3, 3)
    grid = [Fr(k, 4) for k in range(0, 41)]
    x0 = (Fr(1, 2), Fr(-1), Fr(1))
    bound_ok = all(
        max(abs(float(math.exp(-2 * float(tt)) * (x0[1] + float(tt) * x0[2]))), abs(float(math.exp(-2 * float(tt)) * x0[2])))
        <= (1 + float(tt)) * math.exp(-2 * float(tt)) * max(abs(float(v)) for v in x0) + 1e-15 for tt in grid)
    ok = (v1[0] == "attracting" and v1[1] and closed and bound_ok and v2[0] == "not attracting" and v3[0] == "kernel larger than the line")
    check("planted_hyp", ok, f"Jordan field (0,-2,-2): {v1[0]}, Jordan block {v1[1]}, e^(At) closed form {closed}, |x_2,3| <= (1+t) e^(-2t) |x_0| "
          f"on a grid {bound_ok}; diag(0,-2,2): {v2[0]}; diag(0,0,-2): {v3[0]}")


# ---------------------------------------------------------------------------------------------------------------
# M1 n_identity
# ---------------------------------------------------------------------------------------------------------------
def kap_mu_c(sg, b):
    kap, mu, c, _ = cubic_params(sg, b, 0)
    return kap, mu, c


def c_n_identity():
    Nb = build_N(-1 if MUT == "d-sign" else 1)
    ok_a = sp.expand(Nb - N_TH) == 0
    sg = nu ** 2
    kap, mu, c = kap_mu_c(sg, G ** 2 * sg)
    c3 = -2 * G * sg * kap
    c2 = -p * (2 * G * (G ** 2 - 1) * sg ** 2 + p * mu)
    c1 = 4 * G * sg * c * p ** 4
    c0 = p ** 4 * (c ** 2 * p ** 2 - 4 * G ** 2 * sg ** 2)
    ok_b = sp.expand(N_TH - (c3 * q ** 3 + c2 * q ** 2 + c1 * q + c0)) == 0
    om = 2 * G * sg / p
    kap, mu, c, M = cubic_params(sg, G ** 2 * sg, om)
    if MUT == "P-drop":
        M = mu
    Pw = Pcubic(q / p ** 2, G * sg, kap, M, c, om)
    ok_c = sp.cancel(sp.together(N_TH - p ** 6 * Pw)) == 0
    ok_d = sp.expand(N_TH.subs({G: -G, p: -p, q: -q}, simultaneous=True) - N_TH) == 0
    Re = sp.expand(sp.cancel(sp.expand(N_TH.subs({p: nu ** 2 * r_, q: nu ** 2 * s_}, simultaneous=True)) / nu ** 8))
    R = sp.expand(Re.subs(nu, sp.sqrt(sig)))
    kap, mu, c, M = cubic_params(sig, G ** 2 * sig, 2 * G / r_)
    if MUT == "P-drop":
        M = mu
    target = sig ** 2 * r_ ** 6 * Pcubic(s_ / (sig * r_ ** 2), G * sig, kap, M, c, 2 * G / r_)
    ok_e = sp.cancel(sp.together(R - target)) == 0
    check("n_identity", ok_a and ok_b and ok_c and ok_d and ok_e,
          f"(a) N from D (retyped) = A0 + nu^2 A1 + nu^4 A2: {ok_a}; (b) N = c3 q^3 + c2 q^2 + c1 q + c0: {ok_b}; "
          f"(c) N = p^6 P(q/p^2), omega = 2 G sigma/p: {ok_c}; (d) N(-G,-p,-q) = N: {ok_d}; (e) R = N/sigma^4 = sigma^2 r^6 P(s/(sigma r^2)): {ok_e}")
    return R


# ---------------------------------------------------------------------------------------------------------------
# M2 reduction
# ---------------------------------------------------------------------------------------------------------------
def real_R():
    Re = sp.expand(sp.cancel(sp.expand(N_TH.subs({p: nu ** 2 * r_, q: nu ** 2 * s_}, simultaneous=True)) / nu ** 8))
    return sp.expand(Re.subs(nu, sp.sqrt(sig)))


def c_reduction():
    R = real_R()
    # (a) y' from G' = sigma r, r' = 2r + s, s = (G^2 e - r^2)/(2G)
    yy, ee = y_, e_
    rr = G * yy
    ss = (G ** 2 * ee - rr ** 2) / (2 * G)
    Gp = sig * rr
    rp = 2 * rr + ss
    yp = sp.simplify(rp / G - rr * Gp / G ** 2)
    sgn = 1 if MUT == "ehat-sign" else -1
    want = yy * (4 - yy) / 2 + sgn * sig * yy ** 2 + ee / 2
    ok_a = sp.simplify(yp - want) == 0 and sp.simplify(Gp / G - sig * yy) == 0
    # (b) rho = R / G^4 at r = G y, s = G (e - y^2)/2
    rho = sp.expand(R.subs({r_: G * y_, s_: G * (e_ - y_ ** 2) / 2}, simultaneous=True) / G ** 4)
    rho = sp.expand(rho.subs(G, sp.sqrt(g2_)))
    no_odd = not rho.has(sp.sqrt(g2_))
    rho_e = sp.diff(rho, e_).subs({e_: 0, sig: 0})
    rho_s = sp.diff(rho, sig).subs({e_: 0, sig: 0})
    e1 = sp.factor(sp.simplify(-rho_s / rho_e))
    e1_at4 = sp.factor(e1.subs(y_, 4))
    ok_b = no_odd and sp.simplify(rho_e + y_ ** 4 / 4) == 0 and sp.simplify(e1_at4 + 8 * (g2_ + 3)) == 0
    # (c) f = -zeta^2/2 + e/2 - sigma (4 + zeta)^2, e = e1(y) sigma + O(sigma^2)
    zeta = sp.Symbol("zeta")
    e_sig = sp.simplify(e1.subs(y_, 4 + zeta))
    f_sigma = sp.factor(sp.simplify(e_sig.subs(zeta, 0) / 2 - 16))
    f_zeta = sp.simplify(sp.diff(-zeta ** 2 / 2 + e_sig * sig / 2 - sig * (4 + zeta) ** 2, zeta).subs({zeta: 0, sig: 0}))
    rho_y0 = sp.simplify(sp.diff(rho.subs({e_: 0, sig: 0}), y_))
    c_want = -4 * (g2_ + (6 if MUT == "f-sigma" else 7))
    ok_c = sp.simplify(f_sigma - c_want) == 0 and f_zeta == 0 and rho_y0 == 0
    # (d) Jacobian of (G' = sigma r, r' = 2 r + S, sigma' = -2 sigma) at (G, 4G, 0)
    pt = {s_: -8 * G, r_: 4 * G, sig: 0}
    Rs = sp.diff(R, s_).subs(pt)
    Sr = sp.simplify(-sp.diff(R, r_).subs(pt) / Rs)
    SG = sp.simplify(-sp.diff(R, G).subs(pt) / Rs)
    Ss = sp.simplify(-sp.diff(R, sig).subs(pt) / Rs)
    J = sp.Matrix([[0, 0, 4 * G], [SG, 2 + Sr, Ss], [0, 0, -2]])
    eig = J.eigenvals()
    want_eig = {0: 1, -2: 2} if MUT != "eig-plus" else {0: 1, 2: 1, -2: 1}
    minor = sp.factor(sp.Matrix([[2, 4 * G], [8, Ss]]).det())
    ok_d = (sp.simplify(R.subs(pt)) == 0 and eig == want_eig and sp.simplify(2 * 4 * G - 8 * G) == 0
            and sp.simplify(minor + 8 * G * (G ** 2 + 7)) == 0 and Sr == -4 and SG == 8)
    check("reduction", ok_a and ok_b and ok_c and ok_d,
          f"(a) y' = y(4-y)/2 - sigma y^2 + e/2, (ln G)' = sigma y: {ok_a}; (b) d rho/d e = -y^4/4, e = {e1_at4} sigma + ... at y = 4: {ok_b}; "
          f"(c) f_sigma = {f_sigma}, f_zeta(0,0) = {f_zeta}, rho_y(e=0,sigma=0) = {rho_y0}: {ok_c}; (d) eigenvalues {eig} (S_r = {Sr}, S_G = {SG}, "
          f"S_sigma = {Ss}), minor of J + 2I = {minor} (never 0 for real G != 0, so a Jordan block): {ok_d}")
    return Ss


# ---------------------------------------------------------------------------------------------------------------
# M3 hyp_line_real
# ---------------------------------------------------------------------------------------------------------------
def c_hyp_line_real(Ss):
    Sig = sp.lambdify(G, Ss, "sympy")
    res = []
    for g in [-3, sp.Rational(-1, 2), 1, 2, 7]:
        lam = 2 if MUT == "eig-plus" else -2
        J = sp.Matrix([[0, 0, 4 * g], [8, lam, Sig(g)], [0, 0, -2]])
        v = hyp_line(J, sp.Matrix([1, 4, 0]))
        res.append((g, v[0], v[1]))
    ok = all(v == "attracting" and j for _, v, j in res)
    check("hyp_line_real", ok, "hyp_line (analytic, tangent (1,4,0) in ker J, 0 simple, other eigenvalues Re < 0) at G = -3, -1/2, 1, 2, 7: "
          + "; ".join(f"G={g}: {v}, Jordan {j}" for g, v, j in res))


# ---------------------------------------------------------------------------------------------------------------
# P3 / M4 shape lemma
# ---------------------------------------------------------------------------------------------------------------
def shape_hyps_from_coeffs(c3, c2, c1, c0):
    """hypotheses of the shape lemma for P(w) = c3 w^3 + c2 w^2 + c1 w + c0 written as -2 a kappa w^3 - M w^2 + ..., kappa = 1:
    returns {name: bool} (exact Fractions)."""
    a, kap, M = -Fr(c3) / 2, Fr(1), -Fr(c2)
    P = lambda w: c3 * w ** 3 + c2 * w ** 2 + c1 * w + c0
    return {"M>0": M > 0, "4M-3kappa>0": 4 * M - 3 * kap > 0, "P(-1/(4a))<0": P(Fr(-1) / (4 * a)) < 0,
            "P(2)<0": P(Fr(2)) < 0, "P(-2)<0": P(Fr(-2)) < 0}


def classify_roots(c3, c2, c1, c0):
    """real roots of the cubic split into the huge one (<= -1/(4a)) and the others."""
    a = -float(c3) / 2
    roots = sorted(z.real for z in mp.polyroots([float(c3), float(c2), float(c1), float(c0)], maxsteps=200, extraprec=200)
                   if abs(z.imag) < 1e-12)
    return [z for z in roots if z <= -1 / (4 * a)], [z for z in roots if z > -1 / (4 * a)]


def planted_cubic(a, h):
    """-2a (w - h)(w^2 - 1): roots h, 1, -1; M = -2 a h."""
    return -2 * a, 2 * a * h, 2 * a, -2 * a * h


def c_planted_shape():
    a = Fr(1, 100)
    h = Fr(-1, 8) / a if MUT == "plant-shape" else Fr(-1) / (2 * a)
    cf = planted_cubic(a, h)
    hyp = shape_hyps_from_coeffs(*cf)
    huge, small = classify_roots(*cf)
    classified = (len(huge) == 1 and abs(huge[0] - float(Fr(-1) / (2 * a))) < 1e-9 and len(small) == 2
                  and abs(small[0] + 1) < 1e-9 and abs(small[1] - 1) < 1e-9)
    hypb = shape_hyps_from_coeffs(*planted_cubic(a, -h if MUT != "plant-shape" else Fr(1) / (2 * a)))
    fails = (not hypb["M>0"]) and not all(hypb.values())
    check("planted_shape", classified and all(hyp.values()) and fails,
          f"planted -2a (w - h)(w^2 - 1), a = 1/100, h = {h}: huge roots {huge}, other roots {small}, classified as 1 huge + 2 at +-1: {classified}; "
          f"hypotheses {hyp}; with the huge root at +1/(2a) (M < 0) the hypotheses fail: {fails} ({[k for k, v in hypb.items() if not v]})")


def c_shape():
    U = Fr(1, 500)
    S, B, A = param_boxes(U)
    om_range = Fr(300) if MUT == "shape-om" else Fr(100)
    om = Iv(-om_range, om_range)
    kap, mu, c, M = cubic_params(S, B, om)
    # (i)
    ok_i = kap.lo >= Fr(99, 100) and mu.lo >= Fr(99, 100) and c.lo >= 1 and c.hi <= Fr(101, 100) and M.lo >= Fr(79, 100)
    # (ii) P'(-1/(4a)) = (4M - 3 kappa)/(8a) + 4 a c > 0
    ok_ii = (4 * M - 3 * kap).lo > 0
    # (iii) 32 a^2 P(-1/(4a)) = kappa - 2M + 32 a^2 (c^2 - c - omega^2) < 0   (mutation: test point -1/(2a): 8 a^2 P(-1/(2a)) = ...)
    if MUT == "shape-3":
        # a^2 P(-1/(2a)) * 8 = kappa/... : P(-1/(2a)) = kappa/(4a^2) - M/(4a^2) - 2c + c^2 - omega^2, times 4a^2
        val = kap - M + 4 * A ** 2 * (c ** 2 - 2 * c - om ** 2)
    else:
        val = kap - 2 * M + 32 * A ** 2 * (c ** 2 - c - om ** 2)
    ok_iii = val.hi < 0
    # (iv)
    P2 = -16 * A * kap - 4 * M + 8 * A * c + c ** 2 - om ** 2
    Pm2 = 16 * A * kap - 4 * M - 8 * A * c + c ** 2 - om ** 2
    Pp2 = -24 * A * kap - 4 * M + 4 * A * c          # P'(2): with P'(0) = 4ac > 0 it puts the maximum w_M of P in (0, 2)
    ok_iv = P2.hi < 0 and Pm2.hi < 0 and Pp2.hi < 0
    # (v) a real root in [-2, 2] forces omega^2 <= c^2 + 16 a kappa + 8 a c
    bound = c ** 2 + 16 * A * kap + 8 * A * c
    ok_v = bound.hi <= Fr(105, 100) ** 2
    # (vi) |omega| <= 1/2
    om5 = Iv(Fr(-1, 2), Fr(1, 2))
    kap5, mu5, c5, M5 = cubic_params(S, B, om5)
    P0 = c5 ** 2 - om5 ** 2
    Ph = -A * kap5 * Fr(1, 4) - M5 * Fr(1, 4) + c5 ** 2 - om5 ** 2
    Pmh = -M5 * Fr(1, 4) - 2 * A * c5 + c5 ** 2 - om5 ** 2
    ok_vi = P0.lo > 0 and Ph.lo > 0 and Pmh.lo > 0
    # symbolic part: P' = -6 a kappa w^2 - 2 M w + 4 a c, discriminant, P'(0), P'(-1/(4a))
    a_, k_, M_, c_ = sp.symbols("a kappa M c", positive=True)
    Pp = sp.diff(Pcubic(w_, a_, k_, M_, c_, om_), w_)
    disc = sp.expand(sp.discriminant(Pp, w_))
    ok_s = (sp.expand(disc - (4 * M_ ** 2 + 96 * a_ ** 2 * k_ * c_)) == 0 and sp.simplify(Pp.subs(w_, 0) - 4 * a_ * c_) == 0
            and sp.simplify(Pp.subs(w_, -1 / (4 * a_)) - ((4 * M_ - 3 * k_) / (8 * a_) + 4 * a_ * c_)) == 0)
    check("shape", ok_i and ok_ii and ok_iii and ok_iv and ok_v and ok_vi and ok_s,
          f"sigma, b <= 1/500, a <= 1/1000, |omega| <= {om_range}: (i) kappa >= 0.99, mu >= 0.99, c in [1, 1.01], M >= 0.79: {ok_i} "
          f"(kappa {kap}, M {M}); (ii) 4M - 3 kappa > 0: {ok_ii} {(4 * M - 3 * kap)}; (iii) {'P(-1/(2a))' if MUT == 'shape-3' else '32 a^2 P(-1/(4a))'} < 0: "
          f"{ok_iii} {val}; (iv) P(+-2) < 0 and P'(2) < 0: {ok_iv}; (v) omega^2 <= {bound} <= 1.05^2: {ok_v}; (vi) |omega| <= 1/2: P(0), P(1/2), P(-1/2) > 0: {ok_vi}; "
          f"symbolic P' (discriminant 4M^2 + 96 a^2 kappa c, P'(0) = 4ac, P'(-1/(4a))): {ok_s}")


# ---------------------------------------------------------------------------------------------------------------
# M5 ehat
# ---------------------------------------------------------------------------------------------------------------
def c_ehat():
    a_, k_, M_, c_, om, x = sp.symbols("a kappa M c omega x", positive=True)
    # exact identity: kappa x^2 (1 - 2 a x) - ((kappa - M) x^2 - 4 a c x - (omega^2 - c^2)) = -P(-x)
    idt = sp.expand(k_ * x ** 2 * (1 - 2 * a_ * x) - ((k_ - M_) * x ** 2 - 4 * a_ * c_ * x - (om ** 2 - c_ ** 2)) + Pcubic(-x, a_, k_, M_, c_, om))
    s_, b_ = sp.symbols("sigma b")
    kmm = sp.expand(cubic_params(s_, b_, 0)[0] - cubic_params(s_, b_, 0)[1] + (2 * s_ + 2 * s_ ** 2 - 6 * b_ * s_))
    ok_id = idt == 0 and kmm == 0
    U = Fr(1, 500)
    alpha = Fr(2) if MUT == "ehat-alpha" else Fr(51, 25)
    beta = Fr(181, 50)
    f = Iv(0, 1)
    u = Iv(0, U)
    sg, b = u * f, u * (1 - f)       # sigma = u f, b = u (1 - f), u in [0, U] (u -> the value of sigma + b)
    kap = 1 - b + sg
    c = 1 + 2 * sg
    # T2 / u = [2 sigma + 2 sigma^2 + 6 b sigma + 4 u^2 (c + c^2)] / (kappa u) over the simplex (u in (0, U], f in [0, 1])
    uu = Iv(0, U)
    t2 = (2 * f + 2 * uu * f ** 2 + 6 * uu * f * (1 - f) + 4 * uu * (c + c ** 2))
    t2_ok = t2.hi / kap.lo <= alpha
    t1 = 2 * (Iv(-1, 1))   # |b - sigma| <= u: 2|b - sigma|/kappa <= 2u/kappa
    t1_ok = Fr(2) / kap.lo <= Fr(201, 100)
    t0 = 16 * uu   # 64 a^2 <= 16 u^2 (a^2 = sigma b <= u^2/4), so 64 a^2/(kappa u) <= 16 U/kappa
    t0_ok = t0.hi / kap.lo <= Fr(321, 10000)
    y_ok = Fr(321, 10000) <= (beta - Fr(201, 100)) * Fr(1, 50)
    check("ehat", ok_id and t2_ok and t1_ok and t0_ok and y_ok,
          f"exact identity kappa x^2 (1 - 2ax) - ((kappa - M) x^2 - 4acx - (omega^2 - c^2)) = -P(-x) and kappa - mu = -(2 sigma + 2 sigma^2 - 6 b sigma): {ok_id}; "
          f"with 1/x <= 4a and y >= 1/50: T2 = [2 sigma + 2 sigma^2 + 6 b sigma + 4 u^2 (c + c^2)]/(kappa u) <= {float(t2.hi / kap.lo):.5f} <= alpha = {alpha}: {t2_ok}; "
          f"T1 = 2 |b - sigma|/(kappa u) <= {float(Fr(2) / kap.lo):.5f} <= 201/100: {t1_ok}; T0 = 64 a^2/(kappa u) <= {float(t0.hi / kap.lo):.5f} <= 321/10000: {t0_ok}; "
          f"321/10000 <= (beta - 201/100)/50 = {float((beta - Fr(201, 100)) / 50):.5f}: {y_ok}; so |e-hat| <= u (alpha y^2 + beta y), beta = {beta}")


# ---------------------------------------------------------------------------------------------------------------
# P2 planted_gronwall, M6 theorem_a: the bootstrap chain
# ---------------------------------------------------------------------------------------------------------------
def bootstrap_chain(z0, cf, ubar, a_z=Fr(1, 2), kappa=None):
    """Stage 1: |zeta| <= B e^(-kappa t), kappa = 2 - a_z/2, B = |z0| + cf ubar/(2 - kappa).
    Stage 2: |zeta| <= e^H (|z0| + cf ubar t) e^(-2t), H = B/(2 kappa) >= int |zeta|/2 dt. Returns (kappa, B, H, e^H bound)."""
    kappa = 2 - a_z / 2 if kappa is None else Fr(kappa)
    B = abs(Fr(z0)) + Fr(cf) * Fr(ubar) / (2 - (2 - a_z / 2))
    H = B / (2 * kappa)
    return kappa, B, H, exp_upper(H)


def c_planted_gronwall():
    cf = Fr(1, 2) if MUT == "plant-gron" else Fr(1)
    kap_over = Fr(2) if MUT == "A-kappa" else None
    ok_all, rows = True, []
    for z0, ubar in [(Fr(1, 4), Fr(1, 1000)), (Fr(-1, 4), Fr(1, 100)), (Fr(0), Fr(1, 100)), (Fr(1, 10), Fr(1, 100))]:
        kappa, B, H, eH = bootstrap_chain(z0, cf, ubar, kappa=kap_over)
        s1 = s2 = True
        ratio_end = None
        for k in range(0, 201):
            t = k / 4
            exact = abs((float(z0) + float(ubar) * t) * math.exp(-2 * t))     # exact solution of zeta' = -2 zeta + ubar e^(-2t)
            b1 = float(B) * math.exp(-float(kappa) * t)
            b2 = float(eH) * (abs(float(z0)) + float(cf) * float(ubar) * t) * math.exp(-2 * t)
            s1 &= b1 >= exact * (1 - 1e-12)
            s2 &= b2 >= exact * (1 - 1e-12)
            if k == 200 and exact > 0:
                ratio_end = b2 / exact
        sharp = ratio_end is None or z0 < 0 or abs(ratio_end / float(eH) - 1) <= 1e-6      # for zeta_0 >= 0 the rate part is attained: ratio -> e^H
        rows.append((float(z0), float(ubar), s1, s2, ratio_end))
        ok_all &= s1 and s2 and sharp
    check("planted_gronwall", ok_all,
          "zeta' = -2 zeta + ubar e^(-2t), exact (zeta_0 + ubar t) e^(-2t), c_f = 1: stage-1 bound B e^(-7t/4) and the rate bound e^H (|zeta_0| + c_f ubar t) e^(-2t) "
          "are >= the exact solution on t in [0, 50], and for zeta_0 >= 0 the ratio of the rate bound to the exact solution at t = 50 equals e^H to 1e-6: " + "; ".join(f"(z0 {z}, ubar {u}) {a}/{b} ratio at t = 50 {r if r is None else round(r, 4)}" for z, u, a, b, r in rows)
          + f"; c_f used {cf}, kappa override {kap_over}")


def c_theorem_a():
    alpha, beta = Fr(51, 25), Fr(181, 50)
    U = Fr(1, 500)
    u0 = Fr(1, 100) if MUT == "A-u0" else Fr(1, 1000)
    cf = Fr(40) if MUT == "A-cf" else Fr(50)
    ce = alpha * Fr(9, 2) ** 2 + beta * Fr(9, 2)
    ok_ce = ce <= 58
    ok_cf = cf >= Fr(58, 2) + Fr(9, 2) ** 2
    ubar = exp_upper(Fr(9, 2) * u0) * u0              # sigma_0 <= u_0, so e^(9 sigma_0/2) <= e^(9 u_0/2)
    ok_ubar = ubar <= Fr(10046, 10 ** 7) and ubar <= U and u0 <= U
    kappa, B, H, eH = bootstrap_chain(Fr(1, 4), cf, ubar)
    ok_B = B <= Fr(451, 1000) and B < Fr(1, 2)
    ok_H = H <= Fr(129, 1000) and eH <= Fr(114, 100)
    ok_range = Fr(7, 2) >= Fr(1, 50)
    check("theorem_a", ok_ce and ok_cf and ok_ubar and ok_B and ok_H and ok_range,
          f"c_e = alpha (9/2)^2 + beta (9/2) = {float(ce):.4f} <= 58: {ok_ce}; c_f = {cf} >= 58/2 + (9/2)^2 = 49.25: {ok_cf}; ubar = e^(9 u0/2) u0 = {float(ubar):.6g} "
          f"<= 1.0046e-3 and <= 1/500: {ok_ubar}; kappa = {kappa}, B = 1/4 + 4 c_f ubar = {float(B):.5f} <= 0.451 < 1/2: {ok_B}; "
          f"H = B/(2 kappa) = {float(H):.5f} <= 0.129, e^H <= {float(eH):.5f} <= 1.14: {ok_H}; y in [7/2, 9/2] within [1/50, inf): {ok_range}")


# ---------------------------------------------------------------------------------------------------------------
# M7 numeric_a
# ---------------------------------------------------------------------------------------------------------------
def huge_state(Gv, sg, y):
    """huge root x of P(-x) = 0 by Newton, and e-hat by the exact fixed-point form (no cancellation): returns (x, e-hat)."""
    om = 2.0 / y
    b, a = Gv * Gv * sg, Gv * sg
    kap = 1 - b + sg
    mu = 1 + 3 * sg - b + 2 * sg * sg - 6 * b * sg
    c = 1 + 2 * sg
    M = mu + (b - sg) * om
    x = M / (2 * a * kap)
    for _ in range(100):
        f = 2 * a * kap * x ** 3 - M * x ** 2 - 4 * a * c * x + c * c - om * om
        d = 6 * a * kap * x ** 2 - 2 * M * x - 4 * a * c
        dx = f / d
        x -= dx
        if abs(dx) <= 1e-15 * abs(x):
            break
    km = -(2 * sg + 2 * sg * sg - 6 * b * sg) - (b - sg) * om
    eh = y * y * (km - 4 * a * c / x - (om * om - c * c) / (x * x)) / kap
    return x, eh


def rk4(f, t, Y, h):
    k1 = f(t, Y)
    k2 = f(t + h / 2, [y + h / 2 * k for y, k in zip(Y, k1)])
    k3 = f(t + h / 2, [y + h / 2 * k for y, k in zip(Y, k2)])
    k4 = f(t + h, [y + h * k for y, k in zip(Y, k3)])
    return [y + h / 6 * (a + 2 * b + 2 * c + d) for y, a, b, c, d in zip(Y, k1, k2, k3, k4)]


_NFUN = {}


def direct_funs():
    if not _NFUN:
        Ng = sp.lambdify((G, p, q, nu), N_TH, "math")
        Nq = sp.lambdify((G, p, q, nu), sp.diff(N_TH, q), "math")
        _NFUN["N"], _NFUN["Nq"] = Ng, Nq
    return _NFUN["N"], _NFUN["Nq"]


def c_numeric_a():
    alpha, beta = 51 / 25, 181 / 50
    Bc, cf = 0.451, 50.0
    rate = 3.0 if MUT == "num-rate" else 2.0
    h, T, slack = 1 / 400, 20.0, 1e-12
    worst = {"stage1": 1e300, "rate": 1e300, "ehat": 1e300, "G": 1e300}
    ok = True
    ncases = 0
    for G0 in [1e-2, 0.3, 1.0, 3.0, 30.0]:
        for z0 in [-0.25, 0.0, 0.25]:
            u0 = 1 / 1000
            s0 = u0 / (1 + G0 * G0)
            ubar = math.exp(4.5 * s0) * u0

            def rhs(t, Y, s0=s0):
                lnG, y = Y
                sg = s0 * math.exp(-2 * t)
                x, eh = huge_state(math.exp(lnG), sg, y)
                return [sg * y, y * (4 - y) / 2 - sg * y * y + eh / 2]

            Y, t = [math.log(G0), 4 + z0], 0.0
            traj = []
            for k in range(int(T / h)):
                Y = rk4(rhs, t, Y, h)
                t += h
                if k % 10 == 9:
                    traj.append((t, Y[0], Y[1]))
            lnGinf = traj[-1][1]
            for t, lnG, y in traj:
                sg = s0 * math.exp(-2 * t)
                Gv = math.exp(lnG)
                _, eh = huge_state(Gv, sg, y)
                z = abs(y - 4)
                b1 = Bc * math.exp(-7 * t / 4)
                b2 = 1.14 * (abs(z0) + cf * ubar * t) * math.exp(-rate * t)
                b3 = sg * (1 + Gv * Gv) * (alpha * y * y + beta * y)
                b4 = 2.25 * sg
                a4 = abs(lnG - lnGinf)
                ok &= z <= b1 + slack and z <= b2 + slack and abs(eh) <= b3 + slack and a4 <= b4 + slack
                for key, bnd, act in (("stage1", b1, z), ("rate", b2, z), ("ehat", b3, abs(eh)), ("G", b4, a4)):
                    if act > 1e-9:
                        worst[key] = min(worst[key], bnd / act)
            ncases += 1
    # direct integration of N = 0 against the reduced system, at h = 1/800, and the step study 1/800 vs 1/1600 for both integrators
    Nf, Nq = direct_funs()

    def direct_compare(G0, z0, hh):
        s0 = 1 / (1000 * (1 + G0 * G0))

        def rhs_red(t, Y):
            sg = s0 * math.exp(-2 * t)
            x, eh = huge_state(math.exp(Y[0]), sg, Y[1])
            return [sg * Y[1], Y[1] * (4 - Y[1]) / 2 - sg * Y[1] ** 2 + eh / 2]

        def rhs_dir(t, Y):
            Gv, pv = Y
            sg = s0 * math.exp(-2 * t)
            qv = -pv * pv / (2 * Gv * sg)
            for _ in range(60):
                nn = Nf(Gv, pv, qv, math.sqrt(sg))
                dq = nn / Nq(Gv, pv, qv, math.sqrt(sg))
                qv -= dq
                if abs(dq) <= 1e-14 * abs(qv):
                    break
            return [pv, qv]

        Yr, Yd, t = [math.log(G0), 4 + z0], [G0, s0 * G0 * (4 + z0)], 0.0
        for k in range(int(round(6 / hh))):
            Yr = rk4(rhs_red, t, Yr, hh)
            Yd = rk4(rhs_dir, t, Yd, hh)
            t += hh
        sg = s0 * math.exp(-2 * t)
        return math.exp(Yr[0]), Yr[1], Yd[0], Yd[1] / (sg * Yd[0])

    maxdev, maxstep = 0.0, 0.0
    for G0, z0 in [(1.0, 0.2), (0.3, -0.2)]:
        a800 = direct_compare(G0, z0, 1 / 800)
        a1600 = direct_compare(G0, z0, 1 / 1600)
        maxdev = max(maxdev, abs(a800[0] - a800[2]), abs(a800[1] - a800[3]))
        maxstep = max(maxstep, *(abs(u - v) for u, v in zip(a800, a1600)))
    ok_dir = maxdev <= 1e-9 and maxstep <= 1e-10
    check("numeric_a", ok and ok_dir,
          f"{ncases} cases (G_0 in {{1e-2, 0.3, 1, 3, 30}}, zeta_0 in {{-1/4, 0, 1/4}}, u_0 = 1/1000), t to 20 (h = 1/400): every step satisfies |zeta| <= 0.451 e^(-7t/4), "
          f"|zeta| <= 1.14 (|zeta_0| + 50 ubar t) e^(-{rate:g}t), |e-hat| <= u (alpha y^2 + beta y), |ln (G/G_inf)| <= (9/4) sigma (slack 1e-12): {ok}; "
          f"smallest bound/actual ratios (actual > 1e-9): stage1 {worst['stage1']:.3g}, rate {worst['rate']:.3g}, e-hat {worst['ehat']:.3g}, G {worst['G']:.3g}; "
          f"direct N = 0 vs reduced system at t = 6 (h = 1/800), max difference in G and y {maxdev:.2e} <= 1e-9, and h = 1/800 vs 1/1600 for both integrators differ by at most {maxstep:.2e} <= 1e-10: {ok_dir}")


# ---------------------------------------------------------------------------------------------------------------
# P5 planted_series, M8 series, M9 stage_symbol, M10 expansion_link
# ---------------------------------------------------------------------------------------------------------------
def c_planted_series():
    # (a) q + p = 0: F_inf + c nu, no logarithm
    g1, log1 = stage_solve(lambda g: g.theta().theta() + g.theta(), Ser(), range(1, 7), 0, lambda k: k + 1, free=(1, Fr(5, 3)))
    ok1 = all(l[4] == "ok" for l in log1) and g1.d == {(1, 0): Fr(5, 3)} and log1[0][3] == [[1] + [0] * (log1[0][2] - 1)]
    # (b) q + 4 p = 8 nu^4 (theta^2 g + 4 theta g - 8 nu^4): without logs an obstruction at k = 4, with logs a_{4,1} = -2
    forcing = Ser({(4, 0): Fr(0 if MUT == "plant-series" else 8)})
    N2 = lambda g: g.theta().theta() + 4 * g.theta() - forcing
    g2, log2 = stage_solve(N2, Ser(), range(1, 7), 0, lambda k: 0, free=None)
    ob = any(l[4] == "obstruction" and l[0] == 4 for l in log2)
    g3, log3 = stage_solve(N2, Ser(), range(1, 7), 0, lambda k: 2, free=(4, Fr(0)))
    ok3 = all(l[4] == "ok" for l in log3) and g3.d.get((4, 1)) == -2 and set(g3.d) <= {(4, 0), (4, 1)}
    check("planted_series", ok1 and ob and ok3,
          f"q + p = 0: g = {dict(g1.d)} (c nu only, null direction a_(1,0)): {ok1}; q + 4p = 8 nu^4: no-log ansatz stops with an obstruction at k = 4: {ob}; "
          f"with logs a_(4,1) = {g3.d.get((4, 1))} (= -2) and nothing else: {ok3}")


def table_value(key, Fv, Cv):
    F, C = Fv, Cv
    if key == (4, 1):
        return F * (F ** 2 + (8 if MUT == "series-8" else 7))
    if key == (6, 0):
        return -(107 * F ** 6 - 134 * F ** 4 + 587 * F ** 2 + 48 * C * F ** 3 - 384 * C * F + 144 * C ** 2) / (216 * F)
    if key == (6, 1):
        return -2 * (F ** 2 + 7) * (F ** 3 - 8 * F + 6 * C) / 9
    if key == (6, 2):
        return -2 * F * (F ** 2 + 7) ** 2 / 3
    if key == (8, 3):
        return F * (F ** 2 + 7) ** 3 / 2
    if key == (2, 0):
        return -2 * F
    if key == (4, 0):
        return C
    raise KeyError(key)


def c_series():
    ok = True
    notes = []
    for Fv in [Fr(1), Fr(-3, 2), Fr(2, 7), Fr(5)]:
        for Cv in [Fr(0), Fr(3), Fr(-5, 2)]:
            g, log = solve_series(Fv, Cv, a2=(-Fr(11, 5) * Fv if MUT == "series-a2" else None))
            stages_ok = len(log) == 7 and all(l[4] == "ok" for l in log) and all(not l[5] for l in log)
            uniq = stages_ok and all((l[3] == [] if l[0] != 4 else (len(l[3]) == 1 and l[3][0][0] == 1 and all(x == 0 for x in l[3][0][1:]))) for l in log)
            tab = stages_ok and all(g.d.get(k, Fr(0)) == table_value(k, Fv, Cv) for k in [(2, 0), (4, 0), (4, 1), (6, 0), (6, 1), (6, 2), (8, 3)])
            odd = stages_ok and all(v == 0 for (n, j), v in g.d.items() if n % 2 == 1)
            if stages_ok:
                Nr = N_series(Ser.const(Fv) + g)
                low = min((n for (n, j), v in Nr.d.items() if v != 0), default=99)
            else:
                low = None
            res_ok = low == 16
            ok &= stages_ok and uniq and tab and odd and res_ok
            notes.append((Fv, Cv, stages_ok, uniq, tab, odd, low))
    bad = [n for n in notes if not (n[2] and n[3] and n[4] and n[5] and n[6] == 16)]
    check("series", ok,
          f"exact Fraction solver, F_inf in {{1, -3/2, 2/7, 5}} x C in {{0, 3, -5/2}} (12 cases), stages k = 3..9: all solvable with no lower-order residual, full column rank "
          f"except k = 4 (null space spanned by a_(4,0)), table of conical_limit.md (a_(2,0), a_(4,0), a_(4,1), a_(6,0..2), a_(8,3)), odd coefficients zero, first nonzero "
          f"coefficient of N at nu^16: {ok}" + ("" if ok else f"; failing cases (F, C, stages, uniq, table, odd, first order): {bad[:3]}"))


def c_stage_symbol():
    Fs = sp.Symbol("F", real=True, nonzero=True)
    base = {G: Fs - 2 * Fs * nu ** 2, p: 4 * Fs * nu ** 2, q: -8 * Fs * nu ** 2}
    # base: G = F - 2F nu^2, p = theta G = 4 F nu^2, q = theta p = -8 F nu^2 (theta nu = -nu)
    Nq = sp.expand(sp.diff(N_TH, q).subs(base, simultaneous=True))
    Np = sp.expand(sp.diff(N_TH, p).subs(base, simultaneous=True))
    NG = sp.expand(sp.diff(N_TH, G).subs(base, simultaneous=True))
    low = lambda e: min(m[0] for m in sp.Poly(e, nu).monoms())
    A = sp.Poly(Nq, nu).coeff_monomial(nu ** 6)
    Bp = sp.Poly(Np, nu).coeff_monomial(nu ** 6)
    k = sp.Symbol("k")
    I = sp.factor(A * k ** 2 - Bp * k)         # delta G = eps nu^k: delta p = -k eps nu^k, delta q = k^2 eps nu^k
    want = -128 * Fs ** 3 * k * (k - (5 if MUT == "stage-I" else 4))
    ok = low(Nq) == 6 and low(Np) == 6 and low(NG) >= 8 and sp.simplify(I - want) == 0
    check("stage_symbol", ok, f"dN/dq = {A} nu^6 + ..., dN/dp = {Bp} nu^6 + ..., dN/dG of order nu^{low(NG)}: I(k) = {I}, roots {sp.solve(I, k)}; "
          f"so stage k >= 3, k != 4 has a nonzero indicial coefficient for every k: {ok}")


def series_inverse(Fs):
    """1/Fs for a series with constant term F (nu^0 T^0) and the rest of order nu^2 or higher."""
    F0 = Fs.d[(0, 0)]
    h = Ser({k: v / F0 for k, v in Fs.d.items() if k != (0, 0)})
    tot, pw = Ser.const(Fr(1)), Ser.const(Fr(1))
    for _ in range(NM // 2 + 1):
        pw = pw * (-h)
        tot = tot + pw
    return tot * (Fr(1) / F0)


def c_expansion_link():
    ok = True
    rows = []
    for Fv in [Fr(1), Fr(2)]:
        vals = {}
        for Cv in [Fr(0), Fr(3), Fr(-5, 2)]:
            g, _ = solve_series(Fv, Cv)
            Ft = Ser.const(Fv) + g
            pt = g.theta()
            rt = Ser({(n - 2, j): v for (n, j), v in pt.d.items()})           # r~ = p~ / sigma
            zeta = rt * series_inverse(Ft) - 4
            vals[Cv] = zeta
            lnk = g * (Fr(1) / Fv)
            c20 = zeta.d.get((2, 0), Fr(0))
            c21 = zeta.d.get((2, 1), Fr(0))
            want21 = -4 * (Fv ** 2 + (6 if MUT == "link-slope" else 7))
            ok_c = zeta.d.get((0, 0), Fr(0)) == 0 and c21 == want21 and c20 == Fv ** 2 + 15 - 4 * Cv / Fv and lnk.d.get((2, 0)) == -2
            ok &= ok_c
            rows.append((Fv, Cv, c20, c21, ok_c))
        slope = (vals[Fr(3)].d.get((2, 0), 0) - vals[Fr(0)].d.get((2, 0), 0)) / 3
        ok &= slope == -4 / Fv
    check("expansion_link", ok,
          "zeta~ = eta^3 F~'/F~ - 4 from the series (exact): nu^0 coefficient 0, sigma t coefficient -4 (F^2 + 7) = f_sigma, sigma coefficient F^2 + 15 - 4C/F "
          "(affine in C, slope -4/F), ln(F~/F) = -2 sigma + ...: " + "; ".join(f"F={f} C={c}: ({a}, {b})" for f, c, a, b, _ in rows[:3]) + f" ... {ok}")


# ---------------------------------------------------------------------------------------------------------------
# M11 numeric_b
# ---------------------------------------------------------------------------------------------------------------
def c_numeric_b():
    mp.mp.dps = 34
    G0, zeta0, eta0 = mp.mpf(1), mp.mpf("0.2"), mp.mpf(45)
    sig0 = 1 / eta0 ** 2

    def hugeroot(Gv, sg, y):
        om = 2 / y
        b, a = Gv * Gv * sg, Gv * sg
        kap = 1 - b + sg
        mu = 1 + 3 * sg - b + 2 * sg ** 2 - 6 * b * sg
        c = 1 + 2 * sg
        M = mu + (b - sg) * om
        x = M / (2 * a * kap)
        for _ in range(80):
            f = 2 * a * kap * x ** 3 - M * x ** 2 - 4 * a * c * x + c * c - om * om
            d = 6 * a * kap * x ** 2 - 2 * M * x - 4 * a * c
            dx = f / d
            x -= dx
            if abs(dx) < mp.mpf(10) ** -31 * x:
                break
        km = -(2 * sg + 2 * sg ** 2 - 6 * b * sg) - (b - sg) * om
        return (y * y * (km - 4 * a * c / x - (om * om - c * c) / (x * x)) / kap)

    def rhs(t, Y):
        lnG, y = Y
        sg = sig0 * mp.e ** (-2 * t)
        eh = hugeroot(mp.e ** lnG, sg, y)
        return [sg * y, y * (4 - y) / 2 - sg * y * y + eh / 2]

    etaend = mp.mpf(10) ** 7
    f = mp.odefun(rhs, 0, [mp.log(G0), 4 + zeta0], tol=mp.mpf(10) ** -28, degree=12)
    Yend = f(mp.log(etaend / eta0))
    Finf = mp.e ** Yend[0] * mp.e ** (2 / etaend ** 2)           # int_t^inf sigma y dt = 2 sigma + O(sigma^2 t)
    to_fr = lambda x, dig=30: Fr(int(mp.nint(x * mp.mpf(10) ** dig)), 10 ** dig)
    Fi = to_fr(Finf)

    def series_at(g, eta, m):
        T, nv = mp.log(eta), 1 / eta
        tot = Finf
        for (n, j), v in g.d.items():
            if n <= 2 * m:
                coef = mp.mpf(v.numerator) / mp.mpf(v.denominator)
                if MUT == "numB-coef" and (n, j) == (6, 0):
                    coef *= mp.mpf("1.01")
                tot += coef * nv ** n * T ** j
        return tot

    Fref = mp.e ** f(mp.log(mp.mpf(1000) / eta0))[0]

    def fC(Cv):
        g, _ = solve_series(Fi, to_fr(Cv, 20))
        return series_at(g, mp.mpf(1000), 4) - Fref

    c0, c1 = mp.mpf(0), mp.mpf(1)
    f0, f1 = fC(c0), fC(c1)
    for _ in range(6):
        if f1 == f0 or abs(f1) < mp.mpf(10) ** -32:       # fit converged (the secant would divide by zero)
            break
        c2 = c1 - f1 * (c1 - c0) / (f1 - f0)
        c0, f0, c1 = c1, f1, c2
        f1 = fC(c1)
    Cfit = c1
    g, _ = solve_series(Fi, to_fr(Cfit, 20))
    table, ok = [], True
    for eta in [100, 200, 400, 800]:
        Fn = mp.e ** f(mp.log(mp.mpf(eta) / eta0))[0]
        row = []
        for m in [1, 2, 3, 4]:
            err = abs(Fn - series_at(g, mp.mpf(eta), m))
            em = err / ((mp.mpf(1) / eta) ** (2 * m + 2) * (1 + mp.log(eta)) ** m)
            row.append(float(em))
            ok &= em <= 1e5
        table.append((eta, row))
    check("numeric_b", ok,
          f"mpmath (30 digits) from G_0 = 1, zeta_0 = 1/5, eta_0 = 45 to 1e7: F_inf = {mp.nstr(Finf, 20)}, C = {mp.nstr(Cfit, 12)} (residual of the fit at eta = 1000 "
          f"{mp.nstr(abs(f1), 3)}); e_m = |F - F^(2m)|/(nu^(2m+2) (1 + ln eta)^m) <= 1e5 for m = 1..4: " + "; ".join(f"eta {e}: " + ", ".join(f"{x:.3g}" for x in r) for e, r in table))


# ---------------------------------------------------------------------------------------------------------------
# P4 planted_local, M12 local_lemma
# ---------------------------------------------------------------------------------------------------------------
def local_lemma_applies(Npoly):
    """the hypotheses of the local lemma on a polynomial Npoly(G, p, q, nu): (a) the q^2, q^1, q^0 coefficients are divisible by p, p^4, p^4,
    (c) Npoly(G, 0, q, nu) = c3 q^3 with c3 = -2 G nu^2 kappa. Returns (a, c)."""
    Pq = sp.Poly(sp.expand(Npoly), q)
    cs = {m[0]: c for m, c in zip(Pq.monoms(), Pq.coeffs())}
    c3, c2, c1, c0 = (cs.get(k, sp.Integer(0)) for k in (3, 2, 1, 0))
    div = lambda e, k: sp.fraction(sp.cancel(e / p ** k))[1] == 1
    a_ok = div(c2, 1) and div(c1, 4) and div(c0, 4)
    kap = 1 - (G ** 2 - 1) * nu ** 2
    c_ok = sp.expand(Npoly.subs(p, 0) - (-2 * G * nu ** 2 * kap) * q ** 3) == 0
    return a_ok, c_ok


def c_planted_local():
    Nosc = -2 * G * nu ** 2 * q ** 3 * (1 - (G ** 2 - 1) * nu ** 2) if MUT == "plant-local" else nu ** 2 * (q - p) + G
    a_ok, c_ok = local_lemma_applies(Nosc)
    applies = a_ok and c_ok
    # F = cos(eta): eta = pi, G = -1, p = eta F' = 0, q = eta^2 F'' + eta F' = pi^2: N_osc = 0 although q != 0
    ev = sp.simplify(((nu ** 2 * (q - p) + G)).subs({nu: 1 / sp.pi, q: sp.pi ** 2, p: 0, G: -1}))
    counter = ev == 0
    check("planted_local", (not applies) and counter,
          f"N_osc = nu^2 (q - p) + G (F'' + F = 0): divisibility {a_ok}, N(G, 0, q, nu) = c3 q^3 {c_ok}: applies {applies} (it must not); F = cos eta has p = 0 with "
          f"q = pi^2 != 0 at eta = pi and N_osc = {ev}: {counter}")


def c_local_lemma():
    Nuse = N_TH + (q ** 2 if MUT == "local-c2" else 0)
    a_ok, c_ok = local_lemma_applies(Nuse)
    # ratios used by the root bound
    sg = nu ** 2
    kap, mu, c = kap_mu_c(sg, G ** 2 * sg)
    r2 = sp.simplify(((-p * (2 * G * (G ** 2 - 1) * sg ** 2 + p * mu)) / (-2 * G * sg * kap)) - p * ((G ** 2 - 1) * sg / kap + p * mu / (2 * G * sg * kap)))
    r1 = sp.simplify((4 * G * sg * c * p ** 4) / (-2 * G * sg * kap) + 2 * c * p ** 4 / kap)
    ratios_ok = r2 == 0 and r1 == 0
    # numeric test of |q| <= 3 Theta on random roots (complex roots included)
    import numpy as np
    rng = random.Random(7)
    cf = sp.lambdify((G, p, nu), [-2 * G * nu ** 2 * (1 - (G ** 2 - 1) * nu ** 2),
                                 -p * (2 * G * (G ** 2 - 1) * nu ** 4 + p * (1 + 3 * nu ** 2 - G ** 2 * nu ** 2 + 2 * nu ** 4 - 6 * G ** 2 * nu ** 4)),
                                 4 * G * nu ** 2 * (1 + 2 * nu ** 2) * p ** 4,
                                 p ** 4 * ((1 + 2 * nu ** 2) ** 2 * p ** 2 - 4 * G ** 2 * nu ** 4)], "math")
    viol, nroots = 0, 0
    for _ in range(10000):
        Gv = rng.choice([-1, 1]) * 10 ** rng.uniform(-1.5, 1.0)
        nv = 10 ** rng.uniform(-3, -0.5)
        pv = rng.choice([-1, 1]) * 10 ** rng.uniform(-12, 0)
        c3, c2, c1, c0 = cf(Gv, pv, nv)
        if c3 == 0:
            continue
        th = max(abs(c2 / c3), abs(c1 / c3) ** 0.5, abs(c0 / c3) ** (1 / 3))
        for z in np.roots([c3, c2, c1, c0]):
            nroots += 1
            viol += abs(z) > 3 * th * (1 + 1e-9)
    ok = a_ok and c_ok and ratios_ok and viol == 0
    check("local_lemma", ok,
          f"N: c2 divisible by p, c1 and c0 by p^4: {a_ok}; N(G, 0, q, nu) = -2 G nu^2 kappa q^3, c3 != 0 iff G != 0 and kappa != 0: {c_ok}; "
          f"c2/c3 = p [(G^2-1) sigma/kappa + p mu/(2 G sigma kappa)], c1/c3 = -2 c p^4/kappa: {ratios_ok}; |q| <= 3 max(|c2/c3|, |c1/c3|^(1/2), |c0/c3|^(1/3)) on "
          f"{nroots} random roots (1e4 parameter points): {viol} violations")


# ---------------------------------------------------------------------------------------------------------------
# M13 regime_b, M14 small_branch, M15 huge_branch
# ---------------------------------------------------------------------------------------------------------------
def c_regime_b():
    U = Fr(1, 500)
    ymax = Fr(1, 5) if MUT == "regimeb-y" else Fr(1, 50)
    kap_lo, mu_hi, c_hi = Fr(998, 1000), Fr(100601, 100000), Fr(1004, 1000)
    a_hi = U / 2
    th1 = (U + ymax * mu_hi / 2) / kap_lo
    ok2 = 2 * c_hi / kap_lo <= Fr(9, 4)
    th2 = a_hi * ymax * Fr(3, 2)
    tau3 = Fr(35, 10000)
    ok3 = tau3 ** 3 >= 2 * a_hi ** 2 * ymax / kap_lo
    three = 3 * max(th1, th2, tau3)
    ok = ok2 and ok3 and three <= Fr(37, 1000) and Fr(37, 1000) <= Fr(1, 10)
    check("regime_b", ok,
          f"|y| <= {ymax}, a <= 1/1000, u <= 1/500: Theta_1/|p| <= (U + |y| mu/2)/kappa <= {float(th1):.5f}; sqrt(2c/kappa) <= 3/2: {ok2}, Theta_2/|p| <= {float(th2):.2e}; "
          f"Theta_3/|p| <= (2 a^2 |y|/kappa)^(1/3) <= {float(tau3)} (cube {float(tau3 ** 3):.3e} >= {float(2 * a_hi ** 2 * ymax / kap_lo):.3e}: {ok3}); 3 max <= {float(three):.4f} <= 0.037 <= 1/10: {ok}; "
          f"so |s| <= |r|/10 and |r|' >= 1.9 |r|")


def c_small_branch():
    t = sp.Symbol("t")
    Gf, pf, wf = sp.Function("G")(t), sp.Function("p")(t), sp.Function("w")(t)
    om = 2 * Gf * sp.exp(-2 * t) / pf
    dln = sp.diff(sp.log(om), t).subs({sp.Derivative(Gf, t): pf, sp.Derivative(pf, t): pf ** 2 * wf})
    ok_id = sp.simplify(dln - (pf / Gf - 2 - pf * wf)) == 0
    U = Fr(1, 500)
    worst = 4 * U + 8 * (U / 2)       # |p| (1/G + 2) <= (2 G sigma/|omega|)(1/G + 2) <= 4 sigma + 8 a for |omega| >= 1/2
    claim = Fr(199, 100) if MUT == "small-sigma" else Fr(1)
    ok_b = -2 + worst <= -claim
    w0, p0, tt = sp.symbols("w0 p0 tt", positive=True)
    sol = 1 / (1 / p0 - w0 * tt)
    ok_ode = sp.simplify(sp.diff(sol, tt) - w0 * sol ** 2) == 0
    integ = sp.integrate(1 / (1 / p0 + w0 * tt), (tt, 0, sp.oo)) if False else sp.limit(sp.log(1 + w0 * p0 * tt) / w0, tt, sp.oo)
    ok_div = integ == sp.oo
    check("small_branch", ok_id and ok_b and ok_ode and ok_div,
          f"(ln|omega|)' = p/G - 2 - p w: {ok_id}; for |omega| >= 1/2, |w| <= 2, sigma + b <= 1/500: (ln|omega|)' <= -2 + 4 sigma + 8 a <= -2 + {float(worst):.4f} <= -{float(claim)}: {ok_b}; "
          f"p' = w p^2 with constant w: p = 1/(1/p0 - w t) (blow-up at t = 1/(p0 w) for w > 0): {ok_ode}; for w < 0, int p dt = ln(1 + |w| p0 t)/|w| -> infinity: {ok_div}")


def c_huge_branch():
    U = Fr(1, 500)
    alpha, beta = Fr(51, 25), Fr(181, 50)
    # y' = 2y + a y^2 w - sigma y^2 and e-hat = y^2 (1 + 2 a w) give y' = y (4 - y)/2 - sigma y^2 + e-hat/2
    a_, w, y, sg = sp.symbols("a w y sigma")
    ok_id = sp.simplify((2 * y + a_ * y ** 2 * w - sg * y ** 2) - (y * (4 - y) / 2 - sg * y ** 2 + y ** 2 * (1 + 2 * a_ * w) / 2)) == 0
    L1 = Fr(3, 2) - U - U * (alpha + beta) / 2
    L2 = Fr(1, 8) - U * Fr(15, 4) - U * (alpha * Fr(15, 4) + beta) / 2
    L3 = Fr(1, 8) - Fr(17, 4) * U * (1 + alpha / 2) - U * beta / 2
    slope_pos = Fr(1, 2) - U * (1 + alpha / 2) > 0
    c1 = Fr(3, 2) if MUT == "huge-thr" else Fr(149, 100)
    ok = ok_id and L1 >= c1 and L2 >= Fr(1, 10) and L3 >= Fr(104, 1000) and slope_pos
    check("huge_branch", ok,
          f"y' = 2y + a y^2 w - sigma y^2 = y(4-y)/2 - sigma y^2 + e-hat/2 with e-hat = y^2 (1 + 2aw): {ok_id}; with |e-hat| <= u (alpha y^2 + beta y), u, sigma <= 1/500: "
          f"y' >= {float(L1):.4f} y >= {float(c1)} y on [1/50, 1]: {L1 >= c1}; y' >= {float(L2):.4f} y >= 1/10 on [1, 15/4]: {L2 >= Fr(1, 10)}; "
          f"y' <= -{float(L3):.4f} y <= -0.104 y on [17/4, inf) (bracket increasing: {slope_pos}): {L3 >= Fr(104, 1000)}; w <= -1/(4a) gives r' <= 2r - r^2/(4G) (s = sigma r^2 w)")


# ---------------------------------------------------------------------------------------------------------------
# M16 classify_numeric
# ---------------------------------------------------------------------------------------------------------------
def cubic_coeffs(Gv, sg, y):
    om = 2.0 / y
    b, a = Gv * Gv * sg, Gv * sg
    kap = 1 - b + sg
    mu = 1 + 3 * sg - b + 2 * sg * sg - 6 * b * sg
    c = 1 + 2 * sg
    M = mu + (b - sg) * om
    return a, kap, M, c, om


def track_root(w, Gv, sg, y):
    a, kap, M, c, om = cubic_coeffs(Gv, sg, y)
    for _ in range(100):
        f = -2 * a * kap * w ** 3 - M * w ** 2 + 4 * a * c * w + c * c - om * om
        d = -6 * a * kap * w ** 2 - 2 * M * w + 4 * a * c
        dw = f / d
        w -= dw
        if abs(dw) <= 1e-14 * max(abs(w), 1e-300):
            break
    return w


def start_root(kind, Gv, sg, y):
    import numpy as np
    a, kap, M, c, om = cubic_coeffs(Gv, sg, y)
    rts = sorted(z.real for z in np.roots([-2 * a * kap, -M, 4 * a * c, c * c - om * om]) if abs(z.imag) < 1e-9 * max(1, abs(z)))
    if kind == "huge":
        return rts[0]
    return rts[1] if kind == "minus" else rts[-1]


def classify_run(kind, Gv0, y0, T, h=1 / 200, eta0=100.0):
    s0 = 1 / eta0 ** 2
    st = {"w": start_root(kind, Gv0, s0, y0)}

    def rhs(t, Y):
        sg = s0 * math.exp(-2 * t)
        Gv = math.exp(Y[0])
        w = track_root(st["w"], Gv, sg, Y[1])
        st["w"] = w
        a = Gv * sg
        return [sg * Y[1], 2 * Y[1] + a * Y[1] ** 2 * w - sg * Y[1] ** 2]

    Y, t = [math.log(Gv0), y0], 0.0
    ws = []
    for k in range(int(T / h)):
        Y = rk4(rhs, t, Y, h)
        t += h
        ws.append((t, st["w"]))
        if kind == "huge" and abs(Y[1]) > 1e6:        # on the small roots y grows like e^(2t) while p = sigma G y stays bounded: no stop
            return "blowup", t, Y, ws
    return "global", t, Y, ws


def c_classify_numeric():
    ok_all, notes = True, []
    # huge root, p > 0: converges to y = 4
    for G0 in [0.5, 2.0]:
        for y0 in [0.5, 8.0, 30.0]:
            kind, t, Y, _ = classify_run("huge", G0, y0, 15.0)
            good = kind == "global" and abs(Y[1] - 4) <= 1e-6
            ok_all &= good
            notes.append(f"H+ G0={G0} y0={y0}: {kind}, |y-4| = {abs(Y[1] - 4):.1e}")
    # huge root, p < 0: finite-time blow-up (mutation: predicts convergence)
    for G0 in [0.5, 2.0]:
        for y0 in [-0.5, -4.0]:
            kind, t, Y, _ = classify_run("huge", G0, y0, 15.0)
            good = (kind == "global" and abs(Y[1] - 4) <= 1e-6) if MUT == "class-sign" else (kind == "blowup" and t < 15.0)
            ok_all &= good
            notes.append(f"H- G0={G0} y0={y0}: {kind} at t = {t:.2f}")
    # small roots stay in their bands on [3, 15]
    for y0 in [3.0, 6.0, -3.0, -6.0]:
        for kind, lo, hi in (("minus", -2.0, -0.5), ("plus", 0.5, 2.0)):
            res, t, Y, ws = classify_run(kind, 1.0, y0, 15.0)
            band = all(lo < w < hi for tt, w in ws if tt >= 3.0)
            good = res == "global" and band
            ok_all &= good
            notes.append(f"S y0={y0} {kind}: w in ({lo}, {hi}) on [3, 15]: {band}")
    check("classify_numeric", ok_all,
          "eta_0 = 100, RK4 h = 1/200, root tracked by continuity: huge root with p > 0 converges to y = 4 (|y - 4| <= 1e-6 at t = 15), huge root with p < 0 blows up "
          "before t = 15, small roots keep w in (-2, -1/2) or (1/2, 2) on [3, 15]: " + "; ".join(notes))


# ---------------------------------------------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------------------------------------------
SHARED = {}


def run_reduction():
    SHARED["Ss"] = c_reduction()


def run_hyp_line_real():
    if "Ss" not in SHARED:
        SHARED["Ss"] = reduction_Ss_only()
    c_hyp_line_real(SHARED["Ss"])


def reduction_Ss_only():
    R = real_R()
    pt = {s_: -8 * G, r_: 4 * G, sig: 0}
    return sp.simplify(-sp.diff(R, sig).subs(pt) / sp.diff(R, s_).subs(pt))


CHECKS = [
    ("planted_hyp", c_planted_hyp), ("planted_gronwall", c_planted_gronwall), ("planted_shape", c_planted_shape),
    ("planted_local", c_planted_local), ("planted_series", c_planted_series),
    ("n_identity", c_n_identity), ("reduction", run_reduction), ("hyp_line_real", run_hyp_line_real), ("shape", c_shape), ("ehat", c_ehat),
    ("theorem_a", c_theorem_a), ("numeric_a", c_numeric_a), ("series", c_series), ("stage_symbol", c_stage_symbol),
    ("expansion_link", c_expansion_link), ("numeric_b", c_numeric_b), ("local_lemma", c_local_lemma), ("regime_b", c_regime_b),
    ("small_branch", c_small_branch), ("huge_branch", c_huge_branch), ("classify_numeric", c_classify_numeric),
]


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
    if MUT is None and "no-controls" not in OPTS and ONLY is None and "skip-num" not in OPTS:
        print("\n== Controls: each mutation must make its target checks fail, and no others ==", flush=True)
        for name, (want, what) in MUTATIONS.items():
            args = [sys.executable, __file__, f"--mutate={name}"]
            if not any(w in NUMERIC for w in want):
                args.append("--skip-num")
            pr = subprocess.run(args, capture_output=True, text=True)
            if "logdir" in OPTS:
                with open(os.path.join(OPTS["logdir"], f"ctl_{name}.log"), "w") as fh:
                    fh.write(pr.stdout + pr.stderr)
                with open(os.path.join(OPTS["logdir"], f"ctl_{name}.rc"), "w") as fh:
                    fh.write(f"{pr.returncode}\n")
            fails = [ln.split(":")[0][7:] for ln in pr.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = all(w in fails for w in want) and pr.returncode == 1
            extra = [f for f in fails if f not in want and f not in failed]
            controls_ok &= fired and not extra
            print(f"[{'FIRED' if fired and not extra else 'SILENT' if not fired else 'EXTRA'}] --mutate={name} ({what}): exit {pr.returncode}, failed {fails}", flush=True)
            if pr.returncode not in (0, 1):
                print(pr.stderr[-2000:])
        print(f"controls: {'all fired on target only' if controls_ok else 'SOME DID NOT FIRE ON TARGET ONLY'}")
    sys.exit(0 if not failed and controls_ok else 1)

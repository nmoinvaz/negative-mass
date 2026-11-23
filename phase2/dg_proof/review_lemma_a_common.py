#!/usr/bin/env python3
"""Shared code of the #138 review of dg_proof Lemma A, L1, T1 and T2.

Written from RESULTS.md alone (sections 1, 3 and the static Zerilli V of 4.2). Exact arithmetic only (sympy, Fraction).
Run as a script it runs every planted problem PL1..PL7 once:

    python3 review_lemma_a_common.py [--mutate=NAME]

Exit 0: all pass. Exit 1: a real check failed. Exit 2: a planted check failed (the caller stops before any real data).
`--mutate=NAME` switches on one deliberate error (list in MUTATIONS); every script that imports this module accepts it.
"""
import math
import sys
from fractions import Fraction

import sympy as sp

MUTATIONS = [
    "plant_recur", "plant_exp", "plant_rel", "plant_sc", "plant_min", "plant_zer",
    "cn_n4to3", "cn_poch", "cn_sign5", "V_coeff", "P2_claim", "P1_6c_5c", "P0_12to13",
    "rec_12to11", "rec_m4to3", "poch_ratio", "ratio_claim",
    "t1_coef18", "t1_drop_sum2", "t1_bound3", "t2_shift", "t2_flip_d0", "t2_window",
]

MUT = None
for _a in sys.argv[1:]:
    if _a.startswith("--mutate="):
        MUT = _a.split("=", 1)[1]
if MUT is not None and MUT not in MUTATIONS:
    print(f"unknown mutation {MUT!r}")
    sys.exit(3)

# ---------------------------------------------------------------------------------------------------------------------
# harness

FAILS = []          # ids of failed checks, in order
_PLANTED_FAILS = []


def check(cid, ok, msg="", planted=False):
    if ok:
        print(f"PASS {cid}")
    else:
        print(f"FAIL {cid}: {msg}")
        FAILS.append(cid)
        if planted:
            _PLANTED_FAILS.append(cid)
    sys.stdout.flush()
    return ok


def check_range(cid, ls, fn, show=3):
    """fn(l) -> (ok, msg). One PASS line if every l passes, else a FAIL line per failing l (first `show`)."""
    bad = []
    for l in ls:
        ok, msg = fn(l)
        if not ok:
            bad.append((l, msg))
    if not bad:
        print(f"PASS {cid} [l={ls[0]}..{ls[-1]}]")
    else:
        for l, msg in bad[:show]:
            print(f"FAIL {cid}[l={l}]: {msg}")
            FAILS.append(f"{cid}[l={l}]")
        if len(bad) > show:
            print(f"FAIL {cid}: {len(bad) - show} more failing l, last l={bad[-1][0]}")
    sys.stdout.flush()
    return not bad


def gate_planted():
    """Stop with exit 2 if a planted check failed."""
    if _PLANTED_FAILS:
        print(f"STOP: planted check(s) failed: {', '.join(_PLANTED_FAILS)}; no real data was examined")
        sys.exit(2)


def finish():
    n = len(FAILS)
    print(f"{'FAILED' if n else 'ALL PASSED'}: {n} failing check(s)" + (f", first {FAILS[0]}" if n else ""))
    sys.exit(1 if n else 0)


# ---------------------------------------------------------------------------------------------------------------------
# Lemma A closed form (typed from RESULTS section 1):
#   c_n = (n-2)((n-4) l(l+1) + n-1) (l-n+3)_{2n-4} / (2^n ((n-1)!)^2),  chi = sum_{n=1}^{l+2} c_n r^n
# with the Pochhammer convention (x)_{-k} = 1/((x-1)...(x-k)) (needed for c_1).

def _n4():
    return 3 if MUT == "cn_n4to3" else 4


def _pk():
    return 2 if MUT == "cn_poch" else 3


def poch_sym(a, k):
    """(a)_k for integer k (negative k by the convention above); a symbolic or numeric."""
    if k >= 0:
        return sp.Mul(*[a + j for j in range(k)])
    return 1 / sp.Mul(*[a - j for j in range(1, -k + 1)])


def poch_q(a, k):
    """Exact (a)_k as Fraction for integer a, integer k (0 if the product contains 0)."""
    if k >= 0:
        return Fraction(math.prod(range(a, a + k)))
    r = 1
    for j in range(1, -k + 1):
        r *= a - j
    return Fraction(1, r)


def cn_sym(n, l):
    """c_n, n integer, l symbolic."""
    L = l * (l + 1)
    v = (n - 2) * ((n - _n4()) * L + n - 1) * poch_sym(l - n + _pk(), 2 * n - 4) / (2**n * sp.factorial(n - 1) ** 2)
    if MUT == "cn_sign5" and n == 5:
        v = -v
    return sp.together(v)


def cn_q(n, l):
    """c_n as Fraction, n >= 1 and l integers."""
    L = l * (l + 1)
    v = Fraction((n - 2) * ((n - _n4()) * L + n - 1)) * poch_q(l - n + _pk(), 2 * n - 4) / (2**n * _fact(n - 1) ** 2)
    if MUT == "cn_sign5" and n == 5:
        v = -v
    return v


_FACT = {0: 1}


def _fact(n):
    if n not in _FACT:
        _FACT[n] = n * _fact(n - 1)
    return _FACT[n]


def c_of_l(l):
    return (l + 2) * (l - 1)


def rs_q(l):
    return Fraction(6, c_of_l(l))


R = sp.Symbol("r")


def chi_poly(l):
    """chi as sympy Poly in R over QQ, from the closed form, n = 1..l+2 (higher n give 0 and are skipped)."""
    top = l + 2
    d = {}
    for n in range(1, top + 3):
        v = cn_q(n, l)
        if v != 0:
            d[(n,)] = sp.Rational(v.numerator, v.denominator)
    return sp.Poly.from_dict(d, R, domain="QQ")


# ---------------------------------------------------------------------------------------------------------------------
# static Zerilli V of RESULTS 4.2 (M = -1, f = 1 + 2/r), lambda = c/2:
#   V = 2f [lambda^2 r^2 ((lambda+1) r - 3) + 9 (lambda r - 1)] / (r^3 (lambda r - 3)^2)

def V_zerilli(r, c):
    lam = c / sp.Integer(2)
    f = 1 + 2 / r
    nine = 9 * (lam * r - (2 if MUT == "V_coeff" else 1))
    return 2 * f * (lam**2 * r**2 * ((lam + 1) * r - 3) + nine) / (r**3 * (lam * r - 3) ** 2)


def zerilli_residual(psi, V, r):
    """f (f psi_r)_r - V psi with f = 1 + 2/r (d/dx = f d/dr)."""
    f = 1 + (1 if MUT == "plant_zer" else 2) / r
    return sp.cancel(sp.together(f * sp.diff(f * sp.diff(psi, r), r) - V * psi))


# hand-derived coefficients of the chi ODE  P2 chi'' + P1 chi' + P0 chi = 0  (derived from V before any code, see README)

def P_hand(r, c):
    P2 = -(r**2) * (r + 2) * (c * r - 6)
    P1 = 2 * c * r**3 + (5 if MUT == "P1_6c_5c" else 6) * c * r**2 - 12 * r
    P0 = c**2 * r**2 - 6 * c * r + (13 if MUT == "P0_12to13" else 12)
    return P2, P1, P0


def P2_results(r, c):
    """P2 as RESULTS section 1 states it."""
    return -(r**2) * (r + 2) * (c * r - (5 if MUT == "P2_claim" else 6))


def rec_hand(m, c):
    """My recurrence on r^m: B0 c_m + B1 c_{m-1} + B2 c_{m-2} = 0 (typed from the hand derivation)."""
    b0 = (11 if MUT == "rec_12to11" else 12) * (m - 1) ** 2
    b1 = (m - 2) * (6 * (m - 1) - 2 * c * (m - (3 if MUT == "rec_m4to3" else 4)))
    b2 = c * (c - (m - 2) * (m - 5))
    return b0, b1, b2


# ---------------------------------------------------------------------------------------------------------------------
# generic routines (the ones the planted problems test)

def three_term(P2, P1, P0, x):
    """Coefficients B_s(n) such that (P2 y'' + P1 y' + P0 y) applied to x^n is sum_s B_s(n) x^{n+s}.
    A series y = sum a_k x^k solves the ODE iff for every power k: sum_s B_s(k-s) a_{k-s} = 0.
    Returns {s: B_s(n)} with n the sympy Symbol 'n'."""
    n = sp.Symbol("n")
    p2 = sp.Poly(sp.expand(P2), x).all_coeffs()[::-1]
    p1 = sp.Poly(sp.expand(P1), x).all_coeffs()[::-1]
    p0 = sp.Poly(sp.expand(P0), x).all_coeffs()[::-1]
    out = {}
    for j, a in enumerate(p2):
        out[j - 2] = out.get(j - 2, 0) + a * n * (n - 1)
    for j, a in enumerate(p1):
        out[j - 1] = out.get(j - 1, 0) + a * n
    if MUT != "plant_recur":
        for j, a in enumerate(p0):
            out[j] = out.get(j, 0) + a
    out = {s: sp.expand(v) for s, v in out.items()}
    return {s: v for s, v in out.items() if v != 0}      # zero entries (low-degree terms of P2) are not shifts


def series_residual(B, a, k):
    """sum_s B_s(k-s) a(k-s) for a callable a (a(i) = 0 for i < 0 is the caller's duty)."""
    n = sp.Symbol("n")
    tot = 0
    for s, bs in B.items():
        i = k - s
        if i < 0:
            continue
        tot += bs.subs(n, i) * a(i)
    return sp.simplify(tot)


def exponents(P2, P1, P0, x, x0):
    """Indicial roots {root: multiplicity} at x0 of P2 y'' + P1 y' + P0 y = 0 (regular singular point)."""
    t = sp.Symbol("t")
    Q2, Q1, Q0 = (sp.expand(P.subs(x, x0 + t)) for P in (P2, P1, P0))
    p0 = sp.cancel(t * Q1 / Q2).subs(t, 0)
    q0 = sp.cancel(t**2 * Q0 / Q2).subs(t, 0)
    rho = sp.Symbol("rho")
    ind = rho * (rho - 1) + (-p0 if MUT == "plant_exp" else p0) * rho + q0
    return sp.roots(sp.Poly(sp.expand(ind), rho))


def local_relations(P2, P1, P0, x, x0):
    """At a point x0 with P2(x0) = 0: from the ODE E = 0 and its first two derivatives (built by differentiating
    the operator applied to a generic y, then reading the coefficients of y, y', ..., y'''').
    Returns dict: yp (y'/y), ypp (y''/y), res_coef (coefficient of y''' in E''), res_rest (E'' without y''' term, after
    substituting yp, ypp)."""
    y = sp.Function("y")
    Y = sp.symbols("Y0:5")
    E = P2 * y(x).diff(x, 2) + P1 * y(x).diff(x) + P0 * y(x)
    E1 = sp.diff(E, x)
    E2 = sp.diff(E1, x)

    def coeffs(expr):
        for k in range(4, 0, -1):
            expr = expr.subs(y(x).diff(x, k), Y[k])
        expr = expr.subs(y(x), Y[0])
        return {k: sp.simplify(sp.diff(sp.expand(expr), Y[k]).subs(x, x0)) for k in range(5)}

    c0, c1, c2 = coeffs(E), coeffs(E1), coeffs(E2)
    assert sp.simplify(P2.subs(x, x0)) == 0
    yp = sp.simplify(-c0[0] / c0[1])
    a_, b_, c_ = c1[2], c1[1], c1[0]
    if MUT == "plant_rel":
        b_ = b_ - c0[0]          # P1' instead of P1' + P0
    ypp = sp.simplify(-(b_ * yp + c_) / a_)
    res_rest = sp.simplify(c2[2] * ypp + c2[1] * yp + c2[0])
    return {"yp": yp, "ypp": ypp, "res_coef": sp.simplify(c2[3]), "res_rest": res_rest}


def sign_changes(co):
    """Number of sign changes of a coefficient list, zeros skipped."""
    if MUT == "plant_sc":
        return sum(1 for u, v in zip(co, co[1:]) if u * v <= 0)
    s = [1 if u > 0 else -1 for u in co if u != 0]
    return sum(1 for u, v in zip(s, s[1:]) if u != v)


def poly_min(p, x, a, b):
    """Exact minimum of the polynomial p(x) on [a, b]: (value, argument)."""
    dp = sp.diff(p, x)
    cands = [r_ for r_ in sp.real_roots(sp.Poly(dp, x)) if a <= r_ <= b] if dp != 0 else []
    pts = cands if MUT == "plant_min" else cands + [sp.sympify(a), sp.sympify(b)]
    vals = [(sp.simplify(p.subs(x, q)), sp.simplify(q)) for q in pts]
    return min(vals, key=lambda t: t[0])


# ---------------------------------------------------------------------------------------------------------------------
# planted problems

def ode_from_solutions(y1, y2, x):
    """P2 y'' + P1 y' + P0 y = 0 with solutions y1, y2, by the Wronskian determinant, common factor removed."""
    W = y1 * sp.diff(y2, x) - sp.diff(y1, x) * y2
    P1 = -(y1 * sp.diff(y2, x, 2) - sp.diff(y1, x, 2) * y2)
    P0 = sp.diff(y1, x) * sp.diff(y2, x, 2) - sp.diff(y1, x, 2) * sp.diff(y2, x)
    g = sp.gcd(sp.gcd(sp.Poly(sp.expand(W), x), sp.Poly(sp.expand(P1), x)), sp.Poly(sp.expand(P0), x))
    P2, P1, P0 = (sp.cancel(sp.expand(P) / g.as_expr()) for P in (W, P1, P0))
    for y in (y1, y2):
        assert sp.simplify(P2 * sp.diff(y, x, 2) + P1 * sp.diff(y, x) + P0 * y) == 0
    return P2, P1, P0


X = sp.Symbol("x")


def planted_odes():
    a = sp.Rational(3, 2)
    t = X - a
    E2 = ode_from_solutions(1 + 2 * t + 5 * t**2, t**3 * (1 + X), X)
    E3 = ode_from_solutions(1 + t**2, t**7 * (X**2 + 3), X)
    E4 = (X, sp.Rational(1, 2), -sp.Integer(1))
    E5 = (X, sp.Integer(1), sp.Integer(0))        # x y'' + y' = 0, double root
    return a, E2, E3, E4, E5


def run_planted(ids):
    """Run the planted checks named in ids (subset of PL1..PL7)."""
    N_, M_, k_ = sp.symbols("N M k")
    if "PL1" in ids or "PL2" in ids:
        P2 = (1 + X) ** 2
        P1 = -(7 + 4 - 1) * (1 + X)
        P0 = sp.Integer(28)
        B = three_term(P2, P1, P0, X)
        for NN in (7, 4):
            ok = all(series_residual(B, lambda i, NN=NN: sp.binomial(NN, i), k) == 0 for k in range(0, 13))
            if "PL1" in ids:
                check(f"PL1 numeric C({NN},k) k=0..12", ok, "residual not 0", planted=True)
        Bs = three_term((1 + X) ** 2, -(N_ + M_ - 1) * (1 + X), N_ * M_, X)
        n = sp.Symbol("n")
        ratio = {-2: (N_ - k_) * (N_ - k_ - 1) / ((k_ + 1) * (k_ + 2)), -1: (N_ - k_) / (k_ + 1), 0: sp.Integer(1)}
        # sum_s B_s(k-s) a_{k-s} / a_k, with a_{k-s}/a_k = ratio[s] for s in {-2,-1,0}
        sym = sum(Bs[s].subs(n, k_ - s) * ratio[s] for s in Bs)
        if "PL1" in ids:
            check("PL1 symbolic (N,M,k)", sp.simplify(sym) == 0, f"residual {sp.simplify(sym)}", planted=True)
        if "PL2" in ids:
            bad = lambda i: sp.binomial(7, i) + (1 if i == 3 else 0)
            res = [series_residual(B, bad, k) for k in range(0, 13)]
            check("PL2 wrong sequence has nonzero residual", any(v != 0 for v in res), "residual stayed 0", planted=True)
    if "PL3" in ids or "PL4" in ids:
        a, E2, E3, E4, E5 = planted_odes()
    if "PL3" in ids:
        for name, E, x0, want in (("E2", E2, a, {0: 1, 3: 1}), ("E3", E3, a, {0: 1, 7: 1}),
                                  ("E4", E4, sp.Integer(0), {0: 1, sp.Rational(1, 2): 1}),
                                  ("E5", E5, sp.Integer(0), {0: 2})):
            got = exponents(*E, X, x0)
            check(f"PL3 exponents {name}", dict(got) == want, f"got {dict(got)}, want {want}", planted=True)
    if "PL4" in ids:
        a, E2, E3, E4, E5 = planted_odes()
        r2 = local_relations(*E2, X, a)
        check("PL4 E2 y'/y = 2, y''/y = 10", r2["yp"] == 2 and r2["ypp"] == 10, f"{r2}", planted=True)
        check("PL4 E2 resonance (coef 0, rest 0)", r2["res_coef"] == 0 and r2["res_rest"] == 0, f"{r2}", planted=True)
        r3 = local_relations(*E3, X, a)
        check("PL4 E3 y'/y = 0, y''/y = 2", r3["yp"] == 0 and r3["ypp"] == 2, f"{r3}", planted=True)
        r4 = local_relations(*E4, X, sp.Integer(0))
        check("PL4 E4 y'/y = 2, y''/y = 4/3", r4["yp"] == 2 and r4["ypp"] == sp.Rational(4, 3), f"{r4}", planted=True)
    if "PL5" in ids:
        t = sp.Symbol("t")
        for name, p, want_sc, want_roots in (
                ("(t-2)(t+1)(t+3)", (t - 2) * (t + 1) * (t + 3), 1, 1),
                ("(t-1)(t-2)(t+3)", (t - 1) * (t - 2) * (t + 3), 2, 2),
                ("t^2-t+1", t**2 - t + 1, 2, 0),
                ("t^3+2t^2", t**3 + 2 * t**2, 0, 0)):
            co = sp.Poly(sp.expand(p), t).all_coeffs()
            nroots = sp.Poly(sp.expand(p), t).count_roots(sp.Rational(1, 10**9), sp.oo)
            sc = sign_changes(co)
            check(f"PL5 {name}", sc == want_sc and nroots == want_roots, f"changes {sc} (want {want_sc}), roots {nroots} (want {want_roots})",
                  planted=True)
    if "PL6" in ids:
        x = sp.Symbol("x")
        p = x**3 - 3 * x + 2
        for a_, b_, wv, wa in ((0, 3, 0, 1), (2, 3, 4, 2), (-3, 0, -16, -3)):
            try:
                v, arg = poly_min(p, x, a_, b_)
                ok = (v == wv and arg == wa)
                msg = f"got ({v}, {arg}), want ({wv}, {wa})"
            except ValueError as e:
                ok, msg = False, f"error {e}"
            check(f"PL6 min on [{a_},{b_}]", ok, msg, planted=True)
    if "PL7" in ids:
        r = sp.Symbol("r", positive=True)
        V = 2 * (r + 2) / r**3
        check("PL7 psi=r^2 residual 0", zerilli_residual(r**2, V, r) == 0, "residual nonzero", planted=True)
        check("PL7 perturbed V has nonzero residual", zerilli_residual(r**2, sp.Rational(11, 10) * V, r) != 0, "residual stayed 0",
              planted=True)


if __name__ == "__main__":
    print(f"review_lemma_a_common: planted problems PL1..PL7 (mutation: {MUT})")
    run_planted({"PL1", "PL2", "PL3", "PL4", "PL5", "PL6", "PL7"})
    if _PLANTED_FAILS:
        print(f"planted check(s) failed: {', '.join(_PLANTED_FAILS)}")
        sys.exit(2)
    finish()

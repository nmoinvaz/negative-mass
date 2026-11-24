#!/usr/bin/env python3
"""Do conical shifts u = F(x/rho) with det R_3 = 0 have limits as x/rho -> +-infinity? (issue #159)

N(F, F', F''; eta) = 8 rho^6 D, eta = x/rho, is the ODE of det_classify.md section 6. With theta = eta d/d eta,
G = F, p = theta G, q = theta^2 G, nu = 1/eta and T = t = ln eta, N is a polynomial in (G, p, q, nu^2).
  C1  n_build         N rebuilt from the jets of F(x/rho) (D retyped) equals det_classify's typed N, the printed
                      parts of degree 3 and 4, and the jets of a test F at exact points
  C2  theta_parity    the theta-form, nu only as nu^0, nu^2, nu^4; N(F0, -F1, F2; -eta) = N, so eta -> -inf mirrors +inf
  C3  planted         the formal solver of C6, C7 on four ODEs with known transseries
  C4  polygon         Newton polygon of N: balances for F_inf != 0 and F_inf = 0
  C5  a2_zero         the branch a_2 = 0 and every other start k = 1, ..., 8
  C6  series_pure     the series without logs fails at eta^-10 with -512 F_inf^4 (F_inf^2 + 7)
  C7  series_log      the series with logs through eta^-9, its residual, the one free parameter
  C8  reduced_system  the first-order system in (G, r = eta^2 p, nu^2): a line of equilibria, eigenvalues (0, -2, -2)
  C9  linear_modes    exponents of the linearization about the series: 0 and 4 only
  C10 num_inward      mpmath integration of N = 0 inward from the series, compared to 1e-10 at a second eta
  C11 num_outward     the same outward from perturbed data; the limit exists and follows the series
  C12 farfield_scale  eps^3 and eps^4 coefficients of D for u_inf + eps w, w = -2 u_inf rho^2/x^2
Imports: det_classify.conical_N (the typed N, as a comparison target only). D is retyped here.
Usage: conical_limit.py [--no-controls] [--mutate=NAME] [--only=check1,check2] [--skip-num] [--logdir=DIR].
Exit 1 if a check fails or a control does not fire."""
import os
import subprocess
import sys
import time

import mpmath as mp
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import det_classify as dc  # noqa: E402  (only dc.conical_N is used)

# name: (checks that must fail, what the mutation does)
MUTATIONS = {
    "d-sign": (["n_build"], "builds Y as u u_xrho/2 - u_x u_rho"),
    "parity-break": (["theta_parity"], "replaces one term eta^3 of the typed N by eta^2"),
    "theta-sub": (["theta_parity"], "takes F'' = (q + p)/eta^2 in the substitution"),
    "plant-wrong": (["planted"], "P1 uses F'' + F'/eta"),
    "plant-res": (["planted"], "P4 with the forcing 0 in place of 8 nu^4"),
    "poly-drop": (["polygon"], "drops the term -p^2 q^2 from N before the polygon"),
    "poly0-edge": (["polygon"], "adds nu^4 q^2 to N (an edge at k = 2 appears for F_inf = 0)"),
    "a2zero-edge": (["a2_zero"], "replaces -2 G nu^2 q^3 by -2 G nu^3 q^3 (an edge at k = 3)"),
    "series-8": (["series_pure"], "expects F_inf^4 (F_inf^2 + 8)"),
    "series-nolog": (["series_log"], "uses jmax = 0 in the log run"),
    "reduced-eig": (["reduced_system"], "expects the eigenvalue +2 in place of -2 for r"),
    "modes-a2": (["linear_modes"], "the series with a_2 = -(11/5) F_inf in place of -2 F_inf"),
    "num-nolog": (["num_inward", "num_outward"], "the series without the term a_{4,1}, as data and as comparison"),
    "ff-drop": (["farfield_scale"], "drops the eps^4 coefficient"),
}
NUMERIC = ("num_inward", "num_outward")

RESULTS = []
OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate")
ONLY = OPTS["only"].split(",") if "only" in OPTS else None

eta = sp.Symbol("eta", positive=True)
F0, F1, F2 = sp.symbols("F0 F1 F2")
G, p, q, nu = sp.symbols("G p q nu", positive=True)
Fi = sp.Symbol("Fi", real=True, nonzero=True)  # F_inf != 0
T = sp.Symbol("T", real=True)  # ln eta


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def wanted(name):
    if "skip-num" in OPTS and name in NUMERIC and not (MUT and name in MUTATIONS[MUT][0]):
        return False
    return ONLY is None or name in ONLY


def zero(e):
    return sp.expand(e) == 0


# ---------------------------------------------------------------------------------------------------------------
# N from D
# ---------------------------------------------------------------------------------------------------------------
def D_of(P, Q, Z, Y, h):
    return h * (h**2 - P**2 - Q**2 + Y**2 + Z**2) - P * (Y**2 - Z**2) + 2 * Q * Y * Z


def D_static(u, x, rho, ysign=1):
    """det R_3 of the static shift u(x, rho), the block of det_classify.md Setup with b = u, b_t = 0."""
    ux, ur = sp.diff(u, x), sp.diff(u, rho)
    uxx, uxr, urr = sp.diff(u, x, 2), sp.diff(u, x, rho), sp.diff(u, rho, 2)
    P = -(u * uxx + ux**2)
    Q = -(urr + ur / rho) / 2
    Z = uxr / 2
    Y = u * uxr / 2 + ysign * ux * ur
    h = ur**2 / 2
    return D_of(P, Q, Z, Y, h)


def build_N(ysign=1):
    """8 rho^6 D for u = F(x/rho) as a polynomial in (F0, F1, F2, eta)."""
    x, rho = sp.symbols("x rho", positive=True)
    F = sp.Function("F")
    D = D_static(F(x / rho), x, rho, ysign).doit().subs(x, eta * rho).doit()
    D = D.replace(lambda e: isinstance(e, sp.Subs), lambda e: e.doit())
    d2, d1 = sp.Derivative(F(eta), (eta, 2)), sp.Derivative(F(eta), eta)
    D = D.subs(d2, F2).subs(d1, F1).subs(F(eta), F0)
    return sp.expand(sp.simplify(8 * rho**6 * D))


N_ETA = None


def N_canonical():
    global N_ETA
    if N_ETA is None:
        N_ETA = build_N()
    return N_ETA


def theta_form(Neta, f2sign=-1):
    """N(F0, F1, F2; eta) with F1 = p nu, F2 = (q - p) nu^2, F0 = G, eta = 1/nu: a polynomial in (G, p, q, nu)."""
    return sp.expand(Neta.subs({F0: G, F1: p / eta, F2: (q + f2sign * p) / eta**2}).subs(eta, 1 / nu))


def N_theta():
    return theta_form(N_canonical())


A0 = p**2 * (p**2 - q) * (p**2 + q)
A1 = G**2 * p**2 * q**2 + 4 * G * p**4 * q - 2 * G * q**3 + 4 * p**6 - 3 * p**2 * q**2
A2 = 2 * (G * (q - p) + p**2) * (G**2 * q**2 + 2 * G * p**3 + 2 * G * p**2 * q + 2 * p**4 - q**2)


def c1_n_build():
    Nb = build_N(-1 if MUT == "d-sign" else 1)
    typed = dc.conical_N()
    a_ok = zero(Nb - typed)
    # (b) lowest parts: degree in (F1, F2) with F0 fixed, then in (F0, F1, F2)
    lam = sp.Symbol("lam")
    w = F1 + eta * F2
    deg3 = 2 * F0 * w**2 * ((F0**2 - 1) * F2 - eta * w)
    deg4 = -w**2 * ((eta**2 + 1) * (eta**2 + 2) * F1**2 + 2 * F0 * (eta * w + F2))
    PL3 = sp.Poly(sp.expand(Nb.subs({F1: lam * F1, F2: lam * F2})), lam)
    PL4 = sp.Poly(sp.expand(Nb.subs({F0: lam * F0, F1: lam * F1, F2: lam * F2})), lam)
    b3 = all(zero(PL3.coeff_monomial(lam**i)) for i in range(3)) and zero(PL3.coeff_monomial(lam**3) - deg3)
    b4 = all(zero(PL4.coeff_monomial(lam**i)) for i in range(4)) and zero(PL4.coeff_monomial(lam**4) - deg4)
    # (c) exact points, direct differentiation of u(x, rho) = F(x/rho)
    x, rho = sp.symbols("x rho", positive=True)
    Fe = lambda s: s**3 - 2 * s + sp.Rational(1, 3)
    u = Fe(x / rho)
    Dd = D_static(u, x, rho, -1 if MUT == "d-sign" else 1)
    pts = [(sp.Rational(3, 2), sp.Rational(2, 3)), (sp.Rational(1, 2), 1), (5, sp.Rational(7, 3)),
           (sp.Rational(1, 7), 4), (sp.Rational(11, 3), sp.Rational(5, 2))]
    c_ok = True
    for xv, rv in pts:
        ev = xv / rv
        Fv = [Fe(eta).subs(eta, ev), sp.diff(Fe(eta), eta).subs(eta, ev), sp.diff(Fe(eta), eta, 2).subs(eta, ev)]
        lhs = Dd.subs({x: xv, rho: rv})
        rhs = Nb.subs({F0: Fv[0], F1: Fv[1], F2: Fv[2], eta: ev}) / (8 * rv**6)
        c_ok &= sp.simplify(lhs - rhs) == 0
    check("n_build", a_ok and b3 and b4 and c_ok,
          f"N = 8 rho^6 D from the jets of F(x/rho) equals det_classify.conical_N(): {a_ok}; lowest part in (F', F'') is the "
          f"printed degree-3 part 2 F w^2((F^2 - 1)F'' - eta w), lower parts 0: {b3}; lowest part in (F, F', F'') is the "
          f"printed degree-4 part -w^2((eta^2+1)(eta^2+2)F'^2 + 2F(eta w + F'')), lower parts 0: {b4}; D from direct "
          f"derivatives of F(x/rho) = (x/rho)^3 - 2x/rho + 1/3 equals N/(8 rho^6) at {len(pts)} exact points: {c_ok}")


def c2_theta_parity():
    Nn = N_canonical()
    if MUT == "parity-break":
        Nn = sp.expand(Nn - 4 * F0 * F1**5 * eta**3 + 4 * F0 * F1**5 * eta**2)
    Nt = theta_form(Nn, 1 if MUT == "theta-sub" else -1)
    P = sp.Poly(Nt, nu)
    powers = sorted(m[0] for m in P.monoms())
    typed = sp.expand(A0 + nu**2 * A1 + nu**4 * A2)
    form_ok = powers == [0, 2, 4] and zero(Nt - typed)
    par_ok = zero(Nn.subs({F1: -F1, eta: -eta}, simultaneous=True) - Nn)
    check("theta_parity", form_ok and par_ok,
          f"N(G, p, q, nu) has the powers nu^{powers} and equals A0 + nu^2 A1 + nu^4 A2 with A0 = p^2 (p^2 - q)(p^2 + q), "
          f"A1 = G^2 p^2 q^2 + 4G p^4 q - 2G q^3 + 4p^6 - 3p^2 q^2, A2 = 2(G(q - p) + p^2)(G^2 q^2 + 2G p^3 + 2G p^2 q + "
          f"2p^4 - q^2): {form_ok}; N(F0, -F1, F2; -eta) = N(F0, F1, F2; eta), so F(-eta) solves N = 0 when F(eta) "
          f"does: {par_ok}")


# ---------------------------------------------------------------------------------------------------------------
# truncated transseries in nu = e^-t with T = t = ln eta
# ---------------------------------------------------------------------------------------------------------------
class Ser:
    """sum c[(n, j)] nu^n T^j, truncated at nu^M."""

    def __init__(self, d=None, M=40):
        self.d = {k: v for k, v in (d or {}).items() if v != 0}
        self.M = M

    def __add__(self, o):
        o = o if isinstance(o, Ser) else Ser({(0, 0): sp.sympify(o)}, self.M)
        r = dict(self.d)
        for k, v in o.d.items():
            r[k] = sp.expand(r.get(k, 0) + v)
        return Ser(r, self.M)

    __radd__ = __add__

    def __mul__(self, o):
        if not isinstance(o, Ser):
            return Ser({k: sp.expand(v * o) for k, v in self.d.items()}, self.M)
        r = {}
        for (n1, j1), v1 in self.d.items():
            for (n2, j2), v2 in o.d.items():
                if n1 + n2 <= self.M:
                    r[(n1 + n2, j1 + j2)] = r.get((n1 + n2, j1 + j2), 0) + v1 * v2
        return Ser({k: sp.expand(v) for k, v in r.items()}, self.M)

    def Dt(self):
        """d/dt, nu = e^-t: d(nu^n T^j) = -n nu^n T^j + j nu^n T^(j-1)."""
        r = {}
        for (n, j), v in self.d.items():
            if n:
                r[(n, j)] = r.get((n, j), 0) - n * v
            if j:
                r[(n, j - 1)] = r.get((n, j - 1), 0) + j * v
        return Ser({k: sp.expand(v) for k, v in r.items()}, self.M)

    def low(self):
        return min((n for (n, j) in self.d), default=None)

    def coeffT(self, n):
        return sum(v * T**j for (m, j), v in self.d.items() if m == n)

    def expr(self):
        return sum(v * nu**n * T**j for (n, j), v in self.d.items())


def eval_terms(terms, fi, g, M):
    """N(G = fi + g, p = g', q = g''; nu) as a Ser truncated at nu^M. terms: [((a, b, c, m), coef)] for G^a p^b q^c nu^m."""
    g = Ser(g.d, M)
    Gs = Ser({(0, 0): fi}, M) + g
    pp = g.Dt()
    qq = pp.Dt()
    cache = {}

    def power(name, base, k):
        if (name, k) not in cache:
            cache[(name, k)] = Ser({(0, 0): 1}, M) if k == 0 else power(name, base, k - 1) * base
        return cache[(name, k)]

    tot = Ser({}, M)
    for (a, b, c, m), co in terms:
        if m > M:
            continue
        t = Ser({(m, 0): co}, M) * power("G", Gs, a) * power("p", pp, b) * power("q", qq, c)
        tot = tot + t
    return tot


def poly_terms(Nthe):
    P = sp.Poly(sp.expand(Nthe), G, p, q, nu)
    return list(zip(P.monoms(), P.coeffs()))


def stage(terms, fi, g, k, js, M, free_prefix="a"):
    """Add nu^k sum_{j in js} u_j T^j to g, take the lowest nonzero order of N, and solve it for the u_j.
    Returns dict(status, order, sol, free, cond): status 'solved' (linear in the u_j), 'obstruction' (the lowest order holds
    no u_j: cond is a condition on the earlier coefficients), 'inconsistent' or 'zero' (N vanishes to order M)."""
    us = {j: sp.Symbol(f"u{j}", real=True) for j in js}
    gk = g + Ser({(k, j): us[j] for j in js})
    E = eval_terms(terms, fi, gk, M)
    o = E.low()
    if o is None:
        return dict(status="zero", order=None)
    ce = sp.expand(E.coeffT(o))
    eqs = [sp.together(c) for c in sp.Poly(ce, T).all_coeffs()]
    if not any(ce.has(u) for u in us.values()):
        return dict(status="obstruction", order=o, cond=sp.factor(sp.together(ce)))
    unk = list(us.values())
    sol = list(sp.linsolve(eqs, unk))
    if not sol:
        return dict(status="inconsistent", order=o, eqs=eqs)
    vals = dict(zip(unk, sol[0]))
    free = sorted(set().union(*[sp.sympify(v).free_symbols for v in vals.values()]) & set(unk), key=str)
    ren = {u: sp.Symbol(f"{free_prefix}{k}_{j}", real=True) for j, u in us.items() if u in free}
    vals = {u: sp.factor(sp.sympify(v).subs(ren)) for u, v in vals.items()}
    vals.update({u: ren[u] for u in ren})
    if not all(sp.simplify(e.subs({u: vals[u] for u in unk})) == 0 for e in eqs):
        return dict(status="inconsistent", order=o, eqs=eqs)
    return dict(status="solved", order=o, sol={j: vals[us[j]] for j in js}, free=[str(v) for v in ren.values()])


def formal_series(Nthe, fi, k0, c0, kmax, jmax_of, pad=8):
    """g = c0 nu^k0 + sum_{k0 < k <= kmax} nu^k sum_{j <= jmax_of(k)} u_j T^j solving N(fi + g, g', g''; nu) = 0 stage by
    stage. Returns (g as Ser, log list); stops at the first stage that is not 'solved'."""
    terms = poly_terms(Nthe)
    g = Ser({(k0, 0): c0}) if c0 != 0 else Ser({})
    log = []
    for k in range(k0 + 1, kmax + 1):
        r = stage(terms, fi, g, k, list(range(jmax_of(k) + 1)), k + pad)
        r["k"] = k
        log.append(r)
        if r["status"] == "zero":
            continue
        if r["status"] != "solved":
            return g, log
        g = g + Ser({(k, j): v for j, v in r["sol"].items()})
    return g, log


def residual_orders(Nthe, fi, g, M):
    """Orders nu^m, m <= M, at which N(fi + g) does not vanish."""
    E = eval_terms(poly_terms(Nthe), fi, g, M)
    return sorted({n for (n, j) in E.d})


def jlog(k):
    return max(1, k // 2)


# ---------------------------------------------------------------------------------------------------------------
# C3 planted truth
# ---------------------------------------------------------------------------------------------------------------
def c3_planted():
    c1 = sp.Symbol("c1", real=True)
    N1 = nu**2 * ((q + p) if MUT != "plant-wrong" else q)  # F'' + 2F'/eta  (mutation: F'' + F'/eta)
    g1, log1 = formal_series(N1, 0, 1, c1, 8, lambda k: 3, pad=4)
    p1 = all(r["status"] == "solved" for r in log1) and len(log1) == 7 and g1.d == {(1, 0): c1}
    t = poly_terms(N1)
    s2 = stage(t, 0, Ser({}), 0, [1, 2, 3], 6)
    p2 = s2["status"] == "solved" and all(s2["sol"][j] == 0 for j in (1, 2, 3))
    N3 = nu**2 * q  # F'' + F'/eta
    t3 = poly_terms(N3)
    s3 = stage(t3, 0, Ser({}), 0, [1, 2, 3], 6)
    p3a = s3["status"] == "solved" and s3["sol"][2] == 0 and s3["sol"][3] == 0 and s3["free"] == ["a0_1"]
    g3, log3 = formal_series(N3, 0, 1, c1, 4, lambda k: 3, pad=4)
    p3b = bool(log3) and log3[0]["status"] == "obstruction" and sp.simplify(log3[0]["cond"] / c1).free_symbols == set()
    forcing = 0 if MUT == "plant-res" else 8
    N4 = q + 4 * p - forcing * nu**4
    g4n, log4n = formal_series(N4, 0, 0, 0, 8, lambda k: 0, pad=4)
    p4a = log4n[-1]["status"] == "obstruction" and log4n[-1]["k"] == 4 and log4n[-1]["order"] == 4
    g4, log4 = formal_series(N4, 0, 0, 0, 8, lambda k: 1, pad=4)
    a40, a41 = sp.Symbol("a4_0", real=True), None
    p4b = all(r["status"] == "solved" for r in log4) and g4.d == {(4, 0): a40, (4, 1): -2}
    check("planted", p1 and p2 and p3a and p3b and p4a and p4b,
          f"P1 F'' + 2F'/eta = 0 gives exactly F_inf + c1/eta through eta^-8, no other coefficient: {p1}; P2 the terms "
          f"c (ln eta)^j, j = 1..3, are rejected: {p2}; P3 F'' + F'/eta = 0 accepts one c ln eta (a0_1 free, j = 2, 3 zero): "
          f"{p3a}, and c1/eta is then an obstruction (c1 = 0): {p3b}; P4 q + 4p = 8 nu^4: the pure series fails at "
          f"nu^4: {p4a}, with logs g = (a4_0 - 2 T) nu^4 and nothing else through nu^8: {p4b}")


# ---------------------------------------------------------------------------------------------------------------
# C4 Newton polygon, C5 the other starts
# ---------------------------------------------------------------------------------------------------------------
def N_for_polygon(fi):
    """N in (G, p, q, nu), with the mutations of the polygon controls."""
    Nt = N_theta()
    if MUT == "poly-drop":
        Nt = sp.expand(Nt + p**2 * q**2)
    if MUT == "poly0-edge" and fi == 0:
        Nt = sp.expand(Nt + nu**4 * q**2)
    return Nt


def polygon(Nt, fi):
    """Points (d, m) of N(fi + g, p, q, nu) in the monomials g^i p^b q^c nu^m, d = i + b + c, each with its polynomial in
    (k, c) for g = c nu^k, p -> -k c, q -> k^2 c. Returns (points, edges, vertex roots) for exponents k > 0."""
    g, k, c = sp.symbols("g k c")
    E = sp.expand(Nt.subs(G, fi + g))
    P = sp.Poly(E, g, p, q, nu)
    pts = {}
    for (i, b, cc, m), co in zip(P.monoms(), P.coeffs()):
        pts[(i + b + cc, m)] = pts.get((i + b + cc, m), 0) + co * c**i * (-k * c) ** b * (k**2 * c) ** cc
    pts = {pt: sp.factor(v) for pt, v in pts.items() if sp.expand(v) != 0}
    ks = set()
    for (d1, m1) in pts:
        for (d2, m2) in pts:
            if d1 < d2 and m1 > m2:
                ks.add(sp.Rational(m1 - m2, d2 - d1))
    edges, vertex = [], []
    for kk in sorted(ks):
        w = {pt: kk * pt[0] + pt[1] for pt in pts}
        wmin = min(w.values())
        S = [pt for pt in pts if w[pt] == wmin]
        if len(S) >= 2:
            ep = sp.factor(sum(pts[pt].subs(k, kk) for pt in S))
            rts = sp.roots(sp.Poly(sp.expand(ep), c))
            alone = {pt: sp.simplify(pts[pt].subs(k, kk) / c ** pt[0]) for pt in S}
            edges.append(dict(k=kk, points=S, poly=ep, roots=rts, alone=alone))
    for pt, v in pts.items():
        V = sp.Poly(sp.expand(sp.cancel(v / c ** pt[0])), k)
        rts = sp.roots(V) if V.degree() > 0 else {}
        assert sum(rts.values()) == V.degree(), f"roots of {V} not all found"
        for r in rts:
            if r.is_real and r > 0:
                w = {q_: r * q_[0] + q_[1] for q_ in pts}
                if w[pt] == min(w.values()):
                    vertex.append((pt, r))
    return pts, edges, vertex


def c4_polygon():
    out = []
    ok = True
    for label, fi in (("F_inf != 0", Fi), ("F_inf = 0", sp.Integer(0))):
        pts, edges, vertex = polygon(N_for_polygon(fi), fi)
        out.append(f"{label}: points (d, m) {sorted(pts)}, edges {[(str(e['k']), e['roots']) for e in edges]}, vertex "
                   f"balances with k > 0 {[(v[0], str(v[1])) for v in vertex]}")
        if fi != 0:
            e1 = len(edges) == 1 and edges[0]["k"] == 2
            if e1:
                e = edges[0]
                nz = {r for r in e["roots"] if r != 0}
                e1 = nz == {-2 * Fi} and e["roots"].get(0) == 3 and all(v != 0 for v in e["alone"].values())
            ok &= e1 and not vertex
        else:
            ok &= not edges and not vertex
    check("polygon", ok, "; ".join(out))


def N_for_a2(fi):
    Nt = N_theta()
    if MUT == "a2zero-edge":
        Nt = sp.expand(Nt + 2 * G * nu**2 * q**3 - 2 * G * nu**3 * q**3)
    return Nt


def c5_a2_zero():
    ok = True
    notes = []
    for label, fi in (("F_inf != 0", Fi), ("F_inf = 0", sp.Integer(0))):
        terms = poly_terms(N_for_a2(fi))
        for k in range(1, 9):
            us = [sp.Symbol(f"u{j}", real=True) for j in range(3)]
            g = Ser({(k, j): us[j] for j in range(3)})
            E = eval_terms(terms, fi, g, 4 * k + 4)
            o = E.low()
            eqs = [sp.together(c) for c in sp.Poly(sp.expand(E.coeffT(o)), T).all_coeffs()]
            sols = sp.solve(eqs, us, dict=True)
            sols = [{u: s.get(u, u) for u in us} for s in sols]
            nontrivial = [s for s in sols if any(v != 0 for v in s.values())]
            want = [{us[0]: -2 * Fi, us[1]: 0, us[2]: 0}] if (fi != 0 and k == 2) else []
            good = len(nontrivial) == len(want) and all(
                any(all(sp.simplify(s[u] - w[u]) == 0 for u in us) for w in want) for s in nontrivial)
            ok &= good
            if not good or (fi != 0 and k in (1, 2, 3)) or (fi == 0 and k in (1, 4)):
                notes.append(f"{label}, k = {k}: lowest order nu^{o}, nonzero solutions {[tuple(s.values()) for s in nontrivial]}")
    check("a2_zero", ok, "g = nu^k (u0 + u1 T + u2 T^2) from g = 0, k = 1..8, F_inf symbolic nonzero and F_inf = 0: only u = 0 except "
          "(F_inf != 0, k = 2) with u = (-2 F_inf, 0, 0); " + "; ".join(notes))


# ---------------------------------------------------------------------------------------------------------------
# C6, C7 the series
# ---------------------------------------------------------------------------------------------------------------
KMAX = 9
SERIES = {}


def series_log_run(jmax_of):
    return formal_series(N_theta(), Fi, 2, -2 * Fi, KMAX, jmax_of)


def get_series():
    if "g" not in SERIES:
        SERIES["g"], SERIES["log"] = series_log_run(jlog)
    return SERIES["g"], SERIES["log"]


def c6_series_pure():
    g, log = formal_series(N_theta(), Fi, 2, -2 * Fi, 6, lambda k: 0)
    want = 8 if MUT == "series-8" else 7
    a3 = log[0]["status"] == "solved" and log[0]["order"] == 9 and log[0]["sol"][0] == 0
    ob = log[1]["status"] == "obstruction" and log[1]["order"] == 10
    ratio = sp.cancel(log[1]["cond"] / (Fi**4 * (Fi**2 + want))) if ob else None
    cond_ok = ob and ratio == -512
    check("series_pure", a3 and cond_ok,
          f"pure series F_inf - 2 F_inf/eta^2 + ...: a_3 = 0 at eta^-9: {a3}; at eta^-10 no unknown is left and the condition "
          f"is {log[1].get('cond')} = 0, which is {ratio} times F_inf^4 (F_inf^2 + {want}): {cond_ok}")


def c7_series_log():
    jf = (lambda k: 0) if MUT == "series-nolog" else jlog
    g, log = series_log_run(jf) if MUT == "series-nolog" else get_series()
    solved = all(r["status"] == "solved" for r in log) and len(log) == KMAX - 2
    free = [f for r in log for f in r.get("free", [])]
    one_free = free == ["a4_0"]
    a41 = g.d.get((4, 1))
    b_ok = a41 is not None and sp.simplify(a41 - Fi * (Fi**2 + 7)) == 0 and not sp.real_roots(sp.Poly(Fi**2 + 7, Fi))
    odd = all(v == 0 for (n, j), v in g.d.items() if n % 2 == 1)
    res = residual_orders(N_theta(), Fi, g, KMAX + 6) if solved else None
    res_ok = res == []
    a4 = sp.Symbol("a4_0", real=True)
    check("series_log", solved and one_free and b_ok and odd and res_ok,
          f"with logs every stage k = 3..{KMAX} is solvable: {solved}; free parameters {free}; a_{{4,1}} = F_inf (F_inf^2 + 7), "
          f"nonzero for every real F_inf != 0: {b_ok}; odd coefficients zero: {odd}; the series through eta^-{KMAX} leaves N "
          f"nonzero only from nu^{KMAX + 7} on (orders <= nu^{KMAX + 6} that do not vanish: {res})")
    print("  series (g = F - F_inf, nu = 1/eta, T = ln eta):")
    for (n, j), v in sorted(g.d.items()):
        print(f"    nu^{n} T^{j}: {sp.factor(v)}")


# ---------------------------------------------------------------------------------------------------------------
# C8 reduced system, C9 linear modes
# ---------------------------------------------------------------------------------------------------------------
def c8_reduced_system():
    r, s, sig = sp.symbols("r s sigma")
    Nt = N_theta()
    # p = nu^2 r, q = nu^2 s with s = r' - 2r, sigma = nu^2
    E = sp.expand(Nt.subs({p: nu**2 * r, q: nu**2 * s}))
    Pn = sp.Poly(E, nu)
    low = min(m[0] for m in Pn.monoms())
    R = sp.expand(sp.cancel(E / nu**low)).subs(nu, sp.sqrt(sig))
    R = sp.expand(R)
    R0 = R.subs(sig, 0)
    lead_ok = low == 8 and zero(R0 + s**2 * (r**2 + 2 * G * s))
    # solve R = 0 for s near (s, r, sigma) = (-8 G, 4 G, 0): IFT
    pt = {s: -8 * G, r: 4 * G, sig: 0}
    dRs = sp.simplify(sp.diff(R, s).subs(pt))
    ift_ok = sp.simplify(R.subs(pt)) == 0 and dRs != 0
    # Jacobian of (G' = sigma r, r' = 2 r + s(G, r, sigma), sigma' = -2 sigma) at the equilibrium
    dsdr = sp.simplify(-sp.diff(R, r).subs(pt) / dRs)
    dsdG = sp.simplify(-sp.diff(R, G).subs(pt) / dRs)
    dsds = sp.simplify(-sp.diff(R, sig).subs(pt) / dRs)
    J = sp.Matrix([[0, sig, r], [dsdG, 2 + dsdr, dsds], [0, 0, -2]]).subs({sig: 0, r: 4 * G})
    eig = J.eigenvals()
    want = {0: 1, -2: 2} if MUT != "reduced-eig" else {0: 1, 2: 1, -2: 1}
    eq_ok = eig == want and sp.simplify(2 * 4 * G + (-8 * G)) == 0
    check("reduced_system", lead_ok and ift_ok and eq_ok,
          f"N/nu^8 = -s^2 (r^2 + 2 G s) + O(nu^2) with s = r' - 2r: {lead_ok}; at (r, s, nu^2) = (4G, -8G, 0) the equation "
          f"holds and its derivative in s is {dRs} != 0, so r' = 2r + s(G, r, nu^2) is analytic there: {ift_ok}; the line "
          f"r = 4G, nu^2 = 0 is a line of equilibria of (G' = nu^2 r, r', (nu^2)' = -2 nu^2) with eigenvalues {eig} "
          f"(d s/d r = {dsdr}, d s/d G = {dsdG}): {eq_ok}")


def c9_linear_modes():
    g, log = get_series()
    a2 = -sp.Rational(11, 5) * Fi if MUT == "modes-a2" else -2 * Fi
    g = Ser({**g.d, (2, 0): a2})
    Nt = N_theta()
    k = sp.Symbol("k")
    lead = {}
    for name, var in (("A", q), ("B", p), ("C", G)):
        terms = poly_terms(sp.diff(Nt, var))
        E = eval_terms(terms, Fi, g, 12)
        o = E.low()
        lead[name] = (o, sp.factor(E.coeffT(o)))
    oA, A = lead["A"]
    oB, B = lead["B"]
    oC, C = lead["C"]
    ind = sp.factor(A * k**2 - B * k) if oA == oB else None
    roots = sp.solve(ind, k) if ind is not None else None
    ok = ind is not None and sorted(roots) == [0, 4] and A != 0 and oC > oA
    check("linear_modes", ok,
          f"dN/dq = {A} nu^{oA}, dN/dp = {B} nu^{oB}, dN/dG is of order nu^{oC}, smaller; for delta ~ nu^k the leading "
          f"characteristic polynomial A k^2 - B k = {ind} has the roots {roots} (k = 0 shifts F_inf, k = 4 shifts a_4); "
          f"A != 0, so the linearization is a regular second-order equation with no exp(-kappa eta) modes: {ok}")


# ---------------------------------------------------------------------------------------------------------------
# C10, C11 numerics
# ---------------------------------------------------------------------------------------------------------------
mp.mp.dps = 30
ODE_TOL, ODE_DEG = mp.mpf(10) ** -20, 10  # odefun raises the working precision by (1 + degree)
NUM = {}


def num_setup():
    if NUM:
        return NUM
    g, _ = get_series()
    C = sp.Symbol("a4_0", real=True)
    d = dict(g.d)
    if MUT == "num-nolog":
        d.pop((4, 1))
    syms = (Fi, C)
    nus, Ts = sp.symbols("nus Ts")
    def build(kmax):
        e = Fi + sum(v * nus**n * Ts**j for (n, j), v in d.items() if n <= kmax)
        pe = sum(v * (-n * nus**n * Ts**j + (j * nus**n * Ts**(j - 1) if j else 0)) for (n, j), v in d.items() if n <= kmax)
        return (sp.lambdify((Fi, C, nus, Ts), e, "mpmath"), sp.lambdify((Fi, C, nus, Ts), pe, "mpmath"))
    NUM["G9"], NUM["p9"] = build(KMAX)
    NUM["G7"], NUM["p7"] = build(KMAX - 2)  # without the nu^8 terms (the odd coefficients vanish)
    Nt = N_theta()
    NUM["N"] = sp.lambdify((G, p, q, nu), Nt, "mpmath")
    NUM["dN"] = sp.lambdify((G, p, q, nu), sp.diff(Nt, q), "mpmath")
    return NUM


def q_root(Gv, pv, n):
    """Root of N(G, p, q, nu) = 0 in q near the leading branch q = -p^2/(2 G nu^2), by Newton."""
    Nf, dN = NUM["N"], NUM["dN"]
    q0 = -pv**2 / (2 * Gv * n**2)
    x = q0
    tol = mp.mpf(10) ** (-mp.mp.dps + 8)
    for _ in range(100):
        dx = Nf(Gv, pv, x, n) / dN(Gv, pv, x, n)
        x -= dx
        if abs(dx) <= tol * abs(x):
            if abs(x - q0) > abs(q0) / 2:
                raise RuntimeError("left the leading branch")
            return x
    raise RuntimeError("Newton did not converge")


def rhs(t, y):
    return [y[1], q_root(y[0], y[1], mp.e ** (-t))]


def c10_num_inward():
    n_ = num_setup()
    cases = [(mp.mpf(1), mp.mpf(0)), (mp.mpf(-1), mp.mpf(2)), (mp.mpf(1) / 3, mp.mpf(-1)), (mp.mpf(2), mp.mpf(1) / 2)]
    e1 = mp.mpf(1000)
    t1 = mp.log(e1)
    rows, ok = [], True
    for Fv, Cv in cases:
        nv = 1 / e1
        y0 = [n_["G9"](Fv, Cv, nv, t1), n_["p9"](Fv, Cv, nv, t1)]
        sol = mp.odefun(lambda tau, y: [-v for v in rhs(t1 - tau, y)], 0, y0, tol=ODE_TOL, degree=ODE_DEG)
        for e2 in (200, 100):
            t2 = mp.log(e2)
            y = sol(t1 - t2)
            nv2 = mp.mpf(1) / e2
            dG = abs(y[0] - n_["G9"](Fv, Cv, nv2, t2))
            dp = abs(y[1] - n_["p9"](Fv, Cv, nv2, t2))
            spread = abs(n_["G9"](Fv, Cv, nv2, t2) - n_["G7"](Fv, Cv, nv2, t2))
            good = dG <= mp.mpf("1e-10") and dp <= mp.mpf("1e-10")
            ok &= good
            rows.append(f"(F_inf, a_4) = ({mp.nstr(Fv, 4)}, {mp.nstr(Cv, 4)}), eta_2 = {e2}: |dG| = {mp.nstr(dG, 3)}, |dp| = "
                        f"{mp.nstr(dp, 3)}, size of the nu^8 terms {mp.nstr(spread, 3)}, G = {mp.nstr(y[0], 12)}: "
                        f"{'ok' if good else 'FAIL'}")
    check("num_inward", ok, "N = 0 integrated inward from eta_1 = 1000 (series data through eta^-9) agrees with the K = 9 series "
          "to 1e-10 in G and p. " + " | ".join(rows))


def c11_num_outward():
    n_ = num_setup()
    Fv, Cv = mp.mpf(1), mp.mpf(0)
    e0 = mp.mpf(100)
    t0 = mp.log(e0)
    y0 = [n_["G9"](Fv, Cv, 1 / e0, t0) + mp.mpf(1) / 100, n_["p9"](Fv, Cv, 1 / e0, t0) * mp.mpf(6) / 5]
    sol = mp.odefun(rhs, t0, y0, tol=ODE_TOL, degree=ODE_DEG)
    e_end = mp.mpf(10) ** 5
    yend = sol(mp.log(e_end))
    ne = 1 / e_end
    te = mp.log(e_end)
    fit = mp.findroot([lambda a, b: n_["G9"](a, b, ne, te) - yend[0], lambda a, b: n_["p9"](a, b, ne, te) - yend[1]],
                      (yend[0], mp.mpf(0)))
    Fp, Cp = fit[0], fit[1]
    rows, ok = [], True
    for e2 in (1000, 10000):
        t2 = mp.log(e2)
        y = sol(t2)
        nv = mp.mpf(1) / e2
        dG = abs(y[0] - n_["G9"](Fp, Cp, nv, t2))
        dp = abs(y[1] - n_["p9"](Fp, Cp, nv, t2))
        good = dG <= mp.mpf("1e-10") and dp <= mp.mpf("1e-10")
        ok &= good
        rows.append(f"eta = {e2}: |dG| = {mp.nstr(dG, 3)}, |dp| = {mp.nstr(dp, 3)}: {'ok' if good else 'FAIL'}")
    moved = abs(Fp - 1) > mp.mpf("1e-3")
    check("num_outward", ok and moved,
          f"from eta_0 = 100 with G + 1/100 and p x 6/5 of the (1, 0) series, N = 0 integrated outward to eta = 1e5 tends to the "
          f"limit F_inf' = {mp.nstr(Fp, 12)} (a_4' = {mp.nstr(Cp, 8)}), different from 1: {moved}; at eta = 1e3, 1e4 it equals "
          f"the K = 9 series with these two parameters to 1e-10: " + " | ".join(rows))


# ---------------------------------------------------------------------------------------------------------------
# C12 far field
# ---------------------------------------------------------------------------------------------------------------
def c12_farfield_scale():
    x, rho = sp.symbols("x rho", positive=True)
    ep = sp.Symbol("eps")
    U = sp.Symbol("U", real=True, nonzero=True)
    w = -2 * U * rho**2 / x**2
    u = U + ep * w
    D = sp.together(D_static(u, x, rho))
    num, den = sp.fraction(D)
    Pe = sp.Poly(sp.expand(num), ep)
    c = {k: sp.factor(Pe.coeff_monomial(ep**k) / den) for k in range(0, 8)}
    low_zero = all(c[k] == 0 for k in range(3))
    wxx, wxr = sp.diff(w, x, 2), sp.diff(w, x, rho)
    dperp = sp.diff(w, rho, 2) + sp.diff(w, rho) / rho
    sec7 = sp.simplify(c[3] - (-U * (wxr / 2) ** 2 * ((1 - U**2) * wxx + dperp))) == 0
    lim = {k: sp.limit(sp.simplify(c[k] * x**8 / rho**2), x, sp.oo) for k in (3, 4, 5, 6)}
    c4 = 0 if MUT == "ff-drop" else lim[4]
    ok = low_zero and sec7 and lim[3] != 0 and c4 != 0 and sp.simplify(lim[3] + c4) == 0 and lim[5] == 0 and lim[6] == 0
    check("farfield_scale", ok,
          f"u = U + eps w, w = -2 U rho^2/x^2: the eps^0, eps^1, eps^2 coefficients of D vanish: {low_zero}; the eps^3 coefficient "
          f"is section 7's -U Z^2[(1 - U^2) w_xx + Delta_perp w]: {sec7}; at fixed rho and x -> infinity (x^8/rho^2) times the "
          f"eps^3, eps^4, eps^5, eps^6 coefficients tends to {lim[3]}, {lim[4]}, {lim[5]}, {lim[6]}, so the eps^3 and eps^4 "
          f"terms are of the same order in 1/x and cancel at that order: {ok}")


CHECKS = [("n_build", c1_n_build), ("theta_parity", c2_theta_parity), ("planted", c3_planted), ("polygon", c4_polygon),
          ("a2_zero", c5_a2_zero), ("series_pure", c6_series_pure), ("series_log", c7_series_log),
          ("reduced_system", c8_reduced_system), ("linear_modes", c9_linear_modes), ("num_inward", c10_num_inward),
          ("num_outward", c11_num_outward), ("farfield_scale", c12_farfield_scale)]


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
            p_ = subprocess.run(args, capture_output=True, text=True)
            if "logdir" in OPTS:
                with open(os.path.join(OPTS["logdir"], f"ctl_{name}.log"), "w") as fh:
                    fh.write(p_.stdout + p_.stderr)
                with open(os.path.join(OPTS["logdir"], f"ctl_{name}.rc"), "w") as fh:
                    fh.write(f"{p_.returncode}\n")
            fails = [ln.split(":")[0][7:] for ln in p_.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = all(w in fails for w in want) and p_.returncode == 1
            extra = [f for f in fails if f not in want and f not in failed]
            controls_ok &= fired and not extra
            print(f"[{'FIRED' if fired and not extra else 'SILENT' if not fired else 'EXTRA'}] --mutate={name} ({what}): "
                  f"exit {p_.returncode}, failed {fails}", flush=True)
            if p_.returncode not in (0, 1):
                print(p_.stderr[-2000:])
        print(f"controls: {'all fired on target only' if controls_ok else 'SOME DID NOT FIRE ON TARGET ONLY'}")
    sys.exit(0 if not failed and controls_ok else 1)

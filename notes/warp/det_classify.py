#!/usr/bin/env python3
"""Which shifts b(t, x, rho) have det R_3 = 0 on an open set, and does any of them localize in x? (issue #80)

Metric ds^2 = -dt^2 + (dx + b dt)^2 + dy^2 + dz^2, rho^2 = y^2 + z^2. D = det R_3 on the frame (n, d_x, e_rho), with
n = d_t - b d_x, is #65's D = h (h^2 - P^2 - Q^2 + Y^2 + Z^2) - P (Y^2 - Z^2) + 2 Q Y Z. By #65's lemma D = 0 is the
algebraic stage of the normal + phantom scalar test, with lambda = 0, for every axisymmetric shift. This script
  C1  re-derives the block of R_3 for a generic b, with a = b_t - b b_x: P = a_x, Y = (b_x b_rho - a_rho)/2;
  C2  reduces the traveling shift B(x - v t, rho) to the static shift u = B + v;
  C3  the planar shifts b(t, x): D = 0 identically, and the algebraic stage at exact points;
  C4  closure for planar shifts: d_p d_q ln|R_pq| = 0, and Theta'/2 for a static one;
  C5  the inputs of the theorem for traveling planar shifts;
  C6  the family u = c x + g(rho) (Bessel), which contains #65's b(rho);
  C7  additive separable b = X(x) + R(rho): the classification;
  C8  multiplicative b = X(xi) R(rho) + v: the reduced equation and its first compatibility relation;
  C9  the conical shifts F(x/rho): the ODE and its parts of degree 3 and 4;
  C10 the far field of a localized shift;
  C11 the constant-speed bubble, as planted truth: no open set.
Imports: phantom_source (ricci, ricci_at_point, TX, y, z) and two_scalars (closure, dform, two_scalar_form).
Usage: det_classify.py [--no-controls] [--mutate=NAME] [--only=check1,check2]. Exit 1 if a check fails or a control
does not fire."""
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
    "d-sign": (["d_block"], "builds Y as (b_x b_rho + a_rho)/2"),
    "travel-offset": (["travel"], "takes the comoving shift u = B, without + v"),
    "planar-leak": (["planar"], "adds rho^2/10 to the planar b in the symbolic check D = 0"),
    "planar-point": (["planar"], "adds rho^2/10 to the planar shifts in the 'form' check"),
    "closure-deform": (["planar_closure"], "adds p q^2/10 to sigma of the planted passes"),
    "theta-factor": (["planar_closure"], "expects Theta'/4 in place of Theta'/2"),
    "planted-bubble": (["planar_closure"], "replaces the bubble of (iii) by the jet with Theta'(0) = 0"),
    "pg-null": (["planar_traveling"], "takes A' = 1/(1 + u) in the null coordinates"),
    "bessel-freq": (["bessel"], "uses J0(c rho), Y0(c rho) in place of J0(sqrt2 c rho), Y0(sqrt2 c rho)"),
    "bessel-closed": (["bessel"], "replaces (a, b) by the closed pair (d rho, dt + dx)"),
    "add-drop": (["additive"], "replaces (X + R)^2 by X^2 in E"),
    "add-quad": (["additive"], "drops the y^2 coefficient equation of the elimination step"),
    "mult-drop-v": (["multiplicative"], "builds u = X R without + v in D"),
    "conical-term": (["conical"], "changes the coefficient of F^3 F''^3 eta^2 in N from 2 to 3"),
    "conical-hand": (["conical_parts"], "expects the original wrong degree-3 part 2 F w^2 (F^2 F'' - eta w)"),
    "conical-plant": (["conical_planted"], "replaces the planted F'' + 2 F'/eta by F'' + F'/eta (no limit solution)"),
    "far-sign": (["farfield"], "expects + u_inf in front of the e^3 coefficient"),
    "bubble-w2": (["bubble"], "drops the shift's square u u_xx from P"),
    "bubble-ode": (["bubble"], "takes r w'' = +w' in place of r w'' = -w' in the ODE step"),
}

RESULTS = []
OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate")
ONLY = OPTS["only"].split(",") if "only" in OPTS else None


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def wanted(name):
    return ONLY is None or name in ONLY


# ---------------------------------------------------------------------------------------------------------------
# D from jets
# ---------------------------------------------------------------------------------------------------------------
def blocks_jets(b, bx, br, bt, btx, bxx, bxr, btr, brr, rho):
    """#65's P, Q, Z, Y, h from the jets of b(t, x, rho), n = d_t - b d_x."""
    P = btx - b * bxx - bx**2
    Q = -(brr + br / rho) / 2
    Z = bxr / 2
    Y = -(btr - b * bxr) / 2 + bx * br
    return P, Q, Z, Y, br**2 / 2


def D_of(P, Q, Z, Y, h):
    return h * (h**2 - P**2 - Q**2 + Y**2 + Z**2) - P * (Y**2 - Z**2) + 2 * Q * Y * Z


def D_static(u, ux, ur, uxx, uxr, urr, rho):
    return D_of(*blocks_jets(u, ux, ur, 0, 0, uxx, uxr, 0, urr, rho))


def D_explicit(b, t_, x_, r_):
    d = sp.diff
    return D_of(*blocks_jets(b, d(b, x_), d(b, r_), d(b, t_), d(b, t_, x_), d(b, x_, 2), d(b, x_, r_), d(b, t_, r_),
                             d(b, r_, 2), r_))


def zero(e):
    if isinstance(e, sp.MatrixBase):
        return all(sp.simplify(x) == 0 for x in e)
    return sp.simplify(e) == 0


def theta_prime(u, xi):
    """Theta' for a static planar shift u(xi): Theta = (1 - u^2) Phi'/Phi, Phi = (u u')' (1 - u^2)."""
    K = sp.diff(u * sp.diff(u, xi), xi)
    Phi = K * (1 - u**2)
    return sp.diff((1 - u**2) * sp.diff(Phi, xi) / Phi, xi)


# ---------------------------------------------------------------------------------------------------------------
# C1 d_block
# ---------------------------------------------------------------------------------------------------------------
def c1_d_block():
    t, x, rho, phi = sp.symbols("t x rho phi", positive=True)
    b = sp.Function("b")(t, x, rho)
    g = sp.Matrix([[-1 + b**2, b, 0, 0], [b, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, rho**2]])
    X4 = [t, x, rho, phi]
    Rc = ps.ricci(g, X4, lambda E, k: sp.diff(E, X4[k]), sp.simplify)
    E = sp.Matrix([[1, 0, 0, 0], [-b, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1 / rho]])
    Rf = (E.T * Rc * E).applyfunc(sp.simplify)
    a = sp.diff(b, t) - b * sp.diff(b, x)
    h = sp.diff(b, rho)**2 / 2
    Pn = sp.diff(a, x)
    sgn = 1 if MUT == "d-sign" else -1
    Yn = (sp.diff(b, x) * sp.diff(b, rho) + sgn * sp.diff(a, rho)) / 2
    Qn = -(sp.diff(b, rho, 2) + sp.diff(b, rho) / rho) / 2
    Zn = sp.diff(b, x, rho) / 2
    P65, Q65, Z65, Y65, h65 = blocks_jets(b, sp.diff(b, x), sp.diff(b, rho), sp.diff(b, t), sp.diff(b, t, x),
                                          sp.diff(b, x, 2), sp.diff(b, x, rho), sp.diff(b, t, rho),
                                          sp.diff(b, rho, 2), rho)
    diffs = {"R00": Rf[0, 0] - (Pn - h), "R11": Rf[1, 1] - (-Pn - h), "R22": Rf[2, 2] - h, "R01": Rf[0, 1] - Qn,
             "R02": Rf[0, 2] - Zn, "R12": Rf[1, 2] - Yn, "Rphi": sum(abs(Rf[3, j]) for j in range(4)),
             "P=a_x vs #65": Pn - P65, "Y vs #65": Yn - Y65}
    blk = sp.Matrix([[Pn - h, Qn, Zn], [Qn, -Pn - h, Yn], [Zn, Yn, h]])
    diffs["det vs D"] = blk.det() - D_of(P65, Q65, Z65, Y65, h65)
    bad = [k for k, e in diffs.items() if not zero(e)]
    check("d_block", not bad, f"coordinate Ricci of a generic b(t, x, rho), frame (n, d_x, d_rho, d_phi/rho): block = "
          f"[[P-h, Q, Z], [Q, -P-h, Y], [Z, Y, h]], P = a_x, Y = (b_x b_rho - a_rho)/2, a = b_t - b b_x, R(e_phi, .) = 0, "
          f"det = D; nonzero differences: {bad or 'none'}")


# ---------------------------------------------------------------------------------------------------------------
# C2 travel
# ---------------------------------------------------------------------------------------------------------------
def c2_travel():
    B0, Bx, Br, Bxx, Bxr, Brr, v, rho, c = sp.symbols("B0 Bx Br Bxx Bxr Brr v rho c")
    Dlab = D_of(*blocks_jets(B0, Bx, Br, -v * Bx, -v * Bxx, Bxx, Bxr, -v * Bxr, Brr, rho))
    u = B0 if MUT == "travel-offset" else B0 + v
    Dst = D_static(u, Bx, Br, Bxx, Bxr, Brr, rho)
    d1 = sp.expand(Dlab - Dst)
    b0, bx, br, bt, btx, bxx, bxr, btr, brr = sp.symbols("b0 bx br bt btx bxx bxr btr brr")
    D0 = D_of(*blocks_jets(b0, bx, br, bt, btx, bxx, bxr, btr, brr, rho))
    D1 = D_of(*blocks_jets(b0 + c, bx, br, bt + c * bx, btx + c * bxx, bxx, bxr, btr + c * bxr, brr, rho))
    d2 = sp.expand(D1 - D0)
    check("travel", d1 == 0 and d2 == 0, f"D of b = B(x - v t, rho) minus D of the static u = B + v: {d1}; D of "
          f"b(t, x + c t, rho) + c minus D of b: {d2}. So a traveling shift is a static one, with u -> v at infinity for "
          f"a bubble")


# ---------------------------------------------------------------------------------------------------------------
# C3 planar
# ---------------------------------------------------------------------------------------------------------------
def c3_planar():
    b0, bx, bt, btx, bxx, rho = sp.symbols("b0 bx bt btx bxx rho")
    zero_ = sp.Integer(0)
    if MUT == "planar-leak":
        br, brr, bb = rho / 5, sp.Rational(1, 5), b0 + rho**2 / 10
    else:
        br, brr, bb = zero_, zero_, b0
    P, Q, Z, Y, h = blocks_jets(bb, bx, br, bt, btx, bxx, zero_, zero_, brr, rho)
    R3 = sp.Matrix([[P - h, Q, Z], [Q, -P - h, Y], [Z, Y, h]])
    sym_ok = zero(R3 - sp.diag(P, -P, 0)) and zero(R3.det())
    # exact points
    tx = ps.TX
    t_, x_ = tx
    shifts = [(-2 / (1 + (x_ - t_)**2), [(0, 1), (sp.Rational(1, 2), sp.Rational(3, 2)), (-1, sp.Rational(1, 3))]),
              ((t_**2 + 2 * x_) / (1 + x_**2), [(1, sp.Rational(1, 2)), (sp.Rational(-1, 2), 2),
                                                 (2, sp.Rational(-1, 3))]),
              (x_ / (1 + t_**2), [(sp.Rational(1, 3), 1), (2, -1)])]
    yv, zv = sp.Rational(3, 7), sp.Rational(4, 7)
    rv = sp.Rational(5, 7)
    eta = sp.diag(-1, 1, 1, 1)
    rr = sp.Symbol("rr", positive=True)
    forms, algs = [], []
    for bex, pts in shifts:
        bfull = bex + ((ps.y**2 + ps.z**2) / 10 if MUT == "planar-point" else 0)
        bleak = bex + rr**2 / 10
        for (tv, xv) in pts:
            Pv = sp.simplify((sp.diff(bex, t_, x_) - bex * sp.diff(bex, x_, 2) - sp.diff(bex, x_)**2)
                             .subs({t_: tv, x_: xv}))
            res = ts.two_scalar_form(ps.ricci_at_point(bfull, (tv, xv, yv, zv), False), eta)
            forms.append(Pv != 0 and res["status"] == "form" and res["lam"] == 0)
            res2 = ts.two_scalar_form(ps.ricci_at_point(bex + (ps.y**2 + ps.z**2) / 10, (tv, xv, yv, zv), False), eta)
            Dv = D_explicit(bleak, t_, x_, rr).subs({t_: tv, x_: xv, rr: rv})
            algs.append(res2["status"] == "algebraic" and sp.simplify(Dv) != 0)
    check("planar", sym_ok and all(forms) and all(algs), f"b(t, x): R_3 = diag(P, -P, 0), D = 0 symbolically: {sym_ok}; "
          f"{sum(forms)} of {len(forms)} exact points of 3 planar shifts have P != 0 and pass the algebraic stage with "
          f"lam = 0; {sum(algs)} of {len(algs)} points of the same shifts + rho^2/10 fail it with D != 0")


# ---------------------------------------------------------------------------------------------------------------
# C4 planar_closure
# ---------------------------------------------------------------------------------------------------------------
def null_metric_closure(sigma):
    p, q, y, z = sp.symbols("p q y z", real=True)
    X = [p, q, y, z]
    s = sigma(p, q)
    g = sp.Matrix([[0, -sp.exp(2 * s), 0, 0], [-sp.exp(2 * s), 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    R = ps.ricci(g, X, lambda E, k: sp.diff(E, X[k]), sp.simplify)
    Rpq = R[0, 1]
    k = sp.Matrix([1, 0, 0, 0])
    m = sp.Matrix([0, 2 * Rpq, 0, 0])
    a, b = (k + m) / 2, (k - m) / 2
    exact = zero((a * a.T - b * b.T) - R)
    ok, om, why = ts.closure(a, b, X, sp.simplify)
    crit = zero(sp.diff(sp.log(Rpq**2), p, q))
    return exact, ok, crit


def c4_planar_closure():
    deform = (lambda p, q: p * q**2 / 10) if MUT == "closure-deform" else (lambda p, q: 0)
    passes = [("sin(p) e^(q/2)/3", lambda p, q: sp.sin(p) * sp.exp(q / 2) / 3 + deform(p, q)),
              ("p q + p^2 q", lambda p, q: p * q + p**2 * q + deform(p, q)),
              ("(p + q)^2/5", lambda p, q: (p + q)**2 / 5 + deform(p, q))]
    fails = [("p^2 q/2 + p q^2/2", lambda p, q: p**2 * q / 2 + p * q**2 / 2),
             ("e^(p q)", lambda p, q: sp.exp(p * q)),
             ("p^3 q^2/6 + p q", lambda p, q: p**3 * q**2 / 6 + p * q)]
    okall, rows = True, []
    for lab, f in passes + fails:
        want = (lab, f) in passes
        exact, ok, crit = null_metric_closure(f)
        agree = exact and ok == want and crit == want
        okall &= agree
        rows.append(f"sigma = {lab}: R = (a a - b b) {exact}, closure {ok}, criterion {crit}, expected {want}")
    for r in rows:
        print("  " + r)
    # (ii) static planar shift
    t, xi, y, z = sp.symbols("t xi y z", real=True)
    u = sp.Function("u")(xi)
    K = sp.diff(u * sp.diff(u, xi), xi)
    P = -K
    th0, th1 = sp.Matrix([1, 0, 0, 0]), sp.Matrix([u, 1, 0, 0])
    k, m = th0 - th1, P * (th0 + th1)
    a, b = (k + m) / 2, (k - m) / 2
    ok2, om, why = ts.closure(a, b, [t, xi, y, z], sp.simplify)
    obstruction = ts.dform(om, [t, xi, y, z])[0][1]
    fac = sp.Rational(1, 4) if MUT == "theta-factor" else sp.Rational(1, 2)
    ratio = zero(obstruction - fac * theta_prime(u, xi))
    # (iii) jets
    u0, u1, u2, u3 = sp.Rational(1, 3), sp.Rational(1, 2), sp.Rational(1, 5), sp.Rational(2, 7)
    u4 = sp.Symbol("u4")
    poly = u0 + u1 * xi + u2 * xi**2 / 2 + u3 * xi**3 / 6 + u4 * xi**4 / 24
    tp0 = sp.simplify(theta_prime(poly, xi).subs(xi, 0))
    u4v = sp.solve(tp0, u4)
    planted = poly.subs(u4, u4v[0]) if u4v else None
    pass3 = False
    if planted is not None:
        Kp = sp.diff(planted * sp.diff(planted, xi), xi)
        th0, th1 = sp.Matrix([1, 0, 0, 0]), sp.Matrix([planted, 1, 0, 0])
        kk, mm = th0 - th1, -Kp * (th0 + th1)
        aa, bb = (kk + mm) / 2, (kk - mm) / 2
        _, omp, _ = ts.closure(aa, bb, [t, xi, y, z], sp.simplify)
        dw0 = sp.simplify(ts.dform(omp, [t, xi, y, z])[0][1].subs(xi, 0))
        pass3 = dw0 == 0 and sp.simplify(theta_prime(planted, xi).subs(xi, 0)) == 0
    v = sp.Rational(1, 2)
    ub = v * (1 - 1 / (1 + xi**2)**2)
    pts = [sp.Rational(1, 3), 1, sp.Rational(3, 2), 2, sp.Rational(5, 2)]
    if MUT == "planted-bubble" and planted is not None:
        ub, pts = planted, [0] * 5
    vals = [sp.simplify(theta_prime(ub, xi).subs(xi, pt)) for pt in pts]
    Phi_ok = all(sp.simplify((sp.diff(ub * sp.diff(ub, xi), xi) * (1 - ub**2)).subs(xi, pt)) != 0 for pt in pts)
    nz = Phi_ok and all(vv != 0 for vv in vals)
    check("planar_closure", okall and ratio and pass3 and nz,
          f"(i) 3 planted passes and 3 planted failures of the null-coordinate criterion d_p d_q ln|R_pq| = 0 agree with "
          f"`closure`: {okall}; (ii) d omega_(t xi) = {fac} Theta' for an undefined u: {ratio}; (iii) a jet with "
          f"Theta'(0) = 0 has d omega_(t xi)(0) = 0: {pass3}; the 1+1 bubble u = (1/2)(1 - 1/(1 + xi^2)^2) has "
          f"Theta' != 0 at {len(pts)} exact points: {nz}")


# ---------------------------------------------------------------------------------------------------------------
# C5 planar_traveling
# ---------------------------------------------------------------------------------------------------------------
def c5_planar_traveling():
    t, xi = sp.symbols("t xi", real=True)
    u = sp.Function("u")(xi)
    Ap = 1 / (1 + u) if MUT == "pg-null" else 1 / (1 - u)
    Bp = 1 / (1 + u)
    dtt = sp.Matrix([1, 0])
    # dp = dt - A' dxi, dq = dt + B' dxi ; -(1 - u^2) dp dq as a symmetric matrix
    dp, dq = sp.Matrix([1, -Ap]), sp.Matrix([1, Bp])
    gnull = -(1 - u**2) * (dp * dq.T + dq * dp.T) / 2
    gpg = sp.Matrix([[-1 + u**2, u], [u, 1]])
    a_ok = zero(gnull - gpg)
    Tt = [t, xi]
    R2 = ps.ricci(gpg, Tt, lambda E, k: sp.diff(E, Tt[k]), sp.simplify)
    K = sp.diff(u * sp.diff(u, xi), xi)
    b_ok = zero(R2 - K * gpg)
    p, q = sp.symbols("p q")
    Lam = sp.Function("Lam")
    c_ok = zero(sp.diff(Lam(q - p), p, q) + sp.diff(Lam(q - p), q, 2))
    check("planar_traveling", a_ok and b_ok and c_ok, f"for a static planar shift u(xi): the null coordinates p = t - A, "
          f"q = t + B, A' = 1/(1 - u), B' = 1/(1 + u) give -(1 - u^2) dp dq = the metric: {a_ok}; R_ab = (u u')' g_ab, "
          f"so R_pq = -Phi/2 with Phi = (u u')' (1 - u^2): {b_ok}; d_p d_q Lam(q - p) = -Lam'': {c_ok}")


# ---------------------------------------------------------------------------------------------------------------
# C6 bessel
# ---------------------------------------------------------------------------------------------------------------
def c6_bessel():
    t, x, rho, phi = sp.symbols("t x rho phi", positive=True)
    c = sp.Symbol("c", positive=True)
    g = sp.Function("g")(rho)
    g1, g2 = sp.diff(g, rho), sp.diff(g, rho, 2)
    S = g2 + g1 / rho
    Dfam = D_static(c * x + g, c, g1, 0, 0, g2, rho)
    expected = g1**2 * (2 * c**2 + g1**2 - S) * (2 * c**2 + g1**2 + S) / 8
    a_ok = zero(Dfam - expected)
    # (b) numerics
    mp.mp.dps = 130  # working precision: D is a difference of terms of size ~1, 60 correct digits need more
    random.seed(80)
    worst, worst_pert, npt = mp.mpf(0), mp.inf, 0
    where, pert_rows = None, []
    for cc in (sp.Rational(1, 2), sp.Rational(3, 2), 0):
        for eps in (1, -1):
            if cc == 0:
                C1, C2 = sp.Rational(random.randint(2, 9), 3), sp.Rational(random.randint(1, 5), 4)
                psi = C1 + C2 * sp.log(rho)
            else:
                al, be = sp.Rational(random.randint(1, 9), 3), sp.Rational(random.randint(-9, -1), 5)
                kap = cc if MUT == "bessel-freq" else sp.sqrt(2) * cc
                psi = al * sp.besselj(0, kap * rho) + be * sp.bessely(0, kap * rho)
            for pert in (False, True):
                ps_ = psi + (rho**2 / 10 if pert else 0)
                gg = -eps * sp.log(ps_)
                Dx = D_static(cc * x + gg, cc, sp.diff(gg, rho), 0, 0, sp.diff(gg, rho, 2), rho)
                Df = sp.lambdify(rho, Dx.subs(x, sp.Rational(1, 3)), modules="mpmath")
                for rv in (sp.Rational(3, 7), sp.Rational(5, 7), sp.Rational(9, 7), sp.Rational(13, 7),
                           sp.Rational(17, 7), sp.Rational(23, 7)):
                    if abs(sp.N(psi.subs(rho, rv), 30)) < sp.Rational(1, 1000):
                        continue
                    val = abs(Df(mp.mpf(rv.p) / rv.q))
                    if pert:
                        pert_rows.append((cc, eps, rv, val))
                        if val < worst_pert:
                            worst_pert, where = val, (cc, eps, rv)
                    else:
                        worst = max(worst, val)
                        npt += 1
    b_ok = worst < mp.mpf(10)**-50
    p_ok = worst_pert > mp.mpf(10)**-3
    print("  perturbed psi + rho^2/10, |D| at each point (c, eps, rho): " +
          "; ".join(f"({r[0]}, {r[1]:+d}, {r[2]}): {mp.nstr(r[3], 3)}" for r in pert_rows))
    print(f"  observation (not a criterion): unperturbed |D| < 1e-50 at all {npt} points; smallest perturbed |D| = "
          f"{mp.nstr(worst_pert, 3)} at c = {where[0]}, eps = {where[1]:+d}, rho = {where[2]}")
    # (c) R = a a - b b on the ODE (exact), for eps = +-1
    X4 = [t, x, rho, phi]
    u = c * x + g
    gm = sp.Matrix([[-1 + u**2, u, 0, 0], [u, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, rho**2]])
    Rc = ps.ricci(gm, X4, lambda E, k: sp.diff(E, X4[k]), sp.cancel)
    th0, th1, th2 = sp.Matrix([1, 0, 0, 0]), sp.Matrix([u, 1, 0, 0]), sp.Matrix([0, 0, 1, 0])
    c_ok = True
    for eps in (1, -1):
        rule2 = eps * (g1**2 + 2 * c**2) - g1 / rho
        a = (g1 * th2 + 2 * c * th1) / sp.sqrt(2)
        nu = th0 + eps * th1
        if MUT == "bessel-closed":
            a, nu = th2, th0 + th1
            bb = nu * nu.T
        else:
            bb = (g1**2 + 2 * c**2) / 2 * nu * nu.T
        res = (Rc - (a * a.T - bb)).applyfunc(lambda e: sp.simplify(e.subs(sp.Derivative(g, (rho, 2)), rule2)))
        c_ok &= all(e == 0 for e in res)
    # (d) closure at exact jets: the 2-form equations of the omega construction must have no solution
    G0, G1 = sp.symbols("G0 G1")
    rows, d_ok = [], True
    for cc, gv in ((2, 1), (6, 17), (0, sp.Rational(1, 2))):
        for eps in (1, -1):
            rule2 = eps * (g1**2 + 2 * cc**2) - g1 / rho
            rule3 = sp.diff(rule2, rho).subs(sp.Derivative(g, (rho, 2)), rule2)
            uu = cc * x + g
            th0, th1, th2 = sp.Matrix([1, 0, 0, 0]), sp.Matrix([uu, 1, 0, 0]), sp.Matrix([0, 0, 1, 0])
            W = sp.sqrt(g1**2 + 2 * cc**2)
            a = (g1 * th2 + 2 * cc * th1) / sp.sqrt(2)
            b = W * (th0 + eps * th1) / sp.sqrt(2)
            if MUT == "bessel-closed":
                a, b = th2, th0 + th1
            pt = {rho: sp.Rational(5, 7), x: sp.Rational(1, 3), t: sp.Rational(1, 10)}

            def simp(e, rule2=rule2, rule3=rule3, gv=gv, pt=pt):
                e = e.subs(sp.Derivative(g, (rho, 3)), rule3).subs(sp.Derivative(g, (rho, 2)), rule2)
                e = e.subs(sp.Derivative(g, rho), G1).subs(g, G0)
                return sp.simplify(e.subs({G1: gv, G0: sp.Rational(1, 5)}).subs(pt))
            ok, om, why = ts.closure(a, b, X4, simp)
            fails = (not ok) and why.startswith("no chi")
            d_ok &= fails
            rows.append((cc, eps, ok, why[:20]))
    print(f"  closure at exact jets (c, eps, ok, why): {rows}")
    check("bessel", a_ok and b_ok and c_ok and d_ok,
          f"(a) D = g'^2 (2c^2 + g'^2 - S)(2c^2 + g'^2 + S)/8 for u = c x + g(rho): {a_ok}; (b, first part) D = 0 to 1e-50 at "
          f"{npt} points (largest |D| {mp.nstr(worst, 3)}); (c) R = a a - b b on the ODE for eps = +-1: {c_ok}; (d) the omega "
          f"equations have no solution at 6 exact jets (c = 2, 6, 0): {d_ok}")
    check("bessel_pert", p_ok, f"(b, second part) for psi + rho^2/10, |D| > 1e-3 at the same points: smallest |D| = "
          f"{mp.nstr(worst_pert, 3)} at c = {where[0]}, eps = {where[1]:+d}, rho = {where[2]}")


# ---------------------------------------------------------------------------------------------------------------
# C7 additive
# ---------------------------------------------------------------------------------------------------------------
def c7_additive():
    X0, X1, X2, R0, R1, R2, rho, c = sp.symbols("X0 X1 X2 R0 R1 R2 rho c")
    Dadd = D_static(X0 + R0, X1, R1, X2, 0, R2, rho)
    h = R1**2 / 2
    S = R2 + R1 / rho
    shift = X0 if MUT == "add-drop" else X0 + R0
    E = (h + X1**2)**2 - shift**2 * X2**2 - S**2 / 4
    a_ok = sp.expand(Dadd - h * E) == 0
    A, B = R1**2 / 2, S**2 / 4
    f = [1, X1**2, X2**2, X0 * X2**2, X1**4 - X0**2 * X2**2]
    gg = [A**2 - B, 2 * A, -R0**2, -2 * R0, 1]
    Etrue = (h + X1**2)**2 - (X0 + R0)**2 * X2**2 - S**2 / 4
    b_ok = sp.expand(Etrue - sum(fi * gi for fi, gi in zip(f, gg))) == 0
    rr = sp.Symbol("rr")
    wr = sp.simplify(sp.wronskian([sp.Integer(1), sp.Function("R")(rr), sp.Function("R")(rr)**2], rr))
    w_ok = zero(wr - 2 * sp.diff(sp.Function("R")(rr), rr)**3)
    # elimination step
    y, e0, b5, g0, b3, d0, b4, m, s_, c1, c0, xx = sp.symbols("y e0 b5 g0 b3 d0 b4 m s c1 c0 xx")
    ident = sp.expand((y**2 - e0 - b5 * y) * (g0 + b3 * y) - (d0 + b4 * y)**2)
    co = sp.Poly(ident, y).all_coeffs()  # y^3, y^2, y^1, y^0
    eqs = [co[0], co[2], co[3]] if MUT == "add-quad" else list(co)
    # X'' = s m const (s = +-1), X = s m xx^2/2 + c1 xx + c0, y = X'^2: the slope of m^2 X against y is b4
    Xq = s_ * m * xx**2 / 2 + c1 * xx + c0
    yq = sp.diff(Xq, xx)**2
    slope = sp.simplify(sp.diff(m**2 * Xq, xx) / sp.diff(yq, xx))
    tm = sp.Symbol("tm")
    unit = []
    for sv in (1, -1):
        eqs2 = [e.subs(s_, sv) for e in eqs] + [b4 - slope.subs(s_, sv), g0 - m**2, b3, tm * m - 1]
        Gb = sp.groebner(eqs2, tm, e0, b5, g0, b3, d0, b4, m, order="grevlex")
        unit.append(list(Gb.exprs) == [1])
    e_ok = all(unit)
    # f5 affine in y fails for the quadratic X: the y^2 coefficient of X'^4 - X^2 X''^2 is 3/4 (information)
    # X'' = 0: E reduces to ((R'^2 + 2c^2)^2 - S^2)/4
    d_ok = zero(E.subs({X2: 0, X1: c}) - ((R1**2 + 2 * c**2)**2 - S**2) / 4)
    check("additive", a_ok and b_ok and w_ok and e_ok and d_ok,
          f"b = X + R: D = h [(h + X'^2)^2 - ((X + R) X'')^2 - Q^2]: {a_ok}; E = sum f_i g_i with the five f_i, g_i: {b_ok}; "
          f"Wronskian of (1, R, R^2) = 2 R'^3: {w_ok}; elimination step with X'' = s m, s = +-1: the equations and m t = 1 "
          f"generate the unit ideal ({unit}): {e_ok}; X'' = 0 gives ((R'^2 + 2c^2)^2 - S^2)/4: {d_ok}")


# ---------------------------------------------------------------------------------------------------------------
# C8 multiplicative
# ---------------------------------------------------------------------------------------------------------------
def E_mult_poly(X0, X1, X2, R0, R1, S, v):
    """8D/R'^2 for u = X(x) R(rho) + v, as a polynomial in the jets (X0, X1, X2) = (X, X', X''), (R0, R1) = (R, R') and
    S = R'' + R'/rho (derived in the scratch session, checked against D in C8 (a))."""
    return (-4 * R0**4 * X0**4 * X2**2 + 10 * R0**4 * X0**3 * X1**2 * X2 + 14 * R0**4 * X0**2 * X1**4
            - 8 * R0**3 * X0**3 * X2**2 * v + 22 * R0**3 * X0**2 * X1**2 * X2 * v + 12 * R0**3 * X0 * X1**4 * v
            + 9 * R0**2 * R1**2 * X0**4 * X1**2 - 4 * R0**2 * X0**2 * X2**2 * v**2
            + 14 * R0**2 * X0 * X1**2 * X2 * v**2 - 2 * R0**2 * X0 * X1**2 * X2 + 2 * R0**2 * X1**4 * v**2
            - 2 * R0**2 * X1**4 + 6 * R0 * R1**2 * X0**3 * X1**2 * v - 6 * R0 * S * X0**2 * X1**2
            + 2 * R0 * X1**2 * X2 * v**3 - 2 * R0 * X1**2 * X2 * v + R1**4 * X0**6 + R1**2 * X0**2 * X1**2 * v**2
            + R1**2 * X0**2 * X1**2 - S**2 * X0**4 - 2 * S * X0 * X1**2 * v)


def E_add_poly(X0, X1, X2, R0, R1, S):
    """E of C7 for u = X + R."""
    return (R1**2 / 2 + X1**2)**2 - (X0 + R0)**2 * X2**2 - S**2 / 4


def c8_multiplicative():
    X0, X1, X2, R0, R1, S, v, rho = sp.symbols("X0 X1 X2 R0 R1 S v rho")
    vv = 0 if MUT == "mult-drop-v" else v
    u = X0 * R0 + vv
    ux, ur, uxx, uxr = X1 * R0, X0 * R1, X2 * R0, X1 * R1
    # u_rho,rho + u_rho/rho = X0 S with S = R'' + R'/rho
    P, Q, Z, Y, h = blocks_jets(u, ux, ur, 0, 0, uxx, uxr, 0, X0 * S - ur / rho, rho)
    D = D_of(P, Q, Z, Y, h)
    Eprint = E_mult_poly(X0, X1, X2, R0, R1, S, v)
    a_ok = sp.expand(sp.cancel(8 * D / R1**2) - Eprint) == 0
    # (b) v = 0 in log derivatives
    p1, p2, c1, s2 = sp.symbols("p1 p2 c1 s2")
    T = 5 * p1**2 + c1**2 - 2 * p2
    A = T * (9 * p1**2 + 2 * c1**2 - T)
    Bq = p1**2 * T - (s2 + 3 * p1**2)**2
    sub = {X1: X0 * p1, X2: X0 * (p2 + p1**2), R1: R0 * c1, S: R0 * s2}
    E0 = sp.expand(Eprint.subs(v, 0).subs(sub))
    b_ok = sp.expand(E0 - R0**2 * X0**4 * (R0**2 * X0**2 * A + Bq)) == 0
    # (c) subcases
    c1_ok = sp.expand(D.subs({R1: 0, S: 0})) == 0
    Dx = sp.expand(sp.cancel(D.subs({X1: 0, X2: 0})))
    hh = X0**2 * R1**2 / 2
    Qq = -X0 * S / 2
    c2_ok = sp.expand(Dx - hh * (hh**2 - Qq**2)) == 0
    c_ok = c1_ok and c2_ok
    # (d) compatibility
    p3, s2p = sp.symbols("p3 s2p")
    Af = c1**4 + 9 * c1**2 * p1**2 + 20 * p1**4 + 2 * p1**2 * p2 - 4 * p2**2
    Bf = c1**2 * p1**2 - 4 * p1**4 - 2 * p1**2 * p2 - 6 * p1**2 * s2 - s2**2
    form_ok = sp.expand(Af - A) == 0 and sp.expand(Bf - Bq) == 0
    c2s = s2 - c1**2 - c1 / rho
    Ax = sp.diff(Af, p1) * p2 + sp.diff(Af, p2) * p3
    Bx = sp.diff(Bf, p1) * p2 + sp.diff(Bf, p2) * p3
    eqI = sp.expand(Af * Bx - Ax * Bf - 2 * p1 * Af * Bf)
    Ar = sp.diff(Af, c1) * c2s
    Br = sp.diff(Bf, c1) * c2s + sp.diff(Bf, s2) * s2p
    eqII = sp.expand(Af * Br - Ar * Bf - 2 * c1 * Af * Bf)
    p3s = sp.solve(eqI, p3)[0]
    s2s = sp.solve(eqII, s2p)[0]
    C1 = sp.together(sp.diff(p3s, c1) * c2s + sp.diff(p3s, s2) * s2s)
    C2 = sp.together(sp.diff(s2s, p1) * p2 + sp.diff(s2s, p2) * p3s)
    n1, n2 = sp.fraction(C1)[0], sp.fraction(C2)[0]
    d_ok = form_ok and sp.expand(n1 - n2) == 0 and n1 != 0
    degs = sp.Poly(sp.cancel(n1 / (2 * c1 * p1 * (c1**2 + 4 * p1**2 + 2 * p2) * (c1**2 + 5 * p1**2 - 2 * p2))),
                   p1, p2, c1, s2, rho).total_degree()
    check("multiplicative", a_ok and b_ok and c_ok and d_ok,
          f"u = X R + v: 8D/R'^2 = the printed polynomial: {a_ok}; v = 0: E0 = R^2 X^4 [X^2 R^2 A + B]: {b_ok}; "
          f"R' = 0 gives D = 0 and X' = 0 gives h (h^2 - Q^2): {c_ok}; the first compatibility conditions have equal "
          f"numerators (relation Gamma of degree {degs} after removing the factors c1 p1 T+ T-): {d_ok}. "
          f"A classification is NOT claimed")


# ---------------------------------------------------------------------------------------------------------------
# C9 conical
# ---------------------------------------------------------------------------------------------------------------
def c9_conical():
    x, rho = sp.symbols("x rho", positive=True)
    eta = ETA
    F = sp.Function("F")
    u = F(x / rho)
    ux, ur = sp.diff(u, x), sp.diff(u, rho)
    Dc = D_static(u, ux, ur, sp.diff(u, x, 2), sp.diff(u, x, rho), sp.diff(u, rho, 2), rho)
    F0, F1, F2 = F0_, F1_, F2_
    Dn = sp.simplify(Dc.doit().subs(x, eta * rho).doit())
    Dn = Dn.replace(lambda e: isinstance(e, sp.Subs), lambda e: e.doit())
    d2, d1 = sp.Derivative(F(eta), (eta, 2)), sp.Derivative(F(eta), eta)
    Dn = sp.simplify(Dn.subs(d2, F2).subs(d1, F1).subs(F(eta), F0))
    N = conical_N(3 if MUT == "conical-term" else 2)
    a_ok = zero(8 * rho**6 * Dn - N)
    check("conical", a_ok, f"b = F(x/rho): 8 rho^6 D = N(F, F', F''; eta), cubic in F'', N typed in the script: {a_ok}. "
          f"(C9 (b), (c) as first worded failed in dev_mid.run1.log, see the README, A1; C9b below replaces them)")


F0_, F1_, F2_ = sp.symbols("F0 F1 F2")
ETA = sp.Symbol("eta", positive=True)


def conical_N(c33=2):
    F0, F1, F2, eta = F0_, F1_, F2_, ETA
    return (2 * F0**3 * F1**2 * F2 + 4 * F0**3 * F1 * F2**2 * eta + c33 * F0**3 * F2**3 * eta**2 + F0**2 * F1**4 * eta**2
            + 2 * F0**2 * F1**4 + 2 * F0**2 * F1**3 * F2 * eta**3 + 12 * F0**2 * F1**3 * F2 * eta
            + F0**2 * F1**2 * F2**2 * eta**4 + 6 * F0**2 * F1**2 * F2**2 * eta**2 + 4 * F0 * F1**5 * eta**3
            + 8 * F0 * F1**5 * eta + 4 * F0 * F1**4 * F2 * eta**4 + 8 * F0 * F1**4 * F2 * eta**2 - 2 * F0 * F1**3 * eta
            - 6 * F0 * F1**2 * F2 * eta**2 - 2 * F0 * F1**2 * F2 - 6 * F0 * F1 * F2**2 * eta**3
            - 4 * F0 * F1 * F2**2 * eta - 2 * F0 * F2**3 * eta**4 - 2 * F0 * F2**3 * eta**2 + F1**6 * eta**6
            + 4 * F1**6 * eta**4 + 4 * F1**6 * eta**2 - F1**4 * eta**4 - 3 * F1**4 * eta**2 - 2 * F1**4
            - 2 * F1**3 * F2 * eta**5 - 6 * F1**3 * F2 * eta**3 - 4 * F1**3 * F2 * eta - F1**2 * F2**2 * eta**6
            - 3 * F1**2 * F2**2 * eta**4 - 2 * F1**2 * F2**2 * eta**2)


def c9b_parts():
    """C9b (i): the parts of degree 3 and 4 of N from the hand-derived entries, compared with the typed N."""
    F0, F1, F2, eta = F0_, F1_, F2_, ETA
    w = F1 + eta * F2
    hand = {"P": -(F1**2 + F0 * F2), "Q": -eta * w / 2, "Z": -w / 2, "Y": -(eta * F1**2 + F0 * w / 2),
            "h": eta**2 * F1**2 / 2}
    # the same entries from sympy's jets of F(x/rho), times rho^2
    x, rho = sp.symbols("x rho", positive=True)
    F = sp.Function("F")
    u = F(x / rho)
    P, Q, Z, Y, h = blocks_jets(u, sp.diff(u, x), sp.diff(u, rho), 0, 0, sp.diff(u, x, 2), sp.diff(u, x, rho), 0,
                                sp.diff(u, rho, 2), rho)
    d2, d1 = sp.Derivative(F(eta), (eta, 2)), sp.Derivative(F(eta), eta)

    def conv(e):
        e = (rho**2 * e).doit().subs(x, eta * rho).doit()
        e = e.replace(lambda q: isinstance(q, sp.Subs), lambda q: q.doit())
        return sp.simplify(e.subs(d2, F2).subs(d1, F1).subs(F(eta), F0))
    jets = {"P": conv(P), "Q": conv(Q), "Z": conv(Z), "Y": conv(Y), "h": conv(h)}
    ent_ok = all(zero(jets[k] - hand[k]) for k in hand)
    Dhat = D_of(hand["P"], hand["Q"], hand["Z"], hand["Y"], hand["h"])
    N_hand = sp.expand(8 * Dhat)
    N_typed = conical_N()
    n_ok = sp.expand(N_hand - N_typed) == 0
    tt = sp.Symbol("tt")
    exp3 = 2 * F0 * w**2 * ((F0**2 * F2 if MUT == "conical-hand" else (F0**2 - 1) * F2) - eta * w)
    exp4 = -w**2 * ((eta**2 + 1) * (eta**2 + 2) * F1**2 + 2 * F0 * (eta * w + F2))
    low3 = sp.expand(N_hand.subs({F1: tt * F1, F2: tt * F2}))
    low4 = sp.expand(N_hand.subs({F0: tt * F0, F1: tt * F1, F2: tt * F2}))
    typ3 = sp.expand(N_typed.subs({F1: tt * F1, F2: tt * F2}))
    typ4 = sp.expand(N_typed.subs({F0: tt * F0, F1: tt * F1, F2: tt * F2}))
    lowest3 = min(m[0] for m in sp.Poly(low3, tt).monoms())
    lowest4 = min(m[0] for m in sp.Poly(low4, tt).monoms())
    p3 = lowest3 == 3 and sp.expand(low3.coeff(tt, 3) - exp3) == 0 and sp.expand(typ3.coeff(tt, 3) - exp3) == 0
    p4 = lowest4 == 4 and sp.expand(low4.coeff(tt, 4) - exp4) == 0 and sp.expand(typ4.coeff(tt, 4) - exp4) == 0
    check("conical_parts", ent_ok and n_ok and p3 and p4,
          f"hand entries P = -(F'^2 + F F'')/rho^2, Q = -eta w/(2 rho^2), Z = -w/(2 rho^2), Y = -(eta F'^2 + F w/2)/rho^2, "
          f"h = eta^2 F'^2/(2 rho^2) equal sympy's jets of F(x/rho): {ent_ok}; 8 rho^6 D from them equals the typed N: "
          f"{n_ok}; degree 3 in (F', F'') is {sp.factor(exp3)}: {p3}; degree 4 with F, F', F'' small is "
          f"{sp.factor(exp4)}: {p4}")


def limit_scan(N, finfs, ks, kappas):
    """Can N = 0 have F -> F_inf by a power law c eta^-k or an exponential c e^(-kappa eta)? Returns the list of
    (kind, F_inf, parameter, reason) where it cannot be excluded at leading order."""
    F0, F1, F2, eta = F0_, F1_, F2_, ETA
    c = sp.Symbol("c", real=True)
    s = sp.Symbol("s", positive=True)  # eta = 1/s^2, so that half-integer k are integer powers
    uu = sp.Symbol("uu")
    open_ = []
    for fi in finfs:
        for k in ks:
            sub = {F0: fi + c * s**(2 * k), F1: -k * c * s**(2 * k + 2), F2: k * (k + 1) * c * s**(2 * k + 4),
                   eta: 1 / s**2}
            e = sp.expand(N.subs(sub))
            if e == 0:
                open_.append(("power", fi, k, "N vanishes identically"))
                continue
            e = sp.expand(e * s**200)
            low = sp.Poly(e, s)
            lo = min(m[0] for m in low.monoms())
            coef = sp.Poly(low.coeff_monomial(s**lo), c)
            q = sp.Poly(sp.cancel(coef.as_expr() / c**min(m[0] for m in coef.monoms())), c)
            roots = [r for r in sp.real_roots(q)] if q.degree() > 0 else []
            if roots:
                open_.append(("power", fi, k, f"real roots c = {[sp.N(r, 6) for r in roots]}"))
        for kap in kappas:
            sub = {F0: fi + uu, F1: -kap * uu, F2: kap**2 * uu}
            e = sp.Poly(sp.expand(N.subs(sub)), uu)
            if e.is_zero:
                open_.append(("exp", fi, kap, "N vanishes identically"))
                continue
            lo = min(m[0] for m in e.monoms())
            if sp.expand(e.coeff_monomial(uu**lo)) == 0:
                open_.append(("exp", fi, kap, "lowest coefficient is zero"))
    return open_


def formal_series_obstruction(order=5):
    """Observation, not a criterion: F = F_inf + a2/eta^2 + a3/eta^3 + ... in N = 0, the first conditions on the a_j."""
    F0, F1, F2, eta = F0_, F1_, F2_, ETA
    Fi, e = sp.Symbol("Finf"), sp.Symbol("e", positive=True)
    a = sp.symbols(f"a2:{order + 1}")
    Fs = Fi + sum(a[j - 2] * e**j for j in range(2, order + 1))
    d = lambda f: -e**2 * sp.diff(f, e)
    ex = sp.expand(conical_N().subs({F0: Fs, F1: d(Fs), F2: d(d(Fs)), eta: 1 / e}) * e**12)
    poly = sp.Poly(ex, e)
    sols, notes = {}, []
    for (m,) in sorted(poly.monoms())[:6]:
        cf = sp.factor(poly.coeff_monomial(e**m).subs(sols))
        if cf == 0:
            continue
        free = [q for q in a if cf.has(q) and q not in sols]
        if not free:
            notes.append(f"order e^{m - 12}: no unknown left, condition {cf} = 0")
            break
        var = min(free, key=lambda q: int(str(q)[1:]))
        r = sp.solve(cf, var)
        notes.append(f"order e^{m - 12}: {cf} = 0 gives {var} = {r}")
        if not r:
            break
        sols[var] = r[-1]
    return notes


def c9b_limit():
    eta = ETA
    F0, F1, F2 = F0_, F1_, F2_
    finfs = [0, sp.Rational(1, 3), 1, -1, 2]
    ks = [sp.Rational(1, 2), 1, sp.Rational(3, 2), 2, sp.Rational(5, 2), 3, 4, 5, 6]
    kappas = [sp.Rational(1, 2), 1, 2]
    open_main = limit_scan(conical_N(), finfs, ks, kappas)
    check("conical_limit", not open_main,
          f"N = 0 has no limit at eta -> +inf by c eta^-k (k in {[str(k) for k in ks]}) or c e^(-kappa eta) (kappa in "
          f"{[str(k) for k in kappas]}) for F_inf in {[str(f) for f in finfs]}: cases not excluded at leading order: "
          f"{open_main or 'none'}")
    for n in formal_series_obstruction():
        print(f"  observation (not a criterion), formal series F_inf + sum a_j/eta^j, a2 = -2 F_inf branch: {n}")


def c9b_planted():
    eta = ETA
    F0, F1, F2 = F0_, F1_, F2_
    finfs = [0, sp.Rational(1, 3), 1, -1, 2]
    Npow = F2 + F1 / eta if MUT == "conical-plant" else F2 + 2 * F1 / eta
    Nexp = F2 + F1
    pl = limit_scan(Npow, finfs, [1, 2], [])
    pe = limit_scan(Nexp, finfs, [], [1, 2])
    pow_hit = sorted({(r[1], r[2]) for r in pl})
    exp_hit = sorted({(r[1], r[2]) for r in pe})
    want = sorted({(fi, 1) for fi in finfs})
    check("conical_planted", pow_hit == want and exp_hit == want,
          f"the same scan on the planted N = F'' + 2F'/eta (solutions F_inf + c/eta) detects a limit at k = 1 for every "
          f"F_inf and not at k = 2: {pow_hit == want}; on N = F'' + F' (solutions F_inf + c e^(-eta)) at kappa = 1 and not at "
          f"kappa = 2: {exp_hit == want}")


# ---------------------------------------------------------------------------------------------------------------
# C10 farfield
# ---------------------------------------------------------------------------------------------------------------
def c10_farfield():
    ui, e, rho = sp.symbols("ui e rho")
    w0, wx, wr, wxx, wxr, wrr = sp.symbols("w0 wx wr wxx wxr wrr")
    Dff = sp.expand(D_static(ui + e * w0, e * wx, e * wr, e * wxx, e * wxr, e * wrr, rho))
    co = [sp.expand(Dff.coeff(e, k)) for k in range(4)]
    Z = wxr / 2
    sgn = 1 if MUT == "far-sign" else -1
    pred = sgn * ui * Z**2 * ((1 - ui**2) * wxx + wrr + wr / rho)
    ok = co[0] == 0 and co[1] == 0 and co[2] == 0 and sp.expand(co[3] - pred) == 0
    check("farfield", ok, f"u = u_inf + e w: the e^0, e^1, e^2 coefficients of D vanish "
          f"({[c == 0 for c in co[:3]]}) and the e^3 coefficient is {'-' if sgn < 0 else '+'}u_inf Z^2 "
          f"[(1 - u_inf^2) w_xx + Delta_perp w], Z = w_x,rho/2: {sp.expand(co[3] - pred) == 0}")


# ---------------------------------------------------------------------------------------------------------------
# C11 bubble
# ---------------------------------------------------------------------------------------------------------------
def c11_bubble():
    r, c, s, z = sp.symbols("r c s z")
    w0, w1, w2 = sp.symbols("w0 w1 w2")
    ux, ur = w1 * c, w1 * s
    uxx = w2 * c**2 + w1 * s**2 / r
    uxr = (w2 - w1 / r) * c * s
    urr = w2 * s**2 + w1 * c**2 / r
    rho = r * s
    P, Q, Zq, Y, h = blocks_jets(w0, ux, ur, 0, 0, uxx, uxr, 0, urr, rho)
    if MUT == "bubble-w2":
        P = -(ux**2)
        P = sp.expand(P)
    D = D_of(P, Q, Zq, Y, h)
    num, den = sp.fraction(sp.together(sp.expand(16 * r**3 * D)))
    num = sp.expand(sp.rem(sp.Poly(sp.expand(num), s).as_expr(), s**2 + c**2 - 1, s))
    P0 = sp.Poly(sp.expand(num / den), c)
    coef = [sp.factor(P0.coeff_monomial(c**k)) for k in range(P0.degree() + 1)]
    odd = any(P0.coeff_monomial(c**k) != 0 for k in (1, 3, 5)) or sp.expand(num).has(s)
    c0 = 2 * w1**2 * r * (w1**4 * r**2 - 4 * w0**2 * w1**2 - (w1 + r * w2)**2)
    c6 = -2 * (2 * r * w0 * w2 - 2 * w0 * w1 + r * w1**2) * (r * w0 * w2 - w0 * w1 + r * w1**2)**2
    a_ok = P0.degree() == 6 and not odd and sp.expand(coef[0] - c0) == 0 and sp.expand(coef[6] - c6) == 0
    A, B = 2 * w0 - r * w1, w1 + r * w2
    zz = sp.Symbol("zz")
    gens, sat = (zz, w0, w2, w1, r), 1 - zz * w1 * r
    want = [A, B**2]
    nz = [cf for cf in coef if cf != 0]
    Gw = sp.groebner(want + [sat], *gens, order="grevlex")
    Gc = sp.groebner(nz + [sat], *gens, order="grevlex")
    sub = [Gw.contains(p) for p in nz]
    sup = [Gc.contains(p) for p in want]
    b_ok = all(sub) and all(sup)
    # ODE step: differentiate 2 w = r w' (-> w' - r w'' = 0), and w' + r w'' = 0 (or r w'' = +w' for the control)
    rel2 = (w1 - r * w2) if MUT == "bubble-ode" else (w1 + r * w2)
    sol = sp.solve([w1 - r * w2, rel2], [w1, w2], dict=True)
    c_ok = sol == [{w1: 0, w2: 0}]
    check("bubble", a_ok and b_ok and c_ok, f"u = w(r): 16 r^3 D is a cubic in c^2 with the printed c^0 and c^6 "
          f"coefficients: {a_ok}; its coefficients and (2w - r w', (w' + r w'')^2) generate the same ideal after inverting "
          f"r w' (inclusions {all(sub)}, {all(sup)}); differentiating 2 w = r w' and adding w' + r w'' = 0 gives "
          f"w' = 0: {c_ok}. So D = 0 on no open set of the constant-speed wall")


CHECKS = [("d_block", c1_d_block), ("travel", c2_travel), ("planar", c3_planar), ("planar_closure", c4_planar_closure),
          ("planar_traveling", c5_planar_traveling), ("bessel", c6_bessel), ("additive", c7_additive),
          ("multiplicative", c8_multiplicative), ("conical", c9_conical), ("conical_parts", c9b_parts),
          ("conical_limit", c9b_limit), ("conical_planted", c9b_planted), ("farfield", c10_farfield),
          ("bubble", c11_bubble)]


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

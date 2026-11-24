#!/usr/bin/env python3
"""Closure for planar shifts that cross |u| = 1 or do not travel (issue #161).

Planar shift b(t, x): R = K g2 + 0 with K = -a_x, a = b_t - b b_x (det_classify.md section 2). Closure of the normal +
phantom scalar test holds iff g2 = -2 exp(2 H(p) J(q)) dp dq, i.e. C[b] = box ln|K| - 2 K = 0 on {K != 0}, box the 2D
d'Alembertian. This script
  C1  the criterion in null coordinates and its planted passes and failures;
  C2  the criterion in a Painleve-Gullstrand chart, C[u] = Theta' for a static shift, planted failures (Theorem T's
      case, #80's planted failure) and a formal planted pass;
  C3  the crossing analysis for a static (traveling) shift: horizons, the map s, Phi = K (1 - u^2)  (Theorem T');
  C4  Gaussian normal chart = PG chart, Jacobi equation A_tt = K A, b_x = -A_t/A;
  C5  exterior of H' J' >= 0 lumps, geodesic derivation of A_tt = K A, total curvature 2 h j;
  C6  necessary conditions for equal limits b_+ = b_-, linear order;
  E1-E5  numerical evidence for the localized PG chart of sigma = H(p) J(q) with H', J' >= 0 (Theorem E).
Imports: phantom_source (ricci), two_scalars (closure, dform). Not det_classify (it reads --mutate at import).
Usage: planar_closure.py [--no-controls] [--mutate=NAME] [--only=check1,check2] [--logdir=DIR]. Exit 1 if a check fails
or a control does not fire on its target."""
import os
import subprocess
import sys
import time

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phantom_source as ps  # noqa: E402
import two_scalars as ts  # noqa: E402

# name: (target check, what the mutation does)
MUTATIONS = {
    "crit-factor": ("crit_null", "expects d_p d_q ln|R_pq| = -e^(2 sigma) C"),
    "crit-deform": ("crit_null", "adds p q^2/10 to sigma of the planted passes"),
    "crit-static-g": ("crit_pg", "g^xx = 1 + u^2 in box"),
    "crit-plant-jet": ("crit_pg", "the planted pass jets get c6 + 1/7"),
    "ivb-flat": ("crit_pg_ivb", "u1 replaced by the constant 1/3 (K = 0 at every point of the list)"),
    "hor-orient": ("horizon", "s' = 2/(1 + u^2)"),
    "hor-alpha": ("horizon", "alpha = +u u'(xi0)"),
    "jacobi-sign": ("jacobi", "K = -A_tt/A"),
    "jacobi-theta": ("jacobi", "b_x = +A_t/A"),
    "ext-hol": ("exterior", "expects integral of K dA = h j"),
    "ext-null": ("exterior", "Gamma^p_pp = sigma_p"),
    "lin-order": ("limits_necessary", "K1 = -b1_tt"),
    "conv-zero": ("limits_necessary", "G = 0 in the convolution"),
    "E-norm": ("E_norm", "initial velocity (1, 1)"),
    "E-kneg": ("E_focus", "A'' = -K A"),
    "E-ext": ("E_ext", "support of J' is [0, 2]"),
    "E-hol": ("E_hol", "expects L'(10) = h j"),
    "E-geo": ("E_geo", "finite-difference spacing 0.5"),
}

RESULTS = []
OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate")
ONLY = OPTS["only"].split(",") if "only" in OPTS else None


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def zero(e):
    if isinstance(e, sp.MatrixBase):
        return all(sp.simplify(x) == 0 for x in e)
    return sp.simplify(e) == 0


def christoffel(g, X):
    gi = g.inv()
    n = len(X)
    return [[[sum(gi[a, d] * (sp.diff(g[d, c], X[b]) + sp.diff(g[d, b], X[c]) - sp.diff(g[b, c], X[d]))
                  for d in range(n)) / 2 for c in range(n)] for b in range(n)] for a in range(n)]


def box(f, g, X, sq):
    """d'Alembertian (1/sq) d_a (sq g^{ab} d_b f), sq = sqrt|g| given explicitly."""
    gi = g.inv()
    return sum(sp.diff(sq * gi[a, b] * sp.diff(f, X[b]), X[a]) for a in range(len(X)) for b in range(len(X))) / sq


# ---------------------------------------------------------------------------------------------------------------
# PG chart helpers: g = -dt^2 + (dx + b dt)^2, g^{-1} = -n n + d_x d_x, n = (1, -b), sqrt|g| = 1
# ---------------------------------------------------------------------------------------------------------------
def pg_metric(b):
    return sp.Matrix([[-1 + b**2, b], [b, 1]])


def K_of(b, t, x):
    """K = -a_x, a = b_t - b b_x."""
    a = sp.diff(b, t) - b * sp.diff(b, x)
    return -sp.diff(a, x)


def box_pg(f, b, t, x, gxx=None):
    """box f in a PG chart. Mutation crit-static-g uses g^{xx} = 1 + b^2."""
    gxx = (1 + b**2) if MUT == "crit-static-g" else (1 - b**2)
    gi = sp.Matrix([[-1, b], [b, gxx]])
    X = [t, x]
    return sum(sp.diff(gi[i, j] * sp.diff(f, X[j]), X[i]) for i in range(2) for j in range(2))


def C_static(u, xi):
    """C[u] = box ln|K| - 2 K for the static shift b = u(x); ln|K| as (1/2) ln K^2."""
    t = sp.Symbol("t", real=True)
    K = sp.diff(u * sp.diff(u, xi), xi)
    return box_pg(sp.log(K**2) / 2, u, t, xi) - 2 * K


def theta_prime(u, xi):
    K = sp.diff(u * sp.diff(u, xi), xi)
    Th = (1 - u**2) * sp.diff(K, xi) / K - 2 * u * sp.diff(u, xi)
    return sp.diff(Th, xi)


# ---------------------------------------------------------------------------------------------------------------
# C1 crit_null
# ---------------------------------------------------------------------------------------------------------------
def null_planted(sigma):
    """exact R = a a - b b, `closure` ok, and C = box ln|K| - 2K (null chart, covariant formula) for sigma(p, q)."""
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
    g2 = g[:2, :2]
    K = 2 * sp.diff(s, p, q) * sp.exp(-2 * s)
    C = box(sp.log(K**2) / 2, g2, [p, q], sp.exp(2 * s)) - 2 * K
    return exact, ok, zero(C)


def c1_crit_null():
    p, q = sp.symbols("p q", real=True)
    sig = sp.Function("sigma")(p, q)
    g = sp.Matrix([[0, -sp.exp(2 * sig)], [-sp.exp(2 * sig), 0]])
    X = [p, q]
    R = ps.ricci(g, X, lambda E, k: sp.diff(E, X[k]), sp.simplify)
    K = 2 * sp.diff(sig, p, q) * sp.exp(-2 * sig)
    i_a = zero(R - K * g) and zero(R[0, 1] + 2 * sp.diff(sig, p, q))
    f = sp.Function("f")(p, q)
    i_b = zero(box(f, g, X, sp.exp(2 * sig)) + 2 * sp.exp(-2 * sig) * sp.diff(f, p, q))
    Cexpr = box(sp.log(K**2) / 2, g, X, sp.exp(2 * sig)) - 2 * K
    fac = -sp.exp(2 * sig) if MUT == "crit-factor" else -sp.exp(2 * sig) / 2
    i_c = zero(sp.diff(sp.log(R[0, 1]**2) / 2, p, q) - fac * Cexpr)
    deform = (lambda P, Q: P * Q**2 / 10) if MUT == "crit-deform" else (lambda P, Q: 0)
    passes = [("sin(p) e^(q/2)/3", lambda P, Q: sp.sin(P) * sp.exp(Q / 2) / 3 + deform(P, Q)),
              ("p q + p^2 q", lambda P, Q: P * Q + P**2 * Q + deform(P, Q)),
              ("(p + q)^2/5", lambda P, Q: (P + Q)**2 / 5 + deform(P, Q))]
    fails = [("p^2 q/2 + p q^2/2", lambda P, Q: P**2 * Q / 2 + P * Q**2 / 2),
             ("e^(p q)", lambda P, Q: sp.exp(P * Q)),
             ("p^3 q^2/6 + p q", lambda P, Q: P**3 * Q**2 / 6 + P * Q)]
    okall, rows = True, []
    for lab, fn in passes + fails:
        want = (lab, fn) in passes
        exact, ok, c0 = null_planted(fn)
        agree = exact and ok == want and c0 == want
        okall &= agree
        rows.append(f"sigma = {lab}: R = a a - b b {exact}, closure {ok}, C = 0 {c0}, expected {want}")
    for r in rows:
        print("  " + r)
    check("crit_null", i_a and i_b and i_c and okall,
          f"(i) R = K g, K = 2 sigma_pq e^(-2 sigma), R_pq = -2 sigma_pq: {i_a}; box f = -2 e^(-2 sigma) f_pq: {i_b}; "
          f"d_p d_q ln|R_pq| = {'-e^(2 sigma)' if MUT == 'crit-factor' else '-e^(2 sigma)/2'} C: {i_c}; "
          f"(ii) 3 planted passes and 3 planted failures agree with `closure` and with C = 0: {okall}")


# ---------------------------------------------------------------------------------------------------------------
# C2 crit_pg
# ---------------------------------------------------------------------------------------------------------------
def taylor_C_coeffs(cs, xi, nco):
    """Taylor coefficients 0..nco-1 of C[u] for u = sum cs[k] xi^k / k!."""
    u = sum(sp.Rational(1) * c * xi**k / sp.factorial(k) for k, c in enumerate(cs))
    return sp.series(C_static(u, xi), xi, 0, nco).removeO()


def planted_jet(perturb):
    """c0..c3 given, c4..c8 solved so that C[u] = O(xi^5). perturb adds to c6."""
    xi = sp.Symbol("xi", real=True)
    cs = [sp.Rational(1, 3), sp.Rational(1, 2), sp.Rational(1, 5), sp.Rational(2, 7)]
    for m in range(5):
        c = sp.Symbol("c")
        ser = taylor_C_coeffs(cs + [c], xi, m + 1)
        sol = sp.solve(sp.expand(ser).coeff(xi, m), c)
        assert len(sol) == 1, (m, sol)
        cs.append(sol[0])
    if perturb:
        cs[6] += perturb
    return xi, cs


def c2_crit_pg():
    t, x, y, z = sp.symbols("t x y z", real=True)
    b = sp.Function("b")(t, x)
    X2 = [t, x]
    R = ps.ricci(pg_metric(b), X2, lambda E, k: sp.diff(E, X2[k]), sp.simplify)
    Kb = -(sp.diff(b, t, x) - b * sp.diff(b, x, 2) - sp.diff(b, x)**2)
    i_a = zero(R - Kb * pg_metric(b)) and zero(Kb - K_of(b, t, x))
    xi = sp.Symbol("xi", real=True)
    u = sp.Function("u")(xi)
    i_b = zero(C_static(u, xi) - theta_prime(u, xi))
    # (iii) the obstruction of `closure` on the 4D PG metric of a static undefined u
    Ku = sp.diff(u * sp.diff(u, xi), xi)
    th0, th1 = sp.Matrix([1, 0, 0, 0]), sp.Matrix([u, 1, 0, 0])
    kk, mm = th0 - th1, -Ku * (th0 + th1)
    aa, bb = (kk + mm) / 2, (kk - mm) / 2
    ok2, om, why = ts.closure(aa, bb, [t, xi, y, z], sp.simplify)
    obstruction = ts.dform(om, [t, xi, y, z])[0][1]
    i_c = zero(obstruction - sp.Rational(1, 2) * theta_prime(u, xi))
    # (iv) Theorem T's case, as worded in the README: 5 fixed points, K != 0 and C != 0 (FAILS for u3 at xi = 1/3, K = 0)
    pts = [sp.Rational(1, 3), 1, sp.Rational(3, 2), 2, sp.Rational(5, 2)]
    profiles = [("u1 = (1 - 1/(1 + xi^2)^2)/2", (1 - 1 / (1 + xi**2)**2) / 2),
                ("u2 = 1/(2(1 + xi^2))", 1 / (2 * (1 + xi**2))),
                ("u3 = 1/(3(1 + xi^2)^2)", 1 / (3 * (1 + xi**2)**2))]
    i_d = True
    for lab, uu in profiles:
        Kx = sp.diff(uu * sp.diff(uu, xi), xi)
        Cx = C_static(uu, xi)
        vals = [(sp.simplify(Kx.subs(xi, pt)), sp.simplify(Cx.subs(xi, pt))) for pt in pts]
        bad = [str(pt) for pt, (kv, cv) in zip(pts, vals) if kv == 0 or cv == 0]
        good = not bad
        i_d &= good
        print(f"  (iv) {lab}: K != 0 and C != 0 at the {len(pts)} points: {good}"
              f"{'' if good else '; K = 0 or C = 0 at xi = ' + ', '.join(bad)}", flush=True)
    # (iv-b) amendment A3: fixed list, only points with K != 0 verified exactly, first 5 taken
    ell = [sp.Rational(1, 5), sp.Rational(1, 3), sp.Rational(1, 2), sp.Rational(2, 3), 1, sp.Rational(3, 2), 2,
           sp.Rational(5, 2), 3, 4]
    i_f = True
    for lab, uu in profiles:
        if MUT == "ivb-flat" and lab.startswith("u1"):
            uu = sp.Rational(1, 3) + 0 * xi
        Kx = sp.diff(uu * sp.diff(uu, xi), xi)
        taken = [pt for pt in ell if sp.simplify(Kx.subs(xi, pt)) != 0][:5]
        Cx = C_static(uu, xi) if taken else None
        good = len(taken) == 5 and all(sp.simplify(Cx.subs(xi, pt)) != 0 for pt in taken)
        i_f &= good
        print(f"  (iv-b) {lab}: points taken {[str(pt) for pt in taken]}, C != 0 at all: {good}", flush=True)
    # (v) formal planted pass and its perturbation
    perturb = sp.Rational(1, 7) if MUT == "crit-plant-jet" else 0
    xi2, cs = planted_jet(perturb)
    ser = sp.expand(taylor_C_coeffs(cs, xi2, 5))
    pass_ok = all(ser.coeff(xi2, m) == 0 for m in range(5))
    xi3, cs2 = planted_jet(0)
    cs2[6] += sp.Rational(1, 7)
    ser2 = sp.expand(taylor_C_coeffs(cs2, xi3, 5))
    fail_ok = ser2.coeff(xi3, 2) != 0
    print(f"  planted jets c0..c8 = {[str(c) for c in cs]}")
    i_e = pass_ok and fail_ok
    check("crit_pg", i_a and i_b and i_c and i_e,
          f"(i) R = K g, K = -(b_tx - b b_xx - b_x^2): {i_a}; (ii) C[u] = Theta' for an undefined static u: {i_b}; (iii) the "
          f"obstruction of `closure` is Theta'/2: {i_c}; (v) a formal solution of C = 0 with C = O(xi^5): {pass_ok}, the same "
          f"with c6 + 1/7 has a nonzero xi^2 coefficient: {fail_ok}")
    check("crit_pg_iv", i_d,
          f"(iv) 3 subluminal traveling bubbles have K != 0 and C != 0 at 5 fixed exact points each: {i_d} (as worded in the "
          f"README; coordinator ruling A keeps it as worded)")
    check("crit_pg_ivb", i_f,
          f"(iv-b) the same 3 profiles at the first 5 points of a fixed list where K != 0 is verified exactly: C != 0 at all "
          f"15 points: {i_f}")


# ---------------------------------------------------------------------------------------------------------------
# C3 horizon
# ---------------------------------------------------------------------------------------------------------------
def sprime(u):
    return 2 / (1 + u**2) if MUT == "hor-orient" else 2 / (1 - u**2)


def pole_residue(expr, xi, e):
    """residue of a rational expr at a simple pole e; 0 if e is not a pole."""
    num, den = sp.fraction(sp.together(expr))
    if sp.simplify(den.subs(xi, e)) != 0:
        return 0
    d1 = sp.diff(den, xi).subs(xi, e)
    assert sp.simplify(d1) != 0, "double pole"
    return sp.simplify(num.subs(xi, e) / d1)


def end_behavior(uu, xi, e):
    """sign (+1/-1) of the divergence of s towards the end e of a component; 0 = no divergence."""
    sp_ = sprime(uu)
    if e in (sp.oo, -sp.oo):
        lim = sp.limit(sp_, xi, e)
        # s' -> L != 0 at e: s ~ L xi, so s -> sign(L) * sign(e) * oo
        return 0 if lim == 0 else int(sp.sign(lim) * (1 if e == sp.oo else -1))
    r = pole_residue(sp_, xi, e)
    # s ~ r ln|xi - e| -> -sign(r) * oo
    return 0 if r == 0 else int(-sp.sign(r))


def c3_horizon():
    xi, eta = sp.symbols("xi eta", real=True)
    u1, u2, u3, al, gam = sp.symbols("u1 u2 u3 alpha gamma", real=True)
    # (i) the ODE of closure, Theta = 2 alpha: (1 - u^2) K' = (2 alpha + 2 u u') K; with K = kappa eta^gamma the ODE divided
    # by kappa eta^gamma reads (1 - u^2) gamma / eta = 2 alpha + 2 u u'
    uu = 1 + u1 * eta + u2 * eta**2 + u3 * eta**3
    E = sp.simplify((1 - uu**2) * gam / eta - (2 * al + 2 * uu * sp.diff(uu, eta)))
    lead = sp.limit(E, eta, 0)
    gsol = sp.solve(lead, gam)
    i_a = len(gsol) == 1 and sp.simplify(gsol[0] + (al + u1) / u1) == 0
    print("  (i) leading order: gamma =", gsol, flush=True)
    # at a horizon point with K(xi0) != 0 the ODE gives 2 alpha + 2 u u' = 0
    uv, up = sp.symbols("uv up", real=True)
    alpha_sol = sp.solve(sp.Eq(2 * al + 2 * uv * up, 0), al)[0]
    claim = uv * up if MUT == "hor-alpha" else -uv * up
    i_b = sp.simplify(alpha_sol - claim) == 0
    # (ii) the superluminal bubble u_s = 3/(1 + xi^2): two horizons, opposite alpha
    us = 3 / (1 + xi**2)
    Ks = sp.diff(us * sp.diff(us, xi), xi)
    h1, h2 = sp.sqrt(2), -sp.sqrt(2)
    al1 = sp.simplify(-(us * sp.diff(us, xi)).subs(xi, h1))
    al2 = sp.simplify(-(us * sp.diff(us, xi)).subs(xi, h2))
    horizon_K = sp.simplify(Ks.subs(xi, h1)) != 0 and sp.simplify(Ks.subs(xi, h2)) != 0
    i_c = horizon_K and sp.sign(al1) == -sp.sign(al2) and al1 != 0
    print(f"  (ii) u_s = 3/(1 + xi^2): |u_s| = 1 at xi = +-sqrt2, K(+-sqrt2) != 0: {horizon_K}; alpha = {al1} and {al2}",
          flush=True)
    # (iii) components, orientation of s, Phi -> 0 at both ends
    u2p = 1 / (2 * (1 + xi**2))
    r2 = sp.sqrt(2)
    comps = [("u_s", us, -sp.oo, -r2, -3), ("u_s", us, -r2, r2, sp.Rational(1, 7)), ("u_s", us, r2, sp.oo, 3),
             ("u2", u2p, -sp.oo, sp.oo, sp.Rational(1, 7))]
    i_d = True
    for lab, uu_, a, b_, mid in comps:
        sgn_mid = sp.sign(sp.simplify(sprime(uu_).subs(xi, mid)))
        e_a, e_b = end_behavior(uu_, xi, a), end_behavior(uu_, xi, b_)
        Kc = sp.diff(uu_ * sp.diff(uu_, xi), xi)
        Phi = sp.together(Kc * (1 - uu_**2))
        zeros = []
        for e in (a, b_):
            zeros.append(sp.limit(Phi, xi, e) == 0 if e in (sp.oo, -sp.oo) else sp.simplify(Phi.subs(xi, e)) == 0)
        Kmid = sp.simplify(Kc.subs(xi, mid)) != 0
        good = sgn_mid != 0 and e_a * e_b == -1 and all(zeros) and Kmid
        i_d &= good
        print(f"  (iii) component ({a}, {b_}) of {lab}: sign s' = {sgn_mid}, s diverges with sign {e_a} at the left end and "
              f"{e_b} at the right, Phi -> 0 at both: {all(zeros)}, K != 0 inside: {Kmid} -> excludes Phi = C e^(alpha s), "
              f"C != 0: {good}", flush=True)
    check("horizon", i_a and i_b and i_c and i_d,
          f"(i) Theta' = 0 at a simple crossing: K ~ |eta|^gamma, gamma = -(alpha + u1)/u1: {i_a}; alpha = "
          f"{'+' if MUT == 'hor-alpha' else '-'}u u'(xi0) when K(xi0) != 0: {i_b}; (ii) u_s has opposite alpha at its two "
          f"horizons: {i_c}; (iii) on each of 4 components s runs from one infinity to the other and Phi -> 0 at both ends, so "
          f"Phi = C e^(alpha s) forces C = 0: {i_d}")


CHECKS = [(("crit_null",), c1_crit_null), (("crit_pg", "crit_pg_iv", "crit_pg_ivb"), c2_crit_pg),
          (("horizon",), c3_horizon)]


def wanted(name):
    return ONLY is None or name in ONLY


def main():
    t0 = time.time()
    if MUT is not None and MUT not in MUTATIONS:
        print(f"unknown mutation {MUT}")
        sys.exit(2)
    for names, fn in CHECKS:
        if any(wanted(n) for n in names):
            t1 = time.time()
            fn()
            print(f"  ({time.time() - t1:.0f} s)", flush=True)
    failed = [n for n, ok in RESULTS if not ok]
    print(f"\nchecks: {len(RESULTS) - len(failed)} passed, {len(failed)} failed {failed}. {time.time() - t0:.0f} s")
    if MUT is not None:
        want = MUTATIONS[MUT][0]
        print(f"control {MUT}: target check {want}, fired {want in failed}")
    return failed


if __name__ == "__main__":
    failed = main()
    controls_ok = True
    if MUT is None and "no-controls" not in OPTS and ONLY is None:
        print("\n== Controls: each mutation must make its target check fail ==", flush=True)
        for name, (want, what) in MUTATIONS.items():
            p = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True)
            if "logdir" in OPTS:
                with open(os.path.join(OPTS["logdir"], f"mut_{name}.log"), "w") as fh:
                    fh.write(p.stdout + p.stderr)
            fails = [ln.split(":")[0][7:] for ln in p.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = want in fails and p.returncode == 1
            extra = [f for f in fails if f != want and f not in failed]
            controls_ok &= fired
            print(f"[{'FIRED' if fired and not extra else 'SILENT' if not fired else 'FIRED+EXTRA'}] --mutate={name} "
                  f"({what}): exit {p.returncode}, failed {fails}", flush=True)
            if p.returncode not in (0, 1):
                print(p.stderr[-2000:])
        print(f"controls: {'all fired on target' if controls_ok else 'SOME DID NOT FIRE ON TARGET'}")
    sys.exit(0 if not failed and controls_ok else 1)

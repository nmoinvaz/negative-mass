#!/usr/bin/env python3
"""Can a normal plus a phantom scalar source the Alcubierre wall? (issue #65)

Two minimally coupled scalars, phi_1 normal and phi_2 phantom, with any potential V(phi_1, phi_2), have
T_ab = d_a phi_1 d_b phi_1 - d_a phi_2 d_b phi_2 - g_ab ((d phi_1)^2 / 2 - (d phi_2)^2 / 2 + V), and G_ab = 8 pi T_ab
gives R_ab = lam g_ab + a_a a_b - b_a b_b with lam = 8 pi V, a = sqrt(8 pi) d phi_1 and b = sqrt(8 pi) d phi_2. The test
has four stages.

    algebraic   some real lam makes R - lam g a real form with at most one positive and one negative square
    closure     some field-space boost (a, b) -> (cosh chi a + sinh chi b, sinh chi a + cosh chi b) makes a and b closed
    potential   d lam ^ a ^ b = 0, so that V is a function of (phi_1, phi_2)
    residual    G_ab - 8 pi T_ab = 0 with the fields substituted back

The script applies the test to Alcubierre's metric (the frame Ricci tensor of phantom_source.py, issues #57 and #61)
and to a control, the Ellis-Bronnikov wormhole with an added normal scalar.

    python3 notes/warp/two_scalars.py                 all checks, then every control (exit 1 if any check fails)
    python3 notes/warp/two_scalars.py --mutate=NAME   one control (its target checks must fail, exit 1)
    python3 notes/warp/two_scalars.py --no-controls   the checks only

Shared with phantom_source.py, imported read-only: `ricci`, the comoving spherical frame Ricci tensor
`sph_tetrad_ricci`, the lab frame Ricci tensor `lab_tetrad_ricci`, the section 6 route `shift_tetrad_ricci` with its
jet rules, the profile jets and the sample points. The two-scalar test, the determinant analysis and the control are
new here.
"""
import os
import random
import subprocess
import sys
import time

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phantom_source as ps  # noqa: E402

# name: (checks that must fail, what the mutation does)
MUTATIONS = {
    "planted-sign": (["planted_pass"], "builds the planted two-scalar forms as lam eta + a a + b b"),
    "planted-fail-sign": (["planted_fail"], "builds the planted failures with one normal and one phantom square"),
    "lemma-weight": (["lemma_samples"], "builds the lemma's samples with R_00 - R_11 + 2 R_22 = 0 in place of (I2)"),
    "alc-drop-w2": (["det_axis", "det_jets"], "drops the shift's square from g_tt (as in phantom_source.py)"),
    "ode-sign": (["no_open_set"], "flips the sign of the W2 relation before the ODE step"),
    "alc-planted": (["alc_points"], "replaces R at the off-axis sample points by a planted two-scalar form"),
    "eb2-ads": (["eb2_algebraic"], "takes the two-phantom metric, cos(2 kappa l) -> cosh(2 kappa l), beta > 1"),
    "eb2-deform": (["eb2_closed"], "adds (1 - beta) cos(4 kappa l) / (20 kappa^2) to the control's r(l)^2"),
    "eb2-unclosed": (["eb2_closed"], "multiplies the extracted b by 1 + cos(theta)^2"),
    "eb2-potential": (["eb2_potential"], "adds kappa^2 cos(theta) / 10 to lam in the potential stage"),
    "eb2-residual": (["eb2_residual"], "scales the phantom gradient by 11/10 before substituting it back"),
}

RESULTS = []


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# The two-scalar test
# ---------------------------------------------------------------------------------------------------------------
lam = sp.Symbol("lam")


def sign_of(e, sampler):
    """+1, -1 or 0 for an expression whose sign is fixed on the sampled domain, None if it is not."""
    if e == 0:
        return 0
    if e.is_positive:
        return 1
    if e.is_negative:
        return -1
    if sampler is None:
        return None
    vals = [sp.N(e.subs(sampler()), 30) for _ in range(8)]
    if all(x > 0 for x in vals):
        return 1
    if all(x < 0 for x in vals):
        return -1
    return None


def ssqrt(e):
    """A square root of e that is smooth where e > 0: even powers of the factors of e leave the root unsigned."""
    num, den = sp.fraction(sp.together(e))

    def root(p):
        c, fl = sp.factor_list(sp.expand(p))
        out, inside = sp.Integer(1), c
        for f, m in fl:
            out *= f**(m // 2)
            inside *= f**(m % 2)
        return out * sp.sqrt(inside)
    return root(num) / root(den)


def split_form(M, sampler, simp):
    """Write a symmetric M of rank <= 2 as a a^T - b b^T. Returns (a, b) or None if M has two squares of one sign.

    With X a vector where c = M(X, X) != 0 and q = M X, M - q q^T / c has rank one less than M. A zero a or b
    stands for a missing square."""
    n = M.shape[0]
    zero = sp.zeros(n, 1)
    if all(simp(x) == 0 for x in M):
        return zero, zero
    i = next((k for k in range(n) if simp(M[k, k]) != 0), None)
    X = sp.zeros(n, 1)
    if i is None:
        i, j = next((p, q) for p in range(n) for q in range(p + 1, n) if simp(M[p, q]) != 0)
        X[i] = X[j] = 1
    else:
        X[i] = 1
    q = M * X
    c = simp((X.T * M * X)[0])
    s1 = sign_of(c, sampler)
    if s1 is None:
        return None
    first = q / ssqrt(s1 * c)
    M1 = (M - s1 * first * first.T).applyfunc(simp)
    if all(x == 0 for x in M1):
        return (first, zero) if s1 > 0 else (zero, first)
    j = next(k for k in range(n) if M1[k, k] != 0)
    s2 = sign_of(M1[j, j], sampler)
    if s2 is None or s2 == s1:
        return None
    second = M1[:, j] / ssqrt(s2 * M1[j, j])
    if any(simp(x) != 0 for x in (M1 - s2 * second * second.T)):
        return None  # rank of M1 > 1
    return (first, second) if s1 > 0 else (second, first)


def candidates(Rm, gm):
    """Real lam for which R - lam g has rank <= 2: roots of the gcd over its 3 x 3 minors, each cubic in lam."""
    n = Rm.shape[0]
    A = Rm - lam * gm
    polys = []
    rows = [(i, j, k) for i in range(n) for j in range(i + 1, n) for k in range(j + 1, n)]
    for r in rows:
        for c in rows:
            num, _ = sp.fraction(sp.together(A.extract(list(r), list(c)).det()))
            num = sp.expand(num)
            if num != 0:
                polys.append(sp.Poly(num, lam))
    gcd = polys[0]
    for p in polys[1:]:
        gcd = sp.gcd(gcd, p)
    if gcd.degree() < 1:
        return []
    return [lr for lr in sp.roots(gcd.as_expr(), lam) if lr.is_real is not False]


def dform(k, X):
    """(dk)_{mu nu} = d_mu k_nu - d_nu k_mu."""
    n = len(X)
    return [[sp.diff(k[v], X[u]) - sp.diff(k[u], X[v]) for v in range(n)] for u in range(n)]


def closure(a, b, X, simp):
    """Is there chi with cosh chi a + sinh chi b and sinh chi a + cosh chi b both closed?

    With k = a + b and m = a - b this asks for e^chi k and e^-chi m closed, i.e. omega = d chi with
    omega ^ k = -dk and omega ^ m = dm. Contracting with K and Mv, where k(K) = m(Mv) = 1 and k(Mv) = m(K) = 0,
    gives the only candidate omega = dm(K, Mv) k + i_K dk. Returns (ok, omega, why)."""
    n = len(X)
    nz = lambda v: any(simp(x) != 0 for x in v)
    if not (nz(a) and nz(b)):  # one scalar: d of the nonzero covector must vanish, chi is then irrelevant
        c = a if nz(a) else b
        dc = dform(c, X)
        ok = all(simp(dc[u][v]) == 0 for u in range(n) for v in range(u + 1, n))
        return ok, sp.zeros(n, 1), "one square, d of it" + ("" if ok else " != 0")
    k, m = a + b, a - b
    kk, km, mm = (k.T * k)[0], (k.T * m)[0], (m.T * m)[0]
    det = kk * mm - km**2
    K, Mv = (mm * k - km * m) / det, (kk * m - km * k) / det
    dk, dm = dform(k, X), dform(m, X)
    om = sum((K[u] * Mv[v] * dm[u][v] for u in range(n) for v in range(n)), sp.Integer(0)) * k
    om += sp.Matrix([sum(K[u] * dk[u][v] for u in range(n)) for v in range(n)])
    om = om.applyfunc(simp)
    eqs = [om[u] * k[v] - om[v] * k[u] + dk[u][v] for u in range(n) for v in range(u + 1, n)]
    eqs += [om[u] * m[v] - om[v] * m[u] - dm[u][v] for u in range(n) for v in range(u + 1, n)]
    if any(simp(e) != 0 for e in eqs):
        return False, om, "no chi: omega ^ k = -dk and omega ^ m = dm have no common solution"
    dom = dform(om, X)
    ok = all(simp(dom[u][v]) == 0 for u in range(n) for v in range(u + 1, n))
    return ok, om, "omega closed" if ok else "omega = d chi needed but d omega != 0"


def wedge3(p, a, b, X):
    """Components (p ^ a ^ b)_{mu nu rho}, mu < nu < rho."""
    n = len(X)
    return [sp.Matrix([[p[i], p[j], p[k]], [a[i], a[j], a[k]], [b[i], b[j], b[k]]]).det()
            for i in range(n) for j in range(i + 1, n) for k in range(j + 1, n)]


def canon_trig(x):
    """A normalizer for expressions in the trig or hyperbolic functions of x: multiple angles are expanded and
    sin(x)^2 -> 1 - cos(x)^2, sinh(x)^2 -> cosh(x)^2 - 1, so that equal expressions have one polynomial form."""
    s, co, sh, ch = sp.sin(x), sp.cos(x), sp.sinh(x), sp.cosh(x)

    def red(p):
        p = sp.expand(p)
        if p.has(s):
            p = sp.rem(p, s**2 + co**2 - 1, s)
        if p.has(sh):
            p = sp.rem(p, sh**2 - ch**2 + 1, sh)
        return sp.expand(p)

    def f(e):
        n, d = sp.fraction(sp.together(sp.expand_trig(e)))
        return sp.factor(red(n)) / sp.factor(red(d))
    return f


def two_scalar_form(Rm, gm, X=None, simp=sp.simplify, sampler=None, lam_shift=0, canon=None):
    """R_ab = lam g_ab + a_a a_b - b_a b_b, with a = sqrt(8 pi) d phi_1 (normal) and b = sqrt(8 pi) d phi_2 (phantom).

    Returns a dict with 'status' 'algebraic' (no real lam gives the form), 'form' (the form holds, X not given),
    'closure', 'potential' or 'pass', and lam, a, b, omega where reached. Without X only the algebraic stage runs.
    lam_shift is added to lam in the potential stage only (a control). canon, if given, normalizes the entries."""
    out = {"status": "algebraic", "lam": None, "a": None, "b": None}
    if canon is not None:
        Rm, gm = Rm.applyfunc(canon), gm.applyfunc(canon)
    for lr in candidates(Rm, gm):
        M = (Rm - lr * gm).applyfunc(simp)
        ab = split_form(M, sampler, simp)
        if ab is None:
            continue
        a, b = ab
        if any(simp(x) != 0 for x in M - (a * a.T - b * b.T)):
            continue
        out.update(status="form", lam=lr, a=a, b=b)
        break
    if out["status"] != "form" or X is None:
        return out
    a, b = out["a"], out["b"]
    ok, om, why = closure(a, b, X, simp)
    out.update(omega=om, closure_why=why)
    if not ok:
        out["status"] = "closure"
        return out
    lp = out["lam"] + lam_shift
    dl = sp.Matrix([sp.diff(lp, s) for s in X])
    w3 = [simp(e) for e in wedge3(dl, a, b, X)]
    out["dlam_ab"] = w3
    if any(e != 0 for e in w3):
        out["status"] = "potential"
        return out
    out["status"] = "pass"
    return out


def interval_profiles(rels, W, r):
    """Can relations W0 = g(W1, r), W2 = h(W1, r) between the jets (w, w', w'') hold on an interval of r with w' != 0?

    Along the interval W0' = W1 and W1' = W2. Differentiating W0 = g gives W1 = g_r + g_W1 h, an equation in W1 and r
    alone. If it holds identically, every solution of W1' = h(W1, r) is a profile. Otherwise W1 is one of its roots,
    and a nonzero root must also satisfy W1' = h. Returns (found, notes): the relations and branches that give a
    profile with w' != 0, and for each relation the equation and its nonzero roots."""
    W0, W1, W2 = W
    found, notes = [], []
    for s in rels:
        g, h = s[W0], s[W2]
        e = sp.factor(sp.simplify(W1 - (sp.diff(g, r) + sp.diff(g, W1) * h)))
        if e == 0:
            found.append((s, "every solution of W1' = h"))
            notes.append((e, "any"))
            continue
        br = [b for b in sp.solve(e, W1) if sp.simplify(b) != 0]
        notes.append((e, br))
        found += [(s, b) for b in br if sp.simplify(sp.diff(b, r) - h.subs(W1, b)) == 0]
    return found, notes


# ---------------------------------------------------------------------------------------------------------------
# Random rational Lorentz transformations and planted forms
# ---------------------------------------------------------------------------------------------------------------
eta = sp.diag(-1, 1, 1, 1)


def rq(lo=-9, hi=9):
    return sp.Rational(random.randint(lo, hi), random.randint(1, 9))


def lorentz():
    """A random rational Lorentz matrix: a boost along each spatial axis and a rotation in each spatial plane."""
    L = sp.eye(4)
    for i in (1, 2, 3):
        u = sp.Rational(random.randint(-4, 4), 9)
        ch, sh = (1 + u**2) / (1 - u**2), 2 * u / (1 - u**2)
        B = sp.eye(4)
        B[0, 0] = B[i, i] = ch
        B[0, i] = B[i, 0] = sh
        L = L * B
    for i, j in ((1, 2), (2, 3), (1, 3)):
        u = sp.Rational(random.randint(-9, 9), 7)
        co, si = (1 - u**2) / (1 + u**2), 2 * u / (1 + u**2)
        Q = sp.eye(4)
        Q[i, i] = Q[j, j] = co
        Q[i, j], Q[j, i] = -si, si
        L = L * Q
    assert L.T * eta * L == eta
    return L


def vec4():
    return sp.Matrix([rq() for _ in range(4)])


# ---------------------------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------------------------
def main():
    opts = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
    mut = opts.get("mutate")
    if mut is not None and mut not in MUTATIONS:
        sys.exit(f"unknown mutation {mut!r}. Known: {', '.join(MUTATIONS)}")
    random.seed(65)
    t0 = time.time()
    print(f"two_scalars.py, issue #65. sympy {sp.__version__}, python {sys.version.split()[0]}. mutation: "
          f"{mut or 'none'}" + (f" ({MUTATIONS[mut][1]})" if mut else ""))

    # 1. Planted truth for the algebraic stage, in random Lorentz frames.
    print("\n== 1. Planted truth: R = lam eta + a a - b b passes, other signatures fail ==")
    bad = []
    for _ in range(12):
        L, lv, a, b = lorentz(), rq(), vec4(), vec4()
        sb = 1 if mut == "planted-sign" else -1
        Rm = L.T * (lv * eta + a * a.T + sb * b * b.T) * L
        res = two_scalar_form(Rm, eta)
        ok = res["status"] == "form" and res["lam"] == lv
        bad += [] if ok else [(lv, res["status"], res["lam"])]
    check("planted_pass", not bad, f"12 random lam eta + a a - b b, transformed: lam recovered and the form split "
          f"in all but {len(bad)} {bad[:2]}")
    bad = []
    kinds = [(1, 1, 0), (-1, -1, 0), (1, -1, 1), (1, -1, -1)]  # signs of the squares a a, b b, c c
    for kind in kinds * 3:
        L, lv, a, b, c = lorentz(), rq(), vec4(), vec4(), vec4()
        sa, sb, sc = (1, -1, 0) if mut == "planted-fail-sign" else kind
        Rm = L.T * (lv * eta + sa * a * a.T + sb * b * b.T + sc * c * c.T) * L
        res = two_scalar_form(Rm, eta)
        bad += [] if res["status"] == "algebraic" else [(kind, res["status"], res["lam"])]
    check("planted_fail", not bad, f"12 random forms with signs (+,+), (-,-), (+,-,+), (+,-,-): all fail the "
          f"algebraic stage but {len(bad)} {bad[:2]}")

    # 2. The lemma's prediction on random tensors that satisfy (I1) and (I2).
    print("\n== 2. Lemma: with (I1), (I2), the form holds iff lam = 0 and det R_3 = 0 (R_3 the (0, 1, 2) block) ==")
    wt = (1, -1, 2) if mut == "lemma-weight" else (1, 1, 2)
    wn = lambda u: wt[0] * u[0]**2 + wt[1] * u[1]**2 + wt[2] * u[2]**2
    mism, npass, nfail = [], 0, 0
    for trial in range(24):
        if trial % 2 == 0:  # det R_3 = 0: R_3 = a a - t b b with t fixed by the identity
            a, b = sp.Matrix([rq() for _ in range(3)]), sp.Matrix([rq() for _ in range(3)])
            if wn(b) == 0:
                continue
            R3 = a * a.T - (wn(a) / wn(b)) * b * b.T
        else:  # a generic R_3 with the identity, R_22 fixed by it
            R3 = sp.zeros(3, 3)
            for i in range(3):
                for j in range(i, 3):
                    R3[i, j] = R3[j, i] = rq()
            R3[2, 2] = -(wt[0] * R3[0, 0] + wt[1] * R3[1, 1]) / wt[2]
        Rm = sp.zeros(4, 4)
        Rm[:3, :3] = R3
        L = lorentz()
        res = two_scalar_form(L.T * Rm * L, eta)
        want = R3.det() == 0
        got = res["status"] == "form" and res["lam"] == 0
        npass += got
        nfail += res["status"] == "algebraic"
        if want != got or (not got and res["status"] != "algebraic"):
            mism.append((trial, R3.det(), res["status"], res["lam"]))
    check("lemma_samples", not mism, f"{npass} pass with lam = 0 (det R_3 = 0), {nfail} fail (det R_3 != 0), "
          f"{len(mism)} disagree with the lemma {mism[:2]}")

    # 3. Alcubierre, constant v: det R_3 in the comoving spherical frame, as a cubic in c^2 = cos(theta)^2.
    print("\n== 3. Alcubierre, constant v: D = det R_3 on (E_0 = n, E_1 = d_x, e_rho), general profile ==")
    Rsph = ps.sph_tetrad_ricci("alc-drop-w2" if mut == "alc-drop-w2" else None)
    W0, W1, W2 = ps.W
    r, th = ps.r_, ps.th
    c = sp.Symbol("c", real=True)
    R3 = Rsph[:3, :3]
    print(f"  ({time.time() - t0:.0f} s) identities from phantom_source.py: R(e_varphi, .) = "
          f"{[sp.simplify(Rsph[i, 3]) for i in range(4)]}, R_00 + R_11 + 2 R_rho,rho = "
          f"{sp.simplify(R3[0, 0] + R3[1, 1] + 2 * R3[2, 2])}")
    S = sp.Symbol("S", real=True)  # S = sin(theta), reduced by S^2 = 1 - c^2
    Dth = sp.expand(sp.simplify(R3.det()) * 16 * r**3)
    D = sp.expand(sp.expand_trig(Dth).subs({sp.sin(th): S, sp.cos(th): c}))
    D = sp.expand(sp.rem(sp.Poly(D, S).as_expr(), S**2 + c**2 - 1, S))
    P0, P1 = sp.Poly(D.coeff(S, 0), c), sp.Poly(D.coeff(S, 1), c)
    odd = [f"c^{k}" for (k,) in P0.monoms() if k % 2] + [f"S c^{k}" for (k,) in P1.monoms() if not P1.is_zero]
    coef = {f"c^{k}": sp.factor(P0.coeff_monomial(c**k)) for k in range(0, P0.degree() + 1)
            if P0.coeff_monomial(c**k) != 0}
    coef.update({f"S c^{k}": sp.factor(P1.coeff_monomial(c**k)) for (k,) in P1.monoms() if not P1.is_zero})
    print(f"  16 r^3 D is a polynomial of degree {P0.degree()} in c = cos(theta), odd terms {odd}. Coefficients:")
    for k, v in coef.items():
        print(f"    {k}: {v}")
    on_axis = [sp.factor(D.subs({S: 0, c: cv})) for cv in (1, -1)]
    check("det_axis", on_axis == [0, 0] and not odd, f"D on the axis (c = 1 and c = -1) = {on_axis}, odd terms {odd}")
    zz = sp.Symbol("zz")
    G = sp.groebner(list(coef.values()) + [1 - zz * W1 * r], zz, W0, W2, W1, r, order="lex")
    Gs = [g for g in G.exprs if not g.has(zz)]
    print(f"  Groebner basis (lex, zz W1 r = 1 enforces W1 != 0, r != 0), elements free of zz: "
          f"{[sp.factor(g) for g in Gs]}")
    want = [2 * W0 - r * W1, (W1 + r * W2)**2]
    contains = [G.contains(p) for p in want]
    sol = {W0: r * W1 / 2, W2: -W1 / r}
    kills = [sp.simplify(v.subs(sol)) for v in coef.values()]
    check("det_jets", all(contains) and all(x == 0 for x in kills),
          f"2 W0 - r W1 and (W1 + r W2)^2 lie in the ideal {contains}, and W0 = r W1 / 2, W2 = -W1 / r kill every "
          f"coefficient {kills}. So D(r, .) = 0 for all theta, with w'(r) != 0, only where 2 w = r w' and "
          f"r w'' = -w'")

    # The ODE step: on an interval of r, those two relations force w' = 0.
    rel = sp.solve(Gs, [W0, W2], dict=True)
    print(f"  relations solved from the basis: {rel}")
    if mut == "ode-sign":
        rel = [{W0: s[W0], W2: -s[W2]} for s in rel]
    found, notes = interval_profiles(rel, (W0, W1, W2), r)
    for e, br in notes:
        print(f"  differentiating W0 = g(W1, r) along r with W2 = h(W1, r) leaves {e} = 0, so W1 is one of {br}")
    check("no_open_set", not found, f"relations that hold on an interval of r with w' != 0: {found or 'none'}. "
          "So D vanishes on no open set of the wall")

    # 4. The test itself at exact points of Alcubierre's wall, through the same code as the control.
    print("\n== 4. Two-scalar test at exact points (lab route of phantom_source.py, same code as the control) ==")
    Rt = ps.lab_tetrad_ricci(None)
    random.seed(57)  # phantom_source.py's sample points; its random jets differ, since its seed feeds more draws first
    samples = []
    for (x0, y0, z0, r0) in [(12, 3, 4, 13), (-12, 4, 3, 13), (13, 0, 0, 13), (-13, 0, 0, 13), (0, 5, 12, 13)]:
        p = [sp.Rational(q, 13) for q in (x0, y0, z0, r0)]
        for vv in (sp.Integer(1), sp.Rational(1, 2), sp.Integer(3)):
            samples.append(("tanh sig=8 Rb=1", p, vv, ps.tanh_jets(p[3], 8, 1)))
    for (x0, y0, z0, r0) in [(4, 3, 0, 5), (-3, 0, 4, 5), (5, 0, 0, 5), (0, 3, 4, 5)]:
        p = [sp.Rational(q, 5) for q in (x0, y0, z0, r0)]
        samples.append(("bump Rb=2", p, sp.Rational(7, 4), ps.bump_jets(p[3], 2)))
    for (x0, y0, z0, r0) in [(12, 3, 4, 13), (-12, 4, 3, 13), (4, 3, 0, 5), (-3, 0, 4, 5)]:
        p = [sp.Rational(q, 7) for q in (x0, y0, z0, r0)]
        samples.append(("random jets", p, sp.Rational(random.randint(1, 9), 3),
                        [sp.Rational(random.randint(-9, 9), random.randint(1, 9)) for _ in range(3)]))
    random.seed(650)
    out, axis_info = [], []
    for (lab_, p, vv, jets) in samples:
        sub = {ps.X: p[0], ps.y: p[1], ps.z: p[2], ps.R: p[3], ps.v: vv, ps.F[0]: jets[0], ps.F[1]: jets[1],
               ps.F[2]: jets[2]}
        Rm = Rt.subs(sub)
        axis = p[1] == 0 and p[2] == 0
        if mut == "alc-planted" and not axis:
            a, b = vec4(), vec4()
            Rm = a * a.T - b * b.T
        Dv = D.subs({W0: vv * (1 - jets[0]), W1: -vv * jets[1], W2: -vv * jets[2], r: p[3], c: p[0] / p[3],
                     S: sp.sqrt(p[1]**2 + p[2]**2) / p[3]})
        res = two_scalar_form(Rm, eta)
        ok = (res["status"] == "form" and res["lam"] == 0) if axis else res["status"] == "algebraic"
        ok = ok and ((Dv == 0) == axis)
        out.append((lab_, axis, res["status"], ok))
        if axis and res["status"] == "form" and len(axis_info) < 2:
            ab = (res["a"] * res["a"].T - res["b"] * res["b"].T)
            ev = [sp.N(e, 6) for e in (eta * ab).eigenvals()]
            axis_info.append(f"x - v t = {p[0]}, v = {vv} ({lab_}): a = {[sp.N(x, 5) for x in res['a']]}, "
                             f"b = {[sp.N(x, 5) for x in res['b']]}, eigenvalues of eta (a a - b b): {ev}")
    for key in sorted(set((o[0], o[1]) for o in out)):
        sts = [o[2] for o in out if (o[0], o[1]) == key]
        print(f"  {key[0]:>16}, {'on the axis' if key[1] else 'off the axis'}: {len(sts)} points, outcomes "
              f"{sorted(set(sts))}")
    for s in axis_info:
        print(f"  info: {s}")
    nax = sum(o[1] for o in out)
    check("alc_points", all(o[3] for o in out), f"{len(out) - nax} off-axis points fail the algebraic stage with "
          f"D != 0, {nax} on-axis points pass it with lam = 0 and D = 0; {sum(not o[3] for o in out)} disagree")

    # 5. Control: Ellis-Bronnikov with an added normal scalar, ultrastatic, with a constant potential.
    print("\n== 5. Control: ds^2 = -dt^2 + dl^2 + r(l)^2 dOmega^2, r^2 = (1 - beta cos 2 kappa l) / (2 kappa^2) ==")
    t, l, the, phv = sp.symbols("t l theta phi", real=True)
    ka, be = sp.symbols("kappa beta", positive=True)
    Xe = [t, l, the, phv]
    if mut == "eb2-ads":
        r2 = (-1 + be * sp.cosh(2 * ka * l)) / (2 * ka**2)
        brange = (1, 3)
    else:
        r2 = (1 - be * sp.cos(2 * ka * l)) / (2 * ka**2)
        brange = (0, sp.Rational(1, 2))
        if mut == "eb2-deform":
            r2 += (1 - be) * sp.cos(4 * ka * l) / (20 * ka**2)

    def sampler():
        u = random.randint(1, 99)
        return {ka: sp.Rational(random.randint(1, 30), 10), be: brange[0] + (brange[1] - brange[0]) * sp.Rational(u, 100),
                l: sp.Rational(random.randint(-50, 50), 10), the: sp.Rational(random.randint(1, 30), 10),
                t: sp.Rational(random.randint(-50, 50), 10)}
    ge = sp.diag(-1, 1, r2, r2 * sp.sin(the)**2)
    Re = ps.ricci(ge, Xe, lambda E, k: sp.diff(E, Xe[k]), sp.simplify)
    print(f"  nonzero R_ab: " + ", ".join(f"R_{Xe[i]}{Xe[j]} = {sp.factor(Re[i, j])}" for i in range(4)
                                          for j in range(i, 4) if Re[i, j] != 0))
    res = two_scalar_form(Re, ge, Xe, sampler=sampler, canon=canon_trig(ka * l),
                          lam_shift=ka**2 * sp.cos(the) / 10 if mut == "eb2-potential" else 0)
    if mut == "eb2-unclosed" and res["status"] != "algebraic":
        res["b"] = res["b"] * (1 + sp.cos(the)**2)
        ok, om, why = closure(res["a"], res["b"], Xe, sp.simplify)
        res.update(status="pass" if ok else "closure", omega=om, closure_why=why)
    print(f"  status {res['status']}, lam = {res['lam']}, a = {list(res['a']) if res['a'] is not None else None}, "
          f"b = {list(res['b']) if res['b'] is not None else None}")
    check("eb2_algebraic", res["status"] != "algebraic",
          f"lam = {res['lam']} (2 kappa^2 for the solution), one normal and one phantom square: "
          f"{res['status'] != 'algebraic'}")
    if res["status"] != "algebraic":
        check("eb2_closed", res["status"] not in ("closure",), f"omega = {list(res['omega'])}: {res['closure_why']}")
    if res["status"] in ("potential", "pass"):
        check("eb2_potential", res["status"] == "pass", f"d lam ^ a ^ b = {res['dlam_ab']}")
    if res["status"] == "pass":
        # The fields. omega = 0 here, so chi is constant and a, b are themselves the gradients (times sqrt(8 pi)).
        q = sp.sqrt((1 + be) / (1 - be))
        phi1, phi2 = ka * t / (2 * sp.sqrt(sp.pi)), -sp.atan(q * sp.tan(ka * l)) / (2 * sp.sqrt(sp.pi))
        u1, u2 = res["a"] / sp.sqrt(8 * sp.pi), res["b"] / sp.sqrt(8 * sp.pi)
        # a and b are fixed up to sign (phi_i -> -phi_i), so the sign of each phi_i is matched to them. beta = 1 - om
        # with om > 0 lets sympy take the square roots of 1 - beta.
        om = sp.Symbol("omega", positive=True)
        zero = lambda e: sp.simplify(sp.simplify(e).subs(be, 1 - om)) == 0
        sg = lambda ph, u: 1 if all(zero(sp.diff(ph, s) - u[i]) for i, s in enumerate(Xe)) else -1
        phi1, phi2 = sg(phi1, u1) * phi1, sg(phi2, u2) * phi2
        grad = [0 if zero(sp.diff(phi1, s) - u1[i]) else sp.diff(phi1, s) - u1[i] for i, s in enumerate(Xe)]
        grad += [0 if zero(sp.diff(phi2, s) - u2[i]) else sp.diff(phi2, s) - u2[i] for i, s in enumerate(Xe)]
        V = res["lam"] / (8 * sp.pi)
        print(f"  phi_1 = {phi1}, phi_2 = {phi2} (|kappa l| < pi/2), V = {V}; d phi_i - u_i = {grad}")
        if mut == "eb2-residual":
            u2 = u2 * sp.Rational(11, 10)
        gi = ge.inv()
        kin = lambda u: sp.simplify((u.T * gi * u)[0])
        T = u1 * u1.T - u2 * u2.T - ge * (kin(u1) / 2 - kin(u2) / 2 + V)
        Rsc = sp.simplify(sum(gi[i, j] * Re[i, j] for i in range(4) for j in range(4)))
        resid = (Re - ge * Rsc / 2 - 8 * sp.pi * T).applyfunc(sp.simplify)
        check("eb2_residual", all(x == 0 for x in grad) and all(x == 0 for x in resid),
              f"d phi_i = u_i: {all(x == 0 for x in grad)}; G_ab - 8 pi T_ab = "
              f"{'0' if all(x == 0 for x in resid) else list(resid)}")
    print(f"  ({time.time() - t0:.0f} s)")

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

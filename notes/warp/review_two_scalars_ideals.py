#!/usr/bin/env python3
"""Independent review of #65, part I, for issue #79: the constant-speed coefficients, the ideals, the step to "no open
set", the accelerating bubble and a shift b(rho).

Notation. A bubble of constant speed has b = w(r) - v with w = v (1 - f), r the distance to the centre, c = cos theta,
s = sin theta, w1 = w', w2 = w''. P, Q, Z, Y, h are the blocks of R_3 in the frame (n, d_x, e_rho), from
review_two_scalars_cartan.py. Here they are derived again by differentiating b(xi, y, z) with my own jet code, and
compared with the imported Cartan routine of review_cartan.py (Jets, alcubierre_shift, cartan_ricci). Ideal statements
are decided by explicit cofactors found by linear algebra over Q(r, w1, ...) (no Groebner basis) and verified by exact
expansion, and by a resultant.

    python3 notes/warp/review_two_scalars_ideals.py                  checks, then every control
    python3 notes/warp/review_two_scalars_ideals.py --no-controls    checks only
    python3 notes/warp/review_two_scalars_ideals.py --mutate=NAME    one control (its target checks must fail, exit 1)
"""
import importlib.util
import os
import subprocess
import sys
import time

import sympy as sp
from sympy.polys.matrices import DomainMatrix

HERE = os.path.dirname(os.path.abspath(__file__))
OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate", "")
CONTROLS = {
    "alc-lapse": ("imported Cartan with lapse 1 + b", ["r3_const", "acc_R3"]),
    "coef-c0-sign": ("expected c_0 with the opposite sign", ["coeffs_const"]),
    "ideal-cofactor": ("one cofactor of the certificate altered", ["ideal_const"]),
    "planted-ideal": ("planted jets with A = 1, expecting D = 0", ["zero_set_const"]),
    "ode-sign": ("second relation r w2 = +w1", ["no_open_set"]),
    "acc-no-A": ("my comoving R_3 without the A terms", ["acc_R3"]),
    "acc-drop-odd": ("odd coefficients dropped from the certificate", ["acc_unit"]),
    "brho-h2": ("expected h doubled", ["brho_D"]),
}
if MUT and MUT not in CONTROLS:
    sys.exit(f"unknown mutation {MUT}")

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def load_review_cartan():
    saved = sys.argv
    sys.argv = [saved[0]]
    try:
        spec = importlib.util.spec_from_file_location("review_cartan", os.path.join(HERE, "review_cartan.py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    finally:
        sys.argv = saved
    return mod


rc = load_review_cartan()

# symbols
r, c, s = sp.symbols("r c s", positive=True)
w, w1, w2, aa = sp.symbols("w w1 w2 aa")
Asym, Bsym, om = sp.symbols("A B omega")
u = sp.Symbol("u")


# ---------------------------------------------------------------------------------------------------------------
# My own jets for a function of the distance, b = W(r_s) - const: symbols W_k for the k-th derivative at r, and the
# total derivative in (xi, y, z), with dr/dx_i = x_i / r. r is an independent symbol with that derivative rule.
# ---------------------------------------------------------------------------------------------------------------
xi, y, z = sp.symbols("xi y z", real=True)
rs = sp.Symbol("r_s", positive=True)
Wk = [sp.Symbol("W%d" % k) for k in range(5)]
COORD = (xi, y, z)


def Dtot(e, i):
    """d/dx_i of an expression in (xi, y, z, r_s, W_k)."""
    out = sp.diff(e, COORD[i]) + sp.diff(e, rs) * COORD[i] / rs
    for k in range(4):
        out += sp.diff(e, Wk[k]) * Wk[k + 1] * COORD[i] / rs
    return out


def blocks_comoving(Aval):
    """P, Q, Z, Y, h at the slice y = rho, z = 0 (rho a symbol), for b = W(r_s) - X1 with shift speed W0 and
    d_t|_xi g = Aval * g for every spatial derivative g of b, so n(g) = Aval * g - W0 * d_xi g."""
    bs = Wk[0]  # b up to the constant X1, whose spatial derivatives vanish
    nn = lambda g: Aval * g - Wk[0] * Dtot(g, 0)
    bx = Dtot(bs, 0)
    by, bz = Dtot(bs, 1), Dtot(bs, 2)
    beta = y * by + z * bz  # rho b_rho
    P = nn(bx) - bx**2
    lap = Dtot(by, 1) + Dtot(bz, 2)
    Q = -lap / 2
    Zr = Dtot(beta, 0) / 2
    Yr = -nn(beta) / 2 + bx * beta
    q = y**2 + z**2
    h = beta**2 / (2 * q)
    return P, Q, Zr, Yr, h


def to_slice(e, rho):
    """substitute the slice y = rho, z = 0, xi = r c, rho = r s, r_s = r, W_k = w_k and reduce s^2 = 1 - c^2."""
    e = e.subs({z: 0}).subs({y: rho})
    e = e.subs({xi: r * c, rs: r})
    e = e.subs({rho: r * s})
    e = e.subs({Wk[0]: w, Wk[1]: w1, Wk[2]: w2})
    return sp.together(e)


def reduce_s(e):
    """rewrite a rational function of r, c, s, polynomial in s, with s^2 = 1 - c^2."""
    e = sp.expand(sp.cancel(e))
    num, den = sp.fraction(sp.together(e))
    num = sp.expand(num)
    pn = sp.Poly(num, s)
    out = 0
    for (k,), co in pn.terms():
        out += co * (1 - c**2) ** (k // 2) * s ** (k % 2)
    return sp.expand(out) / den


def D_of(P, Q, Z, Y, h):
    return h * (h**2 - P**2 - Q**2 + Y**2 + Z**2) - P * (Y**2 - Z**2) + 2 * Q * Y * Z


def main():
    t0 = time.time()
    print("review_two_scalars_ideals.py, issue #79. sympy %s, python %s. mutation: %s"
          % (sp.__version__, sys.version.split()[0], MUT or "none"))
    print("imports from review_cartan.py: Jets, alcubierre_shift, cartan_ricci")
    rho = sp.Symbol("rho", positive=True)

    # ---- 1. R_3 for constant speed: my derivation against the imported Cartan Ricci and the note's table --------
    print("\n== 1. R_3 for constant v: blocks from my jet code, the imported Cartan Ricci, and the note's table ==")
    P0, Q0, Zr0, Yr0, h0 = blocks_comoving(0)
    # slice blocks as functions of r, c, s: Z = Zr / rho, Y = Yr / rho
    Pm, Qm, hm = (reduce_s(to_slice(e, rho)) for e in (P0, Q0, h0))
    Zm = reduce_s(to_slice(Zr0, rho) / (r * s))
    Ym = reduce_s(to_slice(Yr0, rho) / (r * s))
    table = {  # phantom_source.md, "The Ricci tensor in the frame"
        "00": -(2 * w * w1 * s**2 + 2 * r * w * w2 * c**2 + r * w1**2 * (1 + c**2)) / (2 * r),
        "01": -(w1 * (1 + c**2) + r * w2 * s**2) / (2 * r),
        "0r": (r * w2 - w1) * s * c / (2 * r),
        "11": (2 * w * w1 * s**2 + 2 * r * w * w2 * c**2 + r * w1**2 * (3 * c**2 - 1)) / (2 * r),
        "1r": (w * (r * w2 - w1) + 2 * r * w1**2) * s * c / (2 * r),
        "rr": w1**2 * s**2 / 2,
    }
    mine = {"00": Pm - hm, "01": Qm, "0r": Zm, "11": -Pm - hm, "1r": Ym, "rr": hm}
    d_table = {k: sp.expand(reduce_s(mine[k] - table[k])) for k in table}
    ok_table = all(v == 0 for v in d_table.values())

    # Cartan, lab coordinates, b = -v H(q), q = (x - v t)^2 + y^2 + z^2 (review_cartan.alcubierre_shift)
    t, x, yy, zz = sp.symbols("t x y z", real=True)
    J = rc.Jets((t, x, yy, zz))
    v = sp.Symbol("v", positive=True)
    b, Hs = rc.alcubierre_shift(J, v)
    lapse = 1 + b if MUT == "alc-lapse" else 1
    E = sp.Matrix([[lapse, 0, 0, 0], [b, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    Ric, bad, _ = rc.cartan_ricci(E, J.total)
    slc = {t: 0, x: xi, yy: rho, zz: 0}
    Rs = Ric.subs(slc)
    H0, H1, H2 = Hs(0), Hs(1), Hs(2)
    # f(r) = H(r^2): w = v(1 - f), w1 = -v f', w2 = -v f''; r = sqrt(xi^2 + rho^2), c = xi / r, s = rho / r
    rr_ = sp.sqrt(xi**2 + rho**2)
    jets = {w: v * (1 - H0), w1: -2 * v * rr_ * H1, w2: -v * (2 * H1 + 4 * rr_**2 * H2)}
    back = {r: rr_, c: xi / rr_, s: rho / rr_}
    ent = {"00": Rs[0, 0], "01": Rs[0, 1], "0r": Rs[0, 2], "11": Rs[1, 1], "1r": Rs[1, 2], "rr": Rs[2, 2]}
    d_cartan = {k: sp.simplify(sp.expand(sp.expand(mine[k]).subs(jets).subs(back) - ent[k])) for k in ent}
    ok_cartan = all(v_ == 0 for v_ in d_cartan.values()) and not bad
    check("r3_const", ok_table and ok_cartan, "R_3 entries from my jet code equal the note's table (nonzero: %s) and the "
          "imported Cartan Ricci of b = -v H((x - vt)^2 + y^2 + z^2) (nonzero: %s; structure-equation violations %d)"
          % ([k for k, v_ in d_table.items() if v_ != 0], [k for k, v_ in d_cartan.items() if v_ != 0], len(bad)))

    # ---- 2. The coefficients of 16 r^3 D ------------------------------------------------------------------------
    print("\n== 2. 16 r^3 D as a polynomial in c^2 for constant v ==")
    Dfull = reduce_s(D_of(Pm, Qm, Zm, Ym, hm))
    D16 = sp.expand(sp.cancel(16 * r**3 * Dfull))
    pc = sp.Poly(D16, c)
    cks = [sp.factor(pc.coeff_monomial(c ** (2 * k))) for k in range(4)]
    odd = [pc.coeff_monomial(c ** (2 * k + 1)) for k in range(3)]
    printed_c0 = 2 * w1**2 * r * (w1**4 * r**2 - 4 * w**2 * w1**2 - (w1 + r * w2) ** 2)
    printed_c6 = -2 * (2 * r * w * w2 - 2 * w * w1 + r * w1**2) * (r * w * w2 - w * w1 + r * w1**2) ** 2
    sgn = -1 if MUT == "coef-c0-sign" else 1
    ok = (sp.expand(cks[0] - sgn * printed_c0) == 0 and sp.expand(cks[3] - printed_c6) == 0 and pc.degree() == 6
          and all(o == 0 for o in odd))
    check("coeffs_const", ok, "16 r^3 D has degree %d in c, no odd powers; c^0 and c^6 equal the printed ones" % pc.degree())
    print("  c^2 coefficient:", cks[1])
    print("  c^4 coefficient:", cks[2])
    ck = [sp.expand(e) for e in cks]

    # ---- 3. The ideal of the four coefficients: A = 2w - r w1, B = w1 + r w2 -------------------------------------
    print("\n== 3. ideal(c_0, c_1, c_2, c_3) = ideal(A, B^2) in Q[r, w, w1, w2][1/(r w1)], by cofactors ==")
    sub = {w: (Asym + r * w1) / 2, w2: (Bsym - w1) / r}
    cAB = [sp.expand(sp.cancel(e.subs(sub))) for e in ck]
    # (a) c_k in ideal(A, B^2): split monomials
    inside, lost = [], []
    for k, e in enumerate(cAB):
        alpha, beta = 0, 0
        for (i, j), co in sp.Poly(e, Asym, Bsym).terms():
            if i >= 1:
                alpha += co * Asym ** (i - 1) * Bsym**j
            elif j >= 2:
                beta += co * Bsym ** (j - 2)
            else:
                lost.append((k, i, j))
        inside.append((alpha, beta))
    ok_in = not lost and all(sp.expand(inside[k][0] * Asym + inside[k][1] * Bsym**2 - cAB[k]) == 0 for k in range(4))
    # (b) A and B^2 in ideal(c_k): cofactors by linear algebra over Q(r, w1)
    K = sp.QQ.frac_field(r, w1)

    def find_cert(polys, target, d):
        mons = [(i, j) for i in range(d + 1) for j in range(d + 1 - i)]
        unknown = [(k, m) for k in range(len(polys)) for m in mons]
        cols, allm = [], {(0, 0)}
        for (k, m) in unknown:
            P_ = sp.Poly(sp.expand(Asym ** m[0] * Bsym ** m[1] * polys[k]), Asym, Bsym)
            dct = dict(P_.terms())
            cols.append(dct)
            allm |= set(dct)
        tdct = dict(sp.Poly(target, Asym, Bsym).terms())
        allm |= set(tdct)
        allm = sorted(allm)
        rows = [[K.from_sympy(sp.sympify(cdict.get(mon, 0))) for cdict in cols] + [K.from_sympy(sp.sympify(tdct.get(mon, 0)))]
                for mon in allm]
        R_, piv = DomainMatrix(rows, (len(allm), len(unknown) + 1), K).rref()
        if len(unknown) in piv:
            return None
        Rl = R_.to_Matrix()
        sol = {}
        for row, pcol in enumerate(piv):
            val = Rl[row, len(unknown)]
            if val != 0:
                sol[unknown[pcol]] = sp.factor(val)
        return sol

    def verify(polys, target, sol):
        tot = 0
        for (k, m), co in sol.items():
            tot += co * Asym ** m[0] * Bsym ** m[1] * polys[k]
        return sp.expand(sp.cancel(sp.together(tot - target))) == 0

    certs, ok_out, dens_ok, used_d = {}, True, True, {}
    for name, tgt in (("A", Asym), ("B^2", Bsym**2)):
        sol = None
        for d in range(3, 7):
            sol = find_cert(cAB, tgt, d)
            if sol is not None:
                used_d[name] = d
                break
        if sol is None:
            ok_out = False
            continue
        if MUT == "ideal-cofactor":  # control: alter one cofactor
            key = sorted(sol)[0]
            sol[key] = sol[key] + 1
        ok_out &= verify(cAB, tgt, sol)
        dens_ok &= all(sp.Poly(sp.denom(sp.together(co)), r, w1).is_monomial for co in sol.values())
        certs[name] = sol
    check("ideal_const", ok_in and ok_out and dens_ok,
          "each c_k = alpha_k A + beta_k B^2 (lost monomials %s); A and B^2 are combinations of the c_k with cofactors of degree "
          "%s in (A, B), %s nonzero cofactor coefficients, denominators monomials in r, w1: %s, verified by exact expansion"
          % (lost, used_d, {k: len(v) for k, v in certs.items()}, dens_ok))
    print("  c_0 = %s" % sp.factor(cAB[0]))
    if "certs" in OPTS:
        print("  CERTIFICATES (constant speed). c_k in (A, B, r, w1); sum_k cofactor(k, A^i B^j) A^i B^j c_k = target")
        print("  c_k:", [str(e) for e in cAB])
        for name, sol in certs.items():
            print("  target %s:" % name)
            for key in sorted(sol):
                print("    (k=%d, A^%d B^%d): %s" % (key[0], key[1][0], key[1][1], sol[key]))

    # ---- 4. Planted jets ------------------------------------------------------------------------------------------
    print("\n== 4. Planted jets: A = B = 0 gives D = 0 for every c; A = 0, B != 0 or B = 0, A != 0 gives D != 0 ==")
    import random
    rng = random.Random(7916)
    Dpoly = sp.expand(D16 / 1)
    cvals = [sp.Rational(rng.randint(-9, 9), 10) for _ in range(10)]

    def jets_AB(Av, Bv, rv=sp.Integer(2), w1v=sp.Integer(3)):
        return {r: rv, w1: w1v, w: (Av + rv * w1v) / 2, w2: (Bv - w1v) / rv}

    Aval_pass = 1 if MUT == "planted-ideal" else 0
    j00 = jets_AB(sp.Integer(Aval_pass), sp.Integer(0))
    z_ok = all(Dpoly.subs(j00).subs(c, cv) == 0 for cv in cvals)
    j_a = jets_AB(sp.Integer(0), sp.Rational(3, 2))
    j_b = jets_AB(sp.Rational(5, 2), sp.Integer(0))
    n_a = sum(Dpoly.subs(j_a).subs(c, cv) != 0 for cv in cvals)
    n_b = sum(Dpoly.subs(j_b).subs(c, cv) != 0 for cv in cvals)
    check("zero_set_const", z_ok and n_a >= 8 and n_b >= 8,
          "A = %d, B = 0 at (r, w1) = (2, 3): D = 0 at 10 random c: %s; A = 0, B = 3/2: D != 0 at %d of 10; B = 0, A = 5/2: "
          "D != 0 at %d of 10" % (Aval_pass, z_ok, n_a, n_b))

    # ---- 5. The step to "no open set" -----------------------------------------------------------------------------
    print("\n== 5. A = 0 and B = 0 on an interval force w' = 0 ==")
    Wf = sp.Function("W")
    rr = sp.Symbol("r", positive=True)
    Aexpr = 2 * Wf(rr) - rr * Wf(rr).diff(rr)
    Bexpr = Wf(rr).diff(rr) + rr * Wf(rr).diff(rr, 2)
    dA = sp.simplify(sp.diff(Aexpr, rr) - (2 * Wf(rr).diff(rr) - Bexpr))
    relB = (Wf(rr).diff(rr) - rr * Wf(rr).diff(rr, 2)) if MUT == "ode-sign" else Bexpr
    solA = sp.dsolve(sp.Eq(Aexpr, 0), Wf(rr))
    solB = sp.dsolve(sp.Eq(relB, 0), Wf(rr))
    print("  dsolve(A = 0):", solA, "  dsolve(second relation = 0):", solB)
    C1, C2 = sp.symbols("C1 C2")
    fa, fb = solA.rhs, solB.rhs
    cc, kk, ll = sp.symbols("cc kk ll")
    # intersection: W = cc * f_A-shape = kk * (ln r or r^2) + ll as functions of r; compare via the Wronskian
    shape_A = sp.simplify(fa.subs(C1, 1))
    free_B = [sp.simplify(fb.diff(sy)) for sy in (C1, C2) if fb.has(sy)]
    fam = [shape_A] + free_B
    wr = sp.simplify(sp.wronskian(fam, rr)) if len(fam) > 1 else 1
    # the functions of A's family and B's family are linearly independent iff the Wronskian is nonzero; then the only
    # common member is 0
    lin_indep = wr != 0 and len(set(fam)) == len(fam)
    check("no_open_set", dA == 0 and lin_indep,
          "dA/dr - (2 w1 - B) = %s; A = 0 gives W = %s, the second relation gives W = %s; the functions %s have Wronskian %s, "
          "so the only common solution is W = 0 and then w' = 0" % (dA, fa, fb, fam, wr))

    # ---- 6. The accelerating bubble: my comoving R_3 against the imported Cartan Ricci in lab coordinates ---------
    print("\n== 6. Accelerating bubble b = -V(t) H((x - X(t))^2 + y^2 + z^2): comoving R_3 with the A terms ==")
    J2 = rc.Jets((t, x, yy, zz))
    jc = {}

    def sym(name, k):
        key = (name, k)
        if key not in jc:
            jc[key] = sp.Symbol("%s_%d" % (name, k), real=True)
            jc[jc[key]] = key
        return jc[key]

    qlab = (x - sym("X", 0)) ** 2 + yy**2 + zz**2

    def deriv(sm):
        name, k = jc[sm]
        if name in ("V", "X"):
            return [sym(name, k + 1), 0, 0, 0]
        return [sym("H", k + 1) * J2.total(qlab, mu) for mu in range(4)]

    J2.rules.append((lambda sm: sm in jc, deriv))
    sym("V", 0), sym("X", 0), sym("H", 0)
    b2 = -sym("V", 0) * sym("H", 0)
    E2 = sp.Matrix([[1 + b2 if MUT == "alc-lapse" else 1, 0, 0, 0], [b2, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    Ric2, bad2, _ = rc.cartan_ricci(E2, J2.total)
    V0, V1, X1 = sym("V", 0), sym("V", 1), sym("X", 1)
    H0_, H1_, H2_ = sym("H", 0), sym("H", 1), sym("H", 2)
    slc2 = {t: 0, x: xi + sym("X", 0), yy: rho, zz: 0}
    Rs2 = Ric2.subs(slc2)
    # mine, with A = V1 / V0: comoving shift speed w~ = X1 - V0 H0, w~1 = -2 V0 r H1, w~2 = -V0 (2 H1 + 4 r^2 H2)
    AAs = sp.Symbol("AA")
    Aacc = 0 if MUT == "acc-no-A" else V1 / V0
    Pa, Qa, Zra, Yra, ha = blocks_comoving(AAs)
    jets2 = {Wk[0]: X1 - V0 * H0_, Wk[1]: -2 * V0 * rr_ * H1_, Wk[2]: -V0 * (2 * H1_ + 4 * rr_**2 * H2_), AAs: Aacc}

    def at_slice(e):
        e = e.subs({z: 0}).subs({y: rho}).subs({xi: xi, rs: rr_})
        return e.subs(jets2)

    Pn, Qn, hn = (at_slice(e) for e in (Pa, Qa, ha))
    Zn, Yn = at_slice(Zra) / rho, at_slice(Yra) / rho
    R3mine = {"00": Pn - hn, "01": Qn, "0r": Zn, "11": -Pn - hn, "1r": Yn, "rr": hn}
    ent2 = {"00": Rs2[0, 0], "01": Rs2[0, 1], "0r": Rs2[0, 2], "11": Rs2[1, 1], "1r": Rs2[1, 2], "rr": Rs2[2, 2]}
    nz2 = [k for k in ent2 if sp.simplify(sp.expand(R3mine[k] - ent2[k])) != 0]
    check("acc_R3", not nz2 and not bad2, "comoving R_3 (P and Y shifted by +A w1 c and -A w1 s / 2, A = V'/V) equals the imported "
          "Cartan Ricci in lab coordinates, generic jets of V, X, H; entries that differ: %s" % nz2)

    # ---- 7. The seven coefficients -----------------------------------------------------------------------------
    print("\n== 7. 16 r^3 D for the accelerating bubble: seven coefficients ==")
    Pa_, Qa_, Za_, Ya_, ha_ = blocks_comoving(aa)
    Pm2, Qm2, hm2 = (reduce_s(to_slice(e, rho)) for e in (Pa_, Qa_, ha_))
    Zm2 = reduce_s(to_slice(Za_, rho) / (r * s))
    Ym2 = reduce_s(to_slice(Ya_, rho) / (r * s))
    D16a = sp.expand(sp.cancel(16 * r**3 * reduce_s(D_of(Pm2, Qm2, Zm2, Ym2, hm2))))
    pa = sp.Poly(D16a, c)
    ca = [sp.expand(pa.coeff_monomial(c**k)) for k in range(7)]
    odd_fac = [sp.cancel(ca[k] / (4 * aa * r * w1)) for k in (1, 3, 5)]
    odd_ok = all(sp.denom(sp.together(f)) == 1 for f in odd_fac)
    at0 = [sp.expand(ca[k].subs(aa, 0) - (ck[k // 2] if k % 2 == 0 else 0)) == 0 for k in range(7)]
    check("acc_coeffs", pa.degree() == 6 and odd_ok and all(at0),
          "16 r^3 D has degree %d in c; the odd coefficients c^1, c^3, c^5 are divisible by 4 aa r w1: %s; at aa = 0 the seven "
          "reduce to section 2's four" % (pa.degree(), odd_ok))
    print("  c^6 = %s" % sp.factor(ca[6]))
    print("  c^5 / (4 aa r w1) = %s" % sp.factor(odd_fac[2]))

    # ---- 8. The unit ideal, with aa r w1 inverted ---------------------------------------------------------------
    print("\n== 8. The seven coefficients generate the unit ideal in Q[r, w1, aa, 1/(aa r w1)][w, w2] ==")
    ek = [ca[k] if k % 2 == 0 else sp.expand(sp.cancel(ca[k] / (4 * aa * r * w1))) for k in range(7)]
    w2om = (w1 - r * w1**2 * om) / r  # omega = (w1 - r w2)/(r w1^2)
    ekw = [sp.expand(sp.cancel(p.subs(w2, w2om))) for p in ek]
    Lam = r * w * w2 + r * w1**2 - w * w1
    rho1 = r * w1**2
    chain = [
        sp.expand(ek[5] - Lam * (3 * Lam - rho1)) == 0,
        sp.expand(ek[6] + 2 * Lam**2 * (2 * Lam - rho1)) == 0,
        sp.expand(ek[6] - (2 * rho1**2 / 9 * Lam + ek[5] * (2 * rho1 / 9 - 4 * Lam / 3))) == 0,
        sp.expand(sp.cancel(Lam.subs(w2, w2om)) - r * w1**2 * (1 - w * om)) == 0,  # Lambda = r w1^2 (1 - w omega)
    ]
    # modulo Lambda (w2 = (w w1 - r w1^2)/(r w), w a unit): F from c_4, G0 from c_0, then a resultant
    w2s = (w * w1 - r * w1**2) / (r * w)
    Fw = sp.expand(sp.cancel(ek[4].subs(w2, w2s) * w**2 / (2 * r**2 * w1**4)))
    G0w = sp.expand(sp.cancel(ek[0].subs(w2, w2s) * w**2 / (2 * r * w1**3)))
    chain.append(sp.Poly(Fw, w).degree() == 3 and sp.Poly(G0w, w).degree() == 4)
    res = sp.factor(sp.resultant(Fw, G0w, w))
    chain.append(res == 512 * aa**2 * r**4 * w1**9)
    # certificate: cofactors in (w, omega) over Q(r, w1, aa), support c_0, c_4, c_5, c_6 (odd ones dropped for the control)
    support = (0, 4, 6) if MUT == "acc-drop-odd" else (0, 4, 5, 6)
    Kf = sp.QQ.frac_field(r, w1, aa)

    def find_cert2(polys, d):
        # a certificate with monomial denominators specializes to numbers, so a failure at a numeric point rules it out
        if not numeric_has_cert(polys, d):
            return None
        mons = [(i, j) for i in range(d + 1) for j in range(d + 1 - i)]
        unknown = [(k, m) for k in range(len(polys)) for m in mons]
        cols, allm = [], {(0, 0)}
        for (k, m) in unknown:
            dct = dict(sp.Poly(sp.expand(w ** m[0] * om ** m[1] * polys[k]), w, om).terms())
            cols.append(dct)
            allm |= set(dct)
        allm = sorted(allm)
        rows = [[Kf.from_sympy(sp.sympify(cd.get(mon, 0))) for cd in cols] + [Kf.from_sympy(sp.Integer(1 if mon == (0, 0) else 0))]
                for mon in allm]
        R_, piv = DomainMatrix(rows, (len(allm), len(unknown) + 1), Kf).rref()
        if len(unknown) in piv:
            return None
        Rl = R_.to_Matrix()
        sol = {}
        for row, pcol in enumerate(piv):
            val = Rl[row, len(unknown)]
            if val != 0:
                sol[unknown[pcol]] = sp.factor(val)
        return sol

    def numeric_has_cert(polys, d, pt={r: 2, w1: 3, aa: 5}):
        mons = [(i, j) for i in range(d + 1) for j in range(d + 1 - i)]
        unknown = [(k, m) for k in range(len(polys)) for m in mons]
        cols, allm = [], {(0, 0)}
        for (k, m) in unknown:
            dct = dict(sp.Poly(sp.expand((w ** m[0] * om ** m[1] * polys[k]).subs(pt)), w, om).terms())
            cols.append(dct)
            allm |= set(dct)
        allm = sorted(allm)
        rows = [[sp.QQ.from_sympy(sp.sympify(cd.get(mon, 0))) for cd in cols] + [sp.QQ.from_sympy(sp.Integer(1 if mon == (0, 0) else 0))]
                for mon in allm]
        _, piv = DomainMatrix(rows, (len(allm), len(unknown) + 1), sp.QQ).rref()
        return len(unknown) not in piv

    polys = [ekw[k] for k in support]
    sol, dused = None, None
    for d in range(2, 6):
        sol = find_cert2(polys, d)
        if sol is not None:
            dused = d
            break
    ok_cert = False
    if sol is not None:
        tot = sum(co * w ** m[0] * om ** m[1] * polys[k] for (k, m), co in sol.items())
        ok_cert = sp.expand(sp.cancel(sp.together(tot - 1))) == 0
        mono = all(sp.Poly(sp.denom(sp.together(co)), r, w1, aa).is_monomial for co in sol.values())
        # back in (w, w2): omega -> (w1 - r w2)/(r w1^2), the identity holds in the original variables
        omb = (w1 - r * w2) / (r * w1**2)
        tot2 = sum(co * w ** m[0] * omb ** m[1] * ek[support[k]] for (k, m), co in sol.items())
        ok_back = sp.expand(sp.cancel(sp.together(tot2 - 1))) == 0
        ok_cert = ok_cert and mono and ok_back
    check("acc_unit", ok_cert, "support c_%s: sum g_k c_k = 1 with cofactors of degree %s in (w, omega), %s nonzero coefficients, "
          "monomial denominators in r, w1, aa, verified by exact expansion also in (w, w2)%s"
          % (",".join(str(k) for k in support), dused, None if sol is None else len(sol),
             "" if sol is not None else " [NO CERTIFICATE up to degree 5]"))
    if "certs" in OPTS and sol is not None:
        print("  CERTIFICATE (accelerating bubble). e_k = c_k for k even, c_k/(4 aa r w1) for k odd, in (w, omega), omega = (w1 - r w2)/(r w1^2); "
              "sum_k cofactor(k, w^i omega^j) w^i omega^j e_k = 1, support %s" % (support,))
        print("  e_k:", [str(ekw[k]) for k in support])
        for key in sorted(sol):
            print("    (k=%d, w^%d omega^%d): %s" % (support[key[0]], key[1][0], key[1][1], sol[key]))
    print("  chain: c^5/(4 aa r w1) = Lam (3 Lam - r w1^2), c^6 = -2 Lam^2 (2 Lam - r w1^2), Lam = r w1^2 (1 - w omega), "
          "c^6 = (2 (r w1^2)^2 / 9) Lam + e_5 (2 r w1^2/9 - 4 Lam/3), F and G0 of degree 3 and 4: %s; Res_w(F, G0) = %s"
          % (chain[:4], res))

    # ---- 9. A shift b(rho) -------------------------------------------------------------------------------------------
    print("\n== 9. A shift b(rho): D = h (h^2 - Q^2); the solutions of b'' + b'/rho = eps b'^2 ==")
    tt_, xx_, y_, z_ = sp.symbols("tt xx y_ z_", real=True)
    rb = sp.Symbol("rho_b", positive=True)
    Bf = sp.Function("Bf")
    bexp = Bf(sp.sqrt(y_**2 + z_**2))  # a shift depending on rho only: no t or x dependence
    nb = lambda g: sp.diff(g, tt_) - bexp * sp.diff(g, xx_)
    bx_ = sp.diff(bexp, xx_)
    brho_ = (y_ * sp.diff(bexp, y_) + z_ * sp.diff(bexp, z_)) / sp.sqrt(y_**2 + z_**2)
    Pb = nb(bx_) - bx_**2
    Zb = sp.diff(brho_, xx_) / 2
    Yb = -nb(brho_) / 2 + bx_ * brho_
    pzy_zero = (Pb, Zb, Yb) == (0, 0, 0)
    # Delta_perp b(rho) = b'' + b'/rho, tested on three explicit functions with y, z positive
    yp, zp = sp.symbols("yp zp", positive=True)
    rp = sp.sqrt(yp**2 + zp**2)
    lapz = 0
    for fx in (lambda q_: q_**3, sp.sin, lambda q_: sp.log(1 + q_**2)):
        bb_ = fx(rp)
        lhs = sp.diff(bb_, yp, 2) + sp.diff(bb_, zp, 2)
        qq_ = sp.Symbol("qq", positive=True)
        rhs = (sp.diff(fx(qq_), qq_, 2) + sp.diff(fx(qq_), qq_) / qq_).subs(qq_, rp)
        lapz += sp.simplify(lhs - rhs) != 0
    ok_b = []
    for eps in (1, -1):
        for C0 in (0, 3):
            sol_b = -eps * sp.log(C0 - eps * sp.log(rb))
            ode = sp.simplify(sol_b.diff(rb, 2) + sol_b.diff(rb) / rb - eps * sol_b.diff(rb) ** 2)
            Qv = -(sol_b.diff(rb, 2) + sol_b.diff(rb) / rb) / 2
            hv = sol_b.diff(rb) ** 2 / 2
            hv_exp = hv * (2 if MUT == "brho-h2" else 1)
            Dv = sp.simplify(D_of(0, Qv, 0, 0, hv))
            Dexp = sp.simplify(hv_exp * (hv_exp**2 - Qv**2))
            ok_b.append(ode == 0 and Dv == 0 and sp.simplify(Dv - Dexp) == 0)
    check("brho_D", pzy_zero and lapz == 0 and all(ok_b),
          "for b(rho): P = Z = Y = 0 (%s) and Delta_perp b = b'' + b'/rho (3 test functions, %s differ), so D = h (h^2 - Q^2); b = -eps ln(C - eps ln rho) "
          "solves b'' + b'/rho = eps b'^2 and has D = 0 for eps = +-1 and C = 0, 3: %s" % (pzy_zero, lapz, ok_b))

    fails = [n for n, ok in RESULTS if not ok]
    print("\nchecks: %d passed, %d failed %s. %.0f s" % (len(RESULTS) - len(fails), len(fails), fails, time.time() - t0))
    if MUT:
        target = CONTROLS[MUT][1]
        fired = all(n in fails for n in target)
        print("control %s: targets %s %s" % (MUT, target, "FIRED" if fired else "SILENT"))
        sys.exit(1)
    status = 1 if fails else 0
    if "no-controls" not in OPTS:
        print("\n== Controls: each mutation must make its target checks fail ==")
        for name, (what, target) in CONTROLS.items():
            p = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True)
            failed = [ln.split("]")[1].split(":")[0].strip() for ln in p.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = p.returncode == 1 and all(n in failed for n in target)
            status |= 0 if fired else 1
            print("[%s] --mutate=%s (%s): exit %d, failed %s" % ("FIRED" if fired else "SILENT", name, what,
                                                                p.returncode, failed), flush=True)
        print("controls: %s" % ("all fired" if status == 0 else "SOME SILENT OR CHECK FAILED"))
    sys.exit(status)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Independent review of #65, part D, for issue #79: the determinant D of the frame Ricci block R_3.

For ds^2 = -dt^2 + (dx + b dt)^2 + dy^2 + dz^2 with b = b(t, x, rho), #65 prints, with n = d_t - b d_x,
    P = n(b_x) - b_x^2,  Q = -Delta_perp b / 2,  Z = d_x b_rho / 2,  Y = -n(b_rho) / 2 + b_x b_rho,  h = b_rho^2 / 2,
    D = h (h^2 - P^2 - Q^2 + Y^2 + Z^2) - P (Y^2 - Z^2) + 2 Q Y Z = det R_3,
R_3 the block of the Ricci tensor on (n, d_x, e_rho). This script checks that against two routes that share no code with
phantom_source.py or two_scalars.py:
  - the Cartan routine of review_cartan.py (#59), imported: Jets, axisymmetric_shift, cartan_ricci;
  - a coordinate-Christoffel routine written here (coord_ricci_at), at exact points.

    python3 notes/warp/review_two_scalars_cartan.py                  checks, then every control
    python3 notes/warp/review_two_scalars_cartan.py --no-controls    checks only
    python3 notes/warp/review_two_scalars_cartan.py --mutate=NAME    one control (its target checks must fail, exit 1)

Exit status is 1 if a check fails or a control does not fire.
"""
import importlib.util
import os
import subprocess
import sys
import time

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate", "")
CONTROLS = {
    "d-qyz-sign": ("printed D with -2QYZ", ["D_identity", "D_cartan_jets", "D_coordinate", "D_slice_bz0"]),
    "d-drop-z2": ("printed D without the Z^2 terms", ["D_identity", "D_cartan_jets", "D_coordinate", "D_slice_bz0"]),
    "blocks-q-sign": ("expected R_01 = +Delta_perp b / 2", ["blocks_jets", "D_cartan_jets", "D_slice_bz0", "D_slice_offslice"]),
    "cartan-lapse": ("imported routine with lapse 1 + b", ["blocks_jets", "D_cartan_jets"]),
    "cartan-sign": ("imported routine, sign of omega^omega flipped", ["cartan_planted", "blocks_jets", "D_cartan_jets"]),
    "nonaxi": ("shift a function of y^2 + 2 z^2", ["blocks_jets", "planted_D"]),  # README A1: original target list [blocks_jets, D_cartan_jets]; D_cartan_jets stayed silent
    "coord-gamma": ("my Christoffel symbol with a wrong index order", ["coord_planted", "D_coordinate", "planted_D", "D_slice_bz0", "D_slice_offslice"]),
    "planted-h": ("planted expectation with h doubled", ["planted_D"]),
}
if MUT and MUT not in CONTROLS:
    sys.exit(f"unknown mutation {MUT}")

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def load_review_cartan():
    """Import review_cartan.py without letting it read this script's command line."""
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
if MUT == "cartan-sign":
    rc.MUT = "cartan-sign"  # module global read inside rc.cartan_ricci


def D_printed(P, Q, Z, Y, h):
    """#65's D, as printed, with the two mutations of the claim."""
    sq = -1 if MUT == "d-qyz-sign" else 1
    z2 = 0 if MUT == "d-drop-z2" else 1
    return h * (h**2 - P**2 - Q**2 + Y**2 + z2 * Z**2) - P * (Y**2 - z2 * Z**2) + sq * 2 * Q * Y * Z


def Q_sign():
    return +1 if MUT == "blocks-q-sign" else -1  # Q = -Delta_perp b / 2


# ---------------------------------------------------------------------------------------------------------------
# My own coordinate-Christoffel Ricci tensor at an exact point. No inverse is taken symbolically: the jets of g up to
# second order are evaluated at the point, then Gamma, d Gamma and R_bc = d_a G^a_bc - d_c G^a_ba + G^a_ad G^d_bc
# - G^a_cd G^d_ba are computed with exact (or 80-digit) numbers.
# ---------------------------------------------------------------------------------------------------------------
def coord_ricci_at(g, X, pt, prec=None):
    n = 4

    def ev(e):
        e = sp.sympify(e).subs(pt)
        return sp.N(e, prec) if prec else e

    g0 = sp.Matrix(n, n, lambda a, b: ev(g[a, b]))
    dg = [[[ev(sp.diff(g[a, b], X[d])) for b in range(n)] for a in range(n)] for d in range(n)]
    ddg = [[[[ev(sp.diff(g[a, b], X[d], X[e])) for b in range(n)] for a in range(n)] for e in range(n)]
           for d in range(n)]
    gi = g0.inv()
    if MUT == "coord-gamma":
        def low(a, b, c):  # wrong: not Gamma_{a,bc}
            return (dg[b][a][c] + dg[a][b][c] - dg[c][a][b]) / 2

        def dlow(d, a, b, c):
            return (ddg[d][b][a][c] + ddg[d][a][b][c] - ddg[d][c][a][b]) / 2
    else:
        def low(a, b, c):  # Gamma_{a,bc}
            return (dg[b][a][c] + dg[c][a][b] - dg[a][b][c]) / 2

        def dlow(d, a, b, c):  # d_d Gamma_{a,bc}
            return (ddg[d][b][a][c] + ddg[d][c][a][b] - ddg[d][a][b][c]) / 2
    Gam = [[[sum(gi[a, e] * low(e, b, c) for e in range(n)) for c in range(n)] for b in range(n)] for a in range(n)]
    dgi = [[[-sum(gi[a, f] * dg[d][f][h] * gi[h, e] for f in range(n) for h in range(n)) for e in range(n)]
            for a in range(n)] for d in range(n)]  # dgi[d][a][e] = d_d g^{ae}
    dGam = [[[[sum(dgi[d][a][e] * low(e, b, c) + gi[a, e] * dlow(d, e, b, c) for e in range(n))
               for c in range(n)] for b in range(n)] for a in range(n)] for d in range(n)]  # [d][a][b][c]
    Ric = sp.zeros(n, n)
    for b in range(n):
        for c in range(n):
            Ric[b, c] = (sum(dGam[a][a][b][c] - dGam[c][a][b][a] for a in range(n))
                         + sum(Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][b][a]
                               for a in range(n) for d in range(n)))
    return g0, Ric


def shift_metric(b):
    return sp.Matrix([[-1 + b**2, b, 0, 0], [b, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])


def printed_from_b(bexpr, X, pt):
    """P, Q, Z, Y, h of #65 at a point, from sympy derivatives of the explicit b(t, x, y, z) (rho = sqrt(y^2 + z^2))."""
    t, x, y, z = X
    d = sp.diff
    nb = lambda e: d(e, t) - bexpr * d(e, x)  # n = d_t - b d_x
    brho = (y * d(bexpr, y) + z * d(bexpr, z)) / sp.sqrt(y**2 + z**2)
    bx = d(bexpr, x)
    P = nb(bx) - bx**2
    Q = Q_sign() * (d(bexpr, y, y) + d(bexpr, z, z)) / 2
    Z = d(brho, x) / 2
    Y = -nb(brho) / 2 + bx * brho
    h = brho**2 / 2
    return [sp.simplify(sp.sympify(e).subs(pt)) for e in (P, Q, Z, Y, h)]


def frame_R3(g0, Ric, bval, y, z):
    """R(V_a, V_b) for V = (n, d_x, e_rho) at the point; also returns g(V, V) to confirm orthonormality."""
    rho = sp.sqrt(y**2 + z**2)
    V = [sp.Matrix([1, -bval, 0, 0]), sp.Matrix([0, 1, 0, 0]), sp.Matrix([0, 0, y / rho, z / rho])]
    R3 = sp.Matrix(3, 3, lambda i, j: (V[i].T * Ric * V[j])[0])
    G3 = sp.Matrix(3, 3, lambda i, j: (V[i].T * g0 * V[j])[0])
    return R3, G3


def main():
    t0 = time.time()
    print("review_two_scalars_cartan.py, issue #79. sympy %s, python %s. mutation: %s"
          % (sp.__version__, sys.version.split()[0], MUT or "none"))
    print("imports from review_cartan.py: Jets, axisymmetric_shift, cartan_ricci (module MUT = %r)" % (rc.MUT,))

    # ---- 1. Planted truth for both Ricci routines ----------------------------------------------------------------
    print("\n== 1. Planted truth: de Sitter (flat slicing) and Ellis-Bronnikov, through both routines ==")
    tt, r, th, ph = sp.symbols("t r theta phi", positive=True)
    Hs, asym = sp.symbols("H a", positive=True)
    plain = rc.Jets((tt, r, th, ph))
    Ec = sp.diag(1, sp.exp(Hs * tt), sp.exp(Hs * tt), sp.exp(Hs * tt))
    Ric, bad, _ = rc.cartan_ricci(Ec, plain.total)
    ok_ds = (sp.simplify(Ric - 3 * Hs**2 * rc.ETA) == sp.zeros(4)) and not bad
    rr = sp.sqrt(r**2 + asym**2)
    Ec = sp.diag(1, 1, rr, rr * sp.sin(th))
    Ric, bad, _ = rc.cartan_ricci(Ec, plain.total, simp=sp.simplify)
    want = sp.diag(0, -2 * asym**2 / (r**2 + asym**2) ** 2, 0, 0)
    ok_eb = (Ric - want).applyfunc(sp.simplify) == sp.zeros(4) and not bad
    check("cartan_planted", ok_ds and ok_eb, "imported Cartan routine: de Sitter R = 3H^2 eta %s, "
          "Ellis-Bronnikov R_ll = -2a^2/(l^2+a^2)^2 %s" % (ok_ds, ok_eb))

    Hq, aq = sp.Rational(3, 7), sp.Rational(2, 5)
    X = sp.symbols("T L U V", real=True)  # (t, l, theta, phi)
    g_ds = sp.diag(-1, sp.exp(2 * Hq * X[0]), sp.exp(2 * Hq * X[0]), sp.exp(2 * Hq * X[0]))
    g0, Rds = coord_ricci_at(g_ds, X, {X[0]: 0, X[1]: sp.Rational(1, 3), X[2]: sp.Rational(1, 2), X[3]: 0})
    ok_ds2 = sp.simplify(Rds - 3 * Hq**2 * g0) == sp.zeros(4)
    l0 = sp.Rational(7, 5)
    g_eb = sp.diag(-1, 1, X[1] ** 2 + aq**2, (X[1] ** 2 + aq**2) * sp.sin(X[2]) ** 2)
    g0, Reb = coord_ricci_at(g_eb, X, {X[0]: 0, X[1]: l0, X[2]: sp.pi / 2, X[3]: 0})
    want = sp.diag(0, -2 * aq**2 / (l0**2 + aq**2) ** 2, 0, 0)
    ok_eb2 = sp.simplify(Reb - want) == sp.zeros(4)
    check("coord_planted", ok_ds2 and ok_eb2, "my coordinate routine: de Sitter R = 3H^2 g %s, Ellis-Bronnikov "
          "R_ll = %s (expected %s), others 0: %s" % (ok_ds2, Reb[1, 1], want[1, 1], ok_eb2))

    # ---- 2. Cartan Ricci of b = G(t, x, q), q = y^2 + z^2, against the printed blocks and D ------------------------
    print("\n== 2. Cartan routine, b = G(t, x, y^2 + z^2): the blocks of R_3 and D, exact in the jets of G ==")
    t, x, y, z = sp.symbols("t x y z", real=True)
    J = rc.Jets((t, x, y, z))
    b, G = rc.axisymmetric_shift(J, weight_z=2 if MUT == "nonaxi" else 1)
    lapse = 1 + b if MUT == "cartan-lapse" else 1
    E = sp.Matrix([[lapse, 0, 0, 0], [b, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    Ric, bad, _ = rc.cartan_ricci(E, J.total)
    q = y**2 + z**2
    nn = lambda e: J.total(e, 0) - b * J.total(e, 1)  # my n = d_t - b d_x
    by, bz, bx = J.total(b, 2), J.total(b, 3), J.total(b, 1)
    beta = y * by + z * bz  # rho b_rho
    P = nn(bx) - bx**2
    Q = Q_sign() * (J.total(by, 2) + J.total(bz, 3)) / 2
    Zr = J.total(beta, 1) / 2  # rho Z
    Yr = -nn(beta) / 2 + bx * beta  # rho Y
    h = sp.cancel(sp.expand(beta**2) / (2 * q))
    R = Ric
    diffs = {
        "R00-(P-h)": R[0, 0] - (P - h),
        "R11-(-P-h)": R[1, 1] - (-P - h),
        "R01-Q": R[0, 1] - Q,
        "rho R0rho-rho Z": y * R[0, 2] + z * R[0, 3] - Zr,
        "rho R1rho-rho Y": y * R[1, 2] + z * R[1, 3] - Yr,
        "rho^2 Rrr-rho^2 h": y**2 * R[2, 2] + 2 * y * z * R[2, 3] + z**2 * R[3, 3] - q * h,
        "rho R0phi": -z * R[0, 2] + y * R[0, 3],
        "rho R1phi": -z * R[1, 2] + y * R[1, 3],
        "rho^2 Rrphi": -y * z * R[2, 2] + (y**2 - z**2) * R[2, 3] + y * z * R[3, 3],
        "rho^2 Rphiphi": z**2 * R[2, 2] - 2 * y * z * R[2, 3] + y**2 * R[3, 3],
    }
    bad_keys = [k for k, v in diffs.items() if sp.expand(v) != 0]
    check("blocks_jets", not bad_keys and h.is_polynomial(y, z) and not bad,
          "R_00 = P-h, R_11 = -P-h, R_01 = Q, R_0rho = Z, R_1rho = Y, R_rho,rho = h, R_a,phi = 0; "
          "structure-equation violations %d; nonzero differences %s" % (len(bad), bad_keys))

    rho = sp.Symbol("rho", positive=True)
    sl = {y: rho, z: 0}
    Pm, Qm, hm = (sp.expand(e.subs(sl)) for e in (P, Q, h))
    Zm, Ym = sp.cancel(sp.expand(Zr.subs(sl)) / rho), sp.cancel(sp.expand(Yr.subs(sl)) / rho)
    R3 = sp.Matrix([[R[0, 0], R[0, 1], R[0, 2]], [R[0, 1], R[1, 1], R[1, 2]], [R[0, 2], R[1, 2], R[2, 2]]]).subs(sl)
    detR3 = sp.expand(R3.det())
    Dp = sp.expand(D_printed(Pm, Qm, Zm, Ym, hm))
    dd = sp.expand(detR3 - Dp)
    check("D_cartan_jets", dd == 0, "det R_3 (Cartan, in the frame n, d_x, e_rho at y = rho, z = 0) minus D(printed), "
          "exact in the jets of G: %s (%d terms in D)" % (dd if dd != 0 else 0, len(Dp.args)))
    print("  D = " + str(sp.factor(Dp))[:400])

    # ---- 3. The algebra of D, and the axis ----------------------------------------------------------------------
    print("\n== 3. det of the printed blocks, and the axis ==")
    p_, q_, z_, y_, h_ = sp.symbols("P Q Z Y h")
    Mb = sp.Matrix([[p_ - h_, q_, z_], [q_, -p_ - h_, y_], [z_, y_, h_]])
    di = sp.expand(Mb.det() - D_printed(p_, q_, z_, y_, h_))
    check("D_identity", di == 0, "det [[P-h, Q, Z], [Q, -P-h, Y], [Z, Y, h]] - D(printed) = %s" % (di if di != 0 else 0))
    ax = {y: 0, z: 0}
    rows = [sp.expand(R[2, c].subs(ax)) for c in range(4)] + [sp.expand(R[3, c].subs(ax)) for c in range(4)]
    R3a = sp.Matrix([[R[0, 0], R[0, 1], R[0, 2]], [R[0, 1], R[1, 1], R[1, 2]], [R[0, 2], R[1, 2], R[2, 2]]]).subs(ax)
    Da = sp.expand(D_printed(Pm, Qm, Zm, Ym, hm).subs(rho, 0))
    i2_axis = sp.expand((R[0, 0] + R[1, 1]).subs(ax))
    check("axis_D", all(r_ == 0 for r_ in rows) and sp.expand(R3a.det()) == 0 and Da == 0 and i2_axis == 0,
          "on y = z = 0: R(E_2, .), R(E_3, .) = 0, R_00 + R_11 = %s, det R_3 = 0, D(printed)(rho = 0) = %s" % (i2_axis, Da))

    # ---- 4. My coordinate route at exact points, three explicit shifts ------------------------------------------
    print("\n== 4. Coordinate-Christoffel Ricci at exact points: the entries of R_3 and D against the printed formulas ==")
    X4 = (t, x, y, z)
    qq = y**2 + z**2
    shifts = {
        "b1 = (1 + t x)(rho^2 + x)/(2 + t^2 + rho^2)": (1 + t * x) * (qq + x) / (2 + t**2 + qq),
        "b2 = (x^2 - t) rho^2/(1 + rho^2) + t/3": (x**2 - t) * qq / (1 + qq) + t / 3,
        "b3 = rho^2 (3 + t + x^2)/(1 + x^2 + t^2 rho^2)": qq * (3 + t + x**2) / (1 + x**2 + t**2 * qq),
    }
    pts = [(sp.Rational(1, 3), -sp.Rational(2, 5), sp.Rational(3, 7), sp.Rational(4, 7)),
           (-sp.Rational(1, 2), sp.Rational(3, 4), sp.Rational(5, 13), sp.Rational(12, 13)),
           (sp.Rational(2, 3), sp.Rational(1, 5), -sp.Rational(8, 17), sp.Rational(15, 17))]
    worst, nonz, count = 0, 0, 0
    for name, bexpr in shifts.items():
        gm = shift_metric(bexpr)
        for pt in pts:
            ptd = dict(zip(X4, pt))
            g0, Rc = coord_ricci_at(gm, X4, ptd)
            R3c, G3 = frame_R3(g0, Rc, bexpr.subs(ptd), pt[2], pt[3])
            Pp, Qp, Zp, Yp, hp = printed_from_b(bexpr, X4, ptd)
            Rexp = sp.Matrix([[Pp - hp, Qp, Zp], [Qp, -Pp - hp, Yp], [Zp, Yp, hp]])
            e_ent = sp.simplify(R3c - Rexp)
            e_det = sp.simplify(R3c.det() - D_printed(Pp, Qp, Zp, Yp, hp))
            e_og = sp.simplify(G3 - sp.diag(-1, 1, 1))
            count += 1
            if e_ent != sp.zeros(3) or e_det != 0 or e_og != sp.zeros(3):
                worst += 1
            nonz += e_det == 0 and sp.simplify(R3c.det()) != 0
    check("D_coordinate", worst == 0 and nonz == 9, "%d of %d (shift, point) pairs agree exactly in every entry of R_3 and "
          "in D; D != 0 at %d of them; frame orthonormal" % (count - worst, count, nonz))

    # ---- 4b. Non-axisymmetric shifts (README A4: written after the claim "D(printed) = det R(n, d_x, d_y) for every b" failed) ----
    print("\n== 4b. Non-axisymmetric b(t, x, y, z): the block R(n, d_x, d_y) against [[P-h, Q, Z], [Q, -P-h, Y], [Z, Y, h]], h = b_y^2/2 ==")

    def slice_block(bexpr, pt):
        """my coordinate Ricci in the frame (n, d_x, d_y) at the point, the printed block from derivatives of b, b_z, and R_zz."""
        ptd = dict(zip(X4, pt))
        g0, Rc = coord_ricci_at(shift_metric(bexpr), X4, ptd)
        bv = bexpr.subs(ptd)
        V = [sp.Matrix([1, -bv, 0, 0]), sp.Matrix([0, 1, 0, 0]), sp.Matrix([0, 0, 1, 0])]
        R3c = sp.Matrix(3, 3, lambda i, j: (V[i].T * Rc * V[j])[0])
        nb = lambda e: sp.diff(e, t) - bexpr * sp.diff(e, x)
        by_, bx_, bz_ = sp.diff(bexpr, y), sp.diff(bexpr, x), sp.diff(bexpr, z)
        Pp = nb(bx_) - bx_**2
        Qp = Q_sign() * (sp.diff(bexpr, y, y) + sp.diff(bexpr, z, z)) / 2
        Zp = sp.diff(by_, x) / 2
        Yp = -nb(by_) / 2 + bx_ * by_
        hp = by_**2 / 2
        Pp, Qp, Zp, Yp, hp, bzv = (sp.simplify(e.subs(ptd)) for e in (Pp, Qp, Zp, Yp, hp, bz_))
        Rexp = sp.Matrix([[Pp - hp, Qp, Zp], [Qp, -Pp - hp, Yp], [Zp, Yp, hp]])
        return R3c, Rexp, (Pp, Qp, Zp, Yp, hp), bzv

    even_z = {  # not functions of rho, even in z
        "b1 = (y + 2 z^2 + t x)/(1 + y^2 + x^2)": (y + 2 * z**2 + t * x) / (1 + y**2 + x**2),
        "b2 = (1 + t)(x^2 + y z^2)/(2 + y^2 + 3 z^2)": (1 + t) * (x**2 + y * z**2) / (2 + y**2 + 3 * z**2),
        "b3 = t y^3 + x z^2/(1 + z^2) + x y^2": t * y**3 + x * z**2 / (1 + z**2) + x * y**2,
    }
    pts0 = [(sp.Rational(1, 3), -sp.Rational(2, 5), sp.Rational(3, 7), 0), (-sp.Rational(1, 2), sp.Rational(3, 4), sp.Rational(5, 13), 0),
            (sp.Rational(2, 3), sp.Rational(1, 5), -sp.Rational(8, 17), 0)]
    bad0 = 0
    for name, bexpr in even_z.items():
        for pt in pts0:
            R3c, Rexp, (Pp, Qp, Zp, Yp, hp), bzv = slice_block(bexpr, pt)
            ok = sp.simplify(R3c - Rexp) == sp.zeros(3) and sp.simplify(R3c.det() - D_printed(Pp, Qp, Zp, Yp, hp)) == 0 and bzv == 0
            bad0 += not ok
    bzj = sp.expand(J.total(b, 3).subs(z, 0))  # b = G(t, x, y^2 + z^2) jets, imported axisymmetric_shift
    check("D_slice_bz0", bad0 == 0 and bzj == 0,
          "%d of 9 (shift even in z and not a function of rho, point with z = 0): all six entries of R(n, d_x, d_y) equal the printed block "
          "and det = D(printed); for the axisymmetric b = G(t, x, y^2 + z^2) of section 2, b_z on z = 0 is %s exactly in the jets"
          % (9 - bad0, bzj))

    gen = {
        "b1 = (y + 2 z^2 + t x)/(1 + y^2 + x^2)": (y + 2 * z**2 + t * x) / (1 + y**2 + x**2),
        "b2 = (1 + t)(x^2 + y z)/(2 + y^2 + 3 z^2)": (1 + t) * (x**2 + y * z) / (2 + y**2 + 3 * z**2),
        "b3 = t y^3 + x z/(1 + z^2) + x y^2": t * y**3 + x * z / (1 + z**2) + x * y**2,
    }
    npts = [(sp.Rational(1, 3), -sp.Rational(2, 5), sp.Rational(3, 7), sp.Rational(4, 7)),
            (-sp.Rational(1, 2), sp.Rational(3, 4), sp.Rational(5, 13), sp.Rational(12, 13)),
            (sp.Rational(2, 3), sp.Rational(1, 5), -sp.Rational(8, 17), sp.Rational(15, 17))]
    bad1, nbz, ndet = 0, 0, 0
    for name, bexpr in gen.items():
        for pt in npts:
            R3c, Rexp, (Pp, Qp, Zp, Yp, hp), bzv = slice_block(bexpr, pt)
            diff = (R3c - Rexp).applyfunc(sp.simplify)
            half = bzv**2 / 2
            ok = (diff[0, 0] == -half and diff[1, 1] == -half and diff[0, 1] == 0 and diff[0, 2] == 0 and diff[1, 2] == 0
                  and diff[2, 2] == 0)
            bad1 += not ok
            nbz += bzv != 0
            ndet += sp.simplify(R3c.det() - D_printed(Pp, Qp, Zp, Yp, hp)) != 0
    check("D_slice_offslice", bad1 == 0 and nbz == 9 and ndet == 9,
          "%d of 9 (non-axisymmetric shift, point with z != 0): only R_00 and R_11 differ from the block, each by exactly -b_z^2/2 "
          "(R_00 = P - h - b_z^2/2, R_11 = -P - h - b_z^2/2), the other four entries agree; b_z != 0 at %d, "
          "det R(n, d_x, d_y) != D(printed) at %d" % (9 - bad1, nbz, ndet))

    # ---- 5. Planted values of D: b = k rho^2 and the two b(rho) solutions ----------------------------------------
    print("\n== 5. Planted D: b = k rho^2 (D = h (h^2 - 4k^2)) and b = -eps ln|C - eps ln rho| (D = 0), both routes ==")
    k = sp.Rational(4, 9)
    hexp = 2 * k**2 * rho**2 * (2 if MUT == "planted-h" else 1)
    Dexp = hexp * (hexp**2 - 4 * k**2)
    ok = []
    # jet route: b = G(q) with G_000 = k q, G_001 = k, higher zero
    subs_k = {s: 0 for s in R3.free_symbols if s.name.startswith("G_")}
    subs_k[G((0, 0, 0))] = k * rho**2
    subs_k[G((0, 0, 1))] = k
    jet_val = sp.simplify(R3.subs(subs_k).det())
    ok.append(sp.simplify(jet_val - Dexp) == 0)
    # coordinate route at rho = 5/7 and at rho^2 = 1/k = 9/4 (D = 0)
    bk = k * (y**2 + z**2)
    gm = shift_metric(bk)
    for yy in (sp.Rational(5, 7), sp.Rational(3, 2)):
        ptd = {t: sp.Rational(1, 3), x: -sp.Rational(2, 5), y: yy, z: 0}
        g0, Rc = coord_ricci_at(gm, X4, ptd)
        R3c, _ = frame_R3(g0, Rc, bk.subs(ptd), yy, 0)
        ok.append(sp.simplify(R3c.det() - Dexp.subs(rho, yy)) == 0)
    # log solutions: coordinate route at 80 digits, jet route at 70 digits
    Lq = sp.Symbol("q", positive=True)
    for eps, C0, rho0 in ((1, 0, sp.Rational(3, 5)), (-1, 0, sp.Rational(5, 2)), (1, 2, sp.Rational(5, 2)),
                          (-1, 1, sp.Rational(1, 2))):
        sg = 1 if sp.N(C0 - eps * sp.log(rho0)) > 0 else -1  # ln|u| = ln(sg u) where u has the sign sg
        bsol = -eps * sp.log(sg * (C0 - eps * sp.log(sp.sqrt(y**2 + z**2))))
        gm = shift_metric(bsol)
        ptd = {t: 0, x: 0, y: rho0, z: 0}
        g0, Rc = coord_ricci_at(gm, X4, ptd, prec=80)
        R3c, _ = frame_R3(g0, Rc, bsol.subs(ptd).evalf(80), rho0, 0)
        dval = sp.N(R3c.det(), 80)
        scale = max(abs(e) for e in R3c)
        ok.append(abs(dval) < sp.Float(10) ** -60 * max(1, scale**3) and scale > sp.Float(10) ** -5)
        Gq = -eps * sp.log(sg * (C0 - eps * sp.log(Lq) / 2))
        sj = {s: 0 for s in R3.free_symbols if s.name.startswith("G_")}
        sj[G((0, 0, 0))] = sp.N(Gq.subs(Lq, rho0**2), 70)
        sj[G((0, 0, 1))] = sp.N(sp.diff(Gq, Lq).subs(Lq, rho0**2), 70)
        sj[G((0, 0, 2))] = sp.N(sp.diff(Gq, Lq, 2).subs(Lq, rho0**2), 70)
        dj = R3.subs(sj).subs(rho, rho0).det()
        ok.append(abs(sp.N(dj, 70)) < sp.Float(10) ** -50)
    check("planted_D", all(ok), "b = k rho^2: both routes give D = %s at three points; four b(rho) solutions have D = 0 "
          "(coordinate route 80 digits, jet route 60 digits): %s" % (sp.factor(Dexp), ok))

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

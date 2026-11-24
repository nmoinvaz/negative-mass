#!/usr/bin/env python3
"""Does a sheared warp shift admit the scalar form R_ab = lam g_ab + s u_a u_b? (issue #81)

#70 showed that (I1) and (I2) hold in the frame adapted to grad_perp(sigma) for every shift b(t, x, sigma(y, z)), so the
lemma of phantom_source.md applies there. Here the shift is

    ds^2 = -dt^2 + (dx + b dt)^2 + dy^2 + dz^2,   b = -v f(q),   q = X^2 + y^2 + 2 z^2 + X y / 2,   X = x - v t,

for which grad_perp(d_x b) is not along grad_perp(b) wherever v f'(q) z != 0. It reuses phantom_source.py's `scalar_form`,
`ricci_at_point`, `shift_tetrad_ricci` and `JetD`, and nonaxisymmetric.py's `q_rule`, `route_T`, `lemma_residues`, `G`,
`ETA` and jets of profiles.

    python3 notes/warp/sheared.py                       all checks, then every control (exit 1 on any failure)
    python3 notes/warp/sheared.py --mutate=NAME         one control (its target checks must fail, exit 1)
    python3 notes/warp/sheared.py --no-controls         the checks only
    python3 notes/warp/sheared.py --logdir=DIR          also write each control's stdout to DIR/NAME.log (numbered)

Sections. D: Ricci tensor of the shift on the jets of f (route D), its closed form, and the identities. X: route D
against routes E and T at random points. P: the exact sample points. L: planted scalar forms at the same points. R: the
by-hand reduction of sheared.md on planted matrices. Z: a point on the z axis where the sheared Ricci tensor is a scalar
form. A `--mutate` run executes the sections its target checks need.

Frame. E_0 = n = d_t - b d_x, E_1 = d_x, E_2 = d_y, E_3 = d_z. For a transverse vector u = (u_y, u_z), R(u, .) is
u_y R(E_2, .) + u_z R(E_3, .). w = grad_perp q = (q_y, q_z), wp = (-w_1, w_0), so (E_0, E_1, w, wp) is orthogonal.
"""
import os
import random
import subprocess
import sys
import time

import mpmath as mp
import sympy as sp

import nonaxisymmetric as na
import phantom_source as ps

X, y, z, v = ps.X, ps.y, ps.z, ps.v
TX = ps.TX
G = na.G  # G[k] = f^(k)(q)
ETA = na.ETA
QS = X**2 + y**2 + 2 * z**2 + X * y / 2  # the sheared q, used for every expected value

# name: (target checks that must fail, what the mutation does)
MUTATIONS = {
    "static": (["xcheck_E", "xcheck_T", "points_xcheck"], "drops every t-derivative in route D's jets"),
    "jet-g2": (["points_xcheck"], "route E at the sample points uses f'' + 1/7"),
    "lapse": (["xcheck_T", "closed_form", "shear_I1", "shear_I2", "reduction_identity"],
              "gives the slicing the lapse 1 + b in routes D and E"),
    "unsheared": (["shear_cross", "shear_I1", "closed_form", "circular_closed", "reduction_identity",
                   "points_residues"], "drops X y / 2 from the measured q, the expected values keep it"),
    "circular": (["shear_I1", "shear_I2", "points_residues"],
                 "tests (I1), (I2) in the circular frame (e_rho, e_varphi), not the frame adapted to grad_perp q"),
    "circ-swap": (["circular_closed"], "builds the circular frame from (z, y), not (y, z)"),
    "form-drop": (["closed_form"], "drops -grad_perp(d_t b - b d_x b) / 2 from the closed form of R(E_1, E_a)"),
    "red-coef": (["reduction_identity"], "doubles v^2 in the expected value of the reduction identity"),
    "planted-sign": (["planted_pass"], "plants lam eta + 8 pi u u, the wrong sign for eps = -1"),
    "planted-null": (["planted_deform"], "deforms the planted form by zero"),
    "planted-flat": (["planted_premise"], "plants R = 0, which satisfies the lemma's premises"),
    "broken-test": (["points_outcome", "lemma_agree", "planted_sign", "planted_deform", "reduction_agree"],
                    "replaces scalar_form by a test that always returns 'pass'"),
    "planted-in-data": (["points_outcome", "points_residues"],
                        "replaces the Ricci tensor at the sample points by a planted scalar form"),
    "red-no-c1": (["reduction_deform"], "the reduction ignores condition c1, R(w, E_0) = 0"),
    "red-no-c2": (["reduction_deform"], "the reduction ignores condition c2, R(w, E_1) = 0"),
    "red-no-c3": (["reduction_deform"], "the reduction ignores condition c3, on R_00"),
    "red-no-c4": (["reduction_deform"], "the reduction ignores condition c4, on R_11"),
    "red-no-c5": (["reduction_deform"], "the reduction ignores condition c5, on R_01"),
    "red-sign": (["reduction_pass"], "the reduction returns 'pass' for both signs"),
    "zaxis-off": (["zaxis_pass"], "adds 1/10 to f'' at the z-axis point of the pass side"),
}

CHECK_SECTION = {
    "xcheck_E": "X", "xcheck_T": "X",
    "closed_form": "D", "shear_cross": "D", "shear_I1": "D", "shear_I2": "D", "circular_closed": "D",
    "reduction_identity": "D",
    "points_xcheck": "P", "points_residues": "P", "points_outcome": "P", "lemma_agree": "P", "reduction_agree": "P",
    "planted_pass": "L", "planted_sign": "L", "planted_deform": "L", "planted_premise": "L",
    "reduction_pass": "R", "reduction_deform": "R",
    "zaxis_pass": "Z",
}

RESULTS = []


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def rat(a, b=1):
    return sp.Rational(a, b)


def fmt(r):
    return str(r) if sp.Rational(r).q < 10**6 else f"{float(r):.6g}"


def qexpr(mut):
    return X**2 + y**2 + 2 * z**2 if mut == "unsheared" else QS


# ---------------------------------------------------------------------------------------------------------------
# Routes D, E, T
# ---------------------------------------------------------------------------------------------------------------
def route_D(q, mut):
    """Frame Ricci tensor on the jets G_k = f^(k)(q), through phantom_source's `ricci`."""
    return ps.shift_tetrad_ricci(-v * G[0], na.q_rule(q, mut == "static"), mut == "lapse")


def explicit_b(vv, jets, q0, q):
    """b(t, x, y, z) = -v f(q), f(q) = G0 + G1 (q - q0) + G2 (q - q0)^2 / 2, so f has the jets G at q0."""
    qq = q.subs(X, TX[1] - vv * TX[0])
    return -vv * (jets[0] + jets[1] * (qq - q0) + jets[2] * (qq - q0)**2 / 2)


# ---------------------------------------------------------------------------------------------------------------
# Closed forms, from the 3+1 equations by hand (sheared.md), in the derivatives of b
# ---------------------------------------------------------------------------------------------------------------
def closed_form_tensor(mut):
    """The frame Ricci tensor of any b(t, x, y, z) with unit lapse and flat slices, evaluated on b = -v f(QS), as a
    4 x 4 matrix of expressions in X, y, z, v and the jets."""
    D = ps.JetD(na.q_rule(QS, False))
    b = -v * G[0]
    bt, bx, by, bz = (D(b, a) for a in range(4))
    bxx, bxt, bxy, bxz = D(bx, 1), D(bx, 0), D(bx, 2), D(bx, 3)
    byy, bzz = D(by, 2), D(bz, 3)
    lap = [by, bz]
    sq = by**2 + bz**2
    Dtb = bt - b * bx  # d_t b - b d_x b
    R = sp.zeros(4, 4)
    R[0, 0] = -sq / 2 - bx**2 + bxt - b * bxx
    R[1, 1] = bx**2 - sq / 2 - bxt + b * bxx
    R[0, 1] = R[1, 0] = -(byy + bzz) / 2
    for i, bxa in enumerate((bxy, bxz)):
        R[0, 2 + i] = R[2 + i, 0] = bxa / 2
        d_Dtb = D(Dtb, 2 + i)
        second = 0 if mut == "form-drop" else d_Dtb / 2
        R[1, 2 + i] = R[2 + i, 1] = bx * lap[i] / 2 - second
        for j in range(2):
            R[2 + i, 2 + j] = lap[i] * lap[j] / 2
    return R.applyfunc(sp.expand)


# ---------------------------------------------------------------------------------------------------------------
# The by-hand reduction (sheared.md, "The reduction"), as a decision procedure on a frame Ricci tensor at a point
# ---------------------------------------------------------------------------------------------------------------
def reducer(Rm, w, eps, mut):
    """Predicted status of scalar_form for Rm, with w the unnormalised transverse direction of the frame. Valid when
    R(wp, w) = R(wp, wp) = 0. Returns 'n/a' if that, or in the case r0 = r1 = 0 the identity (I2), is not met."""
    if Rm.is_zero_matrix:
        return "trivial"
    wv = sp.Matrix(w)
    wp = sp.Matrix([-w[1], w[0]])
    T = Rm[2:, 2:]
    n2 = w[0]**2 + w[1]**2
    Rww = (wv.T * T * wv)[0]
    if (wp.T * T * wv)[0] != 0 or (wp.T * T * wp)[0] != 0:
        return "n/a"
    r0, r1 = (wp.T * Rm[2:, 0])[0], (wp.T * Rm[2:, 1])[0]
    if r0 == 0 and r1 == 0:  # u_3 = 0: the lemma, with (I2)
        return "algebraic" if n2 * (Rm[0, 0] + Rm[1, 1]) + 2 * Rww == 0 else "n/a"
    lam = Rww / n2
    if lam == 0:  # u_3 != 0 needs lam = R(e_2, e_2) = -s u_3^2 != 0
        return "algebraic"
    conds = {1: (wv.T * Rm[2:, 0])[0], 2: (wv.T * Rm[2:, 1])[0],
             3: (Rm[0, 0] + lam) * lam * n2 + r0**2, 4: (Rm[1, 1] - lam) * lam * n2 + r1**2,
             5: Rm[0, 1] * lam * n2 + r0 * r1}
    for k in conds:
        if mut == f"red-no-c{k}":
            conds[k] = 0
    if any(c != 0 for c in conds.values()):
        return "algebraic"
    if mut == "red-sign":
        return "pass"
    return "pass" if eps * (-lam) > 0 else "sign"


def structured_planted(w, alpha, beta, r, deform=None):
    """R = r (eta - m m^T) in the frame (E_0, E_1, e_2, e_3), m = (alpha, beta, 0, 1), alpha^2 + beta^2 = 2, rotated to
    the coordinate frame with e_2 = w / |w|, e_3 = wp / |w|. deform = (i, j, delta) adds delta to R_ij (and R_ji)."""
    m = sp.Matrix([alpha, beta, 0, 1])
    Rf = r * (ETA - m * m.T)
    if deform:
        i, j, d = deform
        Rf[i, j] += d
        if i != j:
            Rf[j, i] += d
    nw = sp.sqrt(sp.Integer(w[0])**2 + sp.Integer(w[1])**2)
    Rot = sp.eye(4)
    Rot[2:, 2:] = sp.Matrix([[w[0], -w[1]], [w[1], w[0]]]) / nw  # columns e_2, e_3
    return Rot * Rf * Rot.T


# ---------------------------------------------------------------------------------------------------------------
# Sample points: (profile, kind, (X, y, z), v, jets)
# ---------------------------------------------------------------------------------------------------------------
def qparts(p):
    qX, qy, qz = (sp.diff(QS, s).subs(dict(zip((X, y, z), p))) for s in (X, y, z))
    return qX, qy, qz, qy**2 + qz**2


TH, BP, RP = "tanh, sigma = 8, R = 1", "bump (1 - q / 4)^3", "1 / (1 + q^2)"
RJ, RJ1, RJ0, RJC1, RJC2, BEY = ("random jets", "random jets, f = 1", "random jets, f' = 0", "random jets, c1 = 0",
                                 "random jets, c2 = 0", "bump, beyond it (f = 0)")
# (profile, kind, (X, y, z), v)
POINTS = [
    (TH, "axis", (1, 0, 0), 1), (TH, "axis", (-1, 0, 0), 1), (TH, "z = 0 plane", (rat(3, 5), rat(4, 5), 0), rat(1, 2)),
    (TH, "critical line", (1, -rat(1, 4), 0), 3), (TH, "X = 0 plane", (0, rat(1, 3), rat(2, 3)), 1),
    (TH, "y = 0 plane", (rat(1, 3), 0, rat(2, 3)), 3), (TH, "general", (rat(1, 2), rat(1, 2), rat(1, 2)), 3),
    (TH, "general", (rat(4, 5), rat(1, 5), rat(2, 5)), 1), (TH, "general", (-rat(1, 5), rat(4, 5), rat(2, 5)), rat(1, 2)),
    (TH, "general", (rat(5, 6), rat(1, 2), rat(1, 6)), 1), (TH, "general", (-rat(4, 7), rat(1, 7), rat(4, 7)), 3),
    (TH, "z axis", (0, 0, rat(7, 10)), 1), (TH, "z axis", (0, 0, -rat(5, 7)), rat(3, 2)),
    (BP, "axis", (1, 0, 0), rat(7, 4)), (BP, "z = 0 plane", (rat(3, 5), rat(4, 5), 0), rat(7, 4)),
    (BP, "general", (rat(1, 2), rat(1, 2), rat(1, 2)), rat(7, 4)), (BP, "general", (-rat(1, 5), rat(4, 5), rat(2, 5)), rat(7, 4)),
    (BP, "general", (rat(3, 2), rat(1, 3), -rat(1, 2)), rat(7, 4)), (BP, "z axis", (0, 0, rat(1, 2)), rat(7, 4)),
    (RP, "general", (rat(4, 5), rat(1, 5), rat(2, 5)), 2), (RP, "X = 0 plane", (0, rat(1, 3), rat(2, 3)), 2),
    (RP, "general", (rat(1, 3), -rat(2, 5), rat(3, 4)), rat(3, 2)),
    (RJ, "general", (rat(1, 2), rat(1, 2), rat(1, 2)), rat(8, 3)), (RJ, "general", (rat(5, 6), rat(1, 2), rat(1, 6)), rat(5, 3)),
    (RJ1, "general", (rat(4, 5), rat(1, 5), rat(2, 5)), rat(3, 2)), (RJ0, "general", (rat(4, 5), rat(1, 5), rat(2, 5)), rat(5, 2)),
    (RJC1, "general", (rat(1, 2), rat(1, 2), rat(1, 2)), 2), (RJC2, "general", (rat(1, 2), rat(1, 2), rat(1, 2)), 2),
    (BEY, "general", (2, 1, rat(1, 2)), rat(7, 4)),
]


def sample_points():
    """The exact points of the scalar-form test, as (profile, kind, (X, y, z), v, jets in q)."""
    rng = random.Random(81)

    def nz():
        return rat(rng.choice([-1, 1]) * rng.randint(1, 9), rng.randint(1, 9))

    def jets_of(prof, p, q0):
        qX, qy, qz, n2 = qparts(p)
        if prof == TH:
            return na.tanh_jets(q0)
        if prof in (BP, BEY):
            return na.bump_jets(q0)
        if prof == RP:
            return na.rational_jets(q0)
        if prof == RJ:
            return [rat(rng.randint(-9, 9), rng.randint(1, 9)), nz(), nz()]
        if prof == RJ1:
            return [sp.Integer(1), nz(), nz()]
        if prof == RJ0:
            return [rat(2, 3), sp.Integer(0), rat(-3, 5)]
        if prof == RJC1:  # f'' chosen so that R(w, E_0) = 0
            return [rat(1, 2), rat(-3, 2), -rat(-3, 2) * qy / (2 * qX * n2)]
        G0, G1 = rat(1, 2), rat(-3, 2)  # RJC2: f'' chosen so that R(w, E_1) = 0
        return [G0, G1, (2 * G1**2 * qX * n2 - G1 * (1 - G0) * qy / 2) / ((1 - G0) * qX * n2)]
    pts = []
    for prof, kind, p, vv in POINTS:
        p = tuple(sp.Rational(c) for c in p)
        q0 = sp.nsimplify(QS.subs(dict(zip((X, y, z), p))))
        pts.append((prof, kind, p, sp.Rational(vv), jets_of(prof, p, q0)))
    return pts


def adapted_w(p, mut):
    """The frame's transverse direction at p: grad_perp q, or (y, z) for `circular`; (1, 0) where it vanishes."""
    if mut == "circular":
        w0 = (p[1], p[2])
    else:
        w0 = tuple(sp.diff(qexpr(mut), s).subs(dict(zip((X, y, z), p))) for s in (y, z))
    return w0 if w0 != (0, 0) else (sp.Integer(1), sp.Integer(0))


def circular_w(p, mut):
    w0 = (p[2], p[1]) if mut == "circ-swap" else (p[1], p[2])
    return w0 if w0 != (0, 0) else (sp.Integer(1), sp.Integer(0))


def resid(Rm, w0):
    """Exact (I1) and (I2) residues in the frame (E_0, E_1, w0, w0p), and the sizes after normalising to unit vectors."""
    i1, i2 = na.lemma_residues(Rm, w0)
    n2 = w0[0]**2 + w0[1]**2
    nw = sp.sqrt(n2)
    n1 = float(sp.sqrt(sum(e**2 for e in [i1[0] / nw, i1[1] / nw, i1[2] / n2, i1[3] / n2])))
    nF = float(sp.sqrt(sum(e**2 for e in Rm)))
    return i1, i2, n1, abs(float(i2 / n2)), nF


# ---------------------------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------------------------
def main():
    opts = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
    mut = opts.get("mutate")
    if mut is not None and mut not in MUTATIONS:
        sys.exit(f"unknown mutation {mut!r}. Known: {', '.join(MUTATIONS)}")
    secs = set("DXPLRZ")
    if mut is not None:
        secs = {CHECK_SECTION[t] for t in MUTATIONS[mut][0]}
        if secs & set("PZ"):
            secs.add("D")
    random.seed(81)
    t0 = time.time()
    print(f"sheared.py, issue #81. sympy {sp.__version__}, mpmath {mp.__version__}, python {sys.version.split()[0]}. "
          f"mutation: {mut or 'none'}" + (f" ({MUTATIONS[mut][1]})" if mut else "") + f". sections {''.join(sorted(secs))}")
    lapse = mut == "lapse"
    q = qexpr(mut)
    qy, qz = sp.diff(q, y), sp.diff(q, z)
    wad = (y, z) if mut == "circular" else (qy, qz)  # the adapted frame's transverse direction
    wS = (sp.diff(QS, y), sp.diff(QS, z))
    qXS = sp.diff(QS, X)
    Rt = None

    # D. Route D, the closed form and the identities, for general jets
    if "D" in secs:
        print(f"\n== D. Frame Ricci tensor of b = -v f(q), q = {q}, on the jets of f ==")
        Rt = route_D(q, mut)
        for a in range(4):
            for b in range(a, 4):
                print(f"  R_({a},{b}): {len(sp.Add.make_args(Rt[a, b]))} terms")
        print("  (frame E_0 = n, E_1 = d_x, E_2 = d_y, E_3 = d_z. G_k = f^(k)(q). The expressions are in the log of the "
              "run's raw folder, not printed here.)")
        Rc = closed_form_tensor(mut)
        bad = [(a, b) for a in range(4) for b in range(a, 4) if sp.expand(Rt[a, b] - Rc[a, b]) != 0]
        check("closed_form", not bad, "R_ab of route D equals the closed form of any b(t, x, y, z) in its derivatives, "
              "for b = -v f(q) with q the sheared one, as a polynomial identity in X, y, z, v, f, f', f''. Components "
              f"that differ: {bad}")
        D = ps.JetD(na.q_rule(QS, False))
        b0 = -v * G[0]
        bx = D(b0, 1)
        cr = sp.expand(D(bx, 2) * D(b0, 3) - D(bx, 3) * D(b0, 2))
        check("shear_cross", sp.expand(cr - 2 * v**2 * G[1]**2 * z) == 0,
              f"d_y(d_x b) d_z b - d_z(d_x b) d_y b - 2 v^2 f'^2 z = {sp.expand(cr - 2 * v**2 * G[1]**2 * z)}, so "
              "grad_perp d_x b is along grad_perp b iff v f' z = 0")
        i1, i2 = na.lemma_residues(Rt, wad)
        want = [v * G[1] * z, v**2 * G[1] * (1 - G[0]) * z, 0, 0]
        d1 = [sp.expand(a - b) for a, b in zip(i1, want)]
        check("shear_I1", all(e == 0 for e in d1), f"w = {tuple(wad)}: R(wp, E_0), R(wp, E_1), R(wp, w), R(wp, wp) minus "
              f"(v f' z, v^2 f' (1 - f) z, 0, 0) = {d1}")
        check("shear_I2", sp.expand(i2) == 0, f"|w|^2 (R_00 + R_11) + 2 R(w, w) = {sp.factor(i2)}")
        cw = (y, z) if mut != "circ-swap" else (z, y)
        c1, c2 = na.lemma_residues(Rt, cw)
        h, h1 = G[1] * (1 - G[0]), G[2] * (1 - G[0]) - G[1]**2
        cpw = z * (2 * y - X / 2)  # cp . grad_perp q, cp = (-z, y)
        qXs = 2 * X + y / 2
        want_c = [-v / 2 * (G[2] * qXs * cpw - G[1] * z / 2),
                  v**2 / 2 * (G[1]**2 * qXs * cpw - h1 * qXs * cpw + h * z / 2),
                  v**2 * G[1]**2 * cpw * (y * (2 * y + X / 2) + z * 4 * z) / 2,
                  v**2 * G[1]**2 * cpw**2 / 2]
        dc = [sp.expand(a - b) for a, b in zip(c1, want_c)] + [sp.expand(c2 + v**2 * G[1]**2 * cpw**2)]
        check("circular_closed", all(e == 0 for e in dc), f"circular frame w = (y, z): (I1) residues minus "
              f"[-(v/2)(f'' q_X cp.w - f' z/2), ..., (v^2 f'^2/2)(cp.w)^2] and (I2) + v^2 f'^2 (cp.w)^2, "
              f"cp.w = z (2y - X/2): {dc}")
        w = sp.Matrix(wad)
        Rw0, Rw1 = [sp.expand((w.T * Rt[2:, a])[0]) for a in (0, 1)]
        kv = 2 if mut == "red-coef" else 1
        ident = sp.expand(Rw1 - v * (1 - G[0]) * Rw0 - kv * v**2 * G[1]**2 * qXS * (wS[0]**2 + wS[1]**2))
        check("reduction_identity", ident == 0, f"R(w, E_1) - v (1 - f) R(w, E_0) - v^2 f'^2 q_X |w|^2 = {ident}")
        print("  Lemma (by hand, phantom_source.md): in any orthonormal frame, R(e_3, .) = 0 and R_00 + R_11 + 2 R_22 = 0\n"
              "  force R_ab = lam eta_ab + s u_a u_b to have R_ab = 0. Here R(e_3, .) = (v f' z, v^2 f' (1 - f) z, 0, 0) / |w|,\n"
              "  so (I1) fails where v f' z != 0, and the lemma does not apply there.")

    # X. Route D against routes E and T at random exact points
    if "X" in secs:
        print("\n== X. Route D against route E (Christoffel at a point) and route T (3+1 equations) ==")
        if Rt is None:
            Rt = route_D(q, mut)
        rng = random.Random(7)
        badE = badT = done = 0
        while done < 8:
            p = [rat(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(3)]
            vv, tv = rat(rng.randint(1, 9), rng.randint(1, 5)), rat(rng.randint(-5, 5), rng.randint(1, 4))
            jets = [rat(rng.randint(-9, 9), rng.randint(1, 9)) for _ in range(3)]
            if vv * jets[0] == 1:
                continue
            done += 1
            q0 = q.subs(dict(zip((X, y, z), p)))
            bex = explicit_b(vv, jets, q0, q)
            pt = (tv, p[0] + vv * tv, p[1], p[2])
            D_ = Rt.subs({X: p[0], y: p[1], z: p[2], v: vv, G[0]: jets[0], G[1]: jets[1], G[2]: jets[2]})
            badE += any(sp.simplify(e) != 0 for e in D_ - ps.ricci_at_point(bex, pt, lapse))
            badT += any(sp.simplify(e) != 0 for e in D_ - na.route_T(bex, pt))
        check("xcheck_E", badE == 0, f"route D differs from route E at {badE} of 8 exact points (random v, point, jets)")
        check("xcheck_T", badT == 0, f"route D differs from route T at {badT} of 8 exact points")

    # L. Planted scalar forms at the sample points (positions only)
    pts = sample_points() if secs & set("PL") else []

    def test(Rm, eps):
        if mut == "broken-test":
            return "pass", 0, sp.zeros(4, 4), sp.zeros(4, 1)
        return ps.scalar_form(Rm, ETA, eps)

    planted = []
    if "L" in secs or mut == "planted-in-data":
        rng = random.Random(8181)
        for _ in pts:
            lam_p = rat(rng.randint(-9, 9), rng.randint(1, 9))
            u_p = sp.Matrix([rat(rng.choice([-1, 1]) * rng.randint(1, 9), rng.randint(1, 9)) for _ in range(4)])
            if mut == "planted-flat":
                lam_p, u_p = sp.Integer(0), sp.zeros(4, 1)
            sgn = 1 if mut == "planted-sign" else -1
            planted.append((lam_p, u_p, lam_p * ETA + sgn * 8 * sp.pi * u_p * u_p.T))

    def premises(Rm, p):
        i1, i2 = na.lemma_residues(Rm, adapted_w(p, mut))
        return all(e == 0 for e in i1) and i2 == 0

    if "L" in secs:
        print("\n== L. Planted control at the sample points: R = lam eta - 8 pi u u passes, a deformed one fails ==")
        bad_pass = bad_sign = bad_def = bad_prem = 0
        for (prof, kind, p, vv, jets), (lam_p, u_p, Rp) in zip(pts, planted):
            st, lr, M, u = test(Rp, -1)
            ok = st == "pass" and lr == lam_p and (M + 8 * sp.pi * u_p * u_p.T).applyfunc(sp.simplify).is_zero_matrix \
                and (u * u.T - u_p * u_p.T).applyfunc(sp.simplify).is_zero_matrix
            bad_pass += not ok
            bad_sign += test(Rp, +1)[0] != "sign"
            mu = 0 if mut == "planted-null" else rat(1, 10)
            Rd = Rp + mu * sp.Matrix(4, 4, lambda i, j: 1 if (i, j) == (3, 3) else 0)
            bad_def += test(Rd, -1)[0] != "algebraic"
            bad_prem += premises(Rp, p)
        n = len(pts)
        check("planted_pass", bad_pass == 0, f"lam eta - 8 pi u u, eps = -1: status 'pass' with lam and u u recovered "
              f"exactly at {n - bad_pass} of {n} points")
        check("planted_sign", bad_sign == 0, f"the same form with eps = +1: status 'sign' at {n - bad_sign} of {n} points")
        check("planted_deform", bad_def == 0, f"the form plus (1/10) E_3 E_3 (rank 2): status 'algebraic' at "
              f"{n - bad_def} of {n} points")
        check("planted_premise", bad_prem == 0, f"the planted form violates (I1) or (I2) of the adapted frame at "
              f"{n - bad_prem} of {n} points (the lemma's contrapositive)")

    # P. The scalar-form test on the sheared Ricci tensor, point by point
    if "P" in secs:
        print("\n== P. Scalar-form test at exact points, eps = -1 and +1 (phantom_source.scalar_form) ==")
        need_sf = mut is None or any(t in ("points_outcome", "lemma_agree", "reduction_agree") for t in MUTATIONS[mut][0])
        print("  A = frame adapted to grad_perp q (w = (q_y, q_z)); C = circular frame (y, z). Both at points where the "
              "vector vanishes use (E_2, E_3).\n"
              "  |I1|, |I2|: residues of (I1) and (I2) for unit vectors, rel = sqrt(|I1|^2 + |I2|^2) / ||R||_F, "
              "||R||_F over the 16 frame components.")
        print(f"  {'#':>2}  {'profile':<24} {'kind':<13} {'v':>5} {'X':>6} {'y':>6} {'z':>6} {'q':>6} {'f_q':>9} {'R!=0':>5} "
              f"{'A |I1|':>9} {'A |I2|':>8} {'A rel':>8} {'C |I1|':>9} {'C |I2|':>9} {'C rel':>8}  eps=-1     eps=+1")
        outc, badE, bad_res, bad_prem_h = [], 0, 0, 0
        n_fail = n_agree = n_red = n_premh = 0
        n_I1 = 0
        for i, (prof, kind, p, vv, jets) in enumerate(pts, 1):
            sub = {X: p[0], y: p[1], z: p[2], v: vv, G[0]: jets[0], G[1]: jets[1], G[2]: jets[2]}
            Rm = Rt.subs(sub)
            q0 = q.subs(dict(zip((X, y, z), p)))
            jE = [jets[0], jets[1], jets[2] + (rat(1, 7) if mut == "jet-g2" else 0)]
            tv = rat(1, 3)
            Em = ps.ricci_at_point(explicit_b(vv, jE, q0, q), (tv, p[0] + vv * tv, p[1], p[2]), lapse)
            badE += any(e != 0 for e in Rm - Em)
            if mut == "planted-in-data":
                Rm = planted[i - 1][2]
            nz = not Rm.is_zero_matrix
            wA = adapted_w(p, mut)
            iA1, iA2, nA1, nA2, nF = resid(Rm, wA)
            iC1, iC2, nC1, nC2, _ = resid(Rm, circular_w(p, mut))
            # expected residues in the frame (E_0, E_1, wA, wAp), from the closed form with the sheared q
            exp_i1 = [vv * jets[1] * p[2], vv**2 * jets[1] * (1 - jets[0]) * p[2], 0, 0]
            ok_res = [sp.expand(a - b) for a, b in zip(iA1, exp_i1)] == [0] * 4 and iA2 == 0
            bad_res += not ok_res
            prem = all(e == 0 for e in iA1) and iA2 == 0
            n_premh += prem
            n_I1 += any(e != 0 for e in iA1)
            sts = [test(Rm, e)[0] for e in (-1, 1)] if need_sf else ["-", "-"]
            want = "trivial" if not nz else "algebraic"
            n_fail += any(s != want for s in sts)
            n_agree += any(prem and s != want for s in sts)
            if need_sf:
                red = [reducer(Rm, wA, e, mut) for e in (-1, 1)]
                n_red += any(r != s for r, s in zip(red, sts))
            outc.append(sts)
            print(f"  {i:>2}  {prof:<24} {kind:<13} {fmt(vv):>5} {fmt(p[0]):>6} {fmt(p[1]):>6} {fmt(p[2]):>6} "
                  f"{fmt(sp.nsimplify(q0)):>6} {float(jets[1]):>9.4g} {str(nz):>5} {nA1:>9.3g} {nA2:>8.3g} "
                  f"{(nA1**2 + nA2**2)**.5 / nF if nF else 0:>8.3g} {nC1:>9.3g} {nC2:>9.3g} "
                  f"{(nC1**2 + nC2**2)**.5 / nF if nF else 0:>8.3g}  {sts[0]:<10} {sts[1]}", flush=True)
        n = len(pts)
        check("points_xcheck", badE == 0, f"route D equals route E (explicit b, sympy's own derivatives) in all 16 frame "
              f"components at {n - badE} of {n} sample points")
        check("points_residues", bad_res == 0, f"at {n - bad_res} of {n} points the adapted-frame residues equal "
              "(v f' z, v^2 f' (1 - f) z, 0, 0) for (I1) and 0 for (I2). (I1) fails at "
              f"{n_I1} of {n} points, the premises of the lemma hold at {n_premh}")
        nw = sum(1 for (_, _, _, _, j), s in zip(pts, outc) if j[1] != 0 and s == ["algebraic"] * 2)
        nzero = sum(1 for s in outc if s == ["trivial"] * 2)
        if need_sf:
            check("points_outcome", n_fail == 0, f"{n} points, both signs ({2 * n} tests). R != 0 gives 'algebraic' (no "
                  f"lam makes R - lam eta of rank <= 1) and R = 0 gives 'trivial'. {nw} points with f'(q) != 0 are "
                  f"'algebraic' at both signs, {nzero} point is 'trivial', {n_fail} points differ from that")
            check("lemma_agree", n_agree == 0, f"at the {n_premh} points where (I1) and (I2) hold, scalar_form agrees "
                  f"with the lemma ('algebraic' if R != 0, 'trivial' if R = 0) at both signs: {n_agree} differ")
            check("reduction_agree", n_red == 0, f"the by-hand reduction predicts scalar_form's status at {n - n_red} "
                  f"of {n} points, both signs")

    # R. The by-hand reduction on planted matrices
    if "R" in secs:
        print("\n== R. The reduction on planted R = r (eta - m m^T), m = (alpha, beta, 0, 1), alpha^2 + beta^2 = 2 ==")
        cases = [((3, 4), (1, 1), rat(1, 2)), ((5, 12), (rat(7, 5), rat(1, 5)), 3), ((8, 15), (rat(1, 5), rat(7, 5)), rat(2, 7)),
                 ((4, -3), (-1, 1), 5)]
        bp = 0
        for w, (al, be), r in cases:
            Rp = structured_planted(w, al, be, r)
            st = [test(Rp, e) for e in (-1, 1)]
            red = [reducer(Rp, w, e, mut) for e in (-1, 1)]
            ok = st[0][0] == "pass" and st[0][1] == r and st[1][0] == "sign" and red == ["pass", "sign"]
            bp += not ok
        check("reduction_pass", bp == 0, f"scalar_form gives 'pass' (lam = r) at eps = -1 and 'sign' at eps = +1, and the "
              f"reduction predicts both, at {len(cases) - bp} of {len(cases)} planted matrices")
        defs = [(0, 2), (1, 2), (0, 0), (1, 1), (0, 1)]
        bd = 0
        names = {(0, 2): "c1", (1, 2): "c2", (0, 0): "c3", (1, 1): "c4", (0, 1): "c5"}
        det = []
        for w, (al, be), r in cases:
            for (i, j) in defs:
                Rd = structured_planted(w, al, be, r, (i, j, rat(1, 10)))
                s_ = [test(Rd, e)[0] for e in (-1, 1)]
                r_ = [reducer(Rd, w, e, mut) for e in (-1, 1)]
                good = s_ == ["algebraic"] * 2 and r_ == ["algebraic"] * 2
                bd += not good
                if not good:
                    det.append((w, names[(i, j)], s_, r_))
        check("reduction_deform", bd == 0, f"R_ij + 1/10 for (i, j) in (0,2) c1, (1,2) c2, (0,0) c3, (1,1) c4, (0,1) c5: "
              f"scalar_form and the reduction both return 'algebraic' at {len(cases) * len(defs) - bd} of "
              f"{len(cases) * len(defs)} deformed matrices. Differing: {det}")

    # Z. A point on the z axis where the sheared Ricci tensor is a scalar form
    if "Z" in secs:
        print("\n== Z. A point of the z axis (X = y = 0) where the sheared Ricci tensor is a phantom scalar form ==")
        mp.mp.dps = 60
        v0, z0 = mp.mpf(1), mp.mpf(1) / 2
        p0 = (3 * mp.sqrt(2) - 4) / 8  # root of 32 p^2 + 32 p - 1 = 0, p = v^2 (1 - f)^2
        a0 = mp.sqrt(p0) / v0
        g0 = 1 - a0
        g1 = -1 / (256 * v0**2 * a0 * z0**2)
        g2 = (-a0 / (64 * z0**2) - 6 * g1) / (16 * z0**2)
        out = {}
        for tag, d2 in (("pass", mp.mpf(1) / 10 if mut == "zaxis-off" else 0), ("deformed", mp.mpf(1) / 10)):
            sub = {X: 0, y: 0, z: sp.Float(str(z0), 60), v: sp.Float(str(v0), 60), G[0]: sp.Float(str(g0), 60),
                   G[1]: sp.Float(str(g1), 60), G[2]: sp.Float(str(g2 + d2), 60)}
            Rm = Rt.subs(sub).applyfunc(lambda e: sp.N(e, 60))
            lam_ = mp.mpf(str(Rm[3, 3]))  # w = (0, 4 z), so R(w, w) / |w|^2 = R_33
            M = mp.matrix(4, 4)
            for i in range(4):
                for j in range(4):
                    M[i, j] = mp.mpf(str(Rm[i, j])) - lam_ * (mp.mpf(-1 if i == 0 else 1) if i == j else 0)
            minors = [abs(M[i, k] * M[j, m] - M[i, m] * M[j, k]) for i in range(4) for j in range(i + 1, 4)
                      for k in range(4) for m in range(k + 1, 4)]
            diag_ok = all(M[i, i] <= mp.mpf(10)**-40 for i in range(4)) and sum(M[i, i] for i in range(4)) < 0
            out[tag] = (max(minors), diag_ok, lam_)
        print(f"  v = 1, z = 1/2, f = {mp.nstr(g0, 20)}, f' = {mp.nstr(g1, 20)}, f'' = {mp.nstr(g2, 20)}, lam = R(w, w)/|w|^2 = "
              f"{mp.nstr(out['pass'][2], 20)}")
        check("zaxis_pass", out["pass"][0] < mp.mpf(10)**-40 and out["pass"][1] and out["deformed"][0] > mp.mpf(10)**-6,
              f"at 60 digits the largest 2x2 minor of R - lam eta is {mp.nstr(out['pass'][0], 5)} (needs < 1e-40) with "
              f"every diagonal entry <= 1e-40 and trace < 0 ({out['pass'][1]}), so R - lam eta = -8 pi u u with eps = -1; with "
              f"f'' + 1/10 it is {mp.nstr(out['deformed'][0], 5)} (needs > 1e-6)")

    failed = [nm for nm, ok in RESULTS if not ok]
    print(f"\nchecks: {len(RESULTS) - len(failed)} passed, {len(failed)} failed {failed}. {time.time() - t0:.0f} s")
    if mut is not None:
        want = MUTATIONS[mut][0]
        print(f"control {mut}: target checks {want}, fired {[nm for nm in want if nm in failed]}")
    return failed


def numbered(path):
    stem, n, cand = path[:-4], 1, path
    while os.path.exists(cand):
        n += 1
        cand = f"{stem}.{n}.log"
    return cand


if __name__ == "__main__":
    failed = main()
    controls_ok = True
    opts = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
    if "mutate" not in opts and "no-controls" not in opts:
        print("\n== Controls: each mutation must make its target checks fail ==", flush=True)
        for name, (want, what) in MUTATIONS.items():
            p = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True,
                               timeout=3600)
            if "logdir" in opts:
                with open(numbered(f"{opts['logdir']}/{name}.log"), "w") as fh:
                    fh.write(p.stdout + p.stderr)
            fails = [ln.split(":")[0][7:] for ln in p.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = all(w in fails for w in want) and p.returncode == 1
            controls_ok &= fired
            print(f"[{'FIRED' if fired else 'SILENT'}] --mutate={name} ({what}): exit {p.returncode}, failed {fails}",
                  flush=True)
            if p.returncode not in (0, 1):
                print(p.stderr[-2000:])
        print(f"controls: {'all fired' if controls_ok else 'SOME DID NOT FIRE'}")
    sys.exit(0 if not failed and controls_ok else 1)

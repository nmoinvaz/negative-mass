#!/usr/bin/env python3
"""The Hawking-Ellis type on the surface D = 0 off the equator of the Alcubierre wall (issue #72).

he_type.md (#64) decides D > 0 (Type I) and D < 0 (Type IV) and leaves D = 0 to a rank test, and review_he.md (#71) found
200 random points of D = 0 off the equator all Type II. This script decides the type on D = 0 off the equator exactly.

With a = f', b = f'', F = 1 - f, u = 1/r, c = cos(theta), s = sin(theta), and off the axis (s != 0), off the equator
(c != 0), with v != 0 and (a, b) != (0, 0), the point is on D = 0 and Type I exactly when

    T1:  b + a u = 0,                                  that is (r f')' = 0
    T2:  v^2 (r a + F)(r a + c^2 F) + c^2 = 0
    T3:  2 r a + (1 + c^2) F != 0.

Everything else on D = 0 is Type II or III. The rank test of he_type.py's `he_type_exact` runs at exact rational points.

Routes. The frame Ricci tensor at a point comes from review_he_derive.py's route C (Cartan's structure equations, #71),
which shares no code with he_type.py. The rank test is `he_type.he_type_exact` (#64). The closed forms are typed here,
checked equal to route C as polynomial identities, and the proofs below use only the typed forms.

    python3 notes/warp/he_surface.py                 all checks, then every control (exit 1 if any check fails)
    python3 notes/warp/he_surface.py --mutate=NAME   one control (its target checks must fail, exit 1)
    python3 notes/warp/he_surface.py --no-controls   the checks only
"""
import argparse
import math
import os
import random
import subprocess
import sys
import time
from fractions import Fraction as Fr

import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import he_type as ht  # noqa: E402  (#64: he_type_exact, the rank test)
import review_he_derive as rd  # noqa: E402  (#71: route C, the exact frame Ricci tensor)

v, F, a, b, u, c, s = rd.SYMS
lam = sp.Symbol("lam")
phi_, X_, gm_, y_, w_ = sp.symbols("phi X gamma y w")  # effective variables of the elimination

# name: (what it changes, checks that must fail)
CONTROLS = {
    "table-typo": ("R_(0,rho) of the typed table has the wrong factor, -v/2 becomes -v",
                   ["frame_ricci", "null_points"]),
    "tgt-tau": ("the target R'_(r,theta) = -c v F j_theta becomes -v F j_theta", ["rotated_frame"]),
    "tgt-cert": ("the certificate of necessity has (1 - v^2 F^2) for (1 + v^2 F^2)", ["necessity_T1"]),
    "tgt-DT1": ("the target D on T1 has 3 c^2 for 4 c^2", ["D_on_T1"]),
    "tgt-ideal": ("the quadratic of the elimination has w phi^2 for w (phi^2 + 1)", ["elimination"]),
    "loosen-T1": ("the description drops T1", ["elimination", "outside_description"]),
    "loosen-T2": ("the description drops T2", ["elimination", "outside_description"]),
    "loosen-T3": ("the description drops T3", ["outside_description"]),
    "no-rank": ("he_type_exact labels every D = 0 tensor Type II without the rank test",
                ["planted", "locus_points"]),
    "planted-deform": ("the planted Type I tensor is deformed before the check", ["planted"]),
    "null-nonnull": ("the eigenvector of the D = 0 sample is (1, 1/2, 1/2), not null", ["null_points"]),
    "tanh-sign": ("the scan tests F - r f' for F + r f'", ["tanh_scan"]),
}
MUT = set()
RESULTS = []
NSIMP = {"points": 0, "value_changed": 0}


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# Small exact tools
# ---------------------------------------------------------------------------------------------------------------
def pyth(m):
    """(c, s) = ((1 - m^2)/(1 + m^2), 2m/(1 + m^2)), a rational point of the unit circle."""
    return (1 - m * m) / (1 + m * m), 2 * m / (1 + m * m)


def red_sp(e):
    """Expand and reduce s^2 -> 1 - c^2 in a sympy expression."""
    return sp.expand(sp.rem(sp.expand(e), s ** 2 - 1 + c ** 2, s))


def ring_ev(p, vals):
    """Exact value of a ring element of review_he_derive.R at vals = (v, F, a, b, u, c, s), Fractions."""
    tot = Fr(0)
    for mon, co in p.terms():
        t = Fr(int(co.numerator), int(co.denominator))
        for x, e in zip(vals, mon):
            if e:
                t *= x ** e
        tot += t
    return tot


def vals_of(pt):
    return tuple(pt[x] for x in rd.SYMS)


def R3_at(RF, pt):
    """The 3x3 block (n, e_x, e_rho) of the frame Ricci tensor of ring elements RF at the point pt, Fractions."""
    vals = vals_of(pt)
    return [[ring_ev(RF[i][j], vals) for j in range(3)] for i in range(3)]


def to_R4(R3):
    R4 = sp.zeros(4, 4)
    for i in range(3):
        for j in range(3):
            R4[i, j] = sp.Rational(R3[i][j].numerator, R3[i][j].denominator)
    return R4


def M_of(R3):
    """M = eta R on (n, e_x, e_rho), nested lists of Fractions."""
    return ht.mixed_block(R3)


def cubic_D(R3):
    """(c2, c1, c0) and D of the block, by he_type.py's own `cubic` and `disc`."""
    c2, c1, c0 = ht.cubic(M_of(R3))
    return (c2, c1, c0), ht.disc(c2, c1, c0)


def repeated_root(c2, c1, c0):
    """(lam, multiplicity) of the repeated root of a cubic with D = 0, exact. A triple root has 2 c2^2 = 6 c1."""
    den = 2 * c2 ** 2 - 6 * c1
    if den == 0:
        return -c2 / 3, 3
    return (9 * c0 - c1 * c2) / den, 2


def rank_exact(M, lm):
    """Rank of M - lm I by minors in Fractions, with no use of sympy's rank."""
    A = [[M[i][j] - (lm if i == j else 0) for j in range(3)] for i in range(3)]
    if all(x == 0 for row in A for x in row):
        return 0
    minors = [A[i][k] * A[j][l] - A[i][l] * A[j][k] for i in range(3) for j in range(i + 1, 3)
              for k in range(3) for l in range(k + 1, 3)]
    if all(x == 0 for x in minors):
        return 1
    det = (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
           + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))
    return 3 if det != 0 else 2


def type_from_rank(mult, rk):
    """The lemma of he_type.md, for D = 0: Type I if rank = 3 - mult, II for a 2x2 Jordan block, III for a 3x3."""
    if rk == 3 - mult:
        return "I"
    if mult == 3 and rk == 2:
        return "III"
    return "II"


def classify(RF, pt):
    """(type by he_type_exact, D, rank by minors, multiplicity, lam) at an exact point of the Alcubierre family."""
    R3 = R3_at(RF, pt)
    (c2, c1, c0), D = cubic_D(R3)
    typ, Dsp = ht.he_type_exact(to_R4(R3))
    Dq = sp.Rational(D.numerator, D.denominator)
    # he_type_exact returns D after sp.nsimplify, which can replace an exact rational of large height by a nearby
    # radical (and send |D| < 1e-60 to 0). The sign and zero-ness must agree with the exact D, the value need not.
    assert (Dq > 0) == (Dsp > 0) and (Dq < 0) == (Dsp < 0), (Dq, Dsp)
    NSIMP["points"] += 1
    NSIMP["value_changed"] += (Dq != Dsp)
    if D != 0:
        return typ, D, None, None, None
    lm, mult = repeated_root(c2, c1, c0)
    rk = rank_exact(M_of(R3), lm)
    return typ, D, rk, mult, lm


# ---------------------------------------------------------------------------------------------------------------
# The typed closed forms
# ---------------------------------------------------------------------------------------------------------------
def typed_table():
    """Frame Ricci on (n, e_x, e_rho) from review_he.md section 1, typed here (sympy)."""
    R00 = v ** 2 * (F * (b * c ** 2 + a * u * s ** 2) - a ** 2 * (c ** 2 + s ** 2 / 2))
    R0x = v / 2 * (a * u * (1 + c ** 2) + b * s ** 2)
    Rxx = v ** 2 * (-F * (b * c ** 2 + a * u * s ** 2) + a ** 2 * (c ** 2 - s ** 2 / 2))
    Rrr = v ** 2 * a ** 2 * s ** 2 / 2
    R0r = -(v if "table-typo" in MUT else v / 2) * (b - a * u) * c * s
    Rxr = v ** 2 * c * s * (a ** 2 - F * (b - a * u) / 2)
    return sp.Matrix([[R00, R0x, R0r], [R0x, Rxx, Rxr], [R0r, Rxr, Rrr]])


# The frame-rotation O takes (n, e_x, e_rho) to (n, rhat, thetahat), rhat = (c, s), thetahat = (-s, c).
O3 = sp.Matrix([[1, 0, 0], [0, c, s], [0, -s, c]])
# Closed forms in the rotated frame (derived by hand from the table, checked in `rotated_frame`)
JR = v * c * a * u
JTH = -v * s * (b + a * u) / 2
P00 = v ** 2 * (F * (a * u * s ** 2 + b * c ** 2) - a ** 2 * (1 + c ** 2) / 2)
PRR = v ** 2 * (a ** 2 * (1 + c ** 2) / 2 - F * b * c ** 2)
PTT = -v ** 2 * a * s ** 2 * (a + 2 * F * u) / 2
PRT = v ** 2 * F * c * s * (b + a * u) / 2

P2 = v ** 2 * (a / u + F) * (a / u + c ** 2 * F) + c ** 2  # T2, with r = 1/u


def T1_val(p):
    return p[b] + p[a] * p[u]


def T2_val(p):
    r = 1 / p[u]
    return p[v] ** 2 * (r * p[a] + p[F]) * (r * p[a] + p[c] ** 2 * p[F]) + p[c] ** 2


def T3_val(p):
    r = 1 / p[u]
    return 2 * r * p[a] + (1 + p[c] ** 2) * p[F]


def in_description(p):
    """The description of the Type I subset, with the controls that loosen it."""
    ok = True
    if "loosen-T1" not in MUT:
        ok = ok and T1_val(p) == 0
    if "loosen-T2" not in MUT:
        ok = ok and T2_val(p) == 0
    if "loosen-T3" not in MUT:
        ok = ok and T3_val(p) != 0
    return ok


# ---------------------------------------------------------------------------------------------------------------
# Planted tensors, run through he_type_exact (the same code path as every point below)
# ---------------------------------------------------------------------------------------------------------------
def boost_rot(R3):
    """R3 -> L^T R3 L for a boost (3/5), a rotation (3/5, 4/5) and a second boost (5/13)."""
    def boost(beta):
        g = 1 / sp.sqrt(1 - beta ** 2)
        g = sp.Rational(g)
        return sp.Matrix([[g, g * beta, 0], [g * beta, g, 0], [0, 0, 1]])
    rot = sp.Matrix([[1, 0, 0], [0, sp.Rational(3, 5), -sp.Rational(4, 5)], [0, sp.Rational(4, 5), sp.Rational(3, 5)]])
    L = boost(sp.Rational(3, 5)) * rot * boost(sp.Rational(5, 13))
    assert (L.T * sp.diag(-1, 1, 1) * L) == sp.diag(-1, 1, 1)
    return L.T * R3 * L


def planted_tensors():
    eta = sp.diag(-1, 1, 1)
    lm = sp.Rational(3, 2)
    w1 = sp.Matrix([1, 2, 1])  # not null, eta(w, w) = 4
    sig1 = -sp.Rational(1, 2)
    if "planted-deform" in MUT:
        w1 = sp.Matrix([1, 2, 1]) + sp.Matrix([0, 0, sp.Rational(1, 3)])  # still rank 1 but now |w| shifts D
        sig1 = sig1  # (the deformation below also adds a rank 2 piece)
    pI = lm * eta + sig1 * w1 * w1.T
    if "planted-deform" in MUT:
        pI = pI + sp.diag(0, 0, sp.Rational(1, 5))
    p1 = sp.Matrix([0, 4, -3])
    p2 = sp.Matrix([3, -5, 0])
    pII = lm * eta + 2 * p1 * p1.T - p2 * p2.T           # Sigma u = 0 for the null u = (5, 3, 4), rank 2
    ell = sp.Matrix([5, 3, 4])                             # null: -25 + 9 + 16 = 0
    pII3 = lm * eta + ell * ell.T                          # triple root, rank 1
    lv, mv = sp.Matrix([1, 1, 0]), sp.Matrix([0, 0, 1])    # N k = -m, N m = l, N l = 0
    T = lv * mv.T + mv * lv.T
    pIII = lm * eta + eta * T * eta                        # triple root, a 3x3 Jordan block
    out = {"I (rank 1, w not null)": (pI, "I", 1, 2), "II (rank 2, null eigenvector)": (pII, "II", 2, 2),
           "II (triple root, rank 1)": (pII3, "II", 1, 3), "III (3x3 Jordan block)": (pIII, "III", 2, 3)}
    res = {}
    for k, (R3, typ, rk, mult) in out.items():
        res[k] = (boost_rot(R3), typ, rk, mult)
    return res


def check_planted():
    rows, bad = [], []
    for name, (R3, want, rk_want, mult_want) in planted_tensors().items():
        R4 = sp.zeros(4, 4)
        R4[:3, :3] = R3
        typ, D = ht.he_type_exact(R4)
        Rf = [[Fr(int(R3[i, j].p), int(R3[i, j].q)) for j in range(3)] for i in range(3)]
        (c2, c1, c0), Dq = cubic_D(Rf)
        lm, mult = repeated_root(c2, c1, c0) if Dq == 0 else (None, None)
        rk = rank_exact(M_of(Rf), lm) if Dq == 0 else None
        rows.append(f"{name}: type {typ}, D = {D}, rank {rk}, multiplicity {mult}")
        if typ != want or Dq != 0 or rk != rk_want or mult != mult_want:
            bad.append(name)
    # deformed copies must NOT keep the claimed type and rank: they leave D = 0 or change the rank
    deformed = []
    for name, (R3, want, rk_want, mult_want) in planted_tensors().items():
        R3d = R3.copy()
        R3d[0, 0] += sp.Rational(1, 7)
        R3d[2, 2] -= sp.Rational(1, 11)
        R4 = sp.zeros(4, 4)
        R4[:3, :3] = R3d
        typ, D = ht.he_type_exact(R4)
        Rf = [[Fr(int(R3d[i, j].p), int(R3d[i, j].q)) for j in range(3)] for i in range(3)]
        (c2, c1, c0), Dq = cubic_D(Rf)
        keeps = (Dq == 0 and typ == want and rank_exact(M_of(Rf), repeated_root(c2, c1, c0)[0]) == rk_want)
        deformed.append((name, typ, str(D), keeps))
        if keeps:
            bad.append("deformed " + name)
    check("planted", not bad,
          "; ".join(rows) + ". Deformed copies (R_00 + 1/7, R_22 - 1/11): "
          + "; ".join(f"{n} -> type {t}, D = {d}, claim {'holds' if k else 'fails'}" for n, t, d, k in deformed)
          + (f". Wrong: {bad}" if bad else ""))
    return not bad


# ---------------------------------------------------------------------------------------------------------------
# Symbolic checks
# ---------------------------------------------------------------------------------------------------------------
def check_frame_ricci(RF):
    T = typed_table()
    bad = [(i, j) for i in range(3) for j in range(3) if not rd.req(rd.to_ring(T[i, j]), RF[i][j])]
    e_phi = all(RF[3][j] == 0 and RF[j][3] == 0 for j in range(4))
    check("frame_ricci", not bad and e_phi,
          f"the 9 entries of the typed table equal route C's as polynomials in (v, F, a, b, u, c, s), and R(e_phi, .) = 0"
          + (f"; differ at {bad}" if bad else ""))


def check_rotated():
    T = typed_table()
    Rp = (O3 * T * O3.T).applyfunc(red_sp)
    want = {(0, 0): P00, (0, 1): JR, (0, 2): JTH, (1, 1): PRR, (1, 2): PRT, (2, 2): PTT}
    bad = [k for k, e in want.items() if red_sp(Rp[k] - e) != 0]
    tau_coef = -c * v * F if "tgt-tau" not in MUT else -v * F
    ids = {
        "R'_(r,theta) = -c v F j_theta": red_sp(PRT - tau_coef * JTH),
        "R'_00 + R'_rr = v^2 F a u s^2": red_sp(P00 + PRR - v ** 2 * F * a * u * s ** 2),
        "R'_(0,r)^2 = (v c a u)^2": red_sp(Rp[0, 1] ** 2 - JR ** 2),
        "R_00 + R_xx + 2 R_rr = 0": red_sp(T[0, 0] + T[1, 1] + 2 * T[2, 2]),
        "trace of the spatial block is rotation invariant": red_sp(Rp[1, 1] + Rp[2, 2] - T[1, 1] - T[2, 2]),
    }
    bad += [k for k, e in ids.items() if e != 0]
    check("rotated_frame", not bad,
          "in the frame (n, rhat, thetahat) the entries are j_r = v c a u, j_theta = -(v s/2)(b + a u), "
          "R'_(r,theta) = -c v F j_theta, and R'_00, R'_rr, R'_(theta,theta) as typed; the identities "
          + ", ".join(ids) + " hold" + (f"; fails: {bad}" if bad else ""))


def minors_rot(lm):
    """The six distinct 2x2 minors of Sigma' = R' - lm eta in the rotated frame, as sympy expressions."""
    A_, Bp, Cp = P00 + lm, PRR - lm, PTT - lm
    jr, jt, tau = JR, JTH, PRT
    return {"a": A_ * Bp - jr ** 2, "b": A_ * Cp - jt ** 2, "c": Bp * Cp - tau ** 2,
            "d": jr * tau - jt * Bp, "e": A_ * tau - jt * jr, "f": jr * Cp - tau * jt}, (A_, Bp, Cp)


def check_necessity():
    """rank(R - lam eta) <= 1 forces j_theta = 0, by two identities (the hand proof of T1)."""
    m, (A_, Bp, Cp) = minors_rot(lam)
    sign = -1 if "tgt-cert" in MUT else 1
    cert = JTH * v * c * a * u * (1 + sign * v ** 2 * F ** 2) + m["e"] + c * v * F * m["d"]
    e1 = red_sp(cert)
    # a = 0: j_theta^4 + m_a C'^2 + 2 A C' m_b - m_b^2 = 0
    e2 = red_sp((JTH ** 4 + m["a"] * Cp ** 2 + 2 * A_ * Cp * m["b"] - m["b"] ** 2).subs(a, 0))
    # the two relations the proof uses: m_d = -j_theta (B' + c v F j_r), m_e = -j_theta (c v F A + j_r)
    e3 = red_sp(m["d"] + JTH * (Bp + c * v * F * JR))
    e4 = red_sp(m["e"] + JTH * (c * v * F * A_ + JR))
    check("necessity_T1", e1 == 0 and e2 == 0 and e3 == 0 and e4 == 0,
          "for every lam, j_theta v c a u (1 + v^2 F^2) + m_e + c v F m_d = 0 identically, so rank 1 (all minors 0) "
          "gives j_theta v c a u (1 + v^2 F^2) = 0; at a = 0, j_theta^4 + m_a C'^2 + 2 A C' m_b - m_b^2 = 0 gives "
          "j_theta = 0. With c != 0, u != 0, v != 0, 1 + v^2 F^2 > 0, rank <= 1 forces j_theta = 0, that is "
          "b + a u = 0 off the axis")
    return e1 == 0 and e2 == 0


def check_sufficiency(RF):
    """On T1 the point splits into (n, rhat) and thetahat, and the rank-1 condition is the single equation T2."""
    sub = {b: -a * u}
    m, (A_, Bp, Cp) = minors_rot(PTT)  # lam = R'_(theta,theta), so C' = 0
    ids = {}
    ids["j_theta = R'_(r,theta) = 0 on T1"] = red_sp(JTH.subs(sub)) == 0 and red_sp(PRT.subs(sub)) == 0
    ids["A = -v^2 a (a + F u c^2)"] = red_sp(A_.subs(sub) + v ** 2 * a * (a + F * u * c ** 2)) == 0
    ids["B' = v^2 a (a + F u)"] = red_sp(Bp.subs(sub) - v ** 2 * a * (a + F * u)) == 0
    ids["A B' - j_r^2 = -v^2 a^2 u^2 P2"] = red_sp(sp.cancel(sp.together(
        (A_ * Bp - JR ** 2).subs(sub) + v ** 2 * a ** 2 * u ** 2 * P2))) == 0
    # the other five minors vanish identically on T1 at lam = R'_(theta,theta)
    others = [red_sp(m[k].subs(sub)) == 0 for k in "bcdef"]
    ids["the minors b, c, d, e, f vanish on T1 at lam = R'_(theta,theta)"] = all(others)
    # the third eigenvalue: tr M - 3 lam = B' - A = v^2 (a u) (2 r a + (1 + c^2) F)
    trM = -P00 + PRR + PTT
    ids["tr M - 3 lam = (v^2 a u)(2 r a + (1 + c^2) F) on T1"] = red_sp(
        (trM - 3 * PTT).subs(sub) - v ** 2 * a * u * (2 * a / u + (1 + c ** 2) * F)) == 0
    bad = [k for k, ok in ids.items() if not ok]
    check("sufficiency_T1", not bad, "on T1: " + "; ".join(ids) + (f". Fails: {bad}" if bad else ""))



def check_D_on_T1(Rr_ring):
    """D on T1 is v^6 a^6 u^6 (v^2 F^2 s^4 - 4 c^2) P2^2, with D from review_he_derive's ring."""
    Dring = rd.red(rd.Dof(Rr_ring))
    De = Dring.as_expr()
    four = 3 if "tgt-DT1" in MUT else 4
    target = v ** 6 * a ** 6 * u ** 6 * (v ** 2 * F ** 2 * s ** 4 - four * c ** 2) * P2 ** 2
    diff = red_sp(sp.expand(sp.cancel(sp.together(De.subs(b, -a * u) - target))))
    check("D_on_T1", diff == 0,
          "D(b = -a u) = v^6 a^6 u^6 (v^2 F^2 s^4 - 4 c^2) P2^2 as a polynomial identity, so on T1 the point is Type I "
          "where v F s^2 > 2 |c| (D > 0 or P2 = 0 with the double root), Type IV where v F s^2 < 2 |c|, and the "
          "Jordan block of the (n, rhat) pair on v F s^2 = 2 c")
    return diff == 0


def check_axis(RF, Rr_ring):
    """On the axis D = -4 a^2 v^6 Q_ax^2 / r^6 vanishes only for a = 0, where M is diagonal."""
    Dring = rd.red(rd.Dof(Rr_ring)).as_expr()
    sgn = 1 if "tgt-axis" not in MUT else -1
    Qax = a ** 2 + v ** 2 * (F * b - sgn * a ** 2) ** 2 / u ** 2
    target = -4 * a ** 2 * v ** 6 * u ** 6 * Qax ** 2
    d = sp.cancel(sp.together(Dring.subs({s: 0, c: 1}) - target))
    # exact axis points with a = 0: M diagonal, Type I
    rng = random.Random(5)
    bad = []
    for _ in range(20):
        pt = {v: Fr(rng.randint(1, 9), 4), F: Fr(rng.randint(1, 8), 8), a: Fr(0), b: Fr(rng.randint(-8, 8), 4),
              u: Fr(rng.randint(1, 6), rng.randint(1, 6)), c: Fr(1), s: Fr(0)}
        typ, D, rk, mult, lm = classify(RF, pt)
        M = M_of(R3_at(RF, pt))
        diag = all(M[i][j] == 0 for i in range(3) for j in range(3) if i != j)
        if D != 0 or typ != "I" or not diag:
            bad.append(pt)
    # and a != 0 gives D < 0 (Type IV) at exact axis points
    neg = 0
    for _ in range(20):
        pt = {v: Fr(rng.randint(1, 9), 4), F: Fr(rng.randint(1, 8), 8), a: Fr(rng.choice([-1, 1]) * rng.randint(1, 8), 4),
              b: Fr(rng.randint(-8, 8), 4), u: Fr(rng.randint(1, 6), rng.randint(1, 6)), c: Fr(1), s: Fr(0)}
        typ, D, rk, mult, lm = classify(RF, pt)
        neg += (typ == "IV" and D < 0)
    check("axis", d == 0 and not bad and neg == 20,
          "on the axis (c = 1, s = 0) D = -4 a^2 v^6 Q_ax^2/r^6 as a polynomial identity, so D = 0 only where v a = 0 "
          "(Q_ax = 0 forces a = 0), and there M = diag(-v^2 F b, -v^2 F b, 0) is Type I at 20 exact points; "
          f"a != 0 is Type IV at {neg} of 20 exact points")


def check_elimination():
    """Elimination of lam from the six minors in effective variables, by a Groebner basis, an independent route."""
    S = 1 - w_
    kap = gm_ + w_ * y_
    R00 = phi_ * kap - X_ * (1 + w_) / 2
    Rxx = -phi_ * kap + X_ * (3 * w_ - 1) / 2
    Rrr = X_ * S / 2
    p = gm_ + S * y_ / 2
    qq, tt = -y_ / 2, X_ - phi_ * y_ / 2
    cs2 = w_ * S
    A_, B_, C_ = R00 + lam, Rxx - lam, Rrr - lam
    ms = [A_ * B_ - p ** 2, A_ * C_ - cs2 * qq ** 2, B_ * C_ - cs2 * tt ** 2, A_ * tt - p * qq, B_ * qq - p * tt,
          C_ * p - cs2 * qq * tt]
    # effective variables are phi = v F, X = v^2 a^2, gamma = v a u, y = v (b - a u), w = c^2 (R entries per
    # review_he.md); y = 1 is the normalisation (R is linear in X, gamma, y, lam at fixed phi)
    eff = {phi_: v * F, X_: v ** 2 * a ** 2, gm_: v * a * u, y_: v * (b - a * u), w_: c ** 2}
    T = typed_table()
    ok_eff = all(red_sp((e - T[i, j]).subs(eff)) == 0 for e, (i, j) in [
        (R00, (0, 0)), (Rxx, (1, 1)), (Rrr, (2, 2)), (p, (0, 1)), (qq * c * s, (0, 2)),
        (c * s * (X_ - phi_ * y_ / 2), (1, 2))])
    m1 = [sp.expand(e.subs(y_, 1)) for e in ms]
    quad = 4 * X_ ** 2 - 2 * X_ * phi_ * (1 + w_) + w_ * (phi_ ** 2 + (0 if "tgt-ideal" in MUT else 1))
    claim = []
    if "loosen-T1" not in MUT:
        claim.append(2 * gm_ + 1)
    if "loosen-T2" not in MUT:
        claim.append(quad)
    claim.append(lam - (w_ - 1) * (X_ - phi_) / 2)
    order = (lam, gm_, X_, phi_, w_)
    G = sp.groebner(m1, *order, order="lex")
    Gc = sp.groebner(claim, *order, order="lex")
    h = (phi_ ** 2 + 1) * (w_ - 1)
    nec = all(G.contains(sp.expand(h * g)) for g in claim)       # h * claim in the ideal of the minors
    suf = all(Gc.reduce(mm)[1] == 0 for mm in m1)                 # the minors in the ideal of the claim
    # y = 0 branch (q = 0): normalise gamma = 1, the ideal is (lam, X, phi^2 + 1), no real point
    ms0 = [sp.expand(e.subs(y_, 0).subs(gm_, 1)) for e in ms]
    G0 = sp.groebner(ms0, lam, X_, phi_, w_, order="lex")
    y0 = (sorted(map(str, G0.exprs)) == sorted(map(str, [lam, X_, phi_ ** 2 + 1])))
    # translate: 2 gamma + y = v (b + a u) and Q = 4 v^2 a^2 u^2 P2 on T1
    t1 = red_sp((2 * gm_ + y_).subs(eff) - v * (b + a * u)) == 0
    Qy = (4 * X_ ** 2 - 2 * X_ * phi_ * y_ * (1 + w_) + w_ * y_ ** 2 * (phi_ ** 2 + 1)).subs(eff).subs(b, -a * u)
    t2 = red_sp(sp.cancel(sp.together(Qy - 4 * v ** 2 * a ** 2 * u ** 2 * P2))) == 0
    check("elimination", ok_eff and nec and suf and y0 and t1 and t2,
          "with y = 1, the minors' ideal in (lam, gamma, X, phi, w) satisfies h g in I(minors) for each generator g of "
          "<2 gamma + y, 4X^2 - 2 X phi y (1 + w) + w y^2 (phi^2 + 1), lam - (w - 1)(X - phi y)/2> and h = (phi^2 + 1) "
          f"(w - 1), and the minors lie in that ideal (necessity {nec}, sufficiency {suf}); for y = 0 the ideal is "
          f"(lam, X, phi^2 + 1), no real point ({y0}). 2 gamma + y = v (b + a u) and, on T1, the quadratic is "
          f"4 v^2 a^2 u^2 T2 ({t1 and t2}); effective forms equal the typed table ({ok_eff}). So off the axis "
          "(w != 1), the real solutions with y != 0 are exactly T1 and T2")
    return nec and suf


# ---------------------------------------------------------------------------------------------------------------
# Exact points
# ---------------------------------------------------------------------------------------------------------------
def locus_point(mc, vv, k, r, sign):
    """A rational point of T1, T2: F and a r from the conic d^2 = F^2 s^4 - 4 c^2/v^2, d = (c/v)(k - 1/k)."""
    cc, ss = pyth(mc)
    Fv = cc * (k + 1 / k) / (vv * ss * ss)
    d = (cc / vv) * (k - 1 / k)
    x = (-(1 + cc * cc) * Fv + sign * d) / 2           # x = r a
    aa = x / r
    if "locus-T2-shift" in MUT:
        aa = aa * Fr(101, 100)
    return {v: vv, F: Fv, a: aa, b: -aa / r, u: 1 / r, c: cc, s: ss}


def check_locus_points(RF):
    rng = random.Random(72)
    mcs = [Fr(1, 3), Fr(1, 2), Fr(2, 3), Fr(1, 5), Fr(3, 4), Fr(2, 5), Fr(1, 4)]
    n, nI, bad = 0, 0, []
    ranks, mults = {}, {}
    worst_T3 = None
    for _ in range(60):
        mc = rng.choice(mcs)
        k = rng.choice([Fr(2), Fr(3), Fr(3, 2), Fr(5, 2), Fr(1, 2), Fr(4), Fr(1, 3), Fr(7, 3)])
        vv = Fr(rng.randint(1, 6), rng.randint(1, 3))
        r = Fr(rng.randint(1, 8), rng.randint(1, 4))
        pt = locus_point(mc, vv, k, r, rng.choice([1, -1]))
        typ, D, rk, mult, lm = classify(RF, pt)
        n += 1
        t1, t2, t3 = T1_val(pt), T2_val(pt), T3_val(pt)
        R3 = R3_at(RF, pt)
        lam_th = None
        # the double eigenvalue is R'_(theta,theta)
        sub = {v: pt[v], F: pt[F], a: pt[a], b: pt[b], u: pt[u], c: pt[c], s: pt[s]}
        lam_th = Fr(int(sp.Rational(PTT.subs(sub)).p), int(sp.Rational(PTT.subs(sub)).q))
        ok = (typ == "I" and D == 0 and rk == 1 and mult == 2 and t1 == 0 and t2 == 0 and t3 != 0 and lm == lam_th
              and in_description(pt) and type_from_rank(mult, rk) == "I")
        nI += ok
        ranks[rk] = ranks.get(rk, 0) + 1
        mults[mult] = mults.get(mult, 0) + 1
        if not ok:
            bad.append((pt, typ, D, rk, mult))
    check("locus_points", not bad and nI == n,
          f"{nI} of {n} exact rational points with T1 = T2 = 0 and T3 != 0, placed by v F s^2 = c (k + 1/k) and "
          f"r a = [-(1 + c^2) F +- (c/v)(k - 1/k)]/2 (c, s from a Pythagorean pair), have D = 0 exactly, "
          f"he_type_exact says Type I, rank(M - lam) = {sorted(ranks)} by minors with multiplicity {sorted(mults)}, "
          "and the double eigenvalue is R'_(theta,theta)")
    return not bad


def null_point(rng, RF):
    """A rational point on D = 0 placed by a null eigenvector: (R - lam eta) e = 0 with e = (1, m1, m2) null forces a
    Jordan block, so D = 0. The system is linear in (X, gamma, y, lam) at fixed (c, e, phi), a ray, scaled so X is a
    square. Returns a point of (v, F, a, b, u, c, s), or None."""
    mc = Fr(rng.randint(1, 7), rng.randint(2, 9))
    if not 0 < mc < 1:
        return None
    cc, ss = pyth(mc)
    tu = Fr(rng.randint(1, 9), rng.randint(1, 9))
    m1, m2 = (1 - tu * tu) / (1 + tu * tu), 2 * tu / (1 + tu * tu)
    if "null-nonnull" in MUT:
        m1, m2 = Fr(1, 2), Fr(1, 2)
    ph = Fr(rng.randint(-9, 9), rng.randint(1, 5))
    if ph == 0:
        return None
    w = cc * cc
    # rows of (R - lam eta) e = 0 in the unknowns (X, gamma, y, lam)
    # unknowns (X, gamma, y, lam); R_00 = phi (gamma + w y) - X (1 + w)/2, R_0x = gamma + (1 - w) y/2,
    # R_0rho = -y c s/2, R_xx = -phi (gamma + w y) + X (3 w - 1)/2, R_xrho = c s (X - phi y/2), R_rr = X (1 - w)/2
    cs_ = cc * ss
    row0 = [-(1 + w) / 2, ph + m1, ph * w + m1 * (1 - w) / 2 - m2 * cs_ / 2, Fr(1)]
    row1 = [m1 * (3 * w - 1) / 2 + m2 * cs_, 1 - ph * m1, (1 - w) / 2 - ph * w * m1 - m2 * cs_ * ph / 2, -m1]
    row2 = [cs_ * m1 + m2 * (1 - w) / 2, Fr(0), -cs_ / 2 - cs_ * ph * m1 / 2, -m2]
    Mx = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in row] for row in (row0, row1, row2)])
    ns = Mx.nullspace()
    if len(ns) != 1:
        return None
    X0, g0, y0, l0 = [Fr(int(q.p), int(q.q)) for q in ns[0]]
    if X0 == 0 or g0 == 0:
        return None
    if X0 < 0:
        X0, g0, y0, l0 = -X0, -g0, -y0, -l0
    kk = Fr(rng.randint(1, 5), rng.randint(1, 5))
    tau = X0 * kk * kk
    X, gam, yy, lm_ = tau * X0, tau * g0, tau * y0, tau * l0
    alpha = X0 * kk
    if gam / alpha < 0:
        alpha = -alpha
    vv = Fr(rng.randint(1, 12), rng.randint(1, 4))
    uu = gam / alpha
    if uu <= 0:
        return None
    aa = alpha / vv
    return {v: vv, F: ph / vv, a: aa, b: yy / vv + aa * uu, u: uu, c: cc, s: ss}


def check_null_points(RF):
    rng = random.Random(7172)
    cnt, bad = {}, []
    ranks = {}
    n_noT1 = 0
    tries = 0
    nt3 = 0
    while sum(cnt.values()) < 250 and tries < 20000:
        tries += 1
        pt = null_point(rng, RF)
        if pt is None:
            continue
        typ, D, rk, mult, lm = classify(RF, pt)
        key = (typ, rk, mult) if D == 0 else (typ, "D != 0")
        cnt[key] = cnt.get(key, 0) + 1
        if D == 0 and T1_val(pt) != 0:
            n_noT1 += 1
        if D != 0 or typ not in ("II", "III") or in_description(pt):
            bad.append((key, pt))
        if D == 0 and typ == "III":
            nt3 += 1
    total = sum(cnt.values())
    check("null_points", not bad and total >= 250,
          f"{total} exact rational points placed on D = 0 by a null eigenvector (R - lam eta) e = 0, e = (1, m1, m2) null, "
          f"with a ray scaled so X = v^2 a^2 is a square; every one has D = 0 exactly, and (type, rank of M - lam, "
          f"multiplicity): {dict(sorted(cnt.items(), key=str))}; T1 fails at {n_noT1} of them, and the description "
          f"fails at all {total}" + (f". Violations: {bad[:3]}" if bad else ""))
    return not bad


def check_designated(RF):
    """Type II points outside the description that satisfy all but one condition, so loosening it is caught."""
    out = []
    rng = random.Random(3)
    # (a) T1, T3, not T2: on T1 with v F s^2 = 2 c the (n, rhat) pair is a Jordan block, D = 0
    for mc, vv, aa, r in [(Fr(1, 3), Fr(2), Fr(-1, 2), Fr(3, 2)), (Fr(1, 2), Fr(3), Fr(-3, 5), Fr(2)),
                          (Fr(2, 3), Fr(1), Fr(1, 4), Fr(1, 3))]:
        cc, ss = pyth(mc)
        Fv = 2 * cc / (vv * ss * ss)
        pt = {v: vv, F: Fv, a: aa, b: -aa / r, u: 1 / r, c: cc, s: ss}
        out.append(("T1, T3, not T2", pt))
    # (b) T1, T2, not T3: k = 1, the double root of the (n, rhat) pair meets the thetahat eigenvalue
    for mc, vv, r in [(Fr(1, 3), Fr(2), Fr(3, 2)), (Fr(1, 2), Fr(3), Fr(1)), (Fr(2, 3), Fr(1), Fr(5, 2))]:
        pt = locus_point(mc, vv, Fr(1), r, 1)
        out.append(("T1, T2, not T3", pt))
    rows, bad = [], []
    for label, pt in out:
        typ, D, rk, mult, lm = classify(RF, pt)
        t = (T1_val(pt) == 0, T2_val(pt) == 0, T3_val(pt) != 0)
        ok = D == 0 and typ == "II" and in_description(pt) is False
        rows.append((label, typ, rk, mult, t))
        if not ok:
            bad.append((label, typ, D, rk, mult, t, in_description(pt)))
    check("outside_description", not bad,
          "Type II points on D = 0 that satisfy two of T1, T2, T3 and fail the third, exact rational points: "
          + "; ".join(f"[{lb}] type {ty}, rank {rk}, multiplicity {mu}, (T1, T2, T3) = {t}" for lb, ty, rk, mu, t in rows)
          + ". The description rejects each, and fails if the failed condition is dropped"
          + (f". Violations: {bad}" if bad else ""))
    return not bad


# ---------------------------------------------------------------------------------------------------------------
# The tanh wall
# ---------------------------------------------------------------------------------------------------------------
def scan_roots(jets, rgrid):
    """Roots of (r f')' = f' + r f'' on a grid of radii, with (r, r f', F, F + r f') at each."""
    out = []
    g = lambda r: (lambda j: j[1] + r * j[2])(jets(r))
    prev_r, prev = rgrid[0], g(rgrid[0])
    for r in rgrid[1:]:
        cur = g(r)
        if prev * cur < 0:
            r0 = mp.findroot(g, (prev_r, r), solver="illinois", tol=1e-20, verify=False)
            f0, f1, _ = jets(r0)
            E = (1 - f0) + (-(r0 * f1) if "tanh-sign" in MUT else r0 * f1)
            out.append((r0, r0 * f1, 1 - f0, E))
        prev_r, prev = r, cur
    return out


def check_tanh_scan():
    mp.mp.dps = 30
    # the jets of he_type.py against mpmath's numerical derivatives of the profile, at three radii
    sig, Rb = mp.mpf(3), mp.mpf(1)
    prof = lambda r: (mp.tanh(sig * (r + Rb)) - mp.tanh(sig * (r - Rb))) / (2 * mp.tanh(sig * Rb))
    jerr = 0
    for r in (mp.mpf("0.4"), mp.mpf("1.1"), mp.mpf("1.9")):
        f, f1, f2 = ht.tanh_jets(r, sig, Rb, lib=mp)
        jerr = max(jerr, abs(f1 - mp.diff(prof, r)), abs(f2 - mp.diff(prof, r, 2)))
    # planted truth: a profile with a point of T1 and F + r f' > 0 must be found non-empty by the same scan
    eps = mp.mpf(1) / 10
    # f = eps sech^2 r, f' = -2 eps sech^2 tanh, f'' = 2 eps sech^2 (2 tanh^2 - sech^2)... written out below
    def planted_jets(r):
        sh2 = 1 / mp.cosh(r) ** 2
        th = mp.tanh(r)
        return eps * sh2, -2 * eps * sh2 * th, 2 * eps * sh2 * (2 * th ** 2 - sh2)
    pr = scan_roots(planted_jets, [mp.mpf(i) / 100 for i in range(1, 600)])
    planted_ok = bool(pr) and all(E > 0 for _, _, _, E in pr)
    # the scan on the tanh wall, sigma R_b from 0.01 to 316
    sigs = [mp.mpf(10) ** (mp.mpf(i) / 12 - 2) for i in range(0, 12 * 4 + 1)]
    worst, nroots, missing = None, [], []
    for sg in sigs:
        rmax = max(mp.mpf(12), 30 / sg)
        grid = [mp.mpf(10) ** (-3 + mp.mpf(i) * (mp.log10(rmax) + 3) / 1500) for i in range(1501)]
        roots = scan_roots(lambda r: ht.tanh_jets(r, sg, mp.mpf(1), lib=mp), grid)
        nroots.append(len(roots))
        if not roots:
            missing.append(float(sg))
        for r0, rf1, Fv, E in roots:
            if worst is None or E > worst[0]:
                worst = (E, float(sg), float(r0), float(rf1), float(Fv))
    empty = worst is not None and worst[0] < 0 and not missing
    check("tanh_scan", jerr < 1e-15 and planted_ok and empty,
          f"tanh jets of he_type.py agree with mpmath's derivatives to {mp.nstr(jerr, 3)}. Planted: f = eps sech^2 r "
          f"(eps = 1/10) has {len(pr)} root(s) of (r f')' with F + r f' = {[mp.nstr(E, 4) for *_, E in pr]} > 0, so "
          f"the scan finds a non-empty T1 and T2 there. Tanh wall, sigma R_b in {len(sigs)} log-spaced values from 0.01 "
          f"to 100: (r f')' = 0 has {sorted(set(nroots))} root(s) on r > 0, and F + r f' at them is at most "
          f"{mp.nstr(worst[0], 4)} (sigma R_b = {worst[1]:.4g}, r = {worst[2]:.4g}, r f' = {worst[3]:.4g}, "
          f"F = {worst[4]:.4g}), always negative, so T2 (which needs F + r f' > 0) has no real solution"
          if worst else "no roots found")
    return empty


# ---------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mutate", default=None)
    ap.add_argument("--no-controls", action="store_true")
    args = ap.parse_args()
    if args.mutate:
        if args.mutate not in CONTROLS:
            sys.exit(f"unknown mutation {args.mutate}; known: {', '.join(CONTROLS)}")
        MUT.add(args.mutate)
        if args.mutate == "no-rank":
            ht.MUT = "no-rank"
    t0 = time.time()
    Rr = rd.ring_frame(rd.build_C())
    check_planted()
    check_frame_ricci(Rr)
    check_rotated()
    check_necessity()
    check_sufficiency(Rr)
    check_D_on_T1(Rr)
    check_axis(Rr, Rr)
    check_elimination()
    check_locus_points(Rr)
    check_null_points(Rr)
    check_designated(Rr)
    check_tanh_scan()
    bad = sorted(n for n, ok in RESULTS if not ok)
    print(("FAILED: " + ",".join(bad)) if bad else f"all {len(RESULTS)} checks passed in {time.time() - t0:.0f} s")
    if args.mutate or args.no_controls:
        sys.exit(1 if bad else 0)
    code = 1 if bad else 0
    print("\ncontrols (each is the same script with --mutate=NAME and must fail its target checks)")
    for name, (what, targets) in CONTROLS.items():
        p = subprocess.run([sys.executable, "-u", __file__, f"--mutate={name}"], capture_output=True, text=True)
        line = [ln for ln in p.stdout.splitlines() if ln.startswith("FAILED:")]
        got = set(line[-1][len("FAILED: "):].split(",")) if line else set()
        fired = p.returncode == 1 and set(targets) <= got
        print(f"[{'FIRED' if fired else 'DID NOT FIRE'}] --mutate={name}: {what}; exit {p.returncode}; "
              f"checks that fail: {', '.join(sorted(got)) or 'none'}", flush=True)
        code = code or (0 if fired else 1)
    sys.exit(code)


if __name__ == "__main__":
    main()

"""Step 2 of the independent review of #72 (issue #101): the Ricci tensor and the type test, with planted truth first.

Checks
  planted_types   own type test on six planted tensors (perfect fluid, three distinct real roots, Jordan 2x2 plus a simple
                  root, null dust plus pressure, 3x3 Jordan block, complex pair), each boosted, rotated and boosted again
                  with rational Lorentz matrices, in 3x3 and, with a decoupled fourth axis, 4x4, and a deformed copy of each
                  of the first five that must get a different type
  ricci_numeric   the symbolic frame Ricci (jet of the shift) against 50-digit finite differences of the lab-frame metric
                  -dt^2 + (dx - v f(r_s) dt)^2 + dy^2 + dz^2 at t = 0, three points, tanh profile
  frame_table     the typed table of he_surface.md section 1 minus my frame Ricci, modulo c^2 + s^2 = 1
  equator_axis    D of my 3x3 block against the closed forms of he_type.md on the equator and on the axis

Run: python3 notes/warp/review_he_surface_ricci.py [--no-controls | --mutate=NAME]
Exit 0 only if every check passes and every control fires. Written without reading he_surface.py, he_type.py,
review_he_derive.py or phantom_source.py.
"""
import sys
import time
from fractions import Fraction as Fr

import mpmath as mp
import sympy as sp

import review_he_surface_lib as L
from review_he_surface_lib import a, b, c, F, s, u, v, x

MUT = ""
T0 = time.time()
mp.mp.dps = 60  # before any constant below


def set_mut(name):
    global MUT
    MUT = name
    L.MUT = name


# ----------------------------------------------------------------------------------------------------------------------
# planted types
# ----------------------------------------------------------------------------------------------------------------------
def lorentz3():
    """boost along x (beta 3/5), rotation in (x, y) by (5/13, 12/13), boost along y (beta 8/17); rational entries."""
    def boost(i, j, g, gb):
        m = sp.eye(3)
        m[0, 0] = m[i, i] = g
        m[0, i] = m[i, 0] = gb
        return m
    B1 = boost(1, 0, sp.Rational(5, 4), sp.Rational(3, 4))
    cs, sn = sp.Rational(5, 13), sp.Rational(12, 13)
    Rot = sp.Matrix([[1, 0, 0], [0, cs, -sn], [0, sn, cs]])
    B2 = boost(2, 0, sp.Rational(17, 15), sp.Rational(8, 15))
    Lm = B1 * Rot * B2
    eta = sp.diag(-1, 1, 1)
    assert (Lm.T * eta * Lm).applyfunc(sp.expand) == eta
    return Lm


def planted_list():
    """name, covariant R in a rest frame (3x3), perturbation dR (3x3), 4th-axis eigenvalue, expected type."""
    eta = sp.diag(-1, 1, 1)
    mu, p = sp.Rational(3, 2), sp.Rational(1, 5)
    E01 = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    diagp = sp.diag(sp.Rational(1, 7), 0, -sp.Rational(1, 11))
    kk = sp.Matrix([[1, -1, 0], [-1, 1, 0], [0, 0, 0]])  # k_a k_b with k_a = (-1, 1, 0), k null
    III = sp.Matrix([[0, 0, -1], [0, 0, 1], [-1, 1, 0]])  # l_a m_b + m_a l_b
    return [
        ("perfect fluid (3/2, 1/5)", sp.diag(mu, p, p), E01, p, "I"),
        ("three distinct real roots", sp.diag(2, -1, 4), E01, sp.Rational(7, 3), "I"),
        ("Jordan 2x2 plus simple root", eta * sp.diag(2, 2, 5) + 3 * kk, sp.diag(sp.Rational(1, 7), 0, 0), 5, "II"),
        ("null dust plus pressure", p * eta + 2 * kk, diagp, p, "II"),
        ("3x3 Jordan block", sp.Rational(1, 3) * eta + III, diagp, sp.Rational(1, 3), "III"),
        ("complex pair", E01, None, 0, "IV"),
    ]


def to_frac(M):
    return [[Fr(int(e.p), int(e.q)) for e in M.row(i)] for i in range(M.rows)]


def planted_types(h):
    Lm = lorentz3()
    eta = sp.diag(-1, 1, 1)
    ok_all, notes = True, []
    for name, R0, dR, ephi, want in planted_list():
        for variant in ("3x3", "4x4"):
            def build(Rr):
                Rp = (Lm.T * Rr * Lm).applyfunc(sp.expand)  # boosted, rotated, boosted covariant components
                M = eta * Rp
                if variant == "3x3":
                    return to_frac(M)
                M4 = sp.zeros(4, 4)
                M4[:3, :3] = M
                M4[3, 3] = ephi
                return to_frac(M4)
            got = L.classify(build(R0))["type"]
            ok = got == want
            line = "%s %s: got %s want %s" % (name, variant, got, want)
            if dR is not None:
                got2 = L.classify(build(R0 + dR))["type"]
                ok = ok and got2 != want
                line += "; deformed copy %s" % got2
            ok_all = ok_all and ok
            notes.append(line + ("" if ok else "  <-- FAIL"))
    for n in notes:
        print("   " + n)
    h.check("planted_types", ok_all, "six planted tensors, 3x3 and 4x4, each boosted, rotated, boosted; deformed copies change type")


# ----------------------------------------------------------------------------------------------------------------------
# the Ricci tensor against 50-digit finite differences of the lab-frame metric
# ----------------------------------------------------------------------------------------------------------------------
SIG, RB = mp.mpf(5) / 2, mp.mpf(6) / 5


def fprof(r):
    return (mp.tanh(SIG * (r + RB)) - mp.tanh(SIG * (r - RB))) / (2 * mp.tanh(SIG * RB))


def lab_metric(vv, X):
    """g_mu_nu at (t, x, y, z) = X for ds^2 = -dt^2 + (dx - v f dt)^2 + dy^2 + dz^2, r_s = |(x - v t, y, z)|."""
    t_, x_, y_, z_ = X
    rs = mp.sqrt((x_ - vv * t_) ** 2 + y_**2 + z_**2)
    w = vv * fprof(rs)
    g = mp.zeros(4, 4)
    g[0, 0] = -1 + (w**2 if MUT != "drop-w2" else 0)
    g[0, 1] = g[1, 0] = -w
    g[1, 1] = g[2, 2] = g[3, 3] = 1
    return g


D1 = [mp.mpf(1) / 280, -mp.mpf(4) / 105, mp.mpf(1) / 5, -mp.mpf(4) / 5]  # offsets -4..-1; antisymmetric for +1..+4
D2 = [-mp.mpf(1) / 560, mp.mpf(8) / 315, -mp.mpf(1) / 5, mp.mpf(8) / 5]  # offsets -4..-1; symmetric; centre -205/72


def d_one(func, X, i, h):
    out = None
    for k in range(1, 5):
        for sgn in (1, -1):
            Y = list(X)
            Y[i] += sgn * k * h
            val = func(Y)
            co = -sgn * D1[4 - k]
            out = val * co if out is None else out + val * co
    return out / h


def d_two_diag(func, X, i, h):
    out = func(list(X)) * (-mp.mpf(205) / 72)
    for k in range(1, 5):
        for sgn in (1, -1):
            Y = list(X)
            Y[i] += sgn * k * h
            out = out + func(Y) * D2[4 - k]
    return out / h**2


def numeric_frame_ricci(vv, X):
    """Frame components R(E_a, E_b), E = (n, d_x, d_y, d_z), from finite differences of the metric only."""
    h = mp.mpf(10) ** -8
    gfun = lambda Y: lab_metric(vv, Y)
    g = gfun(X)
    gi = g**-1
    dg = [d_one(gfun, X, i, h) for i in range(4)]
    ddg = [[None] * 4 for _ in range(4)]
    for i in range(4):
        ddg[i][i] = d_two_diag(gfun, X, i, h)
        for j in range(i + 1, 4):
            ddg[i][j] = ddg[j][i] = d_one(lambda Y, i=i: d_one(gfun, Y, j, h), X, i, h)
    dgi = [-(gi * dg[m] * gi) for m in range(4)]
    Gam = [[[sum(gi[k, l] * (dg[i][l, j] + dg[j][l, i] - dg[l][i, j]) for l in range(4)) / 2 for j in range(4)]
            for i in range(4)] for k in range(4)]

    def dGam(m, k, i, j):
        t1 = sum(dgi[m][k, l] * (dg[i][l, j] + dg[j][l, i] - dg[l][i, j]) for l in range(4)) / 2
        t2 = sum(gi[k, l] * (ddg[m][i][l, j] + ddg[m][j][l, i] - ddg[m][l][i, j]) for l in range(4)) / 2
        return t1 + t2
    R = mp.zeros(4, 4)
    for i in range(4):
        for j in range(4):
            e = 0
            for k in range(4):
                e += dGam(k, k, i, j) - dGam(j, k, i, k)
                if MUT != "no-gamma2":
                    for d in range(4):
                        e += Gam[k][k][d] * Gam[d][i][j] - Gam[k][j][d] * Gam[d][i][k]
            R[i, j] = e
    # frame: n = d_t + v f d_x (lab), then d_x, d_y, d_z
    r_s = mp.sqrt(X[1] ** 2 + X[2] ** 2 + X[3] ** 2)  # t = 0
    n = mp.matrix([1, vv * fprof(r_s), 0, 0])
    E = [n] + [mp.matrix([1 if k == i else 0 for k in range(4)]) for i in (1, 2, 3)]
    return mp.matrix(4, 4) + mp.matrix([[(E[i].T * R * E[j])[0] for j in range(4)] for i in range(4)])


def ricci_numeric(h):
    mp.mp.dps = 60
    out, _, _ = L.ricci_jet()
    fl = sp.lambdify((v, F, a, b, u, c, s), out, modules="mpmath")
    worst = mp.mpf(0)
    pts = [(mp.mpf("0.7"), mp.mpf("1.3"), mp.mpf("0.6")), (mp.mpf("1.4"), mp.mpf("0.9"), mp.mpf("-0.3")),
           (mp.mpf("2.2"), mp.mpf("2.1"), mp.mpf("0.8"))]
    lines = []
    for vv, r, cc in pts:
        ss = mp.sqrt(1 - cc**2)
        X = [mp.mpf(0), r * cc, r * ss, mp.mpf(0)]
        Rn = numeric_frame_ricci(vv, X)
        f0, f1, f2 = fprof(r), mp.diff(fprof, r, 1), mp.diff(fprof, r, 2)
        Rs = fl(vv, 1 - f0, f1, f2, 1 / r, cc, ss)
        dev = max(abs(Rn[i, j] - Rs[i, j]) / max(1, abs(Rs[i, j])) for i in range(4) for j in range(4))
        ephi = max(abs(Rn[3, j]) for j in range(4))
        worst = max(worst, dev, ephi)
        lines.append("v = %s r = %s c = %s: max rel. deviation %s, max |R(e_phi, .)| %s, R_rr^numeric = %s"
                     % (v_s(vv), v_s(r), v_s(cc), mp.nstr(dev, 3), mp.nstr(ephi, 3), mp.nstr(Rn[2, 2], 12)))
    for ln in lines:
        print("   " + ln)
    h.check("ricci_numeric", worst < mp.mpf(10) ** -30,
            "symbolic frame Ricci equals 50-digit second differences of the lab metric at 3 points, worst %s (tol 1e-30)"
            % mp.nstr(worst, 3))


def v_s(z):
    return mp.nstr(z, 4)


# ----------------------------------------------------------------------------------------------------------------------
# the typed table of he_surface.md section 1
# ----------------------------------------------------------------------------------------------------------------------
def typed_table():
    half = sp.Rational(1, 2) if MUT != "table-typo" else sp.Rational(1, 4)
    return {
        (0, 0): v**2 * (F * (b * c**2 + a * u * s**2) - a**2 * (c**2 + s**2 / 2)),
        (0, 1): v / 2 * (a * u * (1 + c**2) + b * s**2),
        (1, 1): v**2 * (-F * (b * c**2 + a * u * s**2) + a**2 * (c**2 - s**2 / 2)),
        (0, 2): -v / 2 * (b - a * u) * c * s,
        (2, 2): v**2 * a**2 * s**2 / 2,
        (1, 2): v**2 * c * s * (a**2 - half * F * (b - a * u)),
    }


def frame_table(h):
    R3, row = L.frame_block()
    diffs = {k: L.red(R3[k] - e) for k, e in typed_table().items()}
    trace = L.red(R3[0, 0] + R3[1, 1] + 2 * R3[2, 2])
    ok = all(d == 0 for d in diffs.values()) and trace == 0 and all(e == 0 for e in row)
    h.check("frame_table", ok,
            "six typed entries minus my frame Ricci reduce to %s; R_00 + R_xx + 2R_rr = %s; R(e_phi, .) = %s"
            % ([str(d) for d in diffs.values()], trace, row))


# ----------------------------------------------------------------------------------------------------------------------
# D on the equator and on the axis
# ----------------------------------------------------------------------------------------------------------------------
def disc_block(R3):
    M = L.mixed(R3)
    chi = sp.Poly(sp.expand((x * sp.eye(3) - M).det()), x)
    return chi, sp.expand(sp.discriminant(chi.as_expr(), x))


def equator_axis(h):
    R3, _ = L.frame_block()
    res = []
    # equator c = 0, s = 1: D = v^6 [r^2 v^2 a^4 - (a + r b)^2] Q_eq^2 / (16 r^6), Q_eq = 4 v^2 a^2 F (F + r a) + (a + r b)^2
    Re = R3.subs({c: 0, s: 1})
    _, De = disc_block(Re)
    r = 1 / u
    sg = -1 if MUT == "eq-typo" else 1
    Qeq = 4 * v**2 * a**2 * F * (F + sg * r * a) + (a + r * b) ** 2
    closed = v**6 * (r**2 * v**2 * a**4 - (a + r * b) ** 2) * Qeq**2 / (16 * r**6)
    d1 = sp.expand(De - closed)
    res.append(("equator", d1 == 0))
    # axis c = 1, s = 0: D = -4 a^2 v^6 Q_ax^2 / r^6, Q_ax = a^2 + r^2 v^2 (F b - a^2)^2
    Ra = R3.subs({c: 1, s: 0})
    _, Da = disc_block(Ra)
    sa = 1 if MUT == "axis-typo" else -1
    Qax = a**2 + r**2 * v**2 * (F * b + sa * a**2) ** 2
    d2 = sp.expand(Da - (-4 * a**2 * v**6 * Qax**2 / r**6))
    res.append(("axis", d2 == 0))
    h.check("equator_axis", all(ok for _, ok in res),
            "D of my 3x3 block equals he_type.md's closed forms: %s" % ", ".join("%s %s" % (n, "yes" if ok else "NO") for n, ok in res))


CHECKS = {
    "planted_types": planted_types,
    "ricci_numeric": ricci_numeric,
    "frame_table": frame_table,
    "equator_axis": equator_axis,
}
CONTROLS = {
    "no-rank": (["planted_types"], "the type test returns II at every D = 0"),
    "disc-sign": (["planted_types"], "the distinct-roots branch uses -D"),
    "mult-off": (["planted_types"], "rank is compared with n - m + 1 instead of n - m"),
    "drop-w2": (["ricci_numeric"], "the lab-frame metric loses the shift's square in g_tt"),
    "no-gamma2": (["ricci_numeric"], "the numeric Ricci drops the Gamma Gamma terms"),
    "table-typo": (["frame_table"], "the typed table's R_x,rho coefficient 1/2 -> 1/4"),
    "eq-typo": (["equator_axis"], "Q_eq: F + r a -> F - r a"),
    "axis-typo": (["equator_axis"], "Q_ax: (F b - a^2)^2 -> (F b + a^2)^2"),
}

if __name__ == "__main__":
    print("review_he_surface_ricci.py: sympy %s, mpmath %s" % (sp.__version__, mp.__version__), flush=True)
    H = L.Harness(CHECKS, CONTROLS, "ricci")
    H.main(sys.argv[1:], set_mut)

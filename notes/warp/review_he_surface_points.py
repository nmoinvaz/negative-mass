"""Step 3 of the independent review of #72 (issue #101): exact points, the axis, the Type III point, the tanh wall.

All frame Ricci entries come from my own review_he_surface_lib.ricci_jet and are evaluated in exact arithmetic (Fractions,
or Q(sqrt 113) for the Type III point) or at 60 digits. The type comes from my own test (lib.classify: sqf factors of the
characteristic polynomial, the rank of M - lam by Fraction Gaussian elimination, the size of the Jordan block).

  planted_types   the planted tensors of review_he_surface_ricci.py, run first through the same test
  exact_points    Type I points placed by my own search on T1, T2, T3; points failing exactly one of T1, T2, T3
  axis_points     on the axis D = 0 iff v a = 0, and there M is diagonal
  type3_exact     a Type III point in Q(sqrt 113) from the recipe of he_surface.md section 3, built with my own Ricci
  tanh_scan       the scan of F + r f' at the root of (r f')' = 0 for 97 values of sigma R_b, with a planted profile
  type3_tanh      Newton solves for a triple root of rank 2 on the tanh wall at sigma R_b = 1, 4, 8

Run: python3 notes/warp/review_he_surface_points.py [--no-controls | --mutate=NAME]
Written without reading he_surface.py, he_type.py, review_he_derive.py or phantom_source.py.
"""
import sys
from fractions import Fraction as Fr

import mpmath as mp
import sympy as sp

import review_he_surface_lib as L
import review_he_surface_ricci as RR
from review_he_surface_lib import QS, F, a, b, c, lam, s, u, v, x

MUT = ""
mp.mp.dps = 60


def set_mut(name):
    global MUT
    MUT = name
    L.MUT = name
    RR.MUT = name
    QS.d = 112 if name == "nf-modulus" else 113


_TERMS = {}


def entry_terms():
    if not _TERMS:
        R3, _ = L.frame_block()
        for i in range(3):
            for j in range(i, 3):
                _TERMS[i, j] = L.poly_terms(R3[i, j])
    return _TERMS


def R_at(pt):
    T = entry_terms()
    R = [[None] * 3 for _ in range(3)]
    for (i, j), tt in T.items():
        R[i][j] = R[j][i] = L.peval(tt, pt)
    return R


def M_at(pt):
    R = R_at(pt)
    return [[(-1 if i == 0 else 1) * R[i][j] for j in range(3)] for i in range(3)]


PYTH = [(Fr(3, 5), Fr(4, 5)), (Fr(5, 13), Fr(12, 13)), (Fr(8, 17), Fr(15, 17)), (Fr(7, 25), Fr(24, 25)),
        (Fr(20, 29), Fr(21, 29)), (Fr(12, 37), Fr(35, 37)), (Fr(9, 41), Fr(40, 41)), (Fr(-3, 5), Fr(4, 5))]


def conds(vv, FF, aa, bb, uu, cc, ss):
    r = 1 / uu
    t1 = bb + aa * uu
    t2 = vv**2 * (r * aa + FF) * (r * aa + cc**2 * FF) + cc**2
    t3 = 2 * r * aa + (1 + cc**2) * FF
    return t1, t2, t3


def issq(q):
    if q <= 0:
        return None
    n, d = sp.integer_nthroot(q.numerator, 2), sp.integer_nthroot(q.denominator, 2)
    return Fr(int(n[0]), int(d[0])) if n[1] and d[1] else None


def find_typeI(count):
    out = []
    Fs = [Fr(1, 2), Fr(2, 3), Fr(1), Fr(3, 2), Fr(3, 4), Fr(2), Fr(-1, 2), Fr(-3, 4)]
    rs = [Fr(1), Fr(2), Fr(3, 2), Fr(5, 2), Fr(4), Fr(5, 3)]
    ks = [Fr(p, q) for q in range(2, 13) for p in range(1, q)]
    for (cc, ss) in PYTH:
        for FF in Fs:
            for k in ks:
                xx = -FF * (cc**2 + (1 - cc**2) * k)  # between -F and -c^2 F
                den = (xx + FF) * (xx + cc**2 * FF)
                vv = issq(-cc**2 / den) if den != 0 else None
                if vv is None or 2 * xx + (1 + cc**2) * FF == 0:
                    continue
                r = rs[len(out) % len(rs)]
                aa = xx / r
                out.append((vv, FF, aa, -aa / r, 1 / r, cc, ss))
                break
            if len(out) >= count:
                return out
    return out


def classify_pt(pt):
    return L.classify(M_at(pt))


def D_poly_in_b(vv, FF, aa, uu, cc, ss):
    """D as an exact polynomial in b at fixed rational (v, F, a, u, c, s): discriminant of det(x - M(b))."""
    bs = sp.Symbol("bs")
    T = entry_terms()
    Rs = sp.zeros(3, 3)
    for (i, j), tt in T.items():
        val = sum(sp.Rational(cf.numerator, cf.denominator) * sp.Mul(*[sp.Rational(z.numerator, z.denominator) ** e if k != 3 else bs**e
                                                                       for k, (z, e) in enumerate(zip((vv, FF, aa, None, uu, cc, ss), mn))])
                  for mn, cf in tt)
        Rs[i, j] = Rs[j, i] = val
    M = L.mixed(Rs)
    chi = sp.Poly(sp.expand((x * sp.eye(3) - M).det()), x)
    c2, c1, c0 = chi.coeff_monomial(x**2), chi.coeff_monomial(x), chi.coeff_monomial(1)
    D = sp.expand(18 * c2 * c1 * c0 - 4 * c2**3 * c0 + c2**2 * c1**2 - 4 * c1**3 - 27 * c0**2)
    return sp.Poly(D, bs), bs


def numeric_type(pt_mp):
    """Rank structure of M at a 60-digit point with D ~ 0: smallest singular values of M - lam at the near-double root."""
    vv, FF, aa, bb, uu, cc, ss = pt_mp
    M = mp.matrix(frame_fn()(vv, FF, aa, bb, uu, cc, ss))
    tr = sum(M[i, i] for i in range(3))
    c1 = sum(M[i, i] * M[j, j] - M[i, j] * M[j, i] for i in range(3) for j in range(i + 1, 3))
    c0 = -mp.det(M)
    roots = mp.polyroots([1, -tr, c1, c0], maxsteps=500, extraprec=200)
    best = min(((abs(roots[i] - roots[j]), i, j) for i in range(3) for j in range(i + 1, 3)))
    lm = mp.re((roots[best[1]] + roots[best[2]]) / 2)
    sv = sorted(mp.svd_r(M - lm * mp.eye(3), compute_uv=False))
    return best[0], sv


def _jac_exact():
    """Jacobian of (T1, Pu) in (v, F, a, b, u, c) as callables on a point (v, F, a, b, u, c, s) of Fractions."""
    vs = (v, F, a, b, u, c)
    Pu_ = v**2 * (a + F * u) * (a + c**2 * F * u) + c**2 * u**2
    out = []
    for e in (b + a * u, Pu_):
        row = []
        for z in vs:
            terms = L.poly_terms(sp.diff(e, z))
            row.append(lambda *pt, terms=terms: L.peval(terms, pt))
        out.append(row)
    return out


def exact_points(h):
    pts = find_typeI(24)
    res = []
    nneg = nvf_neg = 0
    global JAC
    JAC = _jac_exact()
    for pt in pts:
        if MUT == "locus-shift":
            pt = (pt[0] * Fr(101, 100),) + tuple(pt[1:])
        t1, t2, t3 = conds(*pt)
        info = classify_pt(pt)
        D = info["D"]
        ok = (t1 == 0 and t2 == 0 and t3 != 0 and D == 0 and info["type"] == "I" and info["m"] == 2 and info["rank"] == 1)
        # (v F s^2)^2 >= 4 c^2 on T1 and T2, from Theta^2 = v^2 a^2 u^2 (v^2 F^2 s^4 - 4 c^2)
        ok = ok and (pt[0] * pt[1] * pt[6] ** 2) ** 2 >= 4 * pt[5] ** 2
        # T1 and T2 are independent equations: the Jacobian of (b + a u, Pu) in (v, F, a, b, u, c) has rank 2
        ok = ok and L.rank_frac([[g(*pt) for g in row] for row in JAC]) == 2
        res.append(ok)
        nneg += pt[1] < 0
        nvf_neg += pt[0] * pt[1] * pt[6] ** 2 < 0
    n_typeI = len(pts)
    # points failing exactly one condition
    desc = lambda t1, t2, t3: (MUT == "loosen-T1" or t1 == 0) and (MUT == "loosen-T2" or t2 == 0) and (MUT == "loosen-T3" or t3 != 0)
    outs = []
    # (a) T3 fails: T1, T2 hold, 2 r a + (1 + c^2) F = 0, v F s^2 = 2|c|
    cc, ss, FF, r = Fr(3, 5), Fr(4, 5), Fr(1, 2), Fr(3, 2)
    xx = -(1 + cc**2) * FF / 2
    vv = 2 * cc / (FF * ss**2)
    aa = xx / r
    pa = (vv, FF, aa, -aa / r, 1 / r, cc, ss)
    t = conds(*pa)
    i_a = classify_pt(pa)
    outs.append(("T3 fails", t, desc(*t) or False, i_a["type"] != "I" and i_a["D"] == 0 and i_a["m"] == 3 and i_a["rank"] == 1, i_a))
    # (b) T2 fails: T1, T3 hold, v F s^2 = 2|c| (the other factor of D on T1)
    FF, r = Fr(1), Fr(2)
    vv = 2 * cc / (FF * ss**2)
    xx = Fr(-3, 10)
    aa = xx / r
    pb = (vv, FF, aa, -aa / r, 1 / r, cc, ss)
    t = conds(*pb)
    i_b = classify_pt(pb)
    outs.append(("T2 fails", t, desc(*t) or False, i_b["type"] != "I" and i_b["D"] == 0 and i_b["m"] == 2 and i_b["rank"] == 2, i_b))
    # (c) T1 fails: T2, T3 hold at the first Type I point's parameters, b a real root of D(b) = 0 other than -a u
    vv, FF, aa, bb, uu, cc, ss = pts[0]
    Dp, bs = D_poly_in_b(vv, FF, aa, uu, cc, ss)
    roots = [rt for rt in Dp.nroots(n=60, maxsteps=500) if abs(sp.im(rt)) < sp.Float(10) ** -40]
    bT1 = sp.Rational((-aa * uu).numerator, (-aa * uu).denominator)
    cand = [sp.re(rt) for rt in roots if abs(sp.re(rt) - bT1) > sp.Rational(1, 100)]
    ok_c, line_c = False, "no real root of D(b) other than -a u"
    desc_c = None
    if cand:
        b0 = mp.mpf(str(cand[0]))
        gap, sv = numeric_type((mp.mpf(vv.numerator) / vv.denominator, mp.mpf(FF.numerator) / FF.denominator,
                                mp.mpf(aa.numerator) / aa.denominator, b0, mp.mpf(uu.numerator) / uu.denominator,
                                mp.mpf(cc.numerator) / cc.denominator, mp.mpf(ss.numerator) / ss.denominator))
        t1c = b0 + mp.mpf(aa.numerator) / aa.denominator * mp.mpf(uu.numerator) / uu.denominator
        _, t2c, t3c = conds(vv, FF, aa, Fr(0), uu, cc, ss)
        desc_c = (MUT == "loosen-T1" or abs(t1c) < mp.mpf(10) ** -30) and (MUT == "loosen-T2" or t2c == 0) and (MUT == "loosen-T3" or t3c != 0)
        rank2 = sv[1] > mp.mpf(10) ** -10  # second smallest singular value O(1): rank >= 2, a Jordan block
        ok_c = rank2 and gap < mp.mpf(10) ** -25
        line_c = "b = %s (-a u = %s), D-gap %s, singular values of M - lam: %s" % (
            mp.nstr(b0, 12), mp.nstr(-mp.mpf(aa.numerator) / aa.denominator * mp.mpf(uu.numerator) / uu.denominator, 12), mp.nstr(gap, 3),
            [mp.nstr(z, 3) for z in sv])
        outs.append(("T1 fails", (t1c, t2c, t3c), bool(desc_c), ok_c, {"type": "not I (rank 2)", "line": line_c}))
    else:
        outs.append(("T1 fails", (None, None, None), False, False, {"type": "?"}))
    for nm, t, d, ok, info in outs:
        print("   %s: conditions (T1, T2, T3) = %s, in the description %s, not Type I as required %s %s" % (
            nm, tuple(mp.nstr(z, 6) if not isinstance(z, Fr) else str(z) for z in t), d, ok, info.get("line", "type %s m %s rank %s" % (info.get("type"), info.get("m"), info.get("rank")))))
    # necessity sample: six generic parameter sets, roots of D(b) off -a u are never Type I
    samples = [(Fr(1), Fr(1, 2), Fr(-1, 3), Fr(3, 2), Fr(3, 5), Fr(4, 5)), (Fr(2, 3), Fr(2, 3), Fr(-1, 2), Fr(1), Fr(5, 13), Fr(12, 13)),
               (Fr(3, 2), Fr(1, 4), Fr(1, 3), Fr(2), Fr(-8, 17), Fr(15, 17)), (Fr(5, 4), Fr(3, 4), Fr(-1, 5), Fr(5, 4), Fr(4, 5), Fr(3, 5)),
               (Fr(7, 5), Fr(1, 3), Fr(-2, 7), Fr(3, 4), Fr(7, 25), Fr(24, 25)), (Fr(2), Fr(3, 5), Fr(-1, 4), Fr(5, 3), Fr(20, 29), Fr(21, 29))]
    nsamp, nI = 0, 0
    for vv, FF, aa, uu, cc, ss in samples:
        Dp, bs = D_poly_in_b(vv, FF, aa, uu, cc, ss)
        for rt in Dp.nroots(n=60, maxsteps=500):
            if abs(sp.im(rt)) > sp.Float(10) ** -40:
                continue
            b0 = mp.mpf(str(sp.re(rt)))
            q = lambda z: mp.mpf(z.numerator) / z.denominator
            gap, sv = numeric_type((q(vv), q(FF), q(aa), b0, q(uu), q(cc), q(ss)))
            if gap < mp.mpf(10) ** -25:
                nsamp += 1
                if sv[1] < mp.mpf(10) ** -10:
                    nI += 1
    ok_samp = nsamp >= 6 and nI == 0
    all_out = all(ok and not d for _, _, d, ok, _ in outs)
    print("   Type I points with F < 0 (where -F < r a < -c^2 F fails and v F s^2 >= 2|c| fails): %d, with v F s^2 < 0: %d" % (nneg, nvf_neg))
    h.check("exact_points", all(res) and n_typeI >= 20 and all_out and ok_samp and nneg >= 2,
            "%d points on T1, T2, T3 have D = 0 exactly, a double root, rank 1, (v F s^2)^2 >= 4 c^2 and independent T1, T2 (Type I %d of %d; %d with F < 0); "
            "the three points failing exactly one condition are not Type I and lie outside the description: %s; "
            "%d real roots of D(b) off T1 at six generic parameter sets, %d of them Type I"
            % (n_typeI, sum(res), n_typeI, nneg, all_out, nsamp, nI))


# ----------------------------------------------------------------------------------------------------------------------
def axis_points(h):
    ok_all = True
    notes = []
    # a = 0, v != 0: M = diag(-v^2 F b, -v^2 F b, 0), Type I
    for vv, FF, bb, uu in [(Fr(1), Fr(1, 2), Fr(3), Fr(1, 2)), (Fr(2), Fr(3, 4), Fr(-5, 3), Fr(2)), (Fr(3, 2), Fr(2), Fr(1, 7), Fr(1))]:
        for cc in (Fr(1), Fr(-1)):
            sgn = -1 if MUT == "axis-typo" else 1
            M = M_at((vv, FF, Fr(0), bb, uu, cc, Fr(0)))
            want = [[-sgn * vv**2 * FF * bb if i == j and i < 2 else Fr(0) for j in range(3)] for i in range(3)]
            info = L.classify(M)
            ok = M == want and info["D"] == 0 and info["type"] == "I"
            ok_all = ok_all and ok
    # a != 0, v != 0: D < 0 (Type IV), D = -4 a^2 v^6 Q_ax^2 u^6 exactly
    import random
    random.seed(7)
    nIV = 0
    for _ in range(20):
        vv, FF, aa, bb, uu = [Fr(random.randint(1, 9), random.randint(1, 5)) * random.choice((-1, 1)) for _ in range(5)]
        cc = random.choice((Fr(1), Fr(-1)))
        M = M_at((vv, FF, aa, bb, uu, cc, Fr(0)))
        info = L.classify(M)
        Qax = aa**2 + vv**2 * (FF * bb - aa**2) ** 2 / uu**2
        Dc = -4 * aa**2 * vv**6 * Qax**2 * uu**6
        ok = info["D"] == Dc and info["D"] < 0 and info["type"] == "IV"
        nIV += ok
        ok_all = ok_all and ok
    h.check("axis_points", ok_all,
            "a = 0 (6 points): M = diag(-v^2 F b, -v^2 F b, 0), D = 0, Type I; a != 0 (20 random exact points): D = -4 a^2 v^6 Q_ax^2 u^6 < 0, Type IV (%d of 20)" % nIV)


# ----------------------------------------------------------------------------------------------------------------------
def qs_charpoly(M):
    n = len(M)
    Mk = [[QS(1 if i == j else 0) for j in range(n)] for i in range(n)]
    cs = [QS(1)]
    for k in range(1, n + 1):
        AM = L.qs_matmul(M, Mk)
        ck = QS(0) - sum((AM[i][i] for i in range(n)), QS(0)) / k
        cs.append(ck)
        Mk = [[AM[i][j] + (ck if i == j else QS(0)) for j in range(n)] for i in range(n)]
    return cs


def type3_exact(h):
    d = QS.d
    phi = QS(Fr(-9, 16), Fr(1, 16))
    vv, cc, ss = QS(2), QS(Fr(3, 5)), QS(Fr(4, 5))
    FF = phi / vv
    e = [QS(1), QS(Fr(-7, 25)), QS(Fr(24, 25))]
    # R entries in the effective variables X = v^2 a^2, g = v a u, y = v (b - a u), linear in (X, g, y) at fixed (phi, c)
    def Rlin(X, g, y):
        e0 = phi * (y * cc**2 + g) - X * (cc**2 + ss**2 / 2)
        e2 = X * ss**2 / 2
        return [[e0, g * (1 + cc**2) / 2 + (y + g) * ss**2 / 2 + QS(0) * g, -y * cc * ss / 2],
                [None, QS(0) - e0 - 2 * e2, cc * ss * (X - phi * y / 2)],
                [None, None, e2]]
    # kernel of Sigma(lam) e = 0 in the unknowns (X, g, y, lam): column of each unknown from unit vectors
    def Sig_e(X, g, y, lm):
        R = Rlin(X, g, y)
        for i in range(3):
            for j in range(i):
                R[i][j] = R[j][i]
        eta = (-1, 1, 1)
        return [sum((R[i][j] * e[j] for j in range(3)), QS(0)) - lm * eta[i] * e[i] for i in range(3)]
    zero = QS(0)
    cols = [Sig_e(QS(1), zero, zero, zero), Sig_e(zero, QS(1), zero, zero), Sig_e(zero, zero, QS(1), zero), Sig_e(zero, zero, zero, QS(1))]
    A = [[cols[j][i] for j in range(4)] for i in range(3)]
    # row reduce, solve for the kernel with free variable lam = 1
    rows = [r[:] for r in A]
    piv = []
    ri = 0
    for col in range(4):
        pr = next((i for i in range(ri, 3) if not rows[i][col].iszero()), None)
        if pr is None:
            continue
        rows[ri], rows[pr] = rows[pr], rows[ri]
        rows[ri] = [z / rows[ri][col] for z in rows[ri]]
        for i in range(3):
            if i != ri and not rows[i][col].iszero():
                f = rows[i][col]
                rows[i] = [z1 - f * z2 for z1, z2 in zip(rows[i], rows[ri])]
        piv.append(col)
        ri += 1
    free = [j for j in range(4) if j not in piv]
    sol = [QS(0)] * 4
    sol[free[0]] = QS(1)
    for k, col in enumerate(piv):
        sol[col] = QS(0) - rows[k][free[0]]
    X0, g0, y0, l0 = sol
    ok_kernel = all(z.iszero() for z in Sig_e(X0, g0, y0, l0))
    # scale: X = X0^2, a = X0 / v, u = sigma g0, b = (y + g)/v
    sgn = 1 if float(g0) > 0 else -1
    X, g, y, lm = X0 * X0, X0 * g0, X0 * y0, X0 * l0
    aa = QS(sgn) * X0 / vv
    uu = QS(sgn) * g0
    bb = (y + g) / vv
    consist = (vv**2 * aa**2 == X) and (vv * aa * uu == g) and (vv * (bb - aa * uu) == y)
    pt = (vv, FF, aa, bb, uu, cc, ss)
    R = R_at(pt)
    # my R3 entries equal the effective-variable forms (they are the same polynomials), and Sigma(lam) e = 0
    eta = (-1, 1, 1)
    sig_e = [sum((R[i][j] * e[j] for j in range(3)), QS(0)) - lm * eta[i] * e[i] for i in range(3)]
    ok_sig = all(z.iszero() for z in sig_e)
    M = [[QS(eta[i]) * R[i][j] for j in range(3)] for i in range(3)]
    cs = qs_charpoly(M)
    c2, c1, c0 = cs[1], cs[2], cs[3]
    Dq = QS(18) * c2 * c1 * c0 - QS(4) * c2**3 * c0 + c2**2 * c1**2 - QS(4) * c1**3 - QS(27) * c0**2
    triple = (c2**2 == QS(3) * c1) and (c0 * QS(27) == c2**3)
    lt = QS(0) - c2 / 3
    N = [[M[i][j] - (lt if i == j else QS(0)) for j in range(3)] for i in range(3)]
    N2 = L.qs_matmul(N, N)
    N3 = L.qs_matmul(N2, N)
    n2nz = any(not z.iszero() for row in N2 for z in row)
    n3z = all(z.iszero() for row in N3 for z in row)
    rk = L.qs_rank(N)
    t1 = bb + aa * uu
    posr = float(uu) > 0
    eq8 = QS(8) * phi**2 + QS(9) * phi - QS(1)
    ok = ok_kernel and consist and ok_sig and Dq.iszero() and triple and n2nz and n3z and rk == 2 and not t1.iszero() and posr and eq8.iszero()
    print("   Q(sqrt %d): F = %.4f, a = %.4f, b = %.4f, r = %.4f, b + a u = %.4f, lam = %.5f" % (
        d, float(FF), float(aa), float(bb), 1 / float(uu), float(t1), float(lt)))
    h.check("type3_exact", ok,
            "e null, kernel of Sigma(lam) e: %s, consistent effective variables: %s, D = 0: %s, triple root (c2^2 = 3 c1, 27 c0 = c2^3): %s, N^2 != 0: %s, N^3 = 0: %s, "
            "rank N = %d, b + a u != 0: %s, r > 0: %s, 8 phi^2 + 9 phi - 1 = 0: %s" % (
                ok_kernel, consist, Dq.iszero(), triple, n2nz, n3z, rk, not t1.iszero(), posr, eq8.iszero()))


# ----------------------------------------------------------------------------------------------------------------------
def tanh_jets(sig, r):
    """f, f', f'' of f = [tanh sig (r + 1) - tanh sig (r - 1)] / (2 tanh sig), R_b = 1."""
    sech2 = lambda z: 1 / mp.cosh(z) ** 2
    zp, zm = sig * (r + 1), sig * (r - 1)
    n = 2 * mp.tanh(sig)
    f = (mp.tanh(zp) - mp.tanh(zm)) / n
    f1 = sig * (sech2(zp) - sech2(zm)) / n
    f2 = -2 * sig**2 * (sech2(zp) * mp.tanh(zp) - sech2(zm) * mp.tanh(zm)) / n
    return f, f1, f2


def planted_jets(r):
    sech2 = 1 / mp.cosh(r) ** 2
    th = mp.tanh(r)
    return sech2 / 10, -sech2 * th / 5, (2 * sech2 * th**2 - sech2**2) / 5


def scan_margin(jets, rmin, rmax, npts):
    """roots of g(r) = f' + r f'' on [rmin, rmax] (log grid, bisection to 50 digits); return [(r, F + r f', r f')]."""
    out = []
    rs = [mp.exp(mp.log(rmin) + (mp.log(rmax) - mp.log(rmin)) * k / npts) for k in range(npts + 1)]
    g = lambda r: jets(r)[1] + r * jets(r)[2]
    prev = g(rs[0])
    for k in range(1, len(rs)):
        cur = g(rs[k])
        if prev * cur < 0:
            r0 = mp.findroot(g, (rs[k - 1], rs[k]), solver="anderson", tol=mp.mpf(10) ** -45, maxsteps=200)
            f0, f1, _ = jets(r0)
            out.append((r0, (1 - f0) + r0 * f1, r0 * f1, 1 - f0))
        prev = cur
    return out


def t2_possible(F0, x0):
    """T1 and T2 can hold at a root with r f' = x0, F = F0 iff (x0 + F0)(x0 + c^2 F0) < 0 for some c^2 in (0, 1)."""
    if MUT == "tanh-sign":
        x0 = -x0
    return (F0 > 0 and -F0 < x0 < 0) or (F0 < 0 and 0 < x0 < -F0)


def tanh_scan(h):
    sigs = [mp.exp(mp.log(mp.mpf("0.01")) + (mp.log(100) - mp.log(mp.mpf("0.01"))) * k / 96) for k in range(97)]
    worst, worst_s, allone, anyok = mp.mpf("-inf"), None, True, False
    rows = []
    for sg in sigs:
        roots = scan_margin(lambda r, sg=sg: tanh_jets(sg, r), mp.mpf(10) ** -4, mp.mpf(10) ** 5, 3000)
        allone = allone and len(roots) == 1
        for r0, marg, x0, F0 in roots:
            if marg > worst:
                worst, worst_s = marg, (sg, r0, x0, F0)
            anyok = anyok or t2_possible(F0, x0)
    # planted profile f = sech^2(r)/10: its (r f')' root has F + r f' > 0 and r f' < 0
    pr = scan_margin(planted_jets, mp.mpf("1e-3"), mp.mpf(30), 3000)
    planted_found = any(t2_possible(F0, x0) for _, _, x0, F0 in pr)
    print("   97 values of sigma R_b in [0.01, 100], one root of (r f')' = 0 on r in [1e-4, 1e5] at every value: %s" % allone)
    print("   largest F + r f' at the root: %s at sigma R_b = %s, r = %s (r f' = %s, F = %s)" % (
        mp.nstr(worst, 6), mp.nstr(worst_s[0], 4), mp.nstr(worst_s[1], 6), mp.nstr(worst_s[2], 6), mp.nstr(worst_s[3], 6)))
    print("   planted profile f = sech^2(r)/10: roots %s" % [(mp.nstr(r0, 6), mp.nstr(m, 6)) for r0, m, _, _ in pr])
    h.check("tanh_scan", allone and not anyok and planted_found and worst < 0,
            "T2 needs F + r f' > 0 and r f' < 0 at the root; the tanh scan finds none (max F + r f' = %s < 0), the planted profile is found non-empty: %s"
            % (mp.nstr(worst, 5), planted_found))


# ----------------------------------------------------------------------------------------------------------------------
_FL = {}


def frame_fn():
    if "f" not in _FL:
        R3, _ = L.frame_block()
        _FL["f"] = sp.lambdify(L.VARS7, [[(-1 if i == 0 else 1) * R3[i, j] for j in range(3)] for i in range(3)], modules="mpmath")
    return _FL["f"]


def M_tanh(sig, r, vv, cc):
    f0, f1, f2 = tanh_jets(sig, r)
    ss = mp.sqrt(1 - cc**2)
    return mp.matrix(frame_fn()(vv, 1 - f0, f1, f2, 1 / r, cc, ss))


def cubic_coeffs(M):
    tr = sum(M[i, i] for i in range(3))
    c1 = sum(M[i, i] * M[j, j] - M[i, j] * M[j, i] for i in range(3) for j in range(i + 1, 3))
    return -tr, c1, -mp.det(M)


def type3_tanh(h):
    cases = [(1, mp.mpf("2.0656"), mp.mpf("1.1282"), mp.mpf("0.7768")), (4, mp.mpf("0.898"), mp.mpf("2.2089"), mp.mpf("0.2116")),
             (8, mp.mpf("0.916"), mp.mpf("2.0194"), mp.mpf("0.6437"))]
    ok_all = True
    for sg, r, v0, c0 in cases:
        eqs = lambda vv, cc: (lambda co: (co[0] ** 2 - 3 * co[1], 27 * co[2] - co[0] ** 3))(cubic_coeffs(M_tanh(mp.mpf(sg), r, vv, cc)))
        try:
            vv, cc = mp.findroot(eqs, (v0, c0), tol=mp.mpf(10) ** -45, maxsteps=100)
        except Exception as ex:  # a failed solve is a failed check
            print("   sigma R_b = %d: Newton failed: %s" % (sg, ex))
            ok_all = False
            continue
        if MUT == "type3-off":
            vv = vv * (1 + mp.mpf(10) ** -3)
        M = M_tanh(mp.mpf(sg), r, vv, cc)
        c2, c1, c0_ = cubic_coeffs(M)
        lm = -c2 / 3
        N = M - lm * mp.eye(3)
        sv = sorted(mp.svd_r(N, compute_uv=False), reverse=True)
        n3 = mp.norm(N * N * N, 1)
        n2 = mp.norm(N * N, 1)
        ok = sv[2] < mp.mpf(10) ** -30 and sv[1] > mp.mpf("1e-3") and n3 < mp.mpf(10) ** -30 and n2 > mp.mpf("1e-3")
        ok_all = ok_all and ok
        print("   sigma R_b = %d, r = %s: v = %s, cos theta = %s; singular values of M - lam %s; |N^3| = %s, |N^2| = %s" % (
            sg, mp.nstr(r, 5), mp.nstr(vv, 8), mp.nstr(cc, 8), [mp.nstr(z, 4) for z in sv], mp.nstr(n3, 3), mp.nstr(n2, 3)))
    h.check("type3_tanh", ok_all, "Newton at 60 digits on c2^2 = 3 c1, 27 c0 = c2^3 at fixed r, starting from rounded values, gives a 3x3 Jordan block "
                                  "(two singular values O(1), one < 1e-30, N^3 = 0, N^2 != 0) at sigma R_b = 1, 4, 8")


CHECKS = {
    "planted_types": RR.planted_types,
    "exact_points": exact_points,
    "axis_points": axis_points,
    "type3_exact": type3_exact,
    "tanh_scan": tanh_scan,
    "type3_tanh": type3_tanh,
}
CONTROLS = {
    "no-rank": (["planted_types", "exact_points"], "the type test returns II at every D = 0"),
    "loosen-T1": (["exact_points"], "the description of the Type I set drops T1"),
    "loosen-T2": (["exact_points"], "the description drops T2"),
    "loosen-T3": (["exact_points"], "the description drops T3"),
    "locus-shift": (["exact_points"], "v on the Type I points is scaled by 101/100"),
    "axis-typo": (["axis_points"], "the axis check expects the diagonal with the sign of b flipped"),
    "nf-modulus": (["type3_exact"], "the modulus of Q(sqrt 113) is changed to 112"),
    "tanh-sign": (["tanh_scan"], "F + r f' -> F - r f' in the T2 test of the scan"),
    "type3-off": (["type3_tanh"], "v of each Newton solution is scaled by 1 + 1e-3"),
}

if __name__ == "__main__":
    print("review_he_surface_points.py: sympy %s, mpmath %s" % (sp.__version__, mp.__version__), flush=True)
    H = L.Harness(CHECKS, CONTROLS, "points")
    H.main(sys.argv[1:], set_mut)

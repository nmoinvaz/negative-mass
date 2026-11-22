#!/usr/bin/env python3
"""The exact Type III set on D = 0 off the equator of the Alcubierre wall (issue #102).

Setting. Frame (n, e_x, e_rho) at the point r (c, s, 0) of the comoving metric -dt^2 + (dx + w dt)^2 + dy^2 + dz^2,
w = v (1 - f(r)), with a = f', b = f'', F = 1 - f, u = 1/r, and the 3x3 block M = eta R, whose characteristic polynomial is
x^3 + c2 x^2 + c1 x + c0 (e_phi decouples). A Type III point has a triple root with rank(M - lam) = 2. For v c s u != 0 and
(a, b) != (0, 0) the claim, with T = b + a u, E1 = c2^2 - 3 c1 and E2 = 27 c0 - c2^3, is

    Type III  <=>  E1 = 0,  E2 = 0,  T != 0.

Checks (each can fail and each has a control, `--mutate=NAME`):
  planted_types   Jordan 3x3, Jordan 2x2 + simple root, triple root of rank 1, fluid, three real roots, complex pair, each
                  boosted and rotated, through own `classify` and `he_type_exact`; E1 and E2 vanish exactly on the triple roots
  polynomials     E1, E2 as polynomials in (v, F, a, b, u, c); the frame table typed from he_surface.md section 1; the
                  effective form in (phi, A, B, u, c); E1 = (p - q)^2 and E2 = (p - q)^2 (8p + q) on a double root p, third q
  cert_rank       cofactor certificate: rank Sigma(lam) <= 1 forces T = 0 (the unrotated certificate of #101, restated)
  cert_T0         cofactor certificate: T = 0 with E1 = E2 = 0 forces rank <= 1, so rank 2 at a triple root forces T != 0
  cert_disc       27 D + E2^2 = 27 E1 G in (c2, c1, c0): on D = 0, E1 = 0 forces E2 = 0
  equator_axis    on the equator and the axis M is block diagonal, so no 3x3 Jordan block
  ray_family      a Type III point has a null kernel vector; on the null ray E1 = tau^2, E2 = tau^2 (9 lam + tau) with
                  tau = tr M - 3 lam, which gives the exact points
  exact_points    exact Type III points in Q(sqrt d), #72's point in Q(sqrt 113) among them: D = 0, triple root, rank 2
  type2_failures  Type II points of D = 0 failing exactly one condition, and a Type IV point failing only E2
  description     the description against the type, at exact points of every type
  tanh_jets       own mpmath jets of the tanh profile against he_type.tanh_jets and mp.diff
  tanh_planted    the resultant code finds a planted Type III point
  tanh_points     the Type III points of the tanh wall at Le's and Rodal's speeds, at 60 digits
  tanh_curve      the curve of Type III points over r at sigma R_b = 8 and 20
  tanh_enclosure  interval Newton enclosures of the points at the papers' speeds
  lab_ricci       the same points from second differences of the lab metric (own Christoffel code, no comoving jets)

Imports. From `review_he_surface_lib.py` (#101): `frame_block` (the frame Ricci tensor, by its own Christoffel and Ricci
code), `sigma`, `minors`, `red`, `classify` and `rank_frac` (its own type test), `QS` and `qs_rank` (arithmetic in Q(sqrt d)),
`Harness`. From `he_type.py` (#64, #103): `he_type_exact` for the type of a tensor and `tanh_jets` for the profile cross-check.
Not imported: `he_surface.py`, `review_he_surface_typeI.py`, `review_he_surface_points.py`. The cofactors of `cert_rank` are
#101's, restated and re-verified here; every other certificate is new. `lab_ricci` shares nothing with the library's Ricci code.

    python3 notes/warp/he_type3.py                 all checks, then every control (exit 1 if any check fails)
    python3 notes/warp/he_type3.py --mutate=NAME   one control (its target checks must all fail, exit 1)
    python3 notes/warp/he_type3.py --no-controls   the checks only
"""
import math
import os
import sys
import time
from fractions import Fraction as Fr

import mpmath as mp
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import review_he_surface_lib as L  # noqa: E402
from review_he_surface_lib import QS, a, b, c, F, lam, s, u, v, x  # noqa: E402

import he_type as HT  # noqa: E402

MUT = ""
T0 = time.time()


def set_mut(name):
    global MUT
    MUT = name
    L.MUT = name
    HT.MUT = name


ph, X, ga, yy, w = sp.symbols("phi X gamma y w")
e1s, e2s = sp.symbols("e1 e2")
Ap, Bp = sp.symbols("Ap Bp")
VARS7 = (v, F, a, b, u, c, s)
Tpol = b + a * u


def say(msg):
    print(msg, flush=True)


# ----------------------------------------------------------------------------------------------------------------------
# helpers: exact evaluation of polynomials at points of Q or Q(sqrt d)
# ----------------------------------------------------------------------------------------------------------------------
class PolyEval:
    """A sympy polynomial in `syms`, with Fraction coefficients, evaluated at points of Q or Q(sqrt d)."""

    def __init__(self, expr, syms):
        P = sp.Poly(sp.expand(expr), *syms)
        self.terms = [(mn, Fr(int(cf.p), int(cf.q))) for mn, cf in P.terms()]

    def __call__(self, pt):
        tot = 0
        for mn, cf in self.terms:
            t = cf
            for xx, e in zip(pt, mn):
                if e:
                    t = t * xx**e
            tot = tot + t
        return tot


def cubic_at(M):
    """(c2, c1, c0) of det(x - M) = x^3 + c2 x^2 + c1 x + c0 for a 3x3 nested list over Q or Q(sqrt d)."""
    c2 = -(M[0][0] + M[1][1] + M[2][2])
    c1 = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) + (M[0][0] * M[2][2] - M[0][2] * M[2][0]) + (
        M[1][1] * M[2][2] - M[1][2] * M[2][1])
    det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
           + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    return c2, c1, -det


def E12(c2, c1, c0):
    e2_ = 27 * c0 + c2**3 if MUT == "E2-typo" else 27 * c0 - c2**3
    e1_ = c2**2 - 4 * c1 if MUT == "E1-typo" else c2**2 - 3 * c1
    return e1_, e2_


def disc3(c2, c1, c0):
    return 18 * c2 * c1 * c0 - 4 * c2**3 * c0 + c2**2 * c1**2 - 4 * c1**3 - 27 * c0**2


def is0(z):
    return z == 0


def rank_any(M):
    """Rank over Q (own Fraction elimination) or over Q(sqrt d) (exact zero tests)."""
    if any(isinstance(e, QS) for row in M for e in row):
        return L.qs_rank(M)
    return L.rank_frac(M)


def matmul_any(A, B):
    n = len(A)
    return [[sum((A[i][k] * B[k][j] for k in range(n)), 0 * A[0][0]) for j in range(n)] for i in range(n)]


def shift(M, lm):
    return [[M[i][j] - (lm if i == j else 0) for j in range(3)] for i in range(3)]


def type_own(M):
    """Own Hawking-Ellis typing of a mixed 3x3 over Q or Q(sqrt d), using exact formulas for the repeated root.
    Returns dict(type, D, E1, E2, mult, lam, rank, N2zero, N3zero)."""
    c2, c1, c0 = cubic_at(M)
    D = disc3(c2, c1, c0)
    E1, E2 = E12(c2, c1, c0)
    out = dict(D=D, E1=E1, E2=E2, mult=1, lam=None, rank=None)
    if not is0(D):
        isfield = any(isinstance(z, QS) for z in (c2, c1, c0))
        if isfield:
            out["type"] = "?"
            return out
        out["type"] = "I" if D > 0 else "IV"
        return out
    if is0(E1):
        lm, mult = -c2 / 3, 3
    else:
        lm, mult = (9 * c0 - c2 * c1) / (2 * (c2**2 - 3 * c1)), 2
    N = shift(M, lm)
    rk = rank_any(N)
    N2 = matmul_any(N, N)
    N3 = matmul_any(N2, N)
    out.update(mult=mult, lam=lm, rank=rk, N2zero=all(is0(e) for r_ in N2 for e in r_),
               N3zero=all(is0(e) for r_ in N3 for e in r_))
    if mult == 3:
        out["type"] = {0: "I", 1: "II", 2: "III"}[rk]
    else:
        out["type"] = {1: "I", 2: "II"}[rk]
    return out


# ----------------------------------------------------------------------------------------------------------------------
# the algebra: Ricci block, invariants, effective variables
# ----------------------------------------------------------------------------------------------------------------------
_AL = {}


def typed_table():
    """The six frame entries of he_surface.md section 1 in (n, e_x, e_rho), typed here from the text."""
    sg = -1 if MUT == "table-typo" else 1
    R00 = v**2 * (F * (b * c**2 + a * u * s**2) - a**2 * (c**2 + s**2 / 2))
    R0x = sg * (v / 2) * (a * u * (1 + c**2) + b * s**2)
    Rxx = v**2 * (-F * (b * c**2 + a * u * s**2) + a**2 * (c**2 - s**2 / 2))
    R0r = -(v / 2) * (b - a * u) * c * s
    Rrr = sp.Rational(1, 2) * v**2 * a**2 * s**2
    Rxr = v**2 * c * s * (a**2 - F * (b - a * u) / 2)
    return sp.Matrix([[R00, R0x, R0r], [R0x, Rxx, Rxr], [R0r, Rxr, Rrr]])


def algebra():
    """Block R3, mixed M, (c2, c1, c0), E1, E2 reduced mod s^2 = 1 - c^2. Cached unless a control changes the formulas."""
    key = MUT if MUT in ("E2-typo", "E1-typo") else ""
    if key in _AL:
        return _AL[key]
    R3, _ = L.frame_block()
    M = L.mixed(R3)
    c2 = L.red(-M.trace())
    c1 = L.red(sum(M[i, i] * M[j, j] - M[i, j] * M[j, i] for i, j in [(0, 1), (0, 2), (1, 2)]))
    c0 = L.red(-M.det())
    E1 = L.red(c2**2 - 4 * c1) if MUT == "E1-typo" else L.red(c2**2 - 3 * c1)
    E2 = L.red(27 * c0 + c2**3) if MUT == "E2-typo" else L.red(27 * c0 - c2**3)
    _AL[key] = dict(R3=R3, M=M, c2=c2, c1=c1, c0=c0, E1=E1, E2=E2)
    return _AL[key]


def eff(e):
    """Rewrite a polynomial of (v, F, a, b, u, c) in the effective variables phi = vF, A = vA, B = vb: v drops out."""
    A_, B_ = sp.symbols("A B")
    e = sp.expand(e.subs({F: ph / v, a: A_ / v, b: B_ / v}))
    return e, A_, B_


_EFF = {}


def effective():
    """Entries and invariants in (X = A^2, gamma = A u, y = B - gamma, phi, w = c^2): R is linear in (X, gamma, y)."""
    if "e" in _EFF:
        return _EFF["e"]
    al = algebra()
    A_, B_ = sp.symbols("A B")
    sub = {F: ph / v, a: A_ / v, b: B_ / v}
    Rg = al["R3"].subs(sub).applyfunc(sp.expand)

    def toXgy(e):
        e = sp.expand(e.subs(B_, yy + ga))
        out = 0
        for (i, j), cf in sp.Poly(e, A_, u).terms():
            if (i, j) == (2, 0):
                out += cf * X
            elif (i, j) == (1, 1):
                out += cf * ga
            elif (i, j) == (0, 0):
                out += cf
            else:
                raise ValueError("entry is not linear in (X, gamma, y)")
        return sp.expand(out)

    Rl = Rg.applyfunc(toXgy)
    _EFF["e"] = (Rl, A_, B_)
    return _EFF["e"]


def inv_c(Mx):
    """c2, c1, c0 of a mixed 3x3 sympy matrix, unreduced."""
    c2_ = -Mx.trace()
    c1_ = sum(Mx[i, i] * Mx[j, j] - Mx[i, j] * Mx[j, i] for i, j in [(0, 1), (0, 2), (1, 2)])
    c0_ = -Mx.det()
    return sp.expand(c2_), sp.expand(c1_), sp.expand(c0_)


# ----------------------------------------------------------------------------------------------------------------------
# planted types
# ----------------------------------------------------------------------------------------------------------------------
def lor():
    """A rational Lorentz transformation of 2+1 Minkowski space (boost, rotation, boost)."""
    g1, gb1 = Fr(5, 4), Fr(3, 4)
    g2, gb2 = Fr(13, 12), Fr(5, 12)
    Lb = [[g1, gb1, 0], [gb1, g1, 0], [0, 0, 1]]
    Lr = [[1, 0, 0], [0, Fr(3, 5), Fr(-4, 5)], [0, Fr(4, 5), Fr(3, 5)]]
    Lc = [[g2, 0, gb2], [0, 1, 0], [gb2, 0, g2]]
    Lm = lambda A_, B_: [[sum(Fr(A_[i][k]) * Fr(B_[k][j]) for k in range(3)) for j in range(3)] for i in range(3)]  # noqa: E731
    return Lm(Lm(Lb, Lr), Lc)


def planted_tensors():
    """(name, R lower-index 3x3 as Fractions, expected type). R + lam0 eta shifts every eigenvalue by lam0, so that the
    invariants c2, c1, c0 do not all vanish."""
    l0 = Fr(3, 7)
    ell = [Fr(1), Fr(1), Fr(0)]  # null
    m = [Fr(0), Fr(0), Fr(1)]
    outer = lambda p, q: [[p[i] * q[j] for j in range(3)] for i in range(3)]  # noqa: E731
    add = lambda *Ms: [[sum(Mx[i][j] for Mx in Ms) for j in range(3)] for i in range(3)]  # noqa: E731
    sc = lambda k, Mx: [[k * Mx[i][j] for j in range(3)] for i in range(3)]  # noqa: E731
    eta = [[Fr(-1), 0, 0], [0, Fr(1), 0], [0, 0, Fr(1)]]
    sh = sc(l0, eta)
    jordan3 = add(outer(ell, m), outer(m, ell), sh)
    jordan2 = add(outer(ell, ell), sc(Fr(2, 3), outer(m, m)), sh)
    triple1 = add(outer(ell, ell), sh)
    fluid = [[Fr(3, 2), 0, 0], [0, Fr(1, 5), 0], [0, 0, Fr(1, 5)]]
    distinct = [[Fr(2), 0, 0], [0, Fr(1), 0], [0, 0, Fr(3)]]
    cplx = [[Fr(0), Fr(1), 0], [Fr(1), Fr(0), 0], [0, 0, Fr(1, 2)]]
    Lm = lor()
    LT = [[Lm[j][i] for j in range(3)] for i in range(3)]
    mm = lambda A_, B_: [[sum(A_[i][k] * B_[k][j] for k in range(3)) for j in range(3)] for i in range(3)]  # noqa: E731
    tf = lambda R_: mm(mm(LT, R_), Lm)  # noqa: E731
    base = [("Jordan 3x3", jordan3, "III"), ("Jordan 2x2 + simple root", jordan2, "II"), ("triple root, rank 1", triple1, "II"),
            ("perfect fluid", fluid, "I"), ("three distinct real roots", distinct, "I"), ("complex pair", cplx, "IV")]
    return [(nm, tf(R_), ty) for nm, R_, ty in base]


def mixed_frac(R_):
    return [[(-1 if i == 0 else 1) * Fr(R_[i][j]) for j in range(3)] for i in range(3)]


def deform(R_):
    out = [[Fr(e) for e in row] for row in R_]
    out[0][0] += Fr(1, 7)
    out[2][2] -= Fr(1, 11)
    return out


def planted_types(h):
    rows, ok = [], True
    for nm, R_, want in planted_tensors():
        M_ = mixed_frac(R_)
        R4 = sp.Matrix(4, 4, lambda i, j: sp.Rational(R_[i][j].numerator, R_[i][j].denominator) if i < 3 and j < 3 else 0)
        t_own = L.classify(M_)["type"]
        t_ht, _ = HT.he_type_exact(R4)
        t_fld = type_own(M_)["type"]
        c2_, c1_, c0_ = cubic_at(M_)
        e1_, e2_ = E12(c2_, c1_, c0_)
        trip = is0(e1_) and is0(e2_)
        want_trip = want == "III" or nm == "triple root, rank 1"
        good = t_own == want and t_ht == want and t_fld == want and trip == want_trip
        ok = ok and good
        rows.append("%s: classify %s, he_type_exact %s, own %s, E1 = E2 = 0 %s (want %s, %s)" % (
            nm, t_own, t_ht, t_fld, trip, want, want_trip))
    # a deformed copy changes type
    dtypes = []
    for nm, R_, want in planted_tensors()[:3]:
        Rd = deform(R_)
        td = L.classify(mixed_frac(Rd))["type"]
        dtypes.append((nm, want, td))
        ok = ok and td != want
    h.check("planted_types", ok, "; ".join(rows) + ". Deformed (R_00 + 1/7, R_22 - 1/11): " + ", ".join(
        "%s %s -> %s" % t for t in dtypes))


# ----------------------------------------------------------------------------------------------------------------------
# polynomials
# ----------------------------------------------------------------------------------------------------------------------
def forms():
    """E1 (a quadratic form) and E2 (a cubic form) in (X, gamma, y) = (v^2 a^2, v a u, v (b - a u)) with coefficients in phi = vF and
    w = c^2, from the effective Ricci block: the compact form of the conditions."""
    Rl, _, _ = effective()
    Mx = sp.diag(-1, 1, 1) * Rl
    c2_, c1_, c0_ = inv_c(Mx)

    def rd(e):
        e = sp.rem(sp.Poly(sp.expand(e), s), sp.Poly(s**2 - 1 + c**2, s)).as_expr()
        return sp.expand(sp.expand(e).subs(c, sp.sqrt(w)))

    c2r, c1r, c0r = rd(c2_), rd(c1_), rd(c0_)
    e1_, e2_ = E12(c2r, c1r, c0r)
    return sp.expand(e1_), sp.expand(e2_), c2r


def polynomials(h):
    al = algebra()
    R3, c2, c1, c0, E1, E2 = al["R3"], al["c2"], al["c1"], al["c0"], al["E1"], al["E2"]
    tab = typed_table()
    d_tab = [L.red(R3[i, j] - tab[i, j]) for i in range(3) for j in range(i, 3)]
    ok_tab = all(e == 0 for e in d_tab)
    # chi(x) = det(x - M) = det Sigma(x)
    chi = L.red(L.sigma(R3, x).det())
    ok_chi = L.red(chi - (x**3 + c2 * x**2 + c1 * x + c0)) == 0
    # symmetry, s-free, v-free in the effective variables, weights
    ok_c = (not E1.has(s)) and (not E2.has(s)) and sp.expand(E1.subs(c, -c) - E1) == 0 and sp.expand(E2.subs(c, -c) - E2) == 0
    e1e, A_, B_ = eff(E1)
    e2e, _, _ = eff(E2)
    ok_v = (not e1e.has(v)) and (not e2e.has(v))
    ell = sp.Symbol("ell")
    k1 = 5 if MUT == "weight-typo" else 4
    w1 = sp.expand(e1e.subs({A_: ell * A_, B_: ell**2 * B_, u: ell * u}) - ell**k1 * e1e)
    w2 = sp.expand(e2e.subs({A_: ell * A_, B_: ell**2 * B_, u: ell * u}) - ell**6 * e2e)
    ok_w = w1 == 0 and w2 == 0
    # E1 and E2 on a double root p with third root q
    p_, q_ = sp.symbols("p q")
    ee1, ee2 = E12(-(2 * p_ + q_), p_**2 + 2 * p_ * q_, -p_**2 * q_)
    ok_pq = sp.expand(ee1 - (p_ - q_)**2) == 0 and sp.expand(ee2 - (p_ - q_)**2 * (8 * p_ + q_)) == 0
    # a triple root: p = q gives E1 = E2 = 0 for all p
    ee3 = E12(-3 * p_, 3 * p_**2, -p_**3)
    ok_tr = sp.expand(ee3[0]) == 0 and sp.expand(ee3[1]) == 0
    E1x, E2x, c2x = forms()
    back = {X: v**2 * a**2, ga: v * a * u, yy: v * (b - a * u), ph: v * F, w: c**2}
    ok_forms = sp.expand(E1x.subs(back) - E1) == 0 and sp.expand(E2x.subs(back) - E2) == 0
    cx = lambda e: sp.collect(e, [X, ga, yy], sp.factor)  # noqa: E731
    n1, n2 = len(sp.Poly(E1, v, F, a, b, u, c).terms()), len(sp.Poly(E2, v, F, a, b, u, c).terms())
    f1 = sp.factor_list(E1)
    f2 = sp.factor_list(E2)
    h.check("polynomials", ok_tab and ok_chi and ok_c and ok_v and ok_w and ok_pq and ok_tr and ok_forms,
            "typed table of he_surface.md section 1 minus own Ricci block reduces to 0: %s; det(x - M) = det Sigma(x) = x^3 + c2 x^2 + "
            "c1 x + c0: %s; E1, E2 are s-free and even in c: %s; in phi = vF, A = va, B = vb they are v-free: %s, with weights 4 "
            "and 6 under (A, B, u) -> (l A, l^2 B, l u): %s; E1 = (p - q)^2 and E2 = (p - q)^2 (8p + q) on (p, p, q): %s; "
            "E1 = E2 = 0 on (p, p, p): %s; E1 has %d terms and E2 has %d terms in (v, F, a, b, u, c); factors of E1: %s, of E2: %s. In "
            "X = v^2 a^2, gamma = v a u, y = v (b - a u), phi = v F, w = c^2 (identity with the above: %s) c2 = %s, E1 = %s, E2 = %s" % (
                ok_tab, ok_chi, ok_c, ok_v, ok_w, ok_pq, ok_tr, n1, n2,
                [(sp.Poly(f, v, F, a, b, u, c).total_degree(), k) for f, k in f1[1]],
                [(sp.Poly(f, v, F, a, b, u, c).total_degree(), k) for f, k in f2[1]], ok_forms, sp.factor(c2x), cx(E1x), cx(E2x)))


# ----------------------------------------------------------------------------------------------------------------------
# certificates
# ----------------------------------------------------------------------------------------------------------------------
def lamstar():
    return sp.expand(-sp.Rational(1, 2) * v**2 * a * s**2 * (a + 2 * F * u))


def Theta():
    cc = (1 + c) if MUT == "theta-typo" else (1 + c**2)
    return sp.expand(v**2 * a * (2 * a + cc * F * u))


def Pu():
    return sp.expand(v**2 * (a + F * u) * (a + c**2 * F * u) + c**2 * u**2)


def Kq():
    return sp.expand(v**2 * a**2 * Pu())


def divT(f):
    """f = q (b + a u) + r, r = f(b = -a u): exact division of a polynomial by T in the variable b."""
    q, r = sp.div(sp.Poly(sp.expand(f), b), sp.Poly(b + a * u, b))
    return q.as_expr(), r.as_expr()


def minors_sigma():
    al = algebra()
    S = L.sigma(al["R3"])
    return [L.red(m) for m in L.minors(S)]


def cert_rank(h):
    """C1 of #101, restated: K1 T = sum q_k m_k for every lam, where the m_k are the six 2x2 minors of Sigma(lam)."""
    mins = minors_sigma()
    K1 = sp.Rational(1, 2) * v**2 * a * u * c * s * (1 + v**2 * F**2)
    q = [-c * s, c**2 - s**2, v * F * c**2, c * s, v * F * c * s, 0]
    if MUT == "cof-rank-perturb":
        q[0] = q[0] + 1
    if MUT == "cof-rank-drop":
        q[4] = 0
    lhs = L.red(sum(qk * mk for qk, mk in zip(q, mins)))
    d1 = L.red(lhs - K1 * Tpol)
    S = L.sigma(algebra()["R3"])
    d6 = L.red(mins[0].subs(a, 0) + S[0, 0].subs(a, 0)**2 + S[0, 1].subs(a, 0)**2)
    s01 = L.red(S[0, 1].subs(a, 0) - v * b * s**2 / 2)
    h.check("cert_rank", d1 == 0 and d6 == 0 and s01 == 0,
            "K1 (b + a u) = -cs m01,01 + (c^2 - s^2) m01,02 + vF c^2 m01,12 + cs m02,02 + vF cs m02,12, K1 = v^2 a u c s (1 + v^2 F^2)/2, "
            "for symbolic lam: residual %s; at a = 0, m01,01 + Sigma_00^2 + Sigma_01^2 = %s and Sigma_01 = v b s^2/2: %s. So for real lam and "
            "v c s u != 0, (a, b) != (0, 0): all six minors zero (rank Sigma <= 1) gives a != 0 (else Sigma_00 = Sigma_01 = 0, b = 0), "
            "then K1 != 0 and T = 0" % (d1, d6, s01 == 0))


def chain():
    """All exact quantities of the T = 0 certificate, as sympy expressions reduced modulo c^2 + s^2 = 1."""
    al = algebra()
    R3, E1, E2, c2 = al["R3"], al["E1"], al["E2"], al["c2"]
    ls, Th, K_ = lamstar(), Theta(), Kq()
    chi = L.red(L.sigma(R3, x).det())
    chiT = L.red(chi.subs(b, -a * u))
    target = sp.expand((x - ls) * ((x - ls) * (x - ls - Th) + K_))
    H, rH = divT(L.red(chi - chiT))
    E1T = sp.expand(Th**2 - 3 * K_)
    E2T = sp.expand(Th**3 + 9 * ls * (Th**2 - 3 * K_))
    g1, r1 = divT(L.red(E1 - E1T))
    g2, r2 = divT(L.red(E2 - E2T))
    mins = minors_sigma()
    Gam = [c**2, c * s, 0, s**2, 0, 0]
    A_, B_, r_ = [], [], []
    for mk in mins:
        Ak, rk = divT(mk)
        Bk, remk = sp.div(sp.Poly(L.red(rk), lam), sp.Poly(lam - ls, lam))
        A_.append(Ak)
        B_.append(Bk.as_expr())
        r_.append(remk.as_expr())
    h2 = sp.Poly(H, x).coeff_monomial(x**2)
    return dict(ls=ls, Th=Th, K=K_, chi=chi, chiT=chiT, target=target, H=H, rH=rH, E1T=E1T, E2T=E2T, g1=g1, r1=r1, g2=g2,
                r2=r2, mins=mins, Gam=Gam, A=A_, B=B_, rem=r_, h2=h2, c2=c2, E1=E1, E2=E2)


def cert_T0(h):
    """On T = 0 the characteristic polynomial is (x - l*)[(x - l*)(x - l* - Theta) + K], E1 = Theta^2 - 3K and
    E2 = Theta^3 + 9 l* E1 up to multiples of T, and the six minors at lam0 = tr M / 3 lie in <T, Theta, K>."""
    C = chain()
    ls, Th, K_, g1, g2, h2 = C["ls"], C["Th"], C["K"], C["g1"], C["g2"], C["h2"]
    E1, E2, c2 = C["E1"], C["E2"], C["c2"]
    if MUT == "cof-g1-perturb":
        g1 = g1 + 1
    if MUT == "cof-g2-drop":
        g2 = 0 * g2
    z1 = L.red(C["chiT"] - C["target"])
    z1g = L.red(C["chi"] - C["target"] - Tpol * C["H"])
    z3 = L.red(E1 - C["E1T"] - Tpol * g1)
    z4 = L.red(E2 - C["E2T"] - Tpol * g2)
    z6 = L.red(c2 + 3 * ls + Th - Tpol * h2)
    rem_ok = all(L.red(e) == 0 for e in (C["rH"], C["r1"], C["r2"]))
    nulls = L.red(Th**3 - (E2 - 9 * ls * E1 + Tpol * (9 * ls * g1 - g2)))
    nullk = L.red(3 * K_ - (Th**2 - E1 + Tpol * g1))
    z5, z5t = [], []
    for k_ in range(6):
        Bk = C["B"][k_]
        if MUT == "cof-B-perturb" and k_ == 0:
            Bk = Bk + 1
        resid = L.red(C["mins"][k_] - (Tpol * C["A"][k_] + (lam - ls) * Bk - K_ * C["Gam"][k_]))
        z5.append(resid == 0)
        z5t.append(L.red(C["rem"][k_] + K_ * C["Gam"][k_]) == 0)
    l0 = -c2 / 3
    z7 = []
    for k_ in range(6):
        Bk = C["B"][k_]
        if MUT == "cof-B-perturb" and k_ == 0:
            Bk = Bk + 1
        rhs = (Tpol * (C["A"][k_].subs(lam, l0) - h2 * Bk.subs(lam, l0) / 3) + Th / 3 * Bk.subs(lam, l0) - K_ * C["Gam"][k_])
        z7.append(L.red(C["mins"][k_].subs(lam, l0) - rhs) == 0)
    ok = z1 == 0 and z1g == 0 and z3 == 0 and z4 == 0 and z6 == 0 and rem_ok and nulls == 0 and nullk == 0 and all(z5) and all(z5t) and all(z7)
    h.check("cert_T0", ok,
            "chi(x) at T = 0 minus (x - l*)[(x - l*)(x - l* - Theta) + K]: %s; chi - target = T H: %s; E1 - (Theta^2 - 3K) = T g1: %s "
            "(g1 has %d terms); E2 - (Theta^3 + 9 l* (Theta^2 - 3K)) = T g2: %s (g2 has %d terms); c2 + 3 l* + Theta = T h2: %s; "
            "Theta^3 = E2 - 9 l* E1 + T (9 l* g1 - g2): %s; 3K = Theta^2 - E1 + T g1: %s; the six minors m_k(lam) = T A_k + (lam - l*) B_k "
            "- K Gamma_k, Gamma = (c^2, cs, 0, s^2, 0, 0): %s, and their remainders at l* are -K Gamma_k: %s; at lam0 = -c2/3, "
            "m_k(lam0) = T (A_k - h2 B_k / 3) + (Theta / 3) B_k - K Gamma_k: %s. Hence at T = 0 and E1 = E2 = 0: Theta^3 = 0, "
            "K = 0, lam0 = l*, every m_k(lam0) = 0, rank(M - lam0) <= 1. l* = -v^2 a s^2 (a + 2Fu)/2, Theta = v^2 a (2a + (1 + c^2) F u), "
            "K = v^2 a^2 Pu, Pu = v^2 (a + Fu)(a + c^2 F u) + c^2 u^2" % (
                z1, z1g == 0, z3 == 0, len(sp.Poly(g1, v, F, a, u, c, s).terms()), z4 == 0, len(sp.Poly(g2, v, F, a, u, c, s).terms()),
                z6 == 0, nulls == 0, nullk == 0, all(z5), all(z5t), all(z7)))


def cert_disc(h):
    """27 D + E2^2 = 27 E1 G in Q[c2, c1, c0]: on E1 = 0, D = -E2^2/27, so D = 0 forces E2 = 0."""
    C2, C1, C0 = sp.symbols("C2 C1 C0")
    D = sp.expand(disc3(C2, C1, C0))
    e1_, e2_ = E12(C2, C1, C0)
    k = 26 if MUT == "disc27" else 27
    num = sp.expand(k * D + e2_**2)
    Gq = sp.cancel(num / (27 * e1_))
    num_, den_ = sp.fraction(Gq)
    is_poly = den_.is_number
    G = sp.expand(Gq) if is_poly else Gq
    if MUT == "cof-disc-perturb":
        G = G + C2
    res = sp.expand(k * D + e2_**2 - 27 * e1_ * G)
    at_e1 = sp.expand(D.subs(C1, C2**2 / 3) + E12(C2, C2**2 / 3, C0)[1]**2 / 27)
    h.check("cert_disc", is_poly and res == 0 and at_e1 == 0,
            "27 D + E2^2 - 27 E1 G = %s with G = %s (polynomial: %s); on E1 = 0, D + E2^2/27 = %s. So E1 = 0 and D = 0 give E2 = 0, "
            "and E1 = E2 = 0 give D = 0" % (res, sp.factor(G) if is_poly else "n/a", is_poly, at_e1))


def equator_axis(h):
    """R_0rho and R_x,rho carry the factor c s, so at c = 0 and at s = 0 the block M is (n, e_x) + (e_rho): no 3x3 Jordan block."""
    M = algebra()["M"]
    if MUT == "axis-typo":
        pts = [{c: sp.Rational(1, 2), s: sp.sqrt(3) / 2}, {c: sp.Rational(1, 3), s: 2 * sp.sqrt(2) / 3}]
    else:
        pts = [{c: 0, s: 1}, {c: 1, s: 0}]
    zero = all(sp.expand(M[i, j].subs(p_)) == 0 for p_ in pts for (i, j) in [(0, 2), (1, 2), (2, 0), (2, 1)])
    p_, q_, r_, t_ = sp.symbols("p q r t")
    N2 = sp.Matrix([[p_, q_], [r_, t_]])
    ch = (N2 * N2 - N2.trace() * N2 + N2.det() * sp.eye(2)).applyfunc(sp.expand) == sp.zeros(2, 2)
    h.check("equator_axis", zero and ch,
            "M[0,2], M[1,2], M[2,0], M[2,1] vanish at c = 0 (equator) and at s = 0 (axis): %s. Then M = (2x2 block) + (1x1), a nilpotent "
            "2x2 block N2 has N2^2 = tr(N2) N2 - det(N2) I = 0 (Cayley-Hamilton, identity %s), so N^2 = 0 and rank <= 1 at a triple root: "
            "no Type III point on the equator or the axis" % (zero, ch))


# ----------------------------------------------------------------------------------------------------------------------
# the null ray: a Type III point has a null kernel vector, which gives its exact points
# ----------------------------------------------------------------------------------------------------------------------
def red_ce(ex):
    """Reduce modulo e2^2 = 1 - e1^2 and s^2 = 1 - c^2."""
    ex = sp.expand(ex)
    ex = sp.rem(sp.Poly(ex, e2s), sp.Poly(e2s**2 - 1 + e1s**2, e2s)).as_expr()
    ex = sp.rem(sp.Poly(sp.expand(ex), s), sp.Poly(s**2 - 1 + c**2, s)).as_expr()
    return sp.expand(ex)


_RAY = {}


def ray_symbolic():
    """Cofactors n = (lam, X, gamma, y) of the 3x4 matrix of Sigma(lam) e = 0, e = (1, e1, e2), in the effective variables,
    and the polynomials tr M - kappa lam on the ray as polynomials in phi."""
    key = "r" + ("-typo" if MUT == "ray-typo" else "")
    if key in _RAY:
        return _RAY[key]
    Rl, _, _ = effective()
    eta = sp.diag(-1, 1, 1)
    e = sp.Matrix([1, e1s, e2s])
    eqs = ((Rl - lam * eta) * e).applyfunc(sp.expand)
    unk = [lam, X, ga, yy]
    K = sp.Matrix([[sp.diff(q, z) for z in unk] for q in eqs])
    n = [sp.expand((-1)**i * K[:, [j for j in range(4) if j != i]].det()) for i in range(4)]
    if MUT == "ray-typo":
        n[1] = n[1] + 1
    sub = {lam: n[0], X: n[1], ga: n[2], yy: n[3]}
    trM = (-Rl[0, 0] + Rl[1, 1] + Rl[2, 2]).subs(sub)
    kn = sp.expand(K * sp.Matrix(n))
    taus = {}
    for kap in (3, -6):
        t_ = red_ce(trM - kap * n[0])  # reduced modulo the two circle relations: valid at every exact point used here
        taus[kap] = [PolyEval(cf, (c, s, e1s, e2s)) for cf in sp.Poly(t_, ph).all_coeffs()]
    syms = (c, s, e1s, e2s, ph)
    _RAY[key] = dict(Rl=Rl, K=K, n=n, kn=kn, trM=trM, nev=[PolyEval(red_ce(q), syms) for q in n], taus=taus)
    return _RAY[key]


def ray_family(h):
    """On the null ray, E1 = tau^2 and E2 = tau^2 (9 lam + tau), tau = tr M - 3 lam, symbolically modulo c^2 + s^2 = 1 and
    e1^2 + e2^2 = 1; tau is a quadratic in phi with leading coefficient c s^3 / 2."""
    R = ray_symbolic()
    zero_kn = all(z == 0 for z in R["kn"])
    nr = [red_ce(q) for q in R["n"]]
    sub = {lam: nr[0], X: nr[1], ga: nr[2], yy: nr[3]}
    M = (sp.diag(-1, 1, 1) * R["Rl"]).applyfunc(lambda z: red_ce(z.subs(sub)))
    tau = red_ce(M.trace() - 3 * nr[0])
    c2_, c1_, c0_ = (red_ce(q) for q in inv_c(M))
    E1_, E2_ = E12(c2_, c1_, c0_)
    d1 = red_ce(E1_ - tau**2)
    d2 = red_ce(E2_ - tau**2 * (9 * nr[0] + tau))
    q_ = M.trace() - 2 * nr[0]  # the third root: chi(x) = (x - lam)^2 (x - q)
    dD = [red_ce(c2_ + 2 * nr[0] + q_), red_ce(c1_ - (nr[0]**2 + 2 * nr[0] * q_)), red_ce(c0_ + nr[0]**2 * q_)]
    dD = 0 if all(z == 0 for z in dD) else dD
    P = sp.Poly(tau, ph)
    lead = sp.factor(P.all_coeffs()[0])
    h.check("ray_family", zero_kn and d1 == 0 and d2 == 0 and dD == 0 and P.degree() == 2 and L.red(lead - c * s**3 / 2) == 0,
            "the 3x4 system of Sigma(lam) e = 0, e = (1, e1, e2), is solved by the cofactor vector n(c, s, e1, e2, phi) (K n = 0: %s); on it "
            "E1 - tau^2 = %s, E2 - tau^2 (9 lam + tau) = %s and chi(x) - (x - lam)^2 (x - q) = %s (q = tr M - 2 lam, so D = 0) modulo the two circle relations, so E1 = E2 = 0 exactly when "
            "tau = tr M - 3 lam = 0; tau is a quadratic in phi with leading coefficient %s" % (zero_kn, d1, d2, dD, lead))


# ----------------------------------------------------------------------------------------------------------------------
# exact points in Q and Q(sqrt d)
# ----------------------------------------------------------------------------------------------------------------------
_BLK = {}


def block_evals():
    if "b" not in _BLK:
        R3 = algebra()["R3"]
        _BLK["b"] = [[PolyEval(R3[i, j], VARS7) for j in range(3)] for i in range(3)]
        _BLK["E"] = (PolyEval(algebra()["E1"], VARS7), PolyEval(algebra()["E2"], VARS7))
        Jv = (v, F, a, b, u, c)
        al = algebra()
        _BLK["J"] = [[PolyEval(sp.diff(al[k], z), VARS7) for z in Jv] for k in ("E1", "E2")]
    return _BLK


def block_at(pt):
    B = block_evals()["b"]
    return [[B[i][j](pt) for j in range(3)] for i in range(3)]


def mixed_at(R_):
    return [[(-1 if i == 0 else 1) * R_[i][j] for j in range(3)] for i in range(3)]


def pyth_pairs(maxm):
    out = []
    for m in range(2, maxm):
        for n in range(1, m):
            if math.gcd(m, n) == 1 and (m - n) % 2 == 1:
                p_, q_, r_ = m * m - n * n, 2 * m * n, m * m + n * n
                out += [(Fr(p_, r_), Fr(q_, r_)), (Fr(q_, r_), Fr(p_, r_))]
    return out


def sqfree(n):
    """n = d k^2 with d squarefree, for a positive integer n."""
    d_, k_ = 1, 1
    for pr, ex in sp.factorint(n).items():
        d_ *= pr ** (ex % 2)
        k_ *= pr ** (ex // 2)
    return d_, k_


def quad_roots(al_, be_, ga_):
    """Roots of al x^2 + be x + ga = 0 over Q: list of ('Q', root) or ('S', d, p, q) meaning p + q sqrt d, d squarefree > 1."""
    if al_ == 0:
        return [("Q", -ga_ / be_)] if be_ != 0 else []
    Dl = be_**2 - 4 * al_ * ga_
    if Dl < 0:
        return []
    d_, k_ = sqfree(Dl.numerator * Dl.denominator)
    if d_ == 1:
        r0 = Fr(k_, Dl.denominator)
        return [("Q", (-be_ + sg * r0) / (2 * al_)) for sg in (1, -1)]
    r0 = Fr(k_, Dl.denominator)
    return [("S", d_, -be_ / (2 * al_), sg * r0 / (2 * al_)) for sg in (1, -1)]


def ray_data(c_, s_, e1_, e2_, kappa):
    """Roots (phi) of tau_kappa(phi) = 0 at the exact point (c, s, e1, e2) of the null ray."""
    R = ray_symbolic()
    co = [f((c_, s_, e1_, e2_)) for f in R["taus"][kappa]]
    while len(co) > 3 and co[0] == 0:
        co = co[1:]
    assert len(co) <= 3, "tau is not quadratic in phi"
    co = [Fr(0)] * (3 - len(co)) + co
    return quad_roots(*co)


def make_point(c_, s_, e1_, e2_, root, vv):
    """The point (v, F, a, b, u, c, s) of the null ray at phi = root, scaled to X = A^2 = 1, with u > 0, as elements of Q or
    Q(sqrt d) (QS.d is set here). Returns (point, lam, phi, field)."""
    R = ray_symbolic()
    if root[0] == "Q":
        phi_, field = Fr(root[1]), None
    else:
        QS.d = root[1] + (1 if MUT == "nf-radicand" else 0)
        phi_, field = QS(root[2], root[3]), root[1]
    if MUT == "phi-off":
        phi_ = phi_ + Fr(1, 1000)
    nv = [f((c_, s_, e1_, e2_, phi_)) for f in R["nev"]]
    if is0(nv[1]):
        return None
    inv = 1 / nv[1]
    lam1, ga1, y1 = nv[0] * inv, nv[2] * inv, nv[3] * inv
    if is0(ga1):
        return None
    A_ = 1 if float(ga1) > 0 else -1
    u_ = ga1 * A_
    B_ = y1 + ga1
    pt = (Fr(vv), phi_ / vv, Fr(A_) / vv, B_ / vv, u_, c_, s_)
    return pt, lam1, phi_, field


def point_invariants(pt):
    """Block, mixed matrix, cubic data at an exact point: dict from `type_own` plus T and the block."""
    R_ = block_at(pt)
    M_ = mixed_at(R_)
    out = type_own(M_)
    out["T"] = pt[3] + pt[2] * pt[4]
    out["R"], out["M"] = R_, M_
    return out


def to_mp(z, dps=60):
    mp.mp.dps = dps
    if isinstance(z, (mp.mpf, mp.mpc)):
        return z
    if isinstance(z, QS):
        return mp.mpf(z.p.numerator) / z.p.denominator + mp.mpf(z.q.numerator) / z.q.denominator * mp.sqrt(QS.d)
    z = Fr(z)
    return mp.mpf(z.numerator) / z.denominator


def numeric_type3(M_, dps=60):
    """Singular values of N = M - lam0 and the size of N^3 at dps digits, lam0 = tr M / 3."""
    mp.mp.dps = dps
    Mm = mp.matrix([[to_mp(M_[i][j], dps) for j in range(3)] for i in range(3)])
    l0 = (Mm[0, 0] + Mm[1, 1] + Mm[2, 2]) / 3
    N = Mm - l0 * mp.eye(3)
    sv = mp.svd_r(N, compute_uv=False)
    N3 = N * N * N
    return sorted([abs(z) for z in sv], reverse=True), max(abs(N3[i, j]) for i in range(3) for j in range(3))


def search_seeds(maxm=7):
    """The search that produced TYPE3_SEEDS (`--search-seeds`): over Pythagorean (c, s) and null directions e = (1, e1, e2) it keeps,
    for each squarefree d, the seed with the smallest denominators whose root phi gives a point with F > 0, u > 0, T != 0 and E1 = E2 = 0
    exactly, and counts the seeds whose phi is rational. Returns (seeds, number of seeds tried, rational points found)."""
    ps = pyth_pairs(maxm)
    sg = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    cand, tried, nrat = {}, 0, 0
    for (cc, ss) in ps:
        for (p1, p2) in ps:
            for g1, g2 in sg:
                sd = (cc, ss, p1 * g1, p2 * g2)
                tried += 1
                for root in ray_data(*sd, 3):
                    res = make_point(*sd, root, 2)
                    if res is None:
                        continue
                    pt = res[0]
                    if float(pt[1]) <= 0 or float(pt[4]) <= 0:
                        continue
                    inv = point_invariants(pt)
                    if is0(inv["T"]) or not (is0(inv["E1"]) and is0(inv["E2"])):
                        continue
                    if root[0] == "Q":
                        nrat += 1
                        continue
                    key = max(cc.denominator, sd[2].denominator)
                    if root[1] not in cand or key < cand[root[1]][0]:
                        cand[root[1]] = (key, sd, root)
    return [cand[d_][1:] for d_ in sorted(cand)], tried, nrat


TYPE3_SEEDS = [
    ((Fr(3, 5), Fr(4, 5), Fr(-7, 25), Fr(24, 25)), ("S", 113, Fr(-9, 16), Fr(1, 16))),  # #72's point
    ((Fr(4, 5), Fr(3, 5), Fr(7, 25), Fr(24, 25)), ("S", 463, Fr(-16, 9), Fr(1, 9))),
    ((Fr(3, 5), Fr(4, 5), Fr(24, 25), Fr(7, 25)), ("S", 949, Fr(-11, 24), Fr(1, 24))),
    ((Fr(4, 5), Fr(3, 5), Fr(24, 25), Fr(-7, 25)), ("S", 1649, Fr(-13, 8), Fr(1, 24))),
    ((Fr(7, 25), Fr(24, 25), Fr(3, 5), Fr(-4, 5)), ("S", 7793, Fr(-65, 448), Fr(5, 1344))),
    ((Fr(24, 25), Fr(7, 25), Fr(3, 5), Fr(4, 5)), ("S", 120367, Fr(-1660, 147), Fr(5, 147))),
]


def ht_type(R_):
    """he_type_exact (#64, #103) on the exact 3x3 block R_ of Fractions, decoupled axis added."""
    R4 = sp.Matrix(4, 4, lambda i, j: sp.Rational(R_[i][j].numerator, R_[i][j].denominator) if i < 3 and j < 3 else 0)
    return HT.he_type_exact(R4)[0]


def fmt_num(z, dps=8):
    return mp.nstr(to_mp(z, 30), dps)


def exact_points(h):
    """Type III points of the null ray in Q(sqrt d), #72's in Q(sqrt 113) first: every exact test, then 60 digits."""
    rows, ok = [], True
    seeds = TYPE3_SEEDS
    for idx, (sd, root) in enumerate(seeds):
        ok = ok and (root in ray_data(*sd, 3))
        res = make_point(*sd, root, 2)
        if res is None:
            ok = False
            rows.append("seed %d gives no point" % idx)
            continue
        pt, lam1, phi_, field = res
        inv = point_invariants(pt)
        R_, M_ = inv["R"], inv["M"]
        e_ = (Fr(1), sd[2], sd[3])
        eta_ = (-1, 1, 1)
        ker = all(is0(sum((R_[i][j] * e_[j] for j in range(3)), 0 * R_[0][0]) - lam1 * eta_[i] * e_[i]) for i in range(3))
        Jm = [[f(pt) for f in row] for row in block_evals()["J"]]
        rJ = rank_any(Jm)
        sv, n3 = numeric_type3(M_)
        good = (ker and inv["type"] == "III" and inv["mult"] == 3 and inv["rank"] == 2 and is0(inv["D"]) and is0(inv["E1"])
                and is0(inv["E2"]) and not is0(inv["T"]) and not inv["N2zero"] and inv["N3zero"] and inv["lam"] == lam1
                and rJ == 2 and sv[2] < mp.mpf(10)**-50 and sv[1] > mp.mpf(10)**-6 and n3 < mp.mpf(10)**-50
                and float(pt[1]) > 0 and float(pt[4]) > 0)
        ok = ok and good
        fld = "Q(sqrt %d)" % field if field else "Q"
        extra = ""
        if idx == 0:
            ar, br2 = float(pt[2] / pt[4]), float(pt[3] / pt[4]**2)
            # #72's printed values a = -0.2394, b = -0.1307, r = 7.25 and F = 0.0509 (he_surface.md section 3)
            same = abs(ar - (-0.2394 * 7.25)) < 4e-3 and abs(br2 - (-0.1307 * 7.25**2)) < 2e-2 and abs(float(pt[1]) - 0.0509) < 5e-5
            ok = ok and same
            co = [f((sd[0], sd[1], sd[2], sd[3])) for f in ray_symbolic()["taus"][3]]
            co = [Fr(0)] * (3 - len(co)) + co
            quad = "tau(phi) = (%s)(8 phi^2 + 9 phi - 1)" % (co[0] / 8) \
                if co[1] / co[0] == Fr(9, 8) and co[2] / co[0] == Fr(-1, 8) else "tau(phi) = %s phi^2 + %s phi + %s (NOT #72's quadratic)" % tuple(co)
            ok = ok and co[1] / co[0] == Fr(9, 8) and co[2] / co[0] == Fr(-1, 8)
            extra = " [#72's point: same c, e, phi, v, and %s, the quadratic of he_surface.md section 3; a r = %.4f, b r^2 = %.3f against -1.7357, -6.870 " \
                    "from the printed a, b, r; F = %.5f against 0.0509: %s]" % (quad, ar, br2, float(pt[1]), same)
        rows.append("%s c = %s, e = (1, %s, %s), phi = %s, v = 2: F = %s, a = %s, b = %s, u = %s; T = %s != 0, D = 0, E1 = E2 = 0, triple root "
                    "lam = %s, rank(M - lam) = %d, N^2 != 0, N^3 = 0, Jacobian of (E1, E2) in (v, F, a, b, u, c) has rank %d, singular values "
                    "%s, %s, %s, |N^3| = %s%s" % (
                        fld, sd[0], sd[2], sd[3], mp.nstr(to_mp(phi_, 30), 12), fmt_num(pt[1]), fmt_num(pt[2]), fmt_num(pt[3]),
                        fmt_num(pt[4]), fmt_num(inv["T"]), fmt_num(inv["lam"]), inv["rank"], rJ, mp.nstr(sv[0], 5), mp.nstr(sv[1], 5),
                        mp.nstr(sv[2], 3), mp.nstr(n3, 3), extra))
    # he_type_exact on #72's point in Q(sqrt 113): the sympy matrix with a radical
    t_ht = "not run"
    if seeds:
        sd, root = seeds[0]
        pt, lam1, phi_, field = make_point(*sd, root, 2)
        R_ = block_at(pt)
        sq = sp.sqrt(QS.d)
        conv = lambda z: sp.Rational(z.p.numerator, z.p.denominator) + sp.Rational(z.q.numerator, z.q.denominator) * sq if isinstance(z, QS) else sp.Rational(z.numerator, z.denominator)  # noqa: E731
        R4 = sp.Matrix(4, 4, lambda i, j: conv(R_[i][j]) if i < 3 and j < 3 else 0)
        t_ht, D_ht = HT.he_type_exact(R4)
        ok = ok and t_ht == "III"
    h.check("exact_points", ok and len(seeds) >= 5,
            "%d Type III points of the null ray. %s. he_type_exact (#64, #103) on #72's point in Q(sqrt 113), the matrix with radicals: Type %s" % (
                len(seeds), ". ".join(rows), t_ht))


# ----------------------------------------------------------------------------------------------------------------------
# the description, Type II points failing one condition, and points of every type
# ----------------------------------------------------------------------------------------------------------------------
def in_description(sm):
    """E1 = 0, E2 = 0 and T != 0, from the summary of an exact point; the controls loosen one condition."""
    e1ok, e2ok, tok = sm["E1z"], sm["E2z"], not sm["Tz"]
    if MUT == "loosen-E1":
        e1ok = True
    if MUT == "loosen-E2":
        e2ok = True
    if MUT == "loosen-T":
        tok = True
    return e1ok and e2ok and tok


def summ(pt, with_ht=True):
    inv = point_invariants(pt)
    out = dict(type=inv["type"], D0=is0(inv["D"]), E1z=is0(inv["E1"]), E2z=is0(inv["E2"]), Tz=is0(inv["T"]), rank=inv["rank"],
               mult=inv["mult"], Dsign=(None if is0(inv["D"]) or any(isinstance(z, QS) for z in (inv["D"],)) else (1 if inv["D"] > 0 else -1)))
    if with_ht and not any(isinstance(z, QS) for row in inv["R"] for z in row):
        out["ht"] = ht_type(inv["R"])
        out["cl"] = L.classify(inv["M"])["type"]
    return out


def triple_rank1_points():
    """T1, T2 and not T3: |vF| s^2 = 2|c|, r a = -(1 + c^2) F / 2, b = -a u. A triple root with rank 1."""
    out = []
    for cc, ss in [(Fr(3, 5), Fr(4, 5)), (Fr(5, 13), Fr(12, 13)), (Fr(8, 17), Fr(15, 17))]:
        phi_ = 2 * cc / ss**2
        vv, u_ = Fr(3, 2), Fr(1)
        A_ = -(1 + cc**2) * phi_ / 2 * u_
        out.append(("triple root, rank 1, c = %s" % cc, (vv, phi_ / vv, A_ / vv, -(A_ / vv) * u_, u_, cc, ss)))
    return out


MINUS8_SEEDS = [(Fr(7, 25), Fr(24, 25), Fr(-35, 37), Fr(-12, 37)), (Fr(12, 37), Fr(35, 37), Fr(-24, 25), Fr(-7, 25)),
                (Fr(3, 5), Fr(4, 5), Fr(-7, 25), Fr(24, 25))]


def minus8_points():
    """A double root p with rank 2 and third root -8p: tr M = -6 lam on the null ray, so E2 = 0 and E1 = 81 p^2 != 0. The first two
    seeds give rational points, the third a point in Q(sqrt d)."""
    out = []
    for sd in MINUS8_SEEDS:
        for root in ray_data(*sd, -6):
            res = make_point(*sd, root, 2)
            if res is None:
                continue
            pt, lam1, phi_, field = res
            if float(pt[1]) > 0 and float(pt[4]) > 0 and not is0(point_invariants(pt)["T"]):
                out.append(("double root p, third root -8p, e = (1, %s, %s), %s" % (sd[2], sd[3],
                            ("Q(sqrt %d)" % field) if field else "Q"), pt, field, root, sd))
                break
    return out


def generic_points():
    sd = (Fr(3, 5), Fr(4, 5), Fr(-7, 25), Fr(24, 25))
    out = []
    for ph_ in (Fr(1, 3), Fr(2)):
        pt, lam1, phi_, field = make_point(*sd, ("Q", ph_), 2)
        out.append(("generic double root, phi = %s" % ph_, pt))
    return out


def typeI_points():
    """The Type I points of #72's family, placed here: v F s^2 = c (k + 1/k) and r a = [-(1 + c^2) F + e (c/v)(k - 1/k)]/2."""
    out = []
    for (cc, ss), kk, eps, r_, F_ in [((Fr(3, 5), Fr(4, 5)), Fr(2), 1, Fr(3, 2), Fr(1, 2)), ((Fr(5, 13), Fr(12, 13)), Fr(3), -1, Fr(5, 4), Fr(2, 3)),
                                      ((Fr(8, 17), Fr(15, 17)), Fr(5, 2), 1, Fr(2), Fr(3, 4))]:
        vv = cc * (kk + 1 / kk) / (F_ * ss**2)
        ra = (-(1 + cc**2) * F_ + eps * (cc / vv) * (kk - 1 / kk)) / 2
        a_ = ra / r_
        u_ = 1 / r_
        T2 = vv**2 * (ra + F_) * (ra + cc**2 * F_) + cc**2
        assert T2 == 0
        out.append(("Type I, T1 T2 T3, c = %s, k = %s" % (cc, kk), (vv, F_, a_, -a_ * u_, u_, cc, ss)))
    return out


def random_points():
    return [("random point %d" % i, p_) for i, p_ in enumerate([
        (Fr(1, 2), Fr(1, 2), Fr(1, 3), Fr(1, 5), Fr(2, 3), Fr(3, 5), Fr(4, 5)),
        (Fr(3), Fr(1, 2), Fr(-1, 3), Fr(1, 5), Fr(1), Fr(5, 13), Fr(12, 13)),
        (Fr(10), Fr(3, 4), Fr(1, 2), Fr(-1, 7), Fr(2), Fr(4, 5), Fr(3, 5)),
        (Fr(1, 5), Fr(1), Fr(1, 2), Fr(1, 3), Fr(1, 2), Fr(8, 17), Fr(15, 17))])]


def e1_conic_point():
    """A point with E1 = 0 and E2 != 0, hence D = -E2^2/27 < 0 (Type IV), with T != 0. E1 is a quadratic form in (X, gamma, y) with
    coefficients in (phi, c); at c = 3/5, X = 1 and rational (phi, gamma) on a fixed grid, the first y that solves E1 = 0 rationally."""
    Rl, _, _ = effective()
    Mx = sp.diag(-1, 1, 1) * Rl
    c2_, c1_, _ = inv_c(Mx)
    Qy = sp.Poly(sp.expand(c2_**2 - 3 * c1_), yy)
    cfs = [PolyEval(cf, (X, ga, ph, c, s)) for cf in Qy.all_coeffs()]
    cc, ss = Fr(3, 5), Fr(4, 5)
    for k in range(1, 25):
        for g in range(-16, 17):
            if g == 0:
                continue
            phi_, ga_ = Fr(k, 4), Fr(g, 4)
            q2, q1, q0 = [f((Fr(1), ga_, phi_, cc, ss)) for f in cfs]
            if q2 == 0:
                continue
            Dl = q1 * q1 - 4 * q2 * q0
            if Dl < 0 or math.isqrt(Dl.numerator)**2 != Dl.numerator or math.isqrt(Dl.denominator)**2 != Dl.denominator:
                continue
            r0 = Fr(math.isqrt(Dl.numerator), math.isqrt(Dl.denominator))
            for sg in (1, -1):
                y_ = (-q1 + sg * r0) / (2 * q2)
                if y_ + 2 * ga_ == 0:
                    continue
                A_ = 1 if ga_ > 0 else -1
                u_ = ga_ * A_
                pt = (Fr(1), phi_, Fr(A_), y_ + ga_, u_, cc, ss)
                return ("Type IV with E1 = 0, E2 != 0 (phi = %s, gamma = %s, y = %s)" % (phi_, ga_, y_), pt)
    return None


def type2_failures(h):
    """Each Type II point of D = 0 below fails exactly one condition and is outside the description; a Type IV point fails only E2."""
    ok, rows = True, []
    for lab, pt in triple_rank1_points():
        sm = summ(pt)
        good = (sm["type"] == "II" and sm["ht"] == "II" and sm["cl"] == "II" and sm["D0"] and sm["E1z"] and sm["E2z"] and sm["Tz"]
                and sm["rank"] == 1 and sm["mult"] == 3 and not in_description(sm))
        ok = ok and good
        rows.append("%s: type %s (he_type_exact %s, classify %s), D = 0, triple root, rank %d, E1 = E2 = 0, T = 0, in description: %s" % (
            lab, sm["type"], sm["ht"], sm["cl"], sm["rank"], in_description(sm)))
    for lab, pt, field, root, sd in minus8_points():
        QS.d = field if field else QS.d
        sm = summ(pt)
        inv = point_invariants(pt)
        p_ = inv["lam"]
        q_ = -(inv["M"][0][0] + inv["M"][1][1] + inv["M"][2][2]) - 2 * p_  # sign fixed below
        q_ = (inv["M"][0][0] + inv["M"][1][1] + inv["M"][2][2]) - 2 * p_
        good = (sm["type"] == "II" and sm["D0"] and (not sm["E1z"]) and sm["E2z"] and (not sm["Tz"]) and sm["rank"] == 2 and sm["mult"] == 2
                and q_ == -8 * p_ and not is0(p_) and sm.get("ht", "II") == "II" and sm.get("cl", "II") == "II" and not in_description(sm))
        ok = ok and good
        rows.append("%s: type %s, D = 0, double root p = %s, third root q = %s = -8p, rank %d, E2 = 0, E1 = %s != 0, T != 0, in description: %s" % (
            lab, sm["type"], fmt_num(p_), fmt_num(q_), sm["rank"], fmt_num(inv["E1"]), in_description(sm)))
    for lab, pt in generic_points():
        sm = summ(pt)
        good = (sm["type"] == "II" and sm["ht"] == "II" and sm["cl"] == "II" and sm["D0"] and (not sm["E1z"]) and (not sm["E2z"]) and (not sm["Tz"])
                and sm["rank"] == 2 and not in_description(sm))
        ok = ok and good
        rows.append("%s: type %s (he_type_exact %s), D = 0, rank %d, E1 != 0, E2 != 0, T != 0, in description: %s" % (
            lab, sm["type"], sm["ht"], sm["rank"], in_description(sm)))
    lab, pt = e1_conic_point()
    sm = summ(pt)
    inv = point_invariants(pt)
    good = (sm["type"] == "IV" and sm["ht"] == "IV" and sm["E1z"] and (not sm["E2z"]) and (not sm["Tz"]) and (not sm["D0"])
            and inv["D"] == -inv["E2"]**2 / 27 and not in_description(sm))
    ok = ok and good
    rows.append("%s: type %s (he_type_exact %s), E1 = 0, E2 = %s != 0, T != 0, D = -E2^2/27 = %s < 0, in description: %s" % (
        lab, sm["type"], sm["ht"], fmt_num(inv["E2"], 5), fmt_num(inv["D"], 5), in_description(sm)))
    h.check("type2_failures", ok, "; ".join(rows))


def description(h):
    """Type III <=> E1 = 0, E2 = 0, T != 0, tested at exact points of every type, with the types from two routes."""
    pts = [(lab, pt, None) for lab, pt in triple_rank1_points()] + [(lab, pt, fld) for lab, pt, fld, _, _ in minus8_points()] \
        + [(lab, pt, None) for lab, pt in generic_points() + typeI_points() + random_points() + [e1_conic_point()]]
    ok, rows = True, []
    for lab, pt, fld in pts:
        if fld:
            QS.d = fld
        sm = summ(pt)
        iii = sm["type"] == "III"
        good = in_description(sm) == iii and sm["type"] in ("I", "II", "IV") and (sm.get("ht", sm["type"]) == sm["type"]) \
            and (sm.get("cl", sm["type"]) == sm["type"])
        ok = ok and good
        rows.append("%s -> %s" % (lab, sm["type"] + (" (in description)" if in_description(sm) else "")))
    n3 = 0
    for sd, root in TYPE3_SEEDS:
        pt, lam1, phi_, field = make_point(*sd, root, 2)
        sm = summ(pt)
        good = sm["type"] == "III" and in_description(sm)
        ok = ok and good
        n3 += 1
    h.check("description", ok and len(pts) > 14, "in the description iff Type III at %d exact points of Types I, II and IV (%s) and at the %d Type III points" % (
        len(pts), "; ".join(rows), n3))


# ----------------------------------------------------------------------------------------------------------------------
# the tanh wall
# ----------------------------------------------------------------------------------------------------------------------
WALLS = {  # name -> (sigma, R_b, speed v of the paper, mask)
    "Le": (mp.mpf(8), mp.mpf(1), mp.mpf(1) / 2),
    "Rodal": (mp.mpf(4), mp.mpf(5), mp.mpf(1)),
}


def jets_mp(r, sig, Rb):
    """(F, f', f'', f''') of f = [tanh sigma (r + Rb) - tanh sigma (r - Rb)] / (2 tanh sigma Rb), own mpmath formulas."""
    tp, tm = mp.tanh(sig * (r + Rb)), mp.tanh(sig * (r - Rb))
    N = 2 * mp.tanh(sig * Rb)
    f0 = (tp - tm) / N
    f1 = sig * (tm**2 - tp**2) / N
    f2 = 2 * sig**2 * (tm * (1 - tm**2) - tp * (1 - tp**2)) / N
    f3 = -2 * sig**3 * ((1 - tp**2) * (1 - 3 * tp**2) - (1 - tm**2) * (1 - 3 * tm**2)) / N
    if MUT == "jets-sign":
        f2 = -f2
    return 1 - f0, f1, f2, f3


_SH = {}


def shadow():
    """E1, E2 as polynomials in w = c^2 with coefficients in (phi, A', B') = (vF, v r a, v r^2 b), and the factor P of their
    resultant in w: P = 0 is necessary for a common root, P(phi, A', B') = v^10 (P0 + P1 v^2 + ... + P4 v^8) on the wall."""
    ckey = "s" + ("-P" if MUT == "P-term" else "")
    if ckey in _SH:
        return _SH[ckey]
    al = algebra()
    out = []
    for key in ("E1", "E2"):
        e, A_, B_ = eff(al[key])
        e = sp.expand(e.subs({u: 1, A_: Ap, B_: Bp}).subs(c, sp.sqrt(w)))
        out.append(e)
    e1w, e2w = out
    Res = sp.expand(sp.resultant(e1w, e2w, w))
    fl = sp.factor_list(Res)
    big = max(fl[1], key=lambda fk: len(sp.Poly(fk[0], ph, Ap, Bp).terms()))
    P = sp.expand(big[0])
    rest = sp.expand(Res / P)
    PP = sp.Poly(P, ph, Ap, Bp)
    if MUT == "P-term":
        terms = PP.terms()[1:]
    else:
        terms = PP.terms()
    ptab = [(mn, mp.mpf(int(cf.p)) / int(cf.q)) for mn, cf in terms]
    ecoef = []
    for e in (e1w, e2w):
        Pw = sp.Poly(e, w)
        ecoef.append([[(mn, mp.mpf(int(cf.p)) / int(cf.q)) for mn, cf in sp.Poly(co, ph, Ap, Bp).terms()] for co in Pw.all_coeffs()[::-1]])
    eexact = [[[(mn, int(cf.p), int(cf.q)) for mn, cf in sp.Poly(co, ph, Ap, Bp).terms()] for co in sp.Poly(e, w).all_coeffs()[::-1]]
              for e in (e1w, e2w)]
    _SH[ckey] = dict(e1w=e1w, e2w=e2w, Res=Res, fl=fl, P=P, rest=rest, ptab=ptab, ecoef=ecoef, eexact=eexact, nterms=len(PP.terms()),
                    degs=sorted({sum(mn) for mn, _ in PP.terms()}))
    return _SH[ckey]


def eval_terms(tab, vals):
    tot = mp.mpf(0)
    for mn, cf in tab:
        t = cf
        for xx, e in zip(vals, mn):
            if e:
                t *= xx**e
        tot += t
    return tot


def Pquartic(F_, a_, b_, r_):
    """Coefficients [Q0, ..., Q4] with P(vF, v r a, v r^2 b) = v^10 sum Q_k v^(2k), at the wall point r."""
    sh = shadow()
    base = (F_, a_ * r_, b_ * r_**2)
    co = [mp.mpf(0)] * 5
    for mn, cf in sh["ptab"]:
        d = sum(mn)
        t = cf
        for xx, e in zip(base, mn):
            if e:
                t *= xx**e
        co[(d - 10) // 2] += t
    return co


def Eeval(phi_, Ap_, Bp_, w_):
    """E1 and E2 of the effective point (u = 1), each divided by the sum of the moduli of its terms (a dimensionless residual)."""
    sh = shadow()
    vals = (phi_, Ap_, Bp_)
    out = []
    for co in sh["ecoef"]:
        cs_ = [eval_terms(t, vals) for t in co]
        scale = sum(abs(cs_[k]) * abs(w_)**k for k in range(len(cs_)))
        out.append(sum(cs_[k] * w_**k for k in range(len(cs_))) / (scale if scale != 0 else 1))
    return out


def wcommon(phi_, Ap_, Bp_, tol=mp.mpf(10)**-10):
    """Candidate common roots w of E1 and E2 at the effective point: roots of the quadratic E1 that make E2 (relatively) small."""
    sh = shadow()
    vals = (phi_, Ap_, Bp_)
    c1_ = [eval_terms(t, vals) for t in sh["ecoef"][0]]
    out = []
    if c1_[2] == 0:
        return out
    disc = c1_[1]**2 - 4 * c1_[2] * c1_[0]
    if disc < 0:
        return out
    for sg in (1, -1):
        w0 = (-c1_[1] + sg * mp.sqrt(disc)) / (2 * c1_[2])
        if abs(Eeval(phi_, Ap_, Bp_, w0)[1]) < tol:
            out.append(w0)
    return out


def polish_vw(F_, a_, b_, r_, v0, w0):
    """Newton on (v, w) for E1 = E2 = 0 at the wall point r, from the candidate (v0, w0). Returns (v, w, residual) or None."""
    f = lambda vv, ww: Eeval(vv * F_, vv * a_ * r_, vv * b_ * r_**2, ww)  # noqa: E731
    try:
        sol = mp.findroot(f, (v0, w0), tol=mp.mpf(10)**-(mp.mp.dps - 8), maxsteps=15)
    except (ValueError, ZeroDivisionError):
        return None
    vv, ww = sol[0], sol[1]
    res = max(abs(z) for z in f(vv, ww))
    return vv, ww, res


def curve_at(r_, sig, Rb, dps=60):
    """All (v, c^2, T) with E1 = E2 = 0 at the wall point r (v > 0, 0 < c^2 < 1, T != 0), from the quartic in v^2 and the common
    root w, each polished by Newton to a residual below 10^-(dps - 15)."""
    mp.mp.dps = dps
    r_ = mp.mpf(r_)
    F_, a_, b_, _ = jets_mp(r_, sig, Rb)
    return curve_from_jets(F_, a_, b_, r_, dps)


def admit(F_, a_, b_, r_, vv, w0, dps=60, filter_T=True):
    """Newton polish of a candidate (v, w) and the admission tests: residual below 10^-(dps - 15), 0 < w < 1, v > 0 and T != 0."""
    pl = polish_vw(F_, a_, b_, r_, vv, w0)
    if pl is None or pl[2] > mp.mpf(10)**-(dps - 15) or not (0 < pl[1] < 1) or pl[0] <= 0:
        return None
    T_ = b_ + a_ / r_
    Tabs = abs(T_) / (abs(b_) + abs(a_ / r_))
    if filter_T and MUT != "T-filter-off" and Tabs < mp.mpf(10)**-(dps // 2):
        return None
    return (pl[0], pl[1], T_)


def curve_from_jets(F_, a_, b_, r_, dps=60, filter_T=True):
    mp.mp.dps = dps
    co = Pquartic(F_, a_, b_, r_)
    while co and co[-1] == 0:
        co = co[:-1]
    if len(co) < 2:
        return []
    try:
        rts = mp.polyroots(co[::-1], maxsteps=400, extraprec=300)
    except mp.libmp.libhyper.NoConvergence:
        return []
    out = []
    for z in rts:
        if abs(mp.im(z)) > mp.mpf(10)**-(dps // 2) * max(1, abs(z)) or mp.re(z) <= 0:
            continue
        vv = mp.sqrt(mp.re(z))
        for w0 in wcommon(vv * F_, vv * a_ * r_, vv * b_ * r_**2):
            pnt = admit(F_, a_, b_, r_, vv, w0, dps, filter_T)
            if pnt is not None:
                out.append(pnt)
    out.sort(key=lambda t: t[0])
    return out


# ----------------------------------------------------------------------------------------------------------------------
# the tanh wall: checks
# ----------------------------------------------------------------------------------------------------------------------
_NB = {}


def block_num(vals):
    """The frame block R3 at a point of mp numbers (v, F, a, b, u, c, s), at the current precision."""
    if "b" not in _NB:
        R3 = algebra()["R3"]
        _NB["b"] = [[[(mn, int(cf.p), int(cf.q)) for mn, cf in sp.Poly(sp.expand(R3[i, j]), *VARS7).terms()] for j in range(3)] for i in range(3)]
    out = []
    for i in range(3):
        row = []
        for j in range(3):
            tot = mp.mpf(0)
            for mn, p_, q_ in _NB["b"][i][j]:
                t = mp.mpf(p_) / q_
                for xx, e in zip(vals, mn):
                    if e:
                        t *= xx**e
                tot += t
            row.append(tot)
        out.append(row)
    return out


def point_numeric(sig, Rb, r_, v0, w_):
    """Singular values of N = M - lam0, |N^3|, E1, E2 (frame values) and T at the wall point (r, c^2 = w, v = v0), 60 digits."""
    F_, a_, b_, _ = jets_mp(r_, sig, Rb)
    c_ = mp.sqrt(w_)
    s_ = mp.sqrt(1 - w_)
    pt = (v0, F_, a_, b_, 1 / r_, c_, s_)
    R_ = block_num(pt)
    M_ = mixed_at(R_)
    sv, n3 = numeric_type3(M_)
    c2_, c1_, c0_ = cubic_at(M_)
    e1_, e2_ = E12(c2_, c1_, c0_)
    return dict(sv=sv, n3=n3, E1=e1_, E2=e2_, T=b_ + a_ / r_, F=F_, scale=max(abs(z) for row in R_ for z in row), c=c_, s=s_, R=R_, M=M_)


def tanh_jets_check(h):
    """Own jets against he_type.tanh_jets (floats) and against mp.diff of the profile at 40 digits."""
    ok, rows = True, []
    mp.mp.dps = 40
    for nm, (sig, Rb, v0) in WALLS.items():
        w1 = w2 = mpw = mp.mpf(0)
        fn = lambda xx: (mp.tanh(sig * (xx + Rb)) - mp.tanh(sig * (xx - Rb))) / (2 * mp.tanh(sig * Rb))  # noqa: E731
        for k in range(12):
            r_ = Rb * (mp.mpf(1) / 2 + mp.mpf(k) / 8)
            F_, a_, b_, c3 = jets_mp(r_, sig, Rb)
            f_, f1_, f2_ = HT.tanh_jets(float(r_), float(sig), float(Rb))
            w1 = max(w1, abs(1 - f_ - F_), abs(f1_ - a_) / max(1, abs(a_)), abs(f2_ - b_) / max(1, abs(b_)))
            d = [mp.diff(fn, r_, n) for n in (1, 2, 3)]
            mpw = max(mpw, abs(d[0] - a_) / max(1, abs(a_)), abs(d[1] - b_) / max(1, abs(b_)), abs(d[2] - c3) / max(1, abs(c3)))
        ok = ok and w1 < mp.mpf(10)**-9 and mpw < mp.mpf(10)**-25
        rows.append("%s (sigma = %s, R_b = %s): max |own - he_type.tanh_jets| = %s over 12 radii, max |own - mp.diff| = %s for f', f'', f'''" % (
            nm, sig, Rb, mp.nstr(w1, 3), mp.nstr(mpw, 3)))
    h.check("tanh_jets", ok, "; ".join(rows))


def tanh_planted(h):
    """Planted truth for the numerical path: the solver, given the jets (F, a, b) and r of an exact Type III point, finds its
    (v, c^2); given those of a triple-root rank-1 point (T = 0) it sees the common root and the T filter rejects it."""
    ok, rows = True, []
    for idx in (0, 1):
        sd, root = TYPE3_SEEDS[idx]
        pt, lam1, phi_, field = make_point(*sd, root, 2)
        mp.mp.dps = 60
        F_, a_, b_ = to_mp(pt[1]), to_mp(pt[2]), to_mp(pt[3])
        r_ = 1 / to_mp(pt[4])
        res = curve_from_jets(F_, a_, b_, r_)
        want_w = mp.mpf(sd[0].numerator)**2 / sd[0].denominator**2
        hit = [z for z in res if abs(z[0] - 2) < mp.mpf(10)**-40 and abs(z[1] - want_w) < mp.mpf(10)**-40]
        ok = ok and len(hit) == 1
        rows.append("planted Type III point in Q(sqrt %d), (v, c^2) = (2, %s): the solver returns %d point(s) at r = %s, %s" % (
            root[1], sd[0]**2, len(res), mp.nstr(r_, 8), "recovers it to 1e-40" if hit else "MISSES it"))
    pt = triple_rank1_points()[0][1]
    mp.mp.dps = 60
    F_, a_, b_ = to_mp(pt[1]), to_mp(pt[2]), to_mp(pt[3])
    r_ = 1 / to_mp(pt[4])
    co = Pquartic(F_, a_, b_, r_)
    Pres = abs(sum(co[k] * (mp.mpf(9) / 4)**k for k in range(5))) / sum(abs(co[k]) * (mp.mpf(9) / 4)**k for k in range(5))
    seen = admit(F_, a_, b_, r_, mp.mpf(3) / 2, mp.mpf(9) / 25, filter_T=False)
    gone = admit(F_, a_, b_, r_, mp.mpf(3) / 2, mp.mpf(9) / 25)
    ok = ok and Pres < mp.mpf(10)**-50 and seen is not None and gone is None
    rows.append("planted triple root of rank 1, (v, c^2) = (3/2, 9/25), T = 0: V = 9/4 is a root of the shadow quartic (relative residual %s, a "
                "triple root, so not returned by the root finder), Newton from (3/2, 9/25) converges without the T filter (%s) and the T filter "
                "rejects it (%s)" % (mp.nstr(Pres, 3), seen is not None, gone is None))
    h.check("tanh_planted", ok, "; ".join(rows))


def r_grid(name):
    sig, Rb, v0 = WALLS[name]
    segs = {"Le": [(0.05, 0.3, 0.005), (0.3, 3.0, 0.0005), (3.0, 10.0, 0.01)],
            "Rodal": [(0.5, 3.0, 0.02), (3.0, 9.0, 0.001), (9.0, 25.0, 0.02)]}[name]
    out = []
    for lo, hi, st in segs:
        n = int(round((hi - lo) / st))
        out += [mp.mpf(lo) + (mp.mpf(hi) - mp.mpf(lo)) * k / n for k in range(n)]
    out.append(mp.mpf(segs[-1][1]))
    return out


def gsum(r_, sig, Rb, v0):
    F_, a_, b_, _ = jets_mp(r_, sig, Rb)
    co = Pquartic(F_, a_, b_, r_)
    V = v0**2
    return sum(co[k] * V**k for k in range(5))


def points_at_speed(name):
    """All Type III points of the wall at the paper's speed: sign changes of P(v F, v r a, v r^2 b) over a fine grid of r,
    bisection at 60 digits, the common root w, Newton on (r, w), and the checks. Also near-tangent minima of |P|."""
    sig, Rb, v0 = WALLS[name]
    grid = r_grid(name)
    mp.mp.dps = 30
    vals = [gsum(r_, sig, Rb, v0) for r_ in grid]
    brackets = [(grid[i], grid[i + 1]) for i in range(len(grid) - 1) if vals[i] * vals[i + 1] < 0]
    tang = [grid[i] for i in range(1, len(grid) - 1) if abs(vals[i]) < abs(vals[i - 1]) and abs(vals[i]) < abs(vals[i + 1])
            and vals[i - 1] * vals[i] > 0 and vals[i] * vals[i + 1] > 0 and abs(vals[i]) < mp.mpf(10)**-3 * min(abs(vals[i - 1]), abs(vals[i + 1]))]
    mp.mp.dps = 60
    out, rejected = [], []
    for lo, hi in brackets:
        lo, hi = mp.mpf(lo), mp.mpf(hi)
        glo = gsum(lo, sig, Rb, v0)
        if glo * gsum(hi, sig, Rb, v0) >= 0:
            rejected.append((lo, "no sign change at 60 digits"))
            continue
        for _ in range(170):
            mid = (lo + hi) / 2
            gm = gsum(mid, sig, Rb, v0)
            if gm * glo > 0:
                lo, glo = mid, gm
            else:
                hi = mid
        r0 = (lo + hi) / 2
        F_, a_, b_, _ = jets_mp(r0, sig, Rb)
        cand = wcommon(v0 * F_, v0 * a_ * r0, v0 * b_ * r0**2)
        if not cand:
            rejected.append((r0, "no common root w"))
            continue
        for w0 in cand:
            def fr(rr, ww):
                Fx, ax, bx, _ = jets_mp(rr, sig, Rb)
                return Eeval(v0 * Fx, v0 * ax * rr, v0 * bx * rr**2, ww)
            try:
                sol = mp.findroot(fr, (r0, w0), tol=mp.mpf(10)**-52, maxsteps=40)
            except (ValueError, ZeroDivisionError):
                rejected.append((r0, "Newton failed"))
                continue
            rr, ww = sol[0], sol[1]
            res = max(abs(z) for z in fr(rr, ww))
            Fx, ax, bx, _ = jets_mp(rr, sig, Rb)
            Tv = bx + ax / rr
            if res < mp.mpf(10)**-45 and 0 < ww < 1 and abs(Tv) > mp.mpf(10)**-20:
                out.append(dict(r=rr, w=ww, res=res, T=Tv, F=Fx))
            else:
                rejected.append((r0, "residual %s, w %s" % (mp.nstr(res, 3), mp.nstr(ww, 5))))
    return out, rejected, tang, len(brackets)


def in_mask(name, r_):
    sig, Rb, v0 = WALLS[name]
    if name == "Le":
        f_ = 1 - jets_mp(r_, sig, Rb)[0]
        return 0.1 <= f_ <= 0.9
    return Rb - 3 / sig <= r_ <= Rb + 3 / sig


def tanh_points(h):
    ok, rows = True, []
    for name in ("Le", "Rodal"):
        sig, Rb, v0 = WALLS[name]
        pts, rejected, tang, nbr = points_at_speed(name)
        shift_ = mp.mpf(10)**-6 if MUT == "tanh-off-root" else 0
        good = len(pts) >= 1
        desc = []
        for pnt in pts:
            r_ = pnt["r"]
            nb = point_numeric(sig, Rb, r_ + shift_, v0, pnt["w"])
            sv = nb["sv"]
            ok_pt = (sv[2] < mp.mpf(10)**-45 * max(1, nb["scale"]) and sv[1] > mp.mpf(10)**-8 * max(1, nb["scale"]) and
                     nb["n3"] < mp.mpf(10)**-40 * max(1, nb["scale"])**3 and abs(nb["E1"]) < mp.mpf(10)**-40 * max(1, nb["scale"])**2)
            good = good and ok_pt
            th = mp.degrees(mp.acos(nb["c"]))
            desc.append("r = %s (f = %s, in mask: %s), cos(theta) = %s (theta = %s deg), T = %s, singular values of M - lam: %s, %s, %s, |N^3| = %s" % (
                mp.nstr(r_, 12), mp.nstr(1 - nb["F"], 6), in_mask(name, r_), mp.nstr(nb["c"], 12), mp.nstr(th, 8), mp.nstr(nb["T"], 6),
                mp.nstr(sv[0], 6), mp.nstr(sv[1], 6), mp.nstr(sv[2], 3), mp.nstr(nb["n3"], 3)))
        ok = ok and good
        rows.append("%s (sigma = %s, R_b = %s, v = %s): %d sign changes of P over %d grid radii, %d Type III points of residual < 1e-45, %d rejected %s, "
                    "%d near-tangent minima of |P|; %s" % (
                        name, sig, Rb, v0, nbr, len(r_grid(name)), len(pts), len(rejected), [z[1] for z in rejected], len(tang), "; ".join(desc)))
    h.check("tanh_points", ok, " ... ".join(rows))


_PTS = {}


def points_cached(name):
    key = (name, MUT if MUT in ("jets-sign", "P-term") else "")
    if key not in _PTS:
        _PTS[key] = points_at_speed(name)
    return _PTS[key]


def eq_end(sig, Rb, r0):
    """The equator end of branch A. At c = 0 the resultant of E1 and E2 in B' is 36 A'^6 (A' + 2 phi)^6, so a common root has
    A' = -2 phi, i.e. r a = -2 F, and then E1 = 0 gives B' = 2 phi (1 +- 2 phi), i.e. v = |r^2 b - 2F| / (4 F^2)."""
    def g(rr):
        Fx, ax, bx, _ = jets_mp(rr, sig, Rb)
        return rr * ax + (3 if MUT == "eq-end-typo" else 2) * Fx
    rr = mp.findroot(g, r0, tol=mp.mpf(10)**-52)
    Fx, ax, bx, _ = jets_mp(rr, sig, Rb)
    vv = abs(rr**2 * bx - 2 * Fx) / (4 * Fx**2)
    return rr, vv


def ax_start(sig, Rb, rlo, rhi):
    """The start of branch B: the radius where the leading coefficient Q4 of the shadow quartic changes sign, so that a root V = v^2
    goes to infinity."""
    q4 = lambda rr: Pquartic(*jets_mp(rr, sig, Rb)[:3], rr)[4]  # noqa: E731
    lo, hi = mp.mpf(rlo), mp.mpf(rhi)
    glo = q4(lo)
    assert glo * q4(hi) < 0
    for _ in range(160):
        mid = (lo + hi) / 2
        gm = q4(mid)
        if gm * glo > 0:
            lo, glo = mid, gm
        else:
            hi = mid
    return (lo + hi) / 2


# the three Newton points of #72 on the tanh wall, he_surface.md section 4: (sigma, R_b = 1) -> (r, v, cos theta)
POINTS72 = [(1, mp.mpf("2.0656"), mp.mpf("1.1282"), mp.mpf("0.7768")), (4, mp.mpf("0.898"), mp.mpf("2.2089"), mp.mpf("0.2116")),
            (8, mp.mpf("0.916"), mp.mpf("2.0194"), mp.mpf("0.6437"))]


def solve_near(sig, Rb, r0, v0, c0):
    """Newton on (r, w) for E1 = E2 = 0 at fixed v from (r0, c0^2)."""
    def fr(rr, ww):
        Fx, ax, bx, _ = jets_mp(rr, sig, Rb)
        return Eeval(v0 * Fx, v0 * ax * rr, v0 * bx * rr**2, ww)
    sol = mp.findroot(fr, (r0, c0**2), tol=mp.mpf(10)**-52, maxsteps=60)
    return sol[0], sol[1], max(abs(z) for z in fr(sol[0], sol[1]))


def tanh_curve(h):
    """The curve of Type III points over r: branch A from the centre (v -> infinity as r decreases) to the equator, branch B from
    infinity on the axis side to r -> infinity (v -> 0). Counts per radius, the two ends, the minimum of v on A, a table, the
    crossings of v, and #72's three Newton points."""
    ok, rows = True, []
    mp.mp.dps = 60
    for name, (lo, hi, n) in (("Le", (0.7, 2.0, 650)), ("Rodal", (4.0, 7.0, 600))):
        sig, Rb, v0 = WALLS[name]
        mp.mp.dps = 60
        A, B, counts = [], [], []
        for k in range(n + 1):
            r_ = mp.mpf(lo) + (mp.mpf(hi) - mp.mpf(lo)) * k / n
            pts = sorted(curve_at(r_, sig, Rb), key=lambda z: z[1])
            counts.append(len(pts))
            if len(pts) == 2:
                A.append((r_, pts[0][0], pts[0][1]))
                B.append((r_, pts[1][0], pts[1][1]))
            elif len(pts) == 1:
                (A if not B else B).append((r_, pts[0][0], pts[0][1]))
        pattern = [counts[0]] + [counts[i] for i in range(1, len(counts)) if counts[i] != counts[i - 1]]
        i_ax = next(i for i in range(len(counts)) if counts[i] == 2)
        i_eq = next(i for i in range(i_ax, len(counts)) if counts[i] == 1)
        dr = (mp.mpf(hi) - mp.mpf(lo)) / n
        r_i = lambda i: mp.mpf(lo) + dr * i  # noqa: E731
        okp = pattern == [1, 2, 1]
        rax = ax_start(sig, Rb, r_i(i_ax - 3), r_i(i_ax))
        pB = sorted(curve_at(rax + mp.mpf(10)**-6, sig, Rb), key=lambda z: z[1])
        axis_ok = len(pB) == 2 and pB[1][0] > 100 and pB[1][1] > mp.mpf("0.999") and abs(rax - r_i(i_ax)) < 4 * dr
        req, veq = eq_end(sig, Rb, r_i(i_eq))
        neq = point_numeric(sig, Rb, req, veq, mp.mpf(0))
        Fq, aq, bq, _ = jets_mp(req, sig, Rb)
        Qeq = 4 * veq**2 * aq**2 * Fq * (Fq + req * aq) + (aq + req * bq)**2  # he_type.md: Q_eq = 0
        edge = req**2 * veq**2 * aq**4 - (aq + req * bq)**2  # he_type.md: the edge v r f'^2 = |(r f')'|
        eq_ok = (abs(req - r_i(i_eq)) < 2 * dr and abs(neq["E1"]) < mp.mpf(10)**-40 and abs(neq["E2"]) < mp.mpf(10)**-40
                 and neq["sv"][1] < mp.mpf(10)**-40 * max(1, neq["scale"]) and neq["sv"][0] > mp.mpf(10)**-6
                 and abs(Qeq) < mp.mpf(10)**-40 * max(1, neq["scale"]**2) and abs(edge) < mp.mpf(10)**-40 * max(1, neq["scale"]**2))
        amin = min(A, key=lambda z: z[1])
        rT = mp.findroot(lambda rr: (lambda j: j[2] + j[1] / rr)(jets_mp(rr, sig, Rb)), amin[0])
        pT = [curve_at(rT + sg_ * mp.mpf(10)**-5, sig, Rb) for sg_ in (-1, 1)]
        tT_ok = all(len(z) >= 1 and min(q[0] for q in z) < mp.mpf("0.01") and min(mp.sqrt(q[1]) for q in z) < mp.mpf("0.01") for z in pT)
        cross = 0
        for br in (A, B):
            for i in range(len(br) - 1):
                if (br[i][1] - v0) * (br[i + 1][1] - v0) < 0 and br[i + 1][0] - br[i][0] < 3 * dr:
                    cross += 1
        npts = len(points_cached(name)[0])
        ok = ok and okp and axis_ok and eq_ok and tT_ok and cross == npts
        sel = {"Le": [0.8, 0.9, 0.916, 0.95, 1.0, 1.01, 1.03, 1.05, 1.08, 1.1, 1.13, 1.2, 1.5, 2.0, 3.0],
               "Rodal": [4.5, 4.7, 4.9, 5.0, 5.05, 5.1, 5.2, 5.3, 5.35, 5.5, 6.0, 8.0, 12.0]}[name]
        tab = []
        for r_ in sel:
            pts = curve_at(mp.mpf(r_), sig, Rb)
            tab.append("r = %s (f = %.4f): %s" % (r_, float(1 - jets_mp(mp.mpf(r_), sig, Rb)[0]),
                       ", ".join("(v = %s, cos = %s)" % (mp.nstr(z[0], 5), mp.nstr(mp.sqrt(z[1]), 5)) for z in pts) or "none"))
        rows.append("%s (sigma = %s, R_b = %s): the number of Type III points per radius over r in [%s, %s] at %d radii runs %s. Branch B starts at "
                    "r = %s where the leading coefficient of the shadow quartic changes sign (v -> infinity, c -> 1): at r + 1e-6 it has v = %s, c^2 = %s. "
                    "Branch A ends on the equator at r = %s, v = %s (r a = -2F, v = |r^2 b - 2F| / (4 F^2); there E1 = %s, E2 = %s, singular values "
                    "%s, %s, %s: a triple root of rank 1, Type II; he_type.md's Q_eq = %s and r^2 v^2 a^4 - (a + r b)^2 = %s there, so it is the crossing of "
                    "the two curves of the equator). v on A has its minimum %s on the grid at r = %s (cos = %s), and A passes through v = 0, c = 0 "
                    "at r_T = %s, the root of (r f')' = 0 (T = 0): at r_T -+ 1e-5 the smallest (v, cos) are (%s, %s) and (%s, %s). v crosses the paper's "
                    "v = %s %d times over the branches and tanh_points found %d points. Table: %s" % (
                        name, sig, Rb, lo, hi, n + 1, pattern, mp.nstr(rax, 12), mp.nstr(pB[1][0], 6) if len(pB) == 2 else "none",
                        mp.nstr(pB[1][1], 8) if len(pB) == 2 else "none", mp.nstr(req, 12), mp.nstr(veq, 10), mp.nstr(abs(neq["E1"]), 2),
                        mp.nstr(abs(neq["E2"]), 2), mp.nstr(neq["sv"][0], 5), mp.nstr(neq["sv"][1], 2), mp.nstr(neq["sv"][2], 2),
                        mp.nstr(abs(Qeq), 2), mp.nstr(abs(edge), 2), mp.nstr(amin[1], 6), mp.nstr(amin[0], 6), mp.nstr(mp.sqrt(amin[2]), 6), mp.nstr(rT, 10),
                        mp.nstr(min(q[0] for q in pT[0]), 3), mp.nstr(min(mp.sqrt(q[1]) for q in pT[0]), 3),
                        mp.nstr(min(q[0] for q in pT[1]), 3), mp.nstr(min(mp.sqrt(q[1]) for q in pT[1]), 3), v0, cross, npts, "; ".join(tab)))
    # #72's Newton points at R_b = 1: sigma R_b = 1, 4, 8
    rr72 = []
    ok72 = True
    for sig_, r0, v0_, c0 in POINTS72:
        sol = solve_near(mp.mpf(sig_), mp.mpf(1), r0, v0_, c0)
        near = abs(sol[0] - r0) < mp.mpf("0.005") * max(1, r0) and abs(mp.sqrt(sol[1]) - c0) < mp.mpf("0.002") and sol[2] < mp.mpf(10)**-45
        ok72 = ok72 and near
        rr72.append("sigma R_b = %s: Newton at v = %s from (r, cos) = (%s, %s) gives r = %s, cos = %s (residual %s)" % (
            sig_, mp.nstr(v0_, 5), mp.nstr(r0, 5), mp.nstr(c0, 5), mp.nstr(sol[0], 8), mp.nstr(mp.sqrt(sol[1]), 6), mp.nstr(sol[2], 2)))
    ok = ok and ok72
    h.check("tanh_curve", ok, " ... ".join(rows) + " ... #72's points (he_surface.md section 4) lie on the curve: " + "; ".join(rr72))


# ----------------------------------------------------------------------------------------------------------------------
# interval Newton enclosures
# ----------------------------------------------------------------------------------------------------------------------
def iv_tanh(xx):
    e = mp.iv.exp(2 * xx)
    return (e - 1) / (e + 1)


def iv_jets(r_, sig, Rb):
    sig, Rb = mp.iv.mpf(sig), mp.iv.mpf(Rb)
    tp, tm, tR = iv_tanh(sig * (r_ + Rb)), iv_tanh(sig * (r_ - Rb)), iv_tanh(sig * Rb)
    N = 2 * tR
    f0 = (tp - tm) / N
    f1 = sig * (tm * tm - tp * tp) / N
    f2 = 2 * sig**2 * (tm * (1 - tm * tm) - tp * (1 - tp * tp)) / N
    f3 = -2 * sig**3 * ((1 - tp * tp) * (1 - 3 * tp * tp) - (1 - tm * tm) * (1 - 3 * tm * tm)) / N
    return 1 - f0, f1, f2, f3


def iv_E(r_, w_, v0, sig, Rb):
    """E1, E2 (u = 1 effective form, unscaled) and their partial derivatives in r and w, as intervals, at (r, w)."""
    iv = mp.iv
    v0 = iv.mpf(v0)
    F_, a_, b_, c3 = iv_jets(r_, sig, Rb)
    vals = (v0 * F_, v0 * r_ * a_, v0 * r_ * r_ * b_)
    dv = (-v0 * a_, v0 * (a_ + r_ * b_), v0 * (2 * r_ * b_ + r_ * r_ * c3))
    sh = shadow()
    outE, outR, outW = [], [], []
    for k in range(2):
        E_, dEr, dEw = iv.mpf(0), iv.mpf(0), iv.mpf(0)
        for m, tab in enumerate(sh["eexact"][k]):
            C_, dC = iv.mpf(0), iv.mpf(0)
            for mn, p_, q_ in tab:
                cf = iv.mpf(p_) / q_
                C_ += cf * vals[0]**mn[0] * vals[1]**mn[1] * vals[2]**mn[2]
                for t_ in range(3):
                    if mn[t_]:
                        ex = list(mn)
                        ex[t_] -= 1
                        dC += cf * mn[t_] * vals[0]**ex[0] * vals[1]**ex[1] * vals[2]**ex[2] * dv[t_]
            E_ += C_ * w_**m
            dEr += dC * w_**m
            if m:
                dEw += m * C_ * w_**(m - 1)
        outE.append(E_)
        outR.append(dEr)
        outW.append(dEw)
    return outE, outR, outW


def krawczyk(name, rc, wc, dr, dw, v0):
    """Krawczyk test for G = (E1, E2)(r, w) at fixed v0 on the box rc +- dr, wc +- dw. Returns (contracts, T enclosure, details)."""
    iv = mp.iv
    sig, Rb, _ = WALLS[name]
    iv.dps = mp.mp.dps
    X = [iv.mpf([rc - dr, rc + dr]), iv.mpf([wc - dw, wc + dw])]
    xm = [iv.mpf(rc), iv.mpf(wc)]
    Gm, dRm, dWm = iv_E(xm[0], xm[1], v0, sig, Rb)
    Jm = mp.matrix([[dRm[0].mid, dWm[0].mid], [dRm[1].mid, dWm[1].mid]])
    Y = Jm**-1
    GX, dRX, dWX = iv_E(X[0], X[1], v0, sig, Rb)
    JX = [[dRX[0], dWX[0]], [dRX[1], dWX[1]]]
    Yi = [[iv.mpf(Y[i, j]) for j in range(2)] for i in range(2)]
    corr = [[(1 if i == j else 0) - (Yi[i][0] * JX[0][j] + Yi[i][1] * JX[1][j]) for j in range(2)] for i in range(2)]
    dX = [X[0] - xm[0], X[1] - xm[1]]
    K = []
    for i in range(2):
        K.append(xm[i] - (Yi[i][0] * Gm[0] + Yi[i][1] * Gm[1]) + (corr[i][0] * dX[0] + corr[i][1] * dX[1]))
    inside = all(K[i].a > X[i].a and K[i].b < X[i].b for i in range(2))
    F_, a_, b_, _ = iv_jets(X[0], sig, Rb)
    T = b_ + a_ / X[0]
    Tpos = T.a > 0 or T.b < 0
    cond = (corr[0][0].b - corr[0][0].a, max(abs(corr[i][j].a) + abs(corr[i][j].b) for i in range(2) for j in range(2)))
    return inside and Tpos and F_.a > 0 and a_.b < 0 and X[1].a > 0 and X[1].b < 1, T, (K, X, cond)


def tanh_enclosure(h):
    """Krawczyk enclosures of the Type III points of tanh_points: a box with a unique zero of (E1, E2) at the paper's v, with T != 0
    on it. A box shifted by 20 half-widths must fail (planted negative)."""
    ok, rows = True, []
    mp.mp.dps = 60
    for name in ("Le", "Rodal"):
        sig, Rb, v0 = WALLS[name]
        pts = points_cached(name)[0]
        dr = mp.mpf(10)**-25
        dw = mp.mpf(10)**-25
        sh_ = 20 * dr if MUT == "enclosure-shift" else 0
        for pnt in pts:
            good, T, det = krawczyk(name, pnt["r"] + sh_, pnt["w"], dr, dw, v0)
            bad, _, _ = krawczyk(name, pnt["r"] + 20 * dr, pnt["w"], dr, dw, v0)
            ok = ok and good and not bad
            rows.append("%s r = %s, c^2 = %s: Krawczyk operator inside the box %s (half-widths 1e-25), T = %s with a relative enclosure width %s, "
                        "shifted box rejected %s" % (name, mp.nstr(pnt["r"], 14), mp.nstr(pnt["w"], 14), good, mp.nstr(T.mid, 8),
                                                     mp.nstr((T.b - T.a).mid / abs(T.mid), 2), not bad))
    h.check("tanh_enclosure", ok, "; ".join(rows))


# ----------------------------------------------------------------------------------------------------------------------
# an independent route to the frame Ricci tensor: second differences of the lab metric
# ----------------------------------------------------------------------------------------------------------------------
def lab_frame_ricci(sig, Rb, v0, r_, c_, s_, dps=70):
    """Frame Ricci tensor of the lab metric ds^2 = -dt^2 + (dx - v f(r_s) dt)^2 + dy^2 + dz^2, r_s^2 = (x - v t)^2 + y^2 + z^2,
    at t = 0, (x, y, z) = r (c, s, 0), in the frame E_0 = (1, v f, 0, 0), E_1 = d_x, E_2 = d_y, E_3 = d_z. Own Christoffel and
    Ricci code, with the first and second derivatives of g from mp.diff at `dps` digits; it shares no code with the comoving jets of
    review_he_surface_lib. Returns the 4x4 frame matrix (mp numbers)."""
    mp.mp.dps = dps
    sg = -1 if MUT == "lab-sign" else 1

    def prof(rs):
        return (mp.tanh(sig * (rs + Rb)) - mp.tanh(sig * (rs - Rb))) / (2 * mp.tanh(sig * Rb))

    def fl(t_, x_, y_, z_):
        return prof(mp.sqrt((x_ - v0 * t_)**2 + y_**2 + z_**2))

    comps = {(0, 0): lambda t_, x_, y_, z_: -1 + (v0 * fl(t_, x_, y_, z_))**2,
             (0, 1): lambda t_, x_, y_, z_: -sg * v0 * fl(t_, x_, y_, z_)}
    X0 = (mp.mpf(0), r_ * c_, r_ * s_, mp.mpf(0))
    g = mp.zeros(4, 4)
    for i in range(4):
        g[i, i] = 1
    g[0, 0] = -1
    for (i, j), fn in comps.items():
        g[i, j] = g[j, i] = fn(*X0)
    dg = [[[mp.mpf(0)] * 4 for _ in range(4)] for _ in range(4)]
    ddg = [[[[mp.mpf(0)] * 4 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for (i, j), fn in comps.items():
        for k in range(4):
            od = tuple(1 if q == k else 0 for q in range(4))
            dg[k][i][j] = dg[k][j][i] = mp.diff(fn, X0, od)
            for l in range(k, 4):
                od2 = tuple((1 if q == k else 0) + (1 if q == l else 0) for q in range(4))
                ddg[k][l][i][j] = ddg[k][l][j][i] = ddg[l][k][i][j] = ddg[l][k][j][i] = mp.diff(fn, X0, od2)
    gi = g**-1
    dgi = [-(gi * mp.matrix([[dg[k][i][j] for j in range(4)] for i in range(4)]) * gi) for k in range(4)]
    G = [[[sum(gi[m, n] * (dg[i][n][j] + dg[j][n][i] - dg[n][i][j]) for n in range(4)) / 2 for j in range(4)] for i in range(4)] for m in range(4)]
    dG = [[[[sum(dgi[k][m, n] * (dg[i][n][j] + dg[j][n][i] - dg[n][i][j]) +
                 gi[m, n] * (ddg[k][i][n][j] + ddg[k][j][n][i] - ddg[k][n][i][j]) for n in range(4)) / 2
             for j in range(4)] for i in range(4)] for m in range(4)] for k in range(4)]
    R = mp.zeros(4, 4)
    for i in range(4):
        for j in range(4):
            tot = mp.mpf(0)
            for m in range(4):
                tot += dG[m][m][i][j] - dG[j][m][i][m]
                for n in range(4):
                    tot += G[m][m][n] * G[n][i][j] - G[m][j][n] * G[n][i][m]
            R[i, j] = tot
    fv = prof(r_)
    E = [mp.matrix([1, v0 * fv, 0, 0]), mp.matrix([0, 1, 0, 0]), mp.matrix([0, 0, 1, 0]), mp.matrix([0, 0, 0, 1])]
    return mp.matrix([[(E[a_].T * R * E[b_])[0] for b_ in range(4)] for a_ in range(4)])


def lab_ricci(h):
    """The tanh points of tanh_points, and a generic point, from the lab metric: the frame Ricci block equals the comoving-jet block of
    the library, R(e_phi, .) = 0, and the rank test on the lab block alone gives a triple root of rank 2. Planted: v = 0 is flat."""
    ok, rows = True, []
    mp.mp.dps = 70
    flat = lab_frame_ricci(mp.mpf(8), mp.mpf(1), mp.mpf(0), mp.mpf("0.95"), mp.mpf("0.6"), mp.mpf("0.8"))
    okf = max(abs(flat[i, j]) for i in range(4) for j in range(4)) < mp.mpf(10)**-60
    ok = ok and okf
    rows.append("planted: v = 0 gives max |R| = %s" % mp.nstr(max(abs(flat[i, j]) for i in range(4) for j in range(4)), 2))
    cases = [("generic", mp.mpf(8), mp.mpf(1), mp.mpf("0.7"), mp.mpf("0.9"), mp.mpf("0.6"))]
    for name in ("Le", "Rodal"):
        sig, Rb, v0 = WALLS[name]
        for pnt in points_cached(name)[0]:
            cases.append((name, sig, Rb, v0, pnt["r"], mp.sqrt(pnt["w"])))
    for name, sig, Rb, v0, r_, c_ in cases:
        s_ = mp.sqrt(1 - c_**2)
        Rl4 = lab_frame_ricci(sig, Rb, v0, r_, c_, s_)
        F_, a_, b_, _ = jets_mp(r_, sig, Rb)
        Rb3 = block_num((v0, F_, a_, b_, 1 / r_, c_, s_))
        scale = max(abs(Rl4[i, j]) for i in range(4) for j in range(4))
        dif = max(abs(Rl4[i, j] - Rb3[i][j]) for i in range(3) for j in range(3)) / scale
        dphi = max(abs(Rl4[3, j]) for j in range(4)) / scale
        good = dif < mp.mpf(10)**-25 and dphi < mp.mpf(10)**-25
        info = ""
        if name != "generic":
            M_ = mixed_at([[Rl4[i, j] for j in range(3)] for i in range(3)])
            sv, n3 = numeric_type3(M_, 70)
            good = good and sv[2] < mp.mpf(10)**-20 * max(1, scale) and sv[1] > mp.mpf(10)**-6 and n3 < mp.mpf(10)**-18 * max(1, scale)**3
            info = ", lab block alone: singular values %s, %s, %s, |N^3| = %s" % (mp.nstr(sv[0], 5), mp.nstr(sv[1], 5), mp.nstr(sv[2], 2), mp.nstr(n3, 2))
        ok = ok and good
        rows.append("%s r = %s, cos = %s: max |lab - comoving jets| / max |R| = %s, R(e_phi, .) = %s%s" % (
            name, mp.nstr(r_, 8), mp.nstr(c_, 8), mp.nstr(dif, 2), mp.nstr(dphi, 2), info))
    h.check("lab_ricci", ok, "; ".join(rows))


# ----------------------------------------------------------------------------------------------------------------------
# harness
# ----------------------------------------------------------------------------------------------------------------------
CHECKS = {
    "planted_types": planted_types,
    "polynomials": polynomials,
    "cert_rank": cert_rank,
    "cert_T0": cert_T0,
    "cert_disc": cert_disc,
    "equator_axis": equator_axis,
    "ray_family": ray_family,
    "exact_points": exact_points,
    "type2_failures": type2_failures,
    "description": description,
    "tanh_jets": tanh_jets_check,
    "tanh_planted": tanh_planted,
    "tanh_points": tanh_points,
    "tanh_curve": tanh_curve,
    "tanh_enclosure": tanh_enclosure,
    "lab_ricci": lab_ricci,
}

CONTROLS = {
    "no-rank": (["planted_types", "exact_points", "description"], "every D = 0 tensor is called Type II"),
    "disc-sign": (["planted_types", "type2_failures", "description"], "the sign of D flipped in the type tests"),
    "E2-typo": (["planted_types", "polynomials", "cert_T0", "cert_disc", "ray_family", "exact_points"], "E2 = 27 c0 + c2^3"),
    "E1-typo": (["planted_types", "polynomials", "cert_T0", "cert_disc", "ray_family", "exact_points"], "E1 = c2^2 - 4 c1"),
    "table-typo": (["polynomials"], "the sign of R_0x in the typed table"),
    "weight-typo": (["polynomials"], "the weight of E1 taken as 5"),
    "cof-rank-perturb": (["cert_rank"], "the first cofactor of the rank certificate plus 1"),
    "cof-rank-drop": (["cert_rank"], "the cofactor of m02,12 set to 0"),
    "cof-g1-perturb": (["cert_T0"], "g1 plus 1"),
    "cof-g2-drop": (["cert_T0"], "g2 set to 0"),
    "cof-B-perturb": (["cert_T0"], "B_0 plus 1"),
    "theta-typo": (["cert_T0"], "Theta with (1 + c) for (1 + c^2)"),
    "disc27": (["cert_disc"], "26 D + E2^2"),
    "cof-disc-perturb": (["cert_disc"], "G plus c2"),
    "axis-typo": (["equator_axis"], "the vanishing tested off the equator and the axis"),
    "ray-typo": (["ray_family", "exact_points"], "n_X of the null ray plus 1"),
    "nf-radicand": (["exact_points", "description"], "the radicand of Q(sqrt d) plus 1"),
    "phi-off": (["exact_points", "description"], "phi shifted by 1/1000"),
    "loosen-E1": (["type2_failures", "description"], "the description without E1 = 0"),
    "loosen-E2": (["type2_failures", "description"], "the description without E2 = 0"),
    "loosen-T": (["type2_failures", "description"], "the description without T != 0"),
    "jets-sign": (["tanh_jets"], "the sign of f'' in the own jets"),
    "P-term": (["tanh_planted"], "the first term of P dropped"),
    "T-filter-off": (["tanh_planted"], "the T != 0 test skipped"),
    "tanh-off-root": (["tanh_points"], "r shifted by 1e-6 before the rank test"),
    "enclosure-shift": (["tanh_enclosure"], "the Krawczyk box moved off the root"),
    "eq-end-typo": (["tanh_curve"], "the equator end from r a + 3F = 0"),
    "lab-sign": (["lab_ricci"], "the sign of the shift in the lab metric"),
}


def main():
    argv = sys.argv[1:]
    if "--search-seeds" in argv:
        found, tried, nrat = search_seeds()
        for sd, root in found:
            print([str(z) for z in sd], root)
        say("%d seeds tried, %d squarefree d, %d points with rational phi" % (tried, len(found), nrat))
        return
    import platform
    say("he_type3.py: python %s, sympy %s, mpmath %s" % (platform.python_version(), sp.__version__, mp.__version__))
    Hn = L.Harness(CHECKS, CONTROLS, "he_type3")
    sel = next((x.split("=", 1)[1] for x in argv if x.startswith("--checks=")), "")
    if sel:
        set_mut("")
        Hn.run(sel.split(","))
        say("checks failed: %s; wall %.1f s" % (Hn.failed or "none", time.time() - T0))
        raise SystemExit(1 if Hn.failed else 0)
    try:
        Hn.main(argv, set_mut)
    finally:
        say("wall %.1f s" % (time.time() - T0))


if __name__ == "__main__":
    main()

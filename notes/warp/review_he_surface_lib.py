"""Shared code of the independent review of #72 (issue #101): own Ricci tensor, own Hawking-Ellis type test.

Written from the metric and from the claims of he_surface.md and he_type.md. It does not import, and was written without
reading, he_surface.py, he_type.py, review_he_derive.py or phantom_source.py.

Geometry. Comoving coordinates (t, x, y, z) with x -> x - v t turn ds^2 = -dt^2 + (dx - v f dt)^2 + dy^2 + dz^2 into
-dt^2 + (dx + w dt)^2 + dy^2 + dz^2, w = v (1 - f(r)), r = |(x, y, z)|, stationary. At the point (x, y, z) = r (c, s, 0)
the Ricci tensor depends only on the jet (w, w_i, w_ij), which `ricci_jet` takes as independent symbols. The orthonormal
frame at t = 0 is (n, d_x, d_y, d_z) with n = d_t - w d_x, and e_phi = d_z, e_rho = d_y at z = 0.

Type test. `classify` takes a mixed tensor M = eta R, finds the characteristic polynomial, its repeated real roots, the
rank of M - lambda by Gaussian elimination in Fractions (own routine, not sympy's rank) and the size of the largest Jordan
block, and returns I, II, III or IV.
"""
import itertools
from fractions import Fraction as Fr

import sympy as sp

MUT = ""  # the control in force, set by the scripts

v, F, a, b, u, c, s, lam, x = sp.symbols("v F a b u c s lam x")
ETA3 = (-1, 1, 1)


# ----------------------------------------------------------------------------------------------------------------------
# Ricci tensor from the jet of the shift
# ----------------------------------------------------------------------------------------------------------------------
def ricci(g, gi, X):
    """R_ab = d_c G^c_ab - d_b G^c_ac + G^c_cd G^d_ab - G^c_bd G^d_ac, G the Levi-Civita connection of g."""
    n = len(X)
    dg = [[[sp.diff(g[i, j], X[k]) for k in range(n)] for j in range(n)] for i in range(n)]
    G = [[[sum(gi[k, l] * (dg[l][i][j] + dg[l][j][i] - dg[i][j][l]) for l in range(n)) / 2
           for j in range(n)] for i in range(n)] for k in range(n)]
    R = sp.zeros(n, n)
    for p in range(n):
        for q in range(p, n):
            e = 0
            for k in range(n):
                e += sp.diff(G[k][p][q], X[k]) - sp.diff(G[k][p][k], X[q])
                for d in range(n):
                    e += G[k][k][d] * G[d][p][q] - G[k][q][d] * G[d][p][k]
            R[p, q] = R[q, p] = e
    return R


_CACHE = {}


def ricci_jet():
    """Frame Ricci (n, d_x, d_y, d_z) of the shift w with independent jet symbols, and the substitution dict
    that puts the jet of w = v (1 - f(r)) at the point r (c, s, 0)."""
    if "rj" in _CACHE:
        return _CACHE["rj"]
    t, X1, X2, X3 = X = sp.symbols("t X1 X2 X3")
    W0 = sp.Symbol("W0")
    w1 = sp.symbols("w1 w2 w3")
    wij = {(i, j): sp.Symbol("w%d%d" % (i + 1, j + 1)) for i in range(3) for j in range(i, 3)}
    xs = [X1, X2, X3]
    W = W0 + sum(w1[i] * xs[i] for i in range(3)) + sp.Rational(1, 2) * sum(
        wij[min(i, j), max(i, j)] * xs[i] * xs[j] for i in range(3) for j in range(3))
    g = sp.Matrix([[-1 + W**2, W, 0, 0], [W, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    gi = sp.Matrix([[-1, W, 0, 0], [W, 1 - W**2, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    assert (g * gi).applyfunc(sp.expand) == sp.eye(4)
    R = ricci(g, gi, X).subs({t: 0, X1: 0, X2: 0, X3: 0}).applyfunc(sp.expand)
    E = [sp.Matrix([1, -W0, 0, 0])] + [sp.Matrix([1 if k == i else 0 for k in range(4)]) for i in (1, 2, 3)]
    Rf = sp.Matrix(4, 4, lambda i, j: sp.expand((E[i].T * R * E[j])[0, 0]))
    n = [c, s, 0]
    sub = {W0: v * F}
    for i in range(3):
        sub[w1[i]] = -v * a * n[i]
    for (i, j), sym in wij.items():
        sub[sym] = -v * (b * n[i] * n[j] + a * u * ((1 if i == j else 0) - n[i] * n[j]))
    out = Rf.subs(sub).applyfunc(sp.expand)
    _CACHE["rj"] = (out, Rf, sub)
    return _CACHE["rj"]


def red(e):
    """Reduce modulo s^2 = 1 - c^2."""
    e = sp.expand(e)
    if not e.has(s):
        return e
    return sp.expand(sp.rem(sp.Poly(e, s), sp.Poly(s**2 - (1 - c**2), s)).as_expr())


def frame_block():
    """The 3x3 symmetric R on (n, d_x, d_y) with entries reduced modulo s^2 + c^2 = 1, and the e_phi row."""
    out, _, _ = ricci_jet()
    R3 = sp.Matrix(3, 3, lambda i, j: red(out[i, j]))
    return R3, [red(out[i, 3]) for i in range(4)]


def mixed(R3):
    return sp.Matrix(3, 3, lambda i, j: ETA3[i] * R3[i, j])


def sigma(R3, l=lam):
    """Sigma = R - lam eta, so that M - lam = eta Sigma."""
    return sp.Matrix(3, 3, lambda i, j: R3[i, j] - (l * ETA3[i] if i == j else 0))


PAIRS = [(0, 1), (0, 2), (1, 2)]
MINOR_NAMES = ["m01,01", "m01,02", "m01,12", "m02,02", "m02,12", "m12,12"]
MINOR_IDX = [(p, q) for i, p in enumerate(PAIRS) for q in PAIRS[i:]]


def minors(S):
    """The six distinct 2x2 minors of a symmetric 3x3: rows (r1, r2), columns (c1, c2)."""
    return [S[p[0], q[0]] * S[p[1], q[1]] - S[p[0], q[1]] * S[p[1], q[0]] for p, q in MINOR_IDX]


# ----------------------------------------------------------------------------------------------------------------------
# exact linear algebra in Fractions
# ----------------------------------------------------------------------------------------------------------------------
def rank_frac(M):
    """Rank of a matrix of Fractions by Gaussian elimination with exact pivots."""
    A = [[Fr(e) for e in row] for row in M]
    n, m = len(A), len(A[0])
    r = 0
    for col in range(m):
        piv = next((i for i in range(r, n) if A[i][col] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        for i in range(r + 1, n):
            if A[i][col] != 0:
                f = A[i][col] / A[r][col]
                A[i] = [x1 - f * x2 for x1, x2 in zip(A[i], A[r])]
        r += 1
    return r


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def matsub_scalar(M, lm):
    return [[M[i][j] - (lm if i == j else 0) for j in range(len(M))] for i in range(len(M))]


def eye(n):
    return [[Fr(1 if i == j else 0) for j in range(n)] for i in range(n)]


def matpow(M, k):
    out = eye(len(M))
    for _ in range(k):
        out = matmul(out, M)
    return out


def is_zero(M):
    return all(e == 0 for row in M for e in row)


def charpoly_frac(M):
    """Coefficients [1, c_{n-1}, ..., c_0] of det(x - M) by Faddeev-LeVerrier, exactly."""
    n = len(M)
    Mk = eye(n)
    cs = [Fr(1)]
    for k in range(1, n + 1):
        AM = matmul(M, Mk)
        ck = -sum(AM[i][i] for i in range(n)) / k
        cs.append(ck)
        Mk = [[AM[i][j] + (ck if i == j else 0) for j in range(n)] for i in range(n)]
    return cs


def classify(M):
    """Hawking-Ellis type of the mixed tensor M (list of lists of Fractions, eta-self-adjoint).

    Returns dict(type, D (for n = 3), lam, m, rank, jordan) ; lam is the repeated real root.
    Roots of the characteristic polynomial: the sqf decomposition gives each repeated root as a rational linear factor
    when the entries are rational.
    """
    n = len(M)
    cs = charpoly_frac(M)
    X = sp.Symbol("X")
    chi = sp.Poly(sum(sp.Rational(ci.numerator, ci.denominator) * X ** (n - i) for i, ci in enumerate(cs)), X)
    out = {"type": None, "D": None, "lam": None, "m": 1, "rank": None, "jordan": None}
    if n == 3:
        c2, c1, c0 = cs[1], cs[2], cs[3]
        D = 18 * c2 * c1 * c0 - 4 * c2**3 * c0 + c2**2 * c1**2 - 4 * c1**3 - 27 * c0**2
        out["D"] = D
    sq = sp.factor_list(chi.as_expr(), X)[1]
    sq = [(sp.Poly(f, X), k) for f, k in sq]
    rep = [(f, k) for f, k in sq if k >= 2]
    if not rep:
        nreal = chi.count_roots()
        allreal = nreal == n
        if MUT == "disc-sign":
            allreal = not allreal
        out["type"] = "I" if allreal else "IV"
        return out
    ok_rank = {}
    best = None
    jmax = 1
    for f, k in rep:
        if f.degree() != 1:
            raise NotImplementedError("repeated non-rational root")
        lm = Fr(int(-f.all_coeffs()[1].p), int(f.all_coeffs()[1].q)) / Fr(int(f.all_coeffs()[0].p), int(f.all_coeffs()[0].q))
        mult = k
        N = matsub_scalar(M, lm)
        rk = rank_frac(N)
        # largest Jordan block: smallest p with rank N^p = n - mult
        p, Np = 1, N
        while rank_frac(Np) != n - mult:
            Np = matmul(Np, N)
            p += 1
        ok = (rk == n - mult)
        if MUT == "no-rank":
            ok = False
        if MUT == "mult-off":
            ok = rk == n - mult + 1
        ok_rank[lm] = ok
        jmax = max(jmax, p)
        best = (lm, mult, rk)
        out.update(lam=lm, m=mult, rank=rk, jordan=p)
    if all(ok_rank.values()):
        # all repeated roots diagonalizable. Type I iff the remaining roots are real.
        rest = chi
        for f, k in rep:
            rest = sp.quo(rest, f**k)
        out["type"] = "I" if rest.count_roots() == rest.degree() else "IV"
    else:
        if MUT == "no-rank":
            out["type"] = "II"
        else:
            out["type"] = "III" if jmax >= 3 else "II"
    return out


# ----------------------------------------------------------------------------------------------------------------------
# evaluation of polynomials at exact points, in Fractions or in a number field
# ----------------------------------------------------------------------------------------------------------------------
VARS7 = (v, F, a, b, u, c, s)


def poly_terms(expr):
    P = sp.Poly(sp.expand(expr), *VARS7)
    return [(tuple(mn), Fr(int(cf.p), int(cf.q))) for mn, cf in P.terms()]


def peval(terms, pt):
    tot = 0
    for mn, cf in terms:
        t = cf
        for xx, e in zip(pt, mn):
            if e:
                t = t * xx**e
        tot = tot + t
    return tot


class QS:
    """Element p + q sqrt(d) of Q(sqrt d), Fractions p, q. The radicand is the class attribute d (controls change it)."""
    d = 113

    def __init__(self, p, q=0):
        self.p, self.q = Fr(p), Fr(q)

    def __add__(self, o):
        o = o if isinstance(o, QS) else QS(o)
        return QS(self.p + o.p, self.q + o.q)

    __radd__ = __add__

    def __neg__(self):
        return QS(-self.p, -self.q)

    def __sub__(self, o):
        return self + (-(o if isinstance(o, QS) else QS(o)))

    def __rsub__(self, o):
        return QS(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, QS) else QS(o)
        return QS(self.p * o.p + QS.d * self.q * o.q, self.p * o.q + self.q * o.p)

    __rmul__ = __mul__

    def inv(self):
        den = self.p**2 - QS.d * self.q**2
        return QS(self.p / den, -self.q / den)

    def __truediv__(self, o):
        o = o if isinstance(o, QS) else QS(o)
        return self * o.inv()

    def __rtruediv__(self, o):
        return QS(o) * self.inv()

    def __pow__(self, e):
        out = QS(1)
        for _ in range(e):
            out = out * self
        return out

    def __eq__(self, o):
        o = o if isinstance(o, QS) else QS(o)
        return self.p == o.p and self.q == o.q

    def __hash__(self):
        return hash((self.p, self.q))

    def iszero(self):
        return self.p == 0 and self.q == 0

    def __float__(self):
        return float(self.p) + float(self.q) * QS.d**0.5

    def __repr__(self):
        return "(%s + %s*sqrt(%d))" % (self.p, self.q, QS.d)


def qs_matmul(A, B):
    n = len(A)
    return [[sum((A[i][k] * B[k][j] for k in range(n)), QS(0)) for j in range(n)] for i in range(n)]


def qs_rank(M):
    """Rank over Q(sqrt d) by elimination with exact zero tests."""
    A = [[e if isinstance(e, QS) else QS(e) for e in row] for row in M]
    n, m = len(A), len(A[0])
    r = 0
    for col in range(m):
        piv = next((i for i in range(r, n) if not A[i][col].iszero()), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        for i in range(r + 1, n):
            if not A[i][col].iszero():
                f = A[i][col] / A[r][col]
                A[i] = [x1 - f * x2 for x1, x2 in zip(A[i], A[r])]
        r += 1
    return r


# ----------------------------------------------------------------------------------------------------------------------
# harness
# ----------------------------------------------------------------------------------------------------------------------
class Harness:
    def __init__(self, checks, controls, title):
        """checks: dict name -> function; controls: dict name -> (list of target check names, description)."""
        self.checks, self.controls, self.title = checks, controls, title
        self.failed = []

    def check(self, name, ok, msg):
        print("[%s] %s: %s" % ("PASS" if ok else "FAIL", name, msg), flush=True)
        if not ok:
            self.failed.append(name)

    def run(self, names):
        for nm in names:
            n0 = len(self.failed)
            try:
                self.checks[nm](self)
            except Exception as ex:  # an exception is a failure of the check
                import traceback
                traceback.print_exc()
                self.check(nm, False, "raised %r" % (ex,))
            if len(self.failed) == n0 and False:
                pass

    def main(self, argv, set_mut):
        mut = next((x.split("=", 1)[1] for x in argv if x.startswith("--mutate=")), "")
        if mut:
            targets, desc = self.controls[mut]
            set_mut(mut)
            print("control %s: %s; targets %s" % (mut, desc, targets), flush=True)
            self.run(targets)
            fired = [t for t in targets if t in self.failed]
            print("control %s: failed %s of targets %s -> %s" % (
                mut, self.failed, targets, "FIRED" if fired == targets else "SILENT or partial"), flush=True)
            raise SystemExit(1 if fired == targets else 0)
        set_mut("")
        self.run(list(self.checks))
        base = list(self.failed)
        if "--no-controls" in argv:
            print("checks failed: %s" % (base or "none"), flush=True)
            raise SystemExit(1 if base else 0)
        silent = []
        for cn, (targets, desc) in self.controls.items():
            self.failed = []
            set_mut(cn)
            import io, contextlib
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                self.run(targets)
            fired = [t for t in targets if t in self.failed]
            ok = fired == targets
            print("control %-14s targets %-40s failed %-40s %s" % (cn, targets, self.failed, "FIRED" if ok else "SILENT"),
                  flush=True)
            if not ok:
                silent.append(cn)
        set_mut("")
        print("%d checks, %d failed %s; %d controls, %d silent %s" % (
            len(self.checks), len(base), base, len(self.controls), len(silent), silent), flush=True)
        raise SystemExit(1 if base or silent else 0)

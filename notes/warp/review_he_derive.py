#!/usr/bin/env python3
"""Independent derivation of #64's Hawking-Ellis closed form (issue #71).

Re-derives the four identities of notes/warp/he_type.md from the metric alone:

  T1  D = v^6 (p6 + p8 v^2 + p10 v^4 + p12 v^6) with p6 v^6 = -4|j|^6
  T2  p12 >= 0
  T3  axis:     D = -4 a^2 v^6 Q_ax^2 / r^6
  T4  equator:  D = v^6 [r^2 v^2 a^4 - (a + r b)^2] Q_eq^2 / (16 r^6), Type I iff v r a^2 > |a + r b|

with a = f'(r), b = f''(r), F = 1 - f(r), c = cos(theta), s = sin(theta), u = 1/r.

The metric is ds^2 = -dt^2 + (dx + w dt)^2 + dy^2 + dz^2 with w = v (1 - f(r)), the comoving form of the Alcubierre
metric. D is the discriminant of the cubic of the mixed Ricci tensor R^a_b on its 3-block (n, e_x, e_rho), and is
unchanged by any frame, by the sign convention of Ricci and by the sign of v.

Two routes share no code with each other, nor with notes/warp/he_type.py or phantom_source.py (never imported, see
check `independence`):

  route K  coordinate Ricci tensor of a general metric (own Christoffel code) in Cartesian (t, x, y, z), at the point
           (r c, 0, r s). The frame is read off from the coordinate components.
  route C  Cartan's structure equations (own exterior calculus) in the cylindrical coframe
           (dt, dx + w dt, d rho, rho d phi).

Everything is exact (sympy rationals, sparse polynomial rings, python Fractions, and mpmath for the numerical
eigenvectors of `type_direct`). Every check can fail, and `--mutate=NAME` makes each control break its target.

    python3 notes/warp/review_he_derive.py                 # checks, then every control as a subprocess
    python3 notes/warp/review_he_derive.py --no-controls   # checks alone
    python3 notes/warp/review_he_derive.py --mutate=NAME   # one control; exit 1 means it fired
"""
import argparse
import random
import subprocess
import sys
import time
from fractions import Fraction as Fr

import mpmath as mp
import sympy as sp
from sympy import Function, Matrix, QQ, diff, symbols, zeros
from sympy.polys.rings import ring

MUT = set()


def mut(name):
    return name in MUT


# ------------------------------------------------------------------------------------------------------------------
# controls: name -> (what it changes, checks that must fail)
CONTROLS = {
    'eng-quad': ('drop the Gamma*Gamma terms (route K) and the omega^omega term (route C)', ['engines']),
    'k-drop-w2': ("route K: g_tt = -1, the shift's square dropped", ['alc_energy', 'routes_agree', 'T1_structure',
                                                                      'T3_axis', 'T4_equator']),
    'c-no-nshift': ('route C: the frame derivative E_0 = d_t - w d_x loses its shift, E_0 = d_t',
                    ['routes_agree', 'T1_structure', 'T3_axis', 'T4_equator']),
    'jets-b': ("the radial jets use f'' x 1.001", ['jets', 'T1_structure', 'T2_p12', 'T3_axis', 'T4_equator']),
    'disc-formula': ('the discriminant formula uses 17 c2 c1 c0 for 18 c2 c1 c0 (c0 = 0 on the axis and at lowest '
                     'order in v, so T1 and T3 do not see it)', ['disc_tools', 'T2_p12', 'T4_equator', 'type_direct']),
    'disc-c1': ('the discriminant formula uses 3 c1^3 for 4 c1^3',
                ['disc_tools', 'T1_structure', 'T3_axis', 'T4_equator', 'type_direct']),
    'lemma-diag': ('classification of every D = 0 tensor skips the rank test and answers Type II',
                   ['planted_types', 'Qeq0_and_edges']),
    'crit-factor': ('the equatorial criterion becomes v r a^2 > 2 |a + r b|', ['T4_equator']),
    'tgt-p6': ('the target p6 = -J^6/(16 r^6) becomes -J^6/(8 r^6)', ['T1_structure']),
    'tgt-p12': ('the target Q2 loses its -2 r a^2 term', ['T2_p12']),
    'tgt-axis': ('the target Q_ax = a^2 + r^2 v^2 (F b - a^2)^2 becomes a^2 + r^2 v^2 (F b + a^2)^2', ['T3_axis']),
    'tgt-eq': ('the target Q_eq = 4 v^2 a^2 F (F + r a) + (a + r b)^2 becomes 3 v^2 a^2 F (F + r a) + (a + r b)^2',
               ['T4_equator']),
}

# ------------------------------------------------------------------------------------------------------------------
# polynomial ring: v, F, a, b, u = 1/r, c = cos, s = sin, with s^2 -> 1 - c^2
R, v_, F_, a_, b_, u_, c_, s_ = ring('v,F,a,b,u,c,s', QQ)
IV, IC, IS = 0, 5, 6
v, F, a, b, u, c, s = sp.symbols('v F a b u c s')
SYMS = (v, F, a, b, u, c, s)


class NotPoly(Exception):
    pass


def red(p):
    """Reduce s^2 -> 1 - c^2 so that every polynomial has a unique form."""
    groups = {}
    for mon, co in p.terms():
        q, e = divmod(mon[IS], 2)
        groups.setdefault(q, {})[mon[:IS] + (e,)] = co
    out = R.zero
    for q, d in groups.items():
        out += R.from_dict(d) * (1 - c_ ** 2) ** q
    return out


def to_ring(expr):
    """sympy expression -> reduced ring element, or NotPoly when it is not a polynomial."""
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    if den.free_symbols:
        raise NotPoly(str(den))
    try:
        return red(R.from_expr(sp.expand(num)) * (QQ(1) / QQ(int(den.p), int(den.q)) if den != 1 else 1))
    except Exception as exc:
        raise NotPoly(str(exc))


def req(p, q):
    return red(p - q) == 0


def to_fr(x):
    x = sp.Rational(x)
    return Fr(int(x.p), int(x.q))


def at_point(expr, pt):
    """Exact value of a sympy expression in SYMS at the point pt (a dict symbol -> Fraction)."""
    sub = {k: sp.Rational(x.numerator, x.denominator) for k, x in pt.items()}
    return to_fr(sp.sympify(expr).xreplace(sub))


# ------------------------------------------------------------------------------------------------------------------
# engine 1: coordinate Ricci tensor of a metric
def coord_ricci(g, X):
    n = len(X)
    gi = g.inv().applyfunc(sp.simplify)
    Gam = [[[sum(gi[l, k] * (diff(g[k, m], X[nu]) + diff(g[k, nu], X[m]) - diff(g[m, nu], X[k])) for k in range(n)) / 2
             for nu in range(n)] for m in range(n)] for l in range(n)]
    Ric = zeros(n, n)
    for m in range(n):
        for nu in range(n):
            e = sum(diff(Gam[l][m][nu], X[l]) for l in range(n)) - sum(diff(Gam[l][m][l], X[nu]) for l in range(n))
            if not mut('eng-quad'):
                e += sum(Gam[l][l][k] * Gam[k][m][nu] for l in range(n) for k in range(n))
                e -= sum(Gam[l][nu][k] * Gam[k][m][l] for l in range(n) for k in range(n))
            Ric[m, nu] = e
    return Ric, gi


# engine 2: Cartan's structure equations for an orthonormal coframe
def cartan_ricci(e, X, eta):
    """e[a, mu] = theta^a_mu. Returns the frame Ricci tensor R_ab. Solves d theta^a = -omega^a_b ^ theta^b for omega
    (antisymmetric once lowered) and curvature Omega = d omega + omega ^ omega, and checks the first equation."""
    n = len(X)
    E = e.inv()

    def Ed(k, f):
        if mut('c-no-nshift') and k == 0:
            return diff(f, X[0])
        return sum(E[mu, k] * diff(f, X[mu]) for mu in range(n))

    # d theta^a = 1/2 Xs^a_{bc} theta^b ^ theta^c
    Xs = [[[sum((diff(e[k, nu], X[mu]) - diff(e[k, mu], X[nu])) * E[mu, p] * E[nu, q] for mu in range(n)
                for nu in range(n)) for q in range(n)] for p in range(n)] for k in range(n)]
    Xl = [[[eta[k, k] * Xs[k][p][q] for q in range(n)] for p in range(n)] for k in range(n)]
    om = [[[eta[k, k] * (Xl[k][p][q] + Xl[p][q][k] - Xl[q][k][p]) / 2 for q in range(n)] for p in range(n)]
          for k in range(n)]  # omega^k_{pq}, the 1-form omega^k_p = omega^k_{pq} theta^q
    for k in range(n):
        for p in range(n):
            for q in range(n):
                assert sp.simplify(Xs[k][p][q] - (om[k][p][q] - om[k][q][p])) == 0, 'structure equation'
                assert sp.simplify(eta[k, k] * om[k][p][q] + eta[p, p] * om[p][k][q]) == 0, 'omega antisymmetry'
    Rm = [[[[0] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for k in range(n):
        for p in range(n):
            for d in range(n):
                for q in range(d + 1, n):
                    t1 = Ed(d, om[k][p][q]) - Ed(q, om[k][p][d]) + sum(om[k][p][m] * Xs[m][d][q] for m in range(n))
                    t2 = 0 if mut('eng-quad') else sum(om[k][m][d] * om[m][p][q] - om[k][m][q] * om[m][p][d]
                                                       for m in range(n))
                    Rm[k][p][d][q] = t1 + t2
                    Rm[k][p][q][d] = -(t1 + t2)
    Ric = zeros(n, n)
    for p in range(n):
        for q in range(n):
            Ric[p, q] = sum(Rm[k][p][k][q] for k in range(n))
    return Ric


def jet_subs(expr, func, jets):
    """Replace the derivatives of func by jet symbols and then by the values in jets (dict symbol -> expr)."""
    reps = {}
    for d in expr.atoms(sp.Derivative):
        reps[d] = sp.Symbol('w_' + ''.join(str(vv) * k for vv, k in d.variable_count))
    expr = expr.xreplace(reps).xreplace({func: sp.Symbol('w0')})
    return expr.xreplace(jets)


# ------------------------------------------------------------------------------------------------------------------
# radial jets of w = v (1 - f(r)) at the point (x, y, z) = (r c, 0, r s) and (x, rho) = (r c, r s)
def jets_cart():
    bj = b * sp.Rational(1001, 1000) if mut('jets-b') else b
    nv, nm = [c, 0, s], ['x', 'y', 'z']
    J = {'w0': v * F}
    for i in range(3):
        J['w_' + nm[i]] = -v * a * nv[i]
        for j in range(i, 3):
            J['w_' + nm[i] + nm[j]] = -v * (bj * nv[i] * nv[j] + a * u * ((1 if i == j else 0) - nv[i] * nv[j]))
    return {sp.Symbol(k): sp.sympify(val) for k, val in J.items()}


def jets_cyl():
    bj = b * sp.Rational(1001, 1000) if mut('jets-b') else b
    J = {'w0': v * F, 'w_x': -v * a * c, 'w_rho': -v * a * s, 'w_xx': -v * (bj * c ** 2 + a * u * s ** 2),
         'w_rhorho': -v * (bj * s ** 2 + a * u * c ** 2), 'w_xrho': -v * (bj - a * u) * c * s,
         'w_rhox': -v * (bj - a * u) * c * s}
    return {sp.Symbol(k): sp.sympify(val) for k, val in J.items()}


tt, xx, yy, zz, rho, phi = symbols('t x y z rho phi')


def build_K():
    """Frame Ricci (n, e_x, e_rho = d_z, e_phi = d_y) and coordinate mixed Ricci R^mu_nu of route K, in SYMS."""
    w = Function('w')(xx, yy, zz)
    gtt = -1 if mut('k-drop-w2') else -1 + w ** 2
    g = Matrix([[gtt, w, 0, 0], [w, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    Ric, gi = coord_ricci(g, [tt, xx, yy, zz])
    J = jets_cart()
    ev = lambda e: sp.cancel(sp.together(jet_subs(sp.sympify(e), w, J)))
    Rc = Ric.applyfunc(ev)
    gij = gi.applyfunc(ev)
    W = v * F
    vec = [Matrix([1, -W, 0, 0]), Matrix([0, 1, 0, 0]), Matrix([0, 0, 0, 1]), Matrix([0, 0, 1, 0])]
    Rf = Matrix(4, 4, lambda i, j: sp.cancel(sp.together((vec[i].T * Rc * vec[j])[0, 0])))
    Rmix = (gij * Rc).applyfunc(lambda e: sp.cancel(sp.together(e)))
    return Rf, Rmix


def build_C():
    """Frame Ricci (E_0 = n, E_1 = d_x, E_2 = e_rho, E_3 = e_phi) of route C, in SYMS."""
    w = Function('w')(xx, rho)
    e = Matrix([[1, 0, 0, 0], [w, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, rho]])
    Ric = cartan_ricci(e, [tt, xx, rho, phi], sp.diag(-1, 1, 1, 1))
    J = jets_cyl()
    ev = lambda e_: sp.cancel(sp.together(jet_subs(sp.sympify(e_), w, J).xreplace({rho: s / u})))
    return Ric.applyfunc(ev)


# ------------------------------------------------------------------------------------------------------------------
# numeric-type-generic algebra (ring elements and Fractions)
def disc_formula(c2, c1, c0):
    k = 17 if mut('disc-formula') else 18
    k1 = 3 if mut('disc-c1') else 4
    return k * c2 * c1 * c0 - 4 * c2 ** 3 * c0 + c2 ** 2 * c1 ** 2 - k1 * c1 ** 3 - 27 * c0 ** 2


def cubic(M):
    """Characteristic polynomial l^3 + c2 l^2 + c1 l + c0 of a 3x3 matrix."""
    tr = M[0][0] + M[1][1] + M[2][2]
    c1_ = sum(M[i][i] * M[j][j] - M[i][j] * M[j][i] for i, j in ((0, 1), (0, 2), (1, 2)))
    det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
           + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    return -tr, c1_, -det


def block(Rf, jzero=False):
    """M = eta R on (n, e_x, e_rho) with eta = diag(-1, 1, 1)."""
    eta = (-1, 1, 1)
    M = [[eta[i] * Rf[i][j] for j in range(3)] for i in range(3)]
    if jzero:
        for i in (1, 2):
            M[0][i] = M[0][i] * 0
            M[i][0] = M[i][0] * 0
    return M


def Dof(Rf, jzero=False):
    return disc_formula(*cubic(block(Rf, jzero)))


def closed_forms(ctx):
    """The targets of he_type.md, written as 16 r^6 p12 etc. with u = 1/r cleared (see review_he.md section 3)."""
    v_, F_, a_, b_, u_, c_, s_ = ctx
    Jt = s_ ** 2 * (a_ * u_ + b_) ** 2 + 4 * c_ ** 2 * a_ ** 2 * u_ ** 2            # u^2 J^2
    p6x16 = -Jt ** 3 * (2 if mut('tgt-p6') else 1)                                         # 16 p6
    Q1t = -c_ ** 2 * F_ * (a_ * u_ - b_) ** 2 - 4 * u_ * a_ ** 3                       # u^2 Q1
    inner = F_ * (b_ - a_ * u_) - (0 if mut('tgt-p12') else 2 * a_ ** 2)
    Q2t = a_ ** 2 * (F_ * u_ + a_) ** 2 + c_ ** 2 * F_ * (a_ * u_ + b_) * inner        # u^2 Q2
    p12x16 = s_ ** 4 * F_ ** 2 * Q1t ** 2 * Q2t                                          # 16 p12
    sg = 1 if mut('tgt-axis') else -1
    Qax = a_ ** 2 * u_ ** 2 + v_ ** 2 * (F_ * b_ + sg * a_ ** 2) ** 2                  # u^2 Q_ax
    Dax = -4 * a_ ** 2 * v_ ** 6 * u_ ** 2 * Qax ** 2                                    # D on the axis
    k4 = 3 if mut('tgt-eq') else 4
    Qeq = k4 * v_ ** 2 * a_ ** 2 * F_ * u_ * (F_ * u_ + a_) + (a_ * u_ + b_) ** 2     # u^2 Q_eq
    T = v_ ** 2 * a_ ** 4 - (a_ * u_ + b_) ** 2                                          # u^2 [r^2 v^2 a^4 - (a + r b)^2]
    Deq16 = v_ ** 6 * T * Qeq ** 2                                                       # 16 D on the equator
    return dict(Jt=Jt, p6x16=p6x16, Q1t=Q1t, Q2t=Q2t, p12x16=p12x16, Qax=Qax, Dax=Dax, Qeq=Qeq, T=T, Deq16=Deq16)


def coef_v(p, k):
    """Coefficient of v^k of a ring element, with v removed."""
    return R.from_dict({m[:IV] + (0,) + m[IV + 1:]: co for m, co in p.terms() if m[IV] == k})


def sub_axis(p):
    d = {}
    for m, co in red(p).terms():
        if m[IS] != 0:
            continue
        if m[IC] % 2:
            raise ValueError('odd power of c at s = 0')
        k = m[:IC] + (0,) + m[IC + 1:]
        d[k] = d.get(k, 0) + co
    return R.from_dict(d)


def sub_equator(p):
    d = {}
    for m, co in red(p).terms():
        if m[IC] != 0:
            continue
        k = m[:IS] + (0,)
        d[k] = d.get(k, 0) + co
    return R.from_dict(d)


# ------------------------------------------------------------------------------------------------------------------
# exact type classification of a 3x3 mixed block M = eta R (Fractions)
def spM(M):
    return sp.Matrix(3, 3, lambda i, j: sp.Rational(M[i][j].numerator, M[i][j].denominator))


ETA3 = sp.diag(-1, 1, 1)


def classify_lemma(M):
    """The lemma of he_type.md: the sign of D, then the rank test on D = 0."""
    c2, c1, c0 = cubic(M)
    D = disc_formula(c2, c1, c0)
    if D > 0:
        return 'I'
    if D < 0:
        return 'IV'
    if mut('lemma-diag'):
        return 'II'
    roots = sp.Poly([1, sp.Rational(c2.numerator, c2.denominator), sp.Rational(c1.numerator, c1.denominator),
                     sp.Rational(c0.numerator, c0.denominator)], sp.Symbol('l'), domain=QQ).ground_roots()
    if sum(roots.values()) != 3:
        return '?'
    A = spM(M)
    lam, m = max(roots.items(), key=lambda kv: kv[1])
    rk = (A - lam * sp.eye(3)).rank()
    if rk == 3 - m:
        return 'I'
    return 'III' if (m == 3 and rk == 2) else 'II'


def _psd(G):
    n = G.shape[0]
    import itertools
    for k in range(1, n + 1):
        for idx in itertools.combinations(range(n), k):
            if G.extract(list(idx), list(idx)).det() < 0:
                return False
    return True


def classify_direct_exact(M):
    """Definition of the Hawking-Ellis types for D = 0 (every eigenvalue rational): Type I iff a timelike
    eigenvector exists; else the Jordan structure of the null eigenvectors gives II or III. Uses no discriminant."""
    A = spM(M)
    c2, c1, c0 = cubic(M)
    roots = sp.Poly([1, sp.Rational(c2.numerator, c2.denominator), sp.Rational(c1.numerator, c1.denominator),
                     sp.Rational(c0.numerator, c0.denominator)], sp.Symbol('l'), domain=QQ).ground_roots()
    if sum(roots.values()) != 3:
        return '?'
    jordan = 1
    for lam in roots:
        B = A - lam * sp.eye(3)
        N = B.nullspace()
        Nm = sp.Matrix.hstack(*N)
        if not _psd(Nm.T * ETA3 * Nm):
            return 'I'
        dims = [len(N), len((B ** 2).nullspace()), len((B ** 3).nullspace())]
        if dims[2] > dims[1]:
            jordan = max(jordan, 3)
        elif dims[1] > dims[0]:
            jordan = max(jordan, 2)
    return {3: 'III', 2: 'II'}.get(jordan, 'I?')


def eig_eta_mp(M):
    """For D != 0: the eta-norms (eta(u, u) / |u|^2) of the real eigenvectors, by mpmath at 60 digits."""
    mp.mp.dps = 60
    c2, c1, c0 = cubic(M)
    q = lambda x: mp.mpf(x.numerator) / mp.mpf(x.denominator)
    try:
        roots = mp.polyroots([1, q(c2), q(c1), q(c0)], maxsteps=500, extraprec=300)
    except mp.libmp.libhyper.NoConvergence:
        return []
    A0 = mp.matrix(3, 3)
    for i in range(3):
        for j in range(3):
            A0[i, j] = q(M[i][j])
    out = []
    for lam in roots:
        if abs(mp.im(lam)) > mp.mpf(10) ** -40:
            continue
        lam = mp.re(lam)
        rows = [[A0[i, j] - (lam if i == j else 0) for j in range(3)] for i in range(3)]
        best = None
        for i, j in ((0, 1), (0, 2), (1, 2)):
            cr = [rows[i][1] * rows[j][2] - rows[i][2] * rows[j][1], rows[i][2] * rows[j][0] - rows[i][0] * rows[j][2],
                  rows[i][0] * rows[j][1] - rows[i][1] * rows[j][0]]
            nrm = sum(x * x for x in cr)
            if best is None or nrm > best[0]:
                best = (nrm, cr)
        nrm, cr = best
        out.append((-cr[0] ** 2 + cr[1] ** 2 + cr[2] ** 2) / nrm)
    return out


# ------------------------------------------------------------------------------------------------------------------
class Checks:
    def __init__(self):
        self.rows = []

    def add(self, name, ok, detail=''):
        self.rows.append((name, bool(ok), detail))
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)

    def failed(self):
        return sorted({n for n, ok, _ in self.rows if not ok})


def randpt(rng, c0=None, s0=None):
    """A random exact point of SYMS. (c0, s0) is a rational point of the unit circle, random unless given."""
    vv = Fr(rng.randint(1, 40), 8)
    FF = Fr(rng.randint(1, 16), 16)
    aa = Fr(rng.choice([-1, 1]) * rng.randint(1, 12), 4)
    bb = Fr(rng.choice([-1, 1]) * rng.randint(0, 12), 4)
    uu = Fr(rng.randint(1, 12), rng.randint(2, 12))
    if c0 is None:
        t = Fr(rng.randint(1, 5), rng.randint(1, 5))
        c0, s0 = (1 - t * t) / (1 + t * t), 2 * t / (1 + t * t)
    return {v: vv, F: FF, a: aa, b: bb, u: uu, c: c0, s: s0}


def frame_at(Rf, pt):
    """Exact frame Ricci at a point. Only the 3x3 block of (n, e_x, e_rho) is evaluated, the rest is 0."""
    return [[at_point(Rf[i, j], pt) if (i < 3 and j < 3) else Fr(0) for j in range(4)] for i in range(4)]


def ctx_of(pt):
    return tuple(pt[x] for x in SYMS)


def ring_frame(Rf):
    return [[to_ring(Rf[i, j]) for j in range(4)] for i in range(4)]


RING_CTX = (v_, F_, a_, b_, u_, c_, s_)


def fit_cubic_in_V(Rfs_by_v, pt_rest):
    """For fixed (F, a, b, u, c, s), D(v)/v^6 as a cubic in V = v^2 from 5 exact values of v. Returns
    (p6, p8, p10, p12, residual at the fifth point)."""
    vals = sorted(Rfs_by_v)
    E = [(vv * vv, Dof(Rfs_by_v[vv]) / vv ** 6) for vv in vals]
    xs, ys = [e[0] for e in E[:4]], [e[1] for e in E[:4]]
    X = sp.Symbol('V')
    poly = sp.interpolate([(sp.Rational(x_.numerator, x_.denominator), sp.Rational(y_.numerator, y_.denominator))
                           for x_, y_ in zip(xs, ys)], X)
    P = sp.Poly(poly, X)
    co = [to_fr(P.coeff_monomial(X ** k)) for k in range(4)]
    V5, y5 = E[4]
    res = sum(co[k] * V5 ** k for k in range(4)) - y5
    return co, res


# ------------------------------------------------------------------------------------------------------------------
def main_checks():
    ck = Checks()
    t0 = time.time()
    rng = random.Random(71)

    # ---- independence ------------------------------------------------------------------------------------------
    bad = [m for m in ('he_type', 'phantom_source') if m in sys.modules]
    ck.add('independence', not bad, 'imports: stdlib, sympy %s, mpmath %s; he_type and phantom_source are not '
           'imported' % (sp.__version__, mp.__version__))

    # ---- planted truth for the two engines ---------------------------------------------------------------------
    A = Function('A')(tt)
    Rk, _ = coord_ricci(sp.diag(-1, A ** 2, A ** 2, A ** 2), [tt, xx, yy, zz])
    Rcf = cartan_ricci(sp.diag(1, A, A, A), [tt, xx, yy, zz], sp.diag(-1, 1, 1, 1))
    Ad, Add = A.diff(tt), A.diff(tt, 2)
    frw_k = [-3 * Add / A, A * Add + 2 * Ad ** 2, A * Add + 2 * Ad ** 2, A * Add + 2 * Ad ** 2]
    frw_c = [-3 * Add / A] + [(A * Add + 2 * Ad ** 2) / A ** 2] * 3
    okK = all(sp.simplify(Rk[i, i] - frw_k[i]) == 0 for i in range(4)) and all(
        sp.simplify(Rk[i, j]) == 0 for i in range(4) for j in range(4) if i != j)
    okC = all(sp.simplify(Rcf[i, i] - frw_c[i]) == 0 for i in range(4)) and all(
        sp.simplify(Rcf[i, j]) == 0 for i in range(4) for j in range(4) if i != j)
    # Painleve-Gullstrand Schwarzschild: unit lapse, flat slices, the same family as the Alcubierre metric. R_ab = 0.
    r_, th, ph, k = sp.symbols('r theta phi k', positive=True)
    beta = k / sp.sqrt(r_)                                                      # k^2 = 2 M
    gp = Matrix([[-1 + beta ** 2, beta, 0, 0], [beta, 1, 0, 0], [0, 0, r_ ** 2, 0], [0, 0, 0, (r_ * sp.sin(th)) ** 2]])
    RpK, _ = coord_ricci(gp, [tt, r_, th, ph])
    ep = Matrix([[1, 0, 0, 0], [beta, 1, 0, 0], [0, 0, r_, 0], [0, 0, 0, r_ * sp.sin(th)]])
    RpC = cartan_ricci(ep, [tt, r_, th, ph], sp.diag(-1, 1, 1, 1))
    okPK = all(sp.simplify(RpK[i, j]) == 0 for i in range(4) for j in range(4))
    okPC = all(sp.simplify(RpC[i, j]) == 0 for i in range(4) for j in range(4))
    ck.add('engines', okK and okC and okPK and okPC,
           f'FRW R_tt = -3 A\'\'/A and R_ij = A A\'\' + 2 A\'^2 by route K {okK}, route C {okC}; '
           f'Painleve-Gullstrand Schwarzschild R_ab = 0 by route K {okPK}, route C {okPC}')

    # ---- jets --------------------------------------------------------------------------------------------------
    def jet_ok():
        res = []
        x_, y_, z_, rr = sp.symbols('x_ y_ z_ rr')
        profiles = [lambda r: 1 / (1 + r ** 2), lambda r: r ** -3 + 2 * r ** 2]
        for prof in profiles:
            for (X0, Z0) in ((3, 4), (5, 12), (-8, 15), (5, 0), (-5, 0)):
                r0 = sp.sqrt(X0 ** 2 + Z0 ** 2)
                vv = sp.Rational(7, 3)
                W = vv * (1 - prof(sp.sqrt(x_ ** 2 + y_ ** 2 + z_ ** 2)))
                f = prof(rr)
                pt = {v: to_fr(vv), F: to_fr((1 - f).subs(rr, r0)), a: to_fr(f.diff(rr).subs(rr, r0)),
                      b: to_fr(f.diff(rr, 2).subs(rr, r0)), u: to_fr(1 / r0), c: to_fr(sp.Rational(X0) / r0),
                      s: to_fr(sp.Rational(Z0) / r0)}
                J = jets_cart()
                nm = ['x', 'y', 'z']
                crd = [x_, y_, z_]
                pos = {x_: X0, y_: 0, z_: Z0}
                for i in range(3):
                    d1 = W.diff(crd[i]).subs(pos)
                    res.append(at_point(J[sp.Symbol('w_' + nm[i])], pt) == to_fr(sp.simplify(d1)))
                    for j in range(i, 3):
                        d2 = W.diff(crd[i], crd[j]).subs(pos)
                        res.append(at_point(J[sp.Symbol('w_' + nm[i] + nm[j])], pt) == to_fr(sp.simplify(d2)))
                # cylindrical jets: w(x, rho) = v (1 - f(sqrt(x^2 + rho^2))) at (X0, Z0)
                Wc = vv * (1 - prof(sp.sqrt(x_ ** 2 + z_ ** 2)))
                Jc = jets_cyl()
                for nmj, (cs) in (('w_x', (x_,)), ('w_rho', (z_,)), ('w_xx', (x_, x_)), ('w_rhorho', (z_, z_)),
                                  ('w_xrho', (x_, z_))):
                    dd = Wc.diff(*cs).subs({x_: X0, z_: Z0})
                    res.append(at_point(Jc[sp.Symbol(nmj)], pt) == to_fr(sp.simplify(dd)))
        return all(res), len(res)
    okj, nj = jet_ok()
    ck.add('jets', okj, f'the radial jets w_i, w_ij (Cartesian) and w_x, w_rho, w_xx, w_xrho, w_rhorho equal sympy\'s '
           f'derivatives of v (1 - f(r)) for two profiles at 5 exact points each, {nj} comparisons')

    # ---- the two routes ----------------------------------------------------------------------------------------
    RfK, RmixK = build_K()
    RfC = build_C()
    print(f'  built routes K and C in {time.time() - t0:.1f} s', flush=True)
    try:
        ringK, ringC = ring_frame(RfK), ring_frame(RfC)
        symbolic = True
    except NotPoly as exc:
        ringK = ringC = None
        symbolic = False
        print(f'  entries are not polynomial ({str(exc)[:60]}); only the pointwise checks run', flush=True)

    # Alcubierre's Eulerian energy density, G(n, n) = (1/2)(R_00 + R_11 + R_22 + R_33) = -(1/4) v^2 s^2 a^2
    target = -v ** 2 * s ** 2 * a ** 2 / 4
    res_alc = []
    for nm_, Rf_, rg in (('K', RfK, ringK), ('C', RfC, ringC)):
        Gnn = (Rf_[0, 0] + Rf_[1, 1] + Rf_[2, 2] + Rf_[3, 3]) / 2
        if rg is not None:
            gr = red((rg[0][0] + rg[1][1] + rg[2][2] + rg[3][3]) * QQ(1, 2))
            res_alc.append(req(gr, to_ring(target)))
        for _ in range(20):
            pt = randpt(rng)
            res_alc.append(at_point(Gnn, pt) == at_point(target, pt))
    ck.add('alc_energy', all(res_alc), 'G(n, n) = -(1/4) v^2 s^2 a^2, i.e. rho_E = -v^2 rho^2 f\'^2/(32 pi r^2), '
           'both routes, symbolic and at 20 random exact points each (the Eulerian energy density of Alcubierre\'s '
           'metric, Eq. (19) of the 1994 paper, quoted in he_type.md and review.md and not re-read from the paper)')

    # routes agree entry by entry, (I1), and the coordinate cubic
    res_ag = []
    names = []
    if symbolic:
        res_ag.append(all(req(ringK[i][j], ringC[i][j]) for i in range(4) for j in range(4)))
        names.append('symbolic')
    pts = [randpt(rng) for _ in range(30)]
    for pt in pts:
        res_ag.append(all(at_point(RfK[i, j] - RfC[i, j], pt) == 0 for i in range(4) for j in range(4)))
    I1 = all(at_point(RfK[3, j], pt) == 0 and at_point(RfC[3, j], pt) == 0 for pt in pts for j in range(4))
    # coordinate-basis mixed Ricci of route K: det = 0, and the cubic of the 4x4 equals the cubic of the 3-block
    def coordcubic(pt):
        Mx = Matrix(4, 4, lambda i, j: sp.Rational(at_point(RmixK[i, j], pt).numerator,
                                                   at_point(RmixK[i, j], pt).denominator))
        lam = sp.Symbol('lam')
        P = sp.Poly((lam * sp.eye(4) - Mx).det(), lam)
        return [to_fr(P.coeff_monomial(lam ** k)) for k in range(5)]
    res_cc = []
    for pt in pts[:10]:
        P = coordcubic(pt)
        c2_, c1_, c0_ = cubic(block(frame_at(RfC, pt)))
        res_cc.append(P[0] == 0 and P[1] == c0_ and P[2] == c1_ and P[3] == c2_)
    ck.add('routes_agree', all(res_ag) and all(res_cc) and I1,
           f'R_ab of route K (coordinates) = R_ab of route C (Cartan), {"symbolic and " if symbolic else ""}'
           f'at 30 exact points; R(e_phi, .) = 0 (I1) {I1}; det of the 4x4 mixed Ricci of route K is 0 and its '
           f'cubic equals the 3-block cubic at 10 points {all(res_cc)}')

    # ---- discriminant tools -------------------------------------------------------------------------------------
    A2, B2, C2, X2 = sp.symbols('A2 B2 C2 X2')
    sym_ok = sp.expand(sp.discriminant(X2 ** 3 + A2 * X2 ** 2 + B2 * X2 + C2, X2) - disc_formula(A2, B2, C2)) == 0
    pl = []
    for roots in ((1, 2, 3), (0, 1, -1), (2, 2, 5), (1, sp.I, -sp.I), (3, 1 + 2 * sp.I, 1 - 2 * sp.I)):
        poly = sp.expand((X2 - roots[0]) * (X2 - roots[1]) * (X2 - roots[2]))
        co = sp.Poly(poly, X2).all_coeffs()
        prod = sp.expand((roots[0] - roots[1]) ** 2 * (roots[0] - roots[2]) ** 2 * (roots[1] - roots[2]) ** 2)
        pl.append(sp.simplify(disc_formula(co[1], co[2], co[3]) - prod) == 0)
    sturm = []
    for _ in range(60):
        cc_ = [Fr(rng.randint(-20, 20), rng.randint(1, 6)) for _ in range(3)]
        Dv = disc_formula(*cc_)
        nreal = sp.Poly([1] + [sp.Rational(x_.numerator, x_.denominator) for x_ in cc_], X2).count_roots()
        sturm.append((Dv > 0 and nreal == 3) or (Dv < 0 and nreal == 1) or (Dv == 0 and nreal in (1, 2, 3)))
    ck.add('disc_tools', sym_ok and all(pl) and all(sturm),
           f'disc formula = sympy.discriminant (symbolic) {sym_ok}; = product of squared root differences for 5 '
           f'planted cubics {all(pl)}; sign(D) agrees with Sturm\'s real-root count at 60 random rational cubics '
           f'{all(sturm)}')

    # ---- planted tensors for the type classifier ---------------------------------------------------------------
    def lorentz():
        def boost(i, j, ch, sh):
            L = sp.eye(3)
            L[i, i] = L[j, j] = ch
            L[i, j] = L[j, i] = sh
            return L
        def rot(i, j, co, si):
            L = sp.eye(3)
            L[i, i] = L[j, j] = co
            L[i, j] = -si
            L[j, i] = si
            return L
        return boost(0, 1, sp.Rational(5, 4), sp.Rational(3, 4)) * rot(1, 2, sp.Rational(3, 5), sp.Rational(4, 5)) * \
            boost(0, 2, sp.Rational(13, 12), sp.Rational(5, 12))
    L = lorentz()
    assert sp.simplify(L.T * ETA3 * L - ETA3) == sp.zeros(3, 3)
    ell = sp.Matrix([1, 1, 0])
    m_ = sp.Matrix([0, 0, 1])
    planted = {
        'perfect fluid (mu, p) = (3/2, 1/5)': (sp.diag(sp.Rational(3, 2), sp.Rational(1, 5), sp.Rational(1, 5)), 'I'),
        'dust': (sp.diag(2, 0, 0), 'I'),
        'anisotropic fluid': (sp.diag(sp.Rational(1, 2), -1, 3), 'I'),
        'lambda eta': (5 * ETA3, 'I'),
        'null dust': (ell * ell.T, 'II'),
        'null dust on a fluid, triple root': (3 * ETA3 + 2 * ell * ell.T, 'II'),
        'Jordan 2x2 plus a simple root': (sp.Matrix([[-1 + 1, 1, 0], [1, 1 + 1, 0], [0, 0, 7]]), 'II'),
        'ell m + m ell': (ell * m_.T + m_ * ell.T, 'III'),
        'complex pair, trace 0': (sp.Matrix([[1, 2, 0], [2, -1, 0], [0, 0, 3]]), 'IV'),
        'complex pair, trace != 0': (sp.Matrix([[1, 3, 0], [3, 0, 0], [0, 0, -2]]), 'IV'),
    }
    out = []
    for nm_, (R0, want) in planted.items():
        for Lm in (sp.eye(3), L):
            Rm_ = Lm.T * R0 * Lm
            Mm = (ETA3 * Rm_)
            M = [[to_fr(Mm[i, j]) for j in range(3)] for i in range(3)]
            got = classify_lemma(M)
            dd = disc_formula(*cubic(M))
            if dd == 0:
                got2 = classify_direct_exact(M)
            else:
                ev_ = eig_eta_mp(M)
                got2 = 'I' if any(x_ < 0 for x_ in ev_) else 'IV'
            out.append((nm_, got == want and got2 == want))
    bad = sorted({n for n, ok in out if not ok})
    ck.add('planted_types', not bad, f'{len(planted)} tensors of known type, plain and in a frame boosted twice and '
           f'rotated, get their type from the lemma (D, then rank) and from the definition (timelike or null '
           f'eigenvectors, Jordan sizes){"; wrong: " + ", ".join(bad) if bad else ""}')

    # ---- the four identities ------------------------------------------------------------------------------------
    cf_r = closed_forms(RING_CTX)

    # T1: the powers of v, p6 v^6 = -4 |j|^6, and the weights
    t1 = []
    detail1 = []
    if symbolic:
        for nm_, rg in (('K', ringK), ('C', ringC)):
            D = red(Dof(rg))
            exps = sorted({m[IV] for m in D.keys()})
            wts = all({m[IV] for m in rg[0][j].keys()} <= {1} for j in (1, 2, 3)) and all(
                {m[IV] for m in rg[i][j].keys()} <= {2} for (i, j) in ((0, 0), (1, 1), (1, 2), (2, 2)))
            j2 = red(rg[0][1] ** 2 + rg[0][2] ** 2 + rg[0][3] ** 2)
            t1.append(exps == [6, 8, 10, 12] and wts)
            t1.append(req(4 * j2, v_ ** 2 * cf_r['Jt']))
            t1.append(req(coef_v(D, 6) * v_ ** 6, -4 * j2 ** 3))
            t1.append(req(16 * coef_v(D, 6), cf_r['p6x16']))
            detail1.append(f'route {nm_}: powers of v in D = {exps}, weights ok {wts}, #terms {len(D)}, '
                           f'p6 v^6 = -4|j|^6 and p6 = -J^6/(16 r^6)')
    ptsT1 = []
    for _ in range(12):
        pt0 = randpt(rng)
        byv = {}
        for k_ in range(1, 6):
            pt = dict(pt0)
            pt[v] = Fr(k_, 3)
            byv[pt[v]] = (frame_at(RfK, pt), frame_at(RfC, pt), pt)
        for idx in (0, 1):
            co, res = fit_cubic_in_V({vv: byv[vv][idx] for vv in byv}, pt0)
            cfp = closed_forms(ctx_of(pt0))
            t1.append(res == 0 and 16 * co[0] == cfp['p6x16'])
    ck.add('T1_structure', all(t1), '; '.join(detail1) + ('; ' if detail1 else '') +
           'pointwise: D/v^6 is a cubic in v^2 (5 exact values of v, 12 points, both routes) with p6 as stated')

    # T2: p12 >= 0
    t2 = []
    detail2 = []
    if symbolic:
        for nm_, rg in (('K', ringK), ('C', ringC)):
            D = red(Dof(rg))
            D0 = red(Dof(rg, jzero=True))
            M0 = block(rg, jzero=True)
            S11, S12, S22 = rg[1][1], rg[1][2], rg[2][2]
            lam0 = M0[0][0]
            chi = red((lam0 - S11) * (lam0 - S22) - S12 ** 2)
            Dl = red((S11 - S22) ** 2 + 4 * S12 ** 2)
            t2.append(req(16 * coef_v(D, 12), cf_r['p12x16']))
            t2.append(req(D0, v_ ** 12 * coef_v(D, 12)))
            t2.append(req(D0, chi ** 2 * Dl))
            # constants relating chi and the discriminant of S to Q1 and Q2
            kc = red(chi * 4) - red(v_ ** 4 * s_ ** 2 * F_ * cf_r['Q1t'])
            kd = red(Dl) - red(v_ ** 4 * cf_r['Q2t'])
            t2.append(kc == 0 and kd == 0)
            detail2.append(f'route {nm_}: p12 v^12 = D(j = 0) = chi_S(lam0)^2 [(S11 - S22)^2 + 4 S12^2] = '
                           f'16^-1 s^4 F^2 Q1^2 Q2 r^-6 v^12, chi_S(lam0) = (1/4) v^4 s^2 F u^2 Q1 and '
                           f'(S11 - S22)^2 + 4 S12^2 = v^4 u^2 Q2, both exact')
    for _ in range(40):
        pt = randpt(rng)
        for Rf_ in (RfK, RfC):
            Rfp = frame_at(Rf_, pt)
            M0 = block(Rfp, jzero=True)
            S11, S12, S22 = Rfp[1][1], Rfp[1][2], Rfp[2][2]
            lam0 = M0[0][0]
            chi = (lam0 - S11) * (lam0 - S22) - S12 ** 2
            Dj = Dof(Rfp, jzero=True)
            t2.append(Dj == chi ** 2 * ((S11 - S22) ** 2 + 4 * S12 ** 2) and Dj >= 0)
            cfp = closed_forms(ctx_of(pt))
            t2.append(16 * Dj == pt[v] ** 12 * cfp['p12x16'])
    ck.add('T2_p12', all(t2), '; '.join(detail2) + ('; ' if detail2 else '') +
           'pointwise at 40 points per route: D(j = 0) = chi_S(lam0)^2 Delta_S >= 0 and equals v^12 p12')

    # T3: the axis
    t3 = []
    detail3 = []
    if symbolic:
        for nm_, rg in (('K', ringK), ('C', ringC)):
            D = red(Dof(rg))
            t3.append(req(sub_axis(D), sub_axis(cf_r['Dax'])))
            detail3.append(f'route {nm_}: D(s = 0, c^2 = 1) = -4 a^2 v^6 Q_ax^2 / r^6')
    for sign in (1, -1):
        for _ in range(25):
            pt = randpt(rng, c0=Fr(sign), s0=Fr(0))
            cfp = closed_forms(ctx_of(pt))
            for Rf_ in (RfK, RfC):
                t3.append(Dof(frame_at(Rf_, pt)) == cfp['Dax'])
    # the complex pair on the axis, and D < 0 wherever f' != 0
    pair = []
    for _ in range(20):
        pt = randpt(rng, c0=Fr(1), s0=Fr(0))
        Rfp = frame_at(RfK, pt)
        Dv = Dof(Rfp)
        pair.append(Dv < 0)
        alpha = -Rfp[0][0]
        beta_ = Rfp[0][1]
        c2_, c1_, c0_ = cubic(block(Rfp))
        mu = -c2_ - 2 * alpha
        quad = [Fr(1), -2 * alpha, alpha ** 2 + beta_ ** 2]
        poly3 = [quad[0], quad[1] - mu, quad[2] - mu * quad[1], -mu * quad[2]]
        pair.append(poly3 == [Fr(1), c2_, c1_, c0_])
    ck.add('T3_axis', all(t3) and all(pair), '; '.join(detail3) + ('; ' if detail3 else '') +
           'pointwise at 50 axis points per route; D < 0 at 20 axis points with a != 0 and the cubic factors as '
           '(l - mu)[(l + R_00)^2 + R_01^2] {}'.format(all(pair)))

    # T4: the equator and the criterion
    t4 = []
    detail4 = []
    crit_factor = 2 if mut('crit-factor') else 1
    if symbolic:
        for nm_, rg in (('K', ringK), ('C', ringC)):
            D = red(Dof(rg))
            t4.append(req(16 * sub_equator(D), sub_equator(cf_r['Deq16'])))
            # the block splits (n, e_x) + (e_rho) on the equator
            split = all(sub_equator(rg[i][j]) == 0 for (i, j) in ((0, 2), (1, 2)))
            disc2 = red((rg[0][0] + rg[1][1]) ** 2 - 4 * rg[0][1] ** 2)
            t4.append(split)
            t4.append(req(sub_equator(disc2), v_ ** 2 * sub_equator(cf_r['T'])))
            M = block(rg)
            mu = M[2][2]
            chi2 = red((mu - M[0][0]) * (mu - M[1][1]) - M[0][1] * M[1][0])
            t4.append(req(sub_equator(chi2), sub_equator(v_ ** 2 * cf_r['Qeq'] * QQ(1, 4))))
            detail4.append(f'route {nm_}: D(c = 0) = v^6 [r^2 v^2 a^4 - (a + r b)^2] Q_eq^2 / (16 r^6); the block '
                           f'splits as (n, e_x) + (e_rho) with discriminant (R_00 + R_11)^2 - 4 R_01^2 = '
                           f'v^2 u^2 [..] and chi_2(mu) = v^2 u^2 Q_eq / 4')
    # the criterion, from the exact classification
    crit = []
    for _ in range(300):
        pt = randpt(rng, c0=Fr(0), s0=Fr(1))
        for Rf_ in (RfK, RfC):
            Rfp = frame_at(Rf_, pt)
            M = block(Rfp)
            Dv = disc_formula(*cubic(M))
            cfp = closed_forms(ctx_of(pt))
            t4.append(16 * Dv == cfp['Deq16'])
            lhs = pt[v] * pt[a] ** 2
            rhs = crit_factor * abs(pt[a] * pt[u] + pt[b])
            if cfp['Qeq'] != 0:
                crit.append((Dv > 0) == (lhs > rhs))
    ck.add('T4_equator', all(t4) and all(crit), '; '.join(detail4) + ('; ' if detail4 else '') +
           f'{len(crit)} exact random equatorial points per route: D > 0 iff v r a^2 > |a + r b| ({all(crit)}) '
           f'(Q_eq = 0 and D = 0 points are in Qeq0_and_edges)')

    # ---- the lemma against the definition, at random points of the whole wall ----------------------------------
    td = []
    nI = nIV = 0
    for _ in range(150):
        pt = randpt(rng)
        for Rf_ in (RfK, RfC):
            M = block(frame_at(Rf_, pt))
            Dv = disc_formula(*cubic(M))
            if Dv == 0:
                continue
            ev_ = eig_eta_mp(M)
            if Dv > 0:
                td.append(len(ev_) == 3 and any(x_ < 0 for x_ in ev_) and all(x_ != 0 for x_ in ev_))
                nI += 1
            else:
                td.append(len(ev_) == 1 and ev_[0] > 0)
                nIV += 1
    ck.add('type_direct', all(td) and nI > 0 and nIV > 0,
           f'the sign of D against a direct look at the eigenvectors (mpmath, 60 digits): D > 0 gives three real '
           f'eigenvalues, a timelike eigenvector and no null one ({nI} points); D < 0 gives one real eigenvalue with a '
           f'spacelike eigenvector, so no timelike or null eigenvector ({nIV} points)')

    # ---- D = 0 points: Q_eq = 0 and the edges of the equatorial annulus -----------------------------------------
    e0 = []
    typesq = {}
    for _ in range(80):
        # Q_eq = 0 by construction: F + r a = -F q^2, (a + r b) = +-2 v a F q
        vv = Fr(rng.randint(1, 20), 8)
        FF = Fr(rng.randint(1, 16), 16)
        rr_ = Fr(rng.randint(2, 20), 4)
        q = Fr(rng.randint(1, 12), rng.randint(1, 6))
        if q == 1:
            continue
        aa = -FF * (1 + q * q) / rr_
        bb = (rng.choice([-1, 1]) * 2 * vv * aa * FF * q - aa) / rr_
        pt = {v: vv, F: FF, a: aa, b: bb, u: 1 / rr_, c: Fr(0), s: Fr(1)}
        for Rf_ in (RfK, RfC):
            M = block(frame_at(Rf_, pt))
            Dv = disc_formula(*cubic(M))
            cfp = closed_forms(ctx_of(pt))
            if cfp['Qeq'] != 0 or Dv != 0:
                e0.append(False)
                continue
            got, got2 = classify_lemma(M), classify_direct_exact(M)
            holds = pt[v] * aa ** 2 > crit_factor * abs(aa / rr_ + bb)
            typesq[got] = typesq.get(got, 0) + 1
            e0.append(got == got2 == 'I' and holds)
    edges = {}
    for _ in range(80):
        # boundary of the criterion: v r a^2 = |a + r b|, j != 0
        vv = Fr(rng.randint(1, 20), 8)
        FF = Fr(rng.randint(1, 16), 16)
        rr_ = Fr(rng.randint(2, 20), 4)
        aa = Fr(rng.choice([-1, 1]) * rng.randint(1, 12), 4)
        bb = (rng.choice([-1, 1]) * vv * rr_ * aa ** 2 - aa) / rr_
        pt = {v: vv, F: FF, a: aa, b: bb, u: 1 / rr_, c: Fr(0), s: Fr(1)}
        for Rf_ in (RfK, RfC):
            Rfp = frame_at(Rf_, pt)
            M = block(Rfp)
            Dv = disc_formula(*cubic(M))
            jj = Rfp[0][1]
            if jj == 0:
                continue
            if Dv != 0:
                e0.append(False)
                continue
            got, got2 = classify_lemma(M), classify_direct_exact(M)
            edges[got] = edges.get(got, 0) + 1
            e0.append(got == got2 == 'II')
    # the point of he_type.md section 'Checks': a worked Q_eq = 0 example with exact numbers
    ck.add('Qeq0_and_edges', all(e0) and len(e0) > 100,
           f'D = 0 on the equator: {len(e0)} exact points per route. At Q_eq = 0 (a != 0, r a + 2F != 0) the point is '
           f'Type I and v r a^2 > |a + r b| holds ({typesq}); on the edge v r a^2 = |a + r b| with j != 0 it is '
           f'Type II ({edges}); lemma and definition agree at every point')

    ck.add('timing', True, f'{time.time() - t0:.1f} s')
    if symbolic and not MUT:
        print_results(ringC)
    return ck


def print_results(rg):
    """What the derivation found, for the review: the frame Ricci tensor, the cubic and the coefficients p6 .. p12."""
    r = sp.Symbol('r')
    show = lambda p: sp.factor(p.as_expr().subs(u, 1 / r))
    print('\nderived quantities (frame (n, e_x, e_rho, e_phi), r = 1/u, c = cos, s = sin, F = 1 - f, a = f\', '
          'b = f\'\')')
    for (i, j) in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)):
        print(f'  R_{i}{j} = {show(rg[i][j])}')
    print('  R_3a = 0 for every a')
    M = block(rg)
    c2_, c1_, c0_ = cubic(M)
    print(f'  c2 = {show(red(c2_))}')
    print(f'  c1 = {show(red(c1_))}')
    print(f'  c0 = {show(red(c0_))}')
    D = red(Dof(rg))
    for k_ in (6, 8, 10, 12):
        pk = coef_v(D, k_)
        print(f'  p{k_} ({len(pk)} terms) = {show(pk)}')
    j2 = red(rg[0][1] ** 2 + rg[0][2] ** 2 + rg[0][3] ** 2)
    print(f'  |j|^2 = {show(j2)}')


# ------------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mutate', default=None)
    ap.add_argument('--no-controls', action='store_true')
    args = ap.parse_args()
    if args.mutate:
        names = args.mutate.split(',')
        for nme in names:
            if nme not in CONTROLS:
                sys.exit(f'unknown control {nme}; known: {", ".join(CONTROLS)}')
        MUT.update(names)
        print(f'mutations: {sorted(MUT)}')
    ck = main_checks()
    failed = ck.failed()
    print(('FAILED: ' + ','.join(failed)) if failed else f'all {len(ck.rows)} checks passed', flush=True)
    if args.mutate:
        sys.exit(1 if failed else 0)
    exit_code = 1 if failed else 0
    if not args.no_controls:
        print('\ncontrols (each is the same script with --mutate=NAME and must fail its target checks)')
        rows = []
        for name, (what, targets) in CONTROLS.items():
            t0 = time.time()
            p = subprocess.run([sys.executable, '-u', __file__, f'--mutate={name}'], capture_output=True, text=True)
            line = [ln for ln in p.stdout.splitlines() if ln.startswith('FAILED:') or ln.startswith('all ')]
            got = set(line[-1][len('FAILED: '):].split(',')) if line and line[-1].startswith('FAILED:') else set()
            fired = p.returncode == 1 and set(targets) <= got
            rows.append((name, fired, sorted(got), time.time() - t0, what))
            print(f"[{'FIRED' if fired else 'DID NOT FIRE'}] --mutate={name}: {what}; exit {p.returncode}; "
                  f"checks that fail: {', '.join(sorted(got)) or 'none'}", flush=True)
            if not fired:
                exit_code = 1
        print(f"{sum(r[1] for r in rows)} of {len(rows)} controls fired")
    sys.exit(exit_code)


if __name__ == '__main__':
    main()

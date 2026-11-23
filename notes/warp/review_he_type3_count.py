#!/usr/bin/env python3
"""Independent interval review of #133's count of the Type III points of the tanh wall (issue #136).

Shares nothing with he_type3_count.py or he_type3.py (not opened, not imported). Built from the frame table printed in
he_type3.md section 1 and the statement of he_type3_count.md. Interval arithmetic is exact rational arithmetic
(fractions.Fraction) with outward rounding to a p-bit grid; exp, cosh and sinh are enclosed by Taylor series with a
remainder bound. mpmath is used only at point precision (mp, never iv) as an independent numerical reference.

Usage: review_he_type3_count.py [--checks=a,b,...] [--mutate=NAME] [--bits=N]   (stdout only; exit 1 if any check fails)
"""
import sys
import time
import math
import random
import argparse
import traceback
from fractions import Fraction as Fr

MUT = None        # name of the active mutation (controls), set from --mutate
BITS = 160        # working precision of the interval arithmetic, significant bits


# ---------------------------------------------------------------------------------------------------------------------
# 1. Rational intervals with outward rounding
# ---------------------------------------------------------------------------------------------------------------------

def _scale(x, p):
    """Return s such that |x| 2^s lies in (2^(p-1), 2^(p+1)): the grid of spacing 2^-s has p significant bits."""
    return p - (abs(x.numerator).bit_length() - x.denominator.bit_length())


def rd(x, p=None):
    """Largest number <= x on the p-bit grid."""
    n = x.numerator
    if n == 0:
        return x
    s = (BITS if p is None else p) - (abs(n).bit_length() - x.denominator.bit_length())
    d = x.denominator
    m = (n << s) // d if s >= 0 else n // (d << -s)
    return Fr(m, 1 << s) if s > 0 else Fr(m << -s)


def ru(x, p=None):
    """Smallest number >= x on the p-bit grid."""
    return -rd(-x, p)


def ulp(x, p=None):
    """Spacing of the p-bit grid at x ("one unit in the last place")."""
    s = _scale(Fr(x), BITS if p is None else p)
    return Fr(1, 1 << s) if s > 0 else Fr(1 << -s)


class Iv:
    """Closed interval [lo, hi] with Fraction endpoints. Every operation rounds outward to the BITS-bit grid."""
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        lo = Fr(lo)
        hi = lo if hi is None else Fr(hi)
        if lo > hi:
            raise ValueError('empty interval [%s, %s]' % (float(lo), float(hi)))
        self.lo, self.hi = lo, hi

    @staticmethod
    def _mk(lo, hi):
        if lo > hi:
            raise ValueError('empty interval [%s, %s]' % (float(lo), float(hi)))
        r = Iv.__new__(Iv)
        r.lo, r.hi = lo, hi
        return r

    def __repr__(self):
        return 'Iv[%s, %s]' % (float(self.lo), float(self.hi))

    def mid(self):
        return (self.lo + self.hi) / 2

    def width(self):
        return self.hi - self.lo

    def rad(self):
        return (self.hi - self.lo) / 2

    def contains(self, x):
        return self.lo <= x <= self.hi

    def contains0(self):
        return self.lo <= 0 <= self.hi

    def pos(self):
        return self.lo > 0

    def neg(self):
        return self.hi < 0

    def subset(self, o):
        return o.lo <= self.lo and self.hi <= o.hi

    def interior_in(self, o):
        return o.lo < self.lo and self.hi < o.hi

    def __neg__(self):
        return Iv._mk(-self.hi, -self.lo)

    def __add__(a, b):
        if not isinstance(b, Iv):
            b = Iv(b)
        return Iv._mk(rd(a.lo + b.lo), ru(a.hi + b.hi))

    __radd__ = __add__

    def __sub__(a, b):
        if not isinstance(b, Iv):
            b = Iv(b)
        return Iv._mk(rd(a.lo - b.hi), ru(a.hi - b.lo))

    def __rsub__(a, b):
        return Iv(b) - a

    def __mul__(a, b):
        if not isinstance(b, Iv):
            b = Fr(b)
            if b >= 0:
                return Iv._mk(rd(a.lo * b), ru(a.hi * b))
            return Iv._mk(rd(a.hi * b), ru(a.lo * b))
        al, ah, bl, bh = a.lo, a.hi, b.lo, b.hi
        if al >= 0 and bl >= 0:
            lo, hi = al * bl, ah * bh
        elif ah <= 0 and bh <= 0:
            lo, hi = ah * bh, al * bl
        elif al >= 0 and bh <= 0:
            lo, hi = ah * bl, al * bh
        elif ah <= 0 and bl >= 0:
            lo, hi = al * bh, ah * bl
        else:
            q = (al * bl, al * bh, ah * bl, ah * bh)
            lo, hi = min(q), max(q)
        return Iv._mk(rd(lo), ru(hi))

    __rmul__ = __mul__

    def recip(self):
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError('interval contains 0: %r' % (self,))
        return Iv._mk(rd(1 / self.hi), ru(1 / self.lo))

    def __truediv__(a, b):
        if not isinstance(b, Iv):
            b = Iv(b)
        return a * b.recip()

    def __rtruediv__(a, b):
        return Iv(b) * a.recip()

    def sq(self):
        if self.lo >= 0:
            return Iv._mk(rd(self.lo ** 2), ru(self.hi ** 2))
        if self.hi <= 0:
            return Iv._mk(rd(self.hi ** 2), ru(self.lo ** 2))
        return Iv._mk(Fr(0), ru(max(self.lo ** 2, self.hi ** 2)))

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        if n == 0:
            return Iv(1)
        if n == 1:
            return self
        h = self.__pow__(n // 2)
        r = h.sq()
        return r * self if n % 2 else r

    def hull(a, b):
        return Iv._mk(min(a.lo, b.lo), max(a.hi, b.hi))

    def intersect(a, b):
        lo, hi = max(a.lo, b.lo), min(a.hi, b.hi)
        return None if lo > hi else Iv._mk(lo, hi)


def iv_sum(xs):
    lo = hi = Fr(0)
    for x in xs:
        lo += x.lo
        hi += x.hi
    return Iv._mk(rd(lo), ru(hi))


# ---------------------------------------------------------------------------------------------------------------------
# 2. Rational enclosures of exp, cosh and sinh (Taylor series with a remainder bound)
# ---------------------------------------------------------------------------------------------------------------------

_exp_cache = {}


def exp_pt(x):
    """[lo, hi] with lo <= e^x <= hi, endpoints on the BITS-bit grid, hi - lo at most two grid spacings.

    Method: y = x / 2^k with |y| <= 1/2; e^y in [S - R, S + R] with S the Taylor sum of degree N-1 and
    R = 2 |y|^N / N! (the tail of a series whose term ratio is at most |y|/(N+1) <= 1/4 is at most 4/3 of its first
    term); then k interval squarings, each rounded outward to the working grid of BITS + 40 bits.
    "Last place" below means the spacing of the BITS-bit grid at the endpoint (ulp).
    """
    x = Fr(x)
    key = (x, BITS)
    r = _exp_cache.get(key)
    if r is None:
        r = _exp_pt(x)
        _exp_cache[key] = r
    return r


def _exp_pt(x):
    if x == 0:
        lo, hi = Fr(1), Fr(1)
    else:
        pw = BITS + 40
        k = 0
        while abs(x) > Fr(1, 2) * (1 << k):
            k += 1
        y = x / (1 << k)
        # Taylor degree: 2 |y|^N / N! < 2^-(pw+8)
        N, term = 1, abs(y)
        while True:
            if 2 * term < Fr(1, 1 << (pw + 8)):
                break
            N += 1
            term = term * abs(y) / N
        S, t = Fr(0), Fr(1)
        for n in range(N):
            S += t
            t = t * y / (n + 1)
        R = 2 * term
        lo, hi = rd(S - R, pw), ru(S + R, pw)
        for _ in range(k):
            lo, hi = rd(lo * lo, pw), ru(hi * hi, pw)
        lo, hi = rd(lo), ru(hi)
    if MUT == 'exp-ulp':
        lo, hi = lo + ulp(lo), hi - ulp(hi)       # control: narrowed by one unit in the last place at each end
    return (lo, hi)


def exp_iv(a, b):
    """Enclosure of e^t for t in [a, b] (monotone)."""
    return Iv._mk(exp_pt(a)[0], exp_pt(b)[1])


def cosh_iv(a, b):
    """Enclosure of cosh t for t in [a, b], with 0 <= a (cosh is increasing there)."""
    assert a >= 0
    lo = rd((exp_pt(a)[0] + exp_pt(-a)[0]) / 2)
    hi = ru((exp_pt(b)[1] + exp_pt(-b)[1]) / 2)
    return Iv._mk(lo, hi)


def sinh_iv(a, b):
    """Enclosure of sinh t for t in [a, b] (sinh is increasing)."""
    lo = rd((exp_pt(a)[0] - exp_pt(-a)[1]) / 2)
    hi = ru((exp_pt(b)[1] - exp_pt(-b)[0]) / 2)
    return Iv._mk(lo, hi)


def cosh_const(x):
    """Enclosure of cosh x for a rational x >= 0."""
    return cosh_iv(Fr(x), Fr(x))


# ---------------------------------------------------------------------------------------------------------------------
# 3. Forward-mode derivative (dual numbers over intervals)
# ---------------------------------------------------------------------------------------------------------------------

class Du:
    """(value, derivative) pair of intervals. The derivative is with respect to the one independent variable z."""
    __slots__ = ('v', 'd')

    def __init__(self, v, d=None):
        self.v = v if isinstance(v, Iv) else Iv(v)
        self.d = (Iv(0) if d is None else d) if not isinstance(d, (int, Fr)) else Iv(d)

    def __add__(a, b):
        if not isinstance(b, Du):
            b = Du(b)
        return Du(a.v + b.v, a.d + b.d)

    __radd__ = __add__

    def __neg__(a):
        return Du(-a.v, -a.d)

    def __sub__(a, b):
        if not isinstance(b, Du):
            b = Du(b)
        return Du(a.v - b.v, a.d - b.d)

    def __rsub__(a, b):
        return Du(b) - a

    def __mul__(a, b):
        if not isinstance(b, Du):
            if isinstance(b, Iv):
                return Du(a.v * b, a.d * b)
            b = Fr(b)
            return Du(a.v * b, a.d * b)
        return Du(a.v * b.v, a.v * b.d + a.d * b.v)

    __rmul__ = __mul__

    def recip(a):
        r = a.v.recip()
        return Du(r, -(a.d * r.sq()))

    def __truediv__(a, b):
        if not isinstance(b, Du):
            if isinstance(b, Iv):
                return a * b.recip()
            return a * (1 / Fr(b))
        return a * b.recip()

    def __rtruediv__(a, b):
        return Du(b) * a.recip()

    def sq(a):
        return Du(a.v.sq(), 2 * a.v * a.d)


def zvar(a, b):
    """The independent variable z on the cell [a, b]."""
    return Du(Iv(a, b), Iv(1))


def du_cosh_sinh(z):
    """(cosh z, sinh z) as Du for a Du z on a cell of nonnegative z with a derivative 1 (z the independent variable)."""
    a, b = z.v.lo, z.v.hi
    c = cosh_iv(a, b)
    s = sinh_iv(a, b)
    return Du(c, s * z.d), Du(s, c * z.d)


# ---------------------------------------------------------------------------------------------------------------------
# 4. E1, E2 and P from the frame table of he_type3.md section 1 (typed from the text)
# ---------------------------------------------------------------------------------------------------------------------

class Poly:
    """Sparse polynomial {exponent tuple: Fraction}, evaluated on intervals or dual numbers by grouped Horner sums."""

    def __init__(self, terms):
        self.terms = {e: Fr(c) for e, c in terms.items() if c != 0}
        self.nv = len(next(iter(self.terms)))
        self.maxe = [max(e[v] for e in self.terms) for v in range(self.nv)]
        tree = {}
        for e, c in self.terms.items():
            node = tree
            for x in e[:-1]:
                node = node.setdefault(x, {})
            node.setdefault(e[-1], c)
        self.tree = tree

    def __len__(self):
        return len(self.terms)

    def eval(self, vals, one):
        pw = []
        for v in range(self.nv):
            p = [one]
            if self.maxe[v] > 0:
                p.append(vals[v])
            for _ in range(2, self.maxe[v] + 1):
                p.append(p[-1] * vals[v])
            pw.append(p)

        def rec(node, depth):
            if depth == self.nv - 1:
                acc = None
                for e, c in node.items():
                    t = pw[depth][e] * c
                    acc = t if acc is None else acc + t
                return acc
            acc = None
            for e, sub in node.items():
                t = rec(sub, depth + 1)
                if e > 0:
                    t = pw[depth][e] * t
                acc = t if acc is None else acc + t
            return acc

        return rec(self.tree, 0)


def frame_entries(w, phi, X, g, y, cs):
    """The frame entries of he_type3.md section 1 in the frame (n, e_x, e_rho), cs = c s."""
    R00 = -X * (1 + w) / 2 + phi * (g + w * y)
    R0x = g + (1 - w) * y / 2
    R0r = -cs * y / 2
    Rxx = X * (3 * w - 1) / 2 - phi * (g + w * y)
    Rxr = cs * (X - phi * y / 2)
    Rrr = X * (1 - w) / 2
    if MUT == 'frame-typo':
        Rxr = cs * (X + phi * y / 2)       # control: one sign of the printed table typed wrongly
    return R00, R0x, R0r, Rxx, Rxr, Rrr


_alg = None


def algebra():
    """sympy construction. Returns a dict of symbolic objects and of the term tables used by the interval code."""
    global _alg
    key = MUT if MUT in ('frame-typo', 'drop-P-term') else None
    if _alg is not None and _alg['key'] == key:
        return _alg
    import sympy as sp
    phi, A, B, w, x, cs = sp.symbols('phi A B w x cs')
    X, g, y = A ** 2, A, B - A          # (X, gamma, y) at u = 1 after the weight scaling (a, b, u) -> (l a, l^2 b, l u)
    R00, R0x, R0r, Rxx, Rxr, Rrr = frame_entries(w, phi, X, g, y, cs)
    R = sp.Matrix([[R00, R0x, R0r], [R0x, Rxx, Rxr], [R0r, Rxr, Rrr]])
    M = sp.diag(-1, 1, 1) * R                       # M = eta R, eta = diag(-1, 1, 1)
    cp = sp.Poly(sp.expand((x * sp.eye(3) - M).det()), x)
    c2, c1, c0 = [sp.expand(cp.coeff_monomial(x ** k)) for k in (2, 1, 0)]

    def red(e):                                      # c^2 s^2 = w (1 - w); only even powers of cs occur
        p = sp.Poly(sp.expand(e), cs)
        out = 0
        for (k,), co in p.terms():
            assert k % 2 == 0
            out += co * (w * (1 - w)) ** (k // 2)
        return sp.expand(out)

    c2, c1, c0 = red(c2), red(c1), red(c0)
    E1 = sp.expand(c2 ** 2 - 3 * c1)
    E2 = sp.expand(27 * c0 - c2 ** 3)
    Res = sp.expand(sp.resultant(E1, E2, w))
    Pq = sp.expand(sp.cancel(Res / (sp.Rational(27, 256) * phi ** 3 * (A - B) ** 3)))
    Ptab = {m: Fr(int(c.p), int(c.q)) for m, c in sp.Poly(Pq, phi, A, B).terms()}
    if MUT == 'drop-P-term':
        del Ptab[(3, 0, 9)]                           # control: the term 4 phi^3 B'^9 of P is dropped
    E1tab = {(m[3], m[0], m[1], m[2]): Fr(int(c.p), int(c.q)) for m, c in sp.Poly(E1, phi, A, B, w).terms()}
    E2tab = {(m[3], m[0], m[1], m[2]): Fr(int(c.p), int(c.q)) for m, c in sp.Poly(E2, phi, A, B, w).terms()}
    _alg = dict(key=key, sp=sp, phi=phi, A=A, B=B, w=w, x=x, cs=cs, c2=c2, c1=c1, c0=c0, E1=E1, E2=E2, Res=Res, Pq=Pq,
                P=Poly(Ptab), E1p=Poly(E1tab), E2p=Poly(E2tab), Ptab=Ptab, E1tab=E1tab, E2tab=E2tab)
    return _alg


# ---------------------------------------------------------------------------------------------------------------------
# 5. The tanh wall: f = (1 + C)/(cosh z + C), z = 2 sigma r, C = cosh 2 sigma R_b; (phi, A', B') = (v F, v r f', v r^2 f'')
# ---------------------------------------------------------------------------------------------------------------------

class Wall:
    def __init__(self, name, sigma, Rb, v, zt, radii, ws, Ts, other_w, y_r):
        self.name, self.v, self.zt = name, Fr(v), zt
        self.sigma = Fr(sigma) + (Fr(1, 1000) if MUT == 'sigma' else 0)      # control: sigma + 1e-3
        self.sigma0 = Fr(sigma)
        self.Rb = Fr(Rb)
        self.radii, self.ws, self.Ts, self.other_w, self.y_r = radii, ws, Ts, other_w, y_r
        self.C = cosh_const(2 * self.sigma * self.Rb)
        self.C1 = self.C + 1

    def r_of_z(self, z):
        return z / (2 * self.sigma)


def walls():
    # statement of he_type3_count.md: radii, w, T of the three admissible points; w of the other roots; radius of y's root
    le = Wall('le', 8, 1, Fr(1, 2), 40,
              [Fr('0.990567889806'), Fr('1.02393326111'), Fr('1.09466230932')],
              [Fr('0.0600518'), Fr('0.0409869'), Fr('0.912867')], [Fr('-8.80794'), Fr('7.90635'), Fr('22.0305')],
              [Fr('58.459')], Fr('0.99211'))
    ro = Wall('rodal', 4, 5, 1, 60,
              [Fr('4.71448517823'), Fr('5.07106634302'), Fr('5.10338056385')],
              [Fr('0.957639'), Fr('0.284767'), Fr('0.977030')], [Fr('-4.51919'), Fr('3.7259'), Fr('4.97171')],
              [Fr('4.1007'), Fr('8.7173'), Fr('1.0300')], Fr('4.99374'))
    return {'le': le, 'rodal': ro}


def phiAB(zd, wall):
    """(phi, A', B') as dual numbers in z, from the cosh form of the profile."""
    c, s = du_cosh_sinh(zd)
    C = wall.C
    D = c + C
    iD = D.recip()
    v = wall.v
    phi = (c - 1) * iD * v
    A = zd * s * iD.sq() * wall.C1 * (-v)
    B = zd.sq() * (c.sq() * (-1) + c * C + 2) * iD.sq() * iD * wall.C1 * (-v)
    return phi, A, B


def P_du(zd, wall):
    ph, A, B = phiAB(zd, wall)
    return algebra()['P'].eval((ph, A, B), Du(1))


def y_du(zd, wall):
    ph, A, B = phiAB(zd, wall)
    return B - A


def T_du(zd, wall):
    """B' + A', proportional to T = b + a u = (A' + B') / (v r^2) at r = z / (2 sigma)."""
    ph, A, B = phiAB(zd, wall)
    return B + A


# ---------------------------------------------------------------------------------------------------------------------
# 6. The cover: leaves on [lo, hi] that exclude a root or hold exactly one
# ---------------------------------------------------------------------------------------------------------------------

def sgn(iv):
    return 1 if iv.lo > 0 else (-1 if iv.hi < 0 else 0)


class Cover:
    """Adaptive bisection of [lo, hi] for a function fun(Du) -> Du, with value and derivative enclosures.

    A leaf is one of
      'excl'  : the natural enclosure or the mean-value form excludes 0 on the leaf (no root);
      'mono'  : the derivative enclosure excludes 0 and the end values have the same sign (no root);
      'root'  : the derivative enclosure excludes 0 and the end values have opposite signs (exactly one root);
      'open'  : a leaf narrower than minw that is none of these (unresolved: fails the check).
    """

    def __init__(self, fun, lo, hi, minw=Fr(1, 1 << 24)):
        self.fun, self.lo, self.hi, self.minw = fun, Fr(lo), Fr(hi), minw
        self.pts = {}
        self.leaves = []
        self.nevals = 0

    def pt(self, z):
        r = self.pts.get(z)
        if r is None:
            r = self.fun(zvar(z, z)).v
            self.pts[z] = r
        return r

    def run(self):
        stack = [(self.lo, self.hi)]
        while stack:
            a, b = stack.pop()
            r = self.classify(a, b)
            if r is not None:
                self.leaves.append((a, b, r[0], r[1]))
                continue
            m = self.split(a, b)
            stack.append((m, b))
            stack.append((a, m))
        self.leaves.sort()
        return self.leaves

    def split(self, a, b):
        m = (a + b) / 2
        if self.pt(m).contains0():                    # a point too close to a root to give a sign: move the split
            m = a + (b - a) * Fr(3, 8)
        return m

    def classify(self, a, b):
        """(kind, sign info) for a leaf, None to bisect."""
        F = self.fun(zvar(a, b))
        self.nevals += 1
        if not F.v.contains0():
            return 'excl', sgn(F.v)
        m = (a + b) / 2
        mv = self.pt(m) + F.d * (Iv(a, b) - m)
        if not mv.contains0():
            return 'excl', sgn(mv)
        if not F.d.contains0():
            sa, sb = sgn(self.pt(a)), sgn(self.pt(b))
            if sa and sb:
                return ('root', (sa, sb)) if sa != sb else ('mono', sa)
        if b - a < self.minw:
            return 'open', None
        return None


def tiling_ok(leaves, lo, hi):
    if not leaves or leaves[0][0] != lo or leaves[-1][1] != hi:
        return False
    return all(leaves[i][1] == leaves[i + 1][0] for i in range(len(leaves) - 1))


# ---------------------------------------------------------------------------------------------------------------------
# 7. Root enclosures (interval Newton), Krawczyk boxes in (z, w), exclusion of every other common root
# ---------------------------------------------------------------------------------------------------------------------

class bits_ctx:
    """Temporarily set the working precision."""

    def __init__(self, bits):
        self.bits = bits

    def __enter__(self):
        global BITS
        self.old, BITS = BITS, self.bits

    def __exit__(self, *a):
        global BITS
        BITS = self.old


def dyadic(x, p=64):
    """A short dyadic rational near x (a point, not an enclosure)."""
    x = Fr(x)
    return rd(x, p) if x != 0 else x


def refine_root(fun, a, b, target=Fr(1, 10 ** 30), maxit=40):
    """Enclosure [lo, hi] of the one root of fun in [a, b] (derivative of fixed sign there) by interval Newton."""
    Z = Iv(a, b)
    for _ in range(maxit):
        if Z.width() < target:
            break
        m = dyadic(Z.mid(), 3 * BITS // 2)
        Fm = fun(zvar(m, m)).v
        Fd = fun(zvar(Z.lo, Z.hi)).d
        if Fd.contains0():
            raise ValueError('refine_root: derivative contains 0')
        N = Iv(m) - Fm / Fd
        Zn = N.intersect(Z)
        if Zn is None:
            raise ValueError('refine_root: empty Newton intersection')
        if Zn.lo == Z.lo and Zn.hi == Z.hi:
            break
        Z = Zn
    return Z


def mat_inv2(M):
    (a, b), (c, d) = M
    det = a * d - b * c
    return [[d / det, -b / det], [-c / det, a / det]]


def krawczyk(model, zc, wc, hz, hw):
    """Krawczyk test for F = (E1, E2) on X = [zc +- hz] x [wc +- hw].

    Returns (ok, K) with ok True iff K(X) lies in the interior of X, that is, F has exactly one zero in X.
    K(X) = x0 - Y F(x0) + (I - Y J(X)) (X - x0), Y the inverse of the midpoint Jacobian at x0.
    """
    zc, wc, hz, hw = Fr(zc), Fr(wc), Fr(hz), Fr(hw)
    zp, wp = Iv(zc), Iv(wc)
    F0 = model.E(Du(zp), Du(wp))
    Fz = model.E(Du(zp, Iv(1)), Du(wp))
    Fw = model.E(Du(zp), Du(wp, Iv(1)))
    J0 = [[Fz[0].d.mid(), Fw[0].d.mid()], [Fz[1].d.mid(), Fw[1].d.mid()]]
    Y = [[dyadic(e, 100) for e in row] for row in mat_inv2(J0)]
    Zb, Wb = Iv(zc - hz, zc + hz), Iv(wc - hw, wc + hw)
    Gz = model.E(Du(Zb, Iv(1)), Du(Wb))
    Gw = model.E(Du(Zb), Du(Wb, Iv(1)))
    JX = [[Gz[0].d, Gw[0].d], [Gz[1].d, Gw[1].d]]
    dX = [Iv(-hz, hz), Iv(-hw, hw)]
    K = []
    for i, x0 in enumerate((zc, wc)):
        t = Iv(x0) - (F0[0].v * Y[i][0] + F0[1].v * Y[i][1])
        for j in range(2):
            Mij = Iv(1 if i == j else 0) - (JX[0][j] * Y[i][0] + JX[1][j] * Y[i][1])
            t = t + Mij * dX[j]
        K.append(t)
    X = (Zb, Wb)
    return all(K[i].interior_in(X[i]) for i in range(2)), K


def common_roots(model, Zroot, wlo, whi, cands, hz=None, hw=Fr(1, 10 ** 20), minw=Fr(1, 10 ** 28), nbox=100000):
    """All common zeros (z, w) of F with z in Zroot and w in [wlo, whi].

    cands: numerical candidates (zc, wc), only starting points. Each is tried as a Krawczyk box. The rest of
    Zroot x [wlo, whi] is covered by boxes that exclude a zero of E1 or E2, or lie inside a verified box.
    Returns (boxes, complete): boxes = verified [(Zbox, Wbox)], complete = every point of the range is accounted for.
    """
    if hz is None:
        hz = max(Zroot.width(), Fr(1, 10 ** 26))
    zc = dyadic(Zroot.mid(), 3 * BITS // 2)
    boxes = []
    for (_, wc) in cands:
        wc = dyadic(wc, 3 * BITS // 2)
        if MUT == 'krawczyk-shift':
            wc += 20 * hw                                  # control: the box is centred twenty half-widths from the zero
        ok, K = krawczyk(model, zc, wc, hz, hw)
        if ok:
            Zb, Wb = Iv(zc - hz, zc + hz), Iv(wc - hw, wc + hw)
            if not any(Wb.lo <= o[1].hi and o[1].lo <= Wb.hi for o in boxes):
                boxes.append((Zb, Wb))
    Zr = Zroot
    stack, n, complete = [Iv(wlo, whi)], 0, True
    while stack:
        W = stack.pop()
        n += 1
        E = model.E(Du(Zr), Du(W))
        if not E[0].v.contains0() or not E[1].v.contains0():
            continue
        if any(Zr.subset(b[0]) and W.subset(b[1]) for b in boxes):
            continue
        if W.width() < minw or n > nbox:
            complete = False
            continue
        m = W.mid()
        stack.append(Iv(W.lo, m))
        stack.append(Iv(m, W.hi))
    return boxes, complete


# ---------------------------------------------------------------------------------------------------------------------
# 8. mpmath at point precision (never iv): independent numerical reference and the candidate finder
# ---------------------------------------------------------------------------------------------------------------------

def mpf_of(x):
    import mpmath as mp
    x = Fr(x)
    return mp.mpf(x.numerator) / x.denominator


def mp_cosh_form(z, wall):
    """(phi, A', B') in mp from the cosh form."""
    import mpmath as mp
    v = mpf_of(wall.v)
    C = mp.cosh(mpf_of(2 * wall.sigma * wall.Rb))
    c, s = mp.cosh(z), mp.sinh(z)
    D = c + C
    return v * (c - 1) / D, -v * z * (1 + C) * s / D ** 2, -v * z ** 2 * (1 + C) * (-c * c + C * c + 2) / D ** 3


def mp_poly(tab, vals):
    import mpmath as mp
    return sum(mpf_of(c) * mp.fprod(vals[n] ** e for n, e in enumerate(ex)) for ex, c in tab.items())


def mp_candidates(wall, zroot):
    """Numerical starting points (z, w) of common roots near zroot (an mp number): roots of E1(w), polished in 2D."""
    import mpmath as mp
    al = algebra()
    ph, A, B = mp_cosh_form(zroot, wall)
    e1 = [mp.mpf(0)] * 3
    for (a, i, j, k), c in al['E1tab'].items():
        e1[a] += mpf_of(c) * ph ** i * A ** j * B ** k
    out = []
    disc = e1[1] ** 2 - 4 * e1[2] * e1[0]
    if e1[2] != 0 and disc >= 0:
        for sgn_ in (1, -1):
            w0 = (-e1[1] + sgn_ * mp.sqrt(disc)) / (2 * e1[2])
            out.append((zroot, w0))
    res = []
    for z0, w0 in out:
        def F(z, w):
            p = mp_cosh_form(z, wall)
            return (mp_poly(al['E1tab'], (w,) + tuple(p)), mp_poly(al['E2tab'], (w,) + tuple(p)))
        try:
            sol = mp.findroot(F, (z0, w0), tol=mp.mpf(10) ** (-mp.mp.dps + 15), maxsteps=30)
            res.append((sol[0], sol[1]))
        except Exception:
            pass
    return res


# ---------------------------------------------------------------------------------------------------------------------
# 9. End bounds: P and y do not vanish on (0, zs] (series with remainder) and on [zt, oo) (dominance)
# ---------------------------------------------------------------------------------------------------------------------

def _fact(n):
    return math.factorial(n)


# positive power series in x = z^2:  K1 = (cosh z - 1)/z^2,  S = sinh z / z,  K2 = (S - 2 K1)/z^2,  K3 = (S - 1)/z^2
SERIES = {
    'K1': lambda n: Fr(1, _fact(2 * n + 2)),
    'S': lambda n: Fr(1, _fact(2 * n + 1)),
    'K2': lambda n: Fr(2 * (n + 1), _fact(2 * n + 4)),
    'K3': lambda n: Fr(1, _fact(2 * n + 3)),
}


def series_iv(name, x_lo, x_hi):
    """Enclosure of sum_n c_n x^n for x in [x_lo, x_hi], 0 <= x_lo. All c_n > 0 and c_{n+1}/c_n decreases in n.

    The sum of the first N terms rounded down is a lower bound; with c_{N+1} x_hi / c_N <= 1/2 the omitted tail is at
    most 2 c_N x_hi^N, which is added to the upper bound. N is the smallest with that tail below 2^-(BITS+20) c_0.
    """
    c = SERIES[name]
    N = 1
    while True:
        if c(N + 1) * x_hi <= c(N) / 2 and 2 * c(N) * x_hi ** N < c(0) / (1 << (BITS + 20)):
            break
        N += 1
    lo = sum((c(n) * x_lo ** n for n in range(N)), Fr(0))
    hi = sum((c(n) * x_hi ** n for n in range(N)), Fr(0)) + 2 * c(N) * x_hi ** N
    return Iv._mk(rd(lo), ru(hi))


def small_z_forms(a, b, wall):
    """Enclosures of Phi = phi/z^2, Psi = (A' + 2 phi)/z^4 and Delta = (B' - A')/z^4 on z in [a, b], 0 <= a < b."""
    x_lo, x_hi = a * a, b * b
    K1, S, K2, K3 = [series_iv(n, x_lo, x_hi) for n in ('K1', 'S', 'K2', 'K3')]
    Zsq = Iv(x_lo, x_hi)
    C, C1, v = wall.C, wall.C1, wall.v
    D = C1 + Zsq * K1                                           # cosh z + C
    Phi = K1 * v / D
    Psi = (K1.sq() * 2 - C1 * K2) * v / D.sq()
    brace = -(C1 * K3) + (C - 2) * K1 - S * K1 - Zsq * K1.sq()
    Delta = -(brace * C1 * v) / (D.sq() * D)
    return Zsq, Phi, Psi, Delta


_gtab = None


def G_table():
    """P(phi, psi - 2 phi, delta + psi - 2 phi) as a table in (phi, psi, delta) and the weight-24 factored form of G."""
    global _gtab
    al = algebra()
    if _gtab is not None and _gtab[0] is al:
        return _gtab[1]
    sp = al['sp']
    ph, ps, dl = sp.symbols('ph ps dl')
    expr = 0
    for (i, j, k), c in al['Ptab'].items():
        Aq = ps - 2 * ph
        Bq = dl + Aq
        expr += sp.Rational(c.numerator, c.denominator) * ph ** i * Aq ** j * Bq ** k
    tab = {m: Fr(int(c.p), int(c.q)) for m, c in sp.Poly(sp.expand(expr), ph, ps, dl).terms()}
    wts = {m: 2 * m[0] + 4 * m[1] + 4 * m[2] for m in tab}
    Gt = {(wts[m] // 2 - 12, m[0], m[1], m[2]): c for m, c in tab.items() if wts[m] >= 24}
    res = dict(tab=tab, wts=wts, minw=min(wts.values()), low={m: c for m, c in tab.items() if wts[m] == min(wts.values())},
               G=Poly(Gt))
    _gtab = (al, res)
    return res


def small_z_cover(wall, zs, maxdepth=14):
    """Cells of [0, zs] on which G > 0 and Delta < 0, with the cell count and the smallest G / (leading term)."""
    gt = G_table()
    G = gt['G']
    cells, bad = [], []
    stack = [(Fr(0), Fr(zs), 0)]
    while stack:
        a, b, d = stack.pop()
        Zsq, Phi, Psi, Delta = small_z_forms(a, b, wall)
        Gv = G.eval((Zsq, Phi, Psi, Delta), Iv(1))
        if Gv.lo > 0 and Delta.hi < 0:
            cells.append((a, b, Gv, Delta))
        elif d >= maxdepth:
            bad.append((a, b, Gv, Delta))
        else:
            m = (a + b) / 2
            stack.append((a, m, d + 1))
            stack.append((m, b, d + 1))
    return cells, bad


_stab = None


def tail_table(wall):
    """S = P / (tau^9 z^9) as a table in (s, T, phihat, alphahat, betahat); s = 1/z, T = tau, with nonnegative exponents."""
    al = algebra()
    terms = {}
    for (i, j, k), c in al['Ptab'].items():
        if k > 9 or j + k < 9:
            raise ValueError('term (%d, %d, %d) of P breaks the tail form' % (i, j, k))
        terms[(9 - k, j + k - 9, i, j, k)] = c * wall.v ** i * (-1) ** j
    return Poly(terms)


def tail_bounds(wall, zt, ncell=64):
    """Enclosures that hold for every z >= zt (e = e^-z in (0, e_t], tau in (0, tau_t], s = 1/z in (0, 1/zt])."""
    et = exp_pt(-Fr(zt))[1]
    E = Iv(0, et)
    C, C1, v = wall.C, wall.C1, wall.v
    den = Iv(1) + C * E * 2 + E.sq()
    phih = (Iv(1) - E).sq() / den
    alph = (Iv(1) - E.sq()) / den.sq()
    beth = ((Iv(1) + E.sq()).sq() - C * E * (Iv(1) + E.sq()) * 2 - E.sq() * 8) / (den.sq() * den)
    tau_t = ru(2 * v * C1.hi * Fr(zt) * et)
    Ti = Iv(0, tau_t)
    tab = tail_table(wall)
    cells = []
    for n in range(ncell):
        s_lo, s_hi = Fr(n, ncell * zt), Fr(n + 1, ncell * zt)
        Sv = tab.eval((Iv(s_lo, s_hi), Ti, phih, alph, beth), Iv(1))
        cells.append((s_lo, s_hi, Sv))
    return dict(e_t=et, tau_t=tau_t, phih=phih, alph=alph, beth=beth, cells=cells)


# ---------------------------------------------------------------------------------------------------------------------
# 10. Models (what the cover and the common-root code are run on) and the planted problem
# ---------------------------------------------------------------------------------------------------------------------

class WallModel:
    """E1, E2 of the wall as functions of dual numbers (z, w)."""

    def __init__(self, wall):
        self.wall = wall

    def E(self, zd, wd):
        ph, A, B = phiAB(zd, self.wall)
        al = algebra()
        one = Du(1)
        return al['E1p'].eval((wd, ph, A, B), one), al['E2p'].eval((wd, ph, A, B), one)

    def T(self, zd):
        return T_du(zd, self.wall)


class PlantedModel:
    """E1 = (w - u1)(w - u2), E2 = (w - v1)(w - v2)(w - v3): common roots exactly where some u_i = v_j."""
    EPS = Fr(1, 10 ** 5)

    @staticmethod
    def uv(z):
        u1 = z * Fr(1, 100) + Fr(3, 10) - Fr(3, 100)
        u2 = z * Fr(1, 50) + Fr(7, 2) - Fr(7, 50)
        v1 = u1 - (z - 3)
        v2 = u2 - (z - 7) * Fr(1, 2)
        v3 = u1 - (z - 5) * (z - 5 - PlantedModel.EPS)
        return u1, u2, v1, v2, v3

    def E(self, zd, wd):
        u1, u2, v1, v2, v3 = self.uv(zd)
        return (wd - u1) * (wd - u2), (wd - v1) * (wd - v2) * (wd - v3)

    def P(self, zd):
        u1, u2, v1, v2, v3 = self.uv(zd)
        r = Du(1)
        for u in (u1, u2):
            for v in (v1, v2, v3):
                r = r * (u - v)
        return r

    def T(self, zd):
        return zd - 5 - self.EPS


PLANTED_ZETAS = ['0.9', '2.3', '2.30001', '7.5', '12', '25', '39.9', '0.2', '41']      # the last two lie outside [1/4, 40]


def planted_c():
    import mpmath as mp
    mp.mp.dps = 60
    cs = []
    for z in PLANTED_ZETAS:
        cs.append(Fr(int(mp.nint(mp.cosh(mp.mpf(z)) * 10 ** 30)), 10 ** 30))
    return cs


def planted_g(zd, cs, double=False):
    c, _ = du_cosh_sinh(zd)
    r = Du(1)
    for ci in cs:
        r = r * (c - ci)
    if double:
        r = r * (c - Fr(int(math.cosh(5.5) * 10 ** 12), 10 ** 12)).sq()      # control: a double root at about 5.5
    return r


# ---------------------------------------------------------------------------------------------------------------------
# 11. Checks
# ---------------------------------------------------------------------------------------------------------------------

CHECKS = {}
STATE = {}
LINES = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    LINES.append(s)


def check(name):
    def deco(f):
        CHECKS[name] = f
        return f
    return deco


def fl(x, n=12):
    import mpmath as mp
    return mp.nstr(mpf_of(x), n)


def _mp_enc_contains(iv, val, slack_rel=0):
    import mpmath as mp
    lo, hi = mpf_of(iv.lo), mpf_of(iv.hi)
    tol = abs(val) * mp.mpf(10) ** (-slack_rel) if slack_rel else 0
    return lo - tol <= val <= hi + tol


@check('exp')
def check_exp():
    import mpmath as mp
    mp.mp.dps = 120
    rng = random.Random(136)
    ok = True
    bad = 0
    maxw = 0
    for _ in range(300):
        x = Fr(rng.randint(-70000, 70000), 1000)
        lo, hi = exp_pt(x)
        ref = mp.exp(mpf_of(x))
        if not (mpf_of(lo) <= ref <= mpf_of(hi)):
            bad += 1
        elif lo > 0:
            maxw = max(maxw, float((hi - lo) / ulp(lo)))
    log('  exp(x), 300 random x in [-70, 70]: %d outside the 120-digit mpmath value; widest enclosure %.1f ulp' % (bad, maxw))
    ok &= bad == 0 and maxw <= 2.0
    bad2 = 0
    for _ in range(100):
        x = Fr(rng.randint(-60000, 60000), 1000)
        y = Fr(rng.randint(-60000, 60000), 1000)
        a, b = exp_pt(x), exp_pt(-x)
        prod = Iv._mk(a[0], a[1]) * Iv._mk(b[0], b[1]) if a[0] <= a[1] and b[0] <= b[1] else None
        if prod is None or not prod.contains(1):
            bad2 += 1
        c1, c2, c3 = exp_pt(x), exp_pt(y), exp_pt(x + y)
        if c1[0] <= c1[1] and c2[0] <= c2[1]:
            pr = Iv._mk(c1[0], c1[1]) * Iv._mk(c2[0], c2[1])
            if pr.hi < c3[0] or pr.lo > c3[1]:
                bad2 += 1
        else:
            bad2 += 1
    log('  exp(x) exp(-x) contains 1 and exp(x) exp(y) meets exp(x+y), 100 pairs: %d failures' % bad2)
    ok &= bad2 == 0
    bad3 = 0
    for _ in range(100):
        a = Fr(rng.randint(0, 60000), 1000)
        b = a + Fr(rng.randint(0, 3000), 1000)
        try:
            ch, sh = cosh_iv(a, b), sinh_iv(a, b)
            d = ch.sq() - sh.sq()
        except ValueError:
            bad3 += 1
            continue
        t = a + (b - a) * Fr(rng.randint(0, 100), 100)
        if not (d.contains(1) and _mp_enc_contains(ch, mp.cosh(mpf_of(t))) and _mp_enc_contains(sh, mp.sinh(mpf_of(t)))):
            bad3 += 1
    log('  cosh^2 - sinh^2 contains 1 and cosh, sinh enclose the mpmath values on 100 random cells: %d failures' % bad3)
    ok &= bad3 == 0
    return ok


@check('profile')
def check_profile():
    import mpmath as mp
    import sympy as sp
    ok = True
    q, k = sp.symbols('q k', positive=True)
    ftanh = ((q * k - 1) / (q * k + 1) - (q - k) / (q + k)) * (k + 1) / (2 * (k - 1))   # e^(2 sigma r) = q, e^(2 sigma R_b) = k
    Cs = (k + 1 / k) / 2
    fcosh = (1 + Cs) / ((q + 1 / q) / 2 + Cs)
    e1 = sp.simplify(sp.cancel(sp.together(ftanh - fcosh)))
    log('  tanh form = (1 + C)/(cosh z + C), exact rational identity in (q, k): difference simplifies to', e1)
    ok &= e1 == 0
    z, C = sp.symbols('z C', positive=True)
    f = (1 + C) / (sp.cosh(z) + C)
    c, s = sp.cosh(z), sp.sinh(z)
    D = c + C
    d1 = sp.simplify((sp.diff(f, z) + (1 + C) * s / D ** 2).rewrite(sp.exp))
    d2 = sp.simplify((sp.diff(f, z, 2) + (1 + C) * (-c ** 2 + C * c + 2) / D ** 3).rewrite(sp.exp))
    log('  f_z = -(1+C) sinh z / D^2 and f_zz = -(1+C)(-cosh^2 z + C cosh z + 2)/D^3 (D = cosh z + C): differences', d1, d2)
    ok &= d1 == 0 and d2 == 0
    mp.mp.dps = 60
    r_, sg, Rb, v_ = sp.symbols('r sigma Rb v')
    ft = (sp.tanh(sg * (r_ + Rb)) - sp.tanh(sg * (r_ - Rb))) / (2 * sp.tanh(sg * Rb))
    fa, fb = sp.diff(ft, r_), sp.diff(ft, r_, 2)
    lam = sp.lambdify((r_, sg, Rb), (ft, fa, fb), 'mpmath')
    rng = random.Random(7)
    bad = 0
    for wn, W in walls().items():
        for _ in range(12):
            r = Fr(rng.randint(5, 4000), 1000) * (1 if wn == 'le' else 3)
            z0 = 2 * W.sigma * r
            if z0 > 90:
                continue
            fv, fa_, fb_ = lam(mpf_of(r), mpf_of(W.sigma), mpf_of(W.Rb))
            v = mpf_of(W.v)
            tgt = (v * (1 - fv), v * mpf_of(r) * fa_, v * mpf_of(r) ** 2 * fb_)
            ph, A, B = phiAB(zvar(z0, z0), W)
            for t, e in zip(tgt, (ph.v, A.v, B.v)):
                if not _mp_enc_contains(e, t, 40):
                    bad += 1
            Tt = fb_ + fa_ / mpf_of(r)
            Te = (A.v + B.v) / (Iv(W.v * r * r))
            if not _mp_enc_contains(Te, Tt, 40):
                bad += 1
    log('  (phi, A\', B\') and T from the tanh definition (sympy derivatives, mpmath 60 digits) against the enclosures at 2 x 12 radii: %d misses' % bad)
    ok &= bad == 0
    return ok


@check('frame')
def check_frame():
    al = algebra()
    sp = al['sp']
    phi, A, B, w = al['phi'], al['A'], al['B'], al['w']
    ok = True
    X, g, y = A ** 2, A, B - A
    c2p = sp.expand((4 * phi * (g + w * y) - X * (3 * w + 1)) / 2)
    E1p = X ** 2 - phi * ((3 * w - 1) * g + 2 * w * y) * X + (phi ** 2 - 3) * g ** 2 + (2 * phi ** 2 * w + 3 * w - 3) * g * y \
        + sp.Rational(1, 4) * (phi ** 2 * w ** 2 + 3 * phi ** 2 * w + 3 * w - 3) * y ** 2
    d_c2, d_E1 = sp.expand(al['c2'] - c2p), sp.expand(al['E1'] - sp.expand(E1p))
    log('  c2 and E1 from the frame table minus the printed c2 and E1 of he_type3.md section 1: %s, %s' % (d_c2, d_E1))
    ok &= d_c2 == 0 and d_E1 == 0
    n1, n2 = len(al['E1tab']), len(al['E2tab'])
    log('  E1, E2 have %d and %d terms (printed: 17 and 42)' % (n1, n2))
    ok &= (n1, n2) == (17, 42)
    # E2 as a cubic form in (X, gamma, y) with 10 terms
    Xs, gs, ys, cs = sp.symbols('Xs gs ys cs')
    R00, R0x, R0r, Rxx, Rxr, Rrr = frame_entries(w, phi, Xs, gs, ys, cs)
    Rm = sp.Matrix([[R00, R0x, R0r], [R0x, Rxx, Rxr], [R0r, Rxr, Rrr]])
    Mm = sp.diag(-1, 1, 1) * Rm
    xx = al['x']
    cp = sp.Poly(sp.expand((xx * sp.eye(3) - Mm).det()), xx)
    cc = [sp.expand(cp.coeff_monomial(xx ** k)) for k in (2, 1, 0)]

    def red(e):
        out = 0
        for (kk,), co in sp.Poly(sp.expand(e), cs).terms():
            out += co * (w * (1 - w)) ** (kk // 2)
        return sp.expand(out)
    cc = [red(e) for e in cc]
    E2s = sp.expand(27 * cc[2] - cc[0] ** 3)
    pe = sp.Poly(E2s, Xs, gs, ys)
    degs = {sum(m) for m, _ in pe.terms()}
    log('  E2 in (X, gamma, y): total degrees %s, %d monomials (printed: a cubic form with 10 terms)' % (sorted(degs), len(pe.terms())))
    ok &= degs == {3} and len(pe.terms()) == 10
    # E1 = (p-q)^2, E2 = (p-q)^2 (8p+q) on a double root p, third root q
    p_, q_ = sp.symbols('p q')
    cp2 = sp.Poly(sp.expand((al['x'] - p_) ** 2 * (al['x'] - q_)), al['x'])
    c2_, c1_, c0_ = [cp2.coeff_monomial(al['x'] ** k) for k in (2, 1, 0)]
    dd = sp.expand(c2_ ** 2 - 3 * c1_ - (p_ - q_) ** 2), sp.expand(27 * c0_ - c2_ ** 3 - (p_ - q_) ** 2 * (8 * p_ + q_))
    log('  double root p, third root q: E1 - (p-q)^2 and E2 - (p-q)^2 (8p + q) =', dd)
    ok &= dd == (0, 0)
    # resultant, P
    Pq = sp.Poly(al['Pq'], phi, A, B)
    fac = sp.factor_list(al['Pq'])
    degs = sorted({sum(m) for m, _ in Pq.terms()})
    log('  P: %d terms, total degrees %s, irreducible over Q: %s (printed: 167, 10..18, irreducible)'
        % (len(al['Ptab']), degs, len(fac[1]) == 1 and fac[1][0][1] == 1))
    ok &= len(al['Ptab']) == 167 and degs == [10, 12, 14, 16, 18] and len(fac[1]) == 1 and fac[1][0][1] == 1
    dres = sp.expand(al['Res'] - sp.Rational(27, 256) * phi ** 3 * (A - B) ** 3 * al['Pq'])
    log('  Res_w(E1, E2) - (27/256) phi^3 (A\' - B\')^3 P =', dres)
    ok &= dres == 0
    dlead = (sp.factor(sp.Poly(al['E1'], w).coeff_monomial(w ** 2)), sp.factor(sp.Poly(al['E2'], w).coeff_monomial(w ** 3)))
    log('  leading coefficients in w: E1: %s ; E2: %s' % dlead)
    # the interval code's tables against an exact Sylvester determinant at random rational points
    rng = random.Random(1360)
    badS = 0
    for _ in range(6):
        ph = Fr(rng.randint(1, 40), rng.randint(1, 9))
        Aa = Fr(rng.randint(-30, 30), rng.randint(1, 9))
        Bb = Fr(rng.randint(-30, 30), rng.randint(1, 9))

        def wcoef(tab, deg):
            co = [Fr(0)] * (deg + 1)
            for (a, i, j, k), c in tab.items():
                co[a] += c * ph ** i * Aa ** j * Bb ** k
            return co[::-1]            # highest power first
        f, h = wcoef(al['E1tab'], 2), wcoef(al['E2tab'], 3)
        S = [f + [0, 0], [0] + f + [0], [0, 0] + f, h + [0], [0] + h]
        det = fdet([[Fr(x) for x in row] for row in S])
        val = Fr(27, 256) * ph ** 3 * (Aa - Bb) ** 3 * al['P'].eval((ph, Aa, Bb), Fr(1))
        if det != val:
            badS += 1
    log('  exact 5 x 5 Sylvester determinant of the tables of E1, E2 against (27/256) phi^3 (A\' - B\')^3 P (table of the interval code)'
        ' at 6 rational points: %d mismatches' % badS)
    ok &= badS == 0
    return ok


def fdet(M):
    """Exact determinant of a square matrix of Fractions."""
    M = [row[:] for row in M]
    n, det = len(M), Fr(1)
    for i in range(n):
        piv = next((r for r in range(i, n) if M[r][i] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != i:
            M[i], M[piv] = M[piv], M[i]
            det = -det
        det *= M[i][i]
        for r in range(i + 1, n):
            fct = M[r][i] / M[i][i]
            for c in range(i, n):
                M[r][c] -= fct * M[i][c]
    return det


def sample_ok(iv_tuple, mpvals, slack=40):
    return all(_mp_enc_contains(e, t, slack) for e, t in zip(iv_tuple, mpvals))


@check('planted')
def check_planted():
    import mpmath as mp
    mp.mp.dps = 60
    ok = True
    dbl = MUT == 'planted-double'
    # A: g(z) = prod (cosh z - c_i), c_i rational, root i at arccosh(c_i); two roots 1e-5 apart, two outside [1/4, 40]
    cs = planted_c()
    roots = [mp.acosh(mpf_of(c)) for c in cs]
    lo_c, hi_c = cosh_const(Fr(1, 4)), cosh_const(40)
    inside = [c for c in cs if c > lo_c.hi and c < hi_c.lo]
    exact_n = len(inside)
    gap = abs(roots[2] - roots[1])
    log('  A: g has %d roots in [1/4, 40] (exact: %d of %d rational c_i lie strictly between the enclosures of cosh(1/4) and cosh(40));'
        ' pair separation %s' % (exact_n, exact_n, len(cs), mp.nstr(gap, 6)))
    ok &= exact_n == 7 and Fr(99, 10 ** 7) < Fr(str(mp.nstr(gap, 10))) < Fr(101, 10 ** 7)
    cov = Cover(lambda z: planted_g(z, cs, dbl), Fr(1, 4), Fr(40))
    t0 = time.time()
    leaves = cov.run()
    kinds = {}
    for lf in leaves:
        kinds[lf[2]] = kinds.get(lf[2], 0) + 1
    rl = [lf for lf in leaves if lf[2] == 'root']
    log('  A: %d leaves %s, tiling %s, %.0f s' % (len(leaves), kinds, tiling_ok(leaves, Fr(1, 4), Fr(40)), time.time() - t0))
    ok &= tiling_ok(leaves, Fr(1, 4), Fr(40)) and 'open' not in kinds and len(rl) == exact_n
    exp_in = sorted(r for r in roots if mp.mpf(1) / 4 < r < 40)
    hit = [sum(1 for r in exp_in if mpf_of(lf[0]) < r < mpf_of(lf[1])) for lf in rl]
    log('  A: every root leaf holds exactly one known root: %s' % (hit == [1] * len(rl) and len(rl) == len(exp_in)))
    ok &= hit == [1] * len(rl) and len(rl) == len(exp_in)
    # B: planted (E1, E2) with common roots at known (z, w): classification, Krawczyk boxes, close pair, T = 0, w outside (0, 1)
    pm = PlantedModel()
    cov = Cover(lambda z: pm.P(z), Fr(1, 4), Fr(40))
    leaves = cov.run()
    kinds = {}
    for lf in leaves:
        kinds[lf[2]] = kinds.get(lf[2], 0) + 1
    rl = [lf for lf in leaves if lf[2] == 'root']
    exp_pts = sorted([Fr(3), Fr(659, 49), Fr(5), Fr(5) + PlantedModel.EPS, Fr(7)])
    exp_w = [pm.uv(Du(z))[1 if z == 7 else 0].v.lo for z in exp_pts]         # u1(z) (u2(z) at z = 7): a start for the box
    exp_adm = {Fr(3): True, Fr(659, 49): True, Fr(5): True, Fr(5) + PlantedModel.EPS: False, Fr(7): False}
    log('  B: P_pl has %d leaves %s, %d roots (exact: 5 at z = 3, 5, 5 + 1e-5, 7, 659/49)' % (len(leaves), kinds, len(rl)))
    ok &= tiling_ok(leaves, Fr(1, 4), Fr(40)) and 'open' not in kinds and len(rl) == 5
    adm = {}
    with bits_ctx(320):
        for lf, zx, wx in zip(rl, exp_pts, exp_w):
            Z = refine_root(lambda z: pm.P(z), lf[0], lf[1])
            good = Z.contains(zx)
            boxes, complete = common_roots(pm, Z, Fr(-10 ** 4), Fr(10 ** 4), [(0, wx)])
            good &= complete and len(boxes) == 1 and boxes[0][1].contains(wx)
            if boxes:
                Zb, Wb = boxes[0]
                Tv = pm.T(Du(Zb)).v
                adm[zx] = (Wb.lo > 0 and Wb.hi < 1) and not Tv.contains0()
            log('     root z = %s: enclosure %s, boxes %d, complete %s' % (fl(zx, 12), 'contains it' if Z.contains(zx) else 'MISSES it', len(boxes), complete))
            ok &= good
    log('  B: admissible (0 < w < 1, T != 0) at z = 3, 5, 5 + 1e-5, 7, 659/49: %s; exact expectation [True, True, False (T = 0), False (w = 7/2), True]'
        % [adm.get(z) for z in exp_pts])
    ok &= adm == exp_adm
    return ok


def _wall_checks(prefix, fn):
    for wn in ('le', 'rodal'):
        CHECKS['%s:%s' % (prefix, wn)] = (lambda wn=wn: fn(wn))


def small_zs(wn):
    if MUT == 'small-z-range':
        return Fr(20) if wn == 'le' else Fr(45)           # control: beyond the first roots of y and P
    return Fr(1, 4)


def check_small_z(wn):
    import sympy as sp
    import mpmath as mp
    mp.mp.dps = 100
    W = walls()[wn]
    ok = True
    gt = G_table()
    log('  P(phi, psi - 2 phi, delta + psi - 2 phi): %d terms, least weight %d (phi: 2, psi, delta: 4), weight-%d part %s'
        % (len(gt['tab']), gt['minw'], gt['minw'], {m: int(c) for m, c in gt['low'].items()}))
    ok &= gt['minw'] == 24 and gt['low'] == {(8, 0, 2): Fr(774144)}
    # closed forms
    z, C, v = sp.symbols('z C v', positive=True)
    c, s = sp.symbols('c s')
    K1, S = (c - 1) / z ** 2, s / z
    K2, K3 = (S - 2 * K1) / z ** 2, (S - 1) / z ** 2
    D = 1 + C + z ** 2 * K1
    phi_ = v * (c - 1) / (c + C)
    A_ = -v * z * (1 + C) * s / (c + C) ** 2
    B_ = -v * z ** 2 * (1 + C) * (-c ** 2 + C * c + 2) / (c + C) ** 3
    Phi_f = v * K1 / D
    Psi_f = v * (2 * K1 ** 2 - (1 + C) * K2) / D ** 2
    Del_f = -v * (1 + C) * (-(1 + C) * K3 + (C - 2) * K1 - S * K1 - z ** 2 * K1 ** 2) / D ** 3
    dd = [sp.simplify(sp.together(z ** 2 * Phi_f - phi_)), sp.simplify(sp.together(z ** 4 * Psi_f - (A_ + 2 * phi_))),
          sp.simplify(sp.together(z ** 4 * Del_f - (B_ - A_)))]
    log('  closed forms Phi, Psi, Delta of the weights (exact, sympy): z^2 Phi - phi, z^4 Psi - (A\' + 2 phi), z^4 Delta - (B\' - A\') =', dd)
    ok &= all(d == 0 for d in dd)
    zz = sp.Symbol('zz')
    ser = {'K1': (sp.cosh(zz) - 1) / zz ** 2, 'S': sp.sinh(zz) / zz,
           'K2': (sp.sinh(zz) / zz - 2 * (sp.cosh(zz) - 1) / zz ** 2) / zz ** 2, 'K3': (sp.sinh(zz) / zz - 1) / zz ** 2}
    okser = True
    for name, expr in ser.items():
        tay = sp.series(expr, zz, 0, 22).removeO()
        for n in range(11):
            okser &= sp.Rational(SERIES[name](n).numerator, SERIES[name](n).denominator) == tay.coeff(zz, 2 * n)
    log('  series coefficients of K1, S, K2, K3 against the Taylor expansions of their closed forms (11 coefficients each):', okser)
    ok &= okser
    # numerical cross-check at points
    bad = 0
    for zt_ in (Fr(1, 1000), Fr(1, 20), Fr(1, 8), Fr(1, 4)):
        Zsq, Phi, Psi, Delta = small_z_forms(zt_, zt_, W)
        ph, A, B = mp_cosh_form(mpf_of(zt_), W)
        zm = mpf_of(zt_)
        for e, t in ((Phi, ph / zm ** 2), (Psi, (A + 2 * ph) / zm ** 4), (Delta, (B - A) / zm ** 4)):
            if not _mp_enc_contains(e, t, 50):
                bad += 1
    log('  Phi, Psi, Delta enclosures against mpmath (100 digits) at z = 1/1000, 1/20, 1/8, 1/4: %d misses' % bad)
    ok &= bad == 0
    zs = small_zs(wn)
    t0 = time.time()
    cells, badc = small_z_cover(W, zs, maxdepth=14 if MUT != 'small-z-range' else 6)
    if cells:
        log('  [0, %s]: %d cells, G > 0 and Delta < 0 on each (smallest lower bound of G %s, largest upper bound of Delta %s); %d cells fail; %.1f s'
            % (fl(zs, 4), len(cells), fl(min(c_[2].lo for c_ in cells), 5), fl(max(c_[3].hi for c_ in cells), 5), len(badc), time.time() - t0))
    else:
        log('  [0, %s]: no cell proves G > 0 and Delta < 0; %d cells fail' % (fl(zs, 4), len(badc)))
    ok &= not badc and bool(cells)
    log('  so P = z^24 G > 0 and y = z^4 Delta < 0 on (0, %s]' % fl(zs, 4))
    return ok


def tail_zt(wn):
    if MUT == 'tail-zt':
        return 4                                           # control: a start where the dominance does not hold
    return walls()[wn].zt


def check_tail(wn):
    import sympy as sp
    W = walls()[wn]
    ok = True
    z, C, v, e = sp.symbols('z C v e', positive=True)
    c, s = (1 + e ** 2) / (2 * e), (1 - e ** 2) / (2 * e)
    phi_ = v * (c - 1) / (c + C)
    A_ = -v * z * (1 + C) * s / (c + C) ** 2
    B_ = -v * z ** 2 * (1 + C) * (-c ** 2 + C * c + 2) / (c + C) ** 3
    den = 1 + 2 * C * e + e ** 2
    tau = 2 * v * (1 + C) * z * e
    ph_h = (1 - e) ** 2 / den
    al_h = (1 - e ** 2) / den ** 2
    be_h = ((1 + e ** 2) ** 2 - 2 * C * e * (1 + e ** 2) - 8 * e ** 2) / den ** 3
    dd = [sp.simplify(sp.together(phi_ - v * ph_h)), sp.simplify(sp.together(A_ + tau * al_h)),
          sp.simplify(sp.together(B_ - tau * z * be_h))]
    log('  with e = e^-z, cosh = (1 + e^2)/(2e), sinh = (1 - e^2)/(2e): phi - v phih, A\' + tau alphah, B\' - tau z betah (exact) =', dd)
    ok &= all(d == 0 for d in dd)
    rng = random.Random(1361)
    bad = 0
    tab = tail_table(W)
    for _ in range(5):
        phh, alh, beh = [Fr(rng.randint(5, 20), 10) for _ in range(3)]
        tau_, z_ = Fr(rng.randint(1, 50), 7), Fr(rng.randint(30, 900), 10)
        lhs = algebra()['P'].eval((W.v * phh, -tau_ * alh, tau_ * z_ * beh), Fr(1))
        rhs = (tau_ * z_) ** 9 * tab.eval((1 / z_, tau_, phh, alh, beh), Fr(1))
        bad += lhs != rhs
    log('  P(v phih, -tau alphah, tau z betah) = (tau z)^9 S(1/z, tau, phih, alphah, betah) at 5 random rational points (exact): %d mismatches' % bad)
    ok &= bad == 0
    zt = tail_zt(wn)
    t0 = time.time()
    tb = tail_bounds(W, zt)
    smin = min(c_[2].lo for c_ in tb['cells'])
    log('  z >= %s: e <= %s, tau <= %s, phih in [%s, %s], betah in [%s, %s]'
        % (zt, fl(tb['e_t'], 4), fl(tb['tau_t'], 4), fl(tb['phih'].lo, 8), fl(tb['phih'].hi, 3), fl(tb['beth'].lo, 8), fl(tb['beth'].hi, 3)))
    log('  S = P/(tau z)^9 over %d cells of s = 1/z in (0, 1/%s]: smallest lower bound %s (leading term 4 v^3 = %s); %.1f s'
        % (len(tb['cells']), zt, fl(smin, 8), fl(4 * W.v ** 3, 6), time.time() - t0))
    ok &= smin > 0
    ypos = tb['beth'].lo > 0 and tb['alph'].lo > 0
    log('  y = tau (z betah + alphah) > 0: betah, alphah lower bounds %s, %s -> %s' % (fl(tb['beth'].lo, 6), fl(tb['alph'].lo, 6), ypos))
    ok &= ypos
    # the enclosure of S contains the true S at points z >= zt
    bad = 0
    import mpmath as mp
    mp.mp.dps = 120
    for zk in (zt, 2 * zt, 3 * zt, 5 * zt):
        zm = mpf_of(zk)
        ph, A, B = mp_cosh_form(zm, W)
        Pm = mp_poly(algebra()['Ptab'], (ph, A, B))
        tm = 2 * mpf_of(W.v) * (1 + mp.cosh(mpf_of(2 * W.sigma * W.Rb))) * zm * mp.exp(-zm)
        Sm = Pm / (tm * zm) ** 9
        cell = next(c_ for c_ in tb['cells'] if c_[0] <= Fr(1, zk) <= c_[1])
        if not _mp_enc_contains(cell[2], Sm):
            bad += 1
    log('  the enclosure of S contains the true S (mpmath, 120 digits) at z = 1, 2, 3, 5 times %s: %d misses' % (zt, bad))
    ok &= bad == 0
    return ok


_wall_checks('small_z', check_small_z)
_wall_checks('tail', check_tail)


_tanh_lam = None


def mp_tanh_form(z, wall):
    """(phi, A', B') in mp straight from the tanh definition of f (derivatives by sympy), an independent route."""
    global _tanh_lam
    import sympy as sp
    if _tanh_lam is None:
        r_, sg, Rb = sp.symbols('r sigma Rb')
        ft = (sp.tanh(sg * (r_ + Rb)) - sp.tanh(sg * (r_ - Rb))) / (2 * sp.tanh(sg * Rb))
        _tanh_lam = sp.lambdify((r_, sg, Rb), (ft, sp.diff(ft, r_), sp.diff(ft, r_, 2)), 'mpmath')
    r = z / (2 * mpf_of(wall.sigma))
    f, fa, fb = _tanh_lam(r, mpf_of(wall.sigma), mpf_of(wall.Rb))
    v = mpf_of(wall.v)
    return v * (1 - f), v * r * fa, v * r ** 2 * fb


def mp_tanh_P(z, wall):
    return mp_poly(algebra()['Ptab'], mp_tanh_form(z, wall))


def mp_tanh_y(z, wall):
    _, A, B = mp_tanh_form(z, wall)
    return B - A


def mp_bisect(f, a, b, n=110):
    fa = f(a)
    for _ in range(n):
        m = (a + b) / 2
        fm = f(m)
        if (fa < 0) == (fm < 0):
            a, fa = m, fm
        else:
            b = m
    return (a + b) / 2


def get_cover(wn):
    """Cover of P and of y over [1/4, zt] for a wall, with refined root enclosures. Cached in STATE."""
    key = ('cover', wn)
    if key in STATE:
        return STATE[key]
    W = walls()[wn]
    lo, hi = Fr(1, 4), Fr(W.zt)
    out = {}
    for nm, fun in (('P', lambda z: P_du(z, W)), ('y', lambda z: y_du(z, W))):
        t0 = time.time()
        cv = Cover(fun, lo, hi)
        leaves = cv.run()
        if MUT == 'drop-leaf':
            leaves = leaves[:len(leaves) // 2] + leaves[len(leaves) // 2 + 1:]       # control: one leaf is lost
        roots = [lf for lf in leaves if lf[2] == 'root']
        with bits_ctx(320):
            enc = [refine_root(fun, lf[0], lf[1]) for lf in roots]
        out[nm] = dict(leaves=leaves, roots=roots, enc=enc, time=time.time() - t0, nevals=cv.nevals)
    STATE[key] = out
    return out


def check_cover(wn):
    import mpmath as mp
    mp.mp.dps = 80
    W = walls()[wn]
    lo, hi = Fr(1, 4), Fr(W.zt)
    cv = get_cover(wn)
    ok = True
    nexp = {'le': 4, 'rodal': 6}[wn]
    for nm, ne in (('P', nexp), ('y', 1)):
        d = cv[nm]
        kinds = {}
        for lf in d['leaves']:
            kinds[lf[2]] = kinds.get(lf[2], 0) + 1
        tl = tiling_ok(d['leaves'], lo, hi)
        log('  %s on [1/4, %s]: %d leaves %s; tiling %s; %d roots (statement: %d); %.0f s' % (nm, W.zt, len(d['leaves']), kinds, tl, len(d['roots']), ne, d['time']))
        ok &= tl and 'open' not in kinds and len(d['roots']) == ne
        fun = (lambda z: mp_tanh_P(z, W)) if nm == 'P' else (lambda z: mp_tanh_y(z, W))
        # sign of the tanh-form function (independent route) at the midpoint of every leaf agrees with the leaf's sign
        bad = 0
        for a, b, kind, info in d['leaves']:
            if kind in ('excl', 'mono'):
                sg_ = info
                val = fun(mpf_of((a + b) / 2))
                if (val > 0) != (sg_ > 0):
                    bad += 1
            elif kind == 'root':
                if (fun(mpf_of(a)) > 0) != (info[0] > 0) or (fun(mpf_of(b)) > 0) != (info[1] > 0):
                    bad += 1
        log('     tanh-form function (mpmath, 80 digits): sign at a point of every leaf vs the leaf sign: %d disagreements' % bad)
        ok &= bad == 0
        # root enclosures contain the mpmath root of the tanh form
        miss = 0
        rr = []
        for lf, enc in zip(d['roots'], d['enc']):
            zr = mp_bisect(fun, mpf_of(lf[0]), mpf_of(lf[1]))
            if not (mpf_of(enc.lo) - mp.mpf(10) ** -28 <= zr <= mpf_of(enc.hi) + mp.mpf(10) ** -28):
                miss += 1
            rr.append(W.r_of_z(enc.mid()))
        log('     refined root enclosures (width <= %s) contain the mpmath root of the tanh form: %d misses; r = %s'
            % (fl(max(e.width() for e in d['enc']), 3) if d['enc'] else '-', miss, [fl(r_, 10) for r_ in rr]))
        ok &= miss == 0
    return ok


def check_points(wn):
    import mpmath as mp
    mp.mp.dps = 80
    W = walls()[wn]
    cv = get_cover(wn)
    model = WallModel(W)
    ok = True
    adm_w = (Fr(0), Fr(100) if MUT == 'w-range' else Fr(1))
    pts = []
    t0 = time.time()
    with bits_ctx(320):
        for enc in cv['P']['enc']:
            zm = mpf_of(enc.mid())
            cands = mp_candidates(W, zm)
            boxes, complete = common_roots(model, enc, Fr(-10 ** 4), Fr(10 ** 4), cands)
            ok &= complete
            if not boxes:
                log('     root r = %s: no common root with |w| <= 1e4 verified (complete %s)' % (fl(W.r_of_z(enc.mid()), 10), complete))
            for Zb, Wb in boxes:
                Tv = model.T(Du(Zb)).v
                tz = 4 * W.sigma ** 2 * Tv / (W.v * Zb.sq())              # T = (A' + B') / (v r^2), r = z / (2 sigma)
                inside = Wb.lo > adm_w[0] and Wb.hi < adm_w[1]
                nz = not Tv.contains0()
                pts.append(dict(z=Zb, w=Wb, T=tz, adm=inside and nz, inside=inside, nz=nz, r=W.r_of_z(Zb.mid())))
    adm = [p for p in pts if p['adm']]
    oth = [p for p in pts if not p['adm']]
    log('  %s: %d P roots, %d Krawczyk boxes in (z, w) (K(X) inside X), every other point of z in the root cell x w in [-1e4, 1e4] excluded: %s; %.0f s'
        % (wn, len(cv['P']['enc']), len(pts), ok, time.time() - t0))
    log('  admissible (0 < w < 1, T != 0): %d (statement: 3)' % len(adm))
    ok &= len(adm) == 3
    for p, rs, ws, Ts in zip(adm, W.radii, W.ws, W.Ts):
        dr = abs(p['r'] - rs)
        good = dr <= Fr(1, 10 ** 9) and abs(p['w'].mid() - ws) <= Fr(1, 10 ** 6) and abs(p['T'].mid() - Ts) <= Fr(1, 10 ** 4)
        log('     r = %s (statement %s, |diff| %s), w in [%s, %s], T = %s; box half-widths z %s, w %s -> %s'
            % (fl(p['r'], 12), fl(rs, 12), fl(dr, 2), fl(p['w'].lo, 10), fl(p['w'].hi, 10), fl(p['T'].mid(), 6), fl(p['z'].rad(), 2), fl(p['w'].rad(), 2),
               'ok' if good else 'MISMATCH'))
        ok &= good
    ows = sorted(p['w'].mid() for p in oth)
    exp_ows = sorted(W.other_w)
    log('  the other P roots have their common root at w = %s (statement: %s)' % ([fl(x, 8) for x in ows], [str(x) for x in exp_ows]))
    ok &= len(ows) == len(exp_ows) and all(abs(a - b) <= Fr(1, 1000) for a, b in zip(ows, exp_ows))
    return ok


def check_spurious(wn):
    W = walls()[wn]
    cv = get_cover(wn)
    model = WallModel(W)
    ok = True
    if MUT == 'spurious-on-P-root':
        Z = cv['P']['enc'][0]                                 # control: a root of P, where a common root exists
        what = 'the first root of P'
    else:
        Z = cv['y']['enc'][0] if cv['y']['enc'] else None
        what = 'the root of y'
    if Z is None:
        log('  y has no root')
        return False
    with bits_ctx(320):
        r_ = W.r_of_z(Z.mid())
        ph, A, B = phiAB(Du(Z), W)
        lead = ph.v.sq() * (B - A).v.sq() * Fr(1, 4)
        log('  %s at r = %s (statement %s): leading coefficient of E1 in w, phi^2 y^2 / 4, encloses %s .. %s' % (what, fl(r_, 8), fl(W.y_r, 5), fl(lead.lo, 3), fl(lead.hi, 3)))
        boxes, complete = common_roots(model, Z, Fr(-10 ** 4), Fr(10 ** 4), [])
    log('  common roots with z in the root cell and |w| <= 1e4: %d (complete %s)' % (len(boxes), complete))
    ok &= len(boxes) == 0 and complete and abs(r_ - W.y_r) <= Fr(5, 10 ** 6)
    return ok


_wall_checks('cover', check_cover)
_wall_checks('points', check_points)
_wall_checks('spurious', check_spurious)


# ---------------------------------------------------------------------------------------------------------------------
# 12. Driver
# ---------------------------------------------------------------------------------------------------------------------

ORDER = ['exp', 'profile', 'frame', 'planted', 'small_z:le', 'small_z:rodal', 'tail:le', 'tail:rodal',
         'cover:le', 'points:le', 'spurious:le', 'cover:rodal', 'points:rodal', 'spurious:rodal']


def main():
    global MUT, BITS
    ap = argparse.ArgumentParser()
    ap.add_argument('--checks', default='')
    ap.add_argument('--mutate', default=None)
    ap.add_argument('--bits', type=int, default=BITS)
    args = ap.parse_args()
    MUT, BITS = args.mutate, args.bits
    names = []
    for n in (args.checks.split(',') if args.checks else ORDER):
        names += [m for m in ORDER if m == n or m.split(':')[0] == n] if n else []
    log('review_he_type3_count.py  mutation: %s  bits: %d' % (MUT, BITS))
    fails, t00 = [], time.time()
    for n in names:
        t0 = time.time()
        log('check %s' % n)
        try:
            ok = bool(CHECKS[n]())
        except Exception as ex:
            tb = traceback.extract_tb(sys.exc_info()[2])[-1]
            log('  EXCEPTION %s: %s (line %d)' % (type(ex).__name__, ex, tb.lineno))
            ok = False
        log('  -> %s (%.0f s)' % ('PASS' if ok else 'FAIL', time.time() - t0))
        if not ok:
            fails.append(n)
    log('')
    log('summary: %d checks run, %d failed %s, %.0f s' % (len(names), len(fails), fails, time.time() - t00))
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()

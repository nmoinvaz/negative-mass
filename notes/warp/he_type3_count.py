#!/usr/bin/env python3
"""A rigorous count of the Type III points of the tanh wall at Le's and Rodal's speeds (issue #133).

Setting. he_type3.md: off the equator and axis, with v c s u != 0 and (a, b) != (0, 0), a point of the Alcubierre wall is
Type III exactly when E1 = E2 = 0 and T = b + a u != 0. With (phi, A', B') = (v F, v r f', v r^2 f'') and w = c^2, E1 and E2
are polynomials in w (degrees 2 and 3) with coefficients in (phi, A', B') (u = 1 by the weights 4 and 6). Their resultant in w
is (27/256) phi^3 (A' - B')^3 P with P irreducible (167 terms). So an admissible point at radius r has P = 0 or y = B' - A' = 0
(phi > 0 for r > 0). The script isolates every root of P and of y in r > 0 and decides, at each, the common roots w.

Imports. The tables of E1, E2 and P come from `he_type3.shadow()` (imported, not rebuilt): E1 and E2 from the frame Ricci block
of `review_he_surface_lib.py`, P from sympy's resultant. `tables` re-derives their identity independently: an exact Sylvester
determinant in fractions, and E1, E2 from the characteristic polynomial of the frame table printed in he_type3.md section 1.

Variables. z = 2 sigma r, eps = e^-z, C = cosh 2 sigma R_b, D = 1 + 2 C eps + eps^2. The profile is f = 2 (1 + C) eps / D and

    phi = v (1 - eps)^2 / D,   A' = -t a~,   B' = z t b~,   t = 2 v (1 + C) z eps,
    a~ = (1 - eps^2) / D^2,    b~ = [(1 + eps^2)^2 - 2 C eps (1 + eps^2) - 8 eps^2] / D^3,

so P = t^9 g with g = sum cf (-1)^j phi^i a~^j b~^k z^k t^(j + k - 9) (every term has j + k >= 9), and y = B' - A' = t Yn / D^3 with
Yn = z [(1 + eps^2)^2 - 2 C eps (1 + eps^2) - 8 eps^2] + (1 - eps^2) D. Ranges of z:

  (0, z_s]      P = z^24 G with G = 774144 Phi^8 Delta^2 + (higher weight) in phi = z^2 Phi, psi = A' + 2 phi = z^4 Psi,
                delta = B' - A' = z^4 Delta, each an even series in z with an enclosed remainder; G and Delta exclude 0.
  [z_s, z_t]    interval cover of g and of Yn: excluded leaves, leaves with one root (sign-definite derivative, opposite signs
                at the ends), no unresolved leaf. Forward-mode derivatives on mpmath `iv` intervals.
  [z_t, inf)    dominance of the z^9 coefficient 4 phi^3 b~^9 of the (A', B') degree 9 part; Yn > 0.

At each root cell the common roots (z, w) with w in a wide range are isolated by bisection and the Krawczyk operator, and the
admissible ones (0 < w < 1 and T != 0 on the box) are counted. A planted problem with a known number of admissible roots (sympy
Sturm count, a close pair, roots with w outside (0, 1) and with T = 0) runs through the same cover and classification.

    python3 notes/warp/he_type3_count.py                  all checks, then every control (exit 1 on any failure)
    python3 notes/warp/he_type3_count.py --mutate=NAME    one control, all checks (exit 1 if any check fails)
    python3 notes/warp/he_type3_count.py --no-controls    the checks only
    python3 notes/warp/he_type3_count.py --checks=a,b --walls=Le    a subset of checks, optionally one wall

Trusted: mpmath 1.3.0 `iv` (outward-rounded + - * / and exp), sympy 1.14.0 (resultant, Poly, ring, Sturm count), Python's exact
Fraction, and he_type3's tables (checked in `tables`).
"""
import math
import os
import re
import subprocess
import sys
import time
from fractions import Fraction as Fr

import mpmath as mp
import sympy as sp
from mpmath import iv
from sympy.polys.rings import ring
from sympy import QQ

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import he_type3 as H3  # noqa: E402

iv.dps = 50
mp.mp.dps = 60
T0 = time.time()
MUT = ""
HERE = os.path.dirname(os.path.abspath(__file__))

# name: (checks that must fail, what it changes). The targets were written to the raw README before the first run.
CONTROLS = {
    "sigma-1e-3": (["points", "cross"], "sigma + 1e-3 at both walls (C and the radii; v fixed)"),
    "drop-P-term": (["tables", "small_z", "points", "cross", "assembly"], "the first term of P's table is dropped everywhere"),
    "no-tail-large": (["tail_z", "planted", "assembly"], "the proof on [z_t, infinity) is removed"),
    "no-tail-small": (["small_z", "assembly"], "the proof on (0, z_s] is removed"),
    "cell-gap": (["cover", "assembly"], "one leaf is removed from the tiling of each wall's cover"),
    "eps-typo": (["identity", "cover", "spurious", "points", "cross", "assembly"], "D = 1 + 2 C eps + eps^2 becomes 1 + C eps + eps^2 in the jets"),
    "T-off": (["planted"], "the T != 0 condition is dropped from admissibility"),
    "w-range": (["planted", "points"], "admissible w in (0, 1) becomes (-10^4, 10^4)"),
    "krawczyk-nocorr": (["planted", "points"], "the Krawczyk operator drops its (I - Y J)(X - m) term"),
    "double-root": (["planted"], "the planted close pair becomes a double root"),
    "skip-y": (["planted", "spurious", "assembly"], "the roots of the spurious factor y are not classified"),
}
W_DIAG = 10**4  # the diagnostic range of w for the common roots


def say(msg):
    print(msg, flush=True)


# ----------------------------------------------------------------------------------------------------------------------
# interval helpers and forward-mode derivatives
# ----------------------------------------------------------------------------------------------------------------------
def ivpt(x):
    """The degenerate interval [x, x] (x an mp.mpf; exact when x has at most the interval precision)."""
    return iv.mpf([x, x])


def ivr(a, b):
    return iv.mpf([a, b])


def mid(I):
    return (mp.mpf(I.a) + mp.mpf(I.b)) / 2


def width(I):
    return mp.mpf(I.b) - mp.mpf(I.a)


def exc0(I):
    """True if the interval excludes 0."""
    return I.a > 0 or I.b < 0


class Du:
    """A pair (value, derivative in z) of intervals. Operands: Du with Du, or Du on the left with an iv scalar."""
    __slots__ = ("v", "d")

    def __init__(s, v, d):
        s.v, s.d = v, d

    def __add__(s, o):
        return Du(s.v + o.v, s.d + o.d) if isinstance(o, Du) else Du(s.v + o, s.d)

    def __sub__(s, o):
        return Du(s.v - o.v, s.d - o.d) if isinstance(o, Du) else Du(s.v - o, s.d)

    def __mul__(s, o):
        if isinstance(o, Du):
            return Du(s.v * o.v, s.v * o.d + s.d * o.v)
        return Du(s.v * o, s.d * o)

    def __truediv__(s, o):
        if isinstance(o, Du):
            iq = 1 / o.v
            return Du(s.v * iq, (s.d * o.v - s.v * o.d) * iq * iq)
        return Du(s.v / o, s.d / o)

    def __neg__(s):
        return Du(-s.v, -s.d)


def du_const(x):
    return Du(iv.mpf(x), iv.mpf(0))


def powers(x, n):
    out = [du_const(1)]
    for _ in range(n):
        out.append(out[-1] * x)
    return out


_CF = {}


def cfi(p, q):
    """The rational p/q as an interval."""
    key = (p, q)
    if key not in _CF:
        _CF[key] = iv.mpf(p) / iv.mpf(q)
    return _CF[key]


# ----------------------------------------------------------------------------------------------------------------------
# the walls and the profile
# ----------------------------------------------------------------------------------------------------------------------
class Wall:
    """sigma, R_b, v exact; C = cosh 2 sigma R_b as an interval; z-ranges of the proof."""

    def __init__(self, name, sigma, Rb, v, zs, zt, table, mut_sigma=False):
        self.name = name
        self.sig = Fr(sigma) + (Fr(1, 1000) if mut_sigma else 0)
        self.Rb, self.v0 = Fr(Rb), Fr(v)
        self.zs, self.zt = Fr(zs), Fr(zt)
        self.table = table
        self.v = iv.mpf(self.v0.numerator) / self.v0.denominator
        x = iv.mpf(2 * self.sig.numerator * self.Rb.numerator) / (self.sig.denominator * self.Rb.denominator)
        self.p2 = x  # 2 sigma R_b
        self.C = (iv.exp(x) + iv.exp(-x)) / 2
        self.mps = mp.mpf(self.sig.numerator) / self.sig.denominator

    def Tn_of(self, Z):
        """A' + B' over the interval Z (T = this / (v r^2))."""
        J = jets_real(Z, self)
        return (J["A"] + J["B"]).v

    def r_iv(self, Z):
        """The radius r = z / (2 sigma) of an interval Z."""
        return Z / (2 * iv.mpf(self.sig.numerator) / self.sig.denominator)


def make_walls(names=None):
    ms = MUT == "sigma-1e-3"
    walls = {
        "Le": Wall("Le", 8, 1, Fr(1, 2), Fr(1, 4), 40, [0.990568, 1.023933, 1.094662], ms),
        "Rodal": Wall("Rodal", 4, 5, 1, Fr(1, 4), 60, [4.714485, 5.071066, 5.103381], ms),
    }
    return {k: w for k, w in walls.items() if names is None or k in names}


def jets_real(X, wl):
    """phi, A', B', t, z, (a~, b~), D as Du over the interval X, from the (eps, z) forms."""
    C, v = wl.C, wl.v
    z = Du(X, iv.mpf(1))
    e = iv.exp(-X)
    eps = Du(e, -e)
    one = du_const(1)
    e2 = eps * eps
    Dt = one + eps * (C if MUT == "eps-typo" else 2 * C) + e2
    omE = one - eps
    phi = (omE * omE) * v / Dt
    Dt2 = Dt * Dt
    anum = one - e2
    at = anum / Dt2
    ope2 = one + e2
    bnum = ope2 * ope2 - eps * ope2 * (2 * C) - e2 * 8
    bt = bnum / (Dt2 * Dt)
    t = z * eps * (2 * v * (1 + C))
    A = -(t * at)
    B = z * t * bt
    return dict(phi=phi, A=A, B=B, t=t, z=z, at=at, bt=bt, bnum=bnum, anum=anum, Dt=Dt, eps=eps)


# ----------------------------------------------------------------------------------------------------------------------
# tables
# ----------------------------------------------------------------------------------------------------------------------
_TAB = {}


def tables():
    """P terms (i, j, k, p, q) in (phi, A', B'), E1 and E2 terms per power of w, and the exact sympy objects."""
    if "t" in _TAB:
        return _TAB["t"]
    sh = H3.shadow()
    PP = sp.Poly(sh["P"], H3.ph, H3.Ap, H3.Bp)
    pterms = [(i, j, k, int(cf.p), int(cf.q)) for (i, j, k), cf in PP.terms()]
    Eterms = [[[(mn[0], mn[1], mn[2], p, q) for mn, p, q in co] for co in sh["eexact"][kk]] for kk in range(2)]
    _TAB["t"] = dict(sh=sh, P=pterms, E=Eterms)
    return _TAB["t"]


def P_table():
    pt = tables()["P"]
    if MUT == "drop-P-term":
        return pt[1:]
    return pt


def group_P(terms):
    """Terms grouped by (j, k): {(j, k): [(i, cf interval)]} with the sign (-1)^j folded in, for g."""
    gr = {}
    for i, j, k, p, q in terms:
        sg = -1 if j % 2 else 1
        gr.setdefault((j, k), []).append((i, cfi(sg * p, q)))
    return gr


_GP = {}


def g_fun(wl):
    """The Du function X -> g(X) = P / t^9 for the wall (jets from the (eps, z) forms)."""
    key = (wl.name, MUT)
    if key not in _GP:
        _GP[key] = group_P(P_table())
    gr = _GP[key]
    imax = max(i for lst in gr.values() for i, _ in lst)
    kmax = max(k for _, k in gr)
    jmax = max(j for j, _ in gr)

    def g(X):
        J = jets_real(X, wl)
        pphi = powers(J["phi"], imax)
        pat, pbt = powers(J["at"], jmax), powers(J["bt"], kmax)
        pz, pt = powers(J["z"], kmax), powers(J["t"], max(j + k - 9 for j, k in gr))
        tot = du_const(0)
        for (j, k), lst in gr.items():
            h = du_const(0)
            for i, cf in lst:
                h = h + pphi[i] * cf
            tot = tot + h * pat[j] * pbt[k] * pz[k] * pt[j + k - 9]
        return tot
    return g


def yn_fun(wl):
    """Yn(z) = z bnum + (1 - eps^2) D (the sign of B' - A', which equals t Yn / D^3), as a Du function."""
    def yn(X):
        J = jets_real(X, wl)
        return J["z"] * J["bnum"] + J["anum"] * J["Dt"]
    return yn


# ----------------------------------------------------------------------------------------------------------------------
# generic cover and Newton shrink on an interval Du function
# ----------------------------------------------------------------------------------------------------------------------
def cover(fun, lo_, hi_, maxdepth=64):
    """Adaptive bisection of [lo_, hi_] for a Du function. Leaves (a, b, kind): 'x' no root (the enclosure, or the
    mean-value form, excludes 0, or the derivative is sign-definite and the ends have one sign), 'r' exactly one root (the
    derivative is sign-definite and the ends have opposite signs, each excluding 0), 'u' unresolved at maxdepth."""
    stack = [(mp.mpf(lo_), mp.mpf(hi_), 0)]
    leaves, nev = [], 0
    while stack:
        a, b, dep = stack.pop()
        X = ivr(a, b)
        r = fun(X)
        nev += 1
        G, dG = r.v, r.d
        if exc0(G):
            leaves.append((a, b, "x"))
            continue
        m = (a + b) / 2
        rm = fun(ivpt(m))
        nev += 1
        mv = rm.v + dG * (X - ivpt(m))
        if exc0(mv):
            leaves.append((a, b, "x"))
            continue
        if exc0(dG):
            ga, gb = fun(ivpt(a)).v, fun(ivpt(b)).v
            nev += 2
            if exc0(ga) and exc0(gb):
                if (ga.b < 0 < gb.a) or (ga.a > 0 > gb.b):
                    leaves.append((a, b, "r"))
                else:
                    leaves.append((a, b, "x"))
                continue
        if dep >= maxdepth:
            leaves.append((a, b, "u"))
            continue
        stack.append((m, b, dep + 1))
        stack.append((a, m, dep + 1))
    leaves.sort(key=lambda t_: t_[0])
    return leaves, nev


def tiles(leaves, lo_, hi_):
    """True if the leaves tile [lo_, hi_] exactly."""
    if not leaves or leaves[0][0] != mp.mpf(lo_) or leaves[-1][1] != mp.mpf(hi_):
        return False
    return all(leaves[i][1] == leaves[i + 1][0] for i in range(len(leaves) - 1))


def newton_shrink(fun, a, b, tolw=mp.mpf(10) ** -20, maxit=14):
    """Interval Newton on a cell with one root: the root stays in the result. Raises if the derivative contains 0."""
    for _ in range(maxit):
        if b - a < tolw:
            break
        m = (a + b) / 2
        gm = fun(ivpt(m)).v
        dX = fun(ivr(a, b)).d
        if not exc0(dX):
            raise ValueError("derivative contains 0 in the Newton cell")
        N = ivpt(m) - gm / dX
        na, nb = max(a, mp.mpf(N.a)), min(b, mp.mpf(N.b))
        if na > nb:
            raise ValueError("Newton cell has no root")
        a, b = na, nb
    return a, b


# ----------------------------------------------------------------------------------------------------------------------
# E1, E2 on boxes, Krawczyk, bisection
# ----------------------------------------------------------------------------------------------------------------------
class Prob:
    """What the classification needs: E1, E2 term tables per power of w, jets(X) -> (phi, A', B') as Du, T(X) as Du."""

    def __init__(self, name, E, jets, Tfun):
        self.name, self.E, self.jets, self.Tfun = name, E, jets, Tfun
        self.imax = max(t[0] for kk in E for co in kk for t in co)
        self.jmax = max(t[1] for kk in E for co in kk for t in co)
        self.kmax = max(t[2] for kk in E for co in kk for t in co)
        self.cf = [[[(i, j, k, cfi(p, q)) for i, j, k, p, q in co] for co in kk] for kk in E]


def real_prob(wl):
    def jets(X):
        J = jets_real(X, wl)
        return J["phi"], J["A"], J["B"]

    def Tfun(X):
        J = jets_real(X, wl)
        return J["A"] + J["B"]
    return Prob("tanh " + wl.name, tables()["E"], jets, Tfun)


def E_all(pr, X, W):
    """E1, E2 on z in X, w in W: values, z-derivatives and w-derivatives as intervals."""
    phi, A, B = pr.jets(X)
    pp, pa, pb = powers(phi, pr.imax), powers(A, pr.jmax), powers(B, pr.kmax)
    Wp = [iv.mpf(1)]
    for _ in range(4):
        Wp.append(Wp[-1] * W)
    vals, dz, dw = [], [], []
    for kk in range(2):
        Ev, Ez, Ew = iv.mpf(0), iv.mpf(0), iv.mpf(0)
        for m, terms in enumerate(pr.cf[kk]):
            c = du_const(0)
            for i, j, k, cf in terms:
                c = c + pp[i] * pa[j] * pb[k] * cf
            Ev = Ev + c.v * Wp[m]
            Ez = Ez + c.d * Wp[m]
            if m:
                Ew = Ew + c.v * (m * Wp[m - 1])
        vals.append(Ev)
        dz.append(Ez)
        dw.append(Ew)
    return vals, dz, dw


def krawczyk(pr, Z, W):
    """One Krawczyk step on the box Z x W (intervals). Returns (K_z, K_w, inside)."""
    zm, wm = ivpt(mid(Z)), ivpt(mid(W))
    Fm, Jz, Jw = E_all(pr, zm, wm)
    Fx, Jzx, Jwx = E_all(pr, Z, W)
    # E_all at a point interval returns the derivative at that point too, so (Jz, Jw) is the midpoint Jacobian
    Jmid = mp.matrix([[mid(Jz[0]), mid(Jw[0])], [mid(Jz[1]), mid(Jw[1])]])
    Y = Jmid ** -1
    Yi = [[iv.mpf(Y[a_, b_]) for b_ in range(2)] for a_ in range(2)]
    JX = [[Jzx[0], Jwx[0]], [Jzx[1], Jwx[1]]]
    dX = [Z - zm, W - wm]
    xm = [zm, wm]
    K = []
    for a_ in range(2):
        base = xm[a_] - (Yi[a_][0] * Fm[0] + Yi[a_][1] * Fm[1])
        if MUT != "krawczyk-nocorr":
            for b_ in range(2):
                co = (1 if a_ == b_ else 0) - (Yi[a_][0] * JX[0][b_] + Yi[a_][1] * JX[1][b_])
                base = base + co * dX[b_]
        K.append(base)
    inside = K[0].a > Z.a and K[0].b < Z.b and K[1].a > W.a and K[1].b < W.b
    return K, inside


def prune(pr, za, zb, wa, wb, maxdepth=200, maxboxes=30000):
    """All common roots of E1, E2 in the box [za, zb] x [wa, wb]: boxes with a unique root (Krawczyk), and boxes left
    unresolved. Bisects the wider side."""
    stack = [(mp.mpf(za), mp.mpf(zb), mp.mpf(wa), mp.mpf(wb), 0)]
    roots, unres, nbox = [], [], 0
    while stack:
        a, b, c_, d_, dep = stack.pop()
        nbox += 1
        if nbox > maxboxes:
            unres.append((a, b, c_, d_))
            continue
        Z, W = ivr(a, b), ivr(c_, d_)
        vals, _, _ = E_all(pr, Z, W)
        if exc0(vals[0]) or exc0(vals[1]):
            continue
        try:
            K, inside = krawczyk(pr, Z, W)
        except ZeroDivisionError:
            K, inside = None, False
        if inside:
            roots.append((a, b, c_, d_))
            continue
        if dep >= maxdepth:
            unres.append((a, b, c_, d_))
            continue
        if (b - a) >= (d_ - c_):
            m = (a + b) / 2
            stack += [(a, m, c_, d_, dep + 1), (m, b, c_, d_, dep + 1)]
        else:
            m = (c_ + d_) / 2
            stack += [(a, b, c_, m, dep + 1), (a, b, m, d_, dep + 1)]
    return roots, unres


def refine_box(pr, a, b, c_, d_, n=6):
    """Krawczyk iteration X <- K(X) cap X on a box that already has K inside, for a tight enclosure."""
    for _ in range(n):
        K, inside = krawczyk(pr, ivr(a, b), ivr(c_, d_))
        if not inside:
            break
        a, b, c_, d_ = mp.mpf(K[0].a), mp.mpf(K[0].b), mp.mpf(K[1].a), mp.mpf(K[1].b)
    return a, b, c_, d_


def boxes_disjoint(boxes):
    for x in range(len(boxes)):
        for y in range(x + 1, len(boxes)):
            a, b, c_, d_ = boxes[x]
            e, f, g_, h = boxes[y]
            if not (b < e or f < a or d_ < g_ or h < c_):
                return False
    return True


def classify_cell(pr, za, zb, wlim=(0, 1), diag=W_DIAG):
    """Common roots of E1, E2 on a z-cell [za, zb] (w in [-diag, diag]). Returns dict(points=[(box, T, admissible)],
    unresolved=[...], disjoint=bool). A point is admissible when 0 < w < 1 (or the mutated range) with the whole box
    inside, and T excludes 0 on the box."""
    roots, unres = prune(pr, za, zb, -diag, diag)
    pts = []
    refined = []
    for bx in roots:
        refined.append(refine_box(pr, *bx))
    for bx in refined:
        a, b, c_, d_ = bx
        Tn = pr.Tfun(ivr(a, b)).v
        lo_w, hi_w = (-W_DIAG, W_DIAG) if MUT == "w-range" else wlim
        inw = c_ > lo_w and d_ < hi_w
        outw = d_ <= lo_w or c_ >= hi_w
        tnz = exc0(Tn) or MUT == "T-off"
        if not (inw or outw):
            unres.append(bx)
            continue
        pts.append(dict(box=bx, T=Tn, inw=inw, tnz=tnz, adm=inw and tnz))
    return dict(points=pts, unresolved=unres, disjoint=boxes_disjoint(refined))


# ----------------------------------------------------------------------------------------------------------------------
# the proofs at the two ends
# ----------------------------------------------------------------------------------------------------------------------
def tail_P_lemma(wl):
    """On [z_t, inf): g >= z^9 L with L = c9_lo - sum_{k<9} |c_{9,k}| z_t^(k-9) - sum_{n>=10} M_n s_t^(n-9) > 0, where c_{n,k}
    are the coefficients of z^k in the part of (A', B') degree n of g/t^0 with (phi, a~, b~) enclosed on eps in [0, eps_t],
    M_n = sum_k |c_{n,k}|, s_t = t(z_t) z_t >= t z for z >= z_t (t z = 2 v (1 + C) z^2 e^-z decreases for z >= 2)."""
    if MUT == "no-tail-large":
        return False, "removed"
    zt = iv.mpf(wl.zt.numerator) / wl.zt.denominator
    e_hi = iv.exp(-zt).b
    E = ivr(0, e_hi)
    C, v = wl.C, wl.v
    D = 1 + 2 * C * E + E * E
    phi = v * (1 - E) * (1 - E) / D
    at = (1 - E * E) / (D * D)
    bt = ((1 + E * E) ** 2 - 2 * C * E * (1 + E * E) - 8 * E * E) / (D * D * D)
    c = {}
    for i, j, k, p, q in P_table():
        n = j + k
        sg = -1 if j % 2 else 1
        term = cfi(sg * p, q) * phi ** i * at ** j * bt ** k
        c[(n, k)] = c.get((n, k), iv.mpf(0)) + term
    c99 = c[(9, 9)]
    if not c99.a > 0:
        return False, "c_{9,9} lower bound %s not positive" % mp.nstr(mp.mpf(c99.a), 6)
    sub = iv.mpf(0)
    for (n, k), val in c.items():
        if n == 9 and k < 9:
            sub = sub + max(abs(val.a), abs(val.b)) / zt ** (9 - k)
    st = 2 * v * (1 + C) * zt * zt * iv.exp(-zt)
    st_hi = iv.mpf(st.b)
    rest = iv.mpf(0)
    for n in sorted({nn for nn, _ in c if nn > 9}):
        Mn = iv.mpf(0)
        for (nn, k), val in c.items():
            if nn == n:
                Mn = Mn + max(abs(val.a), abs(val.b))
        rest = rest + Mn * st_hi ** (n - 9)
    if min(nn for nn, _ in c) < 9:
        return False, "a term of (A', B') degree below 9"
    L = iv.mpf(c99.a) - sub - rest
    ok = L.a > 0
    return ok, "z_t = %s, c_{9,9} >= %s, sum_(k<9) |c| z_t^(k-9) <= %s, higher degrees <= %s, L >= %s" % (
        wl.zt, mp.nstr(mp.mpf(c99.a), 8), mp.nstr(mp.mpf(sub.b), 6), mp.nstr(mp.mpf(rest.b), 6), mp.nstr(mp.mpf(L.a), 8))


def tail_Y_lemma(wl):
    """On [z_t, inf): Yn = (1 - eps^4) + 2 C eps (1 - eps^2) + z [(1 + eps^2)^2 - 2 C eps (1 + eps^2) - 8 eps^2] >= z beta > 0
    with beta = 1 - 2 kappa_t (1 + eps_t^2) - 8 eps_t^2, kappa_t = C eps_t (C eps is decreasing for z >= 0)."""
    if MUT == "no-tail-large":
        return False, "removed"
    zt = iv.mpf(wl.zt.numerator) / wl.zt.denominator
    et = iv.exp(-zt)
    et_hi = iv.mpf(et.b)
    kap = wl.C * et_hi
    beta = 1 - 2 * kap * (1 + et_hi * et_hi) - 8 * et_hi * et_hi
    return beta.a > 0, "kappa_t = C e^-z_t <= %s, beta >= %s" % (mp.nstr(mp.mpf(kap.b), 6), mp.nstr(mp.mpf(beta.a), 8))


def _fact(n):
    return math.factorial(n)


def series_cell(wl, a, b, K=16):
    """Enclosures of Phi, Psi, Delta on z in [a, b] subset of [0, 1/2]: partial sums over s = z^2 plus the remainder."""
    C, v = wl.C, wl.v
    A_, B_ = iv.mpf(a), iv.mpf(b)
    s = ivr(mp.mpf(a) ** 2, mp.mpf(b) ** 2)
    zc = ivr(a, b)
    chz = (iv.exp(zc) + iv.exp(-zc)) / 2
    D = C + chz
    # S2 = (cosh z - 1)/z^2 = sum_{k>=1} s^(k-1)/(2k)!
    S2 = iv.mpf(0)
    for k in range(K, 0, -1):
        S2 = S2 * s + iv.mpf(1) / _fact(2 * k)
    S2 = S2 + ivr(0, iv.mpf(1) / _fact(2 * K + 2) * 2 * (iv.mpf(b) ** (2 * K)))
    Phi = v * S2 / D
    # Nn1 = N1/z^4 = sum_{k>=2} n_k s^(k-2), n_k = (4^k + 2C - 2 - 2k(1 + C))/(2k)!
    def nk(k):
        return (iv.mpf(4) ** k + 2 * C - 2 - 2 * k * (1 + C)) / _fact(2 * k)

    def mu(k):
        return (iv.mpf(4) ** k + 2 * C + 2 + 2 * k * (1 + C)) / _fact(2 * k)
    N1 = iv.mpf(0)
    for k in range(K, 1, -1):
        N1 = N1 * s + nk(k)
    rem1 = mu(K + 1) * (iv.mpf(b) ** (2 * K - 2)) * 2
    N1 = N1 + ivr(-rem1.b, rem1.b)
    Psi = v * N1 / (D * D)
    # Nn2 = N2/z^3 = sum_{k>=1} m_k s^(k-1), m_k = (4^k (2k + 3)/2 - 2 k C)/(2k + 1)!
    def mk(k):
        return (iv.mpf(4) ** k * (2 * k + 3) / 2 - 2 * k * C) / _fact(2 * k + 1)

    def nu(k):
        return (iv.mpf(4) ** k * (2 * k + 3) / 2 + 2 * k * C) / _fact(2 * k + 1)
    N2 = iv.mpf(0)
    for k in range(K, 0, -1):
        N2 = N2 * s + mk(k)
    rem2 = nu(K + 1) * (iv.mpf(b) ** (2 * K)) * 2
    N2 = N2 + ivr(-rem2.b, rem2.b)
    Delta = v * (1 + C) * N2 / (D * D * D)
    return Phi, Psi, Delta


_PH = {}


def weighted_P(terms):
    """P(phi, psi - 2 phi, psi - 2 phi + delta) as {(i, j, k): Fraction} in (phi, psi, delta), with weights 2 i + 4 j + 4 k."""
    key = tuple(terms)
    if key in _PH:
        return _PH[key]
    R, phi, psi, dl = ring("phi,psi,delta", QQ)
    Ap_, Bp_ = psi - 2 * phi, psi - 2 * phi + dl
    tot = R.zero
    for i, j, k, p, q in terms:
        tot += QQ(p, q) * phi ** i * Ap_ ** j * Bp_ ** k
    out = {mn: Fr(int(cf.numerator), int(cf.denominator)) for mn, cf in tot.terms()}
    _PH[key] = out
    return out


def small_P_lemma(wl, cells=None):
    """P != 0 on (0, z_s]: G = P / z^24 = sum cf Phi^i Psi^j Delta^k z^(2 i + 4 j + 4 k - 24) has an enclosure excluding 0 on
    each cell of a tiling of [0, z_s], and Delta excludes 0 (so y = z^4 Delta v-factor != 0)."""
    if MUT == "no-tail-small":
        return False, "removed"
    wt = weighted_P(P_table())
    if any(2 * i + 4 * j + 4 * k < 24 for (i, j, k) in wt):
        return False, "a term of weight below 24"
    zs = wl.zs
    n = 8
    edges = [Fr(m, n) * zs for m in range(n + 1)]
    res = []
    for a, b in zip(edges[:-1], edges[1:]):
        af, bf = mp.mpf(a.numerator) / a.denominator, mp.mpf(b.numerator) / b.denominator
        Phi, Psi, Delta = series_cell(wl, af, bf)
        zc = ivr(af, bf)
        mx = max(max(mn) for mn in wt)
        pP, pS, pD = powers_iv(Phi, mx), powers_iv(Psi, mx), powers_iv(Delta, mx)
        tot = iv.mpf(0)
        for (i, j, k), cf in wt.items():
            e = 2 * i + 4 * j + 4 * k - 24
            tot = tot + iv.mpf(cf.numerator) / cf.denominator * pP[i] * pS[j] * pD[k] * zc ** e
        ok = exc0(tot) and exc0(Delta)
        res.append((ok, tot, Delta))
        if not ok:
            return False, "cell [%s, %s]: G in [%s, %s], Delta in [%s, %s]" % (
                a, b, mp.nstr(mp.mpf(tot.a), 6), mp.nstr(mp.mpf(tot.b), 6), mp.nstr(mp.mpf(Delta.a), 6), mp.nstr(mp.mpf(Delta.b), 6))
    lead = res[0][1]
    return True, "%d cells tile [0, %s]; G in [%s, %s] on the first cell, G sign %s on all; Delta < 0 on all" % (
        n, zs, mp.nstr(mp.mpf(lead.a), 6), mp.nstr(mp.mpf(lead.b), 6), "+" if all(r[1].a > 0 for r in res) else "-")


def powers_iv(x, n):
    out = [iv.mpf(1)]
    for _ in range(n):
        out.append(out[-1] * x)
    return out


# ----------------------------------------------------------------------------------------------------------------------
# the cover at a wall, the spurious factor, the classification (cached per run)
# ----------------------------------------------------------------------------------------------------------------------
_COV, _SPU, _PTS = {}, {}, {}


def zlim(wl):
    return mp.mpf(wl.zs.numerator) / wl.zs.denominator, mp.mpf(wl.zt.numerator) / wl.zt.denominator


def cover_wall(wl):
    if wl.name in _COV:
        return _COV[wl.name]
    t = time.time()
    lo_, hi_ = zlim(wl)
    leaves, nev = cover(g_fun(wl), lo_, hi_)
    if MUT == "cell-gap":
        xs = [n for n, lf in enumerate(leaves) if lf[2] == "x"]
        del leaves[xs[len(xs) // 2]]
    _COV[wl.name] = dict(leaves=leaves, nev=nev, lo=lo_, hi=hi_, time=time.time() - t)
    return _COV[wl.name]


def spur_wall(wl):
    if wl.name in _SPU:
        return _SPU[wl.name]
    t = time.time()
    lo_, hi_ = zlim(wl)
    leaves, nev = cover(yn_fun(wl), lo_, hi_)
    _SPU[wl.name] = dict(leaves=leaves, nev=nev, lo=lo_, hi=hi_, time=time.time() - t)
    return _SPU[wl.name]


def tanh_values(wl, z, dps=150):
    """(phi, A', B') from the tanh form of he_type3.jets_mp at z = 2 sigma r, in mp at `dps` digits."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        r = mp.mpf(z) / (2 * wl.mps)
        F_, a_, b_, _ = H3.jets_mp(r, wl.mps, mp.mpf(wl.Rb.numerator) / wl.Rb.denominator)
        v = mp.mpf(wl.v0.numerator) / wl.v0.denominator
        return v * F_, v * r * a_, v * r * r * b_
    finally:
        mp.mp.dps = old


def P_tanh(wl, z):
    old = mp.mp.dps
    mp.mp.dps = 150
    try:
        ph_, A_, B_ = tanh_values(wl, z)
        return sum(mp.mpf(p) / q * ph_ ** i * A_ ** j * B_ ** k for i, j, k, p, q in P_table())
    finally:
        mp.mp.dps = old


def y_tanh(wl, z):
    ph_, A_, B_ = tanh_values(wl, z)
    return B_ - A_


def sign_of(I):
    return 1 if I.a > 0 else (-1 if I.b < 0 else 0)


def sign_agreement(leaves, fun, tanh_fun):
    """The signs of the tanh-form function and of the interval function agree at every leaf boundary, the tanh-form sign is
    constant inside each excluded leaf (checked at its midpoint) and changes inside each root leaf; returns (ok, nchanges)."""
    bounds = [leaves[0][0]] + [lf[1] for lf in leaves]
    sg_i = [sign_of(fun(ivpt(z)).v) for z in bounds]
    sg_t = []
    for z in bounds:
        val = tanh_fun(z)
        sg_t.append(1 if val > 0 else -1)
    ok = sg_i == sg_t and 0 not in sg_i
    kinds = [lf[2] for lf in leaves]
    ok = ok and all((sg_t[n] != sg_t[n + 1]) == (kinds[n] == "r") for n in range(len(leaves)))
    for n, lf in enumerate(leaves):
        if lf[2] == "x":
            val = tanh_fun((lf[0] + lf[1]) / 2)
            ok = ok and (1 if val > 0 else -1) == sg_t[n]
    return ok, sum(1 for n in range(len(leaves)) if sg_t[n] != sg_t[n + 1])


def fmt(x, n=12):
    return mp.nstr(mp.mpf(x), n)


def root_cells(leaves, fun):
    return [newton_shrink(fun, a, b) for a, b, kd in leaves if kd == "r"]


def points_wall(wl):
    if wl.name in _PTS:
        return _PTS[wl.name]
    cov = cover_wall(wl)
    g = g_fun(wl)
    pr = real_prob(wl)
    out = []
    for a, b in root_cells(cov["leaves"], g):
        out.append(dict(cell=(a, b), res=classify_cell(pr, a, b)))
    _PTS[wl.name] = out
    return out


def admissible(wl):
    """Admissible points of the wall, sorted by radius: (r interval, w interval, T interval)."""
    pts = []
    for cell in points_wall(wl):
        for p_ in cell["res"]["points"]:
            if p_["adm"]:
                a, b, c_, d_ = p_["box"]
                rr = wl.r_iv(ivr(a, b))
                pts.append((rr, ivr(c_, d_), p_["T"], T_phys(wl, a, b, rr)))
    pts.sort(key=lambda t_: mp.mpf(t_[0].a))
    return pts


# ----------------------------------------------------------------------------------------------------------------------
# the planted problem
# ----------------------------------------------------------------------------------------------------------------------
_PL = {}


def planted_data():
    """E1 = (w - al)(w - al - 1), E2 = prod_j (w - al + kap_j) on the path (phi, A', B') = (z, z^2, z^3). The common roots of
    E1 and E2 are at kap_j = 0 (w = al) and kap_j = -1 (w = al + 1), so P = prod kap_j (1 + kap_j). The exact roots, w and T are
    from sympy, independent of the interval code."""
    key = MUT == "double-root"
    if key in _PL:
        return _PL[key]
    z = sp.symbols("z")
    phs, As, Bs, ws = H3.ph, H3.Ap, H3.Bp, H3.w
    R = sp.Rational

    def aff(poly):
        c = sp.Poly(poly, z).all_coeffs()[::-1] + [0] * 4
        return c[0] + c[1] * phs + c[2] * As + c[3] * Bs

    sep = 0 if key else R(1, 10 ** 5)
    kz = [2 * (z - R(1, 3)) * (z - R(1, 3) - sep) * (z - R(3, 2)), (z - R(5, 7)) * (z - R(11, 5)) * (z - R(13, 3)), z ** 3 + z + 1]
    alz = R(5, 6) - z / 3
    al = aff(alz)
    kap = [aff(k) for k in kz]
    E1 = sp.expand((ws - al) * (ws - al - 1))
    E2 = sp.expand(sp.prod([ws - al + k for k in kap]))
    Pfull = sp.Poly(sp.expand(sp.prod([k * (1 + k) for k in kz])), z)
    Pz = sp.Poly(sp.cancel(Pfull.as_expr() / (z - R(11, 5))), z)  # the factor (z - 11/5) is carried by the spurious factor
    Yz = sp.expand((z - R(11, 5)) * (z ** 2 + 1))
    Tz = sp.expand((z - R(5, 7)) * (z + 7))

    def table(expr):
        Pw = sp.Poly(expr, ws)
        co = Pw.all_coeffs()[::-1]
        return [[(i, j, k, int(cf.p), int(cf.q)) for (i, j, k), cf in sp.Poly(c, phs, As, Bs).terms()] if c != 0 else [] for c in co]
    Ptab = [(n_, 0, 0, int(c_.p), int(c_.q)) for (n_,), c_ in sp.Poly(Pz, z).terms()]  # P_pl(z) written in phi = z
    roots = [r for r in sp.real_roots(Pfull) if r > 0]
    truth = []
    def polyv(expr, x):
        return mp.polyval([mp.mpf(int(c_.p)) / int(c_.q) for c_ in sp.Poly(expr, z).all_coeffs()], x)
    for rt in roots:
        mp.mp.dps = 60
        rv = mp.mpf(sp.N(rt, 60))
        kv = [polyv(k, rv) for k in kz]
        alv = polyv(alz, rv)
        wv = None
        for kk_ in kv:
            if abs(kk_) < mp.mpf(10) ** -25:
                wv = alv
            if abs(1 + kk_) < mp.mpf(10) ** -25:
                wv = alv + 1
        Tv = polyv(Tz, rv)
        truth.append(dict(z=rv, w=wv, T0=abs(Tv) < mp.mpf(10) ** -25, adm=(wv is not None and 0 < wv < 1 and abs(Tv) > mp.mpf(10) ** -25)))
    cs = Pz.all_coeffs()
    n_ = len(cs) - 1

    def dominant(N):
        return sum(abs(c_) * sp.Integer(N) ** (n_ - m_) for m_, c_ in enumerate(cs[1:], start=1)) < abs(cs[0]) * sp.Integer(N) ** n_
    ZC = 1
    while not dominant(ZC):
        ZC *= 2
    lo_n, hi_n = ZC // 2 + 1 if ZC > 1 else 1, ZC
    while lo_n < hi_n:
        m_ = (lo_n + hi_n) // 2
        if dominant(m_):
            hi_n = m_
        else:
            lo_n = m_ + 1
    ZC = lo_n
    tail_ok = MUT != "no-tail-large"
    _PL[key] = dict(E=[table(E1), table(E2)], P=Ptab, Tz=Tz, truth=truth, ZC=int(ZC), tail_ok=tail_ok, nroots=len(roots), deg=Pfull.degree(), Yz=Yz, z=z)
    return _PL[key]


def planted_prob(pd):
    def jets(X):
        z_ = Du(X, iv.mpf(1))
        return z_, z_ * z_, z_ * z_ * z_

    def Tfun(X):
        z_ = Du(X, iv.mpf(1))
        return (z_ - du_const(iv.mpf(5) / 7)) * (z_ + du_const(7))

    def g(X):
        ph_, A_, B_ = jets(X)
        pp, pa, pb = powers(ph_, 20), powers(A_, 20), powers(B_, 20)
        tot = du_const(0)
        for i, j, k, p, q in pd["P"]:
            tot = tot + pp[i] * pa[j] * pb[k] * cfi(p, q)
        return tot

    def yf(X):
        z_ = Du(X, iv.mpf(1))
        return (z_ - du_const(iv.mpf(11) / 5)) * (z_ * z_ + du_const(1))
    return Prob("planted", pd["E"], jets, Tfun), g, yf


# ----------------------------------------------------------------------------------------------------------------------
# checks
# ----------------------------------------------------------------------------------------------------------------------
STATUS = {}


class Harness:
    def __init__(self):
        self.failed = []

    def check(self, name, ok, msg):
        STATUS[name] = bool(ok)
        say("[%s] %s: %s" % ("PASS" if ok else "FAIL", name, msg))
        if not ok:
            self.failed.append(name)


def check_identity(h, walls):
    ok, rows = True, []
    mp.mp.dps = 100
    for name, wl in walls.items():
        sig, Rb = wl.mps, mp.mpf(wl.Rb.numerator) / wl.Rb.denominator
        v = mp.mpf(wl.v0.numerator) / wl.v0.denominator
        zs_ = [mp.mpf("0.7"), mp.mpf("3.1"), 2 * sig * Rb * mp.mpf("0.99"), 2 * sig * Rb + 3]
        worst, wder = mp.mpf(0), mp.mpf(0)
        for z0 in zs_:
            def tanh_form(zz):
                r = zz / (2 * sig)
                F_, a_, b_, _ = H3.jets_mp(r, sig, Rb)
                return v * F_, v * r * a_, v * r * r * b_
            tv = tanh_form(z0)
            J = jets_real(ivpt(z0), wl)
            for n, key in enumerate(("phi", "A", "B")):
                I = J[key].v
                err = abs(mid(I) - tv[n]) / (1 + abs(tv[n]))
                ok = ok and width(I) < mp.mpf(10) ** -30 * (1 + abs(tv[n])) and err < mp.mpf(10) ** -40
                worst = max(worst, err)
                dd = mp.diff(lambda zz: tanh_form(zz)[n], z0)
                derr = abs(mid(J[key].d) - dd) / (1 + abs(dd))
                ok = ok and derr < mp.mpf(10) ** -25
                wder = max(wder, derr)
        # the series forms of Phi, Psi, Delta against psi = A' + 2 phi and delta = B' - A'
        wser = mp.mpf(0)
        for z0 in (mp.mpf("0.1"), mp.mpf("0.2")):
            ph_, A_, B_ = tanh_values(wl, z0, 100)
            Phi, Psi, Delta = series_cell(wl, z0, z0)
            for I, val, e in ((Phi, ph_, 2), (Psi, A_ + 2 * ph_, 4), (Delta, B_ - A_, 4)):
                err = abs(mid(I) * z0 ** e - val) / abs(val)
                ok = ok and err < mp.mpf(10) ** -40
                wser = max(wser, err)
        rows.append("%s: jets within %s of the tanh form, derivatives within %s of mp.diff, series within %s" % (
            name, mp.nstr(worst, 3), mp.nstr(wder, 3), mp.nstr(wser, 3)))
    mp.mp.dps = 60
    h.check("identity", ok, "; ".join(rows))


def sylvester_det(f, g):
    """Sylvester determinant of f (degree m) and g (degree n), coefficients from the top, over Fraction; rows of f first."""
    m, n = len(f) - 1, len(g) - 1
    N = m + n
    M = [[Fr(0)] * N for _ in range(N)]
    for r_ in range(n):
        for c_, x in enumerate(f):
            M[r_][r_ + c_] = x
    for r_ in range(m):
        for c_, x in enumerate(g):
            M[n + r_][r_ + c_] = x
    det = Fr(1)
    for c_ in range(N):
        piv = next((r_ for r_ in range(c_, N) if M[r_][c_] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != c_:
            M[c_], M[piv] = M[piv], M[c_]
            det = -det
        det *= M[c_][c_]
        for r_ in range(c_ + 1, N):
            fac = M[r_][c_] / M[c_][c_]
            if fac:
                M[r_] = [x - fac * y for x, y in zip(M[r_], M[c_])]
    return det


def poly_at(terms, pt):
    tot = Fr(0)
    for i, j, k, p, q in terms:
        tot += Fr(p, q) * pt[0] ** i * pt[1] ** j * pt[2] ** k
    return tot


def check_tables(h):
    tb = tables()
    sh = tb["sh"]
    ph_, Ap_, Bp_, w_ = H3.ph, H3.Ap, H3.Bp, H3.w
    msgs, ok = [], True
    Pt = P_table()
    # (a) exact: Res = (27/256) phi^3 (A' - B')^3 P
    PPoly = sp.Poly({(i, j, k): sp.Rational(p, q) for i, j, k, p, q in Pt}, ph_, Ap_, Bp_)
    spur = sp.Poly(sp.Rational(27, 256) * ph_ ** 3 * (Ap_ - Bp_) ** 3, ph_, Ap_, Bp_)
    diff = sp.Poly(sh["Res"], ph_, Ap_, Bp_) - spur * PPoly
    a_ok = diff.is_zero
    ok = ok and a_ok
    msgs.append("Res = (27/256) phi^3 (A' - B')^3 P exactly: %s (%d terms of P)" % (a_ok, len(Pt)))
    # (b) an independent Sylvester determinant at rational points
    import random
    rnd = random.Random(133)
    b_ok = True
    for _ in range(5):
        pt = tuple(Fr(rnd.randint(1, 40), rnd.randint(1, 40)) * (1 if rnd.random() < 0.5 else -1) for _ in range(3))
        f = [poly_at(co, pt) for co in tb["E"][0]][::-1]
        g = [poly_at(co, pt) for co in tb["E"][1]][::-1]
        res = sylvester_det(f, g)
        want = Fr(27, 256) * pt[0] ** 3 * (pt[1] - pt[2]) ** 3 * poly_at(Pt, pt)
        b_ok = b_ok and res == want
    # exact on three random lines: a determinant of polynomials in one variable against the product, as polynomials
    s_ = sp.Symbol("s")
    lines_ok = True
    for _ in range(3):
        p0 = [sp.Integer(rnd.randint(-9, 9)) for _ in range(3)]
        dd = [sp.Integer(rnd.randint(-9, 9)) for _ in range(3)]
        vs = [p0[n_] + dd[n_] * s_ for n_ in range(3)]

        def on_line(terms):
            return sp.Poly(sum(sp.Rational(p, q) * vs[0] ** i * vs[1] ** j * vs[2] ** k for i, j, k, p, q in terms), s_)
        fl = [on_line(co) for co in tb["E"][0]][::-1]
        gl = [on_line(co) for co in tb["E"][1]][::-1]
        Mx = sp.zeros(5, 5)
        for r_ in range(3):
            for c_, x in enumerate(fl):
                Mx[r_, r_ + c_] = x.as_expr()
        for r_ in range(2):
            for c_, x in enumerate(gl):
                Mx[3 + r_, r_ + c_] = x.as_expr()
        det_l = sp.Poly(sp.expand(Mx.det(method="berkowitz")), s_)
        want_l = sp.Poly(sp.expand(sp.Rational(27, 256) * vs[0] ** 3 * (vs[1] - vs[2]) ** 3 * on_line(Pt).as_expr()), s_)
        lines_ok = lines_ok and (det_l - want_l).is_zero
    b_ok = b_ok and lines_ok
    ok = ok and b_ok
    msgs.append("Sylvester determinant in fractions equals the same product at 5 rational points, and exactly as polynomials on 3 random lines: %s" % b_ok)
    # (c) E1, E2 from the characteristic polynomial of the frame table printed in he_type3.md section 1
    X_, ga_, y_, q_ = sp.symbols("X gamma y q")
    ph, w = ph_, w_
    R00 = -X_ * (1 + w) / 2 + ph * (ga_ + w * y_)
    R0x = ga_ + (1 - w) * y_ / 2
    R0r = -q_ * y_ / 2
    Rxx = X_ * (3 * w - 1) / 2 - ph * (ga_ + w * y_)
    Rxr = q_ * (X_ - ph * y_ / 2)
    Rrr = X_ * (1 - w) / 2
    Rm = sp.Matrix([[R00, R0x, R0r], [R0x, Rxx, Rxr], [R0r, Rxr, Rrr]])
    M = sp.diag(-1, 1, 1) * Rm
    c2 = -M.trace()
    c1 = sum(M[a, a] * M[b, b] - M[a, b] * M[b, a] for a, b in [(0, 1), (0, 2), (1, 2)])
    c0 = -M.det()
    E1c, E2c = sp.expand(c2 ** 2 - 3 * c1), sp.expand(27 * c0 - c2 ** 3)

    def red(e):
        pq = sp.Poly(e, q_)
        out = 0
        for (m_,), cf in pq.terms():
            if m_ % 2:
                raise ValueError("odd power of c s")
            out += cf * (w * (1 - w)) ** (m_ // 2)
        return sp.expand(out)
    sub = {X_: Ap_ ** 2, ga_: Ap_, y_: Bp_ - Ap_}
    c_ok = True
    for kk, Ec in enumerate((E1c, E2c)):
        tab_e = sum(sp.Rational(p, q) * ph_ ** i * Ap_ ** j * Bp_ ** k * w_ ** m_ for m_, co in enumerate(tb["E"][kk]) for i, j, k, p, q in co)
        c_ok = c_ok and sp.expand(red(Ec).subs(sub) - tab_e) == 0
    ok = ok and c_ok
    msgs.append("E1, E2 equal c2^2 - 3 c1 and 27 c0 - c2^3 of the printed frame table (u = 1): %s" % c_ok)
    # (d) structure
    degs = {j + k for _, j, k, _, _ in Pt}
    tot = {i + j + k for i, j, k, _, _ in Pt}
    d_ok = min(degs) == 9 and max(degs) == 15 and tot == {10, 12, 14, 16, 18}
    ok = ok and d_ok
    msgs.append("(A', B') degrees %d to %d, total degrees %s: %s" % (min(degs), max(degs), sorted(tot), d_ok))
    wt = weighted_P(Pt)
    wmin = min(2 * i + 4 * j + 4 * k for (i, j, k) in wt)
    lead = {mn: cf for mn, cf in wt.items() if 2 * mn[0] + 4 * mn[1] + 4 * mn[2] == wmin}
    e_ok = wmin == 24 and lead == {(8, 0, 2): Fr(774144)}
    ok = ok and e_ok
    msgs.append("weighted expansion in (phi, psi, delta): minimum weight %d, leading part %s: %s" % (wmin, {k_: str(v_) for k_, v_ in lead.items()}, e_ok))
    h.check("tables", ok, "; ".join(msgs))


def check_small_z(h, walls):
    ok, rows = True, []
    for name, wl in walls.items():
        r_, m_ = small_P_lemma(wl)
        ok = ok and r_
        rows.append("%s: %s" % (name, m_))
    h.check("small_z", ok, "; ".join(rows))


def check_tail_z(h, walls):
    ok, rows = True, []
    for name, wl in walls.items():
        r1, m1 = tail_P_lemma(wl)
        r2, m2 = tail_Y_lemma(wl)
        ok = ok and r1 and r2
        rows.append("%s: P %s; y %s" % (name, m1, m2))
    h.check("tail_z", ok, "; ".join(rows))


def check_cover(h, walls):
    ok, rows = True, []
    for name, wl in walls.items():
        cv = cover_wall(wl)
        lv = cv["leaves"]
        t_ok = tiles(lv, cv["lo"], cv["hi"])
        u = sum(1 for x in lv if x[2] == "u")
        nr = sum(1 for x in lv if x[2] == "r")
        s_ok, nch = sign_agreement(lv, g_fun(wl), lambda z: P_tanh(wl, z)) if t_ok and u == 0 else (False, -1)
        good = t_ok and u == 0 and s_ok and nch == nr
        ok = ok and good
        rows.append("%s: %d leaves tile [%s, %s]: %s, %d unresolved, %d root leaves, tanh-form P changes sign %d times at the leaf ends "
                    "(agreement %s), %d evaluations, %.0f s" % (name, len(lv), wl.zs, wl.zt, t_ok, u, nr, nch, s_ok, cv["nev"], cv["time"]))
    h.check("cover", ok, "; ".join(rows))


def check_spurious(h, walls):
    ok, rows = True, []
    for name, wl in walls.items():
        sp_ = spur_wall(wl)
        lv = sp_["leaves"]
        t_ok = tiles(lv, sp_["lo"], sp_["hi"])
        u = sum(1 for x in lv if x[2] == "u")
        nr = sum(1 for x in lv if x[2] == "r")
        s_ok, nch = sign_agreement(lv, yn_fun(wl), lambda z: y_tanh(wl, z)) if t_ok and u == 0 else (False, -1)
        cells = root_cells(lv, yn_fun(wl)) if t_ok and u == 0 else []
        if MUT == "skip-y":
            cells = []
        pr = real_prob(wl)
        gcells = [c_["cell"] for c_ in points_wall(wl)]
        nadm, nunres, disj = 0, 0, True
        desc = []
        for a, b in cells:
            res = classify_cell(pr, a, b)
            nadm += sum(1 for p_ in res["points"] if p_["adm"])
            nunres += len(res["unresolved"])
            disj = disj and res["disjoint"] and all(b < ga or gb < a for ga, gb in gcells)
            desc.append("y = 0 at z in [%s, %s] (r = %s), common roots %d, w in %s" % (
                fmt(a, 10), fmt(b, 10), fmt(wl.r_iv(ivr(a, b)).a, 10), len(res["points"]),
                [fmt(mid(ivr(p_["box"][2], p_["box"][3])), 6) for p_ in res["points"]]))
        good = t_ok and u == 0 and s_ok and nch == nr and len(cells) == nr and nadm == 0 and nunres == 0 and disj
        ok = ok and good
        rows.append("%s: %d leaves tile: %s, %d unresolved, %d root leaves, tanh-form y changes sign %d times (agreement %s); %s; "
                    "admissible points at y-roots %d, unresolved boxes %d, cells disjoint from the P cells %s" % (
                        name, len(lv), t_ok, u, nr, nch, s_ok, "; ".join(desc), nadm, nunres, disj))
    h.check("spurious", ok, " ... ".join(rows))


EXPECT_OUT = {"Le": [58.459], "Rodal": [4.1007, 8.7173, 1.0300]}  # he_type3_output.txt, rejected roots


def T_phys(wl, a_, b_, r_iv):
    """T = b + a u = (A' + B') / (v r^2), as an interval, from the box (z in [a_, b_])."""
    return wl.Tn_of(ivr(a_, b_)) / (wl.v * r_iv * r_iv)


def check_points(h, walls):
    ok, rows = True, []
    for name, wl in walls.items():
        cells = points_wall(wl)
        adm = admissible(wl)
        unres = sum(len(c_["res"]["unresolved"]) for c_ in cells)
        disj = all(c_["res"]["disjoint"] for c_ in cells)
        outs = []
        desc = []
        for c_ in cells:
            za, zb = c_["cell"]
            rr = wl.r_iv(ivr(za, zb))
            items = []
            for p_ in c_["res"]["points"]:
                if not p_["adm"]:
                    outs.append(mid(ivr(p_["box"][2], p_["box"][3])))
                a, b, c2, d2 = p_["box"]
                Tp = T_phys(wl, a, b, wl.r_iv(ivr(a, b)))
                items.append("w = %s %s, T = %s" % (fmt(mid(ivr(c2, d2)), 8), "admissible" if p_["adm"] else (
                    "T = 0" if p_["inw"] and not p_["tnz"] else "outside (0, 1)"), fmt(mid(Tp), 6)))
            desc.append("P-root r = %s (z = %s, width %s): %s" % (fmt(mid(rr), 12), fmt(mid(ivr(za, zb)), 12), mp.nstr(zb - za, 2),
                                                                  "; ".join(items) if items else "no common root, w in [-1e4, 1e4]"))
        exp_out = sorted(EXPECT_OUT[name])
        out_ok = len(outs) == len(exp_out) and all(abs(float(o) - e) < 5e-4 * max(1, e) for o, e in zip(sorted(outs), exp_out))
        tab_ok = len(adm) == 3 and all(abs(float(mid(a_[0])) - t_) < 5e-7 for a_, t_ in zip(adm, wl.table))
        good = unres == 0 and disj and len(adm) == 3 and tab_ok and out_ok and len(cells) == 3 + len(exp_out)
        ok = ok and good
        pts = "; ".join("r = %s, cos(theta) = %s, T = %s" % (
            fmt(mid(a_[0]), 12), fmt(mp.sqrt(mid(a_[1])), 12), fmt(mid(a_[3]), 6)) for a_ in adm)
        rows.append("%s (sigma = %s, R_b = %s, v = %s): %d P-root cells, %d admissible points (table radii within 5e-7: %s), %d with w outside "
                    "[0, 1]: w = %s, unresolved boxes %d, boxes disjoint %s. Cells: %s. Admissible: %s" % (
                        name, float(wl.sig), wl.Rb, wl.v0, len(cells), len(adm), tab_ok, len(outs), [fmt(o, 7) for o in sorted(outs)],
                        unres, disj, " | ".join(desc), pts))
    h.check("points", ok, " ... ".join(rows))


def check_cross(h, walls):
    txt = open(os.path.join(HERE, "he_type3_output.txt")).read()
    line = next(ln for ln in txt.splitlines() if ln.startswith("[PASS] tanh_points"))
    vals = re.findall(r"r = ([0-9.]+) \(f = [0-9.]+, in mask: \w+\), cos\(theta\) = ([0-9.]+)", line)
    ok, rows = True, []
    for wi, (name, wl) in enumerate(walls.items()):
        want = vals[3 * list(make_walls()).index(name): 3 * list(make_walls()).index(name) + 3]
        adm = admissible(wl)
        good = len(adm) == 3
        worst = mp.mpf(0)
        if good:
            for a_, (rv, cv) in zip(adm, want):
                dr = abs(mid(a_[0]) - mp.mpf(rv))
                dc = abs(mp.sqrt(mid(a_[1])) - mp.mpf(cv))
                worst = max(worst, dr, dc)
                good = good and dr < mp.mpf(10) ** -11 and dc < mp.mpf(10) ** -11
        ok = ok and good
        rows.append("%s: radii and cos(theta) of the %d admissible points within %s of he_type3_output.txt (limit 1e-11)" % (
            name, len(adm), mp.nstr(worst, 3)))
    h.check("cross", ok, "; ".join(rows))


def check_assembly(h, walls):
    need = ["small_z", "tail_z", "cover", "spurious"]
    ok = all(STATUS.get(n_) for n_ in need)
    h.check("assembly", ok, "(0, z_s] by small_z, [z_s, z_t] by cover and spurious, [z_t, infinity) by tail_z: %s" % (
        ", ".join("%s %s" % (n_, "PASS" if STATUS.get(n_) else "not passed") for n_ in need)))


def check_planted(h):
    pd = planted_data()
    pr, g, yf = planted_prob(pd)
    ZC = pd["ZC"]
    leaves, nev = cover(g, mp.mpf(0), mp.mpf(ZC))
    t_ok = tiles(leaves, 0, ZC)
    u = sum(1 for x in leaves if x[2] == "u")
    nr = sum(1 for x in leaves if x[2] == "r")
    yleaves, nev2 = cover(yf, mp.mpf(0), mp.mpf(ZC))
    ty_ok = tiles(yleaves, 0, ZC) and not any(x[2] == "u" for x in yleaves)
    nry = sum(1 for x in yleaves if x[2] == "r")
    nev += nev2
    truth = pd["truth"]
    cnt_ok = nr + nry == pd["nroots"] and u == 0 and ty_ok and nry == 1
    pts, unres = [], 0
    if t_ok and u == 0 and ty_ok:
        cells = root_cells(leaves, g) + ([] if MUT == "skip-y" else root_cells(yleaves, yf))
        for a, b in cells:
            res = classify_cell(pr, a, b)
            unres += len(res["unresolved"])
            for p_ in res["points"]:
                if p_["adm"]:
                    pts.append((mid(ivr(p_["box"][0], p_["box"][1])), mid(ivr(p_["box"][2], p_["box"][3]))))
    pts.sort()
    want = sorted((tr["z"], tr["w"]) for tr in truth if tr["adm"])
    match = len(pts) == len(want) and all(abs(a_[0] - b_[0]) < mp.mpf(10) ** -20 and abs(a_[1] - b_[1]) < mp.mpf(10) ** -20
                                          for a_, b_ in zip(pts, want))
    nT0 = sum(1 for tr in truth if tr["T0"])
    nout = sum(1 for tr in truth if tr["w"] is not None and not (0 < tr["w"] < 1))
    ok = t_ok and ty_ok and cnt_ok and match and pd["tail_ok"] and unres == 0
    h.check("planted", ok, "P_pl (the factor z - 11/5 carried by the spurious factor y_pl) and y_pl have %d positive roots by a Sturm sequence (exact); "
            "the cover finds %d root leaves of P_pl and %d of y_pl over [0, %d] (tiles %s %s, %d unresolved, %d evaluations; Cauchy tail proved %s); "
            "admissible points %d found, %d expected (%d roots with w outside (0, 1), %d with T = 0); positions within 1e-20: %s; "
            "unresolved boxes %d" % (
                pd["nroots"], nr, nry, ZC, t_ok, ty_ok, u, nev, pd["tail_ok"],
                len(pts), len(want), nout, nT0, match, unres))


CHECKS = {
    "identity": (check_identity, True),
    "tables": (check_tables, False),
    "planted": (check_planted, False),
    "small_z": (check_small_z, True),
    "tail_z": (check_tail_z, True),
    "cover": (check_cover, True),
    "spurious": (check_spurious, True),
    "points": (check_points, True),
    "cross": (check_cross, True),
    "assembly": (check_assembly, True),
}
DEPS = {"assembly": ["small_z", "tail_z", "cover", "spurious"]}


def run_checks(names, walls):
    h = Harness()
    for nm in names:
        t = time.time()
        try:
            fn, needs = CHECKS[nm]
            fn(h, walls) if needs else fn(h)
        except Exception as ex:
            import traceback
            traceback.print_exc()
            h.check(nm, False, "raised %r" % (ex,))
        say("    (%s: %.0f s)" % (nm, time.time() - t))
    return h


def set_mut(name):
    global MUT
    MUT = name


def main():
    argv = sys.argv[1:]
    import platform
    say("he_type3_count.py: python %s, sympy %s, mpmath %s" % (platform.python_version(), sp.__version__, mp.__version__))
    wsel = next((x.split("=", 1)[1] for x in argv if x.startswith("--walls=")), "")
    csel = next((x.split("=", 1)[1] for x in argv if x.startswith("--checks=")), "")
    mut = next((x.split("=", 1)[1] for x in argv if x.startswith("--mutate=")), "")
    names = csel.split(",") if csel else list(CHECKS)
    if mut:
        if mut not in CONTROLS:
            raise SystemExit("unknown control %s" % mut)
        set_mut(mut)
        targets, desc = CONTROLS[mut]
        if "--only-targets" in argv:
            names = [n_ for n_ in CHECKS if n_ in targets or n_ in sum((DEPS.get(t_, []) for t_ in targets), [])]
        say("control %s: %s; targets %s" % (mut, desc, targets))
    walls = make_walls(wsel.split(",") if wsel else None)
    h = run_checks(names, walls)
    dump = next((x.split("=", 1)[1] for x in argv if x.startswith("--dump-cells=")), "")
    if dump:
        with open(dump, "w") as fh:
            fh.write("# function\twall\tleft end\tright end\tkind (x: no root, r: one root); ends are exact binary fractions\n")
            for fname, cache in (("g", _COV), ("Yn", _SPU)):
                for wn, cv in cache.items():
                    for a, b, kd in cv["leaves"]:
                        fh.write("%s\t%s\t%s\t%s\t%s\n" % (fname, wn, mp.nstr(a, 70), mp.nstr(b, 70), kd))
    if mut:
        fired = [t_ for t_ in targets if t_ in h.failed]
        exact = sorted(h.failed) == sorted(targets)
        say("control %s: failed %s, targets %s -> %s" % (mut, h.failed, targets,
            "FIRED on exactly its targets" if exact else ("FIRED, but on %s" % h.failed if fired else "SILENT")))
        say("wall %.1f s" % (time.time() - T0))
        raise SystemExit(1 if h.failed else 0)
    base = list(h.failed)
    if "--no-controls" in argv or csel:
        say("checks failed: %s; wall %.1f s" % (base or "none", time.time() - T0))
        raise SystemExit(1 if base else 0)
    silent = []
    for cn, (targets, desc) in CONTROLS.items():
        t = time.time()
        out = subprocess.run([sys.executable, os.path.abspath(__file__), "--mutate=" + cn, "--only-targets"], capture_output=True, text=True)
        last = [ln for ln in out.stdout.splitlines() if ln.startswith("control %s: failed" % cn)]
        failed = eval(last[0].split("failed ", 1)[1].split(", targets")[0]) if last else None
        fired = failed is not None and all(t_ in failed for t_ in targets) and out.returncode == 1
        say("control %-16s targets %-60s failed %-60s %s (%.0f s)" % (cn, targets, failed, "FIRED" if fired else "SILENT", time.time() - t))
        if not fired:
            silent.append(cn)
    say("checks failed: %s; controls silent: %s; wall %.1f s" % (base or "none", silent or "none", time.time() - T0))
    raise SystemExit(1 if base or silent else 0)


if __name__ == "__main__":
    main()

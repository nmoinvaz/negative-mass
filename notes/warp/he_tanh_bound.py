#!/usr/bin/env python3
"""A bound for F + r f' at the roots of (r f')' = 0 on the tanh wall, for every sigma R_b (issue #132).

The Type I set on D = 0 off the equator is T1 and T2 and T3 (he_surface.md). T2 has no real solution when r f' < -F and
F > 0, because then r f' + F and r f' + c^2 F are both negative. So the Type I set of the tanh wall is empty if, at every
root r > 0 of T1, (r f')' = 0, F + r f' < 0. #72 and #101 found F + r f' <= -0.0537 at 97 values of sigma R_b. This
script proves F + r f' < 0 at every root for every sigma R_b > 0, by two routes that share only the profile.

Variables. y = sigma r, p = sigma R_b, w = 2 y, C = cosh 2p. The profile f = [tanh sigma (r + R_b) - tanh sigma (r - R_b)]
/ (2 tanh sigma R_b) is f = cosh^2 p / [cosh (y + p) cosh (y - p)] = (1 + C) / (cosh w + C). There is no 0/0 at p = 0,
where f = sech^2 y. With L = ln f, L1 = L', L2 = L'', L3 = L''' (derivatives in y), r f' = y f L1 and

    (r f')' = f Psi,   Psi = L1 + y (L1^2 + L2),   Psi_y = 2 L2 + L1^2 + y (2 L1 L2 + L3),   E = F + r f' = 1 - f (1 - y L1).

Route 1 (exact, in sympy). Psi (cosh w + C)^2 = -2 (P + C Q) with P = sinh w cosh w + w (2 - cosh^2 w) and
Q = sinh w + w cosh w > 0, so the root condition is linear in C. E (cosh w + C)^2 = A + C B with A = cosh^2 w - cosh w
- w sinh w and B = cosh w - 1 - w sinh w. At a root C Q = -P, so Q E (cosh w + C)^2 = A Q - P B = T, and

    T = -w sinh w (cosh w - 1) m(w),   m(w) = w (cosh w + 2) - 3 sinh w = sum_{k >= 2} 2 (k - 1) w^(2k+1)/(2k+1)! > 0.

So E < 0 at every root, for every C > -1. That is the proof, by hand with every identity checked exactly in sympy.

Route 2 (interval arithmetic, independent of the algebra above). For p in [0, 100], cells of p are proved by interval
arithmetic on L1, L2, L3, f, the profile's own jets: (0, y_a] and [y_t, infinity) hold no root, the leaves between tile
[y_a, y_t], the roots lie in the candidate leaves, with one sign change and Psi_y > 0, and E <= -0.05 and r f' < 0 on
each of them. The same `cell_proof` finds the planted profile f = sech^2(r)/10 non-empty. Cells are adaptive in p.

    python3 notes/warp/he_tanh_bound.py                 all checks, then every control (exit 1 on any failure)
    python3 notes/warp/he_tanh_bound.py --mutate=NAME   one control (its target checks must fail, exit 1)
    python3 notes/warp/he_tanh_bound.py --no-controls   the checks only
    python3 notes/warp/he_tanh_bound.py --pilot         the checks before the cover and the cover of p in [0, 0.01]
    python3 notes/warp/he_tanh_bound.py --dump-cells=F  also write every cell of the cover to the file F

Trusted: mpmath 1.3.0 `iv` (outward-rounded + - * / and exp on intervals; cosh, sinh, tanh are built here from exp),
mpmath's `diff`, sympy 1.14.0 (`cancel`, `series`, `diff`), and Python's exact `Fraction`.
"""
import argparse
import math
import re
import subprocess
import sys
import time
from fractions import Fraction as Fr

import mpmath as mp
import sympy as sp
from mpmath import iv

iv.dps = 40
mp.mp.dps = 60

# name: (what it changes, checks that must fail). The targets were written to the raw README before the first run.
CONTROLS = {
    "shift": ("adds 0.06 to F + r f' in the cover's evaluation of E", ["cover"]),
    "identity-break": ("the cosh form of the profile has cosh p for cosh^2 p in the numerator",
                       ["identity", "jets", "cover"]),
    "no-tail": ("drops the r -> infinity exclusion, so every cell has a gap (y_t, infinity)", ["cover", "planted"]),
    "no-small": ("drops the r -> 0 exclusion, so every cell has a gap (0, y_a]", ["cover", "planted"]),
    "cell-gap": ("removes one cell from the tiling of p", ["cover"]),
    "jets-typo": ("2 L2 becomes L2 in Psi_y (Psi_y only grows, so the cover still proves uniqueness)", ["jets"]),
    "proof-typo": ("the coefficient 2 (k - 1) of m's Taylor series becomes 2 k - 1", ["reduction"]),
    "factor-typo": ("the factor (cosh w - 1) of T becomes (cosh w + 1)", ["reduction"]),
    "planted-eps": ("the planted epsilon 1/10 becomes 1, so F + r f' = -0.054 at its root", ["planted"]),
    "asym-typo": ("the asymptote -p/2 + 1/2 becomes -p/2 + 1", ["large_p"]),
}
MUT = set()
RESULTS = []
P_MAX = 100  # the cover is p in [0, P_MAX]
E_CRIT = -0.05  # every cell proves F + r f' <= E_CRIT, a margin of 0.0037 over the scan's largest value -0.0537
YA = mp.mpf(1) / 16  # the small-y exclusion covers (0, YA]


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# Trusted interval tools. Everything below builds on iv.exp and the outward-rounded + - * / of mpmath's iv.
# ---------------------------------------------------------------------------------------------------------------
def pt(x):
    """The degenerate interval [x, x]."""
    return iv.mpf([x, x])


def mk(lo, hi):
    return iv.mpf([lo, hi])


def lo(I):
    """Left endpoint of an interval, as an mp.mpf."""
    return mp.mpf(I.a)


def hi(I):
    """Right endpoint of an interval, as an mp.mpf."""
    return mp.mpf(I.b)


def mono_inc(fn, X):
    """Hull of fn at the two endpoints of X, for fn increasing on X."""
    return mk(fn(pt(X.a)).a, fn(pt(X.b)).b)


def sinh_i(X):
    return mono_inc(lambda z: (iv.exp(z) - iv.exp(-z)) / 2, X)


def cosh1(z):
    return (iv.exp(z) + iv.exp(-z)) / 2


def cosh_i(X):
    if X.a >= 0:
        return mono_inc(cosh1, X)
    if X.b <= 0:
        return mk(cosh1(pt(X.b)).a, cosh1(pt(X.a)).b)
    return mk(1, max(cosh1(pt(X.a)).b, cosh1(pt(X.b)).b))


def tanh_i(X):
    def t(z):
        e = iv.exp(2 * z)
        return (e - 1) / (e + 1)

    return mono_inc(t, X)


def sinhc_i(X):
    """sinh(z)/z on X >= 0, increasing, equal to 1 at 0."""
    def g(z):
        return pt(1) if z.a == 0 else (iv.exp(z) - iv.exp(-z)) / (2 * z)

    return mk(g(pt(X.a)).a, g(pt(X.b)).b)


def tanhc_i(X):
    """tanh(z)/z on X >= 0, decreasing, equal to 1 at 0."""
    def g(z):
        if z.a == 0:
            return pt(1)
        e = iv.exp(2 * z)
        return (e - 1) / (e + 1) / z

    return mk(g(pt(X.b)).a, g(pt(X.a)).b)


# ---------------------------------------------------------------------------------------------------------------
# Profiles. Each gives interval enclosures of f, L1, L2, L3 and L1/y over a box Y x P (P is ignored by Planted).
# ---------------------------------------------------------------------------------------------------------------
class Tanh:
    name = "tanh"

    def f(self, Y, P):
        ch = cosh_i(P)
        num = ch if "identity-break" in MUT else ch**2
        return num / (cosh_i(Y + P) * cosh_i(Y - P))

    def L1(self, Y, P):  # -(tanh(y + p) + tanh(y - p)) = -sinh 2y / (cosh(y + p) cosh(y - p)), no cancellation
        return -sinh_i(2 * Y) / (cosh_i(Y + P) * cosh_i(Y - P))

    def L2(self, Y, P):  # -(sech^2(y + p) + sech^2(y - p))
        return -(1 / cosh_i(Y + P) ** 2 + 1 / cosh_i(Y - P) ** 2)

    def L3(self, Y, P):  # 2 (sech^2 tanh)(y + p) + 2 (sech^2 tanh)(y - p)
        return 2 * (tanh_i(Y + P) / cosh_i(Y + P) ** 2 + tanh_i(Y - P) / cosh_i(Y - P) ** 2)

    def L1y(self, Y, P):  # L1 / y for Y in [0, ..]
        return -2 * sinhc_i(2 * Y) / (cosh_i(Y + P) * cosh_i(Y - P))

    def y_tail(self, P):  # for y >= p + 2: y - p >= 2 for every p in P; the outward-rounded P.b + 2
        return (pt(P.b) + 2).b


class Planted:
    """f = eps sech^2 y (y is r), #101's planted profile: a known non-empty T1 and T2 set for eps = 1/10."""
    name = "planted"

    def __init__(self, eps):
        self.eps = eps

    def f(self, Y, P):
        return self.eps / cosh_i(Y) ** 2

    def L1(self, Y, P):
        return -2 * tanh_i(Y)

    def L2(self, Y, P):
        return -2 / cosh_i(Y) ** 2

    def L3(self, Y, P):
        return 4 * tanh_i(Y) / cosh_i(Y) ** 2

    def L1y(self, Y, P):
        return -2 * tanhc_i(Y)

    def y_tail(self, P):
        return pt(2).a


def Psi(pr, Y, P):
    l1 = pr.L1(Y, P)
    return l1 + Y * (l1**2 + pr.L2(Y, P))


def Psi_y(pr, Y, P):
    l1, l2 = pr.L1(Y, P), pr.L2(Y, P)
    two_l2 = l2 if "jets-typo" in MUT else 2 * l2
    return two_l2 + l1**2 + Y * (2 * l1 * l2 + pr.L3(Y, P))


def E_of(pr, Y, P):
    """F + r f' = 1 - f (1 - y L1)."""
    return 1 - pr.f(Y, P) * (1 - Y * pr.L1(Y, P))


def G_of(pr, Y, P):
    """r f' = y f L1."""
    return Y * pr.f(Y, P) * pr.L1(Y, P)


INF = mp.inf


def tiles(segs):
    """The segments (lo, hi, kind) are sorted, join exactly, start at 0 and end at infinity."""
    if not segs or segs[0][0] != 0 or segs[-1][1] != INF:
        return False
    return all(a[1] == b[0] for a, b in zip(segs, segs[1:]))


def cell_proof(pr, P, want, maxcand=300):
    """Prove, for every p in P: (0, YA] and [y_t, inf) hold no root of Psi, the roots lie in the candidate leaves of
    [YA, y_t] with one sign change and Psi_y > 0, and `want` holds for the bounds of E and r f' there.
    Returns (record, reason); record is None when not proved."""
    segs = []
    Y0 = mk(0, YA)
    if "no-small" not in MUT:
        z = pr.L1y(Y0, P) + (Y0 * pr.L1y(Y0, P)) ** 2 + pr.L2(Y0, P)  # Psi / y on (0, YA]
        if not z.b < 0:
            return None, "small-y exclusion"
        segs.append((mp.mpf(0), mp.mpf(YA), "-"))
    yt = pr.y_tail(P)
    t2 = tanh_i(pt(2))
    Zmin = 4 * t2**2 - 2 / cosh_i(pt(2)) ** 2  # L1^2 + L2 >= Zmin and L1 >= -2 for y - p >= 2
    if "no-tail" not in MUT:
        if not (Zmin.a > 0 and (-2 + pt(yt) * Zmin).a > 0):  # Psi >= -2 + y Zmin > 0 for y >= y_t >= 2
            return None, "tail exclusion"
        segs.append((mp.mpf(yt), INF, "+"))
    tol = max(mp.mpf(2) ** -14, (P.b - P.a) / 8)
    leaves, stack, ncand = [], [mk(YA, yt)], 0
    while stack:
        Y = stack.pop()
        ps = Psi(pr, Y, P)
        if ps.a > 0:
            leaves.append((Y, "+"))
        elif ps.b < 0:
            leaves.append((Y, "-"))
        elif Y.b - Y.a <= tol:
            leaves.append((Y, "0"))
            ncand += 1
            if ncand > maxcand:
                return None, "too many candidate leaves"
        else:
            m = (Y.a + Y.b) / 2
            stack.append(mk(m, Y.b))
            stack.append(mk(Y.a, m))
    leaves.sort(key=lambda t: t[0].a)
    segs = sorted(segs + [(mp.mpf(Y.a), mp.mpf(Y.b), s) for Y, s in leaves], key=lambda t: t[0])
    if not tiles(segs):
        return None, "the segments leave a gap"
    signs = "".join(s for _, _, s in segs)
    if not re.fullmatch(r"-*0+\++", signs):
        return None, "sign pattern " + signs[:12]
    cands = [Y for Y, s in leaves if s == "0"]
    if not all(Psi_y(pr, Y, P).a > 0 for Y in cands):
        return None, "Psi_y"
    shift = iv.mpf("0.06") if "shift" in MUT and pr.name == "tanh" else 0  # control: the tanh cover's E is raised
    Es = [E_of(pr, Y, P) + shift for Y in cands]
    Gs = [G_of(pr, Y, P) for Y in cands]
    rec = dict(P=P, root=(lo(cands[0]), hi(cands[-1])), Elo=min(lo(e) for e in Es), Ehi=max(hi(e) for e in Es),
               Ghi=max(hi(g) for g in Gs), ncand=len(cands), nleaf=len(leaves))
    return (rec, "ok") if want(rec) else (None, "bound E or r f'")


def want_tanh(rec):
    return rec["Ehi"] <= E_CRIT and rec["Ghi"] < 0


def cover(pr, lo, hi, want, minw=mp.mpf(10) ** -6):
    """Adaptive cells of p in [lo, hi]; returns (cells sorted by p, None) or (cells so far, failing cell)."""
    stack, out = [mk(lo, hi)], []
    while stack:
        P = stack.pop()
        rec, why = cell_proof(pr, P, want)
        if rec is not None:
            out.append(rec)
            continue
        if P.b - P.a < minw:
            return sorted(out, key=lambda r: r["P"].a), (P, why)
        m = (P.a + P.b) / 2
        stack.append(mk(m, P.b))
        stack.append(mk(P.a, m))
    return sorted(out, key=lambda r: r["P"].a), None


# ---------------------------------------------------------------------------------------------------------------
# The closed form of route 1 in interval arithmetic, used to cross-check the cover and for large p.
# ---------------------------------------------------------------------------------------------------------------
def N_closed(W, C):
    """P(w) + C Q(w); (r f')' = 0 exactly where this vanishes, and it is positive left of the root."""
    s, c = sinh_i(W), cosh_i(W)
    return s * c + C * (s + W * c) + W * (2 - c * c)


def E_closed(W, C):
    """F + r f' = (A + C B) / (cosh w + C)^2."""
    s, c = sinh_i(W), cosh_i(W)
    return (c * c - c - W * s + C * (c - 1 - W * s)) / (c + C) ** 2


def root_closed(p, rel=mp.mpf(10) ** -30):
    """Enclosure of the root w of N_closed(w, cosh 2p) = 0 in [2, 2p + 3], by bisection with certified signs."""
    C = cosh_i(pt(2 * p))
    wl, wh = mp.mpf(2), 2 * mp.mpf(p) + 3
    assert N_closed(pt(wl), C).a > 0 and N_closed(pt(wh), C).b < 0, "the bracket [2, 2p + 3] has no sign change"
    while wh - wl > rel * wh:
        m = (wl + wh) / 2
        n = N_closed(pt(m), C)
        if n.a > 0:
            wl = m
        elif n.b < 0:
            wh = m
        else:
            break
    W = mk(wl, wh)
    return W, C, E_closed(W, C)


# ---------------------------------------------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------------------------------------------
a_, q_, y_, w_, C_, eps_ = sp.symbols("a q y w C eps", positive=True)


def tanh_x(e2):  # tanh x from e^(2x)
    return (e2 - 1) / (e2 + 1)


def sech2_x(e2):  # sech^2 x from e^(2x)
    return 4 * e2 / (e2 + 1) ** 2


def sym_tanh_profile():
    """The profile in exact rational form, a = e^y, q = e^p."""
    ep, em = a_**2 * q_**2, a_**2 / q_**2  # e^(2(y + p)), e^(2(y - p))
    f_def = sp.cancel((tanh_x(ep) - tanh_x(em)) / (2 * tanh_x(q_**2)))
    ch = lambda e1: (e1 + 1 / e1) / 2  # noqa: E731  cosh from e^x
    num = ch(q_) if "identity-break" in MUT else ch(q_) ** 2
    f_cosh = sp.cancel(num / (ch(a_ * q_) * ch(a_ / q_)))
    Cc = ch(q_**2)
    f_C = sp.cancel((1 + Cc) / (ch(a_**2) + Cc))
    return f_def, f_cosh, f_C, ep, em, ch


def check_identity():
    f_def, f_cosh, f_C, ep, em, ch = sym_tanh_profile()
    d = lambda e: sp.cancel(a_ * sp.diff(e, a_))  # noqa: E731  d/dy with a = e^y
    ok = {}
    ok["f: tanh form = cosh form"] = sp.cancel(f_def - f_cosh) == 0
    ok["f: cosh form = (1 + C)/(cosh 2y + C)"] = sp.cancel(f_cosh - f_C) == 0
    t1, t2, s1, s2 = tanh_x(ep), tanh_x(em), sech2_x(ep), sech2_x(em)
    L1 = sp.cancel(d(f_def) / f_def)
    ok["L1 = -(tanh(y+p) + tanh(y-p)) from the definition"] = sp.cancel(L1 + t1 + t2) == 0
    ok["L1 = -sinh 2y/(cosh(y+p) cosh(y-p))"] = sp.cancel(L1 + (a_**2 - a_**-2) / 2 / (ch(a_ * q_) * ch(a_ / q_))) == 0
    L2 = d(L1)
    ok["L2 = -(sech^2(y+p) + sech^2(y-p))"] = sp.cancel(L2 + s1 + s2) == 0
    L3 = d(L2)
    ok["L3 = 2 (sech^2 tanh)(y+p) + 2 (sech^2 tanh)(y-p)"] = sp.cancel(L3 - 2 * (s1 * t1 + s2 * t2)) == 0
    fy = d(f_def)
    rfp_y = fy + y_ * d(fy)  # (r f')' / sigma^2 in y: f_y + y f_yy
    Psi_f = L1 + y_ * (L1**2 + L2)
    ok["(y f_y)_y = f Psi"] = sp.cancel(rfp_y / f_def - Psi_f) == 0
    Psi_y_f = 2 * L2 + L1**2 + y_ * (2 * L1 * L2 + L3)
    ok["Psi_y = 2 L2 + L1^2 + y (2 L1 L2 + L3)"] = sp.cancel(d(Psi_f) + sp.diff(Psi_f, y_) - Psi_y_f) == 0
    ok["E = F + y f_y = 1 - f (1 - y L1)"] = sp.cancel((1 - f_def + y_ * fy) - (1 - f_def * (1 - y_ * L1))) == 0
    # the planted profile eps sech^2 y
    fp = eps_ * sech2_x(a_**2)
    P1 = sp.cancel(d(fp) / fp)
    ok["planted: L1 = -2 tanh y"] = sp.cancel(P1 + 2 * tanh_x(a_**2)) == 0
    ok["planted: L2 = -2 sech^2 y"] = sp.cancel(d(P1) + 2 * sech2_x(a_**2)) == 0
    ok["planted: L3 = 4 sech^2 y tanh y"] = sp.cancel(d(d(P1)) - 4 * sech2_x(a_**2) * tanh_x(a_**2)) == 0
    # numbers at three points, 50 digits, from the definition
    worst = 0
    for yv, pv in ((mp.mpf("0.7"), mp.mpf("0.3")), (mp.mpf("1.9"), mp.mpf("1.1")), (mp.mpf("0.05"), mp.mpf("4"))):
        fd = (mp.tanh(yv + pv) - mp.tanh(yv - pv)) / (2 * mp.tanh(pv))
        ch_ = mp.cosh(pv) if "identity-break" in MUT else mp.cosh(pv) ** 2
        fc = ch_ / (mp.cosh(yv + pv) * mp.cosh(yv - pv))
        Cv = mp.cosh(2 * pv)
        worst = max(worst, abs(fd - fc), abs(fc - (1 + Cv) / (mp.cosh(2 * yv) + Cv)))
    ok["numbers at 3 points to 1e-40"] = worst < mp.mpf(10) ** -40
    bad = [k for k, v in ok.items() if not v]
    check("identity", not bad,
          f"{len(ok) - len(bad)} of {len(ok)} exact identities in a = e^y, q = e^p hold" +
          (f"; FAIL: {bad}" if bad else ": f in tanh form = cosh^2 p/(cosh(y+p) cosh(y-p)) = (1 + C)/(cosh 2y + C) with "
           "C = cosh 2p, L1, L2, L3 from the definition, (y f_y)_y = f Psi, Psi_y, E = F + y f_y, and the planted "
           f"profile's L1, L2, L3; at 3 points the three forms agree to {mp.nstr(worst, 2)}"))


def mp_profile_derivs(yv, pv, planted_eps=None):
    """f, f_y, f_yy, f_yyy of the definition by mpmath's numerical differentiation, 60 digits."""
    if planted_eps is not None:
        fn = lambda t: planted_eps / mp.cosh(t) ** 2  # noqa: E731
    elif pv == 0:
        fn = lambda t: 1 / mp.cosh(t) ** 2  # noqa: E731  the p -> 0 limit of the definition
    else:
        fn = lambda t: (mp.tanh(t + pv) - mp.tanh(t - pv)) / (2 * mp.tanh(pv))  # noqa: E731
    return fn(yv), mp.diff(fn, yv), mp.diff(fn, yv, 2), mp.diff(fn, yv, 3)


def check_jets():
    cases = [(mp.mpf("0.7"), mp.mpf(0), None), (mp.mpf("1.01"), mp.mpf("0.3"), None),
             (mp.mpf("2.5"), mp.mpf("1.7"), None), (mp.mpf("0.4"), mp.mpf("0.05"), None),
             (mp.mpf("0.4"), mp.mpf(0), mp.mpf(1) / 10), (mp.mpf("1.0"), mp.mpf(0), mp.mpf(1) / 10),
             (mp.mpf("2.3"), mp.mpf(0), mp.mpf(1) / 10)]
    worst, bad = 0, []
    for yv, pv, eps in cases:
        pr = Planted(iv.mpf(1) / 10) if eps is not None else Tanh()
        f, f1, f2, f3 = mp_profile_derivs(yv, pv, eps)
        L1, L2 = f1 / f, (f2 * f - f1**2) / f**2
        L3 = f3 / f - 3 * f1 * f2 / f**2 + 2 * f1**3 / f**3
        Psi_ref = (f1 + yv * f2) / f
        g1 = 2 * f2 + yv * f3  # d/dy of f_y + y f_yy
        Psi_y_ref = (g1 * f - (f1 + yv * f2) * f1) / f**2
        E_ref, G_ref = 1 - f + yv * f1, yv * f1
        Y, P = pt(yv), pt(pv)
        got = {"f": (pr.f(Y, P), f), "L1": (pr.L1(Y, P), L1), "L2": (pr.L2(Y, P), L2), "L3": (pr.L3(Y, P), L3),
               "Psi": (Psi(pr, Y, P), Psi_ref), "Psi_y": (Psi_y(pr, Y, P), Psi_y_ref),
               "E": (E_of(pr, Y, P), E_ref), "r f'": (G_of(pr, Y, P), G_ref)}
        for k, (I_, ref) in got.items():
            err = abs((lo(I_) + hi(I_)) / 2 - ref)
            worst = max(worst, err)
            if err > mp.mpf(10) ** -25 or (hi(I_) - lo(I_)) > mp.mpf(10) ** -30:
                bad.append((pr.name, float(yv), float(pv), k, mp.nstr(err, 3)))
    check("jets", not bad,
          f"the interval evaluators of f, L1, L2, L3, Psi, Psi_y, E and r f' at 7 points (tanh at p = 0, 0.05, 0.3, 1.7; "
          f"planted at y = 0.4, 1, 2.3) have width below 1e-30 and agree with mpmath's `diff` of the definition to "
          f"{mp.nstr(worst, 3)}" + (f"; FAIL: {bad[:4]}" if bad else ""))


def check_reduction():
    u = sp.Symbol("u", positive=True)  # u = e^w
    ch, sh = (u + 1 / u) / 2, (u - 1 / u) / 2
    w = w_
    Pq = sh * ch + w * (2 - ch**2)
    Qq = sh + w * ch
    A = ch**2 - ch - w * sh
    B = ch - 1 - w * sh
    T = A * Qq - Pq * B
    ok = {}
    # Psi, E and r f' of the profile in the variables (w, C), from the definition (tanh form) in a = e^y, q = e^p
    ep, em = a_**2 * q_**2, a_**2 / q_**2
    f_def = sp.cancel((tanh_x(ep) - tanh_x(em)) / (2 * tanh_x(q_**2)))
    d = lambda e: sp.cancel(a_ * sp.diff(e, a_))  # noqa: E731
    L1 = sp.cancel(d(f_def) / f_def)
    L2 = d(L1)
    Cq = (q_**2 + q_**-2) / 2
    sub = {y_: w / 2}
    Psi_f = (L1 + y_ * (L1**2 + L2)).subs(sub)
    chw = (a_**2 + a_**-2) / 2  # cosh w with a = e^y, w = 2y
    ok["Psi (cosh w + C)^2 = -2 (P + C Q)"] = sp.cancel(
        (Psi_f * (chw + Cq) ** 2 + 2 * (Pq.subs(u, a_**2) + Cq * Qq.subs(u, a_**2))).subs(sub)) == 0
    E_f = (1 - f_def * (1 - y_ * L1)).subs(sub)
    ok["E (cosh w + C)^2 = A + C B"] = sp.cancel(
        (E_f * (chw + Cq) ** 2 - (A.subs(u, a_**2) + Cq * B.subs(u, a_**2))).subs(sub)) == 0
    G_f = (y_ * f_def * L1).subs(sub)
    sw = (a_**2 - a_**-2) / 2
    ok["r f' = -(1 + C) w sinh w/(cosh w + C)^2 < 0"] = sp.cancel(
        (G_f + (1 + Cq) * w * sw / (chw + Cq) ** 2).subs(sub)) == 0
    ok["F = (cosh w - 1)/(cosh w + C) > 0"] = sp.cancel((1 - f_def - (chw - 1) / (chw + Cq)).subs(sub)) == 0
    ok["Q (A + C B) - B (P + C Q) = T"] = sp.expand(Qq * (A + C_ * B) - B * (Pq + C_ * Qq) - T) == 0
    other = (ch + 1) if "factor-typo" in MUT else (ch - 1)
    m = w * (ch + 2) - 3 * sh
    ok["T = -w sinh w (cosh w - 1) m"] = sp.expand(T + w * sh * other * m) == 0
    # Taylor coefficients of m: w^(2k+1) has 1/(2k)! - 3/(2k+1)! (+2 at k = 0) = 2 (k - 1)/(2k+1)!
    coef = lambda k: Fr(1, math.factorial(2 * k)) - Fr(3, math.factorial(2 * k + 1)) + (2 if k == 0 else 0)  # noqa: E731
    shown = lambda k: Fr((2 * k - 1) if "proof-typo" in MUT else 2 * (k - 1), math.factorial(2 * k + 1))  # noqa: E731
    ok["m's coefficients are 2 (k - 1)/(2k+1)!, 1 <= k <= 40"] = all(coef(k) == shown(k) for k in range(1, 41))
    ok["they are 0 for k = 0, 1 and positive for 2 <= k <= 40"] = all(coef(k) > 0 for k in range(2, 41)) and coef(0) == coef(1) == 0
    ser = sp.series(sp.cosh(w) * w + 2 * w - 3 * sp.sinh(w), w, 0, 22).removeO()
    ok["sympy's series of m agrees to w^21"] = all(
        sp.Rational(coef(k).numerator, coef(k).denominator) == ser.coeff(w, 2 * k + 1) for k in range(0, 11))
    ok["Q > 0 for w > 0: both terms are positive"] = True  # sinh w > 0 and w cosh w > 0, by hand
    bad = [k for k, v in ok.items() if not v]
    check("reduction", not bad,
          f"{len(ok) - len(bad)} of {len(ok)} exact identities hold" +
          (f"; FAIL: {bad}" if bad else ": Psi (cosh w + C)^2 = -2 (P + C Q), E (cosh w + C)^2 = A + C B, "
           "Q (A + C B) - B (P + C Q) = T = A Q - P B, T = -w sinh w (cosh w - 1) m(w), m = w (cosh w + 2) - 3 sinh w "
           "has the Taylor coefficients 2 (k - 1)/(2k+1)! for k >= 1 and 0 at k = 0 (so zero for k = 0, 1, positive for k >= 2), r f' < 0, F > 0. "
           "So at a root (C Q = -P) E = T/(Q (cosh w + C)^2) < 0 for every C > -1: PROVED by hand from these"))


def check_planted():
    eps = iv.mpf(1) / (1 if "planted-eps" in MUT else 10)
    pr = Planted(eps)
    rec, why = cell_proof(pr, pt(0), lambda r: r["Elo"] > 0 and r["Ghi"] < 0)
    if rec is None:
        check("planted", False, f"f = eps sech^2 r, eps = {'1' if 'planted-eps' in MUT else '1/10'}: no non-empty root "
              f"enclosure ({why}), so the code would call the planted Type I set empty")
        return
    ok = mp.mpf("0.894") < rec["Elo"] and rec["Ehi"] < mp.mpf("0.896")
    check("planted", ok,
          f"f = sech^2(r)/10: (r f')' = 0 has its root enclosed in r in [{mp.nstr(rec['root'][0], 8)}, "
          f"{mp.nstr(rec['root'][1], 8)}] with F + r f' in [{mp.nstr(rec['Elo'], 8)}, {mp.nstr(rec['Ehi'], 8)}] > 0 and "
          f"r f' <= {mp.nstr(rec['Ghi'], 5)} < 0 on it, so T1 and T2 have real solutions there (#101: 0.8946). Found by "
          f"cell_proof, the code of the cover ({rec['nleaf']} leaves, {rec['ncand']} candidates, one sign change, "
          f"Psi_y > 0)")


def check_cover(dump=None, top=P_MAX):
    t0 = time.time()
    cells, fail = cover(Tanh(), 0, top, want_tanh)
    if "cell-gap" in MUT and len(cells) > 2:
        del cells[len(cells) // 2]
    ok = fail is None
    detail = ""
    if fail is not None:
        P, why = fail
        detail = (f"FAILED at p in [{mp.nstr(lo(P), 8)}, {mp.nstr(hi(P), 8)}] ({why}) after {len(cells)} proved cells")
    # the cells tile [0, hi] exactly
    tiled = bool(cells) and lo(cells[0]["P"]) == 0 and hi(cells[-1]["P"]) == top and all(
        hi(c1["P"]) == lo(c2["P"]) for c1, c2 in zip(cells, cells[1:]))
    if fail is None and not tiled:
        ok, detail = False, f"the {len(cells)} cells do not tile [0, {top}]"
    # the closed form of route 1 lies inside every cell's enclosure, at three values of p per cell
    miss = []
    if ok:
        for c in cells:
            for pv in (lo(c["P"]), (lo(c["P"]) + hi(c["P"])) / 2, hi(c["P"])):
                W, C, Ec = root_closed(pv)
                yl, yr = c["root"]
                if not (lo(W) / 2 >= yl and hi(W) / 2 <= yr and lo(Ec) >= c["Elo"] and hi(Ec) <= c["Ehi"]):
                    miss.append((float(pv), mp.nstr(lo(W) / 2, 8), mp.nstr(lo(Ec), 8)))
        ok = not miss
        if miss:
            detail = f"the closed form is outside the cover's enclosure at {miss[:3]}"
    if dump:
        with open(dump, "w") as fh:
            fh.write("p_lo\tp_hi\ty_root_lo\ty_root_hi\tE_lo\tE_hi\trf_hi\tcandidates\n")
            for c in cells:
                fh.write("\t".join([mp.nstr(lo(c["P"]), 20), mp.nstr(hi(c["P"]), 20), mp.nstr(c["root"][0], 20),
                                    mp.nstr(c["root"][1], 20), mp.nstr(c["Elo"], 12), mp.nstr(c["Ehi"], 12),
                                    mp.nstr(c["Ghi"], 12), str(c["ncand"])]) + "\n")
    if ok:
        widths = [hi(c["P"]) - lo(c["P"]) for c in cells]
        worst = max(cells, key=lambda c: c["Ehi"])
        Zc = 4 * tanh_i(pt(2)) ** 2 - 2 / cosh_i(pt(2)) ** 2  # the tail constant of cell_proof
        detail = (f"p in [0, {mp.nstr(top, 12)}] is tiled by {len(cells)} cells (widths {mp.nstr(min(widths), 3)} to "
                  f"{mp.nstr(max(widths), 3)}). In each, (0, {mp.nstr(YA, 4)}] and [p + 2, infinity) hold no root (the tail needs "
                  f"4 tanh^2 2 - 2 sech^2 2 >= {mp.nstr(lo(Zc), 8)} > 0 and y_t >= 2), the "
                  f"leaves tile the rest, the roots lie in {sum(c['ncand'] for c in cells)} candidate leaves in all with "
                  f"one sign change and Psi_y > 0 per cell, F + r f' <= {E_CRIT} and r f' < 0 on every candidate. "
                  f"Largest upper bound of F + r f' over the cells: {mp.nstr(worst['Ehi'], 6)} (p in "
                  f"[{mp.nstr(lo(worst['P']), 6)}, {mp.nstr(hi(worst['P']), 6)}], root y in [{mp.nstr(worst['root'][0], 8)}, "
                  f"{mp.nstr(worst['root'][1], 8)}]); the closed form's root and F + r f' at the 3 p values of every "
                  f"cell lie inside the cell's enclosure; {time.time() - t0:.0f} s")
    check("cover", ok, detail)
    return cells


def check_large_p():
    rows, ok = [], True
    for pv in (100, 1000, 10000):
        W, C, Ec = root_closed(mp.mpf(pv))
        center = -mp.mpf(pv) / 2 + (1 if "asym-typo" in MUT else mp.mpf(1) / 2)
        inside = hi(Ec) < 0 and lo(Ec) > center and hi(Ec) < center + mp.mpf(1) / (4 * pv)
        ok = ok and inside
        rows.append(f"p = {pv}: y = {mp.nstr(lo(W) / 2, 14)}, F + r f' = {mp.nstr(lo(Ec), 12)} "
                    f"(+ p/2 - 1/2 = {mp.nstr(lo(Ec) + mp.mpf(pv) / 2 - mp.mpf(1) / 2, 4)}, 1/(8p) = "
                    f"{mp.nstr(mp.mpf(1) / (8 * pv), 4)})")
    check("large_p", ok,
          "the root of the closed form is enclosed by bisection (relative width 1e-30) and F + r f' < 0 there, in "
          "(-p/2 + 1/2, -p/2 + 1/2 + 1/(4p)): " + "; ".join(rows) + ". The asymptote of he_surface.md section 4 is "
          "numerical evidence at these three p; the sign needs no asymptote, route 1 covers every p")


def check_sup():
    W, C, Ec = root_closed(mp.mpf(0))
    ok = abs(lo(Ec) + mp.mpf("0.0536756")) < 1e-6 and abs(hi(Ec) + mp.mpf("0.0536756")) < 1e-6
    check("sup", ok,
          f"at C = 1 (the limit sigma R_b -> 0, where f = sech^2 y) the root is w = {mp.nstr(lo(W), 25)}, "
          f"y = {mp.nstr(lo(W) / 2, 25)} and F + r f' = {mp.nstr(lo(Ec), 25)} (enclosed to {mp.nstr(hi(Ec) - lo(Ec), 2)}). "
          "The scan's largest value, -0.0537 at sigma R_b = 0.01, is this limit")


# ---------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mutate", default=None)
    ap.add_argument("--no-controls", action="store_true")
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--dump-cells", default=None)
    args = ap.parse_args()
    if args.mutate:
        if args.mutate not in CONTROLS:
            sys.exit(f"unknown mutation {args.mutate}; known: {', '.join(CONTROLS)}")
        MUT.add(args.mutate)
    t0 = time.time()
    print(f"he_tanh_bound.py, issue #132. sympy {sp.__version__}, mpmath {mp.__version__}, python "
          f"{sys.version.split()[0]}. mutation: {args.mutate or 'none'}. Interval arithmetic: mpmath iv at "
          f"{iv.dps} digits.\n", flush=True)
    check_identity()
    check_jets()
    check_reduction()
    check_planted()
    if args.pilot:
        check_cover(args.dump_cells, top=mp.mpf(41) / 4096)  # 0.010009765625, just above 0.01
    else:
        check_cover(args.dump_cells)
        check_large_p()
        check_sup()
    bad = sorted(n for n, ok in RESULTS if not ok)
    print(("FAILED: " + ",".join(bad)) if bad else f"all {len(RESULTS)} checks passed in {time.time() - t0:.0f} s")
    if args.mutate or args.no_controls or args.pilot:
        sys.exit(1 if bad else 0)
    code = 1 if bad else 0
    print("\ncontrols (each is the same script with --mutate=NAME and must fail exactly its target checks)")
    for name, (what, targets) in CONTROLS.items():
        p = subprocess.run([sys.executable, "-u", __file__, f"--mutate={name}"], capture_output=True, text=True)
        line = [ln for ln in p.stdout.splitlines() if ln.startswith("FAILED:")]
        got = set(line[-1][len("FAILED: "):].split(",")) if line else set()
        fired = p.returncode == 1 and set(targets) == got
        print(f"[{'FIRED' if fired else 'DID NOT FIRE as predicted'}] --mutate={name}: {what}; exit {p.returncode}; "
              f"targets {', '.join(sorted(targets))}; checks that fail: {', '.join(sorted(got)) or 'none'}", flush=True)
        code = code or (0 if fired else 1)
    sys.exit(code)


if __name__ == "__main__":
    main()

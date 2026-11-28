#!/usr/bin/env python3
"""Independent review of the exchange identities E4 (k = 6), E2W and E5 (k = 7) of the large-D matching hierarchy (#264).

Own code. Shares nothing with PanelK / order4.py ... order7.py or the coordinator's ind26 / ind208. The order equations are
expanded from README's master equation at run time (sympy); the hierarchy is a Taylor-jet integration in l with fixed-point
Python integers; the theta integrals are Gauss-Legendre panels in v = sqrt(l) (dtheta = v dv / u, u = tan theta).

Notation (derivations.md of raw I264-rpi5): u = tan th, l = ln R = -2 ell, ell = ln cos th, S = sin^2 th = 1 - exp(-l),
x = pi/2 - th, lt = ln(2 cos th) = ell + ln 2, C = Cl2(pi - 2 th), D3 = 2 int_0^th C = (3/4) zeta(3) + Cl3(pi - 2 th).
  dY_k/dl = F_{k-1} = sum_{j<k} L_{k-1-j} l^j / j!,  S L_k' - L_k/2 = h_k = S L_{k-1} + S_k + Y_k/2,  S_k = sum_{j<k} Y_j c_{k-j}.

Parts: sym | hier | exch | bnd | masters | fp | ctlsum.  Each check prints PASS/FAIL, the script exits 1 if any check fails.
--mutate=NAME changes one thing on purpose (see MUTATIONS) and must make the stated checks fail.
"""
import argparse
import os
import pickle
import sys
import time
from math import sumprod

import mpmath as mp
import sympy as sp

FAILS = []
MUT = None
T0 = time.time()

MUTATIONS = {
    "e4_20": "E4's 22 l^2 as 20 l^2 on the right side",
    "e4_Om1": "E4's -2 Om1 as +2 Om1",
    "e4_C": "E4's 16 u C as 8 u C",
    "e2w_P": "P without the D3 term",
    "e2w_g": "g4's -32 l^2 as -30 l^2, P as derived",
    "e5_noD3": "W without 16 D3 cot th",
    "e5_a1": "a1 = 8 Om1 + 8 Om0^2 + 64 ln2 as 8 Om1 + 8 Om0^2 + 32 ln2 in W only",
    "e5_K": "W with the constant +1/100",
    "hyp_G4": "Om2 shifted by 1e-3 in the hierarchy",
    "c_wrong": "c3's -2 ell^3 as -3 ell^3 in the hierarchy",
    "g6r_slope": "G6 reduced form without the Om4 terms",
    "g7r_slope": "G7 reduced form without the Om5 terms",
    "m1_basis": "M1's closed form with 5 zeta(5) as 4 zeta(5)",
    "m1_J2": "M1's closed form with J2 + 0.01",
    "fp_const": "FP rule's constant term replaced by the next term",
    "m2_K2K2": "<K2^2>'s (61/6) zeta(6) as (60/6) zeta(6)",
    "m2_J5": "J5's 1080 zeta(5) as 1008 zeta(5)",
    "m2_M3": "<lt^3 D3>'s (11/4) zeta(6) as (10/4) zeta(6)",
}


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (" : " + detail if detail else ""))
    sys.stdout.flush()
    if not ok:
        FAILS.append(name)
    return ok


def digits(a, b, scale=None):
    """-log10(|a - b| / scale), scale defaults to max(|a|, |b|)."""
    d = abs(a - b)
    s = scale if scale is not None else max(abs(a), abs(b))
    if d == 0:
        return mp.mpf(999)
    if s == 0:
        return -mp.log10(d)
    return -mp.log10(d / s)


def dd(x, n=4):
    return mp.nstr(x, n)


# ----------------------------------------------------------------------------------------------------------------------------
# closed forms of the Omegas (README, #220, #263)
# ----------------------------------------------------------------------------------------------------------------------------
def omegas(upto=5):
    q = mp.log(2)
    pi = mp.pi
    z3, z5 = mp.zeta(3), mp.zeta(5)
    Om = [None] * 7
    Om[0] = -4 * q
    Om[1] = 2 * pi**2 / 3 + 8 * q**2
    Om[2] = 80 * z3 - mp.mpf(8) / 3 * pi**2 * q - mp.mpf(32) / 3 * q**3
    Om[3] = (mp.mpf(22) / 45 * pi**4 + mp.mpf(16) / 3 * pi**2 * q**2 + mp.mpf(32) / 3 * q**4 - 320 * z3 * q - 384 * z3)
    Om[4] = (-10432 * z5 - 160 * pi**2 * z3 + 640 * z3 * q**2 - mp.mpf(64) / 9 * pi**2 * q**3 - mp.mpf(88) / 45 * pi**4 * q
             - mp.mpf(128) / 15 * q**5 - mp.mpf(2432) / 45 * pi**4 + 1536 * z3 * q + 3072 * z3)
    Om[5] = (mp.mpf(256) / 45 * q**6 + mp.mpf(64) / 9 * pi**2 * q**4 - mp.mpf(2560) / 3 * z3 * q**3 + mp.mpf(176) / 45 * pi**4 * q**2
             - 3072 * z3 * q**2 + 640 * pi**2 * z3 * q + mp.mpf(9728) / 45 * pi**4 * q - 12288 * z3 * q + 41728 * z5 * q
             - mp.mpf(38812) / 315 * pi**6 + 3200 * z3**2 - mp.mpf(256) / 3 * pi**2 * z3 + 74240 * z5 + mp.mpf(39424) / 45 * pi**4
             - 36864 * z3)
    Om[6] = mp.mpf(0)
    return Om


# ----------------------------------------------------------------------------------------------------------------------------
# c_m from the master equation, by sympy:  K(eps) = (1 + eps Om) e^{eps l} - e^{-eps l}/2 - eps S e^{-eps l},  S = 1 - Rinv
# ----------------------------------------------------------------------------------------------------------------------------
KMAX = 7
_eps, _l, _Ri = sp.symbols("eps l Rinv")
_Oms = sp.symbols("Om0:7")


def c_table():
    """c_m = A_m(l) + Rinv B_m(l); returns sympy polynomials in l (coefficients in Om0..Om6) for m = 0..7."""
    Omega = sum(_eps**i * _Oms[i] for i in range(7))
    S = 1 - _Ri
    # the master equation's bracket, un-simplified: (1+eps Om)e^{eps l} - (1/2+eps) S e^{-eps l} - e^{-eps l} Rinv/2
    K = (1 + _eps * Omega) * sp.exp(_eps * _l) - (sp.Rational(1, 2) + _eps) * S * sp.exp(-_eps * _l) - sp.exp(-_eps * _l) * _Ri / 2
    ser = sp.series(K, _eps, 0, KMAX + 1).removeO()
    out = []
    for m in range(KMAX + 1):
        cm = sp.expand(ser.coeff(_eps, m))
        A = sp.expand(cm.subs(_Ri, 0))
        B = sp.expand((cm - A) / _Ri)
        assert not B.has(_Ri)
        if MUT == "c_wrong" and m == 3:
            # -2 ell^3 = +l^3/4 in l; the wrong table uses -3 ell^3 = +3 l^3/8: add l^3/8 to A_3 and subtract nothing else
            A = sp.expand(A + _l**3 / 8)
        out.append((A, B))
    return out


_CT = None


def ctab():
    global _CT
    if _CT is None:
        _CT = c_table()
    return _CT


def c_polys(Om):
    """numeric polynomial coefficient lists (mpf, by power of l) of A_m and B_m for m = 0..KMAX."""
    res = []
    for A, B in ctab():
        la = []
        for P_ in (A, B):
            pl = sp.Poly(P_, _l)
            deg = pl.degree() if not pl.is_zero else 0
            cf = [mp.mpf(0)] * (deg + 1)
            for (pw,), co in pl.terms():
                f = sp.lambdify(_Oms, co, "mpmath")
                cf[pw] = f(*Om)
            la.append(cf)
        res.append(tuple(la))
    return res


# ----------------------------------------------------------------------------------------------------------------------------
# fixed-point series (Python ints, scale 2^PB)
# ----------------------------------------------------------------------------------------------------------------------------
class Fix:
    def __init__(self, bits, N, Nthroat):
        self.P = bits
        self.ONE = 1 << bits
        self.N = N
        self.Nt = Nthroat

    def fx(self, x):
        with mp.workprec(self.P + 200):
            return int(mp.nint(mp.mpf(x) * self.ONE))

    def fl(self, i):
        with mp.workprec(self.P + 64):
            return mp.mpf(i) / self.ONE


def smul(a, b, N, P):
    rb = b[::-1]
    return [sumprod(a[: n + 1], rb[N - 1 - n:]) >> P for n in range(N)]


def spoly(a, poly, N, P):
    """series times a short polynomial (list of ints, same scale)."""
    out = [0] * N
    d = len(poly)
    for i in range(d):
        pi_ = poly[i]
        if pi_ == 0:
            continue
        for n in range(i, N):
            out[n] += pi_ * a[n - i]
    return [o >> P for o in out]


def sadd(*xs):
    return [sum(t) for t in zip(*xs)]


def shift_poly(coefs, lc):
    """mpf polynomial in l -> mpf polynomial in delta = l - lc."""
    d = len(coefs) - 1
    out = []
    for n in range(d + 1):
        s = mp.mpf(0)
        for i in range(n, d + 1):
            s += coefs[i] * mp.binomial(i, n) * lc ** (i - n)
        out.append(s)
    return out


class Hier:
    """Taylor-jet hierarchy Y_0..Y_K, L_0..L_{K-1} for given Omegas.  segments[i] covers l in [1 + i, 2 + i), center 1.5 + i;
    the throat series (center 0) covers [0, 1)."""

    def __init__(self, Om, K=7, bits=440, N=160, Nt=230, lmax=345.0, verbose=True):
        self.Om = [mp.mpf(o) for o in Om]
        self.K = K
        self.fix = Fix(bits, N, Nt)
        self.N, self.Nt = N, Nt
        self.lmax = lmax
        self.verbose = verbose
        with mp.workprec(bits + 200):
            self.cp = c_polys(self.Om)
            self._build()

    def _series_set(self, lc, Ncur, Y0, L0, throat):
        fx = self.fix.fx
        P = self.fix.P
        ONE = self.fix.ONE
        N = Ncur
        K = self.K
        lc = mp.mpf(lc)
        # e^{-l} series: E_n = (-1)^n e^{-lc} / n!
        e0 = fx(mp.exp(-lc))
        E = [0] * N
        mag = e0
        for n in range(N):
            E[n] = mag if n % 2 == 0 else -mag
            mag = mag // (n + 1)
        sig = [-e for e in E]
        sig[0] = ONE - E[0]
        # short polynomials: (lc + delta)^j / j!
        PW = []
        for j in range(K + 1):
            pc = [mp.mpf(0)] * (j + 1)
            pc[j] = mp.mpf(1)
            pcs = shift_poly(pc, lc)
            fj = mp.factorial(j)
            PW.append([fx(c / fj) for c in pcs])
        # c_m polynomials shifted
        Ap, Bp = [], []
        for m in range(K + 1):
            a, b = self.cp[m]
            Ap.append([fx(c) for c in shift_poly(a, lc)])
            Bp.append([fx(c) for c in shift_poly(b, lc)])
        Y = [None] * (K + 1)
        L = [None] * K
        for k in range(K + 1):
            if k == 0:
                Y[0] = [ONE] + [0] * (N - 1)
            else:
                f = [0] * N
                for j in range(k):
                    f = sadd(f, spoly(L[k - 1 - j], PW[j], N, P))
                Y[k] = [Y0[k]] + [f[n] // (n + 1) for n in range(N - 1)]
            if k == K:
                break
            # h_k
            if k == 0:
                h = [ONE // 2] + [0] * (N - 1)
            else:
                TA = [0] * N
                TB = [0] * N
                for j in range(k):
                    TA = sadd(TA, spoly(Y[j], Ap[k - j], N, P))
                    TB = sadd(TB, spoly(Y[j], Bp[k - j], N, P))
                inner = sadd(TB, [-t for t in L[k - 1]])
                conv = smul(E, inner, N, P)
                h = sadd(L[k - 1], [y >> 1 for y in Y[k]], TA, conv)
            if True:
                lam = [0] * N
                mu = [0] * N
                if throat:
                    # (n - 1/2) lam_n = h_n - sum_{i=2}^{n} s_i (n - i + 1) lam_{n-i+1},  s_i = sig_i
                    lam[0] = -2 * h[0]
                    for n in range(1, N):
                        acc = sumprod(sig[2: n + 1], mu[n - 1: 0: -1]) >> P if n >= 2 else 0
                        lam[n] = (2 * (h[n] - acc)) // (2 * n - 1)
                        mu[n] = n * lam[n]
                else:
                    lam[0] = L0[k]
                    s0 = sig[0]
                    for n in range(N - 1):
                        acc = (sumprod(sig[1: n + 1], mu[n: 0: -1]) >> P) if n >= 1 else 0
                        num = h[n] + (lam[n] >> 1) - acc
                        lam[n + 1] = (num << P) // (s0 * (n + 1))
                        mu[n + 1] = (n + 1) * lam[n + 1]
                L[k] = lam
            if throat and k >= 1:
                self.h0 = getattr(self, "h0", {})
                self.h0[k] = self.fix.fl(h[0])
        return Y, L, h

    def _eval(self, coefs, delta_pows):
        return sumprod(coefs, delta_pows[: len(coefs)]) >> self.fix.P

    def _build(self):
        t0 = time.time()
        P = self.fix.P
        K = self.K
        fx = self.fix.fx
        # throat
        Y0 = [0] * (K + 1)
        L0 = [0] * K
        Y, L, _ = self._series_set(0, self.Nt, Y0, L0, True)
        self.throat = (Y, L)
        # propagate to 1.5
        def pows(delta, n):
            d = fx(delta)
            pw = [self.fix.ONE]
            for i in range(1, n):
                pw.append((pw[-1] * d) >> P)
            return pw
        pw = pows(mp.mpf(3) / 2, self.Nt)
        Yv = [self._eval(Y[k], pw) for k in range(K + 1)]
        Lv = [self._eval(L[k], pw) for k in range(K)]
        self.seg = []
        nseg = int(self.lmax - 1.0) + 2
        pw1 = pows(mp.mpf(1), self.N)
        for i in range(nseg):
            Ys, Ls, _ = self._series_set(mp.mpf(3) / 2 + i, self.N, Yv, Lv, False)
            self.seg.append((Ys, Ls))
            Yv = [self._eval(Ys[k], pw1) for k in range(K + 1)]
            Lv = [self._eval(Ls[k], pw1) for k in range(K)]
        self.build_time = time.time() - t0

    def at(self, l):
        """values (mpf) of Y_0..Y_K, L_0..L_{K-1} at l (mpf)."""
        P = self.fix.P
        fx = self.fix.fx
        K = self.K
        with mp.workprec(P + 100):
            if l < 1:
                Y, L = self.throat
                N = self.Nt
                delta = l
            else:
                i = int(l - 1)
                Y, L = self.seg[i]
                N = self.N
                delta = l - (mp.mpf(3) / 2 + i)
            d = fx(delta)
            pw = [self.fix.ONE]
            for j in range(1, N):
                pw.append((pw[-1] * d) >> P)
            Yv = [self.fix.fl(self._eval(Y[k], pw)) for k in range(K + 1)]
            Lv = [self.fix.fl(self._eval(L[k], pw)) for k in range(K)]
        return Yv, Lv



# ----------------------------------------------------------------------------------------------------------------------------
# Gauss-Legendre panels in v = sqrt(l); nodes carry the hierarchy's values
# ----------------------------------------------------------------------------------------------------------------------------
_GL = {}


def gl_rule(n):
    """nodes and weights of n-point Gauss-Legendre on [-1, 1] at the current precision."""
    key = (n, mp.mp.prec)
    if key in _GL:
        return _GL[key]
    xs, ws = [], []
    for i in range(1, n + 1):
        x = mp.cos(mp.pi * (i - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < mp.mpf(10) ** (-(mp.mp.dps + 5)):
                break
        p0, p1 = mp.mpf(1), x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = n * (x * p1 - p0) / (x * x - 1)
        xs.append(x)
        ws.append(2 / ((1 - x * x) * dp * dp))
    _GL[key] = (xs, ws)
    return xs, ws


def unode(v):
    """u = tan theta at l = v^2."""
    return mp.sqrt(mp.expm1(v * v))


def make_nodes(H, n=40, width=mp.mpf(1) / 2, vend=mp.mpf(37) / 2):
    """list of node dicts over v in [0, vend] (panels of the given width), with the end node at vend."""
    xs, ws = gl_rule(n)
    npan = int(round(vend / width))
    nodes = []
    # integral of g(v) = ell^2 v / u over each full panel, for Lc2 (partial integrals use the same rule)
    def g2(v):
        return (v**4 / 4) * v / unode(v)

    cum = mp.mpf(0)
    for p in range(npan):
        a, b = p * width, (p + 1) * width
        half = (b - a) / 2
        mid = (a + b) / 2
        for x_, w_ in zip(xs, ws):
            v = mid + half * x_
            l = v * v
            u = unode(v)
            wt = half * w_ * v / u
            Y, L = H.at(l)
            # Lc2(theta) = int_0^v g2 : cumulative + partial panel
            ha = (v - a) / 2
            part = mp.mpf(0)
            for x2, w2 in zip(xs, ws):
                vv = a + ha * (1 + x2)
                part += ha * w2 * g2(vv)
            nodes.append(dict(v=v, l=l, w=wt, Y=Y, L=L, Lc2=cum + part))
        full = mp.mpf(0)
        for x_, w_ in zip(xs, ws):
            full += half * w_ * g2(mid + half * x_)
        cum += full
    Y, L = H.at(vend * vend)
    end = dict(v=vend, l=vend * vend, w=mp.mpf(0), Y=Y, L=L, Lc2=cum)
    return nodes, end


def derive(d, Om, cp):
    """fill the derived quantities of a node (mp numbers at the working precision)."""
    l = d["l"]
    ell = -l / 2
    S = -mp.expm1(-l)
    u = mp.sqrt(mp.expm1(l))
    xang = mp.atan(1 / u)
    th = mp.pi / 2 - xang
    # pi - 2 th = 2 xang
    C = mp.clsin(2, 2 * xang)
    D3 = mp.mpf(3) / 4 * mp.zeta(3) + mp.clcos(3, 2 * xang)
    q = mp.log(2)
    Lc1 = (C - 2 * th * q) / 2
    Rinv = mp.exp(-l)
    c = []
    for A, B in cp:
        a = mp.polyval(A[::-1], l)
        b = mp.polyval(B[::-1], l)
        c.append(a + b * Rinv)
    Y, L = d["Y"], d["L"]
    K = len(Y) - 1
    h = [None] * (K + 1)
    Sk = [None] * (K + 1)
    for k in range(1, K + 1):
        Sk[k] = mp.fsum(Y[j] * c[k - j] for j in range(k))
        h[k] = S * L[k - 1] + Sk[k] + Y[k] / 2
    F = []
    for k in range(len(L)):
        F.append(mp.fsum(L[k - j] * l**j / mp.factorial(j) for j in range(k + 1)))
    d.update(ell=ell, S=S, u=u, th=th, xang=xang, C=C, D3=D3, Lc1=Lc1, c=c, h=h, Sk=Sk, F=F, cot=1 / u, q=q)
    return d


# ----------------------------------------------------------------------------------------------------------------------------
# weights of derivations.md (typed from the derivation); the mutations change the *claimed* weights only
# ----------------------------------------------------------------------------------------------------------------------------
def claim():
    """the claimed constants of E4, E2W, E5 (derivations.md); the mutations change them."""
    return dict(a22=20 if MUT == "e4_20" else 22, sg1=+2 if MUT == "e4_Om1" else -2, cC=8 if MUT == "e4_C" else 16,
                gl=-30 if MUT == "e2w_g" else -32, Pd3=0 if MUT == "e2w_P" else 1, a1q=32 if MUT == "e5_a1" else 64,
                d3=0 if MUT == "e5_noD3" else 16, K=sp.Rational(1, 100) if MUT == "e5_K" else 0)


def weights(d, Om):
    """claimed and chain-derived weights at a node."""
    O0, O1, O2 = Om[0], Om[1], Om[2]
    ell, S, u, cot, C, D3, q = d["ell"], d["S"], d["u"], d["cot"], d["C"], d["D3"], d["q"]
    c1, c2, c3 = d["c"][1], d["c"][2], d["c"][3]
    out = {}
    # chain-derived weights from the c_m (not from the claims)
    out["w4"] = (4 * O0 * ell + 8 * ell * c1 + 2 * c2) / S
    # claim: E4's primitive (times u) and W4
    cl = claim()
    a22, sg1, cC = cl["a22"], cl["sg1"], cl["cC"]
    out["W4u"] = a22 * ell**2 - 8 * O0 * ell + sg1 * O1 + cC * u * C
    out["W4"] = cot * (a22 * ell**2 - 8 * O0 * ell + sg1 * O1) + cC * C
    # claim: g4 and its primitive P
    g4 = cl["gl"] * ell**2 + 16 * O0 * ell + 4 * O1 - 32 * u * C
    out["g4"] = g4
    Ptrue = mp.mpf(32) / 3 * ell**3 - 8 * O0 * ell**2 - 4 * O1 * ell - 32 * (u * C + ell**2 + 2 * q * ell - D3 / 2)
    Pclaim = (mp.mpf(32) / 3 * ell**3 - 8 * O0 * ell**2 - 4 * O1 * ell
              - 32 * (u * C + ell**2 + 2 * q * ell - cl["Pd3"] * D3 / 2))
    out["Pt"], out["P"] = Ptrue, Pclaim
    out["w4p"] = (8 * O0 * ell * c1 + 8 * ell * c2 + 2 * c3 - Ptrue) / S
    # claim: E5's primitive W
    a1 = 8 * O1 + 8 * O0**2 + cl["a1q"] * q
    a2 = 32 - 28 * O0
    a3 = -mp.mpf(20) / 3
    Lc1, Lc2, th = d["Lc1"], d["Lc2"], d["th"]
    W = (-2 * O2 * cot - a1 * (ell * cot + th) - a2 * (ell**2 * cot + 2 * Lc1) - a3 * (ell**3 * cot + 3 * Lc2)
         + 32 * (-C * ell + 2 * Lc2 + 2 * q * Lc1) + cl["d3"] * D3 * cot - 20 * Lc2 - 8 * O0 * Lc1)
    W += mp.mpf(cl["K"].p) / cl["K"].q if cl["K"] != 0 else 0
    out["W"] = W
    return out


def integrands(d, Om):
    """all integrand values at a node, by name."""
    Y, L, F, h, c, Sk = d["Y"], d["L"], d["F"], d["h"], d["c"], d["Sk"]
    ell, S, u, cot, l = d["ell"], d["S"], d["u"], d["cot"], d["l"]
    O0 = Om[0]
    wt = weights(d, Om)
    r = {}
    for k in range(2, 8):
        r["Gdef%d" % k] = 2 * (h[k] - Om[k - 1]) / S
    # E1 at k
    for k in (6, 7):
        r["E1l_%d" % k] = r["Gdef%d" % k]
        r["E1r_%d" % k] = 4 * L[k - 1] + 2 * (F[k - 1] - L[k - 1]) + 2 * (Sk[k] - Om[k - 1]) / S
    # E2 at m
    for m in range(2, 7):
        r["E2l_%d" % m] = L[m]
        r["E2r_%d" % m] = 2 * ell * h[m] / S
    # E3 at m (Y_m weight 4 ell + 2 c1 over S = 2 P_E')
    for m in range(3, 8):
        r["E3l_%d" % m] = Y[m] * (4 * ell + 2 * c[1]) / S
        r["E3r_%d" % m] = -4 * (ell - O0) * F[m - 1]
    # E4 at k = 6 (Y4, F3) and the k = 7 use (Y5, F4)
    r["E4l_6"] = Y[4] * wt["w4"]
    r["E4r_6"] = -2 * wt["W4u"] * F[3]
    r["E4l_7"] = Y[5] * wt["w4"]
    r["E4r_7"] = -2 * wt["W4u"] * F[4]
    # E2W: weight 12 ell^2 (P_a = -4 ell^3), and g4 with P
    r["E2Wal"] = 12 * ell**2 * L[4]
    r["E2War"] = -2 * (-4 * ell**3) * h[4] / S
    r["E2Wbl"] = wt["g4"] * L[4]
    r["E2Wbr"] = -2 * wt["P"] * h[4] / S
    # E5
    r["E5l"] = Y[4] * wt["w4p"]
    r["E5r"] = -2 * wt["W"] * u * F[3]
    # reduced forms
    def Spp(n):
        return mp.fsum(Y[j] * c[n - j] for j in range(4))
    S4, S5, S6, S7 = Spp(4), Spp(5), Spp(6), Spp(7)
    l2, l3, l4, l5, l6 = l**2, l**3, l**4, l**5, l**6
    f6tail = L[3] * l3 / 6 + L[2] * l4 / 24 + L[1] * l5 / 120 + L[0] * l6 / 720   # F6 - L6 - L5 l - L4 l^2/2
    f5t = L[3] * l2 / 2 + L[2] * l3 / 6 + L[1] * l4 / 24 + L[0] * l5 / 120
    f4t = L[3] * l + L[2] * l2 / 2 + L[1] * l3 / 6 + L[0] * l4 / 24
    f3 = F[3]
    # G6 reduced
    t6 = (8 * O0 * ell * L[3] + 8 * O0 * ell * S4 / S + 8 * ell * S5 / S
          + 2 * (L[3] * l2 / 2 + L[2] * l3 / 6 + L[1] * l4 / 24 + L[0] * l5 / 120)
          - 4 * (ell - O0) * (L[3] * l + L[2] * l2 / 2 + L[1] * l3 / 6 + L[0] * l4 / 24)
          + 2 * (S6 - Om[5]) / S - 2 * wt["W4u"] * f3)
    if MUT == "g6r_slope":
        t6 -= 8 * ell * Om[4] / S      # the Om4 of Y0 c5, the only Om4 left in the reduced form
    r["G6Rr"] = t6
    r["G6Rl"] = r["Gdef6"]
    # G7 reduced
    t7 = (8 * O0 * ell * S5 / S + 8 * ell * S6 / S + 2 * (S7 - Om[6]) / S - 2 * wt["P"] * S4 / S - 2 * wt["P"] * L[3]
          + 2 * f6tail - 4 * (ell - O0) * f5t - 2 * wt["W4u"] * f4t - 2 * wt["W"] * u * f3)
    if MUT == "g7r_slope":
        t7 -= 8 * ell * Om[5] / S      # the Om5 of Y0 c6
    r["G7Rr"] = t7
    r["G7Rl"] = r["Gdef7"]
    return r


def boundary_terms(e, Om):
    """predicted boundary terms of every by-parts step at the cutoff node e (derive() applied)."""
    Y, L, ell, cot = e["Y"], e["L"], e["ell"], e["cot"]
    O0 = Om[0]
    wt = weights(e, Om)
    b = {}
    for k in (6, 7):
        b["E1_%d" % k] = -cot * Y[k]
    for m in range(2, 7):
        b["E2_%d" % m] = -ell * L[m] * cot
    for m in range(3, 8):
        b["E3_%d" % m] = 2 * cot * (ell - O0) * Y[m]
    b["E4_6"] = wt["W4"] * Y[4]
    b["E4_7"] = wt["W4"] * Y[5]
    b["E2Wa"] = -4 * ell**3 * L[4] * cot
    b["E2Wb"] = wt["P"] * L[4] * cot
    b["E5"] = wt["W"] * Y[4]
    b["G6R"] = b["E1_6"] + 4 * b["E2_5"] + b["E3_5"] + 4 * O0 * b["E2_4"] + b["E4_6"]
    b["G7R"] = (b["E1_7"] + 4 * b["E2_6"] + b["E3_6"] + 4 * O0 * b["E2_5"] + b["E4_7"] + b["E2Wb"] + b["E5"])
    return b


# ----------------------------------------------------------------------------------------------------------------------------
# runs on the hierarchy
# ----------------------------------------------------------------------------------------------------------------------------
DPS = 100
DEFAULT = dict(bits=440, N=160, Nt=230, lmax=345.0, n=40, width=mp.mpf(1) / 2, vend=mp.mpf(37) / 2)
BAR_RAW = 30      # digits, identities (the issue's acceptance)
BAR_COR = 40      # digits, identities after the predicted boundary term at the cutoff is subtracted (declared before the first run)
REQ = {  # highest Omega index that must be at its root for the identity to apply
    "E1_6": 3, "E1_7": 4, "E2_2": 0, "E2_3": 1, "E2_4": 2, "E2_5": 3, "E2_6": 4,
    "E3_3": 0, "E3_4": 1, "E3_5": 2, "E3_6": 3, "E3_7": 4,
    "E4_6": 1, "E4_7": 2, "E2Wa": 2, "E2Wb": 2, "E5": 1, "G6R": 3, "G7R": 4,
}


def parse_shift(s):
    out = {}
    if s:
        for tok in s.split(","):
            i, v = tok.split(":")
            out[int(i)] = mp.mpf(v)
    return out


def shifted_omegas(shift):
    with mp.workprec(700):
        Om = omegas()
        for i, v in shift.items():
            Om[i] = Om[i] + v
    return Om


def tag_of(shift, cfg):
    s = "_".join("%d_%s" % (i, mp.nstr(v, 6)) for i, v in sorted(shift.items())) or "roots"
    if MUT == "c_wrong":
        s += "_cwrong"
    return "%s_b%d_N%d_n%d_w%s" % (s, cfg["bits"], cfg["N"], cfg["n"], mp.nstr(cfg["width"], 3))


def get_nodes(shift, cfg, cachedir, log=True):
    """(nodes, end, Om, H or None).  Derived nodes are pickled by configuration."""
    Om = shifted_omegas(shift)
    path = os.path.join(cachedir, tag_of(shift, cfg) + ".pkl") if cachedir else None
    with mp.workdps(DPS):
        cp = c_polys(Om)
    if path and os.path.exists(path):
        with open(path, "rb") as f:
            nodes, end = pickle.load(f)
        if log:
            print("# nodes loaded from", path)
        return nodes, end, Om, cp
    t = time.time()
    with mp.workprec(cfg["bits"] + 300):
        H = Hier(Om, K=7, bits=cfg["bits"], N=cfg["N"], Nt=cfg["Nt"], lmax=cfg["lmax"])
    if log:
        print("# hierarchy built in %.1f s (N=%d, %d bits)" % (time.time() - t, cfg["N"], cfg["bits"]))
    with mp.workdps(DPS):
        t = time.time()
        nodes, end = make_nodes(H, n=cfg["n"], width=cfg["width"], vend=cfg["vend"])
        for d in nodes:
            derive(d, Om, cp)
        derive(end, Om, cp)
        if log:
            print("# %d nodes made and derived in %.1f s" % (len(nodes), time.time() - t))
    if path:
        with open(path, "wb") as f:
            pickle.dump((nodes, end), f)
    return nodes, end, Om, cp


def accumulate(nodes, Om):
    acc, absacc = {}, {}
    for d in nodes:
        r = integrands(d, Om)
        for k, v in r.items():
            acc[k] = acc.get(k, 0) + d["w"] * v
            absacc[k] = absacc.get(k, 0) + d["w"] * abs(v)
    return acc, absacc


def applicable(name, shift):
    if MUT == "hyp_G4":
        return True
    j = REQ[name]
    return all(i > j for i in shift)


def run_exch(args):
    cfg = dict(DEFAULT)
    shift = parse_shift(args.shift)
    if MUT == "hyp_G4":
        shift = {2: mp.mpf("0.001")}
    with mp.workdps(DPS):
        nodes, end, Om, cp = get_nodes(shift, cfg, args.cache)
        t = time.time()
        acc, absacc = accumulate(nodes, Om)
        b = boundary_terms(end, Om)
        print("# integrals accumulated in %.1f s; %d nodes; shift %s" % (time.time() - t, len(nodes), {i: mp.nstr(v, 5) for i, v in shift.items()}))
        # ---- G_k from the definition against (L_k + 2 Om_{k-1}) / u at the cutoff (the ODE against the integral)
        for k in range(2, 8):
            g_int = acc["Gdef%d" % k]
            g_end = (end["L"][k] + 2 * Om[k - 1]) / end["u"] if k < 7 else None
            if k < 7:
                dg = digits(g_int, g_end, absacc["Gdef%d" % k])
                check("H3_G%d" % k, dg >= BAR_COR, "G_%d: quadrature %s, (L_k + 2 Om_{k-1})/u at the cutoff %s, digits %s" % (k, dd(g_int, 8), dd(g_end, 8), dd(dg, 3)))
            else:
                print("# G7 quadrature %s (no L_7 in the hierarchy)" % dd(g_int, 10))
        # ---- H2: G_k vanishes at the closed forms; the slope check when a constant is shifted
        for k in range(2, 8):
            g_int = acc["Gdef%d" % k]
            idx = k - 2
            if idx in shift:
                pred = -4 * mp.pi * shift[idx]
                check("SLOPE_G%d" % k, abs(g_int - pred) <= mp.mpf(10) ** -BAR_RAW * abs(pred),
                      "G_%d = %s against -4 pi delta = %s" % (k, dd(g_int, 12), dd(pred, 12)))
            elif all(i > k - 2 for i in shift) or MUT == "hyp_G4":
                tag = "H2_G%d" % k
                if k <= 6:
                    check(tag, abs(g_int) <= mp.mpf(10) ** -BAR_RAW, "|G_%d| = %s at the closed forms (bar 1e-%d)" % (k, dd(abs(g_int), 4), BAR_RAW))
                else:
                    print("# G_7 at #263's Omega_5: %s (reported)" % dd(g_int, 4))
        # ---- identities
        rows = [
            ("E1_6", "E1l_6", "E1r_6", "E1_6"), ("E1_7", "E1l_7", "E1r_7", "E1_7"),
            ("E2_2", "E2l_2", "E2r_2", "E2_2"), ("E2_3", "E2l_3", "E2r_3", "E2_3"), ("E2_4", "E2l_4", "E2r_4", "E2_4"),
            ("E2_5", "E2l_5", "E2r_5", "E2_5"), ("E2_6", "E2l_6", "E2r_6", "E2_6"),
            ("E3_3", "E3l_3", "E3r_3", "E3_3"), ("E3_4", "E3l_4", "E3r_4", "E3_4"), ("E3_5", "E3l_5", "E3r_5", "E3_5"),
            ("E3_6", "E3l_6", "E3r_6", "E3_6"), ("E3_7", "E3l_7", "E3r_7", "E3_7"),
            ("E4_6", "E4l_6", "E4r_6", "E4_6"), ("E4_7", "E4l_7", "E4r_7", "E4_7"),
            ("E2Wa", "E2Wal", "E2War", "E2Wa"), ("E2Wb", "E2Wbl", "E2Wbr", "E2Wb"),
            ("E5", "E5l", "E5r", "E5"), ("G6R", "G6Rl", "G6Rr", "G6R"), ("G7R", "G7Rl", "G7Rr", "G7R"),
        ]
        res = {}
        for name, ln, rn, bn in rows:
            lhs, rhs = acc[ln], acc[rn]
            scale = max(absacc[ln], absacc[rn])
            raw = digits(lhs, rhs, scale)
            cor = digits(lhs - rhs, b[bn], scale)
            res[name] = (raw, cor)
            ok = (raw >= BAR_RAW and cor >= BAR_COR)
            if applicable(name, shift):
                check(name, ok, "lhs %s rhs %s | raw %s digits, with the boundary term %s at the cutoff %s digits (scale %s)"
                      % (dd(lhs, 12), dd(rhs, 12), dd(raw, 3), dd(b[bn], 4), dd(cor, 3), dd(scale, 3)))
            else:
                print("INFO %s (hypothesis not met, no verdict) : lhs %s rhs %s raw %s digits, corrected %s" % (name, dd(lhs, 12), dd(rhs, 12), dd(raw, 3), dd(cor, 3)))
        if set(shift) == {2} and abs(shift[2] - mp.mpf("0.001")) < mp.mpf(10) ** -12 and MUT is None:
            ok = (res["E2Wb"][0] < 8 and res["E2Wb"][1] >= BAR_COR and res["E4_6"][0] >= BAR_RAW and res["E5"][0] >= BAR_RAW)
            check("HYP", ok, "Om2 + 1e-3 (G_4 = -4 pi 1e-3 != 0): E2Wb fails by %s digits only (bar: fails, < 8) yet equals its boundary term at the cutoff to %s digits; E4_6 holds to %s, E5 to %s digits"
                  % (dd(res["E2Wb"][0], 3), dd(res["E2Wb"][1], 3), dd(res["E4_6"][0], 3), dd(res["E5"][0], 3)))
    return acc, absacc, b, Om, nodes, end


def main():
    global MUT
    ap = argparse.ArgumentParser()
    ap.add_argument("parts", nargs="*", default=["sym", "exch"])
    ap.add_argument("--mutate", default=None)
    ap.add_argument("--shift", default="")
    ap.add_argument("--cache", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--ctldir", default=".")
    ap.add_argument("--ctlprefix", default="ctl")
    args = ap.parse_args()
    MUT = args.mutate
    if MUT is not None and MUT not in MUTATIONS:
        sys.exit("unknown mutation " + MUT)
    if args.out:
        class Tee:
            def __init__(self, f):
                self.f, self.o = f, sys.stdout
            def write(self, s):
                self.f.write(s)
                self.o.write(s)
            def flush(self):
                self.f.flush()
                self.o.flush()
        sys.stdout = Tee(open(args.out, "x"))
    print("# identity_review.py parts=%s mutate=%s shift=%s" % (",".join(args.parts), MUT, args.shift))
    for part in args.parts:
        print("\n== part %s ==" % part)
        t = time.time()
        {"sym": run_sym, "exch": run_exch, "shifts": run_shifts, "bnd": run_bnd, "masters": run_masters, "masters2": run_masters2, "fp": run_fp, "ctlsum": run_ctlsum}[part](args)
        print("# part %s: %.1f s" % (part, time.time() - t))
    print("\n%d check(s) failed%s" % (len(FAILS), ": " + ", ".join(FAILS) if FAILS else ""))
    sys.exit(1 if FAILS else 0)


def run_sym(args):
    """symbolic checks (sympy residual 0) of the algebra of derivations.md."""
    th = sp.symbols("theta", positive=True)
    ell = sp.log(sp.cos(th))
    u, S, cot = sp.tan(th), sp.sin(th) ** 2, sp.cot(th)
    q = sp.log(2)
    O = sp.symbols("O0:7")
    Cf, D3f, Lc1f, Lc2f = (sp.Function(n)(th) for n in ("C", "D3", "Lc1", "Lc2"))
    # derivative rules: C' = 2 (ell + ln 2), D3' = 2 C, Lc1' = ell, Lc2' = ell^2
    def D(expr):
        e = sp.diff(expr, th)
        e = e.subs({sp.Derivative(Cf, th): 2 * (ell + q), sp.Derivative(D3f, th): 2 * Cf,
                    sp.Derivative(Lc1f, th): ell, sp.Derivative(Lc2f, th): ell ** 2})
        return e

    def zero(expr):
        # the functions C, D3, Lc1, Lc2, and ell, become plain symbols, then the rest is trigonometric and rational
        e = expr.subs({Cf: sp.Symbol("Cs"), D3f: sp.Symbol("D3s"), Lc1f: sp.Symbol("L1s"), Lc2f: sp.Symbol("L2s")})
        e = sp.simplify(sp.expand_log(e, force=True))
        return sp.simplify(e) == 0

    cl = claim()
    Ro = sp.symbols("Ro")                      # Ro = ell as a symbol
    # --- c_m from the master equation against the closed form in ell
    Sx = sp.Symbol("Sx")
    ok = True
    for m in range(1, KMAX + 1):
        A, B = ctab()[m]
        cm = sp.expand((A + B * (1 - Sx)).subs(_l, -2 * Ro))
        cf = (sum(_Oms[i] * (-2 * Ro) ** (m - 1 - i) / sp.factorial(m - 1 - i) for i in range(m))
              + (1 - sp.Rational(1, 2) * (-1) ** m) * (-2 * Ro) ** m / sp.factorial(m) - Sx * (2 * Ro) ** (m - 1) / sp.factorial(m - 1))
        ok = ok and sp.expand(cm - cf) == 0
    check("SY_c", ok, "c_1 .. c_7 from the expansion of (1+eps Om)e^{eps l} - (1/2+eps)S e^{-eps l} - e^{-eps l}/(2R) equal the closed form in ell")
    c = [sp.expand((ctab()[m][0] + ctab()[m][1] * (1 - Sx)).subs(_l, -2 * Ro)) for m in range(KMAX + 1)]
    ok = (sp.expand(c[1] - (_Oms[0] - 3 * Ro - Sx)) == 0 and sp.expand(c[2] - (_Oms[1] - 2 * _Oms[0] * Ro + Ro**2 - 2 * Sx * Ro)) == 0
          and sp.expand(c[3] - (_Oms[2] - 2 * _Oms[1] * Ro + 2 * _Oms[0] * Ro**2 - 2 * Ro**3 - 2 * Sx * Ro**2)) == 0)
    check("SY_cS", ok, "c_1 = Om0 - 3 ell - S, c_2 = Om1 - 2 Om0 ell + ell^2 - 2 S ell, c_3 = Om2 - 2 Om1 ell + 2 Om0 ell^2 - 2 ell^3 - 2 S ell^2")
    Ri = _Ri
    c6 = (_Oms[0] * _l**5 / 120 + _Oms[1] * _l**4 / 24 + _Oms[2] * _l**3 / 6 + _Oms[3] * _l**2 / 2 + _Oms[4] * _l + _Oms[5]
          + _l**6 / 1440 + _l**5 / 120 - _l**5 * Ri / 120)
    c7 = (_Oms[0] * _l**6 / 720 + _Oms[1] * _l**5 / 120 + _Oms[2] * _l**4 / 24 + _Oms[3] * _l**3 / 6 + _Oms[4] * _l**2 / 2 + _Oms[5] * _l
          + _Oms[6] + _l**7 / 3360 - _l**6 / 720 + _l**6 * Ri / 720)
    ok = all(sp.expand(ctab()[m][0] + ctab()[m][1] * Ri - cm) == 0 for m, cm in ((6, c6), (7, c7)))
    check("SY_c67", ok, "c_6 and c_7 equal the definitions printed in omega4.md / omega5.md (read as definitions)")
    ok = all(sp.expand((ctab()[m][0] + ctab()[m][1]).subs(_l, 0) - _Oms[m - 1]) == 0 for m in range(1, KMAX + 1))
    check("SY_c0", ok, "c_m(l = 0, R = 1) = Om_{m-1}  (h_k(1) = Om_{k-1})")
    # --- (L cot th)' = 2 h / S from the L equation (R-1) dL/dR - L/2 = h
    Lf = sp.Function("Lf")(th)
    hf = sp.Function("hf")(th)
    R = 1 / sp.cos(th) ** 2
    dLdR = sp.diff(Lf, th) / sp.diff(R, th)
    Lprime = sp.solve(sp.Eq((R - 1) * dLdR - Lf / 2, hf), sp.diff(Lf, th))[0]
    res = sp.diff(Lf * cot, th).subs(sp.diff(Lf, th), Lprime) - 2 * hf / S
    check("SY_L", sp.simplify(res) == 0, "(L cot th)' = 2 h / S from (R - 1) dL/dR - L/2 = h")
    # --- weights of the chain
    O0, O1, O2 = O[0], O[1], O[2]
    c1, c2, c3 = (O0 - 3 * ell - S), (O1 - 2 * O0 * ell + ell**2 - 2 * S * ell), (O2 - 2 * O1 * ell + 2 * O0 * ell**2 - 2 * ell**3 - 2 * S * ell**2)
    w4 = (4 * O0 * ell + 8 * ell * c1 + 2 * c2) / S
    check("SY_w4", zero(w4 - ((2 * O1 + 8 * O0 * ell - 22 * ell**2) / S - 12 * ell)), "w4 = (4 Om0 ell + 8 ell c1 + 2 c2)/S = (2 Om1 + 8 Om0 ell - 22 ell^2)/S - 12 ell")
    PE = cot * (ell - O0)
    check("SY_E3", zero((4 * ell + 2 * c1) / S - 2 * D(PE)), "E3: (4 ell + 2 c1)/S = 2 P_E', P_E = cot th (ell - Om0);  u P_E = ell - Om0")
    # E4: W4' = w4 with Om0 = -4 ln 2
    W4 = cot * (cl["a22"] * ell**2 - 8 * O0 * ell + cl["sg1"] * O1) + cl["cC"] * Cf
    d = (D(W4) - w4).subs(O0, -4 * q)
    check("SY_W4", zero(d), "W4 = cot th (22 ell^2 - 8 Om0 ell - 2 Om1) + 16 C has W4' = w4 at Om0 = -4 ln 2")
    # g4 from the chain, and P' = u g4
    lq = -2 * ell
    g4chain = 8 * O0 * ell + lq**2 + 8 * ell * (ell - O0) - 2 * (22 * ell**2 - 8 * O0 * ell - 2 * O1 + 16 * u * Cf)
    g4 = cl["gl"] * ell**2 + 16 * O0 * ell + 4 * O1 - 32 * u * Cf
    check("SY_g4", zero(g4chain - (-32 * ell**2 + 16 * O0 * ell + 4 * O1 - 32 * u * Cf)), "g4 = 8 Om0 ell + l^2 + 8 ell (ell - Om0) - 2 (22 ell^2 - 8 Om0 ell - 2 Om1 + 16 u C)")
    P = sp.Rational(32, 3) * ell**3 - 8 * O0 * ell**2 - 4 * O1 * ell - 32 * (u * Cf + ell**2 + 2 * q * ell - cl["Pd3"] * D3f / 2)
    check("SY_P", zero(D(P) - u * g4), "P' = u g4 for P = (32/3) ell^3 - 8 Om0 ell^2 - 4 Om1 ell - 32 (u C + ell^2 + 2 ln2 ell - D3/2)")
    Pa = -4 * ell**3
    check("SY_Pa", zero(D(Pa) - u * 12 * ell**2), "P_a' = u 12 ell^2, P_a = -4 ell^3")
    # w4' from the chain with the true P, and W' = w4'
    Ptrue = sp.Rational(32, 3) * ell**3 - 8 * O0 * ell**2 - 4 * O1 * ell - 32 * (u * Cf + ell**2 + 2 * q * ell - D3f / 2)
    w4p = (8 * O0 * ell * c1 + 8 * ell * c2 + 2 * c3 - Ptrue) / S
    a1t, a2t, a3t = 8 * O1 + 8 * O0**2 + 64 * q, 32 - 28 * O0, -sp.Rational(20, 3)
    w4p_s = ((2 * O2 + a1t * ell + a2t * ell**2 + a3t * ell**3 + 32 * u * Cf - 16 * D3f) / S - 20 * ell**2 - 8 * O0 * ell)
    check("SY_w4p", zero(w4p - w4p_s), "w4' = (8 Om0 ell c1 + 8 ell c2 + 2 c3 - P)/S = [2 Om2 + a1 ell + a2 ell^2 + a3 ell^3 + 32 u C - 16 D3]/S - 20 ell^2 - 8 Om0 ell")
    a1 = 8 * O1 + 8 * O0**2 + cl["a1q"] * q
    Wf = (-2 * O2 * cot - a1 * (ell * cot + th) - a2t * (ell**2 * cot + 2 * Lc1f) - a3t * (ell**3 * cot + 3 * Lc2f)
          + 32 * (-Cf * ell + 2 * Lc2f + 2 * q * Lc1f) + cl["d3"] * D3f * cot - 20 * Lc2f - 8 * O0 * Lc1f)
    check("SY_W", zero(D(Wf) - w4p_s), "W' = w4'  (with Lc1' = ell, Lc2' = ell^2, C' = 2 (ell + ln2), D3' = 2C; the C cot th terms cancel)")
    # W(pi/2) = K0 + K:  K0 = -a1 pi/2 - 2 a2 J1 + 64 J2 + (64 ln2 - 8 Om0) J1
    J1 = -sp.pi / 2 * q
    J2 = sp.pi / 2 * (q**2 + sp.pi**2 / 12)
    K0 = -a1t * sp.pi / 2 - 2 * a2t * J1 + (-3 * a3t + 44) * J2 + (64 * q - 8 * O0) * J1
    K0r = sp.simplify(K0.subs({O0: -4 * q, O1: 2 * sp.pi**2 / 3 + 8 * q**2}))
    Kfull = sp.simplify(K0.subs(O0, -4 * q))
    check("SY_K", sp.simplify(K0r + cl["K"]) == 0, "W(pi/2) = K0 + K with K0 = %s at the roots of Om0, Om1 (and %s at the root of Om0 only)" % (K0r, Kfull))
    # --- boundary expansions at th = 0: series of C and D3 from C' and D3'
    x = th
    ellser = sp.series(ell, th, 0, 9).removeO()
    Cser = 2 * q * x + 2 * sp.integrate(ellser, x)
    D3ser = 2 * sp.integrate(Cser, x)
    def ser(e, n=3):
        e = e.subs({Cf: Cser, D3f: D3ser})
        return sp.series(sp.expand(e), th, 0, n).removeO()
    L3_0 = -2 * O2
    Y4 = -2 * O2 * (-2 * ellser)         # Y4 = F3(0) l = -2 Om2 l
    Y5 = -2 * O[3] * (-2 * ellser)
    ok = True
    for nm, expr, pred in (("W4 Y4", W4 * Y4, 4 * O1 * O2 * th), ("W4 Y5", W4 * Y5, 4 * O1 * O[3] * th), ("W Y4", Wf.subs({Lc1f: sp.integrate(ellser, x), Lc2f: sp.integrate(ellser**2, x)}) * Y4, 4 * O2**2 * th),
                           ("P L4 cot", P * (-2 * O[3]) * cot, -4 * O1 * O[3] * th)):
        e = sp.series(sp.expand(ser(expr, 2)), th, 0, 2).removeO()
        ok = ok and sp.simplify(e - pred) == 0
    check("SY_bnd0", ok, "throat leading terms: W4 Y4 ~ 4 Om1 Om2 th, W4 Y5 ~ 4 Om1 Om3 th, W Y4 ~ 4 Om2^2 th, P L4 cot ~ -4 Om1 Om3 th")
    # E1 / E2 / E3 boundary leading terms at the throat
    ok = True
    for nm, expr, pred in (("E1", -cot * (-2 * O[2] * (-2 * ellser)), 2 * O[2] * th), ("E2", -ell * (-2 * O[3]) * cot, -O[3] * th),
                           ("E3", 2 * cot * (ell - O0) * (-2 * O[2] * (-2 * ellser)), 4 * O0 * O[2] * th)):
        e = sp.series(sp.expand(expr.subs(ell, ellser)), th, 0, 2).removeO()
        ok = ok and sp.simplify(e - pred) == 0
    check("SY_bnd1", ok, "throat leading terms of E1, E2, E3's boundary terms: -cot Y_k ~ 2 Om_{k-2} th, -ell L cot ~ -Om th, 2 P_E Y ~ 4 Om0 Om_{m-2} th")
    # --- moments: J_n from the Beta function
    a = sp.Symbol("a")
    Bf = sp.sqrt(sp.pi) / 2 * sp.gamma((a + 1) / 2) / sp.gamma(a / 2 + 1)
    j1 = sp.simplify(sp.diff(Bf, a).subs(a, 0))
    j2 = sp.simplify(sp.diff(Bf, a, 2).subs(a, 0))
    ok = sp.simplify(j1 - J1) == 0 and sp.simplify(sp.expand_func(j2) - J2) == 0
    check("SY_J", ok, "J1 = -(pi/2) ln2 and J2 = (pi/2)(ln^2 2 + pi^2/12) from B(a) = int cos^a")
    Bt = sp.pi / 2 * sp.gamma(1 + a) / sp.gamma(1 + a / 2) ** 2
    # B~(a) = 2^a B(a)
    ok = sp.simplify(sp.expand_func(sp.gammasimp(2**a * Bf - Bt))) == 0
    print("# B~(a) = (pi/2) Gamma(1+a)/Gamma(1+a/2)^2 vs 2^a B(a) (duplication): %s" % ("equal" if ok else "not simplified by sympy; checked numerically in masters"))


def run_bnd(args):
    """planted truth of the hierarchy (H0), boundary terms at both ends (BND), growth (GROW), convergence (H1)."""
    cfg = dict(DEFAULT)
    with mp.workdps(DPS):
        nodes, end, Om, cp = get_nodes({}, cfg, args.cache)
        # ---- H0: Y1, L1, Y2, L2 against README's closed forms, 20 nodes with l <= 60
        sel = [nodes[i] for i in range(0, len(nodes), 1) if nodes[i]["l"] <= 60]
        step = max(1, len(sel) // 20)
        sel = sel[::step][:20]
        worst = mp.mpf(0)
        for d in sel:
            l, u = d["l"], d["u"]
            K_ = mp.quad(lambda t: mp.log(1 + t * t) / (1 + t * t), [0, 1, u] if u > 1 else [0, u])
            Y1, L1 = -l, -2 * Om[0] - 2 * l
            Y2 = -2 * Om[0] * l - mp.mpf(3) / 2 * l**2
            L2 = -2 * Om[1] + 2 * Om[0] * l + 4 * l**2 - 8 * Om[0] * u * mp.atan(u) - 16 * u * K_
            for a_, b_ in ((d["Y"][1], Y1), (d["L"][1], L1), (d["Y"][2], Y2), (d["L"][2], L2)):
                worst = max(worst, abs(a_ - b_))
        check("H0", worst <= mp.mpf(10) ** -35, "Y1, L1, Y2, L2 against README's closed forms at %d nodes, l from %s to %s: worst difference %s (bar 1e-35)"
              % (len(sel), dd(sel[0]["l"], 3), dd(sel[-1]["l"], 4), dd(worst, 3)))
        # ---- BND at the throat: the first two nodes, b / (theta * predicted coefficient) -> 1
        O0, O1, O2, O3 = Om[0], Om[1], Om[2], Om[3]
        pred = {"E1_6": 2 * Om[4], "E1_7": 2 * Om[5], "E2_4": -Om[3], "E2_5": -Om[4], "E2_6": -Om[5],
                "E3_5": 4 * O0 * Om[3], "E3_6": 4 * O0 * Om[4], "E4_6": 4 * O1 * O2, "E4_7": 4 * O1 * O3, "E2Wb": -4 * O1 * O3, "E5": 4 * O2**2}
        worst = {}
        first = sorted(nodes, key=lambda z: z["l"])[:2]
        for d in first:
            b = boundary_terms(d, Om)
            for k, pv in pred.items():
                r = b[k] / (d["th"] * pv)
                worst[k] = max(worst.get(k, 0), abs(r - 1))
        ok = all(v < mp.mpf(10) ** -4 for v in worst.values())
        check("BND_throat", ok, "b/(theta * predicted) - 1 at theta = %s, %s: worst over %d terms %s (bar 1e-4); e.g. E4_6 %s, E5 %s, E2Wb %s"
              % (dd(first[0]["th"], 3), dd(first[1]["th"], 3), len(pred), dd(max(worst.values()), 3), dd(worst["E4_6"], 3), dd(worst["E5"], 3), dd(worst["E2Wb"], 3)))
        # ---- BND at the outer end: b / (x l^q) at l near 30, 60, 100
        qexp = {"E1_6": 6, "E1_7": 7, "E2_5": 6, "E2_6": 7, "E3_5": 6, "E3_6": 7, "E4_6": 6, "E4_7": 7, "E2Wb": 7, "E5": 7}
        rows = []
        bad = False
        for target in (30, 60, 100):
            d = min(nodes, key=lambda z: abs(z["l"] - target))
            b = boundary_terms(d, Om)
            line = []
            for k, qq in qexp.items():
                r = b[k] / (d["xang"] * d["l"] ** qq)
                line.append("%s %s" % (k, dd(r, 3)))
                bad = bad or abs(r) > 1000
            rows.append("l = %s: %s" % (dd(d["l"], 4), ", ".join(line)))
        check("BND_outer", not bad, "b/(x l^q) bounded by 1e3 at l = 30, 60, 100:")
        for r in rows:
            print("#   " + r)
        # ---- GROW: Y_k / l^k and L_k / l^k near l = 256 (reported)
        d = min(nodes, key=lambda z: abs(z["l"] - 256))
        print("# GROW at l = %s: Y_k/l^k k=1..7: %s" % (dd(d["l"], 4), ", ".join(dd(d["Y"][k] / d["l"] ** k, 3) for k in range(1, 8))))
        print("# GROW at l = %s: L_k/l^k k=1..6: %s" % (dd(d["l"], 4), ", ".join(dd(d["L"][k] / d["l"] ** k, 3) for k in range(1, 7))))
        # ---- H1: the same integrals on a second setting
        cfg2 = dict(bits=500, N=190, Nt=260, lmax=345.0, n=44, width=mp.mpf(2) / 5, vend=mp.mpf(92) / 5)
        nodes2, end2, Om2, cp2 = get_nodes({}, cfg2, args.cache)
        acc1, ab1 = accumulate(nodes, Om)
        acc2, ab2 = accumulate(nodes2, Om2)
        worst, wname = mp.mpf(0), ""
        for k in ("E4l_6", "E4r_6", "E4l_7", "E2Wbl", "E2Wbr", "E5l", "E5r", "G6Rr", "G7Rr", "E1r_7", "E3l_7"):
            r = abs(acc1[k] - acc2[k]) / ab1[k]
            if r > worst:
                worst, wname = r, k
        check("H1", worst < mp.mpf(10) ** -40, "second setting (N=190, 500 bits, 44 nodes on panels of 0.4, l_end 339): largest relative change %s (%s), bar 1e-40" % (dd(worst, 3), wname))


def run_masters(args):
    mp.mp.dps = 75
    q, pi = mp.log(2), mp.pi
    z3, z5, z2, z4 = mp.zeta(3), mp.zeta(5), mp.zeta(2), mp.zeta(4)
    half = pi / 2

    def Cf(t):
        return mp.clsin(2, pi - 2 * t)

    def D3f(t):
        return mp.mpf(3) / 4 * z3 + mp.clcos(3, pi - 2 * t)

    # M0: evaluators against nested quadrature
    worst = mp.mpf(0)
    for t in (mp.mpf("0.3"), mp.mpf("0.7"), mp.mpf("1.1"), mp.mpf("1.4"), mp.mpf("1.55")):
        c_q = 2 * mp.quad(lambda s: mp.log(2 * mp.cos(s)), [0, t])
        d_q = 4 * mp.quad(lambda s: (t - s) * mp.log(2 * mp.cos(s)), [0, t])
        worst = max(worst, abs(c_q - Cf(t)), abs(d_q - D3f(t)))
    check("M0", worst < mp.mpf(10) ** -40, "C = Cl2(pi - 2 th) and D3 = (3/4) zeta(3) + Cl3(pi - 2 th) against C = 2 int lt, D3 = 4 int (th - t) lt at 5 angles: worst %s (bar 1e-40)" % dd(worst, 3))
    J2 = half * (q**2 + pi**2 / 12)
    # M1
    t = time.time()
    M1 = mp.quad(lambda th: mp.log(mp.cos(th)) ** 2 * D3f(th), [0, half])
    zz = 5 if MUT != "m1_basis" else 4
    J2c = J2 + (mp.mpf("0.01") if MUT == "m1_J2" else 0)
    M1c = mp.mpf(3) / 4 * z3 * J2c + pi / 4 * (2 * q * z4 + zz * z5 - 2 * z2 * z3)
    check("M1", digits(M1, M1c) >= 40, "int l^2 D3 = %s by quadrature (%.0f s) against 3/4 zeta(3) J2 + (pi/4)(2 ln2 zeta(4) + 5 zeta(5) - 2 zeta(2) zeta(3)) = %s: %s digits (bar 40)"
          % (mp.nstr(M1, 30), time.time() - t, mp.nstr(M1c, 30), dd(digits(M1, M1c), 4)))
    t = time.time()
    M2 = mp.quad(lambda th: mp.log(2 * mp.cos(th)) ** 2 * Cf(th) ** 2, [0, half])
    M3 = mp.quad(lambda th: mp.log(2 * mp.cos(th)) ** 3 * D3f(th), [0, half])
    print("# M2 = <lt^2 C^2> = %s ; M3 = <lt^3 D3> = %s  (%.0f s)" % (mp.nstr(M2, 60), mp.nstr(M3, 60), time.time() - t))
    # M4: PSLQ in my own bases
    mp.mp.dps = 70
    b1 = [M1 / pi, z5, pi**2 * z3, pi**4 * q, z3 * q**2, pi**2 * q**3, q**5]
    r1 = mp.pslq(b1, maxcoeff=10**5, maxsteps=10**6)
    print("# M4 M1/pi in {1, zeta(5), pi^2 zeta(3), pi^4 ln2, zeta(3) ln^2 2, pi^2 ln^3 2, ln^5 2}: %s" % (r1,))
    b2 = lambda M: [M / pi, pi**6, z3**2, pi**2 * z3, pi**4, z5, z3, pi**2, mp.mpf(1)]
    r2 = mp.pslq(b2(M2), maxcoeff=10**6, maxsteps=10**6)
    r3 = mp.pslq(b2(M3), maxcoeff=10**6, maxsteps=10**6)
    print("# M4 M2/pi in {pi^6, zeta(3)^2, pi^2 zeta(3), pi^4, zeta(5), zeta(3), pi^2, 1}: %s" % (r2,))
    print("# M4 M3/pi in the same basis: %s" % (r3,))
    check("M4", r1 is not None and r2 is not None and r3 is not None, "PSLQ found a relation for each of M1, M2, M3 in the bases above")
    # the identified forms of M2 and M3, tested at 120 digits (not the precision they were found at)
    mp.mp.dps = 125
    q, pi = mp.log(2), mp.pi
    z3 = mp.zeta(3)
    Cf2 = lambda t: mp.clsin(2, pi - 2 * t)
    D3f2 = lambda t: mp.mpf(3) / 4 * mp.zeta(3) + mp.clcos(3, pi - 2 * t)
    t = time.time()
    M2h = mp.quad(lambda th: mp.log(2 * mp.cos(th)) ** 2 * Cf2(th) ** 2, [0, pi / 2])
    M3h = mp.quad(lambda th: mp.log(2 * mp.cos(th)) ** 3 * D3f2(th), [0, pi / 2])
    M2c = pi * (mp.mpf(3) / 4 * z3**2 - mp.mpf(19) / 18144 * pi**6)
    M3c = pi * (mp.mpf(9) / 16 * z3**2 - mp.mpf(11) / 3780 * pi**6)
    check("M2", digits(M2h, M2c) >= 40, "<lt^2 C^2> = %s by quadrature at 125 digits (%.0f s) against pi (3/4 zeta(3)^2 - 19 pi^6/18144) found by PSLQ: %s digits"
          % (mp.nstr(M2h, 30), time.time() - t, dd(digits(M2h, M2c), 4)))
    check("M3", digits(M3h, M3c) >= 40, "<lt^3 D3> = %s against pi (9/16 zeta(3)^2 - 11 pi^6/3780) found by PSLQ: %s digits" % (mp.nstr(M3h, 30), dd(digits(M3h, M3c), 4)))


def run_masters2(args):
    """every master integral that omega4.md and omega5.md use, by quadrature of the definitions (no Fourier series, no Parseval)
    against the closed forms the notes state (read after the derivations were committed; amendment A2)."""
    mp.mp.dps = 58
    q, pi = mp.log(2), mp.pi
    z2, z3, z4, z5, z6 = (mp.zeta(n) for n in (2, 3, 4, 5, 6))
    half = pi / 2
    ell = lambda t: mp.log(mp.cos(t))
    lt = lambda t: mp.log(2 * mp.cos(t))
    Cf = lambda t: mp.clsin(2, pi - 2 * t)
    D3f = lambda t: mp.mpf(3) / 4 * z3 + mp.clcos(3, pi - 2 * t)
    cacheK2, cacheK3 = {}, {}
    def K2(t):
        if t not in cacheK2:
            cacheK2[t] = mp.quad(lambda s_: lt(s_) ** 2, [0, t])
        return cacheK2[t]
    def K3(t):
        if t not in cacheK3:
            cacheK3[t] = mp.quad(lambda s_: lt(s_) ** 3, [0, t])
        return cacheK3[t]
    def Q(f):
        return mp.quad(f, [0, half])
    def cmp(name, val, ref, what, bar=40):
        d = digits(val, ref) if ref != 0 else -mp.log10(abs(val))
        dtxt = dd(d, 4) if d < 900 else ">= %d (equal at the working precision)" % mp.mp.dps
        check(name, d >= bar, "%s: quadrature %s, closed form %s, %s digits (bar %d)" % (what, mp.nstr(val, 22), mp.nstr(ref, 22), dtxt, bar))
    # moments of ln cos and of lt
    J = [Q(lambda t, n=n: ell(t) ** n) for n in range(0, 6)]
    Bf = lambda a: mp.sqrt(pi) / 2 * mp.gamma((a + 1) / 2) / mp.gamma(a / 2 + 1)
    tay = mp.taylor(Bf, 0, 5)
    worst = max(digits(J[n], tay[n] * mp.factorial(n)) for n in range(1, 6))
    check("MJ", min(digits(J[n], tay[n] * mp.factorial(n)) for n in range(1, 6)) >= 40, "J_n = int ln^n cos for n = 1 .. 5 against n! [a^n] B(a), B(a) = (sqrt(pi)/2) Gamma((a+1)/2)/Gamma(a/2+1): worst %s digits" % dd(min(digits(J[n], tay[n] * mp.factorial(n)) for n in range(1, 6)), 4))
    c1219 = 1080 if MUT != "m2_J5" else 1008
    J5c = pi / 96 * (-19 * pi**4 * q - c1219 * z5 - 60 * pi**2 * z3 - 720 * q**2 * z3 - 40 * pi**2 * q**3 - 48 * q**5)
    cmp("MJ5", J[5], J5c, "J5 of omega4.md")
    Jt = [Q(lambda t, n=n: lt(t) ** n) for n in range(0, 7)]
    Jtc = {2: pi**3 / 24, 3: -3 * pi * z3 / 4, 4: 19 * pi**5 / 480, 5: -(5 * pi / 8) * (18 * z5 + pi**2 * z3),
           6: 5 * pi / 2688 * (6048 * z3**2 + 55 * pi**6)}
    for n in range(2, 7):
        cmp("MJt%d" % n, Jt[n], Jtc[n], "J~_%d = int ln^%d(2 cos)" % (n, n))
    cmp("Mth2", Q(lambda t: t**2), pi**3 / 24, "<theta^2>")
    C2 = Q(lambda t: Cf(t) ** 2)
    cmp("MC2", C2, pi**5 / 360, "<C^2>")
    cmp("MlC2", Q(lambda t: ell(t) * Cf(t) ** 2), -q * pi**5 / 360, "<ell C^2> = -ln2 <C^2>")
    v = Q(lambda t: lt(t) * Cf(t) ** 2)
    check("MltC2", abs(v) < mp.mpf(10) ** -40, "<lt C^2> = %s (closed form 0, bar 1e-40)" % mp.nstr(v, 5))
    cmp("MthltC", Q(lambda t: t * lt(t) * Cf(t)), -C2 / 4, "<theta lt C> = -<C^2>/4")
    cmp("MthC", Q(lambda t: t * Cf(t)), pi * z3 / 4, "<theta C>")
    cmp("MuC", Q(lambda t: mp.tan(t) * Cf(t)), pi**3 / 12, "<u C>")
    l2D3 = Q(lambda t: lt(t) ** 2 * D3f(t))
    cmp("Mlt2D3", l2D3, pi / 4 * (5 * z5 - 2 * z2 * z3) + pi**3 * z3 / 32, "<lt^2 D3>")
    ell2D3 = Q(lambda t: ell(t) ** 2 * D3f(t))
    Lc2 = lambda t: K2(t) - q * Cf(t) + q**2 * t        # Lc2 = K2 - ln2 C + ln^2 2 theta
    cmp("MCLc2", Q(lambda t: Cf(t) * Lc2(t)), mp.mpf(7) / 8 * z3 * J[2] - ell2D3 / 2, "<C Lc2> = (7/8) zeta(3) J2 - (1/2) int ell^2 D3")
    cmp("MCK2", Q(lambda t: Cf(t) * K2(t)), mp.mpf(7) / 8 * z3 * Jt[2] - l2D3 / 2, "<C K2> = (7/8) zeta(3) J~2 - <lt^2 D3>/2")
    lt2C2 = Q(lambda t: lt(t) ** 2 * Cf(t) ** 2)
    lt3D3 = Q(lambda t: lt(t) ** 3 * D3f(t))
    cmp("MltCK2", Q(lambda t: lt(t) * Cf(t) * K2(t)), -lt2C2 / 4, "<lt C K2> = -<lt^2 C^2>/4")
    cmp("MCK3", Q(lambda t: Cf(t) * K3(t)), mp.mpf(7) / 8 * z3 * Jt[3] - lt3D3 / 2, "<C K3> = (7/8) zeta(3) J~3 - <lt^3 D3>/2")
    th2lt2 = Q(lambda t: t**2 * lt(t) ** 2)
    cmp("Mth2lt2", th2lt2, 11 * pi**5 / 1440, "<theta^2 lt^2>")
    cmp("MthK2", Q(lambda t: t * K2(t)), pi**2 / 8 * Jt[2] - th2lt2 / 2, "<theta K2> = (pi^2/8) J~2 - <theta^2 lt^2>/2")
    c61 = mp.mpf(61) / 6 if MUT != "m2_K2K2" else mp.mpf(60) / 6
    cmp("MK2K2", Q(lambda t: K2(t) ** 2), pi**7 / 3456 - pi**7 / 2880 + pi / 16 * (c61 * z6 - 6 * z3**2), "<K2^2>")
    cmp("NM2", lt2C2, pi * (mp.mpf(3) / 4 * z3**2 - mp.mpf(95) / 96 * z6), "<lt^2 C^2> against omega5.md's pi (3/4 zeta(3)^2 - (95/96) zeta(6))")
    c114 = mp.mpf(11) / 4 if MUT != "m2_M3" else mp.mpf(10) / 4
    cmp("NM3", lt3D3, pi * (mp.mpf(9) / 16 * z3**2 - c114 * z6), "<lt^3 D3> against omega5.md's pi (9/16 zeta(3)^2 - (11/4) zeta(6))")


def run_fp(args):
    mp.mp.dps = 60
    q, pi = mp.log(2), mp.pi
    half = pi / 2
    J1 = -half * q
    J2 = half * (q**2 + pi**2 / 12)
    # rule (derivations.md section 9):  FP int u = 0 (FP of -ln sin x), FP int u th = J1,  FP int u Lc1 = J2
    fpu = mp.mpf(0) if MUT != "fp_const" else -q
    def Cf(t):
        return mp.clsin(2, pi - 2 * t)
    I1 = mp.quad(lambda t: mp.tan(t) * (half - t), [0, half])
    R1 = half * fpu - J1
    check("FP1", digits(I1, R1) >= 40, "int u (pi/2 - th) = %s by quadrature; the rule gives (pi/2) FP int u - FP int u th = %s: %s digits" % (mp.nstr(I1, 25), mp.nstr(R1, 25), dd(digits(I1, R1), 3)))
    I2 = mp.quad(lambda t: mp.tan(t) * (J1 - (Cf(t) - 2 * t * q) / 2), [0, half])
    R2 = J1 * fpu - J2
    check("FP2", digits(I2, R2) >= 40, "int u (J1 - Lc1) = %s by quadrature; the rule gives J1 FP int u - FP int u Lc1 = %s: %s digits" % (mp.nstr(I2, 25), mp.nstr(R2, 25), dd(digits(I2, R2), 3)))


def run_shifts(args):
    """E-checks again with one constant moved off its root: Om4 + 1/10 (G6 = -0.4 pi), Om5 + 1/10 (G7 = -0.4 pi), Om2 + 1e-3 (HYP)."""
    import copy
    for sh in ("4:0.1", "5:0.1", "2:0.001"):
        a2 = copy.copy(args)
        a2.shift = sh
        print("\n-- shift %s --" % sh)
        run_exch(a2)


# declared in the raw README (amendment A1) from development runs dev_*; each control runs the parts listed in CTL_PARTS
CTL_PARTS = {
    "e4_20": "sym exch", "e4_Om1": "sym exch", "e4_C": "sym exch", "e2w_P": "sym exch", "e2w_g": "sym exch",
    "e5_noD3": "sym exch", "e5_a1": "sym exch", "e5_K": "sym exch", "hyp_G4": "exch", "c_wrong": "exch bnd",
    "g6r_slope": "exch", "g7r_slope": "exch", "m1_basis": "masters", "m1_J2": "masters", "fp_const": "fp",
    "m2_K2K2": "masters2", "m2_J5": "masters2", "m2_M3": "masters2",
}
TARGETS = {   # name -> the checks that must fail (and no other), from development runs dev_ctl_* / dev2_ctl_* of raw I264-rpi5 (amendment A1)
    "e4_20": ['E4_6', 'E4_7', 'G6R', 'G7R', 'SY_W4'],
    "e4_Om1": ['E4_6', 'E4_7', 'G6R', 'G7R', 'SY_W4', 'SY_bnd0'],
    "e4_C": ['E4_6', 'E4_7', 'G6R', 'G7R', 'SY_W4'],
    "e2w_P": ['E2Wb', 'G7R', 'SY_P', 'SY_bnd0'],
    "e2w_g": ['E2Wb', 'SY_P'],
    "e5_noD3": ['E5', 'G7R', 'SY_W'],
    "e5_a1": ['E5', 'G7R', 'SY_W'],
    "e5_K": ['E5', 'G7R', 'SY_K'],
    "hyp_G4": ['E1_6', 'E1_7', 'E2Wa', 'E2Wb', 'E2_4', 'E2_5', 'E2_6', 'E3_5', 'E3_6', 'E3_7', 'E4_7', 'G6R', 'G7R', 'H2_G5', 'H2_G6'],
    "c_wrong": ['BND_outer', 'E1_6', 'E1_7', 'E2Wa', 'E2Wb', 'E2_3', 'E2_4', 'E2_5', 'E2_6', 'E3_4', 'E3_5', 'E3_6', 'E3_7', 'E4_6', 'E4_7', 'E5', 'G6R', 'G7R', 'H1', 'H2_G3', 'H2_G4', 'H2_G5', 'H2_G6'],
    "g6r_slope": ['G6R'],
    "g7r_slope": ['G7R'],
    "m1_basis": ['M1'],
    "m1_J2": ['M1'],
    "fp_const": ['FP1', 'FP2'],
    "m2_K2K2": ['MK2K2'],
    "m2_J5": ['MJ5'],
    "m2_M3": ['NM3'],
}


def run_ctlsum(args):
    """read ctl_<name>.log in --ctldir: the set of FAIL lines against TARGETS."""
    on, tot = 0, 0
    for name in sorted(TARGETS):
        path = os.path.join(args.ctldir, "%s_%s.log" % (args.ctlprefix, name))
        rc = os.path.join(args.ctldir, "%s_%s.rc" % (args.ctlprefix, name))
        if not os.path.exists(path):
            check("CT_" + name, False, "no log " + path)
            continue
        fired = sorted({ln.split()[1] for ln in open(path, encoding="utf-8") if ln.startswith("FAIL ")})
        want = sorted(TARGETS[name])
        code = open(rc).read().strip() if os.path.exists(rc) else "?"
        ok = fired == want and code == "1"
        tot += 1
        on += ok
        check("CT_" + name, ok, "fired %s, declared %s, exit %s" % (",".join(fired) or "-", ",".join(want) or "-", code))
    print("# controls on target: %d of %d" % (on, tot))



if __name__ == "__main__":
    main()

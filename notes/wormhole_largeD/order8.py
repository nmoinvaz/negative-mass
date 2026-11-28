"""Order 8 of the matched expansion (matched.py, order3.py ... order7.py): Omega6 from G8 = 0 by method B, the structure of F, and the checks of issue #282.

W^2 a^2 = s + Omega, s = 2(D - 3), eps = 1/s, Omega = Omega0 + Omega1 eps + Omega2 eps^2 + ..., R = r^s, l = ln R = -2 ln cos(theta).
The order equations (E0 derives them from the master equation), the same for every k:
  dY_k/dl = F_{k-1},  F_{k-1} = sum_{j<k} L_{k-1-j} l^j / j!,   (R - 1) L_k' - L_k/2 = h_k,
  h_k = (1 - 1/R) L_{k-1} + sum_{j<k} Y_j c_{k-j} + Y_k/2,   L_k = -2 Omega_{k-1} + u C_k,   G_k = C_k(inf) = 0,
with u = sqrt(R - 1). G_k is affine in Omega_{k-2} with slope -4 pi: G2 = 0 gives Omega0, ..., G7 = 0 Omega5, **G8 = 0 gives Omega6**.
One path computes them all (`omega_chain`); Omega0 .. Omega5 are planted (their closed forms of #26, #208, #220, #263 are the truth).

Method B (order4.py's PanelB, shared): v = sqrt(l), panels of width w on [0, V], N Chebyshev nodes per panel, cumulative integrals by the
spectral integration matrix; Y_k from the left, the tail C_right from the right, L_k = -2 h_k(1) - u C_right, G_k = C_right(0). It uses no
closed form of Y_k or L_k and no input value of Omega_0 .. Omega_5. Levels below k - 1 do not depend on Omega_{k-2}, so the chain keeps them
(`cached`; the same operations on the same numbers as the uncached chain, which CACHE compares).

The structure: W^2 a^2 = s 2^(-4/s) F(1/s), F(eps) = exp(4 ln 2 eps) (1 + Omega0 eps + Omega1 eps^2 + ...) = 1 + f2 eps^2 + ... (#263: free of ln 2 through eps^6).
f7 = Omega6 + sum_{j<=5} Omega_j (4 ln 2)^(6-j)/(6-j)! + (4 ln 2)^7/7! (Omega0 .. Omega5 exact), and its ln 2-free basis is PSLQ's target.

    python3 order8.py grid --grid=N,w,V,dps [--kmax=8]      (a development run of one grid; no gate)

Parts (the arguments of a run; the checks of omega6.md and of the raw README of #282 in parentheses):
    sym      E0, the table against the epsilon-expansion of the master equation, k <= 8
    cache    CACHE, the cached chain against the uncached chain on a small grid
    struct   LN2: F's coefficients f2 .. f6 free of ln 2 (sympy), from the closed forms of Omega0 .. Omega5
    evals    the grids (--grids=G1,G2,B1,B2,B3 or --cgrids), --save=FILE.json
    values   PL0 .. PL5, SL, RES, M6, RNG (--load=FILE.json,...)
    pslq     PQ(basis) and the real PSLQ of f7 in the bases F1 .. F6 (ln 2-free) and X3, X0 (with ln 2) (--bases=F6,X0 --no-real --pslq-cap=D --save-rel=FILE)
    remainder R7, R8, X4 against #195's W^2 (--w2raw=inputs_w2_growthD.txt)
    closed6  ID6, ZL: the closed form of Omega6 from a PSLQ relation (--rel=FILE) against the grids, and zeros on the ln 2 monomials
Mutations: --mutate=NAME (MUTATIONS below). Exit status 1 when a check fails.
"""
import json
import math
import os
import random
import re
import sys
import time
import types

import mpmath as mp
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_argv, sys.argv = sys.argv, sys.argv[:1]          # order4.py parses sys.argv when it is imported
from order4 import PanelB  # noqa: E402
import order4 as o4  # noqa: E402
sys.argv = _argv

MUTATIONS = ("wrong_h", "wrong_c6", "wrong_c7", "wrong_c8", "wrong_Om0", "wrong_Om1", "wrong_Om2", "wrong_Om3", "wrong_Om4", "wrong_Om5", "coarse_B",
             "short_tail", "stale_cache", "pslq_wrong_basis", "pslq_wrong_target", "wrong_Om6_series", "wrong_closed5", "ln2_leak", "f7_rel", "ln2_struct", "ln2_zl")
MUTATE = None
W2FILE = "w2_output.txt"
SAVE = None
LOAD = None
GRIDARG = None
GRIDSEL = None
CGRIDS_ON = False
W2RAW = None
KMAX = 8
PARTS = []
REAL = True
BASES_SEL = None
PSLQ_CAP = None
REL = None
SAVE_REL = None


def parse(argv):
    global MUTATE, W2FILE, W2RAW, SAVE, LOAD, GRIDARG, GRIDSEL, CGRIDS_ON, KMAX, PARTS, REAL, BASES_SEL, PSLQ_CAP, REL, SAVE_REL
    parts = []
    for a in argv:
        if a.startswith("--mutate="):
            MUTATE = a.split("=", 1)[1]
            assert MUTATE in MUTATIONS, f"unknown mutation {MUTATE}; known: {MUTATIONS}"
        elif a.startswith("--w2="):
            W2FILE = a.split("=", 1)[1]
        elif a.startswith("--w2raw="):
            W2RAW = a.split("=", 1)[1]
        elif a.startswith("--rel="):
            REL = a.split("=", 1)[1]
        elif a.startswith("--save-rel="):
            SAVE_REL = a.split("=", 1)[1]
        elif a.startswith("--save="):
            SAVE = a.split("=", 1)[1]
        elif a.startswith("--load="):
            LOAD = a.split("=", 1)[1]
        elif a.startswith("--grid="):
            GRIDARG = a.split("=", 1)[1]
        elif a.startswith("--grids="):
            GRIDSEL = a.split("=", 1)[1].split(",")
        elif a.startswith("--bases="):
            BASES_SEL = a.split("=", 1)[1].split(",")
        elif a.startswith("--pslq-cap="):
            PSLQ_CAP = int(a.split("=", 1)[1])
        elif a == "--cgrids":
            CGRIDS_ON = True
        elif a == "--no-real":
            REAL = False
        elif a.startswith("--kmax="):
            KMAX = int(a.split("=", 1)[1])
        else:
            parts.append(a)
    PARTS = parts or ["sym", "struct", "evals", "values", "pslq", "remainder"]


RESULTS = []


def report(cid, ok, text):
    RESULTS.append((cid, ok))
    print(f"{'PASS' if ok else 'FAIL'}    {cid}: {text}", flush=True)


def digits(a, b):
    """-log10 of the relative difference, capped at the working precision in digits"""
    a, b = mp.mpf(a), mp.mpf(b)
    d = abs(a - b)
    return mp.mpf(mp.mp.dps) if d == 0 else min(mp.mpf(mp.mp.dps), -mp.log10(d / abs(a)))


def nd(x, n=40):
    return mp.nstr(x, n)


# ---------------------------------------------------------------------------------------------------------------
# the table of c_m and h_m (integer and rational arithmetic only, so that it also runs on sympy symbols)
# ---------------------------------------------------------------------------------------------------------------
def cc(m, Om, l, R):
    if m == 1:
        return Om[0] + 3 * l / 2 - 1 + 1 / R
    if m == 2:
        return Om[0] * l + Om[1] + l**2 / 4 + l + (l / R if MUTATE == "wrong_h" else -l / R)
    if m == 3:
        return Om[0] * l**2 / 2 + Om[1] * l + Om[2] + l**3 / 4 - l**2 / 2 + l**2 / (2 * R)
    if m == 4:
        return Om[0] * l**3 / 6 + Om[1] * l**2 / 2 + Om[2] * l + Om[3] + l**4 / 48 + l**3 / 6 - l**3 / (6 * R)
    if m == 5:
        return (Om[0] * l**4 / 24 + Om[1] * l**3 / 6 + Om[2] * l**2 / 2 + Om[3] * l + Om[4] + l**5 / 80 - l**4 / 24
                + l**4 / (24 * R))
    if m == 6:
        return (Om[0] * l**5 / 120 + Om[1] * l**4 / 24 + Om[2] * l**3 / 6 + Om[3] * l**2 / 2 + Om[4] * l + Om[5]
                + l**6 / 1440 + l**5 / 120 + (1 if MUTATE == "wrong_c6" else -1) * l**5 / (120 * R))
    if m == 7:           # new at this order: the coefficient of eps^7 of (1 + eps Om) e^{eps l} - (1/2 + eps)(1 - 1/R) e^{-eps l} - e^{-eps l}/(2R)
        return (Om[0] * l**6 / 720 + Om[1] * l**5 / 120 + Om[2] * l**4 / 24 + Om[3] * l**3 / 6 + Om[4] * l**2 / 2 + Om[5] * l + Om[6]
                + l**7 / 3360 - l**6 / 720 + (-1 if MUTATE == "wrong_c7" else 1) * l**6 / (720 * R))
    if m == 8:           # new at this order: the coefficient of eps^8 of the same bracket: Om_j l^(7-j)/(7-j)! (j = 0 .. 7), l^8/80640 + l^7/5040 - l^7/(5040 R)
        return (Om[0] * l**7 / 5040 + Om[1] * l**6 / 720 + Om[2] * l**5 / 120 + Om[3] * l**4 / 24 + Om[4] * l**3 / 6 + Om[5] * l**2 / 2 + Om[6] * l + Om[7]
                + l**8 / 80640 + l**7 / 5040 + (1 if MUTATE == "wrong_c8" else -1) * l**7 / (5040 * R))
    raise ValueError(m)


def hk(m, Om, l, R, Y, L):
    """h_m"""
    if m == 0:
        return Y[0] / 2
    return (1 - 1 / R) * L[m - 1] + sum(Y[j] * cc(m - j, Om, l, R) for j in range(m)) + Y[m] / 2


def Fsum(k, L, l):
    return sum(L[k - 1 - j] * l**j / math.factorial(j) for j in range(k))


# ---------------------------------------------------------------------------------------------------------------
# E0: the epsilon-expansion of the master equation against the table, k = 0 .. 8
# ---------------------------------------------------------------------------------------------------------------
def check_sym():
    t0 = time.time()
    R, eps, l = sp.symbols("R epsilon l", positive=True)
    Om = sp.symbols("Omega0:8")
    Ys, Ls = sp.symbols("Y0:9"), sp.symbols("L0:9")
    kk = 8
    r = sum((eps * l) ** j / sp.factorial(j) for j in range(kk + 1))          # exact through eps^8, which is all that is read
    rinv = sum((-eps * l) ** j / sp.factorial(j) for j in range(kk + 1))
    Y = sum(Ys[k] * eps**k for k in range(kk + 1))
    L = sum(Ls[k] * eps**k for k in range(kk + 1))
    O = sum(Om[k] * eps**k for k in range(kk))
    half = sp.Rational(1, 2)
    rhs = ((1 + eps * O) * r * Y + L / (2 * R) + (half + eps) * L * (R - 1) / R - (half + eps) * Y * (R - 1) / R * rinv
           - Y / (2 * R) * rinv)
    ser = sp.expand(rhs)
    yser = sp.expand(eps * L * r)
    res = []
    for k in range(kk + 1):
        res.append(sp.simplify(ser.coeff(eps, k) - Ls[k] / 2 - hk(k, Om, l, R, Ys, Ls)))
        if k >= 1:
            res.append(sp.simplify(yser.coeff(eps, k) - Fsum(k, Ls, l)))
    report("E0", all(x == 0 for x in res),
           f"h_k (k = 0 .. 8) and dY_k/dl (k = 1 .. 8) from the epsilon-expansion of the master equation against the table: "
           f"residuals {[str(x) for x in res if x != 0] or '0'} ({time.time() - t0:.1f} s)")


# ---------------------------------------------------------------------------------------------------------------
# method B at order k
# ---------------------------------------------------------------------------------------------------------------
class PanelK(PanelB):
    def __init__(self, N, w, V):
        super().__init__(N, w, V)
        self._lp, self._cg, self._Sg = {}, {}, None

    def grid(self, fn):
        return [[fn(p, i) for i in range(self.N)] for p in range(self.P)]

    def lp(self, j):
        """l^j / j! on the nodes"""
        if j not in self._lp:
            self._lp[j] = self.grid(lambda p, i: self.l[p][i] ** j / mp.factorial(j))
        return self._lp[j]

    def cg(self, m, Om):
        """c_m on the nodes (kept while Om[0 .. m - 1] stays the same)"""
        key = tuple(Om[:m])
        hit = self._cg.get(m)
        if hit is None or hit[0] != key:
            hit = (key, self.grid(lambda p, i: cc(m, Om, self.l[p][i], self.R[p][i])))
            self._cg[m] = hit
        return hit[1]

    def Sg(self):
        if self._Sg is None:
            self._Sg = self.grid(lambda p, i: self.Rm1[p][i] / self.R[p][i])                      # 1 - 1/R
        return self._Sg

    def cum_left(self, F):
        N, half = self.N, self.half
        out, acc = [], mp.mpf(0)
        for p in range(self.P):
            out.append([acc + half * mp.fdot(self.S[i], F[p]) for i in range(N)])
            acc += half * mp.fdot(self.wend, F[p])
        return out

    def cum_right(self, F):
        N, half = self.N, self.half
        tot = [half * mp.fdot(self.wend, F[p]) for p in range(self.P)]
        out, acc = [None] * self.P, mp.mpf(0)
        for p in range(self.P - 1, -1, -1):
            out[p] = [acc + tot[p] - half * mp.fdot(self.S[i], F[p]) for i in range(N)]
            acc += tot[p]
        return out, acc

    def yk(self, k, L):
        """Y_k = int_0^l F_{k-1} dl from the left (needs L_0 .. L_{k-1}; k >= 1)"""
        F = self.grid(lambda p, i: 2 * self.v[p][i] * mp.fdot([L[k - 1 - j][p][i] for j in range(k)], [self.lp(j)[p][i] for j in range(k)]))
        return self.cum_left(F)

    def lk(self, k, Om, Y, L, G):
        """append L_k and set G[k], from Y_0 .. Y_k and L_0 .. L_{k-1}; reads Om[0 .. k - 1] (Om[k - 1] is h_k(1), which drops out of G_k)"""
        thr = mp.mpf(1) / 2 if k == 0 else Om[k - 1]
        if k == 0:
            s = self.grid(lambda p, i: Y[0][p][i] / 2 - thr)
        else:
            cg = {m: self.cg(m, Om) for m in range(1, k + 1)}
            Sg = self.Sg()
            s = self.grid(lambda p, i: Sg[p][i] * L[k - 1][p][i] + mp.fdot([Y[j][p][i] for j in range(k)], [cg[k - j][p][i] for j in range(k)])
                          + Y[k][p][i] / 2 - thr)
        Cint = self.grid(lambda p, i: s[p][i] / self.Rm1[p][i] * self.R[p][i] / self.u[p][i] * 2 * self.v[p][i])
        Cr, G[k] = self.cum_right(Cint)
        L.append(self.grid(lambda p, i: -2 * thr - self.u[p][i] * Cr[p][i]))

    def hierarchy(self, Om, kmax):
        """Y_0 .. Y_kmax, L_0 .. L_kmax and G_0 .. G_kmax for the given Omega list, from level 0 (the uncached path)"""
        Y = [self.grid(lambda p, i: mp.mpf(1))]
        L, G = [], {}
        for k in range(kmax + 1):
            if k >= 1:
                Y.append(self.yk(k, L))
            self.lk(k, Om, Y, L, G)
        return Y, L, G


def solve_omega(Gfun, k, Om):
    """Omega_{k-2} from the affine G_k; returns (value, slope)"""
    o0 = list(Om)
    o0[k - 2] = mp.mpf(0)
    o1 = list(Om)
    o1[k - 2] = mp.mpf(1)
    g0, g1 = Gfun(k, o0), Gfun(k, o1)
    return -g0 / (g1 - g0), g1 - g0


def closed_omegas():
    """Omega_0 .. Omega_5 of #26, #208, #220 and #263 (derived there), the planted truth"""
    pi, ln2, z3, z5 = mp.pi, mp.ln2, mp.zeta(3), mp.zeta(5)
    return [-4 * ln2, 2 * pi**2 / 3 + 8 * ln2**2,
            80 * z3 - mp.mpf(8) / 3 * pi**2 * ln2 - mp.mpf(32) / 3 * ln2**3,
            22 * pi**4 / 45 + mp.mpf(16) / 3 * pi**2 * ln2**2 + mp.mpf(32) / 3 * ln2**4 - 320 * z3 * ln2 - 384 * z3,
            -10432 * z5 - 160 * pi**2 * z3 + 640 * ln2**2 * z3 - mp.mpf(64) / 9 * pi**2 * ln2**3 - mp.mpf(88) / 45 * pi**4 * ln2
            - mp.mpf(128) / 15 * ln2**5 - mp.mpf(2432) / 45 * pi**4 + 1536 * ln2 * z3 + 3072 * z3,
            mp.mpf(256) / 45 * ln2**6 + mp.mpf(64) / 9 * pi**2 * ln2**4 - mp.mpf(2560) / 3 * z3 * ln2**3 + mp.mpf(176) / 45 * pi**4 * ln2**2 - 3072 * z3 * ln2**2
            + 640 * pi**2 * z3 * ln2 + mp.mpf(9728) / 45 * pi**4 * ln2 - 12288 * z3 * ln2 + 41728 * z5 * ln2 - mp.mpf(38812) / 315 * pi**6 + 3200 * z3**2
            - mp.mpf(256) / 3 * pi**2 * z3 + 74240 * z5 + mp.mpf(39424) / 45 * pi**4 - 36864 * z3 + (mp.mpf(1) / 1000 if MUTATE == "wrong_closed5" else 0)]


WRONG_OM = ("wrong_Om0", "wrong_Om1", "wrong_Om2", "wrong_Om3", "wrong_Om4", "wrong_Om5")


def omega_chain_ref(B, kmax, log=None):
    """Omega_0 .. Omega_{kmax-2} from G_2 .. G_kmax = 0, in turn, each G_k from the whole hierarchy (the uncached path); returns (Om, slopes, secs)"""
    Om = [mp.mpf(0)] * kmax
    slopes, secs, resid = {}, {}, {}
    cl = closed_omegas()
    for k in range(2, kmax + 1):
        t0 = time.time()
        Om_in = list(Om)
        for i, name in enumerate(WRONG_OM):
            if MUTATE == name and k >= i + 3:
                Om_in[i] = cl[i] + mp.mpf(1) / 100
        val, slope = solve_omega(lambda kk, o: B.hierarchy(o, kk)[2][kk], k, Om_in)
        Om[k - 2] = val
        slopes[k] = slope
        Om_root = list(Om_in)
        Om_root[k - 2] = val
        resid[k] = abs(B.hierarchy(Om_root, k)[2][k]) / abs(slope)           # G_k at the root
        secs[k] = time.time() - t0
        if log:
            print(f"    k = {k}: Omega_{k - 2} = {nd(val, 30)}  slope {nd(slope, 20)}  ({secs[k]:.0f} s)", flush=True)
    return Om, slopes, secs, resid


def omega_chain(B, kmax, log=None):
    """The same chain with the levels below k - 1 kept: Y_j, L_j (j <= k - 2) do not depend on Omega_{k-2}. Per step: Y_{k-1}, two trials of
    L_{k-1}, Y_k, G_k at Omega_{k-2} = 0, 1, and L_{k-1} at the root. Returns (Om, slopes, secs, resid) with resid[k] = G_k at the root / |slope|."""
    if MUTATE in WRONG_OM:
        return omega_chain_ref(B, kmax, log=log)           # a wrong input must reach every level: recompute from level 0
    Om = [mp.mpf(0)] * kmax
    slopes, secs, resid = {}, {}, {}
    Y = [B.grid(lambda p, i: mp.mpf(1))]
    L, G = [], {}
    B.lk(0, Om, Y, L, G)
    stale = None
    for k in range(2, kmax + 1):
        t0 = time.time()
        Y.append(B.yk(k - 1, L))                                  # Y_{k-1}: from L_0 .. L_{k-2}, independent of Omega_{k-2}

        def trial(o):
            Omt = list(Om)
            Omt[k - 2] = o
            Yt, Lt, Gt = list(Y), list(L), {}
            B.lk(k - 1, Omt, Yt, Lt, Gt)
            Yt.append(B.yk(k, Lt))
            B.lk(k, Omt, Yt, Lt, Gt)
            return Gt[k], Lt[k - 1], Gt[k - 1]
        g0, l0, gm0 = trial(mp.mpf(0))
        g1 = trial(mp.mpf(1))[0]
        val, slope = -g0 / (g1 - g0), g1 - g0
        Om[k - 2] = val
        slopes[k] = slope
        if MUTATE == "stale_cache" and k >= 3:
            L.append(l0)                                          # level k - 1 from the first trial, not recomputed at the root
            G[k - 1] = gm0
        else:
            B.lk(k - 1, Om, Y, L, G)                              # level k - 1 at the root (G_{k-1} at the root is read below)
        if k >= 3:
            resid[k - 1] = abs(G[k - 1]) / abs(slopes[k - 1])
        if k == kmax:                                             # the last order: G_kmax at the root, through one more pair of levels
            Yt, Lt, Gt = list(Y), list(L), {}
            Yt.append(B.yk(k, Lt))
            B.lk(k, Om, Yt, Lt, Gt)
            resid[k] = abs(Gt[k]) / abs(slope)
        secs[k] = time.time() - t0
        if log:
            print(f"    k = {k}: Omega_{k - 2} = {nd(val, 30)}  slope {nd(slope, 20)}  ({secs[k]:.0f} s)", flush=True)
    return Om, slopes, secs, resid


def run_B(N, w, V, dps, kmax=8, log=False, cached=True):
    t0 = time.time()
    if MUTATE == "coarse_B":
        N = 10
    if MUTATE == "short_tail":
        V = 6
    with mp.workdps(dps):
        w = mp.mpf(w[0]) / w[1] if isinstance(w, tuple) else mp.mpf(w)
        tinit = time.time()
        B = PanelK(N, w, V)
        tinit = time.time() - tinit
        Om, slopes, secs, resid = (omega_chain if cached else omega_chain_ref)(B, kmax, log=log)
        out = {f"Om{i}": +Om[i] for i in range(kmax - 1)}
        out.update({f"slope{k}": +s for k, s in slopes.items()})
        out.update({f"resid{k}": +r for k, r in resid.items()})
        out["secs"] = time.time() - t0
        out["init_secs"] = tinit
        out["grid"] = f"N={N} w={mp.nstr(w, 6)} V={V} dps={dps} panels={B.P}"
    return out


# declared grids and the control grids are fixed in the raw README before the runs of record (see there)
GRIDS = {"G1": (300, (1, 4), 50, 465), "G2": (310, (3, 10), 49, 455)}             # w as a fraction (num, den), evaluated at the working precision
REPORT_GRIDS = {"B1": (60, (2, 5), 26, 100), "B2": (100, (1, 4), 32, 170), "B3": (170, (1, 4), 40, 285)}
CGRIDS = {"G1": (90, (1, 4), 32, 160), "G2": (85, (3, 10), 31, 150)}


def save_ev(ev, path):
    with open(path, "w") as f:
        json.dump({n: {k: (mp.nstr(v, 600) if isinstance(v, mp.mpf) else v) for k, v in e.items()} for n, e in ev.items()}, f, indent=1)


def load_ev(paths):
    out = {}
    for path in paths.split(","):
        raw = json.load(open(path))
        out.update({n: {k: (mp.mpf(v) if k.startswith(("Om", "slope", "resid")) else v) for k, v in e.items()} for n, e in raw.items()})
    return out


def check_grid():
    """one grid, no gate: the development run"""
    N, w, V, dps = GRIDARG.split(",")
    mp.mp.dps = 60
    e = run_B(int(N), tuple(int(t) for t in w.split("/")) if "/" in w else w, int(V), int(dps), kmax=KMAX, log=True)
    mp.mp.dps = int(dps)
    cl = closed_omegas()
    print(f"  {e['grid']}  {e['secs']:.0f} s (panel setup {e['init_secs']:.0f} s)")
    for i in range(KMAX - 1):
        print(f"    Omega{i} = {nd(e[f'Om{i}'], int(dps))}")
    for i in range(min(6, KMAX - 1)):
        print(f"    Omega{i} against closed form: {float(digits(cl[i], e[f'Om{i}'])):.1f} digits")
    for k in range(2, KMAX + 1):
        print(f"    slope G{k} against -4 pi: {float(digits(-4 * mp.pi, e[f'slope{k}'])):.1f} digits")
    for k, r in sorted((int(n[5:]), e[n]) for n in e if n.startswith("resid")):
        print(f"    |G{k}| at the root / |slope| = {nd(r, 3)}")
    return e


def check_cache():
    """CACHE: the cached chain against the uncached chain, on a small grid: the same Omegas"""
    N, w, V, dps = 24, (1, 2), 12, 60
    mp.mp.dps = 40
    a = run_B(N, w, V, dps, kmax=8)
    b = run_B(N, w, V, dps, kmax=8, cached=False)
    mp.mp.dps = dps
    worst = max(abs(a[f"Om{i}"] - b[f"Om{i}"]) / abs(b[f"Om{i}"]) for i in range(7))
    dig = [float(digits(b[f"Om{i}"], a[f"Om{i}"])) for i in range(7)]
    report("CACHE", worst <= mp.mpf(10) ** -(dps - 5),
           f"cached chain against the uncached chain on (N {N}, w 1/2, V {V}, {dps} digits): Omega0 .. Omega6 equal to "
           + ", ".join(f"{d:.1f}" for d in dig) + f" digits, worst relative difference {nd(worst, 3)} (bar 1e-{dps - 5}); times {a['secs']:.0f} s and {b['secs']:.0f} s")


# ---------------------------------------------------------------------------------------------------------------
# the structure of F: f_n = [eps^n] exp(4 ln 2 eps)(1 + Omega0 eps + Omega1 eps^2 + ...)
# ---------------------------------------------------------------------------------------------------------------
def f_coeff(n, a):
    """f_n = sum_{j<=n} a_j (4 ln 2)^(n-j)/(n-j)!, a_0 = 1, a_j = Omega_{j-1}"""
    q4 = 4 * mp.ln2
    return sum(a[j] * q4 ** (n - j) / mp.factorial(n - j) for j in range(n + 1))


def f7_of(om6):
    """f_7 from Omega6 and the exact Omega0 .. Omega5 (the closed forms)"""
    return f_coeff(7, [mp.mpf(1)] + closed_omegas() + [om6])


def check_struct():
    """LN2: f_2 .. f_6 free of ln 2, symbolically, from the closed forms of Omega0 .. Omega5 (d/d ln 2 of each coefficient is 0)"""
    t0 = time.time()
    q, z3, z5, p2 = sp.symbols("q z3 z5 p2")                      # ln 2, zeta(3), zeta(5), pi^2
    R = sp.Rational
    Om = [-4 * q, 2 * p2 / 3 + 8 * q**2,
          80 * z3 - R(8, 3) * p2 * q - R(32, 3) * q**3,
          22 * p2**2 / 45 + R(16, 3) * p2 * q**2 + R(32, 3) * q**4 - 320 * z3 * q - 384 * z3,
          -10432 * z5 - 160 * p2 * z3 + 640 * q**2 * z3 - R(64, 9) * p2 * q**3 - R(88, 45) * p2**2 * q - R(128, 15) * q**5 - R(2432, 45) * p2**2 + 1536 * q * z3 + 3072 * z3,
          R(256, 45) * q**6 + R(64, 9) * p2 * q**4 - R(2560, 3) * z3 * q**3 + R(176, 45) * p2**2 * q**2 - 3072 * z3 * q**2 + 640 * p2 * z3 * q + R(9728, 45) * p2**2 * q
          - 12288 * z3 * q + 41728 * z5 * q - R(38812, 315) * p2**3 + 3200 * z3**2 - R(256, 3) * p2 * z3 + 74240 * z5 + R(39424, 45) * p2**2 - 36864 * z3]
    if MUTATE == "ln2_struct":
        Om[3] = Om[3] + q * z3                                    # a stray ln 2 term in Omega3
    a = [sp.Integer(1)] + Om
    fs = [sp.expand(sum(a[j] * (4 * q) ** (n - j) / sp.factorial(n - j) for j in range(n + 1))) for n in range(7)]
    dq = [sp.expand(sp.diff(f, q)) for f in fs]
    report("LN2", all(x == 0 for x in dq) and fs[0] == 1 and fs[1] == 0,
           "F = exp(4 ln 2 eps)(1 + Omega0 eps + ... + Omega5 eps^6) from the closed forms: f0 = 1, f1 = 0, and d f_n / d ln 2 = 0 for n = 2 .. 6, sympy residuals "
           f"{[str(x) for x in dq if x != 0] or '0'}; f2 = {fs[2]}, f3 = {fs[3]}, f4 = {fs[4]}, f5 = {fs[5]}, f6 = {fs[6]} ({time.time() - t0:.1f} s)")


# ---------------------------------------------------------------------------------------------------------------
# the evaluations and their checks
# ---------------------------------------------------------------------------------------------------------------
def check_evals():
    mp.mp.dps = 60
    ev = {}
    main_grids = CGRIDS if CGRIDS_ON else GRIDS
    todo = list(main_grids.items()) + ([] if CGRIDS_ON else list(REPORT_GRIDS.items()))
    if GRIDSEL:
        todo = [(n, g) for n, g in todo if n in GRIDSEL]
    for name, (N, w, V, dps) in todo:
        ev[name] = run_B(N, w, V, dps, log=True)
        e = ev[name]
        print(f"  {name}: {e['grid']}  {e['secs']:.0f} s (panel setup {e['init_secs']:.0f} s)\n    Omega0 = {nd(e['Om0'], 60)}\n    Omega1 = {nd(e['Om1'], 60)}"
              f"\n    Omega2 = {nd(e['Om2'], 60)}\n    Omega3 = {nd(e['Om3'], 110)}\n    Omega4 = {nd(e['Om4'], 200)}\n    Omega5 = {nd(e['Om5'], dps - 40)}"
              f"\n    Omega6 = {nd(e['Om6'], dps - 40)}"
              "\n    slopes G2 .. G8 = " + " ".join(nd(e[f"slope{k}"], 20) for k in range(2, 9))
              + "\n    |G_k| at the root / |slope|, k = 2 .. 8 = " + " ".join(nd(e[f"resid{k}"], 3) for k in range(2, 9)), flush=True)
    if SAVE:
        save_ev(ev, SAVE)
    return ev


OM6_REF = "3.67e6"          # a reference for PLW of `wrong_Om5` (the estimate of omega5.md; only a change of 1e-3 is read)


def check_values(ev):
    mp.mp.dps = max(int(e["grid"].split("dps=")[1].split()[0]) for e in ev.values()) + 30
    cl = closed_omegas()
    grids = ["G1", "G2"]
    # PL0 .. PL5: the planted orders by the same code
    for i, cid in enumerate(("PL0", "PL1", "PL2", "PL3", "PL4", "PL5")):
        rows = {n: digits(cl[i], ev[n][f"Om{i}"]) for n in ev}
        ok = all(rows[n] >= 100 for n in grids)
        report(cid, ok, f"planted: Omega{i} against its closed form, digits: " + ", ".join(f"{n} {float(d):.1f}" for n, d in rows.items())
               + " (bar 100 on G1, G2; the others reported)")
    # PLW: a wrong input must move the orders above it (the control's own statement, evaluated when a wrong-input mutation is on)
    if MUTATE in WRONG_OM:
        i = int(MUTATE[-1]) + 1
        ref = cl[i] if i < 6 else mp.mpf(OM6_REF)
        moved = all(abs(ev[n][f"Om{i}"] - ref) > mp.mpf(1) / 1000 for n in grids)
        report("PLW", moved, f"control {MUTATE} moves Omega{i} of " + ", ".join(f"{n} by {nd(ev[n][f'Om{i}'] - ref, 5)}" for n in grids))
    # SL: slopes
    rows, ok = [], True
    for n in ev:
        ds = [digits(-4 * mp.pi, ev[n][f"slope{k}"]) for k in range(2, 9)]
        rows.append(f"{n} " + "/".join(f"{float(d):.1f}" for d in ds))
        ok = ok and (n not in grids or min(ds) >= 100)
    report("SL", ok, "slope of G_k in Omega_(k-2) against -4 pi, k = 2 .. 8, digits: " + ", ".join(rows) + " (bar 100 on G1, G2)")
    # RES: G_k at the root
    rs = {n: max(ev[n][f"resid{k}"] for k in range(2, 9)) for n in ev}
    report("RES", all(rs[n] <= mp.mpf(10) ** -100 for n in grids),
           "max over k = 2 .. 8 of |G_k at the chain's root| / |slope|: " + ", ".join(f"{n} {nd(r, 3)}" for n, r in rs.items()) + " (bar 1e-100 on G1, G2)")
    # M6: the two grids on Omega6
    a, b = grids
    d6 = digits(ev[a]["Om6"], ev[b]["Om6"])
    report("M6", d6 >= 100, f"Omega6 of {a} against {b} (two grids, two precisions): {float(d6):.1f} digits (bar 100); |diff| = "
           f"{nd(abs(ev[a]['Om6'] - ev[b]['Om6']), 3)}")
    # RNG: the range of the estimate
    report("RNG", all(mp.mpf("3.5e6") <= ev[n]["Om6"] <= mp.mpf("3.9e6") for n in grids), "Omega6 between 3.5e6 and 3.9e6: "
           + ", ".join(f"{n} {nd(ev[n]['Om6'], 12)}" for n in grids))
    others = [n for n in ev if n not in grids]
    print("  reported: digits between evaluations, Omega5 / Omega6: " + ", ".join(
        f"{x}-{y} " + "/".join(f"{float(digits(ev[x][f'Om{i}'], ev[y][f'Om{i}'])):.1f}" for i in (5, 6))
        for x, y in [("G1", "G2")] + [(o, "G1") for o in others]))
    return d6


# ---------------------------------------------------------------------------------------------------------------
# PSLQ of f7
# ---------------------------------------------------------------------------------------------------------------
PSLQ_DIGITS = 390             # cap on the working digits; the digits the two grids share, less 32, if fewer
PSLQ_SLACK = 12               # tolerance 10^-(d - slack)
VERIFY_EXTRA = 20             # a relation must hold at this many more digits than were used to find it


def free_constants():
    """the ln 2-free monomials of pi^2, zeta(3), zeta(5), zeta(7) of weight <= 7 (pi^2: 2, zeta(3): 3, zeta(5): 5, zeta(7): 7, products add): (name, weight, value);
    the 11 constants span the non-alternating Euler sums of weight <= 7 (Zagier's dimensions 1, 0, 1, 1, 1, 2, 2, 3 for weights 0 .. 7)"""
    pi, z3, z5, z7 = mp.pi, mp.zeta(3), mp.zeta(5), mp.zeta(7)
    p2 = pi**2
    return [("pi^4 zeta(3)", 7, p2**2 * z3), ("pi^2 zeta(5)", 7, p2 * z5), ("zeta(7)", 7, z7), ("pi^6", 6, p2**3), ("zeta(3)^2", 6, z3**2),
            ("pi^2 zeta(3)", 5, p2 * z3), ("zeta(5)", 5, z5), ("pi^4", 4, p2**2), ("zeta(3)", 3, z3), ("pi^2", 2, p2), ("1", 0, mp.mpf(1))]


# basis name -> (lowest weight of the ln 2-free monomials, highest power of ln 2, maxcoeff H, planted seed); ln 2^j times the free monomials of weight <= 7 - j
PSLQ_BASES = {"F1": (7, 0, 10**16, 282001), "F2": (6, 0, 10**16, 282002), "F3": (5, 0, 10**16, 282003), "F4": (4, 0, 10**16, 282004),
              "F5": (3, 0, 10**16, 282005), "F6": (0, 0, 10**16, 282006),
              "X1": (0, 1, 10**12, 282007), "X2": (0, 2, 10**11, 282008), "X3": (0, 3, 10**10, 282009), "XF": (0, 7, 10**8, 282010)}


def constants(name):
    lo, jmax, _, _ = PSLQ_BASES[name]
    q = mp.ln2
    out = []
    for j in range(jmax + 1):
        for nm, wt, v in free_constants():
            if wt >= lo and wt <= 7 - j:
                pw = "ln 2" if j == 1 else f"ln^{j} 2"
                out.append((nm if j == 0 else f"{pw} {nm}" if nm != "1" else pw, q**j * v))
    return out


def planted(n, seed):
    """the planted relation: denominator and n integer coefficients, drawn once from a fixed seed (declared in the raw README)"""
    r = random.Random(seed)
    den = r.choice([3, 5, 7, 9, 11, 13, 15, 21, 45])
    return den, [r.choice([-1, 1]) * r.randint(100, 9999) for _ in range(n)]


def pslq_run(x, tau, H):
    t0 = time.time()
    r = mp.pslq(x, tol=tau, maxcoeff=H, maxsteps=10**6)
    return r, time.time() - t0


def residual(r, x):
    return abs(sum(c * v for c, v in zip(r, x)))


def f7_target(ev, g):
    """f7 of a grid: its Omega6 and the exact Omega0 .. Omega5; `ln2_leak` adds 1000 ln 2 zeta(3)^2 (a ln 2 term in F at eps^7)"""
    f = f7_of(ev[g]["Om6"])
    if MUTATE == "ln2_leak":
        f = f + 1000 * mp.ln2 * mp.zeta(3) ** 2
    return f


def check_pslq(ev):
    grids = ["G1", "G2"]
    mp.mp.dps = max(int(e["grid"].split("dps=")[1].split()[0]) for e in ev.values()) + 30
    agree = digits(ev[grids[0]]["Om6"], ev[grids[1]]["Om6"])
    cap = PSLQ_CAP or PSLQ_DIGITS
    d = int(max(30, min(cap, math.floor(agree) - (PSLQ_SLACK + VERIFY_EXTRA))))      # PSLQ needs at least 53 bits
    pert = mp.sqrt(2) / 1000        # irrational, so that the shifted constant is outside every basis's rational span (1 is in some bases)
    print(f"  PSLQ at {d} digits (the two grids share {float(agree):.1f}), tolerance 1e-{d - PSLQ_SLACK}, verification at {d + VERIFY_EXTRA} digits", flush=True)
    saved = {}
    for name, (_, _, H, seed) in PSLQ_BASES.items():
        if BASES_SEL and name not in BASES_SEL:
            continue
        with mp.workdps(d):
            tau = mp.mpf(10) ** (-(d - PSLQ_SLACK))
            names, cs = zip(*constants(name))
            cs_true = list(cs)
            cs = list(cs)
            if MUTATE == "pslq_wrong_basis":                  # PSLQ is given a basis that is not the one the target was built from
                cs[0] = cs[0] + pert
            n = len(cs)
            if CGRIDS_ON or PSLQ_CAP:                  # a control run: d is small, so H is lowered until the budget is below d - 25 (at least 10^4, the planted coefficients)
                H = 10 ** max(4, min(int(math.log10(H)), int((d - 25) // (n + 1)) - 1))
            budget = (n + 1) * math.log10(2 * H + 1) + 1
            f7 = +f7_target(ev, grids[1])
            # planted truth: den * x = sum c_i b_i, through the same call
            den, c = planted(n, seed)
            xp = sum(ci * b for ci, b in zip(c, cs_true)) / den
            if MUTATE == "pslq_wrong_target":                 # the planted target is off by 1e-20 relative
                xp = xp * (1 + mp.mpf(10) ** -20)
            rp, tp = pslq_run([xp] + cs, tau, H)
            want = [den] + [-ci for ci in c]
            okp = rp is not None and (list(rp) == want or list(rp) == [-v for v in want])
            # perturbed basis on the planted target: must find nothing
            csp = list(cs)
            csp[1 if n > 1 else 0] = csp[1 if n > 1 else 0] + pert
            rpp, _ = pslq_run([xp] + csp, tau, H)
            report("P" + name, okp and rpp is None,
                   f"{name} = {{f7, {', '.join(names)}}}, n = {n + 1}, H = {H}, d = {d}: planted relation (den {den}, seed {seed}) "
                   f"recovered: {okp} ({tp:.0f} s); perturbed basis on the planted target finds {'nothing' if rpp is None else list(rpp)}; "
                   f"budget (n + 1) log10(2H + 1) + 1 = {budget:.1f}")
            # the real target
            if not REAL:
                print(f"    REAL {name}: skipped (--no-real)", flush=True)
                continue
            rr, tr = pslq_run([f7] + cs, tau, H)
            found, why = False, "not found"
            if rr is not None:
                res = residual(rr, [f7] + cs)
                found = res <= tau and budget < d and rr[0] != 0      # a relation without f7 is among the constants alone
                why = f"relation {list(rr)}, residual {nd(res, 3)} (bar {nd(tau, 2)}), budget with H = {H}: {budget:.1f} < {d}: {budget < d}; f7 coefficient {rr[0]}"
                if found:
                    j = next(i for i in range(1, len(rr)) if rr[i] != 0)
                    csr = list(cs)
                    csr[j - 1] = csr[j - 1] + pert
                    rrp, _ = pslq_run([f7] + csr, tau, H)
                    found = found and rrp is None
                    why += f"; perturbed basis (constant {names[j - 1]} + sqrt(2)/1000) finds {'nothing' if rrp is None else list(rrp)}"
                    # verification at VERIFY_EXTRA more digits, on both grids
                    with mp.workdps(d + VERIFY_EXTRA):
                        csv = [c_ for _, c_ in constants(name)]
                        tv = mp.mpf(10) ** -(d + VERIFY_EXTRA - PSLQ_SLACK)
                        vres = {g: residual(rr, [f7_target(ev, g)] + csv) for g in grids}
                        vok = all(r_ <= tv for r_ in vres.values())
                        found = found and vok
                        why += (f"; verified at {d + VERIFY_EXTRA} digits (bar {nd(tv, 2)}): " + ", ".join(f"{g} {nd(r_, 3)}" for g, r_ in vres.items())
                                + f" -> {vok}")
                    if found:
                        saved[name] = {"names": list(names), "rel": [int(v) for v in rr], "d": d, "H": H}
            print(f"    REAL {name}: " + (("FOUND (counts) " if found else "found but does not count ") + why if rr is not None
                                           else f"not found at H = {H} ({tr:.0f} s)"), flush=True)
    if SAVE_REL:
        with open(SAVE_REL, "w") as f:
            json.dump(saved, f, indent=1)


# ---------------------------------------------------------------------------------------------------------------
# the closed form of Omega6 from a relation of f7, and the structure verdict
# ---------------------------------------------------------------------------------------------------------------
def check_closed6(ev):
    """ST7, ID6, ZL: from the saved relations (--rel=F.json,X.json): the ln 2-free relation of F6 gives f7 and Omega6; it holds on both grids at the verification digits;
    every relation found in a basis with ln 2 has zero coefficients on all ln 2 monomials"""
    grids = ["G1", "G2"]
    saved = {}
    for path in REL.split(","):
        saved.update(json.load(open(path)))
    mp.mp.dps = 400
    have = "F6" in saved
    # the closed form of Omega6, symbolically
    R = sp.Rational
    q, z3, z5, z7, p2 = sp.symbols("q z3 z5 z7 p2")
    sym_of = {"pi^4 zeta(3)": p2**2 * z3, "pi^2 zeta(5)": p2 * z5, "zeta(7)": z7, "pi^6": p2**3, "zeta(3)^2": z3**2, "pi^2 zeta(3)": p2 * z3, "zeta(5)": z5,
              "pi^4": p2**2, "zeta(3)": z3, "pi^2": p2, "1": sp.Integer(1)}
    if have:
        rel = list(saved["F6"]["rel"])
        names = saved["F6"]["names"]
        if MUTATE == "f7_rel":
            rel[4] += 1                                       # a coefficient of the relation off by one
        f7s = -sum(R(rel[1 + i], rel[0]) * sym_of[nm] for i, nm in enumerate(names))
        Om = [-4 * q, 2 * p2 / 3 + 8 * q**2,
              80 * z3 - R(8, 3) * p2 * q - R(32, 3) * q**3,
              22 * p2**2 / 45 + R(16, 3) * p2 * q**2 + R(32, 3) * q**4 - 320 * z3 * q - 384 * z3,
              -10432 * z5 - 160 * p2 * z3 + 640 * q**2 * z3 - R(64, 9) * p2 * q**3 - R(88, 45) * p2**2 * q - R(128, 15) * q**5 - R(2432, 45) * p2**2 + 1536 * q * z3 + 3072 * z3,
              R(256, 45) * q**6 + R(64, 9) * p2 * q**4 - R(2560, 3) * z3 * q**3 + R(176, 45) * p2**2 * q**2 - 3072 * z3 * q**2 + 640 * p2 * z3 * q + R(9728, 45) * p2**2 * q
              - 12288 * z3 * q + 41728 * z5 * q - R(38812, 315) * p2**3 + 3200 * z3**2 - R(256, 3) * p2 * z3 + 74240 * z5 + R(39424, 45) * p2**2 - 36864 * z3]
        a = [sp.Integer(1)] + Om
        om6 = sp.expand(f7s - sum(a[j] * (4 * q) ** (7 - j) / sp.factorial(7 - j) for j in range(7)))
        val = {q: mp.ln2, z3: mp.zeta(3), z5: mp.zeta(5), z7: mp.zeta(7), p2: mp.pi**2}
        f = sp.lambdify((q, z3, z5, z7, p2), om6, "mpmath")
        om6n = f(mp.ln2, mp.zeta(3), mp.zeta(5), mp.zeta(7), mp.pi**2)
        ds = {g: digits(om6n, ev[g]["Om6"]) for g in grids}
        report("ID6", all(d_ >= 380 for d_ in ds.values()),
               f"Omega6 from the ln 2-free relation of F6 (f7 = -sum c_i b_i / c_0, Omega6 = f7 - the exact ln 2 terms) = {nd(om6n, 60)} against the grids, digits: "
               + ", ".join(f"{g} {float(d_):.1f}" for g, d_ in ds.items()) + " (bar 380; the working precision here is 400 digits)")
        print("  Omega6 = " + str(om6))
        wt = {q: 1, p2: 2, z3: 3, z5: 5, z7: 7}
        pol = sp.Poly(om6, q, p2, z3, z5, z7)
        byw = {}
        for mon, co in pol.terms():
            w_ = mon[0] + 2 * mon[1] + 3 * mon[2] + 5 * mon[3] + 7 * mon[4]
            byw.setdefault(w_, []).append((mon, co))
        for w_ in sorted(byw, reverse=True):
            part = sum(mp.mpf(co.p) / co.q * mp.ln2 ** m[0] * mp.pi ** (2 * m[1]) * mp.zeta(3) ** m[2] * mp.zeta(5) ** m[3] * mp.zeta(7) ** m[4] for m, co in byw[w_])
            print(f"    weight {w_}: {len(byw[w_])} terms, part {nd(part, 20)}")
    else:
        report("ID6", False, "no relation of F6 was found (the file has " + (", ".join(saved) or "no relation") + "): f7 is not in the span of the ln 2-free monomials of weight <= 7 at these digits")
    # ZL: every relation found with ln 2 has zero coefficients on the ln 2 monomials, and is the F6 relation
    bad = []
    for nm, r in saved.items():
        if nm.startswith("X"):
            nz = [r["names"][i] for i in range(len(r["names"])) if "ln" in r["names"][i] and r["rel"][1 + i] != 0]
            if nz:
                bad.append(f"{nm}: nonzero on {nz}")
    if MUTATE == "ln2_zl" and not bad:
        bad.append("ln2_zl: a nonzero coefficient on a ln 2 monomial put in by the control")
    ok = have and not bad
    report("ZL", ok, "F6 found: " + str(have) + "; relations found in the bases with ln 2 (" + (", ".join(n_ for n_ in saved if n_.startswith("X")) or "none") + "): "
           + ("all with zero coefficients on every ln 2 monomial" if not bad else "; ".join(bad)))


# ---------------------------------------------------------------------------------------------------------------
# the remainder against #195's W^2
# ---------------------------------------------------------------------------------------------------------------
def read_w2(path):
    """{D: W2 string} of the 64-value list of w2_output.txt (qd values, 60 digits)"""
    lst = {}
    for line in open(path):
        m = re.match(r"\s*D=(\d+)\s+(\S+)\s+digits qd/od", line)
        if m:
            lst[int(m.group(1))] = m.group(2)
    return lst


def read_w2raw(path):
    """{D: W2 string} of the RESULT lines of I195-dh610/run_02 (qd, 60 digits) at D = 200, 500, 1000, 2000"""
    out = {}
    for line in open(path):
        m = re.search(r'"D":(\d+),"W2":"([^"]+)"', line)
        if m:
            out[int(m.group(1))] = m.group(2)
    return out


def refit_windows(lst, om, n, Ks=(14, 18, 22, 26, 30, 36, 44)):
    """Omega_n from the 64 values of W^2 with Omega0 .. Omega_{n-1} fixed (om): y = (W^2 - [s + Omega0 + Omega1/s + ... + Omega_{n-1}/s^(n-1)]) s^n =
    Omega_n + Omega_{n+1}/s + ..., polynomial extrapolation to 1/s = 0 through the K points of largest s; returns {K: estimate}"""
    pts = []
    for D in sorted(lst):
        s = mp.mpf(2 * (D - 3))
        pts.append((1 / s, (mp.mpf(lst[D]) - (s + sum(om[i] / s**i for i in range(n)))) * s**n))
    out = {}
    for K in Ks:
        w = pts[len(pts) - K:]
        est = mp.mpf(0)
        for i, (ti, yi) in enumerate(w):
            Li = mp.mpf(1)
            for j, (tj, _) in enumerate(w):
                if j != i:
                    Li *= (0 - tj) / (ti - tj)
            est += Li * yi
        out[K] = est
    return out


X4_BAR = mp.mpf("1e-18")      # set knowing a development extrapolation (1.1e-23 .. 4.4e-21 on a 148-digit Omega6, dev_07); see the raw README
R7_BAR = mp.mpf("5e-5")       # set knowing dev_07: relative difference 1.47e-6 (an estimate)
R8_BAR = mp.mpf(10**5)        # set knowing dev_07: |x8(1e4) - x8(5103)| = 52722


def check_remainder(ev):
    mp.mp.dps = 80
    grids = ["G1", "G2"]
    om = [ev[grids[1]][f"Om{i}"] for i in range(1, 7)]
    if MUTATE == "wrong_Om6_series":
        om[5] = om[5] + 1000
    lst = read_w2(W2FILE)
    if W2RAW:
        lst.update(read_w2raw(W2RAW))      # D = 1000, 2000 are not in the 64-value list

    def x(D, nterms):
        s = mp.mpf(2 * (D - 3))
        rem = mp.mpf(lst[D]) - (s - 4 * mp.ln2 + sum(om[i] / s ** (i + 1) for i in range(nterms)))
        return rem, rem * s ** (nterms + 1)
    print("  W^2 (qd, w2_output.txt) against the series s - 4 ln 2 + Omega1/s + ... (Omega1 .. Omega6 of G2), s = 2(D - 3):")
    print("  D       s       rem7 s^6 (tends to Omega6)    (rem7 s^6 - Omega6) s (tends to Omega7)   rem8 s^7 (tends to Omega7)")
    xs7, xs8 = {}, {}
    for D in (1000, 2000, 5103, 10000):
        s = mp.mpf(2 * (D - 3))
        _, xs7[D] = x(D, 5)
        _, xs8[D] = x(D, 6)
        print(f"  {D:<6d} {2 * (D - 3):<7d} {nd(xs7[D], 14):<28s} {nd((xs7[D] - om[5]) * s, 12):<40s} {nd(xs8[D], 12)}")
    s1, s2 = mp.mpf(2 * (5103 - 3)), mp.mpf(2 * (10000 - 3))
    est7 = (s2 * xs7[10000] - s1 * xs7[5103]) / (s2 - s1)
    report("R7", abs(est7 - om[5]) <= R7_BAR * abs(om[5]), f"two-point extrapolation (s2 x2 - s1 x1)/(s2 - s1) of rem7 s^6 from D = 5103, 10^4: {nd(est7, 8)} against the computed "
           f"Omega6 {nd(om[5], 8)}, relative difference {nd((est7 - om[5]) / om[5], 3)} (bar {R7_BAR}: an estimate, not a bound)")
    est8 = (s2 * xs8[10000] - s1 * xs8[5103]) / (s2 - s1)
    print(f"  reported: two-point extrapolation of rem8 s^7 from D = 5103, 10^4: Omega7 ~ {nd(est8, 8)} (an estimate, not a bound)")
    report("R8", abs(xs8[10000] - xs8[5103]) <= R8_BAR, f"rem8 s^7 at D = 5103 and 10^4 differs by {nd(abs(xs8[10000] - xs8[5103]), 5)} (bar {R8_BAR}: a wrong Omega6 adds (error) s, about 2e4 per unit)")
    # X4: Omega6 from the crate's 64 values with Omega0 .. Omega5 fixed (the closed forms)
    mp.mp.dps = 120
    est = refit_windows(read_w2(W2FILE), closed_omegas(), 6)
    rows = {n: {K: abs(e - ev[n]["Om6"]) for K, e in est.items()} for n in grids}
    mx = max(max(r.values()) for r in rows.values())
    report("X4", mx <= X4_BAR, "Omega6 against the polynomial extrapolation in 1/s of #195's 64 W^2 values (Omega0 .. Omega5 fixed at the closed forms), "
           "windows of the K largest s, |difference| " + "; ".join(f"{n}: " + ", ".join(f"K = {K} {nd(d, 2)}" for K, d in r.items()) for n, r in rows.items())
           + f" (bar {X4_BAR} on every window; the spread over windows is an estimator range, {nd(mx, 2)}; not a bound)")
    mp.mp.dps = 80
    print("  reported, not gated: rem7 s^6 at the D of the list")
    print("  D       rem7 s^6")
    for D in sorted(lst):
        print(f"  {D:<7d} {nd(x(D, 5)[1], 14)}")


def main():
    parse(sys.argv[1:])
    mp.mp.dps = 60
    print(f"order8.py  mutate={MUTATE}  parts={PARTS}  mpmath {mp.__version__}  sympy {sp.__version__}", flush=True)
    t0 = time.time()
    if "grid" in PARTS:
        check_grid()
    if "cache" in PARTS:
        check_cache()
    if "sym" in PARTS:
        check_sym()
    if "struct" in PARTS:
        check_struct()
    ev = None
    if "evals" in PARTS:
        ev = check_evals()
        if "values" in PARTS or not GRIDSEL:
            check_values(ev)
    elif LOAD and any(p in PARTS for p in ("values", "pslq", "remainder", "closed6")):
        mp.mp.dps = 600
        ev = load_ev(LOAD)
        if "values" in PARTS:
            check_values(ev)
    if "pslq" in PARTS:
        check_pslq(ev)
    if "closed6" in PARTS:
        check_closed6(ev)
    if "remainder" in PARTS:
        check_remainder(ev)
    nfail = sum(1 for _, ok in RESULTS if not ok)
    print(f"\n{len(RESULTS) - nfail} passed, {nfail} failed of {len(RESULTS)} checks; {time.time() - t0:.0f} s")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()

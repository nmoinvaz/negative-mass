"""Order 7 of the matched expansion (matched.py, order3.py ... order6.py): Omega5 from G7 = 0 by method B, and the checks of issue #263.

W^2 a^2 = s + Omega, s = 2(D - 3), eps = 1/s, Omega = Omega0 + Omega1 eps + Omega2 eps^2 + ..., R = r^s, l = ln R = -2 ln cos(theta).
The order equations (E0 derives them from the master equation), the same for every k:
  dY_k/dl = F_{k-1},  F_{k-1} = sum_{j<k} L_{k-1-j} l^j / j!,   (R - 1) L_k' - L_k/2 = h_k,
  h_k = (1 - 1/R) L_{k-1} + sum_{j<k} Y_j c_{k-j} + Y_k/2,   L_k = -2 Omega_{k-1} + u C_k,   G_k = C_k(inf) = 0,
with u = sqrt(R - 1), C_k = int_1^R (h_k - Omega_{k-1}) (R' - 1)^{-3/2} dR'. G_k is affine in Omega_{k-2} with slope -4 pi:
G2 = 0 gives Omega0, ..., G6 = 0 Omega4, **G7 = 0 gives Omega5**. One path computes them all (`omega_chain`).

Method B (order4.py's PanelB, shared): v = sqrt(l), panels of width w on [0, V], N Chebyshev nodes per panel, cumulative integrals by the
spectral integration matrix; Y_k from the left, the tail C_right from the right, L_k = -2 h_k(1) - u C_right, G_k = C_right(0). It uses no
closed form of Y_k or L_k and no input value of Omega_0 .. Omega_4. Levels below k - 1 do not depend on Omega_{k-2}, so the chain keeps them
(`cached`; the same operations on the same numbers as the uncached chain, which CACHE compares).

    python3 order7.py grid --grid=N,w,V,dps [--kmax=7]      (a development run of one grid; no gate)
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

MUTATIONS = ("wrong_h", "wrong_c6", "wrong_c7", "wrong_Om0", "wrong_Om1", "wrong_Om2", "wrong_Om3", "wrong_Om4", "coarse_B", "short_tail",
             "stale_cache", "pslq_wrong_basis", "pslq_wrong_target", "wrong_Om5_series", "wrong_closed4", "e1_sign", "e2_weight", "e3_factor",
             "e4_g", "e4_noC", "e2w_g", "e2w_P", "e5_noD3", "e5_const", "rd_kit", "rd_Dstep", "rd_J5", "rd_master", "rd_B")
MUTATE = None
W2FILE = "w2_output.txt"
SAVE = None
LOAD = None
GRIDARG = None
GRIDSEL = None
CGRIDS_ON = False
W2RAW = None
KMAX = 7
PARTS = []
REAL = True
BASES_SEL = None
PSLQ_CAP = None


def parse(argv):
    global MUTATE, W2FILE, W2RAW, SAVE, LOAD, GRIDARG, GRIDSEL, CGRIDS_ON, KMAX, PARTS, REAL, BASES_SEL, PSLQ_CAP
    parts = []
    for a in argv:
        if a.startswith("--mutate="):
            MUTATE = a.split("=", 1)[1]
            assert MUTATE in MUTATIONS, f"unknown mutation {MUTATE}; known: {MUTATIONS}"
        elif a.startswith("--w2="):
            W2FILE = a.split("=", 1)[1]
        elif a.startswith("--w2raw="):
            W2RAW = a.split("=", 1)[1]
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
    PARTS = parts or ["sym", "evals", "values", "pslq", "remainder"]


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
    raise ValueError(m)


def hk(m, Om, l, R, Y, L):
    """h_m"""
    if m == 0:
        return Y[0] / 2
    return (1 - 1 / R) * L[m - 1] + sum(Y[j] * cc(m - j, Om, l, R) for j in range(m)) + Y[m] / 2


def Fsum(k, L, l):
    return sum(L[k - 1 - j] * l**j / math.factorial(j) for j in range(k))


# ---------------------------------------------------------------------------------------------------------------
# E0: the epsilon-expansion of the master equation against the table, k = 0 .. 7
# ---------------------------------------------------------------------------------------------------------------
def check_sym():
    t0 = time.time()
    R, eps, l = sp.symbols("R epsilon l", positive=True)
    Om = sp.symbols("Omega0:7")
    Ys, Ls = sp.symbols("Y0:8"), sp.symbols("L0:8")
    kk = 7
    r = sum((eps * l) ** j / sp.factorial(j) for j in range(kk + 1))          # exact through eps^7, which is all that is read
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
           f"h_k (k = 0 .. 7) and dY_k/dl (k = 1 .. 7) from the epsilon-expansion of the master equation against the table: "
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
    """Omega_0 .. Omega_4 of #26, #208 and #220 (derived there), the planted truth"""
    pi, ln2, z3, z5 = mp.pi, mp.ln2, mp.zeta(3), mp.zeta(5)
    return [-4 * ln2, 2 * pi**2 / 3 + 8 * ln2**2,
            80 * z3 - mp.mpf(8) / 3 * pi**2 * ln2 - mp.mpf(32) / 3 * ln2**3,
            22 * pi**4 / 45 + mp.mpf(16) / 3 * pi**2 * ln2**2 + mp.mpf(32) / 3 * ln2**4 - 320 * z3 * ln2 - 384 * z3,
            -10432 * z5 - 160 * pi**2 * z3 + 640 * ln2**2 * z3 - mp.mpf(64) / 9 * pi**2 * ln2**3 - mp.mpf(88) / 45 * pi**4 * ln2
            - mp.mpf(128) / 15 * ln2**5 - mp.mpf(2432) / 45 * pi**4 + 1536 * ln2 * z3 + 3072 * z3
            + (mp.mpf(1) / 1000 if MUTATE == "wrong_closed4" else 0)]


WRONG_OM = ("wrong_Om0", "wrong_Om1", "wrong_Om2", "wrong_Om3", "wrong_Om4")


def omega_chain_ref(B, kmax, log=None):
    """Omega_0 .. Omega_{kmax-2} from G_2 .. G_kmax = 0, in turn, each G_k from the whole hierarchy (the uncached path); returns (Om, slopes, secs)"""
    Om = [mp.mpf(0)] * kmax
    slopes, secs = {}, {}
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
        secs[k] = time.time() - t0
        if log:
            print(f"    k = {k}: Omega_{k - 2} = {nd(val, 30)}  slope {nd(slope, 20)}  ({secs[k]:.0f} s)", flush=True)
    return Om, slopes, secs, {}


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
            return Gt[k], Lt[k - 1]
        g0, l0 = trial(mp.mpf(0))
        g1, _ = trial(mp.mpf(1))
        val, slope = -g0 / (g1 - g0), g1 - g0
        Om[k - 2] = val
        slopes[k] = slope
        if MUTATE == "stale_cache" and k >= 3:
            L.append(l0)                                          # level k - 1 from the first trial, not recomputed at the root
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


def run_B(N, w, V, dps, kmax=7, log=False, cached=True):
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
GRIDS = {"G1": (285, (1, 4), 48, 445), "G2": (275, (3, 10), 47, 410)}             # w as a fraction (num, den), evaluated at the working precision
REPORT_GRIDS = {"B1": (60, (2, 5), 24, 100), "B2": (100, (1, 4), 30, 170), "B3": (170, (1, 4), 38, 285)}
CGRIDS = {"G1": (90, (1, 4), 29, 160), "G2": (85, (3, 10), 28, 150)}


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
    for i in range(min(5, KMAX - 1)):
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
    a = run_B(N, w, V, dps, kmax=7)
    b = run_B(N, w, V, dps, kmax=7, cached=False)
    mp.mp.dps = dps
    worst = max(abs(a[f"Om{i}"] - b[f"Om{i}"]) / abs(b[f"Om{i}"]) for i in range(6))
    dig = [float(digits(b[f"Om{i}"], a[f"Om{i}"])) for i in range(6)]
    report("CACHE", worst <= mp.mpf(10) ** -(dps - 5),
           f"cached chain against the uncached chain on (N {N}, w 1/2, V {V}, {dps} digits): Omega0 .. Omega5 equal to "
           + ", ".join(f"{d:.1f}" for d in dig) + f" digits, worst relative difference {nd(worst, 3)} (bar 1e-{dps - 5}); times {a['secs']:.0f} s and {b['secs']:.0f} s")


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
        ev[name] = run_B(N, w, V, dps)
        e = ev[name]
        print(f"  {name}: {e['grid']}  {e['secs']:.0f} s (panel setup {e['init_secs']:.0f} s)\n    Omega0 = {nd(e['Om0'], 60)}\n    Omega1 = {nd(e['Om1'], 60)}"
              f"\n    Omega2 = {nd(e['Om2'], 60)}\n    Omega3 = {nd(e['Om3'], 110)}\n    Omega4 = {nd(e['Om4'], 200)}\n    Omega5 = {nd(e['Om5'], dps - 40)}"
              "\n    slopes G2 .. G7 = " + " ".join(nd(e[f"slope{k}"], 20) for k in range(2, 8))
              + "\n    |G_k| at the root / |slope|, k = 2 .. 7 = " + " ".join(nd(e[f"resid{k}"], 3) for k in range(2, 8)), flush=True)
    if SAVE:
        save_ev(ev, SAVE)
    return ev


OM5_REF = "40860.4553545256873125440648260337"          # a reference for PLW of `wrong_Om4` (a development value of Omega5, 31 digits; only a change of 1e-3 is read)


def check_values(ev):
    mp.mp.dps = max(int(e["grid"].split("dps=")[1].split()[0]) for e in ev.values()) + 30
    cl = closed_omegas()
    grids = ["G1", "G2"]
    # PL0 .. PL4: the planted orders by the same code
    for i, cid in enumerate(("PL0", "PL1", "PL2", "PL3", "PL4")):
        rows = {n: digits(cl[i], ev[n][f"Om{i}"]) for n in ev}
        ok = all(rows[n] >= 100 for n in grids)
        report(cid, ok, f"planted: Omega{i} against its closed form, digits: " + ", ".join(f"{n} {float(d):.1f}" for n, d in rows.items())
               + " (bar 100 on G1, G2; the others reported)")
    # PLW: a wrong input must move the orders above it (the control's own statement, evaluated when a wrong-input mutation is on)
    if MUTATE in WRONG_OM:
        i = int(MUTATE[-1]) + 1
        ref = cl[i] if i < 5 else mp.mpf(OM5_REF)
        moved = all(abs(ev[n][f"Om{i}"] - ref) > mp.mpf(1) / 1000 for n in grids)
        report("PLW", moved, f"control {MUTATE} moves Omega{i} of " + ", ".join(f"{n} by {nd(ev[n][f'Om{i}'] - ref, 5)}" for n in grids))
    # SL: slopes
    rows, ok = [], True
    for n in ev:
        ds = [digits(-4 * mp.pi, ev[n][f"slope{k}"]) for k in range(2, 8)]
        rows.append(f"{n} " + "/".join(f"{float(d):.1f}" for d in ds))
        ok = ok and (n not in grids or min(ds) >= 100)
    report("SL", ok, "slope of G_k in Omega_(k-2) against -4 pi, k = 2 .. 7, digits: " + ", ".join(rows) + " (bar 100 on G1, G2)")
    # RES: G_k at the root
    rs = {n: max(ev[n][f"resid{k}"] for k in range(2, 8)) for n in ev}
    report("RES", all(rs[n] <= mp.mpf(10) ** -100 for n in grids),
           "max over k = 2 .. 7 of |G_k at the chain's root| / |slope|: " + ", ".join(f"{n} {nd(r, 3)}" for n, r in rs.items()) + " (bar 1e-100 on G1, G2)")
    # M5: the two grids on Omega5
    a, b = grids
    d5 = digits(ev[a]["Om5"], ev[b]["Om5"])
    report("M5", d5 >= 100, f"Omega5 of {a} against {b} (two grids, two precisions): {float(d5):.1f} digits (bar 100); |diff| = "
           f"{nd(abs(ev[a]['Om5'] - ev[b]['Om5']), 3)}")
    # RNG: the range of the estimate
    report("RNG", all(mp.mpf("4.0e4") <= ev[n]["Om5"] <= mp.mpf("4.2e4") for n in grids), "Omega5 between 4.0e4 and 4.2e4: "
           + ", ".join(f"{n} {nd(ev[n]['Om5'], 12)}" for n in grids))
    others = [n for n in ev if n not in grids]
    print("  reported: digits between evaluations, Omega4 / Omega5: " + ", ".join(
        f"{x}-{y} " + "/".join(f"{float(digits(ev[x][f'Om{i}'], ev[y][f'Om{i}'])):.1f}" for i in (4, 5))
        for x, y in [("G1", "G2")] + [(o, "G1") for o in others]))
    return d5


# ---------------------------------------------------------------------------------------------------------------
# PSLQ
# ---------------------------------------------------------------------------------------------------------------
PSLQ_DIGITS = 330             # cap on the working digits; the digits the two grids share, less 32, if fewer
PSLQ_SLACK = 12               # tolerance 10^-(d - slack)
VERIFY_EXTRA = 20             # a relation must hold at this many more digits than were used to find it


def weighted_constants():
    """{weight: [(name, value)]}: the monomials ln2^a (pi^2)^b zeta(3)^c zeta(5)^e zeta(7)^f of weight a + 2b + 3c + 5e + 7f, and the
    polylogarithm products Li4(1/2) x {1, ln 2, ln^2 2, pi^2}, Li5(1/2) x {1, ln 2}, Li6(1/2)"""
    pi, q, z3, z5, z7 = mp.pi, mp.ln2, mp.zeta(3), mp.zeta(5), mp.zeta(7)
    p2 = pi**2
    li4, li5, li6 = mp.polylog(4, mp.mpf(1) / 2), mp.polylog(5, mp.mpf(1) / 2), mp.polylog(6, mp.mpf(1) / 2)
    w = {0: [("1", mp.mpf(1))], 1: [("ln 2", q)], 2: [("ln^2 2", q**2), ("pi^2", p2)],
         3: [("ln^3 2", q**3), ("pi^2 ln 2", p2 * q), ("zeta(3)", z3)],
         4: [("ln^4 2", q**4), ("pi^2 ln^2 2", p2 * q**2), ("pi^4", p2**2), ("zeta(3) ln 2", z3 * q)],
         5: [("ln^5 2", q**5), ("pi^2 ln^3 2", p2 * q**3), ("pi^4 ln 2", p2**2 * q), ("zeta(3) ln^2 2", z3 * q**2), ("pi^2 zeta(3)", p2 * z3),
             ("zeta(5)", z5)],
         6: [("ln^6 2", q**6), ("pi^2 ln^4 2", p2 * q**4), ("pi^4 ln^2 2", p2**2 * q**2), ("pi^6", p2**3), ("zeta(3) ln^3 2", z3 * q**3),
             ("pi^2 zeta(3) ln 2", p2 * z3 * q), ("zeta(3)^2", z3**2), ("zeta(5) ln 2", z5 * q)],
         7: [("ln^7 2", q**7), ("pi^2 ln^5 2", p2 * q**5), ("pi^4 ln^3 2", p2**2 * q**3), ("pi^6 ln 2", p2**3 * q), ("zeta(3) ln^4 2", z3 * q**4),
             ("pi^2 zeta(3) ln^2 2", p2 * z3 * q**2), ("pi^4 zeta(3)", p2**2 * z3), ("zeta(3)^2 ln 2", z3**2 * q), ("zeta(5) ln^2 2", z5 * q**2),
             ("pi^2 zeta(5)", p2 * z5), ("zeta(7)", z7)]}
    li = {4: [("Li4(1/2)", li4)], 5: [("Li4(1/2) ln 2", li4 * q), ("Li5(1/2)", li5)], 6: [("Li4(1/2) ln^2 2", li4 * q**2), ("Li4(1/2) pi^2", li4 * p2),
                                                                                           ("Li5(1/2) ln 2", li5 * q), ("Li6(1/2)", li6)]}
    return w, li


# basis name -> (weights of the monomials, weights of the polylogarithm products, maxcoeff H, planted seed); the bases of the raw README
PSLQ_BASES = {"Q1": ((6,), (), 10**12, 263001), "Q2": ((6,), (6,), 10**12, 263002), "Q3": ((6, 5), (), 10**11, 263003),
              "Q4": ((6, 5, 4), (), 10**10, 263004), "Q5": ((6, 5, 4, 3), (), 10**10, 263005), "Q6": ((6, 5, 4, 3, 2, 1, 0), (), 10**9, 263006),
              "Q7": ((6, 5, 4, 3, 2, 1, 0), (4, 5, 6), 10**9, 263007), "Q8": ((6, 5, 4, 3), (4, 5, 6), 10**9, 263008),
              "Q9": ((7, 6, 5, 4, 3, 2, 1, 0), (), 10**8, 263009)}


def constants(name):
    w, li = weighted_constants()
    mons, lis, _, _ = PSLQ_BASES[name]
    return [c for k in mons for c in w[k]] + [c for k in lis for c in li[k]]


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


def check_pslq(ev):
    grids = ["G1", "G2"]
    mp.mp.dps = max(int(e["grid"].split("dps=")[1].split()[0]) for e in ev.values()) + 30
    agree = digits(ev[grids[0]]["Om5"], ev[grids[1]]["Om5"])
    cap = PSLQ_CAP or PSLQ_DIGITS
    d = int(max(30, min(cap, math.floor(agree) - (PSLQ_SLACK + VERIFY_EXTRA))))      # PSLQ needs at least 53 bits
    pert = mp.sqrt(2) / 1000        # irrational, so that the shifted constant is outside every basis's rational span (1 is in some bases)
    print(f"  PSLQ at {d} digits (the two grids share {float(agree):.1f}), tolerance 1e-{d - PSLQ_SLACK}, verification at {d + VERIFY_EXTRA} digits", flush=True)
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
                H = 10 ** max(4, int((d - 25) // (n + 1)) - 1)
            budget = (n + 1) * math.log10(2 * H + 1) + 1
            om5 = +ev[grids[1]]["Om5"]
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
                   f"{name} = {{Omega5, {', '.join(names)}}}, n = {n + 1}, H = {H}, d = {d}: planted relation (den {den}, seed {seed}) "
                   f"recovered: {okp} ({tp:.0f} s); perturbed basis on the planted target finds {'nothing' if rpp is None else list(rpp)}; "
                   f"budget (n + 1) log10(2H + 1) + 1 = {budget:.1f}")
            # the real target
            if not REAL:
                print(f"    REAL {name}: skipped (--no-real)", flush=True)
                continue
            rr, tr = pslq_run([om5] + cs, tau, H)
            found, why = False, "not found"
            if rr is not None:
                res = residual(rr, [om5] + cs)
                found = res <= tau and budget < d and rr[0] != 0      # a relation without Omega5 is among the constants alone
                why = f"relation {list(rr)}, residual {nd(res, 3)} (bar {nd(tau, 2)}), budget with H = {H}: {budget:.1f} < {d}: {budget < d}; Omega5 coefficient {rr[0]}"
                if found:
                    j = next(i for i in range(1, len(rr)) if rr[i] != 0)
                    csr = list(cs)
                    csr[j - 1] = csr[j - 1] + pert
                    rrp, _ = pslq_run([om5] + csr, tau, H)
                    found = found and rrp is None
                    why += f"; perturbed basis (constant {names[j - 1]} + sqrt(2)/1000) finds {'nothing' if rrp is None else list(rrp)}"
                    # verification at VERIFY_EXTRA more digits, on both grids
                    with mp.workdps(d + VERIFY_EXTRA):
                        csv = [c_ for _, c_ in constants(name)]
                        tv = mp.mpf(10) ** -(d + VERIFY_EXTRA - PSLQ_SLACK)
                        vres = {g: residual(rr, [ev[g]["Om5"]] + csv) for g in grids}
                        vok = all(r_ <= tv for r_ in vres.values())
                        found = found and vok
                        why += (f"; verified at {d + VERIFY_EXTRA} digits (bar {nd(tv, 2)}): " + ", ".join(f"{g} {nd(r_, 3)}" for g, r_ in vres.items())
                                + f" -> {vok}")
            print(f"    REAL {name}: " + (("FOUND (counts) " if found else "found but does not count ") + why if rr is not None
                                           else f"not found at H = {H} ({tr:.0f} s)"), flush=True)


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


def refit_windows(lst, om, Ks=(14, 18, 22, 26, 30, 36, 44)):
    """Omega5 from the 64 values of W^2 with Omega0 .. Omega4 fixed: y = (W^2 - [s + Omega0 + Omega1/s + ... + Omega4/s^4]) s^5 =
    Omega5 + Omega6/s + ..., polynomial extrapolation to 1/s = 0 through the K points of largest s; returns {K: estimate}"""
    pts = []
    for D in sorted(lst):
        s = mp.mpf(2 * (D - 3))
        pts.append((1 / s, (mp.mpf(lst[D]) - (s + sum(om[i] / s**i for i in range(5)))) * s**5))
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


X3_BAR = mp.mpf("1e-22")      # set knowing a development extrapolation (7.5e-28 .. 3.3e-25 on a 148-digit Omega5); see the raw README
R7_BAR = mp.mpf(2500)         # set knowing rem7 s^6 on a development Omega5: |x(1e4) - x(5103)| = 1219
R6_BAR = mp.mpf("1e-4")


def check_remainder(ev):
    mp.mp.dps = 80
    grids = ["G1", "G2"]
    om = [ev[grids[1]][f"Om{i}"] for i in range(1, 6)]
    if MUTATE == "wrong_Om5_series":
        om[4] = om[4] + 10
    lst = read_w2(W2FILE)
    if W2RAW:
        lst.update(read_w2raw(W2RAW))      # D = 1000, 2000 are not in the 64-value list

    def x(D, nterms):
        s = mp.mpf(2 * (D - 3))
        rem = mp.mpf(lst[D]) - (s - 4 * mp.ln2 + sum(om[i] / s ** (i + 1) for i in range(nterms)))
        return rem, rem * s ** (nterms + 1)
    print("  W^2 (qd, w2_output.txt) against the series s - 4 ln 2 + Omega1/s + ... (Omega1 .. Omega5 of G2), s = 2(D - 3):")
    print("  D       s       rem6 s^5 (tends to Omega5)    (rem6 s^5 - Omega5) s (tends to Omega6)   rem7 s^6 (tends to Omega6)")
    xs6, xs7 = {}, {}
    for D in (1000, 2000, 5103, 10000):
        s = mp.mpf(2 * (D - 3))
        _, xs6[D] = x(D, 4)
        _, xs7[D] = x(D, 5)
        print(f"  {D:<6d} {2 * (D - 3):<7d} {nd(xs6[D], 14):<28s} {nd((xs6[D] - om[4]) * s, 12):<40s} {nd(xs7[D], 12)}")
    s1, s2 = mp.mpf(2 * (5103 - 3)), mp.mpf(2 * (10000 - 3))
    est6 = (s2 * xs6[10000] - s1 * xs6[5103]) / (s2 - s1)
    report("R6", abs(est6 - om[4]) <= R6_BAR * abs(om[4]), f"two-point extrapolation (s2 x2 - s1 x1)/(s2 - s1) of rem6 s^5 from D = 5103, 10^4: {nd(est6, 8)} against the computed "
           f"Omega5 {nd(om[4], 8)}, relative difference {nd((est6 - om[4]) / om[4], 3)} (bar 1e-4: an estimate, not a bound)")
    falls = all(xs7[a] != 0 for a in xs7) and R7_BAR is not None and abs(xs7[10000] - xs7[5103]) <= R7_BAR
    report("R7", falls, f"rem7 s^6 at D = 5103 and 10^4 differs by {nd(abs(xs7[10000] - xs7[5103]), 5)} (bar {R7_BAR}: a wrong Omega5 adds (error) s, about 10^4 per unit)")
    # X3: Omega5 from the crate's 64 values with Omega0 .. Omega4 fixed (the closed forms)
    mp.mp.dps = 120
    est = refit_windows(read_w2(W2FILE), closed_omegas())
    rows = {n: {K: abs(e - ev[n]["Om5"]) for K, e in est.items()} for n in grids}
    mx = max(max(r.values()) for r in rows.values())
    report("X3", X3_BAR is not None and mx <= X3_BAR, "Omega5 against the polynomial extrapolation in 1/s of #195's 64 W^2 values (Omega0 .. Omega4 fixed at the closed forms), "
           "windows of the K largest s, |difference| " + "; ".join(f"{n}: " + ", ".join(f"K = {K} {nd(d, 2)}" for K, d in r.items()) for n, r in rows.items())
           + f" (bar {X3_BAR} on every window; the spread over windows is an estimator range, {nd(mx, 2)}; not a bound)")
    mp.mp.dps = 80
    print("  reported, not gated: rem6 s^5 at the D of the list")
    print("  D       rem6 s^5")
    for D in sorted(lst):
        print(f"  {D:<7d} {nd(x(D, 4)[1], 14)}")


# ---------------------------------------------------------------------------------------------------------------
# E1, E2, E3 at k = 7 and the steps that remove Y_5 (E4), L_4 (E2W) and Y_4 (E5), on the hierarchy (no closed form, no matching beyond method B)
# ---------------------------------------------------------------------------------------------------------------
EGRID = (40, (1, 2), 18, 100)

# Generators on theta in (0, pi/2): ell = ln cos th, u = tan th, ui = cot th, th, C = Cl_2(pi - 2 th), D3 = Cl_3(pi - 2 th) + 3 zeta(3)/4, Lc2 = int_0^th ell^2.
# Their derivatives close: ell' = -u, u' = 1 + u^2, ui' = -(1 + ui^2), th' = 1, C' = 2 ell + 2 q, D3' = 2 C, Lc2' = ell^2 (q = ln 2), Lc1 = (C - 2 q th)/2.
X_LG, X_U, X_UI, X_TH, X_C, X_D3, X_LC = sp.symbols("ell u ui th C D3 Lc2")
X_Q = sp.Symbol("q", positive=True)
X_W = sp.symbols("w0:3")
X_J2 = sp.Symbol("J2")


def x_dth(e):
    """d/dtheta on the generators"""
    return (sp.diff(e, X_LG) * (-X_U) + sp.diff(e, X_U) * (1 + X_U**2) + sp.diff(e, X_TH) + sp.diff(e, X_C) * (2 * X_LG + 2 * X_Q)
            + sp.diff(e, X_D3) * 2 * X_C + sp.diff(e, X_LC) * X_LG**2 + sp.diff(e, X_UI) * (-(1 + X_UI**2)))


def x_pieces():
    """the weights and primitives of the steps E4, E2W, E5 in the generators and w0, w1, w2 = Omega0, Omega1, Omega2 (sympy)"""
    lg, u, ui, th, C, D3, Lc2, q = X_LG, X_U, X_UI, X_TH, X_C, X_D3, X_LC, X_Q
    w0, w1, w2 = X_W
    S = u**2 / (1 + u**2)
    Lc1 = (C - 2 * q * th) / 2
    c1 = w0 - 3 * lg - S
    c2 = w1 - 2 * w0 * lg + lg**2 - 2 * S * lg
    c3 = 2 * w0 * lg**2 - 2 * w1 * lg + w2 - 2 * lg**3 - 2 * S * lg**2
    g = (20 if MUTATE == "e4_g" else 22) * lg**2 - 8 * w0 * lg - 2 * w1       # E4: int Y_5 w5 = -2 int (g + cC u C) F_4, w5 = (8 ell c1 + 2 c2 + 4 w0 ell)/S
    cC = 0 if MUTATE == "e4_noC" else 16
    w5 = (8 * lg * c1 + 2 * c2 + 4 * w0 * lg) / S
    g4 = 12 * lg**2 - 2 * g - 2 * cC * u * C                                    # L_4's weight (E2 and the collapse of L_5, then E4)
    # P' = g4 u, P(0) = 0: int u^2 C = u C + ell^2 + 2 q ell - D3/2
    P = -sp.integrate(sp.expand(12 * lg**2 - 2 * g), lg) - 2 * cC * (u * C + lg**2 + 2 * q * lg - (0 if MUTATE == "e2w_P" else 1) * D3 / 2)
    w4p = (8 * lg * c2 + 2 * c3 + 8 * w0 * lg * c1 - P) / S                      # Y_4's weight
    # W' = w4p: S w4p = pol(ell) + 32 u C - 16 D3 - S (20 ell^2 + 8 w0 ell); int ell^b/S = -ell^b ui - b Lc_{b-1}; int u C/S = Z + int C tan,
    # int D3/S = -D3 ui + 2 Z (Z = int C cot th cancels); int C tan = -C ell + 2 Lc2 + 2 q Lc1
    a1, a2, a3 = 8 * w1 + 8 * w0**2 + 64 * q, 32 - 28 * w0, -sp.Rational(20, 3)
    Wn = (-2 * w2 * ui + a1 * (-lg * ui - th) + a2 * (-lg**2 * ui - 2 * Lc1) + a3 * (-lg**3 * ui - 3 * Lc2)
          + 32 * (-C * lg + 2 * Lc2 + 2 * q * Lc1) + (0 if MUTATE == "e5_noD3" else 16) * D3 * ui - 20 * Lc2 - 8 * w0 * Lc1)
    lim = sp.expand(Wn.subs({ui: 0, C: 0, D3: 0}).subs({th: sp.pi / 2, Lc2: X_J2}))      # W(pi/2) of the primitive above
    return dict(S=S, g=g, cC=cC, w5=w5, g4=g4, P=P, w4p=w4p, Wn=Wn, lim=lim)


_XP = {}


def x_funcs():
    """lambdified weights: name -> function of (ell, u, ui, th, C, D3, Lc2, q, w0, w1, w2, J2)"""
    if not _XP:
        pc = x_pieces()
        gens = (X_LG, X_U, X_UI, X_TH, X_C, X_D3, X_LC, X_Q) + tuple(X_W) + (X_J2,)
        for k in ("g", "w5", "g4", "P", "w4p", "Wn", "lim"):
            _XP[k] = sp.lambdify(gens, pc[k], "mpmath")
        _XP["_pieces"] = pc
    return _XP


class Fields:
    """Y_0..Y_k, L_0..L_k, G_k of the hierarchy for one Omega list, with the nodes' l, R, ell, sin^2, tan, the generators and the th-quadrature"""

    def __init__(self, B, Om, k):
        self.B, self.Om, self.k = B, list(Om), k
        self.N, self.P = B.N, B.P
        self.Y, self.L, self.G = B.hierarchy(self.Om, k)
        self.l, self.R, self.u = B.l, B.R, B.u
        self.S = self.grid(lambda p, i: B.Rm1[p][i] / B.R[p][i])                 # sin^2 th = 1 - 1/R
        self.ell = self.grid(lambda p, i: -B.l[p][i] / 2)                        # ln cos th = -l/2
        self.W = self.grid(lambda p, i: B.half * B.wend[i] * B.v[p][i] / B.u[p][i])   # th-quadrature weights: dth = (v/u) dv
        cs = self.grid(lambda p, i: mp.exp(-B.l[p][i] / 2))                      # cos th = sin x
        self.th = self.grid(lambda p, i: mp.acos(cs[p][i]))
        self.ui = self.grid(lambda p, i: 1 / B.u[p][i])
        self.C = self.grid(lambda p, i: mp.clsin(2, 2 * mp.asin(cs[p][i])))      # Cl_2(pi - 2 th) = Cl_2(2 x), x = pi/2 - th
        self.D3 = self.grid(lambda p, i: mp.clcos(3, 2 * mp.asin(cs[p][i])) + 3 * mp.zeta(3) / 4)
        self.Lc2 = B.cum_left(self.grid(lambda p, i: self.ell[p][i] ** 2 * B.v[p][i] / B.u[p][i]))      # int_0^th ell^2 dth, by the panels' own integration

    def grid(self, fn):
        return [[fn(p, i) for i in range(self.N)] for p in range(self.P)]

    def integral(self, f):
        """int_0^{pi/2} f dth, f(p, i)"""
        return sum(self.W[p][i] * f(p, i) for p in range(self.P) for i in range(self.N))

    def h(self, m, p, i):
        return hk(m, self.Om, self.l[p][i], self.R[p][i], [y[p][i] for y in self.Y], [lv[p][i] for lv in self.L] if m else None)

    def F(self, j, p, i):
        """F_j = sum_{q <= j} L_{j-q} l^q / q!"""
        return sum(self.L[j - q][p][i] * self.l[p][i] ** q / mp.factorial(q) for q in range(j + 1))

    def Fpart(self, j, lo, hi, p, i):
        """sum_{q = lo}^{hi} L_{j-q} l^q / q!"""
        return sum(self.L[j - q][p][i] * self.l[p][i] ** q / mp.factorial(q) for q in range(lo, hi + 1))

    def Sk(self, k, p, i, upto=None):
        """sum_{j < upto} Y_j c_{k-j}"""
        upto = k if upto is None else upto
        return sum(self.Y[j][p][i] * cc(k - j, self.Om, self.l[p][i], self.R[p][i]) for j in range(upto))

    def gens(self, p, i):
        om = self.Om
        return (self.ell[p][i], self.u[p][i], self.ui[p][i], self.th[p][i], self.C[p][i], self.D3[p][i], self.Lc2[p][i], mp.ln2, om[0], om[1], om[2],
                mp.pi * (12 * mp.ln2**2 + mp.pi**2) / 24)


def e1_sides(fl, k):
    om = fl.Om
    A = 4 * fl.integral(lambda p, i: fl.L[k - 1][p][i])
    Bs = 2 * fl.integral(lambda p, i: fl.F(k - 1, p, i) - fl.L[k - 1][p][i])
    if MUTATE == "e1_sign":
        Bs = -Bs
    Cs = 2 * fl.integral(lambda p, i: (fl.Sk(k, p, i) - om[k - 1]) / fl.S[p][i])
    return fl.G[k], A + Bs + Cs, (A, Bs, Cs)


def e2_sides(fl, m):
    lhs = fl.integral(lambda p, i: fl.L[m][p][i])
    shift = mp.ln2 if MUTATE == "e2_weight" else 0
    rhs = 2 * fl.integral(lambda p, i: fl.h(m, p, i) * (fl.ell[p][i] + shift) / fl.S[p][i])
    return lhs, rhs


def e3_sides(fl, m):
    """E3: int Y_m [(ell - Omega0)/sin^2 + 1] against 2 int (ell - Omega0) F_{m-1}"""
    O0 = fl.Om[0]
    fac = 1 if MUTATE == "e3_factor" else 2
    lhs = fl.integral(lambda p, i: fl.Y[m][p][i] * ((fl.ell[p][i] - O0) / fl.S[p][i] + 1))
    rhs = fac * fl.integral(lambda p, i: (fl.ell[p][i] - O0) * fl.F(m - 1, p, i))
    return lhs, rhs


def e123_sides(fl):
    """G_7 against E1, E2 (m = 6) and E3 (m = 6) together: 8 int H_6 ell/S + 2 int (F_6 - L_6) + 2 int (S'_7 - Omega6)/S - 4 int (ell - Omega0) F_5"""
    om, O0 = fl.Om, fl.Om[0]
    ell, S, l = fl.ell, fl.S, fl.l
    shift = mp.ln2 if MUTATE == "e2_weight" else 0
    f3 = 2 if MUTATE == "e3_factor" else 4
    H6 = lambda p, i: S[p][i] * fl.L[5][p][i] + fl.Sk(6, p, i, upto=6)  # noqa: E731
    t = [8 * fl.integral(lambda p, i: (ell[p][i] + shift) * H6(p, i) / S[p][i]),
         (-2 if MUTATE == "e1_sign" else 2) * fl.integral(lambda p, i: fl.Fpart(6, 1, 6, p, i)),
         2 * fl.integral(lambda p, i: (fl.Sk(7, p, i, upto=6) - om[6]) / S[p][i]),
         -f3 * fl.integral(lambda p, i: (ell[p][i] - O0) * fl.F(5, p, i))]
    return fl.G[7], sum(t), max(abs(x) for x in t)


def e4_sides(fl):
    """E4 (the Y_5 step at k = 7): int Y_5 w5 = -2 int (g + 16 u C) F_4, w5 = (2 Omega1 + 8 Omega0 ell - 22 ell^2)/S - 12 ell"""
    xf = x_funcs()
    S = fl.S
    lhs = fl.integral(lambda p, i: fl.Y[5][p][i] * xf["w5"](*fl.gens(p, i)))
    rhs = -2 * fl.integral(lambda p, i: (xf["g"](*fl.gens(p, i)) + xf["_pieces"]["cC"] * fl.u[p][i] * fl.C[p][i]) * fl.F(4, p, i))
    return lhs, rhs


def e2w_sides(fl, part):
    """E2W: int L_4 g = -2 int P h_4/S for the weight g = g4' of P (a: g = 12 ell^2, P = -4 ell^3 as a primitive of g u; b: the weight after E4)"""
    xf = x_funcs()
    if part == "a":
        lhs = fl.integral(lambda p, i: fl.L[4][p][i] * (10 if MUTATE == "e2w_g" else 12) * fl.ell[p][i] ** 2)
        rhs = 2 * fl.integral(lambda p, i: 4 * fl.ell[p][i] ** 3 * fl.h(4, p, i) / fl.S[p][i])
    else:
        lhs = fl.integral(lambda p, i: fl.L[4][p][i] * xf["g4"](*fl.gens(p, i)))
        rhs = -2 * fl.integral(lambda p, i: xf["P"](*fl.gens(p, i)) * fl.h(4, p, i) / fl.S[p][i])
    return lhs, rhs


def Wfun(fl, p, i):
    xf = x_funcs()
    g = fl.gens(p, i)
    return xf["Wn"](*g) - (0 if MUTATE == "e5_const" else xf["lim"](*g))


def e5_sides(fl):
    """E5 (the Y_4 step): int Y_4 w4' = -2 int W u F_3, W' = w4', W(pi/2) = 0"""
    xf = x_funcs()
    lhs = fl.integral(lambda p, i: fl.Y[4][p][i] * xf["w4p"](*fl.gens(p, i)))
    rhs = -2 * fl.integral(lambda p, i: Wfun(fl, p, i) * fl.u[p][i] * fl.F(3, p, i))
    return lhs, rhs


def g7_reduced(fl):
    """G_7 with no Y_4 .. Y_6, L_4 .. L_6: E1, E2 (m = 6, 5), E3 (m = 6), E4, E2W and E5, on the hierarchy's own Y_0 .. Y_3, L_0 .. L_3"""
    om, O0 = fl.Om, fl.Om[0]
    ell, S, l = fl.ell, fl.S, fl.l
    xf = x_funcs()
    shift = mp.ln2 if MUTATE == "e2_weight" else 0
    f3 = 2 if MUTATE == "e3_factor" else 4
    t = []
    t.append(8 * fl.integral(lambda p, i: (ell[p][i] + shift) * fl.Sk(6, p, i, upto=4) / S[p][i]))
    t.append((-2 if MUTATE == "e1_sign" else 2) * fl.integral(lambda p, i: fl.Fpart(6, 3, 6, p, i)))
    t.append(2 * fl.integral(lambda p, i: (fl.Sk(7, p, i, upto=4) - om[6]) / S[p][i]))
    t.append(-f3 * fl.integral(lambda p, i: (ell[p][i] - O0) * fl.Fpart(5, 2, 5, p, i)))
    t.append(8 * O0 * fl.integral(lambda p, i: (ell[p][i] + shift) * fl.Sk(5, p, i, upto=4) / S[p][i]))
    t.append(-2 * fl.integral(lambda p, i: (xf["g"](*fl.gens(p, i)) + xf["_pieces"]["cC"] * fl.u[p][i] * fl.C[p][i]) * fl.Fpart(4, 1, 4, p, i)))
    t.append(-2 * fl.integral(lambda p, i: xf["P"](*fl.gens(p, i)) * (S[p][i] * fl.L[3][p][i] + fl.Sk(4, p, i, upto=4)) / S[p][i]))
    t.append(-2 * fl.integral(lambda p, i: Wfun(fl, p, i) * fl.u[p][i] * fl.F(3, p, i)))
    return sum(t), t


def check_exch7():
    t0 = time.time()
    N, w, V, dps = EGRID
    mp.mp.dps = 60
    with mp.workdps(dps):
        B = PanelK(N, mp.mpf(w[0]) / w[1], V)
        Om, _, _, _ = omega_chain(B, 7)
        print(f"  E-grid N = {N}, w = {w[0]}/{w[1]}, V = {V}, dps {dps}, {B.P} panels; its Omega0 .. Omega4 meet the closed forms to "
              + ", ".join(f"{float(digits(c, o)):.1f}" for c, o in zip(closed_omegas(), Om)) + " digits", flush=True)
        res = {}
        for label, o5v in (("root", Om[5]), ("0", mp.mpf(0)), ("1/2", mp.mpf(1) / 2)):
            Omv = list(Om)
            Omv[5] = o5v
            fl = Fields(B, Omv, 7)
            g7, rhs, parts = e1_sides(fl, 7)
            scale = max(abs(x) for x in parts)
            e1 = abs(g7 - rhs) / scale
            lhs2, rhs2 = e2_sides(fl, 6)
            lhs2b, rhs2b = e2_sides(fl, 5)
            lhs2c, rhs2c = e2_sides(fl, 4)
            lhs3, rhs3 = e3_sides(fl, 6)
            g123, r123, sc123 = e123_sides(fl)
            lhs4, rhs4 = e4_sides(fl)
            lw_a, rw_a = e2w_sides(fl, "a")
            lw_b, rw_b = e2w_sides(fl, "b")
            lhs5, rhs5 = e5_sides(fl)
            gr, grp = g7_reduced(fl)
            sc = max(abs(x) for x in grp)
            j2 = digits(fl.integral(lambda p, i: fl.ell[p][i] ** 2), mp.pi * (12 * mp.ln2**2 + mp.pi**2) / 24)
            res[label] = dict(
                E1=(-mp.log10(e1) if e1 else mp.mpf(dps), g7, scale),
                E2_6=(digits(lhs2, rhs2), lhs2, rhs2), E2_5=(digits(lhs2b, rhs2b), lhs2b, rhs2b), E2_4=(digits(lhs2c, rhs2c), lhs2c, rhs2c),
                E3_7=(digits(lhs3, rhs3), lhs3, rhs3), E123=(-mp.log10(abs(g123 - r123) / sc123) if g123 != r123 else mp.mpf(dps), g123, r123),
                E4_7=(digits(lhs4, rhs4), lhs4, rhs4), E2Wa=(digits(lw_a, rw_a), lw_a, rw_a), E2Wb=(digits(lw_b, rw_b), lw_b, rw_b),
                E5_7=(digits(lhs5, rhs5), lhs5, rhs5), G7R=(-mp.log10(abs(g7 - gr) / sc) if g7 != gr else mp.mpf(dps), g7, gr), LC2=(j2, 0, 0))
        for key, cid, what in (("E1", "E1_7", "G_7 against E1's right side (relative to the largest of its three terms)"),
                               ("E2_6", "E2_6", "int L_6 = 2 int h_6 ell/sin^2"), ("E2_5", "E2_5", "int L_5 = 2 int h_5 ell/sin^2"),
                               ("E2_4", "E2_4", "int L_4 = 2 int h_4 ell/sin^2"),
                               ("E3_7", "E3_7", "int Y_6 [(ell - Omega0)/sin^2 + 1] = 2 int (ell - Omega0) F_5"),
                               ("E123", "E123_7", "G_7 against E1 + E2 (m = 6) + E3 (m = 6): 8 int H_6 ell/S + 2 int (F_6 - L_6) + 2 int (S'_7 - Omega6)/S - 4 int (ell - Omega0) F_5 (relative to its largest term)"),
                               ("E4_7", "E4_7", "int Y_5 w5 = -2 int (g + 16 u C) F_4 (the order-6 step E4 at one level higher)"),
                               ("E2Wa", "E2Wa", "12 int L_4 ell^2 = 8 int ell^3 h_4/S (E2 with a weight)"),
                               ("E2Wb", "E2Wb", "int L_4 g4 = -2 int P h_4/S, g4 = 12 ell^2 - 2 g - 32 u C, P' = g4 tan th, P(0) = 0"),
                               ("E5_7", "E5_7", "int Y_4 w4' = -2 int W u F_3, W' = w4', W(pi/2) = 0"),
                               ("G7R", "G7R", "G_7 against the reduced form (E1, E2, E3, E4, E2W, E5: Y_0 .. Y_3, L_0 .. L_3 only; relative to its largest term)"),
                               ("LC2", "LC2", "the nodes' int ell^2 dth against J_2 = pi (12 ln^2 2 + pi^2)/24")):
            ds = {lab: r[key][0] for lab, r in res.items()}
            report(cid, all(d >= 30 for d in ds.values()), f"{what}: " + ", ".join(f"{float(d):.1f} digits at Omega5 = {lab}" for lab, d in ds.items())
                   + " (bar 30)")
        print(f"    G_7 at Omega5 = 0 and 1/2: {nd(res['0']['E1'][1], 20)}, {nd(res['1/2']['E1'][1], 20)} ({time.time() - t0:.0f} s)")


def check_xsym():
    """XS: the primitives of the steps, symbolically: P' = g4 u, W' = w4' (with ui = 1/u), and the limit W(pi/2) of W's primitive"""
    pc = x_pieces()
    r1 = sp.simplify(x_dth(pc["P"]) - pc["g4"] * X_U)
    r2 = sp.simplify(x_dth(pc["Wn"]).subs(X_UI, 1 / X_U) - pc["w4p"])
    r3 = sp.simplify(sp.expand(pc["w5"] - ((2 * X_W[1] + 8 * X_W[0] * X_LG - 22 * X_LG**2) / pc["S"] - 12 * X_LG)))
    report("XS", r1 == 0 and r2 == 0 and r3 == 0, f"sympy: P' - g4 u = {r1}; W' - w4' = {r2}; w5 - w_4(order 6) = {r3}; W(pi/2) of the primitive = {pc['lim']}")



# ---------------------------------------------------------------------------------------------------------------
# the reduced integrand of G_7 (after E1 .. E5, from the closed forms of Y_0 .. Y_3, L_0 .. L_3): its decomposition (RD1) and method A (A7)
# ---------------------------------------------------------------------------------------------------------------
def r7_terms(Om, gen, Y, L, l, R, S, u, th, C, D3, Lc2, lg, q, num):
    """the eight terms t1 .. t8 of G_7 (see g7_reduced) from Y_0 .. Y_3, L_0 .. L_3: one function for the symbolic decomposition and the numeric integrand;
    gen(name) gives the weights g, P, W - lim (functions of the generators)"""
    O0 = Om[0]
    fa = math.factorial

    def Fp(j, lo, hi):
        return sum(L[j - i] * l**i / fa(i) for i in range(lo, hi + 1))

    def Sk(k, upto, om=Om):
        return sum(Y[j] * cc(k - j, om, l, R) for j in range(upto))

    g, P, Wd, cC = gen("g"), gen("P"), gen("Wd"), gen("cC")
    t = [8 * lg * Sk(6, 4) / S, 2 * Fp(6, 3, 6), 2 * (Sk(7, 4) - Om[6]) / S, -4 * (lg - O0) * Fp(5, 2, 5), 8 * O0 * lg * Sk(5, 4) / S,
         -2 * (g + cC * u * C) * Fp(4, 1, 4), -2 * P * (S * L[3] + Sk(4, 4)) / S, -2 * Wd * u * Fp(3, 0, 3)]
    return t


def r7_parts():
    """the integrand of G_7 as coefficients of C^a D3^b Lc2^c (u-Laurent polynomials in ell, th, q, w1 .. w5): parts[(a, b, c)]; Omega0 = -4 ln 2"""
    lg, u, ui, th, C, D3, Lc2, q = X_LG, X_U, X_UI, X_TH, X_C, X_D3, X_LC, X_Q
    w = sp.symbols("v0:6")
    S = sp.Symbol("S")
    R = 1 / (1 - S)
    l = -2 * lg
    Om = list(w) + [sp.Integer(0)]
    Y = [sp.Integer(1), 2 * lg, 4 * w[0] * lg - 6 * lg**2, 4 * w[1] * lg - 4 * lg**3 + 32 * lg**2 - 16 * w[0] * lg + 32 * u * C - 16 * D3]
    L = [sp.Integer(-1), -2 * w[0] + 4 * lg, -2 * w[1] - 4 * w[0] * lg + 16 * lg**2 + 16 * u * C]
    h3 = hk(3, Om, l, R, Y, L)
    sc = u * (1 - S)
    B = -32 * (lg + q) * C + 16 * S * C + 32 * Lc2 + 16 * sc * lg * (lg + 2 * q) - 4 * w[1] * th - 2 * w[1] * sc
    L = L + [-2 * h3 + 2 * u * B]
    pc = x_pieces()
    rep = {X_W[0]: w[0], X_W[1]: w[1], X_W[2]: w[2]}
    pcs = {k: pc[k].subs(rep) for k in ("g", "P", "Wn", "lim")}
    lim = pcs["lim"].subs(X_J2, sp.Symbol("J2"))
    gen = {"g": pcs["g"], "P": pcs["P"], "Wd": (pcs["Wn"] - lim).subs(X_UI, 1 / u), "cC": pc["cC"]}
    t = r7_terms(Om, lambda name: gen[name], Y, L, l, R, S, u, th, C, D3, Lc2, lg, q, None)
    tot = sp.expand(sum(t).subs(w[0], -4 * q))
    tot = sp.expand(sp.cancel(tot.subs(S, u**2 / (1 + u**2))))
    Pl = sp.Poly(tot, C, D3, Lc2)
    return {mon: sp.expand(sp.cancel(cf)) for mon, cf in Pl.terms()}, tot


def r7_numeric_integrand(Om):
    """the integrand of G_7 after E1 .. E5, from the closed forms of Y_0 .. Y_3, L_0 .. L_3 (order4.py, order6.py's L_3), Om = [w0 .. w4, w5, 0]"""
    q = mp.ln2
    xf = x_funcs()
    J2 = mp.pi * (12 * q**2 + mp.pi**2) / 24
    lgs = x_pieces()

    def f(th, x, tail):
        if tail:
            th = mp.pi / 2 - x
        c, lc, l, u, R, s2, Y, L, K, T = o4.closed(th, x, tail, Om)
        ell, s2b, ub, Cn, D3n = o4.basics(th, x, tail)
        Lc2 = o4.Lc2_poly(th)
        sc = mp.sin(th) * mp.cos(th)
        B = -32 * (lc + q) * Cn + 16 * s2 * Cn + 32 * Lc2 + 16 * sc * lc * (lc + 2 * q) - 4 * Om[1] * th - 2 * Om[1] * sc
        LL = [L[0], L[1], L[2], -2 * hk(3, Om, l, R, Y, L) + 2 * u * B]
        gens = (lc, u, 1 / u, th, Cn, D3n, Lc2, q, Om[0], Om[1], Om[2], J2)
        W = {"g": xf["g"](*gens), "P": xf["P"](*gens), "Wd": xf["Wn"](*gens) - xf["lim"](*gens), "cC": lgs["cC"]}
        t = r7_terms(Om, lambda name: W[name], Y, LL, l, R, s2, u, th, Cn, D3n, Lc2, lc, q, None)
        return sum(t)

    return o4.guarded(f, 150)


# ---------------------------------------------------------------------------------------------------------------
# the finite-part reduction of the integrand of G_7 to master integrals (RD)
# ---------------------------------------------------------------------------------------------------------------
# Generators on theta in (0, pi/2), in lt = ln(2 cos th) = ell + ln 2 (so that the masters have no ln 2): lt' = -u, u' = 1 + u^2, ui = cot th, C' = 2 lt, D3' = 2 C,
# K2 = int_0^th lt^2 (K2' = lt^2), K3' = lt^3. <f> is the finite part at theta -> pi/2 (cutoff x = pi/2 - th: the constant term of x^0 ln^0 x) and at theta -> 0 of int f;
# for an integrable f it is the integral, <.> is linear, and <F'> = B[F] = FP F(pi/2) - FP F(0). One step is I <- I - dF, const += B[F]; the integral of the final I is a sum of masters.
RS = types.SimpleNamespace(**{n: sp.Symbol(n) for n in ("lt", "u", "ui", "th", "C", "D3", "K2", "K3", "Jt2", "Jt3", "v5")})
RS_X, RS_L, RS_T, RS_S, RS_R = sp.symbols("x L t s r", positive=True)
RS_ORD, RS_SHIFT = 10, 9
RD_STATE = {}


def rd_dth(e):
    lt, u, ui, th, C, D3, K2, K3 = (RS.lt, RS.u, RS.ui, RS.th, RS.C, RS.D3, RS.K2, RS.K3)
    return (sp.diff(e, lt) * (-u) + sp.diff(e, u) * (1 + u**2) + sp.diff(e, ui) * (-(1 + ui**2)) + sp.diff(e, th) + sp.diff(e, C) * 2 * lt
            + sp.diff(e, D3) * 2 * C + sp.diff(e, K2) * lt**2 + sp.diff(e, K3) * lt**3)


def rd_kit():
    """the asymptotic series of the generators at theta -> pi/2 (x = pi/2 - th, L = ln x) and at theta -> 0, for the finite parts"""
    x, L, t, s, r, q = RS_X, RS_L, RS_T, RS_S, RS_R, X_Q
    z3 = sp.zeta(3)
    logs = lambda e: sp.expand(e).subs(sp.log(x), L)  # noqa: E731
    ell1 = sp.series(sp.log(sp.sin(s) / s), s, 0, RS_ORD).removeO()
    sgn = -1 if MUTATE == "rd_kit" else 1
    ell_s = sp.log(s) + q + ell1                                      # ln(2 sin x)
    C_s = sp.expand(-sp.integrate(sp.expand(2 * (sp.log(s) + sgn * q + ell1)), (s, 0, s)))
    D3_s = sp.expand(sp.Rational(7, 4) * z3 - 2 * sp.integrate(C_s.subs(s, r), (r, 0, s)))
    K2_s = sp.expand(RS.Jt2 - sp.integrate(sp.expand(ell_s.subs(s, r) ** 2), (r, 0, s)))
    K3_s = sp.expand(RS.Jt3 - sp.integrate(sp.expand(ell_s.subs(s, r) ** 3), (r, 0, s)))
    pi2 = {"lt": L + q + ell1.subs(s, x), "u": sp.series(sp.cot(x), x, 0, RS_ORD).removeO(), "ui": sp.series(sp.tan(x), x, 0, RS_ORD).removeO(),
           "th": sp.pi / 2 - x, "C": logs(C_s.subs(s, x)), "D3": logs(D3_s.subs(s, x)), "K2": logs(K2_s.subs(s, x)), "K3": logs(K3_s.subs(s, x))}
    lt_t = q + sp.series(sp.log(sp.cos(t)), t, 0, RS_ORD).removeO()
    C_t = sp.expand(sp.integrate(sp.expand(2 * lt_t), (t, 0, t)))
    D3_t = sp.expand(sp.integrate(sp.expand(2 * C_t.subs(t, r)), (r, 0, t)))
    K2_t = sp.expand(sp.integrate(sp.expand(lt_t.subs(t, r) ** 2), (r, 0, t)))
    K3_t = sp.expand(sp.integrate(sp.expand(lt_t.subs(t, r) ** 3), (r, 0, t)))
    zero = {"lt": lt_t, "u": sp.series(sp.tan(t), t, 0, RS_ORD).removeO(), "ui": sp.series(sp.cot(t), t, 0, RS_ORD).removeO(), "th": t, "C": C_t, "D3": D3_t,
            "K2": K2_t, "K3": K3_t}
    return pi2, zero


def rd_const_term(e, v, logsym=None):
    e = sp.expand(e * v**RS_SHIFT)
    gens = (v, logsym) if logsym is not None else (v,)
    P = sp.Poly(e, *gens)
    return P.coeff_monomial((RS_SHIFT,) + ((0,) if logsym is not None else ()))


def rd_bnd(expr):
    """B[F] = FP_{pi/2} F - FP_0 F"""
    if "kit" not in RD_STATE:
        RD_STATE["kit"] = rd_kit()
    pi2, zero = RD_STATE["kit"]
    expr = sp.expand(expr).replace(lambda e: e.is_Pow and e.base == RS.u and e.exp < 0, lambda e: RS.ui ** (-e.exp))
    names = {"lt": RS.lt, "u": RS.u, "ui": RS.ui, "th": RS.th, "C": RS.C, "D3": RS.D3, "K2": RS.K2, "K3": RS.K3}
    syms = {k: v for k, v in names.items() if expr.has(v)}
    a = rd_const_term(expr.subs({syms[k]: pi2[k] for k in syms}, simultaneous=True), RS_X, RS_L)
    b = rd_const_term(expr.subs({syms[k]: zero[k] for k in syms}, simultaneous=True), RS_T)
    return sp.expand(a - b)


def rd_integrand():
    """the integrand of G_7 after E1 .. E5 in the generators lt, u, th, C, D3, K2 (Omega0 = -4 ln 2, Omega1 .. Omega4 at their closed forms, Omega5 = v5)"""
    parts, tot = r7_parts()
    q, pi, z3, z5 = X_Q, sp.pi, sp.zeta(3), sp.zeta(5)
    v = sp.symbols("v0:6")
    om = {v[1]: 2 * pi**2 / 3 + 8 * q**2, v[2]: 80 * z3 - sp.Rational(8, 3) * pi**2 * q - sp.Rational(32, 3) * q**3,
          v[3]: 22 * pi**4 / 45 + sp.Rational(16, 3) * pi**2 * q**2 + sp.Rational(32, 3) * q**4 - 320 * z3 * q - 384 * z3,
          v[4]: -10432 * z5 - 160 * pi**2 * z3 + 640 * q**2 * z3 - sp.Rational(64, 9) * pi**2 * q**3 - sp.Rational(88, 45) * pi**4 * q
          - sp.Rational(128, 15) * q**5 - sp.Rational(2432, 45) * pi**4 + 1536 * q * z3 + 3072 * z3}
    e = tot.subs(om).subs(X_J2, pi**3 / 24 + q**2 * pi / 2).subs(v[5], RS.v5)
    e = e.subs(X_LC, RS.K2 - q * X_C + q**2 * X_TH).subs(X_LG, RS.lt - q)
    return sp.expand(e)


def rd_chain(I):
    """the integration by parts as a list of (name, F); returns (I_final, Fs, notes)"""
    lt, u, ui, th, C, D3, K2, K3 = (RS.lt, RS.u, RS.ui, RS.th, RS.C, RS.D3, RS.K2, RS.K3)
    st = {"I": sp.expand(I), "Fs": []}
    notes = {}

    def coeff(a, b, c):
        return sp.expand(sp.Poly(st["I"], C, D3, K2).coeff_monomial(C**a * D3**b * K2**c))

    def apply(name, F):
        dF = sp.expand(rd_dth(F).subs(ui, 1 / u) if F.has(ui) else rd_dth(F))
        st["I"] = sp.expand(st["I"] - dF)
        st["Fs"].append((name, F))

    def P_b(b):
        return [th, C / 2, K2, RS.K3][b]

    def primitive_c(c):
        """the primitive in th of c = alpha(lt) + beta(lt) u^-2: lt^b -> P_b, lt^b (1 + ui^2) -> -lt^b ui - b P_{b-1}"""
        cu = sp.collect(sp.expand(c), u, evaluate=False)
        assert set(cu) <= {sp.Integer(1), u**-2}, list(cu)
        alpha, beta = sp.expand(cu.get(sp.Integer(1), 0)), sp.expand(cu.get(u**-2, 0))
        c0, c1 = sp.expand(alpha - beta), beta
        A = 0
        sg = 1 if MUTATE == "rd_Dstep" else -1
        if c0 != 0:
            pc0 = sp.Poly(c0, lt)
            A += sum(pc0.coeff_monomial(lt**b) * P_b(b) for b in range(pc0.degree() + 1))
        if c1 != 0:
            pc1 = sp.Poly(c1, lt)
            A += sum(pc1.coeff_monomial(lt**b) * (-lt**b * ui + sg * (b * P_b(b - 1) if b else 0)) for b in range(pc1.degree() + 1))
        return sp.expand(A), sp.simplify(sp.expand(rd_dth(sp.expand(A)).subs(ui, 1 / u) - c))

    def solve_ft(g):
        """f with f - f' = g (a polynomial in lt)"""
        f, d = 0, sp.expand(g)
        while d != 0:
            f += d
            d = sp.diff(d, lt)
        return sp.expand(f)

    # the quadratic terms
    a1 = sp.expand(coeff(0, 0, 2) / u**2)
    assert not a1.has(u)
    apply("LL", a1 * u * K2**2)                                                  # <u^2 K2^2> = B[u K2^2] - 2<u lt^2 K2> - <K2^2>
    a2 = sp.expand(coeff(0, 1, 1) / u)
    assert not a2.has(u)
    apply("DL", -a2 * lt * D3 * K2)                                              # <u D3 K2> = B[-lt D3 K2] + 2<lt C K2> + <lt^3 D3>
    A, res = primitive_c(coeff(0, 2, 0))
    notes["DD primitive residual"] = res
    apply("DD", A * D3**2)                                                       # <c D3^2> = B[A D3^2] - 4<A D3 C>, A = 512 (3 lt + q) cot
    cCD = sp.collect(coeff(1, 1, 0), u, evaluate=False)
    notes["CD: 1/u part after DD"] = sp.expand(cCD.get(1 / u, 0))
    G = sp.integrate(sp.expand(cCD.get(u, 0)), lt)
    apply("CD", -G * C * D3)                                                     # <u g C D3> = B[-G C D3] + <G (2 lt D3 + 2 C^2)>
    cu = sp.collect(coeff(1, 0, 1), u, evaluate=False)
    gt = solve_ft(cu.get(u**2, 0))
    apply("CL", u * gt * C * K2)                                                 # <u^2 (gt - gt') C K2> = B[u gt C K2] - <gt C K2> - ...
    cu = sp.collect(coeff(2, 0, 0), u, evaluate=False)
    f = solve_ft(cu.get(u**2, 0))
    apply("CC", u * f * C**2)                                                    # <u^2 (f - f') C^2> = B[u f C^2] - <f C^2> - 4<u f lt C>
    # the D3 coefficient: its u theta part, then the primitive of the rest
    cu = sp.collect(coeff(0, 1, 0), u, evaluate=False)
    apply("D_u", -lt * sp.expand(cu.get(u, 0)) * D3)
    A, res = primitive_c(coeff(0, 1, 0))
    notes["D primitive residual"] = res
    apply("D", A * D3)                                                           # <c D3> = B[A D3] - 2<A C>
    cK3C = sp.expand(sp.expand(st["I"]).coeff(RS.K3, 1)).coeff(C, 1)
    notes["K3 apart from C K3"] = sp.expand(sp.expand(st["I"]).coeff(RS.K3, 1) - cK3C * C)
    apply("K3C", sp.Rational(1, 2) * cK3C * D3 * RS.K3)                          # <C K3> = B[D3 K3/2] - <lt^3 D3>/2

    def uC():
        c1 = sp.expand(sp.collect(coeff(1, 0, 0), u, evaluate=False).get(u, 0))
        if c1 != 0:
            pc = sp.Poly(c1, lt)
            for b in range(pc.degree() + 1):
                if pc.coeff_monomial(lt**b) != 0:
                    apply(f"C_u_lt{b}", -pc.coeff_monomial(lt**b) * lt**(b + 1) / (b + 1) * C)        # <u lt^b C> = B[-lt^{b+1} C/(b+1)] + 2<lt^{b+2}>/(b+1)
    notes["C: 1/u part after D"] = sp.expand(sp.collect(coeff(1, 0, 0), u, evaluate=False).get(1 / u, 0))
    uC()
    c2 = sp.expand(sp.collect(coeff(1, 0, 0), u, evaluate=False).get(u**2, 0))        # u^2 theta g(lt) C: by parts, top degree first
    for b in range(sp.Poly(c2, lt).degree(), -1, -1):
        cb = sp.expand(sp.Poly(sp.expand(sp.collect(coeff(1, 0, 0), u, evaluate=False).get(u**2, 0)), lt).coeff_monomial(lt**b))
        if cb != 0:
            cth = sp.Poly(cb, th).coeff_monomial(th)
            assert sp.expand(cb - cth * th) == 0, cb
            apply(f"C_u2th_lt{b}", cth * u * th * lt**b * C)
    uC()
    # the K2 coefficient
    def uK():
        c1 = sp.expand(sp.collect(coeff(0, 0, 1), u, evaluate=False).get(u, 0))
        if c1 != 0:
            pc = sp.Poly(c1, lt)
            for b in range(pc.degree() + 1):
                if pc.coeff_monomial(lt**b) != 0:
                    apply(f"K_u_lt{b}", -pc.coeff_monomial(lt**b) * lt**(b + 1) / (b + 1) * K2)       # <u lt^b K2> = B[-lt^{b+1} K2/(b+1)] + <lt^{b+3}>/(b+1)
    uK()
    c2 = sp.expand(sp.collect(coeff(0, 0, 1), u, evaluate=False).get(u**2, 0))
    cth = sp.Poly(c2, th).coeff_monomial(th)
    assert sp.expand(c2 - cth * th) == 0
    apply("K_u2th", cth * u * th * K2)
    uK()
    # the pure terms
    c2 = sp.expand(sp.collect(coeff(0, 0, 0), u, evaluate=False).get(u**2, 0))
    if c2 != 0:
        assert sp.expand(c2 - sp.Poly(c2, th).coeff_monomial(th**2) * th**2) == 0
        apply("P_u2", sp.Poly(c2, th).coeff_monomial(th**2) * u * th**2)                  # <u^2 th^2> = B[u th^2] - 2<u th> - <th^2>
    st["I"] = sp.expand(st["I"].replace(lambda e: e.is_Pow and e.base == u and e.exp < 0, lambda e: ui ** (-e.exp)))
    pcm = sp.Poly(sp.expand(sp.collect(sp.expand(coeff(0, 0, 0)), ui, evaluate=False).get(ui**2, 0)), lt)
    for b in range(pcm.degree(), -1, -1):
        cbl = sp.Poly(sp.expand(sp.collect(sp.expand(coeff(0, 0, 0)), ui, evaluate=False).get(ui**2, 0)), lt).coeff_monomial(lt**b)
        if cbl != 0:
            F = -cbl * lt**b * ui                                                          # d(-lt^b ui) = lt^b (1 + ui^2) + b lt^(b-1)
            dF = sp.expand(rd_dth(F).replace(lambda e: e.is_Mul and e.has(u) and e.has(ui), lambda e: sp.simplify(e.subs(ui, 1 / u))))
            st["I"] = sp.expand(st["I"] - dF)
            st["Fs"].append((f"P_ui_lt{b}", F))
    c1 = sp.expand(sp.collect(sp.expand(coeff(0, 0, 0)), u, evaluate=False).get(u, 0))
    for (b, j), cf in sp.Poly(c1, lt, th).terms():
        apply(f"P_u_th{j}_lt{b}", -cf * th**j * lt**(b + 1) / (b + 1))                    # <u th^j lt^b> = B[-th^j lt^{b+1}/(b+1)] + ...
    I = sp.expand(st["I"])
    notes["u or 1/u left"] = I.has(u) or I.has(ui)
    return I, st["Fs"], notes


def rd_Jt():
    """Jt_n = int_0^{pi/2} ln^n(2 cos th) = n! [a^n] (pi/2) Gamma(a + 1)/Gamma(a/2 + 1)^2, n <= 7 (sympy, closed forms)"""
    a = sp.Symbol("a")
    ser = sp.series((sp.pi / 2) * sp.gamma(a + 1) / sp.gamma(a / 2 + 1)**2, a, 0, 8).removeO()
    Jt = [sp.simplify(sp.expand_func(sp.factorial(n) * ser.coeff(a, n))) for n in range(8)]
    if MUTATE == "rd_J5":
        Jt[5] = Jt[5].subs(sp.zeta(5), sp.Rational(17, 18) * sp.zeta(5))
    return Jt


def rd_masters(Jt):
    """closed forms of the master integrals <C^a D3^b K2^c th^j lt^k>, keyed (a, b, c, j, k); A = <lt^2 C^2> and B = <lt^3 D3> are identified numerically"""
    pi, z3, z5 = sp.pi, sp.zeta(3), sp.zeta(5)
    z2, z6 = pi**2 / 6, pi**6 / 945
    A = pi * (sp.Rational(3, 4) * z3**2 - sp.Rational(95, 96) * z6)
    B = pi * (sp.Rational(9, 16) * z3**2 - sp.Rational(11, 4) * z6)
    if MUTATE == "rd_master":
        A = pi * (sp.Rational(3, 4) * z3**2 - sp.Rational(94, 96) * z6)
    Cm = pi**7 / 3456 - pi**7 / 2880 + (pi / 16) * (sp.Rational(61, 6) * z6 - 6 * z3**2)          # <K2^2>: Parseval; Euler's sums of H_n/n^5 and H_n^2/n^4
    Dm = 11 * pi**5 / 1440                                                                      # <th^2 lt^2>: Fourier; Euler's sum of H_n/n^3
    lt2D3 = (pi / 4) * (5 * z5 - 2 * z2 * z3) + pi**3 * z3 / 32                                 # <lt^2 D3>: #220's integral of ell^2 D3, with lt = ell + ln 2
    m = {(0, 0, 2, 0, 0): Cm, (1, 0, 1, 0, 0): sp.Rational(7, 8) * z3 * Jt[2] - lt2D3 / 2, (1, 0, 1, 0, 1): -A / 4, (2, 0, 0, 0, 0): pi**5 / 360,
         (2, 0, 0, 0, 1): sp.Integer(0), (2, 0, 0, 0, 2): A, (0, 0, 1, 1, 0): pi**2 / 8 * Jt[2] - Dm / 2, (0, 1, 0, 0, 3): B, (1, 0, 0, 1, 0): pi * z3 / 4,
         (1, 0, 0, 1, 1): -pi**5 / 1440, (0, 0, 0, 2, 0): pi**3 / 24}
    for b in range(7):
        m[(0, 0, 0, 0, b)] = Jt[b]
    return m


RD_NEED_K2 = {(0, 0, 2, 0, 0), (1, 0, 1, 0, 0), (1, 0, 1, 0, 1), (0, 0, 1, 1, 0)}


def rd_quad_master(key, dps):
    """<C^a D3^b K2^c th^j lt^k> by quadrature at dps digits"""
    mp.mp.dps = dps
    q_, z3 = mp.ln2, mp.zeta(3)
    a, b, c, j, k = key

    def fn(th, x, tail):
        if tail:
            th = mp.pi / 2 - x
        with mp.extradps(60):
            lt = mp.log(2 * mp.sin(x)) if tail else mp.log(2 * mp.cos(th))
            a2 = 2 * x if tail else mp.pi - 2 * th
            C = mp.clsin(2, a2)
            D3 = mp.clcos(3, a2) + 3 * z3 / 4 if b else 0
            K2 = (o4.Lc2_poly(th) + q_ * C - q_**2 * th) if c else 0
            r = C**a * D3**b * K2**c * th**j * lt**k
        return +r
    return o4.both(fn)


def rd_kit_check():
    """RD4: the series of the generators at pi/2 and 0 against the functions"""
    pi2, zero = RD_STATE["kit"] if "kit" in RD_STATE else (RD_STATE.setdefault("kit", rd_kit()))
    mp.mp.dps = 60
    q_, z3 = mp.ln2, mp.zeta(3)
    Jt2v, Jt3v = mp.pi**3 / 24, -3 * mp.pi * z3 / 4
    xv = mp.mpf(1) / 1000
    th_v = mp.pi / 2 - xv
    K2f = lambda a: o4.Lc2_poly(a) + q_ * mp.clsin(2, mp.pi - 2 * a) - q_**2 * a  # noqa: E731
    K3f = lambda a: mp.quad(lambda p: mp.log(2 * mp.cos(p)) ** 3, [0, a])  # noqa: E731
    ref = {"lt": mp.log(2 * mp.cos(th_v)), "u": mp.tan(th_v), "ui": 1 / mp.tan(th_v), "C": mp.clsin(2, mp.pi - 2 * th_v),
           "D3": mp.clcos(3, mp.pi - 2 * th_v) + 3 * z3 / 4, "K2": K2f(th_v), "K3": K3f(th_v)}
    subs = {RS_X: xv, RS_L: mp.log(xv), X_Q: q_, RS.Jt2: Jt2v, RS.Jt3: Jt3v}
    e1 = max(abs(mp.mpf(sp.N(pi2[k].subs(subs), 70)) - ref[k]) / max(abs(ref[k]), mp.mpf(1)) for k in ref)
    tv0 = mp.mpf(1) / 1000
    ref0 = {"lt": mp.log(2 * mp.cos(tv0)), "u": mp.tan(tv0), "ui": 1 / mp.tan(tv0), "C": mp.clsin(2, mp.pi - 2 * tv0),
            "D3": mp.clcos(3, mp.pi - 2 * tv0) + 3 * z3 / 4, "K2": K2f(tv0), "K3": K3f(tv0)}
    e2 = max(abs(mp.mpf(sp.N(zero[k].subs({RS_T: tv0, X_Q: q_}), 70)) - ref0[k]) / max(abs(ref0[k]), mp.mpf(1e-30)) for k in ref0)
    return e1, e2


def check_reduce7(ev=None):
    t0 = time.time()
    mp.mp.dps = 45
    parts, tot = r7_parts()
    cl = closed_omegas()
    pi, z3 = sp.pi, sp.zeta(3)
    # RD1: the decomposition against the numeric integrand of closed forms at four points
    w5v = mp.mpf("0.3")
    Om = [cl[0], cl[1], cl[2], cl[3], cl[4], w5v, mp.mpf(0)]
    gens = (X_LG, X_U, X_TH, X_C, X_D3, X_LC, X_Q, X_J2) + tuple(sp.symbols("v1:6"))
    tot_f = sp.lambdify(gens, sum(cf * X_C**a * X_D3**b * X_LC**c for (a, b, c), cf in parts.items()), "mpmath")
    I0 = rd_integrand()
    I0_f = sp.lambdify((RS.lt, RS.u, RS.th, RS.C, RS.D3, RS.K2, RS.v5, X_Q), I0, "mpmath")
    fn = r7_numeric_integrand(Om)
    worst, worst2 = mp.mpf(99), mp.mpf(99)
    for tv in ("0.2", "0.7", "1.1", "1.4"):
        tv = mp.mpf(tv)
        with mp.extradps(60):
            ell, s2, u_, Cn, D3n = o4.basics(tv, mp.pi / 2 - tv, False)
            Lc2v = o4.Lc2_poly(tv)
            dec = tot_f(ell, u_, tv, Cn, D3n, Lc2v, mp.ln2, mp.pi * (12 * mp.ln2**2 + mp.pi**2) / 24, cl[1], cl[2], cl[3], cl[4], w5v)
            dec2 = I0_f(ell + mp.ln2, u_, tv, Cn, D3n, Lc2v + mp.ln2 * Cn - mp.ln2**2 * tv, w5v, mp.ln2)
            num = fn(tv, mp.pi / 2 - tv, False)
        worst = min(worst, digits(num, dec))
        worst2 = min(worst2, digits(num, dec2))
    report("RD1", worst >= 30 and worst2 >= 30, f"the decomposition (sympy; coefficients of C^a D3^b Lc2^c: {sorted(parts)}) against the numeric integrand of closed forms at theta = 0.2, 0.7, 1.1, 1.4: "
           f"worst {float(worst):.1f} digits; the same integrand in lt = ln(2 cos th) and K2 = int lt^2 (the masters' variables): {float(worst2):.1f} digits (bar 30)")
    # the chain of integrations by parts
    t1 = time.time()
    I1, Fs, notes = rd_chain(I0)
    dFs = sum(sp.expand(rd_dth(F).subs(RS.ui, 1 / RS.u)) for _, F in Fs)
    r_id = sp.expand(I0 - dFs - I1)
    report("RD2", r_id == 0 and notes["DD primitive residual"] == 0 and notes["D primitive residual"] == 0 and not notes["u or 1/u left"],
           f"{len(Fs)} integrations by parts I <- I - dF (steps {', '.join(n for n, _ in Fs)}): I - sum dF - I_final = {r_id}; the primitives' residuals {notes['DD primitive residual']}, "
           f"{notes['D primitive residual']}; no u or 1/u left in I_final: {not notes['u or 1/u left']} ({time.time() - t1:.0f} s)")
    Bs = [(n, rd_bnd(F)) for n, F in Fs]
    if MUTATE == "rd_B":
        Bs = [(n, -b if n == "CL" else b) for n, b in Bs]
    Jt = rd_Jt()
    Bsum = sp.expand(sum(b for _, b in Bs).subs({RS.Jt2: Jt[2], RS.Jt3: Jt[3]}))
    masters = rd_masters(Jt)
    Pm = sp.Poly(I1, RS.C, RS.D3, RS.K2, RS.th, RS.lt)
    total = Bsum + sum(cf * masters[m] for m, cf in Pm.terms())
    total = sp.expand(total)
    slope = sp.simplify(sp.diff(total, RS.v5))
    cancel_ok = (notes["CD: 1/u part after DD"] == 0 and notes["C: 1/u part after D"] == 0 and notes["K3 apart from C K3"] == 0)
    report("RD3", cancel_ok and sp.simplify(slope + 4 * pi) == 0,
           f"cancellations: the 1/u part of the C D3 coefficient after the D3^2 step: {notes['CD: 1/u part after DD']}; the 1/u part of the C coefficient after the D3 step: "
           f"{notes['C: 1/u part after D']}; K3 only in C K3: {notes['K3 apart from C K3'] == 0}; slope of G_7 in Omega5 = {slope} (= -4 pi)")
    e1, e2 = rd_kit_check()
    report("RD4", max(e1, e2) <= mp.mpf("1e-24"), f"the asymptotic series (order {RS_ORD}) of lt, u, 1/u, C, D3, K2, K3 against the functions at x = pi/2 - theta = 1e-3 (worst relative "
           f"difference {nd(e1, 3)}) and at theta = 1e-3 ({nd(e2, 3)}) (bar 1e-24)")
    # RD5: the masters, numerically (30 digits)
    dm = {}
    t2 = time.time()
    for key, cf in sorted(masters.items(), key=lambda kv: kv[0]):
        mp.mp.dps = 30
        vnum = rd_quad_master(key, 30)
        mp.mp.dps = 30
        vcl = mp.mpf(sp.N(cf, 40))
        dm[key] = digits(vcl, vnum) if vcl != 0 else (mp.mpf(30) if abs(vnum) < mp.mpf("1e-28") else -mp.log10(abs(vnum)))
    report("RD5", min(dm.values()) >= 25, "the master integrals <C^a D3^b K2^c th^j lt^k> against their closed forms (quadrature at 30 digits), digits: "
           + ", ".join(f"({','.join(str(i) for i in k)}) {float(v):.1f}" for k, v in dm.items()) + f" (bar 25; {time.time() - t2:.0f} s)")
    # RD5b: the two masters identified numerically (found at 55 digits) at 100 digits
    t2b = time.time()
    dmb = {}
    for key in ((2, 0, 0, 0, 2), (0, 1, 0, 0, 3)):
        mp.mp.dps = 100
        vnum = rd_quad_master(key, 100)
        mp.mp.dps = 100
        dmb[key] = digits(mp.mpf(sp.N(masters[key], 110)), vnum)
    report("RD5b", min(dmb.values()) >= 90, "the two identified masters <lt^2 C^2> = pi (3 zeta(3)^2/4 - 95 zeta(6)/96) and <lt^3 D3> = pi (9 zeta(3)^2/16 - 11 zeta(6)/4) (found by PSLQ at 55 digits "
           "in {pi^7, pi zeta(3)^2}) against quadrature at 100 digits, digits: " + ", ".join(f"({','.join(str(i) for i in k)}) {float(v):.1f}" for k, v in dmb.items()) + f" (bar 90; {time.time() - t2b:.0f} s)")
    # RD6: the reduction numerically: sum of the boundary terms + quadrature of I_final at Omega5 = 0 against method A
    mp.mp.dps = 45
    Jtn = [mp.mpf(sp.N(j, 60)) for j in Jt]
    q_ = mp.ln2
    Bnum = mp.mpf(0)
    for n, b in Bs:
        Bnum += mp.mpf(str(sp.N(b.subs({RS.Jt2: sp.Float(str(Jtn[2]), 50), RS.Jt3: sp.Float(str(Jtn[3]), 50), X_Q: sp.log(2)}), 50)))
    I1_f = sp.lambdify((RS.lt, RS.th, RS.C, RS.D3, RS.K2, RS.v5, X_Q), I1, "mpmath")

    def fn1(th, x, tail):
        if tail:
            th = mp.pi / 2 - x
        with mp.extradps(80):
            lt = mp.log(2 * mp.sin(x)) if tail else mp.log(2 * mp.cos(th))
            a2 = 2 * x if tail else mp.pi - 2 * th
            Cv = mp.clsin(2, a2)
            r = I1_f(lt, th, Cv, mp.clcos(3, a2) + 3 * mp.zeta(3) / 4, o4.Lc2_poly(th) + q_ * Cv - q_**2 * th, mp.mpf(0), q_)
        return +r
    g_red = Bnum + o4.both(fn1)
    w_ = [cl[0], cl[1], cl[2], cl[3], cl[4]]

    def G(o5):
        f = r7_numeric_integrand(w_ + [o5, mp.mpf(0)])
        return (mp.quad(lambda t: f(t, mp.pi / 2 - t, False), [0, mp.mpf(1) / 2, 1], method="gauss-legendre")
                + mp.quad(lambda x: f(mp.pi / 2 - x, x, True), [0, mp.pi / 2 - 1]))
    t3 = time.time()
    g0n, g1n = G(mp.mpf(0)), G(mp.mpf(1))
    slope_n, om5_a = g1n - g0n, -g0n / (g1n - g0n)
    mp.mp.dps = 60
    dred = digits(g0n, g_red)
    report("RD6", dred >= 30, f"G_7(Omega5 = 0) from the integrations by parts (sum of the boundary terms {nd(Bnum, 25)} + quadrature of I_final) = {nd(g_red, 32)} against the integrand of closed forms by "
           f"quadrature = {nd(g0n, 32)}: {float(dred):.1f} digits (bar 30)")
    # A7: method A
    ds = digits(-4 * mp.pi, slope_n)
    msg = f"method A at order 7 (closed forms, E1 .. E5, quadrature at 40 digits, {time.time() - t3:.0f} s): slope against -4 pi {float(ds):.1f} digits; Omega_5 = {nd(om5_a, 38)}"
    mp.mp.dps = 600
    om5_closed = sp.N(total.subs(RS.v5, 0).subs(X_Q, sp.log(2)) / (4 * sp.pi), 520)
    om5_c = mp.mpf(str(om5_closed))
    RD_STATE["om5_closed_expr"] = sp.expand(total.subs(RS.v5, 0) / (4 * pi))
    if ev is not None:
        dc = {n: digits(ev[n]["Om5"], om5_a) for n in ev}
        report("A7", ds >= 30 and all(dc[n] >= 30 for n in ("G1", "G2") if n in dc), msg + " against the grids' Omega_5, digits: " + ", ".join(f"{n} {float(d):.1f}" for n, d in dc.items()) + " (bar 30)")
        dd = {n: digits(ev[n]["Om5"], om5_c) for n in ev}
        report("RD7", all(dd[n] >= 100 for n in ("G1", "G2") if n in dd), "Omega5 of the closed form of the reduction (boundary terms + masters, A and B identified) against the grids' Omega5, digits: "
               + ", ".join(f"{n} {float(d):.1f}" for n, d in dd.items()) + " (bar 100 on G1, G2; the others reported)")
    else:
        print("  " + msg)
        print(f"  closed form of the reduction: Omega5 = {mp.nstr(om5_c, 60)}")
    mp.mp.dps = 60
    print("  Omega5 (derived, symbolic, in q = ln 2) = " + str(RD_STATE["om5_closed_expr"]))
    print(f"  ({time.time() - t0:.0f} s)")
    return om5_a


def main():
    parse(sys.argv[1:])
    mp.mp.dps = 60
    print(f"order7.py  mutate={MUTATE}  parts={PARTS}  mpmath {mp.__version__}  sympy {sp.__version__}", flush=True)
    t0 = time.time()
    if "grid" in PARTS:
        check_grid()
    if "cache" in PARTS:
        check_cache()
    if "sym" in PARTS:
        check_sym()
    ev = None
    if "evals" in PARTS:
        ev = check_evals()
        if "values" in PARTS or not GRIDSEL:
            check_values(ev)
    elif LOAD and any(p in PARTS for p in ("values", "pslq", "remainder", "reduce7")):
        mp.mp.dps = 600
        ev = load_ev(LOAD)
        if "values" in PARTS:
            check_values(ev)
    if "xsym" in PARTS:
        check_xsym()
    if "exch7" in PARTS:
        check_exch7()
    if "reduce7" in PARTS:
        check_reduce7(ev)
    if "pslq" in PARTS:
        check_pslq(ev)
    if "remainder" in PARTS:
        check_remainder(ev)
    nfail = sum(1 for _, ok in RESULTS if not ok)
    print(f"\n{len(RESULTS) - nfail} passed, {nfail} failed of {len(RESULTS)} checks; {time.time() - t0:.0f} s")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()

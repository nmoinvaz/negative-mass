"""Order 6 of the matched expansion (matched.py, order3.py ... order5.py): Omega4 from G6 = 0 by method B, and the checks of issue #220.

W^2 a^2 = s + Omega, s = 2(D - 3), eps = 1/s, Omega = Omega0 + Omega1 eps + Omega2 eps^2 + ..., R = r^s, l = ln R.
The order equations (E0 derives them from the master equation), the same for every k:
  dY_k/dl = F_{k-1},  F_{k-1} = sum_{j<k} L_{k-1-j} l^j / j!,   (R - 1) L_k' - L_k/2 = h_k,
  h_k = (1 - 1/R) L_{k-1} + sum_{j<k} Y_j c_{k-j} + Y_k/2,   L_k = -2 Omega_{k-1} + u C_k,   G_k = C_k(inf) = 0,
with u = sqrt(R - 1), C_k = int_1^R (h_k - Omega_{k-1}) (R' - 1)^{-3/2} dR'. G_k is affine in Omega_{k-2} with slope -4 pi:
G2 = 0 gives Omega0, G3 = 0 Omega1, G4 = 0 Omega2, G5 = 0 Omega3, **G6 = 0 gives Omega4**. One path computes them all (`omega_chain`).

Method B (order4.py's PanelB, shared): v = sqrt(l), panels of width w on [0, V], N Chebyshev nodes per panel, cumulative integrals by the
spectral integration matrix; Y_k from the left, the tail C_right from the right, L_k = -2 h_k(1) - u C_right, G_k = C_right(0). It uses no
closed form of Y_k or L_k and no input value of Omega_0 .. Omega_3.

    python3 order6.py grid --grid=N,w,V,dps [--kmax=6]      (a development run of one grid; no gate)
"""
import json
import math
import os
import random
import re
import sys
import time

import mpmath as mp
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_argv, sys.argv = sys.argv, sys.argv[:1]          # order4.py parses sys.argv when it is imported
from order4 import PanelB  # noqa: E402
import order4 as o4  # noqa: E402
sys.argv = _argv

MUTATIONS = ("wrong_h", "wrong_c6", "wrong_Om0", "wrong_Om1", "wrong_Om2", "wrong_Om3", "coarse_B", "short_tail", "pslq_wrong_basis",
             "wrong_Om4_series", "wrong_closed3", "pslq_wrong_target", "e1_sign", "e2_weight", "e3_factor", "e4_g", "e4_noC",
             "rd_w4", "rd_M1", "rd_M3", "rd_master", "rd_J5", "rd_kit", "rd_kit2", "rd_Dstep")
MUTATE = None
W2FILE = "w2_output.txt"
SAVE = None
LOAD = None
GRIDARG = None
GRIDSEL = None
CGRIDS_ON = False
W2RAW = None
KMAX = 6
PARTS = []
REAL = True


def parse(argv):
    global MUTATE, W2FILE, W2RAW, SAVE, LOAD, GRIDARG, GRIDSEL, CGRIDS_ON, KMAX, PARTS, REAL
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
        elif a == "--cgrids":
            CGRIDS_ON = True
        elif a == "--no-real":
            REAL = False
        elif a.startswith("--kmax="):
            KMAX = int(a.split("=", 1)[1])
        else:
            parts.append(a)
    PARTS = parts or ["sym", "evals", "pslq", "remainder"]


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
    raise ValueError(m)


def hk(m, Om, l, R, Y, L):
    """h_m"""
    if m == 0:
        return Y[0] / 2
    return (1 - 1 / R) * L[m - 1] + sum(Y[j] * cc(m - j, Om, l, R) for j in range(m)) + Y[m] / 2


def Fsum(k, L, l):
    return sum(L[k - 1 - j] * l**j / math.factorial(j) for j in range(k))


# ---------------------------------------------------------------------------------------------------------------
# E0: the epsilon-expansion of the master equation against the table, k = 0 .. 6
# ---------------------------------------------------------------------------------------------------------------
def check_sym():
    t0 = time.time()
    R, eps, l = sp.symbols("R epsilon l", positive=True)
    Om = sp.symbols("Omega0:6")
    Ys, Ls = sp.symbols("Y0:7"), sp.symbols("L0:7")
    kk = 6
    r = sum((eps * l) ** j / sp.factorial(j) for j in range(kk + 1))          # exact through eps^6, which is all that is read
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
           f"h_k (k = 0 .. 6) and dY_k/dl (k = 1 .. 6) from the epsilon-expansion of the master equation against the table: "
           f"residuals {[str(x) for x in res if x != 0] or '0'} ({time.time() - t0:.1f} s)")


# ---------------------------------------------------------------------------------------------------------------
# method B at order k
# ---------------------------------------------------------------------------------------------------------------
class PanelK(PanelB):
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

    def hierarchy(self, Om, kmax):
        """Y_0 .. Y_kmax, L_0 .. L_kmax and G_0 .. G_kmax for the given Omega list (Om[kmax - 1] is not used: h_kmax(1) = Om[kmax - 1])"""
        N, P = self.N, self.P
        grid = lambda fn: [[fn(p, i) for i in range(N)] for p in range(P)]  # noqa: E731
        lp = [grid(lambda p, i, j=j: self.l[p][i] ** j / mp.factorial(j)) for j in range(kmax + 1)]       # l^j / j!
        cg = {m: grid(lambda p, i, m=m: cc(m, Om, self.l[p][i], self.R[p][i])) for m in range(1, kmax + 1)}
        Sg = grid(lambda p, i: self.Rm1[p][i] / self.R[p][i])                                              # 1 - 1/R
        Y = [grid(lambda p, i: mp.mpf(1))]
        L = []
        G = {}
        for k in range(kmax + 1):
            if k >= 1:
                F = grid(lambda p, i: 2 * self.v[p][i] * mp.fdot([L[k - 1 - j][p][i] for j in range(k)], [lp[j][p][i] for j in range(k)]))
                Y.append(self.cum_left(F))
            thr = mp.mpf(1) / 2 if k == 0 else Om[k - 1]
            if k == 0:
                s = grid(lambda p, i: Y[0][p][i] / 2 - thr)
            else:
                s = grid(lambda p, i: Sg[p][i] * L[k - 1][p][i] + mp.fdot([Y[j][p][i] for j in range(k)], [cg[k - j][p][i] for j in range(k)])
                         + Y[k][p][i] / 2 - thr)
            Cint = grid(lambda p, i: s[p][i] / self.Rm1[p][i] * self.R[p][i] / self.u[p][i] * 2 * self.v[p][i])
            Cr, G[k] = self.cum_right(Cint)
            L.append(grid(lambda p, i: -2 * thr - self.u[p][i] * Cr[p][i]))
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
    """Omega_0 .. Omega_3 of #26 and #208 (derived there), the planted truth"""
    pi, ln2, z3 = mp.pi, mp.ln2, mp.zeta(3)
    return [-4 * ln2, 2 * pi**2 / 3 + 8 * ln2**2,
            80 * z3 - mp.mpf(8) / 3 * pi**2 * ln2 - mp.mpf(32) / 3 * ln2**3,
            22 * pi**4 / 45 + mp.mpf(16) / 3 * pi**2 * ln2**2 + mp.mpf(32) / 3 * ln2**4 - 320 * z3 * ln2 - 384 * z3
            + (mp.mpf(1) / 1000 if MUTATE == "wrong_closed3" else 0)]


def omega_chain(B, kmax, log=None):
    """Omega_0 .. Omega_{kmax-2} from G_2 .. G_kmax = 0, in turn, by the same code at every order; returns (Om, slopes, seconds per order)"""
    Om = [mp.mpf(0)] * kmax
    slopes, secs = {}, {}
    cl = closed_omegas()
    for k in range(2, kmax + 1):
        t0 = time.time()
        Om_in = list(Om)
        for i, name in enumerate(("wrong_Om0", "wrong_Om1", "wrong_Om2", "wrong_Om3")):
            if MUTATE == name and k >= i + 3:
                Om_in[i] = cl[i] + mp.mpf(1) / 100
        val, slope = solve_omega(lambda kk, o: B.hierarchy(o, kk)[2][kk], k, Om_in)
        Om[k - 2] = val
        slopes[k] = slope
        secs[k] = time.time() - t0
        if log:
            print(f"    k = {k}: Omega_{k - 2} = {nd(val, 30)}  slope {nd(slope, 20)}  ({secs[k]:.0f} s)", flush=True)
    return Om, slopes, secs


def run_B(N, w, V, dps, kmax=6, log=False):
    t0 = time.time()
    if MUTATE == "coarse_B":
        N = 10
    if MUTATE == "short_tail":
        V = 6
    with mp.workdps(dps):
        w = mp.mpf(w[0]) / w[1] if isinstance(w, tuple) else mp.mpf(w)
        B = PanelK(N, w, V)
        Om, slopes, secs = omega_chain(B, kmax, log=log)
        out = {f"Om{i}": +Om[i] for i in range(kmax - 1)}
        out.update({f"slope{k}": +s for k, s in slopes.items()})
        out["secs"] = time.time() - t0
        out["grid"] = f"N={N} w={mp.nstr(w, 6)} V={V} dps={dps} panels={B.P}"
    return out


# declared grids: name -> (N, w, V, working digits). G1, G2 are the grids of record; B1 .. B3 are reported only; C1, C2 are the grids of the
# control runs (--cgrids), sized to meet the same bars (100 digits) with a few digits to spare.
GRIDS = {"G1": (170, (1, 4), 38, 400), "G2": (150, (3, 10), 36, 360)}              # w as a fraction (num, den), evaluated at the working precision
REPORT_GRIDS = {"B1": (40, (1, 2), 18, 120), "B2": (60, (2, 5), 24, 200), "B3": (100, (1, 4), 30, 260)}
CGRIDS = {"G1": (84, (1, 4), 28, 200), "G2": (80, (3, 10), 28, 190)}

OM4_FIT = "-12794.499860023"          # #195's fit (w2_output.txt: 14 digits, range over the 6 admitted cells 1.5e-30)
OM4_FIT_HALF_ULP = mp.mpf("5e-10")


def save_ev(ev, path):
    with open(path, "w") as f:
        json.dump({n: {k: (mp.nstr(v, 430) if isinstance(v, mp.mpf) else v) for k, v in e.items()} for n, e in ev.items()}, f, indent=1)


def load_ev(paths):
    out = {}
    for path in paths.split(","):
        raw = json.load(open(path))
        out.update({n: {k: (mp.mpf(v) if k.startswith(("Om", "slope")) else v) for k, v in e.items()} for n, e in raw.items()})
    return out


def check_grid():
    """one grid, no gate: the development run"""
    N, w, V, dps = GRIDARG.split(",")
    mp.mp.dps = 60
    e = run_B(int(N), tuple(int(t) for t in w.split("/")) if "/" in w else w, int(V), int(dps), kmax=KMAX, log=True)
    mp.mp.dps = int(dps)
    cl = closed_omegas()
    print(f"  {e['grid']}  {e['secs']:.0f} s")
    for i in range(KMAX - 1):
        print(f"    Omega{i} = {nd(e[f'Om{i}'], int(dps))}")
    for i in range(min(4, KMAX - 1)):
        print(f"    Omega{i} against closed form: {float(digits(cl[i], e[f'Om{i}'])):.1f} digits")
    if "Om4" in e:
        print(f"    Omega4 against #195's fit {OM4_FIT}: difference {nd(e['Om4'] - mp.mpf(OM4_FIT), 5)}")
    for k in range(2, KMAX + 1):
        print(f"    slope G{k} against -4 pi: {float(digits(-4 * mp.pi, e[f'slope{k}'])):.1f} digits")
    return e


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
        print(f"  {name}: {e['grid']}  {e['secs']:.0f} s\n    Omega0 = {nd(e['Om0'], 60)}\n    Omega1 = {nd(e['Om1'], 60)}"
              f"\n    Omega2 = {nd(e['Om2'], 60)}\n    Omega3 = {nd(e['Om3'], 110)}\n    Omega4 = {nd(e['Om4'], dps - 20)}\n    slopes G2 .. G6 = "
              + " ".join(nd(e[f"slope{k}"], 20) for k in range(2, 7)), flush=True)
    if SAVE:
        save_ev(ev, SAVE)
    return ev


def check_values(ev):
    mp.mp.dps = 430
    cl = closed_omegas()
    grids = ["G1", "G2"]
    # PL0 .. PL3: the planted orders by the same code
    for i, cid in enumerate(("PL0", "PL1", "PL2", "PL3")):
        rows = {n: digits(cl[i], ev[n][f"Om{i}"]) for n in ev}
        ok = all(rows[n] >= 100 for n in grids)
        report(cid, ok, f"planted: Omega{i} against its closed form, digits: " + ", ".join(f"{n} {float(d):.1f}" for n, d in rows.items())
               + " (bar 100 on G1, G2; the others reported)")
    # PLW: a wrong input must move the orders above it (the control's own statement, evaluated when a wrong-input mutation is on)
    if MUTATE in ("wrong_Om0", "wrong_Om1", "wrong_Om2", "wrong_Om3"):
        i = int(MUTATE[-1]) + 1
        moved = all(abs(ev[n][f"Om{i}"] - (cl[i] if i < 4 else mp.mpf(OM4_FIT))) > mp.mpf(1) / 1000 for n in grids)
        report("PLW", moved, f"control {MUTATE} moves Omega{i} of " + ", ".join(
            f"{n} by {nd(ev[n][f'Om{i}'] - (cl[i] if i < 4 else mp.mpf(OM4_FIT)), 5)}" for n in grids))
    # SL: slopes
    rows, ok = [], True
    for n in ev:
        ds = [digits(-4 * mp.pi, ev[n][f"slope{k}"]) for k in range(2, 7)]
        rows.append(f"{n} " + "/".join(f"{float(d):.1f}" for d in ds))
        ok = ok and (n not in grids or min(ds) >= 100)
    report("SL", ok, "slope of G_k in Omega_(k-2) against -4 pi, k = 2 .. 6, digits: " + ", ".join(rows) + " (bar 100 on G1, G2)")
    # M4: the two grids on Omega4
    a, b = grids
    d4 = digits(ev[a]["Om4"], ev[b]["Om4"])
    report("M4", d4 >= 100, f"Omega4 of {a} against {b} (two grids, two precisions): {float(d4):.1f} digits (bar 100); |diff| = "
           f"{nd(abs(ev[a]['Om4'] - ev[b]['Om4']), 3)}")
    others = [n for n in ev if n not in grids]
    print("  reported: digits between evaluations, Omega3 / Omega4: " + ", ".join(
        f"{x}-{y} " + "/".join(f"{float(digits(ev[x][f'Om{i}'], ev[y][f'Om{i}'])):.1f}" for i in (3, 4))
        for x, y in [("G1", "G2")] + [(o, "G1") for o in others]))
    # X1: Omega4 against #195's fit
    fit = mp.mpf(OM4_FIT)
    dd = {n: abs(ev[n]["Om4"] - fit) for n in ev}
    report("X1", all(dd[n] <= OM4_FIT_HALF_ULP for n in grids),
           "Omega4 against #195's fit " + OM4_FIT + " (14 digits, half a unit of its last digit is 5e-10; the fit's range over its cells is 1.5e-30): "
           + ", ".join(f"{n} {nd(d, 3)}" for n, d in dd.items()) + " (bar 5e-10 on G1, G2)")
    return d4


# ---------------------------------------------------------------------------------------------------------------
# PSLQ
# ---------------------------------------------------------------------------------------------------------------
PSLQ_DIGITS = 200             # cap on the working digits; the digits the two grids share if fewer
PSLQ_SLACK = 12               # tolerance 10^-(d - slack)
# basis name -> (maxcoeff H, planted seed)
PSLQ_BASES = {"P1": (10**9, 220001), "P2": (10**8, 220002), "P3": (10**7, 220003), "P4": (10**9, 220004)}
PSLQ_BASES_C = {"P1": (10**9, 220001), "P2": (10**5, 220002), "P3": (10**4, 220003), "P4": (10**9, 220004)}   # --cgrids: d ~ 113 < the budgets above


def constants(name):
    pi, ln2, z3, z5 = mp.pi, mp.ln2, mp.zeta(3), mp.zeta(5)
    li4, li5 = mp.polylog(4, mp.mpf(1) / 2), mp.polylog(5, mp.mpf(1) / 2)
    w5 = [("pi^4 ln 2", pi**4 * ln2), ("pi^2 ln^3 2", pi**2 * ln2**3), ("ln^5 2", ln2**5), ("pi^2 zeta(3)", pi**2 * z3),
          ("zeta(3) ln^2 2", z3 * ln2**2), ("zeta(5)", z5)]
    w4 = [("pi^4", pi**4), ("pi^2 ln^2 2", pi**2 * ln2**2), ("ln^4 2", ln2**4), ("zeta(3) ln 2", z3 * ln2)]
    low = [("zeta(3)", z3), ("pi^2 ln 2", pi**2 * ln2), ("ln^3 2", ln2**3), ("pi^2", pi**2), ("ln^2 2", ln2**2), ("ln 2", ln2), ("1", mp.mpf(1))]
    if name == "P1":
        return w5
    if name == "P2":
        return w5 + w4 + low
    if name == "P3":
        return w5 + w4 + low + [("Li4(1/2)", li4), ("Li4(1/2) ln 2", li4 * ln2), ("Li5(1/2)", li5)]
    if name == "P4":
        return w5 + [("Li4(1/2) ln 2", li4 * ln2), ("Li5(1/2)", li5)]
    raise ValueError(name)


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
    mp.mp.dps = 430
    agree = digits(ev[grids[0]]["Om4"], ev[grids[1]]["Om4"])
    d = int(max(30, min(PSLQ_DIGITS, math.floor(agree))))      # PSLQ needs at least 53 bits
    pert = mp.sqrt(2) / 1000        # irrational, so that the shifted constant is outside every basis's rational span (1 is in P2, P3)
    print(f"  PSLQ at {d} digits (the two grids share {float(agree):.1f}), tolerance 1e-{d - PSLQ_SLACK}")
    for name, (H, seed) in (PSLQ_BASES_C if CGRIDS_ON else PSLQ_BASES).items():
        with mp.workdps(d):
            tau = mp.mpf(10) ** (-(d - PSLQ_SLACK))
            names, cs = zip(*constants(name))
            cs_true = list(cs)
            cs = list(cs)
            if MUTATE == "pslq_wrong_basis":                  # PSLQ is given a basis that is not the one the target was built from
                cs[0] = cs[0] + pert
            n = len(cs)
            budget = (n + 1) * math.log10(2 * H + 1) + 1
            om4 = +ev[grids[1]]["Om4"]
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
            csp[1] = csp[1] + pert
            rpp, _ = pslq_run([xp] + csp, tau, H)
            report("Q" + name[1], okp and rpp is None,
                   f"{name} = {{Omega4, {', '.join(names)}}}, n = {n + 1}, H = {H}, d = {d}: planted relation (den {den}, seed {seed}) "
                   f"recovered: {okp} ({tp:.0f} s); perturbed basis on the planted target finds {'nothing' if rpp is None else list(rpp)}; "
                   f"budget (n + 1) log10(2H + 1) + 1 = {budget:.1f}")
            # the real target
            if not REAL:
                print(f"    REAL {name}: skipped (--no-real)", flush=True)
                continue
            rr, tr = pslq_run([om4] + cs, tau, H)
            found, why = False, "not found"
            if rr is not None:
                res = residual(rr, [om4] + cs)
                found = res <= tau and budget < d and rr[0] != 0      # a relation without Omega4 is among the constants alone
                why = f"relation {list(rr)}, residual {nd(res, 3)} (bar {nd(tau, 2)}), budget with H = {H}: {budget:.1f} < {d}: {budget < d}; Omega4 coefficient {rr[0]}"
                if found:
                    j = next(i for i in range(1, len(rr)) if rr[i] != 0)
                    csr = list(cs)
                    csr[j - 1] = csr[j - 1] + pert
                    rrp, _ = pslq_run([om4] + csr, tau, H)
                    found = found and rrp is None
                    why += f"; perturbed basis (constant {names[j - 1]} + sqrt(2)/1000) finds {'nothing' if rrp is None else list(rrp)}"
            print(f"    REAL {name}: " + (("FOUND (counts) " if found else "found but does not count ") + why if rr is not None
                                           else f"not found at H = {H} ({tr:.0f} s)"), flush=True)


# ---------------------------------------------------------------------------------------------------------------
# E1, E2, E3 at k = 6 and the Y_4 step (E4), on the hierarchy (no closed form, no matching beyond method B)
# ---------------------------------------------------------------------------------------------------------------
EGRID = (36, (1, 2), 16, 80)


class Fields:
    """Y_0..Y_k, L_0..L_k, G_k of the hierarchy for one Omega list, with the nodes' l, R, ell, sin^2, tan and the th-quadrature"""

    def __init__(self, B, Om, k):
        self.B, self.Om, self.k = B, list(Om), k
        self.N, self.P = B.N, B.P
        self.Y, self.L, self.G = B.hierarchy(self.Om, k)
        self.l, self.R, self.u = B.l, B.R, B.u
        self.S = self.grid(lambda p, i: B.Rm1[p][i] / B.R[p][i])                 # sin^2 th = 1 - 1/R
        self.ell = self.grid(lambda p, i: -B.l[p][i] / 2)                        # ln cos th = -l/2
        self.W = self.grid(lambda p, i: B.half * B.wend[i] * B.v[p][i] / B.u[p][i])   # th-quadrature weights: dth = (v/u) dv
        self.C = self.grid(lambda p, i: mp.clsin(2, 2 * mp.asin(mp.exp(-B.l[p][i] / 2))))   # Cl_2(pi - 2 th) = Cl_2(2 x), x = pi/2 - th

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

    def Sk(self, k, p, i, upto=None):
        """sum_{j < upto} Y_j c_{k-j}"""
        upto = k if upto is None else upto
        return sum(self.Y[j][p][i] * cc(k - j, self.Om, self.l[p][i], self.R[p][i]) for j in range(upto))


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


def e4_sides(fl):
    """E4 (the Y_4 step at k = 6): int Y_4 w4 = -2 int (g + 16 u C) F_3, w4 = (2 Omega1 + 8 Omega0 ell - 22 ell^2)/S - 12 ell,
    g = 22 ell^2 - 8 Omega0 ell - 2 Omega1; it needs Omega0 = -4 ln 2"""
    O0, O1 = fl.Om[0], fl.Om[1]
    gg = 20 if MUTATE == "e4_g" else 22
    cfac = 0 if MUTATE == "e4_noC" else 16
    lhs = fl.integral(lambda p, i: fl.Y[4][p][i] * ((2 * O1 + 8 * O0 * fl.ell[p][i] - 22 * fl.ell[p][i] ** 2) / fl.S[p][i] - 12 * fl.ell[p][i]))
    rhs = -2 * fl.integral(lambda p, i: (gg * fl.ell[p][i] ** 2 - 8 * O0 * fl.ell[p][i] - 2 * O1 + cfac * fl.u[p][i] * fl.C[p][i]) * fl.F(3, p, i))
    return lhs, rhs


def g6_reduced(fl):
    """G_6 with no Y_4, Y_5, L_4, L_5: E1, E2 (m = 5 and 4), E3 (m = 5) and E4, on the hierarchy's own Y_0..Y_3, L_0..L_3"""
    om = fl.Om
    O0 = om[0]
    ell, l, S = fl.ell, fl.l, fl.S
    shift = mp.ln2 if MUTATE == "e2_weight" else 0                       # E2's weight ell (used in terms 1 and 5)
    f3 = 2 if MUTATE == "e3_factor" else 4                               # E3's factor
    gg = 20 if MUTATE == "e4_g" else 22
    cfac = 0 if MUTATE == "e4_noC" else 16
    t = []
    t.append(8 * fl.integral(lambda p, i: (ell[p][i] + shift) * fl.Sk(5, p, i, upto=4) / S[p][i]))      # 8 int ell (sum_{j<4} Y_j c_{5-j})/S
    t.append((-2 if MUTATE == "e1_sign" else 2) * fl.integral(lambda p, i: sum(fl.L[5 - j][p][i] * l[p][i] ** j / mp.factorial(j) for j in range(2, 6))))
    t.append(2 * fl.integral(lambda p, i: (fl.Sk(6, p, i, upto=4) - om[5]) / S[p][i]))
    t.append(-f3 * fl.integral(lambda p, i: (ell[p][i] - O0) * sum(fl.L[4 - j][p][i] * l[p][i] ** j / mp.factorial(j) for j in range(1, 5))))
    t.append(8 * O0 * fl.integral(lambda p, i: (ell[p][i] + shift) * (S[p][i] * fl.L[3][p][i] + fl.Sk(4, p, i, upto=4)) / S[p][i]))
    g = lambda p, i: gg * ell[p][i] ** 2 - 8 * O0 * ell[p][i] - 2 * om[1]  # noqa: E731
    t.append(-2 * fl.integral(lambda p, i: (g(p, i) + cfac * fl.u[p][i] * fl.C[p][i]) * fl.F(3, p, i)))
    return sum(t), t


def check_exch6():
    t0 = time.time()
    N, w, V, dps = EGRID
    mp.mp.dps = 60
    with mp.workdps(dps):
        B = PanelK(N, mp.mpf(w[0]) / w[1], V)
        Om, _, _ = omega_chain(B, 6)
        print(f"  E-grid N = {N}, w = {w[0]}/{w[1]}, V = {V}, dps {dps}, {B.P} panels; its Omega0 .. Omega3 meet the closed forms to "
              + ", ".join(f"{float(digits(c, o)):.1f}" for c, o in zip(closed_omegas(), Om)) + " digits")
        res = {}
        for label, o4v in (("root", Om[4]), ("0", mp.mpf(0)), ("1/2", mp.mpf(1) / 2)):
            Omv = list(Om)
            Omv[4] = o4v
            fl = Fields(B, Omv, 6)
            g6, rhs, parts = e1_sides(fl, 6)
            scale = max(abs(x) for x in parts)
            e1 = abs(g6 - rhs) / scale
            lhs3, rhs3 = e3_sides(fl, 5)
            lhs2, rhs2 = e2_sides(fl, 5)
            lhs2b, rhs2b = e2_sides(fl, 4)
            lhs4, rhs4 = e4_sides(fl)
            gr, grp = g6_reduced(fl)
            sc = max(abs(x) for x in grp)
            res[label] = dict(
                E1=(-mp.log10(e1) if e1 else mp.mpf(dps), g6, scale),
                E2_5=(digits(lhs2, rhs2), lhs2, rhs2), E2_4=(digits(lhs2b, rhs2b), lhs2b, rhs2b),
                E3_6=(digits(lhs3, rhs3), lhs3, rhs3), E4=(digits(lhs4, rhs4), lhs4, rhs4),
                G6R=(-mp.log10(abs(g6 - gr) / sc) if g6 != gr else mp.mpf(dps), g6, gr))
        for key, cid, what in (("E1", "E1_6", "G_6 against E1's right side (relative to the largest of its three terms)"),
                               ("E2_5", "E2_5", "int L_5 = 2 int h_5 ell/sin^2"), ("E2_4", "E2_4", "int L_4 = 2 int h_4 ell/sin^2"),
                               ("E3_6", "E3_6", "int Y_5 [(ell - Omega0)/sin^2 + 1] = 2 int (ell - Omega0) F_4"),
                               ("E4", "E4_6", "int Y_4 w4 = -2 int (g + 16 u C) F_3"),
                               ("G6R", "G6R", "G_6 against the reduced form (E1, E2, E3, E4: Y_0 .. Y_3, L_0 .. L_3 only; relative to its largest term)")):
            ds = {lab: r[key][0] for lab, r in res.items()}
            report(cid, all(d >= 30 for d in ds.values()), f"{what}: " + ", ".join(f"{float(d):.1f} digits at Omega4 = {lab}" for lab, d in ds.items())
                   + " (bar 30)")
        print(f"    G_6 at Omega4 = 0 and 1/2: {nd(res['0']['E1'][1], 20)}, {nd(res['1/2']['E1'][1], 20)} ({time.time() - t0:.0f} s)")


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
    """Omega4 from the 64 values of W^2 with Omega0 .. Omega3 fixed: y = (W^2 - [s + Omega0 + Omega1/s + Omega2/s^2 + Omega3/s^3]) s^4 =
    Omega4 + Omega5/s + ..., polynomial extrapolation to 1/s = 0 through the K points of largest s; returns {K: estimate}"""
    pts = []
    for D in sorted(lst):
        s = mp.mpf(2 * (D - 3))
        pts.append((1 / s, (mp.mpf(lst[D]) - (s + om[0] + om[1] / s + om[2] / s**2 + om[3] / s**3)) * s**4))
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


def check_remainder(ev):
    mp.mp.dps = 80
    grids = ["G1", "G2"]
    om = [ev[grids[1]][f"Om{i}"] for i in range(1, 5)]
    if MUTATE == "wrong_Om4_series":
        om[3] = om[3] + 1
    lst = read_w2(W2FILE)
    if W2RAW:
        lst.update(read_w2raw(W2RAW))      # D = 1000, 2000 are not in the 64-value list

    def x(D):
        s = mp.mpf(2 * (D - 3))
        rem = mp.mpf(lst[D]) - (s - 4 * mp.ln2 + om[0] / s + om[1] / s**2 + om[2] / s**3 + om[3] / s**4)
        return rem, rem * s**5
    print("  six-term series s - 4 ln 2 + Omega1/s + ... + Omega4/s^4 against W^2 (qd, w2_output.txt), s = 2(D - 3):")
    print("  D       s       remainder              rem*s^5 (tends to Omega5)")
    xs = {}
    for D in (1000, 2000, 5103, 10000):
        rem, xs[D] = x(D)
        print(f"  {D:<6d} {2 * (D - 3):<7d} {nd(rem, 12):<22s} {nd(xs[D], 12)}")
    falls = xs[1000] > xs[2000] > xs[5103] > xs[10000] > 0
    s1, s2 = mp.mpf(2 * (5103 - 3)), mp.mpf(2 * (10000 - 3))
    est = (s2 * xs[10000] - s1 * xs[5103]) / (s2 - s1)
    report("R5", falls and abs(xs[10000] - xs[5103]) <= 500, f"rem*s^5 at D = 1000, 2000, 5103, 10^4 falls ({falls}); |x(10^4) - x(5103)| = "
           f"{nd(abs(xs[10000] - xs[5103]), 5)} (bar 500: a wrong Omega4 adds (error)*s, 10^4 per unit)")
    report("R6", mp.mpf("4.0e4") <= est <= mp.mpf("4.2e4"), f"two-point extrapolation (s2 x2 - s1 x1)/(s2 - s1) from D = 5103, 10^4: Omega5 ~ {nd(est, 8)} "
           f"(bar 4.0e4 to 4.2e4: an estimate, not a bound)")
    # X2: Omega4 from the crate's 64 values with Omega0 .. Omega3 fixed (the closed forms), not the printed 14 digits of the fit
    mp.mp.dps = 120
    est = refit_windows(read_w2(W2FILE), closed_omegas())
    rows = {n: {K: abs(e - ev[n]["Om4"]) for K, e in est.items()} for n in grids}
    mx = max(max(r.values()) for r in rows.values())
    report("X2", mx <= mp.mpf("1e-25"), "Omega4 against the polynomial extrapolation in 1/s of #195's 64 W^2 values (Omega0 .. Omega3 fixed at the closed forms), "
           "windows of the K largest s, |difference| " + "; ".join(f"{n}: " + ", ".join(f"K = {K} {nd(d, 2)}" for K, d in r.items()) for n, r in rows.items())
           + f" (bar 1e-25 on every window; the spread over windows is an estimator range, {nd(mx, 2)}; not a bound)")
    mp.mp.dps = 80
    print("  reported, not gated: rem*s^5 at the D of the list")
    print("  D       rem*s^5")
    for D in sorted(lst):
        print(f"  {D:<7d} {nd(x(D)[1], 14)}")


# ---------------------------------------------------------------------------------------------------------------
# the reduction of G_6 by hand (RD0 .. RD7) and method A at order 6 (A6)
# ---------------------------------------------------------------------------------------------------------------
# Generators on theta in (0, pi/2): ell = ln cos th, u = tan th, th, C = Cl_2(pi - 2 th), D3 = Cl_3(pi - 2 th) + 3 zeta(3)/4, Lc = int_0^th ell^2.
# Their derivatives close: ell' = -u, u' = 1 + u^2, th' = 1, C' = 2 ell + 2 q, D3' = 2 C, Lc' = ell^2 (q = ln 2).
# Integrals are finite parts at theta -> pi/2 (cutoff x = pi/2 - th = eps: the constant term of the expansion in eps and ln eps) and at theta -> 0;
# every integration by parts holds for them, with the boundary terms the finite parts of F(pi/2) - F(0); the integral of the integrable total is the sum.
R6_LG, R6_U, R6_UI, R6_TH, R6_C, R6_D3, R6_LC = sp.symbols("ell u ui th C D3 Lc2")
R6_Q = sp.Symbol("q", positive=True)
R6_W = sp.symbols("w0:6")
R6_X, R6_L, R6_T, R6_S, R6_R = sp.symbols("x L t s r", positive=True)
R6_Z3, R6_J2 = sp.Symbol("z3"), sp.Symbol("J2")
R6_ORD, R6_SHIFT = 9, 8


def r6_dth(e):
    """d/dtheta on the generators"""
    L_, U_, T_, C_, D_, K_ = R6_LG, R6_U, R6_TH, R6_C, R6_D3, R6_LC
    return (sp.diff(e, L_) * (-U_) + sp.diff(e, U_) * (1 + U_**2) + sp.diff(e, T_) + sp.diff(e, C_) * (2 * L_ + 2 * R6_Q) + sp.diff(e, D_) * 2 * C_
            + sp.diff(e, K_) * L_**2 + sp.diff(e, R6_UI) * (-(1 + R6_UI**2)))


def r6_kit():
    """asymptotic series of the generators at theta -> pi/2 (x = pi/2 - th, L = ln x) and theta -> 0, for the finite parts"""
    x, L, t, s, r, q, ORD = R6_X, R6_L, R6_T, R6_S, R6_R, R6_Q, R6_ORD
    logs = lambda e: sp.expand(e).subs(sp.log(x), L)  # noqa: E731
    ell1 = sp.series(sp.log(sp.sin(s) / s), s, 0, ORD).removeO()
    sgn = -1 if MUTATE == "rd_kit" else 1
    ell_s = sp.log(s) + ell1
    C_s = sp.expand(-sp.integrate(sp.expand(2 * (sp.log(s) + ell1) + 2 * sgn * q), (s, 0, s)))
    D3_s = sp.expand(sp.Rational(7, 4) * R6_Z3 - 2 * sp.integrate(C_s.subs(s, r), (r, 0, s)))
    Lc_s = sp.expand((sp.Rational(101, 100) if MUTATE == "rd_kit2" else 1) * R6_J2 - sp.integrate(sp.expand(ell_s.subs(s, r) ** 2), (r, 0, s)))
    pi2 = {"ell": L + ell1.subs(s, x), "u": sp.series(sp.cot(x), x, 0, ORD).removeO(), "ui": sp.series(sp.tan(x), x, 0, ORD).removeO(),
           "th": sp.pi / 2 - x, "C": logs(C_s.subs(s, x)), "D3": logs(D3_s.subs(s, x)), "Lc": logs(Lc_s.subs(s, x))}
    ell_t = sp.series(sp.log(sp.cos(t)), t, 0, ORD).removeO()
    C_t = sp.expand(sp.integrate(sp.expand(2 * ell_t + 2 * q), (t, 0, t)))
    D3_t = sp.expand(sp.integrate(sp.expand(2 * C_t.subs(t, r)), (r, 0, t)))
    Lc_t = sp.expand(sp.integrate(sp.expand(ell_t.subs(t, r) ** 2), (r, 0, t)))
    zero = {"ell": ell_t, "u": sp.series(sp.tan(t), t, 0, ORD).removeO(), "ui": sp.series(sp.cot(t), t, 0, ORD).removeO(), "th": t, "C": C_t,
            "D3": D3_t, "Lc": Lc_t}
    return pi2, zero


_KIT = {}


def r6_const_term(e, v, logsym=None):
    e = sp.expand(e * v**R6_SHIFT)
    gens = (v, logsym) if logsym is not None else (v,)
    P = sp.Poly(e, *gens)
    return P.coeff_monomial((R6_SHIFT,) + ((0,) if logsym is not None else ()))


def r6_bnd(expr):
    """B[F] = FP_{pi/2} F - FP_0 F"""
    if not _KIT:
        _KIT["pi2"], _KIT["zero"] = r6_kit()
    pi2, zero = _KIT["pi2"], _KIT["zero"]
    expr = sp.expand(expr).replace(lambda e: e.is_Pow and e.base == R6_U and e.exp < 0, lambda e: R6_UI ** (-e.exp))
    names = {"ell": R6_LG, "u": R6_U, "ui": R6_UI, "th": R6_TH, "C": R6_C, "D3": R6_D3, "Lc": R6_LC}
    syms = {k: v for k, v in names.items() if expr.has(v)}
    a = r6_const_term(expr.subs({syms[k]: pi2[k] for k in syms}), R6_X, R6_L)
    b = r6_const_term(expr.subs({syms[k]: zero[k] for k in syms}), R6_T)
    return sp.expand(a - b)


def r6_moments():
    """J_n = int_0^{pi/2} ln^n cos, n <= 6, from the Beta function (as a sympy expression in q = ln 2, pi, zeta)"""
    a = sp.Symbol("a")
    ser = sp.series(sp.sqrt(sp.pi) / 2 * sp.gamma((a + 1) / 2) / sp.gamma(a / 2 + 1), a, 0, 7).removeO()
    out = []
    for n in range(7):
        j = sp.expand_log(sp.simplify(sp.expand_func(sp.factorial(n) * ser.coeff(a, n))), force=True)
        out.append(sp.expand(j.subs(sp.log(2), R6_Q)))
    if MUTATE == "rd_J5":
        out[5] = out[5].subs(sp.zeta(5), sp.Rational(1000, 1080) * sp.zeta(5))
    return out


def r6_parts():
    """the integrand of G_6 after E1, E2, E3, E4 (Y_0 .. Y_3, L_0 .. L_3 only), at Omega_0 = -4 ln 2, as coefficients of the generators C, D3, Lc2
    (u-Laurent polynomials in ell, th, q, w1 .. w4): parts[(a, b, c)] is the coefficient of C^a D3^b Lc2^c"""
    lg, u, th, q = R6_LG, R6_U, R6_TH, R6_Q
    C, D3, Lc, S = R6_C, R6_D3, R6_LC, sp.Symbol("S")
    w = R6_W
    R = 1 / (1 - S)
    l = -2 * lg
    Om = list(w)
    Y = [sp.Integer(1), 2 * lg, 4 * w[0] * lg - 6 * lg**2, 4 * w[1] * lg - 4 * lg**3 + 32 * lg**2 - 16 * w[0] * lg + 32 * u * C - 16 * D3]
    L = [sp.Integer(-1), -2 * w[0] + 4 * lg, -2 * w[1] - 4 * w[0] * lg + 16 * lg**2 + 16 * u * C]
    h3 = hk(3, Om, l, R, Y, L)
    sc = u * (1 - S)
    B = -32 * (lg + q) * C + 16 * S * C + 32 * Lc + 16 * sc * lg * (lg + 2 * q) - 4 * w[1] * th - 2 * w[1] * sc
    L = L + [-2 * h3 + 2 * u * B]
    cx = lambda m, o=Om: cc(m, o, l, R)  # noqa: E731
    o5z = list(Om)
    o5z[5] = 0
    fa = sp.factorial
    T1 = 8 * lg * sum(Y[j] * cx(5 - j) for j in range(4)) / S
    T2 = 2 * sum(L[5 - j] * l**j / fa(j) for j in range(2, 6))
    T3 = 2 * sum(Y[j] * cx(6 - j, o5z) for j in range(4)) / S
    T4 = -4 * (lg - w[0]) * sum(L[4 - j] * l**j / fa(j) for j in range(1, 5))
    T5 = 8 * w[0] * lg * (S * L[3] + sum(Y[j] * cx(4 - j) for j in range(4))) / S
    g = (20 if MUTATE == "rd_w4" else 22) * lg**2 - 8 * w[0] * lg - 2 * w[1]
    F3 = sum(L[3 - j] * l**j / fa(j) for j in range(4))
    T6 = -2 * (g + 16 * u * C) * F3
    tot = (T1 + T2 + T3 + T4 + T5 + T6).subs(w[0], -4 * q)
    tot = sp.expand(sp.cancel(sp.expand(tot).subs(S, u**2 / (1 + u**2))))
    P = sp.Poly(tot, C, D3, Lc)
    return {mon: sp.expand(sp.cancel(cf)) for mon, cf in P.terms()}


R6_RES = {}


def r6_reduce(parts, Jv):
    """the steps Q (three quadratic terms), D (D3 coefficient), C, L and the rest; returns the sum of the finite parts (G_6 at Omega_4 = w4)"""
    lg, u, th, q = R6_LG, R6_U, R6_TH, R6_Q
    C, D3, Lc = R6_C, R6_D3, R6_LC
    z3, z5 = sp.zeta(3), sp.zeta(5)
    pi = sp.pi
    fix = lambda e: sp.expand(e).subs({R6_J2: Jv[2], R6_Z3: z3})  # noqa: E731
    cCC, cCD, cCL = parts[(2, 0, 0)], parts[(1, 1, 0)], parts[(1, 0, 1)]
    cC, cD, cL, cR = parts[(1, 0, 0)], parts[(0, 1, 0)], parts[(0, 0, 1)], parts[(0, 0, 0)]
    z2, z4 = pi**2 / 6, pi**4 / 90
    M_C2, M_lC2, M_thC = pi**5 / 360, -q * pi**5 / 360, pi / 4 * z3
    M_uC = [sp.expand((2 * Jv[b + 2] + 2 * q * Jv[b + 1]) / (b + 1)) for b in range(4)]
    M_l2D3 = sp.Rational(3, 4) * z3 * Jv[2] + pi / 4 * (2 * q * z4 + (4 if MUTATE == "rd_master" else 5) * z5 - 2 * z2 * z3)
    M_CLc = sp.Rational(7, 8) * z3 * Jv[2] - M_l2D3 / 2
    res = R6_RES
    const = sp.Integer(0)
    Ct, Dt, Lt, Rt = sp.expand(cC), sp.expand(cD), sp.expand(cL), sp.expand(cR)
    # Q: M1  <u^2 (f - f_ell) C^2> = B[u f C^2] - <f C^2> - 4 <u f (ell + q) C>
    f = 3072 * lg + 2048 * q + (4095 if MUTATE == "rd_M1" else 4096)
    res["M1 coefficient"] = sp.expand(cCC - u**2 * (f - sp.diff(f, lg)))
    res["M1 derivative"] = sp.expand(r6_dth(u * f * C**2) - (u**2 * (f - sp.diff(f, lg)) * C**2 + f * C**2 + 4 * u * f * (lg + q) * C))
    const += fix(r6_bnd(u * f * C**2)) - (3072 * M_lC2 + (2048 * q + 4096) * M_C2)
    Ct += -4 * u * f * (lg + q)
    # M2  <u C D3> = B[-ell C D3] + <ell (2 ell + 2 q) D3> + 2 <ell C^2>
    res["M2 coefficient"] = sp.expand(cCD + 512 * u)
    res["M2 derivative"] = sp.expand(r6_dth(-lg * C * D3) - (u * C * D3 - lg * (2 * lg + 2 * q) * D3 - 2 * lg * C**2))
    const += -512 * fix(r6_bnd(-lg * C * D3)) - 1024 * M_lC2
    Dt += -512 * lg * (2 * lg + 2 * q)
    # M3  <u^2 C Lc> = B[u C Lc] - <C Lc> - 2 <u (ell + q) Lc> - <u ell^2 C>
    res["M3 coefficient"] = sp.expand(cCL + 2048 * u**2)
    res["M3 derivative"] = sp.expand(r6_dth(u * C * Lc) - (u**2 * C * Lc + C * Lc + 2 * u * (lg + q) * Lc + u * lg**2 * C))
    const += (2048 if MUTATE == "rd_M3" else -2048) * fix(r6_bnd(u * C * Lc)) + 2048 * M_CLc
    Ct += 2048 * u * lg**2
    Lt += 4096 * u * (lg + q)
    # D: <c D3> = B[A D3] - 2 <A C>, A' = c (c0 + c1/sin^2; F_b = int ell^b/sin^2 = -ell^b cot - b P_{b-1}, P_0 = th, P_1 = (C - 2 q th)/2, P_2 = Lc)
    pu = sp.collect(sp.expand(Dt), u, evaluate=False)
    assert set(pu) <= {u**-2, sp.Integer(1)}, list(pu)
    Eq, P0 = pu.get(u**-2, sp.Integer(0)), pu.get(sp.Integer(1), sp.Integer(0))
    c1, c0 = Eq, sp.expand(P0 - Eq)
    a_b = [sp.Poly(c0, lg).coeff_monomial(lg**b) for b in range(3)]
    beta = [sp.Poly(c1, lg).coeff_monomial(lg**b) for b in range(4)]
    assert sp.Poly(c0, lg).degree() <= 2 and sp.Poly(c1, lg).degree() <= 3
    Pf = [th, (C - 2 * q * th) / 2, Lc]
    sg = 1 if MUTATE == "rd_Dstep" else -1
    Fb = lambda b: -lg**b / u + sg * (b * Pf[b - 1] if b else 0)  # noqa: E731
    A = sp.expand(sum(a_b[b] * Pf[b] for b in range(3)) + sum(beta[b] * Fb(b) for b in range(4)))
    res["D derivative"] = sp.expand(r6_dth(A) - (c0 + c1 * (1 + u**-2)))
    const += fix(r6_bnd(A * D3))
    AC = sp.expand(-2 * A * C)
    lin = sum(t_ for t_ in sp.Add.make_args(AC) if t_.has(C) and not t_.has(C**2) and not t_.has(Lc))
    Ct += sp.expand(lin / C)
    rest_AC = sp.expand(AC - lin)
    cC2 = rest_AC.coeff(C, 2)
    cCLc = sp.expand(rest_AC - cC2 * C**2).coeff(Lc, 1).coeff(C, 1)
    assert sp.expand(rest_AC - cC2 * C**2 - cCLc * C * Lc) == 0
    const += cC2 * M_C2 + cCLc * M_CLc
    # C
    pc = sp.collect(sp.expand(Ct), u, evaluate=False)
    assert set(pc) <= {u**-1, u, u**2, sp.Integer(1)}, list(pc)
    res["C: u^-1 part"] = sp.expand(pc.get(u**-1, sp.Integer(0)))
    c1u, c2u, c0u = sp.expand(pc.get(u, 0)), sp.expand(pc.get(u**2, 0)), sp.expand(pc.get(sp.Integer(1), 0))
    assert sp.Poly(c1u, lg).degree() <= 3
    const += sum(sp.Poly(c1u, lg).coeff_monomial(lg**b) * M_uC[b] for b in range(4))
    c2th = sp.Poly(c2u, th).coeff_monomial(th)
    assert sp.expand(c2u - c2th * th) == 0
    res["C: u th derivative"] = sp.expand(r6_dth(u * th * C) - (u**2 * th * C + th * C + u * C + 2 * u * th * (lg + q)))
    res["C: u th (ell + q) derivative"] = sp.expand(r6_dth(-th * (lg + q) ** 2 / 2) - (u * th * (lg + q) - (lg + q) ** 2 / 2))
    uth = fix(r6_bnd(-th * (lg + q) ** 2 / 2)) + sp.Rational(1, 2) * (Jv[2] + 2 * q * Jv[1] + q**2 * Jv[0])
    const += c2th * (fix(r6_bnd(u * th * C)) - M_thC - M_uC[0] - 2 * uth)
    c0t = sp.Poly(c0u, th)
    assert c0t.degree() <= 1 and sp.expand(c0t.coeff_monomial(1)) == 0
    const += c0t.coeff_monomial(th) * M_thC
    # L: <u g(ell) Lc> = B[-G Lc] + <G ell^2>, G = int g d ell, G(0) = 0
    pl = sp.collect(sp.expand(Lt), u, evaluate=False)
    assert list(pl) == [u], list(pl)
    g = pl[u]
    G = sp.integrate(g, lg)
    assert G.subs(lg, 0) == 0
    res["L derivative"] = sp.expand(r6_dth(-G * Lc) - (u * g * Lc - G * lg**2))
    const += fix(r6_bnd(-G * Lc)) + sum(sp.Poly(G, lg).coeff_monomial(lg**b) * Jv[b + 2] for b in range(4))
    # the rest: <l^b / sin^2> = -b J_{b-1}; <u th l^b> = B[-th l^{b+1}/(b+1)] + J_{b+1}/(b+1)
    for term in sp.Add.make_args(sp.expand(Rt)):
        c_, rest_ = term.as_independent(u, th, lg)
        d = rest_.as_powers_dict()
        a, k, b = d.get(u, 0), d.get(th, 0), d.get(lg, 0)
        if a == -2 and k == 0:
            assert b >= 1
            const += c_ * (-b * Jv[b - 1] - Jv[b])
        elif a == 0 and k == 0:
            const += c_ * Jv[b]
        elif a == 1 and k == 1:
            res.setdefault("rest: u th derivative", sp.Integer(0))
            res["rest: u th derivative"] += sp.expand(r6_dth(-th * lg ** (b + 1) / (b + 1)) - (u * th * lg**b - lg ** (b + 1) / (b + 1))) * c_
            const += c_ * (fix(r6_bnd(-th * lg ** (b + 1) / (b + 1))) + Jv[b + 1] / (b + 1))
        else:
            raise ValueError(("unhandled term of the rest", term))
    return sp.expand(const)


def r6_closed_symbolic(const):
    """G_6 at Omega_4 = 0 and Omega_4's closed form from the sum of the finite parts"""
    q, pi = R6_Q, sp.pi
    z3 = sp.zeta(3)
    w = R6_W
    wsub = {w[1]: 2 * pi**2 / 3 + 8 * q**2, w[2]: 80 * z3 - sp.Rational(8, 3) * pi**2 * q - sp.Rational(32, 3) * q**3,
            w[3]: 22 * pi**4 / 45 + sp.Rational(16, 3) * pi**2 * q**2 + sp.Rational(32, 3) * q**4 - 320 * z3 * q - 384 * z3}
    slope = sp.simplify(sp.diff(const, w[4]))
    g0 = sp.expand(const.subs(w[4], 0).subs(wsub))
    return slope, g0, sp.expand(g0 / (4 * pi))


def r6_numeric_integrand(Om):
    """the integrand of G_6 after E1 .. E4, from the closed forms of Y_0 .. Y_3, L_0 .. L_3 (order4.py, order5.py's L_3), Om = [w0 .. w4, 0]"""
    q = mp.ln2

    def f(th, x, tail):
        if tail:
            th = mp.pi / 2 - x
        c, lc, l, u, R, s2, Y, L, K, T = o4.closed(th, x, tail, Om)
        Cn = mp.clsin(2, mp.pi - 2 * th)
        sc = mp.sin(th) * mp.cos(th)
        B = -32 * (lc + q) * Cn + 16 * s2 * Cn + 32 * o4.Lc2_poly(th) + 16 * sc * lc * (lc + 2 * q) - 4 * Om[1] * th - 2 * Om[1] * sc
        LL = [L[0], L[1], L[2], -2 * hk(3, Om, l, R, Y, L) + 2 * u * B]
        o = list(Om)
        o[5] = 0
        fa = math.factorial
        t1 = 8 * lc * sum(Y[j] * cc(5 - j, Om, l, R) for j in range(4)) / s2
        t2 = 2 * sum(LL[5 - j] * l**j / fa(j) for j in range(2, 6))
        t3 = 2 * sum(Y[j] * cc(6 - j, o, l, R) for j in range(4)) / s2
        t4 = -4 * (lc - Om[0]) * sum(LL[4 - j] * l**j / fa(j) for j in range(1, 5))
        t5 = 8 * Om[0] * lc * (s2 * LL[3] + sum(Y[j] * cc(4 - j, Om, l, R) for j in range(4))) / s2
        F3 = sum(LL[3 - j] * l**j / fa(j) for j in range(4))
        t6 = -2 * (22 * lc**2 - 8 * Om[0] * lc - 2 * Om[1] + 16 * u * Cn) * F3
        return t1 + t2 + t3 + t4 + t5 + t6

    return o4.guarded(f, 80)


def check_reduce6(ev=None):
    t0 = time.time()
    lg, u, th, q = R6_LG, R6_U, R6_TH, R6_Q
    C, D3, Lc = R6_C, R6_D3, R6_LC
    mp.mp.dps = 45
    # RD0: the derivative rules of the generators, numerically
    t_ = mp.mpf("0.7")
    ell_f = lambda a: mp.log(mp.cos(a))  # noqa: E731
    C_f = lambda a: mp.clsin(2, mp.pi - 2 * a)  # noqa: E731
    D3_f = lambda a: mp.clcos(3, mp.pi - 2 * a) + 3 * mp.zeta(3) / 4  # noqa: E731
    Lc_f = o4.Lc2_poly
    errs = [abs(mp.diff(ell_f, t_) + mp.tan(t_)), abs(mp.diff(mp.tan, t_) - (1 + mp.tan(t_) ** 2)),
            abs(mp.diff(C_f, t_) - (2 * ell_f(t_) + 2 * mp.ln2)), abs(mp.diff(D3_f, t_) - 2 * C_f(t_)), abs(mp.diff(Lc_f, t_) - ell_f(t_) ** 2)]
    report("RD0", max(errs) <= mp.mpf("1e-30"), "derivative rules ell' = -u, u' = 1 + u^2, C' = 2 ell + 2 ln 2, D3' = 2 C, Lc2' = ell^2 at theta = 0.7 (numerical "
           f"derivatives): worst difference {nd(max(errs), 3)} (bar 1e-30)")
    Jv = r6_moments()
    parts = r6_parts()
    # RD1: the decomposition against the numeric integrand of closed forms at four points
    w1v, w2v, w3v = closed_omegas()[1], closed_omegas()[2], closed_omegas()[3]
    w4v = mp.mpf("0.3")
    Om = [-4 * mp.ln2, w1v, w2v, w3v, w4v, mp.mpf(0)]
    gens = (lg, u, th, C, D3, Lc, q) + tuple(R6_W[1:5])
    tot = sum(cf * C**a * D3**b * Lc**c for (a, b, c), cf in parts.items())
    tot_f = sp.lambdify(gens, tot, "mpmath")
    fn = r6_numeric_integrand(Om)
    worst = mp.mpf(99)
    for tv in ("0.2", "0.7", "1.1", "1.4"):
        tv = mp.mpf(tv)
        with mp.extradps(40):
            ell, s2, u_, Cn, D3n = o4.basics(tv, mp.pi / 2 - tv, False)
            dec = tot_f(ell, u_, tv, Cn, D3n, o4.Lc2_poly(tv), mp.ln2, w1v, w2v, w3v, w4v)
            num = fn(tv, mp.pi / 2 - tv, False)
        worst = min(worst, digits(num, dec))
    report("RD1", worst >= 30, f"the decomposition (sympy; coefficients of C^a D3^b Lc2^c in u, ell, theta: {sorted(parts)}) against the numeric integrand of closed forms "
           f"at theta = 0.2, 0.7, 1.1, 1.4: worst {float(worst):.1f} digits (bar 30)")
    # the reduction (finite parts)
    R6_RES.clear()
    const = r6_reduce(parts, Jv)
    bad = {k: v for k, v in R6_RES.items() if sp.simplify(v) != 0}
    report("RD2", not bad, f"the by-parts identities, symbolically ({len(R6_RES)}: " + "; ".join(sorted(R6_RES)) + f"): nonzero residuals {sorted(bad) or 'none'}")
    slope, g0, om4c = r6_closed_symbolic(const)
    z3sq = sp.expand(const).coeff(sp.zeta(3), 2)
    report("RD3", sp.simplify(slope + 4 * sp.pi) == 0 and sp.simplify(z3sq) == 0 and not bad.get("C: u^-1 part"),
           f"cancellations: slope of G_6 in Omega_4 = {slope} (= -4 pi); the 1/u terms of the C coefficient cancel ({sp.simplify(R6_RES['C: u^-1 part'])}); "
           f"J_6 drops out (the zeta(3)^2 coefficient is {sp.simplify(z3sq)})")
    expr = sp.expand(om4c)
    print("  Omega4 (derived) = " + str(expr))
    cl = sp.lambdify((q,), expr, "mpmath")
    mp.mp.dps = 330
    om4num = cl(mp.ln2)
    print(f"  Omega4 (derived) = {nd(om4num, 60)}")
    # RD4: the series of the finite-part kit against the functions
    kit_pi2, kit_zero = (_KIT["pi2"], _KIT["zero"]) if _KIT else r6_kit()
    mp.mp.dps = 60
    xv = mp.mpf(1) / 1000
    subs = {R6_X: xv, R6_L: mp.log(xv), R6_Q: mp.ln2, R6_Z3: mp.zeta(3), R6_J2: mp.pi * (12 * mp.ln2**2 + mp.pi**2) / 24}
    th_v = mp.pi / 2 - xv
    ref = {"ell": ell_f(th_v), "u": mp.tan(th_v), "ui": 1 / mp.tan(th_v), "C": C_f(th_v), "D3": D3_f(th_v), "Lc": Lc_f(th_v)}
    e1 = max(abs(mp.mpf(sp.N(kit_pi2[k].subs(subs), 70)) - ref[k]) / max(abs(ref[k]), mp.mpf(1)) for k in ref)
    tv0 = mp.mpf(1) / 1000
    subs0 = {R6_T: tv0, R6_Q: mp.ln2}
    ref0 = {"ell": ell_f(tv0), "u": mp.tan(tv0), "ui": 1 / mp.tan(tv0), "C": C_f(tv0), "D3": D3_f(tv0), "Lc": Lc_f(tv0)}
    e2 = max(abs(mp.mpf(sp.N(kit_zero[k].subs(subs0), 70)) - ref0[k]) / max(abs(ref0[k]), mp.mpf(1e-30)) for k in ref0)
    report("RD4", max(e1, e2) <= mp.mpf("1e-24"), f"the asymptotic series of the generators (order {R6_ORD}) against the functions at x = pi/2 - theta = 1e-3 (worst relative "
           f"difference {nd(e1, 3)}) and at theta = 1e-3 ({nd(e2, 3)}) (bar 1e-24)")
    # RD5: the masters, numerically
    mp.mp.dps = 40
    pi_, q_, z3_, z5_ = mp.pi, mp.ln2, mp.zeta(3), mp.zeta(5)

    def integ(fn_):
        return o4.both(o4.guarded(fn_, 40))

    def bs(th_, x, tail):
        if tail:
            th_ = mp.pi / 2 - x
        ell, s2, u_, Cn, D3n = o4.basics(th_, x, tail)
        return th_, ell, u_, Cn, D3n

    Jn = [mp.mpf(sp.N(j.subs({R6_Q: sp.log(2)}), 50)) for j in Jv]
    masters = {
        "C^2 = pi^5/360": (integ(lambda a, x, t: bs(a, x, t)[3] ** 2), pi_**5 / 360),
        "ell C^2 = -ln 2 pi^5/360": (integ(lambda a, x, t: bs(a, x, t)[1] * bs(a, x, t)[3] ** 2), -q_ * pi_**5 / 360),
        "theta C = pi zeta(3)/4": (integ(lambda a, x, t: bs(a, x, t)[0] * bs(a, x, t)[3]), pi_ / 4 * z3_),
        "ell^2 D3": (integ(lambda a, x, t: bs(a, x, t)[1] ** 2 * bs(a, x, t)[4]),
                     mp.mpf(3) / 4 * z3_ * Jn[2] + pi_ / 4 * (2 * q_ * pi_**4 / 90 + (4 if MUTATE == "rd_master" else 5) * z5_ - 2 * pi_**2 / 6 * z3_)),
        "C Lc2": (integ(lambda a, x, t: bs(a, x, t)[3] * o4.Lc2_poly(bs(a, x, t)[0])), mp.mpf(7) / 8 * z3_ * Jn[2]
                  - (mp.mpf(3) / 4 * z3_ * Jn[2] + pi_ / 4 * (2 * q_ * pi_**4 / 90 + (4 if MUTATE == "rd_master" else 5) * z5_ - 2 * pi_**2 / 6 * z3_)) / 2)}
    for b in range(4):
        masters[f"u ell^{b} C = (2 J_{b + 2} + 2 ln2 J_{b + 1})/{b + 1}"] = (integ(lambda a, x, t, b=b: bs(a, x, t)[2] * bs(a, x, t)[1] ** b * bs(a, x, t)[3]),
                                                                              (2 * Jn[b + 2] + 2 * q_ * Jn[b + 1]) / (b + 1))
    for n in range(7):
        masters[f"J_{n}"] = (integ(lambda a, x, t, n=n: bs(a, x, t)[1] ** n), Jn[n])
    dm = {k: digits(a, b) for k, (a, b) in masters.items()}
    report("RD5", min(dm.values()) >= 30, "the master integrals against their closed forms (quadrature at 40 digits), digits: "
           + ", ".join(f"{k}: {float(v):.1f}" for k, v in dm.items()) + " (bar 30)")
    # RD6: method A at order 6, the numeric integrand of closed forms, Omega_4 from the affine G_6, against the derived closed form
    cl = closed_omegas()
    w_ = [cl[0], cl[1], cl[2], cl[3]]

    def G(o4v):
        f = r6_numeric_integrand(w_ + [o4v, mp.mpf(0)])
        return (mp.quad(lambda t: f(t, mp.pi / 2 - t, False), [0, mp.mpf(1) / 2, 1], method="gauss-legendre")
                + mp.quad(lambda x: f(mp.pi / 2 - x, x, True), [0, mp.pi / 2 - 1]))

    t1 = time.time()
    g0n, g1n = G(mp.mpf(0)), G(mp.mpf(1))
    slope_n, om4_a = g1n - g0n, -g0n / (g1n - g0n)
    mp.mp.dps = 60
    ds, dc = digits(-4 * mp.pi, slope_n), digits(om4num, om4_a)
    report("RD6", ds >= 30 and dc >= 30, f"method A at order 6 (closed forms, E1 .. E4, quadrature at 40 digits, {time.time() - t1:.0f} s): G_6(Omega_4 = 0) = {nd(g0n, 30)}, "
           f"slope against -4 pi {float(ds):.1f} digits; Omega_4 = {nd(om4_a, 38)} against the derived closed form {float(dc):.1f} digits (bar 30)")
    # RD7: the derived closed form against the grids of record
    if ev is not None:
        mp.mp.dps = 330
        rows = {n: digits(om4num, ev[n]["Om4"]) for n in ev}
        report("RD7", all(rows[n] >= 100 for n in ("G1", "G2") if n in rows), "Omega_4 of the derived closed form against the grids' Omega_4, digits: "
               + ", ".join(f"{n} {float(d):.1f}" for n, d in rows.items()) + " (bar 100 on G1, G2; the others reported)")
    # RD8: the integer relation PSLQ found in run_05 (raw I220-rpi5, P2 and P3, residual 6.4e-196) is 45 x the derived closed form
    ps, zs3, zs5, qs_ = sp.symbols("ps zs3 zs5 qs_", positive=True)
    poly = sp.Poly(sp.expand(45 * expr.subs({sp.pi: ps, sp.zeta(3): zs3, sp.zeta(5): zs5, R6_Q: qs_})), qs_, ps, zs3, zs5)
    monos = [(1, 4, 0, 0), (3, 2, 0, 0), (5, 0, 0, 0), (0, 2, 1, 0), (2, 0, 1, 0), (0, 0, 0, 1), (0, 4, 0, 0), (2, 2, 0, 0), (4, 0, 0, 0), (1, 0, 1, 0), (0, 0, 1, 0),
             (1, 2, 0, 0), (3, 0, 0, 0), (0, 2, 0, 0), (2, 0, 0, 0), (1, 0, 0, 0), (0, 0, 0, 0)]       # q^a pi^b zeta(3)^c zeta(5)^d of P2's constants
    pslq_rel = [-45, -88, -320, -384, -7200, 28800, -469440, -2432, 0, 0, 69120, 138240, 0, 0, 0, 0, 0, 0]    # run_05, P2: [Omega4, constants...]
    mine = [-45] + [poly.coeff_monomial(m) for m in monos]
    report("RD8", [int(v) for v in mine] == pslq_rel and set(poly.monoms()) <= set(monos),
           f"the integer relation of run_05 (PSLQ, P2 and P3) {pslq_rel} against -45 followed by 45 x the derived closed form's coefficients {[int(v) for v in mine]} "
           f"(every term of the closed form is a constant of P2)")
    print(f"  ({time.time() - t0:.0f} s)")
    return expr


def main():
    parse(sys.argv[1:])
    mp.mp.dps = 60
    print(f"order6.py  mutate={MUTATE}  parts={PARTS}  mpmath {mp.__version__}  sympy {sp.__version__}", flush=True)
    t0 = time.time()
    if "grid" in PARTS:
        check_grid()
    if "sym" in PARTS:
        check_sym()
    ev = None
    if "evals" in PARTS:
        ev = check_evals()
        if "values" in PARTS or not GRIDSEL:
            check_values(ev)
    elif LOAD and any(p in PARTS for p in ("values", "pslq", "remainder", "reduce6")):
        mp.mp.dps = 430
        ev = load_ev(LOAD)
        if "values" in PARTS:
            check_values(ev)
    if "exch6" in PARTS:
        check_exch6()
    if "reduce6" in PARTS:
        check_reduce6(ev)
    if "pslq" in PARTS:
        check_pslq(ev)
    if "remainder" in PARTS:
        check_remainder(ev)
    nfail = sum(1 for _, ok in RESULTS if not ok)
    print(f"\n{len(RESULTS) - nfail} passed, {nfail} failed of {len(RESULTS)} checks; {time.time() - t0:.0f} s")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()

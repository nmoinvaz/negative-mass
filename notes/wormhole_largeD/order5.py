"""Order 5 of the matched expansion (matched.py, order3.py, order4.py): Omega3 from G5 = 0 by method B, and the checks of issue #208.

W^2 a^2 = s + Omega, s = 2(D - 3), eps = 1/s, Omega = Omega0 + Omega1 eps + Omega2 eps^2 + Omega3 eps^3 + ..., R = r^s, l = ln R.
The order equations (E0 derives them from the master equation), the same for every k:
  dY_k/dl = F_{k-1},  F_{k-1} = sum_{j<k} L_{k-1-j} l^j / j!,   (R - 1) L_k' - L_k/2 = h_k,
  h_k = (1 - 1/R) L_{k-1} + sum_{j<k} Y_j c_{k-j} + Y_k/2,   L_k = -2 Omega_{k-1} + u C_k,   G_k = C_k(inf) = 0,
with u = sqrt(R - 1), C_k = int_1^R (h_k - Omega_{k-1}) (R' - 1)^{-3/2} dR'. G_k is affine in Omega_{k-2} with slope -4 pi:
G2 = 0 gives Omega0, G3 = 0 Omega1, G4 = 0 Omega2, **G5 = 0 gives Omega3**. One path computes them all (`omega_chain`).

Method B (order4.py's PanelB, shared): v = sqrt(l), panels of width w on [0, V], N Chebyshev nodes per panel, cumulative integrals by the
spectral integration matrix; Y_k from the left, the tail C_right from the right, L_k = -2 h_k(1) - u C_right, G_k = C_right(0). It uses no
closed form of Y_k or L_k and no input value of Omega_0 .. Omega_2.

    python3 order5.py [sym] [evals] [pslq] [remainder] [closedform] [reduce5] [a5] [--mutate=NAME] [--w2=w2_output.txt] [--w2raw=FILE] [--save=FILE] [--load=FILE] [--no-real]
    python3 order5.py grid --grid=N,w,V,dps [--kmax=5]      (a development run of one grid; no gate)

No argument runs sym, evals, pslq, remainder. The exit status is 1 when a check fails. Status words: sympy residual 0 is "checked
symbolically", an agreement of integrals "numerical evidence", a PSLQ relation an "identification"; a spread between grids or precisions
is an estimator range, not a bound.
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

MUTATIONS = ("wrong_h", "wrong_c5", "wrong_Om0", "wrong_Om1", "wrong_Om2", "coarse_B", "short_tail", "pslq_wrong_basis",
             "wrong_Om3_series", "wrong_closed", "wrong_E3", "wrong_AD", "wrong_C2", "wrong_reg", "wrong_J4")
MUTATE = None
W2FILE = "w2_output.txt"
SAVE = None
LOAD = None
GRIDARG = None
W2RAW = None
KMAX = 5
PARTS = []
REAL = True


def parse(argv):
    global MUTATE, W2FILE, W2RAW, SAVE, LOAD, GRIDARG, KMAX, PARTS, REAL
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
                + (-1 if MUTATE == "wrong_c5" else 1) * l**4 / (24 * R))
    raise ValueError(m)


def hk(m, Om, l, R, Y, L):
    """h_m"""
    if m == 0:
        return Y[0] / 2
    return (1 - 1 / R) * L[m - 1] + sum(Y[j] * cc(m - j, Om, l, R) for j in range(m)) + Y[m] / 2


def Fsum(k, L, l):
    return sum(L[k - 1 - j] * l**j / math.factorial(j) for j in range(k))


# ---------------------------------------------------------------------------------------------------------------
# E0: the epsilon-expansion of the master equation against the table, k = 0 .. 5
# ---------------------------------------------------------------------------------------------------------------
def check_sym():
    t0 = time.time()
    R, eps, l = sp.symbols("R epsilon l", positive=True)
    Om = sp.symbols("Omega0:5")
    Ys, Ls = sp.symbols("Y0:6"), sp.symbols("L0:6")
    kk = 5
    r = sum((eps * l) ** j / sp.factorial(j) for j in range(kk + 1))          # exact through eps^5, which is all that is read
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
           f"h_k (k = 0 .. 5) and dY_k/dl (k = 1 .. 5) from the epsilon-expansion of the master equation against the table: "
           f"residuals {[str(x) for x in res if x != 0] or '0'} ({time.time() - t0:.1f} s)")


# ---------------------------------------------------------------------------------------------------------------
# method B at order k
# ---------------------------------------------------------------------------------------------------------------
class PanelK(PanelB):
    def hierarchy(self, Om, kmax):
        """Y_0 .. Y_kmax, L_0 .. L_kmax and G_0 .. G_kmax for the given Omega list (Om[kmax - 1] is not used: h_kmax(1) = Om[kmax - 1])"""
        N, P = self.N, self.P
        grid = lambda fn: [[fn(p, i) for i in range(N)] for p in range(P)]  # noqa: E731
        Y = [grid(lambda p, i: mp.mpf(1))]
        L = []
        G = {}
        for k in range(kmax + 1):
            if k >= 1:
                F = grid(lambda p, i: 2 * self.v[p][i] * sum(L[k - 1 - j][p][i] * self.l[p][i] ** j / mp.factorial(j)
                                                             for j in range(k)))
                Y.append(self.cum_left(F))
            thr = mp.mpf(1) / 2 if k == 0 else Om[k - 1]
            s = grid(lambda p, i: hk(k, Om, self.l[p][i], self.R[p][i], [y[p][i] for y in Y],
                                     [lv[p][i] for lv in L] if k else None) - thr)
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


def omega_chain(B, kmax, log=None):
    """Omega_0 .. Omega_{kmax-2} from G_2 .. G_kmax = 0, in turn, by the same code at every order; returns (Om, slopes, seconds per order)"""
    Om = [mp.mpf(0)] * kmax
    slopes, secs = {}, {}
    for k in range(2, kmax + 1):
        t0 = time.time()
        Om_in = list(Om)
        if MUTATE == "wrong_Om0" and k >= 3:
            Om_in[0] = -4 * mp.ln2 + mp.mpf(1) / 100
        if MUTATE == "wrong_Om1" and k >= 4:
            Om_in[1] = 2 * mp.pi**2 / 3 + 8 * mp.ln2**2 + mp.mpf(1) / 100
        if MUTATE == "wrong_Om2" and k >= 5:
            Om_in[2] = (80 * mp.zeta(3) - mp.mpf(8) / 3 * mp.pi**2 * mp.ln2 - mp.mpf(32) / 3 * mp.ln2**3) + mp.mpf(1) / 100
        val, slope = solve_omega(lambda kk, o: B.hierarchy(o, kk)[2][kk], k, Om_in)
        Om[k - 2] = val
        slopes[k] = slope
        secs[k] = time.time() - t0
        if log:
            print(f"    k = {k}: Omega_{k - 2} = {nd(val, 30)}  slope {nd(slope, 20)}  ({secs[k]:.0f} s)", flush=True)
    return Om, slopes, secs


def run_B(N, w, V, dps, kmax=5, log=False):
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


# ---------------------------------------------------------------------------------------------------------------
# the closed forms of the planted orders, #195's fit
# ---------------------------------------------------------------------------------------------------------------
def closed_omegas():
    return [-4 * mp.ln2, 2 * mp.pi**2 / 3 + 8 * mp.ln2**2,
            80 * mp.zeta(3) - mp.mpf(8) / 3 * mp.pi**2 * mp.ln2 - mp.mpf(32) / 3 * mp.ln2**3]


# declared grids: name -> (N, w, V, working digits). B3, B4 are the grids of record; B1, B2 are order4.py's (#26), reported only.
GRIDS = {"B3": (80, (3, 10), 28, 180), "B4": (90, (1, 4), 30, 240)}              # w as a fraction (num, den), evaluated at the working precision
REPORT_GRIDS = {"B1": (40, (1, 2), 18, 120), "B2": (50, (2, 5), 20, 140)}

OM3_ISSUE = "-652.8401075700561254352393037996"                      # the issue's string (#195's fit, 28 decimals)
OM3_FIT40 = "-652.8401075700561254352393037996473610805"             # raw I195-dh610/dev_25_pslq_omega3.log, line 1; range 2.74e-35


def read_fit(path):
    """Omega4 of #195's fit (w2_output.txt); Omega3 there has 14 digits only, so the longer strings above are used for it"""
    txt = open(path).read()
    return mp.mpf(re.search(r"Omega4 = (\S+)", txt).group(1))


def save_ev(ev, path):
    with open(path, "w") as f:
        json.dump({n: {k: (mp.nstr(v, 330) if isinstance(v, mp.mpf) else v) for k, v in e.items()} for n, e in ev.items()}, f, indent=1)


def load_ev(path):
    raw = json.load(open(path))
    return {n: {k: (mp.mpf(v) if k.startswith(("Om", "slope")) else v) for k, v in e.items()} for n, e in raw.items()}


def check_grid():
    """one grid, no gate: the development run"""
    N, w, V, dps = GRIDARG.split(",")
    mp.mp.dps = 60
    e = run_B(int(N), tuple(int(t) for t in w.split("/")) if "/" in w else w, int(V), int(dps), kmax=KMAX, log=True)
    mp.mp.dps = int(dps)
    cl = closed_omegas()
    fit = mp.mpf(OM3_FIT40)
    print(f"  {e['grid']}  {e['secs']:.0f} s")
    for i in range(KMAX - 1):
        print(f"    Omega{i} = {nd(e[f'Om{i}'], int(dps))}")
    for i in range(3):
        if f"Om{i}" in e:
            print(f"    Omega{i} against closed form: {float(digits(cl[i], e[f'Om{i}'])):.1f} digits")
    if "Om3" in e:
        print(f"    Omega3 against #195's fit: difference {nd(e['Om3'] - fit, 5)}")
    for k in range(2, KMAX + 1):
        print(f"    slope G{k} against -4 pi: {float(digits(-4 * mp.pi, e[f'slope{k}'])):.1f} digits")
    return e


# ---------------------------------------------------------------------------------------------------------------
# the evaluations and their checks
# ---------------------------------------------------------------------------------------------------------------
def check_evals():
    mp.mp.dps = 60
    ev = {}
    for name, (N, w, V, dps) in list(GRIDS.items()) + list(REPORT_GRIDS.items()):
        ev[name] = run_B(N, w, V, dps)
        e = ev[name]
        print(f"  {name}: {e['grid']}  {e['secs']:.0f} s\n    Omega0 = {nd(e['Om0'], 60)}\n    Omega1 = {nd(e['Om1'], 60)}"
              f"\n    Omega2 = {nd(e['Om2'], 60)}\n    Omega3 = {nd(e['Om3'], 110)}\n    slopes G2 .. G5 = "
              + " ".join(nd(e[f"slope{k}"], 20) for k in range(2, 6)), flush=True)
    if SAVE:
        save_ev(ev, SAVE)
    return ev


def check_values(ev):
    mp.mp.dps = 330
    cl = closed_omegas()
    grids = list(GRIDS)
    # PL0 .. PL2: the planted orders by the same code
    for i, cid in enumerate(("PL0", "PL1", "PL2")):
        rows = {n: digits(cl[i], ev[n][f"Om{i}"]) for n in ev}
        ok = all(rows[n] >= 40 for n in grids)
        report(cid, ok, f"planted: Omega{i} against its closed form, digits: " + ", ".join(f"{n} {float(d):.1f}" for n, d in rows.items())
               + " (bar 40 on B3, B4; B1, B2 reported)")
    # PL3: a wrong input must move them (the control's own statement, evaluated when a wrong-input mutation is on)
    if MUTATE in ("wrong_Om0", "wrong_Om1"):
        i = 1 if MUTATE == "wrong_Om0" else 2
        moved = all(abs(ev[n][f"Om{i}"] - cl[i]) > mp.mpf(1) / 1000 for n in grids)
        report("PL3", moved, f"control {MUTATE} moves Omega{i} of " + ", ".join(f"{n} by {nd(ev[n][f'Om{i}'] - cl[i], 5)}" for n in grids))
    # SL: slopes
    rows, ok = [], True
    for n in ev:
        ds = [digits(-4 * mp.pi, ev[n][f"slope{k}"]) for k in range(2, 6)]
        rows.append(f"{n} " + "/".join(f"{float(d):.1f}" for d in ds))
        ok = ok and (n not in grids or min(ds) >= 40)
    report("SL", ok, "slope of G_k in Omega_(k-2) against -4 pi, k = 2 .. 5, digits: " + ", ".join(rows) + " (bar 40 on B3, B4)")
    # M3: the two grids on Omega3 (and Omega1, Omega2 reported)
    a, b = grids
    d3 = digits(ev[a]["Om3"], ev[b]["Om3"])
    report("M3", d3 >= 30, f"Omega3 of {a} against {b} (two grids, two precisions): {float(d3):.1f} digits (bar 30); |diff| = "
           f"{nd(abs(ev[a]['Om3'] - ev[b]['Om3']), 3)}")
    print("  reported: digits between evaluations, Omega1 / Omega2 / Omega3: " + ", ".join(
        f"{x}-{y} " + "/".join(f"{float(digits(ev[x][f'Om{i}'], ev[y][f'Om{i}'])):.1f}" for i in (1, 2, 3))
        for x, y in (("B3", "B4"), ("B1", "B2"), ("B2", "B3"), ("B2", "B4"))))
    # X1: Omega3 against #195's fit
    issue, fit40 = mp.mpf(OM3_ISSUE), mp.mpf(OM3_FIT40)
    dd = {n: abs(ev[n]["Om3"] - issue) for n in grids}
    report("X1", all(d <= mp.mpf("1e-30") for d in dd.values()),
           "Omega3 against the issue's string " + OM3_ISSUE + " (28 decimals, the truncation of #195's median: no Omega3 can meet the bar): "
           + ", ".join(f"{n} {nd(d, 3)}" for n, d in dd.items()) + " (bar 1e-30)")
    d40 = {n: abs(ev[n]["Om3"] - fit40) for n in ev}
    report("X1b", all(d40[n] <= mp.mpf("1e-30") for n in grids),
           "Omega3 against #195's 40-digit median " + OM3_FIT40 + " (range 2.74e-35), amendment A1: " + ", ".join(f"{n} {nd(d, 3)}" for n, d in d40.items())
           + " (bar 1e-30 on B3, B4; B1, B2 reported)")
    return d3


# ---------------------------------------------------------------------------------------------------------------
# PSLQ
# ---------------------------------------------------------------------------------------------------------------
PSLQ_DIGITS = 110             # cap on the working digits; the digits the two grids share if fewer
PSLQ_SLACK = 12               # tolerance 10^-(d - slack)
# basis name -> (maxcoeff H, planted seed)
PSLQ_BASES = {"P1": (10**9, 208001), "P2": (10**6, 208002), "P3": (10**6, 208003)}


def constants(name):
    pi, ln2, z3 = mp.pi, mp.ln2, mp.zeta(3)
    li4 = mp.polylog(4, mp.mpf(1) / 2)
    G, b4 = mp.catalan, mp.dirichlet(4, [0, 1, 0, -1])
    p1 = [("pi^4", pi**4), ("pi^2 ln^2 2", pi**2 * ln2**2), ("ln^4 2", ln2**4), ("zeta(3) ln 2", z3 * ln2), ("Li4(1/2)", li4)]
    if name == "P1":
        return p1
    if name == "P2":
        return p1 + [("zeta(3)", z3), ("pi^2 ln 2", pi**2 * ln2), ("ln^3 2", ln2**3), ("pi^2", pi**2), ("ln^2 2", ln2**2), ("ln 2", ln2),
                     ("1", mp.mpf(1))]
    if name == "P3":
        return p1 + [("pi zeta(3)", pi * z3), ("pi^3 ln 2", pi**3 * ln2), ("pi ln^3 2", pi * ln2**3), ("pi^2 G", pi**2 * G), ("beta(4)", b4),
                     ("G ln^2 2", G * ln2**2), ("pi G ln 2", pi * G * ln2)]
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
    grids = list(GRIDS)
    mp.mp.dps = 330
    agree = digits(ev[grids[0]]["Om3"], ev[grids[1]]["Om3"])
    d = int(max(30, min(PSLQ_DIGITS, math.floor(agree))))      # PSLQ needs at least 53 bits
    pert = mp.sqrt(2) / 1000        # irrational, so that the shifted constant is outside every basis's rational span (1 is in P2)
    print(f"  PSLQ at {d} digits (the two grids share {float(agree):.1f}), tolerance 1e-{d - PSLQ_SLACK}")
    for name, (H, seed) in PSLQ_BASES.items():
        with mp.workdps(d):
            tau = mp.mpf(10) ** (-(d - PSLQ_SLACK))
            names, cs = zip(*constants(name))
            cs_true = list(cs)
            cs = list(cs)
            if MUTATE == "pslq_wrong_basis":                  # PSLQ is given a basis that is not the one the target was built from
                cs[0] = cs[0] + pert
            n = len(cs)
            budget = (n + 1) * math.log10(2 * H + 1) + 1
            om3 = +ev[grids[1]]["Om3"]
            # planted truth: den * x = sum c_i b_i, through the same call
            den, c = planted(n, seed)
            xp = sum(ci * b for ci, b in zip(c, cs_true)) / den
            rp, tp = pslq_run([xp] + cs, tau, H)
            want = [den] + [-ci for ci in c]
            okp = rp is not None and (list(rp) == want or list(rp) == [-v for v in want])
            # perturbed basis on the planted target: must find nothing
            csp = list(cs)
            csp[1] = csp[1] + pert
            rpp, _ = pslq_run([xp] + csp, tau, H)
            report("Q" + name[1], okp and rpp is None,
                   f"{name} = {{Omega3, {', '.join(names)}}}, n = {n + 1}, H = {H}, d = {d}: planted relation (den {den}, seed {seed}) "
                   f"recovered: {okp} ({tp:.0f} s); perturbed basis on the planted target finds {'nothing' if rpp is None else list(rpp)}; "
                   f"budget n log10(2H + 1) + 1 = {budget:.1f}")
            # the real target
            if not REAL:
                print(f"    REAL {name}: skipped (--no-real)", flush=True)
                continue
            rr, tr = pslq_run([om3] + cs, tau, H)
            found, why = False, "not found"
            if rr is not None:
                res = residual(rr, [om3] + cs)
                found = res <= tau and budget < d and rr[0] != 0      # a relation without Omega3 is among the constants alone
                why = f"relation {list(rr)}, residual {nd(res, 3)} (bar {nd(tau, 2)}), budget with H = {H}: {budget:.1f} < {d}: {budget < d}; Omega3 coefficient {rr[0]}"
                if found:
                    j = next(i for i in range(1, len(rr)) if rr[i] != 0)
                    csr = list(cs)
                    csr[j - 1] = csr[j - 1] + pert
                    rrp, _ = pslq_run([om3] + csr, tau, H)
                    found = found and rrp is None
                    why += f"; perturbed basis (constant {names[j - 1]} + sqrt(2)/1000) finds {'nothing' if rrp is None else list(rrp)}"
            print(f"    REAL {name}: " + (("FOUND (counts) " if found else "found but does not count ") + why if rr is not None
                                           else f"not found at H = {H} ({tr:.0f} s)"), flush=True)


# ---------------------------------------------------------------------------------------------------------------
# the remainder against #195's W^2
# ---------------------------------------------------------------------------------------------------------------
def closed_omega3():
    """the P2 relation of run_02 (A2): 45 Omega3 - 22 pi^4 - 240 pi^2 ln^2 2 - 480 ln^4 2 + 14400 zeta(3) ln 2 + 17280 zeta(3) = 0"""
    pi, ln2, z3 = mp.pi, mp.ln2, mp.zeta(3)
    return ((23 if MUTATE == "wrong_closed" else 22) * pi**4 + 240 * pi**2 * ln2**2 + 480 * ln2**4 - 14400 * z3 * ln2 - 17280 * z3) / 45


def check_closedform(ev):
    """X3 (A2): Omega3 of B3 and B4 against 22 pi^4/45 + 16 pi^2 ln^2 2/3 + 32 ln^4 2/3 - 320 zeta(3) ln 2 - 384 zeta(3)"""
    mp.mp.dps = 330
    cf = closed_omega3()
    rows = {n: digits(cf, ev[n]["Om3"]) for n in ev}
    fit40 = mp.mpf(OM3_FIT40)
    report("X3", all(rows[n] >= 100 for n in GRIDS),
           f"Omega3 against 22 pi^4/45 + 16 pi^2 ln^2 2/3 + 32 ln^4 2/3 - 320 zeta(3) ln 2 - 384 zeta(3) = {nd(cf, 50)}, digits: "
           + ", ".join(f"{n} {float(d):.1f}" for n, d in rows.items()) + " (bar 100 on B3, B4; B1, B2 reported); against #195's 40-digit median: "
           f"{nd(abs(cf - fit40), 3)} (its range 2.74e-35)")


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


def check_remainder(ev):
    mp.mp.dps = 80
    grids = list(GRIDS)
    om1, om2, om3 = ev[grids[1]]["Om1"], ev[grids[1]]["Om2"], ev[grids[1]]["Om3"] + (mp.mpf(1) / 100 if MUTATE == "wrong_Om3_series" else 0)
    om4 = read_fit(W2FILE)
    lst = read_w2(W2FILE)
    if W2RAW:
        lst.update(read_w2raw(W2RAW))      # D = 1000, 2000 are not in the 64-value list
    x = {}
    print("  five-term series s - 4 ln 2 + Omega1/s + Omega2/s^2 + Omega3/s^3 against W^2 (qd, w2_output.txt), s = 2(D - 3):")
    print("  D       s       remainder              rem*s^4           (rem*s^4 - Omega4)")
    for D in (1000, 2000, 5103, 10000):
        s = mp.mpf(2 * (D - 3))
        rem = mp.mpf(lst[D]) - (s - 4 * mp.ln2 + om1 / s + om2 / s**2 + om3 / s**3)
        x[D] = rem * s**4
        print(f"  {D:<6d} {int(s):<7d} {nd(rem, 12):<22s} {nd(x[D], 12):<17s} {nd(x[D] - om4, 6)}")
    dev = [abs(x[D] - om4) for D in (1000, 2000, 5103, 10000)]
    falls = all(dev[i] > dev[i + 1] for i in range(3))
    rel = dev[3] / abs(om4)
    report("R3", falls and rel < mp.mpf("5e-4"), f"|rem*s^4 - Omega4| at D = 1000, 2000, 5103, 10^4: {', '.join(nd(v, 5) for v in dev)} "
           f"(falling: {falls}); at D = 10^4 {nd(100 * rel, 3)} % of |Omega4| = {nd(abs(om4), 11)} (bar 0.05 %)")
    s1, s2 = mp.mpf(2 * (5103 - 3)), mp.mpf(2 * (10000 - 3))
    est = (s2 * x[10000] - s1 * x[5103]) / (s2 - s1)
    report("R4", abs(est - om4) <= mp.mpf("0.1"), f"two-point extrapolation (s2 x2 - s1 x1)/(s2 - s1) from D = 5103, 10^4: {nd(est, 12)}, Omega4 = "
           f"{nd(om4, 14)}: difference {nd(est - om4, 4)} (bar 0.1)")
    print("  reported, not gated: rem*s^4 at the 64 values of D and D = 200 .. 2000 (60-digit qd), and (rem*s^4 - Omega4)*s, which tends to Omega5")
    print("  D       rem*s^4             (rem*s^4 - Omega4)*s")
    for D in sorted(lst):
        s = mp.mpf(2 * (D - 3))
        xx = (mp.mpf(lst[D]) - (s - 4 * mp.ln2 + om1 / s + om2 / s**2 + om3 / s**3)) * s**4
        print(f"  {D:<7d} {nd(xx, 14):<19s} {nd((xx - om4) * s, 8)}")


# ---------------------------------------------------------------------------------------------------------------
# the reduction of G_5 by hand (RD0 .. RD8) and method A at order 5 (A5, RD9)
# ---------------------------------------------------------------------------------------------------------------
def e3f():
    return 2 if MUTATE == "wrong_E3" else 4


def sym_pieces():
    """sympy: the Y_4-free integrand I_5 at Omega0 = -4 q in the generators ell = ln cos, S = sin^2, u = tan, C = Cl2(pi - 2 th),
    D3 = Cl3(pi - 2 th) + 3 zeta(3)/4, Lc2 = int_0^th ln^2 cos, th; I_5 = rest + cC C + cD D3 + cL Lc2"""
    lg, S, u, C, D3, Lc, th, q, K = sp.symbols("ell S u C D3 Lc2 th q K")
    w0, w1, w2, w3, w4 = sp.symbols("w0 w1 w2 w3 w4")
    R = 1 / (1 - S)
    l = -2 * lg
    sc = u * (1 - S)
    L = [-1, -2 * w0 + 4 * lg, -2 * w1 - 4 * w0 * lg + 16 * lg**2 + 16 * u * C]
    Y = [1, 2 * lg, 4 * w0 * lg - 6 * lg**2, 4 * w1 * lg - 4 * lg**3 + 32 * lg**2 - 16 * w0 * lg + 32 * u * C - 16 * D3]
    Om = [w0, w1, w2, w3, w4]
    h3 = hk(3, Om, l, R, Y + [0], L)
    B = -32 * (lg + q) * C + 16 * S * C + 32 * Lc + 16 * sc * lg * (lg + 2 * q) - 4 * w1 * th - 2 * w1 * sc
    LL = L + [-2 * h3 + 2 * u * B]
    h4 = S * LL[3] + sum(Y[j] * cc(4 - j, Om, l, R) for j in range(4))            # Y_4/2 removed by E3
    t1 = 8 * lg * h4 / S
    t2 = 2 * sum(LL[3 - j] * l ** (j + 1) / sp.factorial(j + 1) for j in range(4))
    o = list(Om)
    o[4] = 0
    t3 = 2 * sum(Y[j] * cc(5 - j, o, l, R) for j in range(4)) / S
    F3 = sum(LL[3 - j] * l**j / sp.factorial(j) for j in range(4))
    t4 = -e3f() * (lg - w0) * F3

    def red_u(e):
        p = sp.Poly(sp.expand(e), u)
        return sp.expand(sum(cf * (S / (1 - S)) ** (k // 2) * u ** (k % 2) for (k,), cf in p.terms()))

    I5 = red_u((t1 + t2 + t3 + t4).subs(w0, -4 * q))
    d = dict(lg=lg, S=S, u=u, C=C, D3=D3, Lc=Lc, th=th, q=q, K=K, w0=w0, w1=w1, w2=w2, w3=w3, w4=w4, I5=I5, red_u=red_u,
             cC=sp.together(sp.Poly(I5, C).coeff_monomial(C)), cD=sp.together(sp.Poly(I5, D3).coeff_monomial(D3)),
             cL=sp.together(sp.Poly(I5, Lc).coeff_monomial(Lc)), rest=sp.expand(sp.together(I5.subs({C: 0, D3: 0, Lc: 0}))))
    # D3 by parts: cD = 32 (6 ell - 8 q) + cD1/S; A_D = int cD with K = 2 q th - C, K' = -2 ell, F_n = int ell^n/S
    cD, cot = d["cD"], 1 / u
    cDe = sp.expand(sp.simplify(cD * S))
    cD0, cD1 = cDe.coeff(S, 1), cDe.coeff(S, 0)
    assert sp.simplify(cD - (cD0 + cD1 / S)) == 0 and sp.Poly(cD1, lg).degree() == 2
    cD1c = [sp.Poly(cD1, lg).coeff_monomial(lg**n) for n in range(3)]
    F0, F1, F2 = -cot, -lg * cot - th, -lg**2 * cot + K
    AD = 32 * (-3 * K - 8 * q * th) + cD1c[0] * F0 + cD1c[1] * F1 + (cD1c[2] + (32 if MUTATE == "wrong_AD" else 0)) * F2
    d.update(AD=AD, cD1c=cD1c)
    return d


def check_reduce5():
    t0 = time.time()
    dp = sym_pieces()
    lg, S, u, C, D3, Lc, th, q, K = (dp[k] for k in ("lg", "S", "u", "C", "D3", "Lc", "th", "q", "K"))
    w0, w1, w2, w3 = dp["w0"], dp["w1"], dp["w2"], dp["w3"]
    Su = u**2 / (1 + u**2)
    qv = sp.log(2)
    w1v = 2 * sp.pi**2 / 3 + 8 * qv**2
    w2v = 80 * sp.zeta(3) - sp.Rational(8, 3) * sp.pi**2 * qv - sp.Rational(32, 3) * qv**3

    def Dth(f):                                   # d/dth on the generators: ell' = -u, th' = 1, K' = -2 ell, u' = 1 + u^2, S' = 2 u (1 - S), C' = 2 ell + 2 q
        return (sp.diff(f, lg) * (-u) + sp.diff(f, th) + sp.diff(f, K) * (-2 * lg) + sp.diff(f, u) * (1 + u**2)
                + sp.diff(f, S) * 2 * u * (1 - S) + sp.diff(f, C) * (2 * lg + 2 * q))

    # RD0: c_1 = Omega0 - 3 ell - S, (4 ell + 2 c_1)/S = 2 P', P = cot(th)(ell - Omega0)
    c1 = cc(1, [w0, w1, w2, w3, 0], -2 * lg, 1 / (1 - S))
    P = (lg - w0) / u
    r0 = [sp.simplify(c1 - (w0 - 3 * lg - S)), sp.simplify(((4 * lg + 2 * c1) / S - 2 * Dth(P)).subs(S, Su))]
    report("RD0", all(x == 0 for x in r0), f"c_1 = Omega0 - 3 ell - S and (4 ell + 2 c_1)/S = 2 P', P = cot(th)(ell - Omega0): residuals {[str(x) for x in r0]}")
    # numbers: Omega0, Omega1, Omega2 closed, generator values from order4.basics
    mp.mp.dps = 45
    W0, W1, W2 = -4 * mp.ln2, 2 * mp.pi**2 / 3 + 8 * mp.ln2**2, 80 * mp.zeta(3) - mp.mpf(8) / 3 * mp.pi**2 * mp.ln2 - mp.mpf(32) / 3 * mp.ln2**3
    Om = [W0, W1, W2, mp.mpf(0), mp.mpf(0)]
    # RD1: E3 at j = 3 (the same step as for Y_4): int 2 P' Y_3 = -4 int (ell - Omega0) F_2, closed Y_3, L_0 .. L_2
    def e3_lhs(th_, x, tail):
        ell, s2, u_, Cn, D3n = o4.basics(th_, x, tail)
        Y3 = 4 * W1 * ell - 4 * ell**3 + 32 * ell**2 - 16 * W0 * ell + 32 * u_ * Cn - 16 * D3n
        return 2 * (W0 - ell - s2) / s2 * Y3

    def e3_rhs(th_, x, tail):
        ell, s2, u_, Cn, D3n = o4.basics(th_, x, tail)
        l_ = -2 * ell
        L2 = -2 * W1 - 4 * W0 * ell + 16 * ell**2 + 16 * u_ * Cn
        F2 = L2 + (-2 * W0 + 4 * ell) * l_ - l_**2 / 2
        return -e3f() * (ell - W0) * F2

    lhs, rhs = o4.both(o4.guarded(e3_lhs, 30)), o4.both(o4.guarded(e3_rhs, 30))
    d1 = digits(lhs, rhs)
    report("RD1", d1 >= 30, f"E3 at j = 3: int 2 P' Y_3 dth = {nd(lhs, 30)} against -{e3f()} int (ell - Omega0) F_2 dth = {nd(rhs, 30)}: {float(d1):.1f} digits")
    # RD2: the decomposition rest + cC C + cD D3 + cL Lc2 against the numerical integrand of method A (a5_integrand) at four points
    args = (lg, S, u, C, D3, Lc, th, q, w1, w2, w3)
    parts_f = sp.lambdify(args, dp["rest"] + dp["cC"] * C + dp["cD"] * D3 + dp["cL"] * Lc, "mpmath")
    fA = a5_integrand(Om)
    worst = mp.mpf(99)
    for t in ("0.2", "0.7", "1.1", "1.4"):
        th_ = mp.mpf(t)
        with mp.extradps(40):
            ell, s2, u_, Cn, D3n = o4.basics(th_, mp.pi / 2 - th_, False)
            lc2 = o4.Lc2_poly(th_)
            dec = parts_f(ell, s2, u_, Cn, D3n, lc2, th_, mp.ln2, W1, W2, mp.mpf(0))
            num = fA(th_, mp.pi / 2 - th_, False)
        worst = min(worst, digits(num, dec))
    report("RD2", worst >= 30, f"I_5 = rest + cC C + cD D3 + cL Lc2 (sympy) against method A's integrand with Y_4 removed, th = 0.2, 0.7, 1.1, 1.4: worst {float(worst):.1f} digits")
    # RD3: A_D' = cD (symbolic), A_D(pi/2) = 0, and int cD D3 = -2 int A_D C (numerical)
    AD = dp["AD"]
    r3 = sp.simplify((Dth(AD) - dp["cD"]).subs(S, Su))
    cD1c = dp["cD1c"]
    AD_half = 32 * (-3 * sp.pi * q - 4 * sp.pi * q) + cD1c[1] * (-sp.pi / 2) + cD1c[2] * (sp.pi * q) + (32 if MUTATE == "wrong_AD" else 0) * sp.pi * q
    r3b = sp.simplify(AD_half)
    ADf = sp.lambdify((lg, S, u, C, th, q, w1), AD.subs(K, 2 * q * th - C), "mpmath")
    cDf = sp.lambdify((lg, S, u, C, D3, th, q, w1), dp["cD"], "mpmath")

    def lhs3(th_, x, tail):
        if tail:
            th_ = mp.pi / 2 - x
        ell, s2, u_, Cn, D3n = o4.basics(th_, x, tail)
        return cDf(ell, s2, u_, Cn, D3n, th_, mp.ln2, W1) * D3n

    def rhs3(th_, x, tail):
        if tail:
            th_ = mp.pi / 2 - x
        ell, s2, u_, Cn, D3n = o4.basics(th_, x, tail)
        return -2 * ADf(ell, s2, u_, Cn, th_, mp.ln2, W1) * Cn

    a3, b3 = o4.both(o4.guarded(lhs3, 30)), o4.both(o4.guarded(rhs3, 30))
    d3 = digits(a3, b3)
    report("RD3", r3 == 0 and r3b == 0 and d3 >= 30, f"A_D' - cD = {r3}, A_D(pi/2) = {r3b} (so the boundary term [D3 A_D] vanishes; D3(0) = 0); int cD D3 dth = {nd(a3, 30)} "
           f"against -2 int A_D C dth = {nd(b3, 30)}: {float(d3):.1f} digits")
    # RD4: the C terms after the by parts: no 1/u, no th, = u P(ell) + 512 C
    tot = sp.together((dp["cC"] - 2 * AD).subs(K, 2 * q * th - C).subs(S, Su))
    expect = u * (64 * (16 * q**2 + 8 * q + w1) + 128 * (8 * q - 3) * lg - 512 * lg**2) + 512 * C
    r4 = sp.simplify(sp.expand(sp.cancel(tot - expect)))
    report("RD4", r4 == 0, f"cC - 2 A_D (K = 2 q th - C) = u (64 (16 q^2 + 8 q + w1) + 128 (8 q - 3) ell - 512 ell^2) + 512 C: residual {r4}")
    # RD5: int u ell^n C = (2 J_{n+2} + 2 q J_{n+1})/(n+1), n = 0 .. 2, and int C^2 = pi^5/360
    J = o4.moments()
    if MUTATE == "wrong_J4":
        J[4] = J[4] + sp.pi * sp.zeta(3) / 8
    Jn = [mp.mpf(sp.N(j, 60)) for j in J]

    def uC(n):
        def f(th_, x, tail):
            ell, s2, u_, Cn, D3n = o4.basics(th_, x, tail)
            return u_ * ell**n * Cn
        return f

    dd = [digits(o4.both(o4.guarded(uC(n), 30)), (2 * Jn[n + 2] + 2 * mp.ln2 * Jn[n + 1]) / (n + 1)) for n in range(3)]
    c2num = o4.both(o4.guarded(lambda th_, x, tail: o4.basics(th_, x, tail)[3] ** 2, 30))
    c2cl = mp.pi**5 / (180 if MUTATE == "wrong_C2" else 360)
    dc2 = digits(c2num, c2cl)
    report("RD5", min(dd) >= 30 and dc2 >= 30, f"int u ell^n C dth, n = 0, 1, 2, against (2 J_(n+2) + 2 q J_(n+1))/(n + 1): worst {float(min(dd)):.1f} digits; "
           f"int C^2 dth = {nd(c2num, 30)} against {'pi^5/180' if MUTATE == 'wrong_C2' else 'pi^5/360'} (Parseval): {float(dc2):.1f} digits")
    # RD6: the u (Lc2, th) term: the divergent constants cancel, and the regularized integral
    k_th = 128 * q * w1                                  # the coefficient of th u in rest (read from sympy below)
    kth_sym = sp.simplify(sp.Poly(sp.expand(dp["rest"] * S), th).coeff_monomial(th) / (u * S))
    J2s = sp.pi / 2 * (qv**2 + sp.pi**2 / 12)
    r6 = sp.simplify(sp.expand_log((dp["cL"].coeff(u, 1) * J2s * 1 + kth_sym * sp.pi / 2).subs({q: qv, w1: w1v}), force=True))
    regJ = (Jn[2] if MUTATE == "wrong_reg" else Jn[1])

    def regint(th_, x, tail):
        if tail:
            th_ = mp.pi / 2 - x
        ell, s2, u_, Cn, D3n = o4.basics(th_, x, tail)
        lc2 = o4.Lc2_poly(th_)
        return u_ * (-1024 * mp.ln2 * (lc2 - Jn[2]) + 128 * mp.ln2 * W1 * (th_ - mp.pi / 2))

    reg_num = o4.both(o4.guarded(regint, 90))
    reg_cl = -1024 * mp.ln2 * Jn[3] + 128 * mp.ln2 * W1 * regJ
    d6 = digits(reg_num, reg_cl)
    report("RD6", sp.simplify(kth_sym - k_th.subs({q: q, w1: w1})) == 0 and r6 == 0 and d6 >= 25,
           f"coefficient of th u in rest = 128 q w1, (coefficient of u Lc2) J_2 + 128 q w1 pi/2 = {r6} (the divergence at pi/2 cancels); "
           f"int u [-1024 q (Lc2 - J_2) + 128 q w1 (th - pi/2)] dth = {nd(reg_num, 25)} against -1024 q J_3 + 128 q w1 J_1 = {nd(reg_cl, 25)}: {float(d6):.1f} digits (bar 25)")
    # RD7: the moments
    def ellb(b):
        def f(th_, x, tail):
            ell, s2, _, _, _ = o4.basics(th_, x, tail)
            return ell**b / s2
        return f

    d7a = min(digits(o4.both(o4.guarded(ellb(b), 30)), -b * Jn[b - 1]) for b in range(1, 6))
    d7b = min(digits(o4.both(o4.guarded(lambda th_, x, tail, n=n: o4.basics(th_, x, tail)[0] ** n, 30)), Jn[n]) for n in range(5))
    report("RD7", d7a >= 30 and d7b >= 30, f"int ell^b/S dth = -b J_(b-1), b = 1 .. 5: worst {float(d7a):.1f} digits; J_n, n <= 4, from the Beta function (sympy) against quadrature: "
           f"worst {float(d7b):.1f} digits; J_4 = {sp.simplify(J[4])}")
    # RD8: assembly (sympy)
    rest = dp["rest"]
    rest_S = sp.expand(rest * S)
    polyS = sp.Poly(dp["red_u"](rest_S.subs(th, 0)), S)
    A_, B_ = polyS.coeff_monomial(1), polyS.coeff_monomial(S)
    assert polyS.degree() <= 1
    tot_int = 0
    for (b,), cf in sp.Poly(A_, lg).terms():
        assert b >= 1, "A has a constant term: the integral diverges at the throat"
        tot_int += -b * cf * J[b - 1]
    for (b,), cf in sp.Poly(B_, lg).terms():
        tot_int += cf * J[b]
    pc = {2: -512, 1: 1024 * q - 384, 0: 1024 * q**2 + 512 * q + 64 * w1}
    tot_int += sum(p * (2 * J[n + 2] + 2 * q * J[n + 1]) / (n + 1) for n, p in pc.items())
    tot_int += 512 * (sp.pi**5 / (180 if MUTATE == "wrong_C2" else 360))
    tot_int += -1024 * q * J[3] + 128 * q * w1 * (J[2] if MUTATE == "wrong_reg" else J[1])
    G5 = sp.expand_log(sp.expand(tot_int.subs(q, qv).subs(w0, -4 * qv)), force=True).subs({w1: w1v, w2: w2v})
    sol = sp.solve(sp.expand(G5), w3)
    om3 = sp.simplify(sp.expand(sol[0]))
    target = 22 * sp.pi**4 / 45 + 16 * sp.pi**2 * qv**2 / 3 + 32 * qv**4 / 3 - 320 * sp.zeta(3) * qv - 384 * sp.zeta(3)
    diff8 = sp.simplify(sp.expand(sp.expand_log(om3 - target, force=True)))
    slope8 = sp.simplify(sp.diff(sp.expand(G5), w3) + 4 * sp.pi)
    report("RD8", diff8 == 0 and slope8 == 0, f"G_5 = 0 gives Omega3 = {om3}; difference from 22 pi^4/45 + 16 pi^2 ln^2 2/3 + 32 ln^4 2/3 - 320 zeta(3) ln 2 - 384 zeta(3): {diff8}; "
           f"slope of G_5 in Omega3 + 4 pi = {slope8}")
    G5_0 = mp.mpf(sp.N(G5.subs(w3, 0), 60))
    print(f"  derived G_5(Omega3 = 0) = {nd(G5_0, 40)} ({time.time() - t0:.0f} s)")
    return G5_0


def a5_integrand(Om):
    """integrand of G_5 by method A with Y_4 removed by E3; Om = [w0, w1, w2, w3, 0], closed forms of order4.py"""
    q = mp.ln2

    def f(th, x, tail):
        if tail:
            th = mp.pi / 2 - x                 # at the guarded precision: u ~ 1/x multiplies terms in th
        c, lc, l, u, R, s2, Y, L, K, T = o4.closed(th, x, tail, Om)
        Cn = mp.clsin(2, mp.pi - 2 * th)
        sc = mp.sin(th) * mp.cos(th)
        B = -32 * (lc + q) * Cn + 16 * s2 * Cn + 32 * o4.Lc2_poly(th) + 16 * sc * lc * (lc + 2 * q) - 4 * Om[1] * th - 2 * Om[1] * sc
        LL = [L[0], L[1], L[2], -2 * hk(3, Om, l, R, Y, L) + 2 * u * B]
        h4 = s2 * LL[3] + sum(Y[j] * cc(4 - j, Om, l, R) for j in range(4))
        o = list(Om)
        o[4] = 0
        t3 = 2 * sum(Y[j] * cc(5 - j, o, l, R) for j in range(4)) / s2
        F3 = sum(LL[3 - j] * l**j / math.factorial(j) for j in range(4))
        return 8 * lc * h4 / s2 + 2 * sum(LL[3 - j] * l ** (j + 1) / math.factorial(j + 1) for j in range(4)) + t3 - e3f() * (lc - Om[0]) * F3

    return o4.guarded(f, 80)


def check_a5(ev=None, derived=None):
    """A5: G_5 by method A (closed forms of order4.py, E1, E2, E3) at 40 digits: Omega3 from the affine G_5, against the derived closed form and B"""
    t0 = time.time()
    dps = 40
    mp.mp.dps = dps
    cl = closed_omegas()

    def G(Om):
        f = a5_integrand(Om)
        return (mp.quad(lambda t: f(t, mp.pi / 2 - t, False), [0, mp.mpf(1) / 2, 1], method="gauss-legendre")
                + mp.quad(lambda x: f(mp.pi / 2 - x, x, True), [0, mp.pi / 2 - 1]))

    g0 = G([cl[0], cl[1], cl[2], mp.mpf(0), mp.mpf(0)])
    g1 = G([cl[0], cl[1], cl[2], mp.mpf(1), mp.mpf(0)])
    slope, om3 = g1 - g0, -g0 / (g1 - g0)
    mp.mp.dps = 60
    cf = closed_omega3()
    ds, dc = digits(-4 * mp.pi, slope), digits(cf, om3)
    msg = f"G_5(Omega3 = 0) = {nd(g0, 30)}, slope against -4 pi: {float(ds):.1f} digits; Omega3 = {nd(om3, 38)} against the closed form: {float(dc):.1f} digits"
    ok = ds >= 30 and dc >= 30
    if ev is not None:
        db = digits(ev["B4"]["Om3"], om3)
        msg += f"; against B4: {float(db):.1f} digits"
        ok = ok and db >= 30
    if derived is not None:
        dg = digits(derived, g0)
        msg += f"; against the derived G_5(0) = {nd(derived, 30)}: {float(dg):.1f} digits"
        ok = ok and dg >= 30
    report("A5", ok, msg + f" (bar 30; {time.time() - t0:.0f} s)")


def main():
    parse(sys.argv[1:])
    mp.mp.dps = 60
    print(f"order5.py  mutate={MUTATE}  parts={PARTS}  mpmath {mp.__version__}  sympy {sp.__version__}", flush=True)
    t0 = time.time()
    if "grid" in PARTS:
        check_grid()
    if "sym" in PARTS:
        check_sym()
    ev = None
    if "evals" in PARTS:
        ev = check_evals()
        check_values(ev)
    elif LOAD and any(p in PARTS for p in ("values", "pslq", "remainder", "closedform", "a5")):
        mp.mp.dps = 330
        ev = load_ev(LOAD)
        if "values" in PARTS:
            check_values(ev)
    if "closedform" in PARTS:
        check_closedform(ev)
    derived = None
    if "reduce5" in PARTS:
        derived = check_reduce5()
    if "a5" in PARTS:
        check_a5(ev, derived)
    if "pslq" in PARTS:
        check_pslq(ev)
    if "remainder" in PARTS:
        check_remainder(ev)
    nfail = sum(1 for _, ok in RESULTS if not ok)
    print(f"\n{len(RESULTS) - nfail} passed, {nfail} failed of {len(RESULTS)} checks; {time.time() - t0:.0f} s")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()

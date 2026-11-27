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

    python3 order5.py [sym] [evals] [pslq] [remainder] [--mutate=NAME] [--w2=w2_output.txt] [--w2raw=FILE] [--save=FILE] [--load=FILE] [--no-real]
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
sys.argv = _argv

MUTATIONS = ("wrong_h", "wrong_c5", "wrong_Om0", "wrong_Om1", "wrong_Om2", "coarse_B", "short_tail", "pslq_wrong_basis",
             "wrong_Om3_series")
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
    elif LOAD and any(p in PARTS for p in ("values", "pslq", "remainder")):
        mp.mp.dps = 330
        ev = load_ev(LOAD)
        if "values" in PARTS:
            check_values(ev)
    if "pslq" in PARTS:
        check_pslq(ev)
    if "remainder" in PARTS:
        check_remainder(ev)
    nfail = sum(1 for _, ok in RESULTS if not ok)
    print(f"\n{len(RESULTS) - nfail} passed, {nfail} failed of {len(RESULTS)} checks; {time.time() - t0:.0f} s")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()

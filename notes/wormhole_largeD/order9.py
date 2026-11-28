"""Order 9 of the matched expansion (matched.py, order3.py ... order8.py): Omega7 from G9 = 0 by method B, the structure of F at eps^8, and the checks of issue #315.

W^2 a^2 = s + Omega, s = 2(D - 3), eps = 1/s, Omega = Omega0 + Omega1 eps + Omega2 eps^2 + ..., R = r^s, l = ln R = -2 ln cos(theta).
The order equations are order8.py's, for k = 0 .. 9, with the new
  c_9 = Om0 l^8/8! + Om1 l^7/7! + ... + Om7 l + Om8 + l^9/241920 - l^8/40320 + l^8/(40320 R)
(the coefficient of eps^9 of (1 + eps Om) e^{eps l} - (1/2 + eps)(1 - 1/R) e^{-eps l} - e^{-eps l}/(2R); check E0 derives it from the master equation).
G_k is affine in Omega_{k-2} with slope -4 pi: G2 = 0 gives Omega0, ..., G8 = 0 Omega6, **G9 = 0 gives Omega7**.
order8.py is imported, not copied: its table `cc` (patched here to add m = 9), `PanelK`, the cached chain `omega_chain` and the uncached `omega_chain_ref`, `run_B`, the W^2 readers and
`refit_windows`. Omega0 .. Omega5 are planted at their closed forms (#26, #208, #220, #263), Omega6 against #282's computed value (430 digits) and its identified closed form.

The structure: F(eps) = exp(4 ln 2 eps)(1 + Omega0 eps + Omega1 eps^2 + ...) = 1 + f2 eps^2 + ...; f7 is ln 2-free by #282's identification, and
f8 = Omega7 + sum_{j<=6} Omega_j (4 ln 2)^(7-j)/(7-j)! + (4 ln 2)^8/8! (Omega0 .. Omega6 exact: closed forms) is PSLQ's target, in bases of ln 2-free monomials of weight <= 8
that include zeta(5,3) (weight 8 has the four constants pi^8, pi^2 zeta(3)^2, zeta(3) zeta(5), zeta(5,3): Zagier's dimension 4).

    python3 order9.py grid --grid=N,w,V,dps [--kmax=9]      (a development run of one grid; no gate)

Parts (the arguments of a run; the checks of omega7.md and of the raw README of #315 in parentheses):
    sym      E0, the table against the epsilon-expansion of the master equation, k <= 9
    cache    CACHE, the cached chain against the uncached chain on a small grid
    struct   LN2 (order8.py's, f2 .. f6) and LN2b (f7 from the printed closed form of Omega6, ln 2-free, equal to the PSLQ relation of #282) and TR6 (the printed Omega6 equals the one from the relation)
    zv       ZV1 .. ZV4: zeta(5,3) by an asymptotic tail against the stuffle and shuffle relations, zeta(6,2) and a PSLQ of the weight-8 constants
    evals    the grids (--grids=G1,G2,B1,B2,B3 or --cgrids), --save=FILE.json
    values   PL0 .. PL6, PL6c, SL, RES, M7, RNG (--load=FILE.json,...; --ref282=FILE.json for PL6)
    pslq     PQ(basis) and the real PSLQ of f8 in the bases F1 .. F7, F7m (ln 2-free) and X1, X2, X3, XF (with ln 2) (--bases=F7,X1 --no-real --pslq-cap=D --save-rel=FILE)
    plant    writes --save=FILE.json: a planted f8 (a ln 2-free combination of the 13 monomials of weights 3 .. 8) as Omega7 of two "grids", for the structure controls
    closed7  ID7, ZL: the closed form of Omega7 from a PSLQ relation (--rel=FILE) against the grids, and zeros on the ln 2 monomials
    remainder R8, R9, X5 against #195's W^2 (--w2raw=inputs_w2_growthD.txt)
    table    the table of the bases (n, H, budget) for the raw README; --check=README.md compares it with the block between <!-- table --> and <!-- /table --> (TBL)
Mutations: --mutate=NAME (MUTATIONS below). Exit status 1 when a check fails.
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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import order8 as o8  # noqa: E402

o8_cc = o8.cc
o8_closed_omegas = o8.closed_omegas
report, digits, nd = o8.report, o8.digits, o8.nd

MUTATIONS = ("wrong_h", "wrong_c6", "wrong_c7", "wrong_c8", "wrong_c9",
             "wrong_Om0", "wrong_Om1", "wrong_Om2", "wrong_Om3", "wrong_Om4", "wrong_Om5", "wrong_Om6",
             "coarse_B", "short_tail", "stale_cache", "pslq_wrong_basis", "pslq_wrong_target", "wrong_Om7_series", "wrong_closed5", "wrong_closed6",
             "ln2_leak", "f8_rel", "ln2_struct", "ln2_zl", "zeta53_wrong")
MUTATE = None
W2FILE = "w2_output.txt"
W2RAW = None
SAVE = None
LOAD = None
REF282 = None
GRIDARG = None
GRIDSEL = None
CGRIDS_ON = False
KMAX = 9
PARTS = []
REAL = True
BASES_SEL = None
PSLQ_CAP = None
REL = None
SAVE_REL = None
CHECK_README = None
ZV_DPS = 560


def parse(argv):
    global MUTATE, W2FILE, W2RAW, SAVE, LOAD, REF282, GRIDARG, GRIDSEL, CGRIDS_ON, KMAX, PARTS, REAL, BASES_SEL, PSLQ_CAP, REL, SAVE_REL, CHECK_README, ZV_DPS
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
        elif a.startswith("--ref282="):
            REF282 = a.split("=", 1)[1]
        elif a.startswith("--grid="):
            GRIDARG = a.split("=", 1)[1]
        elif a.startswith("--grids="):
            GRIDSEL = a.split("=", 1)[1].split(",")
        elif a.startswith("--bases="):
            BASES_SEL = a.split("=", 1)[1].split(",")
        elif a.startswith("--pslq-cap="):
            PSLQ_CAP = int(a.split("=", 1)[1])
        elif a.startswith("--check="):
            CHECK_README = a.split("=", 1)[1]
        elif a.startswith("--zv-dps="):
            ZV_DPS = int(a.split("=", 1)[1])
        elif a == "--cgrids":
            CGRIDS_ON = True
        elif a == "--no-real":
            REAL = False
        elif a.startswith("--kmax="):
            KMAX = int(a.split("=", 1)[1])
        else:
            parts.append(a)
    PARTS = parts or ["sym", "struct", "evals", "values", "pslq", "remainder"]
    o8.MUTATE = MUTATE
    o8.W2FILE = W2FILE


# ---------------------------------------------------------------------------------------------------------------
# the table: c_9 added to order8.py's c_1 .. c_8 (the patch is one assignment: order8's PanelK and chains look `cc` up in their module)
# ---------------------------------------------------------------------------------------------------------------
def cc9(m, Om, l, R):
    if m == 9:           # new at this order: Om_j l^(8-j)/(8-j)! (j = 0 .. 8), l^9/241920 - l^8/40320 + l^8/(40320 R)
        return (Om[0] * l**8 / 40320 + Om[1] * l**7 / 5040 + Om[2] * l**6 / 720 + Om[3] * l**5 / 120 + Om[4] * l**4 / 24 + Om[5] * l**3 / 6
                + Om[6] * l**2 / 2 + Om[7] * l + Om[8] + l**9 / 241920 - l**8 / 40320 + (-1 if MUTATE == "wrong_c9" else 1) * l**8 / (40320 * R))
    return o8_cc(m, Om, l, R)


# ---------------------------------------------------------------------------------------------------------------
# closed forms: Omega0 .. Omega5 (#26, #208, #220, #263) and Omega6 identified in #282 (from f7 by exact algebra)
# ---------------------------------------------------------------------------------------------------------------
def f7_closed():
    """f7 of #282's PSLQ relation (F6 at 390 digits; ln 2-free)"""
    pi, z3, z5, z7 = mp.pi, mp.zeta(3), mp.zeta(5), mp.zeta(7)
    p2 = pi**2
    return (mp.mpf(7656704) / 2835 * p2**3 + mp.mpf(328672) / 45 * p2**2 * z3 - mp.mpf(74752) / 5 * p2**2 + 30720 * p2 * z3 - 2176 * p2 * z5
            + 36864 * z3**2 + 565248 * z3 - 307200 * z5 + 974336 * z7)


def om6_closed():
    """Omega6 = f7 - sum_{j<=5} Omega_j (4 ln 2)^(6-j)/(6-j)! - (4 ln 2)^7/7!"""
    a = [mp.mpf(1)] + o8_closed_omegas()
    q4 = 4 * mp.ln2
    return f7_closed() - sum(a[j] * q4 ** (7 - j) / mp.factorial(7 - j) for j in range(7)) + (mp.mpf(1) / 1000 if MUTATE == "wrong_closed6" else 0)


def closed_omegas9():
    return o8_closed_omegas() + [om6_closed()]


def sym_omegas():
    """Omega0 .. Omega5 as sympy expressions of q = ln 2, p2 = pi^2, z3, z5 (the closed forms), and the symbols"""
    q, z3, z5, z7, p2, z53 = sp.symbols("q z3 z5 z7 p2 z53")
    R = sp.Rational
    Om = [-4 * q, 2 * p2 / 3 + 8 * q**2,
          80 * z3 - R(8, 3) * p2 * q - R(32, 3) * q**3,
          22 * p2**2 / 45 + R(16, 3) * p2 * q**2 + R(32, 3) * q**4 - 320 * z3 * q - 384 * z3,
          -10432 * z5 - 160 * p2 * z3 + 640 * q**2 * z3 - R(64, 9) * p2 * q**3 - R(88, 45) * p2**2 * q - R(128, 15) * q**5 - R(2432, 45) * p2**2 + 1536 * q * z3 + 3072 * z3,
          R(256, 45) * q**6 + R(64, 9) * p2 * q**4 - R(2560, 3) * z3 * q**3 + R(176, 45) * p2**2 * q**2 - 3072 * z3 * q**2 + 640 * p2 * z3 * q + R(9728, 45) * p2**2 * q
          - 12288 * z3 * q + 41728 * z5 * q - R(38812, 315) * p2**3 + 3200 * z3**2 - R(256, 3) * p2 * z3 + 74240 * z5 + R(39424, 45) * p2**2 - 36864 * z3]
    if MUTATE == "ln2_struct":
        Om[3] = Om[3] + q * z3
    return Om, (q, z3, z5, z7, p2, z53)


def f7_sym():
    """f7 of the relation, sympy"""
    Om, (q, z3, z5, z7, p2, z53) = sym_omegas()
    R = sp.Rational
    return (R(7656704, 2835) * p2**3 + R(328672, 45) * p2**2 * z3 - R(74752, 5) * p2**2 + 30720 * p2 * z3 - 2176 * p2 * z5 + 36864 * z3**2 + 565248 * z3
            - 307200 * z5 + 974336 * z7)


def om6_printed():
    """Omega6 as printed in omega6.md (24 terms), sympy"""
    Om, (q, z3, z5, z7, p2, z53) = sym_omegas()
    R = sp.Rational
    return (R(155248, 315) * p2**3 * q + R(7656704, 2835) * p2**3 - R(704, 135) * p2**2 * q**3 - R(19456, 45) * p2**2 * q**2 - R(157696, 45) * p2**2 * q
            + R(328672, 45) * p2**2 * z3 - R(74752, 5) * p2**2 - R(256, 45) * p2 * q**5 - 1280 * p2 * z3 * q**2 + R(1024, 3) * p2 * z3 * q + 30720 * p2 * z3 - 2176 * p2 * z5
            - R(1024, 315) * q**7 + R(2560, 3) * z3 * q**4 + 4096 * z3 * q**3 + 24576 * z3 * q**2 - 83456 * z5 * q**2 - 12800 * z3**2 * q + 147456 * z3 * q - 296960 * z5 * q
            + 36864 * z3**2 + 565248 * z3 - 307200 * z5 + 974336 * z7)


def fsym(n, a, q):
    return sp.expand(sum(a[j] * (4 * q) ** (n - j) / sp.factorial(n - j) for j in range(n + 1)))


# ---------------------------------------------------------------------------------------------------------------
# E0 through k = 9, CACHE, LN2, LN2b, TR6
# ---------------------------------------------------------------------------------------------------------------
def check_sym9():
    t0 = time.time()
    R, eps, l = sp.symbols("R epsilon l", positive=True)
    kk = 9
    Om = sp.symbols("Omega0:9")
    Ys, Ls = sp.symbols("Y0:10"), sp.symbols("L0:10")
    r = sum((eps * l) ** j / sp.factorial(j) for j in range(kk + 1))          # exact through eps^9, which is all that is read
    rinv = sum((-eps * l) ** j / sp.factorial(j) for j in range(kk + 1))
    Y = sum(Ys[k] * eps**k for k in range(kk + 1))
    L = sum(Ls[k] * eps**k for k in range(kk + 1))
    O = sum(Om[k] * eps**k for k in range(kk))
    half = sp.Rational(1, 2)
    rhs = ((1 + eps * O) * r * Y + L / (2 * R) + (half + eps) * L * (R - 1) / R - (half + eps) * Y * (R - 1) / R * rinv - Y / (2 * R) * rinv)
    ser = sp.expand(rhs)
    yser = sp.expand(eps * L * r)
    res = []
    for k in range(kk + 1):
        res.append(sp.simplify(ser.coeff(eps, k) - Ls[k] / 2 - o8.hk(k, Om, l, R, Ys, Ls)))
        if k >= 1:
            res.append(sp.simplify(yser.coeff(eps, k) - o8.Fsum(k, Ls, l)))
    report("E0", all(x == 0 for x in res),
           f"h_k (k = 0 .. 9) and dY_k/dl (k = 1 .. 9) from the epsilon-expansion of the master equation against the table: "
           f"residuals {[str(x) for x in res if x != 0] or '0'} ({time.time() - t0:.1f} s)")


def check_cache9():
    N, w, V, dps = 24, (1, 2), 12, 60
    mp.mp.dps = 40
    a = o8.run_B(N, w, V, dps, kmax=9)
    b = o8.run_B(N, w, V, dps, kmax=9, cached=False)
    mp.mp.dps = dps
    worst = max(abs(a[f"Om{i}"] - b[f"Om{i}"]) / abs(b[f"Om{i}"]) for i in range(8))
    dig = [float(digits(b[f"Om{i}"], a[f"Om{i}"])) for i in range(8)]
    report("CACHE", worst <= mp.mpf(10) ** -(dps - 5),
           f"cached chain against the uncached chain on (N {N}, w 1/2, V {V}, {dps} digits): Omega0 .. Omega7 equal to "
           + ", ".join(f"{d:.1f}" for d in dig) + f" digits, worst relative difference {nd(worst, 3)} (bar 1e-{dps - 5}); times {a['secs']:.0f} s and {b['secs']:.0f} s")


def check_struct9():
    o8.check_struct()                                           # LN2: f2 .. f6 from the closed forms of Omega0 .. Omega5
    t0 = time.time()
    Om, (q, z3, z5, z7, p2, z53) = sym_omegas()
    f7s = f7_sym()
    a6 = [sp.Integer(1)] + Om
    om6_rel = sp.expand(f7s - sum(a6[j] * (4 * q) ** (7 - j) / sp.factorial(7 - j) for j in range(7)))
    om6_prn = om6_printed()
    tr = sp.expand(om6_rel - om6_prn)
    report("TR6", tr == 0, f"Omega6 printed in omega6.md (24 terms) minus Omega6 from the f7 relation by exact algebra: {tr if tr != 0 else '0'} ({time.time() - t0:.1f} s)")
    a = [sp.Integer(1)] + Om + [om6_prn]
    f7 = fsym(7, a, q)
    dq = sp.expand(sp.diff(f7, q))
    eq = sp.expand(f7 - f7s)
    report("LN2b", dq == 0 and eq == 0,
           f"f7 = [eps^7] exp(4 ln 2 eps)(1 + Omega0 eps + ... + Omega6 eps^7) with the printed Omega6: d f7/d ln 2 = {dq if dq != 0 else '0'} and f7 minus the relation's f7 = "
           f"{eq if eq != 0 else '0'} (sympy)")


# ---------------------------------------------------------------------------------------------------------------
# the evaluations and their checks
# ---------------------------------------------------------------------------------------------------------------
GRIDS = {"G1": (405, (1, 4), 57, 625), "G2": (428, (3, 10), 57, 625)}                       # w as a fraction (num, den), evaluated at the working precision; declared in the raw README
REPORT_GRIDS = {"B1": (60, (2, 5), 26, 100), "B2": (100, (1, 4), 32, 170), "B3": (170, (1, 4), 40, 285)}
CGRIDS = {"G1": (90, (1, 4), 32, 160), "G2": (85, (3, 10), 31, 150)}


def check_grid9():
    N, w, V, dps = GRIDARG.split(",")
    mp.mp.dps = 60
    e = o8.run_B(int(N), tuple(int(t) for t in w.split("/")) if "/" in w else w, int(V), int(dps), kmax=KMAX, log=True)
    mp.mp.dps = int(dps)
    cl = closed_omegas9()
    print(f"  {e['grid']}  {e['secs']:.0f} s (panel setup {e['init_secs']:.0f} s)")
    for i in range(KMAX - 1):
        print(f"    Omega{i} = {nd(e[f'Om{i}'], int(dps) - 20)}")
    for i in range(min(7, KMAX - 1)):
        print(f"    Omega{i} against its closed form: {float(digits(cl[i], e[f'Om{i}'])):.1f} digits")
    for k in range(2, KMAX + 1):
        print(f"    slope G{k} against -4 pi: {float(digits(-4 * mp.pi, e[f'slope{k}'])):.1f} digits")
    for k, r in sorted((int(n[5:]), e[n]) for n in e if n.startswith("resid")):
        print(f"    |G{k}| at the root / |slope| = {nd(r, 3)}")
    return e


def check_evals9():
    mp.mp.dps = 60
    ev = {}
    main_grids = CGRIDS if CGRIDS_ON else GRIDS
    todo = list(main_grids.items()) + ([] if CGRIDS_ON else list(REPORT_GRIDS.items()))
    if GRIDSEL:
        todo = [(n, g) for n, g in todo if n in GRIDSEL]
    for name, (N, w, V, dps) in todo:
        ev[name] = o8.run_B(N, w, V, dps, kmax=9, log=True)
        e = ev[name]
        print(f"  {name}: {e['grid']}  {e['secs']:.0f} s (panel setup {e['init_secs']:.0f} s)\n    Omega0 = {nd(e['Om0'], 60)}\n    Omega1 = {nd(e['Om1'], 60)}"
              f"\n    Omega2 = {nd(e['Om2'], 60)}\n    Omega3 = {nd(e['Om3'], 110)}\n    Omega4 = {nd(e['Om4'], 200)}\n    Omega5 = {nd(e['Om5'], dps - 40)}"
              f"\n    Omega6 = {nd(e['Om6'], dps - 40)}\n    Omega7 = {nd(e['Om7'], dps - 40)}"
              "\n    slopes G2 .. G9 = " + " ".join(nd(e[f"slope{k}"], 20) for k in range(2, 10))
              + "\n    |G_k| at the root / |slope|, k = 2 .. 9 = " + " ".join(nd(e[f"resid{k}"], 3) for k in range(2, 10)), flush=True)
    if SAVE:
        o8.save_ev(ev, SAVE)
    return ev


OM7_REF = mp.mpf("2.55e7")      # #282's two-point estimate (25545476), the reference of PLW for wrong_Om6; only a change of 1e-3 is read
RNG7 = (mp.mpf("2.4e7"), mp.mpf("2.7e7"))


def load_ref282():
    """Omega6 of #282's G1 and G2 (430.8 digits between them), from the saved JSON"""
    out = {}
    with mp.workdps(700):
        for path in REF282.split(","):
            out.update({n: mp.mpf(e["Om6"]) for n, e in json.load(open(path)).items() if "Om6" in e})
    return out


def check_values9(ev):
    mp.mp.dps = max(int(e["grid"].split("dps=")[1].split()[0]) for e in ev.values()) + 30
    cl = closed_omegas9()
    grids = ["G1", "G2"]
    for i in range(6):
        rows = {n: digits(cl[i], ev[n][f"Om{i}"]) for n in ev}
        ok = all(rows[n] >= 100 for n in grids)
        report(f"PL{i}", ok, f"planted: Omega{i} against its closed form, digits: " + ", ".join(f"{n} {float(d):.1f}" for n, d in rows.items()) + " (bar 100 on G1, G2; the others reported)")
    if REF282:
        ref = load_ref282()
        rows = {n: digits(ref["G1"], ev[n]["Om6"]) for n in ev}
        report("PL6", all(rows[n] >= 100 for n in grids),
               "planted: Omega6 against #282's computed value (its G1; its G2 agrees to " + f"{float(digits(ref['G1'], ref['G2'])):.1f} digits), digits: "
               + ", ".join(f"{n} {float(d):.1f}" for n, d in rows.items()) + " (bar 100 on G1, G2; the others reported)")
    rows = {n: digits(cl[6], ev[n]["Om6"]) for n in ev}
    report("PL6c", all(rows[n] >= 100 for n in grids), "planted: Omega6 against its identified closed form (#282), digits: "
           + ", ".join(f"{n} {float(d):.1f}" for n, d in rows.items()) + " (bar 100 on G1, G2; the others reported)")
    if MUTATE in o8.WRONG_OM:
        i = int(MUTATE[-1]) + 1
        ref_i = cl[i] if i < 7 else OM7_REF
        moved = all(abs(ev[n][f"Om{i}"] - ref_i) > mp.mpf(1) / 1000 for n in grids)
        report("PLW", moved, f"control {MUTATE} moves Omega{i} of " + ", ".join(f"{n} by {nd(ev[n][f'Om{i}'] - ref_i, 5)}" for n in grids))
    rows, ok = [], True
    for n in ev:
        ds = [digits(-4 * mp.pi, ev[n][f"slope{k}"]) for k in range(2, 10)]
        rows.append(f"{n} " + "/".join(f"{float(d):.1f}" for d in ds))
        ok = ok and (n not in grids or min(ds) >= 100)
    report("SL", ok, "slope of G_k in Omega_(k-2) against -4 pi, k = 2 .. 9, digits: " + ", ".join(rows) + " (bar 100 on G1, G2)")
    rs = {n: max(ev[n][f"resid{k}"] for k in range(2, 10)) for n in ev}
    report("RES", all(rs[n] <= mp.mpf(10) ** -100 for n in grids),
           "max over k = 2 .. 9 of |G_k at the chain's root| / |slope|: " + ", ".join(f"{n} {nd(r, 3)}" for n, r in rs.items()) + " (bar 1e-100 on G1, G2)")
    a, b = grids
    d7 = digits(ev[a]["Om7"], ev[b]["Om7"])
    report("M7", d7 >= 100, f"Omega7 of {a} against {b} (two grids, two precisions): {float(d7):.1f} digits (bar 100); |diff| = {nd(abs(ev[a]['Om7'] - ev[b]['Om7']), 3)}")
    report("RNG", all(RNG7[0] <= ev[n]["Om7"] <= RNG7[1] for n in grids), f"Omega7 between {nd(RNG7[0], 2)} and {nd(RNG7[1], 2)}: " + ", ".join(f"{n} {nd(ev[n]['Om7'], 12)}" for n in grids))
    others = [n for n in ev if n not in grids]
    print("  reported: digits between evaluations, Omega6 / Omega7: " + ", ".join(
        f"{x}-{y} " + "/".join(f"{float(digits(ev[x][f'Om{i}'], ev[y][f'Om{i}'])):.1f}" for i in (6, 7)) for x, y in [("G1", "G2")] + [(o, "G1") for o in others]))
    return d7


# ---------------------------------------------------------------------------------------------------------------
# zeta(5,3): the double zeta value, by an exact head and an asymptotic tail
# ---------------------------------------------------------------------------------------------------------------
_MZV = {}
_BD = {}


def _bd(j):
    """B_2j/(2j)! at the working precision"""
    key = (j, mp.mp.prec)
    if key not in _BD:
        _BD[key] = mp.bernoulli(2 * j) / mp.factorial(2 * j)
    return _BD[key]


def hz(s, x, eps):
    """zeta(s, x) for an integer s >= 2 and a large integer x (>= 1000): x^(1-s)/(s-1) + x^-s/2 + sum_j B_2j/(2j)! (s)_(2j-1) x^(1-s-2j), summed until a term is below
    eps relative to the sum (mpmath's zeta(s, a) is good to about 290 digits at a = 2001, not to its working precision). The series is asymptotic and its terms fall
    until 2j is about 2 pi x, far beyond what is used; the remainder is below the first omitted term"""
    xm = mp.mpf(x)
    xi2 = 1 / (xm * xm)
    p = xm ** (-s)
    val = p * xm / (s - 1) + p / 2
    rf = mp.mpf(s)
    pj = p / xm
    for j in range(1, 10000):
        term = _bd(j) * rf * pj
        val += term
        if abs(term) < eps * abs(val):
            return val
        rf *= (s + 2 * j - 1) * (s + 2 * j)
        pj *= xi2
    raise ArithmeticError(f"hz({s}, {x}) did not converge")


def mzv2(a, b, dps=None):
    """zeta(a, b) = sum_{n > m >= 1} n^-a m^-b (a >= 2, b >= 2, integers) to `dps` digits: the head n <= N0 exactly, the tail
    sum_{n > N0} n^-a (zeta(b) - zeta(b, n)) with zeta(b, n) = n^(1-b)/(b-1) + n^-b/2 + sum_j B_2j/(2j)! (b)_(2j-1) n^(1-b-2j) (the remainder of an Euler-Maclaurin
    series of positive x is below its first omitted term), summed over n by the Hurwitz zeta values zeta(a + b + 2j - 1, N0 + 1) (`hz`)"""
    dps = dps or mp.mp.dps
    key = (a, b)
    if key in _MZV and _MZV[key][0] >= dps:
        return +_MZV[key][1]
    with mp.workdps(dps + 30):
        N0 = 2000
        eps = mp.mpf(10) ** -(dps + 25)
        H, S = mp.mpf(0), mp.mpf(0)
        for n in range(2, N0 + 1):
            H += mp.mpf(n - 1) ** (-b)
            S += H * mp.mpf(n) ** (-a)
        tail = mp.zeta(b) * hz(a, N0 + 1, eps) - (hz(a + b - 1, N0 + 1, eps) / (b - 1) + hz(a + b, N0 + 1, eps) / 2)
        last = None
        for j in range(1, 4 * dps):
            term = _bd(j) * mp.rf(b, 2 * j - 1) * hz(a + b + 2 * j - 1, N0 + 1, eps)
            if last is not None and abs(term) > abs(last):
                raise ArithmeticError(f"the asymptotic series of zeta({a},{b}) turned up at j = {j}: N0 = {N0} is too small for {dps} digits")
            tail -= term
            last = term
            if abs(term) < eps:
                break
        else:
            raise ArithmeticError(f"the asymptotic series of zeta({a},{b}) did not reach {dps} digits")
        val = S + tail
    _MZV[key] = (dps, val)
    return +val


def zeta53(dps=None):
    v = mzv2(5, 3, dps)
    if MUTATE == "zeta53_wrong":
        v = v * (1 + mp.mpf(10) ** -30)
    return v


def check_zv():
    """ZV1 .. ZV4"""
    t0 = time.time()
    dps = ZV_DPS
    mp.mp.dps = dps + 10
    z = mp.zeta
    z53, z35, z62, z26 = zeta53(dps + 10), mzv2(3, 5, dps + 10), mzv2(6, 2, dps + 10), mzv2(2, 6, dps + 10)
    z44 = (z(4) ** 2 - z(8)) / 2
    z71 = mp.mpf(7) / 2 * z(8) - (z(2) * z(6) + z(3) * z(5) + z(4) * z(4) + z(5) * z(3) + z(6) * z(2)) / 2            # Euler's formula for zeta(m, 1), m = 7
    bar = mp.mpf(10) ** -(dps - 5)
    r1 = abs(z53 + z35 - (z(3) * z(5) - z(8)))
    report("ZV1", r1 <= bar, f"zeta(5,3) + zeta(3,5) = zeta(3) zeta(5) - zeta(8) (stuffle; zeta(3,5) by the same routine with b = 5): residual {nd(r1, 3)} (bar {nd(bar, 2)}); zeta(5,3) = {nd(z53, 40)}")
    r2 = abs(6 * z53 + 15 * z62 + 30 * z71 + 3 * z44 - z(8))
    report("ZV2", r2 <= bar, f"zeta(3) zeta(5) as a shuffle product minus the stuffle: 6 zeta(5,3) + 15 zeta(6,2) + 30 zeta(7,1) + 3 zeta(4,4) - zeta(8) = 0 "
           f"(zeta(6,2) by the routine, zeta(7,1) Euler's formula, zeta(4,4) = (zeta(4)^2 - zeta(8))/2): residual {nd(r2, 3)} (bar {nd(bar, 2)})")
    r2b = abs(z26 + z62 - (z(2) * z(6) - z(8)))
    report("ZV2b", r2b <= bar, f"zeta(2,6) + zeta(6,2) = zeta(2) zeta(6) - zeta(8): residual {nd(r2b, 3)} (bar {nd(bar, 2)})")
    # ZV3: a low-precision independent value (Richardson and Shanks on the plain partial sums, mpmath's nsum)
    with mp.workdps(30):
        ref = mp.nsum(lambda n: (mp.zeta(3) - mp.zeta(3, n)) / n**5, [2, mp.inf])
        d3 = digits(ref, z53)
    report("ZV3", d3 >= 20, f"zeta(5,3) against mp.nsum (30 digits, Richardson/Shanks on sum n^-5 (zeta(3) - zeta(3, n))): {float(d3):.1f} digits (bar 20)")
    # ZV4: the weight-8 constants are independent at this height, and zeta(6,2) is in their span (the relation of ZV2)
    d = dps - 10
    with mp.workdps(d):
        pi, z3, z5 = mp.pi, mp.zeta(3), mp.zeta(5)
        c8 = [pi**8, pi**2 * z3**2, z3 * z5, +z53]
        tau = mp.mpf(10) ** -(d - 20)
        r_ind = mp.pslq(c8, tol=tau, maxcoeff=10**12, maxsteps=10**6)
        r_62 = mp.pslq([+z62] + c8, tol=tau, maxcoeff=10**12, maxsteps=10**6)
    ok4 = r_ind is None and r_62 is not None and r_62[0] != 0 and r_62[2] == 0 and r_62[1] != 0
    report("ZV4", ok4, f"PSLQ at {d} digits, H = 1e12: {{pi^8, pi^2 zeta(3)^2, zeta(3) zeta(5), zeta(5,3)}} has no relation ({'none' if r_ind is None else list(r_ind)}); "
           f"{{zeta(6,2), those four}} has {'none' if r_62 is None else list(r_62)}")
    print(f"  ZV: {time.time() - t0:.0f} s at {dps} digits")


# ---------------------------------------------------------------------------------------------------------------
# PSLQ of f8
# ---------------------------------------------------------------------------------------------------------------
PSLQ_DIGITS = 560            # cap on the working digits; the digits the two grids share, less 32, if fewer
PSLQ_SLACK = 30              # tolerance 10^-(d - slack); the budget plus the slack must stay below d
VERIFY_EXTRA = 20            # a relation must hold at this many more digits than were used to find it
MAXSTEPS = 4 * 10**6


def free9():
    """the ln 2-free monomials of weight <= 8 (pi^2: 2, zeta(3): 3, zeta(5): 5, zeta(7): 7, zeta(5,3): 8, products add): (name, weight, value); the 15 span the Euler sums
    of weight <= 8 (Zagier's dimensions 1, 0, 1, 1, 1, 2, 2, 3, 4 for weights 0 .. 8; independence is conjectural)"""
    pi, z3, z5, z7 = mp.pi, mp.zeta(3), mp.zeta(5), mp.zeta(7)
    p2 = pi**2
    return [("pi^8", 8, p2**4), ("pi^2 zeta(3)^2", 8, p2 * z3**2), ("zeta(3) zeta(5)", 8, z3 * z5), ("zeta(5,3)", 8, zeta53()),
            ("pi^4 zeta(3)", 7, p2**2 * z3), ("pi^2 zeta(5)", 7, p2 * z5), ("zeta(7)", 7, z7), ("pi^6", 6, p2**3), ("zeta(3)^2", 6, z3**2),
            ("pi^2 zeta(3)", 5, p2 * z3), ("zeta(5)", 5, z5), ("pi^4", 4, p2**2), ("zeta(3)", 3, z3), ("pi^2", 2, p2), ("1", 0, mp.mpf(1))]


# basis -> (lowest weight of the ln 2-free monomials, highest power of ln 2, maxcoeff H, planted seed, drop zeta(5,3)); ln 2^j times the free monomials of weight <= 8 - j
PSLQ_BASES = {"F1": (8, 0, 10**16, 315001, False), "F2": (7, 0, 10**16, 315002, False), "F3": (6, 0, 10**16, 315003, False), "F4": (5, 0, 10**16, 315004, False),
              "F5": (4, 0, 10**16, 315005, False), "F6": (3, 0, 10**16, 315006, False), "F7": (0, 0, 10**16, 315007, False), "F7m": (0, 0, 10**16, 315011, True),
              "X1": (0, 1, 10**15, 315008, False), "X2": (0, 2, 10**14, 315009, False), "X3": (0, 3, 10**12, 315010, False), "XF": (0, 8, 10**9, 315012, False)}


def constants(name):
    lo, jmax, _, _, drop = PSLQ_BASES[name]
    q = mp.ln2
    out = []
    fr = free9()
    for j in range(jmax + 1):
        for nm, wt, v in fr:
            if wt >= lo and wt <= 8 - j and not (drop and nm == "zeta(5,3)"):
                pw = "ln 2" if j == 1 else f"ln^{j} 2"
                out.append((nm if j == 0 else f"{pw} {nm}" if nm != "1" else pw, q**j * v))
    return out


def basis_size(name):
    lo, jmax, _, _, drop = PSLQ_BASES[name]
    names = [("pi^8", 8), ("pi^2 zeta(3)^2", 8), ("zeta(3) zeta(5)", 8), ("zeta(5,3)", 8), ("pi^4 zeta(3)", 7), ("pi^2 zeta(5)", 7), ("zeta(7)", 7), ("pi^6", 6), ("zeta(3)^2", 6),
             ("pi^2 zeta(3)", 5), ("zeta(5)", 5), ("pi^4", 4), ("zeta(3)", 3), ("pi^2", 2), ("1", 0)]
    return sum(1 for j in range(jmax + 1) for nm, wt in names if wt >= lo and wt <= 8 - j and not (drop and nm == "zeta(5,3)"))


def budget_of(n, H):
    """(n + 1) log10(2H + 1) + 1: the digits the search space uses up"""
    return (n + 1) * math.log10(2 * H + 1) + 1


def table_text(d=None):
    rows = ["| basis | constants besides f8 | n | H | budget | budget + slack |", "|---|---|---|---|---|---|"]
    desc = {"F1": "the 4 monomials of weight 8", "F2": "F1 and the 3 of weight 7", "F3": "F2 and the 2 of weight 6", "F4": "F3 and the 2 of weight 5", "F5": "F4 and pi^4",
            "F6": "F5 and zeta(3)", "F7": "F6 and pi^2, 1 (the ln 2-free basis, weight <= 8)", "F7m": "F7 without zeta(5,3)",
            "X1": "F7 and ln 2 x the 11 monomials of weight <= 7", "X2": "X1 and ln^2 2 x the 8 of weight <= 6", "X3": "X2 and ln^3 2 x the 6 of weight <= 5",
            "XF": "X3 and ln^4 2 x 4, ln^5 2 x 3, ln^6 2 x 2, ln^7 2, ln^8 2 (the full basis, weight <= 8)"}
    for name, (lo, jmax, H, seed, drop) in PSLQ_BASES.items():
        n = basis_size(name)
        b = budget_of(n, H)
        rows.append(f"| {name} | {desc[name]} | {n} | 1e{round(math.log10(H))} | {b:.1f} | {b + PSLQ_SLACK:.1f} |")
    return "\n".join(rows)


def check_table():
    t = table_text()
    print(t)
    if CHECK_README:
        txt = open(CHECK_README).read()
        m = re.search(r"<!-- table -->\n(.*?)\n<!-- /table -->", txt, re.S)
        ok = bool(m) and m.group(1).strip() == t.strip()
        report("TBL", ok, f"the table of bases in {CHECK_README} equals the code's PSLQ_BASES (H, n, budgets): {ok}")


def planted(n, seed):
    """the planted relation: denominator and n integer coefficients, drawn once from a fixed seed (declared in the raw README)"""
    r = random.Random(seed)
    den = r.choice([3, 5, 7, 9, 11, 13, 15, 21, 45])
    return den, [r.choice([-1, 1]) * r.randint(100, 9999) for _ in range(n)]


def f8_of(om7):
    """f8 from Omega7 and the exact Omega0 .. Omega6 (closed forms)"""
    return o8.f_coeff(8, [mp.mpf(1)] + closed_omegas9() + [om7])


def f8_target(ev, g):
    """f8 of a grid; `ln2_leak` adds 1000 ln 2 zeta(7) (a ln 2 term in F at eps^8)"""
    f = f8_of(ev[g]["Om7"])
    if MUTATE == "ln2_leak":
        f = f + 1000 * mp.ln2 * mp.zeta(7)
    return f


def check_pslq9(ev):
    grids = ["G1", "G2"]
    mp.mp.dps = max(int(e["grid"].split("dps=")[1].split()[0]) for e in ev.values()) + 30
    agree = digits(ev[grids[0]]["Om7"], ev[grids[1]]["Om7"])
    cap = PSLQ_CAP or PSLQ_DIGITS
    d = int(max(30, min(cap, math.floor(agree) - 32)))
    pert = mp.sqrt(2) / 1000
    tau_digits = d - PSLQ_SLACK
    print(f"  PSLQ at {d} digits (the two grids share {float(agree):.1f}), tolerance 1e-{tau_digits}, verification at {d + VERIFY_EXTRA} digits, maxsteps {MAXSTEPS}", flush=True)
    mp.mp.dps = d + VERIFY_EXTRA + 40
    t0 = time.time()
    zeta53(d + VERIFY_EXTRA + 40)                                  # cached for every basis below
    print(f"  zeta(5,3) at {d + VERIFY_EXTRA + 40} digits: {time.time() - t0:.0f} s", flush=True)
    saved = {}
    for name, (_, _, H, seed, _) in PSLQ_BASES.items():
        if BASES_SEL and name not in BASES_SEL:
            continue
        with mp.workdps(d):
            tau = mp.mpf(10) ** (-tau_digits)
            names, cs = zip(*constants(name))
            cs_true = list(cs)
            cs = list(cs)
            if MUTATE == "pslq_wrong_basis":
                cs[0] = cs[0] + pert
            n = len(cs)
            if CGRIDS_ON or PSLQ_CAP:                  # a control run: d is small, so H is lowered until the budget plus the slack is below d - 5 (at least 10^4)
                H = 10 ** max(4, min(int(math.log10(H)), int((d - PSLQ_SLACK - 5) // (n + 1)) - 1))
            budget = budget_of(n, H)
            f8 = +f8_target(ev, grids[1])
            den, c = planted(n, seed)
            xp = sum(ci * b for ci, b in zip(c, cs_true)) / den
            if MUTATE == "pslq_wrong_target":
                xp = xp * (1 + mp.mpf(10) ** -20)
            rp, tp = o8.pslq_run([xp] + cs, tau, H)
            want = [den] + [-ci for ci in c]
            okp = rp is not None and (list(rp) == want or list(rp) == [-v for v in want])
            csp = list(cs)
            csp[1 if n > 1 else 0] = csp[1 if n > 1 else 0] + pert
            rpp, _ = o8.pslq_run([xp] + csp, tau, H)
            report("P" + name, okp and rpp is None,
                   f"{name} = {{f8, {', '.join(names)}}}, n = {n + 1}, H = {H}, d = {d}: planted relation (den {den}, seed {seed}) recovered: {okp} ({tp:.0f} s); "
                   f"perturbed basis on the planted target finds {'nothing' if rpp is None else list(rpp)}; budget (n + 1) log10(2H + 1) + 1 = {budget:.1f}")
            if not REAL:
                print(f"    REAL {name}: skipped (--no-real)", flush=True)
                continue
            rr, tr = o8.pslq_run([f8] + cs, tau, H)
            found, why = False, "not found"
            if rr is not None:
                res = o8.residual(rr, [f8] + cs)
                found = res <= tau and budget + PSLQ_SLACK < d and rr[0] != 0
                why = (f"relation {list(rr)}, residual {nd(res, 3)} (bar {nd(tau, 2)}), budget with H = {H}: {budget:.1f} + {PSLQ_SLACK} < {d}: {budget + PSLQ_SLACK < d}; "
                       f"f8 coefficient {rr[0]}")
                if found:
                    j = next(i for i in range(1, len(rr)) if rr[i] != 0)
                    csr = list(cs)
                    csr[j - 1] = csr[j - 1] + pert
                    rrp, _ = o8.pslq_run([f8] + csr, tau, H)
                    found = found and rrp is None
                    why += f"; perturbed basis (constant {names[j - 1]} + sqrt(2)/1000) finds {'nothing' if rrp is None else list(rrp)}"
                    with mp.workdps(d + VERIFY_EXTRA):
                        csv = [c_ for _, c_ in constants(name)]
                        tv = mp.mpf(10) ** -(d + VERIFY_EXTRA - PSLQ_SLACK)
                        vres = {g: o8.residual(rr, [f8_target(ev, g)] + csv) for g in grids}
                        vok = all(r_ <= tv for r_ in vres.values())
                        found = found and vok
                        why += f"; verified at {d + VERIFY_EXTRA} digits (bar {nd(tv, 2)}): " + ", ".join(f"{g} {nd(r_, 3)}" for g, r_ in vres.items()) + f" -> {vok}"
                    if found:
                        saved[name] = {"names": list(names), "rel": [int(v) for v in rr], "d": d, "H": H}
            print(f"    REAL {name}: " + (("FOUND (counts) " if found else "found but does not count ") + why if rr is not None else f"not found at H = {H} ({tr:.0f} s)"), flush=True)
    if SAVE_REL:
        with open(SAVE_REL, "w") as f:
            json.dump(saved, f, indent=1)


def check_plant():
    """--save=FILE: two "grids" whose Omega0 .. Omega6 are the closed forms and whose Omega7 is built from a planted f8 (an invented ln 2-free combination of the
    13 monomials of weights 3 .. 8), 700 digits; the structure controls run on it (the real relation does not enter them)"""
    mp.mp.dps = 760
    r = random.Random(315999)
    den = 14175
    fr = [v for nm, wt, v in free9() if wt >= 3]
    coef = [r.choice([-1, 1]) * r.randint(10**4, 10**9) for _ in fr]
    f8 = sum(c * v for c, v in zip(coef, fr)) / den
    a = [mp.mpf(1)] + closed_omegas9()
    q4 = 4 * mp.ln2
    om7 = f8 - sum(a[j] * q4 ** (8 - j) / mp.factorial(8 - j) for j in range(8))
    cl = closed_omegas9()
    e = {f"Om{i}": cl[i] for i in range(7)}
    e["Om7"] = om7
    e["grid"] = "planted dps=700"
    out = {"G1": dict(e), "G2": dict(e)}
    with open(SAVE, "w") as f:
        json.dump({n: ({k: (mp.nstr(v, 740) if isinstance(v, mp.mpf) else v) for k, v in x.items()}) for n, x in out.items()}, f, indent=1)
    with open(SAVE + ".planted.json", "w") as f:
        json.dump({"den": den, "coef": coef, "names": [nm for nm, wt, v in free9() if wt >= 3]}, f, indent=1)
    print(f"  planted f8 = (sum c_i b_i)/{den} with c = {coef}; Omega7 = {nd(om7, 30)} written to {SAVE}")


# ---------------------------------------------------------------------------------------------------------------
# the closed form of Omega7 from a relation of f8, and the structure verdict
# ---------------------------------------------------------------------------------------------------------------
SYM_OF = None


def sym_of():
    q, z3, z5, z7, p2, z53 = sp.symbols("q z3 z5 z7 p2 z53")
    return {"pi^8": p2**4, "pi^2 zeta(3)^2": p2 * z3**2, "zeta(3) zeta(5)": z3 * z5, "zeta(5,3)": z53, "pi^4 zeta(3)": p2**2 * z3, "pi^2 zeta(5)": p2 * z5, "zeta(7)": z7,
            "pi^6": p2**3, "zeta(3)^2": z3**2, "pi^2 zeta(3)": p2 * z3, "zeta(5)": z5, "pi^4": p2**2, "zeta(3)": z3, "pi^2": p2, "1": sp.Integer(1)}


ID7_BAR = 520


def check_closed7(ev):
    grids = ["G1", "G2"]
    saved = {}
    for path in REL.split(","):
        saved.update(json.load(open(path)))
    mp.mp.dps = 600
    have = "F7" in saved
    q, z3, z5, z7, p2, z53 = sp.symbols("q z3 z5 z7 p2 z53")
    sof = sym_of()
    if have:
        rel = list(saved["F7"]["rel"])
        names = saved["F7"]["names"]
        if MUTATE == "f8_rel":
            rel[4] += 1
        f8s = -sum(sp.Rational(rel[1 + i], rel[0]) * sof[nm] for i, nm in enumerate(names))
        Om, _ = sym_omegas()
        a = [sp.Integer(1)] + Om + [om6_printed()]
        om7 = sp.expand(f8s - sum(a[j] * (4 * q) ** (8 - j) / sp.factorial(8 - j) for j in range(8)))
        f = sp.lambdify((q, z3, z5, z7, p2, z53), om7, "mpmath")
        zz = zeta53(600) if True else None
        om7n = f(mp.ln2, mp.zeta(3), mp.zeta(5), mp.zeta(7), mp.pi**2, zz)
        ds = {g: digits(om7n, ev[g]["Om7"]) for g in grids}
        report("ID7", all(d_ >= ID7_BAR for d_ in ds.values()),
               f"Omega7 from the ln 2-free relation of F7 (f8 = -sum c_i b_i / c_0, Omega7 = f8 - the exact ln 2 terms) = {nd(om7n, 60)} against the grids, digits: "
               + ", ".join(f"{g} {float(d_):.1f}" for g, d_ in ds.items()) + f" (bar {ID7_BAR}; the working precision here is 600 digits)")
        print("  Omega7 = " + str(om7))
        pol = sp.Poly(om7, q, p2, z3, z5, z7, z53)
        byw = {}
        for mon, co in pol.terms():
            w_ = mon[0] + 2 * mon[1] + 3 * mon[2] + 5 * mon[3] + 7 * mon[4] + 8 * mon[5]
            byw.setdefault(w_, []).append((mon, co))
        for w_ in sorted(byw, reverse=True):
            part = sum(mp.mpf(co.p) / co.q * mp.ln2 ** m[0] * mp.pi ** (2 * m[1]) * mp.zeta(3) ** m[2] * mp.zeta(5) ** m[3] * mp.zeta(7) ** m[4] * zz ** m[5] for m, co in byw[w_])
            print(f"    weight {w_}: {len(byw[w_])} terms, part {nd(part, 20)}")
        print(f"  coefficient of zeta(5,3) in f8: {sp.Rational(-rel[1 + names.index('zeta(5,3)')], rel[0])}; height of the integer vector {max(abs(v) for v in rel)}")
    else:
        report("ID7", False, "no relation of F7 was found (the file has " + (", ".join(saved) or "no relation") + "): f8 is not in the span of the ln 2-free monomials of weight <= 8 at these digits")
    bad = []
    for nm, r in saved.items():
        if nm.startswith("X"):
            nz = [r["names"][i] for i in range(len(r["names"])) if "ln" in r["names"][i] and r["rel"][1 + i] != 0]
            if nz:
                bad.append(f"{nm}: nonzero on {nz}")
    if MUTATE == "ln2_zl" and not bad:
        bad.append("ln2_zl: a nonzero coefficient on a ln 2 monomial put in by the control")
    ok = have and not bad
    report("ZL", ok, "F7 found: " + str(have) + "; relations found in the bases with ln 2 (" + (", ".join(n_ for n_ in saved if n_.startswith("X")) or "none") + "): "
           + ("all with zero coefficients on every ln 2 monomial" if not bad else "; ".join(bad)))


# ---------------------------------------------------------------------------------------------------------------
# the remainder against #195's W^2
# ---------------------------------------------------------------------------------------------------------------
X5_BAR = mp.mpf("1e-14")      # set knowing a development extrapolation of a 177-digit Omega7: 1.6e-19 .. 5.8e-17 on the windows K = 14 .. 44 (dev_05)
R8_BAR = mp.mpf("1e-4")      # #282 estimated Omega7 as 25545476; the development run (N = 120) has 25545344.6, relative 5.1e-6
R9_BAR = mp.mpf(5 * 10**6)  # set knowing dev_05: |rem9 s^8 (10^4) - rem9 s^8 (5103)| = 1.29e6


def check_remainder9(ev):
    mp.mp.dps = 80
    grids = ["G1", "G2"]
    om = [ev[grids[1]][f"Om{i}"] for i in range(1, 8)]
    if MUTATE == "wrong_Om7_series":
        om[6] = om[6] + 10**6
    lst = o8.read_w2(W2FILE)
    if W2RAW:
        lst.update(o8.read_w2raw(W2RAW))

    def x(D, nterms):
        s = mp.mpf(2 * (D - 3))
        rem = mp.mpf(lst[D]) - (s - 4 * mp.ln2 + sum(om[i] / s ** (i + 1) for i in range(nterms)))
        return rem, rem * s ** (nterms + 1)
    print("  W^2 (qd, w2_output.txt) against the series s - 4 ln 2 + Omega1/s + ... (Omega1 .. Omega7 of G2), s = 2(D - 3):")
    print("  D       s       rem8 s^7 (tends to Omega7)    (rem8 s^7 - Omega7) s (tends to Omega8)   rem9 s^8 (tends to Omega8)")
    xs8, xs9 = {}, {}
    for D in (1000, 2000, 5103, 10000):
        s = mp.mpf(2 * (D - 3))
        _, xs8[D] = x(D, 6)
        _, xs9[D] = x(D, 7)
        print(f"  {D:<6d} {2 * (D - 3):<7d} {nd(xs8[D], 14):<28s} {nd((xs8[D] - om[6]) * s, 12):<40s} {nd(xs9[D], 12)}")
    s1, s2 = mp.mpf(2 * (5103 - 3)), mp.mpf(2 * (10000 - 3))
    est8 = (s2 * xs8[10000] - s1 * xs8[5103]) / (s2 - s1)
    report("R8", abs(est8 - om[6]) <= R8_BAR * abs(om[6]), f"two-point extrapolation (s2 x2 - s1 x1)/(s2 - s1) of rem8 s^7 from D = 5103, 10^4: {nd(est8, 8)} against the computed Omega7 "
           f"{nd(om[6], 8)}, relative difference {nd((est8 - om[6]) / om[6], 3)} (bar {R8_BAR}: an estimate, not a bound)")
    est9 = (s2 * xs9[10000] - s1 * xs9[5103]) / (s2 - s1)
    print(f"  reported: two-point extrapolation of rem9 s^8 from D = 5103, 10^4: Omega8 ~ {nd(est9, 8)} (an estimate, not a bound)")
    report("R9", abs(xs9[10000] - xs9[5103]) <= R9_BAR, f"rem9 s^8 at D = 5103 and 10^4 differs by {nd(abs(xs9[10000] - xs9[5103]), 5)} (bar {R9_BAR}: a wrong Omega7 adds (error)·s, about 2e4 per unit)")
    mp.mp.dps = 120
    est = o8.refit_windows(o8.read_w2(W2FILE), closed_omegas9(), 7)
    rows = {n: {K: abs(e - ev[n]["Om7"]) for K, e in est.items()} for n in grids}
    mx = max(max(r.values()) for r in rows.values())
    report("X5", mx <= X5_BAR, "Omega7 against the polynomial extrapolation in 1/s of #195's 64 W^2 values (Omega0 .. Omega6 fixed at the closed forms), "
           "windows of the K largest s, |difference| " + "; ".join(f"{n}: " + ", ".join(f"K = {K} {nd(d, 2)}" for K, d in r.items()) for n, r in rows.items())
           + f" (bar {X5_BAR} on every window; the spread over windows is an estimator range, {nd(mx, 2)}; not a bound)")


def main():
    parse(sys.argv[1:])
    o8.cc = cc9
    o8.closed_omegas = closed_omegas9
    o8.WRONG_OM = o8.WRONG_OM + ("wrong_Om6",)
    mp.mp.dps = 60
    print(f"order9.py  mutate={MUTATE}  parts={PARTS}  mpmath {mp.__version__}  sympy {sp.__version__}", flush=True)
    t0 = time.time()
    if "grid" in PARTS:
        check_grid9()
    if "cache" in PARTS:
        check_cache9()
    if "sym" in PARTS:
        check_sym9()
    if "struct" in PARTS:
        check_struct9()
    if "zv" in PARTS:
        check_zv()
    if "table" in PARTS:
        check_table()
    if "plant" in PARTS:
        check_plant()
    ev = None
    if "evals" in PARTS:
        ev = check_evals9()
        if "values" in PARTS or not GRIDSEL:
            check_values9(ev)
    elif LOAD and any(p in PARTS for p in ("values", "pslq", "remainder", "closed7")):
        mp.mp.dps = 800
        ev = o8.load_ev(LOAD)
        if "values" in PARTS:
            check_values9(ev)
    if "pslq" in PARTS:
        check_pslq9(ev)
    if "closed7" in PARTS:
        check_closed7(ev)
    if "remainder" in PARTS:
        check_remainder9(ev)
    nfail = sum(1 for _, ok in o8.RESULTS if not ok)
    print(f"\n{len(o8.RESULTS) - nfail} passed, {nfail} failed of {len(o8.RESULTS)} checks; {time.time() - t0:.0f} s")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()

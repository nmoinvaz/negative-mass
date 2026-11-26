"""W^2 a^2 of the Ellis wormhole's growing mode at large D to 30+ digits (#195): the driver.

The crate `w2_crate/` (lane-generic Taylor stepper: f64, dd, qd, od) computes W^2(D) of E2 of master.py. This script
runs the checks of the raw README (negative-mass-raw/I195-dh610), fits Omega_1 and Omega_2 of
    W^2 a^2 = 2(D-3) - 4 ln 2 + Omega_1 eps + Omega_2 eps^2 + ...,   eps = 1/(2(D-3)),
and searches integer relations (PSLQ). Long crate runs (the D list in qd and od, the planted case in od) are separate
processes whose logs are passed with --logs; short ones are run here.

    python w2_highprec.py --logs=RAW/run_01_qd.log,RAW/run_02_od.log [--only=E0,S1,...] [--mutate=NAME]

Exit status 1 when a check fails (or, with --mutate, when the control does not fire on its declared targets).
Needs the crate built: cargo build --release -j 2 in w2_crate/.
"""
import argparse
import json
import os
import random
import subprocess
import sys
import time

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "w2_crate", "target", "release", "w2_crate")
mp.mp.dps = 300

DLIST = [50, 54, 59, 64, 70, 76, 83, 90, 98, 107, 116, 126, 137, 149, 162, 177, 192, 209, 227, 247, 269, 292, 318,
         346, 376, 409, 445, 484, 527, 573, 623, 678, 737, 802, 873, 949, 1032, 1123, 1221, 1329, 1445, 1572, 1710,
         1860, 2023, 2201, 2394, 2604, 2832, 3081, 3351, 3645, 3965, 4313, 4691, 5103, 5550, 6037, 6567, 7143, 7770,
         8452, 9193, 10000]
LANE_D = {"f64": 11, "dd": 27, "qd": 56, "od": 110}       # stated digits of the planted checks (README, lane table)
LN2 = mp.log(2)
OM1_CLOSED = 2 * mp.pi**2 / 3 + 8 * LN2**2

# fit cells (README, "Fit")
KS = (10, 12, 14, 16, 18, 20, 24)
KS_K = KS + (28, 32)        # F3k: the grid extended after dev_21 (Amendment A2)
KS2 = (8, 10, 12, 14)
DLOS = (50, 100, 200, 500)
GATE_FACTOR = 100
DATA_FLOOR = mp.mpf("1e-60")

# controls: name -> checks that must fail under it (computed in development runs dev_14, dev_19; raw README Amendment A2).
# Only checks that pass without the control are targets (S4, S4w, S5 and F3 fail without one, by the rulings).
TARGETS = {
    "wrong_eps": ("F3k",),
    "wrong_B": ("E0", "S1", "S2", "S3", "S5s", "O2", "M1", "F3k"),
    "wrong_A_sign": ("E0", "O2", "M1", "F3k"),
    "recip_f64": ("S1", "S2", "S3", "S5s", "F3k"),
    "short_series": ("S1", "S5s", "F3k"),
    "short_contraction": ("S3", "S5s", "F3k"),
    "step_loose": ("S2", "S3", "S5s", "F3k"),
    "wrong_planted": ("S1", "S2", "S3"),
    "fit_wrong_eps": ("F1", "F2", "F3k", "P1"),
    "fit_drop_term": ("F1", "F2"),       # silent on F3k (3 of 6 cells still admitted): reported, dev_23
    "pslq_wrong_basis": ("P1",),         # P2 flips from found to not found, P3 unchanged (not found)
}
# (control, check): the cases (lane, D) in which the check must fail, from the development run (coordinator's ruling: dd at
# D = 5 is silent under short_contraction, -27.33 against -27.0)
CASE_TARGETS = {("short_contraction", "S5s"): ("dd D=50", "dd D=1000", "dd D=10000", "qd D=5", "qd D=50", "qd D=1000",
                                              "qd D=10000")}

CRATE_MUTATIONS = ("wrong_eps", "wrong_B", "wrong_A_sign", "recip_f64", "short_series", "short_contraction", "step_loose",
                   "wrong_planted")
LAST_CELLS = []
RESULTS = []          # (name, status, detail) with status in PASS / FAIL / SKIPPED
POOL = []             # RESULT dicts parsed from crate logs
MUTATE = None
CRATE_MUT = None      # MUTATE when it is one of the crate's controls


def report(name, ok, detail):
    status = "PASS" if ok else "FAIL"
    RESULTS.append((name, status, detail))
    print(f"{status:7s} {name}: {detail}", flush=True)


def skipped(name, why):
    RESULTS.append((name, "SKIPPED", why))
    print(f"SKIPPED {name}: {why}", flush=True)


def parse_logs(paths):
    for p in paths:
        with open(p) as f:
            for line in f:
                if line.startswith("RESULT "):
                    POOL.append(json.loads(line[7:]))


def crate(mode, lane, ds, extra=()):
    """run the crate and return its RESULT dicts"""
    cmd = [BIN, mode, f"--lane={lane}", "--D=" + ",".join(str(d) for d in ds), *extra]
    if CRATE_MUT:
        cmd.append(f"--mutate={CRATE_MUT}")
    out = subprocess.run(["nice", "-n", "10", *cmd], capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed: {out.stderr[-400:]}")
    return [json.loads(line[7:]) for line in out.stdout.splitlines() if line.startswith("RESULT ")]


def pool_w2(lane, d):
    for r in POOL:
        if r["kind"] == "w2" and r["lane"] == lane and r["D"] == d and r.get("T") is not None:
            return r
    return None


def w2_values(lane, ds, allow_run=True):
    """{D: result dict} from the pool, else from a crate run (not under a mutation: the pool is the real data)"""
    got = {}
    miss = []
    for d in ds:
        r = pool_w2(lane, d)
        if r is not None and not CRATE_MUT:
            got[d] = r
        else:
            miss.append(d)
    if miss and allow_run:
        for r in crate("w2", lane, miss, extra=_guess_args(miss)):
            got[r["D"]] = r
    return got


def _guess_args(ds):
    """growth.csv's values as Newton starts for D < 50 (only a start; the result does not depend on it)"""
    csv = read_csv()
    if all(d in csv for d in ds):
        return ["--W2=" + ",".join(str(csv[d][0]) for d in ds)]
    return []


def read_csv():
    out = {}
    with open(os.path.join(HERE, "growth.csv")) as f:
        next(f)
        for line in f:
            d, w, src = line.strip().split(",")
            out[int(d)] = (w, src)
    return out


def digits_between(a, b):
    """-log10 of the relative difference of two decimal strings (or mpf)"""
    a, b = mp.mpf(a), mp.mpf(b)
    if a == b:
        return mp.inf
    return -mp.log10(abs(a - b) / abs(b))


# ---------------------------------------------------------------------------------------------------------------
# E0: the transcription of E2 into the system of the crate, symbolically
# ---------------------------------------------------------------------------------------------------------------
def check_E0():
    import sympy as sp
    l, eps, Om = sp.symbols("l epsilon Omega", positive=True)
    Y0, Y1, Y2 = sp.symbols("Y0 Y1 Y2")
    s = 1 / eps
    n = s / 2 + 1
    E = sp.exp(-l)
    F = 1 - E
    r = sp.exp(eps * l)
    W2 = s + Om
    q2 = F                                   # r0'^2
    rpp = (s / 2) * E / r                    # r0'' = (n - 1) r0^{1 - 2n}
    rp_over_r = sp.Symbol("x")               # r0'/r0 appears squared or times r0'/q only: handle by q2
    # E2 with Y_l = Y1 s r0'/r0, Y_ll = Y2 (s r0'/r0)^2 + Y1 s (r0''/r0 - r0'^2/r0^2); r0' = q, q^2 = F
    q = sp.Symbol("q", positive=True)
    Yl = Y1 * s * q / r
    Yll = Y2 * (s * q / r)**2 + Y1 * s * (rpp / r - q**2 / r**2)
    e2 = Yll - (2 * rpp / q + n * q / r) * Yl + (n * q**2 / r**2 + rpp / r) * Y0 - W2 * Y0
    e2 = e2.subs(q, sp.sqrt(F))
    A = (1 + eps * Om) * r**2 - (sp.Rational(1, 2) + eps) * F - E / 2
    B = (sp.Rational(1, 2) + 2 * eps) * F + E / 2
    mine = (s**2 / r**2) * (F * Y2 - B * Y1 - eps * A * Y0)
    if MUTATE == "wrong_A_sign":
        mine = (s**2 / r**2) * (F * Y2 - B * Y1 - eps * ((1 + eps * Om) * r**2 - (sp.Rational(1, 2) + eps) * F + E / 2) * Y0)
    if MUTATE == "wrong_B":
        mine = (s**2 / r**2) * (F * Y2 - (sp.Rational(1, 2) + 2 * eps) * F * Y1 - eps * A * Y0)
    d = sp.simplify(sp.expand(e2 - mine))
    pt = {l: sp.Rational(7, 5), eps: sp.Rational(1, 9), Om: sp.Rational(-13, 10), Y0: sp.Rational(3, 7), Y1: sp.Rational(-5, 11), Y2: sp.Rational(2, 3)}
    num = sp.N((e2 - mine).subs(pt), 30)
    report("E0", d == 0 and abs(num) < 1e-25, f"E2 in l minus the crate's system, symbolic residual {d}, numeric {sp.N(num, 5)}")


# ---------------------------------------------------------------------------------------------------------------
# the fit
# ---------------------------------------------------------------------------------------------------------------
def fit_one(pts, K, model):
    """least squares of g/eps on eps^(k-1) (k = 1..K) [and eps^(k-1) ln eps (k = 3..K) in M2]; returns (x, max relative residual)"""
    rows, rhs, rel = [], [], []
    for d, g in pts:
        e = mp.mpf(1) / (2 * (d - 3))
        row = [e**(k - 1) for k in range(1, K + 1)]
        if model == "M2":
            row += [e**(k - 1) * mp.log(e) for k in range(3, K + 1)]
        rows.append(row)
        rhs.append(g / e)
        rel.append(abs(g / e))
    A = mp.matrix(rows)
    b = mp.matrix(rhs)
    x, _ = mp.qr_solve(A, b)
    res = A * x - b
    return x, max(abs(res[i]) / rel[i] for i in range(len(pts)))


def fit_cells(data, delta, model="M1", wrong_eps=False, ks=None):
    """data: {D: g}; returns the list of cells with their Omega_1, Omega_2 and whether the gate admits them"""
    ks = (ks or KS) if model == "M1" else KS2
    gate = GATE_FACTOR * max(delta, DATA_FLOOR)
    cells = []
    for dlo in DLOS:
        pts = [(d, g) for d, g in sorted(data.items()) if d >= dlo]
        for K in ks:
            nparam = K if model == "M1" else K + (K - 2)
            if len(pts) < nparam + 4:
                continue
            if wrong_eps:      # control fit_wrong_eps: eps = 1/(D-3) in the fit
                pts_f = [(2 * d - 3, g) for d, g in pts]
            else:
                pts_f = pts
            x, res = fit_one(pts_f, K - 4 if MUTATE == "fit_drop_term" else K, model)
            cells.append({"model": model, "Dlo": dlo, "K": K, "n": len(pts), "Om1": x[0], "Om2": x[1], "Om3": x[2], "Om4": x[3], "res": res,
                          "admitted": res <= gate})
    return cells, gate


def headline(cells, key):
    vals = sorted(c[key] for c in cells if c["admitted"])
    if not vals:
        return None
    return {"median": vals[len(vals) // 2] if len(vals) % 2 else (vals[len(vals) // 2 - 1] + vals[len(vals) // 2]) / 2,
            "range": vals[-1] - vals[0], "n": len(vals)}


def synthetic(noise, with_log=0.0, seed=195):
    """F1/F2: g from Omega_k = A (k+1)! c^k (truncated at its smallest term) on the D list, relative noise"""
    om1, om2 = mp.mpf("10.4233603787385"), mp.mpf(73)
    c = om2 / om1 / 3
    A = om1 / (2 * c)
    rng = random.Random(seed)
    data = {}
    for d in DLIST:
        e = mp.mpf(1) / (2 * (d - 3))
        s, prev = mp.mpf(0), None
        for k in range(1, 80):
            t = A * mp.factorial(k + 1) * c**k * e**k
            if prev is not None and abs(t) > abs(prev):
                break
            s += t
            prev = t
        s += with_log * e**3 * mp.log(e)
        data[d] = s * (1 + noise * (2 * mp.mpf(rng.random()) - 1))
    return data, om1, om2


def check_fits_planted():
    noise = mp.mpf("1e-40")
    data, om1, om2 = synthetic(noise)
    cells, gate = fit_cells(data, noise, "M1", wrong_eps=MUTATE == "fit_wrong_eps")
    h1, h2 = headline(cells, "Om1"), headline(cells, "Om2")
    if not h1:
        report("F1", False, "no admitted cell on planted data")
    else:
        e1, e2 = abs(h1["median"] - om1), abs(h2["median"] - om2)
        ok = (h1["n"] >= 3 and e1 <= mp.mpf("1e-15") and e2 <= mp.mpf("1e-16") * om2
              and abs(h1["median"] - om1) <= h1["range"] and abs(h2["median"] - om2) <= h2["range"])
        report("F1", ok, f"{h1['n']} cells admitted; Omega1 error {mp.nstr(e1, 3)} (range {mp.nstr(h1['range'], 3)}), "
                         f"Omega2 error {mp.nstr(e2, 3)} (range {mp.nstr(h2['range'], 3)})")
    data2, _, _ = synthetic(noise, with_log=mp.mpf("0.5"))
    c1, _ = fit_cells(data2, noise, "M1", wrong_eps=MUTATE == "fit_wrong_eps")
    c2, _ = fit_cells(data2, noise, "M2", wrong_eps=MUTATE == "fit_wrong_eps")
    h2m = headline(c2, "Om2")
    none1 = not any(c["admitted"] for c in c1)
    ok = none1 and h2m is not None and abs(h2m["median"] - om2) <= mp.mpf("1e-10") * om2
    report("F2", ok, f"M1 admits {sum(c['admitted'] for c in c1)} of {len(c1)} cells on data with a log term; M2 admits "
                     f"{h2m['n'] if h2m else 0} and Omega2 error {mp.nstr(abs(h2m['median'] - om2), 3) if h2m else 'n/a'}")


# ---------------------------------------------------------------------------------------------------------------
# PSLQ
# ---------------------------------------------------------------------------------------------------------------
def pslq_found(x, tau, H, expected=None):
    """(found, relation, residual): a relation with max |coeff| <= H, residual <= tau, the same at 10 tau"""
    r = mp.pslq(x, tol=tau, maxcoeff=H, maxsteps=10**6)
    if r is None:
        return False, None, None
    res = abs(sum(c * v for c, v in zip(r, x)))
    ok = res <= tau and max(abs(c) for c in r) <= H
    r10 = mp.pslq(x, tol=10 * tau, maxcoeff=H, maxsteps=10**6)
    ok = ok and r10 is not None and (list(r10) == list(r) or list(r10) == [-c for c in r])
    if expected is not None:
        ok = ok and (list(r) == expected or list(r) == [-c for c in expected])
    return ok, r, res


def budget(n, H, tau, scale):
    """digits a chance relation would need vs digits available: found only if n log10(2H+1) + 1 < -log10(tau/scale)"""
    need = n * mp.log10(2 * H + 1) + 1
    have = -mp.log10(tau / abs(scale))
    return need < have, need, have


# ---------------------------------------------------------------------------------------------------------------
# checks that run the crate
# ---------------------------------------------------------------------------------------------------------------
def planted_rows(lane, ds=(5, 50, 1000, 10000)):
    rows = [r for r in POOL if r["kind"] == "planted" and r["lane"] == lane and r["D"] in ds]
    if (not rows or CRATE_MUT) and lane != "od":
        rows = crate("planted", lane, ds)
    return rows


def check_planted(lanes):
    """S1 to S5 and the added S4w, S5s (raw README, Amendment A2), per lane and D"""
    fails = {k: [] for k in ("S1", "S2", "S3", "S4", "S4w", "S5", "S5s")}
    seen = {k: 0 for k in fails}
    notconv = []
    for lane in lanes:
        d = LANE_D[lane]
        rows = planted_rows(lane)
        if not rows:
            skipped(f"S* ({lane})", "no planted result for this lane (pass its log with --logs)")
            continue
        for r in rows:
            D = r["D"]
            s = 2 * (D - 3)
            tag = f"{lane} D={D}"
            def bad(name, cond, val):
                seen[name] += 1
                if not cond:
                    fails[name].append(f"{tag} ({val})")
            bad("S1", max(r["S1y"], r["S1w"]) <= -d, max(r["S1y"], r["S1w"]))
            bad("S2", max(r["S2y"], r["S2w"]) <= -d, max(r["S2y"], r["S2w"]))
            bad("S3", r["S3"] <= -d, r["S3"])
            for nm, key, conv in (("a", "S4a", r["S4conv"][0]), ("b", "S4b", r["S4conv"][1])):
                if not conv:
                    notconv.append(f"{tag} start {nm}")
                bad("S4", r[key] <= -d, f"start {nm}: {r[key]}")
                # S4w: |W2 - W2*| <= 10^-d W2*, from the relative error of Omega
                werr = r[key] + float(mp.log10(mp.mpf("2.5") / (s + mp.mpf("2.5"))))
                bad("S4w", werr <= -d, f"start {nm}: {werr:.2f}")
            bad("S5", r["S5"] <= r["S5bound"], f"{r['S5']} vs {r['S5bound']}")
            b5s = float(mp.log10(100 * mp.mpf(s)**2 * mp.exp(-mp.mpf(r["T5"])) + mp.mpf(10)**(-d))) if lane != "f64" else None
            if lane != "f64":
                bad("S5s", r["S5"] <= b5s, f"{r['S5']} vs {b5s:.2f}")
    for k, v in fails.items():
        if seen[k] == 0:
            continue
        report(k, not v, f"{seen[k] - len(v)} of {seen[k]} cases meet it" + ("" if not v else "; NOT MET: " + "; ".join(v)))
    if notconv:
        print(f"        (Newton's own stop test, two steps below 10^-d relative, was not reached in 60 iterations at: {'; '.join(notconv)})")


def check_series():
    ok, det = True, []
    for lane, ds in (("qd", (1000, 10000)),):
        for r in crate("series", lane, ds):
            det.append(f"D={r['D']}: Y - closed form = {r['diff']:.2e} = {r['diff_over_eps3']:.2f} eps^3")
            ok = ok and abs(r["diff_over_eps3"]) <= 100
    report("O2", ok, "throat series at l = 1 against Y0 + eps Y1 + eps^2 Y2 of README.md (|diff| <= 100 eps^3): " + "; ".join(det))


def check_start_robust():
    """R4: Newton from the start shifted by -0.5 and +0.5 reaches the same root to 1e-50, qd"""
    ds = list(range(4, 12)) + [50, 1000]
    bad, det = [], []
    base = w2_values("qd", ds)
    for sh in (-0.5, 0.5):
        got = {}
        for d in ds:
            args = _guess_args([d]) + [f"--shift={sh}"]
            got[d] = crate("w2", "qd", [d], extra=args)[0]
        for d in ds:
            dg = digits_between(got[d]["W2"], base[d]["W2"])
            conv = got[d]["converged"]
            if not conv or dg < 50:
                bad.append(f"D={d} shift {sh:+}: {'not converged, ' if not conv else ''}W2 = {mp.nstr(mp.mpf(got[d]['W2']), 12)} agrees {mp.nstr(dg, 3)} digits")
    report("R4", not bad, "Newton from the start +-0.5, D = 4..11, 50, 1000, qd, same root to 1e-50" +
           ("" if not bad else "; NOT MET: " + "; ".join(bad)))


def check_master():
    """M1: master.py's l-integration (scipy) at D = 5..8 against the crate"""
    sys.path.insert(0, HERE)
    import master
    ds = (5, 6, 7, 8)
    ours = w2_values("qd", ds)
    worst, det = 0.0, []
    for d in ds:
        roots = master.growing_mode(d)
        w = float(mp.mpf(ours[d]["W2"]))
        best = min(roots, key=lambda x: abs(x - w)) if roots else float("nan")
        worst = max(worst, abs(best - w))
        det.append(f"D={d}: master.py {best:.6f} ({len(roots)} roots), crate {w:.8f}")
    report("M1", worst <= 2e-3, f"worst difference {worst:.1e} (bar 2e-3); " + "; ".join(det))


# ---------------------------------------------------------------------------------------------------------------
# growth.csv, the lanes, and the fit of record
# ---------------------------------------------------------------------------------------------------------------
def modes_r_tight():
    """{D: value} of modes_r.py at X 50, rho0/10, rtol 1e-13, read from logs_modes_r.txt (not rerun)"""
    out = {}
    with open(os.path.join(HERE, "logs_modes_r.txt")) as f:
        for line in f:
            if line.startswith("D =") and "rtol 1e-13:" in line:
                d = int(line.split("D =")[1].split(":")[0])
                out[d] = mp.mpf(line.split("rtol 1e-13:")[1].split(",")[0])
    return out


def check_growth():
    csv = read_csv()
    ds = list(csv)
    ours = w2_values("qd", ds)
    tight = modes_r_tight()
    rows, bad = [], []
    for d in ds:
        c = mp.mpf(csv[d][0])
        w = mp.mpf(ours[d]["W2"])
        diff = w - c
        ok = abs(diff) <= mp.mpf("5e-13")
        if not ok:
            bad.append(d)
        rows.append((d, csv[d][1], c, w, diff, ok))
    r12 = lambda x: mp.floor(x * 10**12 + mp.mpf(1) / 2) / 10**12      # noqa: E731
    selftest = all(abs(r12(mp.mpf(ours[d]["W2"])) - mp.mpf(ours[d]["W2"])) <= mp.mpf("5e-13") for d in ds)
    report("G1", not bad, f"{len(ds) - len(bad)} of {len(ds)} rows within 5e-13 of the printed value; NOT MET at D = {bad}; "
                          f"the comparator passes the crate's own values printed to 12 decimals: {selftest}")
    print("\n  growth.csv against the crate (qd), and the tighter modes_r.py run of logs_modes_r.txt (reported, not a gate)")
    print("  D       source          csv                   crate (qd)                crate - csv   tight - csv   moved by fraction of gap")
    for d, src, c, w, diff, ok in rows:
        t = tight.get(d)
        if t is not None:
            mv = t - c
            frac = mv / diff if diff != 0 else mp.nan
            tcol = f"{mp.nstr(mv, 3):>12s}   {mp.nstr(frac, 3):>8s}"
        else:
            tcol = "         n/a        n/a"
        print(f"  {d:<6d}  {src:<14s}  {mp.nstr(c, 16):<21s} {mp.nstr(w, 21):<25s} {mp.nstr(diff, 3):>11s}  {tcol}  {'ok' if ok else 'differs'}")
    print()


def check_lanes():
    q = w2_values("qd", DLIST, allow_run=False)
    o = w2_values("od", DLIST, allow_run=False)
    common = [d for d in DLIST if d in q and d in o]
    if len(common) < len(DLIST):
        skipped("L1", f"qd and od agree on {len(common)} of {len(DLIST)} D in the logs given")
        return None
    worst, rows = mp.inf, []
    delta = mp.mpf(0)
    for d in DLIST:
        dg = digits_between(q[d]["W2"], o[d]["W2"])
        rows.append((d, dg))
        worst = min(worst, dg)
        gj = mp.mpf(q[d]["W2"]) - 2 * (d - 3) + 4 * LN2
        delta = max(delta, abs(mp.mpf(q[d]["W2"]) - mp.mpf(o[d]["W2"])) / abs(gj))
    report("L1", worst >= 30, f"qd against od at {len(DLIST)} values of D: fewest agreeing digits {mp.nstr(worst, 5)} (bar 30)")
    return delta


def table_w2():
    q = w2_values("qd", DLIST, allow_run=False)
    o = w2_values("od", DLIST, allow_run=False)
    print("  W^2 a^2 (qd, 60 digits shown) and the digits shared with od; l_inf and Taylor steps of the qd run")
    for d in DLIST:
        if d not in q:
            continue
        dg = digits_between(q[d]["W2"], o[d]["W2"]) if d in o else None
        print(f"  D={d:<6d} {mp.nstr(mp.mpf(q[d]['W2']), 60):<68s} digits qd/od: {mp.nstr(dg, 4) if dg is not None else 'n/a':>6s}  "
              f"l_inf {q[d]['l_inf']}, {q[d]['steps']} steps, {q[d]['secs']:.1f} s")
    print()


def _fit_summary(label, data, delta, ks, verbose):
    cells, gate = fit_cells(data, delta, "M1", wrong_eps=MUTATE == "fit_wrong_eps", ks=ks)
    model = "M1"
    h1, h2 = headline(cells, "Om1"), headline(cells, "Om2")
    if h1 is None:
        cells, gate = fit_cells(data, delta, "M2")
        model = "M2"
        h1, h2 = headline(cells, "Om1"), headline(cells, "Om2")
    global LAST_CELLS
    LAST_CELLS = cells
    if verbose:
        print(f"  fit cells of {label}, model {model} (gate: relative residual <= {mp.nstr(gate, 3)}; delta_data = {mp.nstr(delta, 3)})")
        print("  model  Dlo    K   n    Omega1                      Omega2                       rel. residual   admitted")
        for c in cells:
            print(f"  {c['model']}     {c['Dlo']:<5d} {c['K']:<3d} {c['n']:<3d} {mp.nstr(c['Om1'], 24):<27s} {mp.nstr(c['Om2'], 24):<28s} {mp.nstr(c['res'], 3):<14s}  {c['admitted']}")
        print()
    if h1 is None:
        report(label, False, "no cell of M1 or M2 is admitted by the gate: the series model does not describe the data")
        return None
    d1 = h1["median"] - OM1_CLOSED
    ok1 = h1["range"] <= mp.mpf("1e-12") and abs(d1) <= mp.mpf("1e-12") and abs(d1) <= 3 * h1["range"]
    ok2 = h2["range"] <= mp.mpf("1e-16") * abs(h2["median"])
    digs2 = -mp.log10(h2["range"] / abs(h2["median"])) if h2["range"] > 0 else mp.inf
    det = (f"{h1['n']} admitted cell(s), model {model}; Omega1 {mp.nstr(h1['median'], 22)}, range {mp.nstr(h1['range'], 2)}, "
           f"Omega1 - (2 pi^2/3 + 8 ln^2 2) = {mp.nstr(d1, 3)} (bars 1e-12, 3 x range); "
           f"Omega2 {mp.nstr(h2['median'], 24)}, range {mp.nstr(h2['range'], 2)} = {mp.nstr(digs2, 3)} digits (bar 16)")
    if h1["n"] == 1:
        det += ("; one admitted cell has range 0, so the '<= 3 x range' clause cannot hold: a defect of the criterion's "
                "wording, not of the data (the difference from the closed form is " + mp.nstr(abs(d1), 3) + ")")
    if model != "M1":
        det += "; MODEL M1 FAILED THE GATE, M2 used"
    report(label, ok1 and ok2 and model == "M1", det)
    return h1, h2


def check_fit_of_record(delta):
    q = w2_values("qd", DLIST, allow_run=bool(CRATE_MUT))
    if len(q) < len(DLIST) or delta is None:
        skipped("F3", "needs the qd values of the D list and L1's lane disagreement")
        skipped("F3k", "needs the qd values of the D list and L1's lane disagreement")
        return None
    data = {d: mp.mpf(q[d]["W2"]) - 2 * (d - 3) + 4 * LN2 for d in DLIST}
    _fit_summary("F3", data, delta, None, False)
    hk = _fit_summary("F3k", data, delta, KS_K, True)
    if hk is not None:
        print("  headline from F3k (an estimator range over admitted cells, not a bound); status: numerical evidence\n")
        print_tail(hk, LAST_CELLS)
    return hk


def print_tail(hk, cells):
    """Omega3, Omega4 from the same cells, and the four-term series' remainder times s^3 at the growth.csv D (#26's third test)"""
    for k in ("Om3", "Om4"):
        h = headline(cells, k)
        print(f"  {k.replace('Om', 'Omega')} = {mp.nstr(h['median'], 14)}  range over the {h['n']} admitted cells {mp.nstr(h['range'], 2)}")
    csv = read_csv()
    ds = [200, 500, 1000, 2000]
    ours = w2_values("qd", ds, allow_run=False)
    if len(ours) < len(ds):
        return
    om1, om2 = hk[0]["median"], hk[1]["median"]
    print("\n  the four-term series 2(D-3) - 4 ln 2 + Omega1/s + Omega2/s^2 (s = 2(D-3)) against W^2: remainder x s^3 (tends to Omega3)")
    print("  D       crate (qd)         growth.csv (modes_r.py)")
    for d in ds:
        s_ = 2 * (d - 3)
        ser = s_ - 4 * LN2 + om1 / s_ + om2 / s_**2
        print(f"  {d:<6d}  {mp.nstr((mp.mpf(ours[d]['W2']) - ser) * s_**3, 10):<18s} {mp.nstr((mp.mpf(csv[d][0]) - ser) * s_**3, 10)}")
    print()


def check_pslq(h):
    if h is None:
        skipped("P1", "no fit of record")
        return
    h1, h2 = h
    om1, om2 = h1["median"], h2["median"]
    pert = mp.mpf("1e-3")
    basis1 = [om1, mp.pi**2, LN2**2 + (pert if MUTATE == "pslq_wrong_basis" else 0)]
    tau, H = mp.mpf("1e-11"), 100
    okb, need, have = budget(3, H, tau, om1)
    found, rel, res = pslq_found(basis1, tau, H, expected=[3, -2, -24])
    ctrl, _, _ = pslq_found([om1, mp.pi**2, LN2**2 + pert], tau, H)
    report("P1", found and okb and not ctrl,
           f"{{Omega1, pi^2, ln^2 2}} tau {mp.nstr(tau, 1)}, H {H}: relation {list(rel) if rel else 'not found'}, residual "
           f"{mp.nstr(res, 3) if res is not None else 'n/a'}, budget {mp.nstr(need, 3)} < {mp.nstr(have, 3)} digits: {okb}; "
           f"perturbed basis finds a relation: {ctrl}; identification, not a derivation")
    z3 = mp.zeta(3)
    pert2 = pert if MUTATE == "pslq_wrong_basis" else 0
    b2 = [om2, z3 + pert2, mp.pi**2 * LN2, LN2**3]
    tau2, H2 = mp.mpf("1e-13"), 1000
    okb2, need2, have2 = budget(4, H2, tau2, om2)
    found2, rel2, res2 = pslq_found(b2, tau2, H2)
    ctrl2, _, _ = pslq_found([om2, z3 + pert, mp.pi**2 * LN2, LN2**3], tau2, H2)
    if found2 and not okb2:
        verdict = "inconclusive (budget)"
    else:
        verdict = "found" if found2 else "not found"
    report("P2", (verdict == "found" or verdict == "not found") and not ctrl2 and not (found2 and MUTATE == "pslq_wrong_basis"),
           f"{{Omega2, zeta(3), pi^2 ln 2, ln^3 2}} tau {mp.nstr(tau2, 1)}, H {H2}: {verdict}" +
           (f", relation {list(rel2)} residual {mp.nstr(res2, 3)}" if rel2 else "") +
           f"; budget {mp.nstr(need2, 3)} vs {mp.nstr(have2, 3)} digits; perturbed basis finds a relation: {ctrl2}")
    tau3, H3 = 10 * h2["range"], 100
    if tau3 == 0:
        tau3 = mp.mpf("1e-30")
    b3 = [om2, z3 + pert2, mp.pi**2 * LN2, LN2**3, mp.pi**2, LN2**2, LN2, mp.mpf(1)]
    okb3, need3, have3 = budget(8, H3, tau3, om2)
    found3, rel3, res3 = pslq_found(b3, tau3, H3)
    ctrl3, _, _ = pslq_found([om2, z3 + pert, mp.pi**2 * LN2, LN2**3, mp.pi**2, LN2**2, LN2, mp.mpf(1)], tau3, H3)
    verdict3 = "inconclusive (budget)" if not okb3 else ("found" if found3 else "not found")
    report("P3", not ctrl3 and not (found3 and okb3 and MUTATE == "pslq_wrong_basis"),
           f"{{Omega2, zeta(3), pi^2 ln 2, ln^3 2, pi^2, ln^2 2, ln 2, 1}} tau {mp.nstr(tau3, 2)}, H {H3}: {verdict3}" +
           (f", relation {list(rel3)} residual {mp.nstr(res3, 3)}" if rel3 else "") +
           f"; budget {mp.nstr(need3, 3)} vs {mp.nstr(have3, 3)} digits; perturbed basis finds a relation: {ctrl3}")


# ---------------------------------------------------------------------------------------------------------------
def main():
    global MUTATE, CRATE_MUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--logs", default="", help="comma-separated crate logs (RESULT lines) of the long runs")
    ap.add_argument("--only", default="", help="comma-separated check names")
    ap.add_argument("--mutate", default=None)
    ap.add_argument("--delta", default=None, help="development runs without an od log: the lane disagreement to use")
    a = ap.parse_args()
    MUTATE = a.mutate
    CRATE_MUT = MUTATE if MUTATE in CRATE_MUTATIONS else None
    if a.logs:
        parse_logs(a.logs.split(","))
    only = set(filter(None, a.only.split(",")))
    want = lambda *names: not only or any(n in only for n in names)       # noqa: E731
    t0 = time.time()
    rev = subprocess.run(["git", "-C", HERE, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    print(f"w2_highprec.py  code {rev}  mutate={MUTATE}  crate={os.path.relpath(BIN, HERE)}\n")
    if want("E0"):
        check_E0()
    if want("S1", "S2", "S3", "S4", "S4w", "S5", "S5s"):
        check_planted(["f64", "dd", "qd", "od"] if not MUTATE else ["f64", "dd", "qd"])
    if want("O2"):
        check_series()
    if want("M1"):
        check_master()
    if want("R4") and not MUTATE:
        check_start_robust()
    if want("G1") and not CRATE_MUT:
        check_growth()
    delta = None
    if want("L1", "F3", "F3k", "P1", "P2", "P3"):
        if CRATE_MUT:
            delta = mp.mpf(a.delta) if a.delta else mp.mpf(0)      # no od run under a control: --delta (the record's value)
        else:
            delta = check_lanes()
            if delta is None and a.delta:
                delta = mp.mpf(a.delta)
            table_w2()
    if want("F1", "F2"):
        check_fits_planted()
    h = None
    if want("F3", "F3k", "P1", "P2", "P3"):
        h = check_fit_of_record(delta)
        if want("P1", "P2", "P3"):
            check_pslq(h)
    npass = sum(r[1] == "PASS" for r in RESULTS)
    nfail = sum(r[1] == "FAIL" for r in RESULTS)
    nskip = sum(r[1] == "SKIPPED" for r in RESULTS)
    print(f"\n{npass} passed, {nfail} failed, {nskip} skipped ({', '.join(r[0] for r in RESULTS if r[1] == 'SKIPPED')}) "
          f"of {len(RESULTS)} checks; {time.time() - t0:.0f} s")
    if MUTATE:
        failed = {r[0] for r in RESULTS if r[1] == "FAIL"}
        tg = set(TARGETS.get(MUTATE, ()))
        cells_ok = all(all(c + " (" in r[2] for c in cells) for (m, ck), cells in CASE_TARGETS.items() if m == MUTATE
                       for r in RESULTS if r[0] == ck)
        fires = bool(tg) and tg <= failed and cells_ok
        print(f"control {MUTATE}: failed checks {sorted(failed)}; declared targets {sorted(tg)}"
              + (f"; cells {CASE_TARGETS[(MUTATE, 'S5s')]}" if (MUTATE, "S5s") in CASE_TARGETS else "") + "; "
              + ("fires" if fires else "DOES NOT FIRE (or no targets declared yet)"))
        sys.exit(0 if fires else 1)
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()

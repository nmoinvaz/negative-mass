#!/usr/bin/env python3
"""The coordinator's independent check of #40 (ind40 of negative-mass-raw I40-mac-review, where its criteria, amendments,
runs and mutations are recorded; this copy differs only in the default of --outtxt). The lifetime law for family S by
the coordinator's own double-null evolution, which shares no code with dn.py, and a refit of #40's points.

    python3 threshold_check.py [--mutate=oracle-slope|fit-pstar|time-unit|sphere-coef] [--only=K0,K1,K2,K3,K4] [--D=4,5,6]
                     [--out=DIR] [--raw=I40-dh610/runs] [--outtxt=threshold_output.txt]
"""
import argparse
import json
import math
import os
import re
import sys
import time

import numba as nb
import numpy as np
from scipy.integrate import solve_ivp

ap = argparse.ArgumentParser()
ap.add_argument("--mutate", default=None, choices=["oracle-slope", "fit-pstar", "time-unit", "sphere-coef"])
ap.add_argument("--only", default="K0,K1,K2,K3,K4")
ap.add_argument("--D", default="4,5,6")
ap.add_argument("--out", default=None, help="directory for the evolution points (jsonl per D and h)")
ap.add_argument("--raw", default=os.path.expanduser("~/Source/negative-mass-raw/I40-dh610/runs"))
ap.add_argument("--outtxt", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "threshold_output.txt"))
OPT = ap.parse_args()
MUT = OPT.mutate
ONLY = OPT.only.split(",")
DS = [int(x) for x in OPT.D.split(",")]

SQ2 = math.sqrt(2.0)
W2 = {4: 1.39705243433350684, 5: 2.98495893035176025, 6: 4.68662054311946479}        # #195's crate, a = 1
GRID = {4: ((0.02, 0.01, 0.005), 28.0), 5: ((0.01 / SQ2, 0.005 / SQ2, 0.0025 / SQ2), 18.0), 6: ((0.005, 0.0025, 0.00125), 15.0)}
ES = [4.0 + 0.5 * k for k in range(13)]                                                # e = 4 ... 10
WIN = (1e-9, 1e-4)
RLO, RHI = 0.5, 2.0
XC, WP = 3.0, 0.5

RESULTS = []


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


# ------------------------------------------------------------------------------------------------------------------
# fit, bisection and series (shared by the oracle and the evolution; none of it is #40's)
# ------------------------------------------------------------------------------------------------------------------
def fit(pts, pstar, win, sides=(-1, 1)):
    """Common slope, one intercept per side, closed-form within-side normal equations. pts: (side, p, fate, T)."""
    by = {}
    for s, p, fate, T in pts:
        d = abs(p - pstar)
        if s in sides and fate == s and T is not None and win[0] * (1 - 1e-6) <= d <= win[1] * (1 + 1e-6):
            by.setdefault(s, []).append((math.log(1.0 / d), T))
    num = den = 0.0
    cnt = 0
    for s, xy in by.items():
        x = np.array([a for a, _ in xy])
        y = np.array([b for _, b in xy])
        num += float(((x - x.mean()) * (y - y.mean())).sum())
        den += float(((x - x.mean()) ** 2).sum())
        cnt += len(xy)
    return (num / den if den > 0 else float("nan")), cnt


def bisect(run, lo, hi, tol=1e-14):
    log = []
    for p, want in ((lo, -1), (hi, 1)):
        r = run(p)
        log.append(r)
        if r["fate"] != want:
            raise RuntimeError(f"bracket does not straddle: fate {r['fate']} at {p}, wanted {want}")
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:
            break
        r = run(mid)
        log.append(r)
        if r["fate"] == 1:
            hi = mid
        elif r["fate"] == -1:
            lo = mid
        else:
            raise RuntimeError(f"no event at p = {mid!r}")
    return lo, hi, log


def series(run, pstar):
    out = []
    for e in ES:
        for s in (-1, 1):
            r = run(pstar + s * 10.0 ** (-e))
            r.update(side=s, e=e)
            out.append(r)
    return out


# ------------------------------------------------------------------------------------------------------------------
# the evolution (the coordinator's own code)
# ------------------------------------------------------------------------------------------------------------------
def ray_data(n, N, h, p):
    """r and phi on the ray x- = 0 at s = i h: background and constraint integrated together, f = 0 (the gauge)."""
    q = math.sqrt(n * (n - 1.0))
    g0 = math.exp(-(XC / WP) ** 2)

    def rhs(s, y):
        rb, rbs, _, r, rs = y
        g = math.exp(-((s - XC) / WP) ** 2)
        phbs = q / (SQ2 * rb ** n)
        phs = phbs + p * (-2.0 * (s - XC) / WP ** 2) * g
        return [rbs, (n - 1.0) / (2.0 * rb ** (2 * n - 1)), phbs, rs, (phs * phs / n) * r]
    s = np.arange(N + 1) * h
    sol = solve_ivp(rhs, [0.0, s[-1]], [1.0, 0.0, 0.0, 1.0, 0.0], method="DOP853", rtol=1e-13, atol=1e-14, t_eval=s)
    if not sol.success:
        raise RuntimeError(sol.message)
    _, _, phb, r, _ = sol.y
    phi = phb + p * (np.exp(-((s - XC) / WP) ** 2) - g0)
    return np.ascontiguousarray(r), np.ascontiguousarray(phi)


@nb.njit(cache=False)
def march(n, h, r0, p0, c1, rlo, rhi):
    """Half domain x+ >= x- of the symmetric member, row by row in x-. Returns the diagonal r and f, the fate and the
    last row. Cells i < cut of the previous row are valid."""
    N = r0.shape[0] - 1
    pr = r0.copy()
    pf = np.zeros(N + 1)
    pp = p0.copy()
    cr = np.zeros(N + 1)
    cf = np.zeros(N + 1)
    cp = np.zeros(N + 1)
    dr = np.zeros(N + 1)
    df = np.zeros(N + 1)
    dr[0] = pr[0]
    df[0] = 0.0
    nn = float(n)
    h2 = h * h
    ih = 0.5 / h
    cut = N + 1
    fate = 0
    jend = 0
    for j in range(1, N + 1):
        newcut = cut
        for i in range(j, cut):
            rE = pr[i]
            fE = pf[i]
            qE = pp[i]
            rS = pr[i - 1]
            fS = pf[i - 1]
            qS = pp[i - 1]
            if i == j:
                rW = rE
                fW = fE
                qW = -qE
            else:
                rW = cr[i - 1]
                fW = cf[i - 1]
                qW = cp[i - 1]
            rN = rE + rW - rS
            fN = fE + fW - fS
            qN = qE + qW - qS
            for it in range(3):
                rc = 0.25 * (rN + rE + rW + rS)
                fc = 0.25 * (fN + fE + fW + fS)
                rp = (rN - rW + rE - rS) * ih
                rm = (rN - rE + rW - rS) * ih
                qp = (qN - qW + qE - qS) * ih
                qm = (qN - qE + qW - qS) * ih
                ruv = -c1 * (rp * rm + 0.5 * math.exp(-fc)) / rc
                fuv = (nn / rc) * ruv - qp * qm
                quv = -nn * (rp * qm + rm * qp) / (2.0 * rc)
                rN = rE + rW - rS + h2 * ruv
                fN = fE + fW - fS + h2 * fuv
                qN = qE + qW - qS + h2 * quv
            if not (np.isfinite(rN) and np.isfinite(fN) and np.isfinite(qN)) or rN < 0.05 or fN < -18.42:
                newcut = i
                break
            cr[i] = rN
            cf[i] = fN
            cp[i] = qN
        jend = j
        if newcut <= j:
            fate = 0
            break
        dr[j] = cr[j]
        df[j] = cf[j]
        if cr[j] <= rlo:
            fate = -1
            break
        if cr[j] >= rhi:
            fate = 1
            break
        cut = newcut
        for i in range(j, cut):
            pr[i] = cr[i]
            pf[i] = cf[i]
            pp[i] = cp[i]
    return dr[:jend + 1], df[:jend + 1], fate, jend


def run_member(n, h, X, p):
    t0 = time.time()
    N = int(round(X / h))
    r0, p0 = ray_data(n, N, h, p)
    c1 = float(n) if MUT == "sphere-coef" else n - 1.0
    dr, df, fate, jend = march(n, h, r0, p0, c1, RLO, RHI)
    out = dict(A=p, fate=int(fate), T=None, tau=None, j=int(jend), secs=round(time.time() - t0, 3))
    if fate != 0:
        k = jend
        target = RLO if fate == -1 else RHI
        dt = SQ2 * h
        fr = (target - dr[k - 1]) / (dr[k] - dr[k - 1])
        T = dt * (k - 1 + fr)
        lap = np.exp(-0.5 * df[:k + 1])
        tau = dt * (0.5 * lap[0] + lap[1:k - 1].sum() + 0.5 * lap[k - 1]) if k > 1 else 0.0
        tau += dt * fr * 0.5 * (lap[k - 1] + (lap[k - 1] + fr * (lap[k] - lap[k - 1])))
        sc = 1.0 / SQ2 if MUT == "time-unit" else 1.0
        out.update(T=T * sc, tau=tau * sc)
    return out


# ------------------------------------------------------------------------------------------------------------------
# K0: planted
# ------------------------------------------------------------------------------------------------------------------
if "K0" in ONLY:
    p0, s0, cc = 3.7e-6, 0.8, {-1: 1.0, 1: 3.0}

    def oracle(p):
        f = 1 if p > p0 else -1
        return dict(A=p, fate=f, T=cc[f] + s0 * math.log(1.0 / abs(p - p0)))
    lo, hi, _ = bisect(oracle, -0.01, 0.01)
    ps = 0.5 * (lo + hi)
    rows = series(oracle, ps)
    s, cnt = fit([(r["side"], r["A"], r["fate"], r["T"]) for r in rows], ps, WIN)
    want = 0.81 if MUT == "oracle-slope" else s0
    check("K0", abs(ps - p0) <= 1e-13 and abs(s / want - 1) <= 1e-9,
          f"p* {ps!r} against {p0} (diff {ps - p0:.1e}); slope {s:.12f} against {want} (rel {s / want - 1:+.1e}, {cnt} points)")

# ------------------------------------------------------------------------------------------------------------------
# K1: refit of #40's points
# ------------------------------------------------------------------------------------------------------------------
PRINTED = {}
for line in open(OPT.outtxt):
    m = re.match(r"\s+SL D=(\d) (\w) h = ([0-9.]+): p\*_h = \S+; primary slope ([0-9.]+)", line)
    if m:
        PRINTED[(int(m[1]), m[2], float(m[3]))] = float(m[4])


def load40(D, fam, h):
    d = os.path.join(OPT.raw, f"D{D}_{fam}_h{h:.6g}")
    summ = json.load(open(os.path.join(d, "summary.json")))
    pts = [json.loads(x) for x in open(os.path.join(d, "points.jsonl"))]
    pts = [(r["side"], r["A"], r["fate"], r["T"]) for r in pts if r.get("phase") == "series"]
    return summ, pts


if "K1" in ONLY:
    worst, notes = 0.0, []
    for (D, fam, h), sp in sorted(PRINTED.items()):
        summ, pts = load40(D, fam, h)
        ps = 0.5 * (summ["lo"] + summ["hi"]) + (1e-9 if MUT == "fit-pstar" else 0.0)
        s, cnt = fit(pts, ps, (1e-11, 1e-3))
        worst = max(worst, abs(s - sp))
        notes.append(f"D{D}{fam} h{h:.6g} {s:.6f}/{sp:.6f}")
    check("K1", len(PRINTED) == 26 and worst <= 1e-6,
          f"{len(PRINTED)} jobs; largest |refit - printed| {worst:.1e}; " + ", ".join(notes))

# ------------------------------------------------------------------------------------------------------------------
# K2 to K4: the evolution
# ------------------------------------------------------------------------------------------------------------------
if any(k in ONLY for k in ("K2", "K3", "K4")):
    EV = {}
    for D in DS:
        n = D - 2
        hs, X = GRID[D]
        Winv = 1.0 / math.sqrt(W2[D])
        for h in hs:
            t0 = time.time()
            lo, hi, blog = bisect(lambda p: run_member(n, h, X, p), -0.01, 0.01)
            ps = 0.5 * (lo + hi)
            rows = series(lambda p: run_member(n, h, X, p), ps)
            pts = [(r["side"], r["A"], r["fate"], r["T"]) for r in rows]
            s, cnt = fit(pts, ps, WIN)
            s5, _ = fit(pts, ps, (1e-9, 1e-5))
            sc, _ = fit(pts, ps, WIN, (-1,))
            se, _ = fit(pts, ps, WIN, (1,))
            stau, _ = fit([(r["side"], r["A"], r["fate"], r["tau"]) for r in rows], ps, WIN)
            fitted = [r for r in rows if WIN[0] * (1 - 1e-6) <= 10 ** (-r["e"]) <= WIN[1] * (1 + 1e-6)]
            own = all(r["fate"] == r["side"] for r in fitted)
            floor = [(r["e"], r["side"], r["fate"], None if r["T"] is None else round(r["T"], 4)) for r in rows if r["e"] > 9.0]
            EV[(D, h)] = dict(pstar=ps, lo=lo, hi=hi, slope=s, cnt=cnt, own=own)
            print(f"  D={D} h={h:.6g}: p*_h {ps:.12e} ({len(blog)} bisection runs); slope {s:.6f} ({cnt} pts; rel {s / Winv - 1:+.2e}); "
                  f"[1e-9,1e-5] {s5:.6f}; collapse {sc:.6f}, expansion {se:.6f}; in tau {stau:.6f}; own side {own}; floor {floor}; "
                  f"{time.time() - t0:.0f} s", flush=True)
            if OPT.out:
                os.makedirs(OPT.out, exist_ok=True)
                with open(os.path.join(OPT.out, f"D{D}_h{h:.6g}.jsonl"), "w") as fh:
                    for r in blog:
                        fh.write(json.dumps(dict(phase="bisect", **r)) + "\n")
                    for r in rows:
                        fh.write(json.dumps(dict(phase="series", **r)) + "\n")
        ph = [EV[(D, h)]["pstar"] for h in hs]
        if all(abs(x) > 0 for x in ph):
            print(f"  D={D}: p*_h {ph[0]:.6e}, {ph[1]:.6e}, {ph[2]:.6e}; orders with p* = 0: {math.log2(ph[0] / ph[1]):.3f}, "
                  f"{math.log2(ph[1] / ph[2]):.3f}; Richardson {math.log2(abs((ph[0] - ph[1]) / (ph[1] - ph[2]))):.3f}", flush=True)
    if "K2" in ONLY:
        ok, notes = True, []
        for D in DS:
            hs, _ = GRID[D]
            Winv = 1.0 / math.sqrt(W2[D])
            a, b = EV[(D, hs[1])], EV[(D, hs[2])]
            okD = abs(a["slope"] / Winv - 1) <= 0.01 and abs(b["slope"] / Winv - 1) <= 0.01 and abs(a["slope"] / b["slope"] - 1) <= 0.005 \
                and a["own"] and b["own"]
            ok &= okD
            notes.append(f"D={D}: {a['slope']:.6f}, {b['slope']:.6f} against 1/W {Winv:.6f} (rel {a['slope'] / Winv - 1:+.2e}, "
                         f"{b['slope'] / Winv - 1:+.2e}; agreement {a['slope'] / b['slope'] - 1:+.2e}; own side {a['own'] and b['own']})")
        check("K2", ok, "; ".join(notes))
    if "K3" in ONLY:
        ok, notes = True, []
        for D in [d for d in DS if d in (4, 5)]:
            hs, _ = GRID[D]
            Winv = 1.0 / math.sqrt(W2[D])
            s = EV[(D, hs[2])]["slope"]
            ok &= abs(s / Winv - 1) <= 2e-3
            notes.append(f"D={D} h={hs[2]:.6g}: {s:.6f}, rel {s / Winv - 1:+.2e}")
        check("K3", ok, "; ".join(notes))
    if "K4" in ONLY:
        ok, notes = True, []
        for D in DS:
            hs, _ = GRID[D]
            summ, pts = load40(D, "S", float(f"{hs[2]:.6g}"))
            s40, cnt40 = fit(pts, 0.5 * (summ["lo"] + summ["hi"]), WIN)
            s = EV[(D, hs[2])]["slope"]
            ok &= abs(s / s40 - 1) <= 3e-3
            notes.append(f"D={D}: ind40 {s:.6f}, #40 refit over [1e-9, 1e-4] {s40:.6f} ({cnt40} pts), rel {s / s40 - 1:+.2e}")
        check("K4", ok, "; ".join(notes))

bad = [nm for nm, ok in RESULTS if not ok]
print(f"\n{len(RESULTS)} criteria, {len(RESULTS) - len(bad)} pass; failed: {bad}; mutation {MUT}", flush=True)
sys.exit(1 if bad else 0)

#!/usr/bin/env python3
"""Lifetime law at the Ellis threshold (#40, #150): bisection to the discrete threshold p*_h, the lifetime against ln|p - p*_h|, and its
slope against 1/W, at D = 4, 5, 6. Uses dn.py and checks.py of #39 unchanged. The criteria are CRIT and MUTATIONS below, as declared in
negative-mass-raw I40-dh610/README.md before the runs of record.

    python3 threshold.py job --D=4 --family=S --h=0.02 --out=DIR [--ed] [--stall] [--cons]   bisect, series, event-definition, stall and constraint runs
    python3 threshold.py g6 --out=DIR                       D = 6 growth control: checks.check_G1 at n = 4, declared h and X
    python3 threshold.py planted                            the planted lifetime series through the same bisection and fit code
    python3 threshold.py analyze --raw=DIR [--mutate=NAME]  every check and table, from the stored runs (planted recomputed)

Families (field 0 is the ghost, 8 pi G = 1, a = 1, t = (x+ + x-)/sqrt2, rays through the throat corner):
    S  symmetric two-sided ghost pulse (A on '+', mirrored on '-')        passes through the static wormhole
    O  one-sided ghost pulse (A on '+')                                   passes through the static wormhole too (not generic)
    G  one-sided ghost pulse A (tuned) plus a fixed normal pulse B on the same ray, eps = (-1, +1)   added by the coordinator's ruling
The lifetime is T = t of the event cell: collapse = first row with a trapped cell (r_+ and r_- below -tol); expansion = first row whose diagonal
cell is anti-trapped (both above tol) with r >= rexp. Fates: collapse for p < p*_h, expansion for p > p*_h (as in checks.bisect_A).
"""
import json
import math
import os
import struct
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import checks   # noqa: E402
import dn       # noqa: E402

SQ2 = math.sqrt(2.0)
W_EXACT = checks.W_EXACT        # #195: sqrt(W^2 a^2) for n = 2, 3, 4

# criteria: declared in the raw README (I40-dh610) after the development runs dev_01 to dev_NN and before the first run of record
SQ2H = 1.0 / SQ2
CRIT = dict(
    trap_tol=1e-2, rexp=2.0, rcut=0.05, fcut=25.0,            # events and excision as #39's runs of record
    pulse=dict(xc=3.0, w=0.5),                                  # Gaussian pulses (kind 'phi'), as checks.fate_run
    G_B=0.1,                                                    # the normal pulse of family G (added by the coordinator's ruling)
    H={2: (0.02, 0.01, 0.005), 3: (0.01 * SQ2H, 0.005 * SQ2H, 0.0025 * SQ2H), 4: (0.005, 0.0025, 0.00125)},   # h, coarse to fine, ratio 2
    X=dict(S={2: 28.0, 3: 18.0, 4: 15.0}, O={2: 28.0, 3: 18.0, 4: 15.0}, G={2: 28.0, 3: 22.0, 4: 28.0}),   # rectangle side: sqrt2 X must exceed the longest lifetime
    bracket=dict(S=(-0.01, 0.01), O=(-0.01, 0.01), G=(0.005, 0.03)),       # fate(lo) = collapse, fate(hi) = expansion
    bisect_tol=1e-14,                                           # hi - lo at the end of the bisection (the brief asks for 1e-10)
    es=[3.0 + 0.5 * k for k in range(19)],                      # series p* -/+ 10^-e, e = 3 .. 12 (e > 11 is the roundoff floor, not fitted)
    windows=[(1e-11, 1e-3), (1e-11, 1e-4), (1e-11, 1e-5), (1e-10, 1e-3), (1e-10, 1e-4)],   # fit windows in delta; the first is primary (8 decades)
    fm_e=11.0,                                                  # the fate-monotone check covers e <= fm_e (above it roundoff decides)
    order_lo=1.75, order_hi=2.25,                               # observed order of p*_h, as #39's order band
    agree=0.005,                                                # two finest h: primary slopes agree within 0.5 % (as G1)
    ed_events={"tol0.1": dict(tol=0.1), "rexp3": dict(rexp=3.0)}, ed_es=[4.0 + k for k in range(8)], ed_window=1,
    stall_e=8.0,                                                # the stall and formation runs at p* -/+ 1e-8
    not_run={(4, "G"): 2},                                      # (n, family): index of the h that is not run (cost about 4 h, over the 2 h stop; accepted)
    planted_seeds=40, planted_cover=0.9,                        # coverage of the planted control over noise seeds
    planted_sigma=0.01,                                         # Gaussian noise of the planted series (plus quantization q = sqrt2 h)
    g6_h=(0.005 / math.sqrt(3.0), 0.0025 / math.sqrt(3.0)), g6_X=12.0,
)

MUTATIONS = {   # name: the checks that must fire; computed in the development runs dev_11 to dev_17 (negative-mass-raw I40-dh610/dev), D = 4 family S
    "time_unit": ["PL", "PX", "SL", "G6"],                  # T scaled by sqrt2 in the fits (checks.fit_growth too); G6 from its own mutated run
    "pstar_off": ["FM", "PL", "PX", "RA", "SL"],            # p*_h shifted by 1e-9 in the fits
    "bisect_loose": ["BR", "FM", "RA"],                     # bisection stopped at 1e-9
    "ed_unit": ["ED"],                                      # the alternative events timed in x- units
    "pstar_h1": ["PS"],                                     # p*_h replaced by p*_h h_1/h (a first-order seed)
    "no_corrector": ["SL"],                                 # dn.py: predictor-only update; fires at h = 0.02, 0.01 and not at 0.005
    "theta_coef": ["JOB"], "ghost_sign": ["JOB"], "finder_sign": ["JOB"], "ellis_charge": ["JOB"], "scalar_coef": ["JOB"],
    "data_constraint_sign": ["JOB"],                        # (the six JOB ones are dn.py's: the bracket no longer straddles the threshold, the job fails)
}
ONLY = None                                                     # analyze --only=4S,5O,planted,g6: restrict the checks

results = []


def check(cid, text, ok, detail=""):
    ok = None if ok is None else bool(ok)
    results.append((cid, ok))
    tag = "SKIP" if ok is None else ("pass" if ok else "FAIL")
    print(f"[{tag}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


# --------------------------------------------------------------------------------------------------------------------
# evolution of one family member
# --------------------------------------------------------------------------------------------------------------------
def family(fam, A):
    """(pulses, symmetric, eps) of the member with tuned amplitude A."""
    pu = CRIT["pulse"]
    ghost = dict(field=0, ray="+", A=A, xc=pu["xc"], w=pu["w"])
    if fam == "S":
        return [ghost], True, (-1.0,)
    if fam == "O":
        return [ghost], False, (-1.0,)
    if fam == "G":
        return [ghost, dict(field=1, ray="+", A=CRIT["G_B"], xc=pu["xc"], w=pu["w"])], False, (-1.0, 1.0)
    raise ValueError(fam)


def evolve_member(n, fam, h, X, A, store=False, stop=True, tol=None, rexp=None):
    pulses, sym, eps = family(fam, A)
    N = int(round(X / h))
    d = dn.ellis_rays(n, N, h, pulses=pulses, symmetric=sym, eps=eps)
    res = dn.evolve(d, stop=stop, store=store, tol=tol or CRIT["trap_tol"], rexp=rexp or CRIT["rexp"], rcut=CRIT["rcut"], fcut=CRIT["fcut"])
    return d, res


def point(n, fam, h, X, A, **kw):
    """One run to its event: fate (-1 collapse, +1 expansion, 0 none through X), T = t of the event cell, event cell (j, i)."""
    t0 = time.time()
    d, res = evolve_member(n, fam, h, X, A, **kw)
    ev = res.events
    out = dict(A=A, fate=0, T=None, j=None, i=None)
    if ev:
        key = next(iter(ev))
        out.update(fate=-1 if key == "collapse" else 1, T=dn.event_time(res, key)[2], j=int(ev[key][0]), i=int(ev[key][1]))
    out["secs"] = round(time.time() - t0, 3)
    return out


# --------------------------------------------------------------------------------------------------------------------
# bisection, series and fit: shared by the evolution and by the planted oracle
# --------------------------------------------------------------------------------------------------------------------
def bisect(run, lo, hi, tol, sink=None):
    """Bisection of the fate between lo (collapse) and hi (expansion) to hi - lo <= tol. run(p) -> dict(fate, ...). Every run goes to sink."""
    for p, want in ((lo, -1), (hi, 1)):
        r = run(p)
        if sink:
            sink(dict(phase="bisect", **r))
        if r["fate"] != want:
            raise RuntimeError(f"bracket [{lo}, {hi}] does not straddle the threshold: fate {r['fate']} at {p}, wanted {want}")
    stopped = None
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:
            stopped = "float"
            break
        r = run(mid)
        if sink:
            sink(dict(phase="bisect", **r))
        if r["fate"] == 1:
            hi = mid
        elif r["fate"] == -1:
            lo = mid
        else:
            stopped = "static through X"
            break
    return lo, hi, stopped


def pstar_of(lo, hi):
    return 0.5 * (lo + hi) + (1e-9 if dn.mut("pstar_off") else 0.0)


def series(run, pstar, es, sink=None, phase="series", fate_tag=True):
    """Runs at p = p* -/+ 10^-e for each e; returns rows (side, e, fate, T, A)."""
    rows = []
    for e in es:
        for s in (-1, 1):
            p = pstar + s * 10.0 ** (-e)
            r = run(p)
            r.update(phase=phase, side=s, e=float(e))
            if sink:
                sink(dict(r))
            rows.append(r)
    return rows


def delta_of(A, pstar):
    return abs(A - pstar)


def fit_joint(rows, pstar, window, sides=(-1, 1)):
    """T = c_side + s ln(1/delta) over the points with side in sides and window[0] <= delta <= window[1]. Common slope s, one intercept per side.
    Returns (s, SE of s, rms residual, n). The unit of T is scaled by sqrt2 under the time_unit mutation."""
    sc = SQ2 if dn.mut("time_unit") else 1.0
    pts = [(r["side"], delta_of(r["A"], pstar), r["T"] * sc) for r in rows
           if r["T"] is not None and r["side"] in sides and r["fate"] == r["side"]]
    lo, hi = window
    pts = [(s, d, T) for s, d, T in pts if lo * (1 - 1e-6) <= d <= hi * (1 + 1e-6)]
    present = sorted({s for s, _, _ in pts})
    if len(pts) < 4 or not present:
        return float("nan"), float("nan"), float("nan"), len(pts)
    M = np.array([[1.0 if s == sd else 0.0 for sd in present] + [math.log(1.0 / d)] for s, d, _ in pts])
    y = np.array([T for _, _, T in pts])
    coef, *_ = np.linalg.lstsq(M, y, rcond=None)
    res = y - M @ coef
    dof = max(1, len(y) - M.shape[1])
    s2 = float(res @ res) / dof
    cov = s2 * np.linalg.inv(M.T @ M)
    return float(coef[-1]), float(math.sqrt(cov[-1, -1])), float(math.sqrt(float(res @ res) / len(y))), len(y)


def fit_set(rows, pstar):
    """The declared fit choices: for each window the joint fit and each side alone. Returns the fits, the primary one (joint, first window) with
    its OLS standard error, and the fit spread: the envelope of every choice's slope +- 2 OLS standard errors (window and side choice, plus noise)."""
    allfits = []
    for w in CRIT["windows"]:
        for sides in ((-1, 1), (-1,), (1,)):
            s, se, rms, n = fit_joint(rows, pstar, w, sides)
            allfits.append(dict(window=w, sides=sides, slope=s, se=se, rms=rms, n=n))
    prim = allfits[0]
    ok = [f for f in allfits if np.isfinite(f["slope"]) and np.isfinite(f["se"])]
    env = (min(f["slope"] - 2 * f["se"] for f in ok), max(f["slope"] + 2 * f["se"] for f in ok)) if ok else (float("nan"), float("nan"))
    return allfits, prim, env


def verdict_inside(env, target):
    return bool(np.isfinite(env[0]) and env[0] <= target <= env[1])


# --------------------------------------------------------------------------------------------------------------------
# job: one (D, family, h)
# --------------------------------------------------------------------------------------------------------------------
class Sink:
    def __init__(self, path):
        self.f = open(path, "a")

    def __call__(self, rec):
        self.f.write(json.dumps(rec) + "\n")
        self.f.flush()


def job(n, fam, h, out, ed=False, stall=False, cons=False, X=None):
    X = X or CRIT["X"][fam][n]
    os.makedirs(out, exist_ok=True)
    sink = Sink(os.path.join(out, "points.jsonl"))
    t0 = time.time()
    run = lambda p: point(n, fam, h, X, p)
    lo, hi = CRIT["bracket"][fam]
    summ = dict(D=n + 2, family=fam, h=h, X=X, N=int(round(X / h)), mutation=dn.MUT, B=CRIT["G_B"] if fam == "G" else None)
    try:
        lo, hi, stopped = bisect(run, lo, hi, 1e-9 if dn.mut("bisect_loose") else CRIT["bisect_tol"], sink)
        summ.update(lo=lo, hi=hi, stopped=stopped, bisect_secs=round(time.time() - t0, 1))
        json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)
        ps = 0.5 * (lo + hi)                      # the experiment's p*_h; the fits use pstar_of (mutable)
        series(run, ps, CRIT["es"], sink)
        summ["series_secs"] = round(time.time() - t0, 1)
        if ed:
            for name, kw in CRIT["ed_events"].items():
                kw = dict(kw)
                alt = lambda p, kw=kw: point(n, fam, h, X, p, **kw)
                if dn.mut("ed_unit"):                                # the alternative events timed in x- units (T / sqrt2), [SH02]'s factor
                    alt = lambda p, kw=kw: (lambda r: dict(r, T=None if r["T"] is None else r["T"] / SQ2))(point(n, fam, h, X, p, **kw))
                series(alt, ps, CRIT["ed_es"], sink, phase="ed:" + name)
            summ["ed_secs"] = round(time.time() - t0, 1)
        if cons:
            json.dump(constraint_runs(n, fam, h, X, lo, hi), open(os.path.join(out, "constraints.json"), "w"), indent=1)
            summ["cons_secs"] = round(time.time() - t0, 1)
        if stall:
            json.dump(stall_runs(n, fam, h, X, ps), open(os.path.join(out, "stall.json"), "w"), indent=1)
            summ["stall_secs"] = round(time.time() - t0, 1)
        summ["status"] = "done"
    except Exception as e:           # recorded; analyze reports the job as failed
        import traceback
        traceback.print_exc()
        summ["status"] = "failed: " + repr(e)
    summ["secs"] = round(time.time() - t0, 1)
    json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print(json.dumps(summ), flush=True)
    return summ["status"] == "done"


# --------------------------------------------------------------------------------------------------------------------
# stored runs (store=True, stop=True: rows up to the event are filled; nothing outside the event's past is read)
# --------------------------------------------------------------------------------------------------------------------
def stored(n, fam, h, X, A):
    d, res = evolve_member(n, fam, h, X, A, store=True, stop=True)
    key = next(iter(res.events), None)
    if key is None:
        return d, res, None, None, None
    j0, i0 = res.events[key]
    return d, res, key, int(j0), int(i0)


def constraint_stats(res, i0, j0, chunk=500):
    """max |C_++|, max |C_--| and the largest sum of the absolute values of their terms (|r_++| + |f_+ r_+| + |(r/n) sum eps phi_+^2| and the minus
    version) over the nodes of the event cell's past: i <= i0 + 1, j <= j0. dn.constraints does the sums; chunks of rows j with one row of overlap
    bound the memory. `repro` is the largest difference between dn.constraints and the same sum assembled here from its terms."""
    from types import SimpleNamespace
    d, G, h = res.d, res.G, res.d.h
    n, eps = d.n, d.eps
    I = i0 + 3                                     # nodes 0 .. i0 + 2: C_++ at i = 1 .. i0 + 1
    out = dict(pp=0.0, mm=0.0, spp=0.0, smm=0.0, repro=0.0, cells=0)
    for a in range(0, j0 + 1, chunk):
        b = min(a + chunk, j0 + 1)                 # this chunk's nodes j = a .. b - 1
        lo, hi = max(0, a - 1), min(d.N + 1, b + 1)
        Gv, dv = G[:, :I, lo:hi], res.dead[:I, lo:hi]
        Cpp, Cmm = dn.constraints(SimpleNamespace(d=d, G=Gv, dead=dv))
        r, f = Gv[0], Gv[1]
        sp = sum(eps[k] * ((Gv[2 + k][2:, :] - Gv[2 + k][:-2, :]) / (2 * h)) ** 2 for k in range(len(eps)))
        t1, t2, t3 = (r[2:, :] - 2 * r[1:-1, :] + r[:-2, :]) / h ** 2, (f[2:, :] - f[:-2, :]) / (2 * h) * (r[2:, :] - r[:-2, :]) / (2 * h), (r[1:-1, :] / n) * sp
        sm = sum(eps[k] * ((Gv[2 + k][:, 2:] - Gv[2 + k][:, :-2]) / (2 * h)) ** 2 for k in range(len(eps)))
        u1, u2, u3 = (r[:, 2:] - 2 * r[:, 1:-1] + r[:, :-2]) / h ** 2, (f[:, 2:] - f[:, :-2]) / (2 * h) * (r[:, 2:] - r[:, :-2]) / (2 * h), (r[:, 1:-1] / n) * sm
        # C_++: rows i = 1 .. i0 + 1 (index i - 1), columns j in [a, b); C_--: j = 1 .. N - 1 (index j - lo - 1), rows i = 0 .. i0 + 1
        cj = np.arange(lo, hi)
        mp = (cj >= a) & (cj < b)
        P, SP = Cpp[: i0 + 1][:, mp], (np.abs(t1) + np.abs(t2) + np.abs(t3))[: i0 + 1][:, mp]
        vp = np.isfinite(P)
        dj = np.arange(lo + 1, hi - 1)
        mm_ = (dj >= a) & (dj < b) & (dj <= j0)
        Mm, SM = Cmm[: i0 + 2][:, mm_], (np.abs(u1) + np.abs(u2) + np.abs(u3))[: i0 + 2][:, mm_]
        vm = np.isfinite(Mm)
        if vp.any():
            out["pp"] = max(out["pp"], float(np.max(np.abs(P[vp]))))
            out["spp"] = max(out["spp"], float(np.max(SP[vp])))
            out["repro"] = max(out["repro"], float(np.max(np.abs((t1 + t2 + t3)[: i0 + 1][:, mp][vp] - P[vp]))))
        if vm.any():
            out["mm"] = max(out["mm"], float(np.max(np.abs(Mm[vm]))))
            out["smm"] = max(out["smm"], float(np.max(SM[vm])))
            out["repro"] = max(out["repro"], float(np.max(np.abs((u1 + u2 + u3)[: i0 + 2][:, mm_][vm] - Mm[vm]))))
        out["cells"] += int(vp.sum() + vm.sum())
    out["rel_pp"] = out["pp"] / out["spp"] if out["spp"] > 0 else float("nan")
    out["rel_mm"] = out["mm"] / out["smm"] if out["smm"] > 0 else float("nan")
    return out


def constraint_runs(n, fam, h, X, lo, hi):
    """The null constraints along the two runs that bracket p*_h at the end of the bisection (reported, not gated)."""
    out = {}
    for name, A in (("lo", lo), ("hi", hi)):
        d, res, key, j0, i0 = stored(n, fam, h, X, A)
        rec = dict(A=A, fate=key, event=[j0, i0])
        if key is not None:
            rec.update(constraint_stats(res, i0, j0))
        out[name] = rec
        del res, d
    return out


def throat_partial(res, jmax):
    """dn.throat_min on the slices t = const, restricted to nodes j <= jmax (the rows a stopped run has filled): (t, r_min, l_min)."""
    d, G, h = res.d, res.G, res.d.h
    R, N = G[0], d.N
    out = []
    for s in range(2, min(2 * N - 1, 2 * jmax + 1)):
        lo, hi = max(0, s - N, s - jmax), min(N, s)
        if hi - lo < 2:
            continue
        i = np.arange(lo, hi + 1)
        v = R[i, s - i]
        if np.any(res.dead[i, s - i]):
            continue
        k = int(np.argmin(v))
        if k == 0 or k == len(v) - 1:
            continue
        y0, y1, y2 = v[k - 1], v[k], v[k + 1]
        den = y0 - 2 * y1 + y2
        dk = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        ip = i[k] + dk
        out.append(((s * h + d.x0p + d.x0m) / SQ2, y1 - 0.25 * (y0 - y2) * dk, (ip - (s - ip)) * h / SQ2))
    return np.array(out)


def stall_runs(n, fam, h, X, pstar):
    """At p* -/+ 10^-stall_e: the event cell's offset from the diagonal, the formation values at the first trapped cell (collapse side) and the
    throat locus during the near-critical stall."""
    out = {}
    for s in (-1, 1):
        A = pstar + s * 10.0 ** (-CRIT["stall_e"])
        d, res, key, j0, i0 = stored(n, fam, h, X, A)
        rec = dict(A=A, side=s, event=None if key is None else [j0, i0])
        if key is not None:
            rec["fate"] = key
            rec["T"] = dn.event_time(res, key)[2]
            rec["event_offset_cells"] = i0 - j0                              # event cell minus diagonal, in cells of x+
            rec["event_offset_l"] = (i0 - j0) * h / SQ2                      # the same in proper distance l = (x+ - x-)/sqrt2
            G = res.G
            S_, E_, W_, N_ = G[:, i0, j0], G[:, i0 + 1, j0], G[:, i0, j0 + 1], G[:, i0 + 1, j0 + 1]
            c = 0.25 * (S_ + E_ + W_ + N_)
            rp, rm = (E_[0] + N_[0] - S_[0] - W_[0]) / (2 * h), (W_[0] + N_[0] - S_[0] - E_[0]) / (2 * h)
            rec["r_cell"] = float(c[0])
            rec["M_cell"] = float(c[0] ** (n - 1) * (1 + 2 * math.exp(c[1]) * rp * rm) / 2)    # Misner-Sharp M at the cell, dn.mass_nodes' formula
            rec["M_hor"] = float(c[0] ** (n - 1) / 2)                                           # its value on r_+ = 0 (F1's definition)
            th = throat_partial(res, j0)
            if len(th):
                rec["stall"] = stall_of(th[:, 0], th[:, 1], th[:, 2], rec["T"])
        out[str(s)] = rec
        del res, d
    return out


def stall_of(t, r, l, Tev):
    """The near-critical stall of the throat locus (minimum of r on the slices t = const): the time in [sqrt2 xc + 2, Tev - 2] where |dr/dt| is
    smallest after a 1-unit running mean; r and l there, the span over which r stays within 2 % of it, and the drift dl/dt of the locus over that span."""
    t0 = SQ2 * CRIT["pulse"]["xc"] + 2.0
    m = (t >= t0) & (t <= Tev - 2.0)
    if m.sum() < 10:
        return None
    dt = float(np.median(np.diff(t)))
    k = max(1, int(round(1.0 / dt)))
    rs = np.convolve(r, np.ones(k) / k, mode="same")
    dr = np.abs(np.gradient(rs, t))
    idx = np.where(m)[0]
    ref = idx[int(np.argmin(dr[idx]))]
    inside = np.abs(r - r[ref]) <= 0.02 * r[ref]
    a = ref
    while a > 0 and inside[a - 1]:
        a -= 1
    b = ref
    while b < len(r) - 1 and inside[b + 1]:
        b += 1
    v = float(np.polyfit(t[a:b + 1], l[a:b + 1], 1)[0]) if b - a >= 5 else float("nan")      # drift of the throat locus, dl/dt, over the stall
    return dict(t=float(t[ref]), r=float(r[ref]), l=float(l[ref]), span=float(t[b] - t[a]), t_from=float(t[a]), t_to=float(t[b]), v=v)


# --------------------------------------------------------------------------------------------------------------------
# D = 6 growth control (checks.check_G1 at n = 4)
# --------------------------------------------------------------------------------------------------------------------
def g6(out):
    os.makedirs(out, exist_ok=True)
    checks.CRIT["G1_h"][4] = CRIT["g6_h"]
    checks.CRIT["G1_X"][4] = CRIT["g6_X"]
    checks.check_G1(cases=(4,))
    json.dump(dict(info=checks.INFO.get("G1_D6"), results=checks.results, h=CRIT["g6_h"], X=CRIT["g6_X"], mutation=dn.MUT),
              open(os.path.join(out, "g6.json"), "w"), indent=1, default=float)
    return all(ok for _, ok in checks.results)


# --------------------------------------------------------------------------------------------------------------------
# planted lifetime series (#150): a synthetic type I law through the same bisection, series and fit code
# --------------------------------------------------------------------------------------------------------------------
def oracle(slope, c, pstar, sigma, q, seed=1):
    """run(p) of a synthetic threshold: fate = side of p*, T = c[side] + slope ln(1/|p - p*|) + Gaussian noise (a function of the bits of p,
    so a run is repeatable), rounded to the quantum q of the event row."""
    def run(p):
        d = p - pstar
        side = -1 if d < 0 else 1
        key = struct.unpack("<Q", struct.pack("<d", p))[0]
        g = float(np.random.default_rng([key, seed]).normal()) * sigma
        T = c[side] + slope * math.log(1.0 / max(abs(d), 1e-300)) + g
        return dict(A=p, fate=side, T=round(T / q) * q, j=None, i=None)
    return run


def planted_case(n, slope, pstar=4.7e-5, seed=1):
    """Bisect and fit the synthetic series for D = n + 2; returns the recovered p*, the fits and the envelope."""
    h = CRIT["H"][n][1]
    q = SQ2 * h
    run = oracle(slope, {-1: -0.79, 1: 4.05}, pstar, CRIT["planted_sigma"], q, seed)
    lo, hi, stopped = bisect(run, -0.01, 0.01, CRIT["bisect_tol"])
    rows = series(run, 0.5 * (lo + hi), CRIT["es"])
    ps = pstar_of(lo, hi)
    fits, prim, env = fit_set(rows, ps)
    return dict(n=n, q=q, lo=lo, hi=hi, pstar_rec=ps, pstar=pstar, fits=fits, prim=prim, env=env, stopped=stopped)


def check_planted():
    seeds = range(1, CRIT["planted_seeds"] + 1)
    for n in (2, 3, 4):
        W = W_EXACT[n]
        cov, cov_p, rej, errs = 0, 0, 0, []
        for sd in seeds:
            c = planted_case(n, 1.0 / W, seed=sd)
            cov += verdict_inside(c["env"], 1.0 / W)
            cov_p += abs(c["pstar_rec"] - c["pstar"]) <= 2 * CRIT["bisect_tol"] and c["stopped"] is None
            x = planted_case(n, 1.0 / (SQ2 * W), seed=sd)
            rej += not verdict_inside(x["env"], 1.0 / W)
            errs.append(c["prim"]["slope"] * W - 1)
            if sd == 1:
                c1, x1 = c, x
        m = len(seeds)
        print(f"    PL D={n + 2}: planted slope {1 / W:.6f}, q = {c1['q']:.5f}, sigma = {CRIT['planted_sigma']}; seed 1: p* recovered {c1['pstar_rec']:.15e} (planted {c1['pstar']:.15e}), "
              f"primary slope {c1['prim']['slope']:.6f} +- {c1['prim']['se']:.6f} (OLS SE), spread [{c1['env'][0]:.6f}, {c1['env'][1]:.6f}]; over {m} seeds: "
              f"1/W inside the spread in {cov}, p* recovered in {cov_p}, primary rel error mean {np.mean(errs):+.1e}, rms {np.std(errs):.1e}", flush=True)
        check(f"PL[D={n + 2}]", f"planted type I series through the same bisection and fit code: p* recovered and 1/W inside the fit spread in at least "
              f"{CRIT['planted_cover']:.0%} of {m} noise seeds", cov >= CRIT["planted_cover"] * m and cov_p == m, f"{cov} of {m} inside, p* {cov_p} of {m}")
        print(f"    PX D={n + 2}: planted slope 1/(sqrt2 W) = {1 / (SQ2 * W):.6f}; seed 1 fitted {x1['prim']['slope']:.6f}, spread [{x1['env'][0]:.6f}, {x1['env'][1]:.6f}] "
              f"against 1/W = {1 / W:.6f}; rejected in {rej} of {m} seeds", flush=True)
        check(f"PX[D={n + 2}]", "a planted series with the time unit scaled by sqrt2 (slope 1/(sqrt2 W)) is rejected by the comparison with 1/W in every seed", rej == m,
              f"{rej} of {m}")


# --------------------------------------------------------------------------------------------------------------------
# analysis of the stored runs
# --------------------------------------------------------------------------------------------------------------------
FAMS = ("S", "O", "G")


def load_jobs(raw):
    jobs = {}
    for name in sorted(os.listdir(raw)):
        p = os.path.join(raw, name)
        sp = os.path.join(p, "summary.json")
        if not (os.path.isdir(p) and os.path.exists(sp)):
            continue
        s = json.load(open(sp))
        pts = [json.loads(l) for l in open(os.path.join(p, "points.jsonl"))] if os.path.exists(os.path.join(p, "points.jsonl")) else []
        s["points"] = pts
        s["stall"] = json.load(open(os.path.join(p, "stall.json"))) if os.path.exists(os.path.join(p, "stall.json")) else None
        s["dir"] = name
        jobs[(s["D"], s["family"], s["h"])] = s
    return jobs


def fmtp(x, f="{:.3e}"):
    return f.format(x) if x is not None and np.isfinite(x) else "n/a"


def analyze(raw):
    global RAW_ROOT
    RAW_ROOT = raw
    jobs = load_jobs(raw)
    print(f"analyze: {raw}: {len(jobs)} jobs; mutation {dn.MUT}")
    if ONLY is None or "planted" in ONLY:
        check_planted()
    # D = 6 growth control
    gp = os.path.join(raw, "g6", "g6.json")
    if ONLY is not None and "g6" not in ONLY:
        pass
    elif os.path.exists(gp):
        g = json.load(open(gp))
        info = g["info"] or []
        W = W_EXACT[4]
        wf = [o["Wfit"] for o in info]
        rel = [w / W - 1 for w in wf]
        ok = len(info) == 2 and all(abs(r) <= checks.CRIT["G1_bar"] for r in rel) and abs(wf[0] / wf[1] - 1) <= checks.CRIT["G1_agree"]
        print(f"    G6 D=6: h = {', '.join(f'{x:.6f}' for x in g['h'])}, X = {g['X']}; W_fit {', '.join(f'{w:.6f}' for w in wf)} against W = {W:.6f} "
              f"(rel {', '.join(f'{r:+.2e}' for r in rel)}); mutation in that run: {g['mutation']}")
        check("G6", f"D = 6 growth control: tuned growth within {checks.CRIT['G1_bar']:.0%} of W at two resolutions that agree within {checks.CRIT['G1_agree']:.1%}",
              ok and all(v for _, v in g["results"]), f"results {g['results']}")
    else:
        check("G6", "D = 6 growth control", False, f"{gp} missing")
    summary = {}
    for n in (2, 3, 4):
        for fam in FAMS:
            if ONLY is not None and f"{n + 2}{fam}" not in ONLY:
                continue
            hs = CRIT["H"][n]
            tag = f"D={n + 2} {fam}"
            sl = []
            pst = []
            for k, h in enumerate(hs):
                if CRIT["not_run"].get((n, fam)) == k:
                    print(f"    {tag} h = {h:.6g}: not run (declared: about 4 h, over the 2 h stop)")
                    sl.append(None)
                    pst.append(None)
                    continue
                j = jobs.get((n + 2, fam, h))
                if j is None or j.get("status") != "done":
                    check(f"JOB[{tag} h={h:.6g}]", "job finished", False, "missing" if j is None else j.get("status"))
                    sl.append(None)
                    pst.append(None)
                    continue
                lo, hi = j["lo"], j["hi"]
                ps = pstar_of(lo, hi)
                allrows = [r for r in j["points"] if r["phase"] == "series"]
                rows = [r for r in allrows if r["e"] <= CRIT["fm_e"]]
                bad = [r for r in rows if r["T"] is None or (1 if r["A"] > ps else -1) != r["fate"]]
                floor = [r for r in allrows if r["e"] > CRIT["fm_e"] and (r["T"] is None or (1 if r["A"] > ps else -1) != r["fate"])]
                check(f"BR[{tag} h={h:.6g}]", f"the bisection ended with hi - lo <= {CRIT['bisect_tol']:g} and no run static through X", hi - lo <= CRIT["bisect_tol"] and j.get("stopped") is None,
                      f"hi - lo = {hi - lo:.2e}, stopped {j.get('stopped')}, {j.get('B') and 'B = ' + str(j['B']) or ''}")
                check(f"FM[{tag} h={h:.6g}]", f"every series point with |p - p*_h| >= 1e-{CRIT['fm_e']:g} lands on its own side of p*_h (fate monotone) and ends in an event", not bad,
                      f"{len(bad)} of {len(rows)} points off; below that, {len(floor)} of {len(allrows) - len(rows)} points off the side (roundoff floor, not fitted)")
                fits, prim, env = fit_set([r for r in allrows if r["T"] is not None], ps)
                W = W_EXACT[n]
                ok = verdict_inside(env, 1.0 / W)
                sc, ex = [f for f in fits if f["window"] == CRIT["windows"][0] and f["sides"] == (-1,)][0], [f for f in fits if f["window"] == CRIT["windows"][0] and f["sides"] == (1,)][0]
                print(f"    SL {tag} h = {h:.6g}: p*_h = {ps:.12e}; primary slope {prim['slope']:.6f} +- {prim['se']:.6f} (OLS SE), n = {prim['n']}, rms {prim['rms']:.4f}; "
                      f"collapse {sc['slope']:.6f}, expansion {ex['slope']:.6f}; spread [{env[0]:.6f}, {env[1]:.6f}]; 1/W = {1 / W:.6f}, rel dev {prim['slope'] * W - 1:+.2e}", flush=True)
                check(f"SL[{tag} h={h:.6g}]", "1/W lies inside the spread of the fitted slopes (windows x {joint, collapse, expansion})", ok,
                      f"primary {prim['slope']:.6f}, spread [{env[0]:.6f}, {env[1]:.6f}], 1/W {1 / W:.6f}")
                sl.append((prim["slope"], env, ok))
                pst.append(ps)
            summary[(n, fam)] = (sl, pst)
            if all(p is not None for p in pst):
                p1, p2, p3 = [p * (hs[0] / h) if dn.mut("pstar_h1") else p for p, h in zip(pst, hs)]    # pstar_h1: a first-order seed, p*_h proportional to h
                q_r = math.log2((p1 - p2) / (p2 - p3)) if (p2 - p3) != 0 and (p1 - p2) / (p2 - p3) > 0 else float("nan")
                lines = f"p*_h = {p1:.10e}, {p2:.10e}, {p3:.10e}; Richardson order {q_r:.3f}"
                ok_o = CRIT["order_lo"] <= q_r <= CRIT["order_hi"]
                if fam == "S":                  # the continuum threshold of the symmetric family is p = 0 (the unperturbed static wormhole)
                    o12, o23 = math.log2(abs(p1 / p2)), math.log2(abs(p2 / p3))
                    lines += f"; with p* = 0: orders {o12:.3f}, {o23:.3f}"
                    ok_o &= CRIT["order_lo"] <= o12 <= CRIT["order_hi"] and CRIT["order_lo"] <= o23 <= CRIT["order_hi"]
                if np.isfinite(q_r) and q_r > 0:
                    lines += f"; Richardson p* = {p3 - (p2 - p3) / (2 ** q_r - 1):.10e}"
                print(f"    PS {tag}: {lines}", flush=True)
                check(f"PS[{tag}]", f"p*_h converges at second order (order in [{CRIT['order_lo']}, {CRIT['order_hi']}])", ok_o, lines)
            fin = [x for x in sl if x is not None][-2:]
            if len(fin) == 2:
                rel = abs(fin[0][0] / fin[1][0] - 1)
                check(f"RA[{tag}]", f"primary slopes at the two finest h run agree within {CRIT['agree']:.1%}", rel <= CRIT["agree"],
                      f"{fin[0][0]:.6f}, {fin[1][0]:.6f}; rel diff {rel:.2e}")
            if any(p is None for p in pst) and sum(p is not None for p in pst) == 2:
                q = [p for p in pst if p is not None]
                print(f"    PS {tag}: two h only: p*_h = {q[0]:.10e}, {q[1]:.10e}, ratio {q[0] / q[1]:.4f} (second order: 4)")
            # event-definition sensitivity at the middle h
            jm = jobs.get((n + 2, fam, hs[1]))
            if jm is not None and jm.get("status") == "done" and any(r["phase"].startswith("ed:") for r in jm["points"]) and sl[1] is not None:
                ps = pstar_of(jm["lo"], jm["hi"])
                env = sl[1][1]
                oks = []
                for name in CRIT["ed_events"]:
                    rows = [r for r in jm["points"] if r["phase"] == "ed:" + name and r["T"] is not None]
                    s, se, rms, m = fit_joint(rows, ps, CRIT["windows"][CRIT["ed_window"]])
                    oks.append(verdict_inside(env, s))
                    print(f"    ED {tag} h = {hs[1]:.6g}: event {name}: slope {s:.6f} +- {se:.6f} (n = {m}), primary spread [{env[0]:.6f}, {env[1]:.6f}]")
                check(f"ED[{tag}]", "the slope does not depend on the event definition: each alternative lies inside the primary spread", all(oks))
    print_constraints(jobs, summary)
    print_tables(jobs)
    bad = [c for c, ok in results if ok is False]
    print(f"\n{sum(1 for _, ok in results if ok)} of {len(results)} checks pass; {len(bad)} fail" + (f": {', '.join(bad)}" if bad else ""))
    return jobs, bad


def print_constraints(jobs, summary):
    """Reported, not gated (the coordinator's addition): the null constraints along the two runs that bracket p*_h, over the event cell's past."""
    print("\nNull constraints along the two runs that bracket p*_h (lo = last collapse, hi = last expansion), over the nodes of the event cell's past,")
    print("max |C| and max |C| / max(sum of the absolute values of its terms); change with h against the change of p*_h and of the primary slope")
    for n in (2, 3, 4):
        for fam in FAMS:
            if (n, fam) not in summary:
                continue
            sl, pst = summary[(n, fam)]
            rows, cm = [], []
            for k, h in enumerate(CRIT["H"][n]):
                j = jobs.get((n + 2, fam, h))
                cj = None
                if j is not None and j.get("dir"):
                    cp = os.path.join(raw_dir_of(j), "constraints.json")
                    cj = json.load(open(cp)) if os.path.exists(cp) else None
                if not cj or any("pp" not in cj[x] for x in ("lo", "hi")):
                    rows.append(None)
                    cm.append(None)
                    continue
                m = max(cj[x][c] for x in ("lo", "hi") for c in ("pp", "mm"))
                rel = max(cj[x][c] for x in ("lo", "hi") for c in ("rel_pp", "rel_mm"))
                cm.append(m)
                rows.append(f"    CS D={n + 2} {fam} h = {h:.6g}: " + "; ".join(
                    f"{x}: |C++| {cj[x]['pp']:.3e} (rel {cj[x]['rel_pp']:.2e}), |C--| {cj[x]['mm']:.3e} (rel {cj[x]['rel_mm']:.2e}), event {cj[x]['event']}, repro {cj[x]['repro']:.1e}" for x in ("lo", "hi"))
                    + f"; max {m:.3e}, max rel {rel:.2e}")
            for r in rows:
                if r:
                    print(r)
            for k in (0, 1):
                if cm[k] and cm[k + 1] and pst[k] is not None and pst[k + 1] is not None and sl[k] and sl[k + 1]:
                    print(f"    CS D={n + 2} {fam} h {CRIT['H'][n][k]:.6g} -> {CRIT['H'][n][k + 1]:.6g}: max|C| ratio {cm[k] / cm[k + 1]:.3f} (second order: 4); "
                          f"change of p*_h {pst[k + 1] - pst[k]:+.3e} (relative {pst[k + 1] / pst[k] - 1 if pst[k] else float('nan'):+.3f}); "
                          f"change of the primary slope {sl[k + 1][0] - sl[k][0]:+.2e} (relative {sl[k + 1][0] / sl[k][0] - 1:+.2e})")


def raw_dir_of(j):
    return os.path.join(RAW_ROOT, j["dir"])


RAW_ROOT = ""


def print_tables(jobs):
    print("\nFormation values at the first trapped cell and the stall, p* -/+ 1e-%g (not settled masses; #41's question stays open)" % CRIT["stall_e"])
    print("D family h    side fate     T       event cell offset (cells, l)   r_cell   M_cell   M_hor=r^(D-3)/2 | stall: t     r      l     span  dl/dt")
    for (D, fam, h), j in sorted(jobs.items()):
        if not j.get("stall"):
            continue
        for s in ("-1", "1"):
            r = j["stall"].get(s)
            if not r or "T" not in r:
                continue
            st = r.get("stall") or {}
            ms = f"{r['r_cell']:.5f}  {r['M_cell']:.5f}  {r['M_hor']:.5f}" if r.get("fate") == "collapse" else "n/a (expansion event: r = rexp)"
            print(f"{D} {fam}  {h:.6g} {s:>3} {r.get('fate', '-'):9} {r['T']:7.3f}  {r['event_offset_cells']:6d} {r['event_offset_l']:+.4f}   "
                  f"{ms} | " + (f"{st['t']:6.2f} {st['r']:.5f} {st['l']:+.4f} {st['span']:.2f} {st.get('v', float('nan')):+.4f}" if st else "n/a"))


# --------------------------------------------------------------------------------------------------------------------
def main():
    args = sys.argv[1:]
    if not args:
        sys.exit(__doc__)
    cmd, opts = args[0], {}
    for a in args[1:]:
        if a.startswith("--mutate="):
            dn.MUT = a.split("=", 1)[1]
        elif a.startswith("--only="):
            global ONLY
            ONLY = a.split("=", 1)[1].split(",")
        elif a.startswith("--") and "=" in a:
            k, v = a[2:].split("=", 1)
            opts[k] = v
        elif a.startswith("--"):
            opts[a[2:]] = True
        else:
            sys.exit(f"unknown argument {a}")
    if dn.MUT is not None and dn.MUT not in MUTATIONS and dn.MUT not in checks.MUTATIONS:
        sys.exit(f"unknown mutation {dn.MUT}; known: {', '.join(sorted(set(MUTATIONS) | set(checks.MUTATIONS)))}")
    t0 = time.time()
    import platform
    la = open("/proc/loadavg").read().split()[:3] if platform.system() == "Linux" else ["?"]
    print(f"threshold.py {cmd} {opts}: numba {dn.nb.__version__ if dn.HAVE_NUMBA else 'absent'}, numpy {np.__version__}, python {platform.python_version()}, "
          f"{platform.node()}, load average at start {' '.join(la)}, mutation {dn.MUT}, {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    if cmd == "job":
        n = int(opts["D"]) - 2
        ok = job(n, opts["family"], float(opts["h"]), opts["out"], ed=bool(opts.get("ed")), stall=bool(opts.get("stall")), cons=bool(opts.get("cons")), X=float(opts["X"]) if "X" in opts else None)
        print(f"job done: {ok}; {time.time() - t0:.1f} s")
        sys.exit(0 if ok else 1)
    if cmd == "g6":
        ok = g6(opts["out"])
        print(f"g6 done: {ok}; {time.time() - t0:.1f} s")
        sys.exit(0 if ok else 1)
    if cmd == "planted":
        check_planted()
        bad = [c for c, ok in results if ok is False]
    elif cmd == "analyze":
        jobs, bad = analyze(opts["raw"])
    else:
        sys.exit(f"unknown command {cmd}")
    if dn.MUT:
        want = set(MUTATIONS.get(dn.MUT, []))
        fired = {c.split("[")[0] for c in bad}
        print(f"mutation {dn.MUT}: expected to fire {sorted(want)}; fired {sorted(fired)}; "
              + ("targets all fired" if want <= fired else f"SILENT: {sorted(want - fired)}"))
    print(f"{time.time() - t0:.1f} s")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

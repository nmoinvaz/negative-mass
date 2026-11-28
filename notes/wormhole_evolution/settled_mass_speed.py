#!/usr/bin/env python3
"""Re-read the throat speed and the member of the settled-mass refit with estimators that share nothing with `speed_at` (#314). D = 4, delta = 1e-8, h = 0.005. Stored runs only.

#295 froze a refit y = M / r_th - M_0 = a gamma_1 + b gamma_1^2 + c rdot_th over fifteen members; #305 found it misses a route chosen for leverage (X1 to X3) by +9.0e-5, +2.1e-4, +2.0e-4. The miss equals a
coefficient 10.5 % below c, which a speed read 12 % too large in magnitude on that route would also give. This driver reads the throat speed and the member again, from the same stored runs, by
estimators that call none of `speed_at`, `pair_arrays`, `rdot_local`, `departure_time`, and refits the fifteen with each reading.

Speed estimators (all on the pair's mean r(tau), tau the proper time of the stored locus, evaluated at the departure time of |odd|/r = 1e-3 recomputed here):
  S1   dr/dtau = (dr/dt) / (dtau/dt), both by local linear least squares over |t - t_dep| <= 0.5 (window in coordinate time)
  S1r  dr/dt alone (the variable changes; the refit has its own coefficients)
  S2   local linear least squares of r against tau over |tau - tau_dep| <= 1.0 (a window twice as wide as #267's)
  S3   derivative of a quintic polynomial fitted to r(tau) over tau_dep - 2.0 <= tau <= tau_dep + 0.5
Member readings (the static member's fit to M_+(r), M_-(r) at r = 1.25 .. 4 r_th, stored in member.json; the slices are interpolated in t to t_dep as #266's `member_at` does):
  WA   s_th from M_+ - M_- with weights (r_th / r)^2 (inner radii emphasised), r_th from the locus as before
  WJ   (r_th, s_th) jointly from M_+ and M_-, equal weights
  WJi  (r_th, s_th) jointly from M_+ and M_-, weights (r_th / r)^2

    python3 settled_mass_speed.py analyze [--part=dev] [--mutate=NAME] [--baseline=FILE]
    python3 settled_mass_speed.py planted [--mutate=NAME]

--part=dev drops X1 to X3 from the points and raises if any new reading touches them (the mutation effects of the raw README were computed so, before X1 to X3 were read by any new estimator).
Units: a = 1 is the throat radius of the initial static wormhole, 8 pi G = 1, D = 4 (n = 2). Status words as in AGENTS.md: every number here is numerical evidence.
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import settled_mass_route4 as s4          # noqa: E402  (#305's driver; imports #295's, #283's, #266's, #267's and everything below them unchanged)
s3, r = s4.s3, s4.r
b, m, g, th, tl = s3.b, s3.m, s3.g, s3.th, s3.tl

H = 0.005
SQ2 = math.sqrt(2.0)
MUT = None
PART = "all"
I305 = os.path.join(s3.RAW_ROOT, "I305-mac")
EPS = 1e-3
M0 = s3.M0
P4 = s3.CRIT3["P4"]

CRIT = dict(
    # windows (declared in the raw README)
    s1_half=0.5, s2_half=1.0, s3_lo=2.0, s3_hi=0.5, s3_deg=5,
    # bars
    bar_pred=6e-5, bar_frz=1e-9,
    bar_dep=3e-6, bar_spd0=5e-6, bar_memref=1e-7,
    bar_pspd={"S1": 1e-5, "S1r": 1e-5, "S2": 3e-5, "S3": 1e-5},
    bar_pmem=1e-7, bar_pmem_w=1e-6, bar_ppred=1e-12,
    # the notes' printed residuals of the reference reading (settled_mass_route3.md, settled_mass_route4.md), half a unit of the last digit
    res_N=(-4.77e-07, -2.15e-06, +1.21e-06), bar_resN=1e-8,
    res_X=(+9.00e-05, +2.13e-04, +1.98e-04), bar_resX=6e-8,
    # the verdict map (#305's recommended task 2)
    map_15=0.02, map_X=0.10,
    planted_noise_seed=314,
)
MEM_MODES = ("WA", "WJ", "WJi")
SPD_KEYS = ("S1", "S1r", "S2", "S3")
READINGS = ("REF", "S1", "S1r", "S2", "S3", "WA", "WJ", "WJi")      # (member, speed) = REF for none, one changed
GATING = ("S1", "S1r", "S2", "S3", "WA", "WJ", "WJi")
results = []
BASE = set()


def mut(name):
    return MUT == name


def check(cid, text, ok, detail="", ph=False):
    ok = None if ok is None else bool(ok)
    results.append((cid, ok, ph))
    print(f"[{'SKIP' if ok is None else ('pass' if ok else 'FAIL')}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def finish(t0):
    nfail = sum(1 for _, ok, ph in results if ok is False and not ph)
    nph = sum(1 for _, ok, ph in results if ok is False and ph)
    print(f"\n{len(results)} checks ({sum(1 for _, _, ph in results if ph)} post hoc, gating nothing), {nfail} gating checks failed, {nph} post hoc checks failed; {time.time() - t0:.0f} s")
    if MUT is not None:
        fired = sorted({c for c, ok, ph in results if ok is False})
        new = sorted(set(fired) - BASE)
        new_dev = [c for c in new if c not in XCHECKS]
        new_x = [c for c in new if c in XCHECKS]
        want = sorted(MUTATIONS.get(MUT, []))
        flipped = sorted(c for c in BASE if c in XCHECKS and c not in fired)
        extra = ""
        ok = new_dev == want
        if MUT == "x_speed":
            ok = ok and PART != "dev" and sorted(c for c in BASE if c in XCHECKS) == flipped and len(flipped) > 0
            extra = f"; baseline X-checks that now pass {flipped}"
        print(f"mutation {MUT}: fired {fired}; failing without a mutation (--baseline) {sorted(BASE)}; newly fired {new_dev}, declared {want}; X-checks newly fired (not declared: X1 to X3 are not in the dev partition) {new_x}{extra}: {'MATCH' if ok else 'DIFFERENT'}")
        return 0 if ok else 1
    return 1 if nfail else 0


# --------------------------------------------------------------------------------------------------------------------
# the stored points
# --------------------------------------------------------------------------------------------------------------------
def load_points(part="all"):
    """The 15 fitted points, #295's N1 to N3 and (part != 'dev') #305's X1 to X3 at h = 0.005, from the stored runs: the run (locus, member rows), the plateau M, and the reference reading
    (#266's member_at at eps = 1e-3, #267's speed_at). Each is a dict(name, grp, e, M, mem, sp)."""
    P = []
    Rb = m.load266(H)
    for B in s3.B_266:
        x = Rb[B]
        P.append(dict(name=f"B{B:g}", grp="B", e=x["e"], M=x["plat"], mem=x["mem"], sp=x["speed"]))
    for e_ in s3.E_267:
        base = os.path.join(s3.I267, "runs", m.TAG, f"h{H:g}") + f"/e{e_:g}"
        d = dict(e=e_, h=H, base=base)
        d["s"] = dict(np.load(base + "/chain/series.npz"))
        d["m"] = json.load(open(base + "/member/member.json"))
        d["loc"] = dict(np.load(base + "/member/locus.npz"))
        mem = b.member_at(d, EPS, "m")
        P.append(dict(name=f"e{e_}", grp="D", e=d, M=g.summary_g(d["s"], H)["plat1"], mem=mem, sp=m.speed_at(d["loc"], mem["t_dep"])))
    for nm, spec in s3.POINTS_283:
        e = r.load_point(os.path.join(s3.I283, "runs"), spec, H)
        x = r.read_point(e)
        P.append(dict(name=nm, grp=nm[0], e=e, M=x["M"], mem=x["mem"], sp=x["sp"]))
    for nm, spec, _ in s3.POINTS3:
        e = r.load_point(os.path.join(s4.I295, "runs"), spec, H)
        x = r.read_point(e)
        P.append(dict(name=nm, grp="N", e=e, M=x["M"], mem=x["mem"], sp=x["sp"]))
    if part != "dev":
        for nm, spec, _ in s4.POINTS4:
            e = r.load_point(os.path.join(I305, "runs"), spec, H)
            x = r.read_point(e)
            P.append(dict(name=nm, grp="X", e=e, M=x["M"], mem=x["mem"], sp=x["sp"]))
    return P


def guard(p):
    """The dev partition: X1 to X3 are not read by any new reading."""
    if PART == "dev" and p["grp"] == "X":
        raise RuntimeError("dev partition: X1 to X3 are not read by a new reading")


# --------------------------------------------------------------------------------------------------------------------
# the pair's mean, the departure time, the four speed estimators (own code: none of #267's / #255's / #266's)
# --------------------------------------------------------------------------------------------------------------------
def pair_mean(loc):
    """On the slices both runs have: t, tau, r (the mean of the two runs, the tuned member to second order in the offset) and odd = (r_b - r_a)/2 (the unstable mode's part)."""
    sa, sb = loc["s_a"], loc["s_b"]
    ia = np.nonzero(np.isin(sa, sb))[0]
    ib = np.searchsorted(sb, sa[ia])
    t = 0.5 * (loc["t_a"][ia] + loc["t_b"][ib])
    tau = 0.5 * (loc["tau_a"][ia] + loc["tau_b"][ib])
    if mut("tau_unit"):
        tau = tau * SQ2
    rr = 0.5 * (loc["r_a"][ia] + loc["r_b"][ib])
    odd = 0.5 * (loc["r_b"][ib] - loc["r_a"][ia])
    return t, tau, rr, odd


def dep_own(t, rr, odd, eps=EPS, tmin=6.0):
    """First time (t >= tmin) at which |odd| / r reaches eps, by linear interpolation of ln(|odd| / r) between the two slices that bracket it (the mode grows exponentially)."""
    eps = eps * 10.0 if mut("dep_eps") else eps
    q = np.abs(odd) / rr
    idx = np.nonzero((t >= tmin) & (q >= eps))[0]
    k = int(idx[0])
    if k == 0 or q[k - 1] <= 0.0:
        return float(t[k])
    w = (math.log(eps) - math.log(q[k - 1])) / (math.log(q[k]) - math.log(q[k - 1]))
    return float(t[k - 1] + w * (t[k] - t[k - 1]))


def lsq_poly(x, y, x0, deg):
    """Least-squares polynomial in (x - x0) of degree deg through (x, y); returns the coefficients [c0, c1, ...] (c1 = the derivative at x0)."""
    dx = np.asarray(x, float) - x0
    V = np.vander(dx, deg + 1, increasing=True)
    return np.linalg.lstsq(V, np.asarray(y, float), rcond=None)[0]


def window(x, x0, lo, hi):
    return np.nonzero((x >= x0 - lo) & (x <= x0 + hi))[0]


def speeds_at(t, tau, rr, t0):
    """The four speed readings at t0 from the arrays (t, tau, r): the same code for the stored points and for the planted cases."""
    tau0 = float(np.interp(t0, t, tau))
    # S1: window in coordinate time
    h1 = 1.0 if mut("s1_win") else CRIT["s1_half"]
    i1 = window(t, t0, h1, h1)
    drdt = lsq_poly(t[i1], rr[i1], t0, 1)[1]
    dtdt = lsq_poly(t[i1], tau[i1], t0, 1)[1]
    if mut("s1_nolapse"):
        s1 = drdt
    elif mut("s1_lapse_inv"):
        s1 = drdt * dtdt
    else:
        s1 = drdt / dtdt
    # S1r: dr/dt, or (mutation) the slope against tau
    if mut("s1r_tau"):
        s1r = lsq_poly(tau[i1], rr[i1], tau0, 1)[1]
    else:
        s1r = drdt
    # S2: window twice as wide in tau
    h2 = 3.0 if mut("s2_win3") else CRIT["s2_half"]
    i2 = window(tau, tau0, h2, h2)
    s2 = lsq_poly(tau[i2], rr[i2], tau0, 1)[1]
    # S3: derivative of a quintic fitted to r(tau) over the late stall
    deg = 2 if mut("s3_deg2") else CRIT["s3_deg"]
    i3 = window(tau, tau0, CRIT["s3_lo"], CRIT["s3_hi"])
    s3_ = lsq_poly(tau[i3], rr[i3], tau0, deg)[1]
    # S0: #267's estimator (linear, +-0.5 in tau) written here, evaluated exactly at t0 (control SPD0 only)
    h0 = 1.0 if mut("spd0_win") else 0.5
    i0 = window(tau, tau0, h0, h0)
    s0 = lsq_poly(tau[i0], rr[i0], tau0, 1)[1]
    return dict(S1=float(s1), S1r=float(s1r), S2=float(s2), S3=float(s3_), S0=float(s0), lapse=float(dtdt))


def point_speeds(p):
    """Departure time and the four speed readings of one stored point."""
    guard(p)
    t, tau, rr, odd = pair_mean(p["e"]["loc"])
    t0 = dep_own(t, rr, odd)
    out = speeds_at(t, tau, rr, t0)
    out["t_dep"] = t0
    return out


# --------------------------------------------------------------------------------------------------------------------
# the member readings (weights on the static member's fit)
# --------------------------------------------------------------------------------------------------------------------
_MFIT = {}


def static_pred(rth, s, rr):
    return tl.member_M(2, rth, s, rr)


def fit_row(row, mode, cache=True):
    """(s_th, r_th) of one stored profile row. 'unit': the stored fit's objective with unit weights (control); 'WA': M_+ - M_- with weights (r_th / r)^2, r_th the locus's; 'WJ', 'WJi': (r_th, s_th) jointly
    from M_+ and M_- with unit weights or weights (r_th / r)^2 (the weights use the locus r_th, fixed)."""
    key = (id(row), mode, MUT if MUT in ("wa_outer", "wj_swap", "memref_w") else None)
    if cache and key in _MFIT:
        return _MFIT[key]
    from scipy.optimize import least_squares, minimize_scalar
    rr = np.array(row["rr_m"], float)
    Mp, Mm = np.array(row["Mp_m"], float), np.array(row["Mm_m"], float)
    rth0 = float(row["r_th"])
    if mode in ("unit", "WA"):
        w = np.ones_like(rr) if mode == "unit" else (rth0 / rr) ** 2
        if mode == "unit" and mut("memref_w"):
            w = (rth0 / rr) ** 2
        if mode == "WA" and mut("wa_outer"):
            w = (rr / rth0) ** 2
        dm = Mp - Mm

        def cost(s):
            a, bm = static_pred(rth0, s, rr)
            return float(np.sum(w * ((a - bm) - dm) ** 2))
        o = minimize_scalar(cost, bounds=(-0.95, 0.95), method="bounded", options=dict(xatol=1e-9))
        out = (float(o.x), rth0)
    else:
        w = np.ones_like(rr) if mode == "WJ" else (rth0 / rr) ** 2
        sw = np.sqrt(w)
        if mut("wj_swap"):
            Mp, Mm = Mm, Mp

        def res(x):
            a, bm = static_pred(x[0], x[1], rr)
            return np.concatenate([sw * (a - Mp), sw * (bm - Mm)])
        o = least_squares(res, [rth0, float(row["s_m"])], x_scale=[1e-3, 1e-3], xtol=1e-14, ftol=1e-14, gtol=1e-14)
        out = (float(o.x[1]), float(o.x[0]))
    if cache:
        _MFIT[key] = out
    return out


def member_read(p, mode, t0):
    """(gamma_1, r_th) at t0 by the member reading `mode`, the two bracketing slices interpolated linearly in t as #266's member_at does."""
    guard(p)
    rows = p["e"]["m"]["rows"]
    ts = np.array([x["t"] for x in rows])
    k = min(max(int(np.searchsorted(ts, t0)), 1), len(rows) - 1)
    r0, r1 = rows[k - 1], rows[k]
    w = (t0 - r0["t"]) / (r1["t"] - r0["t"])
    f0, f1 = fit_row(r0, mode), fit_row(r1, mode)
    s = f0[0] + w * (f1[0] - f0[0])
    rt = f0[1] + w * (f1[1] - f0[1])
    return b.gamma_of_s(s), rt


# --------------------------------------------------------------------------------------------------------------------
# readings, the refit, the prediction
# --------------------------------------------------------------------------------------------------------------------
def readings_of(p):
    """The point's value of every reading: dict(name -> (g1, rth, rdot)) for REF, S1..S3 (reference member), WA, WJ, WJi (reference speed), plus the cross terms (S_k with W_j)."""
    sp = point_speeds(p)
    ref = (p["mem"]["g1"], p["mem"]["r_th"], p["sp"])
    out = dict(REF=ref)
    for k in SPD_KEYS:
        out[k] = (ref[0], ref[1], sp[k])
    t0 = p["mem"]["t_dep"]
    for j in MEM_MODES:
        g1, rt = member_read(p, j, t0)
        out[j] = (g1, rt, ref[2])
        for k in ("S1", "S2", "S3"):
            out[f"{k}+{j}"] = (g1, rt, sp[k])
    out["_sp"] = sp
    return out


def arrays15(P15, rd, key):
    """gamma_1, r_th, rdot_th (reading `key`) and M of the points P15."""
    g1 = np.array([x[key][0] for x in rd[:len(P15)]])
    rth = np.array([x[key][1] for x in rd[:len(P15)]])
    rdot = np.array([x[key][2] for x in rd[:len(P15)]])
    return g1, rth, rdot, np.array([p["M"] for p in P15])


def refit(P15, rd, key, form="P4"):
    g1, rth, rdot, M = arrays15(P15, rd, key)
    m0 = 0.223095 if mut("m0_shift") else M0
    return s3.fit(g1, rth, rdot, M, form, m0=m0)


def pred(co, g1, rth, rdot, m0=None):
    m0 = (0.223095 if mut("m0_shift") else M0) if m0 is None else m0
    c = co[2] * (0.8 if mut("coef_c") else 1.0)
    return rth * (m0 + co[0] * g1 + co[1] * g1 ** 2 + c * rdot)


def refit_residuals(P, rd, key):
    """Refit the fifteen with reading `key`; M - predicted at every point of P (the fifteen first, then N, then X): returns (coefficients, residuals array, readings of the points)."""
    P15 = [p for p in P if p["grp"] not in ("N", "X")]
    n15 = len(P15)
    drop = [k for k, p in enumerate(P15) if not (mut("drop_pt") and p["name"] == "e12")]
    Pd = [P15[k] for k in drop]
    rdd = [rd[k] for k in drop]
    co = refit(Pd, rdd, key)
    res = []
    for p, x in zip(P, rd):
        g1_, rth_, rdot_ = x[key]
        if p["grp"] == "N" and mut("n_speed"):
            rdot_ = rdot_ * 1.3
        if p["grp"] == "X" and mut("x_speed"):
            rdot_ = rdot_ * 0.895
        if p["grp"] == "N" and mut("n_g1"):
            g1_ = g1_ * 1.05
        res.append(p["M"] - float(pred(co, g1_, rth_, rdot_)))
    return co, np.array(res)


# --------------------------------------------------------------------------------------------------------------------
# planted cases
# --------------------------------------------------------------------------------------------------------------------
def planted_arrays(p, noise, rng, stress=1.0):
    """A planted pair mean around t_dep with a known speed. The throat radius is the point's own degree-7 model of r(tau) over [tau_dep - 2.5, tau_dep + 1.2] (stress multiplies every coefficient above the
    first by `stress`), tau(t) is the point's cubic model of the lapse, the grid is the stored slices, Gaussian noise `noise` on r. Returns t, tau, r, t_dep, the true dr/dtau and dr/dt there."""
    t, tau, rr, odd = pair_mean(p["e"]["loc"])
    t0 = dep_own(t, rr, odd)
    tau0 = float(np.interp(t0, t, tau))
    sel = window(tau, tau0, 2.5, 1.2)
    cr = np.polyfit(tau[sel] - tau0, rr[sel], 7)[::-1]
    sel_t = window(t, t0, 3.0, 3.0)
    ct = np.polyfit(t[sel_t] - t0, tau[sel_t] - tau0, 3)[::-1]
    tt = t[sel_t]
    tp = tau0 + np.polyval(ct[::-1], tt - t0)
    cs = cr.copy()
    cs[2:] *= stress
    rp = np.polyval(cs[::-1], tp - tau0) + rng.normal(0.0, noise, len(tt))
    return tt, tp, rp, t0, float(cs[1]), float(cs[1] * ct[1])


def measured_noise(P):
    """The scatter of r(tau) about its smooth model: the median over the points of the residual sd of a degree-9 fit over [tau_dep - 2.5, tau_dep + 0.6]."""
    sds = []
    for p in P:
        if p["grp"] == "X":
            continue
        t, tau, rr, odd = pair_mean(p["e"]["loc"])
        t0 = dep_own(t, rr, odd)
        tau0 = float(np.interp(t0, t, tau))
        sel = window(tau, tau0, 2.5, 0.6)
        c = np.polyfit(tau[sel] - tau0, rr[sel], 9)
        sds.append(float(np.std(rr[sel] - np.polyval(c, tau[sel] - tau0))))
    return float(np.median(sds)), float(np.max(sds))


def check_pspd(P, stress_report=True):
    """PSPD: each speed estimator recovers a planted speed through the same code path (speeds_at): the planted r(tau) of each of the points' own models with the measured noise; truth dr/dtau (dr/dt for S1r)."""
    noise, noise_max = measured_noise(P)
    rng = np.random.default_rng(CRIT["planted_noise_seed"])
    errs = {k: [] for k in SPD_KEYS}
    worst = {}
    for p in P:
        if p["grp"] == "X":
            continue
        tt, tp, rp, t0, v, vt = planted_arrays(p, noise, rng)
        sp = speeds_at(tt, tp, rp, t0)
        for k in SPD_KEYS:
            truth = vt if k == "S1r" else v
            errs[k].append(abs(sp[k] - truth))
            if abs(sp[k] - truth) >= max(errs[k]):
                worst[k] = (p["name"], sp[k] - truth)
    print(f"    PSPD planted: noise sd {noise:.2e} (median of the points; largest {noise_max:.2e}), {len(errs['S1'])} planted cases (each point's own model of r(tau), degree 7 over [-2.5, +1.2])")
    for k in SPD_KEYS:
        print(f"      {k:3s}: max |error| {max(errs[k]):.2e} at {worst[k][0]} (signed {worst[k][1]:+.2e}), median {np.median(errs[k]):.2e}, bar {CRIT['bar_pspd'][k]:g}")
        check(f"PSPD_{k}", f"estimator {k} recovers the planted speed of every planted case ({len(errs[k])}) to the bar", max(errs[k]) <= CRIT["bar_pspd"][k],
              f"max error {max(errs[k]):.2e} (bar {CRIT['bar_pspd'][k]:g})")
    if stress_report:
        rng2 = np.random.default_rng(CRIT["planted_noise_seed"] + 1)
        for k in SPD_KEYS:
            errs[k] = []
        for p in P:
            if p["grp"] == "X":
                continue
            tt, tp, rp, t0, v, vt = planted_arrays(p, noise, rng2, stress=2.0)
            sp = speeds_at(tt, tp, rp, t0)
            for k in SPD_KEYS:
                errs[k].append(abs(sp[k] - (vt if k == "S1r" else v)))
        print("    PSPD stress (curvature coefficients of r(tau) x 2, reported, gates nothing): " + ", ".join(f"{k} {max(errs[k]):.2e}" for k in SPD_KEYS))


def check_pmem():
    """PMEM: planted static members (exact tl.member_M profiles of known gamma_1 and r_th at the stored radii) through the same fits. WJ and WJi recover (r_th, s_th) to the bar from a locus r_th that is 3e-4
    off; WA returns the minimizer of its own objective, compared with a grid search (an independent route); unit weights through `fit_row` give the stored objective's minimizer."""
    from scipy.optimize import minimize_scalar          # noqa: F401
    r_pts = np.array(tl.STALL["r_pts"], float)
    e_j, e_a, e_u = [], [], []
    rng = np.random.default_rng(CRIT["planted_noise_seed"] + 2)
    for g1 in (0.003, 0.012, 0.048, 0.075):
        rth_t = 1.0 + 0.005 * rng.standard_normal()
        s_t = g1 / math.sqrt(1 + g1 * g1)
        rth_loc = rth_t + 3e-4
        rr = r_pts * rth_loc
        Mp, Mm = tl.member_M(2, rth_t, s_t, rr)
        row = dict(s=int(1000 * g1) + 7, r_th=rth_loc, rr_m=rr.tolist(), Mp_m=Mp.tolist(), Mm_m=Mm.tolist(), s_m=s_t * (1 + 1e-3))
        for mode in ("WJ", "WJi"):
            s_f, r_f = fit_row(row, mode, cache=False)
            e_j.append(max(abs(s_f - s_t), abs(r_f - rth_t)))
        # WA against a grid search of the same weighted objective (r_th = the locus's: its minimizer differs from s_t by the 3e-4 offset)
        s_a, _ = fit_row(row, "WA", cache=False)
        w = (rth_loc / rr) ** 2
        dm = Mp - Mm
        grid = np.linspace(s_t - 2e-3, s_t + 2e-3, 801)
        cost = np.array([float(np.sum(w * ((lambda ab: ab[0] - ab[1])(tl.member_M(2, rth_loc, s, rr)) - dm) ** 2)) for s in grid])
        k = int(np.argmin(cost))
        c2 = np.polyfit(grid[max(k - 3, 0):k + 4], cost[max(k - 3, 0):k + 4], 2)
        e_a.append(abs(s_a - (-c2[1] / (2 * c2[0]))))
        s_u, _ = fit_row(row, "unit", cache=False)
        cost_u = np.array([float(np.sum(((lambda ab: ab[0] - ab[1])(tl.member_M(2, rth_loc, s, rr)) - dm) ** 2)) for s in grid])
        ku = int(np.argmin(cost_u))
        cu = np.polyfit(grid[max(ku - 3, 0):ku + 4], cost_u[max(ku - 3, 0):ku + 4], 2)
        e_u.append(abs(s_u - (-cu[1] / (2 * cu[0]))))
    print(f"    PMEM planted: joint fits recover (r_th, s_th) to {max(e_j):.1e}; WA against its grid search {max(e_a):.1e}; unit weights against the grid {max(e_u):.1e}")
    check("PMEM", "planted static members: WJ and WJi recover (r_th, s_th) to 1e-7 from a locus r_th 3e-4 off; WA and the unit-weight fit equal an independent grid search of their objective to 1e-6",
          max(e_j) <= CRIT["bar_pmem"] and max(e_a) <= CRIT["bar_pmem_w"] and max(e_u) <= CRIT["bar_pmem_w"], f"joint {max(e_j):.1e} (bar {CRIT['bar_pmem']:g}), WA {max(e_a):.1e}, unit {max(e_u):.1e} (bar {CRIT['bar_pmem_w']:g})")


def check_ppred(P):
    """PPRED: planted data that follow a refit pass the prediction test; with c x 0.8 (the mutation's factor) an X-like row misses by more than the bar; a point displaced by 1e-4 in M is found by the leave-one-out."""
    g1 = np.array([0.003, 0.012, 0.048, 0.075, 0.048, 0.048, 0.048, 0.048, 0.012, 0.048, 0.075, 0.012, 0.048, 0.075, 0.075])
    rth = np.linspace(0.99, 1.02, 15)
    rdot = np.array([-1e-4, -3e-4, 6e-5, 1.6e-3, -6e-4, 5e-4, 8e-4, 1e-3, -5e-4, -6.6e-4, 5e-4, -4e-4, -3e-4, 1e-3, 3e-4])
    co0 = np.array(P4)
    M = rth * (M0 + co0[0] * g1 + co0[1] * g1 ** 2 + co0[2] * rdot)
    co = s3.fit(g1, rth, rdot, M, "P4")
    xrow = (0.048, 0.99, -2.7e-3)
    Mx = xrow[1] * (M0 + co0[0] * xrow[0] + co0[1] * xrow[0] ** 2 + co0[2] * xrow[2])
    pp = pred(co, *xrow)
    miss = abs(Mx - pp)
    co_bad = co * np.array([1.0, 1.0, 0.8])
    miss_bad = abs(Mx - xrow[1] * (M0 + co_bad[0] * xrow[0] + co_bad[1] * xrow[0] ** 2 + co_bad[2] * xrow[2]))
    Md = M.copy()
    Md[5] += 1e-4
    eM, _ = s3.loo(g1, rth, rdot, Md, "P4")
    ok = float(np.max(np.abs(co - co0))) <= 1e-9 and miss <= CRIT["bar_ppred"] and miss_bad > CRIT["bar_pred"] and abs(abs(eM[5]) - 1e-4) <= 3e-5
    check("PPRED", "planted data that follow the refit return its coefficients (1e-9), pass the prediction at an X-like row (1e-12), fail it with c x 0.8 (more than 6e-5), and the leave-one-out finds a point displaced by 1e-4",
          ok, f"coefficients {float(np.max(np.abs(co - co0))):.1e}; planted miss {miss:.1e}; with c x 0.8 {miss_bad:.2e}; displaced point's LOO error {abs(eM[5]):.2e}")


def check_planted(P):
    check_pspd(P)
    check_pmem()
    check_ppred(P)


# --------------------------------------------------------------------------------------------------------------------
# the analysis
# --------------------------------------------------------------------------------------------------------------------
def analyze_cmd():
    t0 = time.time()
    P = load_points(PART)
    names = [p["name"] for p in P]
    n15 = sum(1 for p in P if p["grp"] not in ("N", "X"))
    idxN = [k for k, p in enumerate(P) if p["grp"] == "N"]
    idxX = [k for k, p in enumerate(P) if p["grp"] == "X"]
    print(f"== points: {len(P)} (the fifteen, N1 to N3" + (", X1 to X3" if idxX else "; dev partition: X1 to X3 not loaded") + f") at h = {H}")
    rd = [readings_of(p) for p in P]
    # ---- controls on the points
    dep_err = max(abs(x["_sp"]["t_dep"] - p["mem"]["t_dep"]) for p, x in zip(P, rd))
    check("DEP", "the departure time recomputed here (ln-linear in |odd|/r = 1e-3) equals #266's to 3e-6 at every point", dep_err <= CRIT["bar_dep"], f"max {dep_err:.1e}")
    pa_err = 0.0
    for p in P:
        t, tau, rr, odd = pair_mean(p["e"]["loc"])
        t2, tau2, r2, odd2 = m.pair_arrays(p["e"]["loc"])
        pa_err = max(pa_err, (0.0 if len(t) == len(t2) else 1.0) + float(np.max(np.abs(t - t2))) + float(np.max(np.abs(tau - tau2))) + float(np.max(np.abs(rr - r2))) + float(np.max(np.abs(odd - odd2))))
    check("PAIR", "the pair's mean written here has the slices, t, tau, r and odd part of #267's `pair_arrays` (data equality, to 1e-12)", pa_err <= 1e-12, f"max difference {pa_err:.1e}")
    s0_err = max(abs(x["_sp"]["S0"] - p["sp"]) for p, x in zip(P, rd))
    check("SPD0", "#267's estimator (linear, +-0.5 in tau) written here and evaluated exactly at the departure time equals `speed_at` to 5e-6 at every point", s0_err <= CRIT["bar_spd0"], f"max {s0_err:.1e}")
    mr_err = 0.0
    for p in P:
        if p["grp"] == "X":
            continue
        for row in p["e"]["m"]["rows"][::3]:
            mr_err = max(mr_err, abs(fit_row(row, "unit")[0] - row["s_m"]))
    check("MEMREF", "the member's fit written here with unit weights equals the stored fit (s_th, #225's `fit_asym`) to 1e-7 on every third stored slice", mr_err <= CRIT["bar_memref"], f"max {mr_err:.1e}")
    # ---- reference refit
    co0, res0 = refit_residuals(P, rd, "REF")
    frz_err = float(np.max(np.abs(co0 - np.array(P4))))
    check("FRZ", "the refit over the fifteen read from the stored runs here returns the frozen (a, b, c) of #295 to 1e-9", frz_err <= CRIT["bar_frz"], f"max {frz_err:.1e}: {co0}")
    rN = [res0[k] for k in idxN]
    check("REFN", "the reference reading's residuals M - P4 at N1 to N3 equal #295's printed values to 1e-8", max(abs(a - b_) for a, b_ in zip(rN, CRIT["res_N"])) <= CRIT["bar_resN"], ", ".join(f"{x:+.2e}" for x in rN))
    if idxX:
        rX = [res0[k] for k in idxX]
        check("REFX", "the reference reading's residuals M - P4 at X1 to X3 equal #305's printed values to 6e-8", max(abs(a - b_) for a, b_ in zip(rX, CRIT["res_X"])) <= CRIT["bar_resX"], ", ".join(f"{x:+.2e}" for x in rX))
    # ---- per reading
    ref15 = np.array([rd[k]["REF"][2] for k in range(n15)])
    rms15 = float(np.sqrt(np.mean(ref15 ** 2)))
    print(f"\n== speeds at t_dep (reference = #267's `speed_at`); rms of the reference speeds over the fifteen {rms15:.3e}")
    print("    point    t_dep    lapse    REF          S1           S1r(dr/dt)   S2           S3           | S1-REF   S2-REF   S3-REF")
    for p, x in zip(P, rd):
        sp = x["_sp"]
        print(f"    {p['name']:6s} {sp['t_dep']:8.4f} {sp['lapse']:7.4f}  {p['sp']:+.4e}  {sp['S1']:+.4e}  {sp['S1r']:+.4e}  {sp['S2']:+.4e}  {sp['S3']:+.4e}  | {sp['S1'] - p['sp']:+.1e} {sp['S2'] - p['sp']:+.1e} {sp['S3'] - p['sp']:+.1e}")
    print("\n== the member at t_dep (gamma_1, r_th): reference, then the change under WA, WJ, WJi")
    print("    point    g1 ref      r_th ref    | dg1 WA     dg1 WJ     dg1 WJi    | drth WJ    drth WJi")
    for p, x in zip(P, rd):
        print(f"    {p['name']:6s} {x['REF'][0]:.6f}  {x['REF'][1]:.6f}  | {x['WA'][0] - x['REF'][0]:+.2e}  {x['WJ'][0] - x['REF'][0]:+.2e}  {x['WJi'][0] - x['REF'][0]:+.2e}  | {x['WJ'][1] - x['REF'][1]:+.2e}  {x['WJi'][1] - x['REF'][1]:+.2e}")
    print("\n== refit over the fifteen with each reading: coefficients, in-sample and leave-one-out errors of the fifteen (M - predicted), and M - predicted at N1 to N3" + (" and X1 to X3" if idxX else ""))
    print("    reading   a            b            c            max in-sample  max LOO     |  N1        N2        N3" + ("        |  X1        X2        X3" if idxX else ""))
    table = {}
    for key in ("REF",) + GATING:
        co, res = refit_residuals(P, rd, key)
        g1_, rth_, rdot_, M_ = arrays15(P[:n15], rd, key)
        eM, _ = s3.loo(g1_, rth_, rdot_, M_, "P4")
        table[key] = (co, res)
        line = f"    {key:6s}  {co[0]:+.6f}  {co[1]:+.6f}  {co[2]:+.6f}  {np.max(np.abs(res[:n15])):.2e}      {np.max(np.abs(eM)):.2e}   | " + "  ".join(f"{res[k]:+.2e}" for k in idxN)
        if idxX:
            line += "  | " + "  ".join(f"{res[k]:+.2e}" for k in idxX)
        print(line)
    print("\n    supplementary (declared, gates nothing): a speed estimator with a member reading together")
    for k in ("S1", "S2", "S3"):
        for j in MEM_MODES:
            co, res = refit_residuals(P, rd, f"{k}+{j}")
            line = f"    {k}+{j:4s}  {co[0]:+.6f}  {co[1]:+.6f}  {co[2]:+.6f}  {np.max(np.abs(res[:n15])):.2e}                   | " + "  ".join(f"{res[kk]:+.2e}" for kk in idxN)
            if idxX:
                line += "  | " + "  ".join(f"{res[kk]:+.2e}" for kk in idxX)
            print(line)
    # ---- the gating predictions
    print()
    for key in ("REF",) + GATING:
        co, res = table[key]
        rN = np.array([res[k] for k in idxN])
        check(f"PN_{key}", f"reading {key}: #295's N1 to N3 within 6e-5 of the refit over the fifteen", float(np.max(np.abs(rN))) <= CRIT["bar_pred"], "M - pred " + ", ".join(f"{x:+.2e}" for x in rN))
        if idxX:
            rX = np.array([res[k] for k in idxX])
            check(f"PX_{key}", f"reading {key}: #305's X1 to X3 within 6e-5 of the refit over the fifteen", float(np.max(np.abs(rX))) <= CRIT["bar_pred"], "M - pred " + ", ".join(f"{x:+.2e}" for x in rX))
    # ---- the verdict map for the speed readings (declared in the README before the run)
    if idxX:
        print("\n== verdict map per speed reading (#305's task 2): D15 = max |change of rdot_th| over the fifteen / (0.02 rms), DX = change of rdot_th at X1, X2, X3 relative to the reference")
        for key in SPD_KEYS:
            d15 = max(abs(rd[k][key][2] - rd[k]["REF"][2]) for k in range(n15))
            dX = [abs(rd[k][key][2] - rd[k]["REF"][2]) / abs(rd[k]["REF"][2]) for k in idxX]
            res = table[key][1]
            inbar = float(np.max(np.abs([res[k] for k in idxX]))) <= CRIT["bar_pred"]
            small15 = d15 <= CRIT["map_15"] * rms15
            bigX = all(x >= CRIT["map_X"] for x in dX)
            if inbar and small15:
                v = "reading (the miss belonged to the reading)"
            elif inbar:
                v = "reading, but the fifteen's readings change by 2 % or more"
            elif bigX and small15:
                v = "route-dependent reading"
            else:
                v = "form (neither)"
            print(f"    {key:3s}: D15 {d15:.2e} ({d15 / rms15 * 100:.2f} % of rms; bar 2 %), DX " + ", ".join(f"{x * 100:.2f} %" for x in dX) + f"; X within 6e-5: {inbar}  ->  {v}")
        print("\n== the shift of the speed that would close each miss under the reference refit: -(M - P4) / (r_th c), against each estimator's shift at X1, X2, X3")
        for k in idxX:
            need = -res0[k] / (rd[k]["REF"][1] * P4[2] * (0.8 if mut("coef_c") else 1.0))
            print(f"    {P[k]['name']}: needed {need:+.3e} ({need / rd[k]['REF'][2] * 100:+.1f} % of {rd[k]['REF'][2]:+.3e}); shifts " + ", ".join(f"{key} {rd[k][key][2] - rd[k]['REF'][2]:+.1e}" for key in SPD_KEYS))
    # ---- planted
    print()
    check_planted(P)
    return finish(t0)


MUTATIONS = {       # mutation -> the checks it newly fires on the dev partition (X1 to X3 not read by any new estimator); X-checks (PX_*, REFX) are not in the dev partition and are listed apart
    "s1_win": ['PSPD_S1', 'PSPD_S1r'],
    "s1_nolapse": ['PSPD_S1'],
    "s1_lapse_inv": ['PSPD_S1'],
    "s1r_tau": ['PSPD_S1r'],
    "s2_win3": ['PSPD_S2'],
    "s3_deg2": ['PSPD_S3'],
    "spd0_win": ['SPD0'],
    "tau_unit": ['PAIR', 'SPD0'],
    "dep_eps": ['DEP', 'PN_S1', 'PN_S1r', 'PN_S2', 'PN_S3', 'PSPD_S1', 'PSPD_S1r', 'PSPD_S2', 'SPD0'],
    "m0_shift": ['FRZ', 'PPRED', 'REFN'],
    "drop_pt": ['FRZ', 'REFN'],
    "coef_c": ['PN_REF', 'PN_S1', 'PN_S1r', 'PN_S2', 'PN_S3', 'PN_WA', 'PN_WJ', 'PN_WJi', 'PPRED', 'REFN'],
    "n_speed": ['PN_REF', 'PN_S1', 'PN_S1r', 'PN_S2', 'PN_S3', 'PN_WA', 'PN_WJ', 'PN_WJi', 'REFN'],
    "n_g1": ['PN_REF', 'PN_S1', 'PN_S1r', 'PN_S2', 'PN_S3', 'PN_WA', 'PN_WJ', 'PN_WJi', 'REFN'],
    "wa_outer": ['PMEM'],
    "wj_swap": ['PMEM'],
    "memref_w": ['MEMREF', 'PMEM'],
    "x_speed": [],      # record only (X1 to X3 speeds x 0.895 = #305's post hoc c_eff / c, from its published output): every PX_rho that fails in the baseline is expected to pass; fires nothing
}
XCHECKS = ("REFX",) + tuple(f"PX_{k}" for k in READINGS)


def main():
    global MUT, PART, BASE
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        return 2
    cmd = args[0]
    for a in args[1:]:
        if a.startswith("--mutate="):
            MUT = a.split("=", 1)[1]
        if a.startswith("--part="):
            PART = a.split("=", 1)[1]
        if a.startswith("--baseline="):
            for ln in open(a.split("=", 1)[1]):
                if ln.startswith("[FAIL] "):
                    BASE.add(ln[7:].split(":", 1)[0])
    print(f"settled_mass_speed.py {cmd} {args[1:]}: numpy {np.__version__}, mutation {MUT}, part {PART}, {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    if cmd == "analyze":
        return analyze_cmd()
    if cmd == "planted":
        t0 = time.time()
        P = load_points("dev")
        check_planted(P)
        return finish(t0)
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main())

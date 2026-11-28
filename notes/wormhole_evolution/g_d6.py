#!/usr/bin/env python3
"""Family G at D = 6 at its third resolution (#228): a thin driver around #40's threshold.py, #225's throat_locus.py and #255's throat_drift.py, which are imported unchanged.

    python3 g_d6.py tail --D=6 --family=G --h=0.00125 --out=DIR     constraint_runs and stall_runs of #40's job, after the job (the same two files, written separately)
    python3 g_d6.py threshold-analyze --raw=DIR [--mutate=NAME]      threshold.py analyze with CRIT["not_run"] = {} (the third h of G at D = 6 is run) and the P246 check
    python3 g_d6.py locus series|stall ...                           throat_locus.py main, with RAW40 = the merged record
    python3 g_d6.py drift bisect|stall ...                           throat_drift.py main, with the same RAW40

RAW below is the raw folder I228-dh610; `runs40/` is the merged record: links to #40's jobs and to this task's job at h = 0.00125. The criteria are in RAW/README.md.
"""
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import checks           # noqa: E402,F401
import dn               # noqa: E402
import threshold as th  # noqa: E402
import throat_locus as tl   # noqa: E402
import throat_drift as td   # noqa: E402

RAW = os.path.expanduser(os.environ.get("G_D6_RAW", "~/Source/negative-mass-raw/I228-dh610"))
MERGED = os.path.join(RAW, "runs40")
P246_JSON = os.path.expanduser("~/Source/negative-mass-raw/I246-dh610/runs/pstar_D6_G_k2.json")
SQ2 = math.sqrt(2.0)
LMUT = None         # this file's mutation


def lmut(name):
    return LMUT == name


def merged():
    """Point #225's and #255's code at the merged record (#40's raw plus the job of this task); the third h of D = 6 is read from the job of record. With G_D6_SEED246=1 the p*_h of that job,
    while it is not there yet, is replaced by #246's: used only to seed the bracket of a bisection (never an amplitude of a run of record)."""
    tl.RAW40 = MERGED
    orig = tl.pstar_raw40

    def pstar_raw40(D, fam, hk):
        if (D, fam, hk) == (6, "G", 2) and os.environ.get("G_D6_SEED246") == "1":
            h = th.CRIT["H"][4][2]
            sp = os.path.join(MERGED, f"D6_G_h{h:g}", "summary.json")
            if not os.path.exists(sp):
                s = json.load(open(P246_JSON))
                assert s["h"] == h, (s["h"], h)
                return s["ps"], s
        return orig(D, fam, hk)
    tl.pstar_raw40 = pstar_raw40
    return tl


def cmd_tail(a):
    """#40's job runs `constraint_runs` and `stall_runs` at its end; here they are run on the job's own bracket after it. Same arguments, same files (constraints.json, stall.json)."""
    n, fam, h, out = int(a["D"]) - 2, a["family"], float(a["h"]), a["out"]
    s = json.load(open(os.path.join(out, "summary.json")))
    assert s.get("status") == "done" and s["h"] == h and s["family"] == fam, s
    X = s["X"]
    t0 = time.time()
    if not os.path.exists(os.path.join(out, "constraints.json")):
        json.dump(th.constraint_runs(n, fam, h, X, s["lo"], s["hi"]), open(os.path.join(out, "constraints.json"), "w"), indent=1)
        print(f"constraints.json written, {time.time() - t0:.0f} s", flush=True)
    if not os.path.exists(os.path.join(out, "stall.json")):
        json.dump(th.stall_runs(n, fam, h, X, 0.5 * (s["lo"] + s["hi"])), open(os.path.join(out, "stall.json"), "w"), indent=1)
        print(f"stall.json written, {time.time() - t0:.0f} s", flush=True)
    s["tail_secs"] = round(time.time() - t0, 1)
    json.dump(s, open(os.path.join(out, "summary.json"), "w"), indent=1)
    return True


# --------------------------------------------------------------------------------------------------------------------
# declared settings of the D = 6 analyses (RAW/README.md, part 2)
# --------------------------------------------------------------------------------------------------------------------
C6 = dict(
    hks=(1, 2),                                   # the two finest h of th.CRIT["H"][4]
    p246_tol=1e-14,                               # P246: p*_h of the job against #246's
    cmp=dict(t_stall=6.0, win=(7.0, 11.0)),       # CMP1/2 at #40's earlier stall time (6.24, expansion side), CMP3 over the slices where the pair at 1e-9 is linear (|odd| <= 3e-4)
    dr=dict(win=(7.0, 11.0)),                     # DR: the throat radius of the locus over the same slices (a drifting scale: growth >= tl.CRIT["dr"]["d5_min"])
    ref_ctl={0: 0.990046, 1: 0.994756},           # #40's primary diagonal slope (threshold.md), 6 digits; the finest h's is the job's own, computed from its points
    delta=1e-10, Bs=(0.04, 0.1, 0.2), Bs_extra=(0.02, 0.07), ts=(8.0, 10.0, 12.0), vwins=(0.25, 0.5, 1.0),   # local rdot over +- vwin in tau (0.5 as #255; the others give the window part of its error)
    fit_ok=dict(kappa=0.05, rel_se=0.1),          # the fit of #255 is read only where it converges: kappa >= 0.05 and SE(rdot_inf) <= 10 % of |rdot_inf|
    odd_max=3e-4,                                 # NO6: as #255
    pl1=dict(c=1.2, hs=(0.005, 0.0025, 0.00125), t=(0.8, 2.7), bar=2e-4, rdot=5e-4, order=(1.75, 2.25)),       # #255's PL1, PL2 bars at n = 4, at the production h of D = 6
    pl2=dict(vs=(0.13, -0.3), hs=(0.005, 0.0025, 0.00125), t=(0.8, 2.7), bar=5e-4, rdot=1e-3),
    pl3=dict(cases=((0.026, 0.48, -0.013), (0.05, 0.15, -0.015), (0.0, 0.30, 0.05)), noise=1e-6, rdot_bar=2e-4, nse=3.0),   # synthetic rdot(tau): local slope at ts and the fit's rdot_inf
)
MUTATIONS = {   # name: the checks whose verdict changes against the unmutated analysis (computed in the development runs, RAW/dev/mut, before they were declared in README part 2)
    "p246_off": None, "ctl_shift": None, "odd_scale6": None, "ra_scale6": None,
}
results = []


def check(cid, text, ok, detail=""):
    ok = None if ok is None else bool(ok)
    results.append((cid, ok))
    tag = "SKIP" if ok is None else ("pass" if ok else "FAIL")
    print(f"[{tag}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def fmt(xs, f="{:.2e}"):
    return "[" + ", ".join(f.format(x) for x in xs) + "]"


def declare():
    """Put the D = 6 entries of part 2 into the dictionaries of #225's and #255's code (their D = 4 and 5 entries are untouched)."""
    tl.CRIT["cmp"]["t_stall"][6] = C6["cmp"]["t_stall"]
    tl.CRIT["cmp"]["win"][6] = C6["cmp"]["win"]
    tl.CRIT["dr"]["win"][6] = C6["dr"]["win"]


def all_results():
    return results + tl.results + td.results


# ---- P246 ------------------------------------------------------------------------------------------------------------
def chk_P246(hks):
    """p*_h of the job at each h against #246's stored value (hk = 2: its own bisection; hk = 0, 1: #40's, which it copied). Agree to the bisection tolerance, or the difference is reported."""
    for hk in hks:
        h = th.CRIT["H"][4][hk]
        sj = json.load(open(os.path.join(MERGED, f"D6_G_h{h:g}", "summary.json")))
        s2 = json.load(open(P246_JSON.replace("_k2", f"_k{hk}")))
        pj, p2 = 0.5 * (sj["lo"] + sj["hi"]), s2["ps"] + (1e-12 if lmut("p246_off") else 0.0)
        check(f"P246[h={h:g}]", f"p*_h of this job agrees with #246's within {C6['p246_tol']:g}", abs(pj - p2) <= C6["p246_tol"] and sj.get("stopped") is None,
              f"{pj!r} against {p2!r}, difference {pj - p2:+.2e}; lo difference {sj['lo'] - s2['lo']:+.1e}, hi difference {sj['hi'] - s2['hi']:+.1e}; job bracket {sj['hi'] - sj['lo']:.1e}")


# ---- the locus series and the member at the stall --------------------------------------------------------------------
def job_slope_ref(hk):
    """The primary diagonal slope of the job's own points (threshold.fit_set, as threshold.analyze), for the control CTL at a h that threshold.md does not quote."""
    h = th.CRIT["H"][4][hk]
    j = th.load_jobs(MERGED)[(6, "G", h)]
    ps = th.pstar_of(j["lo"], j["hi"])
    rows = [r for r in j["points"] if r["phase"] == "series" and r["T"] is not None]
    return th.fit_set(rows, ps)[1]["slope"]


def load_locus(base, hks):
    import glob
    S, stalls = {}, {}
    for hk in hks:
        h = th.CRIT["H"][4][hk]
        recs = []
        for f in sorted(glob.glob(os.path.join(base, "series", f"D6_G_h{h:g}", "points_*.jsonl"))):
            recs += [json.loads(l) for l in open(f)]
        if recs:
            S[(6, hk)] = recs
        a, b = (os.path.join(base, "stall", f"stall_D6_h{h:g}_side{sd}.npz") for sd in (-1, 1))
        if os.path.exists(a) and os.path.exists(b):
            stalls[(6, hk)] = (tl.load_run(a), tl.load_run(b))
    return S, stalls


def analyze_locus(S, stalls):
    """#225's analyze_series unchanged at D = 6 (JOB, CTL, FA, ED, RA), the control made non-vacuous (CTLN), the expansion side alone (the collapse side has no locus event for e >= 4.5, hk = 0),
    and #225's member_report with the D = 6 windows."""
    for (D, hk) in S:
        tl.CRIT["ctl"][(6, hk)] = C6["ref_ctl"].get(hk) or job_slope_ref(hk)
    out = tl.analyze_series(S, full=True)
    for (D, hk), recs in sorted(S.items()):
        h = th.CRIT["H"][4][hk]
        W = tl.W_EXACT[4]
        ps = out[(D, hk)]["ps"]
        old = tl.points40(6, "G", hk)
        shift = (lambda a: a * (1 + 1e-15)) if lmut("ctl_shift") else (lambda a: a)
        nin = sum(shift(r["A"]) in old for r in recs)
        check(f"CTLN[D=6 G h={h:g}]", "every amplitude of the locus series and brackets is a point of #40's job at that h (so that CTL compares them and is not vacuous)", nin == len(recs), f"{nin} of {len(recs)}")
        ser = [r for r in recs if r["phase"] == "series"]
        ncol = sum(1 for r in ser if r["side"] < 0 and r["loc"].get("c:0.5") is not None)
        nexp = sum(1 for r in ser if r["side"] > 0 and r["loc"].get("e:2.0") is not None)
        print(f"    LC D=6 h={h:g}: locus events at the primary level: collapse side {ncol} of {sum(1 for r in ser if r['side'] < 0)} points (r_min = 0.5), expansion side {nexp} of {sum(1 for r in ser if r['side'] > 0)} (r_min = 2)")
        key_p = lambda r: (r["loc"].get("e:2.0") or [None, None])[1]       # noqa: E731
        rows = [r for r in tl.series_rows(recs, key_p) if r["side"] > 0]
        sx = tl.fit_summary(rows, ps, W)
        print(f"    SLx D=6 h={h:g}: expansion side alone, locus tau slope {sx['prim']['slope']:.6f} +- {sx['prim']['se']:.6f} (OLS SE), n = {sx['prim']['n']}; spread [{sx['env'][0]:.6f}, {sx['env'][1]:.6f}]; "
              f"1/W = {1 / W:.6f}, rel {sx['rel']:+.3e}")
        out[(D, hk)]["spx"] = sx
    ks = sorted(out)
    if len(ks) >= 2:
        a, b = out[ks[-2]]["spx"]["prim"]["slope"], out[ks[-1]]["spx"]["prim"]["slope"] * (1.01 if lmut("ra_scale6") else 1.0)
        check("LRAx[D=6 G]", f"the expansion-side locus tau slopes at the two finest h agree within {tl.CRIT['RA']:.1%}", abs(a / b - 1) <= tl.CRIT["RA"],
              f"{a:.6f} (h = {out[ks[-2]]['h']:g}), {b:.6f} (h = {out[ks[-1]]['h']:g}); rel {a / b - 1:+.2e}")
    for (D, hk), (La, Lb) in sorted(stalls.items()):
        tl.member_report(6, hk, La, Lb, sp=out.get((D, hk), {}).get("spx"))
    return out


# ---- the drift -------------------------------------------------------------------------------------------------------
def rdot_at(S, t0, win=None):
    """dr/dtau at t0 by a local linear fit of r against tau over +- win in tau (throat_drift.rdot_local's definition), on the pair's mean."""
    win = td.CRIT["vwin"] if win is None else win
    k = int(np.argmin(np.abs(S["t"] - t0)))
    m = np.abs(S["tau"] - S["tau"][k]) <= win
    if m.sum() < 5 or abs(S["t"][k] - t0) > 0.1:
        return float("nan")
    return float(np.polyfit(S["tau"][m], S["r"][m], 1)[0])


def stopped_pstar(base, B, h):
    """p*_h of a bisection that its 2 h limit stopped: from its points (the last collapse and the first expansion), written to pstar_stopped.json by `stopped-pstar`. None if there is none."""
    f = os.path.join(base, "bisect", f"D6_B{B:g}_h{h:g}", "pstar_stopped.json")
    return json.load(open(f)) if os.path.exists(f) else None


def cmd_stopped_pstar(a):
    """From the points of a stopped bisection: lo = the largest A with fate -1, hi = the smallest A with fate +1; writes pstar_stopped.json (pstar = (lo + hi)/2, the bracket)."""
    d = a["dir"]
    pts = [json.loads(l) for l in open(os.path.join(d, "points.jsonl"))]
    lo = max(p["A"] for p in pts if p.get("fate") == -1 and p.get("phase") == "bisect")
    hi = min(p["A"] for p in pts if p.get("fate") == 1 and p.get("phase") == "bisect")
    out = dict(lo=lo, hi=hi, pstar=0.5 * (lo + hi), bracket=hi - lo, runs=sum(1 for p in pts if p.get("phase") == "bisect"), source="points.jsonl of a bisection stopped by its 2 h limit")
    json.dump(out, open(os.path.join(d, "pstar_stopped.json"), "w"), indent=1)
    print(json.dumps(out))
    return True


def load_drift(base, hks, Bs):
    ents = {}
    for B in Bs:
        for hk in hks:
            h = th.CRIT["H"][4][hk]
            pr = td.load_pair(lambda sd: os.path.join(base, "stall_drift", td.stall_name(6, B, h, C6["delta"], sd)))
            if pr is None:
                continue
            bis = os.path.join(base, "bisect")
            if B != td.CRIT["B_ref"] and not os.path.exists(os.path.join(bis, f"D6_B{B:g}_h{h:g}", "summary.json")):
                sp = stopped_pstar(base, B, h)
                if sp is None:
                    continue
                ps = sp["pstar"]
                print(f"    note: B = {B:g}, h = {h:g}: p*_h from a bisection stopped by its 2 h limit, bracket {sp['bracket']:.1e} after {sp['runs']} runs (declared 1e-12); the pair at delta = {C6['delta']:g} straddles it")
            else:
                ps = td.pstar_of(6, B, hk, bis)
            ents[(B, hk)] = td.pair_entry(6, B, hk, pr[0], pr[1], ps)
            ents[(B, hk)]["delta"] = C6["delta"]
    return ents


def drift_fit(S):
    f = td.fit_rdot(S["t"], S["tau"], S["r"], S["odd"])
    return f


def analyze_drift(ents):
    """For each pair: fates, the odd part (NO6), r_th - 1 and the local rdot at CRIT ts, the fit rdot_inf + c exp(-kappa tau) of #255, the quasi-static ghost density. Then the B exponents
    (p*_h, r_th - 1, local rdot, rdot_inf), each with its propagated error and the difference between the two finest h."""
    n = 4
    tab = {}
    for (B, hk), e in sorted(ents.items()):
        S = e["S"]
        h = e["h"]
        ok = e["fate"] == (-1, 1)
        check(f"JOB6[B={B:g} h={h:g}]", "the pair p*_h -/+ delta has fates (collapse, expansion) on #40's events", ok, f"fates {e['fate']}, on the locus {e['fate_loc']}; p*_h {e['pstar']:.12e}")
        w = (S["t"] >= td.CRIT["win"][0]) & (S["t"] <= td.CRIT["win"][1])
        odd = float(np.max(np.abs(S["odd"][w] * (10.0 if lmut("odd_scale6") else 1.0)))) if w.any() else float("inf")
        check(f"NO6[B={B:g} h={h:g}]", f"the pair's odd part is at most {C6['odd_max']:g} for t in {td.CRIT['win']}: the pair's mean is the tuned member", odd <= C6["odd_max"], f"max |odd| {odd:.1e}")
        f = drift_fit(S)
        conv = bool(np.isfinite(f["inf"]) and f["kappa"] >= C6["fit_ok"]["kappa"] and f["se_inf"] <= C6["fit_ok"]["rel_se"] * abs(f["inf"]))
        R = td.Rt(S)
        rdw = {w: [rdot_at(S, t, w) for t in C6["ts"]] for w in C6["vwins"]}
        rd = rdw[td.CRIT["vwin"]]
        d = td.locus_diag(S, n)
        wq = (S["t"] >= 8.0) & (S["t"] <= 12.0)
        dq = float(np.max(np.abs(d["ratio"][wq] - 1))) if wq.any() else float("nan")
        print(f"    DR B={B:g} h={h:g}: r_th - 1 at t = {C6['ts']}: {fmt(R, '{:.5f}')}; local rdot (tau, +-{td.CRIT['vwin']}) {fmt(rd, '{:+.5f}')} (window +-{C6['vwins'][0]}: {fmt(rdw[C6['vwins'][0]], '{:+.5f}')}, "
              f"+-{C6['vwins'][-1]}: {fmt(rdw[C6['vwins'][-1]], '{:+.5f}')}); fit rdot_inf {f['inf']:+.5f} +- {f['se_inf']:.1e}, kappa {f['kappa']:.4f} +- {f['se_kappa']:.1e}, rms {f['rms']:.1e}, "
              f"n {f['n']} (t {f['t0']:.1f} to {f['t1']:.1f}), {'converged' if conv else 'NOT converged (a slow transient, not an approach to a constant speed)'}; max |rho_g/rho_eq - 1| over t in [8, 12] {dq:.1e} (reported)")
        tab[(B, hk)] = dict(R=R, rd=rd, rdw=rdw, fit=f, conv=conv, ps=e["pstar"], h=h)
    print("    B exponents: a ~ B^alpha between the B values at hand; error = the fit's or the finite-difference's relative errors propagated (sigma_alpha = sqrt((s_a/a)^2 + (s_b/b)^2)/ln(B_a/B_b)),")
    print("    and the difference between the two finest h (the resolution part)")
    Bs = sorted({B for B, _ in tab})
    for hk in sorted({k for _, k in tab}):
        for Ba, Bb in zip(Bs[1:], Bs[:-1]):
            if (Ba, hk) not in tab or (Bb, hk) not in tab:
                continue
            A_, B_ = tab[(Ba, hk)], tab[(Bb, hk)]
            lnB = math.log(Ba / Bb)
            al_p = td.expo(A_["ps"], B_["ps"], Ba, Bb)
            al_R = [td.expo(x, y, Ba, Bb) for x, y in zip(A_["R"], B_["R"])]
            al_d = [td.expo(x, y, Ba, Bb) for x, y in zip(A_["rd"], B_["rd"])]
            al_w = {w: [td.expo(x, y, Ba, Bb) for x, y in zip(A_["rdw"][w], B_["rdw"][w])] for w in C6["vwins"]}
            wspread = [max(al_w[w][i] for w in C6["vwins"]) - min(al_w[w][i] for w in C6["vwins"]) for i in range(len(C6["ts"]))]
            fa, fb = A_["fit"], B_["fit"]
            if A_["conv"] and B_["conv"]:
                al_i = td.expo(fa["inf"], fb["inf"], Ba, Bb)
                err_i = math.sqrt((fa["se_inf"] / fa["inf"]) ** 2 + (fb["se_inf"] / fb["inf"]) ** 2) / lnB
                txt_i = f"{al_i:.3f} +- {err_i:.3f}"
            else:
                al_i, err_i, txt_i = float("nan"), float("nan"), "n/a (a fit did not converge)"
            prev = tab.get((Ba, hk - 1)) is not None and tab.get((Bb, hk - 1)) is not None
            res = ""
            if prev:
                A0, B0 = tab[(Ba, hk - 1)], tab[(Bb, hk - 1)]
                d0 = [td.expo(x, y, Ba, Bb) for x, y in zip(A0["rd"], B0["rd"])]
                res = f"; resolution (h = {A0['h']:g} against {A_['h']:g}): local rdot {fmt([a - b for a, b in zip(al_d, d0)], '{:+.3f}')}, p*_h {al_p - td.expo(A0['ps'], B0['ps'], Ba, Bb):+.3f}"
            print(f"    EXP h={A_['h']:g} B {Bb:g} to {Ba:g}: p*_h {al_p:.3f}; r_th - 1 at ts {fmt(al_R, '{:.3f}')}; local rdot at ts {fmt(al_d, '{:.3f}')} (window part {fmt(wspread, '{:.3f}')}); "
                  f"rdot_inf {txt_i}{res}")
            tab[(Ba, hk)].setdefault("exp", {})[Bb] = dict(p=al_p, R=al_R, d=al_d, wspread=wspread, i=al_i, ei=err_i)
    return tab


def posthoc_departure(hk, La, Lb):
    """Post hoc, gates nothing. The deviation of the run at p* + delta from the tuned throat, in the linear regime, is delta/delta_0 times the pair's odd part (delta_0 = 1e-9, the pair of `stall`):
    its time of departure t_dep(delta) is where |odd| = eps delta_0/delta (eps = 0.1 of the throat radius; 0.03 and 0.3 give the sensitivity). Reported against #40's lifetimes T of the job
    at the same h: the offset T - t_dep, the local slopes per two e-folds of dtau_dep/dln(1/delta) (which is 1/lambda = r_th/W) and of T in t, and the lapse factor dt/dtau on the locus."""
    h = th.CRIT["H"][4][hk]
    W = tl.W_EXACT[4]
    com = np.intersect1d(La["s"], Lb["s"])
    ia, ib = np.searchsorted(La["s"], com), np.searchsorted(Lb["s"], com)
    t, tau = La["t"][ia], 0.5 * (La["tau"][ia] + Lb["tau"][ib])
    odd, r = 0.5 * (Lb["r"][ib] - La["r"][ia]), 0.5 * (La["r"][ia] + Lb["r"][ib])
    ok = np.abs(odd) > 0
    lt, tt, ta, rr = np.log(np.abs(odd[ok])), t[ok], tau[ok], r[ok]
    m = np.concatenate([[True], np.diff(lt) > 0]) & (tt < 17.5)               # the monotone, linear part of the growth
    lt, tt, ta, rr = lt[m], tt[m], ta[m], rr[m]
    d0 = 10.0 ** (-tl.STALL["e"])
    pts = [json.loads(l) for l in open(os.path.join(MERGED, f"D6_G_h{h:g}", "points.jsonl"))]
    ser = {(round(p["e"], 1), p["side"]): p["T"] for p in pts if p["phase"] == "series" and p["T"] is not None}
    print(f"    PH D=6 h={h:g} (post hoc, gates nothing): departure of p* + delta from the tuned throat, |odd(t)| = eps delta_0/delta of the pair at delta_0 = {d0:g}; T of #40's job in t")
    for eps in (0.03, 0.1, 0.3):
        rows = []
        for e in np.arange(3.0, 11.01, 0.5):
            target = math.log(eps * d0 * 10.0 ** e)
            if not (lt.min() <= target <= lt.max()):
                continue
            rows.append(dict(e=float(e), t=float(np.interp(target, lt, tt)), tau=float(np.interp(target, lt, ta)), r=float(np.interp(target, lt, rr))))
        for x in rows:
            x["Tm"], x["Tp"] = ser.get((round(x["e"], 1), -1)), ser.get((round(x["e"], 1), 1))
        if eps == 0.1:
            for x in rows:
                print(f"       eps {eps} e {x['e']:4.1f}: t_dep {x['t']:6.2f} tau_dep {x['tau']:6.2f} r_th {x['r']:.3f} r_th/W {x['r'] / W:.3f}; T(collapse) {x['Tm'] if x['Tm'] is None else round(x['Tm'], 2)} "
                      f"(T - t_dep {x['Tm'] - x['t'] if x['Tm'] else float('nan'):+.2f}), T(expansion) {x['Tp'] if x['Tp'] is None else round(x['Tp'], 2)} (T - t_dep {x['Tp'] - x['t'] if x['Tp'] else float('nan'):+.2f})")
        by = {round(x["e"], 1): x for x in rows}
        out = []
        for a, b in ((3.0, 5.0), (4.0, 6.0), (5.0, 7.0), (6.0, 8.0), (7.0, 9.0)):
            if a in by and b in by:
                L = (b - a) * math.log(10.0)
                sl_tau = (by[b]["tau"] - by[a]["tau"]) / L
                sl_t = (by[b]["t"] - by[a]["t"]) / L
                sl_Tm = (by[b]["Tm"] - by[a]["Tm"]) / L if by[a]["Tm"] and by[b]["Tm"] else float("nan")
                rW = 0.5 * (by[a]["r"] + by[b]["r"]) / W
                out.append(f"e {a:g}-{b:g}: dtau_dep/dln {sl_tau:.3f} (r_th/W at the mid {rW:.3f}), dt_dep/dln {sl_t:.3f}, dT(collapse)/dln {sl_Tm:.3f}")
        offs = [x["Tm"] - x["t"] for x in rows if x["Tm"] and x["e"] >= 4.5]
        print(f"       eps {eps}: " + "; ".join(out) + (f"; T(collapse) - t_dep over e >= 4.5: {min(offs):+.2f} to {max(offs):+.2f}" if offs else ""))


def posthoc_B(tab, base):
    """Post hoc, gates nothing: p*_h(B) at the three h (h = 0.005 from dev/bisect; B = 0.1 from #40's raw), the Richardson order and p* and the exponent of the extrapolated p*, and, with the
    energy E_net = (B^2 - p*_h^2) I_2 that the two pulses leave (#255's definition, I_2 from throat_drift.pulse_integrals(4)), the fitted rdot_inf / E_net^(3/2) where the fit converged."""
    I2, _ = td.pulse_integrals(4)

    def pstar(B, hk):
        h = th.CRIT["H"][4][hk]
        if B == td.CRIT["B_ref"]:
            return tl.pstar_raw40(6, "G", hk)[0]
        for b in (os.path.join(base, "bisect"), os.path.join(RAW, "dev", "bisect")):
            f = os.path.join(b, f"D6_B{B:g}_h{h:g}", "summary.json")
            if os.path.exists(f):
                return json.load(open(f))["pstar"]
        return None
    print("    PHB D=6 (post hoc, gates nothing): p*_h(B) at h = 0.005, 0.0025, 0.00125, the observed order, the Richardson p*; E_net = (B^2 - p*_h^2) I_2 (I_2 = %.4f), rdot_inf / E_net^(3/2) where the fit converged" % I2)
    ext = {}
    for B in sorted({B for B, _ in tab} | {0.02, 0.04, 0.07, 0.1, 0.2}):
        p = [pstar(B, k) for k in (0, 1, 2)]
        if all(x is not None for x in p) and (p[1] - p[2]) != 0 and (p[0] - p[1]) / (p[1] - p[2]) > 0:
            q = math.log2((p[0] - p[1]) / (p[1] - p[2]))
            ext[B] = p[2] - (p[1] - p[2]) / (2 ** q - 1)
            print(f"       B {B:g}: p*_h {fmt(p, '{:.7e}')}, order {q:.3f}, Richardson p* {ext[B]:.7e}")
        else:
            print(f"       B {B:g}: p*_h {fmt([x for x in p if x is not None], '{:.7e}')} (not all three h)")
    Bs = sorted(ext)
    print("       exponent of the extrapolated p*: " + "; ".join(f"B {a:g} to {b:g}: {math.log(ext[b] / ext[a]) / math.log(b / a):.3f}" for a, b in zip(Bs, Bs[1:])) + (f"; B {Bs[0]:g} to {Bs[-1]:g}: {math.log(ext[Bs[-1]] / ext[Bs[0]]) / math.log(Bs[-1] / Bs[0]):.3f}" if len(Bs) > 1 else ""))
    for (B, hk), x in sorted(tab.items()):
        f = x["fit"]
        if x["conv"]:
            A = pstar(B, hk) if pstar(B, hk) is not None else x["ps"]
            E = (B * B - A * A) * I2
            print(f"       B {B:g} h {x['h']:g}: E_net {E:.5f}, rdot_inf {f['inf']:+.5f} +- {f['se_inf']:.1e}, rdot_inf / E_net^(3/2) = {f['inf'] / E ** 1.5:.4f}")


def chk_planted6():
    """The estimators of #255 at n = 4: planted static throat of radius c and planted moving throats through the locus code (PL1, PL2 of #255's bars), and planted rdot(tau) series through the
    local slope and the fit (PL3)."""
    n = 4
    for cid, kind, pars, c in (("PL1_6", "scale", (C6["pl1"]["c"],), C6["pl1"]), ("PL2_6", "boost", C6["pl2"]["vs"], C6["pl2"])):
        ok = True
        for par in pars:
            errs = td.planted_errors(kind, n, par, c["hs"], *c["t"])
            fin = errs[-1]
            ords = td.order_of([e["r"] for e in errs]) if kind == "scale" else []
            ok &= fin["r"] <= c["bar"] and fin["rdot"] <= c["rdot"]
            if kind == "scale":
                ok &= all(c["order"][0] <= x <= c["order"][1] for x in ords)
            print(f"    {cid} D = 6, " + (f"static throat of radius {par:g}" if kind == "scale" else f"static throat moving with v = {par:g} (r_th = 1)") + f", h = {c['hs'][-1]:g}: max |r_th - exact| {fin['r']:.1e}, "
                  f"max |rdot| {fin['rdot']:.1e}" + (f"; orders {fmt(ords, '{:.2f}')}" if ords else ""))
        check(cid, ("planted static throat of known radius c = %g (n = 4)" % c["c"]) if kind == "scale" else "planted static throat moving in the grid (v = %s, n = 4): radius 1, rate 0" % (c["vs"],),
              ok, f"bars: |r_th - exact| <= {c['bar']:g}, |rdot| <= {c['rdot']:g}")
    p3 = C6["pl3"]
    rng = np.random.default_rng(228)
    tau = np.arange(0.0, 14.0, 0.0005)
    worst_loc, worst_se = 0.0, 0.0
    for inf, kap, c in p3["cases"]:
        r = 1.0 + inf * tau + (c / kap) * (1.0 - np.exp(-kap * tau)) + p3["noise"] * rng.standard_normal(len(tau))
        tt = tau * (td.SQ2 if lmut("tau_unit") or td.lmut("tau_unit") else 1.0)       # tau_unit: tau in units of sqrt2 t
        S = dict(t=tau, tau=tt, r=r)
        loc = [abs(rdot_at(S, t) - (inf + c * math.exp(-kap * t))) for t in C6["ts"]]
        f = td.fit_rdot(tau, tt, r, np.zeros_like(tau), win=(6.5, 12.0))
        nse = abs(f["inf"] - inf) / f["se_inf"] if f["se_inf"] > 0 else float("inf")
        worst_loc, worst_se = max(worst_loc, max(loc)), max(worst_se, nse)
        print(f"    PL3_6 planted rdot_inf {inf:+.3f}, kappa {kap}, c {c:+.3f}: local rdot error at ts {fmt(loc)}; fit rdot_inf {f['inf']:+.5f} +- {f['se_inf']:.1e} ({nse:.1f} SE), kappa {f['kappa']:.3f}")
    check("PL3_6", "planted rdot(tau) = rdot_inf + c exp(-kappa tau) (three cases incl. a slow one): the local slope at ts to the bar, the fit's rdot_inf within the declared number of its SE",
          worst_loc <= p3["rdot_bar"] and worst_se <= p3["nse"], f"worst local error {worst_loc:.1e} (bar {p3['rdot_bar']:g}), worst fit error {worst_se:.1f} SE (bar {p3['nse']:g})")


def analyze(base, hks, mutation=None):
    results.clear()
    tl.results.clear()
    td.results.clear()
    declare()
    print(f"g_d6.py analyze: base {base}, hks {hks}, mutation {mutation}, numba {dn.nb.__version__ if dn.HAVE_NUMBA else 'absent'}, numpy {np.__version__}, {time.strftime('%Y-%m-%d %H:%M:%S')}")
    chk_P246(hks)
    S, stalls = load_locus(base, hks)
    analyze_locus(S, stalls)
    for (D, hk), (La, Lb) in sorted(stalls.items()):
        posthoc_departure(hk, La, Lb)
    ents = load_drift(base, hks, C6["Bs"] + C6["Bs_extra"])
    tab = analyze_drift(ents)
    posthoc_B(tab, base)
    chk_planted6()
    bad = [c for c, ok in all_results() if ok is False]
    print(f"\n{sum(1 for _, ok in all_results() if ok)} of {len(all_results())} checks pass; {len(bad)} fail" + (f": {', '.join(bad)}" if bad else ""))
    return dict(all_results())


def set_mutation(m):
    """Route a mutation name to the module that knows it."""
    global LMUT
    LMUT = None
    if m is None:
        return
    if m in MUTATIONS:
        LMUT = m
    elif m in tl.MUTATIONS or m in ("time_unit", "tau_diag", "tau_nolapse", "rho_diag", "w_sign", "asym_sign", "pstar_off", "ra_scale", "locus_diag", "drop_bracket", "ed_t", "xc_tau"):
        tl.LMUT = m
    elif m in td.MUTATIONS:
        td.LMUT = m
    else:
        sys.exit(f"unknown mutation {m}")


def main():
    global LMUT
    argv = sys.argv[1:]
    if not argv:
        sys.exit(__doc__)
    cmd, rest = argv[0], argv[1:]
    if cmd == "locus":
        merged()
        sys.argv = ["throat_locus.py"] + rest
        return tl.main()
    if cmd == "drift":
        merged()
        sys.argv = ["throat_drift.py"] + rest
        return td.main()
    a, mut = {}, None
    for x in rest:
        if x.startswith("--mutate="):
            mut = x.split("=", 1)[1]
        elif x.startswith("--") and "=" in x:
            k, v = x[2:].split("=", 1)
            a[k] = v
    if cmd == "tail":
        sys.exit(0 if cmd_tail(a) else 1)
    if cmd == "stopped-pstar":
        sys.exit(0 if cmd_stopped_pstar(a) else 1)
    if cmd == "threshold-analyze":
        th.CRIT["not_run"] = {}
        if mut:
            dn.MUT = mut
        sys.argv = ["threshold.py", "analyze", f"--raw={a['raw']}"] + ([f"--mutate={mut}"] if mut else [])
        return th.main()
    if cmd in ("analyze", "mutate"):
        merged()
        base = a.get("base", os.path.join(RAW, "runs"))
        hks = tuple(int(x) for x in a.get("hks", "1,2").split(","))
        t0 = time.time()
        if cmd == "mutate":
            set_mutation(None)
            base_r = analyze(base, hks)
            print("\n==== baseline done; now the mutation ====\n", flush=True)
            set_mutation(mut)
            mut_r = analyze(base, hks, mut)
            changed = sorted({c.split("[")[0] for c in set(base_r) | set(mut_r) if base_r.get(c) != mut_r.get(c)})
            print(f"\nmutation {mut}: verdicts changed against the unmutated analysis: {changed}; {time.time() - t0:.1f} s")
            return
        set_mutation(mut)
        r = analyze(base, hks, mut)
        print(f"{time.time() - t0:.1f} s")
        sys.exit(1 if any(v is False for v in r.values()) else 0)
    sys.exit(f"unknown command {cmd}")


if __name__ == "__main__":
    main()

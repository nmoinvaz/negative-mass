#!/usr/bin/env python3
"""Settled black-hole mass of family G against the normal amplitude B (#266): #41's apparent-horizon chain (settled_mass.py, through #253's driver
settled_mass_g.py) at D = 4 for B = 0, 0.025, 0.05, 0.1, 0.2 and for family O, each with its own bisected p*_h, and the static member the data stall
near (throat radius r_th, asymmetry gamma_1) at departure, measured by #225's method (throat_locus.py).

A thin driver. It imports settled_mass.py (#41), settled_mass_g.py (#253), threshold.py (#40), throat_locus.py (#225) and dn.py (#39) unchanged.

    python3 settled_mass_b.py pstar  --fam=G --B=0.1 --h=0.01 --out=DIR              bisect p*_h of the member to hi - lo <= 1e-14 (#40's bisection)
    python3 settled_mass_b.py job    --fam=G --B=0.1 --h=0.01 --out=DIR --pstar=DIR   #41's chain at p*_h - delta (delta = 1e-8), the series of M_AH
    python3 settled_mass_b.py member --fam=G --B=0.1 --h=0.01 --out=DIR --pstar=DIR   the pair p*_h -/+ 1e-8 kept whole: locus, profile slices, the member's gamma_1, r_th
    python3 settled_mass_b.py planted [--mutate=NAME]                                  planted cases (member extraction, departure time)
    python3 settled_mass_b.py analyze --raw=DIR [--mutate=NAME]                        every check and table from the stored runs

Units: a = 1 is the throat radius of the initial static wormhole, 8 pi G = 1, D = 4 (n = 2), M = r/2 on the apparent horizon (r_+ = 0), t = (x+ + x-)/sqrt2.
Status words as in AGENTS.md: every mass and every member parameter is numerical evidence.
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dn                      # noqa: E402
import settled_mass as sm      # noqa: E402
import settled_mass_g as g     # noqa: E402  (replaces settled_mass.track by its wrapper, as in #253)
import threshold as th         # noqa: E402
import throat_locus as tl      # noqa: E402

SQ2 = math.sqrt(2.0)
MUT = None
ND = 2
HS = (0.02, 0.01, 0.005)                       # h of #40 for G at D = 4; 0.02 is added to read the order, the two of the issue are 0.01 and 0.005
H2 = (0.01, 0.005)
DELTA = 1e-8
E_PAIR = 8.0                                   # the member pair is p*_h -/+ 10^-8: its collapse side is the run of record
FAMS = (("G", 0.0), ("G", 0.025), ("G", 0.05), ("G", 0.1), ("G", 0.2), ("O", None))
BRACKET = {("O", None): (-0.01, 0.01), ("G", 0.0): (-0.01, 0.01), ("G", 0.025): (-0.01, 0.01), ("G", 0.05): (0.0, 0.02), ("G", 0.1): (0.005, 0.03), ("G", 0.125): (0.01, 0.05), ("G", 0.15): (0.02, 0.07), ("G", 0.2): (0.03, 0.1)}
S_PLATEAU = 0.22309                            # the issue's value of S's settled mass
BAR_S = 6e-4                                   # the issue's bar for B = 0
I40 = os.path.expanduser("~/Source/negative-mass-raw/I40-dh610/runs")
I253 = os.path.expanduser("~/Source/negative-mass-raw/I253-rpi5")
I225 = os.path.expanduser("~/Source/negative-mass-raw/I225-dh610/runs")


def mut(name):
    return MUT == name


def tag_of(fam, B):
    return "O" if fam == "O" else f"G_B{B:g}"


def set_family(fam, B):
    """Select family and normal amplitude in the imported modules (threshold.CRIT['G_B'] is family G's B)."""
    g.FAM = fam
    if fam == "G":
        th.CRIT["G_B"] = B + (0.05 if mut("b_offset") else 0.0)


def hstr(h):
    return repr(float(h))


# --------------------------------------------------------------------------------------------------------------------
# p*_h by bisection (#40's: threshold.bisect on threshold.point, hi - lo <= 1e-14)
# --------------------------------------------------------------------------------------------------------------------
def pstar_cmd(fam, B, h, out, tol=1e-14):
    set_family(fam, B)
    os.makedirs(out, exist_ok=True)
    X = th.CRIT["X"][fam][ND]
    sink = th.Sink(os.path.join(out, "points.jsonl"))
    n = [0]

    def run(p):
        n[0] += 1
        return th.point(ND, fam, h, X, p)
    t0 = time.time()
    lo, hi = BRACKET[(fam, B)]
    lo, hi, stopped = th.bisect(run, lo, hi, tol, sink)
    summ = dict(fam=fam, B=B, h=h, h_arg=hstr(h), X=X, lo=lo, hi=hi, ps=0.5 * (lo + hi), stopped=stopped, runs=n[0], secs=round(time.time() - t0, 1))
    json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print(json.dumps(summ), flush=True)


def load_pstar(pdir, fam, B, h):
    s = json.load(open(os.path.join(pdir, tag_of(fam, B), f"h{h:g}", "pstar", "summary.json")))
    return s


# --------------------------------------------------------------------------------------------------------------------
# the chain
# --------------------------------------------------------------------------------------------------------------------
def job_cmd(fam, B, h, out, pdir, delta=DELTA):
    import resource
    set_family(fam, B)
    ps = load_pstar(pdir, fam, B, h)["ps"]
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    r = g.run_real_g(h, delta, ps, log=lambda m: print(m, flush=True))
    s = r["ser"]
    np.savez(os.path.join(out, "series.npz"), **s)
    X = th.CRIT["X"][fam][ND]
    t6 = {k: th.point(ND, fam, h, X, ps + sg * 1e-6) for k, sg in (("m6", -1), ("p6", 1))}      # the lifetimes at 1e-6: the stall's growth with ln(1/delta)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0 ** 2 if sys.platform == "darwin" else 1024.0)      # bytes on macOS, KB on Linux
    summ = dict(fam=fam, B=B, h=h, delta=delta, ps=ps, A=r["A"], status=r["status"], event=r["event"], point=r["point"], up=r["up"], pt6=t6["m6"], up6=t6["p6"], Np=r["Np"],
                nm0=r["nm0"], F1=r["F1"], K=r["K"], Nc=r["Nc"], levels=len(r["levels"]), cols=int(len(s["r"])), xp0=float(s["xp"][0]) if len(s["r"]) else None,
                xp1=float(s["xp"][-1]) if len(s["r"]) else None, F_end=float(s["F"][-1]) if len(s["r"]) else None, M_start=float(s["M"][0]) if len(s["r"]) else None,
                M_end=float(s["M"][-1]) if len(s["r"]) else None, min_edge=float(s["edge"].min()) if len(s["r"]) else None, secs=round(time.time() - t0, 1), maxrss_MB=round(rss, 1))
    json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print(json.dumps(summ), flush=True)


# --------------------------------------------------------------------------------------------------------------------
# the member at the stall: the pair p*_h -/+ 1e-8 kept whole (#225's `job_stall` at e = 8), profile slices every PROF_DT
# --------------------------------------------------------------------------------------------------------------------
TCAP_MB = 34.0                                 # triangle t <= 34 (#225's TCAP)
PROF_T0, PROF_DT = 3.0, 0.5                    # profile slices t = 3, 3.5, ... up to the end of the triangle


def pair_runs(fam, B, h, ps, e=E_PAIR, tcap=TCAP_MB, log=None, keep_prof=True):
    """The two runs A = p*_h -/+ 10^-e through throat_locus.sweep (dn.row_step unchanged), whole (no early stop on the expansion side beyond #225's)."""
    set_family(fam, B)
    smax = tl.tcap_smax(h, tcap)
    tk = np.arange(PROF_T0, tcap - 3.0 + 1e-9, PROF_DT)
    prof_s = sorted({int(round(t * SQ2 / h)) for t in tk}) if keep_prof else []
    runs = {}
    for sd in (-1, 1):
        A = ps + sd * 10.0 ** (-e)
        t0 = time.time()
        d = tl.make_data(ND, fam, h, smax, A, bare=True)
        L = tl.sweep(d, smax=smax, stop=True, col_levels=tl.LEVELS["col"], exp_levels=tl.LEVELS["exp"], prof_s=prof_s, log=log)
        tl.locus_series(L)
        runs[sd] = L
        if log:
            log(f"pair side {sd}: A = {A!r}, {time.time() - t0:.0f} s, slices {len(L['s'])}, profiles {len(prof_s)}, cross {L['cross']}, ev40 {L['ev40']}")
    return runs


def mean_profile(La, Lb, s_k):
    """Profile of the pair's mean on the slice s_k (the unstable mode is odd in the pair, so it cancels to first order): r, f, r_+, r_- averaged node by node."""
    pa, pb = La["prof"][s_k], Lb["prof"][s_k]
    out = {k: 0.5 * (pa[k] + pb[k]) for k in ("r", "f", "rp", "rm")}
    out["dead"] = pa["dead"] | pb["dead"]
    return out


def gamma_of_s(s):
    return s / math.sqrt(max(1e-12, 1.0 - s * s))


def member_rows(La, Lb, tmax):
    """For each profile slice up to tmax: throat radius of the pair's mean locus, and the member fitted (#225's slice_masses + fit_asym) on La alone (#225's
    literal reading) and on the pair's mean profile. Returns a list of dicts."""
    rows = []
    r_pts = np.array(tl.STALL["r_pts"], float)
    sa, sb = La["s"], Lb["s"]
    for s_k in sorted(La["prof"]):
        h = La["h"]
        t = s_k * h / SQ2
        if t > tmax:
            break
        ka, kb = np.searchsorted(sa, s_k), np.searchsorted(sb, s_k)
        if ka >= len(sa) or sa[ka] != s_k or kb >= len(sb) or sb[kb] != s_k:
            continue
        rth_a, rth_b = float(La["r"][ka]), float(Lb["r"][kb])
        rth = 0.5 * (rth_a + rth_b)
        row = dict(s=int(s_k), t=float(t), r_th=rth, r_th_a=rth_a, r_th_b=rth_b, v=float(0.5 * (La["v"][ka] + Lb["v"][kb])), f=float(0.5 * (La["f"][ka] + Lb["f"][kb])),
                   l=float(0.5 * (La["l"][ka] + Lb["l"][kb])))
        for name, L1, rr in (("a", La, rth_a), ("m", {"prof": {s_k: mean_profile(La, Lb, s_k)}}, rth)):
            Mp, Mm, r_node, ic, M = tl.slice_masses(L1, s_k, ND, r_pts * rr)
            good = np.isfinite(Mp) & np.isfinite(Mm)
            if good.sum() >= 3:
                sth, rms = tl.fit_asym(ND, rr, (r_pts * rr)[good], Mp[good], Mm[good])
                row[f"s_{name}"], row[f"rms_{name}"], row[f"g1_{name}"] = sth, rms, gamma_of_s(sth)
                row[f"Mp_{name}"], row[f"Mm_{name}"], row[f"rr_{name}"] = Mp[good].tolist(), Mm[good].tolist(), (r_pts * rr)[good].tolist()
            else:
                row[f"s_{name}"] = row[f"rms_{name}"] = row[f"g1_{name}"] = float("nan")
        rows.append(row)
    return rows


def odd_series(La, Lb):
    """On the slices both runs have: t, the tuned (mean) throat radius and the odd part (r_b - r_a)/2 = the unstable mode's part (#255)."""
    sa, sb = La["s"], Lb["s"]
    common = np.intersect1d(sa, sb)
    ia, ib = np.searchsorted(sa, common), np.searchsorted(sb, common)
    t = La["t"][ia]
    ra, rb = La["r"][ia], Lb["r"][ib]
    return t, 0.5 * (ra + rb), 0.5 * (rb - ra)


def departure_time(t, rtun, odd, eps, t_min=6.0):
    """First time (linear interpolation between slices, t >= t_min) at which the unstable mode's part reaches eps of the tuned throat radius: |odd|/rtun >= eps."""
    eps = eps * 10.0 if mut("dep_eps") else eps
    q = np.abs(odd) / rtun
    ok = t >= t_min
    idx = np.where(ok & (q >= eps))[0]
    if len(idx) == 0:
        return float("nan")
    k = int(idx[0])
    if k == 0 or not ok[k - 1]:
        return float(t[k])
    w = (eps - q[k - 1]) / (q[k] - q[k - 1])
    return float(t[k - 1] + w * (t[k] - t[k - 1]))


def member_cmd(fam, B, h, out, pdir, e=E_PAIR):
    set_family(fam, B)
    ps = load_pstar(pdir, fam, B, h)["ps"]
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    runs = pair_runs(fam, B, h, ps, e=e, log=lambda m: print(m, flush=True))
    La, Lb = runs[-1], runs[1]
    t, rtun, odd = odd_series(La, Lb)
    t_dep = {str(eps): departure_time(t, rtun, odd, eps) for eps in (1e-3, 1e-2, 1e-1)}
    tmax = max(v for v in t_dep.values() if np.isfinite(v)) + 1.0 if any(np.isfinite(v) for v in t_dep.values()) else 25.0
    rows = member_rows(La, Lb, tmax)
    json.dump(dict(fam=fam, B=B, h=h, e=e, ps=ps, t_dep=t_dep, rows=rows, cross_a=La["cross"], cross_b=Lb["cross"], ev40_a={k: list(v) for k, v in La["ev40"].items()},
                   ev40_b={k: list(v) for k, v in Lb["ev40"].items()}, secs=round(time.time() - t0, 1)), open(os.path.join(out, "member.json"), "w"))
    np.savez_compressed(os.path.join(out, "locus.npz"), t_a=La["t"], r_a=La["r"], l_a=La["l"], f_a=La["f"], v_a=La["v"], tau_a=La["tau"], s_a=La["s"],
                        t_b=Lb["t"], r_b=Lb["r"], l_b=Lb["l"], f_b=Lb["f"], v_b=Lb["v"], tau_b=Lb["tau"], s_b=Lb["s"])
    print(json.dumps(dict(fam=fam, B=B, h=h, t_dep=t_dep, secs=round(time.time() - t0, 1))), flush=True)


# --------------------------------------------------------------------------------------------------------------------
# planted member: an exact static ghost wormhole of known gamma_1 as a profile slice of the sweep, through the same extraction (slice_masses, fit_asym)
# --------------------------------------------------------------------------------------------------------------------
def planted_slice(h, t, g1, nmax=None):
    """Profile slice s = round(t sqrt2 / h) of the exact static member of GGS09a's family with gamma_1 = g1 (tl.ggs_closed: r_th, s_th, m_+, m_-), in the
    double-null gauge x+- = (t -+ ... ): e^{-f} = N^2, r a function of l~ = (x+ - x-)/sqrt2 (N(0) = 1, dl~ = dl/N). The nodes (i, j = s - i) carry r, f and the
    centred differences r_+ = (r[i+1, j] - r[i-1, j])/2h, r_- = (r[i, j+1] - r[i, j-1])/2h that the sweep forms. Returns (s, profile dict, throat radius, s_th)."""
    rth, sth, mp, mm = tl.ggs_closed(g1)
    s = int(round(t * SQ2 / h))
    from scipy.interpolate import CubicSpline
    lmax = (s + 2) * h / SQ2 + 2.0                       # the slice reaches |l~| = t; r <= l~ + r_th, so r_max = lmax + 3 suffices
    side = tl.static_member(ND, rth, sth, rmax=lmax + 3.0)
    tab = {}
    for sg in (1, -1):
        l, r, p, N = side(sg)
        lt = np.concatenate([[0.0], np.cumsum(0.5 * (1.0 / N[1:] + 1.0 / N[:-1]) * np.diff(l))])      # l~(l) on this side
        tab[sg] = (CubicSpline(lt, r), CubicSpline(lt, N))

    def rf(i, j):
        lt_ = (i - j) * h / SQ2
        sg = 1 if lt_ >= 0 else -1
        return float(tab[sg][0](abs(lt_))), float(tab[sg][1](abs(lt_)))
    n = s + 1
    R, F = np.zeros(n), np.zeros(n)
    rp, rm = np.full(n, np.nan), np.full(n, np.nan)
    for i in range(n):
        j = s - i
        r0, N0 = rf(i, j)
        R[i], F[i] = r0, -2.0 * math.log(N0)
        if 1 <= i <= n - 2:
            rp[i] = (rf(i + 1, j)[0] - rf(i - 1, j)[0]) / (2 * h)
            rm[i] = (rf(i, j + 1)[0] - rf(i, j - 1)[0]) / (2 * h)
    prof = dict(r=R, f=F, rp=rp, rm=rm, dead=np.zeros(n, bool))
    return s, prof, rth, sth


def extract_member(L1, s_k, rth):
    """The extraction of member_rows on one profile: M(r) on both sides at the radii of #225's r_pts, the fitted s_th and gamma_1."""
    r_pts = np.array(tl.STALL["r_pts"], float)
    Mp, Mm, r_node, ic, M = tl.slice_masses(L1, s_k, ND, r_pts * rth)
    good = np.isfinite(Mp) & np.isfinite(Mm)
    sth, rms = tl.fit_asym(ND, rth, (r_pts * rth)[good], Mp[good], Mm[good])
    return sth, gamma_of_s(sth), rms, r_pts * rth, Mp, Mm


def planted_member_case(h, t, g1):
    """Planted static member -> slice -> extraction. Returns dict(known, recovered, errors) for gamma_1, throat radius (parabolic minimum of the nodes, as
    the sweep) and the end masses m_+, m_- of the recovered (r_th, s_th) against the closed forms."""
    s, prof, rth, sth = planted_slice(h, t, g1)
    r = prof["r"]
    ic = int(np.argmin(r))
    y0, y1, y2 = r[ic - 1], r[ic], r[ic + 1]
    dk = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
    rth_rec = float(y1 - 0.25 * (y0 - y2) * dk)
    sth_rec, g1_rec, rms, rr, Mp, Mm = extract_member({"prof": {s: prof}}, s, rth_rec)
    mp_rec, mm_rec = tl.asym_masses(ND, rth_rec, sth_rec)
    _, _, mp, mm = tl.ggs_closed(g1)
    return dict(h=h, t=t, g1=g1, g1_rec=g1_rec, rth=rth, rth_rec=rth_rec, sth=sth, sth_rec=sth_rec, rms=rms, mp=mp, mp_rec=mp_rec, mm=mm, mm_rec=mm_rec)


# --------------------------------------------------------------------------------------------------------------------
# checks
# --------------------------------------------------------------------------------------------------------------------
TAGS = (("O", None), ("G", 0.0), ("G", 0.025), ("G", 0.05), ("G", 0.1), ("G", 0.125), ("G", 0.2))       # B = 0.125 is added by this note (supplementary)
STALL_EXPECTED = {("G", 0.2): False}           # declared from the development scan (dev03): the near-threshold run of G has no stall at B = 0.2, a stall at B <= 0.15
CRITB = dict(
    bar_bracket=1e-14, bar_p40=1e-10, order_band=(1.75, 2.25), stall_dT=2.0,
    bar_flux=2e-2, bar_plat=1e-4, bar_flux05=1e-5, bar_hpair=5e-5,
    bar_b10=1e-4, bar_s=BAR_S,
    bar_g0=2e-5, bar_r0=1e-4, eps_dep=1e-3, eps_alt=1e-2,
    bar_c225_s=1e-4, bar_c225_r=3e-4, bar_hst_g=5e-5, bar_hst_r=1e-4,
    bar_pm_g=1e-5, bar_pm_r=1e-6, bar_pm_m=1e-4, bar_dep=2e-5,
    bar_rel=6e-5, rel_B=(0.025, 0.05, 0.1), rel_B_ph=(0.025, 0.05, 0.1, 0.125),
)
C225 = {12.0: (1.0031, 0.0484), 13.0: (1.0028, 0.0483), 14.0: (1.0028, 0.0482)}      # #225's member table (throat_locus.md): t -> r_th, s_th as printed
results = []


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
        want = sorted(MUTATIONS.get(MUT, []))
        print(f"mutation {MUT}: fired {fired}, declared {want}: {'MATCH' if fired == want else 'DIFFERENT'}")
        return 0 if fired == want else 1
    return 1 if nfail else 0


def fmt_list(xs, f="{:.2e}"):
    return "[" + ", ".join(f.format(x) for x in xs) + "]"


def load_all(raw):
    D = {}
    for fam, B in TAGS:
        tg = tag_of(fam, B)
        for h in HS:
            base = os.path.join(raw, "runs", tg, f"h{h:g}")
            e = dict(fam=fam, B=B, tag=tg, h=h, base=base)
            try:
                e["p"] = json.load(open(base + "/pstar/summary.json"))
                e["pts"] = [json.loads(l) for l in open(base + "/pstar/points.jsonl")]
                e["i"] = json.load(open(base + "/chain/summary.json"))
                e["s"] = dict(np.load(base + "/chain/series.npz"))
                e["m"] = json.load(open(base + "/member/member.json"))
                e["loc"] = dict(np.load(base + "/member/locus.npz"))
            except Exception as ex:
                e["err"] = repr(ex)
            D[(tg, h)] = e
    return D


CHAIN_RERUN = {"mass_exponent", "finder_sign", "hm_swap", "interp_shift", "no_zoom", "short_limit", "family_B", "b_offset", "theta_coef", "col_shift"}


def chain_of(e):
    """The stored series and summary, or, under a mutation that acts on the chain, a fresh run from the stored p*_h."""
    if MUT not in CHAIN_RERUN:
        return e["s"], e["i"]
    set_family(e["fam"], e["B"])
    K = 0 if mut("no_zoom") else (12 if mut("short_limit") else None)
    if mut("family_B") and e["fam"] == "G" and e["B"] == 0.1:
        th.CRIT["G_B"] = 0.11
    r = g.run_real_g(e["h"], DELTA, e["p"]["ps"], K=K)
    return r["ser"], dict(status=r["status"], event=r["event"], point=r["point"], up=r["up"], Nc=r["Nc"], F1=r["F1"], K=r["K"], h=e["h"])


def stalls(e):
    i = e["i"]
    ok = i["point"]["fate"] == -1 and i["pt6"]["fate"] == -1 and i["up"]["fate"] == 1 and i["up6"]["fate"] == 1
    dm = i["point"]["T"] - i["pt6"]["T"] if not mut("stl_flat") else 0.0
    dp = i["up"]["T"] - i["up6"]["T"] if not mut("stl_flat") else 0.0
    return ok and dm >= CRITB["stall_dT"] and dp >= CRITB["stall_dT"], dm, dp


def bracket_after(pts, k):
    """lo, hi after the two ends and k bisection steps, replayed from the recorded points (threshold.bisect's rule)."""
    lo, hi = pts[0]["A"], pts[1]["A"]
    for r in pts[2:2 + k]:
        if r["fate"] == 1:
            hi = r["A"]
        elif r["fate"] == -1:
            lo = r["A"]
    return lo, hi


def odd_from_loc(loc):
    sa, sb = loc["s_a"], loc["s_b"]
    common = np.intersect1d(sa, sb)
    ia, ib = np.searchsorted(sa, common), np.searchsorted(sb, common)
    t = loc["t_a"][ia]
    ra, rb = loc["r_a"][ia], loc["r_b"][ib]
    return t, 0.5 * (ra + rb), 0.5 * (rb - ra)


_FIT = {}


def refit(row, name, rth):
    """s_th of the stored profile (row, 'a' collapse side alone, 'm' pair mean) by #225's fit_asym with the given throat radius (so a mutation of the fit acts)."""
    key = (row["s"], row["t"], name, round(rth, 12), tl.LMUT)
    if key not in _FIT:
        rr = np.array(row[f"rr_{name}"])
        if len(rr) < 3:
            _FIT[key] = float("nan")
        else:
            _FIT[key] = tl.fit_asym(ND, rth, rr, np.array(row[f"Mp_{name}"]), np.array(row[f"Mm_{name}"]))[0]
    return _FIT[key]


def member_at(e, eps, kind="m"):
    """Member at the departure time of the pair (|odd|/r_tuned = eps): throat radius r_th, s_th, gamma_1 by linear interpolation in t of the two profile slices that
    bracket it. kind 'm': throat radius and profile of the pair's mean; 'a': the collapse side alone (#225's reading)."""
    if mut("eps_h") and e["h"] == HS[-1]:
        eps = eps * 100.0
    t, rtun, odd = odd_from_loc(e["loc"])
    td = departure_time(t, rtun, odd, eps)
    rows = e["m"]["rows"]
    ts = np.array([r["t"] for r in rows])
    if not np.isfinite(td) or td < ts[0] or td > ts[-1]:
        return dict(t_dep=td, r_th=float("nan"), s_th=float("nan"), g1=float("nan"))
    k = int(np.searchsorted(ts, td))
    k = min(max(k, 1), len(rows) - 1)
    r0, r1 = rows[k - 1], rows[k]
    w = (td - r0["t"]) / (r1["t"] - r0["t"])
    key_r = "r_th_a" if (kind == "a" or mut("rth_collapse_side")) else "r_th"
    rth = [r0[key_r] + 0.0, r1[key_r] + 0.0]
    st = [refit(r0, kind, r0[key_r]), refit(r1, kind, r1[key_r])]
    rt = rth[0] + w * (rth[1] - rth[0])
    sth = st[0] + w * (st[1] - st[0])
    return dict(t_dep=td, r_th=rt, s_th=sth, g1=gamma_of_s(sth))


def s_at(e, t, kind="a"):
    """(r_th, s_th) at the time t by interpolation of the stored slices (for #225's table)."""
    rows = e["m"]["rows"]
    ts = np.array([r["t"] for r in rows])
    k = min(max(int(np.searchsorted(ts, t)), 1), len(rows) - 1)
    r0, r1 = rows[k - 1], rows[k]
    w = (t - r0["t"]) / (r1["t"] - r0["t"])
    key_r = "r_th_a" if kind == "a" else "r_th"
    return r0[key_r] + w * (r1[key_r] - r0[key_r]), refit(r0, kind, r0[key_r]) + w * (refit(r1, kind, r1[key_r]) - refit(r0, kind, r0[key_r]))


def rel_loo(g1, rth, M, c0, perm=False):
    """y = M / r_th - c0 = a g1 + b g1^2; leave one out: fit (a, b) on the others, predict the held-out one."""
    M = np.array(M, float)
    if perm:
        M = np.roll(M, 1)
    y = M / np.array(rth) - c0
    out = []
    for k in range(len(g1)):
        idx = [j for j in range(len(g1)) if j != k]
        A = np.array([[g1[j], g1[j] ** 2] for j in idx])
        ab = np.linalg.lstsq(A, y[idx], rcond=None)[0]
        pred = ab[0] * g1[k] + ab[1] * g1[k] ** 2
        out.append((pred - y[k], ab, y[k]))
    return out


def check_planted():
    """PMEM: exact static members of known gamma_1 as profile slices through the extraction; PDEP: a planted growing pair through the departure time;
    PREL: a planted relation through the leave-one-out test."""
    errs, rows = [], []
    for g1 in (0.0, 0.00293, 0.0483, 0.0755, 0.3):
        for h in H2:
            r = planted_member_case(h, 13.0, g1)
            errs.append((abs(r["g1_rec"] - g1), abs(r["rth_rec"] - r["rth"]), max(abs(r["mp_rec"] - r["mp"]), abs(r["mm_rec"] - r["mm"]))))
            rows.append(f"g1 {g1:g} h {h:g}: err g1 {errs[-1][0]:.1e}, r_th {errs[-1][1]:.1e}, m+- {errs[-1][2]:.1e}")
    for ln in rows:
        print("    PMEM " + ln)
    eg, er, em = (max(x[i] for x in errs) for i in range(3))
    check("PMEM", "planted static members (gamma_1 = 0 to 0.3, known r_th, s_th, end masses m_+, m_-) as profile slices through slice_masses, fit_asym, the parabolic throat minimum and asym_masses, "
          f"h = {H2}", eg <= CRITB["bar_pm_g"] and er <= CRITB["bar_pm_r"] and em <= CRITB["bar_pm_m"],
          f"max error gamma_1 {eg:.1e} (bar {CRITB['bar_pm_g']:g}), r_th {er:.1e} (bar {CRITB['bar_pm_r']:g}), m_+- {em:.1e} (bar {CRITB['bar_pm_m']:g})")
    lam, q0, h_ = 1.18, 1e-9, 0.01 / SQ2
    t = np.arange(0.0, 30.0, h_)
    q = q0 * np.exp(lam * t)
    rt = 1.0 + 0.003 * np.tanh(t / 5.0)
    odd = q * rt
    worst = 0.0
    for eps in (1e-3, 1e-2, 1e-1):
        tstar = math.log(eps / q0) / lam
        worst = max(worst, abs(departure_time(t, rt, odd, eps) - tstar))
    check("PDEP", "planted growing pair (odd part q0 e^{lambda t}, known crossing times of 1e-3, 1e-2, 1e-1) through departure_time", worst <= CRITB["bar_dep"], f"max |t_dep - exact| {worst:.1e} (bar {CRITB['bar_dep']:g})")
    Bs = np.array(CRITB["rel_B"])
    g1 = 4.8 * Bs ** 2
    rth = 1.0 + 0.1 * Bs - 3.0 * Bs ** 2
    c0, a_, b_ = 0.22309, 0.206, 0.8
    M = rth * (c0 + a_ * g1 + b_ * g1 ** 2)
    lo = rel_loo(g1, rth, M, c0, perm=mut("rel_perm"))
    check("PREL", "planted relation M / r_th = M_0 + a gamma_1 + b gamma_1^2 through the leave-one-out test (every held-out point predicted within the bar)",
          max(abs(x[0]) for x in lo) <= CRITB["bar_rel"], f"max |prediction error| {max(abs(x[0]) for x in lo):.1e} (bar {CRITB['bar_rel']:g})")


# --------------------------------------------------------------------------------------------------------------------
# mutations: name -> the checks that must fire (computed in the development runs, declared in the raw README before the runs of record)
# --------------------------------------------------------------------------------------------------------------------
MUTATIONS = {      # name: the checks that fire (computed on the development set, raw I266-mac dev/mut_b, before being declared in its README; REL fails without a mutation, so it is in most sets)
    "crit_dm": ['CFG', 'REL'],   # dm = 2.0 in the chain constants
    "bisect_early": ['BR', 'BR-ph', 'REL'],   # the bracket of the bisection after 30 steps instead of the last
    "pstar_h_swap": ['P40', 'P40-ph', 'REL'],   # #40's p*_h of the next h
    "go_h": ['GO', 'GO-ph', 'REL'],   # O at the next h against G(B = 0)
    "pstar_h1": ['PO-ph', 'REL'],   # p*_h scaled by h_1/h (a first-order seed)
    "stl_flat": ['REL', 'STL', 'STL-ph'],   # the lifetimes at 1e-6 and 1e-8 taken equal
    "evt_side": ['EVT', 'EVT-ph', 'REL'],   # the pair's expansion-side event in place of the collapse-side one
    "flux_sign": ['RC4', 'RC4-ph', 'REL'],   # the horizon flux with the wrong sign
    "short_limit": ['RB10', 'RB10-ph', 'RC5', 'RC5-ph', 'RC8', 'RC8-ph'],   # K = 12: the chain stops at e^-f = 6144
    "no_zoom": ['RB0', 'RB0-ph', 'RB10', 'RB10-ph', 'RC4', 'RC4-ph', 'RC5', 'RC5-ph', 'RC6-ph', 'RC8', 'RC8-ph', 'REL', 'REL-ph', 'RH', 'RH-ph'],   # K = 0: a uniform grid
    "mass_exponent": ['RB0', 'RB0-ph', 'RB10', 'RB10-ph', 'RC4', 'RC4-ph', 'REL-ph'],   # M = r^n/2 in place of r^(n-1)/2
    "family_B": ['EVT', 'EVT-ph', 'RB10', 'RB10-ph', 'RC5', 'RC5-ph', 'RC6-ph', 'REL', 'RH'],   # B = 0.11 in place of 0.1 (the chain refuses: fates at p*_h -/+ delta)
    "b_offset": ['EVT', 'EVT-ph', 'RB0', 'RB0-ph', 'RB10', 'RB10-ph', 'RC5', 'RC5-ph', 'RC6-ph', 'REL', 'RH', 'RH-ph'],   # B = 0.05 in place of 0 (the chain refuses)
    "h_plat": ['REL', 'RH'],   # the plateau at h = 0.02 in place of 0.005 in RH
    "rth_collapse_side": ['MB0', 'MB0-ph', 'REL', 'REL-ph'],   # r_th of the collapse side alone in place of the pair's mean
    "asym_sign": ['C225', 'PMEM', 'REL'],   # #225's fit_asym with the sides swapped
    "dep_eps": ['MB0', 'MB0-ph', 'PDEP', 'REL', 'REL-ph'],   # epsilon x 10 in departure_time
    "eps_h": ['HST', 'HST-ph', 'MB0', 'REL', 'REL-ph'],   # epsilon x 100 at the finest h only
    "rel_perm": ['PREL', 'REL', 'REL-ph'],   # M of the members rolled by one
}


# --------------------------------------------------------------------------------------------------------------------
# the analysis
# --------------------------------------------------------------------------------------------------------------------
G_B_DEFAULT = th.CRIT["G_B"]


def appl(fam, B):
    return STALL_EXPECTED.get((fam, B), True)


def stallscan_cmd(h, Bs, out):
    """Post hoc scan of the stall against B: p*_h (bisected to hi - lo <= 1e-12) and the lifetime at p*_h -/+ 1e-6 and 1e-8 for each B (family G, D = 4)."""
    os.makedirs(out, exist_ok=True)
    X = th.CRIT["X"]["G"][ND]
    rows = []
    for B in Bs:
        set_family("G", B)
        t0 = time.time()
        run = lambda p: th.point(ND, "G", h, X, p)
        lo, hi, stopped = th.bisect(run, 0.02, 0.1, 1e-12)
        ps = 0.5 * (lo + hi)
        T = {}
        for e in (6, 8):
            for sd, nm in ((-1, "m"), (1, "p")):
                r = th.point(ND, "G", h, X, ps + sd * 10.0 ** (-e))
                T[f"{nm}{e}"] = dict(T=r["T"], fate=r["fate"])
        gm = T["m8"]["T"] - T["m6"]["T"]
        gp = T["p8"]["T"] - T["p6"]["T"]
        rows.append(dict(B=B, h=h, ps=ps, hi_lo=hi - lo, stopped=stopped, T=T, gain_m=gm, gain_p=gp, secs=round(time.time() - t0, 1)))
        print(f"B = {B:g}: p* {ps:.10e} (hi - lo {hi - lo:.1e}, {stopped}); T(-1e-6) {T['m6']['T']:.3f} T(-1e-8) {T['m8']['T']:.3f} gain {gm:+.2f}; T(+1e-6) {T['p6']['T']:.3f} T(+1e-8) {T['p8']['T']:.3f} gain {gp:+.2f}", flush=True)
    json.dump(rows, open(os.path.join(out, f"stallscan_h{h:g}.json"), "w"), indent=1)


def stallscan_report(raw):
    """Prints the stored stall scans (every runs/stallscan*/stallscan_h*.json and runs/stallscan_h*.json, merged per h and sorted by B) and the first B at which the stall
    (gain >= 2.0 per side) is gone, per h. A report: it gates nothing."""
    import glob
    byh = {}
    for f in sorted(glob.glob(os.path.join(raw, "runs", "stallscan_h*.json")) + glob.glob(os.path.join(raw, "runs", "stallscan*", "stallscan_h*.json"))):
        for r in json.load(open(f)):
            byh.setdefault(r["h"], {})[r["B"]] = r
    for h in sorted(byh, reverse=True):
        rows = [byh[h][b] for b in sorted(byh[h])]
        print(f"\n== Post hoc: the stall against B at h = {h:g} (T gain from 1e-6 to 1e-8 on each side; a stall gains about 3.9 to 5; prompt gains 0)")
        first = None
        for r in rows:
            sl = r["gain_m"] >= CRITB["stall_dT"] and r["gain_p"] >= CRITB["stall_dT"]
            if not sl and first is None:
                first = r["B"]
            print(f"    B {r['B']:.4f}: p* {r['ps']:.8e} (hi - lo {r['hi_lo']:.1e}), collapse side T {r['T']['m6']['T']:.3f} -> {r['T']['m8']['T']:.3f} (+{r['gain_m']:.2f}), expansion side {r['T']['p6']['T']:.3f} -> {r['T']['p8']['T']:.3f} (+{r['gain_p']:.2f}): {'stall' if sl else 'prompt'}")
        if first is None:
            print("    the stall persists over the whole scan")
        else:
            prev = [r["B"] for r in rows if r["B"] < first]
            print(f"    first B with no stall (gain < {CRITB['stall_dT']}): {first:g}" + (f"; the last with a stall: {max(prev):g}" if prev else "; no scanned B below it"))


PH_TAGS = {("G", 0.125)}          # post hoc by the coordinator's ruling (2026-11-27): supplementary, gates nothing
PH_H = (0.02,)


def is_ph(fam, B, h):
    return (fam, B) in PH_TAGS or h in PH_H


class Acc:
    """Accumulates one check over (member, h): the issue's members at the issue's h gate; B = 0.125 and h = 0.02 are post hoc and gate nothing."""

    def __init__(self):
        self.g, self.p, self.ng, self.np = True, True, 0, 0

    def add(self, good, fam, B, h):
        if is_ph(fam, B, h):
            self.p, self.np = self.p and bool(good), self.np + 1
        else:
            self.g, self.ng = self.g and bool(good), self.ng + 1

    def add_ph(self, good):
        self.p, self.np = self.p and bool(good), self.np + 1


def emit(cid, text, acc, detail=""):
    if acc.ng:
        check(cid, text, acc.g, detail)
    if acc.np:
        check(cid + "-ph", "POST HOC, gates nothing (B = 0.125 and/or h = 0.02): " + text, acc.p, detail, ph=True)


def analyze_cmd(raw):
    t0 = time.time()
    D = load_all(raw)
    miss = {k: e["err"] for k, e in D.items() if "err" in e}
    print(f"== Runs read from {raw}/runs: {len(D) - len(miss)} of {len(D)} (tag, h) complete" + (f"; missing: {miss}" if miss else ""))
    tags = [(fam, B, tag_of(fam, B)) for fam, B in TAGS]
    ok_tags = [(fam, B, tg) for fam, B, tg in tags if all("err" not in D[(tg, h)] for h in HS)]
    st_tags = [(fam, B, tg) for fam, B, tg in ok_tags if appl(fam, B)]
    bad = {k: (sm.CRIT[k], v) for k, v in g.DECLARED.items() if k != "X" and sm.CRIT[k] != v}
    xg, xo = th.CRIT["X"]["G"][ND], th.CRIT["X"]["O"][ND]
    check("CFG", "settled_mass.CRIT (F1, K, W1, rcut, tol, dxp, dm) equal #41's declared chain, X = 28 for G and O, threshold's G_B = 0.1 at import",
          not bad and xg == g.DECLARED["X"] and xo == g.DECLARED["X"] and G_B_DEFAULT == 0.1, f"differs: {bad}; X(G) {xg}, X(O) {xo}, G_B {G_B_DEFAULT}" if bad or xg != 28.0 or xo != 28.0 or G_B_DEFAULT != 0.1 else "")
    acc = Acc()
    for (tg, h), e in D.items():
        acc.add("err" not in e, e["fam"], e["B"], h)
    emit("RUNS", "every (member, h) has its pstar, chain and member output", acc, f"missing {sorted(miss)}" if miss else "")
    print("\n== Planted cases")
    check_planted()

    print("\n== p*_h of each member: bisection to hi - lo <= 1e-14 (threshold.bisect on threshold.point, X = 28), the h of the job as passed")
    acc, det = Acc(), []
    for fam, B, tg in ok_tags:
        for h in HS:
            e = D[(tg, h)]
            p, pts = e["p"], e["pts"]
            lo, hi = bracket_after(pts, 30) if mut("bisect_early") else (p["lo"], p["hi"])
            good = (hi - lo <= CRITB["bar_bracket"] and pts[0]["fate"] == -1 and pts[1]["fate"] == 1 and p["h_arg"] == repr(float(h)) and float(p["h_arg"]) == h
                    and p["stopped"] is None and bracket_after(pts, len(pts)) == (p["lo"], p["hi"]))
            acc.add(good, fam, B, h)
            det.append(f"{tg} h {h:g}: hi - lo {hi - lo:.1e}, {p['runs']} runs, h_arg '{p['h_arg']}', stopped {p['stopped']}{'' if good else '  <-- NOT OK'}")
    for ln in det:
        print("    " + ln)
    emit("BR", "bracket ends are collapse and expansion, the replayed bisection equals lo and hi, hi - lo <= 1e-14, not stopped by 'static through X', h passed as its float repr", acc)
    acc, det = Acc(), []
    for k, h in enumerate(HS):
        e = D.get(("G_B0.1", h))
        if e is None or "err" in e:
            acc.add(False, "G", 0.1, h)
            continue
        ref = g.i40_g(HS[(k + 1) % len(HS)] if mut("pstar_h_swap") else h)
        d = abs(e["p"]["ps"] - ref["ps"])
        acc.add(d <= CRITB["bar_p40"] and float(ref["h_log"]) == h, "G", 0.1, h)
        det.append(f"h {h:g}: p*_h {e['p']['ps']:.16e}, #40's {ref['ps']:.16e} (job.log h '{ref['h_log']}'), difference {d:.1e}")
    for ln in det:
        print("    " + ln)
    emit("P40", "B = 0.1: this machine's p*_h equals #40's (D4_G_h*, dh610) within 1e-10, the h of #40's job.log equal to the h passed", acc, "bar 1e-10")
    acc, det = Acc(), []
    for k, h in enumerate(HS):
        eg, eo = D.get(("G_B0", h)), D.get(("O", HS[(k + 1) % len(HS)] if mut("go_h") else h))
        if eg is None or eo is None or "err" in eg or "err" in eo:
            acc.add(False, "G", 0.0, h)
            continue
        acc.add(eg["p"]["ps"] == eo["p"]["ps"] and eg["i"]["event"] == eo["i"]["event"], "G", 0.0, h)
        det.append(f"h {h:g}: p*_h G(B = 0) {eg['p']['ps']:.16e}, O {eo['p']['ps']:.16e}, events {eg['i']['event']} {eo['i']['event']}")
    for ln in det:
        print("    " + ln)
    emit("GO", "family G with B = 0 and family O have the same p*_h (bit for bit) and the same first trapped cell at every h", acc)
    acc, det = Acc(), []
    for fam, B, tg in ok_tags:
        ps = [D[(tg, h)]["p"]["ps"] * (HS[0] / h if mut("pstar_h1") else 1.0) for h in HS]
        o = math.log2((ps[0] - ps[1]) / (ps[1] - ps[2]))
        gated = appl(fam, B)
        if gated:
            acc.add_ph(CRITB["order_band"][0] <= o <= CRITB["order_band"][1])
        rich = ps[2] - (ps[1] - ps[2]) / 3.0
        det.append(f"{tg}: p*_h = {ps[0]:.6e}, {ps[1]:.6e}, {ps[2]:.6e}; observed order {o:.3f}; Richardson (second order) p* = {rich:.6e}" + ("" if gated else "  (prompt threshold: reported, not checked)"))
    for ln in det:
        print("    " + ln)
    emit("PO", f"p*_h converges at second order at h = {HS}: observed order in {CRITB['order_band']} for every stalling member (it needs three h, and h = 0.02 is post hoc)", acc)

    print("\n== Does the near-threshold run stall? Lifetime T at p*_h -/+ 1e-8 against 1e-6 (the stall adds about ln(100)/W = 3.9 per side), h = 0.02, 0.01, 0.005")
    acc = Acc()
    for fam, B, tg in ok_tags:
        for h in HS:
            e = D[(tg, h)]
            sl, dm, dp = stalls(e)
            want = appl(fam, B)
            acc.add(sl == want, fam, B, h)
            i = e["i"]
            print(f"    {tg:9s} h {h:5g}: collapse side T {i['pt6']['T']:.3f} (1e-6) -> {i['point']['T']:.3f} (1e-8), +{dm:.2f}; expansion side {i['up6']['T']:.3f} -> {i['up']['T']:.3f}, +{dp:.2f}: "
                  f"{'stall' if sl else 'prompt'} (declared {'stall' if want else 'prompt'})")
    emit("STL", "each member's classification (T gains >= 2.0 per side from 1e-6 to 1e-8: stall; else prompt) equals the declared one (stall for B <= 0.125 and O, prompt at B = 0.2)", acc)
    acc, det = Acc(), []
    for fam, B, tg in ok_tags:
        for h in HS:
            e = D[(tg, h)]
            try:
                s_, i_ = chain_of(e)
            except Exception as ex:
                acc.add(False, fam, B, h)
                det.append(f"{tg} h {h:g}: the chain run failed: {ex!r}")
                continue
            ev = e["m"]["ev40_b"]["expand"] if mut("evt_side") and "expand" in e["m"]["ev40_b"] else e["m"]["ev40_a"].get("collapse")
            good = list(i_["event"]) == list(ev) and i_["point"]["fate"] == -1 and i_["up"]["fate"] == 1
            if tg == "G_B0.1":
                good &= list(i_["event"]) == list(g.i40_g(h)["event"])
            acc.add(good, fam, B, h)
            det.append(f"{tg} h {h:g}: chain event {list(i_['event'])}, pair's collapse-side event {ev}" + (f", #40's {g.i40_g(h)['event']}" if tg == "G_B0.1" else ""))
    for ln in det:
        print("    " + ln)
    emit("EVT", "the chain's first trapped cell equals the pair's collapse-side one (same amplitude p*_h - 1e-8), fates collapse at p*_h - 1e-8 and expansion at +1e-8, and #40's for B = 0.1", acc)

    print("\n== The chain (#41's, unchanged) at delta = 1e-8: status, limit, first law, second order, settled value")
    S, I = {}, {}
    for fam, B, tg in st_tags:
        for h in HS:
            try:
                S[(tg, h)], I[(tg, h)] = chain_of(D[(tg, h)])
                if len(S[(tg, h)]["r"]) == 0:
                    raise RuntimeError("empty series: " + str(I[(tg, h)]["status"]))
            except Exception as ex:
                print(f"    run {tg} h = {h:g} failed: {ex!r}")
    for fam, B, tg in ok_tags:
        if not appl(fam, B):
            for h in HS:
                i = D[(tg, h)]["i"]
                print(f"    {tg} h {h:g}: {i['cols']} columns, status '{i['status']}'; M_AH at the first trapped cell {i['M_start']:.5f} (r {2 * i['M_start']:.4f}), F_end {i['F_end']:.3f}  (no stall: not run through the checks below)")
    Fmax = sm.CRIT["F1"] * 2.0 ** sm.CRIT["K"]
    acc, det = Acc(), []
    for fam, B, tg in st_tags:
        for h in HS:
            if (tg, h) not in S:
                acc.add(False, fam, B, h)
                continue
            s, i = S[(tg, h)], I[(tg, h)]
            fine = np.allclose(np.diff(s["xp"]), h, atol=1e-9)
            good = i["status"] == "limit" and 0.5 * Fmax <= s["F"][-1] < Fmax and s["edge"].min() >= i["Nc"] / 8 and fine
            acc.add(good, fam, B, h)
            det.append(f"{tg} h {h:g}: {i['status']}, F_end {s['F'][-1]:.0f}, edge {s['edge'].min():.0f}/{i['Nc']}")
    emit("RC5", f"every run of a stalling member ends at the declared limit (F = e^-f on the horizon, {0.5 * Fmax:g} <= F_end < {Fmax:g}), horizon at least Nc/8 cells from the box edges, no gaps",
         acc, f"{len(det)} runs; worst F_end {min(float(x.split('F_end ')[1].split(',')[0]) for x in det):.4g}" if det else "")
    for ln in det:
        print("    " + ln)
    acc, parts = Acc(), []
    for (tg, h), s in S.items():
        i1 = int(np.argmax(s["F"] >= 20.0))
        dMs = s["M"][i1] - s["M"][0]
        with np.errstate(all="ignore"):
            rel = abs(dMs - g.flux_integral_g(s, 0, i1)) / abs(dMs)
        rel = rel if np.isfinite(rel) else float("inf")
        fam_, B_ = D[(tg, h)]["fam"], D[(tg, h)]["B"]
        acc.add(rel <= CRITB["bar_flux"], fam_, B_, h)
        parts.append(f"{tg} h {h:g}: dM {dMs:+.4e}, flux {g.flux_integral_g(s, 0, i1):+.4e}, relative {rel:.1e}")
    if not S:
        acc.add(False, "G", 0.1, 0.01)
    emit("RC4", "first law on the horizon, both fields: M(x+_b) - M(x+_a) = integral of (r^2/2) e^f r_- (phi0_+^2 - phi1_+^2), from the first trapped cell to F = 20", acc, f"bar {CRITB['bar_flux']:g}")
    for ln in parts:
        print("    " + ln)
    acc = Acc()
    for fam, B, tg in st_tags:
        if not all((tg, h) in S for h in HS):
            acc.add_ph(False)
            continue
        rows = sm.band_orders({h: S[(tg, h)] for h in HS}, HS, g.WINDOWS, DELTA)
        for a_, b_, Ee, o in rows:
            good = all(g.BAND_R[0] <= x <= g.BAND_R[1] for x in o)
            acc.add_ph(good)
            print(f"    RC6 {tg}, x+ - x+_f in [{a_:g}, {b_:g}]: E12 {Ee[0]:.2e}, E23 {Ee[1]:.2e}, order {o[0]:.2f} {'ok' if good else 'OUT'}")
    emit("RC6", f"M_AH(x+) converges at second order at the three h {HS}, three windows, every stalling member (observed order in {g.BAND_R}; it needs three h, and h = 0.02 is post hoc)", acc)
    SS = {k: g.summary_g(s, k[1]) for k, s in S.items()}
    acc = Acc()
    for (tg, h), q in SS.items():
        acc.add(abs(q["dplat"]) <= CRITB["bar_plat"] and q["flux05"] <= CRITB["bar_flux05"], D[(tg, h)]["fam"], D[(tg, h)]["B"], h)
    if not SS:
        acc.add(False, "G", 0.1, 0.01)
    emit("RC8", "M_AH has settled at the limit: plateau mean over the last 1.0 differs from the previous 1.0 by <= 1e-4 and the mean horizon flux (both fields) over the last 0.5 is <= 1e-5, every run", acc,
         "plateau change " + ", ".join(f"{q['dplat']:+.1e}" for q in SS.values()) + "; flux " + ", ".join(f"{q['flux05']:.1e}" for q in SS.values()))
    acc, dets = Acc(), []
    for fam, B, tg in st_tags:
        if not all((tg, h) in SS for h in HS):
            acc.add(False, fam, B, 0.01)
            continue
        h_a = HS[0] if mut("h_plat") else HS[2]
        d = abs(SS[(tg, HS[1])]["plat1"] - SS[(tg, h_a)]["plat1"])
        acc.add(d <= CRITB["bar_hpair"], fam, B, HS[1])
        dets.append(f"{tg}: |plat(0.01) - plat({h_a:g})| = {d:.1e}")
    emit("RH", "the settled value at two h: |plat(0.01) - plat(0.005)| (plateau mean over the last 1.0 of x+) <= 5e-5 for every stalling member", acc, "; ".join(dets))

    print("\n== M(B) with the error bars of this note, kept separate (plateau mean over the last 1.0 of x+, delta = 1e-8)")
    print("    member     plat(0.02)  plat(0.01)  plat(0.005)  Richardson   h: |0.01 - 0.005|   scatter (#253, 6.6e-6 at h = 0.005, 2.5e-5 at 0.01)   sum")
    for fam, B, tg in st_tags:
        if all((tg, h) in SS for h in HS):
            p1, p2, p3 = (SS[(tg, h)]["plat1"] for h in HS)
            bh = abs(p2 - p3)
            print(f"    {tg:9s}  {p1:.6f}   {p2:.6f}   {p3:.6f}    {p3 - (p2 - p3) / 3:.6f}    {bh:.2e}          6.6e-06                                        {bh + 6.6e-6:.2e}")
    print("    delta dependence (physical, not an error): at B = 0.1 the plateau is +5.7e-4 from delta = 1e-6 to 1e-8 (#253); not measured at the other B")
    print("\n== B = 0 against family S's settled mass, and B = 0.1 against #253")
    acc, det = Acc(), []
    for tg, B_ in (("G_B0", 0.0), ("O", None)):
        for h in HS:
            fam_ = "O" if tg == "O" else "G"
            if (tg, h) not in SS:
                acc.add(False, fam_, B_, h)
                continue
            d = SS[(tg, h)]["plat1"] - S_PLATEAU
            acc.add(abs(d) <= CRITB["bar_s"], fam_, B_, h)
            det.append(f"{tg} h {h:g}: plat {SS[(tg, h)]['plat1']:.6f}, minus {S_PLATEAU} = {d:+.2e}")
    for ln in det:
        print("    " + ln)
    emit("RB0", f"B = 0 (G with a zero normal pulse, and O) reproduces S's 0.22309 within {CRITB['bar_s']:g}: plateau means at h = {HS}", acc, "bar 6e-4 (the issue's)")
    acc, det = Acc(), []
    for h in HS:
        if ("G_B0.1", h) not in SS:
            acc.add(False, "G", 0.1, h)
            continue
        ref = g.summary_g(dict(np.load(f"{I253}/runs/h{h:g}_d1e-08/series.npz")), h)
        d = SS[("G_B0.1", h)]["plat1"] - ref["plat1"]
        acc.add(abs(d) <= CRITB["bar_b10"], "G", 0.1, h)
        det.append(f"h {h:g}: plat {SS[('G_B0.1', h)]['plat1']:.6f}, #253's {ref['plat1']:.6f}, difference {d:+.2e}")
    for ln in det:
        print("    " + ln)
    emit("RB10", "B = 0.1 reproduces #253's plateau means (I253-rpi5, delta = 1e-8) within 1e-4 (the roundoff-realization scatter of #253 is 2.5e-5 at h = 0.01)", acc, "bar 1e-4")
    s41 = dict(np.load(f"{g.I41RAW}/runs/h0.005_d1e-08/series.npz"))
    if ("G_B0", 0.005) in S:
        s0 = S[("G_B0", 0.005)]
        print(f"    M_AH(x+ - x+_f), h = 0.005, delta = 1e-8: G(B = 0) x+_f {s0['xp'][0]:.4f}, S (I41-rpi5) x+_f {s41['xp'][0]:.4f}")
        print("    dx+   M_G(B=0)  M_S       difference")
        for dv in (0, 1, 2, 3, 4, 5, 6, 8, 9.4):
            i = int(np.argmin(np.abs(s0["xp"] - (s0["xp"][0] + dv))))
            j = int(np.argmin(np.abs(s41["xp"] - (s41["xp"][0] + dv))))
            print(f"    {dv:4.1f}  {s0['M'][i]:.6f}  {s41['M'][j]:.6f}  {s0['M'][i] - s41['M'][j]:+.2e}")

    print("\n== The static member at departure (#225's slice_masses and fit_asym; throat radius r_th and gamma_1 = s_th / sqrt(1 - s_th^2); t_dep: the unstable mode's part reaches 1e-3 of r_th)")
    MEM = {}
    for fam, B, tg in st_tags:
        for h in HS:
            e = D[(tg, h)]
            MEM[(tg, h)] = dict(m3=member_at(e, CRITB["eps_dep"], "m"), m2=member_at(e, CRITB["eps_alt"], "m"), a3=member_at(e, CRITB["eps_dep"], "a"))
            q = MEM[(tg, h)]
            print(f"    {tg:9s} h {h:5g}: t_dep {q['m3']['t_dep']:.3f}: r_th {q['m3']['r_th']:.5f}, s_th {q['m3']['s_th']:+.5f}, gamma_1 {q['m3']['g1']:+.5f}   | 1e-2: t_dep {q['m2']['t_dep']:.3f}, r_th {q['m2']['r_th']:.5f}, "
                  f"gamma_1 {q['m2']['g1']:+.5f}   | collapse side alone, 1e-3: r_th {q['a3']['r_th']:.5f}, gamma_1 {q['a3']['g1']:+.5f}")
    acc, det = Acc(), []
    for tg, B_ in (("G_B0", 0.0), ("O", None)):
        for h in HS:
            fam_ = "O" if tg == "O" else "G"
            q = MEM.get((tg, h))
            if q is None:
                acc.add(False, fam_, B_, h)
                continue
            gg, rr = q["m3"]["g1"], q["m3"]["r_th"]
            acc.add(np.isfinite(gg) and abs(gg) <= CRITB["bar_g0"] and abs(rr - 1.0) <= CRITB["bar_r0"], fam_, B_, h)
            det.append(f"{tg} h {h:g}: gamma_1 {gg:+.1e}, r_th - 1 {rr - 1.0:+.1e}")
    emit("MB0", "B = 0 (G and O): the member at departure is the symmetric one, |gamma_1| <= 2e-5 and |r_th - 1| <= 1e-4 (planted by symmetry, proved in #213)", acc, "; ".join(det))
    acc, det = Acc(), []
    for h in H2:
        e = D.get(("G_B0.1", h))
        if e is None or "err" in e:
            acc.add(False, "G", 0.1, h)
            continue
        for t_, (rp, sp) in C225.items():
            r_, s_ = s_at(e, t_, "m")
            acc.add(abs(r_ - rp) <= CRITB["bar_c225_r"] and abs(s_ - sp) <= CRITB["bar_c225_s"], "G", 0.1, h)
            det.append(f"h {h:g} t = {t_:g}: r_th {r_:.5f} (#225 {rp}), s_th {s_:.5f} (#225 {sp})")
    for ln in det:
        print("    " + ln)
    emit("C225", "B = 0.1: r_th and s_th of the pair's mean at t = 12, 13, 14 equal #225's printed member table (stall pair at 1e-9) within 3e-4 and 1e-4", acc, "bars 3e-4 (r_th), 1e-4 (s_th)")
    acc, det = Acc(), []
    for fam, B, tg in st_tags:
        a, b = MEM[(tg, HS[1])]["m3"], MEM[(tg, HS[2])]["m3"]
        dg, dr = abs(a["g1"] - b["g1"]), abs(a["r_th"] - b["r_th"])
        acc.add(dg <= CRITB["bar_hst_g"] and dr <= CRITB["bar_hst_r"], fam, B, HS[1])
        det.append(f"{tg}: |gamma_1(0.01) - gamma_1(0.005)| {dg:.1e}, |r_th| {dr:.1e}")
    for ln in det:
        print("    " + ln)
    emit("HST", "the member at departure at the two finest h agree: |gamma_1| <= 5e-5, |r_th| <= 1e-4, every stalling member", acc, "bars 5e-5, 1e-4")

    print("\n== M(B) against gamma_1(B) and r_th(B) at departure (h = 0.005; M = plateau mean over the last 1.0 of x+; dM' = M / r_th - M(B = 0))")
    fin = HS[2]
    rel_ok = all((f"G_B{b}", fin) in SS for b in CRITB["rel_B_ph"]) and ("G_B0", fin) in SS

    def lexp(x1, x0, B1, B0):
        """local exponent d ln|x| / d ln B between two B (nan when a value is not positive: a mutated run)."""
        return math.log(abs(x1) / abs(x0)) / math.log(B1 / B0) if x1 * x0 > 0 else float("nan")
    if rel_ok:
        c0 = SS[("G_B0", fin)]["plat1"]
        rows = []
        for fam, B, tg in st_tags:
            q = MEM[(tg, fin)]["m3"]
            M = SS[(tg, fin)]["plat1"]
            rows.append((B, tg, q["g1"], q["r_th"], M, M / q["r_th"] - c0))
        print("    B        p*_h(0.005)      gamma_1    r_th      M(0.005)   M(0.01)    M(0.02)    M/r_th     dM'        dM'/gamma_1  dM'/gamma_1^2  M - c0")
        for B, tg, gg, rr, M, dMp in rows:
            ps = D[(tg, fin)]["p"]["ps"]
            m1, m2 = SS[(tg, HS[1])]["plat1"], SS[(tg, HS[0])]["plat1"]
            print(f"    {B if B is not None else 'O':<8} {ps:.8e}  {gg:+.5f}  {rr:.5f}  {M:.6f}  {m1:.6f}  {m2:.6f}  {M / rr:.6f}  {dMp:+.5f}   " + (f"{dMp / gg:+.4f}      {dMp / gg ** 2:+.3f}" if abs(gg) > 1e-3 else "-            -") + f"       {M - c0:+.5f}")
        pos = [r for r in rows if r[0] is not None and r[0] > 0 and abs(r[2]) > 1e-3]
        for a_, b_ in zip(pos[:-1], pos[1:]):
            print(f"    local exponents B {a_[0]:g} -> {b_[0]:g}: gamma_1 {lexp(b_[2], a_[2], b_[0], a_[0]):.2f}, dM' {lexp(b_[5], a_[5], b_[0], a_[0]):.2f}, M - c0 "
                  f"{lexp(b_[4] - c0, a_[4] - c0, b_[0], a_[0]):.2f}, r_th - 1 {lexp(b_[3] - 1, a_[3] - 1, b_[0], a_[0]):.2f}")
        print("    the member's asymptotic masses (asym_masses at the fitted r_th, s_th; a report): B, m_+, m_-, (M - M_0)/m_+, dM'/m_+")
        for B, tg, gg, rr, M, dMp in rows:
            if B is None or B <= 0:
                continue
            try:
                mp_, mm_ = tl.asym_masses(ND, MEM[(tg, fin)]["m3"]["r_th"], MEM[(tg, fin)]["m3"]["s_th"])
                print(f"    {B:<6g} m_+ {mp_:+.5f}  m_- {mm_:+.5f}  (M - M_0)/m_+ {(M - c0) / mp_:.4f}  dM'/m_+ {dMp / mp_:.4f}")
            except Exception as ex:
                print(f"    {B:<6g} asym_masses failed: {ex!r}")
        for name, Bset, gate in (("REL", CRITB["rel_B"], True), ("REL-ph", CRITB["rel_B_ph"], False)):
            use = [r for r in rows if r[0] in Bset]
            lo = rel_loo([r[2] for r in use], [r[3] for r in use], [r[4] for r in use], c0, perm=mut("rel_perm"))
            for r, (err, ab, y) in zip(use, lo):
                print(f"    {name}: leave out B = {r[0]:g}: (a, b) = ({ab[0]:.4f}, {ab[1]:.4f}) from the others; dM' measured {y:+.5f}, predicted {y + err:+.5f}, error {err:+.2e}")
            lo_u = rel_loo([r[2] for r in use], [1.0] * len(use), [r[4] for r in use], c0, perm=False)
            print(f"    {name} without the r_th scaling (M - M_0 = a gamma_1 + b gamma_1^2, a report): leave-one-out errors " + ", ".join(f"{x[0]:+.1e}" for x in lo_u) + "; with it: " + ", ".join(f"{x[0]:+.1e}" for x in lo))
            ok_ = max(abs(x[0]) for x in lo) <= CRITB["bar_rel"]
            txt = (f"M / r_th - M_0 = a gamma_1 + b gamma_1^2, leave one out over B = {', '.join(f'{b:g}' for b in Bset)} (M_0 = B = 0's plateau at h = 0.005): every held-out point within 6e-5")
            det_ = f"max error {max(abs(x[0]) for x in lo):.1e} (bar {CRITB['bar_rel']:g}); the form was chosen after the h = 0.02, 0.01 development runs"
            if gate:
                check("REL", txt, ok_, det_)
            else:
                check("REL-ph", "POST HOC, gates nothing (B = 0.125): " + txt, ok_, det_, ph=True)
        use2 = [r for r in rows if r[0] in CRITB["rel_B_ph"]]
        try:
            lo2 = rel_loo([MEM[(r[1], fin)]["m2"]["g1"] for r in use2], [MEM[(r[1], fin)]["m2"]["r_th"] for r in use2], [r[4] for r in use2], c0)
            print("    the same relation with the member read at epsilon = 1e-2 instead of 1e-3 (a report, four B): leave-one-out errors " + ", ".join(f"{x[0]:+.1e}" for x in lo2)
                  + f"; (a, b) from the first three: ({lo2[0][1][0]:.4f}, {lo2[0][1][1]:.4f})")
        except (np.linalg.LinAlgError, ValueError) as ex:
            print(f"    the relation with the member read at epsilon = 1e-2: not computable ({ex!r}; a mutated run)")
        print("    *the relation is a one-parameter family in B: it constrains the shape of M(gamma_1, r_th) along the curve, and does not test the dependence on gamma_1 and r_th separately*")
    else:
        check("REL", "M / r_th against gamma_1, leave one out", False, "a run is missing")
    return finish(t0)


def main():
    global MUT
    argv = sys.argv[1:]
    cmd = argv[0] if argv else "help"
    opts = dict(a.lstrip("-").split("=", 1) if "=" in a else (a.lstrip("-"), True) for a in argv[1:])
    if "mutate" in opts:
        MUT = opts["mutate"]
        if MUT not in MUTATIONS:
            raise SystemExit(f"unknown mutation {MUT}; known: {', '.join(sorted(MUTATIONS))}")
        g.MUT = MUT
        sm.MUT = MUT
        dn.MUT = MUT if MUT in ("theta_coef", "ghost_sign", "scalar_coef") else None
        tl.LMUT = MUT if MUT in ("asym_sign",) else None
        if MUT == "crit_dm":
            sm.CRIT["dm"] = 2.0
    import platform
    import numba as nb
    print(f"settled_mass_b.py {cmd} {opts}: numba {nb.__version__}, numpy {np.__version__}, python {platform.python_version()}, {platform.node()}, "
          f"mutation {MUT}, {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    fam = opts.get("fam", "G")
    B = None if (fam == "O" or cmd == "stallscan") else float(opts.get("B", 0.1))
    if cmd == "pstar":
        pstar_cmd(fam, B, float(opts["h"]), opts["out"])
    elif cmd == "job":
        job_cmd(fam, B, float(opts["h"]), opts["out"], opts["pstar"], delta=float(opts.get("delta", DELTA)))
    elif cmd == "member":
        member_cmd(fam, B, float(opts["h"]), opts["out"], opts["pstar"])
    elif cmd == "planted":
        t0 = time.time()
        check_planted()
        sys.exit(finish(t0))
    elif cmd == "stallscan":
        stallscan_cmd(float(opts["h"]), [float(x) for x in opts["B"].split(",")], opts["out"])
    elif cmd == "analyze":
        rc = analyze_cmd(opts["raw"])
        stallscan_report(opts["raw"])
        sys.exit(rc)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()

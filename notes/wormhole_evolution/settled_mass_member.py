#!/usr/bin/env python3
"""Does the settled mass of family G follow the static member at departure as delta -> 0 (#267)? D = 4, B = 0.1.

#41's apparent-horizon chain (settled_mass.py, through #253's driver settled_mass_g.py) at A = p*_h - delta for delta = 10^-e, and the static member the data
stall near (throat radius r_th, asymmetry gamma_1) at departure, measured by #225's method (throat_locus.py) on the pair p*_h -/+ delta, through #266's driver
settled_mass_b.py. The settled mass is set against the member by #266's relation (the declared prediction) and by a fit over the delta scan; the drift of the
throat radius along the stall is fitted by #255's estimator (throat_drift.py).

A thin driver. It imports settled_mass_b.py (#266), which imports settled_mass.py (#41), settled_mass_g.py (#253), threshold.py (#40), throat_locus.py (#225) and
dn.py (#39), and throat_drift.py (#255), all unchanged.

    python3 settled_mass_member.py pstar   --h=0.01 --out=DIR                          bisect p*_h of G (B = 0.1) to hi - lo <= 1e-14 (#266's `pstar`)
    python3 settled_mass_member.py job     --h=0.01 --e=8 --out=DIR --pstar=DIR        #41's chain at p*_h - 10^-e: the series of M_AH
    python3 settled_mass_member.py member  --h=0.01 --e=8 --out=DIR --pstar=DIR        the pair p*_h -/+ 10^-e kept whole: locus, profile slices, gamma_1, r_th along the stall
    python3 settled_mass_member.py planted [--mutate=NAME]                             planted cases (the delta scan through the analysis)
    python3 settled_mass_member.py analyze --raw=DIR [--mutate=NAME]                   every check and table from the stored runs

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
import settled_mass_b as b         # noqa: E402  (#266's driver; imports dn, settled_mass, settled_mass_g, threshold, throat_locus unchanged)
import throat_drift as td          # noqa: E402  (#255's estimator of the throat radius's speed)
dn, sm, g, th, tl = b.dn, b.sm, b.g, b.th, b.tl

SQ2 = math.sqrt(2.0)
MUT = None
ND = 2
FAM, B0 = "G", 0.1
TAG = "G_B0.1"
H2 = (0.01, 0.005)                             # the issue's two h
H3 = (0.02, 0.01, 0.005)                       # 0.02 is added for the order in h (post hoc)
E_GATE = (7, 8, 9, 10, 11, 12)                 # the six delta = 10^-e of the acceptance (gating)
E_SUP = (3, 4, 5, 6)                           # the rest of #253's scan to 1e-3 (supplementary: the member is read before the stall has settled; gates nothing)
E_ALL = E_SUP + E_GATE
I266 = os.path.expanduser("~/Source/negative-mass-raw/I266-mac")
I255 = os.path.expanduser("~/Source/negative-mass-raw/I255-dh610")


def mut(name):
    return MUT == name


def delta_of(e):
    return 10.0 ** (-float(e))


# --------------------------------------------------------------------------------------------------------------------
# runs: p*_h (#266's pstar), the chain at delta = 10^-e, the pair p*_h -/+ 10^-e
# --------------------------------------------------------------------------------------------------------------------
def job_cmd(h, e, out, pdir):
    """#41's chain at A = p*_h - 10^-e (g.run_real_g: #40's fates at p*_h -/+ delta, then the chain), the series and the lifetimes at p*_h -/+ delta."""
    import resource
    b.set_family(FAM, B0)
    ps = b.load_pstar(pdir, FAM, B0, h)["ps"]
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    delta = delta_of(e)
    r = g.run_real_g(h, delta, ps, log=lambda m: print(m, flush=True))
    s = r["ser"]
    np.savez(os.path.join(out, "series.npz"), **s)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0 ** 2 if sys.platform == "darwin" else 1024.0)      # bytes on macOS, KB on Linux
    ok = len(s["r"]) > 0
    summ = dict(fam=FAM, B=B0, h=h, e=e, delta=delta, ps=ps, A=r["A"], status=r["status"], event=r["event"], point=r["point"], up=r["up"], Np=r["Np"], nm0=r["nm0"],
                F1=r["F1"], K=r["K"], Nc=r["Nc"], levels=len(r["levels"]), cols=int(len(s["r"])), xp0=float(s["xp"][0]) if ok else None, xp1=float(s["xp"][-1]) if ok else None,
                F_end=float(s["F"][-1]) if ok else None, M_start=float(s["M"][0]) if ok else None, M_end=float(s["M"][-1]) if ok else None,
                min_edge=float(s["edge"].min()) if ok else None, secs=round(time.time() - t0, 1), maxrss_MB=round(rss, 1))
    json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print(json.dumps(summ), flush=True)


def member_cmd(h, e, out, pdir):
    """#266's `member` at the pair p*_h -/+ 10^-e: whole runs through throat_locus.sweep, the locus, the profile slices every 0.5 in t, the member fitted on each."""
    b.member_cmd(FAM, B0, h, out, pdir, e=float(e))


# --------------------------------------------------------------------------------------------------------------------
# declared bars (the raw README, written before the runs of record)
# --------------------------------------------------------------------------------------------------------------------
CRITM = dict(
    bar_bracket=1e-14, bar_p266=1e-13, bar_ev=0.0, bar_flux=2e-2, bar_plat=1e-4, bar_flux05=1e-5, bar_hpair=5e-5, order_band=(1.5, 2.5),
    bar_r266_M=1e-6, bar_r266_mem=1e-6, bar_r253=1e-4,
    tg_band=(1.8, 2.5), tdep_band=(1.9, 2.1),
    lin_t=(10.0, 11.0), lin_e=(8, 9, 10), lin_e_ph=(11, 12), bar_lin=2e-3, bar_lin_ph=5e-2,
    traj_t=(11.0, 12.0, 13.0), bar_traj_g=1e-5, bar_traj_r=2e-5, traj_odd=1e-3,
    eps_dep=1e-3, eps_alt=1e-2, bar_fitq=1e-4, bar_hst_rel=1e-4,
    bar_pred=6e-5, bar_smooth=6e-5,
    dr_win=(6.5, 12.0), dr_odd=3e-4, dr_ext_win=(6.5, 20.5), dr_ext_t=(14.0, 20.0), bar_dr255_inf=2e-5, bar_dr255_kap=1e-2, bar_dr255_r=1e-5, bar_drext=1e-4, bar_drext_inf=2e-4,
    bar_ptd_g=1e-5, bar_ptd_r=5e-5, bar_pdr=1e-5,
)
REL_B = (0.025, 0.05, 0.1, 0.125)             # the four B of #266's relation (REL-ph there)
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


def is_ph(e, h):
    """Gating: the six delta of the issue (E_GATE) at its two h. Post hoc, gates nothing: the rest of #253's scan to 1e-3, and h = 0.02."""
    return (e not in E_GATE) or (h not in H2)


class Acc:
    def __init__(self):
        self.g, self.p, self.ng, self.np = True, True, 0, 0

    def add(self, good, e, h):
        if is_ph(e, h):
            self.p, self.np = self.p and bool(good), self.np + 1
        else:
            self.g, self.ng = self.g and bool(good), self.ng + 1


def emit(cid, text, acc, detail=""):
    if acc.ng:
        check(cid, text, acc.g, detail)
    if acc.np:
        check(cid + "-ph", "POST HOC, gates nothing (delta = 1e-3 .. 1e-6 and/or h = 0.02): " + text, acc.p, detail, ph=True)


# --------------------------------------------------------------------------------------------------------------------
# loading and the member at departure
# --------------------------------------------------------------------------------------------------------------------
def base_of(raw, h):
    return os.path.join(raw, "runs", TAG, f"h{h:g}")


def load_all(raw):
    D = {}
    for h in H3:
        base = base_of(raw, h)
        try:
            p = json.load(open(base + "/pstar/summary.json"))
            pts = [json.loads(ln) for ln in open(base + "/pstar/points.jsonl")]
        except Exception as ex:
            p, pts = None, None
            print(f"    pstar of h = {h:g}: {ex!r}")
        for e in E_ALL:
            d = dict(e=e, h=h, base=f"{base}/e{e:g}", p=p, pts=pts)
            try:
                d["i"] = json.load(open(d["base"] + "/chain/summary.json"))
                d["s"] = dict(np.load(d["base"] + "/chain/series.npz"))
                d["m"] = json.load(open(d["base"] + "/member/member.json"))
                d["loc"] = dict(np.load(d["base"] + "/member/locus.npz"))
                if p is None:
                    raise RuntimeError("no pstar")
            except Exception as ex:
                d["err"] = repr(ex)
            D[(e, h)] = d
    return D


CHAIN_RERUN = {"family_B", "mass_exponent", "short_limit", "delta_off"}
REF = (8, 0.01)                                # under a mutation that acts on the chain, this one run is rerun; every other run is read


def chain_of(d):
    """The stored series and summary, or, under a mutation that acts on the chain, a fresh run of the reference point (e = 8, h = 0.01) from the stored p*_h."""
    if mut("no_settle"):
        keep = d["s"]["xp"] <= d["s"]["xp"][-1] - 3.0
        return {k: (v[keep] if isinstance(v, np.ndarray) and v.ndim == 1 and len(v) == len(keep) else v) for k, v in d["s"].items()}, d["i"]
    if MUT not in CHAIN_RERUN or (d["e"], d["h"]) != REF:
        return d["s"], d["i"]
    b.set_family(FAM, B0)
    if mut("family_B"):
        th.CRIT["G_B"] = 0.11
    e_use = d["e"] + 1 if mut("delta_off") else d["e"]
    r = g.run_real_g(d["h"], delta_of(e_use), d["p"]["ps"], K=12 if mut("short_limit") else None)
    return r["ser"], dict(status=r["status"], event=r["event"], point=r["point"], up=r["up"], Nc=r["Nc"], F1=r["F1"], K=r["K"], h=d["h"])


def pair_arrays(loc):
    """The pair's mean on the slices both runs have (the tuned member to second order in the offset): t, tau, r_th, and odd = (r_+ - r_-)/2, the unstable mode's part."""
    sa, sb = loc["s_a"], loc["s_b"]
    common = np.intersect1d(sa, sb)
    ia, ib = np.searchsorted(sa, common), np.searchsorted(sb, common)
    t = 0.5 * (loc["t_a"][ia] + loc["t_b"][ib])
    tau = 0.5 * (loc["tau_a"][ia] + loc["tau_b"][ib])
    if mut("tau_unit"):
        tau = tau * SQ2
    if mut("rdot_nolapse"):
        tau = t.copy()
    r = 0.5 * (loc["r_a"][ia] + loc["r_b"][ib])
    odd = 0.5 * (loc["r_b"][ib] - loc["r_a"][ia])
    return t, tau, r, odd


def derive(D):
    for d in D.values():
        if "err" in d:
            continue
        try:
            d["S"], d["I"] = chain_of(d)
            d["SS"] = g.summary_g(d["S"], d["h"])
            d["mem"] = dict(m3=b.member_at(d, CRITM["eps_dep"] * (10.0 if (mut("tdep_e9") and d["e"] == 9) else 1.0), "m"), m2=b.member_at(d, CRITM["eps_alt"], "m"), a3=b.member_at(d, CRITM["eps_dep"], "a"))
        except Exception as ex:
            d["err"] = "derive: " + repr(ex)
    return D


def load266(h):
    """#266's G members (stored runs, read only) at this h: the plateau of M_AH and the member at departure (eps = 1e-3), B = 0 and the four B of its relation."""
    out = {}
    for B in (0.0,) + REL_B:
        base = f"{I266}/runs/G_B{B:g}/h{h:g}"
        e = dict(h=h, loc=dict(np.load(base + "/member/locus.npz")), m=json.load(open(base + "/member/member.json")))
        s = dict(np.load(base + "/chain/series.npz"))
        out[B] = dict(e=e, plat=g.summary_g(s, h)["plat1"], mem=b.member_at(e, CRITM["eps_dep"], "m"), mem2=b.member_at(e, CRITM["eps_alt"], "m"), i=json.load(open(base + "/chain/summary.json")))
        out[B]["speed"] = speed_at(e["loc"], out[B]["mem"]["t_dep"])
    return out


def speed_at(loc, t0):
    """dr_th/dtau of the pair's mean at the time t0 (#255's local slope, interpolated)."""
    t, tau, r, odd = pair_arrays(loc)
    rd = td.rdot_local(tau, r)
    ok = np.isfinite(rd)
    return float(np.interp(t0, t[ok], rd[ok])) if np.isfinite(t0) else float("nan")


def fit_relation(R, scale=True, key="mem"):
    """#266's relation by least squares through the four B: M / r_th - M_0 = a gamma_1 + b gamma_1^2 (scale=False: M - M_0, without the division by r_th). key: the member
    reading, 'mem' (eps = 1e-3) or 'mem2' (eps = 1e-2)."""
    c0 = R[0.0]["plat"]
    g1 = np.array([R[B_][key]["g1"] for B_ in REL_B])
    rth = np.array([R[B_][key]["r_th"] for B_ in REL_B]) if scale else np.ones(len(REL_B))
    M = np.array([R[B_]["plat"] for B_ in REL_B])
    if not (np.all(np.isfinite(g1)) and np.all(np.isfinite(rth))):
        return c0, np.array([float("nan"), float("nan")])          # a mutated reading outside the stored slices
    ab = np.linalg.lstsq(np.stack([g1, g1 ** 2], 1), M / rth - c0, rcond=None)[0]
    return c0, ab


def predict(c0, ab, g1, rth):
    return rth * (c0 + ab[0] * g1 + ab[1] * g1 ** 2)


def scan_errors(entries, c0, ab, scale=True):
    """entries: list of (M, g1, r_th) over the delta scan. Prediction errors M_pred - M (the unscaled relation is M - M_0 = a g + b g^2, so M_pred = M_0 + ...)."""
    out = []
    for M, g1, rth in entries:
        pr = rth * (c0 + ab[0] * g1 + ab[1] * g1 ** 2) if scale else c0 + ab[0] * g1 + ab[1] * g1 ** 2
        out.append(pr - M)
    return out


# --------------------------------------------------------------------------------------------------------------------
# the throat radius's speed along the stall (#255's estimator)
# --------------------------------------------------------------------------------------------------------------------
def drift_fit(loc, win=None, odd_max=None):
    t, tau, r, odd = pair_arrays(loc)
    win = ((8.0, 12.0) if mut("drift_win") else CRITM["dr_win"]) if win is None else win
    odd_max = CRITM["dr_odd"] if odd_max is None else odd_max
    f = td.fit_rdot(t, tau, r, odd, win=win, odd_max=odd_max)
    if np.isfinite(f["t0"]):
        f["tau0"] = float(tau[int(np.argmin(np.abs(t - f["t0"])))])
    else:
        f["tau0"] = float("nan")
    return f, (t, tau, r, odd)


def rdot_model(f, tau):
    return f["inf"] + f["c"] * np.exp(-f["kappa"] * (tau - f["tau0"]))


def drift_extrap(loc, win=None, ext_t=None):
    """The fit on #255's window, extrapolated beyond it: max |measured local slope - model| over t in ext_t where the pair's odd part is within the bar; and the slope there."""
    win = ((10.5, 12.0) if mut("drext_win") else CRITM["dr_win"]) if win is None else win
    ext_t = CRITM["dr_ext_t"] if ext_t is None else ext_t
    f, (t, tau, r, odd) = drift_fit(loc, win)
    rd = td.rdot_local(tau, r)
    sel = np.isfinite(rd) & (t >= ext_t[0]) & (t <= ext_t[1]) & (np.abs(odd) <= CRITM["dr_odd"])
    if sel.sum() < 5:
        return f, float("nan"), 0, None
    err = float(np.max(np.abs(rd[sel] - rdot_model(f, tau[sel]))))
    return f, err, int(sel.sum()), (t[sel], rd[sel], rdot_model(f, tau[sel]))


def parse255():
    """#255's D = 4, B = 0.1, h = 0.01 line and its D = 5, B = 0.1, finest-h rdot_inf, from throat_drift_output.txt (the file of record)."""
    import re
    txt = open(os.path.join(HERE, "throat_drift_output.txt")).read()
    m4 = re.search(r"D=4 B=0\.1 h=0\.01: p\*_h ([0-9.e+-]+); r_th - 1 at t = \(8\.0, 10\.0, 12\.0\): \[([^\]]+)\]; rdot_inf ([+-][0-9.]+) \+- ([0-9.e+-]+), kappa ([0-9.]+), n (\d+)", txt)
    m5 = re.search(r"rdot_inf \+0\.00144 \+- [0-9.e+-]+, (\+0\.02602) \+- ([0-9.e+-]+)", txt)
    if not m4 or not m5:
        raise RuntimeError("throat_drift_output.txt: the lines of #255's D = 4 and D = 5 fits were not found")
    R = [float(x) for x in m4.group(2).split(",")]
    return dict(ps=float(m4.group(1)), R=R, inf=float(m4.group(3)), se=float(m4.group(4)), kappa=float(m4.group(5)), n=int(m4.group(6)), inf5=float(m5.group(1)), se5=float(m5.group(2)))


# --------------------------------------------------------------------------------------------------------------------
# planted cases through the same code
# --------------------------------------------------------------------------------------------------------------------
LAM_T = 1.15                                   # growth rate of the unstable mode per unit t (planted; the measured pair rate is 1.18 per unit proper time)
Q0 = 1.6e-2                                   # planted odd part at t = 0 per unit offset (odd = Q0 * delta * exp(LAM_T t) r_th)


def planted_traj(t):
    """A planted member trajectory of the size of the measured one: r_th(t) with a minimum near t = 13.5 and a late speed 1.1e-3, gamma_1(t) falling slowly."""
    r = 1.0029 + 2.5e-4 * (t - 13.5) ** 2 / (1.0 + 0.02 * (t - 13.5) ** 2) + 3.0e-4 * np.clip(t - 13.5, 0, None)
    g1 = 0.0483 - 6.0e-5 * (t - 13.5)
    return r, g1


def planted_scan_entry(h, e, c0, ab, extra=0.0):
    """A planted run of the delta scan: pair loc arrays with the planted odd part, rows every 0.5 in t with the exact static member of (r_th(t), gamma_1(t)),
    and M = r_th (c0 + a gamma_1 + b gamma_1^2) at the planted departure (plus extra * (e - 8))."""
    delta = delta_of(e)
    dt = h / SQ2
    t = np.arange(2.0, 34.0, dt)
    rt, g1t = planted_traj(t)
    q = Q0 * delta * np.exp(LAM_T * t) * rt
    s = np.arange(len(t))
    loc = dict(s_a=s, s_b=s, t_a=t, t_b=t, tau_a=0.97 * t, tau_b=0.97 * t, r_a=rt - q, r_b=rt + q)
    rows = []
    r_pts = np.array(tl.STALL["r_pts"], float)
    for tk in np.arange(3.0, 30.0 + 1e-9, 0.5):
        rth_k, g_k = (float(x) for x in planted_traj(np.array(tk)))
        s_k = g_k / math.sqrt(1 + g_k * g_k)
        Mp, Mm = tl.member_M(ND, rth_k, s_k, r_pts * rth_k)
        rows.append(dict(s=int(round(tk * SQ2 / h)), t=float(tk), r_th=rth_k, r_th_a=rth_k, rr_m=(r_pts * rth_k).tolist(), Mp_m=Mp.tolist(), Mm_m=Mm.tolist(),
                         rr_a=(r_pts * rth_k).tolist(), Mp_a=Mp.tolist(), Mm_a=Mm.tolist()))
    d = dict(e=e, h=h, loc=loc, m=dict(rows=rows))
    mem = b.member_at(d, CRITM["eps_dep"], "m")
    td_true = math.log(CRITM["eps_dep"] / (Q0 * delta)) / (LAM_T * (1.05 if mut("ptd_lambda") else 1.0))
    r_true, g_true = (float(x) for x in planted_traj(np.array(td_true)))
    M = r_true * (c0 + ab[0] * g_true + ab[1] * g_true ** 2) + extra * (e - 8) + (1e-4 if mut("ppos_off") else 0.0)
    return dict(mem=mem, t_true=td_true, r_true=r_true, g_true=g_true, M=M, loc=loc)


def check_planted(ref=None):
    """PMEM, PDEP, PREL of #266 through this driver; PTD: the planted delta scan through the departure time and the member at departure; PPOS, PNEG: the planted relation
    through the prediction test passes, and the same with an extra delta dependence of 3e-4 per decade fails; PDR: a planted throat speed through the drift fit."""
    n0 = len(b.results)
    b.check_planted()
    for cid, ok, ph in b.results[n0:]:
        results.append((cid, ok, ph))
    h = 0.01
    c0, ab = 0.223086, (0.2027, 0.853)
    ents = {e: planted_scan_entry(h, e, c0, ab) for e in E_GATE}
    eg = max(abs(en["mem"]["g1"] - en["g_true"]) for en in ents.values())
    er = max(abs(en["mem"]["r_th"] - en["r_true"]) for en in ents.values())
    et = max(abs(en["mem"]["t_dep"] - en["t_true"]) for en in ents.values())
    for e, en in ents.items():
        print(f"    PTD e = {e}: t_dep {en['mem']['t_dep']:.4f} (planted {en['t_true']:.4f}), r_th {en['mem']['r_th']:.6f} ({en['r_true']:.6f}), gamma_1 {en['mem']['g1']:.6f} ({en['g_true']:.6f})")
    check("PTD", "planted delta scan (odd part Q0 delta e^{lambda t}, member trajectory of the measured size, rows every 0.5 in t): the departure time, r_th and gamma_1 at departure for six delta",
          eg <= CRITM["bar_ptd_g"] and er <= CRITM["bar_ptd_r"] and et <= 2e-5, f"max error gamma_1 {eg:.1e} (bar {CRITM['bar_ptd_g']:g}), r_th {er:.1e} (bar {CRITM['bar_ptd_r']:g}), t_dep {et:.1e} (bar 2e-5)")
    errs = scan_errors([(en["M"], en["mem"]["g1"], en["mem"]["r_th"]) for en in ents.values()], c0, ab)
    check("PPOS", "planted relation through the prediction test: a scan that follows the member passes", max(abs(x) for x in errs) <= CRITM["bar_pred"], f"max |error| {max(abs(x) for x in errs):.1e} (bar {CRITM['bar_pred']:g})")
    ents2 = {e: planted_scan_entry(h, e, c0, ab, extra=3e-4) for e in E_GATE}
    errs2 = scan_errors([(en["M"], en["mem"]["g1"], en["mem"]["r_th"]) for en in ents2.values()], c0, ab)
    check("PNEG", "planted relation plus an extra 3e-4 per decade of delta through the same test: it fires", max(abs(x) for x in errs2) > CRITM["bar_pred"] if not mut("pneg_off") else max(abs(x) for x in errs2) <= CRITM["bar_pred"],
          f"max |error| {max(abs(x) for x in errs2):.1e} (the test fires above {CRITM['bar_pred']:g})")
    t = np.arange(0.0, 30.0, 0.01 / SQ2)
    inf0, c_, kap = 0.0013, -0.0030, 0.244
    tau = 0.975 * t
    r = 1.01 + inf0 * tau - (c_ / kap) * np.exp(-kap * tau) + (c_ / kap)
    loc = dict(s_a=np.arange(len(t)), s_b=np.arange(len(t)), t_a=t, t_b=t, tau_a=tau, tau_b=tau, r_a=r - 1e-12, r_b=r + 1e-12)
    ws = ((6.5, 12.0), (6.5, 20.5))
    worst = 0.0
    for w in ws:
        f, _ = drift_fit(loc, win=w)
        worst = max(worst, abs(f["inf"] - inf0), abs(f["kappa"] - kap) * 1e-1)
        print(f"    PDR window {w}: rdot_inf {f['inf']:+.6f} (planted {inf0:+.6f}), kappa {f['kappa']:.4f} ({kap})")
    fx, err, n, _ = drift_extrap(loc)
    check("PDR", "planted throat speed rdot = rdot_inf + c e^{-kappa tau} (rdot_inf = 0.0013, kappa = 0.244, as the fit of #255 at D = 4) through the drift fit on #255's window and on the extended one, and its extrapolation",
          worst <= CRITM["bar_pdr"] and err <= 1e-5, f"worst |rdot_inf error| {worst:.1e} (bar {CRITM['bar_pdr']:g}), extrapolation error {err:.1e} (bar 1e-5)")


# --------------------------------------------------------------------------------------------------------------------
# the analysis
# --------------------------------------------------------------------------------------------------------------------
def analyze_cmd(raw):
    t0 = time.time()
    D = derive(load_all(raw))
    miss = {k: d["err"] for k, d in D.items() if "err" in d}
    print(f"== Runs read from {raw}/runs: {len(D) - len(miss)} of {len(D)} (delta, h) complete" + (f"; missing: {sorted(miss)}" if miss else ""))
    ok_keys = [k for k, d in D.items() if "err" not in d]
    xg, xo = th.CRIT["X"]["G"][ND], th.CRIT["X"]["O"][ND]
    bad = {k: (sm.CRIT[k], v) for k, v in g.DECLARED.items() if k != "X" and sm.CRIT[k] != v}
    check("CFG", "settled_mass.CRIT (F1, K, W1, rcut, tol, dxp, dm) equal #41's declared chain, X = 28 for G, threshold's G_B = 0.1 at import",
          not bad and xg == g.DECLARED["X"] and b.G_B_DEFAULT == 0.1, f"differs: {bad}; X(G) {xg}, G_B {b.G_B_DEFAULT}" if bad or xg != 28.0 or b.G_B_DEFAULT != 0.1 else "")
    acc = Acc()
    for (e, h), d in D.items():
        acc.add("err" not in d, e, h)
    emit("RUNS", "every (delta, h) has its chain and member output, and every h its p*_h", acc, f"missing {sorted(miss)}" if miss else "")

    print("\n== Planted cases")
    check_planted()

    print("\n== p*_h: bisection to hi - lo <= 1e-14 (threshold.bisect on threshold.point, X = 28), and #266's p*_h of the same member")
    acc, det = Acc(), []
    acc2 = Acc()
    for k, h in enumerate(H3):
        d0 = D.get((E_GATE[0], h))
        if d0 is None or d0["p"] is None:
            acc.add(False, E_GATE[0], h)
            acc2.add(False, E_GATE[0], h)
            continue
        p, pts = d0["p"], d0["pts"]
        lo, hi = b.bracket_after(pts, 30) if mut("bisect_early") else (p["lo"], p["hi"])
        good = (hi - lo <= CRITM["bar_bracket"] and pts[0]["fate"] == -1 and pts[1]["fate"] == 1 and p["h_arg"] == repr(float(h)) and p["stopped"] is None
                and b.bracket_after(pts, len(pts)) == (p["lo"], p["hi"]))
        acc.add(good, E_GATE[0], h)
        h_ref = H3[(k + 1) % len(H3)] if mut("pstar_h_swap") else h
        ref = json.load(open(f"{I266}/runs/{TAG}/h{h_ref:g}/pstar/summary.json"))
        dd = abs(p["ps"] - ref["ps"])
        acc2.add(dd <= CRITM["bar_p266"], E_GATE[0], h)
        det.append(f"h {h:g}: hi - lo {hi - lo:.1e}, {p['runs']} runs, stopped {p['stopped']}; p* {p['ps']:.16e}, #266's {ref['ps']:.16e}, difference {dd:.1e}")
    for ln in det:
        print("    " + ln)
    emit("BR", "bracket ends are collapse and expansion, the replayed bisection equals lo and hi, hi - lo <= 1e-14, h passed as its float repr", acc)
    emit("P266", "this machine's p*_h of G (B = 0.1) equals #266's within 1e-13 (the same bisection, the same code)", acc2, "bar 1e-13")

    print("\n== The chain at p*_h - delta (#41's, unchanged): events, limit, first law, settled value")
    acc, det = Acc(), []
    for (e, h) in ok_keys:
        d = D[(e, h)]
        i, pair_ev = d["I"], d["m"]["ev40_a"].get("collapse")
        good = list(i["event"]) == list(pair_ev) and i["point"]["fate"] == -1 and i["up"]["fate"] == 1
        acc.add(good, e, h)
        if not good:
            det.append(f"e = {e} h {h:g}: chain event {list(i['event'])}, pair's collapse-side event {pair_ev}, fates {i['point']['fate']}, {i['up']['fate']}")
    for ln in det:
        print("    " + ln)
    emit("EVT", "the chain's first trapped cell equals the pair's collapse-side one (same amplitude p*_h - delta), fate collapse at p*_h - delta and expansion at p*_h + delta", acc)
    Fmax = sm.CRIT["F1"] * 2.0 ** sm.CRIT["K"]
    acc, det = Acc(), []
    for (e, h) in ok_keys:
        s, i = D[(e, h)]["S"], D[(e, h)]["I"]
        fine = np.allclose(np.diff(s["xp"]), h, atol=1e-9)
        good = i["status"] == "limit" and 0.5 * Fmax <= s["F"][-1] < Fmax and s["edge"].min() >= i["Nc"] / 8 and fine
        acc.add(good, e, h)
        det.append(f"{i['status']} F_end {s['F'][-1]:.3g}")
    emit("RC5", f"every run ends at the declared limit (F = e^-f on the horizon, {0.5 * Fmax:g} <= F_end < {Fmax:g}), horizon at least Nc/8 cells from the box edges, no gaps", acc)
    acc, parts = Acc(), []
    for (e, h) in ok_keys:
        s = D[(e, h)]["S"]
        i1 = int(np.argmax(s["F"] >= 20.0))
        dMs = s["M"][i1] - s["M"][0]
        with np.errstate(all="ignore"):
            rel = abs(dMs - g.flux_integral_g(s, 0, i1)) / abs(dMs)
        rel = rel if np.isfinite(rel) else float("inf")
        acc.add(rel <= CRITM["bar_flux"], e, h)
        parts.append(rel)
    emit("RC4", "first law on the horizon, both fields: M(x+_b) - M(x+_a) = integral of (r^2/2) e^f r_- (phi0_+^2 - phi1_+^2), from the first trapped cell to F = 20", acc, f"bar {CRITM['bar_flux']:g}; worst {max(parts):.1e}" if parts else "")
    acc = Acc()
    for (e, h) in ok_keys:
        q = D[(e, h)]["SS"]
        acc.add(abs(q["dplat"]) <= CRITM["bar_plat"] and q["flux05"] <= CRITM["bar_flux05"], e, h)
    emit("RC8", "M_AH has settled at the limit: plateau mean over the last 1.0 differs from the previous 1.0 by <= 1e-4 and the mean horizon flux over the last 0.5 is <= 1e-5", acc)
    acc, dets = Acc(), []
    for e in E_ALL:
        if (e, H2[0]) in ok_keys and (e, H2[1]) in ok_keys:
            h_a = H3[0] if mut("h_plat") else H2[1]
            dd = abs(D[(e, H2[0])]["SS"]["plat1"] - D[(e, h_a)]["SS"]["plat1"]) if (e, h_a) in D and "err" not in D[(e, h_a)] else float("inf")
            acc.add(dd <= CRITM["bar_hpair"], e, H2[0])
            dets.append(f"e = {e}: {dd:.1e}")
        else:
            acc.add(False, e, H2[0])
    emit("RH", "the settled value at two h: |plat(0.01) - plat(0.005)| (plateau mean over the last 1.0 of x+) <= 5e-5", acc, "; ".join(dets))
    acc, det = Acc(), []
    for e in E_ALL:
        if all((e, h) in ok_keys for h in H3):
            rows = sm.band_orders({h: D[(e, h)]["S"] for h in H3}, H3, g.WINDOWS, delta_of(e))
            good = all(CRITM["order_band"][0] <= x <= CRITM["order_band"][1] for _, _, _, o in rows for x in o)
            acc.p, acc.np = acc.p and bool(good), acc.np + 1
            det.append(f"e = {e}: orders " + ", ".join(f"{o[0]:.2f}" for _, _, _, o in rows))
    if acc.np:
        check("RC6-ph", f"POST HOC, gates nothing (needs h = 0.02): M_AH(x+) converges at second order at the three h, three windows, every delta (observed order in {CRITM['order_band']})", acc.p, "; ".join(det), ph=True)

    print("\n== Against #266 (reproduction of its G, B = 0.1 at delta = 1e-8) and #253's scan (rpi5)")
    acc, det = Acc(), []
    for h in H3:
        d = D.get((8, h))
        if d is None or "err" in d:
            acc.add(False, 8, h)
            continue
        base = f"{I266}/runs/{TAG}/h{h:g}"
        s266 = dict(np.load(base + "/chain/series.npz"))
        e266 = dict(h=h, loc=dict(np.load(base + "/member/locus.npz")), m=json.load(open(base + "/member/member.json")))
        m266 = b.member_at(e266, CRITM["eps_dep"], "m")
        dM = d["SS"]["plat1"] - g.summary_g(s266, h)["plat1"]
        dg, dr = d["mem"]["m3"]["g1"] - m266["g1"], d["mem"]["m3"]["r_th"] - m266["r_th"]
        acc.add(abs(dM) <= CRITM["bar_r266_M"] and abs(dg) <= CRITM["bar_r266_mem"] and abs(dr) <= CRITM["bar_r266_mem"], 8, h)
        det.append(f"h {h:g}: plateau {d['SS']['plat1']:.6f} (#266 {d['SS']['plat1'] - dM:.6f}), gamma_1 {d['mem']['m3']['g1']:.5f} ({m266['g1']:.5f}), r_th {d['mem']['m3']['r_th']:.5f} ({m266['r_th']:.5f}); differences {dM:+.1e}, {dg:+.1e}, {dr:+.1e}")
    for ln in det:
        print("    " + ln)
    emit("R266", "delta = 1e-8 reproduces #266's G (B = 0.1): the plateau of M_AH to 1e-6 and gamma_1, r_th at departure to 1e-6", acc)
    acc, det = Acc(), []
    for h in H2:
        sc = {r["delta"]: r for r in json.load(open(f"{b.I253}/runs/scan_G_h{h:g}.json"))} if os.path.exists(f"{b.I253}/runs/scan_G_h{h:g}.json") else {}
        for e in E_ALL:
            d = D.get((e, h))
            r = sc.get(delta_of(e))
            if d is None or "err" in d or r is None:
                acc.add(False, e, h)
                continue
            dd = d["SS"]["plat1"] - r["plat1"]
            acc.add(abs(dd) <= CRITM["bar_r253"], e, h)
            det.append(f"h {h:g} e = {e}: plateau {d['SS']['plat1']:.6f}, #253's {r['plat1']:.6f}, {dd:+.1e}")
    for ln in det:
        print("    " + ln)
    emit("R253", "the plateau at each delta of the scan reproduces #253's scan (rpi5, #40's p*_h) within 1e-4 (its roundoff-realization scatter is 2.5e-5 at h = 0.01)", acc, "bar 1e-4")

    print("\n== The stall: lifetimes and the unstable mode")
    print("    lifetimes T of #40's events at p*_h -/+ delta (collapse side, expansion side) and the departure time of the unstable mode (eps = 1e-3)")
    for h in H2:
        print(f"    h = {h:g}: " + "; ".join(f"e = {e}: T {D[(e, h)]['I']['point']['T']:.3f}/{D[(e, h)]['I']['up']['T']:.3f}, t_dep {D[(e, h)]['mem']['m3']['t_dep']:.3f}" for e in E_ALL if (e, h) in ok_keys))
    acc, det = Acc(), []
    for h in H3:
        for e in E_GATE[1:]:
            if (e, h) not in ok_keys or (e - 1, h) not in ok_keys:
                acc.add(False, e, h)
                continue
            Tm = {k: D[(k, h)]["I"]["point"]["T"] for k in (e - 1, e)}
            Tp = {k: D[(k, h)]["I"]["up"]["T"] for k in (e - 1, e)}
            if mut("tg_flat"):
                Tm[e], Tp[e] = Tm[e - 1], Tp[e - 1]
            gm, gp = Tm[e] - Tm[e - 1], Tp[e] - Tp[e - 1]
            acc.add(CRITM["tg_band"][0] <= gm <= CRITM["tg_band"][1] and CRITM["tg_band"][0] <= gp <= CRITM["tg_band"][1], e, h)
            det.append(f"h {h:g} e {e - 1}->{e}: +{gm:.2f}, +{gp:.2f}")
    emit("TG", f"the stall: the lifetime T at p*_h -/+ delta gains {CRITM['tg_band'][0]} to {CRITM['tg_band'][1]} per decade of delta on both sides, from 1e-7 to 1e-12", acc, "; ".join(det[:5]))
    acc, det = Acc(), []
    for h in H3:
        for e in E_GATE[1:]:
            if (e, h) not in ok_keys or (e - 1, h) not in ok_keys:
                acc.add(False, e, h)
                continue
            gd = D[(e, h)]["mem"]["m3"]["t_dep"] - D[(e - 1, h)]["mem"]["m3"]["t_dep"]
            acc.add(CRITM["tdep_band"][0] <= gd <= CRITM["tdep_band"][1], e, h)
            det.append(f"h {h:g} e {e - 1}->{e}: +{gd:.3f}")
    emit("TDEP", f"the departure time (eps = 1e-3) gains {CRITM['tdep_band'][0]} to {CRITM['tdep_band'][1]} per decade of delta, from 1e-7 to 1e-12 (ln 10 over the mode's rate 1.15)", acc, "; ".join(det[:5]))
    acc, acc2, det = Acc(), Acc(), []
    for h in H3:
        vals = {}
        for e in CRITM["lin_e"] + CRITM["lin_e_ph"]:
            if (e, h) not in ok_keys:
                continue
            t, tau, r, odd = pair_arrays(D[(e, h)]["loc"])
            scale = 2.0 if (mut("lin_scale") and e == 10) else (1.2 if (mut("lin_scale12") and e == 12) else 1.0)
            vals[e] = [scale * float(np.interp(tt, t, np.abs(odd))) / delta_of(e) for tt in CRITM["lin_t"]]
        if len(vals) < len(CRITM["lin_e"] + CRITM["lin_e_ph"]):
            acc.add(False, CRITM["lin_e"][0], h)
            acc2.add(False, CRITM["lin_e"][0], h)
            continue
        for j, tt in enumerate(CRITM["lin_t"]):
            col = np.array([vals[e][j] for e in CRITM["lin_e"]])
            spread = float((col.max() - col.min()) / col.mean())
            acc.add(spread <= CRITM["bar_lin"], CRITM["lin_e"][0], h)
            dev = [abs(vals[e][j] / col.mean() - 1.0) for e in CRITM["lin_e_ph"]]
            acc2.add(max(dev) <= CRITM["bar_lin_ph"], CRITM["lin_e"][0], h)
            det.append(f"h {h:g} t = {tt:g}: odd/delta {col.mean():.2f}, spread over e = 8..10 {spread:.1e}; e = 11, 12 differ by " + ", ".join(f"{x:.1e}" for x in dev))
    for ln in det:
        print("    " + ln)
    emit("LIN", f"the unstable mode is linear and the same for every delta: |odd|(t)/delta at t = {CRITM['lin_t']} agrees over delta = 1e-8, 1e-9, 1e-10 to {CRITM['bar_lin']:g} (relative spread)", acc)
    emit("LIN2", f"(roundoff floor) at delta = 1e-11 and 1e-12 the pair's odd part differs from the linear value by at most {CRITM['bar_lin_ph']:g} (the roundoff-realization of the evolution, #253)", acc2)

    print("\n== The member along the stall does not depend on delta before departure (the pair's mean profile, the same t)")
    acc, det = Acc(), []
    for h in H3:
        for tt in CRITM["traj_t"]:
            vals = []
            for e in E_GATE:
                if (e, h) not in ok_keys:
                    continue
                d = D[(e, h)]
                t, tau, r, odd = pair_arrays(d["loc"])
                if abs(np.interp(tt, t, odd)) / np.interp(tt, t, r) > CRITM["traj_odd"] * 0.3:
                    continue
                shift = 0.5 if (mut("traj_shift") and e == 9) else 0.0
                rows = d["m"]["rows"]
                k = int(np.argmin([abs(x["t"] - (tt + shift)) for x in rows]))
                vals.append((e, rows[k]["r_th"], rows[k]["g1_m"]))
            if len(vals) >= 2:
                sr = max(v[1] for v in vals) - min(v[1] for v in vals)
                sg = max(v[2] for v in vals) - min(v[2] for v in vals)
                acc.add(sr <= CRITM["bar_traj_r"] and sg <= CRITM["bar_traj_g"], E_GATE[0], h)
                det.append(f"h {h:g} t = {tt:g}: {len(vals)} delta, spread r_th {sr:.1e}, gamma_1 {sg:.1e}")
            else:
                acc.add(False, E_GATE[0], h)
    for ln in det:
        print("    " + ln)
    emit("TRAJ", f"the member's r_th(t) and gamma_1(t) at t = {CRITM['traj_t']} are the same for every delta whose unstable mode is below 3e-4 there: spread <= {CRITM['bar_traj_r']:g} (r_th), {CRITM['bar_traj_g']:g} (gamma_1)", acc)

    print("\n== The member at departure (#225's slice_masses and fit_asym on the pair's mean profile; t_dep: the unstable mode's part reaches eps = 1e-3 of r_th)")
    print("    h      e    t_dep   r_th      gamma_1   rms(fit)  | eps = 1e-2: t_dep  r_th      gamma_1  | collapse side alone: r_th  gamma_1")
    for h in H3:
        for e in E_ALL:
            if (e, h) not in ok_keys:
                continue
            m3, m2, a3 = (D[(e, h)]["mem"][k] for k in ("m3", "m2", "a3"))
            rms = [x["rms_m"] for x in D[(e, h)]["m"]["rows"] if abs(x["t"] - m3["t_dep"]) <= 0.5] if np.isfinite(m3["t_dep"]) else []
            print(f"    {h:<6g} {e:<4d} {m3['t_dep']:7.3f} {m3['r_th']:.6f}  {m3['g1']:.6f}  {max(rms) if rms else float('nan'):.1e}   | {m2['t_dep']:7.3f}  {m2['r_th']:.6f}  {m2['g1']:.6f} | {a3['r_th']:.6f}  {a3['g1']:.6f}")
    acc, det = Acc(), []
    for (e, h) in ok_keys:
        m3 = D[(e, h)]["mem"]["m3"]
        rms = [x["rms_m"] for x in D[(e, h)]["m"]["rows"] if abs(x["t"] - m3["t_dep"]) <= 0.5] if np.isfinite(m3["t_dep"]) else []
        acc.add(bool(rms) and max(rms) * (10.0 if mut("fitq_scale") else 1.0) <= CRITM["bar_fitq"], e, h)
    emit("FITQ", "the member's fit of M(r) on the two sides has an rms residual of at most 1e-4 at the slices around departure (a static member describes the state there)", acc)
    acc, det = Acc(), []
    for e in E_ALL:
        if all((e, h) in ok_keys for h in H2):
            a_, b_ = D[(e, H2[0])]["mem"]["m3"], D[(e, H2[1])]["mem"]["m3"]
            dg, dr = abs(a_["g1"] - b_["g1"]) / abs(b_["g1"]), abs(a_["r_th"] - b_["r_th"]) / abs(b_["r_th"])
            acc.add(np.isfinite(dg) and np.isfinite(dr) and dg <= CRITM["bar_hst_rel"] and dr <= CRITM["bar_hst_rel"], e, H2[0])
            det.append(f"e = {e}: gamma_1 {dg:.1e}, r_th {dr:.1e}")
        else:
            acc.add(False, e, H2[0])
    for ln in det:
        print("    h = 0.01 against 0.005, relative: " + ln)
    emit("HST", "the member at departure at the two h (0.01, 0.005) agrees to 1e-4 relative in gamma_1 and in r_th, for each delta (the issue's acceptance)", acc, "bar 1e-4 relative")

    print("\n== The settled mass against the member at departure; the relation of #266 as the prediction")
    pred_ok = {}
    for h in H2:
        R = load266(h)
        c0, ab = fit_relation(R)
        c0n, abn = fit_relation(R, scale=False)
        print(f"    h = {h:g}: #266's relation from B = {REL_B}: M_0 = {c0:.6f}, (a, b) = ({ab[0]:.4f}, {ab[1]:.4f}); without the r_th scaling: ({abn[0]:.4f}, {abn[1]:.4f})")
        rows = []
        for e in E_ALL:
            if (e, h) not in ok_keys:
                continue
            q = D[(e, h)]["mem"]["m3"]
            M = D[(e, h)]["SS"]["plat1"]
            Mp = predict(c0, ab, q["g1"], q["r_th"])
            Mn = c0n + abn[0] * q["g1"] + abn[1] * q["g1"] ** 2
            rows.append((e, M, q["g1"], q["r_th"], Mp, Mn))
        pred_ok[h] = (c0, ab, rows)
        print("      e   M(plateau)  gamma_1    r_th       M predicted  M - pred   | without r_th scaling: M pred   M - pred")
        for e, M, g1, rth, Mp, Mn in rows:
            print(f"      {e:<3d} {M:.6f}   {g1:.6f}  {rth:.6f}  {Mp:.6f}    {M - Mp:+.2e}  |                       {Mn:.6f}  {M - Mn:+.2e}")
    acc = Acc()
    for h in H2:
        c0, ab, rows = pred_ok[h]
        for e, M, g1, rth, Mp, Mn in rows:
            acc.add(abs(M - Mp) <= CRITM["bar_pred"], e, h)
    emit("PRED", "the declared prediction: M_AH settled at each delta equals r_th (M_0 + a gamma_1 + b gamma_1^2) at the member at departure, with (M_0, a, b) from #266's four B at the same h, within 6e-5", acc, "bar 6e-5")
    acc, det = Acc(), []
    for h in H2:
        use = [r for r in pred_ok[h][2] if r[0] in E_GATE]
        if len(use) != len(E_GATE):
            acc.add(False, E_GATE[0], h)
            continue
        Ms = np.array([r[1] for r in use])
        if mut("smooth_swap"):
            Ms[[3, 4]] = Ms[[4, 3]]
        A = np.stack([np.ones(len(use)), [r[2] for r in use], [r[3] for r in use]], 1)
        if not np.all(np.isfinite(A)):
            acc.add(False, E_GATE[0], h)
            continue
        co = np.linalg.lstsq(A, Ms, rcond=None)[0]
        res = float(np.max(np.abs(A @ co - Ms)))
        sv = np.linalg.svd(A - A.mean(0) * np.array([0, 1, 1]), compute_uv=False)
        acc.add(res <= CRITM["bar_smooth"], E_GATE[0], h)
        det.append(f"h {h:g}: M = {co[0]:.5f} + {co[1]:+.3f} gamma_1 {co[2]:+.3f} r_th, residual max {res:.1e}")
    for ln in det:
        print("    " + ln)
    emit("SMOOTH", "(a weak test: the six (gamma_1, r_th, M) of the scan lie on one smooth function; they lie on a curve, so any smooth M(delta) passes) M = c_0 + p gamma_1 + q r_th through the six delta, residual within 6e-5", acc, "bar 6e-5")
    c0, ab, rows = pred_ok[H2[1]]
    ref = [r for r in rows if r[0] == 8]
    print("\n    Decomposition at h = 0.005 against delta = 1e-8 (the change of M between the scan's points, and the part the relation assigns to r_th and to gamma_1):")
    print("      e   dM measured   dM predicted   from r_th     from gamma_1   unexplained")
    if ref:
        _, M8, g8, r8, Mp8, _ = ref[0]
        for e, M, g1, rth, Mp, Mn in rows:
            d_r = (M8 / r8) * (rth - r8)
            d_g = r8 * (ab[0] + 2 * ab[1] * g8) * (g1 - g8)
            print(f"      {e:<3d} {M - M8:+.3e}   {Mp - Mp8:+.3e}    {d_r:+.3e}   {d_g:+.3e}   {M - M8 - (Mp - Mp8):+.3e}")
    use = [r for r in rows if r[0] in E_GATE]
    if len(use) >= 4:
        dg_ = np.array([r[2] - g8 for r in use])
        dr_ = np.array([r[3] - r8 for r in use])
        dM_ = np.array([r[1] - M8 for r in use])
        A = np.stack([dg_, dr_], 1)
        sv = np.linalg.svd(A, compute_uv=False)
        pq = np.linalg.lstsq(A, dM_, rcond=None)[0]
        print(f"\n    The scan alone (six delta, h = 0.005): dM = p dgamma_1 + q dr_th by least squares: p = {pq[0]:+.3f}, q = {pq[1]:+.3f}, residual max {np.max(np.abs(A @ pq - dM_)):.1e}, singular values {sv[0]:.2e}, {sv[1]:.2e} (condition {sv[0] / sv[1]:.0f});"
              f" #266's relation gives p = {r8 * (ab[0] + 2 * ab[1] * g8):+.3f}, q = {M8 / r8:+.3f}")
        print("    *the delta scan is a one-parameter curve in (gamma_1, r_th): a smooth function of the member along it is not a test; the cross-family prediction is*")
        un = np.array([r[1] - r[4] for r in use])
        print(f"    PRED errors M - predicted at the six delta: {', '.join(f'{x:+.1e}' for x in un)}")


    print("\n== Post hoc (found in the development runs; declared before the runs of record, gates nothing): the same prediction with the member read at eps = 1e-2, the part of dM not carried by r_th, and a third parameter")
    h = H2[1]
    R = load266(h)
    c0, ab = fit_relation(R)
    c2, ab2 = fit_relation(R, key="mem2")
    print(f"    h = {h:g}: #266's relation with the member read at eps = 1e-2: M_0 = {c2:.6f}, (a, b) = ({ab2[0]:.4f}, {ab2[1]:.4f}); errors M - predicted at e = 7..12 (member of the scan read at eps = 1e-2): "
          + ", ".join(f"{D[(e, h)]['SS']['plat1'] - predict(c2, ab2, D[(e, h)]['mem']['m2']['g1'], D[(e, h)]['mem']['m2']['r_th']):+.1e}" for e in E_GATE if (e, h) in ok_keys))
    print("    the part of dM between neighbouring delta not carried by r_th: excess = dM - (M / r_th) dr_th, and its ratio to the previous decade's")
    kap = None
    if (12, h) in ok_keys:
        fk, _ = drift_fit(D[(12, h)]["loc"])
        kap = fk["kappa"]
    prev = None
    for e in E_GATE[1:]:
        if (e, h) not in ok_keys or (e - 1, h) not in ok_keys:
            continue
        a_, b_ = D[(e - 1, h)], D[(e, h)]
        dM = b_["SS"]["plat1"] - a_["SS"]["plat1"]
        dr = b_["mem"]["m3"]["r_th"] - a_["mem"]["m3"]["r_th"]
        mr = 0.5 * (a_["SS"]["plat1"] / a_["mem"]["m3"]["r_th"] + b_["SS"]["plat1"] / b_["mem"]["m3"]["r_th"])
        ex = dM - mr * dr
        dt_ = b_["mem"]["m3"]["t_dep"] - a_["mem"]["m3"]["t_dep"]
        print(f"      e {e - 1}->{e}: dM {dM:+.3e}, dr_th {dr:+.3e}, M/r_th {mr:.5f}, excess {ex:+.3e}" + (f", ratio {ex / prev:.3f}" if prev else "") + (f"  (exp(-kappa dtau) = {math.exp(-kap * 0.975 * dt_):.3f}, kappa {kap:.3f} of the fit on [6.5, 12], dtau/dt 0.975)" if kap else ""))
        prev = ex
    pts = []
    for B_ in REL_B:
        r_ = R[B_]
        pts.append((f"B {B_:g}", r_["mem"]["g1"], r_["mem"]["r_th"], r_["speed"], r_["plat"]))
    for e in E_GATE:
        if e != 8 and (e, h) in ok_keys:
            d = D[(e, h)]
            pts.append((f"e {e}", d["mem"]["m3"]["g1"], d["mem"]["m3"]["r_th"], speed_at(d["loc"], d["mem"]["m3"]["t_dep"]), d["SS"]["plat1"]))
    print("    the state at departure and y = M / r_th - M_0 (M_0 = B = 0's plateau): point, gamma_1, r_th, dr_th/dtau, y")
    for nm, g1, rth, sp, M in pts:
        print(f"      {nm:7s} {g1:.5f} {rth:.5f} {sp:+.5f} {M / rth - c0:+.6f}")
    y = np.array([p_[4] / p_[2] - c0 for p_ in pts])
    g1v = np.array([p_[1] for p_ in pts])
    spv = np.array([p_[3] for p_ in pts])
    for nm, A in (("y = a g + b g^2", np.stack([g1v, g1v ** 2], 1)), ("y = a g + b g^2 + c rdot", np.stack([g1v, g1v ** 2, spv], 1))):
        if not (np.all(np.isfinite(A)) and np.all(np.isfinite(y))):
            print(f"    fit {nm}: not computable (a mutated run)")
            continue
        loo = []
        for k in range(len(y)):
            idx = [j for j in range(len(y)) if j != k]
            co = np.linalg.lstsq(A[idx], y[idx], rcond=None)[0]
            loo.append(A[k] @ co - y[k])
        co = np.linalg.lstsq(A, y, rcond=None)[0]
        print(f"    fit over these {len(y)} points, {nm}: coefficients {', '.join(f'{x:+.4f}' for x in co)}; residual max {np.max(np.abs(A @ co - y)):.1e}; leave-one-out errors " + ", ".join(f"{x:+.1e}" for x in loo))
    print("    *post hoc: the third parameter was suggested by the scan itself, so these are descriptions, not tests*")

    print("\n== The throat radius along the stall at D = 4: #255's estimator (dr_th/dtau = rdot_inf + c exp(-kappa (tau - tau_0)) over t in [6.5, 12] where |odd| <= 3e-4), on the pair's mean")
    r255 = parse255()
    print(f"    #255 (throat_drift_output.txt), D = 4, B = 0.1, h = 0.01, pair 1e-10: r_th - 1 at t = 8, 10, 12 = {r255['R']}, rdot_inf {r255['inf']:+.5f} +- {r255['se']:.1e}, kappa {r255['kappa']}, n {r255['n']}; D = 5, B = 0.1 (finest h): rdot_inf {r255['inf5']:+.5f} +- {r255['se5']:.1e}")
    acc, det = Acc(), []
    for h in H2:
        d = D.get((10, h))
        if d is None or "err" in d:
            acc.add(False, 10, h)
            continue
        f, (t, tau, r, odd) = drift_fit(d["loc"])
        R = [float(np.interp(tt, t, r)) - 1.0 for tt in (8.0, 10.0, 12.0)]
        good = (abs(f["inf"] - r255["inf"]) <= CRITM["bar_dr255_inf"] and abs(f["kappa"] - r255["kappa"]) <= CRITM["bar_dr255_kap"]
                and max(abs(x - y) for x, y in zip(R, r255["R"])) <= CRITM["bar_dr255_r"])
        acc.add(good, 10, h)
        det.append(f"h {h:g} (pair 1e-10): rdot_inf {f['inf']:+.5f} +- {f['se_inf']:.1e}, kappa {f['kappa']:.3f}, n {f['n']}, r_th - 1 = " + ", ".join(f"{x:+.5f}" for x in R))
    for ln in det:
        print("    " + ln)
    emit("DR255", "#255's D = 4 fit (B = 0.1, h = 0.01, pair 1e-10) is reproduced on this machine: rdot_inf within 2e-5, kappa within 0.01, r_th - 1 at t = 8, 10, 12 within 1e-5 (at h = 0.01 and 0.005)", acc)
    acc, det = Acc(), []
    for h in H2:
        d = D.get((12, h))
        if d is None or "err" in d:
            acc.add(False, 12, h)
            continue
        f, err, n, tr = drift_extrap(d["loc"])
        fe, _ = drift_fit(d["loc"], win=CRITM["dr_ext_win"])
        good = np.isfinite(err) and err <= CRITM["bar_drext"] and abs(fe["inf"] - f["inf"]) <= CRITM["bar_drext_inf"]
        acc.add(good, 12, h)
        det.append(f"h {h:g} (pair 1e-12): fit on [6.5, 12]: rdot_inf {f['inf']:+.5f} +- {f['se_inf']:.1e}, kappa {f['kappa']:.3f}; extended to t = 20.5: rdot_inf {fe['inf']:+.5f} +- {fe['se_inf']:.1e}, kappa {fe['kappa']:.3f}, n {fe['n']}; "
                   f"extrapolated model against the measured slope over t = {CRITM['dr_ext_t']}: max |difference| {err:.1e} ({n} slices)")
    for ln in det:
        print("    " + ln)
    emit("DREXT", f"the drift continues as #255's fit says: the model fitted on [6.5, 12] predicts the measured dr_th/dtau over t = {CRITM['dr_ext_t']} (pair 1e-12) within {CRITM['bar_drext']:g}, and the fit on [6.5, 20.5] gives rdot_inf within {CRITM['bar_drext_inf']:g} of it", acc)
    d = D.get((12, H2[1]))
    if d is not None and "err" not in d:
        f, (t, tau, r, odd) = drift_fit(d["loc"])
        fe, _ = drift_fit(d["loc"], win=CRITM["dr_ext_win"])
        rd = td.rdot_local(tau, r)
        print("    dr_th/dtau along the stall (h = 0.005, pair 1e-12; the pair's odd part below 3e-4 for t <= 20.5):")
        for tt in (8, 10, 12, 14, 16, 18, 20, 22, 24):
            k = int(np.argmin(np.abs(t - tt)))
            near = rd[max(0, k - 20):k + 21]
            v = near[np.isfinite(near)]
            print(f"      t = {tt:2d}: tau {tau[k]:7.3f}, r_th {r[k]:.6f}, dr_th/dtau {v[len(v) // 2] if len(v) else float('nan'):+.5f}, |odd|/r_th {abs(odd[k]) / r[k]:.1e}")
        print(f"    D = 4 late speed rdot_inf = {fe['inf']:+.5f} +- {fe['se_inf']:.1e} (window to t = 20.5), {f['inf']:+.5f} (window to 12); #255's D = 5: {r255['inf5']:+.5f}; ratio D = 5 / D = 4: {r255['inf5'] / fe['inf']:.1f} (window to 20.5), {r255['inf5'] / f['inf']:.1f} (to 12)")
        g1s = [(row["t"], row["g1_m"]) for row in d["m"]["rows"] if 11.0 <= row["t"] <= 22.0]
        pf = np.polyfit([x[0] for x in g1s], [x[1] for x in g1s], 1)
        print(f"    gamma_1(t) along the stall: slope {pf[0]:+.2e} per unit t over t = 11 to 22 (gamma_1 {g1s[0][1]:.5f} to {g1s[-1][1]:.5f}); the speed of r_th/r_th: {fe['inf'] / np.interp(18.0, t, r):+.2e} per unit tau")
    return finish(t0)


# --------------------------------------------------------------------------------------------------------------------
# mutations: name -> the checks that must fire (computed in the development runs, declared in the raw README before the runs of record)
# --------------------------------------------------------------------------------------------------------------------
MUTATIONS = {      # name: the checks that fire (computed on the development set, raw I267-mac dev/mut, before being declared in its README; HST, PRED and the post hoc FITQ-ph, HST-ph, PRED-ph, RC6-ph fail without a mutation, so they are in every set)
    "asym_sign": ['FITQ-ph', 'HST', 'HST-ph', 'PMEM', 'PPOS', 'PRED', 'PRED-ph', 'PTD', 'RC6-ph'],
    "bisect_early": ['BR', 'BR-ph', 'FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'RC6-ph'],
    "crit_dm": ['CFG', 'FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'RC6-ph'],
    "delta_off": ['EVT', 'FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'R253', 'R266', 'RC6-ph', 'RH', 'SMOOTH', 'TG'],
    "dep_eps": ['FITQ-ph', 'HST', 'HST-ph', 'PDEP', 'PPOS', 'PRED', 'PRED-ph', 'PTD', 'RC6-ph'],
    "drext_win": ['DREXT', 'FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'RC6-ph'],
    "drift_win": ['DR255', 'FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'RC6-ph'],
    "eps_h": ['FITQ', 'FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'RC6-ph', 'SMOOTH'],
    "family_B": ['FITQ-ph', 'HST', 'HST-ph', 'LIN', 'LIN2', 'PRED', 'PRED-ph', 'R253', 'R266', 'RC6-ph', 'RH', 'RUNS', 'SMOOTH', 'TDEP', 'TG'],
    "fitq_scale": ['FITQ', 'FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'RC6-ph'],
    "flux_sign": ['FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'RC4', 'RC4-ph', 'RC6-ph'],
    "h_plat": ['FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'RC6-ph', 'RH', 'RH-ph'],
    "lin_scale": ['FITQ-ph', 'HST', 'HST-ph', 'LIN', 'LIN-ph', 'LIN2', 'LIN2-ph', 'PRED', 'PRED-ph', 'RC6-ph'],
    "lin_scale12": ['FITQ-ph', 'HST', 'HST-ph', 'LIN2', 'LIN2-ph', 'PRED', 'PRED-ph', 'RC6-ph'],
    "mass_exponent": ['FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'R253', 'R266', 'RC4', 'RC6-ph', 'RH', 'SMOOTH'],
    "no_settle": ['FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'R253', 'R253-ph', 'R266', 'R266-ph', 'RC5', 'RC5-ph', 'RC6-ph', 'RC8', 'RC8-ph'],
    "pneg_off": ['FITQ-ph', 'HST', 'HST-ph', 'PNEG', 'PRED', 'PRED-ph', 'RC6-ph'],
    "ppos_off": ['FITQ-ph', 'HST', 'HST-ph', 'PPOS', 'PRED', 'PRED-ph', 'RC6-ph'],
    "pstar_h_swap": ['FITQ-ph', 'HST', 'HST-ph', 'P266', 'P266-ph', 'PRED', 'PRED-ph', 'RC6-ph'],
    "ptd_lambda": ['FITQ-ph', 'HST', 'HST-ph', 'PPOS', 'PRED', 'PRED-ph', 'PTD', 'RC6-ph'],
    "rdot_nolapse": ['DR255', 'FITQ-ph', 'HST', 'HST-ph', 'PDR', 'PRED', 'PRED-ph', 'RC6-ph'],
    "rel_perm": ['FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'PREL', 'RC6-ph'],
    "short_limit": ['FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'R253', 'R266', 'RC5', 'RC6-ph', 'RC8', 'RH', 'SMOOTH'],
    "smooth_swap": ['FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'RC6-ph', 'SMOOTH'],
    "tau_unit": ['DR255', 'FITQ-ph', 'HST', 'HST-ph', 'PDR', 'PRED', 'PRED-ph', 'RC6-ph'],
    "tdep_e9": ['FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'RC6-ph', 'SMOOTH', 'TDEP', 'TDEP-ph'],
    "tg_flat": ['FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'RC6-ph', 'TG', 'TG-ph'],
    "traj_shift": ['FITQ-ph', 'HST', 'HST-ph', 'PRED', 'PRED-ph', 'RC6-ph', 'TRAJ', 'TRAJ-ph'],
}


def main():
    global MUT
    argv = sys.argv[1:]
    cmd = argv[0] if argv else "help"
    opts = dict(a.lstrip("-").split("=", 1) if "=" in a else (a.lstrip("-"), True) for a in argv[1:])
    if "mutate" in opts:
        MUT = opts["mutate"]
        if MUT not in MUTATIONS:
            raise SystemExit(f"unknown mutation {MUT}; known: {', '.join(sorted(MUTATIONS))}")
        b.MUT = g.MUT = sm.MUT = MUT
        dn.MUT = MUT if MUT in ("theta_coef", "ghost_sign", "scalar_coef") else None
        tl.LMUT = MUT if MUT in ("asym_sign",) else None
        if MUT == "crit_dm":
            sm.CRIT["dm"] = 2.0
    import platform
    import numba as nb
    print(f"settled_mass_member.py {cmd} {opts}: numba {nb.__version__}, numpy {np.__version__}, python {platform.python_version()}, {platform.node()}, "
          f"mutation {MUT}, {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    if cmd == "pstar":
        b.pstar_cmd(FAM, B0, float(opts["h"]), opts["out"])
    elif cmd == "job":
        job_cmd(float(opts["h"]), float(opts["e"]), opts["out"], opts["pstar"])
    elif cmd == "member":
        member_cmd(float(opts["h"]), float(opts["e"]), opts["out"], opts["pstar"])
    elif cmd == "planted":
        t0 = time.time()
        check_planted()
        sys.exit(finish(t0))
    elif cmd == "analyze":
        sys.exit(analyze_cmd(opts["raw"]))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()

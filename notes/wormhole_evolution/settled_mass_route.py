#!/usr/bin/env python3
"""Does the settled mass of family G follow its member at departure along a second route (#283)? D = 4, delta = 1e-8.

Family G keeps its ghost pulse (field 0, centre 3, width 0.5, amplitude A tuned) and changes only the normal pulse (field 1, same ray): its amplitude B, width w
(route W) or centre xc (route C). The chain (#41, through #253's `run_real_g`), the pair p*_h -/+ 1e-8 and the member at departure (#225's method through #266's
`member_at`) and the throat's speed at departure (#267's `speed_at`) are read exactly as in #266 and #267. The two predictions of #267's post hoc work are fixed in
`CRITR` and tested on the new members.

A thin driver. It imports settled_mass_b.py (#266), settled_mass_member.py (#267), which import settled_mass.py (#41), settled_mass_g.py (#253), threshold.py (#40),
throat_locus.py (#225), dn.py (#39) and throat_drift.py (#255), all unchanged. `threshold.family` is replaced in this process by `route_family` below (the normal pulse of
family G with its own centre and width); with the default pulse (B, 3.0, 0.5) it returns what threshold.family returns (check FAM0).

    python3 settled_mass_route.py bracket --route=B,xc,w --h=0.02                      the first sign change of the fate on a geometric grid of A (development)
    python3 settled_mass_route.py pstar   --route=B,xc,w --h=0.01 --out=DIR --bracket=lo,hi    bisect p*_h to hi - lo <= 1e-14 (#40's bisection)
    python3 settled_mass_route.py job     --route=B,xc,w --h=0.01 --out=DIR --pstar=DIR        #41's chain at p*_h - 1e-8
    python3 settled_mass_route.py member  --route=B,xc,w --h=0.01 --out=DIR --pstar=DIR        the pair p*_h -/+ 1e-8 kept whole: locus, slices, member
    python3 settled_mass_route.py planted [--mutate=NAME]                              planted cases
    python3 settled_mass_route.py analyze --raw=DIR [--mutate=NAME]                    every check and table from the stored runs

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
import settled_mass_member as m    # noqa: E402  (#267's driver; imports throat_drift unchanged)
dn, sm, g, th, tl, td = b.dn, b.sm, b.g, b.th, b.tl, m.td

SQ2 = math.sqrt(2.0)
MUT = None
ND = 2
FAM = "G"
H2 = (0.01, 0.005)                             # the issue's two h
DELTA = 1e-8
E_PAIR = 8.0
I266 = m.I266
GHOST = dict(xc=3.0, w=0.5)                    # the ghost pulse, as threshold.CRIT["pulse"]


def mut(name):
    return MUT == name


# --------------------------------------------------------------------------------------------------------------------
# the routes: family G with the normal pulse's own amplitude, centre and width
# --------------------------------------------------------------------------------------------------------------------
_family_orig = th.family


def route_family(fam, A):
    """(pulses, symmetric, eps) of family G with the normal pulse th.CRIT['G_normal'] = dict(B, xc, w); the ghost pulse (field 0, amplitude A) as threshold.family."""
    if fam != "G":
        return _family_orig(fam, A)
    nrm = th.CRIT["G_normal"]
    pu = th.CRIT["pulse"]
    ghost = dict(field=0, ray="+", A=A, xc=pu["xc"], w=pu["w"])
    return [ghost, dict(field=1, ray="+", A=nrm["B"], xc=nrm["xc"], w=nrm["w"])], False, (-1.0, 1.0)


def set_route(fam, spec):
    """spec = (B, xc, w). Selects family G with that normal pulse in every imported module (they all call threshold.family)."""
    assert fam == "G", fam
    B, xc, w = spec
    th.CRIT["G_normal"] = dict(B=B, xc=xc + (0.1 if mut("xc_offset") else 0.0), w=w * (1.1 if mut("w_scale") else 1.0))
    g.FAM = fam


def tag_of(fam, spec):
    B, xc, w = spec
    return f"R_B{B:g}_xc{xc:g}_w{w:g}"


th.family = route_family
b.set_family = set_route
b.tag_of = tag_of


def parse_route(s):
    B, xc, w = (float(x) for x in s.split(","))
    return (B, xc, w)


def bracket_scan(spec, h, grid=None):
    """First sign change of the fate (collapse -1 below, expansion +1 above) over a geometric grid of A; returns (lo, hi, rows). Development."""
    set_route(FAM, spec)
    X = th.CRIT["X"][FAM][ND]
    grid = grid or [1e-4 * 2.0 ** k for k in range(0, 11)]
    rows = []
    for p in grid:
        r = th.point(ND, FAM, h, X, p)
        rows.append(dict(A=p, fate=r["fate"], T=r["T"]))
        if len(rows) >= 2 and rows[-2]["fate"] == -1 and rows[-1]["fate"] == 1:
            return rows[-2]["A"], rows[-1]["A"], rows
    return None, None, rows


def pstar_cmd(spec, h, out, bracket):
    b.BRACKET[(FAM, spec)] = bracket
    b.pstar_cmd(FAM, spec, h, out)


def job_cmd(spec, h, out, pdir):
    b.job_cmd(FAM, spec, h, out, pdir)


def member_cmd(spec, h, out, pdir):
    b.member_cmd(FAM, spec, h, out, pdir)


def load_point(root, spec, h):
    """The stored run of one route point at h: p*_h, the chain's series and summary, the pair's locus and member rows. 'err' is set when a file is missing."""
    base = os.path.join(root, tag_of(FAM, spec), f"h{h:g}")
    e = dict(fam=FAM, B=spec, tag=tag_of(FAM, spec), h=h, base=base)
    try:
        e["p"] = json.load(open(base + "/pstar/summary.json"))
        e["pts"] = [json.loads(l) for l in open(base + "/pstar/points.jsonl")]
        e["i"] = json.load(open(base + "/chain/summary.json"))
        e["s"] = dict(np.load(base + "/chain/series.npz"))
        e["m"] = json.load(open(base + "/member/member.json"))
        e["loc"] = dict(np.load(base + "/member/locus.npz"))
    except Exception as ex:
        e["err"] = repr(ex)
    return e


def read_point(e, eps=1e-3):
    """Plateau of M_AH, the member at departure (#266's member_at on the pair's mean profile) and the throat's speed there (#267's speed_at)."""
    mem = b.member_at(e, eps, "m")
    return dict(M=g.summary_g(e["s"], e["h"])["plat1"], mem=mem, sp=m.speed_at(e["loc"], mem["t_dep"]))


# --------------------------------------------------------------------------------------------------------------------
# the declared set (the raw README, written before the runs of record)
# --------------------------------------------------------------------------------------------------------------------
POINTS = (          # name, (B, xc, w), target gamma_1; W: width 1.0 at centre 3 (route W), C: centre 4.0 at width 0.5 (route C); B tuned at h = 0.02 to the target
    ("W1", (0.07905, 3.0, 1.0), 0.012), ("W2", (0.1543, 3.0, 1.0), 0.048), ("W3", (0.1909, 3.0, 1.0), 0.075),
    ("C1", (0.0393, 4.0, 0.5), 0.012), ("C2", (0.07726, 4.0, 0.5), 0.048), ("C3", (0.09576, 4.0, 0.5), 0.075),
)
REP = ("REP", (0.1, 3.0, 0.5), None)          # #266's B = 0.1 through this driver
H_FIN = H2[1]
CRITR = dict(
    bar_bracket=1e-14, bar_bracket_X=2e-14, bar_flux=2e-2, bar_plat=1e-4, bar_flux05=1e-5, bar_hpair=5e-5, stall_dT=2.0,
    tgt_rel=0.03, bar_hst_rel=1e-4, bar_hst_rdot=2e-5, bar_fitq=1e-4, eps_dep=1e-3,
    sep_k=5.0, bar_pred=6e-5,
    M0={0.005: 0.223086, 0.01: 0.223095},                    # B = 0's plateau (#266, #267's output)
    P2={0.005: (0.2027, 0.8534), 0.01: (0.2029, 0.8502)},    # #266's relation (a, b) at the same h (settled_mass_member_output.txt, "#266's relation")
    P3=(0.2300, 0.2724, 0.7835),                             # #267's post hoc fit over nine points (a, b, c), h = 0.005 (settled_mass_member_output.txt, "fit over these 9 points")
    bar_rep_M=1e-6, bar_rep_mem=1e-6, bar_rep_p=1e-13,
    bar_pm_g=1e-5, bar_pm_r=1e-6, bar_pm_m=1e-4, bar_dep=2e-5, bar_ptd_g=1e-5, bar_ptd_r=5e-5, bar_ptd_t=2e-5, bar_pspd=5e-6,
)
results = []
BASE = set()                                  # checks that fail without a mutation (from the unmutated analysis, --baseline=FILE): a mutation 'fires' a check that it newly fails
RERUN = {"w_scale", "xc_offset", "short_limit", "no_zoom", "mass_exponent", "finder_sign", "theta_coef"}      # mutations that act on the chain: rerun from the stored p*_h


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
        want = sorted(MUTATIONS.get(MUT, []))
        print(f"mutation {MUT}: fired {fired}; failing without a mutation (--baseline) {sorted(BASE)}; newly fired {new}, declared {want}: {'MATCH' if new == want else 'DIFFERENT'}")
        return 0 if new == want else 1
    return 1 if nfail else 0


def is_ph(name, h):
    """Gating: the six points at both h and the reproduction point at h = 0.005. Post hoc: the reproduction point at h = 0.01."""
    return name == "REP" and h != H_FIN


class Acc:
    def __init__(self):
        self.g, self.p, self.ng, self.np = True, True, 0, 0

    def add(self, good, name, h):
        if is_ph(name, h):
            self.p, self.np = self.p and bool(good), self.np + 1
        else:
            self.g, self.ng = self.g and bool(good), self.ng + 1

    def add_ph(self, good):
        self.p, self.np = self.p and bool(good), self.np + 1


def emit(cid, text, acc, detail=""):
    if acc.ng:
        check(cid, text, acc.g, detail)
    if acc.np:
        check(cid + "-ph", "POST HOC, gates nothing (the reproduction point at h = 0.01): " + text, acc.p, detail, ph=True)


# --------------------------------------------------------------------------------------------------------------------
# the predictions, the separation test, the speed
# --------------------------------------------------------------------------------------------------------------------
def predict(kind, g1, rth, rdot, h, mutated=True):
    """Predicted settled mass r_th (M_0 + y): y = a g + b g^2 (P2, #266's relation) or y = a g + b g^2 + c rdot (P3, #267's post hoc fit). `mutated`: whether the mutations act."""
    M0 = CRITR["M0"][h]
    if kind == "P2":
        a_, b_ = CRITR["P2"][h]
        c_ = 0.0
    else:
        a_, b_, c_ = CRITR["P3"]
    mu = mutated and MUT is not None
    a_ *= 1.05 if (mu and mut("coef_a")) else 1.0
    c_ *= 0.8 if (mu and mut("coef_c")) else 1.0
    y = a_ * g1 + b_ * g1 ** 2 + c_ * rdot
    return (1.0 if (mu and mut("no_rscale")) else rth) * (M0 + y)


def fam_b(h):
    """Family B (#266's B = 0, 0.025, 0.05, 0.1, 0.125, stored runs read only) at h: gamma_1, r_th, rdot_th at departure (eps = 1e-3), and the stored plateaus."""
    R = m.load266(h)
    Bs = (0.0, 0.025, 0.05, 0.1, 0.125)
    g1 = np.array([R[B]["mem"]["g1"] for B in Bs])
    g1[0] = 0.0                                     # B = 0 is symmetric: |gamma_1| <= 3e-8 (#266)
    rth, rdot = np.array([R[B]["mem"]["r_th"] for B in Bs]), np.array([R[B]["speed"] for B in Bs])
    o = np.argsort(g1)                              # ascending in gamma_1 (a mutated asymmetry sign reverses it)
    return dict(R=R, Bs=Bs, g1=g1[o], rth=rth[o], rdot=rdot[o])


def sep_test(g1, v, vh, curves, k=None):
    """Separation of a route point's (r_th, rdot_th) at gamma_1 from family B's at the same gamma_1. v: (r_th, rdot) at the finest h, vh: at h = 0.01 (its own h difference).
    Family B is interpolated in gamma_1 over its five points by a piecewise-linear and by a PCHIP curve; the separation is the smaller of the two. Informative when gamma_1 lies inside
    family B's range and the separation exceeds k times the h difference for r_th or for rdot_th."""
    from scipy.interpolate import PchipInterpolator
    k = CRITR["sep_k"] if k is None else k
    k = k * 1e-3 if mut("sep_k") else k
    x = curves["g1"]
    out = dict(inside=bool(x[0] <= g1 <= x[-1]), q={})
    for j, nm in enumerate(("r_th", "rdot")):
        yB = (curves["rth"], curves["rdot"])[j]
        lin = float(np.interp(g1, x, yB))
        pch = float(PchipInterpolator(x, yB)(g1))
        d = min(abs(v[j] - lin), abs(v[j] - pch))
        hd = abs(vh[j] - v[j])
        out["q"][nm] = dict(d_lin=v[j] - lin, d_pchip=v[j] - pch, d=d, hdiff=hd, ok=bool(d > k * hd))
    out["informative"] = out["inside"] and any(q["ok"] for q in out["q"].values())
    return out


# --------------------------------------------------------------------------------------------------------------------
# planted cases
# --------------------------------------------------------------------------------------------------------------------
def planted_traj(t, g0, r0, v):
    """A planted member trajectory: r_th(t) with a smooth minimum near t = 13.5 and a late speed v per unit t, gamma_1(t) falling slowly (the sizes of the measured ones)."""
    x = np.asarray(t, float) - 13.5
    ramp = 0.5 * (x + np.sqrt(x * x + 4.0))
    return r0 + 2.5e-4 * x * x / (1.0 + 0.02 * x * x) + v * ramp, g0 - 6.0e-5 * x


def planted_entry(h, g0, r0, v, e=8):
    """A planted pair: odd part Q0 delta e^{lambda t} r_th on a locus with tau = 0.97 t, profile rows every 0.5 in t carrying the exact static member of the planted (r_th, gamma_1),
    read through the same member_at and speed_at as the runs. Returns the reading and the known values at the planted departure time."""
    delta = 10.0 ** (-e)
    dt = h / SQ2
    t = np.arange(2.0, 34.0, dt)
    rt, g1t = planted_traj(t, g0, r0, v)
    q = m.Q0 * delta * np.exp(m.LAM_T * t) * rt
    s = np.arange(len(t))
    loc = dict(s_a=s, s_b=s, t_a=t, t_b=t, tau_a=0.97 * t, tau_b=0.97 * t, r_a=rt - q, r_b=rt + q)
    rows = []
    r_pts = np.array(tl.STALL["r_pts"], float)
    for tk in np.arange(3.0, 30.0 + 1e-9, 0.5):
        rth_k, g_k = (float(x) for x in planted_traj(np.array(tk), g0, r0, v))
        s_k = g_k / math.sqrt(1 + g_k * g_k)
        Mp, Mm = tl.member_M(ND, rth_k, s_k, r_pts * rth_k)
        rows.append(dict(s=int(round(tk * SQ2 / h)), t=float(tk), r_th=rth_k, r_th_a=rth_k, rr_m=(r_pts * rth_k).tolist(), Mp_m=Mp.tolist(), Mm_m=Mm.tolist(),
                         rr_a=(r_pts * rth_k).tolist(), Mp_a=Mp.tolist(), Mm_a=Mm.tolist()))
    d = dict(e=e, h=h, loc=loc, m=dict(rows=rows))
    mem = b.member_at(d, CRITR["eps_dep"], "m")
    sp = m.speed_at(loc, mem["t_dep"])
    t_true = math.log(CRITR["eps_dep"] / (m.Q0 * delta)) / m.LAM_T
    r_true, g_true = (float(x) for x in planted_traj(np.array(t_true), g0, r0, v))
    eps_ = 1e-5
    v_true = float((planted_traj(np.array(t_true + eps_), g0, r0, v)[0] - planted_traj(np.array(t_true - eps_), g0, r0, v)[0]) / (2 * eps_) / 0.97)
    return dict(mem=mem, sp=sp, t_true=t_true, r_true=r_true, g_true=g_true, v_true=v_true)


PLANTED_ENT = ((0.012, 0.9970, -1.0e-3), (0.048, 0.9990, 2.0e-4), (0.075, 1.0080, 2.0e-3), (0.012, 1.0010, 1.0e-3))      # (gamma_1, r_th at t = 13.5, late speed per unit t); r_th distinct


def planted_pred_data():
    rng = np.random.default_rng(283)
    g = np.array([0.012, 0.048, 0.075, 0.012, 0.048, 0.075])
    rth = 1.0 + rng.uniform(-6e-3, 6e-3, 6)
    rd = rng.uniform(-6e-4, 1.6e-3, 6)
    return g, rth, rd


def check_planted():
    """PMEM, PDEP: #266's planted static members and growing pair through the extraction and the departure time; PTD: planted pairs through member_at and speed_at (departure time,
    r_th, gamma_1, rdot_th at departure against the known values); PPRED: planted data that follow P2 (P3) exactly pass P2 (P3) and fire P3 (P2); PSEP: planted points separated from a planted
    family B by known amounts are told informative or not."""
    errs, rows = [], []
    for g1 in (0.0, 0.012, 0.048, 0.075):
        for h in H2:
            r = b.planted_member_case(h, 13.0, g1)
            errs.append((abs(r["g1_rec"] - g1), abs(r["rth_rec"] - r["rth"]), max(abs(r["mp_rec"] - r["mp"]), abs(r["mm_rec"] - r["mm"]))))
            rows.append(f"g1 {g1:g} h {h:g}: err g1 {errs[-1][0]:.1e}, r_th {errs[-1][1]:.1e}, m+- {errs[-1][2]:.1e}")
    for ln in rows:
        print("    PMEM " + ln)
    eg, er, em = (max(x[i] for x in errs) for i in range(3))
    check("PMEM", f"planted static members (gamma_1 = 0, 0.012, 0.048, 0.075, known r_th, s_th, end masses) as profile slices through slice_masses, fit_asym, the parabolic throat minimum and asym_masses, h = {H2}",
          eg <= CRITR["bar_pm_g"] and er <= CRITR["bar_pm_r"] and em <= CRITR["bar_pm_m"],
          f"max error gamma_1 {eg:.1e} (bar {CRITR['bar_pm_g']:g}), r_th {er:.1e} (bar {CRITR['bar_pm_r']:g}), m_+- {em:.1e} (bar {CRITR['bar_pm_m']:g})")
    lam, q0, h_ = 1.18, 1e-9, 0.01 / SQ2
    t = np.arange(0.0, 30.0, h_)
    q = q0 * np.exp(lam * t)
    rt = 1.0 + 0.003 * np.tanh(t / 5.0)
    worst = 0.0
    for eps in (1e-3, 1e-2, 1e-1):
        worst = max(worst, abs(b.departure_time(t, rt, q * rt, eps) - math.log(eps / q0) / lam))
    check("PDEP", "planted growing pair (odd part q0 e^{lambda t}, known crossing times of 1e-3, 1e-2, 1e-1) through departure_time", worst <= CRITR["bar_dep"], f"max |t_dep - exact| {worst:.1e} (bar {CRITR['bar_dep']:g})")
    worst = [0.0, 0.0, 0.0, 0.0]
    for g0, r0, v in PLANTED_ENT:
        for h in H2:
            en = planted_entry(h, g0, r0, v)
            err = (abs(en["mem"]["t_dep"] - en["t_true"]), abs(en["mem"]["g1"] - en["g_true"]), abs(en["mem"]["r_th"] - en["r_true"]), abs(en["sp"] - en["v_true"]))
            worst = [max(a_, b_) for a_, b_ in zip(worst, err)]
            print(f"    PTD gamma_1 {g0:g} r_th {r0:g} v {v:+.1e} h {h:g}: t_dep {en['mem']['t_dep']:.4f} ({en['t_true']:.4f}), gamma_1 {en['mem']['g1']:.6f} ({en['g_true']:.6f}), "
                  f"r_th {en['mem']['r_th']:.6f} ({en['r_true']:.6f}), rdot_th {en['sp']:+.3e} ({en['v_true']:+.3e})")
    check("PTD", "planted pairs (odd part Q0 delta e^{lambda t}, trajectories of the measured size, rows every 0.5 in t from the exact static member) through member_at and speed_at: departure time, gamma_1, r_th and rdot_th "
          "at departure", worst[0] <= CRITR["bar_ptd_t"] and worst[1] <= CRITR["bar_ptd_g"] and worst[2] <= CRITR["bar_ptd_r"] and worst[3] <= CRITR["bar_pspd"],
          f"max error t_dep {worst[0]:.1e} (bar {CRITR['bar_ptd_t']:g}), gamma_1 {worst[1]:.1e} ({CRITR['bar_ptd_g']:g}), r_th {worst[2]:.1e} ({CRITR['bar_ptd_r']:g}), rdot_th {worst[3]:.1e} ({CRITR['bar_pspd']:g})")
    g, rth, rd = planted_pred_data()
    out = {}
    for src in ("P2", "P3"):
        M = np.array([predict(src, g[i], rth[i], rd[i], H_FIN, mutated=False) for i in range(len(g))])
        for kind in ("P2", "P3"):
            out[(src, kind)] = max(abs(M[i] - predict(kind, g[i], rth[i], rd[i], H_FIN)) for i in range(len(g)))
    ok = out[("P2", "P2")] <= 1e-12 and out[("P3", "P3")] <= 1e-12 and out[("P3", "P2")] > CRITR["bar_pred"] and out[("P2", "P3")] > CRITR["bar_pred"]
    check("PPRED", "planted data that follow P2 (P3) exactly pass the P2 (P3) test and fire the P3 (P2) test",
          ok, "max |M - predicted|: " + ", ".join(f"data {s_} against {k_}: {v_:.1e}" for (s_, k_), v_ in out.items()) + f" (bar {CRITR['bar_pred']:g}; the matching pairs pass at 1e-12)")
    x = np.array([0.0, 0.0029, 0.0118, 0.0483, 0.0755])
    curves = dict(g1=x, rth=1.0 - 0.4 * x + 6.0 * x ** 2, rdot=-2e-4 + 5e-3 * x + 0.4 * x ** 2)
    cases = (("shift r_th 3e-3", 0.0483, (3e-3, 0.0), (1e-6, 1e-5), True), ("shift rdot 5e-4", 0.0483, (0.0, 5e-4), (1e-6, 1e-5), True), ("shift 1e-7, 1e-6", 0.0483, (1e-7, 1e-6), (1e-6, 1e-5), False),
             ("outside the range", 0.09, (3e-3, 5e-4), (1e-6, 1e-5), False), ("shift 3e-6 against h-diff 1e-6", 0.0118, (3e-6, 0.0), (1e-6, 1e-5), False))
    bad, det = [], []
    for nm, g1, sh, hd, want in cases:
        base = (float(np.interp(g1, x, curves["rth"])), float(np.interp(g1, x, curves["rdot"])))
        v = (base[0] + sh[0], base[1] + sh[1])
        vh = (v[0] + hd[0], v[1] + hd[1])
        r = sep_test(g1, v, vh, curves)
        det.append(f"{nm}: informative {r['informative']} (planted {want})")
        if r["informative"] != want:
            bad.append(nm)
    check("PSEP", "planted points at known separations from a planted family B (shifts 3e-3, 5e-4, 1e-7, 3e-6; h differences 1e-6, 1e-5; one outside the range) through the separation test: informative or not as planted",
          not bad, "; ".join(det))


# --------------------------------------------------------------------------------------------------------------------
# the analysis
# --------------------------------------------------------------------------------------------------------------------
def chain_of(e):
    """The stored series and summary, or, under a mutation that acts on the chain, a fresh run from the stored p*_h."""
    if mut("no_settle"):
        keep = e["s"]["xp"] <= e["s"]["xp"][-1] - 3.0
        return {k: (v[keep] if isinstance(v, np.ndarray) and v.ndim == 1 and len(v) == len(keep) else v) for k, v in e["s"].items()}, e["i"]
    if MUT not in RERUN:
        return e["s"], e["i"]
    set_route(FAM, e["B"])
    K = 0 if mut("no_zoom") else (12 if mut("short_limit") else None)
    r = g.run_real_g(e["h"], DELTA, e["p"]["ps"], K=K)
    return r["ser"], dict(status=r["status"], event=r["event"], point=r["point"], up=r["up"], Nc=r["Nc"], F1=r["F1"], K=r["K"], h=e["h"])


def derive(D):
    for e in D.values():
        if "err" in e:
            continue
        try:
            if mut("no_run") and e["name"] == POINTS[0][0] and e["h"] == H_FIN:
                raise RuntimeError("member output treated as missing")
            e["S"], e["I"] = chain_of(e)
            if len(e["S"]["r"]) == 0:
                raise RuntimeError("empty series: " + str(e["I"]["status"]))
            e["SS"] = g.summary_g(e["S"], e["h"])
            e["mem"] = b.member_at(e, CRITR["eps_dep"], "m")
            e["sp"] = m.speed_at(e["loc"], e["mem"]["t_dep"])
            t_ = e["mem"]["t_dep"]
            rms = [x["rms_m"] for x in e["m"]["rows"] if abs(x["t"] - t_) <= 0.5] if np.isfinite(t_) else []
            e["rms"] = max(rms) if rms else float("nan")
        except Exception as ex:
            e["err"] = "derive: " + repr(ex)
    return D


def declared_spec(nm):
    spec = dict((n, s) for n, s, _ in POINTS + (REP,))[nm]
    return (spec[0], spec[1], spec[2] * (1.1 if mut("spec_w") else 1.0))


def analyze_cmd(raw, rn=None):
    t0 = time.time()
    root = os.path.join(raw, "runs")
    pts = list(POINTS) + [REP]
    D = {}
    for nm, spec, tg in pts:
        for h in H2:
            e = load_point(root, spec, h)
            e["name"], e["target"] = nm, tg
            D[(nm, h)] = e
    miss = {f"{k[0]} h{k[1]:g}": e["err"] for k, e in D.items() if "err" in e}
    print(f"== Runs read from {root}: {len(D) - len(miss)} of {len(D)} (point, h) complete" + (f"; missing: {miss}" if miss else ""))
    derive(D)
    for k, e in D.items():
        if "err" in e and f"{k[0]} h{k[1]:g}" not in miss:
            miss[f"{k[0]} h{k[1]:g}"] = e["err"]
            print(f"    {k[0]} h {k[1]:g}: {e['err']}")
    ok = {k: e for k, e in D.items() if "err" not in e}
    bad = {k: (sm.CRIT[k], v) for k, v in g.DECLARED.items() if k != "X" and sm.CRIT[k] != v}
    xg = th.CRIT["X"]["G"][ND]
    spec_bad = [f"{k[0]} h {k[1]:g}" for k, e in D.items() if "p" in e and "i" in e and (list(e["p"]["B"]) != list(declared_spec(k[0])) or list(e["i"]["B"]) != list(declared_spec(k[0])))]
    check("CFG", "settled_mass.CRIT (F1, K, W1, rcut, tol, dxp, dm) equal #41's declared chain, X = 28 for G, the ghost pulse (xc 3, w 0.5), and the (B, xc, w) stored in each run equal the declared one",
          not bad and xg == g.DECLARED["X"] and th.CRIT["pulse"] == GHOST and not spec_bad, f"differs: {bad}; X(G) {xg}; pulse {th.CRIT['pulse']}; specs {spec_bad}" if (bad or xg != 28.0 or spec_bad) else "")
    set_route(FAM, REP[1])
    A_ = (0.0, 0.0149, 0.02)
    check("FAM0", "route_family with the default normal pulse (B, 3.0, 0.5) returns what threshold.family returns (the pulses, symmetry flag and signs), bit for bit",
          all(route_family("G", a_) == _family_orig("G", a_) for a_ in A_), "")
    acc = Acc()
    for k, e in D.items():
        acc.add("err" not in e, *k)
    emit("RUNS", "every (point, h) has its pstar, chain and member output and reads without error", acc, f"missing {sorted(miss)}" if miss else "")
    print("\n== Planted cases")
    check_planted()

    print("\n== p*_h of each point: bisection to hi - lo <= 1e-14 (threshold.bisect on threshold.point, X = 28), the h of the job as passed")
    acc = Acc()
    for (nm, h), e in ok.items():
        p, pts_ = e["p"], e["pts"]
        lo, hi = b.bracket_after(pts_, 30) if mut("bisect_early") else (p["lo"], p["hi"])
        good = (pts_[0]["fate"] == -1 and pts_[1]["fate"] == 1 and p["h_arg"] == repr(float(h)) and b.bracket_after(pts_, len(pts_)) == (p["lo"], p["hi"])
                and ((hi - lo <= CRITR["bar_bracket"] and p["stopped"] is None) or (p["stopped"] == "static through X" and hi - lo <= CRITR["bar_bracket_X"])))
        acc.add(good, nm, h)
        print(f"    {nm} {tuple(e['B'])} h {h:g}: p*_h {p['ps']:.16e}, hi - lo {hi - lo:.1e}, {p['runs']} runs, stopped {p['stopped']}{'' if good else '  <-- NOT OK'}")
    emit("BR", "bracket ends are collapse and expansion, the replayed bisection equals lo and hi, and either hi - lo <= 1e-14 or the bisection stopped by 'static through X' (a lifetime past the rectangle, as p* is approached) at hi - lo <= 2e-14; h passed as its float repr", acc)

    print("\n== Does the near-threshold run stall? Lifetime T at p*_h -/+ 1e-8 against 1e-6 (a stall adds about 4 per side)")
    acc = Acc()
    for (nm, h), e in ok.items():
        sl, dm, dp = b.stalls(e)
        acc.add(sl, nm, h)
        i = e["i"]
        print(f"    {nm} h {h:g}: collapse side T {i['pt6']['T']:.3f} (1e-6) -> {i['point']['T']:.3f} (1e-8), +{dm:.2f}; expansion side {i['up6']['T']:.3f} -> {i['up']['T']:.3f}, +{dp:.2f}: {'stall' if sl else 'NOT a stall'}")
    emit("STL", "each point stalls: collapse at p*_h - 1e-6 and 1e-8 and expansion at + 1e-6 and 1e-8, and T gains >= 2.0 per side from 1e-6 to 1e-8", acc)
    acc = Acc()
    for (nm, h), e in ok.items():
        ev = e["m"]["ev40_b"].get("expand") if mut("evt_side") else e["m"]["ev40_a"].get("collapse")
        i_ = e["I"]
        acc.add(ev is not None and list(i_["event"]) == list(ev) and i_["point"]["fate"] == -1 and i_["up"]["fate"] == 1, nm, h)
    emit("EVT", "the chain's first trapped cell equals the pair's collapse-side one (same amplitude p*_h - 1e-8), fate collapse at p*_h - 1e-8 and expansion at p*_h + 1e-8", acc)

    print("\n== The chain (#41's, unchanged) at delta = 1e-8: limit, first law, settled value")
    Fmax = sm.CRIT["F1"] * 2.0 ** sm.CRIT["K"]
    acc = Acc()
    for (nm, h), e in ok.items():
        s, i = e["S"], e["I"]
        fine = np.allclose(np.diff(s["xp"]), h, atol=1e-9)
        acc.add(i["status"] == "limit" and 0.5 * Fmax <= s["F"][-1] < Fmax and s["edge"].min() >= i["Nc"] / 8 and fine, nm, h)
        print(f"    {nm} h {h:g}: {i['status']}, F_end {s['F'][-1]:.4g}, edge {s['edge'].min():.0f}/{i['Nc']}")
    emit("RC5", f"every run ends at the declared limit (F = e^-f on the horizon, {0.5 * Fmax:g} <= F_end < {Fmax:g}), horizon at least Nc/8 cells from the box edges, no gaps", acc)
    acc, worst = Acc(), 0.0
    for (nm, h), e in ok.items():
        s = e["S"]
        i1 = int(np.argmax(s["F"] >= 20.0))
        dMs = s["M"][i1] - s["M"][0]
        with np.errstate(all="ignore"):
            rel = abs(dMs - g.flux_integral_g(s, 0, i1)) / abs(dMs)
        rel = rel if np.isfinite(rel) else float("inf")
        worst = max(worst, rel)
        acc.add(rel <= CRITR["bar_flux"], nm, h)
    emit("RC4", "first law on the horizon, both fields: M(x+_b) - M(x+_a) = integral of (r^2/2) e^f r_- (phi0_+^2 - phi1_+^2), from the first trapped cell to F = 20", acc, f"bar {CRITR['bar_flux']:g}; worst {worst:.1e}")
    acc = Acc()
    for (nm, h), e in ok.items():
        q = e["SS"]
        acc.add(abs(q["dplat"]) <= CRITR["bar_plat"] and q["flux05"] <= CRITR["bar_flux05"], nm, h)
    emit("RC8", "M_AH has settled at the limit: plateau mean over the last 1.0 differs from the previous 1.0 by <= 1e-4 and the mean horizon flux (both fields) over the last 0.5 is <= 1e-5", acc,
         "plateau change " + ", ".join(f"{e['SS']['dplat']:+.1e}" for e in ok.values()))
    acc, det = Acc(), []
    names = [nm for nm, _, _ in pts]
    for j, nm in enumerate(names):
        if (nm, H2[0]) not in ok or (nm, H2[1]) not in ok:
            acc.add(False, nm, H_FIN)
            continue
        nm2 = names[(j + 1) % len(names)] if mut("rh_swap") else nm
        d = abs(ok[(nm, H2[0])]["SS"]["plat1"] - ok.get((nm2, H2[1]), ok[(nm, H2[1])])["SS"]["plat1"])
        acc.add(d <= CRITR["bar_hpair"], nm, H_FIN)
        det.append(f"{nm}: {d:.1e}")
    emit("RH", "the settled value at two h: |plat(0.01) - plat(0.005)| (plateau mean over the last 1.0 of x+) <= 5e-5", acc, "; ".join(det))

    print("\n== The member at departure (#225's slice_masses and fit_asym on the pair's mean profile; t_dep: the unstable mode's part reaches eps = 1e-3 of r_th) and the throat's speed there (#267's)")
    print("    point  h       p*_h            t_dep    r_th      gamma_1   rdot_th     rms(fit)  M(plateau)")
    for (nm, h), e in ok.items():
        mem = e["mem"]
        print(f"    {nm:5s} {h:<7g} {e['p']['ps']:.10e}  {mem['t_dep']:7.3f}  {mem['r_th']:.6f}  {mem['g1']:.6f}  {e['sp']:+.3e}  {e['rms']:.1e}  {e['SS']['plat1']:.6f}")
    acc, det = Acc(), []
    for nm, spec, tg in POINTS:
        e = ok.get((nm, H_FIN))
        if e is None:
            acc.add(False, nm, H_FIN)
            continue
        tgv = tg * (1.1 if mut("tgt_scale") else 1.0)
        rel = e["mem"]["g1"] / tgv - 1.0
        acc.add(abs(rel) <= CRITR["tgt_rel"], nm, H_FIN)
        det.append(f"{nm}: {rel:+.1e}")
    emit("TGT", f"gamma_1 at departure (h = {H_FIN:g}) is within {CRITR['tgt_rel']:.0%} of its target (0.012, 0.048, 0.075)", acc, "; ".join(det))
    acc, det = Acc(), []
    for nm in names:
        a_, c_ = ok.get((nm, H2[0])), ok.get((nm, H2[1]))
        if a_ is None or c_ is None:
            acc.add(False, nm, H_FIN)
            continue
        dg = abs(a_["mem"]["g1"] / c_["mem"]["g1"] - 1.0)
        dr = abs(a_["mem"]["r_th"] / c_["mem"]["r_th"] - 1.0)
        dv = abs(a_["sp"] - c_["sp"])
        good = dg <= CRITR["bar_hst_rel"] and dr <= CRITR["bar_hst_rel"] and dv <= CRITR["bar_hst_rdot"]
        acc.add(good, nm, H_FIN)
        det.append(f"{nm}: gamma_1 {dg:.1e}, r_th {dr:.1e}, rdot_th {dv:.1e}")
    for ln in det:
        print("    HST " + ln)
    emit("HST", "the member at departure at the two h (0.01, 0.005) agrees to 1e-4 relative in gamma_1 and r_th, and to 2e-5 absolute in rdot_th (the issue's acceptance)", acc)
    acc = Acc()
    for (nm, h), e in ok.items():
        acc.add(np.isfinite(e["rms"]) and e["rms"] * (10.0 if mut("fitq_scale") else 1.0) <= CRITR["bar_fitq"], nm, h)
    emit("FITQ", "the member's fit of M(r) on the two sides has an rms residual of at most 1e-4 at the slices around departure (a static member describes the state there)", acc)

    print("\n== Separation from family B at the same gamma_1 (#266's B = 0 ... 0.125, stored, h = 0.005; piecewise-linear and PCHIP interpolation; the smaller of the two)")
    curves = fam_b(H_FIN)
    print("    family B: gamma_1 " + ", ".join(f"{x:.5f}" for x in curves["g1"]) + "; r_th " + ", ".join(f"{x:.5f}" for x in curves["rth"]) + "; rdot_th " + ", ".join(f"{x:+.2e}" for x in curves["rdot"]))
    acc, informative, SEP = Acc(), {}, {}
    for nm, spec, tg in POINTS:
        a_, c_ = ok.get((nm, H2[0])), ok.get((nm, H2[1]))
        if a_ is None or c_ is None:
            acc.add(False, nm, H_FIN)
            informative[nm] = False
            continue
        cv = curves
        if mut("sep_self"):
            cv = dict(g1=np.array([0.0, c_["mem"]["g1"] - 1e-4, c_["mem"]["g1"], c_["mem"]["g1"] + 1e-4]), rth=np.array([c_["mem"]["r_th"]] * 4), rdot=np.array([c_["sp"]] * 4))
        r_ = sep_test(c_["mem"]["g1"], (c_["mem"]["r_th"], c_["sp"]), (a_["mem"]["r_th"], a_["sp"]), cv)
        SEP[nm] = r_
        informative[nm] = r_["informative"]
        acc.add(r_["informative"], nm, H_FIN)
        q = r_["q"]
        print(f"    {nm}: gamma_1 {c_['mem']['g1']:.5f} (inside family B's range: {r_['inside']}); r_th - family B: {q['r_th']['d_lin']:+.2e} (linear), {q['r_th']['d_pchip']:+.2e} (PCHIP), its h difference {q['r_th']['hdiff']:.1e}"
              f" -> {'ok' if q['r_th']['ok'] else 'no'}; rdot_th - family B: {q['rdot']['d_lin']:+.2e}, {q['rdot']['d_pchip']:+.2e}, h difference {q['rdot']['hdiff']:.1e} -> {'ok' if q['rdot']['ok'] else 'no'}"
              f": {'informative' if r_['informative'] else 'NOT informative'}")
    emit("SEP", "the route separates the parameters: each point's r_th or rdot_th at departure differs from family B's at the same gamma_1 by more than 5 times its own h difference (a point that fails is not informative and does not count for P2, P3)",
         acc)

    print("\n== The two predictions on the new members: M - predicted, h = 0.005 (gating) and 0.01 (P2 only, post hoc)")
    print("    P2: M = r_th (M_0 + a g + b g^2), M_0 = 0.223086, (a, b) = (0.2027, 0.8534);  P3: M = r_th (M_0 + a g + b g^2 + c rdot_th), (a, b, c) = (0.2300, 0.2724, 0.7835)")
    print("    point   informative   M(plateau)   P2 predicted  M - P2      P3 predicted  M - P3      | h = 0.01: M - P2 (P2's own (a, b) at 0.01)")
    res = {"P2": {}, "P3": {}}
    for nm, spec, tg in POINTS:
        e, e1 = ok.get((nm, H_FIN)), ok.get((nm, H2[0]))
        if e is None or e1 is None:
            continue
        M, g1, rth, rd = e["SS"]["plat1"], e["mem"]["g1"], e["mem"]["r_th"], e["sp"]
        p2, p3 = predict("P2", g1, rth, rd, H_FIN), predict("P3", g1, rth, rd, H_FIN)
        res["P2"][nm], res["P3"][nm] = M - p2, M - p3
        M1 = e1["SS"]["plat1"]
        res.setdefault("P2h", {})[nm] = M1 - predict("P2", e1["mem"]["g1"], e1["mem"]["r_th"], e1["sp"], H2[0])
        print(f"    {nm:6s}  {str(informative.get(nm)):10s}  {M:.6f}   {p2:.6f}      {M - p2:+.2e}   {p3:.6f}      {M - p3:+.2e}   | {res['P2h'][nm]:+.2e}")
    inf = [nm for nm, _, _ in POINTS if informative.get(nm) and nm in res["P2"]]
    for kind in ("P2", "P3"):
        text = (f"{kind}: M_AH settled at the new member equals r_th (M_0 + " + ("a g + b g^2" if kind == "P2" else "a g + b g^2 + c rdot_th") + f") within {CRITR['bar_pred']:g} for every informative point (h = {H_FIN:g}; coefficients "
                "fixed by the issue before the runs)")
        if not inf:
            check(kind, text, None, "no informative point")
        else:
            worst = max(abs(res[kind][nm]) for nm in inf)
            check(kind, text, worst <= CRITR["bar_pred"], f"{len(inf)} informative points ({', '.join(inf)}); max |M - predicted| {worst:.1e} (bar {CRITR['bar_pred']:g}); " +
                  ", ".join(f"{nm} {res[kind][nm]:+.1e}" for nm in inf))
    if inf:
        worst = max(abs(res["P2h"][nm]) for nm in inf)
        check("P2-ph", f"POST HOC, gates nothing: P2 at h = 0.01 with #266's (a, b) at that h", worst <= CRITR["bar_pred"], f"max |M - predicted| {worst:.1e}", ph=True)

    print("\n== Post hoc (gates nothing): where P2's miss goes, and what P3 leaves (h = 0.005)")
    print("    point   M - P2      speed term r_th c rdot_th   (a, b) shift term   M - P3      r_th - family B's (linear)   (M - P3) / that")
    xs_, ys_ = [], []
    for nm, spec, tg in POINTS:
        e = ok.get((nm, H_FIN))
        if e is None or nm not in SEP:
            continue
        g1, rth, rd = e["mem"]["g1"], e["mem"]["r_th"], e["sp"]
        a3, b3, c3 = CRITR["P3"]
        a2, b2 = CRITR["P2"][H_FIN]
        drt = SEP[nm]["q"]["r_th"]["d_lin"]
        r3 = e["SS"]["plat1"] - predict("P3", g1, rth, rd, H_FIN)
        xs_.append(drt)
        ys_.append(r3)
        print(f"    {nm:6s}  {e['SS']['plat1'] - predict('P2', g1, rth, rd, H_FIN):+.2e}   {rth * c3 * rd:+.2e}                  {rth * ((a3 - a2) * g1 + (b3 - b2) * g1 ** 2):+.2e}          {r3:+.2e}   {drt:+.2e}                  {(r3 / drt if abs(drt) > 1e-12 else float('nan')):+.2e}")
    if xs_:
        xs_, ys_ = np.array(xs_), np.array(ys_)
        k_ = float(xs_ @ ys_ / max(xs_ @ xs_, 1e-300))
        print(f"    least squares through the origin: (M - P3) = {k_:+.3e} (r_th - family B's); residual max {np.max(np.abs(ys_ - k_ * xs_)):.1e} (six points, one parameter: a description, not a test)")

    print("\n== Reproduction of #266's B = 0.1 point (the default pulse through this driver; stored runs of #266 read only)")
    acc, det = Acc(), []
    for h in H2:
        e = ok.get(("REP", h))
        if e is None:
            acc.add(False, "REP", h)
            continue
        R = (curves if h == H_FIN else fam_b(h))["R"]
        ref = R[0.125] if mut("rep_off") else R[0.1]
        p266 = json.load(open(f"{I266}/runs/G_B{0.125 if mut('rep_off') else 0.1:g}/h{h:g}/pstar/summary.json"))["ps"]
        dM = e["SS"]["plat1"] - ref["plat"]
        dm_ = (e["mem"]["g1"] - ref["mem"]["g1"], e["mem"]["r_th"] - ref["mem"]["r_th"], e["sp"] - ref["speed"])
        dp = e["p"]["ps"] - p266
        acc.add(abs(dM) <= CRITR["bar_rep_M"] and max(abs(x) for x in dm_) <= CRITR["bar_rep_mem"] and abs(dp) <= CRITR["bar_rep_p"], "REP", h)
        det.append(f"h {h:g}: dM {dM:+.1e}, dgamma_1 {dm_[0]:+.1e}, dr_th {dm_[1]:+.1e}, drdot_th {dm_[2]:+.1e}, dp* {dp:+.1e}")
    for ln in det:
        print("    REP " + ln)
    emit("REP", "the point (B, xc, w) = (0.1, 3.0, 0.5) reproduces #266's B = 0.1 (p*_h to 1e-13, M and the member (gamma_1, r_th, rdot_th) to 1e-6)", acc)

    if rn:
        print(f"\n== Post hoc (gates nothing): a second roundoff realization of the six points at h = {H_FIN:g} (p*_h bisected from another bracket: the same p*_h to within its own floor, a different ulp-noise realization of the chain), {rn}")
        R2 = {}
        for nm, spec, tg in POINTS:
            e2 = load_point(os.path.join(rn, "runs"), spec, H_FIN)
            e2["name"] = nm
            R2[nm] = e2
        derive({(nm, H_FIN): e2 for nm, e2 in R2.items()})
        print("    point  p*_h difference  M(rn) - M    gamma_1, r_th, rdot_th differences      M - P2 (rn)   M - P3 (rn)")
        for nm, spec, tg in POINTS:
            e2, e = R2[nm], ok.get((nm, H_FIN))
            if "err" in e2 or e is None:
                print(f"    {nm}: not available ({e2.get('err')})")
                continue
            M2, g2, r2, d2 = e2["SS"]["plat1"], e2["mem"]["g1"], e2["mem"]["r_th"], e2["sp"]
            print(f"    {nm:5s}  {e2['p']['ps'] - e['p']['ps']:+.1e}       {M2 - e['SS']['plat1']:+.1e}    {g2 - e['mem']['g1']:+.1e}, {r2 - e['mem']['r_th']:+.1e}, {d2 - e['sp']:+.1e}        "
                  f"{M2 - predict('P2', g2, r2, d2, H_FIN):+.2e}    {M2 - predict('P3', g2, r2, d2, H_FIN):+.2e}")
    print("\n== Post hoc (gates nothing): a refit with the new points, and the decomposition of the residuals")
    try:
        R = curves["R"]
        allp = [(f"B {B_:g}", R[B_]["mem"]["g1"], R[B_]["mem"]["r_th"], R[B_]["speed"], R[B_]["plat"]) for B_ in (0.025, 0.05, 0.1, 0.125)]
        allp += [(nm, ok[(nm, H_FIN)]["mem"]["g1"], ok[(nm, H_FIN)]["mem"]["r_th"], ok[(nm, H_FIN)]["sp"], ok[(nm, H_FIN)]["SS"]["plat1"]) for nm, _, _ in POINTS if (nm, H_FIN) in ok]
        c0 = R[0.0]["plat"]
        yv = np.array([p_[4] / p_[2] - c0 for p_ in allp])
        gv = np.array([p_[1] for p_ in allp])
        sv = np.array([p_[3] for p_ in allp])
        print("    y = M / r_th - M_0 at each point (the four B of #266 and the six new): " + ", ".join(f"{p_[0]} {y_:+.6f}" for p_, y_ in zip(allp, yv)))
        for nm_, A_ in (("y = a g + b g^2", np.stack([gv, gv ** 2], 1)), ("y = a g + b g^2 + c rdot", np.stack([gv, gv ** 2, sv], 1))):
            loo = []
            for k_ in range(len(yv)):
                idx = [j for j in range(len(yv)) if j != k_]
                loo.append(A_[k_] @ np.linalg.lstsq(A_[idx], yv[idx], rcond=None)[0] - yv[k_])
            co = np.linalg.lstsq(A_, yv, rcond=None)[0]
            print(f"    fit over these {len(yv)} points, {nm_}: coefficients {', '.join(f'{x:+.4f}' for x in co)}; residual max {np.max(np.abs(A_ @ co - yv)):.1e}; leave-one-out errors "
                  + ", ".join(f"{x:+.1e}" for x in loo))
        print("    *post hoc: fitted on points that include the ones tested; descriptions, not tests*")
    except Exception as ex:
        print(f"    post hoc refit not computable: {ex!r}")
    return finish(t0)


# --------------------------------------------------------------------------------------------------------------------
# mutations: name -> the checks that must fire (computed on the development set, raw I283-mac dev/, before being declared in its README)
# --------------------------------------------------------------------------------------------------------------------
MUTATIONS = {      # name: the checks it newly fails (those failing without a mutation, --baseline, excluded); computed on the development set (raw I283-mac dev/mut_a) before being declared in its README;
                   # HST, P2, P3 and P2-ph fail without a mutation on that set, so no mutation can newly fire them
    "asym_sign": ['PMEM', 'PTD', 'TGT'],
    "bisect_early": ['BR', 'BR-ph'],
    "coef_a": ['PPRED'],
    "coef_c": ['PPRED'],
    "crit_dm": ['CFG'],
    "dep_eps": ['FITQ', 'PDEP', 'PTD'],
    "eps_h": ['FITQ', 'PTD', 'SEP'],
    "evt_side": ['EVT', 'EVT-ph'],
    "fitq_scale": ['FITQ'],
    "flux_sign": ['RC4', 'RC4-ph'],
    "mass_exponent": ['RC4', 'RC4-ph', 'REP', 'REP-ph'],
    "no_rscale": ['PPRED'],
    "no_run": ['RH', 'RUNS', 'SEP', 'TGT'],
    "no_settle": ['RC5', 'RC5-ph', 'RC8', 'RC8-ph', 'REP', 'REP-ph'],
    "normal_drop": [],
    "rdot_nolapse": ['PTD'],
    "rep_off": ['REP', 'REP-ph'],
    "rh_swap": ['RH'],
    "sep_k": ['PSEP'],
    "sep_self": ['SEP'],
    "short_limit": ['RC5', 'RC5-ph', 'RC8', 'RC8-ph', 'REP', 'REP-ph'],
    "spec_w": ['CFG'],
    "stl_flat": ['STL', 'STL-ph'],
    "tgt_scale": ['TGT'],
    "w_scale": ['FAM0', 'REP', 'REP-ph', 'RH', 'RUNS', 'RUNS-ph', 'SEP', 'TGT'],
    "xc_offset": ['FAM0', 'REP', 'REP-ph', 'RH', 'RUNS', 'RUNS-ph', 'SEP', 'TGT'],
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
        if "baseline" in opts:
            BASE.update(ln.split("]", 1)[1].split(":", 1)[0].strip() for ln in open(opts["baseline"]) if ln.startswith("[FAIL]"))
        b.MUT = m.MUT = g.MUT = sm.MUT = MUT
        dn.MUT = MUT if MUT in ("theta_coef", "ghost_sign", "scalar_coef") else None
        tl.LMUT = MUT if MUT in ("asym_sign",) else None
        if MUT == "crit_dm":
            sm.CRIT["dm"] = 2.0
    import platform
    import numba as nb
    print(f"settled_mass_route.py {cmd} {opts}: numba {nb.__version__}, numpy {np.__version__}, python {platform.python_version()}, {platform.node()}, "
          f"mutation {MUT}, {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    if cmd == "tag":
        print(tag_of(FAM, parse_route(opts["route"])))
    elif cmd == "bracket":
        lo, hi, rows = bracket_scan(parse_route(opts["route"]), float(opts["h"]))
        print(json.dumps(dict(lo=lo, hi=hi, rows=rows)))
    elif cmd == "pstar":
        pstar_cmd(parse_route(opts["route"]), float(opts["h"]), opts["out"], tuple(float(x) for x in opts["bracket"].split(",")))
    elif cmd == "job":
        job_cmd(parse_route(opts["route"]), float(opts["h"]), opts["out"], opts["pstar"])
    elif cmd == "member":
        member_cmd(parse_route(opts["route"]), float(opts["h"]), opts["out"], opts["pstar"])
    elif cmd == "planted":
        t0 = time.time()
        check_planted()
        sys.exit(finish(t0))
    elif cmd == "analyze":
        sys.exit(analyze_cmd(opts["raw"], opts.get("rn")))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()

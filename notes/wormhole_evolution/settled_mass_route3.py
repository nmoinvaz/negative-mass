#!/usr/bin/env python3
"""A frozen refit of family G's settled-mass relation, tested blind on a third route to the member (#295). D = 4, delta = 1e-8.

#266 reads the settled mass of family G along the normal amplitude B, #267 along delta, #283 along two other routes (the normal pulse's width w = 1.0, its centre 4.0). This driver
(1) reads the fifteen members at h = 0.005 that exist (the four B of #266, the five delta of #267, the six of #283) with #283's code (settled_mass_route.py, unchanged),
(2) fits y = M / r_th - M_0 = a gamma_1 + b gamma_1^2 + c rdot_th by plain least squares with M_0 = 0.223086 (the coordinator's ruling) and prints the coefficients and the
leave-one-out errors (`refit`), which are copied into the raw README before any run of the third route (they are the constants of `CRIT3["P4"]`, and the analysis recomputes them from the
stored points and compares), and
(3) reaches gamma_1 = 0.012, 0.048, 0.075 by a third route (the normal pulse's width and centre both changed), runs the chain exactly as #283, and tests the frozen refit on the new members.

A thin driver. It imports settled_mass_route.py (#283), which imports settled_mass_b.py (#266) and settled_mass_member.py (#267), which import settled_mass.py (#41), settled_mass_g.py (#253),
threshold.py (#40), throat_locus.py (#225), dn.py (#39) and throat_drift.py (#255), all unchanged.

    python3 settled_mass_route3.py refit   [--raw=DIR]                               the fifteen points and the refit (reads the stored runs of #266, #267, #283)
    python3 settled_mass_route3.py tag|bracket|pstar|job|member ...                  as settled_mass_route.py (#283), for the third route's points
    python3 settled_mass_route3.py planted [--mutate=NAME]                           planted cases
    python3 settled_mass_route3.py analyze --raw=DIR [--mutate=NAME] [--baseline=FILE] [--rn=DIR]

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
import settled_mass_route as r         # noqa: E402  (#283's driver; imports #266's and #267's drivers and everything below them unchanged)
b, m, g, th, tl, sm, dn, td = r.b, r.m, r.g, r.th, r.tl, r.sm, r.dn, r.td

SQ2 = math.sqrt(2.0)
MUT = None
ND = 2
FAM = "G"
H2 = r.H2                                       # (0.01, 0.005)
H_FIN = r.H_FIN
RAW_ROOT = os.path.expanduser("~/Source/negative-mass-raw")
I266 = m.I266
I267 = os.path.join(RAW_ROOT, "I267-mac")
I283 = os.path.join(RAW_ROOT, "I283-mac")
M0 = 0.223086                                   # B = 0's plateau at h = 0.005 (#266), the coordinator's ruling for the refit
EPS_DEP = 1e-3


def mut(name):
    return MUT == name


# --------------------------------------------------------------------------------------------------------------------
# the fifteen points
# --------------------------------------------------------------------------------------------------------------------
B_266 = (0.025, 0.05, 0.1, 0.125)               # #266's four B (B = 0 is the symmetric member, M_0)
E_267 = (7, 9, 10, 11, 12)                      # #267's five delta = 1e-e at B = 0.1 (e = 8 is #266's B = 0.1)
POINTS_283 = (("W1", (0.07905, 3.0, 1.0)), ("W2", (0.1543, 3.0, 1.0)), ("W3", (0.1909, 3.0, 1.0)),
              ("C1", (0.0393, 4.0, 0.5)), ("C2", (0.07726, 4.0, 0.5)), ("C3", (0.09576, 4.0, 0.5)))


def read_267(e, h=H_FIN):
    """#267's run at delta = 10^-e, B = 0.1 (stored, read only): plateau, member at departure (eps = 1e-3), throat speed."""
    base = os.path.join(I267, "runs", m.TAG, f"h{h:g}")
    d = dict(e=e, h=h, base=f"{base}/e{e:g}")
    d["s"] = dict(np.load(d["base"] + "/chain/series.npz"))
    d["m"] = json.load(open(d["base"] + "/member/member.json"))
    d["loc"] = dict(np.load(d["base"] + "/member/locus.npz"))
    mem = b.member_at(d, EPS_DEP, "m")
    return dict(M=g.summary_g(d["s"], h)["plat1"], g1=mem["g1"], rth=mem["r_th"], rdot=m.speed_at(d["loc"], mem["t_dep"]), t_dep=mem["t_dep"])


def points15(h=H_FIN):
    """The fifteen points at h: name, route, gamma_1, r_th, rdot_th, M (plateau). Route: B (#266), D (#267's delta scan), W, C (#283)."""
    Rb = m.load266(h)
    out = []
    for B in B_266:
        x = Rb[B]
        out.append(dict(name=f"B{B:g}", route="B", g1=x["mem"]["g1"], rth=x["mem"]["r_th"], rdot=x["speed"], M=x["plat"]))
    for e in E_267:
        x = read_267(e, h)
        out.append(dict(name=f"e{e}", route="D", g1=x["g1"], rth=x["rth"], rdot=x["rdot"], M=x["M"]))
    for nm, spec in POINTS_283:
        e = r.load_point(os.path.join(I283, "runs"), spec, h)
        x = r.read_point(e)
        out.append(dict(name=nm, route=nm[0], g1=x["mem"]["g1"], rth=x["mem"]["r_th"], rdot=x["sp"], M=x["M"]))
    return out, Rb[0.0]["plat"]


# --------------------------------------------------------------------------------------------------------------------
# the fit: plain least squares, M_0 fixed, y = M / r_th - M_0 = a g + b g^2 (+ c rdot)
# --------------------------------------------------------------------------------------------------------------------
def design(g1, rdot, form):
    g1, rdot = np.asarray(g1, float), np.asarray(rdot, float)
    cols = [g1, g1 ** 2] + ([rdot] if form == "P4" else [])
    return np.stack(cols, 1)


def fit(g1, rth, rdot, M, form="P4", m0=M0):
    """Least-squares coefficients of y = M / r_th - m0 on the form's columns (no intercept, no weights)."""
    y = np.asarray(M, float) / np.asarray(rth, float) - m0
    return np.linalg.lstsq(design(g1, rdot, form), y, rcond=None)[0]


def loo(g1, rth, rdot, M, form="P4", m0=M0):
    """Leave-one-out: fit without point k, predict point k. Returns (M - predicted) in M and (y - predicted) in y, per point."""
    g1, rth, rdot, M = (np.asarray(x, float) for x in (g1, rth, rdot, M))
    A = design(g1, rdot, form)
    y = M / rth - m0
    eM, ey = [], []
    for k in range(len(y)):
        idx = [j for j in range(len(y)) if j != k]
        co = np.linalg.lstsq(A[idx], y[idx], rcond=None)[0]
        ey.append(y[k] - A[k] @ co)
        eM.append(M[k] - rth[k] * (m0 + A[k] @ co))
    return np.array(eM), np.array(ey)


def loo_hat(g1, rth, rdot, M, form="P4", m0=M0):
    """The closed form of the leave-one-out error of a linear least-squares fit: e_k / (1 - H_kk), H the hat matrix (an independent route to loo(), in y)."""
    A = design(g1, rdot, form)
    y = np.asarray(M, float) / np.asarray(rth, float) - m0
    H = A @ np.linalg.pinv(A)
    return (y - H @ y) / (1.0 - np.diag(H))


def refit_cmd():
    pts, c0 = points15()
    print(f"== The fifteen points at h = {H_FIN:g} (B = 0's plateau {c0:.7f}; M_0 = {M0} is the ruling)")
    print("    point  route  gamma_1    r_th       rdot_th      M            y = M / r_th - M_0")
    g1, rth, rdot, M = (np.array([p[k] for p in pts]) for k in ("g1", "rth", "rdot", "M"))
    for p in pts:
        print(f"    {p['name']:6s} {p['route']:5s}  {p['g1']:.6f}  {p['rth']:.6f}  {p['rdot']:+.4e}  {p['M']:.8f}  {p['M'] / p['rth'] - M0:+.8f}")
    for form in ("P4", "P2"):
        co = fit(g1, rth, rdot, M, form)
        A = design(g1, rdot, form)
        y = M / rth - M0
        res = rth * (M0 + A @ co) - M
        eM, ey = loo(g1, rth, rdot, M, form)
        eh = loo_hat(g1, rth, rdot, M, form)
        print(f"\n   form {form}: y = " + (" + ".join(f"{c:+.10f} {n}" for c, n in zip(co, ("g", "g^2", "rdot")))))
        print(f"     in-sample M - predicted:  " + ", ".join(f"{p['name']} {-x:+.2e}" for p, x in zip(pts, res)) + f";  max {np.max(np.abs(res)):.2e}")
        print(f"     leave-one-out M - predicted: " + ", ".join(f"{p['name']} {x:+.2e}" for p, x in zip(pts, eM)) + f";  max {np.max(np.abs(eM)):.2e}")
        print(f"     leave-one-out y errors:     " + ", ".join(f"{p['name']} {x:+.2e}" for p, x in zip(pts, ey)) + f";  hat-matrix form agrees to {np.max(np.abs(ey - eh)):.1e}")
        print(f"     singular values of the design {', '.join(f'{x:.3e}' for x in np.linalg.svd(A, compute_uv=False))}")


# --------------------------------------------------------------------------------------------------------------------
# the declared set (the raw README, written before the runs of record)
# --------------------------------------------------------------------------------------------------------------------
POINTS3 = (("N1", (0.05648, 3.5, 0.75), 0.012), ("N2", (0.1101, 3.5, 0.75), 0.048), ("N3", (0.136, 3.5, 0.75), 0.075))      # (name, (B, xc, w), target gamma_1): the third route, selected in development (README stage 2)
REPW3 = ("RW3", (0.1909, 3.0, 1.0))          # #283's W3, reproduced at h = 0.005 from its stored bracket
TARGETS = (0.012, 0.048, 0.075)
CRIT3 = dict(
    P4=(0.2296318924, 0.2834170383, 0.7657657062),                     # the frozen refit (raw README stage 1, 305e5fee), h = 0.005, M_0 = 0.223086
    P3=(0.2300, 0.2724, 0.7835),                                       # #267's post hoc fit over nine points (comparison only)
    P2p=(0.2177664539, 0.5888595816),                                  # the refit without the speed over the same fifteen points (comparison only)
    loo_max=7.90e-5, ins_max=5.51e-5,                                  # frozen errors of the refit (M - predicted)
    bar_frz=1e-9, bar_frz_err=1e-7,
    bar_pred=6e-5, sep_k=5.0, sep_extrap=0.03, tgt_rel=0.03,
    M0_h={0.005: 0.223086, 0.01: 0.223095},
)
results = []
BASE = set()


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


class Acc:
    def __init__(self):
        self.ok, self.n = True, 0

    def add(self, good):
        self.ok, self.n = self.ok and bool(good), self.n + 1


# --------------------------------------------------------------------------------------------------------------------
# predictions and the frozen constants
# --------------------------------------------------------------------------------------------------------------------
def predict(kind, g1, rth, rdot, h=H_FIN, mutated=True):
    """M_pred = r_th (M_0 + a g + b g^2 + c rdot). kind: P4 (the frozen refit), P3 (#267's), P2p (the refit without the speed). `mutated`: whether the mutations act."""
    m0 = CRIT3["M0_h"][h]
    co = CRIT3[kind]
    a_, b_ = co[0], co[1]
    c_ = co[2] if len(co) > 2 else 0.0
    mu = mutated and MUT is not None
    a_ *= 1.05 if (mu and mut("coef_a")) else 1.0
    c_ *= 0.8 if (mu and mut("coef_c")) else 1.0
    m0 = 0.223095 if (mu and mut("frz_m0") and h == H_FIN) else m0
    return (1.0 if (mu and mut("no_rscale")) else rth) * (m0 + a_ * g1 + b_ * g1 ** 2 + c_ * rdot)


def frozen_check(pts, c0):
    """FRZ: the refit recomputed from the stored fifteen points equals the constants of CRIT3 (the README's frozen numbers), M_0 equals B = 0's plateau, and its errors equal the frozen maxima."""
    P = list(pts)
    if mut("drop_point"):
        P = [p for p in P if p["name"] != "W3"]
    g1, rth, rdot, M = (np.array([p[k] for p in P]) for k in ("g1", "rth", "rdot", "M"))
    m0 = 0.223095 if mut("frz_m0") else M0
    co = fit(g1, rth, rdot, M, "P4", m0)
    eM, _ = loo(g1, rth, rdot, M, "P4", m0)
    if mut("loo_off"):
        eM = M - rth * (m0 + design(g1, rdot, "P4") @ co)
    ins = M - rth * (m0 + design(g1, rdot, "P4") @ co)
    frozen = np.array(CRIT3["P4"]) * (np.array([1.0 + 1e-3, 1.0, 1.0]) if mut("frz_a") else 1.0)
    dco = np.max(np.abs(co - frozen))
    ok = (dco <= CRIT3["bar_frz"] and abs(c0 - m0) <= 1e-6 and abs(np.max(np.abs(eM)) - CRIT3["loo_max"]) <= 1e-7 and abs(np.max(np.abs(ins)) - CRIT3["ins_max"]) <= 1e-7 and len(P) == 15)
    check("FRZ", "the refit recomputed from the stored fifteen points equals the frozen (a, b, c) to 1e-9, M_0 = B = 0's plateau to 1e-6, and the largest leave-one-out and in-sample errors equal the frozen 7.90e-5 and 5.51e-5 to 1e-7",
          ok, f"{len(P)} points; max |co - frozen| {dco:.1e}; |B = 0 plateau - M_0| {abs(c0 - m0):.1e}; max LOO {np.max(np.abs(eM)):.3e}, max in-sample {np.max(np.abs(ins)):.3e}")
    return co


# --------------------------------------------------------------------------------------------------------------------
# the separation test against a route's curve (extrapolation allowed within 3 % of the end gamma_1)
# --------------------------------------------------------------------------------------------------------------------
def sep3(g1, v, vh, cv, k=None):
    """Separation of a new point's (r_th, rdot_th) at gamma_1 from a reference curve cv (gamma_1 ascending; piecewise-linear and PCHIP, the smaller difference), against k times the
    point's own h difference. gamma_1 beyond the curve's last point is extrapolated (linearly from the last two points; PCHIP's own extension) and counts only up to 3 % beyond it."""
    from scipy.interpolate import PchipInterpolator
    k = CRIT3["sep_k"] if k is None else k
    k = k * 1e-3 if mut("sep_k") else k
    x = np.asarray(cv["g1"], float)
    ext = max(0.0, (g1 - x[-1]) / x[-1])
    out = dict(ext=ext, inside=bool(x[0] <= g1 and ext <= CRIT3["sep_extrap"]), q={})
    for j, nm in enumerate(("r_th", "rdot")):
        y = np.asarray((cv["rth"], cv["rdot"])[j], float)
        lin = float(np.interp(g1, x, y)) if g1 <= x[-1] else float(y[-1] + (g1 - x[-1]) * (y[-1] - y[-2]) / (x[-1] - x[-2]))
        pch = float(PchipInterpolator(x, y, extrapolate=True)(g1))
        d = min(abs(v[j] - lin), abs(v[j] - pch))
        hd = abs(vh[j] - v[j])
        out["q"][nm] = dict(d_lin=v[j] - lin, d_pchip=v[j] - pch, d=d, hdiff=hd, ok=bool(d > k * hd))
    out["informative"] = out["inside"] and any(q["ok"] for q in out["q"].values())
    return out


def route_curve(h, names):
    """An earlier route's members at h (#283's stored runs, B = 0 as the anchor at gamma_1 = 0): gamma_1, r_th, rdot_th, ascending in gamma_1."""
    anchor = r.fam_b(h)
    g1, rth, rdot = [float(anchor["g1"][0])], [float(anchor["rth"][0])], [float(anchor["rdot"][0])]
    for nm, spec in POINTS_283:
        if nm in names:
            x = r.read_point(r.load_point(os.path.join(I283, "runs"), spec, h))
            g1.append(x["mem"]["g1"]), rth.append(x["mem"]["r_th"]), rdot.append(x["sp"])
    o = np.argsort(g1)
    return dict(g1=np.array(g1)[o], rth=np.array(rth)[o], rdot=np.array(rdot)[o])


# --------------------------------------------------------------------------------------------------------------------
# planted cases
# --------------------------------------------------------------------------------------------------------------------
def planted_fit_pts():
    """Planted y on the real fifteen points' (gamma_1, r_th, rdot_th): from the data of `points15` (stored members), with the frozen coefficients."""
    pts, _ = points15()
    return pts


def check_planted():
    """PMEM, PDEP, PTD: #283's extraction cases through this analysis's code path. PFIT: planted data of known (a, b, c) through `fit` and `loo`. PPRED4: planted data that follow the frozen refit
    pass it and fire it when a, or c, is wrong. PSEP3: planted points at known separations from three planted route curves."""
    errs, rows = [], []
    for g1 in (0.0, 0.012, 0.048, 0.075):
        for h in H2:
            q = b.planted_member_case(h, 13.0, g1)
            errs.append((abs(q["g1_rec"] - g1), abs(q["rth_rec"] - q["rth"]), max(abs(q["mp_rec"] - q["mp"]), abs(q["mm_rec"] - q["mm"]))))
            rows.append(f"g1 {g1:g} h {h:g}: err g1 {errs[-1][0]:.1e}, r_th {errs[-1][1]:.1e}, m+- {errs[-1][2]:.1e}")
    for ln in rows:
        print("    PMEM " + ln)
    eg, er, em = (max(x[i] for x in errs) for i in range(3))
    C = r.CRITR
    check("PMEM", f"planted static members (gamma_1 = 0, 0.012, 0.048, 0.075, known r_th, s_th, end masses) as profile slices through slice_masses, fit_asym, the parabolic throat minimum and asym_masses, h = {H2}",
          eg <= C["bar_pm_g"] and er <= C["bar_pm_r"] and em <= C["bar_pm_m"], f"max error gamma_1 {eg:.1e} (bar {C['bar_pm_g']:g}), r_th {er:.1e} (bar {C['bar_pm_r']:g}), m_+- {em:.1e} (bar {C['bar_pm_m']:g})")
    lam, q0, h_ = 1.18, 1e-9, 0.01 / SQ2
    t = np.arange(0.0, 30.0, h_)
    rt = 1.0 + 0.003 * np.tanh(t / 5.0)
    worst = max(abs(b.departure_time(t, rt, q0 * np.exp(lam * t) * rt, eps) - math.log(eps / q0) / lam) for eps in (1e-3, 1e-2, 1e-1))
    check("PDEP", "planted growing pair (odd part q0 e^{lambda t}, known crossing times of 1e-3, 1e-2, 1e-1) through departure_time", worst <= C["bar_dep"], f"max |t_dep - exact| {worst:.1e} (bar {C['bar_dep']:g})")
    worst = [0.0, 0.0, 0.0, 0.0]
    for g0, r0, v in r.PLANTED_ENT:
        for h in H2:
            en = r.planted_entry(h, g0, r0, v)
            err = (abs(en["mem"]["t_dep"] - en["t_true"]), abs(en["mem"]["g1"] - en["g_true"]), abs(en["mem"]["r_th"] - en["r_true"]), abs(en["sp"] - en["v_true"]))
            worst = [max(a_, b_) for a_, b_ in zip(worst, err)]
            print(f"    PTD gamma_1 {g0:g} r_th {r0:g} v {v:+.1e} h {h:g}: t_dep {en['mem']['t_dep']:.4f} ({en['t_true']:.4f}), gamma_1 {en['mem']['g1']:.6f} ({en['g_true']:.6f}), "
                  f"r_th {en['mem']['r_th']:.6f} ({en['r_true']:.6f}), rdot_th {en['sp']:+.3e} ({en['v_true']:+.3e})")
    check("PTD", "planted pairs (odd part, trajectories of the measured size, rows every 0.5 in t from the exact static member) through member_at and speed_at: departure time, gamma_1, r_th and rdot_th at departure",
          worst[0] <= C["bar_ptd_t"] and worst[1] <= C["bar_ptd_g"] and worst[2] <= C["bar_ptd_r"] and worst[3] <= C["bar_pspd"],
          f"max error t_dep {worst[0]:.1e} (bar {C['bar_ptd_t']:g}), gamma_1 {worst[1]:.1e} ({C['bar_ptd_g']:g}), r_th {worst[2]:.1e} ({C['bar_ptd_r']:g}), rdot_th {worst[3]:.1e} ({C['bar_pspd']:g})")
    # PFIT: planted y = P4 exactly on the fifteen real (gamma_1, r_th, rdot_th); one point then displaced by a known amount
    pts = planted_fit_pts()
    g1, rth, rdot = (np.array([p[k] for p in pts]) for k in ("g1", "rth", "rdot"))
    a_, b_, c_ = CRIT3["P4"]
    M = rth * (M0 + a_ * g1 + b_ * g1 ** 2 + c_ * rdot)
    co = fit(g1, rth, rdot, M, "P4")
    eM, ey = loo(g1, rth, rdot, M, "P4")
    d_co, d_loo = float(np.max(np.abs(co - np.array(CRIT3["P4"])))), float(np.max(np.abs(eM)))
    k_ = [p["name"] for p in pts].index("W3")                   # W3 displaced by +1e-4 in M
    Md = M.copy()
    Md[k_] += 1e-4
    eMd, eyd = loo(g1, rth, rdot, Md, "P4")
    eh = loo_hat(g1, rth, rdot, Md, "P4")
    d_hat = float(np.max(np.abs(eyd - eh)))
    want = 1e-4 / rth[k_]                                       # the others are exact, so the fit without W3 is the planted one and its leave-one-out error is the displacement
    check("PFIT", "planted data that follow the frozen refit exactly (M = r_th (M_0 + a g + b g^2 + c rdot) on the fifteen real (gamma_1, r_th, rdot_th)) return (a, b, c) and zero errors through fit and loo; with W3 displaced "
          "by 1e-4 in M its leave-one-out error is the displacement (in y) and the explicit and the hat-matrix leave-one-out agree",
          d_co <= 1e-10 and d_loo <= 1e-12 and abs(eyd[k_] - want) <= 1e-12 and d_hat <= 1e-12,
          f"max |co - planted| {d_co:.1e}, max |loo| {d_loo:.1e}; displaced point: loo y error {eyd[k_]:.6e}, displacement {want:.6e}; explicit vs hat-matrix {d_hat:.1e}")
    # PPRED4: planted data following the frozen refit pass it and fire it when a or c of the data is changed
    rng = np.random.default_rng(295)
    pg = np.array([0.012, 0.048, 0.075, 0.012, 0.048, 0.075])
    prth = 1.0 + rng.uniform(-6e-3, 6e-3, 6)
    prd = rng.uniform(-6e-4, 1.6e-3, 6)
    prd[np.abs(prd) < 4e-4] = 8e-4
    Mp = np.array([predict("P4", pg[i], prth[i], prd[i], H_FIN, mutated=False) for i in range(6)])
    own = max(abs(Mp[i] - predict("P4", pg[i], prth[i], prd[i], H_FIN)) for i in range(6))
    Ma = np.array([prth[i] * (M0 + 1.05 * a_ * pg[i] + b_ * pg[i] ** 2 + c_ * prd[i]) for i in range(6)])
    Mc = np.array([prth[i] * (M0 + a_ * pg[i] + b_ * pg[i] ** 2 + 0.8 * c_ * prd[i]) for i in range(6)])
    fa = max(abs(Ma[i] - predict("P4", pg[i], prth[i], prd[i], H_FIN)) for i in range(6))
    fc = max(abs(Mc[i] - predict("P4", pg[i], prth[i], prd[i], H_FIN)) for i in range(6))
    check("PPRED4", "planted data that follow the frozen refit exactly pass the P4 test (1e-12), and planted data with a x 1.05, or c x 0.8, fire it (bar 6e-5)",
          own <= 1e-12 and fa > CRIT3["bar_pred"] and fc > CRIT3["bar_pred"], f"max |M - predicted|: own {own:.1e}; a x 1.05 {fa:.1e}; c x 0.8 {fc:.1e} (bar {CRIT3['bar_pred']:g})")
    # PSEP3: three planted route curves, planted new points at known separations
    x = np.array([0.0, 0.0119, 0.0480, 0.0750])
    cur = {"B": dict(g1=x, rth=1.0 - 0.4 * x + 6.0 * x ** 2, rdot=-2e-4 + 5e-3 * x + 0.4 * x ** 2)}
    cur["W"] = dict(g1=x, rth=cur["B"]["rth"] - 1.0e-2 * x / 0.075, rdot=cur["B"]["rdot"] - 7e-4 * x / 0.075)
    cur["C"] = dict(g1=x, rth=cur["B"]["rth"] + 3e-3 * x / 0.075, rdot=cur["B"]["rdot"] - 3e-4 * x / 0.075)
    hd = (1e-6, 1e-5)

    def at(c, g1, dr=0.0, dv=0.0):
        return (float(np.interp(g1, c["g1"], c["rth"])) + dr, float(np.interp(g1, c["g1"], c["rdot"])) + dv)
    allT, allF = dict(B=True, W=True, C=True), dict(B=False, W=False, C=False)
    cases = (("far from all three", 0.0480, lambda k: at(cur["B"], 0.0480, -4e-3, 2e-4), allT),
             ("on route W's curve", 0.0480, lambda k: at(cur["W"], 0.0480), dict(B=True, W=False, C=True)),
             ("2 % beyond the last point, far", 0.0765, lambda k: at(cur["B"], 0.0750, -4e-3, 2e-4), allT),
             ("10 % beyond the last point, far", 0.0825, lambda k: at(cur["B"], 0.0750, -4e-3, 2e-4), allF),
             ("3e-6 in r_th against an h difference of 1e-6", 0.0119, lambda k: at(cur[k], 0.0119, 3e-6), allF),
             ("6e-6 in r_th against an h difference of 1e-6", 0.0119, lambda k: at(cur[k], 0.0119, 6e-6), allT))
    out, badc = [], []
    for nm, g1, vf, want in cases:
        got = {}
        for k in cur:
            v = vf(k)
            got[k] = sep3(g1, v, (v[0] + hd[0], v[1] + hd[1]), cur[k])["informative"]
        out.append(f"{nm}: " + "".join(f"{k}={'Y' if got[k] else 'n'}" for k in cur))
        if got != want:
            badc.append(nm)
    check("PSEP3", "planted points at known separations from three planted route curves (B, W, C) through the separation test: informative or not as planted (far, on a curve, 2 % and 10 % beyond the last point, "
          "3e-6 and 6e-6 against an h difference of 1e-6)", not badc, "; ".join(out) + (f"; WRONG: {badc}" if badc else ""))


# --------------------------------------------------------------------------------------------------------------------
# the analysis
# --------------------------------------------------------------------------------------------------------------------
def derive3(D):
    for e in D.values():
        if mut("no_run") and e.get("name") == "N1" and e["h"] == H_FIN and "err" not in e:
            e["err"] = "member output treated as missing"
    return r.derive(D)


def declared_spec(nm):
    spec = dict((n, s) for n, s, _ in POINTS3)[nm] if nm != REPW3[0] else REPW3[1]
    return (spec[0], spec[1], spec[2] * (1.1 if mut("spec_w") else 1.0))


def analyze_cmd(raw, rn=None):
    t0 = time.time()
    root = os.path.join(raw, "runs")
    D = {}
    for nm, spec, tg in POINTS3:
        for h in H2:
            e = r.load_point(root, spec, h)
            e["name"], e["target"] = nm, tg
            D[(nm, h)] = e
    rep = r.load_point(os.path.join(raw, "runs_rep"), REPW3[1], H_FIN)
    rep["name"], rep["target"] = REPW3[0], None
    D[(REPW3[0], H_FIN)] = rep
    miss = {f"{k[0]} h{k[1]:g}": e["err"] for k, e in D.items() if "err" in e}
    print(f"== Runs read from {root} (and runs_rep): {len(D) - len(miss)} of {len(D)} (point, h) complete" + (f"; missing: {miss}" if miss else ""))
    derive3(D)
    for k, e in D.items():
        if "err" in e and f"{k[0]} h{k[1]:g}" not in miss:
            miss[f"{k[0]} h{k[1]:g}"] = e["err"]
            print(f"    {k[0]} h {k[1]:g}: {e['err']}")
    ok = {k: e for k, e in D.items() if "err" not in e}
    names = [nm for nm, _, _ in POINTS3]
    C = r.CRITR
    bad = {k: (sm.CRIT[k], v) for k, v in g.DECLARED.items() if k != "X" and sm.CRIT[k] != v}
    xg = th.CRIT["X"]["G"][ND]
    spec_bad = [f"{k[0]} h {k[1]:g}" for k, e in D.items() if "p" in e and "i" in e and (list(e["p"]["B"]) != list(declared_spec(k[0])) or list(e["i"]["B"]) != list(declared_spec(k[0])))]
    check("CFG", "settled_mass.CRIT (F1, K, W1, rcut, tol, dxp, dm) equal #41's declared chain, X = 28 for G, the ghost pulse (xc 3, w 0.5), and the (B, xc, w) stored in each run equal the declared one",
          not bad and xg == g.DECLARED["X"] and th.CRIT["pulse"] == r.GHOST and not spec_bad, f"differs: {bad}; X(G) {xg}; pulse {th.CRIT['pulse']}; specs {spec_bad}" if (bad or xg != 28.0 or spec_bad) else "")
    r.set_route(FAM, (0.1, 3.0, 0.5))
    check("FAM0", "route_family with the default normal pulse (B, 3.0, 0.5) returns what threshold.family returns (the pulses, symmetry flag and signs), bit for bit",
          all(r.route_family("G", a_) == r._family_orig("G", a_) for a_ in (0.0, 0.0149, 0.02)), "")
    check("RUNS", "every (point, h) has its pstar, chain and member output and reads without error", not miss, f"missing {sorted(miss)}" if miss else "")
    print("\n== Planted cases")
    check_planted()

    print("\n== The frozen refit (raw README stage 1) against the refit recomputed from the stored fifteen points")
    pts15, c0 = points15()
    co = frozen_check(pts15, c0)
    print(f"    recomputed (a, b, c) = {', '.join(f'{x:+.10f}' for x in co)}; frozen {', '.join(f'{x:+.10f}' for x in CRIT3['P4'])}")

    print("\n== p*_h of each point: bisection to hi - lo <= 1e-14 (threshold.bisect on threshold.point, X = 28), the h of the job as passed")
    acc = Acc()
    for (nm, h), e in ok.items():
        p, pts_ = e["p"], e["pts"]
        lo, hi = b.bracket_after(pts_, 30) if mut("bisect_early") else (p["lo"], p["hi"])
        good = (pts_[0]["fate"] == -1 and pts_[1]["fate"] == 1 and p["h_arg"] == repr(float(h)) and b.bracket_after(pts_, len(pts_)) == (p["lo"], p["hi"])
                and ((hi - lo <= C["bar_bracket"] and p["stopped"] is None) or (p["stopped"] == "static through X" and hi - lo <= C["bar_bracket_X"])))
        acc.add(good)
        print(f"    {nm} {tuple(e['B'])} h {h:g}: p*_h {p['ps']:.16e}, hi - lo {hi - lo:.1e}, {p['runs']} runs, stopped {p['stopped']}{'' if good else '  <-- NOT OK'}")
    check("BR", "bracket ends are collapse and expansion, the replayed bisection equals lo and hi, and either hi - lo <= 1e-14 or the bisection stopped by 'static through X' at hi - lo <= 2e-14; h passed as its float repr", acc.ok)

    print("\n== Does the near-threshold run stall? Lifetime T at p*_h -/+ 1e-8 against 1e-6 (a stall adds about 4 per side)")
    acc = Acc()
    for (nm, h), e in ok.items():
        sl, dm, dp = b.stalls(e)
        acc.add(sl)
        i = e["i"]
        print(f"    {nm} h {h:g}: collapse side T {i['pt6']['T']:.3f} (1e-6) -> {i['point']['T']:.3f} (1e-8), +{dm:.2f}; expansion side {i['up6']['T']:.3f} -> {i['up']['T']:.3f}, +{dp:.2f}: {'stall' if sl else 'NOT a stall'}")
    check("STL", "each point stalls: collapse at p*_h - 1e-6 and 1e-8 and expansion at + 1e-6 and 1e-8, and T gains >= 2.0 per side from 1e-6 to 1e-8", acc.ok)
    acc = Acc()
    for (nm, h), e in ok.items():
        ev = e["m"]["ev40_b"].get("expand") if mut("evt_side") else e["m"]["ev40_a"].get("collapse")
        i_ = e["I"]
        acc.add(ev is not None and list(i_["event"]) == list(ev) and i_["point"]["fate"] == -1 and i_["up"]["fate"] == 1)
    check("EVT", "the chain's first trapped cell equals the pair's collapse-side one (same amplitude p*_h - 1e-8), fate collapse at p*_h - 1e-8 and expansion at p*_h + 1e-8", acc.ok)

    print("\n== The chain (#41's, unchanged) at delta = 1e-8: limit, first law, settled value")
    Fmax = sm.CRIT["F1"] * 2.0 ** sm.CRIT["K"]
    acc = Acc()
    for (nm, h), e in ok.items():
        s, i = e["S"], e["I"]
        fine = np.allclose(np.diff(s["xp"]), h, atol=1e-9)
        acc.add(i["status"] == "limit" and 0.5 * Fmax <= s["F"][-1] < Fmax and s["edge"].min() >= i["Nc"] / 8 and fine)
        print(f"    {nm} h {h:g}: {i['status']}, F_end {s['F'][-1]:.4g}, edge {s['edge'].min():.0f}/{i['Nc']}")
    check("RC5", f"every run ends at the declared limit (F = e^-f on the horizon, {0.5 * Fmax:g} <= F_end < {Fmax:g}), horizon at least Nc/8 cells from the box edges, no gaps", acc.ok)
    acc, worst = Acc(), 0.0
    for (nm, h), e in ok.items():
        s = e["S"]
        i1 = int(np.argmax(s["F"] >= 20.0))
        dMs = s["M"][i1] - s["M"][0]
        with np.errstate(all="ignore"):
            rel = abs(dMs - g.flux_integral_g(s, 0, i1)) / abs(dMs)
        rel = rel if np.isfinite(rel) else float("inf")
        worst = max(worst, rel)
        acc.add(rel <= C["bar_flux"])
    check("RC4", "first law on the horizon, both fields: M(x+_b) - M(x+_a) = integral of (r^2/2) e^f r_- (phi0_+^2 - phi1_+^2), from the first trapped cell to F = 20", acc.ok, f"bar {C['bar_flux']:g}; worst {worst:.1e}")
    acc = Acc()
    for (nm, h), e in ok.items():
        q = e["SS"]
        acc.add(abs(q["dplat"]) <= C["bar_plat"] and q["flux05"] <= C["bar_flux05"])
    check("RC8", "M_AH has settled at the limit: plateau mean over the last 1.0 differs from the previous 1.0 by <= 1e-4 and the mean horizon flux (both fields) over the last 0.5 is <= 1e-5", acc.ok,
          "plateau change " + ", ".join(f"{e['SS']['dplat']:+.1e}" for e in ok.values()))
    acc, det = Acc(), []
    for j, nm in enumerate(names):
        if (nm, H2[0]) not in ok or (nm, H2[1]) not in ok:
            acc.add(False)
            continue
        nm2 = names[(j + 1) % len(names)] if mut("rh_swap") else nm
        d = abs(ok[(nm, H2[0])]["SS"]["plat1"] - ok.get((nm2, H2[1]), ok[(nm, H2[1])])["SS"]["plat1"])
        acc.add(d <= C["bar_hpair"])
        det.append(f"{nm}: {d:.1e}")
    check("RH", "the settled value at two h: |plat(0.01) - plat(0.005)| (plateau mean over the last 1.0 of x+) <= 5e-5", acc.ok, "; ".join(det))

    print("\n== The member at departure (#225's slice_masses and fit_asym on the pair's mean profile; t_dep: the unstable mode's part reaches eps = 1e-3 of r_th) and the throat's speed there (#267's)")
    print("    point  h       p*_h            t_dep    r_th      gamma_1   rdot_th     rms(fit)  M(plateau)")
    for (nm, h), e in ok.items():
        mem = e["mem"]
        print(f"    {nm:5s} {h:<7g} {e['p']['ps']:.10e}  {mem['t_dep']:7.3f}  {mem['r_th']:.6f}  {mem['g1']:.6f}  {e['sp']:+.3e}  {e['rms']:.1e}  {e['SS']['plat1']:.6f}")
    acc, det = Acc(), []
    for nm, spec, tg in POINTS3:
        e = ok.get((nm, H_FIN))
        if e is None:
            acc.add(False)
            continue
        tgv = tg * (1.1 if mut("tgt_scale") else 1.0)
        rel = e["mem"]["g1"] / tgv - 1.0
        acc.add(abs(rel) <= CRIT3["tgt_rel"])
        det.append(f"{nm}: {rel:+.1e}")
    check("TGT", f"gamma_1 at departure (h = {H_FIN:g}) is within {CRIT3['tgt_rel']:.0%} of its target (0.012, 0.048, 0.075)", acc.ok, "; ".join(det))
    acc, det = Acc(), []
    for nm in names:
        a_, c_ = ok.get((nm, H2[0])), ok.get((nm, H2[1]))
        if a_ is None or c_ is None:
            acc.add(False)
            continue
        dg = abs(a_["mem"]["g1"] / c_["mem"]["g1"] - 1.0)
        dr = abs(a_["mem"]["r_th"] / c_["mem"]["r_th"] - 1.0)
        dv = abs(a_["sp"] - c_["sp"])
        acc.add(dg <= C["bar_hst_rel"] and dr <= C["bar_hst_rel"] and dv <= C["bar_hst_rdot"])
        det.append(f"{nm}: gamma_1 {dg:.1e}, r_th {dr:.1e}, rdot_th {dv:.1e}")
    for ln in det:
        print("    HST " + ln)
    check("HST", "the member at departure at the two h (0.01, 0.005) agrees to 1e-4 relative in gamma_1 and r_th, and to 2e-5 absolute in rdot_th (the issue's acceptance)", acc.ok)
    acc = Acc()
    for (nm, h), e in ok.items():
        acc.add(np.isfinite(e["rms"]) and e["rms"] * (10.0 if mut("fitq_scale") else 1.0) <= C["bar_fitq"])
    check("FITQ", "the member's fit of M(r) on the two sides has an rms residual of at most 1e-4 at the slices around departure (a static member describes the state there)", acc.ok)

    print("\n== Separation from the earlier routes at the same gamma_1 (h = 0.005; B = 0 anchors each curve; piecewise-linear and PCHIP, the smaller difference; gamma_1 up to 3 % beyond a curve's last point is extrapolated)")
    refs = {"B": r.fam_b(H_FIN), "W": route_curve(H_FIN, ("W1", "W2", "W3")), "C": route_curve(H_FIN, ("C1", "C2", "C3"))}
    refs["B"] = dict(g1=refs["B"]["g1"], rth=refs["B"]["rth"], rdot=refs["B"]["rdot"])
    for k_, cv in refs.items():
        print(f"    route {k_}: gamma_1 " + ", ".join(f"{x:.5f}" for x in cv["g1"]) + "; r_th " + ", ".join(f"{x:.5f}" for x in cv["rth"]) + "; rdot_th " + ", ".join(f"{x:+.2e}" for x in cv["rdot"]))
    sepacc, SEP = {k_: Acc() for k_ in refs}, {}
    for nm in names:
        a_, c_ = ok.get((nm, H2[0])), ok.get((nm, H_FIN))
        for k_, cv in refs.items():
            if a_ is None or c_ is None:
                sepacc[k_].add(False)
                continue
            cvu = cv
            if mut("sep_self"):
                gg = c_["mem"]["g1"]
                cvu = dict(g1=np.array([0.0, gg - 1e-4, gg, gg + 1e-4]), rth=np.array([c_["mem"]["r_th"]] * 4), rdot=np.array([c_["sp"]] * 4))
            s_ = sep3(c_["mem"]["g1"], (c_["mem"]["r_th"], c_["sp"]), (a_["mem"]["r_th"], a_["sp"]), cvu)
            SEP[(nm, k_)] = s_
            sepacc[k_].add(s_["informative"])
            q = s_["q"]
            print(f"    {nm} against route {k_}: gamma_1 {c_['mem']['g1']:.5f} (beyond its last point by {s_['ext']:.1%}); r_th - route: {q['r_th']['d_lin']:+.2e} (linear), {q['r_th']['d_pchip']:+.2e} (PCHIP), its h difference {q['r_th']['hdiff']:.1e}"
                  f" -> {'ok' if q['r_th']['ok'] else 'no'}; rdot_th - route: {q['rdot']['d_lin']:+.2e}, {q['rdot']['d_pchip']:+.2e}, h difference {q['rdot']['hdiff']:.1e} -> {'ok' if q['rdot']['ok'] else 'no'}: "
                  f"{'informative' if s_['informative'] else 'NOT informative'}")
    for k_, nmk in (("W", "SEPW"), ("C", "SEPC"), ("B", "SEPB")):
        check(nmk, f"the new members' r_th or rdot_th at departure differs from route {k_}'s at the same gamma_1 by more than 5 times its own h difference, all three points" + (" (the issue's)" if k_ != "B" else " (family B, #283's reference; stricter than the issue)"), sepacc[k_].ok)

    print("\n== The frozen refit on the new members: M - predicted, h = 0.005 (gating)")
    print("    P4: M = r_th (M_0 + a g + b g^2 + c rdot_th), M_0 = 0.223086, (a, b, c) = (0.2296318924, 0.2834170383, 0.7657657062), frozen before any run of the third route")
    print("    point   M(plateau)   P4 predicted  M - P4      (leave-one-out scale 7.9e-5) | comparison, gates nothing: M - P3 (#267's)   M - P2' (refit without the speed)")
    res = {}
    for nm in names:
        e = ok.get((nm, H_FIN))
        if e is None:
            continue
        M, g1, rth, rd = e["SS"]["plat1"], e["mem"]["g1"], e["mem"]["r_th"], e["sp"]
        res[nm] = M - predict("P4", g1, rth, rd)
        print(f"    {nm:6s}  {M:.6f}   {predict('P4', g1, rth, rd):.6f}      {res[nm]:+.2e}                                    | {M - predict('P3', g1, rth, rd):+.2e}                {M - predict('P2p', g1, rth, rd):+.2e}")
    if len(res) == len(names) and names:
        worst = max(abs(v) for v in res.values())
        check("P4", "the frozen refit predicts the settled mass at every new member within 6e-5 (h = 0.005; coefficients frozen in the raw README before any run of the third route)", worst <= CRIT3["bar_pred"],
              f"{len(res)} points; max |M - predicted| {worst:.1e} (bar {CRIT3['bar_pred']:g}); " + ", ".join(f"{nm} {res[nm]:+.1e}" for nm in names))
    else:
        check("P4", "the frozen refit predicts the settled mass at every new member within 6e-5", False, "a point is missing")
    resh = {}
    for nm in names:
        e = ok.get((nm, H2[0]))
        if e is not None:
            resh[nm] = e["SS"]["plat1"] - predict("P4", e["mem"]["g1"], e["mem"]["r_th"], e["sp"], H2[0])
    if len(resh) == len(names) and names:
        check("P4-ph", "POST HOC, gates nothing: the frozen (a, b, c) at h = 0.01 with M_0(0.01) = 0.223095 (the coefficients were fitted at h = 0.005)", max(abs(v) for v in resh.values()) <= CRIT3["bar_pred"],
              ", ".join(f"{nm} {resh[nm]:+.1e}" for nm in names), ph=True)

    print("\n== Reproduction of #283's W3 at h = 0.005 from its stored bracket, bit for bit (stored runs of #283 read only)")
    e = ok.get((REPW3[0], H_FIN))
    if e is None:
        check("REP3", "the W3 point of #283 (B, xc, w) = (0.1909, 3.0, 1.0) at h = 0.005 reproduces #283's stored run bit for bit", False, "run missing")
    else:
        ref = r.load_point(os.path.join(I283, "runs"), (0.1543, 3.0, 1.0) if mut("rep_off") else REPW3[1], H_FIN)
        er = ref if "err" not in ref else None
        if er is None:
            check("REP3", "the W3 point of #283 reproduces #283's stored run bit for bit", False, ref["err"])
        else:
            rm = r.read_point(er)
            same_p = e["p"]["ps"] == er["p"]["ps"]
            same_s = all(np.array_equal(e["s"][k], er["s"][k]) for k in er["s"]) and set(e["s"]) == set(er["s"])
            same_l = all(np.array_equal(e["loc"][k], er["loc"][k]) for k in er["loc"]) and set(e["loc"]) == set(er["loc"])
            same_m = {k: v for k, v in e["m"].items() if k != "secs"} == {k: v for k, v in er["m"].items() if k != "secs"}      # all but the wall time
            dM = e["SS"]["plat1"] - rm["M"]
            dd = (e["mem"]["g1"] - rm["mem"]["g1"], e["mem"]["r_th"] - rm["mem"]["r_th"], e["sp"] - rm["sp"], e["mem"]["t_dep"] - rm["mem"]["t_dep"])
            check("REP3", "the W3 point of #283 at h = 0.005 reproduces #283's stored run bit for bit: p*_h, the chain's series, the locus, the member rows (all but the wall time), the plateau, and (gamma_1, r_th, rdot_th, t_dep)",
                  same_p and same_s and same_l and same_m and dM == 0.0 and all(x == 0.0 for x in dd),
                  f"p*_h equal {same_p}, series {same_s}, locus {same_l}, member rows {same_m}; dM {dM:+.1e}, d(gamma_1, r_th, rdot_th, t_dep) {dd[0]:+.1e}, {dd[1]:+.1e}, {dd[2]:+.1e}, {dd[3]:+.1e}")

    print("\n== Post hoc (gates nothing): leverage, nearest fitted point, and a refit with the new points")
    g1v, rthv, rdv, Mv = (np.array([p[k] for p in pts15]) for k in ("g1", "rth", "rdot", "M"))
    A15 = design(g1v, rdv, "P4")
    G15 = np.linalg.inv(A15.T @ A15)
    lev15 = np.einsum("ij,jk,ik->i", A15, G15, A15)
    print(f"    leverage of the fifteen fitted points: max {lev15.max():.3f} ({pts15[int(np.argmax(lev15))]['name']}), mean {lev15.mean():.3f} (3 parameters / 15 points = {3 / 15:.3f})")
    for nm in names:
        e = ok.get((nm, H_FIN))
        if e is None:
            continue
        a1 = design([e["mem"]["g1"]], [e["sp"]], "P4")[0]
        lv = float(a1 @ G15 @ a1)
        dist = [math.sqrt(((e["mem"]["g1"] - p["g1"]) / 1e-3) ** 2 + ((e["mem"]["r_th"] - p["rth"]) / 1e-3) ** 2 + ((e["sp"] - p["rdot"]) / 1e-4) ** 2) for p in pts15]
        jn = int(np.argmin(dist))
        print(f"    {nm}: leverage {lv:.3f} ({'beyond' if lv > lev15.max() else 'within'} the fitted points' maximum); nearest fitted point {pts15[jn]['name']} at {dist[jn]:.2f} in units of (1e-3 in gamma_1, 1e-3 in r_th, 1e-4 in rdot_th)")
    try:
        allp = list(pts15) + [dict(name=nm, route="N", g1=ok[(nm, H_FIN)]["mem"]["g1"], rth=ok[(nm, H_FIN)]["mem"]["r_th"], rdot=ok[(nm, H_FIN)]["sp"], M=ok[(nm, H_FIN)]["SS"]["plat1"]) for nm in names if (nm, H_FIN) in ok]
        g_, r_, d_, M_ = (np.array([p[k] for p in allp]) for k in ("g1", "rth", "rdot", "M"))
        for form in ("P4", "P2"):
            co2 = fit(g_, r_, d_, M_, form)
            eM, _ = loo(g_, r_, d_, M_, form)
            ins = M_ - r_ * (M0 + design(g_, d_, form) @ co2)
            print(f"    refit of the form {form} over these {len(allp)} points: coefficients {', '.join(f'{x:+.4f}' for x in co2)}; in-sample max {np.max(np.abs(ins)):.1e}; leave-one-out max {np.max(np.abs(eM)):.1e}; "
                  f"new points' leave-one-out: " + ", ".join(f"{p['name']} {x:+.1e}" for p, x in zip(allp, eM) if p["route"] == "N"))
        print("    *post hoc: fitted on points that include the ones tested; descriptions, not tests*")
    except Exception as ex:
        print(f"    post hoc refit not computable: {ex!r}")

    if rn:
        print(f"\n== Post hoc (gates nothing): a second roundoff realization of the new points at h = {H_FIN:g} (p*_h bisected from another bracket), {rn}")
        R2 = {}
        for nm, spec, tg in POINTS3:
            e2 = r.load_point(os.path.join(rn, "runs"), spec, H_FIN)
            e2["name"] = nm
            R2[(nm, H_FIN)] = e2
        derive3(R2)
        for nm, spec, tg in POINTS3:
            e2, e = R2[(nm, H_FIN)], ok.get((nm, H_FIN))
            if "err" in e2 or e is None:
                print(f"    {nm}: not available ({e2.get('err')})")
                continue
            M2, g2, r2, d2 = e2["SS"]["plat1"], e2["mem"]["g1"], e2["mem"]["r_th"], e2["sp"]
            print(f"    {nm:5s} p*_h difference {e2['p']['ps'] - e['p']['ps']:+.1e}  M(rn) - M {M2 - e['SS']['plat1']:+.1e}  gamma_1, r_th, rdot_th differences {g2 - e['mem']['g1']:+.1e}, {r2 - e['mem']['r_th']:+.1e}, {d2 - e['sp']:+.1e}  "
                  f"M - P4 (rn) {M2 - predict('P4', g2, r2, d2, H_FIN):+.2e}")
    return finish(t0)


MUTATIONS = {      # name: the checks it newly fails (those failing without a mutation, --baseline, excluded); computed on the development set (raw I295-mac dev/mut_a, baseline dev/ana_dev2.log: no check fails) before being declared in its README
    "asym_sign": ['FRZ', 'P4', 'P4-ph', 'PMEM', 'PTD', 'SEPC', 'SEPW', 'TGT'],
    "bisect_early": ['BR'],
    "coef_a": ['P4', 'P4-ph', 'PPRED4'],
    "coef_c": ['P4', 'P4-ph', 'PPRED4'],
    "crit_dm": ['CFG'],
    "dep_eps": ['FITQ', 'FRZ', 'P4', 'P4-ph', 'PDEP', 'PTD'],
    "drop_point": ['FRZ'],
    "eps_h": ['FITQ', 'FRZ', 'HST', 'P4', 'PTD', 'SEPB', 'SEPC', 'SEPW'],
    "evt_side": ['EVT'],
    "fitq_scale": ['FITQ'],
    "flux_sign": ['RC4'],
    "frz_a": ['FRZ'],
    "frz_m0": ['FRZ', 'PPRED4'],
    "loo_off": ['FRZ'],
    "mass_exponent": ['P4', 'P4-ph', 'RC4', 'REP3'],
    "no_rscale": ['P4', 'P4-ph', 'PPRED4'],
    "no_run": ['HST', 'P4', 'RH', 'RUNS', 'SEPB', 'SEPC', 'SEPW', 'TGT'],
    "no_settle": ['P4', 'P4-ph', 'RC5', 'RC8', 'REP3'],
    "normal_drop": [],
    "rdot_nolapse": ['FRZ', 'PTD'],
    "rep_off": ['REP3'],
    "rh_swap": ['RH'],
    "sep_k": ['PSEP3'],
    "sep_self": ['SEPB', 'SEPC', 'SEPW'],
    "short_limit": ['P4', 'P4-ph', 'RC5', 'RC8', 'REP3'],
    "spec_w": ['CFG'],
    "stl_flat": ['STL'],
    "tgt_scale": ['TGT'],
    "w_scale": ['FAM0', 'HST', 'P4', 'REP3', 'RH', 'RUNS', 'SEPB', 'SEPC', 'SEPW', 'TGT'],
    "xc_offset": ['FAM0', 'HST', 'P4', 'REP3', 'RH', 'RUNS', 'SEPB', 'SEPC', 'SEPW', 'TGT'],
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
        r.MUT = b.MUT = m.MUT = g.MUT = sm.MUT = MUT
        dn.MUT = MUT if MUT in ("theta_coef", "ghost_sign", "scalar_coef") else None
        tl.LMUT = MUT if MUT in ("asym_sign",) else None
        if MUT == "crit_dm":
            sm.CRIT["dm"] = 2.0
    import platform
    import numba as nb
    print(f"settled_mass_route3.py {cmd} {opts}: numba {nb.__version__}, numpy {np.__version__}, python {platform.python_version()}, {platform.node()}, "
          f"mutation {MUT}, {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    if cmd == "refit":
        refit_cmd()
    elif cmd == "tag":
        print(r.tag_of(FAM, r.parse_route(opts["route"])))
    elif cmd == "bracket":
        lo, hi, rows = r.bracket_scan(r.parse_route(opts["route"]), float(opts["h"]))
        print(json.dumps(dict(lo=lo, hi=hi, rows=rows)))
    elif cmd == "pstar":
        r.pstar_cmd(r.parse_route(opts["route"]), float(opts["h"]), opts["out"], tuple(float(x) for x in opts["bracket"].split(",")))
    elif cmd == "job":
        r.job_cmd(r.parse_route(opts["route"]), float(opts["h"]), opts["out"], opts["pstar"])
    elif cmd == "member":
        r.member_cmd(r.parse_route(opts["route"]), float(opts["h"]), opts["out"], opts["pstar"])
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

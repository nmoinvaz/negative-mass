#!/usr/bin/env python3
"""The frozen refit of family G's settled-mass relation tested blind on a route outside its fitted data, chosen for leverage (#305). D = 4, delta = 1e-8.

#295 refitted y = M / r_th - M_0 = a gamma_1 + b gamma_1^2 + c rdot_th over fifteen members at h = 0.005 and froze (a, b, c) = (0.2296318924, 0.2834170383, 0.7657657062). Its third route lay inside the
fitted data (leverage 0.07 to 0.33 against a fitted maximum of 0.393). This driver tests the same frozen coefficients on a route whose normal pulse (x_c, w) was chosen in development at h = 0.02 for
its leverage in the fit's design (gamma_1, gamma_1^2, rdot_th), with rdot_th outside [-6.7e-4, +1.6e-3], by the rule of the raw README (stage 1, scripts/devselect.py), and runs the chain as #283 and #295.

A thin driver. It imports settled_mass_route3.py (#295), which imports settled_mass_route.py (#283), settled_mass_b.py (#266), settled_mass_member.py (#267) and everything below them, all unchanged.

    python3 settled_mass_route4.py tag|bracket|pstar|job|member ...                  as settled_mass_route.py (#283), for the route's points
    python3 settled_mass_route4.py planted [--mutate=NAME]                           planted cases
    python3 settled_mass_route4.py analyze --raw=DIR [--mutate=NAME] [--baseline=FILE] [--rn=DIR]

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
import settled_mass_route3 as s3           # noqa: E402  (#295's driver: the fifteen points, the frozen refit, the separation test; imports #283's, #266's, #267's unchanged)
r = s3.r
b, m, g, th, tl, sm, dn, td = s3.b, s3.m, s3.g, s3.th, s3.tl, s3.sm, s3.dn, s3.td

FAM = "G"
H2 = r.H2                                       # (0.01, 0.005)
H_FIN = r.H_FIN
MUT = None
ND = 2
I295 = os.path.join(s3.RAW_ROOT, "I295-mac")
I283 = s3.I283
M0 = s3.M0
CRIT3 = s3.CRIT3

POINTS4 = (("X1", (0.09112, 4.0, 1.5), 0.012), ("X2", (0.1695, 4.0, 1.5), 0.048), ("X3", (0.2044, 4.0, 1.5), 0.075))      # (name, (B, xc, w), target gamma_1): the route selected by the rule of the raw README (stage 1b)
REP4 = ("RN3", (0.136, 3.5, 0.75))              # #295's N3, reproduced at h = 0.005 from its stored bracket
TARGETS = (0.012, 0.048, 0.075)
CRIT4 = dict(
    rdot_range=(-6.7e-4, 1.6e-3),               # the issue's range of the fitted ṙ_th
    lev_issue=0.393,                            # the issue's rounding of the fitted maximum leverage (0.3932, W3)
    lev_bar_scale=10.0,                         # mutation lev_bar (the bar x 10: above the largest leverage of the route, 3.06; a x 1.6 first tried did not fire)
)

results = s3.results
check = s3.check
Acc = s3.Acc
BASE = s3.BASE


def mut(name):
    return MUT == name


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


# --------------------------------------------------------------------------------------------------------------------
# leverage in the frozen design
# --------------------------------------------------------------------------------------------------------------------
def fitted_design(pts):
    """Rows (gamma_1, gamma_1^2, rdot_th) of the fitted points, (A^T A)^-1, and the leverage of each fitted point from the hat matrix A pinv(A) (an independent route to the explicit inverse)."""
    g1, rdot = np.array([p["g1"] for p in pts]), np.array([p["rdot"] for p in pts])
    A = s3.design(g1, rdot, "P4")
    G = np.linalg.inv(A.T @ A)
    hat = np.diag(A @ np.linalg.pinv(A))
    return A, G, hat


def leverage(g1, rdot, A, G):
    """L = a^T (A^T A)^-1 a for a = (gamma_1, gamma_1^2, rdot_th). Mutation lev_design: the design without the gamma_1^2 column."""
    a = s3.design([g1], [rdot], "P4")[0]
    if mut("lev_design"):
        cols = [0, 2]
        G2 = np.linalg.inv(A[:, cols].T @ A[:, cols])
        return float(a[cols] @ G2 @ a[cols])
    return float(a @ G @ a)


def leverage_svd(a, A):
    """The same quantity from the singular value decomposition of A: sum_i (v_i . a)^2 / sigma_i^2 (an independent route; (A^T A)^-1 = V S^-2 V^T)."""
    _, S, Vt = np.linalg.svd(A, full_matrices=False)
    return float(np.sum((Vt @ a) ** 2 / S ** 2))


def check_plev(pts):
    """PLEV: leverage planted. (1) the leverage function at each fitted point's (gamma_1, rdot_th) equals the hat-matrix diagonal; (2) it equals the SVD form at random rows of the fitted size; (3) it is
    quadratic in the row (L(s a) = s^2 L(a)); (4) the decision 'above the fitted maximum' is right for planted rows at 0.9 and 1.1 times W3's row (leverage 0.81 and 1.21 times W3's)."""
    A, G, hat = fitted_design(pts)
    e1 = max(abs(leverage(p["g1"], p["rdot"], A, G) - h_) for p, h_ in zip(pts, hat))
    rng = np.random.default_rng(305)
    rows = np.stack([rng.uniform(0.0, 0.08, 6), rng.uniform(0.0, 6e-3, 6), rng.uniform(-8e-4, 2e-3, 6)], 1)
    e2 = max(abs(float(a @ G @ a) - leverage_svd(a, A)) / leverage_svd(a, A) for a in rows)
    k3 = [p["name"] for p in pts].index("W3")
    aw = A[k3]
    e3 = max(abs(float((s * aw) @ G @ (s * aw)) - s * s * hat[k3]) for s in (0.9, 1.1, 2.0))
    Lmax = float(hat.max())
    dec = [float((s * aw) @ G @ (s * aw)) > Lmax * (1.0 + 1e-9) for s in (0.9, 1.0, 1.1)]
    ok = e1 <= 1e-9 and e2 <= 1e-9 and e3 <= 1e-12 and dec == [False, False, True]
    check("PLEV", "planted leverage: the function at the fitted points equals the hat-matrix diagonal (1e-9), at random rows the SVD form (1e-9 relative), L(s a) = s^2 L(a) (1e-12), and 'above the fitted maximum' "
          "is false at 0.9 and 1.0 and true at 1.1 times W3's row", ok, f"max |L - hat| {e1:.1e}; |L - svd|/L {e2:.1e}; scaling {e3:.1e}; decisions at 0.9, 1.0, 1.1 {dec}")


# --------------------------------------------------------------------------------------------------------------------
# curves of the earlier routes (separation test)
# --------------------------------------------------------------------------------------------------------------------
def route_curve_N(h):
    """#295's route N at h (its stored runs, B = 0 as the anchor at gamma_1 = 0): gamma_1, r_th, rdot_th, ascending in gamma_1."""
    anchor = r.fam_b(h)
    g1, rth, rdot = [float(anchor["g1"][0])], [float(anchor["rth"][0])], [float(anchor["rdot"][0])]
    for nm, spec, _ in s3.POINTS3:
        x = r.read_point(r.load_point(os.path.join(I295, "runs"), spec, h))
        g1.append(x["mem"]["g1"]), rth.append(x["mem"]["r_th"]), rdot.append(x["sp"])
    o = np.argsort(g1)
    return dict(g1=np.array(g1)[o], rth=np.array(rth)[o], rdot=np.array(rdot)[o])


def points_N(h=H_FIN):
    """#295's three stored members at h: name, gamma_1, r_th, rdot_th, M."""
    out = []
    for nm, spec, _ in s3.POINTS3:
        e = r.load_point(os.path.join(I295, "runs"), spec, h)
        x = r.read_point(e)
        out.append(dict(name=nm, route="N", g1=x["mem"]["g1"], rth=x["mem"]["r_th"], rdot=x["sp"], M=x["M"]))
    return out


# --------------------------------------------------------------------------------------------------------------------
# where a miss belongs (raw README stage 2a): the verdict map
# --------------------------------------------------------------------------------------------------------------------
def rms_at(d, t_dep):
    """FITQ's quantity: the largest rms of the static-member fit among the stored rows within 0.5 of the departure time."""
    v = [x["rms_m"] for x in d["m"]["rows"] if abs(x["t"] - t_dep) <= 0.5]
    return max(v) if v else float("nan")


def belong_verdict(f, R):
    """f = g_X / g_N (the ratio of the geometric means of |M - P4| / rms over the new and over #295's held-out members), R = per new point the three ratios R_q = A_q / |delta*_q|.
    Returns (M2 label, per-point labels, verdict) by the map of the raw README stage 2a."""
    m2 = "consistent" if 1.0 / 3.0 <= f <= 3.0 else ("apart" if (f > 10.0 or f < 0.1) else "between")
    lab = ["reading-compatible" if max(x) >= 1.0 else ("form-side" if max(x) < 1.0 / 3.0 else "between") for x in R]
    if m2 == "consistent" and all(x == "reading-compatible" for x in lab):
        v = "reading"
    elif m2 == "apart" and all(x == "form-side" for x in lab):
        v = "form"
    else:
        v = "not separable by these data"
    if mut("bel_swap"):
        v = {"reading": "form", "form": "reading"}.get(v, v)
    return f, m2, lab, v


def check_pbel():
    """PBEL: the verdict map on planted inputs: planted answers for reading, form, and three that must come out 'not separable'."""
    cases = (("consistent and compatible", 1.1, [[2.0, 0.1, 0.1]] * 3, "reading"), ("apart and form-side", 25.0, [[0.1, 0.2, 0.3]] * 3, "form"),
             ("apart but compatible", 25.0, [[2.0, 0.1, 0.1]] * 3, "not separable by these data"), ("consistent but form-side", 1.1, [[0.1, 0.1, 0.1]] * 3, "not separable by these data"),
             ("between on M2", 5.0, [[0.1, 0.1, 0.1]] * 3, "not separable by these data"), ("one point between", 25.0, [[0.1, 0.1, 0.1], [0.1, 0.1, 0.1], [0.5, 0.1, 0.1]], "not separable by these data"))
    bad = [nm for nm, f, R, want in cases if belong_verdict(f, R)[3] != want]
    check("PBEL", "planted inputs through the verdict map of the raw README stage 2a: reading, form, and four combinations that must come out 'not separable by these data'", not bad, f"wrong: {bad}" if bad else f"{len(cases)} planted cases")


# --------------------------------------------------------------------------------------------------------------------
# the analysis
# --------------------------------------------------------------------------------------------------------------------
def derive4(D):
    for e in D.values():
        if mut("no_run") and e.get("name") == POINTS4[0][0] and e["h"] == H_FIN and "err" not in e:
            e["err"] = "member output treated as missing"
    return r.derive(D)


def declared_spec(nm):
    spec = dict((n, s) for n, s, _ in POINTS4)[nm] if nm != REP4[0] else REP4[1]
    return (spec[0], spec[1], spec[2] * (1.1 if mut("spec_w") else 1.0))


def analyze_cmd(raw, rn=None):
    t0 = time.time()
    root = os.path.join(raw, "runs")
    D = {}
    for nm, spec, tg in POINTS4:
        for h in H2:
            e = r.load_point(root, spec, h)
            e["name"], e["target"] = nm, tg
            D[(nm, h)] = e
    rep = r.load_point(os.path.join(raw, "runs_rep"), REP4[1], H_FIN)
    rep["name"], rep["target"] = REP4[0], None
    D[(REP4[0], H_FIN)] = rep
    miss = {f"{k[0]} h{k[1]:g}": e["err"] for k, e in D.items() if "err" in e}
    print(f"== Runs read from {root} (and runs_rep): {len(D) - len(miss)} of {len(D)} (point, h) complete" + (f"; missing: {miss}" if miss else ""))
    derive4(D)
    for k, e in D.items():
        if "err" in e and f"{k[0]} h{k[1]:g}" not in miss:
            miss[f"{k[0]} h{k[1]:g}"] = e["err"]
            print(f"    {k[0]} h {k[1]:g}: {e['err']}")
    ok = {k: e for k, e in D.items() if "err" not in e}
    names = [nm for nm, _, _ in POINTS4]
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
    s3.check_planted()
    pts15, c0 = s3.points15()
    check_plev(pts15)
    check_pbel()

    print("\n== The frozen refit (#295's raw README stage 1) against the refit recomputed from the stored fifteen points")
    co = s3.frozen_check(pts15, c0)
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
    for nm, spec, tg in POINTS4:
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
    check("HST", "the member at departure at the two h (0.01, 0.005) agrees to 1e-4 relative in gamma_1 and r_th, and to 2e-5 absolute in rdot_th (the issue's acceptance, as #283)", acc.ok)
    acc = Acc()
    for (nm, h), e in ok.items():
        acc.add(np.isfinite(e["rms"]) and e["rms"] * (10.0 if mut("fitq_scale") else 1.0) <= C["bar_fitq"])
    check("FITQ", "the member's fit of M(r) on the two sides has an rms residual of at most 1e-4 at the slices around departure (a static member describes the state there)", acc.ok)

    print("\n== Leverage in the frozen fit's design and ṙ_th against the fitted range (the issue's first acceptance; h = 0.005)")
    A15, G15, hat15 = fitted_design(pts15)
    Lmax = float(hat15.max())
    kmax = int(np.argmax(hat15))
    lev = {}
    for nm in names:
        e = ok.get((nm, H_FIN))
        if e is None:
            continue
        gg, rd = e["mem"]["g1"], e["sp"] * (0.0 if mut("lev_rdot") else 1.0)
        lev[nm] = leverage(gg, rd, A15, G15)
        print(f"    {nm}: gamma_1 {gg:.6f}, rdot_th {e['sp']:+.4e}: leverage {lev[nm]:.4f} ({'above' if lev[nm] > Lmax else 'within'} the fitted maximum {Lmax:.4f}, {pts15[kmax]['name']}); "
              f"rdot_th {'outside' if (e['sp'] < CRIT4['rdot_range'][0] or e['sp'] > CRIT4['rdot_range'][1]) else 'inside'} [{CRIT4['rdot_range'][0]:+.1e}, {CRIT4['rdot_range'][1]:+.1e}]")
    bar_lev = max(Lmax, CRIT4["lev_issue"]) * (CRIT4["lev_bar_scale"] if mut("lev_bar") else 1.0)
    check("LEV", f"at least one new member (h = 0.005) has leverage above the fitted points' maximum {Lmax:.4f} (W3; the issue's 0.393) in the frozen fit's design (gamma_1, gamma_1^2, rdot_th)",
          len(lev) == len(names) and any(v > bar_lev for v in lev.values()), ", ".join(f"{nm} {v:.4f}" for nm, v in lev.items()) + f"; bar {bar_lev:.4f}")
    lo_, hi_ = CRIT4["rdot_range"] if not mut("outr_range") else (-1e-2, 1e-2)
    check("OUTR", "at least one new member (h = 0.005) has rdot_th outside [-6.7e-4, +1.6e-3] (the issue's range of the fitted points' rdot_th)",
          any(ok[(nm, H_FIN)]["sp"] < lo_ or ok[(nm, H_FIN)]["sp"] > hi_ for nm in names if (nm, H_FIN) in ok), ", ".join(f"{nm} {ok[(nm, H_FIN)]['sp']:+.2e}" for nm in names if (nm, H_FIN) in ok))

    print("\n== Separation from the earlier routes at the same gamma_1 (h = 0.005; B = 0 anchors each curve; piecewise-linear and PCHIP, the smaller difference; gamma_1 up to 3 % beyond a curve's last point is extrapolated)")
    refs = {"B": r.fam_b(H_FIN), "W": s3.route_curve(H_FIN, ("W1", "W2", "W3")), "C": s3.route_curve(H_FIN, ("C1", "C2", "C3")), "N": route_curve_N(H_FIN)}
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
            s_ = s3.sep3(c_["mem"]["g1"], (c_["mem"]["r_th"], c_["sp"]), (a_["mem"]["r_th"], a_["sp"]), cvu)
            SEP[(nm, k_)] = s_
            sepacc[k_].add(s_["informative"])
            q = s_["q"]
            print(f"    {nm} against route {k_}: gamma_1 {c_['mem']['g1']:.5f} (beyond its last point by {s_['ext']:.1%}); r_th - route: {q['r_th']['d_lin']:+.2e} (linear), {q['r_th']['d_pchip']:+.2e} (PCHIP), its h difference {q['r_th']['hdiff']:.1e}"
                  f" -> {'ok' if q['r_th']['ok'] else 'no'}; rdot_th - route: {q['rdot']['d_lin']:+.2e}, {q['rdot']['d_pchip']:+.2e}, h difference {q['rdot']['hdiff']:.1e} -> {'ok' if q['rdot']['ok'] else 'no'}: "
                  f"{'informative' if s_['informative'] else 'NOT informative'}")
    for k_, nmk in (("W", "SEPW"), ("C", "SEPC"), ("B", "SEPB"), ("N", "SEPN")):
        what = {"W": "route W (#283)", "C": "route C (#283)", "B": "family B (#266)", "N": "route N (#295)"}[k_]
        check(nmk, f"the new members' r_th or rdot_th at departure differs from {what}'s at the same gamma_1 by more than 5 times its own h difference, all three points", sepacc[k_].ok)

    print("\n== The frozen refit on the new members: M - predicted, h = 0.005 (gating)")
    print("    P4: M = r_th (M_0 + a g + b g^2 + c rdot_th), M_0 = 0.223086, (a, b, c) = (0.2296318924, 0.2834170383, 0.7657657062), frozen in #295 before any run of its third route")
    print("    point   leverage  M(plateau)   P4 predicted  M - P4      (frozen leave-one-out scale 7.9e-5) | comparison, gates nothing: M - P3 (#267's)   M - P2' (refit without the speed)")
    res = {}
    for nm in names:
        e = ok.get((nm, H_FIN))
        if e is None:
            continue
        M, g1, rth, rd = e["SS"]["plat1"], e["mem"]["g1"], e["mem"]["r_th"], e["sp"]
        res[nm] = M - s3.predict("P4", g1, rth, rd)
        print(f"    {nm:6s}  {lev.get(nm, float('nan')):.4f}   {M:.6f}   {s3.predict('P4', g1, rth, rd):.6f}      {res[nm]:+.2e}                                    | {M - s3.predict('P3', g1, rth, rd):+.2e}                {M - s3.predict('P2p', g1, rth, rd):+.2e}")
    if len(res) == len(names) and names:
        worst = max(abs(v) for v in res.values())
        check("P4", "the frozen refit predicts the settled mass at every new member within 6e-5 (h = 0.005; coefficients frozen in #295's raw README before any run of its third route, and before this issue)", worst <= CRIT3["bar_pred"],
              f"{len(res)} points; max |M - predicted| {worst:.1e} (bar {CRIT3['bar_pred']:g}); " + ", ".join(f"{nm} {res[nm]:+.1e}" for nm in names))
    else:
        check("P4", "the frozen refit predicts the settled mass at every new member within 6e-5", False, "a point is missing")
    resh = {}
    for nm in names:
        e = ok.get((nm, H2[0]))
        if e is not None:
            resh[nm] = e["SS"]["plat1"] - s3.predict("P4", e["mem"]["g1"], e["mem"]["r_th"], e["sp"], H2[0])
    if len(resh) == len(names) and names:
        check("P4-ph", "POST HOC, gates nothing: the frozen (a, b, c) at h = 0.01 with M_0(0.01) = 0.223095 (the coefficients were fitted at h = 0.005)", max(abs(v) for v in resh.values()) <= CRIT3["bar_pred"],
              ", ".join(f"{nm} {resh[nm]:+.1e}" for nm in names), ph=True)

    print("\n== Reproduction of #295's N3 at h = 0.005 from its stored bracket, bit for bit (stored runs of #295 read only)")
    e = ok.get((REP4[0], H_FIN))
    if e is None:
        check("REP4", "the N3 point of #295 (B, xc, w) = (0.136, 3.5, 0.75) at h = 0.005 reproduces #295's stored run bit for bit", False, "run missing")
    else:
        ref = r.load_point(os.path.join(I295, "runs"), (0.1101, 3.5, 0.75) if mut("rep_off") else REP4[1], H_FIN)
        er = ref if "err" not in ref else None
        if er is None:
            check("REP4", "the N3 point of #295 reproduces #295's stored run bit for bit", False, ref["err"])
        else:
            rm = r.read_point(er)
            same_p = e["p"]["ps"] == er["p"]["ps"]
            same_s = all(np.array_equal(e["s"][k], er["s"][k]) for k in er["s"]) and set(e["s"]) == set(er["s"])
            same_l = all(np.array_equal(e["loc"][k], er["loc"][k]) for k in er["loc"]) and set(e["loc"]) == set(er["loc"])
            same_m = {k: v for k, v in e["m"].items() if k != "secs"} == {k: v for k, v in er["m"].items() if k != "secs"}      # all but the wall time
            dM = e["SS"]["plat1"] - rm["M"]
            dd = (e["mem"]["g1"] - rm["mem"]["g1"], e["mem"]["r_th"] - rm["mem"]["r_th"], e["sp"] - rm["sp"], e["mem"]["t_dep"] - rm["mem"]["t_dep"])
            check("REP4", "the N3 point of #295 at h = 0.005 reproduces #295's stored run bit for bit: p*_h, the chain's series, the locus, the member rows (all but the wall time), the plateau, and (gamma_1, r_th, rdot_th, t_dep)",
                  same_p and same_s and same_l and same_m and dM == 0.0 and all(x == 0.0 for x in dd),
                  f"p*_h equal {same_p}, series {same_s}, locus {same_l}, member rows {same_m}; dM {dM:+.1e}, d(gamma_1, r_th, rdot_th, t_dep) {dd[0]:+.1e}, {dd[1]:+.1e}, {dd[2]:+.1e}, {dd[3]:+.1e}")

    print("\n== Post hoc (gates nothing): nearest fitted point, leverage in the eighteen-point design, and refits with the new points")
    ptsN = points_N(H_FIN)
    allN = list(pts15) + ptsN
    A18, G18, hat18 = fitted_design(allN)
    print(f"    leverage of the fifteen fitted points: max {hat15.max():.3f} ({pts15[kmax]['name']}), mean {hat15.mean():.3f}; of the eighteen (with #295's N1, N2, N3): max {hat18.max():.3f} ({allN[int(np.argmax(hat18))]['name']}); "
          f"#295's N1, N2, N3 in the fifteen's design: " + ", ".join(f"{p['name']} {leverage(p['g1'], p['rdot'], A15, G15):.3f}" for p in ptsN))
    for nm in names:
        e = ok.get((nm, H_FIN))
        if e is None:
            continue
        dist = [math.sqrt(((e["mem"]["g1"] - p["g1"]) / 1e-3) ** 2 + ((e["mem"]["r_th"] - p["rth"]) / 1e-3) ** 2 + ((e["sp"] - p["rdot"]) / 1e-4) ** 2) for p in allN]
        jn = int(np.argmin(dist))
        print(f"    {nm}: leverage in the fifteen's design {leverage(e['mem']['g1'], e['sp'], A15, G15):.3f}, in the eighteen's {leverage(e['mem']['g1'], e['sp'], A18, G18):.3f}; nearest of the eighteen {allN[jn]['name']} at {dist[jn]:.2f} "
              f"in units of (1e-3 in gamma_1, 1e-3 in r_th, 1e-4 in rdot_th); gamma_1 {'inside' if e['mem']['g1'] <= max(p['g1'] for p in pts15) else 'beyond'} the fitted range, "
              f"rdot_th {'inside' if min(p['rdot'] for p in pts15) <= e['sp'] <= max(p['rdot'] for p in pts15) else 'outside'} it")
    try:
        newp = [dict(name=nm, route="X", g1=ok[(nm, H_FIN)]["mem"]["g1"], rth=ok[(nm, H_FIN)]["mem"]["r_th"], rdot=ok[(nm, H_FIN)]["sp"], M=ok[(nm, H_FIN)]["SS"]["plat1"]) for nm in names if (nm, H_FIN) in ok]
        for label, base_ in (("the fifteen and the new points", list(pts15)), ("the fifteen, #295's three and the new points", allN)):
            allp = base_ + newp
            g_, r_, d_, M_ = (np.array([p[k] for p in allp]) for k in ("g1", "rth", "rdot", "M"))
            co2 = s3.fit(g_, r_, d_, M_, "P4")
            eM, _ = s3.loo(g_, r_, d_, M_, "P4")
            ins = M_ - r_ * (M0 + s3.design(g_, d_, "P4") @ co2)
            print(f"    refit of P4 over {label} ({len(allp)} points): coefficients {', '.join(f'{x:+.4f}' for x in co2)} (frozen {', '.join(f'{x:+.4f}' for x in CRIT3['P4'])}); in-sample max {np.max(np.abs(ins)):.1e}; "
                  f"leave-one-out max {np.max(np.abs(eM)):.1e}; new points' leave-one-out: " + ", ".join(f"{p['name']} {x:+.1e}" for p, x in zip(allp, eM) if p["route"] == "X"))
        print("    *post hoc: fitted on points that include the ones tested; descriptions, not tests*")
    except Exception as ex:
        print(f"    post hoc refit not computable: {ex!r}")

    print("\n== Post hoc (gates nothing): #283's description of P3's residual on the new members (its constant, K = -3.365e-3 from six points of routes W and C, is not refitted)")
    try:
        K283 = -3.365e-3
        for nm in names:
            e = ok.get((nm, H_FIN))
            if e is None:
                continue
            g1_, rth_, rd_ = e["mem"]["g1"], e["mem"]["r_th"], e["sp"]
            d_ = s3.sep3(g1_, (rth_, rd_), (rth_, rd_), refs["B"], k=0.0)["q"]["r_th"]["d_lin"]
            obs = e["SS"]["plat1"] - s3.predict("P3", g1_, rth_, rd_, H_FIN)
            print(f"    {nm}: r_th - family B's {d_:+.2e}; K x that {K283 * d_:+.2e}; observed M - P3 (#267's) {obs:+.2e}; difference {obs - K283 * d_:+.1e}")
    except Exception as ex:
        print(f"    not computable: {ex!r}")

    print("\n== Post hoc (gates nothing): the frozen-fit residual against the leave-one-out scale at the new members")
    for nm in names:
        e = ok.get((nm, H_FIN))
        if e is None or nm not in lev:
            continue
        pm = s3.predict("P4", e["mem"]["g1"], e["mem"]["r_th"], e["sp"])
        print(f"    {nm}: leverage {lev[nm]:.3f}; M - P4 {res.get(nm, float('nan')):+.2e}; ratio to the largest frozen leave-one-out error 7.9e-5: {abs(res.get(nm, float('nan'))) / 7.9e-5:.2f}; predicted M {pm:.6f}")

    print("\n== Post hoc (gates nothing): the static-member fit's rms at departure (FITQ's quantity) at every stored member at h = 0.005, against |rdot_th|")
    try:
        rows = []
        R266 = m.load266(H_FIN)
        for B_ in (0.0,) + tuple(s3.B_266):
            rows.append((f"B{B_:g}", R266[B_]["mem"]["g1"], R266[B_]["mem"]["r_th"], R266[B_]["speed"], rms_at(R266[B_]["e"], R266[B_]["mem"]["t_dep"])))
        for ee in s3.E_267:
            base = os.path.join(s3.I267, "runs", m.TAG, f"h{H_FIN:g}", f"e{ee:g}")
            d = dict(e=ee, h=H_FIN, m=json.load(open(base + "/member/member.json")), loc=dict(np.load(base + "/member/locus.npz")))
            mem_ = b.member_at(d, CRIT3.get("eps", 1e-3), "m")
            rows.append((f"e{ee}", mem_["g1"], mem_["r_th"], m.speed_at(d["loc"], mem_["t_dep"]), rms_at(d, mem_["t_dep"])))
        for nm_, spec_ in s3.POINTS_283:
            e_ = r.load_point(os.path.join(I283, "runs"), spec_, H_FIN)
            x_ = r.read_point(e_)
            rows.append((nm_, x_["mem"]["g1"], x_["mem"]["r_th"], x_["sp"], rms_at(e_, x_["mem"]["t_dep"])))
        for nm_, spec_, _ in s3.POINTS3:
            e_ = r.load_point(os.path.join(I295, "runs"), spec_, H_FIN)
            x_ = r.read_point(e_)
            rows.append((nm_, x_["mem"]["g1"], x_["mem"]["r_th"], x_["sp"], rms_at(e_, x_["mem"]["t_dep"])))
        for nm_ in names:
            e_ = ok.get((nm_, H_FIN))
            if e_ is not None:
                rows.append((nm_, e_["mem"]["g1"], e_["mem"]["r_th"], e_["sp"], rms_at(e_, e_["mem"]["t_dep"])))
        for nm_, g1_, rth_, sp_, rms_ in rows:
            print(f"    {nm_:6s} gamma_1 {g1_:.5f}  r_th - 1 {rth_ - 1:+.5f}  rdot_th {sp_:+.3e}  rms(fit) {rms_:.2e}{'  (above FITQ bar 1e-4)' if rms_ > 1e-4 else ''}")
        arr = np.array([[abs(x[1]), abs(x[2] - 1), abs(x[3]), x[4]] for x in rows])
        print(f"    correlation of rms with |gamma_1|, |r_th - 1|, |rdot_th| over {len(rows)} members: " + ", ".join(f"{np.corrcoef(arr[:, j], arr[:, 3])[0, 1]:+.3f}" for j in range(3)))
    except Exception as ex:
        print(f"    not computable: {ex!r}")

    print("\n== Post hoc (gates nothing): where the miss belongs (raw README stage 2a: M1, M2, M3 and the verdict map)")
    try:
        arr_v = np.array([abs(x[3]) for x in rows])
        arr_r = np.array([x[4] for x in rows])
        c1 = float(np.corrcoef(arr_v, arr_r)[0, 1])
        print(f"    M1: correlation of the static fit's rms with |rdot_th| over {len(rows)} members {c1:+.3f}: " + ("one coordinate (>= 0.9): a poor member fit and a large throat speed are not separable in these data" if c1 >= 0.9 else "below 0.9"))
        held = []
        for p_ in ptsN:
            held.append(("N", p_["name"], p_["M"] - s3.predict("P4", p_["g1"], p_["rth"], p_["rdot"]), [x[4] for x in rows if x[0] == p_["name"]][0], leverage(p_["g1"], p_["rdot"], A15, G15)))
        for nm_ in names:
            if nm_ in res and (nm_, H_FIN) in ok:
                held.append(("X", nm_, res[nm_], rms_at(ok[(nm_, H_FIN)], ok[(nm_, H_FIN)]["mem"]["t_dep"]), lev.get(nm_, float("nan"))))
        print("    M2: held-out members: name, leverage, static-fit rms, M - P4, |M - P4| / rms")
        for grp, nm_, miss_, rms_, lv_ in held:
            print(f"        {nm_:3s} {lv_:.3f}  {rms_:.2e}  {miss_:+.2e}  {abs(miss_) / rms_:.3f}")
        gm = {k_: math.exp(np.mean([math.log(abs(h_[2]) / h_[3]) for h_ in held if h_[0] == k_])) for k_ in ("N", "X")}
        f_ = gm["X"] / gm["N"]
        print(f"        geometric means g_N {gm['N']:.3f}, g_X {gm['X']:.3f}; g_X / g_N = {f_:.1f}")
        a_, b_, c_ = CRIT3["P4"]

        def m3(e_, miss_):
            """delta*_q, A_q, R_q of one stored member (h = 0.005) with the miss M - P4."""
            mem1, sp1 = e_["mem"], e_["sp"]
            mem2 = b.member_at(e_, 1e-2, "m")
            sp2 = m.speed_at(e_["loc"], mem2["t_dep"])
            g1_, rth_ = mem1["g1"], mem1["r_th"]
            der = (rth_ * (a_ + 2 * b_ * g1_), M0 + a_ * g1_ + b_ * g1_ ** 2 + c_ * sp1, rth_ * c_)
            dstar = [miss_ / x for x in der]
            A_ = (abs(mem2["g1"] - g1_), abs(mem2["r_th"] - rth_), abs(sp2 - sp1))
            return dstar, A_, [A_[i] / abs(dstar[i]) for i in range(3)], mem2["t_dep"]
        Rq, labs = [], []
        print("    M3: per new member: shift that would close the miss, delta*_q = (M - P4) / dP4/dq, the ambiguity A_q = |q(eps = 1e-2) - q(eps = 1e-3)| and R_q = A_q / |delta*_q|, q = gamma_1, r_th, rdot_th")
        for nm_ in names:
            e_ = ok.get((nm_, H_FIN))
            if e_ is None or nm_ not in res:
                continue
            dstar, A_, R_, t2_ = m3(e_, res[nm_])
            Rq.append(R_)
            print(f"        {nm_}: delta*: {dstar[0]:+.2e}, {dstar[1]:+.2e}, {dstar[2]:+.2e}; A: {A_[0]:.2e}, {A_[1]:.2e}, {A_[2]:.2e} (t_dep(1e-2) {t2_:.3f}); R_q: {R_[0]:.2f}, {R_[1]:.2f}, {R_[2]:.2f}")
        print("    M3 calibration (added after the development set showed R_q of order 1 to 7 on the new members; not part of the map): the same on #295's held-out members, whose misses are 1e-6")
        for (nm_, spec_, _), p_ in zip(s3.POINTS3, ptsN):
            dstar, A_, R_, t2_ = m3(r.derive({0: dict(r.load_point(os.path.join(I295, "runs"), spec_, H_FIN), name=nm_, h=H_FIN)})[0], p_["M"] - s3.predict("P4", p_["g1"], p_["rth"], p_["rdot"]))
            print(f"        {nm_}: delta*: {dstar[0]:+.2e}, {dstar[1]:+.2e}, {dstar[2]:+.2e}; A: {A_[0]:.2e}, {A_[1]:.2e}, {A_[2]:.2e}; R_q: {R_[0]:.1f}, {R_[1]:.1f}, {R_[2]:.1f}")
        f2, m2, labs, verdict = belong_verdict(f_, Rq)
        print(f"    M2: {m2}; M3 per point: {', '.join(f'{nm_} {l}' for nm_, l in zip(names, labs))}")
        print(f"    VERDICT by the map of the raw README stage 2a: {verdict}. A reading error larger than M3 measures is not excluded by these data.")
    except Exception as ex:
        print(f"    not computable: {ex!r}")

    if rn:
        print(f"\n== Post hoc (gates nothing): a second roundoff realization of the new points at h = {H_FIN:g} (p*_h bisected from another bracket), {rn}")
        R2 = {}
        for nm, spec, tg in POINTS4:
            e2 = r.load_point(os.path.join(rn, "runs"), spec, H_FIN)
            e2["name"] = nm
            R2[(nm, H_FIN)] = e2
        derive4(R2)
        for nm, spec, tg in POINTS4:
            e2, e = R2[(nm, H_FIN)], ok.get((nm, H_FIN))
            if "err" in e2 or e is None:
                print(f"    {nm}: not available ({e2.get('err')})")
                continue
            M2, g2, r2, d2 = e2["SS"]["plat1"], e2["mem"]["g1"], e2["mem"]["r_th"], e2["sp"]
            print(f"    {nm:5s} p*_h difference {e2['p']['ps'] - e['p']['ps']:+.1e}  M(rn) - M {M2 - e['SS']['plat1']:+.1e}  gamma_1, r_th, rdot_th differences {g2 - e['mem']['g1']:+.1e}, {r2 - e['mem']['r_th']:+.1e}, {d2 - e['sp']:+.1e}  "
                  f"M - P4 (rn) {M2 - s3.predict('P4', g2, r2, d2, H_FIN):+.2e}")
    return finish(t0)


MUTATIONS = {      # name: the checks it newly fails (those failing without a mutation, --baseline, excluded); computed on the development set (raw I305-mac dev/mut_a, baseline dev/ana_dev3.log: FITQ, P4 and P4-ph fail) before being declared in its README
    "asym_sign": ['FRZ', 'PMEM', 'PTD', 'SEPC', 'SEPN', 'SEPW', 'TGT'],
    "bel_swap": ['PBEL'],
    "bisect_early": ['BR'],
    "coef_a": ['PPRED4'],
    "coef_c": ['PPRED4'],
    "crit_dm": ['CFG'],
    "dep_eps": ['FRZ', 'PDEP', 'PTD'],
    "drop_point": ['FRZ'],
    "eps_h": ['FRZ', 'HST', 'LEV', 'PLEV', 'PTD', 'SEPB', 'SEPC', 'SEPN', 'SEPW'],
    "evt_side": ['EVT'],
    "fitq_scale": [],
    "flux_sign": ['RC4'],
    "frz_a": ['FRZ'],
    "frz_m0": ['FRZ', 'PPRED4'],
    "lev_bar": ['LEV'],
    "lev_design": ['PLEV'],
    "lev_rdot": [],
    "loo_off": ['FRZ'],
    "mass_exponent": ['RC4', 'REP4'],
    "no_rscale": ['PPRED4'],
    "no_run": ['HST', 'LEV', 'RH', 'RUNS', 'SEPB', 'SEPC', 'SEPN', 'SEPW', 'TGT'],
    "no_settle": ['RC5', 'RC8', 'REP4'],
    "normal_drop": [],
    "outr_range": ['OUTR'],
    "rdot_nolapse": ['FRZ', 'PTD'],
    "rep_off": ['REP4'],
    "rh_swap": ['RH'],
    "sep_k": ['PSEP3'],
    "sep_self": ['SEPB', 'SEPC', 'SEPN', 'SEPW'],
    "short_limit": ['RC5', 'RC8', 'REP4'],
    "spec_w": ['CFG'],
    "stl_flat": ['STL'],
    "tgt_scale": ['TGT'],
    "w_scale": ['FAM0', 'HST', 'LEV', 'OUTR', 'REP4', 'RH', 'RUNS', 'SEPB', 'SEPC', 'SEPN', 'SEPW', 'TGT'],
    "xc_offset": ['FAM0', 'HST', 'LEV', 'OUTR', 'REP4', 'RH', 'RUNS', 'SEPB', 'SEPC', 'SEPN', 'SEPW', 'TGT'],
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
        s3.MUT = r.MUT = b.MUT = m.MUT = g.MUT = sm.MUT = MUT
        dn.MUT = MUT if MUT in ("theta_coef", "ghost_sign", "scalar_coef") else None
        tl.LMUT = MUT if MUT in ("asym_sign",) else None
        if MUT == "crit_dm":
            sm.CRIT["dm"] = 2.0
    import platform
    import numba as nb
    print(f"settled_mass_route4.py {cmd} {opts}: numba {nb.__version__}, numpy {np.__version__}, python {platform.python_version()}, {platform.node()}, "
          f"mutation {MUT}, {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    if cmd == "tag":
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
        s3.check_planted()
        pts15, _ = s3.points15()
        check_plev(pts15)
        check_pbel()
        sys.exit(finish(t0))
    elif cmd == "analyze":
        sys.exit(analyze_cmd(opts["raw"], opts.get("rn")))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()

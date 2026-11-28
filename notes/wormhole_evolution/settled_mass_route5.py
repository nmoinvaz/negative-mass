#!/usr/bin/env python3
"""The frozen settled-mass refit tested blind on a route of moderate leverage, 0.4 to 1.0 in its design (#313). D = 4, delta = 1e-8.

#295 froze (a, b, c) = (0.2296318924, 0.2834170383, 0.7657657062) of y = M / r_th - M_0 = a gamma_1 + b gamma_1^2 + c rdot_th over fifteen members at h = 0.005. #295's third route lay inside the
fitted data (leverage 0.07 to 0.33 against a fitted maximum of 0.393); #305's route lay far outside it (leverage 0.27 to 3.06, rdot_th down to -2.7e-3) and the refit missed it by 9e-5 to 2.1e-4.
This driver tests the same frozen coefficients on a Gaussian normal pulse (x_c, w) selected in development at h = 0.02 by the issue's rule (raw README stage 1, scripts/devselect5.py): every
leverage at most 1.0, the largest at least 0.43, the largest summed leverage, and runs the chain as #283, #295 and #305.

A thin driver. It imports settled_mass_route4.py (#305), which imports settled_mass_route3.py (#295), settled_mass_route.py (#283), settled_mass_b.py (#266), settled_mass_member.py (#267)
and everything below them, all unchanged.

    python3 settled_mass_route5.py tag|bracket|pstar|job|member ...                  as settled_mass_route.py (#283), for the route's points
    python3 settled_mass_route5.py planted [--mutate=NAME]                           planted cases
    python3 settled_mass_route5.py analyze --raw=DIR [--mutate=NAME] [--baseline=FILE] [--rn=DIR]

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
import settled_mass_route4 as s4           # noqa: E402  (#305's driver: the leverage, the frozen design; imports #295's, #283's, #266's, #267's unchanged)
s3 = s4.s3
r = s3.r
b, m, g, th, tl, sm, dn, td = s3.b, s3.m, s3.g, s3.th, s3.tl, s3.sm, s3.dn, s3.td

FAM = "G"
H2 = r.H2                                       # (0.01, 0.005)
H_FIN = r.H_FIN
MUT = None
ND = 2
I295 = s4.I295
I305 = os.path.join(s3.RAW_ROOT, "I305-mac")
I283 = s3.I283
M0 = s3.M0
CRIT3 = s3.CRIT3

POINTS5 = (("Y1", (0.0394, 5.0, 0.7), 0.012), ("Y2", (0.0768, 5.0, 0.7), 0.048), ("Y3", (0.09431, 5.0, 0.7), 0.075))      # (name, (B, xc, w), target gamma_1): the route selected by the rule of the raw README (stage 1)
POINTS_X = s4.POINTS4                           # #305's route X, stored in I305-mac
REP5 = ("RX2", (0.1695, 4.0, 1.5))              # #305's X2, reproduced at h = 0.005 from its stored bracket
TARGETS = (0.012, 0.048, 0.075)
CRIT5 = dict(
    rdot_range=(-6.7e-4, 1.6e-3),               # the issue's range of the fitted rdot_th (information only here)
    lev_lo=0.40, lev_hi=1.0,                    # LEVB: the issue's 0.4 to 1.0 at h = 0.005 (the fitted maximum is 0.3932)
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
# the leverage band
# --------------------------------------------------------------------------------------------------------------------
def band_ok(L):
    """LEVB's decision: every leverage at most hi and the largest at least lo (mutations lev_hi: hi = 0.5; lev_lo: lo = 1.0)."""
    lo = 1.0 if mut("lev_lo") else CRIT5["lev_lo"]
    hi = 0.5 if mut("lev_hi") else CRIT5["lev_hi"]
    return len(L) > 0 and max(L) <= hi and max(L) >= lo


def check_pband():
    """PBAND: the LEVB decision on planted leverage triples (the raw README stage 1)."""
    cases = (((0.1, 0.2, 0.39), False), ((0.1, 0.2, 0.41), True), ((0.2, 0.5, 1.0), True), ((0.2, 0.5, 1.01), False), ((0.45, 0.2, 0.1), True))
    bad = [L for L, want in cases if band_ok(L) != want]
    check("PBAND", "the LEVB decision (every leverage <= 1.0 and the largest >= 0.40) on planted triples: (0.1, 0.2, 0.39) false, (0.1, 0.2, 0.41) true, (0.2, 0.5, 1.0) true, (0.2, 0.5, 1.01) false, (0.45, 0.2, 0.1) true",
          not bad, f"wrong: {bad}" if bad else f"{len(cases)} planted cases")


# --------------------------------------------------------------------------------------------------------------------
# #305's route X: curve and members (stored runs read only)
# --------------------------------------------------------------------------------------------------------------------
def x_read(h):
    out = []
    for nm, spec, _ in POINTS_X:
        e = r.load_point(os.path.join(I305, "runs"), spec, h)
        x = r.read_point(e)
        out.append((nm, e, x))
    return out


def route_curve_X(h):
    """#305's route X at h (B = 0 as the anchor at gamma_1 = 0): gamma_1, r_th, rdot_th, ascending in gamma_1."""
    anchor = r.fam_b(h)
    g1, rth, rdot = [float(anchor["g1"][0])], [float(anchor["rth"][0])], [float(anchor["rdot"][0])]
    for nm, e, x in x_read(h):
        g1.append(x["mem"]["g1"]), rth.append(x["mem"]["r_th"]), rdot.append(x["sp"])
    o = np.argsort(g1)
    return dict(g1=np.array(g1)[o], rth=np.array(rth)[o], rdot=np.array(rdot)[o])


def points_X(h=H_FIN):
    """#305's three stored members at h: name, gamma_1, r_th, rdot_th, M."""
    return [dict(name=nm, route="X", g1=x["mem"]["g1"], rth=x["mem"]["r_th"], rdot=x["sp"], M=x["M"]) for nm, e, x in x_read(h)]


# --------------------------------------------------------------------------------------------------------------------
# the analysis
# --------------------------------------------------------------------------------------------------------------------
def derive5(D):
    for e in D.values():
        if mut("no_run") and e.get("name") == POINTS5[0][0] and e["h"] == H_FIN and "err" not in e:
            e["err"] = "member output treated as missing"
    return r.derive(D)


def declared_spec(nm):
    spec = dict((n, s) for n, s, _ in POINTS5)[nm] if nm != REP5[0] else REP5[1]
    return (spec[0], spec[1], spec[2] * (1.1 if mut("spec_w") else 1.0))


def analyze_cmd(raw, rn=None):
    t0 = time.time()
    root = os.path.join(raw, "runs")
    D = {}
    for nm, spec, tg in POINTS5:
        for h in H2:
            e = r.load_point(root, spec, h)
            e["name"], e["target"] = nm, tg
            D[(nm, h)] = e
    rep = r.load_point(os.path.join(raw, "runs_rep"), REP5[1], H_FIN)
    rep["name"], rep["target"] = REP5[0], None
    D[(REP5[0], H_FIN)] = rep
    miss = {f"{k[0]} h{k[1]:g}": e["err"] for k, e in D.items() if "err" in e}
    print(f"== Runs read from {root} (and runs_rep): {len(D) - len(miss)} of {len(D)} (point, h) complete" + (f"; missing: {miss}" if miss else ""))
    derive5(D)
    for k, e in D.items():
        if "err" in e and f"{k[0]} h{k[1]:g}" not in miss:
            miss[f"{k[0]} h{k[1]:g}"] = e["err"]
            print(f"    {k[0]} h {k[1]:g}: {e['err']}")
    ok = {k: e for k, e in D.items() if "err" not in e}
    names = [nm for nm, _, _ in POINTS5]
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
    s4.check_plev(pts15)
    check_pband()

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
    for nm, spec, tg in POINTS5:
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

    print("\n== Leverage in the frozen fit's design and rdot_th against the fitted range (h = 0.005)")
    A15, G15, hat15 = s4.fitted_design(pts15)
    Lmax = float(hat15.max())
    kmax = int(np.argmax(hat15))
    lev = {}
    for nm in names:
        e = ok.get((nm, H_FIN))
        if e is None:
            continue
        gg, rd = e["mem"]["g1"], e["sp"] * (0.0 if mut("lev_rdot") else 1.0)
        lev[nm] = s4.leverage(gg, rd, A15, G15)
        print(f"    {nm}: gamma_1 {gg:.6f}, rdot_th {e['sp']:+.4e}: leverage {lev[nm]:.4f} ({'above' if lev[nm] > Lmax else 'within'} the fitted maximum {Lmax:.4f}, {pts15[kmax]['name']}); "
              f"rdot_th {'outside' if (e['sp'] < CRIT5['rdot_range'][0] or e['sp'] > CRIT5['rdot_range'][1]) else 'inside'} [{CRIT5['rdot_range'][0]:+.1e}, {CRIT5['rdot_range'][1]:+.1e}] (information)")
    check("LEVB", "at h = 0.005, in the frozen fit's design (gamma_1, gamma_1^2, rdot_th), every new member's leverage is at most 1.0 and the largest is at least 0.40 (the fitted maximum is 0.3932, W3)",
          len(lev) == len(names) and band_ok(list(lev.values())), ", ".join(f"{nm} {v:.4f}" for nm, v in lev.items()))

    print("\n== Separation from the earlier routes at the same gamma_1 (h = 0.005; B = 0 anchors each curve; piecewise-linear and PCHIP, the smaller difference; gamma_1 up to 3 % beyond a curve's last point is extrapolated)")
    refs = {"B": r.fam_b(H_FIN), "W": s3.route_curve(H_FIN, ("W1", "W2", "W3")), "C": s3.route_curve(H_FIN, ("C1", "C2", "C3")), "N": s4.route_curve_N(H_FIN), "X": route_curve_X(H_FIN)}
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
            if mut("sep_self") or (mut("sepx_self") and k_ == "X"):
                gg = c_["mem"]["g1"]
                cvu = dict(g1=np.array([0.0, gg - 1e-4, gg, gg + 1e-4]), rth=np.array([c_["mem"]["r_th"]] * 4), rdot=np.array([c_["sp"]] * 4))
            s_ = s3.sep3(c_["mem"]["g1"], (c_["mem"]["r_th"], c_["sp"]), (a_["mem"]["r_th"], a_["sp"]), cvu)
            SEP[(nm, k_)] = s_
            sepacc[k_].add(s_["informative"])
            q = s_["q"]
            print(f"    {nm} against route {k_}: gamma_1 {c_['mem']['g1']:.5f} (beyond its last point by {s_['ext']:.1%}); r_th - route: {q['r_th']['d_lin']:+.2e} (linear), {q['r_th']['d_pchip']:+.2e} (PCHIP), its h difference {q['r_th']['hdiff']:.1e}"
                  f" -> {'ok' if q['r_th']['ok'] else 'no'}; rdot_th - route: {q['rdot']['d_lin']:+.2e}, {q['rdot']['d_pchip']:+.2e}, h difference {q['rdot']['hdiff']:.1e} -> {'ok' if q['rdot']['ok'] else 'no'}: "
                  f"{'informative' if s_['informative'] else 'NOT informative'}")
    for k_, nmk in (("W", "SEPW"), ("C", "SEPC"), ("B", "SEPB"), ("N", "SEPN"), ("X", "SEPX")):
        what = {"W": "route W (#283)", "C": "route C (#283)", "B": "family B (#266)", "N": "route N (#295)", "X": "route X (#305)"}[k_]
        check(nmk, f"the new members' r_th or rdot_th at departure differs from {what}'s at the same gamma_1 by more than 5 times its own h difference, all three points", sepacc[k_].ok)

    print("\n== The frozen refit on the new members: M - predicted, h = 0.005 (gating), with the leverage, the static fit's rms and rdot_th beside each residual")
    print("    P4: M = r_th (M_0 + a g + b g^2 + c rdot_th), M_0 = 0.223086, (a, b, c) = (0.2296318924, 0.2834170383, 0.7657657062), frozen in #295 before any run of its third route, and before this issue")
    print("    point   leverage  M(plateau)   P4 predicted  M - P4      rms(fit)  rdot_th     | comparison, gates nothing: M - P3 (#267's)   M - P2' (refit without the speed)")
    res = {}
    for nm in names:
        e = ok.get((nm, H_FIN))
        if e is None:
            continue
        M, g1, rth, rd = e["SS"]["plat1"], e["mem"]["g1"], e["mem"]["r_th"], e["sp"]
        res[nm] = M - s3.predict("P4", g1, rth, rd)
        print(f"    {nm:6s}  {lev.get(nm, float('nan')):.4f}   {M:.6f}   {s3.predict('P4', g1, rth, rd):.6f}      {res[nm]:+.2e}   {e['rms']:.2e}  {rd:+.3e}  | {M - s3.predict('P3', g1, rth, rd):+.2e}                {M - s3.predict('P2p', g1, rth, rd):+.2e}")
    if len(res) == len(names) and names:
        worst = max(abs(v) for v in res.values())
        check("P4", "the frozen refit predicts the settled mass at every new member within 6e-5 (h = 0.005; coefficients frozen in #295's raw README before any run of its third route, and before this issue)", worst <= CRIT3["bar_pred"],
              f"{len(res)} points; max |M - predicted| {worst:.1e} (bar {CRIT3['bar_pred']:g}); " + ", ".join(f"{nm} {res[nm]:+.1e} (leverage {lev.get(nm, float('nan')):.2f}, rms {ok[(nm, H_FIN)]['rms']:.1e})" for nm in names))
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

    print("\n== Reproduction of #305's X2 at h = 0.005 from its stored bracket, bit for bit (stored runs of #305 read only)")
    e = ok.get((REP5[0], H_FIN))
    if e is None:
        check("REP5", "#305's X2 (B, xc, w) = (0.1695, 4.0, 1.5) at h = 0.005 reproduces #305's stored run bit for bit", False, "run missing")
    else:
        ref = r.load_point(os.path.join(I305, "runs"), (0.09112, 4.0, 1.5) if mut("rep_off") else REP5[1], H_FIN)
        er = ref if "err" not in ref else None
        if er is None:
            check("REP5", "#305's X2 reproduces #305's stored run bit for bit", False, ref["err"])
        else:
            rm = r.read_point(er)
            same_p = e["p"]["ps"] == er["p"]["ps"]
            same_s = all(np.array_equal(e["s"][k], er["s"][k]) for k in er["s"]) and set(e["s"]) == set(er["s"])
            same_l = all(np.array_equal(e["loc"][k], er["loc"][k]) for k in er["loc"]) and set(e["loc"]) == set(er["loc"])
            same_m = {k: v for k, v in e["m"].items() if k != "secs"} == {k: v for k, v in er["m"].items() if k != "secs"}      # all but the wall time
            dM = e["SS"]["plat1"] - rm["M"]
            dd = (e["mem"]["g1"] - rm["mem"]["g1"], e["mem"]["r_th"] - rm["mem"]["r_th"], e["sp"] - rm["sp"], e["mem"]["t_dep"] - rm["mem"]["t_dep"])
            check("REP5", "#305's X2 at h = 0.005 reproduces #305's stored run bit for bit: p*_h, the chain's series, the locus, the member rows (all but the wall time), the plateau, and (gamma_1, r_th, rdot_th, t_dep)",
                  same_p and same_s and same_l and same_m and dM == 0.0 and all(x == 0.0 for x in dd),
                  f"p*_h equal {same_p}, series {same_s}, locus {same_l}, member rows {same_m}; dM {dM:+.1e}, d(gamma_1, r_th, rdot_th, t_dep) {dd[0]:+.1e}, {dd[1]:+.1e}, {dd[2]:+.1e}, {dd[3]:+.1e}")

    # ---------------------------------------------------------------- post hoc blocks (declared in the raw README stage 1; they gate nothing)
    print("\n== Post hoc (gates nothing): the nine held-out members, miss against leverage, the static fit's rms and the speed")
    ptsN, ptsX = s4.points_N(H_FIN), points_X(H_FIN)
    rms_of = {}
    for nm_, spec_, _ in s3.POINTS3:
        e_ = r.load_point(os.path.join(I295, "runs"), spec_, H_FIN)
        rms_of[nm_] = s4.rms_at(e_, r.read_point(e_)["mem"]["t_dep"])
    for nm_, e_, x_ in x_read(H_FIN):
        rms_of[nm_] = s4.rms_at(e_, x_["mem"]["t_dep"])
    for nm_ in names:
        if (nm_, H_FIN) in ok:
            rms_of[nm_] = ok[(nm_, H_FIN)]["rms"]
    held = []
    for grp, pl in (("N", ptsN), ("X", ptsX), ("Y", [dict(name=nm_, g1=ok[(nm_, H_FIN)]["mem"]["g1"], rth=ok[(nm_, H_FIN)]["mem"]["r_th"], rdot=ok[(nm_, H_FIN)]["sp"], M=ok[(nm_, H_FIN)]["SS"]["plat1"]) for nm_ in names if (nm_, H_FIN) in ok])):
        for p_ in pl:
            held.append(dict(route=grp, name=p_["name"], g1=p_["g1"], rth=p_["rth"], rdot=p_["rdot"], M=p_["M"], lev=s4.leverage(p_["g1"], p_["rdot"], A15, G15),
                             miss=p_["M"] - s3.predict("P4", p_["g1"], p_["rth"], p_["rdot"]), rms=rms_of[p_["name"]]))
    print("    name  leverage  gamma_1    r_th       rdot_th      M - P4      rms(fit)   |M - P4|/rms  rdot_th vs fitted range")
    for h_ in held:
        inr = CRIT5["rdot_range"][0] <= h_["rdot"] <= CRIT5["rdot_range"][1]
        print(f"    {h_['name']:4s}  {h_['lev']:.3f}    {h_['g1']:.5f}  {h_['rth']:.6f}  {h_['rdot']:+.3e}  {h_['miss']:+.2e}  {h_['rms']:.2e}   {abs(h_['miss']) / h_['rms']:.3f}        {'inside' if inr else 'outside'}")
    byl = sorted(held, key=lambda x: x["lev"])
    print("    ordered by leverage: " + ", ".join(f"{x['name']} {x['lev']:.2f}:{x['miss']:+.1e}" for x in byl))
    print("    the largest |M - P4| among the held-out members whose leverage is within the band [0.40, 1.0]: " + (", ".join(f"{x['name']} {abs(x['miss']):.2e}" for x in byl if 0.40 <= x["lev"] <= 1.0) or "none"))

    try:
        print("\n== Post hoc (gates nothing): the static-member fit's rms at departure (FITQ's quantity) at every stored member at h = 0.005, against |rdot_th| (#305's M1, with the new members)")
        rows = []
        R266 = m.load266(H_FIN)
        for B_ in (0.0,) + tuple(s3.B_266):
            rows.append((f"B{B_:g}", R266[B_]["mem"]["g1"], R266[B_]["mem"]["r_th"], R266[B_]["speed"], s4.rms_at(R266[B_]["e"], R266[B_]["mem"]["t_dep"])))
        for ee in s3.E_267:
            base = os.path.join(s3.I267, "runs", m.TAG, f"h{H_FIN:g}", f"e{ee:g}")
            d = dict(e=ee, h=H_FIN, m=json.load(open(base + "/member/member.json")), loc=dict(np.load(base + "/member/locus.npz")))
            mem_ = b.member_at(d, CRIT3.get("eps", 1e-3), "m")
            rows.append((f"e{ee}", mem_["g1"], mem_["r_th"], m.speed_at(d["loc"], mem_["t_dep"]), s4.rms_at(d, mem_["t_dep"])))
        for nm_, spec_ in s3.POINTS_283:
            e_ = r.load_point(os.path.join(I283, "runs"), spec_, H_FIN)
            x_ = r.read_point(e_)
            rows.append((nm_, x_["mem"]["g1"], x_["mem"]["r_th"], x_["sp"], s4.rms_at(e_, x_["mem"]["t_dep"])))
        for h_ in held:
            rows.append((h_["name"], h_["g1"], h_["rth"], h_["rdot"], h_["rms"]))
        for nm_, g1_, rth_, sp_, rms_ in rows:
            print(f"    {nm_:6s} gamma_1 {g1_:.5f}  r_th - 1 {rth_ - 1:+.5f}  rdot_th {sp_:+.3e}  rms(fit) {rms_:.2e}{'  (above FITQ bar 1e-4)' if rms_ > 1e-4 else ''}")
        arr = np.array([[abs(x[1]), abs(x[2] - 1), abs(x[3]), x[4]] for x in rows])
        print(f"    correlation of rms with |gamma_1|, |r_th - 1|, |rdot_th| over {len(rows)} members: " + ", ".join(f"{np.corrcoef(arr[:, j], arr[:, 3])[0, 1]:+.3f}" for j in range(3)))
        a22 = np.array([[abs(x[1]), abs(x[2] - 1), abs(x[3]), x[4]] for x in rows if not x[0].startswith("Y")])
        print(f"    the same over the {len(a22)} members without Y (#305's note: +0.952 with |rdot_th|, +0.795 with |r_th - 1|, +0.549 with |gamma_1|): " + ", ".join(f"{np.corrcoef(a22[:, j], a22[:, 3])[0, 1]:+.3f}" for j in range(3)))
        print("    rms(fit) at the same gamma_1 target on each route (the static fit's rms against gamma_1 within a route):")
        for rt, nms in (("W", ("W1", "W2", "W3")), ("C", ("C1", "C2", "C3")), ("N", ("N1", "N2", "N3")), ("X", ("X1", "X2", "X3")), ("Y", ("Y1", "Y2", "Y3"))):
            d_ = {x[0]: x for x in rows}
            if all(n_ in d_ for n_ in nms):
                print(f"        route {rt}: " + ", ".join(f"{n_} (gamma_1 {d_[n_][1]:.4f}, rdot_th {d_[n_][3]:+.1e}) {d_[n_][4]:.2e}" for n_ in nms))
    except Exception as ex:
        print(f"    not computable: {ex!r}")

    print("\n== Post hoc (gates nothing): #305's description of its miss (M - P4 = k rdot_th, k the mean of X's three) applied to the new members without refitting k")
    try:
        kx = [x["miss"] / x["rdot"] for x in held if x["route"] == "X"]
        k_ = float(np.mean(kx))
        print(f"    k from X: {', '.join(f'{x:+.4f}' for x in kx)}; mean {k_:+.4f}; c_eff/c on X: " + ", ".join(f"{1.0 + x / CRIT3['P4'][2]:.3f}" for x in kx))
        better = []
        for x in held:
            if x["route"] != "Y":
                continue
            pred = k_ * x["rdot"]
            better.append(abs(x["miss"] - pred) < abs(x["miss"]))
            print(f"    {x['name']}: rdot_th {x['rdot']:+.3e}; predicted miss k rdot_th {pred:+.2e}; observed {x['miss']:+.2e}; error of the description {x['miss'] - pred:+.2e} against {x['miss']:+.2e} with no correction: "
                  f"{'better' if better[-1] else 'not better'}; (M - P4)/rdot_th {x['miss'] / x['rdot']:+.4f}; c_eff/c {1.0 + x['miss'] / x['rdot'] / CRIT3['P4'][2]:.3f}")
        print(f"    the description is better than no correction at {sum(better)} of {len(better)} new members (declared in the raw README stage 1: better if |observed - k rdot_th| < |observed| at all three)")
        ky = [x["miss"] / x["rdot"] for x in held if x["route"] == "Y"]
        print(f"    (M - P4)/rdot_th on Y alone: {', '.join(f'{x:+.4f}' for x in ky)}; mean {np.mean(ky):+.4f}, spread {(max(ky) - min(ky)) / abs(np.mean(ky)):.1%}; on X: mean {k_:+.4f}, spread {(max(kx) - min(kx)) / abs(k_):.1%}; Y's mean over X's {np.mean(ky) / k_:.3f}")
        print("    (M - P4)/rdot_th on N: " + ", ".join(f"{x['name']} {x['miss'] / x['rdot']:+.4f}" for x in held if x["route"] == "N") + " (descriptions of single routes, three points each; not a law)")
    except Exception as ex:
        print(f"    not computable: {ex!r}")

    print("\n== Post hoc (gates nothing): refits including the new points (descriptions, not tests)")
    try:
        newp = [dict(name=x["name"], route="Y", g1=x["g1"], rth=x["rth"], rdot=x["rdot"], M=x["M"]) for x in held if x["route"] == "Y"]
        for label, base_ in (("the fifteen and Y", list(pts15)), ("the fifteen, #295's three and Y", list(pts15) + ptsN), ("the fifteen, #295's three, #305's three and Y", list(pts15) + ptsN + ptsX)):
            allp = base_ + newp
            g_, r_, d_, M_ = (np.array([p[k] for p in allp]) for k in ("g1", "rth", "rdot", "M"))
            co2 = s3.fit(g_, r_, d_, M_, "P4")
            eM, _ = s3.loo(g_, r_, d_, M_, "P4")
            ins = M_ - r_ * (M0 + s3.design(g_, d_, "P4") @ co2)
            print(f"    refit of P4 over {label} ({len(allp)} points): coefficients {', '.join(f'{x:+.4f}' for x in co2)} (frozen {', '.join(f'{x:+.4f}' for x in CRIT3['P4'])}); in-sample max {np.max(np.abs(ins)):.1e}; "
                  f"leave-one-out max {np.max(np.abs(eM)):.1e}; Y's leave-one-out: " + ", ".join(f"{p['name']} {x:+.1e}" for p, x in zip(allp, eM) if p["route"] == "Y"))
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

    if rn:
        print(f"\n== Post hoc (gates nothing): a second roundoff realization of the new points at h = {H_FIN:g} (p*_h bisected from another bracket), {rn}")
        R2 = {}
        for nm, spec, tg in POINTS5:
            e2 = r.load_point(os.path.join(rn, "runs"), spec, H_FIN)
            e2["name"] = nm
            R2[(nm, H_FIN)] = e2
        derive5(R2)
        for nm, spec, tg in POINTS5:
            e2, e = R2[(nm, H_FIN)], ok.get((nm, H_FIN))
            if "err" in e2 or e is None:
                print(f"    {nm}: not available ({e2.get('err')})")
                continue
            M2, g2, r2, d2 = e2["SS"]["plat1"], e2["mem"]["g1"], e2["mem"]["r_th"], e2["sp"]
            print(f"    {nm:5s} p*_h difference {e2['p']['ps'] - e['p']['ps']:+.1e}  M(rn) - M {M2 - e['SS']['plat1']:+.1e}  gamma_1, r_th, rdot_th differences {g2 - e['mem']['g1']:+.1e}, {r2 - e['mem']['r_th']:+.1e}, {d2 - e['sp']:+.1e}  "
                  f"M - P4 (rn) {M2 - s3.predict('P4', g2, r2, d2, H_FIN):+.2e}")
    return finish(t0)


MUTATIONS = {      # name: the checks it newly fails (those failing without a mutation, --baseline, excluded); computed on the development set (raw I313-mac dev/mut_a, baseline dev/ana_dev2.log: FITQ fails) before being declared in its README stage 2
    "asym_sign": ['FRZ', 'P4', 'P4-ph', 'PMEM', 'PTD', 'SEPC', 'SEPN', 'SEPW', 'SEPX', 'TGT'],
    "bisect_early": ['BR'],
    "coef_a": ['P4', 'P4-ph', 'PPRED4'],
    "coef_c": ['P4', 'P4-ph', 'PPRED4'],
    "crit_dm": ['CFG'],
    "dep_eps": ['FRZ', 'P4', 'P4-ph', 'PDEP', 'PTD'],
    "drop_point": ['FRZ'],
    "eps_h": ['FRZ', 'HST', 'LEVB', 'P4', 'PLEV', 'PTD', 'SEPB', 'SEPC', 'SEPN', 'SEPW', 'SEPX'],
    "evt_side": ['EVT'],
    "fitq_scale": [],
    "flux_sign": ['RC4'],
    "frz_a": ['FRZ'],
    "frz_m0": ['FRZ', 'PPRED4'],
    "lev_design": ['PLEV'],
    "lev_hi": ['LEVB', 'PBAND'],
    "lev_lo": ['LEVB', 'PBAND'],
    "lev_rdot": [],
    "loo_off": ['FRZ'],
    "mass_exponent": ['P4', 'P4-ph', 'RC4', 'REP5'],
    "no_rscale": ['P4', 'P4-ph', 'PPRED4'],
    "no_run": ['HST', 'LEVB', 'P4', 'RH', 'RUNS', 'SEPB', 'SEPC', 'SEPN', 'SEPW', 'SEPX', 'TGT'],
    "no_settle": ['P4', 'P4-ph', 'RC5', 'RC8', 'REP5'],
    "normal_drop": [],
    "rdot_nolapse": ['FRZ', 'PTD'],
    "rep_off": ['REP5'],
    "rh_swap": ['RH'],
    "sep_k": ['PSEP3'],
    "sep_self": ['SEPB', 'SEPC', 'SEPN', 'SEPW', 'SEPX'],
    "sepx_self": ['SEPX'],
    "short_limit": ['P4', 'P4-ph', 'RC5', 'RC8', 'REP5'],
    "spec_w": ['CFG'],
    "stl_flat": ['STL'],
    "tgt_scale": ['TGT'],
    "w_scale": ['FAM0', 'HST', 'LEVB', 'P4', 'REP5', 'RH', 'RUNS', 'SEPB', 'SEPC', 'SEPN', 'SEPW', 'SEPX', 'TGT'],
    "xc_offset": ['FAM0', 'HST', 'LEVB', 'P4', 'REP5', 'RH', 'RUNS', 'SEPB', 'SEPC', 'SEPN', 'SEPW', 'SEPX', 'TGT'],
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
        s4.MUT = s3.MUT = r.MUT = b.MUT = m.MUT = g.MUT = sm.MUT = MUT
        dn.MUT = MUT if MUT in ("theta_coef", "ghost_sign", "scalar_coef") else None
        tl.LMUT = MUT if MUT in ("asym_sign",) else None
        if MUT == "crit_dm":
            sm.CRIT["dm"] = 2.0
    import platform
    import numba as nb
    print(f"settled_mass_route5.py {cmd} {opts}: numba {nb.__version__}, numpy {np.__version__}, python {platform.python_version()}, {platform.node()}, "
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
        s4.check_plev(pts15)
        check_pband()
        sys.exit(finish(t0))
    elif cmd == "analyze":
        sys.exit(analyze_cmd(opts["raw"], opts.get("rn")))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()

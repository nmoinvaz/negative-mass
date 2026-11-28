#!/usr/bin/env python3
"""Post hoc reading of #313 (written after the analysis of record; reads the stored runs, changes none, gates nothing).

usage: python3 settled_mass_route5_post.py [RAWDIR]     (RAWDIR defaults to ~/Source/negative-mass-raw/I313-mac)

Prints, for the nine held-out members (#295's N, #305's X, this issue's Y) at h = 0.005: the miss M - P4 of the frozen refit against leverage, |rdot_th|, gamma_1, |r_th - 1| and the static fit's rms
(Spearman rank correlation of |miss|), the geometric mean of |miss|/rms per route, the margin of each new member to the bar against the changes of M with h and with the roundoff realization,
and the smallest separation from the earlier routes in units of the point's own h difference."""
import math
import os
import sys

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
raw = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Source/negative-mass-raw/I313-mac")
sys.argv = ["x"]
import settled_mass_route5 as s5          # noqa: E402

s4, s3, r = s5.s4, s5.s3, s5.r
fin, BAR = s5.H_FIN, s5.CRIT3["bar_pred"]
D = {}
for nm, spec, tg in s5.POINTS5:
    for h in s5.H2:
        e = r.load_point(os.path.join(raw, "runs"), spec, h)
        e["name"] = nm
        D[(nm, h)] = e
r.derive(D)
pts, c0 = s3.points15()
A15, G15, hat15 = s4.fitted_design(pts)
held = []
for kind, pl in (("N", s4.points_N(fin)), ("X", s5.points_X(fin))):
    for p in pl:
        held.append(dict(nm=p["name"], kind=kind, g1=p["g1"], rth=p["rth"], rd=p["rdot"], M=p["M"]))
rms_of = {}
for nm_, spec_, _ in s3.POINTS3:
    e_ = r.load_point(os.path.join(s5.I295, "runs"), spec_, fin)
    rms_of[nm_] = s4.rms_at(e_, r.read_point(e_)["mem"]["t_dep"])
for nm_, e_, x_ in s5.x_read(fin):
    rms_of[nm_] = s4.rms_at(e_, x_["mem"]["t_dep"])
for nm, spec, tg in s5.POINTS5:
    e = D[(nm, fin)]
    held.append(dict(nm=nm, kind="Y", g1=e["mem"]["g1"], rth=e["mem"]["r_th"], rd=e["sp"], M=e["SS"]["plat1"]))
    rms_of[nm] = e["rms"]
for p in held:
    p["lev"] = s4.leverage(p["g1"], p["rd"], A15, G15)
    p["miss"] = p["M"] - s3.predict("P4", p["g1"], p["rth"], p["rd"], fin, mutated=False)
    p["rms"] = rms_of[p["nm"]]

print("== the nine held-out members, h = 0.005 (M - P4 of the frozen refit)")
print("name  leverage  gamma_1   |r_th - 1|  rdot_th      rms(fit)  M - P4      |M - P4|/rms")
for p in held:
    print(f"{p['nm']:4s}  {p['lev']:.3f}    {p['g1']:.5f}  {abs(p['rth'] - 1):.5f}    {p['rd']:+.3e}  {p['rms']:.2e}  {p['miss']:+.2e}  {abs(p['miss']) / p['rms']:.3f}")

print("\n== Spearman rank correlation of |M - P4| with each quantity over the nine (and over the six of X and Y): descriptions of nine points, not tests")
am = np.array([abs(p["miss"]) for p in held])
for label, key in (("leverage", lambda p: p["lev"]), ("|rdot_th|", lambda p: abs(p["rd"])), ("gamma_1", lambda p: p["g1"]), ("|r_th - 1|", lambda p: abs(p["rth"] - 1)), ("static-fit rms", lambda p: p["rms"])):
    v = np.array([key(p) for p in held])
    s9 = spearmanr(am, v)[0]
    k6 = [i for i, p in enumerate(held) if p["kind"] in ("X", "Y")]
    s6 = spearmanr(am[k6], v[k6])[0]
    print(f"    {label:15s} nine {s9:+.3f}   six (X, Y) {s6:+.3f}")

print("\n== geometric mean of |M - P4| / rms per route (as #305's M2)")
gm = {}
for k in ("N", "X", "Y"):
    gm[k] = math.exp(np.mean([math.log(abs(p["miss"]) / p["rms"]) for p in held if p["kind"] == k]))
    print(f"    route {k}: {gm[k]:.3f}")
print(f"    Y over N {gm['Y'] / gm['N']:.1f}, X over Y {gm['X'] / gm['Y']:.1f}, X over N {gm['X'] / gm['N']:.1f}")

print("\n== the margin of each new member to the bar 6e-5 against the changes of M (h = 0.005)")
rn = os.path.join(raw, "runs_rn")
for nm, spec, tg in s5.POINTS5:
    e, e1 = D[(nm, fin)], D[(nm, 0.01)]
    p = [q for q in held if q["nm"] == nm][0]
    dh = abs(e["SS"]["plat1"] - e1["SS"]["plat1"])
    e2 = r.load_point(os.path.join(rn, "runs"), spec, fin)
    e2["name"] = nm
    r.derive({0: e2})
    drn = e2["SS"]["plat1"] - e["SS"]["plat1"]
    margin = BAR - abs(p["miss"])
    print(f"    {nm}: M - P4 {p['miss']:+.2e}; margin to the bar {margin:+.2e}; |M(h = 0.01) - M(h = 0.005)| {dh:.2e}; M(second realization) - M {drn:+.2e}; margin over the h change {margin / dh:.2f}, over the realization change {margin / abs(drn):.1f}")

print("\n== separation from the earlier routes at the same gamma_1 (the smaller of the linear and PCHIP differences) in units of the point's own h difference (the bar is 5), in r_th and in rdot_th; the last line is the smallest over all pairs of the larger of the two")
refs = {"B": r.fam_b(fin), "W": s3.route_curve(fin, ("W1", "W2", "W3")), "C": s3.route_curve(fin, ("C1", "C2", "C3")), "N": s4.route_curve_N(fin), "X": s5.route_curve_X(fin)}
refs["B"] = dict(g1=refs["B"]["g1"], rth=refs["B"]["rth"], rdot=refs["B"]["rdot"])
worst = (float("inf"), None)
for nm, spec, tg in s5.POINTS5:
    a_, c_ = D[(nm, s5.H2[0])], D[(nm, fin)]
    row = []
    for k_, cv in refs.items():
        s_ = s3.sep3(c_["mem"]["g1"], (c_["mem"]["r_th"], c_["sp"]), (a_["mem"]["r_th"], a_["sp"]), cv, k=0.0)
        q = s_["q"]
        ratios = [q["r_th"]["d"] / q["r_th"]["hdiff"], q["rdot"]["d"] / q["rdot"]["hdiff"]]
        row.append(f"{k_}: r_th {ratios[0]:.0f}, rdot_th {ratios[1]:.0f}")
        if max(ratios) < worst[0]:
            worst = (max(ratios), f"{nm} against {k_}")
    print(f"    {nm}: " + "; ".join(row))
print(f"    the smallest, over all (point, route), of the larger of the two ratios: {worst[0]:.0f} ({worst[1]})")

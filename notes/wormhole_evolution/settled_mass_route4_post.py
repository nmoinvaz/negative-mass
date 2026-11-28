#!/usr/bin/env python3
"""Post hoc reading of #305's misses (gates nothing; written after the analysis of record, from the same stored runs and the same functions as settled_mass_route4.py analyze).
usage: python3 settled_mass_route4_post.py [RAWDIR]     (RAWDIR defaults to ~/Source/negative-mass-raw/I305-mac)
Prints, for the three new members and #295's three held-out members at h = 0.005: the miss M - P4 divided by rdot_th, the coefficient c_eff of rdot_th that would remove it with a and b kept,
and the miss divided by its leverage and by the static fit's rms. Descriptions of three points; none is a test."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
raw = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Source/negative-mass-raw/I305-mac")
sys.argv = ["x"]
import settled_mass_route4 as s4          # noqa: E402

s3, r = s4.s3, s4.r
fin = s4.H_FIN
a_, b_, c_ = s3.CRIT3["P4"]
pts15, _ = s3.points15()
A15, G15, _ = s4.fitted_design(pts15)
D = {}
for nm, spec, tg in s4.POINTS4:
    e = r.load_point(os.path.join(raw, "runs"), spec, fin)
    e["name"] = nm
    D[(nm, fin)] = e
r.derive(D)
rows = []
for p in s4.points_N(fin):
    e = r.load_point(os.path.join(s4.I295, "runs"), dict((n, s) for n, s, _ in s3.POINTS3)[p["name"]], fin)
    rows.append(dict(nm=p["name"], g1=p["g1"], rth=p["rth"], rd=p["rdot"], M=p["M"], rms=s4.rms_at(e, r.read_point(e)["mem"]["t_dep"])))
for nm, spec, tg in s4.POINTS4:
    e = D[(nm, fin)]
    rows.append(dict(nm=nm, g1=e["mem"]["g1"], rth=e["mem"]["r_th"], rd=e["sp"], M=e["SS"]["plat1"], rms=s4.rms_at(e, e["mem"]["t_dep"])))
print("post hoc, h = 0.005, frozen P4 (a, b, c) = (%.10f, %.10f, %.10f); c_eff = c + (M - P4) / (r_th rdot_th): the coefficient of rdot_th that removes the miss, a and b kept" % (a_, b_, c_))
print("point  leverage  rdot_th     M - P4      (M - P4)/rdot_th   c_eff     c_eff/c    (M - P4)/leverage   (M - P4)/rms")
for q in rows:
    miss = q["M"] - s3.predict("P4", q["g1"], q["rth"], q["rd"], fin, mutated=False)
    lv = s4.leverage(q["g1"], q["rd"], A15, G15)
    ce = c_ + miss / (q["rth"] * q["rd"])
    q.update(miss=miss, ce=ce)
    print(f"{q['nm']:5s}  {lv:7.3f}  {q['rd']:+.3e}  {miss:+.2e}   {miss / q['rd']:+.4f}          {ce:.4f}    {ce / c_:.4f}    {miss / lv:+.2e}          {miss / q['rms']:+.3f}")
X = [q for q in rows if q["nm"].startswith("X")]
ratio = np.array([q["miss"] / q["rd"] for q in X])
print(f"new members: (M - P4)/rdot_th = {', '.join(f'{x:+.4f}' for x in ratio)}; mean {ratio.mean():+.4f}, spread (max - min)/|mean| {np.ptp(ratio) / abs(ratio.mean()):.3f}; "
      f"c_eff/c = {', '.join(f'{q['ce'] / c_:.3f}' for q in X)}")
print("*post hoc: three points and one constant; a description, not a test. The in-sample errors of the fifteen fitted points are 5.5e-5 at most, so a coefficient c_eff 10 % below c is not seen in them.*")

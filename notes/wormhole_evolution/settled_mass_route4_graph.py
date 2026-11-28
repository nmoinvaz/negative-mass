#!/usr/bin/env python3
"""Figure of #305: the frozen refit's errors against leverage in its design (the fifteen fitted points in-sample and leave-one-out, #295's three held-out members, the three members of the route chosen
for leverage), and where those members sit in the fit's (gamma_1, rdot_th) plane (h = 0.005).
usage: python3 settled_mass_route4_graph.py [RAWDIR [OUT.svg]]     (RAWDIR defaults to ~/Source/negative-mass-raw/I305-mac, OUT to settled_mass_route4.svg here; numbers from the same
functions as settled_mass_route4.py analyze)
The graph primitives are phase2/nonradial/data_graphs.py (the arriba graph-data skill's), unchanged."""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "nonradial"))
import data_graphs as dg                  # noqa: E402
raw = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Source/negative-mass-raw/I305-mac")
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "settled_mass_route4.svg")
sys.argv = ["x"]
import settled_mass_route4 as s4          # noqa: E402

s3, r = s4.s3, s4.r
fin = s4.H_FIN
D = {}
for nm, spec, tg in s4.POINTS4:
    e = r.load_point(os.path.join(raw, "runs"), spec, fin)
    e["name"] = nm
    D[(nm, fin)] = e
r.derive(D)
pts, c0 = s3.points15()
g1, rth, rdot, M = (np.array([p[k] for p in pts]) for k in ("g1", "rth", "rdot", "M"))
co = s3.fit(g1, rth, rdot, M, "P4")
ins = M - rth * (s3.M0 + s3.design(g1, rdot, "P4") @ co)
eloo, _ = s3.loo(g1, rth, rdot, M, "P4")
A15, G15, hat15 = s4.fitted_design(pts)
ptsN = s4.points_N(fin)
held = []
for p in ptsN:
    held.append(dict(nm=p["name"], g1=p["g1"], rd=p["rdot"], lev=s4.leverage(p["g1"], p["rdot"], A15, G15), res=p["M"] - s3.predict("P4", p["g1"], p["rth"], p["rdot"], fin, mutated=False), kind="N"))
for nm, spec, tg in s4.POINTS4:
    e = D[(nm, fin)]
    gn, rd = e["mem"]["g1"], e["sp"]
    held.append(dict(nm=nm, g1=gn, rd=rd, lev=s4.leverage(gn, rd, A15, G15), res=e["SS"]["plat1"] - s3.predict("P4", gn, e["mem"]["r_th"], rd, fin, mutated=False), kind="X"))
print("point  leverage  gamma_1    rdot_th    M - P4")
for p in held:
    print(f"{p['nm']:5s}  {p['lev']:.3f}    {p['g1']:.5f}  {p['rd']:+.3e}  {p['res']:+.2e}")

svg = dg.Svg(760)
X0, PW, PH = 86, 600, 190
ROUTE = {"B": ("square", "#9aa0a6"), "D": ("triangle", dg.SERIES[3]), "W": ("circle", dg.SERIES[0]), "C": ("diamond", dg.SERIES[1])}
y = 28
y = dg.caption(svg, 16, y, "Measured minus the frozen refit, against leverage in the fit's design",
               "h = 0.005; filled: fitted in-sample, open: fitted leave-one-out; grey: #295's held-out; green: this route; lines: bar 6e-5; vertical: fitted maximum")
sq = np.sqrt
xm0, ym = dg.axes(svg, X0, y, PW, PH, (0.0, 1.85), (-3.0e-4, 3.0e-4), "sqrt of the leverage L = a^T (A^T A)^-1 a", "M - predicted", xfmt=lambda v: f"{v:.1f}", yfmt=lambda v: f"{v:+.0e}")
xm = lambda L: xm0(sq(L))             # noqa: E731  (the axis is sqrt(L): L = 0.39, 1, 2, 3 sit at 0.63, 1.0, 1.41, 1.73)
for s_ in (+6e-5, -6e-5):
    svg.line(X0, ym(s_), X0 + PW, ym(s_), dg.SERIES[-1], 1)
svg.line(X0, ym(0.0), X0 + PW, ym(0.0), "#9aa0a6", 1)
svg.line(xm(float(hat15.max())), y, xm(float(hat15.max())), y + PH, "#c8ccd0", 1)
for p, a_, b_, lv in zip(pts, ins, eloo, hat15):
    sh, col = ROUTE[p["route"]]
    dg.marker(svg, sh, xm(lv), ym(a_), col, f"{p['name']} (fitted), in-sample M - P4 = {a_:+.2e}, leverage {lv:.3f}")
    svg.add(f'<circle cx="{xm(lv):.1f}" cy="{ym(b_):.1f}" r="5" fill="none" stroke="{col}" stroke-width="1.5"><title>{dg.esc(p["name"])} (fitted), leave-one-out M - P4 = {b_:+.2e}</title></circle>')
for p in held:
    col = "#9aa0a6" if p["kind"] == "N" else dg.SERIES[2]
    svg.add(f'<circle cx="{xm(p["lev"]):.1f}" cy="{ym(p["res"]):.1f}" r="8" fill="{col}" stroke="{dg.INK}" stroke-width="1.5"><title>{dg.esc(p["nm"])} (held out): leverage {p["lev"]:.3f}, M - P4 = {p["res"]:+.2e}</title></circle>')
    if p["kind"] == "X":
        svg.text(xm(p["lev"]) + 11, ym(p["res"]) + 4, p["nm"], size=11)
dg.legend(svg, X0 + PW, y - 10, [("B", "#9aa0a6"), ("delta scan", dg.SERIES[3]), ("route W", dg.SERIES[0]), ("route C", dg.SERIES[1]), ("new route", dg.SERIES[2])])

y1 = y + PH + 80
y = dg.caption(svg, 16, y1, "Where the members sit in the fit's design: rdot_th against gamma_1 at departure",
               "h = 0.005; lines: the fitted points' rdot_th range [-6.7e-4, +1.6e-3]; the new route lies below it at gamma_1 = 0.012, 0.048 and 0.075")
xm2, ym2 = dg.axes(svg, X0, y, PW, PH, (0.0, 0.085), (-3.2e-3, 2.2e-3), "gamma_1 at departure", "rdot_th at departure", xfmt=lambda v: f"{v:.2f}", yfmt=lambda v: f"{v * 1e3:+.1f}e-3")
for s_ in (-6.7e-4, 1.6e-3):
    svg.line(X0, ym2(s_), X0 + PW, ym2(s_), "#c8ccd0", 1)
for p in pts:
    sh, col = ROUTE[p["route"]]
    dg.marker(svg, sh, xm2(p["g1"]), ym2(p["rdot"]), col, f"{p['name']} (fitted): gamma_1 {p['g1']:.5f}, rdot_th {p['rdot']:+.3e}")
for p in held:
    col = "#9aa0a6" if p["kind"] == "N" else dg.SERIES[2]
    svg.add(f'<circle cx="{xm2(p["g1"]):.1f}" cy="{ym2(p["rd"]):.1f}" r="7" fill="{col}" stroke="{dg.INK}" stroke-width="1.5"><title>{dg.esc(p["nm"])}: gamma_1 {p["g1"]:.5f}, rdot_th {p["rd"]:+.3e}, leverage {p["lev"]:.3f}</title></circle>')
    svg.text(xm2(p["g1"]) + 10, ym2(p["rd"]) + 4, p["nm"], size=11)
dg.legend(svg, X0 + PW, y - 10, [("B", "#9aa0a6"), ("delta scan", dg.SERIES[3]), ("route W", dg.SERIES[0]), ("route C", dg.SERIES[1]), ("new route", dg.SERIES[2])])
try:
    commit = subprocess.check_output(["git", "-C", HERE, "rev-parse", "--short", "HEAD"], text=True).strip()
except (subprocess.SubprocessError, OSError):
    commit = "?"
height = dg.footnote(svg, [f"settled_mass_route4.py at negative-mass {commit}", "chain: #41 via #253", "member: #225 via #266", "speed: #267, #255", "refit: #295, frozen before the runs"], dg.host_specs(),
                     [f"the frozen refit's largest leave-one-out error is {np.max(np.abs(eloo)):.1e}; the new members' largest miss is {max(abs(p['res']) for p in held if p['kind'] == 'X'):.1e}"], y1 + PH + 100,
                     repo_url="https://github.com/nmoinvaz/negative-mass")
open(out, "w").write(svg.finish(height))
print("wrote", out)

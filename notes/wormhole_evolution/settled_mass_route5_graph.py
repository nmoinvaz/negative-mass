#!/usr/bin/env python3
"""Figure of #313: the frozen refit's errors against leverage in its design (the fifteen fitted points in-sample and leave-one-out; the nine held-out members: #295's N, #305's X and this issue's Y),
and the held-out misses against the throat speed at departure (h = 0.005).
usage: python3 settled_mass_route5_graph.py [RAWDIR [OUT.svg]]     (RAWDIR defaults to ~/Source/negative-mass-raw/I313-mac, OUT to settled_mass_route5.svg here; numbers from the same
functions as settled_mass_route5.py analyze)
The graph primitives are phase2/nonradial/data_graphs.py (the arriba graph-data skill's), unchanged."""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "nonradial"))
import data_graphs as dg                  # noqa: E402
raw = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Source/negative-mass-raw/I313-mac")
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "settled_mass_route5.svg")
sys.argv = ["x"]
import settled_mass_route5 as s5          # noqa: E402

s4, s3, r = s5.s4, s5.s3, s5.r
fin = s5.H_FIN
D = {}
for nm, spec, tg in s5.POINTS5:
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
held = []
for kind, pl in (("N", s4.points_N(fin)), ("X", s5.points_X(fin))):
    for p in pl:
        held.append(dict(nm=p["name"], g1=p["g1"], rd=p["rdot"], lev=s4.leverage(p["g1"], p["rdot"], A15, G15), res=p["M"] - s3.predict("P4", p["g1"], p["rth"], p["rdot"], fin, mutated=False), kind=kind))
for nm, spec, tg in s5.POINTS5:
    e = D[(nm, fin)]
    gn, rd = e["mem"]["g1"], e["sp"]
    held.append(dict(nm=nm, g1=gn, rd=rd, lev=s4.leverage(gn, rd, A15, G15), res=e["SS"]["plat1"] - s3.predict("P4", gn, e["mem"]["r_th"], rd, fin, mutated=False), kind="Y"))
print("point  leverage  gamma_1    rdot_th    M - P4")
for p in held:
    print(f"{p['nm']:5s}  {p['lev']:.3f}    {p['g1']:.5f}  {p['rd']:+.3e}  {p['res']:+.2e}")
kx = np.mean([p["res"] / p["rd"] for p in held if p["kind"] == "X"])
ky = np.mean([p["res"] / p["rd"] for p in held if p["kind"] == "Y"])
print(f"(M - P4)/rdot_th: X mean {kx:+.4f}, Y mean {ky:+.4f}")

svg = dg.Svg(760)
X0, PW, PH = 86, 600, 190
ROUTE = {"B": ("square", "#9aa0a6"), "D": ("triangle", dg.SERIES[3]), "W": ("circle", dg.SERIES[0]), "C": ("diamond", dg.SERIES[1])}
HELD = {"N": dg.SERIES[6], "X": dg.SERIES[4], "Y": dg.SERIES[2]}
LEG = [("B", "#9aa0a6"), ("delta scan", dg.SERIES[3]), ("route W", dg.SERIES[0]), ("route C", dg.SERIES[1])]
LEGH = [("#295 N", HELD["N"]), ("#305 X", HELD["X"]), ("this issue Y", HELD["Y"])]
y = 28
y = dg.caption(svg, 16, y, "Measured minus the frozen refit, against leverage in the fit's design",
               "h = 0.005; small markers: the fifteen fitted points (filled: in-sample, open: leave-one-out); large: held out; lines: bar 6e-5; vertical: fitted maximum")
sq = np.sqrt
xm0, ym = dg.axes(svg, X0, y, PW, PH, (0.0, 1.85), (-1.5e-4, 3.0e-4), "sqrt of the leverage L = a^T (A^T A)^-1 a", "M - predicted", xfmt=lambda v: f"{v:.1f}", yfmt=lambda v: f"{v:+.0e}")
xm = lambda L: xm0(sq(L))             # noqa: E731  (the axis is sqrt(L): L = 0.39, 1, 2, 3 sit at 0.63, 1.0, 1.41, 1.73)
for s_ in (+6e-5, -6e-5):
    svg.line(X0, ym(s_), X0 + PW, ym(s_), dg.SERIES[-1], 1)
svg.line(X0, ym(0.0), X0 + PW, ym(0.0), "#9aa0a6", 1)
svg.line(xm(float(hat15.max())), y, xm(float(hat15.max())), y + PH, "#c8ccd0", 1)
svg.line(xm(1.0), y, xm(1.0), y + PH, "#c8ccd0", 1)
for p, a_, b_, lv in zip(pts, ins, eloo, hat15):
    sh, col = ROUTE[p["route"]]
    dg.marker(svg, sh, xm(lv), ym(a_), col, f"{p['name']} (fitted), in-sample M - P4 = {a_:+.2e}, leverage {lv:.3f}")
    svg.add(f'<circle cx="{xm(lv):.1f}" cy="{ym(b_):.1f}" r="5" fill="none" stroke="{col}" stroke-width="1.5"><title>{dg.esc(p["name"])} (fitted), leave-one-out M - P4 = {b_:+.2e}</title></circle>')
for p in held:
    col = HELD[p["kind"]]
    svg.add(f'<circle cx="{xm(p["lev"]):.1f}" cy="{ym(p["res"]):.1f}" r="8" fill="{col}" stroke="{dg.INK}" stroke-width="1.5"><title>{dg.esc(p["nm"])} (held out): leverage {p["lev"]:.3f}, M - P4 = {p["res"]:+.2e}</title></circle>')
    if p["kind"] != "N":
        svg.text(xm(p["lev"]) + 11, ym(p["res"]) + (4 if p["kind"] == "X" else -10), p["nm"], size=11)
dg.legend(svg, X0 + PW, y - 10, LEG + LEGH)
svg.text(xm(1.0) + 4, y + 12, "L = 1", size=10)

y1 = y + PH + 80
y = dg.caption(svg, 16, y1, "The held-out misses against the throat speed at departure",
               "h = 0.005; lines: the post hoc description M - P4 = k rdot_th, k the mean of each route's three (X, Y); a description of single routes, not a law")
xm2, ym2 = dg.axes(svg, X0, y, PW, PH, (-3.2e-3, 1.8e-3), (-1.5e-4, 3.0e-4), "rdot_th at departure", "M - predicted", xfmt=lambda v: f"{v * 1e3:+.1f}e-3", yfmt=lambda v: f"{v:+.0e}")
for s_ in (+6e-5, -6e-5):
    svg.line(X0, ym2(s_), X0 + PW, ym2(s_), dg.SERIES[-1], 1)
svg.line(X0, ym2(0.0), X0 + PW, ym2(0.0), "#9aa0a6", 1)
for s_ in (-6.7e-4, 1.6e-3):
    svg.line(xm2(s_), y, xm2(s_), y + PH, "#c8ccd0", 1)
for k_, col in ((kx, HELD["X"]), (ky, HELD["Y"])):
    svg.line(xm2(-3.0e-3), ym2(k_ * -3.0e-3), xm2(0.0), ym2(0.0), col, 1)
for p, a_ in zip(pts, ins):
    sh, col = ROUTE[p["route"]]
    dg.marker(svg, sh, xm2(p["rdot"]), ym2(a_), col, f"{p['name']} (fitted), in-sample M - P4 = {a_:+.2e}, rdot_th {p['rdot']:+.3e}")
for p in held:
    col = HELD[p["kind"]]
    svg.add(f'<circle cx="{xm2(p["rd"]):.1f}" cy="{ym2(p["res"]):.1f}" r="7" fill="{col}" stroke="{dg.INK}" stroke-width="1.5"><title>{dg.esc(p["nm"])}: rdot_th {p["rd"]:+.3e}, M - P4 = {p["res"]:+.2e}, leverage {p["lev"]:.3f}</title></circle>')
    if p["kind"] != "N":
        svg.text(xm2(p["rd"]) + (10 if p["kind"] == "X" else -4), ym2(p["res"]) + (4 if p["kind"] == "X" else -11), p["nm"], size=11, anchor="start" if p["kind"] == "X" else "end")
dg.legend(svg, X0 + PW, y - 10, LEG + LEGH)
try:
    commit = subprocess.check_output(["git", "-C", HERE, "rev-parse", "--short", "HEAD"], text=True).strip()
except (subprocess.SubprocessError, OSError):
    commit = "?"
warn = []
for nm, spec, tg in s5.POINTS5:
    e1 = r.load_point(os.path.join(raw, "runs"), spec, 0.01)
    e1["name"], e1["h"] = nm, 0.01
    r.derive({0: e1})
    miss_ = [p["res"] for p in held if p["nm"] == nm][0]
    hd = abs(D[(nm, fin)]["SS"]["plat1"] - e1["SS"]["plat1"]) if "SS" in e1 else float("nan")
    margin = 6e-5 - abs(miss_)
    if margin < hd:
        warn.append(f"{nm}: the margin to the bar ({margin:.1e}) is smaller than M's change between h = 0.01 and 0.005 ({hd:.1e})")
height = dg.footnote(svg, [f"settled_mass_route5.py at negative-mass {commit}", "chain: #41 via #253", "member: #225 via #266", "speed: #267, #255", "refit: #295, frozen before the runs"], dg.host_specs(),
                     warn, y1 + PH + 100, repo_url="https://github.com/nmoinvaz/negative-mass")
open(out, "w").write(svg.finish(height))
print("wrote", out)

#!/usr/bin/env python3
"""Figure of #283: the settled mass of family G on two routes to its member, against family B, and the residuals of P2 and P3 (h = 0.005).
usage: python3 settled_mass_route_graph.py [RAWDIR [OUT.svg]]     (RAWDIR defaults to ~/Source/negative-mass-raw/I283-mac, OUT to settled_mass_route.svg here; numbers from the same functions as settled_mass_route.py analyze)
The graph primitives are phase2/nonradial/data_graphs.py (the arriba graph-data skill's), unchanged."""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "nonradial"))
import data_graphs as dg                  # noqa: E402
raw = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Source/negative-mass-raw/I283-mac")
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "settled_mass_route.svg")
sys.argv = ["x"]
import settled_mass_route as sr           # noqa: E402

fin = sr.H_FIN
root = os.path.join(raw, "runs")
D = {}
for nm, spec, tg in sr.POINTS:
    e = sr.load_point(root, spec, fin)
    e["name"] = nm
    D[(nm, fin)] = e
sr.derive(D)
curves = sr.fam_b(fin)
R = curves["R"]
c0 = sr.CRITR["M0"][fin]
fb = [(B, R[B]["mem"]["g1"], R[B]["mem"]["r_th"], R[B]["speed"], R[B]["plat"]) for B in (0.025, 0.05, 0.1, 0.125)]
pts = []
for nm, spec, tg in sr.POINTS:
    e = D[(nm, fin)]
    M, g1, rth, rd = e["SS"]["plat1"], e["mem"]["g1"], e["mem"]["r_th"], e["sp"]
    pts.append(dict(nm=nm, g1=g1, r=rth, rd=rd, M=M, y=M / rth - c0, r2=M - sr.predict("P2", g1, rth, rd, fin), r3=M - sr.predict("P3", g1, rth, rd, fin)))
print("point  gamma_1   r_th      rdot_th    M          M - P2     M - P3")
for p in pts:
    print(f"{p['nm']:5s}  {p['g1']:.5f}  {p['r']:.5f}  {p['rd']:+.2e}  {p['M']:.6f}  {p['r2']:+.2e}  {p['r3']:+.2e}")

svg = dg.Svg(760)
X0, PW, PH = 86, 600, 170
SH = {"W": ("circle", dg.SERIES[0]), "C": ("diamond", dg.SERIES[1])}
y = 28
y = dg.caption(svg, 16, y, "M / r_th - M_0 minus P2's curve, against the member's asymmetry gamma_1",
               "family G, D = 4, delta = 1e-8, h = 0.005; zero is P2 (0.2027 g + 0.8534 g^2); grey: family B (#266), circles: route W, diamonds: route C")
xm, ym = dg.axes(svg, X0, y, PW, PH, (0.0, 0.08), (-9e-4, 3e-4), "gamma_1 at departure", "M / r_th - M_0 - P2", xfmt=lambda v: f"{v:.2f}", yfmt=lambda v: f"{v:+.1e}")
a2, b2 = sr.CRITR["P2"][fin]
svg.line(X0, ym(0.0), X0 + PW, ym(0.0), dg.SERIES[-1], 1)
for B, g1, rth, rd, M in fb:
    dg.marker(svg, "square", xm(g1), ym(M / rth - c0 - a2 * g1 - b2 * g1 * g1), "#9aa0a6", f"family B, B = {B:g}: gamma_1 {g1:.5f}, excess over P2 {M / rth - c0 - a2 * g1 - b2 * g1 * g1:+.2e}")
for p in pts:
    sh, col = SH[p["nm"][0]]
    ex = p["y"] - a2 * p["g1"] - b2 * p["g1"] ** 2
    dg.marker(svg, sh, xm(p["g1"]), ym(ex), col, f"{p['nm']}: gamma_1 {p['g1']:.5f}, r_th {p['r']:.5f}, excess over P2 {ex:+.2e}")
dg.legend(svg, X0 + PW, y - 10, [("family B", "#9aa0a6"), ("route W", dg.SERIES[0]), ("route C", dg.SERIES[1])])
y1 = y + PH + 80
y = dg.caption(svg, 16, y1, "The throat's speed at departure against gamma_1", "the third parameter of P3: dr_th/dtau of the pair's mean at t_dep; grey line: family B")
xm2, ym2 = dg.axes(svg, X0, y, PW, PH, (0.0, 0.08), (-1.2e-3, 1.8e-3), "gamma_1 at departure", "dr_th/dtau at departure", xfmt=lambda v: f"{v:.2f}", yfmt=lambda v: f"{v:+.1e}")
x_, y_ = curves["g1"], curves["rdot"]
svg.add('<polyline fill="none" stroke="#9aa0a6" stroke-width="1.5" points="' + " ".join(f"{xm2(a):.1f},{ym2(b_):.1f}" for a, b_ in zip(x_, y_)) + '"/>')
for B, g1, rth, rd, M in fb:
    dg.marker(svg, "square", xm2(g1), ym2(rd), "#9aa0a6", f"family B, B = {B:g}: gamma_1 {g1:.5f}, rdot_th {rd:+.2e}")
for p in pts:
    sh, col = SH[p["nm"][0]]
    dg.marker(svg, sh, xm2(p["g1"]), ym2(p["rd"]), col, f"{p['nm']}: gamma_1 {p['g1']:.5f}, rdot_th {p['rd']:+.2e}")
y2 = y1 + PH + 112
y = dg.caption(svg, 16, y2, "Measured minus predicted, P2 and P3", "lines: the declared bar +-6e-5; both predictions fixed by the issue before the runs")
lo, hi = -9e-4, 2e-4
xm3, ym3 = dg.axes(svg, X0, y, PW, PH, (-0.5, 5.5), (lo, hi), "point (W: route W, C: route C; 1, 2, 3: gamma_1 near 0.012, 0.048, 0.075)", "M - predicted (a)", xfmt=lambda v: "", yfmt=lambda v: f"{v:+.1e}")
for k in range(6):
    svg.text(xm3(k), y + PH + 16, sr.POINTS[k][0], size=12, anchor="middle")
for s_ in (+6e-5, -6e-5):
    svg.line(X0, ym3(s_), X0 + PW, ym3(s_), dg.SERIES[-1], 1)
for k, p in enumerate(pts):
    dg.marker(svg, "triangle", xm3(k), ym3(p["r2"]), dg.SERIES[3], f"{p['nm']}: M - P2 = {p['r2']:+.2e}")
    dg.marker(svg, "circle", xm3(k), ym3(p["r3"]), dg.SERIES[2], f"{p['nm']}: M - P3 = {p['r3']:+.2e}")
dg.legend(svg, X0 + PW, y - 10, [("P2 (triangles)", dg.SERIES[3]), ("P3 (circles)", dg.SERIES[2])])
try:
    commit = subprocess.check_output(["git", "-C", HERE, "rev-parse", "--short", "HEAD"], text=True).strip()
except (subprocess.SubprocessError, OSError):
    commit = "?"
height = dg.footnote(svg, [f"settled_mass_route.py at negative-mass {commit}", "chain: #41 via #253", "member: #225 via #266", "speed: #267, #255", "P2: #266, P3: #267"], dg.host_specs(),
                     [f"P2 and P3 residuals against the bar 6e-5; W3's P3 residual is {pts[2]['r3']:+.1e}"], y2 + PH + 100, repo_url="https://github.com/nmoinvaz/negative-mass")
open(out, "w").write(svg.finish(height))
print("wrote", out)

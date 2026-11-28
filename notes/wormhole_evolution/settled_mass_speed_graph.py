#!/usr/bin/env python3
"""Figure of #314: the frozen refit's misses at X1 to X3 when the throat speed and the member are read by other estimators (refit over the fifteen with each reading), and how far the speed estimators
move the speed at X1 to X3 against the shift that would close the miss (h = 0.005).
usage: python3 settled_mass_speed_graph.py [OUT.svg]     (numbers from the same functions as settled_mass_speed.py analyze; stored runs read from ~/Source/negative-mass-raw)
The graph primitives are phase2/nonradial/data_graphs.py (the arriba graph-data skill's), unchanged."""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "nonradial"))
import data_graphs as dg                  # noqa: E402
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "settled_mass_speed.svg")
sys.argv = ["x"]
import settled_mass_speed as S            # noqa: E402

P = S.load_points("all")
rd = [S.readings_of(p) for p in P]
idxN = [k for k, p in enumerate(P) if p["grp"] == "N"]
idxX = [k for k, p in enumerate(P) if p["grp"] == "X"]
idxY = [k for k, p in enumerate(P) if p["grp"] == "Y"]
n15 = sum(1 for p in P if p["grp"] not in ("N", "X"))
keys = S.READINGS
res = {k: S.refit_residuals(P, rd, k)[1] for k in keys}
co0, res0 = S.refit_residuals(P, rd, "REF")
need = [res0[k] / (rd[k]["REF"][1] * S.P4[2]) for k in idxX]
print("reading   X1        X2        X3   | Y1        Y2        Y3   | N max")
for k in keys:
    print(f"{k:5s}  " + "  ".join(f"{res[k][i]:+.2e}" for i in idxX) + "  | " + "  ".join(f"{res[k][i]:+.2e}" for i in idxY) + f"  | {max(abs(res[k][i]) for i in idxN):.1e}")

svg = dg.Svg(760)
X0, PW, PH = 86, 600, 200
y = 28
y = dg.caption(svg, 16, y, "The miss at X1 to X3 and Y1 to Y3 under each reading of the speed and the member",
               "h = 0.005; M - predicted after refitting the fifteen with the reading; lines: bar 6e-5; grey: #295's N1 to N3 (largest of the three); amber: #313's Y1 to Y3")
xm, ym = dg.axes(svg, X0, y, PW, PH, (-0.5, len(keys) - 0.5), (-1e-4, 3.0e-4), "reading", "M - predicted", xfmt=lambda v: "", yfmt=lambda v: f"{v:+.0e}")
for j, k in enumerate(keys):
    svg.text(xm(j), y + PH + 16, k, size=12, anchor="middle")
for s_ in (+6e-5, -6e-5):
    svg.line(X0, ym(s_), X0 + PW, ym(s_), dg.SERIES[-1], 1)
svg.line(X0, ym(0.0), X0 + PW, ym(0.0), "#9aa0a6", 1)
SH = ("circle", "square", "triangle")
for j, k in enumerate(keys):
    nmax = max(res[k][i] for i in idxN) if max(abs(res[k][i]) for i in idxN) == max(res[k][i] for i in idxN) else min(res[k][i] for i in idxN)
    svg.add(f'<circle cx="{xm(j) - 32:.1f}" cy="{ym(nmax):.1f}" r="4" fill="#9aa0a6"><title>{dg.esc(k)}: N1 to N3, largest |M - predicted| {nmax:+.2e}</title></circle>')
    for q, i in enumerate(idxX):
        dg.marker(svg, SH[q], xm(j) + (q - 2) * 10 - 8, ym(res[k][i]), dg.SERIES[q], f"{P[i]['name']} under {k}: M - predicted = {res[k][i]:+.2e}")
    for q, i in enumerate(idxY):
        dg.marker(svg, SH[q], xm(j) + q * 10 + 14, ym(res[k][i]), dg.SERIES[3], f"{P[i]['name']} under {k}: M - predicted = {res[k][i]:+.2e}")
dg.legend(svg, X0 + PW, y - 10, [("X1", dg.SERIES[0]), ("X2", dg.SERIES[1]), ("X3", dg.SERIES[2]), ("Y1 to Y3", dg.SERIES[3])])

y1 = y + PH + 80
y = dg.caption(svg, 16, y1, "How far the speed estimators move rdot_th at X1 to X3, against the shift that closes the miss",
               "h = 0.005; shift = estimator minus #267's speed_at; marks: the needed shift +(M - P4)/(r_th c) of the reference refit; S1r reads dr/dt, a different variable, and is not shown")
xm2, ym2 = dg.axes(svg, X0, y, PW, PH, (-0.5, 2.5), (-6e-5, 3.2e-4), "point", "change of rdot_th", xfmt=lambda v: ("X1", "X2", "X3")[int(round(v))] if abs(v - round(v)) < 1e-9 and 0 <= round(v) < 3 else "",
                   yfmt=lambda v: f"{v * 1e4:+.1f}e-4")
svg.line(X0, ym2(0.0), X0 + PW, ym2(0.0), "#9aa0a6", 1)
for q, i in enumerate(idxX):
    svg.add(f'<line x1="{xm2(q) - 45:.1f}" y1="{ym2(need[q]):.1f}" x2="{xm2(q) + 45:.1f}" y2="{ym2(need[q]):.1f}" stroke="{dg.INK}" stroke-width="2.5"><title>{P[i]["name"]}: needed shift {need[q]:+.2e}</title></line>')
    for j, k in enumerate(("S1", "S2", "S3")):
        d = rd[i][k][2] - rd[i]["REF"][2]
        dg.marker(svg, ("circle", "square", "triangle")[j], xm2(q) + (j - 1) * 20, ym2(d), dg.SERIES[j], f"{P[i]['name']}, {k}: change of rdot_th {d:+.2e} (needed {need[q]:+.2e})")
dg.legend(svg, X0 + PW, y - 10, [("S1 (slope in t)", dg.SERIES[0]), ("S2 (+-1.0 in tau)", dg.SERIES[1]), ("S3 (quintic)", dg.SERIES[2])])
try:
    commit = subprocess.check_output(["git", "-C", HERE, "rev-parse", "--short", "HEAD"], text=True).strip()
except (subprocess.SubprocessError, OSError):
    commit = "?"
height = dg.footnote(svg, [f"settled_mass_speed.py at negative-mass {commit}", "chain: #41 via #253", "member: #225 via #266", "speed: #267, #255", "refit: #295, frozen before the runs"], dg.host_specs(),
                     [f"the largest change of rdot_th at X1 to X3 by S1, S2, S3 is {max(abs(rd[i][k][2] - rd[i]['REF'][2]) for i in idxX for k in ('S1', 'S2', 'S3')):.1e}; the needed shift is {min(abs(x) for x in need):.1e} to {max(abs(x) for x in need):.1e}"],
                     y1 + PH + 100, repo_url="https://github.com/nmoinvaz/negative-mass")
open(out, "w").write(svg.finish(height))
print("wrote", out)

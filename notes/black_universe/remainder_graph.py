#!/usr/bin/env python3
"""The late growth at fixed rho of smooth c = -1 data against #262's profile (#322).

    python3 remainder_graph.py SERIES_DIR      (SERIES_DIR = negative-mass-raw/I322-mac/runs/r1/series, from `remainder.py dump`; writes remainder.svg beside this file)

(a) Y_rem / (c2 V R) against V = e^{kappa v} for class G at three fixed rho, with the prediction -J/I0 (dashed).
(b) the growth coefficient a/c2 of every datum at the three rho against the curve -J(rho)/I0 of #262.
Prints the plotted numbers as a table.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "nonradial"))
from data_graphs import INK, SERIES, Svg, axes, callout, caption, footnote, legend, marker  # noqa: E402

_ARGV = sys.argv
sys.argv = sys.argv[:1]
import remainder as rm  # noqa: E402
sys.argv = _ARGV

series_dir = sys.argv[1]
names = [n for n, _ in rm.data_items()]
data = {n: np.load(os.path.join(series_dir, f"series_{n}.npz")) for n in names}
cf = rm.closed_form()
pis = rm.PIS
rhos = [rm.rho_R(p)[0] for p in pis]
Qp = [rm.Q_pred(r) for r in rhos]

width = 780
svg = Svg(width)
svg.text(16, 26, "Late growth at fixed ρ: Y_rem does not decay, it grows at κ with #262's profile", size=14, fill=INK, weight="bold")
px, pw, ph = 96, 600, 230
py = caption(svg, 16, 54, "(a) Y_rem / (c₂ V R) for class G (Richardson, h = 0.000625, 0.0003125)",
             "dashed: −J/I₀ = +0.2961, −0.5737, −0.7327 at ρ = −0.311, 0.682, 1; the ratio would fall to 0 if Y_rem decayed, and sit at 0 if the growth were c₂VR.")
d = data["G"]
c2 = float(d["c2"][2] * 4 - d["c2"][1]) / 3.0
xm, ym = axes(svg, px, py, pw, ph, (2.0, 24.0, "log"), (-0.9, 0.5), "V = e^(κv) (log scale)", "Y_rem / (c₂ V R)", xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:g}")
items = []
for k in range(3):
    V, w = d[f"V_{k}"], d[f"w_{k}"]
    m = V >= 2.0
    pts = " ".join(f"{xm(v):.1f},{ym(y / (c2 * v)):.1f}" for v, y in zip(V[m][::max(1, m.sum() // 200)], w[m][::max(1, m.sum() // 200)]))
    svg.add(f'<polyline points="{pts}" fill="none" stroke="{SERIES[k]}" stroke-width="2"/>')
    svg.add(f'<line x1="{px}" y1="{ym(Qp[k]):.1f}" x2="{px + pw}" y2="{ym(Qp[k]):.1f}" stroke="{SERIES[k]}" stroke-width="1" stroke-dasharray="5 4"/>')
    items.append((f"Π = {pis[k]:g}", SERIES[k]))
svg.add(f'<line x1="{px}" y1="{ym(0):.1f}" x2="{px + pw}" y2="{ym(0):.1f}" stroke="{INK}" stroke-width="1" stroke-dasharray="2 3"/>')
svg.text(px + pw - 6, ym(0) - 5, "0: late growth would be c₂VR", size=11, fill=INK, anchor="end")
legend(svg, px + pw, py - 12, items)

py2 = caption(svg, 16, py + ph + 62, "(b) Growth coefficient a/c₂ of Y_rem/R = aV + b + d/V, all data with c₂ ≠ 0, against −J(ρ)/I₀",
              "fit on V ∈ [8, 24]; one marker per datum and ρ; the curve is #262's Y_∞/Y_ref − 1 (ρ < 0 is inside the horizon).")
ph2 = 230
xm2, ym2 = axes(svg, px, py2, pw, ph2, (-0.4, 1.1), (-0.85, 0.4), "ρ", "a / c₂", xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:g}")
rr = np.linspace(-0.4, 1.1, 120)
svg.add('<polyline points="%s" fill="none" stroke="%s" stroke-width="2"/>' % (" ".join(f"{xm2(r):.1f},{ym2(rm.Q_pred(float(r))):.1f}" for r in rr), SERIES[3]))
worst = 0.0
print("| datum | rho | a/c2 | -J/I0 |")
print("|---|---|---|---|")
for n in names:
    d = data[n]
    c2 = (4.0 * float(d["c2"][2]) - float(d["c2"][1])) / 3.0
    if abs(c2) < rm.C2MIN:
        continue
    for k in range(3):
        V, w = d[f"V_{k}"], d[f"w_{k}"]
        a = rm.fit_poly3(V, w, rm.VA, rm.VB)[0]
        marker(svg, "circle", xm2(rhos[k]), ym2(a / c2), SERIES[0], f"{n}: ρ = {rhos[k]:.4f}, a/c₂ = {a / c2:+.5f}, −J/I₀ = {Qp[k]:+.5f}")
        worst = max(worst, abs(a / c2 - Qp[k]))
        print(f"| {n} | {rhos[k]:.4f} | {a / c2:+.5f} | {Qp[k]:+.5f} |")
legend(svg, px + pw, py2 - 12, [("data (11 at three ρ)", SERIES[0]), ("−J(ρ)/I₀", SERIES[3])])
callout(svg, xm2(rhos[1]), ym2(Qp[1]), f"all data within {worst:.1e} of the curve", dx=40, dy=-50)
h = footnote(svg, ["remainder.py (#322, code 50e1d53)", "raw I322-mac/runs/r1/series", "I₀, Y_∞: #262"], "", [], py2 + ph2 + 44, "https://github.com/nmoinvaz/negative-mass")
path = os.path.join(HERE, "remainder.svg")
open(path, "w").write(svg.finish(h))
print(f"wrote {path}")

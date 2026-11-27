#!/usr/bin/env python3
"""The growth rate of the single radial mode of the black universes against the family parameter c, and its ratio to the
surface gravity of the horizon.

    python3 rate_graph.py        (reads the table in stability_output.txt, writes growth_rate.svg beside this file)

b = 1 units. c = 0 is the Ellis wormhole (W = 1.18197, planted); c = -1 has no growing mode.
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "nonradial"))
from data_graphs import INK, SERIES, Svg, axes, callout, caption, footnote, legend, marker  # noqa: E402

text = open(os.path.join(HERE, "stability_output.txt")).read()
rows = []
for m in re.finditer(r"^\s*(-?[0-9.]+)\s+(-?[0-9.]+)\s+([0-9.]+)\s+\|\s+(\d)\s+(\d)\s+\|\s+(none|[0-9.]+)", text, re.M):
    c, rh, kappa, n1, n2, W = m.groups()
    rows.append((float(c), float(kappa), None if W == "none" else float(W)))
assert len(rows) == 20, len(rows)
ellis = float(re.search(r"c = 0\) R2: growing modes 1 \(declared 1\); W\^2 a\^2 = ([0-9.]+)", text).group(1)) ** 0.5

width = 780
svg = Svg(width)
svg.text(16, 26, "The one growing radial mode of the black universes, W against c (b = 1)", size=14, fill=INK, weight="bold")
px, pw, ph = 96, 600, 230
py = caption(svg, 16, 54, "(a) Growth rate W and the horizon's surface gravity κ",
             "units 1/b; the throat radius b = 1 is the minimum of R; c = −1 puts the horizon at the throat and has no growing mode.")
xm, ym = axes(svg, px, py, pw, ph, (0.008, 12, "log"), (0, 1.4), "−c (log scale)", "rate (1/b)", xfmt=lambda v: f"{v:g}",
              yfmt=lambda v: f"{v:g}")
pts = [(xm(-c), ym(W)) for c, k, W in rows if W is not None]
svg.add('<polyline points="%s" fill="none" stroke="%s" stroke-width="2"/>' % (" ".join(f"{x:.1f},{y:.1f}" for x, y in pts), SERIES[0]))
svg.add('<polyline points="%s" fill="none" stroke="%s" stroke-width="2" stroke-dasharray="5 4"/>' % (
    " ".join(f"{xm(-c):.1f},{ym(k):.1f}" for c, k, W in rows), SERIES[1]))
for c, k, W in rows:
    if W is not None:
        marker(svg, "circle", xm(-c), ym(W), SERIES[0], f"c = {c:g}: W = {W:.10g}, κ = {k:.6g}")
svg.add(f'<line x1="{px}" y1="{ym(ellis):.1f}" x2="{px + pw}" y2="{ym(ellis):.1f}" stroke="{INK}" stroke-width="1" stroke-dasharray="2 3"/>')
svg.text(px + 6, ym(ellis) - 5, f"Ellis wormhole (c = 0, planted): W = {ellis:.5f}", size=11, fill=INK)
x1 = xm(1.0)
svg.add(f'<circle cx="{x1:.1f}" cy="{ym(0):.1f}" r="6" fill="none" stroke="{SERIES[0]}" stroke-width="2"/>')
svg.text(x1, ym(0) - 14, "c = −1: no mode", size=11, fill=INK, anchor="middle")
legend(svg, px + pw, py - 12, [("W", SERIES[0]), ("κ = A′(ρ_h)/2", SERIES[1])])

py2 = caption(svg, 16, py + ph + 62, "(b) W/κ",
              "W/κ → 1 on both sides of c = −1, but the mode is absent at c = −1: it is bound in a Pöschl–Teller well at the horizon (note).")
ph2 = 220
xm2, ym2 = axes(svg, px, py2, pw, ph2, (0.008, 12, "log"), (0, 2.8), "−c (log scale)", "W/κ", xfmt=lambda v: f"{v:g}",
                yfmt=lambda v: f"{v:g}")
svg.add(f'<line x1="{px}" y1="{ym2(1):.1f}" x2="{px + pw}" y2="{ym2(1):.1f}" stroke="{INK}" stroke-width="1" stroke-dasharray="2 3"/>')
pts2 = [(xm2(-c), ym2(min(W / k, 2.78))) for c, k, W in rows if W is not None]
svg.add('<polyline points="%s" fill="none" stroke="%s" stroke-width="2"/>' % (" ".join(f"{x:.1f},{y:.1f}" for x, y in pts2), SERIES[2]))
for c, k, W in rows:
    if W is not None:
        marker(svg, "circle", xm2(-c), ym2(min(W / k, 2.78)), SERIES[2], f"c = {c:g}: W/κ = {W / k:.6f}")
svg.text(xm2(0.0108) + 8, ym2(2.2) + 4, "c = −0.01: W/κ = %.2f (off scale)" % next(W / k for c, k, W in rows if c == -0.01), size=11, fill=INK)
h = footnote(svg, ["stability_output.txt (stability.py, #212)", "c = −0.1, …, −10 and the Ellis plant: raw I212-rpi5"], "", [],
             py2 + ph2 + 44, "https://github.com/nmoinvaz/negative-mass")
path = os.path.join(HERE, "growth_rate.svg")
open(path, "w").write(svg.finish(h))
print(f"wrote {path}")

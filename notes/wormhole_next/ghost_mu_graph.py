#!/usr/bin/env python3
"""mu(D), the least-damped real ghost-sector rate of Choptuik's solution, from ghost_mu.jsonl and ghost_mu2.jsonl.

    python3 ghost_mu_graph.py         (writes ghost_mu.svg beside this file, prints the numbers)

Selection as mu_inf.py's load(): mode-1, unplanted records, one per D (the highest M), whose exact zero at mu = 0 stands
out (pivot(0)/pivot(-h) < 1e-4), here from D = 4. The records at D = 72, 73.27 (marginal) and D >= 74 (singular) fail
that test and are drawn only as the shaded wall. Standard library only.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "nonradial"))
from data_graphs import (GRID, INK, INK_SOFT, SERIES, Svg, axes, better_hint, callout, caption,  # noqa: E402
                         footnote, legend, marker)

pts, rejected = {}, set()
for name in ("ghost_mu.jsonl", "ghost_mu2.jsonl"):
    for line in open(os.path.join(HERE, name)):
        r = json.loads(line)
        if r["mode"] != 1 or r.get("plant") is not None:
            continue
        if r["pivot0"] / r["pivot_h"] >= 1e-4:
            rejected.add(round(r["D"], 4))
            continue
        z = [x["mu"] for x in r["zeros"] if x["ok"]]
        if not z:
            continue
        D = round(r["D"], 4)
        if D not in pts or r["M"] > pts[D][0]:
            pts[D] = (r["M"], max(z))
DS = sorted(pts)
MU = [pts[D][1] for D in DS]
MGM = (-0.450, -0.446)                                     # [MGM26]'s window at D = 5, evolutions (wormhole_next.md §1.4)

print(f"{'D':>9}  {'M':>3}  {'mu':>15}")
for D in DS:
    print(f"{D:>9g}  {pts[D][0]:>3}  {pts[D][1]:>15.10f}")
print(f"rejected by the conditioning test (pivot ratio >= 1e-4): {sorted(rejected)}")

W = 780
svg = Svg(W)
px, pw, ph = 96, 590, 330
py = caption(svg, 16, 26, "The ghost sector of Choptuik's solution decays in every dimension computed",
             "μ = least-damped real rate of a ghost-field perturbation on Choptuik's critical solution, per echo. "
             "Negative = the ghost admixture dies out.")
xm, ym = axes(svg, px, py, pw, ph, (4, 100, "log"), (-1.1, -0.3), "spacetime dimension D",
              "μ (ghost-sector rate)", xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:.1f}")
x74 = xm(74)
svg.add(f'<rect x="{x74:.1f}" y="{py}" width="{px + pw - x74:.1f}" height="{ph}" fill="{GRID}" opacity="0.6"/>')
for i, t in enumerate(("D ≥ 74:", "matrix", "singular", "in double-", "double")):
    svg.text(x74 + 5, py + 16 + 13 * i, t, size=10)
y1 = ym(-1)
svg.add(f'<line x1="{px}" y1="{y1:.1f}" x2="{px + pw}" y2="{y1:.1f}" stroke="{INK_SOFT}" stroke-width="1" '
        f'stroke-dasharray="5 4"/>')
svg.text(px + 8, y1 + 15, "μ = −1: the conjectured large-D limit (not tested: ±0.02 cannot be resolved)", size=10)
xs = xm(5)
svg.line(xs + 12, ym(MGM[0]), xs + 12, ym(MGM[1]), SERIES[1], width=3)
svg.add(f'<g><rect x="{xs + 8:.1f}" y="{ym(MGM[1]):.1f}" width="8" height="{ym(MGM[0]) - ym(MGM[1]):.1f}" '
        f'fill="none" stroke="{SERIES[1]}" stroke-width="1.5"/><title>[MGM26], D = 5, evolutions: '
        f'{MGM[0]} to {MGM[1]}</title></g>')
path = " ".join(f"{xm(D):.1f},{ym(m):.1f}" for D, m in zip(DS, MU))
svg.add(f'<polyline points="{path}" fill="none" stroke="{SERIES[0]}" stroke-width="1.5" opacity="0.6"/>')
for D, m in zip(DS, MU):
    marker(svg, "circle", xm(D), ym(m), SERIES[0], f"D = {D:g}: μ = {m:.10f} (M = {pts[D][0]})")
callout(svg, xm(5), ym(pts[5.0][1]), "D = 5: −0.44304, outside evolutions' −0.450 to −0.446 (open)", dx=36, dy=-46)
callout(svg, xm(4), ym(pts[4.0][1]), "D = 4: −0.3651, within perturbation theory's −0.362 ± 0.012", dx=36, dy=-30)
better_hint(svg, px, py, pw, ph, where="right", label="more damped", reverse=True)
legend(svg, px + pw, py + ph + 54, [("test-field determinant on choptuik-large-d's solutions", SERIES[0]),
                                    ("evolutions at D = 5 [MGM26]", SERIES[1])])
h = footnote(svg, ["μ: ghost_mu.jsonl, ghost_mu2.jsonl (wormhole_next.md §§1.3, 6.3; #30)",
                   "solutions and determinant: choptuik-large-d, read-only"],
             "", [], py + ph + 64, "https://github.com/nmoinvaz/negative-mass")
out = os.path.join(HERE, "ghost_mu.svg")
open(out, "w").write(svg.finish(h))
print(f"wrote {out}")

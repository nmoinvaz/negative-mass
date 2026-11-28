#!/usr/bin/env python3
"""Figure of #253 (settled_mass_g.md): M_AH(x+ - x+_f) of family G against S, its plateau, and the plateau against delta.

    python3 settled_mass_g_graph.py [RAW [OUT.svg]]

RAW is the negative-mass-raw clone (default ~/Source/negative-mass-raw); it reads I253-rpi5/runs (G: h0.005_d1e-0{6,8},
scan_G_h0.01.json, scan_G_h0.005.json) and I41-rpi5/runs (S: h0.005_d1e-08, scan_h0.01.json). OUT defaults to
settled_mass_g.svg beside this file. It replaces I253-rpi5's matplotlib plot_settled_mass_g.py with the same data on
phase2/nonradial/data_graphs.py (the arriba graph-data skill's primitives, unchanged). Prints the plotted numbers.
"""
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "nonradial"))
import data_graphs as dg  # noqa: E402

RAW = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Source/negative-mass-raw")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "settled_mass_g.svg")
G_RUNS, S_RUNS = os.path.join(RAW, "I253-rpi5", "runs"), os.path.join(RAW, "I41-rpi5", "runs")


def series(d):
    z = np.load(os.path.join(d, "series.npz"))
    return z["xp"] - z["xp"][0], z["M"]


CURVES = [("G, delta = 1e-8", series(os.path.join(G_RUNS, "h0.005_d1e-08")), dg.SERIES[0], None),
          ("G, delta = 1e-6", series(os.path.join(G_RUNS, "h0.005_d1e-06")), dg.SERIES[0], "5 4"),
          ("S, delta = 1e-8", series(os.path.join(S_RUNS, "h0.005_d1e-08")), dg.SERIES[1], None)]
S_SCAN = sorted((r["delta"], r["M_end"]) for r in json.load(open(os.path.join(S_RUNS, "scan_h0.01.json")))
                if "M_end" in r and r["delta"] <= 1e-3)
G_SCAN = {h: sorted((r["delta"], r["plat1"]) for r in json.load(open(os.path.join(G_RUNS, f"scan_G_h{h}.json")))
                    if "plat1" in r and r["delta"] <= 1e-3) for h in ("0.01", "0.005")}

print("curve              points  x+ - x+_f end   M start    M end")
for lab, (x, M), _, _ in CURVES:
    print(f"{lab:18s} {len(x):6d}  {x[-1]:13.3f}  {M[0]:.6f}  {M[-1]:.6f}")
print("delta     S (h = 0.01, M_end)   G (h = 0.01, plateau)   G (h = 0.005, plateau)")
sd, g1, g2 = dict(S_SCAN), dict(G_SCAN["0.01"]), dict(G_SCAN["0.005"])
for d in sorted(set(sd) | set(g1) | set(g2), reverse=True):
    f = lambda t: f"{t[d]:.6f}" if d in t else "    -   "  # noqa: E731
    print(f"{d:<9.0e} {f(sd):>12s}          {f(g1):>12s}            {f(g2):>12s}")


def polyline(svg, xm, ym, x, y, color, dash, title, lo=None, hi=None, ylo=None, yhi=None):
    keep = np.ones(len(x), bool)
    if lo is not None:
        keep &= (x >= lo) & (x <= hi) & (y >= ylo) & (y <= yhi)
    idx = np.flatnonzero(keep)[::3]
    pts = " ".join(f"{xm(a):.1f},{ym(b):.1f}" for a, b in zip(x[idx], y[idx]))
    da = f' stroke-dasharray="{dash}"' if dash else ""
    svg.add(f'<polyline fill="none" stroke="{color}" stroke-width="1.8"{da} points="{pts}"><title>{dg.esc(title)}</title></polyline>')


svg = dg.Svg(760)
X0, PW, PH = 86, 600, 180
y = 28
y = dg.caption(svg, 16, y, "Settled black-hole mass of the one-sided family G against the mirrored family S",
               "M_AH = r_AH/2 on the apparent horizon from formation to the chain's limit e^-f = 1.57e6, D = 4, h = 0.005; a = 1")
xm, ym = dg.axes(svg, X0, y, PW, PH, (0, 11), (0.20, 0.52), "x+ - x+ at formation (a)", "M_AH (a)",
                 xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:.2f}")
for lab, (x, M), col, dash in CURVES:
    polyline(svg, xm, ym, x, M, col, dash, f"{lab}: M from {M[0]:.4f} to {M[-1]:.6f}")
dg.legend(svg, X0 + PW, y - 10, [("family G (dashed: delta = 1e-6)", dg.SERIES[0]), ("family S", dg.SERIES[1])])
dg.callout(svg, xm(CURVES[0][1][0][0]), ym(CURVES[0][1][1][0]), f"formation: {CURVES[0][1][1][0]:.3f}", dx=40, dy=10)

y1 = y + PH + 70
y = dg.caption(svg, 16, y1, "The plateau (same runs)", "G settles 0.0125 a above S, and G's value depends on delta")
xm2, ym2 = dg.axes(svg, X0, y, PW, PH, (4, 11), (0.2220, 0.2370), "x+ - x+ at formation (a)", "M_AH (a)",
                   xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:.4f}")
for lab, (x, M), col, dash in CURVES:
    polyline(svg, xm2, ym2, x, M, col, dash, f"{lab}: plateau end {M[-1]:.6f}", 4, 11, 0.2220, 0.2370)
for lab, (x, M), col, dash in (CURVES[0], CURVES[2]):
    dg.callout(svg, xm2(x[-1]), ym2(M[-1]), f"{lab.split(',')[0]}: {M[-1]:.5f}", dx=-110, dy=-22 if lab.startswith("G") else -26)

y2 = y1 + PH + 100
y = dg.caption(svg, 16, y2, "The settled mass against the distance from threshold",
               "S: final M_AH of #41's scan (h = 0.01); G: the plateau mean over the last 1.0 of x+ (#253's scan); delta = p*_h - A, closer to threshold to the right")
e_of = lambda d: -np.log10(d)  # noqa: E731
xm3, ym3 = dg.axes(svg, X0, y, PW, PH, (2.5, 12.5), (0.2220, 0.2390), "distance from threshold delta", "settled mass M_AH (a)",
                   xfmt=lambda v: f"1e-{v:g}", yfmt=lambda v: f"{v:.4f}")
for d, M in S_SCAN:
    dg.marker(svg, "circle", xm3(e_of(d)), ym3(M), dg.SERIES[1], f"S, delta = {d:.0e}: {M:.6f}")
pts = dg.dodge_markers({k: {d: (xm3(e_of(d)), ym3(M)) for d, M in G_SCAN[h]} for k, h in enumerate(("0.01", "0.005"))}, step=7.0)
for k, (h, shape) in enumerate((("0.01", "circle"), ("0.005", "square"))):
    for d, M in G_SCAN[h]:
        dg.marker(svg, shape, pts[k][d][0], pts[k][d][1], dg.SERIES[0], f"G, h = {h}, delta = {d:.0e}: {M:.6f}")
dg.legend(svg, X0 + PW, y - 10, [("G (circle h = 0.01, square h = 0.005)", dg.SERIES[0]), ("S (h = 0.01)", dg.SERIES[1])])
dmin, Mmin = min(G_SCAN["0.005"], key=lambda t: t[1])
dg.callout(svg, xm3(e_of(dmin)), ym3(Mmin), f"minimum near 1e-5: {Mmin:.5f}", dx=30, dy=-30)
d12, M12 = G_SCAN["0.005"][0]
dg.callout(svg, xm3(e_of(d12)), ym3(M12), f"1e-12: {M12:.5f}, still rising", dx=-60, dy=70)

try:
    commit = subprocess.check_output(["git", "-C", HERE, "rev-parse", "--short", "HEAD"], text=True).strip()
except (subprocess.SubprocessError, OSError):
    commit = "?"
height = dg.footnote(svg, [f"settled_mass_g_graph.py at negative-mass {commit}", "G: negative-mass-raw I253-rpi5 (a31d31d4)",
                           "S: I41-rpi5"], "data from runs on rpi5-ubuntu (Raspberry Pi 5)",
                     ["the third panel compares S's final M_AH (#41) with G's plateau mean (#253): two readings of the settled value"],
                     y2 + PH + 100, repo_url="https://github.com/nmoinvaz/negative-mass")
open(OUT, "w").write(svg.finish(height))
print("wrote", OUT)

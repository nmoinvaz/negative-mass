#!/usr/bin/env python3
"""Paper I in one figure: the threshold kappa_c(R) of Theorem 2, and the growth rate of the point mass against finite stars.

    python3 thresholds_graph.py        (writes thresholds.svg beside this file, prints the plotted numbers)

Reads kappa_c.csv and growth_rates.csv (both written by truncated.py) and ../../paper/figures/thresholds.csv (r_s, r_*,
r_neg from the closed forms, make_figures.py). Units |M| = 1. Draws with ../nonradial/data_graphs.py (the arriba
graph-data skill's primitives, unchanged). Exit 1 if the sign change of kappa_c R in kappa_c.csv is not within one grid
step of the closed-form r_*.
"""
import csv
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "nonradial"))
import data_graphs as dg  # noqa: E402

LS = (2, 3, 4, 5)
KC = list(csv.DictReader(open(os.path.join(HERE, "kappa_c.csv"))))
GR = list(csv.DictReader(open(os.path.join(HERE, "growth_rates.csv"))))
TH = {int(r["l"]): r for r in csv.DictReader(open(os.path.join(HERE, "..", "..", "paper", "figures", "thresholds.csv")))}
Rs = [float(r["R"]) for r in KC]

# r_* read from the data: the last sign change of kappa_c R from + to - above r_s, against the closed form
bad = 0
print("l   r_s       r_* (closed form)   sign change in kappa_c.csv   grid step there")
for l in LS:
    y = [float(r[f"kappa_c_x_R_l{l}"]) for r in KC]
    rs, rstar = float(TH[l]["r_s"]), float(TH[l]["r_star"])
    cross = [(Rs[i], Rs[i + 1]) for i in range(len(Rs) - 1) if Rs[i] > rs and y[i] > 0 >= y[i + 1]]
    lo, hi = cross[-1] if cross else (math.nan, math.nan)
    ok = lo <= rstar <= hi
    bad += not ok
    print(f"{l}   {rs:.6f}  {rstar:.6f}            ({lo:.5f}, {hi:.5f})            {hi - lo:.1e}  {'ok' if ok else 'MISMATCH'}")

SER = [("full_DG_k_over_2", "point mass (R -> 0)", "circle"), ("R=0.3_Neumann", "R = 0.3, Neumann", "square"),
       ("R=1.0_Neumann", "R = 1, Neumann", "diamond"), ("R=10.0_kappa=-10", "R = 10, kappa = -10", "triangle")]
print("\nl   " + "   ".join(f"{lab:>20s}" for _, lab, _ in SER))
for d in GR:
    print(f"{d['l']:<3s} " + "   ".join(f"{float(d[k]):20.4f}" for k, _, _ in SER))

svg = dg.Svg(760)
X0, PW, PH = 86, 600, 200
y = 28
y = dg.caption(svg, 16, y, "Theorem 2: the threshold kappa_c(R) of a finite star",
               "l = 2 to 5, |M| = 1. A star with R > r_s grows iff kappa < kappa_c(R); with kappa >= 0 it is stable right of its dot r_*")
xm, ym = dg.axes(svg, X0, y, PW, PH, (0.05, 10, "log"), (-10, 10), "star radius R (|M|)", "kappa_c(R) R",
                 xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:g}")
svg.line(X0, ym(0), X0 + PW, ym(0), dg.INK_SOFT, 1)
for i, l in enumerate(LS):
    ys = [float(r[f"kappa_c_x_R_l{l}"]) for r in KC]
    seg = []
    for k, (R, v) in enumerate(zip(Rs, ys)):
        jump = k > 0 and abs(v - ys[k - 1]) > 4
        if abs(v) > 10 or jump:
            if len(seg) > 1:
                svg.add(f'<polyline fill="none" stroke="{dg.SERIES[i]}" stroke-width="1.8" points="{" ".join(seg)}"/>')
            seg = []
            if abs(v) > 10:
                continue
        seg.append(f"{xm(R):.1f},{ym(v):.1f}")
    if len(seg) > 1:
        svg.add(f'<polyline fill="none" stroke="{dg.SERIES[i]}" stroke-width="1.8" points="{" ".join(seg)}"/>')
    rstar = float(TH[l]["r_star"])
    dg.marker(svg, "circle", xm(rstar), ym(0), dg.SERIES[i], f"l = {l}: r_* = {rstar:.6f} |M| (r_s = {float(TH[l]['r_s']):.6f})")
dg.legend(svg, X0 + PW, y - 10, [(f"l = {l}", dg.SERIES[i]) for i, l in enumerate(LS)])
dg.callout(svg, xm(float(TH[2]["r_star"])), ym(0), f"l = 2: r_* = {float(TH[2]['r_star']):.3f} |M|", dx=40, dy=-60)

y1 = y + PH + 86
y = dg.caption(svg, 16, y1, "Any finite radius caps the point mass's l^4 growth",
               "fastest growth rate against l: k/2 ~ l^4 for the point mass, bounded or zero for a finite star (zeros not drawn)")
xm2, ym2 = dg.axes(svg, X0, y, PW, PH, (1.5, 8.5), (0.5, 1000, "log"), "angular number l", "growth rate gamma |M|",
                   xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:g}")
PTS = {i: {int(d["l"]): float(d[key]) for d in GR if float(d[key]) > 0} for i, (key, _, _) in enumerate(SER)}
POS = dg.dodge_markers({i: {a: (xm2(a), ym2(b)) for a, b in pts.items()} for i, pts in PTS.items()}, step=8.0)
for i, (key, lab, shape) in enumerate(SER):
    col = dg.SERIES[i]
    if len(PTS[i]) > 1:
        svg.add(f'<polyline fill="none" stroke="{col}" stroke-width="1.5" points="'
                + " ".join(f"{xm2(a):.1f},{ym2(b):.1f}" for a, b in PTS[i].items()) + '"/>')
    for a, b in PTS[i].items():
        dg.marker(svg, shape, POS[i][a][0], POS[i][a][1], col, f"{lab}, l = {a}: gamma |M| = {b:.4f}")
dg.legend(svg, X0 + PW, y - 10, [(lab, dg.SERIES[i]) for i, (_, lab, _) in enumerate(SER)])
last = GR[-1]
dg.callout(svg, xm2(int(last["l"])), ym2(float(last["full_DG_k_over_2"])), f"point mass, l = {last['l']}: {float(last['full_DG_k_over_2']):g}",
           dx=-150, dy=10)

try:
    commit = subprocess.check_output(["git", "-C", HERE, "rev-parse", "--short", "HEAD"], text=True).strip()
except (subprocess.SubprocessError, OSError):
    commit = "?"
height = dg.footnote(svg, [f"thresholds_graph.py at negative-mass {commit}", "kappa_c.csv, growth_rates.csv: truncated.py",
                           "r_s, r_*: paper/figures/thresholds.csv"],
                     "closed forms (kappa_c) and Chebyshev spectra (growth rates) committed in phase2/dg_spectrum", [],
                     y1 + PH + 100, repo_url="https://github.com/nmoinvaz/negative-mass")
open(os.path.join(HERE, "thresholds.svg"), "w").write(svg.finish(height))
print("\nwrote thresholds.svg")
sys.exit(1 if bad else 0)

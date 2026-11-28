#!/usr/bin/env python3
"""Figure of #267: the settled mass of family G against delta, the prediction from #266's relation and the residual, and the member's throat radius along the stall (h = 0.005).
usage: python3 settled_mass_member_graph.py [RAWDIR [OUT.svg]]     (RAWDIR defaults to ~/Source/negative-mass-raw/I267-mac, OUT to settled_mass_member.svg here; numbers from the same functions as settled_mass_member.py analyze)
The graph primitives are phase2/nonradial/data_graphs.py (the arriba graph-data skill's), unchanged."""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "nonradial"))
import data_graphs as dg                  # noqa: E402
raw = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Source/negative-mass-raw/I267-mac")
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "settled_mass_member.svg")
sys.argv = ["x"]
import settled_mass_member as mm          # noqa: E402

fin = 0.005
D = mm.derive(mm.load_all(raw))
R = mm.load266(fin)
c0, ab = mm.fit_relation(R)
rows = []
for e in mm.E_ALL:
    d = D[(e, fin)]
    q = d["mem"]["m3"]
    M = d["SS"]["plat1"]
    rows.append(dict(e=e, M=M, g1=q["g1"], r=q["r_th"], t=q["t_dep"], Mp=mm.predict(c0, ab, q["g1"], q["r_th"]), ph=e not in mm.E_GATE))
print("e   t_dep    gamma_1   r_th      M         M predicted  M - predicted")
for r in rows:
    print(f"{r['e']:<3d} {r['t']:7.3f}  {r['g1']:.6f}  {r['r']:.6f}  {r['M']:.6f}  {r['Mp']:.6f}    {r['M'] - r['Mp']:+.2e}")
print(f"relation: M_0 = {c0:.6f}, (a, b) = ({ab[0]:.4f}, {ab[1]:.4f})")

svg = dg.Svg(760)
X0, PW, PH = 86, 600, 170
y = 28
y = dg.caption(svg, 16, y, "Settled mass of family G against delta, and #266's prediction from the member at departure",
               "M_AH plateau at e^-f = 1.57e6, D = 4, B = 0.1, h = 0.005; a = 1. Circles: the six gating delta; diamonds: delta = 1e-3 to 1e-6 (the member is read before the stall settles).")
xm, ym = dg.axes(svg, X0, y, PW, PH, (2.5, 12.5), (0.2330, 0.2380), "distance from threshold delta = p*_h - A", "settled mass M_AH (a)",
                 xfmt=lambda v: f"1e-{v:g}", yfmt=lambda v: f"{v:.4f}")
for r in rows:
    sh = "diamond" if r["ph"] else "circle"
    dg.marker(svg, sh, xm(r["e"]), ym(r["M"]), dg.SERIES[0], f"delta = 1e-{r['e']}: M = {r['M']:.6f}")
    dg.marker(svg, sh, xm(r["e"]), ym(r["Mp"]), dg.SERIES[1], f"delta = 1e-{r['e']}: predicted r_th (M_0 + a gamma_1 + b gamma_1^2) = {r['Mp']:.6f} (gamma_1 {r['g1']:.5f}, r_th {r['r']:.5f})")
dg.legend(svg, X0 + PW, y - 10, [("measured", dg.SERIES[0]), ("predicted from the member at departure", dg.SERIES[1])])
r12 = rows[-1]
dg.callout(svg, xm(12), ym(r12["M"]), f"1e-12: {r12['M']:.5f}, predicted {r12['Mp']:.5f}", dx=-150, dy=28)
y1 = y + PH + 80
y = dg.caption(svg, 16, y1, "Measured minus predicted", "the declared bar is 6e-5 (lines); the relation of #266 is fitted at B = 0.025 to 0.125 at delta = 1e-8, so the point at 1e-8 is on it by construction")
res = [r["M"] - r["Mp"] for r in rows]
xm2, ym2 = dg.axes(svg, X0, y, PW, PH, (2.5, 12.5), (-2.2e-3, 2.4e-3), "distance from threshold delta = p*_h - A", "M - predicted (a)",
                   xfmt=lambda v: f"1e-{v:g}", yfmt=lambda v: f"{v:+.1e}")
for s_ in (+6e-5, -6e-5):
    svg.line(X0, ym2(s_), X0 + PW, ym2(s_), dg.SERIES[-1], 1)
for r, v in zip(rows, res):
    dg.marker(svg, "diamond" if r["ph"] else "circle", xm2(r["e"]), ym2(v), dg.SERIES[0], f"delta = 1e-{r['e']}: M - predicted = {v:+.2e}")
dg.callout(svg, xm2(12), ym2(res[-1]), f"{res[-1]:+.1e}", dx=-40, dy=-26)
y2 = y1 + PH + 112
y = dg.caption(svg, 16, y2, "The member's throat radius along the stall", "mean of the pair p*_h -/+ 1e-12; dots: the member at departure of each delta (t_dep, r_th); the radius falls to 1.0029 at t = 13.5 and then rises")
d12 = D[(12, fin)]
t, tau, rr, odd = mm.pair_arrays(d12["loc"])
sel = (t >= 8.0) & (np.abs(odd) <= 3e-4)
xm3, ym3 = dg.axes(svg, X0, y, PW, PH, (8, 24), (1.000, 1.011), "time t = (x+ + x-)/sqrt(2) (a)", "throat radius r_th (a)", xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:.3f}")
svg.add('<polyline fill="none" stroke="' + dg.SERIES[0] + '" stroke-width="1.5" points="' + " ".join(f"{xm3(a):.1f},{ym3(b_):.1f}" for a, b_ in zip(t[sel][::5], rr[sel][::5])) + '"/>')
for r in rows:
    if not r["ph"]:
        dg.marker(svg, "circle", xm3(r["t"]), ym3(r["r"]), dg.SERIES[1], f"delta = 1e-{r['e']}: t_dep = {r['t']:.3f}, r_th = {r['r']:.5f}, gamma_1 = {r['g1']:.5f}")
dg.legend(svg, X0 + PW, y - 10, [("pair 1e-12 mean", dg.SERIES[0]), ("member at departure", dg.SERIES[1])])
dg.callout(svg, xm3(13.5), ym3(1.002866), "minimum 1.00287 at t = 13.5", dx=-30, dy=-40)
try:
    commit = subprocess.check_output(["git", "-C", HERE, "rev-parse", "--short", "HEAD"], text=True).strip()
except (subprocess.SubprocessError, OSError):
    commit = "?"
height = dg.footnote(svg, [f"settled_mass_member.py at negative-mass {commit}", "chain: #41 via #253", "member: #225 via #266", "prediction: #266's relation"], dg.host_specs(),
                     ["PRED fails: the prediction misses by up to 9e-4 against its bar of 6e-5", "delta = 1e-12 is at the roundoff floor of the pair (1 to 3 % in its odd part)"],
                     y2 + PH + 100, repo_url="https://github.com/nmoinvaz/negative-mass")
open(out, "w").write(svg.finish(height))
print("wrote", out)

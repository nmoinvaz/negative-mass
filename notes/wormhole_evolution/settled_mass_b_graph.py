#!/usr/bin/env python3
"""Figure of #266: the settled mass of family G against B, and M/r_th - M_0 against gamma_1 (h = 0.005, delta = 1e-8).
usage: python3 settled_mass_b_graph.py [RAWDIR [OUT.svg]]     (RAWDIR defaults to ~/Source/negative-mass-raw/I266-mac, OUT to settled_mass_b.svg here; numbers from the same functions as settled_mass_b.py analyze)
The graph primitives are phase2/nonradial/data_graphs.py (the arriba graph-data skill's), unchanged."""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "nonradial"))
import data_graphs as dg                  # noqa: E402
raw = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Source/negative-mass-raw/I266-mac")
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "settled_mass_b.svg")
sys.argv = ["x"]
import settled_mass_b as b                # noqa: E402

D = b.load_all(raw)
fin = 0.005
rows = []
for fam, B in b.TAGS:
    if not b.appl(fam, B) or fam == "O":
        continue
    tg = b.tag_of(fam, B)
    e = D[(tg, fin)]
    s, _ = e["s"], e["i"]
    M = b.g.summary_g(s, fin)["plat1"]
    q = b.member_at(e, 1e-3, "m")
    rows.append(dict(B=B, M=M, g1=q["g1"], r=q["r_th"], ph=(fam, B) in b.PH_TAGS))
M0 = rows[0]["M"]
for r in rows:
    r["dMp"] = r["M"] / r["r"] - M0
    r["dM"] = r["M"] - M0
pos = [r for r in rows if r["B"] > 0]
A = np.array([[r["g1"], r["g1"] ** 2] for r in pos])
a, c = np.linalg.lstsq(A, np.array([r["dMp"] for r in pos]), rcond=None)[0]
print("B      gamma_1   r_th     M         M/r_th    dM'       dM       fit(dM')")
for r in rows:
    print(f"{r['B']:<6g} {r['g1']:+.5f}  {r['r']:.5f}  {r['M']:.6f}  {r['M'] / r['r']:.6f}  {r['dMp']:+.6f}  {r['dM']:+.6f}  {a * r['g1'] + c * r['g1'] ** 2:+.6f}")
print(f"dM' = {a:.4f} gamma_1 + {c:.4f} gamma_1^2 (least squares, four B)")

svg = dg.Svg(760)
W, X0, PW, PH = 760, 86, 600, 190
y = 28
y = dg.caption(svg, 16, y, "Settled mass of family G against the normal amplitude B", "M_AH plateau at e^-f = 1.57e6, D = 4, delta = 1e-8, h = 0.005; a = 1. Diamond: B = 0.125, post hoc. B = 0.2 has no stall.")
xm, ym = dg.axes(svg, X0, y, PW, PH, (0, 0.14), (0.222, 0.250), "normal amplitude B", "settled mass M_AH (a)", xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:.3f}")
for r in rows:
    dg.marker(svg, "diamond" if r["ph"] else "circle", xm(r["B"]), ym(r["M"]), dg.SERIES[0], f"B = {r['B']:g}: M = {r['M']:.6f}")
svg.text(X0 + PW - 6, y + PH - 10, "B = 0.2: no stall, no settled mass (chain stops at 'horizon lost')", size=11, fill=dg.INK, anchor="end")
dg.callout(svg, xm(0), ym(M0), f"B = 0: {M0:.6f} (S: 0.22309)", dx=40, dy=-28)
dg.callout(svg, xm(0.1), ym(rows[3]["M"]), f"B = 0.1: {rows[3]['M']:.6f} (#253)", dx=-30, dy=-36)
y1 = y + PH + 66
y = dg.caption(svg, 16, y1, "M / r_th - M_0 against the asymmetry gamma_1 of the member at departure", "line: dM' = a gamma_1 + b gamma_1^2 through the four B; gray: the same without the r_th scaling")
xm2, ym2 = dg.axes(svg, X0, y, PW, PH, (0, 0.085), (0, 0.025), "asymmetry gamma_1 of the member at departure", "M / r_th - M_0 (a)", xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:.3f}")
gs = np.linspace(0, 0.085, 60)
svg.add('<polyline fill="none" stroke="' + dg.SERIES[0] + '" stroke-width="1.5" points="' + " ".join(f"{xm2(g):.1f},{ym2(a * g + c * g * g):.1f}" for g in gs) + '"/>')
for r in rows[1:]:
    sh = "diamond" if r["ph"] else "circle"
    dg.marker(svg, sh, xm2(r["g1"]), ym2(r["dM"]), dg.SERIES[-1], f"B = {r['B']:g}: M - M_0 = {r['dM']:.6f} (no r_th scaling)")
    dg.marker(svg, sh, xm2(r["g1"]), ym2(r["dMp"]), dg.SERIES[0], f"B = {r['B']:g}: gamma_1 = {r['g1']:.5f}, r_th = {r['r']:.5f}, M/r_th - M_0 = {r['dMp']:.6f}")
dg.legend(svg, X0 + PW, y - 10, [("M / r_th - M_0", dg.SERIES[0]), ("M - M_0, no r_th scaling", dg.SERIES[-1])])
dg.callout(svg, xm2(rows[3]["g1"]), ym2(rows[3]["dMp"]), f"B = 0.1: gamma_1 = {rows[3]['g1']:.4f}, r_th = {rows[3]['r']:.4f}", dx=-40, dy=-40)
try:
    commit = subprocess.check_output(["git", "-C", HERE, "rev-parse", "--short", "HEAD"], text=True).strip()
except (subprocess.SubprocessError, OSError):
    commit = "?"
height = dg.footnote(svg, [f"settled_mass_b.py at negative-mass {commit}", "chain: #41 via #253", "member: #225"], dg.host_specs(),
                     ["post hoc: B = 0.125 (diamond) and h = 0.02 gate nothing", "the gating leave-one-out over three B misses its bar (7.9e-5 against 6e-5)"],
                     y1 + PH + 100, repo_url="https://github.com/nmoinvaz/negative-mass")
open(out, "w").write(svg.finish(height))
print("wrote", out)

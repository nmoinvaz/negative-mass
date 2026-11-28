#!/usr/bin/env python3
"""Figure of #295: the frozen refit's residuals on the fifteen fitted points (in-sample and leave-one-out) and on the three members of the third route, and where the third route
sits between routes W and C (h = 0.005).
usage: python3 settled_mass_route3_graph.py [RAWDIR [OUT.svg]]     (RAWDIR defaults to ~/Source/negative-mass-raw/I295-mac, OUT to settled_mass_route3.svg here; numbers from the same
functions as settled_mass_route3.py analyze)
The graph primitives are phase2/nonradial/data_graphs.py (the arriba graph-data skill's), unchanged."""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "nonradial"))
import data_graphs as dg                  # noqa: E402
raw = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Source/negative-mass-raw/I295-mac")
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "settled_mass_route3.svg")
sys.argv = ["x"]
import settled_mass_route3 as s3          # noqa: E402

r = s3.r
fin = s3.H_FIN
D = {}
for nm, spec, tg in s3.POINTS3:
    e = r.load_point(os.path.join(raw, "runs"), spec, fin)
    e["name"] = nm
    D[(nm, fin)] = e
r.derive(D)
pts, c0 = s3.points15()
g1, rth, rdot, M = (np.array([p[k] for p in pts]) for k in ("g1", "rth", "rdot", "M"))
co = s3.fit(g1, rth, rdot, M, "P4")
ins = M - rth * (s3.M0 + s3.design(g1, rdot, "P4") @ co)
eloo, _ = s3.loo(g1, rth, rdot, M, "P4")
new = []
for nm, spec, tg in s3.POINTS3:
    e = D[(nm, fin)]
    Mn, gn, rn_, dn_ = e["SS"]["plat1"], e["mem"]["g1"], e["mem"]["r_th"], e["sp"]
    new.append(dict(nm=nm, g1=gn, r=rn_, rd=dn_, M=Mn, res=Mn - s3.predict("P4", gn, rn_, dn_, fin, mutated=False)))
print("point  gamma_1   r_th      rdot_th    M - P4 (new)")
for p in new:
    print(f"{p['nm']:5s}  {p['g1']:.5f}  {p['r']:.5f}  {p['rd']:+.2e}  {p['res']:+.2e}")

refs = {"B": r.fam_b(fin), "W": s3.route_curve(fin, ("W1", "W2", "W3")), "C": s3.route_curve(fin, ("C1", "C2", "C3"))}
refs["B"] = dict(g1=refs["B"]["g1"], rth=refs["B"]["rth"], rdot=refs["B"]["rdot"])


def offset(g, r_, v_):
    """(r_th, rdot_th) of a point minus family B's at the same gamma_1 (piecewise linear)."""
    q = s3.sep3(g, (r_, v_), (r_, v_), refs["B"], k=0.0)["q"]
    return q["r_th"]["d_lin"], q["rdot"]["d_lin"]


svg = dg.Svg(760)
X0, PW, PH = 86, 600, 170
ROUTE = {"B": ("square", "#9aa0a6"), "D": ("triangle", dg.SERIES[3]), "W": ("circle", dg.SERIES[0]), "C": ("diamond", dg.SERIES[1])}
y = 28
y = dg.caption(svg, 16, y, "Measured minus the frozen refit, against the member's asymmetry gamma_1",
               "h = 0.005; small filled: in-sample, open: leave-one-out, of the fifteen fitted points; large green: the three new members")
xm, ym = dg.axes(svg, X0, y, PW, PH, (0.0, 0.08), (-1.0e-4, 1.0e-4), "gamma_1 at departure", "M - predicted", xfmt=lambda v: f"{v:.2f}", yfmt=lambda v: f"{v:+.0e}")
for s_ in (+6e-5, -6e-5):
    svg.line(X0, ym(s_), X0 + PW, ym(s_), dg.SERIES[-1], 1)
svg.line(X0, ym(0.0), X0 + PW, ym(0.0), "#9aa0a6", 1)
for p, a_, b_ in zip(pts, ins, eloo):
    sh, col = ROUTE[p["route"]]
    dg.marker(svg, sh, xm(p["g1"]), ym(a_), col, f"{p['name']} (fitted), in-sample M - P4 = {a_:+.2e}")
    svg.add(f'<circle cx="{xm(p["g1"]):.1f}" cy="{ym(b_):.1f}" r="5" fill="none" stroke="{col}" stroke-width="1.5"><title>{dg.esc(p["name"])} (fitted), leave-one-out M - P4 = {b_:+.2e}</title></circle>')
for p in new:
    svg.add(f'<circle cx="{xm(p["g1"]):.1f}" cy="{ym(p["res"]):.1f}" r="8" fill="{dg.SERIES[2]}" stroke="{dg.INK}" stroke-width="1.5"><title>{dg.esc(p["nm"])} (new): gamma_1 {p["g1"]:.5f}, M - P4 = {p["res"]:+.2e}</title></circle>')
dg.legend(svg, X0 + PW, y - 10, [("B", "#9aa0a6"), ("delta scan", dg.SERIES[3]), ("route W", dg.SERIES[0]), ("route C", dg.SERIES[1]), ("new route", dg.SERIES[2])])

y1 = y + PH + 80
y = dg.caption(svg, 16, y1, "Where the third route sits: r_th and rdot_th at departure minus family B's at the same gamma_1",
               "h = 0.005; each route's three points at gamma_1 = 0.012, 0.048, 0.075 (labelled 1, 2, 3); the new route lies near the chord of routes W and C")
rr = [offset(p["g1"], p["r"], p["rd"]) for p in new]
allx, ally = [], []
curves = {}
for key, cv in (("W", refs["W"]), ("C", refs["C"])):
    pp = [(cv["g1"][i], cv["rth"][i], cv["rdot"][i]) for i in range(1, len(cv["g1"]))]
    curves[key] = [offset(*q) for q in pp]
    allx += [q[0] for q in curves[key]]
    ally += [q[1] for q in curves[key]]
allx += [q[0] for q in rr]
ally += [q[1] for q in rr]
xl = (min(allx) * 1.15 - 2e-3, max(allx) * 1.15 + 2e-3)
yl = (min(ally) * 1.15 - 1e-4, max(ally) * 1.15 + 1e-4)
xm2, ym2 = dg.axes(svg, X0, y, PW, PH, xl, yl, "r_th - family B's", "rdot_th - family B's", xfmt=lambda v: f"{v:+.0e}", yfmt=lambda v: f"{v:+.0e}")
for key, sh, col in (("W", "circle", dg.SERIES[0]), ("C", "diamond", dg.SERIES[1])):
    for k, (dx, dy) in enumerate(curves[key]):
        dg.marker(svg, sh, xm2(dx), ym2(dy), col, f"route {key}{k + 1}: r_th - B {dx:+.2e}, rdot_th - B {dy:+.2e}")
        svg.text(xm2(dx) + 8, ym2(dy) - 6, f"{key}{k + 1}", size=11)
for k, (dx, dy) in enumerate(rr):
    dg.marker(svg, "triangle", xm2(dx), ym2(dy), dg.SERIES[2], f"N{k + 1}: r_th - B {dx:+.2e}, rdot_th - B {dy:+.2e}")
    svg.text(xm2(dx) + 8, ym2(dy) - 6, f"N{k + 1}", size=11)
for k in range(3):
    xs = [curves["W"][k][0], curves["C"][k][0]]
    ys = [curves["W"][k][1], curves["C"][k][1]]
    svg.line(xm2(xs[0]), ym2(ys[0]), xm2(xs[1]), ym2(ys[1]), "#c8ccd0", 1)
dg.legend(svg, X0 + PW, y - 10, [("route W", dg.SERIES[0]), ("route C", dg.SERIES[1]), ("new route", dg.SERIES[2])])
try:
    commit = subprocess.check_output(["git", "-C", HERE, "rev-parse", "--short", "HEAD"], text=True).strip()
except (subprocess.SubprocessError, OSError):
    commit = "?"
height = dg.footnote(svg, [f"settled_mass_route3.py at negative-mass {commit}", "chain: #41 via #253", "member: #225 via #266", "speed: #267, #255", "refit: #295, frozen before the runs"], dg.host_specs(),
                     [f"the frozen refit's largest leave-one-out error is {np.max(np.abs(eloo)):.1e}; the new members' largest residual is {max(abs(p['res']) for p in new):.1e}"], y1 + PH + 100,
                     repo_url="https://github.com/nmoinvaz/negative-mass")
open(out, "w").write(svg.finish(height))
print("wrote", out)

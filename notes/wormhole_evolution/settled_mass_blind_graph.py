#!/usr/bin/env python3
"""Figure of #317: the miss of the frozen refits at X1 to X3 under six readings, and the same V3 procedure repeated at other eps (post hoc), h = 0.005.
usage: python3 settled_mass_blind_graph.py [OUT.svg]     (numbers from the same functions as settled_mass_blind.py analyze and post; stored runs read from ~/Source/negative-mass-raw)
The graph primitives are phase2/nonradial/data_graphs.py (the arriba graph-data skill's), unchanged."""
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "nonradial"))
import data_graphs as dg                  # noqa: E402
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "settled_mass_blind.svg")
sys.argv = ["x"]
import settled_mass_blind as B            # noqa: E402

E = B.E
fz = json.load(open(B.frozen_path()))
B.seed_x()
train, test, N3 = B.FIT, B.XN, B.N3
names = (B.REF,) + train + test + N3
P1, P2 = B.readset(names, B.EPS_A), B.readset(names, B.EPS_B)

cols = [("P4, 1e-3 (#305)", None), ("refit, 1e-3", "FIX_1e-3"), ("V3, 1e-3", "V3_1e-3"), ("refit, 1e-2", "FIX_1e-2"), ("V1, 1e-2", "V1_1e-2"), ("V3, 1e-2", "V3_1e-2")]
res = {}
for lab, key in cols:
    if key is None:
        f = dict(eps=B.EPS_A, kind="P3", sub=False, co=list(E.FROZEN_P4), b0=None, m0=B.M0)
        res[lab] = B.resid(f, P1, test)[0]
    else:
        f = fz["fits"][key]
        res[lab] = B.resid(f, P2 if f["eps"] == B.EPS_B else P1, test)[0]
print("reading            X1         X2         X3")
for lab, _ in cols:
    print(f"{lab:16s} " + "  ".join(f"{x:+.3e}" for x in res[lab]))

scan = (1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2)
v3, nmax = [], []
for eps in scan:
    P = B.readset(names, eps)
    Pq = B.apply_offsets(P, **B.offsets_of(P[B.REF]))
    v3.append(E.fit_pred(Pq, train, test, B.M0)["res"])
    rn = E.fit_pred(Pq, train, N3, B.M0)["res"]
    nmax.append(float(rn[int(np.argmax(np.abs(rn)))]))
print("eps      V3: X1         X2         X3      | N largest")
for eps, r_, n_ in zip(scan, v3, nmax):
    print(f"{eps:<7g} " + "  ".join(f"{x:+.3e}" for x in r_) + f" | {n_:+.2e}")

svg = dg.Svg(760)
X0, PW, PH = 86, 600, 200
SH = ("circle", "square", "triangle")
y = 28
y = dg.caption(svg, 16, y, "The miss at X1 to X3 under six readings of the corrected and the earlier fit",
               "h = 0.005; M - predicted of fits frozen over 14 points (P4: #295's, fifteen points); lines: bar 6e-5; V3 at 1e-2 is the gated prediction, the others are reported only")
xm, ym = dg.axes(svg, X0, y, PW, PH, (-0.5, len(cols) - 0.5), (-1.6e-4, 2.4e-4), "reading (fit, eps)", "M - predicted", xfmt=lambda v: "", yfmt=lambda v: f"{v:+.0e}")
for j, (lab, _) in enumerate(cols):
    svg.text(xm(j), y + PH + 16, lab, size=11, anchor="middle")
for s_ in (+6e-5, -6e-5):
    svg.line(X0, ym(s_), X0 + PW, ym(s_), dg.SERIES[-1], 1)
svg.line(X0, ym(0.0), X0 + PW, ym(0.0), "#9aa0a6", 1)
for j, (lab, _) in enumerate(cols):
    for q, n in enumerate(test):
        dg.marker(svg, SH[q], xm(j) + (q - 1) * 12, ym(res[lab][q]), dg.SERIES[q], f"{n} under {lab}: M - predicted = {res[lab][q]:+.2e}")
dg.callout(svg, xm(5) + 0, ym(res["V3, 1e-2"][2]), "gated: X2, X3 beyond the bar", dx=-120, dy=-34)
dg.legend(svg, X0 + PW, y - 10, [("X1", dg.SERIES[0]), ("X2", dg.SERIES[1]), ("X3", dg.SERIES[2])])

y1 = y + PH + 80
y = dg.caption(svg, 16, y1, "The V3 procedure repeated at other eps (post hoc; the frozen eps is 1e-2)",
               "V3 refit over the 14 at each eps; grey: largest of N1 to N3; read after the record, so no eps may be chosen from it")
xm2, ym2 = dg.axes(svg, X0, y, PW, PH, (1e-3, 5e-2, "log"), (-1.0e-4, 2.0e-4), "departure threshold eps (|odd| / r)", "M - predicted",
                   xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:+.0e}")
for s_ in (+6e-5, -6e-5):
    svg.line(X0, ym2(s_), X0 + PW, ym2(s_), dg.SERIES[-1], 1)
svg.line(X0, ym2(0.0), X0 + PW, ym2(0.0), "#9aa0a6", 1)
for q, n in enumerate(test):
    for a, b_ in zip(range(len(scan) - 1), range(1, len(scan))):
        svg.line(xm2(scan[a]), ym2(v3[a][q]), xm2(scan[b_]), ym2(v3[b_][q]), dg.SERIES[q], 1.5)
    for k_, eps in enumerate(scan):
        dg.marker(svg, SH[q], xm2(eps), ym2(v3[k_][q]), dg.SERIES[q], f"{n}, eps {eps:g}: M - predicted = {v3[k_][q]:+.2e}")
for k_, eps in enumerate(scan):
    svg.add(f'<circle cx="{xm2(eps):.1f}" cy="{ym2(nmax[k_]):.1f}" r="4" fill="#9aa0a6"><title>N1 to N3 at eps {eps:g}: largest M - predicted {nmax[k_]:+.2e}</title></circle>')
dg.callout(svg, xm2(1e-2), ym2(v3[3][2]), "frozen eps", dx=40, dy=-30)
dg.legend(svg, X0 + PW, y - 10, [("X1", dg.SERIES[0]), ("X2", dg.SERIES[1]), ("X3", dg.SERIES[2])])
try:
    commit = subprocess.check_output(["git", "-C", HERE, "rev-parse", "--short", "HEAD"], text=True).strip()
except (subprocess.SubprocessError, OSError):
    commit = "?"
height = dg.footnote(svg, [f"settled_mass_blind.py at negative-mass {commit}", "readings: #296", "chain: #41", "member: #266", "speed: #267", "X: #305"], dg.host_specs(),
                     [f"V3 at eps = 1e-2: M - predicted {res['V3, 1e-2'][0]:+.2e}, {res['V3, 1e-2'][1]:+.2e}, {res['V3, 1e-2'][2]:+.2e} against the bar 6e-5"],
                     y1 + PH + 100, repo_url="https://github.com/nmoinvaz/negative-mass")
open(out, "w").write(svg.finish(height))
print("wrote", out)

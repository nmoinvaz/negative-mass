#!/usr/bin/env python3
"""Figure of #291: the throat speed of family G at D = 4, B = 0.02 against proper time, and its local exponent -d ln|v| / d ln tau, with the guides of a power law and of an exponential.
usage: python3 throat_tail_graph.py [RAWDIR [OUT.svg]]     (RAWDIR defaults to ~/Source/negative-mass-raw/I291-mac, OUT to throat_tail.svg here; numbers from the same functions as throat_tail.py analyze)
The graph primitives are phase2/nonradial/data_graphs.py (the arriba graph-data skill's), unchanged."""
import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "nonradial"))
import data_graphs as dg                  # noqa: E402
raw = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Source/negative-mass-raw/I291-mac")
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "throat_tail.svg")
sys.argv = ["x"]
import throat_tail as tt                  # noqa: E402

C = tt.CRIT
M = tt.load_all(raw)
S = {h: {k: tt.series_member(m) for k, m in M[h].items()} for h in M}
SM = {h: tt.ens_series([S[h][k] for k in sorted(S[h])]) for h in tt.H2}
W = {h: tt.window([S[h][k] for k in sorted(S[h])]) for h in tt.H2}
T_end = min(w["T_end"] for w in W.values())
st = {h: tt.tail_stats(SM[h], T_end) for h in tt.H2}
print(f"T_end = {T_end:.2f}; slopes " + ", ".join(f"h = {h:g}: -{st[h]['p']:.3f}" for h in tt.H2))

tg = np.arange(8.0, 20.0 + 1e-9, 0.1)
H0, H1 = 0.005, 0.01
tau_of = lambda Sx, t: np.interp(t, Sx["t"], Sx["tau"])      # noqa: E731
v15 = abs(float(tt.v_on(SM[H0], np.array([15.0]))[0]))
tau15 = float(tau_of(SM[H0], 15.0))

svg = dg.Svg(760)
X0, PW, PH = 86, 600, 230
y = 28
y = dg.caption(svg, 16, y, "Throat speed of family G at D = 4, B = 0.02 against proper time",
               "mean of K = 8 pairs p*_h -/+ 1e-12 (colored), the single pairs at h = 0.005 (grey); guides through the mean at t = 15; log axes")
xm, ym = dg.axes(svg, X0, y, PW, PH, (8, 20, "log"), (1.5e-5, 2.0e-4, "log"), "proper time tau along the throat (a)", "throat speed |dr_th / dtau|",
                 xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:.0e}")


def poly(xs, ys, color, dash=None, w=1.6):
    pts = " ".join(f"{xm(a):.1f},{ym(b_):.1f}" for a, b_ in zip(xs, ys) if np.isfinite(b_) and 1.5e-5 <= b_ <= 2.0e-4)
    svg.add('<polyline fill="none" stroke="' + color + f'" stroke-width="{w}"' + (f' stroke-dasharray="{dash}"' if dash else "") + ' points="' + pts + '"/>')


for k in sorted(S[H0]):
    sel = (tg >= 8.0) & (tg <= 20.0)
    poly(tau_of(S[H0][k], tg[sel]), np.abs(tt.v_on(S[H0][k], tg[sel])), dg.SERIES[-1], w=0.8)
tx = np.linspace(8.0, 20.0, 100)
poly(tx, v15 * (tx / tau15) ** -3.0, dg.SERIES[2], dash="5 4", w=1.4)
poly(tx, v15 * np.exp(-C["kappa255"] * (tx - tau15)), dg.SERIES[3], dash="2 3", w=1.6)
for h, col in ((H1, dg.SERIES[1]), (H0, dg.SERIES[0])):
    sel = tg <= max(W[h]["T_end"], T_end) + 0.4
    poly(tau_of(SM[h], tg[sel]), np.abs(tt.v_on(SM[h], tg[sel])), col, w=1.9)
svg.line(xm(T_end), y, xm(T_end), y + PH, dg.SERIES[-1], 1)
dg.legend(svg, X0 + PW, y - 10, [("mean, h = 0.005", dg.SERIES[0]), ("mean, h = 0.01", dg.SERIES[1]), ("tau^-3 through t = 15", dg.SERIES[2]), ("exp(-0.244 tau)", dg.SERIES[3]), ("single pairs, h = 0.005", dg.SERIES[-1])])
dg.callout(svg, xm(tau15), ym(v15), f"t = 15: {v15:.2e}", dx=-60, dy=-58)
dg.callout(svg, xm(T_end), ym(float(abs(tt.v_on(SM[H0], np.array([T_end]))[0]))), f"window ends at {T_end:.1f}", dx=-70, dy=60)

y1 = y + PH + 84
y = dg.caption(svg, 16, y1, "Local exponent rate x t = -d ln|v| / d ln tau", "central difference over +-0.5 in t; grey lines: the S1 bar 2.9 +- 20 % and the window [15, T_end]; dotted: an exponential's 0.244 t")
xm2, ym2 = dg.axes(svg, X0, y, PW, PH, (8, 20), (1.0, 4.2), "time t = (x+ + x-)/sqrt(2) (a)", "rate x t", xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:g}")
lo, hi = C["p_target"] * (1 - C["p_tol"]), C["p_target"] * (1 + C["p_tol"])
for yy in (lo, hi):
    svg.line(X0, ym2(yy), X0 + PW, ym2(yy), dg.SERIES[-1], 1)
svg.line(xm2(C["t0"]), y, xm2(C["t0"]), y + PH, dg.SERIES[-1], 1)
svg.line(xm2(T_end), y, xm2(T_end), y + PH, dg.SERIES[-1], 1)
svg.add('<polyline fill="none" stroke="' + dg.SERIES[3] + '" stroke-width="1.6" stroke-dasharray="2 3" points="' + " ".join(f"{xm2(a):.1f},{ym2(C['kappa255'] * a):.1f}" for a in (8, 12, 16)) + '"/>')
for h, col in ((H1, dg.SERIES[1]), (H0, dg.SERIES[0])):
    tt_ = np.arange(8.5, 19.0 + 1e-9, 0.25)
    pp = [tt.local_p(SM[h], T) for T in tt_]
    sel = [i for i, (T, p) in enumerate(zip(tt_, pp)) if -0.5 < p < 4.2 and T <= (T_end + 0.5)]
    svg.add('<polyline fill="none" stroke="' + col + '" stroke-width="1.9" points="' + " ".join(f"{xm2(tt_[i]):.1f},{ym2(pp[i]):.1f}" for i in sel) + '"/>')
for T in (10.0, 12.0, 14.0, 16.0):
    dg.marker(svg, "circle", xm2(T), ym2(tt.local_p(SM[H0], T)), dg.SERIES[0], f"h = 0.005, t = {T:g}: rate x t = {tt.local_p(SM[H0], T):.3f}")
dg.marker(svg, "diamond", xm2(15.0), ym2(2.74), dg.SERIES[2], "#272 linear response at t = 15: 2.74")
dg.legend(svg, X0 + PW, y - 10, [("mean, h = 0.005", dg.SERIES[0]), ("mean, h = 0.01", dg.SERIES[1]), ("exponential 0.244 t", dg.SERIES[3]), ("linear response (#272)", dg.SERIES[2])])
dg.callout(svg, xm2(16.0), ym2(tt.local_p(SM[H0], 16.0)), f"2.3 at t = 15 to 16", dx=-90, dy=60)
try:
    commit = subprocess.check_output(["git", "-C", HERE, "rev-parse", "--short", "HEAD"], text=True).strip()
except (subprocess.SubprocessError, OSError):
    commit = "?"
height = dg.footnote(svg, [f"throat_tail.py at {commit}", "evolution #39, #225", "linear response #272"], dg.host_specs(),
[f"S1 slope -{st[H0]['p']:.2f} (h = 0.005), -{st[H1]['p']:.2f} (h = 0.01), bar 2.32 to 3.48", "single pairs: roundoff offset after t = 17"], y1 + PH + 100, repo_url="https://github.com/nmoinvaz/negative-mass")
open(out, "w").write(svg.finish(height))
print("wrote", out)

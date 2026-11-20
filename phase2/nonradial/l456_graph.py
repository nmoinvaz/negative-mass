"""
Graph of the phantom star's mode condition D(gamma) at l = 4, 5, 6 (#4, RESULTS.md section 3a).

Reads scan_eps-1_w0.7603_l{4,5,6}.csv, convergence_l{4,5,6}_output.txt and kappa_eff_l{4,5,6}_output.txt (for |M|),
writes l456_mode_condition.svg and prints the same numbers as a table.  Usage: python l456_graph.py
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import data_graphs
from data_graphs import INK, INK_SOFT, SERIES, Svg, axes, callout, caption, footnote, legend, marker

LS = (4, 5, 6)


def decade_ticks(lo, hi):
    # Whole decades, every second one past 8 decades, so the 15 decades of D stay legible
    e0, e1 = math.ceil(math.log10(lo) - 1e-9), math.floor(math.log10(hi) + 1e-9)
    step = 2 if e1 - e0 > 8 else 1
    return [10.0 ** e for e in range(e1, e0 - 1, -step)][::-1]


data_graphs.nice_log_ticks = decade_ticks


def fmt_dec(v):
    e = round(math.log10(v))
    return f"{v:g}" if -3 <= e <= 3 else f"1e{e}"


def scan(l):
    rows = [ln.split(",") for ln in open(os.path.join(HERE, f"scan_eps-1_w0.7603_l{l}.csv")).read().split("\n")[1:] if ln]
    return [(float(g), float(d)) for g, d in rows]


def conv(l):
    out = []
    for ln in open(os.path.join(HERE, f"convergence_l{l}_output.txt")):
        m = re.match(r"gamma=\s*(\S+): D_base=(\S+)\s+signs=\[(.*?)\]\s+rel\.spread \(same r_m\)=(\S+)", ln)
        if m:
            out.append((float(m.group(1)), float(m.group(2)), m.group(3), float(m.group(4)), "not completed" in ln))
    return out


def mass(l):
    return float(re.search(r"\|M\|=(\S+)", open(os.path.join(HERE, f"kappa_eff_l{l}_output.txt")).read()).group(1))


scans = {l: scan(l) for l in LS}
convs = {l: conv(l) for l in LS}
# DG point-mass rate k/(2|M|), k = (l-1) l (l+1) (l+2)/6 (137 and 686 at l = 2, 3, section 3)
rate = {l: (l - 1) * l * (l + 1) * (l + 2) / 6 / (2 * mass(l)) for l in LS}

print("| l | points | negative | min abs(D) | at gamma | k/(2|M|) | convergence gammas | all variants negative | capped |")
print("|---|---|---|---|---|---|---|---|---|")
for l in LS:
    s = scans[l]
    g0, d0 = min(s, key=lambda p: abs(p[1]))
    neg = all(sg == "np.float64(-1.0)" for _, _, sg, _, _ in convs[l])
    print(f"| {l} | {len(s)} | {sum(d < 0 for _, d in s)} | {abs(d0):.3e} | {g0:.0f} | {rate[l]:.0f} | "
          f"{len(convs[l])} | {neg} | {sum(c[4] for c in convs[l])} |")

svg = Svg(780)
svg.text(16, 26, "No growing polar mode of the phantom star at l = 4, 5, 6", size=15, fill=INK, weight="bold")
top = caption(svg, 16, 50, "Mode condition D(γ) along the real growth-rate axis (#4)",
              "A growing mode needs D = 0. All 300 scan points are negative, so the plot shows −D. "
              "Lines: 100-point scans. Markers: base variant of the convergence runs.")
x, y, w, h = 84, top + 10, 600, 360
xm, ym = axes(svg, x, y, w, h, (1e-3, 2e4, "log"), (1e-15, 1, "log"),
              "growth rate γ (log)", "−D (log)", xfmt=fmt_dec, yfmt=fmt_dec)
for i, l in enumerate(LS):
    c = SERIES[i]
    pts = " ".join(f"{xm(g):.1f},{ym(-d):.1f}" for g, d in scans[l])
    svg.add(f'<polyline points="{pts}" fill="none" stroke="{c}" stroke-width="2"/>')
    rx = xm(rate[l])
    svg.add(f'<line x1="{rx:.1f}" y1="{y}" x2="{rx:.1f}" y2="{y + h}" stroke="{c}" '
            f'stroke-width="1" stroke-dasharray="4 3"/>')
    svg.text(rx - 3, y + 12 + 13 * i, f"point-mass rate, l = {l}", size=10, fill=INK_SOFT, anchor="end")
    for g, d, sg, spread, capped in convs[l]:
        marker(svg, "circle", xm(g), ym(-d), c,
               f"l = {l}, γ = {g:g}: D_base = {d:.6e}, relative spread {spread:.1e}"
               + (", rtol 1e-12 variant not completed" if capped else ""))
for l, dx, dy in ((5, -250, -6), (6, -250, 34)):
    g, d = next((g, d) for g, d, _, _, cap in convs[l] if cap)
    callout(svg, xm(g), ym(-d), f"l = {l}, γ = {g:g}: rtol 1e-12 variant capped at 5 h", dx=dx, dy=dy)
legend(svg, x + w, top - 8, [(f"l = {l}", SERIES[i]) for i, l in enumerate(LS)])

height = footnote(svg, ["scan_growing.py, convergence.py (#4)", "raw I4-win"],
                  "surface-laptop-4, 8 CPUs, 32 GB, Windows 11",
                  ["the top γ of l = 5 and 6 rests on four of the five convergence variants"],
                  y + h + 44, "https://github.com/nmoinvaz/negative-mass")
open(os.path.join(HERE, "l456_mode_condition.svg"), "w").write(svg.finish(height))

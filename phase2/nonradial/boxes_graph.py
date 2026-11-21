"""
Graph of the complex-s boxes counted for the phantom campaign star (#5, RESULTS.md section 5.1), at true scale.

Reads the box table of i5_boxes_output.txt (collect_i5.py), writes complex_s_boxes.svg and prints the boxes as a
table.  Usage: python boxes_graph.py
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import data_graphs
from data_graphs import INK, INK_SOFT, SERIES, Svg, axes, caption, esc, footnote, legend, marker

txt = open(os.path.join(HERE, "i5_boxes_output.txt")).read()
OMEGA = float(re.search(r"omega = ([\d.]+)", txt).group(1))
row = re.compile(r"^\| (\w+?)_(\w+) \| [^|]* \| \[(\S+), (\S+)\] \| (?:abs\(Im s\) ≤ (\S+)|\[(\S+), (\S+)\]) \| "
                 r"([-+\d.]+) \| [\d.]+ \| (\S+) \| (\S+) \| ([-+\d.]+) of (\d) \| (PASS|FAIL) \|$", re.M)
boxes = {}
plants = {}
for m in row.finditer(txt.split("\n# ", 1)[0]):
    name, run, x0, x1, ys, y0, y1, wnd, mind, s0, got, exp, verdict = m.groups()
    sym = ys is not None
    key = (float(x0), float(x1), 0.0 if sym else float(y0), float(ys or y1), sym)
    l = int(run[1]) if run.startswith("l") else 2  # t1_r40 and the like are l = 2
    boxes.setdefault(key, {"name": "stored box" if name == "rep" else name, "runs": []})["runs"].append(
        (f"{name}_{run}", l, float(wnd), float(mind), verdict))
    plants[complex(s0.replace("j", "j"))] = (f"{name}_{run}", float(got), int(exp), verdict)

print("| box | Re s | Im s | l | winding | min abs(D) | planted |")
print("|---|---|---|---|---|---|---|")
for (x0, x1, y0, y1, sym), b in sorted(boxes.items(), key=lambda kv: (kv[0][2], kv[0][0])):
    for run, l, wnd, mind, verdict in b["runs"]:
        im = f"abs(Im s) ≤ {y1:g}" if sym else f"[{y0:g}, {y1:g}]"
        print(f"| {b['name']} ({run}) | [{x0:g}, {x1:g}] | {im} | {l} | {wnd:+.4f} | {mind:.2e} | {verdict} |")
assert all(abs(r[2]) < 0.05 and r[4] == "PASS" for b in boxes.values() for r in b["runs"]), "a box is not 0 and PASS"


def decade_ticks(lo, hi):
    e0, e1 = math.ceil(math.log10(lo) - 1e-9), math.floor(math.log10(hi) + 1e-9)
    return [10.0 ** e for e in range(e0, e1 + 1)]


data_graphs.nice_log_ticks = decade_ticks

svg = Svg(780)
svg.text(16, 26, "Every box counted for the phantom star has winding 0", size=15, fill=INK, weight="bold")
top = caption(svg, 16, 50, "Boxes of the growing half plane at true scale, Im s ≥ 0 (#5)",
              "Each box returned winding 0 and recovered the zero planted in it. Mirror boxes have the same count. "
              "Blank areas were not counted.")
x, y, w, h = 84, top + 10, 600, 330
XLO, YLO = 1e-4, 0.02
xm, ym = axes(svg, x, y, w, h, (XLO, 4, "log"), (YLO, 20, "log"), "Re s (log)", "Im s (log)",
              xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:g}")
C_BOTH, C_L2, C_PLANT = SERIES[0], SERIES[1], SERIES[2]
# Hand-placed labels where boxes overlap, in data coordinates
LPOS = {"n1": (1e-3, 0.06), "s2": (0.0175, 1.0), "s1": (0.055, 0.6), "t1": (0.045, 0.2), "t2": (0.045, 1.76),
        "stored box": (0.45, 0.6)}


def tip_of(b, x0, x1, y0, y1, sym):
    im = f"|Im s| ≤ {y1:g}" if sym else f"Im s ∈ [{y0:g}, {y1:g}]"
    runs = "; ".join(f"{r}: l = {l}, winding {wd:+.4f}, min|D| {md:.2e}" for r, l, wd, md, _ in b["runs"])
    return f"{b['name']}: Re s ∈ [{x0:g}, {x1:g}], {im}; {runs}"


order = sorted(boxes.items(), key=lambda kv: -(math.log(kv[0][1] / max(kv[0][0], XLO)) *
                                                math.log(kv[0][3] / max(kv[0][2], YLO))))
for (x0, x1, y0, y1, sym), b in order:
    if x1 <= XLO:  # n3 lies left of the log axis, the panel below draws it
        continue
    ls = sorted({r[1] for r in b["runs"]})
    color = C_BOTH if ls == [2, 3] else C_L2
    px0, px1, py0, py1 = xm(max(x0, XLO)), xm(x1), ym(max(y0, YLO)), ym(y1)
    tip = tip_of(b, x0, x1, y0, y1, sym)
    if b["name"] == "stored box":
        tip += "; l = 3 was counted in section 5 and recounted with the current rule in section 5.2 (#54): 0"
    svg.add(f'<rect x="{px0:.1f}" y="{py1:.1f}" width="{px1 - px0:.1f}" height="{py0 - py1:.1f}" '
            f'fill="{color}" fill-opacity="0.16" stroke="{color}" stroke-width="1.5"><title>{esc(tip)}</title></rect>')
    lx, ly = LPOS.get(b["name"], (math.sqrt(max(x0, XLO) * x1), math.sqrt(max(y0, YLO) * y1)))
    svg.text(xm(lx), ym(ly) + 4, b["name"], size=11, fill=INK, anchor="middle", weight="bold")
    if b["name"] == "stored box":
        svg.text(xm(lx), ym(ly) + 18, "l = 3: 0, recounted in §5.2", size=10, anchor="middle")
for v, lab in ((1 - OMEGA, "1 − ω"), (1 + OMEGA, "1 + ω")):
    py = ym(v)
    svg.add(f'<line x1="{x}" y1="{py:.1f}" x2="{x + w}" y2="{py:.1f}" stroke="{INK_SOFT}" '
            f'stroke-width="1" stroke-dasharray="2 3"/>')
    svg.text(x + w - 4, py - 4, f"Im s = {lab}", size=10, anchor="end")
for s0, (run, got, exp, verdict) in plants.items():
    if s0.real > XLO:
        marker(svg, "diamond", xm(s0.real), ym(s0.imag), C_PLANT,
               f"planted at {s0.real:g} + {s0.imag:g}i ({run}): returned {got:+.0f} of {exp}, {verdict}")
svg.text(xm(1e-3), ym(3), "not counted", size=11, fill=INK_SOFT, anchor="middle")
legend(svg, x + w, top - 8, [("l = 2 and 3", C_BOTH), ("l = 2 only", C_L2), ("planted zero", C_PLANT)])

# Near the axis: n3 crosses Re s = 0, so it needs a linear Re s
by0 = caption(svg, 16, y + h + 70, "Near the imaginary axis (linear Re s)",
              "n3 crosses the axis into the damped side. n1 starts at Re s = 1e-4 and continues to 0.1.")
ix, iy, iw, ih = 84, by0 + 4, 600, 120
XA, XB = -0.0012, 0.0008
im_, iym = axes(svg, ix, iy, iw, ih, (XA, XB), (0, 0.25), "Re s", "Im s",
                xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:g}")
for (x0, x1, y0, y1, sym), b in boxes.items():
    if x0 < XB:
        ls = sorted({r[1] for r in b["runs"]})
        color = C_BOTH if ls == [2, 3] else C_L2
        px0, px1 = im_(max(x0, XA)), im_(min(x1, XB))
        svg.add(f'<rect x="{px0:.1f}" y="{iym(y1):.1f}" width="{px1 - px0:.1f}" height="{iym(y0) - iym(y1):.1f}" '
                f'fill="{color}" fill-opacity="0.16" stroke="{color}" stroke-width="1.5">'
                f'<title>{esc(tip_of(b, x0, x1, y0, y1, sym))}</title></rect>')
        svg.text((px0 + px1) / 2, iym(0.2) + 4, b["name"], size=11, fill=INK, anchor="middle", weight="bold")
svg.line(im_(0), iy, im_(0), iy + ih, INK)
for s0, (run, got, exp, verdict) in plants.items():
    if XA < s0.real < XB:
        marker(svg, "diamond", im_(s0.real), iym(s0.imag), C_PLANT,
               f"planted at {s0.real:g} + {s0.imag:g}i ({run}): returned {got:+.0f} of {exp}, {verdict}")

height = footnote(svg, ["winding2.py (#5), table by collect_i5.py", "raw I5-rpi5"], "rpi5-ubuntu, Raspberry Pi 5",
                  ["a winding counts zeros minus poles of D; section 5.3 (#121) finds no Zerilli-column pole in "
                   "any box but n3, not counted (#123)"],
                  iy + ih + 44, "https://github.com/nmoinvaz/negative-mass")
open(os.path.join(HERE, "complex_s_boxes.svg"), "w").write(svg.finish(height))

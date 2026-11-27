#!/usr/bin/env python3
"""The coefficients of W^2 a^2 = s + Omega0 + Omega1 eps + Omega2 eps^2 + ... (s = 2(D - 3), eps = 1/s), and the error of
each truncation against #195's quad-double W^2 at 64 values of D.

    python3 omega_graph.py            (writes omega_coeffs.svg beside this file, prints the numbers)

Omega0, Omega1, Omega2 are the closed forms (README.md, omega2.md); Omega3, Omega4 are #195's fit (w2_output.txt), as are
the W^2 values. Standard library only: decimal at 70 digits for the differences, which fall to about 1e-17 of W^2.
"""
import math
import os
import re
import sys
from decimal import Decimal, getcontext

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "nonradial"))
from data_graphs import (INK, SERIES, Svg, axes, better_arrow, callout, caption, footnote,  # noqa: E402
                         legend, marker)

getcontext().prec = 70
D_ = Decimal


def pi_dec():
    """pi by the decimal module's documented recipe"""
    getcontext().prec += 2
    three = D_(3)
    lasts, t, s, n, na, d, da = 0, three, 3, 1, 0, 0, 24
    while s != lasts:
        lasts = s
        n, na = n + na, na + 8
        d, da = d + da, da + 32
        t = (t * n) / d
        s += t
    getcontext().prec -= 2
    return +s


def zeta3_dec():
    """zeta(3) = (5/2) sum_{n>=1} (-1)^{n+1} / (n^3 binom(2n, n))"""
    s, c = D_(0), 1
    for n in range(1, 260):
        c = c * (2 * n) * (2 * n - 1) // (n * n)
        s += D_((-1) ** (n + 1)) / (D_(n) ** 3 * c)
    return D_(5) / 2 * s


PI, Q, Z3 = pi_dec(), D_(2).ln(), zeta3_dec()
out = open(os.path.join(HERE, "w2_output.txt")).read()
OM = [-4 * Q, 2 * PI**2 / 3 + 8 * Q**2, 80 * Z3 - D_(8) / 3 * PI**2 * Q - D_(32) / 3 * Q**3,
      D_(re.search(r"Omega3 = (-?[0-9.]+)", out).group(1)), D_(re.search(r"Omega4 = (-?[0-9.]+)", out).group(1))]
FORM = ["−4 ln 2", "2π²/3 + 8 ln²2", "80ζ(3) − (8/3)π² ln 2 − (32/3) ln³2", "fit (#195)", "fit (#195)"]
CLOSED = [True, True, True, False, False]
ROWS = [(int(m.group(1)), D_(m.group(2)))
        for m in re.finditer(r"^\s+D=(\d+)\s+([0-9.]+)\s+digits qd/od", out, re.M)]
assert len(ROWS) == 64, len(ROWS)


def partial(D, n):
    s = D_(2 * (D - 3))
    return s + OM[0] + sum(OM[k] / s**k for k in range(1, n + 1))


ERR = {n: [(D, abs(W - partial(D, n))) for D, W in ROWS] for n in range(5)}

# ---- stdout: the same numbers ----
print(f"{'k':>2}  {'Omega_k':>28}  form")
for k, o in enumerate(OM):
    print(f"{k:>2}  {float(o):>28.15g}  {FORM[k]}")
print("ratios Omega_{k+1}/Omega_k:", ", ".join(f"{float(OM[k + 1] / OM[k]):.3f}" for k in range(4)))
print(f"\n{'D':>6}  " + "  ".join(f"{'log10 err n=' + str(n):>16}" for n in range(5)))
for i in range(0, 64, 9):
    print(f"{ROWS[i][0]:>6}  " + "  ".join(f"{math.log10(ERR[n][i][1]):>16.3f}" for n in range(5)))
s_last = D_(2 * (ROWS[-1][0] - 3))
o5 = (ROWS[-1][1] - partial(ROWS[-1][0], 4)) * s_last**5
print(f"\nremainder after Omega4, times s^5, at D = {ROWS[-1][0]}: {float(o5):.6g} (an estimate of Omega5)")

# ---- the chart ----
W = 780
svg = Svg(W)
svg.text(16, 26, "The large-D coefficients of the wormhole's growth rate, W²a² = s + Ω₀ + Ω₁/s + Ω₂/s² + …,"
         " s = 2(D − 3)", size=14, fill=INK, weight="bold")

# panel (a): |Omega_k| against k
px, pw, ph = 96, 600, 230
py = caption(svg, 16, 54, "(a) The coefficients",
             "|Ωₖ| on a log scale; ● positive, ▼ negative. Consecutive ratios Ω₁/Ω₀ … Ω₄/Ω₃: "
             + ", ".join(f"{float(OM[k + 1] / OM[k]):.2f}" for k in range(4)) + ".")
xm, ym = axes(svg, px, py, pw, ph, (-0.4, 4.6), (0, 4.6), "order k (the coefficient of 1/sᵏ)",
              "|Ωₖ|", xfmt=lambda v: str(int(v)),
              yfmt=lambda v: {0: "1", 1: "10", 2: "100", 3: "1000", 4: "10⁴"}.get(int(v), ""))
cols = {True: SERIES[0], False: SERIES[1]}
pts = []
for k, o in enumerate(OM):
    x, y = xm(k), ym(math.log10(abs(float(o))))
    pts.append((x, y))
    marker(svg, "circle" if o > 0 else "triangle-down", x, y, cols[CLOSED[k]],
           f"Ω{k} = {float(o):.15g} ({FORM[k]})")
lab = [f"Ω₀ = −4 ln 2 = {float(OM[0]):.6f}", f"Ω₁ = 2π²/3 + 8 ln²2 = {float(OM[1]):.6f}",
       f"Ω₂ = 80ζ(3) − (8/3)π² ln 2 − (32/3) ln³2 = {float(OM[2]):.6f}", f"Ω₃ = {float(OM[3]):.6f} (fit)",
       f"Ω₄ = {float(OM[4]):.1f} (fit)"]
for k in range(4):
    svg.text(pts[k][0] + 16, pts[k][1] + 4, lab[k], size=11, fill=INK)
svg.text(pts[4][0] - 12, pts[4][1] + 4, lab[4], size=11, fill=INK, anchor="end")
callout(svg, *pts[2], "derived in #26", dx=-40, dy=-36)
callout(svg, *pts[3], "#208 computes it from G₅ = 0", dx=-40, dy=-34)
legend(svg, px + pw, py - 12, [("closed form", SERIES[0]), ("fit to W² (#195)", SERIES[1])])

# panel (b): error of each truncation against W^2
py2 = caption(svg, 16, py + ph + 62, "(b) What each coefficient buys",
              "|W²a² − the series through Ωₙ| against #195's 60-digit W² at 64 values of D: each term adds about "
              "one power of 1/D.")
ph2 = 300
xm2, ym2 = axes(svg, px, py2, pw, ph2, (50, 10000, "log"), (-18, 0), "dimension D",
                "error of the truncated series", xfmt=lambda v: f"{v:g}",
                yfmt=lambda v: "1" if v == 0 else f"1e{int(v)}")
order = [(4, SERIES[3]), (3, SERIES[2]), (2, SERIES[4]), (1, SERIES[5]), (0, SERIES[-1])]
name = {0: "s + Ω₀", 1: "+ Ω₁/s", 2: "+ Ω₂/s²", 3: "+ Ω₃/s³", 4: "+ Ω₄/s⁴"}
for n, col in order:
    xy = [(xm2(D), ym2(math.log10(e))) for D, e in ERR[n]]
    path = " ".join(f"{x:.1f},{y:.1f}" for x, y in xy)
    svg.add(f'<polyline points="{path}" fill="none" stroke="{col}" stroke-width="2"/>')
    for i in list(range(0, 64, 9)) + [63]:
        D, e = ERR[n][i]
        marker(svg, "circle", xy[i][0], xy[i][1], col, f"D = {D}, through Ω{n}: {float(e):.3e}")
    svg.text(xy[-1][0] + 16, xy[-1][1] + 4, name[n], size=11, fill=INK)
better_arrow(svg, px + 36, ym2(-9), px + 36, ym2(-15), label="smaller")
xe, ye = xm2(ROWS[-1][0]), ym2(math.log10(ERR[4][-1][1]))
callout(svg, xe, ye, "≈ Ω₅/s⁵, with Ω₅ ≈ 4.1 × 10⁴ (estimate)", dx=-150, dy=-34)

h = footnote(svg, ["closed forms: README.md, omega2.md (#26)",
                   "Ω₃, Ω₄ and W² (quad-double, dh610-fedora): w2_output.txt (#195)"],
             "", [], py2 + ph2 + 44, "https://github.com/nmoinvaz/negative-mass")
path = os.path.join(HERE, "omega_coeffs.svg")
open(path, "w").write(svg.finish(h))
print(f"\nwrote {path}")

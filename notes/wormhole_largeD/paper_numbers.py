"""Every number quoted in Paper III (paper/3-wormhole-largeD), recomputed from the repo's data and the published tables.

    python3 paper_numbers.py            (prints; exit 1 if a stated check fails)

Inputs: growth.csv (profiles.py: Torii & Shinkai's Table I for D <= 11, modes_r.py for D >= 20); logs_order3.txt
(Omega1); Shinkai & Torii 2017 Table I (alpha_GB = 0 rows, transcribed below); Shinkai & Hayward 2002's lifetime
exponents 0.60 (ghost pulses) and 0.61 (normal pulses).
"""
import csv
import math
import os
import re
import sys

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
L2 = math.log(2)
OM0 = -4 * L2
OM1_CLOSED = 2 * math.pi**2 / 3 + 8 * L2**2
bad = 0


def check(cond, what):
    global bad
    if not cond:
        bad += 1
        print("  CHECK FAILED:", what)


# Omega1 from the matching integral (finest grid in logs_order3.txt)
om1 = float(re.findall(r"Omega1 = ([0-9.]+)", open(os.path.join(HERE, "logs_order3.txt")).read())[-1])
print(f"Omega1 (computed) = {om1:.10f};  2 pi^2/3 + 8 ln^2 2 = {OM1_CLOSED:.11f};  diff {om1 - OM1_CLOSED:.1e}")
check(abs(om1 - OM1_CLOSED) < 1e-8, "Omega1 closed form")
mp.mp.dps = 16
rel = mp.pslq([mp.mpf(repr(om1)), mp.pi**2, mp.log(2)**2, mp.log(2), 1], tol=mp.mpf("1e-9"), maxcoeff=60,
              maxsteps=10**6)
print("PSLQ relation [Omega1, pi^2, ln^2 2, ln 2, 1]:", rel)
check(rel == [3, -2, -24, 0, 0], "PSLQ relation")
mp.mp.dps = 30
# J = int_1^inf ln^2 R (R-1)^(-3/2) dR; with R = 1 + u^2 and one integration by parts, J = 8 int_0^inf ln(1+u^2)/(1+u^2) du
J = mp.quad(lambda R: mp.log(R)**2 / (R - 1)**mp.mpf(1.5), [1, 1.5, 2, 5, 10, 100, 1e4, mp.inf])
Ju = 8 * mp.quad(lambda u: mp.log(1 + u * u) / (1 + u * u), [0, 1, 10, mp.inf])
print(f"J = {mp.nstr(J, 15)} (direct), {mp.nstr(Ju, 15)} (u form);  8 pi ln 2 = {mp.nstr(8 * mp.pi * mp.log(2), 15)}")
check(abs(Ju - 8 * mp.pi * mp.log(2)) < 1e-12 and abs(J - Ju) < 1e-8, "J = 8 pi ln 2")

rows = list(csv.DictReader(open(os.path.join(HERE, "growth.csv"))))
print("\n  D    W^2 a^2 (data)      3-term formula   rel. diff    (data - 3 terms) s^2   source")
for r in rows:
    D = int(r["D"])
    w = float(r["W2a2"])
    s = 2 * (D - 3)
    f3 = s + OM0 + OM1_CLOSED / s
    print(f"{D:5d}  {w:18.10f}  {f3:15.8f}  {abs(w - f3) / w:9.2e}   {(w - f3) * s * s:12.2f}       {r['source']}")
    if 7 <= D <= 11:
        check(abs(w - f3) / w < 0.011, f"3-term within 1.1% at D = {D}")
fit = {}
for r in rows:
    D = int(r["D"])
    if D in (200, 500):
        s = 2 * (D - 3)
        fit[D] = (float(r["W2a2"]) - (s + OM0 + OM1_CLOSED / s)) * s * s
check(all(70 < v < 76 for v in fit.values()), "remainder ~ 73 / s^2 at D = 200, 500")
print("  (the r-integration's own uncertainty, logs_modes_r.txt, times s^2: 0.02 at D = 200, 0.6 at 500, 6 at 1000,"
      " 80 at 2000; so the last column means nothing beyond D = 1000)")

# e-folding times and the type I lifetime slope 1/(sqrt 2 W) in the null coordinate x^- = (t - r_*)/sqrt 2
TS = {int(r["D"]): float(r["W2a2"]) for r in rows}
print("\n  D   1/W (e-folding, units of a)   predicted slope 1/(sqrt2 W)")
for D in (4, 5, 6, 100):
    W = math.sqrt(TS[D])
    print(f"{D:4d}   {1 / W:.4f}                        {1 / (math.sqrt(2) * W):.4f}")
# Shinkai & Torii 2017, Table I, alpha_GB = 0 rows: (D, Delta E / a0, x_H^- / a0), two amplitudes per D
ST17 = {4: [(0.03, 2.94), (0.10, 2.09)], 5: [(0.15, 2.26), (0.61, 1.66)], 6: [(0.38, 1.92), (1.50, 1.46)]}
print("\n  D   ST17 two-run slope  (+/- from the table's rounding: dE +/- 0.005, x +/- 0.005)   predicted")
for D, ((e1, x1), (e2, x2)) in ST17.items():
    lr = math.log(e2 / e1)
    sl = (x1 - x2) / lr
    dlr = 0.005 / e1 + 0.005 / e2
    err = 0.01 / lr + sl * dlr / lr
    pred = 1 / (math.sqrt(2) * math.sqrt(TS[D]))
    print(f"{D:4d}   {sl:.3f} +/- {err:.3f}                                                       {pred:.3f}"
          f"   ({(sl - pred) / err:+.1f} sigma)")
    check(abs(sl - pred) < 2 * err, f"ST17 slope within 2 sigma at D = {D}")
print("\n  D = 4, Shinkai & Hayward 2002 small-amplitude runs: 0.60 (ghost pulses), 0.61 (normal pulses);"
      f" predicted {1 / (math.sqrt(2) * math.sqrt(TS[4])):.3f}")
print("\nchecks failed:", bad)
sys.exit(1 if bad else 0)

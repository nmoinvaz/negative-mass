#!/usr/bin/env python3
"""Cost of one run and of one bisection of the planned dual-null code (#39 design), from measured kernel throughput.

    python3 cost_model.py

Inputs: the per-cell cost of a stand-in for the cell update, best of three on N = 2000 (4e6 cells), one thread, nice 10,
dh610-fedora (i7-12700, load average about 9 of 20): negative-mass-raw I39-dh610/bench/run01_kernel_bench_{python,rust}.log.
The run length model follows [SH02]'s lifetime law, x_H^- = x_c + (1/(sqrt2 W)) ln(1/|p - p*|) in null time, with the pulse
reaching the throat at x+ = x- = x_c = 3 and a margin of 3 a after the departure to classify the fate. Time t = (x+ + x-)/sqrt2.
"""
import math

NS = {"numba": 61.7, "rust": 79.6, "numpy": 192.3}                      # ns per cell, from the logs above
W2 = {4: 1.39705243433350684, 5: 2.98495893035176025, 6: 4.68662054311946479}   # W^2 a^2, #195 (qd column)
DELTA = 1e-10
XC, MARGIN = 3.0, 3.0
RUNS_PER_BISECTION = 45      # about 33 halvings from a bracket of order 0.1, plus lifetime points at fixed |p - p*|
H4 = [0.02, 0.01, 0.005, 0.0025]

print("per-cell cost (ns): " + ", ".join(f"{k} {v}" for k, v in NS.items()))
print(f"|p - p*| = {DELTA:g}, x_c = {XC}, margin = {MARGIN} a, {RUNS_PER_BISECTION} runs per bisection (no early exit)\n")
print(f"{'D':>2} {'1/W':>7} {'T_run':>7} {'X':>6} {'h':>8} {'cells':>10} {'numba s':>9} {'rust s':>8} {'numpy s':>9}"
      f" {'bisect numba':>13} {'bisect numpy':>13}")
tot = {k: 0.0 for k in NS}
for D, w2 in W2.items():
    W = math.sqrt(w2)
    T = math.sqrt(2) * XC + math.log(1 / DELTA) / W + MARGIN
    X = T / math.sqrt(2)
    for h4 in H4:
        h = h4 / math.sqrt(D - 3)
        cells = (X / h) ** 2
        t = {k: cells * v * 1e-9 for k, v in NS.items()}
        print(f"{D:>2} {1 / W:7.3f} {T:7.1f} {X:6.1f} {h:8.5f} {cells:10.3g} {t['numba']:9.2f} {t['rust']:8.2f} {t['numpy']:9.2f}"
              f" {RUNS_PER_BISECTION * t['numba'] / 60:11.1f} m {RUNS_PER_BISECTION * t['numpy'] / 60:11.1f} m")
        for k in NS:
            tot[k] += RUNS_PER_BISECTION * t[k]
print("\ncampaign of one bisection per (D, h) in the table, 12 bisections, ghost family: "
      + ", ".join(f"{k} {v / 3600:.2f} h" for k, v in tot.items()) + " of one thread")
print("the same for the normal family: double it; at 4 threads (independent bisections): divide by 4")

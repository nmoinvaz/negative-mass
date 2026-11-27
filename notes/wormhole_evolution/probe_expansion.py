#!/usr/bin/env python3
"""How far the expanding branch can be followed in this gauge (#213's question; a probe, not a check).

    python3 probe_expansion.py

The F1 expansion case (one-sided ghost pulse A = +0.04, derivative form, x_c = 3, w = 1/sqrt3) at D = 4 and D = 5, h = 0.02 and 0.01,
on a long rectangle (X = 30, t up to 42). Prints the diagonal throat r, H = d ln r/d tau-like rate, the maximum of e^{-f} over live cells,
the first dead cell, and the e-folds of the throat radius reached. Nothing is gated; the numbers are read as measurements.
"""
import math
import sys
import time

import numpy as np

import dn

SQ2 = math.sqrt(2.0)
print(f"probe_expansion.py: numba {'present' if dn.HAVE_NUMBA else 'absent'}, {time.strftime('%Y-%m-%d %H:%M:%S')}")
for n, h in ((2, 0.02), (2, 0.01), (3, 0.0141421356), (3, 0.00707106781)):
    X = 30.0
    N = int(round(X / h))
    t0 = time.time()
    d = dn.ellis_rays(n, N, h, pulses=[dict(field=0, ray="+", A=0.04, xc=3.0, w=1 / math.sqrt(3.0), kind="dphi")])
    res = dn.evolve(d, kernel="numba", excise=True, rcut=0.05, fcut=25.0, rexp=1e9)       # rexp huge: run the whole rectangle
    dg = res.diag
    ok = np.isfinite(dg["r"])
    j_last = np.where(ok)[0][-1]
    emf = np.exp(-res.G[1])
    live = ~res.dead.astype(bool)
    emax = np.nanmax(np.where(live, emf, np.nan))
    ic = np.unravel_index(np.nanargmax(np.where(live, emf, np.nan)), emf.shape)
    ndead = int(res.dead.sum())
    print(f"D={n + 2} h={h:.5f} N={N} ({time.time() - t0:.1f} s): fate by events {res.events or 'none'}; last live diagonal cell j = {j_last} of {N - 1} "
          f"(t = {dg['t'][j_last]:.2f}); dead cells {ndead}; max e^-f over live cells {emax:.3f} at x+ = {ic[0] * h:.2f}, x- = {ic[1] * h:.2f}")
    for t in (8, 12, 16, 20, 24, 28, 32, 36, 40):
        j = int(np.argmin(np.abs(dg["t"] - t)))
        if j <= j_last and abs(dg["t"][j] - t) < 0.5:
            r, H = dg["r"][j], dg["H"][j]
            print(f"    t = {dg['t'][j]:6.2f}  tau = {dg['tau'][j]:6.2f}  r_diag = {r:10.4f}  e-folds ln r = {math.log(r):6.3f}  H = {H:+8.4f}  (1/H = {1 / H if H else float('inf'):7.3f})")

"""Exact growth rates What = W_th a of the general ghost wormhole over a list of (D, gamma1, eps2), with a convergence
check (larger end radius). Real growth rates only (a sign change of the real determinant).

    python3 scan_general.py OUT.csv "D,gamma1,eps2,lo,hi" ...
"""
import csv
import sys
import time

import numpy as np
from scipy.optimize import brentq

import modes_general as mg


def root(w, lo, hi, N, **kw):
    f = lambda h: w.det(w.W2_of(h), **kw)[0]  # noqa: E731
    grid = np.linspace(lo, hi, N)
    vals = [f(h) for h in grid]
    for (a, fa), (b, fb) in zip(zip(grid, vals), zip(grid[1:], vals[1:])):
        if np.isfinite(fa) and np.isfinite(fb) and np.sign(fa) != np.sign(fb):
            r = brentq(f, a, b, xtol=1e-10)
            if abs(f(r)) < 1e-6:        # else a jump of the amp form's signs, not a zero: look further
                return r
    return float("nan")


if __name__ == "__main__":
    out = sys.argv[1]
    with open(out, "a", newline="") as fh:
        wr = csv.writer(fh)
        for spec in sys.argv[2:]:
            D, g1, e2, lo, hi = spec.split(",")
            D, g1, e2, lo, hi = int(D), float(g1), float(e2), float(lo), float(hi)
            t = time.time()
            w = mg.Wormhole(D, g1, e2)
            r1 = root(w, lo, hi, 9)
            r2 = root(w, lo, hi, 9, rmax_factor=90.0, X=45.0) if np.isfinite(r1) else float("nan")
            wr.writerow([D, g1, e2, f"{r1:.10f}", f"{r2:.10f}", f"{time.time() - t:.0f}"])
            fh.flush()
            print(f"D = {D}, gamma1 = {g1}, eps2 = {e2}: What = {r1:.10f} (larger end radius {r2:.10f}),"
                  f" {time.time() - t:.0f} s", flush=True)

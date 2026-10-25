"""All real growth rates What = W_th a of the charged (or any) ghost wormhole in a range, for each (D, gamma1, eps2):
the determinant's sign on a grid of step h, then each sign change refined by brentq. Coarse brackets (scan_general.py
takes the first sign change of a 9-point grid) miss pairs of roots: for the symmetric charged wormhole the even
sector has two real roots close together near the charge where they merge into a complex pair (D 6, eps2 0.05:
between What 1.15 and 1.2 and between 1.3 and 1.35).

For each root the solution whose sign index k changes across it is printed: 0 (the even throat solution) or 3 (the
odd one) for a symmetric wormhole; for an asymmetric one both may change.

    python3 charged_scan.py OUT.csv lo hi h "D,gamma1,eps2" ...      (CSV rows: D, gamma1, eps2, What, sector, seconds)
"""
import csv
import sys
import time

import numpy as np
from scipy.optimize import brentq

import modes_general as mg


def main():
    out, lo, hi, h = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    with open(out, "a", newline="") as fh:
        wr = csv.writer(fh)
        for spec in sys.argv[5:]:
            D, g1, e2 = spec.split(",")
            D, g1, e2 = int(D), float(g1), float(e2)
            t = time.time()
            w = mg.Wormhole(D, g1, e2)
            grid = np.arange(lo, hi + 0.5 * h, h)
            vals = [w.det(w.W2_of(x)) for x in grid]
            roots = []
            for (a, (fa, sa)), (b, (fb, sb)) in zip(zip(grid, vals), zip(grid[1:], vals[1:])):
                if np.sign(fa) != np.sign(fb):
                    r = brentq(lambda x: w.det(w.W2_of(x))[0], a, b, xtol=1e-10)
                    moved = [r0 for r0 in (0, 3) if sa["k"][(1, r0)] != sb["k"][(1, r0)]]
                    roots.append((r, "+".join(str(m) for m in moved) or "-"))
            for r, sec in roots:
                wr.writerow([D, g1, e2, f"{r:.10f}", sec, f"{time.time() - t:.0f}"])
            if not roots:
                wr.writerow([D, g1, e2, "nan", "-", f"{time.time() - t:.0f}"])
            fh.flush()
            print(f"D = {D}, gamma1 = {g1}, eps2 = {e2}: real roots in [{lo}, {hi}] " +
                  (", ".join(f"{r:.10f} (sector {s})" for r, s in roots) or "none") + f", {time.time() - t:.0f} s",
                  flush=True)


if __name__ == "__main__":
    main()

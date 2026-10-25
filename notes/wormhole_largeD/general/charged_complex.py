"""The charged symmetric wormhole's even pair of growth rates, followed in the charge at fixed D with modes_complex.py:
two real roots (charged_scan.py) that move toward each other as eps2 grows, meet at eps2*, and continue as a complex
pair W_r +- i W_i (a growth that oscillates). Near the meeting the mismatch is f ~ a (W - m)^2 + b (eps2 - eps2*), so
(W_2 - W_1)^2 is linear in eps2 there: eps2* is where that line reaches 0, and past it the guess is m +- i sqrt(...).

    python3 charged_complex.py OUT.csv D W1 W2 eps2_0 eps2_1 ... (W1, W2: the real pair at eps2_0, from charged_scan.py;
    the continuation solves eps2_1, eps2_2, ...)
    (CSV rows: D, gamma1, eps2, Re What, Im What, seconds; one row per root, the complex ones with Im > 0)
"""
import csv
import sys
import time

import numpy as np

import modes_complex as mc


def main():
    out, D, w1, w2 = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    eps = [float(v) for v in sys.argv[5:]]
    hist = [(eps[0], w1, w2)]         # the given pair at eps[0]; the continuation starts at eps[1]
    root = None                       # the complex root (Im > 0) once the pair has met
    with open(out, "a", newline="") as fh:
        wr = csv.writer(fh)
        for e2 in eps[1:]:
            t = time.time()
            w = mc.ComplexWormhole(D, 0.0, e2)
            if root is None:
                if len(hist) >= 2:
                    (e0, a0, b0), (e1, a1, b1) = hist[-2], hist[-1]
                    d0, d1 = (b0 - a0) ** 2, (b1 - a1) ** 2
                    m = 0.5 * (a1 + b1) + (0.5 * (a1 + b1) - 0.5 * (a0 + b0)) * (e2 - e1) / (e1 - e0)
                    d = d1 + (d1 - d0) * (e2 - e1) / (e1 - e0)
                    guesses = [m - 0.5 * np.sqrt(d + 0j), m + 0.5 * np.sqrt(d + 0j)]
                else:
                    guesses = [w1, w2]
                found = []
                for g in guesses:
                    try:
                        found.append(w.root_c(complex(g) if abs(complex(g).imag) > 1e-9 else complex(g.real, 0.0)))
                    except Exception as exc:               # findroot gave up
                        print(f"  eps2 {e2}: no root from {g:.6f} ({type(exc).__name__})", flush=True)
                found = [complex(f.real, abs(f.imag)) for f in found]
                if len(found) == 2 and max(abs(f.imag) for f in found) < 1e-8 and abs(found[0] - found[1]) > 1e-6:
                    a, b = sorted(f.real for f in found)
                    hist.append((e2, a, b))
                    rows = [(a, 0.0), (b, 0.0)]
                else:
                    cx = [f for f in found if f.imag > 1e-8]
                    if not cx:
                        print(f"  eps2 {e2}: lost the pair: {found}", flush=True)
                        break
                    root = cx[0]
                    rows = [(root.real, root.imag)]
            else:
                root = w.root_c(root)
                root = complex(root.real, abs(root.imag))
                rows = [(root.real, root.imag)]
            for re_, im_ in rows:
                wr.writerow([D, 0.0, e2, f"{re_:.10f}", f"{im_:.10f}", f"{time.time() - t:.0f}"])
            fh.flush()
            print(f"D = {D}, eps2 = {e2}: " + ", ".join(f"{r:.10f}{i:+.10f}i" for r, i in rows) +
                  f" ({time.time() - t:.0f} s)", flush=True)


if __name__ == "__main__":
    main()

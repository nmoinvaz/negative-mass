"""Growth rates What = W_th a of the asymmetric wormhole (gamma1 = 1, eps2 = 0) at several D, for the large-D fit
What = (D - 3) sqrt(omega) + c0 + c1 / (D - 3) + ... (fit_largeD.py). The bracket is the leading guess
1.0293216299 (D - 3) + 0.47, widened until the determinant changes sign; with --far each root is repeated with the
end pushed out (rmax_factor 90, X 45).

    python3 largeD_roots.py OUT.csv D [D ...] [--far]      (CSV rows: D, What, What at the larger end, seconds)
"""
import csv
import sys
import time

from scipy.optimize import brentq

import modes_general as mg

SQRT_OMEGA = 1.0293216299


def root(w, guess, **kw):
    f = lambda h: w.det(w.W2_of(h), **kw)[0]  # noqa: E731
    for hw in (0.05, 0.15, 0.4):
        lo, hi = guess - hw, guess + hw
        if f(lo) * f(hi) < 0:
            return brentq(f, lo, hi, xtol=1e-13)
    return float("nan")


def main():
    out = sys.argv[1]
    far = "--far" in sys.argv
    Ds = [int(a) for a in sys.argv[2:] if not a.startswith("--")]
    with open(out, "a", newline="") as fh:
        wr = csv.writer(fh)
        for D in Ds:
            t = time.time()
            w = mg.Wormhole(D, 1.0, 0.0)
            guess = SQRT_OMEGA * (D - 3) + 0.47
            r = root(w, guess)
            r2 = root(w, r if r == r else guess, rmax_factor=90.0, X=45.0) if far else float("nan")
            wr.writerow([D, f"{r:.13f}", f"{r2:.13f}", f"{time.time() - t:.0f}"])
            fh.flush()
            print(f"D {D}: What {r:.13f}, at the larger end {r2:.13f}, offset {r - SQRT_OMEGA * (D - 3):.12f}, "
                  f"{time.time() - t:.0f} s", flush=True)


if __name__ == "__main__":
    main()

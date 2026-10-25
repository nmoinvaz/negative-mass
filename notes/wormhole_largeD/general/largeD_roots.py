"""Growth rates What = W_th a of the asymmetric wormhole (eps2 = 0; gamma1 = 1 unless --g1) at several D, for the
large-D fit What = (D - 3) sqrt(omega) + c0 + c1 / (D - 2) + ... (fit_largeD.py). The bracket is centred on
sqrt(omega) (D - 3) + c, with sqrt(omega) from leading.py (--slope; 1.0293216299 for gamma1 = 1) and c the offset of
the previous D (0.47 at first), widened until the determinant changes sign; with --far each root is repeated with
the end pushed out (rmax_factor 90, X 45). Give the D in increasing order.

    python3 largeD_roots.py OUT.csv D [D ...] [--far] [--g1=G --slope=S]
    (CSV rows: D, What, What at the larger end, seconds)
"""
import csv
import sys
import time

from scipy.optimize import brentq

import modes_general as mg

SQRT_OMEGA = 1.0293216299


def root(w, guess, **kw):
    f = lambda h: w.det(w.W2_of(h), **kw)[0]  # noqa: E731
    for hw in (0.05, 0.15, 0.4, 1.0):
        lo, hi = guess - hw, guess + hw
        if f(lo) * f(hi) < 0:
            return brentq(f, lo, hi, xtol=1e-13)
    return float("nan")


def main():
    out = sys.argv[1]
    far = "--far" in sys.argv
    opt = dict(a[2:].split("=") for a in sys.argv[2:] if a.startswith("--") and "=" in a)
    g1, slope = float(opt.get("g1", 1.0)), float(opt.get("slope", SQRT_OMEGA))
    Ds = [int(a) for a in sys.argv[2:] if not a.startswith("--")]
    c = 0.47
    with open(out, "a", newline="") as fh:
        wr = csv.writer(fh)
        for D in Ds:
            t = time.time()
            w = mg.Wormhole(D, g1, 0.0)
            guess = slope * (D - 3) + c
            r = root(w, guess)
            r2 = root(w, r if r == r else guess, rmax_factor=90.0, X=45.0) if far else float("nan")
            if r == r:
                c = r - slope * (D - 3)
            wr.writerow([D, f"{r:.13f}", f"{r2:.13f}", f"{time.time() - t:.0f}"])
            fh.flush()
            print(f"D {D}: What {r:.13f}, at the larger end {r2:.13f}, offset {r - slope * (D - 3):.12f}, "
                  f"{time.time() - t:.0f} s", flush=True)


if __name__ == "__main__":
    main()

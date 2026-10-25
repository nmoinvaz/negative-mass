"""Large-D fit of the asymmetric wormhole's growth rate (largeD_roots.py's CSVs): What against the leading result
(D - 3) sqrt(omega) of leading.py.

  1. With the slope fixed at sqrt(omega): What - (D - 3) sqrt(omega) = c0 + c1 / N + ... + c_{K-1} / N^{K-1}, N = D - 3,
     fitted on D in [Dmin, Dfit] and checked on the D above Dfit (held out).
  2. With the slope free: the same fit with an s N term; s must come out at sqrt(omega).
  3. Control: the same free fit to What + 1e-3 N (a wrong slope) must give s off by 1e-3, so check 2 can fail.

    python3 fit_largeD.py CSV [CSV ...] [--dmin=8] [--dfit=20] [--k=5]
"""
import csv
import sys

import numpy as np

SQRT_OMEGA = 1.0293216299


def load(paths):
    rows = {}
    for p in paths:
        with open(p) as fh:
            for r in csv.reader(fh):
                if r and r[1] != "nan":
                    rows[int(r[0])] = (float(r[1]), float(r[2]))
    return dict(sorted(rows.items()))


def basis(D, k, slope):
    """columns N (with slope), 1, 1/N, ..., 1/N^(k-1), N = D - 3"""
    N = np.asarray(D, float) - 3.0
    cols = [N ** -j for j in range(k)]
    return np.array(([N] if slope else []) + cols).T


def fit(D, y, k, slope):
    c, *_ = np.linalg.lstsq(basis(D, k, slope), y, rcond=None)
    return c


def predict(c, Dq, k, slope):
    return float(basis([Dq], k, slope)[0] @ c)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opt = dict(a[2:].split("=") for a in sys.argv[1:] if a.startswith("--") and "=" in a)
    dmin, dfit, k = int(opt.get("dmin", 8)), int(opt.get("dfit", 20)), int(opt.get("k", 5))
    rows = load(args)
    print("D, What, What at the larger end (difference), offset What - (D - 3) sqrt(omega)")
    for D, (w, w2) in rows.items():
        print(f"  {D:3d}  {w:.12f}  {w2:.12f} ({(w2 - w) / w:+.1e})  {w - SQRT_OMEGA * (D - 3):.12f}")
    Ds = np.array([D for D in rows if dmin <= D <= dfit], float)
    held = [D for D in rows if D > dfit]
    W = np.array([rows[int(D)][0] for D in Ds])
    bad = []
    off = W - SQRT_OMEGA * (Ds - 3)
    c = fit(Ds, off, k, False)
    print(f"fixed slope, {k} terms on D {dmin}..{dfit}: c = {', '.join(f'{v:.10g}' for v in c)}")
    print(f"  max fit residual {np.max(np.abs(off - [predict(c, D, k, False) for D in Ds])):.1e}")
    for D in held:
        e = rows[D][0] - SQRT_OMEGA * (D - 3) - predict(c, D, k, False)
        print(f"  held out D {D}: data - fit {e:+.2e} (relative to What {e / rows[D][0]:+.1e})")
    for kk in (k - 1, k):
        cs = fit(Ds, W, kk, True)
        print(f"free slope, {kk} more terms: s = {cs[0]:.10f} (sqrt(omega) {SQRT_OMEGA}, s - sqrt(omega) "
              f"{cs[0] - SQRT_OMEGA:+.1e}), c0 = {cs[1]:.6f}")
    cs = fit(Ds, W, k - 1, True)
    if abs(cs[0] - SQRT_OMEGA) > 1e-5:
        bad.append("the free slope is not sqrt(omega)")
    cm = fit(Ds, W + 1e-3 * (Ds - 3), k - 1, True)
    print(f"control, What + 1e-3 (D - 3): s = {cm[0]:.10f} ({cm[0] - SQRT_OMEGA:+.1e})")
    if abs(cm[0] - SQRT_OMEGA) <= 1e-5:
        bad.append("the control passed")
    print("PASS" if not bad else "FAIL: " + "; ".join(bad))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

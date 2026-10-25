"""Where the charged symmetric wormhole's two real growth rates meet (charged_complex.py's CSVs): eps2*(D).

Near the meeting the mismatch is f ~ a (W - m)^2 + b (eps2 - eps2*), so d = (W_2 - W_1)^2 is linear in eps2 there:
positive on the real side, -(2 Im W)^2 on the complex side. eps2* is the zero of d, from the quadratic through the
three points nearest the sign change (the linear fit through the two bracketing points is printed beside it; their
difference is the error of the reading). Then the large-D trend: eps2* (D - 2)^2 against 1/(D - 2).

Control: d computed with Im W doubled on the complex side (a wrong branch) must move eps2* by more than the
quadratic-linear difference.

    python3 charged_merge.py charged_complex_D4.csv charged_complex_D6.csv ...
"""
import csv
import sys
from collections import defaultdict

import numpy as np


def d_of(path, im_scale=1.0):
    rows = defaultdict(list)
    D = None
    for r in csv.reader(open(path)):
        if not r:
            continue
        D = int(r[0])
        rows[float(r[2])].append((float(r[3]), float(r[4])))
    pts = []
    for e2, ws in sorted(rows.items()):
        if len(ws) == 2 and all(abs(w[1]) < 1e-9 for w in ws):
            pts.append((e2, (ws[1][0] - ws[0][0]) ** 2))
        elif len(ws) == 1 and ws[0][1] > 0:
            pts.append((e2, -(2 * im_scale * ws[0][1]) ** 2))
    return D, pts


def merge(pts):
    i = next(i for i in range(len(pts) - 1) if pts[i][1] > 0 >= pts[i + 1][1])
    (x0, y0), (x1, y1) = pts[i], pts[i + 1]
    lin = x0 - y0 * (x1 - x0) / (y1 - y0)
    j = i - 1 if i > 0 and (i + 2 >= len(pts) or abs(pts[i - 1][0] - lin) < abs(pts[i + 2][0] - lin)) else i
    xs, ys = zip(*pts[j:j + 3])
    c = np.polyfit(xs, ys, 2)
    rts = [r.real for r in np.roots(c) if abs(r.imag) < 1e-12 and x0 <= r.real <= x1]
    quad = rts[0] if rts else float("nan")
    return quad, lin, (x0, x1)


def main():
    res = {}
    bad = []
    for p in sys.argv[1:]:
        D, pts = d_of(p)
        q, l, br = merge(pts)
        _, ptc = d_of(p, im_scale=2.0)
        qc, _, _ = merge(ptc)
        res[D] = (q, l)
        print(f"D {D:2d}: eps2* = {q:.6f} (quadratic), {l:.6f} (linear), bracket {br[0]}..{br[1]}; "
              f"control (Im doubled): {qc:.6f}, moves {qc - q:+.1e}")
        if not abs(qc - q) > abs(q - l):
            bad.append(f"D {D}: the control did not move eps2* beyond the reading error")
    Ds = np.array(sorted(res), float)
    y = np.array([res[int(D)][0] * (D - 2) ** 2 for D in Ds])
    print("eps2* (D - 2)^2: " + ", ".join(f"D {int(D)}: {v:.4f}" for D, v in zip(Ds, y)))
    if len(Ds) >= 3:
        c = np.polyfit(1 / (Ds - 2), y, 2)
        print(f"quadratic in 1/(D - 2): eps2* (D - 2)^2 -> {c[2]:.4f} as D -> infinity (coefficients {c[1]:.4f}, "
              f"{c[0]:.4f}); residual {np.max(np.abs(np.polyval(c, 1 / (Ds - 2)) - y)):.1e}")
        sl = np.polyfit(np.log(Ds - 2), np.log([res[int(D)][0] for D in Ds]), 1)[0]
        print(f"log-log slope of eps2* against D - 2: {sl:.3f}")
        # held out: the same fit without the largest D predicts it
        if len(Ds) >= 4:
            ch = np.polyfit(1 / (Ds[:-1] - 2), y[:-1], 2)
            print(f"held out D {int(Ds[-1])}: predicted {np.polyval(ch, 1 / (Ds[-1] - 2)):.4f}, measured {y[-1]:.4f}; "
                  f"limit without it {ch[2]:.4f}")
        # planted truth through the same fit and the same D: a limit of 0.9 with the fitted 1/(D - 2) and 1/(D - 2)^2
        # terms, plus an unmodelled cubic term of the size of the quadratic one (+-), and noise at 1e-4
        rng = np.random.default_rng(7)
        errs = []
        for c3 in (0.0, abs(c[0]), -abs(c[0])):
            for _ in range(20):
                u = 1 / (Ds - 2)
                yp = 0.9 + c[1] * u + c[0] * u ** 2 + c3 * u ** 3 + 1e-4 * rng.standard_normal(len(Ds))
                errs.append(np.polyfit(u, yp, 2)[2] - 0.9)
        print(f"planted limit 0.9 (cubic term 0, +-{abs(c[0]):.2f}; noise 1e-4): recovered to {min(errs):+.4f} ... "
              f"{max(errs):+.4f}")
    print("PASS" if not bad else "FAIL: " + "; ".join(bad))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

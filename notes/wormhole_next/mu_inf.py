#!/usr/bin/env python3
"""mu_inf, the D = infinity limit of the least-damped real ghost-sector rate, with a planted-truth test of the
extrapolation first (notes/wormhole_next.md section 6.4). Criteria written before the run on the step-2 data:

    python3 mu_inf.py [--records=ghost_mu.jsonl,ghost_mu2.jsonl] [--dmin=16]

Data: the least-damped real zero mu(D) of every mode-1, unplanted record with D >= dmin (one per D, the highest M).

Estimator (fixed here): least-squares fits mu(n) = L + sum_{j=1..k} c_j n^-j, n = D - 2, for k = 2, 3, 4 and for the
data from each of the three lowest D upward (9 fits; a fit needs at least k + 3 points). The estimate is the midpoint of
the fitted L's and its bar is their half-range: a spread over fit choices, not an error bound.

Planted truth (same D grid, same estimator):
  T1  analytic models: the real data refitted with L held at L0 = -1.00, -0.98 and -1.02 (terms 1/n .. 1/n^5), then
      evaluated on the grid. Pass: each estimate's interval contains L0 and its midpoint is within 0.01 of L0.
  T2  a non-analytic model: L0 = -1 + (a + b ln n)/n^2 + c/n^3 + d/n^4 fitted to the real data. Pass: the interval
      contains -1.
  T3  control: the T1 model with L0 = -1.00, shifted by -0.05 at every point; the recovered interval must not contain
      -1.00 (the estimator sees a 0.05 change).
The real-data estimate is reported only if T1-T3 pass; if not, the script says so and exits 1. Whether the real interval
contains -1 is reported, not checked.
"""
import json
import os
import sys

import numpy as np


def load(paths, dmin):
    pts = {}
    for p in paths:
        if not os.path.exists(p):
            continue
        for line in open(p):
            r = json.loads(line)
            if r["mode"] != 1 or r.get("plant") is not None or r["D"] < dmin:
                continue
            z = [x["mu"] for x in r["zeros"] if x["ok"]]
            if not z:
                continue
            D = round(r["D"], 4)
            if D not in pts or r["M"] > pts[D][0]:
                pts[D] = (r["M"], max(z))
    Ds = sorted(pts)
    return np.array(Ds), np.array([pts[D][1] for D in Ds])


def fit(n, y, k, L=None):
    """least squares in 1/n; with L given, it is held"""
    cols = [n ** -j for j in range(1, k + 1)]
    if L is None:
        A = np.stack([np.ones_like(n)] + cols, 1)
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        return c[0], c, float(np.abs(A @ c - y).max())
    A = np.stack(cols, 1)
    c, *_ = np.linalg.lstsq(A, y - L, rcond=None)
    return L, c, float(np.abs(A @ c - (y - L)).max())


def estimate(D, y):
    n = D - 2
    Ls = []
    for k in (2, 3, 4):
        for start in range(3):
            if len(n) - start < k + 3:
                continue
            L, _, res = fit(n[start:], y[start:], k)
            Ls.append((k, D[start], L, res))
    vals = [L for _, _, L, _ in Ls]
    return 0.5 * (max(vals) + min(vals)), 0.5 * (max(vals) - min(vals)), Ls


def main():
    opts = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a)
    here = os.path.dirname(os.path.abspath(__file__))
    recs = [os.path.join(here, f) for f in opts.get("records", "ghost_mu.jsonl,ghost_mu2.jsonl").split(",")]
    D, y = load(recs, float(opts.get("dmin", 16)))
    n = D - 2
    print(f"data: {len(D)} points, D = {[float(d) for d in D]}")
    bad = []

    def show(name, mid, half, Ls):
        print(f"  {name}: midpoint {mid:+.5f}, half-range {half:.5f} over {len(Ls)} fits "
              f"(L from {min(L for *_, L, _ in Ls):+.5f} to {max(L for *_, L, _ in Ls):+.5f})")

    for L0 in (-1.00, -0.98, -1.02):
        _, c, _ = fit(n, y, 5, L=L0)
        syn = L0 + sum(c[j - 1] * n ** -j for j in range(1, 6))
        mid, half, Ls = estimate(D, syn)
        show(f"T1 planted L0 = {L0:+.2f}", mid, half, Ls)
        ok = (mid - half - 1e-12 <= L0 <= mid + half + 1e-12) and abs(mid - L0) < 0.01
        print(f"{'PASS' if ok else 'FAIL'}  T1 L0 = {L0:+.2f}")
        bad += [] if ok else [f"T1 {L0}"]
    A = np.stack([n ** -2, np.log(n) * n ** -2, n ** -3, n ** -4], 1)
    c, *_ = np.linalg.lstsq(A, y + 1, rcond=None)
    syn = -1 + A @ c
    mid, half, Ls = estimate(D, syn)
    show("T2 planted -1 with (a + b ln n)/n^2", mid, half, Ls)
    ok = mid - half <= -1 <= mid + half
    print(f"{'PASS' if ok else 'FAIL'}  T2")
    bad += [] if ok else ["T2"]
    _, c, _ = fit(n, y, 5, L=-1.0)
    syn = -1.0 + sum(c[j - 1] * n ** -j for j in range(1, 6)) - 0.05
    mid, half, Ls = estimate(D, syn)
    show("T3 control, shifted by -0.05", mid, half, Ls)
    ok = not (mid - half <= -1.0 <= mid + half)
    print(f"{'PASS' if ok else 'FAIL'}  T3 control (the shift is seen)")
    bad += [] if ok else ["T3"]
    if bad:
        print(f"planted-truth test failed ({bad}): no estimate of mu_inf is quoted")
        sys.exit(1)
    mid, half, Ls = estimate(D, y)
    show("real data", mid, half, Ls)
    for k, d0, L, res in Ls:
        print(f"    terms to 1/n^{k}, from D = {d0:.2f}: L = {L:+.6f}, max residual {res:.1e}")
    print(f"mu_inf = {mid:+.4f} +- {half:.4f} (spread over fit choices); the interval "
          f"{'contains' if mid - half <= -1 <= mid + half else 'excludes'} -1")


if __name__ == "__main__":
    main()

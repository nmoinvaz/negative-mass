#!/usr/bin/env python3
"""Throughput of a stand-in for the planned dual-null cell update (#39 design; NOT the evolution code).

    python3 kernel_bench.py [N] [reps]

The cell update of the design (second-order predictor-corrector diamond, fields r, f, one ghost and one normal scalar,
n = D - 2 = 2) applied on an N x N grid of spacing h = 0.001 to the exact static Ellis data on the two initial rays. It
measures cells per second only. Nothing is checked or kept; the sweep length is short enough that the static solution has not
begun to leave (the unstable mode has grown by e^(t/0.846) with t <= 2.8). Three implementations of the same arithmetic:
numpy on anti-diagonals, numba on rows, and (kernel_bench.rs) Rust on rows.
"""
import sys
import time

import numpy as np
import numba as nb

NN = 2.0          # n = D - 2
E1, E2 = -1.0, 1.0


def ellis_rays(N, h):
    """Exact static data on x- = 0 (l = x+/sqrt2) and x+ = 0 (l = -x-/sqrt2): r = sqrt(1 + l^2), f = 0, ghost = sqrt2 atan l."""
    x = np.arange(N + 1) * h
    lp = x / np.sqrt(2)
    rp = np.sqrt(1 + lp**2)
    ap = np.sqrt(2) * np.arctan(lp)
    lm = -x / np.sqrt(2)
    rm = np.sqrt(1 + lm**2)
    am = np.sqrt(2) * np.arctan(lm)
    z = np.zeros(N + 1)
    return (rp, z, ap, z), (rm, z, am, z)


def F_np(rc, fc, ac, bc, rp, fp, ap, bp, rm, fm, am, bm):
    emf = np.exp(-fc)
    rpm = -(NN - 1) * (rp * rm + 0.5 * emf) / rc
    spm = E1 * ap * am + E2 * bp * bm
    fpm = (NN / rc) * rpm + spm
    k = -NN / (2 * rc)
    apm = k * (rp * am + rm * ap)
    bpm = k * (rp * bm + rm * bp)
    return rpm, fpm, apm, bpm


def cell_np(S, E, W, h):
    """S, E, W: tuples of 4 arrays. Returns the 4 arrays at N."""
    h2 = h * h
    c = [0.5 * (e + w) for e, w in zip(E, W)]
    up = [(e - s) / h for e, s in zip(E, S)]
    um = [(w - s) / h for w, s in zip(W, S)]
    rpm, fpm, apm, bpm = F_np(c[0], c[1], c[2], c[3], up[0], up[1], up[2], up[3], um[0], um[1], um[2], um[3])
    Nh = [e + w - s + h2 * q for e, w, s, q in zip(E, W, S, (rpm, fpm, apm, bpm))]
    c = [0.25 * (s + e + w + nn) for s, e, w, nn in zip(S, E, W, Nh)]
    up = [(e + nn - s - w) / (2 * h) for s, e, w, nn in zip(S, E, W, Nh)]
    um = [(w + nn - s - e) / (2 * h) for s, e, w, nn in zip(S, E, W, Nh)]
    rpm, fpm, apm, bpm = F_np(c[0], c[1], c[2], c[3], up[0], up[1], up[2], up[3], um[0], um[1], um[2], um[3])
    return [e + w - s + h2 * q for e, w, s, q in zip(E, W, S, (rpm, fpm, apm, bpm))]


def sweep_numpy(N, h):
    (Rp, Rm) = ellis_rays(N, h)
    Dk = [np.zeros(N + 1) for _ in range(4)]
    Dk1 = [np.zeros(N + 1) for _ in range(4)]
    for q in range(4):
        Dk[q][0] = Rp[q][0]
        Dk1[q][1] = Rp[q][1]
        Dk1[q][0] = Rm[q][1]
    for k in range(0, 2 * N - 1):
        Dk2 = [np.zeros(N + 1) for _ in range(4)]
        lo, hi = max(1, k + 2 - N), min(N, k + 1)
        if hi >= lo:
            S = [a[lo - 1:hi] for a in Dk]
            E = [a[lo:hi + 1] for a in Dk1]
            W = [a[lo - 1:hi] for a in Dk1]
            out = cell_np(S, E, W, h)
            for q in range(4):
                Dk2[q][lo:hi + 1] = out[q]
        if k + 2 <= N:
            for q in range(4):
                Dk2[q][k + 2] = Rp[q][k + 2]
                Dk2[q][0] = Rm[q][k + 2]
        Dk, Dk1 = Dk1, Dk2
    return Dk1[0][N]


@nb.njit(cache=True, fastmath=False)
def sweep_numba(N, h, rpR, apR, rmR, amR):
    old = np.zeros((4, N + 1)); new = np.zeros((4, N + 1))
    old[0, :] = rpR; old[2, :] = apR           # row j = 0: ray x- = 0
    h2 = h * h
    for j in range(N):
        new[0, 0] = rmR[j + 1]; new[2, 0] = amR[j + 1]; new[1, 0] = 0.0; new[3, 0] = 0.0
        for i in range(N):
            S0, S1, S2, S3 = old[0, i], old[1, i], old[2, i], old[3, i]
            E0, E1_, E2_, E3 = old[0, i + 1], old[1, i + 1], old[2, i + 1], old[3, i + 1]
            W0, W1, W2, W3 = new[0, i], new[1, i], new[2, i], new[3, i]
            # predictor
            rc = 0.5 * (E0 + W0); fc = 0.5 * (E1_ + W1)
            rp = (E0 - S0) / h; rm = (W0 - S0) / h
            ap = (E2_ - S2) / h; am = (W2 - S2) / h
            bp = (E3 - S3) / h; bm = (W3 - S3) / h
            emf = np.exp(-fc)
            rpm = -(NN - 1) * (rp * rm + 0.5 * emf) / rc
            spm = E1 * ap * am + E2 * bp * bm
            fpm = (NN / rc) * rpm + spm
            k = -NN / (2 * rc)
            apm = k * (rp * am + rm * ap); bpm = k * (rp * bm + rm * bp)
            N0 = E0 + W0 - S0 + h2 * rpm; N1 = E1_ + W1 - S1 + h2 * fpm
            N2 = E2_ + W2 - S2 + h2 * apm; N3 = E3 + W3 - S3 + h2 * bpm
            # corrector
            rc = 0.25 * (S0 + E0 + W0 + N0); fc = 0.25 * (S1 + E1_ + W1 + N1)
            rp = (E0 + N0 - S0 - W0) / (2 * h); rm = (W0 + N0 - S0 - E0) / (2 * h)
            ap = (E2_ + N2 - S2 - W2) / (2 * h); am = (W2 + N2 - S2 - E2_) / (2 * h)
            bp = (E3 + N3 - S3 - W3) / (2 * h); bm = (W3 + N3 - S3 - E3) / (2 * h)
            emf = np.exp(-fc)
            rpm = -(NN - 1) * (rp * rm + 0.5 * emf) / rc
            spm = E1 * ap * am + E2 * bp * bm
            fpm = (NN / rc) * rpm + spm
            k = -NN / (2 * rc)
            apm = k * (rp * am + rm * ap); bpm = k * (rp * bm + rm * bp)
            new[0, i + 1] = E0 + W0 - S0 + h2 * rpm; new[1, i + 1] = E1_ + W1 - S1 + h2 * fpm
            new[2, i + 1] = E2_ + W2 - S2 + h2 * apm; new[3, i + 1] = E3 + W3 - S3 + h2 * bpm
        old, new = new, old
    return old[0, N]


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    h = 0.001
    cells = N * N
    la = open("/proc/loadavg").read().split()[:3]
    print(f"N = {N}, cells = {cells}, h = {h}, load average at start {' '.join(la)}, one thread")
    (Rp, Rm) = ellis_rays(N, h)
    sweep_numba(8, h, Rp[0][:9].copy(), Rp[2][:9].copy(), Rm[0][:9].copy(), Rm[2][:9].copy())    # compile
    for name in ("numba", "numpy"):
        ts = []
        for _ in range(reps):
            t0 = time.perf_counter()
            if name == "numba":
                v = sweep_numba(N, h, Rp[0], Rp[2], Rm[0], Rm[2])
            else:
                v = sweep_numpy(N, h)
            ts.append(time.perf_counter() - t0)
        print(f"{name:6s} best {min(ts):8.3f} s  median {sorted(ts)[len(ts)//2]:8.3f} s  -> {1e9 * min(ts) / cells:7.1f} ns/cell "
              f"({cells / min(ts) / 1e6:6.2f} Mcell/s)  value {v:.6f}")

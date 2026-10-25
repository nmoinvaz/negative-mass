"""Checks of the compiled path (fastcoef.py; Wormhole(fast=True), the default) against the earlier one (fast=False:
lambdified sympy coefficients, the Frobenius step in mpmath at 30 digits).

    python3 fast_check.py [D ...]          (default 7 13; gamma1 = 1, eps2 = 0, the roots of largeD_roots_a.csv)

1. Coefficients A2, A1, A0, r, r_z, r_zz at z = 1e-6 ... 0.7 from both ends, and A2, A1, A0 at complex points near
   the throat: max relative difference <= 1e-13.
2. Frobenius: the series terms c_k z0^(k + r0) of both throat solutions at the det's z0, fast (double) against
   mpmath, relative to the largest term: <= 1e-13.
3. det (amp form) at What = root + 0.01 and root - 0.01: |fast - slow| <= 1e-9 (|det| <= 1).
   Control: slow at What (1 + 1e-7) against fast at What must differ by more than 1e-9.
4. The root by brentq on the fast det (largeD_roots.root) against largeD_roots_a.csv (slow): relative <= 1e-10.
5. modes_complex: det_c at What = root (1 + 0.1 i), fast against slow: relative <= 1e-8; control at What (1 + 1e-7).
Exit 1 if a check fails or a control passes.
"""
import csv
import sys
import time

import numpy as np

import largeD_roots as lr
import modes_complex as mc
import modes_general as mg


def rel(a, b):
    a, b = np.asarray(a, complex), np.asarray(b, complex)
    return float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300)))


def main():
    Ds = [int(a) for a in sys.argv[1:]] or [7, 13]
    csv_roots = {int(r[0]): float(r[1]) for r in csv.reader(open("largeD_roots_a.csv")) if r}
    bad = []
    for D in Ds:
        t = time.time()
        wf, ws = mg.Wormhole(D, 1.0, 0.0), mg.Wormhole(D, 1.0, 0.0, fast=False)
        root = csv_roots[D]
        W = np.sqrt(wf.W2_of(root))
        # 1. coefficients
        e = 0.0
        for side in (-1, 1):
            for zz in (1e-6, 1e-3, 0.05, 0.3, 0.7):
                e = max(e, rel(wf.coeffs_z(zz, W, side), ws.coeffs_z(zz, W, side)))
                e = max(e, rel([wf.r_z(zz, side, k) for k in range(3)], [ws.r_z(zz, side, k) for k in range(3)]))
        p = (wf.n, wf.c, wf.g1, complex(wf.e2))
        for dz in (0.01, 0.01j, -0.02 + 0.005j):
            e = max(e, rel(wf.F.coef_x(wf.yth + dz, *p, complex(W)), [f(wf.yth + dz, *p, W + 0j) for f in ws.f_np]))
        print(f"D {D}: 1. coefficients max relative difference {e:.1e}", flush=True)
        if e > 1e-13:
            bad.append(f"D {D}: check 1")
        # 2. Frobenius
        sf, ss = wf.frobenius(W), ws.frobenius(W)
        # the coefficients enter only through the series at z0 (Wormhole.local): compare c_k z0^(k + r0), relative to
        # the largest term (a high coefficient on its own carries the double-precision error of the Cauchy integral,
        # 1e-16 rho^-k, which z0^k < (0.4 rho)^k suppresses)
        z0 = 0.02 * min(1.0, np.pi / 2 - abs(wf.yth))
        e2 = 0.0
        for r0 in (0, 3):
            tf = np.array(sf[r0]) * z0 ** (np.arange(len(sf[r0])) + r0)
            ts = np.array(ss[r0]) * z0 ** (np.arange(len(ss[r0])) + r0)
            e2 = max(e2, float(np.max(np.abs(tf - ts)) / np.max(np.abs(ts))))
        print(f"   2. Frobenius series terms at z0 max difference {e2:.1e} of the largest (resonance rhs fast "
              f"{sf.get('resonance_rhs'):.1e}, mpmath {ss.get('resonance_rhs'):.1e})", flush=True)
        if e2 > 1e-13:
            bad.append(f"D {D}: check 2")
        # 3. det, with the control
        for h in (root + 0.01, root - 0.01):
            t1 = time.time()
            df = wf.det(wf.W2_of(h))[0]
            t2 = time.time()
            ds = ws.det(ws.W2_of(h))[0]
            t3 = time.time()
            print(f"   3. det at What {h:.4f}: fast {df:+.15f} ({t2 - t1:.1f} s), slow {ds:+.15f} ({t3 - t2:.1f} s), "
                  f"difference {df - ds:+.1e}", flush=True)
            if abs(df - ds) > 1e-9:
                bad.append(f"D {D}: check 3")
        dctl = ws.det(ws.W2_of((root + 0.01) * (1 + 1e-7)))[0]
        df = wf.det(wf.W2_of(root + 0.01))[0]
        print(f"      control: slow at What (1 + 1e-7) differs by {dctl - df:+.1e}", flush=True)
        if abs(dctl - df) <= 1e-9:
            bad.append(f"D {D}: control 3 passed")
        # 4. root
        t1 = time.time()
        r = lr.root(wf, root)
        print(f"   4. root fast {r:.13f} ({time.time() - t1:.1f} s), largeD_roots_a.csv {root:.13f}: relative "
              f"{(r - root) / root:+.1e}", flush=True)
        if not abs(r - root) <= 1e-10 * root:
            bad.append(f"D {D}: check 4")
        # 5. complex solver
        cf, cs = mc.ComplexWormhole(D, 1.0, 0.0), mc.ComplexWormhole(D, 1.0, 0.0, fast=False)
        g = root * (1 + 0.1j)
        Wg = g * np.exp(cf.nu_th) / cf.a
        margins = tuple(cf.margin(Wg, s) for s in (-1, 1))
        t1 = time.time()
        vf = cf.det_c(Wg ** 2, margins)
        t2 = time.time()
        vs = cs.det_c(Wg ** 2, margins)
        t3 = time.time()
        vc = cs.det_c((Wg * (1 + 1e-7)) ** 2, margins)
        e5, e5c = rel(vf, vs), rel(vc, vf)
        print(f"   5. det_c at What {g:.4f}: fast {vf:.12e} ({t2 - t1:.1f} s), slow {vs:.12e} ({t3 - t2:.1f} s), "
              f"relative {e5:.1e}; control (1 + 1e-7) {e5c:.1e}  [{time.time() - t:.0f} s]", flush=True)
        if e5 > 1e-8:
            bad.append(f"D {D}: check 5")
        if e5c <= 1e-8:
            bad.append(f"D {D}: control 5 passed")
    print("PASS" if not bad else "FAIL: " + "; ".join(bad))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

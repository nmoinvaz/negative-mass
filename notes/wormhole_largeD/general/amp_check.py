"""Checks of the amplitude form of Wormhole.det (form="amp") against the earlier e^psi cos(th) form, on the
asymmetric wormhole (gamma1 = 1, eps2 = 0).

  1. Scale: the amplitude form at three What around the root is O(1); the cos form is printed beside it, and the
     two solutions' end angles modulo pi agree to 1e-6 (the amplitude form's assumption).
  2. Agreement: where the cos form is well resolved (D = 8, 10) both forms give the same root to 1e-8 (relative).
  3. Control: psi_shift = 1e-3 on one solution's amplitude must move the amplitude form's root by more than 1e-6,
     so check 2 can fail.
  4. Convergence (--conv): the root with the end pushed out (rmax_factor 90, X 45) and with rtol 1e-13 agrees with
     the default to 1e-8.

    python3 amp_check.py D lo hi [--conv]        (a bracket [lo, hi] around the root; exit 1 if a check fails)
"""
import sys
import time

import numpy as np
from scipy.optimize import brentq

import modes_general as mg


def root(w, lo, hi, **kw):
    f = lambda h: w.det(w.W2_of(h), **kw)[0]  # noqa: E731
    return brentq(f, lo, hi, xtol=1e-12)


def main():
    D, lo, hi = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    w = mg.Wormhole(D, 1.0, 0.0)
    bad = []
    t = time.time()
    mid = 0.5 * (lo + hi)
    for h in (lo, mid, hi):
        (a, sols), c = w.det(w.W2_of(h)), w.det(w.W2_of(h), form="cos")[0]
        print(f"D {D} What {h:.6f}: amp {a:+.6e}   cos {c:+.3e}   th spread {sols['th_spread']:.1e}", flush=True)
        if sols["th_spread"] > 1e-6:
            bad.append(f"th spread {sols['th_spread']:.1e} at What {h}: the solutions are not on one direction")
    r_amp = root(w, lo, hi)
    r_cos = root(w, lo, hi, form="cos")
    print(f"D {D}: root amp {r_amp:.12f}, cos {r_cos:.12f}, rel diff {abs(r_amp - r_cos) / r_amp:.1e}", flush=True)
    if abs(r_amp - r_cos) > 1e-8 * r_amp:
        bad.append("amp and cos roots differ")
    r_mut = root(w, lo, hi, psi_shift=1e-3)
    print(f"D {D}: control psi_shift 1e-3: root {r_mut:.12f}, moved {r_mut - r_amp:+.2e}", flush=True)
    if abs(r_mut - r_amp) <= 1e-6 * r_amp:
        bad.append("the control did not move the root")
    if "--conv" in sys.argv:
        r_far = root(w, lo, hi, rmax_factor=90.0, X=45.0)
        r_tol = root(w, lo, hi, rtol=1e-13)
        print(f"D {D}: rmax 90: {r_far:.12f} ({(r_far - r_amp) / r_amp:+.1e}); rtol 1e-13: {r_tol:.12f} "
              f"({(r_tol - r_amp) / r_amp:+.1e})", flush=True)
        if max(abs(r_far - r_amp), abs(r_tol - r_amp)) > 1e-8 * r_amp:
            bad.append("not converged in the end radius or the tolerance")
    print(f"{time.time() - t:.0f} s; " + ("PASS" if not bad else "FAIL: " + "; ".join(bad)), flush=True)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

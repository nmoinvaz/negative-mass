"""Checks of Wormhole.det on the asymmetric wormhole (gamma1 = 1, eps2 = 0): the amplitude form (form="amp") and the
shooting in z = pi/2 - side y with the phase of (Y, Y'/g) (coord="z"), against the earlier e^psi cos(th) form and
the earlier shooting in y (coord="y").

  1. Scale: the amplitude form at three What around the root is O(1); the cos form is printed beside it, and the
     two solutions' end angles modulo pi agree to 1e-6 (the amplitude form's assumption).
  2. Smoothness: at five What 2e-6 apart (relative) around the root, the determinant is a quadratic in What to 1e-9.
     Control: the same five points with coord="y", whose right end psi carries noise of 1e-7 at D = 8 and 1e-4 at
     D = 10, must fail it. From D = 12 the y shooting takes hours, and the control adds psi_shift = +-1e-8, alternating
     between the five points (a looser rtol is no control: at 1e-6 it moves the determinant but leaves it smooth).
  3. Agreement: the amp and cos forms give the same root to 1e-8 (relative); with --y the coord="y" root must agree
     to 1e-8 too (only at D = 8: its noise moves the root by 5e-6 at D = 10).
  4. Control: psi_shift = 1e-3 on one solution's amplitude must move the amplitude form's root by more than 1e-6,
     so check 3 can fail.
  5. Convergence (--conv): the root with the end pushed out (rmax_factor 90, X 45) and with rtol 1e-13 agrees with
     the default to 1e-8.

    python3 amp_check.py D lo hi [--conv] [--y]   (a bracket [lo, hi] around the root; exit 1 if a check fails)
"""
import sys
import time

import numpy as np
from scipy.optimize import brentq

import modes_general as mg


def root(w, lo, hi, **kw):
    f = lambda h: w.det(w.W2_of(h), **kw)[0]  # noqa: E731
    return brentq(f, lo, hi, xtol=1e-12)


def roughness(w, r, noise=0.0, **kw):
    """largest residual of a quadratic fit to the determinant at five What 2e-6 apart (relative) around r; noise is a
    psi_shift of alternating sign at the five points (for the control)"""
    hs = r * (1 + 2e-6 * np.arange(-2, 3))
    vals = [w.det(w.W2_of(h), psi_shift=noise * (-1) ** i, **kw)[0] for i, h in enumerate(hs)]
    x = (hs - r) / (hs[-1] - hs[0])
    return float(np.max(np.abs(np.polyval(np.polyfit(x, vals, 2), x) - vals))), vals


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
    rough, vals = roughness(w, r_amp)
    ctl = dict(coord="y") if D <= 10 else dict(noise=1e-8)
    rough_y, vals_y = roughness(w, r_amp, **ctl)
    print(f"D {D}: quadratic residual z {rough:.1e} (det {vals[0]:+.6e} .. {vals[-1]:+.6e}); control {ctl} "
          f"{rough_y:.1e} (det {vals_y[0]:+.6e} .. {vals_y[-1]:+.6e})", flush=True)
    if rough > 1e-9:
        bad.append("the z determinant is not smooth")
    if rough_y <= 1e-9:
        bad.append(f"the control {ctl} passed the smoothness check")
    r_mut = root(w, lo, hi, psi_shift=1e-3)
    print(f"D {D}: control psi_shift 1e-3: root {r_mut:.12f}, moved {r_mut - r_amp:+.2e}", flush=True)
    if abs(r_mut - r_amp) <= 1e-6 * r_amp:
        bad.append("the control did not move the root")
    if "--y" in sys.argv:
        r_y = root(w, lo, hi, coord="y")
        print(f"D {D}: root y {r_y:.12f} ({(r_y - r_amp) / r_amp:+.1e})", flush=True)
        if abs(r_y - r_amp) > 1e-8 * r_amp:
            bad.append("the z and y roots differ")
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

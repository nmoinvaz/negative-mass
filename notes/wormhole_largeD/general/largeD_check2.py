"""Exact growth rate of the asymmetric (or charged) wormhole at one D, near a predicted value, for the large-D check:
What = W_th a against (D - 3) sqrt(omega) + offset. A narrow bracket around the prediction keeps the number of
determinant evaluations (each a Frobenius expansion at 30 digits and four shootings) small.

    python3 largeD_check2.py D gamma1 eps2 What_guess [half_width=0.3] [rmax_factor=60]
"""
import sys
import time

import numpy as np
from scipy.optimize import brentq

import modes_general as mg

D, g1, e2, guess = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
hw = float(sys.argv[5]) if len(sys.argv) > 5 else 0.3
rmf = float(sys.argv[6]) if len(sys.argv) > 6 else 60.0
t = time.time()
w = mg.Wormhole(D, g1, e2)
f = lambda h: w.det(w.W2_of(h), rmax_factor=rmf, X=rmf / 2)[0]  # noqa: E731
grid = np.linspace(guess - hw, guess + hw, 4)
vals = [f(h) for h in grid]
root = float("nan")
for (a, fa), (b, fb) in zip(zip(grid, vals), zip(grid[1:], vals[1:])):
    if np.isfinite(fa) and np.isfinite(fb) and np.sign(fa) != np.sign(fb):
        root = brentq(f, a, b, xtol=1e-9)
        break
print(f"D = {D}, gamma1 = {g1}, eps2 = {e2}, rmax_factor {rmf}: What = {root:.10f}, What/(D-3) = {root / (D - 3):.8f}"
      f" ({time.time() - t:.0f} s; bracket values {['%.2e' % v for v in vals]})", flush=True)

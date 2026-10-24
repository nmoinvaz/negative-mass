"""Exact growth rates of the asymmetric wormhole at several D, against the leading large-D result (D-3) sqrt(omega)."""
import sys
import numpy as np
from scipy.optimize import brentq
import modes_general as mg
g1 = float(sys.argv[1]); slope = float(sys.argv[2])
for D in [int(a) for a in sys.argv[3:]]:
    w = mg.Wormhole(D, g1, 0.0)
    guess = slope * (D - 3)
    grid = np.linspace(0.6 * guess, 1.4 * guess + 1.5, 9)
    out = []
    for rmf in (60.0, 90.0):
        f = lambda h: w.det(w.W2_of(h), rmax_factor=rmf, X=30.0 if rmf == 60 else 45.0)[0]
        vals = [f(h) for h in grid]
        root = None
        for (a, fa), (b, fb) in zip(zip(grid, vals), zip(grid[1:], vals[1:])):
            if np.sign(fa) != np.sign(fb):
                root = brentq(f, a, b, xtol=1e-10)
                break
        out.append(root)
    r = out[0]
    print(f"D = {D:3d}: What = {r}  (larger end radius: {out[1]})  What/(D-3) = {r / (D - 3) if r else None}"
          f"  What - (D-3) sqrt(om) = {r - guess if r else None}", flush=True)

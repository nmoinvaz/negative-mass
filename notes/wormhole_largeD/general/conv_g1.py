import sys
import numpy as np
from scipy.optimize import brentq
import modes_general as mg
w = mg.Wormhole(4, 1.0, 0.0)
for rmf, z0 in ((25, None), (40, None), (25, 0.005), (60, 0.005)):
    f = lambda h: w.det(w.W2_of(h), rmax_factor=rmf, z0=z0)[0]
    grid = np.linspace(1.40, 1.60, 6)
    vals = [f(h) for h in grid]
    print(rmf, z0, ["%+.3e" % v for v in vals], flush=True)
    for (a, fa), (b, fb) in zip(zip(grid, vals), zip(grid[1:], vals[1:])):
        if np.sign(fa) != np.sign(fb):
            r = brentq(f, a, b, xtol=1e-10)
            print("   root What =", r, " T =", 1 / r, flush=True)

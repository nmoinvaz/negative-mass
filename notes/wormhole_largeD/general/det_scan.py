"""The determinant of largeD_check2.py at given What, one line each, to tell a root from a pole of its sign change.

    python3 det_scan.py D gamma1 eps2 What1,What2,... [rmax_factor=60]
"""
import sys
import time

import modes_general as mg

D, g1, e2 = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
hs = [float(v) for v in sys.argv[4].split(",")]
rmf = float(sys.argv[5]) if len(sys.argv) > 5 else 60.0
w = mg.Wormhole(D, g1, e2)
for h in hs:
    t = time.time()
    print(f"D = {D}, gamma1 = {g1}, eps2 = {e2}, rmax_factor {rmf}: What = {h:.4f} det = "
          f"{w.det(w.W2_of(h), rmax_factor=rmf, X=rmf / 2)[0]:.6e} ({time.time() - t:.0f} s)", flush=True)

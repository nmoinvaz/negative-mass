"""
Matter-sector diagnostic.  At r_m, a growing mode of the scalar sector alone (e.g. a Q-ball-type
instability, metric weakly involved) would make the regular interior 3-plane's scalar part meet the
decaying scalar 2-plane in more than the generic 1 dimension: the second-smallest singular value of
A = (I - S S^+) Y_scalar (4x3, rank <= 2 generically 2) would go to 0.  We tabulate
sv2/sv1 of A over gamma (sv3 ~ 0 is the generic null vector used by kappa_eff.py).
Also: flat-space l=1 is excluded (l>=2 only).  Usage: python matter_sv.py eps w l
"""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem
eps, w, l = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
P = PolarProblem(B.solve_star(w, eps), l)
gs = np.logspace(-3, np.log10(0.6), 25)      # scalar continuum threshold: bound scalar modes need |s| small vs 1-w
out = []
for g in gs:
    _, sv = P.interior_zerilli(g)
    out.append((g, sv[1] / sv[0], sv[2] / sv[0]))
out = np.array(out)
i = np.argmin(out[:, 1])
print(f"eps={eps:+d} w={w} l={l}: min sv2/sv1 = {out[i,1]:.3e} at gamma={out[i,0]:.3g};  max sv3/sv1 = {out[:,2].max():.1e}")
for row in out[::4]:
    print(f"   gamma={row[0]:.3e}  sv2/sv1={row[1]:.3e}")

"""Sign of D on a deep negative sigma^2 range (growth rates up to |sigma| = 2) for the campaign stars."""
import numpy as np, warnings
warnings.filterwarnings('ignore')
from background import solve_star
from radial_modes import RadialProblem
for eps, w in [(+1, 0.75), (-1, 0.7603)]:
    st = solve_star(0.75, eps)
    if w != 0.75: st = solve_star(w, eps, guess=st)
    P = RadialProblem(st)
    g = -np.logspace(-6, np.log10(4.0), 40)
    v = np.array([P.D(s) for s in g])
    print(f"eps={eps:+d} w={w}: sigma^2 in [-4,-1e-6]: sign changes of D = {int(np.sum(np.diff(np.sign(v)) != 0))}, "
          f"min|D|={np.min(np.abs(v)):.2e}", flush=True)

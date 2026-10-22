"""Canonical (eps=+1) family: M, N and lowest radial sigma^2 vs w."""
import numpy as np, sys, json, warnings
warnings.filterwarnings('ignore')
from background import solve_star
from radial_modes import RadialProblem, lowest_modes
eps = int(sys.argv[1]) if len(sys.argv) > 1 else 1
ws = [float(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else \
     [0.75, 0.80, 0.85, 0.88, 0.90, 0.91, 0.92, 0.925, 0.93, 0.935, 0.94, 0.95, 0.96]
st = solve_star(0.75, eps)
cur = 0.75
rows = []
for w in sorted(ws):
    # continuation from the previous star
    while abs(w - cur) > 1e-12:
        step = np.clip(w - cur, -0.02, 0.02)
        cur = round(cur + step, 6)
        st = solve_star(cur, eps, guess=st)
    P = RadialProblem(st)
    roots, g, vals = lowest_modes(P)
    rows.append(dict(w=w, M=st.M, N=st.N, phic=st.phic, roots=[(float(x), bool(o)) for x, o in roots]))
    print(f"eps={eps:+d} w={w:.4f} M={st.M:+.7f} N={st.N:.7f} sigma2={['%.6e%s' % (x, '' if o else '?') for x, o in roots]}"
          f"  thr={(1-w)**2:.4e} rm={P.rm:.1f} R={P.R:.0f}", flush=True)
json.dump(rows, open(f"scan_eps{eps:+d}.json" if len(sys.argv) < 3 or len(ws) > 10 else f"scan_eps{eps:+d}_fine.json", 'w'), indent=1)

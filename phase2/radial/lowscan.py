"""Radial spectrum on the low-w parts of both families: continue from w=0.5 in particle number N
(w free, passes the phantom w-fold) and compute the lowest sigma^2 at the listed N."""
import numpy as np, warnings, sys, pickle
warnings.filterwarnings('ignore')
from background import solve_star, solve_star_N, StarN
from radial_modes import RadialProblem, lowest_modes
NPOS = int(sys.argv[3]) if len(sys.argv) > 3 else 16
eps = int(sys.argv[1]); Ns = [float(x) for x in sys.argv[2].split(',')]
st = solve_star(0.66, eps, tol=1e-8); w = 0.66; dw = 0.01
while w > 0.5 + 1e-9:
    w2 = round(max(w - dw, 0.5), 6)
    try: st = solve_star(w2, eps, guess=st, tol=1e-8); w = w2; dw = min(dw * 1.5, 0.01)
    except RuntimeError: dw /= 2
lnN = np.log(st.N)
for Nt in Ns:
    dl = 0.05
    while lnN < np.log(Nt) - 1e-12:
        step = min(dl, np.log(Nt) - lnN)
        try: st = solve_star_N(np.exp(lnN + step), eps, st); lnN += step; dl = min(dl * 1.3, 0.05)
        except RuntimeError:
            dl /= 2
            if dl < 1e-4: raise
    P = RadialProblem(st)
    roots, g, vals = lowest_modes(P, s2_min=-0.05, nneg=16, npos=NPOS)
    print(f"eps={eps:+d} N={st.N:.5f} w={st.w:.6f} M={st.M:+.6f} sigma2={['%.5e%s' % (x, '' if o else '?') for x, o in roots]} "
          f"thr={(1-st.w)**2:.4f} rm={P.rm:.1f} R={P.R:.0f}", flush=True)

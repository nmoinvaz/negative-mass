"""Fine scan around the high-w turning point: compare zero of sigma^2 with extrema of M and N,
and test the first law dM = eps * w * dN."""
import numpy as np, sys, warnings, json
warnings.filterwarnings('ignore')
from background import solve_star
from radial_modes import RadialProblem, lowest_modes
from scipy.optimize import brentq
eps = int(sys.argv[1]); w0, w1, nw = 0.9230, 0.9310, 17
ws = np.linspace(w0, w1, nw)
st = solve_star(0.75, eps); cur = 0.75
rows = []
for w in ws:
    while abs(w - cur) > 1e-12:
        cur = cur + np.clip(w - cur, -0.02, 0.02); st = solve_star(cur, eps, guess=st)
    P = RadialProblem(st)
    roots, _, _ = lowest_modes(P, s2_min=-0.01, nneg=8, npos=8, nmax=1)
    s2 = roots[0][0]
    rows.append((w, st.M, st.N, s2)); print("%.5f %+.9f %.9f %+.6e" % rows[-1], flush=True)
a = np.array(rows)
w, M, N, s2 = a.T
cM = np.polyfit(w - 0.927, M, 4); cN = np.polyfit(w - 0.927, N, 4); cS = np.polyfit(w - 0.927, s2, 3)
from numpy.polynomial import polynomial as Pn
def root(c):
    rr = np.roots(c); rr = rr[np.isreal(rr)].real; rr = rr[np.abs(rr) < 0.005]; return rr + 0.927
wM = root(np.polyder(cM)); wN = root(np.polyder(cN)); wS = root(cS)
# first law check via centered differences
dM = np.gradient(M, w); dN = np.gradient(N, w)
print(f"eps={eps:+d}: dM/dw=0 at w={wM}, dN/dw=0 at w={wN}, sigma^2=0 at w={wS}")
print("first law  (dM/dw)/(eps w dN/dw) - 1 at interior points away from extremum:",
      np.array2string(((dM / (eps * w * dN)) - 1)[[1, 2, -3, -2]], precision=2))
json.dump(dict(rows=rows, wM=list(wM), wN=list(wN), wS=list(wS)), open(f"crossing_eps{eps:+d}.json", "w"), indent=1)

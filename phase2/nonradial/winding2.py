"""
Argument-principle count of growing modes with complex s (oscillatory instabilities), done so that the
quantity integrated is HOLOMORPHIC in s up to a zero-free holomorphic factor:
  * scalar decaying plane started at a FIXED radius R_s (polar_solver Rs_fixed), and multiplied by the
    phase of exp[(k+ + k-)(R_s - r_m)], k± = sqrt(1 + (s ∓ i w)^2) (principal branch; its cut lies on Re s = 0,
    outside the contour), which removes the trivial decay of the 2 scalar columns;
  * Zerilli column = (1, q/f) with q the decaying Riccati value at r_m (holomorphic, r_start-independent);
  * interior plane from fixed r0 to fixed r_m.
D(conj s) = conj D(s) and D is real of one sign on the real axis (scan_growing.py), so
  #zeros in {x0<Re s<x1, |Im s|<y1} = 2 * Δarg(D along x1 -> x1+iy1 -> x0+iy1 -> x0) / 2π.
Per-evaluation cache (restart-safe).  Usage: python winding2.py eps w l x0 x1 y1
"""
import sys, os, numpy as np, time, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem
eps, w, l = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
x0, x1, y1 = (float(a) for a in sys.argv[4:7])
st = B.solve_star(w, eps)
P = PolarProblem(st, l, rtol=float(os.environ.get('RTOL', '1e-8')))
# fixed scalar start radius: 40 decay lengths at the slowest-decaying point of the contour
ks = lambda s: (np.sqrt(1 + (s - 1j*w)**2), np.sqrt(1 + (s + 1j*w)**2))
cands = [x0 + 1j*y for y in np.linspace(0, y1, 200)] + [x + 1j*y1 for x in np.linspace(x0, x1, 200)]
kmin = min(min(k.real for k in ks(c)) for c in cands)
P.Rs_fixed = P.rm + 40.0 / kmin
tag = f"eps{eps:+d}_w{w}_l{l}_{x0}_{x1}_{y1}"
cf = os.path.join(HERE, f"winding2_cache_{tag}.pkl")
cache = pickle.load(open(cf, 'rb')) if os.path.exists(cf) else {}
T0 = time.time()
def Dc(s):
    s = complex(round(s.real, 12), round(s.imag, 12))
    if s not in cache:
        kp, km = ks(s)
        # only the PHASE of the trivial factor is removed (its modulus is a positive number: irrelevant for arg)
        # inward integration GROWS the decaying columns by e^{+k(R_s-r_m)}; the regular interior columns grow
        # ~ e^{(k+ + k- + s) r_m}: remove both trivial (zero-free, holomorphic) phases
        cache[s] = P.D(s) * np.exp(-1j * ((kp + km) * (P.Rs_fixed - P.rm) + (kp + km + s) * P.rm).imag)
        pickle.dump(cache, open(cf, 'wb'))
    return cache[s]
print(f"{tag}: r_m={P.rm:.2f}, R_s fixed = {P.Rs_fixed:.1f} (min Re k on contour {kmin:.3f})", flush=True)
corners = [complex(x1, 0), complex(x1, y1), complex(x0, y1), complex(x0, 0)]
total, maxstep, nmin = 0.0, 0.0, np.inf
for a, b in zip(corners[:-1], corners[1:]):
    seg = list(np.linspace(0, 1, 11)); i = 0
    while i < len(seg) - 1:
        da, db = Dc(a + (b - a)*seg[i]), Dc(a + (b - a)*seg[i + 1])
        dph = np.angle(db / da)
        if abs(dph) > 0.4 and seg[i + 1] - seg[i] > 1e-5:
            seg.insert(i + 1, 0.5*(seg[i] + seg[i + 1])); continue
        total += dph; maxstep = max(maxstep, abs(dph)); nmin = min(nmin, abs(da)); i += 1
    print(f"   edge {a}->{b}: cumulative phase {total:+.4f}  ({len(cache)} evals, {time.time()-T0:.0f}s)", flush=True)
print(f"   D(x0)={Dc(complex(x0,0)).real:+.3e}  D(x1)={Dc(complex(x1,0)).real:+.3e} (same sign required)")
print(f"RESULT {tag}: zeros in Re s∈[{x0},{x1}], |Im s|≤{y1}: winding = {2*total/(2*np.pi):+.4f}  "
      f"(max phase step {maxstep:.2f} rad, min|D| {nmin:.2e})")

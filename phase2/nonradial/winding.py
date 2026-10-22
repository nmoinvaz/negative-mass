"""
Complex growth rates: count zeros of the matching determinant in a rectangle of the right half s-plane
by the argument principle.  D(s) = det[Y|Z] x (phases of the QR diagonals) equals the holomorphic
(unnormalised) determinant divided by a POSITIVE function, so arg D = arg of a holomorphic function of s
on Re s > 0 (scalar branch points s = +-i(w +- 1) and the Zerilli exterior are all on/left of Re s = 0).
Winding number of D around the rectangle Re s in [x0, x1], Im s in [-y1, y1] = number of growing modes
inside (oscillatory or not).  Adaptive: the contour is refined until every phase step is < 0.3 rad.
Usage: python winding.py eps w l x0 x1 y1
"""
import sys, os, numpy as np, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem

eps, w, l = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
x0, x1, y1 = (float(a) for a in sys.argv[4:7])
st = B.solve_star(w, eps)
P = PolarProblem(st, l, rtol=float(os.environ.get('RTOL', '1e-8')))
T0 = time.time()
# D(conj s) = conj D(s) and D is real with one sign on the real axis (scan_growing: no sign change), so
# the winding around the symmetric rectangle = 2 x (phase change along the UPPER half: x1 -> x1+iy1 ->
# x0+iy1 -> x0) / 2pi, provided D(x0), D(x1) have the same sign (checked below).
corners = [x1 + 0j, x1 + 1j * y1, x0 + 1j * y1, x0 + 0j]
pts, vals = [], []
cache = {}
import pickle as _pk
_cf = os.path.join(HERE, f"winding_cache_eps{eps:+d}_l{l}.pkl")
if os.path.exists(_cf):
    cache.update(_pk.load(open(_cf, 'rb')))
def Dc(s):
    s = complex(round(s.real, 12), round(s.imag, 12))
    if s not in cache:
        cache[s] = P.D(s)
        _pk.dump(cache, open(_cf, 'wb'))
    return cache[s]
# Strip the zero-free holomorphic factor e^{sL} that comes from normalising solutions at different radii
# (arg D grows ~L Im s). L is fitted once from two nearby points; e^{-sL} has no zeros, so the winding is
# unchanged, but far fewer contour points are needed.
Draw = Dc
_d0, _d1 = Draw(complex(x1, 0.0)), Draw(complex(x1, 0.01))
Lfit = float(np.angle(_d1 / _d0) / 0.01)
for _k in range(3):    # refine L on a longer baseline
    _yb = 0.2 * (_k + 1)
    _db = Draw(complex(x1, _yb))
    Lfit = Lfit + float(np.angle(_db * np.exp(-1j * Lfit * _yb) / _d0) / _yb)
print(f"L (stripped phase rate) = {Lfit:.4f}")
def Dc(s):
    return Draw(s) * np.exp(-complex(s) * Lfit)
total = 0.0
maxstep = 0.0
import json
ckpt = os.path.join(HERE, f"winding_ckpt_eps{eps:+d}_l{l}_{x0}_{x1}_{y1}.json")
done = json.load(open(ckpt)) if os.path.exists(ckpt) else {}
for ie, (a, b) in enumerate(zip(corners[:-1], corners[1:])):
    if str(ie) in done:
        total += done[str(ie)][0]; maxstep = max(maxstep, done[str(ie)][1]); print(f"edge {ie}: from checkpoint {done[str(ie)]}"); continue
    tot_e, max_e = 0.0, 0.0
    seg = list(np.linspace(0, 1, 13))
    i = 0
    while i < len(seg) - 1:
        sa, sb = a + (b - a) * seg[i], a + (b - a) * seg[i + 1]
        da, db = Dc(sa), Dc(sb)
        dphi = np.angle(db / da)
        if abs(dphi) > 0.3 and seg[i + 1] - seg[i] > 1e-4:
            seg.insert(i + 1, 0.5 * (seg[i] + seg[i + 1]))
            continue
        tot_e += dphi; max_e = max(max_e, abs(dphi))
        i += 1
    total += tot_e; maxstep = max(maxstep, max_e)
    done[str(ie)] = [tot_e, max_e, min(abs(v) for v in cache.values())]
    json.dump(done, open(ckpt, 'w')); print(f"edge {ie}: phase {tot_e:+.4f}  ({time.time()-T0:.0f}s)", flush=True)
n = 2 * total / (2 * np.pi)
print("D(x0), D(x1) real parts:", Dc(x0 + 0j).real, Dc(x1 + 0j).real)
print(f"eps={eps:+d} w={w} l={l}: rectangle Re s in [{x0},{x1}], Im s in [-{y1},{y1}] (via upper half x conj. symmetry): winding = {n:+.4f} "
      f"({len(cache)} evaluations, max phase step {maxstep:.3f} rad, min|D| on contour {min(abs(v) for v in cache.values()):.2e})  {time.time()-T0:.0f}s")

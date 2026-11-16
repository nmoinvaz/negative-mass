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
Per-evaluation cache (restart-safe).  Usage: python winding2.py eps w l x0 x1 y1 [options]

Options (the defaults give the behaviour of the first version):
  --y0=Y0         lower edge. Y0 > 0 counts the upper box x0<Re s<x1, Y0<Im s<y1 on its four edges, counter-
                  clockwise, with no doubling (its mirror image has the same count). Default 0, the symmetric box.
  --rsmult=40     R_s = r_m + rsmult / (min Re k on the contour).  Two values give the "two outer radii" check.
  --nseg=10       initial sub-segments per edge (the refinement below 0.4 rad per step follows).
  --cachedir=DIR  directory of the per-evaluation cache (default: beside this script).
  --plant=RE,IM   planted-zero check.  Counts D and D*(s-s0)[*(s-conj s0) if Y0 = 0] on the same cached
                  evaluations, and exits 1 unless the two counts differ by the number of planted zeros inside
                  the box (2 for a pair in the symmetric box, 1 in an upper box, 0 outside) within 0.05.
  --mutate=NAME   breaks the counter so that the --plant check must fail:
                    reverse      traverse the contour clockwise (the sign of the count flips),
                    wrongdouble  doubling applied to an upper box, or left out of a symmetric one,
                    misplant     the planted factor's zero is moved out of the box, the expectation is not,
                    onesided     (s-s0) without (s-conj s0) in the symmetric box (breaks D(conj s) = conj D(s)).
"""
import sys, os, numpy as np, time, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opts = dict(a[2:].split('=', 1) for a in sys.argv[1:] if a.startswith('--') and '=' in a)
unknown = set(opts) - {'y0', 'rsmult', 'nseg', 'cachedir', 'plant', 'mutate'}
if unknown or len(args) != 6:
    sys.exit(__doc__ + f"\nunknown options {sorted(unknown)} or wrong argument count")
eps, w, l = int(args[0]), float(args[1]), int(args[2])
x0, x1, y1 = (float(a) for a in args[3:6])
y0 = float(opts.get('y0', 0.0)); rsmult = float(opts.get('rsmult', 40.0)); nseg = int(opts.get('nseg', 10))
cachedir = opts.get('cachedir', HERE); mutate = opts.get('mutate')
plant = complex(*(float(a) for a in opts['plant'].split(','))) if 'plant' in opts else None
assert mutate in (None, 'reverse', 'wrongdouble', 'misplant', 'onesided'), mutate
assert x0 < x1 and 0 <= y0 < y1
sym = (y0 == 0.0)

st = B.solve_star(w, eps)
P = PolarProblem(st, l, rtol=float(os.environ.get('RTOL', '1e-8')))
# fixed scalar start radius: rsmult decay lengths at the slowest-decaying point of the contour
ks = lambda s: (np.sqrt(1 + (s - 1j*w)**2), np.sqrt(1 + (s + 1j*w)**2))
cands = [x0 + 1j*y for y in np.linspace(y0, y1, 200)] + [x + 1j*y1 for x in np.linspace(x0, x1, 200)]
if not sym:
    cands += [x1 + 1j*y for y in np.linspace(y0, y1, 200)] + [x + 1j*y0 for x in np.linspace(x0, x1, 200)]
kmin = min(min(k.real for k in ks(c)) for c in cands)
P.Rs_fixed = P.rm + rsmult / kmin
tag = f"eps{eps:+d}_w{w}_l{l}_{x0}_{x1}_{y1}" + ("" if sym else f"_y0{y0}") + ("" if rsmult == 40.0 else f"_rs{rsmult:g}")
cf = os.path.join(cachedir, f"winding2_cache_{tag}.pkl")
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
if sym:
    corners = [complex(x1, 0), complex(x1, y1), complex(x0, y1), complex(x0, 0)]
else:
    corners = [complex(x1, y0), complex(x1, y1), complex(x0, y1), complex(x0, y0), complex(x1, y0)]
if mutate == 'reverse':
    corners = corners[::-1]

def count(F, verbose):
    """phase accumulated by F along the contour, adaptive to 0.4 rad per step"""
    total, maxstep, nmin = 0.0, 0.0, np.inf
    for a, b in zip(corners[:-1], corners[1:]):
        seg = list(np.linspace(0, 1, nseg + 1)); i = 0
        while i < len(seg) - 1:
            da, db = F(a + (b - a)*seg[i]), F(a + (b - a)*seg[i + 1])
            dph = np.angle(db / da)
            if abs(dph) > 0.4 and seg[i + 1] - seg[i] > 1e-5:
                seg.insert(i + 1, 0.5*(seg[i] + seg[i + 1])); continue
            total += dph; maxstep = max(maxstep, abs(dph)); nmin = min(nmin, abs(da)); i += 1
        if verbose:
            print(f"   edge {a}->{b}: cumulative phase {total:+.4f}  ({len(cache)} evals, {time.time()-T0:.0f}s)", flush=True)
    doubled = sym != (mutate == 'wrongdouble')
    return (2 if doubled else 1) * total / (2*np.pi), maxstep, nmin

N, maxstep, nmin = count(Dc, True)
if sym:
    print(f"   D(x0)={Dc(complex(x0,0)).real:+.3e}  D(x1)={Dc(complex(x1,0)).real:+.3e} (same sign required)")
print(f"RESULT {tag}: zeros in Re s∈[{x0},{x1}], " + (f"|Im s|≤{y1}" if sym else f"Im s∈[{y0},{y1}]") +
      f": winding = {N:+.4f}  (max phase step {maxstep:.2f} rad, min|D| {nmin:.2e})")

if plant is not None:
    s0 = plant
    inside = x0 < s0.real < x1 and ((0 < abs(s0.imag) < y1) if sym else (y0 < s0.imag < y1))
    expected = (2 if sym else 1) * inside
    s0a = s0 + (10 * (x1 - x0) if mutate == 'misplant' else 0)          # the zero the factor really has
    if sym and mutate != 'onesided':
        pf = lambda s: (s - s0a) * (s - np.conj(s0a))
    else:
        pf = lambda s: (s - s0a)
    Np, maxstep_p, nmin_p = count(lambda s: Dc(s) * pf(s), False)
    got = Np - N
    ok = abs(got - expected) < 0.05 and maxstep_p < 0.45
    print(f"PLANTED {tag}: s0 = {s0}, expected {expected} planted zeros inside; count with the factor {Np:+.4f}, "
          f"without {N:+.4f}, difference {got:+.4f}  (max phase step {maxstep_p:.2f} rad)  "
          f"{'PASS' if ok else 'FAIL'}" + (f"  [mutate={mutate}]" if mutate else ""))
    sys.exit(0 if ok else 1)

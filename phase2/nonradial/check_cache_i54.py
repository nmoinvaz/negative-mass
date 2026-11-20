"""
Cross-machine check of the stored winding2.py caches (#54): the evaluations of #7 (win2, x86 Windows) are reused by
the recount on rpi5, so a few stored points per cache are recomputed here with the same D and compared.

  python check_cache_i54.py CACHEDIR [--mutate=scale]

For every symmetric-box cache winding2_cache_eps*_w*_l*_x0_x1_y1.pkl in CACHEDIR, the points at the 10 %, 50 % and
90 % places of the cache sorted by (Im s, Re s) are recomputed with the setup of winding2.py (rtol 1e-8, fixed R_s,
trivial phase removed) and compared with the stored value.  Pass: the largest relative difference is at most 1e-6
(the printed digits of #7's two machines and #5's recomputation, 3.4e-9, agree far inside this).
Exit status 1 when any point differs by more, or no cache is found.

  --mutate=scale  multiplies each recomputed value by 1.001 (a 1e-3 error the check must catch); exit status 1.
"""
import sys, os, re, glob, numpy as np, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opts = dict(a[2:].split('=', 1) for a in sys.argv[1:] if a.startswith('--') and '=' in a)
if len(args) != 1 or set(opts) - {'mutate'} or opts.get('mutate') not in (None, 'scale'):
    sys.exit(__doc__)
TOL, mutate = 1e-6, opts.get('mutate')
pat = re.compile(r"winding2_cache_eps([+-]\d)_w([\d.]+)_l(\d)_([-\d.]+)_([-\d.]+)_([\d.]+)\.pkl$")
files = sorted(f for f in glob.glob(os.path.join(args[0], 'winding2_cache_*.pkl')) if pat.search(f))
bad = not files
for f in files:
    eps, w, l = int(pat.search(f).group(1)), float(pat.search(f).group(2)), int(pat.search(f).group(3))
    x0, x1, y1 = (float(pat.search(f).group(i)) for i in (4, 5, 6))
    cache = pickle.load(open(f, 'rb'))
    st = B.solve_star(w, eps)
    P = PolarProblem(st, l, rtol=1e-8)
    ks = lambda s: (np.sqrt(1 + (s - 1j*w)**2), np.sqrt(1 + (s + 1j*w)**2))
    cands = [x0 + 1j*y for y in np.linspace(0, y1, 200)] + [x + 1j*y1 for x in np.linspace(x0, x1, 200)]
    kmin = min(min(k.real for k in ks(c)) for c in cands)
    P.Rs_fixed = P.rm + 40.0 / kmin
    keys = sorted(cache, key=lambda s: (s.imag, s.real))
    worst = 0.0
    for q in (0.1, 0.5, 0.9):
        s = keys[int(q * (len(keys) - 1))]
        kp, km = ks(s)
        new = P.D(s) * np.exp(-1j * ((kp + km) * (P.Rs_fixed - P.rm) + (kp + km + s) * P.rm).imag)
        if mutate == 'scale': new *= 1.001
        rel = abs(new - cache[s]) / abs(cache[s]); worst = max(worst, rel)
        print(f"   {os.path.basename(f)}: s = {s}, stored {cache[s]:.9e}, recomputed {new:.9e}, relative {rel:.2e}", flush=True)
    ok = worst <= TOL; bad |= not ok
    print(f"cache eps{eps:+d} w{w} l{l} box [{x0},{x1}] |Im s|<={y1}: {len(cache)} stored, worst relative difference "
          f"of 3 recomputed points {worst:.2e} (limit {TOL:g})  {'PASS' if ok else 'FAIL'}"
          + (f"  [mutate={mutate}]" if mutate else ""), flush=True)
sys.exit(1 if bad else 0)

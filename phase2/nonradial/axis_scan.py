"""
|D(s)| along a vertical line s = x + i y below the continuum thresholds, to locate the zeros of D that lie
close to the imaginary axis (quasi-bound modes: both scalar exponents k± are real there, so the evaluation is as
cheap as on the real axis).  Usage:

  python axis_scan.py eps w l x ymin ymax n [--out=FILE.csv] [--rsmult=40]

Writes FILE.csv (y, Re D, Im D, |D|), appending, restart-safe (points already in the file are skipped).  D is the
same PolarProblem.D as winding2.py (rtol 1e-8 by default, RTOL overrides).  A zero within dx of the line shows as
a dip of |D| of width ~dx in y; local minima of |D| are listed at the end.
"""
import sys, os, numpy as np, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem
args = [a for a in sys.argv[1:] if not a.startswith('--')]
opts = dict(a[2:].split('=', 1) for a in sys.argv[1:] if a.startswith('--') and '=' in a)
eps, w, l = int(args[0]), float(args[1]), int(args[2])
x, ymin, ymax, n = float(args[3]), float(args[4]), float(args[5]), int(args[6])
out = opts.get('out', os.path.join(HERE, f"axis_scan_eps{eps:+d}_w{w}_l{l}_x{x:g}.csv"))
P = PolarProblem(B.solve_star(w, eps), l, rtol=float(os.environ.get('RTOL', '1e-8')))
done = {}
if os.path.exists(out):
    for line in open(out).read().splitlines()[1:]:
        v = line.split(','); done[round(float(v[0]), 12)] = (float(v[1]), float(v[2]))
else:
    open(out, 'w').write('y,ReD,ImD\n')
T0 = time.time()
ys = np.linspace(ymin, ymax, n)
for y in ys:
    if round(y, 12) in done: continue
    d = P.D(complex(x, y))
    done[round(y, 12)] = (d.real, d.imag)
    open(out, 'a').write(f"{y:.12e},{d.real:.12e},{d.imag:.12e}\n")
yy = np.array(sorted(done)); dd = np.array([complex(*done[k]) for k in yy]); a = np.abs(dd)
print(f"x = {x}, y in [{ymin}, {ymax}], {len(yy)} points ({time.time()-T0:.0f}s); |D| range {a.min():.3e} .. {a.max():.3e}")
for i in range(1, len(yy) - 1):
    if a[i] < a[i-1] and a[i] < a[i+1]:
        print(f"   local minimum of |D| at y = {yy[i]:.6f}: |D| = {a[i]:.3e} (neighbours {a[i-1]:.3e}, {a[i+1]:.3e})")

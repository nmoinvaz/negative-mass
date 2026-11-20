"""
|D(s)| along a vertical line s = x + i y below the continuum thresholds, to locate the zeros of D that lie
close to the imaginary axis (quasi-bound modes: both scalar exponents k± are real there, so the evaluation is as
cheap as on the real axis).  Usage:

  python axis_scan.py eps w l x ymin ymax n [--out=FILE.csv] [--mini] [--mutate=move]

--mini   the control: the Mini BS I star (as control_winding2.py), where the known zero s = -5.353e-5 + 0.1195191 i
         must show as a dip of |D| on the line (|D| there below 0.5 of the mean at y +- 0.02, within 0.002 of
         y = 0.1195191); exit 1 otherwise.  --mutate=move looks for the dip at y = 0.1395191 instead (must fail).

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
mini = 'mini' in opts or '--mini' in sys.argv
if mini:
    import mini_bs as mb
    st = mb.star_from_shoot(0.08)
    for _ in range(3): st = B.solve_star(0.853087, +1, guess=st, R=st.R)
    P = PolarProblem(st, 2, LAM=0., MU=0., rtol=float(os.environ.get('RTOL', '1e-8')))
else:
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
def ratio(i, h=0.02):
    """|D_i| over the mean of |D| at y +- h (nearest grid points)"""
    j, k = np.argmin(np.abs(yy - (yy[i] - h))), np.argmin(np.abs(yy - (yy[i] + h)))
    return a[i] / (0.5 * (a[j] + a[k]))
for i in range(1, len(yy) - 1):
    if a[i] < a[i-1] and a[i] < a[i+1]:
        print(f"   local minimum of |D| at y = {yy[i]:.6f}: |D| = {a[i]:.3e}, ratio to the mean at y +- 0.02: {ratio(i):.3f}")
rr = np.array([ratio(i) for i in range(len(yy))])
print(f"   smallest ratio of |D| to its mean at y +- 0.02: {rr.min():.3f} at y = {yy[rr.argmin()]:.4f}; "
      f"max |second difference of ln|D||: {np.abs(np.diff(np.log(a), 2)).max():.3f}")
if mini:
    y0 = 0.1195191 + (0.02 if opts.get('mutate') == 'move' else 0)
    i = int(np.argmin(np.abs(yy - y0))); lo = max(0, i - 3); j = lo + int(np.argmin(a[lo:i + 4]))
    good = abs(yy[j] - y0) < 0.002 and ratio(j) < 0.5
    print(f"CONTROL Mini BS I: dip at y = {yy[j]:.5f} (|D| ratio {ratio(j):.3f}), looked for y = {y0:.5f}  {'PASS' if good else 'FAIL'}")
    sys.exit(0 if good else 1)

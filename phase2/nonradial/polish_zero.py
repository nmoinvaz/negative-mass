"""
Polish zeros of D(s) near the imaginary axis by the complex secant method of mini_bs.py, on the same
PolarProblem.D, and report Re s with its spread over integrator tolerances (a spread over fit choices, not an
enclosure).  Usage:

  python polish_zero.py eps w l RE,IM [RE,IM ...] [--rtol=1e-8,1e-10] [--mini] [--tol=1e-12]

--mini   the control: the Mini BS I star (eps = +1, w = 0.853087, Lambda = mu = 0, as control_winding2.py) with the
         guess 0.0,0.1195 (other guesses RE,IM may be given, only the one near Im 0.1195 is checked); it must return s = -5.353e-5 + 0.1195191 i to 3e-7 (its published value
         sigma = 0.1195 - 5e-5 i; the stored mini_bs_output.txt has -5.35e-5), otherwise exit 1.
--mutate=sign   (with --mini) compares with +5.353e-5 instead: the control must fail.
Exit status 1 when a polish does not converge or the control misses.  A zero is reported as damped (Re s < 0) or
growing (Re s > 0) only when |Re s| exceeds ten times the spread over the tolerances.
"""
import sys, os, numpy as np, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem
args = [a for a in sys.argv[1:] if not a.startswith('--')]
opts = dict((a[2:].split('=', 1) + ['1'])[:2] for a in sys.argv[1:] if a.startswith('--'))
mini = 'mini' in opts
rtols = [float(r) for r in opts.get('rtol', '1e-8,1e-10').split(',')]
tol = float(opts.get('tol', 1e-12))
if mini:
    import mini_bs as mb
    st = mb.star_from_shoot(0.08)
    for _ in range(3): st = B.solve_star(0.853087, +1, guess=st, R=st.R)
    mk = dict(LAM=0., MU=0.); l = 2
    guesses = [complex(*(float(v) for v in a.split(','))) for a in args] or [0.0 + 0.1195j]
else:
    eps, w, l = int(args[0]), float(args[1]), int(args[2])
    st = B.solve_star(w, eps); mk = {}
    guesses = [complex(*(float(v) for v in a.split(','))) for a in args[3:]]
ok = True
for g in guesses:
    sols = []
    for rt in rtols:
        P = PolarProblem(st, l, rtol=rt, **mk)
        x0, x1 = g, g + 1e-4
        f0, f1 = P.D(x0), P.D(x1); T0 = time.time()
        for it in range(40):
            x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
            x0, f0, x1, f1 = x1, f1, x2, P.D(x2)
            if abs(x1 - x0) < tol: break
        conv = abs(x1 - x0) < max(tol, 1e-10)
        sols.append(x1); ok &= conv
        print(f"guess {g}: rtol {rt:g}: s = {x1.real:+.12e} {x1.imag:+.12f} i  |D| = {abs(f1):.2e}  "
              f"last step {abs(x1 - x0):.1e}  {it + 1} iterations  {time.time()-T0:.0f}s  {'converged' if conv else 'NOT CONVERGED'}", flush=True)
    re = [s.real for s in sols]; spread = max(re) - min(re)
    s = sols[-1]
    verdict = ('damped' if s.real < 0 else 'growing') if abs(s.real) > 10 * max(spread, 1e-15) else 'sign not resolved'
    print(f"RESULT zero near {g}: s = {s.real:+.6e} {s.imag:+.8f} i, spread of Re s over rtol {spread:.1e}: {verdict}", flush=True)
    if mini and abs(g.imag - 0.1195) < 0.002:
        target = complex(-5.353e-5 if opts.get('mutate') != 'sign' else 5.353e-5, 0.1195191)
        hit = abs(s - target) < 3e-7 and spread < 1e-8
        print(f"CONTROL Mini BS I: {'PASS' if hit else 'FAIL'}")
        ok &= hit
sys.exit(0 if ok else 1)

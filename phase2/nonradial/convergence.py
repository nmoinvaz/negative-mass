"""
Robustness of the matching determinant D(gamma) (real growth rates) for the phantom campaign star:
integrator tolerance (1e-10 vs 1e-12), matching radius (phi_cut 1e-9 -> r_m ~ 24, 1e-12 -> r_m ~ 34),
background R (60 vs 90), and inner start r0 (1e-3 vs 3e-3).  Reports D for each variant (sign and relative
spread), plus the interior / exterior Zerilli log-derivatives q_int, q_ext at r_m whose equality is the
mode condition (q_int - q_ext is a resolution-independent, well-conditioned diagnostic; see kappa_eff.py).
Usage: python convergence.py l
"""
import sys, os, numpy as np, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem

l = int(sys.argv[1]) if len(sys.argv) > 1 else 2
gams = [1e-3, 0.03, 0.3, 1.0, 3.0, 10.0, 40.0, 137.0, 300.0] + ([1000.0] if l >= 3 else [])
st60 = B.solve_star(0.7603, -1)
st90 = B.solve_star(0.7603, -1, guess=st60, R=90.0)
variants = {
    'base (rtol1e-10, phi_cut1e-9, R60)': dict(star=st60),
    'rtol 1e-12': dict(star=st60, rtol=1e-12),
    'phi_cut 1e-12 (larger r_m)': dict(star=st60, phi_cut=1e-12),
    'background R=90': dict(star=st90),
    'r0 = 3e-3': dict(star=st60, r0=3e-3),
}
probs = {}
for k, v in variants.items():
    v = dict(v); st = v.pop('star')
    probs[k] = PolarProblem(st, l, **v)
    print(f"{k:38s}: r_m={probs[k].rm:.2f}")
print()
for g in gams:
    T = time.time()
    Ds = {k: P.D(g) for k, P in probs.items()}
    base = probs['base (rtol1e-10, phi_cut1e-9, R60)']
    zz, _ = base.interior_zerilli(g)
    f = 1 - 2 * base.Mext / base.rm
    qi = f * zz[1] / zz[0]; qe = base.zerilli_q(g)
    signs = set(np.sign(list(Ds.values())))
    same_rm = [Ds[k] for k in Ds if 'phi_cut' not in k]
    spread = (max(same_rm) - min(same_rm)) / abs(np.mean(same_rm))
    print(f"gamma={g:8.3g}: D_base={Ds['base (rtol1e-10, phi_cut1e-9, R60)']:+.6e}  signs={sorted(signs)}  "
          f"rel.spread (same r_m)={spread:.1e}  D(r_m~34)={Ds['phi_cut 1e-12 (larger r_m)']:+.3e}  "
          f"q_int={qi:+.6e} q_ext={qe:+.6e}  ({time.time()-T:.0f}s)", flush=True)

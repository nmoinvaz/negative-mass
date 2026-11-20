"""
Robustness of the matching determinant D(gamma) (real growth rates) for the phantom campaign star:
integrator tolerance (1e-10 vs 1e-12), matching radius (phi_cut 1e-9 -> r_m ~ 24, 1e-12 -> r_m ~ 34),
background R (60 vs 90), and inner start r0 (1e-3 vs 3e-3).  Reports D for each variant (sign and relative
spread), plus the interior / exterior Zerilli log-derivatives q_int, q_ext at r_m whose equality is the
mode condition (q_int - q_ext is a resolution-independent, well-conditioned diagnostic; see kappa_eff.py).
Usage: python convergence.py l [w]
GAMS (a comma list) replaces the default gammas, NPROC (default 1) runs the (variant, gamma) evaluations in a worker
pool, and each evaluation is appended to convergence_l{l}_w{w}.partial.csv as it finishes, so a stopped run resumes.
"""
import sys, os, json, numpy as np, time
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem

l = int(sys.argv[1]) if len(sys.argv) > 1 else 2
w = float(sys.argv[2]) if len(sys.argv) > 2 else 0.7603
gams = [1e-3, 0.03, 0.3, 1.0, 3.0, 10.0, 40.0, 137.0, 300.0] + ([1000.0] if l >= 3 else [])
if os.environ.get('GAMS'):
    gams = [float(x) for x in os.environ['GAMS'].split(',')]
BASE = 'base (rtol1e-10, phi_cut1e-9, R60)'
VARIANTS = {
    BASE: dict(star='60'),
    'rtol 1e-12': dict(star='60', rtol=1e-12),
    'phi_cut 1e-12 (larger r_m)': dict(star='60', phi_cut=1e-12),
    'background R=90': dict(star='90'),
    'r0 = 3e-3': dict(star='60', r0=3e-3),
}
probs = {}


def init():
    st60 = B.solve_star(w, -1)
    st90 = B.solve_star(w, -1, guess=st60, R=90.0)
    for k, v in VARIANTS.items():
        v = dict(v)
        st = st60 if v.pop('star') == '60' else st90
        probs[k] = PolarProblem(st, l, **v)


def task(kg):
    k, g = kg
    T = time.time()
    if k == 'q':
        base = probs[BASE]
        zz, _ = base.interior_zerilli(g)
        f = 1 - 2 * base.Mext / base.rm
        return k, g, [f * zz[1] / zz[0], base.zerilli_q(g)], time.time() - T
    return k, g, probs[k].D(g), time.time() - T


if __name__ == "__main__":
    init()
    for k in VARIANTS:
        print(f"{k:38s}: r_m={probs[k].rm:.2f}")
    print()
    part = os.path.join(HERE, f"convergence_l{l}_w{w}.partial.csv")
    res = {}
    if os.path.exists(part):
        for line in open(part):
            k, g, v, dt = json.loads(line)
            res[(k, g)] = (v, dt)
    tasks = [(k, g) for g in sorted(gams, reverse=True) for k in list(VARIANTS) + ['q'] if (k, g) not in res]
    with Pool(int(os.environ.get('NPROC', '1')), initializer=init) as pool, open(part, 'a') as fp:
        for k, g, v, dt in pool.imap_unordered(task, tasks, chunksize=1):
            res[(k, g)] = (v, dt)
            fp.write(json.dumps([k, g, v, dt]) + "\n")
            fp.flush()
    for g in gams:
        Ds = {k: res[(k, g)][0] for k in VARIANTS}
        qi, qe = res[('q', g)][0]
        T = sum(res[(k, g)][1] for k in list(VARIANTS) + ['q'])
        signs = set(np.sign(list(Ds.values())))
        same_rm = [Ds[k] for k in Ds if 'phi_cut' not in k]
        spread = (max(same_rm) - min(same_rm)) / abs(np.mean(same_rm))
        print(f"gamma={g:8.3g}: D_base={Ds[BASE]:+.6e}  signs={sorted(signs)}  "
              f"rel.spread (same r_m)={spread:.1e}  D(r_m~34)={Ds['phi_cut 1e-12 (larger r_m)']:+.3e}  "
              f"q_int={qi:+.6e} q_ext={qe:+.6e}  ({T:.0f}s)", flush=True)

"""
Refinement of the real-gamma scan of scan_growing.py until halving the step changes no sign count.

Level 0 is the stored scan: logspace(gmin, gmax, n0) as written by scan_growing.py (same P.D, rtol 1e-10,
phi_cut 1e-9).  Level k adds the log-midpoints of the level k-1 grid, so every earlier point is kept and the step
is halved.  D at a new point comes from the same PolarProblem.D(gamma) as the scan.  Each evaluation is appended to
a csv (restart-safe).  Usage:

  python refine_gamma.py eps w l gmin gmax n0 levels [--csv=SCAN.csv] [--out=DIR] [--check=3] [--procs=1]
                         [--plant=A1,A2] [--mutate=NAME]

Checks (exit 1 when one fails):
  reproduce   D recomputed at --check stored scan points agrees with the stored value to a relative 1e-6.
  halving     the number of sign changes of D is the same at the last two levels (and at the last three, when
              levels >= 2); a change is reported with the bracket and the brentq root.
  planted     with --plant=A1,A2 the factor (g-A1)(g-A2) is put in front of D, the pair is placed to straddle a
              level-1 midpoint, and the check requires that level 0 sees no sign change, and that the finest level
              sees exactly two, in the intervals that contain A1 and A2.
Mutations (the check named must fail):
  --mutate=coarse  stop after level 0 (with --plant: the pair is missed, the planted check fails)
  --mutate=pair    insert an unannounced close pair into D (the halving check fails)
  --mutate=stale   compare the recomputed values with the stored ones scaled by 1.01 (the reproduce check fails)
"""
import sys, os, time, numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opts = dict(a[2:].split('=', 1) for a in sys.argv[1:] if a.startswith('--') and '=' in a)
if set(opts) - {'csv', 'out', 'check', 'procs', 'plant', 'mutate'} or len(args) != 7:
    sys.exit(__doc__)
eps, w, l = int(args[0]), float(args[1]), int(args[2])
gmin, gmax, n0, levels = float(args[3]), float(args[4]), int(args[5]), int(args[6])
mutate = opts.get('mutate'); assert mutate in (None, 'coarse', 'pair', 'stale'), mutate
procs = int(opts.get('procs', 1)); ncheck = int(opts.get('check', 3))
plant = tuple(float(a) for a in opts['plant'].split(',')) if 'plant' in opts else None
stem = f"eps{eps:+d}_w{w}_l{l}"
scan_csv = opts.get('csv', os.path.join(HERE, f"scan_{stem}.csv"))
outdir = opts.get('out', HERE)
cache_csv = os.path.join(outdir, f"refine_{stem}_{gmin:g}_{gmax:g}_{n0}.csv")
rtol = float(os.environ.get('RTOL', '1e-10')); phi_cut = float(os.environ.get('PHICUT', '1e-9'))

_P = None
def init():
    global _P
    _P = PolarProblem(B.solve_star(w, eps), l, rtol=rtol, phi_cut=phi_cut)
def Dg(g):
    return _P.D(g)

def grid(level):
    """log grid of the given level: 2^level (n0 - 1) + 1 points of logspace(gmin, gmax)"""
    return np.logspace(np.log10(gmin), np.log10(gmax), (n0 - 1) * 2**level + 1)

def keyof(g):
    """integer position of g on the log grid in units of 2^-20 of a level-0 step: equal for equal grid points"""
    return int(round((np.log10(g) - np.log10(gmin)) / (np.log10(gmax) - np.log10(gmin)) * (n0 - 1) * 2**20))

def signs_changes(g, d):
    return [i for i in range(len(g) - 1) if np.sign(d[i]) != np.sign(d[i + 1])]

if __name__ == "__main__":
    T0 = time.time()
    st = np.loadtxt(scan_csv, delimiter=',', skiprows=1)
    assert len(st) == n0 and np.allclose(st[:, 0], grid(0), rtol=1e-12), "stored scan is not logspace(gmin, gmax, n0)"
    store = {keyof(g): d for g, d in st}                       # level 0 and every refined value
    if os.path.exists(cache_csv):
        for line in open(cache_csv).read().splitlines()[1:]:
            a, b = line.split(',')
            store[keyof(float(a))] = float(b)
    else:
        with open(cache_csv, 'w') as f: f.write('gamma,D\n')
    ok = True

    # ---- reproduce check: same P.D as the scan, at a few stored points
    if ncheck:
        init()
        idx = np.unique(np.linspace(0, n0 - 1, ncheck).astype(int))
        for i in idx:
            d = Dg(st[i, 0]); ref = st[i, 1] * (1.01 if mutate == 'stale' else 1.0)
            rel = abs(d - ref) / abs(ref); good = rel < 1e-6; ok &= good
            print(f"reproduce gamma={st[i,0]:.6g}: D = {d:+.9e}, stored {ref:+.9e}, relative {rel:.1e}  "
                  f"{'PASS' if good else 'FAIL'}" + (f"  [mutate={mutate}]" if mutate else ""), flush=True)

    # ---- evaluate the missing points of each level
    nlev = 0 if mutate == 'coarse' else levels
    for k in range(1, nlev + 1):
        need = [g for g in grid(k) if keyof(g) not in store]
        print(f"level {k}: {len(grid(k))} points, {len(need)} to evaluate  ({time.time()-T0:.0f}s)", flush=True)
        # sorted by gamma so a partial run covers a contiguous low-gamma range; chunks of 2*procs written as they finish
        if need:
            pool = Pool(procs, initializer=init) if procs > 1 else None
            if pool is None: init()
            for j in range(0, len(need), max(1, 2 * procs)):
                part = need[j:j + 2 * procs]
                vals = pool.map(Dg, part) if pool else [Dg(g) for g in part]
                with open(cache_csv, 'a') as f:
                    for g, d in zip(part, vals):
                        store[keyof(g)] = d; f.write(f"{g:.17e},{d:.17e}\n")
                if (j // max(1, 2 * procs)) % 10 == 0:
                    print(f"   level {k}: {j + len(part)}/{len(need)} done, gamma up to {part[-1]:.4g}  ({time.time()-T0:.0f}s)", flush=True)
            if pool: pool.close()

    # ---- sign counts by level (unplanted data, optionally with a mutated or planted factor)
    def factor(g):
        f = np.ones_like(g)
        if mutate == 'pair':                                      # unannounced close pair inside one level-0 interval
            m = np.sqrt(grid(0)[n0 // 2] * grid(0)[n0 // 2 + 1]); f = f * (g - m * 0.99) * (g - m * 1.01)
        if plant:
            f = f * (g - plant[0]) * (g - plant[1])
        return f
    counts, mind = [], []
    for k in range(0, nlev + 1):
        g = grid(k)
        d = np.array([store[keyof(x)] for x in g]) * factor(g)
        ch = signs_changes(g, d); counts.append(len(ch))
        # dip diagnostic: |D_mid| / sqrt(|D_a D_b|) over each interval of the level (1 for log-linear D)
        rr = np.abs(d[1:-1]) / np.sqrt(np.abs(d[:-2] * d[2:])) if k > 0 else np.array([np.nan])
        print(f"level {k}: step ratio {g[1]/g[0]:.5f}, {len(g)} points, sign changes {len(ch)}"
              + ("" if k == 0 else f", min |D_i|/sqrt|D_(i-1) D_(i+1)| {np.nanmin(rr):.3f}, max {np.nanmax(rr):.3f}")
              + f", min|D| {np.abs(d).min():.3e}" + (f"  brackets {[(float(f'{g[i]:.5g}'), float(f'{g[i+1]:.5g}')) for i in ch]}" if ch else ""),
              flush=True)
    if plant:
        a1, a2 = sorted(plant); g = grid(nlev); d = np.array([store[keyof(x)] for x in g]) * factor(g)
        ch = signs_changes(g, d)
        found = len(ch) == 2 and g[ch[0]] < a1 < g[ch[0] + 1] and g[ch[1]] < a2 < g[ch[1] + 1] and counts[0] == 0
        print(f"PLANTED pair at {a1:.6g}, {a2:.6g}: level 0 sees {counts[0]} sign changes, finest level sees {counts[-1]}"
              f", intervals {[(float(f'{g[i]:.6g}'), float(f'{g[i+1]:.6g}')) for i in ch]}  {'PASS' if found else 'FAIL'}"
              + (f"  [mutate={mutate}]" if mutate else ""))
        ok &= found
    else:
        same = counts[-1] == counts[-2] and (len(counts) < 3 or counts[-2] == counts[-3]) if len(counts) > 1 else False
        print(f"HALVING sign counts by level {counts}: halving the step changes no sign count  {'PASS' if same else 'FAIL'}"
              + (f"  [mutate={mutate}]" if mutate else ""))
        ok &= same
    print(f"{'ALL CHECKS PASS' if ok else 'A CHECK FAILED'}  ({time.time()-T0:.0f}s)")
    sys.exit(0 if ok else 1)

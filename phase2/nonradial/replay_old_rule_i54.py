"""
Replays the FIRST step rule of winding2.py (git 81cc13d, the one #7 ran) on the stored evaluations of a box, and
finds the steps it accepted that lost 2 pi n (#54).  No evaluation of D is made.  Usage, from this directory:

  python replay_old_rule_i54.py ../../../negative-mass-raw [--only=w085_l3,...] [--mutate=thr]

For each of #7's four boxes (Re s in [0.1, 2], |Im s| <= 2.5), with I54-rpi5/cache holding #7's evaluations plus the
ones the current rule asked for:
  1. the old rule, replayed on the cache, must reproduce #7's logged result (I7-win2/wind_*.log): winding, max
     phase step, min|D|, and it must ask for exactly the points of #7's original cache (I7-win2/winding2_cache_*.pkl);
  2. the cached points on each edge, sorted along the edge, form a dense chain whose steps are all small; its total
     phase is the true count.  It must equal the count of the current rule (I54-rpi5/<job>.log);
  3. for each step the old rule accepted, true change - measured angle must be a multiple of 2 pi (within 0.05).
     Each step with a non-zero multiple is printed: this is where 2 pi was lost.
Exit status 1 when a check fails or a run is missing.

  --mutate=thr   the old rule's threshold 0.4 -> 0.3: it asks for other points and no longer matches (must exit 1).
"""
import re, sys, os, pickle, numpy as np

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opts = dict((a[2:].split('=', 1) + ['1'])[:2] for a in sys.argv[1:] if a.startswith('--'))
if len(args) != 1 or set(opts) - {'only', 'mutate'} or opts.get('mutate') not in (None, 'thr'):
    sys.exit(__doc__)
ROOT, THR = args[0], (0.3 if opts.get('mutate') == 'thr' else 0.4)
mutate = opts.get('mutate')
BOXES = {  # job: (tag, I7 log)
    'w085_l3': ('eps-1_w0.85_l3_0.1_2.0_2.5', 'wind_085_l3.log'), 'w085_l2': ('eps-1_w0.85_l2_0.1_2.0_2.5', 'wind_085_l2.log'),
    'w070_l3': ('eps-1_w0.7_l3_0.1_2.0_2.5', 'wind_070_l3.log'), 'w070_l2': ('eps-1_w0.7_l2_0.1_2.0_2.5', 'wind_070_l2.log')}
jobs = opts['only'].split(',') if 'only' in opts else list(BOXES)
res = re.compile(r"RESULT (\S+): .*winding = ([+-]\d+\.\d+)  \(max phase step ([\d.]+) rad, min\|D\| (\S+)\)")
x0, x1, y1 = 0.1, 2.0, 2.5
corners = [complex(x1, 0), complex(x1, y1), complex(x0, y1), complex(x0, 0)]
names = ['right', 'top', 'left']
rk = lambda s: complex(round(s.real, 12), round(s.imag, 12))
bad = False

def logged(path, tag):
    m = [x for x in res.finditer(open(path).read()) if x.group(1) == tag] if os.path.exists(path) else []
    return (float(m[0].group(2)), m[0].group(3), m[0].group(4)) if m else None

for job in jobs:
    tag, oldlog = BOXES[job]
    cache = pickle.load(open(os.path.join(ROOT, 'I54-rpi5', 'cache', f'winding2_cache_{tag}.pkl'), 'rb'))
    orig = pickle.load(open(os.path.join(ROOT, 'I7-win2', f'winding2_cache_{tag}.pkl'), 'rb'))
    old, new = logged(os.path.join(ROOT, 'I7-win2', oldlog), tag), logged(os.path.join(ROOT, 'I54-rpi5', job + '.log'), tag)
    if old is None or new is None:
        print(f"{tag}: log missing (I7 {old}, I54 {new})  FAIL"); bad = True; continue
    asked = set()
    def Dc(s):
        asked.add(rk(s)); return cache[rk(s)]
    class Missing(Exception): pass
    # the first rule, as in 81cc13d
    total, maxstep, nmin, steps = 0.0, 0.0, np.inf, []
    try:
        for e, (a, b) in enumerate(zip(corners[:-1], corners[1:])):
            seg = list(np.linspace(0, 1, 11)); i = 0
            while i < len(seg) - 1:
                pa, pb = a + (b - a)*seg[i], a + (b - a)*seg[i + 1]
                da, db = Dc(pa), Dc(pb)
                dph = np.angle(db / da)
                if abs(dph) > THR and seg[i + 1] - seg[i] > 1e-5:
                    seg.insert(i + 1, 0.5*(seg[i] + seg[i + 1])); continue
                total += dph; maxstep = max(maxstep, abs(dph)); nmin = min(nmin, abs(da)); steps.append((e, pa, pb, dph)); i += 1
    except KeyError as err:
        print(f"{tag}: the old rule asks for a point that is not stored: {err}  FAIL" + (f"  [mutate={mutate}]" if mutate else ""))
        bad = True; continue
    wind = 2 * total / (2 * np.pi)
    same_pts = asked == set(orig)
    ok1 = abs(wind - old[0]) < 5e-5 and f"{maxstep:.2f}" == old[1] and f"{nmin:.2e}" == old[2] and same_pts
    print(f"{tag}: old rule replayed: winding {wind:+.4f}, max step {maxstep:.2f}, min|D| {nmin:.2e}, {len(asked)} points "
          f"(I7's cache {len(orig)}, same set {same_pts}); I7 log {old[0]:+.4f}, {old[1]}, {old[2]}  {'MATCH' if ok1 else 'FAIL'}"
          + (f"  [mutate={mutate}]" if mutate else ""))
    # dense chain along each edge: all cached points on it, sorted along it
    phi, ctot, cmax, ok3, chains = {}, 0.0, 0.0, True, {}
    for e, (a, b) in enumerate(zip(corners[:-1], corners[1:])):
        pts = []
        for s in cache:
            t = (s - a) / (b - a)
            if abs(t.imag) < 1e-9 and -1e-9 <= t.real <= 1 + 1e-9: pts.append((t.real, s))
        pts.sort(key=lambda x: x[0]); chains[e] = pts
        acc = ctot
        for (t0, s0), (t1, s1) in zip(pts[:-1], pts[1:]):
            st = np.angle(cache[s1] / cache[s0]); cmax = max(cmax, abs(st)); acc += st
            phi[(e, s1)] = acc
        phi[(e, pts[0][1])] = ctot; ctot = acc
    cwind = 2 * ctot / (2 * np.pi)
    ok2 = abs(cwind - new[0]) < 0.05 and cmax <= 0.45
    print(f"   dense chain of cached points on the edges: {len(cache)} points, max step {cmax:.2f}, winding {cwind:+.4f}; "
          f"current rule's log {new[0]:+.4f}  {'MATCH' if ok2 else 'FAIL'}")
    lost = []
    for e, pa, pb, dph in steps:
        tru = phi[(e, rk(pb))] - phi[(e, rk(pa))]
        k = (tru - dph) / (2 * np.pi)
        if abs(k - round(k)) * 2 * np.pi > 0.05: ok3 = False
        if round(k) != 0: lost.append((names[e], pa, pb, dph, tru, round(k), e))
    for nm, pa, pb, dph, tru, k, e in lost:
        print(f"   LOST: {nm} edge, s = {pa:.6f} -> {pb:.6f}: measured angle {dph:+.4f}, true change {tru:+.4f} ({k:+d} x 2 pi)")
        a, b = corners[e], corners[e + 1]
        ta, tb = sorted([((pa - a) / (b - a)).real, ((pb - a) / (b - a)).real])
        inside = [s for t, s in chains[e] if ta - 1e-9 <= t <= tb + 1e-9]
        sts = [np.angle(cache[s1] / cache[s0]) for s0, s1 in zip(inside[:-1], inside[1:])]
        print(f"      {len(inside)} stored points on it, |D| in [{min(abs(cache[s]) for s in inside):.2e}, "
              f"{max(abs(cache[s]) for s in inside):.2e}], step angles {min(sts):+.3f} to {max(sts):+.3f} rad "
              f"({'all one sign' if min(sts) * max(sts) > 0 else 'mixed signs'})")
    if not lost: print("   no step of the old rule lost a multiple of 2 pi")
    ok4 = ok3 and abs(wind + 2 * sum(l_[5] for l_ in lost) - cwind) < 0.05
    print(f"   old winding {wind:+.4f} + 2 x (lost multiples of 2 pi) {2 * sum(l_[5] for l_ in lost):+d} = dense-chain winding "
          f"{cwind:+.4f}  {'consistent' if ok4 else 'FAIL'}")
    bad |= not (ok1 and ok2 and ok4)
print("RESULT: " + ("FAIL" if bad else "every replay matches #7's log and every difference is a multiple of 2 pi")
      + (f"  [mutate={mutate}]" if mutate else ""))
sys.exit(1 if bad else 0)

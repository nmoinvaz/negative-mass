"""
Replay of winding2.py's step rule on a box's cache, to find a step that lost 2 pi (#115 s1 at l = 5, #126). Reads only the
pickle cache of the symmetric box, no evaluation of D. Usage:

  python diagnose_plant.py CACHE_PKL x0 x1 y1 RE,IM [--guard=1|all] [--expect=none|MEASURED,CHAIN] [--mutate=NAME]

Replays the adaptive counter of winding2.py (copied, same rule) on the cached values of D and of the planted count
D (s - s0)(s - conj s0), then compares each accepted step with the dense chain of every cached point on the edge, sorted
along it. A step whose measured phase differs from the chain's over the same stretch by more than 0.5 rad lost or gained a
multiple of 2 pi. The chain shares D and its cache with the count, so it tests the step rule and not D.
  --guard=1       deepest bisection level with the half-step test, as winding2.py's --guard (`all` or an integer).
  --expect=...    makes the replay a check that exits 1 when it fails. `none`: no step differs from the chain.
                  MEASURED,CHAIN (rad): exactly one step differs, in D times the factor, with these phases to 4
                  decimals. Both also need the chain's largest increment below 1 rad, so that the chain can be trusted.
  --mutate=NAME   sparse: the chain keeps only the end points of the accepted steps, so it loses the 2 pi as the step does;
                  wrongexpect: MEASURED and CHAIN are shifted by 1e-3.
A point the counter asks for and the cache lacks stops the replay with exit 1: the cache is not the one of the run.
Without --expect the output is a diagnostic and the exit status is 0.
"""
import sys, pickle, numpy as np

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opts = dict(a[2:].split('=', 1) for a in sys.argv[1:] if a.startswith('--') and '=' in a)
if len(args) != 5 or set(opts) - {'guard', 'expect', 'mutate'}:
    sys.exit(__doc__)
cf, (x0, x1, y1) = args[0], (float(a) for a in args[1:4])
s0 = complex(*(float(a) for a in args[4].split(',')))
guard = np.inf if opts.get('guard', '1') == 'all' else int(opts.get('guard', '1'))
expect, mutate = opts.get('expect'), opts.get('mutate')
assert mutate in (None, 'sparse', 'wrongexpect'), mutate
if mutate == 'wrongexpect' and (expect is None or expect == 'none'):
    sys.exit("--mutate=wrongexpect needs --expect=MEASURED,CHAIN")
cache = pickle.load(open(cf, 'rb'))
corners = [complex(x1, 0), complex(x1, y1), complex(x0, y1), complex(x0, 0)]
key = lambda s: complex(round(s.real, 12), round(s.imag, 12))


class Missing(Exception):
    pass


def D(s):
    if key(s) not in cache:
        raise Missing(key(s))
    return cache[key(s)]


pf = lambda s: (s - s0) * (s - np.conj(s0))
F = {'D': D, 'D*factor': lambda s: D(s) * pf(s)}


def count(f, a, b, nseg=10):
    """accepted leaf steps (t0, t1, measured phase) of winding2.count on one edge"""
    out = []
    def step(t0, t1, depth):
        da, db = f(a + (b - a) * t0), f(a + (b - a) * t1)
        dph = np.angle(db / da)
        fine = t1 - t0 <= 1e-5
        if abs(dph) <= 0.4 and depth <= guard and not fine:
            dm = f(a + (b - a) * 0.5 * (t0 + t1))
            h1, h2 = np.angle(dm / da), np.angle(db / dm)
            if max(abs(h1), abs(h2)) > 0.4 or abs(h1 + h2 - dph) > 0.1:
                step(t0, 0.5 * (t0 + t1), depth + 1); step(0.5 * (t0 + t1), t1, depth + 1); return
        if abs(dph) > 0.4 and not fine:
            step(t0, 0.5 * (t0 + t1), depth + 1); step(0.5 * (t0 + t1), t1, depth + 1); return
        out.append((t0, t1, dph))
    ts = np.linspace(0, 1, nseg + 1)
    for t0, t1 in zip(ts[:-1], ts[1:]):
        step(t0, t1, 0)
    return out


def chain(f, a, b, leaves):
    """all cached points on the edge a->b sorted along it, with the principal increments between neighbours"""
    pts = []
    for s in cache:
        t = ((s - a) / (b - a))
        if abs(t.imag) < 1e-9 and -1e-9 <= t.real <= 1 + 1e-9:
            if mutate != 'sparse' or any(abs(t.real - e) < 1e-9 for l in leaves for e in l[:2]):
                pts.append((t.real, s))
    pts.sort()
    v = [f(s) for _, s in pts]
    inc = [np.angle(v[i + 1] / v[i]) for i in range(len(v) - 1)]
    return [p[0] for p in pts], inc


print(f"cache {cf}: {len(cache)} points; planted s0 = {s0}; guard {opts.get('guard', '1')}"
      + (f"; mutate={mutate}" if mutate else ""))
lost, chainmax = [], 0.0
try:
    for name, f in F.items():
        tot_ad, tot_ch = 0.0, 0.0
        print(f"== {name}")
        for a, b in zip(corners[:-1], corners[1:]):
            leaves = count(f, a, b)
            ad = sum(x[2] for x in leaves)
            ts, inc = chain(f, a, b, leaves)
            ch = sum(inc)
            tot_ad += ad; tot_ch += ch; chainmax = max(chainmax, max(abs(x) for x in inc))
            print(f"  edge {a}->{b}: adaptive {ad:+.4f} ({len(leaves)} steps, max {max(abs(x[2]) for x in leaves):.3f}), "
                  f"dense chain {ch:+.4f} ({len(inc)} steps, max {max(abs(x) for x in inc):.3f}), "
                  f"difference/2pi {(ad - ch) / (2 * np.pi):+.4f}")
            # leaf steps whose measured phase differs from the dense chain over the same stretch by a multiple of 2 pi
            for t0, t1, dph in leaves:
                sel = [i for i in range(len(inc)) if ts[i] >= t0 - 1e-9 and ts[i + 1] <= t1 + 1e-9]
                if sel:
                    c = sum(inc[i] for i in sel)
                    if abs(dph - c) > 0.5:
                        lost.append((name, dph, c))
                        print(f"    step t in [{t0:.5f}, {t1:.5f}] (s = {a + (b - a) * t0:.5f} .. {a + (b - a) * t1:.5f}): "
                              f"measured {dph:+.4f}, dense chain {c:+.4f} over {len(sel)} cached steps "
                              f"(max {max(abs(inc[i]) for i in sel):.3f}), difference/2pi {(dph - c) / (2 * np.pi):+.4f}")
        print(f"  total adaptive {tot_ad:+.4f} -> count {2 * tot_ad / (2 * np.pi):+.4f};  dense chain {tot_ch:+.4f} -> "
              f"count {2 * tot_ch / (2 * np.pi):+.4f}")
except Missing as e:
    print(f"CHECK FAIL: the cache lacks the point {e.args[0]} that the counter asks for (guard {opts.get('guard', '1')})")
    sys.exit(1)

if expect is None:
    print("no --expect: diagnostic only")
    sys.exit(0)
why = []
if chainmax >= 1.0:
    why.append(f"the dense chain's largest increment {chainmax:.3f} rad is not below 1")
if expect == 'none':
    want = 'no step differs from the chain'
    if lost:
        why.append(f"{len(lost)} step(s) differ from the chain")
else:
    m, c = (float(v) for v in expect.split(','))
    if mutate == 'wrongexpect':
        m, c = m + 1e-3, c + 1e-3
    want = f"one step, in D*factor, measured {m:+.4f} against the chain's {c:+.4f}"
    if len(lost) != 1 or lost[0][0] != 'D*factor' or abs(lost[0][1] - m) > 5e-5 + 1e-9 or abs(lost[0][2] - c) > 5e-5 + 1e-9:
        why.append("found " + (", ".join(f"{n}: measured {d:+.4f}, chain {k:+.4f}" for n, d, k in lost) or "no step that differs"))
print(f"CHECK expected {want}; chain's largest increment {chainmax:.3f} rad  " + ("PASS" if not why else "FAIL: " + "; ".join(why))
      + (f"  [mutate={mutate}]" if mutate else ""))
sys.exit(1 if why else 0)

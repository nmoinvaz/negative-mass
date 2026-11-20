"""
Collects the recount of #54 (stage 1) from the raw logs into one table, old count next to new, so that every number
quoted in RESULTS.md comes from a file in this repository.  Usage, from this directory:

  python collect_i54.py ../../../negative-mass-raw [--mutate=NAME] > i54_recount_output.txt

Reads the RESULT and PLANTED lines and the per-edge phase lines of the winding2.py logs:
  new   I54-rpi5/<job>.log (the current step rule), and I5-rpi5/rep_l2.log for phantom omega = 0.7603, l = 2 (#5);
  old   I7-win2/wind_0{70,85}_l{2,3}.log (#7) and winding2_output.txt (this directory), both made with the first rule.
Exit status 1 when a new run has no result or no PLANTED line, when its planted zero is not returned or its max phase
step exceeds 0.45 rad, when its winding is not within 0.05 of an integer, or when an old value is missing.

  --mutate=noplant      ignores the PLANTED lines (every row must then fail),
  --mutate=wrongexpect  takes one planted zero as the expectation for every row, instead of the logged one (it must fail).
"""
import re, sys, os

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opts = dict((a[2:].split('=', 1) + ['1'])[:2] for a in sys.argv[1:] if a.startswith('--'))
mutate = opts.get('mutate')
if len(args) != 1 or set(opts) - {'mutate'} or mutate not in (None, 'noplant', 'wrongexpect'):
    sys.exit(__doc__)
ROOT = args[0]
HERE = os.path.dirname(os.path.abspath(__file__))
NEW, I5, I7 = (os.path.join(ROOT, d) for d in ('I54-rpi5', 'I5-rpi5', 'I7-win2'))
OLD0 = os.path.join(HERE, 'winding2_output.txt')

res = re.compile(r"RESULT (\S+): zeros in Re s∈\[(\S+),(\S+)\], \|Im s\|≤(\S+): "
                 r"winding = ([+-]\d+\.\d+)  \(max phase step ([\d.]+) rad, min\|D\| (\S+)\)")
pla = re.compile(r"PLANTED (\S+): s0 = \((\S+)\), expected (\d) planted zeros inside; count with the factor ([+-]\d+\.\d+), "
                 r"without ([+-]\d+\.\d+), difference ([+-]\d+\.\d+)  \(max phase step ([\d.]+) rad\)  (PASS|FAIL)")
edge = re.compile(r"edge \(\S+\)->\(\S+\): cumulative phase ([+-]\d+\.\d+)  \((\d+) evals")
# (star, l, tag, new log, old log).  The old log of the first three comes from winding2_output.txt.
BOX = "0.1_2.0_2.5"
ROWS = [
    ('phantom 0.85', 3, f"eps-1_w0.85_l3_{BOX}", 'w085_l3', os.path.join(I7, 'wind_085_l3.log')),
    ('phantom 0.85', 2, f"eps-1_w0.85_l2_{BOX}", 'w085_l2', os.path.join(I7, 'wind_085_l2.log')),
    ('phantom 0.70', 3, f"eps-1_w0.7_l3_{BOX}", 'w070_l3', os.path.join(I7, 'wind_070_l3.log')),
    ('phantom 0.70', 2, f"eps-1_w0.7_l2_{BOX}", 'w070_l2', os.path.join(I7, 'wind_070_l2.log')),
    ('phantom 0.7603', 3, f"eps-1_w0.7603_l3_{BOX}", 'w07603_l3', OLD0),
    ('canonical 0.75', 2, f"eps+1_w0.75_l2_{BOX}", 'c075_l2', OLD0),
    ('phantom 0.7603', 2, f"eps-1_w0.7603_l2_{BOX}", None, OLD0),         # new: #5's rep_l2
]

def read(path):
    return open(path).read() if os.path.exists(path) else ''

def result(txt, tag):
    """the RESULT line of tag, with the three edge phases printed before it (cumulative -> increments)"""
    for m in res.finditer(txt):
        if m.group(1) == tag:
            prev = max([x.end() for x in res.finditer(txt) if x.end() <= m.start()], default=0)
            cum = [float(x.group(1)) for x in edge.finditer(txt[prev:m.start()])][-3:]
            inc = [b - a for a, b in zip([0.0] + cum[:-1], cum)] if len(cum) == 3 else None
            return dict(wind=float(m.group(5)), step=m.group(6), mind=m.group(7), edges=inc)
    return None

bad = False
rows, edges_out = [], []
for star, l, tag, job, oldlog in ROWS:
    old = result(read(oldlog), tag)
    if job is None:
        newtxt = read(os.path.join(I5, 'rep_l2.log')).split('## mutate=reverse')[0]       # the unmutated run
        src = 'I5-rpi5/rep_l2.log'
    else:
        newtxt, src = read(os.path.join(NEW, job + '.log')), f'I54-rpi5/{job}.log'
    new = result(newtxt, tag)
    p = None if mutate == 'noplant' else next((m for m in pla.finditer(newtxt) if m.group(1) == tag), None)
    if old is None:
        bad = True
    if new is None:
        rows.append(f"| {star} | {l} | {old['wind']:+.4f} | {old['step']} | {old['mind']} | not run | | | | | | | {src} |" if old else
                    f"| {star} | {l} | missing | | | not run | | | | | | | {src} |")
        bad = True; continue
    ok = abs(new['wind'] - round(new['wind'])) < 0.05 and float(new['step']) <= 0.45
    if p:
        expected = 1 if mutate == 'wrongexpect' else int(p.group(3))
        pok = p.group(8) == 'PASS' and abs(float(p.group(6)) - expected) < 0.05 and float(p.group(7)) <= 0.45
        pl = f"{p.group(2)} | {float(p.group(6)):+.4f} of {expected} | {'PASS' if pok else 'FAIL'}"
        ok &= pok
    else:
        pl, ok = " | | no PLANTED line", False
    bad |= not ok
    ow = f"{old['wind']:+.4f} | {old['step']} | {old['mind']}" if old else "missing | |"
    same = 'unchanged' if old and abs(new['wind'] - old['wind']) < 0.05 else 'CHANGED'
    rows.append(f"| {star} | {l} | {ow} | {new['wind']:+.4f} | {new['step']} | {new['mind']} | {pl} | {same} | {src} |")
    if old and old['edges'] and new['edges']:
        f3 = lambda e: ' | '.join(f"{x:+.4f}" for x in e)
        edges_out.append(f"| {star} | {l} | {f3(old['edges'])} | {f3(new['edges'])} |")

print("# Recount of #54, stage 1: Re s in [0.1, 2], abs(Im s) <= 2.5, old step rule against the current one")
print("| star | l | old winding | old max step | old min abs(D) | new winding | new max step | new min abs(D) | planted s0 | "
      "planted found | check | new against old | new log |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
print("\n".join(rows))
print()
print("# Phase gained along each edge of the upper half contour (rad): right (x1, 0 to y1), top (x1 to x0), left (y1 to 0)")
print("| star | l | old right | old top | old left | new right | new top | new left |")
print("|---|---|---|---|---|---|---|---|")
print("\n".join(edges_out))
print()
print("old: first step rule (I7-win2 logs for 0.70 and 0.85, winding2_output.txt for the others); new: current rule, "
      "planted zero 0.3 + 1.76i and its conjugate.")
print("RESULT: " + ("FAIL" if bad else "all rows complete and every planted zero returned")
      + (f"  [mutate={mutate}]" if mutate else ""))
sys.exit(1 if bad else 0)

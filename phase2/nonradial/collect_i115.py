"""
Collects the complex-s counts of #115 (phantom star omega = 0.7603, l = 4, 5, 6) from the raw logs into one table, so
that every number quoted in RESULTS.md comes from a file in this repository.  Usage, from this directory:

  python collect_i115.py ../../../negative-mass-raw/I115-rpi5 [--mutate=NAME] > i115_boxes_output.txt

Reads the RESULT and PLANTED lines of the winding2.py logs <box>_l<L>.log (box = rep, s1, n1) and of s1_l5_guardall.log
(s1 at l = 5 with --guard=all, the same box and plant), the control lines of control.log (control_winding2.py), and the
start and end times of queue.log.  Each log must hold the box it is named for (the tag of its RESULT line is compared
with the expected one) and the planted zero named for that box.  The guard of a run is read from its command in
jobs.txt (1 unless it has --guard=...) and must be 1 for the nine default runs and all for the guarded one.
Exit status 1 when
  - a box has no RESULT line or no PLANTED line, or its RESULT is for another box,
  - its planted zero is not the expected one, is not returned (PASS, expected count, difference within 0.05), or its max
    phase step exceeds 0.45 rad,
  - its winding is not within 0.05 of an integer, or its max phase step exceeds 0.45 rad,
  - its job exit status is missing or not 0,
  - its guard is not the expected one,
  - the control has no PASS on both lines or no "exit 0".

  --mutate=noplant      ignores the PLANTED lines (every box must then fail),
  --mutate=wrongexpect  takes 1 as the expected planted count of every box instead of the logged one,
  --mutate=wrongbox     expects the box of l + 1 in each log (every box must fail),
  --mutate=shift        adds 0.5 to every winding (the integer check must fail),
  --mutate=nocontrol    ignores control.log (the control check must fail),
  --mutate=wrongguard   expects guard all in every run (every default run must then fail).
"""
import re, sys, os
from datetime import datetime

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opts = dict((a[2:].split('=', 1) + ['1'])[:2] for a in sys.argv[1:] if a.startswith('--'))
mutate = opts.get('mutate')
if len(args) != 1 or set(opts) - {'mutate'} or mutate not in (None, 'noplant', 'wrongexpect', 'wrongbox', 'shift', 'nocontrol', 'wrongguard'):
    sys.exit(__doc__)
RAW = args[0]

res = re.compile(r"RESULT (\S+): zeros in Re s∈\[(\S+),(\S+)\], \|Im s\|≤(\S+): "
                 r"winding = ([+-]\d+\.\d+)  \(max phase step ([\d.]+) rad, min\|D\| (\S+)\)")
pla = re.compile(r"PLANTED (\S+): s0 = \((\S+)\), expected (\d) planted zeros inside; count with the factor ([+-]\d+\.\d+), "
                 r"without ([+-]\d+\.\d+), difference ([+-]\d+\.\d+)  \(max phase step ([\d.]+) rad\)  (PASS|FAIL)")
edge = re.compile(r"edge \S+->\S+: cumulative phase [+-]\d+\.\d+  \((\d+) evals")

# box -> (what, Re s lower, Re s upper, Im s bound, planted zero as printed)
BOX = {
    'rep': ('stored box', '0.1', '2.0', '2.5', '0.3+1.76j'),
    's1': ('strip', '0.03', '0.1', '2.5', '0.04+0.25j'),
    'n1': ('below the threshold 1 - w, right of the axis', '0.0001', '0.1', '0.23', '0.0005+0.12j'),
}
LS = (4, 5, 6)
# (box, l, job, expected guard): the nine default runs, then s1 at l = 5 again with --guard=all (same cache and tag)
JOBS = [(k, l, f"{k}_l{l}", '1') for l in LS for k in BOX] + [('s1', 5, 's1_l5_guardall', 'all')]

def read(name):
    f = os.path.join(RAW, name)
    return open(f).read() if os.path.exists(f) else None

def rc_of(job):
    t = read(os.path.join('done', job + '.rc'))
    return t.strip() if t is not None else None

def guard_of(job):
    """the --guard value in the job's command in jobs.txt, 1 when it has none"""
    for line in (read('jobs.txt') or '').splitlines():
        if line.startswith(job + '|'):
            m = re.search(r"--guard=(\S+)", line)
            return m.group(1) if m else '1'
    return None

def times(job):
    """first start and last end of a job in queue.log, in seconds (a restarted job has several start lines)"""
    t = read('queue.log') or ''
    st = [datetime.fromisoformat(l.split()[0]) for l in t.splitlines() if re.match(r"\S+ start " + re.escape(job) + "$", l)]
    en = [datetime.fromisoformat(l.split()[0]) for l in t.splitlines() if re.match(r"\S+ end " + re.escape(job) + r" rc=", l)]
    return (st[0], en[-1], len(st)) if st and en else None

bad = False
failed = []                     # (job, reason) of every box that fails
rows, runs = [], []
for key, l, job, eguard in JOBS:
    what, x0, x1, y1, s0 = BOX[key]
    lw = l + 1 if mutate == 'wrongbox' else l
    tag = f"eps-1_w0.7603_l{lw}_{x0}_{x1}_{y1}"
    guard = guard_of(job)
    txt = read(job + '.log')
    if txt is None:
        rows.append(f"| {key} | {l} | {guard} | {what} | | | not run | | | | | |"); bad = True; failed.append((job, 'not run')); continue
    r = next((m for m in res.finditer(txt) if m.group(1) == tag), None)
    if not r:
        rows.append(f"| {key} | {l} | {guard} | {what} | [{x0}, {x1}] | abs(Im s) ≤ {y1} | no result for {tag} | | | | | |")
        bad = True; failed.append((job, 'no RESULT line')); continue
    wind = float(r.group(5)) + (0.5 if mutate == 'shift' else 0.0)
    winds = f"{wind:+.4f}" if mutate == 'shift' else r.group(5)          # as logged, with the sign of a rounded 0
    step, mind = r.group(6), r.group(7)
    why = []
    want = 'all' if mutate == 'wrongguard' else eguard
    if guard != want: why.append(f"guard {guard}, expected {want}")
    if abs(wind - round(wind)) >= 0.05: why.append('winding not an integer')
    if float(step) > 0.45: why.append('max step above 0.45 rad')
    p = None if mutate == 'noplant' else next((m for m in pla.finditer(txt) if m.group(1) == tag), None)
    if p:
        expected = 1 if mutate == 'wrongexpect' else int(p.group(3))
        # the symmetric box doubles: a pair is planted, 2 zeros
        pok = (p.group(8) == 'PASS' and p.group(2) == s0 and expected == 2
               and abs(float(p.group(6)) - expected) < 0.05 and float(p.group(7)) <= 0.45)
        pl = f"{p.group(2)} | {float(p.group(6)):+.4f} of {expected} | {'PASS' if pok else 'FAIL'}"
        if not pok: why.append('planted zero not returned')
    else:
        pl = " | | no PLANTED line"; why.append('no PLANTED line')
    rc = rc_of(job)
    if rc != '0': why.append(f"job exit {rc}")
    if why: bad = True; failed.append((job, ', '.join(why)))
    ev = [int(m.group(1)) for m in edge.finditer(txt)]
    rows.append(f"| {key} | {l} | {guard} | {what} | [{x0}, {x1}] | abs(Im s) ≤ {y1} | {winds} | {step} | {mind} | {pl} |")
    t = times(job)
    wall = f"{(t[1] - t[0]).total_seconds() / 3600:.2f}" if t else "no times"
    runs.append(f"| {job} | {rc if rc is not None else 'missing'} | {ev[-1] if ev else ''} | {wall} | {t[2] if t else ''} |")

print("# winding2.py boxes (phantom, omega = 0.7603, l = 4, 5, 6); winding = zeros minus poles of D, a symmetric box doubled")
print("| box | l | guard | what | Re s | Im s | winding | max step | min abs(D) | planted s0 | planted found | check |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|")
print("\n".join(rows))
print()
print("# Runs: exit status, evaluations of D in the cache at the last edge, wall hours from the first start to the last end in")
print("# queue.log (includes any wait for load), starts (more than 1: restarted from the cache)")
print("| job | exit | evals | wall h | starts |")
print("|---|---|---|---|---|")
print("\n".join(runs))
print()
print("# control_winding2.py (Mini BS I, l = 2), same machine")
ctl = None if mutate == 'nocontrol' else read('control.log')
if ctl is None:
    print("  control: not run"); bad = True; failed.append(('control', 'not run'))
else:
    lines = [x for x in ctl.splitlines() if x.startswith('winding2 normalisation control') or x.startswith('exit')]
    for x in lines:
        print("  " + x.replace('winding2 normalisation control, ', ''))
    cl = [x for x in lines if x.startswith('winding2')]
    cok = len(cl) == 2 and all(x.rstrip().endswith('PASS') for x in cl) and 'exit 0' in lines
    bad |= not cok
    if not cok: print("  control: FAIL or incomplete"); failed.append(('control', 'FAIL or incomplete'))
print()
print("Each winding counts zeros minus poles of D (#124 counts the poles). No box is covered by a pole exclusion yet.")
for job, why in failed:
    print(f"FAILED {job}: {why}")
print("RESULT: " + ("FAIL" if bad else "all runs complete, every planted pair returned, control PASS")
      + (f"  [mutate={mutate}]" if mutate else ""))
sys.exit(1 if bad else 0)

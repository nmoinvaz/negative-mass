"""
Collects the results of #5 from the raw logs (negative-mass-raw, folder I5-rpi5) into one table, so that every
number quoted in RESULTS.md comes from a file in this repository.  Usage:

  python collect_i5.py RAW_FOLDER > i5_boxes_output.txt

Reads the RESULT and PLANTED lines of the winding2.py logs, the control_winding2.py lines, the polish_zero.py
RESULT lines and the refine_gamma.py level lines named in JOBS.  A job whose log is missing is listed as "not run".
Exit status 1 when a planted-zero line or a control line says FAIL, or when a winding2 run has no result or no PLANTED line.
"""
import re, sys, os

RAW = sys.argv[1]
W2 = [  # (job, what it counts)
    ('rep_l2', 'stored box, l = 2 (#7 box of the first scans)'),
    ('s1_l2', 'strip, l = 2'), ('s1_l3', 'strip, l = 3'),
    ('n1_l2', 'below threshold, right of the axis, l = 2'), ('n1_l3', 'below threshold, right of the axis, l = 3'),
    ('n3_l2', 'below threshold, across the axis, l = 2'), ('n3_l3', 'below threshold, across the axis, l = 3'),
    ('h1_l2', 'high Im s, l = 2'), ('h1_l3', 'high Im s, l = 3'),
    ('t1_r40', 'threshold 1 - w, outer radius 40, l = 2'), ('t1_r80', 'threshold 1 - w, outer radius 80, l = 2'),
    ('t2_r40', 'threshold 1 + w, outer radius 40, l = 2'), ('t2_r80', 'threshold 1 + w, outer radius 80, l = 2'),
    ('s2_l2', 'strip, closer to the axis, l = 2'), ('h2_l2', 'high Im s, second range, l = 2'),
]
res = re.compile(r"RESULT (\S+): zeros in Re s∈\[(\S+),(\S+)\], (?:\|Im s\|≤(\S+)|Im s∈\[(\S+),(\S+)\]): "
                 r"winding = ([+-]\d+\.\d+)  \(max phase step ([\d.]+) rad, min\|D\| (\S+)\)")
pla = re.compile(r"PLANTED (\S+): s0 = \((\S+)\), expected (\d) planted zeros inside; count with the factor ([+-]\d+\.\d+), "
                 r"without ([+-]\d+\.\d+), difference ([+-]\d+\.\d+)  \(max phase step ([\d.]+) rad\)  (PASS|FAIL)")
bad = False
print("# winding2.py boxes (phantom, omega = 0.7603)")
print("| job | what | Re s | Im s | winding | max step | min abs(D) | planted s0 | planted found | check |")
print("|---|---|---|---|---|---|---|---|---|---|")
for job, what in W2:
    f = os.path.join(RAW, job + '.log')
    if job == 'rep_l2': f = os.path.join(RAW, 'rep_l2.log')
    if not os.path.exists(f):
        print(f"| {job} | {what} | | | not run | | | | | |"); bad = True; continue
    txt = open(f).read()
    if job == 'rep_l2': txt = txt.split('## mutate=reverse')[0]       # the unmutated run
    r, p = res.search(txt), pla.search(txt)
    if not r:
        print(f"| {job} | {what} | | | no result yet | | | | | |"); bad = True; continue
    x0, x1, y1, ya, yb, wn, ms, mn = r.group(2), r.group(3), r.group(4), r.group(5), r.group(6), r.group(7), r.group(8), r.group(9)
    im = f"abs(Im s) ≤ {y1}" if y1 else f"[{ya}, {yb}]"
    if p:
        ok = p.group(8) == 'PASS'; bad |= not ok
        pl = f"{p.group(2)} | {p.group(6)} of {p.group(3)} | {p.group(8)}"
    else:
        bad = True; pl = " | | no PLANTED line"
    print(f"| {job} | {what} | [{x0}, {x1}] | {im} | {wn} | {ms} | {mn} | {pl} |")
print()
print("# Mutations of the planted-zero check on the stored box (each must FAIL)")
f = os.path.join(RAW, 'rep_l2.log')
if os.path.exists(f):
    for blk in open(f).read().split('## mutate=')[1:]:
        name = blk.split('\n')[0]; m = pla.search(blk)
        if m:
            print(f"  mutate={name}: difference {m.group(6)} against {m.group(3)} expected, {m.group(8)}")
            bad |= (name == 'none') != (m.group(8) == 'PASS')
print()
print("# control_winding2.py (Mini BS I, l = 2)")
for job in ['cs_default', 'cs_strip2', 'cs_strip2_r80', 'cs_strip2_shift']:
    f = os.path.join(RAW, job + '.log')
    if not os.path.exists(f): print(f"  {job}: not run"); continue
    for line in open(f).read().splitlines():
        if line.startswith('winding2 normalisation control') or line.startswith('exit'):
            print(f"  {job}: {line.replace('winding2 normalisation control, ', '')}")
            if 'mutate=' not in line and line.rstrip().endswith('FAIL'): bad = True
print()
print("# polish_zero.py --mini (Mini BS I zeros near the axis)")
f = os.path.join(RAW, 'polish_mini.log')
if os.path.exists(f):
    for line in open(f).read().splitlines():
        if line.startswith(('RESULT', 'CONTROL')): print('  ' + line)
print()
print("# refine_gamma.py (real-gamma refinement): the main run, then the planted pair and the three mutations")
reftab = []
for job in ['ref_l2', 'ref_l3']:
    f = os.path.join(RAW, job + '.log')
    if not os.path.exists(f): print(f"  {job}: not run"); continue
    parts = re.split(r"^== (.*)$", open(f).read(), flags=re.M)
    for seg, head in zip(parts[0::2], parts[1::2]):
        if head.startswith('main'):
            for m in re.finditer(r"^level (\d): step ratio (\S+), (\d+) points, sign changes (\d+)(?:, min \|D_i\|/sqrt\|D_\(i-1\) D_\(i\+1\)\| ([\d.]+), max ([\d.]+))?", seg, flags=re.M):
                sm = f"{m.group(5)} to {m.group(6)}" if m.group(5) else "not defined"
                reftab.append(f"| {job[-1]} | {m.group(1)} | {float(m.group(2)):.3f} | {m.group(3)} | {m.group(4)} | {sm} |")
        for line in seg.splitlines():
            if re.match(r"(reproduce|HALVING|PLANTED|level \d: step ratio)", line):
                print(f"  {job}: {line}")
        m = re.search(r"rc=(\d+) \((\d) expected\)", head)
        ok = bool(m) and m.group(1) == m.group(2); bad |= not ok
        print(f"  {job}: {head}  {'as expected' if ok else 'NOT AS EXPECTED'}")
    if len(parts) < 11: print(f"  {job}: only {len(parts) // 2} of 5 runs finished"); bad = True
print()
print("# refine_gamma.py, main runs, as a table")
print("| l | level | step ratio | points | sign changes | ratio of abs(D_i) to the geometric mean of its neighbours |")
print("|---|---|---|---|---|---|")
print("\n".join(reftab))
sys.exit(1 if bad else 0)

"""
Collects #126 from one negative-mass-raw folder: every box whose count RESULTS.md reports at the default step rule
(guard 1), recounted with --guard=all on a copy of its stored cache, next to its guard-1 value. Usage:

  python collect_i126.py --jobs                     the job list, one "job|script|arguments|source folder" per line
  python collect_i126.py RAW_FOLDER [--mutate=NAME] > i126_guard_output.txt

A job's log is RAW_FOLDER/<job>.log, the output of `<script> <arguments> --cachedir=RAW_FOLDER/cache`, where the cache
holds copies of the source folder's caches. RAW_FOLDER/jobs.txt must equal the --jobs output. Jobs `w_*` count the zeros
minus poles of D with winding2.py, jobs `z_*` the zeros of Z(r_m; s) with zerilli_zeros.py. The guard-1 values are read
from i54_recount_output.txt (section 5, "new winding"), i5_boxes_output.txt (5.1), i115_boxes_output.txt (3a, guard 1),
i121_zeros_output.txt and i124_zeros_output.txt (5.3).
Exit status 1 when jobs.txt differs from the job list, a log is missing or has no RESULT line, the box in a RESULT line
is not the job's, a count is not an integer, a planted check fails (winding2: the PLANTED line; zerilli_zeros: count,
holo, column, plant and "all checks passed"), a guard-all count differs from its guard-1 value (CHANGED), or a control
does not end as it must.
  --mutate=noplant   ignores the planted checks (every box must fail)
  --mutate=shift     adds 0.5 to the count of w_s1_l2 (not an integer)
  --mutate=changed   expects a guard-1 value 2 higher for w_s1_l2 (CHANGED)
  --mutate=noguard   expects the job list without --guard=all (jobs.txt differs)
  --mutate=control   expects z_mut_lshift to pass (the control check must fail)
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
W, Z, G = 'winding2.py', 'zerilli_zeros.py', ' --guard=all'
ST = '0.1 2.0 2.5 --plant=0.3,1.76'
# section 5 (#54's recount, I54-rpi5 caches; the phantom 0.7603, l = 2 box is #5's rep_l2, I5-rpi5)
S5 = [('p7603_l2', '-1 0.7603 2', 'phantom 0.7603', 2, 'I5-rpi5'), ('p7603_l3', '-1 0.7603 3', 'phantom 0.7603', 3, 'I54-rpi5'),
      ('c075_l2', '1 0.75 2', 'canonical 0.75', 2, 'I54-rpi5'), ('p070_l2', '-1 0.70 2', 'phantom 0.70', 2, 'I54-rpi5'),
      ('p070_l3', '-1 0.70 3', 'phantom 0.70', 3, 'I54-rpi5'), ('p085_l2', '-1 0.85 2', 'phantom 0.85', 2, 'I54-rpi5'),
      ('p085_l3', '-1 0.85 3', 'phantom 0.85', 3, 'I54-rpi5')]
# section 5.1 (#5, I5-rpi5 jobs.txt); rep_l2 is p7603_l2
S51 = [('n1_l2', '-1 0.7603 2 0.0001 0.1 0.23 --plant=0.0005,0.12'), ('n1_l3', '-1 0.7603 3 0.0001 0.1 0.23 --plant=0.0005,0.12'),
       ('n3_l2', '-1 0.7603 2 -0.001 0.0001 0.23 --y0=0.01 --plant=-0.0005,0.12'),
       ('n3_l3', '-1 0.7603 3 -0.001 0.0001 0.23 --y0=0.01 --plant=-0.0005,0.12'),
       ('s1_l2', '-1 0.7603 2 0.03 0.1 2.5 --plant=0.04,0.25'), ('s1_l3', '-1 0.7603 3 0.03 0.1 2.5 --plant=0.04,0.25'),
       ('h1_l2', '-1 0.7603 2 0.1 2.0 6.0 --y0=2.5 --plant=0.3,2.7'), ('h1_l3', '-1 0.7603 3 0.1 2.0 6.0 --y0=2.5 --plant=0.3,2.7'),
       ('t1_r40', '-1 0.7603 2 0.02 0.1 0.33 --y0=0.15 --plant=0.03,0.25'),
       ('t1_r80', '-1 0.7603 2 0.02 0.1 0.33 --y0=0.15 --rsmult=80 --zfac=2 --plant=0.03,0.25'),
       ('t2_r40', '-1 0.7603 2 0.02 0.1 1.84 --y0=1.68 --plant=0.03,1.77'),
       ('t2_r80', '-1 0.7603 2 0.02 0.1 1.84 --y0=1.68 --rsmult=80 --zfac=2 --plant=0.03,1.77'),
       ('s2_l2', '-1 0.7603 2 0.01 0.03 2.5 --plant=0.015,0.25'), ('h2_l2', '-1 0.7603 2 0.1 2.0 10.0 --y0=6.0 --plant=0.3,6.2'),
       ('s3_l2', '-1 0.7603 2 0.03 0.1 6.0 --y0=2.5 --plant=0.05,3.0'), ('s3_l3', '-1 0.7603 3 0.03 0.1 6.0 --y0=2.5 --plant=0.05,3.0')]
# section 3a (#115, I115-rpi5); s1 at l = 5 was counted with --guard=all there
B3A = [('rep', '0.1 2.0 2.5 --plant=0.3,1.76'), ('s1', '0.03 0.1 2.5 --plant=0.04,0.25'),
       ('n1', '0.0001 0.1 0.23 --plant=0.0005,0.12')]
# section 5.3: #121's 21 boxes (I121-win) and #124's nine (I124-rpi5), the arguments of collect_i121.py
S53 = [('p7603_l2', f'-1 0.7603 2 {ST}'), ('p7603_l3', f'-1 0.7603 3 {ST}'), ('c075_l2', f'1 0.75 2 {ST}'),
       ('p070_l2', f'-1 0.70 2 {ST}'), ('p070_l3', f'-1 0.70 3 {ST}'), ('p085_l2', f'-1 0.85 2 {ST}'),
       ('p085_l3', f'-1 0.85 3 {ST}')] + [(j, a) for j, a in S51 if not j.startswith('n3')]
S53 = [(j.replace('_r40', '_z1').replace('_r80', '_z2'), a.replace(' --rsmult=80', '')) for j, a in S53]

JOBS = [(f'w_{j}', W, f'{a} {ST}{G}', '5', ('i54', st, l), src) for j, a, st, l, src in S5]
JOBS += [(f'w_{j}', W, a + G, '5.1', ('i5', j), 'I5-rpi5') for j, a in S51]
JOBS += [(f'w_{b}_l{l}', W, f'-1 0.7603 {l} {a}{G}', '3a', ('i115', b, l), 'I115-rpi5')
         for l in (4, 5, 6) for b, a in B3A if (b, l) != ('s1', 5)]
JOBS += [(f'z_{j}', Z, a + G, '5.3', ('i121', j), 'I121-win') for j, a in S53]
JOBS += [(f'z_{b}_l{l}', Z, f'-1 0.7603 {l} {a}{G}', '5.3', ('i124', f'{b}_l{l}'), 'I124-rpi5')
         for l in (4, 5, 6) for b, a in B3A]
# (job, script, arguments, source, how it must end, the check that must fail or None)
CONTROLS = [('z_flat', Z, f'-1 0.7603 2 {ST.replace(" --", " --flat --")}{G}', 'I121-win', 'pass', None),
            ('z_mut_reverse', Z, f'-1 0.7603 2 0.1 2.0 2.5 --mutate=reverse --plant=0.3,1.76{G}', 'I121-win', 'fail', 'plant'),
            ('z_mut_nonorm', Z, f'-1 0.7603 2 0.1 2.0 2.5 --mutate=nonorm --plant=0.3,1.76{G}', 'I121-win', 'fail', 'holo'),
            ('z_mut_lshift', Z, f'-1 0.7603 2 0.1 2.0 2.5 --mutate=lshift --plant=0.3,1.76{G}', 'I121-win', 'fail', 'column'),
            ('w_mut_reverse', W, f'-1 0.7603 4 0.0001 0.1 0.23 --mutate=reverse --plant=0.0005,0.12{G}', 'I115-rpi5', 'fail',
             'PLANTED')]
MUT = next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--mutate=')), None)
assert MUT in (None, 'noplant', 'shift', 'changed', 'noguard', 'control'), MUT


def joblines(guard=True):
    return [f"{j}|{s}|{a if guard else a.replace(G, '')}|{src}" for j, s, a, _, _, src in JOBS] + \
           [f"{j}|{s}|{a if guard else a.replace(G, '')}|{src}" for j, s, a, src, _, _ in CONTROLS]


if '--jobs' in sys.argv:
    print('\n'.join(joblines()))
    sys.exit(0)
RAW = sys.argv[1]


def table(name, ncol, key, val):
    out = {}
    for ln in open(os.path.join(HERE, name), encoding='utf-8'):
        c = [x.strip() for x in ln.strip().strip('|').split('|')]
        if len(c) == ncol and re.match(r"[+-]\d+\.\d+$", c[val]):
            k = key(c)
            if k is not None:
                out[k] = float(c[val])
    return out


OLD = {}
OLD.update(table('i54_recount_output.txt', 13, lambda c: ('i54', c[0], int(c[1])) if c[1].isdigit() else None, 5))
OLD.update(table('i5_boxes_output.txt', 10, lambda c: ('i5', c[0]), 4))
OLD.update(table('i115_boxes_output.txt', 12, lambda c: ('i115', c[0], int(c[1])) if c[2] == '1' else None, 6))
OLD.update(table('i121_zeros_output.txt', 17, lambda c: ('i121', c[0]), 6))
OLD.update(table('i124_zeros_output.txt', 17, lambda c: ('i124', c[0]), 6))

resw = re.compile(r"RESULT (\S+): zeros in Re s∈\[(\S+),(\S+)\], (?:\|Im s\|≤(\S+)|Im s∈\[(\S+),(\S+)\]): winding = "
                  r"([+-]\d+\.\d+)  \(max phase step ([\d.]+) rad, min\|D\| (\S+)\)")
plw = re.compile(r"PLANTED .*expected (\d) planted zeros inside; .*difference ([+-]\d+\.\d+) .*\)  (PASS|FAIL)")
resz = re.compile(r"RESULT (\S+): zeros of Z\(r_m; s\) in Re s∈\[(\S+),(\S+)\], (?:\|Im s\|≤(\S+)|Im s∈\[(\S+),(\S+)\]): "
                  r"winding = ([+-]\d+\.\d+)  \(max phase step ([\d.]+) rad, min\|Z\| (\S+), (\d+) evaluations")
chk = re.compile(r"\[(PASS|FAIL)\] (\w+): (.*)")


def box(args):
    a = args.split()
    x0, x1, y1 = (float(v) for v in a[3:6])
    y0 = next((float(v.split('=')[1]) for v in a if v.startswith('--y0=')), 0.0)
    return x0, x1, y0, y1


def read(job):
    f = os.path.join(RAW, job + '.log')
    return open(f, encoding='utf-8', errors='replace').read() if os.path.exists(f) else None


bad = []
jf = os.path.join(RAW, 'jobs.txt')
got = open(jf, encoding='utf-8').read().split('\n') if os.path.exists(jf) else []
if [ln for ln in got if ln.strip()] != joblines(guard=MUT != 'noguard'):
    bad.append('jobs.txt is not the job list')
print("# Recount at --guard=all of the boxes counted at guard 1 (D: winding2.py, zeros minus poles; Z: zerilli_zeros.py)")
print(f"jobs.txt {'equals' if 'jobs.txt is not the job list' not in bad else 'DIFFERS FROM'} the job list "
      f"({len(JOBS)} boxes, {len(CONTROLS)} controls)")
print("| job | counts | section | Re s | Im s | guard 1 | guard all | against guard 1 | max step | min abs | evaluations | "
      "planted | check |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
nchanged = 0
for job, script, args, sec, src, _ in JOBS:
    x0, x1, y0, y1 = box(args)
    im = f"abs(Im s) ≤ {y1:g}" if y0 == 0 else f"[{y0:g}, {y1:g}]"
    old = OLD.get(src)
    if MUT == 'changed' and job == 'w_s1_l2' and old is not None:
        old += 2
    kind = 'D' if script == W else 'Z'
    txt = read(job)
    m = (resw if kind == 'D' else resz).search(txt or '')
    if m is None:
        print(f"| {job} | {kind} | {sec} | [{x0:g}, {x1:g}] | {im} | {old} | {'not run' if txt is None else 'no result'} "
              f"| | | | | | FAIL |")
        bad.append(f"{job}: no result")
        continue
    why = []
    rb = (float(m.group(2)), float(m.group(3)), float(m.group(5) or 0), float(m.group(4) or m.group(6)))
    if any(abs(u - v) > 1e-12 for u, v in zip(rb, (x0, x1, y0, y1))):
        why.append(f"box {rb} is not the job's")
    N = float(m.group(7)) + (0.5 if MUT == 'shift' and job == 'w_s1_l2' else 0.0)
    if abs(N - round(N)) >= 0.05:
        why.append("count not an integer")
    if kind == 'D':
        ev = (re.findall(r"\((\d+) evals", txt) or [''])[-1]
        p = plw.search(txt)
        pl = f"{p.group(2)} of {p.group(1)}" if p else ''
        if MUT == 'noplant' or p is None or p.group(3) != 'PASS':
            why.append("planted check not PASS")
    else:
        ev = m.group(10)
        checks = {c.group(2): c.group(1) for c in chk.finditer(txt)}
        last = [ln for ln in txt.splitlines() if re.match(r"(all checks passed|FAILED: )", ln)][-1:] or ['']
        pm = re.search(r"\[PASS\] plant: .*expected (\d), counted ([+-]\d+\.\d+)", txt)
        pl = f"{pm.group(2)} of {pm.group(1)}" if pm else ''
        for c in ('count', 'holo', 'column', 'plant'):
            if checks.get(c) != 'PASS' or (c == 'plant' and MUT == 'noplant'):
                why.append(f"{c} not PASS")
        if not last[0].startswith('all checks passed'):
            why.append("last line not 'all checks passed'")
    if old is None:
        why.append("no guard-1 value")
        vs = ''
    elif round(N) != round(old):
        why.append(f"CHANGED from {old:+.4f}")
        vs = 'CHANGED'
        nchanged += 1
    else:
        vs = 'unchanged'
    print(f"| {job} | {kind} | {sec} | [{x0:g}, {x1:g}] | {im} | {'' if old is None else f'{old:+.4f}'} | {N:+.4f} | {vs} | "
          f"{m.group(8)} | {m.group(9)} | {ev} | {pl} | {'PASS' if not why else 'FAIL'} |")
    bad += [f"{job}: {w}" for w in why]
print()
print(f"{nchanged} of {len(JOBS)} counts changed at --guard=all.")
print()
print("# Controls at --guard=all")
for job, script, args, src, end, must in CONTROLS:
    if MUT == 'control' and job == 'z_mut_lshift':
        end, must = 'pass', None
    txt = read(job)
    if txt is None:
        print(f"  {job}: not run  FAIL"); bad.append(f"{job}: not run"); continue
    if script == W:
        p = plw.search(txt)
        fl = ['PLANTED'] if p is None or p.group(3) != 'PASS' else []
        last = 'PLANTED ' + (p.group(3) if p else 'missing')
    else:
        fl = sorted(c.group(2) for c in chk.finditer(txt) if c.group(1) == 'FAIL')
        last = ([ln for ln in txt.splitlines() if re.match(r"(all checks passed|FAILED: )", ln)][-1:] or [''])[0].strip()
    ok = (not fl and (script == W or last.startswith('all checks passed'))) if end == 'pass' else (must in fl)
    print(f"  {job}: {last}; failed checks {fl or 'none'}, must {'pass' if end == 'pass' else 'fail ' + must}  "
          f"{'PASS' if ok else 'FAIL'}")
    if not ok:
        bad.append(f"{job}: did not {end}")
print()
print("RESULT: " + ("all rows complete, every planted check passed, no count changed, every control as required"
                    if not bad else "FAILED: " + "; ".join(bad)) + (f"  [mutate={MUT}]" if MUT else ""))
sys.exit(1 if bad else 0)

"""
Collects #117 from one negative-mass-raw folder (I117-rpi5): four boxes beyond Re s = 2 and |Im s| = 10 of the phantom
campaign star (eps = -1, omega = 0.7603), each counted with winding2.py (zeros minus poles of D) and with zerilli_zeros.py
(zeros of Z(r_m; s), the poles of D's Zerilli column), both at --guard=all with a planted zero. Usage:

  python collect_i117.py --jobs                     the job list, one "job|script|arguments" per line
  python collect_i117.py RAW_FOLDER [--mutate=NAME] > i117_boxes_output.txt

A job's log is RAW_FOLDER/<job>.log, the output of `<script> <arguments> --cachedir=RAW_FOLDER/cache`, and its exit
status is RAW_FOLDER/done/<job>.rc. RAW_FOLDER/jobs.txt must equal the --jobs output. The zeros of D in a box are
winding(D) + count(Z), since D has the poles of Z(r_m; s) (RESULTS.md 5.3).
Exit status 1 when jobs.txt differs from the job list, a log or exit-status file is missing, an exit status is not 0, a
log has no RESULT line, the box or planted zero in a log is not the job's, a count is not within 0.05 of an integer, a
planted check does not pass (winding2: the difference of the two counts against the planted zeros inside the box, within
0.05, max step below 0.45 rad, and the word PASS; zerilli_zeros: count, holo, column, plant and "all checks passed").
  --mutate=noplant   reads the count with the planted factor as the count without it (every planted check must fail)
  --mutate=shift     adds 0.5 to the winding of w_h3_l2 (not an integer)
  --mutate=nolog     reads no log for w_h2_l3 (a missing log)
  --mutate=badbox    expects w_b25_l2 on Re s up to 5.5 (the log's box is not the job's)
  --mutate=rc        reads exit status 1 for w_b25_l3
  --mutate=zcheck    reads the column check of z_b25_l2 as FAIL
  --mutate=noguard   expects the job list without --guard=all (jobs.txt differs)
"""
import os, re, sys
from datetime import datetime

W, Z, G = 'winding2.py', 'zerilli_zeros.py', ' --guard=all'
BOXES = [('b25_l2', 2, '2.0 5.0 10.0 --plant=2.5,1.0'), ('b25_l3', 3, '2.0 5.0 10.0 --plant=2.5,1.0'),
         ('h2_l3', 3, '0.1 2.0 10.0 --y0=6.0 --plant=0.3,6.2'), ('h3_l2', 2, '0.1 2.0 20.0 --y0=10.0 --plant=0.3,10.5')]
JOBS = [(f'w_{b}', W, f'-1 0.7603 {l} {a}{G}', b) for b, l, a in BOXES] + \
       [(f'z_{b}', Z, f'-1 0.7603 {l} {a}{G}', b) for b, l, a in BOXES]
MUT = next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--mutate=')), None)
assert MUT in (None, 'noplant', 'shift', 'nolog', 'badbox', 'rc', 'zcheck', 'noguard'), MUT


def joblines(guard=True):
    return [f"{j}|{s}|{a if guard else a.replace(G, '')}" for j, s, a, _ in JOBS]


if '--jobs' in sys.argv:
    print('\n'.join(joblines()))
    sys.exit(0)
RAW = next(a for a in sys.argv[1:] if not a.startswith('--'))

resw = re.compile(r"RESULT (\S+): zeros in Re s∈\[(\S+),(\S+)\], (?:\|Im s\|≤(\S+)|Im s∈\[(\S+),(\S+)\]): winding = "
                  r"([+-]\d+\.\d+)  \(max phase step ([\d.]+) rad, min\|D\| (\S+)\)")
plw = re.compile(r"PLANTED \S+: s0 = \((\S+)\), expected (\d) planted zeros inside; count with the factor ([+-]\d+\.\d+), "
                 r"without ([+-]\d+\.\d+), difference ([+-]\d+\.\d+)  \(max phase step ([\d.]+) rad\)  (PASS|FAIL)")
head = re.compile(r"\S+: r_m=([\d.]+), R_s fixed = ([\d.]+) \(min Re k on contour ([\d.]+)\)")
resz = re.compile(r"RESULT (\S+): zeros of Z\(r_m; s\) in Re s∈\[(\S+),(\S+)\], (?:\|Im s\|≤(\S+)|Im s∈\[(\S+),(\S+)\]): "
                  r"winding = ([+-]\d+\.\d+)  \(max phase step ([\d.]+) rad, min\|Z\| (\S+), (\d+) evaluations, \d+ s\)")
chk = re.compile(r"\[(PASS|FAIL)\] (\w+): (.*)")
plz = re.compile(r"\[(?:PASS|FAIL)\] plant: s0 = \((\S+)\), expected (\d), counted ([+-]\d+\.\d+) \(max step ([\d.]+) rad\)")
holo = re.compile(r"\[(?:PASS|FAIL)\] holo: Cauchy-Riemann defect (.*)")
colq = re.compile(r"q against zerilli_q: (.*)")


def box(args):
    a = args.split()
    x0, x1, y1 = (float(v) for v in a[3:6])
    y0 = next((float(v.split('=')[1]) for v in a if v.startswith('--y0=')), 0.0)
    s0 = complex(*(float(v) for v in next(v for v in a if v.startswith('--plant=')).split('=')[1].split(',')))
    return x0, x1, y0, y1, s0


def expected_planted(x0, x1, y0, y1, s0):
    sym = y0 == 0.0
    inside = x0 < s0.real < x1 and ((0 < abs(s0.imag) < y1) if sym else (y0 < s0.imag < y1))
    return (2 if sym else 1) * inside


def read(path, job):
    if MUT == 'nolog' and job == 'w_h2_l3':
        return None
    return open(path, encoding='utf-8', errors='replace').read() if os.path.exists(path) else None


def walls():
    """seconds per job from queue.log start and end lines (a restarted job sums its runs)"""
    out, start = {}, {}
    p = os.path.join(RAW, 'queue.log')
    for ln in (open(p, encoding='utf-8').read().splitlines() if os.path.exists(p) else []):
        m = re.match(r"(\S+) (?:15-min load [\d.]+ < 4: )?start (\S+?): ", ln)
        if m:
            start[m.group(2)] = datetime.fromisoformat(m.group(1))
        m = re.match(r"(\S+) end (\S+) rc=", ln)
        if m and m.group(2) in start:
            out[m.group(2)] = out.get(m.group(2), 0.0) + (datetime.fromisoformat(m.group(1)) - start.pop(m.group(2))).total_seconds()
    return out


def hm(sec):
    return '' if sec is None else f"{int(sec // 3600)} h {int(sec % 3600 // 60):02d} min"


bad = []
jf = os.path.join(RAW, 'jobs.txt')
got = [ln for ln in open(jf, encoding='utf-8').read().split('\n') if ln.strip()] if os.path.exists(jf) else []
if got != joblines(guard=MUT != 'noguard'):
    bad.append('jobs.txt is not the job list')
print("# #117: four boxes of the phantom campaign star (eps = -1, omega = 0.7603) at --guard=all")
print(f"jobs.txt {'DIFFERS FROM' if 'jobs.txt is not the job list' in bad else 'equals'} the job list ({len(JOBS)} jobs)")

info = {}  # job -> dict of what the log says
for job, script, args, bx in JOBS:
    x0, x1, y0, y1, s0 = box(args)
    if MUT == 'badbox' and job == 'w_b25_l2':
        x1 = 5.5
    kind = 'D' if script == W else 'Z'
    d = info[job] = dict(kind=kind, why=[])
    why = d['why']
    txt = read(os.path.join(RAW, job + '.log'), job)
    rcf = os.path.join(RAW, 'done', job + '.rc')
    rc = open(rcf).read().strip() if os.path.exists(rcf) else None
    if MUT == 'rc' and job == 'w_b25_l3':
        rc = '1'
    d['rc'] = rc
    if txt is None:
        why.append('no log')
        continue
    if rc is None:
        why.append('no exit status')
    elif rc != '0':
        why.append(f'exit status {rc}')
    m = (resw if kind == 'D' else resz).search(txt)
    if m is None:
        why.append('no RESULT line')
        continue
    rb = (float(m.group(2)), float(m.group(3)), float(m.group(5) or 0), float(m.group(4) or m.group(6)))
    if any(abs(u - v) > 1e-12 for u, v in zip(rb, (x0, x1, y0, y1))):
        why.append(f"box {rb} is not the job's")
    N = float(m.group(7)) + (0.5 if MUT == 'shift' and job == 'w_h3_l2' else 0.0)
    d.update(maxstep=m.group(8), minabs=m.group(9))
    exp = expected_planted(x0, x1, y0, y1, s0)
    if kind == 'D':
        hm_ = head.search(txt)
        d.update(rm=hm_.group(1) if hm_ else '', Rs=hm_.group(2) if hm_ else '')
        d['evals'] = (re.findall(r"\((\d+) evals", txt) or [''])[-1]
        p = plw.search(txt)
        if p is None:
            why.append('no PLANTED line')
        else:
            diff = float(p.group(5))
            if MUT == 'noplant':  # the count with the factor is read as the count without it
                N, diff = float(p.group(3)), 0.0
            if abs(complex(p.group(1)) - s0) > 1e-12:
                why.append("planted zero is not the job's")
            if int(p.group(2)) != exp:
                why.append(f"log expects {p.group(2)} planted zeros, the box holds {exp}")
            if abs(float(p.group(3)) - float(p.group(4)) - float(p.group(5))) > 2e-4:
                why.append('difference is not with minus without')
            if not (abs(diff - exp) < 0.05 and float(p.group(6)) < 0.45 and p.group(7) == 'PASS'):
                why.append(f"planted check: difference {diff:+.4f} against {exp}, max step {p.group(6)}, {p.group(7)}")
            d.update(s0=p.group(1), found=f"{diff:+.4f} of {exp}")
    else:
        d['evals'] = m.group(10)
        checks = {c.group(2): c.group(1) for c in chk.finditer(txt)}
        if MUT == 'zcheck' and job == 'z_b25_l2':
            checks['column'] = 'FAIL'
        last = ([ln for ln in txt.splitlines() if re.match(r"(all checks passed|FAILED: )", ln)] or [''])[-1]
        pm = plz.search(txt)
        if pm is None:
            why.append('no plant line')
        else:
            diff = float(pm.group(3))
            if MUT == 'noplant':
                N, diff = N + diff, 0.0
            if abs(complex(pm.group(1)) - s0) > 1e-12:
                why.append("planted zero is not the job's")
            if int(pm.group(2)) != exp:
                why.append(f"log expects {pm.group(2)} planted zeros, the box holds {exp}")
            if not (abs(diff - exp) < 0.05 and float(pm.group(4)) < 0.45):
                why.append(f"planted check: difference {diff:+.4f} against {exp}, max step {pm.group(4)}")
            d.update(s0=pm.group(1), found=f"{diff:+.4f} of {exp}")
        for c in ('count', 'holo', 'column', 'plant'):
            if checks.get(c) != 'PASS':
                why.append(f"{c} not PASS")
        if not last.startswith('all checks passed'):
            why.append("last line is not 'all checks passed'")
        h, q = holo.search(txt), colq.search(txt)
        d.update(holo=max((float(v) for v in re.findall(r"[\d.]+e[+-]\d+", h.group(1))), default=None) if h else None,
                 q=max((float(v) for v in re.findall(r": ([\d.]+e[+-]\d+)", q.group(1))), default=None) if q else None)
    if abs(N - round(N)) >= 0.05:
        why.append('count not an integer')
    d['N'] = N

wl = walls()
print("| box | l | Re s | Im s | winding of D | max step | min abs(D) | evaluations | R_s | planted s0 | planted found (D) | "
      "count of Z | planted found (Z) | zeros of D |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for b, l, a in BOXES:
    x0, x1, y0, y1, s0 = box(f'-1 0.7603 {l} {a}')
    im = f"abs(Im s) ≤ {y1:g}" if y0 == 0 else f"[{y0:g}, {y1:g}]"
    dd, zz = info['w_' + b], info['z_' + b]
    if 'N' not in dd or 'N' not in zz:
        print(f"| {b} | {l} | [{x0:g}, {x1:g}] | {im} | " + ('not run' if 'N' not in dd else f"{dd['N']:+.4f}") + " | | | | | | | " +
              ('not run' if 'N' not in zz else f"{zz['N']:+.4f}") + " | | |")
        continue
    zd = round(dd['N']) + round(zz['N'])
    print(f"| {b} | {l} | [{x0:g}, {x1:g}] | {im} | {dd['N']:+.4f} | {dd['maxstep']} | {dd['minabs']} | {dd['evals']} | "
          f"{dd['Rs']} | {dd.get('s0', '')} | {dd.get('found', '')} | {zz['N']:+.4f} | {zz.get('found', '')} | {zd:d} |")
print()
print("| job | script | exit status | max step | min abs | evaluations | wall | checks |")
print("|---|---|---|---|---|---|---|---|")
for job, script, args, bx in JOBS:
    d = info[job]
    extra = ''
    if d['kind'] == 'Z' and 'holo' in d:
        extra = f"holo defect ≤ {d['holo']:.1e}" if d['holo'] is not None else ''
        extra += f", q against zerilli_q ≤ {d['q']:.1e}" if d.get('q') is not None else ''
    print(f"| {job} | {script} | {d['rc'] if d['rc'] is not None else 'none'} | {d.get('maxstep', '')} | {d.get('minabs', '')} | "
          f"{d.get('evals', '')} | {hm(wl.get(job))} | {'PASS' if not d['why'] else 'FAIL'}{'; ' + extra if extra else ''} |")
    bad += [f"{job}: {w}" for w in d['why']]
print()
print("RESULT: " + ("all eight jobs exit 0, every count is an integer and every planted check passed"
                    if not bad else "FAILED: " + "; ".join(bad)) + (f"  [mutate={MUT}]" if MUT else ""))
sys.exit(1 if bad else 0)

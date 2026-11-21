"""
Collects #121 and #124 from the zerilli_zeros.py logs (negative-mass-raw, one folder) into one table: the zeros of
Z(r_m; s) in each box whose winding of D is recorded, and the zeros of D there, winding(D) + count(Z). Usage:

  python collect_i121.py --jobs [--only=3a]           the job list, one "job|arguments of zerilli_zeros.py" per line
  python collect_i121.py RAW_FOLDER [--mutate=NAME] > i121_zeros_output.txt                 (#121, sections 5 and 5.1)
  python collect_i121.py RAW_FOLDER --only=3a [--mutate=NAME] > i124_zeros_output.txt       (#124, section 3a)

--only=3a selects the nine boxes of section 3a (#115, l = 4, 5, 6) and the controls; without it, the 21 boxes of #121.

A job's log is RAW_FOLDER/<job>.log, the output of `zerilli_zeros.py <arguments> --cachedir=RAW_FOLDER/cache`.
The windings of D are read from i54_recount_output.txt (section 5, "new winding"), i5_boxes_output.txt (section 5.1)
and i115_boxes_output.txt (section 3a: s1 at l = 5 from its --guard=all row, the other eight at guard 1).
Exit status 1 when a log is missing or has no RESULT line, when a box job's count, holo, column or plant check is not
PASS or its last line is not "all checks passed", when the box in a RESULT line is not the job's, when a count is not an
integer, or when a control does not end as it must (flat passes; each mutation fails, on its own check).
  --mutate=box      expects s1_l2 (s1_l4 with --only=3a) at Re s in [0.03, 0.2] (the box check must fail)
  --mutate=control  expects mut_lshift to pass (the control check must fail)
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
STORED = '-1 0.7603 2 0.1 2.0 2.5'
# (job, section, star, arguments, where the winding of D is): section 5 (#54's recount), section 5.1 (#5), section 3a
# (#115)
JOBS = [
    ('p7603_l2', '5', 'phantom 0.7603', f'{STORED} --plant=0.3,1.76', ('i54', 'phantom 0.7603', 2)),
    ('p7603_l3', '5', 'phantom 0.7603', '-1 0.7603 3 0.1 2.0 2.5 --plant=0.3,1.76', ('i54', 'phantom 0.7603', 3)),
    ('c075_l2', '5', 'canonical 0.75', '1 0.75 2 0.1 2.0 2.5 --plant=0.3,1.76', ('i54', 'canonical 0.75', 2)),
    ('p070_l2', '5', 'phantom 0.70', '-1 0.70 2 0.1 2.0 2.5 --plant=0.3,1.76', ('i54', 'phantom 0.70', 2)),
    ('p070_l3', '5', 'phantom 0.70', '-1 0.70 3 0.1 2.0 2.5 --plant=0.3,1.76', ('i54', 'phantom 0.70', 3)),
    ('p085_l2', '5', 'phantom 0.85', '-1 0.85 2 0.1 2.0 2.5 --plant=0.3,1.76', ('i54', 'phantom 0.85', 2)),
    ('p085_l3', '5', 'phantom 0.85', '-1 0.85 3 0.1 2.0 2.5 --plant=0.3,1.76', ('i54', 'phantom 0.85', 3)),
    ('s1_l2', '5.1', 'phantom 0.7603', '-1 0.7603 2 0.03 0.1 2.5 --plant=0.04,0.25', ('i5', 's1_l2')),
    ('s1_l3', '5.1', 'phantom 0.7603', '-1 0.7603 3 0.03 0.1 2.5 --plant=0.04,0.25', ('i5', 's1_l3')),
    ('n1_l2', '5.1', 'phantom 0.7603', '-1 0.7603 2 0.0001 0.1 0.23 --plant=0.0005,0.12', ('i5', 'n1_l2')),
    ('n1_l3', '5.1', 'phantom 0.7603', '-1 0.7603 3 0.0001 0.1 0.23 --plant=0.0005,0.12', ('i5', 'n1_l3')),
    ('h1_l2', '5.1', 'phantom 0.7603', '-1 0.7603 2 0.1 2.0 6.0 --y0=2.5 --plant=0.3,2.7', ('i5', 'h1_l2')),
    ('h1_l3', '5.1', 'phantom 0.7603', '-1 0.7603 3 0.1 2.0 6.0 --y0=2.5 --plant=0.3,2.7', ('i5', 'h1_l3')),
    ('t1_z1', '5.1', 'phantom 0.7603', '-1 0.7603 2 0.02 0.1 0.33 --y0=0.15 --plant=0.03,0.25', ('i5', 't1_r40')),
    ('t1_z2', '5.1', 'phantom 0.7603', '-1 0.7603 2 0.02 0.1 0.33 --y0=0.15 --zfac=2 --plant=0.03,0.25',
     ('i5', 't1_r80')),
    ('t2_z1', '5.1', 'phantom 0.7603', '-1 0.7603 2 0.02 0.1 1.84 --y0=1.68 --plant=0.03,1.77', ('i5', 't2_r40')),
    ('t2_z2', '5.1', 'phantom 0.7603', '-1 0.7603 2 0.02 0.1 1.84 --y0=1.68 --zfac=2 --plant=0.03,1.77',
     ('i5', 't2_r80')),
    ('s2_l2', '5.1', 'phantom 0.7603', '-1 0.7603 2 0.01 0.03 2.5 --plant=0.015,0.25', ('i5', 's2_l2')),
    ('h2_l2', '5.1', 'phantom 0.7603', '-1 0.7603 2 0.1 2.0 10.0 --y0=6.0 --plant=0.3,6.2', ('i5', 'h2_l2')),
    ('s3_l2', '5.1', 'phantom 0.7603', '-1 0.7603 2 0.03 0.1 6.0 --y0=2.5 --plant=0.05,3.0', ('i5', 's3_l2')),
    ('s3_l3', '5.1', 'phantom 0.7603', '-1 0.7603 3 0.03 0.1 6.0 --y0=2.5 --plant=0.05,3.0', ('i5', 's3_l3')),
] + [(f'{b}_l{l}', '3a', 'phantom 0.7603', f'-1 0.7603 {l} {box} --plant={pl}',
      ('i115', b, l, 'all' if (b, l) == ('s1', 5) else '1'))
     for l in (4, 5, 6)
     for b, box, pl in (('rep', '0.1 2.0 2.5', '0.3,1.76'), ('s1', '0.03 0.1 2.5', '0.04,0.25'),
                        ('n1', '0.0001 0.1 0.23', '0.0005,0.12'))]
ONLY = next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--only=')), None)
assert ONLY in (None, '3a'), ONLY
JOBS = [j for j in JOBS if (j[1] == '3a') == (ONLY == '3a')]
# boxes of D not counted here, with the reason: k = sqrt(s^2 + V) has its cut on Re s = 0, inside these boxes
NOT_COVERED = [('rep_l2', 'the same box as p7603_l2'),
               ('n3_l2', 'crosses Re s = 0, where the cut of k makes y1(r_m) not holomorphic'),
               ('n3_l3', 'crosses Re s = 0, where the cut of k makes y1(r_m) not holomorphic')] if ONLY is None else []
# (job, arguments, how it must end, the check that must fail or None)
CONTROLS = [('flat', f'{STORED} --flat --plant=0.3,1.76', 'pass', None),
            ('mut_reverse', f'{STORED} --mutate=reverse --plant=0.3,1.76', 'fail', 'plant'),
            ('mut_nonorm', f'{STORED} --mutate=nonorm --plant=0.3,1.76', 'fail', 'holo'),
            ('mut_lshift', f'{STORED} --mutate=lshift --plant=0.3,1.76', 'fail', 'column')]

if '--jobs' in sys.argv:
    for job, *rest in JOBS:
        print(f"{job}|{rest[2]}")
    for job, a, _, _ in CONTROLS:
        print(f"{job}|{a}")
    sys.exit(0)
RAW = sys.argv[1]
MUT = next((a.split('=', 1)[1] for a in sys.argv[2:] if a.startswith('--mutate=')), None)
assert MUT in (None, 'box', 'control'), MUT

res = re.compile(r"RESULT (\S+): zeros of Z\(r_m; s\) in Re s∈\[(\S+),(\S+)\], (?:\|Im s\|≤(\S+)|Im s∈\[(\S+),(\S+)\]): "
                 r"winding = ([+-]\d+\.\d+)  \(max phase step ([\d.]+) rad, min\|Z\| (\S+), (\d+) evaluations, (\d+) s\)")
chk = re.compile(r"\[(PASS|FAIL)\] (\w+): (.*)")


def winding_D():
    out = {}
    for ln in open(os.path.join(HERE, 'i54_recount_output.txt')):
        c = [x.strip() for x in ln.strip().strip('|').split('|')]
        if len(c) == 13 and c[1].isdigit():
            out[('i54', c[0], int(c[1]))] = c[5]
    for ln in open(os.path.join(HERE, 'i5_boxes_output.txt')):
        c = [x.strip() for x in ln.strip().strip('|').split('|')]
        if len(c) == 10 and re.match(r"[+-]\d+\.\d+$", c[4]):
            out[('i5', c[0])] = c[4]
    for ln in open(os.path.join(HERE, 'i115_boxes_output.txt')):
        c = [x.strip() for x in ln.strip().strip('|').split('|')]
        if len(c) == 12 and c[1].isdigit() and c[2] in ('1', 'all'):
            out[('i115', c[0], int(c[1]), c[2])] = c[6]
    return out


def box(args):
    a = args.split()
    x0, x1, y1 = (float(v) for v in a[3:6])
    y0 = next((float(v.split('=')[1]) for v in a if v.startswith('--y0=')), 0.0)
    if MUT == 'box' and a[3:6] == ['0.03', '0.1', '2.5'] and a[2] == ('4' if ONLY else '2'):
        x1 = 0.2
    return x0, x1, y0, y1


def read(job):
    f = os.path.join(RAW, job + '.log')
    if not os.path.exists(f):
        return None
    txt = open(f, encoding='utf-8', errors='replace').read()
    checks = {m.group(2): (m.group(1), m.group(3)) for m in chk.finditer(txt)}
    last = [ln for ln in txt.splitlines() if re.match(r"(all checks passed|FAILED: )", ln)][-1:] or ['']
    return res.search(txt), checks, last[0], txt


bad = []
WD = winding_D()
print("# Zeros of Z(r_m; s), the poles of D, in the boxes of RESULTS.md "
      + ("section 3a (#115)" if ONLY else "section 5 (#54) and section 5.1 (#5)"))
print("| job | section | star | l | Re s | Im s | zeros of Z | max step | min abs(Z) | evaluations | time (s) | "
      "Cauchy-Riemann, worst | q against zerilli_q, worst (Re s ≥ 0.03) | planted | winding of D | zeros of D | check |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for job, sec, star, args, src in JOBS:
    l = args.split()[2]
    x0, x1, y0, y1 = box(args)
    im = f"abs(Im s) ≤ {y1:g}" if y0 == 0 else f"[{y0:g}, {y1:g}]"
    wd = WD.get(src)
    r = read(job)
    if r is None or r[0] is None:
        print(f"| {job} | {sec} | {star} | {l} | [{x0:g}, {x1:g}] | {im} | {'not run' if r is None else 'no result'} "
              f"| | | | | | | | {wd} | | FAIL |")
        bad.append(f"{job}: no result")
        continue
    m, checks, last, txt = r
    why = []
    rb = (float(m.group(2)), float(m.group(3)), float(m.group(5) or 0), float(m.group(4) or m.group(6)))
    if any(abs(u - v) > 1e-12 for u, v in zip(rb, (x0, x1, y0, y1))):
        why.append(f"box {rb} is not the job's")
    for c in ('count', 'holo', 'column', 'plant'):
        if checks.get(c, ('missing',))[0] != 'PASS':
            why.append(f"{c} {checks.get(c, ('missing',))[0]}")
    if not last.startswith('all checks passed'):
        why.append("last line not 'all checks passed'")
    N = float(m.group(7))
    if abs(N - round(N)) >= 0.05:
        why.append("count not an integer")
    if wd is None:
        why.append("no winding of D")
    cr = checks.get('holo', ('', ''))[1].replace('Cauchy-Riemann defect ', '')
    crw = max((float(v) for v in cr.split(', ') if v), default=float('nan'))
    cols = [(float(a), float(d)) for a, d in re.findall(r"s = ([\d.e-]+)[+-][\d.e-]+i: ([\d.e+-]+)", txt)]
    colw = max((d for a, d in cols if a >= 0.03), default=None)
    pm = re.search(r"expected (\d), counted ([+-]\d+\.\d+)", checks.get('plant', ('', ''))[1])
    pl = f"{pm.group(2)} of {pm.group(1)}" if pm else ""
    zd = f"{round(float(wd)) + round(N):+d}" if wd is not None else ""
    print(f"| {job} | {sec} | {star} | {l} | [{x0:g}, {x1:g}] | {im} | {N:+.4f} | {m.group(8)} | {m.group(9)} | "
          f"{m.group(10)} | {m.group(11)} | {crw:.1e} | {'' if colw is None else f'{colw:.1e}'} | {pl} | {wd} | {zd} | "
          f"{'PASS' if not why else 'FAIL'} |")
    bad += [f"{job}: {w}" for w in why]
print()
print("Not counted here: " + ("; ".join(f"{j} ({w})" for j, w in NOT_COVERED) or "none") + ".")
print()
print("# Controls on the stored box (phantom 0.7603, l = 2)")
for job, args, end, must in CONTROLS:
    if MUT == 'control' and job == 'mut_lshift':
        end, must = 'pass', None
    r = read(job)
    if r is None:
        print(f"  {job}: not run  FAIL"); bad.append(f"{job}: not run"); continue
    m, checks, last, _ = r
    fl = sorted(c for c, (v, _) in checks.items() if v == 'FAIL')
    if end == 'pass':
        ok = last.startswith('all checks passed') and not fl and all(
            checks.get(c, ('',))[0] == 'PASS' for c in (('flat', 'flatq', 'plant') if job == 'flat' else ()))
    else:
        ok = last.startswith('FAILED') and must in fl
    print(f"  {job} ({args.split(' ', 6)[-1]}): {last.strip()}; failed checks {fl or 'none'}, must "
          f"{'pass' if end == 'pass' else 'fail ' + must}  {'PASS' if ok else 'FAIL'}")
    if not ok:
        bad.append(f"{job}: did not {end}")
print()
print("RESULT: " + ("all rows complete, every planted zero returned, every control as required" if not bad else
                    "FAILED: " + "; ".join(bad)) + (f"  [mutate={MUT}]" if MUT else ""))
sys.exit(1 if bad else 0)

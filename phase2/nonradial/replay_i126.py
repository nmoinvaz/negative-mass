"""
Replay of the step rule against the dense chain of every box #126 recounted at --guard=all (#128). Usage:

  python replay_i126.py RAW_FOLDER [--mutate=NAME] > i126_replay_output.txt

RAW_FOLDER is negative-mass-raw's I126-rpi5: its jobs.txt (the 66 lines of `python collect_i126.py --jobs`) and the
pickle caches in RAW_FOLDER/cache, read and never written. For each of the 61 boxes (the 31 `w_*` and 30 `z_*` jobs, not
the five controls, which carry --flat or --mutate) it runs

  python diagnose_plant.py CACHE x0 x1 y1 RE,IM [--y0=Y0] --guard=all --expect=none

with the box, the `--y0` and the planted zero of jobs.txt, and prints one row per box: the counted function F (D or Z),
the accepted steps of F and of F times the planted factor, the chain's points, the replayed counts of F and of F times the
factor, their difference (the planted zeros found), and the chain's largest increment for each. Status words: numerical
evidence. A box PASSes when
  replay    diagnose_plant.py exits 0: the cache holds every point the counter asks for, no accepted step differs from the
            dense chain in F or in F times the factor by more than 0.5 rad, and the chain's largest increment is below 1 rad;
  chain     the chain holds every cached point on the contour: its points (derived from the printed edges) equal the
            points of the cache that lie on the box's contour, counted here from the pickle by comparing coordinates;
  planted   the count of F times the factor exceeds the count of F by the planted zeros inside the box (2 for a pair in a
            symmetric box, 1 in an upper box, 0 outside) within 0.05, as winding2.py's PLANTED line does.
The script also needs jobs.txt to equal the job list of collect_i126.py --jobs. It shares with diagnose_plant.py and with
the collector only the cache files, jobs.txt and the geometry of the contour. Exit status 1 when any box fails, a cache is
missing, a replay stops on a point the cache lacks, or jobs.txt is not the job list.
  --mutate=sparse     passes --mutate=sparse to diagnose_plant.py on w_h1_l2: its chain keeps only the accepted steps' end
                      points. That chain agrees with the steps whatever they are, so `replay` alone passes; `chain` fails.
  --mutate=noy0       replays the 20 upper boxes without --y0 (their contour then starts at (x1, 0), which their caches lack);
  --mutate=misplant   computes the expected planted zeros of w_h1_l2 from a zero moved out of the box (10 box widths to the
                      right) while the replay keeps the real one: 1 are found, 0 expected;
  --mutate=nocache    looks for the cache of w_h1_l2 under a name that does not exist;
  --mutate=nojobs     drops the line of w_p7603_l2 from the jobs.txt that it reads.
"""
import os, pickle, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = sys.argv[1]
MUT = next((a.split('=', 1)[1] for a in sys.argv[2:] if a.startswith('--mutate=')), None)
assert MUT in (None, 'sparse', 'noy0', 'misplant', 'nocache', 'nojobs'), MUT
MUTJOB = 'w_h1_l2'
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')

bad = []
lines = [ln for ln in open(os.path.join(RAW, 'jobs.txt'), encoding='utf-8').read().split('\n') if ln.strip()]
if MUT == 'nojobs':
    lines = [ln for ln in lines if not ln.startswith('w_p7603_l2|')]
want = subprocess.run([sys.executable, os.path.join(HERE, 'collect_i126.py'), '--jobs'], capture_output=True, text=True,
                      check=True).stdout.split('\n')
jobs_ok = lines == [ln for ln in want if ln.strip()]
if not jobs_ok:
    bad.append('jobs.txt is not the job list')
boxes = []
for ln in lines:
    job, script, arguments, _ = ln.split('|')
    if '--mutate=' not in arguments and '--flat' not in arguments:
        boxes.append((job, script, arguments))
nD = sum(s == 'winding2.py' for _, s, _ in boxes)


def tag(script, arguments):
    """the cache file name of winding2.py and zerilli_zeros.py for a job's arguments"""
    pos = [a for a in arguments.split() if not a.startswith('--')]
    opt = dict(a[2:].split('=', 1) for a in arguments.split() if a.startswith('--') and '=' in a)
    eps, w, l = int(pos[0]), float(pos[1]), int(pos[2])
    x0, x1, y1 = (float(a) for a in pos[3:6])
    y0, rs, zf = float(opt.get('y0', 0.0)), float(opt.get('rsmult', 40.0)), float(opt.get('zfac', 1.0))
    t = f"eps{eps:+d}_w{w}_l{l}_{x0}_{x1}_{y1}" + ('' if y0 == 0.0 else f"_y0{y0}")
    if script == 'winding2.py':
        t += '' if rs == 40.0 else f"_rs{rs:g}"
    t += '' if zf == 1.0 else f"_z{zf:g}"
    return ('winding2_cache_' if script == 'winding2.py' else 'zerilli_zeros_cache_') + t + '.pkl', (x0, x1, y0, y1), opt


def contour_points(keys, x0, x1, y0, y1):
    """cached points on the box's contour: Re s = x0 or x1 with Im s in [y0, y1], Im s = y1 or y0 (upper box) with Re s in
    [x0, x1]; y0 = 0 is the symmetric box, whose contour has no lower edge"""
    tol, up = 1e-9, y0 > 0
    n = 0
    for s in keys:
        vert = (abs(s.real - x0) < tol or abs(s.real - x1) < tol) and y0 - tol <= s.imag <= y1 + tol
        horz = (abs(s.imag - y1) < tol or (up and abs(s.imag - y0) < tol)) and x0 - tol <= s.real <= x1 + tol
        n += vert or horz
    return n


EDGE = re.compile(r"^  edge .*: adaptive ([+-]\d+\.\d+) \((\d+) steps, max (\d+\.\d+)\), dense chain ([+-]\d+\.\d+) "
                  r"\((\d+) steps, max (\d+\.\d+)\), difference/2pi ([+-]\d+\.\d+)$", re.M)
TOT = re.compile(r"^  total adaptive ([+-]\d+\.\d+) -> count ([+-]\d+\.\d+);  dense chain ([+-]\d+\.\d+) -> count "
                 r"([+-]\d+\.\d+)$", re.M)
CHK = re.compile(r"^CHECK expected .*; chain's largest increment (\d+\.\d+) rad  (PASS|FAIL.*?)(?:  \[mutate=\w+\])?$", re.M)
MISS = re.compile(r"^CHECK FAIL: the cache lacks the point (\S+)", re.M)


def parse(out, sym):
    """per function (F, F times the factor): accepted steps, chain steps, count of the steps, largest chain increment"""
    secs = out.split('\n== ')[1:]
    if len(secs) != 2:
        return None
    res = []
    for sec in secs:
        e, t = EDGE.findall(sec), TOT.findall(sec)
        if len(e) != (3 if sym else 4) or len(t) != 1:
            return None
        res.append((sum(int(x[1]) for x in e), sum(int(x[4]) for x in e), float(t[0][1]) + 0.0,
                    max(float(x[5]) for x in e)))
    return res


print("# Replay at --guard=all of every box recounted in #126 against the dense chain of its cache "
      "(D: winding2.py caches, Z: zerilli_zeros.py caches)")
print(f"jobs.txt {'equals' if jobs_ok else 'DIFFERS FROM'} the job list ({len(boxes)} boxes, {nD} of D and "
      f"{len(boxes) - nD} of Z, {len(lines) - len(boxes)} controls skipped)")
print("| job | F | steps of F | steps of F*factor | chain points | count of F | count of F*factor | planted found | "
      "largest increment F | largest increment F*factor | check |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
rows, biggest, nonzero = 0, (0.0, ''), []
for job, script, arguments in boxes:
    kind = 'D' if script == 'winding2.py' else 'Z'
    name, (x0, x1, y0, y1), opt = tag(script, arguments)
    if MUT == 'nocache' and job == MUTJOB:
        name += '.missing'
    cf = os.path.join(RAW, 'cache', name)
    why, cells = [], [''] * 8
    if not os.path.exists(cf):
        why.append(f"cache {name} is missing")
    else:
        sym = y0 == 0.0
        s0 = complex(*(float(v) for v in opt['plant'].split(',')))
        cmd = [sys.executable, os.path.join(HERE, 'diagnose_plant.py'), cf] + arguments.split()[3:6] + [opt['plant']]
        if not sym and MUT != 'noy0':
            cmd.append(f"--y0={opt['y0']}")
        cmd += ['--guard=all', '--expect=none']
        if MUT == 'sparse' and job == MUTJOB:
            cmd.append('--mutate=sparse')
        r = subprocess.run(cmd, capture_output=True, text=True, env=ENV)
        m, chk, res = MISS.search(r.stdout), CHK.search(r.stdout), parse(r.stdout, sym)
        if m:
            why.append(f"the replay stops on a point the cache lacks, {m.group(1)}")
        elif res is None or chk is None:
            why.append(f"diagnose_plant.py exit {r.returncode}, output not parsed: "
                       f"{(r.stdout + r.stderr).strip().splitlines()[-1:]}")
        else:
            if r.returncode != 0 or chk.group(2) != 'PASS':
                why.append(f"replay {chk.group(2)} (exit {r.returncode})")
            (sf, cf_, nf, mf), (sp, cp, npf, mp) = res
            keys = pickle.load(open(cf, 'rb'))
            pts = cf_ + (1 if sym else 0)                      # chain points: edges share their corners
            if cf_ != cp or pts != contour_points(keys, x0, x1, y0, y1):
                why.append(f"chain has {pts} points, the cache {contour_points(keys, x0, x1, y0, y1)} on the contour")
            if MUT == 'misplant' and job == MUTJOB:
                s0 += 10 * (x1 - x0)
            inside = x0 < s0.real < x1 and ((0 < abs(s0.imag) < y1) if sym else (y0 < s0.imag < y1))
            found = round(npf - nf, 4) + 0.0
            if abs(found - (2 if sym else 1) * inside) >= 0.05:
                why.append(f"planted found {found:+.4f}, expected {(2 if sym else 1) * inside}")
            cells = [sf, sp, pts, f"{nf:+.4f}", f"{npf:+.4f}", f"{found:+.4f}", f"{mf:.3f}", f"{mp:.3f}"]
            if max(mf, mp) > biggest[0]:
                biggest = (max(mf, mp), job)
            if abs(nf) >= 0.5:
                nonzero.append(job)
    rows += 1
    print(f"| {job} | {kind} | " + ' | '.join(str(c) for c in cells) + f" | {'PASS' if not why else 'FAIL'} |")
    bad += [f"{job}: {w}" for w in why]
print()
print(f"largest increment of any chain: {biggest[0]:.3f} rad ({biggest[1]}); boxes whose count of F is not 0: "
      f"{', '.join(nonzero) or 'none'}")
print()
print("RESULT: " + (f"all {rows} boxes pass: no accepted step differs from its dense chain, every chain holds all the cached "
                    f"points of its contour, every planted zero is found" if not bad and rows == 61 else
                    "FAILED: " + "; ".join(bad or [f"{rows} boxes, not 61"])) + (f"  [mutate={MUT}]" if MUT else ""))
sys.exit(1 if bad or rows != 61 else 0)

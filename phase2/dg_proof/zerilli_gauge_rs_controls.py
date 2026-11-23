"""#146: runs zerilli_gauge_rs.py unmutated (must exit 0, no FAIL) and once per mutation with --mutate=NAME --failfast (must exit with the
listed code and the first FAIL id equal to the target written in I146-rpi5/README.md before the first run of the script).

Run: python zerilli_gauge_rs_controls.py > zerilli_gauge_rs_controls_output.txt
Environment: CTRL_LOGDIR = directory for the per-run logs (optional). Exit 0 iff the baseline and every pair is OK.
Pairs whose target lies after PL.G and PL.X run with --skip=PL.G,PL.X (those two planted checks take about 3 minutes and the mutation acts
later; the unmutated baseline runs them).
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
LOGDIR = os.environ.get('CTRL_LOGDIR')
SKIP = ['--skip=PL.G,PL.X']

# (mutation, extra args, target first-FAIL id, expected exit code): exit 2 = a planted check, 1 = a real check
PAIRS = [
    ('dual_inv', SKIP, 'PL.0a', 2),
    ('ds_typo', SKIP, 'PL.0b', 2),
    ('bg_nonvac', SKIP, 'PL.0c', 2),
    ('leg_sign', SKIP, 'PL.0d', 2),
    ('lie_drop', SKIP, 'PL.L', 2),
    ('drop_trace', ['--skip=PL.X'], 'PL.G', 2),
    ('route2_drop', ['--skip=PL.G'], 'PL.X', 2),
    ('ext_G', SKIP, 'PL.E', 2),
    ('solve_sign', SKIP, 'PL.P', 2),
    ('solve_drop_s', SKIP, 'PL.P', 2),
    ('l2_struct', SKIP, 'GA.1', 1),
    ('gs_a', SKIP, 'GS.1', 1),
    ('hand_typo', SKIP, 'GS.2', 1),
    ('gdiv_set', SKIP, 'GS.3', 1),
    ('gdiv_s', SKIP, 'GS.3', 1),
    ('l2', SKIP, 'GU.1', 1),
    ('fam_typo', SKIP, 'GU.2', 1),
    ('st_Edep', SKIP, 'GST.1', 1),
    ('st_l2', SKIP, 'GST.2', 1),
    ('kill_dropG', SKIP, 'GST.3', 1),
    ('flow_typo', SKIP, 'R.1', 1),
    ('rs_typo', SKIP, 'R.2', 1),
    ('pole_inject', SKIP, 'R.3', 1),
    ('flow_h', SKIP, 'R.4', 1),
    ('cz_coef', SKIP, 'R.5', 1),
    ('cz_sign', SKIP, 'R.5', 1),
    ('rw17_typo', SKIP, 'R.6', 1),
    ('nv1_h1', SKIP, 'R.7', 1),
    ('nv1_h2', SKIP, 'R.7', 1),
    ('im_coef', SKIP, 'R.8', 1),
    ('n2_sign', SKIP, 'R.9', 1),
    ('log_inject', SKIP, 'R.10', 1),
    ('rs6c_typo', SKIP, 'R.11', 1),
    ('sflow_typo', SKIP, 'R.S.1', 1),
    ('sflow_h', SKIP, 'R.S.2', 1),
]


def run_one(extra):
    t = time.time()
    p = subprocess.run([sys.executable, '-I', os.path.join(HERE, 'zerilli_gauge_rs.py')] + extra,
                       cwd=HERE, capture_output=True, text=True, timeout=3600, env=dict(os.environ, OMP_NUM_THREADS='1'))
    return p.returncode, p.stdout, p.stderr, time.time() - t


def first_fail(out):
    for ln in out.splitlines():
        if ln.startswith('FAIL '):
            return ln[5:].split(':', 1)[0]
    return None


nok = nmiss = 0
print('# baseline: unmutated must exit 0 with no FAIL line and empty stderr', flush=True)
code, out, err, dt = run_one([])
nf = sum(1 for ln in out.splitlines() if ln.startswith('FAIL '))
npass = sum(1 for ln in out.splitlines() if ln.startswith('PASS '))
ok = code == 0 and nf == 0 and not err.strip()
print('%s baseline exit=%d PASS=%d FAIL=%d %.0fs' % ('OK  ' if ok else 'MISS', code, npass, nf, dt), flush=True)
if LOGDIR:
    open(os.path.join(LOGDIR, 'baseline.log'), 'w').write(out + err)
nok += ok
nmiss += not ok
print('# mutated: --mutate=NAME --failfast; target = first FAIL id, exit code as listed', flush=True)
for mu, extra, target, want_exit in PAIRS:
    code, out, err, dt = run_one(['--mutate=' + mu, '--failfast'] + extra)
    ff = first_fail(out)
    ok = code == want_exit and ff == target
    print('%s %-12s target=%-7s exit=%d (want %d) first-FAIL=%-7s %.0fs' % ('OK  ' if ok else 'MISS', mu, target, code, want_exit, ff, dt), flush=True)
    if LOGDIR:
        open(os.path.join(LOGDIR, '%s.log' % mu), 'w').write(out + err)
    nok += ok
    nmiss += not ok
print('# summary: %d OK, %d MISS (of 1 baseline + %d pairs)' % (nok, nmiss, len(PAIRS)), flush=True)
sys.exit(1 if nmiss else 0)

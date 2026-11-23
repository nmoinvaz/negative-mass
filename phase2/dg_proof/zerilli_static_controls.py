"""#145: runs zerilli_static.py unmutated (must exit 0, no FAIL) and once per mutation with --mutate=NAME --failfast (must exit with the
listed code and the first FAIL id equal to the target written in I145-rpi5/README.md before the first run of the script).

Run: python zerilli_static_controls.py > zerilli_static_controls_output.txt
Environment: CTRL_LOGDIR = directory for the per-run logs (optional). Exit 0 iff the baseline and every pair is OK.
Pairs whose target lies after PL.G and PL.X run with --skip=PL.G,PL.X (those two planted checks take about 2.7 minutes and the mutation acts
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
    ('dual_inv', [], 'PL.0a', 2),
    ('ds_typo', [], 'PL.0b', 2),
    ('bg_nonvac', [], 'PL.0c', 2),
    ('leg_sign', [], 'PL.0d', 2),
    ('drop_trace', [], 'PL.G', 2),
    ('route2_drop', ['--skip=PL.G'], 'PL.X', 2),
    ('ang_norm_odd', SKIP, 'PL.O.1', 2),
    ('inj_rphi', SKIP, 'PL.O.2', 2),
    ('rw_typo', SKIP, 'PL.O.3', 2),
    ('fsign_ans', SKIP, 'PL.O.3', 2),
    ('cpm_typo', SKIP, 'PL.O.3', 2),
    ('cmp_blind', SKIP, 'PL.O.4', 2),
    ('l1', SKIP, 'A.1', 1),
    ('ang_norm', SKIP, 'A.2', 1),
    ('s_leak', SKIP, 'S.1', 1),
    ('inj_tr', SKIP, 'H.1', 1),
    ('inj_tx', SKIP, 'H.2', 1),
    ('inj_h1', SKIP, 'H.3', 1),
    ('ans_rr', SKIP, 'E.1', 1),
    ('grr_drop', SKIP, 'E.2', 1),
    ('tt_drop_H', SKIP, 'E.4', 1),
    ('psi_zero', SKIP, 'E.5a', 1),
    ('mv_coef', SKIP, 'E.5b', 1),
    ('fsign_mv', SKIP, 'E.5b', 1),
    ('mv_lam', SKIP, 'E.5b', 1),
    ('v_typo', SKIP, 'E.5d', 1),
    ('res_typo', SKIP, 'E.5e', 1),
    ('dg3_typo', SKIP, 'E.6a', 1),
    ('inv_typo', SKIP, 'E.6b', 1),
    ('dep_typo', SKIP, 'S.0', 1),
    ('tr_typo', SKIP, 'S.R1', 1),
    ('h1t_typo', SKIP, 'S.R2', 1),
    ('sr_typo', SKIP, 'S.R3', 1),
    ('h1p_typo', SKIP, 'S.R4', 1),
    ('div_set', SKIP, 'D.1', 1),
    ('d2_zero', SKIP, 'D.2', 1),
    ('div_s', SKIP, 'D.3', 1),
    ('div_l2', SKIP, 'D.4', 1),
    ('d5_zero', SKIP, 'D.5', 1),
]


def run_one(extra):
    t = time.time()
    p = subprocess.run([sys.executable, '-I', os.path.join(HERE, 'zerilli_static.py')] + extra,
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

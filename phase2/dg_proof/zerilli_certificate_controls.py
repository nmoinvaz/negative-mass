"""#148: runs zerilli_certificate.py unmutated (must exit 0, no FAIL) and once per mutation with --mutate=NAME --failfast (must exit with the
listed code and the first FAIL id equal to the target written in negative-mass-raw/I148-rpi5/README.md before the first control run).
The exploratory drop_* family is reported separately: an accepted outcome is exit 0 (the equation is not needed) or exit 1 with the first FAIL E.cert or E.den.

Run: python zerilli_certificate_controls.py > zerilli_certificate_controls_output.txt
Environment: CTRL_LOGDIR = directory for the per-run logs (optional). Exit 0 iff the baseline, every required pair and every exploratory run is as listed.
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
LOGDIR = os.environ.get('CTRL_LOGDIR')

# (mutation, extra args, target first-FAIL id, expected exit code): exit 2 = a planted check, 1 = a real check
PAIRS = [
    ('leg_sign', ['--part=planted'], 'PL.0d', 2),
    ('den_allow_M', ['--part=planted'], 'PL.D.1', 2),
    ('odd_v', ['--part=odd'], 'PL.O.2', 2),
    ('odd_psi', ['--part=odd'], 'PL.O.2', 2),
    ('odd_cert', ['--part=odd'], 'PL.O.3', 2),
    ('odd_den', ['--part=odd', '--attempts=1'], 'PL.O.4', 2),
    ('l1', ['--part=even'], 'A.1', 1),
    ('expl_norm', ['--part=even'], 'S.explicit', 1),
    ('v_coef', ['--part=even'], 'E.cert', 1),
    ('psi_coef', ['--part=even'], 'E.cert', 1),
    ('v_s2', ['--part=even'], 'E.cert', 1),
    ('drop_traceless', ['--part=even'], 'E.cert', 1),
    ('flow_target', ['--part=even', '--attempts=1'], 'E.den', 1),
    ('cert_coef', ['--part=even'], 'E.id', 1),
    ('l_spec', ['--part=even'], 'E.L', 1),
    ('m0_v', ['--part=even'], 'E.M0', 1),
    ('s0_s1', ['--part=even'], 'E.S0', 1),
]
EXPLORATORY = ['drop_tt', 'drop_tr', 'drop_rr', 'drop_tx', 'drop_rx', 'drop_trace']
ACCEPT = {(0, None), (1, 'E.cert'), (1, 'E.den')}


def run_one(extra):
    t = time.time()
    p = subprocess.run([sys.executable, '-I', os.path.join(HERE, 'zerilli_certificate.py')] + extra,
                       cwd=HERE, capture_output=True, text=True, timeout=3600, env=dict(os.environ, OMP_NUM_THREADS='1'))
    return p.returncode, p.stdout, p.stderr, time.time() - t


def first_fail(out):
    for ln in out.splitlines():
        if ln.startswith('FAIL '):
            return ln[5:].split(':', 1)[0]
    return None


def save(name, out, err):
    if LOGDIR:
        os.makedirs(LOGDIR, exist_ok=True)
        with open(os.path.join(LOGDIR, name + '.log'), 'w') as fh:
            fh.write(out)
            if err:
                fh.write('\n# stderr\n' + err)


print('# zerilli_certificate_controls: baseline, %d required pairs, %d exploratory runs' % (len(PAIRS), len(EXPLORATORY)), flush=True)
bad = 0
rc, out, err, dt = run_one([])
ff = first_fail(out)
npass = sum(1 for ln in out.splitlines() if ln.startswith('PASS '))
ok = rc == 0 and ff is None
save('baseline', out, err)
print('%s baseline: exit %d, %d PASS, first FAIL %s (%.0f s)' % ('OK  ' if ok else 'BAD ', rc, npass, ff, dt), flush=True)
bad += 0 if ok else 1

nok = 0
for name, extra, target, want in PAIRS:
    rc, out, err, dt = run_one(['--mutate=' + name, '--failfast'] + extra)
    ff = first_fail(out)
    ok = rc == want and ff == target
    save(name, out, err)
    print('%s %-16s exit %d (want %d), first FAIL %s (target %s) (%.0f s)' % ('OK  ' if ok else 'MISS', name, rc, want, ff, target, dt), flush=True)
    nok += ok
    bad += 0 if ok else 1
print('# required pairs on target: %d of %d' % (nok, len(PAIRS)), flush=True)

nex = 0
for name in EXPLORATORY:
    rc, out, err, dt = run_one(['--mutate=' + name, '--failfast', '--part=even'])
    ff = first_fail(out)
    ok = (rc, ff) in ACCEPT
    save(name, out, err)
    kind = 'silent (not needed)' if rc == 0 else 'fires at %s' % ff
    print('%s %-16s exit %d, first FAIL %s: %s (%.0f s)' % ('OK  ' if ok else 'UNEXPECTED', name, rc, ff, kind, dt), flush=True)
    nex += ok
    bad += 0 if ok else 1
print('# exploratory runs with an accepted outcome: %d of %d' % (nex, len(EXPLORATORY)), flush=True)
print('# summary: %d problem(s)' % bad, flush=True)
sys.exit(1 if bad else 0)

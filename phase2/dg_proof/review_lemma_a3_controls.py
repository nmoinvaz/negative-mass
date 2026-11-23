"""#142: runs every review_lemma_a3 check script unmutated (must exit 0) and once per (mutation, script) pair with
--mutate=NAME --failfast (must exit non-zero, first FAIL id = the target written in I142-rpi5/README.md before the first run).

Run: python review_lemma_a3_controls.py > review_lemma_a3_controls_output.txt
Environment: CTRL_LOGDIR = directory for the per-run logs (optional).
Exit 0 iff every baseline and every pair is OK.
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
LOGDIR = os.environ.get('CTRL_LOGDIR')
PLANTED = 2

# (script, mutation, target first-FAIL id, expected exit code), as in "Targets" of the raw README
PAIRS = [
    ('a3', 'pl_crit_root', 'PLA.1', 2), ('a3', 'pl_sturm_sign', 'PLS.1', 2), ('a3', 'a3_c5', 'A3.1', 1), ('a3', 'cn_n4to3', 'A3.1', 1),
    ('a3', 'a3_K120', 'A3.2', 1), ('a3', 'a3_qsign', 'A3.3[l=3]', 1), ('a3', 'a3_s', 'A3.4', 1), ('a3', 'a3_qclosed', 'A3.5', 1),
    ('a3', 'a3_ident', 'A3.6', 1), ('a3', 'a3_l2_min', 'A3.8', 1), ('a3', 'a3_sturm_neg', 'A3.9[l=2]', 1),
    ('ends', 'pl_indicial', 'PLE.1a', 2), ('ends', 'pl_order', 'PLE.3a', 2), ('ends', 'pl_taylor', 'PLE.4a', 2), ('ends', 'pl_vhat', 'PLE.2[n=1]', 2),
    ('ends', 'e_const12', 'E2.1', 1), ('ends', 'e_linear', 'E2.1', 1), ('ends', 'cn_n4to3', 'E2.2[l=2]', 1), ('ends', 'e_den', 'E2.3[l=2]', 1),
    ('ends', 'e_exp', 'E2.4', 1), ('ends', 'e_phi', 'E2.7[l=2]', 1), ('ends', 'e_nn1', 'E3.1', 1), ('ends', 'e_lead_inf', 'E3.2[l=2]', 1),
    ('dg61', 'pl_taylor', 'PLE.4a', 2), ('dg61', 'e_dg61b', 'E2.6b', 1), ('dg61', 'cn_n4to3', 'E2.6c', 1),
    ('l6', 'pl_wronsk', 'PLW.1', 2), ('l6', 'pl_order', 'PLE.3a', 2), ('l6', 'pl_sturm_sign', 'PLS.1', 2), ('l6', 'l6_gen', 'L6.1', 1),
    ('l6', 'l6_w', 'L6.2', 1), ('l6', 'l6_c4', 'L6.3', 1), ('l6', 'cn_n4to3', 'L6.4[l=2]', 1), ('l6', 'l6_k', 'L6.4[l=2]', 1),
    ('l6', 'l6_rs', 'L6.5', 1), ('l6', 'l6_chip', 'L6.6[l=2]', 1),
    ('l4', 'pl_vhat', 'PLL.1', 2), ('l4', 'pl_eig', 'PLL.2', 2), ('l4', 'l4_alg', 'L4.1', 1), ('l4', 'v_typo', 'L4.2', 1),
    ('l4', 'l4_k', 'L4.2', 1), ('l4', 'cn_n4to3', 'L4.4[l=2]', 1), ('l4', 'l4_sign', 'L4.5[l=2]', 1),
    ('v', 'pl_cmp', 'PLV.1b', 2), ('v', 'v_typo', 'V1.1[DG06]', 1), ('v', 'v_dg06', 'V1.1[DG06]', 1), ('v', 'v_ks', 'V1.1[KS]', 1),
    ('v', 'v_x', 'V1.3', 1), ('v', 'v_chand', 'V1.4', 1), ('v', 'v_dg54', 'V1.5', 1), ('v', 'v_dg53', 'V1.6', 1), ('v', 'v_cn', 'V1.8', 1),
    ('v', 'cn_n4to3', 'V1.8', 1), ('v', 'v_m2', 'V1.9', 1),
]
SCRIPTS = ['a3', 'ends', 'dg61', 'l6', 'l4', 'v']


def run_one(script, extra):
    t = time.time()
    p = subprocess.run([sys.executable, '-I', os.path.join(HERE, 'review_lemma_a3_%s.py' % script)] + extra,
                       cwd=HERE, capture_output=True, text=True, timeout=7200, env=dict(os.environ, OMP_NUM_THREADS='1'))
    return p.returncode, p.stdout, p.stderr, time.time() - t


def first_fail(out):
    for ln in out.splitlines():
        if ln.startswith('FAIL '):
            return ln[5:].split(':', 1)[0]
    return None


nok = nmiss = 0
print('# baselines: each script unmutated must exit 0 with no FAIL line', flush=True)
for sc in SCRIPTS:
    code, out, err, dt = run_one(sc, [])
    nf = sum(1 for ln in out.splitlines() if ln.startswith('FAIL '))
    nd = sum(1 for ln in out.splitlines() if ln.startswith('DIFFERS '))
    want_d = 1 if sc == 'dg61' else 0
    ok = code == 0 and nf == 0 and nd == want_d and not err.strip()
    print('%s baseline %-5s exit=%d FAIL=%d DIFFERS=%d (want %d) %.0fs' % ('OK  ' if ok else 'MISS', sc, code, nf, nd, want_d, dt), flush=True)
    if LOGDIR:
        open(os.path.join(LOGDIR, 'baseline_%s.log' % sc), 'w').write(out + err)
    nok += ok
    nmiss += not ok
print('# mutated: --mutate=NAME --failfast; target = first FAIL id, exit code as listed', flush=True)
for sc, mu, target, want_exit in PAIRS:
    code, out, err, dt = run_one(sc, ['--mutate=' + mu, '--failfast'])
    ff = first_fail(out)
    ok = code == want_exit and ff == target
    print('%s %-5s %-14s target=%-12s exit=%d (want %d) first-FAIL=%-12s %.0fs' % ('OK  ' if ok else 'MISS', sc, mu, target, code, want_exit, ff, dt), flush=True)
    if LOGDIR:
        open(os.path.join(LOGDIR, '%s_%s.log' % (sc, mu)), 'w').write(out + err)
    nok += ok
    nmiss += not ok
print('# summary: %d OK, %d MISS (of %d baselines + %d pairs)' % (nok, nmiss, len(SCRIPTS), len(PAIRS)), flush=True)
sys.exit(1 if nmiss else 0)

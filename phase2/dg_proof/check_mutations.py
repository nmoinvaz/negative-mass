"""
Mutation driver for the #2 scripts.  For each script, run it unmutated (must exit 0) and once per mutation (must exit
non-zero), one process at a time.  Every check in those scripts can therefore fail.
Usage: python3 check_mutations.py [script ...]     (default: all)
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

PLAN = {
    # script: (extra args for the unmutated run, [mutations])
    'vhat_structure.py': ([], ['vcoef', 'qsign', 'ssign', 'toy']),
    'threshold_isolation.py': (['--lmax', '60'], ['coef', 'sterm', 'ntail', 'mono', 'planted']),
    'large_l_proof.py': ([], ['beta', 'eta', 'dsign', 'poly', 'window']),
    'monotone_proof.py': ([], ['omega', 'n1sign', 's1sign', 'eta2', 'window']),
}

want = sys.argv[1:] or list(PLAN)
bad = []
for script in want:
    if script not in PLAN or not os.path.exists(os.path.join(HERE, script)):
        print(f"skip {script} (not present)")
        continue
    extra, muts = PLAN[script]
    env = dict(os.environ, OMP_NUM_THREADS='1')
    rc = subprocess.run(['nice', '-n', '10', PY, os.path.join(HERE, script)] + extra, capture_output=True, text=True, env=env).returncode
    print(f"{script:28s} unmutated: exit {rc}  {'OK' if rc == 0 else 'UNEXPECTED'}", flush=True)
    if rc != 0:
        bad.append((script, 'unmutated'))
    for mu in muts:
        rc = subprocess.run(['nice', '-n', '10', PY, os.path.join(HERE, script)] + extra + [f'--mutate={mu}'],
                            capture_output=True, text=True, env=env).returncode
        print(f"{script:28s} --mutate={mu:8s}: exit {rc}  {'fires' if rc != 0 else 'DOES NOT FIRE'}", flush=True)
        if rc == 0:
            bad.append((script, mu))
print("all mutations fire and every unmutated run passes" if not bad else f"PROBLEMS: {bad}")
sys.exit(1 if bad else 0)

"""
Mutation driver for the #3 scripts.  Each script is run unmutated (must exit 0) and once per mutation (must exit non-zero), one process
at a time, nice -n 10, one thread.  Every check in these scripts can therefore fail.
Usage: python3 check_mutations.py [script ...]     (default: all)
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

PLAN = {
    # script: (extra args, [mutations])
    'check_modes.py': ([], ['chi', 'omega', 'jac', 'vhat', 'nodes', 'pow']),
    'check_fastmodes.py': ([], ['coef', 'xr']),
    'planted.py': ([], ['omega', 'signV', 'label', 'kappa', 'dt']),
    'hadamard_sequence.py': (['--lfd', '14'], ['omega', 'vhat', 'sign', 'data', 'dt']),
    'star_control.py': (['--lmax', '24', '--lpm', '16'], ['omega0', 'kappa', 'rs', 'sign', 'bound', 'nodes', 'kl']),
    'causal.py': ([], ['x1', 'omega', 'bound', 'thr']),
}

want = sys.argv[1:] or list(PLAN)
bad = []
env = dict(os.environ, OMP_NUM_THREADS='1')
for script in want:
    extra, muts = PLAN[script]
    t0 = time.time()
    rc = subprocess.run(['nice', '-n', '10', PY, os.path.join(HERE, script)] + extra, capture_output=True, text=True, env=env).returncode
    print(f"{script:24s} unmutated: exit {rc}  {'OK' if rc == 0 else 'UNEXPECTED'}  ({time.time() - t0:.0f} s)", flush=True)
    if rc != 0:
        bad.append((script, 'unmutated'))
    for mu in muts:
        t0 = time.time()
        rc = subprocess.run(['nice', '-n', '10', PY, os.path.join(HERE, script)] + extra + [f'--mutate={mu}'],
                            capture_output=True, text=True, env=env).returncode
        print(f"{script:24s} --mutate={mu:8s}: exit {rc}  {'fires' if rc != 0 else 'DOES NOT FIRE'}  ({time.time() - t0:.0f} s)", flush=True)
        if rc == 0:
            bad.append((script, mu))
print("all mutations fire and every unmutated run passes" if not bad else f"PROBLEMS: {bad}")
sys.exit(1 if bad else 0)

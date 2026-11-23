#!/usr/bin/env python3
"""Controls of the independent review (issue 104): runs every review_large_l_*.py once per `--mutate=NAME`.
Each mutated run must exit non-zero and its FIRST failing check must be the one named in TARGETS (written before the first control run,
copied into the raw README).  A silent mutation (exit 0), a crash, or a different first failure is reported as such and makes this script exit 1.
Usage: review_large_l_controls.py [--only=SCRIPTSUBSTRING] [--logdir=DIR]
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

TARGETS = {
    "review_large_l_common.py": {
        "coef_q": "coefficients a_n(eps) in Lemma A's closed form, l=2",
        "coef_lin": "coefficients a_n(eps) in Lemma A's closed form, l=2",
        "nz_coef": "NZ(rho, eps) matches sympy at a point",
        "accept_zero_enclosure": "planted: endpoint exactly on a root",
        "no_deriv": "planted: three roots in one bracket",
        "alpha_n": "alpha_1 = 3/2, alpha_2 = 0",
        "drop_last": "coefficients a_n(eps) in Lemma A's closed form, l=2",
        "psip_term": "Psi' formula equals sympy",
    },
    "review_large_l_roots.py": {
        "bessel_coef": "Bessel closed form: coefficient of s^n equals",
        "alpha_n": "series of X_inf = Bessel closed form at rho=3.3",
        "nz_coef": "rho_neg reproduces RESULTS.md's 9.635989442586144",
        "shift_enclosure": "rho_*: iv enclosure",
        "no_deriv": "planted: three roots in one bracket",
        "accept_zero_enclosure": "planted: endpoint exactly on a root",
    },
    "review_large_l_theorem3.py": {
        "mut_Z": "(Z) psi_0,xx = V psi_0 exactly, l=2",
        "mut_V_NZ": "V(r = rho/c, lambda = c/2) = 8c^4",
        "mut_Vhat": "2 g^2 - V = 2c^4 (2+eps rho) Psi",
        "mut_S": "Vhat = S c^3 (2c + rho)",
        "mut_P1": "(P1) numerator of V",
        "mut_P2": "(P2) c^5 Q(tau, m) has all coefficients positive",
        "mut_lead": "lead coefficient of Psi",
        "mut_count_psi": "l=2: Psi has exactly one root",
        "mut_count_N": "l=2: N has exactly one root in (6, oo)",
        "mut_NZpos": "N_Z(6 + tau) has all coefficients positive",
    },
    "review_large_l_brackets.py": {
        "coef_q": "coefficients a_n(eps), n <= 40, inside Lemma A's exact closed form, l=29",
        "coef_lin": "coefficients a_n(eps), n <= 40, inside Lemma A's exact closed form, l=29",
        "drop_last": "coefficients a_n(eps), n <= 40, inside Lemma A's exact closed form, l=29",
        "nz_coef": "S shows the same sign change (-,+) on the Psi bracket (does not use (Z)), l=29",
        "S_flip": "S shows the same sign change (-,+) on the Psi bracket (does not use (Z)), l=29",
        "accept_zero_enclosure": "planted: endpoint exactly on a root",
        "no_deriv": "planted: three roots in one bracket",
        "shift_bracket": "N: bracket of c r_* certified (signs, derivative), l=29",
        "rho_limit": "c r_* above its limit for every certified l",
        "l2r_reverse": "l^2 r with c r_*: increases along all sampled l >= 5",
    },
    "review_large_l_thm_l.py": {
        "beta_scale": "step 1: |a_n - alpha_n| <= eps B_n for 3 <= n <= l+10, l=138",
        "slope_sign": "step 3: -N_inf' >= mu_N > 0 on the N window",
        "window_small": "step 4: eps0 E_N < mu_N delta_N",
        "psi_identity": "Psi_l - Psi_inf = eps [rho^3 N_l^2",
        "bracket_k": "Theorem L against the certified brackets, l=138",
        "eta_scale": "step 2 at 11 points of each window, l=138",
    },
    "review_large_l_thm_m.py": {
        "nu_scale": "L2: expansions at 5 points of each window against the exact N_l and Psi_l, l=138",
        "omega_scale": "L2: expansions at 5 points of each window against the exact N_l and Psi_l, l=138",
        "cmp_reverse": "M: (2l+2) nu_min >",
        "l2r_s": "M: l^2 r_* increases: s_N < rho_*",
        "window_narrow": "L2: s_N eps_max <= window radius",
        "brk_s": "Theorem M against my certified brackets, l=500 and 501",
    },
    "review_large_l_table.py": {
        "table_digit": "RESULTS table row l=20",
    },
    "review_large_l_sweep.py": {
        "coef_q": "coefficients a_n(eps), n <= 40, inside Lemma A's exact closed form, l=2",
        "shift_bracket": "N: bracket of c r_* certified (signs, derivative), l=2",
        "S_flip": "S shows the same sign change (-,+) on the Psi bracket (does not use (Z)), l=2",
        "rho_limit": "bracket above its limit, l=",
        "l2r_reverse": "l^2 r increases (N), l=5 -> 6",
    },
}


CSV = os.path.expanduser("~/Source/negative-mass-raw/I104-rpi5/sweep_brackets.csv")


def run(script, mut, logdir):
    t0 = time.time()
    extra = [f"--csv={CSV}"] if script == "review_large_l_table.py" else []
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    p = subprocess.run(["nice", "-n", "10", PY, os.path.join(HERE, script), f"--mutate={mut}"] + extra, cwd=HERE, env=env,
                       capture_output=True, text=True)
    out = p.stdout + ("\n[stderr]\n" + p.stderr if p.stderr.strip() else "")
    if logdir:
        with open(os.path.join(logdir, f"{script[:-3]}__{mut}.log"), "w") as f:
            f.write(out)
    first = next((ln for ln in p.stdout.splitlines() if ln.startswith("FAIL ")), None)
    if first is None:
        first = next((ln for ln in p.stdout.splitlines() if "PLANTED PROBLEM FAILED" in ln), None)
    if first is None and p.returncode != 0:
        first = "(no FAIL line) " + (p.stderr.strip().splitlines()[-1] if p.stderr.strip() else "")
    return p.returncode, first, time.time() - t0


def main():
    only = None
    logdir = None
    for a in sys.argv[1:]:
        if a.startswith("--only="):
            only = a.split("=", 1)[1]
        elif a.startswith("--logdir="):
            logdir = a.split("=", 1)[1]
            os.makedirs(logdir, exist_ok=True)
    print("review_large_l_controls: every mutation must exit non-zero and fire at its target")
    rows = []
    bad = 0
    for script, muts in TARGETS.items():
        if only and only not in script:
            continue
        for mut, target in muts.items():
            rc, first, dt = run(script, mut, logdir)
            fired = rc != 0
            on_target = fired and first is not None and target in first
            verdict = "OK" if on_target else ("SILENT" if not fired else "OFF-TARGET")
            bad += 0 if on_target else 1
            rows.append((script, mut, rc, verdict, target, first))
            print(f"{verdict:10s} {script:30s} --mutate={mut:24s} exit={rc} {dt:6.1f}s", flush=True)
            print(f"           target: {target}")
            print(f"           first FAIL: {(first or '')[:170]}", flush=True)
    print(f"SUMMARY mutations={len(rows)} on-target={len(rows) - bad} not-on-target={bad}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

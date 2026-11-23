#!/usr/bin/env python3
"""Controls of the #138 review: every script unmutated (must exit 0) and each (mutation, script) pair of the raw README (must fail on target).

    python3 review_lemma_a_controls.py [--logdir=DIR] [--only=MUTATION[,MUTATION]]

The target table below is copied from the raw README (negative-mass-raw/I138-rpi5/README.md, "Mutations and targets"), written before the first control run.
A pair passes if the mutated run exits non-zero without a traceback and its first FAIL token matches the target. Exit 1 if any pair misses.
"""
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["common", "signs", "recur", "zerilli", "l1", "t1", "t2"]

TARGETS = {
    "plant_recur": [("common", "PL1"), ("recur", "PL1")],
    "plant_exp": [("common", "PL3"), ("l1", "PL3")],
    "plant_rel": [("common", "PL4"), ("l1", "PL4")],
    "plant_sc": [("common", "PL5"), ("t2", "PL5")],
    "plant_min": [("common", "PL6"), ("t1", "PL6")],
    "plant_zer": [("common", "PL7"), ("zerilli", "PL7")],
    "cn_n4to3": [("signs", "A1.1"), ("recur", "A2.4"), ("zerilli", "Z1[l=2]"), ("l1", "L1.6[l=2]"), ("t1", "T1.2[l=2]"), ("t2", "T2.3[l=2]")],
    "cn_poch": [("signs", "A1.1")],
    "cn_sign5": [("signs", "A1.4[l=3]"), ("recur", "A2.7[l=3]"), ("zerilli", "Z1[l=3]"), ("t2", "T2.3[l=3]")],
    "V_coeff": [("recur", "A2.1"), ("zerilli", "Z1[l=2]")],
    "P2_claim": [("recur", "A2.2")],
    "P1_6c_5c": [("recur", "A2.1"), ("l1", "L1.1a")],
    "P0_12to13": [("recur", "A2.1"), ("l1", "L1.1b")],
    "rec_12to11": [("recur", "A2.3")],
    "rec_m4to3": [("recur", "A2.3")],
    "poch_ratio": [("recur", "A2.5")],
    "ratio_claim": [("l1", "L1.2")],
    "t1_coef18": [("t1", "T1.1")],
    "t1_drop_sum2": [("t1", "T1.1")],
    "t1_bound3": [("t1", "T1.3")],
    "t2_shift": [("t2", "T2.1")],
    "t2_flip_d0": [("t2", "T2.5[l=2]")],
    "t2_window": [("t2", "T2.6[l=2]")],
}


def matches(want, tok):
    if "[" in want:
        return tok == want
    base = tok.split("[")[0]
    return base == want or (base.startswith(want) and re.fullmatch(r"[a-z]", base[len(want):] or "?") is not None)


def run(script, mut, logdir):
    cmd = [sys.executable, "-u", os.path.join(HERE, f"review_lemma_a_{script}.py")] + ([f"--mutate={mut}"] if mut else [])
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=HERE, env=dict(os.environ, OMP_NUM_THREADS="1"))
    dt = time.time() - t0
    if logdir:
        with open(os.path.join(logdir, f"{script}_{mut or 'unmutated'}.log"), "w") as fh:
            fh.write(p.stdout + ("\n--- stderr ---\n" + p.stderr if p.stderr else "") + f"\nexit {p.returncode}\n")
    first = None
    for line in p.stdout.splitlines():
        if line.startswith("FAIL "):
            # the id is the first word, plus the [l=N] suffix that check_range appends at the end of the (spaced) id
            first = re.split(r"[ :]", line[5:], maxsplit=1)[0]
            br = re.search(r"\[l=\d+\](?=: )", line)
            if br and "[" not in first:
                first += br.group(0)
            break
    return p.returncode, first, "Traceback" in p.stderr, dt


def main():
    logdir = None
    only = None
    for a in sys.argv[1:]:
        if a.startswith("--logdir="):
            logdir = a.split("=", 1)[1]
            os.makedirs(logdir, exist_ok=True)
        if a.startswith("--only="):
            only = set(a.split("=", 1)[1].split(","))
    misses = []
    print("unmutated runs (target: exit 0, no FAIL)")
    base_ok = True
    for s in SCRIPTS:
        if only:
            break
        rc, first, crash, dt = run(s, None, logdir)
        ok = rc == 0 and first is None and not crash
        base_ok &= ok
        print(f"  {'OK  ' if ok else 'MISS'} {s:8s} exit {rc} first FAIL {first} ({dt:.0f} s)")
        if not ok:
            misses.append((s, "unmutated"))
    print("mutated runs (target: exit != 0, first FAIL as in the raw README)")
    n = nfire = 0
    for mut, pairs in TARGETS.items():
        if only and mut not in only:
            continue
        for s, want in pairs:
            rc, first, crash, dt = run(s, mut, logdir)
            n += 1
            ok = rc != 0 and not crash and first is not None and matches(want, first)
            nfire += ok
            why = "" if ok else (" SILENT" if rc == 0 else " CRASH" if crash else " WRONG TARGET")
            print(f"  {'OK  ' if ok else 'MISS'} {mut:13s} {s:8s} exit {rc} first FAIL {first} want {want}{why} ({dt:.0f} s)")
            if not ok:
                misses.append((s, mut))
    print(f"summary: {nfire} of {n} mutated runs fired on target; unmutated runs {'all exit 0' if base_ok else 'NOT all clean'}")
    if misses:
        print("misses:", ", ".join(f"{m}@{s}" for s, m in misses))
    sys.exit(1 if misses else 0)


if __name__ == "__main__":
    main()

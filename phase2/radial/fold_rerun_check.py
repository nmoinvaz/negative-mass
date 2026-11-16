#!/usr/bin/env python3
"""Check that a rerun of fold_scan.py reproduces the stored values of the families, before the family is extended.

    python fold_rerun_check.py p fold_p.txt rerun_p.log          # canonical: row by row
    python fold_rerun_check.py m rerun_m.log strong_m.txt low_m.txt   # phantom: stored rows against the rerun
    python fold_rerun_check.py --selftest
    python fold_rerun_check.py p fold_p.txt rerun_p.log --mutate=shift

p: stored fold_p.txt is the output of `fold_scan.py 1 0.33 <Nmax>`. Every row of the rerun is compared with the
row of the same index (w, N, M, phic, alpha_c, max|2m/r|), and the end line ("stuck ...") must agree.
m: the stored phantom values (strong_m.txt from strongfield.py in alpha_c, low_m.txt from lowscan.py in N) were not
made by fold_scan.py. Each stored row (N, w, M) is compared with the rerun of `fold_scan.py -1 0.5 <Nmax>`
interpolated in ln N (PCHIP) on N, so the check is on w and M at the stored N.

The tolerance is TOL, relative to each value, a little above the printed precision (5 to 6 digits). Exit 0 only if
every compared value agrees. --mutate=shift scales the stored M by 1 + 1e-3, --mutate=drop deletes a stored row (p),
--mutate=shiftw scales the stored w by 1 + 1e-4 (m), --mutate=end replaces the stored end line, and each must make the check fail. --selftest runs planted rows
through the same code, once as planted (must pass) and once per mutation (must fail).
"""
import math
import re
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator

TOL = 2e-5
KEYS = ("w", "N", "M", "phic", "alpha_c", "max2mr")
FOLD = re.compile(r"eps=([+-]\d) w=([\d.]+) N=([\d.]+) M=([+-][\d.]+) phic=([\d.]+) alpha_c=([\d.]+) "
                  r"max\|2m/r\|=([\d.]+) mesh=(\d+) R=(\d+)")
STRONG = re.compile(r"eps=([+-]\d) alpha_c=([\d.]+) w=([\d.]+) N=([\d.]+) M=([+-][\d.]+)")
LOW = re.compile(r"eps=([+-]\d) N=([\d.]+) w=([\d.]+) M=([+-][\d.]+)")


def parse_fold(text):
    rows, end = [], None
    for ln in text.splitlines():
        m = FOLD.match(ln)
        if m:
            rows.append({"eps": int(m[1]), "w": float(m[2]), "N": float(m[3]), "M": float(m[4]),
                         "phic": float(m[5]), "alpha_c": float(m[6]), "max2mr": float(m[7]),
                         "mesh": int(m[8]), "R": int(m[9])})
        elif ln.startswith("stuck"):
            end = ln.strip()
    return rows, end


def parse_stored(text):
    """Stored phantom rows from strongfield.py (alpha_c, w, N, M) or lowscan.py (N, w, M)."""
    out = []
    for ln in text.splitlines():
        m = STRONG.match(ln)
        if m:
            out.append({"w": float(m[3]), "N": float(m[4]), "M": float(m[5]), "alpha_c": float(m[2])})
            continue
        m = LOW.match(ln)
        if m:
            out.append({"w": float(m[3]), "N": float(m[2]), "M": float(m[4])})
    return out


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def check_p(stored_text, rerun_text, mutate=None):
    s, s_end = parse_fold(stored_text)
    r, r_end = parse_fold(rerun_text)
    if mutate == "shift":
        for x in s:
            x["M"] *= 1 + 1e-3
    if mutate == "drop":
        del s[len(s) // 2]
    if mutate == "end":
        s_end = "stuck bvp(N) failed N=42.7326858 eps=1: A singular Jacobian encountered when solving the collocation system."
    worst, where, mesh = 0.0, None, 0
    fails = []
    if len(s) != len(r):
        fails.append(f"row count {len(s)} stored vs {len(r)} rerun")
    for i, (a, b) in enumerate(zip(s, r)):
        for k in KEYS:
            d = rel(a[k], b[k])
            if d > worst:
                worst, where = d, (i, k)
            if d > TOL:
                fails.append(f"row {i} {k}: stored {a[k]} rerun {b[k]} rel {d:.2e}")
        mesh += a["mesh"] != b["mesh"]
    def end_key(e):
        m = re.match(r"stuck bvp\(N\) failed N=([\d.]+) eps=(\S+): (.*)", e or "")
        return (float(m[1]), m[2], m[3]) if m else None
    a, b = end_key(s_end), end_key(r_end)
    end_same = (a is None and b is None) or (a is not None and b is not None and rel(a[0], b[0]) <= TOL
                                              and a[1:] == b[1:])
    if not end_same:
        fails.append(f"end line differs: stored {s_end!r} rerun {r_end!r}")
    info = {"rows compared": min(len(s), len(r)), "worst rel": worst, "worst at": where,
            "rows with a different mesh count": mesh}
    return not fails, info, fails


def check_m(rerun_text, stored_texts, mutate=None):
    r, _ = parse_fold(rerun_text)
    stored = [x for t in stored_texts for x in parse_stored(t)]
    if mutate == "shift":
        for x in stored:
            x["M"] *= 1 + 1e-3
    if mutate == "shiftw":
        for x in stored:
            x["w"] *= 1 + 1e-4
    lnN = np.log([x["N"] for x in r])
    fw = PchipInterpolator(lnN, [x["w"] for x in r])
    fM = PchipInterpolator(lnN, [x["M"] for x in r])
    fails, worst, where = [], 0.0, None
    inside = 0
    for x in stored:
        if not (lnN[0] <= math.log(x["N"]) <= lnN[-1]):
            continue
        inside += 1
        for k, f in (("w", fw), ("M", fM)):
            v = float(f(math.log(x["N"])))
            d = rel(x[k], v)
            if d > worst:
                worst, where = d, (x["N"], k)
            if d > TOL:
                fails.append(f"N={x['N']}: {k} stored {x[k]} rerun {v:.6g} rel {d:.2e}")
    if inside != len(stored):
        fails.append(f"{len(stored) - inside} stored rows lie outside the rerun's N range")
    info = {"stored rows": len(stored), "inside the rerun's N range": inside, "worst rel": worst,
            "worst at": where}
    if inside < 8:
        fails.append(f"only {inside} stored rows inside the rerun's N range")
    return not fails, info, fails


def report(name, ok, info, fails):
    print(("PASS" if ok else "FAIL") + f" {name}: " + ", ".join(f"{k} {v:.3g}" if isinstance(v, float) else f"{k} {v}"
                                                              for k, v in info.items()))
    for f in fails[:10]:
        print("   ", f)
    return ok


def planted_rows(eps=-1, n=30):
    """Synthetic fold_scan lines with a known smooth w(N), M(N), written in the real print format."""
    rows = []
    for i in range(n):
        N = 0.3 * math.exp(0.05 * i * 1.3)
        w = 0.5 - 0.1 * (1 - math.exp(-N / 8)) + 0.002 * N
        M = -0.5 * N - 0.01 * N * N
        rows.append(f"eps={eps:+d} w={w:.6f} N={N:.5f} M={M:+.6f} phic=0.019000 alpha_c={1 + 0.01 * i:.4f} "
                    f"max|2m/r|={0.01 * i:.3f} mesh={1500 + i} R=60")
    return "\n".join(rows) + "\nstuck bvp(N) failed N=42 eps=1: The maximum number of mesh nodes is exceeded.\n"


def planted_stored(rerun_text):
    """Stored lines as lowscan.py would print them, taken from the rerun's own curve between its rows."""
    r, _ = parse_fold(rerun_text)
    lnN = np.log([x["N"] for x in r])
    fw = PchipInterpolator(lnN, [x["w"] for x in r]); fM = PchipInterpolator(lnN, [x["M"] for x in r])
    lines = []
    for N in np.exp(np.linspace(lnN[1], lnN[-2], 11)):
        lines.append(f"eps=-1 N={N:.5f} w={float(fw(math.log(N))):.6f} M={float(fM(math.log(N))):+.6f} sigma2=['1e-2']")
    return "\n".join(lines)


def selftest():
    text = planted_rows()
    ok = True
    p_ok, info, f = check_p(text, text)
    ok &= report("selftest p as planted (must pass)", p_ok, info, f)
    for mu in ("shift", "drop", "end"):
        fired = not check_p(text, text, mu)[0]
        print(f"selftest p --mutate={mu}: {'FIRED' if fired else 'DID NOT FIRE'}")
        ok &= fired
    stored = planted_stored(text)
    m_ok, info, f = check_m(text, [stored])
    ok &= report("selftest m as planted (must pass)", m_ok, info, f)
    for mu in ("shift", "shiftw"):
        fired = not check_m(text, [stored], mu)[0]
        print(f"selftest m --mutate={mu}: {'FIRED' if fired else 'DID NOT FIRE'}")
        ok &= fired
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    mutate = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mutate=")), None)
    if "--selftest" in sys.argv:
        if not selftest():
            sys.exit(1)
        if not args:
            return
    if args[0] == "p":
        ok, info, fails = check_p(open(args[1]).read(), open(args[2]).read(), mutate)
    else:
        ok, info, fails = check_m(open(args[1]).read(), [open(a).read() for a in args[2:]], mutate)
    report(" ".join(args) + (f" --mutate={mutate}" if mutate else ""), ok, info, fails)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

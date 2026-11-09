#!/usr/bin/env python3
"""Compare the roots of two lowscan.py logs row by row, matched by N.

    python cmp_lowscan.py COARSE.log FINE.log [--rtol=2e-5] [--selftest]

A row agrees when both logs list the same number of roots, with the same '?' flags, each within rtol relative.
The default 2e-5 allows a rounding step in the printed 6th digit.
Exit 1 if any N is missing from either log or any row disagrees.
--selftest drops the lowest root of FINE's first row and requires the check to fail on it.
"""
import re
import sys


def rows(path):
    out = {}
    for line in open(path, encoding="utf-8", errors="replace").read().splitlines():
        m = re.match(r"eps=([+-]\d) N=([\d.]+) .*sigma2=\[([^\]]*)\]", line)
        if m:
            roots = [(float(x.rstrip("?")), x.endswith("?")) for x in re.findall(r"'([^']*)'", m.group(3))]
            out[(m.group(1), m.group(2))] = roots
    return out


def compare(a, b, rtol):
    bad = []
    for key in sorted(set(a) | set(b), key=lambda k: (k[0], float(k[1]))):
        if key not in a or key not in b:
            bad.append(f"eps={key[0]} N={key[1]}: only in {'coarse' if key in a else 'fine'}")
            continue
        ra, rb = a[key], b[key]
        same = len(ra) == len(rb) and all(fa == fb and abs(xa - xb) <= rtol * abs(xb) for (xa, fa), (xb, fb) in zip(ra, rb))
        if not same:
            bad.append(f"eps={key[0]} N={key[1]}: coarse {ra}, fine {rb}")
    return not bad, bad


def main():
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    opts = dict(x[2:].split("=", 1) if "=" in x else (x[2:], "1") for x in sys.argv[1:] if x.startswith("--"))
    rtol = float(opts.get("rtol", 2e-5))
    a, b = rows(args[0]), rows(args[1])
    if not a or not b:
        print(f"FAIL: no rows in {'COARSE' if not a else 'FINE'}")
        sys.exit(1)
    if "selftest" in opts:
        bent = dict(b)
        k = sorted(bent, key=lambda k: (k[0], float(k[1])))[0]
        bent[k] = bent[k][1:]
        ok, _ = compare(a, bent, rtol)
        print(f"selftest (lowest root of eps={k[0]} N={k[1]} dropped): {'FIRED' if not ok else 'DID NOT FIRE'}")
        if ok:
            sys.exit(1)
    ok, bad = compare(a, b, rtol)
    print(f"{'PASS' if ok else 'FAIL'}: {len(set(a) | set(b))} rows, {len(bad)} disagree, rtol {rtol:.0e}")
    for line in bad:
        print("  " + line)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

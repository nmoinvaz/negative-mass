"""
Compare #2's certified brackets of c r_* and c r_neg for every l = 2..3000 with #104's (#139).

  python compare_brackets.py --ours=PATH --theirs=PATH [--lmax=3000] [--mutate=NAME]

  ours    the csv of `threshold_isolation.py --csv=PATH` (ends as fractions p/q)
  theirs  `I104-rpi5/sweep_brackets.csv` of negative-mass-raw (ends as decimal strings)
Output: compare_brackets_output.txt, from this directory with negative-mass-raw cloned beside the checkout:
  python compare_brackets.py --ours=../../../negative-mass-raw/I139-rpi5/brackets_i2.csv \
      --theirs=../../../negative-mass-raw/I104-rpi5/sweep_brackets.csv > compare_brackets_output.txt
Both have the header l,c,rstar_lo,rstar_hi,rneg_lo,rneg_hi and are read with fractions.Fraction, so no float enters a check.
Each bracket is a closed interval holding its root, so two brackets of one root must intersect; touching counts.

Checks (each prints PASS or FAIL; the exit status is 1 if any fails, or if a file does not parse):
  C1  l = 2..lmax, each once, in each file
  C2  c = (l+2)(l-1) in every row of each file
  C3  lo < hi in both brackets of every row of each file
  C4  for every l the c r_* brackets of the two files intersect
  C5  for every l the c r_neg brackets of the two files intersect
Then it prints the largest bracket width of each file, the smallest overlap (min(hi) - max(lo), negative for a gap) and the
largest distance between the midpoints of the two files' brackets, with the column and the first l where each occurs.

Controls: --mutate=shift | drop | swap change ours after it is read (MUT_L gives the l), and the run must exit 1:
  shift  adds 1e-12 to the four ends of the row
  drop   removes the row
  swap   exchanges the c r_* and c r_neg columns of the row
"""
import argparse
import csv
import hashlib
import os
import sys
from collections import namedtuple
from fractions import Fraction as Fr

HEADER = ['l', 'c', 'rstar_lo', 'rstar_hi', 'rneg_lo', 'rneg_hi']
Row = namedtuple('Row', 'c rs_lo rs_hi rn_lo rn_hi')
COLS = ('c r_*', 'c r_neg')
MUT_L = {'shift': 1500, 'drop': 777, 'swap': 2000}
SHIFT = Fr(1, 10 ** 12)

failures = []


def check(name, ok, msg_pass, msg_fail, detail=()):
    print(("PASS  " if ok else "FAIL  ") + name + ": " + (msg_pass if ok else msg_fail), flush=True)
    if not ok:
        failures.append(name)
        for line in detail:
            print("      " + line, flush=True)


def brackets(row):
    """[(column name, lo, hi)] of a row."""
    return [(COLS[0], row.rs_lo, row.rs_hi), (COLS[1], row.rn_lo, row.rn_hi)]


def sci(x):
    """d.dddddde+xx of a Fraction, rounded half up in exact arithmetic."""
    if x == 0:
        return "0.000000e+00"
    sign = "-" if x < 0 else ""
    x = abs(x)
    e = len(str(x.numerator)) - len(str(x.denominator))
    while x >= Fr(10) ** (e + 1):
        e += 1
    while x < Fr(10) ** e:
        e -= 1
    m = x / Fr(10) ** e
    q = (m.numerator * 2 * 10 ** 6 + m.denominator) // (2 * m.denominator)
    if q == 10 ** 7:
        q, e = 10 ** 6, e + 1
    return f"{sign}{q // 10 ** 6}.{q % 10 ** 6:06d}e{e:+03d}"


def dec(x, nd=22):
    """decimal string of a Fraction to nd places, truncated toward zero (display only)."""
    sign = "-" if x < 0 else ""
    x = abs(x)
    ip = x.numerator // x.denominator
    fp = (x - ip) * 10 ** nd
    return f"{sign}{ip}." + str(fp.numerator // fp.denominator).zfill(nd)


def lst(items, kind="l", limit=10):
    """'l = a, b, c' (first `limit`, then the count), or 'none'."""
    items = list(items)
    if not items:
        return "none"
    s = f"{kind} = " + ", ".join(str(i) for i in items[:limit])
    return s + (f" ... ({len(items)} in all)" if len(items) > limit else "")


def read(path):
    """-> (dict l -> Row of Fractions, list of l that occur more than once). ValueError on a bad header, field or number."""
    rows, dups = {}, []
    with open(path, newline='') as fh:
        rd = csv.reader(fh)
        head = next(rd, None)
        if head != HEADER:
            raise ValueError(f"{os.path.basename(path)}: header {head} is not {HEADER}")
        for n, rec in enumerate(rd, 2):
            if len(rec) != len(HEADER):
                raise ValueError(f"{os.path.basename(path)}, line {n}: {len(rec)} fields")
            try:
                l = int(rec[0])
                row = Row(*[Fr(x) for x in rec[1:]])
            except ValueError as e:
                raise ValueError(f"{os.path.basename(path)}, line {n}: {e}") from None
            if l in rows:
                dups.append(l)
            else:
                rows[l] = row
    return rows, dups


def sha(path):
    with open(path, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ours', required=True)
    ap.add_argument('--theirs', required=True)
    ap.add_argument('--lmax', type=int, default=3000)
    ap.add_argument('--mutate', default='', choices=[''] + sorted(MUT_L))
    args = ap.parse_args()

    try:
        print(f"ours   {os.path.basename(args.ours)}  sha256 {sha(args.ours)}")
        print(f"theirs {os.path.basename(args.theirs)}  sha256 {sha(args.theirs)}")
        ours, dup_o = read(args.ours)
        theirs, dup_t = read(args.theirs)
    except (OSError, ValueError) as e:
        print(f"FAIL  parse: {e}")
        return 1

    if args.mutate:
        L = MUT_L[args.mutate]
        if L not in ours:
            print(f"mutation {args.mutate}: row l = {L} is not in ours")
            return 2
        r = ours[L]
        if args.mutate == 'shift':
            ours[L] = Row(r.c, r.rs_lo + SHIFT, r.rs_hi + SHIFT, r.rn_lo + SHIFT, r.rn_hi + SHIFT)
        elif args.mutate == 'swap':
            ours[L] = Row(r.c, r.rn_lo, r.rn_hi, r.rs_lo, r.rs_hi)
        else:
            del ours[L]
        print(f"mutation: {args.mutate} at l = {L}")

    files = [('ours', ours, dup_o), ('theirs', theirs, dup_t)]
    expected = set(range(2, args.lmax + 1))

    miss = {n: sorted(expected - set(rows)) for n, rows, _ in files}
    extra = {n: sorted(set(rows) - expected) for n, rows, _ in files}
    dups = {n: sorted(set(d)) for n, _, d in files}
    ok1 = all(not (miss[n] or extra[n] or dups[n]) for n, _, _ in files)
    check("C1", ok1,
          f"l = 2..{args.lmax} once each in both files (ours {len(ours)} rows, theirs {len(theirs)} rows)",
          f"ours has {len(ours)} rows, theirs {len(theirs)}; "
          + "; ".join(f"{k} in {n}: {lst(d[n])}" for n, _, _ in files
                      for k, d in (("missing", miss), ("extra", extra), ("duplicate", dups))))

    for n, rows, _ in files:
        bad = [l for l in sorted(rows) if rows[l].c != (l + 2) * (l - 1)]
        check("C2", not bad, f"c = (l+2)(l-1) in every row of {n}",
              f"c != (l+2)(l-1) in {n} at {len(bad)} of {len(rows)} rows: {lst(bad)}")
    for n, rows, _ in files:
        bad = [f"{l} ({name})" for l in sorted(rows) for name, lo, hi in brackets(rows[l]) if not lo < hi]
        check("C3", not bad, f"lo < hi in both brackets of every row of {n}",
              f"lo >= hi in {n} at {len(bad)} brackets: {lst(bad)}")

    common = sorted(set(ours) & set(theirs))
    overlap = {}
    for l in common:
        for i, ((name, a, b), (_, c, d)) in enumerate(zip(brackets(ours[l]), brackets(theirs[l]))):
            overlap[l, i] = min(b, d) - max(a, c)
    for i, cid in enumerate(("C4", "C5")):
        bad = [l for l in common if overlap[l, i] < 0]
        detail = []
        for l in bad[:5]:
            (_, a, b), (_, c, d) = brackets(ours[l])[i], brackets(theirs[l])[i]
            detail.append(f"l = {l}: ours [{dec(a)}, {dec(b)}]  theirs [{dec(c)}, {dec(d)}]  overlap {sci(overlap[l, i])}")
        check(cid, not bad, f"{COLS[i]} brackets intersect for every l ({len(common)} l)",
              f"{COLS[i]} brackets do not intersect: fails at {len(bad)} of {len(common)} l: {lst(bad)}", detail)

    def extreme(items, sign):
        """(value, column, l) of the largest (sign +1) or smallest (sign -1) value; the first of equals."""
        best = None
        for v, i, l in items:
            if best is None or sign * v > sign * best[0]:
                best = (v, i, l)
        return best

    def show(label, best):
        print(f"{label:34s}" + ("n/a" if best is None else f"{sci(best[0])}  ({COLS[best[1]]}, l = {best[2]})"))

    print()
    for n, rows, _ in files:
        show(f"largest bracket width, {n}:",
             extreme(((hi - lo, i, l) for l in sorted(rows) for i, (_, lo, hi) in enumerate(brackets(rows[l]))), +1))
    show("smallest overlap:", extreme(((overlap[l, i], i, l) for l in common for i in (0, 1)), -1))
    gaps = ((abs((a + b) / 2 - (c + d) / 2), i, l) for l in common
            for i, ((_, a, b), (_, c, d)) in enumerate(zip(brackets(ours[l]), brackets(theirs[l]))))
    show("largest midpoint gap:", extreme(gaps, +1))

    print(f"\n{'all checks pass' if not failures else str(len(failures)) + ' check(s) FAIL: ' + ' '.join(failures)}")
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())

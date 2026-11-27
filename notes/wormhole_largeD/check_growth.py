"""Checks growth.csv against w2_output.txt (#200), and recomputes the remainder of the large-D series from its rows.

    python3 check_growth.py [--mutate=NAME] [--old-rev=REV]       (exit 1 if a gated check fails)

Rows (G1-G6): the file's rows against the crate's quad-double table in w2_output.txt (D >= 20, rounded to 12 decimals), against
profiles.py's TS13 values (D <= 11, read from its source text), and against REV's csv for D <= 11 (default 613e866).
Remainder (R0-R2): W^2 a^2 minus the series through Omega_k, times s^(k+1), s = 2(D - 3), from the csv rows; Omega_0 ... Omega_3 are
the closed forms of #26 and #208, Omega_4 is #195's fit. An enclosure is the rounding of the 12th decimal given the crate's value
(5e-13 s^(k+1) for D >= 20); it is not the crate's own error and not a fit spread.

Shares with profiles.py: the input file w2_output.txt and TS13's values (the dict is the only record of its 15 digits). It has its
own parser and its own rounding (exact rationals, not Decimal). --mutate=NAME breaks the data as read, in memory:
  last_digit:D   last digit of row D's W2a2 + 1 unit of the 12th decimal      drop_row:D   row D removed
  source:D       row D's source replaced by modes_r.py                         ts13_swap    W2a2 of rows 4 and 5 swapped
  crate_scale    every crate value times 1 + 1e-12                              round_trunc  the reference truncates, not rounds
  list50         the 64-value list's D = 50 value + 2e-18                       wrong_om2    Omega_2's pi^2 ln 2 term flipped
  wrong_om3      Omega_3's 384 zeta(3) term flipped                             s_off        s = 2(D - 2)
The last line prints the gated checks that failed.
"""
import ast
import csv
import math
import os
import re
import subprocess
import sys
from fractions import Fraction

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
W2_OUT = os.path.join(HERE, "w2_output.txt")
O5_OUT = os.path.join(HERE, "order5_output.txt")
D_ORDER = [4, 5, 6, 7, 8, 9, 10, 11, 20, 30, 50, 100, 200, 500, 1000, 2000]
SRC_CRATE, SRC_TS13 = "w2_crate qd (#195)", "TS13 Table I"
R_DS = (200, 500, 1000, 2000)
HALF = Fraction(1, 2 * 10**12)          # half a unit of the 12th decimal
failed = []


def check(name, ok, what):
    print(f"{'PASS' if ok else 'FAIL'}    {name}: {what}")
    if not ok:
        failed.append(name)


# ---------------------------------------------------------------------------------------------------------------
# exact decimal helpers (rationals; no float, no Decimal)
# ---------------------------------------------------------------------------------------------------------------
def fix(q, nd):
    """the rational q rounded half to even to nd decimals, as a string of exactly nd decimals"""
    n = round(q * 10**nd)
    t = str(abs(n)).rjust(nd + 1, "0")
    return ("-" if n < 0 else "") + t[:-nd] + "." + t[-nd:]


def trunc(q, nd):
    n = math.floor(abs(q) * 10**nd) * (1 if q >= 0 else -1)
    t = str(abs(n)).rjust(nd + 1, "0")
    return ("-" if n < 0 else "") + t[:-nd] + "." + t[-nd:]


def ndec(s):
    return len(s.split(".")[1])


# ---------------------------------------------------------------------------------------------------------------
# readers
# ---------------------------------------------------------------------------------------------------------------
def read_csv_text(text):
    rd = list(csv.reader(text.splitlines()))
    return rd[0], rd[1:]


def read_crate_table(path):
    """{D: (old label, old csv string, crate string)} of the table "growth.csv against the crate (qd)" of w2_output.txt"""
    out, inside = {}, False
    pat = re.compile(r"^\s*(\d+)\s+(\S.*?)\s{2,}(\d+\.\d+)\s+(\d+\.\d{15,})\s+(-?\d\.\d+e[+-]\d+)\s")  # any source label
    for line in open(path):
        if "growth.csv against the crate (qd)" in line:
            inside = True
        elif inside and not line.strip() and out:
            break
        elif inside:
            m = pat.match(line)
            if m:
                out[int(m.group(1))] = (m.group(2), m.group(3), m.group(4))
    return out


def read_list(path):
    out = {}
    for line in open(path):
        m = re.match(r"^\s*D=(\d+)\s+(\d+\.\d+)\s+digits qd/od", line)
        if m:
            out[int(m.group(1))] = m.group(2)
    return out


def read_ts13_text(path):
    """{D: source text of the number} of the TS13 dict in profiles.py"""
    src = open(path).read()
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "TS13":
            return {k.value: ast.get_source_segment(src, v) for k, v in zip(node.value.keys, node.value.values)}
    raise SystemExit("TS13 not found in profiles.py")


def read_w2_remainder_table(path):
    """{D: printed rem*s^3 (crate qd column)} of w2_output.txt's four-term remainder table"""
    out, inside = {}, False
    for line in open(path):
        if "remainder x s^3" in line:
            inside = True
        elif inside:
            m = re.match(r"^\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", line)
            if m:
                out[int(m.group(1))] = m.group(2)
    return out


def read_o5_rem4(path):
    """{D: printed rem*s^4} of order5_output.txt, the 'reported, not gated' table (60-digit qd W^2)"""
    out, inside = {}, False
    for line in open(path):
        if "reported, not gated: rem*s^4" in line:
            inside = True
        elif inside:
            m = re.match(r"^\s+(\d+)\s+(-\d+\.\d+)\s+(-?\d+\.\d+)\s*$", line)
            if m:
                out.setdefault(int(m.group(1)), m.group(2))
    return out


def read_o5_omegas(path):
    """Omega0 ... Omega3 of the first block after the line 'B4:' of order5_output.txt"""
    out, inside = {}, False
    for line in open(path):
        if line.strip().startswith("B4:"):
            inside = True
        elif inside:
            m = re.match(r"^\s+Omega(\d) = (\S+)\s*$", line)
            if m:
                out[int(m.group(1))] = m.group(2)
            if len(out) == 4:
                break
    return out


# ---------------------------------------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------------------------------------
def main(argv):
    mut, old_rev = None, "613e866"
    for a in argv:
        if a.startswith("--mutate="):
            mut = a.split("=", 1)[1]
        elif a.startswith("--old-rev="):
            old_rev = a.split("=", 1)[1]
    mname, _, marg = (mut or "").partition(":")
    marg = int(marg) if marg else None

    header, rows = read_csv_text(open(os.path.join(HERE, "growth.csv")).read())
    csvrows = {}
    order = []
    for r in rows:
        order.append(int(r[0]))
        csvrows[int(r[0])] = r
    tab = read_crate_table(W2_OUT)
    crate = {D: v[2] for D, v in tab.items()}
    lst = read_list(W2_OUT)
    ts13 = read_ts13_text(os.path.join(HERE, "profiles.py"))
    rnd = trunc if mname == "round_trunc" else fix
    print(f"check_growth.py  mutate={mut}  old-rev={old_rev}")

    # ----- mutations of the data as read
    if mname == "last_digit":
        r = csvrows[marg]
        r[1] = fix(Fraction(r[1]) + Fraction(1, 10**12), 12)
    elif mname == "drop_row":
        order.remove(marg)
        del csvrows[marg]
    elif mname == "source":
        csvrows[marg][2] = "modes_r.py"
    elif mname == "ts13_swap":
        csvrows[4][1], csvrows[5][1] = csvrows[5][1], csvrows[4][1]
    elif mname == "crate_scale":
        crate = {D: fix(Fraction(v) * (1 + Fraction(1, 10**12)), 25) for D, v in crate.items()}
    elif mname == "list50":
        lst[50] = fix(Fraction(lst[50]) + Fraction(2, 10**18), 60)
    elif mname not in ("", "round_trunc", "wrong_om2", "wrong_om3", "s_off"):
        raise SystemExit(f"unknown mutation {mut}")

    # ----- rows
    check("G1", header == ["D", "W2a2", "source"] and order == D_ORDER and all(len(r) == 3 for r in csvrows.values()),
          f"header {header}, D = {order}")
    bad2 = []
    for D in D_ORDER:
        if D >= 20:
            r = csvrows.get(D)
            if r is None or r[1] != rnd(Fraction(crate[D]), 12):
                bad2.append(D)
    check("G2", not bad2, f"{8 - len(bad2)} of 8 rows D >= 20 equal the crate's qd value rounded to 12 decimals" +
          (f"; NOT MET at D = {bad2}" if bad2 else ""))
    bad3 = [D for D in D_ORDER if D in csvrows and csvrows[D][2] != (SRC_CRATE if D >= 20 else SRC_TS13)] \
        + [D for D in D_ORDER if D not in csvrows]
    check("G3", not bad3, f"sources: `{SRC_CRATE}` for D >= 20, `{SRC_TS13}` for D <= 11" +
          (f"; NOT MET at D = {bad3}" if bad3 else ""))
    bad4 = []
    for D in range(4, 12):
        r = csvrows.get(D)
        if r is None or r[1] != rnd(Fraction(ts13[D]), 12):
            bad4.append(D)
    check("G4", not bad4, f"{8 - len(bad4)} of 8 rows D <= 11 equal profiles.py's TS13 value at 12 decimals" +
          (f"; NOT MET at D = {bad4}" if bad4 else ""))
    nd = ndec(crate[50])
    check("G5", fix(Fraction(lst[50]), nd) == crate[50],
          f"the table's D = 50 value {crate[50]} and the list's {lst[50][:nd + 6]}... rounded to {nd} decimals "
          f"{fix(Fraction(lst[50]), nd)}")
    oldrows = {}
    try:
        old = subprocess.run(["git", "-C", HERE, "show", f"{old_rev}:notes/wormhole_largeD/growth.csv"], capture_output=True,
                             text=True, check=True).stdout
        _, orows = read_csv_text(old)
        oldrows = {int(r[0]): r for r in orows}
        bad6 = [D for D in range(4, 12) if csvrows.get(D) != oldrows[D]]
        check("G6", not bad6, f"{8 - len(bad6)} of 8 rows D <= 11 equal {old_rev}'s rows" +
              (f"; NOT MET at D = {bad6}" if bad6 else ""))
    except (subprocess.CalledProcessError, FileNotFoundError, KeyError) as exc:
        print(f"SKIP    G6: git show {old_rev} failed ({exc}); the old column below is absent")

    print("\n  reported, not gated: the csv against the crate (exact), and the distance of the crate's value from a rounding tie")
    print("  D       csv                crate (w2_output.txt)         crate - csv      tie margin (units of the 12th decimal)")
    for D in D_ORDER:
        if D not in csvrows:
            continue
        q = Fraction(crate[D])
        diff = q - Fraction(csvrows[D][1])
        frac = (q * 10**12) % 1
        marg_tie = abs(frac - Fraction(1, 2))
        print(f"  {D:<6d}  {csvrows[D][1]:<18s} {crate[D]:<28s} {float(diff):>+11.2e}      {float(marg_tie):.3f}")

    print("\n  reported, not gated: TS13's printed values (profiles.py's dict) against the crate and against the old and new csv rows")
    print("  D       TS13                 crate - TS13  (rel.)         old row - TS13  (rel.)       new row - TS13  (rel.)")
    for D in sorted(ts13):
        if D not in csvrows or D not in crate:
            continue
        t, q = Fraction(ts13[D]), Fraction(crate[D])
        o = Fraction(oldrows[D][1]) if D in oldrows else None
        n = Fraction(csvrows[D][1])
        ocol = f"{float(o - t):>+10.2e} ({float(abs(o - t) / t):.1e})" if o is not None else f"{'n/a':>24s}"
        print(f"  {D:<6d}  {ts13[D]:<20s} {float(q - t):>+10.2e} ({float(abs(q - t) / t):.1e})     {ocol}     "
              f"{float(n - t):>+10.2e} ({float(abs(n - t) / t):.1e})")

    # ----- remainder
    mp.mp.dps = 120
    L2, z3, pi = mp.log(2), mp.zeta(3), mp.pi
    om = {0: -4 * L2, 1: 2 * pi**2 / 3 + 8 * L2**2,
          2: 80 * z3 - (mp.mpf(8) / 3) * pi**2 * L2 * (-1 if mname == "wrong_om2" else 1) - (mp.mpf(32) / 3) * L2**3,
          3: 22 * pi**4 / 45 + (mp.mpf(16) / 3) * pi**2 * L2**2 + (mp.mpf(32) / 3) * L2**4 - 320 * z3 * L2
             - 384 * z3 * (-1 if mname == "wrong_om3" else 1)}
    txt = open(W2_OUT).read()
    om4 = mp.mpf(re.search(r"Omega4 = (\S+)", txt).group(1))
    pr = read_o5_omegas(O5_OUT)
    dev = {k: abs(om[k] - mp.mpf(pr[k])) for k in range(4)}
    check("R0", all(v < mp.mpf("1e-50") for v in dev.values()),
          "closed forms against order5_output.txt's B4 values: " + ", ".join(f"Omega{k} {mp.nstr(dev[k], 2)}" for k in range(4)) +
          " (bar 1e-50)")

    def sdef(D):
        return mp.mpf(2 * (D - (2 if mname == "s_off" else 3)))

    def column(D, W, k):
        s = sdef(D)
        series = s + sum(om[j] / s**j for j in range(k + 1))
        return (mp.mpf(W) - series) * s**(k + 1)

    def enc(D, k):
        s = sdef(D)
        e = mp.mpf("5e-13") if D >= 20 else abs(mp.mpf(crate[D]) - mp.mpf(csvrows[D][1]))
        return e * s**(k + 1)

    p4, p5 = read_w2_remainder_table(W2_OUT), read_o5_rem4(O5_OUT)
    miss = [D for D in R_DS if D not in p4 or D not in p5 or D not in csvrows]
    if miss:
        check("R1", False, f"printed values or rows missing at D = {miss}")
        check("R2", False, f"printed values or rows missing at D = {miss}")
    else:
        d1 = {D: abs(column(D, csvrows[D][1], 2) - mp.mpf(p4[D])) for D in R_DS}
        t1 = {D: mp.mpf("5e-13") * sdef(D)**3 + mp.mpf("5e-8") for D in R_DS}
        check("R1", all(d1[D] <= t1[D] for D in R_DS), "four-term remainder x s^3 from the csv rows against w2_output.txt's printed "
              "column: " + ", ".join(f"D = {D}: {mp.nstr(d1[D], 2)} (bar {mp.nstr(t1[D], 2)})" for D in R_DS))
        d2 = {D: abs(column(D, csvrows[D][1], 3) - mp.mpf(p5[D])) for D in R_DS}
        t2 = {D: mp.mpf("5e-13") * sdef(D)**4 + mp.mpf("5e-8") for D in R_DS}
        check("R2", all(d2[D] <= t2[D] for D in R_DS), "five-term remainder x s^4 from the csv rows against order5_output.txt's "
              "printed column: " + ", ".join(f"D = {D}: {mp.nstr(d2[D], 2)} (bar {mp.nstr(t2[D], 2)})" for D in R_DS))

    print("\n  reported, not gated. s = 2(D - 3); k-term remainder = W^2 a^2 - [s + Om0 + Om1/s + ...]; columns are remainder x s^2 "
          "(3 terms, through Om1; tends to Om2 + Om3/s),")
    print("  x s^3 (4 terms, through Om2; tends to Om3), x s^4 (5 terms, through Om3; tends to Om4). +- is the enclosure "
          "(5e-13 s^p for D >= 20; |crate - csv| s^p for D <= 11).")
    print(f"  Omega0 = {mp.nstr(om[0], 15)}, Omega1 = {mp.nstr(om[1], 15)}, Omega2 = {mp.nstr(om[2], 15)}, "
          f"Omega3 = {mp.nstr(om[3], 15)}, Omega4 = {mp.nstr(om4, 14)} (fit, #195)")
    print("\n  D      s     W2a2 (csv)           rel. diff    3 terms x s^2 (old csv)   (new csv) +-        4 terms x s^3  +-         "
          "5 terms x s^4  +-")
    for D in D_ORDER:
        if D not in csvrows:
            continue
        s = sdef(D)
        W = csvrows[D][1]
        f3 = s + om[0] + om[1] / s
        rel = abs(mp.mpf(W) - f3) / mp.mpf(W)
        c1, c2, c3 = column(D, W, 1), column(D, W, 2), column(D, W, 3)
        oldc = f"{float(column(D, oldrows[D][1], 1)):>12.4f}" if D in oldrows else f"{'n/a':>12s}"
        print(f"  {D:<6d}{int(s):<6d}{W:<20s} {float(rel):10.3e}   {oldc}            {float(c1):>12.4f} {float(enc(D, 1)):9.1e}   "
              f"{float(c2):>12.4f} {float(enc(D, 2)):9.1e}   {float(c3):>12.2f} {float(enc(D, 3)):9.1e}")
    print("\n  the four-term series' relative difference from W^2 a^2, and the old and new csv rows' W2a2 where they differ")
    for D in D_ORDER:
        if D in csvrows and D in oldrows and csvrows[D][1] != oldrows[D][1]:
            s = sdef(D)
            f4 = s + om[0] + om[1] / s + om[2] / s**2
            print(f"  {D:<6d} old {oldrows[D][1]:<18s} new {csvrows[D][1]:<18s} new - old {float(mp.mpf(csvrows[D][1]) - mp.mpf(oldrows[D][1])):>+10.2e}"
                  f"   rel. diff of the four-term series {float(abs(mp.mpf(csvrows[D][1]) - f4) / mp.mpf(csvrows[D][1])):.2e}")

    print("\nfailed:", " ".join(failed) if failed else "none")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

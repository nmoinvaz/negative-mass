#!/usr/bin/env python3
"""Torii and Shinkai's Table I (TS13) against growth.csv and #195's crate (#202).

    python3 ts13_digits.py                 table and checks; exit 1 if a check fails
    python3 ts13_digits.py --mutate=NAME   one input broken on purpose (names below); exit 1 means a check fired
    python3 ts13_digits.py --selftest      every mutation must fire exactly its declared checks; exit 1 if not
    python3 ts13_digits.py --tex=PATH      TeX source of arXiv:1309.2058v1 for T3 (default: the raw folder I202-rpi5 beside this
                                           checkout, ../../../negative-mass-raw/I202-rpi5/ts13_download/src/NdimWH_130910rev2.tex)

Source of the printed values: T. Torii, H. Shinkai, PRD 88, 064027 (2013), arXiv:1309.2058v1, Table I, "The negative
eigenvalues omega^2", n = D, throat radius a = 1 (their note 1: t, r, R scaled by a), perturbations e^{i omega t}. Read from the
rendered page, the PDF text layer and the TeX source (raw I202-rpi5/ts13_table_I.txt); the journal version was not read.
The n = 100 row is typeset "- 191.283017729717" with a space after the minus sign.

Conversion: a mode e^{W t} is e^{i omega t} with omega^2 = -W^2, so W^2 a^2 = -omega^2 (a = 1), and the printed sign is
negative for every row.

Checks (each prints PASS, FAIL or SKIP):
  T1  the typed table has the 11 rows D = 4..11, 20, 50, 100, each with 15 significant digits
  T2  it equals profiles.py's TS13 dict (which holds -omega^2), read from its source text (an earlier transcription), 11 of 11
  T3  it equals the TeX source's table as strings, 11 of 11 (SKIP if the TeX file is absent)
  R1  README.md's quoted [TS13] values are the printed ones rounded to the digits quoted
  C1  the old growth.csv's D = 4..11 rows (OLD_ROWS: main's file at 4bced01, before #221) equal -omega^2 rounded to 12 decimals,
      and 12 is below the decimals printed: those rows held no digit beyond the printed ones. Does not read the current csv.
  C2  the current growth.csv's D = 4..11 rows equal the crate's value rounded to 12 decimals (check_growth.py's G2 for these rows,
      with this script's own parser and Decimal rounding)
  K1  the conversion against the crate: relative difference below 1e-8 (the digits README.md quotes at D = 4 and 10) at all 11 rows
  P1  control: the comparator "agrees with the reference to the last printed decimal, half a unit" accepts the crate's own values
      rounded to 15 significant digits, 11 of 11
  O1  one_mode.log's N4 values (a code sharing none of the crate's, the equation of W1) equal the crate rounded to the digits
      that log prints, 4 of 4
  W1  TS13's regularized potential (3.18) equals one_mode.py's U exactly, at 4 radii for each of the 11 D
  W2  one_mode.py's source holds that U verbatim; W1 evaluates the text found there
Reported, not gated: the digits TS13's printed values share with the crate (F1), and the old rows (OLD_ROWS) and the current rows
of growth.csv against the crate (#221 took option B of #202: the crate rounded to 12 decimals).
"""
import argparse
import ast
import csv
import math
import os
import re
import sys
from decimal import ROUND_HALF_EVEN, Decimal, getcontext
from fractions import Fraction

getcontext().prec = 80
HERE = os.path.dirname(os.path.abspath(__file__))
TEX_DEFAULT = os.path.normpath(os.path.join(HERE, "../../../negative-mass-raw/I202-rpi5/ts13_download/src/NdimWH_130910rev2.tex"))

# TS13 Table I as printed: omega^2, n = D
PRINTED = {4: "-1.39705243371511", 5: "-2.98495893027790", 6: "-4.68662054299460", 7: "-6.46258414126318",
           8: "-8.28975936306259", 9: "-10.1535530451867", 10: "-12.0442650147438", 11: "-13.9552091676647",
           20: "-31.5751101285105", 50: "-91.3457759137153", 100: "-191.283017729717"}
D_CSV = range(4, 12)
# growth.csv's D = 4..11 rows on main before #221 (TS13's printed values rounded to 12 decimals, source "TS13 Table I"):
# git show 4bced01:notes/wormhole_largeD/growth.csv, lines 2-9; byte-identical at 613e866. Kept here for C1; the file now holds the crate's.
OLD_ROWS = {4: "1.397052433715", 5: "2.984958930278", 6: "4.686620542995", 7: "6.462584141263", 8: "8.289759363063",
            9: "10.153553045187", 10: "12.044265014744", 11: "13.955209167665"}
BAR_K1 = Decimal("1e-8")
U_TEXT = "(n - 2) * (n - 4) / (4 * r ** 2) - 3 * n ** 2 / (4 * r ** (2 * n))"

# name: (what is broken, the checks it must fire and no others; a SKIPped check is dropped from the set)
MUTATIONS = {
    "digit_last": ("TS13 D = 4: 14th decimal 1 -> 2", {"T2", "T3"}),
    "digit_8": ("TS13 D = 4: 8th decimal 3 -> 4", {"T2", "T3", "R1", "C1"}),
    "digit_12": ("TS13 D = 9: 12th decimal 6 -> 7", {"T2", "T3", "C1"}),
    "digit_6": ("TS13 D = 7: 6th decimal 4 -> 5", {"T2", "T3", "C1", "K1"}),
    "drop_zero": ("TS13 D = 5: the printed trailing 0 dropped", {"T1", "T3"}),
    "sign": ("conversion W^2 = +omega^2", {"C1", "K1"}),
    "sqrt": ("conversion W = sqrt(-omega^2)", {"C1", "K1"}),
    "csv": ("growth.csv D = 6: 12th decimal + 1", {"C2"}),
    "old_rows": ("OLD_ROWS D = 6: 12th decimal + 1", {"C1"}),
    "crate_lowD": ("crate D = 6: 12th decimal + 1 (2e-12 relative, under K1's bar)", {"C2"}),
    "profiles": ("profiles.py's dict D = 11: last digit + 1", {"T2"}),
    "crate": ("crate D = 20: 6th decimal + 1 (3e-8 relative, over K1's bar)", {"K1", "O1"}),
    "crate_7": ("crate D = 20: 7th decimal + 1 (3e-9 relative, under K1's bar)", {"O1"}),
    "onemode": ("one_mode.log D = 4: last digit + 1", {"O1"}),
    "wformula": ("W(r) with 3 (n-2)^2 -> 4 (n-2)^2", {"W1"}),
    "upot": ("one_mode.py's U text: n ** 2 -> n ** 3", {"W1", "W2"}),
    "halfunit": ("the comparator's tolerance 0.5 -> 0.1 unit", {"P1"}),
}


def dec(s):
    return Decimal(s.replace(" ", ""))


def rnd(x, nd):
    return x.quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_EVEN)


def rnd_sig(x, ns):
    return x.quantize(Decimal(1).scaleb(x.adjusted() - ns + 1), rounding=ROUND_HALF_EVEN)


def ndec(s):
    return len(s.split(".")[1])


def nsig(s):
    return len(s.lstrip("-").replace(".", "").lstrip("0"))


def bump(s, k):
    """the k-th decimal digit of s (k = 1 first) plus one, mod 10"""
    i = s.index(".") + k
    return s[:i] + str((int(s[i]) + 1) % 10) + s[i + 1:]


def read_csv():
    with open(os.path.join(HERE, "growth.csv")) as f:
        rows = list(csv.reader(f))[1:]
    return {int(d): (w, src) for d, w, src in rows}


def read_crate():
    """{D: crate (qd) value} from the table of w2_output.txt headed 'growth.csv against the crate'"""
    pat = re.compile(r"^\s*(\d+)\s+(\S.*?)\s{2,}(\d+\.\d+)\s+(\d+\.\d{15,})\s+(-?\d\.\d+e[+-]\d+)\s")
    return {int(m.group(1)): m.group(4) for m in map(pat.match, open(os.path.join(HERE, "w2_output.txt"))) if m}


def read_profiles():
    """{D: source text of the number} of the TS13 dict in profiles.py"""
    src = open(os.path.join(HERE, "profiles.py")).read()
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "TS13":
            return {k.value: ast.get_source_segment(src, v) for k, v in zip(node.value.keys, node.value.values)}
    raise SystemExit("TS13 not found in profiles.py")


def read_one_mode():
    pat = re.compile(r"^PASS\s+N4: D = (\d+): even bound state -E = (\S+) against W\^2 = \S+ \(relative \S+\)$")
    return {int(m.group(1)): m.group(2) for m in map(pat.match, open(os.path.join(HERE, "one_mode.log"))) if m}


def read_readme():
    out = {}
    for line in open(os.path.join(HERE, "README.md")):
        m = re.match(r"^\|\s*(\d+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|", line)
        if m and int(m.group(1)) in PRINTED:
            out[int(m.group(1))] = m.group(3)
    return out


def read_tex(path):
    if not os.path.exists(path):
        return None
    out = {}
    for line in open(path, encoding="latin-1"):  # the file holds non-UTF-8 bytes in a comment
        m = re.match(r"^(?:~~)?(\d+)&\$(?:~~)?(-\d+\.\d+)", line)
        if m:
            out[int(m.group(1))] = m.group(2)
    return out


def load(mut, tex_path):
    d = dict(printed=dict(PRINTED), old=dict(OLD_ROWS), csv=read_csv(), crate=read_crate(), prof=read_profiles(), one=read_one_mode(),
             readme=read_readme(), tex=read_tex(tex_path), wcoef=3, half=Decimal("0.5"), convert="neg")
    one_src = open(os.path.join(HERE, "one_mode.py")).read()
    if mut == "upot":
        one_src = one_src.replace("3 * n ** 2 / (4", "3 * n ** 3 / (4")
    m = re.search(r"U = Ufun or \(lambda r: (.+)\)\s*$", one_src, re.M)
    d["u_expr"] = m.group(1) if m else None
    if mut == "digit_last":
        d["printed"][4] = PRINTED[4][:-1] + "2"
    elif mut == "digit_8":
        d["printed"][4] = bump(PRINTED[4], 8)
    elif mut == "digit_12":
        d["printed"][9] = bump(PRINTED[9], 12)
    elif mut == "digit_6":
        d["printed"][7] = bump(PRINTED[7], 6)
    elif mut == "drop_zero":
        d["printed"][5] = PRINTED[5][:-1]
    elif mut in ("sign", "sqrt"):
        d["convert"] = mut
    elif mut == "csv":
        w, src = d["csv"][6]
        d["csv"][6] = (w[:-1] + str((int(w[-1]) + 1) % 10), src)
    elif mut == "old_rows":
        d["old"][6] = d["old"][6][:-1] + str((int(d["old"][6][-1]) + 1) % 10)
    elif mut == "crate_lowD":
        d["crate"][6] = bump(d["crate"][6], 12)
    elif mut == "profiles":
        d["prof"][11] = d["prof"][11][:-1] + str((int(d["prof"][11][-1]) + 1) % 10)
    elif mut in ("crate", "crate_7"):
        d["crate"][20] = bump(d["crate"][20], 6 if mut == "crate" else 7)
    elif mut == "onemode":
        d["one"][4] = d["one"][4][:-1] + str((int(d["one"][4][-1]) + 1) % 10)
    elif mut == "wformula":
        d["wcoef"] = 4
    elif mut == "halfunit":
        d["half"] = Decimal("0.1")
    return d


def w2_of(printed, convert="neg"):
    """W^2 a^2 from the printed omega^2 (a = 1)"""
    x = dec(printed)
    return {"neg": -x, "sign": x, "sqrt": (-x).sqrt()}[convert]


def agrees(x_str, ref, half):
    """x_str equals ref to the last printed decimal: |x - ref| <= half a unit of that decimal"""
    return abs(dec(x_str) - ref) <= half * Decimal(1).scaleb(-ndec(x_str))


def w_ts13(D, R, coef):
    """TS13 (3.18), n = D, a = 1: W(r) = -[coef (n-2)^2 R^(-2(n-3)) - (n-4)(n-6)] / (4 R^2), R = R_0(r)"""
    return -(coef * (D - 2) ** 2 * R ** (-2 * (D - 3)) - (D - 4) * (D - 6)) / (4 * R * R)


def u_one_mode(expr, D, R):
    """one_mode.py N4's U (its source text, evaluated), n = D - 2, r = R"""
    return eval(expr, {"__builtins__": {}}, {"n": D - 2, "r": R})


def run(mut, tex_path, verbose):
    d = load(mut, tex_path)
    pr, res = d["printed"], []

    def check(name, ok, msg):
        res.append((name, ok, msg))
        if verbose:
            print(f"{'PASS' if ok else 'FAIL'}    {name}: {msg}")

    def skip(name, msg):
        res.append((name, None, msg))
        if verbose:
            print(f"SKIP    {name}: {msg}")

    Ds = list(PRINTED)
    crate = {D: Decimal(d["crate"][D]) for D in Ds}
    w2 = {D: w2_of(pr[D], d["convert"]) for D in Ds}
    bad1 = [D for D in Ds if D not in pr or nsig(pr[D]) != 15 or not pr[D].startswith("-")]
    check("T1", not bad1 and sorted(pr) == sorted(PRINTED), f"{len(Ds) - len(bad1)} of {len(Ds)} rows hold 15 significant digits"
          + (f"; NOT MET at D = {bad1}" if bad1 else ""))
    bad2 = [D for D in Ds if -dec(pr[D]) != dec(d["prof"][D])]
    check("T2", not bad2, f"{len(Ds) - len(bad2)} of {len(Ds)} rows equal -(profiles.py's TS13 dict), which holds -omega^2" + (f"; differ at D = {bad2}" if bad2 else ""))
    if d["tex"] is None:
        skip("T3", f"{tex_path} not found")
    else:
        bad3 = [D for D in Ds if d["tex"].get(D) != pr[D]]
        check("T3", not bad3 and len(d["tex"]) == len(Ds), f"{len(Ds) - len(bad3)} of {len(Ds)} rows equal the TeX source's table"
              + (f"; differ at D = {bad3}" if bad3 else ""))
    if len(d["readme"]) < 5:
        skip("R1", f"README.md's [TS13] column not found ({sorted(d['readme'])})")
    else:
        bad_r = [D for D, q in d["readme"].items() if str(rnd(-dec(pr[D]), ndec(q))) != q]
        check("R1", not bad_r, f"{len(d['readme']) - len(bad_r)} of {len(d['readme'])} quoted values (D = {sorted(d['readme'])}) are the printed"
              f" ones rounded to their digits" + (f"; differ at D = {bad_r}" if bad_r else ""))
    bad_c = [D for D in D_CSV if d["old"][D] != str(rnd(w2[D], 12)) or ndec(pr[D]) <= 12]
    check("C1", not bad_c, f"{len(D_CSV) - len(bad_c)} of {len(D_CSV)} old rows D = 4..11 (growth.csv at 4bced01) equal -omega^2 at 12 decimals; "
          f"TS13 prints {min(ndec(pr[D]) for D in D_CSV)} to {max(ndec(pr[D]) for D in D_CSV)} decimals there" + (f"; NOT MET at D = {bad_c}" if bad_c else ""))
    bad_c2 = [D for D in D_CSV if d["csv"][D][0] != str(rnd(crate[D], 12))]
    check("C2", not bad_c2, f"{len(D_CSV) - len(bad_c2)} of {len(D_CSV)} rows D = 4..11 of growth.csv equal the crate's value at 12 decimals"
          + (f"; NOT MET at D = {bad_c2}" if bad_c2 else ""))
    rel = {D: abs(w2[D] - crate[D]) / crate[D] for D in Ds}
    bad_k = [D for D in Ds if rel[D] >= BAR_K1]
    check("K1", not bad_k, f"W^2 a^2 = -omega^2 against the crate: relative difference at most {max(rel.values()):.2e} (bar 1e-8); D = 4: "
          f"{w2[4]} against {str(crate[4])[:17]}, {rel[4]:.2e}" + (f"; NOT MET at D = {bad_k}" if bad_k else ""))
    ctl = {D: str(rnd_sig(crate[D], 15)) for D in Ds}
    bad_p = [D for D in Ds if not agrees(ctl[D], crate[D], d["half"])]
    check("P1", not bad_p, f"the comparator accepts the crate's own values at 15 significant digits: {len(Ds) - len(bad_p)} of {len(Ds)}"
          + (f"; NOT MET at D = {bad_p}" if bad_p else ""))
    bad_o = []
    for D, s in sorted(d["one"].items()):
        if s != str(rnd_sig(crate[D], nsig(s))):
            bad_o.append(D)
    check("O1", len(d["one"]) == 4 and not bad_o, f"{len(d['one']) - len(bad_o)} of {len(d['one'])} one_mode.log values (D = {sorted(d['one'])}) equal"
          f" the crate rounded to their {nsig(next(iter(d['one'].values())))} digits" + (f"; differ at D = {bad_o}" if bad_o else ""))
    radii = [Fraction(1), Fraction(3, 2), Fraction(2), Fraction(7, 3)]
    bad_w = [D for D in Ds if d["u_expr"] is None or any(w_ts13(D, R, d["wcoef"]) != u_one_mode(d["u_expr"], D, R) for R in radii)]
    check("W1", not bad_w, f"(3.18) equals one_mode.py's U at {len(radii)} radii for {len(Ds) - len(bad_w)} of {len(Ds)} D"
          + (f"; differs at D = {bad_w}" if bad_w else ""))
    check("W2", d["u_expr"] == U_TEXT, "one_mode.py's source line 'U = Ufun or (lambda r: ...)' holds U: " + U_TEXT)
    if verbose:
        report(d, w2, crate, rel)
    return res


def report(d, w2, crate, rel):
    pr, csvr = d["printed"], d["csv"]
    print("\n  TS13 Table I as printed (omega^2, a = 1), its conversion, and growth.csv (as in this checkout; D = 4..11 are the crate's since #221)")
    print("  D    printed omega^2        digits  W^2 a^2 = -omega^2     growth.csv value (source)                  csv - round(-omega^2, 12)")
    for D in PRINTED:
        w, src = csvr[D]
        print(f"  {D:<4d} {pr[D]:<22s} {nsig(pr[D]):>3d}    {str(w2[D]):<22s} {w + ' (' + src + ')':<42s} {float(Decimal(w) - rnd(w2[D], 12)):+.1e}")
    print("\n  against the crate (qd, w2_output.txt), and F1: the digits TS13's printed values share with it")
    print("  D    crate (qd)                       crate - TS13    relative    agreeing digits   last printed decimal")
    for D in PRINTED:
        diff = crate[D] - w2[D]
        print(f"  {D:<4d} {str(crate[D])[:30]:<32s} {float(diff):+.2e}    {float(rel[D]):.1e}     {-math.log10(rel[D]):6.2f}            "
              f"{ndec(pr[D])}th: {'agrees' if agrees(str(w2[D]), crate[D], Decimal('0.5')) else 'differs by ' + format(float(abs(diff) / Decimal(1).scaleb(-ndec(pr[D]))), '.3g') + ' units'}")
    print(f"  F1: TS13's printed digits agree with the crate to the last printed decimal at "
          f"{sum(agrees(str(w2[D]), crate[D], Decimal('0.5')) for D in PRINTED)} of {len(PRINTED)} rows; they agree to "
          f"{min(-math.log10(r) for r in rel.values()):.2f} .. {max(-math.log10(r) for r in rel.values()):.2f} significant digits")
    print("\n  one_mode.log N4 (TS13's own W(r), shooting; shares no code with the crate) against TS13 and the crate")
    print("  D    one_mode -E         crate rounded    TS13 - one_mode    crate - one_mode")
    for D, s in sorted(d["one"].items()):
        print(f"  {D:<4d} {s:<19s} {str(rnd_sig(crate[D], nsig(s))):<16s} {float(w2[D] - Decimal(s)):+.2e}         {float(crate[D] - Decimal(s)):+.2e}")
    print("\n  growth.csv's D = 4..11 rows: A the old rows (OLD_ROWS: TS13 printed, rounded), B the current rows (the crate rounded to 12 decimals)")
    print("  D    A - crate    B - crate    A = B at 8 decimals   A = B at 10 decimals")
    for D in D_CSV:
        a, b = Decimal(d["old"][D]), Decimal(csvr[D][0])
        print(f"  {D:<4d} {float(a - crate[D]):+.2e}   {float(b - crate[D]):+.2e}   {'yes' if rnd(a, 8) == rnd(b, 8) else 'NO ':<20s}  {'yes' if rnd(a, 10) == rnd(b, 10) else 'NO'}")
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mutate", choices=sorted(MUTATIONS))
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--tex", default=TEX_DEFAULT)
    a = ap.parse_args()
    if a.selftest:
        base = run(None, a.tex, False)
        fired0 = {n for n, ok, _ in base if ok is False}
        skipped = {n for n, ok, _ in base if ok is None}
        print(f"{'PASS' if not fired0 else 'FAIL'}    baseline: failing checks {sorted(fired0) or 'none'}; skipped {sorted(skipped) or 'none'}")
        bad = bool(fired0)
        for name, (what, target) in MUTATIONS.items():
            res = run(name, a.tex, False)
            fired = {n for n, ok, _ in res if ok is False}
            want = target - skipped
            ok = fired == want
            bad |= not ok
            print(f"{'PASS' if ok else 'FAIL'}    --mutate={name}: {what}; fires {sorted(fired)}, declared {sorted(want)}")
        print(f"{'selftest failed' if bad else 'selftest passed: every mutation fires its declared checks and no others'}")
        sys.exit(1 if bad else 0)
    print(f"ts13_digits.py  mutate={a.mutate}  tex={'found' if os.path.exists(a.tex) else 'absent'}\n")
    res = run(a.mutate, a.tex, True)
    failed = [n for n, ok, _ in res if ok is False]
    print(f"{len(failed)} failed: {failed}" if failed else "all checks pass")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()

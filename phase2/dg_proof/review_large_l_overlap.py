#!/usr/bin/env python3
"""Step 7 of the review (issue 104): this review's brackets against #2's stored output files.
Reads (only now) threshold_isolation_output.txt, large_l_proof_output.txt, monotone_proof_output.txt, which this review did not open before.
  (a) every row of #2's printed table (c r_* and c r_neg, 12 decimals, printed "lo" of brackets of width <= 2e-13) against my certified bracket of the same l
      (the sweep csv holds every l = 2..3000): the two intervals must overlap, allowing the 12-digit rounding of the printed value;
  (b) each of my 200 random brackets against #2's Theorem L (l >= 138: |c r - rho| <= K/c) and Theorem M (l >= 500: 0 < c r - rho <= s/c) as stated in
      their output files.
Usage: review_large_l_overlap.py --csv=PATH [--rho-star=...]   `--mutate=overlap_shift` shifts #2's printed values by 1e-9 (the overlap must then fail).
"""
import os
import re
import sys
from fractions import Fraction

from mpmath import mp, mpf

import review_large_l_common as C
import review_large_l_brackets as B
import review_large_l_roots as RT

KNOWN = ["overlap_shift", "bound_scale"]
HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    mut, rest = C.parse_args(sys.argv, KNOWN)
    path = next((a.split("=", 1)[1] for a in rest if a.startswith("--csv=")), None)
    if path is None:
        print("usage: review_large_l_overlap.py --csv=PATH")
        sys.exit(3)
    chk = C.Checker(failfast=mut is not None)
    rows = {}
    with open(path) as f:
        next(f)
        for ln in f:
            l, c, a, b, d, e = ln.strip().split(",")
            rows[int(l)] = (int(c), Fraction(a), Fraction(b), Fraction(d), Fraction(e))
    ti = open(os.path.join(HERE, "threshold_isolation_output.txt")).read()
    lp = open(os.path.join(HERE, "large_l_proof_output.txt")).read()
    mo = open(os.path.join(HERE, "monotone_proof_output.txt")).read()
    # (a) table of section E: lines "  l  c_r_*  c_r_neg  l^2r_*  l^2r_neg" with 12 decimals
    theirs = {}
    in_e = False
    for ln in ti.splitlines():
        if ln.startswith("== E:"):
            in_e = True
        elif ln.startswith("== F:"):
            in_e = False
        elif in_e:
            m = re.match(r"\s*(\d+)\s+(\d+\.\d{12})\s+(\d+\.\d{12})\s+\d+\.\d+\s+\d+\.\d+\s*$", ln)
            if m:
                theirs[int(m.group(1))] = (Fraction(m.group(2)), Fraction(m.group(3)))
    print(f"review_large_l_overlap, mutation={mut}: {len(theirs)} rows of #2's table E, l = {sorted(theirs)}")
    chk.check("parsed 16 rows of #2's table E", len(theirs) == 16)
    shift = Fraction(1, 10 ** 9) if mut == "overlap_shift" else 0
    width2 = Fraction(2, 10 ** 13)  # #2's widest bracket
    unit = Fraction(1, 10 ** 12)  # the 12th decimal of the printed value
    # first version of this script assumed round-to-nearest printing and failed on 15 of 32 values (step7_overlap_try1_wrong_rounding_guess.log).
    # The data show truncation: the printed value p is at or below my bracket's middle by less than one unit of the 12th decimal, for all 32 values.
    diffs = []
    for l in sorted(theirs):
        c, nlo, nhi, plo, phi = rows[l]
        for (lo, hi), p in (((nlo, nhi), theirs[l][0] + shift), ((plo, phi), theirs[l][1] + shift)):
            diffs.append((lo + hi) / 2 - p)
    # #2's root lies in [x, x + w] with p <= x < p + 1e-12 and w <= 2e-13, so (root - p) is in [0, 1.2e-12); rounding to nearest would give [-5e-13, 7e-13)
    chk.check("printing convention of #2's table: 0 <= (middle of my bracket) - (printed value) < 1.2e-12 for all 32 values (truncation to 12 decimals)",
              all(0 <= d < unit + width2 for d in diffs), f"range [{float(min(diffs)):.3e}, {float(max(diffs)):.3e}]")
    chk.check("the same 32 differences are NOT all within round-to-nearest's range [-5e-13, 7e-13)", not all(-5 * unit / 10 <= d < 5 * unit / 10 + width2 for d in diffs),
              "so the table is truncated, not rounded")
    for l in sorted(theirs):
        c, nlo, nhi, plo, phi = rows[l]
        for name, (lo, hi), p in (("c r_*", (nlo, nhi), theirs[l][0] + shift), ("c r_neg", (plo, phi), theirs[l][1] + shift)):
            # #2's bracket is [x, x + w], 0 < w <= 2e-13, with p <= x < p + 1e-12 (truncation); the roots in both brackets must coincide
            ok = lo <= p + unit + width2 and hi >= p
            chk.check(f"l={l}: my bracket of {name} overlaps #2's printed bracket", ok, f"mine [{float(lo):.14f}, {float(hi):.14f}], #2 printed lo {float(p):.12f}", quiet=True)
    chk.check("all 32 brackets (16 l, c r_* and c r_neg) overlap #2's printed brackets", not any(n.startswith("l=") and "overlaps" in n for n in chk.fails), "")
    in_mine = [l for l in sorted(theirs) if l in B.MAIN]
    print(f"  of #2's printed l, {in_mine} are among my 200 random l; the others come from my sweep of every l = 2..3000")
    # (b) my 200 against #2's stated theorems
    mK = re.search(r"THEOREM L: for every integer l >= (\d+):.*?<= ([\d.]+)/c\s+and.*?<= ([\d.]+)/c", lp)
    mM = re.search(r"THEOREM M holds for every l >= (\d+):.*?<= ([\d.]+)/c,.*?<= ([\d.]+)/c", mo)
    lmin_L, KN, KP = int(mK.group(1)), Fraction(mK.group(2)), Fraction(mK.group(3))
    lmin_M, sN, sP = int(mM.group(1)), Fraction(mM.group(2)), Fraction(mM.group(3))
    print(f"  #2's Theorem L: l >= {lmin_L}, K_* = {KN}, K_neg = {KP}; Theorem M: l >= {lmin_M}, s_N = {sN}, s_P = {sP}")
    if mut == "bound_scale":
        KN, KP, sN, sP = KN / 1000, KP / 1000, sN / 1000, sP / 1000
    mp.dps = 40
    r_star = mp.findroot(RT.N_inf_mp, mpf("8.4857"), tol=mpf(10) ** -35)
    r_neg = mp.findroot(RT.Psi_inf_mp, mpf("9.636"), tol=mpf(10) ** -35)
    rs = Fraction(str(mp.nstr(r_star, 35)))
    rn = Fraction(str(mp.nstr(r_neg, 35)))
    badL = badM = 0
    nL = nM = 0
    for l in B.MAIN:
        c, nlo, nhi, plo, phi = rows[l]
        if l >= lmin_L:
            nL += 1
            if not (max(abs(nlo - rs), abs(nhi - rs)) * c <= KN and max(abs(plo - rn), abs(phi - rn)) * c <= KP):
                badL += 1
        if l >= lmin_M:
            nM += 1
            if not (nlo > rs and nhi < rs + sN / c and plo > rn and phi < rn + sP / c):
                badM += 1
    chk.check(f"my 200 brackets satisfy #2's Theorem L bound for the {nL} sampled l >= {lmin_L}", badL == 0, f"{badL} violations")
    chk.check(f"my 200 brackets satisfy #2's Theorem M bound 0 < c r - rho <= s/c for the {nM} sampled l >= {lmin_M}", badM == 0, f"{badM} violations")
    chk.finish()


if __name__ == "__main__":
    main()

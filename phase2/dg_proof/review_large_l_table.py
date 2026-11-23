#!/usr/bin/env python3
"""Review of RESULTS.md 4.3's table (issue 104): the rows of the table against this review's certified brackets (the sweep csv).
Usage: review_large_l_table.py --csv=PATH   (PATH = I104-rpi5/sweep_brackets.csv in the raw repository, written by review_large_l_sweep.py)
Each quoted digit must agree with the bracket, rounded to the digits RESULTS prints.  `--mutate=table_digit` changes one quoted digit.
"""
import sys
from fractions import Fraction

import review_large_l_common as C

KNOWN = ["table_digit"]

# RESULTS.md 4.3: l, c r_*, c r_neg, l^2 r_*, l^2 r_neg
TABLE = [
    (2, "9.115122", "10.696286", "9.115122", "10.696286"),
    (3, "8.789041", "10.143519", "7.910137", "9.129167"),
    (5, "8.604503", "9.833910", "7.682592", "8.780277"),
    (10, "8.517779", "9.689316", "7.886832", "8.971589"),
    (20, "8.494079", "9.649909", "8.128306", "9.234363"),
    (80, "8.486244", "9.636891", "8.384063", "9.520855"),
    (100, "8.486049", "9.636568", "8.403693", "9.543046"),
    (500, "8.485715", "9.636013", "8.468845", "9.616856"),
    (1000, "8.485704", "9.635995", "8.477244", "9.626388"),
    (3000, "8.485701", "9.635990", "8.482876", "9.632781"),
]
LIMITS = ("8.485700967", "9.635989443")  # RESULTS table, row infinity (9 decimals)


def rounds_to(lo, hi, text):
    """True if the whole bracket [lo, hi] rounds to the decimal string `text` (digits as printed)."""
    d = len(text.split(".")[1])
    half = Fraction(1, 2 * 10 ** d)
    v = Fraction(text)
    return v - half <= lo and hi <= v + half


def main():
    mut, rest = C.parse_args(sys.argv, KNOWN)
    path = None
    for a in rest:
        if a.startswith("--csv="):
            path = a.split("=", 1)[1]
    if path is None:
        print("usage: review_large_l_table.py --csv=PATH")
        sys.exit(3)
    chk = C.Checker(failfast=mut is not None)
    rows = {}
    with open(path) as f:
        next(f)
        for ln in f:
            l, c, a, b, d, e = ln.strip().split(",")
            rows[int(l)] = (int(c), Fraction(a), Fraction(b), Fraction(d), Fraction(e))
    print(f"review_large_l_table, mutation={mut}, {len(rows)} brackets read from {path}")
    for l, s1, s2, s3, s4 in TABLE:
        if mut == "table_digit" and l == 20:
            s1 = "8.494078"
        c, slo, shi, nlo, nhi = rows[l]
        # l^2 r = l^2 rho / c is increasing in rho, so the bracket of l^2 r is [l^2 lo/c, l^2 hi/c]
        k = Fraction(l * l, c)
        ok1 = rounds_to(slo, shi, s1)
        ok2 = rounds_to(nlo, nhi, s2)
        ok3 = rounds_to(k * slo, k * shi, s3)
        ok4 = rounds_to(k * nlo, k * nhi, s4)
        chk.check(f"RESULTS table row l={l}: c r_*, c r_neg, l^2 r_*, l^2 r_neg", ok1 and ok2 and ok3 and ok4,
                  f"mine c r_* = {float(slo):.9f}, c r_neg = {float(nlo):.9f}, l^2 r_* = {float(k * slo):.7f}, l^2 r_neg = {float(k * nlo):.7f}; flags {ok1, ok2, ok3, ok4}")
    # the limits row: the brackets for l = 3000 differ from the limit by (3.5155, 5.8394)/c ~ 4e-7, so compare with the review's own roots instead
    import review_large_l_roots as RT
    from mpmath import mp, mpf
    r_star = mp.findroot(RT.N_inf_mp, mpf("8.4857"), tol=mpf(10) ** -50)
    r_neg = mp.findroot(RT.Psi_inf_mp, mpf("9.636"), tol=mpf(10) ** -50)
    chk.check("RESULTS table, row infinity: rho_* and rho_neg to 9 decimals",
              abs(r_star - mpf(LIMITS[0])) < mpf("5e-10") and abs(r_neg - mpf(LIMITS[1])) < mpf("5e-10"), f"{mp.nstr(r_star, 12)}, {mp.nstr(r_neg, 12)}")
    chk.finish()


if __name__ == "__main__":
    main()

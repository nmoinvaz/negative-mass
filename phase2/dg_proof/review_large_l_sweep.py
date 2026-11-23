#!/usr/bin/env python3
"""Review of RESULTS.md 4.3 (issue 104), supplementary: certified brackets of c r_* and c r_neg for EVERY l = 2..3000 (same code as the
200-bracket run: review_large_l_brackets.certify_l), then the consecutive-l claims of RESULTS 4.3:
  r_* < r_neg (the Psi bracket lies beyond the N bracket), c r_* and c r_neg decrease strictly with l and stay above their limits,
  l^2 r_* and l^2 r_neg increase strictly from l = 5 (the minimum is at l = 5).
Usage: review_large_l_sweep.py [--lmax=3000] [--csv=PATH]; the full table of brackets goes to the csv (raw output, not committed here).
`--mutate=NAME` breaks one ingredient (exit non-zero).
"""
import sys
import time
from fractions import Fraction

from mpmath import mp, mpf

import review_large_l_common as C
import review_large_l_brackets as B
import review_large_l_roots as RT

KNOWN = ["coef_q", "shift_bracket", "rho_limit", "l2r_reverse", "S_flip"]


def main():
    mut, rest = C.parse_args(sys.argv, KNOWN)
    lmax = 3000
    csv = None
    for a in rest:
        if a.startswith("--lmax="):
            lmax = int(a.split("=", 1)[1])
        elif a.startswith("--csv="):
            csv = a.split("=", 1)[1]
    chk = C.Checker(failfast=mut is not None)
    t0 = time.time()
    print(f"review_large_l_sweep, mutation={mut}, l = 2..{lmax}")
    if not C.run_planted(chk):
        print("PLANTED PROBLEM FAILED, stop before real data")
        sys.exit(2)
    r_star = mp.findroot(RT.N_inf_mp, mpf("8.4857"), tol=mpf(10) ** -50)
    r_neg = mp.findroot(RT.Psi_inf_mp, mpf("9.636"), tol=mpf(10) ** -50)
    if mut == "rho_limit":
        r_star = r_star + mpf("0.02")
    limhi = {w: Fraction(int(mp.floor(r * mpf(10) ** 25)) + 1, 10 ** 25) for w, r in (("N", r_star), ("P", r_neg))}
    print(f"rho_* = {mp.nstr(r_star, 22)}  rho_neg = {mp.nstr(r_neg, 22)}")

    def F(nn):
        return Fraction(nn, 10 ** B.DIG)

    counters = {}
    prev = None
    bad = {"dec": [], "above": [], "l2r": [], "sep": []}
    out = open(csv, "w") if csv else None
    if out:
        out.write("l,c,rstar_lo,rstar_hi,rneg_lo,rneg_hi\n")
    for l in range(2, lmax + 1):
        res = B.certify_l(l, r_star, r_neg, counters, chk, mut)
        if out:
            out.write(f"{l},{C.c_of(l)},{B.dec(res['N'][0])},{B.dec(res['N'][1])},{B.dec(res['P'][0])},{B.dec(res['P'][1])}\n")
            out.flush()
        for w in ("N", "P"):
            if not F(res[w][0]) > limhi[w]:
                bad["above"].append((l, w))
                chk.check(f"{w} bracket above its limit, l={l}", False)
        if prev is not None:
            pl, pr = prev
            for w in ("N", "P"):
                if not F(res[w][1]) < F(pr[w][0]):
                    bad["dec"].append((l, w))
                    chk.check(f"c r decreases ({w}), l={pl} -> {l}", False)
                if pl >= 5:
                    if mut == "l2r_reverse":
                        ok = l * l * F(res[w][1]) / C.c_of(l) < pl * pl * F(pr[w][0]) / C.c_of(pl)
                    else:
                        ok = l * l * F(res[w][0]) / C.c_of(l) > pl * pl * F(pr[w][1]) / C.c_of(pl)
                    if not ok:
                        bad["l2r"].append((l, w))
                        chk.check(f"l^2 r increases ({w}), l={pl} -> {l}", False)
        prev = (l, res)
        if l % 250 == 0 or l == lmax:
            print(f"  l = {l}: c r_* in [{B.dec(res['N'][0])}, {B.dec(res['N'][1])}], c r_neg in [{B.dec(res['P'][0])}, {B.dec(res['P'][1])}]  ({time.time() - t0:.0f} s)", flush=True)
    if out:
        out.close()
    for name, (p, f) in counters.items():
        chk.check(f"{name}: all {p + f} l", f == 0, f"{p} pass, {f} fail")
    chk.check(f"c r_* and c r_neg strictly decrease for every consecutive l = 2..{lmax} (hi(l+1) < lo(l))", not bad["dec"], f"violations {bad['dec'][:3]}")
    chk.check(f"c r_* and c r_neg stay above their limits for every l = 2..{lmax}", not bad["above"], f"violations {bad['above'][:3]}")
    chk.check(f"l^2 r_* and l^2 r_neg increase strictly for every consecutive l = 5..{lmax} (lo(l+1)(l+1)^2/c(l+1) > hi(l) l^2/c(l))", not bad["l2r"], f"violations {bad['l2r'][:3]}")
    print(f"total time {time.time() - t0:.0f} s")
    chk.finish()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Review of RESULTS.md 4.3 (issue 104): 200 independently certified brackets of c r_* and c r_neg, l drawn at random from 2..3000.

Method (own): for each l all l+2 terms of X, X', X'' are summed in mpmath iv (no truncation, no tail bound), eps an enclosure of 1/c.
The root of N (resp. Psi) is guessed in mp, then a bracket of width 2e-14 with 20-digit decimal endpoints is certified:
N (Psi) enclosures exclude 0 with opposite signs at the ends, and N' (Psi') over the whole bracket excludes 0, so each bracket holds exactly one root.
That the N root is r_* c rests on T2, and that the Psi root is r_neg c on Theorem 3 (both assumed here; Theorem 3's steps are in
review_large_l_theorem3.py).  S (RESULTS 4.3) must show the same sign change on the Psi bracket (a certificate that does not use (Z)).
Seed 104: sorted(random.Random(104).sample(range(2, 3001), 200)).
`--mutate=NAME` breaks one ingredient (exit non-zero).
"""
import random
import sys
import time
from fractions import Fraction

import sympy as sp
from mpmath import iv, mp, mpf

import review_large_l_common as C
import review_large_l_roots as RT

KNOWN = ["coef_q", "coef_lin", "nz_coef", "drop_last", "S_flip", "accept_zero_enclosure", "no_deriv",
         "shift_bracket", "rho_limit", "l2r_reverse"]

MAIN = [29, 34, 42, 47, 48, 70, 78, 82, 95, 104, 108, 116, 122, 131, 134, 155, 158, 162, 165, 168, 169, 225, 231, 251, 256, 294, 295, 296,
        297, 298, 303, 314, 342, 348, 353, 358, 361, 384, 389, 412, 421, 431, 485, 507, 509, 513, 514, 531, 543, 548, 596, 604, 610,
        632, 644, 658, 668, 670, 684, 687, 704, 720, 730, 736, 738, 749, 764, 772, 791, 806, 812, 817, 821, 833, 866, 867, 870, 890,
        900, 926, 934, 940, 945, 975, 981, 1000, 1001, 1023, 1035, 1041, 1054, 1060, 1075, 1108, 1131, 1191, 1227, 1257, 1272, 1322,
        1352, 1356, 1358, 1402, 1430, 1435, 1459, 1512, 1552, 1560, 1584, 1598, 1603, 1604, 1616, 1621, 1679, 1700, 1712, 1713, 1751,
        1769, 1811, 1813, 1826, 1854, 1859, 1882, 1890, 1895, 1912, 1934, 1957, 1961, 1966, 1986, 1993, 2000, 2026, 2064, 2075, 2085,
        2102, 2128, 2141, 2142, 2164, 2176, 2191, 2198, 2205, 2219, 2220, 2244, 2254, 2299, 2301, 2302, 2306, 2315, 2330, 2370, 2371,
        2373, 2383, 2419, 2423, 2450, 2509, 2512, 2513, 2545, 2547, 2554, 2561, 2563, 2575, 2584, 2588, 2593, 2609, 2645, 2674, 2684,
        2686, 2736, 2819, 2826, 2844, 2850, 2865, 2892, 2906, 2925, 2926, 2941, 2949, 2959, 2962, 2999]
SUPP = sorted(set(list(range(2, 31)) + [138, 139, 140, 499, 500, 501, 1000, 2000, 2999, 3000]))

DIG = 20  # decimal digits of the bracket ends
HALF = 10 ** 6  # half width in units of 1e-20, so the width is 2e-14


def dec(n, d=DIG):
    return f"{n // 10 ** d}.{n % 10 ** d:0{d}d}"


class Evaluator:
    """Jets of X at points or intervals for one l, memoized."""

    def __init__(self, l):
        self.l = l
        self.c = C.c_of(l)
        self.eps = iv.mpf(1) / self.c
        self.a = C.coeffs(l)
        self.cache = {}

    def q(self, x):
        key = x._mpi_
        if key not in self.cache:
            self.cache[key] = C.quantities(C.jet_from_coeffs(self.a, x, 2), x, self.eps)
        return self.cache[key]


def guess_root(l, which, x0):
    """Root of N (which='N') or Psi ('P') in mp, starting at x0 (secant)."""
    c = C.c_of(l)
    eps = mpf(1) / c
    a = C.coeffs(l, K=mp.mpf)

    def F(x):
        q = C.quantities(C.jet_from_coeffs(a, x, 2), x, eps)
        return q["N"] if which == "N" else q["Psi"]

    return mp.findroot(F, (x0, x0 + mpf("0.0005")), solver="secant", tol=mpf(10) ** -35, maxsteps=60)


def certify_l(l, lim_star, lim_neg, counters, chk, mut):
    """Returns dict with the two brackets as integers (units 1e-20), or None on failure."""
    ev = Evaluator(l)
    c = ev.c
    # coefficients against Lemma A (exact Fractions) for n <= 40
    nn = min(l + 2, 40)
    ok = all(C.frac_in_iv(C.a_exact(l, n), ev.a[n]) for n in range(1, nn + 1))
    rec(counters, chk, "coefficients a_n(eps), n <= 40, inside Lemma A's exact closed form", ok, l)
    out = {}
    for which, x0 in (("N", lim_star + mpf("3.5") / c), ("P", lim_neg + mpf("5.84") / c)):
        root = guess_root(l, which, x0)
        n0 = int(mp.floor(root * mpf(10) ** DIG))
        shift = 10 ** 8 if mut == "shift_bracket" else 0  # 1e-12, far outside a width 2e-14 bracket
        nlo, nhi = n0 - HALF + shift, n0 + 1 + HALF + shift
        lo = C.endpoint(nlo, 10 ** DIG)
        hi = C.endpoint(nhi, 10 ** DIG)
        if which == "N":
            F = lambda x: ev.q(x)["N"]
            dF = lambda x: ev.q(x)["Np"]
            orient = (1, -1)
        else:
            F = lambda x: ev.q(x)["Psi"]
            dF = lambda x: ev.q(x)["Psip"]
            orient = (-1, 1)
        ok, why = C.certify_bracket(F, dF, lo, hi)
        name = "N: bracket of c r_* certified (signs, derivative)" if which == "N" else "Psi: bracket of c r_neg certified (signs, derivative)"
        rec(counters, chk, name, ok, l, why)
        s1, s2 = C.sgn(F(lo)), C.sgn(F(hi))
        rec(counters, chk, f"{'N' if which == 'N' else 'Psi'}: signs (lo, hi) = {orient}", (s1, s2) == orient, l, f"signs {(s1, s2)}")
        if which == "P":
            sl, sh = C.sgn(ev.q(lo)["S"]), C.sgn(ev.q(hi)["S"])
            rec(counters, chk, "S shows the same sign change (-,+) on the Psi bracket (does not use (Z))", (sl, sh) == (-1, 1), l, f"S signs {(sl, sh)}")
            xs = C.sgn(ev.q(lo)["X"]) + C.sgn(ev.q(hi)["X"])
            rec(counters, chk, "X > 0 at both ends of the Psi bracket (so sign V = sign Psi = sign S)", xs == 2, l)
        out[which] = (nlo, nhi)
    rec(counters, chk, "lo(c r_neg) > hi(c r_*) (the Psi bracket lies beyond r_*)", out["P"][0] > out["N"][1], l)
    return out


def rec(counters, chk, name, ok, l, detail=""):
    p, f = counters.setdefault(name, [0, 0])
    counters[name] = [p + (1 if ok else 0), f + (0 if ok else 1)]
    chk.check(f"{name}, l={l}", ok, detail, quiet=True)


def main():
    mut, _ = C.parse_args(sys.argv, KNOWN)
    chk = C.Checker(failfast=mut is not None)
    t0 = time.time()
    print(f"review_large_l_brackets, mutation={mut}")
    if not C.run_planted(chk):
        print("PLANTED PROBLEM FAILED, stop before real data")
        sys.exit(2)
    # the limits, by own root finding (review_large_l_roots)
    r_star = mp.findroot(RT.N_inf_mp, mpf("8.4857"), tol=mpf(10) ** -50)
    r_neg = mp.findroot(RT.Psi_inf_mp, mpf("9.636"), tol=mpf(10) ** -50)
    if mut == "rho_limit":
        r_star = r_star + mpf("0.02")  # a limit above the true one: "stay above the limit" must fail
    print(f"rho_* = {mp.nstr(r_star, 22)}  rho_neg = {mp.nstr(r_neg, 22)}")
    drawn = sorted(random.Random(104).sample(range(2, 3001), 200))
    chk.check("the list of 200 l equals sorted(random.Random(104).sample(range(2, 3001), 200)) and the list in the raw README", drawn == MAIN and len(set(MAIN)) == 200)
    results = {}
    counters = {}
    todo = list(MAIN) + [l for l in SUPP if l not in MAIN]
    print(f"{'l':>5} {'c':>9}  c r_* in [lo, hi] (20 digits)                         c r_neg in [lo, hi]   (width 2e-14)   seconds")
    for l in todo:
        t1 = time.time()
        res = certify_l(l, r_star, r_neg, counters, chk, mut)
        results[l] = res
        tag = "" if l in MAIN else " (supplementary)"
        print(f"{l:>5} {C.c_of(l):>9}  [{dec(res['N'][0])}, {dec(res['N'][1])}]  [{dec(res['P'][0])}, {dec(res['P'][1])}]  {time.time() - t1:5.1f}{tag}", flush=True)
    for name, (p, f) in counters.items():
        chk.check(f"{name}: all {p + f} l", f == 0, f"{p} pass, {f} fail")

    # ordering checks over the sampled l, as exact Fractions of the decimal endpoints
    def F(nn):
        return Fraction(nn, 10 ** DIG)

    for label, which in (("c r_*", "N"), ("c r_neg", "P")):
        order = sorted(results)
        main_order = [l for l in order if l in MAIN]
        bad = [(a, b) for a, b in zip(main_order, main_order[1:]) if not F(results[b][which][1]) < F(results[a][which][0])]
        chk.check(f"{label} strictly decreases along the 200 sampled l (hi(l2) < lo(l1))", not bad, f"violations: {bad[:3]}")
        lim = r_star if which == "N" else r_neg
        limhi = Fraction(int(mp.floor(lim * mpf(10) ** 25)) + 1, 10 ** 25)
        above = all(F(results[l][which][0]) > limhi for l in results)
        chk.check(f"{label} above its limit for every certified l (lo > rho)", above)
        # l^2 r = l^2 rho/c increases for l >= 5 (sampled pairs, consecutive among the sorted list)
        seq = [l for l in order if l >= 5]
        pairs = list(zip(seq, seq[1:]))
        if mut == "l2r_reverse":
            bad2 = [(a, b) for a, b in pairs if not (b * b * F(results[b][which][1]) / C.c_of(b) < a * a * F(results[a][which][0]) / C.c_of(a))]
        else:
            bad2 = [(a, b) for a, b in pairs if not (b * b * F(results[b][which][0]) / C.c_of(b) > a * a * F(results[a][which][1]) / C.c_of(a))]
        chk.check(f"l^2 r with {label}: increases along all sampled l >= 5 (lo(l2) l2^2/c2 > hi(l1) l1^2/c1)", not bad2, f"violations: {bad2[:3]}")
    # consecutive l (supplementary 2..30): strict decrease of c r and increase of l^2 r from l = 5
    cons = list(range(2, 31))
    for label, which in (("c r_*", "N"), ("c r_neg", "P")):
        bad = [l for l in cons[:-1] if not F(results[l + 1][which][1]) < F(results[l][which][0])]
        chk.check(f"{label} strictly decreasing for consecutive l = 2..30", not bad, f"violations at {bad}")
        bad2 = [l for l in cons[:-1] if l >= 5 and not ((l + 1) ** 2 * F(results[l + 1][which][0]) / C.c_of(l + 1) > l ** 2 * F(results[l][which][1]) / C.c_of(l))]
        chk.check(f"l^2 r with {label}: increasing for consecutive l = 5..30", not bad2, f"violations at {bad2}")
    # l = 2..4: l^2 r is non-monotone below l = 5 (RESULTS: minimum at l = 5)
    chk.check("l^2 r_* and l^2 r_neg have their minimum over l = 2..30 at l = 5 (RESULTS 4.3)",
              all(min((F(results[l][w][0]) * l * l / C.c_of(l), l) for l in range(2, 31))[1] == 5 for w in ("N", "P")))
    # sympy exact root isolation for l = 2..12 (the polynomials built from Lemma A's exact coefficients)
    rho = sp.Symbol("rho")
    for l in range(2, 13):
        c = C.c_of(l)
        e_ = sp.Rational(1, c)
        Xp = sum(sp.Rational(C.a_exact(l, n).numerator, C.a_exact(l, n).denominator) * rho ** n for n in range(1, l + 3))
        Np = sp.Poly(sp.expand(Xp - (rho - 6) * sp.diff(Xp, rho)), rho)
        (a_, b_), _m = Np.intervals(inf=6, eps=sp.Rational(1, 10 ** 18))[0]
        nlo, nhi = results[l]["N"]
        ok1 = Fraction(int(a_.p), int(a_.q)) <= F(nhi) and F(nlo) <= Fraction(int(b_.p), int(b_.q))
        chk.check(f"l={l}: sympy's isolating interval of the N root overlaps the bracket of c r_*", ok1)
        nz = (1 + 2 * e_) * rho ** 3 / 8 - sp.Rational(3, 4) * rho ** 2 + sp.Rational(9, 2) * rho - 9
        Pp = sp.Poly(sp.expand(rho ** 2 * (2 + e_ * rho) * Np.as_expr() ** 2 - 4 * nz * Xp ** 2), rho)
        lo_in = sp.Rational(results[l]["N"][1], 10 ** DIG)
        (a2, b2), _m2 = Pp.intervals(inf=lo_in, eps=sp.Rational(1, 10 ** 18))[0]
        nlo2, nhi2 = results[l]["P"]
        ok2 = Fraction(int(a2.p), int(a2.q)) <= F(nhi2) and F(nlo2) <= Fraction(int(b2.p), int(b2.q))
        chk.check(f"l={l}: sympy's isolating interval of the Psi root beyond r_* overlaps the bracket of c r_neg", ok2)
    # constants K(l) = (c r - rho_inf) c at the sampled l, for the record (numerical evidence of RESULTS 4.3's 3.5155 and 5.8394)
    print("  (c r_* - rho_*) c and (c r_neg - rho_neg) c at the bracket centers, l >= 1000:")
    for l in [x for x in sorted(results) if x in (1000, 2000, 2999, 3000)]:
        cs = F(sum(results[l]["N"])) / 2
        cn = F(sum(results[l]["P"])) / 2
        ks = (mpf(cs.numerator) / mpf(cs.denominator) - r_star) * C.c_of(l)
        kn = (mpf(cn.numerator) / mpf(cn.denominator) - r_neg) * C.c_of(l)
        print(f"    l={l}: {mp.nstr(ks, 8)}  {mp.nstr(kn, 8)}")
    print(f"total time {time.time() - t0:.0f} s")
    chk.finish()


if __name__ == "__main__":
    main()

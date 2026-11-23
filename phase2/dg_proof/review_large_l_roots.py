#!/usr/bin/env python3
"""Review of RESULTS.md 4.1/4.4 (issue 104): rho_* and rho_neg, the roots of N_inf and Psi_inf, to 25 digits.

Order: planted problem (exit 2 if any verdict is wrong), then the Bessel closed form of X_inf against the series,
then root finding in mp at 60 digits, then iv enclosures of width 2e-30 with a certified sign change and a derivative
enclosure that excludes 0, then comparison with the 16 digits quoted in RESULTS.md.
`--mutate=NAME` breaks one ingredient (exit non-zero).
"""
import sys

import sympy as sp
from mpmath import besseli, iv, mp, mpf, sqrt

import review_large_l_common as C

KNOWN = ["alpha_n", "nz_coef", "shift_enclosure", "no_deriv", "accept_zero_enclosure"]

iv.prec = 200
mp.dps = 60

QUOTED_STAR = mpf("8.485700967373530")
QUOTED_NEG = mpf("9.635989442586144")


def bessel_X(rho):
    s = rho / 2
    return s * ((s + 3) * besseli(0, 2 * sqrt(s)) - 4 * sqrt(s) * besseli(1, 2 * sqrt(s)))


def N_inf_mp(rho):
    X = C.jet_inf_mp(rho, 1, nmax=150)
    return X[0] - (rho - 6) * X[1]


def Psi_inf_mp(rho):
    X = C.jet_inf_mp(rho, 2, nmax=150)
    return C.quantities(X, rho, mpf(0))["Psi"]


def decimal_endpoints(root, digits, halfwidth_units):
    """Integers n_lo < n_hi with n/10^digits enclosing root, width 2*halfwidth_units units."""
    n = int(mp.floor(root * mpf(10) ** digits))
    return n - halfwidth_units, n + 1 + halfwidth_units


def dec(n, d):
    return f"{n // 10 ** d}.{n % 10 ** d:0{d}d}"


def main():
    mut, _ = C.parse_args(sys.argv, KNOWN)
    chk = C.Checker(failfast=mut is not None)
    print(f"review_large_l_roots, mutation={mut}")
    # 1. planted problem first
    if not C.run_planted(chk):
        print("PLANTED PROBLEM FAILED, stop before real data")
        sys.exit(2)
    # 2. X_inf: sympy termwise identity of the Bessel closed form, then numerical agreement with the series
    n_, s_ = sp.symbols("n s", positive=True)
    # coefficient of s^n in s[(s+3) I0(2 sqrt s) - 4 sqrt(s) I1(2 sqrt s)], I0(2 sqrt s) = sum s^k/k!^2, sqrt(s) I1(2 sqrt s) = sum s^(k+1)/(k!(k+1)!)
    coef = (1 / sp.factorial(n_ - 2) ** 2 + 3 / sp.factorial(n_ - 1) ** 2 - 4 / (sp.factorial(n_ - 2) * sp.factorial(n_ - 1)))
    target = (n_ - 2) * (n_ - 4) / sp.factorial(n_ - 1) ** 2
    okc = all(sp.simplify((coef - target).subs(n_, k)) == 0 for k in range(3, 14))
    chk.check("Bessel closed form: coefficient of s^n equals (n-2)(n-4)/(n-1)!^2 for n = 3..13", okc)
    for rho in (mpf("3.3"), mpf(6), mpf("8.4857"), mpf("9.636"), mpf(11)):
        X = C.jet_inf_mp(rho, 0, nmax=150)[0]
        d = abs(X - bessel_X(rho))
        chk.check(f"series of X_inf = Bessel closed form at rho={rho}", d < mpf(10) ** -40, f"diff={mp.nstr(d, 3)}")
    # 3. root finding at 60 digits
    r_star = mp.findroot(N_inf_mp, mpf("8.4857"), tol=mpf(10) ** -55)
    r_neg = mp.findroot(Psi_inf_mp, mpf("9.636"), tol=mpf(10) ** -55)
    print(f"rho_*   = {mp.nstr(r_star, 28)}")
    print(f"rho_neg = {mp.nstr(r_neg, 28)}")
    chk.check("residual of N_inf at rho_* below 1e-50", abs(N_inf_mp(r_star)) < mpf(10) ** -50)
    chk.check("residual of Psi_inf at rho_neg below 1e-45", abs(Psi_inf_mp(r_neg)) < mpf(10) ** -45)
    # 4. iv enclosures of width 2e-30 with certified sign change and derivative enclosure
    D = 10 ** 32
    for name, root, F in (("rho_*", r_star, "N"), ("rho_neg", r_neg, "Psi")):
        nlo, nhi = decimal_endpoints(root, 32, 100)
        if mut == "shift_enclosure":
            nlo, nhi = nlo + 10 ** 4, nhi + 10 ** 4  # the root is now outside the enclosure
        lo = C.endpoint(nlo, D)
        hi = C.endpoint(nhi, D)

        def Fpt(x, F=F):
            q = C.quantities(C.jet_inf(x, 2), x, iv.mpf(0))
            return q["N"] if F == "N" else q["Psi"]

        def dF(x, F=F):
            q = C.quantities(C.jet_inf(x, 2), x, iv.mpf(0))
            return q["Np"] if F == "N" else q["Psip"]

        ok, why = C.certify_bracket(Fpt, dF, lo, hi)
        chk.check(f"{name}: iv enclosure [{nlo}e-32, {nhi}e-32] certified (sign change, derivative excludes 0)", ok, why)
        print(f"  {name} in [{dec(nlo, 32)}, {dec(nhi, 32)}]  (decimal endpoints, width {nhi - nlo}e-32)")
        print(f"  {name}: F' enclosure over it = {dF(C.hull(lo, hi))}")
    # 5. comparison with RESULTS.md (16 digits)
    d1 = abs(r_star - QUOTED_STAR)
    d2 = abs(r_neg - QUOTED_NEG)
    chk.check("rho_* reproduces RESULTS.md's 8.485700967373530 to 1e-15", d1 < mpf(10) ** -15, f"diff={mp.nstr(d1, 3)}")
    chk.check("rho_neg reproduces RESULTS.md's 9.635989442586144 to 1e-15", d2 < mpf(10) ** -15, f"diff={mp.nstr(d2, 3)}")
    chk.finish()


if __name__ == "__main__":
    main()

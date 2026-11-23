#!/usr/bin/env python3
"""Review of RESULTS.md 4.2, Theorem 3 (issue 104): the identities and positivity steps, derived in sympy from the definitions.

Definitions used (RESULTS.md section 1 and 4): c = (l+2)(l-1), f = 1 + 2/r, d/dx = f d/dr, chi from Lemma A's closed form,
psi_0 = chi/(c r - 6), g = psi_0'/psi_0, Vhat = g^2 - g', V = 2 f [lam^2 r^2 ((lam+1) r - 3) + 9 (lam r - 1)]/(r^3 (lam r - 3)^2), lam = c/2.
Checked here: (Z), V in (rho, c), Vhat = 2g^2 - V, the Psi identity and the S identity, (P1), (P2), N_Z > 0, the sign of Vhat at infinity,
exact root counts for l = 2..16 (N beyond rho = 6, Psi beyond rho_*).
`--mutate=NAME` breaks one ingredient of one check (exit non-zero).
"""
import sys
from math import factorial

import sympy as sp
from sympy import Rational as R
from sympy import Symbol, rf, symbols

import review_large_l_common as C

KNOWN = ["mut_Z", "mut_V_NZ", "mut_Vhat", "mut_S", "mut_P1", "mut_P2", "mut_lead", "mut_count_psi", "mut_count_N"]

r = Symbol("r", positive=True)
rho = Symbol("rho", positive=True)
cc = Symbol("c", positive=True)


def chi_poly(l, mut):
    """Lemma A: chi = sum_{n=1}^{l+2} c_n r^n, c_n = (n-2)((n-4) l(l+1) + n-1) (l-n+3)_{2n-4} / (2^n (n-1)!^2); c_1 = 3/2."""
    L = l * (l + 1)
    s = 0
    for n in range(1, l + 3):
        if n == 1:
            cn = R(3, 2)
        else:
            last = n if mut == "mut_Z" else n - 1
            cn = (n - 2) * ((n - 4) * L + last) * rf(l - n + 3, 2 * n - 4) / (2 ** n * factorial(n - 1) ** 2)
        s += cn * r ** n
    return sp.expand(s)


def Vz(l, lam=None):
    f = 1 + 2 / r
    lam = R((l + 2) * (l - 1), 2) if lam is None else lam
    return 2 * f * (lam ** 2 * r ** 2 * ((lam + 1) * r - 3) + 9 * (lam * r - 1)) / (r ** 3 * (lam * r - 3) ** 2)


def dx(u):
    return (1 + 2 / r) * sp.diff(u, r)


def NZ_sym(rho_, eps_, mut):
    a2 = R(3, 4) if mut != "mut_V_NZ" else 1
    return (1 + 2 * eps_) * rho_ ** 3 / 8 - a2 * rho_ ** 2 + R(9, 2) * rho_ - 9


def main():
    mut, _ = C.parse_args(sys.argv, KNOWN)
    chk = C.Checker(failfast=mut is not None)
    print(f"review_large_l_theorem3, mutation={mut}")
    # 0. planted problem of the common bracket code (not used here, kept so every script starts with it)
    if not C.run_planted(chk):
        print("PLANTED PROBLEM FAILED, stop")
        sys.exit(2)

    # (Z): psi_0,xx = V psi_0, exact, l = 2..14
    for l in range(2, 15):
        c = (l + 2) * (l - 1)
        psi = chi_poly(l, mut) / (c * r - 6)
        res = sp.cancel(sp.together(dx(dx(psi)) - Vz(l) * psi))
        chk.check(f"(Z) psi_0,xx = V psi_0 exactly, l={l}", res == 0)

    # V in (rho, c) variables
    eps = 1 / cc
    f_rho = cc * (2 + eps * rho) / rho  # f = 1 + 2/r with r = rho/c
    V_rho = Vz(0, lam=cc / 2).subs(r, rho / cc)
    claim = 8 * cc ** 4 * (2 + eps * rho) * NZ_sym(rho, eps, mut) / (rho ** 4 * (rho - 6) ** 2)
    chk.check("V(r = rho/c, lambda = c/2) = 8c^4 (2+eps rho) N_Z / (rho^4 (rho-6)^2), symbolic in c", sp.simplify(V_rho - claim) == 0)

    # generic X: Vhat = g^2 - g' = 2 g^2 - psi_xx/psi, and the S identity
    Xf = sp.Function("X")
    X = Xf(rho)
    X1 = sp.diff(X, rho)
    N = X - (rho - 6) * X1
    g = -cc * f_rho * N / ((rho - 6) * X)
    dxr = lambda u: f_rho * cc * sp.diff(u, rho)  # d/dx = f d/dr = f c d/drho
    psi0 = X / (cc * (rho - 6))  # chi = X(c r)/c, so psi_0 = X/(c (rho-6))
    g_def = dxr(psi0) / psi0
    chk.check("g = psi_0'/psi_0 equals -c f N/((rho-6) X) for generic X", sp.simplify(g_def - g) == 0)
    Vhat = g ** 2 - dxr(g)
    chk.check("Vhat = g^2 - g' = 2 g^2 - psi_xx/psi_0 for generic X",
              sp.simplify(Vhat - (2 * g ** 2 - dxr(dxr(psi0)) / psi0)) == 0)
    sgn_eps = -1 if mut == "mut_S" else 1
    Np = sp.diff(N, rho)
    S = rho * X * (sgn_eps * eps * N + (2 + eps * rho) * Np) - (2 + eps * rho) * N * (X + 2 * rho * X1)
    ratio = sp.factor(sp.simplify(Vhat / S))
    chk.check("Vhat = S c^3 (2c + rho)/(rho^3 (rho-6) X^2) for generic X (S of RESULTS 4.3)",
              sp.simplify(ratio - cc ** 3 * (2 * cc + rho) / (rho ** 3 * (rho - 6) * X ** 2)) == 0, f"ratio={ratio}")
    print("  so sign(Vhat) = sign(S) for rho > 6 where X > 0")

    # Psi identity: Vhat = 2g^2 - V and V = 8c^4(2+eps rho) N_Z/(...) give Vhat = 2c^4 (2+eps rho) Psi/(rho^4 (rho-6)^2 X^2)
    rho_pow = 1 if mut == "mut_Vhat" else 2
    Psi = rho ** rho_pow * (2 + eps * rho) * N ** 2 - 4 * NZ_sym(rho, eps, None) * X ** 2
    form = 2 * cc ** 4 * (2 + eps * rho) * Psi / (rho ** 4 * (rho - 6) ** 2 * X ** 2)
    chk.check("2 g^2 - V = 2c^4 (2+eps rho) Psi / (rho^4 (rho-6)^2 X^2), generic X, symbolic in c",
              sp.simplify(2 * g ** 2 - V_rho - form) == 0)

    # Vhat identity exact from the definition g^2 - g' (no use of (Z)), l = 2..9
    for l in range(2, 10):
        c = (l + 2) * (l - 1)
        chi = chi_poly(l, None)
        psi = chi / (c * r - 6)
        gl = dx(psi) / psi
        Vh = gl ** 2 - dx(gl)
        Xr = c * chi  # X(rho) = c chi(rho/c) = c chi(r) at rho = c r
        Xp = sp.diff(c * chi.subs(r, rho / c), rho).subs(rho, c * r)
        rr = c * r
        Nn = Xr - (rr - 6) * Xp
        nz = NZ_sym(rr, R(1, c), None)
        Psil = rho_pow_apply(rr, mut) * (2 + rr / c) * Nn ** 2 - 4 * nz * Xr ** 2
        formula = 2 * c ** 4 * (2 + rr / c) * Psil / (rr ** 4 * (rr - 6) ** 2 * Xr ** 2)
        chk.check(f"Vhat from the definition = 2c^4 (2+eps rho) Psi / (rho^4 (rho-6)^2 X^2), exact, l={l}", sp.cancel(sp.together(Vh - formula)) == 0)

    # (P1): numerator of V is a polynomial in (sigma, mu) with positive coefficients, lam r = 3 + sigma, mu = 1/lam
    sig, mu = symbols("sigma mu", positive=True)
    lam = Symbol("lambda", positive=True)
    sgn9 = -1 if mut == "mut_P1" else 1
    numV = lam ** 2 * r ** 2 * ((lam + 1) * r - 3) + sgn9 * 9 * (lam * r - 1)
    nv = sp.expand(sp.simplify(numV.subs(lam, 1 / mu).subs(r, (3 + sig) * mu)))
    Pv = sp.Poly(nv, sig, mu)
    chk.check("(P1) numerator of V = polynomial in (sigma, mu) with all positive coefficients", all(co > 0 for _, co in Pv.terms()),
              f"{nv}")
    print("  (P1) other factors: f = 1 + 2/r > 0, r^3 (lam r - 3)^2 = r^3 sigma^2 > 0 for sigma > 0, lam = c/2 >= 2 > 0; hence V > 0 on (r_s, oo)")

    # (P2): V_r = -2Q/(r^5 (cr-6)^3), Q quintic; c^5 Q(tau, m) with r = (6+tau)/c, c = (m+4)(m+1)
    Vr = sp.diff(Vz(0, lam=cc / 2), r)
    Q = sp.cancel(sp.together(-Vr * r ** 5 * (cc * r - 6) ** 3 / 2))
    Qp = sp.Poly(Q, r)
    chk.check("(P2) Q = -V_r r^5 (cr-6)^3/2 is a polynomial of degree 5 in r", Qp.degree() == 5 and Q.is_polynomial(r), f"Q = {sp.factor(Q)}")
    tau, m = symbols("tau m", nonnegative=True)
    shift = 5 if mut == "mut_P2" else 6
    Qtm = sp.expand(sp.cancel((cc ** 5 * Q).subs(r, (shift + tau) / cc).subs(cc, (m + 4) * (m + 1))))
    P2 = sp.Poly(Qtm, tau, m)
    terms = P2.terms()
    chk.check("(P2) c^5 Q(tau, m) has all coefficients positive", all(co > 0 for _, co in terms), f"{len(terms)} coefficients, RESULTS says 64")
    chk.check("(P2) count of coefficients = 64 as stated in RESULTS", len(terms) == 64)
    print("  (P2) hence Q > 0 for tau >= 0, m >= 0 (all real l >= 2), V_r < 0 on (r_s, oo), and V_x = f V_r < 0")

    # N_Z(6 + tau)
    tau2, e2 = symbols("tau2 e2", nonnegative=True)
    nzs = sp.expand(NZ_sym(6 + tau2, e2, None))
    chk.check("N_Z(6 + tau) has all coefficients positive in (tau, eps)", all(co > 0 for _, co in sp.Poly(nzs, tau2, e2).terms()), str(nzs))

    # leading coefficient of Psi at infinity, symbolic in l: X ~ a rho^(l+2), N ~ -(l+1) a rho^(l+2)
    ls = Symbol("l", positive=True)
    e_l = 1 / ((ls + 2) * (ls - 1))
    lead = e_l * (ls + 1) ** 2 - (1 + 2 * e_l) / 2
    if mut == "mut_lead":
        lead = e_l * (ls + 1) ** 2 - (1 + 3 * e_l) / 2
    chk.check("lead coefficient of Psi / a^2 = eps (l+1)^2 - (1+2 eps)/2 = (l+1)/(2(l-1)) > 0", sp.simplify(lead - (ls + 1) / (2 * (ls - 1))) == 0)
    cl = (ls + 2) * (ls - 1)
    # Vhat ~ 2 c^4 (eps rho) lead rho^(2l+7) a^2 / (rho^6 a^2 rho^(2l+4)) = 2 c^4 eps lead / rho^2 with rho = c r
    vr2 = sp.simplify(2 * cl ** 4 * e_l * lead / cl ** 2)
    chk.check("Vhat r^2 -> (l+1)(l+2) at infinity (L3 of RESULTS), symbolic in l", sp.simplify(vr2 - (ls + 1) * (ls + 2)) == 0, f"{sp.factor(vr2)}")

    # exact root counts, l = 2..16
    for l in range(2, 17):
        c = (l + 2) * (l - 1)
        chi = chi_poly(l, None)
        Xp = sp.expand(c * chi.subs(r, rho / c))
        Np_ = sp.expand(Xp - (rho - 6) * sp.diff(Xp, rho))
        Npoly = sp.Poly(Np_, rho)
        cnt_N = Npoly.count_roots(6, None)
        ints = Npoly.intervals(inf=6, eps=R(1, 10 ** 8))
        (a_, b_), mult = ints[0]
        e_ = R(1, c)
        Psip = sp.Poly(sp.expand(rho ** 2 * (2 + e_ * rho) * Np_ ** 2 - 4 * NZ_sym(rho, e_, None) * Xp ** 2), rho)
        lo_for_psi = 6 if mut == "mut_count_psi" else b_
        cnt_Psi = Psip.count_roots(lo_for_psi, None)
        lo_for_N = (b_ + 1) if mut == "mut_count_N" else 6
        cnt_N2 = Npoly.count_roots(lo_for_N, None)
        chk.check(f"l={l}: N has exactly one root in (6, oo) (T2, exact Sturm count), root in [{float(a_):.8f}, {float(b_):.8f}]", cnt_N2 == 1 and len(ints) == 1)
        chk.check(f"l={l}: Psi has exactly one root in [{float(b_):.8f}, oo) (Theorem 3, exact Sturm count) and Psi < 0 there at the left end",
                  cnt_Psi == 1 and Psip.eval(b_) < 0, f"degree {Psip.degree()}")
    chk.finish()


def rho_pow_apply(x, mut):
    return x if mut == "mut_Vhat" else x ** 2


if __name__ == "__main__":
    main()

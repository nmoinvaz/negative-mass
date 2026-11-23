#!/usr/bin/env python3
"""Review of RESULTS.md 4.4, Theorem L (issue 104), with constants of my own, all inequalities in mpmath iv.

Statement reviewed: for every integer l >= 138, |c r_* - rho_*| <= K_*/c and |c r_neg - rho_neg| <= K_neg/c.
Steps redone: (1) the coefficient bound |a_n - alpha_n| <= eps B_n, (2) the perturbation bounds |N_l - N_inf| <= eps E_N and
|Psi_l - Psi_inf| <= eps eta on a window, (3) the slopes of N_inf and Psi_inf, (4) the localisation lemma.
Assumed, not re-proved here: T2 (N has exactly one root beyond rho = 6), Theorem 3 (V-hat has exactly one zero beyond r_*),
the identity V-hat <-> Psi.  Hand lemmas (stated in review_large_l.md): C1 (coefficient bound), C2 (majorants), C3 (localisation).
Design choices fixed before the first run: eps0 = 1/(140*137) (l >= 138), windows rho_* +- 0.01 and rho_neg +- 0.03, 40 subintervals,
60 explicit terms then a tail majorant.  `--mutate=NAME` breaks one ingredient (exit non-zero).
"""
import sys
from fractions import Fraction
from math import factorial

import sympy as sp
from mpmath import iv, mp, mpf

import review_large_l_common as C
import review_large_l_brackets as B
import review_large_l_roots as RT

KNOWN = ["beta_scale", "slope_sign", "window_small", "psi_identity"]

LMIN = 138
EPS0 = iv.mpf(1) / (140 * 137)  # 1/c(138)
N1 = 60
PARTS = 40
DELTA_N = iv.mpf("0.01")
DELTA_P = iv.mpf("0.03")
RESULTS_L = (iv.mpf("91.51"), iv.mpf("229.87"))  # #2's constants, for the comparison only


def e1(n):
    """sum_{a=0}^{n-3} q_a = (n-2)(n^2-4n-3)/3, an integer."""
    v = (n - 2) * (n * n - 4 * n - 3)
    assert v % 3 == 0
    return v // 3


def g(n):
    return iv.mpf(1) / iv.mpf(2 ** n * factorial(n - 1) ** 2)


def B_n(n, eps0, uniform=False):
    """Upper bound (point) of the coefficient bound B_n in |a_n - alpha_n| <= eps B_n:
    beta_n = alpha'_n (1+2 eps0) + |alpha_n| max(2, S_n (1+2 eps0)), S_n = e1(n) + 2, valid for n <= l+2;
    for n >= l+3 the bound is n^2 |alpha_n|.  uniform=True returns the maximum of the two (valid for every l >= 138 and every n);
    uniform=False returns beta_n (valid for n <= 60 <= l+2, the explicit terms of the sums)."""
    kap = 1 + 2 * eps0
    a = abs(C.alpha(n))
    ap = C.alpha_p(n)
    S = iv.mpf(e1(n) + 2)
    mx = iv.mpf(max(2, (S * kap).b))
    beta = ap * kap + a * mx
    other = n * n * a
    val = iv.mpf(max(beta.b, other.b)) if uniform else iv.mpf(beta.b)
    if C.MUT == "beta_scale":
        val = val * iv.mpf("0.001")
    return val


def sum_EX(rho_hi, n1, eps0):
    s = iv.mpf(0)
    for n in range(3, n1 + 1):
        s = s + B_n(n, eps0) * rho_hi ** n
    s = s + C.tail_majorant(5, 0, rho_hi, n1)  # B_n <= n^5 g_n for n > n1 (lemma C2)
    return C.ub(s)


def sum_EN(rho_lo, rho_hi, n1, eps0):
    s = iv.mpf(0)
    for n in range(3, n1 + 1):
        m = iv.mpf(max(abs((6 * n - (n - 1) * rho_lo).a), abs((6 * n - (n - 1) * rho_lo).b),
                       abs((6 * n - (n - 1) * rho_hi).a), abs((6 * n - (n - 1) * rho_hi).b)))
        s = s + B_n(n, eps0) * rho_hi ** (n - 1) * m
    s = s + C.tail_majorant(6, 0, rho_hi, n1)  # B_n rho^(n-1) |6n-(n-1)rho| <= n^6 g_n rho_hi^n
    return C.ub(s)


def exact_pair(l, rho):
    """X_l, X_l', X_l'' (iv, all l+2 terms) and the derived quantities at an iv point."""
    ev = B.Evaluator(l)
    return ev.q(rho)


def derive(chk, mut, verbose=True):
    """Steps 1 to 4.  Returns a dict with the windows, the roots and the constants K_*, K_neg."""
    # roots and windows
    r_star = mp.findroot(RT.N_inf_mp, mpf("8.4857"), tol=mpf(10) ** -50)
    r_neg = mp.findroot(RT.Psi_inf_mp, mpf("9.636"), tol=mpf(10) ** -50)
    zs = [C.endpoint(int(mp.floor(r_star * mpf(10) ** 32)) + d, 10 ** 32) for d in (-100, 101)]
    zn = [C.endpoint(int(mp.floor(r_neg * mpf(10) ** 32)) + d, 10 ** 32) for d in (-100, 101)]
    dN = DELTA_N if mut != "window_small" else iv.mpf("0.0001")
    WN = (zs[0] - dN, zs[1] + dN)
    WP = (zn[0] - DELTA_P, zn[1] + DELTA_P)
    print(f"rho_* = {mp.nstr(r_star, 22)}  rho_neg = {mp.nstr(r_neg, 22)};  windows [{C.fmt(WN[0], 12)}, {C.fmtu(WN[1], 12)}] and [{C.fmt(WP[0], 12)}, {C.fmtu(WP[1], 12)}]")
    eps0b = C.ub(EPS0)

    # ---- hand lemma C2: majorants B_n <= n^5 g_n, |alpha^1_n| <= n^5 g_n (used for tails), polynomial positivity for n >= 5 + explicit iv for 3 <= n <= 400
    n_, t_ = sp.symbols("n t")
    kap = sp.Rational(1002, 1000)
    e1s = (n_ - 2) * (n_ ** 2 - 4 * n_ - 3) / 3
    P_B = n_ ** 5 - (3 * (n_ - 2) * (n_ - 3) * kap + (n_ - 2) * (n_ - 4) * (2 + (e1s + 2) * kap))
    co = sp.Poly(sp.expand(P_B.subs(n_, 5 + t_)), t_).all_coeffs()
    chk.check("C2: n^5 - [alpha'_n (1+2eps0) + alpha_n (2 + S_n (1+2eps0))] g_n^-1 has all coefficients >= 0 in n = 5 + t (kappa = 1.002 >= 1+2 eps0)", all(c >= 0 for c in co))
    ok = True
    for n in range(3, 401):
        ok &= (iv.mpf(n) ** 5 * g(n) - B_n(n, eps0b, uniform=True)).a > 0
    chk.check("C2: B_n <= n^5 g_n for 3 <= n <= 400 (iv, with eps0 = 1/(140*137))", ok)

    # ---- step 1: the coefficient bound |a_n - alpha_n| <= eps B_n, checked in iv for many l and all n up to l + 10
    for l in (138, 139, 140, 145, 200, 276, 500, 1000):
        eps = iv.mpf(1) / C.c_of(l)
        a = C.coeffs(l)
        worst = 0
        ok = True
        for n in range(3, l + 11):
            an = a[n] if n < len(a) else iv.mpf(0)
            diff = abs(an - C.alpha(n))
            bound = eps * B_n(n, eps0b, uniform=(n > l + 2))
            # the bound is an equality at n = 4 for l = 138 (eps = eps0), so the test is
            # "not definitely violated": violated means the enclosure of diff - bound lies above 0
            ok &= not (diff - bound).a > 0
            worst = max(worst, float((diff / bound).b))
        chk.check(f"step 1: |a_n - alpha_n| <= eps B_n for 3 <= n <= l+10, l={l}", ok, f"largest ratio {worst:.6f} (1 only at n = 4, l = 138 where eps = eps0)")

    # ---- step 2: perturbation bounds
    # algebraic identity Psi_l - Psi_inf = eps T
    N0, DN, X0, DX, e, rr, NZ0 = sp.symbols("N0 DN X0 DX e rho NZ0")
    Nl, Xl = N0 + e * DN, X0 + e * DX
    Psil = rr ** 2 * (2 + e * rr) * Nl ** 2 - 4 * (NZ0 + e * rr ** 3 / 4) * Xl ** 2
    Psii = 2 * rr ** 2 * N0 ** 2 - 4 * NZ0 * X0 ** 2
    T = rr ** 3 * Nl ** 2 + 2 * rr ** 2 * DN * (Nl + N0) - rr ** 3 * Xl ** 2
    if mut != "psi_identity":
        T = T - 4 * NZ0 * DX * (Xl + X0)
    chk.check("Psi_l - Psi_inf = eps [rho^3 N_l^2 + 2 rho^2 D_N (N_l + N_inf) - rho^3 X_l^2 - 4 N_Z^0 D_X (X_l + X_inf)] (symbolic)",
              sp.expand(Psil - Psii - e * T) == 0)
    rho_s, e_s = sp.symbols("rho_s e_s")
    nz = lambda e_: (1 + 2 * e_) * rho_s ** 3 / 8 - sp.Rational(3, 4) * rho_s ** 2 + sp.Rational(9, 2) * rho_s - 9
    chk.check("N_Z(rho, eps) - N_Z(rho, 0) = eps rho^3/4 (symbolic)", sp.expand(nz(e_s) - nz(0) - e_s * rho_s ** 3 / 4) == 0)

    EX_N = sum_EX(WN[1], N1, eps0b)
    EN_N = sum_EN(WN[0], WN[1], N1, eps0b)
    EX_P = sum_EX(WP[1], N1, eps0b)
    EN_P = sum_EN(WP[0], WP[1], N1, eps0b)
    print(f"E_X, E_N on the N window: {C.fmtu(EX_N, 10)}, {C.fmtu(EN_N, 10)}")
    print(f"E_X, E_N on the Psi window: {C.fmtu(EX_P, 10)}, {C.fmtu(EN_P, 10)}")
    # Psi: eta = sup |T| over subintervals
    eta = iv.mpf(0)
    mu_P = None
    for w in C.subintervals(WP[0], WP[1], PARTS):
        q = C.quantities(C.jet_inf(w, 2), w, iv.mpf(0))
        n0 = q["N"]
        x0 = q["X"]
        nz0 = C.NZ(w, 0)
        dni = iv.mpf([-EN_P.b, EN_P.b])
        dxi = iv.mpf([-EX_P.b, EX_P.b])
        nl = n0 + eps0b * dni
        xl = x0 + eps0b * dxi
        Tv = w ** 3 * nl ** 2 + 2 * w ** 2 * dni * (nl + n0) - w ** 3 * xl ** 2 - 4 * nz0 * dxi * (xl + x0)
        eta = iv.mpf(max(eta.b, abs(Tv.a), abs(Tv.b)))
        pa = q["Psip"].a
        mu_P = pa if mu_P is None or pa < mu_P else mu_P
    # N: slope -N_inf' = (rho-6) X_inf''
    mu_N = None
    for w in C.subintervals(WN[0], WN[1], PARTS):
        q = C.quantities(C.jet_inf(w, 2), w, iv.mpf(0))
        s = -q["Np"] if mut != "slope_sign" else q["Np"]
        sa = s.a
        mu_N = sa if mu_N is None or sa < mu_N else mu_N
    print(f"eta (Psi) <= {C.fmtu(eta)}; E_N <= {C.fmtu(EN_N)}; slopes: -N_inf' >= {C.fmt(mu_N)}, Psi_inf' >= {C.fmt(mu_P)}")
    chk.check("step 3: -N_inf' >= mu_N > 0 on the N window", mu_N > 0, f"mu_N = {C.fmt(mu_N)}")
    chk.check("step 3: Psi_inf' >= mu_P > 0 on the Psi window", mu_P > 0, f"mu_P = {C.fmt(mu_P)}")
    # ---- step 4: localisation, F' >= mu on [z-delta, z+delta], |F_l - F| <= eta eps < mu delta  =>  F_l has a zero within eta eps/mu of z
    K_N = iv.mpf((iv.mpf(EN_N.b) / mu_N).b) if mu_N > 0 else None
    K_P = iv.mpf((iv.mpf(eta.b) / mu_P).b) if mu_P > 0 else None
    if K_N is None or K_P is None:
        chk.check("localisation constants defined", False)
        chk.finish()
    print(f"K_* = E_N/mu_N = {C.fmtu(K_N)}   K_neg = eta/mu_P = {C.fmtu(K_P)}")
    chk.check("step 4: eps0 E_N < mu_N delta_N (the zero of N_l lies in the window)", (EN_N * eps0b).b < (mu_N * dN).a)
    chk.check("step 4: eps0 eta < mu_P delta_P (the zero of Psi_l lies in the window)", (eta * eps0b).b < (mu_P * DELTA_P).a)
    chk.check("step 4: the windows are disjoint, rho_* + delta_N < rho_neg - delta_P", WN[1].b < WP[0].a)
    chk.check("step 4: the N root of N_l (at most rho_* + K_* eps0) lies below the Psi window, so the Psi zero in the window is beyond r_*",
              (zs[1] + K_N * eps0b).b < WP[0].a)
    print("")
    print(f"THEOREM L (own constants, l >= {LMIN}):")
    print(f"  |c r_*   - rho_*  | <= {C.fmtu(K_N, 6)} / c    (RESULTS: 91.51)")
    print(f"  |c r_neg - rho_neg| <= {C.fmtu(K_P, 6)} / c    (RESULTS: 229.87)")
    print(f"  comparison (not a pass/fail criterion): own K_* {'<=' if K_N.b <= RESULTS_L[0].a else '>'} 91.51, own K_neg {'<=' if K_P.b <= RESULTS_L[1].a else '>'} 229.87")
    return dict(r_star=r_star, r_neg=r_neg, zs=zs, zn=zn, K_N=K_N, K_P=K_P, EN_N=EN_N, eta=eta, WN=WN, WP=WP, mu_N=mu_N, mu_P=mu_P, eps0b=eps0b, dN=dN)


def main():
    mut, _ = C.parse_args(sys.argv, KNOWN)
    chk = C.Checker(failfast=mut is not None)
    print(f"review_large_l_thm_l, mutation={mut}")
    if not C.run_planted(chk):
        print("PLANTED PROBLEM FAILED, stop before real data")
        sys.exit(2)
    D = derive(chk, mut)
    r_star, r_neg, zs, zn = D["r_star"], D["r_neg"], D["zs"], D["zn"]
    K_N, K_P, EN_N, eta, WN, WP, eps0b = D["K_N"], D["K_P"], D["EN_N"], D["eta"], D["WN"], D["WP"], D["eps0b"]

    # ---- numbers quoted in RESULTS 4.4: -N_inf' >= 2.836 on [8, 9]; Psi_inf'(rho_neg) = 33523, |Psi_inf''| <= 1.31e5 on [9.61, 9.66]
    rho_s2 = sp.Symbol("rho")
    X_f = sp.Function("X")(rho_s2)
    N_f = X_f - (rho_s2 - 6) * sp.diff(X_f, rho_s2)
    NZ0_f = rho_s2 ** 3 / 8 - sp.Rational(3, 4) * rho_s2 ** 2 + sp.Rational(9, 2) * rho_s2 - 9
    Psi_f = 2 * rho_s2 ** 2 * N_f ** 2 - 4 * NZ0_f * X_f ** 2
    P2 = sp.diff(Psi_f, rho_s2, 2)
    x0_, x1_, x2_, x3_ = sp.symbols("x0 x1 x2 x3")
    P2s = P2.subs(sp.Derivative(X_f, (rho_s2, 3)), x3_).subs(sp.Derivative(X_f, (rho_s2, 2)), x2_).subs(sp.Derivative(X_f, rho_s2), x1_).subs(X_f, x0_)
    P2fun = sp.lambdify((rho_s2, x0_, x1_, x2_, x3_), P2s, modules=[{}])
    # Psi'' against a difference quotient of Psi' at a point, in mp
    pt = mpf("9.636")
    h = mpf(10) ** -20
    dq = (C.quantities(C.jet_inf_mp(pt + h, 2), pt + h, mpf(0))["Psip"] - C.quantities(C.jet_inf_mp(pt - h, 2), pt - h, mpf(0))["Psip"]) / (2 * h)
    j3 = C.jet_inf_mp(pt, 3)
    chk.check("Psi_inf'' closed form (sympy d^2/d rho^2 with generic X) agrees with a central difference of Psi_inf' at 9.636",
              abs(P2fun(pt, *j3) - dq) < mpf(10) ** -8 * abs(dq), f"{mp.nstr(P2fun(pt, *j3), 10)} vs {mp.nstr(dq, 10)}")
    sl = None
    for w in C.subintervals(iv.mpf(8), iv.mpf(9), 100):
        qq = C.quantities(C.jet_inf(w, 2), w, iv.mpf(0))
        a_ = (-qq["Np"]).a
        sl = a_ if sl is None or a_ < sl else sl
    chk.check("RESULTS 4.4: -N_inf' = (rho-6) X_inf'' >= 2.836 on [8, 9]", sl >= 2.836, f"inf = {C.fmt(sl)}")
    qz = C.quantities(C.jet_inf(C.hull(zn[0], zn[1]), 2), C.hull(zn[0], zn[1]), iv.mpf(0))
    chk.check("RESULTS 4.4: Psi_inf'(rho_neg) = 33523 to the digits quoted", abs(qz["Psip"].a - 33523) < 1, f"{C.fmt(qz['Psip'], 10)}")
    mx = iv.mpf(0)
    for w in C.subintervals(iv.mpf("9.61"), iv.mpf("9.66"), 100):
        jj = C.jet_inf(w, 3)
        v = P2fun(w, *jj)
        mx = iv.mpf(max(mx.b, abs(v.a), abs(v.b)))
    chk.check("RESULTS 4.4: |Psi_inf''| <= 1.31e5 on [9.61, 9.66]", mx.b <= 1.31e5, f"sup |Psi''| <= {C.fmtu(mx, 7)}")
    E108 = sum_EN(iv.mpf(8), iv.mpf(9), N1, C.ub(iv.mpf(1) / 108))
    print(f"  for the record: my E_N on [8, 9] with eps0 = 1/108 is {C.fmtu(E108, 6)}; RESULTS quotes 259.5 (different bounding, not a criterion)")

    # ---- the final inequalities against my certified brackets at the l RESULTS names, and the sample bounds of step 2
    cnt = {}
    for l in (138, 139, 145, 276, 690, 1000, 2000, 3000):
        res = B.certify_l(l, r_star, r_neg, cnt, chk, mut)
        c = C.c_of(l)
        F = lambda n_: iv.mpf(n_) / iv.mpf(10 ** B.DIG)
        okN = max(abs((F(res["N"][0]) - zs[1]).a), abs((F(res["N"][1]) - zs[0]).b)) * c <= K_N.b
        okP = max(abs((F(res["P"][0]) - zn[1]).a), abs((F(res["P"][1]) - zn[0]).b)) * c <= K_P.b
        kN = (F(res["N"][1]) - zs[0]).b * c
        kP = (F(res["P"][1]) - zn[0]).b * c
        chk.check(f"Theorem L against the certified brackets, l={l}", okN and okP, f"(c r_* - rho_*) c <= {C.fmtu(kN, 6)}, (c r_neg - rho_neg) c <= {C.fmtu(kP, 6)}")
    # step 2 sample check: the bounds |N_l - N_inf| <= eps E_N, |Psi_l - Psi_inf| <= eps eta at points of the windows, exact l
    for l in (138, 139, 276, 1000):
        eps = iv.mpf(1) / C.c_of(l)
        ok = True
        worstN = worstP = 0
        for k in range(0, 11):
            xN = WN[0] + (WN[1] - WN[0]) * k / 10
            xP = WP[0] + (WP[1] - WP[0]) * k / 10
            xN = iv.mpf(xN.mid)
            xP = iv.mpf(xP.mid)
            ql = exact_pair(l, xN)
            qi = C.quantities(C.jet_inf(xN, 2), xN, iv.mpf(0))
            d = abs(ql["N"] - qi["N"]) / eps
            ok &= (EN_N.b - d.b) >= 0
            worstN = max(worstN, float(d.b))
            ql = exact_pair(l, xP)
            qi = C.quantities(C.jet_inf(xP, 2), xP, iv.mpf(0))
            d = abs(ql["Psi"] - qi["Psi"]) / eps
            ok &= (eta.b - d.b) >= 0
            worstP = max(worstP, float(d.b))
        chk.check(f"step 2 at 11 points of each window, l={l}: |N_l - N_inf|/eps <= E_N and |Psi_l - Psi_inf|/eps <= eta", ok, f"largest observed {worstN:.4g} (bound {float(EN_N.b):.4g}), {worstP:.4g} (bound {float(eta.b):.4g})")
    chk.finish()


if __name__ == "__main__":
    main()

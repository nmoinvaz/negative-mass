#!/usr/bin/env python3
"""Review of RESULTS.md 4.5, Theorem M (issue 104), with constants of my own, all inequalities in mpmath iv.

Statement reviewed: for every integer l >= 500, c r_*(l+1) < c r_*(l) and c r_neg(l+1) < c r_neg(l), with
0 < c r_*(l) - rho_* <= 3.709/c, 0 < c r_neg(l) - rho_neg <= 6.101/c, and (l+1)^2 r(l+1) > l^2 r(l) for both ("l^2 r increases").
Second order expansion a_n = alpha_n + eps alpha^1_n + eps^2 omega_n, N_l = N_inf + eps nu + eps^2 zeta,
Psi_l = Psi_inf + eps Psi_1 + eps^2 E, |zeta| <= Z, |E| <= E_max on windows that Theorem L (own constants) gives.
The same machinery with eps_max = 1/c(138) refines Theorem L (configuration "L2").
Assumed, not re-proved here: T2, Theorem 3, the identity V-hat <-> Psi, and Theorem L (own constants, review_large_l_thm_l.py).
`--mutate=NAME` breaks one ingredient (exit non-zero).
"""
import sys
from fractions import Fraction
from math import factorial

import sympy as sp
from mpmath import iv, mp, mpf

import review_large_l_common as C
import review_large_l_brackets as B
import review_large_l_thm_l as TL

KNOWN = ["nu_sign", "omega_scale", "cmp_reverse", "l2r_s", "window_narrow"]

N1 = 60
RESULTS_M = (iv.mpf("3.709"), iv.mpf("6.101"))
RESULTS_L = (iv.mpf("91.51"), iv.mpf("229.87"))


def g(n):
    return iv.mpf(1) / iv.mpf(2 ** n * factorial(n - 1) ** 2)


def alpha1(n):
    """alpha^1_n = alpha'_n - alpha_n e1(n)  (the coefficient of eps in a_n)."""
    if C.MUT == "nu_sign":
        return C.alpha_p(n) + C.alpha(n) * TL.e1(n)
    return C.alpha_p(n) - C.alpha(n) * TL.e1(n)


def Omega_n(n, epsb, uniform=False):
    """Upper bound (point) of |omega_n|.  Regime n <= l+2: Oa = |e1| alpha'_n + (|alpha_n| + eps alpha'_n)(e1+4)^2 (1+2eps)/2.
    Regime n >= l+3: Ob = n^4 |alpha_n| + n^2 alpha'_n + n^2 |alpha_n||e1|.  uniform=True is the maximum (every n, every l >= lmin);
    uniform=False is Oa (the explicit terms n <= 60 <= lmin + 2)."""
    a = abs(C.alpha(n))
    ap = C.alpha_p(n)
    E = abs(iv.mpf(TL.e1(n)))
    kap = 1 + 2 * epsb
    sab = iv.mpf(TL.e1(n) + 4)
    Oa = E * ap + (a + epsb * ap) * sab ** 2 * kap / 2
    Ob = n ** 4 * a + n * n * ap + n * n * a * E
    val = iv.mpf(max(Oa.b, Ob.b)) if uniform else iv.mpf(Oa.b)
    if C.MUT == "omega_scale":
        val = val * iv.mpf("0.001")
    return val


def nu_enc(w, n1=N1):
    s = iv.mpf(0)
    for n in range(3, n1 + 1):
        s = s + alpha1(n) * w ** (n - 1) * (6 * n - (n - 1) * w)
    T = C.tail_majorant(6, 0, w.b, n1)
    return s + iv.mpf([-T.b, T.b])


def xi_enc(w, n1=N1):
    s = iv.mpf(0)
    for n in range(3, n1 + 1):
        s = s + alpha1(n) * w ** n
    T = C.tail_majorant(5, 0, w.b, n1)
    return s + iv.mpf([-T.b, T.b])


def Z_sum(wlo, whi, epsb, n1=N1):
    s = iv.mpf(0)
    for n in range(3, n1 + 1):
        m = iv.mpf(max(abs((6 * n - (n - 1) * wlo).a), abs((6 * n - (n - 1) * wlo).b),
                       abs((6 * n - (n - 1) * whi).a), abs((6 * n - (n - 1) * whi).b)))
        s = s + Omega_n(n, epsb) * whi.b ** (n - 1) * m
    return C.ub(s + C.tail_majorant(9, 0, whi.b, n1))


def Th_sum(whi, epsb, n1=N1):
    s = iv.mpf(0)
    for n in range(3, n1 + 1):
        s = s + Omega_n(n, epsb) * whi.b ** n
    return C.ub(s + C.tail_majorant(8, 0, whi.b, n1))


def psi1(u, v, x, y, w, nz0):
    """Coefficient of eps in Psi_l: 4 rho^2 u v + rho^3 u^2 - rho^3 x^2 - 8 N_Z^0 x y."""
    return 4 * w ** 2 * u * v + w ** 3 * u ** 2 - w ** 3 * x ** 2 - 8 * nz0 * x * y


def E_rem(u, v, x, y, zt, th, w, nz0, eps):
    """(Psi_l - Psi_inf - eps Psi_1)/eps^2 with N_l = u + eps v + eps^2 zt, X_l = x + eps y + eps^2 th."""
    c2 = 2 * w ** 2 * (v ** 2 + 2 * u * zt) + 2 * w ** 3 * u * v - 4 * nz0 * (y ** 2 + 2 * x * th) - 2 * w ** 3 * x * y
    c3 = 4 * w ** 2 * v * zt + w ** 3 * (v ** 2 + 2 * u * zt) - 8 * nz0 * y * th - w ** 3 * (y ** 2 + 2 * x * th)
    c4 = 2 * w ** 2 * zt ** 2 + 2 * w ** 3 * v * zt - 4 * nz0 * th ** 2 - 2 * w ** 3 * y * th
    c5 = w ** 3 * zt ** 2 - w ** 3 * th ** 2
    return c2 + eps * c3 + eps ** 2 * c4 + eps ** 3 * c5


def symbolic_checks(chk):
    u, v, zt, x, y, th, e, r, nz0 = sp.symbols("u v zt x y th e r nz0")
    Nl = u + e * v + e ** 2 * zt
    Xl = x + e * y + e ** 2 * th
    Psil = r ** 2 * (2 + e * r) * Nl ** 2 - 4 * (nz0 + e * r ** 3 / 4) * Xl ** 2
    P0 = 2 * r ** 2 * u ** 2 - 4 * nz0 * x ** 2
    P1 = 4 * r ** 2 * u * v + r ** 3 * u ** 2 - r ** 3 * x ** 2 - 8 * nz0 * x * y
    Er = (2 * r ** 2 * (v ** 2 + 2 * u * zt) + 2 * r ** 3 * u * v - 4 * nz0 * (y ** 2 + 2 * x * th) - 2 * r ** 3 * x * y
          + e * (4 * r ** 2 * v * zt + r ** 3 * (v ** 2 + 2 * u * zt) - 8 * nz0 * y * th - r ** 3 * (y ** 2 + 2 * x * th))
          + e ** 2 * (2 * r ** 2 * zt ** 2 + 2 * r ** 3 * v * zt - 4 * nz0 * th ** 2 - 2 * r ** 3 * y * th)
          + e ** 3 * (r ** 3 * zt ** 2 - r ** 3 * th ** 2))
    chk.check("Psi_l = Psi_inf + eps Psi_1 + eps^2 E(eps) with the typed Psi_1 and E (symbolic)", sp.expand(Psil - P0 - e * P1 - e ** 2 * Er) == 0)
    l = sp.Symbol("l", positive=True)
    c = lambda k: (k + 2) * (k - 1)
    chk.check("(l+1)^2 c(l) - l^2 c(l+1) = l^2 - 3l - 2 (symbolic)", sp.expand((l + 1) ** 2 * c(l) - l ** 2 * c(l + 1) - (l ** 2 - 3 * l - 2)) == 0)
    chk.check("c(l+1) - c(l) = 2l + 2 (symbolic)", sp.expand(c(l + 1) - c(l) - (2 * l + 2)) == 0)
    chk.check("(l^2 - 3l - 2)/(l+1)^2 is increasing for l >= 1: derivative has the sign of 5l + 1 > 0 (symbolic)",
              sp.simplify(sp.diff((l ** 2 - 3 * l - 2) / (l + 1) ** 2, l) * (l + 1) ** 3 - (5 * l + 1)) == 0)
    n, t = sp.symbols("n t")
    e1s = (n - 2) * (n ** 2 - 4 * n - 3) / 3
    eps0 = sp.Rational(1, 1000)
    Pa = n ** 8 - (e1s * 3 * (n - 2) * (n - 3) + ((n - 2) * (n - 4) + eps0 * 3 * (n - 2) * (n - 3)) * (e1s + 4) ** 2 * (1 + 2 * eps0) / 2)
    Pb = n ** 8 - (n ** 4 * (n - 2) * (n - 4) + n ** 2 * 3 * (n - 2) * (n - 3) + n ** 2 * (n - 2) * (n - 4) * e1s)
    Pc = n ** 5 - (3 * (n - 2) * (n - 3) + (n - 2) * (n - 4) * e1s)
    for name, P in (("Omega_a", Pa), ("Omega_b", Pb), ("alpha^1", Pc)):
        co = sp.Poly(sp.expand(P.subs(n, 5 + t)), t).all_coeffs()
        chk.check(f"C2: majorant polynomial for {name} has all coefficients >= 0 in n = 5 + t (eps_max <= 1/1000)", all(cc >= 0 for cc in co))


def majorant_checks(chk, epsb_all):
    ok1 = ok2 = True
    for n in range(3, 401):
        ok1 &= (iv.mpf(n) ** 8 * g(n) - Omega_n(n, epsb_all, uniform=True)).a > 0
        ok2 &= (iv.mpf(n) ** 5 * g(n) - abs(alpha1(n))).a > 0
    chk.check("C2: Omega_n <= n^8 g_n and |alpha^1_n| <= n^5 g_n for 3 <= n <= 400 (iv)", ok1 and ok2)


def run_config(chk, D, name, lmin, KpN, KpP, parts, do_monotone, samples, mut):
    print(f"\n=== configuration {name}: l >= {lmin} ===")
    eps_max = iv.mpf(1) / C.c_of(lmin)
    epsb = C.ub(eps_max)
    chk.check(f"{name}: window radii K' = {KpN}, {KpP} are at least Theorem L's own constants K_*, K_neg", KpN >= D["K_N"].b and KpP >= D["K_P"].b,
              f"{C.fmtu(D['K_N'], 6)}, {C.fmtu(D['K_P'], 6)}")
    scaleN = iv.mpf("0.01") if mut == "window_narrow" else iv.mpf(1)
    halfN = KpN * epsb * scaleN
    halfP = KpP * epsb * scaleN
    zs, zn = D["zs"], D["zn"]
    WN = (zs[0] - halfN, zs[1] + halfN)
    WP = (zn[0] - halfP, zn[1] + halfP)
    print(f"  windows: N [{C.fmt(WN[0], 12)}, {C.fmtu(WN[1], 12)}], Psi [{C.fmt(WP[0], 12)}, {C.fmtu(WP[1], 12)}], eps_max = {C.fmtu(epsb, 6)}")
    ZN = Z_sum(WN[0], WN[1], epsb)
    ZP = Z_sum(WP[0], WP[1], epsb)
    ThP = Th_sum(WP[1], epsb)
    print(f"  Z (N window) <= {C.fmtu(ZN, 6)};  Z, Theta (Psi window) <= {C.fmtu(ZP, 6)}, {C.fmtu(ThP, 6)}")
    # N side
    nus = [nu_enc(w) for w in C.subintervals(WN[0], WN[1], parts)]
    nu_min = min(x.a for x in nus)
    nu_max = max(x.b for x in nus)
    nu_root = nu_enc(C.hull(zs[0], zs[1]))
    mu_N = None
    for w in C.subintervals(WN[0], WN[1], parts):
        q = C.quantities(C.jet_inf(w, 2), w, iv.mpf(0))
        a = (-q["Np"]).a
        mu_N = a if mu_N is None or a < mu_N else mu_N
    print(f"  nu in [{C.fmt(nu_min)}, {C.fmtu(nu_max)}]  nu(rho_*) in [{C.fmt(nu_root.a)}, {C.fmtu(nu_root.b)}]  -N_inf' >= {C.fmt(mu_N)}")
    chk.check(f"{name}: -N_inf' >= mu_N > 0 on the N window and nu_min > 0", mu_N > 0 and nu_min > 0, f"mu_N = {C.fmt(mu_N)}, nu_min = {C.fmt(nu_min)}")
    chk.check(f"{name}: nu(rho_*) - eps_max Z > 0, so N_l(rho_*) > 0 and c r_* > rho_*", (nu_root.a - epsb * ZN.b) > 0)
    sN = iv.mpf((((nu_max + eps_max * ZN)) / mu_N).b)
    chk.check(f"{name}: s_N eps_max <= window radius, so rho_* + s_N eps lies in the N window", (sN * eps_max).b <= halfN.a)
    # Psi side
    p1s = []
    mu_P = None
    Emax = iv.mpf(0)
    for w in C.subintervals(WP[0], WP[1], parts):
        q = C.quantities(C.jet_inf(w, 2), w, iv.mpf(0))
        nz0 = C.NZ(w, 0)
        v = nu_enc(w)
        y = xi_enc(w)
        p1 = psi1(q["N"], v, q["X"], y, w, nz0)
        p1s.append(p1)
        a = q["Psip"].a
        mu_P = a if mu_P is None or a < mu_P else mu_P
        zt = iv.mpf([-ZP.b, ZP.b])
        th = iv.mpf([-ThP.b, ThP.b])
        Ev = E_rem(q["N"], v, q["X"], y, zt, th, w, nz0, iv.mpf([0, epsb.b]))
        Emax = iv.mpf(max(Emax.b, abs(Ev.a), abs(Ev.b)))
    p1_min = min(x.a for x in p1s)
    p1_max = max(x.b for x in p1s)
    wr = C.hull(zn[0], zn[1])
    qr = C.quantities(C.jet_inf(wr, 2), wr, iv.mpf(0))
    p1_root = psi1(qr["N"], nu_enc(wr), qr["X"], xi_enc(wr), wr, C.NZ(wr, 0))
    print(f"  Psi_1 in [{C.fmt(p1_min)}, {C.fmtu(p1_max)}]  Psi_1(rho_neg) in [{C.fmt(p1_root.a)}, {C.fmtu(p1_root.b)}]  Psi_inf' >= {C.fmt(mu_P)}  E_max <= {C.fmtu(Emax)}")
    chk.check(f"{name}: Psi_inf' >= mu_P > 0 on the Psi window", mu_P > 0, f"mu_P = {C.fmt(mu_P)}")
    chk.check(f"{name}: Psi_1 < 0 on the whole Psi window", p1_max < 0)
    chk.check(f"{name}: Psi_1(rho_neg) + eps_max E_max < 0, so Psi_l(rho_neg) < 0 and c r_neg > rho_neg", (p1_root.b + epsb * Emax.b) < 0)
    sP = iv.mpf((((-p1_min) + eps_max * Emax) / mu_P).b)
    chk.check(f"{name}: s_P eps_max <= window radius, so rho_neg + s_P eps lies in the Psi window", (sP * eps_max).b <= halfP.a)
    print(f"  s_N = {C.fmtu(sN, 8)}, s_P = {C.fmtu(sP, 8)}")
    # sample comparison with the exact functions, the expansion and its remainder bounds
    for l in samples:
        eps = iv.mpf(1) / C.c_of(l)
        ok = True
        worstZ = worstE = 0
        for k in range(5):
            xN = iv.mpf((WN[0] + (WN[1] - WN[0]) * k / 4).mid)
            xP = iv.mpf((WP[0] + (WP[1] - WP[0]) * k / 4).mid)
            ql = TL.exact_pair(l, xN)
            qi = C.quantities(C.jet_inf(xN, 2), xN, iv.mpf(0))
            resid = (ql["N"] - qi["N"] - eps * nu_enc(xN)) / (eps * eps)
            ok &= not (abs(resid).a > ZN.b)
            worstZ = max(worstZ, float(abs(resid).b))
            ql = TL.exact_pair(l, xP)
            qi = C.quantities(C.jet_inf(xP, 2), xP, iv.mpf(0))
            nz0 = C.NZ(xP, 0)
            p1 = psi1(qi["N"], nu_enc(xP), qi["X"], xi_enc(xP), xP, nz0)
            resid = (ql["Psi"] - qi["Psi"] - eps * p1) / (eps * eps)
            ok &= not (abs(resid).a > Emax.b)
            worstE = max(worstE, float(abs(resid).b))
        chk.check(f"{name}: expansions at 5 points of each window against the exact N_l and Psi_l, l={l}", ok,
                  f"|zeta| seen <= {worstZ:.4g} (bound {float(ZN.b):.4g}), |E| seen <= {worstE:.4g} (bound {float(Emax.b):.4g})")
    out = dict(sN=sN, sP=sP)
    if do_monotone:
        l0 = lmin
        c0, c1 = C.c_of(l0), C.c_of(l0 + 1)
        ratio = C.ub(iv.mpf(c1) / c0 + iv.mpf(c0) / c1)
        lhsN = (2 * l0 + 2) * nu_min
        rhsN = ratio * ZN
        lhsP = (2 * l0 + 2) * (-p1_max)
        rhsP = ratio * Emax
        if mut == "cmp_reverse":
            ok_n, ok_p = lhsN.b < rhsN.a, lhsP.b < rhsP.a
        else:
            ok_n, ok_p = lhsN.a > rhsN.b, lhsP.a > rhsP.b
        chk.check(f"{name}: (2l+2) nu_min > (c(l+1)/c(l) + c(l)/c(l+1)) Z at l = {l0}, so N_(l+1)(rho_l) < 0 and rho_*(l+1) < rho_*(l) for all l >= {l0}", ok_n,
                  f"{C.fmt(lhsN)} vs {C.fmtu(rhsN)}")
        chk.check(f"{name}: (2l+2) |Psi_1|_min > (...) E_max at l = {l0}, so Psi_(l+1)(rho_l) > 0 and rho_neg(l+1) < rho_neg(l) for all l >= {l0}", ok_p,
                  f"{C.fmt(lhsP)} vs {C.fmtu(rhsP)}")
        # l^2 r criterion: s < rho_inf (l^2 - 3l - 2)/(l+1)^2 at l0, increasing in l
        fac = iv.mpf(l0 * l0 - 3 * l0 - 2) / ((l0 + 1) ** 2)
        s_N, s_P = (sN, sP) if mut != "l2r_s" else (iv.mpf(20), iv.mpf(20))
        chk.check(f"{name}: l^2 r_* increases: s_N < rho_* (l^2-3l-2)/(l+1)^2 at l = {l0}", s_N.b < (zs[0] * fac).a, f"{C.fmtu(s_N)} < {C.fmt(zs[0] * fac)}")
        chk.check(f"{name}: l^2 r_neg increases: s_P < rho_neg (l^2-3l-2)/(l+1)^2 at l = {l0}", s_P.b < (zn[0] * fac).a, f"{C.fmtu(s_P)} < {C.fmt(zn[0] * fac)}")
        print("")
        print(f"THEOREM M (own constants, l >= {lmin}):")
        print(f"  c r_*(l+1) < c r_*(l) and c r_neg(l+1) < c r_neg(l); 0 < c r_* - rho_* <= {C.fmtu(sN, 6)}/c (RESULTS 3.709), 0 < c r_neg - rho_neg <= {C.fmtu(sP, 6)}/c (RESULTS 6.101);")
        print("  (l+1)^2 r_*(l+1) > l^2 r_*(l) and the same for r_neg.")
        print(f"  comparison (not a pass/fail criterion): own s_N {'<=' if sN.b <= RESULTS_M[0].a else '>'} 3.709, own s_P {'<=' if sP.b <= RESULTS_M[1].a else '>'} 6.101")
    else:
        print("")
        print(f"THEOREM L, second order refinement (own constants, l >= {lmin}):")
        print(f"  0 < c r_* - rho_* <= {C.fmtu(sN, 6)}/c and 0 < c r_neg - rho_neg <= {C.fmtu(sP, 6)}/c (first order: {C.fmtu(D['K_N'], 6)}, {C.fmtu(D['K_P'], 6)}; RESULTS: 91.51, 229.87)")
        print(f"  comparison (not a pass/fail criterion): own {'<=' if sN.b <= RESULTS_L[0].a else '>'} 91.51 and {'<=' if sP.b <= RESULTS_L[1].a else '>'} 229.87")
    return out


def main():
    mut, _ = C.parse_args(sys.argv, KNOWN)
    chk = C.Checker(failfast=mut is not None)
    print(f"review_large_l_thm_m, mutation={mut}")
    if not C.run_planted(chk):
        print("PLANTED PROBLEM FAILED, stop before real data")
        sys.exit(2)
    symbolic_checks(chk)
    # Theorem L's own constants (first order), the windows below rest on them
    quiet = C.Checker(failfast=False)
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        D = TL.derive(quiet, None)
    chk.check("first order Theorem L (own constants) re-derived with every check passing", not quiet.fails,
              f"K_* = {C.fmtu(D['K_N'], 6)}, K_neg = {C.fmtu(D['K_P'], 6)}")
    majorant_checks(chk, C.ub(iv.mpf(1) / C.c_of(138)))
    KpN, KpP = 46, 339  # integers above the own first-order constants
    res = {}
    res["L2"] = run_config(chk, D, "L2", 138, KpN, KpP, 40, False, (138, 139, 276), mut)
    res["M"] = run_config(chk, D, "M", 500, KpN, KpP, 8, True, (500, 501, 502, 1000, 1001), mut)

    # exact brackets: the consecutive-l claims of Theorem M
    r_star, r_neg = D["r_star"], D["r_neg"]
    cnt = {}
    F = lambda n_: iv.mpf(n_) / iv.mpf(10 ** B.DIG)
    sN, sP = res["M"]["sN"], res["M"]["sP"]
    for l in (500, 501, 502, 1000, 1001, 2000, 2001):
        r0 = B.certify_l(l, r_star, r_neg, cnt, chk, mut)
        r1 = B.certify_l(l + 1, r_star, r_neg, cnt, chk, mut)
        c0 = C.c_of(l)
        dec_ok = F(r1["N"][1]) < F(r0["N"][0]) and F(r1["P"][1]) < F(r0["P"][0])
        up_ok = ((F(r0["N"][0]) - D["zs"][1]).a > 0 and (F(r0["N"][1]) - D["zs"][0]).b * c0 <= sN.b and
                 (F(r0["P"][0]) - D["zn"][1]).a > 0 and (F(r0["P"][1]) - D["zn"][0]).b * c0 <= sP.b)
        l2 = all((l + 1) ** 2 * F(r1[w][0]) / C.c_of(l + 1) > l ** 2 * F(r0[w][1]) / c0 for w in ("N", "P"))
        # N_(l+1)(lo of r_*(l)) < 0 and Psi_(l+1)(lo of r_neg(l)) > 0
        ev1 = B.Evaluator(l + 1)
        lo_n = C.endpoint(r0["N"][0], 10 ** B.DIG)
        lo_p = C.endpoint(r0["P"][0], 10 ** B.DIG)
        cmp_ok = C.sgn(ev1.q(lo_n)["N"]) == -1 and C.sgn(ev1.q(lo_p)["Psi"]) == 1
        chk.check(f"Theorem M against my certified brackets, l={l} and {l + 1}: decrease, 0 < (c r - rho) c <= s, l^2 r increases, N_(l+1) < 0 and Psi_(l+1) > 0 at the lower ends",
                  dec_ok and up_ok and l2 and cmp_ok, f"(c r_* - rho_*) c <= {C.fmtu((F(r0['N'][1]) - D['zs'][0]) * c0, 5)}, (c r_neg - rho_neg) c <= {C.fmtu((F(r0['P'][1]) - D['zn'][0]) * c0, 5)}")
    chk.finish()


if __name__ == "__main__":
    main()

"""
Theorem M (#2): there is an explicit L1 such that for every integer l >= L1
    rho_*(l+1) < rho_*(l),   rho_neg(l+1) < rho_neg(l),                        (rho = c r)
    rho_*(l) > rho_*,        rho_neg(l) > rho_neg,
    rho_*(l) - rho_* <= s_N / c,   rho_neg(l) - rho_neg <= s_P / c,   with s_N, s_P below,
    l^2 r_*(l+1) > l^2 r_*(l)   and   l^2 r_neg(l+1) > l^2 r_neg(l)   (since s < rho (l^2-3l-2)/(l+1)^2).
Together with threshold_isolation.py (certified brackets, consecutive l <= 3000) this gives strict monotonicity of
c r_* and c r_neg for ALL l >= 2 and of l^2 r_*, l^2 r_neg for all l >= 5, when L1 <= 3000.

Second-order expansion: X_l = x + eps u + eps^2 xi, N_l = n + eps nu + eps^2 zeta, with u = sum alpha1_n rho^n,
alpha1_n = alpha'_n - alpha_n e1(n) (the d/d eps coefficient at 0), |xi|, |zeta| bounded coefficientwise (omega_n below):
  a_n - alpha_n - eps alpha1_n = eps^2 omega_n,  omega_n = alpha_n R_n - alpha'_n e1 + eps alpha'_n R_n,
  |R_n| <= (1/2)(e1 + 4)^2 (1 + 2 eps), because P_n(eps) = prod (1 - q_a eps) = 1 - e1 eps + eps^2 R_n (Taylor, P_n'' bounded by
  (sum |q_a|)^2 (1 + 2 eps), the factors with a >= 1 lying in [0,1]); for n >= l+3, a_n = 0 and c < n^2.
Then  N_{l+1}(rho_l) = (eps' - eps) nu + eps'^2 zeta' - eps^2 zeta < 0  when (l+1) nu_min / 2 > eta2_N,  and the same for
Psi_l = Psi_inf + eps Psi_1 + eps^2 E with Psi_1 < 0.
Assumed from large_l_proof.py (Theorem L): l >= 138, K_* <= 92, K_neg <= 230 (checked against its output file).
Mutations: --mutate=omega | n1sign | s1sign | eta2 | window.
"""
import argparse
import math
import os
import re
import sys
import time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rescaled as R

ap = argparse.ArgumentParser()
ap.add_argument('--mutate', default='')
args = ap.parse_args()
MUT = set(m for m in args.mutate.split(',') if m)
failures = []
t0 = time.time()


def check(name, ok):
    print(("PASS  " if ok else "FAIL  ") + name, flush=True)
    if not ok:
        failures.append(name)


def fmt(x, nd=6):
    return f"{float(x):.{nd}g}"


fact = math.factorial
N1 = 60
K_STAR, K_NEG, L0_THM = Fr(92), Fr(230), 138

# the constants of Theorem L must dominate those printed by large_l_proof.py
with open(os.path.join(HERE, 'large_l_proof_output.txt')) as fh:
    txt = fh.read()
m = re.search(r"THEOREM L: for every integer l >= (\d+):.*?<= ([\d.]+)/c\s+and.*?<= ([\d.]+)/c", txt)
check("Theorem L constants used here (l >= 138, K_* <= 92, K_neg <= 230) dominate large_l_proof_output.txt",
      m is not None and int(m.group(1)) <= L0_THM and Fr(m.group(2)) <= K_STAR and Fr(m.group(3)) <= K_NEG)


def alpha(n):
    return Fr((n - 2) * (n - 4), 2 ** n * fact(n - 1) ** 2)


def alphap(n):
    return Fr(3 * (n - 2) * (n - 3), 2 ** n * fact(n - 1) ** 2)


def e1(n):
    return Fr((n - 2) * (n * n - 4 * n - 3), 3)


def alpha1(n):
    return alphap(n) - alpha(n) * e1(n) if n >= 3 else Fr(0)


def tail_K(K, p, jshift, n1, rho):
    rho = Fr(rho)
    ratio = Fr(n1 + 1, n1) ** p * rho / (2 * n1 ** 2)
    assert ratio <= Fr(1, 2), ("ratio", ratio)
    n = n1 + 1
    return 2 * Fr(K * n ** p) * rho ** (n - jshift) / (2 ** n * Fr(fact(n - 1)) ** 2)


def omega_t(n, L1):
    cap = Fr(1, (L1 + 2) * (L1 - 1))
    Rn = Fr(1, 2) * (e1(n) + 4) ** 2 * (1 + 2 * cap)
    om = (abs(alpha(n)) + cap * alphap(n)) * Rn + alphap(n) * abs(e1(n))
    if n > L1 + 2:
        om = max(om, n ** 4 * abs(alpha(n)) + n ** 2 * abs(alpha1(n)))
    return om / 1000 if 'omega' in MUT else om


def xrange(coef, j, a, b, tailK, tailp):
    lo = hi = Fr(0)
    for n, v in coef.items():
        if n < j or v == 0:
            continue
        k = v * math.perm(n, j)
        ta, tb = k * a ** (n - j), k * b ** (n - j)
        lo += min(ta, tb)
        hi += max(ta, tb)
    t = tail_K(tailK, tailp + j, j, N1, b)
    return R.Iv(lo - t, hi + t)


lim = R.make_series(None, n1=N1)
A0 = {n: v for n, v in lim['a'].items()}                       # alpha_n (n <= N1)
A1 = {n: alpha1(n) for n in range(3, N1 + 1)}                  # alpha1_n


def NZ0(rho):
    return rho ** 3 / 8 - 3 * rho ** 2 / 4 + 9 * rho / 2 - 9


print("== 1. second-order coefficient bound |a_n - alpha_n - eps alpha1_n| <= eps^2 omega_n ==")
L_CHECK = 300
ok1 = True
worst = Fr(0)
for l in range(L_CHECK, L_CHECK + 40):
    a, eps = R.a_dict(l)
    cc = (l + 2) * (l - 1)
    for n in range(3, l + 13):
        an = a.get(n, Fr(0))
        diff = abs(an - alpha(n) - eps * alpha1(n))
        bound = eps ** 2 * omega_t(n, L_CHECK)
        if diff > bound:
            ok1 = False
        if bound > 0:
            worst = max(worst, diff / bound)
for n in range(3, 401):
    if omega_t(n, L_CHECK) * 2 ** n * fact(n - 1) ** 2 > n ** 8 and 'omega' not in MUT:
        ok1 = False
    if abs(alpha1(n)) * 2 ** n * fact(n - 1) ** 2 > n ** 5:
        ok1 = False
check(f"1: second-order coefficient bound for l = {L_CHECK}..{L_CHECK + 39}, 3 <= n <= l+12 (with omega built for L1 = {L_CHECK}); "
      "omega_n <= n^8/(2^n (n-1)!^2) and |alpha1_n| <= n^5/(2^n (n-1)!^2) for n <= 400", ok1)
print(f"   largest actual/bound ratio seen: {fmt(worst, 4)}")

print("\n== 2. algebra of the second-order expansion of Psi (exact, random rational points) ==")
import random
rnd = random.Random(5)


def Psi_l_exact(n, nu, ze, x, u, xi, e, rho):
    N = n + e * nu + e * e * ze
    X = x + e * u + e * e * xi
    return rho * rho * (2 + e * rho) * N * N - 4 * ((1 + 2 * e) * rho ** 3 / 8 - 3 * rho ** 2 / 4 + 9 * rho / 2 - 9) * X * X


def Psi1(n, nu, x, u, rho):
    return rho ** 3 * n * n + 4 * rho * rho * n * nu - 8 * NZ0(rho) * x * u - rho ** 3 * x * x


def Ehat(n, nu, ze, x, u, xi, e, rho):
    e2 = 2 * (nu * nu + 2 * n * ze) + 2 * rho * n * nu
    e3 = 4 * nu * ze + rho * (nu * nu + 2 * n * ze)
    e4 = 2 * ze * ze + 2 * rho * nu * ze
    e5 = rho * ze * ze
    f2 = NZ0(rho) * (u * u + 2 * x * xi) + rho ** 3 * x * u / 2
    f3 = 2 * NZ0(rho) * u * xi + rho ** 3 * (u * u + 2 * x * xi) / 4
    f4 = NZ0(rho) * xi * xi + rho ** 3 * u * xi / 2
    f5 = rho ** 3 * xi * xi / 4
    return rho * rho * (e2 + e * e3 + e * e * e4 + e ** 3 * e5) - 4 * (f2 + e * f3 + e * e * f4 + e ** 3 * f5)


okE = True
for _ in range(25):
    rr = Fr(rnd.randint(60, 110), 10)
    ee = Fr(rnd.randint(1, 40), 1000)
    n_, nu_, ze_, x_, u_, xi_ = [Fr(rnd.randint(-60, 60), rnd.randint(3, 13)) for _ in range(6)]
    lhs = Psi_l_exact(n_, nu_, ze_, x_, u_, xi_, ee, rr) - (2 * rr * rr * n_ * n_ - 4 * NZ0(rr) * x_ * x_) \
        - ee * Psi1(n_, nu_, x_, u_, rr)
    if lhs != ee * ee * Ehat(n_, nu_, ze_, x_, u_, xi_, ee, rr):
        okE = False
check("2: Psi_l - Psi_inf - eps Psi_1 = eps^2 E exactly (E, Psi_1 as coded), 25 random points", okE)

print("\n== 3. certified limit roots ==")
rho_star_lo, rho_star_hi = R.refine(lim, 'N', *R.certify_root(lim, 'N', Fr(8485700967, 10 ** 9), +1), +1, 80)
rho_neg_lo, rho_neg_hi = R.refine(lim, 'S', *R.certify_root(lim, 'S', Fr(9635989443, 10 ** 9), -1), -1, 80)
rho_star = (rho_star_lo + rho_star_hi) / 2
rho_neg = (rho_neg_lo + rho_neg_hi) / 2
print(f"   rho_* = {float(rho_star):.15f}   rho_neg = {float(rho_neg):.15f}")


def outward(x, up, den=10 ** 12):
    return Fr(math.ceil(x * den) if up else math.floor(x * den), den)


def attempt(L1):
    """All inequalities of Theorem M for l >= L1; returns (ok, dict of numbers)."""
    eps1 = Fr(1, (L1 + 2) * (L1 - 1))
    wN = K_STAR * eps1
    wP = K_NEG * eps1
    if 'window' in MUT:
        rho_star_c = rho_star + 3 * wN          # mutation: window not containing rho_*
    else:
        rho_star_c = rho_star
    WN = (outward(rho_star_c - wN, False), outward(rho_star_c + wN, True))
    WP = (outward(rho_neg - wP, False), outward(rho_neg + wP, True))
    out = {}
    # ---- N ----
    bN = WN[1]
    x = [xrange(A0, j, WN[0], WN[1], 3, 2) for j in range(3)]
    u = [xrange(A1, j, WN[0], WN[1], 1, 5) for j in range(2)]
    rhoI = R.Iv(*WN)
    nu = u[0] - (rhoI - 6) * u[1]
    if 'n1sign' in MUT:
        nu = -nu
    mu_Np = ((rhoI - 6) * x[2]).lo                       # -N_inf' >= mu_Np
    mu1N, nuN_max = nu.lo, nu.hi
    eta2N = sum(omega_t(n, L1) * bN ** (n - 1) * max(abs(6 * n - (n - 1) * WN[0]), abs(6 * n - (n - 1) * WN[1])) for n in range(3, N1 + 1)) \
        + tail_K(1 * (6 + bN), 9, 1, N1, bN)
    if 'eta2' in MUT:
        eta2N /= 1000
    out.update(mu1N=mu1N, nuN_max=nuN_max, mu_Np=mu_Np, eta2N=eta2N, eps1=eps1)
    condN = (mu1N > 0) and (mu_Np > 0) and ((L1 + 1) * mu1N / 2 > eta2N) and (eps1 * eta2N < mu1N)
    sN = Fr(105, 100) * (nuN_max + eps1 * eta2N) / mu_Np if (mu1N > 0 and mu_Np > 0) else Fr(10 ** 9)
    out['sN'] = sN
    tN = rho_star * Fr(L1 * L1 - 3 * L1 - 2, (L1 + 1) ** 2)
    condN = condN and sN < tN and sN < K_STAR and WN[0] <= rho_star and rho_star + sN * eps1 <= WN[1]
    out['tN'] = tN
    # ---- Psi ----
    bP = WP[1]
    xs = [xrange(A0, j, WP[0], WP[1], 3, 2) for j in range(3)]
    us = [xrange(A1, j, WP[0], WP[1], 1, 5) for j in range(2)]
    rhoP = R.Iv(*WP)
    nP = xs[0] - (rhoP - 6) * xs[1]
    nuP = us[0] - (rhoP - 6) * us[1]
    P1 = Psi1(nP, nuP, xs[0], us[0], rhoP)
    if 's1sign' in MUT:
        P1 = -P1
    mu1P = -P1.hi                                        # Psi_1 <= -mu1P
    P1abs = max(abs(P1.lo), abs(P1.hi))
    # eta2 for zeta on the Psi window and C2_0 for xi
    eta2Nz = sum(omega_t(n, L1) * bP ** (n - 1) * max(abs(6 * n - (n - 1) * WP[0]), abs(6 * n - (n - 1) * WP[1])) for n in range(3, N1 + 1)) \
        + tail_K(1 * (6 + bP), 9, 1, N1, bP)
    C20 = sum(omega_t(n, L1) * bP ** n for n in range(3, N1 + 1)) + tail_K(1, 8, 0, N1, bP)
    ze = R.Iv(-eta2Nz, eta2Nz)
    xi = R.Iv(-C20, C20)
    E = Ehat(nP, nuP, ze, xs[0], us[0], xi, R.Iv(0, eps1), rhoP)
    eta2P = max(abs(E.lo), abs(E.hi)) / (1000 if 'eta2' in MUT else 1)
    # Psi_inf' at rho_neg (point) and Psi_inf'' sup on window, as in large_l_proof.py
    pt = R.xder(lim, rho_neg, 3, 20)
    n0 = pt[0] - (R.Iv(rho_neg) - 6) * pt[1]
    n0p = -(R.Iv(rho_neg) - 6) * pt[2]
    rr_ = R.Iv(rho_neg)
    P0 = 4 * rr_ * n0 * n0 + 4 * rr_ * rr_ * n0 * n0p - 4 * (3 * rr_ ** 2 / 8 - 3 * rr_ / 2 + Fr(9, 2)) * pt[0] * pt[0] - 8 * NZ0(rr_) * pt[0] * pt[1]
    xx = [xrange(A0, j, WP[0], WP[1], 3, 2) for j in range(4)]
    n_ = xx[0] - (rhoP - 6) * xx[1]
    n1_ = -(rhoP - 6) * xx[2]
    n2_ = -xx[2] - (rhoP - 6) * xx[3]
    Pdd = (4 * n_ * n_ + 16 * rhoP * n_ * n1_ + 4 * rhoP * rhoP * (n1_ * n1_ + n_ * n2_) - 4 * (3 * rhoP / 4 - Fr(3, 2)) * xx[0] * xx[0]
           - 16 * (3 * rhoP ** 2 / 8 - 3 * rhoP / 2 + Fr(9, 2)) * xx[0] * xx[1] - 8 * NZ0(rhoP) * (xx[1] * xx[1] + xx[0] * xx[2]))
    Mdd = max(abs(Pdd.lo), abs(Pdd.hi))
    out.update(mu1P=mu1P, P1abs=P1abs, eta2P=eta2P, P0lo=P0.lo, Mdd=Mdd)
    condP = (mu1P > 0) and ((L1 + 1) * mu1P / 2 > eta2P) and (eps1 * eta2P < mu1P)
    # upper localisation: Psi_l(rho_neg + s eps) >= (P0lo - M s eps/2) s eps - eps P1abs - eps^2 eta2P > 0
    sP = Fr(101, 100) * (P1abs + eps1 * eta2P) / (P0.lo - Mdd * 20 * eps1 / 2) if P0.lo > 0 else Fr(10 ** 9)
    out['sP'] = sP
    tP = rho_neg * Fr(L1 * L1 - 3 * L1 - 2, (L1 + 1) ** 2)
    condP = condP and sP < tP and sP < 20 and sP < K_NEG and WP[0] <= rho_neg and rho_neg + sP * eps1 <= WP[1]
    out['tP'] = tP
    return condN and condP, out, (condN, condP)


print("\n== 4. smallest L1 (grid) for which every inequality of Theorem M holds ==")
print("      L1    mu1_N   eta2_N   (L1+1)mu1/2   s_N / t_N      mu1_P      eta2_P   (L1+1)mu1/2     s_P / t_P")
chosen = None
for L1 in [138, 200, 300, 400, 500, 600, 700, 800, 1000, 1200, 1500, 2000, 2500, 3000]:
    ok, o, (cN, cP) = attempt(L1)
    print(f"{L1:8d}  {fmt(o['mu1N'], 4):>7} {fmt(o['eta2N'], 4):>8} {fmt((L1 + 1) * o['mu1N'] / 2, 4):>12}   {fmt(o['sN'], 3)}/{fmt(o['tN'], 3)}"
          f"   {fmt(o['mu1P'], 4):>10} {fmt(o['eta2P'], 4):>10} {fmt((L1 + 1) * o['mu1P'] / 2, 4):>12}   {fmt(o['sP'], 3)}/{fmt(o['tP'], 3)}"
          f"   N:{'ok' if cN else 'no'} Psi:{'ok' if cP else 'no'}", flush=True)
    if ok and chosen is None:
        chosen = (L1, o)
        break
check("4: every inequality of Theorem M holds for some L1 <= 3000 (the range of the certified brackets)", chosen is not None)
if chosen is not None:
    L1, o = chosen
    print(f"\n   constants at L1 = {L1}: -N_inf' >= {fmt(o['mu_Np'], 5)}, nu in [{fmt(o['mu1N'], 5)}, {fmt(o['nuN_max'], 5)}], eta2_N = {fmt(o['eta2N'], 5)};"
          f" Psi_inf'(rho_neg) >= {fmt(o['P0lo'], 6)}, |Psi_inf''| <= {fmt(o['Mdd'], 4)}, |Psi_1| <= {fmt(o['P1abs'], 5)}, mu1_Psi = {fmt(o['mu1P'], 5)}, eta2_Psi = {fmt(o['eta2P'], 5)}")
    print(f"\n   THEOREM M holds for every l >= {L1}:  c r_*(l+1) < c r_*(l),  c r_neg(l+1) < c r_neg(l),  "
          f"0 < c r_* (l) - rho_* <= {fmt(o['sN'], 4)}/c,  0 < c r_neg(l) - rho_neg <= {fmt(o['sP'], 4)}/c,  and l^2 r_*, l^2 r_neg increase")

    print("\n== 5. spot checks of the expansion and of the key inequalities ==")
    eps1 = o['eps1']
    ok5 = ok5b = ok5c = ok5d = True
    wN = K_STAR * eps1
    wP = K_NEG * eps1
    worstN = worstP = Fr(0)
    sers = {l: R.make_series(l) for l in [L1, L1 + 1, L1 + 2, L1 + 3, 2 * L1, 2 * L1 + 1]}
    for l, ser in sers.items():
        eps = ser['eps']
        for rr in [rho_star + (k - 4) * wN / 4 for k in range(9)]:
            xv = R.xder(ser, rr, 1, 20)
            Nl = xv[0] - (rr - 6) * xv[1]
            x0 = R.xder(lim, rr, 1, 20)
            n0 = x0[0] - (rr - 6) * x0[1]
            u0 = [xrange(A1, j, rr, rr, 1, 5) for j in range(2)]
            nu0 = u0[0] - (rr - 6) * u0[1]
            rem = (Nl - n0 - eps * nu0) * (1 / (eps * eps))
            worstN = max(worstN, abs(rem.lo), abs(rem.hi))
            if max(abs(rem.lo), abs(rem.hi)) > o['eta2N']:
                ok5 = False
        for rr in [rho_neg + (k - 4) * wP / 4 for k in range(9)]:
            xv = R.xder(ser, rr, 1, 20)
            nl = xv[0] - (rr - 6) * xv[1]
            Pl = rr * rr * (2 + eps * rr) * nl * nl - 4 * ((1 + 2 * eps) * rr ** 3 / 8 - 3 * rr ** 2 / 4 + 9 * rr / 2 - 9) * xv[0] * xv[0]
            x0 = R.xder(lim, rr, 1, 20)
            n0 = x0[0] - (rr - 6) * x0[1]
            u0 = [xrange(A1, j, rr, rr, 1, 5) for j in range(2)]
            nu0 = u0[0] - (rr - 6) * u0[1]
            P1v = Psi1(n0, nu0, x0[0], u0[0], R.Iv(rr))
            Pinf = 2 * rr * rr * n0 * n0 - 4 * NZ0(rr) * x0[0] * x0[0]
            rem = (Pl - Pinf - eps * P1v) * (1 / (eps * eps))
            worstP = max(worstP, abs(rem.lo), abs(rem.hi))
            if max(abs(rem.lo), abs(rem.hi)) > o['eta2P']:
                ok5b = False
    check("5: |(N_l - N_inf - eps nu)/eps^2| <= eta2_N on the N window, l in {L1, L1+1, L1+2, 2 L1, 2 L1 + 1}", ok5)
    print(f"   largest seen {fmt(worstN, 4)} vs eta2_N = {fmt(o['eta2N'], 4)}")
    check("5: |(Psi_l - Psi_inf - eps Psi_1)/eps^2| <= eta2_Psi on the Psi window, same l", ok5b)
    print(f"   largest seen {fmt(worstP, 4)} vs eta2_Psi = {fmt(o['eta2P'], 4)}")
    # certified brackets: monotone step and localisation
    brk = {}
    for l, ser in sers.items():
        rn_, rs_ = R.locate(ser)
        brk[l] = (R.certify_root(ser, 'N', R.frac_from_float(rn_), +1), R.certify_root(ser, 'S', R.frac_from_float(rs_), -1))
    for l in [L1, L1 + 1]:
        (lo_n, hi_n), (lo_s, hi_s) = brk[l]
        cc = (l + 2) * (l - 1)
        if not (lo_n > rho_star and hi_n - rho_star <= o['sN'] / cc and lo_s > rho_neg and hi_s - rho_neg <= o['sP'] / cc):
            ok5c = False
    check("5: certified brackets satisfy the localisations 0 < rho_l - rho_inf <= s/c at l = L1, L1 + 1", ok5c)
    # the key inequality: N_{l+1} < 0 at the lower end of r_*(l)'s bracket and Psi_{l+1} > 0 at that of r_neg(l)'s bracket (exact)
    for l in [L1, L1 + 1, L1 + 2]:
        sl, s1 = sers[l], sers[l + 1]
        (lo_n, hi_n), (lo_s, hi_s) = brk[l]
        Nn = R.N_iv(s1, lo_n)
        eps_ = s1['eps']
        xv = R.xder(s1, lo_s, 1, 20)
        nl = xv[0] - (lo_s - 6) * xv[1]
        Pl = lo_s * lo_s * (2 + eps_ * lo_s) * nl * nl - 4 * ((1 + 2 * eps_) * lo_s ** 3 / 8 - 3 * lo_s ** 2 / 4 + 9 * lo_s / 2 - 9) * xv[0] * xv[0]
        if not (Nn.sign() == -1 and Pl.sign() == 1):
            ok5d = False
    check("5: N_{l+1}(lo of r_*(l)) < 0 and Psi_{l+1}(lo of r_neg(l)) > 0 for l = L1, L1+1, L1+2, so both roots move left", ok5d)

print(f"\n{time.time() - t0:.0f} s")
sys.exit(1 if failures else 0)

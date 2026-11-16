"""
Theorem L (#2): explicit constants l0, K_*, K_neg with, for every integer l >= l0 and c = (l+2)(l-1),
      |c r_*(l)   - rho_*|   <= K_*   / c,
      |c r_neg(l) - rho_neg| <= K_neg / c,
where rho_* = 8.48570096737... and rho_neg = 9.63598944258... are the simple roots of N_inf and Psi_inf (eps = 0 members).

Objects (rho = c r, eps = 1/c, X = c chi(rho/c), N = X - (rho-6) X'):
  Psi(rho) = rho^2 (2 + eps rho) N^2 - 4 N_Z X^2,   N_Z = (1+2 eps) rho^3/8 - 3 rho^2/4 + 9 rho/2 - 9 > 0 (rho >= 6).
  V^ = 2 c^4 (2 + eps rho) Psi / (rho^4 (rho-6)^2 X^2), so sign V^ = sign Psi on rho > 6 (this is 2 g^2 - V with
  g = -kappa_c = -c f N / ((rho-6) X) and the explicit Zerilli V; checked against chi_poly in vhat_structure.py).
Steps, each a check below:
  1. |a_n(eps) - alpha_n| <= eps beta_n for every n and every l >= LMIN (a_n = A_n(eps) P_n(eps), P_n in (0, 1+2eps],
     Weierstrass product inequality, and a_n = 0, c < n^2 for n >= l+3).
  2. N_l - N_inf = sum (a_n - alpha_n) rho^(n-1) (6n - (n-1) rho), bounded coefficientwise; X_l - X_inf likewise.
  3. (Psi_l - Psi_inf)/eps is an explicit polynomial in (x, n, d_X, d_N, eps, rho), enclosed by interval arithmetic.
  4. on a window around each limit root, N_inf' <= -mu_N and Psi_inf' >= mu_Psi (enclosures of entire series).
  5. an elementary lemma gives |root - limit root| <= eta/(mu c); uniqueness (T2, Theorem 3) makes the root found the
     root r_* / r_neg.
Controls: step 1 against the exact a_n; step 3's algebra against exact rational points; the final inequalities against
certified brackets of rescaled.py; Psi against S (same sign on the certified r_neg bracket).
Mutations: --mutate=beta | eta | dsign | poly | window.
"""
import argparse
import math
import os
import sys
import time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rescaled as R

ap = argparse.ArgumentParser()
ap.add_argument('--mutate', default='')
ap.add_argument('--lmin', type=int, default=10)
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
LMIN = args.lmin
EPS_CAP = Fr(1, (LMIN + 2) * (LMIN - 1))     # eps = 1/c(l) <= EPS_CAP for every l >= LMIN
N1 = 60                                       # terms summed exactly in the constants


def alpha(n):
    return Fr((n - 2) * (n - 4), 2 ** n * fact(n - 1) ** 2)


def alphap(n):
    return Fr(3 * (n - 2) * (n - 3), 2 ** n * fact(n - 1) ** 2)


def Sn(n):
    return Fr((n - 2) * (n * n - 4 * n - 3), 3) + 2


def beta_t(n):
    """beta~_n with |a_n(eps) - alpha_n| <= eps beta~_n for n >= 3, every l >= LMIN (a_n = 0 for n >= l+3 included)."""
    sig = max(Fr(2), Sn(n) * (1 + 2 * EPS_CAP))
    b = alphap(n) * (1 + 2 * EPS_CAP) + abs(alpha(n)) * sig
    out = b if n <= LMIN + 2 else max(b, n * n * abs(alpha(n)))     # n <= l+2 for every l >= LMIN when n <= LMIN+2
    return out / 1000 if 'beta' in MUT else out


def tail_K(K, p, jshift, n1, rho):
    """sum_{n>n1} K n^p rho^(n-jshift) / (2^n (n-1)!^2) <= 2 * first term, ratio <= 1/2 asserted."""
    rho = Fr(rho)
    ratio = Fr(n1 + 1, n1) ** p * rho / (2 * n1 ** 2)
    assert ratio <= Fr(1, 2), ("ratio", ratio)
    n = n1 + 1
    return 2 * Fr(K * n ** p) * rho ** (n - jshift) / (2 ** n * Fr(fact(n - 1)) ** 2)


print(f"== 1. coefficient bound |a_n(eps) - alpha_n| <= eps beta_n, l >= LMIN = {LMIN} ==")
ok1 = True
worst = Fr(0)
for n in range(3, 401):
    if beta_t(n) * 2 ** n * fact(n - 1) ** 2 > 3 * n ** 5 and 'beta' not in MUT:
        ok1 = False
for l in range(LMIN, 241):
    a, eps = R.a_dict(l)
    cc = (l + 2) * (l - 1)
    for n in range(3, l + 13):
        an = a.get(n, Fr(0))
        diff = abs(an - alpha(n))
        bound = eps * beta_t(n)
        if diff > bound:
            ok1 = False
        if bound > 0:
            worst = max(worst, diff / bound)
        if n >= l + 3 and not cc < n * n:
            ok1 = False
check(f"1: |a_n(1/c) - alpha_n| <= beta_n/c for l = {LMIN}..240, 3 <= n <= l+12; beta_n <= 3 n^5/(2^n (n-1)!^2) for n <= 400", ok1)
print(f"   largest actual/bound ratio seen: {fmt(worst, 4)}")

print("\n== 2. limit roots ==")
lim = R.make_series(None, n1=N1)
rho_star_lo, rho_star_hi = R.refine(lim, 'N', *R.certify_root(lim, 'N', Fr(8485700967, 10 ** 9), +1), +1, 80)
rho_neg_lo, rho_neg_hi = R.refine(lim, 'S', *R.certify_root(lim, 'S', Fr(9635989443, 10 ** 9), -1), -1, 80)
rho_star = (rho_star_lo + rho_star_hi) / 2
rho_neg = (rho_neg_lo + rho_neg_hi) / 2
print(f"   rho_*   = {float(rho_star):.15f}   rho_neg = {float(rho_neg):.15f}  (certified brackets of width {float(rho_star_hi - rho_star_lo):.1e})")

WN = (Fr(8), Fr(9)) if 'window' not in MUT else (Fr(9), Fr(10))    # window for rho_*; the mutation misses it
WS = (Fr(961, 100), Fr(966, 100))                                  # [9.61, 9.66] around rho_neg = 9.636


def xrange(coef, j, a, b, tailK, tailp):
    """Iv enclosing {F^(j)(rho) : a <= rho <= b} for F = sum coef[n] rho^n (each term monotone in rho > 0),
    plus a symmetric tail enclosure +-T from tail_K(tailK, tailp + j, j, N1, b)."""
    lo = Fr(0)
    hi = Fr(0)
    for n, v in coef.items():
        if n < j or v == 0:
            continue
        k = v * math.perm(n, j)
        ta, tb = k * a ** (n - j), k * b ** (n - j)
        lo += min(ta, tb)
        hi += max(ta, tb)
    t = tail_K(tailK, tailp + j, j, N1, b)
    return R.Iv(lo - t, hi + t)


def lim_range(j, a, b):
    return xrange(lim['a'], j, a, b, 3, 2)         # 0 <= alpha_n <= 3 n^2 /(2^n (n-1)!^2) for n >= 5


print("\n== 3. N window: mu_N and eta_N ==")


def eta_N_of(a, b):
    """sup over [a,b] of sum beta~_n rho^(n-1) |6n - (n-1) rho|  (bound of |N_l - N_inf|/eps) plus tail."""
    s = Fr(0)
    for n in range(3, N1 + 1):
        s += beta_t(n) * b ** (n - 1) * max(abs(6 * n - (n - 1) * a), abs(6 * n - (n - 1) * b))
    s += tail_K(3 * (6 + b), 6, 1, N1, b)          # beta~_n <= 3 n^5/(..), weight <= rho^(n-1) n (6 + rho)
    return s


def C0_of(b):
    s = sum(beta_t(n) * b ** n for n in range(3, N1 + 1))
    return s + tail_K(3, 5, 0, N1, b)


XN = [lim_range(j, WN[0], WN[1]) for j in range(3)]
rhoN = R.Iv(WN[0], WN[1])
negNp = (rhoN - 6) * XN[2]                       # -N_inf' = (rho-6) X''
mu_N = negNp.lo
eta_N = eta_N_of(*WN)
print(f"   -N_inf' >= {fmt(mu_N)} on [{fmt(WN[0])}, {fmt(WN[1])}],  eta_N = {fmt(eta_N)},  K_* = eta_N/mu_N = {fmt(eta_N / mu_N)}")
check("3: -N_inf' > 0 on the N window and rho_* inside it", mu_N > 0 and WN[0] < rho_star < WN[1])
dN = min(rho_star - WN[0], WN[1] - rho_star)
eps_N = mu_N * dN / eta_N
print(f"   eps_N = mu_N delta / eta_N = {fmt(eps_N)}")

print("\n== 4. Psi window: mu_Psi and eta_Psi ==")


def NZ0(rho):
    return rho ** 3 / 8 - 3 * rho ** 2 / 4 + 9 * rho / 2 - 9


def NZ0p(rho):
    return 3 * rho ** 2 / 8 - 3 * rho / 2 + Fr(9, 2)


def Psi_inf_prime(x, rho):
    """Psi_inf' with x = [X, X', X''] at rho (Fractions or Iv)."""
    n = x[0] - (rho - 6) * x[1]
    npr = -(rho - 6) * x[2]
    return 4 * rho * n * n + 4 * rho * rho * n * npr - 4 * NZ0p(rho) * x[0] * x[0] - 8 * NZ0(rho) * x[0] * x[1]


def Psi_inf_dd(x, rho):
    """Psi_inf'' with x = [X, X', X'', X'''] at rho (Fractions or Iv)."""
    n = x[0] - (rho - 6) * x[1]
    n1 = -(rho - 6) * x[2]
    n2 = -x[2] - (rho - 6) * x[3]
    return (4 * n * n + 16 * rho * n * n1 + 4 * rho * rho * (n1 * n1 + n * n2) - 4 * (3 * rho / 4 - Fr(3, 2)) * x[0] * x[0]
            - 16 * NZ0p(rho) * x[0] * x[1] - 8 * NZ0(rho) * (x[1] * x[1] + x[0] * x[2]))


XS = [lim_range(j, WS[0], WS[1]) for j in range(4)]
rhoS = R.Iv(WS[0], WS[1])
pt = lambda r_: R.xder(lim, r_, 3, 20)
P0 = Psi_inf_prime(pt(rho_neg), R.Iv(rho_neg))          # enclosure of Psi_inf'(rho_neg)
Pdd = Psi_inf_dd(XS, rhoS)                              # enclosure of Psi_inf'' on the window
Mdd = max(abs(Pdd.lo), abs(Pdd.hi))
wS = max(rho_neg - WS[0], WS[1] - rho_neg)
mu_P = P0.lo - Mdd * wS                                  # mean value: Psi_inf' >= Psi_inf'(rho_neg) - M |rho - rho_neg|
if 'dsign' in MUT:
    mu_P = -mu_P
print(f"   Psi_inf'(rho_neg) in [{float(P0.lo):.6g}, {float(P0.hi):.6g}],  |Psi_inf''| <= {fmt(Mdd)} on the window (half-width {fmt(wS)})")
print(f"   mean-value lower bound: Psi_inf' >= {fmt(mu_P)} on [{fmt(WS[0])}, {fmt(WS[1])}]")


def Psi_inf(x, rho):
    n = x[0] - (rho - 6) * x[1]
    return 2 * rho * rho * n * n - 4 * NZ0(rho) * x[0] * x[0]


# Psi_inf' closed form against a difference quotient at rho_neg (exact rationals)
h = Fr(1, 10 ** 5)
mid = lambda iv: (iv.lo + iv.hi) / 2
Pph, Pmh = Psi_inf(pt(rho_neg + h), R.Iv(rho_neg + h)), Psi_inf(pt(rho_neg - h), R.Iv(rho_neg - h))
dq = float((mid(Pph) - mid(Pmh)) / (2 * h))
form = float(mid(P0))
check("4: closed form of Psi_inf' agrees with a symmetric difference quotient of Psi_inf at rho_neg (relative 1e-6)",
      abs(form - dq) < 1e-6 * abs(dq))
ddq = float((mid(Psi_inf_prime(pt(rho_neg + h), R.Iv(rho_neg + h))) - mid(Psi_inf_prime(pt(rho_neg - h), R.Iv(rho_neg - h)))) / (2 * h))
Pdd0 = Psi_inf_dd(pt(rho_neg), R.Iv(rho_neg))
check("4: closed form of Psi_inf'' agrees with a difference quotient of Psi_inf' at rho_neg (relative 1e-6)",
      abs(float(mid(Pdd0)) - ddq) < 1e-6 * abs(ddq))
check("4: the mean-value lower bound of Psi_inf' is positive on the Psi window, rho_neg inside it", mu_P > 0 and WS[0] < rho_neg < WS[1])
Pz = Psi_inf(pt(rho_neg), R.Iv(rho_neg))
print(f"   Psi_inf(rho_neg) enclosure: [{float(Pz.lo):.3e}, {float(Pz.hi):.3e}] (rho_neg from the S bracket)")
check("4: Psi_inf vanishes at the S-root rho_neg (|Psi_inf| < 1e-8 there; Psi and S have the same root)", abs(Pz.lo) < 1e-8 and abs(Pz.hi) < 1e-8)


def DPsi(x, n, dX, dN, e, rho):
    """(Psi_l - Psi_inf)/eps with X = x + eps dX, N = n + eps dN; works for Fractions and Iv."""
    t1 = rho * rho * (4 * n * dN + 2 * e * dN * dN + rho * n * n + 2 * e * rho * n * dN + e * e * rho * dN * dN)
    t2 = 4 * NZ0(rho) * (2 * x * dX + e * dX * dX)
    t3 = rho ** 3 * (x + e * dX) * (x + e * dX)
    return t1 - t2 - (-t3 if 'poly' in MUT else t3)


import random
rnd = random.Random(11)
okT = True
for _ in range(20):
    rr = Fr(rnd.randint(60, 110), 10)
    ee = Fr(rnd.randint(1, 40), 1000)
    xx, nn_ = Fr(rnd.randint(-50, 50), 7), Fr(rnd.randint(-50, 50), 9)
    dX_, dN_ = Fr(rnd.randint(-50, 50), 11), Fr(rnd.randint(-50, 50), 13)
    Psi_l = rr * rr * (2 + ee * rr) * (nn_ + ee * dN_) ** 2 - 4 * ((1 + 2 * ee) * rr ** 3 / 8 - 3 * rr ** 2 / 4 + 9 * rr / 2 - 9) * (xx + ee * dX_) ** 2
    Psi_i = 2 * rr * rr * nn_ ** 2 - 4 * NZ0(rr) * xx ** 2
    if DPsi(xx, nn_, dX_, dN_, ee, rr) != (Psi_l - Psi_i) / ee:
        okT = False
check("4: DPsi reproduces (Psi_l - Psi_inf)/eps exactly at 20 random rational points", okT)

C0S = C0_of(WS[1])
etaN_S = eta_N_of(*WS)
eps_T = EPS_CAP
Pn = XS[0] - (rhoS - 6) * XS[1]
Tenc = DPsi(XS[0], Pn, R.Iv(-C0S, C0S), R.Iv(-etaN_S, etaN_S), R.Iv(0, eps_T), rhoS)
eta_P = max(abs(Tenc.lo), abs(Tenc.hi)) / (1000 if 'eta' in MUT else 1)
dS = min(rho_neg - WS[0], WS[1] - rho_neg)
eps_P = mu_P * dS / eta_P
print(f"   C_0 = {fmt(C0S)},  eta_N(Psi window) = {fmt(etaN_S)},  eta_Psi = {fmt(eta_P)},  K_neg = eta_Psi/mu_Psi = {fmt(eta_P / mu_P)},  eps_Psi = {fmt(eps_P)}")

eps0 = min(eps_N, eps_P, eps_T) / 2          # strict margin
need_c = max(1, math.ceil(1 / eps0))
l0 = max(LMIN, math.isqrt(need_c) - 3)
while (l0 + 2) * (l0 - 1) < need_c:
    l0 += 1
K_star = eta_N / mu_N
K_neg = eta_P / mu_P
print(f"\n   eps0 = {fmt(eps0)}  ->  l0 = {l0}  (c(l0) = {(l0 + 2) * (l0 - 1)}, LMIN = {LMIN})")
sepok = rho_neg - K_neg * eps0 > rho_star + K_star * eps0
check("5: the Psi window shrunk to rho_neg +- K_neg eps0 lies beyond rho_* + K_* eps0 (so it lies in (r_*, oo))", sepok)
print(f"   THEOREM L: for every integer l >= {l0}:  |c r_* - {float(rho_star):.10f}| <= {fmt(K_star)}/c   and   "
      f"|c r_neg - {float(rho_neg):.10f}| <= {fmt(K_neg)}/c")

print("\n== 5. the theorem against certified brackets ==")
ok5 = True
for l in sorted(set([l0, l0 + 1, l0 + 7, 2 * l0, 5 * l0, 1000, 2000, 3000])):
    ser = R.make_series(l)
    rn_, rs_ = R.locate(ser)
    lo_n, hi_n = R.certify_root(ser, 'N', R.frac_from_float(rn_), +1)
    lo_s, hi_s = R.certify_root(ser, 'S', R.frac_from_float(rs_), -1)
    cc = (l + 2) * (l - 1)
    okn = max(abs(lo_n - rho_star), abs(hi_n - rho_star)) <= K_star / cc
    oks = max(abs(lo_s - rho_neg), abs(hi_s - rho_neg)) <= K_neg / cc
    print(f"   l = {l:5d}:  c r_* - rho_* = {float(lo_n - rho_star):+.3e} (bound {float(K_star / cc):.3e})   "
          f"c r_neg - rho_neg = {float(lo_s - rho_neg):+.3e} (bound {float(K_neg / cc):.3e})")
    ok5 = ok5 and okn and oks
check("5: certified brackets of r_*, r_neg satisfy the Theorem L bounds at the l listed", ok5)

print("\n== 6. spot checks of the perturbation constants and of Psi against S ==")
ok6 = ok6b = ok6c = True
worstP = worstN = Fr(0)
for l in [l0, l0 + 3, 2 * l0, 1000]:
    ser = R.make_series(l)
    eps = ser['eps']
    for rr in [WS[0] + (WS[1] - WS[0]) * k / 8 for k in range(9)]:
        xl = R.xder(ser, rr, 2, 20)
        xi = R.xder(lim, rr, 2, 20)
        nl = xl[0] - (rr - 6) * xl[1]
        Pl = rr * rr * (2 + eps * rr) * nl * nl - 4 * ((1 + 2 * eps) * rr ** 3 / 8 - 3 * rr ** 2 / 4 + 9 * rr / 2 - 9) * xl[0] * xl[0]
        Pi = Psi_inf(xi, R.Iv(rr))
        dl = max(abs(Pl.lo - Pi.hi), abs(Pl.hi - Pi.lo))
        worstP = max(worstP, dl / eps)
        if dl / eps > eta_P:
            ok6 = False
    for rr in [WN[0] + (WN[1] - WN[0]) * k / 8 for k in range(9)]:
        Nl, Ni = R.N_iv(ser, rr), R.N_iv(lim, rr)
        dl = max(abs(Nl.lo - Ni.hi), abs(Nl.hi - Ni.lo))
        worstN = max(worstN, dl / eps)
        if dl / eps > eta_N:
            ok6b = False
    rn_, rs_ = R.locate(ser)
    lo_s, hi_s = R.certify_root(ser, 'S', R.frac_from_float(rs_), -1)

    def Psi_l(xv, rr):
        nl = xv[0] - (rr - 6) * xv[1]
        return rr * rr * (2 + eps * rr) * nl * nl - 4 * ((1 + 2 * eps) * rr ** 3 / 8 - 3 * rr ** 2 / 4 + 9 * rr / 2 - 9) * xv[0] * xv[0]
    if not (Psi_l(R.xder(ser, lo_s, 2, 20), lo_s).sign() == -1 and Psi_l(R.xder(ser, hi_s, 2, 20), hi_s).sign() == 1):
        ok6c = False
check("6: |Psi_l - Psi_inf|/eps <= eta_Psi at 9 points of the Psi window for l in {l0, l0+3, 2 l0, 1000}", ok6)
print(f"   largest |Psi_l - Psi_inf|/eps seen: {fmt(worstP, 4)}  vs eta_Psi = {fmt(eta_P, 4)}")
check("6: |N_l - N_inf|/eps <= eta_N at 9 points of the N window for the same l", ok6b)
print(f"   largest |N_l - N_inf|/eps seen: {fmt(worstN, 4)}  vs eta_N = {fmt(eta_N, 4)}")
check("6: Psi_l < 0 at the lower end and > 0 at the upper end of the certified S bracket of r_neg (Psi and S share the root)", ok6c)

print(f"\n{time.time() - t0:.0f} s")
sys.exit(1 if failures else 0)

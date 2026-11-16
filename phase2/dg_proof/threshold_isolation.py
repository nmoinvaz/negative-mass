"""
Certified isolation of c r_* and c r_neg for l = 2..LMAX (#2), by exact rational arithmetic.

  c r_*   : the unique root of N(rho) = X - (rho-6) X' on (6, oo)            (T2 of RESULTS.md, all l)
  c r_neg : the unique zero of V^ on (r_*, oo), where sign V^ = sign S(rho)    (Theorem 3, vhat_structure.py)
For each l a bracket [lo, hi] is certified with sign F(lo) = lo_sign, sign F(hi) = -lo_sign, F evaluated as an exact
rational interval (truncated series at 90 terms plus a rigorous tail enclosure, see rescaled.py).  Uniqueness of the
root is the theorem, so the bracket contains THE root.  No floating point enters a certificate.

Checks (each prints PASS/FAIL, the exit status is non-zero if any fails):
  P  planted truth: X = rho^2 has N = rho(12 - rho) and S = 6 rho^3 (rho - 16), so the certified brackets must contain 12 and 16
  A  the rescaled coefficient formula equals common.chi_poly (exact), l = 2..30
  B  the eps = 0 roots agree with the Bessel closed form X = s[(s+3)I0(2 sqrt s) - 4 sqrt(s) I1(2 sqrt s)], s = rho/2
  C  independent root isolation (sympy Sturm on chi_poly and on the numerator of phi0_xx), l = 2..24
  D  agreement with the stored large_l_limit_output.txt table, l = 2, 3, 5, 10, 20, 40, 80
  E  r_* < r_neg, strict monotone decrease of c r_* and c r_neg in l, approach from above to the eps = 0 roots,
     l^2 r_* and l^2 r_neg increasing from some l and below the limit, for every l = 2..LMAX (interval separation)
  F  a sparse set of larger l
Mutations: --mutate=coef | sterm | ntail | mono | planted.
"""
import argparse
import os
import sys
import time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'dg_spectrum'))
sys.path.insert(0, HERE)
import rescaled as R

ap = argparse.ArgumentParser()
ap.add_argument('--lmax', type=int, default=3000)
ap.add_argument('--mutate', default='')
args = ap.parse_args()
R.MUT = set(m for m in args.mutate.split(',') if m)
MONO_MUT = 'mono' in R.MUT

failures = []


def check(name, ok):
    print(("PASS  " if ok else "FAIL  ") + name, flush=True)
    if not ok:
        failures.append(name)


def fmt(x, nd=20):
    """decimal string of a Fraction to nd places (truncated toward -inf)."""
    x = Fr(x)
    neg = x < 0
    x = abs(x)
    ip = x.numerator // x.denominator
    fp = (x - ip) * 10 ** nd
    return ("-" if neg else "") + f"{ip}." + str(fp.numerator // fp.denominator).zfill(nd)


def isolate(l, refine_steps=0):
    ser = R.make_series(l)
    rn, rsg = R.locate(ser)
    lo_n, hi_n = R.certify_root(ser, 'N', R.frac_from_float(rn), +1)
    lo_s, hi_s = R.certify_root(ser, 'S', R.frac_from_float(rsg), -1)
    if refine_steps:
        lo_n, hi_n = R.refine(ser, 'N', lo_n, hi_n, +1, refine_steps)
        lo_s, hi_s = R.refine(ser, 'S', lo_s, hi_s, -1, refine_steps)
    return dict(lo_n=lo_n, hi_n=hi_n, lo_s=lo_s, hi_s=hi_s)


t0 = time.time()
print("== P: planted truth through the same certification code (X = rho^2, eps = 0) ==")
planted = dict(a={2: Fr(1), 3: (Fr(1, 1000) if 'planted' in R.MUT else Fr(0))}, eps=Fr(0), n1=90, has_tail=False, l=None)
try:
    lo12, hi12 = R.certify_root(planted, 'N', Fr(12), +1)
    okP1 = lo12 <= 12 <= hi12
except RuntimeError:
    okP1 = False
try:
    lo16, hi16 = R.certify_root(planted, 'S', Fr(16), -1)
    okP2 = lo16 <= 16 <= hi16
except RuntimeError:
    okP2 = False
check("P: the certified bracket of the root of N for X = rho^2 contains 12", okP1)
check("P: the certified bracket of the root of S for X = rho^2 contains 16", okP2)

print("\n== A: rescaled coefficients against common.chi_poly (exact) ==")
import sympy as sp
from common import chi_poly, rs_
okA = True
for l in range(2, 31):
    a, eps = R.a_dict(l)
    c = (l + 2) * (l - 1)
    chi = sp.Poly(chi_poly(l), rs_)
    for n in range(1, l + 3):
        want = chi.coeff_monomial(rs_ ** n) * c / sp.Integer(c) ** n
        got = sp.Rational(a[n].numerator, a[n].denominator)
        if want != got:
            okA = False
check("A: a_n(1/c) = c * [r^n]chi / c^n for l = 2..30, every n", okA)

print("\n== B: eps = 0 roots against the Bessel closed form ==")
lim = isolate(None, refine_steps=60)
import mpmath as mp
mp.mp.dps = 40


def Xb(rho):
    s = rho / 2
    z = 2 * mp.sqrt(s)
    return s * ((s + 3) * mp.besseli(0, z) - 4 * mp.sqrt(s) * mp.besseli(1, z))


def Nb(rho):
    return Xb(rho) - (rho - 6) * mp.diff(Xb, rho)


def Sb(rho):
    x0, x1, x2 = Xb(rho), mp.diff(Xb, rho), mp.diff(Xb, rho, 2)
    n = x0 - (rho - 6) * x1
    return -2 * (rho * (rho - 6) * x0 * x2 + n * (x0 + 2 * rho * x1))


rho_star_b = mp.findroot(Nb, mp.mpf('8.4857'))
rho_neg_b = mp.findroot(Sb, mp.mpf('9.636'))
rho_star_lo, rho_star_hi = lim['lo_n'], lim['hi_n']
rho_neg_lo, rho_neg_hi = lim['lo_s'], lim['hi_s']
print(f"  rho_*   in [{fmt(rho_star_lo, 24)}, {fmt(rho_star_hi, 24)}]   Bessel {mp.nstr(rho_star_b, 28)}")
print(f"  rho_neg in [{fmt(rho_neg_lo, 24)}, {fmt(rho_neg_hi, 24)}]   Bessel {mp.nstr(rho_neg_b, 28)}")
tolB = mp.mpf(10) ** (-18)
check("B: certified eps = 0 bracket of rho_* contains the Bessel-form root (to 1e-18)",
      mp.mpf(rho_star_lo.numerator) / rho_star_lo.denominator - tolB <= rho_star_b <= mp.mpf(rho_star_hi.numerator) / rho_star_hi.denominator + tolB)
check("B: certified eps = 0 bracket of rho_neg contains the Bessel-form root (to 1e-18)",
      mp.mpf(rho_neg_lo.numerator) / rho_neg_lo.denominator - tolB <= rho_neg_b <= mp.mpf(rho_neg_hi.numerator) / rho_neg_hi.denominator + tolB)

print("\n== C: independent isolation by sympy Sturm sequences, l = 2..24 ==")
r = rs_
okC = True
for l in range(2, 25):
    c = (l + 2) * (l - 1)
    chi = chi_poly(l)
    Npoly = sp.Poly(sp.expand(c * chi - (c * r - 6) * sp.diff(chi, r)), r)
    isoN = [iv for iv in Npoly.intervals(eps=sp.Rational(1, 10 ** 15)) if iv[0][0] >= sp.Rational(6, c)]
    phi0 = (c * r - 6) / chi
    f = 1 + 2 / r
    num = sp.numer(sp.together(sp.diff(f * sp.diff(phi0, r), r)))
    Pn = sp.Poly(sp.expand(num), r)
    isoS = [iv for iv in Pn.intervals(eps=sp.Rational(1, 10 ** 15)) if iv[0][0] > 0]
    mine = isolate(l, refine_steps=40)
    # sympy: exactly one root of N beyond r_s; largest root of the phi0_xx numerator is r_neg
    a_n, b_n = isoN[0][0]
    a_s, b_s = isoS[-1][0]
    tosp = lambda q_: sp.Rational(q_.numerator, q_.denominator)
    okN = len(isoN) == 1 and a_n <= tosp(mine['lo_n']) / c and tosp(mine['hi_n']) / c <= b_n
    okS = a_s <= tosp(mine['lo_s']) / c and tosp(mine['hi_s']) / c <= b_s
    if not (okN and okS):
        okC = False
        print(f"  l = {l}: N {okN}  S {okS}", flush=True)
check("C: certified brackets lie inside sympy's isolating intervals (N from chi_poly, S as phi0_xx numerator), l = 2..24", okC)

print("\n== D: stored large_l_limit_output.txt ==")
stored = {}
with open(os.path.join(HERE, 'large_l_limit_output.txt')) as fh:
    for line in fh:
        parts = line.split()
        if len(parts) >= 4 and parts[0].isdigit():
            stored[int(parts[0])] = (float(parts[1]), float(parts[2]))
okD = True
for l in sorted(stored):
    got = isolate(l)
    dn = abs(float(got['lo_n']) - stored[l][0])
    ds = abs(float(got['lo_s']) - stored[l][1])
    if dn > 6e-6 or ds > 6e-6:
        okD = False
        print(f"  l = {l}: {dn:.2e} {ds:.2e}")
check(f"D: c r_* and c r_neg agree with the stored table to its 5 printed decimals, l in {sorted(stored)}", okD)

print(f"\n== E: l = 2..{args.lmax}, certified brackets (c r_* [lo, hi], c r_neg [lo, hi]) ==")
rows = {}
for l in range(2, args.lmax + 1):
    rows[l] = isolate(l)
    if l % 500 == 0:
        print(f"   ... l = {l} done, {time.time() - t0:.0f}s", flush=True)
sel = [2, 3, 5, 10, 20, 40, 80, 81, 100, 150, 200, 300, 500, 1000, 2000, 3000]
print("   l     c r_*                 c r_neg                l^2 r_*       l^2 r_neg   (brackets certified, printed lo)")
for l in sel:
    if l in rows:
        row = rows[l]
        c = (l + 2) * (l - 1)
        sc = Fr(l * l, c)
        print(f"{l:5d}  {fmt(row['lo_n'], 12)}  {fmt(row['lo_s'], 12)}  {fmt(sc * row['lo_n'], 8)}  {fmt(sc * row['lo_s'], 8)}")

ok_order = all(rows[l]['hi_n'] < rows[l]['lo_s'] for l in rows)
check(f"E1: r_* < r_neg for every l = 2..{args.lmax} (hi of r_* below lo of r_neg)", ok_order)
def dec(key_lo, key_hi):
    ok = True
    for l in range(2, args.lmax):
        a, b = rows[l], rows[l + 1]
        cond = (b[key_hi] < a[key_lo]) if not MONO_MUT else (b[key_lo] > a[key_hi])
        if not cond:
            ok = False
            print(f"   non-separation at l = {l} for {key_lo}")
            break
    return ok
check(f"E2: c r_* strictly decreasing for l = 2..{args.lmax}", dec('lo_n', 'hi_n'))
check(f"E3: c r_neg strictly decreasing for l = 2..{args.lmax}", dec('lo_s', 'hi_s'))
check("E4: c r_* > rho_* (limit bracket) for every l, c r_neg > rho_neg for every l (approach from above)",
      all(rows[l]['lo_n'] > rho_star_hi and rows[l]['lo_s'] > rho_neg_hi for l in rows))

def l2_increasing_from(key_lo, key_hi):
    """smallest l1 such that l^2 r(l+1) > l^2 r(l) for all l1 <= l < lmax (interval separation), None if none."""
    good = None
    for l in range(args.lmax - 1, 1, -1):
        c, c1 = (l + 2) * (l - 1), (l + 3) * l
        up = Fr((l + 1) ** 2, c1) * rows[l + 1][key_lo] > Fr(l * l, c) * rows[l][key_hi]
        if up:
            good = l
        else:
            break
    return good
l1n = l2_increasing_from('lo_n', 'hi_n')
l1s = l2_increasing_from('lo_s', 'hi_s')
print(f"   l^2 r_* increasing for l >= {l1n} (checked to {args.lmax}); l^2 r_neg increasing for l >= {l1s}")
check("E5: l^2 r_* and l^2 r_neg are strictly increasing in l from l = 5 on, up to LMAX (both have their minimum at l = 5)",
      l1n is not None and l1s is not None and l1n <= 5 and l1s <= 5)
check("E6: l^2 r_* < rho_* and l^2 r_neg < rho_neg for every l = 3..LMAX (approach from below)",
      all(Fr(l * l, (l + 2) * (l - 1)) * rows[l]['hi_n'] < rho_star_lo and Fr(l * l, (l + 2) * (l - 1)) * rows[l]['hi_s'] < rho_neg_lo
          for l in range(3, args.lmax + 1)))
dmax = max((rows[l]['hi_n'] - rows[l]['lo_n'] for l in rows), default=0)
print(f"   widest bracket in the table: {float(dmax):.3g}")

print("\n== F: sparse larger l (brackets refined by exact bisection to 1e-22) ==")
print("      l          c r_* - rho_*            (c r_* - rho_*) c      c r_neg - rho_neg       (c r_neg - rho_neg) c")
prev = None
okF = True
for l in [4000, 10 ** 4, 3 * 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]:
    row = isolate(l, refine_steps=40)
    c = (l + 2) * (l - 1)
    dn_lo, dn_hi = row['lo_n'] - rho_star_hi, row['hi_n'] - rho_star_lo
    ds_lo, ds_hi = row['lo_s'] - rho_neg_hi, row['hi_s'] - rho_neg_lo
    print(f"{l:9d}   {float(dn_lo):.6e}    {float(dn_lo) * c:12.6f}      {float(ds_lo):.6e}     {float(ds_lo) * c:12.6f}")
    if dn_lo <= 0 or ds_lo <= 0:
        okF = False
check("F: c r_* > rho_* and c r_neg > rho_neg also at l = 4000 .. 10^7", okF)

print(f"\n{time.time() - t0:.0f} s")
sys.exit(1 if failures else 0)

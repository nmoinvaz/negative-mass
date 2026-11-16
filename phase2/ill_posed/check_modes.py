"""
Checks of the exact-mode library (modes.py) before it is used for anything else.  Each check can fail.

  A  chi coefficients against dg_spectrum/common.chi_poly (sympy, gamma form), l = 2..8, exact.
  B  psi0 = chi/(c r - 6) solves the static Zerilli equation psi0_xx = V psi0 (V written out, M = -1), l = 2..12, at 6 points.
  C  Vhat = g^2 - g_x equals 2 g^2 - V_Zerilli, l = 2..12.
  D  closed form Psihat = 2^k e^{-omega x} r^4 Pi/chi equals (-d/dx + g) Psi_DG evaluated directly, l = 2..12.
  E  H Psihat = -omega^2 Psihat (residual of the ODE, mp, 6 points incl. next to r_s), l = 2..12, omega = k/2 exact.
  F  Pi and chi/r have no positive root (exact Sturm count), l = 2..40: Psihat is nodeless, chi > 0.
  G  Psihat / r^3 tends to a finite non-zero limit at r -> 0.
  H  ||Psihat||_{L^2(dx)} by two different quadratures (in r and in x).
Mutations: --mutate=chi (c_3 -> 1.01 c_3) | omega (omega -> omega + 1) | jac (drop dx/dr in the norm) | vhat (Vhat * (1 + 1e-6)) | nodes (Pi * (r-1))
  | pow (r^4 -> r^5 in G).
"""
import argparse
import os
import sys

import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'dg_spectrum'))
import modes as M
from chk import Checker

ap = argparse.ArgumentParser()
ap.add_argument('--mutate', default='')
args = ap.parse_args()
if args.mutate in ('chi', 'omega', 'jac', 'vhat'):
    M.MUT.add(args.mutate)
ck = Checker()
mp.mp.dps = 40
LS = range(2, 13)


def sample_points(l):
    rs = mp.mpf(6) / M.cc(l)
    return [rs * s for s in ('0.3', '0.9', '1.001', '1.7', '4', '25')] if False else \
           [rs * mp.mpf(s) for s in ('0.3', '0.9', '1.001', '1.7', '4', '25')]


try:
    # A
    from common import chi_poly, rs_
    ok = True
    for l in range(2, 9):
        ref = sp.Poly(chi_poly(l), rs_)
        mine = M.chi_coeffs(l)
        for n in range(l + 3):
            if sp.Rational(mine[n].numerator, mine[n].denominator) != ref.coeff_monomial(rs_**n):
                ok = False
    ck("A chi coefficients = common.chi_poly, l = 2..8 (exact)", ok)

    # B, C
    worstB = worstC = mp.mpf(0)
    for l in LS:
        c = M.cc(l)
        for r in sample_points(l):
            psi0 = lambda s: M.chi(l, s, 0)[0] / (c * s - 6)
            d1 = lambda s: M.f_of(s) * mp.diff(psi0, s)
            pxx = M.f_of(r) * mp.diff(d1, r)
            worstB = max(worstB, abs(pxx - M.Vzer(l, r) * psi0(r)) / abs(M.Vzer(l, r) * psi0(r)))
            a, b = M.Vhat(l, r), M.Vhat2(l, r)
            worstC = max(worstC, abs(a - b) / (1 + abs(a)))
    ck("B psi0_xx = V_Zerilli psi0, l = 2..12", worstB < 1e-20, f"worst rel {mp.nstr(worstB, 3)}")
    ck("C g^2 - g_x = 2 g^2 - V_Zerilli, l = 2..12", worstC < 1e-20, f"worst rel {mp.nstr(worstC, 3)}")

    # D
    worstD = mp.mpf(0)
    for l in LS:
        for r in sample_points(l):
            ratio = M.psihat_direct(l, r) / (2**M.kk(l) * M.psihat_unnorm(l, r))
            worstD = max(worstD, abs(ratio - 1))
    ck("D (-d/dx + g) Psi_DG = 2^k e^{-omega x} r^4 Pi/chi, l = 2..12", worstD < 1e-15, f"worst |ratio-1| {mp.nstr(worstD, 3)}")

    # E
    worstE = mp.mpf(0)
    for l in LS:
        w2 = float(M.omega(l))**2
        for r in sample_points(l):
            worstE = max(worstE, abs(M.residual(l, r)) / w2)
    ck("E (H + omega^2) Psihat = 0, l = 2..12", worstE < 1e-15, f"worst |res|/omega^2 {mp.nstr(worstE, 3)}")

    # F
    x = rs_
    bad = []
    for l in range(2, 41):
        pi = sp.Poly(sum(sp.Rational(a.numerator, a.denominator) * x**i for i, a in enumerate(M.pi_coeffs(l))), x)
        if args.mutate == 'nodes':
            pi = pi * sp.Poly(x - 1, x)
        ch = sp.Poly(sum(sp.Rational(a.numerator, a.denominator) * x**(i - 1) for i, a in enumerate(M.chi_coeffs(l)) if i), x)
        if pi.count_roots(0, None) != 0 or ch.count_roots(0, None) != 0:
            bad.append(l)
    ck("F Pi and chi/r have no root in [0, oo), l = 2..40 (exact Sturm)", not bad, f"bad l: {bad}")

    # G
    okG = True
    for l in LS:
        extra = (lambda r: r) if args.mutate == 'pow' else (lambda r: 1)
        a = M.psihat_unnorm(l, mp.mpf('1e-9')) * extra(mp.mpf('1e-9')) / mp.mpf('1e-9')**3
        b = M.psihat_unnorm(l, mp.mpf('1e-10')) * extra(mp.mpf('1e-10')) / mp.mpf('1e-10')**3
        if not (abs(a / b - 1) < 1e-5 and abs(a) > 1e-30):
            okG = False
    ck("G Psihat ~ r^3 at 0 (Psihat/r^3 at 1e-9, 1e-10 agree to 1e-5, and are not 0), l = 2..12", okG)

    # H
    worstH = mp.mpf(0)
    for l in (2, 3, 5, 8):
        n1 = M.norm_psihat(l)
        w = float(M.omega(l))
        xs = [0] + [s / w for s in (0.5, 1, 2, 4, 8, 16, 32, 64, 128)] + [mp.inf]
        n2 = mp.sqrt(mp.quad(lambda xx: M.psihat_unnorm(l, M.r_of_x(xx))**2, xs))
        worstH = max(worstH, abs(n1 / n2 - 1))
    ck("H L^2(dx) norm, quadrature in r vs in x, l = 2, 3, 5, 8", worstH < 1e-10, f"worst {mp.nstr(worstH, 3)}")
except AssertionError as e:
    ck("exact identities inside modes.py (assertions)", False, str(e))

ck.done()

"""
Planted truth for the time-domain evolver (evolve.py), run before any real data.  Controls and pass criteria, written first:

  P1  ill-posed toy.  V = 3/(4x^2) - 3 w/x has the exact eigenpair psi = x^{3/2} e^{-w x}, E = -w^2 (the same x^{3/2} branch at the
      singular end as the real problem).  Data (psi, 0) must grow as cosh(w t): |ln G - ln cosh(w t)| <= TOL(w t) at every recorded
      time, for w_n = 5 * 2^n, n = 0..5.  TOL(y) = 0.01 + 5e-4 y.  (Calibrated on the convergence run below: with 40 points per 1/w the
      error is 2e-4 y; second-order convergence is checked in P1b.)
  P1b convergence: halving h divides the error at w t = 20 by at least 3.
  P2  well-posed toy.  V = 3/(4x^2) >= 0.  For every w_n, G(t) = ||u(t)||/||u(0)|| <= 1 + 1e-3 on 0 <= t <= 20/w (v0 = 0, H >= 0).
  P3  Robin end.  V = 0 on x >= x_R = 0, u_x = kappa u.  kappa = -w: data e^{-w x} grows as cosh(w t) (within TOL).  kappa = +w:
      G <= 1 + 1e-3.  This is the sign convention used for the truncated star.
  P4  sequence diagnosis at fixed t0 = 0.1.  For the P1 family the data norm eps_n = 1/cosh(w_n t0) decreases strictly to < 1e-6 and
      the solution norm eps_n G_n(t0) stays within [0.9, 1.1]; for the P2 family G_n(t0) <= 1.001 for all n.  The ill-posed family is
      told apart from the well-posed one by sup_n G_n(t0) (> 1e3 vs <= 1.001).
Mutations: --mutate=omega (predict with 1.05 w) | signV (V -> -V in P1) | label (P2 uses the P1 potential) | kappa (P3 kappa -> -kappa)
  | dt (leapfrog step 2.2 x the stability limit).
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from chk import Checker
from evolve import Op, evolve, ln_cosh

ap = argparse.ArgumentParser()
ap.add_argument('--mutate', default='')
args = ap.parse_args()
MU = args.mutate
ck = Checker()
NPW = 40
TOL = lambda y: 0.01 + 5e-4 * np.asarray(y)
DTF = 2.2 if MU == 'dt' else 0.5


def grid(w, T, x0=0.0, npw=NPW):
    h = 1.0 / (npw * w)
    n = int(round((30.0 / w + T) / h))
    return x0 + h * np.arange(n + 1)


def pot_coulomb(x, w, sign=1.0):
    V = np.zeros_like(x)
    V[1:] = 3.0 / (4 * x[1:]**2) - sign * 3 * w / x[1:]
    return V


def run_dirichlet(w, T, V_of, npw=NPW, nrec=4):
    x = grid(w, T, npw=npw)
    op = Op(x, V_of(x, w), 'dirichlet')
    u0 = op.xs**1.5 * np.exp(-w * op.xs)
    t, ln, _, _ = evolve(op, u0, None, T=T, nrec=nrec, dt_frac=DTF)
    return t, ln


wpred = lambda w: 1.05 * w if MU == 'omega' else w
Vill = lambda x, w: pot_coulomb(x, w, -1.0 if MU == 'signV' else 1.0)
Vwell = (lambda x, w: pot_coulomb(x, w)) if MU == 'label' else (lambda x, w: pot_coulomb(x, 0.0))

WS = [5.0 * 2**n for n in range(6)]
print("== P1  ill-posed toy V = 3/(4x^2) - 3w/x, data (x^1.5 e^{-wx}, 0), npw =", NPW)
print("   w      w t    ln G_num   ln cosh(w t)   diff")
worst = 0.0
okP1 = True
for w in WS:
    T = 20.0 / w
    t, ln = run_dirichlet(w, T, Vill)
    ex = ln_cosh(wpred(w) * t)
    for ti, a, b in zip(t, ln, ex):
        print(f"{w:6.1f} {w * ti:7.2f} {a:10.4f} {b:12.4f} {a - b:+9.4f}")
        okP1 &= abs(a - b) <= TOL(w * ti)
        worst = max(worst, abs(a - b) / TOL(w * ti))
ck("P1 growth of the exact toy mode = cosh(w t), w = 5..160, 4 times each", okP1, f"worst |diff|/TOL = {worst:.2f}")

print("\n== P1b convergence at w = 1, w t = 20")
errs = []
for npw in (20, 40, 80):
    t, ln = run_dirichlet(1.0, 20.0, Vill, npw=npw, nrec=1)
    errs.append(abs(ln[-1] - ln_cosh(wpred(1.0) * t[-1])))
    print(f"   npw = {npw:3d}   |ln G - ln cosh| = {errs[-1]:.3e}")
ck("P1b second-order convergence (error ratio >= 3 per halving of h)", errs[0] / errs[1] >= 3 and errs[1] / errs[2] >= 3,
   f"ratios {errs[0] / errs[1]:.2f}, {errs[1] / errs[2]:.2f}")

print("\n== P2  well-posed toy V = 3/(4x^2)")
mx = -1e9
for w in WS:
    t, ln = run_dirichlet(w, 20.0 / w, Vwell, nrec=8)
    mx = max(mx, ln.max())
    print(f"   w = {w:6.1f}   max ln G = {ln.max():+.2e}")
ck("P2 G <= 1 + 1e-3 for every w", mx <= np.log(1.001), f"max ln G = {mx:+.2e}")

print("\n== P3  Robin end, V = 0, x >= 0, u_x = kappa u, data e^{-w x}")
okP3a = True
okP3b = True
for w in WS[:4]:
    T = 20.0 / w
    h = 1.0 / (NPW * w)
    n = int(round((30.0 / w + T) / h))
    x = h * np.arange(n + 1)
    for sgn in (-1.0, 1.0):
        kap = sgn * w * (-1 if MU == 'kappa' else 1)
        op = Op(x, np.zeros_like(x), 'robin', kap)
        u0 = np.exp(-w * op.xs)
        t, ln, _, _ = evolve(op, u0, None, T=T, nrec=4, dt_frac=DTF)
        if sgn < 0:
            d = np.abs(ln - ln_cosh(w * t))
            okP3a &= bool(np.all(d <= TOL(w * t)))
            print(f"   w = {w:5.1f} kappa = -w  worst |ln G - ln cosh| = {d.max():.4f}")
        else:
            okP3b &= bool(ln.max() <= np.log(1.001))
            print(f"   w = {w:5.1f} kappa = +w  max ln G = {ln.max():+.2e}")
ck("P3a kappa = -w grows as cosh(w t)", okP3a)
ck("P3b kappa = +w stays bounded (G <= 1 + 1e-3)", okP3b)

print("\n== P4  sequence at fixed t0 = 0.1")
t0 = 0.1
eps, sol, Gw = [], [], []
for w in WS:
    t, ln = run_dirichlet(w, t0, Vill, nrec=1)
    e = -ln_cosh(w * t0)
    eps.append(e)
    sol.append(e + ln[-1])
    t, lnw = run_dirichlet(w, t0, Vwell, nrec=1)
    Gw.append(lnw[-1])
    print(f"   w = {w:6.1f}  ln(data norm) = {e:10.3f}  ln(solution norm) = {e + ln[-1]:+8.4f}   well-posed family: ln G = {lnw[-1]:+.2e}")
ck("P4a ill-posed family: data norm -> 0 (strictly decreasing, < 1e-6), solution norm stays in [0.9, 1.1]",
   all(np.diff(eps) < 0) and eps[-1] < np.log(1e-6) and all(np.log(0.9) <= s <= np.log(1.1) for s in sol),
   f"solution norms {[round(float(np.exp(s)), 3) for s in sol]}")
ck("P4b well-posed family: G_n(t0) <= 1.001 for all n", max(Gw) <= np.log(1.001))
supill = max(s - e for s, e in zip(sol, eps))
ck("P4c the two families are told apart by sup G (>1e3 vs <= 1.001)", supill > np.log(1e3) and max(Gw) <= np.log(1.001),
   f"sup ln G ill = {supill:.2f}, well = {max(Gw):+.1e}")
ck.done()

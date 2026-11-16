"""
The sequence of data built from Theorem 1's exact modes (dg_proof), advanced in time.

For l = 2..LFD the data are u_l(0) = eps_l psi_l, v_l(0) = 0, psi_l = Psihat_l/||Psihat_l||_{L^2(dx)} (exact mode, modes.py),
eps_l = 1/cosh(omega_l t0), t0 = 0.01.  By Theorem 1 the solution of u_tt = -H u is u(t) = eps_l cosh(omega_l t) psi_l, so
||u_l(t0)||_{L^2(dx)} = 1 for every l while ||u_l(0)||_{L^2(dx)} = eps_l -> 0.  Here the PDE is advanced numerically (evolve.py,
the planted-checked scheme) and the growth G_l(t) = ||u(t)||/||u(0)|| is compared with cosh(omega_l t).

Thresholds, frozen from the planted runs (planted.py) before this script was run:
  growth      |ln G - ln cosh(w t)| <= 0.01 + 5e-4 w t        (planted worst was 0.2 of this)
  Rayleigh    |<u, A u>/<u, u> / (-w^2) - 1| <= 1.5e-3        (planted at the same resolution: 4e-4)
  resolution  40 points per 1/omega, dt = 0.5 * stability limit
Checks:
  S1  Rayleigh quotient of the exact mode in the discretised operator = -omega_l^2, l = 2..LFD.
  S2  growth at the recorded times in (0, t0], l = 2..LFD.
  S3  growth over T = 20/omega_l (l <= 12), which tests omega_l itself with growth e^{20}.
  S4  the sequence: ln||data|| strictly decreasing and below ln(1e-6) at l = LFD; ln||u(t0)|| within TOL of 0.
  S5  the grid norm of the data against the mp norm (modes.norm_psihat), l = 2..8, relative 1e-2.
Mutations: --mutate=omega (predict with 1.05 w) | vhat (V^ * 1.01 in the evolved operator) | sign (V^ -> -V^) | data (r^4 -> r^5 in the
  data) | dt (step 2.2 x the stability limit).
"""
import argparse
import os
import sys
import time

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import modes as M
import fastmodes as FM
from chk import Checker
from evolve import Op, evolve, ln_cosh

ap = argparse.ArgumentParser()
ap.add_argument('--mutate', default='')
ap.add_argument('--lfd', type=int, default=32)
args = ap.parse_args()
MU = args.mutate
LFD = args.lfd
T0 = 0.01
NPW = 40
TOL = lambda y: 0.01 + 5e-4 * np.asarray(y)
RQTOL = 1.5e-3
DTF = 2.2 if MU == 'dt' else 0.5
ck = Checker()


def build(l, T):
    w = float(M.omega(l))
    h = 1.0 / (NPW * w)
    n = int(round((30.0 / w + T) / h))
    x = h * np.arange(n + 1)
    r = FM.r_of_x(x[1:])
    V = np.zeros_like(x)
    V[1:] = FM.vhat(l, r)
    if MU == 'sign':
        V = -V
    if MU == 'vhat':
        V = V * 1.01
    op = Op(x, V, 'dirichlet')
    rr = FM.r_of_x(op.xs)
    u = FM.psihat(l, rr) * (rr if MU == 'data' else 1.0)
    u = u / np.sqrt(op.norm2(u))
    return op, u, w


wpred = lambda w: 1.05 * w if MU == 'omega' else w
rows = []
okS1 = okS2 = True
worstRQ = worstG = 0.0
t_start = time.time()
print("   l      omega   w*t0   log10 eps_l   ln G(t0)  ln cosh   diff     RQ/(-w^2)-1   N      steps   sec")
for l in range(2, LFD + 1):
    ts = time.time()
    op, u, w = build(l, T0)
    rq = float(np.sum(op.m * u * op.apply(u)) / np.sum(op.m * u * u)) / (-w * w) - 1
    t, ln, _, _ = evolve(op, u, None, T=T0, nrec=4, dt_frac=DTF)
    ex = ln_cosh(wpred(w) * t)
    d = ln - ex
    okS1 &= abs(rq) <= RQTOL
    okS2 &= bool(np.all(np.abs(d) <= TOL(w * t)))
    worstRQ = max(worstRQ, abs(rq) / RQTOL)
    worstG = max(worstG, float(np.max(np.abs(d) / TOL(w * t))))
    rows.append((l, w, -ln_cosh(w * T0), ln[-1]))
    print(f"{l:4d} {w:10.1f} {w * T0:7.2f} {-ln_cosh(w * T0) / np.log(10):12.3f} {ln[-1]:11.3f} {ex[-1]:9.3f} {d[-1]:+8.4f} {rq:+12.2e} "
          f"{len(op.xs):7d} {int(T0 / (DTF * op.dt)):7d} {time.time() - ts:5.1f}", flush=True)
ck("S1 Rayleigh quotient of the exact mode = -omega_l^2 in the discretised operator", okS1, f"worst |RQ-err|/tol = {worstRQ:.2f}")
ck("S2 growth G_l(t) = cosh(omega_l t), 4 times in (0, t0], l = 2..%d" % LFD, okS2, f"worst |diff|/TOL = {worstG:.2f}")

print("\n== S3 growth over T = 20/omega_l")
okS3 = True
w3 = 0.0
print("   l      omega   w*T    ln G(T)   ln cosh(wT)   diff")
for l in range(2, 13):
    w = float(M.omega(l))
    T = 20.0 / w
    op, u, w = build(l, T)
    t, ln, _, _ = evolve(op, u, None, T=T, nrec=4, dt_frac=DTF)
    ex = ln_cosh(wpred(w) * t)
    d = ln - ex
    okS3 &= bool(np.all(np.abs(d) <= TOL(w * t)))
    w3 = max(w3, float(np.max(np.abs(d) / TOL(w * t))))
    print(f"{l:4d} {w:10.1f} {w * T:6.1f} {ln[-1]:10.3f} {ex[-1]:12.3f} {d[-1]:+9.4f}")
ck("S3 growth over T = 20/omega_l, l = 2..12 (ln G = 20 - ln 2)", okS3, f"worst |diff|/TOL = {w3:.2f}")

leps = np.array([r[2] for r in rows])
lu = np.array([r[2] + r[3] for r in rows])
ck("S4a ln||data_l|| = -ln cosh(omega_l t0) strictly decreasing in l, below ln(1e-6) at l = %d" % LFD,
   bool(np.all(np.diff(leps) < 0)) and leps[-1] < -np.log(1e6), f"ln||data|| from {leps[0]:.3f} to {leps[-1]:.1f}")
ck("S4b ln||u_l(t0)||_{L^2(dx)} = 0 within TOL for every l (the solution does not tend to 0)",
   all(abs(v) <= TOL(r[1] * T0) for v, r in zip(lu, rows)), f"max |ln||u(t0)||| = {np.max(np.abs(lu)):.3f}")

print("\n== S5 grid norm of the data against the mp norm")
worst5 = 0.0
for l in range(2, 9):
    w = float(M.omega(l))
    h = 1.0 / (NPW * w)
    x = h * np.arange(1, int(round(30.0 / w / h)) + 1)
    r = FM.r_of_x(x)
    u = FM.psihat(l, r)
    gridn = np.sqrt(np.sum(u * u) * h)
    mpn = float(M.norm_psihat(l))
    worst5 = max(worst5, abs(gridn / mpn - 1))
    print(f"   l = {l}  grid {gridn:.8e}  mp {mpn:.8e}  rel {gridn / mpn - 1:+.1e}")
ck("S5 grid L^2(dx) norm = mp norm of Psihat_l, l = 2..8, relative 1e-2", worst5 < 1e-2, f"worst {worst5:.1e}")
print(f"\ntotal {time.time() - t_start:.0f} s")
ck.done()

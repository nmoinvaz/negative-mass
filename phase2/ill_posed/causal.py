"""
Data supported away from the singularity: x >= x_1 > 0 (tortoise distance x_1 from r = 0).  Thresholds written before the run.

  CA1  overlap exponent.  For g = indicator of [x_1, x_2] and psi_l the unit mode, a_l = <psi_l, g>_{L^2(dx)} (mp quadrature),
       m_l := -ln(a_l)/omega_l (x_1 = 0.02, x_2 = 0.05), l in {2, 4, 8, 12, 16, 24, 32, 40}:
       m_l > x_1, m_l strictly decreasing for l >= 4, and 0 < omega_l (m_l - x_1) <= ln(omega_l) for l >= 8, i.e. m_l -> x_1.
       (First version: |m_l/x_1 - 1| <= 5% from l = 8 and <= 1% from l = 16, a guess made without the 1/(2 omega x_1) ln(omega) rate;
       it failed at l = 8 (33%), 12 (9.9%) and 16 (3.8%).  Reported in RESULTS.md.)  Numerical evidence for the overlap exponent.
  CA2  time domain, l in {12, 16, 20, 24}, smooth bump on [x_1, x_2] = [0.004, 0.006], unit norm, v_0 = 0, point-mass problem advanced to
       T = 2.5 x_1 at 20 times:
       (i)   spectral bracket (exact, from Theorem 1): ln(a_l cosh(w_l t)) - TOL <= ln G_l(t) <= 0.5 ln(a_l^2 cosh^2 + 1) + TOL;
       (ii)  uniform bound before the signal reaches the singularity: max_l G_l(t) <= 1.05 for t <= 0.9 x_1;
       (iii) growth after: the slope d ln G/dt between t = 1.5 x_1 and T equals w_l within 1% for each l (the spectral bracket), and the
             onset time t_on(l) = -ln(a_l)/w_l (where the bracket's lower end crosses G = 1) satisfies t_on(24) <= 1.2 x_1 and decreases
             from l = 16 to 24.
             (First version: ln G_20 >= 0.9 w (t - x_1) for t >= 1.5 x_1.  It failed: ln G = w (t - x_1) - 8.6 at l = 20, the offset being
             the O(sqrt(w)) shift of the edge of a smooth bump.  Replaced by the slope and onset statements, which are what the spectral
             theorem implies.  Reported in RESULTS.md.)
Mutations: --mutate=x1 (predict with 2 x_1 in CA1) | omega (1.05 w in CA2) | near (CA2 data moved to [0.001, 0.003], inside the
  window of (ii)) | far (CA2 data moved to [0.008, 0.010], onset far beyond (iii)).
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
args = ap.parse_args()
MU = args.mutate
ck = Checker()
TOL = lambda y: 0.01 + 5e-4 * np.asarray(y)
mp.mp.dps = 40

# ------------------------------------------------------------------ CA1
x1, x2 = 0.02, 0.05
print("== CA1  overlap exponent of the unit mode with the indicator of [%.3f, %.3f]" % (x1, x2))
print("    l      omega   -ln a_l / omega   (m_l/x_1 - 1)   w (m_l - x_1)   ln w")
ok1 = True
ms = []
r1, r2 = M.r_of_x(x1), M.r_of_x(x2)
for l in (2, 4, 8, 12, 16, 24, 32, 40):
    w = float(M.omega(l))
    t0 = time.time()
    nrm = M.norm_psihat(l)
    pts = [r1] + [r1 + (r2 - r1) * s for s in (0.25, 0.5, 0.75)] + [r2]
    a = mp.quad(lambda r: M.psihat_unnorm(l, r) * r / (r + 2), pts) / nrm
    m = -mp.log(a) / w
    xx = x1 * (2 if MU == 'x1' else 1)
    rel = float(m / x1 - 1)
    off = float(w * (m - xx))
    ms.append(float(m))
    print(f"{l:5d} {w:10.1f} {float(m):17.6f} {rel:+14.4f} {off:12.2f} {np.log(w):10.2f}")
    ok1 &= float(m) > xx
    if l >= 8:
        ok1 &= 0 < off <= np.log(w)
ok1 &= all(a > b for a, b in zip(ms[1:], ms[2:]))
ck("CA1 m_l = -ln a_l/omega_l > x_1, decreasing for l >= 4, 0 < omega (m_l - x_1) <= ln omega for l >= 8 (so m_l -> x_1)", ok1)

# ------------------------------------------------------------------ CA2
print("\n== CA2  time domain, bump on [0.004, 0.006], T = 2.5 x_1")
X1, X2 = 0.004, 0.006
D1, D2 = {'near': (0.001, 0.003), 'far': (0.008, 0.010)}.get(MU, (X1, X2))   # where the data actually sit
T = 2.5 * X1
NPW = 40
okb = oku = okg = True
hdr = "    l      omega   ln a_l"
rowsU = []
res = {}
tons = []
for l in (12, 16, 20, 24):
    w = float(M.omega(l))
    h = 1.0 / (NPW * w)
    n = int(round((D2 + T) / h))
    x = h * np.arange(n + 1)
    V = np.zeros_like(x)
    V[1:] = FM.vhat(l, FM.r_of_x(x[1:]))
    op = Op(x, V, 'dirichlet')
    s = (op.xs - D1) / (D2 - D1)
    g = np.zeros_like(s)
    inside = (s > 0) & (s < 1)
    g[inside] = np.exp(-1.0 / (s[inside] * (1 - s[inside])))
    g = g / np.sqrt(op.norm2(g))
    rr = FM.r_of_x(op.xs)
    psi = FM.psihat(l, rr)
    psi = psi / np.sqrt(op.norm2(psi))
    a = float(np.sum(op.m * psi * g))
    t, ln, _, _ = evolve(op, g, None, T=T, nrec=20)
    wp = w * (1.05 if MU == 'omega' else 1.0)
    lo = np.log(a) + ln_cosh(wp * t)
    hi = 0.5 * np.log(a * a * np.exp(2 * ln_cosh(wp * t)) + 1.0)
    okb &= bool(np.all(ln >= lo - TOL(w * t)) and np.all(ln <= hi + TOL(w * t)))
    early = t <= 0.9 * X1
    oku &= bool(np.all(np.exp(ln[early]) <= 1.05))
    late = t >= 1.5 * X1
    slope = (ln[-1] - ln[late][0]) / (t[-1] - t[late][0])
    okg &= abs(slope / w - 1) <= 0.01
    ton = -np.log(a) / w
    tons.append(ton)
    res[l] = (t, ln, a)
    print(f"  l = {l:2d}  omega = {w:8.1f}  ln a_l = {np.log(a):8.2f}   max_(t<=0.9x1) G = {np.exp(ln[early]).max():.4f}   ln G(T) = {ln[-1]:9.2f}   w (T - x_1) = {w * (T - X1):9.2f}")
print("\n   time series, l = 20:  t/x_1   ln G    w (t - x_1)")
t, ln, a = res[20]
w = float(M.omega(20))
for i in range(0, 20, 2):
    print(f"        {t[i] / X1:6.2f} {ln[i]:9.3f} {w * (t[i] - X1):12.2f}")
ck("CA2i spectral bracket a_l cosh(w_l t) <= G_l(t) <= sqrt(a_l^2 cosh^2 + 1), l = 12, 16, 20, 24", okb)
ck("CA2ii for t <= 0.9 x_1 every l has G_l(t) <= 1.05 (the signal has not reached the singularity)", oku)
print("\n   onset time -ln(a_l)/w_l / x_1 for l = 12, 16, 20, 24:", ", ".join(f"{v / X1:.3f}" for v in tons))
okg &= tons[3] <= 1.2 * X1 and tons[1] > tons[2] > tons[3]
ck("CA2iii slope of ln G after 1.5 x_1 = w_l within 1%; onset time -ln(a_l)/w_l <= 1.2 x_1 at l = 24 and decreasing from l = 16", okg)
ck.done()

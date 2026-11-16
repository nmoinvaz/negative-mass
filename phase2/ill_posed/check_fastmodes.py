"""
fastmodes.py (float64, vectorised) against modes.py (mp, reference).  Each check can fail.
  A  r(x) round trip and x(r) against mp, r from 1e-6 to 1e3.
  B  Vhat(float) against Vhat(mp) at 300 log-spaced points between 1e-3 r_s and 200 r_s, l = 2, 5, 12, 30, 60, relative 1e-7 of
     max(|Vhat|, 1).
  C  Psihat(float)/Psihat(mp) at the same points, relative 1e-9.
  D  kappa_c(float) against mp, away from r_s.
Mutations: --mutate=coef (c_3 -> 1.001 c_3 in the float coefficients) | xr (x(r) off by 1e-6 relative).
"""
import argparse
import os
import sys

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import modes as M
import fastmodes as FM
from chk import Checker

ap = argparse.ArgumentParser()
ap.add_argument('--mutate', default='')
args = ap.parse_args()
ck = Checker()
if args.mutate == 'coef':
    old = FM._coef.__wrapped__

    def patched(l):
        ch, pi = old(l)
        ch = ch.copy()
        ch[3] *= 1.001
        return ch, pi
    FM._coef = patched
if args.mutate == 'xr':
    oldx = FM.x_of_r
    FM.x_of_r = lambda r: oldx(r) * (1 + 1e-6)

r = np.logspace(-6, 3, 400)
xr = np.array([float(M.x_of(mp.mpf(float(v)))) for v in r])
ck("A x(r) float against mp, relative 1e-12", np.max(np.abs(FM.x_of_r(r) / xr - 1)) < 1e-12, f"worst {np.max(np.abs(FM.x_of_r(r) / xr - 1)):.1e}")
back = FM.r_of_x(xr)
ck("A r(x(r)) = r, relative 1e-11", np.max(np.abs(back / r - 1)) < 1e-11, f"worst {np.max(np.abs(back / r - 1)):.1e}")

wV = wP = wK = 0.0
for l in (2, 5, 12, 30, 60):
    rs = 6.0 / M.cc(l)
    rr = rs * np.logspace(-3, np.log10(200), 300)
    Vf = FM.vhat(l, rr)
    Pf = FM.psihat(l, rr)
    Kf = FM.kappa_c(l, rr)
    for i, v in enumerate(rr):
        Vm = float(M.Vhat(l, mp.mpf(float(v))))
        Pm = float(M.psihat_unnorm(l, mp.mpf(float(v))))
        Km = float(M.kappa_c(l, mp.mpf(float(v))))
        wV = max(wV, abs(Vf[i] - Vm) / max(abs(Vm), 1.0))
        if Pm != 0.0 and abs(Pm) > 1e-290:
            wP = max(wP, abs(Pf[i] / Pm - 1))
        if abs(v - rs) > 0.05 * rs:
            wK = max(wK, abs(Kf[i] - Km) / max(abs(Km), 1.0))
ck("B Vhat float vs mp, l = 2, 5, 12, 30, 60", wV < 1e-7, f"worst {wV:.1e}")
ck("C Psihat float vs mp", wP < 1e-9, f"worst {wP:.1e}")
ck("D kappa_c float vs mp (|r - r_s| > 5%)", wK < 1e-8, f"worst {wK:.1e}")
ck.done()

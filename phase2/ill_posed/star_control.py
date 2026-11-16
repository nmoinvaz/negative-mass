"""
Control: the same construction on the truncated star (dg_proof Theorem 2: x >= x_R = x(R), Robin u_x = kappa u, Dirichlet = kappa = oo).

Thresholds written before the run:
  SC0  R = 0.2/c (0.02/c was tried first: the tridiagonal eigensolver loses the eigenvalue, 1e-16 x 2/h^2 ~ omega^2 at l >= 24; the
       count was 1 throughout, only E_0 was off): Dirichlet and Neumann each have exactly one negative eigenvalue, within 1e-2 of -omega_l^2 (the limit R -> 0 of Theorem 1,
       by a code that shares nothing with dg_spectrum), l = 2..60.
  SC1  the number of negative eigenvalues of the discretised Robin problem (Sylvester inertia of the tridiagonal stiffness matrix; grid
       geometric in r, 25 points per e-fold of r / (l+1)) equals Theorem 2's  [R < r_s] + [kappa < kappa_c(R)]  in every case with
       |R - r_s|/r_s > 1e-3 and |kappa - kappa_c| > 2e-3 (|kappa_c| + 1);  R in {0.05, 0.12, 0.3, 1, 3}, kappa in {Dirichlet, 0, 1, 10},
       l = 2..60.  For each (R, kappa) the unstable l form an initial segment {2..l_U} within the range, and
       Omega(R, kappa) = max_l sqrt(-E_0) is recorded.
  SC2  Vhat >= 0 on [R, 50] (the sufficient condition of the proof) for every l with c r_neg-bound 10.7/c < R, up to l = 60.
  SC3  r_* and r_neg from float root finding against the certified brackets of dg_proof/threshold_isolation_output.txt, relative 1e-8
       (l = 2, 3, 5, 10, 20, 40, 80, 100).
  SC4  same data, two problems (R = 0.12, t0 = 0.01, l = 2..24): data = the exact mode of the point mass restricted to x >= x_R, unit
       norm.  (a) point mass, evolved on [0, X]; (b) star, evolved on [x_R, X] with Robin kappa.  Spectral-theorem brackets:
         star:       ln G_l(t) <= ln cosh(Omega t) + TOL            (v0 = 0, all l; Omega from SC1)
         point mass: ln(a_l cosh(w_l t)) - TOL <= ln G_l(t) <= 0.5 ln(a_l^2 cosh^2 + 1) + TOL,   a_l = <psi_l, data>
       and the uniformity statement: max_l G^star_l(t0) <= cosh(Omega t0) while G^pm_l(t0) exceeds it for large l.
  SC5  caveat: kappa_l = -l^2 gives growth sqrt(-E_0) tracking |kappa_l| (ratio >= 0.8 for l = 20..60, >= 0.98 at l = 60, increasing in l).
       (First version of this threshold was 0.9 at every l, a guess; at l = 20 the ratio is 0.83.  The statement was changed to the one made in
       RESULTS.md; reported there.)
Mutations: --mutate=omega0 (SC0 compared with 1.1 omega) | kappa (kappa - 1e4 in SC1) | rs (r_s -> 1.1 r_s in the prediction) | sign (Robin sign flipped in SC4) |
  bound (SC4 bound with cosh(0.5 Omega t)) | kl (SC5 kappa = -1) | nodes (SC2 with Vhat - 1 at large r).
"""
import argparse
import os
import re
import sys
import time

import mpmath as mp
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import modes as M
import fastmodes as FM
from chk import Checker
from evolve import Op, evolve, ln_cosh

ap = argparse.ArgumentParser()
ap.add_argument('--mutate', default='')
ap.add_argument('--lmax', type=int, default=60)
ap.add_argument('--lpm', type=int, default=24)
args = ap.parse_args()
MU = args.mutate
ck = Checker()
LMAX = args.lmax
RS = [0.05, 0.12, 0.3, 1.0, 3.0]
KS = ['D', 0.0, 1.0, 10.0]
TOL = lambda y: 0.01 + 5e-4 * np.asarray(y)
t_all = time.time()


def star_op(l, R, kappa, rmax=50.0, dsf=0.04):
    ds = dsf / (l + 1)
    n = int(np.ceil(np.log(rmax / R) / ds))
    r = R * np.exp(ds * np.arange(n + 1))
    x = FM.x_of_r(r)
    V = FM.vhat(l, r)
    if kappa == 'D':
        return Op(x, V, 'dirichlet')
    return Op(x, V, 'robin', kappa)


# ---------------------------------------------------------------- SC0
print("== SC0  R -> 0 limit: Dirichlet / Neumann at c r = 0.2 has one negative eigenvalue, -omega_l^2 to 1e-2, l = 2..%d" % LMAX)
ok0 = True
worst0 = 0.0
for l in range(2, LMAX + 1):
    Rl = 0.2 / M.cc(l)
    w2 = float(M.omega(l))**2
    for K in ('D', 0.0):
        op = star_op(l, Rl, K)
        cnt = op.count_negative()
        E0 = op.lowest(1)[0]
        e = abs(E0 / (-w2) - 1)
        worst0 = max(worst0, e)
        ok0 &= (cnt == 1 and abs(E0 / (-w2 * (1.21 if MU == 'omega0' else 1.0)) - 1) < 1e-2)
ck("SC0 one negative eigenvalue, E_0 = -omega_l^2 (1e-2), Dirichlet and Neumann at c r = 0.2, l = 2..%d" % LMAX, ok0, f"worst rel {worst0:.1e}")

# ---------------------------------------------------------------- SC1
print("== SC1  count of negative eigenvalues vs Theorem 2, l = 2..%d" % LMAX)
bad = []
skipped = 0
ncase = 0
Omega = {}
Ounst = {}
for R in RS:
    for K in KS:
        unst = []
        om = 0.0
        for l in range(2, LMAX + 1):
            c = M.cc(l)
            rs = 6.0 / c * (1.1 if MU == 'rs' else 1.0)
            if abs(R - 6.0 / c) / (6.0 / c) < 1e-3:
                skipped += 1
                continue
            kc = float(FM.kappa_c(l, np.array([R]))[0])
            kap = K if K == 'D' else K - (1e4 if MU == 'kappa' else 0.0)
            if K != 'D' and abs(K - kc) < 2e-3 * (abs(kc) + 1):
                skipped += 1
                continue
            pred = int(R < rs) + (0 if K == 'D' else int(K < kc))
            op = star_op(l, R, kap)
            cnt = op.count_negative()
            ncase += 1
            if cnt != pred:
                bad.append((R, K, l, cnt, pred))
            if cnt > 0:
                unst.append(l)
                E0 = op.lowest(1)[0]
                if E0 < 0 and np.sqrt(-E0) > om:
                    om = float(np.sqrt(-E0))
        Omega[(R, K)] = om
        Ounst[(R, K)] = unst
        segs = f"{unst[0]}..{unst[-1]}" if unst else "none"
        contiguous = (not unst) or unst == list(range(unst[0], unst[-1] + 1))
        print(f"   R = {R:5.2f}  kappa = {str(K):>4s}  unstable l: {segs:8s} ({len(unst)}), Omega = {om:10.3f}, contiguous: {contiguous}", flush=True)
ck("SC1 negative-eigenvalue count = [R < r_s] + [kappa < kappa_c(R)] in every case", not bad,
   f"{ncase} cases, {skipped} skipped near thresholds; first mismatches {bad[:4]}")
ck("SC1b for kappa >= 0 or Dirichlet the unstable l are finitely many and form {l <= l_U} inside the tested range (l_U < %d)" % LMAX,
   all((not Ounst[k]) or (Ounst[k][-1] < LMAX and Ounst[k] == list(range(Ounst[k][0], Ounst[k][-1] + 1))) for k in Ounst))
print("\n   Omega(R, kappa = 0) * R^2 (the cap scales like 1/R^2):",
      ", ".join(f"R={R}: {Omega[(R, 0.0)] * R * R:.2f}" for R in RS))

# ---------------------------------------------------------------- SC2
print("\n== SC2  Vhat >= 0 on [R, 50] where 10.7/c < R (the sufficient condition)")
ok2 = True
n2 = 0
worst = np.inf
for R in RS:
    for l in range(2, LMAX + 1):
        if 10.7 / M.cc(l) < R:
            r = R * np.exp(np.linspace(0, np.log(50.0 / R), 4000))
            V = FM.vhat(l, r)
            if MU == 'nodes':
                V = V - 1.0 * (r > 10)
            n2 += 1
            mn = float(V.min())
            worst = min(worst, mn)
            if mn < 0:
                ok2 = False
ck("SC2 min Vhat on [R, 50] >= 0 whenever 10.7/c < R (R in %s, l <= %d, %d cases)" % (RS, LMAX, n2), ok2, f"smallest min = {worst:.3e}")

# ---------------------------------------------------------------- SC3
print("\n== SC3  r_*, r_neg (float) against the certified brackets")
tbl = {}
for line in open(os.path.join(HERE, '..', 'dg_proof', 'threshold_isolation_output.txt')):
    m = re.match(r"\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", line)
    if m:
        tbl[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
worst3 = 0.0
n3 = 0
for l in (2, 3, 5, 10, 20, 40, 80, 100):
    c = M.cc(l)
    ru = np.exp(np.linspace(np.log(0.5 / c), np.log(40.0 / c), 6000))
    V = FM.vhat(l, ru)
    s = np.sign(V)
    idx = np.nonzero((s[:-1] < 0) & (s[1:] > 0))[0][-1]
    rneg = brentq(lambda r: FM.vhat(l, np.array([r]))[0], ru[idx], ru[idx + 1], xtol=1e-16, rtol=1e-14)
    kc = FM.kappa_c(l, ru)
    rsl = 6.0 / c
    j = np.nonzero((ru[:-1] > rsl * 1.001) & (kc[:-1] > 0) & (kc[1:] < 0))[0][0]
    rstar = brentq(lambda r: FM.kappa_c(l, np.array([r]))[0], ru[j], ru[j + 1], xtol=1e-16, rtol=1e-14)
    ref = tbl[l]
    e = max(abs(c * rstar / ref[0] - 1), abs(c * rneg / ref[1] - 1))
    worst3 = max(worst3, e)
    n3 += 1
    print(f"   l = {l:4d}  c r_* = {c * rstar:.10f} (table {ref[0]:.10f})   c r_neg = {c * rneg:.10f} (table {ref[1]:.10f})")
ck("SC3 float r_*, r_neg = certified table, l = 2..100", worst3 < 1e-8 and n3 == 8, f"worst rel {worst3:.1e}")

# ---------------------------------------------------------------- SC4
print("\n== SC4  same data, point mass vs star, R = 0.12, t0 = 0.01")
R4 = 0.12
T0 = 0.01
NPW = 40
xR = float(M.x_of(mp.mpf(R4)))
print(f"   x_R = {xR:.6f}  (t0 - x_R = {T0 - xR:.6f})")
for K in ('D', 0.0, 10.0):
    Om = Omega[(R4 if R4 in RS else 0.12, K)]
    okstar = okpm = True
    gmax = -1.0
    gpm_big = []
    rows = []
    for l in range(2, args.lpm + 1):
        w = float(M.omega(l))
        h = 1.0 / (NPW * w)
        Lx = 30.0 / w + 0.5 * T0
        n = int(round((xR + Lx) / h))
        # point mass: grid from 0
        x = h * np.arange(n + 1)
        V = np.zeros_like(x)
        V[1:] = FM.vhat(l, FM.r_of_x(x[1:]))
        opp = Op(x, V, 'dirichlet')
        rr = FM.r_of_x(opp.xs)
        u = FM.psihat(l, rr)
        psi = u / np.sqrt(opp.norm2(u))
        u = np.where(opp.xs >= xR, u, 0.0)
        u = u / np.sqrt(opp.norm2(u))
        a = float(np.sum(opp.m * psi * u))
        tp, lnp, _, _ = evolve(opp, u, None, T=T0, nrec=4)
        # star: grid x_R + j h
        xs_ = xR + h * np.arange(int(round(Lx / h)) + 1)
        Vs = FM.vhat(l, FM.r_of_x(xs_))
        kap = K
        if K != 'D' and MU == 'sign':
            kap = -K - 1.0
        ops = Op(xs_, Vs, 'dirichlet') if K == 'D' else Op(xs_, Vs, 'robin', kap)
        us = FM.psihat(l, FM.r_of_x(ops.xs))
        us = us / np.sqrt(ops.norm2(us))
        ts, lns, _, _ = evolve(ops, us, None, T=T0, nrec=4)
        bnd = np.log(np.cosh(min(Om * (0.5 if MU == 'bound' else 1.0) * ts[-1], 700))) if Om > 0 else 0.0
        okstar &= bool(np.all(lns <= np.log(np.cosh(np.minimum((0.5 if MU == 'bound' else 1.0) * Om * ts, 700))) + TOL(Om * ts) + 1e-3))
        lo = np.log(a) + ln_cosh(w * tp)
        hi = 0.5 * np.log(a * a * np.exp(2 * ln_cosh(w * tp)) + 1.0) if w * tp[-1] < 300 else lo + 1e-9
        okpm &= bool(np.all(lnp >= lo - TOL(w * tp)) and np.all(lnp <= hi + TOL(w * tp)))
        gmax = max(gmax, lns[-1])
        rows.append((l, a, lnp[-1], lo[-1], hi[-1], lns[-1]))
    print(f"   kappa = {str(K):>4s}   Omega = {Om:.2f}   star: max_l ln G(t0) = {gmax:.3f}   bound ln cosh(Omega t0) = {np.log(np.cosh(Om * T0)):.3f}")
    print("      l    ln a_l   ln G_pm(t0)   [lower, upper]            ln G_star(t0)")
    for (l, a, g, lo, hi, gs) in rows:
        if l in (2, 4, 6, 8, 10, 12, 16, 20, 24) or l == args.lpm:
            print(f"     {l:3d} {np.log(a):9.2f} {g:12.3f}   [{lo:9.3f}, {hi:9.3f}]   {gs:10.3f}")
    ck(f"SC4a star (kappa = {K}): ln G_l(t) <= ln cosh(Omega t) for all l = 2..{args.lpm}, 4 times", okstar)
    ck(f"SC4b point mass (kappa = {K} data): spectral bracket a_l cosh(w_l t) <= G_l(t) <= sqrt(a_l^2 cosh^2 + 1), all l", okpm)
    ck(f"SC4c uniform in l for the star, unbounded for the point mass (kappa = {K}): max_l G^star <= cosh(Omega t0), G^pm_{args.lpm} > 1e6 cosh(Omega t0)",
       gmax <= np.log(np.cosh(Om * T0)) + 0.05 and rows[-1][2] > np.log(1e6) + np.log(np.cosh(Om * T0)),
       f"max ln G^star = {gmax:.3f}, ln G^pm_{args.lpm} = {rows[-1][2]:.1f}")

# ---------------------------------------------------------------- SC5
print("\n== SC5  caveat: kappa_l = -l^2 (l-dependent, negative), R = 0.3")
gr = []
for l in (20, 30, 40, 60):
    kl = -float(l * l) if MU != 'kl' else -1.0
    op = star_op(l, 0.3, kl)
    E0 = op.lowest(1)[0]
    gr.append((l, float(np.sqrt(max(-E0, 0))), abs(kl)))
    print(f"   l = {l}  kappa = {kl:.0f}  sqrt(-E_0) = {gr[-1][1]:.2f}   ratio to |kappa| = {gr[-1][1] / abs(kl):.3f}")
ck("SC5 growth tracks |kappa_l|: ratio >= 0.8 at l = 20..60, ratio at l = 60 >= 0.98, growth strictly increasing",
   all(g[1] >= 0.8 * g[2] for g in gr) and gr[-1][1] >= 0.98 * gr[-1][2] and all(a[1] < b[1] for a, b in zip(gr, gr[1:])))
print(f"\ntotal {time.time() - t_all:.0f} s")
ck.done()

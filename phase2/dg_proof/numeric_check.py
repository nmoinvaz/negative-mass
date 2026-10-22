"""Numerical cross-checks at l = 9, 10 (beyond the dg_spectrum runs): full-problem negative count and
truncated count formula, reusing the dg_spectrum discretisations."""
import sys, numpy as np
sys.path.insert(0, '/home/claude/negative-mass/phase2/dg_spectrum')
import common as cm
from full_spectrum import full_eigs
from truncated import trunc_eigs, n_neg
for l in (9, 10):
    k = cm.kval(l)
    for N, a in ((200, 0.1), (240, 0.3)):
        w = full_eigs(l, N, a)[0]
        w = np.real(w[np.isfinite(w)])
        neg = np.sort(w[w < -1e-6])
        print(f"full l={l} N={N} a={a}: #neg={len(neg)}, E0={neg[0]:.10g}, exact={-(k/2)**2}", flush=True)
bad = 0; tot = 0
for l in (6, 8, 10):
    rs = cm.r_s(l)
    for R in (0.3*rs, 0.8*rs, 1.2*rs, 2*rs, 5*rs, 20*rs):
        x_R = float(cm.x_of_r(R)); kc = cm.numeric(l)['kappa_c'](R)
        for kx in (-10, -1, -0.3, 0, 1, 10):
            kap = kx/x_R
            pred = int(R < rs) + int(kap < kc)
            got = n_neg(trunc_eigs(l, R, kap, N=200))
            tot += 1; bad += (pred != got)
            if pred != got: print("MISMATCH", l, R, kap, kc, pred, got)
print(f"truncated count formula: {tot - bad}/{tot} agree (l = 6, 8, 10)")

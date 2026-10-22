"""Signed momentum residual vs time, and whether it is just (mass mismatch + velocity mismatch). Data only."""
import sys, numpy as np
C = sys.argv[1]
for n in (128,192,256):
    D = np.loadtxt(f'{C}/runaway_pair_d10_L64_N{n}_lev0/sector_dynamics.dat', comments='#'); t = D[:,0]
    pc, pp, pt = D[:,11], D[:,12], D[:,13]
    vc = np.gradient(D[:,1], t); vp = np.gradient(D[:,5], t)
    k = t >= 20
    rel = pt[k]/pc[k]
    print(f'N={n}: |P_tot/P+| over t>=20: max {np.abs(rel).max():.3%} median {np.median(np.abs(rel)):.3%}; at end {pt[-1]/pc[-1]:+.3%}')
    # ratio of raw momenta vs ratio of velocities (smoothed over 20 time units)
    def sm(x, w=25): return np.convolve(x, np.ones(w)/w, mode='same')
    for tt in (60,100,150,190):
        i = np.argmin(abs(t-tt))
        print(f'   t={t[i]:.0f}: p-/p+ = {pp[i]/pc[i]:.4f}   v-/v+ (smoothed) = {sm(vp)[i]/sm(vc)[i]:.4f}   p+/v+ = {pc[i]/sm(vc)[i]:.5f}  p-/v- = {pp[i]/sm(vp)[i]:.5f}')

"""Probe time series of headline cells for sampling irregularities / kinks. Data only."""
import sys, numpy as np
C = sys.argv[1]
for n in (128, 192, 256):
    d = np.loadtxt(f'{C}/runaway_pair_d10_L64_N{n}_lev0/sector_barycenters.dat', comments='#')
    t = d[:,0]; m = 0.5*(d[:,2]+d[:,7]); dt = np.diff(t)
    print(f'N={n}: rows {len(t)} t0 {t[0]} t_end {t[-1]:.3f} dt min {dt.min():.3f} max {dt.max():.3f} unique dts {np.unique(np.round(dt,3))[:8]}; nonmonotone t: {(dt<=0).sum()}')
    # velocity by finite differences, smoothed
    v = np.gradient(m, t)
    for tt in (5,25,50,75,100,125,150,175,199):
        i = np.argmin(abs(t-tt)); print(f'   t={t[i]:6.1f} x-x0={m[i]-m[0]:+.5f} v={v[i]:+.5e}')
    # residual of global quadratic fit
    k = t >= 5; c = np.polyfit(t[k], m[k], 2); r = m[k]-np.polyval(c, t[k])
    print(f'   global quad a={2*c[0]:.4e}, rms resid {np.sqrt(np.mean(r**2)):.2e}, max {np.abs(r).max():.2e}')
    # cubic term significance
    c3 = np.polyfit(t[k], m[k], 3); print(f'   cubic coeff {c3[0]:+.3e} -> da/dt {6*c3[0]:+.3e}; a(t=5)={2*c3[1]+6*c3[0]*5:.4e} a(t=200)={2*c3[1]+6*c3[0]*200:.4e}')
    # windowed fits, centred
    for a0, a1 in ((5,100),(100,200),(5,50),(50,100),(100,150),(150,200)):
        q = (t>=a0)&(t<=a1); cc = np.polyfit(t[q]-t[q].mean(), m[q], 2); print(f'     [{a0},{a1}] a={2*cc[0]:.4e} n={q.sum()}')

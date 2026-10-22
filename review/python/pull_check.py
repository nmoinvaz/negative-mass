"""Independent separation-corrected pull fits (x_k = c0 + c1 t + C * I2[1/d^2]) on the core tracker,
mass-ratio cells; plus last-two-thirds fits of the ladder. Data only."""
import sys, numpy as np
C = sys.argv[1]
def D(c): return np.loadtxt(f'{C}/{c}/sector_dynamics.dat', comments='#')
def B(c): return np.loadtxt(f'{C}/{c}/sector_barycenters.dat', comments='#')
def cumtrapz(y, t): return np.r_[0, np.cumsum(0.5*(y[1:]+y[:-1])*np.diff(t))]
def pull(c, sepcol=9):
    d = D(c); t = d[:,0]; g = np.abs(d[:,sepcol]); I2 = cumtrapz(cumtrapz(1/g**2, t), t)
    m = t >= 5; A = np.vstack([np.ones(m.sum()), t[m], I2[m]]).T
    return {k: np.linalg.lstsq(A, d[m, col], rcond=None)[0][2] for k, col in (('c',1),('p',5))}
base = pull('runaway_pair_d10_L64_N128_lev0')
print('matched constants', base)
cells = {'w0804':'massscale_pair_d10_w0804_L64_N128_lev0','w088':'massratio_w088_r060_d10_L64_N128_lev0','heavy':'massratio_heavyphantom_d10_L64_N128_lev0'}
pred_p = {'w0804':0.8026,'w088':0.5997,'heavy':0.7471}; pred_t = {'w0804':0.8059,'w088':0.6079,'heavy':0.7516}
for k, c in cells.items():
    f = pull(c); ch = 'p' if k=='heavy' else 'c'; ot = 'c' if ch=='p' else 'p'
    r = f[ch]/base[ch]
    print(f'{k}: changed pull ratio {r:.4f}  pred paper {pred_p[k]:.4f} (dev {r/pred_p[k]-1:+.2%})  pred true {pred_t[k]:.4f} (dev {r/pred_t[k]-1:+.2%}); control {f[ot]/base[ot]:.4f}')
# gap change at 200 (bary) vs mass ratio
for k, c in list(cells.items()) + [('matched','runaway_pair_d10_L64_N128_lev0')]:
    b = B(c); g = np.abs(b[:,2]-b[:,7]); print(f'  {k}: bary gap change {g[-1]-g[0]:+.4f}')
def fit(t, x, a0, a1):
    m=(t>=a0)&(t<=a1); return 2*np.polyfit(t[m]-t[m].mean(), x[m], 2)[0]
for n in (128,192,256):
    b = B(f'runaway_pair_d10_L64_N{n}_lev0'); t=b[:,0]; x=0.5*(b[:,2]+b[:,7])
    print(f'N={n}: bary last-two-thirds [{t[-1]/3:.1f},{t[-1]:.0f}] a={fit(t,x,t[-1]/3,t[-1]):.4e}')

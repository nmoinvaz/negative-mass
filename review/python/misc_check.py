"""Misc reproductions: 4-pt slope without d=20, per-cell window systematic in the d-scan,
exact disjoint windows on d=10, shift translational part. Data only."""
import sys, numpy as np
C = sys.argv[1]
def L(c, f): return np.loadtxt(f'{C}/{c}/{f}', comments='#')
def fit(t, x, a0, a1):
    m=(t>=a0)&(t<=a1); return 2*np.polyfit(t[m]-t[m].mean(), x[m], 2)[0]
ds = np.array([8,10,12,16,20]); a=[]; lo=[]; hi=[]
for d in ds:
    b = L(f'runaway_pair_d{d:02d}_L64_N128_lev0','sector_barycenters.dat'); t=b[:,0]; x=0.5*(b[:,2]+b[:,7])
    a.append(fit(t,x,5,t[-1])); w=[fit(t,x,p,q) for p,q in ((5,50),(50,100),(100,150),(150,200))]
    Dd = L(f'runaway_pair_d{d:02d}_L64_N128_lev0','sector_dynamics.dat'); tc=Dd[:,0]; xc=0.5*(Dd[:,1]+Dd[:,5])
    wc=[fit(tc,xc,p,q) for p,q in ((5,50),(50,100),(100,150),(150,200))]
    print(f'd={d}: bary windows ' + '/'.join(f'{v*1e4:.3f}' for v in w) + f'  half-spread {50*(max(w)-min(w))/a[-1]:.1f}%;  core windows ' + '/'.join(f'{v*1e4:.3f}' for v in wc) + f' half-spread {50*(max(wc)-min(wc))/np.mean(wc):.1f}%')
a=np.array(a); print('4-pt slope (no d=20):', np.polyfit(np.log(ds[:4]),np.log(a[:4]),1)[0])
for n in (128,256):
    D = L(f'runaway_pair_d10_L64_N{n}_lev0','sector_dynamics.dat'); t=D[:,0]
    sym=0.5*(D[:,14]+D[:,15]); xm=0.5*(D[:,1]+D[:,5]); v=np.gradient(xm,t)[-1]
    print(f'N={n}: max|mean shift| {np.abs(sym).max():.2e}, v_end {v:.3e}, int|sym|dt {np.trapezoid(np.abs(sym),t):.2e}, drift {xm[-1]-xm[0]:.4f}')
for c in ('runaway_pair_d10_L64_N128_lev0','wavezone_pair_d10_L128_N256_lev0','runaway_pair_d10_L64_N256_lev0'):
    b=L(c,'sector_barycenters.dat'); D=L(c,'sector_dynamics.dat')
    xb=0.5*(b[:,2]+b[:,7]); xc=0.5*(D[:,1]+D[:,5])
    print(f'{c}: bary drift {xb[-1]-xb[0]:+.4f} core drift {xc[-1]-xc[0]:+.4f} a_bary {fit(b[:,0],xb,5,200):.4e} a_core {fit(D[:,0],xc,5,200):.4e}')

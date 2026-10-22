"""Independent recomputation of headline numbers from the packed campaign data.
Reads only numeric .dat files. Usage: python3 -I recompute.py CAMPAIGN_DIR"""
import sys, numpy as np
C = sys.argv[1]
def L(cell, f): return np.loadtxt(f'{C}/{cell}/{f}', comments='#')
def bary(cell):
    d = L(cell, 'sector_barycenters.dat'); return d[:,0], d[:,2], d[:,7]
def dyn(cell): return L(cell, 'sector_dynamics.dat')
def fit(t, x, t0=5., t1=None):
    t1 = t[-1] if t1 is None else t1
    m = (t>=t0)&(t<=t1); c = np.polyfit(t[m]-t[m].mean(), x[m], 2); return 2*c[0]

M_P, M_M = 0.014349964, 0.014294692          # paper (5%-radius) masses
MT_P, MT_M = 0.0146265, 0.0145803             # true ADM (iso_star far-field, phase2 mass_check)
MB, MBT = 0.5*(M_P+M_M), 0.5*(MT_P+MT_M)
print(f'Mbar paper {MB:.7f}  true {MBT:.7f}  ratio paper/true {MB/MBT:.5f}')

HEAD = {n: f'runaway_pair_d10_L64_N{n}_lev0' for n in (128,192,256)}
print('\n== headline ladder (barycentre midpoint, t>=5) ==')
A = {}
for n, c in HEAD.items():
    t, xc, xp = bary(c); m = 0.5*(xc+xp); A[n] = fit(t, m)
    D = dyn(c); mc = 0.5*(D[:,1]+D[:,5]); acore = fit(D[:,0], mc)
    print(f' N={n}: drift {m[-1]-m[0]:+.4f} a_bary {A[n]:.4e} a_core {acore:.4e}  '
          f'a d^2/Mbar paper {A[n]*100/MB:.4f} true {A[n]*100/MBT:.4f}  bary gap {xc[0]-xp[0]:.4f}->{xc[-1]-xp[-1]:.4f}')
# Richardson with dx = 0.5, 1/3, 0.25
h = np.array([0.5, 1/3, 0.25]); a = np.array([A[128], A[192], A[256]])
p_obs = np.log(abs((a[1]-a[0])/(a[2]-a[1])))/np.log(1.5)  # crude (unequal ratios)
print(f' differences: 192-128 {a[1]-a[0]:+.3e}, 256-192 {a[2]-a[1]:+.3e} (sign change => non-monotone)')
for p in (2, 4):
    # fit a = a_inf + k h^p to the two finest
    k = (a[1]-a[2])/(h[1]**p-h[2]**p); ainf = a[2]-k*h[2]**p
    print(f'  Richardson p={p} (192,256): a_inf={ainf:.4e} -> a d^2/Mbar paper {ainf*100/MB:.4f} true {ainf*100/MBT:.4f}')

print('\n== separation scan (N=128, barycentre midpoint, t>=5) ==')
ds = np.array([8,10,12,16,20]); acc=[]; accc=[]
for d in ds:
    c = f'runaway_pair_d{d:02d}_L64_N128_lev0'
    t, xc, xp = bary(c); m = 0.5*(xc+xp); acc.append(fit(t, m))
    D = dyn(c); accc.append(fit(D[:,0], 0.5*(D[:,1]+D[:,5])))
    print(f' d={d:2d}: a_bary {acc[-1]:.4e} a_core {accc[-1]:.4e}  a d2/Mbar paper {acc[-1]*d*d/MB:.4f} true {acc[-1]*d*d/MBT:.4f}'
          f'  core-based paper {accc[-1]*d*d/MB:.4f}; gap {xc[0]-xp[0]:.3f}->{xc[-1]-xp[-1]:.3f} core sep {D[0,9]:.3f}->{D[-1,9]:.3f}')
acc = np.array(acc); accc=np.array(accc)
for lab, aa in (('bary', acc), ('core', accc)):
    s, o = np.polyfit(np.log(ds), np.log(aa), 1); res = np.log(aa)-(s*np.log(ds)+o)
    se = np.sqrt(np.sum(res**2)/3/np.sum((np.log(ds)-np.log(ds).mean())**2))
    s4 = np.polyfit(np.log(ds[1:]), np.log(aa[1:]), 1)[0]
    s3 = np.polyfit(np.log(ds[2:]), np.log(aa[2:]), 1)[0]
    print(f' {lab}: slope {s:.4f} +- {se:.4f}; drop d=8: {s4:.4f}; d>=12 only: {s3:.4f}')
# finite-size: fit a d^2 = M (1 + c/d^2)
for lab, aa in (('bary', acc),):
    X = np.vstack([np.ones(5), 1/ds**2]).T; co,*_ = np.linalg.lstsq(X, aa*ds**2, rcond=None)
    print(f' fit a d^2 = A(1+B/d^2): A={co[0]:.6f} (A/Mbar paper {co[0]/MB:.4f}, true {co[0]/MBT:.4f})')
# sensitivity of the slope to a d=10-like resolution correction (N256/N128 = r10) if it were d-dependent
r10 = A[256]/A[128]
print(f' N256/N128 at d=10: {r10:.4f}')

print('\n== window spread of the headline fit ==')
for n in (128,256):
    t, xc, xp = bary(HEAD[n]); m=0.5*(xc+xp)
    print(f' N={n}: ' + ' '.join(f'[{a0},{a1}]={fit(t,m,a0,a1):.4e}' for a0,a1 in ((5,50),(50,100),(100,150),(150,200),(5,100),(100,200))))

print('\n== mirror ==')
t0, xc0, xp0 = bary(HEAD[128]); tm, xcm, xpm = bary('control_mirror_mp_d10_L64_N128_lev0')
m0 = 0.5*(xc0+xp0); mm = 0.5*(xcm+xpm)
print(f' drift ratio {(mm[-1]-mm[0])/(m0[-1]-m0[0]):+.7f}  a ratio {fit(tm,mm)/fit(t0,m0):+.7f}')
# also sum of trajectories about centre: mirror symmetry pointwise
n = min(len(m0), len(mm)); dev = (m0[:n]-m0[0]) + (mm[:n]-mm[0])
print(f' max |drift_head + drift_mirror| {np.max(np.abs(dev)):.2e} vs drift {m0[-1]-m0[0]:.3f}  (rel {np.max(np.abs(dev))/abs(m0[-1]-m0[0]):.1e})')
# which lump is where in the mirror?
print(f'  head x_c0,x_p0 = {xc0[0]:.3f},{xp0[0]:.3f}; mirror x_c0,x_p0 = {xcm[0]:.3f},{xpm[0]:.3f}')

print('\n== signed matter momentum (sector_dynamics px columns) ==')
for n in (128,192,256):
    D = dyn(HEAD[n]); t = D[:,0]
    chk = np.max(np.abs(D[:,13]-(D[:,11]-D[:,12])))
    for tt in (50,100,200):
        i = np.argmin(abs(t-tt))
        print(f' N={n} t={t[i]:6.1f}: p+={D[i,11]:+.4e} p-raw={D[i,12]:+.4e} total={D[i,13]:+.3e} ({100*abs(D[i,13]/D[i,11]):.2f}% of p+)  [col check {chk:.1e}]')
    # momentum vs M v: v from core midpoint
    v = np.gradient(0.5*(D[:,1]+D[:,5]), t)
    i = np.argmin(abs(t-150))
    print(f'   p+/(M_P v) at t~150: paper-M {D[i,11]/(M_P*v[i]):.4f}  true-M {D[i,11]/(MT_P*v[i]):.4f};  p-/(|M_M| v): paper {D[i,12]/(M_M*v[i]):.4f} true {D[i,12]/(MT_M*v[i]):.4f}')

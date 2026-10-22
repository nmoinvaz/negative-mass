"""(1) verify mass_law.csv a / a_err / a_bary from data; (2) pull fits with alternative inputs to locate
the source of campaign/README's 0.597/0.739 and 0.0161/0.0119; (3) psi4 radial falloff on wavezone. Data only."""
import sys, numpy as np
C = sys.argv[1]
def L(c, f): return np.loadtxt(f'{C}/{c}/{f}', comments='#')
def q(t, x, a0, a1=None):
    a1 = t[-1] if a1 is None else a1; m=(t>=a0)&(t<=a1); return 2*np.polyfit(t[m], x[m], 2)[0]
for c in ('masslaw_pair_d20_w0850_L64_N128_lev0','runaway_pair_d20_L64_N128_lev0','masslaw_pair_d20_w0700_L64_N128_lev0','masslaw_pair_d20_w0675_L64_N128_lev0'):
    D=L(c,'sector_dynamics.dat'); t=D[:,0]; x=0.5*(D[:,1]+D[:,5]); B=L(c,'sector_barycenters.dat')
    w=[q(t,x,a,b) for a,b in ((5,50),(50,100),(100,150),(150,200))]
    print(f'{c[:36]}: a {q(t,x,5):.6e} a_err {0.5*(max(w)-min(w)):.4e} a_bary {q(B[:,0],0.5*(B[:,2]+B[:,7]),5):.6e}')
def cum(y,t): return np.r_[0,np.cumsum(0.5*(y[1:]+y[:-1])*np.diff(t))]
def pull(c, mode):
    D=L(c,'sector_dynamics.dat'); B=L(c,'sector_barycenters.dat')
    if mode=='core': t=D[:,0]; xc=D[:,1]; xp=D[:,5]; g=np.abs(D[:,9])
    else: t=B[:,0]; xc=B[:,2]; xp=B[:,7]; g=np.abs(xc-xp)
    I2=cum(cum(1/g**2,t),t); m=t>=5; A=np.vstack([np.ones(m.sum()),t[m],I2[m]]).T
    return {k:np.linalg.lstsq(A,x[m],rcond=None)[0][2] for k,x in (('c',xc),('p',xp))}
for mode in ('core','bary'):
    b=pull('runaway_pair_d10_L64_N128_lev0',mode); s=f'{mode}: matched C {b["c"]:.4f}/{b["p"]:.4f};'
    for k,c in (('w0804','massscale_pair_d10_w0804_L64_N128_lev0'),('w088','massratio_w088_r060_d10_L64_N128_lev0'),('heavy','massratio_heavyphantom_d10_L64_N128_lev0')):
        f=pull(c,mode); ch='p' if k=='heavy' else 'c'; ot='c' if ch=='p' else 'p'
        s+=f' {k} {f[ch]/b[ch]:.3f} (ctrl {f[ot]/b[ot]:.3f})'
    print(s)
d=L('wavezone_pair_d10_L128_N256_lev0','psi4_mode_l2_all.dat'); t=d[:,0]; R=np.array([16,24,32,40])
A=np.array([np.sqrt(sum(d[:,1+i*10+2*m]**2+d[:,2+i*10+2*m]**2 for m in range(5))) for i in range(4)])  # columns are r*psi4? see header
for tt in (50,100,150,190):
    j=np.argmin(abs(t-tt)); s=np.polyfit(np.log(R),np.log(A[:,j]),1)[0]; print(f'psi4 l2 t={t[j]:.0f}: amps '+' '.join(f'{v:.2e}' for v in A[:,j])+f' slope {s:.2f}, ratio R16/R40 {A[0,j]/A[3,j]:.1f}')

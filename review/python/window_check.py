"""Window dependence of the fitted acceleration (bary and core midpoints), and mass-law slope by tracker. Data only."""
import sys, math, numpy as np
C = sys.argv[1]
def L(c, f): return np.loadtxt(f'{C}/{c}/{f}', comments='#')
def fit(t, x, a0, a1):
    m = (t>=a0)&(t<=a1); return 2*np.polyfit(t[m]-t[m].mean(), x[m], 2)[0]
MBT = 0.0146034; MB = 0.0143223
for c in ['runaway_pair_d10_L64_N128_lev0','runaway_pair_d10_L64_N192_lev0','runaway_pair_d10_L64_N256_lev0','wavezone_pair_d10_L128_N256_lev0']:
    b = L(c,'sector_barycenters.dat'); D = L(c,'sector_dynamics.dat')
    for lab, t, x in (('bary', b[:,0], 0.5*(b[:,2]+b[:,7])), ('core', D[:,0], 0.5*(D[:,1]+D[:,5]))):
        edges = np.linspace(5, 200, 5)
        q = [fit(t, x, edges[i], edges[i+1]) for i in range(4)]
        h = [fit(t, x, 5, 102.5), fit(t, x, 102.5, 200)]
        # local acceleration from velocity slope: fit v(t) linearly in sliding 40-unit windows
        print(f'{c[:34]:34s} {lab}: quarters ' + '/'.join(f'{v*1e4:.3f}' for v in q) + f'  halves {h[0]*1e4:.3f}/{h[1]*1e4:.3f} e-4'
              f'  -> a d2/Mbar_true halves {h[0]*100/MBT:.3f}/{h[1]*100/MBT:.3f}')
# longrun windows quoted
b = L('longrun_pair_d10_t400_L64_N128_lev0','sector_barycenters.dat'); t=b[:,0]; x=0.5*(b[:,2]+b[:,7])
print('longrun bary: ' + ' '.join(f'[{a},{z}]={fit(t,x,a,z):.4e}' for a,z in ((5,400),(133,266),(300,400),(350,400),(5,200),(200,400))))
# mass ladder slope by tracker (paper csv masses and true masses)
cells = {'masslaw_pair_d20_w0850_L64_N128_lev0':(0.0091570754,0.0094372),'runaway_pair_d20_L64_N128_lev0':(0.0143216395,0.0146032),
         'masslaw_pair_d20_w0700_L64_N128_lev0':(0.0191728524,0.0194673),'masslaw_pair_d20_w0675_L64_N128_lev0':(0.0225412121,0.0228512)}
for lab in ('core','bary'):
    m1=[];m2=[];a=[]
    for c,(mp,mt) in cells.items():
        if lab=='core': D=L(c,'sector_dynamics.dat'); t=D[:,0]; x=0.5*(D[:,1]+D[:,5])
        else: B=L(c,'sector_barycenters.dat'); t=B[:,0]; x=0.5*(B[:,2]+B[:,7])
        a.append(fit(t,x,5,t[-1])); m1.append(mp); m2.append(mt)
    a=np.array(a)
    for nm, m in (('csv M',np.array(m1)),('true M',np.array(m2))):
        p=np.polyfit(np.log(m),np.log(a),1)[0]
        print(f'mass ladder {lab} {nm}: unweighted slope {p:.4f}; a d2/M = ' + ' '.join(f'{v:.4f}' for v in a*400/m))

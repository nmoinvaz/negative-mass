"""Time series probe of headline d=10 cells. Data only (np.loadtxt)."""
import sys, numpy as np
C = sys.argv[1]
def L(c,f): return np.loadtxt(f'{C}/{c}/{f}', comments='#')
cells = ['runaway_pair_d10_L64_N128_lev0','runaway_pair_d10_L64_N192_lev0','runaway_pair_d10_L64_N256_lev0',
         'longrun_pair_d10_t400_L64_N128_lev0','wavezone_pair_d10_L128_N256_lev0','deepsolve_pair_d10_L64_N128_lev0','gaugetwin_pair_d10_eta2_L64_N128_lev0','amrcheck_pair_d10_L64_N128_lev1']
for c in cells:
    D=L(c,'sector_dynamics.dat'); B=L(c,'sector_barycenters.dat'); K=L(c,'confinement.dat')
    print('==',c)
    for tt in (0,25,50,75,100,125,150,175,200,300,400):
        if tt> D[-1,0]+1: break
        i=np.argmin(abs(D[:,0]-tt)); j=np.argmin(abs(B[:,0]-tt))
        mid=0.5*(D[i,1]+D[i,5])
        print(f't={D[i,0]:6.1f} pk_c {D[i,4]:.5f} pk_p {D[i,8]:.5f} tot_c {B[j,1]:.4f} tot_p {B[j,6]:.4f} rms_c {B[j,5]:.3f} cfr_c {K[j,13]:.3f} mid {mid:8.4f} csep {D[i,9]:.4f} psep {D[i,10]:.4f} pxc {D[i,11]:.3e} pxp {D[i,12]:.3e} shc {D[i,14]:.2e} shp {D[i,15]:.2e} maxsh {D[i,16]:.2e}')

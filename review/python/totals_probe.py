"""Sector 'total' integrals, rms radii and core-gap / acceleration-vs-time for headline cells. Data only."""
import sys, numpy as np
C = sys.argv[1]
cells = ['runaway_pair_d10_L64_N128_lev0','runaway_pair_d10_L64_N192_lev0','runaway_pair_d10_L64_N256_lev0',
         'wavezone_pair_d10_L128_N256_lev0','longrun_pair_d10_t400_L64_N128_lev0','control_lone_canonical_L64_N128_lev0',
         'control_lone_phantom_L64_N128_lev0','runaway_pair_d20_L64_N128_lev0']
for c in cells:
    b = np.loadtxt(f'{C}/{c}/sector_barycenters.dat', comments='#'); t = b[:,0]
    D = np.loadtxt(f'{C}/{c}/sector_dynamics.dat', comments='#')
    j = [np.argmin(abs(t-x)) for x in (0,50,100,200)]
    print(c)
    print('   tot_c ' + ' '.join(f'{b[i,1]:.4f}' for i in j) + ' | tot_p ' + ' '.join(f'{b[i,6]:.4f}' for i in j)
          + ' | rms_c ' + ' '.join(f'{b[i,5]:.3f}' for i in j) + ' | rms_p ' + ' '.join(f'{b[i,10]:.3f}' for i in j))
    k = [np.argmin(abs(D[:,0]-x)) for x in (0,50,100,200)]
    print('   peak_c ' + ' '.join(f'{D[i,4]:.5f}' for i in k) + ' | peak_p ' + ' '.join(f'{D[i,8]:.5f}' for i in k) +
          ' | core sep ' + ' '.join(f'{D[i,9]:.4f}' for i in k) + ' | proper ' + ' '.join(f'{D[i,10]:.4f}' for i in k))

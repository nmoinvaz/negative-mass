"""(1) Symmetric (midpoint) vs antisymmetric ((a_c - a_p)/2) acceleration per quarter; (2) lone-star metric mass proxy
1 - chi (shell means at R=8,16) vs time, compared with omega- and content-implied growth. Data only; python3 -I.
usage: python3 -I sym_anti_chi.py <campaign_dir> <out_dir>"""
import sys, numpy as np
C, OUT = sys.argv[1], sys.argv[2]
def L(c, f): return np.loadtxt(f'{C}/{c}/{f}', comments='#')
def fit2(t, x, a0, a1):
    m = (t >= a0) & (t <= a1); return 2 * np.polyfit(t[m] - t[m].mean(), x[m], 2)[0]
def at(t, x, t0): return x[np.argmin(abs(t - t0))]
E = np.linspace(5, 200, 5)
out = ['grid,quarter,a_sym,a_anti,anti_over_sym']
print('### symmetric / antisymmetric parts (x1e-4)')
for g, c in {'N128': 'runaway_pair_d10_L64_N128_lev0', 'N192': 'runaway_pair_d10_L64_N192_lev0',
             'N256': 'runaway_pair_d10_L64_N256_lev0'}.items():
    D = L(c, 'sector_dynamics.dat'); t = D[:, 0]
    s = [];
    for i in range(4):
        ac, ap = fit2(t, D[:, 1], E[i], E[i + 1]), fit2(t, D[:, 5], E[i], E[i + 1])
        sy, an = 0.5 * (ac + ap), 0.5 * (ac - ap); s.append((sy, an))
        out.append(f'{g},{i+1},{sy:.6e},{an:.6e},{an/sy:.5f}')
    print(f'{g}: sym ' + ' '.join(f'{x[0]*1e4:.3f}' for x in s) + f'  (Q4/Q1 {s[3][0]/s[0][0]:.4f}) | anti ' +
          ' '.join(f'{x[1]*1e4:+.4f}' for x in s))
open(f'{OUT}/sym_anti.csv', 'w').write('\n'.join(out) + '\n')
print('\n### lone stars: 1-chi shell means (|chi-1| for phantom) at R8 and R16, t=40 and 200, and t=1000 for t1000 cell')
for g, c in {'lone_canon': 'control_lone_canonical_L64_N128_lev0', 'lone_phantom': 'control_lone_phantom_L64_N128_lev0',
             'lone_phantom_t1000': 'control_lone_phantom_t1000_L64_N128_lev0'}.items():
    S = L(c, 'shell_profiles.dat'); t = S[:, 0]
    r8 = np.abs(1 - S[:, 1]); r16 = np.abs(1 - S[:, 10])
    s = f'{g:20s} R8 {at(t,r8,40):.4e}->{at(t,r8,200):.4e} ({(at(t,r8,200)/at(t,r8,40)-1)*100:+.2f}%)  R16 {at(t,r16,40):.4e}->{at(t,r16,200):.4e} ({(at(t,r16,200)/at(t,r16,40)-1)*100:+.2f}%)'
    if t[-1] > 900: s += f'  t1000: R8 {(at(t,r8,1000)/at(t,r8,40)-1)*100:+.1f}% R16 {(at(t,r16,1000)/at(t,r16,40)-1)*100:+.1f}%'
    print(s)
print('note: canonical star is box-centred (shells centred on it); phantom cells are off-centre at x=37 (shells not centred).')

"""Per-star sliding accelerations (window 40) on the d=10 ladder + Richardson on the quarter rise. Data only; python3 -I.
usage: python3 -I per_star.py <campaign_dir> <out_dir>"""
import sys, numpy as np
C, OUT = sys.argv[1], sys.argv[2]
def L(c, f): return np.loadtxt(f'{C}/{c}/{f}', comments='#')
def fit2(t, x, a0, a1):
    m = (t >= a0) & (t <= a1); return 2 * np.polyfit(t[m] - t[m].mean(), x[m], 2)[0]
cells = {'N128': 'runaway_pair_d10_L64_N128_lev0', 'N192': 'runaway_pair_d10_L64_N192_lev0', 'N256': 'runaway_pair_d10_L64_N256_lev0'}
cen = np.arange(25, 181, 5)
hdr = ['t_center'] + [f'{g}_{s}' for g in cells for s in ('canon', 'phantom', 'mid')]
rows = [[f'{x}'] for x in cen]
for g, c in cells.items():
    D = L(c, 'sector_dynamics.dat'); t = D[:, 0]
    for k, tc in enumerate(cen):
        ac, ap = fit2(t, D[:, 1], tc - 20, tc + 20), fit2(t, D[:, 5], tc - 20, tc + 20)
        rows[k] += [f'{ac:.6e}', f'{ap:.6e}', f'{0.5*(ac+ap):.6e}']
open(f'{OUT}/per_star_sliding.csv', 'w').write(','.join(hdr) + '\n' + '\n'.join(','.join(r) for r in rows) + '\n')
E = np.linspace(5, 200, 5); R = {}
for g, c in cells.items():
    D = L(c, 'sector_dynamics.dat'); t = D[:, 0]; mid = 0.5 * (D[:, 1] + D[:, 5])
    q = [fit2(t, mid, E[i], E[i + 1]) for i in range(4)]; R[g] = (q[0], q[3], q[3] / q[0])
    print(f'{g}: Q1 {q[0]*1e4:.4f} Q4 {q[3]*1e4:.4f} Q4/Q1 {q[3]/q[0]:.4f}')
for p in (1, 2, 3, 4):
    f = (256 / 192) ** p
    rr = R['N256'][2] + (R['N256'][2] - R['N192'][2]) / (f - 1)
    q4 = R['N256'][1] + (R['N256'][1] - R['N192'][1]) / (f - 1)
    print(f'Richardson N192/N256, p={p}: Q4/Q1 -> {rr:.4f}; Q4 a -> {q4*1e4:.4f}e-4 (a d2/Mtrue at d=10: {q4*100/0.0146034:.4f})')
# order needed for the rise to vanish: (r192-1)/(r256-1) = (256/192)^p
r1, r2 = R['N192'][2] - 1, R['N256'][2] - 1
print(f'order p at which the rise would extrapolate to zero: {np.log(r1/r2)/np.log(256/192):.2f}')

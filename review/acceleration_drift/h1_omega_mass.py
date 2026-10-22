"""H1 with an omega-based mass proxy: M(t) = M_true + dM/domega * (omega(t) - omega(early)), omega from zero crossings of
(max_phi+min_phi) and of (max_Pi+min_Pi) in collapse_diagnostics.dat (canonical sector only). Data only; python3 -I.
usage: python3 -I h1_omega_mass.py <campaign_dir> <out_dir>"""
import sys, numpy as np
C, OUT = sys.argv[1], sys.argv[2]
def L(c, f): return np.loadtxt(f'{C}/{c}/{f}', comments='#')
def fit2(t, x, a0, a1):
    m = (t >= a0) & (t <= a1); return 2 * np.polyfit(t[m] - t[m].mean(), x[m], 2)[0]
dMdw = -0.0661; MT = 0.0146265  # canonical family slope near omega=0.75 (phase2 boson_star_output.txt)
def omega_series(c, cols):
    K = L(c, 'collapse_diagnostics.dat'); t = K[:, 0]; s = K[:, cols[0]] + K[:, cols[1]]
    i = np.where(np.sign(s[:-1]) * np.sign(s[1:]) < 0)[0]
    tz = t[i] - s[i] * (t[i + 1] - t[i]) / (s[i + 1] - s[i])
    return tz
def omega_in(tz, a0, a1):
    m = (tz >= a0) & (tz < a1); return np.pi * np.polyfit(tz[m], np.arange(m.sum()), 1)[0]
cells = {'N128': 'runaway_pair_d10_L64_N128_lev0', 'N192': 'runaway_pair_d10_L64_N192_lev0',
         'N256': 'runaway_pair_d10_L64_N256_lev0', 'L128N256(dx0.5)': 'wavezone_pair_d10_L128_N256_lev0',
         'lone_c_N128': 'control_lone_canonical_L64_N128_lev0',
         'stab_w075': 'stability_canonical_w075_L64_N128_lev0', 'stab_w080': 'stability_canonical_w080_L64_N128_lev0',
         'stab_w085': 'stability_canonical_w085_L64_N128_lev0', 'stab_w090': 'stability_canonical_w090_L64_N128_lev0'}
edges = np.linspace(5, 200, 5)
out = ['grid,quarter,omega_phi,omega_Pi,M_proxy_rel,a_mid,a_over_Mproxy_rel']
for g, c in cells.items():
    tzf, tzp = omega_series(c, (10, 11)), omega_series(c, (12, 13))
    T = min(tzf[-1], tzp[-1])
    E = edges if T > 190 else np.linspace(5, 115, 3)
    w = [(omega_in(tzf, E[i], E[i + 1]), omega_in(tzp, E[i], E[i + 1])) for i in range(len(E) - 1)]
    s = f'{g:16s} omega(phi/Pi) ' + ' '.join(f'{a:.4f}/{b:.4f}' for a, b in w)
    if g in ('N128', 'N192', 'N256', 'L128N256(dx0.5)'):
        D = L(c, 'sector_dynamics.dat'); t = D[:, 0]; mid = 0.5 * (D[:, 1] + D[:, 5])
        a = [fit2(t, mid, E[i], E[i + 1]) for i in range(4)]
        Mr = [1 + dMdw * (0.5 * (x[0] + x[1]) - 0.5 * (w[0][0] + w[0][1])) / MT for x in w]
        r = [a[i] / a[0] / Mr[i] for i in range(4)]
        s += ' | M/M_Q1 ' + ' '.join(f'{x:.4f}' for x in Mr) + ' | a/a_Q1 ' + ' '.join(f'{a[i]/a[0]:.4f}' for i in range(4)) + \
             ' | (a/M)/(a/M)_Q1 ' + ' '.join(f'{x:.4f}' for x in r)
        for i in range(4):
            out.append(f'{g},{i+1},{w[i][0]:.6f},{w[i][1]:.6f},{Mr[i]:.5f},{a[i]:.6e},{r[i]:.5f}')
    print(s)
open(f'{OUT}/h1_omega_mass.csv', 'w').write('\n'.join(out) + '\n')

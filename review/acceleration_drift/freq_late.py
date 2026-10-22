"""(a) Canonical-star oscillation frequency vs time from collapse_diagnostics max_phi+min_phi (proxy for mass via M(omega));
(b) late-time cells: sliding midpoint/individual accelerations, sponge distance. Data only; python3 -I.
usage: python3 -I freq_late.py <campaign_dir> <out_dir>"""
import sys, numpy as np
C, OUT = sys.argv[1], sys.argv[2]
def L(c, f): return np.loadtxt(f'{C}/{c}/{f}', comments='#')
def fit2(t, x, a0, a1):
    m = (t >= a0) & (t <= a1); return 2 * np.polyfit(t[m] - t[m].mean(), x[m], 2)[0]
# M(omega_phys) canonical family from phase2/boson_star_output.txt (our solver)
W = np.array([0.6916, 0.7426, 0.7500, 0.7934, 0.8441]); M = np.array([0.020506, 0.015218, 0.014626, 0.011823, 0.009624])
dMdw = np.polyfit(W[1:4], M[1:4], 1)[0]
print(f'local dM/domega near 0.75 = {dMdw:.4f}  -> 11% mass rise needs domega = {0.11*0.014626/dMdw:+.4f}')
cells = {'N128': 'runaway_pair_d10_L64_N128_lev0', 'N192': 'runaway_pair_d10_L64_N192_lev0',
         'N256': 'runaway_pair_d10_L64_N256_lev0', 'N128_t400': 'longrun_pair_d10_t400_L64_N128_lev0',
         'lone_c_N128': 'control_lone_canonical_L64_N128_lev0'}
out = ['grid,t0,t1,omega_coord,n_crossings']
for g, c in cells.items():
    K = L(c, 'collapse_diagnostics.dat'); t = K[:, 0]; s = K[:, 10] + K[:, 11]
    # zero crossings by linear interpolation
    i = np.where(np.sign(s[:-1]) * np.sign(s[1:]) < 0)[0]
    tz = t[i] - s[i] * (t[i + 1] - t[i]) / (s[i + 1] - s[i])
    res = []
    for a0 in np.arange(0, t[-1] - 49, 50):
        m = (tz >= a0) & (tz < a0 + 50)
        if m.sum() < 4: continue
        k = np.arange(m.sum()); p = np.polyfit(tz[m], k, 1)[0]  # crossings per time
        om = np.pi * p; res.append(om); out.append(f'{g},{a0},{a0+50},{om:.6f},{m.sum()}')
    print(f'{g:12s} coordinate omega per 50-window: ' + ' '.join(f'{x:.5f}' for x in res) +
          f'  | change first->last {res[-1]-res[0]:+.5f}  => implied dM/M {(res[-1]-res[0])*dMdw/0.014626*100:+.2f}%')
open(f'{OUT}/omega_windows.csv', 'w').write('\n'.join(out) + '\n')

print('\n### late-time cells: 50-unit windows of a_mid, a_canon, a_phantom (x1e-4), core sep, distance of leading star to sponge (r=24 from 32)')
late = {'N128_t400': 'longrun_pair_d10_t400_L64_N128_lev0', 'nomill_left': 'nomill_left_pair_d10_L64_N128_lev0',
        'chase03c(treadmill)': 'chase03c_pair_d10_L64_N128_lev0'}
out = ['cell,t0,t1,a_mid,a_canon,a_phantom,core_sep,lead_r_from_centre,trail_r_from_centre']
for g, c in late.items():
    D = L(c, 'sector_dynamics.dat'); t = D[:, 0]; xc, xp = D[:, 1], D[:, 5]; mid = 0.5 * (xc + xp)
    if 'chase' in g:
        # treadmill shifts the grid; use treadmill midpoint_true if available
        T = L(c, 'treadmill.dat')
        print('  chase03c: grid recentred; coordinates jump at shifts -> only gap and midpoint_true used')
        tt, mt = T[:, 0], T[:, 8]
        for a0 in np.arange(5, tt[-1] - 49, 50):
            a = fit2(tt, mt, a0, a0 + 50)
            m = (t >= a0) & (t <= a0 + 50)
            print(f'  {g:20s} [{a0:5.0f},{a0+50:5.0f}] a_mid_true {a*1e4:.3f}  sep {D[m,9].mean():.3f}')
            out.append(f'{g},{a0},{a0+50},{a:.6e},,,{D[m,9].mean():.5f},,')
        continue
    for a0 in np.arange(5, t[-1] - 49, 50):
        m = (t >= a0) & (t <= a0 + 50)
        am, ac, ap = fit2(t, mid, a0, a0 + 50), fit2(t, xc, a0, a0 + 50), fit2(t, xp, a0, a0 + 50)
        lead = np.abs(xc[m] - 32).max(); trail = np.abs(xp[m] - 32).max()
        print(f'  {g:20s} [{a0:5.0f},{a0+50:5.0f}] mid {am*1e4:.3f} can {ac*1e4:.3f} ph {ap*1e4:.3f} sep {D[m,9].mean():.3f} '
              f'lead r {lead:.2f} trail r {trail:.2f}')
        out.append(f'{g},{a0},{a0+50},{am:.6e},{ac:.6e},{ap:.6e},{D[m,9].mean():.5f},{lead:.3f},{trail:.3f}')
open(f'{OUT}/late_windows.csv', 'w').write('\n'.join(out) + '\n')

# Same-time comparison longrun vs nomill (identical physics, pair offset by -10): differences show box/sponge influence
D1 = L(late['N128_t400'], 'sector_dynamics.dat'); D2 = L(late['nomill_left'], 'sector_dynamics.dat')
print('\n### longrun vs nomill (pair 10 left): gap and midpoint displacement at same t')
for t0 in (100, 200, 300, 400):
    i = np.argmin(abs(D1[:, 0] - t0)); j = np.argmin(abs(D2[:, 0] - t0))
    d1 = 0.5 * (D1[i, 1] + D1[i, 5]) - 32.0002; d2 = 0.5 * (D2[j, 1] + D2[j, 5]) - 22.0002
    print(f'  t={t0}: drift longrun {d1:.4f} nomill {d2:.4f} (ratio {d2/d1:.4f}); sep {D1[i,9]:.4f} vs {D2[j,9]:.4f}')

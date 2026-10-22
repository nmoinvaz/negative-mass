"""Acceleration-drift hypotheses H1-H4 on the d=10 cells. Data only (np.loadtxt); run with python3 -I.
usage: python3 -I drift_analysis.py <campaign_dir> <out_dir>
"""
import sys, numpy as np
C, OUT = sys.argv[1], sys.argv[2]
def L(c, f): return np.loadtxt(f'{C}/{c}/{f}', comments='#')
MBT = 0.0146034  # true ADM mean mass
def fit2(t, x, a0, a1):
    m = (t >= a0) & (t <= a1)
    tt = t[m] - t[m].mean(); p = np.polyfit(tt, x[m], 2)
    return 2 * p[0]
def fit1(t, x, a0, a1):
    m = (t >= a0) & (t <= a1); return np.polyfit(t[m] - t[m].mean(), x[m], 1)[0]
def at(t, x, t0): return x[np.argmin(abs(t - t0))]

cells = {'N128': 'runaway_pair_d10_L64_N128_lev0', 'N192': 'runaway_pair_d10_L64_N192_lev0',
         'N256': 'runaway_pair_d10_L64_N256_lev0', 'N128_t400': 'longrun_pair_d10_t400_L64_N128_lev0',
         'L128N256(dx0.5)': 'wavezone_pair_d10_L128_N256_lev0', 'N128_eta2': 'gaugetwin_pair_d10_eta2_L64_N128_lev0',
         'N128_amr1': 'amrcheck_pair_d10_L64_N128_lev1'}
lone = {'lone_c_N128': 'control_lone_canonical_L64_N128_lev0', 'lone_p_N128': 'control_lone_phantom_L64_N128_lev0'}
rows = []
print('### Quarter fits on [5,200] (x1e-4): midpoint a, canonical a, phantom a; a*d^2/Mbar_true with instantaneous d')
edges = np.linspace(5, 200, 5)
csv_q = ['grid,quarter,t0,t1,a_mid,a_canon,a_phantom,a_mid_bary,dbar_core,dbar_proper,a_d2_M_fixed10,a_d2_M_core,a_d2_M_proper,dpdt_c,dpdt_p,meff_c,meff_p,content_c,content_p,peak_c,peak_p']
for g, c in cells.items():
    D = L(c, 'sector_dynamics.dat'); B = L(c, 'sector_barycenters.dat')
    t = D[:, 0]; xc = D[:, 1]; xp = D[:, 5]; mid = 0.5 * (xc + xp)
    tb = B[:, 0]; midb = 0.5 * (B[:, 2] + B[:, 7])
    T = t[-1]
    eds = edges if T < 250 else np.linspace(5, 400, 9) if T < 450 else edges
    for i in range(len(eds) - 1):
        a0, a1 = eds[i], eds[i + 1]
        m = (t >= a0) & (t <= a1)
        am, ac, ap = fit2(t, mid, a0, a1), fit2(t, xc, a0, a1), fit2(t, xp, a0, a1)
        ab = fit2(tb, midb, a0, a1)
        dco = D[m, 9].mean(); dpr = D[m, 10].mean()
        # matter-momentum rate and effective inertial mass p/v (v from linear fit of x in window)
        dpc, dpp = fit1(t, D[:, 11], a0, a1), fit1(t, D[:, 12], a0, a1)
        vc, vp = fit1(t, xc, a0, a1), fit1(t, xp, a0, a1)
        mc, mp = D[m, 11].mean() / vc, D[m, 12].mean() / vp
        mb = (tb >= a0) & (tb <= a1)
        cc, cp = B[mb, 1].mean(), B[mb, 6].mean()
        pc, pp = D[m, 4].mean(), D[m, 8].mean()
        csv_q.append(f'{g},{i+1},{a0:.1f},{a1:.1f},{am:.6e},{ac:.6e},{ap:.6e},{ab:.6e},{dco:.5f},{dpr:.5f},'
                     f'{am*100/MBT:.5f},{am*dco**2/MBT:.5f},{am*dpr**2/MBT:.5f},{dpc:.6e},{dpp:.6e},{mc:.6e},{mp:.6e},{cc:.5f},{cp:.5f},{pc:.6e},{pp:.6e}')
        print(f'{g:16s} [{a0:5.1f},{a1:5.1f}] mid {am*1e4:.3f} can {ac*1e4:.3f} ph {ap*1e4:.3f} bary {ab*1e4:.3f} | '
              f'd_core {dco:.4f} d_prop {dpr:.4f} | ad2/M (10,core,prop) {am*100/MBT:.4f} {am*dco**2/MBT:.4f} {am*dpr**2/MBT:.4f} | '
              f'dp/dt c,p {dpc*1e4:.3f} {dpp*1e4:.3f} | p/v c,p {mc:.5f} {mp:.5f} | content c,p {cc:.3f} {cp:.3f} | peak {pc:.5f} {pp:.5f}')
open(f'{OUT}/quarters.csv', 'w').write('\n'.join(csv_q) + '\n')

# Sliding local acceleration (window 40) for plotting
print('\n### sliding a_mid(t), window 40 (x1e-4)')
lines = ['t_center,' + ','.join(cells)]
centers = np.arange(25, 381, 5)
series = {}
for g, c in cells.items():
    D = L(c, 'sector_dynamics.dat'); t = D[:, 0]; mid = 0.5 * (D[:, 1] + D[:, 5])
    series[g] = [fit2(t, mid, tc - 20, tc + 20) if tc + 20 <= t[-1] + 0.1 else np.nan for tc in centers]
for k, tc in enumerate(centers):
    lines.append(f'{tc}' + ''.join(f',{series[g][k]:.6e}' for g in cells))
open(f'{OUT}/a_mid_sliding.csv', 'w').write('\n'.join(lines) + '\n')
for k, tc in enumerate(centers):
    if tc % 25 == 0: print(f't={tc:4d} ' + ' '.join(f'{g}:{series[g][k]*1e4:.3f}' for g in cells))
# linear trend of the sliding a in [25,180]
print('\n### linear trend of sliding a_mid over centres 25..180: slope (per 100 time) / mean')
for g in cells:
    s = np.array(series[g]); m = (centers >= 25) & (centers <= 180)
    p = np.polyfit(centers[m], s[m], 1); print(f'{g:16s} mean {s[m].mean()*1e4:.4f}e-4  rel slope {p[0]*100/s[m].mean()*100:+.2f}% per 100')

# H2: gauge numbers
print('\n### H2 gauge: v_mid(200), gamma^3-1, shift at cores, mean core shift, lapse')
for g, c in cells.items():
    D = L(c, 'sector_dynamics.dat'); t = D[:, 0]; mid = 0.5 * (D[:, 1] + D[:, 5])
    K = L(c, 'collapse_diagnostics.dat')
    v = fit1(t, mid, 190, 200) if t[-1] >= 200 else np.nan
    sm = 0.5 * (D[:, 14] + D[:, 15])
    i2 = np.argmin(abs(t - 200))
    dsm = fit1(t, sm, 100, 200)
    print(f'{g:16s} v(200)={v:.4f} gamma^3-1={1.5*v*v*100:.3f}% | shift_c {D[i2,14]:+.2e} shift_p {D[i2,15]:+.2e} mean {sm[i2]:+.2e} '
          f'd(mean shift)/dt {dsm:+.2e} (= {dsm/1.5e-4*100:+.2f}% of a) | min lapse t=5,200: {at(K[:,0],K[:,1],5):.5f} {at(K[:,0],K[:,1],200):.5f}')

# H1: content growth vs lone stars
print('\n### H1 field content (int |Phi|-modulus dV) and peak; lone-star controls')
for g, c in {**cells, **lone}.items():
    B = L(c, 'sector_barycenters.dat'); D = L(c, 'sector_dynamics.dat'); K = L(c, 'confinement.dat')
    tb = B[:, 0]; t = D[:, 0]
    s = ''
    for nm, col, colD, colK in (('c', 1, 4, 13), ('p', 6, 8, 15)):
        x0, x2 = at(tb, B[:, col], 0), at(tb, B[:, col], 200)
        if not np.isfinite(x0) or x0 == 0: continue
        p0, p2 = at(t, D[:, colD], 0), at(t, D[:, colD], 200)
        f0, f2 = at(K[:, 0], K[:, colK], 0), at(K[:, 0], K[:, colK], 200)
        s += (f' {nm}: content {x0:.3f}->{x2:.3f} ({(x2/x0-1)*100:+.1f}%) confined-content {(x2*f2/(x0*f0)-1)*100:+.1f}% '
              f'peak {(p2/p0-1)*100:+.2f}%')
    print(f'{g:16s}{s}')

# H4: positions vs sponge
print('\n### H4 geometry: core |x - 32| and leading-edge distance to sponge inner radius 24 (L64 cells)')
for g in ('N128', 'N192', 'N256', 'N128_t400'):
    D = L(cells[g], 'sector_dynamics.dat'); t = D[:, 0]
    for t0 in (0, 100, 200, 300, 400):
        if t0 > t[-1] + 0.1: continue
        xc, xp = at(t, D[:, 1], t0), at(t, D[:, 5], t0)
        print(f'{g:10s} t={t0:3d} canon x={xc:.3f} (r={abs(xc-32):.2f}, gap to sponge {24-abs(xc-32):.2f})  phantom x={xp:.3f} (r={abs(xp-32):.2f})')
# constraint growth
print('\n### H4 constraints L2_Ham at t=5,100,200 and growth')
for g, c in cells.items():
    Q = L(c, 'constraint_norms.dat')
    h = [at(Q[:, 0], Q[:, 1], x) for x in (5, 100, 200)]
    print(f'{g:16s} Ham {h[0]:.3e} {h[1]:.3e} {h[2]:.3e}  ratio200/5 {h[2]/h[0]:.2f}; Mom {at(Q[:,0],Q[:,2],5):.2e} {at(Q[:,0],Q[:,2],200):.2e}')

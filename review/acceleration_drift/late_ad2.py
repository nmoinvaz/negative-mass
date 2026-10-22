"""Late-time a*d^2/M_true against the omega mass proxy for longrun, nomill_left, chase03c (N128). Data only; python3 -I.
usage: python3 -I late_ad2.py <campaign_dir> <out_dir>"""
import sys, numpy as np
C, OUT = sys.argv[1], sys.argv[2]
def L(c, f): return np.loadtxt(f'{C}/{c}/{f}', comments='#')
def fit2(t, x, a0, a1):
    m = (t >= a0) & (t <= a1); return 2 * np.polyfit(t[m] - t[m].mean(), x[m], 2)[0]
MBT, dMdw, MT = 0.0146034, -0.0661, 0.0146265
def tz_of(c):
    K = L(c, 'collapse_diagnostics.dat'); t = K[:, 0]; s = K[:, 10] + K[:, 11]
    i = np.where(np.sign(s[:-1]) * np.sign(s[1:]) < 0)[0]
    return t[i] - s[i] * (t[i + 1] - t[i]) / (s[i + 1] - s[i])
def om(tz, a0, a1):
    m = (tz >= a0) & (tz < a1); return np.pi * np.polyfit(tz[m], np.arange(m.sum()), 1)[0]
out = ['cell,t0,t1,a_mid,d_core,a_d2_over_Mtrue,omega,Mproxy_rel,a_d2_over_Mproxy']
for g, c in {'longrun': 'longrun_pair_d10_t400_L64_N128_lev0', 'nomill_left': 'nomill_left_pair_d10_L64_N128_lev0',
             'chase03c': 'chase03c_pair_d10_L64_N128_lev0'}.items():
    D = L(c, 'sector_dynamics.dat'); t = D[:, 0]
    if g == 'chase03c':
        T = L(c, 'treadmill.dat'); tm, mid = T[:, 0], T[:, 8]
    else:
        tm, mid = t, 0.5 * (D[:, 1] + D[:, 5])
    tz = tz_of(c); w0 = om(tz, 5, 55)
    for a0 in np.arange(5, min(tm[-1], tz[-1]) - 49, 50):
        a = fit2(tm, mid, a0, a0 + 50); m = (t >= a0) & (t <= a0 + 50); d = D[m, 9].mean()
        w = om(tz, a0, a0 + 50); Mr = 1 + dMdw * (w - w0) / MT
        r = a * d * d / MBT
        out.append(f'{g},{a0},{a0+50},{a:.6e},{d:.4f},{r:.5f},{w:.5f},{Mr:.4f},{r/Mr:.5f}')
        print(f'{g:12s} [{a0:4.0f},{a0+50:4.0f}] a {a*1e4:.3f} d {d:7.3f} a d2/Mtrue {r:.4f} omega {w:.4f} Mproxy {Mr:.3f} a d2/Mproxy {r/Mr:.4f}')
open(f'{OUT}/late_ad2.csv', 'w').write('\n'.join(out) + '\n')

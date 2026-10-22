"""Seam test: compare paper-style seeded star (Yukawa tail + Schwarzschild lapse with the 5%-radius
mass beyond r_m, lapse jump) with the exact BVP seed, same grid.  Reports constraint-solved initial
m, N; omega; early breathing amplitude; and the long-time drift; plus a prediction of the omega offset
from the N offset via dN/domega of the static family."""
import json, numpy as np
from analyze import load
fam = json.load(open('family_slopes.json'))
def om(d, w):
    s = (d['t'] >= w[0]) & (d['t'] <= w[1]); return float(-np.polyfit(d['t'][s], np.unwrap(d['phase'])[s], 1)[0])
for ex, se, se0, key in [('c750_h0.1', 'c750_seam_h0.1', 'c750_seam0_h0.1', 'c0750'),
                         ('p7603_h0.1', 'p7603_seam_h0.1', 'p7603_seam0_h0.1', 'p07603'),
                         ('c750_h0.05', 'c750_seam_h0.05', None, 'c0750'),
                         ('p7603_h0.05', 'p7603_seam_h0.05', None, 'p07603')]:
    try:
        _, a = load(ex); _, b = load(se)
    except Exception as e:
        print('missing', ex, se); continue
    if b['t'][-1] < 1499: print('unfinished', se); continue
    f = fam[key]
    dN = b['N50'][0] / a['N50'][0] - 1; dM = b['m50'][0] / a['m50'][0] - 1
    dw_pred = dN / f['dlnNdw']
    w_ex = om(a, (500, 1500)); w_se = om(b, (500, 1500))
    def breathe(d, t0, t1):
        s = (d['t'] >= t0) & (d['t'] <= t1); x = d['amp'][s]; return float((x.max() - x.min()) / 2 / d['amp'][0])
    print(f"{se}: seeded M(50) rel to exact {dM:+.2e}, N {dN:+.2e}  (exact-seed M={a['m50'][0]:+.8f}, seam M={b['m50'][0]:+.8f})")
    print(f"   omega: exact {w_ex:.8f}, seam {w_se:.8f}, diff {w_se-w_ex:+.2e}; predicted from dN via family {dw_pred:+.2e}")
    print(f"   breathing half-amplitude of |Phi_c|: t 0-100 exact {breathe(a,0,100):.1e} seam {breathe(b,0,100):.1e}; t 1400-1500 exact {breathe(a,1400,1500):.1e} seam {breathe(b,1400,1500):.1e}")
    print(f"   drift [0,100]->[1400,1500]: omega {om(b,(1400,1500))-om(b,(0,100)):+.2e}, amp {np.mean(b['amp'][b['t']>1400])/np.mean(b['amp'][b['t']<100])-1:+.2e}, m50 {b['m50'][-1]/b['m50'][0]-1:+.2e}, N50 {b['N50'][-1]/b['N50'][0]-1:+.2e}")
    if se0:
        try:
            _, c = load(se0)
            if c['t'][-1] >= 1499:
                print(f"   no-jump variant: omega {om(c,(500,1500)):.8f} (diff to exact {om(c,(500,1500))-w_ex:+.2e}), breathing t<100 {breathe(c,0,100):.1e}")
        except Exception:
            pass

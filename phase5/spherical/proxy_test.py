"""Test the 'frequency drift implies mass drift' proxy inside our own evolutions:
predicted Delta M/M = (dlnM/domega)_family * Delta omega  vs measured Misner-Sharp Delta m(50)/m,
using windows [0,100] -> [1400,1500] (drift runs) or [200,300] -> [1400,1500] (perturbed runs, after
the initial transient).  Also writes mass_vs_time.csv (Misner-Sharp m at several radii, N, omega(t),
central amplitude for the production runs)."""
import json, numpy as np
from analyze import load, omega_series
fam = json.load(open('family_slopes.json'))
def om(d, w):
    s = (d['t'] >= w[0]) & (d['t'] <= w[1]); return float(-np.polyfit(d['t'][s], np.unwrap(d['phase'])[s], 1)[0])
def mean(d, c, w):
    s = (d['t'] >= w[0]) & (d['t'] <= w[1]); return float(np.mean(d[c][s]))
runs = ['c750_h0.5', 'c750_h0.5_ko1', 'c750_h0.333', 'c750_h0.25', 'c750_h0.2', 'p7603_h0.5', 'p7603_h0.5_ko1',
        'p7603_h0.333', 'c750_seam_h0.1', 'p7603_seam_h0.1', 'c750_scale3_h0.1', 'p7603_scale3_h0.1',
        'c750_scalem3_h0.1', 'p7603_scalem3_h0.1']
print('| run | window | Δω | predicted ΔM/M | measured Δm(50)/m | measured Δm(20)/m | ratio |')
print('|---|---|---|---|---|---|---|')
for r in runs:
    try:
        h, d = load(r)
    except Exception:
        continue
    if d['t'][-1] < 1499: continue
    W0 = (0, 100) if h['kind'] == 'exact' else (200, 300)
    W1 = (1400, 1500)
    fs = fam['c0750' if h['eps'] > 0 else 'p07603']
    dw = om(d, W1) - om(d, W0)
    pred = fs['dlnMdw'] * dw
    m50 = mean(d, 'm50', W1) / mean(d, 'm50', W0) - 1
    m20 = mean(d, 'm20', W1) / mean(d, 'm20', W0) - 1
    print(f"| {r} | {W0}→{W1} | {dw:+.2e} | {pred:+.2e} | {m50:+.2e} | {m20:+.2e} | {m50/pred if pred else float('nan'):.3f} |")
# time series export
prod = ['c750_h0.1', 'p7603_h0.1', 'c750_h0.05', 'p7603_h0.05', 'c750_h0.5', 'p7603_h0.5', 'c750_seam_h0.1',
        'p7603_seam_h0.1']
with open('mass_vs_time.csv', 'w') as f:
    f.write('run,t,omega_50win,amp,m8,m12,m20,m30,m50,mR,N20,N50\n')
    for r in prod:
        try:
            h, d = load(r)
        except Exception:
            continue
        tc, wc = omega_series(d, 50)
        for i in range(0, len(d), 10):
            t = d['t'][i]
            w = np.interp(t, tc, wc) if tc[0] <= t <= tc[-1] else np.nan
            f.write(f"{r},{t:.1f},{w:.9f},{d['amp'][i]:.10e}," + ','.join(f"{d[c][i]:.11e}" for c in
                    ('m8', 'm12', 'm20', 'm30', 'm50', 'mR', 'N20', 'N50')) + '\n')

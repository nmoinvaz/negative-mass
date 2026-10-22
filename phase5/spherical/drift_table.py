"""Drift table for the lone campaign stars: per run, omega and normalised amp/m/N/content in windows,
drifts between [0,100] and [1400,1500], and observed convergence orders of the drifts.
Writes drift_table.csv and prints a markdown table."""
import sys, numpy as np
from analyze import load

W0, W1 = (0, 100), (1400, 1500)

def win(d, col, w):
    s = (d['t'] >= w[0]) & (d['t'] <= w[1])
    return float(np.mean(d[col][s]))

def omega(d, w):
    s = (d['t'] >= w[0]) & (d['t'] <= w[1])
    return float(-np.polyfit(d['t'][s], np.unwrap(d['phase'])[s], 1)[0])

def row(name):
    h, d = load(name)
    r = dict(name=name, h=h['h'], R=h['R'], sko=h['sko'], kind=h['kind'], w_star=h['w'], M=h['M_bvp'], N=h['N_bvp'])
    r['w0'] = omega(d, W0); r['w1'] = omega(d, W1)
    r['dw'] = r['w1'] - r['w0']
    r['w0-wstar'] = r['w0'] - h['w']
    for c in ('amp', 'm20', 'm50', 'N50', 'content20'):
        a, b = win(d, c, W0), win(d, c, W1)
        r[c + '_0'] = a; r['d' + c + '_rel'] = b / a - 1
    r['m50_0/M-1'] = r['m50_0'] / h['M_bvp'] - 1
    r['N50_0/N-1'] = r['N50_0'] / h['N_bvp'] - 1
    r['amp0/phic-1'] = d['amp'][0] / 1.0
    return r

if __name__ == '__main__':
    names = sys.argv[1:]
    rows = [row(n) for n in names]
    keys = list(rows[0].keys())
    with open('drift_table.csv', 'w') as f:
        f.write(','.join(keys) + '\n')
        for r in rows:
            f.write(','.join(str(r[k]) for k in keys) + '\n')
    print('| run | h | ω[0,100] − ω_star | Δω (t 50→1450) | Δamp/amp | Δm(50)/m | ΔN(50)/N | m(50)/M_ADM − 1 |')
    print('|---|---|---|---|---|---|---|---|')
    for r in rows:
        print(f"| {r['name']} | {r['h']:.3g} | {r['w0-wstar']:+.2e} | {r['dw']:+.2e} | {r['damp_rel']:+.2e} | "
              f"{r['dm50_rel']:+.2e} | {r['dN50_rel']:+.2e} | {r['m50_0/M-1']:+.2e} |")
    # convergence orders for successive pairs with the same star/kind
    print()
    for i in range(len(rows) - 1):
        a, b = rows[i], rows[i + 1]
        if a['w_star'] != b['w_star'] or a['kind'] != b['kind'] or a['R'] != b['R']:
            continue
        q = np.log(a['h'] / b['h'])
        for k in ('dw', 'damp_rel', 'dm50_rel', 'w0-wstar', 'm50_0/M-1'):
            if b[k] != 0 and a[k] / b[k] > 0:
                print(f"order {k} {a['name']}→{b['name']}: {np.log(a[k] / b[k]) / q:.2f}")

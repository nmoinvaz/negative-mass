"""Analyse run CSVs: omega(t) (from unwrapped central phase, divided by the asymptotic lapse factor,
which is 1 by construction since alpha*a = 1 at R), relative drifts, windows."""
import sys, json, glob, os
import numpy as np

def load(name):
    p = f'runs/{name}.csv'
    hdr = json.loads(open(p).readline()[2:])
    d = np.genfromtxt(p, delimiter=',', names=True, skip_header=1)
    return hdr, d

def omega_series(d, win=50.0):
    t = d['t']; ph = np.unwrap(d['phase'])
    dt = t[1] - t[0]
    k = max(2, int(round(win / dt)))
    tc, wc = [], []
    for i in range(0, len(t) - k, max(1, k // 2)):
        sl = slice(i, i + k + 1)
        A = np.vstack([t[sl] - t[sl].mean(), np.ones(k + 1)]).T
        c = np.linalg.lstsq(A, ph[sl], rcond=None)[0]
        tc.append(t[sl].mean()); wc.append(-c[0])
    return np.array(tc), np.array(wc)

def window_mean(d, col, t0, t1):
    s = (d['t'] >= t0) & (d['t'] <= t1)
    return float(np.mean(d[col][s]))

def summary(name, wins=((0, 100), (450, 550), (900, 1000), (1400, 1500))):
    hdr, d = load(name)
    out = dict(name=name, h=hdr['h'], R=hdr['R'], kind=hdr['kind'], tend=float(d['t'][-1]))
    a0 = d['amp'][0]
    out['amp0'] = float(a0)
    for (t0, t1) in wins:
        if d['t'][-1] < t1: continue
        s = (d['t'] >= t0) & (d['t'] <= t1)
        ph = np.unwrap(d['phase'])[s]; tt = d['t'][s]
        w = -np.polyfit(tt, ph, 1)[0]
        out[f'w[{t0},{t1}]'] = float(w)
        out[f'amp/amp0-1[{t0},{t1}]'] = float(np.mean(d['amp'][s]) / a0 - 1)
        for c in ('m12', 'm20', 'm50', 'N20', 'N50', 'content20'):
            out[f'{c}[{t0},{t1}]'] = float(np.mean(d[c][s]))
    return out

if __name__ == '__main__':
    for n in sys.argv[1:]:
        print(json.dumps(summary(n), indent=1))

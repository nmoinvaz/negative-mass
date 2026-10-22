"""Growth rate of the radial instability: odd part of the +/-delta perturbation pair,
D(t) = (amp_+ - amp_-)/2 (cancels the discretisation seed and O(delta^2)); fit ln|D| on windows."""
import sys, numpy as np
from analyze import load
def fit(pp, pm, wins):
    _, a = load(pp); _, b = load(pm)
    n = min(len(a), len(b)); t = a['t'][:n]
    D = 0.5 * (a['amp'][:n] - b['amp'][:n]) / a['amp'][0]
    S = 0.5 * (a['amp'][:n] + b['amp'][:n]) / a['amp'][0] - 1
    res = []
    for t0, t1 in wins:
        s = (t >= t0) & (t <= t1)
        sig = np.polyfit(t[s], np.log(np.abs(D[s])), 1)[0]
        res.append((t0, t1, sig))
    return t, D, S, res
if __name__ == '__main__':
    pp, pm = sys.argv[1], sys.argv[2]
    t, D, S, res = fit(pp, pm, [(20, 60), (30, 70), (40, 80), (20, 80)])
    for t0, t1, s in res:
        print(f"{pp}: window [{t0},{t1}] sigma = {s:.5f}")
    for i in range(0, 121, 10):
        print(f"  t={t[i]:.0f} D={D[i]:+.4e} S={S[i]:+.4e}")

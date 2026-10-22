"""Fit the odd part D(t) of a +/-delta pair to A e^{s t} + B e^{-s t} + C (the unstable radial mode
comes as a +/-s pair; C absorbs a static shift onto a neighbouring star), on several windows."""
import sys, numpy as np
from scipy.optimize import curve_fit
from growth import fit
f = lambda t, A, B, C, s: A * np.exp(s * t) + B * np.exp(-s * t) + C
for pp, pm, ref in [(sys.argv[i], sys.argv[i + 1], float(sys.argv[i + 2])) for i in range(1, len(sys.argv), 3)]:
    t, D, S, res = fit(pp, pm, [(40, 80)])
    for t0, t1 in [(0, 60), (0, 80), (10, 90), (20, 100)]:
        s = (t >= t0) & (t <= t1)
        p, _ = curve_fit(f, t[s], D[s], p0=[1e-6, 1e-6, 0, 0.07], maxfev=20000)
        print(f"{pp}/{pm} [{t0},{t1}]: sigma = {p[3]:.5f}  (radial {ref:.5f}, rel {p[3]/ref-1:+.2%})  A={p[0]:.2e} B={p[1]:.2e} C={p[2]:.2e}")

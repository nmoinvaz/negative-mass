"""
Task 5 support: (i) how far the 1PN dipole trajectory tracks exact hyperbolic motion, and where
the truncation visibly fails; (ii) mass multipoles of the Newtonian uniformly accelerated dipole
(radiation estimate). Units G = c = 1.
"""
import numpy as np
import sympy as sp
from integrate_1pn import run

m, d = 0.0143, 10.0
# rigid (boost-symmetric) 1PN configuration: |mu2| = m (1 + m/d); proper acceleration of the midpoint
# at t = 0 is alpha = m/d^2 (eih.py). Exact hyperbolic motion: v = alpha t / sqrt(1 + alpha^2 t^2).
# exact 1PN midpoint acceleration of the rigid pair at rest is (m/d^2)(1 + m/(2d)) (collinear.py)
alpha = m / d**2 * (1 + m / (2 * d))
tmax = 8000.0
o = run(m, -m * (1 + m / d), d, tmax, "exact", npts=801)
print(f"rigid dipole, m = {m}, d = {d}, alpha = (m/d^2)(1+m/2d) = {alpha:.5e}, 1/alpha = {1/alpha:.0f}")
print("   t     v_mid(1PN)  v_hyperbola  rel.err   v_Newton   d*gamma(1PN)")
for tt in (1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000):
    i = np.argmin(abs(o['t'] - tt))
    if i >= len(o['t']):
        break
    vm = o['vmid'][i]
    vh = alpha * tt / np.sqrt(1 + (alpha * tt)**2)
    g = 1 / np.sqrt(1 - vm * vm) if abs(vm) < 1 else np.nan
    print(f"{tt:6.0f}   {vm:.5f}     {vh:.5f}     {vm/vh-1:+.2e}  {alpha*tt:.4f}    {o['sep'][i]*g:.5f}")
print("integration reached t =", o['t'][-1])

# Newtonian multipoles of the uniformly accelerated dipole: x1 = X + d/2, x2 = X - d/2, X = a t^2 / 2
t, a, dd, mm = sp.symbols('t a d m', positive=True)
X = a * t**2 / 2
x1, x2 = X + dd / 2, X - dd / 2
for l in (1, 2, 3, 4):
    Il = mm * x1**l - mm * x2**l
    derivs = [sp.simplify(sp.diff(Il, t, l + 1)), sp.simplify(sp.diff(Il, t, l + 2))]
    print(f"mass {2**l}-pole  I_l = sum mu x^l:  d^(l+1)/dt^(l+1) = {derivs[0]},  d^(l+2) = {derivs[1]}")

"""
Newtonian finite-size check for Shirokov's dipole at d = 10: the two stars' profiles (Phase 2 code, read-only
import) are not compact at d = 10 (tails reach r ~ 20), so the shell theorem does not apply. Compute the
Newtonian force between two rigid spherical distributions with the actual mass profiles m(r):

    W(d) = -(2/pi) int_0^inf rho1^(k) rho2^(k) j0(kd) dk,   rho^(k) = int m'(r) j0(kr) dr,
    F(d) = -dW/dd.

For point masses this gives W = -m1 m2/d. We report F d^2 / (m1 m2) (= 1 for points). Also reports what
fraction of each star's mass lies beyond d/2 from its centre.
"""
import sys, os
import numpy as np
from scipy.integrate import trapezoid
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'phase2'))
import boson_star as bs

def profile(w_t, eps):
    s = bs.star_at_physical_frequency(w_t, eps)
    _, sol = bs.shoot(s['phic'], w_t * 0 + s['w_phys'] * 0 + 1, eps) if False else (None, None)
    return s

def mass_profile(w_t, eps):
    # recompute the cut solution exactly as in star_properties
    w = w_t
    for _ in range(4):
        p = bs.find_star(w, eps)
        s = bs.star_properties(p, w, eps)
        w = w * w_t / s['w_phys']
    _, sol = bs.shoot(s['phic'], w, eps)
    r, m = sol.t, sol.y[2]
    i = np.searchsorted(r, s['rcut'])
    r, m = r[:i + 1], m[:i + 1]
    return r, m, s['M']

r1, m1r, M1 = mass_profile(0.750, +1)
r2, m2r, M2 = mass_profile(0.7603, -1)
print(f"canonical M = {M1:+.6f}, phantom M = {M2:+.6f}, |M2|/M1 = {abs(M2)/M1:.5f}")

def ft(r, m, k):
    dm = np.gradient(m, r)
    return np.array([trapezoid(dm * np.sinc(kk * r / np.pi), r) for kk in k])

k = np.linspace(1e-5, 8.0, 40001)
R1, R2 = ft(r1, m1r, k), ft(r2, m2r, k)
for d in (8.0, 10.0, 12.0, 15.0, 20.0, 25.0):
    # dj0(kd)/dd = k j0'(kd), j0'(x) = (x cos x - sin x)/x^2
    x = k * d
    dj0 = k * (x * np.cos(x) - np.sin(x)) / x**2
    F = -(-(2 / np.pi) * trapezoid(R1 * R2 * dj0, k))      # F = -dW/dd ;  >0 means repulsive (W decreasing)
    W = -(2 / np.pi) * trapezoid(R1 * R2 * np.sinc(x / np.pi), k)
    f1 = 1 - np.interp(d / 2, r1, m1r) / M1
    f2 = 1 - np.interp(d / 2, r2, m2r) / M2
    print(f"d = {d:5.1f}:  F d^2/(M1 M2) = {-F * d * d / (M1 * M2):.5f}   W d/(M1 M2) = {-W * d / (M1 * M2):.5f}"
          f"   mass beyond d/2: {f1:.2e} (canonical), {f2:.2e} (phantom)")
print("(1.00000 = point-particle Newtonian force; < 1 means overlap weakens the interaction)")

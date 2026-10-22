"""Toy adiabatic radiation-reaction model in the pair's instantaneous rest frame (heuristic, weak field).
Phase 3 (1PN): proper separation obeys d'' = -E_tot/d^2. Radiation: dE_tot/dtau = -L0(d) = -k/d^6, k = 24 m^6/35.
Rapidity eta = int A dtau, A = m/d^2. Lab-frame radiated energy dE_lab = L0 cosh(eta) dtau (dP^mu = L0 u^mu dtau).
Scaled: x = d/d0, s = tau/T, T = d0^3/k^(1/3), e = E_tot T^2/d0^3:  x'' = -e/x^2, e' = -1/x^6.
"""
import numpy as np
from scipy.integrate import solve_ivp
def rhs(s, y):
    x, v, e, I = y
    return [v, -e/x**2, -1/x**6, 1/x**2]
for e0 in [0.0, -1e-3, 1e-3]:
    sol = solve_ivp(rhs, [0, 1e6], [1, 0, e0, 0], rtol=1e-11, atol=1e-13, dense_output=True,
                    events=lambda s, y: y[0]-0.05)
    for s in [1, 3, 10, 100, 1e4, 1e6]:
        if s <= sol.t[-1]:
            x, v, e, I = sol.sol(s)
            print(f"e0={e0:+.0e} s={s:8.0e}: x={x:.4g} x'={v:.4g} e={e:.6g} int ds/x^2={I:.6f}")
    print("  terminated at s =", sol.t[-1], "(collision)" if sol.t_events[0].size else "")
m, d0 = 0.0143, 10.0
k = 24*m**6/35; T = d0**3/k**(1/3)
print(f"Shirokov parameters: T = {T:.4g}, eta scale m T/d0^2 = {m*T/d0**2:.4g}, E_tot scale d0^3/T^2 = {d0**3/T**2:.3g}")

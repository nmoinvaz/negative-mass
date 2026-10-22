"""
Phase 4, Lemma B check: a priori bound for time-symmetric data from Witten + refined Kato + Sobolev.

  M >= - E_- / (1 - kappa),   kappa = (8 pi / 3) ||mu_-||_{L^{3/2}(g)} / S_g,   valid if kappa < 1 (kappa0 = 4 pi ||mu_-||/S_g printed for reference),
  E_- = int mu_- dV_g, S_g = Sobolev constant of g (int|dh|^2 >= S_g ||h||_6^2), S_g >= S_0 (inf psi / sup psi)^2
  for g = psi^4 delta, S_0 = 3 (pi/2)^{4/3}. (Lower bound on S_g -> upper bound on kappa -> weaker but rigorous.)
"""
import sys, os
import numpy as np
from scipy.integrate import cumulative_trapezoid as ctz
sys.path.insert(0, os.path.dirname(__file__))
import importlib.util
S0 = 3 * (np.pi / 2) ** (4 / 3)
trap = np.trapezoid

# phantom boson star (reuse loader from witten_bound without re-running its prints)
spec = importlib.util.spec_from_file_location("wb_funcs", os.path.join(os.path.dirname(__file__), "witten_bound.py"))
src = open(spec.origin).read().split("print(\"(1)")[0]
ns = {'__file__': spec.origin}; exec(src, ns)

def report(name, M, E, L32, psimin, psimax):
    Sg = S0 * (psimin / psimax) ** 2
    k0 = 4 * np.pi * L32 / Sg
    k = 2 * k0 / 3
    b = -E / (1 - k) if k < 1 else float('nan')
    print(f"  {name:22s} M={M:+.6f} E_-={E:.6f} ||mu_-||_3/2={L32:.3e} S_g>={Sg:.3f} kappa0={k0:.3e} "
          f"bound={b:+.6f} kappa={k:.3e} holds={M >= b if k < 1 else 'n/a'} ratio M/bound={M/b if k<1 else float('nan'):.5f}")

print("Lemma B: M >= -E_-/(1-kappa)")
for w_t in [0.7603, 0.9050]:
    rA, a, rho, M = ns['star_profile'](w_t, -1)
    rc = rA[-1]
    riso_c = 0.5 * (rc - M + np.sqrt(rc * rc - 2 * M * rc))
    I = -ctz(((a - 1.0) / rA)[::-1], rA[::-1], initial=0.0)[::-1]
    riso = np.exp(np.log(riso_c) - (np.log(rc) - np.log(rA)) - I)
    psi = np.sqrt(rA / riso)
    dV = 4 * np.pi * rA**2 * a
    E = trap(rho * dV, rA)
    L32 = trap(rho**1.5 * dV, rA) ** (2 / 3)
    # psi outside: isotropic Schwarzschild with M<0 -> psi = 1 + M/(2r) increases to 1; sup psi = 1
    report(f"phantom star w={w_t}", M, E, L32, psi.min(), max(psi.max(), 1.0))

r = np.linspace(1e-6, 60, 600001)
for Mabs, wdt in [(0.1, 0.05), (0.1, 0.3), (0.5, 0.3), (1.0, 0.5)]:
    mt = np.exp(-0.5 * ((r - 1.0) / wdt) ** 2); mt *= Mabs / trap(4*np.pi*r**2*mt, r)
    Q = ctz(4*np.pi*r**2*mt, r, initial=0.0)
    dpsi = Q / (2 * r**2)                        # phantom: Delta psi = +2 pi mt psi^... sign: psi increases outward
    psi = 1.0 - Mabs/(2*r[-1]) + (ctz(dpsi, r, initial=0.0) - ctz(dpsi, r, initial=0.0)[-1])
    mu = mt * psi**-5
    dV = 4*np.pi*r**2*psi**6
    E = trap(mu*dV, r); L32 = trap(mu**1.5*dV, r)**(2/3)
    report(f"shell |M|={Mabs} w={wdt}", -Mabs, E, L32, psi.min(), 1.0)

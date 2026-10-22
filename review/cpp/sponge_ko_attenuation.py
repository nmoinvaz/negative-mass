# Independent estimate: the "sponge" is extra 6th-difference Kreiss-Oliger
# dissipation (Source/Grids/SpongeZone.hpp). For a 1D Fourier mode e^{ikx}
# the KO term used (weights 1/64*[1,-6,15,-20,15,-6,1] / dx) gives a damping
# rate  r(k) = sigma/dx * sin^6(k dx / 2).  Estimate the amplitude attenuation
# of a mode crossing the sponge shell once (thickness W, speed ~1).
import numpy as np
w = np.array([1.5625e-2,-9.375e-2,2.34375e-1,-3.125e-1,2.34375e-1,-9.375e-2,1.5625e-2])
assert abs(w.sum()) < 1e-15
for dx in (0.5, 0.25):
    for lam in (16.0, 8.0, 4.0, 2.0, 1.0):
        k = 2*np.pi/lam
        sym = sum(w[j]*np.cos(k*dx*(j-3)) for j in range(7))  # transfer fn
        rate_unit = -sym/dx          # per unit sigma
        # sponge ramp ((r-24)/8)^4 * 4 + base sigma 2, crossing W=8 at v=1
        r = np.linspace(24, 32, 2001)
        sig_extra = 4.0*((r-24)/8)**4
        att_extra = np.exp(-np.trapezoid(sig_extra*rate_unit, r))
        print(f"dx={dx:5.3f} lambda={lam:5.1f} kdx={k*dx:5.3f} "
              f"rate/sigma={rate_unit:9.3e} sponge-only one-pass amplitude factor={att_extra:.6f}")

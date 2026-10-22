"""
Axial (odd-parity) sector.  derive_axial.py shows (sympy, eps symbolic, all three equations used):
    Q = e^{(nu-lam)/2} h1 / r,     d^2Q/dr*^2 = (s^2 + V_ax) Q,    dr*/dr = e^{(lam-nu)/2},
    V_ax = e^nu [ l(l+1)/r^2 - 6 m(r)/r^3 + 4 pi eps (rho - p_r) ],   rho - p_r = 2 V(phi0^2)
(eps = +1: the anisotropic-fluid / boson-star Regge-Wheeler potential of Kojima et al. and Macedo et al.
eq. for V_RW with kappa/2 (rho - p_rad) = 4 pi (rho - p_r); vacuum: Regge-Wheeler).

A growing mode e^{st}, s > 0, needs a normalisable negative eigenvalue -s^2 of  -d^2/dr*^2 + V_ax  on
r* in (0, oo) with the regular condition Q ~ r^{l+1}.  If V_ax > 0 everywhere this is impossible
(multiply by Q and integrate).  Here we evaluate V_ax and its pieces on the backgrounds, and give
r^2 V_ax e^{-nu} / [l(l+1)] (the ratio to the centrifugal term), whose minimum certifies positivity.
"""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B

if __name__ == "__main__":
    for eps, w in [(-1, 0.7603), (+1, 0.75), (-1, 0.70), (-1, 0.85), (-1, 0.90)]:
        st = B.solve_star(w, eps)
        r = np.linspace(st.r0, st.R, 200001)
        al, a, phi, u = st.bg(r)
        m = r * (1 - 1 / a**2) / 2
        Vpot = B.V(phi**2)
        for l in (2, 3):
            LL = l * (l + 1)
            Vax = al**2 * (LL / r**2 - 6 * m / r**3 + 8 * np.pi * eps * Vpot)
            ratio = Vax * r**2 / (al**2 * LL)
            i = np.argmin(ratio)
            print(f"eps={eps:+d} w={w} M={st.M:+.6f} l={l}: min V_ax={Vax.min():.4e} (r={r[np.argmin(Vax)]:.2f}); "
                  f"min r^2 V_ax/(alpha^2 l(l+1)) = {ratio[i]:.6f} at r={r[i]:.3f};  "
                  f"max |6m/r^3|/(l(l+1)/r^2) = {np.max(np.abs(6*m/r)/LL):.2e}, max |8 pi V| r^2/(l(l+1)) = {np.max(np.abs(8*np.pi*Vpot)*r**2/LL):.2e}")

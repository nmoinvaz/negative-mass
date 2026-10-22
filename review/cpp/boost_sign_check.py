# Independent check of the sign convention between the Python painter
# (grtresna/fields/boson_star.py: Pi = -dt_phi / alpha) and the evolution
# (BiComplexScalarField.impl.hpp:132: d_t phi = alpha*Pi + beta.d phi).
# Flat space, alpha=1, beta=0, free massive field m=1, 1D plane-wave packet
# painted with "velocity" +v using the painter's convention; evolve with the
# code's convention (phi_t = Pi, Pi_t = phi_xx - m^2 phi) and report the
# direction the packet centroid moves.
import numpy as np
L, N = 200.0, 4000; x = np.linspace(-L/2, L/2, N, endpoint=False); dx = x[1]-x[0]
m, om, v = 1.0, 1.0, 0.3
g = 1/np.sqrt(1-v*v); k = om*g*v; w = np.sqrt(k*k+m*m)
env = np.exp(-x**2/(2*8.0**2))
# packet moving +x: phi = env(x - v t) cos(k x - w t)  -> dt phi at t=0
phi = env*np.cos(k*x)
dtphi = (v*x/8.0**2)*env*np.cos(k*x) + w*env*np.sin(k*x)
for label, Pi in (("painter Pi=-dt_phi", -dtphi), ("code Pi=+dt_phi", dtphi)):
    p, P = phi.copy(), Pi.copy(); dt = 0.2*dx
    def lap(f): return (np.roll(f,-1)-2*f+np.roll(f,1))/dx**2
    for _ in range(int(40/dt)):           # leapfrog-ish RK2
        k1p, k1P = P, lap(p)-m*m*p
        p2, P2 = p+0.5*dt*k1p, P+0.5*dt*k1P
        p, P = p+dt*P2, P+dt*(lap(p2)-m*m*p2)
    e = P**2 + np.gradient(p,dx)**2 + m*m*p**2
    print(f"{label:22s} centroid after t=40: {np.sum(e*x)/np.sum(e):+.2f} (expected +{v*40*w/w:.1f}-ish for group velocity {k/w:.2f})")

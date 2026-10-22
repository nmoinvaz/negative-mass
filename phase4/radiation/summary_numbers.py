"""Numbers for RESULTS.md. Linearized: dP_rad^mu/dtau = L0 u^mu, L0 = 24 m^6/(35 d^6), A = m/d^2 (rigid pair, released at rest at t=0).
P_sys(t) = P0 - L0 (t, (gamma-1)/A).  Phase-4 floor: W_- ~ gamma(1+v) m = e^eta m."""
import mpmath as mp
m, d = mp.mpf('0.0143'), mp.mpf(10)
L0 = 24*m**6/(35*d**6); A = m/d**2
print("L0 =", mp.nstr(L0, 4), " A =", mp.nstr(A, 4))
for t in [1e3, 1e6, 1e12, 1e16]:
    t = mp.mpf(t); eta = mp.asinh(A*t); g = mp.cosh(eta)
    E = -L0*t; Pz = -L0*(g-1)/A; inv = mp.sqrt(E**2 - Pz**2)
    W = mp.e**eta*m
    print(f"t={mp.nstr(t,3)}: eta={mp.nstr(eta,4)}  E_B-E0={mp.nstr(E,4)}  |P_sys| (past-directed)={mp.nstr(inv,4)}  "
          f"approx L0*sqrt(2t/A)={mp.nstr(L0*mp.sqrt(2*t/A),4)}  W_-={mp.nstr(W,4)}  E_rad/W_-={mp.nstr(-E/W,4)}")
print("asymptotic ratio E_rad/W_- -> L0/(2 A m) = (12/35)(m/d)^4 =", mp.nstr(L0/(2*A*m), 4))
print("lab time for E_B to fall by m:", mp.nstr(m/L0, 4))
# near-scri check of BS luminosity singularity L ~ 3/(8u^2) (units D=a=1)
cu = lambda u, x: 3*u*(1-x**2)/(u**2+1-x**2)**2.5
for u in [1e-2, 1e-3]:
    u = mp.mpf(u); L = 0.5*mp.quad(lambda x: cu(u, x)**2, [-1, 0, 1-u, 1-u**2, 1])
    print(f"u={mp.nstr(u,2)}: u^2 L(u) = {mp.nstr(u**2*L, 8)} (3/8 = 0.375)")

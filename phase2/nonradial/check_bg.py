"""Quick look at the backgrounds: tail of phi0, and the 'Zerilli-pole' denominator
D = 3LL e^{lam+nu} + ... (common denominator of the reduced polar system) along the star."""
import numpy as np, background as B
for eps, w in [(+1, 0.75), (-1, 0.7603)]:
    st = B.solve_star(w, eps)
    rr = np.array([0.01, 1, 2, 4, 6, 8, 10, 15, 20, 25, 30, 40, 50, 59])
    al, a, phi, u = st.bg(rr)
    print(f"eps={eps:+d} w={w} M={st.M:+.6f} R={st.R} phic={st.phic:.5f} alpha_c={st.alpha_c:.6f}")
    for x, p in zip(rr, phi): print(f"   r={x:5.2f} phi0={p:.3e}")
    r = np.linspace(st.r0, st.R, 20000); al, a, phi, u = st.bg(r)
    LAM, MU = B.LAM, B.MU; pi = np.pi
    for l in (2, 3):
        LL = l*(l+1); el, en = a**2, al**2; dphi = u*phi
        Dn = (3*LL*el*en + 6*pi*LAM*eps*r**2*phi**4*el*en - 4*pi*eps*MU*r**2*phi**6*el*en
              - 12*pi*eps*w**2*r**2*phi**2*el - 12*pi*eps*r**2*phi**2*el*en + 36*pi*eps*r**2*en*dphi**2
              + 3*el*en - 9*en)
        print(f"   l={l}: min D/(3(LL-2)) = {np.min(Dn/(3*(LL-2))):.6f}  max = {np.max(Dn/(3*(LL-2))):.6f}")

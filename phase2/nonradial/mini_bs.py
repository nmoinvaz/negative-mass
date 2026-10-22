"""
Validation against a published boson-star polar QNM:
Macedo, Pani, Cardoso & Crispino, PRD 88, 064046 (2013) [arXiv:1307.4812], "Mini BS I":
  omega = 0.853087 m,  M = 0.63300 / m  (the maximum-mass mini star),
  fundamental l = 2 polar mode (direct integration)  sigma = 0.1195 - 5e-5 i  (units of m; e^{-i sigma t}).
Mini star: V = m^2 |Phi|^2 (their normalisation) = our V(s) = s/2 with Lambda = mu = 0.  The model is
identified by (omega, M), which do not depend on the field normalisation.

Steps: shoot a mini star (alpha(0) = 1, code frequency), polish with background.solve_star (alpha(oo) = 1,
continuation in w to the target), then solve D(s) = 0 for s = -i sigma near the published value
(complex secant on the matching determinant of polar_solver, which is holomorphic up to a positive factor;
exterior Zerilli solution = outgoing ~ e^{i sigma r*}).
"""
import sys, os, numpy as np, time
from scipy.integrate import solve_ivp, solve_bvp
from scipy.interpolate import interp1d
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import boson_star as bs
bs.LAM = 0.0; bs.MU = 0.0            # mini boson star (background.py reads these through bs.dVds2 / bs.V)
import background as B
from polar_solver import PolarProblem

def shoot(phic, w, rmax=30.0):
    c2 = (bs.dVds2(phic**2) - w * w) * phic / 6.0
    r0 = 1e-6
    y0 = [phic + c2 * r0**2, 2 * c2 * r0, 0.0, 0.0]
    def cross(r, y, *a): return y[0]
    cross.terminal = True
    def up(r, y, *a): return y[1] if y[0] < 0.5 * phic else -1
    up.terminal = True; up.direction = 1
    sol = solve_ivp(bs.rhs, (r0, rmax), y0, args=(w, +1), events=(cross, up), rtol=1e-11, atol=1e-15, method='DOP853')
    return (-1 if sol.t_events[0].size else +1), sol

def find_w(phic, lo=0.8, hi=2.5):
    grid = np.linspace(lo, hi, 18)
    sg = [shoot(phic, x)[0] for x in grid]
    k = next(i for i in range(len(grid) - 1) if sg[i] == +1 and sg[i + 1] == -1)
    lo, hi = grid[k], grid[k + 1]
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        if shoot(phic, mid)[0] == +1: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)

def star_from_shoot(phic, R=50.0, r0=1e-3):
    wc = find_w(phic)
    _, sol = shoot(phic, wc)
    r, y = sol.t, sol.y
    i = int(np.argmin(np.abs(y[0]) + (r < 2) * 1e9))
    rc = r[i] * 0.8
    al_inf = np.exp(np.interp(rc, r, y[3])) / np.sqrt(1 - 2 * np.interp(rc, r, y[2]) / rc)
    w = wc / al_inf
    k = np.sqrt(1 - w * w)
    x = np.concatenate([np.linspace(r0, rc, 400), np.linspace(rc + 1e-3, R, 200)])
    f = interp1d(r, y, axis=1)
    yi = f(np.minimum(x, rc))
    psi = np.log(np.abs(yi[0])); u = yi[1] / yi[0]; m = yi[2].copy(); nu = 2 * (yi[3] - np.log(al_inf))
    out = x > rc
    psi[out] = psi[~out][-1] - k * (x[out] - rc); u[out] = -k - 1 / x[out]
    m[out] = m[~out][-1]; nu[out] = np.log(1 - 2 * m[out] / x[out])
    solb = solve_bvp(lambda rr, yy: B.odes(rr, yy, w, +1), lambda ya, yb: B.bcs(ya, yb, r0, R, w, +1),
                     x, np.vstack([psi, u, m, nu]), tol=1e-9, max_nodes=200000, bc_tol=1e-10)
    assert solb.status == 0, solb.message
    return B.Star(solb, w, +1, r0, R)

if __name__ == "__main__":
    T0 = time.time()
    # bracket the target frequency on the stable (low phi_c) side of the maximum-mass star
    out = []
    for pc in [0.08, 0.10, 0.12, 0.14]:
        st = star_from_shoot(pc)
        out.append((pc, st.w, st.M)); print(f"phi_c={pc}: w={st.w:.6f} M={st.M:.6f}", flush=True)
    wt = 0.853087
    # continue in w with background.solve_star
    st = star_from_shoot(min(out, key=lambda t: abs(t[1] - wt))[0])
    for _ in range(3):
        st = B.solve_star(wt, +1, guess=st, R=st.R)
    print(f"target star: w={st.w:.6f}  M={st.M:.6f}  (paper: 0.853087, 0.63300)  phi_c={st.phic:.5f}  "
          f"[paper phi_c~ = 0.1916 in their normalisation; ours should be sqrt(2)*0.1916/sqrt(4 pi)?]  {time.time()-T0:.0f}s", flush=True)
    P = PolarProblem(st, 2, LAM=0.0, MU=0.0, phi_cut=1e-9)
    print(f"r_m={P.rm:.2f}  M_ext={P.Mext:.6f}")
    # scan |D| along sigma real near the published value, then complex secant
    def Ds(sig): return P.D(-1j * complex(sig))
    for sg in [0.10, 0.11, 0.115, 0.1195, 0.125, 0.13, 0.14]:
        d = Ds(sg); print(f"   sigma={sg}: D={d:.4e} |D|={abs(d):.3e}", flush=True)
    x0, x1 = 0.118 - 1e-4j, 0.121 - 1e-4j
    f0, f1 = Ds(x0), Ds(x1)
    for it in range(30):
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0, x1, f1 = x1, f1, x2, Ds(x2)
        print(f"   secant {it}: sigma = {x1.real:.7f} {x1.imag:+.3e} i   |D|={abs(f1):.2e}", flush=True)
        if abs(x1 - x0) < 1e-9: break
    print(f"RESULT Mini BS I l=2 fundamental polar: sigma = {x1.real:.6f} {x1.imag:+.2e} i   (paper 0.1195 - 5e-5 i)  {time.time()-T0:.0f}s")

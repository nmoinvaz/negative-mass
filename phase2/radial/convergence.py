"""Convergence of the shooting eigenvalue w.r.t. integrator tolerance, matching radius, outer radius,
inner radius and background BVP tolerance, at two stars near the turning point and the campaign stars."""
import numpy as np, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import brentq
from background import solve_star
from radial_modes import RadialProblem

def root_near(P, s2a, s2b):
    return brentq(P.D, s2a, s2b, xtol=1e-14, rtol=1e-12)

def cont(eps, w, tol=1e-9, R=None):
    st = solve_star(0.75, eps, tol=tol, R=R); cur = 0.75
    while abs(w - cur) > 1e-12:
        cur = cur + np.clip(w - cur, -0.02, 0.02); st = solve_star(cur, eps, guess=st, tol=tol, R=R)
    return st

for eps, w, br in [(+1, 0.92, (1e-3, 4e-3)), (+1, 0.94, (-8e-3, -3e-3)),
                   (-1, 0.9275, (2e-4, 1e-3)), (-1, 0.94, (-6e-3, -3e-3))]:
    st = cont(eps, w)
    base = RadialProblem(st)
    s_ref = root_near(base, *br)
    print(f"eps={eps:+d} w={w}: reference sigma^2 = {s_ref:.10e}  (rm={base.rm:.2f}, R={base.R:.1f}, r0={base.r0})")
    for label, P in [("rtol 1e-12", RadialProblem(st, rtol=1e-12)),
                     ("rm x0.7", RadialProblem(st, rm=0.7 * base.rm)),
                     ("rm x1.4", RadialProblem(st, rm=1.4 * base.rm)),
                     ("R_eig = R-15", RadialProblem(st, R=base.R - 15))]:
        print(f"    {label:14s}: {root_near(P, *br):.10e}")
    st2 = cont(eps, w, tol=1e-10)
    print(f"    {'bg tol 1e-10':14s}: {root_near(RadialProblem(st2), *br):.10e}")
    st3 = cont(eps, w, R=base.R + 30)
    print(f"    {'bg R +30':14s}: {root_near(RadialProblem(st3), *br):.10e}   (M={st3.M:+.9f} vs {st.M:+.9f})", flush=True)

# campaign stars: D has no zero on [-0.05, (1-w)^2) -- show min |D| and sign constancy at two tolerances
for eps, w in [(+1, 0.75), (-1, 0.7603)]:
    st = cont(eps, w)
    for rt in (1e-10, 1e-12):
        P = RadialProblem(st, rtol=rt)
        g = np.linspace(-0.05, 0.999 * (1 - w)**2, 60)
        v = np.array([P.D(s) for s in g])
        print(f"eps={eps:+d} w={w} rtol={rt}: sign changes of D on [-0.05,(1-w)^2): {int(np.sum(np.diff(np.sign(v)) != 0))}; "
              f"min|D| = {np.min(np.abs(v)):.2e}, max|D| = {np.max(np.abs(v)):.2e}", flush=True)

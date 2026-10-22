"""Solve the static stars with phase2/radial/background.py (BVP, alpha(inf)=1, true ADM masses) and save
dense profiles (r, phi0, alpha, a, m) to stars/*.npz for the evolution code.  Read-only import."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'phase2', 'radial'))
import background as bg

def save(w, eps, tag):
    out = f"stars/{tag}.npz"
    if os.path.exists(out):
        return
    st = bg.solve_star(w, eps)
    r = np.linspace(st.r0, st.R, 200001)
    y = st.sol.sol(r)
    np.savez(out, r=r, phi=np.exp(y[0]), u=y[1], m=y[2], nu=y[3], w=w, eps=eps, M=st.M, N=st.N,
             phic=st.phic, R=st.R)
    print(f"{tag}: w={w} eps={eps:+d} M={st.M:+.8f} N={st.N:.8f} phic={st.phic:.8f} alpha_c={st.alpha_c:.8f}", flush=True)

if __name__ == "__main__":
    os.makedirs("stars", exist_ok=True)
    jobs = [(0.75, 1, 'c0750'), (0.7603, -1, 'p07603'), (0.94, 1, 'c0940'), (0.94, -1, 'p0940')]
    # neighbours for dM/domega, dN/domega at the campaign stars
    jobs += [(0.74, 1, 'c0740'), (0.76, 1, 'c0760'), (0.7503, -1, 'p07503'), (0.7703, -1, 'p07703')]
    for w, e, t in jobs:
        save(w, e, t)

"""Continue a family in N (w free) to locate folds of w(N) and the behaviour of M there."""
import numpy as np, warnings, sys, pickle
warnings.filterwarnings('ignore')
from background import solve_star, solve_star_N
eps = int(sys.argv[1]); w_start = float(sys.argv[2]); Nmax = float(sys.argv[3])
st = solve_star(0.66, eps, tol=1e-8); w = 0.66; dw = 0.01
while w > w_start + 1e-9:
    w2 = round(max(w - dw, w_start), 6)
    try: st = solve_star(w2, eps, guess=st, tol=1e-8); w = w2; dw = min(dw * 1.5, 0.01)
    except RuntimeError: dw /= 2
print(f"start eps={eps:+d} w={st.w:.5f} N={st.N:.5f} M={st.M:+.5f}", flush=True)
lnN = np.log(st.N); dl = 0.05; out = []
while np.exp(lnN) < Nmax:
    try:
        st2 = solve_star_N(np.exp(lnN + dl), eps, st)
    except RuntimeError as e:
        dl /= 2
        if dl < 1e-4: print("stuck", e); break
        continue
    st = st2; lnN += dl; dl = min(dl * 1.3, 0.05)
    r = np.linspace(st.r0, st.R, 6000); al, a, phi, u = st.bg(r)
    out.append((st.w, st.N, st.M, st.phic, al[0], st.sol))
    print(f"eps={eps:+d} w={st.w:.6f} N={st.N:.5f} M={st.M:+.6f} phic={st.phic:.6f} alpha_c={al[0]:.4f} "
          f"max|2m/r|={np.max(np.abs(1-1/a**2)):.3f} mesh={st.nmesh} R={st.R:.0f}", flush=True)


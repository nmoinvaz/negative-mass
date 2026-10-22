"""Strong-field end of the families: continue in the central lapse alpha_c through the N/M maximum
(canonical) and the w-fold (phantom); lowest radial sigma^2 at each step."""
import numpy as np, warnings, sys, json
warnings.filterwarnings('ignore')
from background import solve_star, solve_star_N, solve_star_lapse
from radial_modes import RadialProblem, lowest_modes
eps = int(sys.argv[1]); targets = [float(x) for x in sys.argv[2].split(',')]
st = solve_star(0.66, eps, tol=1e-8); w = 0.66; dw = 0.01
while w > 0.45 + 1e-9:
    w2 = round(max(w - dw, 0.45), 6)
    try: st = solve_star(w2, eps, guess=st, tol=1e-8); w = w2; dw = min(dw * 1.5, 0.01)
    except RuntimeError: dw /= 2
ac = st.alpha_c; rows = []
for t in targets:
    d = 0.01 * np.sign(t - ac)
    while abs(t - ac) > 1e-12:
        step = d if abs(d) < abs(t - ac) else t - ac
        try: st = solve_star_lapse(ac + step, eps, st); ac += step; d *= 1.3
        except RuntimeError:
            d /= 2
            if abs(d) < 1e-6: raise
        d = np.sign(d) * min(abs(d), 0.01)
    P = RadialProblem(st)
    roots, g, vals = lowest_modes(P, s2_min=-0.2, nneg=24, npos=40, nmax=2, s2_max=0.1)
    rows.append(dict(alpha_c=ac, w=st.w, N=st.N, M=st.M, roots=[float(x) for x, o in roots]))
    print(f"eps={eps:+d} alpha_c={ac:.4f} w={st.w:.6f} N={st.N:.5f} M={st.M:+.6f} sigma2={['%.4e%s' % (x, '' if o else '?') for x, o in roots]} rm={P.rm:.1f}", flush=True)
json.dump(rows, open(f"strongfield_eps{eps:+d}.json", "w"), indent=1)

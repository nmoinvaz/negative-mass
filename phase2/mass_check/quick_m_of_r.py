"""Our Schwarzschild-gauge star: Misner-Sharp mass m(R) at the radius where phi first
drops below 5% of phi_c (the paper's match/extraction radius) vs the full mass."""
import sys, numpy as np
sys.path.insert(0, '/home/claude/negative-mass/phase2')
import boson_star as bs

for w_t, eps, Mp in [(0.750, +1, 0.014350), (0.7603, -1, -0.014295)]:
    w = w_t
    for _ in range(4):
        p = bs.find_star(w, eps); s = bs.star_properties(p, w, eps)
        w_used = w; w = w * w_t / s['w_phys']
    _, sol = bs.shoot(p, w_used, eps)
    r, phi, m = sol.t, sol.y[0], sol.y[2]
    i5 = np.where((phi < 0.05 * p) & (r > 1e-3))[0][0]
    # interpolate m at exact crossing
    R5 = np.interp(0.05 * p, phi[i5:i5-2:-1], r[i5:i5-2:-1])
    m5 = np.interp(R5, r, m)
    print(f"eps={eps:+d} w_phys={s['w_phys']:.6f} phic={p:.8f} M_full={s['M']:+.6f}  "
          f"R(phi=5%)={R5:.4f} m(R5)={m5:+.6f}  paper={Mp:+.6f}  m5/paper={m5/Mp:.5f}  M/paper={s['M']/Mp:.5f}")

"""Compare relaxed end states of the 5%-bump runs with the static BVP star at the measured late-time
frequency (is the remnant a member of the static family?)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'phase2', 'radial'))
import background as bg
from analyze import load
for name, eps in [('c750_pert5_h0.1', 1), ('p7603_pert5_h0.1', -1)]:
    h, d = load(name)
    s = d['t'] > 1000
    w = float(-np.polyfit(d['t'][s], np.unwrap(d['phase'])[s], 1)[0])
    st = bg.solve_star(w, eps)
    print(f"{name}: late w={w:.6f}; evolved amp={np.mean(d['amp'][s]):.7e} m20={np.mean(d['m20'][s]):+.8f} N50={np.mean(d['N50'][s]):.8f}")
    print(f"   static BVP at that w: phic={st.phic:.7e} M={st.M:+.8f} N={st.N:.8f};  rel diffs amp {np.mean(d['amp'][s])/st.phic-1:+.1e} "
          f"m20 {np.mean(d['m20'][s])/st.M-1:+.1e} N {np.mean(d['N50'][s])/st.N-1:+.1e}")
    print(f"   initial (t=0) m(50)={d['m50'][0]:+.8f} N(50)={d['N50'][0]:.8f}: radiated fraction of N {1-np.mean(d['N50'][s])/d['N50'][0]:.2e}, of M {1-np.mean(d['m50'][s])/d['m50'][0]:.2e}")

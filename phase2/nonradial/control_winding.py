"""
Positive control for the argument-principle counter (winding.py logic): a small rectangle around the
Mini BS I l = 2 polar QNM found by mini_bs.py (s = -i sigma, sigma = 0.1195191 - 5.353e-5 i, i.e.
s = -5.353e-5 - 0.1195191 i) must give winding +1, and a displaced rectangle 0.
"""
import sys, os, numpy as np, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mini_bs as mb                     # sets LAM = MU = 0 (mini star)
import background as B
from polar_solver import PolarProblem

def winding(P, c, hx, hy):
    corners = [c - hx - 1j*hy, c + hx - 1j*hy, c + hx + 1j*hy, c - hx + 1j*hy, c - hx - 1j*hy]
    cache = {}
    def Dc(s):
        if s not in cache: cache[s] = P.D(complex(s))
        return cache[s]
    tot = 0.0
    for a, b in zip(corners[:-1], corners[1:]):
        seg = list(np.linspace(0, 1, 9)); i = 0
        while i < len(seg) - 1:
            da, db = Dc(a + (b - a)*seg[i]), Dc(a + (b - a)*seg[i + 1])
            dph = np.angle(db/da)
            if abs(dph) > 0.3 and seg[i + 1] - seg[i] > 1e-4:
                seg.insert(i + 1, 0.5*(seg[i] + seg[i + 1])); continue
            tot += dph; i += 1
    return tot/(2*np.pi), len(cache)

if __name__ == "__main__":
    T0 = time.time()
    st = mb.star_from_shoot(0.08)
    for _ in range(3):
        st = B.solve_star(0.853087, +1, guess=st, R=st.R)
    P = PolarProblem(st, 2, LAM=0.0, MU=0.0)
    s_qnm = -5.353e-5 - 0.1195191j
    for c, lab in [(s_qnm, 'around QNM'), (s_qnm + 0.004j, 'displaced by +0.004i')]:
        n, ne = winding(P, c, 2e-4, 1.5e-3)
        print(f"{lab}: winding = {n:+.4f}  ({ne} evaluations, {time.time()-T0:.0f}s)", flush=True)

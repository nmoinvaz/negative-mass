"""
Static (s = 0) limit: tidal response of the stars and the sign that controls kappa_eff - kappa_c.

Identity (exact, any s): with Psi^ = (d_x - g) Z,  Z_xx = (s^2 + V_Z) Z,  V_Z = g^2 + g_x,  q = Z_x/Z:
      kappa = Psi^_x/Psi^ = -g + s^2/(q - g),       and  kappa_c = d_x ln phi0 = -g.
So kappa_eff - kappa_c = s^2/(q_int - g): it vanishes like s^2 for ANY star (the intertwiner maps the whole
s = 0 solution space onto phi0), and its sign is sign(q_int - g) as s -> 0.  Since psi0 is the (unique up to
an O((M/r)^{2l+1}) admixture) static solution growing like r^{l+1}, q_int - g measures the decaying
r^{-l} admixture in the star's static tidal solution, i.e. a Love-number-type coefficient:
      Z_int ∝ psi0 (1 + c_l / r^{2l+1} + ...)  =>  q_int - g ≈ -(2l+1) c_l f / R^{2l+2},
      c_l := -(q_int - g) R^{2l+2} / ((2l+1) f(R))         (should be independent of R outside the matter).
psi0 = chi(r)/(6M + c r) with the DG polynomial chi written for general M (checked: static Zerilli residual 0
for either sign of M), so the same construction applies to the canonical star.
"""
import sys, os, numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem

rs, Ms = sp.symbols('r M', real=True)
def psi0_general(l):
    s_ = 0
    for n in range(1, l + 3):
        c = sp.Integer(n - 2) * ((n - 4) * l * (l + 1) + n - 1) * sp.gamma(l + n - 1) \
            / (2**n * sp.gamma(n)**2 * sp.gamma(l - n + 3) * (-Ms)**n)
        s_ += c * rs**n
    return s_ / (6 * Ms + (l + 2) * (l - 1) * rs)

def zerilli_res(psi, l):
    lam = sp.Rational((l - 1) * (l + 2), 2); f = 1 - 2 * Ms / rs
    V = 2 * f * (lam**2 * rs**2 * ((lam + 1) * rs + 3 * Ms) + 9 * Ms**2 * (lam * rs + Ms)) / (rs**3 * (lam * rs + 3 * Ms)**2)
    return sp.simplify(-f * sp.diff(f * sp.diff(psi, rs), rs) + V * psi)

if __name__ == "__main__":
    for l in (2, 3):
        p0 = psi0_general(l)
        print(f"l={l}: static Zerilli residual of psi0(M) = {zerilli_res(p0, l)}")
        gfun = sp.lambdify((rs, Ms), sp.cancel(sp.diff(p0, rs) / p0), 'mpmath')   # d ln psi0 / dr
        for eps, w in [(-1, 0.7603), (+1, 0.75), (-1, 0.85), (+1, 0.85)]:
            st = B.solve_star(w, eps)
            for phi_cut in (1e-9, 1e-12):
                P = PolarProblem(st, l, phi_cut=phi_cut)
                zz, sv = P.interior_zerilli(0.0)
                M = P.Mext; R = P.rm; f = 1 - 2 * M / R
                q_int = zz[1] / zz[0]                 # d ln Z / dr
                g = float(gfun(R, M))
                c = -(q_int - g) * f * R**(2 * l + 2) / ((2 * l + 1) * f)
                print(f"   eps={eps:+d} w={w} M={M:+.6f} r_m={R:6.2f} (phi_cut {phi_cut:.0e}): q_int={q_int:.12e} "
                      f"g={g:.12e}  (q_int-g)/g={(q_int-g)/g:+.3e}  c_l={c:+.5e}  sign(kappa_eff-kappa_c)={'+' if q_int > g else '-'}"
                      f"   [rtol check: ", end='')
                P2 = PolarProblem(st, l, phi_cut=phi_cut, rtol=1e-12)
                zz2, _ = P2.interior_zerilli(0.0)
                print(f"{(zz2[1]/zz2[0]-g)/g:+.3e}]", flush=True)

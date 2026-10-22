"""
Effective surface (Robin) condition that the phantom star imposes on the Dotti-Gleiser variable.

At r_m (just outside the matter, phi0 < 1e-9) the exterior is vacuum Schwarzschild with M < 0.
For growth rate s = gamma > 0:
  * interior: the unique (up to scale) regular interior solution whose scalar part matches the
    decaying scalar exterior (polar_solver.interior_zerilli) -> Zerilli data (Z, dZ/dr) at r_m,
    q_int = f dZ/dr / Z = d ln Z / dx;
  * exterior: decaying Zerilli solution, q_ext (Riccati, polar_solver.zerilli_q).
A growing mode <=> q_int = q_ext (equivalently D = 0).
Mapped to the DG intertwined variable  Psi^ = (d_x - g) Z,  g = d_x ln psi0  (units |M| = 1):
    kappa = Psi^_x / Psi^ = [(s^2 + V_Z - g_x) - g q] / (q - g)          (Z_xx = (s^2 + V_Z) Z)
and compared with kappa_c(R) = d_x ln phi0 (R)  (dg_spectrum: constant-kappa Robin truncation is unstable
iff kappa < kappa_c, for R > r_s).  All quantities in |M| = 1 units (r -> r/|M|, s -> s |M|, q -> q |M|).
Usage: python kappa_eff.py l [w]  -> kappa_eff_l{l}.csv
"""
import sys, os, numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'dg_spectrum'))
import common as C
import background as B
from polar_solver import PolarProblem

l = int(sys.argv[1]); w = float(sys.argv[2]) if len(sys.argv) > 2 else 0.7603
rr = C.rs_
cc = (l + 2) * (l - 1)
f = 1 + 2 / rr
psi0 = C.chi_poly(l) / (-6 + cc * rr)
g = sp.cancel(f * sp.diff(psi0, rr) / psi0)
gx = sp.cancel(f * sp.diff(g, rr))
lam = sp.Rational(cc, 2); M = -1
VZ = 2 * f * (lam**2 * rr**2 * ((lam + 1) * rr + 3 * M) + 9 * M**2 * (lam * rr + M)) / (rr**3 * (lam * rr + 3 * M)**2)
# consistency: psi0 is a static Zerilli solution  =>  V_Z = g^2 + g_x
print("check V_Z - g^2 - g_x = 0:", sp.simplify(VZ - g**2 - gx))
gN, gxN, VN = (sp.lambdify(rr, z, 'mpmath') for z in (g, gx, VZ))
kc = C.numeric(l)['kappa_c']

def kap(q, s, R):
    G, Gx, V = float(gN(R)), float(gxN(R)), float(VN(R))
    return ((s * s + V - Gx) - G * q) / (q - G)

st = B.solve_star(w, -1)
P = PolarProblem(st, l)
Mabs = abs(P.Mext); R = P.rm / Mabs
print(f"phantom w={w} l={l}: |M|={Mabs:.6f} r_m={P.rm:.3f} R=r_m/|M|={R:.2f}  kappa_c(R)={kc(R):.6e}  -(l+1)/R={-(l+1)/R:.6e}")
rows = []
for gam in np.logspace(-4, np.log10(150), 31):
    zz, sv = P.interior_zerilli(gam)
    fm = 1 - 2 * P.Mext / P.rm
    q_int = fm * zz[1] / zz[0] * Mabs
    q_ext = P.zerilli_q(gam) * Mabs
    sM = gam * Mabs
    ki, ke = kap(q_int, sM, R), kap(q_ext, sM, R)
    rows.append((gam, sM, q_int, q_ext, ki, ke, kc(R), sv[-1] / sv[0]))
    print(f"  gamma={gam:9.3e} (s|M|={sM:.3e})  q_int={q_int:+.6e} q_ext={q_ext:+.6e}  kappa_eff={ki:+.6e}  "
          f"kappa_ext={ke:+.6e}  kappa_c={kc(R):+.4e}  (scalar-match sv ratio {sv[-1]/sv[0]:.1e})", flush=True)
np.savetxt(os.path.join(HERE, f'kappa_eff_l{l}_w{w}.csv'), np.array(rows), delimiter=',',
           header='gamma,s_absM,q_int,q_ext,kappa_eff,kappa_ext,kappa_c,svratio', comments='')

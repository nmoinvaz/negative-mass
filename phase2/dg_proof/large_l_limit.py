"""
Large-l limit of phi0, r_* and r_neg.  Scaling rho = c r, c = (l+2)(l-1)  (so r_s c = 6 exactly).
Leading order of c_n r^n with r = rho/c:  c_n ~ (n-2)(n-4) l^(2n-2)/(2^n (n-1)!^2)  [n != 4], c_4 r^4 = O(l^-4)
=> c*chi(rho/c) -> X(rho) = sum_{n>=1} (n-2)(n-4) (rho/2)^n/(n-1)!^2 = s[(s+3) I0(2 sqrt s) - 4 sqrt s I1(2 sqrt s)], s = rho/2.
phi0 -> (rho - 6)/X(rho);  near the singularity x = r^2/4 + O(r^3) = y/c^2, y = rho^2/4, H^ = c^4 [-d_y^2 + V~(y)].
Limit eigenvalue: E c^-4 -> -(1/12)^2 (k = c(c+2)/6).  r_* c -> rho_*, r_neg c -> rho_neg.
"""
import sys
import mpmath as mp
import sympy as sp
sys.path.insert(0, '/home/claude/negative-mass/phase2/dg_spectrum')
from common import chi_poly, rs_ as r
mp.mp.dps = 30

def X(rho):
    s = rho/2; z = 2*mp.sqrt(s)
    return s*((s + 3)*mp.besseli(0, z) - 4*mp.sqrt(s)*mp.besseli(1, z))

# series check
rho = sp.Symbol('rho')
ser = sum(sp.Integer((n - 2)*(n - 4))*(rho/2)**n/sp.factorial(n - 1)**2 for n in range(1, 30))
for t in [0.5, 3, 10]:
    print(f"X({t}) closed form {mp.nstr(X(t), 15)}  series {float(ser.subs(rho, t)):.15g}")
print("X > 0 on (0, 200]:", all(X(mp.mpf(j)/10) > 0 for j in range(1, 2001)))
q_ = sp.Symbol('q', positive=True)
s_ = q_/2; z_ = 2*sp.sqrt(s_)
Xs = s_*((s_ + 3)*sp.besseli(0, z_) - 4*sp.sqrt(s_)*sp.besseli(1, z_))
phis = (q_ - 6)/Xs
phiy = 2/q_*sp.diff(phis, q_)
phiyy = 2/q_*sp.diff(phiy, q_)
phi = sp.lambdify(q_, phis, 'mpmath'); dphi = sp.lambdify(q_, sp.diff(phis, q_), 'mpmath')
phi_yy = sp.lambdify(q_, phiyy, 'mpmath')
rho_star = mp.findroot(dphi, (7.5, 9.5), solver='anderson')
rho_neg = mp.findroot(phi_yy, (9, 10), solver='anderson')
rho_1 = mp.findroot(lambda q: phi_yy(q)/phi(q), (5, 6.2), solver='anderson')
print(f"rho_* = {mp.nstr(rho_star, 12)}   rho_neg = {mp.nstr(rho_neg, 12)}   (lower edge of V~<0: {mp.nstr(rho_1, 12)})")
# confirm phi_yy sign structure on (0, 60)
grid = [mp.mpf(j)/100 for j in range(1, 6000)]
sg = [mp.sign(phi_yy(q)*phi(q)) for q in grid if abs(q - 6) > 1e-3]
print("sign changes of V~ = phi_yy/phi on (0,60):", sum(1 for a, b in zip(sg, sg[1:]) if a != b))

# finite-l values (exact polynomials): r_* from N = c chi - (c r - 6) chi', r_neg from V^ numerator
print("\n  l    c*r_*    c*r_neg    l^2 r_neg")
for l in [2, 3, 5, 10, 20, 40, 80]:
    c = (l + 2)*(l - 1)
    chi = chi_poly(l)
    N = sp.Poly(sp.expand(c*chi - (c*r - 6)*sp.diff(chi, r)), r)
    rstar = max(float(iv[0][1]) for iv in N.intervals(eps=1e-13) if iv[0][0] >= 0)
    # V^ = phi0_xx / phi0 ; zero of phi0_xx : f (f phi0_r)_r = 0  -> numerator
    phi0 = (c*r - 6)/chi
    f = 1 + 2/r
    num = sp.numer(sp.together(sp.diff(f*sp.diff(phi0, r), r)))
    P = sp.Poly(sp.expand(num), r)
    roots = sorted(float(iv[0][1]) for iv in P.intervals(eps=1e-13) if iv[0][0] > 0)
    rneg = roots[-1]
    print(f"{l:3d}  {c*rstar:9.5f}  {c*rneg:9.5f}  {l*l*rneg:9.5f}   (all positive roots of phi0_xx: {[round(x*c,4) for x in roots]})", flush=True)

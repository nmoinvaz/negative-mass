"""
Exact (sympy) verification for l = 2..LMAX of every ingredient of the proofs in RESULTS.md, plus the
closed-form admissible zero-energy solution and its unique zero.
"""
import sys, time
import sympy as sp
import mpmath as mp
sys.path.insert(0, '/home/claude/negative-mass/phase2/dg_spectrum')
from common import chi_poly, rs_ as r

LMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10
mp.mp.dps = 40
zeros_dg = {2: 2.7831, 3: 1.0525, 4: 0.5719, 5: 0.3636, 6: 0.2529, 7: 0.1866, 8: 0.1436}  # dg_spectrum ODE Sturm count

def riccati_E(y_logder_r, Vh, f):
    """(H y)/y for y with d ln y/dr = y_logder_r:  -(sigma_x + sigma^2) + V^, sigma = f * dlny/dr."""
    sig = f*y_logder_r
    return sp.cancel(-(f*sp.diff(sig, r) + sig**2) + Vh)

for l in range(2, LMAX + 1):
    t0 = time.time()
    c = (l + 2)*(l - 1); k = (l - 1)*l*(l + 1)*(l + 2)//6
    f = 1 + 2/r
    chi = chi_poly(l)
    D = c*r - 6
    rs = sp.Rational(6, c)
    psi0 = chi/D
    g = sp.cancel(f*sp.diff(psi0, r)/psi0)
    Vh = sp.cancel(g**2 - f*sp.diff(g, r))
    num, den = sp.fraction(Vh)
    # 1. chi > 0 on r > 0 (exact root count) and V^ poles
    nchi = sp.Poly(chi/r, r).count_roots(0, None)
    den_f = sp.factor_list(den)
    den_ok = sp.cancel(den/(r**4*sp.expand(chi/r)**2)).is_number
    # 2. phi0 = D/chi is an exact zero mode: (H phi0)/phi0 = 0
    phi0 = D/chi
    res_phi0 = riccati_E(sp.cancel(sp.diff(phi0, r)/phi0), Vh, f)
    # 3. admissible E = 0 solution  v = phi0 * F,  F(r) = int_0^r psi0^2 (dx/dr) dr  (zero residue at r_s)
    integrand = sp.cancel(chi**2*r/(D**2*(r + 2)))
    resid = sp.residue(integrand, r, rs)
    ap = sp.apart(integrand, r)
    F = sp.integrate(ap, r)
    F = sp.expand(F - F.subs(r, 0))
    # H(phi0 F) = F H phi0 - 2 phi0_x F_x - phi0 F_xx = -2 psi0 (phi0 psi0)_x = 0 once F_x = psi0^2, i.e. dF/dr = integrand:
    res_v = sp.cancel(sp.diff(F, r) - integrand)
    v = phi0*F
    v0 = sp.nsimplify(sp.series(sp.cancel(phi0*r)*sp.expand(F), r, 0, 4).removeO().coeff(r, 3)) if l <= 4 else None
    v0 = sp.series(F, r, 0, 5).removeO()   # F = (chi'(0)^2/ (36*2)) r^4 + ...  -> v = phi0 F ~ (-4/r) * C r^4
    # zero of v beyond r_s: F(r_z) = 0
    Ff = sp.lambdify(r, F, 'mpmath')
    # F < 0 just beyond r_s (from -inf); F -> +inf; find root by bisection on a log grid
    grid = [rs*(1 + mp.mpf(10)**(-8))] + [rs*mp.mpf(1.0001)**0 * (1 + mp.mpf(j)/50) for j in range(1, 2000)]
    rz = None
    prev = Ff(grid[0])
    for a_, b_ in zip(grid[:-1], grid[1:]):
        if Ff(b_) > 0 and Ff(a_) < 0:
            rz = mp.findroot(Ff, (a_, b_), solver='anderson'); break
    sign_changes = sum(1 for a_, b_ in zip(grid[:-1], grid[1:]) if Ff(a_)*Ff(b_) < 0)
    # F > 0 on (0, r_s) ?
    Fpos = all(Ff(rs*mp.mpf(j)/200) > 0 for j in range(1, 200))
    # 4. DG eigenfunction mapped by A*: Psi^ = f (h chi' - h' chi)/(D chi),  h = r (r+2)^k e^{-k r/2}
    hl = 1/r + sp.Integer(k)/(r + 2) - sp.Rational(k, 2)          # h'/h
    Q = sp.cancel(f*(sp.diff(chi, r)/chi - hl)/D)                 # Psi^ = h * Q * chi / chi  ->  Psi^ = h*Q
    Q = sp.cancel(f*(sp.diff(chi, r) - hl*chi)/(D*chi))
    QN, QD = sp.fraction(Q)
    pole_rs = sp.Poly(QD, r).eval(rs) == 0
    EPsi = riccati_E(sp.cancel(hl + sp.diff(Q, r)/Q), Vh, f)
    ord0 = sp.Poly(QN, r).monoms()[-1][0] - sp.Poly(QD, r).monoms()[-1][0] + 1   # order of h*Q at r=0
    # 5. r_*: roots of N = c chi - D chi' on (r_s, oo);  r_neg: sign changes of V^ numerator on (0, oo)
    N = sp.Poly(sp.expand(c*chi - D*sp.diff(chi, r)), r)
    nrs = N.count_roots(rs, None)
    rstar = [float(x[0][0]) for x in sp.Poly(N, r).intervals(eps=1e-12) if x[0][0] >= rs] if nrs else []
    vroots = sp.Poly(num, r).count_roots(0, None)
    print(f"l={l:2d} k={k:5d}: #roots chi/r on (0,oo)={nchi}; den(V^)=r^4 (chi/r)^2*const: {den_ok}; "
          f"H phi0/phi0={res_phi0}; Res_rs={resid}; dF/dr-psi0^2 dx/dr={res_v}; F near 0 = {v0}; F>0 on (0,r_s): {Fpos}; "
          f"#zeros v beyond r_s (grid)={sign_changes}, r_z={float(rz):.6f} (ODE: {zeros_dg.get(l,"-")}); r_*={rstar}; "
          f"Psi^ pole at r_s: {pole_rs}; H Psi^/Psi^={EPsi} (=-(k/2)^2={-sp.Rational(k,2)**2}); Psi^~r^{ord0}; "
          f"#roots N on (r_s,oo)={nrs}; #pos roots num V^={vroots}  [{time.time()-t0:.1f}s]", flush=True)
    if l <= 3:
        print("     F(r) =", F)

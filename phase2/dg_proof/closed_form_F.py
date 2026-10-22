"""Structure of the admissible zero-energy solution v = phi0*F, F(r) = f.p. int_0^r psi0^2 (dx/dr) dr:
partial fractions of psi0^2 dx/dr = r chi^2/((c r-6)^2 (r+2)):  polynomial + A/(r-r_s)^2 + B/(r-r_s) + C/(r+2).
B = 0 (Lemma L1); check C = -1/2, i.e. chi(-2)^2 = (c+3)^2 ... and record A."""
import sys, sympy as sp
sys.path.insert(0, '/home/claude/negative-mass/phase2/dg_spectrum')
from common import chi_poly, rs_ as r
for l in range(2, 16):
    c = (l + 2)*(l - 1); chi = chi_poly(l); rs = sp.Rational(6, c)
    ig = r*chi**2/((c*r - 6)**2*(r + 2))
    B = sp.residue(ig, r, rs); C = sp.residue(ig, r, -2)
    A = sp.limit(ig*(r - rs)**2, r, rs)
    print(f"l={l:2d}: B (1/(r-r_s) coeff) = {B}, C (1/(r+2) coeff) = {C}, chi(-2) = {chi.subs(r,-2)}, A = {A} = r_s chi(r_s)^2/(c^2 (r_s+2)): {sp.simplify(A - rs*chi.subs(r,rs)**2/(c**2*(rs+2)))==0}")

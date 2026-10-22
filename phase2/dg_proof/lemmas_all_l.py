"""
All-l symbolic lemmas (l symbolic wherever possible).  M = -1, c = (l+2)(l-1) = 2 lam, f = 1 + 2/r.
 L1  chi'(r_s)/chi(r_s) = -1/(r_s(r_s+2))  (from the chi-ODE, whose leading coefficient vanishes at r_s)
     => V^ regular at r_s, and Res_{r_s} psi0^2 dx/dr = 0 (no log in the 2nd zero-energy solution)
 L2  r^4 V^ analytic at r = 0 with value 12 (Frobenius exponents 3, -1 in r for every l)
 L3  V^ ~ (l+1)(l+2)/r^2 at infinity
 L4  DG unstable mode Psi_DG = r (r+2)^k e^{-k r/2}/(c r - 6) solves Zerilli with E = -(k/2)^2, k = c(c+2)/6
     (identically in l, with k kept symbolic until the end)
 L5  DG positivity bound: q(r0) = 4(K-18)(K s - 18(K-16))/den,  (K s)^2 - (18(K-16))^2 = 6 (K-24)^3
"""
import sympy as sp

l, r, k = sp.symbols('l r k', positive=True)
c = (l + 2)*(l - 1)
lam = c/2
M = -1
f = 1 + 2/r
V = 2*f*(lam**2*r**2*((lam + 1)*r + 3*M) + 9*M**2*(lam*r + M))/(r**3*(lam*r + 3*M)**2)
P2 = -r**2*(r + 2)*(c*r - 6)
P1 = 2*r*(l**2*r**2 + 3*l**2*r + l*r**2 + 3*l*r - 2*r**2 - 6*r - 6)
P0 = l**4*r**2 + 2*l**3*r**2 - 3*l**2*r**2 - 6*l**2*r - 4*l*r**2 - 6*l*r + 4*r**2 + 12*r + 12
rs = 6/c

print("== L1 ==")
print("P2 = -r^2 (r+2)(c r - 6): P2(r_s) =", sp.simplify(P2.subs(r, rs)))
ratio = sp.simplify(-P0.subs(r, rs)/P1.subs(r, rs))
print("chi'/chi at r_s = -P0/P1 =", sp.factor(ratio), "; -1/(r_s(r_s+2)) =", sp.factor(-1/(rs*(rs + 2))),
      "; diff =", sp.simplify(ratio + 1/(rs*(rs + 2))))
print("P1(r_s) =", sp.factor(P1.subs(r, rs)), "(nonzero for l>=2, so the relation is forced; chi(r_s)>0 by Lemma A)")
# V^ = phi0_xx/phi0 regular at r_s iff phi0_xx(r_s) = 0: phi0_xx = f^2 phi0_rr + f f_r phi0_r
# with phi0 = (c r - 6)/chi:  phi0_r(r_s) = c/chi, phi0_rr(r_s) = -2 c chi'/chi^2
X, Xp = sp.symbols('X Xp')   # chi(r_s), chi'(r_s)
fr = sp.diff(f, r)
phixx = (f**2*(-2*c*Xp/X**2) + f*fr*(c/X)).subs(r, rs)
print("phi0_xx(r_s) with chi' = ratio*chi:", sp.simplify(phixx.subs(Xp, ratio*X)))
# residue of psi0^2 dx/dr = chi^2 r /((c r-6)^2 (r+2)) at r_s: (1/c^2) d/dr[chi^2 r/(r+2)]
res = (2*X*Xp*r/(r + 2) + X**2*2/(r + 2)**2).subs(r, rs)/c**2
print("Res_{r_s} psi0^2 dx/dr :", sp.simplify(res.subs(Xp, ratio*X)))

print("\n== L2: behaviour at r = 0 (needs chi = (3/2) r - a r^3 + b r^4 + e r^5 + O(r^6)) ==")
a, b, e, e6 = sp.symbols('a b e e6')
chi_s = sp.Rational(3, 2)*r - a*r**3 + b*r**4 + e*r**5 + e6*r**6
psi0 = chi_s/(c*r - 6)
g = f*sp.diff(psi0, r)/psi0
Vh = g**2 - f*sp.diff(g, r)
ser = sp.series(sp.simplify(r**4*Vh), r, 0, 3).removeO()
print("r^4 V^ = ", sp.simplify(ser.coeff(r, 0)), "+ (", sp.factor(ser.coeff(r, 1)), ") r + O(r^2)")
print("-> r^4 V^ analytic at 0 (rational with chi = r*(3/2+...)), value 12; indicial eq. in r:")
s = sp.Symbol('s')
# Psi_rr + (f_r/f) Psi_r - (V^-E)/f^2 Psi = 0; r f_r/f -> -1, r^2 V^/f^2 -> 12/4 = 3
print("   s(s-1) - s - 3 =", sp.factor(s*(s - 1) - s - 3), " -> exponents 3, -1 for every l")

print("\n== L3: infinity ==")
n_ = sp.Symbol('n', positive=True)
# leading behaviour: psi0 ~ C r^{l+1}(1 + d/r + ...): g = f psi0_r/psi0
d1 = sp.Symbol('d1')
psiinf = r**(l + 1)*(1 + d1/r)
ginf = f*sp.diff(psiinf, r)/psiinf
Vhinf = sp.simplify(2*ginf**2 - V)
lim = sp.limit(sp.simplify(r**2*Vhinf), r, sp.oo)
print("lim r^2 V^ =", sp.factor(lim), " (V^ = 2 g^2 - V)")

print("\n== L4: DG unstable mode, all l ==")
Psi = r*(r + 2)**k*sp.exp(-k*r/2)/(c*r - 6)
Ez = -f*sp.diff(f*sp.diff(Psi, r), r) + V*Psi
q = sp.simplify(sp.expand(sp.simplify(Ez/Psi)))
q = sp.factor(sp.simplify(q.subs(k, c*(c + 2)/6)))
print("(H_Z Psi)/Psi with k = c(c+2)/6 :", q, " ; -(k/2)^2 =", sp.factor((c*(c + 2)/12)**2*-1))
print("difference:", sp.simplify(q + (c*(c + 2)/12)**2))
print("k = c(c+2)/6 = (l-1)l(l+1)(l+2)/6 :", sp.expand(c*(c + 2)/6 - (l - 1)*l*(l + 1)*(l + 2)/6) == 0)

print("\n== L5: DG eq. 59 inequality ==")
K, s_ = sp.symbols('K s', positive=True)
print("(K s)^2 - 18^2 (K-16)^2 with s^2 = 6(K-18):", sp.factor(sp.expand(K**2*6*(K - 18) - 324*(K - 16)**2)))
print("K = (l-1)l(l+1)(l+2) >= 120 for l >= 3  ->  q(r0) > 0;  K = 24 (l = 2) is exactly the degenerate case.")

print("\n== L6: DG mode mapped by A* = -d/dx + g is regular at r_s and admissible at r = 0 (all l) ==")
# Psi^ = A* Psi_DG = phi0 * W_x(Psi_DG, psi0) = f (h chi' - h' chi)/(chi (c r - 6)),  h = r (r+2)^k e^{-kr/2}
kk = c*(c + 2)/6
lnh_p = 1/r + kk/(r + 2) - kk/2
print("(ln h)'(r_s) - (ln chi)'(r_s) =", sp.simplify(lnh_p.subs(r, rs) + 1/(rs*(rs + 2))),
      " -> numerator h chi' - h' chi vanishes at r_s: Psi^ regular there")
a_ = (l - 1)*l*(l + 1)*(l + 2)/32; b_ = (l - 1)*l*(l + 1)*(l + 2)/96
X1 = sp.Rational(3, 2) - a_*r**2 + b_*r**3 + sp.Symbol('e5')*r**4
H1p_H1 = kk/(r + 2) - kk/2
diff_ser = sp.series(sp.diff(X1, r)/X1 - H1p_H1, r, 0, 3).removeO()
print("X1'/X1 - H1'/H1 = ", sp.simplify(diff_ser.coeff(r, 0)), "+", sp.simplify(diff_ser.coeff(r, 1)), "* r +", sp.simplify(diff_ser.coeff(r, 2)), "* r^2 + O(r^3)")
print(" -> h chi' - h' chi = r^2 X1 H1 (X1'/X1 - H1'/H1) = O(r^5), Psi^ = O(r^3): admissible branch (b0 = 0).")

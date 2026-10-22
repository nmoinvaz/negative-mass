"""
Linearized RADIAL (l = 0) perturbations of a signed-source complex-scalar boson star.

    G_ab = 8 pi eps T_ab,   Box Phi = 2 V'(|Phi|^2) Phi,
    ds^2 = -e^nu dt^2 + e^lam dr^2 + r^2 dOmega^2      (alpha^2 = e^nu, a^2 = e^lam),
    T_ab = Re(d_a Phi^* d_b Phi) - g_ab [ 1/2 |dPhi|^2 + V ].

Background: nu0(r), lam0(r), Phi = phi0(r) e^{-i w t}.
Perturbation (eps_p -> 0):
    nu  = nu0  + eps_p n(t,r),   lam = lam0 + eps_p l(t,r),
    Phi = [phi0 + eps_p (F(t,r) + i G(t,r))] e^{-i w t}.

Everything is derived here from the metric with eps kept symbolic; the script

  1. builds G_ab - 8 pi eps T_ab and the KG equation, expands to first order;
  2. checks that the background equations used in boson_star.py are reproduced;
  3. goes to the frequency domain  X(t,r) = X(r) e^{i sigma t}  and performs the reduction
         (tr)  -> l algebraic,
         (KG_I) = current conservation -> charge-perturbation function Psi,
         gauge-invariant  Gt = i sigma G + w phi0 n / 2,
         (rr)  -> n' algebraic,
     leaving a closed system for (F, Gt, Psi) that contains sigma only as sigma^2;
  4. verifies symbolically that the reduced system reproduces ALL five linearized equations
     (tt, tr, rr, KG_R, KG_I) identically (residual 0), so nothing was lost or mis-signed;
  5. prints the final equations, and writes them (as sympy srepr) to radial_eqs.txt for reference.

Run:  python derive_radial.py
"""
import sympy as sp

t, r, th = sp.symbols('t r theta', real=True)
w, eps, sig = sp.symbols('omega epsilon sigma')
LAM, MU = sp.symbols('Lambda mu', positive=True)
e = sp.Symbol('e_p')                      # bookkeeping parameter for linearization

nu0 = sp.Function('nu0')(r); lam0 = sp.Function('lam0')(r); p0 = sp.Function('phi0')(r)
n = sp.Function('n')(t, r); l = sp.Function('l')(t, r)
F = sp.Function('F')(t, r); G = sp.Function('G')(t, r)

def Vf(s):  return s/2 - LAM*s**2/4 + MU*s**3/6
def dV2(s): return 1 - LAM*s + MU*s**2          # 2 dV/ds

nu = nu0 + e*n; lam = lam0 + e*l
U = p0 + e*F; W = e*G                              # Z = U + i W,  Phi = Z e^{-i w t}

x = [t, r, th, sp.Symbol('ph')]
g = sp.diag(-sp.exp(nu), sp.exp(lam), r**2, r**2*sp.sin(th)**2)
ginv = g.inv()

def christoffel(g, ginv):
    N = 4
    return [[[sum(ginv[a, d]*(sp.diff(g[d, b], x[c]) + sp.diff(g[d, c], x[b]) - sp.diff(g[b, c], x[d]))
                  for d in range(N))/2 for c in range(N)] for b in range(N)] for a in range(N)]

Gam = christoffel(g, ginv)
def ricci(a, b):
    return sp.simplify(sum(sp.diff(Gam[c][a][b], x[c]) - sp.diff(Gam[c][a][c], x[b])
                           + sum(Gam[c][c][d]*Gam[d][a][b] - Gam[c][b][d]*Gam[d][a][c] for d in range(4))
                           for c in range(4)))
Ric = sp.zeros(4, 4)
for a in range(4):
    for b in range(a, 4):
        Ric[a, b] = ricci(a, b); Ric[b, a] = Ric[a, b]
Rs = sp.simplify(sum(ginv[a, b]*Ric[a, b] for a in range(4) for b in range(4)))
Ein = sp.simplify(Ric - Rs*g/2)

# covariant derivative of Z with the phase removed: D_t = d_t - i w, D_r = d_r
DZ = {0: (sp.diff(U, t) + w*W, sp.diff(W, t) - w*U),   # (Re, Im)
      1: (sp.diff(U, r), sp.diff(W, r))}
def ReDD(a, b):
    if a > 1 or b > 1: return 0
    return DZ[a][0]*DZ[b][0] + DZ[a][1]*DZ[b][1]
kin = sum(ginv[a, b]*ReDD(a, b) for a in range(2) for b in range(2))   # |dPhi|^2
s = U**2 + W**2
T = sp.zeros(4, 4)
for a in range(4):
    for b in range(4):
        T[a, b] = ReDD(a, b) - g[a, b]*(kin/2 + Vf(s))

sqrtg = sp.exp((nu + lam)/2)*r**2
# Box Phi e^{iwt} - 2V' Z  (real and imaginary parts)
def box_part(k):   # k = 0 real, 1 imag
    # D_t acting on (X_R + i X_I): Re -> dX_R/dt + w X_I ; Im -> dX_I/dt - w X_R
    XtR = sqrtg*ginv[0, 0]*DZ[0][0]; XtI = sqrtg*ginv[0, 0]*DZ[0][1]
    if k == 0: tt_ = sp.diff(XtR, t) + w*XtI
    else:      tt_ = sp.diff(XtI, t) - w*XtR
    rr_ = sp.diff(sqrtg*ginv[1, 1]*DZ[1][k], r)
    return (tt_ + rr_)/sqrtg - dV2(s)*(U if k == 0 else W)

eqs = {'tt': Ein[0, 0] - 8*sp.pi*eps*T[0, 0], 'tr': Ein[0, 1] - 8*sp.pi*eps*T[0, 1],
       'rr': Ein[1, 1] - 8*sp.pi*eps*T[1, 1], 'thth': Ein[2, 2] - 8*sp.pi*eps*T[2, 2],
       'KGR': box_part(0), 'KGI': box_part(1)}

# ---------------- background ----------------
bg = {k: sp.simplify(v.subs(e, 0).doit()) for k, v in eqs.items()}
# Solve tt, rr, KGR for lam0', nu0', phi0''
sol_bg = sp.solve([bg['tt'], bg['rr'], bg['KGR']], [sp.diff(lam0, r), sp.diff(nu0, r), sp.diff(p0, r, 2)], dict=True)[0]
sol_bg = {k: sp.simplify(v) for k, v in sol_bg.items()}
# second derivatives of nu0, lam0 (needed after substitution)
dnu0 = sol_bg[sp.diff(nu0, r)]; dlam0 = sol_bg[sp.diff(lam0, r)]
ddp0 = sol_bg[sp.diff(p0, r, 2)]

def bgsub(ex):
    """Replace nu0'', lam0'', phi0'', nu0', lam0' using the background equations."""
    for _ in range(3):
        ex = ex.subs({sp.diff(nu0, r, 2): sp.diff(dnu0, r), sp.diff(lam0, r, 2): sp.diff(dlam0, r)})
        ex = ex.subs(sp.diff(p0, r, 2), ddp0)
        ex = ex.subs({sp.diff(nu0, r): dnu0, sp.diff(lam0, r): dlam0})
    return ex

# check against boson_star.py: m' = 4 pi r^2 eps rho, alpha'/alpha = (a^2-1)/(2r) + 4 pi r a^2 eps p_r, etc.
a2 = sp.exp(lam0); al2 = sp.exp(nu0)
rho = (w**2*p0**2/al2 + sp.diff(p0, r)**2/a2)/2 + Vf(p0**2)
pr = (w**2*p0**2/al2 + sp.diff(p0, r)**2/a2)/2 - Vf(p0**2)
chk_nu = sp.simplify(dnu0/2 - ((a2 - 1)/(2*r) + 4*sp.pi*r*a2*eps*pr))
chk_lam = sp.simplify(dlam0/2 - ((1 - a2)/(2*r) + 4*sp.pi*r*a2*eps*rho))
chk_phi = sp.simplify(ddp0 - (-(2/r + dnu0/2 - dlam0/2)*sp.diff(p0, r) + a2*(dV2(p0**2) - w**2/al2)*p0))
print("background vs boson_star.py (should be 0 0 0):", chk_nu, chk_lam, chk_phi)
print("theta-theta background consistent with Bianchi:", sp.simplify(bgsub(bg['thth'])))

# ---------------- linearization ----------------
lin = {k: sp.expand(sp.diff(v, e).subs(e, 0).doit()) for k, v in eqs.items()}

# frequency domain
A, B, N, L = [sp.Function(nm)(r) for nm in ('A', 'B', 'N', 'L')]
ph = sp.exp(sp.I*sig*t)
fd = {F: A*ph, G: B*ph, n: N*ph, l: L*ph}
def to_fd(ex):
    ex = ex.subs(fd).doit()
    return sp.expand(sp.simplify(ex/ph))
linfd = {k: bgsub(to_fd(v)) for k, v in lin.items()}

# ---- reduction ----
# (tr): solve for L
Lsol = sp.solve(linfd['tr'], L)[0]
print("\n(tr) gives  L =", sp.simplify(Lsol))

# Current conservation: KG_I * phi0 * sqrt(-g) should be d_t(J^t) + d_r(J^r). Define
#   i sig Psi = r^2 e^{(nu0-lam0)/2} (phi0 B' - phi0' B)          [J^r part]
# and the gauge invariant  Gt = i sig B + w phi0 N / 2.
Psi, Gt = sp.Function('Psi')(r), sp.Function('Gt')(r)
E = sp.exp((nu0 - lam0)/2)
Bexpr = (Gt - w*p0*N/2)/(sp.I*sig)

# Express L in terms of Psi: the B-dependence of L is through (phi0 B' - phi0' B) = i sig Psi / (r^2 E)
Bp_from_Psi = sp.solve(sp.Eq(sp.I*sig*Psi, r**2*E*(p0*sp.diff(B, r) - sp.diff(p0, r)*B)), sp.diff(B, r))[0]
L_red = sp.simplify(Lsol.subs(sp.diff(B, r), Bp_from_Psi))
print("\nL in reduced variables:  L =", L_red, "   (has B?", L_red.has(B), ")")

# Current conservation (KG_I x phi0 sqrt(-g)) integrated once:
Psi_p_claim = r**2/E*p0*(Gt - w*p0*L_red/2 - 2*w*A)

def red(ex, Np=None):
    """Substitute L -> L_red, B -> Bexpr, Psi' -> Psi_p_claim (and N' -> Np if given)."""
    ex = ex.subs(L, L_red).doit()
    ex = ex.subs(B, Bexpr).doit()
    for _ in range(2):
        ex = ex.subs(sp.diff(Psi, r, 2), sp.diff(Psi_p_claim, r)).subs(sp.diff(Psi, r), Psi_p_claim)
        if Np is not None:
            ex = ex.subs(sp.diff(N, r, 2), sp.diff(Np, r)).subs(sp.diff(N, r), Np)
        ex = bgsub(ex)
    return sp.expand(ex)

rr_sub = red(linfd['rr'])
print("\n(rr) contains Gt', B?", rr_sub.has(sp.diff(Gt, r)), rr_sub.has(B))
Np = sp.simplify(sp.solve(rr_sub, sp.diff(N, r))[0])
print("(rr) gives  N' =", Np)
print("   N appears undifferentiated in N'? (must be False: gauge)", Np.has(N))

# Psi definition, rewritten with Gt: -sig^2 Psi = r^2 E (phi0 Gt' - phi0' Gt - w phi0^2 N'/2)
Gtp_claim = sp.simplify(sp.solve(sp.Eq(-sig**2*Psi, r**2*E*(p0*sp.diff(Gt, r) - sp.diff(p0, r)*Gt
                                                         - w*p0**2*Np/2)), sp.diff(Gt, r))[0])
KGR_sub = red(linfd['KGR'], Np)
KGR_sub = sp.expand(KGR_sub.subs(sp.diff(Gt, r), Gtp_claim))
print("\nKG_R reduced: has B?", KGR_sub.has(B), " has N?", KGR_sub.has(N), " has Gt'?", KGR_sub.has(sp.diff(Gt, r)))
App = sp.simplify(sp.solve(KGR_sub, sp.diff(A, r, 2))[0])
print("   A'' even in sigma?", sp.simplify(App.subs(sig, -sig) - App) == 0)

def impose(ex):
    ex = red(ex, Np)
    for _ in range(3):
        ex = ex.subs(sp.diff(Gt, r, 2), sp.diff(Gtp_claim, r)).subs(sp.diff(Gt, r), Gtp_claim)
        ex = ex.subs(sp.diff(Psi, r), Psi_p_claim)
        ex = ex.subs(sp.diff(N, r), Np)
        ex = ex.subs(sp.diff(A, r, 3), sp.diff(App, r)).subs(sp.diff(A, r, 2), App)
        ex = bgsub(ex)
    return sp.simplify(ex)

print("\nVerification: residuals of ALL linearized field equations on solutions of the reduced system")
for k in ('tr', 'rr', 'KGR', 'KGI', 'tt', 'thth'):
    print(f"   residual[{k}] =", impose(linfd[k]))

# ---- print the final system ----
print("\n================ FINAL RADIAL SYSTEM (sigma enters only as sigma^2) ================")
print("Variables: A = delta Re(phi) amplitude, Gt = i sigma B + w phi0 N/2 (gauge-invariant),")
print("           Psi = perturbation of enclosed Noether charge (i sig Psi = r^2 e^{(nu0-lam0)/2}(phi0 B' - phi0' B)).")
print("Notation: E = e^{(nu0-lam0)/2} = alpha/a.\n")
print("L   =", L_red)
print("\nN'  =", sp.simplify(Np.subs(L, L_red)))
print("\nPsi' =", sp.simplify(Psi_p_claim))
print("\n-sigma^2 Psi = r^2 E (phi0 Gt' - phi0' Gt - w phi0^2 N'/2)")
print("\nA'' =", sp.simplify(App))

with open(__file__.replace('derive_radial.py', 'radial_eqs.txt'), 'w') as fh:
    fh.write("# sympy srepr of the reduced radial system; see derive_radial.py\n")
    for nm, ex in [('L', L_red), ('Np', sp.simplify(Np.subs(L, L_red))), ('Psip', Psi_p_claim),
                   ('Gtp', Gtp_claim), ('App', App)]:
        fh.write(f"{nm} = {sp.simplify(ex)}\n")

# ================= code generation for the numerical solver =================
# Use chi = Psi/phi0 (Psi decays like phi0 * B, so chi stays O(B) where phi0 is exponentially small)
# and u = phi0'/phi0. State y = (A, A', Gt, chi), linear in sigma^2.
chi = sp.Function('chi')(r)
rs, al, aa, P0, P0p, Uu, Ay, Apy, Gty, chy, s2 = sp.symbols('r alpha a phi0 dphi0 u A dA Gt chi s2')
def numerize(ex):
    ex = ex.subs(Psi, p0*chi).doit()
    ex = ex.subs(sp.diff(A, r), Apy).subs(A, Ay).subs(Gt, Gty).subs(chi, chy)
    ex = ex.subs(sp.diff(p0, r), Uu*p0)                 # phi0' = u phi0
    ex = ex.subs({sp.exp(nu0): al**2, sp.exp(lam0): aa**2})
    ex = ex.subs({nu0: 2*sp.log(al), lam0: 2*sp.log(aa)})
    ex = ex.subs(sig**2, s2).subs(p0, P0)
    return sp.simplify(sp.powsimp(sp.expand_power_exp(ex), force=True))
Np_n = numerize(Np)
Psip_over_phi = numerize(Psi_p_claim/p0)
chip = sp.simplify(Psip_over_phi - Uu*chy)
Gtp = sp.simplify(Uu*Gty + w*P0*Np_n/2 - s2*chy/(rs**2*al/aa))
# check against the Gt' relation derived above
chk = numerize(Gtp_claim) - Gtp
print("\nGt' codegen check (0):", sp.simplify(chk))
Appn = numerize(App)
print("dependence check: A'' has sigma^2 only linearly:", sp.degree(sp.expand(Appn), s2))
from sympy.printing.numpy import NumPyPrinter
pr = NumPyPrinter({'fully_qualified_modules': False, 'allow_unknown_functions': True})
with open(__file__.replace('derive_radial.py', 'radial_rhs_generated.py'), 'w') as fh:
    fh.write('"""AUTO-GENERATED by derive_radial.py -- do not edit.\n'
             'Reduced l=0 system, y=(A, dA, Gt, chi), chi = Psi/phi0, u = phi0\'/phi0.\n'
             'Background: alpha, a (metric), phi0, u; omega code frequency (alpha(0)=1); s2 = sigma^2."""\n')
    fh.write('from numpy import pi, exp, sqrt, log\n\n')
    args = 'r, alpha, a, phi0, u, A, dA, Gt, chi, omega, epsilon, s2, Lambda, mu'
    for nm, ex in [('Np', Np_n), ('ddA', Appn), ('dGt', Gtp), ('dchi', chip),
                   ('Lpert', numerize(L_red))]:
        fh.write(f'def {nm}({args}):\n    return {pr.doprint(ex)}\n\n')
print("wrote radial_rhs_generated.py")

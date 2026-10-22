"""
Phase 0 checks for negative-mass Schwarzschild even-parity (Zerilli) perturbations.

Sources checked:
  GHI  = Gibbons, Hartnoll, Ishibashi, hep-th/0409307  (claim: stable for a chosen BC, q_C = 2 + log(98/9))
  GD06 = Gleiser, Dotti, gr-qc/0604021                 (claim: unstable)
  DG08 = Dotti, Gleiser, arXiv:0809.3615               (IVP via intertwiner; instability excitable)

Conventions (GD06/DG08): f = 1 - 2M/r, dr*/dr = 1/f, lam = (l-1)(l+2)/2,
Zerilli:  d^2 Psi/dt^2 - d^2 Psi/dr*^2 + V Psi = 0,
V = 2 f [lam^2 r^2((lam+1) r + 3M) + 9 M^2 (lam r + M)] / [r^3 (lam r + 3M)^2].
Negative mass: M < 0, mu = -2M > 0.
"""
import sympy as sp

r, t, M, mu = sp.symbols('r t M mu', real=True)
l = sp.Symbol('l', positive=True, integer=True)

def zerilli_V(lam, Mv):
    f = 1 - 2*Mv/r
    return 2*f*(lam**2*r**2*((lam+1)*r + 3*Mv) + 9*Mv**2*(lam*r + Mv)) / (r**3*(lam*r + 3*Mv)**2)

def zerilli_residual(psi, lam, Mv, omega2):
    """-f d/dr(f dpsi/dr) + V psi - omega2 psi  (should vanish)."""
    f = 1 - 2*Mv/r
    return sp.simplify(-f*sp.diff(f*sp.diff(psi, r), r) + zerilli_V(lam, Mv)*psi - omega2*psi)

print("=== A. Do the published exact solutions solve the Zerilli equation? ===")

# A1: GHI eq.(19), l=2 static (omega^2=0) solution, written with mu = -2M.
lam2 = sp.Integer(2)  # l=2 -> lam = (1)(4)/2 = 2
Mneg = -mu/2
sol3 = r*(3*mu**3 - 6*r**2*mu + 4*r**3) / (mu**3*(4*r - 3*mu))
sol4 = (r*(13*mu**3 - 24*mu**2*r + 12*mu*r**2)) / (3*mu**3*(4*r - 3*mu)) \
       - sol3*sp.log((r + mu)/r)
for name, s in [("GHI C3 piece", sol3), ("GHI C4 piece", sol4)]:
    res = zerilli_residual(s, lam2, Mneg, 0)
    print(f"  {name:14s}: residual = {res}")

# A2: DG08 eq.(25) psi0 for l=2 (general M) -- the intertwiner seed
psi0 = r*(r**3 + 3*M*r**2 - 6*M**3) / (8*M**4*(3*M + 2*r))
print(f"  DG08 psi0 (l=2): residual = {zerilli_residual(psi0, lam2, M, 0)}")

# A3: GD06/DG08 unstable mode, general l, with M<0 (|M| = -M).
#     Psi = r (r-2M)^k / (2 lam r + 6M) * exp(k (t - r) / (2|M|)),  k = (l-1)l(l+1)(l+2)/6
for lv in [2, 3, 4, 5]:
    lamv = sp.Rational((lv-1)*(lv+2), 2)
    kv = sp.Rational((lv-1)*lv*(lv+1)*(lv+2), 6)
    Mn = -sp.Symbol('m', positive=True)  # M = -m, m>0
    absM = -Mn
    spatial = r*(r - 2*Mn)**kv / (2*lamv*r + 6*Mn) * sp.exp(-kv*r/(2*absM))
    growth = kv/(2*absM)
    # e^{growth t}: d^2/dt^2 -> growth^2, so Zerilli gives -d2/dr*2 + V = -growth^2  => omega^2 = -growth^2
    res = zerilli_residual(spatial, lamv, Mn, -growth**2)
    print(f"  GD unstable mode l={lv}: k={kv}, growth rate = {growth}, residual = {res}")

print()
print("=== B. GHI critical boundary condition q_C ===")
# GHI pick C3 = C4 [log(7/3) - 4/9] so the mode is regular at r_C = 3mu/4 (for mu=1).
# Check: is the combination finite at r_C = 3 mu / 4?
muv = sp.Integer(1)
s3 = sol3.subs(mu, muv); s4 = sol4.subs(mu, muv)
C3_over_C4 = sp.log(sp.Rational(7, 3)) - sp.Rational(4, 9)
marg = C3_over_C4*s3 + s4
rC = sp.Rational(3, 4)
pole = sp.limit(sp.simplify(marg*(4*r - 3)), r, rC)
print(f"  residue of (C3*sol3 + C4*sol4) at r_C=3/4 (should be 0): {sp.nsimplify(sp.N(pole, 30), rational=False)}  ~ {sp.N(pole, 15)}")
# Leading behaviour at r_C: GHI say Phi ~ (r - r_C)^2
ser = sp.series(marg, r, rC, 4).removeO()
print("  Taylor at r_C:", [sp.N(ser.coeff(r - rC, n) if n else ser.subs(r, rC), 10) for n in range(0, 4)])

# Expand near r=0 in terms of r* (mu=1): r* = r - log(1+r)
rs = sp.Symbol('rs', positive=True)
# r*(r) ~ r^2/2 - r^3/3 + ... ; invert numerically via asymptotic matching below.
# Behaviour Phi ~ a0 r*^{1/2} log r* + b0 r*^{1/2};  q = b0/a0 with log(r*/mu).
import mpmath as mp
mp.mp.dps = 40
margf = sp.lambdify(r, marg, 'mpmath')
def rstar(rv):
    return rv - mp.log(1 + rv)
def fit_q(r1, r2):
    # Solve Phi(r) = a0 sqrt(x) log x + b0 sqrt(x) at two small radii
    x1, x2 = rstar(r1), rstar(r2)
    A = mp.matrix([[mp.sqrt(x1)*mp.log(x1), mp.sqrt(x1)], [mp.sqrt(x2)*mp.log(x2), mp.sqrt(x2)]])
    b = mp.matrix([margf(r1), margf(r2)])
    a0, b0 = mp.lu_solve(A, b)
    return b0/a0
for (r1, r2) in [(mp.mpf('1e-4'), mp.mpf('2e-4')), (mp.mpf('1e-7'), mp.mpf('2e-7')), (mp.mpf('1e-10'), mp.mpf('2e-10'))]:
    print(f"  fitted q at r~{mp.nstr(r1,3)}: {mp.nstr(fit_q(r1, r2), 15)}")
print(f"  GHI q_C = 2 + log(98/9)       : {mp.nstr(2 + mp.log(mp.mpf(98)/9), 15)}")

print()
print("=== C. Is GHI's truncated marginal mode a smooth metric perturbation? ===")
# GHI: Phi = marg for r < r_C, Phi = 0 for r > r_C (static, so H1 = 0).
# Metric functions from Zerilli function (DG08 eq.6), M = -1/2 (mu=1), lam=2:
Mv = -sp.Rational(1, 2); lamv = 2
f = 1 - 2*Mv/r
qf = (lamv*(lamv+1)*r**2 + 3*lamv*Mv*r + 6*Mv**2) / (r**2*(lamv*r + 3*Mv))
hf = (lamv*r**2 - 3*lamv*r*Mv - 3*Mv**2) / ((r - 2*Mv)*(lamv*r + 3*Mv))
def metric_from(P):
    K = qf*P + f*sp.diff(P, r)
    H2 = sp.diff(f*(hf*P + r*sp.diff(P, r)), r) - K
    return K, H2
K_in, H2_in = metric_from(marg)
# Left limits at r_C (right side is identically zero):
for name, expr in [("K", K_in), ("dK/dr", sp.diff(K_in, r)), ("H2", H2_in), ("dH2/dr", sp.diff(H2_in, r))]:
    val = sp.limit(expr, r, rC, '-')
    print(f"  left limit of {name:7s} at r_C: {sp.N(val, 12)}   (right side = 0)")

print()
print("=== D. Lemma 1 (DG08 eq.18) on the unstable l=2 mode: smooth metric <=> Laurent c0 = lam^2 c_{-1} / (3M(3+2 lam)) ===")
Mn = sp.Integer(-1); lamv = 2; kv = 4
spatial = r*(r - 2*Mn)**kv / (2*lamv*r + 6*Mn) * sp.exp(-kv*r/(2*(-Mn)))
rsing = -3*Mn/lamv
eps = sp.Symbol('eps')
ls = sp.series(spatial.subs(r, rsing + eps), eps, 0, 2).removeO()
cm1 = ls.coeff(eps, -1); c0 = ls.coeff(eps, 0)
print(f"  c_-1 = {sp.N(cm1, 12)}, c_0 = {sp.N(c0, 12)}, predicted c_0 = {sp.N(lamv**2*cm1/(3*Mn*(3+2*lamv)), 12)}")

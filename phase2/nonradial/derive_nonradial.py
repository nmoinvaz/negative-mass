"""
Linear NON-RADIAL (l >= 2) perturbations of a signed-source complex-scalar boson star.

    G_ab = 8 pi eps T_ab,   Box Phi = 2 V'(|Phi|^2) Phi,   V(s) = s/2 - Lam s^2/4 + mu s^3/6,
    T_ab = Re(D_a Z^* D_b Z) - g_ab [ 1/2 g^cd Re(D_c Z^* D_d Z) + V ],   Phi = Z e^{-i w t},
    background  ds^2 = -e^nu dt^2 + e^lam dr^2 + r^2 dOmega^2,  Z0 = phi0(r)  (same as ../radial).

Time dependence e^{s t} (s = -i sigma; s = gamma > 0 real is a purely growing mode; s complex allowed).
Angular dependence: axisymmetric m = 0, Y = P_l(cos theta) -- enough by spherical symmetry.

POLAR (Regge-Wheeler gauge):
    h_tt = e^nu H0 Y,  h_tr = H1 Y,  h_rr = e^lam H2 Y,  h_thth = r^2 K Y,  h_phph = r^2 sin^2 K Y,
    delta Z = (F + i G) Y         (F, G complex amplitudes if s is complex; real if s is real)
AXIAL (Regge-Wheeler gauge):
    h_tph = h0 sin(th) dY/dth,   h_rph = h1 sin(th) dY/dth,   delta Z = 0 (a scalar has no axial part).

The script
  1. builds delta(G_ab - 8 pi eps T_ab) and delta(KG) to first order, eps symbolic;
  2. separates the angular dependence;
  3. polar: (thth - phph) gives H2 = H0; (t th), (t r), (r th) give H1', K', H0';
     saves the separated equations E[...] to polar_eqs.pkl for reduce_polar.py.
  (axial: derive_axial.py)
Run:  python derive_nonradial.py   (takes a few minutes)
"""
import sympy as sp
import pickle, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()

t, r, th, ph = sp.symbols('t r theta phi', real=True)
w, eps, s = sp.symbols('omega epsilon s')
LAM, MU = sp.symbols('Lambda mu', positive=True)
LL = sp.Symbol('LL')               # l(l+1)
e = sp.Symbol('e_p')
X = [t, r, th, ph]

nu0 = sp.Function('nu0')(r); lam0 = sp.Function('lam0')(r); p0 = sp.Function('phi0')(r)
Yf = sp.Function('Y')(th)
Ysym, Ypsym = sp.symbols('Ysym Ypsym')

def Vf(x):  return x/2 - LAM*x**2/4 + MU*x**3/6
def dV2(x): return 1 - LAM*x + MU*x**2

ex = sp.exp(s*t)

def linearize(hmat, dZR, dZI):
    """Return dict of linearized equations (multiplied out, with e^{st} removed)."""
    g0 = sp.diag(-sp.exp(nu0), sp.exp(lam0), r**2, r**2*sp.sin(th)**2)
    g = g0 + e*hmat
    g0i = g0.inv()
    gi = g0i - e*g0i*hmat*g0i                      # exact to first order
    N = 4
    def Gam(gm, gim):
        return [[[sum(gim[a, d]*(sp.diff(gm[d, b], X[c]) + sp.diff(gm[d, c], X[b]) - sp.diff(gm[b, c], X[d]))
                      for d in range(N))/2 for c in range(N)] for b in range(N)] for a in range(N)]
    Gfull = Gam(g, gi)
    G0 = [[[sp.simplify(Gfull[a][b][c].subs(e, 0)) for c in range(N)] for b in range(N)] for a in range(N)]
    dG = [[[sp.expand(sp.diff(Gfull[a][b][c], e).subs(e, 0)) for c in range(N)] for b in range(N)] for a in range(N)]
    def ric(Gm):
        R = sp.zeros(N, N)
        for a in range(N):
            for b in range(N):
                R[a, b] = sum(sp.diff(Gm[c][a][b], X[c]) - sp.diff(Gm[c][a][c], X[b])
                              + sum(Gm[c][c][d]*Gm[d][a][b] - Gm[c][b][d]*Gm[d][a][c] for d in range(N))
                              for c in range(N))
        return R
    R0 = ric(G0).applyfunc(sp.simplify)
    # delta R_ab (Palatini-like, written out with partials)
    dR = sp.zeros(N, N)
    for a in range(N):
        for b in range(N):
            dR[a, b] = sum(sp.diff(dG[c][a][b], X[c]) - sp.diff(dG[c][a][c], X[b])
                           + sum(dG[c][c][d]*G0[d][a][b] + G0[c][c][d]*dG[d][a][b]
                                 - dG[c][b][d]*G0[d][a][c] - G0[c][b][d]*dG[d][a][c] for d in range(N))
                           for c in range(N))
    # matter: Z = phi0 + e (dZR + i dZI);  D_t = d_t - i w
    U = p0 + e*dZR; W = e*dZI
    DZ = [(sp.diff(U, t) + w*W, sp.diff(W, t) - w*U)] + [(sp.diff(U, X[k]), sp.diff(W, X[k])) for k in (1, 2, 3)]
    def ReDD(a, b): return DZ[a][0]*DZ[b][0] + DZ[a][1]*DZ[b][1]
    kin = sum(gi[a, b]*ReDD(a, b) for a in range(N) for b in range(N))
    sq = U**2 + W**2
    T = sp.Matrix(N, N, lambda a, b: ReDD(a, b) - g[a, b]*(kin/2 + Vf(sq)))
    trT = sum(gi[a, b]*T[a, b] for a in range(N) for b in range(N))
    Sb = T - g*trT/2                               # trace-reversed source: R_ab = 8 pi eps (T_ab - g T/2)
    dS = Sb.applyfunc(lambda z: sp.expand(sp.diff(z, e).subs(e, 0)))
    S0 = Sb.applyfunc(lambda z: sp.simplify(z.subs(e, 0)))
    # KG: (1/sqrt-g) (d_a - i w delta_a^t)(sqrt-g g^ab D_b Z) - 2V' Z
    sqrtg0 = sp.exp((nu0 + lam0)/2)*r**2*sp.sin(th)
    # sqrt(-g) to first order: sqrtg0 (1 + e/2 g0^ab h_ab)
    trh = sum(g0i[a, b]*hmat[a, b] for a in range(N) for b in range(N))
    sqrtg = sqrtg0*(1 + e*trh/2)
    def KG(k):
        out = 0
        for a in range(N):
            XR = sqrtg*sum(gi[a, b]*DZ[b][0] for b in range(N))
            XI = sqrtg*sum(gi[a, b]*DZ[b][1] for b in range(N))
            if a == 0:
                out += (sp.diff(XR, t) + w*XI) if k == 0 else (sp.diff(XI, t) - w*XR)
            else:
                out += sp.diff(XR if k == 0 else XI, X[a])
        return out/sqrtg - dV2(sq)*(U if k == 0 else W)
    eqs = {}
    for (a, b) in [(0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2), (3, 3), (0, 3), (1, 3), (2, 3)]:
        eqs[(a, b)] = dR[a, b] - 8*sp.pi*eps*dS[a, b]
    for k in (0, 1):
        eqs['KG%d' % k] = sp.expand(sp.diff(KG(k), e).subs(e, 0))
    bgeq = {'R': R0 - 8*sp.pi*eps*S0, 'KG': sp.simplify(KG(0).subs(e, 0))}
    return eqs, bgeq

# ---------- background equations (solve for lam0', nu0', phi0'') ----------
def background_rules():
    a2 = sp.exp(lam0); al2 = sp.exp(nu0)
    rho = (w**2*p0**2/al2 + sp.diff(p0, r)**2/a2)/2 + Vf(p0**2)
    pr = (w**2*p0**2/al2 + sp.diff(p0, r)**2/a2)/2 - Vf(p0**2)
    dnu = (a2 - 1)/r + 8*sp.pi*r*a2*eps*pr
    dlam = (1 - a2)/r + 8*sp.pi*r*a2*eps*rho
    ddp = -(2/r + dnu/2 - dlam/2)*sp.diff(p0, r) + a2*(dV2(p0**2) - w**2/al2)*p0
    return dnu, dlam, ddp
dnu0, dlam0, ddp0 = background_rules()

def bgsub(z):
    for _ in range(3):
        z = z.subs({sp.diff(nu0, r, 3): sp.diff(dnu0, r, 2), sp.diff(lam0, r, 3): sp.diff(dlam0, r, 2)})
        z = z.subs({sp.diff(nu0, r, 2): sp.diff(dnu0, r), sp.diff(lam0, r, 2): sp.diff(dlam0, r)})
        z = z.subs(sp.diff(p0, r, 3), sp.diff(ddp0, r)).subs(sp.diff(p0, r, 2), ddp0)
        z = z.subs({sp.diff(nu0, r): dnu0, sp.diff(lam0, r): dlam0})
    return z

def angsub(z):
    """Y'' -> -cot Y' - LL Y (Legendre), then Y, Y' -> symbols."""
    Ypp = -sp.cos(th)/sp.sin(th)*sp.diff(Yf, th) - LL*Yf
    z = z.subs(sp.diff(Yf, th, 3), sp.diff(Ypp, th)).subs(sp.diff(Yf, th, 2), Ypp)
    z = z.subs(sp.diff(Yf, th, 3), sp.diff(Ypp, th)).subs(sp.diff(Yf, th, 2), Ypp)
    z = z.subs(sp.diff(Yf, th), Ypsym).subs(Yf, Ysym)
    return sp.expand(z)

# =====================================================================================
# POLAR
# =====================================================================================
H0, H1, H2, K, F, G = [sp.Function(n)(r) for n in ('H0', 'H1', 'H2', 'K', 'F', 'G')]
h = sp.zeros(4, 4)
h[0, 0] = sp.exp(nu0)*H0*Yf*ex
h[0, 1] = h[1, 0] = H1*Yf*ex
h[1, 1] = sp.exp(lam0)*H2*Yf*ex
h[2, 2] = r**2*K*Yf*ex
h[3, 3] = r**2*sp.sin(th)**2*K*Yf*ex
print("linearizing polar ...", flush=True)
eqs, bg = linearize(h, F*Yf*ex, G*Yf*ex)
print("  done %.0fs" % (time.time() - T0), flush=True)

# background check
bgR = bg['R'].applyfunc(lambda z: sp.simplify(bgsub(z)))
print("background Einstein residual (all zero?):", all(z == 0 for z in bgR), " KG:", sp.simplify(bgsub(bg['KG'])))

def clean(z):
    z = sp.expand(sp.simplify(z/ex))
    return angsub(bgsub(z))

P = {}
for k, v in eqs.items():
    if k in [(0, 3), (1, 3), (2, 3)]:
        continue
    P[k] = clean(v)
    print("  cleaned", k, "%.0fs" % (time.time() - T0), flush=True)

def coeff(z, which):
    z = sp.expand(z)
    return sp.simplify(z.coeff(which))

# separate angular structure
E = {}
E['tt'] = coeff(P[(0, 0)], Ysym); E['tr'] = coeff(P[(0, 1)], Ysym); E['rr'] = coeff(P[(1, 1)], Ysym)
E['tth'] = coeff(P[(0, 2)], Ypsym); E['rth'] = coeff(P[(1, 2)], Ypsym)
E['KGR'] = coeff(P['KG0'], Ysym); E['KGI'] = coeff(P['KG1'], Ysym)
# angular block: R_thth and R_phph/sin^2 : a Y + b cot Y'
ct = sp.cos(th)/sp.sin(th)
Rthth = P[(2, 2)]; Rphph = sp.expand(sp.simplify(P[(3, 3)]/sp.sin(th)**2))
# traceless part: coefficient of cot*Y' in (thth) (it equals minus that in phph)
Rthth_s = sp.expand(Rthth.subs(sp.cos(th), ct*sp.sin(th)))
tl = sp.simplify(sp.expand(Rthth - Rphph).coeff(Ypsym))
E['ang_tl'] = sp.simplify(tl/ct) if tl != 0 else 0
E['ang_tr'] = sp.simplify(sp.expand(Rthth + Rphph).coeff(Ysym))
chk_angle = sp.simplify(sp.expand(Rthth + Rphph).coeff(Ypsym))
print("trace part free of Y'? ", chk_angle == 0)
for k in E:
    print(f"  E[{k}] has theta?", E[k].has(th) if E[k] != 0 else False)

print("\ntraceless angular eq:", sp.factor(E['ang_tl']))
# expect proportional to (H0 - H2)
sol_H2 = sp.solve(E['ang_tl'], H2)
print("  => H2 =", sol_H2)
H2sol = sol_H2[0]

def subH2(z):
    return sp.expand(z.subs(sp.diff(H2, r, 2), sp.diff(H2sol, r, 2)).subs(sp.diff(H2, r), sp.diff(H2sol, r)).subs(H2, H2sol).doit())

for k in list(E):
    E[k] = subH2(E[k])

# first-order equations
dH1 = sp.diff(H1, r); dK = sp.diff(K, r); dH0 = sp.diff(H0, r)
print("\nderivative content:")
for k in ('tth', 'tr', 'rth', 'rr', 'tt', 'ang_tr', 'KGR', 'KGI'):
    z = E[k]
    print(f"  {k}: H1'={z.has(dH1)} K'={z.has(dK)} H0'={z.has(dH0)} K''={z.has(sp.diff(K,r,2))} "
          f"H0''={z.has(sp.diff(H0,r,2))} H1''={z.has(sp.diff(H1,r,2))} F''={z.has(sp.diff(F,r,2))} G''={z.has(sp.diff(G,r,2))}")

sol1 = sp.solve([E['tth'], E['tr'], E['rth']], [dH1, dK, dH0], dict=True)
assert len(sol1) == 1
sol1 = {k: sp.simplify(v) for k, v in sol1[0].items()}
print("  solved (t th),(t r),(r th) for H1', K', H0'  %.0fs" % (time.time() - T0), flush=True)

def sub1(z):
    for _ in range(2):
        z = z.subs(sp.diff(H1, r, 2), sp.diff(sol1[dH1], r)).subs(sp.diff(K, r, 2), sp.diff(sol1[dK], r))
        z = z.subs(sp.diff(H0, r, 2), sp.diff(sol1[dH0], r))
        z = z.subs({dH1: sol1[dH1], dK: sol1[dK], dH0: sol1[dH0]})
        z = bgsub(z)
    return sp.expand(z)

rr_red = sp.simplify(sub1(E['rr']))
print("\n(rr) after substitution: contains H0?", rr_red.has(H0), "  derivatives?",
      any(rr_red.has(sp.diff(q, r)) for q in (H0, H1, K)))
# WARNING: rr_red still contains G'' (sub1 does not eliminate it), so this H0 / "red" is NOT the correct
# reduction; it is kept only for the pickle. Use reduce_polar.py.
H0sol = sp.solve(rr_red, H0)
assert len(H0sol) == 1
H0sol = sp.simplify(H0sol[0])
print("  algebraic H0 obtained  %.0fs" % (time.time() - T0), flush=True)

# KG equations: solve for F'', G''
KGsol = sp.solve([E['KGR'], E['KGI']], [sp.diff(F, r, 2), sp.diff(G, r, 2)], dict=True)[0]
KGsol = {k: sp.simplify(sub1(v)) for k, v in KGsol.items()}

# reduced first-order system:  H1', K' (with H0 -> H0sol), F'', G''
def subH0(z):
    return z.subs(sp.diff(H0, r), sp.diff(H0sol, r)).subs(H0, H0sol)

red = {'dH1': sp.simplify(subH0(sol1[dH1])), 'dK': sp.simplify(subH0(sol1[dK])),
       'ddF': sp.simplify(subH0(KGsol[sp.diff(F, r, 2)])), 'ddG': sp.simplify(subH0(KGsol[sp.diff(G, r, 2)])),
       'H0': H0sol}
print("  reduced system built  %.0fs" % (time.time() - T0), flush=True)

with open(os.path.join(HERE, 'polar_eqs.pkl'), 'wb') as fh:
    pickle.dump({'E': E, 'red': red, 'H2sol': H2sol, 'sol1': sol1}, fh)

# Reduction to the final first-order system, code generation and verification are done in
#   reduce_polar.py  (joint solve of (rr), KG_R, KG_I for H0, F'', G''; vacuum Zerilli/KG checks)
#   verify_polar.py  (all 8 separated equations satisfied on solutions of the generated system)
# NOTE: an earlier version of this script solved H0 from (rr) alone (it still contains G'') and generated
# code from that; that output was wrong and has been deleted.
print("saved polar_eqs.pkl  %.0fs" % (time.time() - T0))

"""Shared symbolic linearization machinery (extracted verbatim from derive_nonradial.py)."""
import sympy as sp
import time
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


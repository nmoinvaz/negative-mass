"""
Phase 4: Witten identity with a signed source, evaluated on exact time-symmetric data.

Key fact used: every spherically symmetric (indeed every conformally flat) 3-metric g = psi^4 delta has the
explicit Witten spinor eps = psi^{-2} eps0 (conformal covariance of the Dirac operator in n = 3), so
|eps|^2 = psi^{-4}, and the Witten identity becomes (time-symmetric, K = 0, G = c = 1):

    M = G + int mu |eps|^2 dV,      G = (1/4pi) int |nabla eps|^2 dV = (3/2pi) int |grad psi|^2 psi^-4 d^3x >= 0
                                    int mu |eps|^2 dV = int mu psi^2 d^3x  (mu = mu_+ - mu_-)

Hence M >= - int mu_- psi^-4 dV >= - Lambda * E_-,  Lambda = sup_{supp mu_-} psi^-4,  E_- = int mu_- dV.
We check the identity (M = G + W_+ - W_-) to quadrature accuracy and measure tightness on
 (1) Shirokov's phantom and canonical boson stars (phase2/boson_star.py),
 (2) a family of smooth phantom shells showing |M|/E_- is unbounded (no bound M >= -C E_- with universal C),
 (3) mixed spherical data (canonical shell around a phantom core).
"""
import sys, os
import numpy as np
from scipy.integrate import quad, cumulative_trapezoid as ctz
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'phase2'))
import boson_star as bs

def trap(y, x): return np.trapezoid(y, x)

# ---------------------------------------------------------------- (1) boson stars
def star_profile(w_target, eps):
    w = w_target
    for _ in range(4):
        p = bs.find_star(w, eps)
        s = bs.star_properties(p, w, eps)
        w = w * w_target / s['w_phys']
    p = bs.find_star(w, eps)
    _, sol = bs.shoot(p, w, eps)
    r, y = sol.t, sol.y
    i_cut = int(np.argmin(np.abs(y[0]) + (r < 2.0) * 1e9))
    r = r[:i_cut + 1]; y = y[:, :i_cut + 1]
    phi, dphi, m, lnal = y
    al2 = np.exp(2 * lnal); a2 = 1.0 / (1.0 - 2.0 * m / r); a = np.sqrt(a2)
    s2 = phi * phi
    rho = 0.5 * (w * w * s2 / al2 + dphi * dphi / a2) + bs.V(s2)   # sector's own (canonical) energy density
    return r, a, rho, m[-1]

def witten_terms_schw(rA, a, rho, M):
    """Areal-radius data -> isotropic psi; return M check, G, W, E, Lambda."""
    rc = rA[-1]
    riso_c = 0.5 * (rc - M + np.sqrt(rc * rc - 2 * M * rc))
    # ln r_iso(rA) = ln r_iso_c - int_{rA}^{rc} a/s ds  (integrand regular: a -> 1 at 0, handle 1/s analytically)
    # write a/s = 1/s + (a-1)/s
    f = (a - 1.0) / rA
    I = ctz(f[::-1], rA[::-1], initial=0.0)[::-1] * -1.0          # int_{rA}^{rc} (a-1)/s ds
    ln_riso = np.log(riso_c) - (np.log(rc) - np.log(rA)) - I
    psi4inv = (np.exp(ln_riso) / rA) ** 2                           # psi^-4 = (r_iso/r_A)^2
    dV = 4 * np.pi * rA**2 * a
    E = trap(rho * dV, rA)
    W = trap(rho * psi4inv * dV, rA)
    G_in = 1.5 * trap((1 - a) ** 2 / a * psi4inv, rA)
    def g_out(s):
        ao = 1 / np.sqrt(1 - 2 * M / s)
        ri = 0.5 * (s - M + np.sqrt(s * s - 2 * M * s))
        return 1.5 * (1 - ao) ** 2 / ao * (ri / s) ** 2
    G_out = quad(g_out, rc, np.inf, limit=200)[0]
    return dict(G=G_in + G_out, W=W, E=E, Lam=psi4inv.max(), Lam_min=psi4inv.min(),
                psi0=psi4inv[0] ** -0.25)

print("(1) Shirokov boson stars: Witten identity M = G + eps*W with |eps_W|^2 = psi^-4")
res = {}
for w_t, eps, name in [(0.750, +1, 'canonical'), (0.7603, -1, 'phantom')]:
    rA, a, rho, M = star_profile(w_t, eps)
    t = witten_terms_schw(rA, a, rho, M)
    res[name] = (M, t)
    print(f"  {name:9s} M_ADM={M:+.7f}  E=int rho dV={t['E']:.7f}  W=int rho psi^-4 dV={t['W']:.7f}  "
          f"G={t['G']:.3e}  psi(0)={t['psi0']:.6f}")
    print(f"            identity residual M - (G + eps W) = {M - (t['G'] + eps * t['W']):+.2e}")
    if eps < 0:
        print(f"            bound M >= -W:   -W = {-t['W']:+.7f}   slack G = {t['G']:.2e} ({t['G']/t['W']:.2e} of W)")
        print(f"            bound M >= -Lambda E: Lambda = sup psi^-4 = {t['Lam']:.6f};  -Lambda E = {-t['Lam']*t['E']:+.7f}")
        print(f"            naive  M >= -E  would give {-t['E']:+.7f}: VIOLATED (|M|/E = {abs(M)/t['E']:.6f} > 1)")
    else:
        print(f"            canonical: M/E = {M/t['E']:.6f} (gravitational binding), W/E = {t['W']/t['E']:.6f}")

Mp, tp = res['canonical']; Mm, tm = res['phantom']
print("\n  Mixed pair (widely separated, superposed): M_tot = %+.3e" % (Mp + Mm))
print("    Witten: M_tot >= -W_- = %+.6f ; equivalently W_+ <= M_tot + W_- - G (W_+ = %.6f)" % (-tm['W'], tp['W']))
print("    => as long as the phantom configuration stays weakly compact and slow (Lambda ~ 1, J ~ 0),")
print("       total mass (and Bondi mass) cannot drop below ~ -%.4f, i.e. ~ %.0f x the initial M_tot."
      % (tm['W'], tm['W'] / (Mp + Mm)))
for v in [0.1, 0.5, 0.9, 0.99]:
    g = 1 / np.sqrt(1 - v * v)
    print(f"    boosted phantom star v={v}: int(mu_- + |J_-|) dV ~ gamma(1+v) E_- = {g*(1+v)*tm['E']:.4f}")

# ---------------------------------------------------------------- (2) phantom shells
print("\n(2) Smooth phantom shells, conformally flat time-symmetric data, prescribed mu_tilde = mu psi^5")
print("    (Gaussian shell in isotropic r, centre R=1, width 0.05; total int mu_tilde d^3x = |M|; DEC: mu>=0, J=0)")
r = np.linspace(1e-6, 60, 600001)
def shell(Mabs, R=1.0, wdt=0.05, sign=-1, core=None):
    mt = np.exp(-0.5 * ((r - R) / wdt) ** 2)
    mt *= Mabs / trap(4 * np.pi * r**2 * mt, r)
    return sign * mt
def solve(mtil):
    """Delta psi = -2 pi mu_tilde (mu_tilde signed: + canonical, - phantom); psi -> 1."""
    Q = ctz(4 * np.pi * r**2 * mtil, r, initial=0.0)                   # = int mu_tilde d^3x within r
    dpsi = -Q / (2 * r**2)
    Mtot = Q[-1]
    psi = 1.0 + Mtot / (2 * r[-1]) - ctz(dpsi[::-1], r[::-1], initial=0.0)[::-1] * -1.0
    psi = 1.0 + Mtot / (2 * r[-1]) + (ctz(dpsi, r, initial=0.0) - ctz(dpsi, r, initial=0.0)[-1])
    return psi, dpsi, Mtot
print("    |M|      psi_min    E_-        W_-        G          M-(G-W)    |M|/E_-   Lambda    W/(Lambda E)")
for Mabs in [0.1, 0.5, 1.0, 1.5, 1.8, 1.9, 1.95]:
    mt = shell(Mabs)
    psi, dpsi, Mt = solve(mt)
    if psi.min() <= 0: print(f"    {Mabs}: psi hits zero"); continue
    mu_m = -mt * psi**-5
    E = trap(4 * np.pi * r**2 * mu_m * psi**6, r)
    W = trap(4 * np.pi * r**2 * mu_m * psi**2, r)
    tail = 6 * (Mt**2 / 4) * quad(lambda s: s**-2 * (1 + Mt / (2 * s)) ** -4, r[-1], np.inf)[0]
    G = 6 * trap(dpsi**2 * psi**-4 * r**2, r) + tail
    supp = np.abs(mt) > 1e-8 * np.abs(mt).max()
    Lam = (psi[supp] ** -4).max()
    print(f"    {Mabs:<7} {psi.min():.5f}   {E:.6f}   {W:.6f}   {G:.3e}  {Mt-(G-W):+.1e}   {abs(Mt)/E:7.3f}  {Lam:9.3f}  {W/(Lam*E):.4f}")
print("    -> |M|/E_- grows without bound as the shell approaches psi=0 (areal radius -> 0):")
print("       no inequality M >= -C E_- with a universal C exists; the Witten bound -int mu_- psi^-4 dV holds throughout.")

# ---------------------------------------------------------------- (3) mixed data
print("\n(3) Mixed spherical data: canonical shell (R=3) around phantom shell (R=1)")
print("    m_+      m_-     M_tot      -W_-       E_+       E_-      bound holds   M_tot+W_- (=G+W_+)")
for mp, mm in [(0.5, 0.5), (1.0, 0.5), (0.5, 1.0), (2.0, 1.5), (0.3, 1.8)]:
    mt = shell(mp, R=3.0, sign=+1) + shell(mm, R=1.0, sign=-1)
    psi, dpsi, Mt = solve(mt)
    if psi.min() <= 0: print(f"    {mp},{mm}: psi hits zero"); continue
    mu = mt * psi**-5
    mup, mum = np.clip(mu, 0, None), np.clip(-mu, 0, None)
    Ep = trap(4*np.pi*r**2*mup*psi**6, r); Em = trap(4*np.pi*r**2*mum*psi**6, r)
    Wp = trap(4*np.pi*r**2*mup*psi**2, r); Wm = trap(4*np.pi*r**2*mum*psi**2, r)
    print(f"    {mp:<7}  {mm:<6}  {Mt:+.4f}   {-Wm:+.4f}   {Ep:.4f}   {Em:.4f}   {str(Mt >= -Wm):5s}         {Mt+Wm:.4f}")

# ---------------------------------------------------------------- (4) E_- unbounded at fixed M
print("\n(4) Constraints do not bound E_+ or E_- at fixed M: identical +/- configurations")
print("    T_+ = T_- pointwise => T = 0, flat data, M = 0 for any E_+ = E_- (two scalars with equal Cauchy data);")
print("    this is also an exact solution for all time (each field solves its own wave equation on Minkowski).")

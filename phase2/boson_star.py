"""
Static sextic boson stars with a signed gravitational source (Shirokov, arXiv:2608.24577, model).

Matter (both sectors): complex scalar, L = -1/2 |dPhi|^2 - V(|Phi|^2),
    V(s) = 1/2 m^2 s - 1/4 lam s^2 + 1/6 mu s^3,       Box Phi = 2 dV/ds Phi.
Gravity: G_ab = 8 pi eps T_ab,  eps = +1 canonical star, eps = -1 phantom star (negative active / ADM mass,
positive inertial and passive mass: the field equation is unchanged).

Metric (Schwarzschild-like): ds^2 = -alpha^2 dt^2 + a^2 dr^2 + r^2 dOmega^2, a^-2 = 1 - 2 m(r)/r.
Ansatz Phi = phi(r) e^{-i w t}. With rho = 1/2 (w^2 phi^2/alpha^2 + phi'^2/a^2) + V,
p_r = 1/2 (w^2 phi^2/alpha^2 + phi'^2/a^2) - V:
    m'          = 4 pi r^2 eps rho
    alpha'/alpha = (a^2 - 1)/(2r) + 4 pi r a^2 eps p_r
    phi''        = -(2/r + alpha'/alpha - a'/a) phi' + a^2 (2 dV/ds - w^2/alpha^2) phi
Units m = 1, G = c = 1. Couplings as in Shirokov: lam = 10240, lam^2/mu = 4.8.
The physical frequency is w_phys = w / alpha(inf) (we set alpha(0) = 1 and rescale).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

LAM = 10240.0
MU = LAM**2 / 4.8

def dVds2(s):            # 2 dV/ds
    return 1.0 - LAM * s + MU * s * s

def V(s):
    return 0.5 * s - 0.25 * LAM * s * s + MU * s**3 / 6.0

def rhs(r, y, w, eps):
    phi, dphi, m, lnal = y
    al2 = np.exp(2 * lnal)
    a2 = 1.0 / (1.0 - 2.0 * m / r)
    s = phi * phi
    kin = 0.5 * (w * w * s / al2 + dphi * dphi / a2)
    rho = kin + V(s)
    pr = kin - V(s)
    dm = 4 * np.pi * r * r * eps * rho
    dlnal = (a2 - 1.0) / (2 * r) + 4 * np.pi * r * a2 * eps * pr
    dlna = (1.0 - a2) / (2 * r) + 4 * np.pi * r * a2 * eps * rho
    ddphi = -(2.0 / r + dlnal - dlna) * dphi + a2 * (dVds2(s) - w * w / al2) * phi
    return [dphi, ddphi, dm, dlnal]

def shoot(phic, w, eps, rmax=80.0):
    """Integrate outward; return +1 if phi turns back up (undershoot), -1 if it crosses zero (overshoot)."""
    r0 = 1e-6
    # series start: phi'' at origin = (1/3) (2V'(s) - w^2) phi
    c2 = (dVds2(phic * phic) - w * w) * phic / 6.0
    y0 = [phic + c2 * r0**2, 2 * c2 * r0, 0.0, 0.0]
    def cross(r, y, *a): return y[0]
    cross.terminal = True; cross.direction = -1
    # genuine undershoot: phi has fallen below phi_c and turns back up
    def turn(r, y, *a): return y[1] if y[0] < 0.99 * phic else -1.0
    turn.terminal = True; turn.direction = 1
    # runaway past the hill top: counts as "too large", same side as overshoot
    def runaway(r, y, *a): return y[0] - 1.5 * phic
    runaway.terminal = True; runaway.direction = 1
    sol = solve_ivp(rhs, (r0, rmax), y0, args=(w, eps), events=(cross, turn, runaway),
                    rtol=1e-12, atol=1e-16, method='DOP853', max_step=0.05)
    if sol.t_events[0].size or sol.t_events[2].size or not sol.success: return -1, sol
    if sol.t_events[1].size: return +1, sol
    return 0, sol

def find_star(w, eps, phic_lo=0.0, phic_hi=0.0, tol=1e-13):
    """Bisect the central amplitude at fixed (alpha(0)=1) frequency w."""
    have_bracket = phic_lo > 0 and shoot(phic_lo, w, eps)[0] == +1 and shoot(phic_hi, w, eps)[0] == -1
    if not have_bracket:
        # Scan for the physical bracket: undershoot (+1) just below overshoot (-1) as phi_c increases.
        # (Very small phi_c also reads -1: the linear tail grows and runs away; skip that branch.)
        grid = np.linspace(0.008, 0.03, 23)
        signs = [shoot(p, w, eps, rmax=30.0)[0] for p in grid]
        for i in range(len(grid) - 2, -1, -1):
            if signs[i] == +1 and signs[i + 1] == -1:
                phic_lo, phic_hi = grid[i], grid[i + 1]
                break
        else:
            raise RuntimeError(f"no bracket at w={w}, eps={eps}: {signs}")
    flo = +1
    lo, hi = phic_lo, phic_hi
    while hi - lo > tol * max(1.0, hi):
        mid = 0.5 * (lo + hi)
        f, _ = shoot(mid, w, eps)
        if f == flo: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)

def star_properties(phic, w, eps):
    """Integrate the bracketed solution as far as it stays bound, then attach the vacuum exterior."""
    _, sol = shoot(phic, w, eps)
    r, y = sol.t, sol.y
    phi = y[0]
    # cut where |phi| reaches its minimum before the solution peels away
    i_cut = int(np.argmin(np.abs(phi) + (r < 2.0) * 1e9))
    m_cut = y[2, i_cut]
    rc = r[i_cut]
    # exterior: alpha^2 -> (1 - 2M/r) up to the constant;   alpha(inf) from matching at rc
    lnal_c = y[3, i_cut]
    al_inf = np.exp(lnal_c) / np.sqrt(1.0 - 2.0 * m_cut / rc)
    w_phys = w / al_inf
    # rms radius weighted by |phi|^2 (proper volume)
    a = 1.0 / np.sqrt(1.0 - 2.0 * y[2, :i_cut] / r[:i_cut])
    wgt = phi[:i_cut]**2 * a * r[:i_cut]**2
    rrms = np.sqrt(np.trapezoid(wgt * r[:i_cut]**2, r[:i_cut]) / np.trapezoid(wgt, r[:i_cut]))
    # Shirokov's tracker weight is (|Phi|^2 + |Pi|^2)^{1/2}, i.e. proportional to |phi|
    wgt1 = np.abs(phi[:i_cut]) * a * r[:i_cut]**2
    rrms1 = np.sqrt(np.trapezoid(wgt1 * r[:i_cut]**2, r[:i_cut]) / np.trapezoid(wgt1, r[:i_cut]))
    return dict(M=m_cut, w_phys=w_phys, rrms=rrms, rrms_abs=rrms1, rcut=rc, phic=phic,
                phi_cut=phi[i_cut], compact=abs(m_cut) / rrms1)

def star_at_physical_frequency(w_target, eps, iters=4):
    """Fixed-point iteration w <- w_target * alpha(inf); converges fast in weak field."""
    w = w_target
    for _ in range(iters):
        p = find_star(w, eps)
        s = star_properties(p, w, eps)
        w = w * w_target / s['w_phys']
    return s

if __name__ == "__main__":
    print("Shirokov campaign stars (paper: canonical w=0.750 -> M=+0.014350; phantom w=0.7603 -> M=-0.014295)")
    for w_t, eps in [(0.750, +1), (0.7603, -1)]:
        s = star_at_physical_frequency(w_t, eps)
        print(f"  eps={eps:+d} w_phys={s['w_phys']:.6f}  M_ADM={s['M']:+.6f}  phi_c={s['phic']:.6f}  "
              f"r_rms(phi^2)={s['rrms']:.3f} r_rms(|phi|)={s['rrms_abs']:.3f}  |M|/r_rms={s['compact']:.2e}  (cut r={s['rcut']:.1f}, phi={s['phi_cut']:.1e})")
    print()
    print("M(w) families (alpha(0)=1 frequency scanned, physical frequency reported)")
    for eps in (+1, -1):
        print(f"  eps={eps:+d}")
        for w in [0.70, 0.75, 0.80, 0.85, 0.90]:
            try:
                p = find_star(w, eps)
                s = star_properties(p, w, eps)
                print(f"    w_phys={s['w_phys']:.4f}  M={s['M']:+.6f}  phi_c={p:.5f}  r_rms(|phi|)={s['rrms_abs']:.3f}")
            except RuntimeError as e:
                print(f"    w={w}: {e}")

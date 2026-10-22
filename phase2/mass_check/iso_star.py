"""Independent re-implementation (written from the equations in the paper repo's
boson_star_ode.py docstring/MATTER_MODEL.md; none of their code is executed) of
Shirokov's isotropic, maximally sliced dressed star, plus BOTH mass extractions:

  paper:  read at r_m = first point with phi < 5% phi_c (still descending):
          psi_inf = psi + r psi',  alpha_inf = alpha + r alpha'   (A + B/r fit at r_m)
          omega_phys = omega_int/alpha_inf, r_phys = psi_inf^2 r
          M_paper = 2 r_phys (psi(r_m)/psi_inf - 1)
  true:   same A+B/r reads but far outside the matter (phi ~ 1e-8 phi_c), where
          psi, alpha are exactly harmonic:  M_true = 2 A B,  B = -r^2 psi'.
"""
import sys, math, numpy as np
from scipy.integrate import solve_ivp

LAM = 10240.0
MU_EXACT = LAM**2 / 4.8           # 21845333.33
MU_PAPER = 21845333.0             # literal in their scripts & table headers

def make_rhs(w, eps, lam, mu):
    def rhs(r, y):
        phi, dphi, psi, dpsi, al, dal = y
        p4 = psi**4; s = phi*phi
        U = 0.5*s - 0.25*lam*s*s + mu*s**3/6.0
        dU = phi - lam*s*phi + mu*s*s*phi
        oa2 = (w/al)**2
        rho = eps*(0.5*oa2*s + 0.5*dphi*dphi/p4 + U)
        rps = eps*(2*oa2*s - 2*U)
        d2psi = -2/r*dpsi - 2*math.pi*p4*psi*rho
        d2al = -(2/r + 2*dpsi/psi)*dal + 4*math.pi*p4*al*rps
        d2phi = -(2/r + dal/al + 2*dpsi/psi)*dphi + p4*(dU - oa2*phi)
        return [dphi, d2phi, dpsi, d2psi, dal, d2al]
    return rhs

def shoot(phic, w, eps, lam, mu, rmax=60.0):
    rhs = make_rhs(w, eps, lam, mu)
    r0 = 1e-5
    s = phic*phic
    dU = phic - lam*s*phic + mu*s*s*phic
    c = (dU - w*w*phic)/6.0
    y0 = [phic + c*r0*r0, 2*c*r0, 1.0, 0.0, 1.0, 0.0]
    def cross(r, y): return y[0]
    cross.terminal = True; cross.direction = -1
    def turn(r, y): return y[1] if y[0] < 0.99*phic else -1.0
    turn.terminal = True; turn.direction = 1
    def run(r, y): return y[0] - 1.5*phic
    run.terminal = True; run.direction = 1
    sol = solve_ivp(rhs, (r0, rmax), y0, method='DOP853', rtol=1e-12, atol=1e-16,
                    max_step=0.05, events=(cross, turn, run), dense_output=True)
    if sol.t_events[0].size or sol.t_events[2].size or not sol.success: return -1, sol
    if sol.t_events[1].size: return +1, sol
    return 0, sol

def find(w, eps, lam, mu):
    grid = np.linspace(0.008, 0.03, 45)
    sg = [shoot(p, w, eps, lam, mu, rmax=30)[0] for p in grid]
    for i in range(len(grid)-2, -1, -1):
        if sg[i] == +1 and sg[i+1] == -1: lo, hi = grid[i], grid[i+1]; break
    else: raise RuntimeError(sg)
    while hi - lo > 1e-14:
        mid = 0.5*(lo+hi)
        if shoot(mid, w, eps, lam, mu)[0] == +1: lo = mid
        else: hi = mid
    return lo

def extract(phic, w, eps, lam, mu):
    _, sol = shoot(phic, w, eps, lam, mu)
    r = sol.t; phi, dphi, psi, dpsi, al, dal = sol.y
    # --- paper: first point with phi < 5% phi_c (descending)
    i = np.where((phi < 0.05*phic) & (dphi < 0) & (r > 1e-3))[0][0]
    rm = r[i-1] + (r[i]-r[i-1]) * (phi[i-1]-0.05*phic)/(phi[i-1]-phi[i])
    Y = sol.sol(rm)
    psi_inf = Y[2] + rm*Y[3]; al_inf = Y[4] + rm*Y[5]
    M_paper = 2*psi_inf**2*rm*(Y[2]/psi_inf - 1)
    w_paper = w/al_inf
    # --- true: far point where |phi| is smallest before the separatrix peel-off
    j = int(np.argmin(np.abs(phi) + (r < 2)*1e9))
    A = psi[j] + r[j]*dpsi[j]; B = -r[j]**2*dpsi[j]
    Aal = al[j] + r[j]*dal[j]
    M_true = 2*A*B
    w_true = w/Aal
    return dict(phic=phic, w_int=w, rm_iso_phys=psi_inf**2*rm, R_areal_m=Y[2]**2*rm,
                M_paper=M_paper, w_paper=w_paper, M_true=M_true, w_true=w_true,
                r_far=r[j], phi_far=phi[j], sol=sol, psi_inf=psi_inf, al_inf=al_inf, A=A, Aal=Aal)

def at_omega(w_target, eps, lam=LAM, mu=MU_PAPER, which='w_paper', iters=5):
    """Outer-iterate omega_int so that the chosen physical-frequency definition hits w_target."""
    w = w_target
    for _ in range(iters):
        p = find(w, eps, lam, mu)
        e = extract(p, w, eps, lam, mu)
        if abs(e[which] - w_target) < 1e-9: break
        w *= w_target / e[which]
    return e

if __name__ == "__main__":
    for w_t, eps, Mp in [(0.750, +1, 0.014350), (0.7603, -1, -0.014295)]:
        for which in ('w_paper', 'w_true'):
            e = at_omega(w_t, eps, which=which)
            print(f"eps={eps:+d} target {w_t} as {which:7s}: phic={e['phic']:.9f} "
                  f"w_paper={e['w_paper']:.7f} w_true={e['w_true']:.7f}  "
                  f"M_paper-style={e['M_paper']:+.7f}  M_true={e['M_true']:+.7f}  "
                  f"(paper {Mp:+.6f}; r_m iso={e['rm_iso_phys']:.3f}, areal={e['R_areal_m']:.3f}; far r={e['r_far']:.1f} phi={e['phi_far']:.1e})", flush=True)

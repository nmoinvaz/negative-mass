"""COPY of ../radial/background.py (frozen 2026-10-21 for reproducibility; original by the l=0 agent).

Static sextic boson stars (same model as ../boson_star.py) re-solved as a two-point BVP so that the
background is accurate all the way to large r (needed for the perturbation problem).

Variables  y = (psi, u, m, nu)  with  phi = e^psi,  u = phi'/phi,  a^-2 = 1 - 2m/r,  alpha^2 = e^nu.
Normalization alpha(inf) = 1, so omega IS the physical frequency.
    psi' = u
    u'   = -u^2 - (2/r + nu'/2 - lam'/2) u + a^2 (2V'(phi^2) - w^2/alpha^2)
    m'   = 4 pi r^2 eps rho
    nu'  = (a^2-1)/r + 8 pi r a^2 eps p_r
BCs: u(r0) = r0 * c/3 (regular series), m(r0) = 4pi/3 eps rho_c r0^3,
     u(R) = -k - 1/R (decaying Yukawa tail, k = sqrt(1 - w^2)), nu(R) = ln(1 - 2m(R)/R).
Using psi = ln phi removes the trivial solution phi = 0 and keeps the exponentially small tail
well-conditioned (u is O(1) everywhere). The tail BC error is O(1/R^2) at R but is damped like
e^{-2k(R-r)} inward (the decaying solution is the attractor of the inward Riccati flow).
"""
import sys, os
import numpy as np
from scipy.integrate import solve_bvp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import boson_star as bs

LAM, MU = bs.LAM, bs.MU
dV2, V = bs.dVds2, bs.V

def odes(r, y, w, eps):
    psi, u, m, nu = y
    phi2 = np.exp(2 * psi)
    a2 = 1.0 / (1.0 - 2.0 * m / r)
    al2 = np.exp(nu)
    dphi2 = (u * u) * phi2
    kin = 0.5 * (w * w * phi2 / al2 + dphi2 / a2)
    rho = kin + V(phi2); pr = kin - V(phi2)
    dm = 4 * np.pi * r * r * eps * rho
    dnu = (a2 - 1.0) / r + 8 * np.pi * r * a2 * eps * pr
    dlam = (1.0 - a2) / r + 8 * np.pi * r * a2 * eps * rho
    du = -u * u - (2.0 / r + 0.5 * dnu - 0.5 * dlam) * u + a2 * (dV2(phi2) - w * w / al2)
    return np.vstack([u, du, dm, dnu])

def bcs(ya, yb, r0, R, w, eps):
    psi0, u0, m0, nu0 = ya
    phi2 = np.exp(2 * psi0); al2 = np.exp(nu0)
    c = dV2(phi2) - w * w / al2
    rho_c = 0.5 * w * w * phi2 / al2 + V(phi2)
    k = np.sqrt(1.0 - w * w)
    return np.array([u0 - r0 * c / 3.0,
                     m0 - 4 * np.pi / 3 * eps * rho_c * r0**3,
                     yb[1] + k + 1.0 / R,
                     yb[3] - np.log(1.0 - 2.0 * yb[2] / R)])

class Star:
    """Holds a converged background and evaluates (alpha, a, phi0, u) anywhere on [r0, R]."""
    def __init__(self, sol, w, eps, r0, R):
        self.sol, self.w, self.eps, self.r0, self.R = sol, w, eps, r0, R
        rr = np.linspace(r0, R, 20001)
        y = sol.sol(rr)
        self.M = float(sol.y[2, -1])
        a = 1.0 / np.sqrt(1.0 - 2 * y[2] / rr); al = np.exp(y[3] / 2)
        phi2 = np.exp(2 * y[0])
        self.N = float(np.trapezoid(4 * np.pi * rr**2 * w * phi2 * a / al, rr))
        self.phic = float(np.exp(sol.y[0, 0]))
        self.alpha_c = float(np.exp(sol.y[3, 0] / 2))
        self.max_res = float(np.max(sol.rms_residuals))
        self.nmesh = sol.x.size
    def bg(self, r):
        y = self.sol.sol(r)
        a = 1.0 / np.sqrt(1.0 - 2 * y[2] / r)
        al = np.exp(y[3] / 2)
        phi0 = np.exp(y[0])
        return al, a, phi0, y[1]

def solve_star(w, eps, guess=None, r0=1e-3, R=None, tol=1e-9, max_nodes=200000):
    """guess: a previous Star (continuation) or None (use boson_star.py shooting)."""
    k = np.sqrt(1.0 - w * w)
    if R is None:
        R = max(60.0, 30.0 / k)
    if guess is None:
        # shooting guess at code frequency: find code w giving physical w
        wc = w
        for _ in range(4):
            p = bs.find_star(wc, eps)
            s = bs.star_properties(p, wc, eps)
            wc = wc * w / s['w_phys']
        p = bs.find_star(wc, eps)
        s = bs.star_properties(p, wc, eps)
        _, sh = bs.shoot(p, wc, eps)
        rc = s['rcut']
        x = np.concatenate([np.linspace(r0, rc * 0.7, 300), np.linspace(rc * 0.7 + 1e-3, R, 200)])
        # shoot() keeps no dense output; interpolate its steps
        from scipy.interpolate import interp1d
        f = interp1d(sh.t, sh.y, axis=1, fill_value='extrapolate')
        yi = f(np.minimum(x, rc * 0.7))
        lnal_inf = np.log(wc / w)                  # alpha_inf (code) = w_code / w_phys
        psi = np.log(np.abs(yi[0]))
        u = yi[1] / yi[0]
        m = yi[2]
        nu = 2 * (yi[3] - lnal_inf)
        out = x > rc * 0.7
        psi[out] = psi[~out][-1] - k * (x[out] - rc * 0.7)
        u[out] = -k - 1.0 / x[out]
        m[out] = s['M']
        nu[out] = np.log(1 - 2 * s['M'] / x[out])
        y0 = np.vstack([psi, u, m, nu])
    else:
        # resample the previous solution on a moderate mesh (re-using an over-refined mesh
        # makes solve_bvp's node count explode during continuation)
        xo = guess.sol.x
        idx = np.unique(np.linspace(0, xo.size - 1, min(xo.size, 1500)).astype(int))
        x = xo[idx]
        x = x[x <= R]
        if x[-1] < R:
            x = np.append(x, R)
        y0 = guess.sol.sol(np.minimum(x, guess.R))
        # extend tail if R grew
        far = x > guess.R
        if far.any():
            y0[1, far] = -k - 1.0 / x[far]
            y0[0, far] = y0[0, ~far][-1] - k * (x[far] - guess.R)
            y0[2, far] = y0[2, ~far][-1]
            y0[3, far] = np.log(1 - 2 * y0[2, far] / x[far])
    sol = solve_bvp(lambda r, y: odes(r, y, w, eps), lambda ya, yb: bcs(ya, yb, r0, R, w, eps),
                    x, y0, tol=tol, max_nodes=max_nodes, bc_tol=1e-10)
    if sol.status != 0:
        raise RuntimeError(f"bvp failed w={w} eps={eps}: {sol.message}")
    return Star(sol, w, eps, r0, R)

if __name__ == "__main__":
    for eps in (+1, -1):
        w = 0.75 if eps > 0 else 0.7603
        st = solve_star(w, eps)
        print(f"eps={eps:+d} w={w}: M={st.M:+.6f} N={st.N:.6f} phi_c={st.phic:.6f} alpha_c={st.alpha_c:.6f} "
              f"mesh={st.nmesh} max rms residual={st.max_res:.1e} R={st.R}")


# ---------------- continuation in particle number N (w is an unknown parameter) ----------------
# Needed to pass folds w(N) where fixed-w continuation must fail. Extra variable q(r) = enclosed N.
def odes_N(r, y, p, eps):
    w = p[0]
    f = odes(r, y[:4], w, eps)
    phi2 = np.exp(2 * y[0]); a = 1.0 / np.sqrt(1.0 - 2.0 * y[2] / r); al = np.exp(y[3] / 2)
    dq = 4 * np.pi * r * r * w * phi2 * a / al
    return np.vstack([f, dq])

def bcs_N(ya, yb, p, r0, R, eps, Nt):
    return np.concatenate([bcs(ya[:4], yb[:4], r0, R, p[0], eps), [ya[4], yb[4] - Nt]])

class StarN(Star):
    def __init__(self, sol, eps, r0, R):
        self.w = float(sol.p[0])
        Star.__init__(self, sol, self.w, eps, r0, R)
        self.N = float(sol.y[4, -1])

def solve_star_N(Nt, eps, guess, r0=1e-3, R=None, tol=1e-8, max_nodes=200000):
    """guess: a Star / StarN close by (continuation)."""
    R = guess.R if R is None else R
    xo = guess.sol.x
    idx = np.unique(np.linspace(0, xo.size - 1, min(xo.size, 1500)).astype(int))
    x = xo[idx]
    y0 = guess.sol.sol(x)[:4]
    al = np.exp(y0[3] / 2); a = 1 / np.sqrt(1 - 2 * y0[2] / x)
    from scipy.integrate import cumulative_trapezoid
    q = cumulative_trapezoid(4 * np.pi * x**2 * guess.w * np.exp(2 * y0[0]) * a / al, x, initial=0)
    q *= Nt / q[-1]
    sol = solve_bvp(lambda r, y, p: odes_N(r, y, p, eps), lambda ya, yb, p: bcs_N(ya, yb, p, r0, R, eps, Nt),
                    x, np.vstack([y0, q]), p=[guess.w], tol=tol, max_nodes=max_nodes, bc_tol=1e-10)
    if sol.status != 0:
        raise RuntimeError(f"bvp(N) failed N={Nt} eps={eps}: {sol.message}")
    return StarN(sol, eps, r0, R)


# ------------- continuation in the central lapse alpha_c (w free): passes folds in N and in w -------------
def bcs_lapse(ya, yb, p, r0, R, eps, nuc):
    return np.concatenate([bcs(ya, yb, r0, R, p[0], eps), [ya[3] - nuc]])

class StarP(Star):
    def __init__(self, sol, eps, r0, R):
        Star.__init__(self, sol, float(sol.p[0]), eps, r0, R)

def solve_star_lapse(alpha_c, eps, guess, r0=1e-3, R=None, tol=1e-8, max_nodes=200000):
    R = guess.R if R is None else R
    xo = guess.sol.x
    idx = np.unique(np.linspace(0, xo.size - 1, min(xo.size, 1500)).astype(int))
    x = xo[idx]
    y0 = guess.sol.sol(x)[:4]
    sol = solve_bvp(lambda r, y, p: odes(r, y, p[0], eps),
                    lambda ya, yb, p: bcs_lapse(ya, yb, p, r0, R, eps, 2 * np.log(alpha_c)),
                    x, y0, p=[guess.w], tol=tol, max_nodes=max_nodes, bc_tol=1e-10)
    if sol.status != 0:
        raise RuntimeError(f"bvp(alpha_c) failed alpha_c={alpha_c} eps={eps}: {sol.message}")
    return StarP(sol, eps, r0, R)

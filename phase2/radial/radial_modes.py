"""
Radial (l = 0) eigenvalue solver for signed-source sextic boson stars.

Reduced system (derived and verified symbolically in derive_radial.py; RHS in radial_rhs_generated.py):
    y = (A, A', Gt, chi)          A   = delta Re(phi) amplitude
                                  Gt  = i sigma B + w phi0 N/2   (gauge invariant; B = delta Im(phi))
                                  chi = Psi/phi0, Psi = perturbation of the enclosed Noether charge
    A''  = f_A(A, A', Gt, chi; sigma^2)
    Gt'  = u Gt + w phi0 N'/2 - sigma^2 (a/alpha) chi / r^2
    chi' = Psi'/phi0 - u chi
sigma^2 enters linearly; the metric perturbations are eliminated by the (tr) and (rr) constraints.

Boundary conditions
  r -> 0 : A = A0 + O(r^2), Gt = G0 + O(r^2), chi = O(r^3)  -> 2-parameter regular family.
  r -> oo: background phi0 ~ e^{-kr}/r; the system decouples into the free Klein-Gordon pair with
           wave numbers kappa_pm = sqrt(1 - (w +- sigma)^2); 2 decaying + 2 growing solutions.
           We impose decay by taking the 2 contracting eigenvectors of the local (frozen-coefficient)
           matrix at r = R, in variables (A, A', Gt, q = chi/r^2) in which the far-field matrix is
           asymptotically constant. Contamination by the growing pair decays inward like
           e^{-2 kappa (R - r)}, so for R - r_match >~ 15/kappa it is negligible (checked).
Method: shooting/matching. Integrate the 2 regular solutions out to r_m, the 2 decaying solutions in
to r_m, D(sigma^2) = det[Y1 Y2 Z1 Z2] / norms.  Eigenvalues = zeros of D.  sigma^2 < 0 => unstable.
Bound modes exist only for sigma < 1 - w (below the continuum of the field with frequency w+sigma).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import radial_rhs_generated as G
from background import LAM, MU

class RadialProblem:
    def __init__(self, star, r0=None, rm=None, R=None, rtol=1e-10):
        self.rtol = rtol
        self.st = star
        self.w, self.eps = star.w, star.eps
        self.r0 = star.r0 if r0 is None else r0
        self.R = star.R if R is None else R
        if rm is None:
            # matching radius: where phi0 has dropped to 1e-3 of its central value
            rr = np.linspace(self.r0, self.R, 4000)
            phi = self.st.bg(rr)[2]
            rm = rr[np.argmax(phi < 1e-3 * phi[0])]
        self.rm = rm

    def rhs(self, r, y, s2):
        al, a, phi0, u = self.st.bg(r)
        A, dA, Gt, chi = y
        args = (r, al, a, phi0, u, A, dA, Gt, chi, self.w, self.eps, s2, LAM, MU)
        return np.array([dA, G.ddA(*args), G.dGt(*args), G.dchi(*args)])

    # --- in variables (A, A', Gt, q = chi/r^2) for the far boundary ---
    def rhs_q(self, r, z, s2):
        y = np.array([z[0], z[1], z[2], z[3] * r * r])
        f = self.rhs(r, y, s2)
        return np.array([f[0], f[1], f[2], f[3] / (r * r) - 2 * z[3] / r])

    # ---- compound-matrix (bivector) formulation -----------------------------------------
    # For y' = J y, the 2-plane spanned by two solutions is carried by the antisymmetric
    # W = y1 y2^T - y2 y1^T, which obeys W' = J W + W J^T.  We integrate the normalized W
    # (W' = JW + WJ^T - c W, c chosen to keep |W| = 1), so the plane never collapses onto the
    # dominant direction (the failure of plain shooting with 2 decaying solutions of different
    # decay rates). D = det[y1 y2 z1 z2] = (1/4) eps_ijkl W_ij Z_kl.
    def J(self, r, s2, q=False):
        f = self.rhs_q if q else self.rhs
        Jm = np.empty((4, 4))
        for i in range(4):
            e = np.zeros(4); e[i] = 1.0
            Jm[:, i] = f(r, e, s2)          # the system is linear in y
        return Jm

    def _wrhs(self, r, wv, s2, q):
        W = wv.reshape(4, 4)
        Jm = self.J(r, s2, q)
        dW = Jm @ W + W @ Jm.T
        c = np.sum(dW * W) / np.sum(W * W)
        return (dW - c * W).ravel()

    def regular_plane(self, s2):
        r0 = self.r0
        al, a, phi0, u = self.st.bg(r0)
        ys = []
        for A0, G0 in ((1.0, 0.0), (0.0, 1.0)):
            chi0 = r0**3 * (G0 - 2 * self.w * A0) * a / (3 * al)   # Psi ~ r^3 phi0 (Gt-2wA) a/(3 alpha)
            f0 = G.ddA(r0, al, a, phi0, u, A0, 0.0, G0, chi0, self.w, self.eps, s2, LAM, MU)
            ys.append(np.array([A0 + f0 * r0**2 / 6, f0 * r0 / 3, G0, chi0]))
        W = np.outer(ys[0], ys[1]) - np.outer(ys[1], ys[0]); W /= np.linalg.norm(W)
        sol = solve_ivp(self._wrhs, (r0, self.rm), W.ravel(), args=(s2, False), method='DOP853',
                        rtol=self.rtol, atol=1e-14)
        return sol.y[:, -1].reshape(4, 4)

    def decaying_plane(self, s2):
        R = self.R
        # Real spectral projector onto the 2 contracting directions; P e_A, P e_Gt span the plane
        # with an orientation that varies smoothly with sigma^2 (raw eigenvectors have random signs).
        lam, vec = np.linalg.eig(self.J(R, s2, q=True))
        idx = np.argsort(lam.real)[:2]
        Vinv = np.linalg.inv(vec)
        P = (vec[:, idx] @ Vinv[idx, :]).real
        v1, v2 = P[:, 0], P[:, 2]
        W = np.outer(v1, v2) - np.outer(v2, v1); W /= np.linalg.norm(W)
        sol = solve_ivp(self._wrhs, (R, self.rm), W.ravel(), args=(s2, True), method='DOP853',
                        rtol=self.rtol, atol=1e-14)
        Wq = sol.y[:, -1].reshape(4, 4)
        T = np.diag([1.0, 1.0, 1.0, self.rm**2])      # q -> chi = r^2 q
        return T @ Wq @ T.T

    def D(self, s2):
        Y = self.regular_plane(s2)
        Z = self.decaying_plane(s2)
        Y /= np.linalg.norm(Y); Z /= np.linalg.norm(Z)
        # (1/4) eps_ijkl Y_ij Z_kl = Y01 Z23 - Y02 Z13 + Y03 Z12 + Y12 Z03 - Y13 Z02 + Y23 Z01
        return (Y[0, 1] * Z[2, 3] - Y[0, 2] * Z[1, 3] + Y[0, 3] * Z[1, 2]
                + Y[1, 2] * Z[0, 3] - Y[1, 3] * Z[0, 2] + Y[2, 3] * Z[0, 1])

    def eigenvalues(self, s2_lo, s2_hi, n=60, grid=None):
        if grid is None:
            grid = np.linspace(s2_lo, s2_hi, n)
        vals = np.array([self.D(s) for s in grid])
        roots = []
        for i in range(len(grid) - 1):
            if np.sign(vals[i]) != np.sign(vals[i + 1]):
                x = brentq(self.D, grid[i], grid[i + 1], xtol=1e-14, rtol=1e-12)
                # reject poles (sign change through infinity): |D| must be small at the root
                if abs(self.D(x)) < 1e-3 * max(abs(vals[i]), abs(vals[i + 1])):
                    roots.append(x)
        return roots, grid, vals


def lowest_modes(P, s2_min=-0.05, nneg=16, npos=16, nmax=3, s2_max=None):
    """Bound radial modes sigma^2 in [s2_min, (1-w)^2): sign changes of D on a grid that is
    log-dense near sigma^2 = 0 (where the fundamental mode crosses), refined by brentq.
    A sign change is accepted as a root only if |D| at the refined point is < 1e-6 of the bracket
    values (rejects jumps)."""
    thr = (1.0 - P.w)**2 if s2_max is None else s2_max
    g = np.concatenate([-np.logspace(np.log10(-s2_min), -7, nneg), [0.0],
                        np.logspace(-7, np.log10(0.2 * thr), npos // 2),
                        np.linspace(0.25 * thr, 0.995 * thr, npos // 2)])
    g = np.unique(g)
    vals = np.array([P.D(s) for s in g])
    roots = []
    for i in range(len(g) - 1):
        if np.sign(vals[i]) != np.sign(vals[i + 1]):
            x = brentq(P.D, g[i], g[i + 1], xtol=1e-13, rtol=1e-11)
            ok = abs(P.D(x)) < 1e-6 * max(abs(vals[i]), abs(vals[i + 1]))
            roots.append((x, ok))
            if len(roots) >= nmax: break
    return roots, g, vals

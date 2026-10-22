"""
Polar (even-parity) l >= 2 bound-mode solver for signed-source boson stars.

System (derive_nonradial.py -> gen_matrix2.py):  y' = J(r; s^2) y,  y = (H1/s, K, F, F', G/s, G'/s),
time dependence e^{st}; J is polynomial (linear) in s^2 and real for real s.  A growing mode is a
solution with s = gamma > 0 that is regular at r = 0 and decays at r -> oo.

Matching function (shooting):
  * interior: 3 regular solutions launched at r0 with the flat-space regular data r^l, integrated out
    to r_m with continuous QR re-orthonormalisation (the irregular r^{-l-1} branch is suppressed by
    (r0/r)^{2l+1}, and orthonormalisation keeps the 3-plane well conditioned even when the solutions
    grow like e^{gamma r});
  * exterior (r > r_m, where phi0 < ~1e-9 so the matter couplings ~ eps phi0^2 are < 1e-17): the metric
    and scalar blocks decouple exactly in vacuum (checked symbolically in gen_matrix2.py):
      - metric: decaying Zerilli solution Z ~ e^{-s r*}, obtained by integrating the Riccati variable
        q = f Z'/Z inward from r_far = r_m + 40/s (+ power-law zone) -- the decaying solution is the
        inward attractor; mapped to (H1/s, K) by the Moncrief map (verified symbolically);
      - scalar: 2 decaying free-KG solutions, from the contracting eigenvectors of the (frozen) scalar
        block at R_s = r_m + 40/Re(k), integrated inward with QR.
  * D(s) = det[Y1 Y2 Y3 Z1 Z2 Z3] at r_m with every 3-plane orthonormalised (positive-diagonal R) so D is
    a continuous real function of real s; zeros = eigenvalues.
"""
import numpy as np
from scipy.integrate import solve_ivp
from polar_matrix_generated import J_and_h
import background as B

LAM, MU = B.LAM, B.MU


class PolarProblem:
    def __init__(self, star, l, rm=None, r0=None, rtol=1e-10, LAM=LAM, MU=MU, phi_cut=1e-9):
        self.st, self.l, self.LL = star, l, l * (l + 1)
        self.w, self.eps = star.w, star.eps
        self.rtol = rtol
        self.Lam, self.Mu = LAM, MU
        self.r0 = star.r0 if r0 is None else r0
        if rm is None:
            rr = np.linspace(self.r0, star.R, 20000)
            phi = star.bg(rr)[2]
            rm = rr[np.argmax(phi < phi_cut)]
        self.rm = rm
        al, a, phi, u = star.bg(np.array([rm]))
        self.Mext = float(rm * (1 - 1 / a[0] ** 2) / 2)       # m(r_m) = exterior mass
        self.alpha_m = float(al[0])
        self.lapse_ratio = float(al[0] / np.sqrt(1 - 2 * self.Mext / rm))  # should be 1 (alpha(oo)=1)

    # ---------------------------------------------------------------- interior
    def J(self, r, s):
        al, a, phi, u = self.st.bg(np.atleast_1d(r))
        Jm, h = J_and_h(r, al[0], a[0], phi[0], u[0], self.w, self.eps, s, self.LL, self.Lam, self.Mu)
        return Jm

    def _rhs_plane(self, r, yv, s, k):
        Y = yv.reshape(6, k)
        return (self.J(r, s) @ Y).ravel()

    @staticmethod
    def _qr(Y):
        """QR with the phase of the R diagonal returned separately: Y = Q R, ph = prod(R_ii/|R_ii|).
        Keeping ph makes det[...] carry exactly the phase of the holomorphic (unnormalised) determinant."""
        Q, R = np.linalg.qr(Y)
        d = np.diag(R)
        return Q, np.prod(d / np.abs(d))

    def _integrate_plane(self, Y, ra, rb, s, nchunks=None, Jfun=None):
        k = Y.shape[1]
        Jfun = self.J if Jfun is None else Jfun
        if nchunks is None:
            nchunks = int(max(20, abs(rb - ra) * max(1.0, abs(s)) / 2.0))
        grid = np.linspace(ra, rb, nchunks + 1)
        Q, ph = self._qr(Y)
        for i in range(nchunks):
            sol = solve_ivp(lambda r, yv: (Jfun(r, s) @ yv.reshape(6, k)).ravel(), (grid[i], grid[i + 1]),
                            Q.ravel(), method='DOP853', rtol=self.rtol, atol=1e-14)
            Q, p = self._qr(sol.y[:, -1].reshape(6, k))
            ph *= p
        return Q, ph

    def regular_plane(self, s):
        l, r0 = self.l, self.r0
        # flat-space regular data at small r: metric (K ~ r^l, H1/s ~ r^{l+1}) and two scalar branches
        Y = np.zeros((6, 3))
        Y[0, 0] = r0 ** (l + 1) / l; Y[1, 0] = r0 ** l
        Y[2, 1] = r0 ** l; Y[3, 1] = l * r0 ** (l - 1)
        Y[4, 2] = r0 ** l; Y[5, 2] = l * r0 ** (l - 1)
        Y /= np.linalg.norm(Y, axis=0)
        return self._integrate_plane(Y.astype(complex) if np.iscomplexobj(s) else Y, r0, self.rm, s)

    # ---------------------------------------------------------------- exterior: metric (Zerilli)
    def zerilli_q(self, s, r_start=None):
        """q = f dZ/dr / Z = dlnZ/dr* at r_m for the decaying Zerilli solution (mass Mext)."""
        M, lam = self.Mext, (self.l - 1) * (self.l + 2) / 2.0
        def V(r):
            f = 1 - 2 * M / r
            return 2 * f * (lam ** 2 * r ** 2 * ((lam + 1) * r + 3 * M) + 9 * M ** 2 * (lam * r + M)) / (r ** 3 * (lam * r + 3 * M) ** 2)
        # integrate in rho = ln r:  dq/drho = r (s^2 + V - q^2)/f
        if r_start is None:
            r_start = self.rm + 40.0 / abs(s) + 20.0 * self.l / abs(s) + 50.0
        # start on the asymptotic decaying solution, including the leading power-law correction
        # Z ~ e^{-s r*} (1 + c/r): q ≈ -s - LL/(2 r^2 s)... use the 2nd-order WKB/Riccati asymptote
        LLv = self.LL
        q0 = -np.sqrt(s * s + V(r_start))
        def rhs(rho, q):
            r = np.exp(rho); f = 1 - 2 * M / r
            return r * (s * s + V(r) - q * q) / f
        if np.iscomplexobj(q0) or np.iscomplexobj(s):
            q0 = complex(q0); q0 = q0 if q0.real <= 0 else -q0   # decaying / outgoing branch q ~ -s
            q0 = -np.sqrt(complex(s * s + V(r_start)))
            if (q0 / s).real > 0: q0 = -q0
        sol = solve_ivp(rhs, (np.log(r_start), np.log(self.rm)), [q0], method='DOP853', rtol=1e-11, atol=1e-14)
        return sol.y[0, -1] if np.iscomplexobj(sol.y) else float(sol.y[0, -1])

    def moncrief(self, r):
        """Matrix T with (H1/s, K)^T = T (Z, dZ/dr)^T in Schwarzschild of mass Mext (gen_matrix2 check (b))."""
        M, lam = self.Mext, (self.l - 1) * (self.l + 2) / 2.0
        f = 1 - 2 * M / r
        return np.array([[(lam * r * r - 3 * lam * M * r - 3 * M * M) / ((r - 2 * M) * (lam * r + 3 * M)), r],
                         [(lam * (lam + 1) * r * r + 3 * lam * M * r + 6 * M * M) / (r * r * (lam * r + 3 * M)), f]])

    # ---------------------------------------------------------------- exterior: scalar
    def scalar_plane(self, s):
        """Decaying scalar 2-plane at r_m (4x2) and its holomorphic phase."""
        w = self.w
        kr = min(np.real(np.sqrt(1 + (s - 1j * w) ** 2)), np.real(np.sqrt(1 + (s + 1j * w) ** 2)))
        Rs = self.rm + 40.0 / kr
        if getattr(self, 'Rs_fixed', None) is not None:      # for contour integrals: s-independent radius
            Rs = self.Rs_fixed
        Js = self._Jvac(Rs, s)[2:, 2:]
        lam, vec = np.linalg.eig(Js)
        idx = np.argsort(lam.real)[:2]
        Vinv = np.linalg.inv(vec)
        P = vec[:, idx] @ Vinv[idx, :]                     # spectral projector: holomorphic in s
        if not np.iscomplexobj(s):
            P = P.real
        Y = np.zeros((6, 2), dtype=P.dtype); Y[2:, 0] = P[:, 0]; Y[2:, 1] = P[:, 2]
        jf = lambda r, s_: self._Jvac(r, s_) if r > self.st.R else self.J(r, s_)
        Q, ph = self._integrate_plane(Y, Rs, self.rm, s, Jfun=jf)
        return Q[2:, :], ph

    def _Jvac(self, r, s):
        f = 1 - 2 * self.Mext / r
        Jm, h = J_and_h(r, np.sqrt(f), 1 / np.sqrt(f), 0.0, 0.0, self.w, self.eps, s, self.LL, self.Lam, self.Mu)
        return Jm

    # ---------------------------------------------------------------- matching
    def decaying_plane(self, s):
        q = self.zerilli_q(s)
        f = 1 - 2 * self.Mext / self.rm
        zv = np.array([1.0, q / f])                       # (Z, dZ/dr), holomorphic in s
        m = self.moncrief(self.rm) @ zv
        S, ph = self.scalar_plane(s)
        Z = np.zeros((6, 3), dtype=np.result_type(m, S))
        Z[0:2, 0] = m
        Z[2:, 1:] = S
        Qz, p = self._qr(Z)
        return Qz, ph * p

    def D(self, s):
        """det[regular 3-plane | decaying 3-plane] at r_m, orthonormalised, times the exact phase of the
        holomorphic determinant (so arg D is the arg of an analytic function of s on Re s > 0; for real
        s, D is real and continuous)."""
        Y, py = self.regular_plane(s)
        Z, pz = self.decaying_plane(s)
        d = np.linalg.det(np.hstack([Y, Z])) * py * pz
        return d if np.iscomplexobj(d) and not np.isrealobj(s) else float(np.real(d))

    # ---------------------------------------------------------------- effective surface condition
    def interior_zerilli(self, s):
        """Regular interior solution that also matches the decaying scalar exterior: returns
        (Z, dZ/dr) at r_m (up to normalisation) from the interior 3-plane. Generic when the scalar
        sector is not itself at an eigenvalue."""
        Y, _ = self.regular_plane(s)
        S, _ = self.scalar_plane(s)                       # 4x2
        # find c in R^3 with Y[2:] c in span(S):  (I - S S^+) Y[2:] c = 0
        Ps = np.eye(4) - S @ np.linalg.pinv(S)
        A = Ps @ Y[2:, :]
        _, sv, vt = np.linalg.svd(A)
        c = vt[-1]
        met = Y[:2, :] @ c
        zz = np.linalg.solve(self.moncrief(self.rm), met)
        return zz, sv

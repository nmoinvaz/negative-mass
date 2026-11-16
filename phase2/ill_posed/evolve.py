"""
Time-domain evolution of  u_tt = u_xx - V(x) u  on a grid x_0 < x_1 < ... < x_n (Dirichlet at x_n), by P1 Galerkin with a lumped
mass (a symmetric second-order scheme whose quadratic form is  sum (du)^2/h + sum V m u^2 + kappa u_0^2)  and leapfrog.

  left = 'dirichlet'  u(x_0) = 0, x_0 may be 0 (the point mass, V ~ 3/(4 x^2): the regular branch x^{3/2} is selected by the
                      discretisation, which is the Friedrichs extension, the DG choice);
  left = 'robin'      u_x = kappa u at x_0 (the truncated star of dg_proof Theorem 2).

The solution is rescaled during the run so that growth factors of e^{1000} are representable; only log norms are reported.
Independent of dg_spectrum (no eigensolve, no Chebyshev): the data are advanced in time, which is the question asked.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal


class Op:
    def __init__(self, x, V, left='dirichlet', kappa=0.0):
        x = np.asarray(x, float)
        V = np.asarray(V, float)
        n = len(x) - 1
        h = np.diff(x)
        lo = 0 if left == 'robin' else 1
        self.lo, self.x, self.h, self.left, self.kappa = lo, x, h, left, kappa
        idx = np.arange(lo, n)                       # unknown nodes lo .. n-1
        hl = np.where(idx >= 1, h[np.maximum(idx - 1, 0)], np.nan)
        hr = h[idx]
        d = np.where(idx >= 1, 1.0 / hl, 0.0) + 1.0 / hr
        m = np.where(idx >= 1, 0.5 * (hl + hr), 0.5 * hr)
        if left == 'robin':
            d[0] += kappa
        self.m = m
        self.d = d + V[idx] * m                      # diagonal of K (with potential)
        self.e = -1.0 / h[idx[:-1]]                  # K_{j,j+1}
        self.invm = 1.0 / m
        self.xs = x[idx]
        self.V = V[idx]
        lam = (self.d + np.abs(np.r_[self.e, 0.0]) + np.abs(np.r_[0.0, self.e])) * self.invm
        self.lam_max = float(lam.max())              # Gershgorin bound on the largest eigenvalue of M^-1 K
        self.dt = 1.0 / np.sqrt(self.lam_max)        # dt^2 lam_max <= 1 (< 4: leapfrog stable)

    def apply(self, u):
        Ku = self.d * u
        Ku[1:] += self.e * u[:-1]
        Ku[:-1] += self.e * u[1:]
        return Ku * self.invm

    def norm2(self, u):
        return float(np.sum(self.m * u * u))

    def lowest(self, k=1):
        dd = self.d * self.invm
        ee = self.e / np.sqrt(self.m[:-1] * self.m[1:])
        return eigh_tridiagonal(dd, ee, eigvals_only=True, select='i', select_range=(0, k - 1))

    def count_negative(self):
        """number of negative eigenvalues (Sylvester inertia of K, pivots of the tridiagonal LDL^T)."""
        d, e = self.d, self.e
        p = d[0]
        cnt = int(p < 0)
        for j in range(1, len(d)):
            p = d[j] - e[j - 1]**2 / p
            cnt += int(p < 0)
        return cnt


def evolve(op, u0, v0=None, T=1.0, nrec=8, dt=None, dt_frac=0.5):
    """advance (u0, v0) to time T; return times, ln||u(t)||_{L^2(dx)} at nrec equally spaced times (and u at the end, rescaled)."""
    dt = (dt if dt is not None else dt_frac * op.dt)
    nsteps = max(int(np.ceil(T / dt)), nrec)
    dt = T / nsteps
    u0 = np.asarray(u0, float).copy()
    prev = u0
    cur = u0 + (0.0 if v0 is None else dt * np.asarray(v0, float)) - 0.5 * dt * dt * op.apply(u0)
    logscale = 0.0
    rec_at = {int(round(nsteps * (i + 1) / nrec)): i for i in range(nrec)}
    times = np.zeros(nrec)
    logn = np.zeros(nrec)
    n0 = 0.5 * np.log(op.norm2(u0))
    for step in range(1, nsteps + 1):
        if step in rec_at:
            i = rec_at[step]
            times[i] = step * dt
            logn[i] = 0.5 * np.log(op.norm2(cur)) + logscale
        if step == nsteps:
            break
        nxt = 2.0 * cur - prev - dt * dt * op.apply(cur)
        prev, cur = cur, nxt
        if step % 64 == 0:
            s = np.max(np.abs(cur))
            if s > 1e100:
                cur = cur / s
                prev = prev / s
                logscale += np.log(s)
    return times, logn - n0, cur, logscale


def ln_cosh(y):
    y = abs(y)
    return y - np.log(2.0) + np.log1p(np.exp(-2.0 * y))

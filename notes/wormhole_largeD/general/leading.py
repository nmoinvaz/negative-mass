"""Leading large-D eigenvalue problem for the asymmetric / charged ghost wormhole (Lam^2 = gamma1^2 - eps2 != 0).

With a = 1, N1 = D - 3, W_th = N1 sqrt(omega) + ..., L = N1 Lt (L the metric perturbation delta lambda), the near zone
-pi/2 < y < pi/2 gives (nearzone.py; Liouville nu'' = eps2 e^{2 nu}, nu'^2 - nu'' = Lam^2):
  Y'  = nu' Y + B' Lt,                         B' = tan y - nu'
  cos^2(y) B' Lt' = (omega e^{2(nu_th - nu)} - Lam^2 cos^2 y) Y + (1 - sin y cos y nu' + Lam^2 cos^2 y) Lt.
At each end (delta = distance to +-pi/2) the far zone's decaying solution requires Y ~ delta^p,
p = (-1 + sqrt(1 + 4 kappa^2))/2, kappa^2 = omega e^{2(nu_th - nu_end)}. The throat (B' = 0) is a regular singular
point; the two solutions are integrated from the ends, joined along a semicircle around the throat in the complex y
plane (upper and lower paths are compared: equal results mean no logarithm), and the eigenvalue is a zero of the
2 x 2 determinant.

    python3 leading.py gamma1 eps2 [omega_guess]
"""
import sys

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

mp.mp.dps = 25


class Leading:
    def __init__(self, g1, e2):
        self.g1, self.e2 = g1, e2
        self.Lam2 = g1 * g1 - e2
        self.Lam = np.sqrt(complex(self.Lam2))
        self.yth = float(mp.findroot(lambda yy: mp.tan(yy) - mp.re(self.dnu_mp(yy)), mp.atan(g1)))
        self.nu_th = float(np.real(self.nu(self.yth)))
        self.nu_end = [float(np.real(self.nu(s * (np.pi / 2)))) for s in (-1, 1)]

    def ch(self, yy):
        L = self.Lam
        if abs(L) < 1e-12:
            return 1 - self.g1 * yy, -self.g1
        return np.cosh(L * yy) - self.g1 * np.sinh(L * yy) / L, L * np.sinh(L * yy) - self.g1 * np.cosh(L * yy)

    def nu(self, yy):
        return -np.log(self.ch(yy)[0])

    def dnu(self, yy):
        c, cp = self.ch(yy)
        return -cp / c

    def dnu_mp(self, yy):
        L = mp.sqrt(mp.mpc(self.Lam2))
        if abs(L) < 1e-12:
            return self.g1 / (1 - self.g1 * yy)
        c = mp.cosh(L * yy) - self.g1 * mp.sinh(L * yy) / L
        cp = L * mp.sinh(L * yy) - self.g1 * mp.cosh(L * yy)
        return -cp / c

    def rhs(self, yy, X, om):
        Y, Lt = X
        nup = self.dnu(yy)
        Bp = np.tan(yy) - nup
        c2 = np.cos(yy) ** 2
        A = om * np.exp(2 * (self.nu_th - self.nu(yy))) - self.Lam2 * c2
        C = 1 - np.sin(yy) * np.cos(yy) * nup + self.Lam2 * c2
        return np.array([nup * Y + Bp * Lt, (A * Y + C * Lt) / (c2 * Bp)])

    def end_start(self, side, om, d0):
        """decaying solution near the end y = side pi/2, at distance d0: Y = d^p, Lt from Y' = nu' Y + B' Lt"""
        kap2 = om * np.exp(2 * (self.nu_th - self.nu_end[0 if side < 0 else 1]))
        p = (-1 + np.sqrt(1 + 4 * kap2 + 0j)) / 2
        yy = side * (np.pi / 2 - d0)
        Y = 1.0 + 0j                                 # the scale is arbitrary; d0^p can underflow
        dYdy = -side * p / d0                        # d/dy = -side d/d(delta), Y ~ delta^p
        Lt = (dYdy - self.dnu(yy) * Y) / (np.tan(yy) - self.dnu(yy))
        return yy, np.array([Y, Lt], dtype=complex)

    def integrate(self, y0, X0, y1, om, path_sign=+1, rho=None):
        """real segment y0 -> y1. For real omega: Pruefer angle theta (Y = R cos theta, Lt = R sin theta), which stays
        finite however much the solution grows; for complex omega: linear system in renormalized chunks."""
        if abs(np.imag(om)) < 1e-14 and abs(np.imag(X0[0] * np.conj(X0[1]))) < 1e-12 * abs(X0[0] * X0[1]) + 1e-300:
            omr = float(np.real(om))
            th0 = np.arctan2(np.real(X0[1] / X0[0]) if abs(X0[0]) > 0 else 1.0, 1.0)

            def g(yy, th):
                c, s_ = np.cos(th[0]), np.sin(th[0])
                Yp, Lp = self.rhs(yy, np.array([c, s_]), omr).real
                return [-s_ * Yp + c * Lp]
            sol = solve_ivp(g, (y0, y1), [th0], method="DOP853", rtol=1e-11, atol=1e-13)
            th = sol.y[0, -1]
            return np.array([np.cos(th), np.sin(th)], dtype=complex)
        f = lambda yy, X: self.rhs(yy, X, om)  # noqa: E731
        edges = np.linspace(y0, y1, 400)
        X = np.array(X0, dtype=complex)
        for a, b in zip(edges[:-1], edges[1:]):
            sol = solve_ivp(lambda t, Z: f(t, Z), (a, b), X, method="DOP853", rtol=1e-11, atol=1e-30)
            X = sol.y[:, -1]
            X = X / np.linalg.norm(X)
        return X

    def around(self, X0, om, rho, path_sign):
        """from y_th - rho to y_th + rho along a semicircle (upper for path_sign = +1)"""
        def g(th, X):
            yy = self.yth + rho * np.exp(1j * th)
            return self.rhs(yy, X, om) * (1j * rho * np.exp(1j * th))
        th0, th1 = (np.pi, 0.0) if path_sign > 0 else (-np.pi, 0.0)
        sol = solve_ivp(g, (th0, th1), X0, method="DOP853", rtol=1e-11, atol=1e-30)
        return sol.y[:, -1]

    def det(self, om, d0=1e-7, rho=None, path_sign=+1):
        rho = rho or 0.25 * min(np.pi / 2 - abs(self.yth), 0.5)
        yl, Xl = self.end_start(-1, om, d0)
        Xl = self.integrate(yl, Xl, self.yth - rho, om)
        Xl = self.around(Xl, om, rho, path_sign)
        yr, Xr = self.end_start(+1, om, d0)
        Xr = self.integrate(yr, Xr, self.yth + rho, om)
        nl = np.linalg.norm(Xl)
        nr = np.linalg.norm(Xr)
        return (Xl[0] * Xr[1] - Xl[1] * Xr[0]) / (nl * nr)

    def eigen_real(self, lo=1e-4, hi=10.0, N=60):
        """real eigenvalues: sign changes of the (real) determinant on a log grid, refined by brentq"""
        from scipy.optimize import brentq
        grid = np.geomspace(lo, hi, N)
        vals = [float(np.real(self.det(om))) for om in grid]
        roots = []
        for (a, fa), (b, fb) in zip(zip(grid, vals), zip(grid[1:], vals[1:])):
            if np.isfinite(fa) and np.isfinite(fb) and np.sign(fa) != np.sign(fb):
                roots.append(brentq(lambda om: float(np.real(self.det(om))), a, b, xtol=1e-14, rtol=1e-13))
        return roots

    def eigen_complex(self, re=(-2.0, 6.0, 0.25), im=(0.05, 4.0, 0.25), n=4):
        """complex eigenvalues (Im omega > 0; their conjugates are roots too): local minima of |det| on a grid,
        refined by eigen. For the symmetric charged wormhole the leading-order eigenvalue is complex (gamma1 = 0,
        eps2 = 0.05: omega = -0.0749 + 0.2201i), so eigen_real finds none."""
        rs, ims = np.arange(re[0], re[1] + 1e-9, re[2]), np.arange(im[0], im[1] + 1e-9, im[2])
        A = np.array([[abs(self.det(complex(a, b))) for a in rs] for b in ims])
        mins = sorted((A[i, j], complex(rs[j], ims[i])) for i in range(len(ims)) for j in range(len(rs))
                      if A[i, j] == A[max(i - 1, 0):i + 2, max(j - 1, 0):j + 2].min() and A[i, j] < 0.5 * np.median(A))
        roots = []
        for _, z in mins[:n]:
            try:
                om = self.eigen(z)
            except Exception:                   # findroot did not converge from this minimum
                continue
            if om.imag > 0 and all(abs(om - r) > 1e-8 for r in roots):
                roots.append(om)
        return roots

    def eigen(self, guess):
        f = lambda om: complex(self.det(complex(om)))  # noqa: E731
        return complex(mp.findroot(f, mp.mpc(guess), tol=1e-20, maxsteps=100))


if __name__ == "__main__":
    g1, e2 = float(sys.argv[1]), float(sys.argv[2])
    guess = complex(sys.argv[3]) if len(sys.argv) > 3 else None
    L = Leading(g1, e2)
    if guess is not None:
        oms = [L.eigen(guess)]
    else:
        oms = L.eigen_real()
        print("real roots:", oms)
        if not oms:
            oms = L.eigen_complex()
            print("no real root; complex roots (Im omega > 0):", oms)
    print(f"gamma1 = {g1}, eps2 = {e2}: y_th = {L.yth:.6f}, nu_th = {L.nu_th:.6f}, nu_ends = {L.nu_end}")
    for om in oms or [float("nan")]:
        w = np.sqrt(complex(om))
        print(f"  omega = {om:.12f}  -> W_th a ~ (D-3) sqrt(omega) = (D-3) * {w:.10f}"
              + (f" (growth {w.real:.10f}, oscillation {w.imag:.10f})" if abs(w.imag) > 1e-12 else ""))
        if om == om:
            print(f"  upper/lower path determinant at the root: {abs(L.det(om, path_sign=1)):.1e}, "
                  f"{abs(L.det(om, path_sign=-1)):.1e}")

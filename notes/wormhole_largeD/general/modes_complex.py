"""Complex growth rates of the general ghost wormhole (modes_general.py's Wormhole): W = W_r + i W_i, a growth e^{W t}
that oscillates when W_i != 0. The charged symmetric wormhole's growth turns oscillating as the charge or D grows
(charged_scan.py: its two real even roots meet and leave the real axis), and the real determinant of
Wormhole.det cannot follow it.

Method: shoot inward from each end and match at the throat, so that everything is analytic in W.
  1. At each end, at z = margin (z = pi/2 - side y; r = rmax_factor / Re(W_loc) there, as in Wormhole.det), start on the
     solution that decays outward: the eigenvector, with the larger real eigenvalue, of the local first-order system in
     t = log z for (Y, V), V = Y'/g, g = sqrt(1 + r'^2) (Wormhole.shoot), normalised to Y = 1. An error in that
     direction is the solution that grows outward, which decays inward by about e^{-2 rmax_factor}.
  2. Integrate the linear system inward (t increasing) to z0 from the throat, complex, unnormalised: the solution grows
     inward by e^{Re(W) r}, about e^{60..150}, well inside double range.
  3. Write it there as c0 Y_0 + c3 Y_3 in the throat's Frobenius solutions (exponents 0 and 3; Wormhole.frobenius with
     complex W). The two decaying solutions are the same solution when c0_L c3_R - c3_L c0_R = 0; the root search uses
     the equivalent mismatch c3_R / c0_R - c3_L / c0_L, from which the inward growth e^{W r} (which makes the
     determinant itself vary by 10x per 0.1 in W_hat) has cancelled. (For a symmetric wormhole the even modes are
     c3 = 0, where this mismatch vanishes; an odd mode, c0 = 0, would need the inverse ratios.)
The end points (margins) are fixed during a root search (from the starting W), so the mismatch is analytic in W.

    python3 modes_complex.py D gamma1 eps2 W_hat_guess        (W_hat = W_th a, complex, e.g. 2.0+2.5j)
"""
import sys

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

import modes_general as mg


class ComplexWormhole(mg.Wormhole):
    def frobenius_c(self, W, K=14, rho=None):
        """Wormhole.frobenius for complex W: the Taylor coefficients stay complex"""
        if self.fast:
            return self.frobenius_fast(W, K, rho, real=False)
        p = (self.n, self.c, self.g1, self.e2)
        rho = rho or 0.05 * min(1.0, np.pi / 2 - abs(self.yth))
        M = 64
        W = mp.mpc(W)
        z = [rho * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M)]
        Pv, Qv = [], []
        for zz in z:
            A2, A1, A0 = (f(self.yth + zz, *p, W) for f in self.f_mp)
            Pv.append(zz * A1 / A2)
            Qv.append(zz * zz * A0 / A2)
        Pk = [mp.fsum(Pv[j] * z[j] ** (-k) for j in range(M)) / M for k in range(K + 1)]
        Qk = [mp.fsum(Qv[j] * z[j] ** (-k) for j in range(M)) / M for k in range(K + 1)]
        sols = {}
        for r0 in (0, 3):
            cks = [mp.mpc(1)]
            for k in range(1, K + 1 - r0):
                ind = (k + r0) * (k + r0 - 1) + Pk[0] * (k + r0) + Qk[0]
                rhs = -mp.fsum((Pk[j] * (k - j + r0) + Qk[j]) * cks[k - j] for j in range(1, k + 1))
                if abs(ind) < 1e-12:            # resonance (r0 = 0, k = 3): log-free when rhs = 0
                    sols["resonance_rhs"] = complex(rhs)
                    cks.append(mp.mpc(0))
                else:
                    cks.append(rhs / ind)
            sols[r0] = [complex(v) for v in cks]
        return sols

    def coeffs_zc(self, zz, W, side):
        p = (self.n, self.c, self.g1, complex(self.e2))
        if self.fast:
            return list((self.F.end_L if side < 0 else self.F.end_R)(complex(zz), *p, complex(W))[:3])
        return [complex(np.asarray(f(zz + 0j, *p, complex(W)))) for f in mg.end_functions(side)[0]]

    def system(self, zz, W, side):
        """d(Y, V)/dt = z M(z) (Y, V), t = log z (Wormhole.shoot's system, linear and complex)"""
        if self.fast:
            p = (self.n, self.c, self.g1, complex(self.e2))
            return (self.F.lin_L if side < 0 else self.F.lin_R)(complex(zz), *p, complex(W))
        A2, A1, A0 = self.coeffs_zc(zz, W, side)
        rz, rzz = float(self.r_z(zz, side, 1)), float(self.r_z(zz, side, 2))
        g2 = 1.0 + rz * rz
        g = np.sqrt(g2)
        return zz * np.array([[0.0, -side * g], [side * A0 / (A2 * g), side * A1 / A2 - rz * rzz / g2]], dtype=complex)

    def margin(self, W, side, rmax_factor=60.0, X=30.0):
        Wloc = abs(complex(W).real) * np.exp(-self.nu_end[0 if side < 0 else 1])
        target = max(rmax_factor / Wloc, self.a * (1 + X / (2 * (self.n - 1))))
        lo, hi = -700.0, np.log(np.pi / 2 - side * self.yth)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if float(self.r_z(np.exp(mid), side)) > target:
                lo = mid
            else:
                hi = mid
        return np.exp(lo)

    def side_coeffs(self, W, side, margin, sols, z0, rtol=1e-11):
        """(c0, c3) of the solution decaying at the end side, in the Frobenius basis at y_th + side z0"""
        M = self.system(margin, W, side) / margin
        ev, V = np.linalg.eig(M)
        j = int(np.argmax(ev.real))
        X0 = V[:, j] / V[0, j]
        zthr = np.pi / 2 - side * self.yth - z0
        sol = solve_ivp(lambda t, X: self.system(np.exp(t), W, side) @ X, (np.log(margin), np.log(zthr)), X0,
                        method="DOP853", rtol=rtol, atol=1e-300)
        Y, Vv = sol.y[:, -1]
        g0 = np.sqrt(1.0 + float(self.r_z(zthr, side, 1)) ** 2)
        B = np.empty((2, 2), dtype=complex)
        for col, r0 in enumerate((0, 3)):
            Yk, Ypk = self.local(sols, r0, side * z0)
            B[:, col] = (Yk, Ypk / g0)
        return np.linalg.solve(B, np.array([Y, Vv])), sol.y[0, -1]

    def det_c(self, W2, margins=None, z0=None, rtol=1e-11, rmax_factor=60.0, X=30.0):
        """the mismatch c3_R / c0_R - c3_L / c0_L at W = sqrt(W2) (complex); margins (left, right) fixed if given"""
        W = np.sqrt(complex(W2))
        sols = self.frobenius_c(W)
        z0 = z0 or 0.02 * min(1.0, np.pi / 2 - abs(self.yth))
        if margins is None:
            margins = tuple(self.margin(W, s, rmax_factor, X) for s in (-1, 1))
        (cL, _), (cR, _) = (self.side_coeffs(W, s, m, sols, z0, rtol) for s, m in zip((-1, 1), margins))
        return cR[1] / cR[0] - cL[1] / cL[0]

    def root_c(self, What_guess, rmax_factor=60.0, X=30.0, rtol=1e-11, tol=1e-24):
        """complex root What = W_th a near What_guess (mp.findroot, secant), margins fixed from the guess"""
        W = What_guess * np.exp(self.nu_th) / self.a
        margins = tuple(self.margin(W, s, rmax_factor, X) for s in (-1, 1))
        scale = self.a * np.exp(-self.nu_th)
        f = lambda h: self.det_c(complex(h / scale) ** 2, margins, rtol=rtol, rmax_factor=rmax_factor, X=X)  # noqa
        g = complex(What_guess)
        r = mp.findroot(lambda h: mp.mpc(f(complex(h))), (mp.mpc(g), mp.mpc(g * (1 + 1e-4))), solver="secant",
                        tol=tol, maxsteps=60)
        return complex(r)


if __name__ == "__main__":
    D, g1, e2, guess = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), complex(sys.argv[4])
    w = ComplexWormhole(D, g1, e2)
    r = w.root_c(guess)
    print(f"D = {D}, gamma1 = {g1}, eps2 = {e2}: What = {r.real:.10f} {r.imag:+.10f}i "
          f"(growth {r.real:.6f}, oscillation {r.imag:.6f} in throat units)")

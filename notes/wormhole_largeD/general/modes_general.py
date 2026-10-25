"""Growing radial mode of the general static ghost wormhole (asymmetry gamma1, charge eps2) in D = n + 2 dimensions.

Master equation A2 Y'' + A1 Y' + A0 Y = 0 in the harmonic coordinate y in (-pi/2, pi/2) (master_general.py). The throat
y_th (tan y = nu'(y)) is a regular singular point with exponents 0 and 3. Near it the two Frobenius solutions are built
from Taylor coefficients of z A1/A2 and z^2 A0/A2 (z = y - y_th), computed by Cauchy integrals on a small circle; the
exponent-0 solution has a resonance at order 3, whose consistency (no logarithm) is checked and printed. Both
solutions are integrated to each end; the growing part there is measured by Y e^{-W_loc r} (W_loc = W e^{-nu(end)}, the
growth rate in that end's asymptotic time), and the eigenvalue is the zero of the 2 x 2 determinant.

Output: What = W_th a, the growth rate in proper time at the throat, times the throat's areal radius a (as in
Gonzalez, Guzman & Sarbach 2009, whose T = 1/What).

    python3 modes_general.py D gamma1 eps2 [What_lo What_hi]
"""
import functools
import sys

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import derive_general as d
import master_general as mg

mp.mp.dps = 30


@functools.lru_cache(maxsize=None)
def coefficient_functions():
    a2, a1, a0, *_ = mg.master_coeffs()
    nu, lam, r, P, qq = mg.background_exprs()
    bg = {d.nu0: nu, d.lam0: lam, d.r0: r, d.P: P, d.q: qq}
    args = (d.x, mg.n, mg.c, mg.g1, mg.e2, d.W)
    exprs = [sp.together((v.subs(bg).doit())) for v in (a2, a1, a0)]
    f_np = [sp.lambdify(args, e, "numpy") for e in exprs]
    f_mp = [sp.lambdify(args, e, "mpmath") for e in exprs]
    bgf_np = [sp.lambdify((d.x, mg.n, mg.c, mg.g1, mg.e2), e, "numpy") for e in (nu, r, sp.diff(r, d.x))]
    dnu = sp.lambdify((d.x, mg.n, mg.c, mg.g1, mg.e2), sp.diff(nu, d.x), "mpmath")
    nu_mp = sp.lambdify((d.x, mg.n, mg.c, mg.g1, mg.e2), nu, "mpmath")
    r_mp = sp.lambdify((d.x, mg.n, mg.c, mg.g1, mg.e2), r, "mpmath")
    return f_np, f_mp, bgf_np, dnu, nu_mp, r_mp


@functools.lru_cache(maxsize=None)
def end_functions(side):
    """A2, A1, A0 and r as functions of z = pi/2 - side y, the distance to the end y = side pi/2. y enters the
    coefficients through cos y, which becomes sin z exactly, so z keeps full relative precision however close the
    integration gets to the end (in y, pi/2 - y = 7e-10 at D = 10 carries only 3e-7 relative precision)."""
    a2, a1, a0, *_ = mg.master_coeffs()
    nu, lam, r, P, qq = mg.background_exprs()
    bg = {d.nu0: nu, d.lam0: lam, d.r0: r, d.P: P, d.q: qq}
    z = sp.Symbol("z", positive=True)
    sub = {d.x: side * (sp.pi / 2 - z)}
    exprs = [sp.together(v.subs(bg).doit()).subs(sub) for v in (a2, a1, a0)]
    f_np = [sp.lambdify((z, mg.n, mg.c, mg.g1, mg.e2, d.W), e, "numpy") for e in exprs]
    rz = r.subs(sub)
    r_np = [sp.lambdify((z, mg.n, mg.c, mg.g1, mg.e2), e, "numpy") for e in (rz, sp.diff(rz, z), sp.diff(rz, z, 2))]
    return f_np, r_np


class Wormhole:
    def __init__(self, D, g1, e2, c=0.0):
        self.n = D - 2
        self.g1, self.e2, self.c = g1, e2, c
        if abs(g1**2 - e2) < 1e-14:            # Lam = 0: step off it (the formulas divide by Lam)
            self.e2 = e2 - 1e-12
        f_np, f_mp, bgf, dnu, nu_mp, r_mp = coefficient_functions()
        self.f_np, self.f_mp, self.bgf = f_np, f_mp, bgf
        p = (self.n, self.c, self.g1, self.e2)
        # throat: tan y = nu'(y)
        self.yth = float(mp.findroot(lambda yy: mp.tan(yy) - mp.re(dnu(yy, *p)), mp.atan(g1)))
        self.nu_th = float(mp.re(nu_mp(self.yth, *p)))
        self.a = float(mp.re(r_mp(self.yth, *p)))
        self.nu_end = [float(mp.re(nu_mp(s * (mp.pi / 2 - mp.mpf("1e-25")), *p))) for s in (-1, 1)]

    def coeffs(self, yy, W):
        p = (self.n, self.c, self.g1, complex(self.e2))   # complex: Lam = sqrt(g1^2 - e2) is imaginary when e2 > g1^2
        return [np.real(np.asarray(f(yy + 0j, *p, W + 0j))) for f in self.f_np]

    def bg(self, yy):
        p = (self.n, self.c, self.g1, complex(self.e2))
        return [np.real(np.asarray(f(yy + 0j, *p))) for f in self.bgf]

    def coeffs_z(self, zz, W, side):
        p = (self.n, self.c, self.g1, complex(self.e2))
        return [np.real(np.asarray(f(zz + 0j, *p, W + 0j))) for f in end_functions(side)[0]]

    def r_z(self, zz, side, k=0):
        """r (k = 0), dr/dz (k = 1) or d^2r/dz^2 (k = 2) at z = pi/2 - side y"""
        p = (self.n, self.c, self.g1, complex(self.e2))
        return np.real(np.asarray(end_functions(side)[1][k](zz + 0j, *p)))

    def frobenius(self, W, K=14, rho=None):
        """Taylor coefficients of P = z A1/A2 and Q = z^2 A0/A2 at z = 0; the two series solutions"""
        p = (self.n, self.c, self.g1, self.e2)
        rho = rho or 0.05 * min(1.0, np.pi / 2 - abs(self.yth))
        M = 64
        z = [rho * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M)]
        Pv, Qv = [], []
        for zz in z:
            A2, A1, A0 = (f(self.yth + zz, *p, W) for f in self.f_mp)
            Pv.append(zz * A1 / A2)
            Qv.append(zz * zz * A0 / A2)
        Pk = [mp.fsum(Pv[j] * z[j] ** (-k) for j in range(M)) / M for k in range(K + 1)]
        Qk = [mp.fsum(Qv[j] * z[j] ** (-k) for j in range(M)) / M for k in range(K + 1)]
        Pk = [mp.re(v) for v in Pk]
        Qk = [mp.re(v) for v in Qk]
        sols = {}
        for r0 in (0, 3):
            cks = [mp.mpf(1)]
            for k in range(1, K + 1 - r0):
                ind = (k + r0) * (k + r0 - 1) + Pk[0] * (k + r0) + Qk[0]
                rhs = -mp.fsum((Pk[j] * (k - j + r0) + Qk[j]) * cks[k - j] for j in range(1, k + 1))
                if abs(ind) < 1e-12:            # resonance (r0 = 0, k = 3): the log-free condition is rhs = 0
                    sols["resonance_rhs"] = float(rhs)
                    cks.append(mp.mpf(0))
                else:
                    cks.append(rhs / ind)
            sols[r0] = [float(v) for v in cks]
        sols["P0"], sols["Q0"] = float(Pk[0]), float(Qk[0])
        return sols

    def local(self, sols, r0, zz):
        cks = sols[r0]
        Y = sum(ck * zz ** (k + r0) for k, ck in enumerate(cks))
        Yp = sum(ck * (k + r0) * zz ** (k + r0 - 1) for k, ck in enumerate(cks) if k + r0 > 0)
        return Y, Yp

    def shoot(self, W, side, Y0, Yp0, z0, ymax_margin, rtol=1e-11, coord="z"):
        """Integrate the throat solution (Y0, Yp0 = dY/dy at y_th + side z0) to the end, in amplitude-phase form, and
        return psi, th at the end (no overflow however much Y grows).

        coord="z" (the default) integrates in t = log z, z = pi/2 - side y, down to z = ymax_margin, with the
        coefficients as functions of z (end_functions), and the phase of (Y, V), V = Y'/g with g = sqrt(1 + r'^2)
        (' = d/dy): Y = e^psi cos(th), V = e^psi sin(th). Near the throat g = 1 and V = Y'; near the end V = dY/dr up
        to sign, so th stays O(1) from the growing direction's multiples of pi. The equations are
            Y' = g V,   V' = -(A1/A2 + g'/g) V - A0/(A2 g) Y,
        and with d/dz = -side d/dy, a = Y_z e^-psi, b = V_z e^-psi: psi_z = cos a + sin b, th_z = cos b - sin a.
        coord="y" is the earlier integration in y with the phase of (Y, Y'). Its th goes to pi/2 mod pi, and cos(th),
        which is 1/(r' growth rate) there, keeps only 2e-16/cos(th) relative precision (2e-5 at D = 10): that error
        multiplies the A0 term of psi' and gives the right end's psi its noise (1e-7 at D = 8, 1e-4 at D = 10)."""
        if coord == "z":
            g0 = np.sqrt(1.0 + self.r_z(np.pi / 2 - side * self.yth - z0, side, 1) ** 2)
            s0 = [0.5 * np.log(Y0 * Y0 + (Yp0 / g0) ** 2), np.arctan2(Yp0 / g0, Y0)]

            def ft(t, s_):
                zz = np.exp(t)
                psi, th = s_
                A2, A1, A0 = self.coeffs_z(zz, W, side)
                rz, rzz = self.r_z(zz, side, 1), self.r_z(zz, side, 2)
                g2 = 1.0 + rz * rz
                g = np.sqrt(g2)
                c, sn = np.cos(th), np.sin(th)
                a = -side * g * sn
                b = (side * A1 / A2 - rz * rzz / g2) * sn + side * A0 * c / (A2 * g)
                return [zz * (c * a + sn * b), zz * (c * b - sn * a)]
            zstart = np.pi / 2 - side * self.yth - z0
            sol = solve_ivp(ft, (np.log(zstart), np.log(ymax_margin)), s0, method="DOP853", rtol=rtol, atol=rtol * 0.1)
            return sol.y[0, -1], sol.y[1, -1], side * (np.pi / 2 - ymax_margin)

        def f(yy, s_):
            psi, th = s_
            A2, A1, A0 = self.coeffs(yy, W)
            c, sn = np.cos(th), np.sin(th)
            F = -(A1 * sn + A0 * c) / A2
            return [sn * c + F * sn, -sn * sn + F * c]
        y0 = self.yth + side * z0
        yend = side * (np.pi / 2 - ymax_margin)
        s0 = [0.5 * np.log(Y0 * Y0 + Yp0 * Yp0), np.arctan2(Yp0, Y0)]
        sol = solve_ivp(f, (y0, yend), s0, method="DOP853", rtol=rtol, atol=rtol * 0.1)
        return sol.y[0, -1], sol.y[1, -1], yend

    def det(self, W2, rmax_factor=60.0, z0=None, X=30.0, form="amp", rtol=1e-11, psi_shift=0.0, coord="z"):
        """2 x 2 determinant of the growing parts at the two ends; its zero in W is the eigenvalue.

        At an end the growing exponential dominates (the decaying one is e^{-2 rmax_factor} smaller), so both throat
        solutions leave along the same direction in the (Y, V) plane of shoot: th = th_end + k pi, with th_end common
        to the end and k an integer per solution. The growing amplitude of a solution is then A = (-1)^k e^psi times a
        factor common to the end, so the determinant A_L0 A_R3 - A_L3 A_R0 is, up to a nonzero factor,
            s_L0 s_R3 e^{psi_L0 + psi_R3} - s_L3 s_R0 e^{psi_L3 + psi_R0},   s = (-1)^k, k = round((th - th_c) / pi),
        divided by the sum of the two exponentials (form="amp", the default), i.e. (P - Q)/(|P| + |Q|) for the two
        products P = A_L0 A_R3, Q = A_L3 A_R0: continuous wherever P and Q are not both 0, so each sign change is a
        zero. For a symmetric wormhole (gamma1 = 0) |P| = |Q| (the two throat solutions are even and odd), the form is
        +-1, and the zero (A_R0 or A_R3 = 0, where both vanish) is its jump, located to e^{-2 rmax_factor}. Otherwise
        it is analytic in W, in [-1, 1] at every D (divided by the larger exponential instead, it has a jump in its
        second derivative at the root, where the larger one changes over, of (e1 - e2)^2: 1e-8 across 1e-5 in What).
        th_c is a center that
        th_end cannot be pi/2 away from: with coord="z", V = Y'/g with g > 0, and Y'/Y has the sign of side on the
        growing direction (Y grows toward the end), so th_end is in (0, pi/2) at the right end and (pi/2, pi) at the
        left, and th_c = pi/2 - side pi/4; with coord="y", th_end is pi/2 to within 1/(Y'/Y), and th_c = pi/2.
        form="cos" is the earlier e^psi cos(th) form (Y at the ends), scaled by the larger exponential. Its zeros are
        the same; with coord="y" its cos(th_end) is 1e-4 to 1e-9 at D = 8 and falls with D, so that form is 1e-13 at
        D = 8, 3e-15 at D = 10 and at the precision of th from D = 12 (with coord="z" cos(th_end) is O(1)).
        rtol is the shooting's relative tolerance (psi ~ 1e2, so e^psi is good to about 1e2 rtol relative).
        psi_shift is added to psi_L0, for the controls of amp_check.py (a wrong amplitude must move the root).
        coord is the shooting's coordinate and phase (shoot): "z" (the default), or "y", the earlier one, whose right end
        amplitudes carry noise that grows fast with D (1e-7 in psi at D = 8, 1e-4 at D = 10).
        sols["th_spread"] is the larger difference, over the two ends, between the two solutions' th modulo pi: the amp
        form assumes it is 0 (both solutions on the growing direction), so it should stay near the shooting's error.
        sols["k"][(side, r0)] and sols["th"][(side, r0)] are each solution's k (its growing amplitude changes sign
        where k changes by 1) and end angle."""
        W = np.sqrt(W2)
        sols = self.frobenius(W)
        z0 = z0 or 0.02 * min(1.0, np.pi / 2 - abs(self.yth))
        rows = []
        for side in (-1, 1):
            Wloc = W * np.exp(-self.nu_end[0 if side < 0 else 1])
            target = max(rmax_factor / Wloc, self.a * (1 + X / (2 * (self.n - 1))))
            if coord == "z":
                # r(z) falls from infinity at z = 0 to the throat's radius: bisect in log z for r = target
                lo, hi = -700.0, np.log(np.pi / 2 - side * self.yth)
                for _ in range(200):
                    mid = 0.5 * (lo + hi)
                    if float(self.r_z(np.exp(mid), side)) > target:
                        lo = mid
                    else:
                        hi = mid
                margin = np.exp(lo)
            else:
                yl, yh = self.yth, side * (np.pi / 2 - 1e-14)
                for _ in range(200):
                    ym = 0.5 * (yl + yh)
                    if self.bg(ym)[1] < target:
                        yl = ym
                    else:
                        yh = ym
                margin = np.pi / 2 - abs(yh)
            row, ths = [], []
            th_c = np.pi / 2 - side * np.pi / 4 if coord == "z" else np.pi / 2
            for r0 in (0, 3):
                Y0, Yp0 = self.local(sols, r0, side * z0)
                psi, th, _ = self.shoot(W, side, Y0, Yp0, z0, margin, rtol=rtol, coord=coord)
                ths.append(th)
                sols.setdefault("k", {})[(side, r0)] = round((th - th_c) / np.pi)
                sols.setdefault("th", {})[(side, r0)] = th
                if form == "cos":
                    row.append((psi, np.cos(th)))
                else:
                    row.append((psi, -1.0 if round((th - th_c) / np.pi) % 2 else 1.0))
            rows.append(row)
            dth = (ths[0] - ths[1]) % np.pi
            sols["th_spread"] = max(sols.get("th_spread", 0.0), min(dth, np.pi - dth))
        (pL0, cL0), (pL3, cL3) = rows[0]
        (pR0, cR0), (pR3, cR3) = rows[1]
        e1, e2 = pL0 + psi_shift + pR3, pL3 + pR0
        m = np.logaddexp(e1, e2) if form == "amp" else max(e1, e2)
        sols["e1"], sols["e2"] = e1, e2
        return np.exp(e1 - m) * cL0 * cR3 - np.exp(e2 - m) * cL3 * cR0, sols

    def W2_of(self, What):
        """coordinate-time W^2 for a growth rate What = W_th a in throat proper time and throat radius units"""
        return (What * np.exp(self.nu_th) / self.a) ** 2

    def growth(self, lo, hi, N=24):
        """search What = W_th a in [lo, hi]"""
        grid = np.linspace(lo, hi, N)
        vals = [self.det(self.W2_of(h))[0] for h in grid]
        for (a, fa), (b, fb) in zip(zip(grid, vals), zip(grid[1:], vals[1:])):
            if np.sign(fa) != np.sign(fb):
                What = brentq(lambda h: self.det(self.W2_of(h))[0], a, b, xtol=1e-12)
                return self.W2_of(What), What
        return None, None


if __name__ == "__main__":
    D, g1, e2 = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    lo, hi = (float(sys.argv[4]), float(sys.argv[5])) if len(sys.argv) > 5 else (0.3, 3.0)
    w = Wormhole(D, g1, e2)
    print(f"D = {D}, gamma1 = {g1}, eps2 = {e2}: y_th = {w.yth:.6f}, a = {w.a:.6f}, nu_th = {w.nu_th:.6f}, "
          f"nu_ends = {w.nu_end}")
    sols = w.frobenius(1.0)
    print(f"  Frobenius: P0 = {sols['P0']:.6f}, Q0 = {sols['Q0']:.2e}, resonance rhs = {sols.get('resonance_rhs')}")
    W2, What = w.growth(lo, hi)
    print(f"  W^2 = {W2}, What = W_th a = {What}, T = 1/What = {1 / What if What else None}")

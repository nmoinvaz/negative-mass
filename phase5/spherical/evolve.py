"""
1+1 (spherical) nonlinear evolution of a lone sextic boson star with signed gravitational source.

Model (Shirokov arXiv:2608.24577, one sector): G_ab = 8 pi eps T_ab, eps = +1 canonical, -1 phantom;
field equation canonical, Box Phi = 2 V'(|Phi|^2) Phi, V(s) = s/2 - lam s^2/4 + mu s^3/6.

Gauge: polar-areal slicing  ds^2 = -alpha^2 dt^2 + a^2 dr^2 + r^2 dOmega^2.
Variables: complex Phi(t,r), Pi = (a/alpha) d_t Phi.
    d_t Phi = (alpha/a) Pi
    d_t Pi  = r^-2 d_r( r^2 (alpha/a) d_r Phi ) - alpha a 2V'(|Phi|^2) Phi
Metric from the constraints on every RK stage (no metric evolution):
    Hamiltonian:   m' = 4 pi r^2 eps [ K/2 (1 - 2m/r) + V ],  K = |Pi|^2 + |Phi_r|^2,  a^-2 = 1 - 2m/r
                   (m = Misner-Sharp mass, exact quasi-local mass in spherical symmetry)
    slicing:       (ln alpha)' = m/(r(r-2m)) + 4 pi r eps ( K/2 - a^2 V )
                   normalised so that alpha a = 1 at the outer boundary (t = asymptotic time up to
                   the matter outside R_out, which is ~0).
    momentum constraint (diagnostic only): d_t m = 4 pi eps r^2 (alpha/a^3) Re(Pi^* Phi_r)

Discretisation: cell-centred grid r_i = (i+1/2) h, 4th-order centred differences, parity ghosts at
the origin (Phi, Pi even), RK4 in time, 6th-order Kreiss-Oliger dissipation (eps_KO * h^5 term),
quadratic sponge on [r_sp, R].  The radial ODEs are integrated with 4th-order accuracy (RK4 / Simpson
with 4-point midpoint interpolation, parity extension at r = 0).
"""
import numpy as np
import numba as nb
import os, sys, time, json

LAM = 10240.0
MU = LAM**2 / 4.8
NG = 3  # ghost points each side

@nb.njit(cache=True, inline='always')
def V_(s):
    return 0.5 * s - 0.25 * LAM * s * s + MU * s * s * s / 6.0

@nb.njit(cache=True, inline='always')
def dV2_(s):
    return 1.0 - LAM * s + MU * s * s

@nb.njit(cache=True)
def fill_ghosts(f, n):
    # f has length n + 2*NG; interior f[NG:NG+n]; even parity at origin, zero outside
    for g in range(NG):
        f[NG - 1 - g] = f[NG + g]
        f[NG + n + g] = 0.0

@nb.njit(cache=True)
def metric(Phi, Pi, n, h, eps, m, lnal, Phr, Kc, Vc):
    """Given Phi, Pi (with ghosts filled) compute Phi_r, m(r), ln alpha(r) on interior points."""
    for i in range(n):
        j = i + NG
        Phr[i] = (Phi[j - 2] - 8.0 * Phi[j - 1] + 8.0 * Phi[j + 1] - Phi[j + 2]) / (12.0 * h)
        Kc[i] = abs(Pi[j])**2 + abs(Phr[i])**2
        Vc[i] = V_(abs(Phi[j])**2)
    # extended arrays for interpolation: even ghosts at origin, zero beyond R (fields vanish there)
    fourpi = 4.0 * np.pi
    # --- m: RK4 with midpoint interpolation ---
    # half step 0 -> r0
    K0, K1 = Kc[0], Kc[1]
    V0, V1 = Vc[0], Vc[1]
    Kq = 1.09375 * K0 - 0.09375 * K1
    Vq = 1.09375 * V0 - 0.09375 * V1
    r0 = 0.5 * h
    rq = 0.25 * h
    hh = 0.5 * h
    # f(r,m) = 4pi eps r^2 (K/2 (1-2m/r) + V) = A - B m, A = 4pi eps r^2 (K/2 + V), B = 4pi eps r K
    mm = 0.0
    k1 = 0.0
    Aq = fourpi * eps * rq * rq * (0.5 * Kq + Vq); Bq = fourpi * eps * rq * Kq
    k2 = Aq - Bq * (mm + 0.5 * hh * k1)
    k3 = Aq - Bq * (mm + 0.5 * hh * k2)
    A0 = fourpi * eps * r0 * r0 * (0.5 * K0 + V0); B0 = fourpi * eps * r0 * K0
    k4 = A0 - B0 * (mm + hh * k3)
    m[0] = mm + hh / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    for i in range(n - 1):
        ri = (i + 0.5) * h
        rh = ri + 0.5 * h
        rn = ri + h
        Km1 = Kc[i - 1] if i >= 1 else Kc[0]
        Vm1 = Vc[i - 1] if i >= 1 else Vc[0]
        Kp2 = Kc[i + 2] if i + 2 < n else 0.0
        Vp2 = Vc[i + 2] if i + 2 < n else 0.0
        Kh = (-Km1 + 9.0 * Kc[i] + 9.0 * Kc[i + 1] - Kp2) / 16.0
        Vh = (-Vm1 + 9.0 * Vc[i] + 9.0 * Vc[i + 1] - Vp2) / 16.0
        Ai = fourpi * eps * ri * ri * (0.5 * Kc[i] + Vc[i]); Bi = fourpi * eps * ri * Kc[i]
        Ah = fourpi * eps * rh * rh * (0.5 * Kh + Vh); Bh = fourpi * eps * rh * Kh
        An = fourpi * eps * rn * rn * (0.5 * Kc[i + 1] + Vc[i + 1]); Bn = fourpi * eps * rn * Kc[i + 1]
        mm = m[i]
        k1 = Ai - Bi * mm
        k2 = Ah - Bh * (mm + 0.5 * h * k1)
        k3 = Ah - Bh * (mm + 0.5 * h * k2)
        k4 = An - Bn * (mm + h * k3)
        m[i + 1] = mm + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    # --- ln alpha: Simpson with interpolated midpoints ---
    # g(r) = m/(r(r-2m)) + 4 pi eps r (K/2 - a^2 V),  a^2 = r/(r-2m)
    mq = rq**3 * (1.09375 * m[0] / r0**3 - 0.09375 * m[1] / (3 * r0)**3)
    gq = mq / (rq * (rq - 2 * mq)) + fourpi * eps * rq * (0.5 * Kq - rq / (rq - 2 * mq) * Vq)
    g0 = m[0] / (r0 * (r0 - 2 * m[0])) + fourpi * eps * r0 * (0.5 * K0 - r0 / (r0 - 2 * m[0]) * V0)
    lnal[0] = hh / 6.0 * (0.0 + 4.0 * gq + g0)
    gi = g0
    for i in range(n - 1):
        ri = (i + 0.5) * h
        rh = ri + 0.5 * h
        rn = ri + h
        Km1 = Kc[i - 1] if i >= 1 else Kc[0]
        Vm1 = Vc[i - 1] if i >= 1 else Vc[0]
        mm1 = m[i - 1] if i >= 1 else -m[0]
        Kp2 = Kc[i + 2] if i + 2 < n else 0.0
        Vp2 = Vc[i + 2] if i + 2 < n else 0.0
        mp2 = m[i + 2] if i + 2 < n else m[n - 1]
        Kh = (-Km1 + 9.0 * Kc[i] + 9.0 * Kc[i + 1] - Kp2) / 16.0
        Vh = (-Vm1 + 9.0 * Vc[i] + 9.0 * Vc[i + 1] - Vp2) / 16.0
        mh = (-mm1 + 9.0 * m[i] + 9.0 * m[i + 1] - mp2) / 16.0
        gh = mh / (rh * (rh - 2 * mh)) + fourpi * eps * rh * (0.5 * Kh - rh / (rh - 2 * mh) * Vh)
        gn = m[i + 1] / (rn * (rn - 2 * m[i + 1])) + fourpi * eps * rn * (0.5 * Kc[i + 1] - rn / (rn - 2 * m[i + 1]) * Vc[i + 1])
        lnal[i + 1] = lnal[i] + h / 6.0 * (gi + 4.0 * gh + gn)
        gi = gn
    # normalise alpha * a = 1 at R
    R = (n - 0.5) * h
    lna_R = -0.5 * np.log(1.0 - 2.0 * m[n - 1] / R)
    shift = -lna_R - lnal[n - 1]
    for i in range(n):
        lnal[i] += shift

@nb.njit(cache=True)
def rhs(Phi, Pi, n, h, eps, sko, sp, dPhi, dPi, m, lnal, Phr, Kc, Vc, F):
    fill_ghosts(Phi, n); fill_ghosts(Pi, n)
    metric(Phi, Pi, n, h, eps, m, lnal, Phr, Kc, Vc)
    # F = alpha/a with ghosts (even at origin, linear extrapolation outside)
    for i in range(n):
        r = (i + 0.5) * h
        a = 1.0 / np.sqrt(1.0 - 2.0 * m[i] / r)
        F[i + NG] = np.exp(lnal[i]) / a
    for g in range(NG):
        F[NG - 1 - g] = F[NG + g]
        F[NG + n + g] = F[NG + n - 1] + (g + 1) * (F[NG + n - 1] - F[NG + n - 2])
    c = sko / (64.0 * h)
    for i in range(n):
        j = i + NG
        r = (i + 0.5) * h
        a = 1.0 / np.sqrt(1.0 - 2.0 * m[i] / r)
        al = np.exp(lnal[i])
        Fp = (F[j - 2] - 8.0 * F[j - 1] + 8.0 * F[j + 1] - F[j + 2]) / (12.0 * h)
        Phrr = (-Phi[j - 2] + 16.0 * Phi[j - 1] - 30.0 * Phi[j] + 16.0 * Phi[j + 1] - Phi[j + 2]) / (12.0 * h * h)
        s = abs(Phi[j])**2
        dPhi[i] = F[j] * Pi[j]
        dPi[i] = F[j] * Phrr + Fp * Phr[i] + 2.0 * F[j] * Phr[i] / r - al * a * dV2_(s) * Phi[j]
        # Kreiss-Oliger
        dPhi[i] += c * (Phi[j - 3] - 6 * Phi[j - 2] + 15 * Phi[j - 1] - 20 * Phi[j] + 15 * Phi[j + 1] - 6 * Phi[j + 2] + Phi[j + 3])
        dPi[i] += c * (Pi[j - 3] - 6 * Pi[j - 2] + 15 * Pi[j - 1] - 20 * Pi[j] + 15 * Pi[j + 1] - 6 * Pi[j + 2] + Pi[j + 3])
        # sponge
        dPhi[i] -= sp[i] * Phi[j]
        dPi[i] -= sp[i] * Pi[j]

@nb.njit(cache=True)
def steps(Phi, Pi, n, h, dt, eps, sko, sp, nsteps, work_m, lnal, Phr, Kc, Vc, F):
    N = n + 2 * NG
    k1a = np.zeros(n, np.complex128); k1b = np.zeros(n, np.complex128)
    k2a = np.zeros(n, np.complex128); k2b = np.zeros(n, np.complex128)
    k3a = np.zeros(n, np.complex128); k3b = np.zeros(n, np.complex128)
    k4a = np.zeros(n, np.complex128); k4b = np.zeros(n, np.complex128)
    P1 = np.zeros(N, np.complex128); Q1 = np.zeros(N, np.complex128)
    for _ in range(nsteps):
        rhs(Phi, Pi, n, h, eps, sko, sp, k1a, k1b, work_m, lnal, Phr, Kc, Vc, F)
        for i in range(n):
            P1[i + NG] = Phi[i + NG] + 0.5 * dt * k1a[i]; Q1[i + NG] = Pi[i + NG] + 0.5 * dt * k1b[i]
        rhs(P1, Q1, n, h, eps, sko, sp, k2a, k2b, work_m, lnal, Phr, Kc, Vc, F)
        for i in range(n):
            P1[i + NG] = Phi[i + NG] + 0.5 * dt * k2a[i]; Q1[i + NG] = Pi[i + NG] + 0.5 * dt * k2b[i]
        rhs(P1, Q1, n, h, eps, sko, sp, k3a, k3b, work_m, lnal, Phr, Kc, Vc, F)
        for i in range(n):
            P1[i + NG] = Phi[i + NG] + dt * k3a[i]; Q1[i + NG] = Pi[i + NG] + dt * k3b[i]
        rhs(P1, Q1, n, h, eps, sko, sp, k4a, k4b, work_m, lnal, Phr, Kc, Vc, F)
        for i in range(n):
            Phi[i + NG] += dt / 6.0 * (k1a[i] + 2 * k2a[i] + 2 * k3a[i] + k4a[i])
            Pi[i + NG] += dt / 6.0 * (k1b[i] + 2 * k2b[i] + 2 * k3b[i] + k4b[i])
    fill_ghosts(Phi, n); fill_ghosts(Pi, n)
    metric(Phi, Pi, n, h, eps, work_m, lnal, Phr, Kc, Vc)

@nb.njit(cache=True)
def cum_even(q, n, h, p):
    """4th-order cumulative integral of r^p q(r) from 0 to r_i (q even about r=0)."""
    out = np.zeros(n)
    r0 = 0.5 * h; rq = 0.25 * h
    qq = 1.09375 * q[0] - 0.09375 * q[1]
    out[0] = 0.5 * h / 6.0 * (4.0 * rq**p * qq + r0**p * q[0])
    for i in range(n - 1):
        ri = (i + 0.5) * h
        qm1 = q[i - 1] if i >= 1 else q[0]
        qp2 = q[i + 2] if i + 2 < n else 0.0
        qh = (-qm1 + 9.0 * q[i] + 9.0 * q[i + 1] - qp2) / 16.0
        out[i + 1] = out[i] + h / 6.0 * (ri**p * q[i] + 4.0 * (ri + 0.5 * h)**p * qh + (ri + h)**p * q[i + 1])
    return out


# ------------------------------------------------------------------------------------------------
class Sim:
    def __init__(self, star, h, R, eps, kind='exact', delta=0.0, sko=0.2, rsp_frac=0.7, sp0=0.5,
                 cfl=0.25, seam_M=None):
        self.h, self.R, self.eps = h, R, eps
        n = int(round(R / h)); self.n = n
        self.r = (np.arange(n) + 0.5) * h
        self.dt = cfl * h
        self.sko = sko
        rsp = rsp_frac * R
        x = np.clip((self.r - rsp) / (R - rsp), 0, None)
        self.sp = sp0 * x * x
        Phi0, Pi0 = initial_data(star, self.r, kind, delta, seam_M)
        N = n + 2 * NG
        self.Phi = np.zeros(N, np.complex128); self.Pi = np.zeros(N, np.complex128)
        self.Phi[NG:NG + n] = Phi0; self.Pi[NG:NG + n] = Pi0
        self.m = np.zeros(n); self.lnal = np.zeros(n); self.Phr = np.zeros(n, np.complex128)
        self.Kc = np.zeros(n); self.Vc = np.zeros(n); self.F = np.zeros(N)
        self.t = 0.0
        steps(self.Phi, self.Pi, n, h, self.dt, float(eps), sko, self.sp, 0, self.m, self.lnal,
              self.Phr, self.Kc, self.Vc, self.F)

    def advance(self, nsteps):
        steps(self.Phi, self.Pi, self.n, self.h, self.dt, float(self.eps), self.sko, self.sp, nsteps,
              self.m, self.lnal, self.Phr, self.Kc, self.Vc, self.F)
        self.t += nsteps * self.dt

    def m_at(self, r):
        return np.interp(r, self.r, self.m)

    def diagnostics(self, robs):
        n = self.n
        Phi = self.Phi[NG:NG + n]; Pi = self.Pi[NG:NG + n]
        q = -4 * np.pi * np.imag(np.conj(Phi) * Pi)
        Ncum = cum_even(q, n, self.h, 2)
        # central value by 4th-order extrapolation to r=0 from cell centres (even function)
        # even function: quadratic in r^2 through r0, r1, r2 -> value at r = 0
        Phic = 1.171875 * Phi[0] - 0.1953125 * Phi[1] + 0.0234375 * Phi[2]
        cont = cum_even(np.sqrt(np.abs(Phi)**2 + np.abs(Pi)**2), n, self.h, 2)
        return dict(content=float(4 * np.pi * np.interp(20.0, self.r, cont)), t=self.t, Phic=Phic, amp=abs(Phic), phase=np.angle(Phic),
                    alc=np.exp(self.lnal[0]),
                    m=[float(self.m_at(r)) for r in robs], mR=float(self.m[-1]),
                    N=[float(np.interp(r, self.r, Ncum)) for r in robs],
                    maxphi=float(np.max(np.abs(Phi))))

    def mom_residual(self):
        """|d_t m - 4 pi eps r^2 (alpha/a^3) Re(Pi* Phi_r)| using a one-step finite difference."""
        n = self.n
        m0 = self.m.copy()
        Pi = self.Pi[NG:NG + n].copy(); Phr = self.Phr.copy()
        a = 1 / np.sqrt(1 - 2 * m0 / self.r); al = np.exp(self.lnal)
        pred0 = 4 * np.pi * self.eps * self.r**2 * al / a**3 * np.real(np.conj(Pi) * Phr)
        save = (self.Phi.copy(), self.Pi.copy(), self.t)
        self.advance(1)
        Pi1 = self.Pi[NG:NG + n]
        a1 = 1 / np.sqrt(1 - 2 * self.m / self.r); al1 = np.exp(self.lnal)
        pred1 = 4 * np.pi * self.eps * self.r**2 * al1 / a1**3 * np.real(np.conj(Pi1) * self.Phr)
        dm = (self.m - m0) / self.dt
        res = dm - 0.5 * (pred0 + pred1)
        self.Phi, self.Pi, self.t = save[0], save[1], save[2]
        self.advance(0)
        return res


def load_star(tag):
    d = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'stars', tag + '.npz'))
    return {k: d[k] for k in d.files}


def initial_data(st, r, kind, delta, seam_M):
    """kind: 'exact'  -> Phi = phi0, Pi = -i w (a/alpha) phi0 from the BVP star
             'seam'   -> paper-style exterior: for r > r_m (phi = 5% phi_c) phi is replaced by the
                         Yukawa tail phi_m (r_m/r) exp(-kappa (r-r_m)) and the lapse in Pi by the
                         Schwarzschild lapse of mass seam_M (enclosed mass at r_m), as in
                         GRTeclyn boson_star_ode.py
             'pert'   -> exact * (1 + delta exp(-r^2/4))  (real, i.e. amplitude perturbation)
             'scale'  -> exact * (1 + delta)              (uniform amplitude rescale)"""
    rs, phis, us, ms, nus = st['r'], st['phi'], st['u'], st['m'], st['nu']
    w = float(st['w']); Rs = float(st['R'])
    k = np.sqrt(1 - w * w)
    phi = np.interp(r, rs, phis)
    # high-accuracy: use cubic interpolation from the dense table
    from scipy.interpolate import CubicSpline
    lnphi = CubicSpline(rs, np.log(phis))(np.clip(r, rs[0], Rs))
    phi = np.exp(lnphi)
    inner = r < rs[0]
    if inner.any():
        phi[inner] = phis[0]
    far = r > Rs
    phi[far] = phis[-1] * (Rs / r[far]) * np.exp(-k * (r[far] - Rs))
    mm = CubicSpline(rs, ms)(np.clip(r, rs[0], Rs)); mm[far] = ms[-1]
    nu = CubicSpline(rs, nus)(np.clip(r, rs[0], Rs)); nu[far] = np.log(1 - 2 * ms[-1] / r[far])
    a = 1 / np.sqrt(1 - 2 * mm / r); al = np.exp(nu / 2)
    if kind == 'seam':
        # paper (GRTeclyn boson_star_ode.py): core profile up to r_m (first phi < 5% phi_c), then
        # analytic tail phi_m (r_m/r) e^{-kappa (r-r_m)} and Schwarzschild lapse with the enclosed mass
        # at r_m (M_tail = seam_M[0]); the core lapse is normalised independently, leaving a jump
        # alpha_core(r_m) - alpha_tail(r_m) = seam_M[1] (paper: -4.3e-5 canonical, +5.1e-5 phantom).
        # The lapse only enters the seed through Pi = -i w (a/alpha) phi; the metric is re-solved.
        M_tail, jump = seam_M
        phic = phis[0]
        im = np.where(phis < 0.05 * phic)[0][0]
        rm = rs[im]; phim = phis[im]
        ext = r > rm
        phi[ext] = phim * (rm / r[ext]) * np.exp(-k * (r[ext] - rm))
        al_tail = lambda x: np.sqrt(1 - 2 * M_tail / x)
        al_core_rm = np.exp(np.interp(rm, rs, nus) / 2)
        scale = (al_tail(rm) + jump) / al_core_rm
        al = al * scale
        al[ext] = al_tail(r[ext])
    Phi = phi.astype(np.complex128)
    Pi = -1j * w * a / al * phi
    if kind == 'pert':
        f = 1 + delta * np.exp(-r * r / 4)
        Phi *= f; Pi *= f
    elif kind == 'scale':
        Phi *= (1 + delta); Pi *= (1 + delta)
    return Phi, Pi

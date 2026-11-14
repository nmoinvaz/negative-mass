#!/usr/bin/env python3
"""Issue #75: the spectrum of T1 = D^-1 B1 for the tanh profile without a periodic box.

D~ = D + v B1 with B1 = (1/2) tau1 g0 and tau1 = -f'(r) cos(theta). In the coordinates (psi, chi) = (a + ib, a - ib)
of frank_loss.md section 2, D = diag(S, S) and B1 = diag(-M, +M) with S = sigma.grad and M = tau1/2, so the spectrum of
T1 is +-Sigma, Sigma = spec(S^-1 M), and a real zero mode of D~ at speed v is a real eigenvalue mu = 1/v of S^-1 M.
tau1 is axisymmetric about the x axis (taken as the z axis of the 2-spinors, a cyclic relabelling of the Pauli matrices),
so J_z = -i d_phi + sigma_z/2 is conserved. Each half-integer m = J_z is solved separately.

Method (no FFT, no box, no truncation of space). S^-1 M is compact on Hdot^1(R^3) with the inner product <S a, S b>.
Its Galerkin matrix on a finite space V_h in Hdot^1 is K^-1 N, K_ij = <S phi_i, S phi_j>, N_ij = <S phi_i, M phi_j>, and a
compact operator's Galerkin eigenvalues converge without spurious ones. V_h is spanned by f(r) Om_kappa(theta, phi), with
Om_kappa the spin-angular functions (sigma.L Om = -(1 + kappa) Om, J_z = m) and
    f = y^l (1 - y)^2 q(y),   y = r/(r + Lmap),   q a Jacobi polynomial,   l = kappa (kappa > 0) or -kappa - 1,
so f ~ r^l at 0 (regular) and f ~ r^-2 at infinity (the free decaying solution), on the whole of [0, infinity).
S (f Om_kappa) = (f' + (1 + kappa) f / r) (sigma.rhat Om_kappa), so K is block diagonal in kappa and every radial integral
is a polynomial or a smooth weight in y. The matrices are real, so eigenvalues are real or in conjugate pairs.

    python notes/warp/tanh_spectrum.py [--part=reduce,angular,radial,vc,planted,shift,tanh,real] [--sigmas=2,4]
                                       [--ms=0.5,1.5,...] [--procs=4] [--raw=DIR] [--quick]
    python notes/warp/tanh_spectrum.py --controls          every --mutate control, each must make its checks fail

Parts:
  reduce   the 4-to-2 reduction and T1's block form (symbolic, written here from the gammas of mass_bounds.md).
  angular  the spin-angular basis: orthonormal, sigma.rhat is a signed permutation kappa -> -kappa, cos(theta) against the
           Gaunt formula, and sigma.grad (f Om) against finite differences in Cartesian coordinates.
  radial   the radial operator against finite differences, and the y-integrals against quadrature in r.
  vc       v_c = 2 sqrt(S) / ||tau1||_3 by quadrature, against the values printed in mass_bounds.md.
  planted  tau1 = 6/(1 + r^2): mu = +-i/l_n, l_n = 1 + 2n/3, with multiplicities. Passes first or nothing else runs.
  shift    the same weight translated and dilated (axisymmetric, not spherical): the same spectrum.
  tanh     the tanh profile per m and sigma R: rho(T1) v_c, the leading eigenvalues, convergence in every parameter,
           the real candidates, the identity Re(v) int tau1 |psi|^2 = 0 and the residual of each eigenfunction.
  real     the converged window and the real eigenvalues in it (part of tanh's output).

Mutations (each must make its checks fail): g0block, ladder, dterm, half, lscale, sobolev, rot, trunc.
Exits 1 when any check fails.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import math
import resource
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import scipy.linalg as sla
import scipy.sparse as ssp
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.sparse.csgraph import connected_components
from scipy.special import eval_jacobi, eval_legendre, sph_harm_y

MUT = ""
FAILED = []
D = 2  # the radial basis decays as r^-D
VC58 = {2.0: 3.74490, 4.0: 2.47097, 8.0: 1.56243, 16.0: 0.98503}  # mass_bounds.md item 5
REF63 = {2.0: 0.4094, 4.0: 0.3674}  # frank_loss.md section 9, rho(T1) v_c
MUTATIONS = {  # name: (what it breaks, the run's options, the checks that must fail)
    "g0block": ("flips the sign of the lower block of g0", "--part=reduce", ["reduction"]),
    "ladder": ("flips the sign of the sigma.L ladder term", "--part=angular", ["angular_dirac", "angular_u"]),
    "dterm": ("drops (1 + kappa) f/r from the radial operator", "--part=radial", ["radial_op", "radial_int"]),
    "sobolev": ("uses 21/20 of the sharp Sobolev constant", "--part=vc", ["v_c"]),
    "half": ("drops the 1/2 in M = tau1/2", "--part=planted,tanh --quick --sigmas=2", ["planted", "ref63"]),
    "lscale": ("scales the l = 1 Legendre term of the shifted weight by 11/10", "--part=shift --quick", ["shift"]),
    "rot": ("solves with N times i, so the planted eigenvalues are real", "--part=planted --quick",
            ["real_identity"]),
    "trunc": ("keeps 2 angular levels and 8 radial functions", "--part=tanh --quick --sigmas=2",
              ["converged"]),
}


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    if not ok:
        FAILED.append(name)


def peak_gb():
    m = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss  # kB on Linux, bytes on macOS
    return m / 2**30 if sys.platform == "darwin" else m / 2**20


def sobolev():
    return 3 * (math.pi / 2) ** (4 / 3) * (21 / 20 if MUT == "sobolev" else 1.0)


# ------------------------------------------------------------------------------------------------ profile, v_c
def sech2(x):
    e = np.exp(-2 * np.abs(x))
    return 4 * e / (1 + e) ** 2


def tanh_fp(sig, R=1.0):
    """f'(r) for f = [tanh sig(r + R) - tanh sig(r - R)] / (2 tanh sig R)."""
    return lambda r: sig * (sech2(sig * (r + R)) - sech2(sig * (r - R))) / (2 * math.tanh(sig * R))


def v_c(sig):
    """kappa_K(1) = ||tau1||_3 / (2 sqrt S), with int |cos|^3 dOmega = pi; v_c = 1/kappa_K(1)."""
    fp = tanh_fp(sig)
    cut = 1.0 + 40 / sig
    I = sum(quad(lambda r: abs(fp(r)) ** 3 * r**2, a, b, limit=400, epsabs=0, epsrel=1e-13)[0]
            for a, b in ((0.0, 1.0), (1.0, cut)))
    return 2 * math.sqrt(sobolev()) / (math.pi * I) ** (1 / 3)


def part_vc():
    print(f"== vc: v_c by quadrature, S = {sobolev():.9f}")
    worst = 0.0
    for sig, ref in VC58.items():
        vc = v_c(sig)
        worst = max(worst, abs(vc - ref))
        print(f"  sigma R = {sig:4.1f}: v_c = {vc:.6f}  (mass_bounds.md {ref:.5f})")
    check("v_c", worst < 6e-6, f"equal to mass_bounds.md's printed values, max difference {worst:.1e}")


# ------------------------------------------------------------------------------------------------ reduce
def part_reduce():
    print("== reduce: D~ = D + v B1 in the coordinates (psi, chi) = (a + ib, a - ib) (symbolic)")
    d = sp.symbols("d1 d2 d3")  # commuting stand-ins for the partial derivatives
    tau, v = sp.symbols("tau1 v", real=True)
    I2, Z2 = sp.eye(2), sp.zeros(2, 2)
    s = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
    g0 = sp.I * sp.diag(I2, I2 if MUT == "g0block" else -I2)
    g = [sp.Matrix(sp.BlockMatrix([[Z2, sj], [sj, Z2]])) for sj in s]
    Dop = sum((g[j] * d[j] for j in range(3)), sp.zeros(4, 4))
    B1 = tau / 2 * g0
    W = sp.Matrix(sp.BlockMatrix([[I2, sp.I * I2], [I2, -sp.I * I2]]))  # (a, b) -> (psi, chi)
    Wo = sp.Matrix(sp.BlockMatrix([[sp.I * I2, I2], [-sp.I * I2, I2]]))  # (u, w) -> (w + iu, w - iu)
    Sd = sum((s[j] * d[j] for j in range(3)), sp.zeros(2, 2))
    want_D = sp.Matrix(sp.BlockMatrix([[Sd, Z2], [Z2, Sd]]))
    want_B = sp.Matrix(sp.BlockMatrix([[-tau / 2 * I2, Z2], [Z2, tau / 2 * I2]]))
    ok_D = sp.simplify(Wo * Dop * W.inv() - want_D) == sp.zeros(4, 4)
    ok_B = sp.simplify(Wo * B1 * W.inv() - want_B) == sp.zeros(4, 4)
    check("reduction", ok_D and ok_B, "W' D W^-1 = diag(S, S) and W' B1 W^-1 = diag(-M, +M), M = tau1/2")
    # T1 e = lam e  <=>  B1 e = lam D e, so psi: -M psi = lam S psi and chi: M chi = lam S chi
    print("  so spec(T1) = Sigma union -Sigma with Sigma = spec(S^-1 M), and a real zero mode at speed v is a real "
          "mu = 1/v")


# ------------------------------------------------------------------------------------------------ angular
class Chan:
    """Om_kappa with J_z = m: (a(theta) e^{i mu_u phi}, b(theta) e^{i (mu_u + 1) phi}), real amplitudes at phi = 0."""

    def __init__(self, m, kappa):
        self.m, self.kappa = m, kappa
        self.l = kappa if kappa > 0 else -kappa - 1
        self.mu_u, self.mu_d = int(round(m - 0.5)), int(round(m + 0.5))
        l, mu_u = self.l, self.mu_u
        present = [i for i, mu in enumerate((self.mu_u, self.mu_d)) if abs(mu) <= l]
        # sigma.L in the basis (Y_l^{mu_u} up, Y_l^{mu_u + 1} down): [[mu, s], [s, -(mu + 1)]], eigenvalues l, -l - 1
        M = np.array([[mu_u, 0.0], [0.0, -(mu_u + 1)]])
        if len(present) == 2:
            sq = math.sqrt(l * (l + 1) - mu_u * (mu_u + 1))
            M[0, 1] = M[1, 0] = -sq if MUT == "ladder" else sq
        w, V = np.linalg.eigh(M[np.ix_(present, present)])
        k = int(np.argmin(np.abs(w + 1 + kappa)))
        assert abs(w[k] + 1 + kappa) < 1e-12, (m, kappa, w)
        self.c = np.zeros(2)
        self.c[present] = V[:, k]

    def amp(self, theta):
        """(a, b) at phi = 0."""
        out = []
        for c, mu in zip(self.c, (self.mu_u, self.mu_d)):
            out.append(c * sph_harm_y(self.l, mu, theta, 0.0).real if c != 0 else 0 * theta)
        return out[0], out[1]


def kappas(m, dK):
    k0 = int(round(m + 0.5))
    ks = list(range(k0, k0 + dK + 1))
    return [-k for k in reversed(ks)] + ks


class Ang:
    """Angular matrices of the channels of one m on Gauss-Legendre nodes in u = cos(theta)."""

    def __init__(self, m, dK, nu=128):
        self.u, self.wu = leggauss(nu)
        self.s = np.sqrt(1 - self.u**2)
        th = np.arccos(self.u)
        self.ch = [Chan(m, k) for k in kappas(m, dK)]
        ab = [c.amp(th) for c in self.ch]
        self.a = np.array([x[0] for x in ab])
        self.b = np.array([x[1] for x in ab])
        self.ra = self.u * self.a + self.s * self.b  # sigma.rhat Om, up
        self.rb = self.s * self.a - self.u * self.b  # down

    def mats(self, g):
        """O_ij = <Om_i| g |Om_j>, A_ij = <sigma.rhat Om_i| g |Om_j>, for g(u) sampled on the nodes."""
        w = 2 * np.pi * self.wu * g
        O = (self.a * w) @ self.a.T + (self.b * w) @ self.b.T
        A = (self.ra * w) @ self.a.T + (self.rb * w) @ self.b.T
        return O, A


def part_angular():
    print("== angular: spin-angular basis of J_z = m")
    worst_on, worst_perm, worst_u, worst_dirac = 0.0, 0.0, 0.0, 0.0
    for m in (0.5, 1.5, 2.5):
        ang = Ang(m, 6)
        O1, _ = ang.mats(np.ones_like(ang.u))
        worst_on = max(worst_on, np.abs(O1 - np.eye(len(ang.ch))).max())
        R = ang.mats(np.ones_like(ang.u))[1]  # <sigma.rhat Om_i | Om_j> = <Om_i| sigma.rhat |Om_j>
        # signed permutation kappa -> -kappa
        P = np.zeros_like(R)
        for i, ci in enumerate(ang.ch):
            j = [jj for jj, cj in enumerate(ang.ch) if cj.kappa == -ci.kappa][0]
            P[i, j] = 1.0
        worst_perm = max(worst_perm, np.abs(np.abs(R) - P).max())
        # cos(theta) against the Gaunt formula <l + 1, mu| u |l, mu> = sqrt(((l + 1)^2 - mu^2) / ((2l + 1)(2l + 3)))
        Cq, _ = ang.mats(ang.u)

        def gaunt(l1, l2, mu):
            if abs(mu) > min(l1, l2) or abs(l1 - l2) != 1:
                return 0.0
            lo = min(l1, l2)
            return math.sqrt(((lo + 1) ** 2 - mu**2) / ((2 * lo + 1) * (2 * lo + 3)))

        Cf = np.zeros_like(Cq)
        for i, ci in enumerate(ang.ch):
            for j, cj in enumerate(ang.ch):
                Cf[i, j] = ci.c[0] * cj.c[0] * gaunt(ci.l, cj.l, ci.mu_u) + ci.c[1] * cj.c[1] * gaunt(ci.l, cj.l, ci.mu_d)
        worst_u = max(worst_u, np.abs(Cq - Cf).max())
        # sigma.grad (f Om_kappa) = (f' + (1 + kappa) f / r) sigma.rhat Om_kappa, f = r^l e^{-r^2/2}, by finite differences
        rng = np.random.default_rng(m.__hash__() % 1000)
        sg = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]

        def spinor(ch, X):
            r = np.linalg.norm(X)
            th, ph = math.acos(X[2] / r), math.atan2(X[1], X[0])
            a, b = ch.amp(np.array([th]))
            return r**ch.l * math.exp(-r * r / 2) * np.array([a[0] * np.exp(1j * ch.mu_u * ph),
                                                              b[0] * np.exp(1j * ch.mu_d * ph)])

        for ch in ang.ch[:3] + ang.ch[-3:]:
            X = rng.normal(size=3) * 0.7 + np.array([0.3, -0.2, 0.5])
            r = np.linalg.norm(X)
            h = 1e-3
            sgrad = np.zeros(2, complex)
            for j in range(3):
                e = np.zeros(3)
                e[j] = h
                d = (-spinor(ch, X + 2 * e) + 8 * spinor(ch, X + e) - 8 * spinor(ch, X - e) + spinor(ch, X - 2 * e)) / (12 * h)
                sgrad += sg[j] @ d
            f = r**ch.l * math.exp(-r * r / 2)
            fp = (ch.l / r - r) * f
            th, ph = math.acos(X[2] / r), math.atan2(X[1], X[0])
            a, b = ch.amp(np.array([th]))
            om = np.array([a[0] * np.exp(1j * ch.mu_u * ph), b[0] * np.exp(1j * ch.mu_d * ph)])
            rh = X / r
            srh = sum(sg[j] * rh[j] for j in range(3))
            want = (fp + (1 + ch.kappa) * f / r) * (srh @ om)
            worst_dirac = max(worst_dirac, np.abs(sgrad - want).max() / (np.abs(want).max() + 1e-300))
    check("angular", worst_on < 1e-12 and worst_perm < 1e-12,
          f"Om orthonormal (worst {worst_on:.1e}), |<Om_i|sigma.rhat|Om_j>| = delta_(i, -j) (worst {worst_perm:.1e})")
    check("angular_u", worst_u < 1e-12, f"<Om_i|cos theta|Om_j> against the Gaunt formula, worst {worst_u:.1e}")
    check("angular_dirac", worst_dirac < 1e-8,
          f"sigma.grad (f Om_kappa) = (f' + (1 + kappa) f/r) sigma.rhat Om_kappa, relative worst {worst_dirac:.1e}")


# ------------------------------------------------------------------------------------------------ radial basis
def jac(n, a, b, y):
    """Jacobi polynomial P_n^(a,b)(2y - 1) and its derivative in y."""
    x = 2 * y - 1
    q = eval_jacobi(n, a, b, x)
    if n == 0:
        return q, 0 * q
    return q, (n + a + b + 1) * eval_jacobi(n - 1, a + 1, b + 1, x)


def rbasis(kap, p, y):
    """Q_n(y) = q_n, H_n(y) with D f_n = ((1 - y)^(D + 1) / L) y^e H_n, f_n = y^l (1 - y)^D q_n; also e and l."""
    l = kap if kap > 0 else -kap - 1
    e = l - 1 if kap > 0 else l
    Q = np.empty((p, y.size))
    H = np.empty((p, y.size))
    for n in range(p):
        q, dq = jac(n, 2 * D - 2, 2 * e + 2, y)
        Q[n] = q
        if kap > 0:
            c0 = l if MUT == "dterm" else 2 * l + 1  # (l + 1 + kappa)
            H[n] = (c0 - (l + D) * y) * q + y * (1 - y) * dq
        else:
            H[n] = -(l + D) * q + (1 - y) * dq
    return Q, H, e, l


def rad_K(Hi, ei, Hj, ej, y, wy, L):
    return L * (Hi * (wy * y ** (ei + ej + 2) * (1 - y) ** (2 * D - 2))) @ Hj.T


def rad_N(Hi, ei, Qj, lj, y, wy, L, w):
    return L**2 * (Hi * (wy * y ** (2 + ei + lj) * (1 - y) ** (2 * D - 3) * w)) @ Qj.T


def rad_G(Qi, li, Qj, lj, y, wy, L, w):
    return L**3 * (Qi * (wy * y ** (2 + li + lj) * (1 - y) ** (2 * D - 4) * w)) @ Qj.T


def gl_nodes(npan, n=16):
    x, w = leggauss(n)
    edges = np.linspace(0.0, 1.0, npan + 1)
    lo, hi = edges[:-1, None], edges[1:, None]
    return ((hi - lo) / 2 * x + (hi + lo) / 2).ravel(), ((hi - lo) / 2 * w).ravel()


def part_radial():
    print("== radial: the radial operator and the y-integrals against direct quadrature in r")
    L, p = 2.0, 5
    rng = np.random.default_rng(1)
    worst_op, worst_int = 0.0, 0.0
    fpw = tanh_fp(2.0)
    wfun = lambda r: -fpw(r) / 2  # noqa: E731
    for kap in (-3, -1, 1, 2, 4):
        Q, H, e, l = rbasis(kap, p, np.array([0.5]))
        c = rng.normal(size=p)
        y_of = lambda r: r / (r + L)  # noqa: E731

        def f(r, kap=kap, c=c, l=l):
            y = y_of(r)
            q = sum(c[n] * jac(n, 2 * D - 2, 2 * (l - 1 if kap > 0 else l) + 2, y)[0] for n in range(p))
            return y**l * (1 - y) ** D * q

        for r in (0.3, 1.1, 2.9, 7.0):
            h = 1e-3 * r
            fr = f(r)
            fprime = (-f(r + 2 * h) + 8 * f(r + h) - 8 * f(r - h) + f(r - 2 * h)) / (12 * h)
            direct = fprime + (1 + kap) * fr / r
            y = np.array([y_of(r)])
            _, Hy, e_, l_ = rbasis(kap, p, y)
            formula = (1 - y[0]) ** (D + 1) / L * y[0] ** e_ * (c @ Hy[:, 0])
            worst_op = max(worst_op, abs(direct - formula) / (abs(direct) + 1e-12))
    check("radial_op", worst_op < 1e-7, f"D f = ((1 - y)^(D+1)/L) y^e H against finite differences, worst {worst_op:.1e}")
    # the integrals K, N, G of pairs of channels against composite Gauss-Legendre in r, with D f by finite differences
    y, wy = gl_nodes(120)
    rr = L * y / (1 - y)
    xg, wg = leggauss(48)
    edges = [0.0, 0.25, 0.75, 1.5, 3.0, 6.0, 12.0, 30.0, 100.0, 400.0, 3000.0]
    rq = np.concatenate([(hi - lo) / 2 * xg + (hi + lo) / 2 for lo, hi in zip(edges[:-1], edges[1:])])
    wq = np.concatenate([(hi - lo) / 2 * wg for lo, hi in zip(edges[:-1], edges[1:])])

    def fvals(kap, n, r):
        yy = r / (r + L)
        q = rbasis_q(kap, n, yy)
        l = kap if kap > 0 else -kap - 1
        return yy**l * (1 - yy) ** D * q

    def dfvals(kap, n, r):
        h = 1e-3 * r
        d = (-fvals(kap, n, r + 2 * h) + 8 * fvals(kap, n, r + h) - 8 * fvals(kap, n, r - h)
             + fvals(kap, n, r - 2 * h)) / (12 * h)
        return d + (1 + kap) * fvals(kap, n, r) / r

    for ki, kj in ((-2, 3), (1, -1), (3, -4)):
        Qi, Hi, ei, li = rbasis(ki, 3, y)
        Qj, Hj, ej, lj = rbasis(kj, 3, y)
        Kmat = rad_K(Hi, ei, Hi, ei, y, wy, L)
        Nmat = rad_N(Hi, ei, Qj, lj, y, wy, L, wfun(rr))
        Gmat = rad_G(Qi, li, Qj, lj, y, wy, L, wfun(rr))
        got = {"K": [], "N": [], "G": []}
        want = {"K": [], "N": [], "G": []}
        for (a_, b_) in ((0, 0), (1, 2), (2, 1), (2, 2)):
            want["K"].append(np.sum(wq * rq**2 * dfvals(ki, a_, rq) * dfvals(ki, b_, rq)))
            want["N"].append(np.sum(wq * rq**2 * dfvals(ki, a_, rq) * wfun(rq) * fvals(kj, b_, rq)))
            want["G"].append(np.sum(wq * rq**2 * fvals(ki, a_, rq) * wfun(rq) * fvals(kj, b_, rq)))
            got["K"].append(Kmat[a_, b_])
            got["N"].append(Nmat[a_, b_])
            got["G"].append(Gmat[a_, b_])
        for k_ in "KNG":  # error against the largest entry of the same integral, since some entries vanish
            worst_int = max(worst_int, np.abs(np.array(got[k_]) - want[k_]).max() / np.abs(want[k_]).max())
    check("radial_int", worst_int < 1e-6, f"K, N, G integrals in y against quadrature in r, worst relative {worst_int:.1e}")


def rbasis_q(kap, n, y):
    """q_n(y) alone, for the direct evaluation of f."""
    l = kap if kap > 0 else -kap - 1
    e = l - 1 if kap > 0 else l
    return eval_jacobi(n, 2 * D - 2, 2 * e + 2, 2 * y - 1)


# ------------------------------------------------------------------------------------------------ weights, assembly
class Weight:
    """M = sum_t w_t(r) g_t(u) with the name and the data the checks need."""

    def __init__(self, kind, **kw):
        self.kind, self.kw = kind, kw

    def terms(self, r, ang):
        """[(w_t(r) on the radial nodes, g_t(u) on the angular nodes)]."""
        h = 2.0 if MUT == "half" else 1.0
        if self.kind == "tanh":
            return [(-h * tanh_fp(self.kw["sig"])(r) / 2, ang.u)]
        if self.kind == "fl":
            return [(h * 3 / (1 + r**2), np.ones_like(ang.u))]
        s, a, Lw = self.kw["s"], self.kw["a"], self.kw["Lw"]
        un, wn = leggauss(96)
        wv = 3 * s / (1 + s**2 * (r[:, None] ** 2 - 2 * a * r[:, None] * un[None, :] + a**2))  # tau/2, tau = 6 s/(1 + ..)
        out = []
        for l in range(Lw + 1):
            Pl = eval_legendre(l, un)
            c = (2 * l + 1) / 2 * (wv * (wn * Pl)).sum(axis=1)
            if MUT == "lscale" and l == 1:
                c = c * 1.1
            if np.abs(c).max() > 1e-17:
                out.append((h * c, eval_legendre(l, ang.u)))
        return out

    def absterm(self, r, ang):
        """|M| for the normalisation of the identity: only the tanh weight (one term, |w| |u|)."""
        if self.kind != "tanh":
            return None
        return [(np.abs(self.terms(r, ang)[0][0]), np.abs(ang.u))]

    def label(self):
        return self.kind + str(self.kw)


def assemble(wt, m, dK, p, Lmap, npan, diag=False):
    ang = Ang(m, dK)
    nch = len(ang.ch)
    y, wy = gl_nodes(npan)
    r = Lmap * y / (1 - y)
    bas = [rbasis(c.kappa, p, y) for c in ang.ch]
    terms = wt.terms(r, ang)
    n = nch * p
    K = np.zeros((n, n))
    N = np.zeros((n, n))
    pat = np.zeros((nch, nch), bool)
    for i in range(nch):
        Q, H, e, l = bas[i]
        K[i * p:(i + 1) * p, i * p:(i + 1) * p] = rad_K(H, e, H, e, y, wy, Lmap)
    for w, g in terms:
        O, A = ang.mats(g)
        for i in range(nch):
            for j in range(nch):
                if abs(A[i, j]) > 1e-13:
                    pat[i, j] = True
                    N[i * p:(i + 1) * p, j * p:(j + 1) * p] += A[i, j] * rad_N(bas[i][1], bas[i][2], bas[j][0],
                                                                                bas[j][3], y, wy, Lmap, w)
    out = {"ang": ang, "K": K, "N": N, "pat": pat, "p": p, "nch": nch}
    if diag:  # Gram matrices of M and of |M|, and of M^2, for the diagnostics of tanh eigenfunctions
        def gram(tl, power):
            Gm = np.zeros((n, n))
            for w, g in tl:
                O, _ = ang.mats(g * g if power == 2 else g)
                for i in range(nch):
                    for j in range(nch):
                        if abs(O[i, j]) > 1e-13:
                            ww = w * w if power == 2 else w
                            Gm[i * p:(i + 1) * p, j * p:(j + 1) * p] += O[i, j] * rad_G(
                                bas[i][0], bas[i][3], bas[j][0], bas[j][3], y, wy, Lmap, ww)
            return Gm
        out["G"] = gram(terms, 1)
        out["Gabs"] = gram(wt.absterm(r, ang), 1)
        if wt.kind == "tanh":
            out["W2"] = gram(terms, 2)
    return out


def solve_sector(wt, m, dK, p, Lmap, npan, diag=False, nsel=24):
    """All eigenvalues mu of N x = mu K x, block by block, with diagnostics of the nsel largest of each block."""
    A = assemble(wt, m, dK, p, Lmap, npan, diag)
    ang, K, N, pat, nch = A["ang"], A["K"], A["N"], A["pat"], A["nch"]
    ncomp, lab = connected_components(ssp.csr_matrix(pat | pat.T), directed=False)
    rec = []
    for cidx in range(ncomp):
        chs = [i for i in range(nch) if lab[i] == cidx]
        idx = np.concatenate([np.arange(i * p, (i + 1) * p) for i in chs])
        Kb, Nb = K[np.ix_(idx, idx)], (N[np.ix_(idx, idx)] * (1j if MUT == "rot" else 1.0))
        Lc = np.linalg.cholesky(Kb)
        C = sla.solve_triangular(Lc, Nb, lower=True)
        C = sla.solve_triangular(Lc, C.conj().T, lower=True).conj().T
        mu, Z = np.linalg.eig(C)
        order = np.argsort(-np.abs(mu))
        mu, Z = mu[order], Z[:, order]
        sel = list(range(min(nsel, len(mu))))
        d = {"mu": mu, "chs": [ang.ch[i].kappa for i in chs], "diag": []}
        if diag:
            for k in sel:
                x = sla.solve_triangular(Lc.conj().T, Z[:, k], lower=False)
                Kx = (x.conj() @ Kb @ x).real
                Gm = (x.conj() @ A["G"][np.ix_(idx, idx)] @ x).real
                Ga = (x.conj() @ A["Gabs"][np.ix_(idx, idx)] @ x).real
                qrel = Gm / Ga if Ga > 0 else float("nan")
                v = 1 / mu[k]
                res = float("nan")
                if "W2" in A:
                    W2 = (x.conj() @ A["W2"][np.ix_(idx, idx)] @ x).real
                    res = math.sqrt(max(0.0, abs(v) ** 2 * W2 / Kx - 1.0))
                d["diag"].append({"k": k, "qrel": qrel, "res": res})
        rec.append(d)
    return rec


def real_cands(mu, tol=1e-9):
    return mu[np.abs(mu.imag) < tol * np.abs(mu)]


# ------------------------------------------------------------------------------------------------ planted
def fl_expected(m, nmax):
    """mu = +-i/l_n, l_n = 1 + 2n/3, each sign with multiplicity n + 3/2 - m at fixed m > 0 (frank_loss.md section 9)."""
    out = []
    for n in range(nmax + 1):
        c = int(round(n + 1.5 - m))
        if c > 0:
            out += [1j / (1 + 2 * n / 3)] * c + [-1j / (1 + 2 * n / 3)] * c
    return np.array(sorted(out, key=lambda z: (-abs(z), z.imag)))


def compare_levels(muall, m, nmax):
    """Pair the largest eigenvalues with the expected ones, level by level, and return the worst distance."""
    exp = fl_expected(m, nmax)
    pool = list(muall[np.argsort(-np.abs(muall))][:len(exp)])
    worst = 0.0
    for e in exp:
        j = int(np.argmin([abs(z - e) for z in pool]))
        worst = max(worst, abs(pool[j] - e))
        pool.pop(j)
    return worst, np.array(sorted(muall, key=lambda z: -abs(z))[:len(exp)])


def run_planted(kind, cfg, kw, label):
    wt = Weight(kind, **kw)
    worst, worst_re, ev_rho = 0.0, 0.0, 0.0
    cands = 0
    qworst = 0.0
    for m in cfg["planted_ms"]:
        rec = solve_sector(wt, m, cfg["pdK"] if kind == "fl" else cfg["sdK"], cfg["pp"], 2.0, cfg["pnpan"], diag=False)
        muall = np.concatenate([d["mu"] for d in rec])
        err, got = compare_levels(muall, m, cfg["pn"])
        worst = max(worst, err)
        worst_re = max(worst_re, (np.abs(got.real) / np.abs(got)).max())
        ev_rho = max(ev_rho, np.abs(muall).max())
        cands += len(real_cands(got))
        print(f"  {label} m = {m}: top |mu| = {np.abs(muall).max():.10f}, worst error over levels n <= {cfg['pn']} "
              f"{err:.1e}, max |Re mu|/|mu| {worst_re:.1e}", flush=True)
    return worst, worst_re, ev_rho, cands


def part_planted(cfg):
    print("== planted: tau1 = 6/(1 + r^2), mu = +-i/l_n, l_n = 1 + 2n/3, the first levels of each m > 0")
    worst, worst_re, rho, cands = run_planted("fl", cfg, {}, "FL")
    check("planted", worst < cfg["ptol"] and abs(rho - 1) < cfg["ptol"],
          f"rho(T1) = {rho:.12f} (exact 1), worst error over the levels {worst:.1e}")
    check("real_identity", worst_re < 1e-8 and cands == 0,
          f"Re(v) int tau1 |psi|^2 = 0 with int tau1 |psi|^2 > 0 forces Re mu = 0: max |Re mu|/|mu| {worst_re:.1e}, "
          f"real candidates {cands}")


def part_shift(cfg):
    print("== shift: tau1 = 6 s/(1 + s^2 |x - a e_z|^2), the Mobius image of the planted weight, same spectrum")
    kw = {"s": 1.2, "a": 0.5, "Lw": cfg["Lw"]}
    worst, worst_re, rho, cands = run_planted("flshift", cfg, kw, f"s = {kw['s']}, a = {kw['a']}")
    check("shift", worst < cfg["stol"],
          f"worst error over the levels {worst:.1e}, rho = {rho:.10f}, max |Re mu|/|mu| {worst_re:.1e}")


# ------------------------------------------------------------------------------------------------ tanh
def sector_task(t):
    global MUT
    MUT = t["mut"]
    t0 = time.time()
    wt = Weight("tanh", sig=t["sig"])
    out = {}
    for tag in ("A", "B"):
        c = t[tag]
        rec = solve_sector(wt, t["m"], c["dK"], c["p"], c["L"], c["npan"], diag=(tag == "A"))
        out[tag] = rec
    return t["sig"], t["m"], out, time.time() - t0, peak_gb()


def match_stable(muA, muB, tol=1e-5, lo=1e-3):
    """Largest |mu| above which every eigenvalue of A matches one of B to relative tol (A sorted by |mu|)."""
    cut = 0.0
    for z in muA:
        if abs(z) < lo:
            break
        if np.abs(muB - z).min() > tol * abs(z):
            cut = abs(z)
            break
    return cut


def part_tanh(cfg, args):
    sigmas = [float(s_) for s_ in args.sigmas.split(",")]
    ms = [float(x) for x in (args.ms or cfg["ms"]).split(",")]
    tasks = []
    for sig in sigmas:
        for m in ms:
            tasks.append({"sig": sig, "m": m, "mut": MUT,
                          "A": {"dK": cfg["dK"], "p": cfg["p"], "L": 2.0, "npan": cfg["npan"]},
                          "B": {"dK": cfg["dK"] + cfg["dKB"], "p": cfg["p"] + cfg["pB"], "L": 3.0,
                                "npan": int(cfg["npan"] * 1.5)}})
    print(f"== tanh: sigmas {sigmas}, m = {ms}, discretisation A dK = {cfg['dK']}, p = {cfg['p']}, Lmap = 2, "
          f"npan = {cfg['npan']}; B dK + {cfg['dKB']}, p + {cfg['pB']}, Lmap = 3, npan x 1.5", flush=True)
    res = {}
    with ProcessPoolExecutor(max_workers=args.procs) as ex:
        for sig, m, out, dt, pk in ex.map(sector_task, tasks):
            res[(sig, m)] = out
            print(f"  sigma R = {sig:g}, m = {m:4.1f}: done in {dt:.0f} s, peak resident {pk:.2f} GiB", flush=True)
    rho_all = {}
    for sig in sigmas:
        vc = v_c(sig)
        print(f"-- sigma R = {sig:g}, v_c = {vc:.6f}")
        rho, rhoB = 0.0, 0.0
        print("   m    | leading eigenvalues of S^-1 M times v_c (a quadruple +-a +-ib is two conjugate pairs of "
              "blocks)    | |mu| v_c A / B          | stable above | eigenfunction diagnostics")
        for m in ms:
            A, B = res[(sig, m)]["A"], res[(sig, m)]["B"]
            muA = np.concatenate([d["mu"] for d in A]) * vc
            muB = np.concatenate([d["mu"] for d in B]) * vc
            muA = muA[np.argsort(-np.abs(muA))]
            muB = muB[np.argsort(-np.abs(muB))]
            rho = max(rho, abs(muA[0]))
            rhoB = max(rhoB, abs(muB[0]))
            cut = match_stable(muA, muB)
            lead = " ".join(f"{z.real:+.5f}{z.imag:+.5f}i" for z in muA[:4])
            dg = [x for d in A for x in d["diag"][:4]]
            qmax = max((abs(x["qrel"]) for x in dg if np.isfinite(x["qrel"])), default=float("nan"))
            rmax = max((x["res"] for x in dg if np.isfinite(x["res"])), default=float("nan"))
            print(f"   {m:4.1f} | {lead} | {abs(muA[0]):.8f} / {abs(muB[0]):.8f} | {cut:11.4f}  | "
                  f"max |q_rel| {qmax:.1e}, max residual {rmax:.1e}")
            res[(sig, m)]["cut"] = cut
        rho_all[sig] = (rho, rhoB)
        print(f"   rho(T1) v_c = {rho:.8f} (A), {rhoB:.8f} (B); sqrt(3) rho v_c = {math.sqrt(3) * rho:.4f} "
              f"(FL: <= 1), the quadruple is +-Re +-Im i")
    for sig in sigmas:
        rho, rhoB = rho_all[sig]
        check("converged", abs(rho - rhoB) < 1e-6 * rho,
              f"sigma R = {sig:g}: rho(T1) v_c = {rho:.8f}, A and B differ by {abs(rho - rhoB):.1e}")
        if sig in REF63:
            ref = REF63[sig]
            rel = abs(rho - ref) / ref
            check("ref63", rel < 1e-3,
                  f"sigma R = {sig:g}: {rho:.5f} against frank_loss.md {ref} (relative difference {rel:.1e}; "
                  f"three-digit roundings {rho:.3f} and {ref:.3f})")
    # the real eigenvalues, and the identity Re(v) q = 0 on every one
    print("== real: real eigenvalues of S^-1 M (real zero-mode speeds v = 1/mu) and the identity")
    for sig in sigmas:
        vc = v_c(sig)
        for m in ms:
            r = res[(sig, m)]
            muA = np.concatenate([d["mu"] for d in r["A"]]) * vc
            muB = np.concatenate([d["mu"] for d in r["B"]]) * vc
            cut = r["cut"]
            winA, winB = muA[np.abs(muA) >= max(cut, 1e-3)], muB[np.abs(muB) >= max(cut, 1e-3)]
            realA, realB = real_cands(winA), real_cands(winB)
            print(f"   sigma R = {sig:g}, m = {m:4.1f}: eigenvalues with |mu| v_c >= {max(cut, 1e-3):.4f}: "
                  f"{len(winA)} (A), {len(winB)} (B); real: {len(realA)} (A), {len(realB)} (B); "
                  f"all real below that, any size: {len(real_cands(muA[np.abs(muA) > 1e-8]))} (A), "
                  f"{len(real_cands(muB[np.abs(muB) > 1e-8]))} (B)")
            qs = [abs(x["qrel"]) * abs(d["mu"][x["k"]].real) / abs(d["mu"][x["k"]])
                  for d in r["A"] for x in d["diag"] if np.isfinite(x["qrel"])]
            check("real_identity", len(realA) == 0 and len(realB) == 0 and (max(qs, default=0.0) < 1e-8),
                  f"sigma R = {sig:g}, m = {m}: real candidates {len(realA)} (A) {len(realB)} (B); "
                  f"max |Re mu|/|mu| |q_rel| over the leading eigenfunctions {max(qs, default=0.0):.1e}")
    return res


# ------------------------------------------------------------------------------------------------ controls
def run_controls():
    print("== controls: each --mutate makes its target checks fail; the same run without it passes them")
    bad, base = [], {}
    for name, (what, opts, targets) in MUTATIONS.items():
        res = {}
        for mut in ("", name):
            if mut == "" and opts in base:
                res[mut] = base[opts]
                continue
            cmd = [sys.executable, __file__] + ([f"--mutate={mut}"] if mut else []) + opts.split()
            t0 = time.time()
            pr = subprocess.run(cmd, capture_output=True, text=True)
            fl = [ln.split("] ")[1].split(":")[0] for ln in pr.stdout.splitlines() if ln.strip().startswith("[FAIL]")]
            res[mut] = (pr.returncode, fl, time.time() - t0)
            if mut == "":
                base[opts] = res[mut]
        hit = lambda fl: [t_ for t_ in targets if t_ in fl]  # noqa: E731
        fired = res[name][0] == 1 and len(hit(res[name][1])) == len(targets) and not hit(res[""][1])
        print(f"  {name:8s} ({what}; {opts}): baseline exit {res[''][0]}, mutated exit {res[name][0]}, failed "
              f"{res[name][1]}, {res[name][2]:.0f} s -> {'fired' if fired else 'DID NOT FIRE'}", flush=True)
        if not fired:
            bad.append(name)
    check("controls", not bad, "every control fired" if not bad else f"not fired: {bad}")


# ------------------------------------------------------------------------------------------------ main
FULL = {"dK": 20, "p": 100, "npan": 160, "dKB": 6, "pB": 30, "ms": "0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5",
        "planted_ms": [0.5, 1.5, 2.5], "pdK": 12, "sdK": 22, "pp": 60, "pnpan": 120, "pn": 4, "ptol": 1e-9,
        "stol": 1e-9, "Lw": 26}
QUICK = {"dK": 8, "p": 40, "npan": 60, "dKB": 4, "pB": 12, "ms": "0.5,1.5",
         "planted_ms": [0.5], "pdK": 6, "sdK": 14, "pp": 30, "pnpan": 60, "pn": 2, "ptol": 1e-6,
         "stol": 1e-5, "Lw": 16}


def main():
    global MUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="reduce,angular,radial,vc,planted,shift,tanh")
    ap.add_argument("--mutate", default="")
    ap.add_argument("--sigmas", default="2,4")
    ap.add_argument("--ms", default="")
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--controls", action="store_true")
    args = ap.parse_args()
    MUT = args.mutate
    if MUT and MUT not in MUTATIONS:
        sys.exit(f"unknown mutation {MUT}")
    cfg = dict(QUICK if args.quick else FULL)
    if MUT == "trunc":
        cfg.update({"dK": 2, "p": 8})
    parts = args.part.split(",")
    t0 = time.time()
    if args.controls:
        run_controls()
    else:
        order = [("reduce", part_reduce), ("angular", part_angular), ("radial", part_radial), ("vc", part_vc),
                 ("planted", lambda: part_planted(cfg)), ("shift", lambda: part_shift(cfg))]
        for name, fn in order:
            if name in parts:
                fn()
        if "tanh" in parts:
            if any(f in FAILED for f in ("planted", "shift", "angular", "radial_op", "radial_int")):
                print("== tanh skipped: a planted-truth or basis check failed, so real data are not computed")
            else:
                part_tanh(cfg, args)
    print(f"peak resident memory {peak_gb():.2f} GiB, wall {time.time() - t0:.0f} s")
    print("FAILED: " + ", ".join(FAILED) if FAILED else "all checks passed")
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()

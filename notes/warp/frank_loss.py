#!/usr/bin/env python3
"""Issue #63: the Witten-spinor existence threshold of mass_bounds.md item 5, sharpened with Frank-Loss.

Frank and Loss, arXiv:2201.03610 [FL22], prove that a nontrivial psi in Hdot^1 and L^3 with |sigma.grad psi| <= |A| |psi|
needs ||A||_3^2 >= 3 S (their Theorem 1, Remark 4 and Section 3). This script checks the steps that carry that bound to
the flat-slice Witten operator D~ = D + (1/2)(tr K) g0 of mass_bounds.md, and tests the result numerically.

    python notes/warp/frank_loss.py [--part=reduce,fl,vc,planted,krylov,scan] [--sigmas=2,4] [--Ns=64,96,128]
                                    [--vs=...] [--L=...]
    python notes/warp/frank_loss.py --controls          every --mutate control, each must make its check fail

Parts:
  reduce   the 4-to-2 reduction: D~ eps = 0 with eps = (a, b) iff sigma.grad psi = +(1/2) tr K psi and
           sigma.grad chi = -(1/2) tr K chi, psi = a + ib, chi = a - ib (symbolic); FL's Penrose formula for 2x2 and
           4x4 gammas (symbolic).
  fl       FL's constants: S_3, the Loss-Yau zero mode with ||A||_3^2 = 3 S, and the scalar extremal pair of FL's
           appendix written as a zero mode of D~ at imaginary v (symbolic).
  vc       v_c (kappa_K = 1) and sqrt(3) v_c for the tanh profile at sigma R = 2, 4, 8, 16.
  planted  the FFT eigen-solver on the planted extremal tr K = 6/(1 + r^2), whose spectral radius is 1 = 1/(sqrt3 v_c).
  krylov   for the tanh profile: the top eigenvalues of T1 = D^-1 B1, its norm, then GMRES solves of the Witten equation
           (I + v T1) psi = -v T1 eps0 per v, with residuals, the smallest singular value of D~ D^-1 and the identity.
  scan     at N = 64: the top 24 eigenvalues of T1 and s_min(D~ D^-1) for v/v_c = 0.25, 0.5, ..., 4.

Mutations (each must make its check fail): g0block, gamma, ly, sobolev, half, vscale, Jsign.
Exits 1 when any check fails. Peak memory at N = 128 is a few GB; run under nice.
"""
import argparse
import math
import resource
import subprocess
import sys
import time

import numpy as np
import scipy.fft as sfft
import sympy as sp
from scipy.sparse.linalg import LinearOperator, eigs, eigsh, gmres

ap = argparse.ArgumentParser()
ap.add_argument("--part", default="reduce,fl,vc,planted,krylov")
ap.add_argument("--mutate", default="")
ap.add_argument("--sigmas", default="2,4")
ap.add_argument("--Ns", default="64,96,128")
ap.add_argument("--vs", default="")
ap.add_argument("--L", type=float, default=0.0, help="box half-width; 0 uses #58's (8 at sigma R = 2, 5 at 4)")
ap.add_argument("--controls", action="store_true")
args = ap.parse_args()
MUT = args.mutate
PARTS = args.part.split(",")
FAILED = []
KRY = "--part=krylov --sigmas=2 --Ns=96 --vs=4.5"
MUTATIONS = {  # name: (what it breaks, the run's options, the checks that must fail)
    "g0block": ("flips the sign of the lower block of g0", "--part=reduce", ["reduction", "divxi"]),
    "gamma": ("scales gamma^1 and sigma_1 by 11/10", "--part=reduce", ["clifford", "penrose"]),
    "ly": ("scales the Loss-Yau vector potential by 11/10", "--part=fl,vc", ["loss_yau"]),
    "sobolev": ("uses 21/20 of the sharp Sobolev constant", "--part=fl,vc", ["S_3", "loss_yau_norm", "v_c"]),
    "half": ("drops the 1/2 in B1 = (1/2) tr K g0 in the FFT operator", "--part=planted,krylov --sigmas=2 --Ns=96 "
             "--vs=4.5", ["planted", "fl_bound"]),
    "vscale": ("solves with 1.001 v and checks the residual with v", KRY, ["gmres"]),
    "Jsign": ("flips the sign of J.xi in the Witten identity", KRY, ["identity"]),
}
if MUT and MUT not in MUTATIONS:
    sys.exit(f"unknown mutation {MUT}")


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    if not ok:
        FAILED.append(name)


def peak_gb():
    # ru_maxrss is in bytes on macOS, kB on Linux
    m = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return m / 2**30 if sys.platform == "darwin" else m / 2**20


S_SOB = 3 * (np.pi / 2) ** (4 / 3) * (21 / 20 if MUT == "sobolev" else 1)


# ------------------------------------------------------------------------------------------------ gammas
def paulis():
    s = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
    if MUT == "gamma":
        s[0] = s[0] * sp.Rational(11, 10)
    return s


def gammas_sym():
    """mass_bounds.md: g0 = i diag(1, 1, -1, -1), g^j = [[0, s_j], [s_j, 0]]."""
    I2, Z2 = sp.eye(2), sp.zeros(2, 2)
    g0 = sp.I * sp.diag(I2, I2 if MUT == "g0block" else -I2)
    g = [sp.Matrix(sp.BlockMatrix([[Z2, s], [s, Z2]])) for s in paulis()]
    return g0, g


# ------------------------------------------------------------------------------------------------ reduce
def part_reduce():
    print("== reduce: D~ eps = 0 on 4-spinors against two 2-spinor equations (symbolic)")
    x, y, z = X = sp.symbols("x y z", real=True)
    tau = sp.Function("tau", real=True)(x, y, z)  # tr K, any real function
    a = sp.Matrix([sp.Function(f"a{k}")(x, y, z) for k in range(2)])  # complex-valued
    b = sp.Matrix([sp.Function(f"b{k}")(x, y, z) for k in range(2)])
    g0, g = gammas_sym()
    s = paulis()
    eps = sp.Matrix.vstack(a, b)
    Dt = sum((g[j] * eps.diff(X[j]) for j in range(3)), sp.zeros(4, 1)) + tau / 2 * g0 * eps
    u, w = Dt[:2, :], Dt[2:, :]
    sdot = lambda p: sum((s[j] * p.diff(X[j]) for j in range(3)), sp.zeros(2, 1))  # noqa: E731
    psi, chi = a + sp.I * b, a - sp.I * b
    R1 = sdot(psi) - tau / 2 * psi
    R2 = sdot(chi) + tau / 2 * chi
    zero = lambda m: all(sp.expand(e) == 0 for e in m)  # noqa: E731
    # forward: (w + iu, w - iu) = (R1, R2); backward: (u, w) = ((R1 - R2)/(2i), (R1 + R2)/2). The map is invertible.
    fwd = zero(w + sp.I * u - R1) and zero(w - sp.I * u - R2)
    bwd = zero(u - (R1 - R2) / (2 * sp.I)) and zero(w - (R1 + R2) / 2)
    check("reduction", fwd and bwd, "D~eps = 0  <=>  sigma.grad psi = +tau psi/2 and sigma.grad chi = -tau chi/2")
    n2 = sp.expand((psi.H * psi + chi.H * chi)[0] - 2 * (a.H * a + b.H * b)[0])
    check("norms", n2 == 0, "|psi|^2 + |chi|^2 = 2 |eps|^2, so psi = chi = 0 iff eps = 0")
    # 2 Re<eps, g0 D~ eps> = div xi - tau |eps|^2 with xi^j = <eps, g0 g^j eps>: a kernel element has int tau |eps|^2 = 0
    xi = [(eps.H * g0 * g[j] * eps)[0] for j in range(3)]
    q = (eps.H * g0 * Dt)[0]
    dv = sp.expand(q + sp.conjugate(q) - sum(sp.diff(xi[j], X[j]) for j in range(3)) + tau * (eps.H * eps)[0])
    check("divxi", dv == 0, "2 Re<eps, g0 D~eps> = d_j xi^j - tr K |eps|^2, pointwise")
    # the Clifford relations FL assume, for sigma (their N = 2 in d = 3) and for the 4x4 gammas
    cl = all(sp.expand(m[i] * m[j] + m[j] * m[i] - 2 * int(i == j) * sp.eye(m[0].shape[0])) == sp.zeros(*m[0].shape)
             for m in (s, g) for i in range(3) for j in range(3))
    cl = cl and all(m.H == m for m in s + g)
    cl = cl and g0 * g0 == -sp.eye(4) and g0.H == -g0 and all(g0 * gj + gj * g0 == sp.zeros(4, 4) for gj in g)
    check("clifford", cl, "sigma_j and g^j Hermitian and anticommuting; g0^2 = -1, g0 anti-Hermitian, {g0, g^j} = 0")
    # FL eq. (penroseformula): sum_j |p_j - (1/3) g_j (g.p)|^2 = sum_j |p_j|^2 - (1/3) |g.p|^2 for any p_j = d_j phi
    ok = True
    for m in (s, g):
        n = m[0].shape[0]
        p = [sp.Matrix([sp.Symbol(f"p{j}{c}r", real=True) + sp.I * sp.Symbol(f"p{j}{c}i", real=True)
                        for c in range(n)]) for j in range(3)]
        gp = sum((m[j] * p[j] for j in range(3)), sp.zeros(n, 1))
        lhs = sum(((p[j] - m[j] * gp / 3).H * (p[j] - m[j] * gp / 3))[0] for j in range(3))
        rhs = sum((p[j].H * p[j])[0] for j in range(3)) - (gp.H * gp)[0] / 3
        ok = ok and sp.expand(lhs - rhs) == 0
    check("penrose", ok, "FL's pointwise Penrose formula holds for the 2x2 and the 4x4 gammas")


# ------------------------------------------------------------------------------------------------ fl
def part_fl():
    print("== fl: Frank-Loss constants and extremal pairs in d = 3 (symbolic)")
    x, y, z = X = sp.symbols("x y z", real=True)
    r = sp.Symbol("r", positive=True)
    S3 = sp.Rational(3, 4) * (2 * sp.pi**2) ** sp.Rational(2, 3)  # FL: S_d = d(d-2)/4 |S^d|^(2/d), |S^3| = 2 pi^2
    S_used = sp.nsimplify(S_SOB / (3 * (np.pi / 2) ** (4 / 3))) * 3 * (sp.pi / 2) ** sp.Rational(4, 3)
    check("S_3", sp.simplify(S3 - S_used) == 0, f"FL's S_3 = {float(S3):.12f}, item 5's S = {float(S_used):.12f}")
    s = paulis()
    r2 = x**2 + y**2 + z**2
    sx = s[0] * x + s[1] * y + s[2] * z
    # Loss-Yau: A = 3 (1 + r^2)^-2 ((1 - r^2) e1 + 2 x1 x + 2 Sigma x), Psi = (1 + r^2)^(-3/2) (1 + i s sigma.x) Psi0,
    # with (sigma_2 + i sigma_3) Psi0 = 0 and sigma_1 Psi0 = s Psi0 (FL Section 1, d = 3)
    Psi0 = sp.Matrix([1, 1]) / sp.sqrt(2)
    sgn = 1
    assert (s[1] + sp.I * s[2]) * Psi0 == sp.zeros(2, 1) or MUT == "gamma"
    Sig = sp.Matrix([[0, 0, 0], [0, 0, -1], [0, 1, 0]])
    Xv = sp.Matrix([x, y, z])
    num = (1 - r2) * sp.Matrix([1, 0, 0]) + 2 * x * Xv + 2 * Sig * Xv  # A = 3 num / (1 + r^2)^2
    sc = sp.Rational(11, 10) if MUT == "ly" else 1
    P = (sp.eye(2) + sp.I * sgn * sx) * Psi0  # Psi = (1 + r^2)^(-3/2) P
    # (1 + r^2)^(7/2) sigma.(-i grad - A) Psi, a polynomial
    dP = [(1 + r2) ** 2 * P.diff(X[j]) - 3 * X[j] * (1 + r2) * P for j in range(3)]
    eq = -sp.I * sum((s[j] * dP[j] for j in range(3)), sp.zeros(2, 1)) \
        - 3 * sc * sum((s[j] * num[j] for j in range(3)), sp.zeros(2, 2)) * P
    absA2 = sp.factor(sp.expand(9 * sc**2 * (num.T * num)[0]))
    check("loss_yau", all(sp.expand(e) == 0 for e in eq), f"sigma.(-i grad - A) Psi = 0, |A|^2 = {absA2} / (1 + r^2)^4")
    nA3 = 4 * sp.pi * sp.integrate(sp.Integer(27) * sc**3 * r**2 / (1 + r**2) ** 3, (r, 0, sp.oo))
    nA2 = sp.simplify(nA3 ** sp.Rational(2, 3))
    check("loss_yau_norm", sp.simplify(nA2 - 3 * S_used) == 0, f"||A||_3^2 = {nA2} = {float(nA2):.10f}, 3 S = "
          f"{float(3 * S_used):.10f}")
    # FL appendix: sigma.(-i grad) Psi~ = Lambda Psi~, Lambda = 3 s/(1 + r^2), Psi~ = (1 + r^2)^(-3/2)(1 + i s sigma.x) phi0
    # In D~'s reduced form sigma.grad psi = (v/2) tau1 psi this is tau1 = 6/(1 + r^2) at v = i s.
    phi0 = sp.Matrix([sp.Symbol("c0"), sp.Symbol("c1")])
    ok = True
    for sg in (1, -1):
        P = (sp.eye(2) + sp.I * sg * sx) * phi0
        dP = [(1 + r2) * P.diff(X[j]) - 3 * X[j] * P for j in range(3)]
        lhs = sum((s[j] * dP[j] for j in range(3)), sp.zeros(2, 1))  # (1 + r^2)^(5/2) sigma.grad Psi~
        v_ex = sp.I * sg
        ok = ok and all(sp.expand(e) == 0 for e in lhs - v_ex / 2 * 6 * P)
    tau_n3 = sp.simplify((4 * sp.pi * sp.integrate(216 * r**2 / (1 + r**2) ** 3, (r, 0, sp.oo))) ** sp.Rational(1, 3))
    kap1 = sp.simplify(tau_n3 / (2 * sp.sqrt(S_used)))
    check("fl_scalar_extremal", ok and sp.simplify(kap1 - sp.sqrt(3)) == 0,
          f"tau1 = 6/(1 + r^2) has zero modes at v = +-i; kappa_K(1) = {kap1}, so sqrt(3) v_c = |v| = 1")


# ------------------------------------------------------------------------------------------------ vc
def tanh_profile(sig, R=1.0):
    T = math.tanh(sig * R)

    def f(r):
        ta, tb = np.tanh(sig * (r + R)), np.tanh(sig * (r - R))
        sa, sb = 1 - ta**2, 1 - tb**2
        return (ta - tb) / (2 * T), sig * (sa - sb) / (2 * T), sig**2 * (-2 * sa * ta + 2 * sb * tb) / (2 * T)
    return f


def radial_nodes(sig, R=1.0, npan=80, nper=24, outer=40.0):
    a, b = max(0.0, R - 14 / sig), R + 14 / sig
    edges = list(np.linspace(0, a, 6)[:-1]) if a > 0 else []
    edges += list(np.linspace(a, b, npan + 1)) + list(np.linspace(b, R + outer / sig, 9)[1:])
    x, w = np.polynomial.legendre.leggauss(nper)
    rs, ws = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        rs.append((hi - lo) / 2 * x + (hi + lo) / 2)
        ws.append((hi - lo) / 2 * w)
    return np.concatenate(rs), np.concatenate(ws)


def kappa1(sig):
    """kappa_K at v = 1: ||tr K||_3 / (2 sqrt S), tr K = -v f'(r) cos(theta), int |cos|^3 dOmega = pi."""
    r, wr = radial_nodes(sig)
    F1 = tanh_profile(sig)(r)[1]
    return (np.pi * np.sum(wr * r**2 * np.abs(F1) ** 3)) ** (1 / 3) / (2 * math.sqrt(S_SOB))


VC58 = {2.0: 3.74490, 4.0: 2.47097, 8.0: 1.56243, 16.0: 0.98503}  # mass_bounds.md item 5, alcubierre_bounds_output.txt


def part_vc():
    print(f"== vc: v_c = 1/kappa_K(1) (contraction) and sqrt(3) v_c (Frank-Loss), S = {S_SOB:.9f}")
    print("  sigma R | v_c (#58)  | v_c here   | sqrt3 v_c  | thin wall sqrt3 v_c")
    worst = 0.0
    for sig in (2.0, 4.0, 8.0, 16.0):
        vc = 1 / kappa1(sig)
        thin = math.sqrt(3) * 2 * math.sqrt(S_SOB) * (15 / (2 * np.pi * sig**2)) ** (1 / 3)
        worst = max(worst, abs(vc - VC58[sig]))
        print(f"  {sig:7.1f} | {VC58[sig]:10.5f} | {vc:10.5f} | {math.sqrt(3) * vc:10.5f} | {thin:10.5f}")
    check("v_c", worst < 6e-6, f"equal to #58's printed values, max difference {worst:.1e}")


# ------------------------------------------------------------------------------------------------ FFT operators
class Grid:
    """Periodic box [-L, L)^3, N points per side; spinor fields are (4, N, N, N) complex arrays."""

    def __init__(self, L, N):
        self.L, self.N = L, N
        self.x1 = -L + 2 * L / N * np.arange(N)
        k1 = 2 * np.pi * np.fft.fftfreq(N, d=2 * L / N)
        self.k = (k1[:, None, None], k1[None, :, None], k1[None, None, :])
        k2 = self.k[0] ** 2 + self.k[1] ** 2 + self.k[2] ** 2
        k2[0, 0, 0] = 1.0
        self.ik2 = 1 / k2
        self.ik2[0, 0, 0] = 0.0
        self.dV = (2 * L / N) ** 3
        self.shape = (4, N, N, N)
        self.n = 4 * N**3

    def coords(self):
        return np.meshgrid(self.x1, self.x1, self.x1, indexing="ij")

    def fft(self, u):
        return sfft.fftn(u, axes=(1, 2, 3), workers=1)

    def ifft(self, u):
        return sfft.ifftn(u, axes=(1, 2, 3), workers=1)

    def sigk(self, u0, u1):
        """(sigma.k) (u0, u1)."""
        k1, k2, k3 = self.k
        return k3 * u0 + (k1 - 1j * k2) * u1, (k1 + 1j * k2) * u0 - k3 * u1

    def gk(self, h):
        """(g.k) h on a Fourier-space 4-spinor: g^j = [[0, s_j], [s_j, 0]]."""
        out = np.empty_like(h)
        out[0], out[1] = self.sigk(h[2], h[3])
        out[2], out[3] = self.sigk(h[0], h[1])
        return out

    def D(self, u):
        return self.ifft(1j * self.gk(self.fft(u)))

    def Dinv(self, u):
        h = self.gk(self.fft(u))
        h *= -1j * self.ik2  # (i g.k)^-1 = -i g.k/|k|^2; k = 0 set to zero
        return self.ifft(h)

    def dj(self, u, j):
        return self.ifft(1j * self.k[j] * self.fft(u))


G0D = np.array([1, 1, -1, -1])[:, None, None, None]  # g0 = i diag(1, 1, -1, -1)


def gam(j, e):
    """g^j e for j = 1, 2, 3 (x, y, z) and g0 for j = 0, on (4, ...) arrays."""
    if j == 0:
        return 1j * G0D * e
    a0, a1, b0, b1 = e
    if j == 1:
        return np.stack([b1, b0, a1, a0])
    if j == 2:
        return np.stack([-1j * b1, 1j * b0, -1j * a1, 1j * a0])
    return np.stack([b0, -b1, a0, -a1])


class Witten:
    """D~(v) = D + v B1 with B1 = (1/2) tau1 g0, tau1 = tr K at v = 1."""

    def __init__(self, grid, tau1):
        self.g, self.tau1 = grid, tau1
        self.half = 1.0 if MUT == "half" else 0.5
        self.napply = 0

    def B1(self, e):
        return (self.half * 1j) * self.tau1 * G0D * e

    def T1(self, e):  # D^-1 B1
        self.napply += 1
        return self.g.Dinv(self.B1(e))

    def M1(self, e):  # B1 D^-1, the L^2 form of T1
        self.napply += 1
        return self.B1(self.g.Dinv(e))

    def lin(self, fn):
        sh = self.g.shape
        return LinearOperator((self.g.n, self.g.n), matvec=lambda x: fn(x.reshape(sh)).ravel(), dtype=complex)


def top_eigs(W, k=10, ncv=24, tol=1e-7):
    lam = eigs(W.lin(W.T1), k=k, which="LM", ncv=ncv, tol=tol, return_eigenvectors=False)
    return lam[np.argsort(-np.abs(lam))]


def norm_T1(W):
    """||T1|| in Hdot^1 = ||B1 D^-1|| in L^2, the largest singular value."""
    # D^-1 and B1 are anti-Hermitian, so (B1 D^-1)^H = D^-1 B1 = T1 and the Gram operator is T1 M1
    H = W.lin(lambda e: W.T1(W.M1(e)))
    s2 = eigsh(H, k=1, which="LA", ncv=16, tol=1e-8, return_eigenvectors=False)
    return math.sqrt(abs(s2[0]))


def smin_Dt(W, v):
    """Smallest singular value of D~ D^-1 = I + v B1 D^-1 on mean-free L^2 fields (= inf ||D~psi|| / ||grad psi||)."""

    def P0(e):
        h = e - e.mean(axis=(1, 2, 3), keepdims=True)
        return h

    def H(e):
        c = e.mean(axis=(1, 2, 3), keepdims=True)
        h = e - c
        Ah = h + v * W.M1(h)  # A h
        out = Ah + v * W.T1(Ah)  # A^H = I + v D^-1 B1
        return P0(out) + c

    s2 = eigsh(W.lin(H), k=1, which="SA", ncv=20, tol=1e-9, maxiter=4000, return_eigenvectors=False)
    return math.sqrt(max(s2[0], 0.0))


def solve_gmres(W, v, eps0, rtol=1e-13, restart=20, maxcycles=40):
    g = W.g
    vs = v * (1.001 if MUT == "vscale" else 1.0)
    e0 = np.asarray(eps0, dtype=complex)[:, None, None, None] * np.ones((1, g.N, g.N, g.N))
    b = (-vs * W.T1(e0)).ravel()
    hist = []
    A = W.lin(lambda e: e + vs * W.T1(e))
    x, info = gmres(A, b, rtol=rtol, atol=0.0, restart=restart, maxiter=maxcycles,
                    callback=lambda rk: hist.append(rk), callback_type="pr_norm")
    psi = x.reshape(g.shape)
    r_pre = np.linalg.norm(b - A.matvec(x)) / np.linalg.norm(b)
    # the Witten equation itself, with D by spectral derivative and B1 by multiplication: D psi + v B1 (eps0 + psi)
    res = g.D(psi) + v * W.B1(e0 + psi)
    nb = math.sqrt(np.mean(np.abs(v * W.B1(e0)) ** 2))
    r_phys = math.sqrt(np.mean(np.abs(res) ** 2)) / nb
    m = res.mean(axis=(1, 2, 3), keepdims=True)  # the k = 0 mode, which a periodic box cannot set to zero
    r_mean = math.sqrt(np.mean(np.abs(m) ** 2)) / nb
    r_free = math.sqrt(np.mean(np.abs(res - m) ** 2)) / nb
    return e0 + psi, {"info": info, "its": len(hist), "r_pre": r_pre, "r_phys": r_phys, "r_mean": r_mean,
                      "r_free": r_free}


def picard(W, v, eps0, tol=1e-11, maxit=600):
    g = W.g
    e0 = np.asarray(eps0, dtype=complex)[:, None, None, None] * np.ones((1, g.N, g.N, g.N))
    psi = np.zeros_like(e0)
    for it in range(maxit):
        new = -v * W.T1(e0 + psi)
        d = math.sqrt(np.mean(np.abs(new - psi) ** 2))
        psi = new
        if d < tol or not np.isfinite(d) or d > 1e6:
            break
    return it + 1, d


def identity_terms(g, f, v, eps):
    """G = (1/4 pi) int |grad~ eps|^2 and int T(xi, n) = int (mu |eps|^2 - J.xi), mass_bounds.md item 6."""
    X = g.coords()
    r = np.sqrt(X[0] ** 2 + X[1] ** 2 + X[2] ** 2)
    rs = np.where(r == 0, 1.0, r)
    F0, F1, F2 = f(rs)
    n = [Xi / rs for Xi in X]
    del X
    F1r = np.where(r == 0, f(np.array([1e-8]))[2][0], F1 / rs)
    g0e = gam(0, eps)
    # grad~_i eps = d_i eps + (1/2) K_ij g^j g0 eps, K_ij = -(v/2) F1 (d_i1 n_j + d_j1 n_i)
    G2 = np.zeros(r.shape)
    gjg0 = [gam(j, g0e) for j in (1, 2, 3)]
    for i in range(3):
        Kg = -(v / 2) * F1 * n[i] * gjg0[0]
        if i == 0:
            Kg = Kg - (v / 2) * F1 * sum(n[j] * gjg0[j] for j in range(3))
        A = g.dj(eps, i) + 0.5 * Kg
        G2 += np.sum(np.abs(A) ** 2, axis=0)
        del A, Kg
    G = np.sum(G2) * g.dV / (4 * np.pi)
    trK = -v * F1 * n[0]
    nn = n[0] ** 2 + n[1] ** 2 + n[2] ** 2  # 1, except 0 at the grid point r = 0, where F1 and F2 are taken at r = 1
    KK = 2 * (v / 2) ** 2 * F1**2 * (n[0] ** 2 + nn)  # K_ij K_ij
    mu = (trK**2 - KK) / (16 * np.pi)
    hess0 = [F2 * n[0] * n[j] + F1r * ((1 if j == 0 else 0) - n[0] * n[j]) for j in range(3)]
    lap = F2 * nn + F1r * (3 - nn)
    J = [(v / (16 * np.pi)) * (hess0[j] - (lap if j == 0 else 0)) for j in range(3)]
    e2 = np.sum(np.abs(eps) ** 2, axis=0)
    xi = [np.real(np.sum(eps.conj() * gam(0, gam(j, eps)), axis=0)) for j in (1, 2, 3)]
    sJ = 1.0 if MUT == "Jsign" else -1.0
    T = mu * e2 + sJ * sum(J[j] * xi[j] for j in range(3))
    return G, np.sum(T) * g.dV, float(e2.min()), float(e2.max())


# ------------------------------------------------------------------------------------------------ planted
def part_planted():
    print("== planted: tau1 = 6/(1 + r^2) has zero modes at v = +-i (FL appendix), so rho(T1) = 1 = 1/(sqrt3 v_c)")
    L, N = 10.0, 64
    g = Grid(L, N)
    X = g.coords()
    tau1 = 6 / (1 + X[0] ** 2 + X[1] ** 2 + X[2] ** 2)
    del X
    W = Witten(g, tau1)
    t0 = time.time()
    lam = top_eigs(W, k=10, ncv=32)
    print(f"  L = {L}, N = {N}: top |lambda| = " + ", ".join(f"{abs(z_):.4f}" for z_ in lam[:10]))
    print("  top eigenvalues: " + ", ".join(f"{z_.real:+.4f}{z_.imag:+.4f}i" for z_ in lam[:8]))
    rho = abs(lam[0])
    im = max(abs(z_.real) / abs(z_) for z_ in lam[:8])
    spread = abs(lam[7]) / rho - 1
    print(f"  {W.napply} applications, {time.time() - t0:.0f} s, peak {peak_gb():.2f} GB")
    # By stereographic projection these are the Dirac eigenvalues +-(3/2 + k) on the round S^3, multiplicity
    # (k + 1)(k + 2) per sign, twice for the psi and chi sectors, so |lambda| = (3/2)/(3/2 + k) = 1, 3/5, ...
    check("planted", abs(rho - 1) < 0.03 and im < 0.03 and abs(spread) < 1e-3 and abs(abs(lam[8]) - 0.6) < 0.02,
          f"rho(T1) = {rho:.4f} (exact 1) eightfold, max |Re lambda|/|lambda| over the top 8 = {im:.1e} (exact 0), "
          f"ninth |lambda| = {abs(lam[8]):.4f} (exact 3/5)")


# ------------------------------------------------------------------------------------------------ krylov
BOX = {2.0: 8.0, 4.0: 5.0}  # #58's boxes
VS = {2.0: [1.0, 4.5, 5.5, 6.3, 13.0], 4.0: [0.5, 3.0, 3.6, 4.2, 8.5]}  # #58's row, three in (v_c, sqrt3 v_c), one
# near 2 sqrt3 v_c
REF58 = {(2.0, 1.0): (6.50304e-02, 0.9622, 1.0343), (4.0, 0.5): (2.84619e-02, 0.9840, 1.0085)}


def part_krylov():
    sigmas = [float(s_) for s_ in args.sigmas.split(",")]
    Ns = [int(n_) for n_ in args.Ns.split(",")]
    worst_res, worst_id, fl_ratio = 0.0, 0.0, 0.0
    for sig in sigmas:
        f = tanh_profile(sig)
        vc = 1 / kappa1(sig)
        L = args.L or BOX[sig]
        vs = [float(v_) for v_ in args.vs.split(",")] if args.vs else VS[sig]
        print(f"== krylov: sigma R = {sig}, L = {L}, v_c = {vc:.5f}, sqrt3 v_c = {math.sqrt(3) * vc:.5f}", flush=True)
        for N in Ns:
            g = Grid(L, N)
            X = g.coords()
            r = np.sqrt(X[0] ** 2 + X[1] ** 2 + X[2] ** 2)
            tau1 = -f(np.where(r == 0, 1.0, r))[1] * X[0] / np.where(r == 0, 1.0, r)  # tr K / v = -f' cos(theta)
            del X, r
            W = Witten(g, tau1)
            t0 = time.time()
            lam = top_eigs(W, k=8, ncv=20 if N >= 128 else 24)
            nT = norm_T1(W)
            rho = abs(lam[0])
            reals = [z_ for z_ in lam if abs(z_.imag) < 1e-6 * abs(z_)]
            print(f"  N = {N}: rho(T1) v_c = {rho * vc:.5f} (FL: <= 1/sqrt3 = {1 / math.sqrt(3):.5f}), "
                  f"||T1|| v_c = {nT * vc:.5f} (contraction: <= 1), {time.time() - t0:.0f} s", flush=True)
            print("    top eigenvalues of T1 times v_c: " + ", ".join(f"{z_.real * vc:+.4f}{z_.imag * vc:+.4f}i"
                                                                    for z_ in lam), flush=True)
            if reals:
                print(f"    largest real eigenvalue: first real zero-mode speed v = {1 / abs(reals[0]):.4f}")
            else:
                print(f"    no real eigenvalue among the top {len(lam)}: any real zero-mode speed is above "
                      f"{1 / abs(lam[-1]):.4f}")
            if N == Ns[-1]:
                fl_ratio = max(fl_ratio, rho * vc * math.sqrt(3))
            print("     v    | v/v_c | v rho(T1) | GMRES its | r_pre    | r_phys   | k=0 part | k!=0 part | "
                  "s_min(D~D^-1) | G           | rel. identity | |eps|^2 range   | Picard its", flush=True)
            for v in vs:
                t0 = time.time()
                eps, o = solve_gmres(W, v, [1, 0, 0, 0])
                sm = smin_Dt(W, v)
                G, iT, e2min, e2max = identity_terms(g, f, v, eps)
                del eps
                rel = abs(G + iT) / G
                pic = picard(W, v, [1, 0, 0, 0]) if N == Ns[0] else None
                if N == Ns[-1]:
                    worst_res = max(worst_res, o["r_phys"], o["r_pre"])
                    worst_id = max(worst_id, rel)
                pics = (f"{pic[0]} ({pic[1]:.0e})" if pic else "-")
                print(f"    {v:5.2f} | {v / vc:5.3f} | {v * rho:9.4f} | {o['its']:9d} | {o['r_pre']:.2e} | "
                      f"{o['r_phys']:.2e} | {o['r_mean']:.2e} | {o['r_free']:.2e}  | {sm:13.6f} | {G:.5e} | "
                      f"{rel:13.2e} | {e2min:.4f}-{e2max:.4f} | {pics}"
                      f"   [{time.time() - t0:.0f} s, peak {peak_gb():.2f} GB]", flush=True)
                if (sig, v) in REF58 and N == 128:
                    ref = REF58[(sig, v)]
                    check(f"#58 Picard value of G at ({sig:g}, {v:g}), N = 128", abs(G - ref[0]) < 1e-6 * ref[0] * 10,
                          f"{G:.6e} against {ref[0]:.5e}")
            del W, tau1
    check("gmres", worst_res < 1e-10, f"largest residual at the finest grid {worst_res:.1e} (target 1e-10)")
    check("identity", worst_id < 1e-6, f"worst relative |G + int T| at the finest grid {worst_id:.1e}")
    check("fl_bound", fl_ratio <= 1.0 + 1e-3, f"rho(T1) sqrt3 v_c = {fl_ratio:.4f} at the finest grid (FL: <= 1)")


def part_scan():
    """s_min(D~ D^-1) against v, and the top eigenvalues of T1, at one grid: where could a real zero mode sit?"""
    N = 64
    for sig in [float(s_) for s_ in args.sigmas.split(",")]:
        f = tanh_profile(sig)
        vc = 1 / kappa1(sig)
        g = Grid(BOX[sig], N)
        X = g.coords()
        r = np.sqrt(X[0] ** 2 + X[1] ** 2 + X[2] ** 2)
        tau1 = -f(np.where(r == 0, 1.0, r))[1] * X[0] / np.where(r == 0, 1.0, r)
        del X, r
        W = Witten(g, tau1)
        t0 = time.time()
        lam = top_eigs(W, k=24, ncv=56)
        reals = [z_ for z_ in lam if abs(z_.imag) < 1e-6 * abs(z_)]
        print(f"== scan: sigma R = {sig}, L = {BOX[sig]}, N = {N}, v_c = {vc:.5f}; top 24 eigenvalues of T1 times v_c "
              f"({time.time() - t0:.0f} s):")
        for i in range(0, 24, 6):
            print("    " + ", ".join(f"{z_.real * vc:+.4f}{z_.imag * vc:+.4f}i" for z_ in lam[i:i + 6]))
        if reals:
            print(f"  real eigenvalues among them: " + ", ".join(f"{z_.real * vc:+.4f}" for z_ in reals)
                  + f"; the first real zero-mode speed is v = {1 / abs(reals[0]):.4f} = {1 / abs(reals[0]) / vc:.3f} v_c")
        else:
            print(f"  no real eigenvalue among the top 24, so no real zero-mode speed below {1 / abs(lam[-1]):.4f} = "
                  f"{1 / abs(lam[-1]) / vc:.3f} v_c")
        print("  v/v_c | v       | s_min(D~ D^-1)")
        for t in np.arange(0.25, 4.01, 0.25):
            print(f"  {t:5.2f} | {t * vc:7.4f} | {smin_Dt(W, t * vc):.6f}", flush=True)
        print(f"  peak {peak_gb():.2f} GB")


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
            p = subprocess.run(cmd, capture_output=True, text=True)
            fl = [ln.split("] ")[1].split(":")[0] for ln in p.stdout.splitlines() if ln.strip().startswith("[FAIL]")]
            res[mut] = (p.returncode, fl, time.time() - t0)
            if mut == "":
                base[opts] = res[mut]
        hit = lambda fl: [t_ for t_ in targets if any(f_ == t_ or f_.startswith(t_) for f_ in fl)]  # noqa: E731
        fired = res[name][0] == 1 and len(hit(res[name][1])) == len(targets) and not hit(res[""][1])
        print(f"  {name:8s} ({what}; {opts}): baseline exit {res[''][0]}, mutated exit {res[name][0]}, failed "
              f"{res[name][1]}, {res[name][2]:.0f} s -> {'fired' if fired else 'DID NOT FIRE'}", flush=True)
        if not fired:
            bad.append(name)
    check("controls", not bad, "every control fired" if not bad else f"not fired: {bad}")


if args.controls:
    run_controls()
else:
    for name, fn in (("reduce", part_reduce), ("fl", part_fl), ("vc", part_vc), ("planted", part_planted),
                     ("krylov", part_krylov), ("scan", part_scan)):
        if name in PARTS:
            fn()
print(f"peak resident memory {peak_gb():.2f} GB")
print("FAILED: " + ", ".join(FAILED) if FAILED else "all checks passed")
sys.exit(1 if FAILED else 0)

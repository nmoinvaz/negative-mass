#!/usr/bin/env python3
"""Polar l = 1, 2 modes of the phantom black universes and the Ellis wormhole; radial modes of the M-AdS wormholes.  See polar.md.

    polar.py planted [--mutate=M]                       Ellis polar quasinormal modes (l = 2, 3, 4) against [ABCKY23] Table 1; the Ellis l = 1 problem;
                                                        the Ellis bound-state plant; the c -> 0+ limit of the M-AdS radial mode
    polar.py count ELL c [--res=R1|R2] [--out=F.json]   growing polar modes, l = ELL = 1 or 2, member c <= 0 (c = 0: Ellis), shooting, real W
    polar.py mads c [--res=R1|R2] [--out=F.json]        growing radial modes of the M-AdS wormhole c > 0 (Dirichlet at the AdS end)
    polar.py spectral DIR                               Chebyshev check of the M-AdS counts in DIR/mads_*_R2.json
    polar.py report DIR                                 tables, R1 against R2
    polar.py controls                                   runs every --mutate=M and asserts that its target fires

Units b = min R = 1.  Background and linearization: polar_derive.py (checked there); the coefficient table at the end of this file is its output.
Exit status 1 when a check fails (a mutated run is meant to exit 1).
"""
import glob
import json
import math
import os
import re
import subprocess
import sys
import time

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/nm_i218_numba")
import numpy as np
from numba import njit
from scipy.integrate import ode, solve_ivp
from scipy.optimize import brentq, minimize_scalar, root

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stability as st

# ------------------------------------------------------------------------------------------------ declared settings
W_LO, W_HI = 0.01, 8.0                                          # range searched, units 1/b
RES = {
    "R1": dict(rtol=1e-9, rmax0=150.0, rmaxW=30.0, nW=60, nWm=70, z0=0.04, nfrob=10, delta=0.3, ads=1.0),
    "R2": dict(rtol=1e-12, rmax0=300.0, rmaxW=45.0, nW=150, nWm=176, z0=0.03, nfrob=12, delta=0.2, ads=2.0),
}
W_HI_M = 24.0                                                   # M-AdS: the root grows like c^(0.7) (W = 5.2 at c = 3)
MEMBERS = st.MEMBERS + [0.0]                                    # #212's 20 black universes and the Ellis wormhole
ADS_C = [1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0]
PROBES = [0.01, 0.1, 1.0, 3.0, 8.0]                              # W at which R1 and R2 determinants are compared
DIP_RATIO = 1e-3                                                # an interior minimum of |D| below this fraction of its neighbors counts as a double root
MUT = set()


def mutating(name):
    return name in MUT

# ------------------------------------------------------------------------------------------------ coefficient table (polar_derive.py emit)
import math
from math import sqrt

def m4_entries(rho, lam, L, As, A1):
    x0 = As*rho
    x1 = -x0
    x2 = rho**2
    x3 = A1*x2 + A1
    x4 = x2 + 1
    x5 = 1/x4
    x6 = 1/As
    x7 = x5*x6
    x8 = (1/2)*L
    x9 = x7*((1/2)*A1*x4 + x1)
    x10 = math.sqrt(2)
    x11 = x10*x5
    x12 = 2*As
    x13 = x10*x6
    x14 = 4*lam
    x15 = A1**2
    return (1.0*(-x7*(x1 + x3)), 1.0*(x7*(As*x8 + lam*x4)), 1.0*(x9), 1.0*(-x11), 1.0*(x6), 1.0*(-A1*x6), 1.0*(x6), 1.0*(0), 1.0*(rho*x5), 1.0*(x5*x8), 1.0*(x9), 1.0*(x11), 1.0*(-1/4*x13*(A1*rho + L - x12)), 1.0*(-1/8*x13*(A1*L + rho*x14)), 1.0*((1/8)*x10*(2*A1*x0 - 4*As + L*x12 + x14*x2 + x14 - x15*x2 - x15)/As**2), 1.0*((1/2)*x7*(-4*x0 + x3)))

def m3_entries(rho, lam, As, A1):
    x0 = As*rho
    x1 = rho**2
    x2 = A1*x1 + A1
    x3 = x1 + 1
    x4 = 1/x3
    x5 = 1/As
    x6 = x4*x5
    x7 = (1/2)*x6
    x8 = math.sqrt(2)*x4
    x9 = -x8
    x10 = -4*x0
    x11 = A1*rho
    x12 = 2*x1
    x13 = As**2
    x14 = (1/2)*lam
    x15 = rho**4
    return (1.0*(-x7*(-2*x0 + x2)), 1.0*(-x6*(-x0 + x2)), 1.0*(x9), 1.0*(0), 1.0*(x7*(x10 + x2)), 1.0*(x9), 1.0*(-1/4*x5*x8*(A1*rho**3 - As*x12 + x11 + x12 + 2)), 1.0*(-x8*(-1/8*A1**2*(x15 + 1) + lam*x1 + (1/4)*x11*(As*x1 + As - x11) - 1/2*x13 + x14*x15 + x14)/x13), 1.0*(x7*(A1*x3 + x10)))

def g3_entries(rho, lam, As, A1):
    x0 = math.sqrt(As)
    x1 = rho**2
    x2 = x1 + 1
    x3 = As*rho
    x4 = A1*x2
    x5 = 2*x1
    x6 = x0/x2**(3/2)
    return (1.0*((1/2)*(A1*(-2*x3 + x4) - 4*lam*x2)/(x0*math.sqrt(x2))), 1.0*(x6*(A1*rho**3 + A1*rho - As*x5 + x5 + 2)), 1.0*(-math.sqrt(2)*x6*(-x3 + (1/2)*x4)))

def m4_entries_c(rho, lam, L, As, A1):
    x0 = As*rho
    x1 = -x0
    x2 = rho**2
    x3 = A1*x2 + A1
    x4 = x2 + 1
    x5 = 1/x4
    x6 = 1/As
    x7 = x5*x6
    x8 = (1/2)*L
    x9 = x7*((1/2)*A1*x4 + x1)
    x10 = math.sqrt(2)
    x11 = x10*x5
    x12 = 2*As
    x13 = x10*x6
    x14 = 4*lam
    x15 = A1**2
    return ((-x7*(x1 + x3))+0j, (x7*(As*x8 + lam*x4))+0j, (x9)+0j, (-x11)+0j, (x6)+0j, (-A1*x6)+0j, (x6)+0j, (0)+0j, (rho*x5)+0j, (x5*x8)+0j, (x9)+0j, (x11)+0j, (-1/4*x13*(A1*rho + L - x12))+0j, (-1/8*x13*(A1*L + rho*x14))+0j, ((1/8)*x10*(2*A1*x0 - 4*As + L*x12 + x14*x2 + x14 - x15*x2 - x15)/As**2)+0j, ((1/2)*x7*(-4*x0 + x3))+0j)

def m3_entries_c(rho, lam, As, A1):
    x0 = As*rho
    x1 = rho**2
    x2 = A1*x1 + A1
    x3 = x1 + 1
    x4 = 1/x3
    x5 = 1/As
    x6 = x4*x5
    x7 = (1/2)*x6
    x8 = math.sqrt(2)*x4
    x9 = -x8
    x10 = -4*x0
    x11 = A1*rho
    x12 = 2*x1
    x13 = As**2
    x14 = (1/2)*lam
    x15 = rho**4
    return ((-x7*(-2*x0 + x2))+0j, (-x6*(-x0 + x2))+0j, (x9)+0j, (0)+0j, (x7*(x10 + x2))+0j, (x9)+0j, (-1/4*x5*x8*(A1*rho**3 - As*x12 + x11 + x12 + 2))+0j, (-x8*(-1/8*A1**2*(x15 + 1) + lam*x1 + (1/4)*x11*(As*x1 + As - x11) - 1/2*x13 + x14*x15 + x14)/x13)+0j, (x7*(A1*x3 + x10))+0j)

def g3_entries_c(rho, lam, As, A1):
    x0 = math.sqrt(As)
    x1 = rho**2
    x2 = x1 + 1
    x3 = As*rho
    x4 = A1*x2
    x5 = 2*x1
    x6 = x0/x2**(3/2)
    return (((1/2)*(A1*(-2*x3 + x4) - 4*lam*x2)/(x0*math.sqrt(x2)))+0j, (x6*(A1*rho**3 + A1*rho - As*x5 + x5 + 2))+0j, (-math.sqrt(2)*x6*(-x3 + (1/2)*x4))+0j)


# ------------------------------------------------------------------------------------------------ planted values (read from their sources)
# [ABCKY23] = Azad, Blazquez-Salcedo, Chew, Kunz, Yeom, PRD 107, 084024 (2023), arXiv:2212.12601, Table 1, row M/r0 = 0 (massless Ellis wormhole), r0 = 1:
# polar quasinormal modes omega r0 = omega_R + i omega_I, e^{-i omega t}; three printed digits.  Tolerances are half a unit of the last printed digit.
PUB_QNM = {2: (1.74, -0.305), 3: (2.95, -0.410), 4: (4.08, -0.449)}
PUB_TOL = (5e-3, 5e-4)
QNM_THETA = 0.5                                                 # contour rotation for the outgoing conditions, rad (needs |arg omega| < theta)
QNM_THETA_L1 = 0.9                                              # the l = 1 mode has arg omega = -0.73
MA = 1.0                                                        # M-AdS/Ellis radial plant: stability.planted_value() ([TS13] table, W^2 a^2 at D = 4)

# ------------------------------------------------------------------------------------------------ kernels
m4n, m3n, g3n = njit(m4_entries), njit(m3_entries), njit(g3_entries)
m4c, m3c, g3c = njit(m4_entries_c), njit(m3_entries_c), njit(g3_entries_c)


@njit
def Aval(r, c, K):
    """A(rho) of the family; for rho > 5 the series 1 - K sum (-1)^k 2 u^(2k+1)/((2k+1)(2k+3)), u = 1/rho, which has no cancellation"""
    if r > 5.0:
        u = 1.0 / r
        u2 = u * u
        s = 0.0
        term = u
        sg = 1.0
        for k in range(14):
            s += sg * 2.0 * term / ((2 * k + 1) * (2 * k + 3))
            term *= u2
            sg = -sg
        return 1.0 - K * s
    return 1.0 + c * (r * r + 1.0) + K * (r + (r * r + 1.0) * math.atan(r))


@njit
def Apval(r, c, K):
    if r > 5.0:
        u = 1.0 / r
        u2 = u * u
        s = 0.0
        term = u2
        sg = 1.0
        for k in range(14):
            s += sg * 2.0 * term / (2 * k + 3)
            term *= u2
            sg = -sg
        return K * s
    return 2.0 * c * r + 2.0 * K * (1.0 + r * math.atan(r))


@njit
def rhs4(rho, y, lam, L, c, K):
    As = Aval(rho, c, K)
    A1 = Apval(rho, c, K)
    m = m4n(rho, lam, L, As, A1)
    out = np.empty(8)
    for col in range(2):
        for i in range(4):
            s = 0.0
            for j in range(4):
                s += m[i * 4 + j] * y[j * 2 + col]
            out[i * 2 + col] = s
    return out


@njit
def rhs4_hz(z, y, lam, L, c, K, rh):
    """variables (H0, Pi, K, Ft), F = Ft/z: Fuchsian at the horizon"""
    rho = rh + z
    As = Aval(rho, c, K)
    A1 = Apval(rho, c, K)
    m = m4n(rho, lam, L, As, A1)
    out = np.empty(8)
    t = np.array([1.0, 1.0, 1.0, 1.0 / z])
    for col in range(2):
        for i in range(4):
            s = 0.0
            for j in range(4):
                s += m[i * 4 + j] * t[j] / t[i] * y[j * 2 + col]
            if i == 3:
                s += y[i * 2 + col] / z
            out[i * 2 + col] = s
    return out


@njit
def rhs3(rho, y, lam, c, K):
    As = Aval(rho, c, K)
    A1 = Apval(rho, c, K)
    m = m3n(rho, lam, As, A1)
    out = np.empty(3)
    for i in range(3):
        s = 0.0
        for j in range(3):
            s += m[i * 3 + j] * y[j]
        out[i] = s
    return out


@njit
def rhs3_hz(z, y, lam, c, K, rh):
    """variables (H0, Ht, F), H2 = z Ht: Fuchsian at the horizon"""
    rho = rh + z
    As = Aval(rho, c, K)
    A1 = Apval(rho, c, K)
    m = m3n(rho, lam, As, A1)
    out = np.empty(3)
    t = np.array([1.0, z, 1.0])
    for i in range(3):
        s = 0.0
        for j in range(3):
            s += m[i * 3 + j] * t[j] / t[i] * y[j]
        if i == 1:
            s -= y[i] / z
        out[i] = s
    return out


@njit
def rhs4c(s, y, e, lam, L):
    """Ellis wormhole (A = 1) along rho = s e, e = exp(i theta): complex W"""
    r = s * e
    m = m4c(r, lam, L, 1.0, 0.0)
    out = np.empty(8, dtype=np.complex128)
    for col in range(2):
        for i in range(4):
            acc = 0j
            for j in range(4):
                acc += m[i * 4 + j] * y[j * 2 + col]
            out[i * 2 + col] = e * acc
    return out


def integrate_frame(f, args, s0, s1, Y0, chunk, rtol, atol=None):
    """dim x k frame from s0 to s1 in chunks (a float length, or a list of nodes); Gram-Schmidt (QR with positive diagonal) between chunks keeps the span
    and the sign continuous in W"""
    dim, k = Y0.shape
    atol = 0.1 * rtol if atol is None else atol            # the frame is orthonormal at the start of each chunk
    Q, Rm = np.linalg.qr(Y0)
    Q = Q * np.sign(np.diag(Rm))[None, :]
    if isinstance(chunk, (list, tuple)):
        nodes = list(chunk)
    else:
        d = 1.0 if s1 > s0 else -1.0
        n = int(math.ceil(abs(s1 - s0) / chunk - 1e-12))
        nodes = [s0 + d * min(i * chunk, abs(s1 - s0)) for i in range(n + 1)]
    solver = ode(f).set_integrator("dop853", rtol=rtol, atol=atol, nsteps=500000)
    for s, s2 in zip(nodes[:-1], nodes[1:]):
        solver.set_initial_value(Q.reshape(-1), s).set_f_params(*args)
        y = solver.integrate(s2)
        if not solver.successful():
            raise RuntimeError("dop853 failed at s = %g" % s)
        Q, Rm = np.linalg.qr(y.reshape(dim, k))
        Q = Q * np.sign(np.diag(Rm))[None, :]
    return Q


def hz_nodes(z0, zend, expo):
    """chunk boundaries z0 g^k near the horizon, with g^expo <= e^20 so that the frame never spans more than 20 e-folds of growth between normalizations"""
    g = math.exp(min(2.0, 20.0 / max(expo, 1.0)))
    nodes = [z0]
    while nodes[-1] * g < zend:
        nodes.append(nodes[-1] * g)
    nodes.append(zend)
    return nodes


# ------------------------------------------------------------------------------------------------ horizon: ingoing Frobenius solutions
def A_c(r, c, K):
    return 1 + c * (r * r + 1) + K * (r + (r * r + 1) * np.arctan(r))


def Ap_c(r, c, K):
    return 2 * c * r + 2 * K * (1 + r * np.arctan(r))


def S4(M, z):
    t = np.array([1, 1, 1, 1 / z])
    M2 = M * t[None, :] / t[:, None]
    M2[3, 3] += 1 / z
    return z * M2


def S3(M, z):
    t = np.array([1, z, 1])
    M2 = M * t[None, :] / t[:, None]
    M2[1, 1] -= 1 / z
    return z * M2


def taylor_S(entries, nd, tfun, c, rh, lam, extra, nterm=16, rad=0.35, npts=128):
    """Taylor coefficients C_n of S(z) = z M'(z) at the horizon, M' = T^-1 M T - T^-1 T' (T = diag(1,1,1,1/z) for l = 2, diag(1,z,1) for l = 1),
    by FFT on a circle of radius rad around rho_h; z M' is analytic at z = 0 (Fuchsian)"""
    K = -2 * c / math.pi
    zs = rad * np.exp(2j * np.pi * np.arange(npts) / npts)
    Sv = np.empty((npts, nd, nd), dtype=complex)
    for i, z in enumerate(zs):
        r = rh + z
        Sv[i] = tfun(np.array(entries(r, lam, *extra, A_c(r, c, K), Ap_c(r, c, K)), dtype=complex).reshape(nd, nd), z)
    C = np.fft.fft(Sv, axis=0) / npts
    return [C[n] / rad**n for n in range(nterm)]


def frob_vec(C, s, z0, nmax, nrm_idx):
    """solution z^s sum a_n z^n at z0: (s + n - C0) a_n = sum_{k=1..n} C_k a_{n-k}; a_0 the eigenvector of C0 for the exponent s with a_0[nrm_idx] = 1;
    at a resonance (s + n an exponent) the singular system is solved in the least-squares sense; returns the vector and the worst relative
    inconsistency of the recursion (a log term would show there)"""
    nd = C[0].shape[0]
    w, V = np.linalg.eig(C[0])
    j = np.argmin(abs(w - s))
    if abs(w[j] - s) > 1e-6:
        raise RuntimeError("Frobenius exponent %g not an eigenvalue of the residue matrix (%s)" % (s, np.round(w, 6)))
    a = [V[:, j] / V[nrm_idx, j]]
    incons = 0.0
    for n in range(1, nmax + 1):
        rhs = sum(C[k] @ a[n - k] for k in range(1, n + 1))
        Mn = (s + n) * np.eye(nd) - C[0]
        sol = np.linalg.lstsq(Mn, rhs, rcond=1e-12)[0]
        incons = max(incons, np.linalg.norm(Mn @ sol - rhs) / (np.linalg.norm(rhs) + 1e-300))
        a.append(sol)
    v = sum(a[n] * z0**n for n in range(nmax + 1)) * z0**s
    return v.real, incons


# ------------------------------------------------------------------------------------------------ polar determinants (real W)
def match_point(fam):
    if fam.rh is None:
        return 0.0
    return 1.0 if fam.rh < 0.5 else fam.rh + 1.0


def far_frame_l2(W, L, fam, rm, p, side):
    lam = W * W
    rfar = max(p["rmax0"], p["rmaxW"] / W)
    Y0 = np.random.RandomState(3 if side > 0 else 4).randn(4, 2)
    return integrate_frame(rhs4, (lam, L, fam.c, fam.K), side * rfar, rm, Y0, 8.0 / W, p["rtol"])


def D_l2(W, fam, p, L=6.0):
    """det[QL | QR] at the matching point: QL the two ingoing (c < 0: horizon; c = 0: left infinity) solutions, QR the two decaying at +infinity"""
    if mutating("L_sign"):
        L = -L                                             # control: attractive angular term
    if mutating("rmax"):
        p = dict(p, rmax0=4.0, rmaxW=0.0)                  # control: far boundary at rho = 4
    rm = match_point(fam)
    QR = far_frame_l2(W, L, fam, rm, p, +1)
    if fam.rh is None:
        QL = far_frame_l2(W, L, fam, rm, p, -1)
    else:
        lam = W * W
        mu = W / fam.A1
        C = taylor_S(m4_entries, 4, S4, fam.c, fam.rh, lam, (L,))
        sg = -1.0 if mutating("hz_out") else 1.0           # control: the outgoing exponents -mu instead of +mu
        v1, i1 = frob_vec(C, 1 + sg * mu, p["z0"], p["nfrob"], 2)
        v2, i2 = frob_vec(C, -1 + sg * mu, p["z0"], p["nfrob"], 2)
        if max(i1, i2) > 1e-6:
            raise RuntimeError("Frobenius recursion inconsistent: %g" % max(i1, i2))
        QL = integrate_frame(rhs4_hz, (lam, L, fam.c, fam.K, fam.rh), p["z0"], rm - fam.rh, np.stack([v1, v2], axis=1), hz_nodes(p["z0"], rm - fam.rh, 1 + abs(mu)), p["rtol"])
        QL = QL.copy()
        QL[3, :] = QL[3, :] / (rm - fam.rh)                # F = Ft / z
    QL, _ = np.linalg.qr(QL)
    QR, _ = np.linalg.qr(QR)
    return float(np.linalg.det(np.hstack([QL, QR])))


def gauge_vec(rho, lam, fam):
    g = np.array(g3n(rho, lam, float(Aval(rho, fam.c, fam.K)), float(Apval(rho, fam.c, fam.K))))
    if mutating("g1_wrong"):
        g = np.array([0.0, 1.0, 0.0])                      # control: a wrong third column (H2 only) in place of the pure-gauge solution
    return g / np.linalg.norm(g)


def D_l1(W, fam, p):
    """det[QL | QR | g] for l = 1: g the pure-gauge solution (xi^theta = sqrt(A)/R), QL, QR the ingoing and decaying solutions"""
    if mutating("rmax"):
        p = dict(p, rmax0=4.0, rmaxW=0.0)
    lam = W * W
    rm = match_point(fam)
    rfar = max(p["rmax0"], p["rmaxW"] / W)
    QR = integrate_frame(rhs3, (lam, fam.c, fam.K), rfar, rm, np.random.RandomState(3).randn(3, 1), 8.0 / W, p["rtol"])
    if fam.rh is None:
        QL = integrate_frame(rhs3, (lam, fam.c, fam.K), -rfar, rm, np.random.RandomState(4).randn(3, 1), 8.0 / W, p["rtol"])
    else:
        mu = W / fam.A1
        C = taylor_S(m3_entries, 3, S3, fam.c, fam.rh, lam, ())
        sg = -1.0 if mutating("hz_out") else 1.0
        v, inc = frob_vec(C, sg * mu, p["z0"], p["nfrob"], 2)
        if inc > 1e-6:
            raise RuntimeError("Frobenius recursion inconsistent: %g" % inc)
        QL = integrate_frame(rhs3_hz, (lam, fam.c, fam.K, fam.rh), p["z0"], rm - fam.rh, v.reshape(3, 1), hz_nodes(p["z0"], rm - fam.rh, 1 + abs(mu)), p["rtol"])
        QL = QL.copy()
        QL[1, :] *= rm - fam.rh                            # H2 = z Ht
        QL = QL / np.linalg.norm(QL)
    return float(np.linalg.det(np.column_stack([QL[:, 0], QR[:, 0], gauge_vec(rm, lam, fam)])))


def find_roots(Dfun, Ws):
    """sign changes of D on the grid Ws (odd-multiplicity roots, refined by Brent) and interior minima of |D| refined by a bounded search: a minimum
    below DIP_RATIO of its neighbors is a root of even multiplicity (a double root, as at the Ellis wormhole, where the two polar fields are degenerate)"""
    D = np.array([Dfun(W) for W in Ws])
    roots = [brentq(Dfun, Ws[i], Ws[i + 1], xtol=1e-14, rtol=1e-13) for i in range(len(Ws) - 1) if np.sign(D[i]) != np.sign(D[i + 1])]
    dips = []
    for i in range(1, len(Ws) - 1):
        if abs(D[i]) < abs(D[i - 1]) and abs(D[i]) < abs(D[i + 1]):
            r = minimize_scalar(lambda W: abs(Dfun(W)), bounds=(Ws[i - 1], Ws[i + 1]), method="bounded", options=dict(xatol=1e-12))
            ratio = r.fun / min(abs(D[i - 1]), abs(D[i + 1]))
            dips.append(dict(W=r.x, Dmin=r.fun, ratio=ratio, flagged=bool(ratio < DIP_RATIO)))
    return roots, dips, D


def count_polar(ell, c, res, L=None):
    fam = st.Fam(c, "R2")
    p = RES[res]
    Dfun = (lambda W: D_l2(W, fam, p, 6.0 if L is None else L)) if ell == 2 else (lambda W: D_l1(W, fam, p))
    lo = 3.0 if mutating("scan_range") else W_LO
    Ws = np.exp(np.linspace(math.log(lo), math.log(W_HI), p["nW"]))
    t0 = time.time()
    roots, dips, D = find_roots(Dfun, Ws)
    probes = {("%g" % W): Dfun(W) for W in PROBES}
    return dict(ell=ell, c=c, res=res, rh=fam.rh, A1=fam.A1, n_roots=len(roots) + 2 * sum(d["flagged"] for d in dips), roots=roots, dips=dips,
                sign_changes=len(roots), D_first=float(D[0]), D_last=float(D[-1]), D_absmin=float(np.min(np.abs(D))), probes=probes,
                W_lo=lo, W_hi=W_HI, nW=p["nW"], seconds=time.time() - t0, mutate=sorted(MUT))


# ------------------------------------------------------------------------------------------------ polar quasinormal modes (complex W), Ellis wormhole
def frame_c(rhs, args, s0, s1, dim, k, chunk, seed):
    rng = np.random.RandomState(seed)
    Q, _ = np.linalg.qr(rng.randn(dim, k) + 1j * rng.randn(dim, k))
    s = s0
    d = np.sign(s1 - s0)
    while abs(s - s1) > 1e-12:
        s2 = s + d * min(abs(s1 - s), chunk)
        sol = solve_ivp(lambda t, y: rhs(t, y, *args), [s, s2], Q.reshape(-1), method="DOP853", rtol=1e-11, atol=1e-30)
        Q, Rm = np.linalg.qr(sol.y[:, -1].reshape(dim, k))
        Q = Q * (np.diag(Rm) / np.abs(np.diag(Rm)))[None, :]
        s = s2
    return Q


def D_qnm_l2(Wc, L, theta=QNM_THETA):
    """Ellis polar determinant at complex W = -i omega; the radial coordinate is rotated, rho = s e^{i theta}, so that both outgoing conditions
    (decay of exp(-+W rho)) are decaying solutions along the contour; matching at rho = 0"""
    if mutating("qnm_L"):
        L = L + 1.0                                        # control: wrong angular eigenvalue
    e = np.exp(1j * theta)
    mu = (Wc * e).real
    S = 22.0 / abs(mu)
    chunk = 6.0 / abs(mu)
    args = (e, Wc * Wc, L)
    QR = frame_c(rhs4c, args, S, 0.0, 4, 2, chunk, 3)
    QL = frame_c(rhs4c, args, -S, 0.0, 4, 2, chunk, 4)
    return complex(np.linalg.det(np.hstack([QL, QR])))


def qnm_l2(ell, start=None):
    L = ell * (ell + 1.0)
    om0 = complex(*PUB_QNM[ell]) if start is None else start
    f = lambda v: (lambda d: [d.real, d.imag])(D_qnm_l2(-1j * (v[0] + 1j * v[1]), L))
    r = root(f, [om0.real, om0.imag], tol=1e-11)
    return complex(r.x[0], r.x[1]), bool(r.success), float(np.hypot(*f(r.x)))


@njit
def rhs3c(s, y, e, lam):
    """Ellis wormhole, polar l = 1 system (H0, H2, F) along rho = s e"""
    r = s * e
    m = m3c(r, lam, 1.0, 0.0)
    out = np.empty(3, dtype=np.complex128)
    for i in range(3):
        acc = 0j
        for j in range(3):
            acc += m[i * 3 + j] * y[j]
        out[i] = e * acc
    return out


@njit
def rhs_schr(s, y, e, lam, kind):
    """Z'' = (V + lam) Z along rho = s e; kind 0: V = 2/R^2 - 3/R^4 ([ABCKY23] (53) at L = 2), kind 1: V = 2/R^2 + 1/R^4 (test scalar, l = 1)"""
    r = s * e
    R2 = r * r + 1.0
    V = 2.0 / R2 - 3.0 / (R2 * R2) if kind == 0 else 2.0 / R2 + 1.0 / (R2 * R2)
    return np.array([e * y[1], e * (V + lam) * y[0]], dtype=np.complex128)


def D_qnm_l1(Wc, theta=QNM_THETA_L1):
    """Ellis polar l = 1 determinant det[QL | QR | g] at complex W, rotated contour as in D_qnm_l2"""
    e = np.exp(1j * theta)
    mu = (Wc * e).real
    S, chunk = 22.0 / abs(mu), 6.0 / abs(mu)
    args = (e, Wc * Wc)
    QR = frame_c(rhs3c, args, S, 0.0, 3, 1, chunk, 3)
    QL = frame_c(rhs3c, args, -S, 0.0, 3, 1, chunk, 4)
    g = np.array(g3c(0.0, Wc * Wc, 1.0, 0.0), dtype=complex)
    return complex(np.linalg.det(np.column_stack([QL[:, 0], QR[:, 0], g / np.linalg.norm(g)])))


def D_qnm_schr(Wc, kind, theta=QNM_THETA_L1):
    """the independent 1D code: Z'' = (V + W^2) Z, outgoing at both ends; shares nothing with the polar system"""
    if mutating("l1_pot"):
        kind = 1                                           # control: the test-scalar l = 1 potential in place of the master equation at L = 2
    e = np.exp(1j * theta)
    mu = (Wc * e).real
    S, chunk = 22.0 / abs(mu), 6.0 / abs(mu)
    args = (e, Wc * Wc, kind)
    QR = frame_c(rhs_schr, args, S, 0.0, 2, 1, chunk, 3)
    QL = frame_c(rhs_schr, args, -S, 0.0, 2, 1, chunk, 4)
    return complex(QL[0, 0] * QR[1, 0] - QL[1, 0] * QR[0, 0])


def qnm_root(Dfun, om0):
    f = lambda v: (lambda d: [d.real, d.imag])(Dfun(-1j * (v[0] + 1j * v[1])))
    r = root(f, [om0.real, om0.imag], tol=1e-11)
    return complex(r.x[0], r.x[1]), bool(r.success)


def fd_bound_states(Lv, X=400.0, N=160001, k=6):
    """bound states W of -Z'' + V Z = -W^2 Z, V = L/R^2 - 3/R^4, finite differences on two grids with Richardson extrapolation in h^2 (the other 1D code)"""
    from scipy.linalg import eigh_tridiagonal
    out = []
    for n in (N, 2 * N - 1):
        x = np.linspace(-X, X, n)
        h = x[1] - x[0]
        V = Lv / (x**2 + 1) - 3 / (x**2 + 1)**2
        out.append(eigh_tridiagonal(2 / h**2 + V, -np.ones(n - 1) / h**2, select="i", select_range=(0, k - 1), eigvals_only=True))
    lam = (4 * out[1] - out[0]) / 3
    return [math.sqrt(-e) for e in lam if e < -1e-3]


def cscale_qnm(Vf, near, theta=1.1, X=60.0, N=20000):
    """the quasinormal mode of Z'' + (omega^2 - V) Z = 0 nearest to `near`, as a discrete eigenvalue of the complex-scaled operator
    -e^{-2 i theta} d^2/ds^2 + V(s e^{i theta}) on s in [-X, X] (finite differences, Dirichlet ends; shift-invert Arnoldi), Richardson in h^2 from N and 2N;
    the other independent 1D code (an eigenvalue problem, no shooting)"""
    import scipy.sparse as sps
    from scipy.sparse.linalg import eigs
    lam = []
    for n in (N, 2 * N):
        s = np.linspace(-X, X, n + 2)[1:-1]
        h = s[1] - s[0]
        ph = np.exp(-2j * theta)
        H = sps.diags([-np.ones(n - 1) / h**2 * ph, 2 / h**2 * ph + Vf(s * np.exp(1j * theta)), -np.ones(n - 1) / h**2 * ph], [-1, 0, 1], format="csc")
        w, _ = eigs(H, k=1, sigma=near**2, which="LM")
        lam.append(w[0])
    return np.sqrt((4 * lam[1] - lam[0]) / 3)


# ------------------------------------------------------------------------------------------------ M-AdS wormholes (c > 0): radial modes
# The radial equation, the gauge delta phi = 0, the right-hand condition and the analytic throat crossing are those of stability.py (#212).  New: the
# AdS end rho -> -infinity, where A ~ 2 c rho^2 and Y -> alpha + beta rho; Dirichlet (no source, Psi = -sqrt2 Y/rho -> 0) is beta = 0, the series
# Y = 1 + y2/rho^2 + y4/rho^4 + ... (sympy, polar.md), integrated toward the throat.
def ads_start(fam, W, p):
    c = fam.c
    r0 = -p["ads"] * max(100.0, 400.0 * W / c, 100.0 / math.sqrt(c))
    if mutating("ads_bc"):
        return r0, r0, 1.0                                 # control: the non-normalizable solution Y = rho
    y2 = (W**2 - 4 * c**2 - 4 * c) / (24 * c**2)
    y4 = (W**4 - 40 * W**2 * c + 272 * c**4 + 416 * c**3 + 24 * c**2 * (4 - 3 * W**2)) / (1920 * c**4)
    return r0, 1 + y2 / r0**2 + y4 / r0**4, -2 * y2 / r0**3 - 4 * y4 / r0**5


def rhs_rho(fam, W):
    e = fam.eps

    def f(r, s):
        Y, Yr = s
        A = fam.A(r)
        Ap = fam.Ap(r)
        R2 = r * r + 1
        P1 = -2 * e / (r * R2) + 2 * r / R2 - Ap / A
        P0 = -e / R2**2 + e * Ap / (A * r * R2) - 2 / (A * R2)
        return [Yr, P1 * Yr + P0 * Y + W * W * Y / A**2]
    return f


def left_ads(fam, W, rm, p):
    r0, Y, Yp = ads_start(fam, W, p)
    d, rtol = p["delta"], p["rtol"]
    s = np.array([Y, Yp])
    s /= np.hypot(*s)
    f = rhs_rho(fam, W)
    r = r0
    while r < -d - 1e-12:
        r2 = min(-d, r + max(5.0, abs(r) * 0.5))
        sol = solve_ivp(f, [r, r2], s, method="DOP853", rtol=rtol, atol=[1e-30, 1e-30])
        s = sol.y[:, -1]
        s /= np.hypot(*s)
        r = r2
    sol = solve_ivp(lambda th, y: st.rhs_arc(fam, W, d)(th, y), [np.pi, 0.0], [s[0] + 0j, s[1] + 0j], method="DOP853", rtol=rtol, atol=1e-30)
    Yc, Yrc = sol.y[:, -1]
    flag = abs(Yc.imag) / abs(Yc)                         # monodromy test: no logarithm at the throat
    s0 = np.array([d, Yc.real, (Yrc * fam.A(d)).real])
    s0[1:] /= np.hypot(s0[1], s0[2])
    stop = st.make_ev(lambda x, ss: ss[0] - rm, 1, True)
    s, _, _ = st.run_chunks(st.rhs_x(fam, W), s0, +1, W, rtol, [1e-12, 1e-30, 1e-30], stop, None)
    return s[0], s[1], s[2], flag


def mads_mismatch(fam, W, p, flags=None):
    rm = 1.0
    _, Yl, Ylx, flag = left_ads(fam, W, rm, p)
    if flags is not None:
        flags.append(flag)
    fam.res = dict(rtol=p["rtol"], rmax0=p["rmax0"], rmaxW=p["rmaxW"], delta=p["delta"])
    _, Yr, Yrx, _ = st.right_solution(fam, W, rm)
    return Yl * Yrx - Ylx * Yr


def count_mads(c, res):
    fam = st.Fam(c, "R2")
    p = RES[res]
    lo = 3.0 if mutating("scan_range") else W_LO
    Ws = np.exp(np.linspace(math.log(lo), math.log(W_HI_M), p["nWm"]))
    t0 = time.time()
    flags = []
    D = np.array([mads_mismatch(fam, W, p, flags) for W in Ws])
    roots = [brentq(lambda W: mads_mismatch(fam, W, p, flags), Ws[i], Ws[i + 1], xtol=1e-13, rtol=1e-12) for i in range(len(Ws) - 1) if np.sign(D[i]) != np.sign(D[i + 1])]
    return dict(c=c, res=res, n_roots=len(roots), roots=[dict(W=w, W2=w * w) for w in roots], sign_changes=len(roots), D_first=float(D[0]), D_last=float(D[-1]),
                monodromy_imag_max=float(max(flags)), W_lo=lo, W_hi=W_HI_M, nW=p["nWm"], seconds=time.time() - t0, mutate=sorted(MUT))


def spectral_mads(c, N, X):
    """Chebyshev collocation of the Y-form Y_xx = Cx Y_x + Dx Y + lam Y on x in [x_L, X] (x_L the finite AdS end): Neumann at x_L (the Dirichlet
    solution Y = 1 + O(s^2)), Y = 0 at X; rho(x) by inversion of x(rho).  Shares with the shooting code the background A and the coefficients of the Y-form,
    nothing of the AdS series, the throat crossing, the matching or the root search."""
    from scipy.integrate import quad
    from scipy.linalg import eig
    fam = st.Fam(c, "R2")
    xL = -quad(lambda r: 1 / fam.A(r), -np.inf, 0, epsabs=1e-14, epsrel=1e-13, limit=400)[0]
    D, xc = st.cheb(N)
    xs = X + (xc - 1.0) * (X - xL) / 2                    # xc = 1 -> X, xc = -1 -> xL
    rho = np.empty(N + 1)
    for j in range(N + 1):
        xj = xs[j]
        if j == N:
            rho[j] = -np.inf
        elif xj >= 0:
            rho[j] = brentq(lambda r: quad(lambda q: 1 / fam.A(q), 0, r, epsabs=1e-14, epsrel=1e-13)[0] - xj, 0, xj + 10, xtol=1e-14) if xj > 0 else 0.0
        else:
            s = xj - xL
            rho[j] = brentq(lambda r: quad(lambda q: 1 / fam.A(q), -np.inf, r, epsabs=1e-16, epsrel=1e-12, limit=400)[0] - s, -1e12, -1e-9, xtol=1e-14)
    rho[np.abs(rho) < 1e-6] = 1e-6                        # a node at the throat: its Y-form coefficients have the removable pole 1/rho
    A = np.empty(N + 1)
    Ap = np.empty(N + 1)
    for j in range(N):
        A[j], Ap[j] = fam.A(rho[j]), fam.Ap(rho[j])
    A[N], Ap[N] = np.inf, np.inf
    Cx = np.zeros(N + 1)
    Dx = np.zeros(N + 1)
    for j in range(N):
        Cx[j], Dx[j] = st.coef_x(fam, rho[j], A[j], A[j] / rho[j], Ap[j])
    D1 = D * 2 / (X - xL) * -1.0                           # d/dx = (2/(X - xL)) d/dxc, xc decreasing with x
    D1 = D / ((X - xL) / 2)
    D2 = D1 @ D1
    sgn = -1.0 if mutating("spec_sign") else 1.0
    Aop = np.zeros((N + 1, N + 1))
    Bop = np.zeros((N + 1, N + 1))
    for j in range(1, N):
        Aop[j, :] = D2[j, :] - Cx[j] * D1[j, :] - sgn * Dx[j] * np.eye(N + 1)[j, :]
        Bop[j, j] = 1.0
    Aop[0, 0] = 1.0                                        # Y(X) = 0
    Aop[N, :] = D1[N, :]                                   # Y_x(x_L) = 0
    ev = eig(Aop, Bop, right=False)
    ev = ev[np.isfinite(ev)]
    return ev[(ev.real > 1e-4) & (ev.real < 100.0)]


SPEC = {"S-A": dict(N=240, X=40.0, tol=2e-3), "S-B": dict(N=360, X=60.0, tol=2e-4)}
SPEC_CMIN = 0.01                                           # below it the AdS end is at |rho| > 70 and the finite domain is not resolved (development: 6e-2 at c = 1e-4):
                                                           # there the count only is compared, the rates are reported


def spectral_command(d):
    ok_all = True
    rows = {}
    for fn in sorted(glob.glob(d + "/mads_*_R2.json")):
        o = json.load(open(fn))
        rows[o["c"]] = o
    for c in sorted(rows):
        ref = sorted(r["W"] for r in rows[c]["roots"])
        for which in ("S-A", "S-B"):
            ev = spectral_mads(c, SPEC[which]["N"], SPEC[which]["X"])
            realpos = sorted(e.real for e in ev if abs(e.imag) < 1e-6 * abs(e) and e.real < W_HI_M**2)
            cplx = [e for e in ev if abs(e.imag) >= 1e-6 * abs(e)]
            rel = [abs(math.sqrt(a) - b) / b for a, b in zip(realpos, ref)]
            tol = SPEC[which]["tol"] if c >= SPEC_CMIN else float("inf")
            ok = len(realpos) == len(ref) and not cplx and all(r <= tol for r in rel)
            ok_all &= ok
            print(("PASS  " if ok else "FAIL  ") + f"S6 spectral {which} c = {c:g}: eigenvalues with Re lam > 1e-4: real {len(realpos)} (shooting {len(ref)}), complex {len(cplx)}, "
                  f"|dW|/W = {['%.2e' % r for r in rel]} (declared <= {SPEC[which]['tol'] if c >= SPEC_CMIN else 'n/a, c < %g' % SPEC_CMIN})", flush=True)
    return ok_all


# ------------------------------------------------------------------------------------------------ planted cases
# P3 (ruling of the coordinator, 2026-11-27): no published number for polar l = 1.  The derived l = 1 system at c = 0 has the quasinormal mode of
# [ABCKY23]'s master equation (52), (53) at L = 2 (the two are isospectral, Darboux partners, polar_derive.py D6e, D6f), computed by an independent 1D code.
L1_START = complex(0.6814, -0.6177)                             # starting point of the root searches (complex-scaling development run: 0.681429 - 0.617763i)
L1_TOL = 1e-7                                                   # |omega_derived - omega_master|, both by shooting
L1_TOL_FD = 1e-8                                                # against the finite-difference eigenvalue code (theta = 1.1, X = 60, N = 20000, Richardson)
# P4: the Ellis wormhole at the unphysical L = -1: the decoupled master equation has bound states of known W (finite differences); each is a double
# root of the polar determinant (two identical fields), found as a flagged dip of |D|.
P4_L, P4_WLO, P4_TOL = -1.0, 0.2, 5e-5                          # relative tolerance on W
# P5: the M-AdS radial mode at c -> 0+: W^2(c) - W^2_Ellis = a1 c + a2 c^2 + ..., quadratic fit through c = 1e-4, 1e-3, 1e-2
P5_C, P5_TOL = (1e-4, 1e-3, 1e-2), 1e-8


def planted_command():
    ok_all = True

    def line(tag, ok, msg):
        nonlocal ok_all
        ok_all &= ok
        print(("PASS  " if ok else "FAIL  ") + tag + ": " + msg, flush=True)
    t0 = time.time()
    for ell in (2, 3, 4):
        om, conv, res_ = qnm_l2(ell)
        pr, pi = PUB_QNM[ell]
        ok = conv and abs(om.real - pr) <= PUB_TOL[0] and abs(om.imag - pi) <= PUB_TOL[1]
        line("P2 Ellis polar l = %d quasinormal mode against [ABCKY23] Table 1 (M = 0)" % ell, ok,
             "omega = %.6f %+.6fi (printed %.2f %+.3fi), |dRe| = %.1e, |dIm| = %.1e (declared <= %g, %g), |D| at the root %.1e (%.0f s)" %
             (om.real, om.imag, pr, pi, abs(om.real - pr), abs(om.imag - pi), PUB_TOL[0], PUB_TOL[1], res_, time.time() - t0))
    om1, c1 = qnm_root(D_qnm_l1, L1_START)
    om2, c2 = qnm_root(lambda W: D_qnm_schr(W, 0), L1_START)
    line("P3a Ellis polar l = 1 quasinormal mode of the derived system against the master equation (52), (53) of [ABCKY23] at L = 2 (1D shooting)",
         c1 and c2 and abs(om1 - om2) <= L1_TOL, "derived %.9f %+.9fi, master %.9f %+.9fi, |difference| = %.1e (declared <= %g) (%.0f s)" %
         (om1.real, om1.imag, om2.real, om2.imag, abs(om1 - om2), L1_TOL, time.time() - t0))
    Vm = (lambda r: 2 / (r * r + 1) + 1 / (r * r + 1)**2) if mutating("l1_pot") else (lambda r: 2 / (r * r + 1) - 3 / (r * r + 1)**2)
    om3 = cscale_qnm(Vm, L1_START)
    line("P3b the same against the complex-scaled eigenvalue problem for the master equation at L = 2 (finite differences, Richardson in h^2)", abs(om1 - om3) <= L1_TOL_FD,
         "derived %.9f %+.9fi, eigenvalue code %.9f %+.9fi, |difference| = %.1e (declared <= %g)" % (om1.real, om1.imag, om3.real, om3.imag, abs(om1 - om3), L1_TOL_FD))
    ref = fd_bound_states(P4_L)
    ref = [w for w in ref if w >= P4_WLO]
    fam = st.Fam(0.0, "R2")
    p = RES["R1"]
    Ws = np.exp(np.linspace(math.log(P4_WLO), math.log(W_HI), p["nW"]))
    roots, dips, _ = find_roots(lambda W: D_l2(W, fam, p, P4_L), Ws)
    found = sorted(d["W"] for d in dips if d["flagged"])
    ok = len(found) == len(ref) and not roots and all(abs(a - b) / b <= P4_TOL for a, b in zip(found, sorted(ref)))
    line("P4 Ellis, L = %g (unphysical): double roots of the polar determinant against the bound states of the master equation (finite differences)" % P4_L, ok,
         "found %s, finite differences %s, sign changes %d (declared relative <= %g)" % (["%.6f" % w for w in found], ["%.6f" % w for w in sorted(ref)], len(roots), P4_TOL))
    ref2 = float(st.planted_value())
    ds = []
    for c in P5_C:
        o = count_mads(c, "R1")
        ds.append(o["roots"][0]["W2"] - ref2 if o["n_roots"] == 1 else float("nan"))
    a0 = np.linalg.solve(np.vander(np.array(P5_C), 3, increasing=True), np.array(ds))[0]
    line("P5 M-AdS radial mode as c -> 0+ against the Ellis value %s (stability.py planted value)" % st.planted_value(), all(np.isfinite(ds)) and abs(a0) <= P5_TOL,
         "W^2(c) - W^2_Ellis = %s at c = %s; quadratic extrapolation to c = 0: %.2e (declared <= %g) (%.0f s)" % (["%.6e" % x for x in ds], P5_C, a0, P5_TOL, time.time() - t0))
    return ok_all


# ------------------------------------------------------------------------------------------------ report
ROOT_AGREE = 1e-6                                               # relative, M-AdS rates R1 against R2
PROBE_TOL = 1e-6                                                # relative, polar determinants at PROBES, R1 against R2
MONO_TOL = 1e-6                                                 # relative imaginary part of Y after the throat arc (no logarithm)


def load(d, pattern):
    rows = {}
    for fn in glob.glob(d + "/" + pattern):
        o = json.load(open(fn))
        rows.setdefault((o.get("ell"), o["c"]), {})[o["res"]] = o
    return rows


def report_command(d):
    ok_all = True
    rows = load(d, "polar_l*_*.json")
    for ell in (2, 1):
        print("Polar l = %d, real W in [%g, %g] (units 1/b), shooting; n = sign changes + 2 (flagged dips)" % (ell, W_LO, W_HI))
        print("%8s %10s %9s | %-3s %-3s | %-5s %-12s | %-24s %s" % ("c", "rho_h", "kappa", "nR1", "nR2", "dips", "max|dD|/|D|", "min|D| R1, R2", "status"))
        for (e, c) in sorted(rows):
            if e != ell:
                continue
            r1, r2 = rows[(e, c)].get("R1"), rows[(e, c)].get("R2")
            n1 = r1["n_roots"] if r1 else None
            n2 = r2["n_roots"] if r2 else None
            rel = max((abs(r1["probes"][k] - r2["probes"][k]) / max(abs(r1["probes"][k]), abs(r2["probes"][k])) for k in r1["probes"]), default=0.0) if r1 and r2 else float("nan")
            nd = sum(len(x["dips"]) for x in (r1, r2) if x)
            ok = r1 is not None and r2 is not None and n1 == n2 and rel <= PROBE_TOL
            ok_all &= ok
            ref = r2 or r1
            print("%8g %10s %9s | %-3s %-3s | %-5d %-12.2e | %-10.2e %-12.2e %s" % (
                c, "-" if ref["rh"] is None else "%.6f" % ref["rh"], "-" if ref["A1"] is None else "%.6f" % (ref["A1"] / 2), n1, n2, nd, rel,
                r1["D_absmin"] if r1 else float("nan"), r2["D_absmin"] if r2 else float("nan"), "ok" if ok else "FAIL"))
        print()
    mrows = {}
    for fn in glob.glob(d + "/mads_*_R*.json"):
        o = json.load(open(fn))
        mrows.setdefault(o["c"], {})[o["res"]] = o
    print("M-AdS wormholes, radial modes, Dirichlet at the AdS end, W in [%g, %g]" % (W_LO, W_HI))
    print("%8s | %-3s %-3s | %-18s %-18s | %-10s %-10s %s" % ("c", "nR1", "nR2", "W R1", "W R2", "R1 v R2", "monodromy", "status"))
    for c in sorted(mrows):
        r1, r2 = mrows[c].get("R1"), mrows[c].get("R2")
        w1 = [x["W"] for x in r1["roots"]] if r1 else []
        w2 = [x["W"] for x in r2["roots"]] if r2 else []
        rel = max([abs(a - b) / b for a, b in zip(w1, w2)] or [0.0])
        mono = max(x["monodromy_imag_max"] for x in (r1, r2) if x)
        ok = r1 is not None and r2 is not None and len(w1) == len(w2) and rel <= ROOT_AGREE and mono <= MONO_TOL
        ok_all &= ok
        print("%8g | %-3s %-3s | %-18s %-18s | %-10.1e %-10.1e %s" % (c, len(w1), len(w2), ",".join("%.10f" % w for w in w1) or "none", ",".join("%.10f" % w for w in w2) or "none", rel, mono, "ok" if ok else "FAIL"))
    return ok_all


# ------------------------------------------------------------------------------------------------ controls
# kind: count = the count of the control member must differ from the unmutated count; planted/spectral = the command must exit 1 with a FAIL line of the target
CONTROLS = [
    ("L_sign", "count", ("l2", -0.3)), ("hz_out", "count", ("l2", -1.5)), ("rmax", "count", ("l1", -0.3)), ("g1_wrong", "count", ("l1", -1.5)),
    ("scan_range", "count", ("mads", 0.5)), ("ads_bc", "specmads", "S6 spectral"), ("qnm_L", "planted", "P2 Ellis polar l = 2"),
    ("l1_pot", "planted", "P3a Ellis polar l = 1"), ("spec_sign", "spectral", "S6 spectral"),
]


def run_json(args):
    r = subprocess.run([sys.executable, __file__] + args, capture_output=True, text=True)
    return json.loads(r.stdout[r.stdout.index("{"):])


def controls_command():
    ok_all = True
    import tempfile
    for m, kind, tgt in CONTROLS:
        if kind == "count":
            sector, c = tgt
            args = {"l1": ["count", "1", str(c)], "l2": ["count", "2", str(c)], "mads": ["mads", str(c)]}[sector] + ["--res=R1"]
            base = run_json(args)
            mut = run_json(args + ["--mutate=" + m])
            fired = mut["n_roots"] != base["n_roots"]
            detail = f"{sector} c = {c}: n = {mut['n_roots']} against {base['n_roots']} unmutated"
            ret = 0
        else:
            if kind == "planted":
                r = subprocess.run([sys.executable, __file__, "planted", "--mutate=" + m], capture_output=True, text=True)
            elif kind == "specmads":                         # the shooting with the control, the spectral check without it
                tmp = tempfile.mkdtemp(prefix="nm_i218_ctl_")
                open(tmp + "/mads_0.5_R2.json", "w").write(json.dumps(run_json(["mads", "0.5", "--res=R1", "--mutate=" + m]) | {"res": "R2"}))
                r = subprocess.run([sys.executable, __file__, "spectral", tmp], capture_output=True, text=True)
            else:
                tmp = tempfile.mkdtemp(prefix="nm_i218_ctl_")
                open(tmp + "/mads_0.5_R2.json", "w").write(json.dumps(count_mads(0.5, "R1") | {"res": "R2"}))
                r = subprocess.run([sys.executable, __file__, "spectral", tmp, "--mutate=" + m], capture_output=True, text=True)
            txt = r.stdout + r.stderr
            fired = r.returncode == 1 and any(tgt in l and l.startswith("FAIL") for l in txt.splitlines())
            detail = "; ".join(l[:150] for l in txt.splitlines() if l.startswith("FAIL"))[:300]
            ret = r.returncode
        print(("fires   " if fired else "SILENT  ") + f"--mutate={m} (target {tgt}): exit {ret}; {detail}", flush=True)
        ok_all &= fired
    return ok_all


# ------------------------------------------------------------------------------------------------ driver
def main():
    global MUT
    argv = [a for a in sys.argv[1:] if not a.startswith("--mutate=")]
    for a in sys.argv[1:]:
        if a.startswith("--mutate="):
            MUT |= set(a.split("=", 1)[1].split(","))
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if a.startswith("--") and "=" in a}
    pos = [a for a in argv if not a.startswith("--")]
    cmd = pos[0] if pos else None
    if cmd == "planted":
        sys.exit(0 if planted_command() else 1)
    if cmd == "count":
        o = count_polar(int(pos[1]), float(pos[2]), opts.get("--res", "R2"))
        s = json.dumps(o, indent=1)
        if "--out" in opts:
            open(opts["--out"], "w").write(s)
        print(s)
        return
    if cmd == "mads":
        o = count_mads(float(pos[1]), opts.get("--res", "R2"))
        s = json.dumps(o, indent=1)
        if "--out" in opts:
            open(opts["--out"], "w").write(s)
        print(s)
        return
    if cmd == "spectral":
        sys.exit(0 if spectral_command(pos[1]) else 1)
    if cmd == "report":
        sys.exit(0 if report_command(pos[1]) else 1)
    if cmd == "controls":
        sys.exit(0 if controls_command() else 1)
    raise SystemExit(__doc__)


if __name__ == "__main__":
    main()

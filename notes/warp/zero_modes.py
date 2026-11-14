#!/usr/bin/env python3
"""Issue #69: a search for real zero modes of the Alcubierre Witten operator D~ = D + (1/2)(tr K) g0, up to 10 v_c.

#63 proved that D~(v) fails to be invertible only on a discrete set E of real speeds, |v| >= sqrt(3) v_c, and that E is
empty when tr K has one sign. Alcubierre's tr K = -v f'(r) cos(theta) has both signs, so E empty is a conjecture. This
script scans s_min(v) = sigma_min(D~ D^-1) for real v up to 10 v_c at two boxes and two grids, and refines any minimum
below 1e-3. s_min(v) = 0 on the grid exactly when -1/v is an eigenvalue of M = B D^-1, so a zero mode is a dip to 0.

    python notes/warp/zero_modes.py [--part=adjoint,dense,vc,regress63,sectors,lanczos,planted,ly]   planted truth first
    python notes/warp/zero_modes.py --part=scan --sigma=2 --L=8 --N=64 [--tmax=10] [--out=raw/scan.jsonl]
    python notes/warp/zero_modes.py --part=table --raw=DIR       merge the jsonl of the scan runs into the tables
    python notes/warp/zero_modes.py --controls [--raw=DIR] [--only=a,b]   every --mutate control, each must make its check
                                                                 fail; gridperturb needs --raw

Operators (grid: periodic box [-L, L)^3, spectral D^-1 with k = 0 set to zero, A(v) = I + v M, M = B D^-1):
  witten4   4-spinors, D = g.grad, B = (1/2) tau1 g0, tau1 = tr K / v = -f'(r) cos(theta). The literal D~ D^-1.
  witten2   2-spinors of the psi = a + i b sector of #63 section 2, D = sigma.grad, B = -(1/2) tau1. The chi sector is
            v -> -v. D~ D^-1 is unitarily the direct sum of the two, so its s_min is the smaller of the two; `sectors`
            checks it against witten4 on the grid. It halves the work.
  planted   witten2 with tau1 = 6/(1 + r^2): zero modes at v = +-i (FL appendix of #63).
  ly        2-spinors, B = -i sigma.A with the Loss-Yau potential A: a zero mode at the real v = 1, to test the real scan.
s_min is taken on the full space (a zero mode of the grid operator is an eigenvector of M). #63's smin_Dt restricts to
mean-free fields, which a periodic box needs for the k = 0 mode of B psi and which lifts s_min near a zero; `regress63`
reproduces #63's numbers with that restriction (free=True).

Certificate. s_min is Lipschitz in v with constant ||M||, so on [v_a, v_b] it is at least
(s_a + s_b - ||M|| (v_b - v_a))/2. An interval where that is not positive is bisected. Intervals that stay uncertified
down to 1/16 of the step are minimised, and a minimum below 1e-3 is a zero-mode candidate.

Mutations (each must make its check fail): half, conj, adj, bscale, rot, lip, lyscale, lowiter, lipbig, mem, gridperturb.
Exits 1 when any check fails. Peak memory is a few GB at N = 128; run under nice with OMP_NUM_THREADS=1.
"""
import argparse
import glob
import json
import math
import os
import resource
import subprocess
import sys
import time

import numpy as np
import scipy.fft as sfft
from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, eigs, eigsh

ap = argparse.ArgumentParser()
ap.add_argument("--part", default="adjoint,dense,vc,regress63,sectors,lanczos,planted,ly")
ap.add_argument("--mutate", default="")
ap.add_argument("--sigma", type=float, default=2.0)
ap.add_argument("--spinor", type=int, default=2, choices=(2, 4))
ap.add_argument("--L", type=float, default=0.0, help="box half-width; 0 uses the part's default")
ap.add_argument("--N", type=int, default=0, help="points per side; 0 uses the part's default")
ap.add_argument("--tmax", type=float, default=10.0, help="scan up to tmax units of v_c")
ap.add_argument("--step", type=float, default=0.25, help="scan step in units of v_c")
ap.add_argument("--out", default="", help="jsonl file the scan appends its records to")
ap.add_argument("--raw", default="", help="directory of scan jsonl files, for --part=table")
ap.add_argument("--confirm", action="store_true", help="planted and ly: also the box and grid convergence of the zero")
ap.add_argument("--controls", action="store_true")
ap.add_argument("--only", default="", help="with --controls: run only these mutations, comma separated")
args = ap.parse_args()
MUT = args.mutate
PARTS = args.part.split(",")
FAILED = []
MEM_LIMIT = 0.05 if MUT == "mem" else 8.0  # GB, the acceptance limit of #69
LANCZOS = {"runs": 0, "unconverged": 0}
PLANTED_CFG = "--part=planted --L=10 --N=40"
LY_CFG = "--part=ly --L=10 --N=32"
MUTATIONS = {  # name: (what it breaks, the run's options, the checks that must fail)
    "half": ("drops the 1/2 in B (witten2 and witten4) in the FFT operator",
             "--part=dense,regress63,planted --L=10 --N=40", ["dense", "regress63", "planted_pm_i"]),
    "conj": ("Gram operator uses v instead of conj(v); only dense sees it, a zero of A is a zero of any (I + w T) A",
             "--part=dense", ["dense"]),
    "adj": ("treats the anti-Hermitian B of ly as Hermitian in M^H", "--part=adjoint", ["adjoint"]),
    "bscale": ("scales B of witten2 by 11/10", "--part=sectors", ["sectors"]),
    "rot": ("scans the planted tau along the imaginary axis", "--part=planted --L=10 --N=40 --tmax=3", ["planted_real"]),
    "lip": ("uses a tenth of ||M|| in the certificate", "--part=ly --L=10 --N=32", ["ly_real"]),
    "lyscale": ("scales the Loss-Yau potential A by 11/10 in the operator", "--part=ly --L=10 --N=32", ["ly_real"]),
    "lowiter": ("stops Lanczos after two restarts", "--part=lanczos", ["lanczos"]),
    "lipbig": ("uses 50 ||M|| in the certificate", "--part=scan --sigma=2 --L=8 --N=32 --tmax=1 --step=1", ["scan_resolved"]),
    "mem": ("sets the memory limit to 0.05 GB", "--part=vc", ["memory"]),
    "gridperturb": ("scales the second grid's s_min by 1.01 in the table", "--part=table", ["agree_N"]),
}
if MUT and MUT not in MUTATIONS:
    sys.exit(f"unknown mutation {MUT}")


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    if not ok:
        FAILED.append(name)


def fin(x):
    """x, or infinity for a NaN, so that max() cannot drop a failed solve"""
    return x if x == x else float("inf")


def peak_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20  # kB on Linux


def emit(rec):
    if args.out:
        with open(args.out, "a") as fh:
            fh.write(json.dumps(rec) + "\n")


# ------------------------------------------------------------------------------------------------ profile, v_c
S_SOB = 3 * (np.pi / 2) ** (4 / 3)  # Sobolev's S of mass_bounds.md item 5, = FL's S_3
VC58 = {2.0: 3.74490, 4.0: 2.47097, 8.0: 1.56243, 16.0: 0.98503}  # mass_bounds.md item 5
BOX = {2.0: 8.0, 4.0: 5.0}  # #58's boxes


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


# ------------------------------------------------------------------------------------------------ FFT operators
G0D = np.array([1, 1, -1, -1])[:, None, None, None]  # g0 = i diag(1, 1, -1, -1)


class Grid:
    """Periodic box [-L, L)^3, N points per side; fields are (nc, N, N, N) complex arrays, nc = 2 or 4."""

    def __init__(self, L, N, nc):
        self.L, self.N, self.nc = L, N, nc
        self.x1 = -L + 2 * L / N * np.arange(N)
        k1 = 2 * np.pi * np.fft.fftfreq(N, d=2 * L / N)
        self.k = (k1[:, None, None], k1[None, :, None], k1[None, None, :])
        k2 = self.k[0] ** 2 + self.k[1] ** 2 + self.k[2] ** 2
        k2[0, 0, 0] = 1.0
        self.ik2 = 1 / k2
        self.ik2[0, 0, 0] = 0.0  # D^-1 sets k = 0 to zero
        self.shape = (nc, N, N, N)
        self.n = nc * N**3

    def coords(self):
        return np.meshgrid(self.x1, self.x1, self.x1, indexing="ij")

    def radius(self):
        X = self.coords()
        r = np.sqrt(X[0] ** 2 + X[1] ** 2 + X[2] ** 2)
        return X, r, np.where(r == 0, 1.0, r)

    def sigk(self, u0, u1):
        k1, k2, k3 = self.k
        return k3 * u0 + (k1 - 1j * k2) * u1, (k1 + 1j * k2) * u0 - k3 * u1

    def Dinv(self, u):
        """D^-1 = (i g.k)^-1 = -i g.k/|k|^2 on 4-spinors with g^j = [[0, s_j], [s_j, 0]], -i s.k/|k|^2 on 2-spinors."""
        h = sfft.fftn(u, axes=(1, 2, 3), workers=1)
        out = np.empty_like(h)
        if self.nc == 4:
            out[0], out[1] = self.sigk(h[2], h[3])
            out[2], out[3] = self.sigk(h[0], h[1])
        else:
            out[0], out[1] = self.sigk(h[0], h[1])
        out *= -1j * self.ik2
        return sfft.ifftn(out, axes=(1, 2, 3), workers=1)


class Op:
    """A(v) = I + v M with M = B D^-1, B a pointwise matrix with pointwise adjoint BH. T = M^H = -D^-1 BH."""

    def __init__(self, grid, B, BH):
        self.g, self.B, self.BH, self.napply = grid, B, BH, 0
        self.shape, self.n = grid.shape, grid.n

    def M(self, e):
        self.napply += 1
        return self.B(self.g.Dinv(e))

    def T(self, e):
        self.napply += 1
        return -self.g.Dinv(self.BH(e))

    def lin(self, fn):
        sh = self.shape
        return LinearOperator((self.n, self.n), matvec=lambda x: fn(x.reshape(sh)).ravel(), dtype=complex)


HALF = 1.0 if MUT == "half" else 0.5


def witten4(g, tau1):
    def B(e):
        return (HALF * 1j) * tau1 * G0D * e
    return Op(g, B, lambda e: -B(e))  # g0 is anti-Hermitian


def witten2(g, tau1, sgn=1):
    s = sgn * HALF * (1.1 if MUT == "bscale" else 1.0)

    def B(e):
        return -s * tau1 * e
    return Op(g, B, B)  # real multiplication, Hermitian


def ly_potential(g, scale=1.0):
    """Loss-Yau A = 3 (1 + r^2)^-2 ((1 - r^2) e1 + 2 x1 x + 2 Sigma x), Sigma x = (0, -z, y) (FL Section 1, d = 3)."""
    (x, y, z), r, _ = g.radius()
    f = 3 * scale / (1 + r**2) ** 2
    return f * (1 - r**2 + 2 * x * x), f * (2 * x * y - 2 * z), f * (2 * x * z + 2 * y)


def ly_op(g):
    A1, A2, A3 = ly_potential(g, 1.1 if MUT == "lyscale" else 1.0)

    def sA(e):  # sigma.A e
        return np.stack([A3 * e[0] + (A1 - 1j * A2) * e[1], (A1 + 1j * A2) * e[0] - A3 * e[1]])

    def B(e):
        return -1j * sA(e)
    BH = B if MUT == "adj" else (lambda e: 1j * sA(e))  # -i sigma.A is anti-Hermitian: B^H = -B
    return Op(g, B, BH)


def tanh_tau1(g, sig):
    X, r, rr = g.radius()
    return -tanh_profile(sig)(rr)[1] * X[0] / rr


def planted_tau1(g):
    _, r, _ = g.radius()
    return 6 / (1 + r**2)


def make_op(kind, L, N, sig=2.0):
    """The operator, its unit of v (v_c, and 1/sqrt(3) for the planted cases whose kappa_K(1) is sqrt(3))."""
    if kind == "witten4":
        g = Grid(L, N, 4)
        return witten4(g, tanh_tau1(g, sig)), 1 / kappa1(sig)
    g = Grid(L, N, 2)
    if kind == "witten2":
        return witten2(g, tanh_tau1(g, sig)), 1 / kappa1(sig)
    if kind == "planted":
        return witten2(g, planted_tau1(g)), 1 / math.sqrt(3)
    return ly_op(g), 1 / math.sqrt(3)


# ------------------------------------------------------------------------------------------------ s_min and the scan
def smin(W, v, free=False, tol=1e-6, ncv=24, seed=1):
    """sigma_min(I + v M) by Lanczos on the Gram operator, on the full space or (free) on mean-free fields as #63."""
    vc = v if MUT == "conj" else np.conj(v)
    ax = (1, 2, 3)

    def H(e):
        if not free:
            Ae = e + v * W.M(e)
            return Ae + vc * W.T(Ae)
        c = e.mean(axis=ax, keepdims=True)
        Ah = (e - c) + v * W.M(e - c)
        out = Ah + vc * W.T(Ah)
        return out - out.mean(axis=ax, keepdims=True) + c

    v0 = np.random.default_rng(seed).standard_normal(W.n) + 0j
    n0 = W.napply
    LANCZOS["runs"] += 1
    ok = True
    try:
        s2 = eigsh(W.lin(H), k=1, which="SA", ncv=ncv, tol=tol, maxiter=2 if MUT == "lowiter" else 400, v0=v0,
                   return_eigenvectors=False)
    except ArpackNoConvergence as e:
        s2, ok = e.eigenvalues, False
        LANCZOS["unconverged"] += 1
        if len(s2) == 0:
            return float("nan"), W.napply - n0, False
    return math.sqrt(max(s2[0], 0.0)), W.napply - n0, ok


def op_norm(W):
    """||M|| = sqrt(lambda_max(M^H M)), from below by Lanczos; the 2 % margin is the certificate's safety."""
    H = W.lin(lambda e: W.T(W.M(e)))
    try:
        r = eigsh(H, k=1, which="LA", ncv=24, tol=1e-5, maxiter=300, return_eigenvectors=False)
    except ArpackNoConvergence as e:
        r = e.eigenvalues if len(e.eigenvalues) else np.array([np.nan])
    return 1.02 * math.sqrt(abs(r[0]))


def refine(W, direction, lo, hi, tol=1e-4, maxit=40):
    """Golden-section minimum of s_min(direction * t) on [lo, hi] (v in absolute units)."""
    gr = (math.sqrt(5) - 1) / 2
    a, b = lo, hi
    c, d = b - gr * (b - a), a + gr * (b - a)
    fc, fd = smin(W, direction * c)[0], smin(W, direction * d)[0]
    for _ in range(maxit):
        if (b - a) < tol * max(abs(b), 1e-9) or min(fc, fd) < 1e-7:
            break
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = smin(W, direction * c)[0]
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = smin(W, direction * d)[0]
    t, s = (c, fc) if fc < fd else (d, fd)
    return t, s


def confirm_zero(kind, direction, v0, configs, sig=2.0):
    """Refine the dip near |v| = v0 at each (L, N); a zero mode of the continuum operator gives a grid zero whose v*
    settles as the box grows and the grid refines. Returns (rows, converged); rows are (L, N, v*, s_min at v*)."""
    rows = []
    for L, N in configs:
        W, _ = make_op(kind, L, N, sig)
        t, s = refine(W, direction, 0.9 * v0, 1.1 * v0)
        rows.append((L, N, t, s))
        print(f"    (L, N) = ({L:g}, {N}): v* = {t:.5f}, s_min there {s:.1e}", flush=True)
    ts = [r_[2] for r_ in rows]
    spread = (max(ts) - min(ts)) / np.mean(ts)
    ok = all(r_[3] < 1e-3 for r_ in rows) and spread < 0.02
    print(f"    spread of v* {spread:.1e} (limit 2e-2), all s_min < 1e-3: {all(r_[3] < 1e-3 for r_ in rows)}")
    return rows, ok


def distinct(cands, tol=0.02):
    out = []
    for c in sorted(cands):
        if not out or abs(c[0] - out[-1][0]) > tol * abs(c[0]):
            out.append(c)
        elif c[1] < out[-1][1]:
            out[-1] = c
    return out


def scan(W, unit, tmax, step, direction=1.0, tag="", log=print):
    """s_min along v = direction * t * unit, t = step, 2 step, ..., tmax. Returns the records and the interval verdicts."""
    lip = op_norm(W) * (0.1 if MUT == "lip" else 50.0 if MUT == "lipbig" else 1.0)
    log(f"  ||M|| = {lip / 1.02:.5f}, ||M|| v_c = {lip / 1.02 * unit:.5f} (the certificate uses 1.02 ||M||)")
    ts = list(np.round(np.arange(step, tmax + 1e-9, step), 10))
    recs = {0.0: (1.0, 0, True)}  # s_min(0) = 1
    for t in ts:
        t0 = time.time()
        s, na, ok = smin(W, direction * t * unit)
        recs[t] = (s, na, ok)
        row = {"kind": tag, "t": float(t), "v": float(t * unit), "smin": s, "napply": na, "ok": ok,
               "secs": round(time.time() - t0, 1), "peak_gb": round(peak_gb(), 3)}
        emit({"rec": "sample", **row})
        log(f"  t = v/v_c = {t:6.3f}  v = {t * unit:9.4f}  s_min = {s:.6f}  [{na} apps, {row['secs']} s, "
            f"peak {row['peak_gb']:.2f} GB]")
    verdicts, cands = [], []

    def resolve(ta, sa, tb, sb, depth):
        lb = (sa + sb - lip * (tb - ta) * unit) / 2
        if lb > 0:
            verdicts.append((ta, tb, lb, "certified"))
            return
        if depth >= 4 and len(cands) >= 6:  # a budget of refinements, so an uncertifiable scan ends
            verdicts.append((ta, tb, lb, "unresolved"))
            return
        if depth < 4:
            tm = (ta + tb) / 2
            sm = smin(W, direction * tm * unit)[0]
            emit({"rec": "bisect", "kind": tag, "t": float(tm), "smin": sm})
            log(f"    bisect [{ta:.4f}, {tb:.4f}] at {tm:.4f}: s_min = {sm:.6f}")
            resolve(ta, sa, tm, sm, depth + 1)
            resolve(tm, sm, tb, sb, depth + 1)
            return
        tt, ss = refine(W, direction, ta * unit, tb * unit)
        verdicts.append((ta, tb, lb, "candidate" if ss < 1e-3 else "dip"))
        cands.append((tt, ss, ta, tb))
        log(f"    refined [{ta:.4f}, {tb:.4f}]: minimum s_min = {ss:.3e} at v = {tt:.6f} ({tt / unit:.5f} units)")

    for a, b in zip([0.0] + ts[:-1], ts):
        resolve(a, recs[a][0], b, recs[b][0], 0)
    nc = sum(1 for v_ in verdicts if v_[3] == "candidate")
    unresolved = sum(1 for v_ in verdicts if v_[3] != "certified")
    smin_all = min((recs[t][0], t) for t in ts)
    emit({"rec": "summary", "kind": tag, "lip": lip / 1.02, "unit": unit, "min_s": smin_all[0], "min_t": smin_all[1],
          "intervals": len(verdicts), "unresolved": unresolved, "candidates": nc, "peak_gb": round(peak_gb(), 3),
          "napply": W.napply})
    return recs, verdicts, cands, lip


# ------------------------------------------------------------------------------------------------ adjoint, dense
def part_adjoint():
    print("== adjoint: T = M^H for the three operators (dot test with random fields)")
    rng = np.random.default_rng(5)
    worst = 0.0
    for kind in ("witten4", "witten2", "ly"):
        W, _ = make_op(kind, 6.0, 16)
        x = rng.standard_normal(W.shape) + 1j * rng.standard_normal(W.shape)
        y = rng.standard_normal(W.shape) + 1j * rng.standard_normal(W.shape)
        err = abs(np.vdot(x, W.M(y)) - np.vdot(W.T(x), y)) / (np.linalg.norm(x) * np.linalg.norm(W.M(y)))
        print(f"  {kind:8s} |<x, My> - <Tx, y>| / (|x| |My|) = {err:.1e}")
        worst = max(worst, err)
    check("adjoint", worst < 1e-12, f"worst {worst:.1e}")


PAULI = [np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], dtype=complex)]


def dense_smin(L, N, nc, tau1, v):
    """Independent s_min: explicit DFT matrices, explicit gammas, pseudo-inverse of D, dense SVD. Returns (full, free)."""
    n3 = N**3
    F = np.fft.fft(np.eye(N), axis=0)
    F3 = np.kron(np.kron(F, F), F)
    F3i = np.linalg.inv(F3)
    k1 = 2 * np.pi * np.fft.fftfreq(N, d=2 * L / N)
    kk = [m.ravel() for m in np.meshgrid(k1, k1, k1, indexing="ij")]
    d = [F3i @ np.diag(1j * kj) @ F3 for kj in kk]
    if nc == 2:
        gam = PAULI
    else:
        Z = np.zeros((2, 2))
        gam = [np.block([[Z, s], [s, Z]]) for s in PAULI]
    D = sum(np.kron(gam[j], d[j]) for j in range(3))
    tau = tau1.ravel()
    if nc == 4:
        B = np.kron(0.5 * 1j * np.diag([1, 1, -1, -1]), np.diag(tau))
    else:
        B = np.kron(-0.5 * np.eye(2), np.diag(tau))
    A = np.eye(nc * n3) + v * (B @ np.linalg.pinv(D, rcond=1e-9))
    full = np.linalg.svd(A, compute_uv=False)[-1]
    P0 = np.kron(np.eye(nc), np.eye(n3) - np.ones((n3, n3)) / n3)
    sv = np.linalg.svd(A @ P0, compute_uv=False)  # nc zero values for the constant fields
    return full, np.sort(sv)[nc]


def part_dense():
    print("== dense: Lanczos s_min against a dense SVD built from explicit DFT matrices and gammas (N = 8)")
    L, N = 6.0, 8
    worst = 0.0
    for kind, nc in (("witten2", 2), ("witten4", 4)):
        g = Grid(L, N, nc)
        X, r, _ = g.radius()
        tau1 = (X[0] + 0.3 * X[1] - 0.2 * X[2] + 0.5 * X[0] * X[1]) * np.exp(-r**2 / 10)
        W = witten2(g, tau1) if nc == 2 else witten4(g, tau1)
        for v in (3.0, 0.7 + 0.4j, -2.0):
            full, free = dense_smin(L, N, nc, tau1, v)
            sf = smin(W, v, tol=1e-10)[0]
            sr = smin(W, v, free=True, tol=1e-10)[0]
            print(f"  {kind:8s} v = {v}: full dense {full:.9f} Lanczos {sf:.9f}; free dense {free:.9f} Lanczos {sr:.9f}")
            worst = max(worst, fin(abs(full - sf)), fin(abs(free - sr)))
    check("dense", worst < 1e-7, f"largest |dense - Lanczos| {worst:.1e}")


def part_vc():
    print(f"== vc: v_c = 1/kappa_K(1), S = {S_SOB:.9f}")
    worst = 0.0
    for sig, ref in VC58.items():
        vc = 1 / kappa1(sig)
        worst = max(worst, abs(vc - ref))
        print(f"  sigma R = {sig:4.1f}: v_c = {vc:.5f} (#58: {ref:.5f})")
    check("vc", worst < 6e-6, f"equal to #58's printed values, max difference {worst:.1e}")


# #63's frank_loss_output.txt, --part=scan, L = 8 (sigma R = 2) and 5 (sigma R = 4), N = 64, mean-free s_min
REF63 = {(2.0, 1.0): 0.691615, (2.0, 2.0): 0.531532, (2.0, 4.0): 0.389574,
         (4.0, 1.0): 0.729964, (4.0, 2.0): 0.580202, (4.0, 4.0): 0.440720}


def part_regress63():
    print("== regress63: the mean-free s_min at N = 64 equals #63's scan table, and the full-space value is below it")
    worst = 0.0
    for (sig, t), ref in REF63.items():
        W, vc = make_op("witten2", BOX[sig], 64, sig)
        sr = smin(W, t * vc, free=True)[0]
        sf = smin(W, t * vc)[0]
        worst = max(worst, fin(abs(sr - ref)))
        print(f"  sigma R = {sig:g}, v/v_c = {t:g}: mean-free {sr:.6f} (#63: {ref:.6f}), full space {sf:.6f}")
    check("regress63", worst < 2e-6, f"largest difference from #63's table {worst:.1e}")


def part_sectors():
    print("== sectors: s_min(witten4) = s_min(witten2, +v) = s_min(witten2, -v) (sigma R = 2, L = 8, N = 64)")
    W4, vc = make_op("witten4", 8.0, 64)
    W2, _ = make_op("witten2", 8.0, 64)
    worst = 0.0
    for t in (2.0, 4.0):
        s4 = smin(W4, t * vc)[0]
        sp, sm = smin(W2, t * vc)[0], smin(W2, -t * vc)[0]
        print(f"  v/v_c = {t:g}: witten4 {s4:.8f}, witten2 +v {sp:.8f}, witten2 -v {sm:.8f}")
        worst = max(worst, fin(abs(s4 - sp)), fin(abs(s4 - sm)))
    check("sectors", worst < 1e-7, f"largest difference {worst:.1e}")


def part_lanczos():
    print("== lanczos: s_min is the same from three random starting vectors, and every run converged")
    W, vc = make_op("witten2", 8.0, 64)
    worst = 0.0
    for t in (3.0, 7.0, 10.0):
        ss = [smin(W, t * vc, seed=sd)[0] for sd in (1, 2, 3)]
        print(f"  v/v_c = {t:g}: " + ", ".join(f"{s:.9f}" for s in ss))
        worst = max(worst, fin(max(ss) - min(ss)))
    check("lanczos", worst < 1e-8 and LANCZOS["unconverged"] == 0,
          f"seed spread {worst:.1e}, unconverged runs {LANCZOS['unconverged']} of {LANCZOS['runs']}")


# ------------------------------------------------------------------------------------------------ planted
def part_planted():
    L, N = args.L or 10.0, args.N or 64
    print(f"== planted: tau1 = 6/(1 + r^2), zero modes at v = +-i; L = {L}, N = {N}")
    W, unit = make_op("planted", L, N)
    t0 = time.time()
    lam = eigs(W.lin(W.T), k=8, which="LM", ncv=32, tol=1e-7, return_eigenvectors=False)
    lam = lam[np.argsort(-np.abs(lam))]  # eigenvalues of M^H, conjugate to M's; the set is symmetric
    print("  top |lambda| = " + ", ".join(f"{abs(z):.4f}" for z in lam) + f"  [{time.time() - t0:.0f} s]")
    rho, im = abs(lam[0]), max(abs(z.real) / abs(z) for z in lam[:4])
    vstar = [-1 / next(z for z in lam if z.imag > 0), -1 / next(z for z in lam if z.imag < 0)]
    sv = [smin(W, vs)[0] for vs in vstar]
    spm = [smin(W, 1j)[0], smin(W, -1j)[0]]
    print(f"  zero-mode speeds v* = -1/lambda: " + ", ".join(f"{z.real:+.4f}{z.imag:+.4f}i" for z in vstar)
          + "; s_min there " + ", ".join(f"{s:.1e}" for s in sv) + "; s_min(+i), s_min(-i) = "
          + ", ".join(f"{s:.4f}" for s in spm))
    tt, ss = refine(W, 1j, 0.9, 1.1)
    print(f"  refine along +i t: minimum s_min = {ss:.1e} at t = {tt:.5f}")
    ok = abs(rho - 1) < 0.03 and im < 1e-3 and max(sv) < 1e-3 and max(spm) < 0.05 and abs(tt - 1) < 0.03 and ss < 1e-3
    check("planted_pm_i", ok, f"|lambda| = {rho:.4f} (exact 1), |Re|/|lambda| = {im:.1e}, s_min at v* <= {max(sv):.1e}, "
          f"s_min(+-i) <= {max(spm):.4f}, refined |v*| = {tt:.4f}")
    if args.confirm:
        print("  the zero-mode speed 1/|lambda| against box and grid (eigenvalue route):")
        vs = []
        for Lc, Nc in ((10.0, 64), (14.0, 64), (20.0, 64)):
            Wc, _ = make_op("planted", Lc, Nc)
            lc = eigs(Wc.lin(Wc.T), k=8, which="LM", ncv=32, tol=1e-7, return_eigenvectors=False)
            vs.append(1 / max(abs(lc)))
            print(f"    (L, N) = ({Lc:g}, {Nc}): 1/|lambda| = {vs[-1]:.5f}", flush=True)
        check("planted_conv", vs[0] > vs[1] > vs[2] > 1 and vs[2] - 1 < 0.003,
              f"|v*| = {vs[0]:.4f}, {vs[1]:.4f}, {vs[2]:.4f} at L = 10, 14, 20, tending to the exact 1")
    # the real scan of the same operator, through the same code path; rot scans the imaginary axis instead
    direction = 1j if MUT == "rot" else 1.0
    tmax = min(args.tmax, 10.0)
    print(f"  real scan of the planted tau up to {tmax:g} v_c, v_c = 1/sqrt(3):")
    recs, verdicts, cands, _ = scan(W, unit, tmax, args.step, direction, "planted")
    ts = [t for t in recs if t > 0]
    sm_, tm_ = min((recs[t][0], t) for t in ts)
    check("planted_real", not cands and all(v_[3] == "certified" for v_ in verdicts),
          f"no zero-mode candidate, every one of {len(verdicts)} intervals certified, smallest s_min {sm_:.4f} at "
          f"{tm_:g} v_c")


def part_ly():
    L, N = args.L or 10.0, args.N or 32
    print(f"== ly: Loss-Yau potential, a real zero mode at v = 1 (to test the real scan); L = {L}, N = {N}")
    W, unit = make_op("ly", L, N)
    recs, verdicts, cands, _ = scan(W, unit, 3.0, args.step, 1.0, "ly")
    near = [c for c in cands if abs(c[0] - 1) < 0.05]
    ok = bool(near) and all(c[1] < 1e-3 for c in near)
    if near:
        print(f"  refined zero-mode speed v* = {near[0][0]:.5f} with s_min = {near[0][1]:.1e}")
    check("ly_real", ok, f"{len(distinct(cands))} distinct candidate(s); the one within 5 % of the exact v = 1 has "
          "s_min < 1e-3")
    if args.confirm and near:
        print("  box and grid convergence of the zero near v = 1:")
        rows, conv = confirm_zero("ly", 1.0, near[0][0], [(10.0, 32), (10.0, 48), (14.0, 48), (14.0, 64)])
        toward = abs(rows[-1][2] - 1) < abs(rows[0][2] - 1) and abs(rows[-1][2] - 1) < 0.01
        check("ly_converged", conv and toward, f"v* = {rows[0][2]:.5f} at L = 10 and {rows[-1][2]:.5f} at L = 14, "
              "settling toward the exact 1; every s_min below 1e-3")


# ------------------------------------------------------------------------------------------------ real scan, table
def part_scan():
    sig = args.sigma
    L = args.L or BOX[sig]
    N = args.N or 64
    kind = "witten2" if args.spinor == 2 else "witten4"
    W, vc = make_op(kind, L, N, sig)
    print(f"== scan: {kind}, sigma R = {sig}, L = {L}, N = {N}, v_c = {vc:.5f}, up to {args.tmax} v_c in steps of "
          f"{args.step}")
    tag = f"{kind}-s{sig:g}-L{L:g}-N{N}"
    emit({"rec": "start", "kind": tag, "sigma": sig, "L": L, "N": N, "spinor": args.spinor, "vc": vc,
          "tmax": args.tmax, "step": args.step})
    recs, verdicts, cands, lip = scan(W, vc, args.tmax, args.step, 1.0, tag)
    ts = [t for t in recs if t > 0]
    sm_, tm_ = min((recs[t][0], t) for t in ts)
    print(f"  smallest s_min {sm_:.6f} at v/v_c = {tm_:g}; certificate: {sum(v_[3] == 'certified' for v_ in verdicts)} of "
          f"{len(verdicts)} intervals certified, smallest lower bound {min(v_[2] for v_ in verdicts):.4f}")
    for c in distinct([c_ for c_ in cands if c_[1] < 1e-3]):
        print(f"  CANDIDATE: s_min = {c[1]:.3e} at v = {c[0]:.6f} = {c[0] / vc:.4f} v_c; box and grid test:")
        rows, conv = confirm_zero(kind, 1.0, c[0], GRIDS[sig], sig)
        print(f"  CANDIDATE at {c[0] / vc:.4f} v_c: " + ("converges as a zero mode" if conv else "does not converge"))
    nres = sum(v_[3] == "unresolved" for v_ in verdicts)  # certified, or bisected and refined to a located minimum
    check("scan_resolved", nres == 0, f"{sum(v_[3] == 'certified' for v_ in verdicts)} of {len(verdicts)} intervals "
          f"certified, {sum(v_[3] in ('candidate', 'dip') for v_ in verdicts)} refined to a minimum, {nres} unresolved")
    check("lanczos_converged", LANCZOS["unconverged"] == 0, f"{LANCZOS['unconverged']} of {LANCZOS['runs']} runs")


GRIDS = {2.0: [(8.0, 64), (8.0, 96), (16.0, 96), (16.0, 128)], 4.0: [(5.0, 64), (5.0, 96), (10.0, 96), (10.0, 128)]}


def load_raw(path):
    out = {}
    for fn in sorted(glob.glob(os.path.join(path, "*.jsonl"))):
        for ln in open(fn):
            rec = json.loads(ln)
            if rec.get("rec") == "start":
                cur = out.setdefault((rec["sigma"], rec["L"], rec["N"]), {"start": rec, "s": {}, "sum": None, "file": fn})
            elif rec.get("rec") == "sample":
                cur["s"][round(rec["t"], 6)] = rec
            elif rec.get("rec") == "summary":
                cur["sum"] = rec
    return out


def part_table():
    print(f"== table: s_min(v) from the scan runs in {args.raw}")
    data = load_raw(args.raw)
    worst, missing, cand = 0.0, [], 0
    for sig in (2.0, 4.0):
        cols = GRIDS[sig]
        have = [c for c in cols if (sig,) + c in data]
        missing += [(sig,) + c for c in cols if c not in have]
        if not have:
            continue
        ts = sorted(set.intersection(*[set(data[(sig,) + c]["s"]) for c in have]))
        vc = data[(sig,) + have[0]]["start"]["vc"]
        print(f"-- sigma R = {sig:g}, v_c = {vc:.5f}; columns are (L, N), dx = 2L/N; extrap is 2 s(2L1, N2) - s(L1, N1)")
        names = [f"L={c[0]:g},N={c[1]}" for c in have]
        print("  v/v_c | " + " | ".join(f"{n:>11s}" for n in names) + " |  extrap  | s sqrt(v/v_c)")
        s = {}
        for c in have:
            pert = 1.01 if MUT == "gridperturb" and c == cols[1] else 1.0
            s[c] = {t: data[(sig,) + c]["s"][t]["smin"] * pert for t in ts}
        for t in ts:
            ext = ""
            if cols[0] in have and cols[3] in have:  # equal dx
                ext = f"{(cols[3][0] * s[cols[3]][t] - cols[0][0] * s[cols[0]][t]) / (cols[3][0] - cols[0][0]):.6f}"
            print(f"  {t:5.2f} | " + " | ".join(f"{s[c][t]:11.6f}" for c in have) + f" | {ext:>8s} | "
                  f"{s[have[-1]][t] * math.sqrt(t):.4f}")
        for box in sorted({c[0] for c in have}):
            ns = [c for c in have if c[0] == box]
            if len(ns) == 2:
                d = max(abs(s[ns[0]][t] - s[ns[1]][t]) for t in ts)
                worst = max(worst, d)
                print(f"  grid change at L = {box:g}, N = {ns[0][1]} to {ns[1][1]}: max |ds| = {d:.1e}")
        for c in have:
            sm = data[(sig,) + c]["sum"]
            if sm is None:
                missing.append((sig,) + c)
                continue
            cand += sm["candidates"]
            print(f"  (L, N) = ({c[0]:g}, {c[1]}): ||M|| v_c = {sm['lip'] * vc:.5f}, smallest s_min {sm['min_s']:.6f} at "
                  f"{sm['min_t']:g} v_c, intervals {sm['intervals']}, not certified {sm['unresolved']}, candidates "
                  f"{sm['candidates']}, peak {sm['peak_gb']:.2f} GB, {sm['napply']} applications")
    check("complete", not missing, "all eight runs present with their summaries" if not missing else f"missing {missing}")
    check("agree_N", worst < 2e-4, f"largest change in s_min between the two grids of one box {worst:.1e}")
    check("zero_candidates", cand == 0, f"{cand} zero-mode candidate(s) in the scans")


def run_controls():
    print("== controls: each --mutate makes its target checks fail; the same run without it passes them")
    bad, base = [], {}
    for name, (what, opts, targets) in MUTATIONS.items():
        if args.only and name not in args.only.split(","):
            continue
        if name == "gridperturb" and not args.raw:
            print(f"  {name:11s} skipped, needs --raw")
            continue
        extra = [f"--raw={args.raw}"] if args.raw else []
        res = {}
        for mut in ("", name):
            if mut == "" and opts in base:
                res[mut] = base[opts]
                continue
            cmd = [sys.executable, __file__] + ([f"--mutate={mut}"] if mut else []) + opts.split() + extra
            t0 = time.time()
            p = subprocess.run(cmd, capture_output=True, text=True)
            fl = [ln.split("] ")[1].split(":")[0] for ln in p.stdout.splitlines() if ln.strip().startswith("[FAIL]")]
            res[mut] = (p.returncode, fl, time.time() - t0)
            if mut == "":
                base[opts] = res[mut]
        hit = lambda fl: [t_ for t_ in targets if t_ in fl]  # noqa: E731
        fired = res[name][0] == 1 and len(hit(res[name][1])) == len(targets) and not hit(res[""][1])
        print(f"  {name:11s} ({what}; {opts}): baseline exit {res[''][0]}, mutated exit {res[name][0]}, failed "
              f"{res[name][1]}, {res[name][2]:.0f} s -> {'fired' if fired else 'DID NOT FIRE'}", flush=True)
        if not fired:
            bad.append(name)
    check("controls", not bad, "every control fired" if not bad else f"not fired: {bad}")


if args.controls:
    run_controls()
else:
    for name, fn in (("adjoint", part_adjoint), ("dense", part_dense), ("vc", part_vc), ("regress63", part_regress63),
                     ("sectors", part_sectors), ("lanczos", part_lanczos), ("planted", part_planted), ("ly", part_ly),
                     ("scan", part_scan), ("table", part_table)):
        if name in PARTS:
            fn()
print(f"peak resident memory {peak_gb():.2f} GB")
if not args.controls:
    check("memory", peak_gb() < MEM_LIMIT, f"peak {peak_gb():.2f} GB, limit {MEM_LIMIT} GB")
print("FAILED: " + ", ".join(FAILED) if FAILED else "all checks passed")
sys.exit(1 if FAILED else 0)

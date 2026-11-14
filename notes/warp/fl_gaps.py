#!/usr/bin/env python3
"""Checks for the two Frank-Loss gaps that #68 found, issue #74 (notes/warp/fl_gaps.md).

G1. FL22 Remark 4: does |sigma.grad psi| <= |A||psi| with psi in L^3, A in L^3 give psi in H^1-dot? The regularity proof of
    FL21 is run with M = g psi^+/|psi|^2 in place of sigma.A. The checks in section G1 evaluate each operator-level step of
    that proof with a general complex 2x2 matrix field, and find every place where the coefficient is required to be
    Hermitian.
G2. FL22 Proposition 7, Step 3 takes f = chi |psi|_eps^(-2/(d-1)) and says grad f is in L^d, which fails for general psi in
    H^1-dot. The identity is extended to H^1-dot by density. Section G2 evaluates both sides of Proposition 7, of its
    Step 3 and of its Step 4 along smooth psi_n -> psi in H^1-dot, for psi that are not smooth.

    python3 notes/warp/fl_gaps.py                  checks, then every control
    python3 notes/warp/fl_gaps.py --no-controls    checks only
    python3 notes/warp/fl_gaps.py --mutate=NAME    one control (its target checks must fail, exit 1)
    python3 notes/warp/fl_gaps.py --only=a,b       selected sections, no controls (names in SECTIONS)

Exit status is 1 if a check fails or a control does not fire. One thread, Normal priority.

Shares no code with review_fl_check.py (#68). That script evaluates Proposition 7 on a Cartesian grid for smooth spinors;
this one integrates in polar coordinates, radially in t = ln r (so the singularity at the origin is resolved) and over the
sphere by a product Gauss rule, with analytic gradients of the spinor fields. The conventions (Pauli matrices, <u, v> =
u^+ v, Prop. 7 at d = 3) are those of frank_loss.md and FL22.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import platform
import subprocess
import sys
import time
from fractions import Fraction

import mpmath as mp
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
except Exception:
    pass


def peak_mb():
    """Peak resident memory of this process in MB (reported, never acted on)."""
    try:
        if os.name == "nt":
            import ctypes
            from ctypes import wintypes

            class PMC(ctypes.Structure):
                _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                            ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                            ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]

            pmc = PMC()
            pmc.cb = ctypes.sizeof(PMC)
            k32 = ctypes.windll.kernel32
            k32.GetCurrentProcess.restype = wintypes.HANDLE
            ctypes.windll.psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(PMC), wintypes.DWORD]
            ctypes.windll.psapi.GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb)
            return pmc.PeakWorkingSetSize / 2**20
        import resource
        r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return r / 2**20 if sys.platform == "darwin" else r / 2**10
    except Exception:
        return float("nan")


OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate", "")
# name: (what it breaks, checks that must fail, sections it runs)
CONTROLS = {
    # ---- G1
    "mtranspose": ("M = g psi^T/|psi|^2 (transpose, no conjugate) in place of g psi^+/|psi|^2",
                   ["m_matrix", "integral_equation"], ["m_matrix", "integral_equation"]),
    "adjoint_drop": ("chi_eps := G (Gamma*phi) (G taken Hermitian) in place of G^+ (Gamma*phi)",
                     ["duality", "audit"], ["duality", "audit"]),
    "hermpart": ("bounds |G^+ v| by the norm of the Hermitian part of G", ["norm_bound"], ["norm_bound"]),
    "gamma_nonherm": ("replaces sigma_1 by the non-Hermitian sigma_1 + 0.4 i sigma_3", ["kernel", "duality", "h1_step"],
                      ["kernel", "duality", "h1_step"]),
    "hls_exp": ("uses 1/s = 1/r + 1/2 in place of 1/r + 1/3", ["exponents"], ["exponents"]),
    "gamma_i": ("drops the factor i of the Green's function Gamma = (i/4 pi) sigma.x/|x|^3", ["integral_equation"],
                ["integral_equation"]),
    # ---- G2
    "prop7coef": ("uses 4.1 in place of 4 in the bracket of Proposition 7",
                  ["planted_extremal", "prop7_smooth", "prop7_nonsmooth"], ["planted", "prop7"]),
    "step3coef": ("uses 1.05 in place of 1 on the cross term of Step 3",
                  ["step3_smooth", "step3_nonsmooth"], ["prop7"]),
    "dropt5": ("drops the cut-off term (sum over j != k of ... d_j chi) of Step 3", ["step3_smooth", "step3_nonsmooth"], ["prop7"]),
    "weakonly": ("approximants converge to psi in L^p and weakly, not in H^1-dot (an added oscillation of fixed energy)",
                 ["prop7_limit"], ["prop7"]),
    "smoothf": ("evaluates grad f for the smooth psi_delta in place of the singular psi", ["gradf_Ld"], ["gradf"]),
    # planted_extremal compares both sides of an identity that is 0 = 0 for the extremal, so a uniform scale error cannot show there
    "quadbug": ("weights of the sphere rule off by 1 %", ["planted_norms"], ["planted"]),
    "coarse": ("4 Gauss nodes per radial panel in place of 16", ["quad_converged"], ["planted"]),
}
if MUT and MUT not in CONTROLS:
    sys.exit(f"unknown mutation {MUT}")
RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


# ---------------------------------------------------------------------------------------------------------------------
# Common: Pauli matrices (sigma_1 mutated by the gamma_nonherm control), the spinor inner product <u, v> = u^+ v.
# ---------------------------------------------------------------------------------------------------------------------
S1 = np.array([[0, 1], [1, 0]], complex)
S2 = np.array([[0, -1j], [1j, 0]], complex)
S3 = np.array([[1, 0], [0, -1]], complex)
if MUT == "gamma_nonherm":
    S1 = S1 + 0.4j * S3
SIG = np.array([S1, S2, S3])
I2 = np.eye(2, dtype=complex)
HERM_TRUE = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], complex)  # the true Pauli matrices


def opnorm(m):
    return np.linalg.norm(m, 2)


def sdot(v):
    """sigma . v for a real 3-vector v (or an array (..., 3))."""
    return np.einsum("...j,jab->...ab", v, SIG)


def cx(rng, *shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def greens(z, imag=True):
    """Gamma(z) = (i / 4 pi) sigma.z / |z|^3, the kernel of the inverse of sigma.(-i grad), FL21 section 2."""
    n = np.linalg.norm(z, axis=-1)
    pref = (1j if (imag and MUT != "gamma_i") else 1.0) / (4 * np.pi)
    return pref * sdot(z) / n[..., None, None] ** 3


# ---------------------------------------------------------------------------------------------------------------------
# Fields psi(x) = P(r) a + Q(r) sum_j x_j b_j with analytic gradients, evaluated at pairs (r, xhat).
# ---------------------------------------------------------------------------------------------------------------------
class Field:
    """psi = P(r) a + Q(r) sum_j x_j b_j, a in C^2, b in C^(3 x 2). d_k psi = xhat_k (P' a + r Q' B) + Q b_k, B = sum xhat_j b_j."""

    def __init__(self, P, dP, Q, dQ, a, b):
        self.P, self.dP, self.Q, self.dQ, self.a, self.b = P, dP, Q, dQ, np.asarray(a, complex), np.asarray(b, complex)

    def eval(self, r, xh):
        B = xh @ self.b  # (n, 2)
        Pv, dPv, Qv, dQv = self.P(r), self.dP(r), self.Q(r), self.dQ(r)
        psi = Pv[:, None] * self.a[None, :] + (Qv * r)[:, None] * B
        inner = dPv[:, None] * self.a[None, :] + (dQv * r)[:, None] * B
        dpsi = xh.T[:, :, None] * inner[None, :, :] + Qv[None, :, None] * self.b[:, None, :]
        return psi, dpsi  # (n, 2), (3, n, 2)


def theta(s):
    """Smooth cut-off, 1 for s <= 1 and 0 for s >= 2."""
    s = np.asarray(s, float)
    h = lambda x: np.where(x > 0, np.exp(-1.0 / np.where(x > 0, x, 1.0)), 0.0)
    num, den = h(2 - s), h(2 - s) + h(s - 1)
    return np.where(s <= 1, 1.0, np.where(s >= 2, 0.0, num / np.where(den > 0, den, 1.0)))


def dtheta(s):
    s = np.asarray(s, float)
    h = lambda x: np.where(x > 0, np.exp(-1.0 / np.where(x > 0, x, 1.0)), 0.0)
    dh = lambda x: np.where(x > 0, h(x) / np.where(x > 0, x, 1.0) ** 2, 0.0)
    den = h(2 - s) + h(s - 1)
    val = (-dh(2 - s) * h(s - 1) - h(2 - s) * dh(s - 1)) / np.where(den > 0, den, 1.0) ** 2
    return np.where((s <= 1) | (s >= 2), 0.0, val)


def profile_u(delta, R):
    """u(r) = (r^2 + delta^2)^(-1/8) (1 + r^2)^(-3/2) theta(r/R); delta = 0 and R = inf is |x|^(-1/4)(1 + |x|^2)^(-3/2)."""
    def u(r):
        c = theta(r / R) if np.isfinite(R) else 1.0
        return (r * r + delta * delta) ** -0.125 * (1 + r * r) ** -1.5 * c

    def du(r):
        c = theta(r / R) if np.isfinite(R) else 1.0
        dc = dtheta(r / R) / R if np.isfinite(R) else 0.0
        base = (r * r + delta * delta) ** -0.125 * (1 + r * r) ** -1.5
        return base * dc + base * c * (-0.25 * r / (r * r + delta * delta) - 3 * r / (1 + r * r))
    return u, du


ZERO = lambda r: np.zeros_like(r)
PHI0 = np.array([0.6, 0.8j])


def make_case(name, delta, R):
    """The four test spinors. T1 is the issue's psi = |x|^(-1/4)(1 + |x|^2)^(-3/2) phi_0 (delta > 0 and R < inf: smooth, compact support)."""
    if name in ("T1", "T2", "T3"):
        u, du = profile_u(delta, R)
        if name == "T1":
            return Field(u, du, ZERO, ZERO, PHI0, np.zeros((3, 2)))
        if name == "T2":  # (1 + i s sigma.x) phi_0 times the profile, s = +1: Frank-Loss's extremal direction structure
            return Field(u, du, u, du, PHI0, np.array([1j * HERM_TRUE[j] @ PHI0 for j in range(3)]))
        rng = np.random.default_rng(74)
        return Field(u, du, u, du, cx(rng, 2), 0.8 * cx(rng, 3, 2))
    # T4: w(r)(a + sum xhat_j b_j), w = (1 + r^2)^(-3/2); degree-0 homogeneous at the origin, grad psi in L^2 and not in L^3
    rng = np.random.default_rng(7474)
    a, b = cx(rng, 2), 0.8 * cx(rng, 3, 2)
    cut = (lambda r: theta(r / R)) if np.isfinite(R) else (lambda r: 1.0)
    dcut = (lambda r: dtheta(r / R) / R) if np.isfinite(R) else (lambda r: 0.0)
    w = lambda r: (1 + r * r) ** -1.5 * cut(r)
    dw = lambda r: -3 * r * (1 + r * r) ** -2.5 * cut(r) + (1 + r * r) ** -1.5 * dcut(r)
    rho = lambda r: np.sqrt(r * r + delta * delta)
    Q = lambda r: w(r) / rho(r)
    dQ = lambda r: dw(r) / rho(r) - w(r) * r / rho(r) ** 3
    return Field(w, dw, Q, dQ, a, b)


# ---------------------------------------------------------------------------------------------------------------------
# Quadrature: sphere by a product Gauss rule, radius in t = ln r (or in r) by composite Gauss-Legendre with a refinement check.
# ---------------------------------------------------------------------------------------------------------------------
def sphere_rule(nmu=16, nph=24):
    x, w = np.polynomial.legendre.leggauss(nmu)
    ph = 2 * np.pi * (np.arange(nph) + 0.5) / nph
    C, PH = np.meshgrid(x, ph, indexing="ij")
    st = np.sqrt(1 - C * C)
    xh = np.stack([st * np.cos(PH), st * np.sin(PH), C], -1).reshape(-1, 3)
    W = (w[:, None] * np.full(nph, 2 * np.pi / nph)[None, :]).reshape(-1)
    if MUT == "quadbug":
        W = W * 1.01
    return xh, W


XH, WS = sphere_rule()


def sphere_integral(fun, r, xh=None, ws=None):
    """int dOmega fun(r, xhat) for radii r (n,), fun(rr, xx) -> (m, K) on pairs. Returns (n, K)."""
    xh = XH if xh is None else xh
    ws = WS if ws is None else ws
    na = len(ws)
    out = []
    for c0 in range(0, len(r), 48):
        rc = r[c0:c0 + 48]
        rr = np.repeat(rc, na)
        xx = np.tile(xh, (len(rc), 1))
        v = fun(rr, xx)
        out.append((v.reshape(len(rc), na, -1) * ws[None, :, None]).sum(1))
    return np.vstack(out)


M_NODES = 4 if MUT == "coarse" else 16  # Gauss nodes per radial panel; the coarse control uses 4


def radial_quad(fun, segs, m=None, refine=False, rule=None):
    """Integrate int_0^inf dr r^2 int dOmega fun. segs = [(kind, lo, hi, npan)] with kind 't' (r = e^t, weight r^3) or
    'r' (weight r^2). Returns (I, err): I on the panels given, err = |I(2x panels) - I(x panels)| per component when
    refine is set (zeros otherwise). rule = (xhat, weights) of the sphere rule, default the module's."""
    xg, wg = np.polynomial.legendre.leggauss(M_NODES if m is None else m)
    xh_, ws_ = rule if rule is not None else (XH, WS)

    def one(factor):
        tot = 0.0
        for kind, lo, hi, npan in segs:
            e = np.linspace(lo, hi, npan * factor + 1)
            a, b = e[:-1], e[1:]
            x = ((a + b) / 2)[:, None] + ((b - a) / 2)[:, None] * xg[None, :]
            w = ((b - a) / 2)[:, None] * wg[None, :]
            x, w = x.ravel(), w.ravel()
            r = np.exp(x) if kind == "t" else x
            wt = w * (r**3 if kind == "t" else r**2)
            tot = tot + (wt[:, None] * sphere_integral(fun, r, xh_, ws_)).sum(0)
        return tot

    i1 = one(1)
    if not refine:
        return i1, np.zeros_like(i1)
    i2 = one(2)
    return i2, np.abs(i2 - i1)


# ---------------------------------------------------------------------------------------------------------------------
# G1. Regularity of psi with sigma.(-i grad) psi = M psi, M = g psi^+/|psi|^2.
# ---------------------------------------------------------------------------------------------------------------------
def section_m_matrix():
    print("\n== G1.1 The matrix M = g psi^+/|psi|^2: M psi = g, |M| = |g|/|psi| <= |A|, for complex psi and g ==")
    rng = np.random.default_rng(1)
    worst_app, worst_norm, worst_herm = 0.0, 0.0, 0.0
    herm_defect = []
    for _ in range(2000):
        psi, g = cx(rng, 2), cx(rng, 2)
        if MUT == "mtranspose":
            M = np.outer(g, psi) / np.vdot(psi, psi).real  # transpose, not conjugate transpose
        else:
            M = np.outer(g, psi.conj()) / np.vdot(psi, psi).real  # g psi^+/|psi|^2
        worst_app = max(worst_app, np.linalg.norm(M @ psi - g) / np.linalg.norm(g))
        worst_norm = max(worst_norm, abs(opnorm(M) - np.linalg.norm(g) / np.linalg.norm(psi)) / (np.linalg.norm(g) / np.linalg.norm(psi)))
        herm_defect.append(np.linalg.norm(M - M.conj().T) / np.linalg.norm(M))
    # the zero set: psi = 0 forces g = 0 a.e. since |g| <= |A||psi|, and M := 0 there gives M psi = g = 0
    check("m_matrix", worst_app < 1e-13 and worst_norm < 1e-13 and np.median(herm_defect) > 0.3,
          "M psi = g to %.1e and |M|_op = |g|/|psi| to %.1e over 2000 complex samples; M is not Hermitian (median "
          "|M - M^+|/|M| = %.2f), so the proof cannot use M = M^+" % (worst_app, worst_norm, np.median(herm_defect)))


def section_exponents():
    print("\n== G1.2 The exponent bookkeeping of FL21's regularity proof (d = 3), exact rational arithmetic ==")
    third = Fraction(1, 3)
    half_or_third = Fraction(1, 2) if MUT == "hls_exp" else third
    bad = []
    n = 0
    ps = [Fraction(k, 8) for k in range(13, 97)]  # 13/8 ... 12, includes 3/2 < p < infinity samples
    rs = [Fraction(k, 8) for k in range(13, 97)]
    for p in ps:
        for r in rs:
            if p < 3 and r > 3 * p / (3 - p):  # admissible r of the proof
                continue
            n += 1
            s = 1 / (1 / r + half_or_third)  # 1/s = 1/r + 1/3
            pp, rp = p / (p - 1), r / (r - 1)  # dual exponents
            t = 1 / (1 / pp - third)  # Gamma*: L^{p'} -> L^t, kernel |x|^-2 in R^3 (HLS: 1/q* = 1/q - 1/3)
            sp = 1 / (1 / rp - third)  # Gamma*: L^{r'} -> L^{s'}
            ok = (s > 1 and s <= p                 # F_eps psi in L^s from L^p with F bounded, compact support
                  and 1 < pp < 3 and 1 < rp < 3    # HLS for Gamma* on L^{p'} and on L^{r'}
                  and 1 / s + 1 / sp == 1          # s' is the dual of s
                  and third + 1 / t == 1 / pp      # Holder: chi in L^{p'} from G in L^3 and Gamma*phi in L^t
                  and third + 1 / sp == 1 / rp     # Holder: chi in L^{r'} from G in L^3 and Gamma*phi in L^{s'}
                  and 1 / r == 1 / s - third)      # HLS for Gamma*(F psi): L^s -> L^r
            if not ok:
                bad.append((p, r))
    # the cases used: p = 3 (Remark 4), then r = 6, then g in L^{3/2} and L^2 by Holder with A in L^3
    p3 = all(Fraction(1, 3) + 1 / p == Fraction(1, 2) for p in (Fraction(6),)) and Fraction(1, 3) + Fraction(1, 3) == Fraction(2, 3)
    check("exponents", not bad and p3 and n > 5000,
          "%d admissible (p, r) pairs, p and r in (3/2, 12]: s > 1, s <= p, both HLS pairs admissible, the duals consistent; "
          "violations %d. At p = 3, r = 6 are admissible, and 1/2 = 1/3 + 1/6, 2/3 = 1/3 + 1/3 give g in L^2 and L^(3/2)"
          % (n, len(bad)))


def section_kernel():
    print("\n== G1.3 Hermiticity of the Pauli matrices: where it enters (kernel and symbol) ==")
    rng = np.random.default_rng(2)
    w_adj, w_nrm, w_sym = 0.0, 0.0, 0.0
    for _ in range(500):
        z = rng.normal(size=3)
        k = rng.normal(size=3)
        v = cx(rng, 2)
        G = greens(z)
        w_adj = max(w_adj, opnorm(G.conj().T - greens(-z)))  # Gamma(z)^+ = Gamma(-z), so Gamma* is self-adjoint
        w_nrm = max(w_nrm, abs(opnorm(G) * 4 * np.pi * np.dot(z, z) - 1))  # |Gamma(z)| = 1/(4 pi |z|^2) for HLS
        w_sym = max(w_sym, abs(np.linalg.norm(sdot(k) @ v) - np.linalg.norm(k) * np.linalg.norm(v)))  # |sigma.k v| = |k||v|
    check("kernel", max(w_adj, w_nrm, w_sym) < 1e-12,
          "Gamma(z)^+ = Gamma(-z) (%.1e), |Gamma(z)| = 1/(4 pi |z|^2) (%.1e), |sigma.k v| = |k||v| (%.1e): three uses of "
          "sigma_j^+ = sigma_j, none of them involves the coefficient" % (w_adj, w_nrm, w_sym))


def duality_defect(G, rng, hermitian_formula=False):
    """|<phi, Gamma*(G psi)> - <chi, psi>| / scale on a point cloud with weights, chi = G^+ (Gamma* phi), or G (Gamma* phi)."""
    n = 40
    x = rng.normal(size=(n, 3))
    w = rng.uniform(0.2, 1.0, n)
    phi, psi = cx(rng, n, 2), cx(rng, n, 2)
    diff = x[:, None, :] - x[None, :, :]
    K = np.zeros((n, n, 2, 2), complex)
    for i in range(n):
        for j in range(n):
            if i != j:
                K[i, j] = greens(diff[i, j])
    gam = lambda h: np.einsum("ijab,j,jb->ia", K, w, h)  # (Gamma* h)_i = sum_j w_j Gamma(x_i - x_j) h_j
    Gpsi = np.einsum("iab,ib->ia", G, psi)
    lhs = np.sum(w * np.einsum("ia,ia->i", phi.conj(), gam(Gpsi)))
    Gt = G if hermitian_formula else np.conj(np.transpose(G, (0, 2, 1)))
    chi = np.einsum("iab,ib->ia", Gt, gam(phi))
    rhs = np.sum(w * np.einsum("ia,ia->i", chi.conj(), psi))
    return abs(lhs - rhs) / (abs(lhs) + 1e-30)


def section_duality():
    print("\n== G1.4 The duality step: int <phi, Gamma*(G psi)> = int <G^+(Gamma* phi), psi>, G a general complex matrix field ==")
    rng = np.random.default_rng(3)
    n = 40
    Ggen = cx(rng, n, 2, 2)
    Gher = Ggen + np.conj(np.transpose(Ggen, (0, 2, 1)))
    drop = MUT == "adjoint_drop"
    d_gen = duality_defect(Ggen, rng, drop)
    d_her = duality_defect(Gher, rng, drop)
    d_wrong = duality_defect(Ggen, rng, True)  # the formula for Hermitian G, applied to a general G
    check("duality", d_gen < 1e-12 and d_her < 1e-12,
          "chi := %s: relative defect %.1e for a general G, %.1e for Hermitian G (point cloud, 40 nodes, random weights). "
          "With G in place of G^+ the general G gives %.1e: the adjoint is what the step needs"
          % ("G (Gamma* phi)" if drop else "G^+ (Gamma* phi)", d_gen, d_her, d_wrong))


def section_norm_bound():
    print("\n== G1.5 The norm bound |chi| <= |G||Gamma* phi| with chi = G^+ v ==")
    rng = np.random.default_rng(4)
    worst, sup = -1.0, 0.0
    for _ in range(3000):
        G = cx(rng, 2, 2)
        v = cx(rng, 2)
        if MUT == "hermpart":
            bound = opnorm((G + G.conj().T) / 2) * np.linalg.norm(v)  # |G| read as the norm of its Hermitian part
        else:
            bound = opnorm(G) * np.linalg.norm(v)
        lhs = np.linalg.norm(G.conj().T @ v)
        worst = max(worst, lhs / bound - 1)
        u, s, vh = np.linalg.svd(G)
        sup = max(sup, np.linalg.norm(G.conj().T @ u[:, 0]) / (opnorm(G) * 1.0))  # attained at the top left singular vector
    check("norm_bound", worst <= 1e-12 and abs(sup - 1) < 1e-12,
          "|G^+ v| <= |G|_op |v| over 3000 complex matrices (worst ratio - 1 = %.1e), attained at the top singular vector "
          "(%.15f). |G^+| = |G| for the operator norm, so only |G_eps| <= |M| <= |A| enters" % (worst, sup))


def smooth_field(kind):
    """psi decaying like r^-2: kind 'extremal' is Frank-Loss's (1 + r^2)^(-3/2)(1 + i sigma.x) phi_0, for which g = Lambda psi,
    Lambda = 3/(1 + r^2), so M = Lambda psi psi^+/|psi|^2 is Hermitian; kind 'generic' has random a, b_j, so g is not parallel
    to psi and M is a non-Hermitian rank-one matrix."""
    u = lambda r: (1 + r * r) ** -1.5
    du = lambda r: -3 * r * (1 + r * r) ** -2.5
    if kind == "extremal":
        return Field(u, du, u, du, PHI0, np.array([1j * HERM_TRUE[j] @ PHI0 for j in range(3)]))
    rng = np.random.default_rng(5)
    return Field(u, du, u, du, cx(rng, 2), 0.8 * cx(rng, 3, 2))


def section_integral_equation():
    print("\n== G1.6 The integral equation psi = Gamma*(M psi), planted: psi smooth with g = sigma.(-i grad)psi known ==")
    x0 = np.array([0.3, -0.2, 0.5])
    # polar coordinates about x0: (Gamma*h)(x0) = (i/4 pi) int dOmega int_0^inf drho sigma.omega h(x0 - rho omega)
    ug, wu = np.polynomial.legendre.leggauss(200)
    uu, wu = (ug + 1) / 2, wu / 2
    rho, drho = uu / (1 - uu), 1 / (1 - uu) ** 2
    pref = (1.0 if MUT == "gamma_i" else 1j) / (4 * np.pi)
    out = []
    for label, kind in (("extremal, g = Lambda psi", "extremal"), ("generic, M not Hermitian", "generic")):
        case = smooth_field(kind)
        acc = np.zeros(2, complex)
        defect = 0.0
        for om, wo in zip(XH, WS):
            y = x0[None, :] - rho[:, None] * om[None, :]
            r = np.linalg.norm(y, axis=1)
            psi, dpsi = case.eval(r, y / r[:, None])
            g = -1j * np.einsum("jab,jnb->na", SIG, dpsi)  # g = sigma.(-i grad) psi
            nrm = np.einsum("na,na->n", psi.conj(), psi).real
            M = np.einsum("na,nb->nab", g, psi if MUT == "mtranspose" else psi.conj()) / nrm[:, None, None]  # g psi^+/|psi|^2
            Mpsi = np.einsum("nab,nb->na", M, psi)
            acc += wo * np.einsum("ab,r,rb->a", pref * sdot(om), wu * drho, Mpsi)
            defect = max(defect, np.linalg.norm(M - np.conj(np.transpose(M, (0, 2, 1))), axis=(1, 2)).max())
        r0 = np.linalg.norm(x0)
        psi0, _ = case.eval(np.array([r0]), (x0 / r0)[None, :])
        out.append((label, np.linalg.norm(acc - psi0[0]) / np.linalg.norm(psi0[0]), defect))
    check("integral_equation", all(e < 1e-8 for _, e, _ in out) and out[1][2] > 0.1,
          "(Gamma*(M psi))(x0) = psi(x0) at x0 = (0.3, -0.2, 0.5), relative error: %s. For the generic spinor M is "
          "visibly non-Hermitian (max |M - M^+| = %.2f); psi ~ r^-2 is in L^p for every p > 3/2"
          % ("; ".join("%s %.1e" % (lab, e) for lab, e, _ in out), out[1][2]))


def section_h1_step():
    print("\n== G1.7 The last step: psi in L^6 and g in L^2 give grad psi in L^2, since |sigma.k v| = |k||v| ==")
    rows = []
    for name, delta, R in (("T3", 0.0, np.inf), ("T4", 0.0, np.inf), ("T3", 1e-2, 8.0)):
        f = make_case(name, delta, R)

        def comps(r, xh):
            psi, dpsi = f.eval(r, xh)
            sd = np.einsum("jab,jnb->na", SIG, dpsi)
            return np.stack([np.sum(np.abs(dpsi) ** 2, axis=(0, 2)), np.sum(np.abs(sd) ** 2, axis=1)], axis=1)

        v, err = radial_quad(comps, segs_for(delta, R), rule=sphere_rule(8, 8))  # polynomial in xhat: exact
        rows.append((name, delta, v[0], v[1], abs(v[0] - v[1]) / v[0]))
    check("h1_step", all(r[4] < 1e-9 for r in rows),
          "int |sigma.grad psi|^2 = int |grad psi|^2 (name, delta, int|grad psi|^2, int|sigma.grad psi|^2, relative difference): %s"
          % "; ".join("%s, %g, %.10f, %.10f, %.1e" % r for r in rows))


def section_audit():
    print("\n== G1.8 Audit: which steps of the regularity proof need the coefficient to be Hermitian ==")
    rng = np.random.default_rng(8)
    n = 40
    Gg = cx(rng, n, 2, 2)  # a general complex matrix field
    drop = MUT == "adjoint_drop"
    steps = []  # (step, passes for a general matrix field, passes for a general field when the Hermitian-case formula is used)
    # eq. (decomp): int <phi, Gamma*(G psi)> = int <chi, psi>, chi = G^+ (Gamma* phi) (as proved), or G (Gamma* phi) (FL21, G Hermitian)
    steps.append(("duality identity, eq. (decomp)", duality_defect(Gg, rng, drop) < 1e-12, duality_defect(Gg, rng, True) < 1e-12))
    # |chi| <= |G||Gamma* phi| for chi = G^+ v or G v: operator norms of G and G^+ agree
    w_adj = max(np.linalg.norm(G.conj().T @ v) / (opnorm(G) * np.linalg.norm(v)) - 1 for G, v in ((cx(rng, 2, 2), cx(rng, 2)) for _ in range(500)))
    w_her = max(np.linalg.norm(G @ v) / (opnorm(G) * np.linalg.norm(v)) - 1 for G, v in ((cx(rng, 2, 2), cx(rng, 2)) for _ in range(500)))
    steps.append(("pointwise bound |chi| <= |G||Gamma* phi|", w_adj <= 1e-12, w_her <= 1e-12))
    # the matrix M and its split M = F + G: no adjoint, no Hermiticity (G1.1 tests M psi = g and |M| = |g|/|psi|)
    psi, g = cx(rng, 2), cx(rng, 2)
    M = np.outer(g, psi.conj()) / np.vdot(psi, psi).real
    keep = rng.random(50) < 0.5  # F = M on a set E, G = M off it: |G| <= |M| pointwise, F + G = M
    Fs, Gs = np.where(keep[:, None, None], M, 0), np.where(keep[:, None, None], 0, M)
    ok = np.allclose(Fs + Gs, M) and max(opnorm(x) for x in Gs) <= opnorm(M) + 1e-15
    steps.append(("split M = F + G with |G| <= |M| (F bounded, compact support)", ok, ok))
    # properties of the Pauli matrices only (G1.3): |Gamma(x)| = 1/(4 pi |x|^2), Gamma(x)^+ = Gamma(-x), |sigma.k v| = |k||v|
    z = rng.normal(size=3)
    ok = abs(opnorm(greens(z)) * 4 * np.pi * (z @ z) - 1) < 1e-12 and opnorm(greens(z).conj().T - greens(-z)) < 1e-12
    steps.append(("kernel norm, Gamma(x)^+ = Gamma(-x), symbol norm (gammas only)", ok, ok))
    need_more = [s[0] for s in steps if not s[1]]
    needs_adjoint = [s[0] for s in steps if s[1] and not s[2]]
    print("  %-62s %-24s %s" % ("step", "general M, adjoint used", "general M, Hermitian-case formula"))
    for nm, gen, her in steps:
        print("  %-62s %-24s %s" % (nm, "passes" if gen else "FAILS", "passes" if her else "fails"))
    ok = not need_more and needs_adjoint == ["duality identity, eq. (decomp)"]
    check("audit", ok,
          "steps that fail for a general matrix field: %s. Steps that pass only with the adjoint: %s. %s"
          % (need_more or "none", needs_adjoint or "none",
             "So the coefficient is used through its adjoint in exactly one step, the duality identity; the norm bound sees "
             "G and G^+ alike" if ok else "Expected: none fail, and the duality identity alone needs the adjoint."))


# ---------------------------------------------------------------------------------------------------------------------
# G2. Proposition 7 and its Steps 3 and 4 at d = 3, evaluated pointwise on analytic fields and integrated in polar form.
# ---------------------------------------------------------------------------------------------------------------------
COLS = ["P7L", "P7R", "C3L", "C3R", "S3L", "S3R", "gradpsi2", "sigmagradpsi2", "gs12", "sdn_over_s", "gradf3"]
BRACKET = 4.1 if MUT == "prop7coef" else 4.0
CROSS = 1.05 if MUT == "step3coef" else 1.0
RHO_CHI = 1.0  # chi(x) = theta(|x|/RHO_CHI): 1 on |x| < 1 and 0 beyond |x| = 2, a smooth radial cut-off for Step 3


def components(f, eps):
    """Integrand columns of COLS for the spinor field f and regularisation eps (d = 3, N = 2).

    P7L, P7R: the two sides of Proposition 7. S3L, S3R: the two sides of Step 3 (eq. (step3)) with chi = theta(r/RHO_CHI);
    C3L, C3R: the two sides of eq. (comp3), Step 3 with chi = 1 and no cut-off term. gradf3 = |grad f|^3 for f = |psi|_eps^(-1)
    (d = 3: -2/(d-1) = -1), the L^d norm of grad f that Step 3 as printed requires to be finite.
    """
    def fun(r, xh):
        psi, dpsi = f.eval(r, xh)
        p2 = np.einsum("na,na->n", psi.conj(), psi).real
        s = np.sqrt(p2 + eps * eps)  # |psi|_eps
        gs = np.stack([np.einsum("na,na->n", psi.conj(), dpsi[k]).real for k in range(3)]) / s  # grad |psi|_eps
        gs2 = (gs**2).sum(0)
        gpsi2 = np.sum(np.abs(dpsi) ** 2, axis=(0, 2))
        sd = np.einsum("jab,jnb->na", SIG, dpsi)  # sigma . grad psi
        sdn = np.sum(np.abs(sd) ** 2, axis=1)
        # phi = psi |psi|_eps^(-3/2); the bracket [-i d_j - (1/3) sigma_j sigma.(-i grad)] phi = -i [d_j phi - (1/3) sigma_j sigma.grad phi]
        dphi = s[None, :, None] ** -1.5 * dpsi - 1.5 * s[None, :, None] ** -2.5 * gs[:, :, None] * psi[None]
        sdphi = np.einsum("jab,jnb->na", SIG, dphi)
        lhs = s**2 * sum(np.sum(np.abs(dphi[j] - np.einsum("ab,nb->na", SIG[j], sdphi) / 3) ** 2, axis=1) for j in range(3))
        gs12 = gs2 / (4 * s)  # |grad |psi|_eps^(1/2)|^2
        rhs = (2 / 3) * sdn / s - 2 * gs12 * (BRACKET - 3 * p2 / s**2)
        sgs_psi = np.einsum("jab,jn,nb->na", SIG, gs, psi)  # sigma.(grad |psi|_eps) psi
        cross = np.einsum("na,na->n", sgs_psi.conj(), sd).real  # Re <sigma.(grad|psi|_eps) psi, sigma.grad psi>
        c3l = gpsi2 / s
        c3r = 4 * gs12 + sdn / s - CROSS * cross / s**2
        # Step 3's cut-off term: sum_{j != k} s^-1 Re <sigma_j psi, sigma_k d_k psi> d_j chi, d_j chi = chi'(r) xhat_j
        sj_psi = np.einsum("jab,nb->jna", SIG, psi)
        sk_dk = np.einsum("kab,knb->kna", SIG, dpsi)
        E = np.einsum("jna,kna->jkn", sj_psi.conj(), sk_dk).real
        off = sum(xh[:, j] * E[j, k] for j in range(3) for k in range(3) if j != k)
        chi, chip = theta(r / RHO_CHI), dtheta(r / RHO_CHI) / RHO_CHI
        t5 = chip * off / s * (0.0 if MUT == "dropt5" else 1.0)
        return np.stack([lhs, rhs, c3l, c3r, c3l * chi, c3r * chi + t5, gpsi2, sdn, gs12, sdn / s,
                         gs2**1.5 / s**6], axis=1)
    return fun


def dist_fun(fa, fb):
    """Integrand of ||grad(psi_a - psi_b)||_2^2."""
    def fun(r, xh):
        _, da = fa.eval(r, xh)
        _, db = fb.eval(r, xh)
        return np.sum(np.abs(da - db) ** 2, axis=(0, 2))[:, None]
    return fun


def segs_for(delta, R, tail=False):
    """Panels in t = ln r. Width 1 where the integrand is smooth on that scale, 0.1 across the cut-offs chi (r in [1, 2]) and
    theta(r/R) (r in [R, 2R]). With tail set the range runs on to r = e^14 past 2R (the non-smooth psi has no cut-off)."""
    lo = max(-70.0, np.log(delta) - 35.0) if delta > 0 else -70.0
    cuts = [0.0, np.log(2.0)] + ([np.log(R), np.log(2 * R)] if np.isfinite(R) else [])
    hi = 14.0 if (not np.isfinite(R) or tail) else np.log(2 * R)
    pts = sorted(set([lo] + cuts + [hi]))
    segs = []
    for a, b in zip(pts[:-1], pts[1:]):
        fine = (abs(a - 0.0) < 1e-12 and abs(b - np.log(2.0)) < 1e-12) or (np.isfinite(R) and abs(a - np.log(R)) < 1e-12)
        segs.append(("t", a, b, max(1, int(np.ceil((b - a) / (0.1 if fine else 1.0))))))
    return segs


# Angular rules: T1 and T2 have |psi| radial, so the integrands are polynomials of low degree in xhat and an 8 x 8 rule is exact;
# T3 and T4 need 32 x 48 (checked by quad_converged and by the sphere-rule comparison there).
RULE = {"T1": sphere_rule(8, 8), "T2": sphere_rule(8, 8), "T3": sphere_rule(32, 48), "T4": sphere_rule(32, 48)}
RULE_POLY = sphere_rule(8, 8)  # for |grad psi|^2-type integrands, polynomial in xhat for every case
EPS = {"T1": (0.5, 0.05), "T2": (0.5, 0.05), "T3": (0.5, 0.3), "T4": (0.5, 0.05)}
SEQ = [(1e-2, 8.0), (1e-4, 16.0), (1e-6, 32.0), (1e-8, 64.0), (1e-10, 128.0), (1e-12, 256.0)]  # (delta_m, R_m)


def evaluate(name, delta, R, eps, refine=False, field=None, segs=None):
    f = make_case(name, delta, R) if field is None else field
    v, err = radial_quad(components(f, eps), segs_for(delta, R) if segs is None else segs, rule=RULE[name], refine=refine)
    return dict(zip(COLS, v)), dict(zip(COLS, err))


def distance(name, delta, R, fa=None, segs=None):
    """||grad(psi_{delta,R} - psi)||_2 for the case `name`, psi = the non-smooth field (delta = 0, R = inf)."""
    fa = make_case(name, delta, R) if fa is None else fa
    fb = make_case(name, 0.0, np.inf)
    v, _ = radial_quad(dist_fun(fa, fb), segs_for(max(delta, 0.0), np.inf if delta == 0 else R, tail=True) if segs is None else segs,
                       rule=RULE_POLY)
    return float(np.sqrt(v[0]))


def weak_field(m):
    """control: smooth psi_{delta_m, R_m} + (1/n) e^(-r^2) cos(n r) phi_0, n = 4 * 2^m. It tends to psi in sup norm and in L^p,
    and weakly in H^1-dot, but the added term has gradient of fixed L^2 norm (about 0.99)."""
    delta, R = SEQ[m - 1]
    n = 4.0 * 2**m
    u, du = profile_u(delta, R)
    P = lambda r: u(r) + np.exp(-r * r) * np.cos(n * r) / n
    dP = lambda r: du(r) - 2 * r * np.exp(-r * r) * np.cos(n * r) / n - np.exp(-r * r) * np.sin(n * r)
    return Field(P, dP, ZERO, ZERO, PHI0, np.zeros((3, 2))), n


def weak_segs(m):
    delta, R = SEQ[m - 1]
    n = 4.0 * 2**m
    lo = np.log(delta) - 35.0
    cuts = [np.log(0.05), None, np.log(8.0), np.log(R), np.log(2 * R)]
    segs = [("t", lo, np.log(0.05), int(np.ceil((np.log(0.05) - lo)))),
            ("r", 0.05, 8.0, int(np.ceil((8.0 - 0.05) * n / np.pi))),
            ("t", np.log(8.0), np.log(R), max(1, int(np.ceil(np.log(R) - np.log(8.0))))),
            ("t", np.log(R), np.log(2 * R), 10)]
    return segs


# ---- planted truth ------------------------------------------------------------------------------------------------------
EXT_SEGS = [("t", -30.0, 40.0, 70)]  # the extremal decays like r^-2 in |grad s^(1/2)|^2: the tail beyond e^t is 4 pi e^-t


def extremal_field():
    return smooth_field("extremal")  # Frank-Loss's Psi = (1 + r^2)^(-3/2)(1 + i sigma.x) phi_0, |Psi| = 1/(1 + r^2)


def section_planted():
    print("\n== G2.0 Planted truth before real data: quadrature against closed forms, and the Frank-Loss extremal ==")
    # (a) T1: int |grad psi|^2 = 4 pi int u'(r)^2 r^2 dr with u = r^(-1/4)(1 + r^2)^(-3/2), by mpmath, 30 digits
    mp.mp.dps = 30
    du = lambda r: r ** mp.mpf(-0.25) * (1 + r * r) ** mp.mpf(-1.5) * (-1 / (4 * r) - 3 * r / (1 + r * r))
    exact_t1 = 4 * mp.pi * mp.quad(lambda r: du(r) ** 2 * r * r, [0, 1, mp.inf])
    v, _ = evaluate("T1", 0.0, np.inf, 0.5)
    e1 = abs(v["gradpsi2"] - float(exact_t1)) / float(exact_t1)
    # (b) the extremal at eps = 0 (|Psi| > 0 everywhere): int|grad s^(1/2)|^2 = 3 pi^2/4, int |sigma.grad Psi|^2/|Psi| = 9 pi^2/4,
    #     int |grad Psi|^2 = int |sigma.grad Psi|^2 = 9 pi^2/8, so both sides of Proposition 7 are (2/3)(9 pi^2/4) - 2 (3 pi^2/4)(4 - 3) = 0
    ext = extremal_field()
    v0, _ = radial_quad(components(ext, 0.0), EXT_SEGS, rule=RULE["T2"])
    exact = {"gs12": 3 * np.pi**2 / 4, "sdn_over_s": 9 * np.pi**2 / 4, "gradpsi2": 9 * np.pi**2 / 8, "sigmagradpsi2": 9 * np.pi**2 / 8}
    e2 = max(abs(v0[COLS.index(k)] - x) / x for k, x in exact.items())
    check("planted_norms", e1 < 1e-11 and e2 < 1e-11,
          "T1: int|grad psi|^2 = %.12f, mpmath %.12f (relative %.1e). Extremal at eps = 0: int|grad s^(1/2)|^2 = %.10f (3pi^2/4 = %.10f), "
          "int|sigma.grad Psi|^2/s = %.10f (9pi^2/4 = %.10f), int|grad Psi|^2 = %.10f (9pi^2/8 = %.10f), worst relative %.1e"
          % (v["gradpsi2"], float(exact_t1), e1, v0[8], exact["gs12"], v0[9], exact["sdn_over_s"], v0[6], exact["gradpsi2"], e2))
    # (c) both sides of Proposition 7 vanish for the extremal at eps = 0 (its phi is a twistor spinor) and agree at eps > 0
    v3, _ = radial_quad(components(ext, 0.3), EXT_SEGS, rule=RULE["T2"])
    ok = abs(v0[0]) < 1e-10 and abs(v0[1]) < 1e-10 and abs(v3[0] - v3[1]) < 1e-10
    check("planted_extremal", ok,
          "extremal, eps = 0: LHS = %.1e, RHS = %.1e (both 0 by the twistor equation and the integrals above); eps = 0.3: LHS = %.10f, "
          "RHS = %.10f" % (v0[0], v0[1], v3[0], v3[1]))
    # (d) the quadrature is converged: 2x panels change nothing, and a second sphere rule gives the same numbers
    worst, worst_sph = 0.0, 0.0
    keep = [i for i, c in enumerate(COLS) if c != "gradf3"]  # |grad s|^3 is not smooth where grad s = 0 (see gradf)
    big = sphere_rule(48, 64)
    for name in ("T1", "T2", "T3", "T4"):
        for eps in EPS[name]:
            vv, ee = evaluate(name, 0.0, np.inf, eps, refine=True)
            worst = max(worst, max(ee[k] / abs(vv[k]) for k in COLS if abs(vv[k]) > 0))
            if name in ("T3", "T4"):
                f = make_case(name, 0.0, np.inf)
                a, _ = radial_quad(components(f, eps), segs_for(0.0, np.inf), rule=RULE[name])
                b, _ = radial_quad(components(f, eps), segs_for(0.0, np.inf), rule=big)
                worst_sph = max(worst_sph, np.max(np.abs(a - b)[keep] / np.abs(b)[keep]))
    check("quad_converged", worst < IDENT_TOL / 10 and worst_sph < IDENT_TOL / 10,  # a tenth of the identity checks' tolerance
          "doubling the radial panels changes every integral of T1 to T4 (both values of eps, non-smooth psi) by at most %.1e "
          "relative; 32x48 against 48x64 sphere rule for T3 and T4 at both values of eps: %.1e" % (worst, worst_sph))
    # information, not a check: why T3 does not use eps = 0.05
    f = make_case("T3", 0.0, np.inf)
    a, _ = radial_quad(components(f, 0.05), segs_for(0.0, np.inf), rule=RULE["T3"])
    b, _ = radial_quad(components(f, 0.05), segs_for(0.0, np.inf), rule=big)
    print("  info: T3 at eps = 0.05, 32x48 against 48x64 sphere rule, largest relative difference %.1e (not used in any check)"
          % np.max(np.abs(a - b)[keep] / np.abs(b)[keep]))


# ---- the identity on smooth psi_n and on the non-smooth psi, and the limit along psi_n -> psi in H^1-dot --------------------
SIDES = [("P7L", "P7R"), ("S3L", "S3R"), ("C3L", "C3R")]
ALL_SIDES = ["P7L", "P7R", "S3L", "S3R", "C3L", "C3R"]
IDENT_TOL = 1e-9  # relative: quad_converged requires the quadrature to be good to a tenth of this
PROP7_DATA = {}


def section_prop7():
    print("\n== G2.1 Proposition 7, Step 3 (with chi) and Step 4 (chi = 1), along smooth psi_n -> psi in H^1-dot ==")
    cases = ("T1",) if MUT == "weakonly" else ("T1", "T2", "T3", "T4")
    rel = lambda a, b, scale: abs(a - b) / max(1.0, abs(scale))
    worst = {"prop7_smooth": 0.0, "step3_smooth": 0.0, "prop7_nonsmooth": 0.0, "step3_nonsmooth": 0.0}
    limit_rows, ok_limit = [], True
    lim_detail = []
    for name in cases:
        for eps in ((0.5,) if MUT == "weakonly" else EPS[name]):
            base, _ = evaluate(name, 0.0, np.inf, eps)
            worst["prop7_nonsmooth"] = max(worst["prop7_nonsmooth"], rel(base["P7L"], base["P7R"], base["P7L"]))
            worst["step3_nonsmooth"] = max(worst["step3_nonsmooth"], rel(base["S3L"], base["S3R"], base["S3L"]),
                                           rel(base["C3L"], base["C3R"], base["C3L"]))
            norm0 = np.sqrt(base["gradpsi2"])
            ds, gaps, rows = [], [], []
            for m, (delta, R) in enumerate(SEQ, start=1):
                if MUT == "weakonly":
                    f, nosc = weak_field(m)
                    sg = weak_segs(m)
                    v, _ = evaluate(name, delta, R, eps, field=f, segs=sg)
                    d = distance(name, delta, R, fa=f, segs=sg)
                else:
                    v, _ = evaluate(name, delta, R, eps)
                    d = distance(name, delta, R)
                worst["prop7_smooth"] = max(worst["prop7_smooth"], rel(v["P7L"], v["P7R"], v["P7L"]))
                worst["step3_smooth"] = max(worst["step3_smooth"], rel(v["S3L"], v["S3R"], v["S3L"]), rel(v["C3L"], v["C3R"], v["C3L"]))
                gap = max(rel(v[k], base[k], base[k]) for k in ALL_SIDES)
                ds.append(d)
                gaps.append(gap)
                rows.append((delta, R, d, v["P7L"], v["P7R"], gap))
            # (a) distance decreases strictly and ends small; (b) the sides converge: gaps decrease, end small, and gap/distance stays bounded
            ok_d = all(a > b for a, b in zip(ds, ds[1:])) and ds[-1] < 0.02 * norm0
            ok_g = all(a >= b * (1 - 1e-9) for a, b in zip(gaps, gaps[1:])) and gaps[-1] < 0.02 and gaps[-1] <= 0.05 * gaps[0]
            lip = max(g / d for g, d in zip(gaps, ds))
            ok_l = lip < 1.0
            ok_limit &= ok_d and ok_g and ok_l
            limit_rows.append((name, eps, base, rows, lip, ok_d, ok_g, ok_l))
    for name, eps, base, rows, lip, ok_d, ok_g, ok_l in limit_rows:
        print("  %s, eps = %g: psi (non-smooth) P7L = %.10f, P7R = %.10f, ||grad psi||_2 = %.6f; Lipschitz ratio max gap/dist = %.3f"
              % (name, eps, base["P7L"], base["P7R"], np.sqrt(base["gradpsi2"]), lip))
        print("      %8s %6s %12s %16s %16s %10s" % ("delta_n", "R_n", "||psi_n-psi||", "P7 LHS(psi_n)", "P7 RHS(psi_n)", "max gap"))
        for delta, R, d, a, b, gap in rows:
            print("      %8.0e %6.0f %12.4e %16.10f %16.10f %10.2e" % (delta, R, d, a, b, gap))
    check("prop7_smooth", worst["prop7_smooth"] < IDENT_TOL,
          "Proposition 7: LHS = RHS on every smooth psi_n (4 cases, 2 values of eps, 6 approximants each), worst relative difference %.1e"
          % worst["prop7_smooth"])
    check("prop7_nonsmooth", worst["prop7_nonsmooth"] < IDENT_TOL,
          "Proposition 7: LHS = RHS on the non-smooth psi itself, worst relative difference %.1e" % worst["prop7_nonsmooth"])
    check("step3_smooth", worst["step3_smooth"] < IDENT_TOL,
          "Step 3 (chi = theta(r), with the cut-off term) and its limit chi = 1 (comp3): LHS = RHS on every smooth psi_n, worst %.1e"
          % worst["step3_smooth"])
    check("step3_nonsmooth", worst["step3_nonsmooth"] < IDENT_TOL,
          "Step 3 and comp3: LHS = RHS on the non-smooth psi, worst %.1e" % worst["step3_nonsmooth"])
    check("prop7_limit", ok_limit,
          "psi_n -> psi in H^1-dot (distance decreasing to under 2 %% of ||grad psi||), and all six sides of Proposition 7, Step 3 "
          "and comp3 at psi_n converge to their values at psi (max relative gap decreasing, ending under 2 %%, and bounded by "
          "%.2f x the distance, a Lipschitz ratio)" % max(r[4] for r in limit_rows))


# ---- Step 3 as printed needs grad f in L^d; for psi in H^1-dot it fails ---------------------------------------------------
def section_gradf():
    print("\n== G2.2 Step 3 as printed: f = chi |psi|_eps^(-2/(d-1)) needs grad f in L^d. T4 is in H^1-dot and its grad f is not in L^3 ==")
    eps = 0.5
    smooth = MUT == "smoothf"
    f = make_case("T4", 5e-2 if smooth else 0.0, np.inf)  # control: the smooth field psi_delta, delta = 0.05
    incs = []
    for i in range(1, 9):  # decades r in [10^-(i+1), 10^-i]
        segs = [("t", np.log(10.0 ** -(i + 1)), np.log(10.0 ** -i), 8)]
        v, _ = radial_quad(components(f, eps), segs, rule=RULE["T4"])
        incs.append(v[COLS.index("gradf3")])
    base, _ = evaluate("T4", 0.0, np.inf, eps)
    inc = np.array(incs)
    ratio = inc[2:] / inc[2:].mean()  # decades 10^-3 ... 10^-9
    grows = inc[-1] > 0.5 * inc[2] and inc[-1] > 0.05
    const = np.max(np.abs(ratio - 1)) < 1e-3
    check("gradf_Ld", grows and const,
          "int_{10^-(i+1) < r < 10^-i} |grad f|^3 dx for i = 1..8: %s. Equal increments per decade mean int_{s<r<1} |grad f|^3 ~ "
          "(%.4f) ln(1/s) diverges, so grad f is not in L^3, while ||grad psi||_2^2 = %.6f is finite (psi is homogeneous of degree 0 at the origin)"
          % (["%.5f" % x for x in inc], inc[-1] / np.log(10.0), base["gradpsi2"]))


SECTIONS = {"m_matrix": section_m_matrix, "exponents": section_exponents, "kernel": section_kernel, "duality": section_duality,
            "norm_bound": section_norm_bound, "integral_equation": section_integral_equation, "h1_step": section_h1_step,
            "audit": section_audit, "planted": section_planted, "prop7": section_prop7, "gradf": section_gradf}


def main():
    t0 = time.time()
    print("# fl_gaps.py, issue #74. %s, %s, python %s, numpy %s, mpmath %s. mutation: %s"
          % (platform.node(), time.strftime("%Y-%m-%d"), sys.version.split()[0], np.__version__, mp.__version__, MUT or "none"))
    names = CONTROLS[MUT][2] if MUT else OPTS["only"].split(",") if "only" in OPTS else list(SECTIONS)
    for name in names:
        SECTIONS[name]()
    fails = [n for n, ok in RESULTS if not ok]
    print("\nchecks: %d passed, %d failed %s. %.0f s, peak resident memory %.0f MB"
          % (len(RESULTS) - len(fails), len(fails), fails, time.time() - t0, peak_mb()))
    if MUT:
        target = CONTROLS[MUT][1]
        print("control %s: targets %s %s" % (MUT, target, "FIRED" if all(n in fails for n in target) else "SILENT"))
        sys.exit(1)
    status = 1 if fails else 0
    if "no-controls" not in OPTS and "only" not in OPTS:
        print("\n== Controls: each mutation must make its target checks fail ==")
        for name, (what, target, _) in CONTROLS.items():
            p = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True, encoding="utf-8")
            failed = [ln.split("]")[1].split(":")[0].strip() for ln in p.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = p.returncode == 1 and all(n in failed for n in target)
            status |= 0 if fired else 1
            print("[%s] --mutate=%s (%s): exit %d, failed %s" % ("FIRED" if fired else "SILENT", name, what, p.returncode, failed), flush=True)
            for ln in p.stdout.splitlines():  # what the broken input produced, for the record
                if ln.startswith("[FAIL]") or (name == "weakonly" and ln.startswith("  ")):
                    print("      " + ln[:330])
        print("controls: %s" % ("all fired" if status == 0 else "SOME SILENT OR A CHECK FAILED"))
    print("total %.0f s, peak resident memory of the main process %.0f MB" % (time.time() - t0, peak_mb()))
    sys.exit(status)


if __name__ == "__main__":
    main()

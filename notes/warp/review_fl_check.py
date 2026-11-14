#!/usr/bin/env python3
"""Independent checks of #63's Frank-Loss extension (frank_loss.md sections 2-8), for issue #68.

Written for the review. It shares no code with frank_loss.py, only numpy, scipy and sympy, and it reuses none of its
discretisation: the spectrum of the planted case is computed from the radial ODE (no grid, no FFT), and the
Green's-function step of section 3 is evaluated by quadrature in polar coordinates about the evaluation point.

    python3 notes/warp/review_fl_check.py                  checks, then every control
    python3 notes/warp/review_fl_check.py --no-controls    checks only
    python3 notes/warp/review_fl_check.py --mutate=NAME    one control (its target checks must fail, exit 1)
    python3 notes/warp/review_fl_check.py --only=a,b       selected sections, no controls (names in SECTIONS)

Exit status is 1 if a check fails or a control does not fire. One thread, below-normal priority, under 1 GB.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import subprocess
import sys
import time

import numpy as np
import sympy as sp
from scipy import integrate, optimize

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def lower_priority():
    """Below-normal on Windows, nice 10 elsewhere. This machine's budget is shared with another project."""
    try:
        if os.name == "nt":
            import ctypes
            ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x4000)
        else:
            os.nice(10)
    except Exception as e:  # not fatal, but say so
        print("could not lower priority:", e)


lower_priority()

OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate", "")
# name: (what it breaks, checks that must fail, sections it runs)
CONTROLS = {
    "g0block": ("flips the lower block of gamma^0", ["reduction", "divxi"], ["reduction", "identities"]),
    "half": ("drops the 1/2 of B = (1/2) tau gamma^0", ["planted_ode", "planted_greens"], ["planted_ode", "planted_greens"]),
    "ode": ("uses (2 + k)/r in place of (1 + k)/r in the radial system", ["planted_ode"], ["planted_ode"]),
    "sobolev": ("uses 21/20 of the sharp Sobolev constant", ["sobolev_const", "fl_norm"], ["fl_pair"]),
    "greens": ("uses 1/(2 pi) in the kernel gamma.x/(4 pi |x|^3)", ["greens_kernel", "planted_greens"],
               ["identities", "planted_greens"]),
    "prop7": ("replaces 4 - 3|psi|^2/|psi|_eps^2 by 4.1 - 3|psi|^2/|psi|_eps^2", ["prop7"], ["prop7"]),
    "vcthree": ("uses 2 pi int r^2 |f'|^3 in place of pi int r^2 |f'|^3", ["v_c"], ["table"]),
    "extremal": ("drops the derivative of (1 + i sigma.x) from the extremal's gradient", ["fl_local_min"], ["local_min"]),
}
if MUT and MUT not in CONTROLS:
    sys.exit(f"unknown mutation {MUT}")
RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


S_SHARP = 3 * (np.pi / 2) ** (4 / 3)  # sharp Sobolev constant in R^3, S_3 = (3/4)|S^3|^(2/3)
S_USED = S_SHARP * (21 / 20 if MUT == "sobolev" else 1)
SIG = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]], complex),
       np.array([[1, 0], [0, -1]], complex)]
I2, Z2 = np.eye(2, dtype=complex), np.zeros((2, 2), complex)
# sign of the lower block of gamma^0 = i diag(1, 1, -1, -1); the g0block control flips it
LOWER = 1.0 if MUT == "g0block" else -1.0


def gammas4():
    """gamma^0 = i diag(1,1,-1,-1), gamma^j = [[0, sigma_j], [sigma_j, 0]] (mass_bounds.md, Conventions)."""
    g0 = 1j * np.diag([1, 1, LOWER, LOWER]).astype(complex)
    return g0, [np.block([[Z2, s], [s, Z2]]) for s in SIG]


def gammas4_sym():
    s = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
    z = sp.zeros(2)
    g0 = sp.I * sp.diag(1, 1, int(LOWER), int(LOWER))
    return g0, [sp.Matrix(sp.BlockMatrix([[z, sj], [sj, z]])) for sj in s], s


def radial(fn, a=0.0, b=np.inf, **kw):
    return integrate.quad(fn, a, b, limit=500, epsabs=0, epsrel=1e-13, **kw)[0]


# ---- 1. The 4-to-2 reduction (section 2), by random numbers, and the block form of T_1 ---------------------------
def section_reduction():
    print("\n== 1. Reduction of D~ on C^4 to two Pauli equations (section 2), random complex data ==")
    g0, gj = gammas4()
    rng = np.random.default_rng(68)
    worst, worst_inv, worst_norm, worst_blk = 0.0, 0.0, 0.0, 0.0
    cplx = lambda *s: rng.normal(size=s) + 1j * rng.normal(size=s)
    sdot = lambda grad: sum(SIG[j] @ grad[j] for j in range(3))
    for trial in range(400):
        a, b = cplx(2), cplx(2)
        da, db = cplx(3, 2), cplx(3, 2)  # d_j a, d_j b as independent numbers
        tau = rng.normal() + (1j * rng.normal() if trial % 2 else 0.0)  # real, then complex amplitude
        eps = np.concatenate([a, b])
        deps = [np.concatenate([da[j], db[j]]) for j in range(3)]
        Deps = sum(gj[j] @ deps[j] for j in range(3))
        Beps = 0.5 * tau * (g0 @ eps)
        Dt = Deps + Beps  # D~ eps, 4 numbers
        u, w = Dt[:2], Dt[2:]
        psi, chi = a + 1j * b, a - 1j * b
        dpsi = [da[j] + 1j * db[j] for j in range(3)]
        dchi = [da[j] - 1j * db[j] for j in range(3)]
        R1, R2 = sdot(dpsi) - 0.5 * tau * psi, sdot(dchi) + 0.5 * tau * chi
        worst = max(worst, np.abs(w + 1j * u - R1).max(), np.abs(w - 1j * u - R2).max())
        worst_inv = max(worst_inv, np.abs((R1 - R2) / (2j) - u).max(), np.abs((R1 + R2) / 2 - w).max())
        worst_norm = max(worst_norm, abs(psi @ psi.conj() + chi @ chi.conj() - 2 * (eps @ eps.conj())))
        # D and B_1 separately, in the (w + iu, w - iu) coordinates
        cD = np.concatenate([Deps[2:] + 1j * Deps[:2], Deps[2:] - 1j * Deps[:2]])
        cB = np.concatenate([Beps[2:] + 1j * Beps[:2], Beps[2:] - 1j * Beps[:2]])
        worst_blk = max(worst_blk, np.abs(cD - np.concatenate([sdot(dpsi), sdot(dchi)])).max(),
                        np.abs(cB - np.concatenate([-0.5 * tau * psi, 0.5 * tau * chi])).max())
    check("reduction", worst < 1e-13 and worst_inv < 1e-13,
          "w + iu = sigma.grad psi - tau psi/2 and w - iu = sigma.grad chi + tau chi/2, and the inverse map, "
          "400 samples, real and complex tau (worst %.1e, %.1e)" % (worst, worst_inv))
    check("reduction_norms", worst_norm < 1e-12, "|psi|^2 + |chi|^2 = 2|eps|^2 (worst %.1e)" % worst_norm)
    check("block_form", worst_blk < 1e-13,
          "in the coordinates (w + iu, w - iu), D = diag(S, S) and B_1 = diag(-M, +M) with S = sigma.grad and "
          "M = tau_1/2, so T_1 = D^-1 B_1 is -S^-1 M on psi and +S^-1 M on chi (worst %.1e)" % worst_blk)


# ---- 2. Clifford algebra, the identity of section 8, and the Green's function of D (section 3) -----------------------
def section_identities():
    print("\n== 2. Pointwise identity of section 8, and D E = delta for E = gamma.x/(4 pi |x|^3) (section 3) ==")
    g0, gj = gammas4()
    herm = (np.allclose(g0.conj().T, -g0) and np.allclose(g0 @ g0, -np.eye(4))
            and all(np.allclose(g.conj().T, g) for g in gj) and all(np.allclose(g0 @ g + g @ g0, 0) for g in gj)
            and all(np.allclose(gj[i] @ gj[j] + gj[j] @ gj[i], 2 * (i == j) * np.eye(4)) for i in range(3) for j in range(3)))
    check("clifford", herm and all(np.allclose((g0 @ g).conj().T, g0 @ g) for g in gj),
          "g0 anti-Hermitian with g0^2 = -1, g^j Hermitian, {g0, g^j} = 0, {g^i, g^j} = 2 delta, g0 g^j Hermitian")
    # divxi: 2 Re<eps, g0 D~ eps> = d_j xi^j - tau |eps|^2, xi^j = <eps, g0 g^j eps>; explicit smooth fields, symbolic
    # derivatives, evaluated at random points. Real parts are taken after differentiating (x, y, z are real).
    x, y, z = sp.symbols("x y z", real=True)
    X = [x, y, z]
    G0, GJ, _ = gammas4_sym()
    rng = np.random.default_rng(1)
    gauss = sp.exp(-(x**2 + 2 * y**2 + z**2) / 3)
    comp = lambda: (sp.Rational(int(rng.integers(1, 9)), 5) + x * sp.Rational(int(rng.integers(1, 9)), 7)
                    + sp.I * (sp.Rational(int(rng.integers(1, 9)), 3) - y * z * sp.Rational(int(rng.integers(1, 9)), 4))) * gauss
    eps = sp.Matrix([comp() for _ in range(4)])
    tau = sp.cos(x) * y + z**2 / (1 + x**2)  # real, either sign
    Dt = sum((GJ[j] * eps.diff(X[j]) for j in range(3)), sp.zeros(4, 1)) + tau * G0 * eps / 2
    lhs_c = 2 * (eps.H * G0 * Dt)[0]
    xi_c = [(eps.H * G0 * GJ[j] * eps)[0] for j in range(3)]
    rhs_c = sum(sp.diff(xi_c[j], X[j]) for j in range(3)) - tau * (eps.H * eps)[0]
    f_l, f_r = sp.lambdify(X, lhs_c, "numpy"), sp.lambdify(X, rhs_c, "numpy")
    pts = rng.normal(size=(40, 3))
    err = max(abs(complex(f_l(*p)).real - complex(f_r(*p)).real) for p in pts)
    scale = max(abs(complex(f_l(*p)).real) for p in pts)
    check("divxi", err < 1e-10 * max(scale, 1), "2 Re<eps, g0 D~ eps> = d_j xi^j - tau|eps|^2 at 40 points, "
          "residual %.1e of size %.1e" % (err, scale))
    # D E = delta. gamma^j gamma^k d_j(x_k/r^3) keeps only the symmetric part, i.e. d_j(x_j/r^3) = 0 off the origin;
    # the flux of x/(4 pi r^3) through a sphere is 1 (4 pi / (4 pi)); the control uses 2 pi in the kernel.
    r = sp.sqrt(x**2 + y**2 + z**2)
    Mjk = sp.Matrix(3, 3, lambda j, k: sp.diff(X[k] / r**3, X[j]))
    sym = sp.simplify(Mjk - Mjk.T) == sp.zeros(3)
    div0 = sp.simplify(sum(sp.diff(X[j] / r**3, X[j]) for j in range(3))) == 0
    flux = 4 * np.pi / (2 * np.pi if MUT == "greens" else 4 * np.pi)
    check("greens_kernel", sym and div0 and abs(flux - 1) < 1e-14,
          "d_j(x_k/r^3) is symmetric and div(x/r^3) = 0 off the origin; the flux of the kernel through a sphere is "
          "%.3f, so D E = delta (D^2 = Laplacian by the Clifford relations)" % flux)
    xv = rng.normal(size=3)
    gx = sum(gj[j] * xv[j] for j in range(3))
    check("kernel_norm", abs(np.linalg.norm(gx, 2) - np.linalg.norm(xv)) < 1e-13 and np.allclose(gx @ gx, xv @ xv * np.eye(4)),
          "the operator norm of gamma.x is |x|, so |E(x)| = 1/(4 pi |x|^2)")


# ---- 3. Sobolev constant, the Frank-Loss extremal pair, kappa_K = sqrt(3) -------------------------------------------
def section_fl_pair():
    print("\n== 3. S_3, the Frank-Loss appendix pair, ||A||_3^2 = 3 S, kappa_K(tau_1 = 6/(1+r^2)) = sqrt(3) ==")
    s_fl = 0.75 * (2 * np.pi**2) ** (2 / 3)  # S_3 = (3/4)|S^3|^(2/3), |S^3| = 2 pi^2 (FL, Introduction)
    u, du = (lambda r: (1 + r * r) ** -0.5), (lambda r: -r * (1 + r * r) ** -1.5)
    grad2 = 4 * np.pi * radial(lambda r: du(r) ** 2 * r * r)
    l6 = (4 * np.pi * radial(lambda r: u(r) ** 6 * r * r)) ** (1 / 3)
    others = []
    for fn, dfn in ((lambda r: np.exp(-r * r), lambda r: -2 * r * np.exp(-r * r)),
                    (lambda r: 1 / (1 + r * r), lambda r: -2 * r / (1 + r * r) ** 2),
                    (lambda r: (1 + r * r) ** -1.5, lambda r: -3 * r * (1 + r * r) ** -2.5)):
        others.append(4 * np.pi * radial(lambda r: dfn(r) ** 2 * r * r)
                      / (4 * np.pi * radial(lambda r: fn(r) ** 6 * r * r)) ** (1 / 3))
    check("sobolev_const", abs(grad2 / l6 / S_USED - 1) < 1e-10 and abs(s_fl / S_USED - 1) < 1e-12 and min(others) > S_USED,
          "Aubin-Talenti ratio %.12f, (3/4)|S^3|^(2/3) = %.12f, 3 (pi/2)^(4/3) used = %.12f; other trial functions give %s"
          % (grad2 / l6, s_fl, S_USED, ["%.4f" % o for o in others]))
    # Frank-Loss appendix pair, symbolic with phi0 arbitrary:
    #   Psi(x) = (1+r^2)^(-3/2) (1 + i s sigma.x) phi0,  sigma.(-i grad) Psi = Lambda Psi,  Lambda = 3 s / (1+r^2)
    x, y, z = sp.symbols("x y z", real=True)
    p = sp.symbols("p1:5", real=True)
    phi0 = sp.Matrix([p[0] + sp.I * p[1], p[2] + sp.I * p[3]])
    _, _, sig = gammas4_sym()
    r2 = x**2 + y**2 + z**2
    sx = x * sig[0] + y * sig[1] + z * sig[2]
    ok = True
    for s in (1, -1):
        Psi = (1 + r2) ** sp.Rational(-3, 2) * (sp.eye(2) + sp.I * s * sx) * phi0
        lhs = sum((-sp.I * sig[j] * sp.diff(Psi, v) for j, v in enumerate((x, y, z))), sp.zeros(2, 1))
        ok &= sp.simplify(lhs - sp.Rational(3 * s) / (1 + r2) * Psi) == sp.zeros(2, 1)
    check("fl_pair", ok,
          "sigma.(-i grad) Psi = (3 s/(1+r^2)) Psi for Psi = (1+r^2)^(-3/2)(1 + i s sigma.x) phi0, s = +-1, phi0 arbitrary: "
          "so sigma.grad Psi = (i s tau_1/2) Psi, a zero mode of the psi equation at v = i s with tau_1 = 6/(1+r^2)")
    a3 = (4 * np.pi * radial(lambda r: (3 / (1 + r * r)) ** 3 * r * r)) ** (2 / 3)
    t3 = (4 * np.pi * radial(lambda r: (6 / (1 + r * r)) ** 3 * r * r)) ** (1 / 3)
    kap = t3 / (2 * np.sqrt(S_USED))
    check("fl_norm", abs(a3 / (3 * S_USED) - 1) < 1e-10 and abs(kap - np.sqrt(3)) < 1e-9,
          "||3/(1+r^2)||_3^2 = %.10f, 3 S = %.10f; kappa_K(v = 1) = %.10f, sqrt 3 = %.10f, so sqrt3 v_c = 1"
          % (a3, 3 * S_USED, kap, np.sqrt(3)))


# ---- 4. Frank-Loss Proposition 7 (d = 3) on a test spinor, quadrature --------------------------------------------
def section_prop7():
    print("\n== 4. Frank-Loss Proposition 7 at d = 3, as quoted in section 4, for 2- and 4-component spinors ==")
    rng = np.random.default_rng(7)
    bracket = 4.1 if MUT == "prop7" else 4.0
    for label, gl in (("Pauli, N = 2", SIG), ("4x4 gammas of mass_bounds.md, N = 4", gammas4()[1])):
        n = gl[0].shape[0]
        cen = [np.array([0.5, -0.3, 0.2]), np.array([-0.6, 0.4, -0.5])]
        amp = [rng.normal(size=n) + 1j * rng.normal(size=n) for _ in cen]
        wav = [np.array([0.9, -0.4, 0.6]), np.array([-0.5, 0.8, 0.3])]
        wid = [0.9, 0.7]
        out = []
        for N, L in ((64, 6.0), (88, 7.0)):
            h = 2 * L / N
            ax = -L + h * np.arange(N)
            X = np.meshgrid(ax, ax, ax, indexing="ij")
            psi = np.zeros((n, N, N, N), complex)
            dpsi = [np.zeros((n, N, N, N), complex) for _ in range(3)]
            for c, a, q, wd in zip(cen, amp, wav, wid):
                d = [X[i] - c[i] for i in range(3)]
                g = np.exp(-(d[0] ** 2 + d[1] ** 2 + d[2] ** 2) / (2 * wd**2) + 1j * (q[0] * X[0] + q[1] * X[1] + q[2] * X[2]))
                psi += a[:, None, None, None] * g
                for i in range(3):
                    dpsi[i] += a[:, None, None, None] * ((1j * q[i] - d[i] / wd**2) * g)
            eps = 1.0
            p2 = np.sum(np.abs(psi) ** 2, axis=0)
            s = np.sqrt(p2 + eps**2)
            ds = [np.real(np.sum(psi.conj() * dpsi[i], axis=0)) / s for i in range(3)]  # d_i |psi|_eps
            mat = lambda M, v: np.einsum("ab,bxyz->axyz", M, v)
            dphi = [s ** -1.5 * dpsi[i] - 1.5 * s ** -2.5 * ds[i] * psi for i in range(3)]  # phi = psi |psi|_eps^(-3/2)
            sdphi = sum(mat(gl[k], dphi[k]) for k in range(3))
            # [-i d_j - (1/3) g_j g.(-i grad)] phi = -i [d_j phi - (1/3) g_j g.grad phi]
            lhs_pt = sum(np.sum(np.abs(dphi[j] - mat(gl[j], sdphi) / 3) ** 2, axis=0) for j in range(3)) * s**2
            sdpsi = sum(mat(gl[k], dpsi[k]) for k in range(3))
            grad_sqrt2 = sum((0.5 * s**-0.5 * ds[i]) ** 2 for i in range(3))  # |grad |psi|_eps^(1/2)|^2
            rhs_pt = (2 / 3) * np.sum(np.abs(sdpsi) ** 2, axis=0) / s - 2 * grad_sqrt2 * (bracket - 3 * p2 / s**2)
            lhs, rhs = lhs_pt.sum() * h**3, rhs_pt.sum() * h**3
            out.append((N, lhs, rhs, (lhs - rhs) / abs(lhs)))
            del X, psi, dpsi, dphi, sdphi, sdpsi
        check("prop7" if n == 2 else "prop7_4x4", abs(out[-1][3]) < 1e-8,
              "%s: LHS = %.10f, RHS = %.10f, relative difference at N = 64 and 88: %s"
              % (label, out[-1][1], out[-1][2], ["%.1e" % o[3] for o in out]))


# ---- 5. The planted case rho(T_1) = 1, from the radial ODE (no grid) ---------------------------------------------------
def radial_M(ell, k, r0=1e-3, R=1e3, rm=1.0):
    """Matching determinant at r = rm of the regular-at-0 and the decaying-at-infinity solutions, sector k >= 1.

    sigma.grad psi = mu W psi, W = 3/(1+r^2), psi = f Omega_k + g Omega_{-k}, mu = i ell, f = i F gives
        g' + (1-k) g/r = ell W F,    F' + (1+k) F/r = -ell W g.
    Regular at 0: g ~ r^(k-1), F ~ r^k. Decaying at infinity: F ~ r^(-k-1), g ~ r^(-k-2); the other branch is
    g ~ r^(k-1), which grows. In t = ln r both integrations run in the stable direction.
    """
    e = 2 * ell if MUT == "half" else ell  # control: B_1 without its 1/2 doubles the eigenvalue problem
    kc = k + 1 if MUT == "ode" else k  # control: wrong centrifugal term in the F equation

    def rhs(t, y):
        r = np.exp(t)
        w = 3 * r / (1 + r * r)
        return [(k - 1) * y[0] + e * w * y[1], -(1 + kc) * y[1] - e * w * y[0]]

    dt = complex if np.iscomplexobj(ell) else float
    y0 = np.array([r0 ** (k - 1), -3 * e * r0**k / (2 * k + 1)], dtype=dt)
    reg = integrate.solve_ivp(rhs, [np.log(r0), np.log(rm)], y0, method="DOP853", rtol=1e-12, atol=1e-30).y[:, -1]
    y1 = np.array([-3 * e / (2 * k + 1) * R ** (-k - 2), R ** (-k - 1)], dtype=dt)
    dec = integrate.solve_ivp(rhs, [np.log(R), np.log(rm)], y1, method="DOP853", rtol=1e-12, atol=1e-30).y[:, -1]
    return reg[0] * dec[1] - reg[1] * dec[0]


def winding(k, rho, npts=360):
    th = np.linspace(0, 2 * np.pi, npts + 1)
    ph = np.unwrap(np.angle(np.array([radial_M(rho * np.exp(1j * t), k) for t in th])))
    return (ph[-1] - ph[0]) / (2 * np.pi)


def section_planted_ode():
    print("\n== 5. Planted case tau_1 = 6/(1+r^2): eigenvalues of T_1 from the radial ODE, no grid ==")
    print("   T_1 is -S^-1 M on psi and +S^-1 M on chi (section 1), so lambda = i/ell for sigma.(-i grad) psi = ell W psi,")
    print("   W = tau_1/2 = 3/(1+r^2), and |lambda| = 1/|ell|. Roots ell of the matching determinant, sectors k = j + 1/2.")
    found = {}
    for k in (1, 2, 3):
        grid = np.linspace(0.05, 4.2, 168)
        vals = [radial_M(l, k) for l in grid]
        roots = [optimize.brentq(lambda l: radial_M(l, k), a, b, xtol=1e-13)
                 for a, b, va, vb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]) if va * vb < 0]
        found[k] = roots
        print("  k = %d: ell = %s" % (k, ", ".join("%.10f" % r for r in roots)))
        print("        (2k + 1 + 2m)/3 = %s" % ", ".join("%.10f" % ((2 * k + 1 + 2 * m) / 3) for m in range(len(roots))))
    # higher sectors: only the lowest root is sought, and none may lie below (2k + 1)/3
    for k in (4, 5, 6):
        top = (2 * k + 1) / 3 + 0.4
        grid = np.linspace(0.05, top, int((top - 0.05) / 0.025))
        vals = [radial_M(l, k) for l in grid]
        roots = [optimize.brentq(lambda l: radial_M(l, k), a, b, xtol=1e-13)
                 for a, b, va, vb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]) if va * vb < 0]
        found[k] = roots
        print("  k = %d: ell = %s   ((2k + 1)/3 = %.10f)" % (k, ", ".join("%.10f" % r for r in roots), (2 * k + 1) / 3))
    err = max((abs(r - (2 * k + 1 + 2 * m) / 3) for k, rs in found.items() for m, r in enumerate(rs)), default=np.inf)
    n_ok = all(len(found[k]) >= 3 for k in (1, 2, 3)) and all(len(found[k]) >= 1 for k in (4, 5, 6))
    lowest = min((rs[0] for rs in found.values() if rs), default=np.nan)
    check("planted_ode", err < 1e-9 and n_ok and abs(1 / lowest - 1) < 1e-9,
          "real roots ell_{k,m} = (2k + 1 + 2m)/3 in sectors k = 1..6 to %.1e (the Dirac eigenvalues +-(3/2 + n) of "
          "the round S^3, ell = (2/3)(3/2 + n)); smallest |ell| = %.12f, so rho(T_1) = 1/min|ell| = %.12f"
          % (err, lowest, 1 / lowest))
    if MUT:
        return
    # no root off the real axis: argument principle on circles
    w1 = [winding(1, rho) for rho in (0.95, 1.3)]
    w2 = [winding(2, rho) for rho in (1.4, 1.95)]
    print("  winding numbers of the matching determinant: k = 1, |ell| = 0.95 -> %.4f, 1.3 -> %.4f; "
          "k = 2, |ell| = 1.4 -> %.4f, 1.95 -> %.4f" % (w1[0], w1[1], w2[0], w2[1]))
    check("planted_nocomplex", abs(w1[0]) < 0.02 and abs(w1[1] - 2) < 0.02 and abs(w2[0]) < 0.02 and abs(w2[1] - 2) < 0.02,
          "argument principle: no root of the k = 1 determinant in |ell| < 0.95 and two (ell = +-1) in |ell| < 1.3; "
          "none of the k = 2 determinant in |ell| < 1.4 and two (+-5/3) in |ell| < 1.95, so no complex eigenvalue "
          "in these discs")


# ---- 6. T_1 eps = lambda eps through the Green's function, for the planted zero mode ------------------------------------
def section_planted_greens():
    print("\n== 6. T_1 eps = lambda eps for the planted zero mode, D^-1 = convolution with E = gamma.x/(4 pi|x|^3) ==")
    g0, gj = gammas4()
    gstack = np.array(gj)
    s = 1.0
    phi0 = np.array([0.6, 0.8j])
    sigs = np.array(SIG)

    def eps4(P):  # P: (n, 3) points. psi sector only (chi = 0): a = psi/2, b = psi/(2i), eps = (a, b)
        sx = np.einsum("pj,jab->pab", P, sigs)
        w = (1 + np.sum(P * P, axis=1)) ** -1.5
        psi = w[:, None] * (np.einsum("ab,b->a", np.eye(2), phi0)[None, :] + 1j * s * np.einsum("pab,b->pa", sx, phi0))
        return np.concatenate([psi / 2, psi / 2j], axis=1)

    bfac = 1.0 if MUT == "half" else 0.5  # B_1 = (1/2) tau_1 gamma^0; the control drops the 1/2

    def B1eps(P):
        tau1 = 6.0 / (1 + np.sum(P * P, axis=1))
        return bfac * tau1[:, None] * (eps4(P) @ g0.T)

    cnorm = 2 * np.pi if MUT == "greens" else 4 * np.pi
    x0 = np.array([0.3, -0.2, 0.5])
    nth, nph, nr = 24, 32, 80
    xg, wg = np.polynomial.legendre.leggauss(nth)  # cos(theta)
    ph = 2 * np.pi * np.arange(nph) / nph
    C, PH = np.meshgrid(xg, ph, indexing="ij")
    st = np.sqrt(1 - C * C)
    yh = np.stack([st * np.cos(PH), st * np.sin(PH), C], -1).reshape(-1, 3)
    Wf = (wg[:, None] * np.full(nph, 2 * np.pi / nph)[None, :]).reshape(-1)
    gy = np.einsum("pj,jab->pab", yh, gstack)
    tg, wt = np.polynomial.legendre.leggauss(nr)
    tt, wt = (tg + 1) / 2, wt / 2
    rr, drr = tt / (1 - tt), 1 / (1 - tt) ** 2
    # (E*h)(x0) = int_0^inf dr int dOmega (gamma.yhat / c) h(x0 - r yhat), c = 4 pi: the 1/r^2 of E cancels the Jacobian
    acc = np.zeros(4, complex)
    for r_, w_ in zip(rr, wt * drr):
        hv = B1eps(x0[None, :] - r_ * yh)
        acc += w_ * np.einsum("p,pab,pb->a", Wf, gy, hv)
    got = acc / cnorm
    e0 = eps4(x0[None, :])[0]
    lam = (got @ e0.conj()) / (e0 @ e0.conj())
    resid = np.abs(got - lam * e0).max() / np.abs(e0).max()
    want = 1j * s  # D~(v) eps = 0 with v = i s  <=>  T_1 eps = -eps/v = i s eps
    check("planted_greens", abs(lam - want) < 1e-8 and resid < 1e-8,
          "D^-1(B_1 eps)(x0) = E*(B_1 eps)(x0) = lambda eps(x0), lambda = %.10f%+.10fi (expected %+.0fi), residual %.1e; "
          "eps = Frank-Loss zero mode, in H^1-dot and L^3, x0 = %s" % (lam.real, lam.imag, want.imag, resid, x0))


# ---- 7. Section 8 sector by sector, and the numbers of the table ------------------------------------------------------
def section_sector():
    print("\n== 7. Section 8 sector by sector ==")
    x, y, z, v = sp.symbols("x y z v", real=True)
    X = [x, y, z]
    _, _, sig = gammas4_sym()
    rng = np.random.default_rng(3)
    gauss = sp.exp(-(x**2 + y**2 + 2 * z**2) / 3)
    psi = sp.Matrix([(1 + x * sp.Rational(int(rng.integers(1, 6)), 3) + sp.I * y) * gauss,
                     (sp.Rational(int(rng.integers(1, 6)), 4) - sp.I * z * x) * gauss])
    tau1 = sp.sin(y) * x + z / (1 + x**2 + y**2)
    sd = sum((sig[j] * psi.diff(X[j]) for j in range(3)), sp.zeros(2, 1))
    lhs_c = 2 * (psi.H * (sd - v * tau1 * psi / 2))[0]
    jv = [(psi.H * sig[j] * psi)[0] for j in range(3)]
    rhs_c = sum(sp.diff(jv[j], X[j]) for j in range(3)) - v * tau1 * (psi.H * psi)[0]
    fl, fr = sp.lambdify([x, y, z, v], lhs_c, "numpy"), sp.lambdify([x, y, z, v], rhs_c, "numpy")
    pts = rng.normal(size=(30, 3))
    err = max(abs(complex(fl(*p, 1.7)).real - complex(fr(*p, 1.7)).real) for p in pts)
    check("sector_identity", err < 1e-10,
          "for real v and tau_1: 2 Re<psi, (sigma.grad - v tau_1/2) psi> = d_j(psi^+ sigma_j psi) - v tau_1 |psi|^2, so a "
          "Pauli zero mode with psi in L^3 has v int tau_1 |psi|^2 = 0, for psi and for chi separately (residual %.1e)" % err)


def section_table():
    print("\n== 8. The v_c table of frank_loss.md with the sqrt 3 extension ==")
    note_vc = {2: 3.74490, 4: 2.47097, 8: 1.56243, 16: 0.98503}
    note_s3 = {2: 6.48636, 4: 4.27985, 8: 2.70621, 16: 1.70612}
    note_thin = {2: 6.82624, 4: 4.30026, 8: 2.70900, 16: 1.70656}
    ang = 2 * np.pi if MUT == "vcthree" else np.pi  # 4 pi <|cos theta|^3> = 4 pi / 4
    worst = 0.0
    for sR, wvc in note_vc.items():
        T = np.tanh(sR)
        f1 = lambda r: sR * (np.cosh(sR * (r + 1)) ** -2 - np.cosh(sR * (r - 1)) ** -2) / (2 * T)
        I3 = sum(radial(lambda r: abs(f1(r)) ** 3 * r * r, a, b) for a, b in ((0, 0.5), (0.5, 1), (1, 1.5), (1.5, 3), (3, 30)))
        vc = 2 * np.sqrt(S_SHARP) / (ang * I3) ** (1 / 3)
        thin = np.sqrt(3) * 2 * np.sqrt(S_SHARP) * (15 / (2 * np.pi * sR**2)) ** (1 / 3)
        worst = max(worst, abs(vc - wvc), abs(np.sqrt(3) * vc - note_s3[sR]), abs(thin - note_thin[sR]))
        print("  sigma R = %2d: v_c = %.5f (note %.5f), sqrt3 v_c = %.5f (note %.5f), thin wall x sqrt3 = %.5f (note %.5f)"
              % (sR, vc, wvc, np.sqrt(3) * vc, note_s3[sR], thin, note_thin[sR]))
    check("v_c", worst < 6e-6, "v_c, sqrt3 v_c and the thin-wall column agree with the table to its printed digits "
          "(worst difference %.1e)" % worst)


# ---- 9. FL's inequality near its extremal: ||a||_3^2 >= 3 S with a = |sigma.grad psi| / |psi| ------------------------
def section_local_min():
    print("\n== 9. Frank-Loss inequality near the extremal (Remark 4 form): ||a||_3^2 vs 3 S, a = |sigma.grad psi|/|psi| ==")
    rng = np.random.default_rng(11)
    phi0 = np.array([0.6, 0.8j])
    N, L = 96, 12.0
    h = 2 * L / N
    ax = -L + h * (np.arange(N) + 0.5)
    X = np.meshgrid(ax, ax, ax, indexing="ij")
    r2 = X[0] ** 2 + X[1] ** 2 + X[2] ** 2
    sxm = sum(SIG[j][:, :, None, None, None] * X[j] for j in range(3))  # sigma.x as (2, 2, N, N, N)
    M0 = np.eye(2)[:, :, None, None, None] + 1j * sxm  # 1 + i sigma.x
    w = (1 + r2) ** -1.5
    dw = [-3 * X[i] * (1 + r2) ** -2.5 for i in range(3)]
    M0phi = np.einsum("abxyz,b->axyz", M0, phi0)
    psi0 = w * M0phi
    keep = 0.0 if MUT == "extremal" else 1.0  # control: drop d_i(1 + i sigma.x) phi0 = i sigma_i phi0
    dpsi0 = [dw[i] * M0phi + keep * w * (1j * (SIG[i] @ phi0))[:, None, None, None] for i in range(3)]
    amp = rng.normal(size=2) + 1j * rng.normal(size=2)
    c = np.array([0.4, -0.3, 0.2])
    g = np.exp(-((X[0] - c[0]) ** 2 + (X[1] - c[1]) ** 2 + (X[2] - c[2]) ** 2) / (2 * 0.7**2))
    eta = amp[:, None, None, None] * g
    deta = [amp[:, None, None, None] * (-(X[i] - c[i]) / 0.7**2 * g) for i in range(3)]
    base = 27 * np.pi**2 / 4  # int (3/(1+r^2))^3 dx for the extremal, exact
    a0 = 3 / (1 + r2)
    res = {}
    for delta in (0.0, 0.05, -0.05, 0.1, -0.1):
        psi = psi0 + delta * eta
        dpsi = [dpsi0[i] + delta * deta[i] for i in range(3)]
        sd = sum(np.einsum("ab,bxyz->axyz", SIG[i], dpsi[i]) for i in range(3))
        a = np.sqrt(np.sum(np.abs(sd) ** 2, axis=0) / np.sum(np.abs(psi) ** 2, axis=0))
        res[delta] = (base + np.sum(a**3 - a0**3) * h**3) ** (2 / 3)  # a differs from a0 only where eta lives
    e = {d: res[d] - 3 * S_USED for d in res}
    ratio_q = e[0.1] / e[0.05] if abs(e[0.05]) > 0 else np.nan
    ratio_s = e[0.05] / e[-0.05] if abs(e[-0.05]) > 0 else np.nan
    ok = (abs(e[0.0]) < 1e-6 and all(e[d] > 0 for d in (0.05, -0.05, 0.1, -0.1))
          and 3.5 < ratio_q < 4.5 and 0.8 < ratio_s < 1.25)
    check("fl_local_min", ok,
          "||a||_3^2 - 3 S at delta = 0, +-0.05, +-0.1: %s (3 S = %.6f). Zero at the extremal, positive on both sides, "
          "quadratic: excess(0.1)/excess(0.05) = %.2f, excess(0.05)/excess(-0.05) = %.2f"
          % (["%+.2e" % e[d] for d in (0.0, 0.05, -0.05, 0.1, -0.1)], 3 * S_USED, ratio_q, ratio_s))


SECTIONS = {"reduction": section_reduction, "identities": section_identities, "fl_pair": section_fl_pair,
            "prop7": section_prop7, "planted_ode": section_planted_ode, "planted_greens": section_planted_greens,
            "sector": section_sector, "table": section_table, "local_min": section_local_min}


def main():
    t0 = time.time()
    print("review_fl_check.py, issue #68. numpy %s, scipy %s, sympy %s, python %s. mutation: %s"
          % (np.__version__, __import__("scipy").__version__, sp.__version__, sys.version.split()[0], MUT or "none"))
    names = CONTROLS[MUT][2] if MUT else OPTS["only"].split(",") if "only" in OPTS else list(SECTIONS)
    for name in names:
        SECTIONS[name]()
    fails = [n for n, ok in RESULTS if not ok]
    print("\nchecks: %d passed, %d failed %s. %.0f s" % (len(RESULTS) - len(fails), len(fails), fails, time.time() - t0))
    if MUT:
        target = CONTROLS[MUT][1]
        print("control %s: targets %s %s" % (MUT, target, "FIRED" if all(n in fails for n in target) else "SILENT"))
        sys.exit(1)
    status = 1 if fails else 0
    if "no-controls" not in OPTS and "only" not in OPTS:
        print("\n== Controls: each mutation must make its target checks fail ==")
        for name, (what, target, _) in CONTROLS.items():
            p = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True,
                               encoding="utf-8")
            failed = [ln.split("]")[1].split(":")[0].strip() for ln in p.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = p.returncode == 1 and all(n in failed for n in target)
            status |= 0 if fired else 1
            print("[%s] --mutate=%s (%s): exit %d, failed %s" % ("FIRED" if fired else "SILENT", name, what,
                                                                p.returncode, failed), flush=True)
        print("controls: %s" % ("all fired" if status == 0 else "SOME SILENT OR CHECK FAILED"))
    sys.exit(status)


if __name__ == "__main__":
    main()

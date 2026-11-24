#!/usr/bin/env python3
"""Linear e-folding time of the asymmetric Ellis family (#151, #163), against [GGS09a] Table I and [SN26] Sec. IV B.

    python3 ellis_efold.py                  all checks; exit 1 if any check fails
    python3 ellis_efold.py --gamma=0        the gamma_1 = 0 control and its supporting checks only
    python3 ellis_efold.py --section=old    only #151's checks (first row of Table I); --section=new only #163's
    python3 ellis_efold.py --mutate=NAME    break one input on purpose (the control for NAME must fire; see MUTATIONS)

Units are [GGS09a]'s rescaled ones (b = 1, gamma_0 = 0): g = -e^{2 g atan x} dt^2 + e^{-2 g atan x}[dx^2 + (1+x^2) dOmega^2],
rho(x) = int_0^x e^{-2 g atan y} dy, d/drho = e^{2 g atan x} d/dx. The unstable mode is chi e^{beta t} with
-chi_rhorho + W chi = -beta^2 chi, W = [GGS09a] Eq. (32), and the dimensionless time scale (proper time at the throat over the
areal throat radius, [GGS09a] p. 10 and [SN26] Sec. IV B) is
    T = e^{2 g atan g} / (sqrt(1 + g^2) beta).
g is gamma_1; [SN26]'s m/a = gamma_1. The criteria, tolerances and control targets are in negative-mass-raw I151-rpi5/README.md
(written before the first run); every check below is coded as worded there.

Methods for beta:
  A  shooting on the W equation in rho (scipy DOP853), W = the printed Eq. (32);
  B  sixth-order finite differences in rho (scipy eig_banded), W = -V + 2 (dPsi0/Psi0)^2 derived with sympy from the
     definitions (mpmath, 40 digits);
  C  gauge dPhi = 0: the second-order ODE for C = delta c derived here from the Ricci tensor (sympy), decaying solutions
     continued across the throat (singular point, exponents 0 and 3) on semicircles in the complex x plane. No W.

Checks: U1 U2 (units), S1 (background), L1 (linearization), W1 (printed W), P-A P-B P-C (planted truth), G0 .. G7, D1.

#163 (second row of Table I, gamma_1 = 2.5 ... 10, and the large-gamma_1 limit equation): the section "second row" below.
The legacy path above (domain centred at x = 0, rho~ = lam rho in [-30, 30]) is unchanged and still produces the first row. For
gamma_1 >= 2.5 the throat x = gamma_1 lies at rho~ = 85 ... 3e11, outside that domain; the new path puts the origin of the
variable u = (rho~ - rho~_throat)/sc, sc = max(1, gamma_1), at the throat. Checks: N0 N1a N1b N2 N3a N3b U3 LE0, P-A2 P-B2 P-C2
P-LA P-LB P-LC, R1 R2 R3, G5b G1b G2b G7b D1b, LE3 LE4 LE5 LE6; the README is negative-mass-raw I163-rpi5/README.md.
"""
import argparse
import cmath
import math
import platform
import sys
import time

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.linalg import eig_banded, eigh_tridiagonal
from scipy.optimize import brentq

MUT = None
FAILED = []
RESULTS = {}
T_START = time.time()

# mutation -> (what it changes, the checks that must fire on it); the same table as README section 5
MUTATIONS = {
    "potential": ("W -> W + 0.1/(1+x^2)^2 in A and B", ["G0", "G2", "G3", "G4", "G5", "G6", "R1"]),
    "printedW": ("the -4g/(1+x^2) Fx/F term of the printed W gets 3g",
                 ["W1", "G1", "G2", "G3", "G4", "G5", "R1", "N1a", "N1b", "N3a", "R2", "G5b", "G1b", "G2b"]),
    "matrix": ("B's centre weight -49/18 -> -2.7", ["P-B", "G1", "G0", "P-B2", "G1b", "P-LB", "LE3"]),
    "planted": ("planted exact beta = ell -> 1.01 ell", ["P-A", "P-B", "P-A2", "P-B2", "P-LA", "P-LB", "P-LC"]),
    "xmap": ("A's x_rho = e^{-2 g atan x}", ["P-A", "G1", "G2", "G3", "G4", "G5", "R1"]),
    "units": ("gamma_0 = +pi g/2 in the [SN26] -> [GGS09a] map", ["U1"]),
    "tfactor": ("T drops e^{2 g atan g}", ["U2", "G3", "G4", "G5", "G6", "U3"]),
    "count": ("G7 uses W + a second well", ["G7", "G7b", "LE4"]),
    "domain": ("A and B use L = 4; the new path uses L' = 2", ["D1", "G0", "G1", "D1b", "G5b", "R2", "P-A2", "P-B2", "G1b", "LE3"]),
    "deltaphi": ("the 2 in 2 e^{2(d-c)}/c' of method C becomes 1", ["L1", "G2", "G2b"]),
    "cpath": ("method C's semicircle radius -> 0", ["P-C", "G2", "P-C2", "G2b"]),
    "background": ("kappa Phi1^2 = -2(1+g^2) -> -2(1+2 g^2) in S1", ["S1"]),
    "largex": ("far-field bracket (x > 2 max(1, g)): 3(g^2-1-2gx) -> 3(g^2-1-gx)",
               ["N1a", "N2", "R2", "R3", "G5b", "G1b", "G2b"]),
    "centre": ("the centred domain's origin x0 = gamma_1 -> 0, the legacy origin", ["G5b", "G1b", "G7b", "D1b"]),
    "scale": ("T of the new path drops the scale factor sc = max(1, gamma_1)", ["U3", "G5b"]),
    "lelimit": ("the -10 of the limit potential -> -9", ["LE0", "N3b", "LE6"]),
    "lefv": ("the finite-volume mass weight z/a-hat -> z", ["LE3", "P-LC"]),
}


def mut(name):
    return MUT == name


def report(cid, ok, msg):
    ok = bool(ok)
    RESULTS.setdefault(cid, []).append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {cid}: {msg}", flush=True)
    if not ok and cid not in FAILED:
        FAILED.append(cid)


def rel(x, y):
    return abs(x - y) / abs(y)


# ---------------------------------------------------------------------------------------------------------------------
# the potential
# ---------------------------------------------------------------------------------------------------------------------
def W_printed(x, g, ns=math):
    """[GGS09a] Eq. (32), written as printed (ns = math for floats, mpmath for 30-digit values)."""
    at = ns.atan(x)
    q = 1 + x * x
    F = 1 + g / (1 + g * g) * (1 + g * x) * at
    Fx = g / (1 + g * g) * ((1 + g * x) / q + g * at)
    r = Fx / F
    k = 3 if mut("printedW") else 4
    br = -3 / q + 3 * (x - g) ** 2 / q ** 2 + 2 * r * r - k * g / q * r + 4 * g / (1 + g * g) * (x - g) / (F * q * q)
    return ns.exp(4 * g * at) * br


def W_A(x, g):
    w = W_printed(x, g)
    return w + 0.1 / (1 + x * x) ** 2 if mut("potential") else w


_DERIVED = {}


def derived_W():
    """W = -V + 2 (dPsi0/Psi0)^2 from the definitions: V of Eq. (24) with the background of Eqs. (11)-(14), Psi0 of Eq. (28)."""
    if "f" not in _DERIVED:
        x, g = sp.symbols("x g", real=True)
        at = sp.atan(x)
        d = g * at
        a = -d
        c = sp.log(1 + x ** 2) / 2 - g * at
        cx = sp.diff(c, x)
        kphi = -2 * (1 + g ** 2)  # kappa Phi_1^2, Eq. (14); Phi_x = Phi_1/(1+x^2) for b = 1
        V = sp.exp(2 * (d - c)) - sp.exp(2 * (d - a)) * cx ** 2 - kphi / (1 + x ** 2) ** 2 * sp.exp(2 * (d - c)) / cx ** 2
        Psi0 = sp.exp(c) / (x - g) * (1 + g / (1 + g ** 2) * (1 + g * x) * at)
        dPsi = sp.exp(2 * d) * sp.diff(Psi0, x) / Psi0
        W = -V + 2 * dPsi ** 2
        _DERIVED["V"] = sp.lambdify((x, g), V, modules="mpmath", cse=True)
        _DERIVED["f"] = sp.lambdify((x, g), W, modules="mpmath", cse=True)
    return _DERIVED["f"], _DERIVED["V"]


def W_B_vec(xs, g):
    """derived W on an array, evaluated with 40 digits (the two terms cancel like 1/(x-g)^2 at the throat)."""
    f, _ = derived_W()
    old = mp.mp.dps
    mp.mp.dps = 40
    try:
        gm = mp.mpf(g)

        def one(v):
            xv = mp.mpf(float(v))
            if abs(xv - gm) < mp.mpf("1e-14"):  # removable 0/0 exactly at the throat: symmetric average of the two sides
                dl = mp.mpf("1e-8")
                return float((f(xv + dl, gm) + f(xv - dl, gm)) / 2)
            return float(f(xv, gm))

        out = np.array([one(v) for v in xs])
    finally:
        mp.mp.dps = old
    if mut("potential"):
        out = out + 0.1 / (1 + xs * xs) ** 2
    return out


def T_of(beta, g):
    f = math.exp(2 * g * math.atan(g)) if not mut("tfactor") else 1.0
    return f / (math.sqrt(1 + g * g) * beta)


# ---------------------------------------------------------------------------------------------------------------------
# the map x(rho), planted potentials
# ---------------------------------------------------------------------------------------------------------------------
def x_of_rho_end(g, L, lam=1.0):
    """x at rho~ = +L and -L from x_rho~ = e^{2 g atan x}/lam, x(0) = 0 (A's own integration); rho~ = lam rho."""
    s = -1.0 if mut("xmap") else 1.0
    out = []
    for sg in (1.0, -1.0):
        sol = solve_ivp(lambda r, y: [math.exp(s * 2 * g * math.atan(y[0])) / lam], (0.0, sg * L), [0.0], method="DOP853",
                        rtol=1e-13, atol=1e-14)
        out.append(float(sol.y[0, -1]))
    return out[0], out[1]


def rho_spline(g, xmax=1e5, n=2001):
    """rho(x) = int_0^x e^{-2 g atan y} dy by composite Gauss-Legendre panels, with Hermite interpolation. Independent of
    the ODE for x(rho) used by A and B."""
    xs_pos = np.sinh(np.linspace(0.0, math.asinh(xmax), n))
    gl_x, gl_w = np.polynomial.legendre.leggauss(24)

    def cum(sgn):
        vals = [0.0]
        for k in range(1, n):
            lo, hi = xs_pos[k - 1], xs_pos[k]
            y = 0.5 * (hi - lo) * gl_x + 0.5 * (hi + lo)
            vals.append(vals[-1] + 0.5 * (hi - lo) * float(np.sum(gl_w * np.exp(-sgn * 2 * g * np.arctan(y)))))
        return np.array(vals)

    rp = cum(+1.0)  # rho(x) for x >= 0
    rn = -cum(-1.0)  # rho(-x) = -int_0^x e^{+2 g atan u} du
    xs = np.concatenate([-xs_pos[:0:-1], xs_pos])
    rs = np.concatenate([rn[:0:-1], rp])
    return CubicHermiteSpline(xs, rs, np.exp(-2 * g * np.arctan(xs)))


def planted_W(ell, g):
    """W_pl(x) = -ell(ell+1) sech^2(rho(x)); exact ground state beta = ell (bound states ell - n > 0)."""
    spl = rho_spline(g)
    lam = ell * (ell + 1)

    def sech2(r):
        e = np.exp(-2.0 * np.abs(r))
        return 4.0 * e / (1.0 + e) ** 2

    def W(x, g_=None):
        return -lam * float(sech2(spl(x)))

    def Wv(xs, g_=None):
        return -lam * sech2(spl(xs))

    return W, Wv


# ---------------------------------------------------------------------------------------------------------------------
# method A: shooting in rho
# ---------------------------------------------------------------------------------------------------------------------
def A_match(beta, g, Wsc, L, lam, xR, xL, rtol):
    """normalized Wronskian of the right- and left-decaying solutions at rho~ = 0. beta is in the rescaled units rho~ = lam rho
    (beta~ = beta/lam, W~ = W/lam^2). The unit vectors (chi, P) are renormalized every 2.5 in rho~ (the system is linear)."""
    s = -1.0 if mut("xmap") else 1.0
    gb = beta + (1.0 / (beta * L * L)) / (1.0 + 1.0 / (beta * L))
    nseg = max(1, int(math.ceil(L / 2.5)))
    vecs = []

    def rhs(r, y):
        return [y[1], (Wsc(y[2], g) / (lam * lam) + beta * beta) * y[0], math.exp(s * 2 * g * math.atan(y[2])) / lam]

    for sg, x0 in ((1.0, xR), (-1.0, xL)):
        y = [1.0, -sg * gb, x0]
        edges = np.linspace(sg * L, 0.0, nseg + 1)
        for lo, hi in zip(edges[:-1], edges[1:]):
            sol = solve_ivp(rhs, (lo, hi), y, method="DOP853", rtol=rtol, atol=[1e-16, 1e-16, 1e-14])
            y = [float(v) for v in sol.y[:, -1]]
            h = math.hypot(y[0], y[1])
            y[0] /= h
            y[1] /= h
        vecs.append((y[0], y[1]))
    (cR, pR), (cL, pL) = vecs
    return cR * pL - cL * pR


def beta_hi_of(g, Wsc):
    grid = np.linspace(-60, 60, 2401)
    return math.sqrt(-min(Wsc(float(v), g) for v in grid)) + 0.01


def method_A(g, Wsc, L=30.0, rtol=1e-13, scan_n=40, beta_hi=None, beta_lo=0.05, lam=1.0):
    """all roots (descending, in the original units of rho) of the matching function; the first is the ground state. The search
    is done in beta~ = beta/lam over [beta_lo, beta_hi/lam]."""
    xR, xL = x_of_rho_end(g, L, lam)
    if beta_hi is None:
        beta_hi = beta_hi_of(g, Wsc)
    bs = np.linspace(beta_lo, beta_hi / lam, scan_n)
    fs = [A_match(b, g, Wsc, L, lam, xR, xL, 1e-9) for b in bs]
    roots = []
    for i in range(scan_n - 1):
        if fs[i] == 0.0 or fs[i] * fs[i + 1] < 0:
            roots.append(lam * brentq(lambda b: A_match(b, g, Wsc, L, lam, xR, xL, rtol), bs[i], bs[i + 1], xtol=1e-14,
                                      rtol=1e-14, maxiter=200))
    return sorted(roots, reverse=True)


# ---------------------------------------------------------------------------------------------------------------------
# method B: finite differences in rho
# ---------------------------------------------------------------------------------------------------------------------
def method_B(g, Wvec, L=30.0, N=3000, lam=1.0):
    """(beta_0, number of negative eigenvalues, eigenvalues) of -d^2/drho~^2 + W/lam^2 on [-L, L] in rho~ = lam rho, Dirichlet,
    sixth order. beta_0 and the eigenvalues are returned in the original units (beta = lam beta~, E = lam^2 E~)."""
    h = 2.0 * L / N
    rho = h * (np.arange(N + 1) - N // 2)  # exact zero at the centre
    npos = rho[N // 2:]
    nneg = rho[:N // 2 + 1][::-1]
    sp_ = solve_ivp(lambda r, y: [math.exp(2 * g * math.atan(y[0])) / lam], (0.0, npos[-1]), [0.0], method="DOP853", rtol=1e-13,
                    atol=1e-14, t_eval=npos)
    sn_ = solve_ivp(lambda r, y: [math.exp(2 * g * math.atan(y[0])) / lam], (0.0, nneg[-1]), [0.0], method="DOP853", rtol=1e-13,
                    atol=1e-14, t_eval=nneg)
    xs = np.concatenate([sn_.y[0][::-1][:-1], sp_.y[0]])
    xi = xs[1:-1]
    M = N - 1
    Wi = Wvec(xi, g) / (lam * lam)
    c0 = -2.7 if mut("matrix") else -49.0 / 18.0
    c1, c2, c3 = 3.0 / 2.0, -3.0 / 20.0, 1.0 / 90.0
    ab = np.zeros((4, M))
    ab[3, :] = -c0 / h ** 2 + Wi
    ab[2, 1:] = -c1 / h ** 2
    ab[1, 2:] = -c2 / h ** 2
    ab[0, 3:] = -c3 / h ** 2
    ev = eig_banded(ab, lower=False, eigvals_only=True, select="v", select_range=(-80.0, 0.0))
    ev = np.sort(ev)
    if len(ev) == 0:
        return float("nan"), 0, ev
    return lam * math.sqrt(-ev[0]), len(ev), ev * lam * lam


# ---------------------------------------------------------------------------------------------------------------------
# method C: gauge dPhi = 0, no W
# ---------------------------------------------------------------------------------------------------------------------
def odeC_rhs(x, C, Cp, beta, g, ns=cmath):
    """C'' = -2x/(1+x^2) C' + beta^2 e^{-4 g atan x} C + 2 (C' - g C/(1+x^2))/(x - g); derived in check L1."""
    k = 1.0 if mut("deltaphi") else 2.0
    return -2 * x / (1 + x * x) * Cp + beta * beta * ns.exp(-4 * g * ns.atan(x)) * C + k * (Cp - g * C / (1 + x * x)) / (x - g)


def C_det(beta, g, L=25.0, rfrac=0.5, side=+1, rtol=1e-12, lam=1.0):
    """(Re, |Im|) of the normalized determinant of (C, C') of the right- and left-decaying solutions at x_t - r, the right one
    continued on the upper (side=+1) or lower (side=-1) semicircle of radius r about the throat x_t = g."""
    xt = g
    r = rfrac * math.sqrt(1 + g * g)
    if mut("cpath"):
        r = 0.0
    xR, xL = x_of_rho_end(g, L, lam)

    def real_rhs(x, y):
        v = [y[1], odeC_rhs(x, y[0], y[1], beta, g, ns=math)]
        if not (math.isfinite(v[0]) and math.isfinite(v[1])):
            raise ValueError("non-finite right-hand side")
        return v

    sR = solve_ivp(real_rhs, (xR, xt + r), [1.0, -beta * math.exp(-2 * g * math.atan(xR))], method="DOP853", rtol=rtol,
                   atol=[1e-16, 1e-16])
    sL = solve_ivp(real_rhs, (xL, xt - r), [1.0, beta * math.exp(-2 * g * math.atan(xL))], method="DOP853", rtol=rtol,
                   atol=[1e-16, 1e-16])

    def arc_rhs(th, y):
        x = xt + r * cmath.exp(1j * th)
        v = [1j * r * cmath.exp(1j * th) * y[1], 1j * r * cmath.exp(1j * th) * odeC_rhs(x, y[0], y[1], beta, g)]
        if not (cmath.isfinite(v[0]) and cmath.isfinite(v[1])):
            raise ValueError("non-finite right-hand side")
        return v

    sA = solve_ivp(arc_rhs, (0.0, side * math.pi), [complex(sR.y[0, -1]), complex(sR.y[1, -1])], method="DOP853", rtol=rtol,
                   atol=[1e-16, 1e-16])
    cR, pR = sA.y[0, -1], sA.y[1, -1]
    cL, pL = sL.y[0, -1], sL.y[1, -1]
    nR, nL = math.hypot(abs(cR), abs(pR)), math.hypot(cL, pL)
    det = (cR / nR) * (pL / nL) - (pR / nR) * (cL / nL)
    return det.real, abs(det.imag)


def method_C(g, L=25.0, rfrac=0.5, side=+1, scan=(0.2, 4.0, 60), rtol=1e-12, centre=None, lam=1.0):
    """roots (descending) of Re det over the scan; with centre, only the bracket centre*(1 +- 0.1)."""
    if centre is not None:
        bs = np.linspace(centre * 0.9, centre * 1.1, 5)
    else:
        bs = np.linspace(scan[0], scan[1], scan[2])
    fs = [C_det(b, g, L, rfrac, side, rtol, lam)[0] for b in bs]
    roots = []
    for i in range(len(bs) - 1):
        if fs[i] * fs[i + 1] < 0:
            roots.append(brentq(lambda b: C_det(b, g, L, rfrac, side, rtol, lam)[0], bs[i], bs[i + 1], xtol=1e-14, rtol=1e-14,
                                maxiter=200))
    return sorted(roots, reverse=True)


# ---------------------------------------------------------------------------------------------------------------------
# symbolic checks
# ---------------------------------------------------------------------------------------------------------------------
def check_units():
    a, g = sp.symbols("a g", positive=True)
    m = g * a
    atn = sp.atan(g)
    # atan(g) + atan(1/g) = pi/2 for g > 0
    ident = sp.simplify(sp.diff(sp.atan(g) + sp.atan(1 / g), g))
    ok0 = ident == 0 and sp.simplify((sp.atan(g) + sp.atan(1 / g)).subs(g, 1) - sp.pi / 2) == 0
    # [SN26] (p. 2): u(X) = (m/a)(atan X - pi/2); throat at X = m/a; R_min = e^{-u(r_t)} sqrt(a^2+m^2); alpha_th = e^{u(r_t)}
    u_t = g * (atn - sp.pi / 2)
    Rmin_sn = sp.exp(-u_t) * sp.sqrt(a ** 2 + m ** 2)
    Rmin_sn_closed = sp.sqrt(a ** 2 + m ** 2) * sp.exp(g * (sp.pi / 2 - atn))  # = sqrt(a^2+m^2) exp[(m/a) atan(a/m)]
    # [GGS09a]: b = a, gamma_0 = -pi g/2 (the mutation flips its sign), B = b e^{-gamma_0}
    g0 = sp.pi * g / 2 if mut("units") else -sp.pi * g / 2
    B = a * sp.exp(-g0)
    r_throat = B * sp.sqrt(1 + g ** 2) * sp.exp(-g * atn)  # p. 4
    lapse_throat = sp.exp(g * atn + g0)  # e^{d(x_throat)} with d = g atan(x/b) + gamma_0, x_throat = g b
    mass_plus = B * g * sp.exp(-g * sp.pi / 2)  # p. 4
    res = {
        "R_min closed form": sp.simplify(Rmin_sn - Rmin_sn_closed),
        "r_throat([GGS09a]) - R_star": sp.simplify(r_throat - Rmin_sn),
        "lapse at the throat - alpha_th": sp.simplify(lapse_throat - sp.exp(u_t)),
        "m_{+inf} - m": sp.simplify(mass_plus - m),
    }
    # the metric: [SN26] -e^{2u} dt^2 + e^{-2u}[drho_E^2 + (rho_E^2+a^2) dOmega^2], X = rho_E/a  vs  [GGS09a] Eq. (13)
    X = sp.symbols("X", real=True)
    gtt_sn = sp.exp(2 * g * (sp.atan(X) - sp.pi / 2))
    gtt_ggs = sp.exp(2 * g * sp.atan(X) + 2 * g0)
    res["g_tt"] = sp.simplify(gtt_sn - gtt_ggs)
    grr_sn = sp.exp(-2 * g * (sp.atan(X) - sp.pi / 2))  # coefficient of drho_E^2 (and of (rho_E^2 + a^2) dOmega^2)
    grr_ggs = sp.exp(-2 * g * sp.atan(X) - 2 * g0)
    res["g_rr"] = sp.simplify(grr_sn - grr_ggs)
    # isotropic radius: rho_E = r - a^2/(4r): d rho_E = Omega dr, rho_E^2 + a^2 = r^2 Omega^2
    r = sp.symbols("r", positive=True)
    rhoE = r - a ** 2 / (4 * r)
    Om = 1 + a ** 2 / (4 * r ** 2)
    res["d rho_E = Omega dr"] = sp.simplify(sp.diff(rhoE, r) - Om)
    res["rho_E^2 + a^2 = r^2 Omega^2"] = sp.simplify(rhoE ** 2 + a ** 2 - r ** 2 * Om ** 2)
    # time: coordinate-time e-folding of [SN26] is a e^{pi g}/beta (derived in README sec. 2): T = tau alpha_th / R_star
    beta = sp.symbols("beta", positive=True)
    tau_sn = a * sp.exp(sp.pi * g) / beta
    Tsn = tau_sn * sp.exp(u_t) / Rmin_sn
    Tggs = sp.exp(2 * g * atn) / (sp.sqrt(1 + g ** 2) * beta)
    res["T([SN26]) - T([GGS09a])"] = sp.simplify(Tsn - Tggs)
    allzero = ok0 and all(v == 0 for v in res.values())
    bad = {k: v for k, v in res.items() if v != 0}
    report("U1", allzero, f"atan g + atan(1/g) = pi/2 {ok0}; identities (R_min, r_throat, lapse, mass, g_tt, g_rr, isotropic map, "
                          f"T) all 0: {allzero}" + (f"; nonzero: {bad}" if bad else ""))
    # the printed numbers at a = 2, m = 1
    Rs = math.sqrt(5) * math.exp(0.5 * math.atan(2))
    al = math.exp(-0.5 * math.atan(2))
    T5 = 5.13 * al / Rs
    ok = abs(Rs - 3.8895) < 6e-5 and abs(al - 0.575) < 6e-4 and abs(T5 - 0.7584) < 1.5e-3
    if mut("units"):
        ok = ok and allzero  # the mutated map cannot be hidden by the printed numbers
    report("U1", ok, f"printed numbers: R_star(2,1) = {Rs:.5f} (3.8895), alpha_th = {al:.5f} (0.575), 5.13 alpha_th/R_star = "
                     f"{T5:.5f} (0.7584)")


def check_T_formula():
    worst = 0.0
    for g in (0.05, 0.3, 0.5, 1.0, 2.0, 3.0):
        beta = 1.7
        a = 1.0
        Rs = a * math.sqrt(1 + g * g) * math.exp(g * math.atan(1 / g))  # [SN26] closed form, m = g a
        al = math.exp(-g * math.atan(1 / g))
        tau = a * math.exp(math.pi * g) / beta
        Tsn = tau * al / Rs
        worst = max(worst, rel(T_of(beta, g), Tsn))
    report("U2", worst < 1e-13, f"code T(beta, g) against [SN26]'s tau alpha_th/R_star, worst rel {worst:.2e}")


# Ricci tensor of a diagonal 4-metric by Christoffel symbols
def ricci(gm, X):
    n = 4
    gi = sp.diag(*[1 / gm[i, i] for i in range(n)])
    Gam = [[[sum(gi[i, l] * (sp.diff(gm[l, j], X[k]) + sp.diff(gm[l, k], X[j]) - sp.diff(gm[j, k], X[l])) for l in range(n)) / 2
             for k in range(n)] for j in range(n)] for i in range(n)]

    def R(j, k):
        return (sum(sp.diff(Gam[i][j][k], X[i]) for i in range(n)) - sum(sp.diff(Gam[i][j][i], X[k]) for i in range(n))
                + sum(Gam[i][i][l] * Gam[l][j][k] for i in range(n) for l in range(n))
                - sum(Gam[i][k][l] * Gam[l][j][i] for i in range(n) for l in range(n)))

    return R


def check_background():
    t, x, th, ph = sp.symbols("t x theta phi")
    g = sp.symbols("g", real=True)
    at = sp.atan(x)
    gm = sp.diag(-sp.exp(2 * g * at), sp.exp(-2 * g * at), sp.exp(-2 * g * at) * (1 + x ** 2),
                 sp.exp(-2 * g * at) * (1 + x ** 2) * sp.sin(th) ** 2)
    R = ricci(gm, (t, x, th, ph))
    kphi = -2 * (1 + 2 * g ** 2) if mut("background") else -2 * (1 + g ** 2)
    ok = True
    worst = {}
    for (i, j) in ((0, 0), (0, 1), (1, 1), (2, 2), (3, 3), (1, 2)):
        src = kphi / (1 + x ** 2) ** 2 if (i, j) == (1, 1) else 0  # R_mn = kappa Phi_m Phi_n, Phi = Phi(x)
        e = sp.simplify(R(i, j) - src)
        worst[(i, j)] = e
        ok = ok and e == 0
    report("S1", ok, "R_mn - kappa dPhi_m dPhi_n = 0 for tt, tx, xx, thth, phph, xth with kappa Phi_1^2 = "
                     f"{'-2(1+2g^2)' if mut('background') else '-2(1+g^2)'}" + ("" if ok else f"; residuals {worst}"))


def check_linearization():
    t, x, th, ph, eps, beta, g, k = sp.symbols("t x theta phi epsilon beta g k")
    d0f, a0f, c0f, Df, Af, Cf = [sp.Function(n) for n in ("d0", "a0", "c0", "D", "A", "C")]
    d0, a0, c0, D, A, C = d0f(x), a0f(x), c0f(x), Df(x), Af(x), Cf(x)
    E = sp.exp(beta * t)
    gm = sp.diag(-sp.exp(2 * (d0 + eps * D * E)), sp.exp(2 * (a0 + eps * A * E)), sp.exp(2 * (c0 + eps * C * E)),
                 sp.exp(2 * (c0 + eps * C * E)) * sp.sin(th) ** 2)
    R = ricci(gm, (t, x, th, ph))
    lin = {}
    for nm, (i, j) in {"tt": (0, 0), "tx": (0, 1), "xx": (1, 1), "thth": (2, 2)}.items():
        lin[nm] = sp.diff(R(i, j), eps).subs(eps, 0)
    # scalar equation box Phi = 0 for Phi = Phi(x) (gauge dPhi = 0): sqrt(-g) = e^{d+a+2c} sin(theta), g^{xx} = e^{-2a}; the
    # t term d_t(sqrt(-g) g^{tt} Phi_t) vanishes because Phi_t = 0
    Phi = sp.Function("Phi")(x)
    sg = sp.exp((d0 + eps * D * E) + (a0 + eps * A * E) + 2 * (c0 + eps * C * E))
    box = sp.diff(sg * sp.exp(-2 * (a0 + eps * A * E)) * sp.diff(Phi, x), x) / sg
    lin["scalar"] = sp.diff(box, eps).subs(eps, 0)
    # explicit background (b = 1, gamma_0 = 0), gauge dPhi = 0
    at = sp.atan(x)
    d0e, a0e, c0e = g * at, -g * at, sp.log(1 + x ** 2) / 2 - g * at
    Cs, Cp = sp.symbols("Cs Cp")
    dd, cc = sp.diff(d0e, x), sp.diff(c0e, x)
    R2 = odeC_rhs(x, Cs, Cp, beta, g, ns=sp)
    Aexpr = (Cp + (cc - dd) * Cs) / cc  # from tx: c' A = C' + (c' - d') C
    Dexpr = Aexpr - 2 * Cs + k  # from the scalar equation: (D - A + 2C)' = 0
    # first part: the tx equation is solved by Aexpr (checked below with the others), and the scalar equation by Dexpr
    def Dx(e):
        return sp.diff(e, x) + Cp * sp.diff(e, Cs) + R2 * sp.diff(e, Cp)

    A1 = Dx(Aexpr)
    A2 = Dx(A1)
    D1 = Dx(Dexpr)
    D2 = Dx(D1)
    mp.mp.dps = 30
    rng = np.random.default_rng(7)
    worst = 0.0
    P1 = sp.symbols("P1")  # Phi = Phi_1 atan x
    for nm, e in lin.items():
        e = e.subs({d0f(x): d0e, a0f(x): a0e, c0f(x): c0e}).doit()
        e = e.subs(Phi, P1 * at).doit()
        e = e.subs({sp.Derivative(Af(x), (x, 2)): A2, sp.Derivative(Df(x), (x, 2)): D2, sp.Derivative(Cf(x), (x, 2)): R2})
        e = e.subs({sp.Derivative(Af(x), x): A1, sp.Derivative(Df(x), x): D1, sp.Derivative(Cf(x), x): Cp})
        e = e.subs({Af(x): Aexpr, Df(x): Dexpr, Cf(x): Cs}).subs(t, 0)
        f = sp.lambdify((x, g, beta, k, Cs, Cp, P1), e, "mpmath")
        w = 0.0
        for _ in range(6):
            args = [mp.mpf(float(v)) for v in (rng.uniform(-3, 3), rng.uniform(0, 2), rng.uniform(0.3, 3), rng.uniform(-1, 1),
                                               rng.uniform(-1, 1), rng.uniform(-1, 1), rng.uniform(0.5, 2))]
            w = max(w, float(abs(f(*args))))
        print(f"      L1 {nm}: max residual {w:.2e}")
        worst = max(worst, w)
    report("L1", worst < 1e-15, f"the five linearized equations (tt, tx, xx, thth, scalar) hold for the reconstructed (dc, da, dd) "
                                f"of the ODE of method C, worst residual {worst:.2e}")


def check_W1():
    mp.mp.dps = 30
    f, Vf = derived_W()
    xs = [-6.0, -3.0, -1.0, 0.3, 1.0, 2.5, 5.0, 9.0]
    for g in (0.0, 0.5, 2.0):
        gm = mp.mpf(g)
        worst = 0.0
        for x in xs:
            xm = mp.mpf(x)
            wp = W_printed(xm, gm, mp)
            wd = f(xm, gm)
            worst = max(worst, float(abs(wp - wd) / abs(wd)))
        report("W1", worst < 1e-12, f"gamma_1 = {g}: printed Eq. (32) against -V + 2(dPsi0/Psi0)^2, worst rel {worst:.2e}")
    # static V of p. 5 against Eq. (24) built from the background
    worst = 0.0
    for g in (0.5, 2.0):
        for x in xs:
            xm, gm = mp.mpf(x), mp.mpf(g)
            Vp = mp.exp(4 * gm * mp.atan(xm)) / (xm * xm + 1) * (1 - (xm - gm) ** 2 / (xm * xm + 1) + 2 * (1 + gm * gm) / (xm - gm) ** 2)
            worst = max(worst, float(abs(Vp - Vf(xm, gm)) / abs(Vp)))
    report("W1", worst < 1e-12, f"static V (p. 5) against Eq. (24) with the background, worst rel {worst:.2e}")
    worst = max(float(abs(W_printed(mp.mpf(x), mp.mpf(0), mp) + 3 / (1 + mp.mpf(x) ** 2) ** 2) * (1 + mp.mpf(x) ** 2) ** 2 / 3) for x in xs)
    report("W1", worst < 1e-14, f"gamma_1 = 0: W = -3/(1+x^2)^2, worst rel {worst:.2e}")
    worst = 0.0
    for g in (0.5, 2.0):
        gm = mp.mpf(g)
        closed = mp.exp(4 * gm * mp.atan(gm)) * (-3 / (1 + gm ** 2) - 2 * gm ** 2 / (1 + gm ** 2) ** 2)
        worst = max(worst, float(abs(W_printed(gm, gm, mp) - closed) / abs(closed)))
    report("W1", worst < 1e-12, f"W(gamma_1) closed form of p. 7, worst rel {worst:.2e}")
    worst = 0.0
    for g in (0.5, 2.0):
        for sg in (1, -1):
            xm = mp.mpf(sg * 10 ** 6)
            worst = max(worst, float(abs(xm * xm * W_printed(xm, mp.mpf(g), mp) / (2 * mp.exp(sg * 2 * mp.pi * g)) - 1)))
    report("W1", worst < 1e-4, f"x^2 W -> 2 e^(+-2 pi g) at |x| = 1e6, worst rel {worst:.2e}")


# ---------------------------------------------------------------------------------------------------------------------
# planted truth
# ---------------------------------------------------------------------------------------------------------------------
def check_planted(gammas, L):
    for g in gammas:
        for ell in (0.9, 1.6):
            exact = ell * (1.01 if mut("planted") else 1.0)
            nb = sum(1 for n in range(0, 5) if ell - n > 0)
            Wsc, Wv = planted_W(ell, g)
            try:
                roots = method_A(g, Wsc, L=L, scan_n=30, beta_hi=math.sqrt(ell * (ell + 1)) + 0.01)
                ok = len(roots) == nb and rel(roots[0], exact) < 1e-9
                report("P-A", ok, f"gamma_1 = {g}, ell = {ell}: A beta = {roots[0]:.12f} (exact {exact}), rel {rel(roots[0], exact):.1e}, "
                                  f"{len(roots)} bound states (exact {nb})")
            except Exception as exc:  # noqa: BLE001
                report("P-A", False, f"gamma_1 = {g}, ell = {ell}: A raised {type(exc).__name__}: {exc}")
            try:
                b0, nneg, _ = method_B(g, Wv, L=L)
                ok = nneg == nb and rel(b0, exact) < 1e-8
                report("P-B", ok, f"gamma_1 = {g}, ell = {ell}: B beta = {b0:.12f}, rel {rel(b0, exact):.1e}, {nneg} negative "
                                  f"eigenvalues (exact {nb})")
            except Exception as exc:  # noqa: BLE001
                report("P-B", False, f"gamma_1 = {g}, ell = {ell}: B raised {type(exc).__name__}: {exc}")


def check_planted_C(gammas):
    x, g = sp.symbols("x g")
    beta = 1.3
    at = sp.atan(x)
    Cstar = sp.exp(g * at) * sp.exp(-(x - g) ** 2 / 2)
    Cs_, Cp_ = sp.symbols("Cs Cp")
    R2 = odeC_rhs(x, Cs_, Cp_, beta, g, ns=sp)
    S = sp.diff(Cstar, x, 2) - R2.subs({Cs_: Cstar, Cp_: sp.diff(Cstar, x)})
    fC = sp.lambdify((x, g), Cstar, "numpy")
    fCp = sp.lambdify((x, g), sp.diff(Cstar, x), "numpy")
    fS = sp.lambdify((x, g), S, "numpy")
    for gv in gammas:
        r = 0.0 if mut("cpath") else 0.5 * math.sqrt(1 + gv * gv)
        xt = gv
        worst = 0.0
        try:
            with np.errstate(all="ignore"):
                for side in (+1, -1):
                    def rhs(th, y):
                        xx = xt + r * cmath.exp(1j * th)
                        dx = 1j * r * cmath.exp(1j * th)
                        v = [dx * y[1], dx * (odeC_rhs(xx, y[0], y[1], beta, gv) + complex(fS(xx, gv)))]
                        if not all(cmath.isfinite(z) for z in v):
                            raise ValueError("non-finite right-hand side")
                        return v

                    x0 = xt + r
                    sol = solve_ivp(rhs, (0.0, side * math.pi), [complex(fC(x0, gv)), complex(fCp(x0, gv))], method="DOP853",
                                    rtol=1e-12, atol=[1e-16, 1e-16])
                    x1 = xt - r
                    ex = np.array([complex(fC(x1, gv)), complex(fCp(x1, gv))])
                    got = np.array([sol.y[0, -1], sol.y[1, -1]])
                    e = float(np.linalg.norm(got - ex) / np.linalg.norm(ex))
                    worst = max(worst, e if np.isfinite(e) else float("inf"))
        except Exception as exc:  # noqa: BLE001
            print(f"      P-C integration raised {type(exc).__name__}: {exc}")
            worst = float("inf")
        report("P-C", worst < 1e-8, f"gamma_1 = {gv}: manufactured solution carried across the throat on both semicircles, "
                                    f"worst rel error {worst:.2e}")


# ---------------------------------------------------------------------------------------------------------------------
# the real runs
# ---------------------------------------------------------------------------------------------------------------------
TABLE = {0.0: 0.846, 0.1: 0.841, 0.2: 0.826, 0.3: 0.805, 0.4: 0.782, 0.5: 0.758, 0.6: 0.737, 0.7: 0.718, 0.8: 0.701,
         0.9: 0.687, 1.0: 0.675, 1.2: 0.656, 1.4: 0.642, 1.6: 0.631, 1.8: 0.624, 2.0: 0.618}


def lam_of(g):
    """rho~ = lam rho, lam = e^{2 g atan g}: [GGS09a]'s beta-bar scaling (p. 9), which keeps beta~ = beta/lam and W~ = W/lam^2 of
    order one at the throat for large g. Pure conditioning: T = 1/(sqrt(1+g^2) beta~)."""
    return math.exp(2 * g * math.atan(g))


def run_G(gammas, L, full):
    out = {}
    # ------- A and B at the control points
    ctrl = [g for g in (0.0, 0.5, 1.0) if g in gammas]
    for g in ctrl:
        lam = lam_of(g)
        out[g] = {}
        try:
            out[g]["A"] = method_A(g, W_A, L=L, lam=lam)[0]
        except Exception as exc:  # noqa: BLE001
            report("G1", False, f"gamma_1 = {g}: A raised {type(exc).__name__}: {exc}")
        try:
            out[g]["B"] = method_B(g, W_B_vec, L=L, lam=lam)[0]
        except Exception as exc:  # noqa: BLE001
            report("G1", False, f"gamma_1 = {g}: B raised {type(exc).__name__}: {exc}")
        if "A" in out[g] and "B" in out[g]:
            bA, bB = out[g]["A"], out[g]["B"]
            print(f"      gamma_1 = {g}: beta_A = {bA:.12f}, beta_B = {bB:.12f}, T_A = {T_of(bA, g):.9f}, T_B = {T_of(bB, g):.9f}",
                  flush=True)
            report("G1", rel(bA, bB) < 1e-7, f"gamma_1 = {g}: A against B, rel {rel(bA, bB):.2e}")
    # ------- G0
    if 0.0 in gammas:
        if "A" in out[0.0] and "B" in out[0.0]:
            bA, bB = out[0.0]["A"], out[0.0]["B"]
            for nm, b in (("A", bA), ("B", bB)):
                T = T_of(b, 0.0)
                ok = 0.8455 <= T <= 0.8465 and 1.3965 <= b * b <= 1.3975 and rel(b * b, 1.39705243) < 1e-6
                report("G0", ok, f"gamma_1 = 0, method {nm}: T = {T:.6f} (0.846), beta^2 = {b * b:.8f} (1.397; repo 1.39705243, rel "
                                 f"{rel(b * b, 1.39705243):.1e})")
            Wrepo = lambda x, gg: -3.0 / (1.0 + x * x) ** 2 + (0.1 / (1 + x * x) ** 2 if mut("potential") else 0.0)  # one_mode.py S3, n = 2
            try:
                bR = method_A(0.0, Wrepo, L=L)[0]
                report("G0", rel(bA * bA, bR * bR) < 1e-9, f"gamma_1 = 0: A with the repo's D = 4 potential -3/r^4 (r^2 = 1 + l^2): "
                                                           f"beta^2 = {bR * bR:.10f}, rel to A with W {rel(bA * bA, bR * bR):.1e}")
            except Exception as exc:  # noqa: BLE001
                report("G0", False, f"gamma_1 = 0: A with the repo's potential raised {type(exc).__name__}: {exc}")
        else:
            report("G0", False, "gamma_1 = 0: A or B gave no value")
    # ------- method C
    for g in ctrl:
        lam = lam_of(g)
        try:
            bhi = beta_hi_of(g, W_A)
            rp = method_C(g, side=+1, scan=(0.2, bhi, 60), lam=lam)
            rm = method_C(g, side=-1, scan=(0.2, bhi, 60), lam=lam)
            bC = rp[0]
            bC2 = rm[0]
            dets = [C_det(b, g, 25.0, 0.5, sd, 1e-12, lam)[1] for b in (0.8 * bC, 1.1 * bC) for sd in (+1, -1)]
            bA = out[g]["A"]
            ok = rel(bC, bA) < 1e-6 and rel(bC, bC2) < 1e-9
            report("G2", ok, f"gamma_1 = {g}: C beta = {bC:.12f} (upper) {bC2:.12f} (lower), roots in [0.2, {bhi:.2f}]: "
                             f"{[round(v, 6) for v in rp]}; rel to A {rel(bC, bA):.2e}; max |Im det| {max(dets):.1e}")
            out[g]["C"] = bC
        except Exception as exc:  # noqa: BLE001
            report("G2", False, f"gamma_1 = {g}: method C or its comparison failed: {type(exc).__name__}: {exc}")
    if 0.5 in gammas:
        g = 0.5
        if "A" in out[g]:
            bA = out[g]["A"]
            T = T_of(bA, g)
            a, m = 2.0, 1.0
            Rs = math.sqrt(a * a + m * m) * math.exp((m / a) * math.atan(a / m))
            al = math.exp(-(m / a) * math.atan(a / m))
            tau_sn = a * math.exp(math.pi * g) / bA
            print(f"      gamma_1 = 1/2 (m/a): beta = {bA:.10f}; T = {T:.8f}; tau_[SN26] (evolution time, units of m) = {tau_sn:.5f} "
                  f"(their 5.13); alpha_th = {al:.6f}; R_star = {Rs:.5f}; tau alpha_th/R_star = {tau_sn * al / Rs:.8f}")
            report("G3", abs(T / 0.7584 - 1) <= 0.002, f"T(1/2) = {T:.6f} against [SN26] linear 0.7584: {100 * (T / 0.7584 - 1):+.3f} %")
            report("G4", abs(T - 0.758) <= 0.0006, f"T(1/2) = {T:.6f} against [GGS09a] Table I 0.758: {T - 0.758:+.6f}")
            T3, T4 = 5.12 * al / Rs, 5.23 * al / Rs
            d3, d4 = 100 * (T3 / T - 1), 100 * (T4 / T - 1)
            report("G6", -0.25 <= d3 <= -0.15 and 1.85 <= d4 <= 1.95, f"[SN26] measured tau 5.12 -> T = {T3:.5f} ({d3:+.3f} %, they print "
                                                                      f"0.2 %), 5.23 -> T = {T4:.5f} ({d4:+.3f} %, they print 1.9 %)")
        else:
            for cid in ("G3", "G4", "G6"):
                report(cid, False, "gamma_1 = 1/2: A gave no value")
    # ------- G5: the table, A only
    if full:
        worst = 0.0
        row = []
        try:
            for g, tv in TABLE.items():
                b = out[g]["A"] if g in out and "A" in out[g] else method_A(g, W_A, L=L, lam=lam_of(g))[0]
                LEG[g] = b
                T = T_of(b, g)
                row.append((g, T, tv))
                worst = max(worst, abs(T - tv))
            for g, T, tv in row:
                print(f"      Table I: gamma_1 = {g:3.1f}  T = {T:.6f}  printed {tv:.3f}  diff {T - tv:+.6f}")
            report("G5", worst <= 0.001, f"16 values of Table I, worst |T - printed| = {worst:.6f}")
            e1 = max(rel(LEG[g], FIRST_ROW_BETA[g]) for g in TABLE)
            report("R1", e1 <= 1e-9, f"first row of this script's legacy path against #151's beta_A (main 2455a0d) at 16 values: worst rel {e1:.2e}")
        except Exception as exc:  # noqa: BLE001
            report("G5", False, f"table incomplete after {len(row)} values: {type(exc).__name__}: {exc}")
    # ------- G7: counts
    for g in [v for v in (0.0, 0.5, 1.0, 2.0) if v in gammas]:
        lam = lam_of(g)
        Wc = (lambda x, gg: W_A(x, gg) - 60.0 / math.cosh(x - 8.0) ** 2) if mut("count") else W_A
        Wcv = (lambda xs, gg: W_B_vec(xs, gg) - 60.0 / np.cosh(xs - 8.0) ** 2) if mut("count") else W_B_vec
        try:
            rootsA = method_A(g, Wc, L=L, lam=lam)
            bB, nB, evB = method_B(g, Wcv, L=L, lam=lam)
            ok = len(rootsA) == 1 and nB == 1
            report("G7", ok, f"gamma_1 = {g}: A finds {len(rootsA)} root(s) in the scan, B has {nB} negative eigenvalue(s) "
                             f"(lowest {evB[0] if len(evB) else float('nan'):.6f}; the next is >= 0)")
        except Exception as exc:  # noqa: BLE001
            report("G7", False, f"gamma_1 = {g}: raised {type(exc).__name__}: {exc}")
    return out


def run_D1(gammas, L):
    for g in [v for v in (0.0, 0.5) if v in gammas]:
        lam = lam_of(g)
        try:
            base = method_A(g, W_A, L=L, lam=lam)[0]
            vA = [method_A(g, W_A, L=L, lam=lam, rtol=1e-11)[0], method_A(g, W_A, L=L + 10, lam=lam)[0]]
            sA = max(rel(v, base) for v in vA)
            bB = method_B(g, W_B_vec, L=L, lam=lam)[0]
            vB = [method_B(g, W_B_vec, L=L, N=6000, lam=lam)[0], method_B(g, W_B_vec, L=L + 10, lam=lam)[0]]
            sB = max(rel(v, bB) for v in vB)
            ok = sA < 1e-7 and sB < 1e-7
            sC = float("nan")
            if not mut("domain") and not mut("cpath"):
                try:
                    bhi = beta_hi_of(g, W_A)
                    bC = method_C(g, scan=(0.2, bhi, 60), lam=lam)[0]
                    vC = [method_C(g, rfrac=0.3, centre=bC, lam=lam)[0], method_C(g, L=30.0, rfrac=0.7, centre=bC, lam=lam)[0]]
                    sC = max(rel(v, bC) for v in vC)
                    ok = ok and sC < 1e-7
                except Exception as exc:  # noqa: BLE001
                    ok = False
                    print(f"      method C variants failed: {exc}")
            report("D1", ok, f"gamma_1 = {g}: max relative spread of beta: A (rtol 1e-13 -> 1e-11, L +10) {sA:.1e}; B (N 3000 -> 6000, "
                             f"L +10) {sB:.1e}; C (r/dist 0.5 -> 0.3, 0.7; L 25 -> 30) {sC:.1e}")
            print(f"      gamma_1 = {g}: beta_A = {base:.12f}, beta_B = {bB:.12f}; error bar of T from the spreads: "
                  f"{max(sA, sB, 0.0 if math.isnan(sC) else sC) * T_of(base, g):.1e}")
        except Exception as exc:  # noqa: BLE001
            report("D1", False, f"gamma_1 = {g}: raised {type(exc).__name__}: {exc}")


# ---------------------------------------------------------------------------------------------------------------------
# second row (#163): throat-centred domain, far-field bracket, the large-gamma_1 limit equation
# ---------------------------------------------------------------------------------------------------------------------
# [GGS09a] Table I, second row (phantom_lin.tex lines 970-981), and the two numbers of Sec. IV C
TABLE2 = {2.5: 0.608, 3.0: 0.602, 3.5: 0.599, 4.0: 0.597, 4.5: 0.595, 5.0: 0.594, 5.5: 0.593, 6.0: 0.592, 6.5: 0.592, 7.0: 0.591,
          7.5: 0.591, 8.0: 0.590, 8.5: 0.590, 9.0: 0.590, 9.5: 0.590, 10.0: 0.590}
LE_PAPER = 0.588  # 1/beta-hat from the limit equation, [GGS09a] Sec. IV C ("yields 1/beta-hat = 0.588")
LE_TABLE = 0.590  # the large-gamma_1 value of the table

# beta_A of #151's unmodified notes/ellis_efold.py (main 2455a0d), L = 30, and the same call with L = 60 (negative-mass-raw
# I163-rpi5: baseline_first_row_main2455a0d.log, legacy_Lspread_L60.log; generated by make_refdicts.py). L = 60 differs from
# L = 30 by 2e-13 at every gamma_1 <= 1.8 and by 6.3e-6 at gamma_1 = 2.0 (the L = 30 domain truncates the mode there).
FIRST_ROW_BETA = {0.0: 1.1819697264877345, 0.1: 1.2075044125173486, 0.2: 1.2849623990241523, 0.3: 1.417330106314019,
                  0.4: 1.6106745259201527, 0.5: 1.8751275118689374, 0.6: 2.225817717182347, 0.7: 2.6839541377683713,
                  0.8: 3.2782931151592187, 0.9: 4.0471770915464695, 1.0: 5.041327581055734, 1.2: 7.994103290218043,
                  1.4: 12.968438949655779, 1.6: 21.409547538236673, 1.8: 35.83719399710155, 2.0: 60.662389104966536}
FIRST_ROW_BETA_L60 = {0.0: 1.1819697264877347, 0.1: 1.2075044125173569, 0.2: 1.2849623990241468, 0.3: 1.4173301063140165,
                      0.4: 1.6106745259201645, 0.5: 1.8751275118689392, 0.6: 2.2258177171823645, 0.7: 2.683954137768372,
                      0.8: 3.2782931151592147, 0.9: 4.047177091546469, 1.0: 5.041327581055731, 1.2: 7.99410329021804,
                      1.4: 12.968438949655776, 1.6: 21.40954753823666, 1.8: 35.83719399709431, 2.0: 60.66200729583162}
LEG = {}  # legacy beta_A at the first-row values, filled by run_G


def sc_of(g):
    """scale of the centred variable: u = (rho~ - rho~_throat)/sc, sc = max(1, gamma_1); beta' = sc beta~ = O(1)."""
    return max(1.0, g)


def datan(x, g):
    """atan x - atan g without cancellation near the throat."""
    return math.atan2(x - g, 1.0 + g * x)


def a_of(x, g):
    """a(x) = e^{2 g [atan x - atan g]} ([GGS09a] Sec. IV): a(gamma_1) = 1, d(rho~)/dx = 1/a."""
    return math.exp(2.0 * g * datan(x, g))


def _br(x, g, ns, far):
    """bracket of [GGS09a] Eq. (32) (W-bar = a^2 * bracket). far=False: as printed. far=True: the same with the first two terms
    combined, -3/q + 3(x-g)^2/q^2 = 3(g^2 - 1 - 2gx)/q^2 (exact), the far-field form."""
    at = ns.atan(x)
    q = 1 + x * x
    F = 1 + g / (1 + g * g) * (1 + g * x) * at
    Fx = g / (1 + g * g) * ((1 + g * x) / q + g * at)
    r = Fx / F
    k = 3 if mut("printedW") else 4
    if far:
        c = -g * x if mut("largex") else -2 * g * x
        first = 3 * (g * g - 1 + c) / q ** 2
    else:
        first = -3 / q + 3 * (x - g) ** 2 / q ** 2
    return first + 2 * r * r - k * g / q * r + 4 * g / (1 + g * g) * (x - g) / (F * q * q)


def Wbar(x, g, far=True):
    """W-bar = W/lam^2 = a(x)^2 * bracket; the far-field form beyond x_sw = 2 max(1, gamma_1) (exact algebra)."""
    return a_of(x, g) ** 2 * _br(x, g, math, far and x > 2.0 * max(1.0, g))


def wbar_derived_vec(xs, g):
    """W-bar from the derived W = -V + 2 (dPsi0/Psi0)^2 (40 digits), the independent reference for the printed form."""
    f, _ = derived_W()
    old = mp.mp.dps
    mp.mp.dps = 40
    try:
        gm = mp.mpf(g)
        pref = mp.exp(-4 * gm * mp.atan(gm))

        def one(v):
            xv = mp.mpf(float(v))
            if abs(xv - gm) < mp.mpf("1e-14"):  # removable 0/0 exactly at the throat: symmetric average
                dl = mp.mpf("1e-8")
                return float((f(xv + dl, gm) + f(xv - dl, gm)) / 2 * pref)
            return float(f(xv, gm) * pref)

        return np.array([one(v) for v in xs])
    finally:
        mp.mp.dps = old


def LN():
    return 2.0 if mut("domain") else 20.0


# a system: dx/du = dxdu(x), potential W(x, u) (scalar) and Wvec(xs, us) in the variable u; x0 = x at u = 0
def sysd_ggs(g, far=True, count=False):
    sc = sc_of(g)
    x0 = 0.0 if mut("centre") else g

    def W(x, u):
        w = sc * sc * Wbar(x, g, far)
        return w - 60.0 / math.cosh(u - 8.0) ** 2 if count else w

    def Wvec(xs, us):
        w = sc * sc * wbar_derived_vec(xs, g)
        return w - 60.0 / np.cosh(us - 8.0) ** 2 if count else w

    return dict(x0=x0, dxdu=lambda x: sc * a_of(x, g), W=W, Wvec=Wvec, sc=sc, g=g)


def le_W(z):
    c = 9.0 if mut("lelimit") else 10.0
    return math.exp(4.0 * (1.0 - 1.0 / z)) * (2.0 / z ** 2 - c / z ** 3 + 3.0 / z ** 4)


def le_Wv(z):
    c = 9.0 if mut("lelimit") else 10.0
    return np.exp(4.0 * (1.0 - 1.0 / z)) * (2.0 / z ** 2 - c / z ** 3 + 3.0 / z ** 4)


def sysd_le(count=False):
    """[GGS09a] Eq. (43): -(a-hat d/dz)^2 chi + W-hat chi = -beta-hat^2 chi, a-hat = e^{2(1 - 1/z)}; u = rho-hat = int_1^z dz/a-hat."""
    def W(z, u):
        return le_W(z) - 60.0 / math.cosh(u - 8.0) ** 2 if count else le_W(z)

    def Wvec(zs, us):
        return le_Wv(zs) - 60.0 / np.cosh(us - 8.0) ** 2 if count else le_Wv(zs)

    return dict(x0=1.0, dxdu=lambda z: math.exp(2.0 * (1.0 - 1.0 / z)), W=W, Wvec=Wvec, sc=1.0, g=None)


def x_ends(sd, L):
    out = []
    for sg in (1.0, -1.0):
        sol = solve_ivp(lambda u, y: [sd["dxdu"](y[0])], (0.0, sg * L), [sd["x0"]], method="DOP853", rtol=1e-13, atol=1e-14)
        out.append(float(sol.y[0, -1]))
    return out[0], out[1]


def x_on(sd, us):
    """x(u) on the array us (any order, either sign) from dx/du = dxdu(x), x(0) = x0."""
    us = np.asarray(us, dtype=float)
    xs = np.empty(us.shape)
    pos = us >= 0
    tp, tn = us[pos], us[~pos]
    if tp.size and tp.max() > 0:
        s = solve_ivp(lambda u, y: [sd["dxdu"](y[0])], (0.0, float(tp.max())), [sd["x0"]], method="DOP853", rtol=1e-13, atol=1e-14,
                      t_eval=np.sort(tp))
        xs[pos] = np.interp(tp, np.sort(tp), s.y[0])
    elif tp.size:
        xs[pos] = sd["x0"]
    if tn.size:
        order = np.argsort(-tn)
        s = solve_ivp(lambda u, y: [sd["dxdu"](y[0])], (0.0, float(tn.min())), [sd["x0"]], method="DOP853", rtol=1e-13, atol=1e-14,
                      t_eval=tn[order])
        tmp = np.empty(tn.shape)
        tmp[order] = s.y[0]
        xs[~pos] = tmp
    return xs


def gmatch(sd, b, L, xR, xL, rtol):
    """normalized Wronskian at u = 0 of the solutions decaying to the right (start chi' = -b chi at u = L) and left (+b chi at -L)."""
    nseg = max(1, int(math.ceil(L / 2.5)))

    def rhs(u, y):
        return [y[1], (sd["W"](y[2], u) + b * b) * y[0], sd["dxdu"](y[2])]

    vecs = []
    for sg, x0 in ((1.0, xR), (-1.0, xL)):
        y = [1.0, -sg * b, x0]
        edges = np.linspace(sg * L, 0.0, nseg + 1)
        for lo, hi in zip(edges[:-1], edges[1:]):
            sol = solve_ivp(rhs, (lo, hi), y, method="DOP853", rtol=rtol, atol=[1e-16, 1e-16, 1e-14])
            y = [float(v) for v in sol.y[:, -1]]
            h = math.hypot(y[0], y[1])
            y[0] /= h
            y[1] /= h
        vecs.append((y[0], y[1]))
    (cR, pR), (cL, pL) = vecs
    return cR * pL - cL * pR


def gshoot(sd, L=20.0, rtol=1e-13, scan_n=40, b_lo=0.05, b_hi=None):
    """all roots b (descending) of the matching function on [b_lo, b_hi]; the first is the ground state (b = sc beta~)."""
    xR, xL = x_ends(sd, L)
    if b_hi is None:
        us = np.linspace(-L, L, 2401)
        xs = x_on(sd, us)
        b_hi = math.sqrt(max(0.0, -min(sd["W"](float(x), float(u)) for x, u in zip(xs, us)))) + 0.01
    bs = np.linspace(b_lo, b_hi, scan_n)
    fs = [gmatch(sd, b, L, xR, xL, 1e-9) for b in bs]
    roots = []
    for i in range(scan_n - 1):
        if fs[i] == 0.0 or fs[i] * fs[i + 1] < 0:
            roots.append(brentq(lambda b: gmatch(sd, b, L, xR, xL, rtol), bs[i], bs[i + 1], xtol=1e-14, rtol=1e-14, maxiter=200))
    return sorted(roots, reverse=True)


def gmatrix(sd, L=20.0, N=4000):
    """(b_0, number of negative eigenvalues, eigenvalues) of -d^2/du^2 + W on [-L, L], Dirichlet, sixth order, h = 2L/N."""
    h = 2.0 * L / N
    us = h * (np.arange(N + 1) - N // 2)
    xs = x_on(sd, us)
    Wi = sd["Wvec"](xs[1:-1], us[1:-1])
    M = N - 1
    c0 = -2.7 if mut("matrix") else -49.0 / 18.0
    c1, c2, c3 = 3.0 / 2.0, -3.0 / 20.0, 1.0 / 90.0
    ab = np.zeros((4, M))
    ab[3, :] = -c0 / h ** 2 + Wi
    ab[2, 1:] = -c1 / h ** 2
    ab[1, 2:] = -c2 / h ** 2
    ab[0, 3:] = -c3 / h ** 2
    ev = np.sort(eig_banded(ab, lower=False, eigvals_only=True, select="v", select_range=(-60.0, 0.0)))
    if len(ev) == 0:
        return float("nan"), 0, ev
    return math.sqrt(-ev[0]), len(ev), ev


def le_fv(Wfun, N, zmin=0.15, zmax=150.0):
    """finite volume on t = ln z, uniform: -(p chi_t)_t + (z W/a) chi = -beta^2 (z/a) chi, p = a/z, a = e^{2(1-1/z)}; shares only W
    with the u-variable solvers. Returns (lowest eigenvalue lambda = -beta^2, number of eigenvalues in [-60, 0])."""
    t = np.linspace(math.log(zmin), math.log(zmax), N + 2)
    h = t[1] - t[0]
    z = np.exp(t)
    zi = z[1:-1]
    ai = np.exp(2.0 * (1.0 - 1.0 / zi))
    zm = np.exp(0.5 * (t[:-1] + t[1:]))
    pm = np.exp(2.0 * (1.0 - 1.0 / zm)) / zm
    mass = zi if mut("lefv") else zi / ai
    diag = (pm[:-1] + pm[1:]) / h ** 2 + zi * Wfun(zi) / ai
    off = -pm[1:-1] / h ** 2
    d2 = diag / mass
    o2 = off / np.sqrt(mass[:-1] * mass[1:])
    ev = eigh_tridiagonal(d2, o2, eigvals_only=True, select="v", select_range=(-60.0, 0.0))
    ev = np.sort(ev)
    return (float(ev[0]) if len(ev) else float("nan")), len(ev)


def le_fv_beta(Wfun, N=8000):
    """Richardson in h^2 on the eigenvalue (N, 2N)."""
    l1, n1 = le_fv(Wfun, N)
    l2, n2 = le_fv(Wfun, 2 * N)
    lam = (4.0 * l2 - l1) / 3.0
    return math.sqrt(-lam), n2, math.sqrt(-l1), math.sqrt(-l2)


def u_spline(dudx, x0, nodes_lo, nodes_hi):
    """u(x) = int_{x0}^x dudx by 24-point Gauss-Legendre panels with Hermite interpolation; nodes_lo ascending ending at x0 and
    nodes_hi ascending starting at x0. Independent of the ODE for x(u) used by the solvers."""
    gl_x, gl_w = np.polynomial.legendre.leggauss(24)

    def panels(nodes):
        out = []
        for k in range(1, len(nodes)):
            lo, hi = nodes[k - 1], nodes[k]
            y = 0.5 * (hi - lo) * gl_x + 0.5 * (hi + lo)
            out.append(0.5 * (hi - lo) * float(np.sum(gl_w * dudx(y))))
        return np.array(out)

    up = np.concatenate([[0.0], np.cumsum(panels(nodes_hi))])
    pl = panels(nodes_lo)
    ulo = -np.concatenate([np.cumsum(pl[::-1])[::-1], [0.0]])
    xs = np.concatenate([nodes_lo[:-1], nodes_hi])
    us = np.concatenate([ulo[:-1], up])
    return CubicHermiteSpline(xs, us, dudx(xs))


def u_spline_ggs(g, x0, n=4001):
    sc = sc_of(g)
    s0 = max(1.0, g)
    xlo, xhi = -20.0, 1e5
    lo = x0 - s0 * np.sinh(np.linspace(math.asinh((x0 - xlo) / s0), 0.0, n))
    hi = x0 + s0 * np.sinh(np.linspace(0.0, math.asinh((xhi - x0) / s0), n))
    return u_spline(lambda x: np.exp(-2.0 * g * np.arctan2(x - g, 1.0 + g * x)) / sc, x0, lo, hi)


def u_spline_le(n=4001):
    lo = np.exp(np.linspace(math.log(0.05), 0.0, n))
    hi = np.exp(np.linspace(0.0, math.log(1e4), n))
    return u_spline(lambda z: np.exp(-2.0 * (1.0 - 1.0 / z)), 1.0, lo, hi)


def sysd_planted(base, spl, ell):
    """W_pl(x) = -ell(ell+1) sech^2(u(x)), u(x) from the quadrature spline; exact ground state b = ell."""
    lam = ell * (ell + 1)

    def sech2(u):
        e = np.exp(-2.0 * np.abs(u))
        return 4.0 * e / (1.0 + e) ** 2

    sd = dict(base)
    sd["W"] = lambda x, u: -lam * float(sech2(spl(x)))
    sd["Wvec"] = lambda xs, us: -lam * sech2(spl(xs))
    return sd


def T_of2(b, g):
    """T = e^{2 g atan g}/(sqrt(1+g^2) beta) with beta = lam b/sc, i.e. T = sc/(sqrt(1+g^2) b)."""
    sc = 1.0 if mut("scale") else sc_of(g)
    return sc / (math.sqrt(1 + g * g) * b)


def beta_orig(b, g):
    return lam_of(g) * b / sc_of(g)


def odeC2_rhs(x, C, Cp, bt, g, ns=cmath):
    """method C's ODE with beta^2 e^{-4 g atan x} = bt^2/a(x)^2, bt = beta/lam (no overflow at large gamma_1)."""
    k = 1.0 if mut("deltaphi") else 2.0
    a2 = a_of(x, g) ** 2 if ns is math else ns.exp(4 * g * (ns.atan(x) - math.atan(g) if ns is cmath else ns.atan(x) - ns.atan(g)))
    return -2 * x / (1 + x * x) * Cp + bt * bt / a2 * C + k * (Cp - g * C / (1 + x * x)) / (x - g)


def C2_det(b, g, L=20.0, rfrac=0.5, side=+1, rtol=1e-12):
    sc = sc_of(g)
    bt = b / sc
    sd = sysd_ggs(g)
    xR, xL = x_ends(sd, L)
    xt = g
    r = rfrac * math.sqrt(1 + g * g)
    r = min(r, rfrac * (xt - xL)) if xt - xL > 0 else r
    if mut("cpath"):
        r = 0.0

    def real_rhs(x, y):
        v = [y[1], odeC2_rhs(x, y[0], y[1], bt, g, ns=math)]
        if not (math.isfinite(v[0]) and math.isfinite(v[1])):
            raise ValueError("non-finite right-hand side")
        return v

    sR = solve_ivp(real_rhs, (xR, xt + r), [1.0, -bt / a_of(xR, g)], method="DOP853", rtol=rtol, atol=[1e-16, 1e-16])
    sL = solve_ivp(real_rhs, (xL, xt - r), [1.0, bt / a_of(xL, g)], method="DOP853", rtol=rtol, atol=[1e-16, 1e-16])

    def arc_rhs(th, y):
        x = xt + r * cmath.exp(1j * th)
        d = 1j * r * cmath.exp(1j * th)
        v = [d * y[1], d * odeC2_rhs(x, y[0], y[1], bt, g)]
        if not (cmath.isfinite(v[0]) and cmath.isfinite(v[1])):
            raise ValueError("non-finite right-hand side")
        return v

    sA = solve_ivp(arc_rhs, (0.0, side * math.pi), [complex(sR.y[0, -1]), complex(sR.y[1, -1])], method="DOP853", rtol=rtol,
                   atol=[1e-16, 1e-16])
    cR, pR = sA.y[0, -1], sA.y[1, -1]
    cL, pL = sL.y[0, -1], sL.y[1, -1]
    nR, nL = math.hypot(abs(cR), abs(pR)), math.hypot(cL, pL)
    det = (cR / nR) * (pL / nL) - (pR / nR) * (cL / nL)
    return det.real, abs(det.imag)


def method_C2(g, L=20.0, rfrac=0.5, side=+1, scan=(0.2, 3.0, 60), rtol=1e-12, centre=None):
    bs = np.linspace(centre * 0.9, centre * 1.1, 5) if centre is not None else np.linspace(scan[0], scan[1], scan[2])
    fs = [C2_det(b, g, L, rfrac, side, rtol)[0] for b in bs]
    roots = []
    for i in range(len(bs) - 1):
        if fs[i] * fs[i + 1] < 0:
            roots.append(brentq(lambda b: C2_det(b, g, L, rfrac, side, rtol)[0], bs[i], bs[i + 1], xtol=1e-14, rtol=1e-14, maxiter=200))
    return sorted(roots, reverse=True)


# ------------------------------------------------- checks
def Wscale(x, g):
    """size of the terms of W-bar (sum of the absolute values), the scale the cancellations work at."""
    at = math.atan(x)
    q = 1 + x * x
    F = 1 + g / (1 + g * g) * (1 + g * x) * at
    r = g / (1 + g * g) * ((1 + g * x) / q + g * at) / F
    return a_of(x, g) ** 2 * (3 / q + 3 * (x - g) ** 2 / q ** 2 + 2 * r * r + 4 * g / q * abs(r)
                              + 4 * g / (1 + g * g) * abs(x - g) / (F * q * q))


def check_far_field():
    """N1a, N2, and the diagnostic on the printed bracket at large x. Points where |W-bar| < 1e-3 of the size of its terms (next
    to a zero crossing, where a relative error is not meaningful) are skipped and counted."""
    f, _ = derived_W()
    mp.mp.dps = 60
    worst, wp, skipped = 0.0, 0.0, 0
    for g in (2.5, 10.0):
        for k in (1, 2, 3, 4, 6, 8, 10, 12):
            x = 10.0 ** k
            ref = f(mp.mpf(x), mp.mpf(g)) * mp.exp(-4 * mp.mpf(g) * mp.atan(mp.mpf(g)))
            if abs(float(ref)) < 1e-3 * Wscale(x, g):
                skipped += 1
                continue
            for far in (True, False):
                e = float(abs(mp.mpf(Wbar(x, g, far)) / ref - 1))
                if far:
                    worst = max(worst, e)
                else:
                    wp = max(wp, e)
    report("N1a", worst <= 1e-12, f"W-bar, far-field form (double) against the derived W (60 digits), gamma_1 in {{2.5, 10}}, x = 10 ... 1e12: "
                                  f"worst rel {worst:.2e} ({skipped} points skipped)")
    print(f"      DIAG N1: the printed bracket in double at the same points: worst rel {wp:.2e} (the premise that it loses digits at "
          f"large x is {'true' if wp > 1e-6 else 'false'} for this bracket)")
    x, g = sp.symbols("x g", real=True)
    d = sp.simplify(_br(x, g, sp, True) - _br(x, g, sp, False))
    report("N2", d == 0, f"far-field bracket minus printed bracket, symbolic: {d}")
    mp.mp.dps = 30


def throat_position(g):
    return quad(lambda y: math.exp(2 * g * (math.atan(g) - math.atan(y))), 0.0, g, limit=500)[0]


def TeX_What(z, g):
    """W-hat(z) of [GGS09a] Sec. IV C (phantom_lin.tex lines 1021-1031), retyped, mpmath."""
    z, g = mp.mpf(z), mp.mpf(g)
    ah = mp.exp(2 * g * (mp.atan(g * z) - mp.atan(g)))
    Fz = g * g / (1 + g * g) * ((1 + g * g * z) / (1 + g * g * z * z) + g * mp.atan(g * z))
    F = 1 + g / (1 + g * g) * (1 + g * g * z) * mp.atan(g * z)
    return ah ** 2 * (-3 * g * g / (1 + g * g * z * z) + 3 * g ** 4 * (z - 1) ** 2 / (1 + g * g * z * z) ** 2 + 2 * (Fz / F) ** 2
                      - 4 * g * g / (1 + g * g * z * z) * Fz / F + 4 * g * g / (1 + g * g) * g * g * (z - 1) / (F * (1 + g * g * z * z) ** 2))


def check_limit_potential():
    """N3a: the printed W-hat equals g^2 W-bar(g z); N3b: it tends to the printed limit potential."""
    mp.mp.dps = 50
    worst = 0.0
    for g in (2.5, 10.0, 100.0):
        for z in (0.4, 1.0, 2.0, 5.0):
            x = g * z
            ref = TeX_What(z, g)
            gm, xm = mp.mpf(g), mp.mpf(x)
            a2 = mp.exp(4 * gm * (mp.atan(xm) - mp.atan(gm)))
            mine = gm * gm * a2 * _br(xm, gm, mp, False)
            worst = max(worst, float(abs(mine / ref - 1)))
    report("N3a", worst <= 1e-12, f"[GGS09a] W-hat(z) (TeX, retyped) against gamma_1^2 W-bar(gamma_1 z) of the code, 3 gamma_1 x 4 z, worst rel {worst:.2e}")
    zs = (0.4, 0.8, 1.0, 1.5, 3.0, 8.0)
    c = 9 if mut("lelimit") else 10
    ds = []
    for g in (1e2, 1e3, 1e4, 1e5, 1e6):
        d = 0.0
        for z in zs:
            zm = mp.mpf(z)
            lim = mp.exp(4 * (1 - 1 / zm)) * (2 / zm ** 2 - c / zm ** 3 + 3 / zm ** 4)
            d = max(d, float(abs(TeX_What(z, g) / lim - 1)))
        ds.append(d)
    ok = ds[-1] <= 1e-3 and all(ds[i + 1] <= ds[i] / 5 for i in range(len(ds) - 1))
    report("N3b", ok, f"W-hat(z; gamma_1) against the printed limit potential, z in {zs}: max rel diff at gamma_1 = 1e2 ... 1e6: "
                      f"{', '.join(f'{v:.1e}' for v in ds)} (each decade must fall by 5, last <= 1e-3)")
    mp.mp.dps = 30


def check_LE0():
    z = sp.symbols("z", positive=True)
    ah = sp.exp(2 * (1 - 1 / z))
    c = 9 if mut("lelimit") else 10
    Wh = sp.exp(4 * (1 - 1 / z)) * (2 / z ** 2 - c / z ** 3 + 3 / z ** 4)
    chi0 = (z - 1) / z ** 2 * sp.exp(1 - 1 / z)
    res = sp.simplify(-ah * sp.diff(ah * sp.diff(chi0, z), z) + Wh * chi0)
    report("LE0", res == 0, f"the zero mode (2/pi)(z-1)/z^2 e^(1-1/z) of [GGS09a] solves the printed limit equation with beta-hat = 0, "
                            f"symbolic residual {res}")


def check_U3():
    worst = 0.0
    for g in (0.5, 2.5, 10.0):
        b = 1.7
        worst = max(worst, rel(T_of2(b, g), T_of(beta_orig(b, g), g)))
    report("U3", worst <= 1e-13, f"T of the new path (sc/(sqrt(1+g^2) b)) against T_of(lam b/sc), worst rel {worst:.2e}")


def check_planted2(gammas, L):
    for g in gammas:
        base = sysd_ggs(g)
        spl = u_spline_ggs(g, base["x0"])
        for ell in (0.9, 1.6):
            exact = ell * (1.01 if mut("planted") else 1.0)
            nb = sum(1 for n in range(0, 5) if ell - n > 0)
            sd = sysd_planted(base, spl, ell)
            try:
                roots = gshoot(sd, L=L, scan_n=30, b_hi=math.sqrt(ell * (ell + 1)) + 0.01)
                ok = len(roots) == nb and rel(roots[0], exact) < 1e-9
                report("P-A2", ok, f"gamma_1 = {g}, ell = {ell}: A beta' = {roots[0]:.12f} (exact {exact}), rel {rel(roots[0], exact):.1e}, "
                                   f"{len(roots)} bound states (exact {nb})")
            except Exception as exc:  # noqa: BLE001
                report("P-A2", False, f"gamma_1 = {g}, ell = {ell}: A raised {type(exc).__name__}: {exc}")
            try:
                b0, nneg, _ = gmatrix(sd, L=L)
                ok = nneg == nb and rel(b0, exact) < 1e-8
                report("P-B2", ok, f"gamma_1 = {g}, ell = {ell}: B beta' = {b0:.12f}, rel {rel(b0, exact):.1e}, {nneg} negative "
                                   f"eigenvalues (exact {nb})")
            except Exception as exc:  # noqa: BLE001
                report("P-B2", False, f"gamma_1 = {g}, ell = {ell}: B raised {type(exc).__name__}: {exc}")


def check_planted_C2(gammas):
    x, g = sp.symbols("x g")
    bt = 1.3
    at = sp.atan(x)
    Cstar = sp.exp(g * at) * sp.exp(-(x - g) ** 2 / 2)
    Cs_, Cp_ = sp.symbols("Cs Cp")
    R2 = odeC2_rhs(x, Cs_, Cp_, bt, g, ns=sp)
    S = sp.diff(Cstar, x, 2) - R2.subs({Cs_: Cstar, Cp_: sp.diff(Cstar, x)})
    fC = sp.lambdify((x, g), Cstar, "numpy")
    fCp = sp.lambdify((x, g), sp.diff(Cstar, x), "numpy")
    fS = sp.lambdify((x, g), S, "numpy")
    for gv in gammas:
        r = 0.0 if mut("cpath") else 0.5 * math.sqrt(1 + gv * gv)
        xt = gv
        worst = 0.0
        try:
            with np.errstate(all="ignore"):
                for side in (+1, -1):
                    def rhs(th, y):
                        xx = xt + r * cmath.exp(1j * th)
                        dx = 1j * r * cmath.exp(1j * th)
                        v = [dx * y[1], dx * (odeC2_rhs(xx, y[0], y[1], bt, gv) + complex(fS(xx, gv)))]
                        if not all(cmath.isfinite(z) for z in v):
                            raise ValueError("non-finite right-hand side")
                        return v

                    sol = solve_ivp(rhs, (0.0, side * math.pi), [complex(fC(xt + r, gv)), complex(fCp(xt + r, gv))], method="DOP853",
                                    rtol=1e-12, atol=[1e-16, 1e-16])
                    ex = np.array([complex(fC(xt - r, gv)), complex(fCp(xt - r, gv))])
                    got = np.array([sol.y[0, -1], sol.y[1, -1]])
                    e = float(np.linalg.norm(got - ex) / np.linalg.norm(ex))
                    worst = max(worst, e if np.isfinite(e) else float("inf"))
        except Exception as exc:  # noqa: BLE001
            print(f"      P-C2 integration raised {type(exc).__name__}: {exc}")
            worst = float("inf")
        report("P-C2", worst < 1e-8, f"gamma_1 = {gv}: manufactured solution of the rescaled ODE carried across the throat on both "
                                     f"semicircles, worst rel error {worst:.2e}")


def check_planted_LE(L):
    spl = u_spline_le()
    base = sysd_le()
    for ell in (0.9, 1.6):
        exact = ell * (1.01 if mut("planted") else 1.0)
        nb = sum(1 for n in range(0, 5) if ell - n > 0)
        sd = sysd_planted(base, spl, ell)
        try:
            roots = gshoot(sd, L=L, scan_n=30, b_hi=math.sqrt(ell * (ell + 1)) + 0.01)
            report("P-LA", len(roots) == nb and rel(roots[0], exact) < 1e-9,
                   f"limit-equation geometry, ell = {ell}: A beta = {roots[0]:.12f} (exact {exact}), rel {rel(roots[0], exact):.1e}, "
                   f"{len(roots)} bound states (exact {nb})")
        except Exception as exc:  # noqa: BLE001
            report("P-LA", False, f"ell = {ell}: A raised {type(exc).__name__}: {exc}")
        try:
            b0, nneg, _ = gmatrix(sd, L=L)
            report("P-LB", nneg == nb and rel(b0, exact) < 1e-8,
                   f"limit-equation geometry, ell = {ell}: B beta = {b0:.12f}, rel {rel(b0, exact):.1e}, {nneg} negative eigenvalues (exact {nb})")
        except Exception as exc:  # noqa: BLE001
            report("P-LB", False, f"ell = {ell}: B raised {type(exc).__name__}: {exc}")
        try:
            lam = ell * (ell + 1)

            def Wpl(zs):
                e = np.exp(-2.0 * np.abs(spl(zs)))
                return -lam * 4.0 * e / (1.0 + e) ** 2

            bC, nC, _, _ = le_fv_beta(Wpl)
            report("P-LC", nC == nb and rel(bC, exact) < 1e-6,
                   f"limit-equation geometry, ell = {ell}: finite volume (Richardson) beta = {bC:.10f}, rel {rel(bC, exact):.1e}, "
                   f"{nC} eigenvalues in [-60, 0] (exact {nb})")
        except Exception as exc:  # noqa: BLE001
            report("P-LC", False, f"ell = {ell}: FV raised {type(exc).__name__}: {exc}")


def run_rows(L):
    """A at all 32 values of Table I (both rows), the checks N1b R2 G5b R3, the table lines."""
    BN = {}
    worst_pr, skipped_pr = 0.0, 0
    rows = list(TABLE) + list(TABLE2)
    for g in rows:
        sd = sysd_ggs(g)
        try:
            xR, xL = x_ends(sd, L)
            rhot = throat_position(g) if g > 0 else 0.0
            roots = gshoot(sd, L=L)
            BN[g] = roots[0]
            T = T_of2(roots[0], g)
            row2 = g in TABLE2
            print(f"      gamma_1 = {g:4.1f}: b = beta-hat = {roots[0]:.12f}  beta(lam units) = {beta_orig(roots[0], g):.12g}  T = {T:.9f}  "
                  f"x(-L), x(+L) = {xL:.4g}, {xR:.4g}  rho~_throat = {rhot:.4g} (legacy domain +-30)  roots in scan {len(roots)}", flush=True)
            xs_chk = [xL, xR] + list(x_on(sd, np.array([-L / 2, -2.0, -0.5, 0.5, 2.0, L / 2])))
            ref = wbar_derived_vec(np.array(xs_chk), g)
            for xv, rv in zip(xs_chk, ref):
                if abs(rv) >= 1e-3 * Wscale(float(xv), g):
                    worst_pr = max(worst_pr, abs(Wbar(float(xv), g, False) / rv - 1))
                else:
                    skipped_pr += 1
        except Exception as exc:  # noqa: BLE001
            report("G5b" if g in TABLE2 else "R2", False, f"gamma_1 = {g}: A (centred) raised {type(exc).__name__}: {exc}")
    report("N1b", worst_pr <= 1e-12, f"printed W-bar (double) against the derived W (40 digits) at 8 points of each of the {len(rows)} "
                                    f"centred domains, worst rel {worst_pr:.2e} ({skipped_pr} points skipped near a zero crossing)")
    # R2: the first row through the new path
    worst60, worst30, row = 0.0, 0.0, []
    for g in TABLE:
        if g not in BN:
            row.append((g, float("nan"), float("nan")))
            continue
        b = beta_orig(BN[g], g)
        e60, e30 = rel(b, FIRST_ROW_BETA_L60[g]), rel(b, FIRST_ROW_BETA[g])
        worst60 = max(worst60, e60)
        if g <= 1.8:
            worst30 = max(worst30, e30)
        row.append((g, e60, e30))
    for g, e60, e30 in row:
        print(f"      first row via the new path: gamma_1 = {g:3.1f}  rel to #151 (L=60) {e60:.2e}  rel to #151 as run (L=30) {e30:.2e}")
    ok = len([g for g in TABLE if g in BN]) == len(TABLE)
    report("R2", ok and worst60 <= 1e-9 and worst30 <= 1e-9,
           f"first row, 16 values, new path against #151's beta: converged (L = 60) reference worst rel {worst60:.2e}; as-run (L = 30) reference "
           f"at gamma_1 <= 1.8 worst rel {worst30:.2e}; at 2.0 the as-run value differs by {row[-1][2]:.2e} (its domain truncation)")
    # G5b: the second row against the printed table
    worst, lines = 0.0, []
    for g, tv in TABLE2.items():
        if g in BN:
            T = T_of2(BN[g], g)
            lines.append((g, T, tv))
            worst = max(worst, abs(T - tv))
    for g, T, tv in lines:
        print(f"      Table I, second row: gamma_1 = {g:4.1f}  T = {T:.6f}  printed {tv:.3f}  diff {T - tv:+.6f}")
    report("G5b", len(lines) == len(TABLE2) and worst <= 0.001, f"16 values of Table I second row, worst |T - printed| = {worst:.6f}")
    return BN


def run_R3(L, BN):
    worst = 0.0
    for g in (0.5, 2.0, 2.5, 10.0):
        try:
            b2 = gshoot(sysd_ggs(g, far=False), L=L)[0]
            worst = max(worst, rel(b2, BN[g]) if g in BN else float("inf"))
        except Exception as exc:  # noqa: BLE001
            print(f"      R3 gamma_1 = {g}: raised {type(exc).__name__}: {exc}")
            worst = float("inf")
    report("R3", worst <= 1e-10, f"A with the far-field form against A with the printed bracket only, gamma_1 in (0.5, 2, 2.5, 10): worst rel {worst:.2e}")


def run_G1b(L, BN):
    worst, n = 0.0, 0
    try:
        for g in list(TABLE2) + [0.5, 1.0, 2.0]:
            bB, nneg, ev = gmatrix(sysd_ggs(g), L=L)
            e = rel(bB, BN[g]) if g in BN else float("inf")
            worst = max(worst, e)
            n += 1
            print(f"      gamma_1 = {g:4.1f}: b_A = {BN.get(g, float('nan')):.12f}, b_B = {bB:.12f}, rel {e:.1e}, negative eigenvalues {nneg}", flush=True)
        report("G1b", worst <= 1e-7, f"A against B (derived W) at {n} values (16 of the second row, 0.5, 1, 2), worst rel {worst:.2e}")
    except Exception as exc:  # noqa: BLE001
        report("G1b", False, f"B raised {type(exc).__name__}: {exc}")


def run_G2b(L, BN):
    for g in (2.5, 4.0, 6.0, 8.0, 10.0):
        try:
            sc = sc_of(g)
            bhi = BN[g] * 1.6 if g in BN else 3.0
            rp = method_C2(g, L=L, side=+1, scan=(0.2, bhi, 60))
            rm = method_C2(g, L=L, side=-1, scan=(0.2, bhi, 60))
            dets = [C2_det(b, g, L, 0.5, sd_, 1e-12)[1] for b in (0.8 * rp[0], 1.1 * rp[0]) for sd_ in (+1, -1)]
            ok = rel(rp[0], BN[g]) < 1e-6 and rel(rp[0], rm[0]) < 1e-9
            report("G2b", ok, f"gamma_1 = {g}: C beta' = {rp[0]:.12f} (upper) {rm[0]:.12f} (lower), roots in [0.2, {bhi:.2f}]: "
                              f"{[round(v, 6) for v in rp]}; rel to A {rel(rp[0], BN[g]):.2e}; max |Im det| {max(dets):.1e}")
        except Exception as exc:  # noqa: BLE001
            report("G2b", False, f"gamma_1 = {g}: method C or its comparison failed: {type(exc).__name__}: {exc}")


def run_G7b(L):
    for g in (2.5, 5.0, 10.0):
        try:
            sd = sysd_ggs(g, count=mut("count"))
            roots = gshoot(sd, L=L)
            bB, nB, evB = gmatrix(sd, L=L)
            report("G7b", len(roots) == 1 and nB == 1, f"gamma_1 = {g}: A finds {len(roots)} root(s) in the scan, B has {nB} negative eigenvalue(s) "
                                                       f"(lowest {evB[0] if len(evB) else float('nan'):.6f}; the next is >= 0)")
        except Exception as exc:  # noqa: BLE001
            report("G7b", False, f"gamma_1 = {g}: raised {type(exc).__name__}: {exc}")


def run_D1b(L):
    for g in (2.5, 5.0, 10.0):
        try:
            sd = sysd_ggs(g)
            base = gshoot(sd, L=L)[0]
            vA = [gshoot(sd, L=L, rtol=1e-11)[0], gshoot(sd, L=L + 8)[0]]
            sA = max(rel(v, base) for v in vA)
            bB = gmatrix(sd, L=L)[0]
            vB = [gmatrix(sd, L=L, N=8000)[0], gmatrix(sd, L=L + 8)[0]]
            sB = max(rel(v, bB) for v in vB)
            bC = method_C2(g, L=L, scan=(0.2, base * 1.6, 60))[0]
            vC = [method_C2(g, L=L, rfrac=0.3, centre=bC)[0], method_C2(g, L=L + 8, rfrac=0.7, centre=bC)[0]]
            sC = max(rel(v, bC) for v in vC)
            ok = max(sA, sB, sC) < 1e-7
            report("D1b", ok, f"gamma_1 = {g}: max relative spread of beta': A (rtol 1e-13 -> 1e-11, L' +8) {sA:.1e}; B (N 4000 -> 8000, L' +8) "
                              f"{sB:.1e}; C (r/dist 0.5 -> 0.3, 0.7; L' +8) {sC:.1e}")
            print(f"      gamma_1 = {g}: b_A = {base:.12f}, b_B = {bB:.12f}, b_C = {bC:.12f}; error bar of T from the spreads: "
                  f"{max(sA, sB, sC) * T_of2(base, g):.1e}", flush=True)
        except Exception as exc:  # noqa: BLE001
            report("D1b", False, f"gamma_1 = {g}: raised {type(exc).__name__}: {exc}")


def run_LE(L):
    """LE3, LE4, LE5, LE6 and the value of 1/beta-hat."""
    out = {}
    try:
        sd = sysd_le()
        rA = gshoot(sd, L=L)
        bA = rA[0]
        bB, nB, evB = gmatrix(sd, L=L)
        bC, nC, bC1, bC2 = le_fv_beta(le_Wv)
        out.update(A=bA, B=bB, C=bC)
        ok = rel(bA, bB) < 1e-8 and rel(bA, bC) < 1e-6
        report("LE3", ok, f"limit equation: A beta-hat = {bA:.12f}, B = {bB:.12f} (rel {rel(bA, bB):.1e}), finite volume = {bC:.10f} "
                          f"(N = 8000: {bC1:.10f}, 16000: {bC2:.10f}; Richardson rel to A {rel(bA, bC):.1e})")
        sdc = sysd_le(count=mut("count"))
        rc = gshoot(sdc, L=L)
        bcB, ncB, _ = gmatrix(sdc, L=L)
        report("LE4", len(rc) == 1 and ncB == 1 and nB == 1 and nC == 1,
               f"limit equation: A finds {len(rc)} root(s) in the scan, B has {nB} negative eigenvalue(s) ({ncB} with the extra-well control "
               f"geometry), the finite volume {nC}; the zero mode (one node) forces at least one")
        T = 1.0 / bA
        out["T"] = T
        print(f"      limit equation: 1/beta-hat = {T:.9f} (A), {1 / bB:.9f} (B), {1 / bC:.9f} (finite volume); paper {LE_PAPER} (Sec. IV C), "
              f"table {LE_TABLE}", flush=True)
        report("LE5", abs(T - LE_TABLE) <= 0.001, f"1/beta-hat = {T:.6f} against the table's large-gamma_1 value {LE_TABLE}: {T - LE_TABLE:+.6f} (<= 0.001)")
        report("LE5", abs(T - LE_PAPER) < 0.0005, f"1/beta-hat = {T:.6f} against the paper's 1/beta-hat = {LE_PAPER}: {T - LE_PAPER:+.6f} (must round to it, < 0.0005)")
    except Exception as exc:  # noqa: BLE001
        report("LE3", False, f"limit equation raised {type(exc).__name__}: {exc}")
    return out


def run_LE6(L, BN, LE):
    """the finite-gamma_1 problem approaches the limit equation: beta-hat_g = gamma_1 beta-bar = b."""
    if "A" not in LE:
        report("LE6", False, "no limit-equation value")
        return
    bl = LE["A"]
    gs = (10.0, 20.0, 40.0, 80.0, 160.0)
    vals = {}
    try:
        for g in gs:
            vals[g] = BN[g] if g in BN else gshoot(sysd_ggs(g), L=L)[0]
            print(f"      gamma_1 = {g:5.0f}: beta-hat = {vals[g]:.10f}  (limit {bl:.10f})  T = {T_of2(vals[g], g):.8f}  gap {vals[g] - bl:+.3e}", flush=True)
        gaps = [abs(vals[g] - bl) for g in gs]
        dec = all(gaps[i + 1] < gaps[i] / 1.5 for i in range(len(gs) - 1))
        ok = dec and gaps[-1] / bl <= 3e-3
        T80, T160 = T_of2(vals[80.0], 80.0), T_of2(vals[160.0], 160.0)
        print(f"      Richardson in 1/gamma_1 from 80, 160 (reported only): T_inf = {2 * T160 - T80:.6f}; 1/beta-hat = {1 / bl:.6f}")
        report("LE6", ok, f"beta-hat at gamma_1 = 10, 20, 40, 80, 160 against the limit equation's: gaps {', '.join(f'{v:.2e}' for v in gaps)}; each doubling "
                          f"must cut the gap by 1.5 and the last must be <= 0.3 %: {gaps[-1] / bl:.2e}")
    except Exception as exc:  # noqa: BLE001
        report("LE6", False, f"raised {type(exc).__name__}: {exc}")


def run_diag_legacy():
    for g in (2.5, 3.0):
        try:
            r = method_A(g, W_A, L=30.0, lam=lam_of(g))
            print(f"      DIAG L0: the legacy path at gamma_1 = {g} (domain +-30, origin x = 0, throat at rho~ = {throat_position(g):.1f}): "
                  f"{len(r)} root(s) {[f'{v:.6g}' for v in r]}")
        except Exception as exc:  # noqa: BLE001
            print(f"      DIAG L0: the legacy path at gamma_1 = {g} raised {type(exc).__name__}: {exc}")


def run_new():
    L = LN()
    print("--- second row, limit equation (#163)")
    check_far_field()
    check_limit_potential()
    check_LE0()
    check_U3()
    check_planted2((2.5, 5.0, 10.0), L)
    check_planted_C2((2.5, 5.0, 10.0))
    check_planted_LE(L)
    run_diag_legacy()
    BN = run_rows(L)
    run_R3(L, BN)
    run_G1b(L, BN)
    run_G2b(L, BN)
    run_G7b(L)
    run_D1b(L)
    LE = run_LE(L)
    run_LE6(L, BN, LE)


OLD_IDS = {"U1", "U2", "S1", "L1", "W1", "P-A", "P-B", "P-C", "G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "D1", "R1"}
NEW_IDS = {"N1a", "N1b", "N2", "N3a", "N3b", "U3", "LE0", "P-A2", "P-B2", "P-C2", "P-LA", "P-LB", "P-LC", "R2", "R3", "G5b", "G1b",
           "G2b", "G7b", "D1b", "LE3", "LE4", "LE5", "LE6"}


def main():
    global MUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--mutate", default=None, choices=sorted(MUTATIONS))
    ap.add_argument("--gamma", default=None, help="0: only the gamma_1 = 0 control and its supporting checks")
    ap.add_argument("--section", default="all", choices=["all", "old", "new"], help="old: #151's checks; new: #163's (second row, limit)")
    args = ap.parse_args()
    MUT = args.mutate
    only0 = args.gamma == "0"
    sect = "old" if only0 else args.section
    L = 4.0 if mut("domain") else 30.0
    print(f"ellis_efold.py {' '.join(sys.argv[1:])}  python {platform.python_version()} numpy {np.__version__} "
          f"sympy {sp.__version__} mpmath {mp.__version__} on {platform.machine()}")
    print(f"mutation: {MUT}" + (f" ({MUTATIONS[MUT][0]})" if MUT else ""))
    if sect in ("all", "old"):
        check_units()
        check_T_formula()
        if not only0:
            check_background()
            check_linearization()
        check_W1()
        gam = [0.0] if only0 else [0.0, 0.5, 1.0]
        check_planted(gam if only0 else [0.0, 0.5, 1.0], L)
        if not only0:
            check_planted_C(gam)
        gset = [0.0] if only0 else [0.0, 0.5, 1.0, 2.0]
        run_G(gset, L, full=not only0)
        run_D1(gset, L)
    if sect in ("all", "new"):
        run_new()
    print(f"checks: {sum(len(v) for v in RESULTS.values())} run, failed ids: {FAILED}  ({time.time() - T_START:.0f} s)")
    if MUT:
        target = MUTATIONS[MUT][1]
        in_scope = [c for c in target if (c in OLD_IDS and sect in ("all", "old")) or (c in NEW_IDS and sect in ("all", "new"))]
        fired = [c for c in in_scope if c in FAILED]
        silent = [c for c in in_scope if c not in FAILED]
        outside = [c for c in target if c not in in_scope]
        other = [c for c in FAILED if c not in target]
        print(f"control {MUT}: target {target}: fired {fired}, SILENT {silent} (a target that never ran is silent); not in this "
              f"section {outside}; failed outside the target list {other}")
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()

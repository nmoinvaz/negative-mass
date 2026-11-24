#!/usr/bin/env python3
"""Linear e-folding time of the asymmetric Ellis family (#151), against [GGS09a] Table I and [SN26] Sec. IV B.

    python3 ellis_efold.py                  all checks; exit 1 if any check fails
    python3 ellis_efold.py --gamma=0        the gamma_1 = 0 control and its supporting checks only
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
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.linalg import eig_banded
from scipy.optimize import brentq

MUT = None
FAILED = []
RESULTS = {}
T_START = time.time()

# mutation -> (what it changes, the checks that must fire on it); the same table as README section 5
MUTATIONS = {
    "potential": ("W -> W + 0.1/(1+x^2)^2 in A and B", ["G0", "G2", "G3", "G4", "G5", "G6"]),
    "printedW": ("the -4g/(1+x^2) Fx/F term of the printed W gets 3g", ["W1", "G1", "G2", "G3", "G4", "G5"]),
    "matrix": ("B's centre weight -49/18 -> -2.7", ["P-B", "G1", "G0"]),
    "planted": ("planted exact beta = ell -> 1.01 ell", ["P-A", "P-B"]),
    "xmap": ("A's x_rho = e^{-2 g atan x}", ["P-A", "G1", "G2", "G3", "G4", "G5"]),
    "units": ("gamma_0 = +pi g/2 in the [SN26] -> [GGS09a] map", ["U1"]),
    "tfactor": ("T drops e^{2 g atan g}", ["U2", "G3", "G4", "G5", "G6"]),
    "count": ("G7 uses W + a second well", ["G7"]),
    "domain": ("A and B use L = 4", ["D1", "G0", "G1"]),
    "deltaphi": ("the 2 in 2 e^{2(d-c)}/c' of method C becomes 1", ["L1", "G2"]),
    "cpath": ("method C's semicircle radius -> 0", ["P-C", "G2"]),
    "background": ("kappa Phi1^2 = -2(1+g^2) -> -2(1+2 g^2) in S1", ["S1"]),
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
        return [y[1], odeC_rhs(x, y[0], y[1], beta, g, ns=math)]

    sR = solve_ivp(real_rhs, (xR, xt + r), [1.0, -beta * math.exp(-2 * g * math.atan(xR))], method="DOP853", rtol=rtol,
                   atol=[1e-16, 1e-16])
    sL = solve_ivp(real_rhs, (xL, xt - r), [1.0, beta * math.exp(-2 * g * math.atan(xL))], method="DOP853", rtol=rtol,
                   atol=[1e-16, 1e-16])

    def arc_rhs(th, y):
        x = xt + r * cmath.exp(1j * th)
        return [1j * r * cmath.exp(1j * th) * y[1], 1j * r * cmath.exp(1j * th) * odeC_rhs(x, y[0], y[1], beta, g)]

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
                T = T_of(b, g)
                row.append((g, T, tv))
                worst = max(worst, abs(T - tv))
            for g, T, tv in row:
                print(f"      Table I: gamma_1 = {g:3.1f}  T = {T:.6f}  printed {tv:.3f}  diff {T - tv:+.6f}")
            report("G5", worst <= 0.001, f"16 values of Table I, worst |T - printed| = {worst:.6f}")
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


def main():
    global MUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--mutate", default=None, choices=sorted(MUTATIONS))
    ap.add_argument("--gamma", default=None, help="0: only the gamma_1 = 0 control and its supporting checks")
    args = ap.parse_args()
    MUT = args.mutate
    only0 = args.gamma == "0"
    L = 4.0 if mut("domain") else 30.0
    print(f"ellis_efold.py {' '.join(sys.argv[1:])}  python {platform.python_version()} numpy {np.__version__} "
          f"sympy {sp.__version__} mpmath {mp.__version__} on {platform.machine()}")
    print(f"mutation: {MUT}" + (f" ({MUTATIONS[MUT][0]})" if MUT else ""))
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
    print(f"checks: {sum(len(v) for v in RESULTS.values())} run, failed ids: {FAILED}  ({time.time() - T_START:.0f} s)")
    if MUT:
        target = MUTATIONS[MUT][1]
        fired = [c for c in target if c in FAILED]
        silent = [c for c in target if c not in FAILED]
        other = [c for c in FAILED if c not in target]
        print(f"control {MUT}: target {target}: fired {fired}, SILENT {silent}; failed outside the target list {other}")
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()

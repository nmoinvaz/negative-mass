"""Spherically symmetric Einstein gravity in D = n + 2 dimensions with massless scalars of either sign, in double-null
coordinates (#39). Design: design.md. Equations: derive.py (checked symbolically).

    ds^2 = -2 e^{-f} dx+ dx- + r^2 dOmega_n^2,   R_ab = sum_k eps_k d_a phi_k d_b phi_k   (8 pi G = 1, eps = +1 normal, -1 ghost)
    r_{+-} = -(n-1)(r_+ r_- + e^{-f}/2)/r
    f_{+-} = (n/r) r_{+-} + sum eps phi_+ phi_-
    phi_{+-} = -n (r_+ phi_- + r_- phi_+)/(2 r)

The grid is the rectangle [x0p, x0p + N h] x [x0m, x0m + N h] with data on the two null rays through the corner S (row j = 0 is
x- = x0m, column i = 0 is x+ = x0p). Fields are stored as G[q, i, j], q = 0 r, 1 f, 2.. scalars. Two kernels of the same
second-order predictor-corrector diamond update: `numba` (row sweep, production) and `numpy` (anti-diagonal sweep, reference;
needs no numba). Mutations (`MUT`) break one input on purpose; checks.py names the checks that must fire.
"""
import math

import numpy as np

import os
try:
    if os.environ.get("DN_NO_NUMBA"):
        raise ImportError("numba disabled by DN_NO_NUMBA")
    import numba as nb
    HAVE_NUMBA = True
except Exception:          # the numpy kernel must work without numba
    nb = None
    HAVE_NUMBA = False

SQ2 = math.sqrt(2.0)
MUT = None                 # set by checks.py --mutate=NAME; names in checks.py MUTATIONS


def mut(name):
    return MUT == name


# --------------------------------------------------------------------------------------------------------------------
# kernels
# --------------------------------------------------------------------------------------------------------------------
def kernel_params(n, eps, ncorr=1):
    """[n, th, sc, pms, ncorr]: coefficients of the sphere term, the scalar equation and the f source sign; mutable."""
    th = n if mut("theta_coef") else n - 1.0
    sc = n - 1.0 if mut("scalar_coef") else float(n)
    pms = -1.0 if mut("ghost_sign") else 1.0
    return np.array([float(n), th, sc, pms, 0.0 if mut("no_corrector") else float(ncorr)])


def F_np(c, up, um, par, eps):
    """Vectorized U_{+-} for rows of c, up, um (nf, L)."""
    n, th, sc, pms = par[0], par[1], par[2], par[3]
    nf = c.shape[0]
    emf = np.exp(-c[1])
    rpm = -th * (up[0] * um[0] + 0.5 * emf) / c[0]
    spm = 0.0
    for k in range(nf - 2):
        spm = spm + eps[k] * up[2 + k] * um[2 + k]
    g = np.empty_like(c)
    g[0] = rpm
    g[1] = (n / c[0]) * rpm + pms * spm
    kk = -sc / (2.0 * c[0])
    for k in range(2, nf):
        g[k] = kk * (up[0] * um[k] + um[0] * up[k])
    return g


def cell_np(S, E, W, h, par, eps):
    """Update N from S, E, W (nf, L). The numpy kernel; the same arithmetic as row_step."""
    h2 = h * h
    ncorr = int(par[4])
    c = 0.5 * (E + W)
    up = (E - S) / h
    um = (W - S) / h
    Nn = E + W - S + h2 * F_np(c, up, um, par, eps)
    for _ in range(ncorr):
        c = 0.25 * (S + E + W + Nn)
        up = (E + Nn - S - W) / (2 * h)
        um = (W + Nn - S - E) / (2 * h)
        Nn = E + W - S + h2 * F_np(c, up, um, par, eps)
    return Nn


if HAVE_NUMBA:
    @nb.njit(cache=True)
    def _F(c, up, um, par, eps, g):
        n = par[0]
        th = par[1]
        sc = par[2]
        pms = par[3]
        nf = c.shape[0]
        emf = math.exp(-c[1])
        rpm = -th * (up[0] * um[0] + 0.5 * emf) / c[0]
        spm = 0.0
        for k in range(nf - 2):
            spm = spm + eps[k] * up[2 + k] * um[2 + k]
        g[0] = rpm
        g[1] = (n / c[0]) * rpm + pms * spm
        kk = -sc / (2.0 * c[0])
        for k in range(2, nf):
            g[k] = kk * (up[0] * um[k] + um[0] * up[k])

    @nb.njit(cache=True)
    def row_step(old, new, dead_old, dead_new, h, par, eps, rcut, fcut, excise):
        """new[:, 0] and dead_new[0] are set by the caller; fills new[:, 1:] from old and new[:, :-1]."""
        nf, n1 = old.shape
        h2 = h * h
        ncorr = int(par[4])
        c = np.empty(nf)
        up = np.empty(nf)
        um = np.empty(nf)
        nh = np.empty(nf)
        g = np.empty(nf)
        for i in range(n1 - 1):
            if dead_old[i] or dead_old[i + 1] or dead_new[i]:
                dead_new[i + 1] = 1
                for q in range(nf):
                    new[q, i + 1] = old[q, i + 1]
                continue
            for q in range(nf):
                S = old[q, i]
                E = old[q, i + 1]
                W = new[q, i]
                c[q] = 0.5 * (E + W)
                up[q] = (E - S) / h
                um[q] = (W - S) / h
            _F(c, up, um, par, eps, g)
            for q in range(nf):
                nh[q] = old[q, i + 1] + new[q, i] - old[q, i] + h2 * g[q]
            for _ in range(ncorr):
                for q in range(nf):
                    S = old[q, i]
                    E = old[q, i + 1]
                    W = new[q, i]
                    c[q] = 0.25 * (S + E + W + nh[q])
                    up[q] = (E + nh[q] - S - W) / (2 * h)
                    um[q] = (W + nh[q] - S - E) / (2 * h)
                _F(c, up, um, par, eps, g)
                for q in range(nf):
                    nh[q] = old[q, i + 1] + new[q, i] - old[q, i] + h2 * g[q]
            bad = False
            for q in range(nf):
                if not (nh[q] == nh[q]) or abs(nh[q]) > 1e300:
                    bad = True
            if excise and (bad or nh[0] < rcut or math.exp(-nh[1]) > fcut):
                dead_new[i + 1] = 1
                for q in range(nf):
                    new[q, i + 1] = old[q, i + 1]
            else:
                dead_new[i + 1] = 0
                for q in range(nf):
                    new[q, i + 1] = nh[q]


# --------------------------------------------------------------------------------------------------------------------
# initial data
# --------------------------------------------------------------------------------------------------------------------
class Data:
    """Characteristic data: Sp[q, i] on x- = x0m (j = 0), Sm[q, j] on x+ = x0p (i = 0), the model and the grid."""

    def __init__(self, n, eps, N, h, Sp, Sm, x0p=0.0, x0m=0.0, name=""):
        self.n, self.eps, self.N, self.h = n, np.array(eps, float), N, h
        self.Sp, self.Sm, self.x0p, self.x0m, self.name = Sp, Sm, x0p, x0m, name
        self.nf = 2 + len(eps)


_ELLIS = {}


def ellis_dense(n, lmax):
    """Dense solution of the static Ellis wormhole, a = 1: y = (r, r', phi), r'' = (n-1)/r^(2n-1), phi' = q/r^n."""
    from scipy.integrate import solve_ivp
    key = (n, math.ceil(lmax))
    if key in _ELLIS:
        return _ELLIS[key]
    qq = math.sqrt(n * (n - 1.0))

    def rhs(l, y):
        return [y[1], (n - 1.0) / y[0] ** (2 * n - 1), qq / y[0] ** n]
    sol = solve_ivp(rhs, [0.0, float(key[1]) + 2.0], [1.0, 0.0, 0.0], method="DOP853", rtol=1e-13, atol=1e-15, dense_output=True)
    _ELLIS[key] = (sol.sol, qq)
    return _ELLIS[key]


def _pulse(s, A, xc, w, kind="phi"):
    """(value, derivative) of a Gaussian pulse that vanishes at s = 0. kind 'phi': A (g - g(0)), g = exp(-((s-xc)/w)^2);
    kind 'dphi' (the form of [SH02] (3.1), a Gaussian in the derivative): d phi/ds = A g, phi = int_0^s."""
    from scipy.special import erf
    g = np.exp(-((s - xc) / w) ** 2)
    if kind == "dphi":
        return A * w * math.sqrt(math.pi) / 2 * (erf((s - xc) / w) - erf(-xc / w)), A * g
    g0 = math.exp(-((0.0 - xc) / w) ** 2)
    return A * (g - g0), A * (-2.0 * (s - xc) / w ** 2) * g


def ellis_rays(n, N, h, pulses=(), eps=(-1.0,), symmetric=False):
    """Static Ellis wormhole (field 0 is the ghost, eps[0] = -1, a = 1) plus Gaussian pulses on either ray.

    pulses: dicts {field, ray ('+' or '-'), A, xc, w, kind ('phi' default, or 'dphi')}. symmetric=True mirrors each '+' pulse onto '-' with the parity of the
    field (-1 for field 0, which has the odd static profile; +1 for others). The corner is the throat: r = 1, r_s = 0, f = 0.
    The constraint C = r'' + (S/n) r = 0 is solved by RK4 on the grid with the field derivative evaluated analytically."""
    eps = list(eps)
    k = len(eps)
    pulses = [dict(p) for p in pulses]
    if symmetric:
        for p in list(pulses):
            if p["ray"] == "+":
                q = dict(p)
                q["ray"] = "-"
                q["A"] = p["A"] * (-1.0 if p["field"] == 0 else 1.0)
                pulses.append(q)
    sol, q = ellis_dense(n, N * h / SQ2 + 1.0)
    q = q * (1.1 if mut("ellis_charge") else 1.0)       # the data's charge; r(l) does not depend on it
    cs = -1.0 if mut("data_constraint_sign") else 1.0

    def phi_s(s, ray, field):
        """d phi_field / ds on the ray at s (array)."""
        out = np.zeros_like(s)
        if field == 0:
            l = s / SQ2
            r = sol(l)[0]
            out = out + (1.0 if ray == "+" else -1.0) * q / r ** n / SQ2
        for p in pulses:
            if p["ray"] == ray and p["field"] == field:
                out = out + _pulse(s, p["A"], p["xc"], p["w"], p.get("kind", "phi"))[1]
        return out

    def phi_val(s, ray, field):
        out = np.zeros_like(s)
        if field == 0:
            out = out + (1.0 if ray == "+" else -1.0) * sol(s / SQ2)[2]
        for p in pulses:
            if p["ray"] == ray and p["field"] == field:
                out = out + _pulse(s, p["A"], p["xc"], p["w"], p.get("kind", "phi"))[0]
        return out

    rays = []
    nodes = np.arange(N + 1) * h
    for ray in ("+", "-"):
        def S_of(s):
            return sum(eps[f] * phi_s(s, ray, f) ** 2 for f in range(k))
        S0, Sh, S1 = S_of(nodes[:-1]), S_of(nodes[:-1] + 0.5 * h), S_of(nodes[1:])
        r = np.empty(N + 1)
        v = np.empty(N + 1)
        r[0], v[0] = 1.0, 0.0
        for i in range(N):
            y = (r[i], v[i])

            def d(y, S):
                return (y[1], -cs * (S / n) * y[0])
            k1 = d(y, S0[i])
            y2 = (y[0] + 0.5 * h * k1[0], y[1] + 0.5 * h * k1[1])
            k2 = d(y2, Sh[i])
            y3 = (y[0] + 0.5 * h * k2[0], y[1] + 0.5 * h * k2[1])
            k3 = d(y3, Sh[i])
            y4 = (y[0] + h * k3[0], y[1] + h * k3[1])
            k4 = d(y4, S1[i])
            r[i + 1] = y[0] + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            v[i + 1] = y[1] + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        U = np.zeros((2 + k, N + 1))
        U[0] = r
        for f in range(k):
            U[2 + f] = phi_val(nodes, ray, f)
        rays.append(U)
    return Data(n, eps, N, h, rays[0], rays[1], name="ellis")


_KR = {}


def kruskal_tools(n):
    """Tangherlini with r_h = 1 (M = 1/2) in Kruskal null coordinates U, V (UV = G(r) = -(r-1) e^{2 kappa R(r)}, kappa = (n-1)/2):
    R(r) = int_1^r [1/F - 1/(2 kappa (s-1))] ds is regular, F = 1 - r^{-(n-1)}. Returns functions R(r), g(r) = R'(r), and
    the regular ratio F/(r-1)."""
    if n in _KR:
        return _KR[n]
    from scipy.integrate import solve_ivp
    from scipy.interpolate import CubicSpline
    m = n - 1.0
    K = 90
    c = np.zeros(K + 1)
    c[0] = 1.0
    for k in range(1, K + 1):
        c[k] = c[k - 1] * (-m - k + 1) / k          # binom(-m, k)

    def F_over_x(x):
        x = np.asarray(x, float)
        out = np.empty_like(x)
        sm = np.abs(x) < 0.3
        xs = x[sm]
        acc = np.zeros_like(xs)
        for k in range(K, 0, -1):
            acc = acc * xs + (-c[k])
        out[sm] = acc                                # -sum_{k>=1} c_k x^{k-1}
        xl = x[~sm]
        out[~sm] = (1.0 - (1.0 + xl) ** (-m)) / xl
        return out

    def g(r):
        r = np.asarray(r, float)
        x = r - 1.0
        out = np.empty_like(x)
        sm = np.abs(x) < 0.3
        xs = x[sm]
        num = np.zeros_like(xs)
        for k in range(K, 1, -1):
            num = num * xs + c[k]
        num = num * xs * xs                          # sum_{k>=2} c_k x^k
        Fx = F_over_x(xs) * xs
        out[sm] = num / (m * xs * Fx) if np.all(xs != 0) else np.where(xs != 0, num / (m * np.where(xs != 0, xs, 1) * np.where(xs != 0, Fx, 1)), (m + 1) / (2 * m))
        xl = x[~sm]
        Fl = 1.0 - (1.0 + xl) ** (-m)
        out[~sm] = 1.0 / Fl - 1.0 / (m * xl)
        return out

    grid_hi = np.linspace(1.0, 6.0, 5001)
    grid_lo = np.linspace(1.0, 0.02, 5001)
    sh = solve_ivp(lambda r, y: g(np.array([r]))[:1], [1.0, 6.0], [0.0], method="DOP853", rtol=1e-13, atol=1e-15, t_eval=grid_hi)
    sl = solve_ivp(lambda r, y: g(np.array([r]))[:1], [1.0, 0.02], [0.0], method="DOP853", rtol=1e-13, atol=1e-15, t_eval=grid_lo)
    rr = np.concatenate([sl.t[::-1], sh.t[1:]])
    RR = np.concatenate([sl.y[0][::-1], sh.y[0][1:]])
    spl = CubicSpline(rr, RR)
    _KR[n] = (spl, g, F_over_x)
    return _KR[n]


def kruskal_exact(n, U, V, a=1.0):
    """Exact r and f at Kruskal coordinates (U, V) (arrays), x+ = a (V - V0), x- = a (U - U0): UV = G(r) by bisection,
    e^{-f} = [F/((r-1) e^{2 kappa R})]/(2 kappa^2 a^2)."""
    spl, g, F_over_x = kruskal_tools(n)
    kap = (n - 1.0) / 2.0
    Z = np.asarray(U * V, float)
    lo = np.full_like(Z, 0.02)
    hi = np.full_like(Z, 6.0)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        Gm = -(mid - 1.0) * np.exp(2 * kap * spl(mid))
        up = Gm > Z                                   # G decreases with r: G(mid) > Z means r is larger
        lo = np.where(up, mid, lo)
        hi = np.where(up, hi, mid)
    r = 0.5 * (lo + hi)
    emf = F_over_x(r - 1.0) * np.exp(-2 * kap * spl(r)) / (2 * kap ** 2 * a ** 2)
    return r, -np.log(emf)


def tangherlini_rays(n, N, h, U0=-0.5, V0=1.0):
    """Exact Tangherlini (r_h = 1, M = 1/2) in Kruskal null coordinates on the rays through (U0, V0): the horizon is U = 0,
    x- = -U0, a regular null line; r = r_h on it. Vacuum, no scalar."""
    s = np.arange(N + 1) * h
    Sp = np.zeros((2, N + 1))
    Sm = np.zeros((2, N + 1))
    Sp[0], Sp[1] = kruskal_exact(n, U0 + 0 * s, V0 + s)
    Sm[0], Sm[1] = kruskal_exact(n, U0 + s, V0 + 0 * s)
    d = Data(n, [], N, h, Sp, Sm, name="tangherlini")
    d.M0, d.rh, d.U0, d.V0 = 0.5, 1.0, U0, V0
    return d


def tangherlini_exact_nodes(d):
    i = np.arange(d.N + 1)[:, None]
    j = np.arange(d.N + 1)[None, :]
    return kruskal_exact(d.n, d.U0 + j * d.h + 0 * i, d.V0 + i * d.h + 0 * j)


def roberts_exact(p, m, a):
    r = np.sqrt(a * (p ** 2 + m ** 2) - p * m)
    th = math.acos(1.0 / (2.0 * a))
    ph = SQ2 * np.arctan2(m * math.sin(th), p - m * math.cos(th))
    return r, ph


def roberts_rays(N, h, a=0.75, p0=1.0, m0=1.0):
    """n = 2 ghost, exact Roberts-type solution r^2 = a(p^2 + m^2) - p m, f = 0 (derive.py D5), on the rays through (p0, m0)."""
    s = np.arange(N + 1) * h
    Sp = np.zeros((3, N + 1))
    Sm = np.zeros((3, N + 1))
    r, ph = roberts_exact(p0 + s, m0, a)
    Sp[0], Sp[2] = r, ph
    r, ph = roberts_exact(p0, m0 + s, a)
    Sm[0], Sm[2] = r, ph
    d = Data(2, [-1.0], N, h, Sp, Sm, x0p=p0, x0m=m0, name="roberts")
    d.a = a
    return d


# --------------------------------------------------------------------------------------------------------------------
# evolution
# --------------------------------------------------------------------------------------------------------------------
class Result:
    pass


def strip_cells(old, new, h):
    """Cell-centre quantities for the strip between two rows (nf, N+1): r, f, r_+, r_-, and the full centred U."""
    S, E, W, Nn = old[:, :-1], old[:, 1:], new[:, :-1], new[:, 1:]
    c = 0.25 * (S + E + W + Nn)
    up = (E + Nn - S - W) / (2 * h)
    um = (W + Nn - S - E) / (2 * h)
    return c, up, um


def evolve(d, kernel="numba", store=True, stop=False, rcut=0.05, fcut=25.0, excise=True, rexp=2.0, ncorr=1, tol=1e-2):
    """Evolve the data d over its whole rectangle (or until the first event if stop). Returns a Result with

    G[q, i, j] (store), dead[i, j], the diagonal series (t, tau, r, f, H, rp, rm at the cell centres x+ = x- along the
    rectangle's diagonal), the events and the first rows (for E0)."""
    N, h, nf, eps = d.N, d.h, d.nf, d.eps
    par = kernel_params(d.n, eps, ncorr)
    res = Result()
    res.d, res.kernel, res.rcut, res.fcut, res.tol = d, kernel, rcut, fcut, tol
    res.events = {}
    res.rows_done = N
    G = np.empty((nf, N + 1, N + 1)) if store else None
    dead = np.zeros((N + 1, N + 1), np.uint8) if store else None
    if kernel == "numpy" or not HAVE_NUMBA:
        if kernel == "numba" and not HAVE_NUMBA:
            raise RuntimeError("numba is not available; use --kernel=numpy")
        G = np.empty((nf, N + 1, N + 1))
        dead = np.zeros((N + 1, N + 1), np.uint8)
        G[:, :, 0] = d.Sp
        G[:, 0, :] = d.Sm
        spert = 1.0 + 1e-8 if mut("numpy_Spert") else 1.0      # a transcription slip in the numpy kernel only
        for k in range(0, 2 * N - 1):
            lo, hi = max(1, k + 2 - N), min(N, k + 1)
            if hi < lo:
                continue
            I = np.arange(lo, hi + 1)
            J = k + 2 - I
            S = G[:, I - 1, J - 1]
            E = G[:, I, J - 1]
            W = G[:, I - 1, J]
            S = S * spert
            dS, dE, dW = dead[I - 1, J - 1], dead[I, J - 1], dead[I - 1, J]
            dd = (dS | dE | dW).astype(bool)
            with np.errstate(all="ignore"):
                out = cell_np(S, E, W, h, par, eps)
                bad = ~np.all(np.isfinite(out), axis=0)
                if excise:
                    dd = dd | bad | (out[0] < rcut) | (np.exp(-out[1]) > fcut)
            out = np.where(dd[None, :], E, out)
            G[:, I, J] = out
            dead[I, J] = dd
        res.G, res.dead = G, dead
        res.rows_done = N
        res.diag = diag_from_grid(res)
        res.events = events_from_grid(res, rexp, tol)
        res.first = G[:, :, :3].copy()
        return res
    # numba path: rows
    old = np.ascontiguousarray(d.Sp.copy())
    new = np.empty_like(old)
    dold = np.zeros(N + 1, np.uint8)
    dnew = np.zeros(N + 1, np.uint8)
    if store:
        G[:, :, 0] = old
    first = [old.copy()]
    diag = {k: np.full(N, np.nan) for k in ("r", "f", "rp", "rm")}
    diag["phi"] = np.full((nf - 2, N), np.nan)
    for j in range(N):
        new[:, 0] = d.Sm[:, j + 1]
        dnew[0] = 0
        row_step(old, new, dold, dnew, h, par, eps, rcut, fcut, int(excise))
        c, up, um = strip_cells(old, new, h)
        alive = ~(dold[:-1].astype(bool) | dold[1:].astype(bool) | dnew[:-1].astype(bool) | dnew[1:].astype(bool))
        if alive[j]:
            diag["r"][j], diag["f"][j] = c[0, j], c[1, j]
            diag["rp"][j], diag["rm"][j] = up[0, j], um[0, j]
            for kk in range(nf - 2):
                diag["phi"][kk, j] = c[2 + kk, j]
        if store:
            G[:, :, j + 1] = new
            dead[:, j + 1] = dnew
        if j < 2:
            first.append(new.copy())
        if stop and "collapse" not in res.events and "expand" not in res.events:
            cc = cell_class(up[0], um[0], tol)
            trapped = alive & (cc == -1)
            if trapped.any():
                res.events["collapse"] = (j, int(np.argmax(trapped)))
            elif alive[j] and cc[j] == 1 and c[0, j] >= rexp:
                res.events["expand"] = (j, j)
            if res.events:
                res.rows_done = j + 1
                break
        old, new = new, old
        dold, dnew = dnew, dold
    res.G, res.dead = G, dead
    res.first = np.stack(first[:3] + [first[-1]] * max(0, 3 - len(first)), axis=2)
    res.diag = finish_diag(res, diag)
    if store and not stop:
        res.events = events_from_grid(res, rexp, tol)
    return res


def finish_diag(res, diag):
    d = res.d
    h = d.h
    m = res.rows_done
    j = np.arange(m)
    t = SQ2 * (j + 0.5) * h + (d.x0p + d.x0m) / SQ2
    r, f, rp, rm = (diag[k][:m] for k in ("r", "f", "rp", "rm"))
    ef2 = np.exp(f / 2)
    with np.errstate(all="ignore"):
        H = ef2 * (rp + rm) / (SQ2 * r)
    dtau = SQ2 * h * np.exp(-f / 2)
    tau = np.nancumsum(dtau) - 0.5 * np.nan_to_num(dtau)
    return dict(t=t, tau=tau, r=r, f=f, rp=rp, rm=rm, H=H, phi=diag["phi"][:, :m])


def diag_from_grid(res):
    d, G = res.d, res.G
    N, h, nf = d.N, d.h, d.nf
    c, up, um, alive = cell_fields(res)
    diag = {k: np.full(N, np.nan) for k in ("r", "f", "rp", "rm")}
    j = np.arange(N)
    ok = alive[j, j]
    diag["r"][ok], diag["f"][ok] = c[0][j, j][ok], c[1][j, j][ok]
    diag["rp"][ok], diag["rm"][ok] = up[0][j, j][ok], um[0][j, j][ok]
    diag["phi"] = np.full((nf - 2, N), np.nan)
    for kk in range(nf - 2):
        diag["phi"][kk, ok] = c[2 + kk][j, j][ok]
    return finish_diag(res, diag)


def cell_fields(res):
    """Cell-centre c, r_+, r_- as (nf, N, N)[q, i, j] and the alive mask [i, j]."""
    G, h = res.G, res.d.h
    S, E = G[:, :-1, :-1], G[:, 1:, :-1]
    W, Nn = G[:, :-1, 1:], G[:, 1:, 1:]
    c = 0.25 * (S + E + W + Nn)
    up = (E + Nn - S - W) / (2 * h)
    um = (W + Nn - S - E) / (2 * h)
    dd = res.dead.astype(bool)
    alive = ~(dd[:-1, :-1] | dd[1:, :-1] | dd[:-1, 1:] | dd[1:, 1:])
    return c, up, um, alive


def cell_class(rp, rm, tol):
    """+1 anti-trapped (r_+ > tol and r_- > tol), -1 trapped (both < -tol), 0 otherwise; the sign test of the events."""
    s = -1.0 if mut("finder_sign") else 1.0
    return np.where((s * rp > tol) & (s * rm > tol), 1, np.where((s * rp < -tol) & (s * rm < -tol), -1, 0))


def events_from_grid(res, rexp=2.0, tol=1e-2):
    """First trapped cell in row order (smallest j, then smallest i), and the first row whose diagonal cell is anti-trapped
    with r >= rexp; the earlier row wins, a tie goes to the collapse."""
    c, up, um, alive = cell_fields(res)
    cc = cell_class(up[0], um[0], tol)
    trapped = alive & (cc == -1)
    ev = {}
    rows = np.where(trapped.any(axis=0))[0]
    jc = None
    if rows.size:
        jc = int(rows[0])
        ev["collapse"] = (jc, int(np.argmax(trapped[:, jc])))
    N = res.d.N
    jj = np.arange(N)
    dg = alive[jj, jj] & (cc[jj, jj] == 1) & (c[0][jj, jj] >= rexp)
    if dg.any():
        je = int(np.argmax(dg))
        if jc is None or je < jc:
            ev = {"expand": (je, je)}
    return ev


def event_time(res, key):
    j, i = res.events[key]
    d = res.d
    xp, xm = d.x0p + (i + 0.5) * d.h, d.x0m + (j + 0.5) * d.h
    return xp, xm, (xp + xm) / SQ2


# --------------------------------------------------------------------------------------------------------------------
# post-hoc diagnostics on a stored grid
# --------------------------------------------------------------------------------------------------------------------
def node_derivs(res):
    """Second-order r_+, r_-, f_+, f_-, phi_+, phi_- at nodes (np.gradient, one-sided at the boundary, edge_order 2)."""
    G, h = res.G, res.d.h
    gp = np.gradient(G, h, axis=1, edge_order=2)
    gm = np.gradient(G, h, axis=2, edge_order=2)
    return gp, gm


def mass_nodes(res):
    d = res.d
    gp, gm = node_derivs(res)
    r, f = res.G[0], res.G[1]
    ex = d.n if mut("mass_exponent") else d.n - 1
    return r ** ex * (1 + 2 * np.exp(f) * gp[0] * gm[0]) / 2


def constraints(res):
    """C_++ at nodes 1..N-1 (i) and C_-- at nodes 1..N-1 (j), centred differences; arrays [i, j] with NaN where dead."""
    d, G, h = res.d, res.G, res.d.h
    n, eps = d.n, d.eps
    S = mut("constraint_sign")
    sg = -1.0 if S else 1.0
    r, f = G[0], G[1]
    rpp = (r[2:, :] - 2 * r[1:-1, :] + r[:-2, :]) / h ** 2
    fp = (f[2:, :] - f[:-2, :]) / (2 * h)
    rp = (r[2:, :] - r[:-2, :]) / (2 * h)
    sp = sum(eps[k] * ((G[2 + k][2:, :] - G[2 + k][:-2, :]) / (2 * h)) ** 2 for k in range(len(eps))) if len(eps) else 0.0
    Cpp = rpp + fp * rp + sg * (r[1:-1, :] / n) * sp
    rmm = (r[:, 2:] - 2 * r[:, 1:-1] + r[:, :-2]) / h ** 2
    fm = (f[:, 2:] - f[:, :-2]) / (2 * h)
    rm = (r[:, 2:] - r[:, :-2]) / (2 * h)
    sm = sum(eps[k] * ((G[2 + k][:, 2:] - G[2 + k][:, :-2]) / (2 * h)) ** 2 for k in range(len(eps))) if len(eps) else 0.0
    Cmm = rmm + fm * rm + sg * (r[:, 1:-1] / n) * sm
    dd = res.dead.astype(bool)
    Cpp = np.where(dd[1:-1, :] | dd[2:, :] | dd[:-2, :], np.nan, Cpp)
    Cmm = np.where(dd[:, 1:-1] | dd[:, 2:] | dd[:, :-2], np.nan, Cmm)
    return Cpp, Cmm


def edge_flux_residual(res):
    """On the edge x+ = X: M(N, j) - M(N, 0) against the integral of -(r^n/n) e^f r_+ sum eps phi_-^2 d x-; returns the max
    abs residual over j and the largest |M change|."""
    d, G, h = res.d, res.G, res.d.h
    n, eps = d.n, d.eps
    gp, gm = node_derivs(res)
    M = mass_nodes(res)[-1, :]
    r, f = G[0, -1, :], G[1, -1, :]
    sm = sum(eps[k] * gm[2 + k][-1, :] ** 2 for k in range(len(eps)))
    sgn = -1.0 if mut("flux_sign") else 1.0
    flux = -sgn * (r ** n / n) * np.exp(f) * gp[0][-1, :] * sm
    I = np.concatenate([[0.0], np.cumsum(0.5 * h * (flux[1:] + flux[:-1]))])
    resid = (M - M[0]) - I
    return float(np.max(np.abs(resid))), float(np.max(np.abs(M - M[0])))


def contour(res, which, axis="i"):
    """Zero contour of r_+ (which = '+') or r_- ('-') by linear interpolation of the cell-centre values. axis 'i' scans along
    x+ in each row j (a crossing per row), axis 'j' scans along x- in each column i. Returns arrays (x+, x-) of the crossings,
    the cell-centre r interpolated to them, and the best-fit slope dx+/dx- (axis 'i') or dx-/dx+ (axis 'j')."""
    d = res.d
    c, up, um, alive = cell_fields(res)
    if mut("finder_swap"):
        which = "-" if which == "+" else "+"
    v = up[0] if which == "+" else um[0]
    rc = c[0]
    xs, ms, rs = [], [], []
    if axis == "i":
        for j in range(d.N):
            col, ok = v[:, j], alive[:, j]
            idx = np.where(ok[:-1] & ok[1:] & (np.sign(col[:-1]) * np.sign(col[1:]) < 0))[0]
            for i in idx:
                fr = col[i] / (col[i] - col[i + 1])
                xs.append(d.x0p + (i + 0.5 + fr) * d.h)
                ms.append(d.x0m + (j + 0.5) * d.h)
                rs.append(rc[i, j] + fr * (rc[i + 1, j] - rc[i, j]))
    else:
        for i in range(d.N):
            col, ok = v[i, :], alive[i, :]
            idx = np.where(ok[:-1] & ok[1:] & (np.sign(col[:-1]) * np.sign(col[1:]) < 0))[0]
            for j in idx:
                fr = col[j] / (col[j] - col[j + 1])
                xs.append(d.x0p + (i + 0.5) * d.h)
                ms.append(d.x0m + (j + 0.5 + fr) * d.h)
                rs.append(rc[i, j] + fr * (rc[i, j + 1] - rc[i, j]))
    xs, ms, rs = np.array(xs), np.array(ms), np.array(rs)
    if len(xs) > 2:
        slope = np.polyfit(ms, xs, 1)[0] if axis == "i" else np.polyfit(xs, ms, 1)[0]
    else:
        slope = np.nan
    return xs, ms, rs, slope


def throat_min(res):
    """On each slice t = const of the ray gauge (anti-diagonal i + j = s of the nodes), the minimum of r with a parabolic
    refinement: returns (t, r_min, l_min) for the slices that are alive near their minimum."""
    d, G, h = res.d, res.G, res.d.h
    R = G[0]
    N = d.N
    out = []
    for s in range(2, 2 * N - 1):
        lo, hi = max(0, s - N), min(N, s)
        i = np.arange(lo, hi + 1)
        v = R[i, s - i]
        if len(v) < 3 or np.any(res.dead[i, s - i]):
            continue
        k = int(np.argmin(v))
        if k == 0 or k == len(v) - 1:
            continue
        y0, y1, y2 = v[k - 1], v[k], v[k + 1]
        den = y0 - 2 * y1 + y2
        dk = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        rmin = y1 - 0.25 * (y0 - y2) * dk
        ip = i[k] + dk
        l = (ip - (s - ip)) * h / SQ2
        t = (s * h + d.x0p + d.x0m) / SQ2
        out.append((t, rmin, l))
    return np.array(out)


def stress(res):
    """T_++, T_--, angular T^theta_theta, rho = T_uu at cell centres [i, j]."""
    d = res.d
    c, up, um, alive = cell_fields(res)
    k = len(d.eps)
    Tpp = sum(d.eps[a] * up[2 + a] ** 2 for a in range(k))
    Tmm = sum(d.eps[a] * um[2 + a] ** 2 for a in range(k))
    ef = np.exp(c[1])
    Tang = ef * sum(d.eps[a] * up[2 + a] * um[2 + a] for a in range(k))
    rho = ef * (Tpp + Tmm) / 2
    return Tpp, Tmm, Tang, rho


def ellis_exact_nodes(d):
    """Exact static Ellis fields at the nodes of d's rectangle: r(l), f = 0, phi(l), l = (x+ - x-)/sqrt2 (corner at the throat)."""
    sol, q = ellis_dense(d.n, d.N * d.h / SQ2 + 1.0)
    i = np.arange(d.N + 1)[:, None]
    j = np.arange(d.N + 1)[None, :]
    l = (i - j) * d.h / SQ2
    y = sol(np.abs(l).ravel())
    r = y[0].reshape(l.shape)
    ph = (np.sign(l) * y[2].reshape(l.shape))
    return r, ph

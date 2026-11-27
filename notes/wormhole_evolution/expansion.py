#!/usr/bin/env python3
"""The Ellis threshold's expanding branch, followed past the linear mode (#213). Uses dn.py and checks.py of #39 unchanged. The observer, the
definitions and the criteria are CRIT and MUTATIONS below, as declared in negative-mass-raw I213-dh610/README.md before the runs of record.

    python3 expansion.py job --kind=S --D=4 --A=0.05 --h=0.02 --X=200 --out=DIR      one run of the symmetric family S, streamed
    python3 expansion.py job --kind=Sp --D=4 --h=0.02 --X=200 --out=DIR              S at p*_h + delta (bisected first), streamed
    python3 expansion.py job --kind=frw --out=DIR                                    the planted homogeneous ghost (exact open FRW), three h
    python3 expansion.py job --kind=ctrl --D=4 --h=0.02 --out=DIR                    positive-energy controls (they collapse)
    python3 expansion.py analyze --raw=DIR                                           every check and table, from the stored runs (exit 1 if a check fails)
    python3 expansion.py mutate --mutate=NAME                                        the suite at reduced size with one input broken; the checks in MUTATIONS must fire
    python3 expansion.py dev                                                         the reduced suite, nothing broken

Units a = 1, 8 pi G = 1, t = (x+ + x-)/sqrt2, tau the proper time along the diagonal x+ = x-, the throat of the mirror-symmetric run.
The streamed evolution is dn.row_step unchanged (the row sweep of dn.evolve); only the diagnostics differ, and `stream` reproduces dn.evolve's
diagonal and edge masses on a stored run (check ST).
"""
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import checks   # noqa: E402
import dn       # noqa: E402

SQ2 = math.sqrt(2.0)
W_EXACT = checks.W_EXACT
EMUT = None         # this file's mutation (dn.MUT holds dn.py's)


def emut(name):
    return EMUT == name


# --------------------------------------------------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------------------------------------------------
PULSE = dict(xc=3.0, w=0.5)      # #40's S family: phi-Gaussian, centre 3, width 0.5 (checks.fate_run)


def s_data(n, h, X, A):
    """Family S: ghost pulse of amplitude A on the '+' ray mirrored onto the '-' ray with the parity of the static field (odd).
    mirror_parity: the wrong (even) parity, so the data are not mirror symmetric."""
    N = int(round(X / h))
    p = dict(field=0, ray="+", A=A, xc=PULSE["xc"], w=PULSE["w"])
    if emut("mirror_parity"):
        pulses = [p, dict(p, ray="-")]
        return dn.ellis_rays(n, N, h, pulses=pulses, symmetric=False)
    return dn.ellis_rays(n, N, h, pulses=[p], symmetric=True)


def normal_data(n, h, X, A):
    """The positive-energy control: the same Gaussian in a normal scalar (eps = +1), mirrored (even), on the ghost wormhole."""
    N = int(round(X / h))
    field = 0 if emut("control_ghost") else 1
    eps = (-1.0,) if field == 0 else (-1.0, 1.0)
    return dn.ellis_rays(n, N, h, pulses=[dict(field=field, ray="+", A=A, xc=PULSE["xc"], w=PULSE["w"])], symmetric=True, eps=eps)


def frw_exact(eta):
    """Open FRW (k = -1) with a homogeneous ghost, D = 4: a^2 = cosh 2 eta, phi = sqrt6 atan tanh eta (conformal time eta; dev_06: checked
    symbolically, residuals 1e-30). Returns a, phi, H = a'/a^2, rho = p = -3/cosh^3 2 eta, and the proper time offset function."""
    a2 = np.cosh(2 * eta)
    return np.sqrt(a2), SQ2 * math.sqrt(3.0) * np.arctan(np.tanh(eta)), np.sinh(2 * eta) / a2 ** 1.5, -3.0 / a2 ** 3


def frw_data(h, X, eta0, chi0):
    """Characteristic data of the exact solution on the two rays through the corner u0 = eta0 + chi0, v0 = eta0 - chi0 (x+ = u, x- = v)."""
    N = int(round(X / h))
    u0, v0 = eta0 + chi0, eta0 - chi0
    s = np.arange(N + 1) * h
    Sp = np.zeros((3, N + 1))
    Sm = np.zeros((3, N + 1))
    for S, (u, v) in ((Sp, (u0 + s, v0 + 0 * s)), (Sm, (u0 + 0 * s, v0 + s))):
        eta, chi = (u + v) / 2, (u - v) / 2
        a, ph, _H, _rho = frw_exact(eta)
        S[0], S[1], S[2] = a * np.sinh(chi), -np.log(a * a / 2), ph
    return dn.Data(2, [-1.0], N, h, Sp, Sm, x0p=u0, x0m=v0, name="frw")


# --------------------------------------------------------------------------------------------------------------------
# the streamed evolution
# --------------------------------------------------------------------------------------------------------------------
def _ms(c, up, um, n):
    """Misner-Sharp mass at cell centres (the code's normalization)."""
    ex = n if dn.mut("mass_exponent") else n - 1
    return c[0] ** ex * (1 + 2 * np.exp(c[1]) * up[0] * um[0]) / 2


def stream(d, tol=1e-2, rcut=0.05, fcut=25.0, excise=True, colspec=None, log=None):
    """dn.evolve's numba row sweep with streaming diagnostics instead of the stored grid. Returns a dict of arrays:
    diagonal cell series, the minimum of r on each slice t = const (node level), the anti-trapped cells and the two horizon contours per row,
    the ray identity from the diagonal to the last cell, the constraints per row, the columns i0 in colspec (default N/4, N/2, N) with their
    neighbours, and the last three rows."""
    N, h, nf, eps, n = d.N, d.h, d.nf, d.eps, d.n
    if not dn.HAVE_NUMBA:
        raise RuntimeError("numba is not available")
    par = dn.kernel_params(n, eps, 1)
    colspec = list(colspec or (N // 4, N // 2, N))
    kk = nf - 2
    R = {k: np.full(N, np.nan) for k in ("r", "f", "rp", "rm", "rll", "asym_r", "asym_f", "asym_p", "emf_max", "emf_min",
                                         "xL", "rL", "xR", "rR", "m_diag", "m_edge", "I_tot", "I_pos", "I_neg",
                                         "Cpp", "Cpp_rel", "Cmm", "Cmm_rel", "x_anti_lo", "x_anti_hi")}
    R["phi"] = np.full((kk, N), np.nan)
    R["phip"] = np.full((kk, N), np.nan)
    R["phim"] = np.full((kk, N), np.nan)
    R["n_anti"] = np.zeros(N, np.int64)
    R["n_trap"] = np.zeros(N, np.int64)
    R["n_dead"] = np.zeros(N, np.int64)
    R["alive_diag"] = np.zeros(N, bool)
    R["cls_diag"] = np.zeros(N, np.int8)
    minr = np.full(2 * N + 1, np.inf)
    argi = np.zeros(2 * N + 1, np.int64)
    cols = {}
    for i0 in colspec:
        lo, hi = max(0, i0 - 2), min(N, i0 + 2)
        cols[i0] = (lo, hi, np.empty((nf, hi - lo + 1, N + 1)))
    last = {}
    ii = np.arange(N + 1)
    first_dead = None

    old = np.ascontiguousarray(d.Sp.copy())
    new = np.empty_like(old)
    dold = np.zeros(N + 1, np.uint8)
    dnew = np.zeros(N + 1, np.uint8)
    prev = None
    dprev = None

    def node_row(j, row):
        """per-row bookkeeping of the nodes of row j (r at s = i + j, edge columns)"""
        seg = minr[j:j + N + 1]
        msk = row[0] < seg
        seg[msk] = row[0][msk]
        argi[j:j + N + 1][msk] = ii[msk]
        for i0, (lo, hi, arr) in cols.items():
            arr[:, :, j] = row[:, lo:hi + 1]
        if j >= N - 2:
            last[j] = row.copy()

    def c_pp(row, drow):
        r, f = row[0], row[1]
        rpp = (r[2:] - 2 * r[1:-1] + r[:-2]) / h ** 2
        fp = (f[2:] - f[:-2]) / (2 * h)
        rp = (r[2:] - r[:-2]) / (2 * h)
        sp = sum(eps[k] * ((row[2 + k][2:] - row[2 + k][:-2]) / (2 * h)) ** 2 for k in range(kk)) if kk else 0.0
        sgc = -1.0 if dn.mut("constraint_sign") else 1.0
        t3 = sgc * (r[1:-1] / n) * sp
        C = rpp + fp * rp + t3
        sc = np.abs(rpp) + np.abs(fp * rp) + np.abs(t3)
        ok = ~(drow[1:-1].astype(bool) | drow[2:].astype(bool) | drow[:-2].astype(bool))
        return C, sc, ok

    node_row(0, old)
    for j in range(N):
        new[:, 0] = d.Sm[:, j + 1]
        dnew[0] = 0
        dn.row_step(old, new, dold, dnew, h, par, eps, rcut, fcut, int(excise))
        c, up, um = dn.strip_cells(old, new, h)
        alive = ~(dold[:-1].astype(bool) | dold[1:].astype(bool) | dnew[:-1].astype(bool) | dnew[1:].astype(bool))
        nd = int(dnew.sum())
        R["n_dead"][j] = nd
        if nd and first_dead is None:
            first_dead = (j + 1, int(np.argmax(dnew)))
        node_row(j + 1, new)
        # --- diagonal cell
        if alive[j]:
            R["alive_diag"][j] = True
            R["cls_diag"][j] = dn.cell_class(up[0, j], um[0, j], tol)
            R["r"][j], R["f"][j], R["rp"][j], R["rm"][j] = c[0, j], c[1, j], up[0, j], um[0, j]
            for k in range(kk):
                R["phi"][k, j], R["phip"][k, j], R["phim"][k, j] = c[2 + k, j], up[2 + k, j], um[2 + k, j]
        # r_ll across the diagonal node (j, j): neighbours (j+1, j-1) in the previous row, (j-1, j+1) in the new row
        if prev is not None and j >= 1 and j + 1 <= N:
            a_, b_, c_ = prev[:, j + 1], new[:, j - 1], old[:, j]
            R["rll"][j] = (a_[0] + b_[0] - 2 * c_[0]) / (2 * h * h)
            R["asym_r"][j] = abs(a_[0] - b_[0]) / c_[0]
            R["asym_f"][j] = abs(a_[1] - b_[1])
            R["asym_p"][j] = max(abs(a_[2 + k] + b_[2 + k]) for k in range(kk)) if kk else 0.0
        # --- gauge
        if alive.any():
            emf = np.exp(-c[1][alive])
            R["emf_max"][j], R["emf_min"][j] = emf.max(), emf.min()
        # --- anti-trapped cells, trapped cells, the two horizon contours (cell-centre r_+, r_-)
        cc = dn.cell_class(up[0], um[0], tol)
        anti = alive & (cc == 1)
        trap = alive & (cc == -1)
        R["n_anti"][j], R["n_trap"][j] = int(anti.sum()), int(trap.sum())
        if anti.any():
            ia = np.where(anti)[0]
            R["x_anti_lo"][j], R["x_anti_hi"][j] = d.x0p + (ia[0] + 0.5) * h, d.x0p + (ia[-1] + 0.5) * h
        sp_, sm_ = up[0], um[0]
        zl = np.where(alive[:-1] & alive[1:] & (sp_[:-1] < 0) & (sp_[1:] >= 0))[0]
        if len(zl):
            i = zl[0]
            fr = sp_[i] / (sp_[i] - sp_[i + 1])
            R["xL"][j] = d.x0p + (i + 0.5 + fr) * h
            R["rL"][j] = c[0, i] + fr * (c[0, i + 1] - c[0, i])
        zr = np.where(alive[:-1] & alive[1:] & (sm_[:-1] > 0) & (sm_[1:] <= 0))[0]
        if len(zr):
            i = zr[0]
            fr = sm_[i] / (sm_[i] - sm_[i + 1])
            R["xR"][j] = d.x0p + (i + 0.5 + fr) * h
            R["rR"][j] = c[0, i] + fr * (c[0, i + 1] - c[0, i])
        # --- the ray identity from the diagonal cell to the last cell of this row
        if j <= N - 2 and alive[j] and alive[N - 1]:
            mc = _ms(c, up, um, n)
            Tpp = sum(eps[k] * up[2 + k] ** 2 for k in range(kk))
            sg = -1.0 if dn.mut("flux_sign") else 1.0
            ig = -sg * (c[0] ** n / n) * np.exp(c[1]) * um[0] * Tpp
            seg = ig[j:]
            R["I_tot"][j] = h * seg.sum() - 0.5 * h * (seg[0] + seg[-1])
            R["I_pos"][j] = h * seg[seg > 0].sum()
            R["I_neg"][j] = h * seg[seg < 0].sum()
            R["m_diag"][j], R["m_edge"][j] = mc[j], mc[N - 1]
        # --- constraints: C++ on the new row, C-- on the old row (needs prev and new)
        C, sc, ok = c_pp(new, dnew)
        if ok.any():
            R["Cpp"][j] = np.abs(C[ok]).max()
            R["Cpp_rel"][j] = (np.abs(C[ok]) / np.maximum(sc[ok], 1e-300)).max()
        if prev is not None:
            r_, f_ = old[0], old[1]
            rmm = (new[0] - 2 * r_ + prev[0]) / h ** 2
            fm = (new[1] - prev[1]) / (2 * h)
            rm = (new[0] - prev[0]) / (2 * h)
            smm = sum(eps[k] * ((new[2 + k] - prev[2 + k]) / (2 * h)) ** 2 for k in range(kk)) if kk else 0.0
            t3 = (-1.0 if dn.mut("constraint_sign") else 1.0) * (r_ / n) * smm
            Cm = rmm + fm * rm + t3
            scm = np.abs(rmm) + np.abs(fm * rm) + np.abs(t3)
            okm = ~(dold.astype(bool) | dnew.astype(bool) | dprev.astype(bool))
            if okm.any():
                R["Cmm"][j] = np.abs(Cm[okm]).max()
                R["Cmm_rel"][j] = (np.abs(Cm[okm]) / np.maximum(scm[okm], 1e-300)).max()
        prev, dprev = old.copy(), dold.copy()
        old, new = new, old
        dold, dnew = dnew, dold
        if log and (j + 1) % max(1, N // 10) == 0:
            log(f"stream row {j + 1} of {N}")
    R["minr"], R["argi"] = minr, argi
    R["cols"] = {i0: (lo, arr) for i0, (lo, hi, arr) in cols.items()}
    R["last"] = last
    R["first_dead"] = first_dead
    R["N"], R["h"], R["n"] = N, h, n
    R["eps"] = np.array(eps, float)
    R["x0p"], R["x0m"] = d.x0p, d.x0m
    R.update(finish(R))
    return R


def finish(R):
    """t, tau, H on the diagonal as dn.finish_diag does."""
    N, h = R["N"], R["h"]
    j = np.arange(N)
    t = SQ2 * (j + 0.5) * h + (R["x0p"] + R["x0m"]) / SQ2
    r, f, rp, rm = R["r"], R["f"], R["rp"], R["rm"]
    with np.errstate(all="ignore"):
        H = np.exp(f / 2) * (rp + rm) / (SQ2 * r)
        dtau = SQ2 * h * np.exp(-f / 2)
    tau = np.nancumsum(dtau) - 0.5 * np.nan_to_num(dtau)
    if emut("time_unit"):
        tau = tau * SQ2
        H = H / SQ2
    return dict(t=t, tau=tau, H=H)


# --------------------------------------------------------------------------------------------------------------------
# the throat observer: rho, p, w on the diagonal
# --------------------------------------------------------------------------------------------------------------------
def throat_stress(R):
    """On the diagonal of the run, in the unit frame u = e^{f/2}(d_+ + d_-)/sqrt2, e_l = e^{f/2}(d_+ - d_-)/sqrt2, angular directions:
    rho = T_uu = e^f (T_++ + T_--)/2 + T_+- stuff, p_par = T_ll, p_perp = T^theta_theta, w = p_bar/rho with p_bar = (p_par + n p_perp)/(n+1).
    T_ab = sum eps (d_a phi d_b phi - g_ab (d phi)^2 / 2); g_{+-} = -e^{-f}. Returns rho, p_par, p_perp, w_bar, T_+- (zero analytically)."""
    n, eps = R["n"], R["eps"]
    f = R["f"]
    kk = len(eps)
    Tpp = sum(eps[k] * R["phip"][k] ** 2 for k in range(kk))
    Tmm = sum(eps[k] * R["phim"][k] ** 2 for k in range(kk))
    ef = np.exp(f)
    # T_+- = sum eps (phi_+ phi_- - g_{+-} (g^{cd} phi_c phi_d)/2) with g^{+-} = -e^f, (d phi)^2 = 2 g^{+-} phi_+ phi_-
    Tpm = sum(eps[k] * (R["phip"][k] * R["phim"][k] - 0.5 * (-np.exp(-f)) * (2 * (-ef) * R["phip"][k] * R["phim"][k])) for k in range(kk))
    rho = ef * (Tpp + Tmm + 2 * Tpm) / 2
    p_par = ef * (Tpp + Tmm - 2 * Tpm) / 2
    sg = -1.0 if emut("w_sign") else 1.0
    p_perp = sg * ef * sum(eps[k] * R["phip"][k] * R["phim"][k] for k in range(kk))   # T^theta_theta = -eps (d phi)^2 / 2
    with np.errstate(all="ignore"):
        w = (p_par + n * p_perp) / ((n + 1) * rho)
    return rho, p_par, p_perp, w, Tpm


def m_diag_cells(R):
    n = R["n"]
    ex = n if dn.mut("mass_exponent") else n - 1
    return R["r"] ** ex * (1 + 2 * np.exp(R["f"]) * R["rp"] * R["rm"]) / 2


# --------------------------------------------------------------------------------------------------------------------
# masses at the ends and the flux through them
# --------------------------------------------------------------------------------------------------------------------
def end_mass(R, i0):
    """Misner-Sharp mass at the nodes (i0, j), j = 0..N, the flux integrand -(r^n/n) e^f r_+ sum eps phi_-^2 and the residual of
    M(j) - M(0) - int flux dx-, as dn.mass_nodes and dn.edge_flux_residual (one-sided d_+ at the edge i0 = N)."""
    n, h, N = R["n"], R["h"], R["N"]
    eps = R["eps"]
    lo, arr = R["cols"][i0]
    k = i0 - lo
    if i0 == N:
        sl = lambda q: (3 * arr[q, k] - 4 * arr[q, k - 1] + arr[q, k - 2]) / (2 * h)
    else:
        sl = lambda q: (arr[q, k + 1] - arr[q, k - 1]) / (2 * h)
    r, f = arr[0, k], arr[1, k]
    gp = sl(0)
    gm = np.gradient(arr[:, k, :], h, axis=1, edge_order=2)
    ex = n if dn.mut("mass_exponent") else n - 1
    M = r ** ex * (1 + 2 * np.exp(f) * gp * gm[0]) / 2
    sm = sum(eps[q] * gm[2 + q] ** 2 for q in range(len(eps)))
    sg = -1.0 if dn.mut("flux_sign") else 1.0
    flux = -sg * (r ** n / n) * np.exp(f) * gp * sm
    I = np.concatenate([[0.0], np.cumsum(0.5 * h * (flux[1:] + flux[:-1]))])
    return M, flux, I, (M - M[0]) - I


def top_mass(R):
    """The same on the edge x- = X (row j = N): M(i, N), the flux -(r^n/n) e^f r_- sum eps phi_+^2, the residual."""
    n, h, N = R["n"], R["h"], R["N"]
    eps = R["eps"]
    rows = R["last"]
    a, b, c = rows[N - 2], rows[N - 1], rows[N]
    gm0 = (3 * c[0] - 4 * b[0] + a[0]) / (2 * h)
    gp = np.gradient(c, h, axis=1, edge_order=2)
    r, f = c[0], c[1]
    ex = n if dn.mut("mass_exponent") else n - 1
    M = r ** ex * (1 + 2 * np.exp(f) * gp[0] * gm0) / 2
    sp = sum(eps[q] * gp[2 + q] ** 2 for q in range(len(eps)))
    sg = -1.0 if dn.mut("flux_sign") else 1.0
    flux = -sg * (r ** n / n) * np.exp(f) * gm0 * sp
    I = np.concatenate([[0.0], np.cumsum(0.5 * h * (flux[1:] + flux[:-1]))])
    return M, flux, I, (M - M[0]) - I


# --------------------------------------------------------------------------------------------------------------------
# declared criteria (negative-mass-raw I213-dh610/README.md, before the runs of record) and mutations
# --------------------------------------------------------------------------------------------------------------------
CRIT = dict(
    trap_tol=1e-2,                                   # events and horizon cells as #39 and #40
    # runs: D (n = D - 2) -> rectangle side X, the two h (coarse, fine), the members (ghost amplitude A of S, or 'tuned' = p*_h + delta)
    X={2: 200.0, 3: 100.0},
    H={2: (0.02, 0.01), 3: (0.01 / SQ2, 0.005 / SQ2)},
    members={2: ("tuned", 0.01, 0.05, 0.2), 3: ("tuned", 0.05)},
    delta=1e-6, bis_X={2: 20.0, 3: 15.0}, bis_bracket=0.01, bis_tol=1e-10,         # p*_h by checks.bisect_A (#39's G1 settings)
    # observer
    OB_sym=1e-9,                                     # mirror asymmetry of the nodes next to the diagonal (relative r; absolute f, phi)
    OB_slice=1,                                      # the minimum of r on every full-length slice lies within this many nodes of the diagonal
    OB0=dict(n=2, h=0.02, X=6.0, t_c=4.0, tol=1e-3), # planted: the static wormhole for t <= t_c (#39's P1), r_diag = 1, minimum on the diagonal
    # w
    W1=dict(h=(0.04, 0.02, 0.01), X=2.5, eta0=-1.0, chi0=3.0, order=(1.75, 2.25), w_tol=1e-6),
    W2_tol=1e-9,                                     # identities on the symmetric diagonal: p_perp = -rho, p_par = rho, T_+- = 0
    # growth rate
    SL0=dict(window=(1e-3, 1e-2), bar=0.01),         # tuned member, |H|/W in the window: slope of ln|H| against tau equals W within 1 %
    SL=dict(r_eval=(10.0, 25.0, 50.0), bar=0.10),    # #211's test: |H/W - 1| <= bar at every r/a in r_eval that the run reaches
    # classification (final window: the last e-fold of r on the diagonal)
    CL=dict(nfold=1.0, Hratio=1.2, rhoratio=1.2, wband=0.2, ft_resid=1e-2, ft_tsmax=1.2, grow=1.2),
    RES=dict(H_rel=1e-2, k_rel=1e-2, ts_rel=1e-2, lnr=1e-3),
    # mass bookkeeping
    MB=dict(flux=1e-3, mirror=1e-6, ray=1e-3, ray_x=(3.0, 4.0, 6.0, 8.0, 12.0, 20.0, 40.0, 80.0)),
    # controls
    CT=dict(X=12.0, A_normal=0.05, A_ghost=-0.05, dt=4.0),
    ST=dict(rel=1e-9),
    order_lo=1.75, order_hi=2.25,
)

MUTATIONS = {   # name: the checks that must fire, exactly; computed on the frozen code in dev_13 (negative-mass-raw I213-dh610/dev/mut13), D = 4 reduced suite
    "mirror_parity": ["MB1", "MB2", "MB3", "OB1", "OB2", "RES", "W2"],          # the '-' pulse mirrored with the wrong parity: not mirror symmetric
    "w_sign": ["W1", "W2"],                                                      # p_perp with the wrong sign
    "time_unit": ["SL0", "ST", "W1"],                                            # tau scaled by sqrt2 (H by 1/sqrt2)
    "class_loose": ["PL"],                                                       # the classification tolerances widened
    "control_ghost": ["CT"],                                                     # the control pulse put in the ghost field instead of a normal one
    "theta_coef": ["BUILD", "CT", "MB1", "MB2", "OB0", "OB1", "OB2", "RES", "W1"],   # dn.py's kernel mutations; BUILD = the bisection of p*_h cannot straddle
    "scalar_coef": ["BUILD", "MB1", "MB2", "OB0", "OB1", "OB2", "RES", "W1"],
    "ghost_sign": ["BUILD", "MB1", "MB2", "OB0", "OB1", "OB2", "RES", "W1"],
    "no_corrector": ["MB1", "MB3", "RES", "W1"],
    "data_constraint_sign": ["BUILD", "CT", "MB1", "MB2", "OB0", "OB1", "OB2", "RES"],
    "ellis_charge": ["BUILD", "CT", "MB1", "MB3", "OB0"],
    "constraint_sign": [],                                                       # silent by design: only the reported constraint monitor reads it
    "flux_sign": ["MB1", "MB3"],
    "mass_exponent": ["MB1", "MB3"],
    "finder_sign": ["BUILD", "CT"],
}
DN_MUT = ("theta_coef", "scalar_coef", "ghost_sign", "no_corrector", "data_constraint_sign", "ellis_charge", "constraint_sign", "flux_sign",
          "mass_exponent", "finder_sign")
ALL_MUT = None

results = []      # (id, ok)
verdicts = []     # (id, text, value)


def check(cid, text, ok, detail=""):
    ok = None if ok is None else bool(ok)
    results.append((cid, ok))
    tag = "SKIP" if ok is None else ("pass" if ok else "FAIL")
    print(f"[{tag}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def verdict(cid, text, value):
    verdicts.append((cid, text, value))
    print(f"[{cid}] {text}: {value}", flush=True)


def fmt(xs, f="{:.3e}"):
    return "[" + ", ".join(f.format(x) for x in xs) + "]"


def orders(errs):
    return [math.log2(errs[i] / errs[i + 1]) if errs[i + 1] > 0 and errs[i] > 0 else float("inf") for i in range(len(errs) - 1)]


# --------------------------------------------------------------------------------------------------------------------
# classification of a throat series
# --------------------------------------------------------------------------------------------------------------------
def final_window(tau, r, nfold):
    """Indices of the last e-fold (nfold) of r on the finite part of the diagonal series."""
    ok = np.where(np.isfinite(r) & np.isfinite(tau))[0]
    rmax = r[ok[-1]]
    return ok[r[ok] >= rmax * math.exp(-nfold)]


def classify(tau, r, H, rho, w, crit=None):
    """The declared classes, evaluated on the last e-fold of r. In order: DS (de Sitter-like: H, |rho| constant within the ratios, |w + 1| <=
    wband), STIFF (|w - 1| <= wband: a stiff ghost-dominated state), FT (H grows by the factor `grow`, 1/H linear in tau within ft_resid,
    divergence time tau_s in (tau_end, ft_tsmax tau_end]), GROW (H grows, not FT), OTHER."""
    c = dict(CRIT["CL"] if crit is None else crit)
    if emut("class_loose"):
        c.update(Hratio=1e3, rhoratio=1e3, wband=10.0)
    ix = final_window(tau, r, c["nfold"])
    tw, rw, Hw, rhow, ww = tau[ix], r[ix], H[ix], rho[ix], w[ix]
    out = dict(npts=len(ix), r0=float(rw[0]), r1=float(rw[-1]), Hratio=float(Hw.max() / Hw.min()), rhoratio=float(np.abs(rhow).max() / np.abs(rhow).min()),
               wmin=float(ww.min()), wmax=float(ww.max()), Hgrow=float(Hw[-1] / Hw[0]))
    b, a = np.polyfit(tw, 1 / Hw, 1)
    res = float(np.max(np.abs(1 / Hw - (a + b * tw))) / np.mean(1 / Hw))
    ts = -a / b if b < 0 else float("inf")
    out.update(k=float(-1 / b) if b < 0 else float("inf"), tau_s=float(ts), ft_resid=res, tau_end=float(tw[-1]))
    if out["Hratio"] <= c["Hratio"] and out["rhoratio"] <= c["rhoratio"] and max(abs(out["wmin"] + 1), abs(out["wmax"] + 1)) <= c["wband"]:
        cls = "DS"
    elif max(abs(out["wmin"] - 1), abs(out["wmax"] - 1)) <= c["wband"]:
        cls = "STIFF"
    elif out["Hgrow"] >= c["grow"] and np.all(np.diff(Hw) > 0):
        cls = "FT" if (b < 0 and res <= c["ft_resid"] and tw[-1] < ts <= c["ft_tsmax"] * tw[-1]) else "GROW"
    else:
        cls = "OTHER"
    out["class"] = cls
    return out


def slope_test(r, H, W, crit=None):
    """#211's test: |H/W - 1| <= bar at each r/a in r_eval that the run reaches (H at the first diagonal point with r >= r_eval)."""
    c = CRIT["SL"] if crit is None else crit
    vals = {}
    ok = np.isfinite(r) & np.isfinite(H)
    for r0 in c["r_eval"]:
        jj = np.where(ok & (r >= r0))[0]
        if len(jj):
            vals[r0] = float(H[jj[0]] / W)
    met = bool(vals) and all(abs(v - 1) <= c["bar"] for v in vals.values())
    return met, vals


def planted_series():
    """Synthetic throat series with known classes, through classify and slope_test: (name, expected class, expected slope verdict, series)."""
    out = []
    W = W_EXACT[2]
    tau = np.linspace(0.0, 12.0, 24001)
    r = np.exp(W * tau)
    H = np.full_like(tau, W)
    out.append(("de Sitter, H = W, rho = const, w = -1", "DS", True, (tau, r, H, np.full_like(tau, 1.5), np.full_like(tau, -1.0))))
    ts = 5.0
    tau = ts - np.logspace(math.log10(5.0), -4.5, 24001)
    r = 1.0 + (ts - tau) ** -2.0
    H = 2.0 / (ts - tau)
    out.append(("finite-time divergence, r = (tau_s - tau)^-2, rho = -4H^2, w = -1/3", "FT", False, (tau, r, H, -4 * H ** 2, np.full_like(tau, -1.0 / 3))))
    tau = np.linspace(0.5, 40.0, 24001)
    r = tau ** (1.0 / 3)
    H = 1.0 / (3 * tau)
    out.append(("flat stiff fluid, a ~ tau^(1/3), w = +1", "STIFF", False, (tau, r, H, 3 * H ** 2, np.full_like(tau, 1.0))))
    tau = np.linspace(0.0, 2.0, 24001)
    r = np.exp(tau ** 2 / 2)
    H = tau + 1e-12
    out.append(("Gaussian growth, r = exp(tau^2/2), H = tau, w = -1/3", "GROW", False, (tau, r, H, -4 * H ** 2 - 1e-3, np.full_like(tau, -1.0 / 3))))
    tau = np.linspace(0.5, 40.0, 24001)
    r = tau ** 0.5
    H = 0.5 / tau
    out.append(("decelerating radiation, a ~ tau^(1/2), w = +1/3", "OTHER", False, (tau, r, H, 3 * H ** 2, np.full_like(tau, 1.0 / 3))))
    return out


# --------------------------------------------------------------------------------------------------------------------
# jobs and stored runs
# --------------------------------------------------------------------------------------------------------------------
def save_run(path, R, meta):
    arrs = {}
    for k, v in R.items():
        if isinstance(v, np.ndarray):
            arrs[k] = v
        elif k == "cols":
            for i0, (lo, arr) in v.items():
                arrs[f"cols_{i0}"] = arr
                meta.setdefault("cols", {})[str(i0)] = int(lo)
        elif k == "last":
            for j, row in v.items():
                arrs[f"last_{j}"] = row
        elif k == "first_dead":
            meta["first_dead"] = v
        elif k in ("N", "h", "n", "x0p", "x0m"):
            meta[k] = v
    np.savez_compressed(path, meta=json.dumps(meta), **arrs)


def load_run(path):
    z = np.load(path, allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    R = {}
    cols, last = {}, {}
    for k in z.files:
        if k == "meta":
            continue
        if k.startswith("cols_"):
            i0 = int(k[5:])
            cols[i0] = (meta["cols"][str(i0)], z[k])
        elif k.startswith("last_"):
            last[int(k[5:])] = z[k]
        else:
            R[k] = z[k]
    R["cols"], R["last"] = cols, last
    for k in ("N", "h", "n", "x0p", "x0m", "first_dead"):
        if k in meta:
            R[k] = meta[k]
    R["meta"] = meta
    return R


def pstar(n, h):
    """The discrete threshold p*_h of family S by checks.bisect_A (#39's G1 settings; fates by checks.fate_run, X = bis_X)."""
    checks.CRIT["trap_tol"] = CRIT["trap_tol"]
    lo, hi, k = checks.bisect_A(n, h, CRIT["bis_X"][n], -CRIT["bis_bracket"], CRIT["bis_bracket"], CRIT["bis_tol"])
    return 0.5 * (lo + hi), k


def run_member(n, h, X, member, log=print):
    meta = dict(kind="S", n=n, h=h, X=X, member=str(member), mut=EMUT, dn_mut=dn.MUT)
    if member == "tuned":
        ps, k = pstar(n, h)
        A = ps + CRIT["delta"]
        meta.update(pstar=ps, bisect_runs=k, delta=CRIT["delta"])
    else:
        A = float(member)
    meta["A"] = A
    t0 = time.time()
    d = s_data(n, h, X, A)
    R = stream(d, tol=CRIT["trap_tol"])
    meta["secs"] = time.time() - t0
    return R, meta


def job(argv):
    a = dict(x.lstrip("-").split("=", 1) for x in argv)
    kind = a["kind"]
    out = a["out"]
    os.makedirs(out, exist_ok=True)
    if kind in ("S", "Sp"):
        n = int(a["D"]) - 2
        h = float(a["h"])
        X = float(a.get("X", CRIT["X"][n]))
        member = "tuned" if kind == "Sp" else a["A"]
        R, meta = run_member(n, h, X, member)
        name = f"D{n + 2}_{'tuned' if kind == 'Sp' else 'A' + a['A']}_h{h:g}"
        save_run(os.path.join(out, name + ".npz"), R, meta)
        print(f"{name}: {meta['secs']:.1f} s, N = {R['N']}, A = {meta['A']:.10e}")
    elif kind == "frw":
        c = CRIT["W1"]
        for h in c["h"]:
            R = stream(frw_data(h, c["X"], c["eta0"], c["chi0"]), excise=False)
            save_run(os.path.join(out, f"frw_h{h:g}.npz"), R, dict(kind="frw", h=h))
            print(f"frw h = {h}: N = {R['N']}")
    elif kind == "ctrl":
        n = int(a["D"]) - 2
        h = float(a["h"])
        c = CRIT["CT"]
        for name, d in (("normal", normal_data(n, h, c["X"], c["A_normal"])), ("ghost_neg", s_data(n, h, c["X"], c["A_ghost"]))):
            res = dn.evolve(d, kernel="numba", tol=CRIT["trap_tol"], stop=False, store=True)
            info = control_info(res)
            info.update(name=name, n=n, h=h, X=c["X"], A=(c["A_normal"] if name == "normal" else c["A_ghost"]))
            json.dump(info, open(os.path.join(out, f"ctrl_D{n + 2}_{name}_h{h:g}.json"), "w"), indent=1, default=float)
            print(f"ctrl D{n + 2} {name} h = {h}: {info}")
    else:
        sys.exit(f"unknown job kind {kind}")


def control_info(res):
    """Fate, event cell and time, the mass at the first trapped cell (r^(n-1)/2 on the horizon), E0 and the edge mass."""
    d = res.d
    ev = res.events
    key = next(iter(ev), None)
    out = dict(fate=key or "none", dead=int(res.dead.sum()))
    M = dn.mass_nodes(res)
    sol, _q = dn.ellis_dense(d.n, d.N * d.h / SQ2 + 1.0)
    out["E0"] = float(M[-1, 0] - 1.0 / (2.0 * sol(d.N * d.h / SQ2)[0] ** (d.n - 1)))
    if key:
        out["t_event"] = dn.event_time(res, key)[2]
        if key == "collapse":
            j0, i0 = ev["collapse"]
            c, up, um, alive = dn.cell_fields(res)
            out["r_first_trapped"] = float(c[0][i0, j0])
            out["M_first_trapped"] = float(c[0][i0, j0] ** (d.n - 1) / 2)
            out["M_cell"] = float(dn.mass_nodes(res)[i0, j0])
    out["M_edge_end"] = float(M[-1, -1])
    return out


# --------------------------------------------------------------------------------------------------------------------
# the checks
# --------------------------------------------------------------------------------------------------------------------
def diag_arrays(R):
    """The finite part of the diagonal series: tau, r, H, rho, p_par, p_perp, w, m, and the index array."""
    ix = np.where(R["alive_diag"] & np.isfinite(R["r"]) & np.isfinite(R["H"]))[0]
    rho, ppar, pperp, w, Tpm = throat_stress(R)
    m = m_diag_cells(R)
    return dict(ix=ix, t=R["t"][ix], tau=R["tau"][ix], r=R["r"][ix], H=R["H"][ix], rho=rho[ix], ppar=ppar[ix], pperp=pperp[ix], w=w[ix],
                m=m[ix], Tpm=Tpm[ix], f=R["f"][ix])


def chk_ST():
    rel = CRIT["ST"]["rel"]
    worst = dict(diag=0.0, M=0.0, resid=0.0, C=0.0, cls=0)
    for name, d in (("S", s_data(2, 0.04, 20.0, 0.05)), ("normal (collapse)", normal_data(2, 0.04, 20.0, 0.05))):
        res = dn.evolve(d, kernel="numba", excise=True, rexp=1e9, tol=CRIT["trap_tol"])
        R = stream(d, tol=CRIT["trap_tol"])
        dg = res.diag
        for k in ("r", "f", "rp", "rm", "H", "tau"):
            a, b = R[k], dg[k]
            worst["diag"] = max(worst["diag"], float(np.nanmax(np.abs(a - b)) / max(1.0, np.nanmax(np.abs(b)))))
        N = d.N
        Mdn = dn.mass_nodes(res)[N, :]
        Mst, _fl, _I, resid = end_mass(R, N)
        worst["M"] = max(worst["M"], float(np.max(np.abs(Mst - Mdn)) / max(1.0, np.max(np.abs(Mdn)))))
        rr, dM = dn.edge_flux_residual(res)
        worst["resid"] = max(worst["resid"], abs(float(np.max(np.abs(resid))) - rr) / max(dM, 1.0))
        Cpp, Cmm = dn.constraints(res)
        worst["C"] = max(worst["C"], abs(float(np.nanmax(R["Cpp"])) - float(np.nanmax(np.abs(Cpp)))) / max(1.0, float(np.nanmax(np.abs(Cpp)))))
        c, up, um, alive = dn.cell_fields(res)
        cc = dn.cell_class(up[0], um[0], CRIT["trap_tol"])
        worst["cls"] += int(np.sum(R["n_anti"] != ((cc == 1) & alive).sum(axis=0)) + np.sum(R["n_trap"] != ((cc == -1) & alive).sum(axis=0)))
    ok = worst["diag"] <= rel and worst["M"] <= rel and worst["resid"] <= rel and worst["C"] <= rel and worst["cls"] == 0
    check("ST", "the streamed run reproduces dn.evolve (diagonal, edge masses, flux residual, constraints, cell classes) on two grids, one with excision",
          ok, f"diagonal {worst['diag']:.1e}, edge mass {worst['M']:.1e}, flux residual {worst['resid']:.1e}, C++ {worst['C']:.1e}, class counts {worst['cls']}; tol {rel:g}")


def slice_offset(R, smax=None):
    """Largest |i - (s - i)| of the node-level minimum of r on the full-length slices s = 2 .. N."""
    N = R["N"]
    s = np.arange(2, (smax or N) + 1)
    return int(np.max(np.abs(2 * R["argi"][s] - s)))


def chk_OB0():
    c = CRIT["OB0"]
    N = int(round(c["X"] / c["h"]))
    R = stream(dn.ellis_rays(c["n"], N, c["h"]), tol=CRIT["trap_tol"])
    off = slice_offset(R)
    er = float(np.nanmax(np.abs(R["r"] - 1.0)[R["t"] <= c["t_c"]]))
    check("OB0", "planted: the static Ellis wormhole has r = 1 on the diagonal (t <= t_c) and its minimum of r on every full-length slice at the diagonal",
          off <= CRIT["OB_slice"] and er <= c["tol"], f"max |r_diag - 1| {er:.2e} (tol {c['tol']:g}), largest slice offset {off} nodes (tol {CRIT['OB_slice']}), t <= {c['t_c']}")


def chk_W1():
    c = CRIT["W1"]
    errs = {k: [] for k in ("r", "f", "phi", "H", "tau", "rho", "p_perp", "p_par")}
    wdev = []
    bounce = True
    from scipy.integrate import quad
    for h in c["h"]:
        R = stream(frw_data(h, c["X"], c["eta0"], c["chi0"]), excise=False)
        N = R["N"]
        eta = c["eta0"] + (np.arange(N) + 0.5) * h
        a, ph, H, rho = frw_exact(eta)
        tau_ex = np.array([quad(lambda e: math.sqrt(math.cosh(2 * e)), c["eta0"], e0, epsabs=1e-13, epsrel=1e-13)[0] for e0 in eta])
        rho_n, ppar, pperp, w, Tpm = throat_stress(R)
        for k, x, y in (("r", R["r"], a * math.sinh(c["chi0"])), ("f", R["f"], -np.log(a * a / 2)), ("phi", R["phi"][0], ph), ("H", R["H"], H),
                        ("tau", R["tau"], tau_ex), ("rho", rho_n, rho), ("p_perp", pperp, rho), ("p_par", ppar, rho)):
            errs[k].append(float(np.nanmax(np.abs(x - y))))
        wdev.append(float(np.nanmax(np.abs(w - 1.0))))
        bounce &= bool(H.min() < 0 < H.max())
    ords = {k: orders(v) for k, v in errs.items()}
    ok_o = all(c["order"][0] <= x <= c["order"][1] for v in ords.values() for x in v)
    for k in errs:
        print(f"    W1 {k:7s} errors {fmt(errs[k])} orders {fmt(ords[k], '{:.3f}')}")
    print(f"    W1 |w_bar - 1| at the three h {fmt(wdev)}; the observer passes the bounce (H changes sign): {bounce}")
    check("W1", "planted homogeneous ghost (exact open FRW, D = 4): (r, f, phi), H(tau), tau, rho, p_perp, p_par converge at second order and w = p_bar/rho = +1",
          ok_o and wdev[-1] <= c["w_tol"] and bounce, f"orders in [{c['order'][0]}, {c['order'][1]}], |w - 1| <= {c['w_tol']:g}")


def chk_PL():
    ok = True
    rows = []
    for name, cls_exp, sl_exp, (tau, r, H, rho, w) in planted_series():
        out = classify(tau, r, H, rho, w)
        met, vals = slope_test(r, H, W_EXACT[2])
        good = out["class"] == cls_exp and met == sl_exp
        ok &= good
        rows.append(f"{name}: class {out['class']} (expected {cls_exp}), slope test {'met' if met else 'not met'} (expected {'met' if sl_exp else 'not met'})")
    for x in rows:
        print("    PL " + x)
    check("PL", "planted series with known answers through classify and slope_test: de Sitter, finite-time divergence, stiff, Gaussian growth, decelerating", ok)


def chk_member(M):
    """Per-run checks OB1, OB2, W2, SL0 for one stored or mini run; M = dict(n, h, member, A, R)."""
    R, n, h = M["R"], M["n"], M["h"]
    lab = f"D{n + 2} h={h:.5f} {M['member']}"
    asym = max(float(np.nanmax(R["asym_r"])), float(np.nanmax(R["asym_f"])), float(np.nanmax(R["asym_p"])))
    off = slice_offset(R)
    rll = float(np.nanmin(R["rll"][np.isfinite(R["r"])]))
    M["asym"], M["off"], M["rll"] = asym, off, rll
    g = diag_arrays(R)
    ok = np.isfinite(g["rho"]) & (np.abs(g["rho"]) > 0)
    dpp = float(np.max(np.abs(g["pperp"][ok] / g["rho"][ok] + 1)))
    dw = float(np.max(np.abs(g["w"][ok] - (1 - n) / (1 + n))))
    dpl = float(np.max(np.abs(g["ppar"][ok] / g["rho"][ok] - 1)))
    dT = float(np.max(np.abs(g["Tpm"][ok] / g["rho"][ok])))
    M["w2"] = (dpp, dw, dpl, dT)
    return lab, asym, off, rll, dpp, dw, dpl, dT


def sl0(M):
    """Slope of ln|H| against tau after the global minimum of |H| with |H|/W in the window (checks.fit_growth)."""
    g = diag_arrays(M["R"])
    W = W_EXACT[M["n"]]
    H, tau = g["H"], g["tau"]
    i0 = int(np.argmin(np.abs(H)))
    lvl = np.abs(H) / W
    c = CRIT["SL0"]["window"]
    m = (np.arange(len(H)) >= i0) & (lvl >= c[0]) & (lvl <= c[1]) & (H != 0)
    if m.sum() < 10:
        return float("nan"), int(m.sum())
    return float(np.polyfit(tau[m], np.log(np.abs(H[m])), 1)[0]), int(m.sum())


def running_exponent(r, tau):
    """d ln(r/a - 1)/d tau, the exponent of [SH02]'s fit form r/a = 1 + b exp(H (tau - b5)), by centred differences on the series."""
    with np.errstate(all="ignore"):
        return np.gradient(np.log(r - 1.0), tau)


def series_table(M, lines):
    R, n = M["R"], M["n"]
    g = diag_arrays(R)
    W = W_EXACT[n]
    sp = running_exponent(g["r"], g["tau"])
    lines.append(f"    {'r/a':>6s} {'t':>8s} {'tau':>9s} {'H':>9s} {'H/W':>7s} {'s/W':>7s} {'rho':>11s} {'|rho|/H^2':>10s} {'p_par/rho':>9s} {'p_perp/rho':>10s} {'w_bar':>8s} {'m/(r^{n+1}H^2/2)':>17s} {'e^-f':>9s}")
    for r0 in (1.1, 1.5, 2.0, 3.0, 5.0, 10.0, 25.0, 50.0, 100.0):
        jj = np.where(g["r"] >= r0)[0]
        if not len(jj):
            continue
        j = jj[0]
        lines.append(f"    {g['r'][j]:6.2f} {g['t'][j]:8.2f} {g['tau'][j]:9.4f} {g['H'][j]:9.4f} {g['H'][j] / W:7.3f} {sp[j] / W:7.3f} {g['rho'][j]:+11.3e} "
                     f"{abs(g['rho'][j]) / g['H'][j] ** 2:10.4f} {g['ppar'][j] / g['rho'][j]:9.4f} {g['pperp'][j] / g['rho'][j]:10.4f} {g['w'][j]:8.5f} "
                     f"{g['m'][j] / (g['r'][j] ** (n + 1) * g['H'][j] ** 2 / 2):17.5f} {math.exp(-g['f'][j]):9.2e}")
    fits = []
    for rmax in (1.5, 2.0, 3.0, 5.0, 10.0, 25.0):
        m = (g["r"] - 1 >= 1e-2) & (g["r"] <= rmax)
        if m.sum() >= 10 and g["r"].max() > rmax:
            fits.append(f"r/a <= {rmax:g}: {np.polyfit(g['tau'][m], np.log(g['r'][m] - 1), 1)[0]:.4f} a^-1 = {np.polyfit(g['tau'][m], np.log(g['r'][m] - 1), 1)[0] / W:.3f} W")
    lines.append("    least-squares exponent of ln(r/a - 1) against tau over 1e-2 <= r/a - 1 and r/a <= r_max (the form of [SH02]'s fit; 1.1 a^-1 = " + f"{1.1 / W:.3f} W, W = {W:.4f} a^-1): " + "; ".join(fits))
    j = len(g["r"]) - 1
    lines.append(f"    last live diagonal cell: t = {g['t'][j]:.2f}, tau = {g['tau'][j]:.4f}, r/a = {g['r'][j]:.2f} (ln r/a = {math.log(g['r'][j]):.3f}), H = {g['H'][j]:.4f}")


def horizon_table(M, lines):
    R, n, h = M["R"], M["n"], M["h"]
    N = R["N"]
    rows = []
    for xm in (3.0, 4.0, 6.0, 10.0, 20.0, 50.0, 100.0, 150.0, 199.0):
        j = int(round(xm / h - 0.5))
        if j >= N - 1:
            continue
        xl, rl, xr, rr = R["xL"][j], R["rL"][j], R["xR"][j], R["rR"][j]
        rows.append(f"x- = {xm:6.1f}: r_+ = 0 at x+ = {xl:8.3f} (r = {rl:8.3f})" + (f"; r_- = 0 at x+ = {xr:9.3f} (r = {rr:8.3f})" if np.isfinite(xr) else "; r_- = 0 not inside the grid (x+ > X)"))
    jr = np.where(np.isfinite(R["xR"]))[0]
    lines.append(f"    horizon contours at fixed x-: " + ("; ".join(rows[:3]) if rows else "none"))
    for x in rows[3:]:
        lines.append("      " + x)
    if len(jr):
        lines.append(f"    the contour r_- = 0 stays inside the grid up to x- = {(jr[-1] + 0.5) * h:.3f} (x+ = {R['xR'][jr[-1]]:.1f} there, r = {R['rR'][jr[-1]]:.1f}); beyond it the whole row x+ <= X is anti-trapped to the right of the contour r_+ = 0")


def mass_table(M, lines):
    R, n, h = M["R"], M["n"], M["h"]
    N = R["N"]
    cols = sorted(R["cols"])
    ends = {i0: end_mass(R, i0) for i0 in cols}
    Mt, ft, It, rt = top_mass(R)
    lines.append(f"    {'x-':>6s} {'t at x+=X':>9s} " + " ".join(f"{'M(x+=' + format(i0 * h, '.0f') + ')':>13s}" for i0 in cols) + f" {'M(x-=X, x+)':>13s} {'m_diag':>12s} {'m_edge':>12s} {'I_pos':>12s} {'I_neg':>12s}")
    for xm in CRIT["MB"]["ray_x"]:
        j = int(round(xm / h))
        jc = int(round(xm / h - 0.5))
        if j >= N - 1 or not np.isfinite(R["m_diag"][jc]):
            continue
        lines.append(f"    {xm:6.1f} {(N * h + xm) / SQ2:9.2f} " + " ".join(f"{ends[i0][0][j]:13.5e}" for i0 in cols) + f" {Mt[j]:13.5e} {R['m_diag'][jc]:12.4e} {R['m_edge'][jc]:12.4e} {R['I_pos'][jc]:12.4e} {R['I_neg'][jc]:12.4e}")
    M0 = ends[N][0][0]
    sol = dn.ellis_dense(n, N * h / SQ2 + 1.0)[0]
    lines.append(f"    M at the start of the edge (x- = 0), x+ = X: {M0:+.6f}; the static wormhole's 1/(2 r^(n-1)) there is {1 / (2 * float(sol(N * h / SQ2)[0]) ** (n - 1)):.6f}, so the pulse's E0 = {M0 - 1 / (2 * float(sol(N * h / SQ2)[0]) ** (n - 1)):+.6f}")


def analyze_suite(suite, scale):
    runs, ctrls = suite["runs"], suite["ctrl"]
    print(f"expansion.py analyze: scale {scale}, mutation {EMUT or dn.MUT}, numba {dn.nb.__version__ if dn.HAVE_NUMBA else 'absent'}, numpy {np.__version__}, {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"runs: " + ", ".join(f"D{M['n'] + 2} h={M['h']:.5f} {M['member']} (N={M['R']['N']})" for M in runs))
    if suite.get("build_error"):
        check("BUILD", "the runs of the suite complete", False, suite["build_error"])
    for fn in (chk_ST, chk_OB0, chk_W1, chk_PL):
        try:
            fn()
        except Exception as e:                      # an exception inside a check is that check's failure
            import traceback
            traceback.print_exc()
            check(fn.__name__[4:], "check raised an exception", False, repr(e))
    # ---- per-run observer, identities, linear-phase slope
    rows = [chk_member(M) for M in runs]
    for lab, asym, off, rll, dpp, dw, dpl, dT in rows:
        print(f"    {lab}: mirror asymmetry {asym:.2e}; slice offset {off} nodes; min r_ll on the diagonal {rll:+.3e}; |p_perp/rho + 1| {dpp:.1e}, |w_bar - (1-n)/(1+n)| {dw:.1e}, |p_par/rho - 1| {dpl:.1e}, |T_+-/rho| {dT:.1e}")
    if runs:
        check("OB1", "the evolved fields are mirror symmetric at the diagonal (r, f, phi odd) through the whole run",
              max(r_[1] for r_ in rows) <= CRIT["OB_sym"], f"max asymmetry {max(r_[1] for r_ in rows):.2e}, tol {CRIT['OB_sym']:g}")
        check("OB2", "the diagonal is the minimum of r on every full-length slice t = const (node level) and a local minimum (r_ll > 0) through the whole run",
              max(r_[2] for r_ in rows) <= CRIT["OB_slice"] and min(r_[3] for r_ in rows) > 0,
              f"largest offset {max(r_[2] for r_ in rows)} nodes (tol {CRIT['OB_slice']}), smallest r_ll {min(r_[3] for r_ in rows):+.3e}")
        wmax = max(max(r_[4:]) for r_ in rows)
        check("W2", "on the symmetric diagonal p_perp = -rho, p_par = rho, T_+- = 0, so w_bar = (1-n)/(1+n) (identities of the mirror symmetry)", wmax <= CRIT["W2_tol"],
              f"largest deviation {wmax:.2e}, tol {CRIT['W2_tol']:g}; w_bar = {(1 - runs[0]['n']) / (1 + runs[0]['n']):.4f} at D = {runs[0]['n'] + 2}")
    tuned = [M for M in runs if M["member"] == "tuned"]
    if tuned:
        fits = []
        for M in tuned:
            s, npts = sl0(M)
            fits.append((M, s, npts))
            print(f"    SL0 D{M['n'] + 2} h={M['h']:.5f}: p*_h = {M['meta'].get('pstar', float('nan')):.6e}, A = p*_h + {CRIT['delta']:g}; slope of ln|H| {s:.6f} ({npts} pts) against W = {W_EXACT[M['n']]:.6f}: {s / W_EXACT[M['n']] - 1:+.2e}")
        check("SL0", "control: the tuned member departs from the throat at the linear rate W (slope of ln|H| over |H|/W in [1e-3, 1e-2]) at every h",
              all(abs(s / W_EXACT[M["n"]] - 1) <= CRIT["SL0"]["bar"] for M, s, _ in fits), f"bar {CRIT['SL0']['bar']:.0%}")
    # ---- classification, slope test, anti-trapped region, series, masses
    classes = {}
    for M in runs:
        g = diag_arrays(M["R"])
        W = W_EXACT[M["n"]]
        cl = classify(g["tau"], g["r"], g["H"], g["rho"], g["w"])
        met, vals = slope_test(g["r"], g["H"], W)
        M["cls"], M["slope"] = cl, (met, vals)
        classes[(M["n"], M["member"], M["h"])] = cl["class"]
        print(f"  --- D{M['n'] + 2}, h = {M['h']:.5f}, member {M['member']} (A = {M['A']:.6e}) ---")
        lines = []
        series_table(M, lines)
        print("\n".join(lines))
        print(f"    final window (last e-fold of r, {cl['npts']} pts, r/a = {cl['r0']:.1f} to {cl['r1']:.1f}): H max/min {cl['Hratio']:.3f}, H end/start {cl['Hgrow']:.3f}, |rho| max/min {cl['rhoratio']:.3f}, "
              f"w_bar in [{cl['wmin']:.4f}, {cl['wmax']:.4f}]; 1/H linear in tau: k = {cl['k']:.4f}, tau_s = {cl['tau_s']:.4f} (last point {cl['tau_end']:.4f}), max residual {cl['ft_resid']:.1e}  -> class {cl['class']}")
        print(f"    slope test (#211): H/W at r/a = " + ", ".join(f"{r0:g}: {v:.3f}" for r0, v in vals.items()) + f"; bar {CRIT['SL']['bar']:.0%}: {'met' if met else 'not met'}")
        R = M["R"]
        ja = np.where(R["cls_diag"] == 1)[0]
        if len(ja):
            j0 = ja[0]
            alive = np.where(R["alive_diag"])[0]
            pers = bool(np.all(R["cls_diag"][alive[alive >= j0]] == 1))
            print(f"    anti-trapped region: first anti-trapped diagonal cell at t = {R['t'][j0]:.3f}, tau = {R['tau'][j0]:.3f}, r/a = {R['r'][j0]:.3f}; every later live diagonal cell anti-trapped: {pers}; trapped cells anywhere: {int(R['n_trap'].sum())}; "
                  f"excised cells: {int(R['n_dead'].max())}; first excised cell: {R['first_dead']}")
        else:
            print("    anti-trapped region: no anti-trapped diagonal cell")
        print(f"    gauge: max e^-f over live cells {np.nanmax(R['emf_max']):.4f}, min e^-f {np.nanmin(R['emf_min']):.2e}, ln r_th max {math.log(np.nanmax(R['r'])):.3f}; "
              f"constraints: max |C++| {np.nanmax(R['Cpp']):.2e}, relative to its terms {np.nanmax(R['Cpp_rel']):.2e}; max |C--| {np.nanmax(R['Cmm']):.2e}, relative {np.nanmax(R['Cmm_rel']):.2e}")
        lines = []
        horizon_table(M, lines)
        print("\n".join(lines))
        lines = []
        mass_table(M, lines)
        print("\n".join(lines))
    # ---- masses
    worst = dict(flux=0.0, mirror=0.0, ray=0.0)
    for M in runs:
        R = M["R"]
        for i0 in R["cols"]:
            Mj, flux, I, resid = end_mass(R, i0)
            worst["flux"] = max(worst["flux"], float(np.max(np.abs(resid)) / np.max(np.abs(Mj - Mj[0]))))
        Mt, ft, It, rt = top_mass(R)
        worst["flux"] = max(worst["flux"], float(np.max(np.abs(rt)) / np.max(np.abs(Mt - Mt[0]))))
        Me = end_mass(R, R["N"])[0]
        worst["mirror"] = max(worst["mirror"], float(np.max(np.abs(Mt - Me)) / np.max(np.abs(Me))))
        h = M["h"]
        for xm in CRIT["MB"]["ray_x"]:
            jc = int(round(xm / h - 0.5))
            if jc < R["N"] - 1 and np.isfinite(R["m_diag"][jc]):
                scale_ = abs(R["I_pos"][jc]) + abs(R["I_neg"][jc])
                worst["ray"] = max(worst["ray"], abs(R["m_edge"][jc] - R["m_diag"][jc] - R["I_tot"][jc]) / scale_)
    byk = {}
    for M in runs:
        Rr = M["R"]
        Mj, flux, I, resid = end_mass(Rr, Rr["N"])
        byk.setdefault((M["n"], M["member"]), []).append((M["h"], float(np.max(np.abs(resid))), float(np.max(np.abs(Mj - Mj[0])))))
    for key, v in sorted(byk.items(), key=lambda kv: (kv[0][0], str(kv[0][1]))):
        v = sorted(v, key=lambda x: -x[0])
        print(f"    MB1 report D{key[0] + 2} {key[1]}, edge x+ = X: residual max|M(j) - M(0) - int flux| " + ", ".join(f"{x[1]:.3e} (h = {x[0]:.5f}, largest M change {x[2]:.3e})" for x in v)
              + (f"; ratio {v[0][1] / v[1][1]:.2f}" if len(v) > 1 and v[1][1] > 0 else ""))
    if runs:
        c = CRIT["MB"]
        check("MB1", "the change of the mass at each end equals the flux through it: M(j) - M(0) = int flux dx-, on the edge x+ = X and the columns X/2, X/4, and on the edge x- = X",
              worst["flux"] <= c["flux"], f"largest residual / largest mass change {worst['flux']:.2e}, tol {c['flux']:g}")
        check("MB2", "the two ends have the same mass history (mirror symmetry of the masses): M(i, N) = M(N, i)", worst["mirror"] <= c["mirror"], f"largest |M_top - M_edge| / max|M| {worst['mirror']:.2e}, tol {c['mirror']:g}")
        check("MB3", "the mass on the diagonal equals the mass at the edge cell minus the flux integral along the null ray between them, at x- in the declared list",
              worst["ray"] <= c["ray"], f"largest residual / (|I_pos| + |I_neg|) {worst['ray']:.2e}, tol {c['ray']:g}")
    # ---- two resolutions
    groups = {}
    for M in runs:
        groups.setdefault((M["n"], M["member"]), []).append(M)
    cres = CRIT["RES"]
    allok, lines = True, []
    for key, ms in sorted(groups.items(), key=lambda kv: (kv[0][0], str(kv[0][1]))):
        ms = sorted(ms, key=lambda M: -M["h"])
        for a, b in zip(ms[:-1], ms[1:]):
            ga, gb = diag_arrays(a["R"]), diag_arrays(b["R"])
            ca, cb = a["cls"], b["cls"]
            dH = []
            for r0 in (10.0, 25.0, 50.0):
                ja, jb = np.where(ga["r"] >= r0)[0], np.where(gb["r"] >= r0)[0]
                if len(ja) and len(jb):
                    dH.append(abs(np.interp(r0, ga["r"], ga["H"]) / np.interp(r0, gb["r"], gb["H"]) - 1))
            dk, dts = abs(ca["k"] / cb["k"] - 1), abs(ca["tau_s"] / cb["tau_s"] - 1)
            dl = abs(math.log(ga["r"].max()) - math.log(gb["r"].max()))
            ok = ca["class"] == cb["class"] and (max(dH) if dH else 0) <= cres["H_rel"] and dk <= cres["k_rel"] and dts <= cres["ts_rel"] and dl <= cres["lnr"]
            allok &= ok
            lines.append(f"D{key[0] + 2} {key[1]}: h {a['h']:.5f} vs {b['h']:.5f}: class {ca['class']} / {cb['class']}; H at fixed r/a (10, 25, 50) differs by {fmt(dH)}; k {dk:.1e}; tau_s {dts:.1e}; |d ln r_max| {dl:.1e}")
    for x in lines:
        print("    RES " + x)
    if lines:
        check("RES", "two resolutions agree: same class, H at fixed r/a within 1 %, k and tau_s within 1 %, e-folds within 1e-3", allok, f"bars {cres}")
    # ---- controls
    if ctrls:
        ok, lines = True, []
        by = {}
        for c_ in ctrls:
            by.setdefault((c_["n"], c_["name"]), []).append(c_)
        for key, cs in sorted(by.items()):
            cs = sorted(cs, key=lambda c_: -c_["h"])
            col = all(c_["fate"] == "collapse" for c_ in cs)
            dt = abs(cs[0].get("t_event", float("nan")) - cs[-1].get("t_event", float("nan"))) if len(cs) > 1 else 0.0
            tol = CRIT["CT"]["dt"] * cs[0]["h"] * SQ2
            pos = (cs[0]["E0"] > 0) if key[1] == "normal" else True
            ok &= col and dt <= tol and pos
            lines.append(f"D{key[0] + 2} {key[1]} (A = {cs[0]['A']}): " + "; ".join(f"h = {c_['h']:.5f}: {c_['fate']} at t = {c_.get('t_event', float('nan')):.3f}, r = {c_.get('r_first_trapped', float('nan')):.4f}, M = {c_.get('M_first_trapped', float('nan')):.4f}, E0 = {c_['E0']:+.5f}" for c_ in cs) + f"; |dt| {dt:.3e} (tol {tol:.3e})")
        for x in lines:
            print("    CT " + x)
        check("CT", "control: a positive-energy pulse (a normal field, same Gaussian, mirrored) and the S ghost pulse of the other sign collapse to a black hole, at the same time at two resolutions",
              ok, f"first trapped cell, event times within {CRIT['CT']['dt']} sqrt2 h; the normal pulse has E0 > 0")
    bad = [c for c, ok in results if ok is False]
    skipped = [c for c, ok in results if ok is None]
    print(f"\n{sum(1 for _, ok in results if ok)} of {len(results)} checks pass; {len(bad)} fail" + (f": {', '.join(bad)}" if bad else "") + (f"; skipped: {', '.join(skipped)}" if skipped else ""))
    print("classes: " + "; ".join(f"D{k[0] + 2} {k[1]} h={k[2]:.5f}: {v}" for k, v in sorted(classes.items(), key=lambda kv: (kv[0][0], str(kv[0][1]), -kv[0][2]))))
    return bad


def load_suite(raw):
    runs, ctrls = [], []
    d = os.path.join(raw, "runs")
    for f in sorted(os.listdir(d)):
        path = os.path.join(d, f)
        if f.endswith(".npz") and f.startswith("D"):
            R = load_run(path)
            m = R["meta"]
            runs.append(dict(n=m["n"], h=m["h"], member=m["member"], A=m["A"], R=R, meta=m))
        elif f.startswith("ctrl_") and f.endswith(".json"):
            ctrls.append(json.load(open(path)))
    return dict(runs=runs, ctrl=ctrls)


def mini_suite():
    """The suite at reduced size for the mutation runs: D = 4, X = 40, h = 0.04, 0.02; members tuned and A = 0.05; controls at h = 0.04, 0.02, X = CT X.
    A failure to build part of it (a mutation can break the bisection) is recorded and reported as the check BUILD."""
    runs, ctrls, errs = [], [], []
    for h in (0.04, 0.02):
        for member in ("tuned", 0.05):
            try:
                R, meta = run_member(2, h, 40.0, member)
                runs.append(dict(n=2, h=h, member=str(member), A=meta["A"], R=R, meta=meta))
            except Exception as e:
                errs.append(f"member {member} h={h}: {e!r}")
    c = CRIT["CT"]
    for h in (0.04, 0.02):
        for name, d in (("normal", normal_data(2, h, c["X"], c["A_normal"])), ("ghost_neg", s_data(2, h, c["X"], c["A_ghost"]))):
            try:
                res = dn.evolve(d, kernel="numba", tol=CRIT["trap_tol"], stop=False, store=True)
                info = control_info(res)
                info.update(name=name, n=2, h=h, X=c["X"], A=(c["A_normal"] if name == "normal" else c["A_ghost"]))
                ctrls.append(info)
            except Exception as e:
                errs.append(f"control {name} h={h}: {e!r}")
    return dict(runs=runs, ctrl=ctrls, build_error="; ".join(errs) or None)


def main():
    global EMUT
    argv = sys.argv[1:]
    if not argv:
        sys.exit(__doc__)
    cmd, rest = argv[0], argv[1:]
    mut = None
    raw = None
    for x in rest:
        if x.startswith("--mutate="):
            mut = x.split("=", 1)[1]
        elif x.startswith("--raw="):
            raw = x.split("=", 1)[1]
    if mut:
        known = set(MUTATIONS) | set(DN_MUT) | {"mirror_parity", "w_sign", "time_unit", "class_loose", "control_ghost"}
        if mut not in known:
            sys.exit(f"unknown mutation {mut}; known: {', '.join(sorted(known))}")
        if mut in DN_MUT:
            dn.MUT = mut
        else:
            EMUT = mut
    if cmd == "job":
        job([x for x in rest if not x.startswith("--mutate=")])
    elif cmd == "analyze":
        bad = analyze_suite(load_suite(raw), "record")
        sys.exit(1 if bad else 0)
    elif cmd == "mutate":
        suite = mini_suite()
        bad = analyze_suite(suite, "mini")
        want = set(MUTATIONS.get(mut, []))
        fired = {c.split("[")[0] for c in bad}
        print(f"mutation {mut}: declared target {sorted(want)}; fired {sorted(fired)}; " + ("targets all fired" if want <= fired else f"SILENT: {sorted(want - fired)}") + ("" if fired <= want else f"; extra: {sorted(fired - want)}"))
    elif cmd == "dev":
        suite = mini_suite()
        analyze_suite(suite, "mini")
    else:
        sys.exit(f"unknown command {cmd}")


if __name__ == "__main__":
    main()

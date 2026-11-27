#!/usr/bin/env python3
"""Development prototype of the #217 evolution (design stage): the #39 kernel with the potential term, data of a static member on the
rays in the tortoise null gauge, a one-sided pulse, the event finders. Not the build; the build is evolve.py after the design is approved.

    python3 evolve_dev.py static C [--h=H] [--X=X] [--rho_c=R]    static check against the exact member (error at the nodes)
    python3 evolve_dev.py pulse C --A=A [--h=H] [--X=X] [--ray=+|-] [--rho_c=R]    one run with a one-sided pulse; events and timing
"""
import math
import os
import sys
import time

import numba as nb
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "wormhole_evolution"))
sys.path.insert(0, HERE)
import dn                      # noqa: E402
import evolve_potential as ep  # noqa: E402

SQ2 = math.sqrt(2.0)


# --------------------------------------------------------------------------------------------------------------------
# kernel: dn.py's with the potential term
# --------------------------------------------------------------------------------------------------------------------
@nb.njit(cache=True)
def _F(c, up, um, par, eps, g):
    n = par[0]
    th = par[1]
    sc = par[2]
    pms = par[3]
    nf = c.shape[0]
    emf = math.exp(-c[1])
    rpm = -th * (up[0] * um[0] + 0.5 * emf) / c[0]
    hasV = par[5] != 0.0
    Vv = 0.0
    Vp = 0.0
    if hasV:
        Vv, Vp = ep.pot(c[2], par[6], par[7])
        rpm = rpm + 0.5 * Vv * c[0] * emf
    spm = 0.0
    for k in range(nf - 2):
        spm = spm + eps[k] * up[2 + k] * um[2 + k]
    g[0] = rpm
    g1 = (n / c[0]) * rpm + pms * spm
    if hasV:
        g1 = g1 - Vv * emf
    g[1] = g1
    kk = -sc / (2.0 * c[0])
    for k in range(2, nf):
        g[k] = kk * (up[0] * um[k] + um[0] * up[k])
    if hasV:
        g[2] = g[2] - 0.5 * eps[0] * emf * Vp


@nb.njit(cache=True)
def row_step(old, new, dead_old, dead_new, h, par, eps, rcut, fcut, excise):
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


def params(mem, V=True):
    par = np.zeros(8)
    par[:5] = dn.kernel_params(2, [-1.0])
    par[5] = 1.0 if V else 0.0
    par[6], par[7] = mem.c, mem.K
    return par


# --------------------------------------------------------------------------------------------------------------------
# data: the member on the two rays through the corner, tortoise null gauge (e^{-f} = A, l = (x+ - x-)/sqrt2 = tortoise coordinate)
# --------------------------------------------------------------------------------------------------------------------
def tortoise_rho(mem, rho_c, xs):
    """rho, z, A, f at the tortoise coordinates xs (a sorted array containing 0 or starting at 0, either direction) with rho(0) = rho_c.
    Near a horizon the variable z = rho - rho_h is integrated through zeta = ln z, so A and f keep full relative accuracy."""
    from scipy.integrate import solve_ivp
    xs = np.asarray(xs, float)
    out_rho = np.empty_like(xs)
    out_A = np.empty_like(xs)
    out_f = np.empty_like(xs)
    zc = None
    if mem.rh is not None and rho_c > mem.rh:
        zc = rho_c - mem.rh
    for sgn in (+1, -1):
        sel = np.where(sgn * xs > 0)[0]
        if len(sel) == 0:
            continue
        x = xs[sel]
        order = np.argsort(sgn * x)
        x = x[order]
        use_zeta = zc is not None and sgn < 0
        if use_zeta:
            sol = solve_ivp(lambda t, y: mem.a_over_z(np.exp(y)), [0.0, x[-1]], [math.log(zc)], method="DOP853", rtol=1e-13, atol=1e-14, t_eval=x)
            zeta = sol.y[0]
            z = np.exp(zeta)
            a = mem.a_over_z(z)
            rho = mem.rh + z
            A = z * a
            f = -zeta - np.log(a)
        else:
            sol = solve_ivp(lambda t, y: mem.A(y), [0.0, x[-1]], [rho_c], method="DOP853", rtol=1e-13, atol=1e-14, t_eval=x)
            rho = sol.y[0]
            A = mem.A(rho)
            f = -np.log(A)
        if not sol.success or len(rho) != len(x):
            raise RuntimeError("tortoise integration failed: " + sol.message)
        idx = sel[order]
        out_rho[idx], out_A[idx], out_f[idx] = rho, A, f
    z0 = np.where(xs == 0)[0]
    if len(z0):
        out_rho[z0] = rho_c
        out_A[z0] = float(mem.A(rho_c))
        out_f[z0] = -math.log(float(mem.A(rho_c)))
    return out_rho, out_A, out_f


def member_rays(mem, rho_c, N, h, pulse=None):
    """dn.Data of the static member with the corner at rho_c (x+ = x- = 0). pulse: dict(A, xc, w, ray) adds A (g - g0) to phi on a ray
    and solves the constraint for the change of r (RK4 on the deviation, which is exactly zero at A = 0)."""
    s = np.arange(N + 1) * h
    rays = {}
    for ray, sg in (("+", +1.0), ("-", -1.0)):
        rho, A, f = tortoise_rho(mem, rho_c, sg * s / SQ2)
        U = np.zeros((3, N + 1))
        U[0] = mem.R(rho)
        U[1] = f
        U[2] = mem.phi(rho)
        if pulse is not None and pulse["ray"] == ray and pulse["A"] != 0.0:
            # d phi0/ds along the ray: phi0'(rho) A sg / sqrt2, phi0' = sqrt2/(rho^2+1); f_s = -A'(rho) sg/sqrt2
            gc = lambda ss: np.exp(-((ss - pulse["xc"]) / pulse["w"]) ** 2)
            g0 = math.exp(-(pulse["xc"] / pulse["w"]) ** 2)
            U[2] = U[2] + pulse["A"] * (gc(s) - g0)

            def at(ss):
                rh_, A_, _ = tortoise_rho(mem, rho_c, sg * ss / SQ2)
                ph0s = (SQ2 / (rh_ * rh_ + 1.0)) * A_ * sg / SQ2
                ps = ph0s + pulse["A"] * (-2.0 * (ss - pulse["xc"]) / pulse["w"] ** 2) * gc(ss)
                return mem.R(rh_), ph0s, ps, -mem.Ap(rh_) * sg / SQ2
            Rn, ph0n, psn, fsn = at(s)
            sh = s[:-1] + 0.5 * h
            Rh, ph0h, psh, fsh = at(sh)
            # deviation dl = r - R(rho(s)): dl'' + f_s dl' - (1/2) phi_s^2 dl = (1/2) R (phi_s^2 - phi0_s^2)   (n = 2, eps = -1)
            S0, Sh = 0.5 * Rn * (psn ** 2 - ph0n ** 2), 0.5 * Rh * (psh ** 2 - ph0h ** 2)
            Q0, Qh = 0.5 * psn ** 2, 0.5 * psh ** 2
            dl = np.zeros(N + 1)
            dv = np.zeros(N + 1)

            def d(y, S, Q, fq):
                return (y[1], S + Q * y[0] - fq * y[1])
            for i in range(N):
                y = (dl[i], dv[i])
                k1 = d(y, S0[i], Q0[i], fsn[i])
                k2 = d((y[0] + 0.5 * h * k1[0], y[1] + 0.5 * h * k1[1]), Sh[i], Qh[i], fsh[i])
                k3 = d((y[0] + 0.5 * h * k2[0], y[1] + 0.5 * h * k2[1]), Sh[i], Qh[i], fsh[i])
                k4 = d((y[0] + h * k3[0], y[1] + h * k3[1]), S0[i + 1], Q0[i + 1], fsn[i + 1])
                dl[i + 1] = y[0] + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
                dv[i + 1] = y[1] + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            U[0] = U[0] + dl
        rays[ray] = U
    return dn.Data(2, [-1.0], N, h, rays["+"], rays["-"], name="member")


# --------------------------------------------------------------------------------------------------------------------
# evolution with row diagnostics and events
# --------------------------------------------------------------------------------------------------------------------
def evolve(d, par, rcut=0.05, fcut=25.0, tol=1e-2, rexp=2.0, stop=True, store=False, track=None):
    """Row sweep. Events: first row with a trapped cell (collapse) or an anti-trapped cell with r >= rexp (expansion); with stop the run
    ends at the first of them. Per row: r and phi at the cell l = (i - j) h/sqrt2 for the offsets in track, min r over alive cells."""
    N, h, nf, eps = d.N, d.h, d.nf, d.eps
    old = np.ascontiguousarray(d.Sp.copy())
    new = np.empty_like(old)
    dold = np.zeros(N + 1, np.uint8)
    dnew = np.zeros(N + 1, np.uint8)
    G = np.empty((nf, N + 1, N + 1)) if store else None
    dead = np.zeros((N + 1, N + 1), np.uint8) if store else None
    if store:
        G[:, :, 0] = old
    ev = {}
    rows = 0
    rmin = np.full(N, np.nan)
    trk = {k: np.full((N, 3), np.nan) for k in (track or [])}
    maxpsi = 0.0
    for j in range(N):
        new[:, 0] = d.Sm[:, j + 1]
        dnew[0] = 0
        row_step(old, new, dold, dnew, h, par, eps, rcut, fcut, 1)
        c, up, um = dn.strip_cells(old, new, h)
        alive = ~(dold[:-1].astype(bool) | dold[1:].astype(bool) | dnew[:-1].astype(bool) | dnew[1:].astype(bool))
        if store:
            G[:, :, j + 1] = new
            dead[:, j + 1] = dnew
        if alive.any():
            rmin[j] = c[0][alive].min()
            maxpsi = max(maxpsi, float(np.abs(c[2][alive]).max()) / SQ2)
        for k in trk:
            i = j + k
            if 0 <= i < N and alive[i]:
                trk[k][j] = (c[0, i], c[1, i], c[2, i])
        rows = j + 1
        if not ev:
            cc = dn.cell_class(up[0], um[0], tol)
            trapped = alive & (cc == -1)
            anti = alive & (cc == 1) & (c[0] >= rexp)
            if trapped.any():
                ev["collapse"] = (j, int(np.argmax(trapped)))
            elif anti.any():
                ev["expand"] = (j, int(np.argmax(anti)))
            if ev and stop:
                break
        old, new = new, old
        dold, dnew = dnew, dold
    res = dict(events=ev, rows=rows, rmin=rmin[:rows], track={k: v[:rows] for k, v in trk.items()}, maxpsi=maxpsi, G=G, dead=dead)
    return res


class Static:
    """exact static member at the cell centres, indexed by k = i - j + N (cell-averaged like the code's centre values), and the weight of Y_inv."""

    def __init__(self, mem, rho_c, N, h):
        ks = np.arange(-N - 1, N + 2)
        rho, A, f = tortoise_rho(mem, rho_c, ks * h / SQ2)
        Rn, Pn = mem.R(rho), mem.phi(rho)
        self.N = N
        sl = slice(1, -1)
        self.R = 0.25 * (2 * Rn[sl] + Rn[2:] + Rn[:-2])
        self.P = 0.25 * (2 * Pn[sl] + Pn[2:] + Pn[:-2])
        self.Q = (rho * Rn / SQ2)[sl]                         # R'/phi0' = rho R / sqrt2
        self.f = 0.25 * (2 * f[sl] + f[2:] + f[:-2])


def evolve_dep(d, par, st, Y0=1e-2, rcut=0.05, fcut=25.0, tol=1e-2, stop_dep=True, run_events=False, wtol=0.02, off=0):
    """Row sweep that follows the gauge-invariant departure Y = (r - r_s) - (rho R/sqrt2)(phi - phi_s) (cell centres, static values cell-averaged):
    the first row with max |Y| >= Y0 gives the departure time (t of the diagonal cell centre of that row) and the sign of Y at the maximum.
    With run_events the f-weighted finders also run: trapped = r_+ < 0, r_- < 0, 2 e^f r_+ r_- >= wtol; anti-trapped likewise with both > 0."""
    N, h, nf, eps = d.N, d.h, d.nf, d.eps
    old = np.ascontiguousarray(d.Sp.copy())
    new = np.empty_like(old)
    dold = np.zeros(N + 1, np.uint8)
    dnew = np.zeros(N + 1, np.uint8)
    dep = None
    ev = {}
    ymax = np.zeros(N)
    i_idx = np.arange(N)
    rows = 0
    for j in range(N):
        new[:, 0] = d.Sm[:, j + 1]
        dnew[0] = 0
        row_step(old, new, dold, dnew, h, par, eps, rcut, fcut, 1)
        c, up, um = dn.strip_cells(old, new, h)
        alive = ~(dold[:-1].astype(bool) | dold[1:].astype(bool) | dnew[:-1].astype(bool) | dnew[1:].astype(bool))
        k = i_idx - j + N
        Y = (c[0] - st.R[k]) - st.Q[k] * (c[2] - st.P[k])
        ia = j + off
        ok = alive[ia] if 0 <= ia < N else False
        Ya = np.zeros(N)
        if ok:
            Ya[ia] = Y[ia]
        ymax[j] = abs(Ya[ia]) if ok else 0.0
        rows = j + 1
        if dep is None and ymax[j] >= Y0:
            # time of the cell centre (x+ + x-)/sqrt2 and position l
            dep = dict(j=j, i=ia, T=(d.x0p + d.x0m + (ia + j + 1.0) * h) / SQ2, l=(ia - j) * h / SQ2, sign=1 if Ya[ia] > 0 else -1, Y=float(Ya[ia]))
            if stop_dep and not run_events:
                break
        if run_events and not ev:
            wp = np.exp(c[1]) * up[0] * um[0] * 2.0
            trapped = alive & (up[0] < 0) & (um[0] < 0) & (wp >= wtol)
            anti = alive & (up[0] > 0) & (um[0] > 0) & (wp >= wtol)
            if trapped.any():
                ev["collapse"] = (j, int(np.argmax(trapped)))
            elif anti.any():
                ev["expand"] = (j, int(np.argmax(anti)))
        if dep is not None and (not run_events or ev):
            break
        if run_events and ev and not stop_dep:
            break
        old, new = new, old
        dold, dnew = dnew, dold
    return dict(dep=dep, events=ev, rows=rows, ymax=ymax[:rows])


def extract_mode(mem, rho_c, h, X0, A_p=0.0, Y_ext=2e-3, tmin=6.0, ray="+"):
    """Run the member (A_p = 0: no pulse, the departure is seeded by the truncation error; or with a small pulse A_p) until the departure |Y| at the diagonal first reaches Y_ext after tmin, then return the data on
    the two rays through the diagonal corner of that row (fields minus the static member's), normalized to Y_inv = 1 at the corner node.
    Returns (dR, dL, info): arrays (3, K+1) for the right ray (x+ increasing, j = j_s) and the left ray (x- increasing, i = i_s)."""
    N = int(round(X0 / h))
    d = member_rays(mem, rho_c, N, h, pulse=dict(A=A_p, xc=3.0, w=0.5, ray=ray) if A_p != 0.0 else None)
    st = Static(mem, rho_c, N, h)
    par = params(mem)
    nf = d.nf
    old = np.ascontiguousarray(d.Sp.copy())
    new = np.empty_like(old)
    dold = np.zeros(N + 1, np.uint8)
    dnew = np.zeros(N + 1, np.uint8)
    i_idx = np.arange(N)
    rows = {}
    col = []
    js = None
    for j in range(N):
        new[:, 0] = d.Sm[:, j + 1]
        dnew[0] = 0
        row_step(old, new, dold, dnew, h, par, d.eps, 0.05, 25.0, 1)
        if js is None:
            c, up, um = dn.strip_cells(old, new, h)
            tj = (2 * j + 1.0) * h / SQ2
            Y = (c[0][j] - st.R[N]) - st.Q[N] * (c[2][j] - st.P[N])
            if tj >= tmin and abs(Y) >= Y_ext:
                js = j
                rowS = old.copy()                  # row j (x- = j h): nodes i, includes the corner node (j, j)
                rowN = new.copy()                  # row j + 1
        if js is not None:
            col.append(new[:, js].copy())          # node (i = js, j + 1) for j >= js
        old, new = new, old
        dold, dnew = dnew, dold
    if js is None:
        raise RuntimeError("no departure reached")
    K = N - js
    # right ray: row js, i = js .. N ; left ray: column js, j = js .. N  (corner node (js, js) is in rowS[:, js])
    Gs = np.array([rowS[:, js]] + col).T                  # (nf, N - js + 1): column i = js at rows j = js .. N
    right = rowS[:, js:js + K + 1].copy()
    left = Gs[:, :K + 1].copy()
    # static member's data on rays through a corner at l = 0
    ds = member_rays(mem, rho_c, K, h)
    dR, dL = right - ds.Sp, left - ds.Sm
    ks = np.arange(K + 1)
    Q0 = st.Q[N]
    Yc = dR[0, 0] - Q0 * dR[2, 0]
    info = dict(js=js, t_s=(2 * js) * h / SQ2, K=K, Yc=float(Yc), maxdR=np.abs(dR).max(axis=1), maxdL=np.abs(dL).max(axis=1))
    return dR / Yc, dL / Yc, info


def mode_rays(mem, rho_c, N, h, A, dR, dL):
    """static member + A * (extracted mode): r, f, phi on both rays (N <= K)."""
    ds = member_rays(mem, rho_c, N, h)
    return dn.Data(2, [-1.0], N, h, ds.Sp + A * dR[:, :N + 1], ds.Sm + A * dL[:, :N + 1], name="mode")


def event_T(d, ev):
    key = next(iter(ev))
    j, i = ev[key]
    xp, xm = d.x0p + (i + 0.5) * d.h, d.x0m + (j + 0.5) * d.h
    return key, (xp + xm) / SQ2, (xp - xm) / SQ2


# --------------------------------------------------------------------------------------------------------------------
def static_check(c, rho_c, h, X):
    mem = ep.Member(c)
    N = int(round(X / h))
    d = member_rays(mem, rho_c, N, h)
    par = params(mem)
    t0 = time.time()
    res = evolve(d, par, stop=False, store=True, track=[0])
    el = time.time() - t0
    # exact nodes: l = (i - j) h / sqrt2
    ks = np.arange(-N, N + 1)
    rho, A, f = tortoise_rho(mem, rho_c, ks * h / SQ2) if mem.rh is None or rho_c > mem.rh else None
    i = np.arange(N + 1)[:, None]
    j = np.arange(N + 1)[None, :]
    kk = (i - j) + N
    G = res["G"]
    er = np.abs(G[0] - mem.R(rho)[kk])
    ef = np.abs(G[1] - f[kk])
    ep_ = np.abs(G[2] - mem.phi(rho)[kk])
    t = (i + j) * h / SQ2 + 0 * kk
    return dict(h=h, N=N, secs=el, er=er, ef=ef, ep=ep_, t=t, dead=res["dead"])


def main():
    args = sys.argv[1:]
    kw = {a.split("=")[0][2:]: a.split("=")[1] for a in args if a.startswith("--")}
    pos = [a for a in args if not a.startswith("--")]
    cmd, c = pos[0], float(pos[1])
    h = float(kw.get("h", 0.02))
    X = float(kw.get("X", 10.0))
    rho_c = float(kw.get("rho_c", 0.0))
    if cmd == "static":
        r = static_check(c, rho_c, h, X)
        t = r["t"]
        for tc in (2.0, 4.0, 8.0, 12.0, 1e9):
            m = (t <= tc) & (r["dead"] == 0)
            if m.any():
                print(f"c={c} h={h} t<={min(tc, t.max()):.1f}: max err r {r['er'][m].max():.3e} f {r['ef'][m].max():.3e} phi {r['ep'][m].max():.3e}")
        print(f"{r['secs']:.2f} s for N = {r['N']}")
    elif cmd == "pulse":
        A = float(kw["A"])
        mem = ep.Member(c)
        N = int(round(X / h))
        d = member_rays(mem, rho_c, N, h, pulse=dict(A=A, xc=float(kw.get("xc", 3.0)), w=float(kw.get("w", 0.5)), ray=kw.get("ray", "+")))
        t0 = time.time()
        res = evolve(d, params(mem), stop=True, tol=1e-2, rexp=2.0)
        print(c, A, res["events"], event_T(d, res["events"]) if res["events"] else None, f"rows {res['rows']}/{N}  {time.time() - t0:.2f} s  max|psi| {res['maxpsi']:.3f}")


if __name__ == "__main__":
    main()

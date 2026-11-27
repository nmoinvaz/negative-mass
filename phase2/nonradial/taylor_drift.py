#!/usr/bin/env python
"""#196: where scipy's polar D(gamma) drifts from #186's Rust Taylor crate as gamma grows.

    python taylor_drift.py star-steps OUT          the crate's star as Taylor steps (crate command `starsteps`), one JSON line per step
    python taylor_drift.py scipy --out=F.jsonl ... scipy's D (PolarProblem) for a list of cases, one JSON line each, resumable
    python taylor_drift.py crate --out-dir=DIR ... the crate's reference D (dd, f64, order 24) at (l, gamma, variant)
    python taylor_drift.py planted --out=F.jsonl   scipy's integrator on a planted system with a known answer
    python taylor_drift.py check [--mutate=NAME]   the checks (each can fail), on the files named with --dir=
    python taylor_drift.py report --dir=DIR        the tables

The scipy path is `polar_solver.PolarProblem` on `background.solve_star`'s star. `DriftProblem` below changes nothing in it
except what a case names: the star it reads (scipy's, or the crate's Taylor steps), M_ext, whether the exterior scalar plane
integrates in vacuum (as the crate does) and atol / the chunk count of `_integrate_plane`. With none named it is
`PolarProblem`, checked bit for bit (check IDENT).
"""
import argparse
import json
import math
import os
import signal
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

OMEGA, EPS = 0.7603, -1.0
VARIANTS = {"base": dict(phi_cut=1e-9), "phi": dict(phi_cut=1e-12)}   # convergence.py's base and `phi_cut 1e-12 (larger r_m)`
# r_m of PolarProblem on scipy's star (taylor_polar.py R_M); the crate's rule gives the same to the last bit (#187, RM)
RM = {"base": 24.10980364018201, "phi": 34.20514015700785}
NS = 30


# ---------------------------------------------------------------------------------------------------------------
# the crate's star as Taylor steps, in the interface PolarProblem reads
# ---------------------------------------------------------------------------------------------------------------
class CrateStar:
    """(alpha, a, phi0, u) from the crate's star: per step the coefficients 0..30 of (P, u, m, Z) = (phi0^2, phi0'/phi0, m, alpha^-2)
    about the step's x0 (`starsteps`), summed in f64. The 20-digit coefficients were produced in double-double."""

    def __init__(self, path, w=OMEGA, eps=EPS, r0=1e-3, R=60.0):
        x0, x1, C = [], [], []
        for line in open(path):
            o = json.loads(line)
            if o["cmd"] == "starstep":
                x0.append(o["x0"]); x1.append(o["x1"])
                C.append([[float(c) for c in f] for f in o["coef"]])
        self.x0, self.x1, self.C = np.array(x0), np.array(x1), np.array(C)
        self.w, self.eps, self.r0, self.R = w, eps, r0, R
        self.k = np.arange(NS + 1)
        assert self.C.shape[1:] == (4, NS + 1) and self.x1[-1] >= R, (self.C.shape, self.x1[-1])

    def state(self, r):
        r = np.atleast_1d(np.asarray(r, dtype=float))
        idx = np.clip(np.searchsorted(self.x0, r, side="right") - 1, 0, len(self.x0) - 1)
        h = r - self.x0[idx]
        return np.einsum("nfk,nk->fn", self.C[idx], h[:, None] ** self.k[None, :])

    def bg(self, r):
        r = np.atleast_1d(np.asarray(r, dtype=float))
        P, u, m, Z = self.state(r)
        return 1.0 / np.sqrt(Z), 1.0 / np.sqrt(1.0 - 2.0 * m / r), np.sqrt(P), u

    def m(self, r):
        return float(self.state(r)[2][0])


class CrateJ:
    """J(r; s) = J0(r) + s^2 J1(r) from the crate's tables (`jtable`): per table step the Taylor coefficients 0..29 of each nonzero entry about the
    step's x0, double-double rounded to 20 digits, summed in f64. The generated J_and_h of polar_matrix_generated.py (the same symbolic source)
    evaluated in double-double on the crate's star, rounded once per coefficient."""

    def __init__(self, path):
        x0, x1, c0, c1 = [], [], [], []
        for line in open(path):
            o = json.loads(line)
            if o["cmd"] == "jtable_head":
                self.pattern = [tuple(p) for p in o["pattern"]]
                self.nb = o["nb"]
            elif o["cmd"] == "jstep":
                x0.append(o["x0"]); x1.append(o["x1"])
                c0.append([[float(c) for c in e] for e in o["c0"]])
                c1.append([[float(c) for c in e] for e in o["c1"]])
        self.x0, self.x1, self.c0, self.c1 = np.array(x0), np.array(x1), np.array(c0), np.array(c1)
        self.k = np.arange(self.nb)
        self.ii = np.array([p[0] for p in self.pattern]); self.jj = np.array([p[1] for p in self.pattern])
        self._s = {}
        self.shift_end = False       # the control `cj_shift`: h measured from the step's end

    def J(self, r, s):
        key = float(s)
        if key not in self._s:
            self._s[key] = self.c0 + (key * key) * self.c1
        i = max(int(np.searchsorted(self.x0, r, side="right")) - 1, 0)
        h = r - (self.x1[i] if self.shift_end else self.x0[i])
        Jm = np.zeros((6, 6))
        Jm[self.ii, self.jj] = self._s[key][i] @ (h ** self.k)
        return Jm


# ---------------------------------------------------------------------------------------------------------------
# PolarProblem with the named changes
# ---------------------------------------------------------------------------------------------------------------
def make_problem_class():
    from scipy.integrate import solve_ivp
    from polar_solver import PolarProblem

    class DriftProblem(PolarProblem):
        """PolarProblem, with (all default: identical code path)
        atol         solve_ivp's atol in `_integrate_plane` (hard-coded 1e-14 there)
        extra_chunks added to `_integrate_plane`'s chunk count (the same problem, another chunk grid: a null control)
        vac_scalar   the exterior scalar plane integrates in vacuum (Schwarzschild of M_ext) at every r, as the crate does;
                     PolarProblem uses the star's J up to r = star.R (= 60), which is every r of that plane
        mext         M_ext replaced after construction (alpha_m and the lapse ratio stay the star's)
        noise_ulp    every entry of every J the integrator evaluates is multiplied by 1 + noise_ulp 2^-53 xi, xi standard normal
                     (numpy default_rng(noise_seed)): rounding noise relative to each entry, planted at a stated amplitude
        jtab         J(r; s) read from the crate's tables (`CrateJ`) in place of the generated J_and_h on the star's values (the interior and
                     the exterior scalar plane's J; the vacuum J of r > R is not replaced)
        chunk_factor the interior plane (the regular plane, from r0) is orthonormalized every 2/gamma divided by chunk_factor: the number of
                     renormalizations is chunk_factor times PolarProblem's
        riccati_rtol, riccati_atol   the tolerances of the Riccati integration of the exterior Zerilli q (`zerilli_q` hard-codes 1e-11 and 1e-14 and
                     ignores `rtol`)
        q_eps        q is multiplied by 1 + q_eps after the integration (a planted error of the exterior column, for D's sensitivity to it)
        noise_row    every nonzero entry J_ij gets the absolute noise noise_row 2^-53 xi max_j |J_ij| (xi standard normal, the row's largest entry):
                     the size `jround` measures for the generated J in f64 (rms 1.0 ulp of the row maximum)"""

        def __init__(self, star, l, rm=None, atol=1e-14, extra_chunks=0, vac_scalar=False, mext=None, noise_ulp=0.0, noise_seed=0, noise_row=0.0, jtab=None,
                     riccati_rtol=1e-11, riccati_atol=1e-14, q_eps=0.0, chunk_factor=1.0, **kw):
            super().__init__(star, l, rm=rm, **kw)
            self.chunk_factor = chunk_factor
            self.riccati_rtol, self.riccati_atol, self.q_eps = riccati_rtol, riccati_atol, q_eps
            self.atol, self.extra_chunks, self.vac_scalar = atol, extra_chunks, vac_scalar
            self.noise_ulp, self.noise_row, self.rng = noise_ulp, noise_row, np.random.default_rng(noise_seed)
            self.jtab = jtab
            if mext is not None:
                self.Mext = float(mext)

        def J(self, r, s):
            Jm = self.jtab.J(r, s) if self.jtab is not None else super().J(r, s)
            if self.noise_ulp:
                Jm = Jm * (1.0 + self.noise_ulp * 2.0 ** -53 * self.rng.standard_normal(Jm.shape))
            if self.noise_row:
                scale = np.abs(Jm).max(axis=1, keepdims=True)
                Jm = Jm + np.where(Jm != 0.0, self.noise_row * 2.0 ** -53 * scale * self.rng.standard_normal(Jm.shape), 0.0)
            return Jm

        def _integrate_plane(self, Y, ra, rb, s, nchunks=None, Jfun=None):
            k = Y.shape[1]
            Jfun = self.J if Jfun is None else Jfun
            if nchunks is None:
                nchunks = int(max(20, abs(rb - ra) * max(1.0, abs(s)) / 2.0))
            if self.chunk_factor != 1.0 and ra == self.r0:
                nchunks = max(1, int(nchunks * self.chunk_factor))
            nchunks += self.extra_chunks
            grid = np.linspace(ra, rb, nchunks + 1)
            Q, ph = self._qr(Y)
            for i in range(nchunks):
                sol = solve_ivp(lambda r, yv: (Jfun(r, s) @ yv.reshape(6, k)).ravel(), (grid[i], grid[i + 1]),
                                Q.ravel(), method="DOP853", rtol=self.rtol, atol=self.atol)
                Q, p = self._qr(sol.y[:, -1].reshape(6, k))
                ph *= p
            return Q, ph

        def zerilli_q(self, s, r_start=None):
            """PolarProblem.zerilli_q, line for line, with the integrator's tolerances and a planted multiplicative error as knobs (defaults: PolarProblem's)"""
            M, lam = self.Mext, (self.l - 1) * (self.l + 2) / 2.0

            def V(r):
                f = 1 - 2 * M / r
                return 2 * f * (lam ** 2 * r ** 2 * ((lam + 1) * r + 3 * M) + 9 * M ** 2 * (lam * r + M)) / (r ** 3 * (lam * r + 3 * M) ** 2)
            if r_start is None:
                r_start = self.rm + 40.0 / abs(s) + 20.0 * self.l / abs(s) + 50.0
            q0 = -np.sqrt(s * s + V(r_start))

            def rhs(rho, q):
                r = np.exp(rho); f = 1 - 2 * M / r
                return r * (s * s + V(r) - q * q) / f
            if np.iscomplexobj(q0) or np.iscomplexobj(s):
                q0 = complex(q0); q0 = q0 if q0.real <= 0 else -q0
                q0 = -np.sqrt(complex(s * s + V(r_start)))
                if (q0 / s).real > 0: q0 = -q0
            sol = solve_ivp(rhs, (np.log(r_start), np.log(self.rm)), [q0], method="DOP853", rtol=self.riccati_rtol, atol=self.riccati_atol)
            q = sol.y[0, -1] if np.iscomplexobj(sol.y) else float(sol.y[0, -1])
            return q * (1.0 + self.q_eps)

        def D(self, s):
            """PolarProblem.D, the same operations in the same order; the planes are kept (self.planes) for the principal angles"""
            Y, py = self.regular_plane(s)
            Z, pz = self.decaying_plane(s)
            self.planes = (Y, Z)
            d = np.linalg.det(np.hstack([Y, Z])) * py * pz
            return d if np.iscomplexobj(d) and not np.isrealobj(s) else float(np.real(d))

        def sines(self):
            """sin of the three principal angles between the regular and the decaying 3-plane at r_m (|D| is their product):
            the singular values of (I - Y Y^T) Z, which is accurate for small angles"""
            Y, Z = self.planes
            R = Z - Y @ (Y.T @ Z)
            return [float(x) for x in np.linalg.svd(R, compute_uv=False)]

        def scalar_plane(self, s):
            w = self.w
            kr = min(np.real(np.sqrt(1 + (s - 1j * w) ** 2)), np.real(np.sqrt(1 + (s + 1j * w) ** 2)))
            Rs = self.rm + 40.0 / kr
            if getattr(self, "Rs_fixed", None) is not None:
                Rs = self.Rs_fixed
            Js = self._Jvac(Rs, s)[2:, 2:]
            lam, vec = np.linalg.eig(Js)
            idx = np.argsort(lam.real)[:2]
            Vinv = np.linalg.inv(vec)
            P = vec[:, idx] @ Vinv[idx, :]
            if not np.iscomplexobj(s):
                P = P.real
            Y = np.zeros((6, 2), dtype=P.dtype); Y[2:, 0] = P[:, 0]; Y[2:, 1] = P[:, 2]
            if self.vac_scalar:
                jf = lambda r, s_: self._Jvac(r, s_)
            else:
                jf = lambda r, s_: self._Jvac(r, s_) if r > self.st.R else self.J(r, s_)
            Q, ph = self._integrate_plane(Y, Rs, self.rm, s, Jfun=jf)
            return Q[2:, :], ph

    return DriftProblem


_cache = {}


def scipy_star():
    if "sci" not in _cache:
        import background as B
        t = time.time()
        _cache["sci"] = B.solve_star(OMEGA, EPS)
        _cache["t_star"] = time.time() - t
    return _cache["sci"]


def crate_j(path):
    if ("j", path) not in _cache:
        _cache[("j", path)] = CrateJ(path)
    return _cache[("j", path)]


def crate_star(path):
    if path not in _cache:
        _cache[path] = CrateStar(path)
    return _cache[path]


# the plants: which star the interior reads, where M_ext comes from, how the exterior scalar plane integrates
PLANTS = {
    "sci":       dict(star="sci",   mext="star",  vac=False),   # PolarProblem as convergence.py builds it
    "crate":     dict(star="crate", mext="star",  vac=False),   # the crate's star and its M_ext in place of scipy's (the issue's planted input)
    "sci_vac":   dict(star="sci",   mext="star",  vac=True),    # scipy's star, the exterior scalar plane in vacuum as the crate integrates it
    "crate_vac": dict(star="crate", mext="star",  vac=True),
    "int_crate": dict(star="crate", mext="sci",   vac=False),   # crate's star in the interior, scipy's M_ext
    "ext_crate": dict(star="sci",   mext="crate", vac=False),   # scipy's star in the interior, the crate's M_ext
    "crate_j":   dict(star="crate", mext="star",  vac=False, j=True),   # the crate's star, M_ext and J tables: J rounded once from double-double
    "crate_j_r13": dict(star="crate", mext="star", vac=False, j=True, rq=(1e-13, 1e-16)),     # as crate_j, the Riccati of q at rtol 1e-13, atol 1e-16 (A5)
    "crate_j_cf025": dict(star="crate", mext="star", vac=False, j=True, cf=0.25),            # as crate_j, 1/4, 1/2, 2, 4, 8 times the renormalizations of the interior plane (A6)
    "crate_j_cf05": dict(star="crate", mext="star", vac=False, j=True, cf=0.5),
    "crate_j_cf2": dict(star="crate", mext="star", vac=False, j=True, cf=2.0),
    "crate_j_cf4": dict(star="crate", mext="star", vac=False, j=True, cf=4.0),
    "crate_j_cf8": dict(star="crate", mext="star", vac=False, j=True, cf=8.0),
    "crate_j_q13": dict(star="crate", mext="star", vac=False, j=True, qeps=1e-13),            # as crate_j, q planted 1e-13 too large (A5)
    "crate_j_q12": dict(star="crate", mext="star", vac=False, j=True, qeps=1e-12),
    "crate_j_vac": dict(star="crate", mext="star", vac=True, j=True),   # as crate_j, the exterior scalar plane in vacuum as the crate integrates it (A4)
}
JTABLE_PATH = [None]


def build(l, variant, plant, rtol, steps_path, atol=1e-14, extra_chunks=0, rm_from=None, noise_ulp=0.0, noise_seed=0, noise_row=0.0):
    """a DriftProblem for the case. rm_from: None (PolarProblem's own grid rule on the plant's star), or 'crate' (the crate's r_m, passed in)"""
    DP = make_problem_class()
    pl = PLANTS[plant]
    star = scipy_star() if pl["star"] == "sci" else crate_star(steps_path)
    kw = dict(rtol=rtol, phi_cut=VARIANTS[variant]["phi_cut"], atol=atol, extra_chunks=extra_chunks, vac_scalar=pl["vac"],
              noise_ulp=noise_ulp, noise_seed=noise_seed, noise_row=noise_row)
    rm = RM[variant] if rm_from == "crate" else None
    if pl.get("j"):
        kw["jtab"] = crate_j(JTABLE_PATH[0])
    if pl.get("rq"):
        kw["riccati_rtol"], kw["riccati_atol"] = pl["rq"]
    if pl.get("qeps"):
        kw["q_eps"] = pl["qeps"]
    if pl.get("cf"):
        kw["chunk_factor"] = pl["cf"]
    pp = DP(star, l, rm=rm, **kw)
    if pl["mext"] != "star":
        other = scipy_star() if pl["mext"] == "sci" else crate_star(steps_path)
        al, a, phi, u = other.bg(np.array([pp.rm]))
        pp.Mext = float(pp.rm * (1 - 1 / a[0] ** 2) / 2)
    return pp


class Stopped(Exception):
    pass


def _alarm(signum, frame):
    raise Stopped()


def run_cases(args):
    """one JSON line per (l, gamma, variant, rtol, plant, atol, extra_chunks) not yet in --out"""
    out = args.out
    done = set()
    if os.path.exists(out):
        for line in open(out):
            if line.startswith("{"):
                o = json.loads(line)
                done.add(case_key(o))
    cases = []
    for rtol in args.rtols:
        for g in args.gammas:
            for v in args.variants:
                for pl in args.plants:
                    for sd in args.noise_seeds:
                        for ex in args.extra_chunks:
                            cases.append(dict(l=args.l, gamma=g, variant=v, rtol=rtol, plant=pl, atol=args.atol, extra_chunks=ex,
                                              rm_from=args.rm_from, noise_ulp=args.noise_ulp, noise_seed=sd, noise_row=args.noise_row))
    signal.signal(signal.SIGALRM, _alarm)
    for c in cases:
        if case_key(c) in done:
            continue
        t0 = time.time()
        rec = dict(c, host=os.uname().nodename, load=os.getloadavg()[0])
        try:
            signal.alarm(int(args.max_seconds))
            pp = build(c["l"], c["variant"], c["plant"], c["rtol"], args.steps, c["atol"], c["extra_chunks"], c["rm_from"], c["noise_ulp"], c["noise_seed"], c["noise_row"])
            t1 = time.time()
            D = pp.D(c["gamma"])
            signal.alarm(0)
            rec.update(status="done", D=float(D), sines=pp.sines(), rm=float(pp.rm), mext=float(pp.Mext), alpha_m=float(pp.alpha_m), lapse_ratio=float(pp.lapse_ratio),
                       wall=time.time() - t1, t_setup=t1 - t0)
        except Stopped:
            rec.update(status="stopped", wall=time.time() - t0, note=f"not finished in {args.max_seconds} s")
        print(json.dumps(rec), flush=True)
        with open(out, "a") as fp:
            fp.write(json.dumps(rec) + "\n")


# ---------------------------------------------------------------------------------------------------------------
# the crate's reference values (taylor_polar.py's `rows` command, #187), through taylor_polar.get_rows
# ---------------------------------------------------------------------------------------------------------------
CRATE_PLAN = [("dd", ["base", "phi", "tol"], ""), ("dd", ["o24", "phio24"], "_o24"), ("f64", ["base", "phi"], "")]


def run_crate(args):
    import taylor_polar as TP
    os.makedirs(args.out_dir, exist_ok=True)
    TP.JDIR = args.out_dir
    plan = list(CRATE_PLAN) + ([("qd", ["base", "phi"], "")] if args.qd else [])
    for lane, variants, tag in plan:
        t = time.time()
        TP.get_rows(args.l, args.qd_gammas if lane == "qd" else args.gammas, lane, variants, tag=tag)
        print(f"crate {lane} {variants}{tag}: {time.time() - t:.1f} s", flush=True)


def load_crate(dirn, l=5):
    """{(lane, variant, gamma): record} from the crate's rows_l{l}_*.jsonl in dirn"""
    out = {}
    for fn in sorted(os.listdir(dirn)):
        if fn.startswith(f"rows_l{l}_") and fn.endswith(".jsonl"):
            for line in open(os.path.join(dirn, fn)):
                if line.startswith("{"):
                    o = json.loads(line)
                    if o["cmd"] == "row":
                        out[(o["lane"], o["variant"], o["gamma"])] = o
    return out


# ---------------------------------------------------------------------------------------------------------------
# a planted system through scipy's integrator: y' = A y, A = taylor_polar.planted_A(gamma), exact D from the invariant subspaces
# ---------------------------------------------------------------------------------------------------------------
def planted_exact(gamma, mp_dps=90):
    """D of the planted system at the converged limit (mpmath, 90 digits; shares no code with scipy's integrator)"""
    import mpmath as mp
    import taylor_polar as TP
    mp.mp.dps = mp_dps
    A = TP.planted_A(gamma)
    QY, sy, _, _ = TP.converged_plane(A, TP.Y0, mp.mpf(24), True)
    QZ, sz, _, _ = TP.converged_plane(A, TP.Z0, -mp.mpf(24), False)
    M = mp.matrix(6, 6)
    for c in range(3):
        for i in range(6):
            M[i, c], M[i, 3 + c] = QY[c][i], QZ[c][i]
    return mp.det(M) * sy * sz


def planted_scipy(gamma, rtol, atol, L, Lz, noise_ulp=0.0, seed=0):
    """scipy's chunked DOP853 plane integration (PolarProblem._integrate_plane, unchanged) on the planted system"""
    import taylor_polar as TP
    DP = make_problem_class()

    class _Star:     # PolarProblem.__init__ wants a star; it is not read by _integrate_plane with a Jfun
        w, eps, r0, R = OMEGA, EPS, 1e-3, 60.0

        def bg(self, r):
            r = np.atleast_1d(np.asarray(r, dtype=float))
            return np.ones_like(r), np.ones_like(r), np.zeros_like(r) + 1e-30, np.zeros_like(r)

    pp = DP(_Star(), 5, rm=L, r0=1e-3, rtol=rtol, atol=atol)
    A = np.array([[float(TP.planted_A(gamma)[i, j]) for j in range(6)] for i in range(6)])
    rng = np.random.default_rng(seed)
    if noise_ulp:
        Jf = lambda r, s: A * (1.0 + noise_ulp * 2.0 ** -53 * rng.standard_normal(A.shape))   # as DriftProblem.J's planted noise
    else:
        Jf = lambda r, s: A
    Y, py = pp._integrate_plane(np.array(TP.Y0, dtype=float).T, 0.0, L, gamma, Jfun=Jf)
    Z, pz = pp._integrate_plane(np.array(TP.Z0, dtype=float).T, L + Lz, L, gamma, Jfun=Jf)
    R = Z - Y @ (Y.T @ Z)
    planted_scipy.sines = [float(x) for x in np.linalg.svd(R, compute_uv=False)]    # the principal angles of the planes scipy produced
    return float(np.linalg.det(np.hstack([Y, Z])) * py * pz)


def run_planted(args):
    import mpmath as mp
    done = set()
    if os.path.exists(args.out):
        for line in open(args.out):
            if line.startswith("{"):
                o = json.loads(line)
                done.add((o["gamma"], o["rtol"], o["atol"], o["L"], o.get("noise_ulp", 0.0), o.get("noise_seed", 0)))
    for rtol in args.rtols:
        for g in args.gammas:
            for sd in args.noise_seeds:
                if (g, rtol, args.atol, args.L, args.noise_ulp, sd) in done:
                    continue
                t = time.time()
                ex = planted_exact(g)
                Lz = 30.0 / g
                D = planted_scipy(g, rtol, args.atol, args.L, Lz, args.noise_ulp, sd)
                rec = dict(gamma=g, rtol=rtol, atol=args.atol, L=args.L, Lz=Lz, D=D, D_exact=mp.nstr(ex, 30), err=float(abs(mp.mpf(D) / ex - 1)),
                           noise_ulp=args.noise_ulp, noise_seed=sd, sines=planted_scipy.sines, wall=time.time() - t, host=os.uname().nodename)
                print(json.dumps(rec), flush=True)
                with open(args.out, "a") as fp:
                    fp.write(json.dumps(rec) + "\n")


# ---------------------------------------------------------------------------------------------------------------
# the declared criteria (negative-mass-raw I196-mac README.md has the same numbers, the reasons, and the computed effect of each mutation)
# ---------------------------------------------------------------------------------------------------------------
CRIT = dict(
    L=5, GAMMAS=[1000.0, 2000.0, 3000.0], RTOLS=[1e-10, 1e-11, 1e-12, 1e-13], VARIANTS=["base", "phi"],
    IDENT_POINT=(2, 10.0),                      # DriftProblem with no knob set is PolarProblem, bit for bit, here
    CJ_RMS_ULP=10.0, CJ_MAX_ULP=100.0,         # the crate's J tables in f64 against the generated J at 60 digits on the exact star, per entry, in ulps
    CS_BAR=1e-14, CS_RADII=[0.002, 0.05, 0.5, 2.0, 5.0, 8.0, 12.0, 16.0, 20.0, 24.10980364018201, 30.0, 34.20514015700785],
    REP_BAR=1e-8,                               # this Mac's scipy D against the win laptop's stored one at (5, 1000)
    REF_BAR=1e-20,                              # the crate's dd D against order 24 at 1e-22 and against quad-double
    PLN0_GAMMAS=[10.0, 100.0], PLN0_BAR=1e-6,   # the planted harness, integrator at rtol 1e-10: gross errors only
    TRUNC_RATIO=0.1, FLOOR_RATIO=0.5,           # gap(rtol 1e-13) / gap(rtol 1e-10)
    PLANT_REMOVES=0.1, PLANT_NONE=0.5,          # gap with the plant / gap without, at the same rtol
    NULL_FACTOR=3.0, MOVE_FLOOR=1e-12,          # a plant moves D if |D_plant / D - 1| > max(3 eta, 1e-12), eta = the null control's change
    PLN_ACCOUNTS=0.3, PLN_NOT=0.03,             # planted-system error / real gap
    NOISE_ACCOUNTS=0.3, NOISE_NOT=0.03,         # planted-noise scatter (1 ulp) / the null scatter eta of the same cell
    NOISE_LINEAR=0.5,                           # scatter at 0.1 ulp <= 0.5 x scatter at 1 ulp
    KINK_REMOVES=0.1, KINK_NONE=0.5,            # eta with the crate's star (N1) / eta with scipy's (N0), in 5 of 6 cells
)
MUTATE = None
RESULTS = []
PLANES_FILES = []     # the planes*.jsonl files of the report's folder (set by the report)


def report(name, ok, text):
    ok = None if ok is None else bool(ok)
    RESULTS.append((name, ok))
    print(f"CHECK {name}: {'SKIPPED' if ok is None else 'PASS' if ok else 'FAIL'}  {text}", flush=True)


def fmt(x):
    return f"{float(x):.2e}"


def classify_rtol(g10, g13):
    """one cell: the gap at rtol 1e-10 and 1e-13"""
    r = g13 / g10
    return "truncation" if r <= CRIT["TRUNC_RATIO"] else ("floor" if r >= CRIT["FLOOR_RATIO"] else "mixed")


def classify_plant(g_without, g_with):
    r = g_with / g_without
    return "removes" if r <= CRIT["PLANT_REMOVES"] else ("none" if r >= CRIT["PLANT_NONE"] else "partly")


def classify_pln(err, gap):
    r = err / gap
    return "accounts" if r >= CRIT["PLN_ACCOUNTS"] else ("not" if r <= CRIT["PLN_NOT"] else "partly")


def load_scipy(d):
    import glob
    out = {}
    for fn in sorted(glob.glob(os.path.join(d, "scipy*.jsonl"))):
        for line in open(fn):
            if line.startswith("{"):
                o = json.loads(line)
                if o.get("status") == "done":
                    out[case_key(o)] = o
    return out


def sget(sc, g, v, rtol, plant="sci", atol=1e-14, extra=0, rm_from=None, l=5, noise=0.0, seed=0, row=0.0):
    return sc.get((l, g, v, rtol, plant, atol, extra, rm_from, noise, seed, row))


def crate_ref(crate, lane, v, g):
    import mpmath as mp
    mp.mp.dps = 60
    return mp.mpf(crate[(lane, v, g)]["d"])


def gap(D, ref):
    import mpmath as mp
    mp.mp.dps = 60
    return float(abs(mp.mpf(D) / ref - 1))


STORED = {"base": ("base (rtol1e-10, phi_cut1e-9, R60)", 1e-10), "tol": ("rtol 1e-12", 1e-12), "phi": ("phi_cut 1e-12 (larger r_m)", 1e-10)}


def stored_win(raw, l=5, g=1000.0):
    """win laptop's scipy D at (l, g): {name: D} from I4-win's partial file"""
    out = {}
    for line in open(os.path.join(raw, "I4-win", f"convergence_l{l}_w0.7603.partial.csv")):
        if line.startswith("["):
            k, gg, v, _ = json.loads(line)
            if gg == g and isinstance(v, float):
                out[k] = v
    return out


def check_IDENT(steps):
    from polar_solver import PolarProblem
    l, g = CRIT["IDENT_POINT"]
    DP = make_problem_class()
    st = scipy_star()
    d0 = PolarProblem(st, l).D(g)
    kw = dict(atol=1e-6) if MUTATE == "ident_atol" else {}
    d1 = DP(st, l, **kw).D(g)
    d2 = DP(st, l, atol=1e-6).D(g)
    d3 = DP(st, l, extra_chunks=1).D(g)
    d4 = DP(st, l, noise_ulp=1.0, noise_seed=1).D(g)
    d5 = DP(st, l, noise_ulp=1.0, noise_seed=1).D(g)
    d6 = DP(st, l, noise_row=1.0, noise_seed=1).D(g)
    d7 = DP(st, l, noise_row=1.0, noise_seed=1).D(g)
    d8 = DP(st, l, q_eps=1e-9).D(g)
    d9 = DP(st, l, riccati_rtol=1e-5).D(g)
    ok = (d0 == d1) and (d2 != d0) and (d3 != d0) and (d4 != d0) and (d4 == d5) and (d6 != d0) and (d6 == d7) and (d8 != d0) and (d9 != d0)
    report("IDENT", ok, f"(l, gamma) = ({l}, {g:g}): PolarProblem D {d0!r}, DriftProblem with no knob {d1!r} (equal: {d0 == d1}); atol 1e-6 moves D by "
           f"{fmt(abs(d2 / d0 - 1))}, one extra chunk by {fmt(abs(d3 / d0 - 1))}, 1 ulp of planted noise by {fmt(abs(d4 / d0 - 1))}, "
           f"1 ulp of the row maximum by {fmt(abs(d6 / d0 - 1))}, q planted 1e-9 off by {fmt(abs(d8 / d0 - 1))}, the Riccati at rtol 1e-5 by {fmt(abs(d9 / d0 - 1))} (each must move it), the same seed again gives the same D: {d4 == d5 and d6 == d7}")


def check_CS(steps):
    import mpmath as mp
    import taylor_polar as TP
    mp.mp.dps = 40
    cs = crate_star(steps)
    if MUTATE == "cs_order":
        cs = CrateStar(steps)
        cs.k = np.arange(NS + 1)
        cs.C = cs.C.copy()
        cs.C[:, :, 9:] = 0.0          # the series cut at order 8
    radii = CRIT["CS_RADII"]
    sv = TP.star_values(TP.star_kw(), radii)
    worst, at = 0.0, None
    for r in radii:
        o = sv[r]
        al, a, phi, u = (x[0] for x in cs.bg(np.array([r])))
        for name, val in (("alpha", al), ("a", a), ("phi0", phi), ("u", u)):
            e = float(abs(mp.mpf(float(val)) / mp.mpf(o[name]) - 1))
            if e > worst:
                worst, at = e, (r, name)
    report("CS", worst <= CRIT["CS_BAR"], f"the crate's star from its dumped Taylor steps in f64 against `starcheck` (double-double) at {len(radii)} radii, 4 functions: "
           f"worst relative {fmt(worst)} at (r, f) = {at} (bar {CRIT['CS_BAR']:.0e})")


def check_CJ(steps, jt):
    """the crate's J tables (CrateJ) against the generated J_and_h at 60 digits on the crate's star at 12 radii and gamma = 1000, 3000, per entry
    in units of 2^-53; J_and_h in f64 (scipy's) against the same at its own f64 inputs is printed beside it"""
    import mpmath as mp
    import taylor_polar as TP
    import background as B
    from polar_matrix_generated import J_and_h
    mp.mp.dps = 60
    f = TP.mp_J()
    cj = CrateJ(jt)
    cj.shift_end = MUTATE == "cj_shift"
    radii = CRIT["CS_RADII"]
    sv = TP.star_values(TP.star_kw(), radii)
    cs = crate_star(steps)
    l, LL = 5, 30
    e_tab, e_gen = [], []
    for g in (1000.0, 3000.0):
        for r in radii:
            o = sv[r]
            Jm, _ = f(mp.mpf(r), mp.mpf(o["alpha"]), mp.mpf(o["a"]), mp.mpf(o["phi0"]), mp.mpf(o["u"]), mp.mpf(OMEGA), mp.mpf(EPS), mp.mpf(g), mp.mpf(LL), mp.mpf(B.LAM), mp.mpf(B.MU))
            Jt = cj.J(r, g)
            al, a, ph, u = (float(x[0]) for x in cs.bg(np.array([r])))
            Jg, _ = J_and_h(r, al, a, ph, u, OMEGA, EPS, g, LL, B.LAM, B.MU)
            Jm2, _ = f(mp.mpf(r), mp.mpf(al), mp.mpf(a), mp.mpf(ph), mp.mpf(u), mp.mpf(OMEGA), mp.mpf(EPS), mp.mpf(g), mp.mpf(LL), mp.mpf(B.LAM), mp.mpf(B.MU))
            for i in range(6):
                for j in range(6):
                    if Jm[i, j] != 0:
                        e_tab.append(float(abs(mp.mpf(float(Jt[i, j])) - Jm[i, j]) / abs(Jm[i, j])) / 2.0 ** -53)
                    if Jm2[i, j] != 0:
                        e_gen.append(float(abs(mp.mpf(float(Jg[i, j])) - Jm2[i, j]) / abs(Jm2[i, j])) / 2.0 ** -53)
    e_tab, e_gen = np.array(e_tab), np.array(e_gen)
    rms, mx = float(np.sqrt(np.mean(e_tab ** 2))), float(e_tab.max())
    report("CJ", rms <= CRIT["CJ_RMS_ULP"] and mx <= CRIT["CJ_MAX_ULP"],
           f"the crate's J tables in f64 against the generated J at 60 digits on the exact star, {e_tab.size} entries (12 radii, gamma 1000 and 3000): "
           f"rms {rms:.2f} ulp, max {mx:.1f} ulp (bars {CRIT['CJ_RMS_ULP']:g}, {CRIT['CJ_MAX_ULP']:g}); the generated J in f64 on its own inputs, same entries: "
           f"rms {float(np.sqrt(np.mean(e_gen ** 2))):.0f} ulp, max {float(e_gen.max()):.0f} ulp")


def check_RMX(steps):
    from polar_solver import PolarProblem
    st, cs = scipy_star(), crate_star(steps)
    rows, ok = [], True
    for v in CRIT["VARIANTS"]:
        cut = 1e-11 if (MUTATE == "rm_cut" and v == "phi") else VARIANTS[v]["phi_cut"]
        r_s = PolarProblem(st, 2, phi_cut=VARIANTS[v]["phi_cut"]).rm
        r_c = PolarProblem(cs, 2, phi_cut=cut).rm
        ok &= (r_s == RM[v]) and (r_c == RM[v])
        rows.append(f"{v}: scipy's star {r_s!r}, the crate's star {r_c!r}, taylor_polar R_M {RM[v]!r}")
    report("RMX", ok, "PolarProblem's r_m on scipy's star and on the crate's star equal to the last bit: " + "; ".join(rows))


def check_REP(sc, raw):
    stored = stored_win(raw)
    rows, worst = [], 0.0
    for key, (name, rtol) in STORED.items():
        v = "phi" if key == "phi" else "base"
        rec = sget(sc, 1000.0, v, rtol)
        if rec is None:
            rows.append(f"{name}: not run")
            continue
        ref = stored[STORED["phi" if (MUTATE == "rep_variant" and key == "base") else key][0]]
        e = abs(rec["D"] / ref - 1)
        worst = max(worst, e)
        rows.append(f"{name}: here {rec['D']:.15e}, win laptop {ref:.15e}, {fmt(e)}")
    if not any("here" in r for r in rows):
        report("REP", None, "no (5, 1000) case present")
        return
    report("REP", worst <= CRIT["REP_BAR"], f"scipy on this Mac against the win laptop's stored scipy D at (5, 1000): worst {fmt(worst)} (bar {CRIT['REP_BAR']:.0e}); " + "; ".join(rows))


def check_REF(crate):
    import mpmath as mp
    mp.mp.dps = 60
    rows, worst, ok = [], 0.0, True
    for g in CRIT["GAMMAS"]:
        for v in CRIT["VARIANTS"]:
            d = crate_ref(crate, "dd", v, g)
            if MUTATE == "ref_tol" and v == "base":
                d = crate_ref(crate, "dd", "tol", g)
            alt = [("o24", crate_ref(crate, "dd", "o24" if v == "base" else "phio24", g))]
            if ("qd", v, g) in crate:
                alt.append(("qd", crate_ref(crate, "qd", v, g)))
            for n, a in alt:
                e = float(abs(d / a - 1))
                worst = max(worst, e)
                ok &= e <= CRIT["REF_BAR"]
                rows.append(f"({g:g}, {v}) dd vs {n} {fmt(e)}")
    report("REF", ok, f"the crate's dd D against order 24 at tol 1e-22 (and quad-double where run): worst {fmt(worst)} (bar {CRIT['REF_BAR']:.0e}); " + "; ".join(rows))


def name_cause(st):
    """st: {'H1': ..., 'H2': ..., 'H3': ..., 'H4': ...} each 'supported', 'ruled out' or 'undecided'. A cause is named only if exactly
    one hypothesis is supported and the others are ruled out; otherwise 'not established'"""
    sup = [h for h, v in st.items() if v == "supported"]
    if len(sup) == 1 and all(v == "ruled out" for h, v in st.items() if h != sup[0]):
        return sup[0]
    return "not established"


def cje_stats(sc, crate, g, v, rtol, atol):
    """signed gaps D / D_crate - 1 of the CJE_N realizations (extra chunks 0 .. CJE_N - 1) of one cell and setting: (list, mean, std ddof 1)"""
    import mpmath as mp
    mp.mp.dps = 60
    ref = crate_ref(crate, "dd", v, g)
    xs = [float(mp.mpf(sget(sc, g, v, rtol, "crate_j", atol, extra=e)["D"]) / ref - 1) for e in range(CJE_N)]
    return xs, float(np.mean(xs)), float(np.std(xs, ddof=1))


def verdicts(sc, crate, plc):
    """the five hypotheses against the numbers, by the declared rules (README, A3 for H1, H3, H5). Returns (statuses, lines)"""
    G, V = CRIT["GAMMAS"], CRIT["VARIANTS"]
    ref = lambda v, g: crate_ref(crate, "dd", v, g)
    cells = [(g, v) for g in G for v in V]
    lines, st = [], {}
    eta0 = {c: abs(sget(sc, c[0], c[1], 1e-10, extra=1)["D"] / sget(sc, c[0], c[1], 1e-10)["D"] - 1) for c in cells}
    eta1 = {c: abs(sget(sc, c[0], c[1], 1e-10, "crate", extra=1)["D"] / sget(sc, c[0], c[1], 1e-10, "crate")["D"] - 1) for c in cells}
    # H1 (A3): truncation controlled by rtol, in the system whose J is correctly rounded (CJE): sigma(rtol 1e-13) / sigma(rtol 1e-10) at atol 1e-14
    try:
        sg = {c: [cje_stats(sc, crate, c[0], c[1], r, a)[2] for (r, a) in CJE_SETTINGS] for c in CJ_CELLS}
        r1 = [sg[c][1] / sg[c][0] for c in CJ_CELLS]
        r5 = [sg[c][2] / sg[c][1] for c in CJ_CELLS]
        st["H1"] = "supported" if all(x <= CRIT["TRUNC_RATIO"] for x in r1) else ("ruled out" if all(x >= CRIT["FLOOR_RATIO"] for x in r1) else "undecided")
        st["H5"] = "supported" if all(x <= CRIT["TRUNC_RATIO"] for x in r5) else ("ruled out" if all(x >= CRIT["FLOOR_RATIO"] for x in r5) else "undecided")
    except TypeError:
        r1 = r5 = "not run"
        st["H1"] = st["H5"] = "undecided"
    f2 = lambda x: ["%.2f" % y for y in x] if x != "not run" else x
    lines.append(f"H1 truncation by rtol (exact-J system, 8 realizations per cell): std(rtol 1e-13) / std(rtol 1e-10) per cell {f2(r1)}; {st['H1']}")
    lines.append(f"H5 truncation by atol: std(atol 1e-17) / std(atol 1e-14) at rtol 1e-13 per cell {f2(r5)}; {st['H5']}")
    # H2: the crate's star and M_ext move D beyond 3 eta, and remove the gap, in the issue's grid (rtol 1e-10)
    mv = [abs(sget(sc, g, v, 1e-10, "crate")["D"] / sget(sc, g, v, 1e-10)["D"] - 1) > max(CRIT["NULL_FACTOR"] * eta0[(g, v)], CRIT["MOVE_FLOOR"]) for g, v in cells]
    pl = [classify_plant(gap(sget(sc, g, v, 1e-10)["D"], ref(v, g)), gap(sget(sc, g, v, 1e-10, "crate")["D"], ref(v, g))) for g, v in cells]
    st["H2"] = "supported" if all(mv) and all(x == "removes" for x in pl) else ("ruled out" if sum(mv) <= 1 else "undecided")
    lines.append(f"H2 the star and M_ext as inputs: the crate's star moves D beyond 3 eta in {sum(mv)} of 6 cells; the gap at rtol 1e-10 with the plant {pl}; {st['H2']}")
    # H3 (A3): rounding in the generated J_and_h: the crate's J tables against J_and_h on the crate's star, the realization scatter eta of the same cells
    try:
        etaj = {c: abs(sget(sc, c[0], c[1], 1e-10, "crate_j", extra=1)["D"] / sget(sc, c[0], c[1], 1e-10, "crate_j")["D"] - 1) for c in CJ_CELLS}
        rj = [etaj[c] / eta1[c] for c in CJ_CELLS]
        st["H3"] = "supported" if all(x <= CRIT["KINK_REMOVES"] for x in rj) else ("ruled out" if all(x >= CRIT["KINK_NONE"] for x in rj) else "undecided")
    except TypeError:
        rj, st["H3"] = "not run", "undecided"
    lines.append(f"H3 rounding in the generated J_and_h: eta(crate J tables) / eta(J_and_h), both on the crate's star (cells 1000 base, 1000 phi, 3000 base) {f2(rj)}; {st['H3']}")
    # H6 (A4): the exterior scalar plane through the star's J (couplings proportional to phi0 kept) against vacuum, at (3000, base), rtol 1e-13, atol 1e-14
    try:
        import mpmath as mp3
        mp3.mp.dps = 60
        rf = crate_ref(crate, "dd", "base", 3000.0)
        mv_ = float(np.mean([float(mp3.mpf(sget(sc, 3000.0, "base", 1e-13, "crate_j_vac", extra=e)["D"]) / rf - 1) for e in range(CJE_N)]))
        mn_ = cje_stats(sc, crate, 3000.0, "base", 1e-13, 1e-14)[1]
        q6 = abs(mv_) / abs(mn_)
        st["H6"] = "supported" if q6 <= 0.3 else ("ruled out" if q6 >= 0.7 else "undecided")
        lines.append(f"H6 the exterior scalar plane through the star's J against vacuum, (3000, base) rtol 1e-13: mean signed gap {mn_:+.2e} through the star's J, {mv_:+.2e} in vacuum, ratio {q6:.2f}; {st['H6']}")
    except TypeError:
        st["H6"] = "undecided"
        lines.append("H6 the exterior scalar plane through the star's J against vacuum: not run; undecided")
    # H7 (A5, rebased in A6 on CJQ): the Riccati integration of the exterior q (error 6e-14 relative, read from the stored files) times D's measured sensitivity to q
    # (planted q errors 1e-13 and 1e-12, (3000, base)) against the bias of CJE at (3000, base), rtol 1e-13, atol 1e-14
    try:
        d0q = sget(sc, 3000.0, "base", 1e-10, "crate_j", extra=0)["D"]
        sens = max(abs(sget(sc, 3000.0, "base", 1e-10, pl)["D"] / d0q - 1) / eps for pl, eps in (("crate_j_q13", 1e-13), ("crate_j_q12", 1e-12)))
        bias = abs(cje_stats(sc, crate, 3000.0, "base", 1e-13, 1e-14)[1])
        t7 = sens * 6e-14 / bias
        st["H7"] = "supported" if t7 >= 0.3 else ("ruled out" if t7 <= 0.03 else "undecided")
        lines.append(f"H7 the Riccati of the exterior q (hard-coded rtol 1e-11): D's sensitivity to q's relative error {sens:.3g}, times q's error 6e-14, over CJE's bias at (3000, base) {bias:.2e}: {t7:.1e}; {st['H7']}")
    except (TypeError, ZeroDivisionError):
        st["H7"] = "undecided"
        lines.append("H7 the Riccati of the exterior q: not run; undecided")
    # H8 (A6): roundoff accumulation of the interior plane's renormalizations: median eps_Y against the number of renormalizations (factor 1/4 .. 8) at (1000, base), rtol 1e-13
    try:
        med = {}
        for fn in PLANES_FILES:
            for x in open(fn):
                if x.startswith("{"):
                    o = json.loads(x)
                    if o["gamma"] == 1000.0 and o["variant"] == "base" and o["rtol"] == 1e-13:
                        f = {"crate_j": 1.0, "crate_j_cf025": 0.25, "crate_j_cf05": 0.5, "crate_j_cf2": 2.0, "crate_j_cf4": 4.0, "crate_j_cf8": 8.0}.get(o["plant"])
                        if f is not None:
                            med.setdefault(f, []).append(o["eps_y"])
        m = {f: float(np.median(v)) for f, v in med.items()}
        up = m[8.0] / m[1.0]
        down = m[0.25] / m[1.0]
        flat = max(m[0.25], m[1.0], m[8.0]) / min(m[0.25], m[1.0], m[8.0]) < 1.3
        st["H8"] = "supported" if up >= 2.0 and down <= 0.7 else ("ruled out" if flat else "undecided")
        lines.append(f"H8 roundoff of the interior renormalizations: median eps_Y at 8x the renormalizations over 1x {up:.2f}, at 1/4x over 1x {down:.2f}; {st['H8']}")
    except (KeyError, TypeError, FileNotFoundError):
        st["H8"] = "undecided"
        lines.append("H8 roundoff of the interior renormalizations: not run; undecided")
    # H4: the BVP star's kinks make the integrator's error erratic: eta falls when the star is smooth
    rr = [eta1[c] / eta0[c] for c in cells]
    nlow = sum(1 for x in rr if x <= CRIT["KINK_REMOVES"])
    nhigh = sum(1 for x in rr if x >= CRIT["KINK_NONE"])
    st["H4"] = "supported" if nlow >= 5 else ("ruled out" if nhigh >= 5 else "undecided")
    lines.append(f"H4 the BVP star's kinks: eta(crate star) / eta(scipy star) per cell {f2(rr)}; {st['H4']}")
    return st, lines


def check_CLS():
    syn_trunc = (1e-8, 1e-11)      # gap falls with rtol
    syn_floor = (1e-8, 0.9e-8)
    syn_mixed = (1e-8, 2e-9)
    a, b, c = (classify_rtol(*x) for x in (syn_trunc, syn_floor, syn_mixed))
    if MUTATE == "cls_swap":
        a, b = b, a
    p = (classify_plant(1e-8, 1e-10), classify_plant(1e-8, 0.9e-8), classify_plant(1e-8, 2e-9))
    q = (classify_pln(1e-8, 1e-8), classify_pln(1e-11, 1e-8), classify_pln(1e-9, 1e-8))
    S, Rr, U = "supported", "ruled out", "undecided"
    nc = (name_cause(dict(H1=S, H2=Rr, H3=Rr, H4=Rr)), name_cause(dict(H1=U, H2=Rr, H3=Rr, H4=Rr)), name_cause(dict(H1=S, H2=S, H3=Rr, H4=Rr)),
          name_cause(dict(H1=Rr, H2=Rr, H3=S, H4=Rr)))
    if MUTATE == "cls_swap":
        nc = (nc[1],) + nc[1:]
    ok = ((a, b, c) == ("truncation", "floor", "mixed") and p == ("removes", "none", "partly") and q == ("accounts", "not", "partly")
          and nc == ("H1", "not established", "not established", "H3"))
    report("CLS", ok, f"the classifiers and the naming rule on synthetic tables with known labels: rtol curve {(a, b, c)}, plant {p}, planted system {q}, cause {nc}")


def check_PLN0():
    import mpmath as mp
    worst = 0.0
    for g in CRIT["PLN0_GAMMAS"]:
        ex = planted_exact(g)
        if MUTATE == "pln_sign":
            ex = -ex
        D = planted_scipy(g, 1e-10, 1e-14, 24.1, 30.0 / g)
        worst = max(worst, float(abs(mp.mpf(D) / ex - 1)))
    report("PLN0", worst <= CRIT["PLN0_BAR"], f"scipy's integrator on the planted system (taylor_polar's P1, L = 24.1) at gamma = {CRIT['PLN0_GAMMAS']} against the invariant subspaces at 90 digits: "
           f"worst {fmt(worst)} (bar {CRIT['PLN0_BAR']:.0e})")


def check_COMPLETE(sc, crate, need):
    have = dict(sc)
    if MUTATE == "drop_cell":
        have.pop(next(iter(need)), None)
    miss = [k for k in need if k not in have]
    report("COMPLETE", not miss, f"{len(need) - len(miss)} of {len(need)} declared scipy cases present and finished" + (f"; missing {miss[:6]}" if miss else ""))


# ---------------------------------------------------------------------------------------------------------------
# the declared grid and the checks and tables over it
# ---------------------------------------------------------------------------------------------------------------
def K(g, v, rtol, plant="sci", atol=1e-14, extra=0, noise=0.0, seed=0, row=0.0):
    """a scipy case's key (case_key's tuple), l = 5, r_m by PolarProblem's own rule"""
    return (5, g, v, rtol, plant, atol, extra, None, noise, seed, row)


CJ_CELLS = [(1000.0, "base"), (1000.0, "phi"), (3000.0, "base")]
CJE_SETTINGS = [(1e-10, 1e-14), (1e-13, 1e-14), (1e-13, 1e-17)]     # (rtol, atol)
CJE_N = 8                                                           # realizations: extra_chunks 0 .. 7
CJN_NOISE = [("noise", u) for u in (1.0, 10.0, 100.0)] + [("row", u) for u in (0.01, 0.1)]


def grid_cases():
    """every scipy case the issue and this note's additions declare: {name: [case_key, ...]}. CORE is what COMPLETE requires; the others
    (X...) are run after it, in the order of this table, and reported whatever they show. A3 moved NOISE and NOISEROW out of the core
    (not run) and added CJE and CJN."""
    G, V, R = CRIT["GAMMAS"], CRIT["VARIANTS"], CRIT["RTOLS"]
    out = {}
    out["P0"] = [K(g, v, r) for r in R for g in G for v in V]                                    # the issue: scipy's star, 4 rtol, 3 gamma, 2 variants
    out["P1"] = [K(g, v, 1e-10, "crate") for g in G for v in V]                                  # the issue: the crate's star and M_ext (A3: rtol 1e-10 only in the core)
    out["N0"] = [K(g, v, 1e-10, extra=1) for g in G for v in V]                                  # null control: one extra chunk
    out["N1"] = [K(g, v, 1e-10, "crate", extra=1) for g in G for v in V]
    out["CJE"] = [K(g, v, r, "crate_j", a, extra=e) for (g, v) in CJ_CELLS for (r, a) in CJE_SETTINGS for e in range(CJE_N)]      # A3
    out["CJN"] = [K(1000.0, "base", 1e-10, "crate_j", noise=u, seed=sd) for (kind, u) in CJN_NOISE if kind == "noise" for sd in (1, 2, 3, 4)]   # A3, A6: per entry only
    out["X_CJN_ROW"] = [K(1000.0, "base", 1e-10, "crate_j", row=u, seed=sd) for (kind, u) in CJN_NOISE if kind == "row" for sd in (1, 2, 3, 4)]   # A6: 0.01 ulp took 500 s a case, not core
    out["CJV"] = [K(g, v, 1e-13, "crate_j_vac", extra=e) for (g, v) in CJ_CELLS for e in range(CJE_N)]                                      # A4
    out["X_CJR"] = [K(g, v, 1e-13, "crate_j_r13", extra=e) for (g, v) in CJ_CELLS for e in range(CJE_N)]                                     # A5, A6: not core
    out["CJQ"] = [K(g, "base", 1e-10, pl) for g in (1000.0, 3000.0) for pl in ("crate_j_q13", "crate_j_q12")]                                 # A5
    out["X_P1_13"] = [K(g, v, 1e-13, "crate") for g in G for v in V]
    out["X_ATOL"] = [K(1000.0, v, 1e-13, atol=1e-17) for v in V] + [K(3000.0, "base", 1e-13, atol=1e-17)]
    out["X_N0R13"] = [K(1000.0, v, 1e-13, extra=1) for v in V] + [K(3000.0, "base", 1e-13, extra=1)]
    out["X_P2"] = [K(g, v, 1e-10, "sci_vac") for g in (1000.0, 3000.0) for v in V]
    out["X_P3"] = [K(3000.0, "base", 1e-10, pl) for pl in ("int_crate", "ext_crate")]
    out["X_P1R"] = [K(g, v, r, "crate") for r in (1e-11, 1e-12) for g in G for v in V]
    out["X_NOISE"] = ([K(1000.0, "base", 1e-10, noise=1.0, seed=s) for s in (1, 2, 3)] + [K(1000.0, "base", 1e-10, noise=0.1, seed=s) for s in (1, 2)]
                      + [K(3000.0, "base", 1e-10, noise=1.0, seed=s) for s in (1, 2)])
    return out


CORE = ["P0", "P1", "N0", "N1", "CJE", "CJN", "CJV", "CJQ"]


def run_checks(args):
    global MUTATE
    MUTATE = args.mutate
    only = set(args.only.split(",")) if args.only else None
    d = args.dir
    steps = args.steps
    sc = load_scipy(d)
    crate = load_crate(args.crate_dir or d)
    want = lambda n: only is None or n in only
    if want("IDENT"):
        check_IDENT(steps)
    if want("CS"):
        check_CS(steps)
    if want("RMX"):
        check_RMX(steps)
    if want("CJ"):
        check_CJ(steps, args.jtable)
    if want("REP"):
        check_REP(sc, args.raw)
    if want("REF"):
        check_REF(crate)
    if want("CLS"):
        check_CLS()
    if want("PLN0"):
        check_PLN0()
    if want("COMPLETE"):
        need = [k for n, cases in grid_cases().items() if n in CORE for k in cases]
        if args.synthetic:      # a complete set of placeholder records: the check on data it must pass, and fail under drop_cell
            sc = {k: dict(status="done") for k in need}
        check_COMPLETE(sc, crate, need)
    npass = sum(1 for _, ok in RESULTS if ok is True)
    nfail = sum(1 for _, ok in RESULTS if ok is False)
    nskip = [n for n, ok in RESULTS if ok is None]
    print(f"SUMMARY: {npass} passed, {nfail} failed, {len(nskip)} skipped{(' (' + ', '.join(nskip) + ')') if nskip else ''} of {len(RESULTS)}"
          + (f"  [mutation {MUTATE}]" if MUTATE else ""))
    sys.exit(1 if nfail else 0)


def run_report(args):
    import mpmath as mp
    mp.mp.dps = 60
    sc = load_scipy(args.dir)
    crate = load_crate(args.crate_dir or args.dir)
    G, V, R = CRIT["GAMMAS"], CRIT["VARIANTS"], CRIT["RTOLS"]
    ref = lambda v, g: crate_ref(crate, "dd", v, g)
    P = lambda *a, **k: print(*a, **k)

    P("# crate reference D, l = 5 (taylor_polar `rows`, double-double, order 30, tol 1.2e-30), and its own spread")
    P(f"{'gamma':>6s} {'variant':>7s} {'D (dd)':>22s} {'dd/o24-1':>10s} {'dd/qd-1':>10s} {'f64/dd-1':>10s} {'tol1e-12/dd-1':>14s}")
    for g in G:
        for v in V:
            d = ref(v, g)
            o = crate_ref(crate, "dd", "o24" if v == "base" else "phio24", g) if ("dd", "o24" if v == "base" else "phio24", g) in crate else None
            q = crate_ref(crate, "qd", v, g) if ("qd", v, g) in crate else None
            f = crate_ref(crate, "f64", v, g) if ("f64", v, g) in crate else None
            t = crate_ref(crate, "dd", "tol", g) if v == "base" and ("dd", "tol", g) in crate else None
            fm = lambda x: f"{float(abs(x / d - 1)):10.1e}" if x is not None else f"{'-':>10s}"
            P(f"{g:6g} {v:>7s} {mp.nstr(d, 17):>22s} {fm(o)} {fm(q)} {fm(f)} {fm(t) if t is not None else '-':>14s}")

    def block(title, plant, rtols, extra=0, atol=1e-14):
        P(f"\n# {title}")
        P(f"{'gamma':>6s} {'variant':>7s} " + " ".join(f"{'rtol %g' % r:>22s}" for r in rtols))
        for g in G:
            for v in V:
                cells = []
                for r in rtols:
                    rec = sget(sc, g, v, r, plant, atol, extra)
                    cells.append(f"{rec['D']:+.12e} {gap(rec['D'], ref(v, g)):9.1e}"[:22] if rec else f"{'-':>22s}")
                row = []
                for r in rtols:
                    rec = sget(sc, g, v, r, plant, atol, extra)
                    row.append(f"{gap(rec['D'], ref(v, g)):22.2e}" if rec else f"{'-':>22s}")
                P(f"{g:6g} {v:>7s} " + " ".join(row))

    block("P0: scipy's gap |D_scipy / D_crate(dd) - 1| against rtol (scipy's star, atol 1e-14, one realization per cell)", "sci", R)
    block("P1: the same with the crate's star and its M_ext in PolarProblem (r_m the same grid point)", "crate", [1e-10, 1e-13])
    block("NULL: P0 with one extra chunk (the same problem on another chunk grid), rtol 1e-10", "sci", [1e-10], extra=1)
    block("N1: P1 with one extra chunk, rtol 1e-10", "crate", [1e-10], extra=1)
    block("CJ: the crate's star, M_ext and J tables (J rounded once from double-double), rtol 1e-10", "crate_j", [1e-10])
    block("CJ null: the same with one extra chunk", "crate_j", [1e-10], extra=1)
    block("P2: scipy's star, the exterior scalar plane in vacuum as the crate integrates it, rtol 1e-10", "sci_vac", [1e-10])

    P("\n# every realization per cell, signed D/D_crate - 1 (atol 1e-14, no planted noise): scipy's star (P0 at each rtol, N0, X_N0R13) and the crate's star (P1, N1, X_P1R)."
      "\n# n, mean, standard deviation (ddof 1), and mean / (std / sqrt n): a mean beyond about 2 says the realizations share an offset")
    for plant, title in (("sci", "scipy's star"), ("crate", "the crate's star")):
        for g in G:
            for v in V:
                xs = []
                for k, r in sorted(sc.items(), key=lambda kv: (kv[0][3], kv[0][6])):
                    if k[1] == g and k[2] == v and k[4] == plant and k[5] == 1e-14 and k[8] == 0.0 and k[10] == 0.0:
                        xs.append(float(mp.mpf(r["D"]) / ref(v, g) - 1))
                if len(xs) >= 2:
                    m, sd = float(np.mean(xs)), float(np.std(xs, ddof=1))
                    P(f"  {title:>16s} gamma {g:g} {v:>4s}: n {len(xs)}  mean {m:+.2e}  std {sd:.2e}  mean/(std/sqrt n) {m / (sd / np.sqrt(len(xs))):+.1f}   values " + " ".join(f"{x:+.1e}" for x in xs))

    P("\n# the rtol curve per cell: gap(1e-13) / gap(1e-10) and its label (truncation <= %g, floor >= %g)" % (CRIT["TRUNC_RATIO"], CRIT["FLOOR_RATIO"]))
    labels = []
    for g in G:
        for v in V:
            a, b = sget(sc, g, v, 1e-10), sget(sc, g, v, 1e-13)
            if a and b:
                g10, g13 = gap(a["D"], ref(v, g)), gap(b["D"], ref(v, g))
                lab = classify_rtol(g10, g13)
                labels.append(lab)
                P(f"  gamma {g:g} {v}: {g13:.1e} / {g10:.1e} = {g13 / g10:.2f}  {lab}")
    P(f"  cells: {len(labels)}; truncation {labels.count('truncation')}, floor {labels.count('floor')}, mixed {labels.count('mixed')}")

    P("\n# realization noise: |D(extra chunk) / D - 1| at rtol 1e-10, against the gap")
    eta = {}
    for g in G:
        for v in V:
            a, n = sget(sc, g, v, 1e-10), sget(sc, g, v, 1e-10, extra=1)
            if a and n:
                eta[(g, v)] = abs(n["D"] / a["D"] - 1)
                P(f"  gamma {g:g} {v}: eta {eta[(g, v)]:.1e}, gap {gap(a['D'], ref(v, g)):.1e}, gap with the extra chunk {gap(n['D'], ref(v, g)):.1e}")

    P("\n# planted inputs: D_plant / D_scipy - 1 against 3 eta (and 1e-12), and the gap with the plant against the gap without")
    for plant, rtols in (("crate", (1e-10, 1e-13)), ("crate_j", (1e-10,)), ("sci_vac", (1e-10,)), ("int_crate", (1e-10,)), ("ext_crate", (1e-10,))):
        for r in rtols:
            for g in G:
                for v in V:
                    a, b = sget(sc, g, v, r), sget(sc, g, v, r, plant)
                    if a and b:
                        mv = abs(b["D"] / a["D"] - 1)
                        e = eta.get((g, v))
                        moves = (mv > max(CRIT["NULL_FACTOR"] * e, CRIT["MOVE_FLOOR"])) if e is not None else None
                        g0, g1 = gap(a["D"], ref(v, g)), gap(b["D"], ref(v, g))
                        P(f"  {plant:>9s} rtol {r:g} gamma {g:g} {v}: D moves {mv:.1e} (3 eta {CRIT['NULL_FACTOR'] * e if e is not None else float('nan'):.1e}; moves D: {moves}); "
                          f"gap {g0:.1e} -> {g1:.1e} ({g1 / g0:.2f}, {classify_plant(g0, g1)})")

    P("\n# conditioning: sines of the principal angles between the regular and the decaying plane (scipy, P0 rtol 1e-10), |D| = product; "
      "plane error = gap * min sine")
    for g in G:
        for v in V:
            a = sget(sc, g, v, 1e-10)
            if a and "sines" in a:
                sn = a["sines"]
                P(f"  gamma {g:g} {v}: sines {sn[0]:.2e} {sn[1]:.2e} {sn[2]:.2e}, product/|D| {sn[0] * sn[1] * sn[2] / abs(a['D']):.6f}, "
                  f"gap {gap(a['D'], ref(v, g)):.1e}, gap * min sine {gap(a['D'], ref(v, g)) * min(sn):.1e}, lapse_ratio - 1 {a['lapse_ratio'] - 1:.1e}")

    P("\n# atol: scipy's star, rtol 1e-13, atol 1e-14 (P0) against atol 1e-17")
    for g, v in [(1000.0, "base"), (1000.0, "phi"), (3000.0, "base")]:
        a, b = sget(sc, g, v, 1e-13), sget(sc, g, v, 1e-13, atol=1e-17)
        if a and b:
            P(f"  gamma {g:g} {v}: gap {gap(a['D'], ref(v, g)):.1e} -> {gap(b['D'], ref(v, g)):.1e}; D moves {abs(b['D'] / a['D'] - 1):.1e}")

    import glob as _glob
    prs = []
    PLANES_FILES[:] = sorted(_glob.glob(os.path.join(args.dir, "planes*.jsonl")))
    for fn in PLANES_FILES:
        prs += [json.loads(x) for x in open(fn) if x.startswith("{")]
    if prs:
        P("\n# the planes: sine of the largest principal angle between scipy's plane and the crate's double-double plane at r_m (interior Y, exterior Z, and Z's metric direction and scalar plane)")
        P(f"{'gamma':>6s} {'variant':>7s} {'plant':>14s} {'rtol':>6s} {'extra':>5s} {'eps_Y':>9s} {'eps_Z':>9s} {'eps_met':>9s} {'eps_sca':>9s} {'sin min':>9s} {'signed gap':>11s}")
        for o in prs:
            rfv = ref(o["variant"], o["gamma"]) if ("dd", o["variant"], o["gamma"]) in crate else None
            P(f"{o['gamma']:6g} {o['variant']:>7s} {o['plant']:>14s} {o['rtol']:6.0e} {o['extra_chunks']:5d} {o['eps_y']:9.1e} {o['eps_z']:9.1e} {o['eps_metric']:9.1e} {o['eps_scalar']:9.1e} {min(o['sines']):9.1e} "
              + (f"{float(mp.mpf(o['D']) / rfv - 1):+11.2e}" if rfv is not None else f"{'-':>11s}"))
        P("# medians of eps_Y per (gamma, variant, plant, rtol)")
        grp = {}
        for o in prs:
            grp.setdefault((o["gamma"], o["variant"], o["plant"], o["rtol"]), []).append(o["eps_y"])
        for k, v in grp.items():
            P(f"  {k}: median eps_Y {float(np.median(v)):.2e} over {len(v)} realizations")
    pl = os.path.join(args.dir, "planted.jsonl")
    plc = {}
    if os.path.exists(pl):
        P("\n# scipy's integrator on the planted system (taylor_polar P1's A, L = 24.1), error against the invariant subspaces at 90 digits")
        for line in open(pl):
            o = json.loads(line)
            nz = o.get("noise_ulp", 0.0)
            P(f"  gamma {o['gamma']:g} rtol {o['rtol']:g} atol {o['atol']:g}" + (f" planted noise {nz:g} ulp seed {o['noise_seed']}" if nz else "")
              + f": D {o['D']:+.10e} exact {o['D_exact'][:16]}  error {o['err']:.1e}"
              + (f"  sines {o['sines'][0]:.2e} {o['sines'][1]:.2e} {o['sines'][2]:.2e}, error * min sine {o['err'] * min(o['sines']):.1e}" if "sines" in o else "")
              + f"  ({o['wall']:.0f} s)")
            if not nz:
                plc[(o["gamma"], o["rtol"])] = o["err"]

    P("\n# CJE: the exact-J system (crate's star, M_ext, J tables), 8 realizations (extra chunks 0 .. 7) per cell and setting: signed D/D_crate - 1")
    P(f"{'cell':>14s} {'rtol':>7s} {'atol':>7s} {'mean':>10s} {'std':>10s} {'t=mean/(std/sqrt n)':>20s}   realizations")
    for (g, v) in CJ_CELLS:
        for (r, a) in CJE_SETTINGS:
            try:
                xs, m, sd = cje_stats(sc, crate, g, v, r, a)
                P(f"{'%g %s' % (g, v):>14s} {r:7.0e} {a:7.0e} {m:+10.2e} {sd:10.2e} {m / (sd / np.sqrt(len(xs))):+20.1f}   " + " ".join(f"{x:+.1e}" for x in xs))
            except TypeError:
                P(f"{'%g %s' % (g, v):>14s} {r:7.0e} {a:7.0e}   not run")
    P("\n# CJV: as CJE with the exterior scalar plane in vacuum (crate_j_vac), rtol 1e-13, atol 1e-14")
    import mpmath as mp4
    for (g, v) in CJ_CELLS:
        xs = []
        for e in range(CJE_N):
            rec = sget(sc, g, v, 1e-13, "crate_j_vac", extra=e)
            if rec:
                xs.append(float(mp4.mpf(rec["D"]) / ref(v, g) - 1))
        if len(xs) >= 2:
            m, sd = float(np.mean(xs)), float(np.std(xs, ddof=1))
            P(f"  {'%g %s' % (g, v):>14s}: n {len(xs)} mean {m:+.2e} std {sd:.2e} t {m / (sd / np.sqrt(len(xs))):+.1f}   " + " ".join(f"{x:+.1e}" for x in xs))
    P("\n# CJR: as CJE with the Riccati of the exterior q at rtol 1e-13, atol 1e-16 (crate_j_r13), interior rtol 1e-13, atol 1e-14")
    import mpmath as mp6
    for (g, v) in CJ_CELLS:
        xs = []
        for e in range(CJE_N):
            rec = sget(sc, g, v, 1e-13, "crate_j_r13", extra=e)
            if rec:
                xs.append(float(mp6.mpf(rec["D"]) / ref(v, g) - 1))
        if len(xs) >= 2:
            m, sd = float(np.mean(xs)), float(np.std(xs, ddof=1))
            P(f"  {'%g %s' % (g, v):>14s}: n {len(xs)} mean {m:+.2e} std {sd:.2e} t {m / (sd / np.sqrt(len(xs))):+.1f}   " + " ".join(f"{x:+.1e}" for x in xs))
    P("\n# CJQ: D's sensitivity to a planted error of q (the same chunk grid, so no realization noise): (D(q (1 + eps)) / D(q) - 1) / eps, rtol 1e-10")
    for g in (1000.0, 3000.0):
        d0 = sget(sc, g, "base", 1e-10, "crate_j")
        for pl, eps in (("crate_j_q13", 1e-13), ("crate_j_q12", 1e-12)):
            d1 = sget(sc, g, "base", 1e-10, pl)
            if d0 and d1:
                P(f"  gamma {g:g} eps {eps:g}: D moves {abs(d1['D'] / d0['D'] - 1):.2e}, sensitivity {(d1['D'] / d0['D'] - 1) / eps:+.3g}")
    P("\n# CJN: planted noise on the exact-J system (1000 base, rtol 1e-10, atol 1e-14): std of the signed gap over seeds 1 .. 4 (noise 0: the 8 realizations of CJE)")
    try:
        xs0, m0, s0 = cje_stats(sc, crate, 1000.0, "base", 1e-10, 1e-14)
        P(f"  amplitude 0: std {s0:.2e} (mean {m0:+.2e})")
    except TypeError:
        pass
    import mpmath as mp2
    for kind, u in CJN_NOISE:
        xs = []
        for sd in (1, 2, 3, 4):
            rec = sget(sc, 1000.0, "base", 1e-10, "crate_j", noise=(u if kind == "noise" else 0.0), seed=sd, row=(u if kind == "row" else 0.0))
            if rec:
                xs.append(float(mp2.mpf(rec["D"]) / ref("base", 1000.0) - 1))
        if len(xs) >= 2:
            P(f"  {'per-entry' if kind == 'noise' else 'row-max'} {u:g} ulp: std {float(np.std(xs, ddof=1)):.2e}, mean {float(np.mean(xs)):+.2e}, n {len(xs)}: " + " ".join(f"{x:+.1e}" for x in xs))
    try:
        st, lines = verdicts(sc, crate, plc)
        P("\n# verdicts by the declared rules")
        for x in lines:
            P("  " + x)
        P(f"  cause: {name_cause(st)}   ({st})")
    except (TypeError, KeyError) as e:
        P(f"\n# verdicts: not computed, a case is missing ({type(e).__name__}: {e})")


# ---------------------------------------------------------------------------------------------------------------
# the rounding error of the generated J in f64: J_and_h at the inputs the integrator gives it, against the same source at 60 digits
# ---------------------------------------------------------------------------------------------------------------
def run_jround(args):
    """for each gamma: over `n` radii in [0.5, r_m] and the 36 entries of J (l = 5, scipy's star): |J_f64 - J_60digits| / |J_60digits| in units of 2^-53,
    rms and max over the entries that are not exactly zero, the worst entry, and the same error relative to the largest entry of its row"""
    import mpmath as mp
    import taylor_polar as TP
    from polar_matrix_generated import J_and_h
    mp.mp.dps = 60
    f = TP.mp_J()
    st = scipy_star()
    l = 5
    LL = l * (l + 1)
    import background as B
    for g in args.gammas:
        ulps, rowrel, worst = [], [], (0.0, None)
        for r in np.linspace(0.5, RM["phi"], args.n):
            al, a, ph, u = (float(x[0]) for x in st.bg(np.array([r])))
            Jf, _ = J_and_h(r, al, a, ph, u, OMEGA, EPS, g, LL, B.LAM, B.MU)
            Jm, _ = f(mp.mpf(r), mp.mpf(al), mp.mpf(a), mp.mpf(ph), mp.mpf(u), mp.mpf(OMEGA), mp.mpf(EPS), mp.mpf(g), mp.mpf(LL), mp.mpf(B.LAM), mp.mpf(B.MU))
            for i in range(6):
                rmax = max(abs(Jm[i, j]) for j in range(6))
                for j in range(6):
                    if Jm[i, j] == 0:
                        continue
                    e = float(abs(mp.mpf(float(Jf[i, j])) - Jm[i, j]) / abs(Jm[i, j])) / 2.0 ** -53
                    ulps.append(e)
                    rowrel.append(float(abs(mp.mpf(float(Jf[i, j])) - Jm[i, j]) / rmax) / 2.0 ** -53)
                    if e > worst[0]:
                        worst = (e, (round(float(r), 3), i, j))
        ulps, rowrel = np.array(ulps), np.array(rowrel)
        rec = dict(cmd="jround", gamma=g, l=l, radii=args.n, entries=int(ulps.size), rms_ulp=float(np.sqrt(np.mean(ulps ** 2))), max_ulp=float(ulps.max()),
                   rms_ulp_of_row_max=float(np.sqrt(np.mean(rowrel ** 2))), max_ulp_of_row_max=float(rowrel.max()), worst_entry=worst[1])
        print(json.dumps(rec), flush=True)


# ---------------------------------------------------------------------------------------------------------------
# the planes themselves: scipy's two orthonormal 3-planes at r_m against the crate's (double-double, `rows planes=1`)
# ---------------------------------------------------------------------------------------------------------------
def crate_planes(l, gamma, variant):
    """the crate's orthonormal interior and exterior 3-planes at r_m as float64 (6 x 3 each), from the rows command with planes=1"""
    import subprocess
    import taylor_polar as TP
    args = [TP.BIN, "rows", f"l={l}", "lane=dd", f"variants={variant}", f"gammas={gamma!r}", "planes=1"] + [f"{k}={v}" for k, v in TP.star_kw().items()]
    p = subprocess.run(["nice", "-n", "10"] + args, capture_output=True, text=True, env=dict(os.environ, RAYON_NUM_THREADS="2"), cwd=HERE)
    if p.returncode:
        raise RuntimeError(p.stderr[-400:])
    for line in p.stdout.splitlines():
        if line.startswith("{"):
            o = json.loads(line)
            if o["cmd"] == "planes":
                f = lambda m: np.array([[float(x) for x in c] for c in m]).T
                return f(o["y"]), f(o["z"])
    raise RuntimeError("no planes line")


def sin_max(A, B):
    """sine of the largest principal angle between the column spaces of A and B (equal dimensions, A orthonormal columns)"""
    Qa, _ = np.linalg.qr(A)
    Qb, _ = np.linalg.qr(B)
    R = Qb - Qa @ (Qa.T @ Qb)
    return float(np.linalg.svd(R, compute_uv=False)[0])


def split_exterior(Z):
    """the exterior plane = span{(m, 0), (0, S)}: the metric direction (2-vector, from the rank-1 block Z[:2]) and the scalar plane (4 x 2, from Z[2:])"""
    u, sv, _ = np.linalg.svd(Z[:2, :], full_matrices=False)
    us, ss, _ = np.linalg.svd(Z[2:, :], full_matrices=False)
    return u[:, 0], us[:, :2], sv, ss


def run_planes(args):
    DP_plant = args.plant
    for g in args.gammas:
        yc, zc = crate_planes(args.l, g, args.variant)
        mc, Sc, _, _ = split_exterior(zc)
        for e in args.extra_chunks:
            t = time.time()
            pp = build(args.l, args.variant, DP_plant, args.rtol, args.steps, args.atol, e, None, 0.0, 0, 0.0)
            D = pp.D(g)
            ys, zs = pp.planes
            ms, Ss, _, _ = split_exterior(zs)
            ang_m = float(abs(np.sqrt(max(0.0, 1.0 - (ms @ mc) ** 2))))
            rec = dict(gamma=g, variant=args.variant, plant=DP_plant, rtol=args.rtol, atol=args.atol, extra_chunks=e, D=float(D),
                       eps_y=sin_max(ys, yc), eps_z=sin_max(zs, zc), eps_metric=ang_m, eps_scalar=sin_max(Ss, Sc),
                       sines=pp.sines(), wall=time.time() - t, host=os.uname().nodename)
            print(json.dumps(rec), flush=True)
            with open(args.out, "a") as fp:
                fp.write(json.dumps(rec) + "\n")


def case_key(o):
    return (o["l"], o["gamma"], o["variant"], o["rtol"], o["plant"], o.get("atol", 1e-14), o.get("extra_chunks", 0), o.get("rm_from"),
            o.get("noise_ulp", 0.0), o.get("noise_seed", 0), o.get("noise_row", 0.0))


def floats(s):
    return [float(x) for x in s.split(",")]


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("scipy")
    s.add_argument("--out", required=True)
    s.add_argument("--steps", required=True, help="crate star steps (star-steps)")
    s.add_argument("--l", type=int, default=5)
    s.add_argument("--gammas", type=floats, default=[1000.0])
    s.add_argument("--variants", type=lambda x: x.split(","), default=["base"])
    s.add_argument("--rtols", type=floats, default=[1e-10])
    s.add_argument("--plants", type=lambda x: x.split(","), default=["sci"])
    s.add_argument("--atol", type=float, default=1e-14)
    s.add_argument("--extra-chunks", dest="extra_chunks", type=lambda x: [int(y) for y in x.split(",")], default=[0], help="a list: one case per value")
    s.add_argument("--rm-from", dest="rm_from", default=None)
    s.add_argument("--jtable", default=None, help="crate J tables (j-table), for the plant crate_j")
    s.add_argument("--noise-ulp", dest="noise_ulp", type=float, default=0.0)
    s.add_argument("--noise-row", dest="noise_row", type=float, default=0.0)
    s.add_argument("--noise-seeds", dest="noise_seeds", type=lambda x: [int(y) for y in x.split(",")], default=[0])
    s.add_argument("--max-seconds", dest="max_seconds", type=float, default=7200)
    c = sub.add_parser("crate")
    c.add_argument("--out-dir", dest="out_dir", required=True)
    c.add_argument("--l", type=int, default=5)
    c.add_argument("--gammas", type=floats, default=[1000.0, 2000.0, 3000.0])
    c.add_argument("--qd", action="store_true", help="also the quad-double lane, at --qd-gammas only")
    c.add_argument("--qd-gammas", dest="qd_gammas", type=floats, default=[1000.0])
    q = sub.add_parser("planted")
    q.add_argument("--out", required=True)
    q.add_argument("--gammas", type=floats, default=[1000.0, 2000.0, 3000.0])
    q.add_argument("--rtols", type=floats, default=[1e-10, 1e-13])
    q.add_argument("--atol", type=float, default=1e-14)
    q.add_argument("--L", type=float, default=24.1)
    q.add_argument("--noise-ulp", dest="noise_ulp", type=float, default=0.0)
    q.add_argument("--noise-seeds", dest="noise_seeds", type=lambda x: [int(y) for y in x.split(",")], default=[0])
    k = sub.add_parser("check")
    k.add_argument("--dir", required=True, help="folder of the scipy*.jsonl files")
    k.add_argument("--crate-dir", dest="crate_dir", default=None, help="folder of the crate's rows_l5_*.jsonl (default --dir)")
    k.add_argument("--steps", required=True)
    k.add_argument("--raw", default=os.path.normpath(os.path.join(HERE, "..", "..", "..", "negative-mass-raw")))
    k.add_argument("--only", default=None)
    k.add_argument("--mutate", default=None)
    k.add_argument("--synthetic", action="store_true")
    k.add_argument("--jtable", default=None)
    pn = sub.add_parser("planes")
    pn.add_argument("--out", required=True)
    pn.add_argument("--steps", required=True)
    pn.add_argument("--jtable", default=None)
    pn.add_argument("--l", type=int, default=5)
    pn.add_argument("--gammas", type=floats, default=[3000.0])
    pn.add_argument("--variant", default="base")
    pn.add_argument("--plant", default="crate_j")
    pn.add_argument("--rtol", type=float, default=1e-13)
    pn.add_argument("--atol", type=float, default=1e-14)
    pn.add_argument("--extra-chunks", dest="extra_chunks", type=lambda x: [int(y) for y in x.split(",")], default=[0])
    jr = sub.add_parser("jround")
    jr.add_argument("--gammas", type=floats, default=[1000.0, 3000.0])
    jr.add_argument("--n", type=int, default=40)
    rp = sub.add_parser("report")
    rp.add_argument("--dir", required=True)
    rp.add_argument("--crate-dir", dest="crate_dir", default=None)
    t = sub.add_parser("star-steps")
    t.add_argument("out")
    jt = sub.add_parser("j-table")
    jt.add_argument("out")
    jt.add_argument("--l", type=int, default=5)
    a = ap.parse_args()
    if a.cmd == "scipy":
        JTABLE_PATH[0] = a.jtable
        run_cases(a)
    elif a.cmd == "planes":
        JTABLE_PATH[0] = a.jtable
        run_planes(a)
    elif a.cmd == "jround":
        run_jround(a)
    elif a.cmd == "report":
        run_report(a)
    elif a.cmd == "check":
        run_checks(a)
    elif a.cmd == "crate":
        run_crate(a)
    elif a.cmd == "planted":
        run_planted(a)
    elif a.cmd == "j-table":
        import taylor_polar as TP
        res = TP.rust("jtable", l=a.l, rmin=5e-4, rmax=34.7, **TP.star_kw())
        with open(a.out, "w") as fp:
            for o in res:
                fp.write(json.dumps(o) + "\n")
        print(f"wrote {a.out}: {res[0]['steps']} table steps, {len(res[0]['pattern'])} entries, linear defect {res[0]['linear_defect']:.1e}")
    elif a.cmd == "star-steps":
        import taylor_polar as TP
        kw = TP.star_kw()
        res = TP.rust("starsteps", rend=62.0, **kw)
        with open(a.out, "w") as fp:
            for o in res:
                fp.write(json.dumps(o) + "\n")
        print(f"wrote {a.out}: {len(res) - 1} steps to r = {res[0]['r_end']}")


if __name__ == "__main__":
    main()

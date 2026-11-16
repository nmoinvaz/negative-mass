#!/usr/bin/env python3
"""Follow the unfolded phantom family (eps = -1) in the central lapse alpha_c beyond 1.5, with ω free.

    python phantom_follow.py --out follow.jsonl [--alpha-stop 20] [--hmax 0.01] [--resume ckpt.npz]
    python phantom_follow.py --selftest
    python phantom_follow.py --check-R follow.jsonl [--mutate=matching]

The background equations and boundary conditions are background.py's (odes, bcs_lapse); only the outer radius R is
moved as the star grows. The BVP of background.solve_star_lapse keeps R of its guess, which is a fixed 60 in
strongfield.py, and the star's wall reaches r = 46 at alpha_c = 1.5. Here R = rs + MARGIN/k at every star (rs is where
phi = 1e-3 phi_c, radial_modes' matching radius), and the mesh and the tail are extended when R grows.

Each accepted star is one JSON line (flushed), and a checkpoint (the BVP mesh and solution) is rewritten every
CKPT stars so a restart continues from it. When the step falls below HMIN the run stops, writes the failed solve's
message and mesh size, and tries the diagnostics in `diagnose`.
"""
import argparse
import json
import math
import os
import resource
import sys
import time
import types
import warnings

import numpy as np
from scipy.integrate import solve_bvp
from scipy.interpolate import CubicHermiteSpline

warnings.filterwarnings("ignore")
import background as B  # noqa: E402

EPS = -1
MARGIN = 30.0            # R - rs >= MARGIN / k, with k = sqrt(1 - w^2) (rs: phi = 1e-3 phi_c)
HMIN = 1e-6


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def matching_radius(st):
    """radial_modes.RadialProblem's rm: first r where phi0 < 1e-3 phi_c on a 4000-point grid."""
    rr = np.linspace(st.r0, st.R, 4000)
    phi = st.bg(rr)[2]
    return float(rr[np.argmax(phi < 1e-3 * phi[0])])


def half_radius(st):
    rr = np.linspace(st.r0, st.R, 8000)
    phi = st.bg(rr)[2]
    i = int(np.argmax(phi < 0.5 * phi[0]))
    return float(rr[i])


def required_R(st):
    k = math.sqrt(1.0 - st.w ** 2)
    return matching_radius(st) + MARGIN / k


class Guess:
    """A converged BVP solution restored from a checkpoint, with what solve_lapse needs of a Star."""
    def __init__(self, x, y, yp, p, R):
        self.sol = types.SimpleNamespace(x=x, sol=CubicHermiteSpline(x, y, yp, axis=1))
        self.w, self.R = float(p[0]), float(R)


def save_checkpoint(path, st, alpha_c):
    s = st.sol
    np.savez(path, x=s.x, y=s.y, yp=s.yp, p=s.p, R=st.R, alpha_c=alpha_c)


def load_checkpoint(path):
    z = np.load(path)
    return Guess(z["x"], z["y"], z["yp"], z["p"], float(z["R"])), float(z["alpha_c"])


def bcs_lapse(ya, yb, p, r0, R, eps, nuc, mutate=None, tail="yukawa"):
    res = B.bcs_lapse(ya, yb, p, r0, R, eps, nuc)
    if mutate == "matching":                       # drop the Schwarzschild matching of nu(R)
        res = res.copy()
        res[3] = yb[3]
    if tail == "wkb":
        # diagnostic only (phantom_stop.py): u(R) = -kappa - 1/R with the local kappa = a sqrt(1 - w^2/alpha^2), which
        # has no k = sqrt(1 - w^2) in it and so no branch point at w = 1 (alpha^2 = 1 - 2m/R > 1 for m < 0)
        res = res.copy()
        a2 = 1.0 / (1.0 - 2.0 * yb[2] / R)
        kap = np.sqrt(a2 * (1.0 - p[0] ** 2 / np.exp(yb[3])))
        res[2] = yb[1] + kap + 1.0 / R
    return res


def solve_lapse(alpha_c, eps, guess, R, r0=1e-3, tol=1e-8, max_nodes=200000, mutate=None, tail="yukawa"):
    """Like background.solve_star_lapse, with the outer radius R free. A larger R than the guess's extends the
    mesh and the Yukawa tail (psi falling as e^{-k r}, u = -k - 1/r, m constant, nu = ln(1 - 2m/r))."""
    xo = guess.sol.x
    idx = np.unique(np.linspace(0, xo.size - 1, min(xo.size, 1500)).astype(int))
    x = xo[idx]
    x = x[x <= R]
    gR = min(guess.R, R)
    if x[-1] < gR:
        x = np.append(x, gR)
    if R > guess.R:
        x = np.concatenate([x, np.linspace(guess.R, R, max(3, int((R - guess.R) / 2)) + 1)[1:]])
    elif x[-1] < R:
        x = np.append(x, R)
    y0 = guess.sol.sol(np.minimum(x, guess.R))[:4]
    k = math.sqrt(1.0 - guess.w ** 2) if guess.w < 1 else 1e-3
    far = x > guess.R
    if far.any():
        y0[1, far] = -k - 1.0 / x[far]
        y0[0, far] = y0[0, ~far][-1] - k * (x[far] - guess.R)
        y0[2, far] = y0[2, ~far][-1]
        y0[3, far] = np.log(1 - 2 * y0[2, far] / x[far])
    sol = solve_bvp(lambda r, y, p: B.odes(r, y, p[0], eps),
                    lambda ya, yb, p: bcs_lapse(ya, yb, p, r0, R, eps, 2 * np.log(alpha_c), mutate, tail),
                    x, y0, p=[guess.w], tol=tol, max_nodes=max_nodes, bc_tol=1e-10)
    if sol.status != 0:
        err = RuntimeError(f"bvp(alpha_c) failed alpha_c={alpha_c} R={R:.1f}: {sol.message}")
        err.status, err.nmesh = sol.status, sol.x.size
        raise err
    return B.StarP(sol, eps, r0, R)


def row_of(st, alpha_c, h, dt):
    r = np.linspace(st.r0, st.R, 6000)
    al, a, phi, u = st.bg(r)
    m = st.sol.sol(r)[2]
    return dict(alpha_c=alpha_c, w=st.w, k=math.sqrt(1.0 - st.w ** 2), N=st.N, M=st.M, phic=st.phic, rs=matching_radius(st),
                r_half=half_radius(st), max2mr=float(np.max(np.abs(2 * m / r))), alpha_wall=float(al[np.argmax(phi < 0.5 * phi[0])]),
                R=st.R, nmesh=st.nmesh, max_rms=st.max_res, h=h, t_solve=dt, rss_mb=rss_mb())


def start_star(alpha_c_start, eps=EPS):
    """The same path as strongfield.py: w-continuation from 0.66 to 0.45, then alpha_c-continuation."""
    st = B.solve_star(0.66, eps, tol=1e-8)
    w, dw = 0.66, 0.01
    while w > 0.45 + 1e-9:
        w2 = round(max(w - dw, 0.45), 6)
        try:
            st = B.solve_star(w2, eps, guess=st, tol=1e-8); w = w2; dw = min(dw * 1.5, 0.01)
        except RuntimeError:
            dw /= 2
    ac, d = st.alpha_c, 0.01
    while abs(alpha_c_start - ac) > 1e-12:
        step = min(d, alpha_c_start - ac)
        try:
            st = B.solve_star_lapse(ac + step, eps, st); ac += step; d = min(d * 1.3, 0.01)
        except RuntimeError:
            d /= 2
    return st, ac


def follow(args):
    mode = "a" if args.resume else "w"
    targets = sorted(float(x) for x in args.targets.split(",")) if args.targets else []
    if args.resume:
        guess, ac = load_checkpoint(args.resume)
        R = guess.R
        st = solve_lapse(ac, EPS, guess, R)
    else:
        st, ac = start_star(args.alpha_start)
        R = st.R
    out = open(args.out, mode)
    if args.snap_dir:
        os.makedirs(args.snap_dir, exist_ok=True)
    if not args.resume:
        row = row_of(st, ac, 0.0, 0.0)
        row["target"] = any(abs(ac - t) < 1e-12 for t in targets)
        out.write(json.dumps(row) + "\n"); out.flush()
        print(f"start alpha_c={ac:.6f} w={st.w:.6f} N={st.N:.5f} M={st.M:+.6f} R={st.R:.0f}", flush=True)
        if row["target"] and args.snap_dir:
            save_checkpoint(os.path.join(args.snap_dir, f"star_a{ac:.4f}.npz"), st, ac)
    h, nstars = args.h0, 0
    while ac < args.alpha_stop - 1e-12 and nstars < args.max_stars:
        goal = min([t for t in targets if t > ac + 1e-12] + [args.alpha_stop])
        step = min(h, goal - ac)
        t0 = time.time()
        try:
            R_new = max(R, required_R(st))
            st2 = solve_lapse(ac + step, EPS, st, R_new)
        except RuntimeError as e:
            h = step / 2
            if h < HMIN:
                rec = dict(stop=True, alpha_c=ac, step_tried=step, message=str(e), status=getattr(e, "status", None),
                           nmesh=getattr(e, "nmesh", None), R=R, k=math.sqrt(1.0 - st.w ** 2), rss_mb=rss_mb())
                out.write(json.dumps(rec) + "\n"); out.flush()
                print("STOP", rec, flush=True)
                save_checkpoint(args.ckpt_path, st, ac)
                return st, ac, rec
            continue
        landed = abs(ac + step - goal) < 1e-12
        ac = goal if landed else ac + step
        st, R = st2, st2.R
        h = min(h * 1.3, args.hmax)
        nstars += 1
        row = row_of(st, ac, step, time.time() - t0)
        row["target"] = bool(landed and goal in targets)
        out.write(json.dumps(row) + "\n"); out.flush()
        print(f"alpha_c={ac:.4f} w={st.w:.6f} N={st.N:.5f} M={st.M:+.6f} phic={st.phic:.6f} rs={row['rs']:.1f} "
              f"r_half={row['r_half']:.1f} max|2m/r|={row['max2mr']:.3f} R={R:.0f} mesh={st.nmesh} h={step:.4f} "
              f"t={row['t_solve']:.1f}s rss={rss_mb():.0f}MB", flush=True)
        if row["target"] and args.snap_dir:
            save_checkpoint(os.path.join(args.snap_dir, f"star_a{ac:.4f}.npz"), st, ac)
        if nstars % args.ckpt == 0:
            save_checkpoint(args.ckpt_path, st, ac)
    save_checkpoint(args.ckpt_path, st, ac)
    return st, ac, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="follow.jsonl")
    ap.add_argument("--alpha-start", type=float, default=1.1)
    ap.add_argument("--alpha-stop", type=float, default=20.0)
    ap.add_argument("--hmax", type=float, default=0.01)
    ap.add_argument("--h0", type=float, default=0.01)
    ap.add_argument("--max-stars", type=int, default=10 ** 9)
    ap.add_argument("--resume", default=None)
    ap.add_argument("--ckpt", type=int, default=10)
    ap.add_argument("--ckpt-path", default="follow_ckpt.npz")
    ap.add_argument("--targets", default="", help="alpha_c values the steps land on; a star snapshot is saved at each")
    ap.add_argument("--snap-dir", default=None)
    args = ap.parse_args()
    follow(args)


if __name__ == "__main__":
    main()

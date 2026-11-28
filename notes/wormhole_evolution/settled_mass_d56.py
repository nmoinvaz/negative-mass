#!/usr/bin/env python3
"""Settled black-hole mass at D = 5 and 6 (#246): #41's apparent-horizon chain (settled_mass.py) for families S and G at n = 3, 4, with
M = r^(n-1)/2 on the apparent horizon, and a planted exact Tangherlini through the same chain at each D.

A thin driver. It imports settled_mass.py (#41), settled_mass_g.py (#253), threshold.py (#40) and dn.py (#39) unchanged and adds only:
  - #40's p*_h of S and G at D = 5, 6 read from its raw (the h of each job.log, as the log records it);
  - the general-n run (`run_real_d`, #253's `run_real_g` with n and the family as arguments) and the general-n first law (r^n/n);
  - the planted exact Tangherlini, D = 5, 6, in a gauge in which e^{-f} on the horizon grows as e^{kp (x+ - xs)} (the stretching of the real run).

    python3 settled_mass_d56.py pstar   --n=N --fam=F --k=K --out=DIR       read p*_h (K = 0, 1, 2: coarse to fine) from #40's raw
    python3 settled_mass_d56.py job     --n=N --fam=F --k=K --delta=D --out=DIR --pstar=DIR [--noise=U --seed=S]
    python3 settled_mass_d56.py planted [--n=N] [--mutate=NAME]

Units: a = 1, 8 pi G = 1, D = n + 2. M = r^(n-1)/2 (1 + 2 e^f r_+ r_-); on r_+ = 0 it is r^(n-1)/2.
Status words as in AGENTS.md: every mass is numerical evidence.
"""
import hashlib
import json
import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dn                      # noqa: E402
import settled_mass as sm      # noqa: E402
import settled_mass_g as smg   # noqa: E402  (its import sets sm.track = smg.track_g: the normal field's phi_+ on the horizon)
import threshold as th         # noqa: E402

SQ2 = math.sqrt(2.0)
MUT = None
I40 = os.path.expanduser("~/Source/negative-mass-raw/I40-dh610/runs")
FAMS = ("S", "G")
NS = (3, 4)                                   # D = 5, 6


def mut(name):
    return MUT == name


# --------------------------------------------------------------------------------------------------------------------
# #40's p*_h and its jobs
# --------------------------------------------------------------------------------------------------------------------
def i40_dirs(n, fam):
    """The directories of #40's jobs for (n, fam), coarse to fine (h falls by 2): D<D>_<fam>_h<h:g>."""
    D = n + 2
    ds = [d for d in os.listdir(I40) if re.fullmatch(rf"D{D}_{fam}_h[0-9.]+", d)]
    return sorted(ds, key=lambda d: -float(d.split("_h")[1]))


def i40_job(n, fam, k):
    """#40's bracket, p*_h, formation values and the h of its job.log (as a string) of job k (0 = coarsest) at (n, fam)."""
    p = f"{I40}/{i40_dirs(n, fam)[k]}"
    s = json.load(open(p + "/summary.json"))
    st = json.load(open(p + "/stall.json"))["-1"]
    hstr = re.search(r"'h': '([^']+)'", open(p + "/job.log").readline()).group(1)
    ps = 0.5 * (s["lo"] + s["hi"]) + (3e-10 if mut("pstar_shift") else 0.0)
    return dict(ps=ps, lo=s["lo"], hi=s["hi"], stopped=s["stopped"], h=float(hstr), h_log=hstr, h_summary=s["h"], X=s["X"], event=st["event"], T=st["T"], M_hor=st["M_hor"],
                M_cell=st["M_cell"], r_cell=st["r_cell"], A_stall=st["A"], path=p, job_secs=s.get("secs"), bisect_secs=s.get("bisect_secs"))


def n_h(n, fam):
    return len(i40_dirs(n, fam))


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def pstar_cmd(n, fam, k, out):
    """p*_h = (lo + hi)/2 of #40's job, read from its raw folder (not bisected again), with the hashes of what was read."""
    os.makedirs(out, exist_ok=True)
    r = i40_job(n, fam, k)
    res = dict(n=n, D=n + 2, fam=fam, k=k, h=r["h"], h_job_log=r["h_log"], lo=r["lo"], hi=r["hi"], ps=0.5 * (r["lo"] + r["hi"]), stopped=r["stopped"], X=r["X"],
               event=r["event"], T=r["T"], M_hor=r["M_hor"], source=r["path"].replace(os.path.expanduser("~"), "~"),
               sha256={f: sha(os.path.join(r["path"], f)) for f in ("summary.json", "stall.json", "job.log", "points.jsonl")})
    json.dump(res, open(os.path.join(out, f"pstar_D{n + 2}_{fam}_k{k}.json"), "w"), indent=1)
    print(json.dumps(res), flush=True)


# --------------------------------------------------------------------------------------------------------------------
# the chain's constants: #41's, with D-dependent entries declared in the raw README
# --------------------------------------------------------------------------------------------------------------------
CRITD = dict(F1=1.5, K=20, W1=2.4, rcut=0.05, tol=1e-2, dxp=12.0, dm=3.0)       # #41's declared chain; development runs may change entries per D


def nc_of(h, W1=None):
    """Cells per box: level 1 is W1 wide in x- (the horizon moves 0.1 W1 inside it) at width h/2, so Nc = 2 W1 / h rounded to a multiple of 4."""
    W1 = CRITD["W1"] if W1 is None else W1
    return int(round(2.0 * W1 / h / 4.0)) * 4


def run_real_d(n, fam, h, delta, ps, K=None, Nc=None, dxp=None, dm=None, log=lambda *a: None, need_up=True, noise=0.0, seed=0, forced=None):
    """Family `fam` at D = n + 2, p*_h - delta (collapse side): #40's point() for the fate and the first trapped cell, then #41's chain from there."""
    F1 = CRITD["F1"]
    K = CRITD["K"] if K is None else K
    Nc = nc_of(h) if Nc is None else Nc
    dxp = CRITD["dxp"] if dxp is None else dxp
    dm = CRITD["dm"] if dm is None else dm
    A = ps - delta
    pulses, sym, eps = th.family(fam, A)
    X = th.CRIT["X"][fam][n]
    pt = th.point(n, fam, h, X, A)
    up = th.point(n, fam, h, X, ps + delta) if need_up else dict(fate=1)
    if pt["fate"] != -1 or up["fate"] != 1:
        raise RuntimeError(f"fates at p*_h -/+ delta: {pt['fate']}, {up['fate']} (wanted -1, +1)")
    j_ev, i_ev = pt["j"], pt["i"]
    Np = i_ev + int(round(dxp / h))
    nm0 = j_ev + int(round((dm + Nc * h / 4) / h))
    d = dn.ellis_rays(n, max(Np, nm0), h, pulses=pulses, symmetric=sym, eps=eps)
    smg.add_noise(d, noise * smg.ULP, seed)
    out = sm.chain(n, eps, h, d.Sp[:, :Np + 1], d.Sm, 0.0, 0.0, nm0, None, None, F1, F1 * 2.0 ** K, Nc, forced=forced, tol=CRITD["tol"], rcut=CRITD["rcut"], log=log)
    out.update(point=pt, up=up, A=A, Np=Np, nm0=nm0, F1=F1, K=K, Nc=Nc, h=h, delta=delta, ps=ps, data=d, n=n, fam=fam, dxp=dxp, dm=dm)
    return out


def job_cmd(n, fam, k, delta, out, pstar_dir, noise=0.0, seed=0, **kw):
    import resource
    pj = json.load(open(os.path.join(pstar_dir, f"pstar_D{n + 2}_{fam}_k{k}.json")))
    h, ps = pj["h"], pj["ps"]
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    r = run_real_d(n, fam, h, delta, ps, log=lambda m: print(m, flush=True), noise=noise, seed=seed, **kw)
    s = r["ser"]
    np.savez(os.path.join(out, "series.npz"), **s)
    summ = dict(n=n, D=n + 2, fam=fam, k=k, h=h, delta=delta, ps=ps, A=r["A"], noise_ulp=noise, seed=seed, status=r["status"], event=r["event"], point=r["point"], up=r["up"],
                Np=r["Np"], nm0=r["nm0"], F1=r["F1"], K=r["K"], Nc=r["Nc"], dxp=r["dxp"], dm=r["dm"], levels=len(r["levels"]), cols=int(len(s["r"])),
                xp0=float(s["xp"][0]) if len(s["r"]) else None, xp1=float(s["xp"][-1]) if len(s["r"]) else None, F_end=float(s["F"][-1]) if len(s["r"]) else None,
                M_start=float(s["M"][0]) if len(s["r"]) else None, M_end=float(s["M"][-1]) if len(s["r"]) else None, min_edge=float(s["edge"].min()) if len(s["r"]) else None,
                secs=round(time.time() - t0, 1), maxrss_MB=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 1))
    json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print(json.dumps(summ), flush=True)


def load_run(raw, tag):
    p = os.path.join(raw, "runs", tag)
    return dict(np.load(os.path.join(p, "series.npz"))), json.load(open(os.path.join(p, "summary.json")))


def tag_of(n, fam, k, delta):
    return f"D{n + 2}_{fam}_k{k}_d{delta:g}"


# --------------------------------------------------------------------------------------------------------------------
# planted exact Tangherlini, D = 5, 6
# --------------------------------------------------------------------------------------------------------------------
def tang_Y(n, r):
    """Y(r) = -U V = (r - 1) e^{2 kappa R(r)} of the unit Tangherlini (r_h = 1, kappa = (n-1)/2, F = 1 - r^-(n-1)), in closed form for n = 3, 4:
    R(r) = int_1^r [1/F - 1/(2 kappa (s - 1))] ds, regular at r = 1. Y > 0 outside, < 0 inside, increasing on r > 0."""
    if n == 3:
        return 2.0 * (r - 1.0) * np.exp(2.0 * (r - 1.0)) / (r + 1.0)
    if n == 4:
        return (r - 1.0) * np.exp(3.0 * (r - 1.0)) * np.sqrt(3.0 / (r * r + r + 1.0)) * np.exp(-math.sqrt(3.0) * (np.arctan((2.0 * r + 1.0) / math.sqrt(3.0)) - math.pi / 3.0))
    raise ValueError(n)


def tang_exact_unit(n, Y):
    """r and e^{-f'} of the unit Tangherlini at Y = -UV (arrays): Y(r) = Y by bisection on (0.3, 80) (Y is increasing there), and
    e^{-f'} = F / (2 kappa^2 Y) (the metric -F du dv = -2 e^{-f'} dU dV with U = -e^{-kappa u}, V = e^{kappa v}); regular at r = 1, where it is 1/kappa."""
    kap = (n - 1.0) / 2.0
    Y = np.asarray(Y, float)
    lo = np.full_like(Y, 0.3)
    hi = np.full_like(Y, 80.0)
    for _ in range(110):
        mid = 0.5 * (lo + hi)
        up = tang_Y(n, mid) > Y
        hi = np.where(up, mid, hi)
        lo = np.where(up, lo, mid)
    r = 0.5 * (lo + hi)
    x = r - 1.0
    # F / Y near the horizon by the ratio of the regular parts: F = (1 - r^-(n-1)), Y = x e^{2 kappa R}; use the exact ratio away from x = 0
    F = 1.0 - r ** (-(n - 1.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        emf = np.where(np.abs(x) > 1e-7, F / (2.0 * kap ** 2 * tang_Y(n, r)), 1.0 / kap)
    return r, emf


def tang_exact(n, lam, kp, xh, xp, xm, xs=3.0):
    """Exact Tangherlini, horizon radius lam (M = lam^(n-1)/2), in the gauge x- = xh + lam U', V' = kappa e^{kp (x+ - xs)} / (lam kp): the horizon is
    x- = xh and e^{-f} on it is e^{kp (x+ - xs)}. Returns r, f at (xp, xm) arrays. (D = 4 is #41's `schw_exact`, with Lambert W.)"""
    kap = (n - 1.0) / 2.0
    xp, xm = np.broadcast_arrays(np.asarray(xp, float), np.asarray(xm, float))
    Up = (xm - xh) / lam
    lnVp = math.log(kap / (kp * lam)) + kp * (xp - xs)
    Y = -Up * np.exp(lnVp)
    r1, emf1 = tang_exact_unit(n, Y)
    lnemf = np.log(emf1) + math.log(kp * lam) + lnVp
    return lam * r1, -lnemf


def tang_rays(n, lam, kp, xh, hp, Np, nm0, xs=3.0):
    xs_ = np.arange(Np + 1) * hp
    ms = np.arange(nm0 + 1) * hp
    Sp = np.zeros((2, Np + 1))
    Sm = np.zeros((2, nm0 + 1))
    Sp[0], Sp[1] = tang_exact(n, lam, kp, xh, xs_, 0.0 * xs_, xs)
    Sm[0], Sm[1] = tang_exact(n, lam, kp, xh, 0.0 * ms, ms, xs)
    return Sp, Sm


def tang_planted(n, h, lam, kp, xh, Nc, F1=None, K=None, xs=3.0):
    """Exact Tangherlini in the stretching gauge through the same chain (n = 3, 4). Returns the series and the errors of M and of the horizon position."""
    F1 = CRITD["F1"] if F1 is None else F1
    K = CRITD["K"] if K is None else K
    Fmax = F1 * 2.0 ** K
    Xf = xs + math.log(Fmax) / kp + 1.0
    Np = int(round(Xf / h))
    nm0 = int(round(xh / h)) + Nc // 4 + 6
    Sp, Sm = tang_rays(n, lam, kp, xh, h, Np, nm0, xs)
    out = sm.chain(n, [], h, Sp, Sm, 0.0, 0.0, nm0, 0, xh / h, F1, Fmax, Nc, tol=0.0, log=lambda *x: None)
    s = out["ser"]
    Mex = lam ** (n - 1) / 2.0
    ok = len(s["M"]) > 0
    return out, dict(M_exact=Mex, e_M=np.abs(s["M"] - Mex) if ok else np.array([np.inf]), e_x=np.abs(s["xm"] - xh) if ok else np.array([np.inf]))


# --------------------------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------------------------
def main():
    global MUT
    argv = sys.argv[1:]
    opts = dict(a.lstrip("-").split("=", 1) if "=" in a else (a.lstrip("-"), True) for a in argv[1:])
    cmd = argv[0] if argv else "help"
    if "mutate" in opts:
        MUT = opts["mutate"]
        sm.MUT = MUT
        smg.MUT = MUT
        dn.MUT = MUT if MUT in ("theta_coef", "ghost_sign", "scalar_coef") else None
    import platform
    import numba as nb
    la = open("/proc/loadavg").read().split()[:3] if platform.system() == "Linux" else ["?"]
    print(f"settled_mass_d56.py {cmd} {opts}: numba {nb.__version__}, numpy {np.__version__}, python {platform.python_version()}, {platform.node()}, "
          f"load average at start {' '.join(la)}, mutation {MUT}, {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    if cmd == "pstar":
        pstar_cmd(int(opts["n"]), opts["fam"], int(opts["k"]), opts["out"])
    elif cmd == "job":
        kw = {k: float(opts[k]) for k in ("dxp", "dm") if k in opts}
        if "K" in opts:
            kw["K"] = int(opts["K"])
        job_cmd(int(opts["n"]), opts["fam"], int(opts["k"]), float(opts["delta"]), opts["out"], opts["pstar"], noise=float(opts.get("noise", 0.0)), seed=int(opts.get("seed", 0)), **kw)
    elif cmd == "tang_test":
        tang_selftest()
    else:
        raise SystemExit(__doc__)


def tang_selftest():
    """Development: the closed-form Tangherlini against dn.kruskal_exact (r < 6) and against the field equation by finite differences."""
    for n in NS:
        spl, g, Fx = dn.kruskal_tools(n)
        kap = (n - 1.0) / 2.0
        r = np.linspace(0.4, 5.9, 12)
        # dn: UV = G(r) = -(r-1) e^{2 kap R(r)}; here Y = -UV
        Yd = (r - 1.0) * np.exp(2 * kap * spl(r))
        print(f"n = {n}: max |Y closed form - Y of dn.kruskal_tools| over r in [0.4, 5.9]: {np.max(np.abs(tang_Y(n, r) - Yd) / np.maximum(1e-300, np.abs(Yd))):.2e}")
        r_, f_ = dn.kruskal_exact(n, np.array([-0.5, -0.5, -0.3]), np.array([0.2, 1.0, 0.05]))
        r2, emf2 = tang_exact_unit(n, -np.array([-0.5, -0.5, -0.3]) * np.array([0.2, 1.0, 0.05]))
        print(f"   dn.kruskal_exact r, f: {r_}, {f_}; closed form r, f: {r2}, {-np.log(emf2)}")


if __name__ == "__main__":
    main()

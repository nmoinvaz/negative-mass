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
CRITD = dict(F1=1.5, K=28, W1=2.4, rcut=0.05, tol=1e-2, dxp=14.0, dm=3.0)       # #41's chain with K = 28 (e^{-f} = 4.0e8) and dxp = 14 at D = 5, 6 (README: the x+ span after formation shrinks with D)
K41 = 20                                                                            # #41's limit: e^{-f} = F1 2^20 = 1 572 864, read from the same series (it is a bit-identical prefix)


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


def dstr(delta):
    """delta as the lane script writes it: 1e-6, 1e-8."""
    return f"{delta:g}".replace("e-0", "e-")


def tag_of(n, fam, k, delta, suffix=""):
    return f"D{n + 2}_{fam}_k{k}_d{dstr(delta)}{suffix}"


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
    kap = (n - 1.0) / 2.0 * (1.02 if mut("tang_kappa") else 1.0)
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
# criteria and constants (declared in negative-mass-raw I246-dh610/README.md before the runs of record)
# --------------------------------------------------------------------------------------------------------------------
DECLARED = dict(F1=1.5, K=28, W1=2.4, rcut=0.05, tol=1e-2, dxp=14.0, dm=3.0)       # #41's chain with K = 28 and dxp = 14 (limit e^{-f} = F1 2^K = 4.03e8; #41's 1 572 864 is read from the same series)
DELTAS = (1e-6, 1e-8)
WINDOWS = smg.WINDOWS                                                              # x+ - x+_f in [0, 1.5], [1.5, 3], [3, end]
BAND_P, BAND_R = sm.BAND_P, sm.BAND_R                                              # planted (1.75, 2.25), real series (1.5, 2.5)
BAR = dict(
    bracket=1e-11,          # RC1: hi - lo of p*_h (the floor of p*_h is 1e-11, #40; every delta is at least 1e3 times larger)
    rc3=2e-3,               # RC3: M at the first trapped cell against #40's M_hor
    rc4=2e-2,               # RC4: first law, relative
    rc7=1e-4,               # RC7: delta independence of S (plateau means)
    rc8_plat=1e-4, rc8_flux=1e-5,      # RC8: settled
    dep=0.1,                # DEP: the delta dependence of the two finest h agree to this fraction (or both are below rc7)
    machine=smg.BAR_MACHINE,           # #41's review: another machine reproduces the settled value to about 2e-5 (measured at D = 4)
    pt0=1e-9,               # PT0: closed form against dn.kruskal_exact
    pt1_rel=3e-4,           # PT1: planted M error at the third (finest real) h, relative to M_exact
)
PLANTED = {3: [dict(name="S", lam=0.54, kp=2.7, xh=5.0), dict(name="G", lam=0.73, kp=1.9, xh=5.0)],
           4: [dict(name="S", lam=0.60, kp=3.6, xh=5.0), dict(name="G", lam=1.19, kp=1.75, xh=5.0)]}
NOISE_RUNS = smg.NOISE_RUNS + ((1048576.0, 1),)                                    # #253's five and one at 2^20 x 2^-52 (the dose at which #253's scatter responded)
VOL = {2: math.pi, 3: 4.0 * math.pi / 3.0, 4: math.pi ** 2 / 2.0}                  # volumes of the unit balls B^2, B^3, B^4 ([ST17]'s "A_{n-2}" at Omega = 1)
ST17_EF = {4: 3.14, 5: 6.28, 6: 9.87}                                              # [ST17] Table I, GR rows, E_f/a_0 (read from the TeX source)
ST17_A4_PRINTED = math.pi ** 2                                                     # the printed A_4 = pi^2/Omega^4 (the unit 4-ball has pi^2/2)

results = []


def check(cid, text, ok, detail=""):
    ok = None if ok is None else bool(ok)
    results.append((cid, ok))
    print(f"[{'SKIP' if ok is None else ('pass' if ok else 'FAIL')}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def finish(t0):
    nfail = sum(1 for _, ok in results if ok is False)
    print(f"\n{len(results)} checks, {nfail} failed; {time.time() - t0:.0f} s")
    if MUT is not None:
        fired = sorted({c for c, ok in results if ok is False})
        want = sorted(MUTATIONS.get(MUT, []))
        print(f"mutation {MUT}: fired {fired}, declared {want}: {'MATCH' if fired == want else 'DIFFERENT'}")
        return 0 if fired == want else 1
    return 1 if nfail else 0


def cid(base, n, fam=None):
    return f"{base}-{n + 2}" + (fam or "")


# --------------------------------------------------------------------------------------------------------------------
# general-n first law and late-time summary
# --------------------------------------------------------------------------------------------------------------------
def flux_n(s, n):
    """dM/dx+ along the horizon, ghost + normal field: (r^n/n) e^f r_- (phi0_+^2 - phi1_+^2). Mutations: flux_sign, normal_drop, normal_sign, flux_exp."""
    sg = -1.0 if mut("flux_sign") else 1.0
    cn = 0.0 if mut("normal_drop") else (1.0 if mut("normal_sign") else -1.0)
    ex = n + 1 if mut("flux_exp") else n
    return sg * s["r"] ** ex / n / s["F"] * s["rm"] * (s["phip"] ** 2 + cn * s["phip1"] ** 2)


def flux_integral_n(s, n, i0, i1):
    f = flux_n(s, n)
    x = s["xp"][i0:i1 + 1]
    return float(np.sum(0.5 * (f[i0:i1] + f[i0 + 1:i1 + 1]) * np.diff(x)))


def summary_n(s, h, n):
    q = sm.settled_summary(s, h)
    v = s["xp"]
    q["flux05"] = float(np.abs(flux_n(s, n))[v > v[-1] - 0.5].mean())
    return q


def flux_planted():
    """General-n first-law integrand on planted horizon series with a known integral (r = 0.5, e^-f = 2, r_- = -0.4, phi_+ = 1, x+ in [0, 1]):
    (r^n/n) e^f r_- = -0.4 * 0.5^n / (2 n) per unit x+, for ghost only; the normal field has the opposite sign; equal fluxes cancel. n = 3, 4."""
    x = np.linspace(0.0, 1.0, 11)
    errs = []
    for n in NS:
        base = dict(xp=x, r=np.full(11, 0.5), F=np.full(11, 2.0), rm=np.full(11, -0.4))
        g0 = -0.4 * 0.5 ** n / (2.0 * n)
        for g, nm_, want in ((1.0, 0.0, g0), (0.0, 1.0, -g0), (1.0, 1.0, 0.0)):
            s = dict(base, phip=np.full(11, g), phip1=np.full(11, nm_))
            errs.append(abs(flux_integral_n(s, n, 0, 10) - want))
    check("FLX", "general-n first-law integrand on planted series: ghost only (r^n/n) e^f r_-, normal only the opposite sign, equal fluxes 0 (n = 3, 4)", max(errs) <= 1e-12,
          "errors " + ", ".join(f"{e:.1e}" for e in errs))


# --------------------------------------------------------------------------------------------------------------------
# runs: stored, or fresh under a mutation
# --------------------------------------------------------------------------------------------------------------------
def pstar_json(raw, n, fam, k):
    return json.load(open(os.path.join(raw, "runs", f"pstar_D{n + 2}_{fam}_k{k}.json")))


def get_run(raw, n, fam, k, delta, tag=None, noise=0.0, seed=0):
    """The stored series, or, under a mutation, a fresh run from the stored p*_h."""
    if MUT is None or tag is not None:
        return load_run(raw, tag or tag_of(n, fam, k, delta))
    pj = pstar_json(raw, n, fam, k)
    K = 0 if mut("no_zoom") else (12 if mut("short_limit") else None)
    ps = pj["ps"] + 2e-8 if mut("pstar_off") else pj["ps"]
    r = run_real_d(n, fam, pj["h"], delta, ps, K=K)
    return r["ser"], dict(status=r["status"], event=r["event"], point=r["point"], up=r["up"], Nc=r["Nc"], F1=r["F1"], K=r["K"], h=pj["h"], delta=delta, ps=ps)


# --------------------------------------------------------------------------------------------------------------------
# planted Tangherlini
# --------------------------------------------------------------------------------------------------------------------
def planted_hs(n):
    hr = tuple(th.CRIT["H"][n])
    return (2.0 * hr[0],) + hr                          # one coarser than the first real h, then the three real ones


def pt0_check():
    """PT0: the closed-form exact Tangherlini against dn.kruskal_exact (#39's, a different code: spline of the regular R(r), bisection) at points with
    r in [0.4, 5.9], and e^{-f} on the horizon of the stretching gauge = e^{kp (x+ - xs)}."""
    rng = np.random.default_rng(7)
    worst, worst_h = 0.0, 0.0
    for n in NS:
        U = -rng.uniform(0.05, 0.9, 200)
        V = rng.uniform(0.02, 1.5, 200)
        r_dn, f_dn = dn.kruskal_exact(n, U, V)
        m = (r_dn > 0.4) & (r_dn < 5.9)
        r_c, emf_c = tang_exact_unit(n, -U * V)
        worst = max(worst, float(np.max(np.abs(r_c[m] - r_dn[m]))), float(np.max(np.abs(-np.log(emf_c[m]) - f_dn[m]))))
        for c in PLANTED[n]:
            xp = np.linspace(0.0, 9.0, 9)
            r_h, f_h = tang_exact(n, c["lam"], c["kp"], c["xh"], xp, c["xh"] + 0.0 * xp, 3.0)
            worst_h = max(worst_h, float(np.max(np.abs(np.exp(-f_h) - np.exp(c["kp"] * (xp - 3.0))) / np.exp(c["kp"] * (xp - 3.0)))), float(np.max(np.abs(r_h - c["lam"]))))
    check("PT0", "the planted exact Tangherlini (closed form, n = 3, 4) equals dn.kruskal_exact's r and f at 400 points (r in [0.4, 5.9]); on the horizon of the stretching gauge r = lam and e^{-f} = e^{kp (x+ - xs)}",
          worst <= BAR["pt0"] and worst_h <= BAR["pt0"], f"max |r|, |f| difference {worst:.1e}; horizon {worst_h:.1e} (bar {BAR['pt0']:.0e})")


def planted_checks(log=print, ns=NS, cases=None):
    """PT1: the exact Tangherlini through the chain at h = 2 h_1, h_1, h_2, h_3 of #40's set, two planted cases per D (the real S and G runs' horizon radius and
    stretching): M = lam^(n-1)/2 recovered with error <= pt1_rel M at h_3, observed order in BAND_P at the three halvings, limit reached."""
    F1, K = CRITD["F1"], CRITD["K"]
    if mut("no_zoom"):
        K = 0
    for n in ns:
        hs = planted_hs(n)
        for c in PLANTED[n]:
            if cases is not None and (n, c["name"]) not in cases:
                continue
            eM, st, Fend = [], [], []
            for h in hs:
                out, e = tang_planted(n, h, c["lam"], c["kp"], c["xh"], nc_of(h), K=K)
                s = out["ser"]
                eM.append(float(np.max(e["e_M"])))
                st.append(out["status"])
                Fend.append(float(s["F"][-1]) if len(s["F"]) else 0.0)
                log(f"    PT1 D = {n + 2}, {c['name']}-like (lam {c['lam']}, kp {c['kp']}) h = {h:g}: status {out['status']}, columns {len(s['M'])}, F to {Fend[-1]:.4g}, "
                    f"max |M - {e['M_exact']:.5f}| {eM[-1]:.3e}, max |x- - x_h| {np.max(e['e_x']):.2e}")
            o = [sm.order(eM[i], eM[i + 1]) for i in range(len(hs) - 1)]
            Mex = c["lam"] ** (n - 1) / 2.0
            check(f"PT1-{n + 2}{c['name']}", f"exact Tangherlini D = {n + 2} ({c['name']}-like, lam {c['lam']}, kp {c['kp']}) through the chain: M recovered, second order, limit reached",
                  eM[3] <= BAR["pt1_rel"] * Mex and all(BAND_P[0] <= x <= BAND_P[1] for x in o) and all(x == "limit" for x in st) and all(f >= 0.5 * F1 * 2.0 ** CRITD["K"] for f in Fend),
                  f"max |M - M_exact| at h = {hs[3]:g}: {eM[3]:.2e} (bar {BAR['pt1_rel'] * Mex:.1e}); errors {['%.2e' % x for x in eM]}; orders {[round(x, 3) for x in o]} in {BAND_P}; status {st}; F reached {[round(f) for f in Fend]}")


def level0_identity(n, fam, raw):
    """L0ID: level 0 of the chain (forced = [], no zoom) against dn.evolve on the same data (bit-identical nodes and dead flags, same first trapped cell), k = 0, delta = 1e-6."""
    pj = pstar_json(raw, n, fam, 0)
    r = run_real_d(n, fam, pj["h"], 1e-6, pj["ps"], forced=[])
    d, L0 = r["data"], r["levels"][0]
    res = dn.evolve(d, store=True, stop=False, rcut=CRITD["rcut"], fcut=1e300, tol=CRITD["tol"])
    same_G = np.array_equal(res.G[:, :L0.Np + 1, :L0.Nm + 1], L0.G)
    same_d = np.array_equal(res.dead[:L0.Np + 1, :L0.Nm + 1], L0.dead)
    ev = res.events.get("collapse")
    check(cid("L0ID", n, fam), f"D = {n + 2}, family {fam}: level 0 of the chain equals dn.evolve bit for bit (nodes, dead flags) and finds the same first trapped cell",
          same_G and same_d and ev == r["event"], f"h = {pj['h']:g}, delta = 1e-6: nodes equal {same_G}, dead equal {same_d}, events {ev} and {r['event']}")
    del res, r, d, L0


# --------------------------------------------------------------------------------------------------------------------
# analysis of one (D, family)
# --------------------------------------------------------------------------------------------------------------------
def nks(n, fam):
    return 3


def at_limit(s, Fmax):
    """The series up to the last column with e^{-f} < Fmax (the run's own convention: the stop column belongs to the next level). A run to a larger limit
    contains the run to Fmax as a bit-identical prefix (check PFX)."""
    over = s["F"] >= Fmax
    idx = int(np.argmax(over)) if over.any() else len(s["F"])
    return {k: v[:idx] for k, v in s.items()}


def limit_checks(n, fam, S, I, hs, ks, Fmax, sfx, label):
    """RC5, RC6, RC8 at one limit (sfx '' = this driver's K = 28; 'a' = #41's 1.57e6, read from the same series). Returns the late-time summaries."""
    D = n + 2
    ok5, det = True, []
    for key, s in S.items():
        fine = np.allclose(np.diff(s["xp"]), hs[key[0]], atol=1e-9)
        okk = 0.5 * Fmax <= s["F"][-1] < Fmax and s["edge"].min() >= I[key]["Nc"] / 8 and fine and I[key]["status"] == "limit"
        ok5 &= okk
        det.append(f"k={key[0]},d={key[1]:g}: F_end {s['F'][-1]:.3g}, edge {s['edge'].min():.0f}/{I[key]['Nc']}")
    check(cid("RC5" + sfx, n, fam), f"D = {D}, {fam}, limit {label}: every run ends at the limit (e^-f on the horizon, {0.5 * Fmax:.4g} <= F_end < {Fmax:.4g}), horizon at least Nc/8 cells from the box edges, no gaps", ok5, "; ".join(det))
    print(f"\n== D = {D}, {fam}, limit {label}: second-order convergence of M_AH(x+): max |M_h - M_h/2| over windows of x+ - x+_f (cubic splines), h = {tuple(hs[k] for k in ks)}")
    ok6 = True
    for dl in DELTAS:
        rows = sm.band_orders({hs[k]: S[(k, dl)] for k in ks}, tuple(hs[k] for k in ks), WINDOWS, dl)
        for a_, b_, Ee, o in rows:
            good = all(BAND_R[0] <= x <= BAND_R[1] for x in o)
            ok6 &= good
            print(f"    RC6{sfx} D = {D}, {fam}, delta = {dl:g}, x+ - x+_f in [{a_:g}, {b_:g}]: E12 {Ee[0]:.2e}, E23 {Ee[1]:.2e}, order {o[0]:.2f} {'ok' if good else 'OUT'}")
    check(cid("RC6" + sfx, n, fam), f"D = {D}, {fam}, limit {label}: M_AH(x+) converges at second order at the three h {tuple(hs[k] for k in ks)}, three windows, both deltas (observed order in {BAND_R})", ok6)
    dr = {k: summary_n(s, hs[k[0]], n) for k, s in S.items()}
    check(cid("RC8" + sfx, n, fam), f"D = {D}, {fam}, limit {label}: M_AH has settled at the limit: plateau mean over the last 1.0 of x+ differs from the previous 1.0 by <= {BAR['rc8_plat']:g} and the mean horizon flux (both fields) over the last 0.5 is <= {BAR['rc8_flux']:g}, every run",
          all(abs(q["dplat"]) <= BAR["rc8_plat"] and q["flux05"] <= BAR["rc8_flux"] for q in dr.values()),
          "plateau change " + ", ".join(f"{q['dplat']:+.1e}" for q in dr.values()) + "; flux " + ", ".join(f"{q['flux05']:.1e}" for q in dr.values()))
    return dr


def analyze_case(raw, n, fam):
    """RC1 to RC4, RC5 to RC8 at both limits, RC7 or DEP, TRK for one (n, fam): the three h of the record, both delta. Returns the data for the reports."""
    D = n + 2
    ks = list(range(nks(n, fam)))
    hs = tuple(th.CRIT["H"][n])
    F28 = CRITD["F1"] * 2.0 ** CRITD["K"]
    F20 = CRITD["F1"] * 2.0 ** K41
    X = th.CRIT["X"][fam][n]
    ok1, det1 = True, []
    ref = {}
    for k in ks:
        pj = pstar_json(raw, n, fam, k)
        ref[k] = i40_job(n, fam, k) if k < n_h(n, fam) else None
        okk = (pj["hi"] - pj["lo"] <= BAR["bracket"] and pj["h"] == hs[k] and (ref[k] is None or (float(ref[k]["h_log"]) == hs[k] and ref[k]["h_summary"] == hs[k])))
        pts = []
        for sgn in (-1, 1):
            q = th.point(n, fam, hs[k], X, pj["ps"] + sgn * 1e-10)
            pts.append(q["fate"])
        okk &= pts == [-1, 1]
        ok1 &= okk
        det1.append(f"k={k} h={hs[k]:g}: hi-lo {pj['hi'] - pj['lo']:.1e}{' (' + str(pj['stopped']) + ')' if pj['stopped'] else ''}, fates {pts}")
    check(cid("RC1", n, fam), f"D = {D}, {fam}: p*_h (#40's, or bisected here for G at D = 6, h = 0.00125) with hi - lo <= {BAR['bracket']:g}, h as the job takes it, fate(p*_h -/+ 1e-10) = collapse/expansion on this machine", ok1, "; ".join(det1))
    S, I, fail = {}, {}, {}
    for dl in DELTAS:
        for k in ks:
            try:
                S[(k, dl)], I[(k, dl)] = get_run(raw, n, fam, k, dl)
                if len(S[(k, dl)]["r"]) == 0:
                    raise RuntimeError("empty series: " + str(I[(k, dl)]["status"]))
            except Exception as e:
                fail[(k, dl)] = repr(e)
                print(f"    run D = {D}, {fam}, k = {k}, delta = {dl:g} failed: {e!r}")
    if fail:
        for b in ("RC2", "RC3", "RC4", "RC5", "RC6", "RC8", "RC5a", "RC6a", "RC8a", "RC7" if fam == "S" else "DEP"):
            check(cid(b, n, fam), "a run failed", False, "; ".join(f"{kk}: {v}" for kk, v in fail.items()))
        return None
    kref = [k for k in ks if ref[k] is not None]
    ev_ok = all(list(I[(k, 1e-8)]["event"]) == list(ref[k]["event"]) for k in kref)
    t_ok = all(abs(sm.t_of((I[(k, 1e-8)]["event"][1] + 0.5) * hs[k], (I[(k, 1e-8)]["event"][0] + 0.5) * hs[k]) - ref[k]["T"]) <= SQ2 * hs[k] for k in kref)
    fate_ok = all(I[(k, dl)]["point"]["fate"] == -1 and I[(k, dl)]["up"]["fate"] == 1 for k in ks for dl in DELTAS)
    check(cid("RC2", n, fam), f"D = {D}, {fam}: fates at p*_h -/+ delta (all h, both delta); first trapped cell and its time against #40's at delta = 1e-8 (where #40 has the job)", fate_ok and ev_ok and t_ok,
          "events " + str({k: list(I[(k, 1e-8)]["event"]) for k in ks}) + " against #40's " + str({k: ref[k]["event"] for k in kref}))
    dm = [abs(S[(k, 1e-8)]["M"][0] - ref[k]["M_hor"]) for k in kref]
    check(cid("RC3", n, fam), f"D = {D}, {fam}: M_AH at the first trapped cell against #40's formation value at p*_h - 1e-8", max(dm) <= BAR["rc3"],
          "|dM| " + ", ".join(f"{x:.1e}" for x in dm) + f" (bar {BAR['rc3']:g}); M_AH {[round(float(S[(k, 1e-8)]['M'][0]), 5) for k in ks]}, #40 {[round(ref[k]['M_hor'], 5) for k in kref]}")
    worst, parts = 0.0, []
    for key, s in S.items():
        i1 = int(np.argmax(s["F"] >= 20.0))
        dMs = s["M"][i1] - s["M"][0]
        rel = abs(dMs - flux_integral_n(s, n, 0, i1)) / abs(dMs)
        worst = max(worst, rel)
        parts.append(f"k={key[0]},d={key[1]:g}: dM {dMs:+.4e}, flux {flux_integral_n(s, n, 0, i1):+.4e}")
    check(cid("RC4", n, fam), f"D = {D}, {fam}: first law on the horizon, both fields: M(x+_b) - M(x+_a) = integral of (r^n/n) e^f r_- (phi0_+^2 - phi1_+^2) along it, first trapped cell to F = 20", worst <= BAR["rc4"],
          f"max relative difference {worst:.2e} over all runs (bar {BAR['rc4']:g})")
    for p_ in parts:
        print("    " + p_)
    S20 = {key: at_limit(s, F20) for key, s in S.items()}
    dr = limit_checks(n, fam, S, I, hs, ks, F28, "", f"e^-f = {F28:.3g} (K = {CRITD['K']})")
    dr20 = limit_checks(n, fam, S20, I, hs, ks, F20, "a", f"e^-f = {F20:.3g} (#41's)")
    kf = ks[-1]
    dd = [abs(dr[(k, 1e-6)]["plat1"] - dr[(k, 1e-8)]["plat1"]) for k in ks]
    ordd = [math.log2(dd[i] / dd[i + 1]) if dd[i] > 0 and dd[i + 1] > 0 else float("nan") for i in range(len(ks) - 1)]
    print(f"    delta difference |plat(1e-6) - plat(1e-8)| at the three h {[f'{x:.2e}' for x in dd]}, observed orders {[round(x, 2) for x in ordd]}")
    if fam == "S":
        check(cid("RC7", n, fam), f"D = {D}, S: the settled value does not depend on delta: plateau means over the last 1.0 at 1e-6 and 1e-8 agree at the finest h", dd[-1] <= BAR["rc7"],
              "|dM| " + ", ".join(f"{x:.1e}" for x in dd) + f" at the three h, observed orders {[round(x, 2) for x in ordd]} (bar {BAR['rc7']:g} at the finest h)")
    else:
        d_f, d_p = dd[-1], dd[-2]
        check(cid("DEP", n, fam), f"D = {D}, G: the dependence on delta, plateau(1e-8) - plateau(1e-6), is the same at the two finest h to {BAR['dep']:.0%} (a property of the solution, not of the resolution), or both are below {BAR['rc7']:g}",
              (abs(d_f - d_p) <= BAR["dep"] * d_f) or (max(d_f, d_p) <= BAR["rc7"]),
              "|dM| " + ", ".join(f"{x:.2e}" for x in dd) + f" at the three h; signed at the finest h {dr[(kf, 1e-8)]['plat1'] - dr[(kf, 1e-6)]['plat1']:+.2e}")
    if fam == "G":
        trk = min(float(np.max(np.abs(s["phip1"]))) for s in S.values())
        wor = max(float(np.max(np.abs(s["phip0_re"] - s["phip"]))) / max(float(np.max(np.abs(s["phip"]))), 1e-300) for s in S.values())
        check(cid("TRK", n, fam), f"D = {D}, G: the track wrapper's recomputed phi_+ of the ghost equals settled_mass.track's (relative 1e-12) and the normal field's phi_+ is nonzero on the horizon", wor <= 1e-12 and trk > 0,
              f"max relative difference {wor:.1e}, smallest max|phi_+(normal)| {trk:.2e}")
    return dict(S=S, I=I, dr=dr, dr20=dr20, S20=S20, hs=hs, ks=ks, ref=ref)


def prefix_check(raw):
    """PFX: a run to a smaller limit is a bit-identical prefix of the run to K = 28 (the case D = 5, 6 x S, G, k = 0, delta = 1e-8): the K = 20 and K = 24 runs
    (same dxp) against the K = 28 series cut at their limits (M, x+, x-, r, e^-f equal bit for bit), and the end values M_end of the three."""
    rows, ok, found = [], True, 0
    for n in NS:
        for fam in FAMS:
            try:
                s28, _ = load_run(raw, tag_of(n, fam, 0, 1e-8))
            except Exception:
                continue
            ends = [("28", float(s28["M"][-1]))]
            for K in (20, 24):
                try:
                    sk, _ = load_run(raw, tag_of(n, fam, 0, 1e-8, f"_K{K}"))
                except Exception:
                    ok = False
                    rows.append(f"D = {n + 2} {fam} K = {K}: run missing")
                    continue
                found += 1
                cut = at_limit(s28, CRITD["F1"] * 2.0 ** K)
                same = all(np.array_equal(cut[q], sk[q]) for q in ("M", "xp", "xm", "r", "F"))
                ok &= same
                ends.append((str(K), float(sk["M"][-1])))
                rows.append(f"D = {n + 2} {fam} K = {K}: {len(sk['M'])} columns, bit-identical to the K = 28 series cut at F < {CRITD['F1'] * 2.0 ** K:.4g}: {same}")
            print(f"    M_end at the limits K = 20, 24, 28: D = {n + 2}, {fam}: " + ", ".join(f"K = {K}: {m:.6f}" for K, m in sorted(ends)))
    check("PFX", "a run to a smaller limit (K = 20, 24) equals the K = 28 series cut at its limit bit for bit (M, x+, x-, r, e^-f), D = 5, 6, S and G, h = 0.00707 / 0.005, delta = 1e-8", ok and found > 0, "; ".join(rows))


# --------------------------------------------------------------------------------------------------------------------
# G against S, [ST17], reports
# --------------------------------------------------------------------------------------------------------------------
def c_of(D, a4_printed=False):
    """[ST17]'s energy per unit Misner-Sharp mass: E = (D - 2) V_{D-2} M with M = r^(D-3) (1 - (grad r)^2)/2 and kappa^2 = 8 pi G = 1 (both here and in [ST17]);
    V is the volume of the unit ball, as the printed A_2, A_3, A_5 are; the printed A_4 = pi^2 is twice the unit 4-ball's pi^2/2."""
    V = dict(VOL)
    if a4_printed:
        V[4] = ST17_A4_PRINTED
    if mut("conv_vol"):
        V[4] = ST17_A4_PRINTED
    return (D - 2.0) * V[D - 2]


def conversion_checks():
    """CONV: the conversion reproduces [ST17]'s static-throat energy (their E_f of every GR row equals it) in D = 4, 5, 6: c_D / 2 against 3.14, 6.28, 9.87.
    A4P: the printed A_4 = pi^2 does not (it gives 2 pi^2 = 19.7, twice the table's 9.87)."""
    rows = [(D, c_of(D) / 2.0, ST17_EF[D]) for D in (4, 5, 6)]
    check("CONV", "the conversion E = (D - 2) V_{D-2} M (kappa^2 = 1, unit-ball volumes) gives the static throat (M = 1/2) E = pi, 2 pi, pi^2 = [ST17]'s E_f/a_0 = 3.14, 6.28, 9.87 to their digits, D = 4, 5, 6",
          all(abs(e - t) <= 0.005 for _, e, t in rows), "; ".join(f"D = {D}: {e:.4f} against {t}" for D, e, t in rows))
    e6p = c_of(6, a4_printed=True) / 2.0
    check("A4P", "[ST17]'s printed A_4 = pi^2/Omega^4 gives E_throat = 2 pi^2 = 19.74 at n = 6, twice the table's 9.87 (a discrepancy inside the paper)", abs(e6p / ST17_EF[6] - 2.0) <= 0.01,
          f"E_throat from the printed A_4: {e6p:.4f}; table {ST17_EF[6]}; ratio {e6p / ST17_EF[6]:.4f}")


def error_bars(A):
    """Kinds of error bar of a case's settled value at delta = 1e-8, finest h: h (finest pair), delta (1e-6 against 1e-8), machine (declared)."""
    dr, ks = A["dr"], A["ks"]
    kf, kp = ks[-1], ks[-2]
    pl = dr[(kf, 1e-8)]["plat1"]
    Bh = abs(dr[(kp, 1e-8)]["plat1"] - pl)
    Bd = abs(dr[(kf, 1e-6)]["plat1"] - pl)
    rich = pl - (dr[(kp, 1e-8)]["plat1"] - pl) / 3.0
    spread = max(dr[(k, 1e-8)]["plat1"] for k in ks) - min(dr[(k, 1e-8)]["plat1"] for k in ks)
    return dict(plat=pl, Bh=Bh, Bd=Bd, rich=rich, spread=spread)


def compare_gs(n, AS, AG):
    D = n + 2
    eS, eG = error_bars(AS), error_bars(AG)
    bar = smg.cmp_bar(eG["Bh"], eG["Bd"], BAR["machine"])
    dM = eG["plat"] - eS["plat"]
    verdict = smg.cmp_verdict(dM, bar)
    print(f"\n== D = {D}: the settled mass of G against S (error bars named, kept separate), h = {AG['hs'][-1]:g}, delta = 1e-8")
    print(f"  S: plateau mean over the last 1.0 of x+ {eS['plat']:.6f} (r_AH = {(2 * eS['plat']) ** (1.0 / (n - 1)):.5f}); h: finest pair {eS['Bh']:.2e}; delta: |plat(1e-6) - plat(1e-8)| {eS['Bd']:.2e}; Richardson {eS['rich']:.6f}; spread over the three h {eS['spread']:.2e}")
    print(f"  G: plateau mean {eG['plat']:.6f} (r_AH = {(2 * eG['plat']) ** (1.0 / (n - 1)):.5f}); h: finest pair {eG['Bh']:.2e}; delta: |plat(1e-6) - plat(1e-8)| {eG['Bd']:.2e}; Richardson {eG['rich']:.6f}; spread over the three h {eG['spread']:.2e}")
    print(f"  bar (declared as #253's): h {eG['Bh']:.2e} + delta {eG['Bd']:.2e} + machine {BAR['machine']:.0e} = {bar:.2e}; G - S = {dM:+.3e} (Richardson {eG['rich'] - eS['rich']:+.3e}); G/S = {eG['plat'] / eS['plat']:.4f}")
    d20 = {fm: A_["dr20"][(A_["ks"][-1], 1e-8)]["plat1"] for fm, A_ in (("S", AS), ("G", AG))}
    print(f"  at #41's limit (e^-f = 1.57e6, same series cut): S {d20['S']:.6f}, G {d20['G']:.6f}, G - S {d20['G'] - d20['S']:+.3e}")
    print(f"  VERDICT D = {D}: G's settled mass {'equals' if verdict else 'DIFFERS FROM'} S's: |G - S| = {abs(dM):.2e} {'<=' if verdict else '>'} bar {bar:.2e}")
    pb = 3e-5
    check(cid("CMPT", n), f"D = {D}: the comparison on planted offsets: half the bar is 'equal', three times the bar is 'differs' (the comparison can fail)",
          smg.cmp_verdict(0.5 * smg.cmp_bar(pb, pb), smg.cmp_bar(pb, pb)) and not smg.cmp_verdict(3.0 * smg.cmp_bar(pb, pb), smg.cmp_bar(pb, pb)), f"planted bar {smg.cmp_bar(pb, pb):.2e}")
    return dict(S=eS, G=eG, dM=dM, bar=bar, verdict=verdict)


def print_series(n, fam, A):
    """M_AH(x+ - x+_f), both delta, finest h: r_AH, M, E in [ST17]'s unit, e^{-f}."""
    D = n + 2
    kf = A["ks"][-1]
    hf = A["hs"][kf]
    for dl in DELTAS:
        s = A["S"][(kf, dl)]
        tf = sm.t_of(s["xp"][0], s["xm"][0])
        print(f"\n== M_AH(x+) at D = {D}, family {fam}, delta = {dl:g}, h = {hf:g}: x+_f = {s['xp'][0]:.4f}, x-_f = {s['xm'][0]:.4f}, t_f = {tf:.4f}; M = r_AH^{n - 1}/2; a = 1, 8 pi G = 1; t = (x+ + x-)/sqrt2")
        print("    dx+     t - t_f     r_AH      M_AH      2 M       e^-f")
        for dv in (0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.4):
            i = int(np.argmin(np.abs(s["xp"] - (s["xp"][0] + dv))))
            if abs(s["xp"][i] - s["xp"][0] - dv) > 2 * hf:
                break
            print(f"    {dv:5.2f} {sm.t_of(s['xp'][i], s['xm'][i]) - tf:9.4f}  {s['r'][i]:8.5f}  {s['M'][i]:8.5f}  {2 * s['M'][i]:8.5f}  {s['F'][i]:10.4g}")


def late_table(n, fam, A):
    for key, lab in (("dr", f"limit e^-f = {CRITD['F1'] * 2.0 ** CRITD['K']:.3g} (K = {CRITD['K']}, the quoted value)"), ("dr20", f"limit e^-f = {CRITD['F1'] * 2.0 ** K41:.3g} (#41's, same series cut)")):
        print(f"\n== Late-time values, D = {n + 2}, {fam}, {lab} (reports)")
        print("  delta   h        M_end     plat(1.0)  range(1.0) plat(0.5)  drift(1.0)  plat change  flux(0.5)")
        for dl in DELTAS:
            for k in A["ks"]:
                q = A[key][(k, dl)]
                print(f"  {dl:6g} {A['hs'][k]:9g}  {q['M_end']:.6f}  {q['plat1']:.6f}  {q['rng1']:.1e}  {q['plat0.5']:.6f}  {q['drift1']:+.2e}  {q['dplat']:+.2e}  {q['flux05']:.2e}")
        for dl in DELTAS:
            p = [A[key][(k, dl)]["plat1"] for k in A["ks"]]
            a = [A[key][(k, dl)]["M_end"] for k in A["ks"]]
            print(f"  delta = {dl:g}: plateau(1.0) {[round(x, 6) for x in p]}, differences {p[0] - p[1]:+.2e}, {p[1] - p[2]:+.2e}, Richardson (second order) {p[2] - (p[1] - p[2]) / 3:.6f}; "
                  f"M_end {[round(x, 6) for x in a]}, Richardson {a[2] - (a[1] - a[2]) / 3:.6f}")


def st17_table(cases):
    """The settled values against [ST17]'s rows in one unit: E = c_D M (c_D = (D - 2) V_{D-2}, kappa^2 = 1), E_throat = c_D / 2; their E_f/a_0 and E_f/E_throat."""
    print("\n== The comparison with [ST17] (Table I, GR rows, E_f/a_0 = 3.14, 6.28, 9.87 for n = 4, 5, 6) in one unit")
    print("  unit: E = c_D M, c_D = (D - 2) V_{D-2} (V_2 = pi, V_3 = 4 pi/3, V_4 = pi^2/2, unit-ball volumes as [ST17]'s printed A_2, A_3, A_5; kappa^2 = 1), M = r^(D-3)(1 - (grad r)^2)/2 with a_0 = 1")
    print(f"  c_4 = {c_of(4):.4f}, c_5 = {c_of(5):.4f}, c_6 = {c_of(6):.4f}; static throat E = c_D/2 = {c_of(4) / 2:.4f}, {c_of(5) / 2:.4f}, {c_of(6) / 2:.4f}; [ST17]'s printed A_4 would give c_6 = {c_of(6, True):.4f}")
    print("  D  quantity                         M          E = c_D M    E/E_throat   [ST17] E_f   [ST17] E_f/E_throat   here / [ST17]")
    for D in (5, 6):
        n = D - 2
        for fam in FAMS:
            if (n, fam) not in cases:
                continue
            A = cases[(n, fam)]
            kf = A["ks"][-1]
            s = A["S"][(kf, 1e-8)]
            for lab, M in (("formation (first trapped cell)", float(s["M"][0])), ("settled (plateau, delta = 1e-8)", error_bars(A)["plat"])):
                E = c_of(D) * M
                print(f"  {D}  {fam} {lab:31s} {M:9.5f}  {E:11.4f}  {E / (c_of(D) / 2.0):10.4f}   {ST17_EF[D]:9.2f}   {ST17_EF[D] / (c_of(D) / 2.0):12.4f}   {E / ST17_EF[D]:10.4f}")
    print("  D = 4 (#41, #253): S settles at M = 0.2231 (E/E_throat 0.4462), G at 0.2355; [ST17]'s 3.14 is E_throat")


def noise_report(raw, cases_data, ns_fams, t0):
    """The ulp-noise runs: how far M_AH(x+) moves under relative noise in the ray data (delta = 1e-8), against the h-differences (as #253)."""
    for (n, fam) in ns_fams:
        D = n + 2
        A = cases_data.get((n, fam))
        if A is None:
            continue
        ks_n = A["ks"][1:]
        have = all(os.path.exists(os.path.join(raw, "runs", f"D{D}_{fam}_k{k}_d1e-8_rn_u{u:g}_s{sd}", "series.npz")) for k in ks_n for u, sd in NOISE_RUNS)
        print(f"\n== Roundoff floor, D = {D}, {fam}: M_AH(x+) under relative noise in the ray data (delta = 1e-8) against the h-differences")
        if not have:
            check(cid("RN0", n, fam), "noise runs present", False, "missing runs")
            check(cid("RNR", n, fam), "noise runs present", False, "missing runs")
            continue
        N = {}
        for k in ks_n:
            for u, sd in NOISE_RUNS:
                N[(k, u, sd)], _ = load_run(raw, f"D{D}_{fam}_k{k}_d1e-8_rn_u{u:g}_s{sd}")
        ref = {k: (A["S"][(k, 1e-6)] if mut("rn_ref") else A["S"][(k, 1e-8)]) for k in ks_n}
        z = [np.array_equal(N[(k, 0.0, 1)]["M"], ref[k]["M"]) and np.array_equal(N[(k, 0.0, 1)]["xp"], ref[k]["xp"]) for k in ks_n]
        check(cid("RN0", n, fam), f"D = {D}, {fam}: noise 0 reproduces the run of record bit for bit (the injection is inert at zero, the job is deterministic)", all(z), f"k = {ks_n}: {z}")
        sc = {}
        for k in ks_n:
            for u in sorted({u for u, _ in NOISE_RUNS if u > 0}):
                vals = [smg.window_scatter(N[(k, u, sd)], ref[k], WINDOWS) for uu, sd in NOISE_RUNS if uu == u]
                sc[(k, u)] = [max(v[w] for v in vals) for w in range(len(WINDOWS))]
                print(f"    h = {A['hs'][k]:g}, noise {u:g} x 2^-52, {len(vals)} seed(s): max |M_noisy - M| in windows {WINDOWS}: " + ", ".join(f"{x:.2e}" for x in sc[(k, u)]))
        big = max(u for u, _ in NOISE_RUNS)
        resp = [all(sc[(k, 1.0)][w] > 0.0 and sc[(k, big)][w] >= 10.0 * sc[(k, 1.0)][w] for w in range(len(WINDOWS))) for k in ks_n]
        check(cid("RNR", n, fam), f"D = {D}, {fam}: the scatter responds when the noise is large: nonzero at 2^-52 and at {big:g} x 2^-52 at least 10 times that, every window and h", all(resp),
              "; ".join(f"h = {A['hs'][k]:g}: " + ", ".join(f"{sc[(k, big)][w] / max(sc[(k, 1.0)][w], 1e-300):.1f}" for w in range(len(WINDOWS))) for k in ks_n))
        print("  scatter at 2^-52 against the h-difference E between the h and its coarser neighbour (same windows, delta = 1e-8):")
        rows = sm.band_orders({A["hs"][k]: A["S"][(k, 1e-8)] for k in A["ks"]}, tuple(A["hs"][k] for k in A["ks"]), WINDOWS, 1e-8)
        for kk, k in enumerate(ks_n):
            for w, (a_, b_, Ee, o) in enumerate(rows):
                ratio = sc[(k, 1.0)][w] / Ee[kk]
                word = "reached" if ratio >= 1.0 else ("approached" if ratio >= 0.1 else "not reached")
                print(f"    D = {D}, {fam}, h = {A['hs'][k]:g}, x+ - x+_f in [{a_:g}, {b_:g}]: noise scatter {sc[(k, 1.0)][w]:.2e}, E({A['hs'][k - 1]:g} vs {A['hs'][k]:g}) {Ee[kk]:.2e}, ratio {ratio:.3f}: floor {word}")


# --------------------------------------------------------------------------------------------------------------------
# the analysis
# --------------------------------------------------------------------------------------------------------------------
MUTATIONS = {}
CASES_ALL = [(n, f) for n in NS for f in FAMS]


def analyze_cmd(raw, cases=None, planted_only=False, skip=()):
    t0 = time.time()
    cases = CASES_ALL if cases is None else cases
    print("\n== Declared chain (#41's constants), the exact planted solution, and the planted cases")
    bad = {k: (CRITD[k], v) for k, v in DECLARED.items() if CRITD[k] != v}
    xs = {(f, n): th.CRIT["X"][f][n] for f in FAMS for n in NS}
    xdecl = dict(S={3: 18.0, 4: 15.0}, G={3: 22.0, 4: 28.0})
    xbad = {k: v for k, v in xs.items() if v != xdecl[k[0]][k[1]]}
    check("CFG", "the chain constants of this driver equal the declared (#41's) and #40's rectangles X of S and G at D = 5, 6 equal the declared", not bad and not xbad,
          f"differs: {bad} {xbad}" if bad or xbad else f"F1 {CRITD['F1']}, K {CRITD['K']}, W1 {CRITD['W1']}, dxp {CRITD['dxp']}, dm {CRITD['dm']}; X {xs}")
    flux_planted()
    pt0_check()
    planted_checks(cases={(n, f) for n, f in cases} if cases != CASES_ALL else None)
    conversion_checks()
    if planted_only:
        return finish(t0)
    if "pl41" not in skip:
        try:
            smg.planted_through_driver()
        except Exception as e:
            check("PL1", "#41's PL1 through the imported chain", False, "failed: " + repr(e))
    for (n, fam) in cases:
        try:
            level0_identity(n, fam, raw)
        except Exception as e:
            check(cid("L0ID", n, fam), "level 0 against dn.evolve", False, "failed: " + repr(e))
    if "pfx" not in skip:
        print("\n== Runs to a smaller limit against the run to K = 28")
        prefix_check(raw)
    data = {}
    for (n, fam) in cases:
        print(f"\n== D = {n + 2}, family {fam}: p*_h, runs and checks")
        data[(n, fam)] = analyze_case(raw, n, fam)
    for (n, fam), A in data.items():
        if A is not None:
            print_series(n, fam, A)
            late_table(n, fam, A)
    for n in NS:
        if data.get((n, "S")) is not None and data.get((n, "G")) is not None:
            compare_gs(n, data[(n, "S")], data[(n, "G")])
    st17_table({k: v for k, v in data.items() if v is not None})
    noise_report(raw, data, cases, t0)
    return finish(t0)


# --------------------------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------------------------
def parse_cases(txt):
    """'5S,5G,6S' -> [(3, 'S'), (3, 'G'), (4, 'S')] (D and family)."""
    return [(int(c[:-1]) - 2, c[-1]) for c in txt.split(",")]


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
        if MUT == "family_B":
            th.CRIT["G_B"] = 0.11
        if MUT == "crit_dm":
            CRITD["dm"] = 2.0
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
    elif cmd == "planted":
        t0 = time.time()
        flux_planted()
        pt0_check()
        planted_checks(cases=parse_cases(opts["cases"]) if "cases" in opts else None)
        conversion_checks()
        sys.exit(finish(t0))
    elif cmd == "analyze":
        sys.exit(analyze_cmd(opts["raw"], cases=parse_cases(opts["cases"]) if "cases" in opts else None, skip=tuple(opts.get("skip", "").split(","))))
    elif cmd == "tang_test":
        tang_selftest()
    else:
        raise SystemExit(__doc__)


def tang_selftest():
    """Development: the closed-form Tangherlini against dn.kruskal_exact (r < 6)."""
    for n in NS:
        spl, g, Fx = dn.kruskal_tools(n)
        kap = (n - 1.0) / 2.0
        r = np.linspace(0.4, 5.9, 12)
        Yd = (r - 1.0) * np.exp(2 * kap * spl(r))
        print(f"n = {n}: max relative |Y closed form - Y of dn.kruskal_tools| over r in [0.4, 5.9]: {np.max(np.abs(tang_Y(n, r) - Yd) / np.maximum(1e-300, np.abs(Yd))):.2e}")


if __name__ == "__main__":
    main()

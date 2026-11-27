#!/usr/bin/env python3
"""Settled black-hole mass of the one-sided family G (#253): #41's apparent-horizon chain (settled_mass.py) for #40's family G at D = 4.

A thin driver. It imports settled_mass.py (#41), threshold.py (#40) and dn.py (#39) unchanged and adds only what family G needs:
  - the family: ghost pulse A = p*_h - delta plus the fixed normal pulse B = 0.1 on the same ray (threshold.family("G", A));
  - p*_h read from #40's raw (D4_G_h*/summary.json), not bisected again;
  - `track` wrapped to record the normal field's phi_+ on the horizon (the wrapper calls settled_mass.track and adds a column), so
    the first law and the horizon flux carry both fields;
  - the checks of #41 re-posed for G, #41's planted cases PL1 to PL3 through the same code, a comparison with family S's stored
    series (negative-mass-raw I41-rpi5), and an ulp-noise run that measures how far M_AH(t) is from the roundoff floor of #226.

    python3 settled_mass_g.py pstar   --h=H --out=DIR            read p*_h of G from #40's raw, write DIR/pstar_G_h<H>.json
    python3 settled_mass_g.py job     --h=H --delta=D --out=DIR --pstar=DIR [--noise=U --seed=S]   one run (noise U in units of 2^-52)
    python3 settled_mass_g.py planted [--mutate=NAME]            PL1 to PL3 (and the G level-0 identity) through this driver
    python3 settled_mass_g.py analyze --raw=DIR [--mutate=NAME]  all checks and reports
    python3 settled_mass_g.py scan    --h=H --out=DIR --pstar=DIR  M at the limit against the distance from threshold (report)

Units: a = 1, 8 pi G = 1, D = 4 (n = 2). M = r^(n-1)/2 on r_+ = 0. Status words as in AGENTS.md: every mass is numerical evidence.
"""
import contextlib
import hashlib
import io
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
import threshold as th         # noqa: E402

SQ2 = math.sqrt(2.0)
MUT = None
FAM, ND = "G", 2
HS = (0.02, 0.01, 0.005)                  # the three h of #40 for G at D = 4
DELTAS = (1e-6, 1e-8)
NOISE_HS = (0.01, 0.005)                  # h of the ulp-noise runs
NOISE_RUNS = ((0.0, 1), (1.0, 1), (1.0, 2), (1.0, 3), (64.0, 1))     # (noise in units of 2^-52, seed)
ULP = 2.0 ** -52
DELTA_SWAP = 1e-1                         # the delta of the mutation delta_swap (set from the development scan)
BAR_MACHINE = 2e-5                        # #41's review: the same code on another machine reproduces the settled value to about 2e-5
S_PLATEAU = 0.2231                        # #41's quoted settled mass of family S
I40 = os.path.expanduser("~/Source/negative-mass-raw/I40-dh610/runs")
I41RAW = os.path.expanduser("~/Source/negative-mass-raw/I41-rpi5")
I41OUT = os.path.join(HERE, "settled_mass_output.txt")
DECLARED = dict(n=2, F1=1.5, K=20, W1=2.4, rcut=0.05, tol=1e-2, dxp=12.0, dm=3.0, X=28.0)     # #41's declared chain, unchanged


def mut(name):
    return MUT == name


# --------------------------------------------------------------------------------------------------------------------
# track wrapper: sm.track unchanged, plus the normal field's phi_+ on the horizon
# --------------------------------------------------------------------------------------------------------------------
_track0 = sm.track


def _phip_cols(L, out):
    """phi_+ of field 0 (recomputed, for the control) and of field 1 at the recorded horizon points, by the interpolation track uses."""
    n = len(out["c"])
    p0, p1 = np.zeros(n), np.zeros(n)
    for k in range(n):
        row = out["row"][k]
        j = int(math.floor(row - 0.5))
        fr = row - 0.5 - j
        if mut("trk_shift"):
            j += 1
        _, up, _, _ = sm.cell_window(L, int(out["c"][k]), j, j + 1)
        if L.nf > 2:
            p0[k] = up[2][0] + fr * (up[2][1] - up[2][0])
        if L.nf > 3:
            p1[k] = up[3][0] + fr * (up[3][1] - up[3][0])
    return p0, p1


def track_g(L, n, c_start, jguess, **kw):
    out, stop, exhausted = _track0(L, n, c_start, jguess, **kw)
    p0, p1 = _phip_cols(L, out)
    out["phip1"] = p1
    out["phip0_re"] = p0
    return out, stop, exhausted


sm.track = track_g       # sm.chain looks `track` up in settled_mass's namespace


# --------------------------------------------------------------------------------------------------------------------
# #40's p*_h of G, and the run
# --------------------------------------------------------------------------------------------------------------------
def i40_g(h):
    """#40's bracket, p*_h, and the formation values of D4_G_h<h>; the h of its job.log as a string."""
    p = f"{I40}/D4_G_h{h:g}"
    s = json.load(open(p + "/summary.json"))
    st = json.load(open(p + "/stall.json"))["-1"]
    hstr = re.search(r"'h': '([^']+)'", open(p + "/job.log").readline()).group(1)
    ps = 0.5 * (s["lo"] + s["hi"]) + (3e-10 if mut("pstar_shift") else 0.0)
    return dict(ps=ps, lo=s["lo"], hi=s["hi"], h_summary=s["h"], h_log=hstr, event=st["event"], T=st["T"], M_hor=st["M_hor"], M_cell=st["M_cell"],
                r_cell=st["r_cell"], A_stall=st["A"], path=p)


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def pstar_cmd(h, out):
    """p*_h = (lo + hi)/2 of #40's job for G at D = 4, read from its raw folder (not bisected again), with the hashes of what was read."""
    os.makedirs(out, exist_ok=True)
    r = i40_g(h)
    res = dict(h=h, h_job_log=r["h_log"], lo=r["lo"], hi=r["hi"], ps=0.5 * (r["lo"] + r["hi"]), event=r["event"], T=r["T"], M_hor=r["M_hor"],
               A_stall=r["A_stall"], source=r["path"].replace(os.path.expanduser("~"), "~"),
               sha256={f: sha(os.path.join(r["path"], f)) for f in ("summary.json", "stall.json", "job.log", "points.jsonl")})
    json.dump(res, open(os.path.join(out, f"pstar_G_h{h:g}.json"), "w"), indent=1)
    print(json.dumps(res), flush=True)


def add_noise(d, u, seed):
    """Relative ulp-level noise u * xi, xi uniform in (-1, 1), on every node but the corner of both rays: r and the fields are scaled by 1 + u xi,
    f is shifted by u xi. u = 0 leaves the data bit for bit."""
    if u == 0.0 or mut("noise_inert"):
        return
    rng = np.random.default_rng(seed)
    for ray in (d.Sp, d.Sm):
        xi = rng.uniform(-1.0, 1.0, ray.shape)
        xi[:, 0] = 0.0
        ray[1] += u * xi[1]
        sel = [0] + list(range(2, ray.shape[0]))
        ray[sel] *= 1.0 + u * xi[sel]


def run_real_g(h, delta, ps, K=None, Nc=None, log=lambda *a: None, need_up=True, noise=0.0, seed=0, forced=None):
    """Family G at p*_h - delta (collapse side): #40's point() for the fate and the first trapped cell, then #41's chain from there."""
    F1 = sm.CRIT["F1"]
    K = sm.CRIT["K"] if K is None else K
    Nc = sm.nc_of(h) if Nc is None else Nc
    A = ps - delta
    pulses, sym, eps = th.family(FAM, A)
    X = th.CRIT["X"][FAM][ND]
    pt = th.point(ND, FAM, h, X, A)
    up = th.point(ND, FAM, h, X, ps + delta) if need_up else dict(fate=1)
    if pt["fate"] != -1 or up["fate"] != 1:
        raise RuntimeError(f"fates at p*_h -/+ delta: {pt['fate']}, {up['fate']} (wanted -1, +1)")
    j_ev, i_ev = pt["j"], pt["i"]
    Np = i_ev + int(round(sm.CRIT["dxp"] / h))
    nm0 = j_ev + int(round((sm.CRIT["dm"] + Nc * h / 4) / h))
    d = dn.ellis_rays(ND, max(Np, nm0), h, pulses=pulses, symmetric=sym, eps=eps)
    add_noise(d, noise * ULP, seed)
    out = sm.chain(ND, eps, h, d.Sp[:, :Np + 1], d.Sm, 0.0, 0.0, nm0, None, None, F1, F1 * 2.0 ** K, Nc, forced=forced, tol=sm.CRIT["tol"],
                   rcut=sm.CRIT["rcut"], log=log)
    out.update(point=pt, up=up, A=A, Np=Np, nm0=nm0, F1=F1, K=K, Nc=Nc, h=h, delta=delta, ps=ps, data=d)
    return out


def job_cmd(h, delta, out, pstar_dir, noise=0.0, seed=0):
    import resource
    ps = json.load(open(os.path.join(pstar_dir, f"pstar_G_h{h:g}.json")))["ps"]
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    r = run_real_g(h, delta, ps, log=lambda m: print(m, flush=True), noise=noise, seed=seed)
    s = r["ser"]
    np.savez(os.path.join(out, "series.npz"), **s)
    summ = dict(fam=FAM, h=h, delta=delta, ps=ps, A=r["A"], noise_ulp=noise, seed=seed, status=r["status"], event=r["event"], point=r["point"], up=r["up"],
                Np=r["Np"], nm0=r["nm0"], F1=r["F1"], K=r["K"], Nc=r["Nc"], levels=len(r["levels"]), cols=int(len(s["r"])), xp0=float(s["xp"][0]),
                xp1=float(s["xp"][-1]), F_end=float(s["F"][-1]), M_start=float(s["M"][0]), M_end=float(s["M"][-1]), min_edge=float(s["edge"].min()),
                secs=round(time.time() - t0, 1), maxrss_MB=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 1))
    json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print(json.dumps(summ), flush=True)


def load_run(raw, h, delta, tag=None):
    p = os.path.join(raw, "runs", tag or f"h{h:g}_d{delta:g}")
    return dict(np.load(os.path.join(p, "series.npz"))), json.load(open(os.path.join(p, "summary.json")))


def get_run(raw, h, delta):
    """The stored series, or, under a mutation, a fresh run from the stored p*_h."""
    if MUT is None:
        return load_run(raw, h, delta)
    ps = json.load(open(os.path.join(raw, "runs", f"pstar_G_h{h:g}.json")))["ps"]
    K = 0 if mut("no_zoom") else (12 if mut("short_limit") else None)
    ps = ps + 2e-8 if mut("pstar_off") else ps
    r = run_real_g(h, DELTA_SWAP if (mut("delta_swap") and delta == 1e-6) else delta, ps, K=K)
    return r["ser"], dict(status=r["status"], event=r["event"], point=r["point"], up=r["up"], Nc=r["Nc"], F1=r["F1"], K=r["K"], h=h, delta=delta, ps=ps)


# --------------------------------------------------------------------------------------------------------------------
# checks
# --------------------------------------------------------------------------------------------------------------------
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


def captured(fn, *a, **k):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        val = fn(*a, **k)
    return buf.getvalue(), val


def lines41(prefixes):
    """The lines of #41's output of record that start with one of the prefixes (the first block, up to the first blank line after PL3)."""
    txt = open(I41OUT).read().split("== p*_h on this machine")[0]
    return [l for l in txt.splitlines() if any(l.startswith(p) for p in prefixes)]


def planted_through_driver():
    """PL1 and PL2 (the exact Schwarzschild and Roberts solutions through the chain) and PL3 (level 0 against dn.evolve), run by #41's functions
    through this driver, with this driver's `track` in the chain. Unchanged means: the printed lines of each case equal #41's output of record
    character for character. The same lines are also produced with settled_mass's own `track`, and must equal."""
    sels = {"PL1": ("    P1 h", "[pass] PL1", "[FAIL] PL1"), "PL2": ("    P2 h", "[pass] PL2", "[FAIL] PL2")}
    try:
        ref = lines41(("    P1 h", "    P2 h", "[pass] PL1", "[pass] PL2"))
        sm.track = _track0
        txt0, _ = captured(sm.planted_checks)
        sm.track = track_g
        txt1, _ = captured(sm.planted_checks)
        for cid, pre in sels.items():
            want = [l for l in ref if l.startswith(pre)]
            got0 = [l for l in txt0.splitlines() if l.startswith(pre)]
            got1 = [l for l in txt1.splitlines() if l.startswith(pre)]
            same41, same_trk = got1 == want, got0 == got1
            check(cid, f"#41's {cid} through this driver: the printed lines equal #41's output of record, and equal those of settled_mass's own track",
                  same41 and same_trk and len(want) == 5 and got1[-1].startswith("[pass]"),
                  f"{len(got1)} lines compared with {len(want)} of #41; equal to #41: {same41}; equal without the wrapper: {same_trk}; verdict line {got1[-1][:6] if got1 else 'missing'}")
            if not same41:
                for a, b in zip(got1, want):
                    if a != b:
                        print("    differs:\n      here: " + a[:260] + "\n      #41:  " + b[:260])
    except Exception as e:
        sm.track = track_g
        for cid in sels:
            if cid not in {c for c, _ in results}:
                check(cid, f"#41's {cid} through this driver", False, "failed: " + repr(e))
    try:
        ref3 = lines41(("[pass] PL3",))
        sm.track = _track0
        txt3, _ = captured(sm.kernel_identity, I41RAW)
        sm.track = track_g
        got3 = [l for l in txt3.splitlines() if "PL3" in l]
        check("PL3", "#41's PL3 through this driver (S, h = 0.02, delta = 1e-6): level 0 equals dn.evolve bit for bit; the printed line equals #41's",
              got3 == ref3 and len(ref3) == 1, f"{got3[0][:230] if got3 else 'no line'}")
    except Exception as e:
        sm.track = track_g
        check("PL3", "#41's PL3 through this driver", False, "failed: " + repr(e))


def level0_identity_g(pdir, h=0.02, delta=1e-6):
    """PL3 for G: level 0 of the chain (no zoom: forced = []) against dn.evolve on the same data; and the control of the track wrapper on the
    level-0 horizon series: its recomputed phi_+ of field 0 equals settled_mass.track's, and field 1 is not identically zero."""
    ps = json.load(open(os.path.join(pdir, f"pstar_G_h{h:g}.json")))["ps"]
    r = run_real_g(h, delta, ps, forced=[])
    d, L0 = r["data"], r["levels"][0]
    res = dn.evolve(d, store=True, stop=False, rcut=sm.CRIT["rcut"], fcut=1e300, tol=sm.CRIT["tol"])
    same_G = np.array_equal(res.G[:, :L0.Np + 1, :L0.Nm + 1], L0.G)
    same_d = np.array_equal(res.dead[:L0.Np + 1, :L0.Nm + 1], L0.dead)
    ev = res.events.get("collapse")
    check("PL3G", "family G: level 0 of the chain equals dn.evolve bit for bit (nodes, dead flags) and finds the same first trapped cell",
          same_G and same_d and ev == r["event"], f"h = {h:g}, delta = {delta:g}: nodes equal {same_G}, dead equal {same_d}, events {ev} and {r['event']}")
    s = r["ser"]
    if len(s["r"]):
        trk_check("TRK", s)
    return r


def trk_check(cid, series_list):
    """The wrapper's recomputed phi_+ (field 0) equals settled_mass.track's on every column of the given series, and the normal field is not
    identically zero on the horizon."""
    ss = series_list if isinstance(series_list, list) else [series_list]
    n = sum(len(s["r"]) for s in ss)
    worst = max((float(np.max(np.abs(s["phip0_re"] - s["phip"]))) / max(float(np.max(np.abs(s["phip"]))), 1e-300) for s in ss if len(s["r"])), default=float("inf"))
    nz = min((float(np.max(np.abs(s["phip1"]))) for s in ss if len(s["r"])), default=0.0)
    check(cid, "track wrapper: recomputed phi_+ of field 0 equals settled_mass.track's (relative 1e-12), and the normal field's phi_+ is nonzero on the horizon",
          n > 0 and worst <= 1e-12 and nz > 0, f"{n} columns, max relative difference {worst:.1e}, smallest max|phi_+(normal)| over series {nz:.2e}")


# --------------------------------------------------------------------------------------------------------------------
# analysis helpers
# --------------------------------------------------------------------------------------------------------------------
BAND_R = sm.BAND_R


def flux_g(s):
    """dM/dx+ along the horizon for ghost + normal field: (r^2/2) e^f r_- (phi0_+^2 - phi1_+^2) (the ghost enters with the sign of settled_mass's
    flux_integral; the normal field has the opposite energy sign). Mutations: flux_sign, normal_drop, normal_sign."""
    sg = -1.0 if mut("flux_sign") else 1.0
    cn = 0.0 if mut("normal_drop") else (1.0 if mut("normal_sign") else -1.0)
    return sg * 0.5 * s["r"] ** 2 / s["F"] * s["rm"] * (s["phip"] ** 2 + cn * s["phip1"] ** 2)


def flux_integral_g(s, i0, i1):
    f = flux_g(s)
    x = s["xp"][i0:i1 + 1]
    return float(np.sum(0.5 * (f[i0:i1] + f[i0 + 1:i1 + 1]) * np.diff(x)))


def summary_g(s, h):
    """sm.settled_summary with the horizon flux over both fields."""
    q = sm.settled_summary(s, h)
    q["flux05_ghost"] = q["flux05"]
    v = s["xp"]
    q["flux05"] = float(np.abs(flux_g(s))[v > v[-1] - 0.5].mean()) if len(v) else float("nan")
    return q


def window_scatter(S1, S0, windows):
    """max |M_1 - M_0| over windows of x+ - x+_f (cubic splines on common x+), for two series of the same h."""
    from scipy.interpolate import CubicSpline
    c1, c0 = CubicSpline(S1["xp"], S1["M"]), CubicSpline(S0["xp"], S0["M"])
    vf = S0["xp"][0]
    lo_all = max(S1["xp"][0], S0["xp"][0])
    hi_all = min(S1["xp"][-1], S0["xp"][-1])
    out = []
    for a, b in windows:
        v = np.linspace(max(vf + a, lo_all), min(vf + b, hi_all), 2000)
        out.append(float(np.max(np.abs(c1(v) - c0(v)))))
    return out


WINDOWS = [(0.0, 1.5), (1.5, 3.0), (3.0, 99.0)]


def cmp_bar(Bh, Bd, Bm=BAR_MACHINE):
    return 0.0 if mut("cmp_bar0") else Bh + Bd + Bm


def cmp_verdict(dM, bar):
    return abs(dM) <= bar


def load_s41(h, delta):
    p = os.path.join(I41RAW, "runs", f"h{h:g}_d{delta:g}")
    return dict(np.load(os.path.join(p, "series.npz")))


# --------------------------------------------------------------------------------------------------------------------
# the analysis
# --------------------------------------------------------------------------------------------------------------------
def analyze_cmd(raw):
    t0 = time.time()
    hs, ds = HS, DELTAS
    Fmax = sm.CRIT["F1"] * 2.0 ** sm.CRIT["K"]
    print("\n== Declared chain (#41's constants) and the planted cases through this driver")
    bad = {k: (sm.CRIT[k], v) for k, v in DECLARED.items() if k != "X" and sm.CRIT[k] != v}
    xg = th.CRIT["X"][FAM][ND]
    check("CFG", "the chain constants of settled_mass.py and the rectangle of #40's G equal the declared values (#41's)", not bad and xg == DECLARED["X"],
          f"settled_mass.CRIT differs from declared: {bad}; X(G, D = 4) = {xg}" if bad or xg != DECLARED["X"] else f"F1 {sm.CRIT['F1']}, K {sm.CRIT['K']}, W1 {sm.CRIT['W1']}, dxp {sm.CRIT['dxp']}, dm {sm.CRIT['dm']}, X {xg}")
    planted_through_driver()
    try:
        level0_identity_g(os.path.join(raw, "runs"))
    except Exception as e:
        check("PL3G", "family G level 0 against dn.evolve", False, "failed: " + repr(e))
    print("\n== p*_h of G: #40's bracket, and the fates on this machine at p*_h -/+ 1e-10")
    ok1 = True
    for h in hs:
        r = i40_g(h)
        pts = []
        for sgn in (-1, 1):
            q = th.point(ND, FAM, h, xg, r["ps"] + sgn * 1e-10)
            pts.append(q["fate"])
        good = (r["hi"] - r["lo"] <= 1e-14 and float(r["h_log"]) == h and r["h_summary"] == h and pts == [-1, 1])
        ok1 &= good
        print(f"    h = {h:g} (job.log h = '{r['h_log']}'): p*_h = {r['ps']:.16e}, hi - lo = {r['hi'] - r['lo']:.2e}; fates at p*_h -/+ 1e-10 on this machine: {pts}")
    check("RC1", "#40's p*_h of G (hi - lo <= 1e-14, h of its job.log exact), and on this machine fate(p*_h - 1e-10) = collapse, fate(p*_h + 1e-10) = expansion", ok1)
    S, I = {}, {}
    fail = {}
    for dl in ds:
        for h in hs:
            try:
                S[(h, dl)], I[(h, dl)] = get_run(raw, h, dl)
                if len(S[(h, dl)]["r"]) == 0:
                    raise RuntimeError("empty series: " + str(I[(h, dl)]["status"]))
            except Exception as e:
                fail[(h, dl)] = repr(e)
                print(f"    run h = {h:g}, delta = {dl:g} failed: {e!r}")
    if fail:
        check("RC2", "fates, first trapped cell and its time", False, "a run failed: " + "; ".join(f"{k}: {v}" for k, v in fail.items()))
        for cid, txt in (("TRK", "track wrapper"), ("RC3", "formation value"), ("RC4", "first law"), ("RC5", "limit"), ("RC6", "convergence"), ("RC8", "settled"),
                         ("RC7", "delta independence")):
            check(cid, txt, False, "a run failed")
        return report_tail(raw, t0, S, I, failed=True)
    ref = {h: i40_g(h) for h in hs}
    ev_ok = all(list(I[(h, 1e-8)]["event"]) == list(ref[h]["event"]) for h in hs)
    t_ok = all(abs(sm.t_of((I[(h, 1e-8)]["event"][1] + 0.5) * h, (I[(h, 1e-8)]["event"][0] + 0.5) * h) - ref[h]["T"]) <= SQ2 * h for h in hs)
    fate_ok = all(I[(h, dl)]["point"]["fate"] == -1 and I[(h, dl)]["up"]["fate"] == 1 for h in hs for dl in ds)
    check("RC2", "fates (collapse at p*_h - delta, expansion at p*_h + delta; all h, both delta); the first trapped cell and its time against #40's at delta = 1e-8",
          fate_ok and ev_ok and t_ok, "events " + str({h: list(I[(h, 1e-8)]["event"]) for h in hs}) + " against #40's " + str({h: ref[h]["event"] for h in hs}))
    trk_check("TRK", [S[k] for k in S])
    dm = [abs(S[(h, 1e-8)]["M"][0] - ref[h]["M_hor"]) for h in hs]
    check("RC3", "M_AH at the first trapped cell against #40's formation value at p*_h - 1e-8", max(dm) <= 2e-3,
          "|dM| " + ", ".join(f"{x:.1e}" for x in dm) + f" (bar 2e-3); M_AH {[round(float(S[(h, 1e-8)]['M'][0]), 5) for h in hs]}, #40 {[round(ref[h]['M_hor'], 5) for h in hs]}")
    worst, parts = 0.0, []
    for key, s in S.items():
        i1 = int(np.argmax(s["F"] >= 20.0))
        dMs = s["M"][i1] - s["M"][0]
        rel = abs(dMs - flux_integral_g(s, 0, i1)) / abs(dMs)
        worst = max(worst, rel)
        gh = sm.flux_integral(s, 0, i1)
        parts.append(f"h={key[0]:g},d={key[1]:g}: dM {dMs:+.4e}, flux {flux_integral_g(s, 0, i1):+.4e} (ghost alone {gh:+.4e})")
    check("RC4", "first law on the horizon, both fields: M(x+_b) - M(x+_a) = integral of (r^2/2) e^f r_- (phi0_+^2 - phi1_+^2) along it, from the first trapped cell to F = 20",
          worst <= 0.02, f"max relative difference {worst:.2e} over all runs (bar 2e-2)")
    for p_ in parts:
        print("    " + p_)
    ok5, det = True, []
    for key, s in S.items():
        h = key[0]
        fine = np.allclose(np.diff(s["xp"]), h, atol=1e-9)
        ok5 &= I[key]["status"] == "limit" and s["F"][-1] >= 0.5 * Fmax and s["F"][-1] < Fmax and s["edge"].min() >= I[key]["Nc"] / 8 and fine
        det.append(f"h={key[0]:g},d={key[1]:g}: {I[key]['status']}, F_end {s['F'][-1]:.0f}, edge {s['edge'].min():.0f}/{I[key]['Nc']}")
    check("RC5", f"every run ends at the declared limit (F = e^-f on the horizon, {0.5 * Fmax:g} <= F_end < {Fmax:g}), horizon at least Nc/8 cells from the box edges, no gaps",
          ok5, "; ".join(det))
    print("\n== Second-order convergence of M_AH(x+): max |M_h - M_h/2| over windows of x+ - x+_f (cubic splines)")
    ok6 = True
    E = {}
    for dl in ds:
        rows = sm.band_orders({h: S[(h, dl)] for h in hs}, hs, WINDOWS, dl)
        for a_, b_, Ee, o in rows:
            good = all(BAND_R[0] <= x <= BAND_R[1] for x in o)
            ok6 &= good
            E[(dl, a_)] = Ee
            print(f"    RC6 h = {hs}, delta = {dl:g}, x+ - x+_f in [{a_:g}, {b_:g}]: E12 {Ee[0]:.2e}, E23 {Ee[1]:.2e}, order {o[0]:.2f} {'ok' if good else 'OUT'}")
    check("RC6", f"M_AH(x+) converges at second order at the three h {hs}, three windows, both deltas (observed order in {BAND_R})", ok6)
    dr = {k: summary_g(s, k[0]) for k, s in S.items()}
    check("RC8", "M_AH has settled at the limit: the plateau mean over the last 1.0 of x+ differs from the previous 1.0 by <= 1e-4 and the mean horizon flux (both fields) over the last 0.5 is <= 1e-5, every run",
          all(abs(q["dplat"]) <= 1e-4 and q["flux05"] <= 1e-5 for q in dr.values()),
          "plateau change " + ", ".join(f"{q['dplat']:+.1e}" for q in dr.values()) + "; flux " + ", ".join(f"{q['flux05']:.1e}" for q in dr.values()))
    dd = [abs(dr[(h, 1e-6)]["plat1"] - dr[(h, 1e-8)]["plat1"]) for h in hs]
    check("RC7", "the settled value does not depend on delta: plateau means over the last 1.0 at 1e-6 and 1e-8 agree", max(dd) <= 1e-4, "|dM| " + ", ".join(f"{x:.1e}" for x in dd) + " (bar 1e-4)")
    return report_tail(raw, t0, S, I)


def report_tail(raw, t0, S, I, failed=False):
    hs, ds = HS, DELTAS
    if failed:
        return finish(t0)
    S41 = {(h, dl): load_s41(h, dl) for h in hs for dl in ds}
    print("\n== M_AH(x+) of G against S (negative-mass-raw I41-rpi5, same h and delta) against x+ - x+_f; a = 1, 8 pi G = 1; M = r_AH / 2; t = (x+ + x-)/sqrt2")
    for dl in ds:
        hf = hs[-1]
        s, s4 = S[(hf, dl)], S41[(hf, dl)]
        tf, tf4 = sm.t_of(s["xp"][0], s["xm"][0]), sm.t_of(s4["xp"][0], s4["xm"][0])
        print(f"  delta = {dl:g}, h = {hf:g}: G: x+_f = {s['xp'][0]:.4f}, x-_f = {s['xm'][0]:.4f}, t_f = {tf:.4f}; S: x+_f = {s4['xp'][0]:.4f}, x-_f = {s4['xm'][0]:.4f}, t_f = {tf4:.4f}")
        print("    dx+    t-t_f(G)   r_AH(G)   M_G       M_S       M_G-M_S     e^-f(G)   e^-f(S)")
        for dv in (0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.4):
            i = int(np.argmin(np.abs(s["xp"] - (s["xp"][0] + dv))))
            i4 = int(np.argmin(np.abs(s4["xp"] - (s4["xp"][0] + dv))))
            if abs(s["xp"][i] - s["xp"][0] - dv) > 2 * hf or abs(s4["xp"][i4] - s4["xp"][0] - dv) > 2 * hf:
                break
            print(f"    {dv:5.2f} {sm.t_of(s['xp'][i], s['xm'][i]) - tf:9.4f}  {s['r'][i]:8.5f}  {s['M'][i]:8.5f}  {s4['M'][i4]:8.5f}  {s['M'][i] - s4['M'][i4]:+10.2e}  {s['F'][i]:9.4g}  {s4['F'][i4]:9.4g}")
    print("\n== Late-time values (reports)")
    print("  delta   h        M_end     plat(1.0)  range(1.0) plat(0.5)  drift(1.0)  plat change  flux(0.5)   flux(0.5) ghost alone")
    SS = {}
    for dl in ds:
        for h in hs:
            q = summary_g(S[(h, dl)], h)
            SS[(h, dl)] = q
            print(f"  {dl:6g} {h:7g}  {q['M_end']:.6f}  {q['plat1']:.6f}  {q['rng1']:.1e}  {q['plat0.5']:.6f}  {q['drift1']:+.2e}  {q['dplat']:+.2e}  {q['flux05']:.2e}  {q['flux05_ghost']:.2e}")
    print("  family S, same quantities from I41-rpi5 (the same function):")
    SS4 = {}
    for dl in ds:
        for h in hs:
            q = sm.settled_summary(S41[(h, dl)], h)
            SS4[(h, dl)] = q
            print(f"  {dl:6g} {h:7g}  {q['M_end']:.6f}  {q['plat1']:.6f}  {q['rng1']:.1e}  {q['plat0.5']:.6f}  {q['drift1']:+.2e}  {q['dplat']:+.2e}  {q['flux05']:.2e}")
    for dl in ds:
        a = [SS[(h, dl)]["M_end"] for h in hs]
        p = [SS[(h, dl)]["plat1"] for h in hs]
        p4 = [SS4[(h, dl)]["plat1"] for h in hs]
        print(f"  delta = {dl:g}, h = {hs}: G M_end {[round(x, 6) for x in a]}, differences {a[0] - a[1]:+.2e}, {a[1] - a[2]:+.2e}, Richardson (second order) M_end(0) = {a[2] - (a[1] - a[2]) / 3:.6f}; "
              f"plateau(1.0) {[round(x, 6) for x in p]}, differences {p[0] - p[1]:+.2e}, {p[1] - p[2]:+.2e}, Richardson {p[2] - (p[1] - p[2]) / 3:.6f}")
        print(f"      S plateau(1.0) {[round(x, 6) for x in p4]}, differences {p4[0] - p4[1]:+.2e}, {p4[1] - p4[2]:+.2e}, Richardson {p4[2] - (p4[1] - p4[2]) / 3:.6f}; G - S at the same h: "
              + ", ".join(f"{p[i] - p4[i]:+.2e}" for i in range(3)))
    print("\n== The comparison of the settled mass of G with S's (error bar named)")
    plG, plS = SS[(0.005, 1e-8)]["plat1"], SS4[(0.005, 1e-8)]["plat1"]
    Bh = abs(SS[(0.01, 1e-8)]["plat1"] - SS[(0.005, 1e-8)]["plat1"])
    Bd = abs(SS[(0.005, 1e-6)]["plat1"] - SS[(0.005, 1e-8)]["plat1"])
    spread = max(SS[(h, 1e-8)]["plat1"] for h in hs) - min(SS[(h, 1e-8)]["plat1"] for h in hs)
    rich = SS[(0.005, 1e-8)]["plat1"] - (SS[(0.01, 1e-8)]["plat1"] - SS[(0.005, 1e-8)]["plat1"]) / 3
    rich4 = SS4[(0.005, 1e-8)]["plat1"] - (SS4[(0.01, 1e-8)]["plat1"] - SS4[(0.005, 1e-8)]["plat1"]) / 3
    bar = cmp_bar(Bh, Bd)
    dM = plG - plS
    print(f"  G: plateau mean over the last 1.0 of x+ at h = 0.005, delta = 1e-8: {plG:.6f} (r_AH = {2 * plG:.5f}); S (I41-rpi5, same h and delta): {plS:.6f}; S quoted by #41: {S_PLATEAU}")
    print(f"  error bar of G's value (kept separate): h: finest pair |plat(0.01) - plat(0.005)| = {Bh:.2e} (all three h: range {spread:.2e}); delta: |plat(1e-6) - plat(1e-8)| at h = 0.005 = {Bd:.2e}; "
          f"machine to machine (#41's review) = {BAR_MACHINE:.0e}; sum = {bar:.2e}")
    print(f"  G - S at h = 0.005, delta = 1e-8: {dM:+.2e}; Richardson (three h, second order) G {rich:.6f}, S {rich4:.6f}, difference {rich - rich4:+.2e}; G - 0.2231 = {plG - S_PLATEAU:+.2e}")
    print(f"  VERDICT: G's settled mass {'equals' if cmp_verdict(dM, bar) else 'DIFFERS FROM'} S's: |G - S| = {abs(dM):.2e} {'<=' if cmp_verdict(dM, bar) else '>'} bar {bar:.2e}")
    pb = 3e-5
    eq_ok = cmp_verdict(0.5 * cmp_bar(pb, pb), cmp_bar(pb, pb)) and not cmp_verdict(3.0 * cmp_bar(pb, pb), cmp_bar(pb, pb))
    check("CMPT", "the comparison on planted offsets: an offset of half the bar is 'equal', three times the bar is 'differs' (the comparison can fail)", eq_ok,
          f"planted bar {cmp_bar(pb, pb):.2e}")
    return noise_report(raw, t0, S, S41)


def noise_report(raw, t0, S, S41):
    """The ulp-noise runs: how far M_AH(x+) moves under relative noise of 2^-52 in the ray data, against the difference between the h."""
    hs = HS
    have = all(os.path.exists(os.path.join(raw, "runs", f"rn_h{h:g}_u{u:g}_s{sd}", "series.npz")) for h in NOISE_HS for u, sd in NOISE_RUNS)
    print("\n== Roundoff floor: M_AH(x+) under ulp noise in the ray data (delta = 1e-8) against the h-differences")
    if not have:
        check("RN0", "noise runs present", False, "missing runs rn_h*_u*_s*")
        check("RNR", "noise runs present", False, "missing runs rn_h*_u*_s*")
        return finish(t0)
    ref = {h: S[(h, 1e-8)] for h in NOISE_HS}
    N = {}
    for h in NOISE_HS:
        for u, sd in NOISE_RUNS:
            N[(h, u, sd)], _ = load_run(raw, h, 1e-8, f"rn_h{h:g}_u{u:g}_s{sd}")
    rfs = {h: (S[(h, 1e-6)] if mut("rn_ref") else ref[h]) for h in NOISE_HS}
    z = [np.array_equal(N[(h, 0.0, 1)]["M"], rfs[h]["M"]) and np.array_equal(N[(h, 0.0, 1)]["xp"], rfs[h]["xp"]) for h in NOISE_HS]
    check("RN0", "noise 0 reproduces the run of record bit for bit (the injection is inert at zero, the job is deterministic)", all(z), f"h = {NOISE_HS}: {z}")
    sc = {}
    for h in NOISE_HS:
        for u in (1.0, 64.0):
            vals = [window_scatter(N[(h, u, sd)], rfs[h], WINDOWS) for uu, sd in NOISE_RUNS if uu == u]
            sc[(h, u)] = [max(v[w] for v in vals) for w in range(len(WINDOWS))]
            print(f"    h = {h:g}, noise {u:g} x 2^-52, {len(vals)} seed(s): max |M_noisy - M| in windows {WINDOWS}: " + ", ".join(f"{x:.2e}" for x in sc[(h, u)]))
    resp = [all(sc[(h, 64.0)][w] >= 3.0 * sc[(h, 1.0)][w] for w in range(len(WINDOWS))) for h in NOISE_HS]
    check("RNR", "the scatter responds to the noise: at 64 x 2^-52 it is at least 3 times that at 2^-52, in every window and h", all(resp),
          "; ".join(f"h = {h:g}: " + ", ".join(f"{sc[(h, 64.0)][w] / max(sc[(h, 1.0)][w], 1e-300):.1f}" for w in range(len(WINDOWS))) for h in NOISE_HS))
    print("  scatter at 2^-52 against the h-difference E between the h and its coarser neighbour (same windows, delta = 1e-8):")
    rows = sm.band_orders({h: S[(h, 1e-8)] for h in HS}, HS, WINDOWS, 1e-8)
    for k, (hf, nm) in enumerate(((0.01, "E12 (0.02 vs 0.01)"), (0.005, "E23 (0.01 vs 0.005)"))):
        for w, (a_, b_, Ee, o) in enumerate(rows):
            ratio = sc[(hf, 1.0)][w] / Ee[k]
            word = "reached" if ratio >= 1.0 else ("approached" if ratio >= 0.1 else "not reached")
            print(f"    h = {hf:g}, x+ - x+_f in [{a_:g}, {b_:g}]: noise scatter {sc[(hf, 1.0)][w]:.2e}, {nm} {Ee[k]:.2e}, ratio {ratio:.3f}: floor {word}")
    return finish(t0)


# --------------------------------------------------------------------------------------------------------------------
# supplementary: M at the limit against the distance from threshold
# --------------------------------------------------------------------------------------------------------------------
def scan_cmd(h, out, pdir):
    ps = json.load(open(os.path.join(pdir, f"pstar_G_h{h:g}.json")))["ps"]
    os.makedirs(out, exist_ok=True)
    rows = []
    print("  delta        A = p*_h - delta     x+_f     M_start   M_end     r_AH(end)   F_end   status")
    for dl in (1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-2, 1e-1):
        try:
            r = run_real_g(h, dl, ps, need_up=False)
            s = r["ser"]
            row = dict(delta=dl, A=r["A"], xpf=float(s["xp"][0]), M_start=float(s["M"][0]), M_end=float(s["M"][-1]), plat1=summary_g(s, h)["plat1"], r_end=float(s["r"][-1]),
                       F_end=float(s["F"][-1]), status=r["status"])
        except Exception as e:
            row = dict(delta=dl, status="failed: " + repr(e))
        rows.append(row)
        print("  " + "  ".join(f"{row[k]:.6g}" if isinstance(row.get(k), float) else str(row.get(k, "")) for k in ("delta", "A", "xpf", "M_start", "M_end", "plat1", "r_end", "F_end", "status")), flush=True)
    json.dump(rows, open(os.path.join(out, f"scan_G_h{h:g}.json"), "w"), indent=1)


# --------------------------------------------------------------------------------------------------------------------
# mutations: name -> the checks that must fire (computed in the development runs, declared in the raw README before the runs of record)
# --------------------------------------------------------------------------------------------------------------------
MUTATIONS = {}


def main():
    global MUT
    argv = sys.argv[1:]
    opts = dict(a.lstrip("-").split("=", 1) if "=" in a else (a.lstrip("-"), True) for a in argv[1:])
    cmd = argv[0] if argv else "help"
    if "mutate" in opts:
        MUT = opts["mutate"]
        sm.MUT = MUT
        dn.MUT = MUT if MUT in ("theta_coef", "ghost_sign", "scalar_coef") else None
        if MUT == "family_B":
            th.CRIT["G_B"] = 0.11
        if MUT == "crit_dm":
            sm.CRIT["dm"] = 2.0
    import platform
    import numba as nb
    la = open("/proc/loadavg").read().split()[:3] if platform.system() == "Linux" else ["?"]
    print(f"settled_mass_g.py {cmd} {opts}: numba {nb.__version__}, numpy {np.__version__}, python {platform.python_version()}, {platform.node()}, "
          f"load average at start {' '.join(la)}, mutation {MUT}, {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    if cmd == "pstar":
        pstar_cmd(float(opts["h"]), opts["out"])
    elif cmd == "job":
        job_cmd(float(opts["h"]), float(opts["delta"]), opts["out"], opts["pstar"], noise=float(opts.get("noise", 0.0)), seed=int(opts.get("seed", 0)))
    elif cmd == "planted":
        t0 = time.time()
        planted_through_driver()
        level0_identity_g(opts["pstar"])
        sys.exit(finish(t0))
    elif cmd == "analyze":
        sys.exit(analyze_cmd(opts["raw"]))
    elif cmd == "scan":
        scan_cmd(float(opts["h"]), opts["out"], opts["pstar"])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()

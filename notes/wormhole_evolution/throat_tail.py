#!/usr/bin/env python3
"""The power-law tail of family G's throat speed at D = 4, B = 0.02 in the nonlinear data (#291).

A thin driver. It imports settled_mass_b.py (#266; through it #41's settled_mass.py, #253's settled_mass_g.py, #40's threshold.py, #225's throat_locus.py,
#39's dn.py) and throat_drift.py (#255) unchanged.

    python3 throat_tail.py pstar  --h=0.01 --out=DIR                       bisect p*_h of G at B = 0.02 to hi - lo <= 1e-14 (#266's `pstar`, bracket (-0.01, 0.01))
    python3 throat_tail.py member --h=0.01 --e=12 --out=DIR --pstar=ROOT   the pair p*_h -/+ 10^-e through throat_locus.sweep, whole, profile slices every 0.5 in t
    python3 throat_tail.py planted [--mutate=NAME]                         planted series through the same estimator
    python3 throat_tail.py analyze --raw=DIR [--mutate=NAME]               every check and table from the stored runs

Units: a = 1 is the throat radius of the initial static wormhole, 8 pi G = 1, D = 4 (n = 2), t = (x+ + x-)/sqrt2, tau the locus's proper time.
Status words as in AGENTS.md: every number is numerical evidence.
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import settled_mass_b as b         # noqa: E402
import throat_drift as td          # noqa: E402
dn, sm, g, th, tl = b.dn, b.sm, b.g, b.th, b.tl

SQ2 = math.sqrt(2.0)
MUT = None
ND = 2
FAM, B0 = "G", 0.02
TAG = "G_B0.02"
# Brackets of the bisections (lo collapses, hi expands; p*_h at B = 0.02 is 5.7e-4 to 6.6e-4). Member 0 is the bracket #266 used for B = 0.025 and is the run of record; members 1 to 7
# are independent realizations of the bisection and the pair (the roundoff floor of p*_h is 3e-14), used for the error bar.
BRACKETS = ((-0.01, 0.01), (-0.006, 0.013), (-0.008, 0.009), (-0.012, 0.011), (-0.005, 0.0105), (-0.0095, 0.0072), (-0.0115, 0.0135), (-0.0075, 0.0085))
BRACKET = BRACKETS[0]
b.BRACKET[(FAM, B0)] = BRACKET


def mut(name):
    return MUT == name


# --------------------------------------------------------------------------------------------------------------------
# runs: p*_h (#266's `pstar` command), the pair p*_h -/+ 10^-e (#266's `member` command with the exponent e)
# --------------------------------------------------------------------------------------------------------------------
def pstar_cmd(h, out):
    b.pstar_cmd(FAM, B0, h, out)


def member_cmd(h, e, out, pdir):
    b.member_cmd(FAM, B0, h, out, pdir, e=float(e))



# --------------------------------------------------------------------------------------------------------------------
# declared criteria (the raw README, written before the runs of record)
# --------------------------------------------------------------------------------------------------------------------
H2 = (0.01, 0.005)                             # the issue's two h (gating)
H3 = (0.02, 0.01, 0.005)                       # 0.02 is added for the order in h (post hoc)
E_PAIR = 12                                    # the pair p*_h -/+ 1e-12
K_ENS = len(BRACKETS)
CRIT = dict(
    t0=15.0, t_cap=25.0, odd_max=3e-4, dt=0.05,          # window start, cap, the unstable mode's bar (#255), grid step of the speed samples
    ens_tol=0.02, min_len=2.0,                           # the window ends where the standard error of the ensemble mean of v reaches 2 % of |v|; not informative below 2.0 long
    shape_win=(10.0, 16.0),                              # the window of the shape test (rate x t of a power law is flat, an exponential's grows as t); inside the clean part of v
    p_target=2.9, p_tol=0.20,                            # issue: slope within 20 % of -2.9
    spread_frac=0.5,                                     # rate x t spread at most half of an exponential's over the same three times
    h_agree=0.10,                                        # issue: the two h agree on the slope within 10 %
    dr_t0=12.0,                                          # the speed keeps one sign and |v| falls on [12, T_end]: no drift dominates
    bar_bracket=1e-14, ord_band=(1.8, 2.2), lam_band=(1.10, 1.25), lam_win=(10.0, 16.0),
    kappa255=0.244, lin_pred=((15.0, 2.74), (20.0, 2.89), (25.0, 2.97)),
    pl_draws=200, pl_seed=291, pl_frac=0.9, pl_tol=0.25, pl_p=(3.0, 2.3), pl_kappa=0.244, pl_stress=3.0,
)
results = []


def check(cid, text, ok, detail="", ph=False):
    ok = None if ok is None else bool(ok)
    results.append((cid, ok, ph))
    print(f"[{'SKIP' if ok is None else ('pass' if ok else 'FAIL')}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def finish(t0):
    nfail = sum(1 for _, ok, ph in results if ok is False and not ph)
    nph = sum(1 for _, ok, ph in results if ok is False and ph)
    print(f"\n{len(results)} checks ({sum(1 for _, _, ph in results if ph)} post hoc, gating nothing), {nfail} gating checks failed, {nph} post hoc checks failed; {time.time() - t0:.0f} s")
    if MUT is not None:
        fired = sorted({c for c, ok, ph in results if ok is False})
        want = sorted(MUTATIONS.get(MUT, []))
        print(f"mutation {MUT}: fired {fired}, declared {want}: {'MATCH' if fired == want else 'DIFFERENT'}")
        return 0 if fired == want else 1
    return 1 if nfail else 0


# The checks each mutation must fire, computed on the development set (RAW dev/mut) before the README declared them. S1-h0.005 fails without any mutation (the dev baseline), so it is in
# every set where S1 is evaluated (not under odd_t2, ens_one, short_win and drift_add, where the slope tests are not informative and are skipped).
_B = ["S1-h0.005"]
MUTATIONS = {
    "bisect_early": ["BR-h0.005", "BR-h0.01", "BR-h0.02"] + _B,
    "evt_swap": ["EVT-h0.005", "EVT-h0.01", "EVT-h0.02"] + _B,
    "odd_t2": ["PL1", "PL2", "PL4", "PL5", "RATE-h0.005", "RATE-h0.01", "RATE-h0.02", "WIN"],
    "one_side": ["S1-h0.005", "S1-h0.01", "S2-h0.005", "S2-h0.01"],
    "ens_one": ["DR-h0.005", "PL1", "PL2", "PL4", "PL5"],
    "short_win": ["PL1", "PL2", "PL4", "PL5", "WIN"],
    "drift_add": ["DR-h0.005", "DR-h0.01"],
    "rate_not_p": ["PL1", "PL2", "PL5", "S1-h0.005", "S1-h0.01"],
    "v_exp": ["ORD-v10", "ORD-v12", "PL1", "PL2", "PL4", "PL5", "S1-h0.005", "S1-h0.01", "S2-h0.005", "S2-h0.01"],
    "s3_scale": ["S3"] + _B,
    "order_h": ["BR-h0.02", "ORD-p", "ORD-v10", "ORD-v12"] + _B,
    "pl_exp": ["PL1", "PL2", "PL3", "PL5"] + _B,
    "pl_eps": ["PL1", "PL2", "PL4", "PL5"] + _B,
}


def fmt(xs, f="{:.3f}"):
    return "[" + ", ".join(f.format(x) for x in xs) + "]"


# --------------------------------------------------------------------------------------------------------------------
# loading, the pair's mean, the speed
# --------------------------------------------------------------------------------------------------------------------
def runs_dir(raw, k):
    return os.path.join(raw, "runs" if k == 0 else f"runs_k{k}")


def load_member(raw, h, k, e=E_PAIR):
    base = os.path.join(runs_dir(raw, k), TAG, f"h{h:g}")
    ps = json.load(open(os.path.join(base, "pstar", "summary.json")))
    mdir = os.path.join(base, f"e{e}", "member")
    return dict(h=h, k=k, e=e, ps=ps, loc=dict(np.load(os.path.join(mdir, "locus.npz"))), m=json.load(open(os.path.join(mdir, "member.json"))))


def pair_arrays(loc):
    """The pair's mean on the slices both runs have (the tuned member to second order in the offset): t, tau, r_th, and odd = (r_+ - r_-)/2, the unstable mode's part."""
    sa, sb = loc["s_a"], loc["s_b"]
    common = np.intersect1d(sa, sb)
    ia, ib = np.searchsorted(sa, common), np.searchsorted(sb, common)
    t = 0.5 * (loc["t_a"][ia] + loc["t_b"][ib])
    tau = 0.5 * (loc["tau_a"][ia] + loc["tau_b"][ib])
    r = 0.5 * (loc["r_a"][ia] + loc["r_b"][ib])
    odd = 0.5 * (loc["r_b"][ib] - loc["r_a"][ia])
    if mut("one_side"):                                  # one_side: the + side alone instead of the mean (the unstable mode is kept)
        r = loc["r_b"][ib]
    return t, tau, r, odd


def series_of(t, tau, r, odd, delta=1.0 * 10.0 ** (-E_PAIR)):
    """The pair's series and its speed dr_th/dtau: #255's local linear fit over +- 0.5 in tau (throat_drift.rdot_local), samples every k/10-th slice."""
    rd = td.rdot_local(tau, r)
    ok = np.isfinite(rd)
    return dict(t=t, tau=tau, r=r, odd=odd, ts=t[ok], vs=rd[ok], delta=delta)


def series_member(mem):
    t, tau, r, odd = pair_arrays(mem["loc"])
    if mut("odd_t2"):                                    # odd_t2: the odd part's growth rate raised by 0.3
        odd = odd * np.exp(0.3 * t)
    if mut("drift_add"):                                 # drift_add: a constant drift speed of +1e-4 added
        r = r + 1e-4 * tau
    if mut("v_exp"):                                     # v_exp: an exponential speed of rate kappa (255) in place of the data's, same amplitude at t = 15
        S0 = series_of(t, tau, r, odd, 10.0 ** (-mem["e"]))
        r = planted_truth("exp", CRIT["kappa255"], tau, abs(float(np.interp(15.0, S0["ts"], S0["vs"]))), float(np.interp(15.0, t, r)))
    return series_of(t, tau, r, odd, 10.0 ** (-mem["e"]))


def v_on(S, G):
    """The speed at the times G by linear interpolation of its samples; nan outside them."""
    v = np.interp(G, S["ts"], S["vs"])
    v[(G < S["ts"][0]) | (G > S["ts"][-1])] = np.nan
    return v


def t_odd(S):
    """First t >= 6 at which |odd|/r_th reaches the bar (linear interpolation), else the last slice."""
    q = np.abs(S["odd"]) / S["r"]
    idx = np.where((S["t"] >= 6.0) & (q >= CRIT["odd_max"]))[0]
    if len(idx) == 0:
        return float(S["t"][-1])
    k = int(idx[0])
    w = (CRIT["odd_max"] - q[k - 1]) / (q[k] - q[k - 1])
    return float(S["t"][k - 1] + w * (S["t"][k] - S["t"][k - 1]))


def ens_series(members):
    """The ensemble's mean series: r_th and tau averaged over the K members on their common slices (the slice times are equal), |odd| the largest of the members' (the bar on the
    unstable mode must hold for every member). K = 1 returns the member."""
    if mut("ens_one"):                                  # ens_one: the primary pair alone (the window rule then sees no spread)
        members = members[:1]
    n = min(len(S["t"]) for S in members)
    assert all(np.max(np.abs(S["t"][:n] - members[0]["t"][:n])) < 1e-9 for S in members)
    r = np.mean([S["r"][:n] for S in members], axis=0)
    tau = np.mean([S["tau"][:n] for S in members], axis=0)
    odd = np.max([np.abs(S["odd"][:n]) for S in members], axis=0)
    return series_of(members[0]["t"][:n], tau, r, odd, members[0]["delta"])


def window(members):
    """The window [t0, T_end] of the mean series: T_end is the earliest of the cap, the odd part of the worst member reaching the bar, the last speed sample, and the first t >= t0 at which
    the standard error of the ensemble mean of v (sample std / sqrt(K)) reaches ens_tol of |mean v|: the roundoff realization of the pair, which no choice of delta removes."""
    c = CRIT
    G = np.arange(c["t0"], c["t_cap"] + 1e-9, c["dt"])
    use = members[:1] if mut("ens_one") else members
    V = np.array([v_on(S, G) for S in use])
    valid = np.all(np.isfinite(V), axis=0)
    t_last = float(G[valid][-1]) if valid.any() else float("nan")
    if len(use) >= 3 and valid.any():
        rel = np.where(valid, np.std(V, axis=0, ddof=1) / math.sqrt(len(use)) / np.abs(np.mean(V, axis=0)), np.nan)
        bad = np.where(valid & (rel > c["ens_tol"]))[0]
        t_ens = float(G[bad[0]]) if len(bad) else t_last
    else:
        rel, t_ens = np.full(len(G), np.nan), t_last
    t_o = t_odd(ens_series(members))
    return dict(T_end=min(c["t_cap"], t_o, t_ens, t_last), t_odd=t_o, t_ens=t_ens, t_last=t_last, G=G, rel=rel)


def eps_spread(members):
    """The offsets of the pairs' means from the ensemble's: (v_k - mean v) / (lambda |odd_k| / delta) at t = 12 to 15 (constant in t: the first-order offset of p*_h), sample std over k."""
    ts = np.array([12.0, 13.0, 14.0, 15.0])
    V = np.array([v_on(S, ts) for S in members])
    lam = 1.181
    Q = np.array([np.interp(ts, S["t"], np.abs(S["odd"])) / S["delta"] * lam * np.sign(np.interp(ts, S["t"], S["odd"])) for S in members])
    eps = (V - V.mean(0)) / Q
    return float(np.mean(eps.std(0, ddof=1))), eps.mean(1)


def tail_stats(S, T_end, t0=None, regress="tau"):
    """Over t in [t0, T_end]: p = -d ln|v| / d ln(tau) by the least-squares line (the log-log slope), and the quadratic fit's p at the window's start, middle and end with its spread
    against the spread an exponential's rate x t gives over the same three times ((t3 - t1) / t2)."""
    c = CRIT
    t0 = c["t0"] if t0 is None else t0
    G = np.arange(t0, T_end + 1e-9, c["dt"])
    v = v_on(S, G)
    if len(G) < 10 or not np.all(np.isfinite(v)) or np.any(v == 0):
        return dict(p=float("nan"), p3=[float("nan")] * 3, spread=float("nan"), exp_spread=float("nan"), n=len(G), same_sign=False)
    x = np.log(np.interp(G, S["t"], S["tau"])) if regress == "tau" else np.log(G)
    y = np.log(np.abs(v))
    p = -float(np.polyfit(np.exp(x) if mut("rate_not_p") else x, y, 1)[0])      # rate_not_p: the slope against tau (a rate) instead of ln tau
    q = np.polyfit(x, y, 2)
    ts3 = (t0, 0.5 * (t0 + T_end), T_end)
    xt = [math.log(np.interp(T, S["t"], S["tau"])) if regress == "tau" else math.log(T) for T in ts3]
    p3 = [-(q[1] + 2 * q[0] * xx) for xx in xt]
    spread = (max(p3) - min(p3)) / abs(np.mean(p3))
    return dict(p=p, p3=p3, spread=spread, exp_spread=(ts3[2] - ts3[0]) / ts3[1], n=len(G), same_sign=bool(np.all(v < 0) or np.all(v > 0)))


def drift_check(S, T_end):
    """No drift dominates before the window: v keeps one sign and |v| falls monotonically on [dr_t0, T_end]."""
    G = np.arange(CRIT["dr_t0"], T_end + 1e-9, CRIT["dt"])
    v = v_on(S, G)
    ok = bool(np.all(np.isfinite(v)) and (np.all(v < 0) or np.all(v > 0)) and np.all(np.diff(np.abs(v)) < 0))
    return ok, float(v[0]), float(v[-1])


# --------------------------------------------------------------------------------------------------------------------
# planted series through the same code: a power law and an exponential speed, the measured noise, the roundoff offsets of the pairs
# --------------------------------------------------------------------------------------------------------------------
def noise_rms(S):
    """rms of r_th around its local quartic fits (+- 1 in tau) over t in [12, 17]: the measured noise of the series."""
    rr = []
    for T in np.arange(12.5, 16.6, 0.5):
        k = np.abs(S["tau"] - T) <= 1.0
        pp = np.polyfit(S["tau"][k] - T, S["r"][k], 4)
        rr.append(np.std(S["r"][k] - np.polyval(pp, S["tau"][k] - T)))
    return float(np.mean(rr))


def planted_truth(kind, par, tau, A15, r15):
    """r_th(tau) whose speed is -A15 (tau/15)^-p (kind 'pow', par = p) or -A15 exp(-kappa (tau - 15)) (kind 'exp', par = kappa), with r_th(15) = r15."""
    if kind == "pow":
        p = par
        rt = lambda x: A15 * 15.0 ** p * x ** (1.0 - p) / (p - 1.0)
    else:
        kap = par
        rt = lambda x: (A15 / kap) * np.exp(-kap * (x - 15.0))
    return r15 + rt(tau) - rt(15.0)


def planted_ensemble(base, kind, par, rng, sig_eps, K=K_ENS):
    """K planted pairs sharing the base series' grid and unstable-mode profile Q = odd / delta: r_k = r_true + eps_k Q + white noise (the measured rms), eps_k ~ N(0, sig_eps^2)."""
    c = CRIT
    msk = (base["t"] >= 8.0) & (base["t"] <= 24.0)
    t, tau, odd = base["t"][msk], base["tau"][msk], base["odd"][msk]
    v15 = float(np.interp(15.0, base["ts"], base["vs"]))
    r15 = float(np.interp(15.0, base["t"], base["r"]))
    rtrue = planted_truth(kind, par, tau, abs(v15), r15)
    if mut("pl_exp") and kind == "pow":                  # pl_exp: an exponential in place of the power law
        rtrue = planted_truth("exp", c["pl_kappa"], tau, abs(v15), r15)
    sr = noise_rms(base)
    Q = odd / base["delta"]
    out = []
    for k in range(K):
        eps = rng.normal(0.0, sig_eps * (30.0 if mut("pl_eps") else 1.0))
        out.append(series_of(t, tau, rtrue + eps * Q + rng.normal(0.0, sr, len(t)), odd, base["delta"]))
    return out


def shape_stats(S, T_end):
    """The shape test's three-point rate x t on the shape window (clipped to T_end)."""
    lo, hi = CRIT["shape_win"]
    return tail_stats(S, min(hi, T_end), t0=lo)


def planted_run(base, kind, par, draws, seed, sig_eps):
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(draws):
        mem = planted_ensemble(base, kind, par, rng, sig_eps)
        w = window(mem)
        Sm = ens_series(mem)
        st = tail_stats(Sm, w["T_end"])
        sh = shape_stats(Sm, w["T_end"])
        rows.append((w["T_end"], st["p"], st["spread"], st["exp_spread"], sh["spread"], sh["exp_spread"]))
    return np.array(rows)


def check_planted(base, sig, ph=False):
    """PL1..PL5 on the base series (the stored primary pair at the finest gating h) with the offsets' measured spread sig."""
    c = CRIT
    out = {}
    for tag, kind, par, mult in (("pow3", "pow", c["pl_p"][0], 1.0), ("pow2.3", "pow", c["pl_p"][1], 1.0), ("exp", "exp", c["pl_kappa"], 1.0), ("pow3x", "pow", c["pl_p"][0], c["pl_stress"])):
        out[tag] = planted_run(base, kind, par, c["pl_draws"], c["pl_seed"] + len(tag), sig * mult)
    n = c["pl_draws"]
    for tag, truth, cid in (("pow3", c["pl_p"][0], "PL1"), ("pow2.3", c["pl_p"][1], "PL2"), ("pow3x", c["pl_p"][0], "PL5")):
        R = out[tag]
        T, p, ssh, esh = R[:, 0], R[:, 1], R[:, 4], R[:, 5]
        inf = (T - c["t0"]) >= c["min_len"]
        okp = inf & (np.abs(p - truth) <= c["pl_tol"])
        sig_ = sig * (c["pl_stress"] if tag == "pow3x" else 1.0)
        check(cid, f"planted t^-{truth:g} speed, the measured noise and offsets ({sig_:.2e} per pair, K = {K_ENS}, {n} draws): the slope over [{c['t0']:g}, T_end] is recovered to {c['pl_tol']} in at least {c['pl_frac']:.0%}",
              okp.mean() >= c["pl_frac"], f"recovered {okp.mean():.0%}; p median {np.median(p):.3f} (5 to 95 %: {np.percentile(p, 5):.3f} to {np.percentile(p, 95):.3f}); T_end median {np.median(T):.2f} (5 %: {np.percentile(T, 5):.2f}); informative {inf.mean():.0%}", ph=ph)
    print("  (report) rate x t spread on [15, T_end] itself, planted: t^-3 median %.3f (95 %%: %.3f), exponential median %.3f, exponential's own %.3f: %s" % (
        np.median(out["pow3"][:, 2]), np.percentile(out["pow3"][:, 2], 95), np.median(out["exp"][:, 2]), np.median(out["exp"][:, 3]),
        "separates them" if np.percentile(out["pow3"][:, 2], 95) < np.median(out["exp"][:, 2]) else "does not separate them: the window is too short for the noise, so the shape test uses the shape window"))
    R = out["pow3"]
    s2 = np.mean(R[:, 4] <= c["spread_frac"] * R[:, 5])
    check("PL3", "planted t^-3: the shape test (rate x t spread on the shape window) passes in at least the declared fraction (a power law is not called an exponential)", s2 >= c["pl_frac"],
          f"passes {s2:.0%}; shape spread median {np.median(R[:, 4]):.3f} against the exponential's {np.median(R[:, 5]):.3f}", ph=ph)
    R = out["exp"]
    inf = (R[:, 0] - c["t0"]) >= c["min_len"]
    f2 = np.mean(R[:, 4] > c["spread_frac"] * R[:, 5])
    f1 = np.mean(inf & (np.abs(R[:, 1] - c["p_target"]) > c["p_tol"] * c["p_target"]))
    check("PL4", f"planted exponential of rate {c['pl_kappa']} with the same noise: the shape test and the slope test both fail (a negative control that must fire) in at least {c['pl_frac']:.0%}",
          f2 >= c["pl_frac"] and f1 >= c["pl_frac"], f"shape test fails {f2:.0%}, slope test fails {f1:.0%}; slope median {np.median(R[:, 1]):.3f}; shape spread median {np.median(R[:, 4]):.3f} against {np.median(R[:, 5]):.3f}; informative {inf.mean():.0%}", ph=ph)
    return out


# --------------------------------------------------------------------------------------------------------------------
# the analysis of the stored runs
# --------------------------------------------------------------------------------------------------------------------
def local_p(S, T, w=0.5):
    """Local exponent -d ln|v| / d ln tau at t = T by the central difference over +- w in t."""
    a, c = v_on(S, np.array([T - w, T + w]))
    ta, tc = np.interp([T - w, T + w], S["t"], S["tau"])
    return float(-(math.log(abs(c)) - math.log(abs(a))) / (math.log(tc) - math.log(ta)))


def order_of(a, b_, c):
    """Observed order from three values at h, h/2, h/4 of the coarse-to-fine differences: log2((a - b)/(b - c))."""
    d1, d2 = a - b_, b_ - c
    return math.log2(d1 / d2) if d1 * d2 > 0 else float("nan")


def load_all(raw):
    """members[h][k] for every h in H3 and k in range(K_ENS) whose files exist (h = 0.02: the primary only is required)."""
    out = {}
    for h in H3:
        out[h] = {}
        for k in range(K_ENS):
            try:
                out[h][k] = load_member(raw, h, k)
            except (FileNotFoundError, OSError):
                pass
    return out


def analyze_cmd(raw):
    t_start = time.time()
    c = CRIT
    if mut("short_win"):                                 # short_win: a window of at least 20 is required
        c["min_len"] = 20.0
    M = load_all(raw)
    if mut("order_h"):                                   # order_h: the h = 0.02 records replaced by the h = 0.01 ones
        M[0.02] = {k: dict(m, h=0.02) for k, m in M[0.01].items()}
    S = {h: {k: series_member(m) for k, m in M[h].items()} for h in M}
    print("records: " + ", ".join(f"h = {h:g}: members {sorted(M[h])}" for h in H3))
    # ---- integrity, per member
    for h in H3:
        ph = h not in H2
        bad = []
        for k, m in M[h].items():
            pts = [json.loads(l) for l in open(os.path.join(runs_dir(raw, k), TAG, f"h{h:g}", "pstar", "points.jsonl"))]
            lo, hi = b.bracket_after(pts, 30 if mut("bisect_early") else len(pts))
            ends = pts[0]["fate"] == -1 and pts[1]["fate"] == 1
            if not (hi - lo <= c["bar_bracket"] and lo == m["ps"]["lo"] and hi == m["ps"]["hi"] and m["ps"]["stopped"] is None and ends):
                bad.append(k)
        gaps = [m["ps"]["hi"] - m["ps"]["lo"] for m in M[h].values()]
        check(f"BR-h{h:g}", f"every bisection ends with hi - lo <= {c['bar_bracket']:.0e}, no stop, the replayed bracket equals the stored one, lo collapses and hi expands", len(M[h]) == K_ENS and not bad, f"{len(M[h])} of {K_ENS} members, max gap {max(gaps):.1e}, bad {bad}", ph=ph)
        first, second = ("ev40_b", "ev40_a") if mut("evt_swap") else ("ev40_a", "ev40_b")
        badf = [k for k, m in M[h].items() if not (len(m["m"][first]) and next(iter(m["m"][first])) == "collapse" and len(m["m"][second]) and next(iter(m["m"][second])) != "collapse")]
        check(f"EVT-h{h:g}", "the pair p*_h -/+ 1e-12 collapses on the - side and expands on the + side, inside the triangle, in every member", len(M[h]) == K_ENS and not badf, f"bad {badf}", ph=ph)
        lam = []
        for k, sr in S[h].items():
            w = (sr["t"] >= c["lam_win"][0]) & (sr["t"] <= c["lam_win"][1])
            lam.append(float(np.polyfit(sr["t"][w], np.log(np.abs(sr["odd"][w])), 1)[0]))
        check(f"RATE-h{h:g}", f"the pair's odd part grows at the unstable mode's rate on t in {c['lam_win']}: {c['lam_band'][0]} to {c['lam_band'][1]}, in every member",
              len(lam) == K_ENS and all(c["lam_band"][0] <= x <= c["lam_band"][1] for x in lam), f"{fmt(lam)}", ph=ph)
    # ---- the mean series and the window
    mem = {h: [S[h][k] for k in sorted(S[h])] for h in H3 if S[h]}
    SM = {h: ens_series(mem[h]) for h in mem}
    W = {h: window(mem[h]) for h in H2 if h in mem}
    sig = {h: eps_spread(mem[h])[0] for h in mem if len(mem[h]) >= 3}
    for h in W:
        w = W[h]
        print(f"  h = {h:g}: window end {w['T_end']:.2f} (a member's odd part reaches {c['odd_max']:.0e} at {w['t_odd']:.2f}; the mean's standard error reaches {c['ens_tol']:.0%} at {w['t_ens']:.2f}; last sample {w['t_last']:.2f}); K = {len(mem[h])}; offset spread of the pairs' means {sig.get(h, float('nan')):.2e}")
    T_end = min(w["T_end"] for w in W.values()) if W else float("nan")
    inform = bool(np.isfinite(T_end) and T_end - c["t0"] >= c["min_len"])
    check("WIN", f"the common window [{c['t0']:g}, T_end] is at least {c['min_len']:g} long, so the slope tests are informative", inform, f"T_end = {T_end:.2f}")
    # ---- no drift before the window
    drift_ok = {}
    for h in W:
        ok, v0, v1 = drift_check(SM[h], T_end) if inform else (None, float("nan"), float("nan"))
        drift_ok[h] = ok
        check(f"DR-h{h:g}", f"the speed keeps one sign and |v| falls on [{c['dr_t0']:g}, T_end] (no drift dominates before the window)", ok, f"v = {v0:+.3e} at {c['dr_t0']:g}, {v1:+.3e} at {T_end:.2f}")
    inf_all = inform and all(drift_ok.values())
    # ---- the tests of the issue, on the mean series of each h
    st, sh = {}, {}
    for h in W:
        st[h] = tail_stats(SM[h], T_end) if inf_all else None
        sh[h] = shape_stats(SM[h], T_end) if inf_all else None
    for h in W:
        r, q = st[h], sh[h]
        if r is None:
            check(f"S1-h{h:g}", "slope within 20 % of -2.9: not informative", None, "window too short or drift present")
            check(f"S2-h{h:g}", "rate x t grows less than half as fast as an exponential's: not informative", None, "window too short or drift present")
            continue
        check(f"S1-h{h:g}", f"the log-log slope of |v| over t in [{c['t0']:g}, {T_end:.2f}] is within {c['p_tol']:.0%} of -{c['p_target']}", abs(r["p"] - c["p_target"]) <= c["p_tol"] * c["p_target"],
              f"slope -{r['p']:.3f} (|p - 2.9| = {abs(r['p'] - c['p_target']):.3f}, bar {c['p_tol'] * c['p_target']:.3f}), {r['n']} samples")
        lo, hi = c["shape_win"]
        check(f"S2-h{h:g}", f"rate x t on the shape window [{lo:g}, {min(hi, T_end):g}] at its start, middle and end: spread at most {c['spread_frac']:.0%} of an exponential's",
              q["spread"] <= c["spread_frac"] * q["exp_spread"], f"rate x t = {fmt(q['p3'])}, spread {q['spread']:.1%} of the mean against {q['exp_spread']:.1%} for a constant rate (bar {c['spread_frac'] * q['exp_spread']:.1%})")
    if inf_all and len(W) == 2:
        ps_ = [st[h]["p"] * (1.25 if (mut("s3_scale") and h == H2[1]) else 1.0) for h in H2]
        d = abs(ps_[0] - ps_[1]) / (0.5 * (ps_[0] + ps_[1]))
        check("S3", f"the two h agree on the slope within {c['h_agree']:.0%}", d <= c["h_agree"], f"-{ps_[0]:.3f} (h = {H2[0]:g}), -{ps_[1]:.3f} (h = {H2[1]:g}): {d:.1%}")
    else:
        check("S3", "the two h agree on the slope: not informative", None)
    # ---- order in h (post hoc: needs h = 0.02)
    if all(h in M and 0 in M[h] for h in H3):
        pst = [M[h][0]["ps"]["ps"] for h in H3]
        o = order_of(*pst)
        check("ORD-p", f"p*_h converges at second order over h = 0.02, 0.01, 0.005 (observed order in {c['ord_band']})", c["ord_band"][0] <= o <= c["ord_band"][1], f"p*_h = {fmt(pst, '{:.7e}')}, order {o:.2f}", ph=True)
        for T in (10.0, 12.0):
            vv = [float(v_on(SM[h], np.array([T]))[0]) for h in H3]
            o = order_of(*vv)
            check(f"ORD-v{T:g}", f"the speed at t = {T:g} converges at second order over the three h", c["ord_band"][0] <= o <= c["ord_band"][1], f"v = {fmt(vv, '{:.5e}')}, order {o:.2f}", ph=True)
    # ---- reports
    report(M, S, SM, W, T_end, st, inf_all, mem)
    # ---- planted: the base series is the stored primary pair at the finest gating h; the offsets' spread is the measured one
    if sig and 0.005 in S and 0 in S[0.005]:
        sg = math.sqrt(np.mean([sig[h] ** 2 for h in H2 if h in sig]))
        print(f"\nplanted cases: offset spread of the pairs' means, pooled over the two h: {sg:.2e}")
        check_planted(S[0.005][0], sg)
    return finish(t_start)


def report(M, S, SM, W, T_end, st, inf_all, mem):
    c = CRIT
    print("\n--- local exponent -d ln|v| / d ln tau (central difference over +-0.5 in t) of the mean series, and the standard error of its v ---")
    print("t      " + "  ".join(f"h={h:g}: p" for h in H3) + "   v (h = 0.005) | SE/|v| of the mean (K) at h = 0.01, 0.005 | primary p (h = 0.005)")
    for T in (8.0, 10.0, 12.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0):
        row = []
        for h in H3:
            try:
                row.append(f"{local_p(SM[h], T):6.3f}")
            except Exception:
                row.append("   n/a")
        se = []
        for h in (0.01, 0.005):
            ens = [v_on(S[h][k], np.array([T]))[0] for k in S[h]]
            ens = [x for x in ens if np.isfinite(x)]
            se.append(f"{np.std(ens, ddof=1) / math.sqrt(len(ens)) / abs(np.mean(ens)):.1e}" if len(ens) >= 3 else "n/a")
        try:
            pp = f"{local_p(S[0.005][0], T):6.3f}"
        except Exception:
            pp = "n/a"
        print(f"{T:5.1f}  " + "   ".join(row) + f"    {v_on(SM[0.005], np.array([T]))[0]:+.4e} | {se[0]}, {se[1]} ({len(mem[0.01])}, {len(mem[0.005])}) | {pp}")
    print("\n--- rate x t on [15, T_end] itself (quadratic fit at the window's start, middle, end), mean series (report, gates nothing) ---")
    for h in H2:
        r_ = st.get(h)
        if r_ is not None:
            print(f"  h = {h:g}: rate x t = {fmt(r_['p3'])}, spread {r_['spread']:.1%} of the mean against {r_['exp_spread']:.1%} for a constant rate; the line fit gives {r_['p']:.3f}")
    print("\n--- the member along the stall (the primary pair's mean profile, #225's fit on each stored slice): t, r_th - 1, gamma_1, rms of the fit of M(r) ---")
    for h in H2:
        rows = M[h][0]["m"]["rows"]
        tt = np.array([r_["t"] for r_ in rows])
        out = []
        for T in (6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0):
            r_ = rows[int(np.argmin(np.abs(tt - T)))]
            out.append(f"t = {r_['t']:.2f}: {r_['r_th'] - 1:+.3e}, {r_['g1_m']:.5f}, {r_['rms_m']:.1e}")
        print(f"  h = {h:g}: " + " | ".join(out))
    print("\n--- linear response of #272 (rate x t at t = 15, 20, 25): " + ", ".join(f"{p:.2f}" for _, p in c["lin_pred"]) + "; the rate of #255's fit kappa = 0.244 gives rate x t = 3.66 at t = 15 ---")
    print("\n--- sensitivity of the slope (mean series, h = 0.005 | 0.01), regressor ln tau | ln t; windows beyond T_end are marked * (contaminated) ---")
    for t0 in (12.0, 13.0, 14.0, 15.0):
        for T in (16.0, 17.0, T_end):
            row = []
            for h in (0.005, 0.01):
                a = tail_stats(SM[h], T, t0=t0, regress="tau")["p"]
                b_ = tail_stats(SM[h], T, t0=t0, regress="t")["p"]
                row.append(f"{a:.3f} | {b_:.3f}")
            print(f"  [{t0:g}, {T:.2f}]{'*' if T > T_end + 1e-9 else ' '}  " + "     ".join(row))
    print("\n--- error bars of the slope over the common window (kinds kept apart) ---")
    for h in H2:
        ps_ = [tail_stats(S[h][k], T_end)["p"] for k in sorted(S[h])]
        print(f"  h = {h:g}: roundoff realization, the K = {len(ps_)} single pairs: {fmt(ps_)}, mean {np.mean(ps_):.3f}, std {np.std(ps_, ddof=1):.3f} (the mean series: {tail_stats(SM[h], T_end)['p']:.3f}; the primary pair: {ps_[0]:.3f})")
    print(f"  truncation: |p(0.01) - p(0.005)| of the mean series = {abs(tail_stats(SM[0.01], T_end)['p'] - tail_stats(SM[0.005], T_end)['p']):.3f}")
    print("  fit choices: the spread over the windows and regressors in the table above")
    print("\n--- tail plus drift fit v = -A tau^-p + b over t in [t_start, T_end] (post hoc, mean series): the drift #255's B^3.16 scaling predicts is 6e-3 of B = 0.1's +1.26e-3, i.e. +8e-6 ---")
    try:
        from scipy.optimize import curve_fit
        f3 = lambda x_, A, p, b_: -A * x_ ** (-p) + b_
        f2 = lambda x_, A, p: -A * x_ ** (-p)
        for h in H2:
            Sx = SM[h]
            for ts in (8.0, 10.0, 12.0):
                G = np.arange(ts, T_end + 1e-9, c["dt"])
                v = v_on(Sx, G)
                x = np.interp(G, Sx["t"], Sx["tau"])
                pp, cov = curve_fit(f3, x, v, p0=[0.1, 2.5, 0.0], maxfev=20000)
                se = np.sqrt(np.diag(cov))
                p2, _ = curve_fit(f2, x, v, p0=[0.1, 2.5], maxfev=20000)
                print(f"  h = {h:g}, [{ts:g}, {T_end:.2f}]: A = {pp[0]:.3e}, p = {pp[1]:.3f} +- {se[1]:.3f}, b = {pp[2]:+.3e} +- {se[2]:.1e}, rms {np.sqrt(np.mean((f3(x, *pp) - v) ** 2)):.1e}; without b: p = {p2[1]:.3f}, rms {np.sqrt(np.mean((f2(x, *p2) - v) ** 2)):.1e}")
    except Exception as ex:
        print(f"  fit failed: {ex!r}")

def main():
    global MUT
    argv = sys.argv[1:]
    cmd = argv[0] if argv else "help"
    opts = dict(a.lstrip("-").split("=", 1) if "=" in a else (a.lstrip("-"), True) for a in argv[1:])
    if "mutate" in opts:
        MUT = opts["mutate"]
        if MUT not in MUTATIONS:
            raise SystemExit(f"unknown mutation {MUT}; known: {', '.join(sorted(MUTATIONS))}")
    import platform
    import numba as nb
    print(f"throat_tail.py {cmd} {opts}: numba {nb.__version__}, numpy {np.__version__}, python {platform.python_version()}, {platform.node()}, "
          f"mutation {MUT}, {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    if "k" in opts:
        b.BRACKET[(FAM, B0)] = BRACKETS[int(opts["k"])]
    if cmd == "pstar":
        pstar_cmd(float(opts["h"]), opts["out"])
    elif cmd == "member":
        member_cmd(float(opts["h"]), float(opts["e"]), opts["out"], opts["pstar"])
    elif cmd == "planted":
        t0 = time.time()
        raw = opts["raw"]
        M = load_all(raw)
        sg = math.sqrt(np.mean([eps_spread([series_member(m) for m in M[h].values()])[0] ** 2 for h in H2]))
        print(f"offset spread of the pairs' means, pooled over the two h: {sg:.2e}")
        check_planted(series_member(M[0.005][0]), sg)
        sys.exit(finish(t0))
    elif cmd == "analyze":
        sys.exit(analyze_cmd(opts["raw"]))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Does the corrected eps = 1e-2 reading predict #305's route, frozen first? Stored runs only. D = 4, delta = 1e-8, h = 0.005 (#317).

#296 read the settled-mass members at eps = 1e-2 and found the pair's mean profile carries a second-order part of the unstable mode (r_th - 1 = -1.41 eps^2, rdot_th = -3.24 eps^2 at the static
member B = 0). Its variant V3 subtracts B = 0's reading from every member (M_0 = 0.223086 fixed); V1 adds a free intercept. Both were declared, not predicted blind. This driver

    freeze   fits P4's form y = M / r_th - M_0 = a g + b g^2 + c rdot_th over the fifteen training members without delta = 1e-12 (14 points) under five readings and writes the coefficients
             (frozen_v3.json). It never reads X1 to X3.
    analyze  reads the frozen file (it must be committed and pushed), then reads #305's X1 to X3 from their stored runs and predicts their settled mass. Gating: V3 at eps = 1e-2, bar 6e-5.
             Reported, not gating: V1 at eps = 1e-2, V3 at eps = 1e-3, the unsubtracted refits at both eps, #305's frozen P4.

    python3 settled_mass_blind.py freeze  [--part=record|dev]
    python3 settled_mass_blind.py analyze [--part=record|dev] [--mutate=NAME] [--baseline=FILE]
    python3 settled_mass_blind.py post    post hoc, after the record, gates nothing: the shift that would close the miss, other eps, refits that include X

`--part=dev` is the development partition: fifteen training members only (the fit is 11 of the 14, the prediction is three held out), X1 to X3 and N1 to N3 are never read (a guard raises).
The mutations' declared effects (MUTATIONS) are computed there before the freeze. Readings are #296's (settled_mass_eps.read, through #266's member_at and #267's speed_at), imported unchanged.

Units: a = 1 is the throat radius of the initial static wormhole, 8 pi G = 1, D = 4. Status words as in AGENTS.md: every number is numerical evidence.
"""
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import settled_mass_eps as E               # noqa: E402  (#296's driver: the readings, the form, V3's subtraction; imports #295's, #283's, #266's, #267's unchanged)
import settled_mass_route4 as s4           # noqa: E402  (#305's route: the declared X1 to X3)
r, b, m, g, tl = E.r, E.b, E.m, E.g, E.tl

MUT = None
PART = "record"
MODE = "analyze"
H = E.H_FIN
EPS_A, EPS_B = 1e-3, 1e-2
BAR = 6e-5
M0 = 0.223086
REF = E.REF
RAW_REPO = E.RAW_ROOT
RAW = os.path.join(RAW_REPO, "I317-mac")
I305 = os.path.join(RAW_REPO, "I305-mac")
FIT = tuple(n for n in E.TRAIN9 if n != "e12") + E.SIX          # the fifteen without delta = 1e-12: 14 points
DEV_TEST = ("B0.05", "C2", "e10")                               # held out in development: one B point, one route-C point, one delta point (chosen by structure, before any run)
DEV_FIT = tuple(n for n in FIT if n not in DEV_TEST)
XN = tuple(n for n, _, _ in s4.POINTS4)
XSPEC = {n: s for n, s, _ in s4.POINTS4}
N3 = E.THREE
FITS = (("V3_1e-2", EPS_B, True, "P3"), ("V1_1e-2", EPS_B, False, "free"), ("V3_1e-3", EPS_A, True, "P3"), ("FIX_1e-2", EPS_B, False, "P3"), ("FIX_1e-3", EPS_A, False, "P3"))
GATING_KEY = "V3_1e-2"
OUT305 = os.path.join(HERE, "settled_mass_route4_output.txt")
OUT296 = os.path.join(HERE, "settled_mass_eps_output.txt")

results = []
BASE = set()


def mut(name):
    return MUT == name


def check(cid, text, ok, detail="", ph=False):
    ok = None if ok is None else bool(ok)
    results.append((cid, ok, ph))
    print(f"[{'SKIP' if ok is None else ('pass' if ok else 'FAIL')}] {cid}: {'NON-GATING (declared): ' if ph else ''}{text}" + (f"  ({detail})" if detail else ""), flush=True)


def finish(t0):
    nfail = sum(1 for _, ok, ph in results if ok is False and not ph)
    nph = sum(1 for _, ok, ph in results if ok is False and ph)
    print(f"\n{len(results)} checks ({sum(1 for _, _, ph in results if ph)} non-gating), {nfail} gating checks failed, {nph} non-gating checks failed; {time.time() - t0:.0f} s")
    if MUT is not None:
        fired = sorted({c for c, ok, ph in results if ok is False})
        new = sorted(set(fired) - BASE)
        want = sorted(set(MUTATIONS.get(MUT, [])) - BASE)
        print(f"mutation {MUT}: fired {fired}; failing without a mutation (--baseline) {sorted(BASE)}; newly fired {new}, declared {want}: {'MATCH' if new == want else 'DIFFERENT'}")
        return 0 if new == want else 1
    return 1 if nfail else 0


# --------------------------------------------------------------------------------------------------------------------
# sets, files, reading
# --------------------------------------------------------------------------------------------------------------------
def split():
    """(training names, test names) of the part."""
    train, test = (DEV_FIT, DEV_TEST) if PART == "dev" else (FIT, XN)
    if mut("keep_e12"):
        train = train + ("e12",)
    return train, test


def frozen_path():
    return os.path.join(RAW, "dev", "frozen_dev.json") if PART == "dev" else os.path.join(RAW, "frozen_v3.json")


def rd(name, eps):
    """The member at departure eps (settled_mass_eps.read: #266's member_at on the pair's mean profile, #267's speed_at). X and N are never read in development or by `freeze`."""
    if name in XN and (PART == "dev" or MODE == "freeze"):
        raise RuntimeError(f"{PART}/{MODE}: {name} (#305's member) must not be read")
    if PART == "dev" and (name in N3 or name not in set(FIT) | {REF, "e12"}):
        raise RuntimeError(f"development partition: {name} is not one of the fifteen training members")
    if MODE == "freeze" and name in N3:
        raise RuntimeError(f"freeze: {name} must not be read")
    if mut("no_run") and name == "W1" and eps == EPS_B:
        raise RuntimeError("treated as missing")
    return E.read(name, eps)


def load_stored(root, spec, h=H):
    """The stored run of one route point (settled_mass_route.load_point), the same call that reads #283's, #295's and #305's points."""
    e = r.load_point(os.path.join(root, "runs"), spec, 0.01 if mut("x_loader_h") else h)
    if "err" in e:
        raise RuntimeError(e["err"])
    return e


def seed_x():
    for n in XN:
        E._E[(n, H)] = load_stored(I305, XSPEC[n])


def readset(names, eps):
    return {n: rd(n, eps) for n in names}


def finite(p):
    return all(np.isfinite(p[k]) for k in ("g1", "rth", "rdot", "M", "t_dep"))


def arr(P, names, k):
    return np.array([P[n][k] for n in names], float)


# --------------------------------------------------------------------------------------------------------------------
# the form, the subtraction, the fits (freeze) and the frozen prediction (analyze)
# --------------------------------------------------------------------------------------------------------------------
def offsets_of(ref):
    """B = 0's reading: r_th - 1 and rdot_th (both zero for the static member)."""
    return dict(dr=float(ref["rth"] - 1.0), dv=float(ref["rdot"]))


def apply_offsets(P, dr, dv):
    """V3: subtract B = 0's offsets from every member's r_th and rdot_th (settled_mass_eps.sub_ref's operation)."""
    s = -1.0 if mut("sub_sign") else 1.0
    if mut("no_subtract"):
        dr = dv = 0.0
    if mut("sub_r_only"):
        dv = 0.0
    if mut("sub_v_only"):
        dr = 0.0
    return {n: dict(p, rth=p["rth"] - s * dr, rdot=p["rdot"] - s * dv) for n, p in P.items()}


def make_fit(P, train, eps, sub, kind, m0, Pref=None):
    """Fit over `train`. P holds the members read at eps (B = 0 included when sub). Returns the frozen record of the fit."""
    b0 = None
    Pv = P
    if sub:
        ref = Pref if Pref is not None else P[REF]
        b0 = offsets_of(ref)
        Pv = apply_offsets(P, b0["dr"], b0["dv"])
    R = E.fit_pred(Pv, train, train[:1], m0, kind)
    return dict(eps=eps, kind=kind, sub=sub, co=[float(x) for x in R["co"]], b0=b0, ins_max=float(np.max(np.abs(R["ins"]))), loo_max=float(np.max(np.abs(R["loo"]))), m0=m0)


def predict(fit, P, names):
    """The frozen prediction M = r_th' (M_0 + design . co) with r_th', rdot_th' the readings less the frozen offsets; the coefficient mutations act here."""
    Pv = apply_offsets(P, fit["b0"]["dr"], fit["b0"]["dv"]) if fit["sub"] else P
    co = np.array(fit["co"], float)
    kind = "nospeed" if len(co) == 2 else fit["kind"]
    if mut("coef_a"):
        co[0] *= 1.05
    if mut("coef_c") and len(co) > 2:
        co[2] *= 0.8
    A = E.design(arr(Pv, names, "g1"), arr(Pv, names, "rdot"), kind)
    rs = np.ones(len(names)) if mut("no_rscale") else arr(Pv, names, "rth")
    return rs * (fit["m0"] + A @ co), Pv


def resid(fit, P, names):
    pred, Pv = predict(fit, P, names)
    M = arr(P, names, "M")
    if mut("test_swap") and len(names) > 1:
        M = np.roll(M, 1)
    return M - pred, pred, Pv


def make_all(P1, P2, train, m0):
    out = {}
    for key, eps, sub, kind in FITS:
        P = P2 if eps == EPS_B else P1
        Pref = P1[REF] if (mut("b0_eps") and eps == EPS_B) else None
        out[key] = make_fit(P, train, eps, sub, kind, m0, Pref)
    return out


def inputs_hash(P1, P2, train):
    h = hashlib.sha256()
    for eps, P in ((EPS_A, P1), (EPS_B, P2)):
        for n in (REF,) + tuple(train):
            h.update(repr((n, eps, P[n]["g1"], P[n]["rth"], P[n]["rdot"], P[n]["M"])).encode())
    return h.hexdigest()


def freeze_cmd():
    t0 = time.time()
    train, _ = split()
    names = (REF,) + train
    P1, P2 = readset(names, EPS_A), readset(names, EPS_B)
    fits = make_all(P1, P2, train, E.m0_const(H))
    out = dict(part=PART, issue=317, m0=M0, bar=BAR, h=H, train=list(train), ref=REF, inputs_sha256=inputs_hash(P1, P2, train), fits=fits,
               script=os.path.basename(__file__), script_sha256=hashlib.sha256(open(__file__, "rb").read()).hexdigest(), written=time.strftime("%Y-%m-%d %H:%M:%S"))
    path = frozen_path()
    if os.path.exists(path):
        raise SystemExit(f"{path} exists: a frozen file is never overwritten")
    json.dump(out, open(path, "w"), indent=1)
    print(f"froze {len(train)} points {train} at h = {H:g}, M_0 = {M0}; B = 0's readings: eps 1e-2 r_th - 1 {P2[REF]['rth'] - 1:+.4e}, rdot_th {P2[REF]['rdot']:+.4e}; eps 1e-3 {P1[REF]['rth'] - 1:+.4e}, {P1[REF]['rdot']:+.4e}")
    for key, f in fits.items():
        print(f"    {key:9s} (a, b, c{', intercept' if len(f['co']) > 3 else ''}) = {', '.join(f'{x:+.12f}' for x in f['co'])}; in-sample max {f['ins_max']:.3e}, leave-one-out max {f['loo_max']:.3e}"
              + (f"; subtracts dr {f['b0']['dr']:+.6e}, dv {f['b0']['dv']:+.6e}" if f["b0"] else ""))
    print(f"inputs sha256 {out['inputs_sha256']}\nwrote {path}; {time.time() - t0:.0f} s")
    return 0


# --------------------------------------------------------------------------------------------------------------------
# printed numbers of the earlier notes
# --------------------------------------------------------------------------------------------------------------------
def parse296_table():
    """#296's output: the members at departure at both eps: t_dep, gamma_1, r_th, rdot_th at 1e-3 and 1e-2, and M. Returns {name: dict(field -> (value, half unit))}."""
    out = {}
    pat = re.compile(r"^\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+ \| (\S+) (\S+) \S+ \| (\S+) (\S+) \S+ \| (\S+) (\S+) \S+ \| (\S+)\s*$")
    for ln in open(OUT296):
        mt = pat.match(ln)
        if mt and re.match(r"^(B0(\.\d+)?|e\d+|[WCN][123])$", mt.group(1)):
            d = {}
            for key, gi in (("t1", 2), ("t2", 3), ("g1_1", 4), ("g1_2", 5), ("r1", 6), ("r2", 7), ("v1", 8), ("v2", 9), ("M", 10)):
                d[key] = E.pnum(mt.group(gi))
            out[mt.group(1)] = d
    return out


def parse296_v3():
    """#296's post hoc 2: max |M - predicted| of T1 and T2 without e12 at eps = 1e-3 and 1e-2 for the fixed M_0, V1 and V3 readings."""
    out = {}
    pat = re.compile(r"eps (\S+), without e12: fixed: T1 (\S+), T2 (\S+); V1: T1 (\S+), T2 (\S+); V3: T1 (\S+), T2 (\S+)\s*$")
    for ln in open(OUT296):
        mt = pat.search(ln)
        if mt:
            eps = float(mt.group(1))
            out[eps] = {("fixed", "T1"): E.pnum(mt.group(2)), ("fixed", "T2"): E.pnum(mt.group(3)), ("V1", "T1"): E.pnum(mt.group(4)), ("V1", "T2"): E.pnum(mt.group(5)),
                        ("V3", "T1"): E.pnum(mt.group(6)), ("V3", "T2"): E.pnum(mt.group(7))}
    return out


def parse305():
    """#305's output: M - P4 at X1 to X3 (the table 'M - P4' lines), h = 0.005."""
    out = {}
    pat = re.compile(r"^\s+(X[123]): leverage (\S+); M - P4 (\S+);")
    for ln in open(OUT296 if mut("parse_wrong") else OUT305):
        mt = pat.match(ln)
        if mt:
            out[mt.group(1)] = dict(lev=E.pnum(mt.group(2).rstrip(";")), res=E.pnum(mt.group(3)))
    return out


def leverage(P, train, names, kind="P3"):
    A = E.design(arr(P, train, "g1"), arr(P, train, "rdot"), kind)
    G = np.linalg.inv(A.T @ A)
    At = E.design(arr(P, names, "g1"), arr(P, names, "rdot"), kind)
    return np.einsum("ij,jk,ik->i", At, G, At), np.einsum("ij,jk,ik->i", A, G, A)


# --------------------------------------------------------------------------------------------------------------------
# the git state of the frozen file
# --------------------------------------------------------------------------------------------------------------------
def git_state(path):
    rel = os.path.relpath(path, RAW_REPO)

    def gt(*a):
        return subprocess.run(["git", "-C", RAW_REPO, *a], capture_output=True, text=True)
    tracked = gt("ls-files", "--error-unmatch", rel).returncode == 0
    clean = tracked and gt("status", "--porcelain", "--", rel).stdout.strip() == ""
    log = gt("log", "-1", "--format=%H %cI", "--", rel).stdout.strip()
    hsh = log.split()[0] if log else ""
    pushed = bool(hsh) and gt("branch", "-r", "--contains", hsh).stdout.strip() != ""
    return tracked, clean, pushed, log


# --------------------------------------------------------------------------------------------------------------------
# planted cases
# --------------------------------------------------------------------------------------------------------------------
def planted_set(P, names, off_test=1.0, a_t=1.0, c_t=1.0, test=()):
    """Members that follow the form exactly at the real (gamma_1, r_th, rdot_th) of P (these are the 'true' values), measured with a common additive offset (dr, dv) at every member and at the
    static member B = 0 (r_th = 1, rdot_th = 0, M = M_0). The test members get the offset times off_test, a x a_t, c x c_t in M."""
    a0, b0_, c0 = E.FROZEN_P4
    dr, dv = -1.41e-4, -3.24e-4
    Q = {}
    for n in names:
        p = P[n]
        isT = n in test
        a, c = (a0 * a_t, c0 * c_t) if isT else (a0, c0)
        k = off_test if isT else 1.0
        M = p["rth"] * (M0 + a * p["g1"] + b0_ * p["g1"] ** 2 + c * p["rdot"])
        Q[n] = dict(p, rth=p["rth"] + k * dr, rdot=p["rdot"] + k * dv, M=M)
    Q[REF] = dict(P[REF], g1=0.0, rth=1.0 + dr, rdot=dv, M=M0)
    return Q


def check_planted(P2, train, test):
    """PL1: planted members that follow the form, with a common offset like B = 0's second-order part, through make_fit and the frozen prediction (the code of freeze and analyze) with V3: residual <= 1e-12 and the
    coefficients recovered to 1e-9. PL2/PL3: a x 1.05, or c x 0.8, in the test members' M fires (bar 6e-5). PL4: an offset at the test members 1.5 times B = 0's (a member dependence of the second-order part) fires.
    PL296: #296's planted cases (departure time, throat speed through member_at and speed_at at eps = 1e-2, the unsubtracted form) pass unchanged."""
    names = tuple(train) + tuple(test)
    m0 = E.m0_const(H)
    a0, b0_, c0 = E.FROZEN_P4

    def run(**kw):
        Q = planted_set(P2, names, test=test, **kw)
        Qr = {REF: Q[REF], **{n: Q[n] for n in names}}
        f = make_fit(Qr, train, EPS_B, True, "P3", m0)
        res, _, _ = resid(f, Qr, test)
        return f, float(np.max(np.abs(res)))
    f1, own = run()
    dco = max(abs(x - y) for x, y in zip(f1["co"], (a0, b0_, c0))) if len(f1["co"]) == 3 else float("inf")
    check("PL1", "planted members that follow the form exactly, offset by B = 0's second-order part (r_th -1.41e-4, rdot_th -3.24e-4), through the fit and the frozen prediction with V3: residual <= 1e-12, coefficients to 1e-9",
          own <= 1e-12 and dco <= 1e-9, f"max |M - predicted| {own:.1e}; max |coefficient - P4's| {dco:.1e}")
    _, va = run(a_t=1.05)
    check("PL2", "the same with a x 1.05 in the test members' M fires (bar 6e-5)", va > BAR, f"{va:.1e}")
    _, vc = run(c_t=0.8)
    check("PL3", "the same with c x 0.8 in the test members' M fires (bar 6e-5)", vc > BAR, f"{vc:.1e}")
    _, vo = run(off_test=1.0 if mut("no_member_dep") else 1.5)
    check("PL4", "the same with the offset at the test members 1.5 times B = 0's (a member dependence of the second-order part) fires (bar 6e-5)", vo > BAR, f"{vo:.1e}")
    E.results.clear()
    old_sets = E.sets
    if PART == "dev":
        E.sets = lambda: (train, test, train, test)
    try:
        E.check_planted()
    finally:
        E.sets = old_sets
    bad = [c for c, ok, _ in E.results if ok is not True]
    check("PL296", "#296's planted cases (PDEP2, PTD2, PFIT2) pass through the same readings at eps = 1e-2 (the sets patched to this part's)", not bad and len(E.results) == 3, f"{len(E.results)} run; not passing {bad}")


# --------------------------------------------------------------------------------------------------------------------
# the analysis
# --------------------------------------------------------------------------------------------------------------------
def analyze_cmd():
    t0 = time.time()
    train, test = split()
    names = (REF,) + train + test + (N3 if PART == "record" else ())
    path = frozen_path()
    fz = json.load(open(path))
    if mut("stale_frozen"):
        fz["fits"][GATING_KEY]["co"][0] += 1e-6
    spec_bad = [n for n in XN if list(XSPEC[n]) != list(dict((nm, s) for nm, s, _ in s4.POINTS4)[n])]
    print(f"== Part {PART}: train {len(train)} {train}; test {test}; frozen file {path}")
    check("CFG", "M_0 = 0.223086 equals B = 0's plateau (h = 0.005) to 1e-6; P4's constants equal #295's frozen refit; the bar is 6e-5; the frozen file is this part's, over this training set, with M_0 and bar as declared; "
          "X1 to X3 are the declared (B, x_c, w) of #305",
          abs(E.m0_const(H) - E.read(REF, EPS_A)["M"]) <= 1e-6 and max(abs(x - y) for x, y in zip(s4.CRIT3["P4"], E.FROZEN_P4)) <= 1e-12 and (1e-4 if mut("bar_shift") else BAR) == 6e-5
          and fz["part"] == PART and fz["train"] == list(train) and fz["m0"] == M0 and fz["bar"] == 6e-5 and XN == ("X1", "X2", "X3")
          and [list(XSPEC[n]) for n in XN] == [[0.09112, 4.0, 1.5], [0.1695, 4.0, 1.5], [0.2044, 4.0, 1.5]] and not spec_bad,
          f"|M_0 - B = 0 plateau| {abs(E.m0_const(H) - E.read(REF, EPS_A)['M']):.1e}; frozen part {fz['part']}, {len(fz['train'])} training points, M_0 {fz['m0']}, bar {fz['bar']}")
    tracked, clean, pushed, log = git_state(path)
    if mut("git_untracked"):
        tracked, clean, pushed, log = git_state(path + ".none")
    check("GIT", "the frozen file is tracked, unmodified and contained in a pushed commit of the raw repository, before X1 to X3 are read", tracked and clean and pushed, f"tracked {tracked}, unmodified {clean}, pushed {pushed}, last commit {log}")
    if PART == "record":
        seed_x()
    miss = {}
    P1, P2 = {}, {}
    for n in names:
        for eps, P in ((EPS_A, P1), (EPS_B, P2)):
            try:
                P[n] = rd(n, eps)
                if not finite(P[n]):
                    miss[f"{n} {eps:g}"] = "not finite"
            except Exception as ex:
                miss[f"{n} {eps:g}"] = repr(ex)
    check("RUNS", "every member and B = 0 reads without error at both eps, with finite (t_dep, gamma_1, r_th, rdot_th, M)", not miss, f"missing {miss}" if miss else "")
    if miss:
        return finish(t0)

    # the loader of the test points, against the loader of #296 on points that both read
    chk = [n for n in ("W2", "C2", "W1") if n in train + test][:2] if PART == "dev" else list(N3)
    bad = []
    for n in chk:
        spec, root = (E.SPEC[n], E.ROOT_OF[n])
        before = {eps: dict(E.read(n, eps)) for eps in (EPS_A, EPS_B)}
        E._E.pop((n, H), None)
        for eps in (EPS_A, EPS_B):
            E._R.pop((n, eps, H), None)
        E._E[(n, H)] = load_stored(root, spec)
        for eps in (EPS_A, EPS_B):
            after = E.read(n, eps)
            if any(after[k] != before[eps][k] for k in ("g1", "rth", "rdot", "M", "t_dep")):
                bad.append(f"{n} {eps:g}")
    check("RDX", "the stored-run loader used for X1 to X3 (settled_mass_route.load_point, the root and spec of the point) returns members bit for bit equal to #296's on points both read", not bad, f"points {chk}; differing {bad}")

    print("\n== The members at departure (h = 0.005)")
    for n in names:
        print(f"    {n:6s} eps 1e-3: g {P1[n]['g1']:+.6f} r_th {P1[n]['rth']:.6f} rdot {P1[n]['rdot']:+.3e} M {P1[n]['M']:.6f} | eps 1e-2: g {P2[n]['g1']:+.6f} r_th {P2[n]['rth']:.6f} rdot {P2[n]['rdot']:+.3e} M {P2[n]['M']:.6f} | t_dep {P1[n]['t_dep']:.3f} {P2[n]['t_dep']:.3f}")

    tab = parse296_table()
    bad, nc = [], 0
    for n in names:
        if n in XN or n not in tab:
            continue
        for key, eps, fld in (("t1", EPS_A, "t_dep"), ("t2", EPS_B, "t_dep"), ("g1_1", EPS_A, "g1"), ("g1_2", EPS_B, "g1"), ("r1", EPS_A, "rth"), ("r2", EPS_B, "rth"), ("v1", EPS_A, "rdot"), ("v2", EPS_B, "rdot"), ("M", EPS_A, "M")):
            v, tol = tab[n][key]
            x = (P1 if eps == EPS_A else P2)[n][fld]
            nc += 1
            if abs(x - v) > tol + 1e-9:
                bad.append(f"{n} {key}: {x:.7g} vs printed {v}")
    check("REP296", "the members' t_dep, gamma_1, r_th, rdot_th at eps = 1e-3 and 1e-2 and M equal #296's printed table to half a unit of the last printed digit", not bad and nc == 9 * len([n for n in names if n in tab]) and nc > 0, f"{nc} numbers compared; differing {bad[:5]}")

    # frozen against recomputed
    m0 = E.m0_const(H)
    rec = make_all(P1, P2, train, m0)
    bad = []
    for key, eps, sub, kind in FITS:
        f, q = fz["fits"][key], rec[key]
        if len(f["co"]) != len(q["co"]) or max(abs(x - y) for x, y in zip(f["co"], q["co"])) > 1e-9:
            bad.append(f"{key} co {q['co']} vs {f['co']}")
        if sub and (abs(f["b0"]["dr"] - q["b0"]["dr"]) > 1e-12 or abs(f["b0"]["dv"] - q["b0"]["dv"]) > 1e-12):
            bad.append(f"{key} offsets")
    hsh_ok = fz["inputs_sha256"] == inputs_hash({n: P1[n] for n in (REF,) + train}, {n: P2[n] for n in (REF,) + train}, train)
    check("FRZ", "the five frozen fits (coefficients, B = 0's offsets, the sha256 of the inputs) equal the refit recomputed from the stored runs now, to 1e-9", not bad and hsh_ok, f"differing {bad[:3]}; inputs hash equal {hsh_ok}")

    # #296's printed V3 numbers (without e12): T1 (fit TRAIN9 less e12, predict the six) and, at the record, T2 (fit that and the six, predict #295's three)
    pv = parse296_v3()
    tr9 = tuple(n for n in E.TRAIN9 if n != "e12")
    bad, out = [], []
    for eps in (EPS_A, EPS_B):
        P = P1 if eps == EPS_A else P2
        for lab, kind, sub in (("fixed", "P3", False), ("V1", "free", False), ("V3", "P3", True)):
            for tn, fn, tt in (("T1", tr9, E.SIX), ("T2", tr9 + E.SIX, N3)):
                if PART == "dev" and tn == "T2":
                    continue
                Pq = {n: P[n] for n in (REF,) + fn + tt if n in P}
                if tn == "T2" and not all(n in Pq for n in tt):
                    continue
                Pv = apply_offsets(Pq, **offsets_of(Pq[REF])) if sub else Pq
                R = E.fit_pred(Pv, fn, tt, m0, kind)
                v, tol = pv[eps][(lab, tn)]
                x = float(np.max(np.abs(R["res"])))
                out.append(f"{lab} {tn} {eps:g}: {x:.3e} ({v:.2e})")
                if abs(x - v) > tol + 1e-12:
                    bad.append(f"{lab} {tn} eps {eps:g}: {x:.3e} vs printed {v:.2e}")
    check("REPV3", "#296's printed max |M - predicted| without e12 (fixed M_0, V1, V3; T1 over the nine less e12, T2 over those and the six at the record) at eps = 1e-3 and 1e-2 equal the recomputed values to half a unit of the last digit",
          not bad and len(out) == (6 if PART == "dev" else 12), f"{len(out)} compared; differing {bad[:4]}")
    print("    " + "; ".join(out[:12]))

    p305 = parse305()
    check("PARSE305", "#305's printed M - P4 (+9.00e-05, +2.13e-04, +1.98e-04) and leverage are read from its output", [round(p305[n]["res"][0] * 1e6) for n in XN] == [90, 213, 198] if len(p305) == 3 else False,
          ", ".join(f"{n} {p305[n]['res'][0]:+.2e}" for n in p305))
    # M - P4 at eps = 1e-3 (the frozen refit of #295) for the test points, against the earlier printed residuals
    P4 = dict(eps=EPS_A, kind="P3", sub=False, co=list(E.FROZEN_P4), b0=None, m0=m0)
    res4, _, _ = resid(P4, P1, test)
    if PART == "record":
        bad = [f"{n} {x:+.3e} vs {p305[n]['res'][0]:+.2e}" for n, x in zip(test, res4) if abs(x - p305[n]["res"][0]) > p305[n]["res"][1] + 1e-12]
        check("REFX", "at eps = 1e-3 the frozen P4 (#295's refit) misses X1 to X3 by #305's printed +9.00e-05, +2.13e-04, +1.98e-04, to half a unit of the last digit", not bad, ", ".join(f"{n} {x:+.3e}" for n, x in zip(test, res4)) + (f"; differing {bad}" if bad else ""))
    else:
        c83 = E.parse283()
        P3f = dict(P4, co=list(E.P3))
        sx = [n for n in test if n in E.SIX]
        resP3, _, _ = resid(P3f, P1, sx)
        bad = [f"{n} {x:+.3e} vs {c83[n]['rP3'][0]:+.2e}" for n, x in zip(sx, resP3) if abs(x - c83[n]["rP3"][0]) > 1e-6]
        check("REFX", "(development stand-in) at eps = 1e-3 the held-out members of #283 (C2) are missed by P3 (#267's rounded constants) by #283's printed M - P3, to 1e-6", not bad and len(sx) >= 1, ", ".join(f"{n} {x:+.3e}" for n, x in zip(sx, resP3)) + (f"; differing {bad}" if bad else ""))

    print("\n== Planted cases")
    check_planted(P2, train, test)

    print(f"\n== The frozen prediction of {', '.join(test)} (h = 0.005, M_0 = {M0}, bar {BAR:g}; M - predicted)")
    lev_t, lev_tr = leverage(apply_offsets(P2, **fz["fits"][GATING_KEY]["b0"]), train, test)
    print(f"    leverage in the frozen V3 design at eps = 1e-2: training max {lev_tr.max():.3f} ({train[int(np.argmax(lev_tr))]}), mean {lev_tr.mean():.3f}; test " + ", ".join(f"{n} {x:.3f}" for n, x in zip(test, lev_t)))
    allres = {}
    for key, eps, sub, kind in FITS:
        P = P2 if eps == EPS_B else P1
        f = fz["fits"][key]
        res, pred, Pv = resid(f, P, test)
        allres[key] = res
        co = f["co"]
        print(f"    {key:9s} (a, b, c{', icpt' if len(co) > 3 else ''}) = {', '.join(f'{x:+.6f}' for x in co)}; in-sample {f['ins_max']:.2e}, loo {f['loo_max']:.2e}: " + ", ".join(f"{n} {x:+.3e}" for n, x in zip(test, res)) + f"; max {np.max(np.abs(res)):.3e}")
    res = allres[GATING_KEY]
    pred = resid(fz["fits"][GATING_KEY], P2, test)[1]
    Pv = apply_offsets(P2, **fz["fits"][GATING_KEY]["b0"])
    for i, n in enumerate(test):
        print(f"      {n}: gamma_1 {P2[n]['g1']:.6f}, r_th' {Pv[n]['rth']:.6f}, rdot_th' {Pv[n]['rdot']:+.4e}, M {P2[n]['M']:.6f}, predicted {pred[i]:.6f}, M - predicted {res[i]:+.3e}"
              f", (M - predicted)/rdot_th' {res[i] / Pv[n]['rdot']:+.4f} (c = {fz['fits'][GATING_KEY]['co'][2]:.4f})")
    check("PX", f"the frozen V3 fit at eps = 1e-2 (B = 0's reading subtracted, M_0 = 0.223086, coefficients over {len(train)} points) predicts {', '.join(test)} within 6e-5", float(np.max(np.abs(res))) <= BAR,
          f"max |M - predicted| {np.max(np.abs(res)):.2e} (bar {BAR:g}); " + ", ".join(f"{n} {x:+.2e}" for n, x in zip(test, res)))
    for key, lab in (("V1_1e-2", "V1 (free intercept) at eps = 1e-2"), ("V3_1e-3", "V3 at eps = 1e-3"), ("FIX_1e-2", "the unsubtracted refit at eps = 1e-2"), ("FIX_1e-3", "the unsubtracted refit at eps = 1e-3 (#305's setting, e12 dropped)")):
        x = allres[key]
        check(f"PX-{key}", f"{lab}: {', '.join(test)} within 6e-5", float(np.max(np.abs(x))) <= BAR, f"max {np.max(np.abs(x)):.2e}; " + ", ".join(f"{n} {y:+.2e}" for n, y in zip(test, x)), ph=True)
    if PART == "record":
        resN, _, _ = resid(fz["fits"][GATING_KEY], P2, N3)
        print(f"    control, not X: the frozen V3 fit at eps = 1e-2 on #295's N1 to N3: " + ", ".join(f"{n} {x:+.2e}" for n, x in zip(N3, resN)) + f" (#296's T2 V3 without e12: 2.12e-06)")
    return finish(t0)


def post_cmd():
    """Post hoc (after the record; gates nothing): the shift of one reading that would close the miss at X, the same procedure at other eps, and refits that include X (in-sample)."""
    t0 = time.time()
    train, test = FIT, XN
    fz = json.load(open(frozen_path()))
    seed_x()
    names = (REF,) + train + test + N3
    P2 = readset(names, EPS_B)
    f = fz["fits"][GATING_KEY]
    res, pred, Pv = resid(f, P2, test)
    a, b_, c = f["co"]
    dr0, dv0 = f["b0"]["dr"], f["b0"]["dv"]
    print(f"== Post hoc 1 (gates nothing): the shift of one reading that would close the frozen V3 miss at X (eps = 1e-2), against B = 0's second-order part (dr {dr0:+.3e}, dv {dv0:+.3e})")
    print("    point  M - pred    leverage | needed shift of rdot_th' (frac of |rdot_th'|, / B=0's dv) | of r_th' (/ B=0's dr) | of gamma_1 (/ gamma_1)")
    lev_t, _ = leverage(Pv, train, test)
    for i, n in enumerate(test):
        y = pred[i] / Pv[n]["rth"]
        dv_s = res[i] / (Pv[n]["rth"] * c)
        dr_s = res[i] / y
        dg_s = res[i] / (Pv[n]["rth"] * (a + 2 * b_ * Pv[n]["g1"]))
        print(f"    {n}   {res[i]:+.3e}  {lev_t[i]:.3f}  | {dv_s:+.3e} ({dv_s / abs(Pv[n]['rdot']):+.1%}, {dv_s / dv0:+.2f}) | {dr_s:+.3e} ({dr_s / dr0:+.2f}) | {dg_s:+.3e} ({dg_s / Pv[n]['g1']:+.1%})")
    print("\n== Post hoc 2 (gates nothing; X was read at these eps after the record, so no eps may be chosen from this table): the frozen procedure repeated at other eps (V3 and unsubtracted refit over the 14, predict X and N)")
    print("    eps     | V3: X1, X2, X3 | max X, max N | unsubtracted: X1, X2, X3 | max X, max N | in-sample (V3, fixed) | rdot_th at X1, X2, X3 | B = 0's dr, dv")
    for eps in (EPS_A, 2e-3, 5e-3, EPS_B, 2e-2, 5e-2):
        try:
            P = readset(names, eps)
            row = []
            ref = P[REF]
            Pq = apply_offsets(P, **offsets_of(ref))
            R3x = E.fit_pred(Pq, train, test, M0)
            R3n = E.fit_pred(Pq, train, N3, M0)
            Rfx = E.fit_pred(P, train, test, M0)
            Rfn = E.fit_pred(P, train, N3, M0)
            print(f"    {eps:<7g} | " + " ".join(f"{x:+.2e}" for x in R3x["res"]) + f" | {np.max(np.abs(R3x['res'])):.2e} {np.max(np.abs(R3n['res'])):.2e} | " + " ".join(f"{x:+.2e}" for x in Rfx["res"])
                  + f" | {np.max(np.abs(Rfx['res'])):.2e} {np.max(np.abs(Rfn['res'])):.2e} | {np.max(np.abs(R3x['ins'])):.1e} {np.max(np.abs(Rfx['ins'])):.1e} | " + " ".join(f"{P[n]['rdot']:+.2e}" for n in test)
                  + f" | {ref['rth'] - 1:+.2e} {ref['rdot']:+.2e}")
        except Exception as ex:
            print(f"    {eps:<7g} | not computable: {ex!r}")
    print("\n== Post hoc 3 (gates nothing): refits that include X (descriptions on points that include the tested ones), V3 at eps = 1e-2")
    Pq = apply_offsets(P2, **offsets_of(P2[REF]))
    for lab, tr in (("14 + X1-X3 (17 points)", train + test), ("14 + N1-N3 (17 points)", train + N3), ("14 + N + X (20 points)", train + N3 + test)):
        R = E.fit_pred(Pq, tr, tr[:1], M0)
        print(f"    {lab}: (a, b, c) = {', '.join(f'{x:+.4f}' for x in R['co'])}; in-sample max {np.max(np.abs(R['ins'])):.2e}, leave-one-out max {np.max(np.abs(R['loo'])):.2e}; in-sample at X: "
              + ", ".join(f"{n} {R['ins'][tr.index(n)]:+.2e}" for n in test if n in tr) + "; leave-one-out at X: " + ", ".join(f"{n} {R['loo'][tr.index(n)]:+.2e}" for n in test if n in tr))
    for lab, tr in (("14 + X1-X3 (17 points), unsubtracted", train + test),):
        P1_ = P2
        R = E.fit_pred(P1_, tr, tr[:1], M0)
        print(f"    {lab}: (a, b, c) = {', '.join(f'{x:+.4f}' for x in R['co'])}; in-sample max {np.max(np.abs(R['ins'])):.2e}, leave-one-out max {np.max(np.abs(R['loo'])):.2e}")
    print(f"\n{time.time() - t0:.0f} s")
    return 0


MUTATIONS = {      # name: the checks it newly fails on the development partition (--part=dev --baseline=the unmutated log); declared in the raw README before the freeze
    "asym_sign": ['FRZ', 'PL296', 'PX', 'PX-FIX_1e-3', 'PX-V1_1e-2', 'PX-V3_1e-3', 'REFX', 'REP296'],
    "b0_eps": ['FRZ'],
    "bar_shift": ['CFG'],
    "coef_a": ['PL1', 'PL2', 'PL296', 'PX', 'PX-FIX_1e-3', 'PX-V1_1e-2', 'PX-V3_1e-3', 'REFX', 'REPV3'],
    "coef_c": ['PL1', 'PL296', 'PL3', 'PX', 'PX-FIX_1e-3', 'PX-V1_1e-2', 'PX-V3_1e-3', 'REFX', 'REPV3'],
    "dep_eps": ['FRZ', 'PL296', 'PX', 'PX-FIX_1e-3', 'PX-V1_1e-2', 'PX-V3_1e-3', 'REFX', 'REP296', 'REPV3'],
    "eps_swap": ['FRZ', 'PX', 'PX-V1_1e-2', 'REP296', 'REPV3'],
    "git_untracked": ['GIT'],
    "keep_e12": ['CFG', 'FRZ'],
    "m0_shift": ['CFG', 'FRZ', 'PL1', 'REFX', 'REPV3'],
    "no_member_dep": ['PL4'],
    "no_rscale": ['PL1', 'PL296', 'PX', 'PX-FIX_1e-3', 'PX-V1_1e-2', 'PX-V3_1e-3', 'REFX', 'REPV3'],
    "no_run": ['RUNS'],
    "no_speed": ['FRZ', 'PL1', 'PL296', 'REPV3'],
    "no_subtract": ['FRZ', 'PL1', 'PX', 'REPV3'],
    "parse_wrong": ['PARSE305'],
    "rdot_nolapse": ['FRZ', 'PL296', 'REFX', 'REP296', 'REPV3'],
    "stale_frozen": ['FRZ'],
    "sub_r_only": ['FRZ', 'PL1', 'PX', 'REPV3'],
    "sub_sign": ['FRZ', 'PL1', 'PX', 'REPV3'],
    "sub_v_only": ['FRZ', 'PL1', 'REPV3'],
    "test_swap": ['PL1', 'PX', 'PX-FIX_1e-3', 'PX-V1_1e-2', 'PX-V3_1e-3'],
    "x_loader_h": ['RDX'],
}


def main():
    global MUT, PART, MODE
    argv = sys.argv[1:]
    cmd = argv[0] if argv else "help"
    opts = dict(a.lstrip("-").split("=", 1) if "=" in a else (a.lstrip("-"), True) for a in argv[1:])
    PART = opts.get("part", "record")
    if PART not in ("record", "dev"):
        raise SystemExit("--part=record|dev")
    MODE = cmd
    if "mutate" in opts:
        MUT = opts["mutate"]
        if MUT not in MUTLIST:
            raise SystemExit(f"unknown mutation {MUT}; known: {', '.join(MUTLIST)}")
        if "baseline" in opts:
            BASE.update(ln.split("]", 1)[1].split(":", 1)[0].strip() for ln in open(opts["baseline"]) if ln.startswith("[FAIL]"))
        if MUT in E.MUTLIST:
            E.MUT = MUT
            r.MUT = b.MUT = m.MUT = g.MUT = E.sm.MUT = E.r3.MUT = MUT
            tl.LMUT = MUT if MUT in ("asym_sign",) else None
    import platform
    import numba as nb
    print(f"settled_mass_blind.py {cmd} {opts}: numba {nb.__version__}, numpy {np.__version__}, python {platform.python_version()}, {platform.node()}, part {PART}, mutation {MUT}, {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    if cmd == "freeze":
        sys.exit(freeze_cmd())
    if cmd == "analyze":
        sys.exit(analyze_cmd())
    if cmd == "post":
        sys.exit(post_cmd())
    raise SystemExit(__doc__)


MUTLIST = ("coef_a", "coef_c", "m0_shift", "no_rscale", "no_speed", "eps_swap", "dep_eps", "rdot_nolapse", "asym_sign", "no_run", "bar_shift", "no_subtract", "sub_sign", "sub_r_only", "sub_v_only",
           "keep_e12", "stale_frozen", "test_swap", "x_loader_h", "git_untracked", "b0_eps", "parse_wrong", "no_member_dep")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Does the settled-mass relation of family G hang on the definition of departure? Stored runs read at eps = 1e-3 and 1e-2 (#296). D = 4, delta = 1e-8, h = 0.005.

#266 reads the settled mass of family G along the normal amplitude B, #267 along delta, #283 along two other routes, #295 along a third. Every member is read where the unstable mode's part
|odd| / r_tuned reaches eps (#266's definition, eps = 1e-3). The stored member rows reach eps = 1e-2 (and 1e-1). This driver reads all eighteen members at both eps with #283's code
(settled_mass_route.py, through #295's settled_mass_route3.py, unchanged), and tests P3's form (y = M / r_th - M_0 = a g + b g^2 + c rdot_th, M_0 = 0.223086, plain least squares) at eps = 1e-2:
T1: refit over the nine points of #267 (four B, five delta), predict #283's six; T2: refit over the fifteen, predict #295's three; bar 6e-5. At eps = 1e-3 it reproduces #283's P3 residuals and
#295's refit and residuals. No new run: analysis of stored runs only.

    python3 settled_mass_eps.py analyze [--part=record|dev] [--mutate=NAME] [--baseline=FILE]

`--part=dev` is the development partition: the nine test points (#283's six, #295's three) are never read at eps != 1e-3 (a guard raises), T1 and T2 are replaced by tests over held-out
training points, and post hoc blocks that need the test points at other eps are skipped. The mutations' declared effects (MUTATIONS) are computed there, before the test points are read at eps = 1e-2.

Units: a = 1 is the throat radius of the initial static wormhole, 8 pi G = 1, D = 4, M = r/2 on the apparent horizon. Status words as in AGENTS.md: every number is numerical evidence.
"""
import json
import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import settled_mass_route3 as r3        # noqa: E402  (#295's driver: imports #283's, #266's, #267's and everything below them unchanged)
r = r3.r
b, m, g, th, tl, sm, dn, td = r.b, r.m, r.g, r.th, r.tl, r.sm, r.dn, r.td

SQ2 = math.sqrt(2.0)
MUT = None
PART = "record"
H_FIN, H_ALT = 0.005, 0.01
EPS_A, EPS_B = 1e-3, 1e-2                       # #266's definition, and the one tested here
EPS_SCAN = (1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2)  # post hoc
M0 = r3.M0                                      # 0.223086 (the coordinator's ruling; B = 0's plateau at h = 0.005)
BAR = 6e-5
BAR_DEV = 3e-4                                      # development partition only: the held-out training points are not predicted to 6e-5 at eps = 1e-2 (development, training points), so mutations are shown at this bar
RAW_ROOT = r3.RAW_ROOT
I266, I267, I283 = r3.I266, r3.I267, r3.I283
I295 = os.path.join(RAW_ROOT, "I295-mac")
OUT266, OUT267 = os.path.join(HERE, "settled_mass_b_output.txt"), os.path.join(HERE, "settled_mass_member_output.txt")
OUT283, OUT295 = os.path.join(HERE, "settled_mass_route_output.txt"), os.path.join(HERE, "settled_mass_route3_output.txt")

TRAIN9 = tuple(f"B{B:g}" for B in r3.B_266) + tuple(f"e{e}" for e in r3.E_267)       # #267's nine: four B, five delta
SIX = tuple(n for n, _ in r3.POINTS_283)                                              # #283's six
THREE = tuple(n for n, _, _ in r3.POINTS3)                                            # #295's three
TEST_NAMES = SIX + THREE
REF = "B0"                                                                            # the symmetric member (M_0's source)
SPEC = {**dict(r3.POINTS_283), **{n: s for n, s, _ in r3.POINTS3}}
ROOT_OF = {n: (I283 if n[0] in "WC" else I295) for n in TEST_NAMES}
DEV_T1 = (("B0.025", "B0.1", "B0.125", "e7", "e10", "e12"), ("B0.05", "e9", "e11"))   # development partition: fit, held out
DEV_T2 = (("B0.025", "B0.05", "B0.1", "B0.125", "e7", "e9", "e11", "e12"), ("e10",))
FROZEN_P4 = (0.2296318924, 0.2834170383, 0.7657657062)                                 # #295's frozen refit (the brief's numbers)
BAR_PTD2_RDOT = 1e-5                                                                  # planted speed at eps = 1e-2 (#283's 5e-6 at eps = 1e-3; measured 8.1e-6 in development)
P3 = r.CRITR["P3"]                                                                    # #267's, rounded as printed (0.2300, 0.2724, 0.7835)


def mut(name):
    return MUT == name


def bar_t():
    return BAR if PART == "record" else BAR_DEV


results = []
BASE = set()


def check(cid, text, ok, detail="", ph=False):
    ok = None if ok is None else bool(ok)
    results.append((cid, ok, ph))
    print(f"[{'SKIP' if ok is None else ('pass' if ok else 'FAIL')}] {cid}: {'POST HOC, gates nothing: ' if ph else ''}{text}" + (f"  ({detail})" if detail else ""), flush=True)


def finish(t0):
    nfail = sum(1 for _, ok, ph in results if ok is False and not ph)
    nph = sum(1 for _, ok, ph in results if ok is False and ph)
    print(f"\n{len(results)} checks ({sum(1 for _, _, ph in results if ph)} post hoc, gating nothing), {nfail} gating checks failed, {nph} post hoc checks failed; {time.time() - t0:.0f} s")
    if MUT is not None:
        fired = sorted({c for c, ok, ph in results if ok is False})
        new = sorted(set(fired) - BASE)
        want = sorted(set(MUTATIONS.get(MUT, [])) - BASE)         # a check that fails without a mutation cannot be newly fired
        print(f"mutation {MUT}: fired {fired}; failing without a mutation (--baseline) {sorted(BASE)}; newly fired {new}, declared {want}: {'MATCH' if new == want else 'DIFFERENT'}")
        return 0 if new == want else 1
    return 1 if nfail else 0


# --------------------------------------------------------------------------------------------------------------------
# reading a member from the stored runs
# --------------------------------------------------------------------------------------------------------------------
_E, _R = {}, {}


def entry(name, h):
    """The stored run of one point at h: the pair's locus, the member rows, the chain's series (read only)."""
    if (name, h) in _E:
        return _E[(name, h)]
    if name[0] == "B":
        base = f"{I266}/runs/G_B{float(name[1:]):g}/h{h:g}"
    elif name[0] == "e":
        base = os.path.join(I267, "runs", m.TAG, f"h{h:g}", f"e{int(name[1:])}")
    else:
        e = r.load_point(os.path.join(ROOT_OF[name], "runs"), SPEC[name], h)
        if "err" in e:
            raise RuntimeError(f"{name} h {h:g}: {e['err']}")
        _E[(name, h)] = e
        return e
    e = dict(h=h, base=base, loc=dict(np.load(base + "/member/locus.npz")), m=json.load(open(base + "/member/member.json")), s=dict(np.load(base + "/chain/series.npz")))
    _E[(name, h)] = e
    return e


def read(name, eps, h=H_FIN):
    """The member at departure (eps), the throat's speed there and the plateau of M_AH. In the development partition the test points are not read at eps != 1e-3."""
    if PART == "dev" and name in TEST_NAMES and eps != EPS_A:
        raise RuntimeError(f"development partition: {name} must not be read at eps = {eps:g}")
    if mut("eps_swap") and eps == EPS_B:
        eps = EPS_A
    key = (name, eps, h)
    if key in _R:
        return _R[key]
    e = entry(name, h)
    mem = b.member_at(e, eps, "m")
    sp = m.speed_at(e["loc"], mem["t_dep"])
    M = g.summary_g(e["s"], h)["plat1"]
    _R[key] = dict(name=name, eps=eps, h=h, g1=mem["g1"], rth=mem["r_th"], rdot=sp, M=M, t_dep=mem["t_dep"])
    return _R[key]


def readall(names, eps, h=H_FIN):
    return {n: read(n, eps, h) for n in names}


def finite(p):
    return all(np.isfinite(p[k]) for k in ("g1", "rth", "rdot", "M", "t_dep"))


# --------------------------------------------------------------------------------------------------------------------
# the form: y = M / r_th - m0 = a g + b g^2 + c rdot (+ an intercept in the variant 'free'); plain least squares
# --------------------------------------------------------------------------------------------------------------------
def m0_const(h):
    return 0.223095 if (mut("m0_shift") and h == H_FIN) else r.CRITR["M0"][h]


def m0_eps(eps, h):
    """Variant: M_0 read at the same eps, B = 0's plateau divided by B = 0's throat radius at eps (= M_0 at eps = 1e-3 to 1e-6)."""
    p = read(REF, eps, h)
    return p["M"] / p["rth"]


def design(g1, rdot, kind):
    cols = [g1, g1 ** 2]
    if kind != "nospeed":
        cols.append(rdot)
    if kind == "free":
        cols.append(np.ones_like(g1))
    return np.stack(cols, 1)


def arr(P, names, k):
    return np.array([P[n][k] for n in names], float)


def fit_pred(P, train, test, m0, kind="P3"):
    """Fit y = M / r_th - m0 over `train`, predict M = r_th (m0 + design . co) at `test`. Returns the coefficients, M - predicted (test), in-sample and leave-one-out M - predicted (train)."""
    kind = "nospeed" if mut("no_speed") else kind
    if mut("drop_point"):
        train = tuple(n for n in train if n != "e12")
    A = design(arr(P, train, "g1"), arr(P, train, "rdot"), kind)
    rth, M = arr(P, train, "rth"), arr(P, train, "M")
    y = M / rth - m0
    co = np.linalg.lstsq(A, y, rcond=None)[0]
    ins = M - rth * (m0 + A @ co)
    loo = []
    for k in range(len(y)):
        idx = [j for j in range(len(y)) if j != k]
        c_k = np.linalg.lstsq(A[idx], y[idx], rcond=None)[0]
        loo.append(M[k] - rth[k] * (m0 + A[k] @ c_k))
    loo = ins.copy() if mut("loo_off") else np.array(loo)
    cp = co.copy()
    if mut("coef_a"):
        cp[0] *= 1.05
    if mut("coef_c") and kind != "nospeed":
        cp[2] *= 0.8
    At = design(arr(P, test, "g1"), arr(P, test, "rdot"), kind)
    rs = np.ones(len(test)) if mut("no_rscale") else arr(P, test, "rth")
    res = arr(P, test, "M") - rs * (m0 + At @ cp)
    return dict(co=co, res=res, ins=ins, loo=loo, train=train, test=test, A=A)


def sub_ref(P, ref):
    """Variant V3: the symmetric member's reading at the same eps (r_th - 1 and rdot_th, which are 0 for the static member) subtracted from every member's."""
    return {n: dict(p, rth=p["rth"] - (ref["rth"] - 1.0), rdot=p["rdot"] - ref["rdot"]) for n, p in P.items()}


def sets():
    """(T1 fit, T1 test, T2 fit, T2 test) of the partition."""
    if PART == "dev":
        return DEV_T1[0], DEV_T1[1], DEV_T2[0], DEV_T2[1]
    return TRAIN9, SIX, TRAIN9 + SIX, THREE


def fmt(x, f="+.2e"):
    return format(x, f)


# --------------------------------------------------------------------------------------------------------------------
# printed numbers of the earlier notes (the reproduction checks read these files)
# --------------------------------------------------------------------------------------------------------------------
def pnum(s):
    """(value, half a unit of the last printed digit) of a printed number."""
    s = s.strip()
    v = float(s)
    if "e" in s.lower():
        mant, ex = re.split("[eE]", s)
        dec = len(mant.split(".")[1]) if "." in mant else 0
        return v, 0.5 * 10.0 ** (int(ex) - dec)
    dec = len(s.split(".")[1]) if "." in s else 0
    return v, 0.5 * 10.0 ** (-dec)


def parse266():
    """#266's output: member at departure at h = 0.005, eps = 1e-3 and 1e-2 (t_dep, r_th, gamma_1) of B = 0, 0.025, 0.05, 0.1, 0.125."""
    out = {}
    pat = re.compile(r"^\s+G_B(\S+)\s+h 0\.005: t_dep (\S+): r_th (\S+), s_th \S+, gamma_1 (\S+)\s+\| 1e-2: t_dep (\S+), r_th (\S+), gamma_1 (\S+)\s+\|")
    for ln in open(OUT266):
        mt = pat.match(ln)
        if mt:
            out["B" + mt.group(1)] = {EPS_A: tuple(pnum(x.rstrip(",:")) for x in mt.group(2, 3, 4)), EPS_B: tuple(pnum(x.rstrip(",:")) for x in mt.group(5, 6, 7))}
    return out


def parse267():
    """#267's output at h = 0.005: delta = 1e-e at eps = 1e-3 and 1e-2 (t_dep, r_th, gamma_1); the post hoc nine-point fit with its coefficients and leave-one-out errors."""
    rows, fit = {}, None
    pat = re.compile(r"^\s+0\.005\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\|\s+(\S+)\s+(\S+)\s+(\S+)\s+\|")
    for ln in open(OUT267):
        mt = pat.match(ln)
        if mt:
            rows[f"e{int(mt.group(1))}"] = {EPS_A: tuple(pnum(x) for x in mt.group(2, 3, 4)), EPS_B: tuple(pnum(x) for x in mt.group(5, 6, 7))}
        mt = re.search(r"fit over these 9 points, y = a g \+ b g\^2 \+ c rdot: coefficients (\S+), (\S+), (\S+); residual max (\S+); leave-one-out errors (.*)$", ln)
        if mt:
            fit = dict(co=[pnum(x.rstrip(",;")) for x in mt.group(1, 2, 3)], loo=[pnum(x.strip()) for x in mt.group(5).split(",")])
    return rows, fit


def parse283():
    """#283's output: the six members' M and M - P3 (h = 0.005)."""
    out = {}
    pat = re.compile(r"^\s+([WC][123])\s+True\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\|")
    for ln in open(OUT283):
        mt = pat.match(ln)
        if mt:
            out[mt.group(1)] = dict(M=pnum(mt.group(2)), rP3=pnum(mt.group(6)))
    return out


def parse295():
    """#295's output: the recomputed refit (a, b, c) and the three new members' M and M - P4 (h = 0.005)."""
    out, co = {}, None
    pat = re.compile(r"^\s+(N[123])\s+(\S+)\s+(\S+)\s+(\S+)\s+\|")
    for ln in open(OUT295):
        mt = pat.match(ln)
        if mt:
            out[mt.group(1)] = dict(M=pnum(mt.group(2)), rP4=pnum(mt.group(4)))
        mt = re.search(r"recomputed \(a, b, c\) = (\S+), (\S+), (\S+);", ln)
        if mt:
            co = [pnum(x) for x in mt.group(1, 2, 3)]
    return out, co


# --------------------------------------------------------------------------------------------------------------------
# planted cases
# --------------------------------------------------------------------------------------------------------------------
def check_planted():
    """PDEP2: a planted growing pair through the departure time at eps = 1e-2. PTD2: planted pairs (the exact static member along a planted trajectory) through member_at and speed_at at eps = 1e-2.
    PFIT2: planted data that follow the refit's form exactly, through the code that tests T1 and T2: they pass (1e-12), and with a x 1.05 or c x 0.8 in the test points they fire (bar 6e-5)."""
    C = r.CRITR
    lam, q0, h_ = 1.18, 1e-9, 0.01 / SQ2
    t = np.arange(0.0, 30.0, h_)
    rt = 1.0 + 0.003 * np.tanh(t / 5.0)
    worst = abs(b.departure_time(t, rt, q0 * np.exp(lam * t) * rt, EPS_B) - math.log(EPS_B / q0) / lam)
    check("PDEP2", "planted growing pair (odd part q0 e^{lambda t}, known crossing time of 1e-2) through departure_time", worst <= C["bar_dep"], f"|t_dep - exact| {worst:.1e} (bar {C['bar_dep']:g})")
    old = C["eps_dep"]
    C["eps_dep"] = EPS_B
    try:
        worst = [0.0, 0.0, 0.0, 0.0]
        for g0, r0, v in r.PLANTED_ENT:
            for h in r.H2:
                en = r.planted_entry(h, g0, r0, v)
                err = (abs(en["mem"]["t_dep"] - en["t_true"]), abs(en["mem"]["g1"] - en["g_true"]), abs(en["mem"]["r_th"] - en["r_true"]), abs(en["sp"] - en["v_true"]))
                worst = [max(a_, b_) for a_, b_ in zip(worst, err)]
                print(f"    PTD2 gamma_1 {g0:g} r_th {r0:g} v {v:+.1e} h {h:g}: t_dep {en['mem']['t_dep']:.4f} ({en['t_true']:.4f}), gamma_1 {en['mem']['g1']:.6f} ({en['g_true']:.6f}), "
                      f"r_th {en['mem']['r_th']:.6f} ({en['r_true']:.6f}), rdot_th {en['sp']:+.3e} ({en['v_true']:+.3e})")
    finally:
        C["eps_dep"] = old
    check("PTD2", "planted pairs (odd part Q0 delta e^{lambda t}, trajectories of the measured size, rows every 0.5 in t from the exact static member) through member_at and speed_at at eps = 1e-2: "
          "departure time, gamma_1, r_th and rdot_th at departure", worst[0] <= C["bar_ptd_t"] and worst[1] <= C["bar_ptd_g"] and worst[2] <= C["bar_ptd_r"] and worst[3] <= BAR_PTD2_RDOT,
          f"max error t_dep {worst[0]:.1e} (bar {C['bar_ptd_t']:g}), gamma_1 {worst[1]:.1e} ({C['bar_ptd_g']:g}), r_th {worst[2]:.1e} ({C['bar_ptd_r']:g}), rdot_th {worst[3]:.1e} ({BAR_PTD2_RDOT:g})")
    f1, t1, f2, t2 = sets()
    names = tuple(dict.fromkeys(f1 + t1 + f2 + t2))
    P = readall(names, EPS_B)
    a_, b_, c_ = FROZEN_P4
    m0 = m0_const(H_FIN)

    def planted(a, c):
        return {n: dict(P[n], M=P[n]["rth"] * (m0 + a * P[n]["g1"] + b_ * P[n]["g1"] ** 2 + c * P[n]["rdot"])) for n in names}
    ok_own, fire_a, fire_c = 0.0, 0.0, 0.0
    for (fit_n, test_n) in ((f1, t1), (f2, t2)):
        own = planted(a_, c_)
        ok_own = max(ok_own, float(np.max(np.abs(fit_pred(own, fit_n, test_n, m0)["res"]))))
        for tag, pa, pc in (("a", 1.05 * a_, c_), ("c", a_, 0.8 * c_)):
            Pw = dict(planted(a_, c_))
            wr = planted(pa, pc)
            for n in test_n:
                Pw[n] = wr[n]
            v = float(np.max(np.abs(fit_pred(Pw, fit_n, test_n, m0)["res"])))
            if tag == "a":
                fire_a = max(fire_a, v)
            else:
                fire_c = max(fire_c, v)
    check("PFIT2", "planted data that follow the form exactly (M = r_th (M_0 + a g + b g^2 + c rdot) on the real eps = 1e-2 (gamma_1, r_th, rdot_th)) pass the T1 and T2 code at 1e-12; with a x 1.05, or c x 0.8, in the test points it fires (bar 6e-5)",
          ok_own <= 1e-12 and fire_a > BAR and fire_c > BAR, f"max |M - predicted|: own {ok_own:.1e}; a x 1.05 {fire_a:.1e}; c x 0.8 {fire_c:.1e} (bar {BAR:g})")


# --------------------------------------------------------------------------------------------------------------------
# the analysis
# --------------------------------------------------------------------------------------------------------------------
def post_hoc(P1, P2, allp, f1, t1, f2, t2):
    """The record-only post hoc blocks (they read the test points at eps != 1e-3): printed, gate nothing."""
    print(f"\n== Post hoc (gates nothing): the eps = 1e-3 constants applied to the eps = 1e-2 members, without a refit (P3 #267's rounded constants on the six; P4 #295's frozen refit on the three; both on the nine and fifteen)")
    for nm, co, nms in (("P3 on the six", P3, SIX), ("P4 on the three", FROZEN_P4, THREE), ("P4 on the nine", FROZEN_P4, TRAIN9), ("P4 on the six", FROZEN_P4, SIX)):
        res = [P2[n]["M"] - P2[n]["rth"] * (M0 + co[0] * P2[n]["g1"] + co[1] * P2[n]["g1"] ** 2 + co[2] * P2[n]["rdot"]) for n in nms]
        res1 = [P1[n]["M"] - P1[n]["rth"] * (M0 + co[0] * P1[n]["g1"] + co[1] * P1[n]["g1"] ** 2 + co[2] * P1[n]["rdot"]) for n in nms]
        print(f"    {nm}: eps = 1e-2: " + ", ".join(f"{n} {x:+.2e}" for n, x in zip(nms, res)) + f"; max {max(abs(x) for x in res):.2e} | eps = 1e-3: max {max(abs(x) for x in res1):.2e}")

    print("\n== Post hoc (gates nothing): the refit's own errors at both eps (fifteen points, then eighteen)")
    for nm, names in (("fifteen", TRAIN9 + SIX), ("eighteen", TRAIN9 + SIX + THREE)):
        for eps, P_ in ((EPS_A, P1), (EPS_B, P2)):
            R = fit_pred(P_, names, names[:1], M0)
            print(f"    {nm} at eps = {eps:g}: (a, b, c) = {', '.join(f'{x:+.6f}' for x in R['co'])}; in-sample max {np.max(np.abs(R['ins'])):.2e}, leave-one-out max {np.max(np.abs(R['loo'])):.2e}")
            print("        in-sample (x1e-6): " + ", ".join(f"{n} {x * 1e6:+.0f}" for n, x in zip(names, R["ins"])))
            print("        leave-one-out (x1e-6): " + ", ".join(f"{n} {x * 1e6:+.0f}" for n, x in zip(names, R["loo"])))

    print("\n== Post hoc (gates nothing): the same moment read from two pairs. The pair at delta = 10^-e read at eps = 1e-2 and the pair at 10^-(e+1) read at eps = 1e-3 reach it together; their means differ by the unstable mode's second-order part")
    pair_rows = []
    for n, nxt in (("e7", "B0.1"), ("B0.1", "e9"), ("e9", "e10"), ("e10", "e11"), ("e11", "e12")):
        a_, b_ = P2[n], P1[nxt]
        pair_rows.append((n, nxt, b_["rth"] - a_["rth"], b_["g1"] - a_["g1"], b_["rdot"] - a_["rdot"], b_["t_dep"] - a_["t_dep"]))
        print(f"    {n}@1e-2 vs {nxt}@1e-3: t_dep {a_['t_dep']:.3f} vs {b_['t_dep']:.3f}; r_th(1e-3 pair) - r_th(1e-2 pair) {pair_rows[-1][2]:+.2e}; gamma_1 {pair_rows[-1][3]:+.2e}; rdot_th {pair_rows[-1][4]:+.2e}")
    print(f"    B = 0 (the symmetric member, r_th = 1): r_th - 1 = {P1[REF]['rth'] - 1:+.2e} at eps = 1e-3, {P2[REF]['rth'] - 1:+.2e} at 1e-2; M_0(eps) = M/r_th = {m0_eps(EPS_A, H_FIN):.7f} and {m0_eps(EPS_B, H_FIN):.7f} against 0.223086")

    print("\n== Post hoc (gates nothing): leverage of the test points in the training design at eps = 1e-2 (design columns g, g^2, rdot)")
    for cid, (fn, tn) in (("T1", (f1, t1)), ("T2", (f2, t2))):
        A = design(arr(P2, fn, "g1"), arr(P2, fn, "rdot"), "P3")
        G = np.linalg.inv(A.T @ A)
        lev = np.einsum("ij,jk,ik->i", A, G, A)
        At = design(arr(P2, tn, "g1"), arr(P2, tn, "rdot"), "P3")
        lt = np.einsum("ij,jk,ik->i", At, G, At)
        print(f"    {cid}: training leverage max {lev.max():.3f} ({fn[int(np.argmax(lev))]}), mean {lev.mean():.3f}; test points: " + ", ".join(f"{n} {x:.3f}" for n, x in zip(tn, lt))
              + f"; singular values {', '.join(f'{x:.3e}' for x in np.linalg.svd(A, compute_uv=False))}")

    print("\n== Post hoc (gates nothing): the tests at other eps (the stored rows reach 1e-1): max |M - predicted| of T1 and T2, three choices of M_0, and the coefficients")
    print("    eps     | T1: M_0 = 0.223086, M_0(eps), free, B = 0 subtracted | T2: the same four | (a, b, c) of the nine | (a, b, c) of the fifteen (M_0 = 0.223086)")
    for eps in EPS_SCAN:
        Pe = readall(allp, eps)
        out = []
        for fn, tn in ((TRAIN9, SIX), (TRAIN9 + SIX, THREE)):
            x1 = fit_pred(Pe, fn, tn, M0)
            x2 = fit_pred(Pe, fn, tn, m0_eps(eps, H_FIN))
            x3 = fit_pred(Pe, fn, tn, M0, "free")
            x4 = fit_pred(sub_ref(Pe, Pe[REF]), fn, tn, M0)
            out.append((x1, x2, x3, x4))
        print(f"    {eps:<7g} | " + " ".join(f"{np.max(np.abs(x['res'])):.2e}" for x in out[0]) + " | " + " ".join(f"{np.max(np.abs(x['res'])):.2e}" for x in out[1])
              + " | " + ", ".join(f"{x:+.4f}" for x in out[0][0]["co"]) + " | " + ", ".join(f"{x:+.4f}" for x in out[1][0]["co"]))
        if eps in (2e-3, 5e-3, 2e-2, 5e-2):
            print(f"            T1 residuals (M_0 = 0.223086): " + ", ".join(f"{n} {x:+.2e}" for n, x in zip(SIX, out[0][0]["res"])) + "; T2: " + ", ".join(f"{n} {x:+.2e}" for n, x in zip(THREE, out[1][0]["res"])))


def analyze_cmd():
    t0 = time.time()
    f1, t1, f2, t2 = sets()
    allp = (REF,) + TRAIN9 + TEST_NAMES
    names_b = allp if PART == "record" else (REF,) + TRAIN9           # points read at eps = 1e-2 (the development partition keeps the test points at eps = 1e-3)
    miss = {}
    for n in allp:
        for eps in (EPS_A, EPS_B):
            if PART == "dev" and eps != EPS_A and n not in names_b:
                continue
            try:
                p = read(n, eps)
                if not finite(p):
                    miss[f"{n} eps {eps:g}"] = "not finite"
            except Exception as ex:
                miss[f"{n} eps {eps:g}"] = repr(ex)
    if mut("no_run"):
        miss["W1 eps 0.01"] = "treated as missing" if PART == "record" else None
        miss["e9 eps 0.01"] = "treated as missing"
    miss = {k: v for k, v in miss.items() if v}
    print(f"== Part {PART}: {len(allp)} points read at eps = {EPS_A:g} and {len(names_b)} at {EPS_B:g}, h = {H_FIN:g}; fit sets T1 {len(f1)} -> {len(t1)}, T2 {len(f2)} -> {len(t2)}")
    if PART == "dev":
        print(f"    development partition: T1 fits {f1} and holds out {t1}; T2 fits {f2} and holds out {t2}; the test points {TEST_NAMES} are read at eps = 1e-3 only")
    spec_bad = [n for n in TEST_NAMES if list(entry(n, H_FIN)["B"]) != list(SPEC[n])]
    ref = read(REF, EPS_A)
    check("CFG", "M_0 = 0.223086 equals B = 0's plateau (h = 0.005) to 1e-6 and 0.223095 at h = 0.01; P4's constants equal the frozen refit of #295; the sets are #267's nine (four B, five delta), #283's six and #295's three; "
          "the (B, xc, w) stored in each run equal the declared one; bar 6e-5",
          abs(m0_const(H_FIN) - ref["M"]) <= 1e-6 and abs(r.CRITR["M0"][H_ALT] - read(REF, EPS_A, H_ALT)["M"]) <= 1e-6 and max(abs(x - y) for x, y in zip(r3.CRIT3["P4"], FROZEN_P4)) <= 1e-12
          and len(TRAIN9) == 9 and len(SIX) == 6 and len(THREE) == 3 and not spec_bad and BAR == 6e-5,
          f"|M_0 - B = 0 plateau| {abs(m0_const(H_FIN) - ref['M']):.1e}; B = 0 plateau {ref['M']:.7f}; specs differing {spec_bad}")
    check("RUNS", "every point reads without error at both eps, with finite (t_dep, gamma_1, r_th, rdot_th, M)", not miss, f"missing {miss}" if miss else "")
    if miss:
        return finish(t0)
    print("\n== Planted cases")
    check_planted()

    print(f"\n== The members at departure at eps = {EPS_A:g} and {EPS_B:g} (h = {H_FIN:g}; #266's reading on the pair's mean profile; rdot_th = #267's speed there)")
    P1 = readall(allp, EPS_A)
    P2 = readall(names_b, EPS_B)
    print("    point   t_dep(1e-3) t_dep(1e-2)  shift | gamma_1(1e-3) gamma_1(1e-2)  d    | r_th(1e-3)  r_th(1e-2)   d         | rdot(1e-3)  rdot(1e-2)  d          | M")
    for n in allp:
        if n not in P2:
            continue
        a_, b_ = P1[n], P2[n]
        print(f"    {n:7s} {a_['t_dep']:9.3f} {b_['t_dep']:9.3f} {b_['t_dep'] - a_['t_dep']:7.3f} | {a_['g1']:+.6f} {b_['g1']:+.6f} {b_['g1'] - a_['g1']:+.1e} | "
              f"{a_['rth']:.6f} {b_['rth']:.6f} {b_['rth'] - a_['rth']:+.1e} | {a_['rdot']:+.3e} {b_['rdot']:+.3e} {b_['rdot'] - a_['rdot']:+.1e} | {a_['M']:.6f}")
    sh = [P2[n]["t_dep"] - P1[n]["t_dep"] for n in names_b]
    check("MEM", "the departure at eps = 1e-2 comes 1.9 to 2.2 later than at 1e-3 at every point (the mode grows by ln 10 / lambda = 1.95 to 2.1: the stored linear growth), and every member is finite",
          all(1.9 <= s_ <= 2.2 for s_ in sh) and all(finite(P2[n]) for n in names_b), f"shift {min(sh):.3f} to {max(sh):.3f}")

    print("\n== The printed members of #266 and #267 at h = 0.005 (training points only), to the half unit of their last printed digit")
    p266, (p267, fit267) = parse266(), parse267()
    bad, n_cmp = [], 0
    for src, tab in (("#266", p266), ("#267", p267)):
        for n, d_ in tab.items():
            if n not in P1 or n not in P2:
                continue
            for eps, P_ in ((EPS_A, P1), (EPS_B, P2)):
                for k, (key, scale) in enumerate((("t_dep", 1.0), ("rth", 1.0), ("g1", 1.0))):
                    v, tol = d_[eps][k]
                    x = P_[n][key]
                    n_cmp += 1
                    if abs(x - v) > tol + 1e-9:
                        bad.append(f"{src} {n} eps {eps:g} {key}: {x:.7f} vs printed {v}")
    check("E2B", "the members at eps = 1e-3 and 1e-2 of B = 0, 0.025, 0.05, 0.1, 0.125 (#266's output) and of delta = 1e-7, 1e-9 .. 1e-12 (#267's output) equal the printed t_dep, r_th and gamma_1 to the half unit of the last digit",
          not bad and n_cmp == 60, f"{n_cmp} numbers compared; differing: {bad[:6]}")

    print(f"\n== At eps = {EPS_A:g}: the earlier notes' numbers (P3 as #267's rounded constants, P4 frozen, the refits) reproduced")
    c83, (c295, co295) = parse283(), parse295()
    d_s, rows = [], []
    for n in SIX:
        p = P1[n]
        pred = p["rth"] * (m0_const(H_FIN) + (1.05 if mut("coef_a") else 1.0) * P3[0] * p["g1"] + P3[1] * p["g1"] ** 2 + (0.8 if mut("coef_c") else 1.0) * P3[2] * p["rdot"])
        res = p["M"] - (pred if not mut("no_rscale") else pred / p["rth"])
        pm, tm = c83[n]["M"], c83[n]["rP3"]
        d_s.append(max(abs(res - tm[0]), abs(p["M"] - pm[0])))
        rows.append(f"{n} {res:+.2e} ({tm[0]:+.2e})")
    check("E3S", "at eps = 1e-3 the six points' plateau M and M - P3 (P3 = #267's (0.2300, 0.2724, 0.7835), M_0 = 0.223086) equal #283's printed values to 1e-6", max(d_s) <= 1e-6, f"max difference {max(d_s):.1e}; " + ", ".join(rows))
    fit15 = fit_pred(P1, TRAIN9 + SIX, THREE, m0_const(H_FIN))
    a15 = fit_pred(P1, TRAIN9 + SIX, TRAIN9 + SIX, m0_const(H_FIN))["co"]
    dco = max(abs(x - y) for x, y in zip(fit15["co"], FROZEN_P4))
    dco_p = max(abs(x - y[0]) for x, y in zip(fit15["co"], co295)) if co295 else float("inf")
    check("E3F", "at eps = 1e-3 the refit over the fifteen points (a, b, c) equals #295's frozen (0.2296318924, 0.2834170383, 0.7657657062) and its printed recomputation to 1e-6",
          dco <= 1e-6 and dco_p <= 1e-6 and len(fit15["train"]) == 15, f"{len(fit15['train'])} points; refit {', '.join(f'{x:+.10f}' for x in fit15['co'])}; max |refit - frozen| {dco:.1e}, vs printed {dco_p:.1e}")
    pr4 = {}
    for n in THREE:
        p = P1[n]
        a_, b_, c_ = (FROZEN_P4[0] * (1.05 if mut("coef_a") else 1.0), FROZEN_P4[1], FROZEN_P4[2] * (0.8 if mut("coef_c") else 1.0))
        pred = (1.0 if mut("no_rscale") else p["rth"]) * (m0_const(H_FIN) + a_ * p["g1"] + b_ * p["g1"] ** 2 + c_ * p["rdot"])
        pr4[n] = p["M"] - pred
    d_n = [max(abs(pr4[n] - c295[n]["rP4"][0]), abs(P1[n]["M"] - c295[n]["M"][0])) for n in THREE]
    check("E3N", "at eps = 1e-3 the three new members' plateau M and M - P4 (frozen refit, M_0 = 0.223086) equal #295's printed values to 1e-6", max(d_n) <= 1e-6,
          f"max difference {max(d_n):.1e}; " + ", ".join(f"{n} {pr4[n]:+.2e} ({c295[n]['rP4'][0]:+.2e})" for n in THREE))
    fit9 = fit_pred(P1, TRAIN9, SIX, m0_const(H_FIN))
    dco9 = [abs(x - y[0]) - y[1] for x, y in zip(fit9["co"], fit267["co"])] if fit267 else [1.0]
    loo_y = np.array([-fit9["loo"][i] / P1[nm]["rth"] for i, nm in enumerate(fit9["train"])]) if len(fit9["train"]) == 9 else np.zeros(9)
    dloo = max(abs(x - y[0]) for x, y in zip(loo_y, fit267["loo"])) if (fit267 and len(loo_y) == len(fit267["loo"])) else float("inf")
    check("E3P", "at eps = 1e-3 the refit over the nine points (a, b, c) equals #267's printed (0.2300, 0.2724, 0.7835) to the half unit of its last digit, and its nine leave-one-out errors (predicted - measured in y, in the printed order) equal #267's to 1e-6",
          max(dco9) <= 1e-9 and dloo <= 1e-6 and len(fit9["train"]) == 9, f"refit {', '.join(f'{x:+.5f}' for x in fit9['co'])}; max leave-one-out difference {dloo:.1e}")

    print(f"\n== The tests at eps = {EPS_B:g} (h = {H_FIN:g}, M_0 = 0.223086, plain least squares, no intercept, no weights; M - predicted, bar {bar_t():g})")
    res_t = {}
    for cid, (fn, tn) in (("T1", (f1, t1)), ("T2", (f2, t2))):
        R2 = fit_pred(P2, fn, tn, m0_const(H_FIN))
        R1 = fit_pred(P1, fn, tn, m0_const(H_FIN))
        res_t[cid] = (R2, R1)
        print(f"    {cid}: fit over {len(R2['train'])} points {R2['train']}; test {tn}")
        print(f"        eps = 1e-2: (a, b, c) = {', '.join(f'{x:+.6f}' for x in R2['co'])}; M - predicted: " + ", ".join(f"{n} {x:+.2e}" for n, x in zip(tn, R2["res"])) + f"; in-sample max {np.max(np.abs(R2['ins'])):.1e}, leave-one-out max {np.max(np.abs(R2['loo'])):.1e}")
        print(f"        eps = 1e-3: (a, b, c) = {', '.join(f'{x:+.6f}' for x in R1['co'])}; M - predicted: " + ", ".join(f"{n} {x:+.2e}" for n, x in zip(tn, R1["res"])) + f"; in-sample max {np.max(np.abs(R1['ins'])):.1e}, leave-one-out max {np.max(np.abs(R1['loo'])):.1e}   (the same refit at #266's definition, unrounded)")
        mx = float(np.max(np.abs(R2["res"])))
        text = ("refit of P3's form over the nine points of #267 (four B, five delta) at eps = 1e-2 predicts #283's six within 6e-5" if cid == "T1" else
                "refit of P3's form over the fifteen points at eps = 1e-2 predicts #295's three within 6e-5")
        if PART == "dev":
            text = f"(development partition) {cid}: fit {len(fn)} training points at eps = 1e-2, predict {tn} within 6e-5"
        check(cid, text, mx <= bar_t(), f"max |M - predicted| {mx:.2e} (bar {bar_t():g}); " + ", ".join(f"{n} {x:+.2e}" for n, x in zip(tn, R2["res"])))

    print(f"\n== Post hoc, declared before the test points were read at eps = 1e-2 (gate nothing): the same tests with M_0 read at the same eps, with a free intercept, and at h = {H_ALT:g}")
    for cid, (fn, tn) in (("T1", (f1, t1)), ("T2", (f2, t2))):
        v2 = fit_pred(P2, fn, tn, m0_eps(EPS_B, H_FIN))
        check(f"{cid}-V2", f"{cid} with M_0 = B = 0's plateau / B = 0's r_th at eps = 1e-2 ({m0_eps(EPS_B, H_FIN):.6f}) instead of 0.223086: M - predicted within 6e-5", float(np.max(np.abs(v2["res"]))) <= bar_t(),
              f"(a, b, c) = {', '.join(f'{x:+.4f}' for x in v2['co'])}; " + ", ".join(f"{n} {x:+.2e}" for n, x in zip(tn, v2["res"])), ph=True)
        v3 = fit_pred(sub_ref(P2, P2[REF]), fn, tn, m0_const(H_FIN))
        check(f"{cid}-V3", f"{cid} with the symmetric member's r_th - 1 ({P2[REF]['rth'] - 1:+.1e}) and rdot_th ({P2[REF]['rdot']:+.1e}) at eps = 1e-2 subtracted from every member's (M_0 = 0.223086): M - predicted within 6e-5",
              float(np.max(np.abs(v3["res"]))) <= bar_t(), f"(a, b, c) = {', '.join(f'{x:+.4f}' for x in v3['co'])}; " + ", ".join(f"{n} {x:+.2e}" for n, x in zip(tn, v3["res"])), ph=True)
        v1 = fit_pred(P2, fn, tn, m0_const(H_FIN), "free")
        check(f"{cid}-V1", f"{cid} with a free intercept (four parameters; M_0 = 0.223086 + intercept): M - predicted within 6e-5", float(np.max(np.abs(v1["res"]))) <= bar_t(),
              f"(a, b, c, intercept) = {', '.join(f'{x:+.4f}' for x in v1['co'])}; " + ", ".join(f"{n} {x:+.2e}" for n, x in zip(tn, v1["res"])), ph=True)
        try:
            Ph = readall(tuple(dict.fromkeys(fn + tn)), EPS_B, H_ALT)
            vh = fit_pred(Ph, fn, tn, m0_const(H_ALT))
            check(f"{cid}-h", f"{cid} at h = {H_ALT:g} (M_0 = 0.223095, refit at h = 0.01 over the same points, eps = 1e-2): M - predicted within 6e-5", float(np.max(np.abs(vh["res"]))) <= bar_t(),
                  f"(a, b, c) = {', '.join(f'{x:+.4f}' for x in vh['co'])}; " + ", ".join(f"{n} {x:+.2e}" for n, x in zip(tn, vh["res"])), ph=True)
        except Exception as ex:
            check(f"{cid}-h", f"{cid} at h = {H_ALT:g}", False, f"not computable: {ex!r}", ph=True)

    if PART != "record":
        print("\n== Post hoc blocks that read the test points at eps != 1e-3 are skipped in the development partition")
        return finish(t0)

    try:
        post_hoc(P1, P2, allp, f1, t1, f2, t2)
    except Exception as ex:
        import traceback
        print(f"\n== Post hoc block not computable: {ex!r}")
        traceback.print_exc()
    return finish(t0)


MUTATIONS = {      # name: the checks it newly fails on the development partition (--part=dev --baseline=the unmutated log, BAR_DEV for the T checks); declared in the raw README before the test points were read at eps = 1e-2
    "asym_sign": ['E2B', 'E3F', 'E3N', 'E3P', 'E3S', 'PTD2'],
    "coef_a": ['E3N', 'E3S', 'PFIT2', 'T1', 'T1-V1', 'T1-V2', 'T1-V3', 'T1-h', 'T2', 'T2-V1', 'T2-V2', 'T2-V3', 'T2-h'],
    "coef_c": ['E3N', 'E3S', 'PFIT2'],
    "dep_eps": ['E2B', 'E3F', 'E3N', 'E3P', 'E3S', 'PDEP2', 'PTD2'],
    "drop_point": ['E3F', 'E3P'],
    "eps_swap": ['E2B', 'MEM'],
    "loo_off": ['E3P'],
    "m0_shift": ['CFG', 'E3F', 'E3N', 'E3P', 'E3S'],
    "no_rscale": ['E3N', 'E3S', 'PFIT2', 'T1', 'T1-V1', 'T1-V2', 'T1-V3', 'T1-h', 'T2', 'T2-V1', 'T2-V2', 'T2-V3', 'T2-h'],
    "no_run": ['RUNS'],
    "no_speed": ['E3F', 'E3P', 'PFIT2'],
    "rdot_nolapse": ['E3F', 'E3N', 'E3P', 'E3S', 'PTD2'],
}


def main():
    global MUT, PART
    argv = sys.argv[1:]
    cmd = argv[0] if argv else "help"
    opts = dict(a.lstrip("-").split("=", 1) if "=" in a else (a.lstrip("-"), True) for a in argv[1:])
    PART = opts.get("part", "record")
    if PART not in ("record", "dev"):
        raise SystemExit("--part=record|dev")
    if "mutate" in opts:
        MUT = opts["mutate"]
        if MUT not in MUTLIST:
            raise SystemExit(f"unknown mutation {MUT}; known: {', '.join(MUTLIST)}")
        if "baseline" in opts:
            BASE.update(ln.split("]", 1)[1].split(":", 1)[0].strip() for ln in open(opts["baseline"]) if ln.startswith("[FAIL]"))
        r.MUT = b.MUT = m.MUT = g.MUT = sm.MUT = r3.MUT = MUT
        tl.LMUT = MUT if MUT in ("asym_sign",) else None
    import platform
    import numba as nb
    print(f"settled_mass_eps.py {cmd} {opts}: numba {nb.__version__}, numpy {np.__version__}, python {platform.python_version()}, {platform.node()}, part {PART}, "
          f"mutation {MUT}, {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    if cmd == "analyze":
        sys.exit(analyze_cmd())
    raise SystemExit(__doc__)


MUTLIST = ("asym_sign", "coef_a", "coef_c", "dep_eps", "drop_point", "eps_swap", "loo_off", "m0_shift", "no_rscale", "no_run", "no_speed", "rdot_nolapse")


if __name__ == "__main__":
    main()

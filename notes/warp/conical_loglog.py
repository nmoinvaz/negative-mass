#!/usr/bin/env python3
"""The conical regime F ~ ln ln eta: its expansion about q = -p^2, its basin, and whether a member reaches a limit (issue #184)

N(F, F', F''; eta) = det_classify.conical_N() is, with t = ln eta, sigma = eta^-2, G = F, p = eta F', q = eta^2 F'' + eta F',
N = A0 + sigma A1 + sigma^2 A2 (conical_limit.md section 1). The exact solution p = 1/(t + T), q = -p^2 of A0 = 0 is the regime of
#167's regular crossings. This script expands about it (the root series), walks #167's forward grid to eta = 1e12 by the Rust crate's
precision-generic walker (conical_global/src/hp.rs, binary `loglog`) and decides whether a member reaches a limit.
  S1 theta_form, S2 root_series, S3 shift_series, S4 invariants    symbolic (sympy), an independent Laurent-polynomial verifier
  W0 regress, W1 basin, W2 approach, W3 no_limit, W4 invariants_conv, W5 n_floor, W6 dd_qd, W7 independent, W8 recip_check, W9 parity
  P1 plant_167, P2 plant_loglog, P3 plant_shift, P4 plant_limit, P5 plant_scale   planted cases through the same code
Usage (from the repository root, in the environment of the repository's setup: sympy 1.14.0, mpmath 1.3.0, cargo):
  python3 notes/warp/conical_loglog.py [--no-controls] [--mutate=NAME] [--only=a,b] [--logdir=DIR] [--nobuild]
Exit 1 if a check fails or a control does not fire. Status words as in AGENTS.md: everything numerical is numerical evidence."""
import hashlib
import json
import os
import subprocess
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
CRATE = os.path.join(HERE, "conical_global")
BIN_LL = os.path.join(CRATE, "target", "release", "loglog")
BIN_OLD = os.path.join(CRATE, "target", "release", "conical_global")
OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate")
ONLY = OPTS["only"].split(",") if "only" in OPTS else None

import mpmath as mp  # noqa: E402
import sympy as sp  # noqa: E402

_argv, sys.argv = sys.argv, [sys.argv[0]]  # det_classify and conical_limit parse sys.argv at import
sys.path.insert(0, HERE)
import det_classify as dc  # noqa: E402
sys.argv = _argv

mp.mp.dps = 80
RESULTS = []
LOGDIR = OPTS.get("logdir")


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def note(msg):
    print(f"  {msg}", flush=True)


def wanted(name):
    return ONLY is None or name in ONLY


# ---------------------------------------------------------------------------------------------------------------
# the equation in the theta form, typed from conical_limit.md section 1
# ---------------------------------------------------------------------------------------------------------------
G, p, q, s = sp.symbols("G p q sigma")
NU = sp.Symbol("nu")


def theta_parts():
    A0 = p**2 * (p**4 - q**2)
    A1 = G**2 * p**2 * q**2 + 4 * G * p**4 * q - 2 * G * q**3 + 4 * p**6 - 3 * p**2 * q**2
    A2 = 2 * (G * (q - p) + p**2) * (G**2 * q**2 + 2 * G * p**3 + 2 * G * p**2 * q + 2 * p**4 - q**2)
    if MUT == "theta-coeff":
        A1 = A1 + p**4 * q * G  # 4 G p^4 q -> 5 G p^4 q
    return A0, A1, A2


def n_theta():
    A0, A1, A2 = theta_parts()
    return A0 + s * A1 + s**2 * A2


# ---------------------------------------------------------------------------------------------------------------
# S1 theta_form
# ---------------------------------------------------------------------------------------------------------------
def n_exact(terms, pt):
    f0, f1, f2, e = pt
    v = sm = Fraction(0)
    for a, b, c, m, co in terms:
        t = co * f0**a * f1**b * f2**c * e**m
        v += t
        sm += abs(t)
    return v, sm


def ntable_points():
    return [(Fraction((7 * i + 3) % 23 - 11, 7), Fraction((11 * i + 5) % 29 - 14, 9), Fraction((13 * i + 7) % 31 - 15, 11),
             Fraction((17 * i + 2) % 37 - 18, 13)) for i in range(40)]


def fr2(x):
    return Fraction(x[0]) + Fraction(x[1])


def c_theta_form():
    N = dc.conical_N()
    sub = N.subs({dc.F0_: G, dc.F1_: p * NU, dc.F2_: (q - p) * NU**2, dc.ETA: 1 / NU})
    A0, A1, A2 = theta_parts()
    diff = sp.expand(sub - (A0 + NU**2 * A1 + NU**4 * A2))
    ok1 = diff == 0
    # the crate's table (#167's ntable, rerun with the old binary) against the typed conical_N
    Pn = sp.Poly(N, dc.F0_, dc.F1_, dc.F2_, dc.ETA)
    terms = [(a, b, c, m, Fraction(int(co.p), int(co.q))) for (a, b, c, m), co in zip(Pn.monoms(), Pn.coeffs())]
    pts = ntable_points()
    path = os.path.join(CRATE, "target", "ntable_points_i184.txt")
    with open(path, "w") as fh:
        for pt in pts:
            fh.write(" ".join(f"{x.numerator}/{x.denominator}" for x in pt) + "\n")
    out = old_crate(["ntable", path])
    worst = Fraction(0)
    for pt, o in zip(pts, out):
        v, sm = n_exact(terms, pt)
        worst = max(worst, abs(fr2(o["n"]) - v) / sm if sm else abs(fr2(o["n"]) - v))
    ok2 = worst <= Fraction(1, 10**26)
    check("theta_form", ok1 and ok2, f"conical_N with F' = p nu, F'' = (q - p) nu^2, eta = 1/nu minus A0 + nu^2 A1 + nu^4 A2 expands to {diff}; "
          f"the crate's N at 40 rational points against conical_N: worst |dN|/sum|terms| {float(worst):.2e} (threshold 1e-26)")


# ---------------------------------------------------------------------------------------------------------------
# S2 root series, S3 its generator on a planted shift, S4 the invariants
# ---------------------------------------------------------------------------------------------------------------
def root_series(Nexpr, K, tag=""):
    """e_1 .. e_K with N(G, p, -p^2 + sum sigma^k e_k, sigma) = O(sigma^(K+1)); N a polynomial in q with N(q = -p^2, sigma = 0) = 0"""
    w = sp.Symbol("w")
    Nw = sp.expand(Nexpr.subs(q, w - p**2))
    es, wser = [], 0
    t0 = time.time()
    for k in range(1, K + 1):
        ek = sp.Symbol("ek")
        r = sp.expand(Nw.subs(w, wser + s**k * ek))
        sol = sp.factor(sp.solve(r.coeff(s, k), ek)[0])
        es.append(sol)
        wser = wser + s**k * sol
    note(f"{tag}root series to order {K}: {time.time() - t0:.0f} s")
    return es


class LP:
    """Laurent polynomial in G, p and polynomial in sigma truncated above order `top`: {(a, b, c): Fraction} for G^a p^b sigma^c"""
    top = 0

    def __init__(self, d=None):
        self.d = {k: v for k, v in (d or {}).items() if v != 0}

    @staticmethod
    def of(expr):
        expr = sp.expand(expr)
        if expr == 0:
            return LP()
        M = 12
        pl = sp.Poly(sp.expand(expr * p**M), G, p, s)
        return LP({(a, b - M, c): Fraction(int(co.p), int(co.q)) for (a, b, c), co in zip(pl.monoms(), pl.coeffs()) if c <= LP.top})

    def __add__(self, o):
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d.get(k, 0) + v
        return LP(d)

    def __mul__(self, o):
        d = {}
        for (a, b, c), v in self.d.items():
            for (a2, b2, c2), v2 in o.d.items():
                if c + c2 <= LP.top:
                    k = (a + a2, b + b2, c + c2)
                    d[k] = d.get(k, 0) + v * v2
        return LP(d)

    def order(self):
        return min((c for (_, _, c) in self.d), default=None)


def residual_order(Nexpr, es, K, bump=False):
    """the lowest order in sigma of N(G, p, Q_K, sigma), by the independent Laurent arithmetic (None: zero through order K + 1)"""
    LP.top = K + 1
    Q = LP.of(-p**2)
    for k, e in enumerate(es, 1):
        ek = LP.of(e * s**k)
        if bump and k == 2:
            ek = ek + LP({(2, 0, 2): Fraction(1, 7)})
        Q = Q + ek
    P = sp.Poly(sp.expand(Nexpr), q)
    deg = P.degree()
    coefs = [LP.of(c) for c in reversed(P.all_coeffs())]
    acc = LP()
    pw = LP.of(sp.Integer(1))
    for j in range(deg + 1):
        acc = acc + coefs[j] * pw
        pw = pw * Q
    return acc.order()


def c_root_series():
    K = 6
    N = n_theta()
    es = root_series(N, K)
    ro = residual_order(N, es, K, bump=(MUT == "root-coef"))
    e1_ok = sp.simplify(es[0] + p**2 * (G - 1) ** 2 / 2) == 0
    a1 = sp.factor(theta_parts()[1].subs(q, -(p**2)))
    # printed structure: the p-range and the G-degree of each e_k
    rows = []
    for k, e in enumerate(es, 1):
        pl = sp.Poly(sp.expand(e * p**12), G, p)
        pw = [m[1] - 12 for m in pl.monoms()]
        gd = max(m[0] for m in pl.monoms())
        rows.append(f"e{k}: p^{min(pw)}..p^{max(pw)}, deg_G {gd}")
    check("root_series", (ro is None or ro >= K + 1) and e1_ok and a1 == p**6 * (G - 1) ** 2,
          f"N(G, p, Q_{K}, sigma) has no term below order sigma^{K + 1} (Laurent verifier: lowest order {ro}); A1(G, p, -p^2) = {a1}; "
          f"e1 = {sp.factor(es[0])}; e2 = {sp.factor(es[1])}; " + "; ".join(rows))
    ROOT["es"] = es
    return es


ROOT = {}


def get_root():
    if "es" not in ROOT:
        N = n_theta()
        ROOT["es"] = root_series(N, 6, tag="(for the checks) ")
    return ROOT["es"]


def shift_poly(c):
    w1 = q + p**2 + 2 * c * s * (1 - p) + c**2 * s**2
    return p**2 * w1 * (2 * p**2 - w1)


def c_shift_series():
    c = sp.Rational(1, 3)
    N = shift_poly(c)
    es = root_series(N, 4, tag="(planted) ")
    want = [-2 * c * (1 - p), -(c**2), 0, 0]
    if MUT == "shift-series":
        want[1] = -(c**2) / 2
    ok = all(sp.simplify(e - w) == 0 for e, w in zip(es, want))
    ro = residual_order(N, es, 4)
    check("shift_series", ok and (ro is None or ro >= 5), f"the generator on N = p^2 w'(2p^2 - w'), w' = q + p^2 + 2c sigma (1 - p) + c^2 sigma^2, c = 1/3: "
          f"e1 = {sp.factor(es[0])}, e2 = {sp.factor(es[1])}, e3 = {es[2]}, e4 = {es[3]}; expected -2c(1 - p), -c^2, 0, 0; residual order {ro}")


def c_invariants():
    t = sp.Symbol("t")
    Gf = sp.Function("Gf")(t)
    pf = sp.Function("pf")(t)
    qf = sp.Function("qf")(t)
    T = 1 / pf - t
    g = Gf + sp.log(pf)
    dT = sp.diff(T, t).subs(sp.Derivative(pf, t), qf)
    dg = sp.diff(g, t).subs(sp.Derivative(pf, t), qf).subs(sp.Derivative(Gf, t), pf)
    sgn = -1 if MUT == "inv-sign" else 1
    ok1 = sp.simplify(dT - sgn * (-(qf + pf**2) / pf**2)) == 0
    ok2 = sp.simplify(dg - (qf + pf**2) / pf) == 0
    ok3 = sp.simplify(dT.subs(qf, -(pf**2))) == 0 and sp.simplify(dg.subs(qf, -(pf**2))) == 0
    Tc, gc = sp.symbols("T g")
    pe = 1 / (t + Tc)
    Ge = gc + sp.log(t + Tc)
    ok4 = sp.simplify(sp.diff(Ge, t) - pe) == 0 and sp.simplify(sp.diff(pe, t) + pe**2) == 0
    check("invariants", ok1 and ok2 and ok3 and ok4, f"T = 1/p - t: T' = -(q + p^2)/p^2 ({ok1}); g = G + ln p: g' = (q + p^2)/p ({ok2}); both vanish at q = -p^2 ({ok3}); "
          f"p = 1/(t + T), G = ln(t + T) + g solves G' = p, p' = -p^2 ({ok4})")


# ---------------------------------------------------------------------------------------------------------------
# the crates
# ---------------------------------------------------------------------------------------------------------------
def old_crate(args):
    env = dict(os.environ, RAYON_NUM_THREADS="2")
    r = subprocess.run([BIN_OLD] + [str(a) for a in args], capture_output=True, text=True, env=env, timeout=3600)
    if r.returncode != 0:
        raise RuntimeError(f"conical_global {args}: {r.stderr[-400:]}")
    return [json.loads(ln) for ln in r.stdout.splitlines() if ln.startswith(("{", "["))]


def build():
    if "nobuild" in OPTS:
        return
    env = dict(os.environ, RAYON_NUM_THREADS="2")
    r = subprocess.run(["cargo", "build", "--release", "-j", "2"], cwd=CRATE, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        print(r.stderr[-3000:])
        sys.exit(2)



# ---------------------------------------------------------------------------------------------------------------
# the loglog binary: runs, states as sums of parts, and the verdict of a walk
# ---------------------------------------------------------------------------------------------------------------
CACHE = {}
BIN_ID = {}
CRATE_MUT = {"walk-eq": "eq-coeff", "recip-dd": "recip-dd", "no-project": "no-project", "plant-c": "plant-c"}


def bin_id():
    if "id" not in BIN_ID:
        BIN_ID["id"] = hashlib.sha256(open(BIN_LL, "rb").read()).hexdigest()[:12]
    return BIN_ID["id"]


def ll(args, mutate=None):
    """one run of `loglog run ...`: the list of JSON records (checkpoints, then the end record); cached by arguments"""
    args = [str(a) for a in args]
    if mutate and mutate in CRATE_MUT:
        args.append(f"--mutate={CRATE_MUT[mutate]}")
    key = bin_id() + " " + " ".join(args)
    if key in CACHE:
        return CACHE[key]
    env = dict(os.environ, RAYON_NUM_THREADS="2")
    r = subprocess.run([BIN_LL, "run"] + args, capture_output=True, text=True, env=env, timeout=3600)
    if r.returncode != 0:
        raise RuntimeError(f"loglog {args}: exit {r.returncode}: {r.stderr[-500:]}")
    recs = [json.loads(ln) for ln in r.stdout.splitlines() if ln.startswith("{")]
    CACHE[key] = recs
    return recs


def cmut():
    return MUT if MUT in CRATE_MUT else None


def parts_of(x, n):
    out, r = [], mp.mpf(x)
    for _ in range(n):
        f = float(r)
        out.append(f)
        r -= mp.mpf(f)
    return out


def pstr(x, n):
    return ",".join(repr(v) for v in parts_of(x, n))


def val(pp):
    return sum(mp.mpf(v) for v in pp[::-1])


NP = {"dd": 2, "qd": 4, "od": 8}


def state_args(exf, x0, prec):
    n = NP[prec]
    return ";".join(pstr(v, n) for v in exf(x0))


def gpq(rec):
    """(eta, G, p, q) of a record, in mpmath"""
    x = mp.mpf(rec.get("cp", rec.get("x")))
    F, F1, F2 = (val(rec[k]) for k in ("F", "F1", "F2"))
    pp = x * F1
    return x, F, pp, x * x * F2 + pp


def records_at(recs):
    return [r for r in recs if "cp" in r or r.get("end")]


def verdict(Gs, ps, qs):
    """the verdict of a walk to 1e12 from G, p, q at 1e6 and 1e12, already mirrored to p > 0:
    LIMIT if |G(1e12) - G(1e6)| <= 1e-3 and |p(1e12)| <= 1e-6; LOGLOG if |q + p^2|/p^2 <= 1e-10 at 1e12, p > 0 and not LIMIT; else OTHER"""
    (G6, G12), (p6, p12), (q12,) = Gs, ps, qs
    small_p = abs(p12) <= 1e-6 if MUT != "limit-class" else abs(p12) >= 1e-6
    if abs(G12 - G6) <= 1e-3 and small_p:
        return "LIMIT"
    sgn = -1 if MUT == "class-sign" else 1
    rho0 = abs((sgn * q12 + p12 * p12) / (p12 * p12))
    return "LOGLOG" if rho0 <= 1e-10 and p12 > 0 else "OTHER"


# ---------------------------------------------------------------------------------------------------------------
# planted cases: exact solutions
# ---------------------------------------------------------------------------------------------------------------
R_ = mp.mpf
T_PL, G_PL = R_(13) / 10, R_(7) / 10


def ex_pl2(eta):
    eta = R_(eta)
    t = mp.log(eta)
    tau = t + T_PL
    pp = 1 / tau
    qq = -pp * pp
    g = G_PL + (R_(1) / 1000 if MUT == "plant-exact" else 0)
    return g + mp.log(tau), pp / eta, (qq - pp) / eta**2


def ex_pl3(eta):
    eta = R_(eta)
    tau = mp.log(eta) + T_PL
    c = R_(1) / 3
    sg = 1 / eta**2
    pp = 1 / tau + c * sg
    qq = -1 / tau**2 - 2 * c * sg
    return G_PL + mp.log(tau) - c * sg / 2, pp / eta, (qq - pp) / eta**2


def ex_pl4(eta):
    e = R_(eta)
    al, a, b = R_(7) / 10, R_(13) / 10, R_(2) / 5
    v = 1 + e * e
    u = a + b * e
    return al + u / v, b / v - 2 * e * u / v**2, -4 * e * b / v**2 - 2 * u / v**2 + 8 * e * e * u / v**3


def ex_pl5(eta):
    e = R_(eta)
    return R_(7) / 10 - 1 / (3 * e**3), 1 / e**4, -4 / e**5


def plant_run(name, exf, x0, to, cps, prec, mode):
    args = ["--poly", name, "--prec", prec, "--mode", mode, "--x0", x0, "--s0", state_args(exf, x0, prec), "--to", repr(float(to)),
            "--cp", ",".join(repr(float(c)) for c in cps)]
    recs = ll(args, mutate=cmut())
    errs = []
    for r in records_at(recs):
        x = mp.mpf(r.get("cp", r.get("x")))
        ex = exf(x)
        errs.append((float(x), [float(abs((val(r[k]) - ex[i]) / ex[i])) for i, k in enumerate(("F", "F1", "F2"))]))
    return recs, errs


CPS_BIG = [200, 1e3, 1e6, 1e9]


def fmt_errs(errs):
    return " ".join(f"{x:.0e}:{max(e):.1e}" for x, e in errs)


def c_plant_167():
    outs = {}
    for kind in ("linear", "cubic"):
        o = old_crate(["plant", kind])[0]
        ok = o.get("converged", False) and all(o[k] is not None and o[k] <= 1e-20 for k in ("err_cm", "err_e", "err_cp")) \
            and o["err_nodes"] is not None and o["err_nodes"] <= 1e-22
        outs[kind] = (ok, o)
    bl = old_crate(["plant", "blowup"])
    okb = len(bl) == 8 and all("err" in o and o["err"] <= 1e-8 and o["bar"] <= 1e-8 and o["kind"] == ("F->inf" if o["case"] == "cube" else "F->0") for o in bl)
    ok = all(v[0] for v in outs.values()) and okb
    check("plant_167", ok, "#167's planted cases through the rebuilt crate: linear (|C- - a|, |e|, |C+ - b|, nodes) "
          + ", ".join(f"{k}: {o['err_cm']:.1e}, {o['err_e']:.1e}, {o['err_cp']:.1e}, {o['err_nodes']:.1e}" for k, (_, o) in outs.items())
          + f" (thresholds 1e-20, 1e-22); blowup: {len(bl)} walks, worst error {max(o.get('err', 1) for o in bl):.1e}, worst bar {max(o.get('bar', 1) for o in bl):.1e} (1e-8)")


def c_plant_loglog():
    res, ok = [], True
    for prec, mode, bar in (("dd", "project", 1e-26), ("dd", "conserve", 1e-26), ("qd", "project", 1e-55), ("qd", "conserve", 1e-55)):
        recs, errs = plant_run("pl2", ex_pl2, "2", 1e12, CPS_BIG, prec, mode)
        worst = max(max(e) for _, e in errs)
        reached = recs[-1]["stop"] == "End"
        ok &= reached and worst <= bar
        res.append(f"{prec} {mode}: worst {worst:.1e} (bar {bar:.0e}), {recs[-1]['steps']} steps")
    recs = ll(["--poly", "pl2", "--prec", "dd", "--mode", "project", "--x0", "2", "--s0", state_args(ex_pl2, 2, "dd"), "--to", "1e12", "--cp", "1e6"], mutate=cmut())
    rr = records_at(recs)
    _, G6, p6, q6 = gpq(rr[0])
    _, G12, p12, q12 = gpq(rr[-1])
    v = verdict((G6, G12), (p6, p12), (q12,))
    ok &= v == "LOGLOG"
    check("plant_loglog", ok, f"q + p^2 = 0 from eta = 2 (T = 13/10, g = 7/10) to 1e12, error of F, F', F'' at 200, 1e3, 1e6, 1e9, 1e12: " + "; ".join(res) + f"; verdict {v}")


def c_plant_shift():
    c = R_(1) / 3
    res, ok = [], True
    cps = [200, 1e3, 1e6, 1e9]
    recs, errs = plant_run("pl3", ex_pl3, "2", 1e12, cps, "dd", "project")
    worst = max(max(e[1:]) for _, e in errs)
    ok &= recs[-1]["stop"] == "End" and worst <= 1e-25
    rr = {r.get("cp", r.get("x")): r for r in records_at(recs)}
    r6 = rr[1e6]
    x, _, pp, _ = gpq(r6)
    tau = mp.log(x) + T_PL
    ratio = (pp - 1 / tau) / (c / x**2)
    ok &= abs(ratio - 1) <= 1e-3
    check("plant_shift", ok, f"N = p^2 w'(2p^2 - w'), w' = q + p^2 + 2c sigma (1 - p) + c^2 sigma^2, c = 1/3, exact p = 1/(t + T) + c sigma: DD project to 1e12, worst error of F', F'' "
          f"{worst:.1e} (bar 1e-25) at 200, 1e3, 1e6, 1e9, 1e12; (p - 1/(t + T))/(c sigma) at 1e6 = {mp.nstr(ratio, 12)} (1 +- 1e-3)")


def c_plant_limit():
    recs, errs = plant_run("pl4", ex_pl4, "10", 1e12, CPS_BIG, "dd", "project")
    worstF = max(e[0] for _, e in errs)
    worstD = max(max(e[1:]) for _, e in errs)
    rr = records_at(recs)
    _, G6, p6, q6 = gpq(rr[-2])
    _, G12, p12, q12 = gpq(rr[-1])
    v = verdict((G6, G12), (p6, p12), (q12,))
    recs_c, errs_c = plant_run("pl4", ex_pl4, "10", 1e12, CPS_BIG, "dd", "conserve")
    worstD_c = max(max(e[1:]) for _, e in errs_c)
    check("plant_limit", worstF <= 1e-24 and v == "LIMIT" and recs[-1]["stop"] == "End",
          f"the planted linear equation (1 + eta^2)F'' + 4 eta F' + 2(F - alpha), alpha = 7/10, F = alpha + (a + b eta)/(1 + eta^2) to 1e12: DD project worst error of F {worstF:.1e} "
          f"(bar 1e-24), verdict {v}; printed: F', F'' errors {worstD:.1e} (project) against {worstD_c:.1e} (conserve): project slaves F'' to the cancellation in F - alpha")


def c_plant_scale():
    cps = [30, 100, 300, 1000]
    out, ok = {}, True
    bars = {("dd", "conserve"): ("fail", 1e-3), ("qd", "conserve"): ("pass", 1e-20), ("dd", "project"): ("pass", 1e-24), ("qd", "project"): ("pass", 1e-50)}
    for (prec, mode), (want, bar) in bars.items():
        recs, errs = plant_run("pl5", ex_pl5, "10", 3000, cps, prec, mode)
        ep = errs[-1][1][1]
        good = (ep >= bar) if want == "fail" else (ep <= bar)
        if MUT == "no-project" and (prec, mode) == ("dd", "project"):
            pass
        ok &= good and recs[-1]["stop"] == "End"
        out[(prec, mode)] = (ep, want, bar, good)
    check("plant_scale", ok, "N = p^4 (q + 3p), F' = eta^-4, p falls by 2.7e7 and the scale of N by 7e-38 from eta = 10 to 3000; relative error of F' at 3000: "
          + "; ".join(f"{k[0]} {k[1]}: {v[0]:.1e} (predicted {v[1]}, bar {v[2]:.0e}: {'as predicted' if v[3] else 'NOT as predicted'})" for k, v in out.items())
          + f"; u (p0/p)^5 predicts 1e5 (DD) and 1e-27 (QD)")


def c_recip_check():
    def final(extra, prec="qd"):
        recs = ll(["--poly", "conical", "--prec", prec, "--mode", "project", "--fwd", "2", "1", "--to", "1e6", "--cp", "1e3"] + extra, mutate=cmut())
        e = recs[-1]
        return e["steps"], gpq(e)

    exact = []
    for hs in ("1", "0.5"):
        exact.append(final(["--hscale", hs]))
    (n1, a), (n2, b) = exact
    d = max(abs(a[1] - b[1]), abs(a[2] - b[2]))
    ok = d <= 1e-55
    check("recip_check", ok, f"QD walk of (F(0), F'(0)) = (2, 1) to 1e6: halving every step ({n1} to {n2} steps) changes G and p by {float(d):.1e} (bar 1e-55); "
          f"the control (--mutate=recip-dd, constants rounded to double-double) is run in the controls and must exceed 1e-45")



# ---------------------------------------------------------------------------------------------------------------
# the forward grid of #167, the members, their walks
# ---------------------------------------------------------------------------------------------------------------
F0S = ["0.1", "0.25", "0.5", "0.75", "1.5", "2", "3", "5", "-0.5", "-2"]
F1S = ["-100", "-10", "-1", "-0.1", "-0.01", "0.01", "0.1", "1", "10", "100"]
CPS_REAL = [200, 1e3, 1e4, 1e6, 1e9]
GRID = {}


def old_forward():
    if "rows" not in GRID:
        rows = old_crate(["forward", ",".join(F0S), ",".join(F1S), "--rstop=1e-8"])
        GRID["rows"] = {(r["f"], r["f1"]): r for r in rows}
    return GRID["rows"]


def cell(sd):
    if sd["class"] == "SINGULAR":
        return f"{sd['kind']}@{sd['xstar']:.6g}"
    if sd["class"] == "NONFINITE":
        return f"kappa@{sd['xstar']:.6g}"
    return sd["class"]


def printed_table():
    """#167's printed forward table: {(f, f1): token}"""
    path = os.path.join(HERE, "conical_global_output.txt")
    tab = {}
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("   f = "):
            f, rest = ln[7:].split(":", 1)
            for tok in rest.split():
                f1, t = tok.split(":", 1)
                tab[(float(f), float(f1))] = t
    return tab


def members():
    rows = old_forward()
    return [k for k in sorted(rows) if rows[k]["plus"]["class"] == "LOGLOG"]


def walk_member(f, f1, prec, mode, to="1e12", cps=None, extra=()):
    cps = CPS_REAL if cps is None else cps
    args = ["--poly", "conical", "--prec", prec, "--mode", mode, "--fwd", repr(f), repr(f1), "--to", to,
            "--cp", ",".join(repr(float(c)) for c in cps)] + list(extra)
    return ll(args, mutate=cmut())


def mirrored(Gv, pv, qv):
    sg = 1 if pv > 0 else -1
    return sg, sg * Gv, sg * pv, sg * qv


def Tg(x, Gv, pv):
    return (1 / pv + mp.log(x) if MUT == "inv-conv" else 1 / pv - mp.log(x)), Gv + mp.log(pv)


EFUN = {}


def efuns(K=6):
    if "f" not in EFUN:
        es = get_root()
        EFUN["f"] = [sp.lambdify((G, p), e, "mpmath") for e in es[:K]]
    return EFUN["f"]


def Qk(K, Gv, pv, sig):
    """Q_K(G, p, sigma) = -p^2 + sum_{k <= K} sigma^k e_k(G, p), mirrored variables"""
    ef = efuns()
    v = -pv * pv
    for k in range(1, K + 1):
        ek = ef[k - 1](Gv, pv)
        if MUT == "root-coef" and k == 2:
            ek = ek * (1 + mp.mpf(10) ** -3)
        v += sig**k * ek
    return v


def rho(K, x, Gv, pv, qv):
    sg, Gm, pm, qm = mirrored(Gv, pv, qv)
    return abs(qm - Qk(K, Gm, pm, 1 / x**2)) / pm**2


def c_regress():
    rows = old_forward()
    bad, states, tiny = [], [], []
    for (f, f1), r in sorted(rows.items()):
        recs = ll(["--poly", "conical", "--prec", "dd", "--mode", "conserve", "--fwd", repr(f), repr(f1), "--to", "200", "--max-steps", "20000"], mutate=cmut())
        e = recs[-1]
        reached = e["stop"] == "End"
        is_ll = r["plus"]["class"] == "LOGLOG"
        if reached != is_ll:
            bad.append((f, f1, e["stop"], r["plus"]["class"]))
        if is_ll and reached:
            o = old_crate(["fwdstate", repr(f), repr(f1), "200"])[0]
            worst = Fraction(0)
            for k, name in enumerate(("F", "F1", "F2")):
                a = Fraction(o["state"][k][0]) + Fraction(o["state"][k][1])
                b = sum((Fraction(v) for v in e[name]), Fraction(0))
                worst = max(worst, abs(a - b) / abs(a))
            states.append(((f, f1), float(worst)))
        elif not is_ll and r["plus"]["xstar"] is not None:
            tiny.append(abs(e["x"] - r["plus"]["xstar"]) / max(r["plus"]["xstar"], 1e-12))
    wmax = max(w for _, w in states) if states else float("nan")
    reached_ok = not bad and len(states) == 15
    status = "met" if reached_ok and wmax <= 1e-24 else "NOT MET (as first worded)"
    print(f"  {status} [W0 regress, does not gate]: 100 plus sides to 200: the new DD walker reaches 200 exactly where #167's class is LOGLOG ({len(states)} of 100; "
          f"disagreements {bad[:4]}); F, F', F'' at 200 against #167's walker for those: worst relative difference {wmax:.1e} (bar 1e-24). Both are N-conserving double-double "
          f"walks: the rounding of N from the first steps stays in N while the scale of N shrinks, so each is good only to about its own n_rel at 200 (up to 3e-17 at F'(0) = 100); "
          f"check dd_qd (DD project against QD) is the accuracy gate. The others stop at TinyStep near #167's eta*: |x_end - eta*|/eta* at most {max(tiny):.1e} over {len(tiny)} (printed, no bar)")


RAW167 = OPTS.get("raw167", os.path.join(HERE, "..", "..", "..", "negative-mass-raw", "I167-mac", "dev_09_forward.log"))


def raw167_rows():
    """#167's raw records of `forward ... --rstop=1e-8`: {(f, f1): row}"""
    rows = {}
    for ln in open(RAW167, encoding="utf-8"):
        if ln.startswith("{"):
            r = json.loads(ln)
            rows[(r["f"], r["f1"])] = r
    return rows


def side_token(sd):
    if sd["class"] == "SINGULAR":
        return "S:" + sd["kind"]
    if sd["class"] == "NONFINITE":
        return "NONFINITE"
    return sd["class"]


def c_basin():
    rows = old_forward()
    ref = raw167_rows()
    tab = printed_table()
    bad, worst, nxs = [], 0.0, 0
    for key, r in rows.items():
        a, b = r["plus"], ref[key]["plus"]
        if cell(a).split("@")[0] != tab[key].split("@")[0]:
            bad.append((key, "label", cell(a), tab[key]))
        if a.get("xstar") is not None and b.get("xstar") is not None:
            nxs += 1
            d = abs(a["xstar"] - b["xstar"])
            worst = max(worst, d)
            if d > 1e-8:
                bad.append((key, "xstar", a["xstar"], b["xstar"], d))
        elif (a.get("xstar") is None) != (b.get("xstar") is None):
            bad.append((key, "eta* present in one only", a.get("xstar"), b.get("xstar")))
    mem = members()
    ref_mem = sorted(k for k, v in tab.items() if v == "LOGLOG")
    ok_set = mem == ref_mem and len(mem) == 15
    reach = []
    for f, f1 in mem:
        recs = walk_member(f, f1, "dd", "project")
        reach.append(recs[-1]["stop"] == "End" and recs[-1]["x"] == 1e12 and member_verdict(f, f1)[0] == "LOGLOG")
    kappa_bad = [b for b in bad if b[1] == "xstar" and rows[b[0]]["plus"]["class"] == "NONFINITE"]
    bad = [b for b in bad if b not in kappa_bad]
    check("basin", not bad and ok_set and all(reach), f"#167's forward grid (100 solutions, plus sides), a rerun of #167's crate here: the class and kind labels equal #167's printed table at {100 - sum(1 for b in bad if b[1] == 'label')} of 100; "
          f"eta* within 1e-8 of #167's raw records (I167-mac/dev_09_forward.log, full digits) at {nxs - sum(1 for b in bad if b[1] == 'xstar') - len(kappa_bad)} of {nxs} points that have one (the {len(kappa_bad)} kappa-locus points below are NOT MET as first worded and do not gate), worst |d eta*| over all {worst:.1e}; "
          f"other disagreements {[b for b in bad if b[1] not in ('label', 'xstar')][:3]}; the LOGLOG set is the printed 15: {ok_set}; "
          f"the DD project walks of the 15 reach 1e12 with verdict LOGLOG: {sum(reach)} of {len(reach)}")
    if bad:
        print("    first disagreements: " + "; ".join(str(b) for b in bad[:6]))
    if kappa_bad:
        print(f"  NOT MET (as first worded) [W1, does not gate]: eta* to 1e-8 at the kappa-locus points, where eta* is the x of an F'' overflow and #167 states it to 1e-6 only: "
              f"{len(kappa_bad)} points differ from the raw records by {min(b[4] for b in kappa_bad):.1e} to {max(b[4] for b in kappa_bad):.1e}: "
              + "; ".join(f"({b[0][0]:g}, {b[0][1]:g})" for b in kappa_bad))
    print("  the 15 members (F(0), F'(0)): " + " ".join(f"({f:g},{f1:g})" for f, f1 in mem))


def reaches(f, f1):
    """whether the DD conserve walk of (F(0), F'(0)) = (f, f1) reaches eta = 200 in the log-log band of #167 (|q/p^2 + sgn p| <= 0.2)"""
    recs = ll(["--poly", "conical", "--prec", "dd", "--mode", "conserve", "--fwd", repr(f), repr(f1), "--to", "200", "--max-steps", "20000"])
    e = recs[-1]
    if e["stop"] != "End":
        return False
    _, Gv, pv, qv = gpq(e)
    sg = 1 if pv > 0 else -1
    return abs(qv / pv**2 + sg) <= 0.2


def c_basin_map():
    nf, n1 = 20, 20
    fs = [float(mp.mpf(0.1) * (mp.mpf(50) ** (mp.mpf(i) / (nf - 1)))) for i in range(nf)]
    f1p = [float(mp.mpf("0.01") * (mp.mpf(10000) ** (mp.mpf(i) / (n1 - 1)))) for i in range(n1)]
    cols = [-v for v in reversed(f1p)] + f1p
    n_l = 0
    print(f"  basin map (not gating): rows F(0) from 5 down to 0.1, columns F'(0) = -100 ... -0.01, 0.01 ... 100 (log-spaced); L: the DD walk reaches 200 in the band q = -p^2 (p > 0) or q = +p^2 (p < 0); . : it stops at a singularity")
    for f in reversed(fs):
        row = "".join("L" if reaches(f, c) else "." for c in cols)
        n_l += row.count("L")
        print(f"    F(0) = {f:7.3f} {row}")
    print(f"    {n_l} of {nf * 2 * n1} walks reach 200")
    # the boundary in F'(0) for fixed F(0) > 1: the smallest F'(0) > 0 that reaches 200, by bisection (monotone along the row of the map)
    out = []
    for f in (1.05, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0):
        lo, hi = 0.01, 100.0
        if reaches(f, lo) or not reaches(f, hi):
            out.append(f"F(0) = {f:g}: not bracketed")
            continue
        for _ in range(45):
            mid = (lo * hi) ** 0.5
            if reaches(f, mid):
                hi = mid
            else:
                lo = mid
        out.append(f"F(0) = {f:g}: F'(0)* = {hi:.8g}")
    print("  basin boundary (smallest F'(0) > 0 that reaches 200, bisection between 0.01 and 100, each end checked): " + "; ".join(out))


def c_approach():
    mem = members()
    K = 1 if MUT == "series-k" else 4
    worst = {1e3: 0.0, 1e6: 0.0, 1e12: 0.0}
    lines = []
    for f, f1 in mem:
        recs = records_at(walk_member(f, f1, "dd", "project"))
        row = []
        for r in recs:
            x, Gv, pv, qv = gpq(r)
            if float(x) in worst:
                worst[float(x)] = max(worst[float(x)], float(rho(K, x, Gv, pv, qv)))
            row.append(f"{float(x):.0e}: rho0 {float(rho(0, x, Gv, pv, qv)):.1e} rho1 {float(rho(1, x, Gv, pv, qv)):.1e} rho2 {float(rho(2, x, Gv, pv, qv)):.1e} "
                       f"rho4 {float(rho(4, x, Gv, pv, qv)):.1e} rho6 {float(rho(6, x, Gv, pv, qv)):.1e}")
        lines.append(f"    ({f:g}, {f1:g}): " + " | ".join(row))
    ok = all(w <= 1e-10 for w in worst.values())
    check("approach", ok, f"15 members, DD project walks: worst over the members of rho_{K} = |q - Q_{K}(G, p, eta^-2)|/p^2 (mirrored to p > 0) at eta = 1e3: {worst[1e3]:.1e}, 1e6: {worst[1e6]:.1e}, "
          f"1e12: {worst[1e12]:.1e} (bar 1e-10 at each)")
    for ln in lines:
        print(ln)


def member_verdict(f, f1, prec="dd", mode="project"):
    recs = records_at(walk_member(f, f1, prec, mode))
    byx = {float(gpq(r)[0]): gpq(r) for r in recs}
    _, G6, p6, q6 = byx[1e6]
    _, G12, p12, q12 = byx[1e12]
    sg = 1 if p12 > 0 else -1
    return verdict((sg * G6, sg * G12), (sg * p6, sg * p12), (sg * q12,)), byx


def c_no_limit():
    mem = members()
    out, ok = [], True
    worst = mp.mpf(0)
    for f, f1 in mem:
        v, byx = member_verdict(f, f1)
        _, G6, p6, _ = byx[1e6]
        _, G12, p12, _ = byx[1e12]
        sg = 1 if p6 > 0 else -1
        t1, t2 = mp.log(mp.mpf(10) ** 6), mp.log(mp.mpf(10) ** 12)
        T1 = 1 / (sg * p6) - t1
        pred = mp.log((t2 + T1) / (t1 + T1))
        d = abs(sg * (G12 - G6) - pred)
        worst = max(worst, d)
        ok &= v == "LOGLOG"
        out.append(f"({f:g},{f1:g}) {v} dG {float(sg * (G12 - G6)):.4f}")
    ok &= worst <= 1e-8
    check("no_limit", ok, f"verdicts at 1e12 of the 15 (LIMIT if |dG| <= 1e-3 and |p| <= 1e-6): " + ", ".join(out[:6]) + " ...; "
          f"G(1e12) - G(1e6) against ln((t2 + T1)/(t1 + T1)) with T1 = 1/p - t at 1e6: worst difference {float(worst):.1e} (bar 1e-8)")


def c_invariants_conv():
    mem = members()
    worstT = worstg = mp.mpf(0)
    tinf = []
    for f, f1 in mem:
        recs = records_at(walk_member(f, f1, "qd", "project"))
        vals = {}
        for r in recs:
            x, Gv, pv, qv = gpq(r)
            sg, Gm, pm, qm = mirrored(Gv, pv, qv)
            vals[float(x)] = Tg(x, Gm, pm)
        T12, g12 = vals[1e12]
        for e in (1e6, 1e9):
            worstT = max(worstT, abs(vals[e][0] - T12))
            worstg = max(worstg, abs(vals[e][1] - g12))
        tinf.append((f, f1, T12, g12))
    ok = worstT <= 1e-8 and worstg <= 1e-8
    check("invariants_conv", ok, f"15 members, QD project: |T(eta) - T(1e12)| and |g(eta) - g(1e12)| at eta = 1e6, 1e9: worst {float(worstT):.1e} and {float(worstg):.1e} (bar 1e-8). "
          f"T_inf, g_inf printed below (T = 1/p - t, g = G + ln p, mirrored to p > 0)")
    for f, f1, T, g in tinf:
        print(f"    ({f:g}, {f1:g}): T_inf = {mp.nstr(T, 15)}, g_inf = {mp.nstr(g, 15)}")


def c_n_floor():
    mem = members()
    wd = wq = 0.0
    cd = 0.0
    for f, f1 in mem:
        wd = max(wd, walk_member(f, f1, "dd", "project")[-1]["n_rel_max"])
        wq = max(wq, walk_member(f, f1, "qd", "project")[-1]["n_rel_max"])
        cd = max(cd, walk_member(f, f1, "dd", "conserve")[-1]["n_rel_max"])
    ok = wd <= 1e-26 and wq <= 1e-54
    check("n_floor", ok, f"max |N|/sum|terms| at the step ends over the 15 walks: DD project {wd:.1e} (bar 1e-26), QD project {wq:.1e} (bar 1e-54); printed: DD conserve {cd:.1e}")


def c_dd_qd():
    mem = members()
    worst = mp.mpf(0)
    for f, f1 in mem:
        a = records_at(walk_member(f, f1, "dd", "project"))
        b = records_at(walk_member(f, f1, "qd", "project"))
        for ra, rb in zip(a, b):
            _, G1, p1, q1 = gpq(ra)
            _, G2, p2, q2 = gpq(rb)
            worst = max(worst, abs(G1 - G2) / abs(G2), abs(p1 - p2) / abs(p2), abs(q1 - q2) / abs(q2))
    check("dd_qd", worst <= 1e-26, f"DD against QD (project) at 200, 1e3, 1e4, 1e6, 1e9, 1e12 for the 15: worst relative difference of G, p, q {float(worst):.1e} (bar 1e-26)")


def c_independent():
    picks = [(2.0, 1.0), (5.0, 100.0), (-2.0, -10.0)]
    worst = mp.mpf(0)
    detail = []
    kq = 0 if MUT == "indep-q" else 6
    for f, f1 in picks:
        recs = records_at(walk_member(f, f1, "qd", "project", cps=[1e3, 1e6, 1e9]))
        byx = {float(gpq(r)[0]): gpq(r) for r in recs}
        sg, G0, p0, _ = mirrored(*byx[1e3][1:])
        t0 = time.time()
        with mp.workdps(30):
            def rhs(t, y):
                sig = mp.e ** (-2 * t)
                return [y[1], Qk(kq, y[0], y[1], sig)]
            sol = mp.odefun(rhs, mp.log(mp.mpf(1000)), [mp.mpf(G0), mp.mpf(p0)], tol=mp.mpf(10) ** -22)
            for e in (1e6, 1e9, 1e12):
                _, Gw, pw, _ = byx[e]
                Gi, pi = sol(mp.log(mp.mpf(e)))
                d = max(abs(Gi - sg * Gw) / abs(Gw), abs(pi - sg * pw) / abs(pw))
                worst = max(worst, d)
                detail.append(f"({f:g},{f1:g}) at {e:.0e}: {float(d):.1e}")
        detail.append(f"[{time.time() - t0:.0f} s]")
    check("independent", worst <= 1e-12, f"own mpmath integration (30 digits, odefun) of G' = p, p' = Q_{kq}(G, p, e^(-2t)) from eta = 1e3 (QD state) against the QD walk: worst relative difference of G, p "
          f"{float(worst):.1e} (bar 1e-12): " + ", ".join(detail))


def c_parity():
    rows = old_forward()
    bad_a, worst_a, n_a = [], 0.0, 0
    for (f, f1), r in rows.items():
        a, b = r["minus"], rows[(f, -f1)]["plus"]
        if side_token(a) != side_token(b):
            bad_a.append(((f, f1), side_token(a), side_token(b)))
        elif a.get("xstar") is not None and b.get("xstar") is not None:
            n_a += 1
            d = abs(-a["xstar"] - b["xstar"])
            worst_a = max(worst_a, d)
            if d > 1e-8:
                bad_a.append(((f, f1), "xstar", a["xstar"], b["xstar"], d))
    worst_b, n_b = 0.0, 0
    sgn = -1 if MUT == "parity-sign" else 1
    for f, f1 in members():
        a = ll(["--poly", "conical", "--prec", "dd", "--mode", "conserve", "--fwd", repr(f), repr(-f1), "--to", "-200", "--max-steps", "20000"], mutate=cmut())[-1]
        b = ll(["--poly", "conical", "--prec", "dd", "--mode", "conserve", "--fwd", repr(f), repr(f1), "--to", "200", "--max-steps", "20000"], mutate=cmut())[-1]
        if a["stop"] != "End" or b["stop"] != "End":
            bad_a.append(((f, f1), "regular side not reached", a["stop"], b["stop"]))
            continue
        Fa, F1a, F2a = (val(a[k]) for k in ("F", "F1", "F2"))
        Fb, F1b, F2b = (val(b[k]) for k in ("F", "F1", "F2"))
        worst_b = max(worst_b, float(abs(Fa - Fb) / abs(Fb)), float(abs(F1a + sgn * F1b) / abs(F1b)), float(abs(F2a - F2b) / abs(F2b)))
        n_b += 1
    check("parity", not bad_a and worst_b <= 1e-24, f"W9': singular sides, the minus side of (f, f1) against the plus side of (f, -f1) in #167's crate rerun here: class, kind and eta* within 1e-8 at "
          f"{100 - len(bad_a)} of 100 ({n_a} with an eta*, worst {worst_a:.1e}; {bad_a[:3]}); regular sides: the walk of (f, -f1) to -200 against that of (f, f1) to +200 for the {n_b} members, "
          f"F, F' (sign) and F'' worst relative difference {worst_b:.1e} (bar 1e-24)")
    print("  NOT MET (as first worded): W9 asks for the walk of (f, f1) to -200 against the walk of (f, -f1) to +200 for three members whose minus sides end singular (eta* 0.5276, 0.1173, 0.0226), "
          "so neither walk reaches 200; W9' above replaces it (README A1, item 3)")



CHECKS = [
    ("theta_form", c_theta_form),
    ("root_series", c_root_series),
    ("shift_series", c_shift_series),
    ("invariants", c_invariants),
]
MUTATIONS = {
    "theta-coeff": (["theta_form", "root_series"], "one coefficient of A1 typed wrong"),
    "root-coef": (["root_series"], "one coefficient of e2 perturbed in the residual"),
    "shift-series": (["shift_series"], "the planted generator expected e2 = -c^2/2"),
    "inv-sign": (["invariants"], "T' expected +w/p^2"),
}


CHECKS += [
    ("regress", c_regress),
    ("basin", c_basin),
    ("basin_map", c_basin_map),
    ("approach", c_approach),
    ("no_limit", c_no_limit),
    ("invariants_conv", c_invariants_conv),
    ("n_floor", c_n_floor),
    ("dd_qd", c_dd_qd),
    ("independent", c_independent),
    ("parity", c_parity),
    ("plant_167", c_plant_167),
    ("plant_loglog", c_plant_loglog),
    ("plant_shift", c_plant_shift),
    ("plant_limit", c_plant_limit),
    ("plant_scale", c_plant_scale),
    ("recip_check", c_recip_check),
]
MUTATIONS.update({
    "class-sign": (["basin", "no_limit"], "the LOGLOG test uses q = +p^2 for p > 0"),
    "inv-conv": (["invariants_conv"], "T computed as 1/p + t"),
    "parity-sign": (["parity"], "F' summed with the wrong sign in the regular-side comparison"),
    "series-k": (["approach"], "the approach uses K = 1"),
    "walk-eq": (["approach", "independent"], "one coefficient of N in the hp walker changed"),
    "indep-q": (["independent"], "the independent integration uses Q_0 = -p^2 only"),
    "plant-exact": (["plant_loglog"], "the exact solution of the log-log plant has g + 0.001"),
    "plant-c": (["plant_shift"], "the planted equation has c = 1003/3000 in place of 1/3"),
    "limit-class": (["plant_limit", "plant_loglog"], "the LIMIT test reads |p| >= 1e-6"),
    "no-project": (["plant_scale", "n_floor", "dd_qd"], "project runs as conserve"),
    "recip-dd": (["recip_check"], "the QD lanes use double-double constants"),
})


def main():
    t0 = time.time()
    if MUT is not None and MUT not in MUTATIONS:
        print(f"unknown mutation {MUT}")
        sys.exit(2)
    build()
    for name, fn in CHECKS:
        if not wanted(name):
            continue
        t1 = time.time()
        fn()
        print(f"  ({time.time() - t1:.0f} s)", flush=True)
    failed = [n for n, ok in RESULTS if not ok]
    print(f"\nchecks: {len(RESULTS) - len(failed)} passed, {len(failed)} failed {failed}. {time.time() - t0:.0f} s")
    if MUT is not None:
        want = MUTATIONS[MUT][0]
        print(f"control {MUT}: target checks {want}, fired {[n for n in want if n in failed]}")
    return failed


if __name__ == "__main__":
    if "list" in OPTS:
        for n, (w, what) in MUTATIONS.items():
            print(n, ",".join(w))
        sys.exit(0)
    os.nice(10)
    failed = main()
    controls_ok = True
    if MUT is None and "no-controls" not in OPTS and ONLY is None:
        print("\n== Controls: each mutation must make its target checks fail, and no others ==", flush=True)
        for name, (want, what) in MUTATIONS.items():
            t1 = time.time()
            args = [sys.executable, os.path.abspath(__file__), f"--mutate={name}", "--nobuild", "--only=" + ",".join(want)]
            p_ = subprocess.run(args, capture_output=True, text=True)
            if LOGDIR:
                with open(os.path.join(LOGDIR, f"ctl_{name}.log"), "w") as fh:
                    fh.write(p_.stdout + p_.stderr)
                with open(os.path.join(LOGDIR, f"ctl_{name}.rc"), "w") as fh:
                    fh.write(f"{p_.returncode}\n")
            fails = [ln.split(":")[0][7:] for ln in p_.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = all(w in fails for w in want) and p_.returncode == 1
            extra = [f for f in fails if f not in want]
            controls_ok &= fired and not extra
            print(f"[{'FIRED' if fired and not extra else 'SILENT' if not fired else 'EXTRA'}] --mutate={name} ({what}): exit {p_.returncode}, failed {fails}, {time.time() - t1:.0f} s", flush=True)
            if p_.returncode not in (0, 1):
                print(p_.stderr[-1500:])
        print(f"controls: {'all fired on target only' if controls_ok else 'SOME DID NOT FIRE ON TARGET ONLY'}")
    sys.exit(0 if not failed and controls_ok else 1)

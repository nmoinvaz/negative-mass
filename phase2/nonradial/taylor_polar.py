#!/usr/bin/env python
"""#186: the polar matching determinant D(gamma) by a Rust Taylor integrator (crate taylor_polar/).

    python taylor_polar.py                    build if needed, run every check, exit 1 if any fails
    python taylor_polar.py --only=T1,P1       run the named checks only
    python taylor_polar.py --mutate=NAME      run with one deliberate error (a control: the check it targets must fail)
    python taylor_polar.py --gen              write taylor_polar/src/jgen.rs from polar_matrix_generated.py
    python taylor_polar.py --build            cargo build --release (nice 10, -j 2); --features=lib_table adds the cargo feature (#192:
                                              J's table from multishoot's table module; lib_plane likewise; both replace target/release's binary)
    python taylor_polar.py --scipy-timing OUT run the scipy path (convergence.py's base variant) at the timing points,
                                              appending one JSON line per evaluation to OUT (long: about 30 min)
    python taylor_polar.py --list             the checks and the controls
    python taylor_polar.py --complex ...      #188: the crate over complex s (see the block "the crate over complex s" below; --complex --list)
    python taylor_polar.py --rows [--l=2,4,5,6] [--only=QP,RA2,CAP5] [--jsonl-dir=DIR] [--mutate=NAME] [--raw=DIR]
                                              #187: convergence.py's rows (every variant, q_int, q_ext) from the crate, printed in
                                              the stored format and checked against convergence_l{2,4,5,6}_output.txt, with the
                                              two capped rows of #114; --print-rows=3 prints rows for an l without a stored file

jgen.rs holds J(r; s) of polar_matrix_generated.py (written by reduce_polar.py, the symbolic source) as generic Rust
over truncated power series: every statement of the Python function is translated, operation by operation, by walking
its syntax tree. The symbolic source is shared with polar_solver.py, not independent of it. The checks below name what
they share with the scipy path and what they do not.
"""
import ast
import json
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
CRATE = os.path.join(HERE, "taylor_polar")
BIN = os.path.join(CRATE, "target", "release", "taylor_polar")
SRC = os.path.join(HERE, "polar_matrix_generated.py")

OMEGA, EPS = 0.7603, -1.0
LAM, MU = 10240.0, 10240.0 ** 2 / 4.8   # as boson_star.py: LAM**2 / 4.8 in float64
MODEL = dict(omega=OMEGA, eps=EPS, lam=LAM, mu=repr(MU))
R_M = {"base": 24.10980364018201, "phi_cut": 34.20514015700785, "r0": 24.111}   # PolarProblem(...).rm of convergence.py
MEAN_SPREAD = 2.1e-4                                                          # convergence_l2_output.txt

# ---------------------------------------------------------------------------------------------------------------
# the declared criteria (negative-mass-raw I186-mac README.md has the same numbers and the reasons)
# ---------------------------------------------------------------------------------------------------------------
CRIT = dict(
    # planted systems through the same integrator, against mpmath at 90 digits (shares no code with the integrator)
    P_GAMMAS=[3.0, 10.0, 300.0, 12000.0],
    P_L=24.0, P_LZ=40.0,
    P_TOL={"f64": dict(D=1e-6, P=1e-9, LN=1e-9), "dd": dict(D=1e-24, P=1e-24, LN=1e-24),
           "qd": dict(D=1e-56, P=1e-56, LN=1e-56)},
    P_TRANSIENT_GL=40.0,
    # translation of J and its series
    T1_TOL=1e-24, T2_TOL=1e-24,
    # the star
    S1_TOL=1e-6, S1_TOL_MATCHED=1e-9, S2_TOL=1e-9, S3_TOL=1e-16, S3_RATIO=(1e-4 / 3.0, 3e-4), T3_TOL=1e-20, T3_RADII=[1e-3 * 1.037 ** i for i in range(150)], S_RADII=[0.002, 0.05, 0.5, 2.0, 5.0, 8.0, 12.0, 16.0, 20.0, 24.10980364018201],
    # reciprocals at lane precision
    R1_TOL={"qd": 8 * 1.505e-64, "od": 8 * 2.3e-128},
    # D against the stored rows
    A_BAR=MEAN_SPREAD,
    L2_GAMMAS=[0.001, 0.03, 0.3, 1.0, 3.0, 10.0, 40.0, 137.0],
    L6_GAMMAS=[0.001, 0.03, 0.3, 1.0, 3.0, 10.0, 40.0, 137.0, 300.0, 1000.0, 12000.0],
    L4_GAMMAS=[0.001, 0.03, 0.3, 1.0, 3.0, 10.0, 40.0, 137.0, 300.0, 1000.0, 3000.0],
    L5_GAMMAS=[0.001, 0.03, 0.3, 1.0, 3.0, 10.0, 40.0, 137.0, 300.0, 1000.0, 6000.0],
    # lanes
    LANE_GAMMAS=[(2, 300.0), (6, 1000.0), (6, 12000.0)],
    LANE_QD=[(2, 300.0), (6, 137.0)], L1_QD_TOL=1e-20,
    # variants of convergence.py: r0 = 3e-3 (spread bar), r_m ~ 34 (the stored rows have 4 digits), order and tolerance
    V1_POINTS=[(2, 1.0), (2, 137.0), (6, 1000.0)],
    V2_POINTS=[(2, 0.001), (2, 1.0), (2, 137.0), (6, 1.0), (6, 137.0), (6, 1000.0), (6, 12000.0)], V2_DIGITS_ALLOWANCE=0.5e-3,
    V3_POINTS=[(2, 10.0), (6, 1000.0)], V3_ORDER=24, V3_TOL=1e-22, V3_BAR=1e-15,
    # the scipy path on this Mac, base variant (taylor_polar.py --scipy-timing): points and the gate on its D
    V4_POINTS=[(2, 1.0), (2, 10.0), (6, 10.0), (2, 137.0)], V4_TOL=1e-7,
    TS_POINTS=[(2, 1.0), (2, 10.0), (2, 40.0), (2, 137.0), (2, 300.0), (6, 0.001), (6, 1.0), (6, 40.0), (6, 137.0), (6, 300.0), (6, 1000.0)],
)


# ---------------------------------------------------------------------------------------------------------------
# generation of jgen.rs
# ---------------------------------------------------------------------------------------------------------------
SCALARS = {"omega": "omega", "epsilon": "epsilon", "s": "s", "LL": "ll", "Lambda": "lambda", "mu": "mu", "pi": "pi"}
SERIES = {"r", "alpha", "a", "phi0", "u"}


class Gen:
    def __init__(self, fn):
        self.defs, self.order, self.entries, self.types = {}, [], [], {}
        for st in fn.body:
            if isinstance(st, ast.Assign) and len(st.targets) == 1:
                t = st.targets[0]
                if isinstance(t, ast.Name) and t.id.startswith("x") and t.id[1:].isdigit():
                    self.defs[t.id] = st.value
                    self.order.append(t.id)
                elif isinstance(t, ast.Subscript) and isinstance(t.value, ast.Name) and t.value.id == "J":
                    i, j = ast.literal_eval(t.slice)
                    self.entries.append((i, j, st.value))

    def needed(self):
        need, stack = set(), []
        for _, _, e in self.entries:
            for n in ast.walk(e):
                if isinstance(n, ast.Name) and n.id in self.defs and n.id not in need:
                    need.add(n.id)
                    stack.append(n.id)
        while stack:
            for n in ast.walk(self.defs[stack.pop()]):
                if isinstance(n, ast.Name) and n.id in self.defs and n.id not in need:
                    need.add(n.id)
                    stack.append(n.id)
        return [n for n in self.order if n in need]

    def emit(self, e):
        """(rust code, 'S' scalar of type T or 'V' series)"""
        if isinstance(e, ast.Constant):
            v = e.value
            if isinstance(v, bool) or not isinstance(v, (int, float)):
                raise ValueError(ast.dump(e))
            return f"c::<T>({float(v)!r})", "S"
        if isinstance(e, ast.Name):
            if e.id in SCALARS:
                return SCALARS[e.id], "S"
            if e.id in SERIES:
                return e.id, "V"
            return e.id, self.types[e.id]
        if isinstance(e, ast.UnaryOp) and isinstance(e.op, ast.USub):
            a, ta = self.emit(e.operand)
            return (f"(-{a})", "S") if ta == "S" else (f"neg(&{a})", "V")
        if isinstance(e, ast.BinOp):
            if isinstance(e.op, ast.Pow):
                a, ta = self.emit(e.left)
                ex = ast.literal_eval(e.right)
                if ex == -1:
                    return (f"(T::one() / {a})", "S") if ta == "S" else (f"recip(&{a})", "V")
                assert isinstance(ex, int) and ex >= 2, ex
                return (f"pows({a}, {ex})", "S") if ta == "S" else (f"powi(&{a}, {ex})", "V")
            a, ta = self.emit(e.left)
            b, tb = self.emit(e.right)
            k = ta + tb
            if isinstance(e.op, ast.Add):
                return {"SS": (f"({a} + {b})", "S"), "VV": (f"add(&{a}, &{b})", "V"),
                        "VS": (f"adds(&{a}, {b})", "V"), "SV": (f"adds(&{b}, {a})", "V")}[k]
            if isinstance(e.op, ast.Sub):
                return {"SS": (f"({a} - {b})", "S"), "VV": (f"sub(&{a}, &{b})", "V"),
                        "VS": (f"adds(&{a}, -{b})", "V"), "SV": (f"adds(&neg(&{b}), {a})", "V")}[k]
            if isinstance(e.op, ast.Mult):
                return {"SS": (f"({a} * {b})", "S"), "VV": (f"mul(&{a}, &{b})", "V"),
                        "VS": (f"scl(&{a}, {b})", "V"), "SV": (f"scl(&{b}, {a})", "V")}[k]
            if isinstance(e.op, ast.Div):
                return {"SS": (f"({a} / {b})", "S"), "VV": (f"mul(&{a}, &recip(&{b}))", "V"),
                        "VS": (f"scl(&{a}, T::one() / {b})", "V"), "SV": (f"scl(&recip(&{b}), {a})", "V")}[k]
        raise ValueError(ast.dump(e))

    def render(self):
        o = ["//! AUTO-GENERATED by `taylor_polar.py --gen` from polar_matrix_generated.py (which reduce_polar.py wrote).",
             "//! J(r; s) of y' = J y, y = (H1/s, K, F, F', G/s, G'/s), as operations on truncated power series.",
             "//! The symbolic source is shared with polar_solver.py. Do not edit.",
             "#![allow(clippy::all, unused_variables, unused_parens)]",
             "use crate::num::Num;", "use crate::ser::*;", "",
             "#[derive(Clone, Copy, Debug)]", "pub struct Pars<T: Num> {"]
        o += [f"    pub {k}: T," for k in ("omega", "epsilon", "s", "ll", "lambda", "mu", "pi")]
        o += ["}", "", "fn c<T: Num>(x: f64) -> T {", "    T::from_f64(x)", "}", "",
              "fn pows<T: Num>(a: T, e: i32) -> T {", "    let mut r = T::one();", "    for _ in 0..e {",
              "        r = r * a;", "    }", "    r", "}", "",
              "pub fn j_series<T: Num>(r: &Ser<T>, alpha: &Ser<T>, a: &Ser<T>, phi0: &Ser<T>, u: &Ser<T>, p: &Pars<T>) -> [[Ser<T>; 6]; 6] {",
              "    let (r, alpha, a, phi0, u) = (*r, *alpha, *a, *phi0, *u);",
              "    let (omega, epsilon, s, ll, lambda, mu, pi) = (p.omega, p.epsilon, p.s, p.ll, p.lambda, p.mu, p.pi);"]
        for n in self.needed():
            code, ty = self.emit(self.defs[n])
            self.types[n] = ty
            o.append(f"    let {n} = {code};")
        o.append("    let mut j = [[Ser::<T>::zero(); 6]; 6];")
        for i, jj, e in self.entries:
            if isinstance(e, ast.Constant) and e.value == 0:
                continue
            code, ty = self.emit(e)
            o.append(f"    j[{i}][{jj}] = " + (f"Ser::cst({code});" if ty == "S" else f"{code};"))
        o += ["    j", "}"]
        return "\n".join(o) + "\n"


def gen():
    tree = ast.parse(open(SRC).read())
    fn = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "J_and_h"][0]
    g = Gen(fn)
    text = g.render()
    path = os.path.join(CRATE, "src", "jgen.rs")
    open(path, "w").write(text)
    print(f"wrote {path}: {len(text.splitlines())} lines, {len(g.needed())} temporaries, {len(g.entries)} entries")


def build(features=None):
    env = dict(os.environ, CARGO_NET_GIT_FETCH_WITH_CLI="true")
    subprocess.check_call(["nice", "-n", "10", "cargo", "build", "--release", "-j", "2"] + (["--features", features] if features else []), cwd=CRATE, env=env)


# ---------------------------------------------------------------------------------------------------------------
# running the crate
# ---------------------------------------------------------------------------------------------------------------
MUTATE = None
SCIPY_FILE = None
_star_cache = {}


def rust(cmd, whole=False, **kw):
    args = [BIN, cmd] + [f"{k}={v}" for k, v in kw.items()]
    if MUTATE:
        args.append(f"mutate={MUTATE}")
    env = dict(os.environ, RAYON_NUM_THREADS="2")
    if MUTATE:
        env["TP_MAX_SECONDS"] = "60"   # a mutated integration that does not finish is a failure
    p = subprocess.run(["nice", "-n", "10"] + args, capture_output=True, text=True, env=env, cwd=HERE)
    if p.returncode:
        raise RuntimeError(f"{' '.join(args)} exited {p.returncode}: {p.stderr[-500:]}")
    if whole:
        return json.loads(p.stdout)
    return [json.loads(x) for x in p.stdout.splitlines() if x.startswith("{")]


def star_data(rs=None):
    """central data of the star by Newton shooting in the crate (nothing from scipy)"""
    key = rs
    if key not in _star_cache:
        kw = dict(MODEL)
        if rs is not None:
            kw["rs"] = rs
        _star_cache[key] = rust("star", **kw)[0]
    return _star_cache[key]


def star_kw(rs=None):
    s = star_data(rs)
    kw = dict(MODEL, pc=s["pc"], zc=s["zc"])
    if rs is not None:
        kw["rs"] = rs
    return kw


RESULTS = []


def report(name, ok, text):
    """ok is True (passed), False (failed) or None (skipped: it did not run, and counts as neither)"""
    ok = None if ok is None else bool(ok)   # a numpy bool (S2's) is neither `is True` nor `is False`
    RESULTS.append((name, ok))
    print(f"CHECK {name}: {'SKIPPED' if ok is None else 'PASS' if ok else 'FAIL'}  {text}", flush=True)


def fmt(x):
    return f"{float(x):.2e}"


# ---------------------------------------------------------------------------------------------------------------
# mpmath evaluation of the generated function (the same source, evaluated at 60 digits: shares the symbolic J,
# shares no arithmetic with the crate)
# ---------------------------------------------------------------------------------------------------------------
_mpJ = {}


def mp_J():
    import mpmath as mp
    import numpy as np
    if "f" in _mpJ:
        return _mpJ["f"]
    mp.mp.dps = 60
    tree = ast.parse(open(SRC).read())
    fn = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "J_and_h"][0]

    class Tr(ast.NodeTransformer):
        def _c(self, n):
            if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)) and not isinstance(n.value, bool):
                return ast.copy_location(ast.Call(ast.Name("MPF", ast.Load()), [ast.Constant(repr(n.value))], []), n)
            return n

        def visit_BinOp(self, n):
            self.generic_visit(n)
            n.left, n.right = self._c(n.left), self._c(n.right)
            return n

        def visit_UnaryOp(self, n):
            self.generic_visit(n)
            n.operand = self._c(n.operand)
            return n

        def visit_Assign(self, n):
            self.generic_visit(n)
            if isinstance(n.targets[0], ast.Name) and n.targets[0].id == "J" and isinstance(n.value, ast.Call):
                return ast.parse("J = np.empty((6, 6), dtype=object)").body[0]
            n.value = self._c(n.value)
            return n

    fn = ast.fix_missing_locations(Tr().visit(fn))
    mod = ast.Module(body=[fn], type_ignores=[])
    ns = {"np": np, "pi": mp.pi, "exp": mp.exp, "sqrt": mp.sqrt, "log": mp.log, "MPF": mp.mpf}
    exec(compile(ast.fix_missing_locations(mod), SRC, "exec"), ns)
    _mpJ["f"] = ns["J_and_h"]
    return _mpJ["f"]


def mpf(x):
    import mpmath as mp
    return mp.mpf(x) if isinstance(x, str) else mp.mpf(float(x))


# ---------------------------------------------------------------------------------------------------------------
# T1, T2: the translation of J and the series arithmetic
# ---------------------------------------------------------------------------------------------------------------
T1_POINTS = [  # l, r, alpha, a, phi0, u, s
    (2, 2.3, 1.01, 1.05, 0.01, -0.05, 0.0), (2, 0.7, 1.009, 1.002, 0.0195, -0.012, 1.0),
    (2, 11.0, 1.0006, 1.0007, 3.1e-6, -0.65, 40.0), (6, 2.3, 1.01, 1.05, 0.01, -0.05, 0.3),
    (6, 5.5, 1.003, 1.01, 2.5e-3, -0.4, 300.0), (6, 24.1, 1.0006, 0.9999, 1.0e-9, -0.65, 12000.0),
    (3, 1.5, 1.008, 1.03, 0.015, -0.1, 3.0), (4, 9.0, 1.001, 1.004, 2.0e-4, -0.62, 10.0),
]


def check_T1():
    import mpmath as mp
    f = mp_J()
    mp.mp.dps = 60
    worst, where = 0.0, None
    for (l, r, al, a, ph, u, s) in T1_POINTS:
        out = rust("jcheck", l=l, r=r, alpha=al, a=a, phi0=ph, u=u, s=s, **MODEL)
        ref, _ = f(mp.mpf(r), mp.mpf(al), mp.mpf(a), mp.mpf(ph), mp.mpf(u), mp.mpf(OMEGA), mp.mpf(EPS), mp.mpf(s),
                   mp.mpf(l * (l + 1)), mp.mpf(LAM), mp.mpf(MU))
        for o in out:
            i, j = o["i"], o["j"]
            v, w = mp.mpf(o["v"]), ref[i, j]
            if w == 0 and v == 0:
                continue
            e = float(abs(v - w) / max(abs(w), mp.mpf("1e-300")))
            if e > worst:
                worst, where = e, (l, r, s, i, j)
    ok = worst <= CRIT["T1_TOL"]
    report("T1", ok, f"J at order 0, crate (dd) against the generated function at 60 digits, {len(T1_POINTS)} points x 36 entries: "
                     f"worst relative {fmt(worst)} at (l, r, s, i, j) = {where} (bar {fmt(CRIT['T1_TOL'])})")


def check_T2():
    import mpmath as mp
    f = mp_J()
    mp.mp.dps = 60
    worst, tot = 0.0, 0
    for (l, r, s) in [(2, 0.9, 1.0), (6, 3.0, 300.0), (6, 17.0, 12000.0), (3, 7.0, 0.3)]:
        o = rust("jjets", whole=True, l=l, r=r, s=s, **star_kw())
        nb = o["nb"]
        h = o["x1"] - o["center"]
        delta = mp.mpf(h) * mp.mpf("0.1")
        horner = lambda c: mp.polyval([mp.mpf(x) for x in c[::-1]], delta)
        vals = [horner(o[k]) for k in ("r", "alpha", "a", "phi0", "u")]
        ref, _ = f(*vals, mp.mpf(OMEGA), mp.mpf(EPS), mp.mpf(s), mp.mpf(l * (l + 1)), mp.mpf(LAM), mp.mpf(MU))
        for key, c in o["J"].items():
            i, j = map(int, key.split(","))
            w, v = ref[i, j], horner(c)
            e = float(abs(v - w) / max(abs(w), mp.mpf("1e-300")))
            worst, tot = max(worst, e), tot + 1
    ok = worst <= CRIT["T2_TOL"]
    report("T2", ok, f"J's series (30 coefficients) at 0.1 of a star step against the generated function at 60 digits on the "
                     f"series' own background, 4 cases, {tot} entries: worst relative {fmt(worst)} (bar {fmt(CRIT['T2_TOL'])})")


def check_R1():
    import mpmath as mp
    mp.mp.dps = 200
    res = []
    ok = True
    for lane in ("qd", "od"):
        o = rust("recipcheck", lane=lane)[0]
        w = max(float(abs(mp.mpf(x) * (k + 1) - 1)) for k, x in enumerate(o["inv"]))
        res.append(f"{lane} {fmt(w)} (bar {fmt(CRIT['R1_TOL'][lane])})")
        ok &= w <= CRIT["R1_TOL"][lane]
    report("R1", ok, "1/k for k = 1..32 formed in the lane, |k (1/k) - 1| at 200 digits: " + ", ".join(res))


# ---------------------------------------------------------------------------------------------------------------
# P: planted systems
# ---------------------------------------------------------------------------------------------------------------
Y0 = [[1, 0, 0, 0.5, 0, 0.125], [0, 1, 0, 0.25, 0.375, 0], [0, 0, 1, 0, 0.125, 0.375]]
Z0 = [[0.25, 0, 0.125, 1, 0, 0], [0.125, 0.375, 0, 0, 1, 0], [0, 0.5, 0.25, 0, 0, 1]]
DELTA = [0, 0.5, 0.75]


def planted_A(gamma, linear_r=None):
    import mpmath as mp
    A = mp.matrix(6, 6)
    for i in range(6):
        for j in range(6):
            A[i, j] = mp.mpf(((7 * i + 3 * j + 1) % 11) - 5) / 8
            if linear_r is not None:
                A[i, j] += mp.mpf(linear_r) * mp.mpf(((5 * i + 2 * j + 3) % 7) - 3) / 16
    for k in range(3):
        g = mp.mpf(gamma) + mp.mpf(DELTA[k])
        A[2 * k, 2 * k + 1] += 1
        A[2 * k + 1, 2 * k] += g * g
    return A


def gs_cols(cols):
    """positive-diagonal Gram-Schmidt of real 6-vectors (mp); returns Q columns and the ln of the diagonal of R"""
    import mpmath as mp
    q, lnd = [], []
    for c in cols:
        v = mp.matrix(c)
        for _ in range(2):
            for qq in q:
                v = v - mp.fsum(qq[i] * v[i] for i in range(6)) * qq
        n = mp.sqrt(mp.fsum(x * x for x in v))
        q.append(v / n)
        lnd.append(mp.log(n))
    return q, lnd


def converged_plane(A, X0, t, want_growing):
    """The plane expm(A t) X0 once it has converged to the invariant subspace of the three eigenvalues of largest
    (want_growing) or smallest real part, exactly: Q (positive-diagonal Gram-Schmidt of the real eigenbasis times the
    orientation sign det M), the projector, and ln vol. The decaying planes of the inward integration pass t = -Lz."""
    import mpmath as mp
    E, V = mp.eig(A)
    Vi = V ** -1
    idx = sorted(range(6), key=lambda k: -mp.re(E[k]) if want_growing else mp.re(E[k]))[:3]
    tol = mp.mpf(10) ** -40
    real = [k for k in idx if abs(mp.im(E[k])) < tol]
    cplx = [k for k in idx if abs(mp.im(E[k])) >= tol and mp.im(E[k]) > 0]
    assert len(real) + 2 * len(cplx) == 3, [E[k] for k in idx]
    B, rows, lnscale = [], [], mp.mpf(0)
    for k in real:
        # a real eigenvalue's eigenvector made real: divide by its largest component (Vi's row takes the inverse)
        p = max(range(6), key=lambda i: abs(V[i, k]))
        s = V[p, k]
        for i in range(6):
            V[i, k] = V[i, k] / s
        for i in range(6):
            Vi[k, i] = Vi[k, i] * s
        B.append([mp.re(V[i, k]) for i in range(6)])
        c0 = [mp.re(sum(Vi[k, i] * mp.mpf(X0[col][i]) for i in range(6))) for col in range(3)]
        rows.append(c0)
        lnscale += mp.re(E[k]) * t
    for k in cplx:
        B.append([mp.re(V[i, k]) for i in range(6)])
        B.append([mp.im(V[i, k]) for i in range(6)])
        cc = [sum(Vi[k, i] * mp.mpf(X0[col][i]) for i in range(6)) for col in range(3)]
        ph = mp.expj(mp.im(E[k]) * t)
        rows.append([2 * mp.re(ph * c) for c in cc])
        rows.append([-2 * mp.im(ph * c) for c in cc])
        lnscale += 2 * mp.re(E[k]) * t
    M = mp.matrix(rows)
    Q, lnd = gs_cols(B)
    sgn = 1 if mp.det(M) > 0 else -1
    lnvol = lnscale + mp.log(abs(mp.det(M))) + mp.fsum(lnd)
    P = mp.zeros(6, 6)
    for q in Q:
        for i in range(6):
            for j in range(6):
                P[i, j] += q[i] * q[j]
    return Q, sgn, P, lnvol


def projector(cols):
    import mpmath as mp
    P = mp.zeros(6, 6)
    for q in cols:
        for i in range(6):
            for j in range(6):
                P[i, j] += mp.mpf(q[i]) * mp.mpf(q[j])
    return P


def planted_compare(rs, ref_D, ref_lny, ref_lnz, ref_PY, ref_PZ):
    import mpmath as mp
    r = rs[0]
    qy = [[mp.mpf(x) for x in o["v"]] for o in rs[1:] if o["name"] == "y"]
    qz = [[mp.mpf(x) for x in o["v"]] for o in rs[1:] if o["name"] == "z"]
    dD = abs(mp.mpf(r["d"]) - ref_D) / abs(ref_D)
    dP = max(max(abs(projector(qy)[i, j] - ref_PY[i, j]) for i in range(6) for j in range(6)),
             max(abs(projector(qz)[i, j] - ref_PZ[i, j]) for i in range(6) for j in range(6)))
    dLN = max(abs(mp.mpf(r["lnvol_y"]) - ref_lny) / abs(ref_lny), abs(mp.mpf(r["lnvol_z"]) - ref_lnz) / abs(ref_lnz))
    return float(dD), float(dP), float(dLN)


def check_P1(lanes=("f64", "dd", "qd")):
    """long L: the planes have converged to the invariant subspaces; D, the planes and the log-volume ~ e^(gamma L)"""
    import mpmath as mp
    mp.mp.dps = 90
    L, Lz = CRIT["P_L"], CRIT["P_LZ"]
    for lane in lanes:
        worst = dict(D=0.0, P=0.0, LN=0.0)
        info = []
        for g in (([10.0] if lane == "qd" else [3.0]) if MUTATE == "noqr" else CRIT["P_GAMMAS"]):
            if lane == "qd" and (g > 300 or g < 10):
                continue
            Lz_eff = (L + Lz) - L
            A = planted_A(g)
            QY, sy, PY, lny = converged_plane(A, Y0, mp.mpf(L), True)
            QZ, sz, PZ, lnz = converged_plane(A, Z0, -mp.mpf(Lz_eff), False)
            M = mp.matrix(6, 6)
            for c in range(3):
                for i in range(6):
                    M[i, c], M[i, 3 + c] = QY[c][i], QZ[c][i]
            D = mp.det(M) * sy * sz
            rs = rust("planted", case="const", lane=lane, gamma=g, L=L, Lz=Lz)
            dD, dP, dLN = planted_compare(rs, D, lny, lnz, PY, PZ)
            for k, v in (("D", dD), ("P", dP), ("LN", dLN)):
                worst[k] = max(worst[k], v)
            info.append(f"g={g:g}: D={mp.nstr(D, 6)} steps {rs[0]['steps_y']}+{rs[0]['steps_z']}")
        tol = CRIT["P_TOL"][lane]
        ok = all(worst[k] <= tol[k] for k in worst)
        report(f"P1[{lane}]", ok, f"y' = A y, A = blockdiag([[0,1],[g^2,0]]) + C, L = {L}, Lz = {Lz} against the invariant subspaces "
                                 f"at 90 digits; worst D {fmt(worst['D'])}, planes {fmt(worst['P'])}, ln vol {fmt(worst['LN'])} "
                                 f"(bars {fmt(tol['D'])}, {fmt(tol['P'])}, {fmt(tol['LN'])}); " + "; ".join(info))


def check_P0(lanes=("f64", "dd")):
    """transient: gamma L = 40, mpmath's expm at 130 digits, so the dependence on the start data is tested"""
    import mpmath as mp
    GL = CRIT["P_TRANSIENT_GL"]
    for lane in lanes:
        worst = dict(D=0.0, P=0.0, LN=0.0)
        for g in [3.0, 300.0, 12000.0]:
            mp.mp.dps = 140
            L = GL / g
            Lz = (L + 0.5 * L) - L
            A = planted_A(g)
            E, V = mp.eig(A)
            Vi = V ** -1
            ex = lambda t: V * mp.diag([mp.exp(e * t) for e in E]) * Vi
            Y = [[mp.re(v) for v in ex(mp.mpf(L)) * mp.matrix([mp.mpf(x) for x in c])] for c in Y0]
            Z = [[mp.re(v) for v in ex(-mp.mpf(Lz)) * mp.matrix([mp.mpf(x) for x in c])] for c in Z0]
            qy, dy = gs_cols(Y)
            qz, dz = gs_cols(Z)
            M = mp.matrix(6, 6)
            for c in range(3):
                for i in range(6):
                    M[i, c], M[i, 3 + c] = qy[c][i], qz[c][i]
            D = mp.det(M)
            rs = rust("planted", case="const", lane=lane, gamma=g, L=L, Lz=0.5 * L)
            dD, dP, dLN = planted_compare(rs, D, mp.fsum(dy), mp.fsum(dz), projector(qy), projector(qz))
            for k, v in (("D", dD), ("P", dP), ("LN", dLN)):
                worst[k] = max(worst[k], v)
        tol = CRIT["P_TOL"][lane]
        ok = all(worst[k] <= tol[k] for k in worst)
        report(f"P0[{lane}]", ok, f"transient, gamma L = {GL:g}, expm at 140 digits, gamma = 3, 300, 12000: worst D {fmt(worst['D'])}, "
                                 f"planes {fmt(worst['P'])}, ln vol {fmt(worst['LN'])} (bars {fmt(tol['D'])}, {fmt(tol['P'])}, {fmt(tol['LN'])})")


def check_P2(lanes=("dd",)):
    """variable coefficients A(r) = A0 + r A1: mpmath's odefun (Taylor series at 40 digits, no shared code)"""
    import mpmath as mp
    mp.mp.dps = 60
    worst = dict(D=0.0, P=0.0, LN=0.0)
    for lane in lanes:
        for g in (3.0, 6.0):
            L, Lz = 1.0, 0.5
            cols = {}
            for name, X0, x0, direction in (("y", Y0, 0.0, 1), ("z", Z0, L + Lz, -1)):
                out = []
                for col in X0:
                    if direction > 0:
                        F = lambda x, y: list(planted_A(g, linear_r=x) * mp.matrix(y))
                        sol = mp.odefun(F, x0, [mp.mpf(v) for v in col])
                        out.append(list(sol(mp.mpf(L))))
                    else:
                        F = lambda tau, y: list(-(planted_A(g, linear_r=(x0 - tau)) * mp.matrix(y)))
                        sol = mp.odefun(F, 0, [mp.mpf(v) for v in col])
                        out.append(list(sol(mp.mpf(Lz))))
                cols[name] = out
            qy, dy = gs_cols(cols["y"])
            qz, dz = gs_cols(cols["z"])
            M = mp.matrix(6, 6)
            for c in range(3):
                for i in range(6):
                    M[i, c], M[i, 3 + c] = qy[c][i], qz[c][i]
            D = mp.det(M)
            rs = rust("planted", case="linear", lane=lane, gamma=g, L=L, Lz=Lz)
            dD, dP, dLN = planted_compare(rs, D, mp.fsum(dy), mp.fsum(dz), projector(qy), projector(qz))
            for k, v in (("D", dD), ("P", dP), ("LN", dLN)):
                worst[k] = max(worst[k], v)
    tol = dict(D=1e-28, P=1e-28, LN=1e-28)
    ok = all(worst[k] <= tol[k] for k in worst)
    report("P2[dd]", ok, f"A(r) = A0 + r A1, gamma = 3, 6, L = 1, mpmath odefun at 60 digits: worst D {fmt(worst['D'])}, planes "
                         f"{fmt(worst['P'])}, ln vol {fmt(worst['LN'])} (bar {fmt(tol['D'])})")


# ---------------------------------------------------------------------------------------------------------------
# S: the star
# ---------------------------------------------------------------------------------------------------------------
_scipy_star = {}


def scipy_star():
    import background as B
    if "st" not in _scipy_star:
        _scipy_star["st"] = B.solve_star(OMEGA, EPS)
    return _scipy_star["st"]


def star_values(kw, radii):
    return {o["r"]: o for o in rust("starcheck", r=",".join(repr(r) for r in radii), **kw)}


def check_S1():
    import numpy as np
    import mpmath as mp
    st = scipy_star()
    radii = CRIT["S_RADII"]
    res = {}
    for label, rs in (("exact", None), ("matched", 1e-3)):
        kw = star_kw(rs)
        vals = star_values(kw, radii)
        worst, where = 0.0, None
        for r in radii:
            al, a, ph, u = st.bg(np.array([r]))
            ref = dict(alpha=float(al[0]), a=float(a[0]), phi0=float(ph[0]), u=float(u[0]))
            for k, v in ref.items():
                e = abs(float(mp.mpf(vals[r][k])) / v - 1.0)
                if e > worst:
                    worst, where = e, (k, r)
        res[label] = (worst, where)
    ok = res["exact"][0] <= CRIT["S1_TOL"] and res["matched"][0] <= CRIT["S1_TOL_MATCHED"]
    report("S1", ok, f"star by shooting against solve_star's mesh at r = {radii[0]} .. {radii[-1]} ({len(radii)} radii, alpha, a, phi0, u): "
                     f"exact regular start worst relative {fmt(res['exact'][0])} at {res['exact'][1]} (bar {fmt(CRIT['S1_TOL'])}); "
                     f"start at r0 = 1e-3 with scipy's leading-order boundary condition {fmt(res['matched'][0])} at {res['matched'][1]} "
                     f"(bar {fmt(CRIT['S1_TOL_MATCHED'])})")


def check_S2():
    """the crate's Taylor recursion of the star's equations against scipy's DOP853 on background.odes, from the crate's state"""
    import numpy as np
    import mpmath as mp
    from scipy.integrate import solve_ivp
    import background as B
    kw = star_kw()
    ra, rb = 0.5, 8.0
    v = star_values(kw, [ra, rb])
    y = lambda o: [float(mp.log(mp.mpf(o["phi0"]))), float(o["u"]), float(o["m"]), 2 * float(mp.log(mp.mpf(o["alpha"])))]
    ya, yb = y(v[ra]), y(v[rb])
    sol = solve_ivp(lambda r, yy: B.odes(r, yy, OMEGA, EPS).ravel(), (ra, rb), ya, method="DOP853", rtol=1e-13, atol=1e-30)
    e = max(abs(sol.y[k, -1] - yb[k]) / max(abs(yb[k]), 1e-3) for k in range(4))
    ok = e <= CRIT["S2_TOL"]
    report("S2", ok, f"star's equations: crate state at r = {rb} against DOP853 (rtol 1e-13) on background.odes from the crate's state at "
                     f"r = {ra}: worst of (ln phi0, u, m, nu) {fmt(e)} (bar {fmt(CRIT['S2_TOL'])})")


def check_S3():
    """The shooting starts from the leading-order regular series at r = rs, whose error is O(rs^2): the central data
    must move by less than S3_TOL between rs = 1e-8 and 1e-9, and the move must scale as rs^2: the ratio of that
    difference to the one between rs = 1e-6 and 1e-7 is 1e-4 within a factor of 3 (a wrong start term breaks the law)."""
    from mpmath import mp, mpf
    mp.dps = 50
    d = lambda r1, r2, k: abs(mpf(star_data(r1)[k]) / mpf(star_data(r2)[k]) - 1)
    near = {k: d(1e-8, 1e-9, k) for k in ("pc", "zc")}
    far = {k: d(1e-6, 1e-7, k) for k in ("pc", "zc")}
    ratio = {k: float(near[k] / far[k]) for k in near}
    ok = all(float(near[k]) <= CRIT["S3_TOL"] for k in near)
    ok &= all(CRIT["S3_RATIO"][0] <= ratio[k] <= CRIT["S3_RATIO"][1] for k in ratio)
    report("S3", ok, f"central data (P_c, Z_c) against the shooting's start radius: rs = 1e-8 vs 1e-9 moves them by {fmt(near['pc'])}, "
                     f"{fmt(near['zc'])} (bar {fmt(CRIT['S3_TOL'])}); rs = 1e-6 vs 1e-7 by {fmt(far['pc'])}, {fmt(far['zc'])}; ratio "
                     f"{fmt(ratio['pc'])}, {fmt(ratio['zc'])} (rs^2 law: 1e-4, bar [{CRIT['S3_RATIO'][0]:.2e}, {CRIT['S3_RATIO'][1]:.1e}]); "
                     f"phi_c = {star_data(1e-8)['phic']:.12g}, alpha_c = {star_data(1e-8)['alpha_c']:.12g} "
                     f"(solve_star: {scipy_star().phic:.12g}, {scipy_star().alpha_c:.12g})")


def check_T3():
    """the table's order-0 jets (J0 + s^2 J1 shifted from the start of the table step that holds r) against J evaluated
    directly from the star at the same r, at 150 radii from 1e-3 to 23. Gated: each entry's error relative to itself.
    Printed, not gated: the error relative to the largest |J| entry of that matrix. A table step wider than J's radius
    of convergence (J has 1/r^2 terms the star's regular series lacks) breaks it."""
    worst_n, worst_e, at_n, at_e = 0.0, 0.0, None, None
    radii = [r for r in CRIT["T3_RADII"] if r < 23.0]
    for l, sv in ((2, 0.3), (6, 300.0)):
        out = rust("tabcheck", l=l, s=sv, r=",".join(repr(r) for r in radii), **star_kw())
        for o in out:
            if o["norm"] > worst_n:
                worst_n, at_n = o["norm"], (l, sv, o["r"])
            if o["worst"] > worst_e:
                worst_e, at_e = o["worst"], (l, sv, o["r"], o["at"])
    ok = worst_e <= CRIT["T3_TOL"]
    report("T3", ok, f"table against direct J at {len(radii)} radii, l = 2 (s = 0.3) and l = 6 (s = 300), dd: worst entry relative to itself "
                     f"{fmt(worst_e)} at (l, s, r, entry) = {at_e} (bar {fmt(CRIT['T3_TOL'])}); info, not gated: worst error relative to the "
                     f"largest |J| entry {fmt(worst_n)} at (l, s, r) = {at_n}")


# ---------------------------------------------------------------------------------------------------------------
# A: D against the stored rows
# ---------------------------------------------------------------------------------------------------------------
def stored_rows(l):
    path = os.path.join(HERE, f"convergence_l{l}_output.txt")
    rows = {}
    for line in open(path):
        if not line.startswith("gamma="):
            continue
        g = float(line.split(":")[0].split("=")[1])
        d = float(line.split("D_base=")[1].split()[0])
        sp = float(line.split("rel.spread (same r_m)=")[1].split()[0])
        d34 = float(line.split("D(r_m~34)=")[1].split()[0])
        rows[g] = (d, sp, d34)
    return rows


def d_values(l, gammas, lane="dd", rm=None, r0=1e-3, **extra):
    kw = star_kw()
    rm = R_M["base"] if rm is None else rm
    out = rust("d", lane=lane, l=l, r0=r0, rm=repr(rm), gammas=",".join(repr(g) for g in gammas), **kw, **extra)
    setup = out[0]
    return setup, {o["gamma"]: o for o in out[1:]}


def check_A(l, gammas, name, lane="dd", extra_scipy=None):
    rows = stored_rows(l)
    if extra_scipy:
        rows = dict(rows, **extra_scipy)
    setup, vals = d_values(l, gammas, lane)
    worst, lines, ok = 0.0, [], True
    for g in gammas:
        if g not in rows:
            continue
        d_ref, sp = rows[g][:2]
        d = vals[g]["d_f64"]
        e = abs(d / d_ref - 1.0)
        worst = max(worst, e)
        ok &= e <= CRIT["A_BAR"]
        lines.append(f"g={g:g} D={d:.7e} stored {d_ref:.7e} rel {fmt(e)} (row spread {sp:.1e}) {vals[g]['wall']:.2f}s")
    report(name, ok, f"l = {l}, base variant, lane {lane}: worst |D/D_stored - 1| = {fmt(worst)} over {len(lines)} rows (bar {CRIT['A_BAR']:.1e}, "
                     f"the file's relative spread); " + "; ".join(lines))


def check_L1():
    """which lane D needs: f64 against dd (gated by the file's spread), dd against qd (gated at L1_QD_TOL)"""
    import mpmath as mp
    mp.mp.dps = 80
    lines, ok = [], True
    for (l, g) in CRIT["LANE_GAMMAS"]:
        dd = d_values(l, [g], "dd")[1][g]
        f = d_values(l, [g], "f64")[1][g]
        e = abs(f["d_f64"] / dd["d_f64"] - 1)
        ok &= e <= CRIT["A_BAR"]
        lines.append(f"(l={l}, g={g:g}): f64 {f['d_f64']:.10e} vs dd {dd['d_f64']:.10e} rel {fmt(e)} (bar {CRIT['A_BAR']:.1e}), {f['wall']:.3f} s vs {dd['wall']:.2f} s")
    for (l, g) in CRIT["LANE_QD"]:
        dd = d_values(l, [g], "dd")[1][g]
        q = d_values(l, [g], "qd")[1][g]
        e = abs(mp.mpf(dd["d"]) / mp.mpf(q["d"]) - 1)
        ok &= float(e) <= CRIT["L1_QD_TOL"]
        lines.append(f"(l={l}, g={g:g}): dd vs qd rel {fmt(e)} (bar {fmt(CRIT['L1_QD_TOL'])}), qd {q['wall']:.1f} s, {q['steps_y']} steps")
    report("L1", ok, "lanes: " + "; ".join(lines))


def check_V1():
    """r0 = 3e-3 (r_m = 24.111): D against the base variant's D from the same code; the file's spread is the bar"""
    lines, ok = [], True
    for (l, g) in CRIT["V1_POINTS"]:
        b = d_values(l, [g], "dd")[1][g]["d_f64"]
        v = d_values(l, [g], "dd", rm=R_M["r0"], r0=3e-3)[1][g]["d_f64"]
        e = abs(v / b - 1)
        ok &= e <= CRIT["A_BAR"]
        lines.append(f"(l={l}, g={g:g}) {fmt(e)}")
    report("V1", ok, f"variant r0 = 3e-3 against base, same code, dd: {'; '.join(lines)} (bar {CRIT['A_BAR']:.1e})")


def check_V2():
    """r_m = 34.2051 (phi_cut 1e-12): D against the stored D(r_m~34) column, which has 4 digits: the bar is the file's spread
    plus the half unit of the printed digits"""
    lines, ok = [], True
    for (l, g) in CRIT["V2_POINTS"]:
        row = stored_rows(l)[g]
        d34 = row[2]
        mant = abs(d34) / 10 ** math.floor(math.log10(abs(d34)))
        bar = CRIT["A_BAR"] + CRIT["V2_DIGITS_ALLOWANCE"] / mant
        v = d_values(l, [g], "dd", rm=R_M["phi_cut"])[1][g]["d_f64"]
        e = abs(v / d34 - 1)
        ok &= e <= bar
        lines.append(f"(l={l}, g={g:g}) D={v:.4e} stored {d34:.4e} rel {fmt(e)} (bar {bar:.1e})")
    report("V2", ok, "variant r_m = 34.2051, dd: " + "; ".join(lines))


def check_V3():
    """the Taylor order and tolerance, the analog of the rtol variants: order 24 at tol 1e-22 against order 30 at 1.2e-30"""
    lines, ok = [], True
    for (l, g) in CRIT["V3_POINTS"]:
        b = d_values(l, [g], "dd")[1][g]
        v = d_values(l, [g], "dd", order=CRIT["V3_ORDER"], tol=CRIT["V3_TOL"])[1][g]
        import mpmath as mp
        mp.mp.dps = 60
        e = abs(mp.mpf(v["d"]) / mp.mpf(b["d"]) - 1)
        ok &= float(e) <= CRIT["V3_BAR"]
        lines.append(f"(l={l}, g={g:g}) {fmt(e)} steps {b['steps_y']} -> {v['steps_y']}")
    report("V3", ok, f"order {CRIT['V3_ORDER']}, tol {CRIT['V3_TOL']:g} against order 30, tol 1.2e-30, dd: {'; '.join(lines)} (bar {CRIT['V3_BAR']:.0e})")


def check_V4():
    """convergence.py's variants through scipy here and through the crate at the same r0 and r_m (r_m is read from each scipy
    problem): R = 90 (r_m = 24.1084), r0 = 3e-3 (r_m = 24.1110) and rtol 1e-12 (r_m = 24.1098). The crate agrees with each, and
    the 2.1e-4 spread of the stored rows is the shift of r_m between them (D is normalized, so it depends on r_m), not
    integration error: rtol 1e-12 moves D by about 1e-12."""
    import background as B
    from polar_solver import PolarProblem
    st60 = scipy_star()
    st90 = B.solve_star(OMEGA, EPS, guess=st60, R=90.0)
    variants = {"R=90": (st90, {}), "r0=3e-3": (st60, dict(r0=3e-3)), "rtol 1e-12": (st60, dict(rtol=1e-12))}
    lines, ok, worst = [], True, 0.0
    base = {}
    for (l, g) in CRIT["V4_POINTS"]:
        base[(l, g)] = PolarProblem(st60, l).D(g)
        for name, (st, kw) in variants.items():
            p = PolarProblem(st, l, **kw)
            d_sc = p.D(g)
            d_cr = d_values(l, [g], "dd", rm=float(p.rm), r0=float(p.r0))[1][g]["d_f64"]
            e = abs(d_cr / d_sc - 1)
            worst = max(worst, e)
            ok &= e <= CRIT["V4_TOL"]
            lines.append(f"(l={l}, g={g:g}) {name}: r_m {float(p.rm):.4f}, scipy vs base {d_sc / base[(l, g)] - 1:+.2e}, crate vs scipy {fmt(e)}")
    report("V4", ok, f"variants through scipy and the crate at the same (r0, r_m): worst crate vs scipy {fmt(worst)} (bar {fmt(CRIT['V4_TOL'])}); " + "; ".join(lines))


def scipy_timing(out_path):
    """the scipy path of convergence.py, base variant only (rtol 1e-10, phi_cut 1e-9, R = 60, r0 = 1e-3): one evaluation
    per (l, gamma), wall and cpu seconds, appended to out_path as it finishes. Stops an evaluation at 1 h of wall time."""
    import signal
    import background as B
    from polar_solver import PolarProblem

    def handler(signum, frame):
        raise TimeoutError

    st = B.solve_star(OMEGA, EPS)
    probs = {}
    done = set()
    if os.path.exists(out_path):
        for line in open(out_path):
            o = json.loads(line)
            done.add((o["l"], o["gamma"]))
    for (l, g) in CRIT["TS_POINTS"]:
        if (l, g) in done:
            continue
        if l not in probs:
            probs[l] = PolarProblem(st, l)
        signal.signal(signal.SIGALRM, handler)
        signal.alarm(3600)
        t0, c0, load = time.time(), time.process_time(), os.getloadavg()[0]
        try:
            d = probs[l].D(g)
            rec = dict(l=l, gamma=g, D=d, wall=time.time() - t0, cpu=time.process_time() - c0, load1=load, status="done")
        except TimeoutError:
            rec = dict(l=l, gamma=g, D=None, wall=time.time() - t0, cpu=time.process_time() - c0, load1=load, status="stopped at 1 h")
        signal.alarm(0)
        with open(out_path, "a") as fp:
            fp.write(json.dumps(rec) + "\n")
        print(rec, flush=True)


def check_TS():
    """timing against the scipy path on this Mac (needs the file from --scipy-timing): D from scipy here against D from the crate
    (gated at the file's spread), and the wall time of each"""
    path = SCIPY_FILE
    if not path or not os.path.exists(path):
        report("TS", None, "no scipy timing file (--scipy-file=PATH)")
        return
    recs = [json.loads(x) for x in open(path)]
    lines, ok = [], True
    for o in recs:
        l, g = o["l"], o["gamma"]
        rows = {lane: d_values(l, [g], lane)[1][g] for lane in ("f64", "dd")}
        if o["status"] != "done":
            lines.append(f"(l={l}, g={g:g}) scipy {o['status']} after {o['wall']:.0f} s")
            continue
        e = abs(rows["dd"]["d_f64"] / o["D"] - 1)
        ok &= e <= CRIT["A_BAR"]
        lines.append(f"(l={l}, g={g:g}) scipy {o['wall']:.1f} s (cpu {o['cpu']:.1f}, load {o['load1']:.1f}) D={o['D']:.7e}; crate dd {rows['dd']['wall']:.3f} s "
                     f"(x{o['wall'] / rows['dd']['wall']:.0f}), f64 {rows['f64']['wall']:.4f} s (x{o['wall'] / rows['f64']['wall']:.0f}); D rel {fmt(e)}")
    report("TS", ok, f"scipy path (base variant, this Mac, one process) against the crate (one thread); load {os.getloadavg()[0]:.1f} now; " + "; ".join(lines))


# ---------------------------------------------------------------------------------------------------------------
# #187: the rows of convergence.py from the crate (python taylor_polar.py --rows)
# ---------------------------------------------------------------------------------------------------------------
RAW = os.path.normpath(os.path.join(HERE, "..", "..", "..", "negative-mass-raw"))
JDIR = None
G8 = [0.001, 0.03, 0.3, 1.0, 3.0, 10.0, 40.0, 137.0]
ROW_GAMMAS = {2: G8 + [300.0], 3: G8 + [300.0, 1000.0], 4: G8 + [300.0, 1000.0, 3000.0],
              5: G8 + [300.0, 1000.0, 6000.0], 6: G8 + [300.0, 1000.0, 12000.0]}
CAPS = {5: 6000.0, 6: 12000.0}
NAMES = {"base": "base (rtol1e-10, phi_cut1e-9, R60)", "tol": "rtol 1e-12", "phi": "phi_cut 1e-12 (larger r_m)",
         "r90": "background R=90", "r0": "r0 = 3e-3"}
VKEYS = list(NAMES)
HDR = dict(NAMES, tol="rtol 1e-12 [crate: Taylor tol 1e-12]")   # the crate has no scipy rtol: said in the header
CAP_BAR = MEAN_SPREAD           # #114: each capped variant within 2.1e-4 of base
CAP_VS_SCIPY = 1e-8             # crate against #114's finished l = 5 values
RP_CASES = [(1e-9, 60.0, 1e-3), (1e-12, 60.0, 1e-3), (1e-9, 90.0, 1e-3), (1e-9, 60.0, 3e-3)]
QP_BAR = {("f64", "well"): 1e-12, ("dd", "well"): 1e-28, ("qd", "well"): 1e-60,
          ("f64", "near"): 1e-5, ("dd", "near"): 1e-20, ("qd", "near"): 1e-50}
ROW_RE = re.compile(r"^gamma=\s*(?P<gamma>\S+): D_base=(?P<D>\S+)  signs=(?P<signs>\[.*?\])  rel\.spread \(same r_m\)=(?P<spread>\S+)  "
                    r"D\(r_m~34\)=(?P<D34>\S+)  q_int=(?P<qi>\S+) q_ext=(?P<qe>\S+)  \((?P<T>\d+)s\)")
HDR_RE = re.compile(r"^(?P<name>.{38}): r_m=(?P<rm>\S+)$")
_rows_cache = {}


def get_rows(l, gammas, lane, variants, tag="", rs=None):
    """The crate's `rows` for (l, gammas, variants), streamed: every JSON line goes to jsonl-dir/rows_l{l}_{lane}{tag}.jsonl as it arrives,
    so a stopped run resumes at the next gamma. Returns (setup, {(gamma, variant): record}). Records already there are reused.
    A run with a mutation reads and writes no file. rs: the start radius of the crate's star shooting (None: its default 1e-8)."""
    path = None if (MUTATE or not JDIR) else os.path.join(JDIR, f"rows_l{l}_{lane}{tag}.jsonl")
    key = (l, lane, tag, MUTATE)
    got, setup = _rows_cache.setdefault(key, ({}, [None]))
    if not got and path and os.path.exists(path):
        for line in open(path):
            if line.startswith("{"):
                o = json.loads(line)
                if o["cmd"] == "rowsetup":
                    setup[0] = o
                elif o["cmd"] == "row":
                    got[(o["gamma"], o["variant"])] = o
    need = [g for g in gammas if any((g, v) not in got for v in variants)]
    if need or setup[0] is None:
        args = [BIN, "rows", f"l={l}", f"lane={lane}", f"variants={','.join(variants)}"]
        if need:
            args.append(f"gammas={','.join(repr(g) for g in need)}")
        args += [f"{k}={v}" for k, v in star_kw(rs).items()]
        if MUTATE:
            args.append(f"mutate={MUTATE}")
        env = dict(os.environ, RAYON_NUM_THREADS="2")
        if MUTATE:
            env["TP_MAX_SECONDS"] = "900"
        fp = open(path, "a") if path else None
        p = subprocess.Popen(["nice", "-n", "10"] + args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env, cwd=HERE)
        errs = []
        for line in p.stdout:
            if not line.startswith("{"):
                errs.append(line)
                continue
            if fp:
                fp.write(line)
                fp.flush()
            o = json.loads(line)
            if o["cmd"] == "rowsetup":
                setup[0] = o
            else:
                got[(o["gamma"], o["variant"])] = o
        rc = p.wait()
        if fp:
            fp.close()
        if rc:
            raise RuntimeError(f"rows l={l} {lane} exited {rc}: {''.join(errs)[-400:]}")
    return setup[0], got


def rm_of(setup):
    return {v["variant"]: v["rm"] for v in setup["variants"]}


def format_rows(gammas, got, setup, notes=None):
    """convergence.py's header and row lines (its f-strings, its names), from the crate's records"""
    import numpy as np
    rm = rm_of(setup)
    out = [f"{HDR[k]:38s}: r_m={rm[k]:.2f}" for k in VKEYS] + [""]
    for g in gammas:
        Ds = {NAMES[k]: float(got[(g, k)]["d_f64"]) for k in VKEYS}
        qi, qe = float(got[(g, "base")]["q_int_f64"]), float(got[(g, "base")]["q_ext_f64"])
        T = sum(got[(g, k)]["wall"] for k in VKEYS)
        signs = set(np.sign(list(Ds.values())))
        same_rm = [Ds[k] for k in Ds if "phi_cut" not in k]
        spread = (max(same_rm) - min(same_rm)) / abs(np.mean(same_rm))
        out.append(f"gamma={g:8.3g}: D_base={Ds[NAMES['base']]:+.6e}  signs={sorted(signs)}  "
                   f"rel.spread (same r_m)={spread:.1e}  D(r_m~34)={Ds[NAMES['phi']]:+.3e}  "
                   f"q_int={qi:+.6e} q_ext={qe:+.6e}  ({T:.0f}s)" + (notes or {}).get(g, ""))
    return out


def stored_block(l):
    hdr, rows = [], {}
    for line in open(os.path.join(HERE, f"convergence_l{l}_output.txt")).read().splitlines():
        m = ROW_RE.match(line)
        if m:
            rows[float(m["gamma"])] = m.groupdict()
            continue
        m = HDR_RE.match(line)
        if m:
            hdr.append((m["name"].strip(), m["rm"]))
    return hdr, rows


def margin_rel(s, digits):
    """distance of the decimal string s to the nearest rounding boundary of its print with `digits` places after the
    leading digit, relative to |s| (0 would sit on a boundary)"""
    import mpmath as mp
    mp.mp.dps = 50
    x = abs(mp.mpf(s))
    e = int(mp.floor(mp.log10(x)))
    sc = x / mp.mpf(10) ** e * mp.mpf(10) ** digits
    return float(abs((sc - mp.floor(sc)) - mp.mpf(1) / 2) / sc)


def scipy_full(l):
    """scipy's full-precision values of every variant at every gamma (read only, negative-mass-raw): {(name, gamma): D or [q_int, q_ext]}"""
    out = {}
    for folder in ("I4-win", "I114-win"):
        p = os.path.join(RAW, folder, f"convergence_l{l}_w0.7603.partial.csv")
        if os.path.exists(p):
            for line in open(p):
                if line.startswith("["):
                    k, g, v, _ = json.loads(line)
                    out[(k, g)] = v
    return out


def check_RA(l):
    gammas = ROW_GAMMAS[l]
    setup, got = get_rows(l, gammas, "dd", VKEYS)
    hdr_s, rows_s = stored_block(l)
    notes = {g: "  [no stored row]" for g in gammas if g not in rows_s}
    if l in CAPS:
        notes[CAPS[l]] = "  [#114's capped row; the rtol 1e-12 column is the crate's Taylor tol 1e-12]"
    block = format_rows(gammas, got, setup, notes)
    print(f"# l = {l}: rows from the crate (lane dd, one process, nice 10, load {os.getloadavg()[0]:.1f}); setup per l: star {setup['t_star']:.1f} s, J table {setup['t_table']:.1f} s")
    for line in block:
        print(line, flush=True)
    new_rows = {float(m["gamma"]): m for m in (ROW_RE.match(x) for x in block) if m}
    mism, ntok, nok = [], 0, 0

    def cmp(where, tok, st, nw):
        nonlocal ntok, nok
        ntok += 1
        if st == nw:
            nok += 1
        else:
            mism.append(f"{where} {tok}: stored {st} crate {nw}")
    rm = rm_of(setup)
    for k, (sname, srm) in zip(VKEYS, hdr_s):
        cmp(k, "name", sname, HDR[k].split(" [")[0])   # the rtol line carries a bracketed note by design
        cmp(k, "r_m", srm, f"{rm[k]:.2f}")
    for g, st in rows_s.items():
        nw = new_rows[g].groupdict()
        for tok in ("gamma", "D", "signs", "spread", "D34", "qi", "qe"):
            cmp(f"g={g:g}", tok, st[tok], nw[tok])
    # nearest rounding boundary of the crate's unrounded values at the stored rows, and scipy's full-precision values
    near = (9.0, None)
    for g in rows_s:
        for tok, s, d in (("D", got[(g, "base")]["d"], 6), ("q_int", got[(g, "base")]["q_int"], 6), ("q_ext", got[(g, "base")]["q_ext"], 6),
                          ("D34", got[(g, "phi")]["d"], 3)):
            m = margin_rel(s, d)
            if m < near[0]:
                near = (m, (g, tok))
    info = f"closest rounding boundary of any D/q token: {near[0]:.2e} relative at (gamma, token) = {near[1]}"
    sf = scipy_full(l)
    if sf:
        worst = {}
        for g in rows_s:
            for k in VKEYS:
                ref = sf.get((NAMES[k], g))
                if ref is not None:
                    e = abs(float(got[(g, k)]["d_f64"]) / ref - 1)
                    if e > worst.get(NAMES[k], (0.0, 0))[0]:
                        worst[NAMES[k]] = (e, g)
            q = sf.get(("q", g))
            if q is not None:
                for j, kk in enumerate(("q_int", "q_ext")):
                    e = abs(float(got[(g, "base")][kk + "_f64"]) / q[j] - 1)
                    if e > worst.get(kk, (0.0, 0))[0]:
                        worst[kk] = (e, g)
        info += "; crate against scipy's full-precision values (I4-win/I114-win partial files), worst relative per column: " + ", ".join(
            f"{k.split(' (')[0]} {v[0]:.1e}@{v[1]:g}" for k, v in worst.items())
    report(f"RA{l}", nok == ntok, f"l = {l}: {nok} of {ntok} tokens of {len(rows_s)} stored rows and the 5 header lines equal as strings"
           + (f" ({len(gammas) - len(rows_s)} extra row, not compared)" if len(gammas) > len(rows_s) else "") + "; "
           + ("; ".join(mism[:12]) + ("; ..." if len(mism) > 12 else "") + "; " if mism else "") + info)


def check_RX():
    p = os.path.join(RAW, "I186-mac", "scipy_timing.jsonl")
    if not os.path.exists(p):
        report("RX", None, f"{p} not found (--raw=DIR)")
        return
    ref = [o for o in map(json.loads, open(p)) if o["l"] == 2 and o["gamma"] == 300.0 and o["status"] == "done"][0]
    setup, got = get_rows(2, [300.0], "dd", VKEYS)
    tc, ts = f"{float(got[(300.0, 'base')]['d_f64']):+.7e}", f"{ref['D']:+.7e}"
    report("RX", tc == ts, f"l = 2, gamma = 300 (no stored row): crate D_base {tc} against #186's fresh scipy D {ts} (I186-mac/scipy_timing.jsonl), 8-digit tokens equal as strings")


def i114_l5(g):
    """#114's finished l = 5 evaluations at gamma = g, all five variants, scipy full precision (read only): {name: D}, or None"""
    p = os.path.join(RAW, "I114-win", "convergence_l5_w0.7603.partial.csv")
    if not os.path.exists(p):
        return None
    sf = {}
    for line in open(p):
        if line.startswith("["):
            k, gg, v, _ = json.loads(line)
            if gg == g and k in NAMES.values():
                sf[k] = v
    return sf


def check_CAP(l):
    """the capped row of #114 at l: (a) D < 0 in the five variants, (b) the crate's rtol analogue within 2.1e-4 of base (error bars
    E1 to E3 printed), and for l = 5 (c) the crate's five D against #114's scipy D at CAP_VS_SCIPY each"""
    import mpmath as mp
    mp.mp.dps = 60
    g = CAPS[l]
    setup, got = get_rows(l, [g], "dd", VKEYS)
    D = {k: mp.mpf(got[(g, k)]["d"]) for k in VKEYS}
    neg = all(D[k] < 0 for k in VKEYS)
    e1 = abs(D["tol"] / D["base"] - 1)
    _, g24 = get_rows(l, [g], "dd", ["o24"], tag="_o24")
    e2 = abs(mp.mpf(g24[(g, "o24")]["d"]) / D["base"] - 1)
    _, g64 = get_rows(l, [g], "f64", ["base"])
    e3 = abs(mp.mpf(g64[(g, "base")]["d"]) / D["base"] - 1)
    head = f"l = {l}, gamma = {g:g}"
    report(f"CAP{l}a", neg, f"{head}: D (dd) " + ", ".join(f"{k} {mp.nstr(D[k], 12)}" for k in VKEYS) + f"; all five D < 0: {neg}")
    report(f"CAP{l}b", float(e1) <= CAP_BAR,
           f"{head}: E1 = |D(tol 1e-12)/D(base) - 1| = {float(e1):.2e} (bar {CAP_BAR:.1e}); error bars, printed not gated: E2 (order 24, tol 1e-22 against "
           f"order 30, 1.2e-30) = {float(e2):.2e}, E3 (f64 lane against dd) = {float(e3):.2e}; crate time per variant (dd, load {os.getloadavg()[0]:.1f}, "
           f"one thread): " + ", ".join(f"{k} {got[(g, k)]['wall']:.1f} s" for k in VKEYS))
    if l == 5:
        sf = i114_l5(g)
        if sf is None:
            report("CAP5c", None, f"{os.path.join(RAW, 'I114-win')} not found (--raw=DIR)")
            return
        e4 = {k: abs(float(D[k]) / sf[NAMES[k]] - 1) for k in VKEYS}
        report("CAP5c", all(v <= CAP_VS_SCIPY for v in e4.values()),
               f"{head}: E4 = crate D against #114's scipy D (its rtol 1e-12 evaluation finished; I114-win partial file), per variant: "
               + ", ".join(f"{k} {v:.1e}" for k, v in e4.items()) + f" (bar {CAP_VS_SCIPY:.0e} each); scipy's rtol 1e-12 value {sf[NAMES['tol']]:.9e} "
               f"against its base {sf[NAMES['base']]:.9e}")


def check_CAP5s():
    """CAP5 (c) again with the crate's star started as scipy's: rs = 1e-3 with the leading-order boundary condition (#186's S1 matched start)"""
    g = CAPS[5]
    sf = i114_l5(g)
    if sf is None:
        report("CAP5s", None, f"{os.path.join(RAW, 'I114-win')} not found (--raw=DIR)")
        return
    setup, got = get_rows(5, [g], "dd", VKEYS, tag="_rs1e-3", rs=1e-3)
    _, ex = get_rows(5, [g], "dd", VKEYS)
    e4 = {k: abs(float(got[(g, k)]["d_f64"]) / sf[NAMES[k]] - 1) for k in VKEYS}
    shift = {k: abs(float(got[(g, k)]["d_f64"]) / float(ex[(g, k)]["d_f64"]) - 1) for k in VKEYS}
    report("CAP5s", all(v <= CAP_VS_SCIPY for v in e4.values()),
           f"l = 5, gamma = {g:g}, the crate's star shot from rs = 1e-3 with scipy's leading-order boundary condition: crate D against #114's scipy D, per variant: "
           + ", ".join(f"{k} {v:.1e}" for k, v in e4.items()) + f" (bar {CAP_VS_SCIPY:.0e} each); printed, not gated: the change of the crate's D from rs = 1e-8 to rs = 1e-3: "
           + ", ".join(f"{k} {v:.1e}" for k, v in shift.items()))


def check_QP():
    import mpmath as mp
    from fractions import Fraction as Fr
    import sympy
    mp.mp.dps = 120
    lines, ok = [], True
    for case in ("well", "near"):
        y = [[Fr(((3 * i + 5 * c + 1) % 7) - 3, 8) for i in range(6)] for c in range(3)]
        if case == "near":
            for i in range(6):
                y[2][i] = y[1][i] + Fr(((i + 2) % 3) - 1, 2 ** 27)
        c0 = [Fr(1, 2), Fr(-1, 4), Fr(3, 8)]
        s = [[sum(y[c][2 + i] * c0[c] for c in range(3)) for i in range(4)], [Fr(((2 * i + 1) % 5) - 2, 4) for i in range(4)]]
        M = sympy.Matrix(4, 5, lambda i, k: sympy.Rational(y[k][2 + i].numerator, y[k][2 + i].denominator) if k < 3
                         else -sympy.Rational(s[k - 3][i].numerator, s[k - 3][i].denominator))
        ns = M.nullspace()
        assert len(ns) == 1, "planted planes are degenerate"
        c = [Fr(int(ns[0][j].p), int(ns[0][j].q)) for j in range(3)]
        met = [sum(y[cc][r] * c[cc] for cc in range(3)) for r in (0, 1)]
        mm, rm, l = Fr(1, 64), Fr(193, 8), 3
        lam = Fr((l - 1) * (l + 2), 2)
        f = 1 - 2 * mm / rm
        T = [[(lam * rm * rm - 3 * lam * mm * rm - 3 * mm * mm) / ((rm - 2 * mm) * (lam * rm + 3 * mm)), rm],
             [(lam * (lam + 1) * rm * rm + 3 * lam * mm * rm + 6 * mm * mm) / (rm * rm * (lam * rm + 3 * mm)), f]]
        det = T[0][0] * T[1][1] - T[0][1] * T[1][0]
        z = (T[1][1] * met[0] - T[0][1] * met[1]) / det
        zp = (T[0][0] * met[1] - T[1][0] * met[0]) / det
        q = f * zp / z
        exact = mp.mpf(q.numerator) / mp.mpf(q.denominator)
        for lane in ("f64", "dd", "qd"):
            o = rust("qplant", whole=True, lane=lane, case=case, l=l)
            e = float(abs(mp.mpf(o["q_int"]) / exact - 1))
            bar = QP_BAR[(lane, case)]
            ok &= e <= bar
            lines.append(f"{case} {lane} {e:.1e} (bar {bar:.0e})")
    report("QP", ok, "planted q_int: dyadic planes with Y[2:] c0 in span S, exact q_int from sympy's rational null space and Fractions, "
                     "l = 3, r_m = 24.125, M = 1/64: relative error " + "; ".join(lines))


def check_RP():
    import numpy as np
    lines, ok = [], True
    for cut, R, g0 in RP_CASES:
        o = rust("rmplant", whole=True, cut=repr(cut), R=repr(R), g0=repr(g0))
        rr = np.linspace(g0, R, 20000)
        phi = 0.02 * np.exp(-0.65 * rr)
        i = int(np.argmax(phi < cut))
        same = o["i"] == i and float(o["rm"]) == float(rr[i])
        ok &= same
        lines.append(f"(cut {cut:g}, R {R:g}, start {g0:g}) crate i = {o['i']} r_m = {float(o['rm']):.17g}, numpy i = {i} r_m = {float(rr[i]):.17g}: {'equal' if same else 'DIFFER'}")
    report("RP", ok, "planted r_m, phi = 0.02 exp(-0.65 r) through the crate's grid rule against numpy's linspace and argmax: " + "; ".join(lines))


def check_RM():
    import background as B
    from polar_solver import PolarProblem
    st60 = scipy_star()
    st90 = B.solve_star(OMEGA, EPS, guess=st60, R=90.0)
    spec = {"base": (st60, {}), "tol": (st60, dict(rtol=1e-12)), "phi": (st60, dict(phi_cut=1e-12)), "r90": (st90, {}), "r0": (st60, dict(r0=3e-3))}
    setup, _ = get_rows(2, [], "dd", VKEYS)
    rm_c = rm_of(setup)
    lines, ok = [], True
    for k, (st, kw) in spec.items():
        rs = float(PolarProblem(st, 2, **kw).rm)
        same = rs == rm_c[k]
        ok &= same
        mg = [v for v in setup["variants"] if v["variant"] == k][0]
        lines.append(f"{k}: scipy {rs:.17g}, crate {rm_c[k]:.17g} ({'equal' if same else 'DIFFER'}; ln(phi/cut) at r_m {mg['ln_phi_over_cut_at_rm']:.1e}, one grid point before {mg['ln_phi_over_cut_before']:.1e})")
    report("RM", ok, "r_m of the five variants, the crate's own star against PolarProblem(...).rm on scipy's stars: " + "; ".join(lines))


ROW_CHECKS = {"QP": check_QP, "RP": check_RP, "RM": check_RM, "RA2": lambda: check_RA(2), "RX": check_RX,
              "RA4": lambda: check_RA(4), "RA5": lambda: check_RA(5), "CAP5": lambda: check_CAP(5), "CAP5s": check_CAP5s,
              "RA6": lambda: check_RA(6), "CAP6": lambda: check_CAP(6)}

# the controls of --rows mode: --rows --mutate=NAME runs the checks named here with one deliberate error; every one must fail
ROW_CONTROLS = {
    "qint_cofsign": ["QP", "RA2"],   # the cofactor signs of q_int's null vector not alternated
    "qint_f": ["QP", "RA2"],         # q_int without the factor f
    "rm_off1": ["RP", "RM", "RA2"],  # the r_m grid index one too low
    "rm_cut": ["RM", "RA2"],         # the phi_cut variant's cut at 1e-11
    "r0_grid": ["RM", "RA2"],        # the r0 variant's r_m grid starting at 1e-3
    "tol_loose": ["CAP5b"],           # the rtol-analogue variant at tol 1e-1 (1e-3 moved D by 3.9e-6 only: dev_33; its entries of CAP5c and CAP5s are read in the log: those two fail on real data already)
    "qstart": ["RA2"],               # #186's: the Riccati started at r_m + 1 (RX at gamma = 300 does not see it: dev_31)
    "moncrief": ["RA2", "RX", "QP"],  # #186's: a sign in the Moncrief map
    "jsign": ["CAP5a"],              # #186's: J[3][2] negated, which flips the sign of D
}


# ---------------------------------------------------------------------------------------------------------------
# #188: the crate over complex s (python taylor_polar.py --complex)
#
#   python taylor_polar.py --complex [--only=R1,P1C] [--mutate=NAME] [--raw=DIR] [--bins=DIR] [--jsonl-dir=DIR] [--ref-dir=DIR] [--lane=dd]
#   python taylor_polar.py --complex --build                  cargo build --release --features complex --target-dir target/complex
#   python taylor_polar.py --complex --list                   the checks and the controls
#   python taylor_polar.py --complex --scipy-ref-make=BOX --rtol=1e-13 [--points=tc] --out=FILE
#                                                             winding2's D (its normalization) by the scipy path at stored points of BOX
#                                                             (b1: 147 points of the box l = 2, Re s in [0.1, 2]; b6: l = 6, Re s in [0.03, 0.1]),
#                                                             one JSON line per point, resumable (long: about 40 min for b1 at 1e-13)
#
# The counter is a copy of winding2.py's count() and --plant check (commit 0d1d376), taking F(s) from a cache of stored values
# (check W0, which shares nothing with the crate) or from the crate's `dserve` (checks CX2 and CXZ).
# ---------------------------------------------------------------------------------------------------------------
CBIN = os.path.join(CRATE, "target", "complex", "release", "taylor_polar")
BINS = os.path.join(CRATE, "target", "bins")        # old_default, new_default, old_libplane, new_libplane (R1; the raw README says how they are made)
REFDIR = None
LANE_EPS = {"f64": 1.1102230246251565e-16, "dd": 1.232595164407831e-32, "qd": 1.504632769052528e-64}
CCRIT = dict(
    P1C_PLANE=1e4, P1C_LOGVOL=1e4, P1C_ORTH=1e2,        # in units of the lane's epsilon, as multishoot's plane tests
    P1C_F=[7.0, -3.0, 5.0, -11.0], P1C_X=12.0,
    X1_BAR=1e-9, XR_BAR=1e-9,
    XR_RS=[200.0, 300.0, 600.0],
    CX2_WIND=0.05, CX2_STEP=0.45,
    CXZ_BAR=1e4,
    R1_D=[("f64", 2, [0.001, 1.0, 40.0]), ("dd", 2, [0.001, 1.0, 40.0]), ("qd", 2, [0.001, 1.0]),
          ("f64", 6, [0.001, 1.0, 40.0]), ("dd", 6, [0.001, 1.0, 40.0]), ("qd", 6, [0.001, 1.0])],
    CXZ_L=2.0, CXZ_LZ=2.0, CXZ_S0=(0.75, 0.5),
)
TIMEKEYS = {"t_y", "t_q", "t_s", "wall", "t_star", "t_table", "seconds"}
CRAWDIR = None   # set in complex_main
CBOXES = {
    # the box of RESULTS.md section 5, l = 2 (#126's w_p7603_l2): winding 0, planted pair 0.3 + 1.76i
    "b1": dict(l=2, x0=0.1, x1=2.0, y1=2.5, y0=0.0, plant=(0.3, 1.76), guard="all", expect=(0.0, 2.0, 0.39, 0.40),
               cache="I126-rpi5/cache/winding2_cache_eps-1_w0.7603_l2_0.1_2.0_2.5.pkl",
               tc=["0.1+0.24j", "0.1+1.76j", "0.1+2.5j", "2+0j", "2+2.5j", "1+1j", "0.1+0j", "0.5+0.5j"]),
    # #117's box b25_l2, Re s in [2, 5], |Im s| <= 10 (large |s|, small R_s): winding 0, planted 2.5 + 1i
    "b2": dict(l=2, x0=2.0, x1=5.0, y1=10.0, y0=0.0, plant=(2.5, 1.0), guard="all", expect=(0.0, 2.0, 0.39, 0.40),
               cache="I117-rpi5/cache/winding2_cache_eps-1_w0.7603_l2_2.0_5.0_10.0.pkl", tc=[]),
    # #115's strip at l = 6, Re s in [0.03, 0.1] (R_s = 1293): values only (X1b)
    "b6": dict(l=6, x0=0.03, x1=0.1, y1=2.5, y0=0.0, plant=(0.04, 0.25), guard=1, expect=None,
               cache="I115-rpi5/cache/winding2_cache_eps-1_w0.7603_l6_0.03_0.1_2.5.pkl",
               tc=["0.03+0.24j", "0.03+1.76j", "0.03+2.5j", "0.1+0j", "0.1+2.5j", "0.05+1j", "0.03+0j", "0.07+0.5j"]),
}


def cbox_cache(name):
    import pickle
    return pickle.load(open(os.path.join(CRAWDIR, CBOXES[name]["cache"]), "rb"))


def cbox_kmin(bx):
    """min Re k over the contour, as winding2.py computes it (200 points per edge, the edges of the box)"""
    import numpy as np
    ks = lambda s: (np.sqrt(1 + (s - 1j * OMEGA) ** 2), np.sqrt(1 + (s + 1j * OMEGA) ** 2))
    x0, x1, y0, y1 = bx["x0"], bx["x1"], bx["y0"], bx["y1"]
    cands = [x0 + 1j * y for y in np.linspace(y0, y1, 200)] + [x + 1j * y1 for x in np.linspace(x0, x1, 200)]
    if y0 != 0.0:
        cands += [x1 + 1j * y for y in np.linspace(y0, y1, 200)] + [x + 1j * y0 for x in np.linspace(x0, x1, 200)]
    return float(min(min(k.real for k in ks(c)) for c in cands))


def cbox_rs(name, rsmult=40.0):
    return R_M["base"] + rsmult / cbox_kmin(CBOXES[name])


def cstrip(s, rs):
    """winding2.py's removal of the trivial phases: exp(-i Im[(k+ + k-)(R_s - r_m) + (k+ + k- + s) r_m])"""
    import numpy as np
    kp, km = np.sqrt(1 + (s - 1j * OMEGA) ** 2), np.sqrt(1 + (s + 1j * OMEGA) ** 2)
    f = -1.0 if MUTATE == "strip_sign" else 1.0
    return complex(np.exp(-1j * f * ((kp + km) * (rs - R_M["base"]) + (kp + km + s) * R_M["base"]).imag))


def ckey(s):
    return complex(round(s.real, 12), round(s.imag, 12))


def tc_keys(name, cache):
    ks = []
    for t in CBOXES[name]["tc"]:
        k = min(cache, key=lambda c: abs(c - complex(t)))
        if k not in ks:
            ks.append(k)
    return ks


class Server:
    """a crate process that reads `re im` lines and answers one JSON line each (dserve, czserve); evaluations are kept in memory and,
    with a jsonl file, on disk so that a stopped run resumes"""

    def __init__(self, cmd, jsonl=None, binpath=None, **kw):
        args = [binpath or CBIN, cmd] + [f"{k}={v}" for k, v in kw.items()]
        if MUTATE:
            args.append(f"mutate={MUTATE}")
        env = dict(os.environ, RAYON_NUM_THREADS="2")
        if MUTATE:
            env["TP_MAX_SECONDS"] = "120"
        self.p = subprocess.Popen(["nice", "-n", "10"] + args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, env=env, cwd=HERE)
        self.setup = json.loads(self.p.stdout.readline())
        self.mem, self.jsonl, self.n = {}, jsonl, 0
        if jsonl and os.path.exists(jsonl):
            for line in open(jsonl):
                o = json.loads(line)
                self.mem[(o["re"], o["im"])] = o

    def eval(self, s):
        k = (s.real, s.imag)
        if k not in self.mem:
            self.p.stdin.write(f"{s.real!r} {s.imag!r}\n")
            self.p.stdin.flush()
            line = self.p.stdout.readline()
            if not line:
                raise RuntimeError("the crate process ended")
            o = json.loads(line)
            self.mem[k] = o
            self.n += 1
            if self.jsonl:
                with open(self.jsonl, "a") as fp:
                    fp.write(line)
        return self.mem[k]

    def close(self):
        try:
            self.p.stdin.close()
            self.p.wait(timeout=30)
        except Exception:
            self.p.kill()


def cjsonl(tag):
    if not JDIR:
        return None
    # the counter's own mutations leave the crate's values as they are, so they reuse the unmutated file
    return os.path.join(JDIR, tag + (f"__{MUTATE}" if MUTATE and MUTATE not in COUNTER_MUTS else "") + ".jsonl")


def count_winding(F, x0, x1, y1, y0=0.0, nseg=10, guard=1, plant=None, mutate=None):
    """winding2.py's count() and its --plant check as a function of F: returns the winding, the largest phase step, min |D|,
    and with a plant the count with the factor, the difference and the verdict"""
    import numpy as np
    sym = (y0 == 0.0)
    guard = np.inf if guard == "all" else int(guard)
    if sym:
        corners = [complex(x1, 0), complex(x1, y1), complex(x0, y1), complex(x0, 0)]
    else:
        corners = [complex(x1, y0), complex(x1, y1), complex(x0, y1), complex(x0, y0), complex(x1, y0)]
    if mutate == "reverse":
        corners = corners[::-1]

    margin = [np.inf]   # the smallest distance of any quantity the counter compares with a threshold (0.4, 0.4, 0.1) from that threshold

    def count(Fn):
        total, maxstep, nmin = 0.0, 0.0, np.inf
        for a, b in zip(corners[:-1], corners[1:]):
            def step(t0, t1, depth):
                nonlocal maxstep, nmin
                da, db = Fn(a + (b - a) * t0), Fn(a + (b - a) * t1)
                dph = np.angle(db / da)
                fine = t1 - t0 <= 1e-5
                if not fine:
                    margin[0] = min(margin[0], abs(abs(dph) - 0.4))
                if abs(dph) <= 0.4 and depth <= guard and not fine:
                    dm = Fn(a + (b - a) * 0.5 * (t0 + t1))
                    h1, h2 = np.angle(dm / da), np.angle(db / dm)
                    margin[0] = min(margin[0], abs(max(abs(h1), abs(h2)) - 0.4), abs(abs(h1 + h2 - dph) - 0.1))
                    if max(abs(h1), abs(h2)) > 0.4 or abs(h1 + h2 - dph) > 0.1:
                        return step(t0, 0.5 * (t0 + t1), depth + 1) + step(0.5 * (t0 + t1), t1, depth + 1)
                if abs(dph) > 0.4 and not fine:
                    return step(t0, 0.5 * (t0 + t1), depth + 1) + step(0.5 * (t0 + t1), t1, depth + 1)
                maxstep = max(maxstep, abs(dph))
                nmin = min(nmin, abs(da))
                return dph
            for t0, t1 in zip(np.linspace(0, 1, nseg + 1)[:-1], np.linspace(0, 1, nseg + 1)[1:]):
                total += step(t0, t1, 0)
        doubled = sym != (mutate == "wrongdouble")
        return (2 if doubled else 1) * total / (2 * np.pi), maxstep, nmin

    N, maxstep, nmin = count(F)
    out = dict(winding=float(N), maxstep=float(maxstep), nmin=float(nmin), margin=float(margin[0]))
    if plant is not None:
        s0 = complex(*plant)
        inside = x0 < s0.real < x1 and ((0 < abs(s0.imag) < y1) if sym else (y0 < s0.imag < y1))
        expected = (2 if sym else 1) * inside
        s0a = s0 + (10 * (x1 - x0) if mutate == "misplant" else 0)
        pf = (lambda s: (s - s0a) * (s - np.conj(s0a))) if (sym and mutate != "onesided") else (lambda s: (s - s0a))
        Np, maxstep_p, _ = count(lambda s: F(s) * pf(s))
        out.update(plant_count=float(Np), plant_diff=float(Np - N), plant_expected=expected, plant_maxstep=float(maxstep_p), plant_margin=float(margin[0]),
                   plant_ok=bool(abs((Np - N) - expected) < 0.05 and maxstep_p < 0.45))
    return out


COUNTER_MUTS = ("reverse", "wrongdouble", "misplant", "onesided")


def _stored_F(cache, asked, missing):
    def F(s):
        s = ckey(s)
        asked.add(s)
        if s not in cache:
            missing.add(s)
            return complex(1.0)
        return cache[s]
    return F


def check_W0():
    """the counter port on the stored values alone: it must give each stored box's result and ask for no point the cache lacks"""
    parts, ok = [], True
    for name in ("b1", "b2"):
        bx = CBOXES[name]
        cache = cbox_cache(name)
        asked, missing = set(), set()
        r = count_winding(_stored_F(cache, asked, missing), bx["x0"], bx["x1"], bx["y1"], bx["y0"], guard=bx["guard"], plant=bx["plant"],
                          mutate=MUTATE if MUTATE in COUNTER_MUTS else None)
        w_exp, p_exp, ms_exp, msp_exp = bx["expect"]
        good = (abs(r["winding"] - w_exp) < CCRIT["CX2_WIND"] and r["plant_ok"] and not missing
                and round(r["maxstep"], 2) == ms_exp and round(r["plant_maxstep"], 2) == msp_exp)
        ok &= good
        parts.append(f"{name}: winding {r['winding']:+.4f} (stored {w_exp:+.4f}), max step {r['maxstep']:.2f} (stored {ms_exp}), plant difference {r['plant_diff']:+.4f} "
                     f"(expected {r['plant_expected']}, step {r['plant_maxstep']:.2f}, stored {msp_exp}), asked {len(asked)} of {len(cache)} stored points, {len(missing)} missing; "
                     f"smallest distance of a compared quantity from its threshold (0.4, 0.4, 0.1 rad): {r['margin']:.2e} without the plant, {r['plant_margin']:.2e} with it, "
                     f"and the plant's largest step is {CCRIT['CX2_STEP'] - r['plant_maxstep']:.2f} rad under 0.45")
    report("W0", ok, "winding2's counter ported, on the stored values alone: " + "; ".join(parts))


def _crate_F(name, lane, jsonl=None, rsmult=40.0):
    """F(s) = winding2's D at s from the crate (dserve, one star and one J table), with the phase strip; F.walls has the wall time of each"""
    bx = CBOXES[name]
    rs = cbox_rs(name, rsmult)
    rs_int = rs + (1.0 if MUTATE == "rs_shift" else 0.0)    # mutation rs_shift: the scalar plane started at R_s + 1, the strip at R_s
    srv = Server("dserve", jsonl, lane=lane, l=bx["l"], rsfix=repr(rs_int), **star_kw())
    walls = {}

    def F(s):
        s = ckey(s)
        o = srv.eval(s.conjugate() if MUTATE == "conj_s" else s)    # mutation conj_s: D evaluated at conj(s)
        walls[s] = o["wall"]
        d = complex(o["d_f64"][0], o["d_f64"][1])
        return d * (1.0 if MUTATE == "strip_off" else cstrip(s, rs))   # mutation strip_off: winding2's phase strip left out
    F.walls = walls
    return F, srv


def check_CX2():
    """the counter on the crate's values: the stored boxes' counts, with the planted pair, as winding2 gave them"""
    lane = LANE_C
    parts, ok = [], True
    for name in ("b1", "b2"):
        bx = CBOXES[name]
        F, srv = _crate_F(name, lane, cjsonl(f"cx2_{name}_{lane}"))
        t = time.time()
        try:
            stored = cbox_cache(name)
            asked = set()
            Fa = lambda s, F=F, asked=asked: (asked.add(ckey(s)), F(s))[1]
            r = count_winding(Fa, bx["x0"], bx["x1"], bx["y1"], bx["y0"], guard=bx["guard"], plant=bx["plant"],
                              mutate=MUTATE if MUTATE in COUNTER_MUTS else None)
        finally:
            srv.close()
        w_exp = bx["expect"][0]
        good = abs(r["winding"] - w_exp) < CCRIT["CX2_WIND"] and r["plant_ok"] and r["maxstep"] < CCRIT["CX2_STEP"] and r["plant_maxstep"] < CCRIT["CX2_STEP"]
        ok &= good
        parts.append(f"{name} ({lane}): winding {r['winding']:+.4f} (stored {w_exp:+.4f}), max step {r['maxstep']:.2f} rad, min|D| {r['nmin']:.2e}, plant difference "
                     f"{r['plant_diff']:+.4f} (expected {r['plant_expected']}, step {r['plant_maxstep']:.2f}); {len(asked)} points asked, {len(asked - set(stored))} not in the "
                     f"stored cache, {len(set(stored) - asked)} stored points not asked; {srv.n} evaluated now, {time.time() - t:.0f} s")
    report("CX2", ok, "winding2's counter with the crate as the evaluator: " + "; ".join(parts))


def _x1_points(name, cache):
    if MUTATE:   # a control runs on the stratified points only
        return tc_keys(name, cache)
    return sorted(cache, key=lambda c: (c.real, c.imag))


def _load_ref(name, rtol):
    path = os.path.join(REFDIR or "", f"scipy_{name}_rtol{rtol:g}.jsonl")
    if not os.path.exists(path):
        return None
    d = {}
    for line in open(path):
        o = json.loads(line)
        d[complex(o["key"][0], o["key"][1])] = o
    return d


def _check_X(label, name, lane):
    import numpy as np
    bx = CBOXES[name]
    cache = cbox_cache(name)
    ref = _load_ref(name, 1e-13)
    if ref is None:
        report(label, None, f"no scipy reference file scipy_{name}_rtol1e-13.jsonl in --ref-dir")
        return
    keys = [k for k in _x1_points(name, cache) if k in ref]
    F, srv = _crate_F(name, lane, cjsonl(f"{label.lower()}_{name}_{lane}"))
    try:
        rows = []
        for k in keys:
            d = F(k)
            r = complex(ref[k]["d"][0], ref[k]["d"][1])
            st = cache[k]
            rows.append((k, abs(abs(d) / abs(r) - 1), abs(np.angle(d / r)), abs(abs(d) / abs(st) - 1), abs(np.angle(d / st)),
                         abs(abs(r) / abs(st) - 1), abs(np.angle(r / st)), F.walls[ckey(k)]))
    finally:
        srv.close()
    bar = CCRIT["X1_BAR"]
    mod, ph = max(r[1] for r in rows), max(r[2] for r in rows)
    nbad = sum(1 for r in rows if r[1] > bar or r[2] > bar)
    ok = nbad == 0
    wall = sorted(r[7] for r in rows)
    report(label, ok, f"box {name} (l = {bx['l']}, Re s in [{bx['x0']}, {bx['x1']}], |Im s| <= {bx['y1']}, R_s {cbox_rs(name):.4f}), crate {lane} against winding2's D by the scipy path at rtol 1e-13 "
                      f"at {len(rows)} stored points: worst modulus ratio {mod:.1e}, worst phase {ph:.1e} rad (bar {bar:g} each), {nbad} points over the bar. Reported without a bar: crate against the stored values "
                      f"(rtol 1e-8) modulus {max(r[3] for r in rows):.1e}, phase {max(r[4] for r in rows):.1e}; the stored values against scipy at 1e-13 modulus {max(r[5] for r in rows):.1e}, phase {max(r[6] for r in rows):.1e}; "
                      f"crate wall per D: median {wall[len(wall) // 2]:.2f} s, max {wall[-1]:.2f} s")
    return rows


def check_X1():
    _check_X("X1", "b1", LANE_C)
    if not MUTATE:
        try:
            rows = _check_X("X1-f64", "b1", "f64")   # the f64 lane against the same reference, reported: counts as a check of its own
        except Exception as e:
            report("X1-f64", False, f"raised {type(e).__name__}: {str(e)[:200]}")


def check_X1B():
    _check_X("X1b", "b6", LANE_C)


def check_XR():
    """winding2's D with the phase strip is independent of the scalar plane's start radius: the crate at four R_s, eight points"""
    import numpy as np
    name = "b1"
    cache = cbox_cache(name)
    keys = tc_keys(name, cache)
    rs0 = cbox_rs(name)
    rss = [rs0] + CCRIT["XR_RS"]
    vals = {}
    for rs in rss:
        srv = Server("dserve", None, lane=LANE_C, l=CBOXES[name]["l"], rsfix=repr(rs), **star_kw())
        try:
            vals[rs] = [complex(*srv.eval(ckey(k))["d_f64"]) * cstrip(ckey(k), rs) for k in keys]
        finally:
            srv.close()
    worst_m = worst_p = 0.0
    for rs in rss[1:]:
        for a, b in zip(vals[rs], vals[rs0]):
            worst_m = max(worst_m, abs(abs(a) / abs(b) - 1))
            worst_p = max(worst_p, abs(np.angle(a / b)))
    ok = worst_m <= CCRIT["XR_BAR"] and worst_p <= CCRIT["XR_BAR"]
    report("XR", ok, f"stripped D at R_s = {', '.join(f'{r:.2f}' for r in CCRIT['XR_RS'])} against R_s = {rs0:.2f} ({LANE_C}, {len(keys)} points of b1): worst modulus ratio {worst_m:.1e}, "
                     f"worst phase {worst_p:.1e} rad (bar {CCRIT['XR_BAR']:g} each)"
                     + ("" if ok else "; the premise is false: the phase carries tails of the exterior solution that winding2's strip does not remove (1/R_s from the "
                        "start vectors and the angular term, ln R_s from the mass); they are holomorphic and zero-free at a fixed R_s, so a count is unaffected (XRc)"))


def check_XRC():
    """the counts do not depend on R_s: b1 counted with the scalar plane started at R_s for rsmult 40 and 80 (winding2's --rsmult), winding 0, planted pair 2"""
    bx = CBOXES["b1"]
    parts, ok = [], True
    for rsm in (40.0, 80.0):
        F, srv = _crate_F("b1", LANE_C, cjsonl(f"xrc_b1_{LANE_C}_rsm{rsm:g}"), rsmult=rsm)
        try:
            r = count_winding(F, bx["x0"], bx["x1"], bx["y1"], bx["y0"], guard=bx["guard"], plant=bx["plant"], mutate=MUTATE if MUTATE in COUNTER_MUTS else None)
        finally:
            srv.close()
        good = abs(r["winding"]) < CCRIT["CX2_WIND"] and r["plant_ok"] and r["maxstep"] < CCRIT["CX2_STEP"] and r["plant_maxstep"] < CCRIT["CX2_STEP"]
        ok &= good
        parts.append(f"rsmult {rsm:g} (R_s {cbox_rs('b1', rsm):.1f}): winding {r['winding']:+.4f}, max step {r['maxstep']:.2f}, plant difference {r['plant_diff']:+.4f} (step {r['plant_maxstep']:.2f})")
    report("XRc", ok, f"b1 counted with the scalar plane started at two outer radii ({LANE_C}): " + "; ".join(parts))


# ---------------------------------------------------------------------------------------------------------------
# R1, R2: the real lanes
# ---------------------------------------------------------------------------------------------------------------
def _collect_real(binpath, star):
    out = []
    for lane, l, gs in CCRIT["R1_D"]:
        out += crun("d", binpath, False, lane=lane, l=l, gammas=",".join(repr(g) for g in gs), **star)
    for case in ("const", "linear"):
        for lane in ("f64", "dd", "qd"):
            out += crun("planted", binpath, False, case=case, lane=lane, gamma=10.0, L=1.0, Lz=1.0)
    out += crun("star", binpath, False, **MODEL)
    return out


def crun(cmd, binpath=None, mutate=True, stdin=None, **kw):
    args = [binpath or CBIN, cmd] + [f"{k}={v}" for k, v in kw.items()]
    if MUTATE and mutate:
        args.append(f"mutate={MUTATE}")
    env = dict(os.environ, RAYON_NUM_THREADS="2")
    if MUTATE:
        env["TP_MAX_SECONDS"] = "120"
    p = subprocess.run(["nice", "-n", "10"] + args, capture_output=True, text=True, env=env, cwd=HERE, input=stdin)
    if p.returncode:
        raise RuntimeError(f"{' '.join(args)[:200]} exited {p.returncode}: {p.stderr[-400:]}")
    return [json.loads(x) for x in p.stdout.splitlines() if x.startswith("{")]


def check_R1():
    """the real lanes are unchanged: the binary of main 0d1d376 against the new build, default and lib_plane, every field but the times"""
    star = star_kw()
    parts, ok = [], True
    for old, new in (("old_default", "new_default"), ("old_libplane", "new_libplane")):
        bo, bn = os.path.join(BINS, old), os.path.join(BINS, new)
        if not (os.path.exists(bo) and os.path.exists(bn)):
            report("R1", None, f"{old} or {new} is missing in {BINS}")
            return
        lo, ln = _collect_real(bo, star), _collect_real(bn, star)
        if MUTATE == "r1_ulp":   # the last digit of one D in the new log changed
            d = next(o for o in ln if o["cmd"] == "d")
            d["d"] = d["d"][:-1] + ("1" if d["d"][-1] != "1" else "2")
        diff = [i for i, (a, b) in enumerate(zip(lo, ln)) if {k: v for k, v in a.items() if k not in TIMEKEYS} != {k: v for k, v in b.items() if k not in TIMEKEYS}]
        same = len(lo) == len(ln) and not diff
        ok &= same
        parts.append(f"{old} against {new}: {len(lo)} records, {len(diff)} differ" + ("" if same else f" (first: {lo[diff[0]]['cmd'] if diff else 'length'})"))
    report("R1", ok, "real lanes, every field except times as strings (d at l = 2, 6, gamma = 0.001, 1, 40 in f64, dd, qd; planted const and linear; star): " + "; ".join(parts))


def check_R2():
    """the complex lane at Im s = 0 is the real lane: D (re as a string, im = 0), steps, q, r_s, r_far; and P1c with f = 0"""
    star = star_kw()
    bo = os.path.join(BINS, "old_default")
    if not os.path.exists(bo):
        report("R2", None, f"{bo} is missing")
        return
    n, bad = 0, []
    for lane, l, gs in CCRIT["R1_D"]:
        real = [o for o in crun("d", bo, False, lane=lane, l=l, gammas=",".join(repr(g) for g in gs), **star) if o["cmd"] == "d"]
        cplx = [o for o in crun("dc", CBIN, True, lane=lane, l=l, re=",".join(repr(g) for g in gs), im=",".join("0.0" for _ in gs), **star) if o["cmd"] == "dc"]
        for a, b in zip(real, cplx):
            n += 1
            same = (a["d"] == b["d_re"] and float(b["d_im"]) == 0.0 and (a["steps_y"], a["steps_q"], a["steps_s"]) == (b["steps_y"], b["steps_q"], b["steps_s"])
                    and a["q"] == b["q_f64"][0] and b["q_f64"][1] == 0.0 and a["r_s"] == b["r_s"] and a["r_far"] == b["r_far"])
            if not same:
                bad.append(f"{lane} l={l} g={a['gamma']:g}")
    for lane in ("f64", "dd", "qd"):
        r = crun("cp1", CBIN, True, lane=lane, cx=0, f="0,0,0,0")
        c = crun("cp1", CBIN, True, lane=lane, cx=1, f="0,0,0,0")
        n += 1
        same = (r[0]["lnvol_re"] == c[0]["lnvol_re"] and float(c[0]["lnvol_im"]) == 0.0 and r[0]["steps"] == c[0]["steps"]
                and all(a["re"] == b["re"] and all(float(x) == 0.0 for x in b["im"]) for a, b in zip(r[1:], c[1:])))
        if not same:
            bad.append(f"P1c f = 0 {lane}")
    report("R2", not bad, f"complex lane at real s against the real lane, {n} cases (the d list of R1 against `dc`, and P1c with f = 0 on cx = 0 and 1): "
                          + ("all equal as strings, Im D = 0" if not bad else f"{len(bad)} differ: " + ", ".join(bad[:6])))


# ---------------------------------------------------------------------------------------------------------------
# P1c: the planted complex plane
# ---------------------------------------------------------------------------------------------------------------
def check_P1C():
    import mpmath as mp
    mp.mp.dps = 140
    P = [[1, 1, 0, 1], [1, 2, 1, 1], [0, 1, 2, 1], [1, 1, 1, 3]]
    G = [24, 8, -8, -24]
    f_true = CCRIT["P1C_F"]
    x = mp.mpf(CCRIT["P1C_X"])
    cols = [[mp.mpf(P[i][c]) for i in range(4)] for c in range(2)]
    q0, lnr = [], mp.mpf(0)
    for c in cols:   # positive-diagonal Gram-Schmidt of P[:, :2] at 140 digits (twice, as the stepper does)
        v = mp.matrix(c)
        for _ in range(2):
            for q in q0:
                v = v - mp.fsum(q[i] * v[i] for i in range(4)) * q
        nrm = mp.sqrt(mp.fsum(t * t for t in v))
        q0.append(v / nrm)
        lnr += mp.log(nrm)
    lam_re = G[0] + G[1]
    lam_im = f_true[0] + f_true[1]
    lam_exact = mp.mpc(lam_re * x + lnr, lam_im * x)   # Lambda = (g1 + g2) x + ln det R0 + i (f1 + f2) x
    f_run = list(f_true)
    if MUTATE == "f_scaled":
        f_run[1] = f_run[1] * (1 + 1e-6)
    msgs, ok = [], True
    for lane in ("f64", "dd", "qd"):
        eps = LANE_EPS[lane]
        out = crun("cp1", CBIN, True, lane=lane, cx=1, x=repr(CCRIT["P1C_X"]), f=",".join(repr(t) for t in f_run))
        head, qs = out[0], out[1:]
        Y = [[mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b in zip(o["re"], o["im"])] for o in qs]
        # plane: ||(I - Pi) Y||_F with Pi the projector on span P[:, :2] (exact, from the 140-digit basis)
        err = mp.mpf(0)
        for y in Y:
            proj = [mp.mpc(0)] * 4
            for q in q0:
                cq = mp.fsum(q[i] * y[i] for i in range(4))
                proj = [proj[i] + cq * q[i] for i in range(4)]
            err += mp.fsum(abs(y[i] - proj[i]) ** 2 for i in range(4))
        plane = mp.sqrt(err)
        orth = max(abs(mp.fsum(mp.conj(Y[a][i]) * Y[b][i] for i in range(4)) - (1 if a == b else 0)) for a in range(2) for b in range(2))
        # Lambda = ln_vol + ln det(Q0^T Q)
        m = mp.matrix(2, 2)
        for a in range(2):
            for b in range(2):
                m[a, b] = mp.fsum(q0[a][i] * Y[b][i] for i in range(4))
        lam_got = mp.mpf(head["lnvol_re"]) + mp.log(mp.det(m))
        dre = abs(lam_got.real - lam_exact.real)
        dim = abs(mp.mpf((lam_got.imag - lam_exact.imag + mp.pi) % (2 * mp.pi) - mp.pi))
        dlam = max(dre, dim)
        bars = (CCRIT["P1C_PLANE"] * eps, CCRIT["P1C_LOGVOL"] * eps * max(1.0, float(abs(lam_exact))), CCRIT["P1C_ORTH"] * eps)
        good = float(plane) <= bars[0] and float(dlam) <= bars[1] and float(orth) <= bars[2]
        ok &= good
        msgs.append(f"{lane}: plane {fmt(plane)} (bar {fmt(bars[0])}), Lambda {fmt(dlam)} (re {fmt(dre)}, im {fmt(dim)}; bar {fmt(bars[1])}), orthonormality {fmt(orth)} (bar {fmt(bars[2])}), {head['steps']} steps")
    report("P1C", ok, f"J = P diag(24 + 7i, 8 - 3i, -8 + 5i, -24 - 11i) P^-1, plane P[:, :2], x = 12; exact Lambda = ln_vol + ln det(Q0^T Q) = {mp.nstr(lam_exact, 12)} at 140 digits; " + "; ".join(msgs))


# ---------------------------------------------------------------------------------------------------------------
# CXZ: a planted determinant with zeros
# ---------------------------------------------------------------------------------------------------------------
def _cz_J(mp):
    C = [2.0, -3.0, 1.0, 4.0, -2.0, 3.0]
    J = mp.matrix(6, 6)
    for i in range(6):
        for j in range(6):
            re = (((7 * i + 3 * j + 1) % 11) - 5.0) / 8.0
            im = (((3 * i + 5 * j + 2) % 9) - 4.0) / 4.0
            if i == j:
                im += C[i]
            J[i, j] = mp.mpc(re, im)
    return J


def check_CXZ():
    """a planted determinant: D(s) = c p(s), p(s) = (s - s0)(s - conj s0), through the complex stepper; D against mpmath's expm at 20 contour points
    and the count of the box Re s in [0.2, 1.2], |Im s| <= 1 (2 zeros) and of Re s in [1.3, 2.0] (none)"""
    import mpmath as mp
    import numpy as np
    mp.mp.dps = 140
    L, Lz, s0 = CCRIT["CXZ_L"], CCRIT["CXZ_LZ"], CCRIT["CXZ_S0"]
    J = _cz_J(mp)
    z0 = mp.matrix([mp.mpf(v) for v in Z0[0]])
    zt = mp.expm(-J * mp.mpf(L + Lz)) * z0
    lane = LANE_C
    eps = LANE_EPS[lane]
    srv = Server("czserve", None, lane=lane, L=L, Lz=Lz, s0re=s0[0], s0im=s0[1],
                 zt_re=",".join(mp.nstr(zt[i].real, 135, strip_zeros=False) for i in range(6)), zt_im=",".join(mp.nstr(zt[i].imag, 135, strip_zeros=False) for i in range(6)))
    try:
        # J as the crate has it
        jr = srv.setup["j_re"]
        jm = srv.setup["j_im"]
        same_j = all(abs(J[i, j].real - jr[6 * i + j]) == 0 and abs(J[i, j].imag - jm[6 * i + j]) == 0 for i in range(6) for j in range(6))
        eL, eZ = mp.expm(J * mp.mpf(L)), mp.expm(-J * mp.mpf(Lz))
        s0c = mp.mpc(*s0)

        def exact(s):
            sc = mp.mpc(s.real, s.imag)
            p = (sc - s0c) * (sc - mp.conj(s0c))
            Y = mp.matrix(6, 3)
            for c in range(2):
                for i in range(6):
                    Y[i, c] = mp.mpf(Y0[c][i])
            for i in range(6):
                Y[i, 2] = zt[i] + p * mp.mpf(Y0[2][i])
            Zc = mp.matrix(6, 3)
            for c in range(3):
                for i in range(6):
                    Zc[i, c] = mp.mpf(Z0[c][i])
            YL, ZL = eL * Y, eZ * Zc
            M = mp.matrix(6, 6)
            for c in range(3):
                for i in range(6):
                    M[i, c], M[i, 3 + c] = YL[i, c], ZL[i, c]
            return mp.det(M)

        worst = 0.0
        pts = [complex(1.2, 0), complex(1.2, 0.5), complex(1.2, 1.0)] + [complex(1.2 - 1.0 * t, 1.0) for t in np.linspace(0.05, 1.0, 10)] + [complex(0.2, 1.0 - t) for t in np.linspace(0.1, 1.0, 7)]
        for s in pts:
            o = srv.eval(s)
            d_got = mp.mpc(mp.mpf(o["d_re"]), mp.mpf(o["d_im"])) * mp.exp(mp.mpf(o["ln_y"]) + mp.mpf(o["ln_z"]))
            d_ex = exact(s)
            scale = max(1.0, float(abs(mp.mpf(o["ln_y"])) + abs(mp.mpf(o["ln_z"]))))
            worst = max(worst, float(abs(d_got / d_ex - 1)) / (eps * scale))
        # the count
        res = []
        for (x0, x1, expect) in ((0.2, 1.2, 2.0), (1.3, 2.0, 0.0)):
            F = lambda s: complex(*srv.eval(ckey(s))["d_f64"])
            r = count_winding(F, x0, x1, 1.0, 0.0, mutate=MUTATE if MUTATE in COUNTER_MUTS else None)
            res.append((x0, x1, expect, r))
    finally:
        srv.close()
    bar = CCRIT["CXZ_BAR"]
    ok = same_j and worst <= bar and all(abs(r["winding"] - e) < CCRIT["CX2_WIND"] and r["maxstep"] < CCRIT["CX2_STEP"] for (_, _, e, r) in res)
    report("CXZ", ok, f"planted determinant, {lane}, L = {L}, Lz = {Lz}, zeros at {s0[0]} +- {s0[1]}i; J equal to the crate's: {same_j}; D against mpmath expm at {len(pts)} contour points: "
                      f"worst |D/D_exact - 1| / (eps max(1, ln growth)) = {worst:.2e} (bar {bar:g}); " + "; ".join(
                          f"box Re s in [{x0}, {x1}]: winding {r['winding']:+.4f} (planted {e:.0f}), max step {r['maxstep']:.2f}" for (x0, x1, e, r) in res))


def check_L1C():
    """the lane: the crate's D at three points of b1 in qd against dd (bar 1e-20 in modulus ratio and phase, as #186's L1), f64 against dd reported"""
    import numpy as np
    name = "b1"
    cache = cbox_cache(name)
    keys = [min(cache, key=lambda c, t=t: abs(c - complex(t))) for t in ("0.1+1.76j", "2+2.5j", "0.5+0.5j")]
    rs = cbox_rs(name)
    star = star_kw()
    vals = {}
    for lane in ("f64", "dd", "qd"):
        out = crun("dc", CBIN, True, lane=lane, l=CBOXES[name]["l"], rsfix=repr(rs), re=",".join(repr(k.real) for k in keys), im=",".join(repr(k.imag) for k in keys), **star)
        vals[lane] = [o for o in out if o["cmd"] == "dc"]
    import mpmath as mp
    mp.mp.dps = 80
    def dm(a, b):
        za, zb = mp.mpc(mp.mpf(a["d_re"]), mp.mpf(a["d_im"])), mp.mpc(mp.mpf(b["d_re"]), mp.mpf(b["d_im"]))
        return float(abs(abs(za) / abs(zb) - 1)), float(abs(mp.arg(za / zb)))
    q = [dm(a, b) for a, b in zip(vals["qd"], vals["dd"])]
    f = [dm(a, b) for a, b in zip(vals["f64"], vals["dd"])]
    bar = 1e-20
    ok = all(m <= bar and p <= bar for m, p in q)
    report("L1C", ok, f"crate D at {len(keys)} points of b1: qd against dd worst modulus ratio {max(m for m, _ in q):.1e}, phase {max(p for _, p in q):.1e} (bar {bar:g}); "
                      f"f64 against dd (reported) {max(m for m, _ in f):.1e}, {max(p for _, p in f):.1e}; wall per D dd {sum(o['wall'] for o in vals['dd']) / 3:.2f} s, qd {sum(o['wall'] for o in vals['qd']) / 3:.1f} s")


# ---------------------------------------------------------------------------------------------------------------
# TC: time per D
# ---------------------------------------------------------------------------------------------------------------
def check_TC():
    parts = []
    for name in ("b1", "b6"):
        bx = CBOXES[name]
        ref8 = _load_ref(name, 1e-8)
        ref13 = _load_ref(name, 1e-13)
        if ref8 is None:
            continue
        keys = [k for k in tc_keys(name, cbox_cache(name)) if k in ref8]
        for lane in ("f64", LANE_C) if LANE_C != "f64" else ("f64",):
            srv = Server("dserve", None, lane=lane, l=bx["l"], rsfix=repr(cbox_rs(name)), **star_kw())
            try:
                w = [srv.eval(ckey(k))["wall"] for k in keys]
            finally:
                srv.close()
            ts8 = [ref8[k]["t"] for k in keys]
            ts13 = [ref13[k]["t"] for k in keys] if ref13 else []
            parts.append(f"{name} (l = {bx['l']}, R_s {cbox_rs(name):.1f}, {len(keys)} points) crate {lane}: {sum(w) / len(w):.3f} s per D (range {min(w):.3f} to {max(w):.3f}); scipy rtol 1e-8: "
                         f"{sum(ts8) / len(ts8):.2f} s (range {min(ts8):.2f} to {max(ts8):.2f})" + (f", rtol 1e-13: {sum(ts13) / len(ts13):.2f} s" if ts13 else "")
                         + f"; speed-up over rtol 1e-8 x{sum(ts8) / sum(w):.0f}")
    report("TC", True if parts else None, "time per D, one process each, nice 10, load " + f"{os.getloadavg()[0]:.1f}: " + ("; ".join(parts) if parts else "no scipy rtol 1e-8 files in --ref-dir"))


# ---------------------------------------------------------------------------------------------------------------
# the scipy reference
# ---------------------------------------------------------------------------------------------------------------
def scipy_ref_make(name, out_path, rtol, which):
    import pickle
    import numpy as np
    import background as B
    from polar_solver import PolarProblem
    bx = CBOXES[name]
    cache = cbox_cache(name)
    keys = tc_keys(name, cache) if which == "tc" else sorted(cache, key=lambda c: (c.real, c.imag))
    st = B.solve_star(OMEGA, EPS)
    P = PolarProblem(st, bx["l"], rtol=rtol)
    assert float(P.rm) == R_M["base"], (P.rm, R_M["base"])
    rs = cbox_rs(name)
    P.Rs_fixed = rs
    ks = lambda s: (np.sqrt(1 + (s - 1j * OMEGA) ** 2), np.sqrt(1 + (s + 1j * OMEGA) ** 2))
    done = set()
    if os.path.exists(out_path):
        for line in open(out_path):
            o = json.loads(line)
            done.add(complex(o["key"][0], o["key"][1]))
    print(f"# scipy path, box {name}, l = {bx['l']}, rtol {rtol:g}, R_s fixed {rs!r}, {len(keys)} points ({len(done)} done), load {os.getloadavg()[0]:.2f}", flush=True)
    for k in keys:
        if k in done:
            continue
        kp, km = ks(k)
        t0, c0 = time.time(), time.process_time()
        d = P.D(k) * np.exp(-1j * ((kp + km) * (P.Rs_fixed - P.rm) + (kp + km + k) * P.rm).imag)
        rec = dict(key=[k.real, k.imag], d=[float(d.real), float(d.imag)], t=time.time() - t0, cpu=time.process_time() - c0, load1=os.getloadavg()[0], rtol=rtol, rs=rs, box=name)
        with open(out_path, "a") as fp:
            fp.write(json.dumps(rec) + "\n")
        print(f"{k} {d} {rec['t']:.1f} s", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# the mode
# ---------------------------------------------------------------------------------------------------------------
LANE_C = "dd"
C_CHECKS = {"R1": check_R1, "R2": check_R2, "P1C": check_P1C, "XR": check_XR, "XRc": check_XRC, "L1C": check_L1C, "X1": check_X1, "X1b": check_X1B, "W0": check_W0, "CX2": check_CX2, "CXZ": check_CXZ, "TC": check_TC}
C_CONTROLS = {
    "r1_ulp": ["R1"],                   # the last digit of one D of the new build's log changed
    "cx_recip_general": ["R2"],         # Cx's reciprocal and quotient by the general formula at im = 0 (the real lanes' bits lost)
    "bilinear": ["P1C"],                # Gram-Schmidt without the conjugate (and the squared norm without |.|^2)
    "redot": ["P1C"],                   # the imaginary part of the projector dropped from the inner product
    "f_scaled": ["P1C"],                # f2 scaled by 1 + 1e-6 in the system, not in the reference
    "cx_mul_sign": ["P1C", "X1"],       # re = ac + bd in the complex product
    "strip_off": ["X1", "X1b"],                # winding2's phase strip left out
    "conj_s": ["X1", "X1b", "CXZ"],            # D evaluated at conj(s)
    "rs_shift": ["X1", "X1b"],                 # the scalar plane started at R_s + 1, the strip at R_s
    "lowsum": ["L1C"],                  # #186's: the Taylor sum cut at order 8 (the two lanes err differently)
    "strip_sign": ["X1"],               # exp(+i ...) in the strip (XR cannot serve: it fails on the unmutated crate)
    "jsign": ["X1"], "moncrief": ["X1"], "swapcols": ["X1"], "mext0": ["X1"], "star_eps": ["X1"], "shift_off": ["X1"],   # #186's controls on D
    "reverse": ["W0", "CX2", "XRc", "CXZ"], "wrongdouble": ["W0", "CX2", "XRc", "CXZ"],     # the counter's own (winding2's)
    "misplant": ["W0", "CX2", "XRc"], "onesided": ["W0", "CX2", "XRc"],
}


def complex_main(args):
    global MUTATE, JDIR, RAW, CRAWDIR, REFDIR, BINS, LANE_C, SCIPY_FILE
    if "--build" in args:
        env = dict(os.environ, CARGO_NET_GIT_FETCH_WITH_CLI="true")
        return subprocess.check_call(["nice", "-n", "10", "cargo", "build", "--release", "-j", "2", "--features", "complex", "--target-dir", "target/complex"], cwd=CRATE, env=env)
    if "--list" in args:
        print("--complex checks:", " ".join(C_CHECKS))
        print("--complex controls:", " ".join(C_CONTROLS))
        return
    only, ref_make, out, rtol, which = None, None, None, 1e-13, "all"
    CRAWDIR = RAW
    for a in args:
        if a.startswith("--mutate="):
            MUTATE = a.split("=", 1)[1]
        if a.startswith("--only="):
            only = [x for x in a.split("=", 1)[1].split(",") if x]
        if a.startswith("--raw="):
            CRAWDIR = a.split("=", 1)[1]
        if a.startswith("--bins="):
            BINS = a.split("=", 1)[1]
        if a.startswith("--jsonl-dir="):
            JDIR = a.split("=", 1)[1]
        if a.startswith("--ref-dir="):
            REFDIR = a.split("=", 1)[1]
        if a.startswith("--lane="):
            LANE_C = a.split("=", 1)[1]
        if a.startswith("--scipy-ref-make="):
            ref_make = a.split("=", 1)[1]
        if a.startswith("--out="):
            out = a.split("=", 1)[1]
        if a.startswith("--rtol="):
            rtol = float(a.split("=", 1)[1])
        if a.startswith("--points="):
            which = a.split("=", 1)[1]
    if ref_make:
        return scipy_ref_make(ref_make, out, rtol, which)
    if not os.path.exists(CBIN):
        sys.exit("build first: python taylor_polar.py --complex --build")
    if not os.path.exists(BIN):
        build()
    if MUTATE and MUTATE not in C_CONTROLS:
        sys.exit(f"unknown control {MUTATE}; known: {' '.join(C_CONTROLS)}")
    names = only if only is not None else (C_CONTROLS[MUTATE] if MUTATE else list(C_CHECKS))
    names = list(dict.fromkeys(names))
    if JDIR:
        os.makedirs(JDIR, exist_ok=True)
    t0 = time.time()
    print(f"# taylor_polar.py --complex: start {time.strftime('%Y-%m-%d %H:%M:%S')}, load average {os.getloadavg()[0]:.2f}, one process, nice 10; lane {LANE_C}; mutate = {MUTATE}", flush=True)
    for n in names:
        try:
            C_CHECKS[n]()
        except Exception as e:
            report(n, False, f"raised {type(e).__name__}: {str(e)[:300]}")
    passed = [n for n, s in RESULTS if s is True]
    failed = [n for n, s in RESULTS if s is False]
    skipped = [n for n, s in RESULTS if s is None]
    print(f"SUMMARY: {len(passed)} passed, {len(failed)} failed, {len(skipped)} skipped" + (f" ({', '.join(skipped)}: not run)" if skipped else "")
          + f" of {len(RESULTS)} checks" + (f"; failed: {', '.join(failed)}" if failed else "") + f"; mutate = {MUTATE}; {time.time() - t0:.0f} s")
    if MUTATE:
        not_fired = [n for n in C_CONTROLS[MUTATE] if not any(r.startswith(n) and s is False for r, s in RESULTS)]
        print(f"CONTROL {MUTATE}: " + ("every targeted check failed, as it must" if not not_fired else f"NOT FIRED on {', '.join(not_fired)}"))
    sys.exit(1 if failed else 0)


# ---------------------------------------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------------------------------------
CHECKS = {
    "T1": check_T1, "T2": check_T2, "R1": check_R1,
    "P1": check_P1, "P0": check_P0, "P2": check_P2,
    "S1": check_S1, "S2": check_S2, "S3": check_S3, "T3": check_T3,
    "A2": lambda: check_A(2, CRIT["L2_GAMMAS"], "A2[l=2]"),
    "A6": lambda: check_A(6, CRIT["L6_GAMMAS"], "A6[l=6]"),
    "A4": lambda: check_A(4, CRIT["L4_GAMMAS"], "A4[l=4]"),
    "A5": lambda: check_A(5, CRIT["L5_GAMMAS"], "A5[l=5]"),
    "L1": check_L1, "V1": check_V1, "V2": check_V2, "V3": check_V3, "V4": check_V4, "TS": check_TS,
}


# the controls: --mutate=NAME runs the checks named here with one deliberate error; every one must fail
CONTROLS = {
    "jsign": ["T1", "T2", "A2"],          # J[3][2], the largest entry, negated (jcheck, jjets, the table)
    "lowsum": ["P0", "P1", "P2"],          # the Taylor sum cut at order 8, the step sized for order 30 (D moves 2e-6: A2 is not a target)
    "lin_drop": ["P2"],                   # planted matrix without its linear term (the reference keeps it)
    "recip_dd": ["R1", "P1"],             # 1/k formed in double-double and widened (quad- and octo-double lanes)
    "noqr": ["P1", "A2"],                 # no orthonormalization between steps (P1 at gamma = 3 only: the others never finish)
    "star_eps": ["S1", "S2", "A2"],       # the star solved with the sign of eps flipped
    "zc1": ["S1", "A2"],                  # the lapse normalization skipped (Z_c = 1)
    "start_term": ["S1", "S3"],           # the shooting's start u = r c / 2 instead of r c / 3
    "shift_off": ["T3", "A2"],            # the table's jets not shifted to the step start
    "wide_table": ["T3", "A2"],           # table steps as wide as the star's, not cut by J's convergence
    "moncrief": ["A2", "V4"],             # a sign in the Moncrief map
    "swapcols": ["A2"],                   # the scalar plane's two start columns swapped (orientation)
    "mext0": ["A2"],                      # the exterior mass set to zero
    "qstart": ["A2"],                     # the Riccati started at r_m + 1 instead of r_m + 40/s + 20 l/s + 50
}


def main():
    global MUTATE, JDIR, RAW
    args = sys.argv[1:]
    if "--complex" in args:
        return complex_main(args)
    if "--gen" in args:
        return gen()
    if "--build" in args:
        return build(next((a.split("=", 1)[1] for a in args if a.startswith("--features=")), None))
    if "--scipy-timing" in args:
        return scipy_timing(args[args.index("--scipy-timing") + 1])
    rows_mode = "--rows" in args
    if "--list" in args:
        print("checks:", " ".join(CHECKS))
        print("--rows checks:", " ".join(ROW_CHECKS))
        print("--rows controls:", " ".join(ROW_CONTROLS))
        return
    only, ls, print_rows = None, [2, 4, 5, 6], []
    for a in args:
        if a.startswith("--mutate="):
            MUTATE = a.split("=", 1)[1]
        if a.startswith("--scipy-file="):
            global SCIPY_FILE
            SCIPY_FILE = a.split("=", 1)[1]
        if a.startswith("--only="):
            only = [x for x in a.split("=", 1)[1].split(",") if x]   # --only= alone: no check
        if a.startswith("--l="):
            ls = [int(x) for x in a.split("=", 1)[1].split(",")]
        if a.startswith("--print-rows="):
            print_rows = [int(x) for x in a.split("=", 1)[1].split(",")]
        if a.startswith("--jsonl-dir="):
            JDIR = a.split("=", 1)[1]
        if a.startswith("--raw="):
            RAW = a.split("=", 1)[1]
    if not os.path.exists(BIN):
        build()
    checks, controls = (ROW_CHECKS, ROW_CONTROLS) if rows_mode else (CHECKS, CONTROLS)
    if MUTATE and MUTATE not in controls:
        sys.exit(f"unknown control {MUTATE}; known: {' '.join(controls)}")
    if rows_mode:
        default = [n for n in checks if re.sub(r"\D", "", n) == "" or int(re.sub(r"\D", "", n)) in ls or n == "RX" and 2 in ls]
    else:
        default = list(checks)
    names = only if only is not None else (controls[MUTATE] if MUTATE else default)
    # a control's targets are report names (CAP5b is a clause of the check CAP5): run each check once
    names = list(dict.fromkeys(n if n in checks else re.sub(r"[a-z]$", "", n) for n in names))
    if JDIR:
        os.makedirs(JDIR, exist_ok=True)
    t0 = time.time()
    if rows_mode:
        print(f"# taylor_polar.py --rows: start {time.strftime('%Y-%m-%d %H:%M:%S')}, load average {os.getloadavg()[0]:.2f}, one process, nice 10; mutate = {MUTATE}", flush=True)
    for n in names:
        try:
            checks[n]()
        except Exception as e:  # a crash under a mutation is a failure of the check, and is reported as one
            report(n, False, f"raised {type(e).__name__}: {str(e)[:300]}")
    for l in print_rows:
        gm = ROW_GAMMAS[l]
        setup, got = get_rows(l, gm, "dd", VKEYS)
        print(f"# l = {l}: rows from the crate (lane dd), no stored file, not compared")
        print("\n".join(format_rows(gm, got, setup)), flush=True)
    passed = [n for n, s in RESULTS if s is True]
    failed = [n for n, s in RESULTS if s is False]
    skipped = [n for n, s in RESULTS if s is None]
    print(f"SUMMARY: {len(passed)} passed, {len(failed)} failed, {len(skipped)} skipped" + (f" ({', '.join(skipped)}: not run)" if skipped else "")
          + f" of {len(RESULTS)} checks" + (f"; failed: {', '.join(failed)}" if failed else "") + f"; mutate = {MUTATE}; {time.time() - t0:.0f} s")
    if MUTATE:
        not_fired = [n for n in controls[MUTATE] if not any(r.startswith(n) and s is False for r, s in RESULTS)]
        print(f"CONTROL {MUTATE}: " + ("every targeted check failed, as it must" if not not_fired else f"NOT FIRED on {', '.join(not_fired)}"))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()

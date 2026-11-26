#!/usr/bin/env python
"""#186: the polar matching determinant D(gamma) by a Rust Taylor integrator (crate taylor_polar/).

    python taylor_polar.py                    build if needed, run every check, exit 1 if any fails
    python taylor_polar.py --only=T1,P1       run the named checks only
    python taylor_polar.py --mutate=NAME      run with one deliberate error (a control: the check it targets must fail)
    python taylor_polar.py --gen              write taylor_polar/src/jgen.rs from polar_matrix_generated.py
    python taylor_polar.py --build            cargo build --release (nice 10, -j 2)
    python taylor_polar.py --scipy-timing OUT run the scipy path (convergence.py's base variant) at the timing points,
                                              appending one JSON line per evaluation to OUT (long: about 30 min)
    python taylor_polar.py --list             the checks and the controls

jgen.rs holds J(r; s) of polar_matrix_generated.py (written by reduce_polar.py, the symbolic source) as generic Rust
over truncated power series: every statement of the Python function is translated, operation by operation, by walking
its syntax tree. The symbolic source is shared with polar_solver.py, not independent of it. The checks below name what
they share with the scipy path and what they do not.
"""
import ast
import json
import math
import os
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


def build():
    env = dict(os.environ, CARGO_NET_GIT_FETCH_WITH_CLI="true")
    subprocess.check_call(["nice", "-n", "10", "cargo", "build", "--release", "-j", "2"], cwd=CRATE, env=env)


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
    RESULTS.append((name, ok))
    print(f"CHECK {name}: {'PASS' if ok else 'FAIL'}  {text}", flush=True)


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
        print("CHECK TS: SKIPPED  no scipy timing file (--scipy-file=PATH)")
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
    "L1": check_L1, "V1": check_V1, "V2": check_V2, "V3": check_V3, "TS": check_TS,
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
    "moncrief": ["A2"],                   # a sign in the Moncrief map
    "swapcols": ["A2"],                   # the scalar plane's two start columns swapped (orientation)
    "mext0": ["A2"],                      # the exterior mass set to zero
    "qstart": ["A2"],                     # the Riccati started at r_m + 1 instead of r_m + 40/s + 20 l/s + 50
}


def main():
    global MUTATE
    args = sys.argv[1:]
    if "--gen" in args:
        return gen()
    if "--build" in args:
        return build()
    if "--scipy-timing" in args:
        return scipy_timing(args[args.index("--scipy-timing") + 1])
    if "--list" in args:
        print("checks:", " ".join(CHECKS))
        return
    only = None
    for a in args:
        if a.startswith("--mutate="):
            MUTATE = a.split("=", 1)[1]
        if a.startswith("--scipy-file="):
            global SCIPY_FILE
            SCIPY_FILE = a.split("=", 1)[1]
        if a.startswith("--only="):
            only = a.split("=", 1)[1].split(",")
    if not os.path.exists(BIN):
        build()
    if MUTATE and MUTATE not in CONTROLS:
        sys.exit(f"unknown control {MUTATE}; known: {' '.join(CONTROLS)}")
    names = only or (CONTROLS[MUTATE] if MUTATE else list(CHECKS))
    t0 = time.time()
    for n in names:
        try:
            CHECKS[n]()
        except Exception as e:  # a crash under a mutation is a failure of the check, and is reported as one
            report(n, False, f"raised {type(e).__name__}: {str(e)[:300]}")
    bad = [n for n, ok in RESULTS if not ok]
    print(f"SUMMARY: {len(RESULTS) - len(bad)} of {len(RESULTS)} checks passed" + (f"; failed: {', '.join(bad)}" if bad else "")
          + f"; mutate = {MUTATE}; {time.time() - t0:.0f} s")
    if MUTATE:
        target_ok = [n for n in CONTROLS[MUTATE] if not any(r.startswith(n) and not ok for r, ok in RESULTS if r.startswith(n))]
        print(f"CONTROL {MUTATE}: " + ("every targeted check failed, as it must" if not target_ok else f"NOT FIRED on {', '.join(target_ok)}"))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

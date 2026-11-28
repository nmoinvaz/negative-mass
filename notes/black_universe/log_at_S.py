#!/usr/bin/env python3
"""Independent rederivation of the logarithm of the c = -1 growing mode at the bifurcation sphere S (#321).

Claim of #302: the perturbation with gauge-invariant Y = c2 V R(rho(Pi)) is not smooth at S; d_V delta phi ~ -(pi/sqrt2) c2 ln V.
Nothing here imports or reuses horizon_basis.py, full_slice.py, the slice constraint M or the Kruskal solvers.
Shared with #302: the background closed forms (A, R, phi = sqrt2 arctan rho; checked against my own Ricci tensor in D1),
Kruskal V = e^{kappa v}, U = -e^{-kappa u}, Pi = -VU, and, in B only, #269's exact family (checked against my own Ricci in B1).

Conventions: units theta = phi/sqrt2, so R_ab = -2 theta_a theta_b + Vp g_ab (eps = -1) and -2 box theta = Vp_theta.
Double-null metric  ds^2 = -2 Om dV dU + r^2 dOmega^2 + hv dV^2 + hu dU^2  (hv = hu = 0 is double-null gauge); Om = e^{-F}.
Y = dr - g dth, g = dr0/dTh0 (gauge invariant). The static background has horizon at rho = 0, A'(0) = 2 kappa, A''(0) = -2/b^2.

    log_at_S.py derive   [--mutate=NAME]   D1-D4   sympy and exact series arithmetic, kappa symbolic
    log_at_S.py series   [--mutate=NAME]   S1-S6   Route A: formal series of the linearized equations, with ln V terms
    log_at_S.py check    [--only=B2,..] [--mutate=NAME]   B1-B5, X1   Route B: exact family, outgoing rays, 30 digits
    log_at_S.py mutations DIR PREFIX       fired sets of the stored logs DIR/PREFIX<name>.log against the declared ones

Status words as in AGENTS.md. Exit 1 when a check fails.
"""
import sys
import time

import mpmath as mp
import sympy as sp
from sympy import QQ
from sympy.polys.matrices import DomainMatrix

MUT = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mutate=")), None)
ONLY = next((a.split("=", 1)[1].split(",") for a in sys.argv[1:] if a.startswith("--only=")), None)
RESULTS = []
T0 = time.time()


def want(cid):
    return ONLY is None or cid in ONLY


def check(cid, text, ok, detail=""):
    RESULTS.append((cid, bool(ok)))
    print(f"[{'pass' if ok else 'FAIL'}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def info(text):
    print(f"       {text}", flush=True)


# ----------------------------------------------------------------------------------------------------------------------
# exact arithmetic: truncated bivariate series in (V, U) over Q(kappa)
# ----------------------------------------------------------------------------------------------------------------------
KAP = sp.Symbol("kappa", positive=True)
DOM = QQ.frac_field(KAP)
TRUNC = 11


def cv(x):
    return DOM.from_sympy(sp.sympify(x))


class Ser:
    """sum c_ij V^i U^j, total degree <= TRUNC."""
    __slots__ = ("d",)

    def __init__(self, d=None):
        self.d = {k: v for k, v in (d or {}).items() if v != DOM.zero and sum(k) <= TRUNC}

    def _o(self, o):
        if isinstance(o, Ser):
            return o
        return Ser({(0, 0): o if DOM.of_type(o) else cv(o)})

    def __add__(self, o):
        o = self._o(o)
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d.get(k, DOM.zero) + v
        return Ser(d)

    __radd__ = __add__

    def __neg__(self):
        return Ser({k: -v for k, v in self.d.items()})

    def __sub__(self, o):
        return self + (-self._o(o))

    def __rsub__(self, o):
        return self._o(o) - self

    def __mul__(self, o):
        o = self._o(o)
        d = {}
        for (i, j), a in self.d.items():
            for (p, q), b in o.d.items():
                if i + j + p + q <= TRUNC:
                    k = (i + p, j + q)
                    d[k] = d.get(k, DOM.zero) + a * b
        return Ser(d)

    __rmul__ = __mul__

    def inv(self):
        c0 = self.d.get((0, 0), DOM.zero)
        assert c0 != DOM.zero
        ic = DOM.one / c0
        rest = Ser({k: v for k, v in self.d.items() if k != (0, 0)})
        u = Ser({(0, 0): ic})
        for _ in range(TRUNC + 1):
            u = Ser({(0, 0): ic}) - rest * u * ic
        return u

    def __truediv__(self, o):
        return self * self._o(o).inv()

    def __rtruediv__(self, o):
        return self._o(o) * self.inv()

    def __pow__(self, n):
        if n >= 0:
            r = Ser({(0, 0): DOM.one})
            for _ in range(n):
                r = r * self
            return r
        return self.inv() ** (-n)

    def dV(self):
        return Ser({(i - 1, j): v * i for (i, j), v in self.d.items() if i > 0})

    def dU(self):
        return Ser({(i, j - 1): v * j for (i, j), v in self.d.items() if j > 0})

    def jet(self, a, b):
        s = self
        for _ in range(a):
            s = s.dV()
        for _ in range(b):
            s = s.dU()
        return s

    def mindeg(self):
        return min([sum(k) for k in self.d], default=None)


def pi_to_ser(coefs):
    """series in Pi (list over DOM) -> bivariate, Pi = -V U"""
    return Ser({(k, k): c * (DOM.one if k % 2 == 0 else -DOM.one) for k, c in enumerate(coefs)})


def u_mul(a, b, n):
    out = [DOM.zero] * (n + 1)
    for i, x in enumerate(a):
        if x == DOM.zero:
            continue
        for j, y in enumerate(b):
            if i + j > n:
                break
            out[i + j] += x * y
    return out


def u_compose(coefs_rho, rho_pi, n):
    out = [DOM.zero] * (n + 1)
    for c in reversed(coefs_rho):
        out = u_mul(out, rho_pi, n)
        out[0] += c
    return out


def series_in_rho(expr, rho, n):
    p = sp.Poly(sp.expand(sp.series(expr, rho, 0, n + 1).removeO()), rho)
    co = [DOM.zero] * (n + 1)
    for (k,), c in p.terms():
        co[k] = cv(c)
    return co


def build_background(b, G0=1):
    """series in Pi of the member (b, kappa) of the static family with the horizon at the throat rho = 0:
    R^2 = b^2 + rho^2, theta = arctan(rho/b), A = 1 + c R^2 + kappa Z, c = -1/b^2, Z = rho + R^2 arctan(rho/b)/b,
    Pi = e^{2 kappa x}, dx = drho/A, normalized by rho = G0 Pi + ...; Om0 = A/(2 kappa^2 Pi) from -A du dv."""
    rho = sp.Symbol("rho")
    b, G0 = sp.Rational(b), sp.Rational(G0)
    R2 = b ** 2 + rho ** 2
    Z = rho + R2 / b * sp.atan(rho / b)
    A = 1 - R2 / b ** 2 + KAP * Z
    nPi = TRUNC // 2 + 3
    nr = 2 * nPi + 4
    Aco = series_in_rho(A, rho, nr)
    if MUT == "a2":
        Aco[2] = Aco[2] * cv(sp.Rational(101, 100))
    kd = cv(KAP)
    g = [DOM.zero, cv(G0)] + [DOM.zero] * nPi
    for k in range(2, nPi + 2):
        rp = g[:k]
        pw = u_mul(rp, rp, k)
        tot = DOM.zero
        for n in range(2, k + 1):
            tot += Aco[n] * pw[k]
            pw = u_mul(pw, rp, k)
        g[k] = tot / (2 * kd) / (k - 1)
    rho_pi = g[:nPi + 2]
    N1 = nPi + 1
    out = {"rho": rho_pi, "Aco": Aco, "G0": G0}
    comp = lambda co: u_compose(co, rho_pi, N1)  # noqa: E731
    out["A"] = comp(Aco)
    assert out["A"][0] == DOM.zero
    out["Om"] = [x / ((1 if MUT == "om" else 2) * kd ** 2) for x in out["A"][1:]]
    out["r"] = comp(series_in_rho(sp.sqrt(R2), rho, nr))
    out["th"] = comp(series_in_rho(sp.atan(rho / b), rho, nr))
    Vexpr = (1 - sp.diff(rho * A, rho)) / R2                      # from E_thth of the static metric (D1)
    th_rho = b / R2
    V1 = sp.diff(Vexpr, rho) / th_rho
    V2 = sp.diff(V1, rho) / th_rho
    for nm, e in (("P0", Vexpr), ("P1", V1), ("P2", V2)):
        out[nm] = comp(series_in_rho(e, rho, nr))
    return out


# ----------------------------------------------------------------------------------------------------------------------
# the field equations by a brute-force Ricci tensor (sympy), and their linearization
# ----------------------------------------------------------------------------------------------------------------------
EPS_SIGN = -1 if MUT == "eps" else 1          # +1: phantom (eps = -1), the sign of 2 theta_a theta_b in E_ab


def ricci(g, gi, X):
    n = len(X)
    Gam = [[[sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(n)) / 2
             for c in range(n)] for b in range(n)] for a in range(n)]

    def ric(b, d):
        s = sum(sp.diff(Gam[a][b][d], X[a]) - sp.diff(Gam[a][b][a], X[d]) for a in range(n))
        s += sum(Gam[a][a][k] * Gam[k][b][d] - Gam[a][d][k] * Gam[k][b][a] for a in range(n) for k in range(n))
        return s

    return ric


def dn_equations():
    """E_ab = R_ab + 2 s theta_a theta_b - Vp g_ab (VV, VU, UU, thth) and -2 s box theta - Vp_theta, with the fields
    r0 + eta dr, Om0 + eta dOm, Th0 + eta dth, hv = eta hv, hu = eta hu; returns (E at eta = 0, d/d eta at 0)."""
    V, U, ang, ph, eta = sp.symbols("V U ang ph eta")
    X = [V, U, ang, ph]
    fn = lambda n: sp.Function(n)(V, U)  # noqa: E731
    r0, Om0, Th0, dr, dOm, dth, hv, hu = [fn(n) for n in ("r0", "Om0", "Th0", "dr", "dOm", "dth", "hv", "hu")]
    P0, P1, P2 = sp.symbols("P0 P1 P2")
    r, Om, Th = r0 + eta * dr, Om0 + eta * dOm, Th0 + eta * dth
    g = sp.zeros(4, 4)
    g[0, 0], g[1, 1], g[0, 1], g[1, 0] = eta * hv, eta * hu, -Om, -Om
    g[2, 2], g[3, 3] = r ** 2, r ** 2 * sp.sin(ang) ** 2
    gi = sp.zeros(4, 4)
    det2 = g[0, 0] * g[1, 1] - g[0, 1] ** 2
    gi[0, 0], gi[1, 1], gi[0, 1], gi[1, 0] = g[1, 1] / det2, g[0, 0] / det2, -g[0, 1] / det2, -g[0, 1] / det2
    gi[2, 2], gi[3, 3] = 1 / g[2, 2], 1 / g[3, 3]
    ric = ricci(g, gi, X)
    Vp, Vtp = P0 + eta * P1 * dth, P1 + eta * P2 * dth
    dT = [sp.diff(Th, x) for x in X]
    E = {}
    for nm, (i, j) in {"VV": (0, 0), "VU": (0, 1), "UU": (1, 1), "thth": (2, 2)}.items():
        E[nm] = ric(i, j) + 2 * EPS_SIGN * dT[i] * dT[j] - Vp * g[i, j]
    sq = sp.sqrt(-g.det() / sp.sin(ang) ** 2) * sp.sin(ang)
    box = sum(sp.diff(sq * gi[a, b] * dT[b], X[a]) for a in range(4) for b in range(4)) / sq
    E["sc"] = -2 * EPS_SIGN * box - Vtp
    E0 = {k: e.subs(eta, 0) for k, e in E.items()}
    E1 = {k: sp.diff(e, eta).subs(eta, 0) for k, e in E.items()}
    for k in ("thth", "sc"):
        E0[k], E1[k] = sp.simplify(E0[k]), sp.simplify(E1[k])
    return E0, E1, (V, U)


FIELDS = ["dr", "dOm", "dth", "hv", "hu"]
BGF = ["r0", "Om0", "Th0"]


def to_jets(expr, VU):
    V, U = VU
    reps = {}
    for d in expr.atoms(sp.Derivative):
        a = b = 0
        for var, n in d.variable_count:
            if var == V:
                a += n
            else:
                b += n
        reps[d] = sp.Symbol(f"{d.expr.func.__name__}_{a}_{b}")
    for f in FIELDS + BGF:
        reps[sp.Function(f)(V, U)] = sp.Symbol(f"{f}_0_0")
    return expr.xreplace(reps)


def operators():
    E0, E1, VU = dn_equations()
    E0j = {k: to_jets(e, VU) for k, e in E0.items()}
    ops = {}
    for k, e in E1.items():
        ej = to_jets(e, VU)
        ops[k] = []
        for f in FIELDS:
            for sym in sorted([s for s in ej.free_symbols if s.name.startswith(f + "_")], key=lambda s: s.name):
                a, b = map(int, sym.name.split("_")[1:])
                ops[k].append((f, a, b, sp.diff(ej, sym)))
    return E0j, ops


def evaluate(expr, env, memo):
    def ev(e):
        if e in memo:
            return memo[e]
        if e.is_Symbol:
            out = env[e.name]
        elif e.is_Rational:
            out = Ser({(0, 0): cv(e)})
        elif e.is_Add:
            out = ev(e.args[0])
            for a in e.args[1:]:
                out = out + ev(a)
        elif e.is_Mul:
            out = ev(e.args[0])
            for a in e.args[1:]:
                out = out * ev(a)
        elif e.is_Pow and e.args[1].is_Integer:
            out = ev(e.args[0]) ** int(e.args[1])
        else:
            raise ValueError(e)
        memo[e] = out
        return out

    return ev(expr)


def make_env(bg):
    env = {}
    for nm, key in (("r0", "r"), ("Om0", "Om"), ("Th0", "th")):
        s = pi_to_ser(bg[key])
        for a in range(4):
            for b in range(4):
                env[f"{nm}_{a}_{b}"] = s.jet(a, b)
    for nm in ("P0", "P1", "P2"):
        env[nm] = pi_to_ser(bg[nm])
    return env


def g_theta(bg):
    """dr0/dTh0 = r0_Pi / Th0_Pi as a series in Pi (list), lifted to (V, U)."""
    rp = [bg["r"][i] * i for i in range(1, len(bg["r"]))]
    tp = [bg["th"][i] * i for i in range(1, len(bg["th"]))]
    n = len(tp)
    q, rem = [], rp[:] + [DOM.zero] * 3
    for i in range(n - 1):
        qi = rem[i] / tp[0]
        q.append(qi)
        for j in range(n):
            if i + j < len(rem):
                rem[i + j] -= qi * tp[j]
    if MUT == "gth":
        q = [x * cv(sp.Rational(101, 100)) for x in q]
    return pi_to_ser(q)


class Problem:
    def __init__(self, b, N, G0=1, ops=None):
        global TRUNC
        TRUNC = N + 4
        self.b, self.N = sp.Rational(b), N
        self.bg = build_background(b, G0)
        self.env = make_env(self.bg)
        self.E0j, self.ops = ops if ops is not None else operators()
        memo = {}
        self.coef = {k: [(f, a, bb, evaluate(co, self.env, memo)) for f, a, bb, co in lst] for k, lst in self.ops.items()}
        self.r0, self.Om0, self.Th0 = [pi_to_ser(self.bg[x]) for x in ("r", "Om", "th")]
        self.g = g_theta(self.bg)
        self.memo = memo

    def FPi(self):
        """dF/dPi at S, F = -ln Om0, read from my Om0 series."""
        return -self.bg["Om"][1] / self.bg["Om"][0]

    def A2(self):
        return 2 * self.bg["Aco"][2]            # A''(0)


def mono_deriv(i, j, k, a, b):
    """d_V^a d_U^b of V^i U^j (ln V)^k as {(i, j, k): coefficient}"""
    terms = {(i, j, k): 1}
    for _ in range(a):
        new = {}
        for (ii, jj, kk), c in terms.items():
            if ii != 0:
                new[(ii - 1, jj, kk)] = new.get((ii - 1, jj, kk), 0) + c * ii
            if kk > 0:
                new[(ii - 1, jj, kk - 1)] = new.get((ii - 1, jj, kk - 1), 0) + c * kk
        terms = new
    for _ in range(b):
        new = {}
        for (ii, jj, kk), c in terms.items():
            if jj != 0:
                new[(ii, jj - 1, kk)] = new.get((ii, jj - 1, kk), 0) + c * jj
        terms = new
    return terms


def build_system(P, fields, log, N=None):
    """columns (field, i, j, k); rows (equation, m, n, k) for outputs of total degree m + n <= N - 2."""
    N = N or P.N
    cols = [(f, i, j, k) for f in fields for i in range(N + 1) for j in range(N + 1 - i) for k in ((0, 1) if log else (0,))]
    cidx = {c: n for n, c in enumerate(cols)}
    rows = {}
    for eq, lst in P.coef.items():
        for fn, a, b, C in lst:
            if fn not in fields:
                continue
            for (f, i, j, k), ci in cidx.items():
                if f != fn:
                    continue
                for (ii, jj, kk), cc in mono_deriv(i, j, k, a, b).items():
                    for (p, q), val in C.d.items():
                        m, n = ii + p, jj + q
                        if m + n > N - 2:
                            continue
                        r = rows.setdefault((eq, m, n, kk), {})
                        r[ci] = r.get(ci, DOM.zero) + val * DOM.convert(cc)
    return cols, cidx, rows


def y_rows(P, cidx, log, ymax):
    out = {}
    for i in range(ymax + 1):
        for j in range(ymax + 1 - i):
            for k in ((0, 1) if log else (0,)):
                r = {}
                c = cidx.get(("dr", i, j, k))
                if c is not None:
                    r[c] = DOM.one
                for (p, q), gv in P.g.d.items():
                    c2 = cidx.get(("dth", i - p, j - q, k))
                    if c2 is not None:
                        r[c2] = r.get(c2, DOM.zero) - gv
                out[(i, j, k)] = r
    return out


def to_dm(rows, ncols, keys=None):
    keys = keys or list(rows.keys())
    dok = {}
    for n, kk in enumerate(keys):
        for c, v in rows[kk].items():
            if v != DOM.zero:
                dok[(n, c)] = v
    return DomainMatrix.from_dok(dok, (len(keys), ncols), DOM), keys


def dod(M):
    return M.to_sdm().to_dod() if hasattr(M.to_sdm(), "to_dod") else dict(M.to_sdm())


def relations(P, fields, ymax):
    """linear relations among the Y coefficients (i + j <= ymax) of all smooth solutions of the truncated system."""
    cols, cidx, rows = build_system(P, fields, False)
    A, _ = to_dm(rows, len(cols))
    Nb = A.nullspace()
    yr = y_rows(P, cidx, False, ymax)
    Pm, ykeys = to_dm(yr, len(cols))
    rel = (Pm * Nb.transpose()).transpose().nullspace().rref()[0]
    return rel, [k[:2] for k in ykeys], Nb.shape[0]


def in_rowspace(rel, ykeys, form):
    vec = [form.get(k, DOM.zero) for k in ykeys]
    r1 = rel.rank()
    return DomainMatrix.vstack(rel, DomainMatrix([vec], (1, len(ykeys)), DOM)).rank() == r1


def target_series(P, kind, N):
    """Y coefficients of c V R / R(0), c U R / R(0) (horizon value c V resp. c U) or R."""
    R = P.bg["r"]
    out = {}
    for n, rn in enumerate(R):
        if 2 * n + 1 > N + 1:
            break
        sgn = DOM.one if n % 2 == 0 else -DOM.one
        if kind == "R":
            out[(n, n)] = rn * sgn
        else:
            val = rn * sgn / R[0]
            out[(n + 1, n) if kind == "VR" else (n, n + 1)] = val
    return {k: v for k, v in out.items() if v != DOM.zero}


def solve_affine(P, fields, log, target, ymax=None):
    N = P.N
    ymax = N if ymax is None else ymax
    cols, cidx, rows = build_system(P, fields, log)
    yr = y_rows(P, cidx, log, ymax)
    allrows, rhs = dict(rows), {}
    for key, r in yr.items():
        allrows[("Y",) + key] = r
        rhs[("Y",) + key] = target.get(key[:2], DOM.zero) if key[2] == 0 else DOM.zero
    keys = list(allrows.keys())
    A, _ = to_dm(allrows, len(cols), keys)
    bvec = DomainMatrix([[rhs.get(kk, DOM.zero)] for kk in keys], (len(keys), 1), DOM)
    rref, piv = A.hstack(bvec).rref()
    nc = len(cols)
    out = dict(cols=cols, cidx=cidx, consistent=nc not in piv, nullity=nc - len([p for p in piv if p < nc]))
    if out["consistent"]:
        sd = dod(rref)
        x = [DOM.zero] * nc
        for row, pc in enumerate(piv):
            x[pc] = sd.get(row, {}).get(nc, DOM.zero)
        out["xp"] = x
        Nb = A.nullspace()
        out["Nb"] = dod(Nb)
    return out


def log_coeff(res, field="dth", i=1, j=0):
    """(particular value, spread over the solution space) of the coefficient of V^i U^j ln V in a field"""
    c = res["cidx"].get((field, i, j, 1))
    if c is None:
        return None, None
    spread = any(row.get(c, DOM.zero) != DOM.zero for row in res["Nb"].values())
    return res["xp"][c], spread


# ----------------------------------------------------------------------------------------------------------------------
# derive: D1 - D4
# ----------------------------------------------------------------------------------------------------------------------
_OPS = []


def get_ops():
    if not _OPS:
        _OPS.append(operators())
    return _OPS[0]


def derive_D1():
    """static chart, my Ricci: the member (b, kappa) solves E_tt, E_rr, scalar with Vp from E_thth."""
    t, rho, ang, ph = sp.symbols("t rho ang ph", real=True)
    b = sp.Symbol("b", positive=True)
    X = [t, rho, ang, ph]
    R2 = b ** 2 + rho ** 2
    th0 = sp.atan(rho / b)
    A = 1 - R2 / b ** 2 + KAP * (rho + R2 / b * sp.atan(rho / b))
    Af = sp.Function("A")(rho)
    g = sp.diag(-Af, 1 / Af, R2, R2 * sp.sin(ang) ** 2)
    gi = g.inv()
    ric = ricci(g, gi, X)
    Rtt, Rrr, Rthth = [sp.simplify(ric(i, i)) for i in (0, 1, 2)]
    Vsol = sp.simplify(Rthth / R2)
    dth = sp.diff(th0, rho)
    Ett = sp.simplify((Rtt - Vsol * g[0, 0]).subs(Af, A).doit())
    Err = sp.simplify((Rrr + 2 * EPS_SIGN * dth ** 2 - Vsol * g[1, 1]).subs(Af, A).doit())
    Vc = sp.simplify(Vsol.subs(Af, A).doit())
    box = sp.diff(R2 * A * dth, rho) / R2
    Esc = sp.simplify(-2 * EPS_SIGN * box - sp.diff(Vc, rho) / dth)
    Vser = (1 - sp.diff(rho * A, rho)) / R2
    dV = sp.simplify(Vc - Vser)
    A0, A1, A2 = [sp.simplify(sp.diff(A, rho, n).subs(rho, 0)) for n in (0, 1, 2)]
    ok = Ett == 0 and Err == 0 and Esc == 0 and dV == 0 and A0 == 0 and sp.simplify(A1 - 2 * KAP) == 0 and sp.simplify(A2 + 2 / b ** 2) == 0
    check("D1", "static chart, own Ricci: E_tt, E_rr and the scalar equation vanish with Vp from E_thth; Vp = (1 - (rho A)')/R^2",
          ok, f"residuals {Ett}, {Err}, {Esc}, {dV}; A(0) = {A0}, A'(0) = {A1}, A''(0) = {A2}")


def lie_gauge(P, xiV, xiU):
    f = {"dr": xiV * P.r0.dV() + xiU * P.r0.dU(), "dth": xiV * P.Th0.dV() + xiU * P.Th0.dU()}
    sgn = -1 if MUT == "gauge" else 1
    f["dOm"] = xiV * P.Om0.dV() + xiU * P.Om0.dU() + sgn * P.Om0 * (xiV.dV() + xiU.dU())
    f["hv"] = P.Om0 * xiU.dV() * (-2)
    f["hu"] = P.Om0 * xiV.dU() * (-2)
    return f


def det_poly(deg, seed):
    """deterministic rational polynomial in (V, U)"""
    d = {}
    for i in range(deg + 1):
        for j in range(deg + 1 - i):
            n = (7 * i + 3 * j + 5 * seed) % 11 - 5
            if n:
                d[(i, j)] = cv(sp.Rational(n, 1 + (i * j + seed) % 3))
    return Ser(d)


def residual_degree(P, f):
    out = {}
    for k, lst in P.coef.items():
        tot = Ser()
        for fn, a, b, C in lst:
            tot = tot + C * f[fn].jet(a, b)
        out[k] = tot.mindeg()
    return out


def derive_D2_D3_D4():
    ops = get_ops()
    rows2, rows3, rows4 = [], [], []
    for b, G0 in ((1, 1), (2, 1), (sp.Rational(3, 2), 3)):
        P = Problem(b, 7, G0, ops)
        lim = TRUNC - 3
        r0 = {k: (M if (M := P.E0j[k]) is None else evaluate(M, P.env, P.memo).mindeg()) for k in P.E0j}
        rows2.append((b, G0, r0))
        xiV, xiU = det_poly(5, 1), det_poly(5, 2)
        res = residual_degree(P, lie_gauge(P, xiV, xiU))
        rows3.append((b, G0, res))
        xv1 = Ser({(i, 0): c for (i, j), c in xiV.d.items() if j == 0})
        xu1 = Ser({(0, j): c for (i, j), c in xiU.d.items() if i == 0})
        f = lie_gauge(P, xv1, xu1)
        res_dn = residual_degree(P, f)
        y = f["dr"] - P.g * f["dth"]
        rows4.append((b, G0, y.mindeg(), res_dn))
    lim = TRUNC - 3
    ok2 = all((m is None or m > lim) for _, _, r in rows2 for m in r.values())
    ok3 = all((m is None or m > lim) for _, _, r in rows3 for m in r.values()) and all((m is None or m > lim) for *_, r in rows4 for m in r.values())
    check("D2", "Kruskal chart: the background series (Om0 = A/(2 kappa^2 Pi), r0, Th0, Vp) solves VV, VU, UU, thth and the scalar equation",
          ok2, "lowest degree of a nonzero residual (b, G0, equation: degree; bar > %d): " % lim + "; ".join(f"b={b} G0={G0}: " + ",".join(f"{k}:{m}" for k, m in r.items()) for b, G0, r in rows2))
    check("D3", "linearized operator: pure gauge (xi^V, xi^U arbitrary in (V, U), and xi^V(V), xi^U(U)) solves all five linearized equations",
          ok3, "lowest degree of a nonzero residual (bar > %d): " % lim + "; ".join(f"b={b} G0={G0}: " + ",".join(f"{k}:{m}" for k, m in r.items()) for b, G0, r in rows3))
    okY = all(m is None or m > TRUNC - 3 for _, _, m, _ in rows4)
    check("D4", "Y = dr - g dth vanishes on pure gauge (gauge invariant)", okY,
          "lowest degree of Y(pure gauge): " + ", ".join(f"b={b}: {m}" for b, _, m, _ in rows4) + f"; bar > {TRUNC - 3}")


def run_derive():
    if MUT in ("eps", None, "a2", "om", "gauge", "gth"):
        pass
    derive_D1()
    derive_D2_D3_D4()


# ----------------------------------------------------------------------------------------------------------------------
# series: S1 - S6 (Route A)
# ----------------------------------------------------------------------------------------------------------------------
F3 = ["dr", "dOm", "dth"]
F5 = F3 + ["hv", "hu"]
KD = cv(KAP)


def fmt(x):
    return str(DOM.to_sympy(x))


def forms(FPi):
    """linear forms in the Y coefficients y_ij of the smooth relations; Y_N = (d_V + d_U) Y on V = s, U = -s"""
    m = 2 if MUT == "fpi" else 1
    h = FPi * m / 2
    return {
        "y21 = (F_Pi/2) y10": {(2, 1): DOM.one, (1, 0): -h},
        "y12 = (F_Pi/2) y01": {(1, 2): DOM.one, (0, 1): -h},
        "Y_N''(0) = -F_Pi Y_N(0)": {(3, 0): cv(6), (2, 1): cv(-2), (1, 2): cv(-2), (0, 3): cv(6), (1, 0): FPi * m, (0, 1): FPi * m},
    }


def series_S1():
    for tag, fields, N, ymax in (("double-null gauge (dr, dOm, dth)", F3, 7, 4), ("general smooth gauge (dr, dOm, dth, hv, hu)", F5, 5, 3)):
        P = Problem(1, N, 1, get_ops())
        rel, ykeys, nul = relations(P, fields, ymax)
        FPi = P.FPi()
        ok, det = True, []
        for nm, form in forms(FPi).items():
            r = in_rowspace(rel, ykeys, form)
            ok &= r
            det.append(f"{nm}: {r}")
        hor = [in_rowspace(rel, ykeys, {(i, 0): DOM.one}) and in_rowspace(rel, ykeys, {(0, i): DOM.one}) for i in range(2, ymax + 1)]
        ok &= all(hor)
        neg = not in_rowspace(rel, ykeys, {(2, 1): DOM.one})            # y21 = 0 (what V R has) is excluded
        ok &= neg
        check("S1" if tag.startswith("double") else "S1g", f"{tag}: smooth solutions (nullity {nul} of the truncated system, N = {N}) obey y21 = (F_Pi/2) y10, y12 = (F_Pi/2) y01, "
                    f"y_i0 = y_0i = 0 (i >= 2); y21 = 0 is excluded", ok,
              f"F_Pi(0) = {fmt(FPi)} from my Om0; " + "; ".join(det) + f"; horizon laws {all(hor)}; y21 = 0 in the span: {not neg}")
        if tag.startswith("double"):
            rl = dod(rel)
            for row in range(rel.shape[0]):
                terms = [f"({fmt(v)}) y{ykeys[c][0]}{ykeys[c][1]}" for c, v in sorted(rl.get(row, {}).items())]
                info("relation: " + " + ".join(terms) + " = 0")
            P6 = P
            tg = {(1, 0): DOM.one}
            if MUT == "hor":
                tg[(2, 0)] = DOM.one                    # contradicts the horizon law y20 = 0
            res = solve_affine(P6, F3, False, tg, ymax=2 if MUT == "hor" else 1)
            check("S1b", "a smooth solution with delta r = c2 V on H+ and 0 on H- (y00 = y01 = 0, y10 = 1, nothing else imposed) exists", res["consistent"],
                  f"nullity {res['nullity']}; its Y = V + (F_Pi/2) V^2 U + ... = V + {fmt(FPi / 2)} V^2 U + ...; V R = V + O(V^3 U^2)")


def series_S2():
    P = Problem(1, 7, 1, get_ops())
    out = []
    for kind, ymax in (("VR", 7), ("UR", 7), ("VR", 3), ("R", 7)):
        tgt = target_series(P, kind, 7)
        if MUT == "s2t" and kind == "VR" and ymax == 3:
            tgt[(2, 1)] = tgt.get((2, 1), DOM.zero) + P.FPi() / 2        # the V^2 U term smooth solutions need
        res = solve_affine(P, F3, False, tgt, ymax=ymax)
        out.append((kind, ymax, res["consistent"]))
    ok = (not out[0][2]) and (not out[1][2]) and (not out[2][2]) and out[3][2]
    check("S2", "smooth ansatz: no solution with Y = c2 V R, none with Y = c2' U R (also when only the degree <= 3 coefficients of Y are imposed), a solution with Y = c1 R",
          ok, "consistent: " + ", ".join(f"Y={k} (degree <= {m}): {c}" for k, m, c in out))


def log_solve(P, fields, kind, ymax=None, log=True, tgt=None):
    if MUT == "nolog":
        log = False
    if MUT == "rhs" and kind == "VR":
        tgt = {(1, 0): DOM.one}
    res = solve_affine(P, fields, log, tgt if tgt is not None else target_series(P, kind, P.N), ymax=ymax)
    return res


def series_S3():
    rows = []
    ok = True
    for tag, fields, N in (("double-null gauge", F3, 7), ("general smooth gauge, logs also in hv, hu", F5, 5)):
        P = Problem(1, N, 1, get_ops())
        res = log_solve(P, fields, "VR")
        if not res["consistent"]:
            check("S3", f"{tag}: log ansatz has a solution with Y = c2 V R", False, "inconsistent")
            return
        q, spread = log_coeff(res)
        want_q = -DOM.one / KD
        good = q is not None and q == want_q and not spread
        low = [(c[0], c[1], c[2], fmt(res["xp"][n])) for n, c in enumerate(res["cols"]) if c[3] == 1 and c[0] == "dth" and c[1] + c[2] <= 1 and res["xp"][n] != DOM.zero]
        ok &= good
        rows.append(f"{tag}: coefficient of V ln V in delta theta = {fmt(q)}, spread over the solution space {spread}, nullity {res['nullity']}, log terms of delta theta of degree <= 1: {low}")
    check("S3", "log ansatz: Y = c2 V R has a solution; the coefficient of V ln V in delta theta is -c2/kappa for every solution, i.e. d_N delta phi = -sqrt2/kappa c2 ln|s| = -(pi/sqrt2) c2 ln|s| at kappa = 2/pi",
          ok, " | ".join(rows))


def series_S4():
    P = Problem(1, 7, 1, get_ops())
    res = log_solve(P, F3, "R")
    q, spread = (None, None) if not res["consistent"] else log_coeff(res)
    ok = res["consistent"] and q is not None and q == DOM.zero and not spread
    lows = [] if not res["consistent"] else [(c, fmt(res["xp"][n])) for n, c in enumerate(res["cols"]) if c[3] == 1 and c[0] == "dth" and res["xp"][n] != DOM.zero]
    check("S4", "planted: Y = c1 R (the static mass mode) through the same log ansatz: the coefficient of V ln V in delta theta is 0 for every solution",
          ok, f"q = {None if q is None else fmt(q)}, spread {spread}, log terms of delta theta in the particular solution: {lows}")


def series_S5():
    rows, ok = [], True
    for b, G0 in ((1, 1), (2, 1), (sp.Rational(3, 2), 1), (sp.Rational(1, 2), 1), (1, 3), (2, sp.Rational(5, 2))):
        P = Problem(b, 6, G0, get_ops())
        res = log_solve(P, F3, "VR")
        if not res["consistent"]:
            ok = False
            rows.append(f"b={b}, G0={G0}: inconsistent")
            continue
        q, spread = log_coeff(res)
        pred = P.A2() * (cv(sp.Rational(101, 100)) if MUT == "pred" else DOM.one) / (2 * KD)
        good = q is not None and q == pred and not spread
        ok &= good
        rows.append(f"b={b}, G0={G0}: A''(0) = {fmt(P.A2())}, q = {None if q is None else fmt(q)}, A''(0)/(2 kappa) = {fmt(pred)}: {'equal' if good else 'DIFFER'}")
    check("S5", "control: members (b, kappa) with A''(0) = -2/b^2 and kappa independent, and G0 changed: q = A''(0) c_h / (2 kappa), c_h = d_V delta r on H+", ok, " | ".join(rows))


def series_S6():
    P = Problem(1, 6, 1, get_ops())
    FPi = P.FPi()
    h = FPi * (2 if MUT == "fpi" else 1) / 2
    rows, ok = [], True
    for num, den in ((0, 1), (1, 2), (1, 1), (2, 1)):
        t = h * cv(sp.Rational(num, den))
        tg = {(1, 0): DOM.one, (2, 1): t}
        sm = solve_affine(P, F3, False, tg, ymax=3)["consistent"]
        res = log_solve(P, F3, "VR", ymax=3, tgt=tg)
        q, spread = log_coeff(res) if res["consistent"] else (None, None)
        if q is None:
            ok = False
            rows.append(f"t = {num}/{den} (F_Pi/2): no log solution")
            continue
        wantq = 2 * t - FPi
        good = (q == wantq) and not spread and (sm == (num == den))
        ok &= good
        rows.append(f"y21 = {num}/{den} (F_Pi/2) y10: smooth solution {sm}, q = {fmt(q)} (2 y21 - F_Pi y10 = {fmt(wantq)})")
    check("S6", "the log is the violation of the smooth relation: Y = V + t V^2 U + ... has q = 2 t - F_Pi, and a smooth solution iff t = F_Pi/2", ok, " | ".join(rows))


def run_series():
    for cid, fn in (("S1", series_S1), ("S2", series_S2), ("S3", series_S3), ("S4", series_S4), ("S5", series_S5), ("S6", series_S6)):
        if want(cid):
            fn()


# ----------------------------------------------------------------------------------------------------------------------
# check: B1 - B5, X1 (Route B): the exact family of #269 in a double-null chart built from the outgoing rays
# ----------------------------------------------------------------------------------------------------------------------
class Fam:
    """b = 1, c = -1, K = kappa = 2/pi; L = 1. Coordinates: w = ln V = kappa v. F_in(w, y, beta) is g_vv = -F in the ingoing chart."""

    def __init__(self):
        mp.mp.dps = 30
        self.kap = 2 / mp.pi
        self.Kr = mp.mpf("1.01") if MUT == "rate" else mp.mpf(1)        # K / kappa
        self.phis = mp.mpf("1.01") if MUT == "phi" else mp.mpf(1)       # scale of dtheta/drho

    def A(self, y):
        return -(y * y) + self.kap * (y + (y * y + 1) * mp.atan(y))

    def Fin(self, w, y, beta):
        if MUT == "beta":
            beta = -beta
        s = beta * mp.e ** (self.Kr * w)
        Kc = self.Kr * self.kap
        if MUT == "h":
            return self.A(y) - Kc * y * s / (1 + s)
        return self.A(y) - (1 + Kc * y) * s / (1 + s)

    def lam(self, w, beta):
        return mp.sqrt(1 + beta * mp.e ** (self.Kr * w))

    def dphi_dy(self, y):              # d theta / d rho times sqrt2 (phi = sqrt2 theta)
        return self.phis * mp.sqrt(2) / (1 + y * y)

    def ray(self, wb, rb, beta, kc=None):
        kc = kc or self.kap
        return mp.odefun(lambda tau, y: -self.Fin(wb - tau, y, beta) / (2 * kc), 0, rb)


def fam_ricci():
    v, rho, ang, ph = sp.symbols("v rho ang ph", real=True)
    X = [v, rho, ang, ph]
    Ff, rf, pf = sp.Function("F")(v, rho), sp.Function("r")(v, rho), sp.Function("p")(rho)
    g = sp.zeros(4, 4)
    g[0, 0], g[0, 1], g[1, 0], g[2, 2], g[3, 3] = -Ff, 1, 1, rf ** 2, rf ** 2 * sp.sin(ang) ** 2
    gi = g.inv()
    ric = ricci(g, gi, X)
    Vs, V1s = sp.symbols("Vs V1s")
    dT = [sp.diff(pf, x) for x in X]
    E = {}
    for nm, (i, j) in {"vv": (0, 0), "vr": (0, 1), "rr": (1, 1), "thth": (2, 2)}.items():
        E[nm] = ric(i, j) + 2 * EPS_SIGN * dT[i] * dT[j] - Vs * g[i, j]
    sq = rf ** 2 * sp.sin(ang)
    E["sc"] = -2 * EPS_SIGN * sum(sp.diff(sq * gi[a, b] * dT[b], X[a]) for a in range(4) for b in range(4)) / sq - V1s
    E = {k: sp.simplify(e.subs(ang, sp.pi / 3)) for k, e in E.items()}
    jets = {}
    for k, e in E.items():
        reps, keys = {}, []
        for d in sorted(e.atoms(sp.Derivative), key=lambda z: -z.derivative_count):
            cnt = {v: 0, rho: 0}
            for var, n in d.variable_count:
                cnt[var] += n
            sym = sp.Symbol(f"{d.expr.func.__name__}_{cnt[v]}_{cnt[rho]}")
            reps[d] = sym
            keys.append((d.expr.func.__name__, cnt[v], cnt[rho], sym))
        for f in (Ff, rf, pf):
            if f in e.atoms(sp.Function):
                sym = sp.Symbol(f"{f.func.__name__}_0_0")
                reps[f] = sym
                keys.append((f.func.__name__, 0, 0, sym))
        ee = e.xreplace(reps)
        jets[k] = (sp.lambdify([s for *_, s in keys] + [Vs, V1s], ee, "mpmath"), keys)
    return jets


def check_B1(fam):
    jets = fam_ricci()
    kap = fam.kap

    def member(beta):
        lam = lambda s: fam.lam(kap * s, beta)  # noqa: E731
        dlam = lambda s: mp.diff(lam, s)  # noqa: E731
        hterm = (lambda s: 0) if MUT == "h" else (lambda s: 1 / lam(s) ** 2 - 1)
        bsg = -1 if MUT == "beta" else 1
        lamb = lambda s: fam.lam(kap * s, bsg * beta)  # noqa: E731
        dlamb = lambda s: mp.diff(lamb, s)  # noqa: E731
        F = lambda s, y: fam.A(y) - 2 * y * dlamb(s) / lamb(s) + (0 if MUT == "h" else 1 / lamb(s) ** 2 - 1)  # noqa: E731
        r = lambda s, y: lam(s) * mp.sqrt(1 + y * y)  # noqa: E731
        p = lambda y: fam.phis * mp.atan(y)  # noqa: E731
        return F, r, p

    def Vp(y):
        Vf = lambda z: (1 - mp.diff(lambda u: u * fam.A(u), z)) / (1 + z * z)  # noqa: E731
        return Vf(y), (1 + y * y) * mp.diff(Vf, y)

    worst, rows = mp.mpf(0), []
    for beta in (mp.mpf("0.3"), mp.mpf("-0.4")):
        F, r, p = member(beta)
        for s, y in ((mp.mpf("-1"), mp.mpf("0.4")), (mp.mpf("0.2"), mp.mpf("-0.7")), (mp.mpf("-0.5"), mp.mpf("2.5"))):
            V0, V1 = Vp(y)
            for k, (f, keys) in jets.items():
                vals = []
                for nm, i, j, _ in keys:
                    fn = {"F": F, "r": r}.get(nm)
                    if nm == "p":
                        vals.append(mp.diff(p, y, j) if i == 0 else mp.mpf(0))
                    else:
                        vals.append(mp.diff(fn, (s, y), (i, j)))
                val = f(*vals, V0, V1)
                scale = max([abs(x) for x in vals] + [abs(V0), abs(V1), mp.mpf(1)])
                worst = max(worst, abs(val) / scale)
    check("B1", "the family r = lam R, F = A - 2 rho lam'/lam + 1/lam^2 - 1 solves vv, v-rho, rho-rho, thth and the scalar equation of my Ricci in the ingoing chart",
          worst < mp.mpf("1e-20"), f"worst relative residual {mp.nstr(worst, 3)} over 2 values of beta and 3 points; bar 1e-20")


def amplitudes(fam, eps):
    """c1, c2 of the family on the ray of label 0 (rho = 0 at V_b = 1): delta r = d_beta (lam R) at V = 1e-3 ... 1e-7."""
    sp_, sm_ = fam.ray(0, 0, eps), fam.ray(0, 0, -eps)
    out = []
    for k in (3, 5, 7):
        V = mp.mpf(10) ** (-k)
        w = mp.log(V)
        rp = fam.lam(w, eps) * mp.sqrt(1 + sp_(-w) ** 2)
        rm = fam.lam(w, -eps) * mp.sqrt(1 + sm_(-w) ** 2)
        out.append((V, (rp - rm) / (2 * eps)))
    return out


def growing_values(fam, Vb, U, eps, Vs):
    rb = -U * Vb
    wb = mp.log(Vb)
    sol = {s: fam.ray(wb, rb, s * eps) for s in (1, -1)}
    res = []
    for V in Vs:
        w = mp.log(V)
        d = []
        for s in (1, -1):
            y = sol[s](wb - w)
            d.append(fam.dphi_dy(y) * fam.Fin(w, y, s * eps) / (2 * fam.kap * V))
        res.append((d[0] - d[1]) / (2 * eps))
    return res


def slopes(vals, Vs):
    return [(vals[i + 1] - vals[i]) / (mp.log(Vs[i + 1]) - mp.log(Vs[i])) for i in range(len(vals) - 1)]


def Pi_of(fam, y):
    """Kruskal Pi(rho) = rho exp(2 kappa int_0^rho (1/A - 1/(2 kappa t)) dt), normalized by G0 = 1 (the same as the series)"""
    kap = fam.kap
    f = lambda t: (1 / fam.A(t) - 1 / (2 * kap * t)) if t != 0 else mp.mpf(0)  # noqa: E731
    return y * mp.e ** (2 * kap * mp.quad(f, [0, y]))


def static_values(fam, Lv, U, Vb, Vs, own_chart):
    """static mass mode r = sqrt(L) R, F = A + 1/L - 1: ray label U measured from the horizon rho_h(L); V = e^{kappa_L v} (own chart)."""
    dl = 1 / Lv - 1
    rh = mp.findroot(lambda y: fam.A(y) + dl, -dl / (2 * fam.kap))
    kL = mp.diff(fam.A, rh) / 2
    kc = kL if own_chart else fam.kap
    rb = rh - U * Vb
    wb = mp.log(Vb)
    sol = mp.odefun(lambda tau, y: -(fam.A(y) + dl) / (2 * kc), 0, rb)
    out = []
    for V in Vs:
        y = sol(wb - mp.log(V))
        out.append(fam.dphi_dy(y) * (fam.A(y) + dl) / (2 * kc * V))
    return out


def run_check():
    fam = Fam()
    t0 = time.time()
    if want("B1"):
        check_B1(fam)
        info(f"{time.time() - t0:.0f} s")
    eps = mp.mpf("0.3") if MUT == "step" else mp.mpf("1e-9")
    Vs = [mp.mpf(10) ** (-k) for k in (3, 4, 5, 6, 7, 8)]
    c2 = None
    if want("B5") or want("B2") or want("X1"):
        am = amplitudes(fam, eps)
        c2 = am[-1][1] / am[-1][0]
        c1 = am[-1][1] - c2 * am[-1][0]
        spread = max(abs(a[1] / a[0] / c2 - 1) for a in am)
        check("B5", "amplitude of the family: delta r = c1 + c2 V on the ray U = 0 with c1 = 0, c2 = 1/2 per unit beta", abs(c2 - mp.mpf(1) / 2) < mp.mpf("1e-8") and abs(c1) < mp.mpf("1e-12") and spread < mp.mpf("1e-6"),
              f"c2 = {mp.nstr(c2, 12)} (V = 1e-7), delta r / V over V = 1e-3, 1e-5, 1e-7 differs by {mp.nstr(spread, 3)}; c1 = {mp.nstr(c1, 3)}")
    if want("B2") or want("X1"):
        wantc = -mp.pi / mp.sqrt(2) * c2
        ok, rows, last = True, [], {}
        for Vb, U in ((mp.mpf(1), mp.mpf("-0.5")), (mp.mpf(1), mp.mpf(-1)), (mp.mpf("0.05"), mp.mpf("-0.5")), (mp.mpf("0.05"), mp.mpf(-1))):
            vals = growing_values(fam, Vb, U, eps, Vs)
            r = [s / wantc for s in slopes(vals, Vs)]
            last[(Vb, U)] = r[-1]
            good = abs(r[-1] - 1) < mp.mpf("1e-4") and abs(r[-1] - 1) < abs(r[-3] - 1)
            ok &= good
            rows.append(f"V_b = {mp.nstr(Vb, 3)}, U = {mp.nstr(U, 3)}: " + ", ".join(mp.nstr(x, 9) for x in r))
        check("B2", "d_V delta phi of the exact family at fixed U: slope against ln V over -(pi/sqrt2) c2, at V = 1e-4 ... 1e-8 (bar |ratio - 1| < 1e-4 at 1e-8 and falling)", ok, "; ".join(rows))
        info(f"{time.time() - t0:.0f} s")
        if not (ONLY and "B2" not in ONLY and "X1" not in ONLY):
            v2 = growing_values(fam, mp.mpf(1), mp.mpf("-0.5"), eps / 10, Vs)
            r2 = slopes(v2, Vs)[-1] / wantc
            check("X1", "cross-check of the two routes: the numeric slope of B2 equals Route A's exact -(sqrt2/kappa) c2 (S3) to 1e-5, and does not depend on the step in beta (eps/10)",
                  abs(last[(mp.mpf(1), mp.mpf("-0.5"))] - 1) < mp.mpf("1e-5") and abs(r2 - last[(mp.mpf(1), mp.mpf("-0.5"))]) < mp.mpf("1e-8"),
                  f"ratio at 1e-8: {mp.nstr(last[(mp.mpf(1), mp.mpf('-0.5'))], 9)} (eps) and {mp.nstr(r2, 9)} (eps/10); -(sqrt2/kappa) c2 = {mp.nstr(wantc, 12)}")
            info(f"{time.time() - t0:.0f} s")
    if want("B3"):
        own = MUT != "chart"
        rows, ok = [], True
        for U in (mp.mpf("-0.5"), mp.mpf(-1)):
            d = [(a - b) / (2 * mp.mpf("1e-9")) for a, b in zip(static_values(fam, 1 + mp.mpf("1e-9"), U, mp.mpf(1), Vs, own), static_values(fam, 1 - mp.mpf("1e-9"), U, mp.mpf(1), Vs, own))]
            sl = slopes(d, Vs)
            scale = mp.pi / mp.sqrt(2) / 2
            good = abs(sl[-1]) / scale < mp.mpf("1e-6")
            ok &= good
            rows.append(f"U = {mp.nstr(U, 3)}: slope / ((pi/sqrt2)/2) = " + ", ".join(mp.nstr(x / scale, 3) for x in sl))
        check("B3", "planted: the static mass mode (r = sqrt(L) R, F = A + 1/L - 1), in its own Kruskal chart V = e^{kappa_L v}, has no log: slope / ((pi/sqrt2) c) < 1e-6 at V = 1e-8", ok, "; ".join(rows))
        info(f"{time.time() - t0:.0f} s")
    if want("B4"):
        rows, ok = [], True
        for U in (mp.mpf("-0.5"), mp.mpf(-1)):
            d = [(a - b) / (2 * mp.mpf("1e-9")) for a, b in zip(static_values(fam, 1 + mp.mpf("1e-9"), U, mp.mpf(1), Vs, False), static_values(fam, 1 - mp.mpf("1e-9"), U, mp.mpf(1), Vs, False))]
            sl = slopes(d, Vs)[-1]
            sgn = -1 if MUT == "pgauge" else 1
            pred = -sgn * mp.sqrt(2) / (2 * fam.kap ** 2) * Pi_of(fam, -U)
            good = abs(sl / pred - 1) < mp.mpf("1e-3")
            ok &= good
            rows.append(f"U = {mp.nstr(U, 3)}: slope {mp.nstr(sl, 8)}, predicted sqrt2 U_K/(2 kappa^2) = {mp.nstr(pred, 8)}")
        check("B4", "the same mode in the chart with V = e^{kappa v} of the background (kappa not kappa_L) has the log of the gauge xi^V = (d kappa/kappa) V ln V, proportional to U_K = -Pi(rho_b)/V_b: sqrt2 U_K/(2 kappa^2) per unit dL, to 1e-3", ok, "; ".join(rows))
        info(f"{time.time() - t0:.0f} s")


# ----------------------------------------------------------------------------------------------------------------------
# mutations: name -> (command, the checks it is declared to fire); the effects were computed in development runs (raw README)
# ----------------------------------------------------------------------------------------------------------------------
MUTATIONS = {
    "eps": ("derive", ["D1", "D2", "D3"]), "a2": ("derive", ["D2", "D3"]), "om": ("derive", ["D2", "D3"]),
    "gauge": ("derive", ["D3"]), "gth": ("derive", ["D4"]),
    "fpi": ("series", ["S1", "S1g", "S6"]), "nolog": ("series", ["S3", "S4", "S5", "S6"]), "pred": ("series", ["S5"]),
    "rhs": ("series", ["S3", "S5", "S6"]), "hor": ("series", ["S1b"]), "s2t": ("series", ["S2"]),
    "rate": ("check", ["B1", "B2", "B5", "X1"]), "h": ("check", ["B1", "B2", "X1"]), "beta": ("check", ["B1", "B2", "X1"]),
    "phi": ("check", ["B1", "B2", "B4", "X1"]), "chart": ("check", ["B3"]), "pgauge": ("check", ["B4"]), "step": ("check", ["B5", "X1"]),
}


def report_mutations(raw, prefix):
    import os
    import re
    good = 0
    for name, (cmd, decl) in MUTATIONS.items():
        path = os.path.join(raw, f"{prefix}{name}.log")
        if not os.path.exists(path):
            print(f"{name}: no log {path}")
            continue
        fired = sorted(set(re.findall(r"^\[FAIL\] (\w+):", open(path).read(), re.M)))
        ok = decl is not None and fired == sorted(decl)
        good += ok
        print(f"{name:7s} ({cmd}): fired {','.join(fired) or '-'}; declared {','.join(sorted(decl)) if decl is not None else '?'}: {'equal' if ok else 'DIFFERENT'}")
    print(f"{good} of {len(MUTATIONS)} fired sets equal the declared ones")


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "mutations":
        report_mutations(sys.argv[2], sys.argv[3])
        return
    if MUT is not None and (MUT not in MUTATIONS or MUTATIONS[MUT][0] != cmd):
        sys.exit(f"mutation {MUT} is not declared for {cmd}")
    if cmd == "derive":
        run_derive()
    elif cmd == "series":
        run_series()
    elif cmd == "check":
        run_check()
    else:
        sys.exit(__doc__)
    bad = [c for c, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS)} checks, {len(RESULTS) - len(bad)} pass; failed: {bad}; mutation {MUT}; {time.time() - T0:.0f} s", flush=True)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""The boundary term of item 6's identity against E - P.xi_0 (#100), checked symbolically.

The identity of #87 (weitzenbock_conf.py), on g = psi^4 delta in the frame e_i = psi^-2 d_i, is
    |nabla~ eps|^2 - |D~ eps|^2 + 4 pi (mu |eps|^2 - J.xi) = div V,   div V = psi^-6 d_a(psi^4 V^a),
    V^i = Re<eps, nabla^i eps - g^i D eps> - (1/2) P_ij xi^j,   P_ij = K_ij - delta_ij tr K,   xi^j = <eps, g^0 g^j eps>,
with the gammas and <a, b> = a^dagger b of mass_bounds.md, nabla~_i = nabla_i + (1/2) K_ij g^j g^0 and D = g^i nabla_i. With
D~ eps = 0, integrated over the ball |x| < r with r -> infinity, it is item 6: (1/4 pi) int |nabla~ eps|^2 dV + int (mu
|eps|^2 - J.xi) dV = B, with the boundary term B = (1/4 pi) lim of the flux of V through the sphere, which item 6 and
phase4/RESULTS.md section 1 state to be E - P.xi_0 (|eps_0| = 1, xi_0^j = <eps_0, g^0 g^j eps_0>). This script computes B.

Data (asymptotics, with s = 1/r and n = x/r, every coefficient an independent symbol): psi = 1 + (m/2) s + psi_2(n) s^2 +
psi_3(n) s^3, psi_2 a polynomial of degree 2 in n (a dipole d.n, a quadrupole and a constant), psi_3 of degree 1; K_ij (frame
components) = psi^-6 K-hat_ij with K-hat_ij = (3/2) s^2 [P_i n_j + P_j n_i - (delta_ij - n_i n_j) P.n] (Bowen-York, tr = 0)
+ s^3 K3_ij(n); eps = eps_0 + eps_1(n) s + eps_2(n) s^2 with eps_0 a constant spinor of eight real symbols, eps_1 four
components of polynomials of degree <= 2 in n with independent complex coefficients, eps_2 of degree <= 1. Two extensions beyond
the issue: a K-hat term s^2 T_ij(n) with a trace (`bdy_vs_adm_T`, so that -delta_ij tr K is exercised), and an O(1) term
eta(n) in eps (the control `eps-O1`).

The flux. The covariant divergence of a vector V = V^a e_a is psi^-6 d_a(psi^4 V^a) (`div_formula` checks that against the
Christoffel symbols of psi^4 delta for arbitrary psi and V^a). Its volume element is psi^6 d^3x, so the integral of div V over
a shell is the integral of d_a(psi^4 V^a) d^3x, which is the flux of psi^4 V^a through the two coordinate spheres by Gauss in
flat R^3 (by hand): int over |x| = r of psi^4 V^a n_a r^2 dOmega, n the Euclidean unit normal. It is the covariant flux
too: the unit normal has frame components n_a and the area element is psi^4 r^2 dOmega. `planted_B_schwarzschild` checks the
weight psi^4 on an exact flux through a finite sphere. B is the coefficient of s^0 of (1/4 pi) int psi^4 V^a n_a r^2 dOmega.
E = (1/16 pi) int (d_j g_ij - d_i g_jj) n_i r^2 dOmega and P_j = (1/8 pi) int (K_ij - delta_ij tr K) n_i r^2 dOmega are the
ADM integrals of the same data, with the coordinate components g_ij = psi^4 delta_ij, K_ij = psi^4 (frame components), tr K =
delta^ij K_ij, written separately from V's code (`adm_E_series`, `adm_P_series`).

How. Every quantity is a polynomial in n_i, s and the symbols, truncated at s^3 (V to s^3, so the flux to order 1/r). The
Cartesian derivative at x = n/s is D_i f = s (d_ni f - n_i (n.d_n f + s d_s f)), which is consistent on the sphere. The sphere
integral is by exact moments, (1/4 pi) int n_x^a n_y^b n_z^c dOmega = (a-1)!!(b-1)!!(c-1)!!/(a+b+c+1)!! for even a, b, c.
Data of exactly these polynomial forms have a flux that is rational in s with a denominator a power of psi, analytic at s = 0,
so the limit r -> infinity is the s^0 coefficient. The V code is one function (`Vvec`) of a field representation, run in
three: the s-series, exact Taylor jets at a point (for the identity at exact points) and exact generic functions (for the
checks on arbitrary psi). By hand, for eps = eps_0 constant, D eps_0 = psi^-2 phislash eps_0 with phi = 2 d ln psi, Re<eps_0,
Om_i eps_0> = 0 and V^i = -2 psi^-3 d_i psi |eps_0|^2 - (1/2) P_ij xi_0^j, whose flux is m (1 + m/2r)|eps_0|^2 - (1 +
m/2r)^-2 P.xi_0 for psi = 1 + m/2r and K = psi^-6 Bowen-York; this is the planted case.

The eps_1 terms. With eta = eps - eps_0, the terms of V linear in eta at order r^-2 are Re<eps_0, (delta_ik - g^i g^k) d_k eta>.
The flux of the first term, n_i d_i eta, is not zero, but together they are -Re<eps_0, nslash dslash_tan eta>, because n_i
(delta_ik - g^i g^k) n_k = 0, and nslash dslash_tan = sum_{i<k} g^i g^k (n_i d_k - n_k d_i) (`eps1_hand` checks this for eps_1). It
is a sum of rotation generators with the antisymmetric g^i g^k. The generators are divergence-free on the sphere, so the
integral over the sphere is 0 for every smooth angular dependence, and eps_1 drops out at leading order (`eps1_drops`,
`bdy_vs_adm`, `free_symbols`). Decay needed (power counting, by hand, not computed): for eps - eps_0 = O(r^-d) with
d(eps) = O(r^-1-d), the linear term is zero for every d by the antisymmetry (the curved-space factors psi^2 - 1 times it are
O(r^-d)), the quadratic terms Re<eta, (delta_ik - g^i g^k) n_i d_k eta> r^2 are O(r^(1-2d)), and the cross terms with the spin
connection and K are O(r^-d). So d > 1/2 is enough, and the quadratic term is the one that can survive below it. For a power law,
d(eps) in L^2 at infinity needs 2 (1 + d) > 3, the same d > 1/2. With d = 0, eta = O(1) angular, the flux grows like r
(control `eps-O1`). The coefficients of psi_2, psi_3, K3 and eps_2 drop out by power counting, and `free_symbols` shows it.

Checks, with the controls that must fail them in CONTROLS. Planted first: the gammas (`clifford`, `dagger`), the derivative
and the sphere moments against sympy's own derivative and integral (`planted_derivation`, `planted_sphere`), the spin
connection against Geo's (`omega_vs_conf`), the flux formula (`div_formula`), Schwarzschild: E = m and B = m|eps_0|^2 with
the finite-r series m(1 + m/2r)^3 and m(1 + m/2r)|eps_0|^2 (`planted_E`, `planted_B_schwarzschild`), a dipole in psi
(`planted_dipole`), P_ADM = P (1 + m/2r)^-2 and the Bowen-York flux (`planted_P`, `planted_BY`), the eps_0 with xi_0 = +P/|P|
that gives B = m - |P| (`aligned_eps0`). Then the general data: B = E |eps_0|^2 - P.xi_0 (`bdy_vs_adm`, with E, P and xi_0
computed from the same data), with the trace term (`bdy_vs_adm_T`), no term of order r^2 or r pointwise or integrated
(`no_divergence`, `limit_exists`), and the eps_1 terms (`eps1_drops`, `eps1_hand`, `free_symbols`). Last, V against the
identity itself on the same data, at three exact points with exact rational values of every symbol (`identity_points`), whose
left side uses weitzenbock_conf's Geo and constraints and not this script's spin connection, and the other sign of nabla~
(`flip_consistent`): with nabla_i - (1/2) K_ij g^j g^0 the identity holds with +J.xi and V's K term +(1/2) P_ij xi^j, and the
boundary term is E + P.xi_0. The controls `Vsign` and `sign` both change V's K term from -(1/2) P_ij xi^j to +(1/2) P_ij xi^j
and nothing else in V; `Vsign` keeps nabla~ and the identity as stated, `sign` also takes nabla~ with the other sign and
leaves the identity's -J.xi, so each fails `identity_points` and each makes the P.xi_0 term disagree (B = E + P.xi_0).

What is shared with weitzenbock_conf.py, loaded read-only with its argv hidden: `note_gammas`, `inner` and `re` (so <a, b> and
the conjugation, and the controls `gamma` and `no-conj` act on them), `Geo` and `constraints` (the left side of
`identity_points`, `omega_vs_conf`), `christoffel`, `fn`, `sym_K`, `is_zero`, `ETA`, `X3` and `I3`. V's code, the series,
the jets, the ADM integrals and the sphere moments are new here, typed from the formulas above. V and the ADM integrals share
the derivative D_i and the sphere moments, which `planted_derivation` and `planted_sphere` check against independent code, and
the hand values of the Schwarzschild and Bowen-York cases do not use the code. The 4-metric route of #87, which ties mu and J
to G_ab, is not repeated.

What the code proves (*checked symbolically*): for the data above, with m, P, eps_0 and every angular coefficient arbitrary
symbols, the flux has no term of order r^2 or r, and its order r^0 term is, as a polynomial identity in the symbols, E
|eps_0|^2 - P.xi_0, with E and P the ADM integrals of the same data; V satisfies #87's identity on these data at three exact
points, in exact arithmetic (for all psi, K and eps the identity is #87's, *checked symbolically*). Assumed, not
computed: that data with remainders (psi - 1 - m/2r - psi_2 s^2 = O(r^-3) and so on) contribute only terms of higher order,
by power counting with each derivative gaining a power of 1/r; the decay statement above; and the integrated identity, whose
D~ eps = 0 solution and existence (item 5) are not touched. Not covered: a frame or slice other than psi^4 delta, K-hat other
than Bowen-York plus the extensions, and angular dependence of eps_1 that is not polynomial.

    python3 notes/warp/weitzenbock_boundary.py                  checks, then every control
    python3 notes/warp/weitzenbock_boundary.py --no-controls    checks only
    python3 notes/warp/weitzenbock_boundary.py --mutate=NAME    one control (its target checks must fail, exit 1)
    python3 notes/warp/weitzenbock_boundary.py --no-controls --only=NAME,...   (or --skip=NAME,...) some checks, for staged runs

Exit status is 1 if a check fails or a control does not fire.
"""
import importlib.util
import pathlib
import subprocess
import sys
import time

import sympy as sp

OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate", "")
ONLY = set(OPTS["only"].split(",")) if "only" in OPTS else None
SKIP = set(OPTS.get("skip", "").split(","))
HERE = pathlib.Path(__file__).resolve().parent
RESULTS = []
T0 = time.time()
SPHERE_LIMIT = 4 * 3600  # wall seconds allowed to any single sphere integral

# Every control is a mutation of the code that must make its target checks fail. The target lists are written before the first
# control run and are not changed after it.
VSPIN_T = ["identity_points", "planted_B_schwarzschild", "planted_dipole", "bdy_vs_adm", "eps1_drops", "eps1_hand", "free_symbols"]
CONTROLS = {
    "gamma": ("scales g^1 by 11/10", ["clifford", "identity_points"]),
    "no-conj": ("takes <a, b> = a^T b, without the complex conjugate", ["dagger", "identity_points"]),
    "deriv": ("drops the radial part of d_i, taking D_i f = s d_ni f", ["planted_derivation", "planted_E", "planted_B_schwarzschild"]),
    "sphere": ("takes the moment as a!! b!! c!!/(a+b+c+1)!!, not (a-1)!! (b-1)!! (c-1)!!/(a+b+c+1)!!",
               ["planted_sphere", "planted_BY", "planted_P"]),
    "omega-sign": ("takes the spin connection with -(1/4) omega_jk g^j g^k", ["omega_vs_conf", "planted_B_schwarzschild", "identity_points"]),
    "weight": ("takes the flux and the divergence without the weight psi^4 (and psi^-6)",
               ["div_formula", "planted_B_schwarzschild", "identity_points"]),
    "Vsign": ("flips V's K term, +(1/2) P_ij xi^j, with nabla~ and the identity as stated",
              ["identity_points", "planted_BY", "bdy_vs_adm", "bdy_vs_adm_T", "flip_consistent", "aligned_eps0"]),
    "sign": ("takes nabla~ with the other sign of its K term, nabla_i - (1/2) K_ij g^j g^0, and V's K term derived from it, "
             "-(1/2) sK P_ij xi^j with sK = -1, so +(1/2) P_ij xi^j; the identity keeps -J.xi",
             ["identity_points", "planted_BY", "bdy_vs_adm", "bdy_vs_adm_T", "aligned_eps0"]),
    "Jsign": ("writes +J.xi in the pointwise identity", ["identity_points"]),
    "Kdiv": ("drops V's K term", ["identity_points", "planted_BY", "bdy_vs_adm", "bdy_vs_adm_T", "flip_consistent", "aligned_eps0"]),
    "Vspin": ("drops -g^i D eps from V", VSPIN_T),
    "Dhalf": ("takes -(1/2) g^i D eps in V, so g^i g^k is no longer antisymmetrized against delta_ik", VSPIN_T),
    "trace": ("takes P_ij = K_ij in V, without -delta_ij tr K", ["identity_points", "bdy_vs_adm_T"]),
    "adm-E": ("drops -d_i g_jj from the ADM energy integrand", ["planted_E", "planted_dipole", "bdy_vs_adm"]),
    "adm-P-sign": ("takes the ADM momentum with the other sign", ["planted_P", "bdy_vs_adm", "bdy_vs_adm_T"]),
    "adm-frame": ("takes the frame components of K as its coordinate components in the ADM momentum", ["planted_P"]),
    "adm-trace": ("drops -delta_ij tr K from the ADM momentum", ["bdy_vs_adm_T"]),
    "eps-O1": ("adds an O(1) angular term eta(n) to eps, so eps - eps_0 is not O(1/r)", ["no_divergence", "limit_exists"]),
}


def load_conf():
    """weitzenbock_conf.py, loaded read-only with its argv hidden, as that script loads conventions.py."""
    spec = importlib.util.spec_from_file_location("weitzenbock_conf", HERE / "weitzenbock_conf.py")
    mod = importlib.util.module_from_spec(spec)
    saved, sys.argv = sys.argv, sys.argv[:1]
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.argv = saved
    return mod


WC = load_conf()
WC.MUT = MUT if MUT in ("no-conj", "gamma") else ""  # the two mutations of its helpers that this script passes through

nx, ny, nz = N3 = sp.symbols("nx ny nz", real=True)
s = sp.Symbol("s", positive=True)
NMAX = 3  # series in s are kept to s^3: V to s^3, so the flux (s^-2 times V) to s^1


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def todo(name):
    """Whether the check `name` runs: not when staged away, and under a mutation only if the mutation targets it."""
    if ONLY is not None and name not in ONLY:
        return False
    return name not in SKIP and (not MUT or name in CONTROLS[MUT][1])


def zero(e):
    """e == 0 after expansion, for a scalar or every entry of a matrix"""
    if isinstance(e, sp.MatrixBase):
        return all(sp.expand(x) == 0 for x in e)
    return sp.expand(e) == 0


def gammas():
    return WC.note_gammas()


# ---- Sphere averages, exact --------------------------------------------------------------------------------------------------
def dfact(k):
    return 1 if k <= 0 else k * dfact(k - 2)


def moment(a, b, c):
    """(1/4 pi) of the integral over the unit sphere of nx^a ny^b nz^c: 0 unless all even, else (a-1)!!(b-1)!!(c-1)!!/(a+b+c+1)!!"""
    if a % 2 or b % 2 or c % 2:
        return sp.Integer(0)
    if MUT == "sphere":
        return sp.Rational(dfact(a + 1) * dfact(b + 1) * dfact(c + 1), dfact(a + b + c + 1))
    return sp.Rational(dfact(a - 1) * dfact(b - 1) * dfact(c - 1), dfact(a + b + c + 1))


def savg(e, t0=None):
    """(1/4 pi) times the integral over the unit sphere of the polynomial e(n), by exact moments."""
    t0 = time.time() if t0 is None else t0
    tot = 0
    for t in sp.Add.make_args(sp.expand(e)):
        if time.time() - t0 > SPHERE_LIMIT:
            sys.exit("sphere integral not finished in 4 h: stopped")
        c, mon = t.as_independent(*N3, as_Add=False)
        pw = mon.as_powers_dict()
        tot += c * moment(*(int(pw.get(v, 0)) for v in N3))
    return sp.expand(tot)


# ---- Series in s: split by power -------------------------------------------------------------------------------------------
def split_s(e):
    """dict k -> coefficient of s^k in the scalar e."""
    out = {}
    for t in sp.Add.make_args(sp.expand(e)):
        c, k = t.as_coeff_exponent(s)
        out[int(k)] = out.get(int(k), 0) + c
    return out


def split_vec(v):
    """dict k -> column vector of the coefficients of s^k of the entries of the column vector v."""
    parts = [split_s(x) for x in v]
    ks = sorted({k for p in parts for k in p})
    return {k: sp.Matrix([p.get(k, 0) for p in parts]) for k in ks}


# ---- Representations: a field is a value (scalar or matrix) that the representation can multiply, differentiate and invert -----
class Fld:
    """A scalar or matrix valued field in one representation `rep`, with the operators the V code needs."""
    __slots__ = ("rep", "v", "ord")

    def __init__(self, rep, v, order=None):
        self.rep, self.v, self.ord = rep, v, order

    def lift(self, o):
        return o if isinstance(o, Fld) else Fld(self.rep, o, None)

    def __add__(self, o):
        o = self.lift(o)
        return self.rep.mk(self.v + o.v, self.rep.ordmin(self, o))

    __radd__ = __add__

    def __neg__(self):
        return Fld(self.rep, -self.v, self.ord)

    def __sub__(self, o):
        return self + (-self.lift(o))

    def __rsub__(self, o):
        return self.lift(o) - self

    def __mul__(self, o):
        return self.rep.mul(self, self.lift(o))

    def __rmul__(self, o):
        return self.rep.mul(self.lift(o), self)

    def d(self, i):
        return self.rep.deriv(self, i)

    def at(self, i, j):
        return Fld(self.rep, self.v[i, j], self.ord)

    def col(self, i):
        return Fld(self.rep, self.v[i, 0], self.ord)


def fsum(items, zero_):
    out = zero_
    for it in items:
        out = out + it
    return out


class Series:
    """Representation 1: series in s = 1/r at the direction n, coefficients polynomial in n_i and the symbols. d_i is the
    Cartesian derivative at x = n/s on the sphere: D_i f = s (d_ni f - n_i (n.d_n f + s d_s f)), which raises the power of s."""
    name = "series"

    def elementwise(self, f, v):
        return v.applyfunc(f) if isinstance(v, sp.MatrixBase) else f(v)

    def trunc(self, v, order=NMAX):
        def one(e):
            return sp.Add(*[t for t in sp.Add.make_args(sp.expand(e)) if t.as_coeff_exponent(s)[1] <= order])
        return self.elementwise(one, v)

    def minorder(self, v):
        ents = list(v) if isinstance(v, sp.MatrixBase) else [v]
        ks = [min(split_s(e)) for e in ents if sp.expand(e) != 0]
        return min(ks) if ks else 10**6

    def ordmin(self, a, b):
        return None

    def mk(self, v, order=None):
        return Fld(self, v, None)

    def const(self, M):
        return Fld(self, M, None)

    def mul(self, a, b):
        ma, mb = self.minorder(a.v), self.minorder(b.v)
        av = self.trunc(a.v, NMAX - min(mb, NMAX + 1))
        bv = self.trunc(b.v, NMAX - min(ma, NMAX + 1))
        return Fld(self, self.trunc(av * bv), None)

    def deriv(self, a, i):
        def one(f):
            radial = sum(N3[j] * sp.diff(f, N3[j]) for j in range(3)) + s * sp.diff(f, s)
            if MUT == "deriv":
                return sp.expand(s * sp.diff(f, N3[i]))
            return sp.expand(s * (sp.diff(f, N3[i]) - N3[i] * radial))
        return Fld(self, self.elementwise(one, a.v), None)

    def ipow(self, psi, p):
        """psi^-p as a series, psi = 1 + u with u = O(s)."""
        u = psi.v - 1
        out, term = sp.Integer(1), sp.Integer(1)
        for k in range(1, NMAX + 1):
            term = self.trunc(term * u)
            out = out + sp.binomial(-p, k) * term
        return Fld(self, sp.expand(out), None)

    def inner(self, a, b):
        """<a, b> of column vectors, to s^NMAX, using weitzenbock_conf.inner on the coefficient vectors."""
        pa, pb = split_vec(a.v), split_vec(b.v)
        tot = 0
        for i, x in pa.items():
            for j, y in pb.items():
                if i + j <= NMAX:
                    tot += s ** (i + j) * WC.inner(x, y)
        return Fld(self, sp.expand(tot), None)


SER = Series()


class Sym:
    """Representation 2: exact, generic functions of (x, y, z), no truncation. For checks on arbitrary psi and V."""
    name = "symbolic"

    def ordmin(self, a, b):
        return None

    def mk(self, v, order=None):
        return Fld(self, v, None)

    def const(self, M):
        return Fld(self, M, None)

    def mul(self, a, b):
        return Fld(self, a.v * b.v, None)

    def deriv(self, a, i):
        return Fld(self, a.v.diff(WC.X3[i]), None)

    def ipow(self, psi, p):
        return Fld(self, psi.v ** (-p), None)

    def inner(self, a, b):
        return Fld(self, WC.inner(a.v, b.v), None)


SYM = Sym()


class Jet:
    """Representation 3: Taylor jets at the exact point `pt`, to total degree 2 in the offsets (dx, dy, dz), exact rational
    arithmetic. A field carries its valid order, which a derivative lowers by one."""
    name = "jet"
    DV = sp.symbols("dx dy dz", real=True)

    def __init__(self, pt, order=2):
        self.pt, self.order = pt, order

    def tdeg(self, t):
        pw = t.as_powers_dict()
        return sum(int(pw.get(v, 0)) for v in self.DV)

    def trunc(self, v, order):
        if order is None:
            return v

        def one(e):
            return sp.Add(*[t for t in sp.Add.make_args(sp.expand(e)) if self.tdeg(t) <= order])
        return v.applyfunc(one) if isinstance(v, sp.MatrixBase) else one(v)

    def ordmin(self, a, b):
        os_ = [o for o in (a.ord, b.ord) if o is not None]
        return min(os_) if os_ else None

    def mk(self, v, order=None):
        return Fld(self, v, order)

    def const(self, M):
        return Fld(self, M, None)

    def mul(self, a, b):
        o = self.ordmin(a, b)
        return Fld(self, self.trunc(a.v * b.v, o), o)

    def deriv(self, a, i):
        o = None if a.ord is None else a.ord - 1
        f = lambda e: sp.expand(sp.diff(e, self.DV[i]))
        v = a.v.applyfunc(f) if isinstance(a.v, sp.MatrixBase) else f(a.v)
        return Fld(self, v, o)

    def ipow(self, psi, p):
        c0 = psi.v.subs({v: 0 for v in self.DV})
        u = sp.expand(psi.v / c0 - 1)
        out, term = sp.Integer(1), sp.Integer(1)
        for k in range(1, psi.ord + 1):
            term = self.trunc(sp.expand(term * u), psi.ord)
            out = out + sp.binomial(-p, k) * term
        return Fld(self, sp.expand(c0 ** (-p) * out), psi.ord)

    def inner(self, a, b):
        o = self.ordmin(a, b)
        return Fld(self, self.trunc(WC.inner(a.v, b.v), o), o)

    def field(self, expr):
        """The jet of an expression in (n, s) at the point, as a field of valid order `order`."""
        px, py, pz = self.pt
        rr = sp.sqrt((px + self.DV[0])**2 + (py + self.DV[1])**2 + (pz + self.DV[2])**2)
        jn = [self.taylor((self.pt[i] + self.DV[i]) / rr) for i in range(3)]
        js = self.taylor(1 / rr)

        def one(e):
            e = sp.expand(e)
            if e == 0:
                return sp.Integer(0)
            P = sp.Poly(e, *N3, s)
            out = 0
            for mon, c in P.terms():
                t = c
                for var, ex in zip(jn + [js], mon):
                    for _ in range(ex):
                        t = self.trunc(sp.expand(t * var), self.order)
                out += t
            return sp.expand(out)
        v = expr.applyfunc(one) if isinstance(expr, sp.MatrixBase) else one(expr)
        return Fld(self, v, self.order)

    def taylor(self, f):
        """Taylor polynomial of f(dx, dy, dz) at 0, to total degree `order`, exact."""
        zero0 = {v: 0 for v in self.DV}
        out = 0
        import itertools
        for a in itertools.product(range(self.order + 1), repeat=3):
            if sum(a) > self.order:
                continue
            g = f
            for v, k in zip(self.DV, a):
                if k:
                    g = sp.diff(g, v, k)
            val = sp.simplify(g.subs(zero0))
            fact = 1
            for k in a:
                fact *= sp.factorial(k)
            out += val / fact * self.DV[0]**a[0] * self.DV[1]**a[1] * self.DV[2]**a[2]
        return sp.expand(out)


# ---- V of #87 and the surface integrals ------------------------------------------------------------------------------------
def spin_conn(psi):
    """Om_i = (1/4) omega_jk(e_i) g^j g^k = (1/4) psi^-2 (g^i phislash - phislash g^i), phi_k = 2 d_k psi / psi, which is
    omega_jk(e_i) = psi^-2 (delta_ij phi_k - delta_ik phi_j) of weitzenbock_conf.omega_hand, as a list of matrix fields."""
    rep = psi.rep
    g0, gj = gammas()
    G = [rep.const(g) for g in gj]
    ipsi, ipsi2 = rep.ipow(psi, 1), rep.ipow(psi, 2)
    phi = [2 * (psi.d(k) * ipsi) for k in range(3)]
    ps = fsum([phi[k] * G[k] for k in range(3)], rep.const(sp.zeros(4)))
    sg = -1 if MUT == "omega-sign" else 1
    return [sg * (ipsi2 * (G[i] * ps - ps * G[i])) * sp.Rational(1, 4) for i in range(3)]


def Vvec(psi, Kf, eps, sK=None):
    """V^i = Re<eps, nabla_i eps - g^i D eps> - (1/2) P_ij xi^j, P_ij = K_ij - delta_ij tr K, D = g^i nabla_i, the frame
    components of #87's V on psi^4 delta in e_i = psi^-2 d_i, as scalar fields. Kf is the 3x3 matrix field of frame components."""
    rep = psi.rep
    g0, gj = gammas()
    G0, G = rep.const(g0), [rep.const(g) for g in gj]
    ipsi2 = rep.ipow(psi, 2)
    Om = spin_conn(psi)
    nab = [ipsi2 * eps.d(i) + Om[i] * eps for i in range(3)]
    zero4 = rep.const(sp.zeros(4, 1))
    Deps = fsum([G[i] * nab[i] for i in range(3)], zero4)
    cD = {"Vspin": 0, "Dhalf": sp.Rational(1, 2)}.get(MUT, 1)
    xi = [rep.inner(eps, G0 * (G[j] * eps)) for j in range(3)]
    trK = Kf.at(0, 0) + Kf.at(1, 1) + Kf.at(2, 2)
    ctr = 0 if MUT == "trace" else 1
    if sK is None:  # nabla~'s K term is sK (1/2) K_ij g^j g^0, and V's K term is -(1/2) sK P_ij xi^j
        sK = -1 if MUT == "sign" else 1
    kc = {"Vsign": -1, "Kdiv": 0}.get(MUT, 1) * sK
    V = []
    for i in range(3):
        spin = nab[i] - cD * (G[i] * Deps)
        vs = WC.re(rep.inner(eps, spin).v)
        Pxi = fsum([(Kf.at(i, j) - (trK * ctr if i == j else 0)) * xi[j] for j in range(3)], rep.const(sp.Integer(0)))
        V.append(Fld(rep, vs, None) - sp.Rational(1, 2) * kc * Pxi)
    return V


def wfac(psi):
    """The weight psi^4 of the flux through a coordinate sphere, psi^-6 d_a (psi^4 V^a) being the divergence."""
    return psi * psi * psi * psi if MUT != "weight" else psi.rep.const(sp.Integer(1))


def flux_series(psi, Kf, eps, sK=None):
    """dict j -> (1/4 pi) of the integral of psi^4 V^a n_a r^2 dOmega at order r^-j, i.e. the coefficient of s^j, for j = -2..1."""
    V = Vvec(psi, Kf, eps, sK)
    nV = fsum([V[a] * Fld(SER, N3[a]) for a in range(3)], Fld(SER, sp.Integer(0)))
    F = wfac(psi) * nV
    return {k - 2: savg(c) for k, c in sorted(split_s(F.v).items()) if k <= NMAX}, F


def adm_E_series(psi):
    """E at order r^-j: (1/16 pi) of the integral of (d_j g_ij - d_i g_jj) n_i r^2 dOmega for g = psi^4 delta, j = -2..2"""
    p4 = psi * psi * psi * psi
    z = Fld(SER, sp.Integer(0))
    f = []
    for i in range(3):
        dj = p4.d(i)                                  # d_j g_ij, only j = i
        di_tr = 0 if MUT == "adm-E" else 3 * p4.d(i)  # d_i g_jj
        f.append(dj - di_tr)
    F = fsum([f[i] * Fld(SER, N3[i]) for i in range(3)], z)
    return {k - 2: sp.Rational(1, 4) * savg(c) for k, c in sorted(split_s(F.v).items()) if k <= NMAX}


def adm_P_series(psi, Kf):
    """P_j at order r^-k: (1/8 pi) of the integral of (K_ij - delta_ij tr K) n_i r^2 dOmega with the coordinate components
    K_ij = psi^4 (frame components), tr K = delta^ij K_ij. A list of three dicts."""
    p4 = psi * psi * psi * psi
    Kc = [[p4 * Kf.at(i, j) if MUT != "adm-frame" else Kf.at(i, j) for j in range(3)] for i in range(3)]
    trc = 0 if MUT == "adm-trace" else Kc[0][0] + Kc[1][1] + Kc[2][2]
    sg = -1 if MUT == "adm-P-sign" else 1
    out = []
    z = Fld(SER, sp.Integer(0))
    for j in range(3):
        F = fsum([(Kc[i][j] - (trc if i == j else 0)) * Fld(SER, N3[i]) for i in range(3)], z)
        out.append({k - 2: sg * sp.Rational(1, 2) * savg(c) for k, c in sorted(split_s(F.v).items()) if k <= NMAX})
    return out


# ---- The data ---------------------------------------------------------------------------------------------------------------
SYMS = {}


def sym(name):
    return SYMS.setdefault(name, sp.Symbol(name, real=True))


def monos(dmax):
    return [(a, b, d - a - b) for d in range(dmax + 1) for a in range(d + 1) for b in range(d - a + 1)]


def poly_n(tag, dmax, cplx=False):
    """A polynomial in n_i of degree <= dmax with independent symbolic coefficients (complex: real and imaginary parts)."""
    out = 0
    for a, b, c in monos(dmax):
        co = sym(f"{tag}_{a}{b}{c}") + (sp.I * sym(f"{tag}i_{a}{b}{c}") if cplx else 0)
        out += co * nx**a * ny**b * nz**c
    return out


def data(psi2=None, P=False, K3=False, T=False, eps1=False, eps2=False, eta=False):
    """Expressions in (n, s, symbols): psi, the Bowen-York K-hat plus optional terms, eps. Returns a dict.
    psi = 1 + (m/2) s + psi_2(n) s^2 + psi_3(n) s^3, psi_2 = None (none), 'dipole' (d.n) or 'general' (degree <= 2).
    K-hat_ij = (3/2) s^2 [P_i n_j + P_j n_i - (delta_ij - n_i n_j) P.n] + s^2 T_ij(n) + s^3 K3_ij(n), K = psi^-6 K-hat (frame).
    eps = eps_0 + eta(n) + eps_1(n) s + eps_2(n) s^2, eps_0 constant, eps_1 and eps_2 with independent coefficients, eta O(1)."""
    m = sym("m")
    n = sp.Matrix(N3)
    psi = 1 + m * s / 2
    if psi2 == "dipole":
        psi += sum(sym(f"d{c}") * n[i] for i, c in enumerate("xyz")) * s**2
    elif psi2 == "general":
        psi += poly_n("p2", 2) * s**2 + poly_n("p3", 1) * s**3
    Pv = sp.Matrix([sym(f"P{c}") for c in "xyz"]) if P else sp.zeros(3, 1)
    Pn = (Pv.T * n)[0]
    Kh = sp.Matrix(3, 3, lambda i, j: sp.Rational(3, 2) * s**2 * (Pv[i] * n[j] + Pv[j] * n[i] - ((1 if i == j else 0) - n[i] * n[j]) * Pn))
    if T:
        Kh += sp.Matrix(3, 3, lambda i, j: s**2 * poly_n(f"T{min(i, j)}{max(i, j)}", 1))
    if K3:
        Kh += sp.Matrix(3, 3, lambda i, j: s**3 * poly_n(f"k{min(i, j)}{max(i, j)}", 1))
    e0 = sp.Matrix([sym(f"u{a}") + sp.I * sym(f"v{a}") for a in range(4)])
    eps = e0
    if eta:
        eps = eps + sp.Matrix([poly_n(f"h{a}", 1, True) for a in range(4)])
    if eps1:
        eps = eps + s * sp.Matrix([poly_n(f"f{a}", 2, True) for a in range(4)])
    if eps2:
        eps = eps + s**2 * sp.Matrix([poly_n(f"g{a}", 1, True) for a in range(4)])
    return dict(psi=sp.expand(psi), Kh=Kh.applyfunc(sp.expand), eps=eps.applyfunc(sp.expand), e0=e0, P=Pv, m=m)


def fields(rep, D, vals=None):
    """The fields psi, K (frame, psi^-6 K-hat) and eps of a data dict in the representation rep."""
    if rep is SER:
        mk = lambda e: Fld(rep, e, None)
    else:
        mk = lambda e: rep.field(e.subs(vals))
    psi, Kh, eps = mk(D["psi"]), mk(D["Kh"]), mk(D["eps"])
    return psi, rep.ipow(psi, 6) * Kh, eps


def xi0(D):
    g0, gj = gammas()
    e0 = D["e0"]
    return [WC.inner(e0, g0 * gj[j] * e0) for j in range(3)]


def e2(D):
    return WC.inner(D["e0"], D["e0"])


# ---- Runs of the series code, cached ---------------------------------------------------------------------------------------
_CACHE = {}


def series_run(tag, **kw):
    """flux, E and P series of one data family, once per process."""
    if tag not in _CACHE:
        t0 = time.time()
        D = data(**kw)
        psi, Kf, eps = fields(SER, D)
        flux, F = flux_series(psi, Kf, eps)
        _CACHE[tag] = dict(D=D, flux=flux, F=F, E=adm_E_series(psi), P=adm_P_series(psi, Kf), xi=xi0(D), e2=e2(D))
        _CACHE[tag]["secs"] = time.time() - t0
        nsym = len({v for e in (D["psi"], *D["Kh"], *D["eps"]) for v in e.free_symbols} - set(N3) - {s})
        print(f"    [{tag}: {nsym} symbols, flux, E and P series in {time.time() - t0:.1f} s]", flush=True)
    return _CACHE[tag]


def hand_series(expr, kmax=1):
    """coefficients of s^0 .. s^kmax of expr(s)"""
    ser = sp.series(expr, s, 0, kmax + 1).removeO()
    c = split_s(ser)
    return {j: c.get(j, 0) for j in range(kmax + 1)}


def agree(got, want, js):
    """got[j] == want[j] for each j in js (a missing entry is 0)"""
    return all(zero(got.get(j, 0) - want.get(j, 0)) for j in js)


# ---- 1. Conventions --------------------------------------------------------------------------------------------------------
def section_conventions():
    print("\n== 1. The gammas and the inner product of weitzenbock_conf.py ==")
    g0, gj = gammas()
    G = [g0] + gj
    if todo("clifford"):
        anti = all(WC.is_zero(G[a] * G[b] + G[b] * G[a] - 2 * WC.ETA[a, b] * sp.eye(4)) for a in range(4) for b in range(4))
        herm = g0.H == -g0 and all(g.H == g for g in gj)
        check("clifford", anti and herm, "{g^a, g^b} = 2 eta^ab, g^0 anti-Hermitian, g^j Hermitian: the note's gammas")
    if todo("dagger"):
        e = sp.Matrix([sym(f"u{a}") + sp.I * sym(f"v{a}") for a in range(4)])
        norm = zero(WC.inner(e, e) - sum(sym(f"u{a}")**2 + sym(f"v{a}")**2 for a in range(4)))
        real = all(zero(WC.inner(e, g0 * g * e) - WC.conj(WC.inner(e, g0 * g * e))) for g in gj)
        check("dagger", norm and real, "<a, b> = a^dagger b: |eps|^2 = sum u^2 + v^2 and xi^j = <eps, g^0 g^j eps> real")


# ---- 2. The machinery ------------------------------------------------------------------------------------------------------
def section_machinery():
    print("\n== 2. The machinery: derivative, sphere moments, spin connection, flux formula ==")
    if todo("planted_derivation"):
        # D_i f against the Cartesian derivative of f(x/r, y/r, z/r, 1/r), exactly, at points with rational r; once and twice
        x, y, z = WC.X3
        rr = sp.sqrt(x**2 + y**2 + z**2)
        fs = [nx * ny * s**2 + 3 * nz * s + ny**2 * s**3, nx**3 * s**2 - 2 * nz * ny * s**3 + s, (nx + 2 * ny * nz) * s**2]
        pts = [(1, 2, 2), (2, 3, 6), (-4, 4, 7)]
        ok = True
        for f in fs:
            fc = f.subs({nx: x / rr, ny: y / rr, nz: z / rr, s: 1 / rr}, simultaneous=True)
            for i in range(3):
                for j in range(-1, 3):  # j = -1: first derivative only
                    got = SER.deriv(Fld(SER, f), i).v if j < 0 else SER.deriv(SER.deriv(Fld(SER, f), i), j).v
                    want = sp.diff(fc, WC.X3[i]) if j < 0 else sp.diff(fc, WC.X3[i], WC.X3[j])
                    for p in pts:
                        r0 = sp.sqrt(sum(c * c for c in p))
                        a = got.subs({nx: sp.Rational(p[0]) / r0, ny: sp.Rational(p[1]) / r0, nz: sp.Rational(p[2]) / r0, s: 1 / r0})
                        b = want.subs(dict(zip(WC.X3, p)))
                        ok = ok and sp.simplify(a - b) == 0
        check("planted_derivation", ok, "D_i f = s (d_ni f - n_i (n.d_n f + s d_s f)) equals d/dx_i of f(x/r, y/r, z/r, 1/r), once and "
              "twice, for three polynomials at three exact points on the sphere")
    if todo("planted_sphere"):
        th, ph = sp.symbols("theta phi", real=True)
        mons = [(0, 0, 0), (2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (2, 2, 0), (4, 0, 0), (0, 2, 2), (2, 2, 2), (6, 0, 0), (2, 4, 0),
                (3, 0, 0), (1, 2, 1), (4, 2, 0), (0, 0, 4)]
        ok = True
        for a, b, c in mons:
            f = (sp.sin(th) * sp.cos(ph))**a * (sp.sin(th) * sp.sin(ph))**b * sp.cos(th)**c * sp.sin(th)
            direct = sp.integrate(sp.integrate(f, (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)) / (4 * sp.pi)
            ok = ok and sp.simplify(direct - moment(a, b, c)) == 0
        check("planted_sphere", ok, f"(a-1)!!(b-1)!!(c-1)!!/(a+b+c+1)!! is (1/4 pi) of the integral of nx^a ny^b nz^c, against sympy's "
              f"integral in theta and phi for {len(mons)} monomials")
    if todo("omega_vs_conf"):
        psi_f = WC.fn("psi")
        Om = spin_conn(Fld(SYM, psi_f))
        geo = WC.Geo(psi_f)
        ok = all(zero(Om[i].v - geo.Om[i]) for i in range(3))
        check("omega_vs_conf", ok, "Om_i = (1/4) psi^-2 (g^i phislash - phislash g^i), phi = 2 d ln psi, equals weitzenbock_conf's "
              "Geo(psi).Om_i computed from the Christoffel symbols, for arbitrary psi(x, y, z)")
    if todo("div_formula"):
        psi_f = WC.fn("psi")
        W = [WC.fn(f"W{a}") for a in range(3)]  # frame components V^a
        w = psi_f**4 if MUT != "weight" else 1
        pre = psi_f**-6 if MUT != "weight" else 1
        mine = sp.expand(pre * sum(sp.diff(w * W[a], WC.X3[a]) for a in range(3)))
        g, gi = psi_f**4 * WC.I3, psi_f**-4 * WC.I3
        Chr = WC.christoffel(g, gi, WC.X3)
        U = [W[a] / psi_f**2 for a in range(3)]  # coordinate components of V = V^a e_a, e_a = psi^-2 d_a
        cov = sp.expand(sum(sp.diff(U[a], WC.X3[a]) for a in range(3)) + sum(Chr[a][a][k] * U[k] for a in range(3) for k in range(3)))
        check("div_formula", zero(mine - cov), "psi^-6 d_a(psi^4 V^a) equals the covariant divergence d_a U^a + Gamma^a_ak U^k of "
              "U = psi^-2 V, for arbitrary psi and V; so the flux of V through a coordinate sphere is the integral of "
              "psi^4 V^a n_a r^2 dOmega by Gauss's theorem in flat R^3")


# ---- 3. Planted boundary values --------------------------------------------------------------------------------------------
def section_planted():
    print("\n== 3. Planted boundary values: Schwarzschild, a dipole in psi, Bowen-York ==")
    m = sym("m")
    names = ("planted_E", "planted_B_schwarzschild", "planted_dipole", "planted_P", "planted_BY")
    if not any(todo(n) for n in names):
        return
    R = series_run("schwarzschild")
    e2_ = R["e2"]
    if todo("planted_E"):
        hand = hand_series(sp.expand(m * (1 + m * s / 2)**3), 1)
        check("planted_E", agree(R["E"], hand, (-2, -1, 0, 1)), "psi = 1 + m/2r, K = 0: the ADM integral of g = psi^4 delta is "
              f"m (1 + m/2r)^3, E = {R['E'].get(0)} at order 0 and {R['E'].get(1)} at order 1/r")
    if todo("planted_B_schwarzschild"):
        hand = {0: m * e2_, 1: m**2 * e2_ / 2}
        check("planted_B_schwarzschild", agree(R["flux"], hand, (-2, -1, 0, 1)), "psi = 1 + m/2r, K = 0, eps = eps_0: V^i = "
              "-2 psi^-3 d_i psi |eps_0|^2 by hand, flux m (1 + m/2r) |eps_0|^2, so the boundary term is m |eps_0|^2 = E |eps_0|^2")
    if todo("planted_dipole"):
        Rd = series_run("dipole", psi2="dipole")
        ok = zero(Rd["E"].get(0, 0) - m) and zero(Rd["flux"].get(0, 0) - m * Rd["e2"])
        ok = ok and all(zero(Rd["E"].get(j, 0)) and zero(Rd["flux"].get(j, 0)) for j in (-2, -1))
        check("planted_dipole", ok, "psi = 1 + m/2r + (d.n)/r^2: E and the boundary term are still m and m |eps_0|^2, with no "
              "terms of order r^2 or r")
    if todo("planted_P") or todo("planted_BY"):
        Rb = series_run("BY", P=True)
        Pv = [sym(f"P{c}") for c in "xyz"]
        if todo("planted_P"):
            ok = all(agree(Rb["P"][j], {0: Pv[j], 1: -m * Pv[j]}, (-2, -1, 0, 1)) for j in range(3))
            check("planted_P", ok, "psi = 1 + m/2r, K = psi^-6 Bowen-York: the ADM momentum integral is P_j (1 + m/2r)^-2, "
                  "P_j = " + str([Rb["P"][j].get(0) for j in range(3)]) + " at order 0, -m P_j at order 1/r")
        if todo("planted_BY"):
            xi = Rb["xi"]
            Pxi = sum(Pv[j] * xi[j] for j in range(3))
            hand = {0: m * Rb["e2"] - Pxi, 1: m**2 * Rb["e2"] / 2 + m * Pxi}
            check("planted_BY", agree(Rb["flux"], hand, (-2, -1, 0, 1)), "same data: the flux is m (1 + m/2r) |eps_0|^2 - "
                  "(1 + m/2r)^-2 P.xi_0 by hand (tr K-hat = 0 and the sphere mean of n_i K-hat_ij r^2 is 2 P_j), so the boundary term is m |eps_0|^2 - P.xi_0")


def section_aligned():
    """E - |P| is reached by the eps_0 with xi_0 = P/|P|, and E + |P| by xi_0 = -P/|P| (with the sign of nabla~ of the note)."""
    if not todo("aligned_eps0"):
        return
    print("\n== 3b. The eps_0 that gives E - |P| ==")
    R = series_run("BY", P=True)
    Pn = [1, 2, 2]  # |P| = 3, so that the projector and the spinors are exact
    g0, gj = gammas()
    al = [g0 * g for g in gj]  # alpha_j = g^0 g^j, Hermitian, alpha_i alpha_j + alpha_j alpha_i = 2 delta_ij
    Pa = sum((sp.Rational(Pn[j], 3) * al[j] for j in range(3)), sp.zeros(4))
    out = []
    ok = True
    for sgn, name in ((1, "xi_0 = +P/|P|"), (-1, "xi_0 = -P/|P|")):
        Pi = (sp.eye(4) + sgn * Pa) / 2  # projector onto the +1 (or -1) eigenspace of P.alpha/|P|
        c = next(j for j in range(4) if any(Pi[:, j]))
        e0 = (sp.sqrt(2) * Pi[:, c]).applyfunc(sp.expand)  # |eps_0|^2 = 2 Pi_cc = 1 for a projector column
        vals = {sym(f"u{a}"): sp.re(e0[a]) for a in range(4)}
        vals.update({sym(f"v{a}"): sp.im(e0[a]) for a in range(4)})
        vals.update({sym(f"P{c_}"): Pn[j] for j, c_ in enumerate("xyz")})
        xi = [sp.expand(x.subs(vals)) for x in R["xi"]]
        e2v = sp.expand(R["e2"].subs(vals))
        B = sp.expand(R["flux"].get(0, 0).subs(vals))
        want = sym("m") - sgn * 3
        good = e2v == 1 and xi == [sp.Rational(sgn * Pn[j], 3) for j in range(3)] and sp.expand(B - want) == 0
        ok = ok and good
        out.append(f"{name}: |eps_0|^2 = {e2v}, xi_0 = {xi}, boundary term {B}")
    check("aligned_eps0", ok, "P = (1, 2, 2), |P| = 3, eps_0 = sqrt(2) times a column of (1 +- P.alpha/|P|)/2: " + "; ".join(out)
          + ". So E - |P| is the boundary term of the eps_0 with xi_0 = +P/|P|, and E + |P| that of xi_0 = -P/|P|")


# ---- 4. The general boundary term ------------------------------------------------------------------------------------------
GEN = dict(psi2="general", P=True, K3=True, eps1=True, eps2=True)


def gen_run():
    # eps-O1 adds eta(n) to eps_0 with the Bowen-York K and no other angular terms: with all of the terms the flux series was not
    # finished after 2 h 14 min (1.2 GiB) in step 4's first run, and the divergence of the O(1) term does not depend on them
    kw = dict(P=True, eta=True) if MUT == "eps-O1" else GEN
    return series_run("general", **kw)


def eps1_hand_check(R):
    """The eps_1-dependent part of the order r^0 integrand against the hand form -Re<eps_0, nslash dslash_tan eps_1>, modulo n^2 = 1."""
    g0, gj = gammas()
    F0 = split_s(R["F"].v).get(2, 0)
    fs = {v for k, v in SYMS.items() if k.startswith("f")}
    part = sp.Add(*[t for t in sp.Add.make_args(sp.expand(F0)) if t.free_symbols & fs])
    e1 = sp.Matrix([poly_n(f"f{a}", 2, True) for a in range(4)])
    dtan = [e1.diff(N3[k]) - N3[k] * sum((N3[j] * e1.diff(N3[j]) for j in range(3)), sp.zeros(4, 1)) for k in range(3)]
    ns = sum((N3[i] * gj[i] for i in range(3)), sp.zeros(4))
    vec = sum((ns * gj[k] * dtan[k] for k in range(3)), sp.zeros(4, 1))
    hand = -WC.re(WC.inner(R["D"]["e0"], vec))
    diff = sp.expand(part - hand)
    rem = sp.rem(sp.Poly(diff, nz), sp.Poly(nz**2 + nx**2 + ny**2 - 1, nz)).as_expr()
    return sp.expand(rem) == 0, len(sp.Add.make_args(part))


def section_general():
    print("\n== 4. The boundary term of the general data against E |eps_0|^2 - P.xi_0 ==")
    names = ("bdy_vs_adm", "bdy_vs_adm_T", "no_divergence", "limit_exists", "eps1_drops", "eps1_hand", "free_symbols")
    if not any(todo(n) for n in names):
        return
    if any(todo(n) for n in ("bdy_vs_adm", "no_divergence", "limit_exists", "eps1_drops", "eps1_hand", "free_symbols")):
        R = gen_run()
        flux, E, P, xi, e2_ = R["flux"], R["E"], R["P"], R["xi"], R["e2"]
        eps1_syms = {v for k, v in SYMS.items() if k.startswith(("f", "g", "h"))}
        outside = {v for k, v in SYMS.items() if k[0] in "fghkpdT" and k not in ("m",)}
        F = split_s(R["F"].v)
        if todo("no_divergence"):
            nz_ = [k - 2 for k in (0, 1) if not zero(F.get(k, 0))]
            check("no_divergence", not nz_, "pointwise, psi^4 V^a n_a r^2 has no term of order r^2 or r: its first term is the order "
                  f"r^0 one, which has {len(sp.Add.make_args(F.get(2, 0)))} terms in n")
        if todo("limit_exists"):
            nt = [len(sp.Add.make_args(sp.expand(flux.get(j, 0)))) if not zero(flux.get(j, 0)) else 0 for j in (-2, -1)]
            check("limit_exists", nt == [0, 0], "the sphere integrals of the r^2 and r terms vanish (they have "
                  f"{nt[0]} and {nt[1]} terms): the limit r -> infinity of the boundary term exists")
        if todo("bdy_vs_adm"):
            Pdot = sum(P[j].get(0, 0) * xi[j] for j in range(3))
            want = E.get(0, 0) * e2_ - Pdot
            res = sp.expand(flux.get(0, 0) - want)
            check("bdy_vs_adm", res == 0, f"B = (1/4 pi) lim flux of V = E |eps_0|^2 - P.xi_0 with E = {E.get(0, 0)}, "
                  f"P_ADM = {[P[j].get(0, 0) for j in range(3)]}, xi_0^j = <eps_0, g^0 g^j eps_0> from the same data; "
                  f"{len(sp.Add.make_args(sp.expand(flux.get(0, 0))))} terms in B, {len(sp.Add.make_args(res)) if res != 0 else 0} in B - (E |eps_0|^2 - P.xi_0)")
        if todo("eps1_drops"):
            F0 = F.get(2, 0)
            dep_pointwise = bool(F0.free_symbols & eps1_syms)
            dep_integral = bool(sp.expand(flux.get(0, 0)).free_symbols & eps1_syms)
            check("eps1_drops", dep_pointwise and not dep_integral, "the eps_1 and eps_2 coefficients enter the order r^0 integrand "
                  f"(pointwise it has {len(sp.Add.make_args(F0))} terms) but not its sphere integral: the spin part's eps_1 term is "
                  "-Re<eps_0, nslash dslash_tan eps_1> (eps1_hand), a rotation derivative with antisymmetric g^i g^k whose sphere integral is 0")
        if todo("eps1_hand"):
            okh, nth = eps1_hand_check(R)
            check("eps1_hand", okh, f"the {nth} terms of the order r^0 integrand that contain eps_1 equal -Re<eps_0, nslash dslash_tan eps_1> "
                  "modulo n^2 = 1, with dslash_tan = g^k (d_nk - n_k n.d_n): nslash dslash_tan = sum_{i<k} g^i g^k (n_i d_k - n_k d_i) is a "
                  "sum of rotation generators with antisymmetric g^i g^k, whose sphere integral is 0")
        if todo("free_symbols"):
            ok = not (sp.expand(flux.get(0, 0)).free_symbols & outside)
            check("free_symbols", ok, "the boundary term depends on m, P and eps_0 only: not on the angular parts of psi, the r^-3 terms "
                  "of K-hat, or eps_1 and eps_2 (" + str(len(outside)) + " symbols)")
    if todo("bdy_vs_adm_T"):
        R = series_run("T", psi2="dipole", P=True, T=True)
        flux, E, P, xi, e2_ = R["flux"], R["E"], R["P"], R["xi"], R["e2"]
        Pdot = sum(P[j].get(0, 0) * xi[j] for j in range(3))
        res = sp.expand(flux.get(0, 0) - (E.get(0, 0) * e2_ - Pdot))
        Pv = [sym(f"P{c}") for c in "xyz"]
        shifted = any(not zero(P[j].get(0, 0) - Pv[j]) for j in range(3))
        check("bdy_vs_adm_T", res == 0 and shifted, "with an added r^-2 term T_ij(n) in K-hat that has a trace (not Bowen-York), "
              f"the boundary term is still E |eps_0|^2 - P.xi_0, P_ADM = P + {[sp.expand(P[j].get(0, 0) - Pv[j]) for j in range(3)][:1]}... "
              "so the trace term of the ADM integral and of V is exercised")


# ---- 5. The identity on these data, at exact points ------------------------------------------------------------------------
def values():
    """Exact rational values for every symbol of the fullest data family, fixed by the sorted names."""
    D = data(psi2="general", P=True, K3=True, T=True, eps1=True, eps2=True, eta=True)
    syms = sorted({v for e in (D["psi"], *D["Kh"], *D["eps"]) for v in e.free_symbols} - set(N3) - {s}, key=str)
    return {v: (sp.Rational(3, 5) if v.name == "m" else sp.Rational(((5 * i + 2) % 11) - 5, 10)) for i, v in enumerate(syms)}


def atom_values(exprs, jets):
    """Values at the point of the generic jets (psi, K_ij and their derivatives) in exprs, from the Taylor polynomials `jets`."""
    from sympy.core.function import AppliedUndef
    DV = Jet.DV
    zero0 = {v: 0 for v in DV}
    atoms = set()
    for e in exprs:
        atoms |= e.atoms(sp.Derivative, AppliedUndef)
    out = {}
    for a in atoms:
        if isinstance(a, sp.Derivative):
            poly, vc = jets[a.expr.func.__name__], a.variable_count
            args = [DV[WC.X3.index(v)] for v, c in vc for _ in range(c)]
            out[a] = sp.diff(poly, *args).subs(zero0)
        else:
            out[a] = jets[a.func.__name__].subs(zero0)
    return out


def identity_residual(pt, vals, sK, sJ):
    """lhs - rhs of |nabla~ eps|^2 - |D~ eps|^2 + 4 pi (mu |eps|^2 + sJ J.xi) = psi^-6 d_a(psi^4 V^a) on the data at the point
    pt, exactly. The left side uses weitzenbock_conf's Geo (spin connection) and constraints (mu, J), not this script's."""
    rep = Jet(pt)
    D = data(psi2="general", P=True, K3=True, T=True, eps1=True, eps2=True, eta=True)
    psi, Kf, eps = fields(rep, D, vals)
    V = Vvec(psi, Kf, eps, sK)
    div = fsum([(wfac(psi) * V[a]).d(a) for a in range(3)], Fld(rep, sp.Integer(0), None))
    z0 = {v: 0 for v in Jet.DV}
    pre = rep.ipow(psi, 6) if MUT != "weight" else rep.const(sp.Integer(1))
    rhs = sp.expand((pre * div).v.subs(z0))
    # left side
    if "generic" not in _CACHE:
        psi_f, Kg = WC.fn("psi"), WC.sym_K()
        geo = WC.Geo(psi_f)
        _CACHE["generic"] = (geo,) + tuple(WC.constraints(psi_f, Kg, geo))
    geo, mu_g, J_g = _CACHE["generic"]
    jets = {"psi": psi.v}
    for i in range(3):
        for j in range(i, 3):
            jets[f"K{i + 1}{j + 1}"] = Kf.v[i, j]
    exprs = [mu_g] + list(J_g) + [M[i, j] for M in geo.Om for i in range(4) for j in range(4)]
    av = atom_values(exprs, jets)
    mu = sp.expand(mu_g.xreplace(av))
    J = [sp.expand(e.xreplace(av)) for e in J_g]
    Om = [M.xreplace(av).applyfunc(sp.expand) for M in geo.Om]
    g0, gj = gammas()
    e = eps.v.subs(z0)
    ps = psi.v.subs(z0)
    Kn = Kf.v.subs(z0)
    nabt = []
    for i in range(3):
        de = eps.v.applyfunc(lambda f: sp.diff(f, Jet.DV[i]).subs(z0))
        nab = de / ps**2 + Om[i] * e
        nabt.append(sp.expand(nab + sK * sum((Kn[i, j] * gj[j] * g0 for j in range(3)), sp.zeros(4)) * e / 2))
    Dt = sum((gj[i] * nabt[i] for i in range(3)), sp.zeros(4, 1))
    xi = [WC.inner(e, g0 * gj[j] * e) for j in range(3)]
    lhs = sum(WC.inner(nabt[i], nabt[i]) for i in range(3)) - WC.inner(Dt, Dt) \
        + 4 * sp.pi * (mu * WC.inner(e, e) + sJ * sum(J[j] * xi[j] for j in range(3)))
    lhs = sp.expand(lhs)
    return sp.expand(lhs - rhs), lhs, rhs


PTS = [(1, 2, 2), (-2, 3, 6), (4, -4, 7)]


def section_identity():
    print("\n== 5. V against the pointwise identity, on the same data, at exact points ==")
    names = ("identity_points", "flip_consistent")
    if not any(todo(n) for n in names):
        return
    vals = values()
    sJ = 1 if MUT == "Jsign" else -1
    if todo("identity_points"):
        t0, ok, msgs = time.time(), True, []
        for pt in PTS:
            res, lhs, rhs = identity_residual(pt, vals, -1 if MUT == "sign" else 1, sJ)
            ok = ok and res == 0
            msgs.append(f"{pt}: lhs {sp.N(lhs, 8)}, rhs {sp.N(rhs, 8)}, lhs - rhs {sp.N(res, 3)}")
        check("identity_points", ok, "|nabla~ eps|^2 - |D~ eps|^2 + 4 pi (mu |eps|^2 - J.xi) = psi^-6 d_a(psi^4 V^a), with V as used in "
              "the flux, on data with every symbol group (psi_2, psi_3, T, K3, eps_0, eta, eps_1, eps_2) at exact rational values: "
              + "; ".join(msgs) + f" ({time.time() - t0:.0f} s)")
    if todo("flip_consistent"):
        # The other sign of nabla~, nabla_i - (1/2) K_ij g^j g^0, with +J.xi: the identity holds with V's K term +(1/2) P_ij xi^j,
        # and the boundary term of that V is E |eps_0|^2 + P.xi_0.
        res, lhs, rhs = identity_residual(PTS[0], vals, -1, 1)
        R = series_run("BY", P=True)
        psi, Kf, eps = fields(SER, R["D"])
        fl, _ = flux_series(psi, Kf, eps, sK=-1)
        Pxi = sum(sym(f"P{c}") * R["xi"][j] for j, c in enumerate("xyz"))
        ok = res == 0 and agree(fl, {0: sym("m") * R["e2"] + Pxi}, (-2, -1, 0))
        check("flip_consistent", ok, "with nabla~ = nabla - (1/2) K g g^0 the identity holds with +J.xi and V_K = +(1/2) P_ij xi^j "
              "(residual " + str(res) + f" at {PTS[0]}), and the boundary term is m |eps_0|^2 + P.xi_0: the sign of P.xi_0 follows nabla~")


def main():
    print("weitzenbock_boundary.py, issue #100. sympy %s, python %s. mutation: %s" % (sp.__version__, sys.version.split()[0], MUT or "none"))
    section_conventions()
    section_machinery()
    section_planted()
    section_aligned()
    section_general()
    section_identity()
    fails = [n for n, ok in RESULTS if not ok]
    print("\nchecks: %d passed, %d failed %s. %.0f s" % (len(RESULTS) - len(fails), len(fails), fails, time.time() - T0))
    if MUT:
        target = CONTROLS[MUT][1]
        print("control %s: targets %s %s" % (MUT, target, "FIRED" if all(n in fails for n in target) else "SILENT"))
        sys.exit(1)
    status = 1 if fails else 0
    if "no-controls" not in OPTS and ONLY is None and not SKIP - {""}:
        print("\n== Controls: each mutation must make its target checks fail ==")
        for name, (what, target) in CONTROLS.items():
            t0 = time.time()
            p = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True)
            failed = [ln.split("]")[1].split(":")[0].strip() for ln in p.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = p.returncode == 1 and all(n in failed for n in target)
            status |= 0 if fired else 1
            print("[%s] --mutate=%s (%s): exit %d, failed %s (%.0f s)" % ("FIRED" if fired else "SILENT", name, what, p.returncode, failed,
                                                                         time.time() - t0), flush=True)
        print("controls: %s" % ("all fired" if status == 0 else "SOME SILENT OR CHECK FAILED"))
    sys.exit(status)


if MUT and MUT not in CONTROLS:
    sys.exit(f"unknown mutation {MUT}")

if __name__ == "__main__":
    main()

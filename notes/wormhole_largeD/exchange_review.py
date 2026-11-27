"""Issue #209: the exchange identities E1, E2, E3 and the Clausen integrals behind Omega_2 (#26) and Omega_3 (#208), re-derived.

The derivations are in exchange_review.md, written before omega2.md and omega3.md were opened. This script checks every step.
Notation: R = r^s = sec^2 th, l = ln R = -2 ell with ell = ln cos th, u = tan th, S = sin^2 th, C = Cl2(pi - 2 th).
The order equations (order4.py's docstring, README.md "Method"):
  dY_k/dl = F_{k-1} = sum_{j<k} L_{k-1-j} l^j/j!,   (R - 1) L_k' - L_k/2 = h_k,
  h_k = (1 - 1/R) L_{k-1} + sum_{j<k} Y_j c_{k-j} + Y_k/2,   L_k = -2 Omega_{k-1} + u C_k,   G_k = C_k(pi/2) = 0.
The c_m, the panel hierarchy (method B) and its Omega chain are imported from order5.py (cc, hk, PanelK, run_B); nothing else of
order4.py and order5.py is used.

Parts (run in this order; no argument runs sym, exact, clausen, exch):
  sym      sympy: the algebra of E1, E2, E3, the boundary limits, the moments, the Clausen derivative and the sine-series integrals
  exact    E1 at k = 2 and E2 at m = 1, 2 from closed forms, by mpmath quadrature, with no panel code (planted truth)
  clausen  the two Clausen integrals and the moments of ln cos to 40 digits, and the controls
  exch     E1 (k = 3, 4, 5), E2 (m = 2, 3, 4), E3 (k = 3, 4, 5), their boundary terms, on the panel hierarchy of the E-grid
  grid     Omega_2 and Omega_3 on the third method-B grid, against #26's A70 and #208's B4 (about 100 s on rpi5)

    python3 exchange_review.py [sym] [exact] [clausen] [exch] [grid] [--mutate=NAME] [--egrid=N,w,V,dps] [--grid=N,w,V,dps]

The exit status is 1 when a check fails. Status words (AGENTS.md): sympy residual 0 is "checked symbolically", an agreement of
integrals with each other or with a constant "numerical evidence". A digit count is -log10 of a relative difference between two
evaluations; it is not an error bar.
"""
import math
import os
import re
import sys
import time

import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
_argv, sys.argv = sys.argv, sys.argv[:1]          # order4.py parses sys.argv when it is imported
import order5 as o5  # noqa: E402
sys.argv = _argv

MUTATIONS = ("sym_S1", "sym_S2", "sym_c1", "sym_P", "sym_S6", "sym_S7", "sym_S8", "sym_S9", "sym_S10", "sym_S11",
             "exact_K", "exact_L1", "exact_h2", "e1_sign", "e1_factor", "e2_weight", "bnd_G", "bnd_Y0", "e3_Om0", "e3_factor",
             "cl_moment", "cl_arg", "cl_arg2", "grid_coarse", "grid_short", "grid_wrongh", "grid_wrongOm1")
MUTATE = None
PARTS = []
EGRID = (36, (1, 2), 16, 80)          # N, w as (num, den), V, dps: the grid of the E checks (declared in the raw README)
GRID3 = (64, (1, 3), 24, 130)         # the third method-B grid (declared in the raw README)
ORDER4_OUT = os.path.join(HERE, "order4_output.txt")
ORDER5_OUT = os.path.join(HERE, "order5_output.txt")
OUTFILE = None
CTLDIR = None
CTLPREFIX = "ctl"

for a in sys.argv[1:]:
    if a.startswith("--mutate="):
        MUTATE = a.split("=", 1)[1]
        assert MUTATE in MUTATIONS, f"unknown mutation {MUTATE}; known: {MUTATIONS}"
    elif a.startswith("--egrid="):
        n, w, v, d = a.split("=", 1)[1].split(",")
        EGRID = (int(n), tuple(int(t) for t in w.split("/")) if "/" in w else (int(w), 1), int(v), int(d))
    elif a.startswith("--out="):
        OUTFILE = a.split("=", 1)[1]
    elif a.startswith("--ctldir="):
        CTLDIR = a.split("=", 1)[1]
    elif a.startswith("--ctlprefix="):
        CTLPREFIX = a.split("=", 1)[1]
    elif a.startswith("--grid="):
        n, w, v, d = a.split("=", 1)[1].split(",")
        GRID3 = (int(n), tuple(int(t) for t in w.split("/")) if "/" in w else (int(w), 1), int(v), int(d))
    else:
        PARTS.append(a)
PARTS = PARTS or ["sym", "exact", "clausen", "exch"]

RESULTS = []


def report(cid, ok, text):
    RESULTS.append((cid, ok))
    print(f"{'PASS' if ok else 'FAIL'}    {cid}: {text}", flush=True)


def digits(a, b):
    """-log10 of the relative difference of two numbers (to the larger one), capped at the working precision in digits"""
    a, b = mp.mpf(a), mp.mpf(b)
    d = abs(a - b)
    return mp.mpf(mp.mp.dps) if d == 0 else min(mp.mpf(mp.mp.dps), -mp.log10(d / max(abs(a), abs(b))))


def dg(a, b):
    return f"{float(digits(a, b)):.1f}"


def nd(x, n=30):
    return mp.nstr(x, n)


def cl2(x):
    """Cl_2(x) = Im Li_2(e^{ix}), real; mpmath's clsin returns an mpc with a rounding-size imaginary part for some arguments"""
    return mp.re(mp.clsin(2, x))


# ---------------------------------------------------------------------------------------------------------------
# sym: sympy
# ---------------------------------------------------------------------------------------------------------------
def check_sym():
    t0 = time.time()
    th = sp.Symbol("theta", positive=True)
    Om0 = sp.Symbol("Omega0")
    n, m = sp.symbols("n m", positive=True, integer=True)
    ell = sp.log(sp.cos(th))
    S, cth, tth = sp.sin(th) ** 2, sp.cos(th), sp.tan(th)
    R = 1 / cth ** 2
    lq = -2 * ell
    zero = lambda e: sp.simplify(e) == 0  # noqa: E731

    # S1: the substitution R = sec^2 th turns (R - 1)^(-3/2) dR into 2 dth / sin^2 th (u = sqrt(R - 1) = tan th)
    ok = zero(R - 1 - tth ** 2) and zero(sp.diff(R, th) / tth ** 3 - (2 if MUTATE != "sym_S1" else 1) / S)
    report("S1", ok, "R - 1 = tan^2 and (R - 1)^(-3/2) dR/dth = 2/sin^2: residuals 0 (sympy)")

    # S2: 1 - 1/R = sin^2, l = ln R = -2 ln cos, dl/dth = 2 tan
    ok = zero(1 - 1 / R - S) and zero(sp.expand_log(sp.log(R) - lq, force=True)) and zero(sp.diff(lq, th) - (2 if MUTATE != "sym_S2" else 1) * tth)
    report("S2", ok, "1 - 1/R = sin^2, ln R = -2 ln cos, d(ln R)/dth = 2 tan: residuals 0 (sympy)")

    # S3: c_1 of the table, divided by sin^2: (Omega0 - 3 ell)/sin^2 - 1
    c1 = o5.cc(1, [Om0], lq, R)
    if MUTATE == "sym_c1":
        c1 = c1 - lq / 2                                      # 3l/2 taken as l
    ok3 = zero(c1 / S - ((Om0 - 3 * ell) / S - 1))
    report("S3", ok3, "c_1/sin^2 = (Omega0 - 3 ell)/sin^2 - 1 with c_1 = Omega0 + 3l/2 - 1 + 1/R, l = -2 ell: residual 0 (sympy)")

    # S4: P = cot (ell - Omega0): P' = -[(ell - Omega0)/sin^2 + 1] and P dl/dth = 2 (ell - Omega0)
    P = (cth / sp.sin(th)) * (ell - Om0)
    sgn = 1 if MUTATE != "sym_P" else -1
    ok4 = zero(sp.diff(P, th) + sgn * ((ell - Om0) / S + 1)) and zero(P * sp.diff(lq, th) - 2 * (ell - Om0))
    report("S4", ok4, "P' = -[(ell - Omega0)/sin^2 + 1] and P dl/dth = 2 (ell - Omega0), P = cot (ell - Omega0): residuals 0 (sympy)")

    # S5: the Y_{k-1} terms of G_k add up to -2 [(ell - Omega0)/sin^2 + 1] = 2 P'
    ok5 = zero(4 * ell / S + 2 * c1 / S + 2 * ((ell - Om0) / S + 1))
    report("S5", ok5, "4 ell/sin^2 (from E2) + 2 c_1/sin^2 (from h_k) = -2[(ell - Omega0)/sin^2 + 1] = 2 P': residual 0 (sympy)")

    # S6: ell/sin^2 = d/dth (-cot ell - th); the boundary values of -cot ell - th
    A = -(cth / sp.sin(th)) * ell - (th if MUTATE != "sym_S6" else -th)
    lim0, lim1 = sp.limit(A, th, 0, "+"), sp.limit(A, th, sp.pi / 2, "-")
    ok = zero(sp.diff(A, th) - ell / S) and lim0 == 0 and sp.simplify(lim1 + sp.pi / 2) == 0
    report("S6", ok, f"d/dth(-cot ell - th) = ell/sin^2; boundary values {lim0} at 0 and {lim1} at pi/2, so int ell/sin^2 = -pi/2 (sympy)")

    # S7: boundary limits used by the by-parts steps (x = pi/2 - th at the outer end)
    x = sp.Symbol("x", positive=True)
    lims = {
        "cot ell at 0": sp.limit((cth / sp.sin(th)) * ell, th, 0, "+"),
        "cot ell at pi/2": sp.limit((cth / sp.sin(th)) * ell if MUTATE != "sym_S7" else ell, th, sp.pi / 2, "-"),
        "ell th^-1 ...: ell * th at 0 (C = O(th))": sp.limit(ell * th, th, 0, "+"),
        "x ln^6 x at 0 (Y_4 P, Y_4 = O(l^5), P ~ x (ln x - Omega0))": sp.limit(x * sp.log(x) ** 6, x, 0, "+"),
        "x ln^2 x at 0 (ell C_Cl, C ~ x(1 - ln x))": sp.limit(x * sp.log(x) ** 2, x, 0, "+"),
    }
    ok = all(v == 0 for v in lims.values())
    series_throat = sp.series(ell, th, 0, 6).removeO()
    report("S7", ok, f"boundary limits {', '.join(f'{k} = {v}' for k, v in lims.items())}; ell = {series_throat} + O(th^6) (sympy)")

    # S8: Beta-function moments: B(a) = int_0^{pi/2} cos^a = sqrt(pi) Gamma((a+1)/2) / (2 Gamma(a/2 + 1))
    a = sp.Symbol("a")
    Bf = sp.sqrt(sp.pi) * sp.gamma((a + 1) / 2) / (2 * sp.gamma(a / 2 + 1))
    M0, M1, M2 = (sp.simplify(sp.diff(Bf, a, j).subs(a, 0)) for j in range(3))
    M1c, M2c = -sp.pi / 2 * sp.log(2), sp.pi / 2 * (sp.log(2) ** 2 + sp.pi ** 2 / 12)
    if MUTATE == "sym_S8":
        M1c = -M1c
    ok = (M0 == sp.pi / 2 and sp.simplify(sp.expand_func(M1) - M1c) == 0 and sp.simplify(sp.expand_func(M2) - M2c) == 0)
    report("S8", ok, f"moments from B(a): int ell = {M1}, int ell^2 = {sp.simplify(sp.expand_func(M2))} (sympy)")

    # S9: C' = d/dth Cl2(pi - 2 th) = 2 ln(2 cos th), from Cl2'(x) = -ln(2 sin(x/2)); and the Clausen value from the moments
    x_ = sp.Symbol("x_", positive=True)
    Cl2p = lambda z: -sp.log(2 * sp.sin(z / 2))  # noqa: E731
    dC = sp.diff(sp.pi - 2 * th, th) * Cl2p(sp.pi - 2 * th)
    ok = zero(sp.expand_log(sp.simplify(dC - 2 * sp.log(2 * (cth if MUTATE != "sym_S9" else sp.sin(th)))), force=True))
    val = sp.simplify(sp.expand_func(2 * (M2 + sp.log(2) * M1)))
    ok = ok and sp.simplify(val - sp.pi ** 3 / 12) == 0
    report("S9", ok, f"C' = 2 ln(2 cos th); 2 (int ell^2 + ln2 int ell) = {val} = pi^3/12 (sympy)")

    # S10: orthogonality of sin(2n th) on [0, pi/2] and int tan sin(2n th) = (-1)^(n+1) pi/2, by product-to-sum identities
    k = sp.Symbol("k", integer=True, nonzero=True)
    int_cos = sp.integrate(sp.cos(2 * k * th), (th, 0, sp.pi / 2))            # sin(pi k)/(2k) = 0 for a nonzero integer k
    prod = sp.expand_trig(sp.sin(2 * n * th) * sp.sin(2 * m * th)) - (sp.expand_trig(sp.cos(2 * (n - m) * th)) - sp.expand_trig(sp.cos(2 * (n + m) * th))) / 2
    ok = int_cos == 0 and zero(prod)                                          # sin A sin B = [cos(A - B) - cos(A + B)]/2; n +- m integers
    diag = sp.integrate(sp.sin(2 * n * th) ** 2, (th, 0, sp.pi / 2))
    ok = ok and sp.simplify(diag.subs(n, 3) - sp.pi / 4) == 0
    # tan [sin 2(n+1)th - sin 2(n-1)th] = 2 cos 2n th - cos 2(n+1)th - cos 2(n-1)th, an identity in n
    lhs = sp.tan(th) * sp.expand_trig(sp.sin(2 * n * th + 2 * th) - sp.sin(2 * n * th - 2 * th))
    rhs = 2 * sp.cos(2 * n * th) - sp.expand_trig(sp.cos(2 * n * th + 2 * th) + sp.cos(2 * n * th - 2 * th))
    ok = ok and zero(sp.simplify(lhs) - sp.simplify(rhs))
    I1 = sp.integrate(sp.expand(2 * sp.sin(th) ** 2), (th, 0, sp.pi / 2))     # tan sin 2th = 2 sin^2
    I2 = sp.integrate(sp.expand(4 * sp.sin(th) ** 2 * sp.cos(2 * th)), (th, 0, sp.pi / 2))   # tan sin 4th = 4 sin^2 cos 2th
    ok = ok and sp.simplify(I1 - sp.pi / 2) == 0 and sp.simplify(I2 + (1 if MUTATE != "sym_S10" else -1) * sp.pi / 2) == 0
    report("S10", ok, "int sin(2n)sin(2m) = (pi/4) delta_nm; tan [sin 2(n+1)th - sin 2(n-1)th] = 2 cos 2n th - cos 2(n+1)th - cos 2(n-1)th "
                      "(residual 0), whose integral gives I_(n+1) = I_(n-1) for n >= 2 with I_1 = pi/2, I_2 = -pi/2, "
                      "I_n = int tan sin(2n th) = (-1)^(n+1) pi/2 (sympy)")

    # S11: the series sums of the second route: sum (pi/2)/n^2 = pi^3/12, sum pi/(4 n^4) = pi^5/360 (Parseval)
    nn = sp.Symbol("nn", positive=True, integer=True)
    s2, s4 = sp.summation(sp.pi / 2 / nn ** 2, (nn, 1, sp.oo)), sp.summation(sp.pi / 4 / nn ** 4, (nn, 1, sp.oo))
    ok = sp.simplify(s2 - (sp.pi ** 3 / 12 if MUTATE != "sym_S11" else sp.pi ** 3 / 6)) == 0 and sp.simplify(s4 - sp.pi ** 5 / 360) == 0
    report("S11", ok, f"sum_n (pi/2)/n^2 = {s2} and sum_n (pi/4)/n^4 = {s4} (sympy)")
    print(f"  sym: {time.time() - t0:.0f} s", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# exact: E1 at k = 2 and E2 at m = 1, 2 from the closed forms of README.md "Method", no panel code
# ---------------------------------------------------------------------------------------------------------------
def closed_fields(th, Om0, Om1):
    """Y_0..Y_2, L_0..L_2 at th from README.md's closed forms (L_2 with K = 2 th ln 2 - Cl2(pi - 2 th)), and l, R, 1 - 1/R"""
    c, s = mp.cos(th), mp.sin(th)
    ell = mp.log(c)
    l, R, S, u = -2 * ell, 1 / c ** 2, s ** 2, mp.tan(th)
    K = (2 if MUTATE != "exact_K" else 1) * th * mp.ln2 - cl2(mp.pi - 2 * th)
    Y = [mp.mpf(1), -l, -2 * Om0 * l - mp.mpf(3) / 2 * l ** 2]
    L1 = -2 * Om0 - 2 * l if MUTATE != "exact_L1" else -2 * Om0 - l
    L = [mp.mpf(-1), L1, -2 * Om1 + 2 * Om0 * l + 4 * l ** 2 - 8 * Om0 * u * th - 16 * u * K]
    return Y, L, ell, l, R, S


def qthird(f, dps):
    """int_0^{pi/2} f dth by tanh-sinh at dps digits, as int_0^{pi/4} f(th) dth + int_0^{pi/4} f(pi/2 - x) dx; f is evaluated at 3 dps,
    so that (h - Omega)/sin^2 near th = 0 does not cancel and a node x near 0 does not land beyond pi/2"""
    def g_in(th):
        with mp.workdps(3 * dps):
            v = f(+th)
        return +v

    def g_out(x):
        with mp.workdps(3 * dps):
            v = f(mp.pi / 2 - x)
        return +v
    with mp.workdps(dps):
        a, ea = mp.quad(g_in, [0, mp.pi / 4], error=True)
        b, eb = mp.quad(g_out, [0, mp.pi / 4], error=True)
        return a + b, max(ea, eb)


def check_exact():
    t0 = time.time()
    dps = 50
    mp.mp.dps = dps
    # K(u) against its definition, the planted pieces of L_2
    worst = mp.inf
    for th in (mp.mpf("0.3"), mp.mpf("0.9"), mp.mpf("1.4")):
        u = mp.tan(th)
        Kq = mp.quad(lambda t: mp.log(1 + t * t) / (1 + t * t), [0, u])
        Kc = (2 if MUTATE != "exact_K" else 1) * th * mp.ln2 - cl2(mp.pi - 2 * th)
        worst = min(worst, digits(Kq, Kc)) if worst != mp.inf else digits(Kq, Kc)
    report("P0b", worst > 40, f"K(u) = int_0^u ln(1 + t^2)/(1 + t^2) dt against 2 th ln 2 - Cl2(pi - 2 th), th = 0.3, 0.9, 1.4: worst "
                             f"{float(worst):.1f} digits (bar 40)")
    # the closed forms against the table: h_1 and h_2 of README
    Om0, Om1 = mp.mpf(3) / 10, mp.mpf(1) / 7
    worst = mp.inf
    for th in (mp.mpf("0.2"), mp.mpf("0.8"), mp.mpf("1.3")):
        Y, L, ell, l, R, S = closed_fields(th, Om0, Om1)
        h1 = o5.hk(1, [Om0, Om1], l, R, Y, L)
        h2 = o5.hk(2, [Om0, Om1], l, R, Y, L)
        h1r = Om0 + l - 2 * S
        h2r = Om1 - 2 * Om0 * S - Om0 * l - 2 * l ** 2
        if MUTATE == "exact_h2":
            h2r = h2r + l / 10
        worst = min(worst, digits(h1, h1r), digits(h2, h2r)) if worst != mp.inf else min(digits(h1, h1r), digits(h2, h2r))
    report("P0c", worst > 40, f"h_1, h_2 from the table (cc, hk) against README's h_1 = Omega0 + l - 2(1 - 1/R), h_2 = Omega1 - 2 Omega0 (1 - 1/R) "
                             f"- Omega0 l - 2 l^2, th = 0.2, 0.8, 1.3: worst {float(worst):.1f} digits (bar 40)")

    # E1 at k = 2, Omega0 not at its root: G_2 = -4 pi Omega0 - 2 J, J = 8 pi ln 2
    Om0, Om1 = mp.mpf(3) / 10, mp.mpf(0)
    Om = [Om0, Om1]

    def e1_parts(th):
        Y, L, ell, l, R, S = closed_fields(th, Om0, Om1)
        S2 = Y[0] * o5.cc(2, Om, l, R) + Y[1] * o5.cc(1, Om, l, R)
        sg = -1 if MUTATE == "e1_sign" else 1
        fac = 2 if MUTATE == "e1_factor" else 4
        return fac * L[1] + sg * 2 * L[0] * l + 2 * (S2 - Om1) / S
    def g_def(th):
        Y, L, ell, l, R, S = closed_fields(th, Om0, Om1)
        return 2 * (o5.hk(2, Om, l, R, Y, L) - Om1) / S
    (rhs, e1), (lhs, e2) = qthird(e1_parts, dps), qthird(g_def, dps)
    target = -4 * mp.pi * Om0 - 16 * mp.pi * mp.ln2
    ok = digits(lhs, target) > 40 and digits(rhs, target) > 40
    report("P1", ok, f"E1 at k = 2, Omega0 = 3/10: G_2 from the definition {nd(lhs, 35)}, E1's right side {nd(rhs, 35)}, against "
                     f"-4 pi Omega0 - 16 pi ln 2 = {nd(target, 35)}: {dg(target, lhs)} and {dg(target, rhs)} digits (bar 40; quadrature "
                     f"error estimates {mp.nstr(e1, 2)}, {mp.nstr(e2, 2)})")

    # E2 at m = 1 and m = 2, Omega0 = -4 ln 2 (G_2 = 0 needs it for m = 2; G_1 = 0 holds for any Omega0)
    Om0, Om1 = -4 * mp.ln2, mp.mpf(1) / 3
    Om = [Om0, Om1]

    def parts(m):
        def lhs_f(th):
            Y, L, ell, l, R, S = closed_fields(th, Om0, Om1)
            return L[m]

        def rhs_f(th):
            Y, L, ell, l, R, S = closed_fields(th, Om0, Om1)
            w = ell if MUTATE != "e2_weight" else ell + mp.ln2
            return 2 * o5.hk(m, Om, l, R, Y, L) * w / S
        return lhs_f, rhs_f
    for m in (1, 2):
        lf, rf = parts(m)
        (a, ea), (b, eb) = qthird(lf, dps), qthird(rf, dps)
        extra = ""
        if m == 1:
            tgt = -mp.pi * Om0 - 2 * mp.pi * mp.ln2
            extra = f"; closed form -pi Omega0 - 2 pi ln 2 = {nd(tgt, 30)}: {dg(tgt, a)} digits"
            ok = digits(a, b) > 40 and digits(tgt, a) > 40
        else:
            ok = digits(a, b) > 40
        report(f"P2_{m}", ok, f"E2 at m = {m}: int L_{m} dth = {nd(a, 35)} against 2 int h_{m} ln cos/sin^2 dth = {nd(b, 35)}: "
                              f"{dg(a, b)} digits (bar 40){extra}")
    print(f"  exact: {time.time() - t0:.0f} s", flush=True)

# ---------------------------------------------------------------------------------------------------------------
# clausen: the two Clausen integrals and the moments of ln cos, to 40 digits
# ---------------------------------------------------------------------------------------------------------------
def Cfun(th):
    arg = mp.pi - 2 * th
    if MUTATE == "cl_arg":
        arg = mp.pi - th
    elif MUTATE == "cl_arg2":
        arg = 2 * th
    return cl2(arg)


def check_clausen():
    t0 = time.time()
    dps = 60
    mp.mp.dps = dps
    pi = mp.pi
    sh = mp.ln2 if MUTATE == "cl_moment" else 0
    m1, e1 = qthird(lambda th: mp.log(mp.cos(th)) + sh, dps)
    m2, e2 = qthird(lambda th: (mp.log(mp.cos(th)) + sh) ** 2, dps)
    t1, t2 = -pi / 2 * mp.ln2, pi / 2 * (mp.ln2 ** 2 + pi ** 2 / 12)
    report("CM", min(digits(m1, t1), digits(m2, t2)) > 40,
           f"int_0^(pi/2) ln cos = {nd(m1, 30)}, int ln^2 cos = {nd(m2, 30)} against -(pi/2) ln 2 and (pi/2)(ln^2 2 + pi^2/12): "
           f"{dg(t1, m1)} and {dg(t2, m2)} digits (bar 40)")
    # CL1: int tan C = pi^3/12
    a, ea = qthird(lambda th: mp.tan(th) * Cfun(th), dps)
    report("CL1", digits(a, pi ** 3 / 12) > 40,
           f"int_0^(pi/2) tan th Cl2(pi - 2 th) dth = {nd(a, 40)} against pi^3/12 = {nd(pi ** 3 / 12, 40)}: {dg(pi ** 3 / 12, a)} digits (bar 40; "
           f"quadrature error estimate {mp.nstr(ea, 2)})")
    # CL1a: the by-parts step, int tan C = 2 int ln cos ln(2 cos)
    b, eb = qthird(lambda th: 2 * mp.log(mp.cos(th)) * mp.log(2 * mp.cos(th)), dps)
    report("CL1a", digits(a, b) > 40 and digits(b, pi ** 3 / 12) > 40,
           f"the by-parts step: int tan C dth = {nd(a, 30)} against 2 int ln cos ln(2 cos) dth = {nd(b, 30)}: {dg(a, b)} digits (bar 40)")
    # CL2: int C^2 = pi^5/360
    c, ec = qthird(lambda th: Cfun(th) ** 2, dps)
    report("CL2", digits(c, pi ** 5 / 360) > 40,
           f"int_0^(pi/2) Cl2(pi - 2 th)^2 dth = {nd(c, 40)} against pi^5/360 = {nd(pi ** 5 / 360, 40)}: {dg(pi ** 5 / 360, c)} digits (bar 40; "
           f"quadrature error estimate {mp.nstr(ec, 2)})")
    print(f"  clausen: {time.time() - t0:.0f} s", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# exch: E1, E2, E3 and their boundary terms on the panel hierarchy (method B's PanelK, imported)
# ---------------------------------------------------------------------------------------------------------------
class Fields:
    """Y_0..Y_k, L_0..L_k, G_k of the hierarchy for one Omega list, with the nodes' l, R, ell, sin^2 and the th-quadrature"""

    def __init__(self, B, Om, k):
        self.B, self.Om, self.k = B, list(Om), k
        self.N, self.P = B.N, B.P
        self.Y, self.L, self.G = B.hierarchy(self.Om, k)
        self.l, self.R, self.u = B.l, B.R, B.u
        self.S = self.grid(lambda p, i: B.Rm1[p][i] / B.R[p][i])                 # sin^2 th = 1 - 1/R
        self.ell = self.grid(lambda p, i: -B.l[p][i] / 2)                        # ln cos th = -l/2
        self.W = self.grid(lambda p, i: B.half * B.wend[i] * B.v[p][i] / B.u[p][i])   # th-quadrature weights: dth = (v/u) dv
        if MUTATE == "bnd_Y0":
            self.Y = [self.Y[0]] + [self.grid(lambda p, i, y=y: y[p][i] + mp.mpf(1) / 1000) for y in self.Y[1:]]

    def grid(self, fn):
        return [[fn(p, i) for i in range(self.N)] for p in range(self.P)]

    def node(self, p, i):
        return ([y[p][i] for y in self.Y], [lv[p][i] for lv in self.L])

    def integral(self, f):
        """int_0^{pi/2} f dth, f(p, i)"""
        return sum(self.W[p][i] * f(p, i) for p in range(self.P) for i in range(self.N))

    def h(self, m, p, i):
        Y, L = self.node(p, i)
        return o5.hk(m, self.Om, self.l[p][i], self.R[p][i], Y, L if m else None)

    def F(self, j, p, i):
        """F_j = sum_{q <= j} L_{j-q} l^q / q!"""
        return sum(self.L[j - q][p][i] * self.l[p][i] ** q / mp.factorial(q) for q in range(j + 1))

    def Sk(self, k, p, i, upto=None):
        """sum_{j < upto} Y_j c_{k-j}"""
        upto = k if upto is None else upto
        return sum(self.Y[j][p][i] * o5.cc(k - j, self.Om, self.l[p][i], self.R[p][i]) for j in range(upto))


def e1_residual(fl, k):
    """G_k from the hierarchy against E1's right side"""
    om = fl.Om
    A = 4 * fl.integral(lambda p, i: fl.L[k - 1][p][i])
    if MUTATE == "e1_factor":
        A = A / 2
    Bs = 2 * fl.integral(lambda p, i: fl.F(k - 1, p, i) - fl.L[k - 1][p][i])
    if MUTATE == "e1_sign":
        Bs = -Bs
    Cs = 2 * fl.integral(lambda p, i: (fl.Sk(k, p, i) - om[k - 1]) / fl.S[p][i])
    return fl.G[k], A + Bs + Cs, (A, Bs, Cs)


def e2_sides(fl, m):
    lhs = fl.integral(lambda p, i: fl.L[m][p][i])
    shift = mp.ln2 if MUTATE == "e2_weight" else 0
    rhs = 2 * fl.integral(lambda p, i: fl.h(m, p, i) * (fl.ell[p][i] + shift) / fl.S[p][i])
    return lhs, rhs


def e3_sides(fl, m):
    """E3 at k = m + 1: int Y_m [(ell - Omega0)/sin^2 + 1] against 2 int (ell - Omega0) F_{m-1}; both need Omega0 = fl.Om[0]"""
    O0 = fl.Om[0] if MUTATE != "e3_Om0" else 0
    lhs = fl.integral(lambda p, i: fl.Y[m][p][i] * ((fl.ell[p][i] - O0) / fl.S[p][i] + 1))
    rhs = 2 * fl.integral(lambda p, i: (fl.ell[p][i] - O0) * fl.F(m - 1, p, i))
    return lhs, rhs


def g_reduced(fl, k):
    """G_k with no Y_{k-1}, no Y_k: 8 int (h_m - Y_m/2) ell/S + 2 int (F_m - L_m) + 2 int (sum_{j<m} Y_j c_{k-j} - Omega_{k-1})/S
    - 4 int (ell - Omega0) F_{m-1}, m = k - 1"""
    m, om = k - 1, fl.Om
    O0 = om[0] if MUTATE != "e3_Om0" else 0
    fac = 4 if MUTATE != "e3_factor" else 2

    def hm_noY(p, i):
        return fl.h(m, p, i) - fl.Y[m][p][i] / 2
    t1 = 8 * fl.integral(lambda p, i: hm_noY(p, i) * fl.ell[p][i] / fl.S[p][i])
    t2 = 2 * fl.integral(lambda p, i: fl.F(m, p, i) - fl.L[m][p][i])
    t3 = 2 * fl.integral(lambda p, i: (fl.Sk(k, p, i, upto=m) - om[k - 1]) / fl.S[p][i])
    t4 = -fac * fl.integral(lambda p, i: (fl.ell[p][i] - O0) * fl.F(m - 1, p, i))
    return t1 + t2 + t3 + t4, (t1, t2, t3, t4)


def check_exch():
    t0 = time.time()
    N, w, V, dps = EGRID
    mp.mp.dps = dps
    with mp.workdps(dps):
        wv = mp.mpf(w[0]) / w[1]
        B = o5.PanelK(N, wv, V)
        Om, slopes, secs = o5.omega_chain(B, 5)
    cl = o5.closed_omegas()
    print(f"  E-grid: N = {N}, w = {w[0]}/{w[1]}, V = {V}, dps = {dps}, panels = {B.P}; chain {sum(secs.values()):.0f} s; "
          f"Omega_0..2 against closed forms {', '.join(dg(cl[j], Om[j]) for j in range(3))} digits; Omega_3 = {nd(Om[3], 40)}", flush=True)
    bar = 30
    base = list(Om)
    if MUTATE == "bnd_G":
        base[0] = base[0] + mp.mpf(1) / 100                 # G_2 != 0: L_m grows like u
    # planted: the hierarchy against README's closed forms through order 2
    with mp.workdps(dps):
        fl = Fields(B, base, 3)
        worst = {"Y1": 0, "L1": 0, "Y2": 0, "L2": 0}
        lmax = V * V // 4           # L_k = -2 Omega - u C_right loses the tail cut at V, and amplifies rounding by u = e^(l/2), at large l
        used = 0
        for p in range(B.P):
            for i in range(B.N):
                th, u, l = mp.atan(B.u[p][i]), B.u[p][i], B.l[p][i]
                if l > lmax:
                    continue
                used += 1
                K = 2 * th * mp.ln2 - cl2(mp.pi - 2 * th)
                ref = {"Y1": -l, "L1": -2 * base[0] - 2 * l, "Y2": -2 * base[0] * l - mp.mpf(3) / 2 * l ** 2,
                       "L2": -2 * base[1] + 2 * base[0] * l + 4 * l ** 2 - 8 * base[0] * u * th - 16 * u * K}
                got = {"Y1": fl.Y[1][p][i], "L1": fl.L[1][p][i], "Y2": fl.Y[2][p][i], "L2": fl.L[2][p][i]}
                for key in worst:
                    worst[key] = max(worst[key], abs(got[key] - ref[key]) / (1 + abs(ref[key])))
        report("P0", max(worst.values()) < mp.mpf(10) ** -bar,
               f"hierarchy on the E-grid against README's Y_1, L_1, Y_2, L_2 at the {used} of {B.N * B.P} nodes with l <= V^2/4: worst |difference|/(1 + |value|) "
               f"{', '.join(f'{k} {mp.nstr(v, 2)}' for k, v in worst.items())} (bar 1e-{bar}; the closed L_2 needs Cl2 and K)")

    # E1 at k = 3, 4, 5: Omega_{k-2} at its root (G_k = 0, so the right side must be 0 to 1e-bar) and shifted by 1/2 (G_k = -2 pi)
    for k in (3, 4, 5):
        res = []
        for shift in (0, mp.mpf(1) / 2):
            om = list(base)
            om[k - 2] = om[k - 2] + shift
            with mp.workdps(dps):
                fl = Fields(B, om, k)
                lhs, rhs, parts = e1_residual(fl, k)
            res.append((shift, lhs, rhs, parts))
        (s0, l0, r0, p0), (s1, l1, r1, p1) = res
        scale = max(abs(x) for x in p0)
        ok = abs(l0 - r0) < mp.mpf(10) ** -bar * scale and digits(l1, r1) > bar
        report(f"E1_{k}", ok,
               f"G_{k} against E1's right side: at the root {mp.nstr(l0, 3)} against {mp.nstr(r0, 3)} (terms of size {mp.nstr(scale, 3)}); "
               f"Omega_{k - 2} + 1/2: {nd(l1, 30)} against {nd(r1, 30)}: {dg(l1, r1)} digits (bar {bar})")

    # E2 at m = 2, 3, 4, Omega_{m-1} at the solved value, at 0 and at 1/2 (an arbitrary constant)
    for m in (2, 3, 4):
        outs = []
        for om_m in (base[m - 1], mp.mpf(0), mp.mpf(1) / 2):
            om = list(base)
            om[m - 1] = om_m
            with mp.workdps(dps):
                fl = Fields(B, om, m)
                lhs, rhs = e2_sides(fl, m)
            outs.append((lhs, rhs))
        ds = [digits(a, b) for a, b in outs]
        report(f"E2_{m}", min(ds) > bar,
               f"int L_{m} dth = {nd(outs[0][0], 30)} against 2 int h_{m} ln cos/sin^2 dth = {nd(outs[0][1], 30)}: "
               f"{', '.join(f'{float(d):.1f}' for d in ds)} digits at Omega_{m - 1} = solved, 0, 1/2 (bar {bar})")

    # E3: the by-parts identity for Y_m, m = 2, 3, 4 (k = m + 1), and G_k with no Y_{k-1}, k = 3, 4, 5
    for k in (3, 4, 5):
        m = k - 1
        om = list(base)
        om[k - 2] = om[k - 2] + mp.mpf(1) / 2
        with mp.workdps(dps):
            fl = Fields(B, om, k)
            lhs, rhs = e3_sides(fl, m)
            gr, parts = g_reduced(fl, k)
        report(f"E3_{k}", digits(lhs, rhs) > bar,
               f"int Y_{m} [(ell - Omega0)/sin^2 + 1] dth = {nd(lhs, 30)} against 2 int (ell - Omega0) F_{m - 1} dth = {nd(rhs, 30)}: "
               f"{dg(lhs, rhs)} digits (bar {bar})")
        report(f"E3G_{k}", digits(fl.G[k], gr) > bar,
               f"G_{k} (Omega_{k - 2} + 1/2) = {nd(fl.G[k], 30)} against the reduced form with no Y_{m}, no Y_{k}: {nd(gr, 30)}: "
               f"{dg(fl.G[k], gr)} digits (bar {bar}); its four terms {', '.join(mp.nstr(t, 8) for t in parts)}")

    # boundary terms, at the first and the last node of the grid
    with mp.workdps(dps):
        fl = Fields(B, base, 5)
        p1, i1 = 0, 0
        pN, iN = B.P - 1, B.N - 1
        th1, thN = mp.atan(B.u[p1][i1]), mp.atan(B.u[pN][iN])
        xN = mp.atan(1 / B.u[pN][iN])
        rows = []

        def Cm(m, p, i):
            return -(fl.L[m][p][i] + 2 * (base[m - 1] if m else mp.mpf(1) / 2)) / fl.u[p][i]    # L_m = -2 Omega_{m-1} + u C_m

        def P(p, i):
            return (fl.ell[p][i] - base[0]) / fl.u[p][i]
        for k in (3, 4, 5):
            rows.append((f"Y_{k} cot", lambda p, i, k=k: fl.Y[k][p][i] / fl.u[p][i], 1, -2 * base[k - 2], k + 1))
        for m in (2, 3, 4):
            rows.append((f"ell C_{m}", lambda p, i, m=m: fl.ell[p][i] * Cm(m, p, i), 3, None, m + 1))
        rows.append(("cot ell", lambda p, i: fl.ell[p][i] / fl.u[p][i], 1, mp.mpf(-1) / 2, 1))
        rows.append(("Y_4 P", lambda p, i: fl.Y[4][p][i] * P(p, i), 1, -base[0] * (-2 * base[2]), 6))
        ok_t, ok_o = True, True
        lines = []
        th2 = mp.atan(B.u[0][1])
        flat = [(p, i) for p in range(B.P) for i in range(B.N)]
        def near(lt):
            return min(flat, key=lambda pi_: abs(B.l[pi_[0]][pi_[1]] - lt))
        sel = [near(lt) for lt in (30, 60, 100)]
        for name, fn, pw, lead, q in rows:
            b1, b2, bN = fn(0, 0), fn(0, 1), fn(pN, iN)
            r1, r2 = b1 / th1 ** pw, b2 / th2 ** pw
            stable = abs(r1 - r2) < mp.mpf(1) / 100 * abs(r1)         # b / th^pw tends to a constant at the throat
            lead_ok = True if lead is None else abs(r1 - lead) < mp.mpf(1) / 100 * abs(lead)
            outer = abs(bN) < mp.mpf(10) ** -bar
            rat = [fn(p, i) / (mp.atan(1 / B.u[p][i]) * B.l[p][i] ** q) for (p, i) in sel]   # b / (x l^q), bounded if b = O(x l^q)
            bounded = max(abs(r) for r in rat) < 10 * (min(abs(r) for r in rat) + mp.mpf(1) / 100)
            ok_t, ok_o = ok_t and stable and lead_ok, ok_o and outer and bounded
            lines.append(f"{name}: /th^{pw} {mp.nstr(r1, 6)}, {mp.nstr(r2, 6)}; outer {mp.nstr(bN, 2)}, b/(x l^{q}) at l = 30, 60, 100: "
                         f"{', '.join(mp.nstr(r, 3) for r in rat)}")
        report("BND", ok_t and ok_o, f"boundary terms: throat end b/th^p at the first two nodes th = {mp.nstr(th1, 3)}, {mp.nstr(th2, 3)} "
                                     f"(constant to 1%, and = the leading coefficient -2 Omega_(k-2) th where known), pi/2 end at the last node "
                                     f"x = pi/2 - th = {mp.nstr(xN, 3)} (bar |b| < 1e-{bar}, and b/(x l^q) bounded along the grid): " + "; ".join(lines))
        # growth: local exponents d ln|f| / d ln l between l = 30 and l = 60 (Y_k, L_k) and Y_k at the last node over l^k
        n30, n60 = near(30), near(60)
        l30, l60 = B.l[n30[0]][n30[1]], B.l[n60[0]][n60[1]]
        ex = []
        for k in range(1, 6):
            y30, y60 = fl.Y[k][n30[0]][n30[1]], fl.Y[k][n60[0]][n60[1]]
            ex.append(f"Y_{k} {mp.nstr(mp.log(abs(y60 / y30)) / mp.log(l60 / l30), 3)}")
        for k in range(1, 5):
            y30, y60 = fl.L[k][n30[0]][n30[1]], fl.L[k][n60[0]][n60[1]]
            ex.append(f"L_{k} {mp.nstr(mp.log(abs(y60 / y30)) / mp.log(l60 / l30), 3)}")
        lN = B.l[pN][iN]
        yl = ', '.join(f"Y_{k}/l^{k} = {mp.nstr(fl.Y[k][pN][iN] / lN ** k, 4)}" for k in range(1, 6))
        report("GROW", True, f"growth (reported, no bar): local exponents between l = {mp.nstr(l30, 4)} and {mp.nstr(l60, 4)}: " + ", ".join(ex)
                             + f"; at the last node l = {mp.nstr(lN, 5)}: " + yl)
    print(f"  exch: {time.time() - t0:.0f} s", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# grid: Omega_2 and Omega_3 on the third method-B grid
# ---------------------------------------------------------------------------------------------------------------
def read_a70_om2():
    txt = open(ORDER4_OUT).read()
    return mp.mpf(re.search(r"Omega2 \(A70\) = (\d+\.\d+)", txt).group(1))


def read_b4_om3():
    lines = open(ORDER5_OUT).read().splitlines()
    for j, line in enumerate(lines):
        if line.strip().startswith("B4:"):
            for t in lines[j + 1:j + 8]:
                m = re.match(r"\s+Omega3 = (\S+)", t)
                if m:
                    return mp.mpf(m.group(1))
    raise RuntimeError("B4's Omega3 not found in " + ORDER5_OUT)


def check_grid():
    t0 = time.time()
    N, w, V, dps = GRID3
    mp.mp.dps = dps
    o5.MUTATE = {"grid_coarse": "coarse_B", "grid_short": "short_tail", "grid_wrongh": "wrong_h", "grid_wrongOm1": "wrong_Om1"}.get(MUTATE)
    a70, b4 = read_a70_om2(), read_b4_om3()
    e = o5.run_B(N, w, V, dps, kmax=5, log=True)
    mp.mp.dps = dps
    cl = o5.closed_omegas()
    print(f"  {e['grid']}  {e['secs']:.0f} s", flush=True)
    for j in range(4):
        print(f"    Omega{j} = {nd(e[f'Om{j}'], min(dps, 110))}")
    bar = 40
    d_a70, d_b4 = digits(a70, e["Om2"]), digits(b4, e["Om3"])
    report("GR1", d_a70 > bar, f"Omega_2 on the third grid {nd(e['Om2'], 50)} against #26's A70 {nd(a70, 50)} (order4_output.txt, 48 printed digits): "
                               f"{float(d_a70):.1f} digits (bar {bar})")
    report("GR2", d_b4 > bar, f"Omega_3 on the third grid {nd(e['Om3'], 50)} against #208's B4 {nd(b4, 50)} (order5_output.txt, 113 printed digits): "
                              f"{float(d_b4):.1f} digits (bar {bar})")
    d0, d1 = digits(cl[0], e["Om0"]), digits(cl[1], e["Om1"])
    report("GR3", min(d0, d1) > bar, f"planted: Omega_0, Omega_1 against -4 ln 2 and 2 pi^2/3 + 8 ln^2 2: {float(d0):.1f} and {float(d1):.1f} digits (bar {bar})")
    sl = [digits(-4 * mp.pi, e[f"slope{k}"]) for k in range(2, 6)]
    report("GR4", min(sl) > bar, f"slope of G_k in Omega_(k-2), k = 2..5, against -4 pi: {', '.join(f'{float(x):.1f}' for x in sl)} digits (bar {bar})")
    d2 = digits(cl[2], e["Om2"])
    om3 = (22 * mp.pi ** 4 / 45 + mp.mpf(16) / 3 * mp.pi ** 2 * mp.ln2 ** 2 + mp.mpf(32) / 3 * mp.ln2 ** 4 - 320 * mp.zeta(3) * mp.ln2 - 384 * mp.zeta(3))
    d3 = digits(om3, e["Om3"])
    print(f"  reported: Omega_2 against 80 zeta(3) - (8/3) pi^2 ln 2 - (32/3) ln^3 2: {float(d2):.1f} digits; "
          f"Omega_3 against #208's closed form: {float(d3):.1f} digits", flush=True)
    o5.MUTATE = None
    print(f"  grid: {time.time() - t0:.0f} s", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# ctlsum: the controls' logs against their declared targets
# ---------------------------------------------------------------------------------------------------------------
# name -> (checks that must FAIL, checks declared silent (must PASS: a failure there is reported as off target)).
# Declared from the development runs dev_ctl_*, dev2_ctl_*, dev3_ctl_* (raw I209-rpi5), before the declared runs.
_E1 = {"E1_3", "E1_4", "E1_5"}
_E2 = {"E2_2", "E2_3", "E2_4"}
_E3 = {"E3_3", "E3_4", "E3_5"}
_E3G = {"E3G_3", "E3G_4", "E3G_5"}
TARGETS = {
    "sym_S1": ({"S1"}, set()), "sym_S2": ({"S2"}, set()), "sym_c1": ({"S3", "S5"}, set()), "sym_P": ({"S4"}, set()),
    "sym_S6": ({"S6"}, set()), "sym_S7": ({"S7"}, set()), "sym_S8": ({"S8"}, {"S9"}), "sym_S9": ({"S9"}, set()),
    "sym_S10": ({"S10"}, set()), "sym_S11": ({"S11"}, set()),
    "exact_K": ({"P0b", "P2_2"}, {"P0c", "P1", "P2_1"}),
    "exact_L1": ({"P0c", "P1", "P2_1", "P2_2"}, {"P0b"}),
    "exact_h2": ({"P0c"}, {"P1", "P2_1", "P2_2"}),
    "e1_sign": ({"P1"} | _E1, _E2 | _E3G), "e1_factor": ({"P1"} | _E1, _E2 | _E3G),
    "e2_weight": ({"P2_1", "P2_2"} | _E2, _E1 | _E3G),
    # the grid's L_m = -2 Omega - u (tail) has no growing term for any Omega_0, so the identities stay true for it; the shift
    # shows at the throat (C_2(0) = -G_2) and in the closed L_2
    "bnd_G": ({"P0", "BND"}, _E1 | _E2 | _E3 | _E3G),
    "bnd_Y0": ({"P0", "BND"} | _E1 | _E2 | _E3 | _E3G, set()),
    # the by-parts identity holds for any constant in place of Omega_0; only the combination with c_1 fixes it
    "e3_Om0": (_E3G, _E3 | _E1 | _E2), "e3_factor": (_E3G, _E3 | _E1 | _E2),
    "cl_moment": ({"CM"}, {"CL1", "CL1a", "CL2"}),
    "cl_arg": ({"CL1", "CL1a", "CL2"}, {"CM"}),
    "cl_arg2": ({"CL1", "CL1a"}, {"CM", "CL2"}),          # sin(2n th) orthogonality does not see the sign (-1)^(n+1)
    "grid_coarse": ({"GR1", "GR2", "GR3", "GR4"}, set()), "grid_short": ({"GR1", "GR2", "GR3", "GR4"}, set()),
    "grid_wrongh": ({"GR1", "GR2", "GR3"}, {"GR4"}), "grid_wrongOm1": ({"GR1", "GR2"}, {"GR3", "GR4"}),
}
assert set(TARGETS) == set(MUTATIONS)


def check_ctlsum():
    n_on = 0
    for name in MUTATIONS:
        path = os.path.join(CTLDIR, f"{CTLPREFIX}_{name}.log")
        if not os.path.exists(path):
            report(f"CT_{name}", False, f"no log {path}")
            continue
        txt = open(path).read()
        fired = set(re.findall(r"^FAIL\s+(\S+):", txt, re.M))
        ran = re.findall(r"^(?:PASS|FAIL)\s+(\S+):", txt, re.M)
        fire, silent = TARGETS[name]
        missing = fire - fired
        off = silent & fired
        extra = fired - fire - silent
        ok = not missing and not off and "exit=1" in txt
        n_on += ok
        report(f"CT_{name}", ok, f"fired {sorted(fired)} of {len(ran)} checks; declared targets {sorted(fire)}"
                                 + (f"; MISSING {sorted(missing)}" if missing else "")
                                 + (f"; declared silent but fired {sorted(off)}" if off else "")
                                 + (f"; also fired (not declared) {sorted(extra)}" if extra else "")
                                 + (f"; declared silent and silent {sorted(silent)}" if silent and not off else ""))
    print(f"  controls on target: {n_on} of {len(MUTATIONS)}", flush=True)


class Tee:
    def __init__(self, f):
        self.f, self.o = f, sys.__stdout__

    def write(self, t):
        self.o.write(t)
        self.f.write(t)

    def flush(self):
        self.o.flush()
        self.f.flush()


if __name__ == "__main__":
    if OUTFILE:
        sys.stdout = Tee(open(OUTFILE, "w"))
    mp.mp.dps = 60
    print(f"exchange_review.py  mutate={MUTATE}  parts={PARTS}  mpmath {mp.__version__}  sympy {sp.__version__}", flush=True)
    t0 = time.time()
    if "sym" in PARTS:
        check_sym()
    if "exact" in PARTS:
        check_exact()
    if "clausen" in PARTS:
        check_clausen()
    if "exch" in PARTS:
        check_exch()
    if "grid" in PARTS:
        check_grid()
    if "ctlsum" in PARTS:
        check_ctlsum()
    nfail = sum(1 for _, ok in RESULTS if not ok)
    print(f"\n{len(RESULTS) - nfail} passed, {nfail} failed of {len(RESULTS)} checks; {time.time() - t0:.0f} s")
    sys.exit(1 if nfail else 0)

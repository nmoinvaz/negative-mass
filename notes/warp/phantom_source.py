#!/usr/bin/env python3
"""Can a minimally coupled scalar with the phantom sign source the Alcubierre metric? (issues #57 and #61)

A scalar with T_ab = eps (d_a phi d_b phi - g_ab (d phi)^2 / 2) - g_ab V(phi) and G_ab = 8 pi T_ab gives
R_ab = 8 pi (eps d_a phi d_b phi + V g_ab). So R_ab - lam g_ab must be 8 pi eps u_a u_b for some lam (= 8 pi V), the
covector u must be closed (u = d phi), and dV ^ d phi = 0. The script tests these conditions on Alcubierre's metric

    ds^2 = -dt^2 + (dx - v f(r_s) dt)^2 + dy^2 + dz^2,   r_s^2 = (x - v t)^2 + y^2 + z^2,   v constant,

for a general profile f, and on the D = 4 Ellis-Bronnikov wormhole of paper/3-wormhole-largeD (the control). Section 6
(issue #61) checks the identities behind the test for an accelerating bubble, b = -v_s(t) f(r_s) with
r_s^2 = (x - x_s(t))^2 + y^2 + z^2 and v_s, x_s independent free functions, and for any shift b(t, x, rho).

    python3 notes/warp/phantom_source.py                 all checks, then every control (exit 1 if any check fails)
    python3 notes/warp/phantom_source.py --mutate=NAME   one control (its target check must fail, exit 1)
    python3 notes/warp/phantom_source.py --no-controls   the checks only

Routes. (A) The Einstein tensor in the lab coordinates (t, x, y, z), with f(r_s) differentiated by the chain rule
on its jets F_k = f^(k)(r_s) as free symbols. (B) The 3+1 Hamiltonian constraint with K_ij = (d_i b_j + d_j b_i)/2.
(C) The Ricci tensor in comoving spherical coordinates (t, r, theta, varphi) around the axis of motion, with the
shift w(r) = v (1 - f(r)) a sympy Function. A and C share the Christoffel/Ricci routine `ricci`, B shares only the
chain-rule derivative `D` with A. All three use sympy.

Section 6 has two more. (D) `ricci` in the lab coordinates, with the jets of v_s, x_s and f, or of b(t, x, rho), as
free symbols and the total derivative `JetD`. (E) `ricci_at_point`, the Christoffel formula in exact rationals at a
point, with g and its first and second derivatives from sympy's own diff of explicit functions. D and E share only the
metric's definition `shift_metric`. E is the planted truth for D's time derivatives.

Expressions in (X = x - v t, y, z, R = r_s) are reduced modulo R^2 - X^2 - y^2 - z^2. That relation is irreducible,
so a remainder of zero is exact vanishing on the spacetime. In section 6, X = x - x_s(t), and expressions in rho are
reduced modulo rho^2 - y^2 - z^2.
"""
import random
import subprocess
import sys
import time

import mpmath as mp
import sympy as sp

# name: (checks that must fail, what the mutation does)
MUTATIONS = {
    "rho-v2": (["rho_einstein", "rho_3p1"], "doubles v in the expected energy density"),
    "xcheck-v": (["ricci_xcheck"], "doubles v in the spherical-coordinate Ricci tensor only"),
    "alc-drop-w2": (["identity_sph", "identity_lab", "identity_axis"],
                    "drops the shift's square from g_tt in both coordinate systems"),
    "wall-v2": (["wall_curvature"], "doubles v in the expected wall curvature"),
    "alc-planted": (["alc_points_fail"], "replaces the Ricci tensor at the sample points by a planted scalar form"),
    "eb-eps+1": (["eb_sign"], "tests Ellis-Bronnikov with eps = +1"),
    "eb-deform": (["eb_algebraic"], "deforms Ellis-Bronnikov's r(l)^2 = l^2 + a^2 by a^4/(l^2 + a^2)"),
    "eb-unclosed": (["eb_closed"], "multiplies the extracted gradient u by 1 + cos(theta)^2"),
    "eb-potential": (["eb_potential"], "adds t/(8 pi a^3) to V"),
    "eb-residual": (["eb_residual"], "scales phi by 11/10 before substituting it back"),
    "accel-static": (["gen_xcheck", "accel_xcheck"], "drops every t-derivative of b, v_s and x_s in section 6's jets"),
    "accel-lapse": (["gen_I1", "gen_I2", "gen_Rrr", "accel_I1", "accel_I2", "accel_axis", "accel_wall"],
                    "gives section 6's slicing the lapse 1 + b, not 1"),
    "accel-rotate": (["gen_I1", "gen_I2", "gen_Rrr", "accel_I1", "accel_I2", "accel_wall"],
                     "rotates section 6's (e_rho, e_varphi) by alpha, tan(alpha/2) = 1/20"),
    "accel-wall-v2": (["gen_Rrr", "accel_wall"], "doubles d_rho b and v_s in section 6's expected wall curvature"),
}

RESULTS = []


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def ricci(g, X, Dfun, simp):
    """Ricci tensor R_bc of the metric g in coordinates X, with Dfun(E, a) = d_a E."""
    n = len(X)
    gi = g.inv().applyfunc(simp)
    dg = [[[Dfun(g[b, c], a) for c in range(n)] for b in range(n)] for a in range(n)]  # d_a g_bc
    Gam = [[[simp(sum(gi[a, d] * (dg[b][d][c] + dg[c][d][b] - dg[d][b][c]) for d in range(n)) / 2)
             for c in range(n)] for b in range(n)] for a in range(n)]
    Ric = sp.zeros(n, n)
    for b in range(n):
        for c in range(b, n):
            s = 0
            for a in range(n):
                s += Dfun(Gam[a][b][c], a) - Dfun(Gam[a][a][b], c)
                for d in range(n):
                    s += Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][a][b]
            Ric[b, c] = Ric[c, b] = simp(s)
    return Ric


# ---------------------------------------------------------------------------------------------------------------
# Alcubierre in lab coordinates (route A) and its 3+1 energy density (route B)
# ---------------------------------------------------------------------------------------------------------------
X, y, z, R, v = sp.symbols("X y z R v", real=True)
F = sp.symbols("F0:5", real=True)  # F[k] = f^(k)(r_s)
Q = X**2 + y**2 + z**2  # = R^2


def D(E, a):
    """d/dx^a in lab coordinates (t, x, y, z). E depends on t and x only through X = x - v t."""
    if a == 0:
        return sp.expand(-v * D(E, 1))
    xa = (None, X, y, z)[a]
    dR = sp.diff(E, R) + sum(F[k + 1] * sp.diff(E, F[k]) for k in range(4))
    return sp.expand(sp.diff(E, xa) + xa / R * dR)


def reduce_R(E, rel=Q, var=R):
    """Remainder of E's numerator modulo var^2 - rel, as an expression (zero iff E vanishes)."""
    num, _ = sp.fraction(sp.together(sp.expand(E)))
    return sp.expand(sp.rem(sp.expand(num), var**2 - rel, var))


def lab_tetrad_ricci(mut):
    gtt = -1 if mut == "alc-drop-w2" else -1 + v**2 * F[0]**2
    g = sp.Matrix([[gtt, -v * F[0], 0, 0], [-v * F[0], 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    Ric = ricci(g, [0, 1, 2, 3], D, sp.expand)
    # Orthonormal frame dual to e^0 = dt, e^1 = dx - v f dt, e^2 = dy, e^3 = dz. E_0 is the Eulerian normal n.
    E = sp.Matrix([[1, 0, 0, 0], [v * F[0], 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    return (E.T * Ric * E).applyfunc(sp.expand)


def rho_3p1():
    """Eulerian energy density from the Hamiltonian constraint 16 pi rho = K^2 - K_ij K^ij (flat slices, lapse 1)."""
    beta = [-v * F[0], 0, 0]  # beta_i = beta^i, gamma_ij = delta_ij
    K = sp.Matrix(3, 3, lambda i, j: (D(beta[j], i + 1) + D(beta[i], j + 1)) / 2)
    trK = K.trace()
    return sp.expand((trK**2 - sum(K[i, j]**2 for i in range(3) for j in range(3))) / (16 * sp.pi))


# ---------------------------------------------------------------------------------------------------------------
# Alcubierre in comoving spherical coordinates (route C)
# ---------------------------------------------------------------------------------------------------------------
t_, r_, th, vp = sp.symbols("t r theta varphi", real=True)
W = sp.symbols("W0:3", real=True)  # W[k] = w^(k)(r), w = v (1 - f)


def sph_tetrad_ricci(mut):
    w = sp.Function("w")(r_)
    e = sp.zeros(4, 4)  # e[a, mu]: e^0 = dt, e^1 = d xi + w dt, e^2 = d rho_cyl, e^3 = rho_cyl d varphi
    e[0, 0] = 1
    e[1, :] = sp.Matrix([[w, sp.cos(th), -r_ * sp.sin(th), 0]])
    e[2, :] = sp.Matrix([[0, sp.sin(th), r_ * sp.cos(th), 0]])
    e[3, :] = sp.Matrix([[0, 0, 0, r_ * sp.sin(th)]])
    eta = sp.diag(-1, 1, 1, 1)
    g = (e.T * eta * e).applyfunc(sp.simplify)
    if mut == "alc-drop-w2":
        g[0, 0] = -1
    Xs = [t_, r_, th, vp]
    Ric = ricci(g, Xs, lambda E, a: sp.diff(E, Xs[a]), sp.simplify)
    Rt = (e.inv().T * Ric * e.inv()).applyfunc(sp.simplify)
    jets = [(w.diff(r_, 2), W[2]), (w.diff(r_), W[1]), (w, W[0])]
    return Rt.applyfunc(lambda E: sp.expand(sp.simplify(E.subs(jets))))


# ---------------------------------------------------------------------------------------------------------------
# The scalar-form test, shared by Alcubierre's sample points and the Ellis-Bronnikov control
# ---------------------------------------------------------------------------------------------------------------
lam = sp.Symbol("lam")


def scalar_form(Rm, gm, eps, simp=sp.simplify):
    """Test R_ab = lam g_ab + 8 pi eps u_a u_b at the algebraic level.

    Returns (status, lam*, M, u) with status 'algebraic' if no lam makes R - lam g of rank <= 1, 'sign' if the rank-1
    remainder has the other sign, 'trivial' if R = lam g (u = 0), or 'pass'. The candidates lam* are the roots of the
    gcd over all 2x2 minors of R - lam g, each a polynomial of degree <= 2 in lam.
    """
    n = Rm.shape[0]
    A = Rm - lam * gm
    polys = []
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(n):
                for m in range(k + 1, n):
                    num, _ = sp.fraction(sp.together(A[i, k] * A[j, m] - A[i, m] * A[j, k]))
                    num = sp.expand(num)
                    if num != 0:
                        polys.append(sp.Poly(num, lam))
    if not polys:
        return "trivial", 0, sp.zeros(n, n), None
    gcd = polys[0]
    for p in polys[1:]:
        gcd = sp.gcd(gcd, p)
    if gcd.degree() < 1:
        return "algebraic", None, None, None
    roots = sp.roots(gcd.as_expr(), lam)
    for lr in roots:
        if lr.is_real is False:
            continue
        M = (Rm - lr * gm).applyfunc(simp)
        if all(sp.simplify(x) == 0 for x in M):
            return "trivial", lr, M, None
        k = next(i for i in range(n) if sp.simplify(M[i, i]) != 0)
        c = simp(M[k, k] / (8 * sp.pi * eps))
        pos = c.is_positive
        if pos is None:  # undecided by assumptions: sample the free symbols
            free = sorted(c.free_symbols, key=str)
            vals = [c.subs({s: sp.Rational(random.randint(1, 9), random.randint(1, 9)) for s in free}) for _ in range(5)]
            pos = all(sp.N(x) > 0 for x in vals)
        if not pos:
            return "sign", lr, M, None
        u = (M[:, k] / sp.sqrt(8 * sp.pi * eps * M[k, k])).applyfunc(simp)
        return "pass", lr, M, u
    return "algebraic", None, None, None


# ---------------------------------------------------------------------------------------------------------------
# Section 6 (issue #61): an accelerating bubble and a general axisymmetric shift, lab coordinates (t, x, y, z)
# ---------------------------------------------------------------------------------------------------------------
VS = sp.symbols("vs0:4", real=True)  # VS[k] = v_s^(k)(t)
XS = sp.symbols("xs0:5", real=True)  # XS[k] = x_s^(k)(t) for k >= 1. In section 6, X = x - x_s(t).
RHO = sp.Symbol("rho", positive=True)  # rho = (y^2 + z^2)^(1/2), in the route for a general b(t, x, rho)
BJ = {}  # BJ[(i, j, k)] is the jet d_t^i d_x^j d_rho^k b, and BJ[that symbol] = (i, j, k)
TX = sp.symbols("t x", real=True)  # lab t and x, for the explicit functions of route E


def bjet(i, j, k):
    if (i, j, k) not in BJ:
        BJ[(i, j, k)] = s = sp.Symbol(f"b{i}{j}{k}", real=True)
        BJ[s] = (i, j, k)
    return BJ[(i, j, k)]


class JetD:
    """Total derivative d/dx^a in (t, x, y, z) of an expression in y, z and jet symbols. rule(s) lists d_a s."""

    def __init__(self, rule):
        self.rule = rule

    def __call__(self, E, a):
        out = sp.diff(E, (y, z)[a - 2]) if a >= 2 else 0
        for s in E.free_symbols:
            ds = self.rule(s)
            if ds is not None and ds[a] != 0:
                out += sp.diff(E, s) * ds[a]
        return sp.expand(out)


def bubble_rule(static):
    """Jets of b = -v_s(t) f(r_s), r_s^2 = X^2 + y^2 + z^2, X = x - x_s(t), with v_s and x_s independent."""
    dt = 0 if static else 1

    def rule(s):
        if s == X:
            return [-dt * XS[1], 1, 0, 0]
        if s == R:
            return [-dt * X * XS[1] / R, X / R, y / R, z / R]
        if s in F:
            return [F[F.index(s) + 1] * d for d in rule(R)]
        if s in VS:
            return [dt * VS[VS.index(s) + 1], 0, 0, 0]
        if s in XS[1:]:
            return [dt * XS[XS.index(s) + 1], 0, 0, 0]
        return None
    return rule


def general_rule(static):
    """Jets of a general b(t, x, rho), with d_y rho = y / rho and d_z rho = z / rho."""
    def rule(s):
        if s == RHO:
            return [0, 0, y / RHO, z / RHO]
        if s in BJ:
            i, j, k = BJ[s]
            up = bjet(i, j, k + 1)
            return [0 if static else bjet(i + 1, j, k), bjet(i, j + 1, k), up * y / RHO, up * z / RHO]
        return None
    return rule


def shift_metric(b, lapse):
    """ds^2 = -N^2 dt^2 + (dx + b dt)^2 + dy^2 + dz^2, N = 1 (or 1 + b for a control), and the orthonormal frame
    E_0 = (d_t - b d_x) / N, E_1 = d_x, E_2 = d_y, E_3 = d_z as the columns of E. E_0 is the Eulerian normal n."""
    N = 1 + b if lapse else sp.Integer(1)
    g = sp.Matrix([[-N**2 + b**2, b, 0, 0], [b, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    E = sp.Matrix([[1 / N, 0, 0, 0], [-b / N, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    return g, E


def shift_tetrad_ricci(b, rule, lapse):
    """Route D: frame Ricci tensor of shift_metric, with b a jet symbol expression."""
    g, E = shift_metric(b, lapse)
    Ric = ricci(g, [0, 1, 2, 3], JetD(rule), sp.expand)
    return (E.T * Ric * E).applyfunc(sp.expand)


def ricci_at_point(bexpr, p, lapse):
    """Route E: frame Ricci tensor of shift_metric at the point p = (t, x, y, z), for an explicit b(t, x, y, z)."""
    cs = (TX[0], TX[1], y, z)
    at = dict(zip(cs, p))

    def val(e):
        r = sp.cancel(e.subs(at))
        assert r.is_Rational, r
        return r
    n = 4
    g, E = shift_metric(bexpr, lapse)
    dg = [[[val(sp.diff(g[b, c], cs[a])) for c in range(n)] for b in range(n)] for a in range(n)]  # d_a g_bc
    ddg = [[[[val(sp.diff(g[b, c], cs[a], cs[e])) for c in range(n)] for b in range(n)] for e in range(n)]
           for a in range(n)]  # d_a d_e g_bc
    gi = g.applyfunc(val).inv()
    dgi = [-gi * sp.Matrix(n, n, lambda b, c: dg[a][b][c]) * gi for a in range(n)]  # d_a g^-1
    Gam = [[[sum(gi[k, d] * (dg[b][d][c] + dg[c][d][b] - dg[d][b][c]) for d in range(n)) / 2
             for c in range(n)] for b in range(n)] for k in range(n)]  # Gam[k][b][c] = Gamma^k_bc
    dGam = [[[[sum(dgi[a][k, d] * (dg[b][d][c] + dg[c][d][b] - dg[d][b][c])
                   + gi[k, d] * (ddg[a][b][d][c] + ddg[a][c][d][b] - ddg[a][d][b][c]) for d in range(n)) / 2
               for c in range(n)] for b in range(n)] for k in range(n)] for a in range(n)]  # d_a Gamma^k_bc
    Ric = sp.Matrix(n, n, lambda b, c: sum(dGam[a][a][b][c] - dGam[c][a][a][b] for a in range(n))
                    + sum(Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][a][b]
                          for a in range(n) for d in range(n)))
    E0 = E.applyfunc(val)
    return (E0.T * Ric * E0).applyfunc(sp.cancel)


def cylinder(Rt, tau):
    """rho^2 R_rho,rho, rho^2 R_rho,varphi, rho^2 R_varphi,varphi and [rho R(e_varphi, E_0), rho R(e_varphi, E_1)],
    with (e_rho, e_varphi) rotated by alpha, tan(alpha/2) = tau. tau = 0 is the true frame."""
    tau = sp.Rational(tau)  # exact, so tau = 0 gives the integers 1 and 0, not floats
    ca, sa = (1 - tau**2) / (1 + tau**2), 2 * tau / (1 + tau**2)
    P = sp.Matrix([[ca, sa], [-sa, ca]]) * sp.Matrix([[y, z], [-z, y]])  # rows: rho e_rho, rho e_varphi
    T = Rt[2:, 2:]
    rr, rp, pp = [sp.expand((P[i, :] * T * P[j, :].T)[0]) for i, j in ((0, 0), (0, 1), (1, 1))]
    return rr, rp, pp, [sp.expand(P[1, 0] * Rt[c, 2] + P[1, 1] * Rt[c, 3]) for c in (0, 1)]


# ---------------------------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------------------------
def tanh_jets(r0, sig, Rb):
    """Jets (f, f', f'') of Alcubierre's profile f = [tanh(sig (r + Rb)) - tanh(sig (r - Rb))] / (2 tanh(sig Rb))."""
    mp.mp.dps = 50
    fn = lambda r: (mp.tanh(sig * (r + Rb)) - mp.tanh(sig * (r - Rb))) / (2 * mp.tanh(sig * Rb))
    return [sp.Rational(mp.nstr(mp.diff(fn, mp.mpf(r0.p) / r0.q, k), 45)) for k in range(3)]


def bump_jets(r0, Rb):
    """Jets of the C^2 bump f = (1 - (r/Rb)^2)^3 for r < Rb, 0 beyond. Exact."""
    rr = sp.Symbol("rr")
    fb = (1 - (rr / Rb)**2)**3
    return [sp.diff(fb, rr, k).subs(rr, r0) if r0 < Rb else 0 for k in range(3)]


def main():
    opts = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
    mut = opts.get("mutate")
    if mut is not None and mut not in MUTATIONS:
        sys.exit(f"unknown mutation {mut!r}. Known: {', '.join(MUTATIONS)}")
    random.seed(57)
    t0 = time.time()
    print(f"phantom_source.py, issues #57 and #61. sympy {sp.__version__}, mpmath {mp.__version__}, python "
          f"{sys.version.split()[0]}. mutation: {mut or 'none'}" + (f" ({MUTATIONS[mut][1]})" if mut else ""))
    vexp = 2 * v if mut in ("rho-v2",) else v

    # 1. Eulerian energy density, routes A and B, against -(v^2 / 32 pi) ((y^2 + z^2) / r_s^2) f'^2.
    print("\n== 1. Eulerian energy density ==")
    Rt = lab_tetrad_ricci(mut)
    Rs = sp.expand(-Rt[0, 0] + Rt[1, 1] + Rt[2, 2] + Rt[3, 3])
    rho_E = sp.expand((Rt[0, 0] + Rs / 2) / (8 * sp.pi))  # G_ab n^a n^b / 8 pi, n = E_0
    rho_exp = -vexp**2 / (32 * sp.pi) * (y**2 + z**2) / R**2 * F[1]**2
    rem = reduce_R(rho_E - rho_exp)
    check("rho_einstein", rem == 0, f"G_ab n^a n^b / 8 pi - expected reduces to {sp.factor(rem)}"
          f" ({time.time() - t0:.0f} s)")
    rem = reduce_R(rho_3p1() - rho_exp)
    check("rho_3p1", rem == 0, f"(K^2 - K_ij K^ij) / 16 pi - expected reduces to {sp.factor(rem)}")

    # 2. The Ricci tensor in comoving spherical coordinates, and its cross-check against the lab frame.
    print("\n== 2. Ricci tensor in the tetrad (E_0 = n, E_1 = d_x, e_rho, e_varphi), comoving spherical ==")
    Rsph = sph_tetrad_ricci(mut)
    names = ["0", "1", "rho", "varphi"]
    for a in range(4):
        for b in range(a, 4):
            print(f"  R_({names[a]},{names[b]}) = {sp.factor(Rsph[a, b])}")
    vsph = 2 * v if mut == "xcheck-v" else v
    pts = [(12, 3, 4, 13), (-12, 4, 3, 13), (4, 3, 0, 5), (-3, 0, 4, 5), (8, 9, 12, 17), (0, 7, 24, 25),
           (7, 24, 0, 25), (-15, 12, 16, 25)]
    bad = 0
    ntest = 0
    for (x0, y0, z0, r0) in pts:
        for scale in (sp.Rational(1, 13), sp.Rational(3, 5), 1):
            Xv, yv, zv, Rv = [sp.Rational(c) * scale for c in (x0, y0, z0, r0)]
            rho_c = sp.sqrt(yv**2 + zv**2)
            jets = [sp.Rational(random.randint(-9, 9), random.randint(1, 9)) for _ in range(3)]
            vv = sp.Rational(random.randint(1, 9), random.randint(1, 9))
            lab = Rt.subs({X: Xv, y: yv, z: zv, R: Rv, v: vv, F[0]: jets[0], F[1]: jets[1], F[2]: jets[2]})
            Rot = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, yv / rho_c, -zv / rho_c], [0, 0, zv / rho_c, yv / rho_c]])
            lab = Rot.T * lab * Rot
            ws = {W[0]: vsph.subs(v, vv) * (1 - jets[0]), W[1]: -vsph.subs(v, vv) * jets[1],
                  W[2]: -vsph.subs(v, vv) * jets[2]}
            sph = Rsph.subs(ws).subs({r_: Rv, th: sp.atan2(rho_c, Xv)})
            diff = (lab - sph).applyfunc(sp.simplify)
            ntest += 1
            bad += any(x != 0 for x in diff)
    check("ricci_xcheck", bad == 0, f"lab-frame and spherical tetrad Ricci differ at {bad} of {ntest} exact points"
          " (random rational v and jets)")

    # 3. The two curvature identities that block the scalar form, for a general profile.
    print("\n== 3. Identities: R(e_varphi, .) = 0 and R_00 + R_11 + 2 R_rho,rho = 0 ==")
    ids = [sp.simplify(Rsph[a, 3]) for a in range(4)] + [sp.simplify(Rsph[0, 0] + Rsph[1, 1] + 2 * Rsph[2, 2])]
    check("identity_sph", all(x == 0 for x in ids), f"comoving spherical, general w(r): {ids}")
    rho2 = y**2 + z**2
    Rpp = y**2 * Rt[2, 2] + 2 * y * z * Rt[2, 3] + z**2 * Rt[3, 3]  # rho^2 R(e_rho, e_rho)
    ids = [reduce_R(-z * Rt[2, a] + y * Rt[3, a]) for a in range(4)] + [reduce_R(rho2 * (Rt[0, 0] + Rt[1, 1]) + 2 * Rpp)]
    check("identity_lab", all(x == 0 for x in ids), f"lab Cartesian, general f, residues mod r_s^2: {ids}")
    ax = Rt.subs({y: 0, z: 0})
    ids = [reduce_R(ax[a, b], X**2) for a in (2, 3) for b in range(4)] + [reduce_R(ax[0, 0] + ax[1, 1], X**2)]
    check("identity_axis", all(x == 0 for x in ids), f"on the axis y = z = 0: R(E_2, .), R(E_3, .), R_00 + R_11 "
          f"reduce to {ids}")
    vw = 2 * v if mut == "wall-v2" else v
    rem1 = reduce_R(Rpp - rho2 * vw**2 * F[1]**2 * rho2 / (2 * R**2))
    rem2 = reduce_R(ax[0, 1] - vw * F[1] / R, X**2)
    check("wall_curvature", rem1 == 0 and rem2 == 0,
          f"R_rho,rho = v^2 f'^2 (y^2 + z^2) / (2 r_s^2) [residue {rem1}] and, on the axis, R_01 = v f' / r_s "
          f"[residue {rem2}]. So R_ab != 0 wherever f'(r_s) != 0.")
    print("  Lemma (by hand, notes/warp/phantom_source.md): a symmetric R_ab with R(e_varphi, .) = 0 and\n"
          "  R_00 + R_11 + 2 R_rho,rho = 0 is of the form lam eta_ab + s u_a u_b (s real, nonzero) only if R_ab = 0.\n"
          "  With wall_curvature, the scalar form fails at every point where f'(r_s) != 0, for either eps and any V.")

    # 4. The same scalar-form test at exact sample points, for explicit profiles and free jets.
    print("\n== 4. Scalar-form test at exact points (same code as the Ellis-Bronnikov control) ==")
    samples = []
    wall = [(12, 3, 4, 13), (-12, 4, 3, 13), (13, 0, 0, 13), (-13, 0, 0, 13), (0, 5, 12, 13)]
    for (x0, y0, z0, r0) in wall:  # points on r_s = 1, inside Alcubierre's wall for sig = 8, Rb = 1
        p = [sp.Rational(c, 13) for c in (x0, y0, z0, r0)]
        for vv in (sp.Integer(1), sp.Rational(1, 2), sp.Integer(3)):
            samples.append(("tanh sig=8 Rb=1", p, vv, tanh_jets(p[3], 8, 1)))
    for (x0, y0, z0, r0) in [(4, 3, 0, 5), (-3, 0, 4, 5), (5, 0, 0, 5), (0, 3, 4, 5)]:
        p = [sp.Rational(c, 5) for c in (x0, y0, z0, r0)]  # r_s = 1 < Rb = 2
        samples.append(("bump Rb=2", p, sp.Rational(7, 4), bump_jets(p[3], 2)))
    for (x0, y0, z0, r0) in pts[:4]:
        p = [sp.Rational(c, 7) for c in (x0, y0, z0, r0)]
        samples.append(("random jets", p, sp.Rational(random.randint(1, 9), 3),
                        [sp.Rational(random.randint(-9, 9), random.randint(1, 9)) for _ in range(3)]))
    eta = sp.diag(-1, 1, 1, 1)
    statuses = {}
    for (lab_, p, vv, jets) in samples:
        sub = {X: p[0], y: p[1], z: p[2], R: p[3], v: vv, F[0]: jets[0], F[1]: jets[1], F[2]: jets[2]}
        Rm = Rt.subs(sub)
        if mut == "alc-planted":  # lam eta + 8 pi eps u u, eps = -1, with random rational lam and u
            uu = sp.Matrix([sp.Rational(random.randint(-9, 9), random.randint(1, 9)) for _ in range(4)])
            Rm = sp.Rational(random.randint(-9, 9), 7) * eta - 8 * sp.pi * uu * uu.T
        for eps in (-1, 1):
            st = scalar_form(Rm, eta, eps)[0]
            statuses[(lab_, eps)] = statuses.get((lab_, eps), []) + [st]
    for k, sts in statuses.items():
        print(f"  {k[0]:>16}, eps = {k[1]:+d}: {len(sts)} points, outcomes {sorted(set(sts))}")
    # On the axis R^a_b restricted to (E_0, E_1) is [[-R_00, -R_01], [R_01, -R_00]], eigenvalues -R_00 +- i R_01.
    for (lab_, p, vv, jets) in samples[6:12:3]:
        sub = {X: p[0], y: p[1], z: p[2], R: p[3], v: vv, F[0]: jets[0], F[1]: jets[1], F[2]: jets[2]}
        ev = [sp.N(e, 6) for e in (eta * Rt.subs(sub)).eigenvals()]
        print(f"  info: eigenvalues of R^a_b on the axis at x - v t = {p[0]}, v = {vv} ({lab_}): {ev}")
    allfail = all(st == "algebraic" for sts in statuses.values() for st in sts)
    check("alc_points_fail", allfail, f"every sample point fails the algebraic condition (no lam gives rank <= 1), "
          f"{sum(len(s) for s in statuses.values())} tests")

    # 5. Control: D = 4 Ellis-Bronnikov with a massless phantom scalar (paper/3-wormhole-largeD, Eq. (1), n = 2).
    print("\n== 5. Control: Ellis-Bronnikov, ds^2 = -dt^2 + dl^2 + r(l)^2 dOmega^2, r^2 = l^2 + a^2 ==")
    t, l, the, ph = sp.symbols("t l theta phi", real=True)
    a = sp.Symbol("a", positive=True)
    Xe = [t, l, the, ph]
    r2 = l**2 + a**2 + (a**4 / (l**2 + a**2) if mut == "eb-deform" else 0)
    ge = sp.diag(-1, 1, r2, r2 * sp.sin(the)**2)
    Re = ricci(ge, Xe, lambda E, k: sp.diff(E, Xe[k]), sp.simplify)
    eps = 1 if mut == "eb-eps+1" else -1
    print(f"  eps = {eps:+d}. Nonzero R_ab: " + ", ".join(f"R_{Xe[i]}{Xe[j]} = {sp.factor(Re[i, j])}"
                                                         for i in range(4) for j in range(i, 4) if Re[i, j] != 0))
    status, lr, M, u = scalar_form(Re, ge, eps)
    check("eb_algebraic", status in ("sign", "pass"), f"triple eigenvalue: status {status}, lam = {lr}")
    if status in ("sign", "pass"):
        check("eb_sign", status == "pass", f"rank-1 remainder M = 8 pi eps u u needs M_kk / (8 pi eps) > 0: {status}")
    if status == "pass":
        if mut == "eb-unclosed":
            u = u * (1 + sp.cos(the)**2)
        du = [sp.simplify(sp.diff(u[j], Xe[i]) - sp.diff(u[i], Xe[j])) for i in range(4) for j in range(i + 1, 4)]
        check("eb_closed", all(x == 0 for x in du), f"u = {list(u)}, du = {du}")
        if all(x == 0 for x in du):
            phi = sp.integrate(u[1], l)  # u has only an l-component here, checked next
            grad_ok = all(sp.simplify(sp.diff(phi, Xe[i]) - u[i]) == 0 for i in range(4))
            V = lr / (8 * sp.pi) + (t / (8 * sp.pi * a**3) if mut == "eb-potential" else 0)
            dVdphi = [sp.simplify(sp.diff(V, Xe[i]) * sp.diff(phi, Xe[j]) - sp.diff(V, Xe[j]) * sp.diff(phi, Xe[i]))
                      for i in range(4) for j in range(i + 1, 4)]
            check("eb_potential", grad_ok and all(x == 0 for x in dVdphi),
                  f"phi = {phi}, d phi = u: {grad_ok}, V = {V}, dV ^ d phi = {dVdphi}")
            phis = phi * sp.Rational(11, 10) if mut == "eb-residual" else phi
            gi = ge.inv()
            dphi = sp.Matrix([sp.diff(phis, s) for s in Xe])
            kin = sp.simplify((dphi.T * gi * dphi)[0])
            T = sp.Matrix(4, 4, lambda i, j: eps * (dphi[i] * dphi[j] - ge[i, j] * kin / 2) - ge[i, j] * V)
            Rsc = sp.simplify(sum(gi[i, j] * Re[i, j] for i in range(4) for j in range(4)))
            res = (Re - ge * Rsc / 2 - 8 * sp.pi * T).applyfunc(sp.simplify)
            check("eb_residual", all(x == 0 for x in res), f"G_ab - 8 pi T_ab with (phi, V) substituted back: "
                  f"{'0' if all(x == 0 for x in res) else list(res)}")

    # 6. Issue #61: the identities for an accelerating bubble and for any axisymmetric shift b(t, x, rho).
    print("\n== 6. Accelerating bubble b = -v_s(t) f(r_s), r_s^2 = (x - x_s(t))^2 + rho^2, and any b(t, x, rho) ==")
    print("  ds^2 = -dt^2 + (dx + b dt)^2 + dy^2 + dz^2, frame (E_0 = n, E_1 = d_x, e_rho, e_varphi). v_s and x_s are\n"
          "  independent free functions, so v_s = x_s' is a special case.")
    static, lapse = mut == "accel-static", mut == "accel-lapse"
    tau = sp.Rational(1, 20) if mut == "accel-rotate" else 0
    Rg = shift_tetrad_ricci(bjet(0, 0, 0), general_rule(static), lapse)
    Ra = shift_tetrad_ricci(-VS[0] * F[0], bubble_rule(static), lapse)
    jets_of = lambda E: sorted((str(s) for s in E.free_symbols if s in BJ or s in VS or s in XS), key=str)
    print(f"  info: jets in R_ab, general b: {jets_of(Rg)}; bubble: {jets_of(Ra)}")
    print(f"  info: jets in R_00, general b: {jets_of(Rg[0, 0])}; bubble: {jets_of(Ra[0, 0])}")

    # Planted truth for route D: explicit functions differentiated by sympy itself (route E), at exact points.
    tt, xx = TX
    rr = sp.Symbol("rr", positive=True)
    f_ex, vs_ex, xs_ex = 1 / (1 + rr**3), (1 + tt) / (1 + tt**2), tt**3 / 3 + tt / 2
    b_ex = (tt**2 + xx) * rr / (2 + rr**2 + tt * xx**2)  # a general b(t, x, rho), odd in rho, with rho -> rr
    print(f"  route E: f = {f_ex.subs(rr, sp.Symbol('r'))}, v_s = {vs_ex}, x_s = {xs_ex}; "
          f"b(t, x, rho) = {b_ex.subs(rr, RHO)}")
    times = [sp.Rational(1, 3), sp.Rational(-2, 5), sp.Rational(3, 2), sp.Rational(1, 7), sp.Rational(-5, 4), 2]
    bad_a = bad_g = 0
    for (x0, y0, z0, r0), tv in zip(pts[:6], times):
        Xv, yv, zv, Rv = [sp.Rational(c, 7) for c in (x0, y0, z0, r0)]
        sub = {X: Xv, y: yv, z: zv, R: Rv}
        sub.update({F[k]: sp.diff(f_ex, rr, k).subs(rr, Rv) for k in range(5)})
        sub.update({VS[k]: sp.diff(vs_ex, tt, k).subs(tt, tv) for k in range(4)})
        sub.update({XS[k]: sp.diff(xs_ex, tt, k).subs(tt, tv) for k in range(1, 5)})
        ex = ricci_at_point(-vs_ex * f_ex.subs(rr, sp.sqrt((xx - xs_ex)**2 + y**2 + z**2)),
                            (tv, xs_ex.subs(tt, tv) + Xv, yv, zv), lapse)
        bad_a += any(sp.simplify(e) != 0 for e in Ra.subs(sub) - ex)
        xv, rv = Xv + 1, sp.sqrt(yv**2 + zv**2)
        sub = {y: yv, z: zv, RHO: rv}
        for key in [k for k in BJ if isinstance(k, tuple)]:
            d = b_ex
            for c, m in zip((tt, xx, rr), key):
                d = sp.diff(d, c, m) if m else d
            sub[BJ[key]] = d.subs({tt: tv, xx: xv, rr: rv})
        ex = ricci_at_point(b_ex.subs(rr, sp.sqrt(y**2 + z**2)), (tv, xv, yv, zv), lapse)
        bad_g += any(sp.simplify(e) != 0 for e in Rg.subs(sub) - ex)
    check("gen_xcheck", bad_g == 0, f"route D for b(t, x, rho) differs from route E at {bad_g} of 6 exact points")
    check("accel_xcheck", bad_a == 0, f"route D for the accelerating bubble differs from route E at {bad_a} of 6 "
          "exact points")

    # (I1), (I2) and R_rho,rho = (d_rho b)^2 / 2 for a general b(t, x, rho), off the axis.
    rho2 = y**2 + z**2
    kw = 2 if mut == "accel-wall-v2" else 1
    rr2, rp2, pp2, ap = cylinder(Rg, tau)
    ids = [reduce_R(e, rho2, RHO) for e in ap + [rp2, pp2]]
    check("gen_I1", all(e == 0 for e in ids), f"b(t, x, rho): rho R(e_varphi, E_0), rho R(e_varphi, E_1), "
          f"rho^2 R_rho,varphi, rho^2 R_varphi,varphi reduce to {ids}")
    i2 = reduce_R(rho2 * (Rg[0, 0] + Rg[1, 1]) + 2 * rr2, rho2, RHO)
    check("gen_I2", i2 == 0, f"b(t, x, rho): rho^2 (R_00 + R_11 + 2 R_rho,rho) reduces to {sp.factor(i2)}")
    rem = reduce_R(rr2 - rho2 * (kw * bjet(0, 0, 1))**2 / 2, rho2, RHO)
    check("gen_Rrr", rem == 0, f"b(t, x, rho): rho^2 (R_rho,rho - (d_rho b)^2 / 2) reduces to {sp.factor(rem)}")

    # The same, the axis and the wall curvature for the accelerating bubble.
    rr2, rp2, pp2, ap = cylinder(Ra, tau)
    ids = [reduce_R(e) for e in ap + [rp2, pp2]]
    check("accel_I1", all(e == 0 for e in ids), f"bubble: rho R(e_varphi, E_0), rho R(e_varphi, E_1), "
          f"rho^2 R_rho,varphi, rho^2 R_varphi,varphi reduce to {ids}")
    i2 = reduce_R(rho2 * (Ra[0, 0] + Ra[1, 1]) + 2 * rr2)
    check("accel_I2", i2 == 0, f"bubble: rho^2 (R_00 + R_11 + 2 R_rho,rho) reduces to {sp.factor(i2)}")
    ax = Ra.subs({y: 0, z: 0})
    ids = [reduce_R(ax[a, b], X**2) for a in (2, 3) for b in range(4)] + [reduce_R(ax[0, 0] + ax[1, 1], X**2)]
    check("accel_axis", all(e == 0 for e in ids), f"bubble on the axis y = z = 0: R(E_2, .), R(E_3, .), "
          f"R_00 + R_11 reduce to {ids}")
    vw = kw * VS[0]
    rem1 = reduce_R(rr2 - rho2 * vw**2 * F[1]**2 * rho2 / (2 * R**2))
    rem2 = reduce_R(ax[0, 1] - vw * F[1] / R, X**2)
    check("accel_wall", rem1 == 0 and rem2 == 0,
          f"R_rho,rho = v_s^2 f'^2 rho^2 / (2 r_s^2) [residue {sp.factor(rem1)}] and, on the axis, R_01 = v_s f' / r_s "
          f"[residue {sp.factor(rem2)}]. So R_ab != 0 wherever v_s(t) f'(r_s) != 0.")

    failed = [n for n, ok in RESULTS if not ok]
    print(f"\nchecks: {len(RESULTS) - len(failed)} passed, {len(failed)} failed {failed}. {time.time() - t0:.0f} s")
    if mut is not None:
        want = MUTATIONS[mut][0]
        fired = [n for n in want if n in failed]
        print(f"control {mut}: target checks {want}, fired {fired}")
    return failed


if __name__ == "__main__":
    failed = main()
    controls_ok = True
    if "--mutate" not in " ".join(sys.argv) and "--no-controls" not in sys.argv:
        print("\n== Controls: each mutation must make its target checks fail ==", flush=True)
        for name, (want, what) in MUTATIONS.items():
            p = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True)
            fails = [ln.split(":")[0][7:] for ln in p.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = all(w in fails for w in want) and p.returncode == 1
            controls_ok &= fired
            print(f"[{'FIRED' if fired else 'SILENT'}] --mutate={name} ({what}): exit {p.returncode}, failed {fails}",
                  flush=True)
            if p.returncode not in (0, 1):
                print(p.stderr[-2000:])
        print(f"controls: {'all fired' if controls_ok else 'SOME DID NOT FIRE'}")
    sys.exit(0 if not failed and controls_ok else 1)

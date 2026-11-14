#!/usr/bin/env python3
"""Can a scalar source a warp shift that is not axisymmetric? (issue #70)

#57 and #61 showed that R_ab = lam g_ab + s u_a u_b, the Ricci tensor of a minimally coupled scalar (s = 8 pi eps,
lam = 8 pi V, eps = +-1), fails for every axisymmetric shift b(t, x, rho). Here the shift is

    ds^2 = -dt^2 + (dx + b dt)^2 + dy^2 + dz^2,   b = -v f(q),   q = (x - v t)^2 + y^2 + 2 z^2,

with f a function of q, so f(q) is not a function of r_s. The script reuses phantom_source.py's `scalar_form` (the same
function that passes the Ellis-Bronnikov control there) on the frame Ricci tensor at exact rational points, and its
Christoffel routines `ricci` and `ricci_at_point`.

    python3 notes/warp/nonaxisymmetric.py                 all checks, then every control (exit 1 on any failure)
    python3 notes/warp/nonaxisymmetric.py --mutate=NAME   one control (its target checks must fail, exit 1)
    python3 notes/warp/nonaxisymmetric.py --no-controls   the checks only

Routes for the frame Ricci tensor (E_0 = n = d_t - b d_x, E_i = d_i). (D) phantom_source's `ricci` in the lab
coordinates, with f(q) differentiated on its jets G_k = f^(k)(q) as free symbols. (E) phantom_source's
`ricci_at_point`, the Christoffel formula in exact rationals at a point, with sympy's own diff of an explicit b. (T) The
3+1 equations (momentum constraint, Hamiltonian constraint and the evolution equation for K_ij), also at a point, with
sympy's own diff of the same explicit b. D and E share `shift_metric`, E and T share the explicit b, and T shares
nothing else with D or E.

The lemma of phantom_source.md holds in any orthonormal frame, so it is tested in a frame adapted to the shift, e_2
along grad_perp q = (q_y, q_z) and e_3 along (-q_z, q_y), in place of the circular (e_rho, e_varphi). Its premises (I1)
and (I2) are polynomial identities in the jets, with no relation to reduce modulo.

`--mutate=shear` adds X y / 2 to q. That shift is not of the form b(t, x, sigma(y, z)), the premises fail for it, and
the controls section echoes the scalar-form outcome there.
"""
import random
import subprocess
import sys
import time

import mpmath as mp
import sympy as sp

import phantom_source as ps

X, y, z, v = ps.X, ps.y, ps.z, ps.v
TX = ps.TX  # lab t and x, the symbols of the explicit functions of routes E and T
G = sp.symbols("G0:4", real=True)  # G[k] = f^(k)(q)
ETA = sp.diag(-1, 1, 1, 1)

# name: (checks that must fail, what the mutation does)
MUTATIONS = {
    "static": (["xcheck_E", "xcheck_T", "sigma_xcheck"], "drops every t-derivative in route D's jets"),
    "lapse": (["xcheck_T", "identity_I1", "identity_I2", "identity_axis", "wall_curvature", "lemma_premises",
               "sigma_I1", "sigma_I2", "sigma_Rperp"], "gives the slicing the lapse 1 + b in routes D and E"),
    "circular": (["identity_I1", "identity_I2", "lemma_premises", "sigma_I1", "sigma_I2"],
                 "tests (I1), (I2) in the circular frame (e_rho, e_varphi), not the frame adapted to the shift"),
    "shear": (["identity_I1", "identity_axis", "lemma_premises"],
              "adds X y / 2 to q, so grad_perp d_x b is no longer along grad_perp b"),
    "wall-v2": (["wall_curvature"], "doubles v in the expected wall curvature"),
    "planted-sign": (["planted_pass"], "plants lam eta + 8 pi u u, the wrong sign for eps = -1"),
    "planted-null": (["planted_deform"], "deforms the planted form by zero"),
    "planted-flat": (["planted_premise"], "plants R = 0, which satisfies the lemma's premises"),
    "broken-test": (["points_fail", "lemma_agree", "planted_sign", "planted_deform"],
                    "replaces scalar_form by a test that always returns 'pass'"),
    "planted-in-data": (["points_fail", "lemma_premises"],
                        "replaces the Ricci tensor at the sample points by a planted scalar form"),
}

ECHO = {"shear": ["points_fail", "lemma_premises"]}  # checks whose line of a mutation's run the controls section echoes

RESULTS = []


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def qexpr(mut):
    return X**2 + y**2 + 2 * z**2 + (X * y / 2 if mut == "shear" else 0)


def lemma_residues(Rt, w):
    """The lemma's premises in the frame (E_0, E_1, w, wp), wp = (-w_1, w_0), both unnormalised. Returns the residues
    of (I1), R(wp, E_0), R(wp, E_1), R(wp, w), R(wp, wp), and of (I2), |w|^2 (R_00 + R_11) + 2 R(w, w)."""
    wv = sp.Matrix(w)
    wp = sp.Matrix([-w[1], w[0]])
    T = Rt[2:, 2:]
    I1 = [sp.expand((wp.T * Rt[2:, a])[0]) for a in (0, 1)] + [sp.expand((wp.T * T * wv)[0]),
                                                                  sp.expand((wp.T * T * wp)[0])]
    I2 = sp.expand((w[0]**2 + w[1]**2) * (Rt[0, 0] + Rt[1, 1]) + 2 * (wv.T * T * wv)[0])
    return I1, I2


# ---------------------------------------------------------------------------------------------------------------
# Route D: the frame Ricci tensor with f(q) on its jets
# ---------------------------------------------------------------------------------------------------------------
def q_rule(q, static):
    """Jets of b = -v G_0, G_k = f^(k)(q), q(X, y, z), X = x - v t. d_t = -v d_X on functions of (X, y, z)."""
    qX, qy, qz = (sp.diff(q, s) for s in (X, y, z))
    dq = [0 if static else -v * qX, qX, qy, qz]

    def rule(s):
        if s == X:
            return [0 if static else -v, 1, 0, 0]
        if s in G:
            k = G.index(s)
            assert k < 3, "the Ricci tensor has no third jet of f"
            return [G[k + 1] * d for d in dq]
        return None
    return rule


def route_D(mut):
    q = qexpr(mut)
    return ps.shift_tetrad_ricci(-v * G[0], q_rule(q, mut == "static"), mut == "lapse")


# ---------------------------------------------------------------------------------------------------------------
# Routes E and T: explicit b at an exact point. f is the quadratic in q with the point's jets.
# ---------------------------------------------------------------------------------------------------------------
def explicit_b(vv, G3, q0, mut):
    """b(t, x, y, z) = -v f(q), f(q) = G0 + G1 (q - q0) + G2 (q - q0)^2 / 2, so f has the jets G at q0."""
    qq = qexpr(mut).subs(X, TX[1] - vv * TX[0])
    return -vv * (G3[0] + G3[1] * (qq - q0) + G3[2] * (qq - q0)**2 / 2)


def route_T(b, p):
    """Frame Ricci (E_0 = n, E_i = d_i) from the 3+1 equations, lapse 1, flat slices, beta^i = (b, 0, 0):
    K_ij = (d_i b_j + d_j b_i) / 2, G_nn = (K^2 - K_ij K^ij) / 2, R_ni = d_i K - d_j K_ji,
    R_ij = K K_ij - 2 K_ik K_kj - (d_t - L_beta) K_ij, R_nn = 2 G_nn - R_ii."""
    cs = (TX[0], TX[1], y, z)
    at = dict(zip(cs, p))
    sp_ = cs[1:]  # x, y, z

    def val(e):
        r = sp.cancel(e.subs(at))
        assert r.is_Rational, r
        return r
    db = [sp.diff(b, s) for s in sp_]
    K = sp.Matrix(3, 3, lambda i, j: (db[j] if i == 0 else 0) / 2 + (db[i] if j == 0 else 0) / 2)
    trK = K.trace()
    Gnn = (trK**2 - sum(K[i, j]**2 for i in range(3) for j in range(3))) / 2
    Rni = [sp.diff(trK, sp_[i]) - sum(sp.diff(K[j, i], sp_[j]) for j in range(3)) for i in range(3)]

    def Rij(i, j):
        LK = b * sp.diff(K[i, j], sp_[0]) + K[0, j] * db[i] + K[i, 0] * db[j]
        return trK * K[i, j] - 2 * sum(K[i, k] * K[k, j] for k in range(3)) - sp.diff(K[i, j], cs[0]) + LK
    Rs = sp.Matrix(3, 3, lambda i, j: Rij(i, j))
    Rnn = 2 * Gnn - Rs.trace()
    out = sp.zeros(4, 4)
    out[0, 0] = val(Rnn)
    for i in range(3):
        out[0, i + 1] = out[i + 1, 0] = val(Rni[i])
        for j in range(3):
            out[i + 1, j + 1] = val(Rs[i, j])
    return out


# ---------------------------------------------------------------------------------------------------------------
# A general shift b(t, x, sigma(y, z)): jets B_ijk = d_t^i d_x^j d_sigma^k B and sigma_y, sigma_z, ... (route D again)
# ---------------------------------------------------------------------------------------------------------------
SJ = sp.symbols("sy sz syy syz szz", real=True)
BJ = {}


def bjet(i, j, k):
    if (i, j, k) not in BJ:
        BJ[(i, j, k)] = s = sp.Symbol(f"B{i}{j}{k}", real=True)
        BJ[s] = (i, j, k)
    return BJ[(i, j, k)]


def sigma_rule(static):
    sy, sz, syy, syz, szz = SJ

    def rule(s):
        if s == sy:
            return [0, 0, syy, syz]
        if s == sz:
            return [0, 0, syz, szz]
        if s in SJ:
            raise ValueError("the Ricci tensor has no third jet of sigma")
        if s in BJ:
            i, j, k = BJ[s]
            up = bjet(i, j, k + 1)
            return [0 if static else bjet(i + 1, j, k), bjet(i, j + 1, k), up * sy, up * sz]
        return None
    return rule


# ---------------------------------------------------------------------------------------------------------------
# Profiles (jets in q) and sample points
# ---------------------------------------------------------------------------------------------------------------
def tanh_jets(q0, sig=8, Rb=1):
    """Jets in q of Alcubierre's f = [tanh(sig (r + Rb)) - tanh(sig (r - Rb))] / (2 tanh(sig Rb)), r = q^(1/2)."""
    mp.mp.dps = 50
    fn = lambda qq: (mp.tanh(sig * (mp.sqrt(qq) + Rb)) - mp.tanh(sig * (mp.sqrt(qq) - Rb))) / (2 * mp.tanh(sig * Rb))
    return [sp.Rational(mp.nstr(mp.diff(fn, mp.mpf(q0.p) / q0.q, k), 45)) for k in range(3)]


def bump_jets(q0, Q=4):
    """Jets in q of the C^2 bump (1 - q / Q)^3 for q < Q, 0 beyond. Exact."""
    if q0 >= Q:
        return [sp.Integer(0)] * 3
    s = 1 - q0 / sp.Integer(Q)
    return [s**3, -3 * s**2 / Q, 6 * s / Q**2]


def rational_jets(q0):
    """Jets in q of f = 1 / (1 + q^2). Exact."""
    return [1 / (1 + q0**2), -2 * q0 / (1 + q0**2)**2, (6 * q0**2 - 2) / (1 + q0**2)**3]


def rat(a, b=1):
    return sp.Rational(a, b)


# (label, (X, y, z), v). Points on q = 1 (inside Alcubierre's wall for sigma = 8, R = 1, and inside the bump).
ON_Q1 = [("axis", (1, 0, 0)), ("axis", (-1, 0, 0)), ("y-axis", (0, 1, 0)), ("z = 0 plane", (rat(3, 5), rat(4, 5), 0)),
         ("y = 0 plane", (rat(1, 3), 0, rat(2, 3))), ("X = 0 plane", (0, rat(1, 3), rat(2, 3))),
         ("general", (rat(1, 2), rat(1, 2), rat(1, 2))), ("general", (rat(4, 5), rat(1, 5), rat(2, 5))),
         ("general", (rat(-1, 5), rat(4, 5), rat(2, 5))), ("general", (rat(5, 6), rat(1, 2), rat(1, 6))),
         ("general", (rat(-4, 7), rat(1, 7), rat(4, 7)))]


def sample_points(mut):
    """The exact points of the scalar-form test, as (profile, kind, (X, y, z), v, jets in q)."""
    rng = random.Random(70)
    q = qexpr(mut)
    pts = []

    def add(prof, kind, p, vv, jetfn):
        q0 = sp.nsimplify(q.subs(dict(zip((X, y, z), p))))
        pts.append((prof, kind, tuple(sp.Rational(c) for c in p), sp.Rational(vv), jetfn(q0)))
    for i, (kind, p) in enumerate(ON_Q1):
        add("tanh, sigma = 8, R = 1", kind, p, [1, 1, 1, rat(1, 2), 3, 1, 3, 1, rat(1, 2), 1, 3][i], tanh_jets)
    for kind, p in [ON_Q1[0], ON_Q1[3], ON_Q1[4], ON_Q1[6], ON_Q1[8]]:
        add("bump (1 - q / 4)^3", kind, p, rat(7, 4), bump_jets)
    add("bump (1 - q / 4)^3", "general", (rat(3, 2), rat(1, 3), rat(-1, 2)), rat(7, 4), bump_jets)
    for kind, p in [ON_Q1[1], ON_Q1[5], ON_Q1[7]]:
        add("1 / (1 + q^2)", kind, p, 2, rational_jets)
    add("1 / (1 + q^2)", "general", (rat(1, 3), rat(-2, 5), rat(3, 4)), rat(3, 2), rational_jets)
    def random_jets(q0):  # f' and f'' nonzero, f arbitrary, so no profile in particular
        nz = lambda: rat(rng.choice([-1, 1]) * rng.randint(1, 9), rng.randint(1, 9))
        return [rat(rng.randint(-9, 9), rng.randint(1, 9)), nz(), nz()]
    for kind, p in [ON_Q1[0], ON_Q1[3], ON_Q1[6], ON_Q1[9]]:
        add("random jets", kind, p, rat(rng.randint(1, 9), 3), random_jets)
    add("random jets, f' = 0", "general", ON_Q1[7][1], rat(5, 2), lambda q0: [rat(2, 3), 0, rat(-3, 5)])
    add("bump, beyond it (f = 0)", "general", (2, 1, rat(1, 2)), rat(7, 4), bump_jets)
    return pts


# ---------------------------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------------------------
def fmt(r):
    return str(r) if sp.Rational(r).q < 10**6 else f"{float(r):.6g}"


def main():
    opts = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
    mut = opts.get("mutate")
    if mut is not None and mut not in MUTATIONS:
        sys.exit(f"unknown mutation {mut!r}. Known: {', '.join(MUTATIONS)}")
    random.seed(70)
    t0 = time.time()
    print(f"nonaxisymmetric.py, issue #70. sympy {sp.__version__}, mpmath {mp.__version__}, python "
          f"{sys.version.split()[0]}. mutation: {mut or 'none'}" + (f" ({MUTATIONS[mut][1]})" if mut else ""))
    lapse = mut == "lapse"
    q = qexpr(mut)
    qy, qz = sp.diff(q, y), sp.diff(q, z)
    wfr = (y, z) if mut == "circular" else (qy, qz)  # the frame's e_2 direction, in the (y, z) plane

    # 1. Route D, and its cross-checks against routes E and T, which share none of its machinery.
    print(f"\n== 1. Frame Ricci tensor of b = -v f(q), q = {q}, and two independent routes ==")
    Rt = route_D(mut)
    names = ["0", "1", "2", "3"]
    for a in range(4):
        for b in range(a, 4):
            print(f"  R_({names[a]},{names[b]}) = {sp.factor(Rt[a, b])}")
    print("  (frame: E_0 = n, E_1 = d_x, E_2 = d_y, E_3 = d_z. G_k = f^(k)(q).)")
    rng = random.Random(7)
    badE = badT = nx_done = 0
    nx = 8
    while nx_done < nx:
        p = [rat(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(3)]
        vv, tv = rat(rng.randint(1, 9), rng.randint(1, 5)), rat(rng.randint(-5, 5), rng.randint(1, 4))
        jets = [rat(rng.randint(-9, 9), rng.randint(1, 9)) for _ in range(3)]
        if vv * jets[0] == 1:  # the lapse 1 + b of the `lapse` control would vanish
            continue
        nx_done += 1
        q0 = q.subs(dict(zip((X, y, z), p)))
        bex = explicit_b(vv, jets, q0, mut)
        pt = (tv, p[0] + vv * tv, p[1], p[2])
        D_ = Rt.subs({X: p[0], y: p[1], z: p[2], v: vv, G[0]: jets[0], G[1]: jets[1], G[2]: jets[2]})
        E_ = ps.ricci_at_point(bex, pt, lapse)
        T_ = route_T(bex, pt)
        badE += any(sp.simplify(e) != 0 for e in D_ - E_)
        badT += any(sp.simplify(e) != 0 for e in D_ - T_)
    check("xcheck_E", badE == 0, f"route D differs from route E (Christoffel at a point) at {badE} of {nx} exact "
          "points (random v, point and jets)")
    check("xcheck_T", badT == 0, f"route D differs from route T (3+1 equations) at {badT} of {nx} exact points")

    # 2. The identities of the lemma in the frame adapted to the shift, for general jets.
    print("\n== 2. (I1) R(e_3, .) = 0 and (I2) R_00 + R_11 + 2 R_22 = 0, e_2 along grad_perp q, e_3 normal to it ==")
    ax = Rt.subs({y: 0, z: 0})
    I1, I2 = lemma_residues(Rt, wfr)
    check("identity_I1", all(e == 0 for e in I1),
          f"R(wp, E_0), R(wp, E_1), R(wp, w), R(wp, wp) with w = {tuple(wfr)}, wp = (-w_1, w_0), general jets: {I1}")
    check("identity_I2", I2 == 0, f"|w|^2 (R_00 + R_11) + 2 R(w, w) = {sp.factor(I2)}")
    ids = [sp.expand(ax[a, b]) for a in (2, 3) for b in range(4)] + [sp.expand(ax[0, 0] + ax[1, 1])]
    check("identity_axis", all(e == 0 for e in ids), f"on the axis y = z = 0, where w = 0: R(E_2, .), R(E_3, .), "
          f"R_00 + R_11 = {ids}")
    kw = 2 if mut == "wall-v2" else 1
    w = sp.Matrix((qy, qz))
    Rww = sp.expand((w.T * Rt[2:, 2:] * w)[0])
    Rww_exp = sp.expand(sp.Rational(1, 2) * (kw * v * G[1])**2 * (qy**2 + qz**2)**2)
    R01 = sp.expand(Rt[0, 1])
    R01_exp = sp.expand(kw * v * (G[1] * (sp.diff(q, y, 2) + sp.diff(q, z, 2)) + G[2] * (qy**2 + qz**2)) / 2)
    check("wall_curvature", sp.expand(Rww - Rww_exp) == 0 and sp.expand(R01 - R01_exp) == 0,
          f"R(w, w) - v^2 f'^2 |w|^4 / 2 = {sp.expand(Rww - Rww_exp)} and R_01 - v (f' lap_perp q + f'' |w|^2) / 2 = "
          f"{sp.expand(R01 - R01_exp)}, with w = grad_perp q")
    cc = lemma_residues(Rt, (y, z))
    print(f"  info: in the circular frame w = (y, z) the same residues are R(wp, .) = {[sp.factor(e) for e in cc[0]]}"
          f" and (I2) = {sp.factor(cc[1])}, so (I1) and (I2) as stated for (e_rho, e_varphi) fail for this shift.")
    print("  Lemma (by hand, phantom_source.md): in any orthonormal frame, R(e_3, .) = 0 and R_00 + R_11 + 2 R_22 = 0\n"
          "  force R_ab = lam eta_ab + s u_a u_b to have R_ab = 0. Here R_ab != 0 wherever v f'(q) != 0, since off\n"
          f"  the axis R(w, w) = {sp.factor(Rww)} and on the axis R_01 = {sp.factor(ax[0, 1])}.")

    # 3. A general shift b(t, x, sigma(y, z)), of which the shift above is the case sigma = y^2 + 2 z^2.
    print("\n== 3. The same identities for any b(t, x, sigma(y, z)), in the frame adapted to grad sigma ==")
    static = mut == "static"
    bs = bjet(0, 0, 0)
    Rs = ps.shift_tetrad_ricci(bs, sigma_rule(static), lapse)
    sy, sz = SJ[0], SJ[1]
    tt, xx, ss = sp.symbols("tt xx ss", real=True)
    B_ex = (tt**2 + xx) * ss / (2 + ss**2 + tt * xx**2)  # a general B(t, x, sigma)
    sig_ex = y**2 + 2 * z**2 + y * z / 2 + y / 3  # a general sigma(y, z), not even in y
    bad = 0
    tvs = [rat(1, 3), rat(-2, 5), rat(3, 2), rat(1, 7), rat(-5, 4), 2]
    spts = [(1, rat(1, 2), rat(1, 3)), (rat(-1, 2), 1, rat(-1, 2)), (rat(1, 4), rat(1, 4), 1),
            (2, rat(-1, 3), rat(1, 5)), (0, rat(2, 3), rat(-1, 2)), (rat(-3, 2), rat(1, 6), rat(1, 4))]
    for tv, (x0, y0, z0) in zip(tvs, spts):
        at = {y: y0, z: z0}
        sv = sig_ex.subs(at)
        sub = {y: y0, z: z0, SJ[0]: sp.diff(sig_ex, y).subs(at), SJ[1]: sp.diff(sig_ex, z).subs(at),
               SJ[2]: sp.diff(sig_ex, y, 2), SJ[3]: sp.diff(sig_ex, y, z), SJ[4]: sp.diff(sig_ex, z, 2)}
        for key in [k for k in BJ if isinstance(k, tuple)]:
            d = B_ex
            for c, m in zip((tt, xx, ss), key):
                d = sp.diff(d, c, m) if m else d
            sub[BJ[key]] = d.subs({tt: tv, xx: x0, ss: sv})
        bex = B_ex.subs({tt: TX[0], xx: TX[1], ss: sig_ex})
        E_ = ps.ricci_at_point(bex, (tv, x0, y0, z0), lapse)
        bad += any(sp.simplify(e) != 0 for e in Rs.subs(sub) - E_)
    check("sigma_xcheck", bad == 0, f"route D for b(t, x, sigma(y, z)) differs from route E at {bad} of 6 exact points")
    ws = (y, z) if mut == "circular" else (sy, sz)
    I1s, I2s = lemma_residues(Rs, ws)
    check("sigma_I1", all(e == 0 for e in I1s), f"R(wp, .) with w = grad sigma = {tuple(ws)}: {I1s}")
    check("sigma_I2", I2s == 0, f"|w|^2 (R_00 + R_11) + 2 R(w, w) = {sp.factor(I2s)}")
    Rw = sp.expand((sp.Matrix(ws).T * Rs[2:, 2:] * sp.Matrix(ws))[0])
    Rw_exp = sp.expand(bjet(0, 0, 1)**2 * (sy**2 + sz**2)**2 / 2)
    check("sigma_Rperp", sp.expand(Rw - Rw_exp) == 0, f"R(w, w) - (d_sigma B)^2 |w|^4 / 2 = {sp.expand(Rw - Rw_exp)}")

    # 4. Planted truth, before real data: lam eta - 8 pi u u at the same points, through the same scalar_form.
    print("\n== 4. Planted control at the sample points: R = lam eta - 8 pi u u passes, a deformed one fails ==")
    pts = sample_points(mut)

    def premises(Rm, p):
        """The lemma's premises at an exact point: (I1) and (I2) in the adapted frame, from Rm."""
        w0 = (qy.subs(dict(zip((X, y, z), p))), qz.subs(dict(zip((X, y, z), p))))
        if mut == "circular":
            w0 = (p[1], p[2])
        if w0 == (0, 0):
            w0 = (1, 0)  # on the axis the frame is (E_2, E_3), and R(E_2, .) = 0 there too
        i1, i2 = lemma_residues(Rm, w0)
        return all(e == 0 for e in i1) and i2 == 0

    def test(Rm, eps):
        """The test, as one function for the planted and the real data."""
        if mut == "broken-test":
            return "pass", 0, sp.zeros(4, 4), sp.zeros(4, 1)
        return ps.scalar_form(Rm, ETA, eps)

    planted = []
    for (prof, kind, p, vv, jets) in pts:
        lam_p = rat(random.randint(-9, 9), random.randint(1, 9))
        u_p = sp.Matrix([rat(random.choice([-1, 1]) * random.randint(1, 9), random.randint(1, 9)) for _ in range(4)])
        if mut == "planted-flat":
            lam_p, u_p = sp.Integer(0), sp.zeros(4, 1)
        sgn = 1 if mut == "planted-sign" else -1
        planted.append((lam_p, u_p, lam_p * ETA + sgn * 8 * sp.pi * u_p * u_p.T))
    bad_pass = bad_sign = bad_def = bad_prem = 0
    for (prof, kind, p, vv, jets), (lam_p, u_p, Rp) in zip(pts, planted):
        st, lr, M, u = test(Rp, -1)
        ok = st == "pass" and lr == lam_p and (M + 8 * sp.pi * u_p * u_p.T).applyfunc(sp.simplify).is_zero_matrix \
            and (u * u.T - u_p * u_p.T).applyfunc(sp.simplify).is_zero_matrix
        bad_pass += not ok
        bad_sign += test(Rp, +1)[0] != "sign"
        mu = 0 if mut == "planted-null" else rat(1, 10)
        Rd = Rp + mu * sp.Matrix(4, 4, lambda i, j: 1 if (i, j) == (3, 3) else 0)
        bad_def += test(Rd, -1)[0] != "algebraic"
        bad_prem += premises(Rp, p)  # a planted scalar form with R != 0 must violate the lemma's premises
    n = len(pts)
    check("planted_pass", bad_pass == 0, f"lam eta - 8 pi u u, eps = -1: status 'pass' with lam and u u recovered "
          f"exactly at {n - bad_pass} of {n} points")
    check("planted_sign", bad_sign == 0, f"the same form with eps = +1: status 'sign' at {n - bad_sign} of {n} points")
    check("planted_deform", bad_def == 0, f"the form plus mu E_3 E_3 (rank 2): status 'algebraic' at {n - bad_def} of "
          f"{n} points")
    check("planted_premise", bad_prem == 0, f"the planted form violates (I1) or (I2) of the adapted frame at "
          f"{n - bad_prem} of {n} points (the lemma's contrapositive)")

    # 5. The scalar-form test on the real Ricci tensor, point by point.
    print("\n== 5. Scalar-form test at exact points, eps = -1 and +1 (phantom_source.scalar_form) ==")
    print(f"  {'#':>2}  {'profile':<24} {'kind':<12} {'v':>5} {'X':>6} {'y':>6} {'z':>6} {'q':>5} {'f_q':>11} "
          f"{'R!=0':>5} {'(I1),(I2)':>10}  eps=-1     eps=+1")
    outcomes = []
    n_fail = n_bad_prem = n_bad_agree = 0
    for i, (prof, kind, p, vv, jets) in enumerate(pts, 1):
        sub = {X: p[0], y: p[1], z: p[2], v: vv, G[0]: jets[0], G[1]: jets[1], G[2]: jets[2]}
        Rm = Rt.subs(sub)
        if mut == "planted-in-data":
            Rm = planted[i - 1][2]
        nz = not Rm.is_zero_matrix
        prem = premises(Rm, p)
        sts = [test(Rm, e)[0] for e in (-1, 1)]
        want = "trivial" if not nz else "algebraic"
        n_fail += any(s != want for s in sts)
        n_bad_prem += not prem
        n_bad_agree += any(prem and s != want for s in sts)
        outcomes.append(sts)
        q0 = sp.nsimplify(q.subs(dict(zip((X, y, z), p))))
        print(f"  {i:>2}  {prof:<24} {kind:<12} {fmt(vv):>5} {fmt(p[0]):>6} {fmt(p[1]):>6} {fmt(p[2]):>6} {fmt(q0):>5} "
              f"{float(jets[1]):>11.4g} {str(nz):>5} {str(prem):>10}  {sts[0]:<10} {sts[1]}")
    nw = sum(1 for (_, _, _, _, j), s in zip(pts, outcomes) if j[1] != 0 and s == ["algebraic"] * 2)
    check("points_fail", n_fail == 0, f"{n} points, both signs ({2 * n} tests). R != 0 gives 'algebraic' (no lam makes "
          f"R - lam eta of rank <= 1) and R = 0 gives 'trivial'. {nw} points with f'(q) != 0 fail algebraically at "
          f"both signs, {n_fail} points differ from that")
    check("lemma_premises", n_bad_prem == 0,
          f"(I1) and (I2) of the adapted frame hold at {n - n_bad_prem} of {n} points")
    check("lemma_agree", n_bad_agree == 0, f"where the premises hold the lemma's prediction ('algebraic' if R != 0, "
          f"'trivial' if R = 0) agrees with scalar_form at {n - n_bad_agree} of {n} points, both signs")

    failed = [nm for nm, ok in RESULTS if not ok]
    print(f"\nchecks: {len(RESULTS) - len(failed)} passed, {len(failed)} failed {failed}. {time.time() - t0:.0f} s")
    if mut is not None:
        want = MUTATIONS[mut][0]
        fired = [nm for nm in want if nm in failed]
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
            for ln in p.stdout.splitlines():
                if any(ln.startswith(f"[PASS] {c}:") or ln.startswith(f"[FAIL] {c}:") for c in ECHO.get(name, [])):
                    print(f"  info: under --mutate={name}, {ln[7:]}", flush=True)
            if p.returncode not in (0, 1):
                print(p.stderr[-2000:])
        print(f"controls: {'all fired' if controls_ok else 'SOME DID NOT FIRE'}")
    sys.exit(0 if not failed and controls_ok else 1)

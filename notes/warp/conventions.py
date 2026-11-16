#!/usr/bin/env python3
"""The spinor conventions of mass_bounds.md, checked symbolically, for issues #62, #67 and #86.

On the flat slice of ds^2 = -dt^2 + sum_i (dx^i + b^i dt)^2 (lapse 1, any shift b^i(t, x, y, z)), the note's
nabla~_i = d_i + (1/2) K_ij g^j g^0, with the note's gammas g^a, equals the restriction to the slice of the spacetime
spin connection (1/4) omega_ab G^a G^b in the gammas G^0 = -g^0, G^j = g^j. The other sign of the K term does not (#62).
Section 4 repeats this for ds^2 = -N^2 dt^2 + psi^4 sum_i (dx^i + b^i dt)^2, with N, psi and b^i all arbitrary
functions of (t, x, y, z), where nabla~_i = nabla_i + (1/2) K_ij g^j g^0 and nabla_i is the spin connection of the
slice itself (#67). Planted metrics with a known spin connection pass first.
Sections 6 to 8 repeat this for a diagonal slice g_ij = diag(A^2, B^2, C^2)(t, x, y, z) with e_i = d_i / A_i, and for the
same slice in the orthonormal frame rotated by an angle alpha(t, x, y, z) about z (#86). Planted Bianchi I, Kasner and
Cartan-formula metrics pass first.

The Levi-Civita connection comes from the Christoffel symbols of the coordinate metric, with
omega_ab(X) = g(e_a, nabla_X e_b), e_0 = n and e_i = psi^-2 d_i. It shares no code with review_cartan.py, which solves
Cartan's structure equations instead, or with alcubierre_bounds.py and review_spinor.py. It shares sympy 1.14.0 and the
gammas of mass_bounds.md, "Conventions", which are typed in again here. On the general slice, the slice's own spin
connection comes from the 3-metric psi^4 delta alone, and K_ij from the ADM formula
(1/2N)(-d_t g_ij + D_i b_j + D_j b_i) with the 3-metric's Christoffel symbols. Neither uses the 4-metric. Both routes
share the helper that turns a metric into Christoffel symbols and the helper that differentiates a vector field. In the
rotated frame both routes rotate the same frame by the same cos and sin, and the identity cos^2 + sin^2 = 1 is the only
relation used beyond expanding Laurent polynomials in the jets of the unknown functions.

    python3 notes/warp/conventions.py                  checks, then every control
    python3 notes/warp/conventions.py --no-controls    checks only
    python3 notes/warp/conventions.py --mutate=NAME    one control (its target checks must fail, exit 1)

Exit status is 1 if a check fails or a control does not fire. A minute or two, one process.
"""
import subprocess
import sys
import time

import sympy as sp

OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate", "")
CONTROLS = {
    "sign": ("writes the stated nabla~ with the other sign, d_i - (1/2) K_ij g^j g^0",
             ["restriction", "dirac", "restriction_general", "dirac_general", "restriction_diag", "dirac_diag",
              "restriction_rot", "dirac_rot"]),
    "christoffel": ("flips the sign of the Christoffel term in nabla_X Y",
                    ["planted_desitter", "omega_slice", "planted_ds_rescaled", "planted_conf", "planted_static",
                     "omega_general", "restriction_general", "planted_bianchi", "planted_kasner",
                     "planted_diag_cartan", "omega_diag", "restriction_diag", "omega_rot", "restriction_rot"]),
    "spin-sign": ("uses -(1/4) omega_ab G^a G^b for the spin connection",
                  ["clifford_parallel", "restriction", "note_gammas", "clifford_parallel_general",
                   "restriction_general", "note_gammas_general", "clifford_parallel_diag", "restriction_diag",
                   "note_gammas_diag", "clifford_parallel_rot", "restriction_rot", "note_gammas_rot"]),
    "frame": ("takes e_0 = d_t + b^i d_i, the wrong sign of the shift",
              ["orthonormal", "orthonormal_general", "restriction_general", "orthonormal_diag", "restriction_diag",
               "orthonormal_rot", "restriction_rot"]),
    "gamma": ("scales g^1 by 11/10", ["clifford", "causal", "clifford_parallel_general", "dirac_general",
                                      "clifford_parallel_diag", "dirac_diag", "clifford_parallel_rot", "dirac_rot"]),
    "adm-sign": ("takes the ADM formula with the other sign, K_ij = (1/2N)(d_t g_ij - D_i b_j - D_j b_i)",
                 ["planted_ds_rescaled", "planted_conf", "planted_flrw", "omega_general", "restriction_general",
                  "note_gammas_general", "planted_bianchi", "planted_kasner", "planted_diag_cartan", "omega_diag",
                  "restriction_diag", "note_gammas_diag", "omega_rot", "restriction_rot", "note_gammas_rot"]),
    "lapse-K": ("drops the 1/N of the ADM formula for K_ij",
                ["planted_ds_rescaled", "planted_conf", "planted_flrw", "omega_general", "restriction_general",
                 "planted_bianchi", "planted_diag_cartan", "omega_diag", "restriction_diag", "omega_rot",
                 "restriction_rot"]),
    "lapse-frame": ("takes e_0 = d_t - b^i d_i without the 1/N",
                    ["orthonormal_general", "restriction_general", "planted_bianchi", "planted_diag_cartan",
                     "orthonormal_diag", "restriction_diag", "orthonormal_rot", "restriction_rot"]),
    "psi-frame": ("takes e_i = psi^-1 d_i instead of psi^-2 d_i", ["orthonormal_general", "restriction_general"]),
    "psi-slice": ("gives the slice's own spin connection the metric A_i instead of A_i^2 (psi^2 for psi^4)",
                  ["omega_general", "restriction_general", "omega_diag", "restriction_diag", "omega_rot",
                   "restriction_rot"]),
    "aniso-frame": ("takes e_1 = d_1 / A_2, the scale of the wrong axis",
                    ["planted_bianchi", "planted_kasner", "planted_diag_cartan", "orthonormal_diag",
                     "restriction_diag", "orthonormal_rot", "restriction_rot"]),
    "kasner": ("takes the exponents (1/2, 1/2, 1/2), which are not vacuum, for the Kasner metric",
               ["planted_kasner"]),
    "rot-ortho": ("rotates e_1 only and leaves e_2, so the frame is not orthonormal", ["orthonormal_rot"]),
    "rot-k": ("takes K_ij in the unrotated frame on the slice route",
              ["omega_rot", "restriction_rot", "note_gammas_rot"]),
    "rot-dalpha": ("drops the derivatives of alpha in the slice's own spin connection",
                   ["omega_rot", "restriction_rot", "note_gammas_rot"]),
}
if MUT and MUT not in CONTROLS:
    sys.exit(f"unknown mutation {MUT}")
RESULTS = []

T, X, Y, Z = sp.symbols("t x y z", real=True)
CRD = (T, X, Y, Z)
ETA = sp.diag(-1, 1, 1, 1)


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


ALPHA = sp.Function("alpha")(*CRD)  # the angle of the rotated frame, section 8


def trig_reduce(e):
    """Reduces the expanded e modulo cos(alpha)^2 + sin(alpha)^2 = 1. In c = cos(alpha) and s = sin(alpha) the normal
    form has degree at most 1 in c, so it is zero exactly when e is in the ideal."""
    c, s = sp.Dummy("c"), sp.Dummy("s")
    e = e.xreplace({sp.cos(ALPHA): c, sp.sin(ALPHA): s})
    out = 0
    for key, co in sp.collect(e, c, evaluate=False).items():
        k = 0 if key == 1 else key.as_base_exp()[1]
        out += co * c ** (k % 2) * (1 - s**2) ** (k // 2)
    return sp.expand(out)


def zero(e):
    """True if e vanishes. Every expression here is a Laurent polynomial in the jets of the unknown functions, and in
    cos(alpha) and sin(alpha), for which expand and trig_reduce are canonical, so simplify only runs on a residue that
    they did not clear."""
    e = sp.expand(e)
    if e != 0 and e.has(sp.cos(ALPHA), sp.sin(ALPHA)):
        e = trig_reduce(e)
    return e == 0 or sp.simplify(e) == 0


def is_zero(M):
    M = sp.Matrix(M)
    return all(zero(e) for e in M)


# ---- Gammas of mass_bounds.md, "Conventions" ------------------------------------------------------------------------
def note_gammas():
    s = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
    z2 = sp.zeros(2)
    g0 = sp.I * sp.diag(1, 1, -1, -1)
    gj = [sp.Matrix(sp.BlockMatrix([[z2, sj], [sj, z2]])) for sj in s]
    if MUT == "gamma":
        gj[0] = gj[0] * sp.Rational(11, 10)
    return g0, gj


def section_gammas():
    print("\n== 1. The gammas ==")
    g0, gj = note_gammas()
    G = [g0] + gj
    anti = all(is_zero(G[a] * G[b] + G[b] * G[a] - 2 * ETA[a, b] * sp.eye(4)) for a in range(4) for b in range(4))
    Gs = [-g0] + gj  # spacetime gammas, G^0 = -g^0
    anti_s = all(is_zero(Gs[a] * Gs[b] + Gs[b] * Gs[a] - 2 * ETA[a, b] * sp.eye(4)) for a in range(4) for b in range(4))
    herm = g0.H == -g0 and g0.H * g0 == sp.eye(4) and all(g.H == g for g in gj)
    check("clifford", anti and anti_s and herm,
          "{g^a, g^b} = 2 eta^ab for the note's g^a and for G = (-g^0, g^j); g^0 anti-Hermitian and unitary, "
          "g^j Hermitian")
    # The causal vector. a_j = g^0 g^j Hermitian with {a_i, a_j} = 2 delta_ij gives |xi| <= |eps|^2 for
    # xi^j = <eps, g^0 g^j eps>, since s.a has eigenvalues +-|s|. And <eps, g^0 g^0 eps> = -|eps|^2.
    e = sp.Matrix(sp.symbols("e0:4", complex=True))
    al = [g0 * g for g in gj]
    alpha_ok = all(a_.H == a_ for a_ in al) and all(
        is_zero(al[i] * al[j] + al[j] * al[i] - 2 * int(i == j) * sp.eye(4)) for i in range(3) for j in range(3))
    xi0 = sp.expand((e.H * g0 * g0 * e)[0] + (e.H * e)[0]) == 0
    check("causal", alpha_ok and xi0,
          "a_j = g^0 g^j Hermitian, {a_i, a_j} = 2 delta_ij (so |xi| <= |eps|^2), and <eps, g^0 g^0 eps> = -|eps|^2")


# ---- The Levi-Civita connection in the frame (n, e_i) -----------------------------------------------------------------
def christoffel(g, gi, crd):
    """Chr[m][n][l] = Gamma^m_nl of the metric g with typed-in inverse gi, in the coordinates crd."""
    n = len(crd)
    dg = [[[sp.diff(g[m, k], crd[l]) for l in range(n)] for k in range(n)] for m in range(n)]
    return [[[sp.expand(sum(gi[m, s] * (dg[s][l][k] + dg[s][k][l] - dg[k][l][s]) for s in range(n)) / 2)
              for l in range(n)] for k in range(n)] for m in range(n)]


def nabla(U, V, Chr, crd):
    """(nabla_U V)^m = U^n d_n V^m + Gamma^m_nl U^n V^l, with the sign of the Gamma term mutable."""
    sgn = -1 if MUT == "christoffel" else 1
    n = len(crd)
    return sp.Matrix([sp.expand(sum(U[k] * sp.diff(V[m], crd[k]) for k in range(n))
                                + sgn * sum(Chr[m][k][l] * U[k] * V[l] for k in range(n) for l in range(n)))
                      for m in range(n)])


def metric4(b, N, A):
    """The metric ds^2 = -N^2 dt^2 + sum_i A_i^2 (dx^i + b^i dt)^2 and its inverse, typed in."""
    a2 = [Ai**2 for Ai in A]
    g = sp.zeros(4)
    g[0, 0] = -N**2 + sum(a2[i] * b[i]**2 for i in range(3))
    for i in range(3):
        g[0, i + 1] = g[i + 1, 0] = a2[i] * b[i]
        g[i + 1, i + 1] = a2[i]
    gi = sp.zeros(4)
    gi[0, 0] = -1 / N**2
    for i in range(3):
        gi[0, i + 1] = gi[i + 1, 0] = b[i] / N**2
        for j in range(3):
            gi[i + 1, j + 1] = int(i == j) / a2[i] - b[i] * b[j] / N**2
    return g, gi


def rotate(Es, al):
    """The frame (Es[0], Es[1], Es[2]) rotated by the angle al about the third axis."""
    c, s = sp.cos(al), sp.sin(al)
    if MUT == "rot-ortho":
        return [c * Es[0] + s * Es[1], Es[1], Es[2]]
    return [c * Es[0] + s * Es[1], -s * Es[0] + c * Es[1], Es[2]]


def ricci(Chr, crd):
    """R_ab = d_c Gamma^c_ab - d_b Gamma^c_ac + Gamma^c_cd Gamma^d_ab - Gamma^c_bd Gamma^d_ac."""
    n = len(crd)
    return sp.Matrix(n, n, lambda a, b: sp.expand(sum(
        sp.diff(Chr[c][a][b], crd[c]) - sp.diff(Chr[c][a][c], crd[b])
        + sum(Chr[c][c][e] * Chr[e][a][b] - Chr[c][b][e] * Chr[e][a][c] for e in range(n)) for c in range(n))))


def connection(b, N=sp.Integer(1), psi=sp.Integer(1), A=None, alpha=None):
    """omega[a][b][c] = g(e_a, nabla_{e_c} e_b) for ds^2 = -N^2 dt^2 + sum A_i^2 (dx^i + b^i dt)^2, with
    e_0 = (d_t - b^i d_i)/N and e_i = d_i / A_i, A_i = psi^2 if A is not given (psi^4 delta, #67). With alpha, e_1 and
    e_2 are rotated by alpha about z (#86)."""
    conf = A is None
    A = [psi**2] * 3 if conf else A
    g, gi = metric4(b, N, A)
    inv_ok = is_zero(g * gi - sp.eye(4))  # inverse, typed in and verified
    Chr = christoffel(g, gi, CRD)
    sh = 1 if MUT == "frame" else -1
    lapse = 1 if MUT == "lapse-frame" else N
    if conf:
        sc = [psi**-1 if MUT == "psi-frame" else psi**-2] * 3
    else:
        sc = [1 / Ai for Ai in A]
        if MUT == "aniso-frame":
            sc[0] = 1 / A[1]
    Es = [sp.Matrix([int(m == i + 1) for m in range(4)]) * sc[i] for i in range(3)]
    if alpha is not None:
        Es = rotate(Es, alpha)
    E = [sp.Matrix([1] + [sh * bi for bi in b]) / lapse] + Es
    frame = sp.Matrix(4, 4, lambda a, c: sp.expand((E[a].T * g * E[c])[0]))
    dE = [[g * nabla(E[c], E[bb], Chr, CRD) for bb in range(4)] for c in range(4)]
    om = [[[sp.expand((E[a].T * dE[c][bb])[0]) for c in range(4)] for bb in range(4)] for a in range(4)]
    return inv_ok, frame, om


def section_planted():
    print("\n== 2. Planted truth: de Sitter in flat slicing ==")
    # x^i = a(t) X^i with a = exp(H t) turns -dt^2 + a^2 dX^2 into -dt^2 + sum (dx^i - H x^i dt)^2. The Eulerian
    # observers are comoving, so nabla_{e_i} n = H e_i and nabla_{e_i} e_j = H delta_ij n. Hence
    # omega_j0(e_i) = H delta_ij, omega_0j(e_i) = -H delta_ij and omega_jk(e_i) = 0, with K_ij = -H delta_ij.
    H = sp.symbols("H", positive=True)
    inv_ok, frame, om = connection([-H * X, -H * Y, -H * Z])
    want = all(sp.simplify(om[0][j][i] + H * int(i == j)) == 0 and sp.simplify(om[j][0][i] - H * int(i == j)) == 0
               and all(sp.simplify(om[j][k][i]) == 0 for k in range(1, 4))
               for i in range(1, 4) for j in range(1, 4))
    check("planted_desitter", inv_ok and want,
          "omega_0j(e_i) = -H delta_ij, omega_j0(e_i) = H delta_ij, omega_jk(e_i) = 0, as for comoving observers")


def spin(om, c, Gam, ss=1):
    """(1/4) omega_ab(e_c) G^a G^b, with the sign of the whole mutable."""
    return sp.expand(ss * sum((om[a][bb][c] * Gam[a] * Gam[bb] for a in range(4) for bb in range(4)),
                              sp.zeros(4)) / 4)


def section_slice():
    print("\n== 3. The flat slice with lapse 1 and a general shift b^i(t, x, y, z) ==")
    b = [sp.Function(f"b{i}")(*CRD) for i in (1, 2, 3)]
    inv_ok, frame, om = connection(b)
    check("orthonormal", inv_ok and is_zero(frame - ETA), "g(e_a, e_b) = eta_ab for e_0 = d_t - b^i d_i, e_i = d_i")
    anti = all(sp.expand(om[a][bb][c] + om[bb][a][c]) == 0 for a in range(4) for bb in range(4) for c in range(4))
    # K_ij = -g(nabla_{e_i} n, e_j), which for lapse 1 and g_ij = delta_ij is (d_i b_j + d_j b_i)/2.
    K = sp.Matrix(3, 3, lambda i, j: (sp.diff(b[j], CRD[i + 1]) + sp.diff(b[i], CRD[j + 1])) / 2)
    w0 = all(sp.expand(om[0][j + 1][i + 1] - K[i, j]) == 0 for i in range(3) for j in range(3))
    wjk = all(sp.expand(om[j][k][i]) == 0 for i in range(1, 4) for j in range(1, 4) for k in range(1, 4))
    check("omega_slice", anti and w0 and wjk,
          "omega_ab antisymmetric; on the slice omega_0j(e_i) = K_ij = (d_i b_j + d_j b_i)/2 and omega_jk(e_i) = 0")

    g0, gj = note_gammas()
    Gs = [-g0] + gj  # G^a, the gammas in which the spacetime connection is written
    Gl = [ETA[a, a] * Gs[a] for a in range(4)]  # G_a
    ss = -1 if MUT == "spin-sign" else 1

    # nabla_X (G(Y) eps) = G(nabla_X Y) eps + G(Y) nabla_X eps, i.e. [Gamma_X, G_d] = omega^a_d(X) G_a, with
    # omega^a_d = eta^aa omega_ad. This fixes the sign and the 1/4 of the spin connection.
    par = all(is_zero(spin(om, c, Gs, ss) * Gl[d] - Gl[d] * spin(om, c, Gs, ss)
                      - sum((ETA[a, a] * om[a][d][c] * Gl[a] for a in range(4)), sp.zeros(4)))
              for c in range(4) for d in range(4))
    check("clifford_parallel", par, "(1/4) omega_ab G^a G^b keeps Clifford multiplication parallel, on every e_c")

    sK = -1 if MUT == "sign" else 1
    stated = [sK * sum((K[i, j] * gj[j] * g0 for j in range(3)), sp.zeros(4)) / 2 for i in range(3)]
    res = all(is_zero(spin(om, i + 1, Gs, ss) - stated[i]) for i in range(3))
    check("restriction", res,
          "nabla_{e_i} = d_i + (1/4) omega_ab(e_i) G^a G^b with G^0 = -g^0 equals the stated "
          "nabla~_i = d_i %s (1/2) K_ij g^j g^0" % ("-" if sK < 0 else "+"))
    # In the note's own gammas, the same connection has the other sign (review.md section 3, gap 1).
    other = all(is_zero(spin(om, i + 1, [g0] + gj, ss) + sum((K[i, j] * gj[j] * g0 for j in range(3)), sp.zeros(4)) / 2)
                for i in range(3))
    check("note_gammas", other, "in the note's gammas themselves it is d_i - (1/2) K_ij g^j g^0")
    # The Witten operator of item 5, D~ = g^i nabla~_i = D + (1/2)(tr K) g^0.
    dirac = is_zero(sum((gj[i] * stated[i] for i in range(3)), sp.zeros(4)) - K.trace() * g0 / 2)
    check("dirac", dirac, "g^i (stated K term)_i = (1/2)(tr K) g^0, so D~ = D + (1/2)(tr K) g^0")


# ---- A lapse N and a conformally flat metric psi^4 delta -------------------------------------------------------------
def slice_data(N, psi, b, A=None, alpha=None):
    """The slice's own data from the 3-metric g_ij = A_i^2 delta_ij alone (A_i = psi^2 if A is not given), in the frame
    e_i = d_i / A_i, rotated by alpha about z if given. Returns the spin-connection components
    om3[j][k][i] = g(e_j, nabla_{e_i} e_k) and the ADM extrinsic curvature K_ij = (1/2N)(-d_t g_ij + D_i b_j + D_j b_i)
    in that frame, K_ij = K(e_i, e_j), with b_j = g_jk b^k."""
    X3 = CRD[1:]
    A = [psi**2] * 3 if A is None else A
    a2 = [Ai**2 for Ai in A]
    g3 = sp.diag(*a2)
    gi3 = sp.diag(*[1 / x for x in a2])  # typed in
    if MUT == "psi-slice":
        g3, gi3 = sp.diag(*A), sp.diag(*[1 / Ai for Ai in A])  # g_ij of the slice's own connection with the wrong power
    Chr = christoffel(g3, gi3, X3)
    E3 = [sp.Matrix([int(m == i) for m in range(3)]) / A[i] for i in range(3)]
    Er = E3 if alpha is None else rotate(E3, alpha)  # the frame of the result
    a0 = sp.Symbol("a0")  # a frozen angle, whose derivatives vanish
    Es = rotate(E3, a0) if alpha is not None and MUT == "rot-dalpha" else Er  # the frame nabla differentiates
    dE = [[g3 * nabla(Es[i], Es[k], Chr, X3) for k in range(3)] for i in range(3)]
    back = {a0: alpha} if alpha is not None else {}
    om3 = [[[sp.expand((Es[j].T * dE[i][k])[0]).xreplace(back) for i in range(3)] for k in range(3)]
           for j in range(3)]
    bl = [a2[i] * b[i] for i in range(3)]
    lapse = 1 if MUT == "lapse-K" else N

    def D(i, j):  # D_i b_j
        return sp.diff(bl[j], X3[i]) - sum(Chr[k][i][j] * bl[k] for k in range(3))

    sA = -1 if MUT == "adm-sign" else 1
    Kc = sp.Matrix(3, 3, lambda i, j: sp.expand(sA * (-sp.diff(a2[i] * int(i == j), T) + D(i, j) + D(j, i)) / (2 * lapse)))
    Ek = E3 if MUT == "rot-k" else Er
    return om3, sp.Matrix(3, 3, lambda i, j: sp.expand((Ek[i].T * Kc * Ek[j])[0]))


def table(nz):
    """A 4x4x4 table [a][b][c] of omega_ab(e_c), zero except where nz[(a, b, c)] is given."""
    return [[[nz.get((a, bb, c), 0) for c in range(4)] for bb in range(4)] for a in range(4)]


def same(om, want):
    return all(zero(om[a][bb][c] - want[a][bb][c]) for a in range(4) for bb in range(4) for c in range(4))


def planted(name, N, psi, b, want, wantK, text, A=None):
    """Passes if the 4-metric's omega equals the known table and the ADM K equals the known K."""
    inv_ok, frame, om = connection(b, N, psi, A)
    om3, K = slice_data(N, psi, b, A)
    check(name, inv_ok and is_zero(frame - ETA) and same(om, want) and is_zero(K - wantK), text)


def section_planted_general():
    print("\n== 4. Planted truth with a lapse and a conformal factor ==")
    # Each metric has a spin connection known without the code below. The checks use the code path of section 5.
    # (a) de Sitter in flat slicing again, in the coordinates x^i = s(u) x'^i, t = tau(u), which keep the slices and the
    #     comoving observers. Then N = tau', psi^4 = s^2 and b^i = (s'/s - H tau') x^i, and the frame (n, e_i) is the
    #     one of section 2, so omega_0j(e_i) = -H delta_ij, omega_j0(e_i) = H delta_ij, the rest zero and
    #     K_ij = -H delta_ij.
    H = sp.symbols("H", positive=True)
    tau, s = sp.Function("tau")(T), sp.Function("s")(T)
    N, psi = sp.diff(tau, T), sp.sqrt(s)
    nz = {}
    for i in range(1, 4):
        nz[(0, i, i)], nz[(i, 0, i)] = -H, H
    planted("planted_ds_rescaled", N, psi, [(sp.diff(s, T) / s - H * N) * v for v in (X, Y, Z)], table(nz),
            -H * sp.eye(3), "de Sitter in x = s(t) x', t = tau(u): omega_0j(e_i) = -H delta_ij, omega_jk = 0, "
            "K_ij = -H delta_ij, for arbitrary tau, s")
    # (b) g = Omega^2 eta for any Omega(t, x, y, z): N = Omega, psi^2 = Omega, no shift. With g = e^(2 phi) eta,
    #     theta^a = e^phi dx^a gives omega_ab = phi_b dx_a - phi_a dx_b, so omega_ab(e_c) = e^-phi (phi_b eta_ac -
    #     phi_a eta_bc) with phi = ln Omega. On the slice omega_0j(e_i) = -e^-phi phi_t delta_ij = K_ij
    #     = -(Om_t/Om^2) delta_ij.
    Om = sp.Function("Omega")(*CRD)
    ph = [sp.diff(Om, v) / Om for v in CRD]
    nz = {(a, bb, c): (ph[bb] * ETA[a, c] - ph[a] * ETA[bb, c]) / Om for a in range(4) for bb in range(4)
          for c in range(4)}
    planted("planted_conf", Om, sp.sqrt(Om), [0, 0, 0], table(nz), -sp.diff(Om, T) / Om**2 * sp.eye(3),
            "g = Omega^2 eta: omega_ab(e_c) = Omega^-1 (phi_b eta_ac - phi_a eta_bc), phi = ln Omega, "
            "K_ij = -(d_t Omega / Omega^2) delta_ij")
    # (c) Static, N(x, y, z) and psi(x, y, z) independent, no shift. n is the Killing direction normalised, so
    #     nabla_n n = D ln N, nabla_{e_i} n = 0 and K = 0, giving omega_0j(e_0) = -e_j ln N and omega_jk(e_0) = 0.
    #     The slice's omega_jk(e_i) = psi^-2 (delta_ij phi_k - delta_ik phi_j) with phi = 2 ln psi is the conformally
    #     flat one, from omega_jk = phi_k dx^j - phi_j dx^k.
    Ns, ps_ = sp.Function("N")(X, Y, Z), sp.Function("psi")(X, Y, Z)
    ph = [2 * sp.diff(ps_, v) / ps_ for v in (X, Y, Z)]
    nz = {}
    for j in range(1, 4):
        nz[(0, j, 0)] = -sp.diff(Ns, CRD[j]) / (Ns * ps_**2)
        nz[(j, 0, 0)] = -nz[(0, j, 0)]
        for k in range(1, 4):
            for i in range(1, 4):
                nz[(j, k, i)] = (int(i == j) * ph[k - 1] - int(i == k) * ph[j - 1]) / ps_**2
    planted("planted_static", Ns, ps_, [0, 0, 0], table(nz), sp.zeros(3),
            "static N(x), psi(x): omega_0j(e_0) = -e_j ln N, omega_jk(e_i) = psi^-2 (delta_ij phi_k - delta_ik phi_j), "
            "K = 0")
    # (d) Spatially homogeneous, N(t) and psi(t) independent. Then e_0 is geodesic, e_i is parallel on the slice and
    #     K_ij = -(1/N) d_t(psi^2) / psi^2 delta_ij = -(2 psi'/(N psi)) delta_ij.
    Nt, pt = sp.Function("N")(T), sp.Function("psi")(T)
    Kf = -2 * sp.diff(pt, T) / (Nt * pt)
    nz = {}
    for i in range(1, 4):
        nz[(0, i, i)], nz[(i, 0, i)] = Kf, -Kf
    planted("planted_flrw", Nt, pt, [0, 0, 0], table(nz), Kf * sp.eye(3),
            "N(t), psi(t): omega_0j(e_i) = K_ij = -(2 psi'/(N psi)) delta_ij, every other omega zero")


def check_general(sfx, frame_text, who, inv_ok, frame, om, om3, K):
    """The checks of sections 5, 7 and 8. om is the 4-metric's omega_ab(e_c) in the frame (n, e_i), and om3 and K are the
    slice's own spin connection and the ADM K in the frame e_i, from the 3-metric alone."""
    check("orthonormal_" + sfx, inv_ok and is_zero(frame - ETA), frame_text)
    # K_ij = -g(nabla_{e_i} n, e_j) = -omega_j0(e_i), from the 4-metric. The ADM formula and the slice's 3-metric give
    # the same K, and the slice's omega_jk(e_i) from the 3-metric alone is the 4-metric's.
    anti = all(zero(om[a][bb][c] + om[bb][a][c]) for a in range(4) for bb in range(4) for c in range(4))
    Kdef = sp.Matrix(3, 3, lambda i, j: -om[j + 1][0][i + 1])
    k_ok = is_zero(Kdef - K) and is_zero(K - K.T)
    w0 = all(zero(om[0][j + 1][i + 1] - K[i, j]) for i in range(3) for j in range(3))
    wjk = all(zero(om[j + 1][k + 1][i + 1] - om3[j][k][i]) for i in range(3) for j in range(3) for k in range(3))
    check("omega_" + sfx, anti and k_ok and w0 and wjk,
          "omega_ab antisymmetric; K_ij = -g(nabla_{e_i} n, e_j) equals the ADM K; omega_0j(e_i) = K_ij; "
          "omega_jk(e_i) is the slice's own")

    g0, gj = note_gammas()
    Gs = [-g0] + gj
    Gl = [ETA[a, a] * Gs[a] for a in range(4)]
    ss = -1 if MUT == "spin-sign" else 1
    sp_c = [spin(om, c, Gs, ss) for c in range(4)]
    par = all(is_zero(sp_c[c] * Gl[d] - Gl[d] * sp_c[c] - sum((ETA[a, a] * om[a][d][c] * Gl[a] for a in range(4)),
                                                              sp.zeros(4)))
              for c in range(4) for d in range(4))
    check("clifford_parallel_" + sfx, par,
          "(1/4) omega_ab G^a G^b keeps Clifford multiplication parallel on every e_c, for " + who)

    # nabla_i of the slice, in the note's gammas: d_i + (1/4) omega3_jk(e_i) g^j g^k.
    sK = -1 if MUT == "sign" else 1
    slice_i = [sum((om3[j][k][i] * gj[j] * gj[k] for j in range(3) for k in range(3)), sp.zeros(4)) / 4
               for i in range(3)]
    Kterm = [sum((K[i, j] * gj[j] * g0 for j in range(3)), sp.zeros(4)) / 2 for i in range(3)]
    stated = [slice_i[i] + sK * Kterm[i] for i in range(3)]
    res = all(is_zero(sp_c[i + 1] - stated[i]) for i in range(3))
    check("restriction_" + sfx, res,
          "nabla_{e_i} = d_i + (1/4) omega_ab(e_i) G^a G^b with G^0 = -g^0 equals nabla_i %s (1/2) K_ij g^j g^0, "
          "nabla_i the slice's spin connection" % ("-" if sK < 0 else "+"))
    other = all(is_zero(spin(om, i + 1, [g0] + gj, ss) - (slice_i[i] - Kterm[i])) for i in range(3))
    check("note_gammas_" + sfx, other, "in the note's gammas themselves it is nabla_i - (1/2) K_ij g^j g^0")
    # D~ = g^i nabla~_i = D + (1/2)(tr K) g^0, with D = g^i nabla_i the slice's Dirac operator.
    dirac = is_zero(sum((gj[i] * (stated[i] - slice_i[i]) for i in range(3)), sp.zeros(4)) - K.trace() * g0 / 2)
    check("dirac_" + sfx, dirac, "g^i (stated K term)_i = (1/2)(tr K) g^0, so D~ = D + (1/2)(tr K) g^0")


def section_general():
    print("\n== 5. A general slice: lapse N, conformal factor psi, shift b^i, all of (t, x, y, z) ==")
    N, psi = sp.Function("N")(*CRD), sp.Function("psi")(*CRD)
    b = [sp.Function(f"b{i}")(*CRD) for i in (1, 2, 3)]
    inv_ok, frame, om = connection(b, N, psi)
    om3, K = slice_data(N, psi, b)
    check_general("general", "g(e_a, e_b) = eta_ab for e_0 = (d_t - b^i d_i)/N, e_i = psi^-2 d_i",
                  "the general metric", inv_ok, frame, om, om3, K)


def section_planted_diag():
    print("\n== 6. Planted truth with a diagonal slice ==")
    # (a) Bianchi I, N(t) and A_i = a_i(t) independent, no shift. The slice is flat in the frame e_i = d_i / a_i, which is
    #     parallel on it, and n is geodesic. With H_i = a_i'/(N a_i), nabla_{e_i} n = H_i e_i, so omega_j0(e_i) = H_i
    #     delta_ij = -omega_0j(e_i), omega_jk = 0 and K_ij = -H_i delta_ij.
    a = [sp.Function(f"a{i}")(T) for i in (1, 2, 3)]
    Nt = sp.Function("N")(T)
    H = [sp.diff(ai, T) / (Nt * ai) for ai in a]
    nz = {}
    for i in range(1, 4):
        nz[(0, i, i)], nz[(i, 0, i)] = -H[i - 1], H[i - 1]
    planted("planted_bianchi", Nt, 1, [0, 0, 0], table(nz), -sp.diag(*H),
            "Bianchi I, N(t), a_i(t): omega_0j(e_i) = K_ij = -H_i delta_ij with H_i = a_i'/(N a_i), every other omega zero",
            A=a)
    # (b) Kasner, N = 1, a_i = t^p_i with sum p = sum p^2 = 1, here p = (-2/7, 3/7, 6/7), three different exponents.
    #     It is the Bianchi I of (a) with H_i = p_i/t, and it must be vacuum, R_ab = 0.
    p = [sp.Rational(1, 2)] * 3 if MUT == "kasner" else [sp.Rational(-2, 7), sp.Rational(3, 7), sp.Rational(6, 7)]
    A = [T**pi for pi in p]
    nz = {}
    for i in range(1, 4):
        nz[(0, i, i)], nz[(i, 0, i)] = -p[i - 1] / T, p[i - 1] / T
    g, gi = metric4([0, 0, 0], sp.Integer(1), A)
    vac = is_zero(ricci(christoffel(g, gi, CRD), CRD))
    inv_ok, frame, om = connection([0, 0, 0], A=A)
    om3, K = slice_data(sp.Integer(1), 1, [0, 0, 0], A)
    check("planted_kasner", vac and inv_ok and is_zero(frame - ETA) and same(om, table(nz))
          and is_zero(K + sp.diag(*[pi / T for pi in p])),
          "Kasner t^p_i, p = (-2/7, 3/7, 6/7): R_ab = 0, omega_0j(e_i) = K_ij = -(p_i/t) delta_ij, every other omega zero")
    # (c) A diagonal metric with N and A_i all of (t, x, y, z), no shift. By Cartan, theta^j = A_j dx^j gives
    #     omega_jk(e_i) = delta_ij d_k A_j / (A_j A_k) - delta_ik d_j A_k / (A_j A_k) for j != k, which is the polar frame
    #     of the flat plane, omega_21(e_2) = 1/r, for A = (1, r, 1). The Eulerian observers have acceleration D ln N and
    #     omega_jk(e_0) = 0, and as in (a) omega_0j(e_i) = K_ij = -d_t A_i/(N A_i) delta_ij.
    Ns = sp.Function("N")(*CRD)
    A = [sp.Function(n)(*CRD) for n in "ABC"]
    nz = {}
    for j in range(1, 4):
        nz[(0, j, 0)] = -sp.diff(Ns, CRD[j]) / (Ns * A[j - 1])
        nz[(j, 0, 0)] = -nz[(0, j, 0)]
        Hj = sp.diff(A[j - 1], T) / (Ns * A[j - 1])
        nz[(0, j, j)], nz[(j, 0, j)] = -Hj, Hj
        for k in range(1, 4):
            for i in range(1, 4):
                if j != k:
                    nz[(j, k, i)] = (int(i == j) * sp.diff(A[j - 1], CRD[k]) - int(i == k) * sp.diff(A[k - 1], CRD[j])) \
                        / (A[j - 1] * A[k - 1])
    planted("planted_diag_cartan", Ns, 1, [0, 0, 0], table(nz),
            -sp.diag(*[sp.diff(Ai, T) / (Ns * Ai) for Ai in A]),
            "diag(A^2, B^2, C^2), N, all of (t, x, y, z): omega_jk(e_i) by Cartan's formula, omega_0j(e_0) = -e_j ln N, "
            "omega_0j(e_i) = K_ij = -d_t A_i/(N A_i) delta_ij", A=A)


def diag_fields():
    N = sp.Function("N")(*CRD)
    A = [sp.Function(n)(*CRD) for n in "ABC"]
    b = [sp.Function(f"b{i}")(*CRD) for i in (1, 2, 3)]
    return N, A, b


def section_diag():
    print("\n== 7. A diagonal slice: lapse N, scale factors A, B, C and shift b^i, all of (t, x, y, z) ==")
    N, A, b = diag_fields()
    inv_ok, frame, om = connection(b, N, A=A)
    om3, K = slice_data(N, 1, b, A)
    check_general("diag", "g(e_a, e_b) = eta_ab for e_0 = (d_t - b^i d_i)/N, e_i = d_i / A_i", "the diagonal metric",
                  inv_ok, frame, om, om3, K)


def section_rot():
    print("\n== 8. The diagonal slice in the frame rotated by alpha(t, x, y, z) about z ==")
    N, A, b = diag_fields()
    inv_ok, frame, om = connection(b, N, A=A, alpha=ALPHA)
    om3, K = slice_data(N, 1, b, A, alpha=ALPHA)
    check_general("rot", "g(e'_a, e'_b) = eta_ab for e'_0 = n, e'_1 = cos(alpha) e_1 + sin(alpha) e_2, "
                  "e'_2 = -sin(alpha) e_1 + cos(alpha) e_2, e'_3 = e_3", "the rotated frame",
                  inv_ok, frame, om, om3, K)


def main():
    t0 = time.time()
    print("conventions.py, issues #62, #67 and #86. sympy %s, python %s. mutation: %s"
          % (sp.__version__, sys.version.split()[0], MUT or "none"))
    section_gammas()
    section_planted()
    section_slice()
    section_planted_general()
    section_general()
    section_planted_diag()
    section_diag()
    section_rot()
    fails = [n for n, ok in RESULTS if not ok]
    print("\nchecks: %d passed, %d failed %s. %.0f s" % (len(RESULTS) - len(fails), len(fails), fails,
                                                         time.time() - t0))
    if MUT:
        target = CONTROLS[MUT][1]
        print("control %s: targets %s %s" % (MUT, target, "FIRED" if all(n in fails for n in target) else "SILENT"))
        sys.exit(1)
    status = 1 if fails else 0
    if "no-controls" not in OPTS:
        print("\n== Controls: each mutation must make its target checks fail ==")
        for name, (what, target) in CONTROLS.items():
            p = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True)
            failed = [ln.split("]")[1].split(":")[0].strip() for ln in p.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = p.returncode == 1 and all(n in failed for n in target)
            status |= 0 if fired else 1
            print("[%s] --mutate=%s (%s): exit %d, failed %s" % ("FIRED" if fired else "SILENT", name, what,
                                                                p.returncode, failed), flush=True)
        print("controls: %s" % ("all fired" if status == 0 else "SOME SILENT OR CHECK FAILED"))
    sys.exit(status)


if __name__ == "__main__":
    main()

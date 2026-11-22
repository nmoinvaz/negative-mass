#!/usr/bin/env python3
"""Item 6's Weitzenbock identity on a conformally flat slice (#87) and on a diagonal slice (#99), checked symbolically.

Data, #87: g = psi^4 delta on R^3 with psi(x, y, z) arbitrary, an arbitrary symmetric K_ij(x, y, z) in the orthonormal frame
e_i = psi^-2 d_i, and an arbitrary spinor field eps = u + i v with the eight real components u_a, v_a of (x, y, z). Data,
#99: g = diag(A^2, B^2, C^2)(x, y, z) with A, B, C arbitrary positive functions, e_i = d_i / A_i, K_ij = K(e_i, e_j)
arbitrary and symmetric, eps arbitrary, and div V = (1/ABC) d_a(ABC V^a / A_a) for V = V^a e_a. A = B = C = psi^2 is #87's
slice. The note's gammas (mass_bounds.md, "Conventions") act on eps, <a, b> = a^dagger b, xi^j = <eps, g^0 g^j eps>, and
    nabla~_i = nabla_i + (1/2) K_ij g^j g^0,   nabla_i = e_i + (1/4) omega_jk(e_i) g^j g^k,   D~ = g^i nabla~_i,
with the sign of the K term fixed by #67 and #86. mu and J_i are defined by the constraints,
    16 pi mu = R + (tr K)^2 - K_ij K^ij,    8 pi J_i = D^j (K_ij - g_ij tr K),
which with K_ij = -g(nabla_{e_i} n, e_j) and J_i = -n^a T_ai are #58 item 1's conventions.

The identity checked is
    |nabla~ eps|^2 - |D~ eps|^2 + 4 pi (mu |eps|^2 - J.xi) = div V,
    V^i = Re<eps, nabla^i eps - g^i D eps> - (1/2) P_ij xi^j,   P_ij = K_ij - delta_ij tr K,   D = g^i nabla_i,
where div V = psi^-6 d_a(psi^4 V^a) on psi^4 delta, and (1/ABC) d_a(ABC V^a / A_a) on the diagonal slice. By hand (proved) it
is Lichnerowicz's formula |nabla eps|^2 - |D eps|^2 + (R/4)|eps|^2 = div Re<eps, nabla^i eps - g^i D eps> for the spin part,
and for the K part Re<nabla_i eps, P_ij g^j g^0 eps> = (1/2) div(P_ij <eps, g^j g^0 eps>) - (1/2)(nabla_i P_ij)<eps, g^j g^0
eps>, with <eps, g^j g^0 eps> = -xi^j, nabla_i P_ij = 8 pi J_j and the algebraic terms (1/4)(|K|^2 - (tr K)^2)|eps|^2. Using
the Hamiltonian constraint, -R/4 + (|K|^2 - (tr K)^2)/4 = -4 pi mu, which gives the identity. The script checks the identity
as a polynomial identity in the jets of psi (or A, B, C), K_ij, u_a and v_a, with V written out. The difference is expanded
and must be zero term by term. Inside the integral with D~ eps = 0 this is item 6 (the boundary term is div V).

Planted truth comes first. Each of these has an answer known without the code and goes through the code path of the
general case. On psi^4 delta: R = -8 psi^-5 Laplace psi, the spin connection omega_jk(e_i) = psi^-2 (delta_ij phi_k -
delta_ik phi_j) with phi = 2 ln psi, #58 item 1's mu and J for Alcubierre's flat data (with an arbitrary profile
f(x, y, z)), K = 0, K = f delta, a transverse trace-free K (Hessian of a harmonic polynomial, times psi^-6), and
Schwarzschild. Then the constraints are compared with G_ab of the 4-metric -dt^2 + psi^4 (delta_ij - 2 t K_ij) dx^i dx^j,
whose Einstein tensor at t = 0 gives 8 pi mu = G_tt and 8 pi J_i = -G_ti. On the diagonal slice (sections 6 to 9): R of six
metrics with a hand value (flat, dx^2 + dy^2 + c(x,y)^2 dz^2, dx^2 + f^2 (dy^2 + dz^2), dx^2 + f^2 dy^2 + g^2 dz^2,
A(x), B(x), C(x) arbitrary, and A = B = C = psi^2), omega_jk(e_i) = (delta_ij d_k A_j - delta_ik d_j A_k)/(A_j A_k) for
j != k from Cartan's formula (conventions.py section 6(c), #86), Clifford parallelism, mu and J for Alcubierre's data with
A = B = C = 1 (the planted flat case), K = 0 and K = f delta on two metrics with a hand R, the Bowen-York form on
A = B = C = psi^2, a Kasner slice with constant a_i, and the 4-metric -dt^2 + Lam (delta - 2 t K) Lam dx dx, Lam =
diag(A, B, C). The identities then run in the order flat (A = B = C = 1, K = 0), Lichnerowicz (K = 0, A, B, C arbitrary),
flat with K arbitrary, Alcubierre's K, A = B = C = psi^2 with K arbitrary (#87's identity_general through the diagonal
code), K = diag(k1, k2, k3), and last A, B, C, K_ij and eps all arbitrary.

What the checks share. The gammas are retyped from conventions.py (compared with its note_gammas() in `gammas_match`), and
christoffel and ricci are retyped from it too. Nothing else is imported from it. The identity side (matrices, the frame, V)
shares with the constraint side only sympy, christoffel and ricci. The spin connection of the identity comes from the
Christoffel symbols of the metric and is compared with the closed form above. The constraints are computed in coordinates
and compared with the 4-metric's G_ab, which shares christoffel and ricci with them, so a wrong ricci is seen only by the
planted values of R. The diagonal code (GeoDiag, constraints_diag, identity_diag) is written beside #87's Geo, constraints
and identity and calls none of them, so #87's checks run as they were. They share the helpers (christoffel, ricci, inner,
re, zero, sym_K, spinor, note_gammas), and `diag_conf_limit` and `identity_diag_conf` compare the two codes at
A = B = C = psi^2. Complex conjugation is the structural replacement i -> -i on expressions in the real functions u_a, v_a,
psi, A, B, C and K_ij. Zero is tested by sp.expand, which is canonical for Laurent polynomials in the jets, and sp.simplify
only on the planted cases with sqrt.

What the code proves (*checked symbolically*): at every point, for all psi, K_ij and eps on psi^4 delta, and for all A, B,
C, K_ij and eps on diag(A^2, B^2, C^2)(x, y, z), the two sides above are equal as polynomials in the jets, to second order,
of the functions, in the note's gammas with the stated sign of nabla~ and with mu and J_i the constraints' values. It does
not prove the integrated identity, whose boundary term is div V and needs decay (review.md gap 2), or that D~ eps = 0 has a
solution (item 5). The hand derivation above uses only the Levi-Civita connection, the Clifford relations and Hermiticity,
so it holds on any slice (*proved* by hand). Still covered only by it: slices that are not diagonal in the coordinates, and
frames other than e_i = d_i / A_i (such as one rotated about an axis). Assumed: sp.expand is canonical on these Laurent
polynomials, and the constraints define mu and J_i with K_ij = -g(nabla_{e_i} n, e_j), which `constraints_4d` and
`constraints_4d_diag` tie to G_ab on these classes of slices (N = 1, no shift).

    python3 notes/warp/weitzenbock_conf.py                  checks, then every control
    python3 notes/warp/weitzenbock_conf.py --no-controls    checks only
    python3 notes/warp/weitzenbock_conf.py --mutate=NAME    one control (its target checks must fail, exit 1)
    python3 notes/warp/weitzenbock_conf.py --no-controls --only=NAME,...   (or --skip=NAME,...) some diagonal checks of
        sections 6 to 9 only, for staged runs. Under --mutate a diagonal check runs only if the control targets it.

Exit status is 1 if a check fails or a control does not fire. One process, 20 s for the checks and about 10 minutes with
the controls, 75 MB.
"""
import importlib.util
import pathlib
import subprocess
import sys
import time

import sympy as sp

OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate", "")
ONLY = set(OPTS["only"].split(",")) if "only" in OPTS else None  # staging: run only these diagonal checks
SKIP = set(OPTS.get("skip", "").split(","))  # staging: leave these diagonal checks out

IDENT = ["identity_flat_K0", "identity_lichnerowicz", "identity_flat_K", "identity_alcubierre", "identity_general"]
# The diagonal slice's checks (#99). IDENT_D are its identities, KID those with K != 0 and ANI those with A, B, C not constant.
IDENT_D = ["identity_diag_flat_K0", "identity_diag_lichnerowicz", "identity_diag_flat_K", "identity_diag_alcubierre",
           "identity_diag_conf", "identity_diag_Kdiag", "identity_diag_general"]
KID = ["identity_diag_flat_K", "identity_diag_alcubierre", "identity_diag_conf", "identity_diag_Kdiag",
       "identity_diag_general"]
ANI = ["identity_diag_lichnerowicz", "identity_diag_conf", "identity_diag_Kdiag", "identity_diag_general"]
CONTROLS = {
    "sign": ("takes nabla~ with the other sign, d_i - (1/2) K_ij g^j g^0, and D~ from it",
             ["dirac", "identity_flat_K", "identity_alcubierre", "identity_general"] + KID),
    "Jsign": ("writes +J.xi in the identity, instead of -J.xi",
              ["identity_flat_K", "identity_alcubierre", "identity_general"] + KID),
    "musign": ("writes -mu |eps|^2 in the identity", ["identity_lichnerowicz", "identity_flat_K",
                                                      "identity_alcubierre", "identity_general",
                                                      "identity_diag_lichnerowicz"] + KID),
    "Kdiv": ("drops the term -(1/2) P_ij xi^j from V",
             ["identity_flat_K", "identity_alcubierre", "identity_general", "V_reduced"] + KID + ["V_reduced_diag"]),
    "Vspin": ("drops -g^i D eps from V", IDENT + ["V_reduced"] + IDENT_D + ["V_reduced_diag"]),
    "weight": ("takes div V = d_a V^a, without psi^-6 d_a psi^4 (diagonal: without (1/ABC) d_a ABC)",
               ["identity_lichnerowicz", "identity_general"] + ANI),
    "omega-sign": ("takes the slice's spin connection with -(1/4) omega_jk g^j g^k",
                   ["clifford_parallel", "identity_lichnerowicz", "identity_general", "clifford_parallel_diag"] + ANI),
    "christoffel": ("flips the sign of the Christoffel term in nabla_X Y, which gives omega_jk",
                    ["planted_omega", "clifford_parallel", "identity_lichnerowicz", "identity_general",
                     "planted_omega_diag", "clifford_parallel_diag"] + ANI),
    "psi-frame": ("takes e_i = psi^-1 d_i instead of psi^-2 d_i",
                  ["orthonormal", "planted_omega", "identity_lichnerowicz", "identity_general", "diag_conf_limit"]),
    "trace": ("takes P_ij = K_ij + g_ij tr K in the momentum constraint",
              ["planted_alcubierre", "planted_trace", "constraints_4d", "identity_flat_K", "identity_alcubierre",
               "identity_general", "planted_alcubierre_diag", "planted_trace_diag", "constraints_4d_diag"] + KID),
    "Ksq": ("takes 16 pi mu = R + (tr K)^2 + K_ij K^ij",
            ["planted_alcubierre", "planted_trace", "planted_TT", "constraints_4d", "identity_flat_K",
             "identity_alcubierre", "identity_general", "planted_alcubierre_diag", "planted_trace_diag",
             "planted_TT_diag", "planted_kasner_diag", "constraints_4d_diag"] + KID),
    "Rsign": ("takes 16 pi mu = -R + (tr K)^2 - K_ij K^ij",
              ["planted_K0", "planted_trace", "planted_TT", "constraints_4d", "identity_lichnerowicz",
               "identity_general", "planted_K0_diag", "planted_trace_diag", "planted_TT_diag", "constraints_4d_diag"]
              + ANI),
    "Kframe": ("takes the frame components K_ij as coordinate components in the constraints",
               ["planted_trace", "planted_TT", "constraints_4d", "identity_general", "planted_trace_diag",
                "planted_TT_diag", "planted_kasner_diag", "constraints_4d_diag", "identity_diag_conf",
                "identity_diag_Kdiag", "identity_diag_general"]),
    "J-christoffel": ("drops the Christoffel terms of the covariant divergence in the momentum constraint",
                      ["planted_trace", "planted_TT", "constraints_4d", "identity_general", "planted_trace_diag",
                       "planted_TT_diag", "constraints_4d_diag", "identity_diag_conf", "identity_diag_Kdiag",
                       "identity_diag_general"]),
    "gamma": ("scales g^1 by 11/10", ["clifford", "gammas_match", "causal", "clifford_parallel", "dirac", "V_reduced"]
              + IDENT + ["clifford_parallel_diag"] + IDENT_D + ["V_reduced_diag"]),
    "no-conj": ("takes <a, b> = a^T b, without the complex conjugate", ["dagger"] + IDENT + IDENT_D),
    "adm-K": ("takes the 4-metric with the other sign of K, g_ij = psi^4 (delta_ij + 2 t K_ij)",
              ["constraints_4d", "constraints_4d_diag"]),
    "J4d": ("takes J_i = +G_ti / 8 pi from the 4-metric", ["constraints_4d", "constraints_4d_diag"]),
    "inv-typo": ("types the inverse metric as psi^-2 delta (diagonal: as diag(A^-1, B^-1, C^-1))",
                 ["inverse", "inverse_diag", "planted_R_diag", "planted_omega_diag", "clifford_parallel_diag",
                  "planted_K0_diag", "planted_trace_diag", "planted_TT_diag", "constraints_4d_diag"] + ANI),
    "ricci-sign": ("takes the Gamma Gamma terms of R_ab with the other sign",
                   ["planted_R", "planted_K0", "planted_schwarzschild", "constraints_4d", "planted_R_diag",
                    "planted_K0_diag", "planted_trace_diag", "planted_TT_diag", "constraints_4d_diag"] + ANI),
    "schwarzschild-power": ("takes psi = 1 + m1/2|x|^2 + m2/2|x - c|^2, which is not harmonic",
                            ["planted_schwarzschild"]),
    "diag-frame": ("takes e_i = d_i / A_i^2 instead of d_i / A_i",
                   ["orthonormal_diag", "planted_omega_diag"] + ANI),
    "aniso-frame": ("takes e_1 = d_1 / B, the scale of the wrong axis",
                    ["orthonormal_diag", "planted_omega_diag", "identity_diag_lichnerowicz", "identity_diag_Kdiag",
                     "identity_diag_general"]),
    "abc": ("takes sqrt(det g) = A^3 in the divergence, not ABC (equal for A = B = C)",
            ["identity_diag_lichnerowicz", "identity_diag_Kdiag", "identity_diag_general"]),
    "omega-hand": ("takes the second term of the closed form of omega_jk(e_i) with the other sign", ["planted_omega_diag"]),
    "R-hand": ("takes -2 f'^2/f^2 as +2 f'^2/f^2 in the hand value of R for dx^2 + f^2 (dy^2 + dz^2)", ["planted_R_diag"]),
    "kasner-exponents": ("takes the exponents (1/2, 1/2, 1/2), which are not vacuum, for the Kasner slice",
                         ["planted_kasner_diag"]),
}
if MUT and MUT not in CONTROLS:
    sys.exit(f"unknown mutation {MUT}")
RESULTS = []

T, X, Y, Z = sp.symbols("t x y z", real=True)
CRD = (T, X, Y, Z)
X3 = (X, Y, Z)
ETA = sp.diag(-1, 1, 1, 1)
I3 = sp.eye(3)


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def fn(name):
    return sp.Function(name)(*X3)


def zero(e):
    return sp.expand(e) == 0


def is_zero(M):
    return all(zero(e) for e in sp.Matrix(M))


def nterms(e):
    return len(sp.Add.make_args(sp.expand(e)))


def sym_K():
    """An arbitrary symmetric K_ij(x, y, z), six functions."""
    return sp.Matrix(3, 3, lambda i, j: fn("K%d%d" % (min(i, j) + 1, max(i, j) + 1)))


# ---- Complex conjugation and the inner product --------------------------------------------------------------------------
def conj(e):
    """The complex conjugate of an expression in the real functions u_a, v_a, psi, K_ij and real coordinates."""
    e = sp.expand(e)
    return e if MUT == "no-conj" else e.xreplace({sp.I: -sp.I})


def inner(a, b):
    """<a, b> = a^dagger b for column vectors a and b."""
    return sp.expand((a.applyfunc(conj).T * b)[0])


def re(e):
    return sp.expand((e + conj(e)) / 2)


def spinor():
    """eps = u + i v, with u_a and v_a real functions of (x, y, z)."""
    return sp.Matrix([fn(f"u{a}") + sp.I * fn(f"v{a}") for a in range(4)])


# ---- The gammas of mass_bounds.md, "Conventions" -------------------------------------------------------------------------
def note_gammas():
    s = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
    z2 = sp.zeros(2)
    g0 = sp.I * sp.diag(1, 1, -1, -1)
    gj = [sp.Matrix(sp.BlockMatrix([[z2, sj], [sj, z2]])) for sj in s]
    if MUT == "gamma":
        gj[0] = gj[0] * sp.Rational(11, 10)
    return g0, gj


def section_gammas():
    print("\n== 1. The gammas, the conjugation and the causal vector ==")
    g0, gj = note_gammas()
    G = [g0] + gj
    anti = all(is_zero(G[a] * G[b] + G[b] * G[a] - 2 * ETA[a, b] * sp.eye(4)) for a in range(4) for b in range(4))
    herm = g0.H == -g0 and g0.H * g0 == sp.eye(4) and all(g.H == g for g in gj)
    check("clifford", anti and herm,
          "{g^a, g^b} = 2 eta^ab; g^0 anti-Hermitian and unitary, g^j Hermitian")
    # The same gammas as conventions.py, where #67 and #86 fixed the sign of nabla~. Its argv is hidden while it loads.
    spec = importlib.util.spec_from_file_location("conventions", pathlib.Path(__file__).with_name("conventions.py"))
    mod = importlib.util.module_from_spec(spec)
    saved, sys.argv = sys.argv, sys.argv[:1]
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.argv = saved
    c0, cj = mod.note_gammas()
    check("gammas_match", c0 == g0 and all(a == b for a, b in zip(cj, gj)),
          "the gammas typed here are conventions.note_gammas() entry by entry")
    # a_j = g^0 g^j is Hermitian with {a_i, a_j} = 2 delta_ij, so |xi| <= |eps|^2 for xi^j = <eps, a_j eps>. And
    # g^j g^0 = -g^0 g^j is Hermitian too, so <eps, g^j g^0 eps> = -xi^j. And <eps, g^0 g^0 eps> = -|eps|^2.
    al = [g0 * g for g in gj]
    be = [g * g0 for g in gj]
    ok = all(a_.H == a_ for a_ in al) and all(b_.H == b_ for b_ in be) and all(is_zero(a_ + b_) for a_, b_ in zip(al, be))
    ok = ok and all(is_zero(al[i] * al[j] + al[j] * al[i] - 2 * int(i == j) * sp.eye(4)) for i in range(3) for j in range(3))
    ok = ok and is_zero(g0 * g0 + sp.eye(4))
    check("causal", ok, "g^0 g^j Hermitian with {a_i, a_j} = 2 delta_ij (so |xi| <= |eps|^2), g^j g^0 = -g^0 g^j, "
          "g^0 g^0 = -1")
    # The structural conjugation is the Hermitian conjugate on the spinor field: |eps|^2 = sum u^2 + v^2, and the Hermitian
    # forms <eps, a_j eps> and <eps, b_j eps> are real.
    eps = spinor()
    e2 = inner(eps, eps)
    real = zero(e2 - sum(fn(f"u{a}")**2 + fn(f"v{a}")**2 for a in range(4)))
    real = real and all(zero(inner(eps, a_ * eps) - conj(inner(eps, a_ * eps))) for a_ in al + be)
    check("dagger", real, "|eps|^2 = sum_a (u_a^2 + v_a^2) and <eps, g^0 g^j eps>, <eps, g^j g^0 eps> are real")


# ---- Geometry of g = psi^4 delta ----------------------------------------------------------------------------------------
def christoffel(g, gi, crd):
    """Chr[m][n][l] = Gamma^m_nl of the metric g with typed-in inverse gi, in the coordinates crd."""
    n = len(crd)
    dg = [[[sp.diff(g[m, k], crd[l]) for l in range(n)] for k in range(n)] for m in range(n)]
    return [[[sp.expand(sum(gi[m, s] * (dg[s][l][k] + dg[s][k][l] - dg[k][l][s]) for s in range(n)) / 2)
              for l in range(n)] for k in range(n)] for m in range(n)]


def ricci(Chr, crd):
    """R_ab = d_c Gamma^c_ab - d_b Gamma^c_ac + Gamma^c_cd Gamma^d_ab - Gamma^c_bd Gamma^d_ac."""
    n = len(crd)
    sg = -1 if MUT == "ricci-sign" else 1
    return sp.Matrix(n, n, lambda a, b: sp.expand(sum(
        sp.diff(Chr[c][a][b], crd[c]) - sp.diff(Chr[c][a][c], crd[b])
        + sg * sum(Chr[c][c][e] * Chr[e][a][b] - Chr[c][b][e] * Chr[e][a][c] for e in range(n)) for c in range(n))))


class Geo:
    """The slice g = psi^4 delta: Christoffel symbols, Ricci scalar, the frame, its spin coefficients omega_jk(e_i) = g(e_j,
    nabla_{e_i} e_k) and the spinor connection matrices Om_i = (1/4) omega_jk(e_i) g^j g^k."""

    def __init__(self, psi):
        self.psi = psi
        self.g, self.gi = psi**4 * I3, (psi**-2 if MUT == "inv-typo" else psi**-4) * I3  # typed in
        self.inv_ok = is_zero(self.g * self.gi - I3)
        self.Chr = christoffel(self.g, self.gi, X3)
        Rab = ricci(self.Chr, X3)
        self.R = sp.expand(sum(self.gi[a, b] * Rab[a, b] for a in range(3) for b in range(3)))
        self.pw = 1 if MUT == "psi-frame" else 2
        E = [sp.Matrix([int(m == i) for m in range(3)]) * psi**-self.pw for i in range(3)]
        self.frame = sp.Matrix(3, 3, lambda i, j: sp.expand((E[i].T * self.g * E[j])[0]))
        sg = -1 if MUT == "christoffel" else 1

        def nabla(U, V):  # (nabla_U V)^m = U^k d_k V^m + Gamma^m_kl U^k V^l
            return sp.Matrix([sp.expand(sum(U[k] * sp.diff(V[m], X3[k]) for k in range(3))
                                        + sg * sum(self.Chr[m][k][l] * U[k] * V[l] for k in range(3) for l in range(3)))
                              for m in range(3)])
        self.om = [[[sp.expand((E[j].T * self.g * nabla(E[i], E[k]))[0]) for i in range(3)] for k in range(3)]
                   for j in range(3)]
        g0, gj = note_gammas()
        ss = -1 if MUT == "omega-sign" else 1
        self.Om = [(ss * sum((self.om[j][k][i] * gj[j] * gj[k] for j in range(3) for k in range(3)),
                             sp.zeros(4)) / 4).applyfunc(sp.expand) for i in range(3)]


def omega_hand(psi):
    """omega_jk(e_i) = psi^-2 (delta_ij phi_k - delta_ik phi_j), phi = 2 ln psi, from omega_jk = phi_k dx^j - phi_j dx^k
    for g = e^(2 phi) delta, by hand."""
    ph = [2 * sp.diff(psi, v) / psi for v in X3]
    return [[[(int(i == j) * ph[k] - int(i == k) * ph[j]) / psi**2 for i in range(3)] for k in range(3)] for j in range(3)]


def constraints(psi, K, geo):
    """mu and the frame components J_i = J(e_i) from the constraints, computed in coordinates from the frame components K_ij,
    whose coordinate components are psi^4 K_ij."""
    g, gi, Chr = geo.g, geo.gi, geo.Chr
    Kc = K if MUT == "Kframe" else psi**4 * K
    trK = sp.expand(sum(gi[a, b] * Kc[a, b] for a in range(3) for b in range(3)))
    K2 = sp.expand(sum(gi[a, c] * gi[b, d] * Kc[a, b] * Kc[c, d] for a in range(3) for b in range(3) for c in range(3)
                       for d in range(3)))
    sR = -1 if MUT == "Rsign" else 1
    sq = 1 if MUT == "Ksq" else -1
    mu = sp.expand((sR * geo.R + trK**2 + sq * K2) / (16 * sp.pi))
    tr = 1 if MUT == "trace" else -1
    Pc = (Kc + tr * g * trK).applyfunc(sp.expand)  # P_ab = K_ab - g_ab tr K

    def DP(c, a, b):  # nabla_c P_ab
        ch = 0 if MUT == "J-christoffel" else 1
        return sp.diff(Pc[a, b], X3[c]) - ch * sum(Chr[d][c][a] * Pc[d, b] + Chr[d][c][b] * Pc[a, d] for d in range(3))
    J = sp.Matrix([sp.expand(sum(gi[a, c] * DP(c, a, b) for a in range(3) for c in range(3)) / (8 * sp.pi) / psi**2)
                   for b in range(3)])
    return mu, J


# ---- Planted truth for the geometry and the constraints -------------------------------------------------------------------
def lap(f):
    return sum(sp.diff(f, v, 2) for v in X3)


def section_planted():
    print("\n== 2. The slice psi^4 delta: planted R, omega and Clifford parallelism ==")
    psi = fn("psi")
    geo = Geo(psi)
    check("inverse", geo.inv_ok, "g_ab g^bc = delta, with g = psi^4 delta and the inverse typed in")
    check("orthonormal", is_zero(geo.frame - I3), "g(e_i, e_j) = delta_ij for e_i = psi^-2 d_i")
    check("planted_R", zero(geo.R + 8 * lap(psi) / psi**5), "R = -8 psi^-5 Laplace psi, the Ricci scalar of psi^4 delta by hand")
    hand = omega_hand(psi)
    check("planted_omega", all(zero(geo.om[j][k][i] - hand[j][k][i]) for i in range(3) for j in range(3) for k in range(3)),
          "omega_jk(e_i) = g(e_j, nabla_{e_i} e_k) = psi^-2 (delta_ij phi_k - delta_ik phi_j), phi = 2 ln psi")
    # nabla_i (g^l eps) = g^l nabla_i eps + omega_jl(e_i) g^j eps, i.e. [Om_i, g^l] = omega_jl(e_i) g^j. This fixes the
    # sign and the 1/4 of the spin connection against omega itself.
    g0, gj = note_gammas()
    par = all(is_zero(geo.Om[i] * gj[l] - gj[l] * geo.Om[i] - sum((geo.om[j][l][i] * gj[j] for j in range(3)), sp.zeros(4)))
              for i in range(3) for l in range(3))
    check("clifford_parallel", par, "[(1/4) omega_jk g^j g^k, g^l] = omega_jl g^j: Clifford multiplication is parallel")
    K = sym_K()
    sK = -1 if MUT == "sign" else 1  # the stated term is mutable, the expected D~ = D + (1/2)(tr K) g^0 is not
    A = [sK * sum((K[i, j] * gj[j] * g0 for j in range(3)), sp.zeros(4)) / 2 for i in range(3)]
    check("dirac", is_zero(sum((gj[i] * A[i] for i in range(3)), sp.zeros(4)) - K.trace() * g0 / 2),
          "g^i (the stated K term)_i = (1/2)(tr K) g^0, so D~ = D + (1/2)(tr K) g^0, with the stated sign")


def section_constraints():
    print("\n== 3. Planted truth for mu and J_i from the constraints ==")
    psi = fn("psi")
    geo = Geo(psi)
    Kz = sp.zeros(3)
    # (a) Alcubierre's flat data, #58 item 1: g = delta, K_ij = (d_i b_j + d_j b_i)/2, b = (-v f, 0, 0), f(x, y, z)
    #     arbitrary. mu = -(v^2/32 pi)(f_y^2 + f_z^2) (the same as -(v^2/32 pi)((y^2+z^2)/r_s^2) f'^2 for f(r_s)) and
    #     J_i = (v/16 pi)(d_x d_i f - delta_ix Laplace f), by hand from the constraints.
    f, v = fn("f"), sp.Symbol("v", real=True)
    b = [-v * f, 0, 0]
    Ka = sp.Matrix(3, 3, lambda i, j: (sp.diff(b[j], X3[i]) + sp.diff(b[i], X3[j])) / 2)
    g1 = Geo(sp.Integer(1))
    mu, J = constraints(sp.Integer(1), Ka, g1)
    wantJ = sp.Matrix([v / (16 * sp.pi) * (sp.diff(f, X, X3[i]) - int(i == 0) * lap(f)) for i in range(3)])
    check("planted_alcubierre", zero(mu + v**2 / (32 * sp.pi) * (sp.diff(f, Y)**2 + sp.diff(f, Z)**2)) and is_zero(J - wantJ),
          "flat slice, K = Alcubierre's, arbitrary f: mu = -(v^2/32 pi)(f_y^2 + f_z^2), "
          "J_i = (v/16 pi)(d_x d_i f - delta_ix Laplace f)")
    # (b) K = 0: mu = R/16 pi = -psi^-5 Laplace psi / 2 pi, J = 0.
    mu, J = constraints(psi, Kz, geo)
    check("planted_K0", zero(mu + lap(psi) / (2 * sp.pi * psi**5)) and is_zero(J),
          "K = 0: mu = -psi^-5 Laplace psi / 2 pi, J = 0")
    # (c) K_ij = f delta_ij in the frame, so K_ab = f g_ab and P_ab = -2 f g_ab: tr K = 3f, K_ij K^ij = 3f^2,
    #     16 pi mu = R + 6 f^2, 8 pi J_a = -2 d_a f, so J_i = -psi^-2 d_i f / 4 pi in the frame.
    Kt = f * I3
    mu, J = constraints(psi, Kt, geo)
    wantJ = sp.Matrix([-sp.diff(f, X3[i]) / (4 * sp.pi * psi**2) for i in range(3)])
    check("planted_trace", zero(mu - (-lap(psi) / (2 * sp.pi * psi**5) + 3 * f**2 / (8 * sp.pi))) and is_zero(J - wantJ),
          "K_ij = f delta_ij: mu = -psi^-5 Laplace psi / 2 pi + 3 f^2 / 8 pi, J_i = -psi^-2 d_i f / 4 pi")
    # (d) K_ij = psi^-6 H_ij with H = Hessian of a harmonic polynomial, trace-free and divergence-free on flat space.
    #     K^ij = psi^-10 H_ij is the Bowen-York form, for which D_j K^ij = psi^-10 d_j H_ij = 0 for every psi. Then J = 0
    #     and 16 pi mu = R - psi^-12 H_ij H_ij.
    phi = X * Y * Z + (X**3 - 3 * X * Y**2) / 6 + sp.Rational(1, 3) * (Y**2 * Z - Z**3 / 3)
    Hh = sp.hessian(phi, X3)
    hh = sp.expand(sum(Hh[i, j]**2 for i in range(3) for j in range(3)))
    mu, J = constraints(psi, Hh / psi**6, geo)
    flat_ok = zero(lap(phi)) and all(zero(sum(sp.diff(Hh[i, j], X3[j]) for j in range(3))) for i in range(3))
    check("planted_TT", flat_ok and zero(mu - (-lap(psi) / (2 * sp.pi * psi**5) - hh / (16 * sp.pi * psi**12))) and is_zero(J),
          "K_ij = psi^-6 H_ij, H the Hessian of a harmonic polynomial: J = 0, 16 pi mu = R - psi^-12 H_ij H_ij")
    # (e) Brill-Lindquist data psi = 1 + m1/2|x| + m2/2|x - c|, K = 0 (Schwarzschild for m2 = 0): Laplace psi = 0 off the
    #     punctures, so mu = 0. The mu of (b), generic in psi, is evaluated at exact points with psi replaced, since
    #     running the constraints on an explicit psi with square roots is far slower than the points are.
    m1, m2 = sp.symbols("m1 m2", positive=True)
    r1, r2 = sp.sqrt(X**2 + Y**2 + Z**2), sp.sqrt((X - 3)**2 + Y**2 + Z**2)
    if MUT == "schwarzschild-power":
        r1, r2 = r1**2, r2**2
    ps = 1 + m1 / (2 * r1) + m2 / (2 * r2)
    mu, J = constraints(psi, Kz, geo)
    mus = mu.subs(psi, ps).doit()
    pts = [(1, 2, 2), (2, 3, 6), (-1, 5, 7), (4, -2, 1)]
    vals = [sp.simplify(mus.subs(dict(zip(X3, p)))) for p in pts]
    check("planted_schwarzschild", all(val == 0 for val in vals) and is_zero(J),
          "psi = 1 + m1/2|x| + m2/2|x - c|, K = 0: mu = 0 (vacuum) at %d exact points, for symbolic m1, m2" % len(pts))


def section_4d():
    print("\n== 4. The constraints against G_ab of the 4-metric ==")
    # g = -dt^2 + psi^4 (delta_ij - 2 t K_ij) dx^i dx^j, with N = 1 and no shift, has n = d_t and -(1/2) d_t g_ij = psi^4 K_ij
    # at t = 0, so its extrinsic curvature K_ij = -g(nabla_{e_i} n, e_j) has the frame components K_ij there. G_ab n^a n^b
    # and G_ab n^a e_i^b carry no second t derivative, so the metric to first order in t is enough, with the inverse
    # carried to second order to keep every first derivative right. Then 8 pi mu = G_tt and 8 pi J_i = -G_ti, since
    # mu = T(n, n) and J_i = -n^a T_ai, and G = 8 pi T.
    psi = fn("psi")
    K = sym_K()
    sg = -1 if MUT == "adm-K" else 1
    K2 = (K * K).applyfunc(sp.expand)
    gs = psi**4 * (I3 - sg * 2 * T * K)
    gis = psi**-4 * (I3 + sg * 2 * T * K + 4 * T**2 * K2)
    prod = (gs * gis).applyfunc(sp.expand)
    inv_ok = all(sp.expand(prod[i, j] - int(i == j)).coeff(T, p) == 0 for i in range(3) for j in range(3) for p in range(3))
    g, gi = sp.zeros(4), sp.zeros(4)
    g[0, 0], gi[0, 0] = -1, -1
    for i in range(3):
        for j in range(3):
            g[i + 1, j + 1], gi[i + 1, j + 1] = gs[i, j], gis[i, j]
    Chr = christoffel(g, gi, CRD)
    Chr = [[[sp.expand(c.subs(T, 0) + T * sp.diff(c, T).subs(T, 0)) for c in row] for row in blk] for blk in Chr]  # O(t)
    Rab = ricci(Chr, CRD).subs(T, 0)
    gi0 = gi.subs(T, 0)
    Rs = sp.expand(sum(gi0[a, b] * Rab[a, b] for a in range(4) for b in range(4)))
    Gtt = sp.expand(Rab[0, 0] + Rs / 2)  # G_tt = R_tt - (1/2) g_tt R, g_tt = -1
    sj = 1 if MUT == "J4d" else -1
    mu4 = Gtt / (8 * sp.pi)
    J4 = sp.Matrix([sj * Rab[0, i + 1] / (8 * sp.pi) / psi**2 for i in range(3)])  # J_i = -G_ti / 8 pi, frame component
    mu, J = constraints(psi, K, Geo(psi))
    check("constraints_4d", inv_ok and zero(mu - mu4) and is_zero(J - J4),
          "the constraints' mu and J_i equal G_tt / 8 pi and -G_ti / 8 pi of -dt^2 + psi^4 (delta - 2 t K) dx dx at t = 0, "
          "for arbitrary psi and K_ij")


# ---- The identity ---------------------------------------------------------------------------------------------------------
def identity(name, psi, K, text, show=False, reduced=False):
    geo = Geo(psi)
    mu, J = constraints(psi, K, geo)
    g0, gj = note_gammas()
    eps = spinor()
    trK = K.trace()
    P = K - trK * I3
    sK = -1 if MUT == "sign" else 1
    A = [sK * sum((K[i, j] * gj[j] * g0 for j in range(3)), sp.zeros(4)) / 2 for i in range(3)]
    nab = [sp.expand(psi**-geo.pw * eps.diff(X3[i]) + geo.Om[i] * eps) for i in range(3)]  # nabla_i eps
    nabt = [sp.expand(nab[i] + A[i] * eps) for i in range(3)]  # nabla~_i eps
    Deps = sum((gj[i] * nab[i] for i in range(3)), sp.zeros(4, 1))
    Dteps = sum((gj[i] * nabt[i] for i in range(3)), sp.zeros(4, 1))  # D~ eps = g^i nabla~_i eps
    e2 = inner(eps, eps)
    xi = [inner(eps, g0 * gj[j] * eps) for j in range(3)]
    Jxi = sum(J[j] * xi[j] for j in range(3))
    sM = -1 if MUT == "musign" else 1
    sJ = 1 if MUT == "Jsign" else -1
    lhs = sum(inner(nabt[i], nabt[i]) for i in range(3)) - inner(Dteps, Dteps) + 4 * sp.pi * (sM * mu * e2 + sJ * Jxi)
    lhs = sp.expand(lhs)
    # V^i = Re<eps, nabla^i eps - g^i D eps> - (1/2) P_ij xi^j, the components of V = V^i e_i
    V = []
    for i in range(3):
        spin = nab[i] if MUT == "Vspin" else nab[i] - gj[i] * Deps
        Pxi = 0 if MUT == "Kdiv" else sum(P[i, j] * xi[j] for j in range(3)) / 2
        V.append(sp.expand(re(inner(eps, spin)) - Pxi))
    if MUT == "weight":
        div = sum(sp.diff(V[a], X3[a]) for a in range(3))
    else:
        div = sum(sp.diff(psi**4 * V[a], X3[a]) for a in range(3)) / psi**6
    div = sp.expand(div)
    res = sp.expand(lhs - div)
    check(name, res == 0, f"{text}; {nterms(lhs)} terms on the left, {nterms(div)} in div V, "
          f"{len(sp.Add.make_args(res)) if res != 0 else 0} in the difference")
    if show:
        for i in range(3):
            print(f"    V^{'xyz'[i]} = {V[i]}")
    if reduced:
        # The j = i term of g^i D eps is g^i g^i nabla_i eps = nabla_i eps, so the same V is
        # V^i = -sum_{j != i} Re<eps, g^i g^j nabla_j eps> - (1/2) P_ij xi^j, written out here without the cancellation.
        Vr = [sp.expand(-sum(re(inner(eps, gj[i] * gj[j] * nab[j])) for j in range(3) if j != i)
                        - sum(P[i, j] * xi[j] for j in range(3)) / 2) for i in range(3)]
        check("V_reduced", all(zero(V[i] - Vr[i]) for i in range(3)),
              "V^i = -sum_{j != i} Re<eps, g^i g^j nabla_j eps> - (1/2) P_ij xi^j, the j = i term having cancelled")


def section_identity():
    print("\n== 5. The identity |nabla~ eps|^2 - |D~ eps|^2 + 4 pi (mu |eps|^2 - J.xi) = div V ==")
    print("    V^i = Re<eps, nabla^i eps - g^i D eps> - (1/2) P_ij xi^j, P_ij = K_ij - delta_ij tr K, "
          "div V = psi^-6 d_a(psi^4 V^a)")
    psi = fn("psi")
    Kg = sym_K()
    one, Kz = sp.Integer(1), sp.zeros(3)
    f, v = fn("f"), sp.Symbol("v", real=True)
    b = [-v * f, 0, 0]
    Ka = sp.Matrix(3, 3, lambda i, j: (sp.diff(b[j], X3[i]) + sp.diff(b[i], X3[j])) / 2)
    print("  Planted first: the flat Weitzenbock identity, with V written out for the arbitrary spinor")
    identity("identity_flat_K0", one, Kz, "psi = 1, K = 0, so mu = J = 0 and V^i = -sum_{j != i} Re<eps, g^i g^j d_j eps>",
             show=True)
    identity("identity_lichnerowicz", psi, Kz, "psi arbitrary, K = 0, mu = R/16 pi, J = 0: Lichnerowicz")
    identity("identity_flat_K", one, Kg, "psi = 1, K_ij arbitrary: mu = ((tr K)^2 - K_ij K_ij)/16 pi, "
             "J_i = (d_j K_ij - d_i tr K)/8 pi")
    identity("identity_alcubierre", one, Ka, "psi = 1, K = Alcubierre's with b = (-v f, 0, 0), f(x, y, z) arbitrary: "
             "#59's weitzenbock, pointwise")
    identity("identity_general", psi, Kg, "psi, K_ij and eps arbitrary", reduced=True)


# ---- The diagonal slice g = diag(A^2, B^2, C^2), issue #99 -----------------------------------------------------------------
# The same checks as above for g = diag(A^2, B^2, C^2)(x, y, z) with A, B, C arbitrary positive functions, the frame
# e_i = d_i / A_i, and div V = (1/ABC) d_a(ABC V^a / A_a). GeoDiag, constraints_diag and identity_diag are written beside Geo,
# constraints and identity and call none of them, so the conformally flat results above are as #87 left them, and
# A = B = C = psi^2 is a planted case for the new code (diag_conf_limit, identity_diag_conf).
def todo(name):
    """Whether the diagonal check `name` runs: not when staged away with --only or --skip, and under a mutation only if the
    mutation targets it, since the general case is expensive and the other checks of a mutated run need not repeat it."""
    if ONLY is not None and name not in ONLY:
        return False
    return name not in SKIP and (not MUT or name in CONTROLS[MUT][1])


class GeoDiag:
    """The slice g = diag(A^2, B^2, C^2): Christoffel symbols, Ricci scalar, the frame e_i = d_i / A_i, its spin coefficients
    omega_jk(e_i) = g(e_j, nabla_{e_i} e_k) and the spinor connection matrices Om_i = (1/4) omega_jk(e_i) g^j g^k. The
    spin connection comes from the Christoffel symbols, as Geo's does."""

    def __init__(self, A):
        self.A = A = list(A)
        self.g = sp.diag(*[a**2 for a in A])
        self.gi = sp.diag(*[a**(-1 if MUT == "inv-typo" else -2) for a in A])  # typed in
        self.inv_ok = is_zero(self.g * self.gi - I3)
        self.Chr = christoffel(self.g, self.gi, X3)
        Rab = ricci(self.Chr, X3)
        self.R = sp.expand(sum(self.gi[a, b] * Rab[a, b] for a in range(3) for b in range(3)))
        self.sc = [a**(-2 if MUT == "diag-frame" else -1) for a in A]  # e_i = sc_i d_i
        if MUT == "aniso-frame":
            self.sc[0] = 1 / A[1]
        E = [sp.Matrix([int(m == i) for m in range(3)]) * self.sc[i] for i in range(3)]
        self.frame = sp.Matrix(3, 3, lambda i, j: sp.expand((E[i].T * self.g * E[j])[0]))
        sg = -1 if MUT == "christoffel" else 1

        def nabla(U, V):  # (nabla_U V)^m = U^k d_k V^m + Gamma^m_kl U^k V^l
            return sp.Matrix([sp.expand(sum(U[k] * sp.diff(V[m], X3[k]) for k in range(3))
                                        + sg * sum(self.Chr[m][k][l] * U[k] * V[l] for k in range(3) for l in range(3)))
                              for m in range(3)])
        self.om = [[[sp.expand((E[j].T * self.g * nabla(E[i], E[k]))[0]) for i in range(3)] for k in range(3)]
                   for j in range(3)]
        g0, gj = note_gammas()
        ss = -1 if MUT == "omega-sign" else 1
        self.Om = [(ss * sum((self.om[j][k][i] * gj[j] * gj[k] for j in range(3) for k in range(3)),
                             sp.zeros(4)) / 4).applyfunc(sp.expand) for i in range(3)]


def omega_hand_diag(A):
    """omega_jk(e_i) = (delta_ij d_k A_j - delta_ik d_j A_k) / (A_j A_k) for j != k and 0 for j = k, from Cartan's first
    structure equation with theta^j = A_j dx^j, as in conventions.py section 6(c) (#86)."""
    s2 = -1 if MUT == "omega-hand" else 1
    return [[[0 if j == k else (int(i == j) * sp.diff(A[j], X3[k]) - s2 * int(i == k) * sp.diff(A[k], X3[j]))
              / (A[j] * A[k]) for i in range(3)] for k in range(3)] for j in range(3)]


def constraints_diag(A, K, geo):
    """mu and the frame components J_i = J(e_i) from the constraints, computed in coordinates from the frame components K_ij,
    whose coordinate components are A_i A_j K_ij."""
    g, gi, Chr = geo.g, geo.gi, geo.Chr
    Kc = sp.Matrix(3, 3, lambda a, b: K[a, b] if MUT == "Kframe" else A[a] * A[b] * K[a, b])
    trK = sp.expand(sum(gi[a, b] * Kc[a, b] for a in range(3) for b in range(3)))
    K2 = sp.expand(sum(gi[a, c] * gi[b, d] * Kc[a, b] * Kc[c, d] for a in range(3) for b in range(3) for c in range(3)
                       for d in range(3)))
    sR = -1 if MUT == "Rsign" else 1
    sq = 1 if MUT == "Ksq" else -1
    mu = sp.expand((sR * geo.R + trK**2 + sq * K2) / (16 * sp.pi))
    tr = 1 if MUT == "trace" else -1
    Pc = (Kc + tr * g * trK).applyfunc(sp.expand)  # P_ab = K_ab - g_ab tr K

    def DP(c, a, b):  # nabla_c P_ab
        ch = 0 if MUT == "J-christoffel" else 1
        return sp.diff(Pc[a, b], X3[c]) - ch * sum(Chr[d][c][a] * Pc[d, b] + Chr[d][c][b] * Pc[a, d] for d in range(3))
    J = sp.Matrix([sp.expand(sum(gi[a, c] * DP(c, a, b) for a in range(3) for c in range(3)) / (8 * sp.pi) / A[b])
                   for b in range(3)])
    return mu, J


def section_planted_diag():
    names = ["inverse_diag", "orthonormal_diag", "planted_R_diag", "planted_omega_diag", "clifford_parallel_diag",
             "diag_conf_limit"]
    if not any(todo(n) for n in names):
        return
    print("\n== 6. The diagonal slice diag(A^2, B^2, C^2): planted R, omega and Clifford parallelism ==")
    A = [fn(n) for n in "ABC"]
    geo = GeoDiag(A)
    if todo("inverse_diag"):
        check("inverse_diag", geo.inv_ok, "g_ab g^bc = delta, with g = diag(A^2, B^2, C^2) and the inverse typed in")
    if todo("orthonormal_diag"):
        check("orthonormal_diag", is_zero(geo.frame - I3), "g(e_i, e_j) = delta_ij for e_i = d_i / A_i")
    if todo("planted_R_diag"):
        # By hand: for a metric h + c^2 dz^2 on a base h, R = R_h - 2 Laplace_h(c) / c, which gives (b), and the Gaussian
        # curvature of dx^2 + f^2 dy^2 is -f''/f, which with (b) gives (d). (c) is (d) with g = f, which is hyperbolic space,
        # R = -6, for f = exp(x). (e) is (d) in the proper length dr = A dx.
        f, h, c2 = sp.Function("f")(X), sp.Function("g")(X), sp.Function("c")(X, Y)
        a_, b_, c_ = [sp.Function(n)(X) for n in "abc"]
        one = sp.Integer(1)
        sf = -1 if MUT == "R-hand" else 1
        cases = [
            ("flat", [one, one, one], sp.Integer(0)),
            ("dx^2 + dy^2 + c(x,y)^2 dz^2", [one, one, c2], -2 * (sp.diff(c2, X, 2) + sp.diff(c2, Y, 2)) / c2),
            ("dx^2 + f(x)^2 (dy^2 + dz^2)", [one, f, f], -4 * sp.diff(f, X, 2) / f - sf * 2 * sp.diff(f, X)**2 / f**2),
            ("dx^2 + f(x)^2 dy^2 + g(x)^2 dz^2", [one, f, h],
             -2 * (sp.diff(f, X, 2) / f + sp.diff(h, X, 2) / h + sp.diff(f, X) * sp.diff(h, X) / (f * h))),
            ("A(x)^2 dx^2 + B(x)^2 dy^2 + C(x)^2 dz^2", [a_, b_, c_],
             -2 * (sp.diff(sp.diff(b_, X) / a_, X) / (a_ * b_) + sp.diff(sp.diff(c_, X) / a_, X) / (a_ * c_)
                   + sp.diff(b_, X) * sp.diff(c_, X) / (a_**2 * b_ * c_))),
        ]
        ok = True
        for _, Ac, want in cases:
            ok = ok and zero(GeoDiag(Ac).R - want)
        psi = fn("psi")
        ok = ok and zero(GeoDiag([psi**2] * 3).R + 8 * lap(psi) / psi**5)
        check("planted_R_diag", ok, "R of %d metrics with a hand value: flat, dx^2 + dy^2 + c(x,y)^2 dz^2 (R = -2 Laplace c / c), "
              "dx^2 + f^2 (dy^2 + dz^2) (R = -4 f''/f - 2 f'^2/f^2), dx^2 + f^2 dy^2 + g^2 dz^2, A(x), B(x), C(x) arbitrary, "
              "and psi^2 for A, B, C (R = -8 psi^-5 Laplace psi)" % (len(cases) + 1))
    if todo("planted_omega_diag"):
        hand = omega_hand_diag(A)
        check("planted_omega_diag", all(zero(geo.om[j][k][i] - hand[j][k][i]) for i in range(3) for j in range(3)
                                        for k in range(3)),
              "omega_jk(e_i) = g(e_j, nabla_{e_i} e_k) = (delta_ij d_k A_j - delta_ik d_j A_k)/(A_j A_k), j != k, 0 for j = k, "
              "Cartan's formula of conventions.py section 6(c)")
    if todo("clifford_parallel_diag"):
        g0, gj = note_gammas()
        par = all(is_zero(geo.Om[i] * gj[l] - gj[l] * geo.Om[i] - sum((geo.om[j][l][i] * gj[j] for j in range(3)),
                                                                      sp.zeros(4)))
                  for i in range(3) for l in range(3))
        check("clifford_parallel_diag", par, "[(1/4) omega_jk g^j g^k, g^l] = omega_jl g^j on the diagonal slice")
    if todo("diag_conf_limit"):
        psi = fn("psi")
        gd, gc = GeoDiag([psi**2] * 3), Geo(psi)
        Kg = sym_K()
        mud, Jd = constraints_diag([psi**2] * 3, Kg, gd)
        muc, Jc = constraints(psi, Kg, gc)
        same = is_zero(gd.g - gc.g) and zero(gd.R - gc.R) and zero(mud - muc) and is_zero(Jd - Jc)
        same = same and all(zero(gd.om[j][k][i] - gc.om[j][k][i]) for i in range(3) for j in range(3) for k in range(3))
        same = same and all(is_zero(gd.Om[i] - gc.Om[i]) for i in range(3))
        check("diag_conf_limit", same, "A = B = C = psi^2 gives #87's g, R, omega_jk, Om_i and, for K_ij arbitrary, mu and J_i")


def section_constraints_diag():
    names = ["planted_alcubierre_diag", "planted_K0_diag", "planted_trace_diag", "planted_TT_diag", "planted_kasner_diag"]
    if not any(todo(n) for n in names):
        return
    print("\n== 7. Planted truth for mu and J_i from the constraints, on diagonal slices ==")
    one = sp.Integer(1)
    f, v = fn("f"), sp.Symbol("v", real=True)
    c2 = sp.Function("c")(X, Y)
    a_, w_ = sp.Function("a")(X), sp.Function("w")(X)
    # Two metrics with a hand value of R. (i) dx^2 + dy^2 + c(x, y)^2 dz^2, R = -2 Laplace_xy c / c, whose Gamma Gamma terms
    # cancel in R and whose g_bb does not depend on x_b. (ii) a(x)^2 dx^2 + w(x)^2 (dy^2 + dz^2), which is the metric of
    # planted_R_diag's third case in the proper length dr = a dx, R = -4 (1/a)(w'/a)'/w - 2 w'^2/(a w)^2, with Gamma Gamma terms
    # in R and d_x g_xx != 0, so that the Christoffel terms of the divergence in J matter.
    metrics = [([one, one, c2], -2 * (sp.diff(c2, X, 2) + sp.diff(c2, Y, 2)) / c2),
               ([a_, w_, w_], -4 * sp.diff(sp.diff(w_, X) / a_, X) / (a_ * w_) - 2 * sp.diff(w_, X)**2 / (a_ * w_)**2)]
    if todo("planted_alcubierre_diag"):
        # The planted flat case: A = B = C = 1, through constraints_diag, with #58 item 1's values.
        b = [-v * f, 0, 0]
        Ka = sp.Matrix(3, 3, lambda i, j: (sp.diff(b[j], X3[i]) + sp.diff(b[i], X3[j])) / 2)
        mu, J = constraints_diag([one] * 3, Ka, GeoDiag([one] * 3))
        wantJ = sp.Matrix([v / (16 * sp.pi) * (sp.diff(f, X, X3[i]) - int(i == 0) * lap(f)) for i in range(3)])
        check("planted_alcubierre_diag", zero(mu + v**2 / (32 * sp.pi) * (sp.diff(f, Y)**2 + sp.diff(f, Z)**2))
              and is_zero(J - wantJ),
              "A = B = C = 1, K = Alcubierre's, arbitrary f: mu = -(v^2/32 pi)(f_y^2 + f_z^2), "
              "J_i = (v/16 pi)(d_x d_i f - delta_ix Laplace f)")
    if todo("planted_K0_diag"):
        ok = True
        for Ac, Rc in metrics:
            mu, J = constraints_diag(Ac, sp.zeros(3), GeoDiag(Ac))
            ok = ok and zero(mu - Rc / (16 * sp.pi)) and is_zero(J)
        check("planted_K0_diag", ok, "K = 0 on the two metrics above: mu = R/16 pi with R by hand, J = 0")
    if todo("planted_trace_diag"):
        # K_ij = f delta_ij in the frame, so K_ab = f g_ab, P_ab = -2 f g_ab, 16 pi mu = R + 6 f^2, 8 pi J_a = -2 d_a f.
        ok = True
        for Ac, Rc in metrics:
            mu, J = constraints_diag(Ac, f * I3, GeoDiag(Ac))
            wantJ = sp.Matrix([-sp.diff(f, X3[i]) / (4 * sp.pi * Ac[i]) for i in range(3)])
            ok = ok and zero(mu - (Rc / (16 * sp.pi) + 3 * f**2 / (8 * sp.pi))) and is_zero(J - wantJ)
        check("planted_trace_diag", ok, "K_ij = f delta_ij on the two metrics above: mu = R/16 pi + 3 f^2/8 pi, "
              "J_i = -d_i f/(4 pi A_i)")
    if todo("planted_TT_diag"):
        # The Bowen-York form of #87's planted_TT, with A = B = C = psi^2 given to the diagonal code.
        psi = fn("psi")
        phi = X * Y * Z + (X**3 - 3 * X * Y**2) / 6 + sp.Rational(1, 3) * (Y**2 * Z - Z**3 / 3)
        Hh = sp.hessian(phi, X3)
        hh = sp.expand(sum(Hh[i, j]**2 for i in range(3) for j in range(3)))
        mu, J = constraints_diag([psi**2] * 3, Hh / psi**6, GeoDiag([psi**2] * 3))
        check("planted_TT_diag", zero(mu - (-lap(psi) / (2 * sp.pi * psi**5) - hh / (16 * sp.pi * psi**12))) and is_zero(J),
              "A = B = C = psi^2, K_ij = psi^-6 H_ij, H the Hessian of a harmonic polynomial: J = 0, "
              "16 pi mu = R - psi^-12 H_ij H_ij")
    if todo("planted_kasner_diag"):
        # A Kasner slice at t = 1: constant scale factors a_i (so g is flat), K_ij = diag(p_i) in the frame with sum p = sum p^2 = 1,
        # which is vacuum, so (tr K)^2 = K_ij K_ij and mu = J = 0 for every a_i.
        a = [sp.Symbol(f"a{i + 1}", positive=True) for i in range(3)]
        p = [sp.Rational(1, 2)] * 3 if MUT == "kasner-exponents" else [sp.Rational(-2, 7), sp.Rational(3, 7),
                                                                       sp.Rational(6, 7)]
        mu, J = constraints_diag(a, sp.diag(*p), GeoDiag(a))
        check("planted_kasner_diag", zero(mu) and is_zero(J),
              "constant A_i = a_i, K_ij = diag(-2/7, 3/7, 6/7) (Kasner at t = 1: sum p = sum p^2 = 1): mu = J = 0")


def section_4d_diag():
    if not todo("constraints_4d_diag"):
        return
    print("\n== 8. The constraints against G_ab of the 4-metric, diagonal slice ==")
    # g = -dt^2 + Lam (delta - 2 t K) Lam, Lam = diag(A, B, C), with N = 1 and no shift: the coordinate components of the
    # extrinsic curvature K(d_i, d_j) = A_i A_j K_ij at t = 0 are -(1/2) d_t g_ij, as in section 4.
    A = [fn(n) for n in "ABC"]
    K = sym_K()
    sg = -1 if MUT == "adm-K" else 1
    K2 = (K * K).applyfunc(sp.expand)
    Lam, Lami = sp.diag(*A), sp.diag(*[1 / a for a in A])
    gs = (Lam * (I3 - sg * 2 * T * K) * Lam).applyfunc(sp.expand)
    gis = (Lami * (I3 + sg * 2 * T * K + 4 * T**2 * K2) * Lami).applyfunc(sp.expand)
    prod = (gs * gis).applyfunc(sp.expand)
    inv_ok = all(sp.expand(prod[i, j] - int(i == j)).coeff(T, p) == 0 for i in range(3) for j in range(3) for p in range(3))
    g, gi = sp.zeros(4), sp.zeros(4)
    g[0, 0], gi[0, 0] = -1, -1
    for i in range(3):
        for j in range(3):
            g[i + 1, j + 1], gi[i + 1, j + 1] = gs[i, j], gis[i, j]
    Chr = christoffel(g, gi, CRD)
    Chr = [[[sp.expand(c.subs(T, 0) + T * sp.diff(c, T).subs(T, 0)) for c in row] for row in blk] for blk in Chr]  # O(t)
    Rab = ricci(Chr, CRD).subs(T, 0)
    gi0 = gi.subs(T, 0)
    Rs = sp.expand(sum(gi0[a, b] * Rab[a, b] for a in range(4) for b in range(4)))
    Gtt = sp.expand(Rab[0, 0] + Rs / 2)
    sj = 1 if MUT == "J4d" else -1
    mu4 = Gtt / (8 * sp.pi)
    J4 = sp.Matrix([sj * Rab[0, i + 1] / (8 * sp.pi) / A[i] for i in range(3)])  # J_i = -G_ti / 8 pi, frame component
    mu, J = constraints_diag(A, K, GeoDiag(A))
    check("constraints_4d_diag", inv_ok and zero(mu - mu4) and is_zero(J - J4),
          "the constraints' mu and J_i equal G_tt / 8 pi and -G_ti / 8 pi of -dt^2 + Lam (delta - 2 t K) Lam at t = 0, "
          "Lam = diag(A, B, C), for arbitrary A, B, C and K_ij")


def identity_diag(name, A, K, text, reduced=False):
    run_main, run_red = todo(name), reduced and todo("V_reduced_diag")
    if not (run_main or run_red):
        return
    t0 = time.time()
    geo = GeoDiag(A)
    mu, J = constraints_diag(A, K, geo)
    g0, gj = note_gammas()
    eps = spinor()
    trK = K.trace()
    P = K - trK * I3
    sK = -1 if MUT == "sign" else 1
    Aterm = [sK * sum((K[i, j] * gj[j] * g0 for j in range(3)), sp.zeros(4)) / 2 for i in range(3)]
    nab = [sp.expand(geo.sc[i] * eps.diff(X3[i]) + geo.Om[i] * eps) for i in range(3)]  # nabla_i eps
    nabt = [sp.expand(nab[i] + Aterm[i] * eps) for i in range(3)]  # nabla~_i eps
    Deps = sum((gj[i] * nab[i] for i in range(3)), sp.zeros(4, 1))
    Dteps = sum((gj[i] * nabt[i] for i in range(3)), sp.zeros(4, 1))  # D~ eps = g^i nabla~_i eps
    e2 = inner(eps, eps)
    xi = [inner(eps, g0 * gj[j] * eps) for j in range(3)]
    Jxi = sum(J[j] * xi[j] for j in range(3))
    sM = -1 if MUT == "musign" else 1
    sJ = 1 if MUT == "Jsign" else -1
    V = []
    for i in range(3):
        spin = nab[i] if MUT == "Vspin" else nab[i] - gj[i] * Deps
        Pxi = 0 if MUT == "Kdiv" else sum(P[i, j] * xi[j] for j in range(3)) / 2
        V.append(sp.expand(re(inner(eps, spin)) - Pxi))
    if run_main:
        lhs = sp.expand(sum(inner(nabt[i], nabt[i]) for i in range(3)) - inner(Dteps, Dteps)
                        + 4 * sp.pi * (sM * mu * e2 + sJ * Jxi))
        sqg = A[0]**3 if MUT == "abc" else A[0] * A[1] * A[2]  # sqrt(det g) = A B C
        if MUT == "weight":
            div = sum(sp.diff(V[a], X3[a]) for a in range(3))
        else:  # div V = (1/ABC) d_a(ABC V^a / A_a), the coordinate components of V = V^a e_a being V^a sc_a
            div = sum(sp.diff(sqg * geo.sc[a] * V[a], X3[a]) for a in range(3)) / sqg
        div = sp.expand(div)
        res = sp.expand(lhs - div)
        check(name, res == 0, f"{text}; {nterms(lhs)} terms on the left, {nterms(div)} in div V, "
              f"{len(sp.Add.make_args(res)) if res != 0 else 0} in the difference")
        print(f"    ({name}: {time.time() - t0:.0f} s)", file=sys.stderr, flush=True)
    if run_red:
        Vr = [sp.expand(-sum(re(inner(eps, gj[i] * gj[j] * nab[j])) for j in range(3) if j != i)
                        - sum(P[i, j] * xi[j] for j in range(3)) / 2) for i in range(3)]
        check("V_reduced_diag", all(zero(V[i] - Vr[i]) for i in range(3)),
              "V^i = -sum_{j != i} Re<eps, g^i g^j nabla_j eps> - (1/2) P_ij xi^j, the j = i term having cancelled")


def section_identity_diag():
    names = ["identity_diag_flat_K0", "identity_diag_lichnerowicz", "identity_diag_flat_K", "identity_diag_alcubierre",
             "identity_diag_conf", "identity_diag_Kdiag", "identity_diag_general"]
    if not any(todo(n) for n in names):
        return
    print("\n== 9. The identity |nabla~ eps|^2 - |D~ eps|^2 + 4 pi (mu |eps|^2 - J.xi) = div V on diag(A^2, B^2, C^2) ==")
    print("    V^i as in section 5, div V = (1/ABC) d_a(ABC V^a / A_a)")
    one, Kz = sp.Integer(1), sp.zeros(3)
    A = [fn(n) for n in "ABC"]
    Kg = sym_K()
    f, v = fn("f"), sp.Symbol("v", real=True)
    b = [-v * f, 0, 0]
    Ka = sp.Matrix(3, 3, lambda i, j: (sp.diff(b[j], X3[i]) + sp.diff(b[i], X3[j])) / 2)
    Kd = sp.diag(*[fn(f"k{i + 1}") for i in range(3)])
    psi = fn("psi")
    print("  Planted first: A = B = C = 1, then K = 0, through the code of the general case")
    identity_diag("identity_diag_flat_K0", [one] * 3, Kz, "A = B = C = 1, K = 0, so mu = J = 0: the flat Weitzenbock identity")
    identity_diag("identity_diag_lichnerowicz", A, Kz, "A, B, C arbitrary, K = 0, mu = R/16 pi, J = 0: Lichnerowicz")
    identity_diag("identity_diag_flat_K", [one] * 3, Kg, "A = B = C = 1, K_ij arbitrary")
    identity_diag("identity_diag_alcubierre", [one] * 3, Ka, "A = B = C = 1, K = Alcubierre's with b = (-v f, 0, 0), "
                  "f(x, y, z) arbitrary")
    identity_diag("identity_diag_conf", [psi**2] * 3, Kg, "A = B = C = psi^2, K_ij and eps arbitrary: #87's "
                  "identity_general through the diagonal code")
    identity_diag("identity_diag_Kdiag", A, Kd, "A, B, C arbitrary, K_ij = diag(k1, k2, k3) arbitrary functions, eps arbitrary")
    identity_diag("identity_diag_general", A, Kg, "A, B, C, K_ij and eps arbitrary", reduced=True)


def main():
    t0 = time.time()
    print("weitzenbock_conf.py, issues #87 and #99. sympy %s, python %s. mutation: %s" % (sp.__version__, sys.version.split()[0],
                                                                                MUT or "none"))
    section_gammas()
    section_planted()
    section_constraints()
    section_4d()
    section_identity()
    section_planted_diag()
    section_constraints_diag()
    section_4d_diag()
    section_identity_diag()
    fails = [n for n, ok in RESULTS if not ok]
    print("\nchecks: %d passed, %d failed %s. %.0f s" % (len(RESULTS) - len(fails), len(fails), fails, time.time() - t0))
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

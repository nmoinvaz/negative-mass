#!/usr/bin/env python3
"""Item 6's Weitzenbock identity on a conformally flat slice (#87), on a diagonal slice (#99) and in a frame rotated about z
on it (#130), checked symbolically.

Data, #87: g = psi^4 delta on R^3 with psi(x, y, z) arbitrary, an arbitrary symmetric K_ij(x, y, z) in the orthonormal
frame e_i = psi^-2 d_i, and an arbitrary spinor field eps = u + i v with the eight real components u_a, v_a of (x, y, z).
Data, #99: g = diag(A^2, B^2, C^2)(x, y, z) with A, B, C arbitrary positive functions, e_i = d_i / A_i, K_ij = K(e_i, e_j)
arbitrary and symmetric, eps arbitrary, and div V = (1/ABC) d_a(ABC V^a / A_a) for V = V^a e_a. A = B = C = psi^2 is #87's
slice. Data, #130: the same slice in the frame e'_1 = c e_1 + s e_2, e'_2 = -s e_1 + c e_2, e'_3 = e_3, with c =
cos(alpha), s = sin(alpha) and alpha(x, y, z) arbitrary, K'_ij = K(e'_i, e'_j) arbitrary and symmetric, eps arbitrary, V =
V'^a e'_a and div V = (1/ABC) d_m(ABC V^m) with the coordinate components V^m = V'^a e'^m_a. The note's gammas
(mass_bounds.md, "Conventions") act on eps, <a, b> = a^dagger b, xi^j = <eps, g^0 g^j eps>, and
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

In the rotated frame (sections 10 to 13) the planted values are: omega'_jk(e'_i) = R_ip R_jl R_km omega_lm(e_p) +
(delta_j2 delta_k1 - delta_j1 delta_k2) e'_i(alpha) with Cartan's omega_lm of the diagonal frame and R the rotation, which
is the closed form of #86 (omega'_21(X) = X(alpha) on a flat frame), and #86's own slice connection from conventions.py
(`omega_rot_conv`); orthonormality of e'; Clifford parallelism; mu and J'_i = J(e'_i) by hand for K = 0, K' = f delta (the
same in every frame, J' = -e'_i(f)/4 pi), Alcubierre's K rotated, K' = R K R^T, and a Kasner slice rotated, K' = R diag(p)
R^T; mu and J' against #99's constraints_diag for K = R^T K' R; and against G_ab of the 4-metric with that K. Then alpha =
0 is compared with identity_diag_general term for term, and the identities run flat (K = 0) and Lichnerowicz with alpha
arbitrary, flat with K' arbitrary, Alcubierre's K rotated, K' = diag(k1, k2, k3), and last alpha, A, B, C, K'_ij and eps
all arbitrary.

The reduction. In the rotated frame every expression is a polynomial in c = cos(alpha) and s = sin(alpha) whose
coefficients are Laurent polynomials in the jets of A, B, C, alpha, K'_ij, u_a and v_a (sympy differentiates cos(alpha) by
the chain rule, so d_i c = -s d_i alpha). `red` expands and, as conventions.py's trig_reduce (#86) does, divides by c^2 +
s^2 - 1: c^k goes to c^(k mod 2) (1 - s^2)^(k div 2). That polynomial is monic in c, so the remainder is unique and has
degree at most 1 in c, and the remainder of e is 0 exactly when e is a multiple of c^2 + s^2 - 1. That is the right test
for "zero for every alpha": the jets of alpha at a point are independent of alpha itself, so the only relation among c, s
and the jets is c^2 + s^2 = 1, and a polynomial vanishing on the real circle vanishes on the irreducible conic c^2 + s^2 =
1 and is a multiple of its equation. Reduced results are mapped back to cos(alpha) and sin(alpha) before they are
differentiated again.

What the checks share. The gammas are retyped from conventions.py (compared with its note_gammas() in `gammas_match`), and
christoffel and ricci are retyped from it too. Nothing else is imported from it, except in `omega_rot_conv`, which calls
its slice_data. The identity side (matrices, the frame, V) shares with the constraint side only sympy, christoffel and
ricci. The spin connection of the identity comes from the Christoffel symbols of the metric and is compared with the
closed form above. The constraints are computed in coordinates and compared with the 4-metric's G_ab, which shares
christoffel and ricci with them, so a wrong ricci is seen only by the planted values of R. The diagonal code (GeoDiag,
constraints_diag, identity_diag) is written beside #87's Geo, constraints and identity and calls none of them, so #87's
checks run as they were. They share the helpers (christoffel, ricci, inner, re, zero, sym_K, spinor, note_gammas), and
`diag_conf_limit` and `identity_diag_conf` compare the two codes at A = B = C = psi^2. The rotated code (GeoRot,
constraints_rot, rot_terms) is written beside the diagonal code in the same way: GeoRot subclasses GeoDiag, so g, the
inverse, the Christoffel symbols and R are GeoDiag's, but its frame, spin connection, constraints and identity are its
own, and `diag_terms`, which identity_diag now calls, is the diagonal code that `rot_alpha0_limit` and
`identity_rot_alpha0` compare it with at alpha = 0. The spin connection of the rotated frame comes from the Christoffel
symbols and is compared with the closed form, which shares only Cartan's omega_lm of the diagonal frame, and with #86's
slice_data, which shares sympy, the retyped christoffel and the retyped trig_reduce. The 4-metric route builds K = R^T K'
R from the hand matrix R and not from GeoRot's coframe. Complex conjugation is the structural replacement i -> -i on
expressions in the real functions u_a, v_a, psi, A, B, C and K_ij. Zero is tested by sp.expand, which is canonical for
Laurent polynomials in the jets, and sp.simplify only on the planted cases with sqrt, and in the rotated frame by `red`,
which is canonical modulo c^2 + s^2 = 1.

What the code proves (*checked symbolically*): at every point, for all psi, K_ij and eps on psi^4 delta, and for all A, B,
C, K_ij and eps on diag(A^2, B^2, C^2)(x, y, z), and for all alpha, A, B, C, K'_ij and eps in the frame e'_a rotated by
alpha(x, y, z) about z on that slice, the two sides above are equal as polynomials in the jets, to second order, of the
functions (in the rotated frame also modulo c^2 + s^2 = 1), in the note's gammas with the stated sign of nabla~ and with
mu and J_i the constraints' values. It does not prove the integrated identity, whose boundary term is div V and needs
decay (review.md gap 2), or that D~ eps = 0 has a solution (item 5). The hand derivation above uses only the Levi-Civita
connection, the Clifford relations and Hermiticity, so it holds on any slice (*proved* by hand). Still covered only by it:
slices that are not diagonal in the coordinates, and frames other than e_i = d_i / A_i and its rotation about z (such as
one rotated about x or y, or a general SO(3) frame), and a rotation angle that depends on t. Assumed: sp.expand is
canonical on these Laurent polynomials, `red` is canonical modulo c^2 + s^2 = 1 (argued above), and the constraints define
mu and J_i with K_ij = -g(nabla_{e_i} n, e_j), which `constraints_4d` and `constraints_4d_diag` and `constraints_4d_rot`
tie to G_ab on these classes of slices (N = 1, no shift).

    python3 notes/warp/weitzenbock_conf.py                  checks, then every control
    python3 notes/warp/weitzenbock_conf.py --no-controls    checks only
    python3 notes/warp/weitzenbock_conf.py --mutate=NAME    one control (its target checks must fail, exit 1)
    python3 notes/warp/weitzenbock_conf.py --no-controls --only=NAME,...   (or --skip=NAME,...) some diagonal and rotated checks
        of sections 6 to 13 only, for staged runs. Under --mutate such a check runs only if the control targets it.

Exit status is 1 if a check fails or a control does not fire. One process, about 4 minutes for the checks and about 1 hour
with the controls, 110 MiB.
"""
import importlib.util
import pathlib
import subprocess
import sys
import time
import types

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
# The rotated frame's identities (#130): IDENT_R are all of them, KIDR those with K != 0, ANIR those with A, B, C not constant
# and ROTR those with alpha not constant.
IDENT_R = ["identity_rot_alpha0", "identity_rot_flat_K0", "identity_rot_lichnerowicz", "identity_rot_flat_K",
           "identity_rot_alcubierre", "identity_rot_Kdiag", "identity_rot_general"]
KIDR = ["identity_rot_alpha0", "identity_rot_flat_K", "identity_rot_alcubierre", "identity_rot_Kdiag",
        "identity_rot_general"]
ANIR = ["identity_rot_alpha0", "identity_rot_lichnerowicz", "identity_rot_Kdiag", "identity_rot_general"]
ROTR = ["identity_rot_flat_K0", "identity_rot_lichnerowicz", "identity_rot_flat_K", "identity_rot_alcubierre",
        "identity_rot_Kdiag", "identity_rot_general"]
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
# The rotated frame's controls (#130). The old controls keep their targets and gain the rotated checks they must also fail.
CONTROLS.update({
    "rot-dalpha": ("drops the derivatives of alpha in the spin connection of the rotated frame",
                   ["planted_omega_rot", "omega_rot_conv"] + ROTR),
    "rot-reverse": ("rotates the frame by -alpha, where the hand values and conventions.py rotate by +alpha",
                    ["planted_omega_rot", "omega_rot_conv", "planted_trace_rot", "planted_alcubierre_rot",
                     "planted_kasner_rot"]),
    "rot-swap": ("takes c and s swapped, the frame rotated by pi/2 - alpha, which is orthonormal but is not alpha = 0 at "
                 "alpha = 0", ["planted_omega_rot", "omega_rot_conv", "rot_alpha0_limit", "identity_rot_alpha0",
                               "planted_trace_rot", "planted_alcubierre_rot", "planted_kasner_rot"]),
    "rot-ortho": ("rotates e_1 only, so the frame is not orthonormal",
                  ["orthonormal_rot", "planted_omega_rot", "omega_rot_conv", "identity_rot_flat_K0",
                   "identity_rot_general"]),
    "rot-k": ("takes K'_ij as components in the unrotated frame in the constraints",
              ["planted_alcubierre_rot", "planted_kasner_rot", "constraints_cov_rot", "constraints_4d_rot",
               "identity_rot_flat_K", "identity_rot_alcubierre", "identity_rot_Kdiag", "identity_rot_general"]),
    "rot-J-frame": ("takes J(e'_i) as J(e_i), the components in the unrotated frame",
                    ["planted_alcubierre_rot", "planted_trace_rot", "constraints_cov_rot", "constraints_4d_rot",
                     "identity_rot_flat_K", "identity_rot_alcubierre", "identity_rot_Kdiag", "identity_rot_general"]),
    "rot-div": ("takes V = V'^a e_a, along the unrotated frame, in the divergence", ROTR),
    "rot-reduce": ("does not reduce modulo c^2 + s^2 = 1", ["orthonormal_rot"] + ROTR),
})
for _n, _extra in {
    "sign": KIDR, "Jsign": KIDR, "musign": ["identity_rot_lichnerowicz"] + KIDR, "Kdiv": KIDR + ["V_reduced_rot"],
    "Vspin": IDENT_R + ["V_reduced_rot"], "weight": ANIR, "omega-sign": ["clifford_parallel_rot"] + IDENT_R,
    "christoffel": ["planted_omega_rot", "omega_rot_conv", "clifford_parallel_rot"] + ANIR,
    "trace": KIDR + ["planted_alcubierre_rot", "planted_trace_rot", "constraints_4d_rot"],
    "Ksq": KIDR + ["planted_alcubierre_rot", "planted_trace_rot", "planted_kasner_rot", "constraints_4d_rot"],
    "Rsign": ["planted_K0_rot", "planted_trace_rot", "constraints_4d_rot"] + ANIR,
    "Kframe": ["planted_alcubierre_rot", "planted_trace_rot", "planted_kasner_rot", "constraints_cov_rot",
               "constraints_4d_rot"] + KIDR,
    "J-christoffel": ["planted_trace_rot", "constraints_4d_rot", "identity_rot_alpha0", "identity_rot_Kdiag",
                      "identity_rot_general"],
    "gamma": ["clifford_parallel_rot"] + IDENT_R, "no-conj": IDENT_R, "adm-K": ["constraints_4d_rot"],
    "J4d": ["constraints_4d_rot"],
    "inv-typo": ["planted_omega_rot", "omega_rot_conv", "clifford_parallel_rot", "planted_K0_rot", "planted_trace_rot",
                 "constraints_4d_rot"] + ANIR,
    "ricci-sign": ["planted_K0_rot", "planted_trace_rot", "constraints_4d_rot"] + ANIR,
    "diag-frame": ["orthonormal_rot", "planted_omega_rot"] + ANIR,
    "aniso-frame": ["orthonormal_rot", "planted_omega_rot", "identity_rot_alpha0", "identity_rot_lichnerowicz",
                    "identity_rot_Kdiag", "identity_rot_general"],
    "abc": ["identity_rot_alpha0", "identity_rot_lichnerowicz", "identity_rot_Kdiag", "identity_rot_general"],
    "omega-hand": ["planted_omega_rot"], "kasner-exponents": ["planted_kasner_rot"],
}.items():
    CONTROLS[_n] = (CONTROLS[_n][0], CONTROLS[_n][1] + _extra)
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
        check("planted_K0_diag", ok, "K = 0 on dx^2 + dy^2 + c(x,y)^2 dz^2 and on a(x)^2 dx^2 + w(x)^2 (dy^2 + dz^2): "
              "mu = R/16 pi with R by hand, J = 0")
    if todo("planted_trace_diag"):
        # K_ij = f delta_ij in the frame, so K_ab = f g_ab, P_ab = -2 f g_ab, 16 pi mu = R + 6 f^2, 8 pi J_a = -2 d_a f.
        ok = True
        for Ac, Rc in metrics:
            mu, J = constraints_diag(Ac, f * I3, GeoDiag(Ac))
            wantJ = sp.Matrix([-sp.diff(f, X3[i]) / (4 * sp.pi * Ac[i]) for i in range(3)])
            ok = ok and zero(mu - (Rc / (16 * sp.pi) + 3 * f**2 / (8 * sp.pi))) and is_zero(J - wantJ)
        check("planted_trace_diag", ok, "K_ij = f delta_ij on the same two metrics: mu = R/16 pi + 3 f^2/8 pi, "
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


def diag_terms(A, K, main=True):
    """The pieces of the identity on diag(A^2, B^2, C^2) in the frame e_i = d_i / A_i: nabla_i eps, V, and, if main, the left
    side and div V."""
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
    lhs = div = None
    if main:
        lhs = sp.expand(sum(inner(nabt[i], nabt[i]) for i in range(3)) - inner(Dteps, Dteps)
                        + 4 * sp.pi * (sM * mu * e2 + sJ * Jxi))
        sqg = A[0]**3 if MUT == "abc" else A[0] * A[1] * A[2]  # sqrt(det g) = A B C
        if MUT == "weight":
            div = sum(sp.diff(V[a], X3[a]) for a in range(3))
        else:  # div V = (1/ABC) d_a(ABC V^a / A_a), the coordinate components of V = V^a e_a being V^a sc_a
            div = sum(sp.diff(sqg * geo.sc[a] * V[a], X3[a]) for a in range(3)) / sqg
        div = sp.expand(div)
    return types.SimpleNamespace(nab=nab, P=P, xi=xi, V=V, lhs=lhs, div=div)


def identity_diag(name, A, K, text, reduced=False):
    run_main, run_red = todo(name), reduced and todo("V_reduced_diag")
    if not (run_main or run_red):
        return
    t0 = time.time()
    d = diag_terms(A, K, run_main)
    if run_main:
        res = sp.expand(d.lhs - d.div)
        check(name, res == 0, f"{text}; {nterms(d.lhs)} terms on the left, {nterms(d.div)} in div V, "
              f"{len(sp.Add.make_args(res)) if res != 0 else 0} in the difference")
        print(f"    ({name}: {time.time() - t0:.0f} s)", file=sys.stderr, flush=True)
    if run_red:
        gj = note_gammas()[1]
        eps, nab, V, P, xi = spinor(), d.nab, d.V, d.P, d.xi
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


# ---- The diagonal slice in the frame rotated about z, issue #130 -------------------------------------------------------------
# The frame e'_1 = c e_1 + s e_2, e'_2 = -s e_1 + c e_2, e'_3 = e_3, c = cos(alpha), s = sin(alpha), alpha(x, y, z) arbitrary,
# on diag(A^2, B^2, C^2) with e_i = d_i / A_i. K'_ij = K(e'_i, e'_j) and eps are arbitrary, V = V'^a e'_a. GeoRot,
# constraints_rot and rot_terms are written beside GeoDiag, constraints_diag and diag_terms and call none of the old
# constraints and identity code. At alpha = 0 they are compared with the diagonal code (rot_alpha0_limit,
# identity_rot_alpha0). Expressions hold cos(alpha) and sin(alpha), which sympy differentiates by the chain rule, so the
# derivatives of alpha appear as jets, and `red` reduces modulo c^2 + s^2 = 1.
ALPHA = fn("alpha")  # the angle of the rotated frame
def trig_reduce(e):
    """Reduces the expanded e modulo c^2 + s^2 - 1, c = cos(alpha), s = sin(alpha). The polynomial c^2 + s^2 - 1 is monic in c,
    so division by it, with the coefficients (Laurent polynomials in the jets of A, B, C, alpha, K'_ij, u_a, v_a) as the
    ground ring, has a unique remainder of degree at most 1 in c: c^k goes to c^(k mod 2) (1 - s^2)^(k div 2). The remainder
    is zero exactly when e is in the ideal. That ideal is all of the polynomials that vanish for every real alpha: the jets
    of alpha at a point are independent of alpha itself, so the only relation among c, s and the jets is c^2 + s^2 = 1, and
    the real circle is dense in the irreducible conic c^2 + s^2 = 1, so a polynomial vanishing on the circle is a multiple of
    c^2 + s^2 - 1. This is the reduction of conventions.py's trig_reduce (#86) for alpha(x, y, z), grouped term by term
    because sympy's collect rejects mixed derivatives such as f_xy, and mapped back to cos and sin."""
    c, s = sp.Dummy("c"), sp.Dummy("s")
    e = e.xreplace({sp.cos(ALPHA): c, sp.sin(ALPHA): s})
    groups = {}  # the terms of e by their power k of c, each a monomial of the expanded e
    for term in sp.Add.make_args(e):
        k = term.as_powers_dict().get(c, 0)
        groups.setdefault(k, []).append(term * c**-k)
    out = sum(sp.Add(*co) * c ** (k % 2) * (1 - s**2) ** (k // 2) for k, co in groups.items())
    return sp.expand(out).xreplace({c: sp.cos(ALPHA), s: sp.sin(ALPHA)})  # back, since reduced results are differentiated


def red(e):
    """The canonical form of e: expanded, and reduced modulo c^2 + s^2 = 1 if it holds cos(alpha) or sin(alpha)."""
    e = sp.expand(e)
    if MUT != "rot-reduce" and e != 0 and e.has(sp.cos(ALPHA), sp.sin(ALPHA)):
        e = trig_reduce(e)
    return e


def zero_rot(e):
    return red(e) == 0


def is_zero_rot(M):
    return all(zero_rot(e) for e in sp.Matrix(M))


def rot_matrix(alpha):
    """R with e'_i = R_ik e_k, the rotation of the hand values."""
    c, s = sp.cos(alpha), sp.sin(alpha)
    return sp.Matrix([[c, s, 0], [-s, c, 0], [0, 0, 1]])


def rot_frame(E, alpha):
    """The frame E rotated by alpha about the third axis, e'_1 = c e_1 + s e_2, e'_2 = -s e_1 + c e_2, e'_3 = e_3."""
    c, s = sp.cos(alpha), sp.sin(alpha)
    if MUT == "rot-reverse":
        s = -s
    if MUT == "rot-swap":
        c, s = s, c
    if MUT == "rot-ortho":
        return [c * E[0] + s * E[1], E[1], E[2]]
    return [c * E[0] + s * E[1], -s * E[0] + c * E[1], E[2]]


def nabla_vec(Chr, U, V):
    """(nabla_U V)^m = U^k d_k V^m + Gamma^m_kl U^k V^l, for coordinate components U and V."""
    sg = -1 if MUT == "christoffel" else 1
    return sp.Matrix([sum(U[k] * sp.diff(V[m], X3[k]) for k in range(3))
                      + sg * sum(Chr[m][k][l] * U[k] * V[l] for k in range(3) for l in range(3)) for m in range(3)])


class GeoRot(GeoDiag):
    """The slice diag(A^2, B^2, C^2) in the frame e'_a rotated by alpha about z: F[i][m] the coordinate components of e'_i,
    th[i][m] = g(e'_i, d_m) its coframe, omr[j][k][i] = omega'_jk(e'_i) = g(e'_j, nabla_{e'_i} e'_k) from the Christoffel
    symbols, reduced modulo c^2 + s^2 = 1, and Omr_i = (1/4) omega'_jk(e'_i) g^j g^k. The unrotated om and Om of GeoDiag stay."""

    def __init__(self, A, alpha):
        super().__init__(A)
        self.alpha = alpha
        E = [sp.Matrix([int(m == i) for m in range(3)]) * self.sc[i] for i in range(3)]
        Er = rot_frame(E, alpha)
        a0 = sp.Symbol("a0", real=True)  # a frozen angle, whose derivatives vanish
        Es = rot_frame(E, a0) if MUT == "rot-dalpha" else Er  # the frame nabla differentiates
        back = {a0: alpha} if MUT == "rot-dalpha" else {}
        self.F = [list(v) for v in Er]
        self.th = [[self.g[m, m] * self.F[i][m] for m in range(3)] for i in range(3)]
        self.framer = sp.Matrix(3, 3, lambda i, j: red((Er[i].T * self.g * Er[j])[0]))
        dE = [[self.g * nabla_vec(self.Chr, Es[i], Es[k]) for k in range(3)] for i in range(3)]
        self.omr = [[[red(sp.expand((Es[j].T * dE[i][k])[0]).xreplace(back)) for i in range(3)] for k in range(3)]
                    for j in range(3)]
        g0, gj = note_gammas()
        ss = -1 if MUT == "omega-sign" else 1
        self.Omr = [(ss * sum((self.omr[j][k][i] * gj[j] * gj[k] for j in range(3) for k in range(3)),
                              sp.zeros(4)) / 4).applyfunc(red) for i in range(3)]


def omega_hand_rot(A, alpha):
    """omega'_jk(e'_i) by hand: the closed form omega_lm(e_p) of omega_hand_diag rotated, R_ip R_jl R_km omega_lm(e_p), plus
    the inhomogeneous term of the rotation, omega'_21 = e'_i(alpha) = -omega'_12 (#86, conventions.py section 6(d))."""
    R = rot_matrix(alpha)
    w = omega_hand_diag(A)
    ea = [sum(R[i, p] * sp.diff(alpha, X3[p]) / A[p] for p in range(3)) for i in range(3)]
    return [[[red(sum(R[i, p] * R[j, l] * R[k, m] * w[l][m][p] for p in range(3) for l in range(3) for m in range(3))
                  + (int(j == 1 and k == 0) - int(j == 0 and k == 1)) * ea[i]) for i in range(3)] for k in range(3)]
            for j in range(3)]


CONV = None


def conventions():
    """conventions.py, loaded once, with its argv hidden as in section_gammas."""
    global CONV
    if CONV is None:
        spec = importlib.util.spec_from_file_location("conventions_rot", pathlib.Path(__file__).with_name("conventions.py"))
        CONV = importlib.util.module_from_spec(spec)
        saved, sys.argv = sys.argv, sys.argv[:1]
        try:
            spec.loader.exec_module(CONV)
        finally:
            sys.argv = saved
    return CONV


def constraints_rot(A, K, geo):
    """mu and the frame components J'_i = J(e'_i) from the constraints, computed in coordinates from the frame components
    K'_ij = K(e'_i, e'_j), whose coordinate components are K_mn = K'_ij th^i_m th^j_n."""
    g, gi, Chr = geo.g, geo.gi, geo.Chr
    if MUT == "Kframe":
        Kc = K
    elif MUT == "rot-k":  # K'_ij taken in the unrotated frame
        Kc = sp.Matrix(3, 3, lambda a, b: A[a] * A[b] * K[a, b])
    else:
        Kc = sp.Matrix(3, 3, lambda a, b: sum(K[i, j] * geo.th[i][a] * geo.th[j][b] for i in range(3) for j in range(3)))
    Kc = Kc.applyfunc(red)
    trK = red(sum(gi[a, b] * Kc[a, b] for a in range(3) for b in range(3)))
    K2 = red(sum(gi[a, c] * gi[b, d] * Kc[a, b] * Kc[c, d] for a in range(3) for b in range(3) for c in range(3)
                 for d in range(3)))
    sR = -1 if MUT == "Rsign" else 1
    sq = 1 if MUT == "Ksq" else -1
    mu = red((sR * geo.R + trK**2 + sq * K2) / (16 * sp.pi))
    tr = 1 if MUT == "trace" else -1
    Pc = (Kc + tr * g * trK).applyfunc(red)  # P_ab = K_ab - g_ab tr K

    def DP(c, a, b):  # nabla_c P_ab
        ch = 0 if MUT == "J-christoffel" else 1
        return sp.diff(Pc[a, b], X3[c]) - ch * sum(Chr[d][c][a] * Pc[d, b] + Chr[d][c][b] * Pc[a, d] for d in range(3))
    Jc = [red(sum(gi[a, c] * DP(c, a, b) for a in range(3) for c in range(3)) / (8 * sp.pi)) for b in range(3)]  # J(d_b)
    if MUT == "rot-J-frame":  # J(e'_i) taken as J(e_i)
        return mu, sp.Matrix([red(Jc[i] / A[i]) for i in range(3)])
    return mu, sp.Matrix([red(sum(geo.F[i][b] * Jc[b] for b in range(3))) for i in range(3)])


def rot_terms(A, K, alpha, main=True):
    """The pieces of the identity in the rotated frame: nabla_i eps, V and, if main, the left side and div V, all reduced
    modulo c^2 + s^2 = 1."""
    geo = GeoRot(A, alpha)
    mu, J = constraints_rot(A, K, geo)
    g0, gj = note_gammas()
    eps = spinor()
    trK = K.trace()
    P = K - trK * I3
    sK = -1 if MUT == "sign" else 1
    Aterm = [sK * sum((K[i, j] * gj[j] * g0 for j in range(3)), sp.zeros(4)) / 2 for i in range(3)]
    nab = [(sum((geo.F[i][m] * eps.diff(X3[m]) for m in range(3)), sp.zeros(4, 1)) + geo.Omr[i] * eps).applyfunc(red)
           for i in range(3)]  # nabla_i eps = e'_i(eps) + Omr_i eps
    nabt = [(nab[i] + Aterm[i] * eps).applyfunc(red) for i in range(3)]
    Deps = sum((gj[i] * nab[i] for i in range(3)), sp.zeros(4, 1)).applyfunc(red)
    Dteps = sum((gj[i] * nabt[i] for i in range(3)), sp.zeros(4, 1)).applyfunc(red)
    e2 = inner(eps, eps)
    xi = [inner(eps, g0 * gj[j] * eps) for j in range(3)]
    Jxi = red(sum(J[j] * xi[j] for j in range(3)))
    sM = -1 if MUT == "musign" else 1
    sJ = 1 if MUT == "Jsign" else -1
    V = []
    for i in range(3):
        spin = nab[i] if MUT == "Vspin" else (nab[i] - gj[i] * Deps).applyfunc(red)
        Pxi = 0 if MUT == "Kdiv" else sum(P[i, j] * xi[j] for j in range(3)) / 2
        V.append(red(re(inner(eps, spin)) - Pxi))
    lhs = div = None
    if main:
        lhs = red(sum(red(inner(nabt[i], nabt[i])) for i in range(3)) - red(inner(Dteps, Dteps))
                  + 4 * sp.pi * (sM * mu * e2 + sJ * Jxi))
        sqg = A[0]**3 if MUT == "abc" else A[0] * A[1] * A[2]  # sqrt(det g) = A B C
        # the coordinate components of V = V'^a e'_a, which with rot-div are taken along the unrotated e_a
        Fd = [[self_sc * int(m == a) for m in range(3)] for a, self_sc in enumerate(geo.sc)] if MUT == "rot-div" else geo.F
        Vc = [sum(V[a] * Fd[a][m] for a in range(3)) for m in range(3)]
        if MUT == "weight":
            div = sum(sp.diff(Vc[m], X3[m]) for m in range(3))
        else:  # div V = (1/ABC) d_m(ABC V^m)
            div = sum(sp.diff(sqg * Vc[m], X3[m]) for m in range(3)) / sqg
        div = red(div)
    return types.SimpleNamespace(geo=geo, nab=nab, P=P, xi=xi, V=V, lhs=lhs, div=div)


def identity_rot(name, A, K, alpha, text, reduced=False):
    run_main, run_red = todo(name), reduced and todo("V_reduced_rot")
    if not (run_main or run_red):
        return
    t0 = time.time()
    d = rot_terms(A, K, alpha, run_main)
    if run_main:
        res = red(d.lhs - d.div)
        check(name, res == 0, f"{text}; {nterms(d.lhs)} terms on the left, {nterms(d.div)} in div V, "
              f"{len(sp.Add.make_args(res)) if res != 0 else 0} in the difference")
        print(f"    ({name}: {time.time() - t0:.0f} s)", file=sys.stderr, flush=True)
    if run_red:
        gj = note_gammas()[1]
        eps, nab, V, P, xi = spinor(), d.nab, d.V, d.P, d.xi
        Vr = [red(-sum(re(inner(eps, gj[i] * gj[j] * nab[j])) for j in range(3) if j != i)
                  - sum(P[i, j] * xi[j] for j in range(3)) / 2) for i in range(3)]
        check("V_reduced_rot", all(zero_rot(V[i] - Vr[i]) for i in range(3)),
              "V^i = -sum_{j != i} Re<eps, g^i g^j nabla_j eps> - (1/2) P_ij xi^j in the rotated frame")


def hand_metrics():
    """Two metrics with a hand value of R, as in section_constraints_diag: dx^2 + dy^2 + c(x, y)^2 dz^2 and a(x)^2 dx^2 +
    w(x)^2 (dy^2 + dz^2)."""
    one = sp.Integer(1)
    c2 = sp.Function("c")(X, Y)
    a_, w_ = sp.Function("a")(X), sp.Function("w")(X)
    return [([one, one, c2], -2 * (sp.diff(c2, X, 2) + sp.diff(c2, Y, 2)) / c2),
            ([a_, w_, w_], -4 * sp.diff(sp.diff(w_, X) / a_, X) / (a_ * w_) - 2 * sp.diff(w_, X)**2 / (a_ * w_)**2)]


def section_planted_rot():
    names = ["orthonormal_rot", "planted_omega_rot", "omega_rot_conv", "clifford_parallel_rot", "rot_alpha0_limit"]
    if not any(todo(n) for n in names):
        return
    print("\n== 10. The diagonal slice in the frame rotated about z: omega', Clifford parallelism, alpha = 0 ==")
    A = [fn(n) for n in "ABC"]
    geo = GeoRot(A, ALPHA)
    if todo("orthonormal_rot"):
        check("orthonormal_rot", is_zero_rot(geo.framer - I3),
              "g(e'_i, e'_j) = delta_ij for e'_1 = c e_1 + s e_2, e'_2 = -s e_1 + c e_2, e'_3 = e_3, modulo c^2 + s^2 = 1")
    if todo("planted_omega_rot"):
        hand = omega_hand_rot(A, ALPHA)
        check("planted_omega_rot", all(zero_rot(geo.omr[j][k][i] - hand[j][k][i]) for i in range(3) for j in range(3)
                                       for k in range(3)),
              "omega'_jk(e'_i) = g(e'_j, nabla_{e'_i} e'_k) from the Christoffel symbols is R_ip R_jl R_km omega_lm(e_p) "
              "+ (delta_j2 delta_k1 - delta_j1 delta_k2) e'_i(alpha), with Cartan's omega_lm of the unrotated frame")
    if todo("omega_rot_conv"):
        om3, _ = conventions().slice_data(sp.Integer(1), sp.Integer(1), [0, 0, 0], A, alpha=ALPHA)
        check("omega_rot_conv", all(zero_rot(geo.omr[j][k][i] - om3[j][k][i]) for i in range(3) for j in range(3)
                                    for k in range(3)),
              "omega'_jk(e'_i) equals #86's rotated-frame connection, conventions.slice_data(A, alpha=alpha), "
              "modulo c^2 + s^2 = 1")
    if todo("clifford_parallel_rot"):
        g0, gj = note_gammas()
        par = all(is_zero_rot(geo.Omr[i] * gj[l] - gj[l] * geo.Omr[i]
                              - sum((geo.omr[j][l][i] * gj[j] for j in range(3)), sp.zeros(4)))
                  for i in range(3) for l in range(3))
        check("clifford_parallel_rot", par, "[(1/4) omega'_jk g^j g^k, g^l] = omega'_jl g^j in the rotated frame")
    if todo("rot_alpha0_limit"):
        gr, gd = GeoRot(A, sp.Integer(0)), GeoDiag(A)
        Kg = sym_K()
        mur, Jr = constraints_rot(A, Kg, gr)
        mud, Jd = constraints_diag(A, Kg, gd)
        same = is_zero(gr.g - gd.g) and zero(gr.R - gd.R) and zero(mur - mud) and is_zero(Jr - Jd)
        same = same and all(zero(gr.omr[j][k][i] - gd.om[j][k][i]) for i in range(3) for j in range(3) for k in range(3))
        same = same and all(is_zero(gr.Omr[i] - gd.Om[i]) for i in range(3))
        check("rot_alpha0_limit", same, "alpha = 0 gives #99's g, R, omega_jk, Om_i and, for K_ij arbitrary, mu and J_i")


def section_constraints_rot():
    names = ["planted_K0_rot", "planted_trace_rot", "planted_alcubierre_rot", "planted_kasner_rot", "constraints_cov_rot"]
    if not any(todo(n) for n in names):
        return
    print("\n== 11. Planted truth for mu and J'_i from the constraints, in the rotated frame ==")
    one = sp.Integer(1)
    f, v = fn("f"), sp.Symbol("v", real=True)
    R = rot_matrix(ALPHA)
    if todo("planted_K0_rot"):
        ok = True
        for Ac, Rc in hand_metrics():
            mu, J = constraints_rot(Ac, sp.zeros(3), GeoRot(Ac, ALPHA))
            ok = ok and zero_rot(mu - Rc / (16 * sp.pi)) and is_zero_rot(J)
        check("planted_K0_rot", ok, "K = 0 on dx^2 + dy^2 + c(x,y)^2 dz^2 and on a(x)^2 dx^2 + w(x)^2 (dy^2 + dz^2), alpha "
              "arbitrary: mu = R/16 pi with R by hand, J = 0")
    if todo("planted_trace_rot"):
        # K'_ij = f delta_ij is the same in every frame: K_ab = f g_ab, so 16 pi mu = R + 6 f^2 and 8 pi J_a = -2 d_a f, whose
        # frame components are J(e'_i) = -e'_i(f)/4 pi with e'_i = R_ik e_k, e_k = d_k / A_k.
        ok = True
        for Ac, Rc in hand_metrics():
            mu, J = constraints_rot(Ac, f * I3, GeoRot(Ac, ALPHA))
            wantJ = sp.Matrix([-sum(R[i, k] * sp.diff(f, X3[k]) / Ac[k] for k in range(3)) / (4 * sp.pi) for i in range(3)])
            ok = ok and zero_rot(mu - (Rc / (16 * sp.pi) + 3 * f**2 / (8 * sp.pi))) and is_zero_rot(J - wantJ)
        check("planted_trace_rot", ok, "K'_ij = f delta_ij, alpha arbitrary, on the same two metrics: mu = R/16 pi + "
              "3 f^2/8 pi, J'_i = -e'_i(f)/4 pi")
    if todo("planted_alcubierre_rot"):
        # Alcubierre's flat data K_ij = (d_i b_j + d_j b_i)/2, b = (-v f, 0, 0), with K'_ij = R K R^T in the rotated frame,
        # so mu is as in the unrotated frame and J'_i = R_ik J_k, with #58 item 1's J_k.
        b = [-v * f, 0, 0]
        Ka = sp.Matrix(3, 3, lambda i, j: (sp.diff(b[j], X3[i]) + sp.diff(b[i], X3[j])) / 2)
        Kp = (R * Ka * R.T).applyfunc(red)
        mu, J = constraints_rot([one] * 3, Kp, GeoRot([one] * 3, ALPHA))
        J0 = sp.Matrix([v / (16 * sp.pi) * (sp.diff(f, X, X3[i]) - int(i == 0) * lap(f)) for i in range(3)])
        check("planted_alcubierre_rot", zero_rot(mu + v**2 / (32 * sp.pi) * (sp.diff(f, Y)**2 + sp.diff(f, Z)**2))
              and is_zero_rot(J - R * J0),
              "A = B = C = 1, K' = R K R^T with K Alcubierre's, alpha and f arbitrary: mu = -(v^2/32 pi)(f_y^2 + f_z^2), "
              "J'_i = R_ik J_k with J_k = (v/16 pi)(d_x d_k f - delta_kx Laplace f)")
    if todo("planted_kasner_rot"):
        # K = diag(p) of a Kasner slice at t = 1 on constant A_i = a_i, whose frame components in the rotated frame are
        # K' = R diag(p) R^T for every alpha, still vacuum, so mu = J' = 0.
        a = [sp.Symbol(f"a{i + 1}", positive=True) for i in range(3)]
        p = [sp.Rational(1, 2)] * 3 if MUT == "kasner-exponents" else [sp.Rational(-2, 7), sp.Rational(3, 7),
                                                                       sp.Rational(6, 7)]
        Kp = (R * sp.diag(*p) * R.T).applyfunc(red)
        mu, J = constraints_rot(a, Kp, GeoRot(a, ALPHA))
        check("planted_kasner_rot", zero_rot(mu) and is_zero_rot(J),
              "constant A_i = a_i, K' = R diag(-2/7, 3/7, 6/7) R^T, alpha arbitrary: mu = J' = 0")
    if todo("constraints_cov_rot"):
        # Frame covariance through #99's code: mu and J of K'_ij in the rotated frame are mu and R J of K = R^T K' R in the
        # unrotated frame, which constraints_diag takes as it is.
        A = [fn(n) for n in "ABC"]
        Kp = sym_K()
        Ku = (R.T * Kp * R).applyfunc(red)
        mud, Jd = constraints_diag(A, Ku, GeoDiag(A))
        mur, Jr = constraints_rot(A, Kp, GeoRot(A, ALPHA))
        check("constraints_cov_rot", zero_rot(mur - mud) and is_zero_rot(Jr - R * Jd),
              "mu and J'_i of K'_ij in the rotated frame are mu and R_ik J_k of K = R^T K' R through constraints_diag")


def section_4d_rot():
    if not todo("constraints_4d_rot"):
        return
    print("\n== 12. The constraints against G_ab of the 4-metric, rotated frame ==")
    # As section 8, with K = R^T K' R, the components of K' in the unrotated frame, built here from the hand matrix R and not
    # from the coframe of GeoRot: g = -dt^2 + Lam (delta - 2 t K) Lam, Lam = diag(A, B, C), and then J'_i = R_ik J_k.
    A = [fn(n) for n in "ABC"]
    Kp = sym_K()
    R = rot_matrix(ALPHA)
    K = (R.T * Kp * R).applyfunc(red)
    sg = -1 if MUT == "adm-K" else 1
    K2 = (K * K).applyfunc(red)
    Lam, Lami = sp.diag(*A), sp.diag(*[1 / a for a in A])
    gs = (Lam * (I3 - sg * 2 * T * K) * Lam).applyfunc(red)
    gis = (Lami * (I3 + sg * 2 * T * K + 4 * T**2 * K2) * Lami).applyfunc(red)
    prod = (gs * gis).applyfunc(red)
    inv_ok = all(red(prod[i, j] - int(i == j)).coeff(T, p) == 0 for i in range(3) for j in range(3) for p in range(3))
    g, gi = sp.zeros(4), sp.zeros(4)
    g[0, 0], gi[0, 0] = -1, -1
    for i in range(3):
        for j in range(3):
            g[i + 1, j + 1], gi[i + 1, j + 1] = gs[i, j], gis[i, j]
    Chr = christoffel(g, gi, CRD)
    Chr = [[[red(c.subs(T, 0) + T * sp.diff(c, T).subs(T, 0)) for c in row] for row in blk] for blk in Chr]  # O(t)
    Rab = ricci(Chr, CRD).subs(T, 0)
    gi0 = gi.subs(T, 0)
    Rs = red(sum(gi0[a, b] * Rab[a, b] for a in range(4) for b in range(4)))
    Gtt = red(Rab[0, 0] + Rs / 2)
    sj = 1 if MUT == "J4d" else -1
    mu4 = Gtt / (8 * sp.pi)
    J4 = R * sp.Matrix([sj * Rab[0, i + 1] / (8 * sp.pi) / A[i] for i in range(3)])  # J'_i = R_ik J_k, J_k = -G_tk / 8 pi
    mu, J = constraints_rot(A, Kp, GeoRot(A, ALPHA))
    check("constraints_4d_rot", inv_ok and zero_rot(mu - mu4) and is_zero_rot(J - J4),
          "the constraints' mu and J'_i equal G_tt / 8 pi and -R_ik G_tk / 8 pi of -dt^2 + Lam (delta - 2 t K) Lam at t = 0, "
          "K = R^T K' R, for arbitrary alpha, A, B, C and K'_ij")


def identity_rot_alpha0(A, K):
    name = "identity_rot_alpha0"
    if not todo(name):
        return
    t0 = time.time()
    r = rot_terms(A, K, sp.Integer(0))
    d = diag_terms(A, K)
    res = red(r.lhs - r.div)
    same = zero(r.lhs - d.lhs) and zero(r.div - d.div) and all(zero(r.V[i] - d.V[i]) for i in range(3))
    check(name, res == 0 and same, f"alpha = 0, A, B, C, K_ij and eps arbitrary; {nterms(r.lhs)} terms on the left, "
          f"{nterms(r.div)} in div V, the difference 0, and the left side, div V and V equal identity_diag_general's "
          f"term for term")
    print(f"    ({name}: {time.time() - t0:.0f} s)", file=sys.stderr, flush=True)


def section_identity_rot():
    if not any(todo(n) for n in IDENT_R + ["V_reduced_rot"]):
        return
    print("\n== 13. The identity |nabla~ eps|^2 - |D~ eps|^2 + 4 pi (mu |eps|^2 - J.xi) = div V in the frame rotated about z ==")
    print("    V = V'^a e'_a, div V = (1/ABC) d_m(ABC V^m) with V^m = V'^a e'^m_a, everything reduced modulo c^2 + s^2 = 1")
    one, Kz = sp.Integer(1), sp.zeros(3)
    A = [fn(n) for n in "ABC"]
    Kg = sym_K()
    f, v = fn("f"), sp.Symbol("v", real=True)
    b = [-v * f, 0, 0]
    Ka = sp.Matrix(3, 3, lambda i, j: (sp.diff(b[j], X3[i]) + sp.diff(b[i], X3[j])) / 2)
    R = rot_matrix(ALPHA)
    Kar = (R * Ka * R.T).applyfunc(red)
    Kd = sp.diag(*[fn(f"k{i + 1}") for i in range(3)])
    print("  Planted first: alpha = 0 against identity_diag_general, then A = B = C = 1 and K = 0 in the rotated frame")
    identity_rot_alpha0(A, Kg)
    identity_rot("identity_rot_flat_K0", [one] * 3, Kz, ALPHA, "A = B = C = 1, K = 0, alpha arbitrary: the flat identity")
    identity_rot("identity_rot_lichnerowicz", A, Kz, ALPHA, "A, B, C, alpha arbitrary, K = 0: Lichnerowicz")
    identity_rot("identity_rot_flat_K", [one] * 3, Kg, ALPHA, "A = B = C = 1, K'_ij and alpha arbitrary")
    identity_rot("identity_rot_alcubierre", [one] * 3, Kar, ALPHA, "A = B = C = 1, K' = R K R^T with K Alcubierre's with "
                 "b = (-v f, 0, 0), f and alpha arbitrary")
    identity_rot("identity_rot_Kdiag", A, Kd, ALPHA, "A, B, C, alpha arbitrary, K'_ij = diag(k1, k2, k3) arbitrary "
                 "functions, eps arbitrary")
    identity_rot("identity_rot_general", A, Kg, ALPHA, "alpha, A, B, C, K'_ij and eps all arbitrary", reduced=True)


def main():
    t0 = time.time()
    print("weitzenbock_conf.py, issues #87, #99 and #130. sympy %s, python %s. mutation: %s" % (sp.__version__, sys.version.split()[0],
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
    section_planted_rot()
    section_constraints_rot()
    section_4d_rot()
    section_identity_rot()
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

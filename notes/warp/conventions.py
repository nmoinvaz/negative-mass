#!/usr/bin/env python3
"""The spinor conventions of mass_bounds.md, checked symbolically, for issue #62.

On the flat slice of ds^2 = -dt^2 + sum_i (dx^i + b^i dt)^2 (lapse 1, any shift b^i(t, x, y, z)), the note's
nabla~_i = d_i + (1/2) K_ij g^j g^0, with the note's gammas g^a, equals the restriction to the slice of the spacetime
spin connection (1/4) omega_ab G^a G^b in the gammas G^0 = -g^0, G^j = g^j. The other sign of the K term does not.

The Levi-Civita connection comes from the Christoffel symbols of the coordinate metric, with
omega_ab(X) = g(e_a, nabla_X e_b), e_0 = n = d_t - b^i d_i and e_i = d_i. It shares no code with review_cartan.py,
which solves Cartan's structure equations instead, or with alcubierre_bounds.py and review_spinor.py. It shares sympy
1.14.0 and the gammas of mass_bounds.md, "Conventions", which are typed in again here.

    python3 notes/warp/conventions.py                  checks, then every control
    python3 notes/warp/conventions.py --no-controls    checks only
    python3 notes/warp/conventions.py --mutate=NAME    one control (its target checks must fail, exit 1)

Exit status is 1 if a check fails or a control does not fire. A few seconds, one process.
"""
import subprocess
import sys
import time

import sympy as sp

OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate", "")
CONTROLS = {
    "sign": ("writes the stated nabla~ with the other sign, d_i - (1/2) K_ij g^j g^0", ["restriction", "dirac"]),
    "christoffel": ("flips the sign of the Christoffel term in nabla_X Y", ["planted_desitter", "omega_slice"]),
    "spin-sign": ("uses -(1/4) omega_ab G^a G^b for the spin connection", ["clifford_parallel", "restriction",
                                                                           "note_gammas"]),
    "frame": ("takes e_0 = d_t + b^i d_i, the wrong sign of the shift", ["orthonormal"]),
    "gamma": ("scales g^1 by 11/10", ["clifford", "causal"]),
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


def is_zero(M):
    return sp.Matrix(M).applyfunc(lambda e: sp.simplify(sp.expand(e))) == sp.zeros(*sp.Matrix(M).shape)


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


# ---- The Levi-Civita connection in the frame (n, d_i) -----------------------------------------------------------------
def connection(b):
    """omega[a][b][c] = g(e_a, nabla_{e_c} e_b) for ds^2 = -dt^2 + sum (dx^i + b^i dt)^2, e_0 = d_t - b^i d_i."""
    g = sp.zeros(4)
    g[0, 0] = -1 + sum(bi**2 for bi in b)
    for i in range(3):
        g[0, i + 1] = g[i + 1, 0] = b[i]
        g[i + 1, i + 1] = 1
    gi = sp.zeros(4)  # inverse, typed in and verified
    gi[0, 0] = -1
    for i in range(3):
        gi[0, i + 1] = gi[i + 1, 0] = b[i]
        for j in range(3):
            gi[i + 1, j + 1] = int(i == j) - b[i] * b[j]
    inv_ok = is_zero(g * gi - sp.eye(4))
    dg = [[[sp.diff(g[m, n], CRD[l]) for l in range(4)] for n in range(4)] for m in range(4)]
    Chr = [[[sp.expand(sum(gi[m, s] * (dg[s][l][n] + dg[s][n][l] - dg[n][l][s]) for s in range(4)) / 2)
             for l in range(4)] for n in range(4)] for m in range(4)]  # Chr[m][n][l] = Gamma^m_nl
    sgn = -1 if MUT == "christoffel" else 1
    sh = 1 if MUT == "frame" else -1
    E = [sp.Matrix([1] + [sh * bi for bi in b])] + [sp.Matrix([int(m == i + 1) for m in range(4)]) for i in range(3)]

    def nabla(U, V):  # (nabla_U V)^m = U^n d_n V^m + Gamma^m_nl U^n V^l
        return sp.Matrix([sp.expand(sum(U[n] * sp.diff(V[m], CRD[n]) for n in range(4))
                                    + sgn * sum(Chr[m][n][l] * U[n] * V[l] for n in range(4) for l in range(4)))
                          for m in range(4)])

    frame = sp.Matrix(4, 4, lambda a, c: sp.expand((E[a].T * g * E[c])[0]))
    om = [[[sp.expand((E[a].T * g * nabla(E[c], E[bb]))[0]) for c in range(4)] for bb in range(4)] for a in range(4)]
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

    def spin(c, Gam):  # (1/4) omega_ab(e_c) G^a G^b
        return sp.expand(ss * sum((om[a][bb][c] * Gam[a] * Gam[bb] for a in range(4) for bb in range(4)),
                                  sp.zeros(4)) / 4)

    # nabla_X (G(Y) eps) = G(nabla_X Y) eps + G(Y) nabla_X eps, i.e. [Gamma_X, G_d] = omega^a_d(X) G_a, with
    # omega^a_d = eta^aa omega_ad. This fixes the sign and the 1/4 of the spin connection.
    par = all(is_zero(spin(c, Gs) * Gl[d] - Gl[d] * spin(c, Gs)
                      - sum((ETA[a, a] * om[a][d][c] * Gl[a] for a in range(4)), sp.zeros(4)))
              for c in range(4) for d in range(4))
    check("clifford_parallel", par, "(1/4) omega_ab G^a G^b keeps Clifford multiplication parallel, on every e_c")

    sK = -1 if MUT == "sign" else 1
    stated = [sK * sum((K[i, j] * gj[j] * g0 for j in range(3)), sp.zeros(4)) / 2 for i in range(3)]
    res = all(is_zero(spin(i + 1, Gs) - stated[i]) for i in range(3))
    check("restriction", res,
          "nabla_{e_i} = d_i + (1/4) omega_ab(e_i) G^a G^b with G^0 = -g^0 equals the stated "
          "nabla~_i = d_i %s (1/2) K_ij g^j g^0" % ("-" if sK < 0 else "+"))
    # In the note's own gammas, the same connection has the other sign (review.md section 3, gap 1).
    other = all(is_zero(spin(i + 1, [g0] + gj) + sum((K[i, j] * gj[j] * g0 for j in range(3)), sp.zeros(4)) / 2)
                for i in range(3))
    check("note_gammas", other, "in the note's gammas themselves it is d_i - (1/2) K_ij g^j g^0")
    # The Witten operator of item 5, D~ = g^i nabla~_i = D + (1/2)(tr K) g^0.
    dirac = is_zero(sum((gj[i] * stated[i] for i in range(3)), sp.zeros(4)) - K.trace() * g0 / 2)
    check("dirac", dirac, "g^i (stated K term)_i = (1/2)(tr K) g^0, so D~ = D + (1/2)(tr K) g^0")


def main():
    t0 = time.time()
    print("conventions.py, issue #62. sympy %s, python %s. mutation: %s"
          % (sp.__version__, sys.version.split()[0], MUT or "none"))
    section_gammas()
    section_planted()
    section_slice()
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

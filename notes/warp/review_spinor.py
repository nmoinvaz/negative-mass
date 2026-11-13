#!/usr/bin/env python3
"""Independent checks of #58's existence lemma (mass_bounds.md item 5) and identity (item 6), for issue #59.

Written for the review. It shares no code with alcubierre_bounds.py or phantom_source.py, only numpy, scipy and
sympy. No Dirac equation is solved here. The identity is tested as the integrated Weitzenbock formula, which holds
for every spinor eps = eps0 + psi with psi decaying, with analytic derivatives on a grid.

    python3 notes/warp/review_spinor.py                  checks, then every control
    python3 notes/warp/review_spinor.py --no-controls    checks only
    python3 notes/warp/review_spinor.py --mutate=NAME    one control (its target checks must fail, exit 1)

Exit status is 1 if a check fails or a control does not fire. About 1 GB and under a minute at nice 10.
"""
import subprocess
import sys
import time

import numpy as np
import sympy as sp
from scipy import integrate

OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate", "")
CONTROLS = {
    "sobolev": ("uses 21/20 of the sharp Sobolev constant", ["sobolev_const", "contraction_sharp"]),
    "gamma": ("scales gamma^1 by 11/10", ["clifford", "isometry", "weitzenbock"]),
    "Jsign": ("flips the sign of J.xi in the identity", ["weitzenbock", "weitzenbock_flip"]),
    "musign": ("flips the sign of mu in the identity", ["weitzenbock"]),
    "vc-angle": ("uses 2 pi for the angular integral of |cos theta|^3, not pi", ["v_c"]),
}
if MUT and MUT not in CONTROLS:
    sys.exit(f"unknown mutation {MUT}")
RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


S_SHARP = 3 * (np.pi / 2) ** (4 / 3)  # Aubin-Talenti, n = 3
S_USED = S_SHARP * (21 / 20 if MUT == "sobolev" else 1)


# ---- Gamma matrices of mass_bounds.md, "Conventions" -------------------------------------------------------------
def gammas(exact=False):
    I, s = (sp.I, [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]) \
        if exact else (1j, [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])])
    if exact:
        z2 = sp.zeros(2)
        g0 = I * sp.diag(1, 1, -1, -1)
        gj = [sp.Matrix(sp.BlockMatrix([[z2, sj], [sj, z2]])) for sj in s]
        if MUT == "gamma":
            gj[0] = gj[0] * sp.Rational(11, 10)
    else:
        z2 = np.zeros((2, 2))
        g0 = I * np.diag([1, 1, -1, -1]).astype(complex)
        gj = [np.block([[z2, sj], [sj, z2]]).astype(complex) for sj in s]
        if MUT == "gamma":
            gj[0] = gj[0] * 1.1
    return g0, gj


def section_algebra():
    print("\n== 1. Clifford algebra, the reduction K_ij g^i g^j = tr K, and the spin connection's sign ==")
    g0, gj = gammas(exact=True)
    G = [g0] + gj
    eta = sp.diag(-1, 1, 1, 1)
    anti = all(sp.simplify(G[a] * G[b] + G[b] * G[a] - 2 * eta[a, b] * sp.eye(4)) == sp.zeros(4)
               for a in range(4) for b in range(4))
    herm = g0.H == -g0 and all(g.H == g for g in gj) and g0.H * g0 == sp.eye(4)
    al = [g0 * g for g in gj]
    alpha_ok = all(a_.H == a_ for a_ in al) and all(
        sp.simplify(al[i] * al[j] + al[j] * al[i] - 2 * int(i == j) * sp.eye(4)) == sp.zeros(4)
        for i in range(3) for j in range(3))
    k = sp.symbols("k0:6", real=True)
    K = sp.Matrix([[k[0], k[1], k[2]], [k[1], k[3], k[4]], [k[2], k[4], k[5]]])
    red = sp.expand(sum((K[i, j] * gj[i] * gj[j] for i in range(3) for j in range(3)), sp.zeros(4))
                    - K.trace() * sp.eye(4))
    check("clifford", anti and herm and alpha_ok and red == sp.zeros(4),
          "{g^a, g^b} = 2 eta^ab, g0 anti-Hermitian and unitary, g^j Hermitian, a_j = g0 g^j Hermitian with "
          "{a_i, a_j} = 2 delta_ij (so |xi| <= |eps|^2), K_ij g^i g^j = tr K for symmetric K")
    # Spin connection Gamma = (1/4) omega_ab g^a g^b. Clifford multiplication is parallel iff
    # [Gamma, g^c] = -omega^c_a g^a, for omega_ab antisymmetric and d theta^a = -omega^a_b ^ theta^b.
    w = sp.symbols("w0:6", real=True)
    om = sp.Matrix([[0, w[0], w[1], w[2]], [-w[0], 0, w[3], w[4]], [-w[1], -w[3], 0, w[5]],
                    [-w[2], -w[4], -w[5], 0]])  # omega_ab
    Gam = sum((om[a, b] * G[a] * G[b] for a in range(4) for b in range(4)), sp.zeros(4)) / 4
    comp = all(sp.expand(Gam * G[c] - G[c] * Gam + sum((eta[c, c] * om[c, a] * G[a] for a in range(4)),
                                                       sp.zeros(4))) == sp.zeros(4) for c in range(4))
    # On the slice, omega_0j(e_i) = K_ij (review.md, from d theta^1 = -d_j b theta^0 ^ theta^j).
    om_i = sp.zeros(4)
    i = 0  # any spatial direction e_i; K's row i
    for j in range(3):
        om_i[0, j + 1], om_i[j + 1, 0] = K[i, j], -K[i, j]
    Gam_i = sp.expand(sum((om_i[a, b] * G[a] * G[b] for a in range(4) for b in range(4)), sp.zeros(4)) / 4)
    minus = sp.expand(Gam_i + sum((K[i, j] * gj[j] * g0 for j in range(3)), sp.zeros(4)) / 2) == sp.zeros(4)
    check("spin_restriction", comp and minus,
          "Gamma = (1/4) omega_ab g^a g^b keeps Clifford multiplication parallel, and on the slice it is "
          "Gamma_i = -(1/2) K_ij g^j g^0 in these gammas (mass_bounds.md writes +; see review.md)")


# ---- 2. Sharp Sobolev constant and the contraction constant -------------------------------------------------------
def radial(fn, a=0, b=np.inf):
    return integrate.quad(fn, a, b, limit=400, epsabs=0, epsrel=1e-13)[0]


def section_sobolev():
    print("\n== 2. Sobolev constant S = 3 (pi/2)^(4/3) and the one-step contraction bound ==")
    ratios = {}
    sech = lambda r: 2 * np.exp(-r) / (1 + np.exp(-2 * r))  # no overflow at large r
    trial = {  # u(r), u'(r)
        "Aubin-Talenti (1 + r^2)^(-1/2)": (lambda r: (1 + r * r) ** -0.5, lambda r: -r * (1 + r * r) ** -1.5),
        "Gaussian exp(-r^2)": (lambda r: np.exp(-r * r), lambda r: -2 * r * np.exp(-r * r)),
        "(1 + r^2)^(-1)": (lambda r: 1 / (1 + r * r), lambda r: -2 * r / (1 + r * r) ** 2),
        "sech(r)": (sech, lambda r: -np.tanh(r) * sech(r)),
    }
    for name, (u, du) in trial.items():
        grad2 = 4 * np.pi * radial(lambda r: du(r) ** 2 * r * r)
        l6 = (4 * np.pi * radial(lambda r: u(r) ** 6 * r * r)) ** (1 / 3)
        ratios[name] = grad2 / l6
        print("  %-32s ||grad u||^2 / ||u||_6^2 = %.12f" % (name, ratios[name]))
    at = ratios["Aubin-Talenti (1 + r^2)^(-1/2)"]
    others = [v for k, v in ratios.items() if not k.startswith("Aubin")]
    check("sobolev_const", abs(at / S_USED - 1) < 1e-10 and min(others) > S_USED,
          "extremal ratio %.12f vs S = %.12f; other trial functions give more" % (at, S_USED))
    # Hoelder is an equality for |tr K| proportional to |psi|^2, and Sobolev for the Aubin-Talenti bubble. So for
    # psi = u chi and tr K = c u^2, (1/2) ||tr K psi||_2 = kappa_K ||grad psi||_2 with kappa_K = ||tr K||_3 / (2 sqrt S):
    # the one-step bound of item 5 is attained and cannot be improved.
    u = trial["Aubin-Talenti (1 + r^2)^(-1/2)"][0]
    du = trial["Aubin-Talenti (1 + r^2)^(-1/2)"][1]
    c = 0.37
    lhs = 0.5 * np.sqrt(4 * np.pi * radial(lambda r: (c * u(r) ** 2 * u(r)) ** 2 * r * r))
    trk3 = (4 * np.pi * radial(lambda r: abs(c * u(r) ** 2) ** 3 * r * r)) ** (1 / 3)
    rhs = trk3 / (2 * np.sqrt(S_USED)) * np.sqrt(4 * np.pi * radial(lambda r: du(r) ** 2 * r * r))
    check("contraction_sharp", lhs <= rhs * (1 + 1e-10) and abs(lhs / rhs - 1) < 1e-10,
          "(1/2) ||tr K psi||_2 = %.12f, kappa_K ||grad psi||_2 = %.12f (equality expected)" % (lhs, rhs))


# ---- 3. The tanh profile ------------------------------------------------------------------------------------------
def profile(sig, R=1.0):
    T = np.tanh(sig * R)
    def f(r):
        return (np.tanh(sig * (r + R)) - np.tanh(sig * (r - R))) / (2 * T)
    def f1(r):
        return sig * (np.cosh(sig * (r + R)) ** -2 - np.cosh(sig * (r - R)) ** -2) / (2 * T)
    def f2(r):
        a, b = sig * (r + R), sig * (r - R)
        return sig**2 * (-2 * np.tanh(a) * np.cosh(a) ** -2 + 2 * np.tanh(b) * np.cosh(b) ** -2) / (2 * T)
    return f, f1, f2


def section_vc():
    print("\n== 3. v_c where kappa_K = 1, tanh profile, R = 1 ==")
    # tr K = d_x beta_x = -v f'(r) cos(theta), so ||tr K||_3^3 = v^3 (int |cos|^3 dOmega) int |f'|^3 r^2 dr.
    ang = 2 * np.pi if MUT == "vc-angle" else np.pi  # int |cos theta|^3 dOmega = 2 pi * 2 * (1/4) = pi
    note = {2: 3.745, 4: 2.471, 8: 1.562, 16: 0.985}  # mass_bounds.md item 5
    worst = 0
    for sR, want in note.items():
        _, f1, _ = profile(sR)
        I3 = sum(radial(lambda r: abs(f1(r)) ** 3 * r * r, a, b) for a, b in ((0, 1), (1, 3), (3, 25)))
        vc = 2 * np.sqrt(S_SHARP) / (ang * I3) ** (1 / 3)
        thin = 2 * np.sqrt(S_SHARP) * (15 / (2 * np.pi * sR**2)) ** (1 / 3)
        worst = max(worst, abs(round(vc, 3) - want))
        print("  sigma R = %2d: v_c = %.5f (note %.3f), thin-wall formula %.5f" % (sR, vc, want, thin))
    check("v_c", worst < 1e-9, "every v_c agrees with the note to its printed 3 decimals (worst %.0e)" % worst)


# ---- 4. Integrated Weitzenbock identity on a grid -----------------------------------------------------------------
def section_identity():
    print("\n== 4. Integrated Weitzenbock identity on the Alcubierre slice, analytic derivatives ==")
    g0, gj = gammas()
    rng = np.random.default_rng(59)
    eps0 = rng.normal(size=4) + 1j * rng.normal(size=4)
    eps0 /= np.linalg.norm(eps0)
    chis = [(rng.normal(size=4) + 1j * rng.normal(size=4)) * 0.6 for _ in range(2)]
    centres = [np.array([0.9, 0.3, -0.2]), np.array([-0.5, -0.6, 0.7])]
    widths = [0.55, 0.45]
    waves = [np.array([0.7, -0.4, 0.5]), np.array([-0.3, 0.9, 0.2])]
    rows = []
    for (sR, v, L, Ns) in ((2, 1.0, 8.0, (80, 112)), (4, 0.5, 5.0, (80, 112))):
        f, f1, f2 = profile(sR)
        for N in Ns:
            h = 2 * L / N
            ax = -L + h * (np.arange(N) + 0.5)  # cell centres, so r = 0 is not a node
            X = np.meshgrid(ax, ax, ax, indexing="ij")
            r = np.sqrt(X[0] ** 2 + X[1] ** 2 + X[2] ** 2)
            F1, F2 = f1(r), f2(r)
            n_ = [Xi / r for Xi in X]
            # K_ij = (d_i beta_j + d_j beta_i) / 2, beta_x = -v f: K_xx = -v f' n_x, K_xA = -v f' n_A / 2.
            K = [[np.zeros_like(r) for _ in range(3)] for _ in range(3)]
            K[0][0] = -v * F1 * n_[0]
            for A in (1, 2):
                K[0][A] = K[A][0] = -0.5 * v * F1 * n_[A]
            trK = K[0][0]
            # item 1: mu = -(v^2/32 pi) (rho/r)^2 f'^2, J_i = (v/16 pi)(d_x d_i f - delta_ix lap f)
            mu = -(v**2 / (32 * np.pi)) * (n_[1] ** 2 + n_[2] ** 2) * F1**2
            hess = lambda i, j: F2 * n_[i] * n_[j] + F1 * (float(i == j) - n_[i] * n_[j]) / r
            lap = F2 + 2 * F1 / r
            Jv = [(v / (16 * np.pi)) * (hess(0, i) - (i == 0) * lap) for i in range(3)]
            if MUT == "musign":
                mu = -mu
            # eps = eps0 + sum_k g_k chi_k, g_k = exp(-|x - c_k|^2 / (2 w_k^2) + i q_k.x)
            eps = np.broadcast_to(eps0[:, None, None, None], (4, N, N, N)).astype(complex).copy()
            deps = [np.zeros((4, N, N, N), complex) for _ in range(3)]
            for chi, c, wd, q in zip(chis, centres, widths, waves):
                d = [X[i] - c[i] for i in range(3)]
                g = np.exp(-(d[0] ** 2 + d[1] ** 2 + d[2] ** 2) / (2 * wd**2) + 1j * (q[0] * X[0] + q[1] * X[1]
                                                                                    + q[2] * X[2]))
                eps += chi[:, None, None, None] * g
                for i in range(3):
                    deps[i] += chi[:, None, None, None] * ((1j * q[i] - d[i] / wd**2) * g)
            mat = lambda M, s: np.einsum("ab,bxyz->axyz", M, s)
            g0e = mat(g0, eps)
            gjg0e = [mat(gj[j], g0e) for j in range(3)]
            nab = [deps[i] + 0.5 * sum(K[i][j] * gjg0e[j] for j in range(3)) for i in range(3)]
            Dpsi = sum(mat(gj[i], deps[i]) for i in range(3))
            Dt = Dpsi + 0.5 * trK * g0e
            dV = h**3
            n2 = lambda s: np.sum(np.abs(s) ** 2) * dV
            G = sum(n2(nab[i]) for i in range(3))
            Dt2 = n2(Dt)
            xi = [np.real(np.einsum("axyz,axyz->xyz", eps.conj(), mat(g0 @ gj[j], eps))) for j in range(3)]
            Jxi = sum(Jv[j] * xi[j] for j in range(3))
            sgn = 1 if MUT == "Jsign" else -1  # item 6 has mu |eps|^2 - J.xi
            src = np.sum(mu * np.sum(np.abs(eps) ** 2, axis=0) + sgn * Jxi) * dV
            resid = G - Dt2 + 4 * np.pi * src
            # The other sign of the K term, the slice's 4D spin connection in these gammas, gives mu|eps|^2 + J.xi.
            Gf = sum(n2(deps[i] - 0.5 * sum(K[i][j] * gjg0e[j] for j in range(3))) for i in range(3))
            srcf = np.sum(mu * np.sum(np.abs(eps) ** 2, axis=0) - sgn * Jxi) * dV
            flip = (Gf - n2(Dpsi - 0.5 * trK * g0e) + 4 * np.pi * srcf) / Gf
            print("    with -K in nabla~ and D~: int|nabla~ eps|^2 - int|D~ eps|^2 + 4 pi int(mu|eps|^2 + J.xi), "
                  "relative %.1e" % flip)
            iso = (n2(Dpsi) - sum(n2(deps[i]) for i in range(3))) / sum(n2(deps[i]) for i in range(3))
            # one-step bound of item 5 on this psi: (1/2) ||tr K psi||_2 <= kappa_K ||grad psi||_2
            psi = eps - eps0[:, None, None, None]
            lhs = 0.5 * np.sqrt(n2(trK * psi))
            kap = (np.sum(np.abs(trK) ** 3) * dV) ** (1 / 3) / (2 * np.sqrt(S_SHARP))
            rhs = kap * np.sqrt(sum(n2(deps[i]) for i in range(3)))
            rows.append((sR, v, N, G, Dt2, 4 * np.pi * src, resid / G, iso, np.sum(Jxi) * dV, lhs / rhs, flip))
            print("  sigma R = %d, v = %.1f, L = %.0f, N = %3d: int|nabla~ eps|^2 = %.10f, int|D~ eps|^2 = %.10f, "
                  "4 pi int(mu|eps|^2 - J.xi) = %+.10f, rel. resid %.1e; int J.xi = %+.3e; "
                  "(||D psi||^2 - ||grad psi||^2)/||grad psi||^2 = %.1e; one-step bound ratio %.3f"
                  % (sR, v, L, N, G, Dt2, 4 * np.pi * src, resid / G, np.sum(Jxi) * dV, iso, lhs / rhs), flush=True)
            del X, r, K, eps, deps, nab, gjg0e, g0e, Dt, Dpsi, xi, psi
    fine = [row for row in rows if (row[0], row[2]) in ((2, 112), (4, 112))]
    check("isometry", all(abs(row[7]) < 1e-9 for row in rows),
          "||D psi||_2 = ||grad psi||_2 on every grid (worst %.1e)" % max(abs(row[7]) for row in rows))
    check("weitzenbock", all(abs(row[6]) < 1e-8 for row in fine) and all(abs(row[8]) > 1e-3 for row in fine),
          "int|nabla~ eps|^2 - int|D~ eps|^2 + 4 pi int(mu|eps|^2 - J.xi) = 0 at the finest grids, relative %s; "
          "int J.xi is not 0, so the J term's sign is tested" % ["%.1e" % row[6] for row in fine])
    check("weitzenbock_flip", all(abs(row[10]) < 1e-8 for row in fine),
          "with -K in nabla~ and D~ the same integrand holds with + J.xi, relative %s"
          % ["%.1e" % row[10] for row in fine])
    check("one_step_bound", all(row[9] <= 1 for row in rows),
          "(1/2) ||tr K psi|| / (kappa_K ||grad psi||) = %s <= 1" % ["%.3f" % row[9] for row in rows])


def main():
    t0 = time.time()
    print("review_spinor.py, issue #59. numpy %s, scipy %s, sympy %s, python %s. mutation: %s"
          % (np.__version__, __import__("scipy").__version__, sp.__version__, sys.version.split()[0], MUT or "none"))
    section_algebra()
    section_sobolev()
    section_vc()
    section_identity()
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

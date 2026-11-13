#!/usr/bin/env python3
"""Independent check of #57's identities (I1) and (I2), for the review of issue #59.

The Ricci tensor is computed from Cartan's structure equations in an orthonormal coframe, by code written for this
review. It shares no code with phantom_source.py or alcubierre_bounds.py, only sympy itself. Functions of the
coordinates are free jet symbols with a total derivative, so every identity is exact in the jets and the point.

    python3 notes/warp/review_cartan.py                  checks, then every control
    python3 notes/warp/review_cartan.py --no-controls    checks only
    python3 notes/warp/review_cartan.py --mutate=NAME    one control (its target checks must fail, exit 1)

Exit status is 1 if a check fails or a control does not fire.
"""
import subprocess
import sys
import time

import sympy as sp

OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate", "")
CONTROLS = {
    "cartan-sign": ("flips the sign of omega ^ omega in the curvature 2-form", ["planted_desitter", "planted_schw"]),
    "torsion": ("solves for omega with T_bca and T_cab swapped", ["structure_eqs"]),
    "rotate": ("rotates (e_rho, e_varphi) by alpha with tan(alpha/2) = 1/20", ["axi_I1", "axi_I2"]),
    "nonaxi": ("makes the shift depend on y^2 + 2 z^2, not on rho", ["axi_I1", "axi_I2"]),
    "lapse": ("gives the slicing the lapse 1 + b, not 1", ["general_I2", "axi_I2"]),
    "table-sign": ("flips the sign of R_0rho in the note's table", ["alc_table"]),
    "wall-v2": ("doubles v in the expected wall curvature", ["alc_wall"]),
}
if MUT and MUT not in CONTROLS:
    sys.exit(f"unknown mutation {MUT}")

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def zero(e):
    if isinstance(e, sp.MatrixBase):
        return all(zero(c) for c in e)
    return sp.simplify(sp.expand(e)) == 0


# ---------------------------------------------------------------------------------------------------------------
# Jet spaces. A jet symbol stands for one partial derivative of an unknown function. D(g, mu) is the total
# derivative of an expression g in the coordinates and the jets.
# ---------------------------------------------------------------------------------------------------------------
class Jets:
    def __init__(self, coords):
        self.coords = coords
        self.rules = []  # each: (predicate on symbol, function symbol -> list of d/dx^mu)

    def total(self, g, mu):
        out = sp.diff(g, self.coords[mu])
        for s in g.free_symbols:
            for owns, deriv in self.rules:
                if owns(s):
                    out += sp.diff(g, s) * deriv(s)[mu]
        return out


def general_shift(jets):
    """b(t, x, y, z) with no symmetry. Symbol b_ijkl is d_t^i d_x^j d_y^k d_z^l b."""
    t, x, y, z = jets.coords
    cache = {}

    def sym(idx):
        if idx not in cache:
            cache[idx] = sp.Symbol("b_%d%d%d%d" % idx, real=True)
            cache[cache[idx]] = idx
        return cache[idx]

    def deriv(s):
        i = cache[s]
        return [sym(tuple(i[k] + (k == m) for k in range(4))) for m in range(4)]

    jets.rules.append((lambda s: s in cache, deriv))
    return sym((0, 0, 0, 0)), sym


def axisymmetric_shift(jets, weight_z=1):
    """b = G(t, x, q), q = y^2 + weight_z z^2. Symbol G_ijk is d_t^i d_x^j d_q^k G. weight_z = 1 is axisymmetric."""
    t, x, y, z = jets.coords
    cache = {}

    def sym(idx):
        if idx not in cache:
            cache[idx] = sp.Symbol("G_%d%d%d" % idx, real=True)
            cache[cache[idx]] = idx
        return cache[idx]

    def deriv(s):
        i, j, k = cache[s]
        up = sym((i, j, k + 1))
        return [sym((i + 1, j, k)), sym((i, j + 1, k)), 2 * y * up, 2 * weight_z * z * up]

    jets.rules.append((lambda s: s in cache, deriv))
    return sym((0, 0, 0)), sym


def alcubierre_shift(jets, v):
    """b = -v H(q), q = (x - v t)^2 + y^2 + z^2, so f(r_s) = H(r_s^2). Symbol H_k is the k-th derivative of H."""
    t, x, y, z = jets.coords
    cache = {}

    def sym(k):
        if k not in cache:
            cache[k] = sp.Symbol("H_%d" % k, real=True)
            cache[cache[k]] = k
        return cache[k]

    xi = x - v * t

    def deriv(s):
        up = sym(cache[s] + 1)
        return [-2 * v * xi * up, 2 * xi * up, 2 * y * up, 2 * z * up]

    jets.rules.append((lambda s: s in cache, deriv))
    return -v * sym(0), sym


# ---------------------------------------------------------------------------------------------------------------
# Cartan's structure equations. Coframe theta^a = E[a, mu] dx^mu, orthonormal for eta = diag(-1, 1, 1, 1).
# A 2-form alpha = (1/2) alpha_cd theta^c ^ theta^d is stored by its frame components alpha_cd.
#   d theta^a = -omega^a_b ^ theta^b,  omega_ab = -omega_ba,  omega^a_b = omega^a_bc theta^c,
#   Omega^a_b = d omega^a_b + omega^a_c ^ omega^c_b = (1/2) R^a_bcd theta^c ^ theta^d,  R_bd = R^a_bad.
# ---------------------------------------------------------------------------------------------------------------
ETA = sp.diag(-1, 1, 1, 1)


def cartan_ricci(E, D, simp=sp.expand):
    n = 4
    F = sp.simplify(E.inv())  # F[mu, a]: e_a = F[mu, a] d_mu
    assert sp.simplify(E * F - sp.eye(n)) == sp.zeros(n)

    def e(a, g):  # frame derivative e_a(g)
        return simp(sum(F[m, a] * D(g, m) for m in range(n)))

    # (d theta^a)_bc = F^mu_b F^nu_c (d_mu E^a_nu - d_nu E^a_mu)
    dE = [[[D(E[a, nu], mu) for nu in range(n)] for mu in range(n)] for a in range(n)]
    dth = [[[simp(sum(F[m, b] * F[k, c] * (dE[a][m][k] - dE[a][k][m]) for m in range(n) for k in range(n)))
             for c in range(n)] for b in range(n)] for a in range(n)]
    # T_acd = eta_ae (d theta^e)_cd = omega_acd - omega_adc, solved by omega_abc = (T_abc + T_bca - T_cab) / 2
    T = [[[ETA[a, a] * dth[a][c][d] for d in range(n)] for c in range(n)] for a in range(n)]
    if MUT == "torsion":
        low = [[[(T[a][b][c] + T[c][a][b] - T[b][c][a]) / 2 for c in range(n)] for b in range(n)] for a in range(n)]
    else:
        low = [[[(T[a][b][c] + T[b][c][a] - T[c][a][b]) / 2 for c in range(n)] for b in range(n)] for a in range(n)]
    om = [[[simp(ETA[a, a] * low[a][b][c]) for c in range(n)] for b in range(n)] for a in range(n)]  # omega^a_bc

    # The structure equations themselves: metric compatibility and zero torsion.
    bad = []
    for a in range(n):
        for b in range(n):
            for c in range(n):
                if not zero(low[a][b][c] + low[b][a][c]):
                    bad.append(("antisym", a, b, c))
                if not zero(dth[a][b][c] - (om[a][b][c] - om[a][c][b])):
                    bad.append(("torsion", a, b, c))

    sgn = -1 if MUT == "cartan-sign" else 1
    Riem = {}
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    dom = e(c, om[a][b][d]) - e(d, om[a][b][c]) + sum(om[a][b][k] * dth[k][c][d] for k in range(n))
                    ww = sum(om[a][k][c] * om[k][b][d] - om[a][k][d] * om[k][b][c] for k in range(n))
                    Riem[a, b, c, d] = dom + sgn * ww
    Ric = sp.Matrix(n, n, lambda b, d: simp(sum(Riem[a, b, a, d] for a in range(n))))
    return Ric, bad, e


def main():
    t0 = time.time()
    print("review_cartan.py, issue #59. sympy %s, python %s. mutation: %s"
          % (sp.__version__, sys.version.split()[0], MUT or "none"))

    # ---- 1. Planted truth: metrics whose Ricci tensor is known --------------------------------------------------
    print("\n== 1. Planted truth for the Cartan routine ==")
    tt, r, th, ph, H, m, a, l = sp.symbols("t r theta phi H m a l", positive=True)
    plain = Jets((tt, r, th, ph))
    E = sp.diag(1, sp.exp(H * tt), sp.exp(H * tt), sp.exp(H * tt))  # de Sitter, flat slicing, coordinates (t, x, y, z)
    Ric, bad, _ = cartan_ricci(E, plain.total)
    check("planted_desitter", zero(Ric - 3 * H**2 * ETA) and not bad, "R_ab = 3 H^2 eta_ab: %s" % list(Ric.diagonal()))
    N2 = 1 - 2 * m / r
    E = sp.diag(sp.sqrt(N2), 1 / sp.sqrt(N2), r, r * sp.sin(th))
    Ric, bad, _ = cartan_ricci(E, plain.total, simp=sp.simplify)
    check("planted_schw", Ric.applyfunc(sp.simplify) == sp.zeros(4) and not bad, "Schwarzschild R_ab = 0")
    rr = sp.sqrt(r**2 + a**2)  # Ellis-Bronnikov with l -> r as the coordinate name
    E = sp.diag(1, 1, rr, rr * sp.sin(th))
    Ric, bad, _ = cartan_ricci(E, plain.total, simp=sp.simplify)
    want = sp.diag(0, -2 * a**2 / (r**2 + a**2) ** 2, 0, 0)
    check("planted_eb", (Ric - want).applyfunc(sp.simplify) == sp.zeros(4) and not bad,
          "Ellis-Bronnikov R_ll = %s, others 0" % sp.simplify(Ric[1, 1]))

    # ---- 2. Unit lapse, flat slices, shift b(t, x, y, z) d_x, no symmetry ----------------------------------------
    print("\n== 2. General shift: ds^2 = -dt^2 + (dx + b dt)^2 + dy^2 + dz^2, frame (n, d_x, d_y, d_z) ==")
    t, x, y, z = sp.symbols("t x y z", real=True)
    J = Jets((t, x, y, z))
    b, B = general_shift(J)
    lapse = 1 + b if MUT == "lapse" else 1
    E = sp.Matrix([[lapse, 0, 0, 0], [b, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    Ric, bad, e = cartan_ricci(E, J.total)
    check("structure_eqs", not bad, "omega antisymmetric and torsion-free: %d violations" % len(bad))
    # The 3+1 formulas (hand derivation in review.md), K_ij = (d_i b_j + d_j b_i)/2 with b_j = delta_jx b.
    db = [B((0, 1, 0, 0)), B((0, 0, 1, 0)), B((0, 0, 0, 1))]  # d_x b, d_y b, d_z b
    K = sp.Matrix(3, 3, lambda i, j: (int(j == 0) * db[i] + int(i == 0) * db[j]) / 2)
    trK = K.trace()
    KK = sum(K[i, j] ** 2 for i in range(3) for j in range(3))
    nK = lambda g: e(0, g)  # n = e_0 = d_t - b d_x
    want00 = nK(trK) - KK
    want0 = [J.total(trK, i + 1) - sum(J.total(K[j, i], j + 1) for j in range(3)) for i in range(3)]
    K2 = K * K
    wantij = sp.Matrix(3, 3, lambda i, j: -nK(K[i, j]) + K[0, j] * db[i] + K[i, 0] * db[j] + trK * K[i, j]
                       - 2 * K2[i, j])
    ok3 = zero(Ric[0, 0] - want00) and all(zero(Ric[0, i + 1] - want0[i]) for i in range(3)) and \
        all(zero(Ric[i + 1, j + 1] - wantij[i, j]) for i in range(3) for j in range(3))
    check("general_3p1", ok3, "R_00 = n(K) - K_ij K_ij, R_0i = d_i K - d_j K_ji, "
          "R_ij = -n(K_ij) + K_jx d_i b + K_ix d_j b + K K_ij - 2 K_ik K_kj")
    perp = db[1] ** 2 + db[2] ** 2
    i2gen = sp.expand(Ric[0, 0] + Ric[1, 1] + 2 * (Ric[2, 2] + Ric[3, 3]))
    rab = [sp.expand(Ric[A, Bq] - db[A - 1] * db[Bq - 1] / 2) for A in (2, 3) for Bq in (2, 3)]
    g00 = sp.expand((Ric[0, 0] + Ric[1, 1] + Ric[2, 2] + Ric[3, 3]) / 2 + perp / 4)
    check("general_I2", i2gen == 0 and rab == [0] * 4 and g00 == 0,
          "R_00 + R_11 + 2 (R_22 + R_33) = %s, R_AB - b_A b_B / 2 = %s, G_00 + |grad_perp b|^2 / 4 = %s"
          % (i2gen, rab, g00))
    print("  R_00 + R_11 =", sp.factor(sp.expand(Ric[0, 0] + Ric[1, 1])))

    # ---- 3. Axisymmetric shift b = G(t, x, y^2 + z^2): the frame (n, d_x, e_rho, e_varphi) -------------------
    print("\n== 3. Axisymmetric shift b = G(t, x, y^2 + z^2), any t-dependence, frame (n, d_x, e_rho, e_varphi) ==")
    J = Jets((t, x, y, z))
    b, Gs = axisymmetric_shift(J, weight_z=2 if MUT == "nonaxi" else 1)
    lapse = 1 + b if MUT == "lapse" else 1
    E = sp.Matrix([[lapse, 0, 0, 0], [b, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    Ric, bad, e = cartan_ricci(E, J.total)
    # rho e_rho = (y, z), rho e_varphi = (-z, y) in the (d_y, d_z) plane. Rotation by alpha for the control.
    if MUT == "rotate":
        tau = sp.Rational(1, 20)
        ca, sa = (1 - tau**2) / (1 + tau**2), 2 * tau / (1 + tau**2)
    else:
        ca, sa = 1, 0
    P = sp.Matrix([[y, z], [-z, y]])  # rows: rho e_rho, rho e_varphi
    P = sp.Matrix([[ca, sa], [-sa, ca]]) * P
    Rt = Ric[2:, 2:]
    Rrr, Rrp, Rpp = [sp.expand((P[i, :] * Rt * P[j, :].T)[0]) for i, j in ((0, 0), (0, 1), (1, 1))]
    Rap = [sp.expand(P[1, 0] * Ric[c, 2] + P[1, 1] * Ric[c, 3]) for c in (0, 1)]
    rho2 = y**2 + z**2
    check("axi_I1", Rap == [0, 0] and Rrp == 0 and Rpp == 0,
          "rho R(e_varphi, n), rho R(e_varphi, d_x), rho^2 R_rho,varphi, rho^2 R_varphi,varphi = %s"
          % ([Rap[0], Rap[1], Rrp, Rpp] if MUT else "[0, 0, 0, 0]"))
    i2 = sp.expand(rho2 * (Ric[0, 0] + Ric[1, 1]) + 2 * Rrr)
    check("axi_I2", i2 == 0, "rho^2 (R_00 + R_11 + 2 R_rho,rho) = %s" % (sp.factor(i2) if i2 != 0 else 0))
    axis = {y: 0, z: 0}
    on = [sp.expand(Ric[A, c].subs(axis)) for A in (2, 3) for c in range(4)] + [sp.expand((Ric[0, 0] + Ric[1, 1])
                                                                                          .subs(axis))]
    check("axi_axis", on == [0] * 9, "on y = z = 0: R(E_2, .), R(E_3, .), R_00 + R_11 = %s" % on)

    # ---- 4. Alcubierre, constant v, f(r_s) = H(r_s^2): the note's table and wall curvature ------------------------
    print("\n== 4. Alcubierre with constant v: the note's frame components and wall curvature ==")
    v = sp.Symbol("v", positive=True)
    J = Jets((t, x, y, z))
    b, Hs = alcubierre_shift(J, v)
    E = sp.Matrix([[1, 0, 0, 0], [b, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    Ric, bad, e = cartan_ricci(E, J.total)
    xi = x - v * t
    rs = sp.sqrt(xi**2 + y**2 + z**2)
    H0, H1, H2 = Hs(0), Hs(1), Hs(2)
    # f(r) = H(r^2): f = H0, f' = 2 r H1, f'' = 2 H1 + 4 r^2 H2. w = v (1 - f).
    f1 = 2 * rs * H1
    f2 = 2 * H1 + 4 * rs**2 * H2
    w, w1, w2 = v * (1 - H0), -v * f1, -v * f2
    s, c = sp.sqrt(rho2) / rs, xi / rs
    note = {  # phantom_source.md, "The Ricci tensor in the frame"
        "00": -(2 * w * w1 * s**2 + 2 * rs * w * w2 * c**2 + rs * w1**2 * (1 + c**2)) / (2 * rs),
        "01": -(w1 * (1 + c**2) + rs * w2 * s**2) / (2 * rs),
        "0r": (rs * w2 - w1) * s * c / (2 * rs),
        "11": (2 * w * w1 * s**2 + 2 * rs * w * w2 * c**2 + rs * w1**2 * (3 * c**2 - 1)) / (2 * rs),
        "1r": (w * (rs * w2 - w1) + 2 * rs * w1**2) * s * c / (2 * rs),
        "rr": w1**2 * s**2 / 2,
    }
    if MUT == "table-sign":
        note["0r"] = -note["0r"]
    rho = sp.sqrt(rho2)
    mine = {
        "00": Ric[0, 0], "01": Ric[0, 1], "11": Ric[1, 1],
        "0r": (y * Ric[0, 2] + z * Ric[0, 3]) / rho,
        "1r": (y * Ric[1, 2] + z * Ric[1, 3]) / rho,
        "rr": (y**2 * Ric[2, 2] + 2 * y * z * Ric[2, 3] + z**2 * Ric[3, 3]) / rho2,
    }
    diffs = {k: sp.simplify(mine[k] - note[k]) for k in note}
    check("alc_table", all(d == 0 for d in diffs.values()), "code minus the note's table: %s" % diffs)
    vv = 2 * v if MUT == "wall-v2" else v
    wall_rr = sp.simplify(mine["rr"] - vv**2 * f1**2 * rho2 / (2 * rs**2))
    wall_01 = sp.simplify(Ric[0, 1].subs(axis) - (vv * f1 / rs).subs(axis))
    rho_e = sp.simplify((Ric[0, 0] + Ric[1, 1] + Ric[2, 2] + Ric[3, 3]) / 2 / (8 * sp.pi)
                        + v**2 / (32 * sp.pi) * rho2 / rs**2 * f1**2)
    check("alc_wall", wall_rr == 0 and wall_01 == 0 and rho_e == 0,
          "R_rho,rho - v^2 f'^2 rho^2 / (2 r_s^2) = %s, on the axis R_01 - v f' / r_s = %s, "
          "G_00 / 8 pi + (v^2 / 32 pi) (rho / r_s)^2 f'^2 = %s" % (wall_rr, wall_01, rho_e))

    fails = [n for n, ok in RESULTS if not ok]
    print("\nchecks: %d passed, %d failed %s. %.0f s" % (len(RESULTS) - len(fails), len(fails), fails,
                                                         time.time() - t0))
    if MUT:
        target = CONTROLS[MUT][1]
        fired = all(n in fails for n in target)
        print("control %s: targets %s %s" % (MUT, target, "FIRED" if fired else "SILENT"))
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

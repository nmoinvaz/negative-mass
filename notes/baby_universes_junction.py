#!/usr/bin/env python3
"""Junction algebra of the S-brane of [BHR21] and of the Schwarzschild-de Sitter matching of [FMM90]'s setting (#224).

    python3 baby_universes_junction.py                   all checks; exit 1 if any check fails
    python3 baby_universes_junction.py --mutate=NAME     break one input on purpose (the targets of NAME must fire)
    python3 baby_universes_junction.py --mutate=all      run every mutation in a subprocess, compare with MUTATIONS
    python3 baby_universes_junction.py --list            the mutations and their declared targets

Conventions. Signature (-+++), G_mu_nu = kappa T_mu_nu with kappa = 8 pi G, positive energy density has G^tau_tau = -kappa rho.
The paper [BHR21] uses (+---); the mixed extrinsic curvature K^i_j = h^{ik} (1/2) L_n h_kj with the FUTURE-directed unit normal n
is the same in both signatures (K^x_x = adot/a), and so is every physical statement (NEC, WEC). Shell stress: T^mu_nu =
S^mu_nu delta(proper time), S^i_j = diag(p_x, p_th, p_th), S^n_n = 0. An S-brane of tension sigma, action -sigma int sqrt(h) d^3x,
has p_x = p_th = -sigma. Inside Schwarzschild, f = 1 - r_s/r < 0; F := -f(r0) = r_s/r0 - 1 > 0; the black-hole interior has
r decreasing to the future, n^r = -sqrt(F). Units of the checks are free (r0, F, H, kappa are positive symbols).

Method.
  J1  the singular part of G^mu_nu for ds^2 = -dtau^2 + A(tau)^2 dx^2 + B(tau)^2 dOmega^2 with A, B continuous and A', B' jumping:
      read off with int A'' dtau = [A'] (no Israel formula used). This anchors the sign of the junction condition.
  J2  the same junction condition as the Israel formula  [K^i_j] - delta^i_j [K] = kappa S^i_j, K from the metrics.
  Everything else is derived through israel_S() and compared with by-hand values written into the checks.
Checks: S0a S0b S0c S0d, J1 J1b, J2a J2b J2c, P0a P0b, P1 P2 P3, X1 X2 X3, E1 E2, D2 D2b D3 D3b D4 D5 D6 D6b; the criteria, the by-hand
values and the mutation targets are in negative-mass-raw I224-rpi5/README.md (written before the run of record).
"""
import argparse
import platform
import subprocess
import sys
import time

import sympy as sp

MUT = None
CHECKS = []          # (id, ok, text)
T0 = time.time()

# mutation -> (what it changes, the checks that must fire on it); the same table as the README, section 5
MUTATIONS = {
    "jump_sign": ("[K] = K_- - K_+ in israel_S (the sign of the jump of K flipped)",
                  ["J1", "P0a", "P3", "X3", "E2", "D3", "D6", "D6b"]),
    "trace_term": ("the -delta^i_j [K] term dropped from the Israel formula",
                   ["J1", "P0a", "P0b", "P3", "X3", "E2", "D3", "D6", "D6b"]),
    "tension_sign": ("p = +sigma instead of -sigma for the S-brane",
                     ["P0a", "P3", "X3", "E1"]),
    "orientation": ("Schwarzschild normal n^r = +sqrt(F) (white-hole interior)",
                    ["J2a", "S0c", "P1", "X1", "X2", "X3", "E2", "D3", "D5"]),
    "lapse": ("Schwarzschild normal n^r = -F (square root forgotten)",
              ["J2c", "J2a", "S0c", "S0d", "P0b", "P1", "X1", "X2", "X3", "D3", "D4", "D5", "D6b"]),
    "printed_g": ("g in the printed (4), (5) read as -det(g_mu_nu) of the 4-metric, not as -f(r0)",
                  ["P1", "P2"]),
    "ds_swap": ("de Sitter slicing a = cosh(H tau)/H, b = sinh(H tau)/H (the two scale factors exchanged)",
                ["S0b", "D2"]),
    "ds_orientation": ("de Sitter static form with n^r = -sqrt(-f2) (past-directed normal)",
                       ["D2b"]),
}

tau, th, ph, t_, r_ = sp.symbols("tau theta phi t r", real=True)
r0, F, H, kap = sp.symbols("r0 F H kappa", positive=True)
rs = sp.symbols("r_s", positive=True)
tt = sp.symbols("t", positive=True)       # t = sqrt(c^2-1) = sinh(H tau0), c = H r0 = cosh(H tau0)


def check(cid, ok, text):
    CHECKS.append((cid, bool(ok), text))
    print(f"{cid:5s} {'PASS' if ok else 'FAIL'}  {text}")


# ------------------------------------------------------------------------------------------------ geometry
def einstein_mixed(g, X):
    """G^a_b for the metric g (sympy Matrix) in coordinates X, R^a_bcd = d_c Gam^a_db - ..., MTW sign"""
    n = len(X)
    gi = g.inv()
    Gam = [[[sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(n)) / 2
             for c in range(n)] for b in range(n)] for a in range(n)]

    def ric(b, c):
        return sum(sp.diff(Gam[a][b][c], X[a]) - sp.diff(Gam[a][b][a], X[c])
                   + sum(Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][b][a] for d in range(n)) for a in range(n))

    Rm = sp.Matrix(n, n, lambda b, c: sp.simplify(ric(b, c)))
    Rs = sp.simplify(sum(gi[a, b] * Rm[a, b] for a in range(n) for b in range(n)))
    return (gi * Rm - Rs / 2 * sp.eye(n)).applyfunc(sp.simplify)


def K_lie(g, X, n, tang):
    """K^a_b = (1/2) g^{ac} (L_n g)_{cb} on the tangential indices"""
    N = len(X)
    gi = g.inv()
    L = sp.zeros(N, N)
    for m in range(N):
        for v in range(N):
            L[m, v] = (sum(n[a] * sp.diff(g[m, v], X[a]) for a in range(N))
                       + sum(g[a, v] * sp.diff(n[a], X[m]) for a in range(N))
                       + sum(g[m, a] * sp.diff(n[a], X[v]) for a in range(N)))
    K = gi * L / 2
    return sp.Matrix([[sp.simplify(K[a, b]) for b in tang] for a in tang])


def K_cov(g, X, n, tang):
    """K^a_b = g^{ac} nabla_c n_b on the tangential indices (n_b = g_bd n^d)"""
    N = len(X)
    gi = g.inv()
    Gam = [[[sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(N))
              / 2 for c in range(N)] for b in range(N)] for a in range(N)]
    nl = [sum(g[b, d] * n[d] for d in range(N)) for b in range(N)]
    D = sp.Matrix(N, N, lambda c, b: sp.diff(nl[b], X[c]) - sum(Gam[a][c][b] * nl[a] for a in range(N)))
    K = gi * D
    return sp.Matrix([[sp.simplify(K[a, b]) for b in tang] for a in tang])


# ------------------------------------------------------------------------------------------------ the junction formula
def jump(plus, minus):
    """[X] = X_+ - X_-, + the future side"""
    return (minus - plus) if MUT == "jump_sign" else (plus - minus)


def israel_S(Kp, Km, kappa):
    """S^i_j (list of the three diagonal entries) from K^i_j on the + and - sides: [K^i_j] - delta^i_j [K] = kappa S^i_j"""
    dK = [jump(p, m) for p, m in zip(Kp, Km)]
    tr = sum(dK)
    if MUT == "trace_term":
        return [d / kappa for d in dK]
    return [(d - tr) / kappa for d in dK]


def pressure_of_tension(sigma):
    return sigma if MUT == "tension_sign" else -sigma


def tension_of_pressure(p):
    return p if MUT == "tension_sign" else -p


def rho_from_constraint(K3, R3, kappa):
    """kappa rho = (1/2)(R3 + K^2 - K^i_j K^j_i) for K^i_j = diag(K3), spacelike hypersurface, n timelike"""
    tr = sum(K3)
    return (R3 + tr ** 2 - sum(k ** 2 for k in K3)) / (2 * kappa)


def sphere_R3(radius):
    return 2 / radius ** 2          # R x S^2 of radius `radius`


# ------------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mutate", default=None)
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    global MUT
    if args.list:
        for k, (what, tg) in MUTATIONS.items():
            print(f"{k:13s} {what}\n{'':13s} targets: {' '.join(tg)}")
        return 0
    if args.mutate == "all":
        return run_all()
    MUT = args.mutate
    if MUT and MUT not in MUTATIONS:
        sys.exit(f"unknown mutation {MUT}")
    print(f"# baby_universes_junction.py, python {platform.python_version()}, sympy {sp.__version__}, mutation: {MUT}")

    # ---------------------------------------------------------------------------- backgrounds
    f_r = 1 - rs / r_
    gS = sp.diag(-f_r, 1 / f_r, r_ ** 2, r_ ** 2 * sp.sin(th) ** 2)
    XS = [t_, r_, th, ph]
    GS = einstein_mixed(gS, XS)
    check("S0a", GS == sp.zeros(4, 4), "Schwarzschild: G^mu_nu = 0 (checks the Einstein routine's vacuum case)")

    if MUT == "ds_swap":
        Aexp, Bexp = sp.cosh(H * tau) / H, sp.sinh(H * tau) / H
    else:
        Aexp, Bexp = sp.sinh(H * tau) / H, sp.cosh(H * tau) / H
    gD = sp.diag(-1, Aexp ** 2, Bexp ** 2, Bexp ** 2 * sp.sin(th) ** 2)
    xx = sp.symbols("x", real=True)
    GD = einstein_mixed(gD, [tau, xx, th, ph])
    okD = (GD + 3 * H ** 2 * sp.eye(4)).applyfunc(sp.simplify) == sp.zeros(4, 4)
    check("S0b", okD, "de Sitter in R x S^2 slicing, a = sinh(H tau)/H, b = cosh(H tau)/H: G^mu_nu = -3 H^2 delta^mu_nu"
          " (so kappa rho_Lambda = 3 H^2; also fixes the sign convention G^tau_tau = -kappa rho)")
    gFRW = sp.diag(-1, sp.Function("a")(tau) ** 2, sp.Function("a")(tau) ** 2, sp.Function("a")(tau) ** 2)
    aa = sp.Function("a")(tau)
    GF = einstein_mixed(gFRW, [tau, xx, sp.Symbol("y", real=True), sp.Symbol("z", real=True)])
    check("S0a'", sp.simplify(GF[0, 0] + 3 * (aa.diff(tau) / aa) ** 2) == 0
          and sp.simplify(GF[1, 1] + 2 * aa.diff(tau, 2) / aa + (aa.diff(tau) / aa) ** 2) == 0,
          "flat FRW: G^tau_tau = -3 H^2, G^i_j = -(2 addot/a + H^2) delta^i_j (sign anchor)")

    # ---------------------------------------------------------------------------- J1: delta-function extraction
    Af, Bf = sp.Function("A")(tau), sp.Function("B")(tau)
    gKS = sp.diag(-1, Af ** 2, Bf ** 2, Bf ** 2 * sp.sin(th) ** 2)
    GK = einstein_mixed(gKS, [tau, xx, th, ph])
    print("      KS metric, mixed Einstein tensor: G^tau_tau =", GK[0, 0], ";  G^x_x =", GK[1, 1], ";  G^th_th =", GK[2, 2])
    a2, b2 = sp.symbols("a2 b2")
    GKs = GK.subs({Af.diff(tau, 2): a2, Bf.diff(tau, 2): b2})
    no2 = not any(e.has(sp.Derivative(Af, (tau, 2))) or e.has(sp.Derivative(Bf, (tau, 2))) for e in GKs)
    lin = all(sp.diff(GKs[i, i], v, 2) == 0 for i in range(4) for v in (a2, b2)) and \
        all(sp.diff(GKs[i, i], a2, b2) == 0 for i in range(4))
    A0, B0, kxp, kxm, kthp, kthm = sp.symbols("A0 B0 kxp kxm kthp kthm", positive=True)
    # the coefficient of delta(tau) in A'' is [A'] = A0 (kxp - kxm), in B'' it is B0 (kthp - kthm); no jump() here
    coef = {i: sp.simplify((sp.diff(GKs[i, i], a2) * A0 * (kxp - kxm) + sp.diff(GKs[i, i], b2) * B0 * (kthp - kthm))
                           .subs({Af: A0, Bf: B0})) for i in range(4)}
    no_d = not any(sp.diff(GKs[i, i], v).has(sp.Derivative) for i in range(4) for v in (a2, b2))
    check("J1b", coef[0] == 0 and no2 and lin and no_d,
          "G^tau_tau has no A'', B'' (no delta in the energy density: S^n_n = 0 is forced); G^mu_nu is linear in A'', B''")
    Sx, Sth, _ = israel_S([kxp, kthp, kthp], [kxm, kthm, kthm], kap)
    ok1 = sp.simplify(coef[1] - kap * Sx) == 0 and sp.simplify(coef[2] - kap * Sth) == 0 and sp.simplify(coef[3] - kap * Sth) == 0
    print(f"      extracted: kappa S^x_x = {coef[1]},  kappa S^th_th = {coef[2]}")
    check("J1", ok1, "singular part of G^i_j (int A'' = [A']) equals the Israel formula [K^i_j] - delta^i_j [K] = kappa S^i_j,"
          " i.e. kappa p_x = -2[K^th], kappa p_th = -([K^x] + [K^th])")

    # ---------------------------------------------------------------------------- J2: K from the metrics
    sub0 = {r_: r0, rs: r0 * (F + 1)}
    if MUT == "orientation":
        nr = sp.sqrt(rs / r_ - 1)
    elif MUT == "lapse":
        nr = -(rs / r_ - 1)
    else:
        nr = -sp.sqrt(rs / r_ - 1)
    nS = [0, nr, 0, 0]
    norm = sp.simplify(sum(gS[i, j] * nS[i] * nS[j] for i in range(4) for j in range(4)))
    check("J2c", sp.simplify(norm + 1) == 0, f"Schwarzschild normal n^mu is a unit timelike vector, g(n,n) = {norm} (must be -1)")
    KmL = K_lie(gS, XS, nS, [0, 2, 3]).subs(sub0).applyfunc(sp.simplify)
    KmC = K_cov(gS, XS, nS, [0, 2, 3]).subs(sub0).applyfunc(sp.simplify)
    Kx_m, Kth_m = KmL[0, 0], KmL[1, 1]
    print(f"      K^t_t(-) = {Kx_m},  K^th_th(-) = {Kth_m},  K^ph_ph(-) = {KmL[2, 2]}")
    ref_Kx = (F + 1) / (2 * r0 * sp.sqrt(F))
    ref_Kth = -sp.sqrt(F) / r0
    check("J2a", sp.simplify(Kx_m - ref_Kx) == 0 and sp.simplify(Kth_m - ref_Kth) == 0 and sp.simplify(KmL[2, 2] - ref_Kth) == 0
          and KmL[0, 1] == 0,
          "K^a_b(-) from (1/2) L_n g at r = r0, future normal n^r = -sqrt(F): K^x_x = (F+1)/(2 r0 sqrt F) = r_s/(2 r0^2 sqrt F),"
          " K^th_th = -sqrt(F)/r0")
    # de Sitter side: K from the closed form via the Lie derivative in the Gaussian normal metric
    KpL = K_lie(gD, [tau, xx, th, ph], [1, 0, 0, 0], [1, 2, 3])
    check("J2b", sp.simplify(KmL - KmC) == sp.zeros(3, 3) and sp.simplify(KpL - K_cov(gD, [tau, xx, th, ph], [1, 0, 0, 0], [1, 2, 3]))
          == sp.zeros(3, 3), "K from L_n g equals K from nabla_c n_b (Schwarzschild r = r0 and de Sitter tau = const)")

    # ---------------------------------------------------------------------------- S0c / S0d
    # the Schwarzschild interior as a Kantowski-Sachs cosmology: B = r(tau), dr/dtau = -sqrt(F(r)) (black hole: r decreases),
    # A = sqrt(F(r)); orientation fixed here, independent of the mutable normal
    rr = sp.symbols("rr", positive=True)
    Fr = rs / rr - 1
    rdot = -sp.sqrt(Fr)
    Ad = sp.diff(sp.sqrt(Fr), rr) * rdot
    Add = sp.diff(Ad, rr) * rdot
    Bd = rdot
    Bdd = sp.diff(Bd, rr) * rdot
    subKS = [(Af.diff(tau, 2), Add), (Bf.diff(tau, 2), Bdd), (Af.diff(tau), Ad), (Bf.diff(tau), Bd), (Af, sp.sqrt(Fr)), (Bf, rr)]
    GKS_S = GK.subs(subKS).applyfunc(sp.simplify)
    okvac = GKS_S == sp.zeros(4, 4)
    KxS = sp.simplify((Ad / sp.sqrt(Fr)).subs(sub0).subs(rr, r0))
    KthS = sp.simplify((Bd / rr).subs(sub0).subs(rr, r0))
    check("S0c", okvac and sp.simplify(KxS - Kx_m) == 0 and sp.simplify(KthS - Kth_m) == 0,
          f"Schwarzschild interior as a KS cosmology (A = sqrt(F), B = r, dr/dtau = -sqrt F): vacuum equations hold and"
          f" Adot/A = {KxS}, Bdot/B = {KthS} equal the K of J2a (second computation of K(-))")
    rho_m = rho_from_constraint([Kx_m, Kth_m, Kth_m], sphere_R3(r0), kap)
    check("S0d", sp.simplify(rho_m) == 0,
          "Hamiltonian constraint on the Schwarzschild side: (R3 + K^2 - K^i_j K^j_i)/2 = kappa rho = 0, R3 = 2/r0^2")

    # ---------------------------------------------------------------------------- planted truth
    H1, H2 = sp.symbols("H1 H2", positive=True)
    Sf = israel_S([H2] * 3, [H1] * 3, kap)
    p_hand = -2 * (H2 - H1) / kap                       # by hand: Hdot = -(kappa/2)(rho + p), [H] = -(kappa/2) p_s
    sig_code = tension_of_pressure(Sf[0])
    rho_code = -sum(S * (H1 + H2) / 2 for S in Sf)      # -S^i_j Kbar^j_i
    rho_hand = 3 * (H2 ** 2 - H1 ** 2) / kap            # Friedmann on both sides
    ok_p0 = (sp.simplify(Sf[0] - p_hand) == 0 and sp.simplify(sig_code + p_hand) == 0 and sp.simplify(rho_code - rho_hand) == 0
             and sp.simplify(Sf[0] - Sf[1]) == 0)
    check("P0a", ok_p0, "planted: flat FRW H1 -> H2 across a layer: p_s = -2[H]/kappa (tension 2[H]/kappa), [rho] = 3(H2^2-H1^2)/kappa"
          " from the junction code, against Raychaudhuri and Friedmann by hand")
    F1, F2 = sp.symbols("F1 F2", positive=True)
    KvF = lambda Fv: [KmL[0, 0].subs(F, Fv), KmL[1, 1].subs(F, Fv), KmL[2, 2].subs(F, Fv)]   # K(-) of J2a at another F
    K1v, K2v = KvF(F1), KvF(F2)
    Kx1, Kth1 = K1v[0], K1v[1]
    Kx2, Kth2 = K2v[0], K2v[1]
    S12 = israel_S([Kx2, Kth2, Kth2], [Kx1, Kth1, Kth1], kap)
    Kb = [(Kx1 + Kx2) / 2, (Kth1 + Kth2) / 2, (Kth1 + Kth2) / 2]
    ident = sp.simplify(sum(S * k for S, k in zip(S12, Kb)))
    zero_S = israel_S([Kx1, Kth1, Kth1], [Kx1, Kth1, Kth1], kap)
    check("P0b", ident == 0 and all(z == 0 for z in zero_S),
          "planted: Schwarzschild(F1) -> Schwarzschild(F2) at the same r0: S^i_j Kbar^j_i = [rho] = 0 identically; F2 = F1 gives S = 0")

    # ---------------------------------------------------------------------------- printed (4), (5)
    adot, bdot = sp.symbols("adot_a bdot_b")
    lam = sp.symbols("lambda", positive=True)
    f0 = -F
    fp0 = (r0 * (F + 1)) / r0 ** 2                      # f'(r0) = r_s / r0^2
    g_det = r0 ** 4 * sp.sin(th) ** 2                   # -det g_mu_nu of (2); the determinant itself is negative
    gp = g_det if MUT == "printed_g" else F             # printed g := -f(r0) = F, or the determinant reading
    p4 = adot + sp.sqrt(gp) * fp0 / (2 * f0)
    p5 = bdot + sp.sqrt(gp) / r0
    s_ = sp.Symbol("s", real=True)                      # the jump of each K^i_j in the S-brane case
    d4 = adot - Kx_m                                    # derived: K^x_x(+) - K^x_x(-) = s
    d5 = bdot - Kth_m
    ok_p1 = sp.simplify(p4 - d4) == 0 and sp.simplify(p5 - d5) == 0
    check("P1", ok_p1, "printed (4) adot/a + sqrt(g) f'/(2f) = S and (5) bdot/b + sqrt(g)/r0 = S, with g := -f(r0),"
          " equal the derived K^x_x(+) - K^x_x(-) and K^th_th(+) - K^th_th(-)")

    def scale_exponent(gexpr):      # sqrt(g) ~ lam^k when every length is multiplied by lam
        sq = sp.sqrt(gexpr)
        return sp.simplify(sp.log(sp.simplify(sq.subs(r0, lam * r0) / sq)) / sp.log(lam))
    k_id, k_det = scale_exponent(gp), scale_exponent(g_det)
    check("P2", k_id == 0 and k_det == 2,
          f"dimensions: each printed term must scale as 1/length. sqrt(g) f'/(2f) and sqrt(g)/r0 do so only if sqrt(g) is"
          f" dimensionless: the identification g := -f(r0) has exponent {k_id}, the determinant reading (the text after (1)) {k_det}")
    sig = sp.symbols("sigma", positive=True)
    kj = sp.symbols("kj")
    Siso = [pressure_of_tension(sig)] * 3
    solved = sp.solve(sp.Eq(israel_S([kj] * 3, [0] * 3, kap)[0], Siso[0]), kj)
    kj_val = solved[0]
    ratio = sp.simplify((kj_val.subs({kap: 1, sig: 1})))
    check("P3", ratio == sp.Rational(1, 2),
          f"S-brane of tension sigma, action (1) with 1/2 R (kappa = 1): [K^i_j] = {kj_val} = {ratio} sigma; the printed right-hand"
          " side is S = sigma, a factor 2 above the derived one (finding, see the note)")

    # ---------------------------------------------------------------------------- thresholds for expansion
    s_a = sp.solve(sp.Eq(Kx_m + s_, 0), s_)[0]           # adot/a(+) = K^x(-) + s > 0
    s_b = sp.solve(sp.Eq(Kth_m + s_, 0), s_)[0]
    Km_tr = Kx_m + 2 * Kth_m
    s_v = sp.solve(sp.Eq(Km_tr + 3 * s_, 0), s_)[0]
    ref_sa = -(F + 1) / (2 * r0 * sp.sqrt(F))
    ref_sb = sp.sqrt(F) / r0
    ref_sv = (3 * F - 1) / (6 * r0 * sp.sqrt(F))
    check("X1", sp.simplify(s_a - ref_sa) == 0 and sp.simplify(s_b - ref_sb) == 0 and sp.simplify(s_v - ref_sv) == 0,
          "thresholds of s = kappa sigma/2 for the sign of the future rate: adot/a > 0 iff s > s_a = -(F+1)/(2 r0 sqrt F) < 0;"
          " bdot/b > 0 iff s > s_b = sqrt(F)/r0 > 0; volume iff s > s_v = (3F-1)/(6 r0 sqrt F)")
    sgn = {
        "s_a<0 for all F>0": ref_sa.is_negative,
        "s_b>0 for all F>0": ref_sb.is_positive,
        "s_v>0 iff 3F-1>0 iff r0<3 r_s/4": sp.simplify((3 * F - 1 - (3 * rs / r0 - 4)).subs(rs, r0 * (F + 1))) == 0,
    }
    ok_sgn = all(v is True for v in sgn.values())
    check("X2", ok_sgn and sp.simplify(2 * Km_tr + 3 * ref_sb - 1 / (r0 * sp.sqrt(F))) == 0,
          "signs: s_a < 0, s_b > 0, s_v > 0 iff r0 < 3 r_s/4 (r0 = 3 r_s/4 is the maximum of the Schwarzschild interior volume);"
          " expansion of both scale factors needs s > s_b > 0, a positive tension; and 2K(-) + 3 s_b = 1/(r0 sqrt F)")
    print("      s_a, s_b, s_v =", ref_sa, ",", ref_sb, ",", ref_sv)
    # the same thresholds in the tension sigma, through the junction code: tension -> pressure -> jump of K (as in P3)
    sg_ = sp.Symbol("sg", real=True)
    kj_sg = sp.solve(sp.Eq(israel_S([kj] * 3, [0] * 3, kap)[0], pressure_of_tension(sg_)), kj)[0]
    sg_b = sp.solve(sp.Eq(Kth_m + kj_sg, 0), sg_)[0]
    sg_a = sp.solve(sp.Eq(Kx_m + kj_sg, 0), sg_)[0]
    sg_v = sp.solve(sp.Eq(Km_tr + 3 * kj_sg, 0), sg_)[0]
    ok_x3 = (sp.simplify(sg_b - 2 * sp.sqrt(F) / (kap * r0)) == 0 and sp.simplify(sg_a + (F + 1) / (kap * r0 * sp.sqrt(F))) == 0
             and sp.simplify(sg_v - (3 * F - 1) / (3 * kap * r0 * sp.sqrt(F))) == 0 and sg_b.is_positive and sg_a.is_negative)
    check("X3", ok_x3, "the same thresholds in the tension sigma = -p (action (1)): bdot/b > 0 iff sigma > 2 sqrt(F)/(kappa r0) > 0;"
          " adot/a > 0 iff sigma > -(F+1)/(kappa r0 sqrt F) < 0; volume iff sigma > (3F-1)/(3 kappa r0 sqrt F)")

    # ---------------------------------------------------------------------------- energy conditions
    al, be, v1, v2, v3, p_x, p_t = sp.symbols("alpha beta v1 v2 v3 p_x p_t", real=True)
    e = [sp.sin(al) * sp.cos(be), sp.sin(al) * sp.sin(be), sp.cos(al)]          # unit direction in the orthonormal spatial frame
    Tdiag = [p_x, p_t, p_t]
    NECg = sp.simplify(sum(Tdiag[i] * e[i] ** 2 for i in range(3)))              # T_mu_nu k^mu k^nu, k = (1, e)
    gam2 = 1 / (1 - (v1 ** 2 + v2 ** 2 + v3 ** 2))
    WECg = gam2 * sum(Tdiag[i] * [v1, v2, v3][i] ** 2 for i in range(3))         # T_mu_nu u^mu u^nu, u = gamma (1, v)
    sgm = sp.symbols("sigma_s", real=True)
    sub_iso = {p_x: pressure_of_tension(sgm), p_t: pressure_of_tension(sgm)}
    nec_iso = sp.simplify(NECg.subs(sub_iso))
    wec_iso = sp.simplify(WECg.subs(sub_iso))
    stat = {}
    for sv in (1, 0, -1):
        nv = nec_iso.subs(sgm, sv)
        wv = [sp.simplify(wec_iso.subs({sgm: sv, v1: a, v2: 0, v3: 0})) for a in (0, sp.Rational(1, 2), sp.Rational(9, 10))]
        stat[sv] = ("NEC " + ("violated" if nv < 0 else "saturated" if nv == 0 else "satisfied"),
                    "WEC " + ("violated" if any(w < 0 for w in wv) else "saturated" if all(w == 0 for w in wv) else "satisfied")
                    + (", value 0 for the comoving observer" if wv[0] == 0 else ""))
    ok_e1 = (nec_iso == -sgm and stat[1][0] == "NEC violated" and stat[1][1].startswith("WEC violated")
             and stat[0][0] == "NEC saturated" and stat[-1][0] == "NEC satisfied" and stat[-1][1].startswith("WEC satisfied"))
    for sv in (1, 0, -1):
        print(f"      sigma = {sv:+d}: {stat[sv][0]}, {stat[sv][1]}")
    check("E1", ok_e1, "S-brane, S^i_j = -sigma delta^i_j: T_mu_nu k^mu k^nu = -sigma for every null k, T_mu_nu u^mu u^nu ="
          " -sigma gamma^2 |v|^2; sigma > 0 violates NEC and WEC, sigma = 0 saturates, sigma < 0 satisfies both")
    # general shell: NEC <=> p_x >= 0 and p_t >= 0; [K^th] = -kappa p_x/2
    minNEC_identity = sp.simplify(NECg - (p_t + (p_x - p_t) * sp.sin(al) ** 2 * sp.cos(be) ** 2)) == 0
    kxP, kthP, kxM, kthM = sp.symbols("kxP kthP kxM kthM", real=True)
    SS = israel_S([kxP, kthP, kthP], [kxM, kthM, kthM], kap)
    rel1 = sp.simplify((kthP - kthM) + kap * SS[0] / 2)               # by hand: [K^th] = -kappa p_x / 2
    rel2 = sp.simplify((kxP - kxM) - (kap * SS[1] * (-1) - (kthP - kthM)))   # by hand: [K^x] = -kappa p_t - [K^th]
    pnn = sp.symbols("pnn", nonnegative=True)
    fut_b = Kth_m - kap * pnn / 2                                     # bdot/b(+) for a shell with p_x = pnn >= 0
    check("E2", minNEC_identity and rel1 == 0 and rel2 == 0 and fut_b.is_negative,
          "any shell: T_mu_nu k^mu k^nu = p_t + (p_x - p_t) e_x^2, so NEC <=> p_x, p_t >= 0 (and WEC <=> the same, since S^n_n = 0);"
          " [K^th] = -kappa p_x/2, so a NEC-satisfying shell cannot raise bdot/b: bdot/b(+) <= -sqrt(F)/r0 < 0")

    # ---------------------------------------------------------------------------- de Sitter side
    cc = sp.sqrt(1 + tt ** 2)                        # c = H r0 = cosh(H tau0)
    Hh = cc / r0
    Kx_p = Hh * cc / tt
    Kth_p = Hh * tt / cc
    U = cc + tt                                      # exp(H tau0) = cosh + sinh
    ex = sp.exp(H * tau)
    KxD = sp.simplify((sp.diff(Aexp, tau) / Aexp).rewrite(sp.exp).subs(ex, U).subs(H, Hh))
    KthD = sp.simplify((sp.diff(Bexp, tau) / Bexp).rewrite(sp.exp).subs(ex, U).subs(H, Hh))
    b0 = sp.simplify(Bexp.rewrite(sp.exp).subs(ex, U).subs(H, Hh))
    check("D2", sp.simplify(KxD - Kx_p) == 0 and sp.simplify(KthD - Kth_p) == 0 and sp.simplify(b0 - r0) == 0,
          "de Sitter wedge tau0 > 0, b(tau0) = r0 <=> cosh(H tau0) = H r0 =: c >= 1, t := sinh(H tau0) = sqrt(c^2-1):"
          " K^x_x = H c/t = c^2/(r0 t), K^th_th = H t/c = t/r0 (from the closed form, e^{H tau0} = c + t; S^2 radius r0 by construction)")
    # de Sitter in static form for r > 1/H: f2 = 1 - H^2 r^2 < 0, r increases to the future, n^r = +sqrt(-f2); same code path as the
    # Schwarzschild side
    f2 = 1 - H ** 2 * r_ ** 2
    gS2 = sp.diag(-f2, 1 / f2, r_ ** 2, r_ ** 2 * sp.sin(th) ** 2)
    G2 = einstein_mixed(gS2, XS)
    n2 = -sp.sqrt(-f2) if MUT == "ds_orientation" else sp.sqrt(-f2)
    K2 = K_lie(gS2, XS, [0, n2, 0, 0], [0, 2, 3]).subs(r_, r0).subs(H, Hh).applyfunc(sp.simplify)
    check("D2b", (G2 + 3 * H ** 2 * sp.eye(4)).applyfunc(sp.simplify) == sp.zeros(4, 4)
          and sp.simplify(K2[0, 0] - Kx_p) == 0 and sp.simplify(K2[1, 1] - Kth_p) == 0,
          "de Sitter as f2 = 1 - H^2 r^2 (r > 1/H, r increasing to the future): G^mu_nu = -3H^2 delta^mu_nu and K(+) from (1/2) L_n g with"
          " n^r = +sqrt(-f2) equals D2: the sphere radius r decreases in the Schwarzschild past and increases in the de Sitter future")
    Sd = israel_S([Kx_p, Kth_p, Kth_p], [Kx_m, Kth_m, Kth_m], kap)
    ref_kth = (tt + sp.sqrt(F)) / r0
    ref_kx = cc ** 2 / (r0 * tt) - (F + 1) / (2 * r0 * sp.sqrt(F))
    px_ref = -2 * ref_kth / kap
    pt_ref = -(ref_kx + ref_kth) / kap
    ok_d3 = sp.simplify(Sd[0] - px_ref) == 0 and sp.simplify(Sd[1] - pt_ref) == 0 and px_ref.is_negative
    check("D3", ok_d3, "Schwarzschild(F, r0) -> de Sitter(H r0 = sqrt(1+t^2)) shell: kappa r0 p_x = -2(t + sqrt F) < 0 for every F, t > 0;"
          " kappa r0 p_th = -[(2t^2+1)/t + (F-1)/(2 sqrt F)]")
    # p_th sign: (2t^2+1)/t - 2 sqrt2 = (sqrt2 t - 1)^2/t >= 0; y^2 + 4 sqrt2 y - 1 = (y - (3-2 sqrt2))(y + 3 + 2 sqrt2), y = sqrt F
    y = sp.symbols("y", positive=True)
    id1 = sp.simplify((2 * tt ** 2 + 1) / tt - 2 * sp.sqrt(2) - (sp.sqrt(2) * tt - 1) ** 2 / tt) == 0
    id2 = sp.simplify(sp.expand(y ** 2 + 4 * sp.sqrt(2) * y - 1 - (y - (3 - 2 * sp.sqrt(2))) * (y + 3 + 2 * sp.sqrt(2)))) == 0
    Fc_val = float(17 - 12 * sp.sqrt(2))
    ok_d3b = id1 and id2 and sp.simplify(sp.expand((3 - 2 * sp.sqrt(2)) ** 2) - (17 - 12 * sp.sqrt(2))) == 0
    check("D3b", ok_d3b, f"p_th >= 0 is possible only for F < (3 - 2 sqrt2)^2 = 17 - 12 sqrt2 = {Fc_val:.6f} (r0 > 0.9714 r_s); for"
          " every larger F both pressures are negative at every t > 0: NEC and WEC are violated in all directions")
    # no shell-free matching
    prod_m = sp.simplify(Kx_m * Kth_m)
    prod_p = sp.simplify(Kx_p * Kth_p)
    check("D4", sp.simplify(prod_m + (F + 1) / (2 * r0 ** 2)) == 0 and sp.simplify(prod_p - Hh ** 2) == 0 and prod_m.is_negative
          and prod_p.is_positive,
          "no shell-free matching: K^x K^th = -(F+1)/(2 r0^2) < 0 on the Schwarzschild side, +H^2 > 0 on the de Sitter side;"
          " [K^i_j] = 0 would need equal products")
    # isotropic S-brane condition: K^x(+) - K^th(+) = 1/(r0 t); K^x(-) - K^th(-) = (3F+1)/(2 r0 sqrt F)
    d_p = sp.simplify(Kx_p - Kth_p)
    d_m = sp.simplify(Kx_m - Kth_m)
    t_iso = 2 * sp.sqrt(F) / (3 * F + 1)
    s_iso = sp.simplify((Kth_p + 0 * Kx_p).subs(tt, t_iso) - Kth_m)
    ref_s_iso = 3 * sp.sqrt(F) * (F + 1) / (r0 * (3 * F + 1))
    ok_d5 = (sp.simplify(d_p - 1 / (r0 * tt)) == 0 and sp.simplify(d_m - (3 * F + 1) / (2 * r0 * sp.sqrt(F))) == 0
             and sp.simplify(((Kx_p - Kx_m) - (Kth_p - Kth_m)).subs(tt, t_iso)) == 0
             and sp.simplify(s_iso - ref_s_iso) == 0 and sp.simplify(((Kx_p - Kx_m)).subs(tt, t_iso) - ref_s_iso) == 0
             and ref_s_iso.is_positive
             and sp.simplify(sp.diff(4 * F / (1 + 3 * F) ** 2, F) - 4 * (1 - 3 * F) / (1 + 3 * F) ** 3) == 0
             and sp.simplify((4 * F / (1 + 3 * F) ** 2).subs(F, sp.Rational(1, 3)) - sp.Rational(1, 3)) == 0)
    check("D5", ok_d5, "isotropic S-brane (equal jumps of K^x and K^th) exists iff t = 2 sqrt F/(3F+1), i.e. H r0 = sqrt((F+1)(9F+1))/(3F+1);"
          " its jump is s = 3 sqrt F (F+1)/(r0 (3F+1)) > 0, a positive tension; H r0 < 2/sqrt3, the maximum at F = 1/3 (r0 = 3 r_s/4)")
    # Hamiltonian constraint across the shell
    rho_p = rho_from_constraint([Kx_p, Kth_p, Kth_p], sphere_R3(r0), kap)
    ok6a = sp.simplify(rho_p - 3 * Hh ** 2 / kap) == 0
    Kbar = [(Kx_p + Kx_m) / 2, (Kth_p + Kth_m) / 2, (Kth_p + Kth_m) / 2]
    rho_jump_shell = -(Sd[0] * Kbar[0] + Sd[1] * Kbar[1] + Sd[2] * Kbar[2])
    ok6b = sp.simplify(rho_jump_shell - (rho_p - rho_m)) == 0
    check("D6", ok6a and ok6b,
          "Hamiltonian constraint: kappa rho(+) = 3H^2 on the de Sitter side (R3 = 2/r0^2 continuous), and the shell obeys"
          " [rho] = -S^i_j Kbar^j_i (Gauss-Codazzi; the vacuum energy of de Sitter against the shell stress)")
    # post-brane energy density of the isotropic S-brane
    K_m = Kx_m + 2 * Kth_m
    Kp_iso = [Kx_m + s_, Kth_m + s_, Kth_m + s_]
    rho_plus = rho_from_constraint(Kp_iso, sphere_R3(r0), kap)
    Sbr = israel_S(Kp_iso, [Kx_m, Kth_m, Kth_m], kap)
    Kb_iso = [(a + b) / 2 for a, b in zip(Kp_iso, [Kx_m, Kth_m, Kth_m])]
    rho_shell = -sum(S * k for S, k in zip(Sbr, Kb_iso))
    ref_rho = s_ * (2 * K_m + 3 * s_) / kap
    ok6c = sp.simplify(rho_plus - rho_shell) == 0 and sp.simplify(rho_plus - ref_rho) == 0 and \
        sp.simplify(Sbr[0] + 2 * s_ / kap) == 0 and sp.simplify(Sbr[0] - Sbr[1]) == 0
    check("D6b", ok6c, "isotropic S-brane of jump s: p = -2s/kappa (tension sigma = 2s/kappa), and the energy density of the matter"
          " it leaves, kappa rho_+ = s (2K(-) + 3s) [K(-) = K^x + 2K^th], from the constraint on the + side and from [rho] = -S.Kbar")
    print("\n# thresholds of the jump s = kappa sigma/2 (units 1/r0), Schwarzschild past at r0, F = r_s/r0 - 1: expansion needs s above these")
    print(f"{'F':>8s} {'r0/r_s':>8s} {'r0 s_a':>10s} {'r0 s_b':>10s} {'r0 s_v':>10s}   {'volume contracts at s = 0?':s}")
    for Fv in (1 / 9, 1 / 3, 1.0, 9.0, 99.0):
        sa_, sb_, sv_ = (float(e.subs({F: Fv, r0: 1})) for e in (ref_sa, ref_sb, ref_sv))
        print(f"{Fv:8.4f} {1 / (Fv + 1):8.4f} {sa_:10.4f} {sb_:10.4f} {sv_:10.4f}   {'yes' if sv_ > 0 else 'no'}")
    print("\n# shell stress of Schwarzschild(F) -> de Sitter(c = H r0), units kappa r0 p, c = sqrt(1+t^2)")
    print(f"{'F':>8s} {'r0/r_s':>8s} {'t':>6s} {'H r0':>8s} {'k r0 p_x':>10s} {'k r0 p_th':>10s}   NEC/WEC")
    for Fv in (0.01, float(Fc_val), 0.1, 1 / 3, 1.0, 3.0, 10.0, 100.0):
        for tv in (0.1, 0.7071, 2.0):
            px = float((Sd[0] * kap * r0).subs({F: Fv, tt: tv, r0: 1, kap: 1}))
            pt = float((Sd[1] * kap * r0).subs({F: Fv, tt: tv, r0: 1, kap: 1}))
            nec = "violated" if min(px, pt) < 0 else "ok"
            print(f"{Fv:8.4f} {1 / (Fv + 1):8.4f} {tv:6.4f} {(1 + tv * tv) ** 0.5:8.4f} {px:10.4f} {pt:10.4f}   {nec}")
    print(f"\n# isotropic S-brane into de Sitter: t = 2 sqrt F/(3F+1); r0 s = 3 sqrt F (F+1)/(3F+1)")
    for Fv in (0.01, 0.1, 1.0, 10.0, 100.0, 1e4):
        tv = 2 * Fv ** 0.5 / (3 * Fv + 1)
        print(f"F = {Fv:9.2f}  r0/r_s = {1 / (Fv + 1):.6f}  H r0 = {(1 + tv * tv) ** 0.5:.8f}  r0 s = {3 * Fv ** 0.5 * (Fv + 1) / (3 * Fv + 1):.6f}"
              f"  ; thresholds r0 s_b = {Fv ** 0.5:.6f}")

    failed = [c for c, ok, _ in CHECKS if not ok]
    print(f"\n# {len(CHECKS)} checks, {len(failed)} failed{': ' + ' '.join(failed) if failed else ''}; {time.time() - T0:.1f} s")
    if MUT:
        print(f"MUTATION {MUT}: " + (f"FIRED {' '.join(failed)}" if failed else "SILENT"))
    return 1 if failed else 0


def run_all():
    print("# --mutate=all: each mutation in a subprocess; targets must fire, any other check that fires is listed")
    bad = 0
    for name, (what, targets) in MUTATIONS.items():
        out = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True)
        line = [l for l in out.stdout.splitlines() if l.startswith("MUTATION")]
        fired = set(line[0].split("FIRED")[1].split()) if line and "FIRED" in line[0] else set()
        miss = [t for t in targets if t not in fired]
        extra = sorted(fired - set(targets))
        status = "ON TARGET" if fired and not miss else ("SILENT" if not fired else "MISSED " + " ".join(miss))
        print(f"{name:13s} {status:12s} fired: {' '.join(sorted(fired)) or '-'}   extra: {' '.join(extra) or '-'}   [{what}]")
        if status != "ON TARGET":
            bad += 1
            if out.returncode not in (0, 1):
                print(out.stderr[-600:])
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

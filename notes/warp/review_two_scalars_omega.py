#!/usr/bin/env python3
"""Independent review of #65, part W, for issue #79: the omega construction of the closure stage.

Setting. Two covectors k = a + b and m = a - b (a and b the gradients, up to a field-space boost, of the normal and the
phantom scalar). A boost multiplies them by e^{chi} and e^{-chi}, so closure asks for a function chi with e^{chi} k and
e^{-chi} m closed. With omega = d chi that is
    omega ^ k = -dk,   omega ^ m = dm,   d omega = 0.
#65 states that when k ^ m != 0 the vectors K, M with k(K) = m(M) = 1 and k(M) = m(K) = 0 give the only candidate
    omega = dm(K, M) k + i_K dk.
Proof (by hand, in review_two_scalars.md): contract the first equation with K and the second with M.

This script tests the formula on planted data with a known chi, on planted data that cannot be closed, on the
Ellis-Bronnikov metric with an added normal scalar (my own coordinate Ricci), and on the two shifts b(rho) with D = 0.

    python3 notes/warp/review_two_scalars_omega.py                  checks, then every control
    python3 notes/warp/review_two_scalars_omega.py --no-controls    checks only
    python3 notes/warp/review_two_scalars_omega.py --mutate=NAME    one control (its target checks must fail, exit 1)
"""
import random
import subprocess
import sys
import time

import mpmath as mp
import sympy as sp

OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate", "")
CONTROLS = {
    "omega-sign": ("omega = -dm(K, M) k + i_K dk", ["omega_planted", "omega_eb"]),
    "omega-swap": ("omega = dm(K, M) m + i_M dm", ["omega_planted"]),
    "omega-nodk": ("omega = dm(K, M) k", ["omega_planted"]),
    "plant-open": ("planted k replaced by a non-closed form, recovery expected", ["omega_planted"]),
    "eb-lambda": ("lambda = 2 kappa^2 replaced by kappa^2", ["omega_eb"]),
    "brho-eps": ("split with theta^0 - eps theta^1", ["omega_brho"]),
}
if MUT and MUT not in CONTROLS:
    sys.exit(f"unknown mutation {MUT}")

RESULTS = []
mp.mp.dps = 80


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# Forms on four coordinates. A 1-form is a list of four components; a 2-form F = (1/2) F_ij dx^i ^ dx^j is the
# antisymmetric matrix F_ij. (d a)_ij = d_i a_j - d_j a_i, (a ^ b)_ij = a_i b_j - a_j b_i, (i_K F)_j = K^i F_ij.
# ---------------------------------------------------------------------------------------------------------------
def d1(a, X):
    return sp.Matrix(4, 4, lambda i, j: sp.diff(a[j], X[i]) - sp.diff(a[i], X[j]))


def wedge11(a, b):
    return sp.Matrix(4, 4, lambda i, j: a[i] * b[j] - a[j] * b[i])


def interior(K, F):
    return sp.Matrix([sum(K[i] * F[i, j] for i in range(4)) for j in range(4)])


def omega_candidate(k, m, dk, dm, K, M):
    """#65's candidate omega = dm(K, M) k + i_K dk, with the controls' mutations."""
    dmKM = sum(K[i] * M[j] * dm[i, j] for i in range(4) for j in range(4))
    kv, mv = sp.Matrix(k), sp.Matrix(m)
    if MUT == "omega-sign":
        return -dmKM * kv + interior(K, dk)
    if MUT == "omega-swap":
        return dmKM * mv + interior(M, dm)
    if MUT == "omega-nodk":
        return dmKM * kv
    return dmKM * kv + interior(K, dk)


def dual_vectors(k, m, extra=None):
    """K, M with k(K) = m(M) = 1, k(M) = m(K) = 0 (minimum-norm), optionally plus a kernel vector of (k, m)."""
    N = sp.Matrix([list(k), list(m)])
    G = (N * N.T).inv()
    sol = N.T * G
    K, M = sol[:, 0], sol[:, 1]
    if extra is not None:
        ker = N.nullspace()
        K = K + extra[0] * ker[0] + extra[1] * ker[1]
        M = M + extra[2] * ker[0] + extra[3] * ker[1]
    return K, M


def residuals(om, k, m, dk, dm):
    """omega ^ k + dk and omega ^ m - dm, as 4 x 4 matrices."""
    return wedge11(om, k) + dk, wedge11(om, m) - dm


def zero_mat(M):
    return all(sp.simplify(e) == 0 for e in M)


def main():
    t0 = time.time()
    print("review_two_scalars_omega.py, issue #79. sympy %s, mpmath %s, python %s. mutation: %s"
          % (sp.__version__, mp.__version__, sys.version.split()[0], MUT or "none"))
    rng = random.Random(7917)
    X = sp.symbols("t x y z", real=True)
    t, x, y, z = X

    def randpoly(deg=2, c=3):
        mons = [1] + list(X) + [X[i] * X[j] for i in range(4) for j in range(i, 4)]
        mons = mons if deg >= 2 else mons[:5]
        return sum(rng.randint(-c, c) * mn for mn in mons)

    # ---- 1. Planted data with a known chi ------------------------------------------------------------------------
    print("\n== 1. Planted pairs: k = e^{-chi} (d phi1 + d phi2), m = e^{chi} (d phi1 - d phi2), omega = d chi ==")
    nplant = 8
    good = 0
    notes = []
    dclosed_checked = 0
    for n in range(nplant):
        pt = {t: sp.Rational(rng.randint(-5, 5), rng.randint(1, 3)), x: sp.Rational(rng.randint(-5, 5), rng.randint(1, 3)),
              y: sp.Rational(rng.randint(-5, 5), rng.randint(1, 3)), z: sp.Rational(rng.randint(-5, 5), rng.randint(1, 3))}
        phi1, phi2 = randpoly(), randpoly()
        qpoly = sum(rng.randint(-2, 2) * (X[i] - pt[X[i]]) * (X[j] - pt[X[j]]) for i in range(4) for j in range(i, 4)) \
            + sum(rng.randint(-2, 2) * (X[i] - pt[X[i]]) for i in range(4))
        chi = sp.log(2) + qpoly  # chi(p) = ln 2, so that e^{+-chi(p)} = 2^{+-1} is rational at the point
        d1p = [sp.diff(phi1, v) for v in X]
        d2p = [sp.diff(phi2, v) for v in X]
        k = [sp.exp(-chi) * (d1p[i] + d2p[i]) for i in range(4)]
        m = [sp.exp(chi) * (d1p[i] - d2p[i]) for i in range(4)]
        if MUT == "plant-open":  # a non-closed k: add y dz
            k[3] = k[3] + y * sp.exp(-chi)
        dk, dm = d1(k, X), d1(m, X)
        dchi = [sp.diff(chi, v).subs(pt) for v in X]
        kp = [sp.simplify(e.subs(pt)) for e in k]
        mp_ = [sp.simplify(e.subs(pt)) for e in m]
        dkp, dmp = dk.subs(pt).applyfunc(sp.simplify), dm.subs(pt).applyfunc(sp.simplify)
        if sp.Matrix([kp, mp_]).rank() < 2:
            continue
        omegas = []
        for choice in range(2):
            ex = None if choice == 0 else [sp.Rational(rng.randint(-3, 3), 2) for _ in range(4)]
            K, M = dual_vectors(kp, mp_, ex)
            assert sp.simplify(sum(K[i] * kp[i] for i in range(4)) - 1) == 0 and sp.simplify(sum(M[i] * mp_[i] for i in range(4)) - 1) == 0
            omegas.append(omega_candidate(kp, mp_, dkp, dmp, K, M).applyfunc(sp.simplify))
        same = all(sp.simplify(om[i] - dchi[i]) == 0 for om in omegas for i in range(4))
        r1, r2 = residuals(omegas[0], sp.Matrix(kp), sp.Matrix(mp_), dkp, dmp)
        eqs = zero_mat(r1) and zero_mat(r2)
        ok_n = same and eqs
        # d omega = 0 for the symbolic field, on the first two samples
        if n < 2 and ok_n and MUT == "":
            Kf, Mf = dual_vectors(k, m)
            omf = omega_candidate(k, m, dk, dm, Kf, Mf)
            dom = d1(list(omf), X).subs(pt).applyfunc(sp.simplify)
            dclosed_checked += 1
            ok_n = ok_n and zero_mat(dom)
        good += ok_n
        if not ok_n:
            notes.append(n)
    check("omega_planted", good == nplant, "%d of %d planted pairs: the candidate equals d chi at the point for two different "
          "(K, M), solves omega ^ k = -dk and omega ^ m = dm, and d omega = 0 for the symbolic field on %d of them; failing "
          "samples %s" % (good, nplant, dclosed_checked, notes))

    # ---- 2. Planted pairs that cannot be closed -------------------------------------------------------------------
    print("\n== 2. Pairs with k ^ dk != 0: no chi exists, the candidate must violate omega ^ k = -dk ==")
    nbad, nfire = 4, 0
    for n in range(nbad):
        while True:
            kk = [randpoly(1, 3) for _ in range(4)]
            mm = [randpoly(1, 3) for _ in range(4)]
            pt = {t: sp.Rational(1, 2), x: sp.Rational(-1, 3), y: sp.Rational(2, 5), z: sp.Rational(3, 7)}
            dkk = d1(kk, X).subs(pt)
            kp = [e.subs(pt) for e in kk]
            mp_ = [e.subs(pt) for e in mm]
            frob = any(kp[i] * dkk[j, l] + kp[j] * dkk[l, i] + kp[l] * dkk[i, j] != 0
                       for i, j, l in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)])
            if sp.Matrix([kp, mp_]).rank() == 2 and frob:
                break
        dmm = d1(mm, X).subs(pt)
        K, M = dual_vectors(kp, mp_)
        om = omega_candidate(kp, mp_, dkk, dmm, K, M)
        r1, r2 = residuals(om, sp.Matrix(kp), sp.Matrix(mp_), dkk, dmm)
        nfire += not (zero_mat(r1) and zero_mat(r2))
    check("omega_planted_fail", nfire == nbad, "%d of %d planted non-closable pairs: the candidate violates omega ^ k = -dk or "
          "omega ^ m = dm" % (nfire, nbad))

    # ---- 3. Ellis-Bronnikov with an added normal scalar ------------------------------------------------------------
    print("\n== 3. Ellis-Bronnikov with a normal scalar: my coordinate Ricci, lambda = 2 kappa^2, omega = 0 and omega = -d chi ==")
    T_, L_, TH, PH = sp.symbols("T L TH PH", real=True)
    XE = (T_, L_, TH, PH)
    kap, beta = sp.symbols("kappa beta", positive=True)
    R2 = (1 - beta * sp.cos(2 * kap * L_)) / (2 * kap**2)
    gE = sp.diag(-1, 1, R2, R2 * sp.sin(TH) ** 2)
    pt = {T_: 0, TH: sp.pi / 2, PH: 0}
    g0, Ric = coord_ricci_symbolic(gE, XE, pt)
    lam_e = (kap**2 if MUT == "eb-lambda" else 2 * kap**2)
    bl2 = (1 - beta**2) / (2 * kap**2 * R2**2)  # b_l^2 = 8 pi phi2'^2 with phi2' = C / r^2, 4 pi C^2 = (1 - beta^2)/(4 kappa^2)
    a_cov = sp.Matrix([sp.sqrt(2) * kap, 0, 0, 0])
    b_sq = sp.diag(0, bl2, 0, 0)
    resid = (Ric - lam_e * g0 - a_cov * a_cov.T + b_sq).applyfunc(lambda e: reduce_trig(e, kap, L_))
    ok_ricci = all(e == 0 for e in resid)
    # phi2' from the printed arctan against C / r^2, 80 digits at random parameters
    ok_phi = True
    for _ in range(6):
        kv, bv, lv = mp.mpf(rng.randint(1, 9)) / 7, mp.mpf(rng.randint(1, 9)) / 11, mp.mpf(rng.randint(-9, 9)) / 13
        phi2p = mp.diff(lambda ll: mp.atan(mp.sqrt((1 + bv) / (1 - bv)) * mp.tan(kv * ll)) / (2 * mp.sqrt(mp.pi)), lv)
        r2v = (1 - bv * mp.cos(2 * kv * lv)) / (2 * kv**2)
        want_sq = (1 - bv**2) / (16 * mp.pi * kv**2 * r2v**2)
        ok_phi &= abs(phi2p**2 - want_sq) < mp.mpf(10) ** -40 * max(1, abs(want_sq))
    # G - 8 pi T with phi1 = kappa t / (2 sqrt(pi)), phi2' = C / r^2, V = kappa^2 / (4 pi)
    Rscal = sum(g0.inv()[i, i] * Ric[i, i] for i in range(4))
    Gt = Ric - Rscal * g0 / 2
    dphi1 = sp.Matrix([kap / (2 * sp.sqrt(sp.pi)), 0, 0, 0])
    dphi2 = sp.Matrix([0, sp.sqrt(bl2 / (8 * sp.pi)), 0, 0])
    Vv = kap**2 / (4 * sp.pi)
    ginv = g0.inv()
    n1 = (dphi1.T * ginv * dphi1)[0]
    n2 = (dphi2.T * ginv * dphi2)[0]
    Tt = dphi1 * dphi1.T - dphi2 * dphi2.T - g0 * (n1 / 2 - n2 / 2 + Vv)
    res_G = (Gt - 8 * sp.pi * Tt).applyfunc(lambda e: reduce_trig(e, kap, L_))
    ok_G = all(e == 0 for e in res_G)
    # omega: k = a + b, m = a - b with a = sqrt2 kappa dt, b = b_l dl, b_l = bl2^(1/2) (positive branch)
    bl = sp.sqrt(1 - beta**2) / (sp.sqrt(2) * kap * R2)
    kE = [sp.sqrt(2) * kap, bl, 0, 0]
    mE = [sp.sqrt(2) * kap, -bl, 0, 0]
    dkE, dmE = d1(kE, XE), d1(mE, XE)
    KE, ME = dual_vectors(kE, mE)
    omE = omega_candidate(kE, mE, dkE, dmE, KE, ME).applyfunc(sp.simplify)
    ok_om0 = all(e == 0 for e in omE)
    # a planted boost chi(t, l): k -> e^{chi} k, m -> e^{-chi} m; then omega = -d chi
    chi = T_ * L_ + L_**2 / 3
    kB = [sp.exp(chi) * e for e in kE]
    mB = [sp.exp(-chi) * e for e in mE]
    dkB, dmB = d1(kB, XE), d1(mB, XE)
    KB, MB = dual_vectors(kB, mB)
    omB = omega_candidate(kB, mB, dkB, dmB, KB, MB)
    want = sp.Matrix([-sp.diff(chi, v) for v in XE])
    par = {kap: sp.Rational(3, 5), beta: sp.Rational(1, 3)}
    ok_omB = True
    for tv, lv in ((sp.Rational(1, 3), sp.Rational(2, 7)), (-sp.Rational(1, 2), sp.Rational(5, 9)), (sp.Rational(3, 4), -sp.Rational(1, 5))):
        pv = {**par, T_: tv, L_: lv, TH: sp.pi / 2, PH: 0}
        diff = [sp.N((omB[i] - want[i]).subs(pv), 50) for i in range(4)]
        ok_omB &= all(abs(e) < sp.Float(10) ** -35 for e in diff)
    check("omega_eb", ok_ricci and ok_phi and ok_G and ok_om0 and ok_omB,
          "R_ab = lambda g_ab + a a - b b with lambda = %s (symbolic in l, kappa, beta; identically zero residual: %s); phi2' = C/r^2 "
          "against the arctan formula at 6 random points (80 digits): %s; G - 8 pi T = 0: %s; omega = 0 for the closed split: %s; "
          "omega = -d chi after the boost chi = t l + l^2/3 (50 digits at 3 points, kappa = 3/5, beta = 1/3): %s"
          % (lam_e, ok_ricci, ok_phi, ok_G, ok_om0, ok_omB))

    # ---- 4. The two shifts b(rho) with D = 0: no omega ---------------------------------------------------------------
    print("\n== 4. b(rho) with b'' + b'/rho = eps b'^2: the split, the candidate omega and the Frobenius obstruction ==")
    Tb, Xb, Rb, Fb = sp.symbols("Tb Xb Rb Fb", positive=True)
    XB = (Tb, Xb, Rb, Fb)
    bf = sp.Function("bf")(Rb)
    uu, bb_ = sp.symbols("uu bb_")
    ok_all = True
    notes = []
    for eps in (1, -1):
        gB = sp.Matrix([[-1 + bf**2, bf, 0, 0], [bf, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, Rb**2]])
        ptb = {Tb: 0, Xb: 0, Fb: 0}
        g0, RicB = coord_ricci_symbolic(gB, XB, ptb)

        def ode(e, eps=eps):
            """replace b'' by eps b'^2 - b'/rho, then b' and b by the symbols uu and bb_"""
            e = sp.sympify(e).subs(sp.Derivative(bf, (Rb, 2)), eps * sp.Derivative(bf, Rb) ** 2 - sp.Derivative(bf, Rb) / Rb)
            return e.subs(sp.Derivative(bf, Rb), uu).subs(bf, bb_)

        sgn = -eps if MUT == "brho-eps" else eps
        h = uu**2 / 2
        nu0 = sp.Matrix([1 + sgn * bb_, sgn, 0, 0])
        rho_cov = sp.Matrix([0, 0, 1, 0])
        want = h * (rho_cov * rho_cov.T - nu0 * nu0.T)
        split_ok = (RicB.applyfunc(ode) - want).applyfunc(sp.simplify) == sp.zeros(4, 4)
        # k = sqrt(h) (d rho + nu0), m = sqrt(h) (d rho - nu0) with sqrt(h) = sigma u / sqrt 2, as functions of rho via b
        viol = []
        for sigma in (1, -1):
            sh = sigma * sp.Derivative(bf, Rb) / sp.sqrt(2)
            nuf = sp.Matrix([1 + eps * bf, eps, 0, 0])
            kf = [sh * nuf[0], sh * nuf[1], sh, 0]
            mf = [-sh * nuf[0], -sh * nuf[1], sh, 0]
            dkf, dmf = d1(kf, XB), d1(mf, XB)
            Kf, Mf = dual_vectors(kf, mf)
            omf = omega_candidate(kf, mf, dkf, dmf, Kf, Mf)
            r1, r2 = residuals(omf, sp.Matrix(kf), sp.Matrix(mf), dkf, dmf)
            nz = [sp.simplify(ode(e)) for e in list(r1) + list(r2)]
            viol.append(any(e != 0 for e in nz))
        nu = [(1 + eps * bf), eps, 1, 0]
        dnu = d1(nu, XB)
        comp = ode(nu[0] * dnu[1, 2] + nu[1] * dnu[2, 0] + nu[2] * dnu[0, 1])
        frob_ok = sp.simplify(comp - uu) == 0
        ok_all &= split_ok and all(viol) and frob_ok
        notes.append((eps, split_ok, viol, frob_ok))
    check("omega_brho", ok_all, "for b'' + b'/rho = eps b'^2, eps = +-1, exact in (b, b', rho): R_ab = h ((d rho)^2 - ((1 + eps b) dt + "
          "eps dx)^2); for both signs of sqrt(h) the candidate omega violates omega ^ k = -dk or omega ^ m = dm; "
          "(nu ^ d nu)_{t x rho} = b' (eps, split, violation, Frobenius): %s" % notes)

    fails = [n for n, ok in RESULTS if not ok]
    print("\nchecks: %d passed, %d failed %s. %.0f s" % (len(RESULTS) - len(fails), len(fails), fails, time.time() - t0))
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


def coord_ricci_symbolic(g, X, pt):
    """My coordinate-Christoffel Ricci (same routine as review_two_scalars_cartan.coord_ricci_at, without its control),
    with the coordinates left symbolic where pt does not fix them. No symbolic matrix inverse beyond g itself."""
    n = 4
    ev = lambda e: sp.sympify(e).subs(pt)
    g0 = sp.Matrix(n, n, lambda a, b: ev(g[a, b]))
    dg = [[[ev(sp.diff(g[a, b], X[d])) for b in range(n)] for a in range(n)] for d in range(n)]
    ddg = [[[[ev(sp.diff(g[a, b], X[d], X[e])) for b in range(n)] for a in range(n)] for e in range(n)] for d in range(n)]
    gi = g0.inv()
    low = lambda a, b, c: (dg[b][a][c] + dg[c][a][b] - dg[a][b][c]) / 2
    dlow = lambda d, a, b, c: (ddg[d][b][a][c] + ddg[d][c][a][b] - ddg[d][a][b][c]) / 2
    Gam = [[[sum(gi[a, e] * low(e, b, c) for e in range(n)) for c in range(n)] for b in range(n)] for a in range(n)]
    dgi = [[[-sum(gi[a, f] * dg[d][f][h] * gi[h, e] for f in range(n) for h in range(n)) for e in range(n)]
            for a in range(n)] for d in range(n)]
    dGam = [[[[sum(dgi[d][a][e] * low(e, b, c) + gi[a, e] * dlow(d, e, b, c) for e in range(n)) for c in range(n)]
              for b in range(n)] for a in range(n)] for d in range(n)]
    Ric = sp.zeros(n, n)
    for b in range(n):
        for c in range(n):
            Ric[b, c] = (sum(dGam[a][a][b][c] - dGam[c][a][b][a] for a in range(n))
                         + sum(Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][b][a] for a in range(n) for d in range(n)))
    return g0, Ric


def reduce_trig(e, kap, L_):
    """exact zero test for an expression in cos(2 kappa l), sin(2 kappa l) and rational functions of the parameters."""
    C, S = sp.symbols("C S")
    e = sp.sympify(e).subs({sp.cos(2 * kap * L_): C, sp.sin(2 * kap * L_): S})
    e = sp.together(e)
    num, den = sp.fraction(e)
    num = sp.expand(num)
    if num.has(S):
        pn = sp.Poly(num, S)
        out = 0
        for (k_,), co in pn.terms():
            out += co * (1 - C**2) ** (k_ // 2) * S ** (k_ % 2)
        num = sp.expand(out)
    return sp.simplify(num)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Independent review of #65, part L, for issue #79: the lemma on random tensors.

Lemma (#65). Let R be real symmetric on an orthonormal frame (e_0, ..., e_3) with R_a3 = 0 for every a and
R_00 + R_11 + 2 R_22 = 0, and let R_3 be its block on (e_0, e_1, e_2). Then R = lambda eta + a a - b b with lambda real and
a, b real covectors iff lambda = 0 and det R_3 = 0.

The proof is by hand in review_two_scalars.md. Here the statement is tested without the proof's route. The search for
lambda takes the real roots of the quartic det(R - lambda eta) (exact real-root isolation), not the gcd of the 3 x 3
minors, and at each root counts the signs of the eigenvalues of the symmetric matrix R - lambda eta at 80 digits.
A form a a - b b has at most one positive and one negative eigenvalue, and by Sylvester's law that count is what is tested.

    python3 notes/warp/review_two_scalars_lemma.py                  checks, then every control
    python3 notes/warp/review_two_scalars_lemma.py --no-controls    checks only
    python3 notes/warp/review_two_scalars_lemma.py --mutate=NAME    one control (its target checks must fail, exit 1)
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
    "lemma-weight": ("hypothesis R_00 - R_11 + 2 R_22 = 0", ["lemma_if", "weighted_trace"]),
    "lemma-r33": ("R_33 allowed nonzero", ["lemma_if"]),
    "planted-sign": ("planted lambda_0 eta + a a + b b expected to pass", ["search_planted"]),
    "trace-sign": ("trace reversal with the wrong sign", ["ricci_from_T"]),
    "search-eig": ("inertia from the eigenvalues of eta^{-1} M", ["search_planted"]),
}
if MUT and MUT not in CONTROLS:
    sys.exit(f"unknown mutation {MUT}")

RESULTS = []
mp.mp.dps = 80
TOL = mp.mpf(10) ** -40
ETA = sp.diag(-1, 1, 1, 1)
W = (1, -1, 2) if MUT == "lemma-weight" else (1, 1, 2)  # weights of the hypothesis W_0 R_00 + W_1 R_11 + W_2 R_22 = 0


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# Random rational data
# ---------------------------------------------------------------------------------------------------------------
def rr(rng, lo=-9, hi=9):
    return sp.Rational(rng.randint(lo, hi), rng.randint(1, 4))


def rand_lorentz(rng):
    """A rational Lorentz transformation: a boost along a random axis times a rational rotation (Cayley transform)."""
    m, n = rng.choice([(2, 1), (3, 1), (3, 2), (4, 1), (5, 2), (5, 3)])
    gam, gb = sp.Rational(m * m + n * n, m * m - n * n), sp.Rational(2 * m * n, m * m - n * n)
    ax = rng.randint(1, 3)
    B = sp.eye(4)
    B[0, 0] = B[ax, ax] = gam
    B[0, ax] = B[ax, 0] = gb
    S = sp.zeros(3, 3)
    for i in range(3):
        for j in range(i + 1, 3):
            S[i, j] = rr(rng, -3, 3) / 3
            S[j, i] = -S[i, j]
    Rot3 = (sp.eye(3) - S) * (sp.eye(3) + S).inv()
    Rot = sp.diag(1, Rot3)
    L = B * Rot
    assert L.T * ETA * L == ETA
    return L


def rand_R3_ndet(rng):
    """random R_3 with the weighted trace zero (weights W) and det != 0."""
    while True:
        r00, r11, r01, r02, r12 = (rr(rng) for _ in range(5))
        r22 = -(W[0] * r00 + W[1] * r11) / sp.Integer(W[2])
        R3 = sp.Matrix([[r00, r01, r02], [r01, r11, r12], [r02, r12, r22]])
        if R3.det() != 0:
            return R3


def rand_R3_det0(rng, force_zero=False):
    """random R_3 with the weighted trace zero (weights W) and det = 0: u u^T - s w w^T, s fixed by the trace."""
    while True:
        u = sp.Matrix([rng.randint(-4, 4) for _ in range(3)])
        w = u if force_zero else sp.Matrix([rng.randint(-4, 4) for _ in range(3)])
        nu = sum(W[i] * u[i] ** 2 for i in range(3))
        nw = sum(W[i] * w[i] ** 2 for i in range(3))
        if nw == 0 or nu == 0 and not force_zero:
            continue
        s = sp.Rational(nu, nw)
        R3 = u * u.T - s * w * w.T
        if force_zero and R3 != sp.zeros(3, 3):
            continue
        if not force_zero and R3 == sp.zeros(3, 3):
            continue
        return R3


def embed(R3, r33=0):
    R = sp.zeros(4, 4)
    R[:3, :3] = R3
    R[3, 3] = r33
    return R


def lorentz(R, L):
    return L.T * R * L


# ---------------------------------------------------------------------------------------------------------------
# The lambda search: real roots of det(R - lambda eta), then the inertia of R - lambda eta at each root.
# ---------------------------------------------------------------------------------------------------------------
lam = sp.Symbol("lambda")


def inertia(R, lam_value):
    """(n_plus, n_minus, n_zero) of the real symmetric matrix R - lam_value * eta at 80 digits (R: sympy, rational)."""
    lv = mp.mpf(sp.N(lam_value, 90))
    A = mp.matrix(4, 4)
    for i in range(4):
        for j in range(4):
            A[i, j] = mp.mpf(sp.N(R[i, j], 90)) - lv * mp.mpf(int(ETA[i, j]))
    if MUT == "search-eig":
        A = mp.matrix([[(-1 if i == 0 else 1) * A[i, j] for j in range(4)] for i in range(4)])  # eta^{-1} M
        ev = mp.eig(A, left=False, right=False)
        ev = [e.real if abs(e.imag) < TOL else mp.mpf(10) for e in ev]  # a complex eigenvalue counts as positive
    else:
        ev = mp.eigsy(A, eigvals_only=True)
    npos = sum(1 for e in ev if e > TOL)
    nneg = sum(1 for e in ev if e < -TOL)
    return npos, nneg, 4 - npos - nneg


def lambda_set(R):
    """the real lambda with R - lambda eta of inertia (<= 1, <= 1, .), as a list of sympy numbers."""
    p = sp.Poly((R - lam * ETA).det(), lam)
    assert not p.is_zero
    out = []
    for root in sorted(set(sp.real_roots(p)), key=lambda r_: float(r_)):
        npos, nneg, _ = inertia(R, root)
        if npos <= 1 and nneg <= 1:
            out.append(root)
    return out


def split_ab(R3):
    """a, b real with R = a a^T - b b^T for a rank <= 2 R_3 of inertia (<= 1, <= 1), at 80 digits; returns residual."""
    A = mp.matrix(3, 3)
    for i in range(3):
        for j in range(3):
            A[i, j] = mp.mpf(sp.N(R3[i, j], 90))
    ev, Q = mp.eigsy(A)
    a = mp.matrix(3, 1)
    b = mp.matrix(3, 1)
    for k in range(3):
        v = Q[:, k]
        if ev[k] > TOL:
            a = mp.sqrt(ev[k]) * v
        elif ev[k] < -TOL:
            b = mp.sqrt(-ev[k]) * v
    res = A - (a * a.T - b * b.T)
    return max(abs(res[i, j]) for i in range(3) for j in range(3))


def main():
    t0 = time.time()
    print("review_two_scalars_lemma.py, issue #79. sympy %s, mpmath %s, python %s. mutation: %s"
          % (sp.__version__, mp.__version__, sys.version.split()[0], MUT or "none"))
    rng = random.Random(7915)

    # ---- 1. The stress tensor and the trace reversal ------------------------------------------------------------
    print("\n== 1. R_ab from G_ab = 8 pi T_ab for a normal and a phantom scalar ==")
    pv = sp.Matrix(sp.symbols("p0:4"))
    qv = sp.Matrix(sp.symbols("q0:4"))
    Vp = sp.Symbol("V")
    etai = ETA.inv()
    X1 = (pv.T * etai * pv)[0]
    X2 = (qv.T * etai * qv)[0]
    T = pv * pv.T - qv * qv.T - ETA * (X1 / 2 - X2 / 2 + Vp)
    Gm = 8 * sp.pi * T
    trG = sum(etai[a, a] * Gm[a, a] for a in range(4))
    # R_ab - R eta_ab / 2 = G_ab and tr(G) = -R in four dimensions
    sgn = 1 if MUT == "trace-sign" else -1
    Rab = Gm + sgn * sp.Rational(1, 2) * trG * ETA
    want = 8 * sp.pi * (pv * pv.T - qv * qv.T) + 8 * sp.pi * Vp * ETA
    d1 = (Rab - want).applyfunc(sp.expand)
    check("ricci_from_T", d1 == sp.zeros(4, 4), "R_ab = 8 pi (d phi1 d phi1 - d phi2 d phi2) + 8 pi V eta_ab: "
          "difference %s" % (d1 if d1 != sp.zeros(4, 4) else 0))

    # ---- 2. The weighted trace of N --------------------------------------------------------------------------------
    print("\n== 2. N = R_3 - lambda diag(-1, 1, 1): weighted trace under the hypothesis ==")
    r = sp.symbols("r0:3", real=True)
    L_ = sp.Symbol("lam", real=True)
    N0, N1, N2 = r[0] + L_, r[1] - L_, r[2] - L_
    hyp = W[0] * r[0] + W[1] * r[1] + W[2] * r[2]
    wt = sp.expand(W[0] * N0 + W[1] * N1 + W[2] * N2 - (hyp - 2 * L_))
    # the lemma's step: with hyp = 0, N_00 + N_11 + 2 N_22 = -2 lambda (real weights 1, 1, 2)
    wt_hyp = sp.expand(W[0] * N0 + W[1] * N1 + W[2] * N2 + 2 * L_ - hyp)
    check("weighted_trace", wt_hyp == 0, "weights %s: W.N + 2 lambda - (hypothesis) = %s" % (W, wt_hyp))

    # ---- 3. Planted truth for the search --------------------------------------------------------------------------
    print("\n== 3. Planted truth for the lambda search ==")
    ok_pass, ok_fail, notes = 0, 0, []
    npass, nfail = 30, 30
    for i in range(npass):
        lam0 = sp.Integer(0) if i % 3 == 0 else rr(rng, -9, 9)
        if lam0 == 0 and i % 3 != 0:
            lam0 = sp.Integer(1)
        a = sp.Matrix([rr(rng) for _ in range(4)])
        b = sp.Matrix([rr(rng) for _ in range(4)])
        sb = 1 if MUT == "planted-sign" else -1  # control: aa + bb expected to pass
        R = lam0 * ETA + a * a.T + sb * b * b.T
        L = rand_lorentz(rng)
        S = lambda_set(lorentz(R, L))
        ok_pass += any(abs(sp.N(s - lam0, 60)) < 1e-30 for s in S)
    for i in range(nfail):
        lam0 = rr(rng, -9, 9)
        a, b, c = (sp.Matrix([rr(rng) for _ in range(4)]) for _ in range(3))
        kind = i % 3
        if kind == 0:
            R = lam0 * ETA + a * a.T + b * b.T
        elif kind == 1:
            R = lam0 * ETA - a * a.T - b * b.T
        else:
            R = lam0 * ETA + a * a.T - b * b.T + c * c.T
        S = lambda_set(lorentz(R, rand_lorentz(rng)))
        ok_fail += len(S) == 0
        if S:
            notes.append((kind, [sp.N(s, 8) for s in S]))
    ok_gen = 0
    for i in range(20):
        Rg = sp.zeros(4, 4)
        for a_ in range(4):
            for b_ in range(a_, 4):
                Rg[a_, b_] = Rg[b_, a_] = rr(rng)
        ok_gen += len(lambda_set(Rg)) == 0
    check("search_planted", ok_pass == npass and ok_fail == nfail and ok_gen == 20,
          "planted lambda_0 eta + a a - b b recovered %d/%d (lambda_0 = 0 in a third), planted (+,+), (-,-), (+,-,+) "
          "give the empty set %d/%d, generic symmetric 4x4 give the empty set %d/20 %s"
          % (ok_pass, npass, ok_fail, nfail, ok_gen, notes))

    # ---- 4. lambda != 0 directly -------------------------------------------------------------------------------------
    print("\n== 4. The lambda != 0 branch, tested at random lambda (not only at roots) ==")
    bad = 0
    total = 0
    for _ in range(100):
        R = embed(rand_R3_ndet(rng) if rng.random() < 0.5 else rand_R3_det0(rng), 0)
        for lam0 in (rr(rng, 1, 9), -rr(rng, 1, 9), sp.Rational(1, 7), -sp.Rational(1, 7)):
            if lam0 == 0:
                continue
            npos, nneg, _ = inertia(R, lam0)
            total += 1
            bad += npos <= 1 and nneg <= 1
    check("lambda_branch", bad == 0,
          "%d random (R, lambda != 0) pairs satisfying the hypotheses: %d have inertia within (1, 1)" % (total, bad))

    # ---- 5. If: det R_3 = 0 gives lambda = 0 and a, b ---------------------------------------------------------------
    print("\n== 5. det R_3 = 0 with the hypotheses: the lambda-set is {0} and a, b rebuild R ==")
    nsamp, good, worst = 40, 0, mp.mpf(0)
    n_rank2 = 0
    for i in range(nsamp):
        force0 = i < 2
        R3 = rand_R3_det0(rng, force_zero=force0)
        assert R3.det() == 0
        r33 = rr(rng, 1, 9) if MUT == "lemma-r33" else 0
        R = embed(R3, r33)
        S = lambda_set(lorentz(R, rand_lorentz(rng)))
        okset = len(S) == 1 and S[0] == 0
        res = split_ab(R3) if okset else mp.mpf(1)
        worst = max(worst, res)
        n_rank2 += R3 != sp.zeros(3, 3)
        good += okset and res < TOL
    check("lemma_if", good == nsamp, "%d of %d tensors with det R_3 = 0 (R_3 != 0 in %d): lambda-set = {0}, a and b rebuild R, "
          "worst residual %s" % (good, nsamp, n_rank2, mp.nstr(worst, 3)))

    # ---- 6. Only if: det R_3 != 0 gives the empty set ---------------------------------------------------------------
    print("\n== 6. det R_3 != 0 with the hypotheses: the lambda-set is empty ==")
    nsamp, good = 40, 0
    for i in range(nsamp):
        R3 = rand_R3_ndet(rng)
        R = embed(R3, 0)
        S = lambda_set(lorentz(R, rand_lorentz(rng)))
        good += len(S) == 0
    check("lemma_onlyif", good == nsamp, "%d of %d tensors with det R_3 != 0 have the empty lambda-set" % (good, nsamp))

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


if __name__ == "__main__":
    main()

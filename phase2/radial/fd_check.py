"""
Independent discretization: global second-order finite differences (box scheme) for the same
first-order system y' = J(r; sigma^2) y = (J0 + sigma^2 J1) y on [r0, R], giving a sparse generalized
eigenproblem K0 y = -sigma^2 K1 y, solved by shift-invert Arnoldi around a target.

BCs: at r0 the two regularity relations (A' = f r0/3 with f linear in sigma^2,
chi = r0^3 (Gt - 2 w A) a / (3 alpha)); at R Dirichlet A(R) = Gt(R) = 0 (the growing pair is
suppressed; error ~ e^{-2 kappa (R - r_core)}, negligible for bound modes not too near threshold).
Uniform-in-x grid with x = log-ish stretching r = r0 + (R - r0) * s^2-type map concentrated at small r.
This shares only the generated RHS (the derived equations) with radial_modes.py, not the method.
"""
import numpy as np, sys, warnings
import scipy.sparse as sps
from scipy.sparse.linalg import eigs
warnings.filterwarnings('ignore')
from radial_modes import RadialProblem
import radial_rhs_generated as G
from background import LAM, MU

def fd_eigs(P, n, return_vecs=False):
    r0, R = P.r0, P.R
    s = np.linspace(0, 1, n)
    r = r0 + (R - r0) * (0.6 * s**2 + 0.4 * s)       # finer near the core
    rm = 0.5 * (r[1:] + r[:-1]); h = np.diff(r)
    J0 = np.array([P.J(x, 0.0) for x in rm]); J1 = np.array([P.J(x, 1.0) for x in rm]) - J0
    N4 = 4 * n
    rows0, cols0, v0 = [], [], []
    rows1, cols1, v1 = [], [], []
    I = np.eye(4)
    for i in range(n - 1):
        B0p = I / h[i] - 0.5 * J0[i]; B0m = -I / h[i] - 0.5 * J0[i]
        B1 = -0.5 * J1[i]
        for a in range(4):
            for b in range(4):
                rr = 4 * i + a
                rows0 += [rr, rr]; cols0 += [4 * (i + 1) + b, 4 * i + b]; v0 += [B0p[a, b], B0m[a, b]]
                if B1[a, b] != 0:
                    rows1 += [rr, rr]; cols1 += [4 * (i + 1) + b, 4 * i + b]; v1 += [B1[a, b], B1[a, b]]
    # BCs at r0: rows 4(n-1) .. 4(n-1)+1 ; at R: rows +2, +3
    al, a, phi0, u = P.st.bg(r0)
    base = 4 * (n - 1)
    # chi0 - r0^3 a/(3 al) (Gt - 2 w A) = 0
    c = r0**3 * a / (3 * al)
    rows0 += [base, base, base]; cols0 += [3, 2, 0]; v0 += [1.0, -c, 2 * P.w * c]
    # A' - (r0/3) f(A, Gt, chi; s2) = 0, f linear: f = fA A + fG Gt + fC chi + s2 (gA A)
    def f(A0, G0, C0, s2): return G.ddA(r0, al, a, phi0, u, A0, 0.0, G0, C0, P.w, P.eps, s2, LAM, MU)
    fA, fG, fC = f(1, 0, 0, 0), f(0, 1, 0, 0), f(0, 0, 1, 0)
    gA = f(1, 0, 0, 1) - fA
    rows0 += [base + 1] * 4; cols0 += [1, 0, 2, 3]; v0 += [1.0, -r0 / 3 * fA, -r0 / 3 * fG, -r0 / 3 * fC]
    rows1 += [base + 1]; cols1 += [0]; v1 += [-r0 / 3 * gA]
    rows0 += [base + 2, base + 3]; cols0 += [4 * (n - 1) + 0, 4 * (n - 1) + 2]; v0 += [1.0, 1.0]
    K0 = sps.csc_matrix((v0, (rows0, cols0)), shape=(N4, N4))
    K1 = sps.csc_matrix((v1, (rows1, cols1)), shape=(N4, N4))
    # (K0 + s2 K1) y = 0  <=>  K0 y = s2 (-K1) y.  K0 is very ill-conditioned (cond ~ 1e22, from the
    # exponentially small phi0 coefficients), which defeats sparse shift-invert; use dense QZ.
    import scipy.linalg as sl
    vals, vecs = sl.eig(K0.toarray(), -K1.toarray())
    ok = np.isfinite(vals) & (np.abs(vals.imag) < 1e-6) & (np.abs(vals) < 1.0)
    vals, vecs = vals[ok].real, vecs[:, ok].real
    o = np.argsort(vals)
    if return_vecs:
        return vals[o], vecs[:, o], r
    return vals[o]

NS = (250, 400, 550)   # dense QZ cost ~n^3; 1200 timed out on the shared 2-core machine
if __name__ == "__main__":
    from background import solve_star
    cases = [(+1, [0.92, 0.94]), (-1, [0.7603, 0.9275, 0.94])]
    for eps, ws in cases:
        st = solve_star(0.75, eps); cur = 0.75
        for w in ws:
            while abs(w - cur) > 1e-12:
                cur = round(cur + np.clip(w - cur, -0.02, 0.02), 6); st = solve_star(cur, eps, guess=st)
            P = RadialProblem(st)
            out = []
            for n in NS:
                out.append(fd_eigs(P, n))
            print(f"eps={eps:+d} w={w}: lowest FD eigenvalues sigma^2 (|.|<1e-8 dropped: trivial mode, see RESULTS)"
                  f"  continuum threshold (1-w)^2={(1-w)**2:.5f}")
            for n, v in zip(NS, out):
                v = v[np.abs(v) > 1e-8]
                print(f"     n={n}: {np.array2string(v[:3], precision=7)}")
            sys.stdout.flush()

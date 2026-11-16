#!/usr/bin/env python3
"""An independent recomputation of #75's sector spectra of S^-1 M for the tanh profile (issue #94).

Sector J_z = m of sigma.grad psi = v (tau/2) psi with tau = -f'(r) cos(theta), axis z. Components
psi = (A e^{i(m-1/2)phi}, B e^{i(m+1/2)phi}) with A, B real functions of (r, theta), so
    S_m(A, B) = (d_z A + d_rho B + (m + 1/2) B/rho,  d_rho A - (m - 1/2) A/rho - d_z B).
Tensor cubic B-splines in (r, theta): uniform in r on [0, rcore], geometric beyond it to Rout, uniform in theta. Dirichlet
at r = Rout. B vanishes on the axis, and A too for m > 1/2. Galerkin in <S a, S b>: N c = mu K c, mu = 1/v, solved as
the dense eigenproblem of C = L^-1 N L^-T with K = L L^T. The planted tau = 6/(1 + r^2) has mu = +-i/l_n,
l_n = 1 + 2n/3, with multiplicity 2(n + 3/2 - m) in sector m for n >= m - 1/2.

Shares with tanh_spectrum.py (#75): the formulation (Galerkin compression of S^-1 M in H^1-dot) and the profile. Not the
basis (B-splines in (r, theta) against spin-weighted angular functions and mapped radial functions), the coordinates, the
reduction to a sector or the solver. Shares nothing with frank_loss.py beyond the profile. It grew from the review script
in #75's closing comment.

    python notes/warp/tanh_spectrum_check.py --part=planted [--ms=...]
    python notes/warp/tanh_spectrum_check.py --part=solve --sigma=2 --disc=A --out=DIR [--ms=...]
    python notes/warp/tanh_spectrum_check.py --part=report --out=DIR
    python notes/warp/tanh_spectrum_check.py --part=planted --mutate=kterm      (a control, must exit 1)

`solve` writes DIR/s<sigma>_<disc>_m<2m>.npz per sector, one sector at a time, and skips sectors already written.
`report` reads them all, prints the tables and exits 1 if a check fails. Run one process per (sigma, disc) with
OMP_NUM_THREADS=1. Discretisation A is h = 0.1 min(2/sigma, 1/(1 + m/10)) on [0, 4], then q = 1.15 to Rout = 800, and
24 (1 + (m - 1/2)/20) theta intervals. B changes the spline counts and Rout: h x 0.9, q = 1.13, Rout = 1200, and
28 (1 + (m - 1/2)/20) theta intervals.
"""
import argparse
import glob
import os
import sys
import time

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
from scipy.interpolate import BSpline

ap = argparse.ArgumentParser()
ap.add_argument("--part", required=True, choices=["planted", "solve", "report"])
ap.add_argument("--mutate", default="")
ap.add_argument("--sigma", type=float, default=2.0)
ap.add_argument("--disc", default="A", choices=["A", "B"])
ap.add_argument("--ms", default="0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5")
ap.add_argument("--out", default="")
ap.add_argument("--nth", type=int, default=0)       # overrides of the discretisation, for convergence studies
ap.add_argument("--rout", type=float, default=0.0)
ap.add_argument("--hfac", type=float, default=0.0)
ap.add_argument("--q", type=float, default=0.0)
ap.add_argument("--rcore", type=float, default=0.0)
args = ap.parse_args()
MUT = args.mutate
if MUT not in ("", "kterm", "half"):
    sys.exit(f"unknown mutation {MUT}")
K4, NG = 3, 6  # cubic splines, Gauss points per interval
S_SOB = 3 * (np.pi / 2) ** (4 / 3)
MS = [float(x) for x in args.ms.split(",")]
# #75's leading |mu| v_c per sector (tanh_spectrum.md section 4) and its "converged above" values
REF = {2.0: {"vc": 3.744900, "lead": [0.409556, 0.298994, 0.232982, 0.190371, 0.160798, 0.139127, 0.122581, 0.109541,
                                      0.099002, 0.090309, 0.083018],
             "conv": [0.0777, 0.0803, 0.0751, 0.0711, 0.0667, 0.0688, 0.0649, 0.0611, 0.0578, 0.0595, 0.0563]},
       4.0: {"vc": 2.470971, "lead": [0.367657, 0.288531, 0.235378, 0.198480, 0.171512, 0.150969, 0.134809, 0.121766,
                                      0.111020, 0.102014, 0.094358],
             "conv": [0.0969, 0.0910, 0.0858, 0.0806, 0.0840, 0.0790, 0.0745, 0.0705, 0.0727, 0.0688, 0.0653]}}
DISC = {"A": dict(hfac=1.0, q=1.15, rout=800.0, nth=24, rcore=4.0), "B": dict(hfac=0.9, q=1.13, rout=1200.0, nth=28, rcore=4.0)}
for _d in DISC.values():
    _d.update({k: v for k, v in (("nth", args.nth), ("rout", args.rout), ("hfac", args.hfac), ("q", args.q), ("rcore", args.rcore)) if v})
FAILED = []


def check(name, ok, msg):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)
    if not ok:
        FAILED.append(name)


def knots_r(h, rcore, q, rout):
    br = list(np.arange(0.0, rcore + 1e-12, h))
    d = h
    while br[-1] < rout:
        d *= q
        br.append(min(br[-1] + d, rout) if rout - br[-1] - d > 0.3 * d else rout)
    return np.array(br)


def basis1d(br):
    t = np.r_[[br[0]] * K4, br, [br[-1]] * K4]
    n = len(t) - K4 - 1
    g, w = np.polynomial.legendre.leggauss(NG)
    a, b = br[:-1, None], br[1:, None]
    x = (0.5 * (b - a) * g + 0.5 * (a + b)).ravel()
    wx = (0.5 * (b - a) * w).ravel()
    V = np.zeros((len(x), n))
    dV = np.zeros((len(x), n))
    for j in range(n):
        c = np.zeros(n)
        c[j] = 1.0
        s = BSpline(t, c, K4)
        V[:, j] = s(x)
        dV[:, j] = s.derivative()(x)
    return x, wx, V, dV


def sector_C(tau_fn, m, br, nth):
    """C = L^-1 N L^-T for sector m. Its eigenvalues are the Galerkin mu."""
    r, wr, Nr, dNr = basis1d(br)
    th, wt, Tt, dTt = basis1d(np.linspace(0.0, np.pi, nth + 1))
    R, TH = np.meshgrid(r, th, indexing="ij")
    R, TH = R.ravel(), TH.ravel()
    W = np.outer(wr * r**2, wt * np.sin(th)).ravel()
    c, s = np.cos(TH), np.sin(TH)
    rho = R * s
    kr, kt = Nr.shape[1], Tt.shape[1]
    keep_r = np.arange(kr - 1)  # Dirichlet at Rout
    axis0 = np.arange(1, kt - 1)  # zero on the axis
    kA = np.arange(kt) if m == 0.5 else axis0
    blocks = {}
    for name, kt_ in (("A", kA), ("B", axis0)):
        Nr_, dNr_ = sp.csr_matrix(Nr[:, keep_r]), sp.csr_matrix(dNr[:, keep_r])
        T_, dT_ = sp.csr_matrix(Tt[:, kt_]), sp.csr_matrix(dTt[:, kt_])
        F = sp.kron(Nr_, T_).tocsr()
        Fr = sp.kron(dNr_, T_).tocsr()
        Ft = sp.kron(Nr_, dT_).tocsr()
        blocks[name] = (F, sp.diags(c) @ Fr - sp.diags(s / R) @ Ft, sp.diags(s) @ Fr + sp.diags(c / R) @ Ft)
    FA, dzA, drA = blocks["A"]
    FB, dzB, drB = blocks["B"]
    kterm = 0.0 if MUT == "kterm" else (m + 0.5)
    S1 = sp.hstack([dzA, drB + sp.diags(kterm / rho) @ FB]).tocsr()
    S2 = sp.hstack([drA - sp.diags((m - 0.5) / rho) @ FA, -dzB]).tocsr()
    Mv = 0.5 * tau_fn(R, c) * (2.0 if MUT == "half" else 1.0)
    F1 = sp.hstack([FA, sp.csr_matrix(FB.shape)]).tocsr()
    F2 = sp.hstack([sp.csr_matrix(FA.shape), FB]).tocsr()
    Wd, WM = sp.diags(W), sp.diags(W * Mv)
    K = (S1.T @ Wd @ S1 + S2.T @ Wd @ S2).toarray()
    N = (S1.T @ WM @ F1 + S2.T @ WM @ F2).toarray()
    L = sla.cholesky(K, lower=True, overwrite_a=True, check_finite=False)
    del K
    C = sla.solve_triangular(L, N.T, lower=True, check_finite=False)
    del N
    C = sla.solve_triangular(L, C.T, lower=True, overwrite_b=True, check_finite=False)
    return C


def sech2(x):
    """sech^2 without overflow at large |x|."""
    e = np.exp(-2 * np.abs(x))
    return 4 * e / (1 + e) ** 2


def tanh_fp(sig):
    return lambda r: sig * (sech2(sig * (r + 1)) - sech2(sig * (r - 1))) / (2 * np.tanh(sig))


def v_c(sig):
    g, w = np.polynomial.legendre.leggauss(400)
    rr, ww = 15 * (g + 1) / 2, 15 * w / 2
    return 2 * np.sqrt(S_SOB) / (np.pi * np.sum(ww * rr**2 * np.abs(tanh_fp(sig)(rr))**3)) ** (1 / 3)


def sizes(sig, disc, m):
    """The grid of sector m. High sectors need both a finer core and more theta intervals (the pilot: at m = 21/2,
    refining either alone leaves the leading eigenvalue 2e-5 to 5e-5 low against #75)."""
    d = DISC[disc]
    h = 0.1 * d["hfac"] * min(2.0 / sig, 1.0 / (1.0 + m / 10.0))
    nth = int(round(d["nth"] * (1.0 + (m - 0.5) / 20.0)))
    return knots_r(h, d["rcore"], d["q"], d["rout"]), nth, h


def part_planted():
    print(f"== planted tau = 6/(1 + r^2), discretisation A at sigma R = 2 (Rout = {DISC['A']['rout']:g}), "
          f"mutation {MUT or 'none'}")
    for m in MS:
        br, nth, h = sizes(2.0, "A", m)
        t0 = time.time()
        C = sector_C(lambda R, c: 6.0 / (1.0 + R**2), m, br, nth)
        mu = sla.eigvals(C, overwrite_a=True, check_finite=False)
        mu = mu[np.argsort(-np.abs(mu))]
        exact = []
        n0 = int(round(m - 0.5))
        for k in range(n0, n0 + 6):
            exact += [1 / (1 + 2 * k / 3)] * (2 * int(round(k + 1.5 - m)))
        exact = np.array(exact)
        got = np.abs(mu[:len(exact)])
        err = np.max(np.abs(got - exact))
        per = [float(np.max(np.abs(got[exact == lv] - lv))) for lv in sorted(set(exact), reverse=True)]
        re = np.max(np.abs(mu[:len(exact)].real))
        check(f"planted m={m:g}", err < 1e-5 and re < 1e-8,
              f"h = {h:.4f}, nth = {nth}, n = {C.shape[0]}, levels l_{n0} to l_{n0 + 5} with multiplicities ({len(exact)} eigenvalues), worst "
              f"|mu| error {err:.1e} (per level {', '.join(f'{e:.0e}' for e in per)}), max |Re mu| {re:.1e}, {time.time() - t0:.0f} s")


def part_solve():
    sig = args.sigma
    vc = v_c(sig)
    os.makedirs(args.out, exist_ok=True)
    print(f"== solve sigma R = {sig:g}, disc {args.disc} (q = {DISC[args.disc]['q']}, Rout = {DISC[args.disc]['rout']:g}), "
          f"v_c = {vc:.6f}", flush=True)
    for m in MS:
        br, nth, h = sizes(sig, args.disc, m)
        fn = os.path.join(args.out, f"s{sig:g}_{args.disc}_m{int(2 * m)}.npz")
        if os.path.exists(fn):
            print(f"  m = {m:g}: {fn} exists, skipped", flush=True)
            continue
        t0 = time.time()
        C = sector_C(lambda R, c: -tanh_fp(sig)(R) * c, m, br, nth)
        n = C.shape[0]
        mu = sla.eigvals(C, overwrite_a=True, check_finite=False)
        del C
        mu = mu[np.argsort(-np.abs(mu))]
        np.savez(fn, mu=mu, vc=vc, n=n, h=h, nth=nth, rout=DISC[args.disc]["rout"], q=DISC[args.disc]["q"])
        print(f"  m = {m:g}: h = {h:.4f}, nth = {nth}, n = {n}, leading mu v_c = {mu[0] * vc:.6f}, {time.time() - t0:.0f} s", flush=True)


def part_report():
    """Each sector's window at two matching tolerances. TOLS[0] = 1e-5 is #75's: an eigenvalue of A is converged if B has
    one within 1e-5 relative, and the window is |mu| v_c above the largest unconverged one. This basis converges fewer
    eigenvalues than #75's spectral one (the pilot: 0.18 against #75's 0.078 at m = 1/2, sigma R = 2), so the window is
    narrower; TOLS[1] = 1e-4 is reported beside it. Real eigenvalues are searched in the wider of the two windows."""
    TOLS, tol_lead, tol_real = (1e-5, 1e-4), 1e-5, 1e-9
    for sig in (2.0, 4.0):
        vc = v_c(sig)
        ref = REF[sig]
        check(f"v_c sR={sig:g}", abs(vc - ref["vc"]) < 1e-5, f"{vc:.6f} against #75's {ref['vc']}")
        print(f"\n== sigma R = {sig:g}, v_c = {vc:.6f}")
        print("  m    | n A / B       | leading mu v_c (A)  | |mu| v_c A | B        | #75      | rel. to #75 |"
              " converged above, 1e-5 / 1e-4 (#75) | eigs in window | real in window A / B")
        bound = {tol: 0.0 for tol in TOLS}
        lead_last = None
        for m in MS:
            fa = os.path.join(args.out, f"s{sig:g}_A_m{int(2 * m)}.npz")
            fb = os.path.join(args.out, f"s{sig:g}_B_m{int(2 * m)}.npz")
            if not (os.path.exists(fa) and os.path.exists(fb)):
                check(f"present sR={sig:g} m={m:g}", False, f"{fa} or {fb} missing")
                continue
            A, B = np.load(fa), np.load(fb)
            ma, mb = A["mu"] * vc, B["mu"] * vc
            am, bm = np.abs(ma), np.abs(mb)
            top_a, top_b = ma[am > 1e-3], mb[bm > 1e-3]
            rel_a = np.array([np.min(np.abs(mb - z)) / abs(z) for z in top_a])
            rel_b = np.array([np.min(np.abs(ma - z)) / abs(z) for z in top_b])

            def conv(tol):  # the largest |mu| v_c of an eigenvalue without a partner, in either direction
                bad = np.r_[np.abs(top_a[rel_a > tol]), np.abs(top_b[rel_b > tol]), 0.0]
                return bad.max()
            cv = {tol: conv(tol) for tol in TOLS}
            wide = min(cv.values())
            win_a, win_b = top_a[np.abs(top_a) > wide], top_b[np.abs(top_b) > wide]
            real_a = win_a[np.abs(win_a.imag) <= tol_real * np.abs(win_a)]
            real_b = win_b[np.abs(win_b.imag) <= tol_real * np.abs(win_b)]
            lead_a, lead_b = am[0], bm[0]
            k = int(round(m - 0.5))
            r75 = ref["lead"][k] if k < len(ref["lead"]) else np.nan
            c75 = ref["conv"][k] if k < len(ref["conv"]) else np.nan
            rel = abs(lead_a - r75) / r75
            for tol in TOLS:
                bound[tol] = max(bound[tol], cv[tol])
            lead_last = lead_a
            print(f"  {m:4g} | {int(A['n']):5d} / {int(B['n']):5d} | {ma[0].real:+.5f}{ma[0].imag:+.5f}i | {lead_a:.6f}   | "
                  f"{lead_b:.6f} | {r75:.6f} | {rel:.1e}     | {cv[TOLS[0]]:.4f} / {cv[TOLS[1]]:.4f} ({c75:.4f})"
                  f"           | {len(win_a):4d}           | {len(real_a)} / {len(real_b)}", flush=True)
            for z in list(real_a) + list(real_b):
                print(f"         real eigenvalue in a window: mu v_c = {z.real:+.8f}")
            check(f"lead sR={sig:g} m={m:g}", rel < tol_lead, f"leading |mu| v_c {lead_a:.6f} against #75's {r75:.6f}")
            check(f"AB sR={sig:g} m={m:g}", abs(lead_a - lead_b) < 1e-6 * lead_a,
                  f"A and B leading differ by {abs(lead_a - lead_b) / lead_a:.1e} relative")
            sym = max(np.min(np.abs(ma + z)) / abs(z) for z in win_a) if len(win_a) else 0.0
            check(f"sym sR={sig:g} m={m:g}", sym < 1e-4, f"window symmetric under mu -> -mu to {sym:.1e} relative")
            check(f"noreal sR={sig:g} m={m:g}", len(real_a) == 0 and len(real_b) == 0,
                  f"no real eigenvalue with |mu| v_c above {wide:.4f}, in A or B")
        if lead_last is not None:
            for tol in TOLS:
                top = max(bound[tol], lead_last)
                print(f"  tolerance {tol:.0e}: every eigenvalue with |mu| v_c >= {top:.4f} is converged in every sector "
                      f"computed (the leading eigenvalue of the last sector is {lead_last:.4f}), so no real zero mode with "
                      f"|v| <= {1 / top:.2f} v_c, provided the leading eigenvalue keeps falling with m "
                      f"(#75, at 1e-5 with its own basis: {'12.05' if sig == 2 else '10.32'} v_c)")


t_start = time.time()
{"planted": part_planted, "solve": part_solve, "report": part_report}[args.part]()
print(f"wall {time.time() - t_start:.0f} s")
print("FAILED: " + ", ".join(FAILED) if FAILED else "all checks passed")
sys.exit(1 if FAILED else 0)

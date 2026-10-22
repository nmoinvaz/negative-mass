"""
Tasks 2-3: truncated DG problem modelling a regular star of areal radius R (|M| = 1):
    H^ Psi = E Psi on x > x_R = x(R),   Robin  Psi_x(x_R) = kappa Psi(x_R)   (or Dirichlet),
quadratic form  Q[Psi] = int_{x_R}^inf (Psi_x^2 + V^ Psi^2) dx + kappa Psi(x_R)^2.

Discretisation: Psi^ = r^3 u and the regular u-equation of common.py, Chebyshev collocation on
r = R + a (1+t)/(1-t), a = max(R, 1), t = 1 (r = inf) node removed (u = 0); the boundary row at
r = R is replaced by the Robin condition f u_r + (3f/R - kappa) u = 0 (zero mass row), or the node
is removed for Dirichlet.

Analytic prediction (exact): phi0 = 1/psi0 = ((l+2)(l-1) r - 6)/chi(r) solves H^ phi0 = 0, is smooth
for r > 0, decays as r^-(l+1) (normalisable), and has a single zero at r_s = 6/((l+2)(l-1)).
By Sturm oscillation theory for a regular left endpoint and limit-point right endpoint:
    #negative eigenvalues(Dirichlet)   = #zeros of phi0 in (R, inf) = [R < r_s]
    #negative eigenvalues(Robin kappa) = [R < r_s] + [kappa < kappa_c(R)],
    kappa_c(R) = d ln(phi0)/dx at r = R  (kappa = kappa_c gives a zero-energy bound state).
The script checks this against the numerical spectra on a grid of (R, kappa).
"""
import os, csv
import numpy as np, scipy.linalg as sl
import common as c

HERE = os.path.dirname(os.path.abspath(__file__))


def trunc_eigs(l, R, kappa, N=160, a=None):
    """Real eigenvalues (sorted) of the truncated problem; kappa = 'D' for Dirichlet.
    Uses Psi^ = r^3 u (same regularised u-equation as the full problem, see common.py); this is
    far better conditioned than the bare Psi^ form when R << 1 (V^ ~ 12/r^4 there)."""
    a = a or max(R, 1.0)
    D, t = c.cheb(N)
    idx = np.arange(1, N + 1)
    tt = t[idx]
    r = R + a * (1 + tt) / (1 - tt)
    r[-1] = R
    rt = 2 * a / (1 - tt)**2
    rtt = 4 * a / (1 - tt)**3
    Dt = D[np.ix_(idx, idx)]
    D2t = (D @ D)[np.ix_(idx, idx)]
    Dr = Dt / rt[:, None]
    Dr2 = D2t / rt[:, None]**2 - (rtt / rt**3)[:, None] * Dt
    n = c.numeric(l)
    A = -np.diag(r) @ Dr2 + np.diag(-6 + 2 / (r + 2)) @ Dr + np.diag(n['C'](r))
    B = np.diag(r**3 / (r + 2)**2)
    if isinstance(kappa, str):
        A, B = A[:-1, :-1], B[:-1, :-1]
    else:
        # Psi_x = f (r^3 u)_r = kappa r^3 u  <=>  f u_r + (3 f / R - kappa) u = 0
        f = 1 + 2 / R
        A[-1, :] = f * Dr[-1, :]
        A[-1, -1] += 3 * f / R - kappa
        B[-1, :] = 0
    w = sl.eigvals(A, B)
    w = w[np.isfinite(w)]
    w = w[np.abs(w.imag) < 1e-7 * (1 + np.abs(w.real))].real
    return np.sort(w)


def n_neg(w, tol=1e-9):
    return int(np.sum(w < -tol))


def kappa_c(l, R):
    return float(c.numeric(l)['kappa_c'](R))


def kappa_c_numeric(l, R, N=160):
    """Bisection on kappa for the change in the number of negative eigenvalues (purely numerical)."""
    base = n_neg(trunc_eigs(l, R, 'D', N))   # count for kappa -> +inf
    lo, hi = -1.0, 1.0
    while n_neg(trunc_eigs(l, R, lo, N)) == base:
        lo *= 2
        if lo < -1e7: return np.nan
    while n_neg(trunc_eigs(l, R, hi, N)) != base:
        hi = hi * 2 if hi > 0 else 1.0
        if hi > 1e7: return np.nan
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if n_neg(trunc_eigs(l, R, mid, N)) == base: hi = mid
        else: lo = mid
        if hi - lo < 1e-9 * (1 + abs(hi)): break
    return 0.5 * (lo + hi)


Rgrid = np.round(np.geomspace(0.05, 10, 25), 4)
Lset = (2, 3, 4, 5)


def main():
    print("r_s(l) = 6/((l+2)(l-1)) [zero of phi0; Dirichlet unstable iff R < r_s]")
    print("r_*(l) = zero of phi0' beyond r_s [kappa_c(R) < 0 iff R > r_*; Neumann unstable iff R < r_*]")
    import sympy as sp
    rstar = {}
    for l in range(2, 9):
        s = c.symbolic(l)
        num, _ = sp.fraction(sp.cancel(s['dlogphi0_dx']))
        rts = sorted(float(sp.re(z)) for z in sp.Poly(num, c.rs_).nroots() if abs(sp.im(z)) < 1e-12 and sp.re(z) > c.r_s(l))
        rstar[l] = rts[0]
    from vhat_rneg import RNEG
    print(f"{'l':>2} {'r_s':>8} {'r_*':>8} {'r_neg':>8}")
    for l in range(2, 9):
        print(f"{l:>2} {c.r_s(l):>8.4f} {rstar[l]:>8.4f} {RNEG[l]:>8.4f}")

    # ---- kappa_c(R): analytic vs numerical bisection
    print("\nkappa_c(R): analytic d ln phi0/dx vs numerical bisection (N=160 and N=240)")
    print(f"{'l':>2} {'R':>8} {'x_R':>10} {'kc_exact':>14} {'kc_num160':>14} {'kc_num240':>14} {'kc*x_R':>9} {'kc*R':>8}")
    kc_rows = []
    for l in Lset:
        for R in Rgrid:
            if abs(R - c.r_s(l)) < 1e-3: continue
            xR = float(c.x_of_r(R))
            ke = kappa_c(l, R)
            kn1 = kappa_c_numeric(l, R, 160) if R in Rgrid[::3] else np.nan
            kn2 = kappa_c_numeric(l, R, 240) if R in Rgrid[::3] else np.nan
            kc_rows.append([l, R, xR, ke, kn1, kn2, ke * xR, ke * R])
            if R in Rgrid[::3]:
                print(f"{l:>2} {R:>8.4f} {xR:>10.4e} {ke:>14.6e} {kn1:>14.6e} {kn2:>14.6e} {ke*xR:>9.3f} {ke*R:>8.3f}")
    # dense analytic table for plotting
    Rd = np.geomspace(0.05, 10, 400)
    with open(os.path.join(HERE, 'kappa_c.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['R'] + [f'kappa_c_l{l}' for l in Lset] + [f'kappa_c_x_R_l{l}' for l in Lset])
        for R in Rd:
            xR = float(c.x_of_r(R))
            kk = [kappa_c(l, R) if abs(R - c.r_s(l)) > 1e-4 else np.nan for l in Lset]
            w.writerow([R] + kk + [k * xR for k in kk])
    with open(os.path.join(HERE, 'kappa_c_check.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['l', 'R', 'x_R', 'kappa_c_exact', 'kappa_c_num_N160', 'kappa_c_num_N240', 'kc_xR', 'kc_R'])
        w.writerows(kc_rows)

    # ---- full grid (R, kappa x_R) vs prediction
    print("\nGrid scan: kappa*x_R in {-10,-3,-1,-0.3,0,0.3,1,3,10} + Dirichlet, R in 25 log-spaced values [0.05,10]")
    scan = []
    mism = 0
    for l in Lset:
        for R in Rgrid:
            if abs(R - c.r_s(l)) < 1e-3: continue
            xR = float(c.x_of_r(R))
            for kx in (-10, -3, -1, -0.3, 0, 0.3, 1, 3, 10, 'D'):
                kap = 'D' if kx == 'D' else kx / xR
                w1 = trunc_eigs(l, R, kap, 160)
                w2 = trunc_eigs(l, R, kap, 240)
                pred = int(R < c.r_s(l)) + (0 if kx == 'D' else int(kap < kappa_c(l, R)))
                n1, n2 = n_neg(w1), n_neg(w2)
                E0 = w2[0] if n2 else np.nan
                if not (n1 == n2 == pred): mism += 1; print('MISMATCH', l, R, kx, n1, n2, pred, w2[:3])
                scan.append([l, R, xR, kx, pred, n1, n2, E0, np.sqrt(-E0) if n2 else 0.0,
                             abs(w1[0] - w2[0]) / abs(w2[0]) if n2 else np.nan])
    print(f"cases: {len(scan)}, mismatches between N=160, N=240 and Sturm prediction: {mism}")
    with open(os.path.join(HERE, 'truncated_scan.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['l', 'R', 'x_R', 'kappa_x_R(or D)', 'n_neg_predicted', 'n_neg_N160', 'n_neg_N240',
                                        'E0_N240', 'growth_sqrt_minusE0', 'rel_diff_E0_N160_vs_N240'])
        w.writerows(scan)
    rel = [q[-1] for q in scan if np.isfinite(q[-1])]
    print(f"max relative difference of lowest negative eigenvalue N=160 vs N=240: {max(rel):.2e}")

    # ---- analytic-bound check: R > r_neg and kappa >= 0 => no negative eigenvalue
    viol = [q for q in scan if q[1] > RNEG[q[0]] and (q[3] == 'D' or q[3] >= 0) and q[6] > 0]
    print(f"R > r_neg(l), kappa >= 0 (incl. Dirichlet): cases {sum(1 for q in scan if q[1] > RNEG[q[0]] and (q[3]=='D' or q[3]>=0))}, with negative eigenvalue: {len(viol)}")

    # ---- small R: Dirichlet / Neumann reproduce DG growth?
    print("\nSmall R: lowest eigenvalue, Dirichlet and Neumann vs DG -(k/2)^2")
    print(f"{'l':>2} {'R':>7} {'E_D':>18} {'E_N':>18} {'E_DG':>12} {'relD':>9} {'relN':>9}")
    small = []
    for l in Lset:
        Eex = -(c.kval(l) / 2)**2
        for R in (1e-3, 3e-3, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0):
            wd = trunc_eigs(l, R, 'D', 240); wn = trunc_eigs(l, R, 0.0, 240)
            ED = wd[0] if n_neg(wd) else np.nan; EN = wn[0] if n_neg(wn) else np.nan
            small.append([l, R, ED, EN, Eex])
            print(f"{l:>2} {R:>7.3f} {ED:>18.10f} {EN:>18.10f} {Eex:>12.1f} {abs(ED/Eex-1):>9.1e} {abs(EN/Eex-1):>9.1e}")
    with open(os.path.join(HERE, 'small_R.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['l', 'R', 'E_Dirichlet', 'E_Neumann', 'E_DG']); w.writerows(small)

    # ---- growth-rate table vs l
    print("\nGrowth rate sqrt(-E0) (units 1/|M|), 0 = stable.  Full DG problem: k/2.")
    Rs = (0.05, 0.3, 1.0, 3.0, 10.0)
    BCs = (('D', 'Dirichlet'), (0.0, 'Neumann'), (-1.0, 'kappa=-1'), (-10.0, 'kappa=-10'))
    hdr = f"{'l':>2} {'k/2':>8} " + " ".join(f"{name[:5]+'@'+str(R):>12}" for R in Rs for _, name in BCs)
    gr_rows = []
    for l in range(2, 9):
        row = [l, c.kval(l) / 2]
        for R in Rs:
            for kap, _ in BCs:
                w = trunc_eigs(l, R, kap, 240)
                row.append(np.sqrt(-w[0]) if n_neg(w) else 0.0)
        gr_rows.append(row)
    with open(os.path.join(HERE, 'growth_rates.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['l', 'full_DG_k_over_2'] + [f'R={R}_{name}' for R in Rs for _, name in BCs])
        w.writerows(gr_rows)
    for R in Rs:
        print(f"\nR = {R}:   {'l':>2} {'DG k/2':>8} " + " ".join(f"{name:>11}" for _, name in BCs))
        for row in gr_rows:
            j = 2 + 4 * Rs.index(R)
            print(f"           {row[0]:>2} {row[1]:>8.1f} " + " ".join(f"{v:>11.4f}" for v in row[j:j + 4]))

    # ---- l -> large at fixed R: last unstable l under Neumann / Dirichlet  (uses exact r_s, r_*)
    print("\nNeumann unstable iff R < r_*(l); Dirichlet unstable iff R < r_s(l) = 6/((l+2)(l-1)).")

    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        for l in Lset:
            kk = np.array([kappa_c(l, R) if abs(R - c.r_s(l)) > 1e-4 else np.nan for R in Rd])
            xr = c.x_of_r(Rd)
            ax[0].plot(Rd, kk * xr, label=f'l={l}')
            ax[0].axvline(c.r_s(l), ls=':', lw=0.7, color=ax[0].lines[-1].get_color())
        for q in kc_rows:
            if np.isfinite(q[5]): ax[0].plot(q[1], q[5] * q[2], 'k.', ms=4)
        ax[0].set_xscale('log'); ax[0].set_ylim(-12, 12); ax[0].axhline(0, color='gray', lw=0.5)
        ax[0].set_xlabel('R/|M|'); ax[0].set_ylabel('kappa_c x_R'); ax[0].legend(fontsize=8)
        ax[0].set_title('critical Robin parameter (dots: numerical bisection)')
        L = np.array([q[0] for q in gr_rows])
        ax[1].semilogy(L, [q[1] for q in gr_rows], 'k-o', label='full DG (k/2)')
        for R in (0.05, 0.3, 1.0):
            j = 2 + 4 * Rs.index(R)
            for off, mk in ((0, 's'), (1, '^')):
                v = np.array([q[j + off] for q in gr_rows]); m = v > 0
                ax[1].semilogy(L[m], v[m], mk + '--', label=f'R={R} {"Dir" if off==0 else "Neu"}')
        ax[1].set_xlabel('l'); ax[1].set_ylabel('growth rate sqrt(-E0) |M|'); ax[1].legend(fontsize=7)
        fig.tight_layout(); fig.savefig(os.path.join(HERE, 'truncated.png'), dpi=120)
    except Exception as e:
        print('plot failed', e)


if __name__ == '__main__':
    main()

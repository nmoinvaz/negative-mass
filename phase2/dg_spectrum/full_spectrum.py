"""
Task 1: spectrum of the DG operator H^ = -d^2/dx^2 + V^ on x in (0, inf), l = 2..8, |M| = 1.

Method.  Write Psi^ = r^3 u(r). Near x = 0, r ~ 2 x^(1/2) so r^3 ~ 8 x^(3/2): this factors out the
admissible (Friedrichs, unique s.a.) branch. The other branch x^(-1/2) ~ r^-1 becomes u ~ r^-4,
which a polynomial basis cannot represent, so the self-adjoint extension is selected automatically
and u is smooth at r = 0 (the u-equation has a regular singular point there with exponents 0, -7... ->
the regular one is picked). The u-equation (see common.py)
    -r u'' + (-6 + 2/(r+2)) u' + C(r) u = E r^3/(r+2)^2 u
is discretised with Chebyshev collocation on r = a (1+t)/(1-t), t in [-1,1] (algebraic map of
(0, inf); a sets the resolution scale / effective domain size), N+1 Gauss-Lobatto points, the
t = 1 (r = inf) node removed (u = 0 there; bound states decay exponentially, so this only affects
the discretised continuum). Generalised eigenproblem solved densely (scipy.linalg.eig).
Independent check: Sturm oscillation. The number of eigenvalues below 0 equals the number of zeros
in (0, inf) of the E = 0 solution that is admissible at x = 0 (Psi ~ x^(3/2)); we integrate it with
mpmath-precision ODE solver (scipy DOP853, rtol 1e-12) from small r outward and count sign changes.
"""
import sys, csv, os
import numpy as np, scipy.linalg as sl
from scipy.integrate import solve_ivp
import common as c

HERE = os.path.dirname(os.path.abspath(__file__))


def full_eigs(l, N, a):
    D, t = c.cheb(N)
    idx = np.arange(1, N + 1)
    t = t[idx]
    r = a * (1 + t) / (1 - t)
    drdt = 2 * a / (1 - t)**2
    Dt = D[np.ix_(idx, idx)]
    Dr = Dt / drdt[:, None]
    D2t = (D @ D)[np.ix_(idx, idx)]
    # d2/dr2 = (1/r_t^2) d2/dt2 - r_tt/r_t^3 d/dt ; r_tt = 4a/(1-t)^3
    rtt = 4 * a / (1 - t)**3
    Dr2 = D2t / drdt[:, None]**2 - (rtt / drdt**3)[:, None] * Dt
    n = c.numeric(l)
    A = -np.diag(r) @ Dr2 + np.diag(-6 + 2 / (r + 2)) @ Dr + np.diag(n['C'](r))
    B = np.diag(r**3 / (r + 2)**2)
    w, V = sl.eig(A, B)
    ok = np.isfinite(w)
    w, V = w[ok], V[:, ok]
    o = np.argsort(w.real)
    return w[o], V[:, o], r


def sturm_count(l, E=0.0, rmax=None):
    """Zeros of admissible solution at energy E in r>0 (E <= 0)."""
    n = c.numeric(l)
    # first-order system in r for Psi: Psi_x = f Psi_r ; use y = (Psi, Psi_x), d/dr = (1/f) d/dx
    def rhs(r, y):
        f = 1 + 2 / r
        return [y[1] / f, (n['Vh'](r) - E) * y[0] / f]
    r0 = 1e-3
    # series start: Psi ~ r^3 (1 + O(r)) ; use the w-equation local Frobenius to O(r)
    # u = 1 + u1 r with -r u'' + (-6+1)u' ... at leading order: -5 u1 + C(0) = 0 -> u1 = C(0)/5
    # (r->0: -r u''+(-6+2/(r+2))u' -> -5 u' at O(1))
    C0 = n['C'](1e-14)
    u1 = C0 / 5
    psi = r0**3 * (1 + u1 * r0)
    dpsidr = 3 * r0**2 + 4 * u1 * r0**3
    y0 = [psi, (1 + 2 / r0) * dpsidr]
    if rmax is None:
        rmax = 40.0 if E == 0 else 2.0 + 30.0 / np.sqrt(-E)   # beyond the turning point, monotone
    sol = solve_ivp(rhs, (r0, rmax), y0, method='DOP853', rtol=1e-12, atol=1e-30, dense_output=True)
    rr = np.geomspace(r0, rmax, 200001)
    ps = sol.sol(rr)[0]
    z = np.sum(np.sign(ps[1:]) != np.sign(ps[:-1]))
    return int(z), rr[np.where(np.sign(ps[1:]) != np.sign(ps[:-1]))[0]]


def main():
    out = []
    print("Full problem: Chebyshev (u = Psi^/r^3), negative eigenvalues vs exact -(k/2)^2\n")
    print(f"{'l':>2} {'N':>4} {'a':>5} {'#neg':>4} {'E0':>22} {'rel err':>9} {'E1 (lowest >=0)':>16} {'#complex':>8}")
    conv_rows = []
    best = {}
    for l in range(2, 9):
        k = c.kval(l); Eex = -(k / 2)**2
        for N in (60, 100, 160, 240):
            for a in (0.1, 0.3, 1.0, 3.0, 10.0):
                w, V, r = full_eigs(l, N, a)
                cplx = np.abs(w.imag) > 1e-8 * (1 + np.abs(w.real))
                wr = w.real[~cplx]
                neg = wr[wr < 0]
                E0 = neg.min() if neg.size else np.nan
                E1 = wr[wr >= 0].min() if np.any(wr >= 0) else np.nan
                rel = abs(E0 / Eex - 1)
                conv_rows.append([l, N, a, len(neg), E0, Eex, rel, E1, int(cplx.sum()), wr.min()])
                if N in (160, 240) or l == 2:
                    print(f"{l:>2} {N:>4} {a:>5} {len(neg):>4} {E0:>22.14f} {rel:>9.1e} {E1:>16.4e} {int(cplx.sum()):>8}")
        print()
    with open(os.path.join(HERE, 'full_convergence.csv'), 'w', newline='') as fh:
        wtr = csv.writer(fh)
        wtr.writerow(['l', 'N', 'a', 'n_neg', 'E0', 'E_exact', 'rel_err', 'E1_lowest_nonneg', 'n_complex', 'min_real'])
        wtr.writerows(conv_rows)

    print("Sturm oscillation check (zeros of admissible E=0 solution) and phi0 = 1/psi0 node r_s = 6/((l+2)(l-1))")
    for l in range(2, 9):
        z, where = sturm_count(l)
        E0 = -(c.kval(l) / 2)**2
        za, _ = sturm_count(l, E=E0 * 0.999)
        zb, _ = sturm_count(l, E=E0 * 1.001)
        print(f"l={l}: zeros at E=0: {z} at r={np.round(where, 4)};  at E=0.999*E0: {za};  at E=1.001*E0: {zb};  r_s={c.r_s(l):.4f}")

    # l = 2 eigenfunction vs closed form
    w, V, r = full_eigs(2, 200, 1.0)
    u = V[:, 0].real
    psi_num = r**3 * u
    psi_ex = c.exact_mode(2, r)
    m = np.argmax(np.abs(psi_ex)); psi_num *= psi_ex[m] / psi_num[m]
    inner = r < 50
    err = np.max(np.abs(psi_num - psi_ex)[inner]) / np.max(np.abs(psi_ex))
    err_all = np.max(np.abs(psi_num - psi_ex)) / np.max(np.abs(psi_ex))
    print(f"\nl=2 eigenfunction (N=200, a=1): E0 = {w[0].real:.15f}; max|Psi_num-Psi_exact|/max|Psi_exact| over nodes r<50: {err:.2e}"
          f" (all nodes incl. r~1.6e4 where r^3 amplifies u round-off: {err_all:.1e})")
    x = c.x_of_r(r)
    with open(os.path.join(HERE, 'l2_eigenfunction.csv'), 'w', newline='') as fh:
        wtr = csv.writer(fh); wtr.writerow(['r', 'x', 'psi_numeric', 'psi_exact'])
        for row in sorted(zip(r, x, psi_num, psi_ex)):
            if row[0] < 30: wtr.writerow(row)
    # also check with the formula as literally quoted (denominator 4r^3+12Mr^2+24M^3 at M=-1)
    psi_q = r**3 * (r + 2)**4 * np.exp(-2 * r) / (4 * r**3 - 12 * r**2 - 24)
    print("Quoted denominator 4r^3+12Mr^2+24M^3 (M=-1 -> 4r^3-12r^2-24) has a pole at r =",
          np.roots([4, -12, 0, -24])[np.isreal(np.roots([4, -12, 0, -24]))].real, "; correct denominator is r^3-3r^2+6 = r^3+3Mr^2-6M^3 (no positive root)")

    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        o = np.argsort(r); m = r[o] < 8
        ax[0].plot(r[o][m], psi_ex[o][m], 'k-', label='closed form')
        ax[0].plot(r[o][m], psi_num[o][m], 'r.', ms=3, label='Chebyshev N=200')
        ax[0].set_xlabel('r/|M|'); ax[0].set_ylabel('Psi^ (l=2)'); ax[0].legend()
        import collections
        for l in range(2, 9):
            rows = [q for q in conv_rows if q[0] == l and q[2] == 1.0]
            ax[1].semilogy([q[1] for q in rows], [max(q[6], 1e-16) for q in rows], 'o-', label=f'l={l}')
        ax[1].set_xlabel('N'); ax[1].set_ylabel('|E0/E_exact - 1|  (a=1)'); ax[1].legend(fontsize=7)
        fig.tight_layout(); fig.savefig(os.path.join(HERE, 'full_spectrum.png'), dpi=120)
    except Exception as e:
        print('plot failed', e)


if __name__ == '__main__':
    main()

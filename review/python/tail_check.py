"""Compare the paper's tabulated star (data) against our trusted isotropic solver
(/home/claude/negative-mass/phase2/mass_check/iso_star.py) INCLUDING beyond the seam,
and quantify the alpha value jump at the seam and the matter-mass error of the painted
Yukawa tail.  Usage: python3 -I tail_check.py STARS_DIR"""
import sys, math, numpy as np
sys.path.insert(0, '/home/claude/negative-mass/phase2/mass_check')
import iso_star as I

S = sys.argv[1]
cases = [('canonical_omega0.750.dat', 0.750, +1), ('phantom_omega0.7603.dat', 0.7603, -1)]
for fn, w, eps in cases:
    d = np.loadtxt(f'{S}/{fn}', comments='#'); r, phi, al = d.T
    phic = phi[0]
    i = np.where((phi < 0.05*phic) & (np.r_[0, np.diff(phi)] < 0))[0][0]
    rm = r[i]
    # tail extrapolated to r_m from rows i+1.. (Schwarzschild form with header-M)
    M_hdr = float(open(f'{S}/{fn}').readline().split('ADM_mass=')[1])
    a_tail_at_rm = (1 - M_hdr/(2*rm))/(1 + M_hdr/(2*rm))
    print(f'{fn}: alpha jump at seam: core {al[i]:.9f} tail-extrap {a_tail_at_rm:.9f} '
          f'diff {al[i]-a_tail_at_rm:+.2e}; next row tail {al[i+1]:.9f}')
    # trusted solution
    e = I.at_omega(w, eps, which='w_paper')
    sol = e['sol']; A = e['A']; Aal = e['Aal']; psi_inf_m = e['psi_inf']; al_inf_m = e['al_inf']
    # paper frame: r_phys = psi_inf(r_m)^2 r_int, psi/psi_inf(r_m), alpha/alpha_inf(r_m)
    rint_max = sol.t[np.argmin(np.abs(sol.y[0]) + (sol.t < 2)*1e9)]
    rr = np.linspace(0.05, min(rint_max*0.9, 25.0), 4000)
    Y = sol.sol(rr)
    r_phys = psi_inf_m**2*rr
    phi_true = Y[0]; al_true = Y[4]/al_inf_m; psi_true = Y[2]/psi_inf_m
    phi_tab = np.interp(r_phys, r, phi)
    al_tab = np.interp(r_phys, r, al)
    for R in (3, 5, rm, 6, 7, 8, 10, 12, 15):
        k = np.argmin(abs(r_phys-R))
        print(f'   r={r_phys[k]:6.2f}  phi_tab/phi_true={phi_tab[k]/phi_true[k]:.5f}  '
              f'alpha_tab-alpha_true={al_tab[k]-al_true[k]:+.2e}')
    # matter (rho psi^5 r^2) integrated: true vs tabulated phi in the tail region,
    # using the true metric (rough first-order estimate of how much tail matter the paint lacks)
    def rho(ph, dph, ps, a, ww):
        s = ph*ph; U = 0.5*s - 0.25*I.LAM*s*s + I.MU_PAPER*s**3/6
        return 0.5*(ww/a)**2*s + 0.5*dph**2/ps**4 + U
    wphys = e['w_paper']
    dphi_true = np.gradient(phi_true, r_phys); dphi_tab = np.gradient(phi_tab, r_phys)
    psi_p = psi_true
    m_true = 2*math.pi*np.cumsum(rho(phi_true, dphi_true, psi_p, al_true, wphys)*psi_p**5*r_phys**2*np.gradient(r_phys))
    m_tab = 2*math.pi*np.cumsum(rho(phi_tab, dphi_tab, psi_p, al_true, wphys)*psi_p**5*r_phys**2*np.gradient(r_phys))
    # (2 pi int rho psi^5 r^2 dr) = -r^2 psi' at r -> M/2; so M ~ 2*that
    out = r_phys > rm
    print(f'   source integral 4pi*int rho psi^5 r^2 dr beyond r_m: true {2*(m_true[-1]-m_true[out][0]):.3e}, '
          f'painted tail {2*(m_tab[-1]-m_tab[out][0]):.3e}; total (to r={r_phys[-1]:.1f}) true {2*m_true[-1]:.6f} painted {2*m_tab[-1]:.6f}')
    print(f"   M_paper-style {e['M_paper']:+.7f}  M_true {e['M_true']:+.7f}")

"""Compare the paper's tabulated stars (stars/*.dat; parsed as data only) with our
independent isotropic solution, and read back the mass their table encodes."""
import sys, math, numpy as np
sys.path.insert(0, '.')
import iso_star as I
D = sys.argv[1]
for fn, w_t, eps in [('canonical_omega0.750.dat', 0.750, +1), ('phantom_omega0.7603.dat', 0.7603, -1)]:
    hdr = open(f"{D}/{fn}").readline().strip()
    t = np.loadtxt(f"{D}/{fn}")
    r, ph, al = t[:, 0], t[:, 1], t[:, 2]
    phic = ph[0]
    # core/tail seam: grid spacing changes (core dr ~ 1/60*psi_inf^2, tail dr = 1/20)
    dr = np.diff(r); k = np.where(np.abs(dr[2:] - 0.05) < 1e-3)[0][0] + 2
    r_m = r[k]
    # mass encoded in the tail lapse alpha = (1-M/2r)/(1+M/2r)  =>  M = 2r(1-alpha)/(1+alpha)
    M_tail = 2*r[k+5:]*(1-al[k+5:])/(1+al[k+5:])
    print(fn, '\n  header:', hdr)
    print(f"  phi_c={phic:.10f}  seam r_m={r_m:.4f} phi(r_m)/phi_c={ph[k]/phic:.4f}  "
          f"M from tail lapse: {M_tail.min():+.8f}..{M_tail.max():+.8f}")
    e = I.at_omega(w_t, eps, which='w_paper')
    sol = e['sol']; s = e['psi_inf']**2
    # our phi0 on THEIR physical isotropic grid (core only): r_int = r_phys / psi_inf^2
    rr = r[1:k]
    ours = sol.sol(rr/s)[0]
    al_ours = sol.sol(rr/s)[4]/e['al_inf']
    psi_ours = sol.sol(rr/s)[2]/math.sqrt(s)
    rel = (ours - ph[1:k])/phic
    print(f"  ours phi_c={e['phic']:.10f} (diff {(e['phic']-phic)/phic:+.2e}); "
          f"max|dphi|/phi_c={np.abs(rel).max():.2e}; max|dalpha|={np.abs(al_ours-al[1:k]).max():.2e}")
    print(f"  areal radius of seam: psi^2 r = {np.interp(r_m, rr, psi_ours**2*rr):.4f}")
    # rms radius (|phi|-weighted, proper volume psi^6 r^2 dr; areal radius) from their table + our psi
    print(f"  M_paper-style(ours)={e['M_paper']:+.7f}  M_true(ours)={e['M_true']:+.7f}  "
          f"ratio header/true={float(hdr.split('ADM_mass=')[1])/e['M_true']:.5f}")

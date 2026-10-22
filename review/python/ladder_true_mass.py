"""True vs 5%-radius masses for the d=20 mass-ladder stars, with our trusted iso solver,
then refit the mass law a ~ M^p using the paper's own accelerations (mass_law.csv values
copied below as data)."""
import sys, math, numpy as np
sys.path.insert(0, '/home/claude/negative-mass/phase2/mass_check')
import iso_star as I
# (cell, omega_canon, omega_phantom, a_core, a_err, m_mean_paper) from campaign/mass_law.csv
rows = [("w0850", 0.85, 0.8590975998840956, 2.34522366e-05, 1.61353642e-06, 0.0091570754),
        ("d20",   0.75, 0.7603045190007871, 3.57037217e-05, 9.19942616e-07, 0.0143216395),
        ("w0700", 0.70, 0.7104849879354665, 4.73932105e-05, 6.41847605e-07, 0.0191728524),
        ("w0675", 0.675,0.6860720218606876, 5.56011543e-05, 1.44400141e-06, 0.0225412121)]
out = []
for name, wc, wp, a, ae, mpap in rows:
    ec = I.at_omega(wc, +1, which='w_paper'); ep = I.at_omega(wp, -1, which='w_paper')
    mt = 0.5*(ec['M_true'] - ep['M_true']); ms = 0.5*(ec['M_paper'] - ep['M_paper'])
    out.append((name, a, ae, mpap, ms, mt))
    print(f"{name}: M+ paper-style {ec['M_paper']:+.7f} true {ec['M_true']:+.7f} | M- paper-style {ep['M_paper']:+.7f} true {ep['M_true']:+.7f} | "
          f"Mbar csv {mpap:.7f} ours-5% {ms:.7f} true {mt:.7f} ratio csv/true {mpap/mt:.5f}; a d^2/Mbar csv {a*400/mpap:.4f} true {a*400/mt:.4f}", flush=True)
def wfit(m, a, ae):
    x, y, s = np.log(m), np.log(a), ae/a; w = 1/s**2
    W, Sx, Sy, Sxx, Sxy = w.sum(), (w*x).sum(), (w*y).sum(), (w*x*x).sum(), (w*x*y).sum()
    D = W*Sxx - Sx**2; p = (W*Sxy - Sx*Sy)/D; c = (Sxx*Sy - Sx*Sxy)/D; pe = math.sqrt(W/D)
    chi2 = (w*(y-p*x-c)**2).sum()
    if chi2/(len(m)-2) > 1: pe *= math.sqrt(chi2/(len(m)-2))
    return p, pe
a = np.array([r[1] for r in out]); ae = np.array([r[2] for r in out])
print("slope with csv masses   p = %.4f +- %.4f" % wfit(np.array([r[3] for r in out]), a, ae))
print("slope with true masses  p = %.4f +- %.4f" % wfit(np.array([r[5] for r in out]), a, ae))

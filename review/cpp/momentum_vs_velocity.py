# Compare packed per-sector field momentum px with |M_ADM| * v_core (v from the
# packed core track), as a consistency check of the momentum diagnostic.
import sys, numpy as np
d = np.loadtxt(sys.argv[1]); MP, MM = 0.014349964, 0.014294692
t = d[:,0]; xc, xp = d[:,1], d[:,5]
vc = np.gradient(xc, t); vp = np.gradient(xp, t)
def sm(a, n=25): return np.convolve(a, np.ones(n)/n, mode='same')
vc, vp = sm(vc), sm(vp)
for T in (25, 50, 100, 150, 190):
    i = np.argmin(abs(t-T))
    print(f"t={t[i]:6.1f} v_c={vc[i]:.4e} v_p={vp[i]:.4e} px_c/(M+ v_c)={d[i,11]/(MP*vc[i]):.3f} px_p/(|M-| v_p)={d[i,12]/(MM*vp[i]):.3f}")

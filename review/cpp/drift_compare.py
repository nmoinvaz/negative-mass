# Independent check: compare barycentre-midpoint vs core-midpoint drift and
# the per-sector matter momenta in the packed .dat files (data only).
import sys, numpy as np
camp = sys.argv[1]
for cell in sys.argv[2:]:
    b = np.loadtxt(f"{camp}/{cell}/sector_barycenters.dat")
    d = np.loadtxt(f"{camp}/{cell}/sector_dynamics.dat")
    tb, mb = b[:,0], 0.5*(b[:,2]+b[:,7])
    td, md = d[:,0], 0.5*(d[:,1]+d[:,5])
    def fa(t,x):
        m=(t>=5); return 2*np.polyfit(t[m],x[m],2)[0]
    print(cell)
    print("  bary drift %.4f  a=%.4e | core drift %.4f a=%.4e" % (mb[-1]-mb[0], fa(tb,mb), md[-1]-md[0], fa(td,md)))
    print("  t_end %.1f  px_c %.3e px_p %.3e px_tot %.3e  ratio px_c/px_p %.4f" % (td[-1], d[-1,11], d[-1,12], d[-1,13], d[-1,11]/d[-1,12]))
    print("  bary sep %.4f -> %.4f ; core sep %.4f -> %.4f" % (b[0,2]-b[0,7], b[-1,2]-b[-1,7], d[0,9], d[-1,9]))

"""Independent check of the r_m seam in the paper's packed star tables (data only).
Locates the seam (first row where phi < 5% phi_c), and reports jumps in phi, phi',
alpha, alpha' across it, plus the M encoded by the alpha tail."""
import sys, numpy as np
for path in sys.argv[1:]:
    hdr = open(path).readline()
    d = np.loadtxt(path, comments='#'); r, phi, al = d.T
    phic = phi[0]
    i = np.where((phi < 0.05*phic) & (np.r_[0, np.diff(phi)] < 0))[0][0]  # seam = last core row
    # core spacing vs tail spacing changes at the seam
    dr = np.diff(r)
    # left derivative from core rows, right from tail rows
    dphiL = (phi[i]-phi[i-1])/(r[i]-r[i-1]); dphiR = (phi[i+1]-phi[i])/(r[i+1]-r[i])
    dalL = (al[i]-al[i-1])/(r[i]-r[i-1]); dalR = (al[i+1]-al[i])/(r[i+1]-r[i])
    # second-order one-sided derivatives
    dphiL2 = (3*phi[i]-4*phi[i-1]+phi[i-2])/(2*(r[i]-r[i-1]))
    dphiR2 = (-3*phi[i]+4*phi[i+1]-phi[i+2])/(2*(r[i+1]-r[i]))
    dalL2 = (3*al[i]-4*al[i-1]+al[i-2])/(2*(r[i]-r[i-1]))
    dalR2 = (-3*al[i]+4*al[i+1]-al[i+2])/(2*(r[i+1]-r[i]))
    # M from the alpha tail at large r: alpha=(1-M/2r)/(1+M/2r) -> M = 2r(1-a)/(1+a)
    j = len(r)-1
    Mtail = 2*r[j]*(1-al[j])/(1+al[j])
    Mseam_tail = 2*r[i+1]*(1-al[i+1])/(1+al[i+1])
    # what M would the core lapse imply at the seam (Schwarzschild matching of value)
    Mcore_val = 2*r[i]*(1-al[i])/(1+al[i])
    # lapse 'mass' from slope: alpha' = M/r^2 /(1+M/2r)^2 approx
    print(path.split('/')[-1], hdr.strip()[-40:])
    print(f"  seam row {i}: r_m={r[i]:.5f} phi/phic={phi[i]/phic:.5f} dr core={r[i]-r[i-1]:.4g} dr tail={r[i+1]-r[i]:.4g}")
    print(f"  phi'  left {dphiL2:+.5e} right {dphiR2:+.5e}  rel jump {(dphiR2/dphiL2-1):+.3e}")
    print(f"  alpha' left {dalL2:+.5e} right {dalR2:+.5e}  rel jump {(dalR2/dalL2-1):+.3e}")
    print(f"  M from tail alpha (far) {Mtail:+.8f}; from tail at seam {Mseam_tail:+.8f}; from core alpha value at seam {Mcore_val:+.8f}")
    print(f"  r_max={r[-1]:.1f} phi(r_max)={phi[-1]:.2e}")

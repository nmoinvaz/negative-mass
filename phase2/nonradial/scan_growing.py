"""
Search for growing polar modes (s = gamma > 0 real) of the phantom and canonical stars:
sign changes of the real matching determinant D(gamma) on a log grid, refined by brentq; a root is
accepted only if |D| at the refined point is < 1e-6 of the bracket values (rejects jumps).
Usage: python scan_growing.py eps w l [npts] [gmin] [gmax]
Writes scan_eps{eps}_w{w}_l{l}.csv (gamma, D).
"""
import sys, os, numpy as np, time
from scipy.optimize import brentq
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem

eps, w, l = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
npts = int(sys.argv[4]) if len(sys.argv) > 4 else 80
gmin = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-3
gmax = float(sys.argv[6]) if len(sys.argv) > 6 else 300.0
rtol = float(os.environ.get('RTOL', '1e-10'))
phi_cut = float(os.environ.get('PHICUT', '1e-9'))

_P = None
def init():
    global _P
    st = B.solve_star(w, eps)
    _P = PolarProblem(st, l, rtol=rtol, phi_cut=phi_cut)
def Dg(g):
    return _P.D(g)

if __name__ == "__main__":
    T0 = time.time()
    grid = np.logspace(np.log10(gmin), np.log10(gmax), npts)
    with Pool(2, initializer=init) as pool:
        vals = np.array(pool.map(Dg, grid))
    init()
    tag = f"eps{eps:+d}_w{w}_l{l}" + (os.environ.get('TAG', ''))
    np.savetxt(os.path.join(HERE, f"scan_{tag}.csv"), np.c_[grid, vals], delimiter=',', header='gamma,D', comments='')
    print(f"{tag}: rm={_P.rm:.2f} M={_P.Mext:+.6f} rtol={rtol} phi_cut={phi_cut}  {npts} pts in [{gmin},{gmax}]  ({time.time()-T0:.0f}s)")
    nchg = 0
    for i in range(npts - 1):
        if np.sign(vals[i]) != np.sign(vals[i + 1]):
            nchg += 1
            x = brentq(Dg, grid[i], grid[i + 1], xtol=1e-12, rtol=1e-10)
            ok = abs(Dg(x)) < 1e-6 * max(abs(vals[i]), abs(vals[i + 1]))
            print(f"   sign change in [{grid[i]:.4g},{grid[i+1]:.4g}] -> gamma={x:.8g}  root-like={ok}")
    print(f"   sign changes: {nchg};  D range: min|D|={np.min(np.abs(vals)):.3e} at gamma={grid[np.argmin(np.abs(vals))]:.4g}; "
          f"sign at gmin {np.sign(vals[0]):+.0f}, at gmax {np.sign(vals[-1]):+.0f}")

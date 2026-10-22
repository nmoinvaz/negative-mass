"""Full FD spectrum (dense QZ) including COMPLEX sigma^2, which shooting along the real sigma^2 axis cannot
see. For eps = -1 the quadratic (kinetic) form of the reduced problem need not be definite, so complex
sigma^2 (oscillatory instability) is not excluded a priori. Reports all eigenvalues with |sigma^2| < 0.3
and |Im| > 1e-7, and the real ones below the continuum threshold, at two resolutions."""
import numpy as np, warnings, sys
warnings.filterwarnings('ignore')
import scipy.linalg as sl
from background import solve_star
from radial_modes import RadialProblem
import fd_check as F

def full_spectrum(P, n):
    """Capture the dense pencil assembled by fd_check.fd_eigs and return ALL finite eigenvalues."""
    cap = {}; orig = sl.eig
    def grab(a, b, *args, **kw):
        cap['a'] = a; cap['b'] = b
        return orig(a, b, *args, **kw)
    sl.eig = grab                     # fd_eigs calls scipy.linalg.eig through the module attribute
    try:
        F.fd_eigs(P, n)
    finally:
        sl.eig = orig
    v = orig(cap['a'], cap['b'], right=False)
    return v[np.isfinite(v)]

cases = [(+1, 0.75), (+1, 0.94), (-1, 0.7603), (-1, 0.85), (-1, 0.94)]
for eps, w in cases:
    st = solve_star(0.75, eps); cur = 0.75
    while abs(w - cur) > 1e-12:
        cur = round(cur + np.clip(w - cur, -0.02, 0.02), 6); st = solve_star(cur, eps, guess=st)
    P = RadialProblem(st)
    for n in (300, 600):
        v = full_spectrum(P, n)
        v = v[np.abs(v) < 0.3]
        cplx = v[np.abs(v.imag) > 1e-7]
        real_b = np.sort(v[(np.abs(v.imag) <= 1e-7) & (v.real < (1 - w)**2) & (np.abs(v) > 1e-8)].real)
        print(f"eps={eps:+d} w={w} n={n}: complex eigenvalues (|s2|<0.3): {np.array2string(cplx, precision=4) if cplx.size else 'none'}; "
              f"real bound sigma^2 below (1-w)^2: {np.array2string(real_b[:3], precision=6)}", flush=True)

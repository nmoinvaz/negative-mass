"""
Count of the poles of D in a box (#121): the zeros of the decaying Zerilli solution at the matching radius, Z(r_m; s).
polar_solver.decaying_plane builds the Zerilli column as (1, q/f) with q = f Z'/Z at r_m, so D has a pole at each zero
of Z(r_m; s) and winding2.py counts zeros of D minus these.

Z is integrated inward from a radius R_z FIXED for the box, in the WKB-scaled variables y = (Z, f Z') exp(Psi),
Psi(r) = int_{R_z}^r k/f, k = sqrt(s^2 + V) (principal branch, cut on Re s = 0):
    dy1/dr = (y2 + k y1)/f,   dy2/dr = ((s^2 + V) y1 + k y2)/f,   y(R_z) = (k^{-1/2}, -k^{1/2}),
for the decaying solution Z normalised by its WKB amplitude at R_z. The counted function is y1(r_m) = Z(r_m; s)
exp(Psi(r_m)). exp(Psi) is holomorphic and zero-free for Re s > 0, so y1(r_m) has the zeros of Z(r_m; s). It is
holomorphic in s (R_z does not depend on s) and varies slowly, and the other solution decays inward as exp(2 Psi). That
makes the system stiff at large R_z, so it is integrated with the implicit Radau method. R_z is the start radius
zerilli_q takes at the smallest |s| of the contour, r_m + zfac (40 + 20 l)/|s|_min + 50.
Z(conj s) = conj Z(s), so the symmetric box is doubled as in winding2.py, with the same step rule. Usage:

  python zerilli_zeros.py eps w l x0 x1 y1 [--y0=Y0] [--zfac=1] [--nseg=10] [--plant=RE,IM] [--cachedir=DIR]
                          [--flat] [--mutate=NAME]

Checks, each printed with PASS or FAIL (exit 1 if any fails):
  count     the winding is within 0.05 of an integer and every accepted step is at most 0.45 rad;
  plant     with --plant, Z (s - s0)(s - conj s0) adds the planted zeros inside the box (as winding2.py);
  holo      Cauchy-Riemann at three interior points: |dZ/dx + i dZ/dy| <= 1e-6 |dZ/dx|, central differences at h and
            h/2, Richardson-extrapolated so that the h^2 error of a holomorphic Z cancels, h = 1e-3 of the box's longer
            side and at most Re s / 2 (the stencil stays off the cut of k on Re s = 0);
  column    q = y2/y1 at r_m against polar_solver's zerilli_q(s) at the same three points, reported, and within 1e-6
            relative where Re s >= 0.03 (there both integrations reach the decaying solution; nearer the axis
            zerilli_q's per-point start radius leaves a WKB error that the inward integration does not damp).
  --flat    Mext = 0 (flat space, V = l(l+1)/r^2): Z is r k_l(s r) up to a zero-free factor, whose zeros all have
            Re s < 0, so every box must count 0, and q = 1/r + s k_l'(s r)/k_l(s r) at r_m exactly (check `flatq`,
            1e-8 relative at the three points). A planted-truth run of the same code path.
  --mutate=reverse  traverse the contour clockwise (the planted check must fail);
  --mutate=nonorm   start from y = (1, -k) at a radius that depends on |s| (zerilli_q's), so Z is not holomorphic
                    (the holo check must fail);
  --mutate=lshift   l + 1 in the potential (the column check must fail).
"""
import os, pickle, sys, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import spherical_kn

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opts = dict((a[2:].split('=', 1) + ['1'])[:2] for a in sys.argv[1:] if a.startswith('--'))
if set(opts) - {'y0', 'zfac', 'nseg', 'plant', 'cachedir', 'flat', 'mutate'} or len(args) != 6:
    sys.exit(__doc__)
eps, w, l = int(args[0]), float(args[1]), int(args[2])
x0, x1, y1 = (float(a) for a in args[3:6])
y0, zfac, nseg = float(opts.get('y0', 0.0)), float(opts.get('zfac', 1.0)), int(opts.get('nseg', 10))
mutate, flat, cachedir = opts.get('mutate'), 'flat' in opts, opts.get('cachedir', HERE)
plant = complex(*(float(a) for a in opts['plant'].split(','))) if 'plant' in opts else None
assert mutate in (None, 'reverse', 'nonorm', 'lshift'), mutate
assert x0 < x1 and 0 <= y0 < y1
sym = y0 == 0.0

st = B.solve_star(w, eps)
P = PolarProblem(st, l, rtol=1e-8)
rm, M = P.rm, (0.0 if flat else P.Mext)
lam = ((l - 1) * (l + 2) if mutate != 'lshift' else l * (l + 3)) / 2.0


def V(r):
    f = 1 - 2 * M / r
    return 2 * f * (lam**2 * r**2 * ((lam + 1) * r + 3 * M) + 9 * M**2 * (lam * r + M)) / (r**3 * (lam * r + 3 * M)**2)


if sym:
    corners = [complex(x1, 0), complex(x1, y1), complex(x0, y1), complex(x0, 0)]
else:
    corners = [complex(x1, y0), complex(x1, y1), complex(x0, y1), complex(x0, y0), complex(x1, y0)]
if mutate == 'reverse':
    corners = corners[::-1]
smin = min(abs(a + (b - a) * t) for a, b in zip(corners[:-1], corners[1:]) for t in np.linspace(0, 1, 401))
Rz = rm + zfac * (40.0 + 20.0 * l) / smin + 50.0


def Zrm(s):
    """(y1(r_m), q(r_m; s)) of the decaying solution, WKB-normalised at R_z."""
    R0 = (rm + zfac * (40.0 + 20.0 * l) / abs(s) + 50.0) if mutate == 'nonorm' else Rz
    k0 = np.sqrt(s * s + V(R0))
    ya = [1.0 + 0j, -k0] if mutate == 'nonorm' else [k0**-0.5, -k0**0.5]

    def A(rho):  # the complex 2x2 system as a real 4x4 one, (Re y1, Re y2, Im y1, Im y2); Radau is real only
        r = np.exp(rho); f = 1 - 2 * M / r; k = np.sqrt(s * s + V(r))
        a = (r / f) * np.array([[k, 1.0], [s * s + V(r), k]])
        return np.block([[a.real, -a.imag], [a.imag, a.real]])
    ya = np.array(ya, dtype=complex)
    sol = solve_ivp(lambda rho, y: A(rho) @ y, (np.log(R0), np.log(rm)), np.r_[ya.real, ya.imag], method='Radau',
                    jac=lambda rho, y: A(rho), rtol=1e-11, atol=1e-14)
    z, p = sol.y[0, -1] + 1j * sol.y[2, -1], sol.y[1, -1] + 1j * sol.y[3, -1]
    return z, p / z


tag = f"eps{eps:+d}_w{w}_l{l}_{x0}_{x1}_{y1}" + ("" if sym else f"_y0{y0}") + ("" if zfac == 1.0 else f"_z{zfac:g}") + \
      ("_flat" if flat else "") + (f"_m{mutate}" if mutate in ('nonorm', 'lshift') else "")
cf = os.path.join(cachedir, f"zerilli_zeros_cache_{tag}.pkl")
cache = pickle.load(open(cf, 'rb')) if os.path.exists(cf) else {}
T0 = time.time()


def Zc(s):
    s = complex(round(s.real, 12), round(s.imag, 12))
    if s not in cache:
        cache[s] = Zrm(s)[0]
        pickle.dump(cache, open(cf, 'wb'))
    return cache[s]


def count(F):
    """Phase of F along the contour with winding2.py's step rule: at most 0.4 rad per step, and a step of the first
    two bisection levels only if its two half steps are small too and add up to it."""
    total, maxstep, fmin = 0.0, 0.0, np.inf
    for a, b in zip(corners[:-1], corners[1:]):
        def step(t0, t1, depth):
            nonlocal maxstep, fmin
            da, db = F(a + (b - a) * t0), F(a + (b - a) * t1)
            dph = np.angle(db / da)
            fine = t1 - t0 <= 1e-5
            if abs(dph) <= 0.4 and depth <= 1 and not fine:
                dm = F(a + (b - a) * 0.5 * (t0 + t1))
                h1, h2 = np.angle(dm / da), np.angle(db / dm)
                if max(abs(h1), abs(h2)) > 0.4 or abs(h1 + h2 - dph) > 0.1:
                    return step(t0, 0.5 * (t0 + t1), depth + 1) + step(0.5 * (t0 + t1), t1, depth + 1)
            if abs(dph) > 0.4 and not fine:
                return step(t0, 0.5 * (t0 + t1), depth + 1) + step(0.5 * (t0 + t1), t1, depth + 1)
            maxstep, fmin = max(maxstep, abs(dph)), min(fmin, abs(da))
            return dph
        for t0, t1 in zip(np.linspace(0, 1, nseg + 1)[:-1], np.linspace(0, 1, nseg + 1)[1:]):
            total += step(t0, t1, 0)
    return (2 if sym else 1) * total / (2 * np.pi), maxstep, fmin


fails = []


def check(name, ok, msg):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)
    if not ok:
        fails.append(name)


print(f"{tag}: r_m = {rm:.2f}, Mext = {M:+.6f}, R_z = {Rz:.1f} (|s|_min = {smin:.4g})" +
      (f"  [mutate={mutate}]" if mutate else ""), flush=True)
N, maxstep, zmin = count(Zc)
print(f"RESULT {tag}: zeros of Z(r_m; s) in Re s∈[{x0},{x1}], " + (f"|Im s|≤{y1}" if sym else f"Im s∈[{y0},{y1}]") +
      f": winding = {N:+.4f}  (max phase step {maxstep:.2f} rad, min|Z| {zmin:.2e}, {len(cache)} evaluations, "
      f"{time.time() - T0:.0f} s)", flush=True)
if sym:
    print(f"   Z(x0) = {Zc(complex(x0, 0)).real:+.4e}, Z(x1) = {Zc(complex(x1, 0)).real:+.4e} (real, one sign)")
check("count", abs(N - round(N)) < 0.05 and maxstep <= 0.45 and (not sym or
      Zc(complex(x0, 0)).real * Zc(complex(x1, 0)).real > 0), f"winding {N:+.4f}, max step {maxstep:.2f} rad")
if flat:
    check("flat", round(N) == 0, f"flat space has no zero with Re s > 0, counted {N:+.4f}")

# three contour points, at 10 %, 50 % and 90 % of the perimeter
L = [abs(b - a) for a, b in zip(corners[:-1], corners[1:])]
pts = []
for frac in (0.1, 0.5, 0.9):
    t = frac * sum(L)
    for (a, b), Le in zip(zip(corners[:-1], corners[1:]), L):
        if t <= Le:
            pts.append(a + (b - a) * t / Le)
            break
        t -= Le
def cauchy_riemann(s, h):  # dZ/dx + i dZ/dy and dZ/dx by central differences, both O(h^2)
    dx = (Zrm(s + h)[0] - Zrm(s - h)[0]) / (2 * h)
    dy = (Zrm(s + 1j * h)[0] - Zrm(s - 1j * h)[0]) / (2 * h)
    return dx + 1j * dy, dx


cr, col = [], []
for fx, fy in ((0.25, 0.3), (0.5, 0.6), (0.75, 0.9)):
    s = complex(x0 + (x1 - x0) * fx, y0 + (y1 - y0) * fy)
    h = min(1e-3 * max(x1 - x0, y1 - y0), 0.5 * s.real)
    c1, _ = cauchy_riemann(s, h)
    c2, dx = cauchy_riemann(s, h / 2)
    cr.append(abs((4 * c2 - c1) / 3) / abs(dx))  # Richardson: the h^2 term of a holomorphic Z cancels
for s in pts:
    q = Zrm(s)[1]
    if flat:
        x = s * rm
        qz = 1 / rm + s * spherical_kn(l, x, derivative=True) / spherical_kn(l, x)
    else:
        qz = P.zerilli_q(s)
    col.append((s, abs(q - qz) / abs(qz)))
check("holo", max(cr) <= 1e-6, "Cauchy-Riemann defect " + ", ".join(f"{c:.1e}" for c in cr))
print(f"   q against {'r k_l(s r)' if flat else 'zerilli_q'}: " +
      ", ".join(f"s = {s.real:.4g}{s.imag:+.4g}i: {d:.1e}" for s, d in col))
if flat:
    check("flatq", all(d <= 1e-8 for s, d in col), "relative difference from the exact flat-space value")
else:
    check("column", all(d <= 1e-6 for s, d in col if s.real >= 0.03), "relative difference where Re s >= 0.03")

if plant is not None:
    inside = x0 < plant.real < x1 and ((0 < abs(plant.imag) < y1) if sym else (y0 < plant.imag < y1))
    expected = (2 if sym else 1) * inside
    pf = (lambda s: (s - plant) * (s - np.conj(plant))) if sym else (lambda s: s - plant)
    Np, mp, _ = count(lambda s: Zc(s) * pf(s))
    check("plant", abs(Np - N - expected) < 0.05 and mp < 0.45,
          f"s0 = {plant}, expected {expected}, counted {Np - N:+.4f} (max step {mp:.2f} rad)")
print(f"{'FAILED: ' + ', '.join(fails) if fails else 'all checks passed'}  ({time.time() - T0:.0f} s)")
sys.exit(1 if fails else 0)

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
                          [--flat] [--guard=1] [--axis] [--mutate=NAME]

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
  --guard=1 deepest bisection level at which a step is accepted only if its two half steps agree with it; `all`
            applies the half-step test at every level, as winding2.py's --guard. The default 1 is the rule of #54.
  --mutate=reverse  traverse the contour clockwise (the planted check must fail);
  --mutate=nonorm   start from y = (1, -k) at a radius that depends on |s| (zerilli_q's), so Z is not holomorphic
                    (the holo check must fail);
  --mutate=lshift   l + 1 in the potential (the column check must fail).

--axis (#123) counts a box that crosses Re s = 0 (an upper box, y0 > 0: Z(conj s) = conj Z(s) gives its mirror the same
count). The form above has k = sqrt(s^2 + V) on the principal branch, whose cut is on Re s = 0, so exp(Psi) and y1 are
not holomorphic across the axis. Here y = (Z, f Z') exp(s r*) (r* the tortoise coordinate) has no k in its equation,

    dy1/dr = (y2 + s y1)/f,   dy2/dr = ((s^2 + V) y1 + s y2)/f,

and exp(s r*) is entire and zero-free, so y1(r_m) = exp(s r_m*) Z(r_m; s) has the zeros of Z(r_m; s). Z is the Jost
solution Z = exp(-s r*) u, u = sum_n a_n r^-n, u(r) -> 1: the a_n follow from f (f u')' - 2 s f u' = V u by a recursion in
1/r with poles at s = 0 only, so y(R_z) = (u, f u' - s u) is holomorphic in the box, summed until the terms are below
1e-18 (|s| R_z >= 40 + 20 l there). y1(r_m) does not depend on R_z. For Re s < 0 the wanted solution grows outward, the
other decays outward and grows inward by exp(-2 Re s (R_z* - r_m*)), about e^16 at Re s = -0.001, so on the real axis
an error at R_z is amplified without bound and Radau needs a long time. Along the ray r = r_m - i p, p from R_z - r_m
to 0, the wanted solution decays outward as exp(-Im s p) and the other decays inward: the integration is stable for
every Im s > 0 and the result is the continuation of the decaying solution of Re s > 0. R_z is the box's, as above.
  holo      takes its three points at Re s = x0/2, 0 and x1/2 when the box crosses Re s = 0 (stencils on both sides of
            the axis), with h = 1e-3 of the longer side and no limit from Re s; the criterion stays 1e-6;
  column    unchanged: applied where Re s >= 0.03 only, so in a box near the axis it is applied to none of the 3 points
            (the line says how many), and the differences are only reported;
  --flat    the count must equal the number of roots in the box of the reverse Bessel polynomial (flat-space Z is
            exp(-s r) sum_n (l+n)!/(n!(l-n)! 2^n) (s r)^-n), 0 for a box near the axis;
  --mutate=principal  starts from y = (1, -sqrt(s^2 + V(R_z))), the principal branch (the holo check must fail);
  also printed: q against zerilli_q at the holo points, and y1 at the three contour points (the zfac comparison).
"""
import math, os, pickle, sys, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import spherical_kn

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opts = dict((a[2:].split('=', 1) + ['1'])[:2] for a in sys.argv[1:] if a.startswith('--'))
if set(opts) - {'y0', 'zfac', 'nseg', 'plant', 'cachedir', 'flat', 'guard', 'mutate', 'axis'} or len(args) != 6:
    sys.exit(__doc__)
eps, w, l = int(args[0]), float(args[1]), int(args[2])
x0, x1, y1 = (float(a) for a in args[3:6])
y0, zfac, nseg = float(opts.get('y0', 0.0)), float(opts.get('zfac', 1.0)), int(opts.get('nseg', 10))
mutate, flat, cachedir = opts.get('mutate'), 'flat' in opts, opts.get('cachedir', HERE)
plant = complex(*(float(a) for a in opts['plant'].split(','))) if 'plant' in opts else None
guard = np.inf if opts.get('guard', '1') == 'all' else int(opts.get('guard', '1'))
axis = 'axis' in opts
assert mutate in (None, 'reverse', 'nonorm', 'lshift', 'principal'), mutate
assert x0 < x1 and 0 <= y0 < y1
assert (mutate != 'principal' or axis) and (mutate != 'nonorm' or not axis), "--mutate=principal needs --axis, nonorm not"
assert not axis or y0 > 0, "--axis counts upper boxes (y0 > 0)"
if x0 <= 0 and not axis:
    sys.exit("the box reaches Re s <= 0, where k has its cut: use --axis")
sym = y0 == 0.0
crosses = x0 < 0 < x1

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


NSER = 400
if axis:  # V = x^2 f(x) W(x) with x = 1/r and W a power series in x
    Nc = np.zeros(NSER + 1); Nc[:4] = [lam**2 * (lam + 1), 3 * M * lam**2, 9 * M**2 * lam, 9 * M**3]
    Wc = 2 * np.convolve(Nc, [(j + 1) * (-3 * M / lam)**j / lam**2 for j in range(NSER + 1)])[:NSER + 1]


def jost_start(s, R):
    """y(R) = (u, f u' - s u) for Z = exp(-s r*) u, u = sum a_n R^-n, a_0 = 1, from the recursion of
    f (f u')' - 2 s f u' = V u in x = 1/r: 2 s (n+1) a_{n+1} = (W a)_n - n (n+1) a_n + 2 M (n-1)(n+1) a_{n-1}."""
    x = 1.0 / R
    a = np.zeros(NSER + 2, dtype=complex); a[0] = 1.0
    u, ux, small = 1.0 + 0j, 0j, 0
    for n in range(NSER):
        a[n + 1] = (np.dot(Wc[:n + 1], a[n::-1]) - n * (n + 1) * a[n] +
                    (2 * M * (n - 1) * (n + 1) * a[n - 1] if n else 0.0)) / (2 * s * (n + 1))
        term = a[n + 1] * x**(n + 1)
        u, ux = u + term, ux + (n + 1) * a[n + 1] * x**n
        small = small + 1 if abs(term) < 1e-18 * abs(u) else 0
        if small == 2:
            break
    else:
        raise RuntimeError(f"the series of u at R = {R} did not converge for s = {s}")
    return np.array([u, (1 - 2 * M / R) * (-ux / R**2) - s * u])


def Zrm_axis(s):
    """(y1(r_m), q(r_m; s)), y = (Z, f Z') exp(s r*), integrated along r = r_m - i p, p = e^sigma - 1, from
    p = R_z - r_m (Jost start) to p = 0; the wanted solution decays outward along the ray."""
    T = Rz - rm
    R0 = rm - 1j * T
    ya = np.array([1.0 + 0j, -np.sqrt(s * s + V(R0))]) if mutate == 'principal' else jost_start(s, R0)

    def A(sig):
        p = np.exp(sig) - 1.0; r = rm - 1j * p; f = 1 - 2 * M / r
        a = (-1j * (1.0 + p) / f) * np.array([[s, 1.0], [s * s + V(r), s]])
        return np.block([[a.real, -a.imag], [a.imag, a.real]])
    sol = solve_ivp(lambda sg, y: A(sg) @ y, (np.log(1.0 + T), 0.0), np.r_[ya.real, ya.imag], method='Radau',
                    jac=lambda sg, y: A(sg), rtol=1e-11, atol=1e-14)
    assert sol.success, f"Radau failed at s = {s}: {sol.message}"
    z, p = sol.y[0, -1] + 1j * sol.y[2, -1], sol.y[1, -1] + 1j * sol.y[3, -1]
    return z, p / z


Zfun = Zrm_axis if axis else Zrm
tag = f"eps{eps:+d}_w{w}_l{l}_{x0}_{x1}_{y1}" + ("" if sym else f"_y0{y0}") + ("" if zfac == 1.0 else f"_z{zfac:g}") + \
      ("_flat" if flat else "") + (f"_m{mutate}" if mutate in ('nonorm', 'lshift', 'principal') else "") + \
      ("_axis" if axis else "")
cf = os.path.join(cachedir, f"zerilli_zeros_cache_{tag}.pkl")
cache = pickle.load(open(cf, 'rb')) if os.path.exists(cf) else {}
T0 = time.time()


def Zc(s):
    s = complex(round(s.real, 12), round(s.imag, 12))
    if s not in cache:
        cache[s] = Zfun(s)[0]
        pickle.dump(cache, open(cf, 'wb'))
    return cache[s]


def count(F):
    """Phase of F along the contour with winding2.py's step rule: at most 0.4 rad per step, and a step of the first
    two bisection levels (every level with --guard=all) only if its two half steps are small too and add up to it."""
    total, maxstep, fmin = 0.0, 0.0, np.inf
    for a, b in zip(corners[:-1], corners[1:]):
        def step(t0, t1, depth):
            nonlocal maxstep, fmin
            da, db = F(a + (b - a) * t0), F(a + (b - a) * t1)
            dph = np.angle(db / da)
            fine = t1 - t0 <= 1e-5
            if abs(dph) <= 0.4 and depth <= guard and not fine:
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
      ("  [axis: y = (Z, f Z') exp(s r*) along r = r_m - i p]" if axis else "") +
      (f"  [mutate={mutate}]" if mutate else ""), flush=True)
N, maxstep, zmin = count(Zc)
print(f"RESULT {tag}: zeros of Z(r_m; s) in Re s∈[{x0},{x1}], " + (f"|Im s|≤{y1}" if sym else f"Im s∈[{y0},{y1}]") +
      f": winding = {N:+.4f}  (max phase step {maxstep:.2f} rad, min|Z| {zmin:.2e}, {len(cache)} evaluations, "
      f"{time.time() - T0:.0f} s)", flush=True)
if sym:
    print(f"   Z(x0) = {Zc(complex(x0, 0)).real:+.4e}, Z(x1) = {Zc(complex(x1, 0)).real:+.4e} (real, one sign)")
check("count", abs(N - round(N)) < 0.05 and maxstep <= 0.45 and (not sym or
      Zc(complex(x0, 0)).real * Zc(complex(x1, 0)).real > 0), f"winding {N:+.4f}, max step {maxstep:.2f} rad")
if flat and axis:
    # roots in s of sum_n c_n (s r_m)^-n, c_n = (l+n)!/(n!(l-n)! 2^n), the zeros of flat-space Z(r_m; s)
    c = [math.factorial(l + n) / (math.factorial(n) * math.factorial(l - n) * 2.0**n) for n in range(l + 1)]
    sroots = [1.0 / (tr * rm) for tr in np.roots(c[::-1])]
    nin = sum(x0 < z.real < x1 and y0 < z.imag < y1 for z in sroots)
    check("flat", abs(N - nin) < 0.05, f"flat space has {nin} root(s) of the reverse Bessel polynomial in the box "
          f"(of {', '.join(f'{z.real:.4g}{z.imag:+.4g}i' for z in sroots)}), counted {N:+.4f}")
elif flat:
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
    dx = (Zfun(s + h)[0] - Zfun(s - h)[0]) / (2 * h)
    dy = (Zfun(s + 1j * h)[0] - Zfun(s - 1j * h)[0]) / (2 * h)
    return dx + 1j * dy, dx


def qref(s):  # the column's q at r_m from the reference: exact flat-space value, or polar_solver's zerilli_q
    if flat:
        x = s * rm
        return 1 / rm + s * spherical_kn(l, x, derivative=True) / spherical_kn(l, x)
    return P.zerilli_q(s)


cr, col, crpts = [], [], []
for k, (fx, fy) in enumerate(((0.25, 0.3), (0.5, 0.6), (0.75, 0.9))):
    s = complex((0.5 * x0, 0.0, 0.5 * x1)[k] if crosses else x0 + (x1 - x0) * fx, y0 + (y1 - y0) * fy)
    h = 1e-3 * max(x1 - x0, y1 - y0) if crosses else min(1e-3 * max(x1 - x0, y1 - y0), 0.5 * abs(s.real))
    crpts.append(s)
    c1, _ = cauchy_riemann(s, h)
    c2, dx = cauchy_riemann(s, h / 2)
    cr.append(abs((4 * c2 - c1) / 3) / abs(dx))  # Richardson: the h^2 term of a holomorphic Z cancels
qpts = []
for s in pts:
    qz = qref(s)
    qpts.append(Zfun(s)[1])
    col.append((s, abs(qpts[-1] - qz) / abs(qz)))
check("holo", max(cr) <= 1e-6, "Cauchy-Riemann defect " + ", ".join(f"{c:.1e}" for c in cr))
print(f"   q against {'r k_l(s r)' if flat else 'zerilli_q'}: " +
      ", ".join(f"s = {s.real:.4g}{s.imag:+.4g}i: {d:.1e}" for s, d in col))
if flat:
    check("flatq", all(d <= 1e-8 for s, d in col), "relative difference from the exact flat-space value")
else:
    check("column", all(d <= 1e-6 for s, d in col if s.real >= 0.03), "relative difference where Re s >= 0.03" +
          (f" ({sum(s.real >= 0.03 for s, _ in col)} of 3 points)" if axis else ""))
if axis:  # reported only: q at the holo points on both sides of the axis, and y1 and q at the points
    qhol = [Zfun(s)[1] for s in crpts]
    print("   q against " + ('r k_l(s r)' if flat else 'zerilli_q') + " at the holo points: " +
          ", ".join(f"s = {s.real:.4g}{s.imag:+.4g}i: {abs(q - qref(s)) / abs(qref(s)):.1e}" for s, q in zip(crpts, qhol)))
    print("   y1 at the three points: " + ", ".join(f"{Zc(s).real:.17e}{Zc(s).imag:+.17e}j" for s in pts))
    for name, ss, qq in (("contour", pts, qpts), ("holo", crpts, qhol)):
        print(f"   q at the {name} points: " + ", ".join(f"{s!r} {q.real:.17e}{q.imag:+.17e}j" for s, q in zip(ss, qq)))

if plant is not None:
    inside = x0 < plant.real < x1 and ((0 < abs(plant.imag) < y1) if sym else (y0 < plant.imag < y1))
    expected = (2 if sym else 1) * inside
    pf = (lambda s: (s - plant) * (s - np.conj(plant))) if sym else (lambda s: s - plant)
    Np, mp, _ = count(lambda s: Zc(s) * pf(s))
    check("plant", abs(Np - N - expected) < 0.05 and mp < 0.45,
          f"s0 = {plant}, expected {expected}, counted {Np - N:+.4f} (max step {mp:.2f} rad)")
print(f"{'FAILED: ' + ', '.join(fails) if fails else 'all checks passed'}  ({time.time() - T0:.0f} s)")
sys.exit(1 if fails else 0)

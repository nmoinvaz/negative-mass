"""
#129: is the failure of #122's harmonic test of arg D (RESULTS.md §5.5, `d_factors.py --part=harm`) integrator error?

At the three harm points of each of the four boxes with Re s >= 0.1 (p7603_l2, p7603_l3, c075_l2, p085_l3; the middle
point is also `mean`'s), phi = arg D is tested by the mean value property on N = 64 points of circles of radius
rho = 0.05 L and 0.025 L (L = 2.5):
    lap = 4 mean(phi_k - phi_0) / rho^2,   |grad phi| from the first Fourier mode,
    pass  iff  |lap| L^2 <= 1e-6 (1 + |grad phi| L)          (§5.5's criterion, unchanged).
Unlike the five-point stencil, the mean does not amplify a point error eps by 1/h^2: an error eps per point moves
|lap| L^2 by at most 4 eps L^2 / rho^2 (1600 eps at 0.05 L, 6400 eps at 0.025 L).

D is winding2.py's D (`PolarProblem.D` with winding2's fixed R_s = r_m + 40 / min Re k over the box's left and top edges,
times winding2's zero-free phase factor), at three settings:
    S8   rtol 1e-8,  nchunks as polar_solver chooses it (per point), zerilli_q at its own rtol 1e-11  (winding2's D);
    S10  rtol 1e-10, nchunks fixed across each circle at the centre's value,    zerilli_q at rtol 1e-12;
    S12  rtol 1e-12, nchunks fixed across each circle,                          zerilli_q at rtol 1e-13.
The fixed nchunks and the zerilli_q tolerance are set on a subclass of PolarProblem defined here; polar_solver.py is
not changed. zerilli_q is copied from polar_solver.py with its rtol made a parameter.

Controls, at every point and radius, through the same circle code:
    planted holomorphic g = (s - a)(s - c) e^{2s}, a = s0 + 0.9 + 0.4i, c = s0 - 0.7 - 0.6i: |lap| L^2 <= 1e-3 of the limit;
    planted non-harmonic g exp(i 1e-3 (Re s)^2 / L^2): |lap| L^2 above the limit.
Mutations: --mutate=nonholo multiplies D by exp(i 1e-3 (Re s)^2 / L^2), so D must fail at every point;
           --mutate=half takes the mean over half the circle, so the planted holomorphic control must fail.

Usage, from phase2/nonradial:
    python harm_circles.py --box=c075_l2 --setting=S12 --raw=DIR [--points=0,1,2] [--n=64] [--mutate=nonholo|half]
        evaluates (cached in DIR/harm_circles_<box>_<setting>.pkl after every point) and prints one line per point and
        radius. Exit 1 if a control fails or D fails at any point and radius of this run.
    python harm_circles.py --collect --raw=DIR
        reads every cache in DIR (no new evaluation) and prints the table of all boxes, points, settings and radii, and
        which setting first passes each point. Exit 1 if a control fails, a cache is missing, or D fails anywhere at S12.
"""
import os, sys, pickle, time
import numpy as np
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem

BOXES = {'p7603_l2': (-1, 0.7603, 2), 'p7603_l3': (-1, 0.7603, 3), 'c075_l2': (1, 0.75, 2), 'p085_l3': (-1, 0.85, 3)}
X0, X1, Y1, RSMULT = 0.1, 2.0, 2.5, 40.0            # every box: Re s in [0.1, 2], |Im s| <= 2.5 (collect_i126.py --jobs)
L = max(X1 - X0, Y1)
POINTS = [complex(X0 + (X1 - X0) * fx, Y1 * fy) for fx, fy in zip((0.25, 0.5, 0.75), (0.3, 0.6, 0.9))]  # d_factors harm_points
SETTINGS = {'S8': (1e-8, False, None), 'S10': (1e-10, True, 1e-12), 'S12': (1e-12, True, 1e-13)}
FRACS = (0.05, 0.025)
TOL = 1e-6


class FixedPolar(PolarProblem):
    """PolarProblem with nchunks fixed at a reference s and zerilli_q at a chosen rtol (both None: unchanged)"""
    s_ref = None
    zq_rtol = None

    def _integrate_plane(self, Y, ra, rb, s, nchunks=None, Jfun=None):
        if nchunks is None and self.s_ref is not None:
            nchunks = int(max(20, abs(rb - ra) * max(1.0, abs(self.s_ref)) / 2.0))
        return super()._integrate_plane(Y, ra, rb, s, nchunks=nchunks, Jfun=Jfun)

    def zerilli_q(self, s, r_start=None):
        if self.zq_rtol is None:
            return super().zerilli_q(s, r_start)
        M, lam = self.Mext, (self.l - 1) * (self.l + 2) / 2.0            # copied from polar_solver.zerilli_q

        def V(r):
            f = 1 - 2 * M / r
            return 2 * f * (lam ** 2 * r ** 2 * ((lam + 1) * r + 3 * M) + 9 * M ** 2 * (lam * r + M)) / (r ** 3 * (lam * r + 3 * M) ** 2)
        if r_start is None:
            r_start = self.rm + 40.0 / abs(s) + 20.0 * self.l / abs(s) + 50.0
        q0 = -np.sqrt(s * s + V(r_start))

        def rhs(rho, q):
            r = np.exp(rho); f = 1 - 2 * M / r
            return r * (s * s + V(r) - q * q) / f
        if np.iscomplexobj(q0) or np.iscomplexobj(s):
            q0 = -np.sqrt(complex(s * s + V(r_start)))
            if (q0 / s).real > 0:
                q0 = -q0
        sol = solve_ivp(rhs, (np.log(r_start), np.log(self.rm)), [q0], method='DOP853', rtol=self.zq_rtol, atol=1e-14)
        return sol.y[0, -1] if np.iscomplexobj(sol.y) else float(sol.y[0, -1])


def kpm(s, w):
    return np.sqrt(1 + (s - 1j * w) ** 2), np.sqrt(1 + (s + 1j * w) ** 2)


def make_D(name, setting):
    eps, w, l = BOXES[name]
    rtol, fixed, zq = SETTINGS[setting]
    P = FixedPolar(B.solve_star(w, eps), l, rtol=rtol)
    P.zq_rtol = zq
    cands = [X0 + 1j * y for y in np.linspace(0.0, Y1, 200)] + [x + 1j * Y1 for x in np.linspace(X0, X1, 200)]
    P.Rs_fixed = P.rm + RSMULT / min(min(k.real for k in kpm(c, w)) for c in cands)     # winding2's rule, symmetric box

    def D(s, s_ref):
        P.s_ref = s_ref if fixed else None
        a, b = kpm(s, w)
        return P.D(s) * np.exp(-1j * ((a + b) * (P.Rs_fixed - P.rm) + (a + b + s) * P.rm).imag)
    return D, P


def circle(F, s0, rho, n, half):
    th = 2 * np.pi * np.arange(n) / n
    f0 = F(s0)
    ph = np.array([np.angle(F(s0 + rho * np.exp(1j * t)) / f0) for t in th])
    assert np.abs(ph).max() < 2.5, "phase change across the circle too large for the principal branch"
    sel = slice(0, n // 2) if half else slice(None)
    lap = 4 * ph[sel].mean() / rho ** 2
    grad = np.hypot(2 * np.mean(ph * np.cos(th)), 2 * np.mean(ph * np.sin(th))) / rho
    return abs(lap) * L ** 2, TOL * (1 + grad * L), grad


def controls(s0, rho, n, half):
    a, c = s0 + 0.9 + 0.4j, s0 - 0.7 - 0.6j
    good = lambda s: (s - a) * (s - c) * np.exp(2 * s)
    bad = lambda s: good(s) * np.exp(1j * 1e-3 * s.real ** 2 / L ** 2)
    lg, rg, _ = circle(good, s0, rho, n, half)
    lb, rb, _ = circle(bad, s0, rho, n, half)
    return lg, rg, lb, rb


def opt(name, default=None):
    return next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith(f'--{name}=')), default)


def run():
    name, setting, raw = opt('box'), opt('setting'), opt('raw')
    pts = [int(i) for i in opt('points', '0,1,2').split(',')]
    n, mut = int(opt('n', '64')), opt('mutate')
    assert name in BOXES and setting in SETTINGS and raw and mut in (None, 'nonholo', 'half') and n % 4 == 0
    cf = os.path.join(raw, f"harm_circles_{name}_{setting}" + ("" if n == 64 else f"_n{n}") + ".pkl")
    cache = pickle.load(open(cf, 'rb')) if os.path.exists(cf) else {}
    D, P = make_D(name, setting)
    eps, w, l = BOXES[name]
    print(f"# harm_circles {name} (eps {eps:+d}, w {w}, l {l}) {setting} {SETTINGS[setting]}, N = {n}, L = {L}, "
          f"R_s = {P.Rs_fixed:.1f}, r_m = {P.rm:.2f}" + (f", mutate={mut}" if mut else ""), flush=True)
    fails = []
    for i in pts:
        s0 = POINTS[i]
        t0, new = time.time(), 0

        def F(s):
            nonlocal new
            key = (i, complex(round(s.real, 13), round(s.imag, 13)))
            if key not in cache:
                cache[key] = D(s, s0); new += 1
            v = cache[key]
            return v * np.exp(1j * 1e-3 * s.real ** 2 / L ** 2) if mut == 'nonholo' else v
        for frac in FRACS:
            rho = frac * L
            ld, rd, gd = circle(F, s0, rho, n, mut == 'half')
            lg, rg, lb, rb = controls(s0, rho, n, mut == 'half')
            ok = (ld <= rd, lg <= 1e-3 * rg, lb > rb)
            print(f"  point {i} s0 = {s0.real:.4g}{s0.imag:+.4g}i rho = {frac} L: D |lap| L^2 {ld:.3e}, limit {rd:.3e}, "
                  f"ratio {ld / rd:.3e}, |grad| {gd:.3g}: {'PASS' if ok[0] else 'FAIL'}; planted holo {lg:.1e} "
                  f"{'PASS' if ok[1] else 'FAIL'}; planted non-harmonic {lb:.1e} {'PASS' if ok[2] else 'FAIL'}", flush=True)
            for nm, good in zip(('D', 'holo', 'nonharm'), ok):
                if not good:
                    fails.append(f"{nm} at point {i}, rho {frac} L")
        pickle.dump(cache, open(cf + '.tmp', 'wb')); os.replace(cf + '.tmp', cf)
        print(f"  point {i}: {new} new evaluations, {(time.time() - t0) / max(new, 1):.1f} s each", flush=True)
    print("RESULT", "FAIL: " + "; ".join(fails) if fails else "PASS")
    return 1 if fails else 0


def collect():
    raw, n = opt('raw'), 64
    bad, rows, missing = [], [], []
    print(f"# harm_circles collect: N = {n}, L = {L}, criterion |lap| L^2 <= {TOL:g} (1 + |grad| L); points {POINTS}")
    print("| box | point | s0 | setting | rho / L | abs(lap) L^2 | limit | ratio | verdict |")
    print("|---|---|---|---|---|---|---|---|---|")
    first = {}
    for name in BOXES:
        for setting in SETTINGS:
            cf = os.path.join(raw, f"harm_circles_{name}_{setting}.pkl")
            if not os.path.exists(cf):
                missing.append(f"{name} {setting}"); continue
            cache = pickle.load(open(cf, 'rb'))
            for i, s0 in enumerate(POINTS):
                def F(s):
                    return cache[(i, complex(round(s.real, 13), round(s.imag, 13)))]
                try:
                    res = [circle(F, s0, frac * L, n, False) for frac in FRACS]
                except KeyError:
                    missing.append(f"{name} {setting} point {i}"); continue
                ok_all = True
                for frac, (ld, rd, _) in zip(FRACS, res):
                    ok = ld <= rd; ok_all &= ok
                    print(f"| {name} | {i} | {s0.real:.4g}{s0.imag:+.4g}i | {setting} | {frac} | {ld:.3e} | {rd:.3e} | "
                          f"{ld / rd:.3e} | {'pass' if ok else 'FAIL'} |")
                    lg, rg, lb, rb = controls(s0, frac * L, n, False)
                    if not (lg <= 1e-3 * rg and lb > rb):
                        bad.append(f"control at {name} point {i} rho {frac} L")
                    if setting == 'S12' and not ok:
                        bad.append(f"D fails at S12: {name} point {i} rho {frac} L")
                if ok_all and (name, i) not in first:
                    first[(name, i)] = setting
    print("\nfirst setting at which a point passes at both radii:")
    for name in BOXES:
        print(f"  {name}: " + ", ".join(f"point {i} {first.get((name, i), 'none')}" for i in range(len(POINTS))))
    if missing:
        bad.append("missing: " + "; ".join(missing))
    print("RESULT", "FAIL: " + "; ".join(bad) if bad else "PASS")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(collect() if '--collect' in sys.argv else run())

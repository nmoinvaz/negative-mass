"""
Checks that the Zerilli column is D's only source of poles in the counted boxes (#122). D is polar_solver's
PolarProblem.D with winding2.py's phase factor (Dc). Usage, from phase2/nonradial:

  python d_factors.py --part=factors,gap,norm,pole,axis --raw=RAW_FOLDER [--cache=DIR] [--boxes=a,b] [--rtol=1e-8]
                      [--stop=SECONDS] [--maxeval=N] [--mutate=NAME]

RAW_FOLDER (negative-mass-raw/I122-rpi5) holds this script's evaluation caches, one pickle per part and box, so a rerun
replays from them. Its siblings I126-rpi5 (the box caches, jobs.txt, the winding2 logs) and I123-rpi5 (the n3 logs) are
read, never written (--cache=DIR overrides I126-rpi5/cache). The 61 boxes are the lines of `collect_i126.py --jobs`, with
winding2.py's Rs_fixed (rsmult 40, or 80 where the job says so; a Z job takes its D twin's). Every check prints PASS or FAIL
and the script exits 1 when one fails. Criteria, fixed before the runs:

  factors  the table of D's factors, then
           moncrief  |det T| / (|T col 1| |T col 2|) >= 1e-3 at r_m (T is the Moncrief map, constant in s);
           jpoly     J(r; s^2) is linear in s^2 and even in s at radii in [r0, r_m] and in vacuum (defect <= 1e-10);
           harm      phi = arg D (winding2's Dc) is harmonic: at three interior points of each of six boxes (p7603_l2,
                     p7603_l3, c075_l2, p085_l3, n1_l4 and the upper box n3_l2) the five-point Laplacian at h and h/2,
                     Richardson-extrapolated (4 L(h/2) - L(h)) / 3, with h = 1e-3 of the box's longer side L, satisfies
                     |lap phi| L^2 <= 1e-6 (1 + |grad phi| L). Planted: the phase of (s - a)(s - b) e^s, a and b outside
                     the box, passes; that phase plus 1e-3 (Re s)^2 / L^2 fails.
  gap      per box, over the contour (401 points per edge), the cache's points and a 20 x 20 grid, polished by a local
           search inside the box: min over the two selected eigenvalues of the vacuum Jacobian's scalar block at R_s against
           the other two, |difference| / max |eigenvalue| >= 1e-3 (gap); the gap in real parts, relative the same way, >= 1e-3
           (order: the selection by real part is stable); sigma_min / sigma_max of the projector's columns 0 and 2 >= 1e-3
           (cols: the start columns of the plane are independent). Planted: boxes around s = i(1 - w) and i(1 + w) must fail
           gap and order; a projector whose columns 0 and 2 coincide must fail cols. The branch points s = +-i(1 +- w) are
           reported with each box's distance to them.
  norm     D with each eigenvector of np.linalg.eig in scalar_plane multiplied by a random nonzero complex factor at each s
           differs from D by at most 1e-10 relative, at five points in each of the six boxes above.
  pole     in n1_l2 a pole 1/((s - s0)(s - conj s0)), s0 = 0.05 + 0.1i, multiplied onto the Zerilli column (zv in
           decaying_plane), then onto the scalar plane's first start column, on a copy of those methods on the script's own
           instance, recounted with winding2.py's step rule (copied) at --guard=all: the count must change by -2 against
           #126's 0 (pole_z, pole_s), D_pole / D must equal c / |c| (c the factor) at every point the cache holds
           (<= 1e-8), and the copy with factor 1 must reproduce the cached D (<= 1e-12, `clone`). `counter`: the copied
           counter run on the cached D reproduces #126's winding, max step, min |D| and point count and returns the +2
           of the planted check.
  axis     in n3_l2 and n3_l3, D with q from zerilli_zeros.py --axis (its function copied) in place of zerilli_q at every
           cached contour point: axis_q the copied q reproduces the six q values of the I123 logs (<= 1e-10 relative),
           axis_ref the dense-chain winding of the cached D is 0, axis_wind that of D_axis is 0 and axis_plant that of
           (s - s0) D_axis with s0 = -0.0005 + 0.12i is +1 (chain increments below 1 rad), axisharm D_axis is harmonic at
           the three points of harm. Reported: the largest |D_axis - D| / |D| on each edge against min |D|.

  --mutate=NAME  (one part per run) breaks the part's input so that a check must fail:
    factors: nonholo (D times exp(i 1e-3 (Re s)^2 / L^2)), singularT (T row 2 = 0.5 row 1), jnonlin (J + 0.01 s^4)
    gap:     branch (every box moved to s = i(1 - w)), cols (the projector replaced by the planted bad one)
    norm:    noproj (the selected eigenvectors themselves as the plane's start columns)
    pole:    nophase (the R-diagonal factor of the planted column's QR dropped), nogauge (that column's phase fixed by
             its largest component, so no phase is left), reverse (the counter runs clockwise)
    axis:    reverse (the chain runs clockwise), wrongq (q times 1 + 1e-4)
  --boxes=a,b    restricts harm, norm and axis to these jobs (for mutation runs); --skip=static,harm leaves out the
  moncrief and jpoly checks, or harm; --plant=z|s|zs restricts pole to the
  Zerilli column or the scalar plane; --stop=SECONDS ends a run that reaches it (exit 3); --maxeval=N ends a part after N
  new D evaluations (the mutation runs of pole use it).
"""
import os, pickle, re, subprocess, sys, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem

T_START = time.time()
opts = dict((a[2:].split('=', 1) + ['1'])[:2] for a in sys.argv[1:] if a.startswith('--'))
PARTS = [p for p in opts.get('part', '').split(',') if p]
MUTS = {'factors': ('nonholo', 'singularT', 'jnonlin'), 'gap': ('branch', 'cols'), 'norm': ('noproj',),
        'pole': ('nophase', 'nogauge', 'reverse'), 'axis': ('reverse', 'wrongq')}
MUT = opts.get('mutate')
if (set(opts) - {'part', 'raw', 'cache', 'boxes', 'rtol', 'stop', 'maxeval', 'mutate', 'plant', 'skip'} or not PARTS or 'raw' not in opts
        or set(PARTS) - set(MUTS) or (MUT and (len(PARTS) != 1 or MUT not in MUTS[PARTS[0]]))):
    sys.exit(__doc__)
RAW = os.path.abspath(opts['raw']); os.makedirs(RAW, exist_ok=True)
SIB = os.path.dirname(RAW)
I126, I123 = os.path.join(SIB, 'I126-rpi5'), os.path.join(SIB, 'I123-rpi5')
CACHE = os.path.abspath(opts['cache']) if 'cache' in opts else os.path.join(I126, 'cache')
RTOL = float(opts.get('rtol', os.environ.get('RTOL', '1e-8')))
STOP = float(opts['stop']) if 'stop' in opts else None
MAXEVAL = int(opts['maxeval']) if 'maxeval' in opts else None
ONLY = opts['boxes'].split(',') if 'boxes' in opts else None
SKIP = opts.get('skip', '').split(',')              # factors: static (moncrief, jpoly) and/or harm
PLANTS = opts.get('plant', 'zs')                   # pole: which columns get the planted pole, z and/or s
THR = 1e-3                                    # gap, order, cols, moncrief
fails = []


def check(name, ok, msg):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)
    if not ok:
        fails.append(name)


class Stop(Exception):
    pass


STOPPED = []


# ---------------------------------------------------------------- the boxes
CONTROLS = {'z_flat', 'z_mut_reverse', 'z_mut_nonorm', 'z_mut_lshift', 'w_mut_reverse'}


class Box:
    def __init__(self, line):
        self.name, self.script, a, self.src = line.split('|')
        t = a.split()
        self.eps, self.w, self.l = int(t[0]), float(t[1]), int(t[2])
        self.x0, self.x1, self.y1 = float(t[3]), float(t[4]), float(t[5])
        o = dict(x[2:].split('=', 1) for x in t[6:] if x.startswith('--') and '=' in x)
        self.y0, self.zfac = float(o.get('y0', 0.0)), float(o.get('zfac', 1.0))
        self.rsmult = float(o.get('rsmult', 40.0))
        self.plant = complex(*(float(v) for v in o['plant'].split(',')))
        self.sym = self.y0 == 0.0

    def key(self):
        return (self.eps, self.w, self.l, self.x0, self.x1, self.y0, self.y1, self.zfac)

    def corners(self):
        if self.sym:
            return [complex(self.x1, 0), complex(self.x1, self.y1), complex(self.x0, self.y1), complex(self.x0, 0)]
        return [complex(self.x1, self.y0), complex(self.x1, self.y1), complex(self.x0, self.y1),
                complex(self.x0, self.y0), complex(self.x1, self.y0)]

    def cachefile(self):
        t = (f"eps{self.eps:+d}_w{self.w}_l{self.l}_{self.x0}_{self.x1}_{self.y1}" + ("" if self.sym else f"_y0{self.y0}")
             + ("" if self.rsmult == 40.0 else f"_rs{self.rsmult:g}") + ("" if self.zfac == 1.0 else f"_z{self.zfac:g}"))
        return os.path.join(CACHE, f"winding2_cache_{t}.pkl")


def load_boxes():
    out = subprocess.run([sys.executable, '-I', os.path.join(HERE, 'collect_i126.py'), '--jobs'], capture_output=True,
                         text=True, check=True).stdout
    boxes = [Box(ln) for ln in out.splitlines() if ln.strip() and ln.split('|')[0] not in CONTROLS]
    wtwin = {b.key(): b for b in boxes if b.name.startswith('w_')}
    for b in boxes:                           # a Z job's D twin sets Rs_fixed; s1_l5 has none (#115, rsmult 40)
        if b.name.startswith('z_') and b.key() in wtwin:
            b.rsmult = wtwin[b.key()].rsmult
    assert len(boxes) == 61, len(boxes)
    return boxes


BOXES = load_boxes()
BYNAME = {b.name: b for b in BOXES}
HARM6 = ['w_p7603_l2', 'w_p7603_l3', 'w_c075_l2', 'w_p085_l3', 'w_n1_l4', 'w_n3_l2']
NORM_EXTRA = ((0.1, 0.15), (0.9, 0.5))


def jobs(default):
    names = [n if n in BYNAME else 'w_' + n for n in ONLY] if ONLY else default
    return [BYNAME[n] for n in names]


_stars = {}


def star(eps, w):
    if (eps, w) not in _stars:
        _stars[(eps, w)] = B.solve_star(w, eps)
    return _stars[(eps, w)]


class DBox:
    """winding2.py's D for one box: PolarProblem at RTOL, Rs_fixed from the contour as winding2.py sets it, Dc's phase factor"""

    def __init__(self, box, rtol=None):
        self.box, self.rtol, self._P = box, rtol or RTOL, None

    @staticmethod
    def ks(w, s):
        return np.sqrt(1 + (s - 1j * w) ** 2), np.sqrt(1 + (s + 1j * w) ** 2)

    @staticmethod
    def rs_fixed(P, b, w):
        cands = [b.x0 + 1j * y for y in np.linspace(b.y0, b.y1, 200)] + [x + 1j * b.y1 for x in np.linspace(b.x0, b.x1, 200)]
        if not b.sym:
            cands += [b.x1 + 1j * y for y in np.linspace(b.y0, b.y1, 200)] + [x + 1j * b.y0 for x in np.linspace(b.x0, b.x1, 200)]
        kmin = min(min(k.real for k in DBox.ks(w, c)) for c in cands)
        return P.rm + b.rsmult / kmin

    @property
    def P(self):
        if self._P is None:
            b = self.box
            P = PolarProblem(star(b.eps, b.w), b.l, rtol=self.rtol)
            if b.zfac != 1.0:
                _zq = P.zerilli_q
                P.zerilli_q = lambda s, r_start=None: _zq(s, r_start=P.rm + b.zfac * (40.0 + 20.0 * b.l) / abs(s) + 50.0)
            P.Rs_fixed = self.rs_fixed(P, b, b.w)
            self._P = P
        return self._P

    def D(self, s):
        P = self.P
        kp, km = self.ks(self.box.w, s)
        return P.D(s) * np.exp(-1j * ((kp + km) * (P.Rs_fixed - P.rm) + (kp + km + s) * P.rm).imag)


class Store:
    """restart-safe cache of evaluations, one pickle per (part, box) in RAW"""

    def __init__(self, name):
        self.path = os.path.join(RAW, f"d_factors_cache_{name}.pkl")
        self.d = pickle.load(open(self.path, 'rb')) if os.path.exists(self.path) else {}
        self.new = 0

    def get(self, key, fn):
        if key in self.d:
            return self.d[key]
        if STOP is not None and time.time() - T_START > STOP:
            raise Stop(f"stop condition: {STOP:.0f} s reached")
        if MAXEVAL is not None and self.new >= MAXEVAL:
            raise Stop(f"{MAXEVAL} new evaluations reached")
        v = fn()
        self.d[key] = v
        self.new += 1
        pickle.dump(self.d, open(self.path + '.tmp', 'wb')); os.replace(self.path + '.tmp', self.path)
        return v


def box_cache(b):
    return pickle.load(open(b.cachefile(), 'rb'))


ckey = lambda s: complex(round(s.real, 12), round(s.imag, 12))


# ---------------------------------------------------------------- winding2.py's counter (copied)
def count_contour(F, corners, sym, guard=np.inf, nseg=10):
    """phase of F along the contour with winding2.py's step rule; returns (count, max step, min |F|)"""
    total, maxstep, nmin = 0.0, 0.0, np.inf
    for a, b in zip(corners[:-1], corners[1:]):
        def step(t0, t1, depth):
            nonlocal maxstep, nmin
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
            maxstep = max(maxstep, abs(dph)); nmin = min(nmin, abs(da))
            return dph
        for t0, t1 in zip(np.linspace(0, 1, nseg + 1)[:-1], np.linspace(0, 1, nseg + 1)[1:]):
            total += step(t0, t1, 0)
    return (2 if sym else 1) * total / (2 * np.pi), maxstep, nmin


# ---------------------------------------------------------------- the harmonic check
def harm_point(F, s0, L, h):
    """Richardson-extrapolated five-point Laplacian and gradient of phi = arg F at s0, and the criterion"""
    f0 = F(s0)

    def stencil(hh):
        ph = [np.angle(F(s0 + hh * m) / f0) for m in (1, -1, 1j, -1j)]
        return sum(ph) / hh ** 2, np.array([(ph[0] - ph[1]) / (2 * hh), (ph[2] - ph[3]) / (2 * hh)])
    l1, g1 = stencil(h); l2, g2 = stencil(h / 2)
    lap, grad = (4 * l2 - l1) / 3, (4 * g2 - g1) / 3
    lhs, rhs = abs(lap) * L ** 2, 1e-6 * (1 + np.linalg.norm(grad) * L)
    return lap, np.linalg.norm(grad), lhs, rhs


def box_L(b):
    return max(b.x1 - b.x0, b.y1 - b.y0)


def harm_points(b):
    ys = [b.y0 + (b.y1 - b.y0) * f for f in (0.3, 0.6, 0.9)]
    xs = (0.5 * b.x0, 0.0, 0.5 * b.x1) if b.x0 < 0 < b.x1 else [b.x0 + (b.x1 - b.x0) * f for f in (0.25, 0.5, 0.75)]
    return [complex(x, y) for x, y in zip(xs, ys)]


def planted_pair(b):
    """a pair of exact functions through the same stencil: holomorphic (must pass) and with a non-harmonic phase (must fail)"""
    L = box_L(b)
    a = complex(b.x1 + 0.5 * L, b.y1 + 0.5 * L); c = complex(b.x0 - 0.5 * L, b.y0 - 0.3 * L)
    good = lambda s: (s - a) * (s - c) * np.exp(s)
    bad = lambda s: good(s) * np.exp(1j * 1e-3 * s.real ** 2 / L ** 2)
    return good, bad


def harm_run(b, F, label, ref_pts=None):
    """the three points of box b with F, one line each; returns (all pass, worst ratio lhs/rhs)"""
    L = box_L(b); h = 1e-3 * L
    ok, worst = True, 0.0
    for s0 in (ref_pts or harm_points(b)):
        assert b.x0 < 0 < b.x1 or h <= 0.5 * s0.real, f"h = {h} too large at {s0}"
        lap, g, lhs, rhs = harm_point(F, s0, L, h)
        ok &= lhs <= rhs; worst = max(worst, lhs / rhs)
        print(f"    {label} s = {s0.real:+.5g}{s0.imag:+.5g}i: |grad phi| {g:.3g}, |lap phi| L^2 {lhs:.2e} against "
              f"{rhs:.2e}: {'ok' if lhs <= rhs else 'EXCEEDS'}", flush=True)
    return ok, worst


# ---------------------------------------------------------------- part factors
FACTORS = [
    ("1", "interior start data Y0 (regular_plane)", "constant, real, s-independent", "3 independent columns", "structure"),
    ("2", "J(r; s^2), interior and vacuum", "polynomial in s^2 (entire); singularities in r only", "-", "jpoly"),
    ("3", "interior flow Phi(r_m, r0; s), DOP853 in chunks", "entire: linear ODE with entire coefficients",
     "Phi invertible, so the plane keeps rank 3", "harm (numerical)"),
    ("4", "QR of the regular plane: Q, phase of R's diagonal", "D = det[Y, Z] / prod abs(R_ii), a positive divisor; arg D = arg det",
     "abs(R_ii) > 0", "pole (the factor enters as c / abs(c)), structure"),
    ("5", "Moncrief map T(r_m)", "constant in s (poles in r at 2 M_ext and -3 M_ext/lam, not at r_m)", "det T != 0", "moncrief"),
    ("6", "Zerilli column zv = (1, q/f), q = f Z'/Z", "meromorphic; poles at the zeros of Z(r_m; s) (the pole source, "
     "none in the 61 boxes: #121, #123, #124)", "first entry 1", "axis (D_axis against D), harm at n3_l2"),
    ("7", "scalar block Js = Jvac(R_s, s)[2:, 2:], R_s fixed per box", "entire (factor 2)", "-", "jpoly"),
    ("8", "eig(Js), the two eigenvalues of smallest real part", "eigenvectors carry arbitrary normalisations and the "
     "choice is by real part: neither is holomorphic by itself", "-", "norm, order"),
    ("9", "projector P = V_sel V^-1_sel", "holomorphic while the selected pair stays apart from the other two "
     "(Riesz projector); independent of the eigenvector normalisation", "-", "gap, norm"),
    ("10", "start columns P[:, 0], P[:, 2]", "holomorphic", "independent iff P e0 and P e2 are", "cols"),
    ("11", "scalar plane flow R_s -> r_m and its QR phases", "as 3 and 4", "as 3", "harm, pole_s"),
    ("12", "decaying plane Z = [m, S] and its QR", "m = T zv and S as above; positive divisor as 4", "m != 0 (T invertible, "
     "zv[0] = 1); S rank 2 (sigma_min(S) reported in norm)", "pole_z, norm"),
    ("13", "det[Y, Z] times the QR phases", "holomorphic matching determinant divided by a positive function", "zeros are "
     "the modes", "harm"),
    ("14", "winding2.py's exp(-i Im Theta), Theta = (k+ + k-)(R_s - r_m) + (k+ + k- + s) r_m",
     "a pure phase, not holomorphic; arg = -Im Theta = arg exp(-Theta), and exp(-Theta) is holomorphic and zero-free where "
     "k+- = sqrt(1 + (s -+ i w)^2) (principal) is, i.e. off the cut on Re s = 0, abs(Im s -+ w) >= 1", "never zero", "harm, gap"),
]


def linear_defect(Jfun):
    """Jfun(t) for t = s^2: defect from linearity in t, relative to |J|, over four values, plus evenness in s"""
    j0, j1 = Jfun(0.0), Jfun(1.0)
    worst = 0.0
    for t in (-1.0, 2.5 + 0.7j, 0.3 - 1.9j):
        worst = max(worst, np.abs(Jfun(t) - (j0 + t * (j1 - j0))).max() / np.abs(j1 - j0).max())
    return worst


def part_factors():
    print("# D's factors, in the order of polar_solver.D (D = det[Y, Z] py pz exp(-i Im Theta))")
    print("| # | factor | holomorphic in a box? | zero-free? | checked by |\n|---|---|---|---|---|")
    for r in FACTORS:
        print("| " + " | ".join(r) + " |")
    print()
    if 'static' in SKIP:
        return part_harm()
    seen = {}
    for b in BOXES:
        seen.setdefault((b.eps, b.w, b.l), b)
    # moncrief and jpoly, once per (star, l)
    worstT, worstJ = np.inf, 0.0
    for (eps, w, l), b in seen.items():
        P = DBox(b).P
        T = P.moncrief(P.rm)
        if MUT == 'singularT':
            T = T.copy(); T[1] = 0.5 * T[0]
        ratio = abs(np.linalg.det(T)) / (np.linalg.norm(T[:, 0]) * np.linalg.norm(T[:, 1]))
        worstT = min(worstT, ratio)
        J = P.J
        if MUT == 'jnonlin':
            J = lambda r, s, _J=P.J: _J(r, s) + 0.01 * s ** 4 * np.eye(6)
        d = 0.0
        for r in (0.02 * P.rm, 0.1 * P.rm, 0.35 * P.rm, 0.7 * P.rm, P.rm):
            d = max(d, linear_defect(lambda t, r=r: J(r, np.sqrt(complex(t)))))
            d = max(d, np.abs(J(r, 0.3 + 0.2j) - J(r, -0.3 - 0.2j)).max() / np.abs(J(r, 0.3 + 0.2j)).max())
        for r in (1.5 * P.rm, 3 * P.rm, 10 * P.rm):
            d = max(d, linear_defect(lambda t, r=r: P._Jvac(r, np.sqrt(complex(t)))))
        worstJ = max(worstJ, d)
        print(f"   eps {eps:+d}, w {w}, l {l}: r_m {P.rm:.3f}, M_ext {P.Mext:+.6f}, T conditioning {ratio:.3e}, "
              f"J linear-in-s^2 and even defect {d:.1e}", flush=True)
    check("moncrief", worstT >= THR, f"smallest |det T| / (|col 1| |col 2|) over the {len(seen)} (star, l) is {worstT:.3e} "
          f"(limit {THR:g}){'  [mutate=' + MUT + ']' if MUT else ''}")
    M0 = np.array([[1.0, 2.0], [0.5, -1.0]]); M1 = np.array([[0.3, 0.1], [-0.2, 0.7]])
    pg = linear_defect(lambda t: M0 + t * M1); pb = linear_defect(lambda t: M0 + t * M1 + 0.01 * t ** 2 * M1)
    check("jpoly", worstJ <= 1e-10 and pg <= 1e-10 and pb > 1e-10,
          f"worst defect {worstJ:.1e} (limit 1e-10); planted linear {pg:.1e} passes, planted with a t^2 term {pb:.1e} fails"
          f"{'  [mutate=' + MUT + ']' if MUT else ''}")
    if 'harm' not in SKIP and MUT in (None, 'nonholo'):
        part_harm()


def part_harm():
    print("# harm: phi = arg D harmonic, Richardson Laplacian at three interior points per box")
    allok, worst_all = True, 0.0
    for b in jobs(HARM6):
        L = box_L(b)
        good, bad = planted_pair(b)
        pg_ok, pg_w = harm_run(b, good, "planted holomorphic")
        pb_ok, pb_w = harm_run(b, bad, "planted non-harmonic")
        print(f"  {b.name}: L = {L:g}, h = {1e-3 * L:.3g}; planted holomorphic {'passes' if pg_ok else 'FAILS'} (worst "
              f"{pg_w:.2e} of the limit), planted non-harmonic {'fails as required' if not pb_ok else 'PASSES (wrong)'} "
              f"({pb_w:.2e})", flush=True)
        st = Store(f"factors_{b.name}"); db = DBox(b)
        mut = (lambda s, L=L: np.exp(1j * 1e-3 * s.real ** 2 / L ** 2)) if MUT == 'nonholo' else (lambda s: 1.0)
        F = lambda s, st=st, db=db, mut=mut: st.get(('D', s), lambda: db.D(s)) * mut(s)
        try:
            ok, worst = harm_run(b, F, "D")
        except Stop as e:
            print(f"  STOP {e}"); check("harm", False, f"stopped in {b.name}: {e}"); return
        print(f"  {b.name}: D {'passes' if ok else 'FAILS'}, worst |lap phi| L^2 / limit = {worst:.3e} "
              f"({st.new} new evaluations, {len(st.d)} stored)", flush=True)
        allok &= ok and pg_ok and not pb_ok; worst_all = max(worst_all, worst)
    check("harm", allok, f"Delta arg D = 0 in the boxes {', '.join(b.name for b in jobs(HARM6))}; worst ratio to the limit "
          f"{worst_all:.3e}{'  [mutate=' + MUT + ']' if MUT else ''}")


# ---------------------------------------------------------------- part gap
def gap_metrics(P, Rs, s, badcols=False):
    """(separation, gap in real parts, column ratio) of the selected pair of Js's eigenvalues, relative to max |lambda|"""
    lam, vec = np.linalg.eig(P._Jvac(Rs, s)[2:, 2:])
    o = np.argsort(lam.real); idx, oth = o[:2], o[2:]
    size = np.abs(lam).max()
    sep = min(abs(lam[i] - lam[j]) for i in idx for j in oth) / size
    gre = (lam.real[oth].min() - lam.real[idx].max()) / size
    try:
        Pm = vec[:, idx] @ np.linalg.inv(vec)[idx, :]
        cr = col_ratio(PBAD if badcols else Pm)
    except np.linalg.LinAlgError:
        cr = 0.0
    return sep, gre, cr


def col_ratio(Pm):
    sv = np.linalg.svd(Pm[:, [0, 2]], compute_uv=False)
    return sv[-1] / sv[0]


def planted_bad_projector():
    """a rank-2 projector with P e0 = P e2 (e0 - e2 in its kernel), built from a basis"""
    V = np.eye(4, dtype=complex)
    V[:, 2] = np.array([1, 0, -1, 0]); V[:, 3] = np.array([0.3, 1, 0.2, -0.5j])
    return V[:, :2] @ np.linalg.inv(V)[:2, :]


PBAD = planted_bad_projector()


def dist_to_box(z, x0, x1, ylo, yhi):
    return np.hypot(max(x0 - z.real, 0, z.real - x1), max(ylo - z.imag, 0, z.imag - yhi))


def gap_box(P, Rs, x0, x1, y0, y1, sym, cache_pts, badcols=False):
    """metrics (sep, gap in real parts, column ratio) over the contour, the cache points and a 20 x 20 grid; and the
    minima of sep and of the real-part gap polished by a local search inside the box"""
    ylo, yhi = (-y1, y1) if sym else (y0, y1)
    c = [complex(x1, y0), complex(x1, y1), complex(x0, y1), complex(x0, y0), complex(x1, y0)]
    sets = {'contour': [a + (b - a) * t for a, b in zip(c[:-1], c[1:]) for t in np.linspace(0, 1, 401)],
            'cache': list(cache_pts), 'grid': [complex(x, y) for x in np.linspace(x0, x1, 20) for y in np.linspace(ylo, yhi, 20)]}
    met = {k: np.array([gap_metrics(P, Rs, s, badcols) for s in v]) for k, v in sets.items() if v}
    allp = [s for k, v in sets.items() if v for s in v]
    allm = np.vstack(list(met.values()))
    polished = {}
    for j in (0, 1):
        f = lambda p, j=j: gap_metrics(P, Rs, complex(np.clip(p[0], x0, x1), np.clip(p[1], ylo, yhi)), badcols)[j]
        best = (allm[:, j].min(), allp[int(allm[:, j].argmin())])
        for i in np.argsort(allm[:, j])[:3]:
            r = minimize(f, [allp[i].real, allp[i].imag], method='Nelder-Mead',
                         options=dict(xatol=1e-13, fatol=1e-16, maxiter=3000))
            if r.fun < best[0]:
                best = (r.fun, complex(np.clip(r.x[0], x0, x1), np.clip(r.x[1], ylo, yhi)))
        polished[j] = best
    return met, polished


def moved_box(b, x0, x1, y0, y1):
    n = Box.__new__(Box); n.__dict__.update(b.__dict__)
    n.x0, n.x1, n.y0, n.y1, n.sym = x0, x1, y0, y1, False
    return n


def part_gap():
    print("# gap: the selected pair of the vacuum scalar block's eigenvalues against the other two, at R_s, relative to "
          "max abs(eigenvalue). Columns: smallest over the dense contour, the cache's points and the 20 x 20 grid, then the "
          "minimum polished inside the box, the gap in real parts (order, polished) and the smallest column ratio (cols)")
    waves = sorted({b.w for b in BOXES})
    print("   branch points s = +-i(1 +- w): " + "; ".join(f"w = {w}: " + ", ".join(f"{v:+.4f}i" for v in
          (1 + w, 1 - w, -(1 - w), -(1 + w))) for w in waves))
    print("| boxes | star, l | Re s | Im s | R_s | gap contour | gap cache | gap grid | gap polished (at s) | order | cols | "
          "dist to branch points |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
    rows = {}
    for b in BOXES:
        rows.setdefault((b.key()[:-1], b.rsmult), []).append(b)
    worst = [np.inf, np.inf, np.inf]
    for (k, rsm), bl in rows.items():
        b = bl[0]
        cp = list(box_cache(b).keys()) if os.path.exists(b.cachefile()) else []
        if MUT == 'branch':                   # the box moved onto the branch point i(1 - w), same size
            c0 = 1j * (1 - b.w); hx, hy = (b.x1 - b.x0) / 2, (b.y1 - b.y0) / 2
            b, cp = moved_box(b, c0.real - hx, c0.real + hx, c0.imag - hy, c0.imag + hy), []
        P = DBox(b).P
        met, pol = gap_box(P, P.Rs_fixed, b.x0, b.x1, b.y0, b.y1, b.sym, cp, badcols=MUT == 'cols')
        ylo, yhi = (-b.y1, b.y1) if b.sym else (b.y0, b.y1)
        d = min(dist_to_box(1j * v, b.x0, b.x1, ylo, yhi) for v in (1 + b.w, 1 - b.w, -(1 - b.w), -(1 + b.w)))
        mn = lambda k, j: float(met[k][:, j].min()) if k in met else float('nan')
        (pm, ps), (po, _), cols = pol[0], pol[1], min(float(m[:, 2].min()) for m in met.values())
        worst = [min(worst[0], pm), min(worst[1], po), min(worst[2], cols)]
        print(f"| {', '.join(x.name for x in bl)} | {b.w}, {b.l} | [{b.x0:g}, {b.x1:g}] | "
              f"{'abs(Im s) <= ' + format(b.y1, 'g') if b.sym else '[' + format(b.y0, 'g') + ', ' + format(b.y1, 'g') + ']'} | "
              f"{P.Rs_fixed:.0f} | {mn('contour', 0):.3e} | {mn('cache', 0):.3e} | {mn('grid', 0):.3e} | "
              f"{pm:.3e} ({ps.real:+.5g}{ps.imag:+.5g}i) | {po:.3e} | {cols:.3e} | {d:.3f} |", flush=True)
    sfx = f"  [mutate={MUT}]" if MUT else ""
    check("gap", worst[0] >= THR, f"smallest relative separation over the 61 boxes {worst[0]:.3e} (limit {THR:g}){sfx}")
    check("order", worst[1] >= THR, f"smallest gap in real parts over the 61 boxes {worst[1]:.3e} (limit {THR:g}){sfx}")
    check("cols", worst[2] >= THR, f"smallest sigma_min / sigma_max of columns 0 and 2 over the 61 boxes {worst[2]:.3e} "
          f"(limit {THR:g}){sfx}")
    print("# planted controls, which must fail")
    base = BYNAME['w_p7603_l2']
    P = DBox(base).P
    okc = True
    for name, c0 in (("i(1 - w)", 1j * (1 - base.w)), ("i(1 + w)", 1j * (1 + base.w))):
        b = moved_box(base, -0.01, 0.01, c0.imag - 0.01, c0.imag + 0.01); b.rsmult = 40.0
        Rs = DBox.rs_fixed(P, b, b.w)
        met, pol = gap_box(P, Rs, b.x0, b.x1, b.y0, b.y1, False, [])
        (pm, ps), (po, pso) = pol[0], pol[1]
        fails_it = pm < THR and po < THR
        okc &= fails_it
        print(f"   box Re s in [{b.x0:g}, {b.x1:g}], Im s in [{b.y0:.4f}, {b.y1:.4f}] around {name}, R_s {Rs:.0f}: grid min "
              f"{met['grid'][:, 0].min():.3e}, contour min {met['contour'][:, 0].min():.3e}, polished gap {pm:.3e} at "
              f"{ps.real:+.6g}{ps.imag:+.6g}i, polished order {po:.3e}: {'fails as required' if fails_it else 'PASSES (wrong)'}")
    check("gapctl", okc, "the boxes around s = i(1 - w) and i(1 + w) fail gap and order")
    P_ok = np.array([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]], dtype=complex)
    cg, cb = col_ratio(P_ok), col_ratio(PBAD)
    check("colsctl", cg >= THR and cb < THR, f"planted projector with independent columns {cg:.3e} passes, with P e0 = P e2 "
          f"{cb:.3e} fails (limit {THR:g})")


# ---------------------------------------------------------------- part norm
class EigScale:
    """np.linalg.eig with each eigenvector multiplied by a random nonzero complex factor (log-uniform modulus in [0.1, 10],
    uniform phase), seeded by s"""

    def __init__(self, s):
        self.rng = np.random.default_rng(abs(int(round(1e6 * s.real))) * 1000003 + abs(int(round(1e6 * s.imag))) + 7)

    def __enter__(self):
        self.orig = np.linalg.eig

        def eig(a):
            lam, vec = self.orig(a)
            c = np.exp(self.rng.uniform(np.log(0.1), np.log(10), vec.shape[1]) + 2j * np.pi * self.rng.uniform(size=vec.shape[1]))
            return lam, vec * c[None, :]
        np.linalg.eig = eig
        return self

    def __exit__(self, *a):
        np.linalg.eig = self.orig


def scalar_plane_variant(P, pf=None, noproj=False, nophase=False, nogauge=False):
    """polar_solver.scalar_plane, copied (Rs_fixed is set, so its s-dependent radius is unused). pf(s) multiplies the
    first start column (a planted pole); noproj starts from the selected eigenvectors themselves instead of the projector's
    columns; nophase drops the QR phase factor; nogauge fixes the first column's phase by its largest component"""
    def scalar_plane(s):
        Rs = P.Rs_fixed
        lam, vec = np.linalg.eig(P._Jvac(Rs, s)[2:, 2:])
        idx = np.argsort(lam.real)[:2]
        Vinv = np.linalg.inv(vec)
        Pm = vec[:, idx] @ Vinv[idx, :]
        if not np.iscomplexobj(s):
            Pm = Pm.real
        if noproj:
            Y = np.zeros((6, 2), dtype=vec.dtype); Y[2:, 0] = vec[:, idx[0]]; Y[2:, 1] = vec[:, idx[1]]
        else:
            Y = np.zeros((6, 2), dtype=Pm.dtype); Y[2:, 0] = Pm[:, 0]; Y[2:, 1] = Pm[:, 2]
        if pf is not None:
            Y = Y.astype(complex); Y[2:, 0] = Y[2:, 0] * pf(s)
        jf = lambda r, s_: P._Jvac(r, s_) if r > P.st.R else P.J(r, s_)
        Q, ph = P._integrate_plane(Y, Rs, P.rm, s, Jfun=jf)
        S = Q[2:, :]
        if nogauge:
            S = S.copy(); k = np.argmax(abs(S[:, 0])); S[:, 0] *= np.conj(S[k, 0] / abs(S[k, 0])); ph = 1.0
        if nophase:
            ph = 1.0
        return S, ph
    return scalar_plane


def decaying_plane_variant(P, pf=None, nophase=False, nogauge=False):
    """polar_solver.decaying_plane, copied; pf(s) multiplies the Zerilli column zv (a planted pole)"""
    def decaying_plane(s):
        q = P.zerilli_q(s)
        f = 1 - 2 * P.Mext / P.rm
        zv = np.array([1.0, q / f])
        if pf is not None:
            zv = zv * pf(s)
        m = P.moncrief(P.rm) @ zv
        S, ph = P.scalar_plane(s)
        Z = np.zeros((6, 3), dtype=np.result_type(m, S))
        Z[0:2, 0] = m
        Z[2:, 1:] = S
        Qz, p = P._qr(Z)
        if nogauge:
            Qz = Qz.copy(); k = np.argmax(abs(Qz[:, 0])); Qz[:, 0] *= np.conj(Qz[k, 0] / abs(Qz[k, 0])); p = 1.0
        if nophase:
            p = 1.0
        return Qz, ph * p
    return decaying_plane


def norm_points(b):
    p = harm_points(b)
    return p + [complex(b.x0 + (b.x1 - b.x0) * fx, b.y0 + (b.y1 - b.y0) * fy) for fx, fy in NORM_EXTRA]


def record_sigma(P, out):
    """append sigma_min of the scalar plane's 4 x 2 block S at every call"""
    orig = P.scalar_plane

    def scalar_plane(s):
        r = orig(s)
        out.append(float(np.linalg.svd(r[0], compute_uv=False)[-1]))
        return r
    P.scalar_plane = scalar_plane


def part_norm():
    print("# norm: D with each eigenvector of np.linalg.eig in scalar_plane scaled by a random nonzero complex factor")
    worst, wsig = 0.0, np.inf
    for b in jobs(HARM6):
        st = Store(f"norm_{b.name}"); db, dp = DBox(b), DBox(b)
        if MUT == 'noproj':
            db.P.scalar_plane = scalar_plane_variant(db.P, noproj=True)
            dp.P.scalar_plane = scalar_plane_variant(dp.P, noproj=True)
        rb, rp = [], []
        record_sigma(db.P, rb); record_sigma(dp.P, rp)

        def base(s):
            return db.D(s), rb[-1]

        def scaled(s):
            with EigScale(s):
                return dp.D(s), rp[-1]
        wb = 0.0
        try:
            for s in norm_points(b):
                d0, g0 = st.get(('D', MUT, s), lambda s=s: base(s))
                d1, g1 = st.get(('Dscaled', MUT, s), lambda s=s: scaled(s))
                r = abs(d1 - d0) / abs(d0); wb = max(wb, r); wsig = min(wsig, g0, g1)
                print(f"    {b.name} s = {s.real:+.6g}{s.imag:+.6g}i: |D| {abs(d0):.3e}, relative difference {r:.2e}, "
                      f"sigma_min(S) {min(g0, g1):.4f}", flush=True)
        except Stop as e:
            print(f"  STOP {e}"); check("norm", False, f"stopped in {b.name}: {e}"); return
        worst = max(worst, wb)
        print(f"  {b.name}: worst {wb:.2e} ({st.new} new evaluations)", flush=True)
    check("norm", worst <= 1e-10, f"largest relative difference of D {worst:.2e} (limit 1e-10) at 5 points in "
          f"{len(jobs(HARM6))} boxes{'  [mutate=' + MUT + ']' if MUT else ''}")
    print(f"   smallest sigma_min(S) of the scalar plane's block at these points (reported, no criterion): {wsig:.4f}")


# ---------------------------------------------------------------- part pole
def part_pole():
    b = BYNAME['w_n1_l2']
    s0 = 0.05 + 0.1j
    c = lambda s: 1.0 / ((s - s0) * (s - np.conj(s0)))
    cache = box_cache(b)
    Dbase = lambda s: cache[ckey(s)]
    sfx = f"  [mutate={MUT}]" if MUT else ""
    # the copied counter on the cached D, against #126's log
    log = open(os.path.join(I126, 'w_n1_l2.log'), encoding='utf-8').read()
    m = re.search(r"winding = ([+-]\d+\.\d+)  \(max phase step ([\d.]+) rad, min\|D\| (\S+)\)", log)
    nev = int(re.findall(r"\((\d+) evals", log)[-1])
    seen = set()

    def Fc(s):
        seen.add(ckey(s)); return Dbase(s)
    cor = b.corners()[::-1] if MUT == 'reverse' else b.corners()
    N0, ms0, nm0 = count_contour(Fc, cor, True)
    Np, msp, _ = count_contour(lambda s: Dbase(s) * (s - s0) * (s - np.conj(s0)), cor, True)
    check("counter", abs(N0 - float(m.group(1))) < 0.05 and f"{ms0:.2f}" == m.group(2) and f"{nm0:.2e}" == m.group(3)
          and len(seen) == nev and abs(Np - N0 - 2) < 0.05,
          f"copied counter on the cached D: winding {N0:+.4f} (#126 {m.group(1)}), max step {ms0:.2f} ({m.group(2)}), min|D| "
          f"{nm0:.2e} ({m.group(3)}), {len(seen)} points ({nev}), the planted pair on the cache gives {Np - N0:+.4f} of +2{sfx}")
    stc = Store(f"pole_{b.name}_clone")
    st = stc
    kw = dict(nophase=MUT == 'nophase', nogauge=MUT == 'nogauge')
    stage = 'clone'
    try:
        dbc = DBox(b); dbc.P.scalar_plane = scalar_plane_variant(dbc.P); dbc.P.decaying_plane = decaying_plane_variant(dbc.P)
        a0, a1 = cor[0], cor[1]
        dev = 0.0
        for t in (0.3, 0.5, 0.7):
            s = ckey(a0 + (a1 - a0) * t)
            v = stc.get(('clone', s), lambda s=s: dbc.D(s))
            dev = max(dev, abs(v - Dbase(s)) / abs(Dbase(s)))
        check("clone", dev <= 1e-12, f"the copied scalar_plane and decaying_plane reproduce the cached D at 3 points to {dev:.1e} "
              "(limit 1e-12)")
        for stage in PLANTS:
            st = Store(f"pole_{b.name}_{stage}" + (f"_{MUT}" if MUT else ""))
            d = DBox(b)
            if stage == 'z':
                d.P.decaying_plane = decaying_plane_variant(d.P, pf=c, **kw)
            else:
                d.P.scalar_plane = scalar_plane_variant(d.P, pf=c, **kw)
            Fp = lambda s, d=d, tag=stage: st.get((tag, MUT, ckey(s)), lambda: d.D(ckey(s)))
            t0 = time.time(); n0 = st.new
            print(f"# plant in the {'Zerilli column' if stage == 'z' else 'scalar plane first column'}: s0 = {s0}, factor "
                  f"1/((s - s0)(s - conj s0)), guard all{sfx}", flush=True)
            Nn, msn, nmn = count_contour(Fp, cor, True)
            pts = [k[2] for k in st.d if k[0] == stage and k[1] == MUT]
            idn = [abs(st.d[(stage, MUT, s)] / Dbase(s) - c(s) / abs(c(s))) for s in pts if s in cache]
            wid = max(idn) if idn else float('nan')
            print(f"   count {Nn:+.4f}, max step {msn:.2f} rad, min |D| {nmn:.2e}, {len(pts)} points ({st.new - n0} new "
                  f"evaluations, {time.time() - t0:.0f} s); identity D_pole / D = c / |c| at {len(idn)} cached points: largest "
                  f"deviation {wid:.2e}", flush=True)
            check(f"pole_{stage}", abs(Nn - N0 + 2) < 0.05 and msn < 0.45 and wid <= 1e-8,
                  f"winding {Nn:+.4f}, change {Nn - N0:+.4f} of -2, max step {msn:.2f} (limit 0.45), identity deviation "
                  f"{wid:.2e} (limit 1e-8){sfx}")
    except Stop as e:
        print(f"  STOP {e}")
        if stage == 'clone':
            check("clone", False, f"no evaluation: {e}")
        else:
            check(f"pole_{stage}", False, f"no count (the identity was checked at the points reached): {e}{sfx}")
            pts = [k[2] for k in st.d if k[0] == stage and k[1] == MUT and k[2] in cache]
            if pts:
                dv = max(abs(st.d[(stage, MUT, s)] / Dbase(s) - c(s) / abs(c(s))) for s in pts)
                print(f"   identity D_pole / D = c / |c| at the {len(pts)} cached points reached: largest deviation {dv:.2e}")


# ---------------------------------------------------------------- part axis (zerilli_zeros.py --axis, copied)
def make_qaxis(P, b):
    rm, M, l = P.rm, P.Mext, b.l
    lam = (l - 1) * (l + 2) / 2.0
    cor = b.corners()
    smin = min(abs(a + (c - a) * t) for a, c in zip(cor[:-1], cor[1:]) for t in np.linspace(0, 1, 401))
    Rz = rm + 1.0 * (40.0 + 20.0 * l) / smin + 50.0

    def V(r):
        f = 1 - 2 * M / r
        return 2 * f * (lam**2 * r**2 * ((lam + 1) * r + 3 * M) + 9 * M**2 * (lam * r + M)) / (r**3 * (lam * r + 3 * M)**2)
    NSER = 400
    Nc = np.zeros(NSER + 1); Nc[:4] = [lam**2 * (lam + 1), 3 * M * lam**2, 9 * M**2 * lam, 9 * M**3]
    Wc = 2 * np.convolve(Nc, [(j + 1) * (-3 * M / lam)**j / lam**2 for j in range(NSER + 1)])[:NSER + 1]

    def jost_start(s, R):
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

    def qaxis(s):
        T = Rz - rm
        R0 = rm - 1j * T
        ya = jost_start(s, R0)

        def A(sig):
            p = np.exp(sig) - 1.0; r = rm - 1j * p; f = 1 - 2 * M / r
            a = (-1j * (1.0 + p) / f) * np.array([[s, 1.0], [s * s + V(r), s]])
            return np.block([[a.real, -a.imag], [a.imag, a.real]])
        sol = solve_ivp(lambda sg, y: A(sg) @ y, (np.log(1.0 + T), 0.0), np.r_[ya.real, ya.imag], method='Radau',
                        jac=lambda sg, y: A(sg), rtol=1e-11, atol=1e-14)
        assert sol.success, f"Radau failed at s = {s}: {sol.message}"
        z, p = sol.y[0, -1] + 1j * sol.y[2, -1], sol.y[1, -1] + 1j * sol.y[3, -1]
        return p / z
    qaxis.Rz = Rz
    return qaxis


def chain_of(cache_keys, corners):
    """cached points on each edge sorted along it (as diagnose_plant.py); returns the list of edge chains and the number of
    points on no edge"""
    edges, used = [], set()
    for a, c in zip(corners[:-1], corners[1:]):
        pts = []
        for s in cache_keys:
            t = (s - a) / (c - a)
            if abs(t.imag) < 1e-8 and -1e-8 <= t.real <= 1 + 1e-8:
                pts.append((t.real, s)); used.add(s)
        pts.sort()
        edges.append([s for _, s in pts])
    return edges, len(set(cache_keys) - used)


def chain_wind(F, edges, rev=False):
    """winding of F along the chain (no doubling: upper boxes) and the largest increment"""
    tot, mx = 0.0, 0.0
    for e in (edges[::-1] if rev else edges):
        e = e[::-1] if rev else e
        v = [F(s) for s in e]
        inc = [np.angle(v[i + 1] / v[i]) for i in range(len(v) - 1)]
        tot += sum(inc); mx = max([mx] + [abs(x) for x in inc])
    return tot / (2 * np.pi), mx


def part_axis():
    sfx = f"  [mutate={MUT}]" if MUT else ""
    s0 = -0.0005 + 0.12j
    for b in jobs(['w_n3_l2', 'w_n3_l3']):
        name = b.name
        print(f"# {name}: Re s in [{b.x0:g}, {b.x1:g}], Im s in [{b.y0:g}, {b.y1:g}], l = {b.l}")
        db, da = DBox(b), DBox(b)
        qa = make_qaxis(da.P, b)
        st = Store(f"axis_{name}"); qtag = MUT if MUT == 'wrongq' else None
        qf = lambda s: st.get(('q', s), lambda: qa(s)) * (1 + 1e-4 if MUT == 'wrongq' else 1.0)
        da.P.zerilli_q = lambda s, r_start=None: qf(s)
        cache = box_cache(b)
        log = open(os.path.join(I123, f"{name[2:]}.log"), encoding='utf-8').read()
        ql = {}
        for ln in log.splitlines():
            mt = re.match(r"\s+q at the (contour|holo) points: (.*)$", ln)
            if mt:
                for item in mt.group(2).split(', '):
                    ss, qq = item.split(' ')
                    ql[complex(ss.strip('()'))] = complex(qq)
        try:
            dq = max(abs(qf(s) - q) / abs(q) for s, q in ql.items())
            check("axis_q", len(ql) == 6 and dq <= 1e-10, f"q from the copied --axis function against the I123 log's {len(ql)} "
                  f"values (R_z {qa.Rz:.1f}): largest relative difference {dq:.2e} (limit 1e-10){sfx}")
            keys = sorted(cache, key=lambda s: (s.real, s.imag))
            Da = lambda s: st.get(('Da', qtag, s), lambda: da.D(s))
            t0 = time.time()
            for s in keys:
                Da(s)
        except Stop as e:
            print(f"  STOP {e}"); check("axis_q", False, f"stopped: {e}"); return
        print(f"   D_axis at the {len(keys)} cached points: {st.new} new evaluations, {time.time() - t0:.0f} s", flush=True)
        edges, off = chain_of(keys, b.corners())
        names = ('right Re s = x1', 'top Im s = y1', 'left Re s = x0', 'bottom Im s = y0')
        dmin = min(abs(cache[s]) for s in keys)
        print(f"   cached points on the contour: {sum(len(e) for e in edges)} on edges (corners twice), {off} on none; "
              f"min |D| on the contour {dmin:.3e}")
        for nm, e in zip(names, edges):
            rel = [abs(Da(s) - cache[s]) / abs(cache[s]) for s in e]
            ab = [abs(Da(s) - cache[s]) for s in e]
            print(f"   {nm}: {len(e)} points, largest |D_axis - D| / |D| {max(rel):.3e}, largest |D_axis - D| {max(ab):.3e} "
                  f"({max(ab) / dmin:.3e} of min |D|)")
        rev = MUT == 'reverse'
        w_ref, m_ref = chain_wind(lambda s: cache[s], edges, rev)
        w_a, m_a = chain_wind(Da, edges, rev)
        w_pa, m_pa = chain_wind(lambda s: Da(s) * (s - s0), edges, rev)
        w_pr, m_pr = chain_wind(lambda s: cache[s] * (s - s0), edges, rev)
        check("axis_ref", abs(w_ref) < 0.05 and off == 0 and m_ref < 1.0,
              f"dense-chain winding of the cached D {w_ref:+.4f} (#126: 0), planted {w_pr - w_ref:+.4f} of +1, largest increment "
              f"{max(m_ref, m_pr):.3f} rad (limit 1){sfx}")
        check("axis_wind", abs(w_a) < 0.05 and m_a < 1.0, f"dense-chain winding of D_axis {w_a:+.4f} (expected 0), largest "
              f"increment {m_a:.3f} rad (limit 1){sfx}")
        check("axis_plant", abs(w_pa - w_a - 1) < 0.05 and m_pa < 1.0, f"with the planted zero s0 = {s0}: winding {w_pa:+.4f}, "
              f"change {w_pa - w_a:+.4f} of +1, largest increment {m_pa:.3f} rad (limit 1){sfx}")
        # harmonic at the three points of harm: D (zerilli_q) and D_axis
        L = box_L(b); stz = Store(f"axisharm_{name}")
        for lab, F in (("D", lambda s: stz.get(('D', s), lambda: db.D(s))),
                       ("D_axis", lambda s: stz.get(('Da', qtag, s), lambda: da.D(s)))):
            try:
                ok, worst = harm_run(b, F, lab)
            except Stop as e:
                print(f"  STOP {e}"); check("axisharm", False, f"stopped: {e}"); return
            msg = f"{lab} in {name}: worst |lap phi| L^2 / limit {worst:.3e}"
            if lab == "D_axis":
                check("axisharm", ok, msg + sfx)
            else:
                print(f"  [{'PASS' if ok else 'FAIL'}] harm_n3 (reported, the check of harm): {msg}", flush=True)


# ---------------------------------------------------------------- run
if __name__ == '__main__':
    fn = {'factors': part_factors, 'gap': part_gap, 'norm': part_norm, 'pole': part_pole, 'axis': part_axis}
    print(f"d_factors.py parts {','.join(PARTS)}, RTOL {RTOL:g}, raw {RAW}" + (f", mutate={MUT}" if MUT else ""), flush=True)
    for p in PARTS:
        print(f"\n===== part {p}", flush=True)
        try:
            fn[p]()
        except Stop as e:
            print(f"STOP {e}"); STOPPED.append(p); check(p, False, f"stopped: {e}")
    print(f"\n{'FAILED: ' + ', '.join(fails) if fails else 'all checks passed'}  ({time.time() - T_START:.0f} s)")
    sys.exit(3 if STOPPED else (1 if fails else 0))

"""Order 4 of the matched expansion (matched.py, order3.py): Omega2 from G4 = 0, and the checks of issue #26.

W^2 a^2 = s + Omega, s = 2(D - 3), eps = 1/s, Omega = Omega0 + Omega1 eps + Omega2 eps^2 + ..., R = r^s, l = ln R.
The order equations (E0 derives them from the master equation)
  dY_k/dl = F_{k-1},  F_{k-1} = sum_{j<k} L_{k-1-j} l^j / j!,   (R - 1) L_k' - L_k/2 = h_k,
  h_k = (1 - 1/R) L_{k-1} + sum_{j<k} Y_j c_{k-j} + Y_k/2,   L_k = -2 Omega_{k-1} + u C_k,   G_k = C_k(inf) = 0,
with u = sqrt(R - 1) = tan(th), C_k = int_1^R (h_k - Omega_{k-1}) (R' - 1)^{-3/2} dR' = 2 int_0^th (h_k - Omega_{k-1}) / sin^2.
G_k is affine in Omega_{k-2} with slope -4 pi: G2 = 0 gives Omega0, G3 = 0 Omega1, G4 = 0 Omega2.

Method A: closed forms in th (K = 2 th ln 2 - Cl2(pi - 2 th), T = int_0^th p ln cos p dp with Cl2, Cl3; Y1..Y3, L0..L2) and
two integrations by parts (E1, E2 of the raw README) that leave G_k as one integral of closed forms.
Method B: the hierarchy from its definitions on Chebyshev panels in v = sqrt(l); no closed form, no exchange.
Reduction by hand (part `reduce`, H0 .. H8): at Omega0 = -4 ln 2, with ell = ln cos, u = tan, C = Cl2(pi - 2 th),
D3 = Cl3(pi - 2 th) + 3 zeta(3)/4, the integrand of G_k is A(ell)/S + B(ell) + gamma C + delta D3, S = sin^2; the C and D3 terms
integrate to 64 (Omega0 - 1) int tan C = 64 (Omega0 - 1) pi^3/12, the rest to the moments int ell^n, so that G3 = 0 gives
Omega1 = 2 pi^2/3 + 8 ln^2 2 and G4 = 0 gives Omega2 = 80 zeta(3) - (8/3) pi^2 ln 2 - (32/3) ln^3 2.
L3 = -2 h3 + 2 u B(th), B closed in ell, u, th, C and Lc2 = int_0^th ln^2 cos (S4, N5).

    python3 order4.py [sym] [closed] [lc2] [reduce] [n4b] [n4bg] [evals] [remainder] [--mutate=NAME] [--w2=w2_output.txt] [--w2raw=FILE]
(the parts run in this order; evals and remainder together give R1, R2; the nested quadratures n4b, n4bg take 30 to 70 min)

No argument runs everything. The exit status is 1 when a check fails. Status words: sympy residual 0 is "checked symbolically",
an agreement of integrals "numerical evidence", a PSLQ relation an "identification".
"""
import re
import sys
import time

import mpmath as mp
import sympy as sp

MUTATIONS = ("wrong_h", "wrong_Om0", "wrong_Om1", "T_const", "K_sign", "bad_exchange", "coarse_B", "short_tail",
             "pslq_wrong_basis", "wrong_Om2_series", "wrong_J3", "wrong_tanC", "wrong_Y3", "wrong_B", "wrong_P", "wrong_L3", "wrong_Lc2")
MUTATE = None
W2FILE = "w2_output.txt"
W2RAW = None
parts = []
for a in sys.argv[1:]:
    if a.startswith("--mutate="):
        MUTATE = a.split("=", 1)[1]
        assert MUTATE in MUTATIONS, f"unknown mutation {MUTATE}; known: {MUTATIONS}"
    elif a.startswith("--w2="):
        W2FILE = a.split("=", 1)[1]
    elif a.startswith("--w2raw="):
        W2RAW = a.split("=", 1)[1]
    else:
        parts.append(a)
PARTS = parts or ["sym", "closed", "lc2", "reduce", "n4b", "n4bg", "evals", "remainder"]

RESULTS = []


def report(cid, ok, text):
    RESULTS.append((cid, ok))
    print(f"{'PASS' if ok else 'FAIL'}    {cid}: {text}", flush=True)


def digits(a, b):
    """-log10 of the relative difference, capped at the working precision in digits"""
    a, b = mp.mpf(a), mp.mpf(b)
    d = abs(a - b)
    return mp.mpf(mp.mp.dps) if d == 0 else min(mp.mpf(mp.mp.dps), -mp.log10(d / abs(a)))


def nd(x, n=40):
    return mp.nstr(x, n)


# ---------------------------------------------------------------------------------------------------------------
# the table of c_m and h_m (integer arithmetic only, so that it also runs on sympy symbols)
# ---------------------------------------------------------------------------------------------------------------
def cc(m, Om, l, R):
    if m == 1:
        return Om[0] + 3 * l / 2 - 1 + 1 / R
    if m == 2:
        return Om[0] * l + Om[1] + l**2 / 4 + l + (l / R if MUTATE == "wrong_h" else -l / R)
    if m == 3:
        return Om[0] * l**2 / 2 + Om[1] * l + Om[2] + l**3 / 4 - l**2 / 2 + l**2 / (2 * R)
    if m == 4:
        return Om[0] * l**3 / 6 + Om[1] * l**2 / 2 + Om[2] * l + Om[3] + l**4 / 48 + l**3 / 6 - l**3 / (6 * R)
    raise ValueError(m)


def hk(m, Om, l, R, Y, L, om_last=True):
    """h_m; om_last False leaves out the constant Omega_{m-1} (it is h_m(1))"""
    if m == 0:
        return Y[0] / 2
    o = list(Om)
    if not om_last:
        o[m - 1] = 0
    return (1 - 1 / R) * L[m - 1] + sum(Y[j] * cc(m - j, o, l, R) for j in range(m)) + Y[m] / 2


def Fsum(k, L, l):
    return sum(L[k - 1 - j] * l**j / sp.factorial(j) if isinstance(l, sp.Basic) else L[k - 1 - j] * l**j / mp.factorial(j)
               for j in range(k))


# ---------------------------------------------------------------------------------------------------------------
# E0, S1, S2, S3: sympy
# ---------------------------------------------------------------------------------------------------------------
def check_sym():
    t0 = time.time()
    # E0: the exact system's epsilon-expansion against the table
    R, eps, l = sp.symbols("R epsilon l", positive=True)
    Om = sp.symbols("Omega0:4")
    Ys, Ls = sp.symbols("Y0:5"), sp.symbols("L0:5")
    r = sp.exp(eps * l)
    Y = sum(Ys[k] * eps**k for k in range(5))
    L = sum(Ls[k] * eps**k for k in range(5))
    O = sum(Om[k] * eps**k for k in range(4))
    half = sp.Rational(1, 2)
    rhs = ((1 + eps * O) * r * Y + L / (2 * R) + (half + eps) * L * (R - 1) / R - (half + eps) * Y * (R - 1) / (R * r)
           - Y / (2 * r * R))
    ser = sp.expand(sp.series(rhs, eps, 0, 5).removeO())
    yser = sp.expand(sp.series(eps * L * r, eps, 0, 5).removeO())
    res = []
    for k in range(5):
        res.append(sp.simplify(ser.coeff(eps, k) - Ls[k] / 2 - hk(k, Om, l, R, Ys, Ls)))
        if k >= 1:
            res.append(sp.simplify(yser.coeff(eps, k) - Fsum(k, Ls, l)))
    report("E0", all(x == 0 for x in res),
           f"h_k (k = 0 .. 4) and dY_k/dl (k = 1 .. 4) from the epsilon-expansion of the master equation against the table: "
           f"residuals {[str(x) for x in res if x != 0] or '0'} ({time.time() - t0:.1f} s)")

    # S1, S2, S3 in u = tan(th): generators th, l, K, T with the derivation D = d/du
    u, th, lg, K, T, B = sp.symbols("u th l K T B")
    O0, O1, O2 = sp.symbols("Omega0 Omega1 Omega2")
    Rr = 1 + u**2

    def D(f):
        return (sp.diff(f, u) + sp.diff(f, th) / Rr + sp.diff(f, lg) * 2 * u / Rr + sp.diff(f, K) * lg / Rr
                - sp.diff(f, T) * th * lg / (2 * Rr) + sp.diff(f, B) * 0)

    Om3 = [O0, O1, O2, 0]
    Yc = [1, -lg, -2 * O0 * lg - sp.Rational(3, 2) * lg**2,
          -2 * O1 * lg + lg**3 / 2 + 8 * lg**2 + 8 * O0 * lg + 8 * O0 * th**2 - 16 * O0 * u * th - 32 * u * K + (16 if MUTATE == "wrong_Y3" else 32) * K * th
          + 64 * T]
    Lc = [-1, -2 * O0 - 2 * lg, -2 * O1 + 2 * O0 * lg + 4 * lg**2 - 8 * O0 * u * th - 16 * u * K]
    dldu = 2 * u / Rr
    s1 = [sp.simplify(sp.together(D(Yc[k]) / dldu - Fsum(k, Lc, lg))) for k in range(1, 4)]
    report("S1", all(x == 0 for x in s1), f"Y_1, Y_2, Y_3 closed forms: dY_k/dl - F_(k-1) = {[str(x) for x in s1]}")
    s2 = []
    for k in range(3):
        h = hk(k, Om3, lg, Rr, Yc, Lc)
        s2.append(sp.simplify(sp.together((u * D(Lc[k]) - Lc[k]) / 2 - h)))
    report("S2", all(x == 0 for x in s2), f"L_0, L_1, L_2 closed forms: (R-1)L_k' - L_k/2 - h_k = {[str(x) for x in s2]}")
    # S3: L_3 = -2 h_3 + 2 u B with dB/du = h_3,u / u
    h3 = hk(3, Om3, lg, Rr, Yc, Lc)
    L3 = -2 * h3 + 2 * u * B
    DB = D(h3) / (2 * u if MUTATE == "wrong_B" else u)

    def D2(f):
        return D(f) + sp.diff(f, B) * DB

    s3 = sp.simplify(sp.together((u * D2(L3) - L3) / 2 - h3))
    report("S3", s3 == 0, f"L_3 = -2 h_3 + 2 u B, B' = h_3,u/u: (R-1)L_3' - L_3/2 - h_3 = {s3}")
    # S4: the closed form of B in ell, u, th, C = Cl2(pi - 2 th), Lc2 = int_0^th ln^2 cos (derivation in th, d/dR = D/(2u(1 + u^2)))
    q, Cs, D3s, Lc = sp.symbols("q C D3 Lc2")
    w0q = -4 * q
    Rr2 = 1 + u**2
    S_, sc, c2 = u**2 / Rr2, u / Rr2, 1 / Rr2
    Om_q = [w0q, O1, O2, 0]
    L2q = -2 * O1 - 4 * w0q * lg + 16 * lg**2 + 16 * u * Cs
    Y3q = 4 * O1 * lg - 4 * lg**3 + 32 * lg**2 - 16 * w0q * lg + 32 * u * Cs - 16 * D3s
    Lq = [-1, -2 * w0q + 4 * lg, L2q]
    Yq = [1, 2 * lg, 4 * w0q * lg - 6 * lg**2, Y3q]
    h3q = hk(3, Om_q, -2 * lg, Rr2, Yq, Lq)

    def Dth(f):
        return (sp.diff(f, u) * Rr2 - sp.diff(f, lg) * u + sp.diff(f, Cs) * (2 * lg + 2 * q) + sp.diff(f, D3s) * 2 * Cs
                + sp.diff(f, th) + sp.diff(f, Lc) * lg**2)

    kc = 16 if MUTATE == "wrong_L3" else 32
    Bc = -32 * (lg + q) * Cs + 16 * S_ * Cs + kc * Lc + 16 * sc * lg * (lg + 2 * q) - 4 * O1 * th - 2 * O1 * sc
    L3c = -2 * h3q + 2 * u * Bc
    s4 = sp.simplify(sp.together(u * Dth(L3c) / (2 * Rr2) - L3c / 2 - h3q))
    report("S4", s4 == 0, f"L_3 = -2 h_3 + 2 u B with B = -32(ell + ln 2)C + 16 S C + 32 Lc2 + 16 sin cos ell (ell + 2 ln 2) - 4 Omega1 th - 2 Omega1 sin cos "
           f"(Lc2 = int_0^th ln^2 cos, C = Cl2(pi - 2 th)): (R-1)L_3' - L_3/2 - h_3 = {s4}")


# ---------------------------------------------------------------------------------------------------------------
# closed forms in th (numerical)
# ---------------------------------------------------------------------------------------------------------------
def cl_args(th, x, tail):
    """(a, c, s2, u): a = pi - 2 th; in the tail cos th = sin x, sin^2 th = cos^2 x, tan th = cot x, a = 2x"""
    if tail:
        return 2 * x, mp.sin(x), mp.cos(x) ** 2, mp.cot(x)
    return mp.pi - 2 * th, mp.cos(th), mp.sin(th) ** 2, mp.tan(th)


def K_of(th, a):
    return 2 * th * mp.ln2 + (mp.clsin(2, a) if MUTATE == "K_sign" else -mp.clsin(2, a))


def T_of(th, a):
    const = 0 if MUTATE == "T_const" else -mp.mpf(3) / 16 * mp.zeta(3)
    return th / 2 * mp.clsin(2, a) - mp.clcos(3, a) / 4 + const - mp.ln2 * th**2 / 2


def closed(th, x, tail, Om):
    a, c, s2, u = cl_args(th, x, tail)
    lc = mp.log(c)
    l = -2 * lc
    R = 1 / c**2
    K, T = K_of(th, a), T_of(th, a)
    Y = [mp.mpf(1), -l, -2 * Om[0] * l - 3 * l**2 / 2,
         -2 * Om[1] * l + l**3 / 2 + 8 * l**2 + 8 * Om[0] * l + 8 * Om[0] * th**2 - 16 * Om[0] * u * th - 32 * u * K
         + 32 * K * th + 64 * T]
    L = [mp.mpf(-1), -2 * Om[0] - 2 * l,
         -2 * Om[1] + 2 * Om[0] * l + 4 * l**2 - 8 * Om[0] * u * th - 16 * u * K]
    return c, lc, l, u, R, s2, Y, L, K, T


def theta_pt(th):
    th = mp.mpf(th)
    return th, mp.pi / 2 - th, False


# ---------------------------------------------------------------------------------------------------------------
# Method A
# ---------------------------------------------------------------------------------------------------------------
GUARD = 30


def A_integrand(k, Om):
    """integrand of (E1) with (E2) inserted: G_k = int_0^{pi/2} of it; Om[k-1] is dropped (it is h_k(1))"""
    fac = 4 if MUTATE == "bad_exchange" else 8

    def f(th, x, tail):
        c, lc, l, u, R, s2, Y, L, K, T = closed(th, x, tail, Om)
        h_prev = s2 * L[k - 2] + sum(Y[j] * cc(k - 1 - j, Om, l, R) for j in range(k - 1)) + Y[k - 1] / 2
        t1 = fac * h_prev * lc / s2
        t2 = 2 * sum(L[k - 1 - j] * l**j / mp.factorial(j) for j in range(1, k))
        o = list(Om)
        o[k - 1] = 0
        t3 = 2 * (sum(Y[j] * cc(k - j, o, l, R) for j in range(k))) / s2
        return t1 + t2 + t3

    return f


def A_G(k, Om):
    f = A_integrand(k, Om)
    left = mp.quad(lambda t: f(t, mp.pi / 2 - t, False), [0, mp.mpf(1) / 2, 1], method="gauss-legendre")
    right = mp.quad(lambda x: f(mp.pi / 2 - x, x, True), [0, mp.pi / 2 - 1])
    return left + right


def solve_omega(Gfun, k, Om):
    """Omega_{k-2} from the affine G_k; returns (value, slope)"""
    o0 = list(Om)
    o0[k - 2] = mp.mpf(0)
    o1 = list(Om)
    o1[k - 2] = mp.mpf(1)
    g0, g1 = Gfun(k, o0), Gfun(k, o1)
    return -g0 / (g1 - g0), g1 - g0


def run_A(dps):
    t0 = time.time()
    out = {}
    with mp.workdps(dps + GUARD):
        Om = [mp.mpf(0)] * 4
        for k in (2, 3, 4):
            Om_in = list(Om)
            if MUTATE == "wrong_Om0" and k >= 3:
                Om_in[0] = -4 * mp.ln2 + mp.mpf(1) / 100
            if MUTATE == "wrong_Om1" and k >= 4:
                Om_in[1] = 2 * mp.pi**2 / 3 + 8 * mp.ln2**2 + mp.mpf(1) / 100
            val, slope = solve_omega(A_G, k, Om_in)
            Om[k - 2] = val
            out[f"Om{k - 2}"], out[f"slope{k}"] = +val, +slope
        out["secs"] = time.time() - t0
    return out


# ---------------------------------------------------------------------------------------------------------------
# Method B
# ---------------------------------------------------------------------------------------------------------------
class PanelB:
    def __init__(self, N, w, V):
        self.N, self.w = N, mp.mpf(w)
        self.P = int(mp.ceil(V / self.w))
        N = self.N
        xs = [mp.cos(mp.pi * (2 * i + 1) / (2 * N)) for i in range(N)][::-1]
        T = []                                           # T[i][n] = T_n(x_i), n <= N
        for i in range(N):
            t = [mp.mpf(1), xs[i]]
            for n in range(2, N + 1):
                t.append(2 * xs[i] * t[-1] - t[-2])
            T.append(t)
        Vinv = mp.inverse(mp.matrix([[T[i][n] for n in range(N)] for i in range(N)]))

        def IT(n, ti, tm, tp):                           # int_{-1}^x T_n from T_{n-1}, T_{n+1} at x
            if n == 0:
                return ti + 1
            if n == 1:
                return (ti * ti - 1) / 2
            return (tp / (n + 1) - tm / (n - 1)) / 2 - (mp.mpf(-1) ** (n + 1) / (n + 1) - mp.mpf(-1) ** (n - 1) / (n - 1)) / 2

        Q = mp.matrix(N, N)
        qend = [IT(n, mp.mpf(1), mp.mpf(1), mp.mpf(1)) if n >= 2 else IT(n, mp.mpf(1), 0, 0) for n in range(N)]
        for n in range(N):
            for i in range(N):
                Q[i, n] = IT(n, xs[i] if n == 0 else T[i][n], T[i][n - 1] if n else 0, T[i][n + 1])
        S = Q * Vinv
        self.S = [[S[i, j] for j in range(N)] for i in range(N)]
        self.wend = [sum(qend[n] * Vinv[n, j] for n in range(N)) for j in range(N)]
        self.half = self.w / 2
        self.v = [[p * self.w + self.half * (x + 1) for x in xs] for p in range(self.P)]
        self.l = [[t * t for t in row] for row in self.v]
        self.R = [[mp.exp(t) for t in row] for row in self.l]
        self.Rm1 = [[mp.expm1(t) for t in row] for row in self.l]
        self.u = [[mp.sqrt(t) for t in row] for row in self.Rm1]

    def cum_left(self, F):
        N, half = self.N, self.half
        out, acc = [], mp.mpf(0)
        for p in range(self.P):
            out.append([acc + half * sum(self.S[i][j] * F[p][j] for j in range(N)) for i in range(N)])
            acc += half * sum(self.wend[j] * F[p][j] for j in range(N))
        return out

    def cum_right(self, F):
        N, half = self.N, self.half
        tot = [half * sum(self.wend[j] * F[p][j] for j in range(N)) for p in range(self.P)]
        out, acc = [None] * self.P, mp.mpf(0)
        for p in range(self.P - 1, -1, -1):
            out[p] = [acc + tot[p] - half * sum(self.S[i][j] * F[p][j] for j in range(N)) for i in range(N)]
            acc += tot[p]
        return out, acc

    def hierarchy(self, Om, kmax):
        """Y, L of orders < kmax + 1 and G_k, k <= kmax, for the given Omega list"""
        N, P = self.N, self.P
        grid = lambda fn: [[fn(p, i) for i in range(N)] for p in range(P)]  # noqa: E731
        Y = [grid(lambda p, i: mp.mpf(1))]
        L = []
        G = {}
        for k in range(kmax + 1):
            if k >= 1:
                F = grid(lambda p, i: 2 * self.v[p][i] * sum(L[k - 1 - j][p][i] * self.l[p][i] ** j / mp.factorial(j)
                                                             for j in range(k)))
                Y.append(self.cum_left(F))
            thr = mp.mpf(1) / 2 if k == 0 else Om[k - 1]
            s = grid(lambda p, i: hk(k, Om, self.l[p][i], self.R[p][i], [y[p][i] for y in Y],
                                     [lv[p][i] for lv in L] if k else None) - thr)
            Cint = grid(lambda p, i: s[p][i] / self.Rm1[p][i] * self.R[p][i] / self.u[p][i] * 2 * self.v[p][i])
            Cr, G[k] = self.cum_right(Cint)
            L.append(grid(lambda p, i: -2 * thr - self.u[p][i] * Cr[p][i]))
        return Y, L, G


def run_B(N, w, V, dps):
    t0 = time.time()
    if MUTATE == "coarse_B":
        N = 10
    if MUTATE == "short_tail":
        V = 6
    out = {}
    with mp.workdps(dps):
        B = PanelB(N, w, V)
        Om = [mp.mpf(0)] * 4
        for k in (2, 3, 4):
            Om_in = list(Om)
            if MUTATE == "wrong_Om0" and k >= 3:
                Om_in[0] = -4 * mp.ln2 + mp.mpf(1) / 100
            if MUTATE == "wrong_Om1" and k >= 4:
                Om_in[1] = 2 * mp.pi**2 / 3 + 8 * mp.ln2**2 + mp.mpf(1) / 100
            val, slope = solve_omega(lambda kk, o: B.hierarchy(o, kk)[2][kk], k, Om_in)
            Om[k - 2] = val
            out[f"Om{k - 2}"], out[f"slope{k}"] = +val, +slope
        out["secs"] = time.time() - t0
        out["grid"] = f"N={N} w={mp.nstr(mp.mpf(w), 6)} V={V} panels={B.P}"
    return out


# ---------------------------------------------------------------------------------------------------------------
# numerical checks of the closed forms (N1 .. N4)
# ---------------------------------------------------------------------------------------------------------------
def guarded(f, g):
    """f evaluated with g guard digits, rounded back to the working precision"""
    if g == 0:
        return f

    def w(*a):
        with mp.extradps(g):
            r = f(*a)
        return +r

    return w


def both(f):
    """int_0^{pi/2} f(th, x, tail) dth: Gauss-Legendre on [0, 1], tanh-sinh in x = pi/2 - th beyond"""
    return (mp.quad(lambda t: f(t, mp.pi / 2 - t, False), [0, mp.mpf(1) / 2, 1], method="gauss-legendre")
            + mp.quad(lambda x: f(mp.pi / 2 - x, x, True), [0, mp.pi / 2 - 1]))


def check_closed():
    t0 = time.time()
    mp.mp.dps = 45
    pts = [mp.mpf(s) for s in ("0.3", "0.8", "1.2", "1.5")]
    Om = [-4 * mp.ln2, 2 * mp.pi**2 / 3 + 8 * mp.ln2**2, mp.mpf(0), mp.mpf(0)]
    # N1
    worst = mp.mpf(99)
    for th in pts:
        a = mp.pi - 2 * th
        Kq = mp.quad(lambda p: -2 * mp.log(mp.cos(p)), [0, th])
        Tq = mp.quad(lambda p: p * mp.log(mp.cos(p)), [0, th])
        worst = min(worst, digits(Kq, K_of(th, a)), digits(Tq, T_of(th, a)))
    report("N1", worst >= 30, f"K and T closed forms against quadrature at th = 0.3, 0.8, 1.2, 1.5: worst {float(worst):.1f} digits")
    # N2
    worst2 = mp.mpf(99)
    for th in pts:
        Y3 = closed(th, mp.pi / 2 - th, False, Om)[6][3]

        def F2(p):
            c, lc, l, u, R, s2, Y, L, K, T = closed(p, mp.pi / 2 - p, False, Om)
            return (L[2] + L[1] * l + L[0] * l**2 / 2) * 2 * u

        worst2 = min(worst2, digits(mp.quad(F2, [0, th]), Y3))
    report("N2", worst2 >= 30, f"Y_3 closed form against the integral of F_2 at the same points: worst {float(worst2):.1f} digits")
    # N3 as declared (no guard digits) and N3g (A1: 30 guard digits)
    for cid, g in (("N3", 0), ("N3g", 30)):
        worst3 = mp.mpf(99)
        for th in pts:
            def src(p):
                c, lc, l, u, R, s2, Y, L, K, T = closed(p, mp.pi / 2 - p, False, Om)
                return 2 * (hk(2, Om, l, R, Y, L, om_last=False)) / s2

            L2q = -2 * Om[1] + mp.tan(th) * mp.quad(guarded(src, g), [0, th])
            worst3 = min(worst3, digits(L2q, closed(th, mp.pi / 2 - th, False, Om)[7][2]))
        report(cid, worst3 >= 30, f"L_2 closed form against -2 h_2(1) + u C_2 by quadrature, {g} guard digits: worst {float(worst3):.1f} digits")
    # N4a: (E2) at k - 1 = 2, closed L_2 against 2 int h_2 ln cos / sin^2
    def f2(th, x, tail):
        c, lc, l, u, R, s2, Y, L, K, T = closed(th, x, tail, Om)
        return 2 * hk(2, Om, l, R, Y, L) * lc / s2

    lhs, rhs = both(lambda th, x, tail: closed(th, x, tail, Om)[7][2]), both(f2)
    d2 = digits(lhs, rhs)
    report("N4a", d2 >= 30, f"(E2) at k - 1 = 2: int L_2 dth = {nd(lhs, 30)} against 2 int h_2 ln cos/sin^2 = {nd(rhs, 30)}: "
           f"{float(d2):.1f} digits ({time.time() - t0:.0f} s)")
    # N5: the closed form of L_3 (with Lc2 by quadrature) against -2 Omega2 + u C_3 by quadrature, Omega2 = 0, 30 guard digits
    kc = 16 if MUTATE == "wrong_L3" else 32
    q = mp.ln2
    worst5 = mp.mpf(99)
    for th in (mp.mpf("0.4"), mp.mpf("1.0"), mp.mpf("1.4")):
        def h3f(p):
            c, lc, l, u, R, s2, Y, L, K, T = closed(p, mp.pi / 2 - p, False, Om)
            return hk(3, Om, l, R, Y, L), s2, u, lc

        def src(p):
            h, s2, u, lc = h3f(p)
            return 2 * h / s2

        with mp.extradps(30):
            c, lc, l, u, R, s2, Y, L, K, T = closed(th, mp.pi / 2 - th, False, Om)
            Cn = mp.clsin(2, mp.pi - 2 * th)
            Lc2 = mp.quad(lambda p: mp.log(mp.cos(p)) ** 2, [0, th])
            sc = mp.sin(th) * mp.cos(th)
            B = (-32 * (lc + q) * Cn + 16 * s2 * Cn + kc * Lc2 + 16 * sc * lc * (lc + 2 * q) - 4 * Om[1] * th - 2 * Om[1] * sc)
            L3c = -2 * hk(3, Om, l, R, Y, L) + 2 * u * B
            L3q = u * mp.quad(guarded(src, 30), [0, th])
        worst5 = min(worst5, digits(L3c, L3q))
    report("N5", worst5 >= 30, f"L_3 closed form (with Lc2 = int_0^th ln^2 cos by quadrature) against -2 h_3(1) + u C_3 by quadrature at th = 0.4, 1.0, 1.4, "
           f"Omega2 = 0: worst {float(worst5):.1f} digits")


def Lc2_poly(th):
    """Lc2(th) = int_0^th ln^2 cos = (F(t) - F(1))/(2i), t = e^{2 i th}, F built from Li2, Li3 of 1/(1 + t) and -t"""
    q = mp.ln2
    sg = -1 if MUTATE == "wrong_Lc2" else 1

    def F(t):
        L1, Lt = mp.log(1 + t), mp.log(t)
        w = 1 / (1 + t)
        return (Lt * L1**2 - mp.mpf(2) / 3 * L1**3 - 2 * L1 * mp.polylog(2, w) - 2 * mp.polylog(3, w) + q**2 * Lt + Lt**3 / 12
                + 2 * q * mp.polylog(2, -t) + Lt * mp.polylog(2, -t) - sg * mp.polylog(3, -t) + q * Lt**2 / 2)

    return ((F(mp.exp(2j * th)) - F(mp.mpf(1))) / 2j).real


def check_lc2():
    """N6: Lc2 in polylogarithms against quadrature, and L_3 closed (with that Lc2) against -2 h_3(1) + u C_3 by quadrature"""
    t0 = time.time()
    mp.mp.dps = 45
    pts = [mp.mpf(s) for s in ("0.3", "0.9", "1.4")]
    Om = [-4 * mp.ln2, 2 * mp.pi**2 / 3 + 8 * mp.ln2**2, mp.mpf(0), mp.mpf(0)]
    q = mp.ln2
    worst_a, worst_b = mp.mpf(99), mp.mpf(99)
    for th in pts:
        def src(p):
            c2, lc2, l2, u2, R2, s22, Y2, L2, K2, T2 = closed(p, mp.pi / 2 - p, False, Om)
            return 2 * hk(3, Om, l2, R2, Y2, L2) / s22

        Lq = mp.quad(guarded(lambda p: mp.log(mp.cos(p)) ** 2, 30), [0, th])
        L3q = mp.tan(th) * mp.quad(guarded(src, 30), [0, th])
        with mp.extradps(30):
            Lp = Lc2_poly(th)
            c, lc, l, u, R, s2, Y, L, K, T = closed(th, mp.pi / 2 - th, False, Om)
            Cn = mp.clsin(2, mp.pi - 2 * th)
            sc = mp.sin(th) * mp.cos(th)
            B = (-32 * (lc + q) * Cn + 16 * s2 * Cn + 32 * Lp + 16 * sc * lc * (lc + 2 * q) - 4 * Om[1] * th - 2 * Om[1] * sc)
            L3c = -2 * hk(3, Om, l, R, Y, L) + 2 * u * B
        worst_a = min(worst_a, digits(Lq, +Lp))
        worst_b = min(worst_b, digits(+L3c, L3q))
    report("N6", worst_a >= 30 and worst_b >= 30, f"Lc2 = int_0^th ln^2 cos in Li2, Li3 of 1/(1 + e^(2 i th)) and -e^(2 i th) against quadrature at th = 0.3, 0.9, 1.4: worst "
           f"{float(worst_a):.1f} digits; L_3 closed (that Lc2) against -2 h_3(1) + u C_3 by quadrature: worst {float(worst_b):.1f} digits ({time.time() - t0:.0f} s)")


def check_n4b(cid, g):
    """(E2) at k - 1 = 3: int L_3 dth with L_3 = -2 h_3(1) + u C_3 by nested quadrature (left sum up to th = 1, right sum beyond,
    G_3 = 0) against 2 int h_3 ln cos / sin^2; Omega2 = 0; working precision 30 plus g guard digits in the integrands"""
    t0 = time.time()
    mp.mp.dps = 30
    Om30 = [-4 * mp.ln2, 2 * mp.pi**2 / 3 + 8 * mp.ln2**2, mp.mpf(0), mp.mpf(0)]

    def h3_of(th, x, tail):
        c, lc, l, u, R, s2, Y, L, K, T = closed(th, x, tail, Om30)
        return hk(3, Om30, l, R, Y, L), s2, lc

    def src3(th, x, tail):
        h, s2, lc = h3_of(th, x, tail)
        return 2 * h / s2   # h_3(1) = Omega2 = 0 here

    src3g = guarded(src3, g)

    def L3(th, x, tail):
        if not tail:
            return mp.tan(th) * mp.quad(lambda p: src3g(p, mp.pi / 2 - p, False), [0, th])
        return -mp.cot(x) * mp.quad(lambda y: src3g(mp.pi / 2 - y, y, True), [0, x])

    def rhs_f(th, x, tail):
        h, s2, lc = h3_of(th, x, tail)
        return 2 * h * lc / s2

    lhs, rhs = both(L3), both(guarded(rhs_f, g))
    d3 = digits(lhs, rhs)
    report(cid, d3 >= 25, f"(E2) at k - 1 = 3, {g} guard digits: int L_3 dth = {nd(lhs, 25)} (L_3 by nested quadrature, dps 30) against "
           f"2 int h_3 ln cos/sin^2 = {nd(rhs, 25)}: {float(d3):.1f} digits ({time.time() - t0:.0f} s)")


# ---------------------------------------------------------------------------------------------------------------
# the reduction of G_k to log-cosine moments, by hand (H0 .. H8)
# ---------------------------------------------------------------------------------------------------------------
def reduced_pieces(k):
    """At Omega0 = -4 ln 2 (w0) the closed forms need only ell = ln cos, u = tan, C = Cl2(pi - 2 th), D3 = Cl3(pi - 2 th) - Cl3(pi):
         L2 = -2 w1 - 4 w0 ell + 16 ell^2 + 16 u C,   Y3 = 4 w1 ell - 4 ell^3 + 32 ell^2 - 16 w0 ell + 32 u C - 16 D3.
       The integrand of G_k (E1 with E2, Omega_{k-1} = 0) is A(ell)/S + B(ell) + gamma C + delta D3, S = sin^2 th."""
    lg, S, u, C, D3, w0, w1, w2 = sp.symbols("ell S u C D3 w0 w1 w2")
    c2s = 1 - S
    L = [-1, -2 * w0 + 4 * lg, -2 * w1 - 4 * w0 * lg + 16 * lg**2 + 16 * u * C]
    Y = [1, 2 * lg, 4 * w0 * lg - 6 * lg**2, 4 * w1 * lg - 4 * lg**3 + 32 * lg**2 - 16 * w0 * lg + 32 * u * C - 16 * D3]
    o = [w0, w1, w2, 0]
    ell_c2 = c2s
    c = {1: w0 - 3 * lg - 1 + ell_c2,
         2: -2 * w0 * lg + w1 + lg**2 - 2 * lg + 2 * lg * ell_c2,
         3: 2 * w0 * lg**2 - 2 * w1 * lg + w2 - 2 * lg**3 - 2 * lg**2 + 2 * lg**2 * ell_c2,
         4: -sp.Rational(4, 3) * w0 * lg**3 + 2 * w1 * lg**2 - 2 * w2 * lg + lg**4 / 3 - sp.Rational(4, 3) * lg**3
         + sp.Rational(4, 3) * lg**3 * ell_c2}
    h = S * L[k - 2] + sum(Y[j] * c[k - 1 - j] for j in range(1, k - 1)) + c[k - 1] + Y[k - 1] / 2
    I = 8 * lg * h / S + 2 * sum(L[k - 1 - j] * (-2 * lg) ** j / sp.factorial(j) for j in range(1, k)) \
        + 2 * (sum(Y[j] * c[k - j] for j in range(1, k)) + c[k]) / S
    I = sp.expand(I)
    if k == 3:
        I = sp.expand(I.subs(w2, 0))
    ell_part = sp.expand(I.subs({C: 0, D3: 0}) * S)
    return dict(I=I, A=ell_part.coeff(S, 0), B=ell_part.coeff(S, 1), extra=ell_part.coeff(S, 2), gamma=I.coeff(C, 1),
                delta=I.coeff(D3, 1), syms=(lg, S, u, C, D3, w0, w1, w2))


def moments():
    """J_n = int_0^{pi/2} ln^n cos, n <= 4, from sqrt(pi)/2 Gamma((a+1)/2)/Gamma(a/2+1) = sum J_n a^n/n!"""
    a = sp.Symbol("a")
    ser = sp.series(sp.sqrt(sp.pi) / 2 * sp.gamma((a + 1) / 2) / sp.gamma(a / 2 + 1), a, 0, 5).removeO()
    J = [sp.simplify(sp.expand_func(sp.factorial(n) * ser.coeff(a, n))) for n in range(5)]
    if MUTATE == "wrong_J3":
        J[3] = J[3] + sp.pi * sp.zeta(3) / 8
    return J


def moment_integral(A, B, J, lg):
    """int_0^{pi/2} (A(ell)/sin^2 + B(ell)) dth with A(0) = 0, from int ell^b/sin^2 = -b J_{b-1} (H5)"""
    PA, PB = sp.Poly(A, lg), sp.Poly(B, lg)
    tot = 0
    for (b,), cf in PA.terms():
        assert b >= 1, "A has a constant term: the integral diverges at the throat"
        tot += -b * cf * J[b - 1]
    for (b,), cf in PB.terms():
        tot += cf * J[b]
    return tot


def basics(th, x, tail):
    a, c, s2, u = cl_args(th, x, tail)
    return mp.log(c), s2, u, mp.clsin(2, a), mp.clcos(3, a) + 3 * mp.zeta(3) / 4


def check_reduce():
    t0 = time.time()
    q = sp.log(2)
    lg, S, u, C, D3, w0, w1, w2 = sp.symbols("ell S u C D3 w0 w1 w2")
    # H0: the closed forms of the raw README with K = 2 th ln 2 - C, T = th C/2 - D3/4 - (ln 2/2) th^2, l = -2 ell, w0 = -4 ln 2
    th, qq = sp.symbols("th q")
    Ks, Ts, ls = 2 * th * qq - C, th * C / 2 - D3 / 4 - qq * th**2 / 2, -2 * lg
    W0 = -4 * qq
    L2r = -2 * w1 + 2 * W0 * ls + 4 * ls**2 - 8 * W0 * u * th - 16 * u * Ks
    Y3r = (-2 * w1 * ls + ls**3 / 2 + 8 * ls**2 + 8 * W0 * ls + 8 * W0 * th**2 - 16 * W0 * u * th - 32 * u * Ks + (16 if MUTATE == "wrong_Y3" else 32) * Ks * th
           + 64 * Ts)
    L2s = -2 * w1 - 4 * W0 * lg + 16 * lg**2 + 16 * u * C
    Y3s = 4 * w1 * lg - 4 * lg**3 + 32 * lg**2 - 16 * W0 * lg + 32 * u * C - 16 * D3
    r = [sp.expand(L2r - L2s), sp.expand(Y3r - Y3s)]
    report("H0", all(x == 0 for x in r), f"L_2 and Y_3 of the raw README against their forms in ell, u, C, D3 at Omega0 = -4 ln 2: residuals {[str(x) for x in r]}")
    P = {k: reduced_pieces(k) for k in (3, 4)}
    mp.mp.dps = 45
    W0n, W1n = -4 * mp.ln2, 2 * mp.pi**2 / 3 + 8 * mp.ln2**2
    Om = [W0n, W1n, mp.mpf(0), mp.mpf(0)]
    # H1: the decomposition against the integrand of A at four points
    worst = mp.mpf(99)
    for k in (3, 4):
        pc = P[k]
        fn = sp.lambdify(pc["syms"], pc["A"] / S + pc["B"] + pc["gamma"] * C + pc["delta"] * D3, "mpmath")
        fA = A_integrand(k, Om)
        for t in ("0.2", "0.7", "1.1", "1.4"):
            th_ = mp.mpf(t)
            ell, s2, u_, Cn, D3n = basics(th_, mp.pi / 2 - th_, False)
            worst = min(worst, digits(fA(th_, mp.pi / 2 - th_, False), fn(ell, s2, u_, Cn, D3n, W0n, W1n, mp.mpf(0))))
    report("H1", worst >= 30, f"A(ell)/S + B(ell) + gamma C + delta D3 against the integrand of method A, k = 3, 4, four points: worst {float(worst):.1f} digits")
    # H2: the D3 term by parts: delta = -32 P', P = cot (ell - w0): int delta D3 = 64 int C P
    g4 = P[4]["gamma"]
    d4 = P[4]["delta"]
    sg = -1 if MUTATE == "wrong_P" else 1
    rd = sp.simplify(sp.together(d4 - (-32) * (sg * w0 - S - lg) / S))
    report("H2a", rd == 0, f"delta = -32 (w0 - S - ell)/S = -32 dP/dth for P = cot th (ell - w0): residual {rd}")

    def fD(th_, x, tail):
        ell, s2, u_, Cn, D3n = basics(th_, x, tail)
        return (-32 * (W0n - s2 - ell) / s2) * D3n

    def fP(th_, x, tail):
        ell, s2, u_, Cn, D3n = basics(th_, x, tail)
        return 64 * Cn * (ell - (-W0n if MUTATE == "wrong_P" else W0n)) / u_

    lhs, rhs = both(guarded(fD, 30)), both(guarded(fP, 30))
    report("H2b", digits(lhs, rhs) >= 30, f"int delta D3 dth = {nd(lhs, 30)} against 64 int C cot(th) (ell - w0) dth = {nd(rhs, 30)} (boundary terms vanish): {float(digits(lhs, rhs)):.1f} digits")
    # H3: C (gamma + 64 P) = 64 (w0 - 1) u C
    t_ = sp.Symbol("t", positive=True)
    gam = g4.subs(S, t_**2 / (1 + t_**2)).subs(u, t_)
    r3 = sp.simplify(sp.together(gam + 64 * (lg - sg * w0) / t_ - 64 * (w0 - 1) * t_))
    report("H3", r3 == 0, f"gamma + 64 cot(th)(ell - w0) - 64 (w0 - 1) u = {r3} (sympy, S = u^2/(1 + u^2))")
    # H4: int tan C = pi^3/12
    J = moments()
    Jn = [mp.mpf(sp.N(j, 50)) for j in J]
    th0 = mp.mpf("0.7")
    dC = mp.diff(lambda t: mp.clsin(2, mp.pi - 2 * t), th0)
    c_ok = digits(dC, 2 * mp.log(mp.cos(th0)) + 2 * mp.ln2)
    xs = mp.mpf("1e-25")
    end = abs(mp.log(mp.sin(xs)) * mp.clsin(2, 2 * xs))
    tanC = mp.pi**3 / 6 if MUTATE == "wrong_tanC" else mp.pi**3 / 12
    def utimesC(th_, x, tail):
        b = basics(th_, x, tail)
        return b[2] * b[3]

    num = both(utimesC)
    byparts = 2 * Jn[2] + 2 * mp.ln2 * Jn[1]
    report("H4", c_ok >= 30 and end < mp.mpf("1e-20") and digits(num, tanC) >= 30 and digits(byparts, tanC) >= 30,
           f"int tan th Cl2(pi - 2 th) dth = {nd(num, 30)}; C' = 2 ell + 2 ln 2 to {float(c_ok):.1f} digits; boundary |ell C| = {nd(end, 2)} at x = 1e-25; "
           f"by parts 2 J_2 + 2 ln 2 J_1 = {nd(byparts, 30)} against pi^3/12 (as coded {nd(tanC, 12)})")
    # H5, H6: moments
    def ellb(b):
        def f(th_, x, tail):
            ell, s2, _, _, _ = basics(th_, x, tail)
            return ell**b / s2

        return f

    d5 = min(digits(both(guarded(ellb(b), 30)), -b * Jn[b - 1]) for b in range(1, 5))
    report("H5", d5 >= 30, f"int ell^b/sin^2 dth = -b J_(b-1), b = 1 .. 4, by quadrature: worst {float(d5):.1f} digits")
    def ellp(n):
        return lambda th_, x, tail: basics(th_, x, tail)[0] ** n

    d6 = min(digits(both(guarded(ellp(n), 30)), Jn[n]) for n in range(5))
    report("H6", d6 >= 30, f"J_n = int ell^n dth, n <= 4, from the Beta function series (sympy) against quadrature: worst {float(d6):.1f} digits; "
           f"J_3 = {sp.simplify(J[3])}")
    # H7: assembly
    W0s = -4 * q + (sp.Rational(1, 100) if MUTATE == "wrong_Om0" else 0)
    G3 = moment_integral(P[3]["A"], P[3]["B"], J, lg)
    sol = sp.solve(G3.subs(w0, W0s), w1)
    Om1 = sp.simplify(sol[0])
    G4 = moment_integral(P[4]["A"].subs(w2, 0), P[4]["B"], J, lg) + 64 * (w0 - 1) * (sp.pi**3 / 6 if MUTATE == "wrong_tanC" else sp.pi**3 / 12)
    G4 = G4.subs(w0, W0s)
    Om2 = sp.simplify(sp.expand((G4 / (4 * sp.pi)).subs(w1, Om1)))
    t1 = sp.simplify(Om1 - (2 * sp.pi**2 / 3 + 8 * q**2))
    t2 = sp.simplify(sp.expand(Om2 - (80 * sp.zeta(3) - sp.Rational(8, 3) * sp.pi**2 * q - sp.Rational(32, 3) * q**3)))
    report("H7", t1 == 0 and t2 == 0, f"G_3 = 0 gives Omega1 = {Om1}; G_4(Omega2 = 0)/(4 pi) with that Omega1 gives Omega2 = {Om2}; "
           f"against 2 pi^2/3 + 8 ln^2 2 and 80 zeta(3) - (8/3) pi^2 ln 2 - (32/3) ln^3 2: differences {t1} and {t2}")
    # H8: the derived G_3(0), G_4(0) against method A at 40 digits
    with mp.workdps(40 + GUARD):
        Oms = [W0n, mp.mpf(0), mp.mpf(0), mp.mpf(0)]
        gA3 = A_G(3, Oms)
        Oms4 = [W0n, W1n, mp.mpf(0), mp.mpf(0)]
        gA4 = A_G(4, Oms4)
        gd3 = mp.mpf(sp.N(G3.subs({w0: -4 * q, w1: 0}), 60))
        gd4 = mp.mpf(sp.N(G4.subs({w1: 2 * sp.pi**2 / 3 + 8 * q**2}), 60))
        e3, e4 = digits(gA3, gd3), digits(gA4, gd4)
    report("H8", min(e3, e4) >= 30, f"derived G_3(0) = {nd(gd3, 30)}, G_4(0) = {nd(gd4, 30)} against method A's: {float(e3):.1f} and {float(e4):.1f} digits "
           f"({time.time() - t0:.0f} s)")


# ---------------------------------------------------------------------------------------------------------------
# evaluations, comparisons, PSLQ
# ---------------------------------------------------------------------------------------------------------------
def pslq_found(x, tau, H, expected=None):
    r = mp.pslq(x, tol=tau, maxcoeff=H, maxsteps=10**6)
    if r is None:
        return False, None, None
    res = abs(sum(c * v for c, v in zip(r, x)))
    ok = res <= tau and max(abs(c) for c in r) <= H
    r10 = mp.pslq(x, tol=10 * tau, maxcoeff=H, maxsteps=10**6)
    ok = ok and r10 is not None and (list(r10) == list(r) or list(r10) == [-c for c in r])
    if expected is not None:
        ok = ok and (list(r) == expected or list(r) == [-c for c in expected])
    return ok, r, res


def check_evals():
    mp.mp.dps = 80
    closed_Om = [-4 * mp.ln2, 2 * mp.pi**2 / 3 + 8 * mp.ln2**2]
    om2_closed = 80 * mp.zeta(3) - mp.mpf(8) / 3 * mp.pi**2 * mp.ln2 - mp.mpf(32) / 3 * mp.ln2**3
    plan = [("A50", lambda: run_A(50)), ("A70", lambda: run_A(70)),
            ("B1", lambda: run_B(40, mp.mpf(1) / 2, 18, 120)), ("B2", lambda: run_B(50, mp.mpf(2) / 5, 20, 140))]
    ev = {}
    for name, fn in plan:
        ev[name] = fn()
        e = ev[name]
        print(f"  {name}: {e.get('grid', '')} {e['secs']:.0f} s\n    Omega0 = {nd(e['Om0'], 60)}\n    Omega1 = {nd(e['Om1'], 60)}"
              f"\n    Omega2 = {nd(e['Om2'], 60)}\n    slopes G2, G3, G4 = {nd(e['slope2'], 20)} {nd(e['slope3'], 20)} {nd(e['slope4'], 20)}",
              flush=True)
    mp.mp.dps = 80
    # P1
    rows, ok1 = [], True
    for n, e in ev.items():
        d0, d1 = digits(closed_Om[0], e["Om0"]), digits(closed_Om[1], e["Om1"])
        rows.append(f"{n} {float(d0):.1f}/{float(d1):.1f}")
        ok1 = ok1 and d0 >= 30 and d1 >= 30
    report("P1", ok1, "planted: digits of Omega0 and Omega1 against -4 ln 2 and 2 pi^2/3 + 8 ln^2 2: " + ", ".join(rows))
    # P2
    rows, ok2 = [], True
    for n, e in ev.items():
        ds = [digits(-4 * mp.pi, e[f"slope{k}"]) for k in (2, 3, 4)]
        rows.append(f"{n} " + "/".join(f"{float(d):.1f}" for d in ds))
        ok2 = ok2 and min(ds) >= 30
    report("P2", ok2, "slope of G_k (k = 2, 3, 4) against -4 pi, digits: " + ", ".join(rows))
    # P3: the control wrong_Om0 must move Omega1 (reported when the mutation is on; the check is its statement)
    if MUTATE == "wrong_Om0":
        moved = all(abs(ev[n]["Om1"] - closed_Om[1]) > mp.mpf(1) / 1000 for n in ("A50", "B1"))
        report("P3", moved, "control wrong_Om0 moves Omega1 of A50 and B1: " + ", ".join(
            f"{n} {nd(ev[n]['Om1'] - closed_Om[1], 6)}" for n in ("A50", "B1")))
    # M1 .. M3
    pairs = [("M1", "A50", "A70", 20), ("M2", "B1", "B2", 20), ("M3", "A70", "B2", 20)]
    for cid, a, b, bar in pairs:
        d = digits(ev[a]["Om2"], ev[b]["Om2"])
        report(cid, d >= bar, f"Omega2 of {a} against {b}: {float(d):.1f} digits (bar {bar}); |diff| = {nd(abs(ev[a]['Om2'] - ev[b]['Om2']), 3)}")
    cross = {f"{a}-{b}": float(digits(ev[a]["Om2"], ev[b]["Om2"])) for a, b in (("A50", "B1"), ("A50", "B2"), ("A70", "B1"), ("A70", "B2"))}
    print("  reported: digits of Omega2 between the other pairs: " + ", ".join(f"{k} {v:.1f}" for k, v in cross.items()))
    # X1, X2
    om2 = ev["A70"]["Om2"]
    ref24 = mp.mpf("74.3693867279333175953847")
    ref38 = mp.mpf("74.369386727933317595384670295846617267")
    d24, d38 = abs(om2 - ref24), abs(om2 - ref38)
    report("X1", d24 <= mp.mpf("5e-23"), f"Omega2 (A70) = {nd(om2, 50)}; against #195's fit 74.3693867279333175953847: {nd(d24, 3)} "
           f"(bar 5e-23); against its 38-digit median: {nd(d38, 3)} (its range 3.2e-40)")
    bars = {"A70": 40, "A50": 40, "B1": 30, "B2": 30}
    rows, okx = [], True
    for n, bar in bars.items():
        d = digits(om2_closed, ev[n]["Om2"])
        rows.append(f"{n} {float(d):.1f}")
        okx = okx and d >= bar
    report("X2", okx, "digits of Omega2 against 80 zeta(3) - (8/3) pi^2 ln 2 - (32/3) ln^3 2 = " + nd(om2_closed, 50) + ": " + ", ".join(rows))
    # Q1, Q2
    agree = min(60, int(digits(ev["A70"]["Om2"], ev["B2"]["Om2"])))
    pert = mp.mpf(1) / 1000
    with mp.workdps(max(agree, 30)):     # a value that agrees with B in fewer digits is still given to PSLQ, which then fails
        om = +ev["A70"]["Om2"]
        z3 = mp.zeta(3) + (pert if MUTATE == "pslq_wrong_basis" else 0)
        b1 = [om, z3, mp.pi**2 * mp.ln2, mp.ln2**3]
        f1, r1, res1 = pslq_found(b1, mp.mpf("1e-18"), 1000, expected=[3, -240, 8, 32])
        c1, _, _ = pslq_found([om, mp.zeta(3) + pert, mp.pi**2 * mp.ln2, mp.ln2**3], mp.mpf("1e-18"), 1000)
        report("Q1", f1 and not c1, f"PSLQ on {{Omega2, zeta(3), pi^2 ln 2, ln^3 2}} at tol 1e-18, H 1000, Omega2 at {agree} digits: "
               f"relation {r1}, residual {nd(res1, 3) if res1 is not None else None}; perturbed basis finds a relation: {c1}")
        b2 = [om, z3, mp.pi**2 * mp.ln2, mp.ln2**3, mp.pi**2, mp.ln2**2, mp.ln2, mp.mpf(1)]
        f2, r2, res2 = pslq_found(b2, mp.mpf("1e-30"), 1000, expected=[3, -240, 8, 32, 0, 0, 0, 0])
        c2, _, _ = pslq_found([om, mp.zeta(3) + pert, mp.pi**2 * mp.ln2, mp.ln2**3, mp.pi**2, mp.ln2**2, mp.ln2, mp.mpf(1)], mp.mpf("1e-30"), 1000)
        report("Q2", f2 and not c2, f"PSLQ on the eight-element weight <= 3 basis at tol 1e-30, H 1000: relation {r2}, residual "
               f"{nd(res2, 3) if res2 is not None else None}; perturbed basis finds a relation: {c2}")
    return ev


# ---------------------------------------------------------------------------------------------------------------
# the remainder against #195's W^2
# ---------------------------------------------------------------------------------------------------------------
def read_w2(path):
    """{D: W2 string} of the growth.csv-D table (crate qd column) and of the 64-value list of w2_output.txt"""
    growth, lst = {}, {}
    for line in open(path):
        m = re.match(r"\s*(\d+)\s+(?:modes_r\.py|TS13 Table I)\s+(\S+)\s+(\S+)\s", line)
        if m:
            growth[int(m.group(1))] = m.group(3)
        m = re.match(r"\s*D=(\d+)\s+(\S+)\s+digits qd/od", line)
        if m:
            lst[int(m.group(1))] = m.group(2)
    return growth, lst


def check_remainder(ev):
    mp.mp.dps = 60
    om1, om2 = ev["A70"]["Om1"], ev["A70"]["Om2"] + (mp.mpf(1) / 100 if MUTATE == "wrong_Om2_series" else 0)
    txt = open(W2FILE).read()
    om3 = mp.mpf(re.search(r"Omega3 = (\S+)", txt).group(1))
    om4 = mp.mpf(re.search(r"Omega4 = (\S+)", txt).group(1))
    growth, lst = read_w2(W2FILE)
    x = {}
    print("  four-term series s - 4 ln 2 + Omega1/s + Omega2/s^2 against W^2 (qd, w2_output.txt, growth table), s = 2(D - 3):")
    print("  D       s       remainder             rem*s^3           (rem*s^3 - Omega3)")
    for D in (200, 500, 1000, 2000):
        s = mp.mpf(2 * (D - 3))
        rem = mp.mpf(growth[D]) - (s - 4 * mp.ln2 + om1 / s + om2 / s**2)
        x[D] = rem * s**3
        print(f"  {D:<6d} {int(s):<7d} {nd(rem, 12):<21s} {nd(x[D], 10):<17s} {nd(x[D] - om3, 6)}")
    dev = [abs(x[D] - om3) for D in (200, 500, 1000, 2000)]
    falls = all(dev[i] > dev[i + 1] for i in range(3))
    rel = abs(x[2000] - om3) / abs(om3)
    report("R1", falls and rel < mp.mpf("0.01"), f"|rem*s^3 - Omega3| at D = 200, 500, 1000, 2000: "
           f"{', '.join(nd(d, 5) for d in dev)} (falling: {falls}); at D = 2000 {nd(100 * rel, 3)} % of |Omega3| (bar 1 %)")
    s1, s2 = mp.mpf(2 * (1000 - 3)), mp.mpf(2 * (2000 - 3))
    est = (s2 * x[2000] - s1 * x[1000]) / (s2 - s1)
    report("R2", abs(est - om3) <= mp.mpf("0.05"), f"two-point extrapolation from D = 1000, 2000: {nd(est, 10)}, Omega3 = {nd(om3, 14)}: "
           f"difference {nd(est - om3, 4)} (bar 0.05)")
    if W2RAW:
        print("  reported, not gated: rem*s^3 from the 60-digit qd values at D = 200 .. 2000 (raw I195-dh610/run_02, W2 strings)")
        for line in open(W2RAW):
            m = re.search(r'"D":(\d+),"W2":"([^"]+)"', line)
            if m and int(m.group(1)) in (200, 500, 1000, 2000):
                D = int(m.group(1))
                s = mp.mpf(2 * (D - 3))
                xx = (mp.mpf(m.group(2)) - (s - 4 * mp.ln2 + om1 / s + om2 / s**2)) * s**3
                print(f"  {D:<6d} rem*s^3 = {nd(xx, 14)}   (21-digit value above: {nd(x[D], 14)}; difference {nd(xx - x[D], 3)})")
    print("  reported, not gated: rem*s^3 at the 64 values of D (60-digit qd), and (rem*s^3 - Omega3)*s against Omega4 = -12794.5")
    print("  D       rem*s^3            (rem*s^3 - Omega3)*s")
    for D in sorted(lst):
        s = mp.mpf(2 * (D - 3))
        rem = mp.mpf(lst[D]) - (s - 4 * mp.ln2 + om1 / s + om2 / s**2)
        xx = rem * s**3
        print(f"  {D:<7d} {nd(xx, 12):<18s} {nd((xx - om3) * s, 8)}")


def main():
    mp.mp.dps = 60
    print(f"order4.py  mutate={MUTATE}  parts={PARTS}  mpmath {mp.__version__}  sympy {sp.__version__}", flush=True)
    t0 = time.time()
    if "sym" in PARTS:
        check_sym()
    if "closed" in PARTS:
        check_closed()
    if "lc2" in PARTS:
        check_lc2()
    if "reduce" in PARTS:
        check_reduce()
    if "n4b" in PARTS:
        check_n4b("N4b", 0)
    if "n4bg" in PARTS:
        check_n4b("N4bg", 30)
    mp.mp.dps = 60
    ev = None
    if "evals" in PARTS or "remainder" in PARTS:
        ev = check_evals()
    if "remainder" in PARTS:
        check_remainder(ev)
    nfail = sum(1 for _, ok in RESULTS if not ok)
    print(f"\n{len(RESULTS) - nfail} passed, {nfail} failed of {len(RESULTS)} checks; {time.time() - t0:.0f} s")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()

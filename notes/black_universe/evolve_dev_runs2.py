#!/usr/bin/env python3
"""More development runs of #217 (design stage), on the prototype evolve_dev.py; companion of evolve_dev_runs.py. Not the build.

    evolve_dev_runs2.py near C RHOC H X         mode family near threshold: raw and weighted finders and the departure, with the late state of the diagonal
    evolve_dev_runs2.py pulseE                  c = -1.5, pulse on the + ray: where the expansion event stops, and where the mode amplitude changes sign
    evolve_dev_runs2.py profile                 profile of the gauge-invariant Y from the truncation seed, c = -1, -1.01, -1.5
    evolve_dev_runs2.py resp C RHOC             Y on the diagonal after a compact pulse of amplitude +-1e-3 (+ ray) at h = 0.02, 0.01: the part odd in the amplitude (the response) and the even part (the seed)
    evolve_dev_runs2.py constr                  constraint violation C++, C-- (stencil of dn.constraints, minus the exact static member's) against Y, c = -1, -1.01, -1.5, -0.3
    evolve_dev_runs2.py tuned C RHOC H X        the mode family at A = A* (tuned): |Y| on the diagonal against t and its local growth rate, against W and kappa
"""
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "wormhole_evolution"))
import dn                      # noqa: E402
import evolve_dev as ed        # noqa: E402
import evolve_potential as ep  # noqa: E402

SQ2 = math.sqrt(2.0)


def fmt(ev, d):
    if not ev:
        return "none"
    key, t, l = ed.event_T(d, ev)
    return f"{key} t = {t:.2f} l = {l:+.2f}"


def near(c, rc, h, X, Astar=None):
    mem = ep.Member(c)
    dR, dL, info = ed.extract_mode(mem, rc, h, X + 20.0)
    N = int(round(X / h))
    par, st = ed.params(mem), ed.Static(mem, rc, N, h)
    print(f"c = {c}, rho_c = {rc}, h = {h}, X = {X}: mode extracted at t_s = {info['t_s']:.2f}, Y at the corner node {info['Yc']:.3e}")
    # the threshold: bisection of the sign of the departure (|Y| = 3e-3) as in evolve_dev_runs.mode
    def run(A):
        d = ed.mode_rays(mem, rc, N, h, A, dR, dL)
        r = ed.evolve_dep(d, par, st, Y0=3e-3)
        dp = r["dep"]
        return 0 if dp is None else dp["sign"]
    lo, hi = -1e-4, 1e-4
    fl = run(lo)
    while hi - lo > 1e-14:
        mid = 0.5 * (lo + hi)
        if run(mid) == fl:
            lo = mid
        else:
            hi = mid
    ps = 0.5 * (lo + hi)
    print(f"A*_h = {ps:.6e}; fate(A < A*) = {fl:+d}; finders: raw = the code's cell_class with tolerance 1e-2; weighted = 2 e^f r+ r- >= 0.02 with the signs")
    for e in (3.0, 4.0, 5.0, 6.0, 7.0, 8.0):
        for s in (-1, 1):
            A = ps + s * 10.0 ** (-e)
            d = ed.mode_rays(mem, rc, N, h, A, dR, dL)
            r1 = ed.evolve(d, par, stop=True, tol=1e-2, rexp=2.0)
            r2 = ed.evolve_dep(d, par, st, Y0=1e9, stop_dep=False, run_events=True, wtol=0.02)
            r3 = ed.evolve_dep(d, par, st, Y0=3e-3)
            dp = r3["dep"]
            print(f"   A = A* {'-' if s < 0 else '+'} 1e-{e:.0f}: raw {fmt(r1['events'], d)};  weighted {fmt(r2['events'], d)};  departure t = {dp['T']:.2f} sign {dp['sign']:+d}")
    # the late state of the diagonal for the two nearest runs
    for s in (-1, 1):
        A = ps + s * 10.0 ** (-8.0)
        d = ed.mode_rays(mem, rc, N, h, A, dR, dL)
        res = ed.evolve(d, par, stop=False, store=True)
        G = res["G"]
        print(f"   late diagonal state, A = A* {'-' if s < 0 else '+'} 1e-8 (cells (j, j); dead cells {int(res['dead'].sum())}):")
        for tt in (16.0, 20.0, 30.0, 40.0, 50.0):
            j = int(tt / (SQ2 * h))
            if j + 1 >= N:
                continue
            S_, E_, W_, N_ = (G[:, j, j], G[:, j + 1, j], G[:, j, j + 1], G[:, j + 1, j + 1])
            if res["dead"][j:j + 2, j:j + 2].any():
                print(f"      t = {tt}: excised")
                continue
            cc = 0.25 * (S_ + E_ + W_ + N_)
            rp = (E_ + N_ - S_ - W_)[0] / (2 * h)
            rm = (W_ + N_ - S_ - E_)[0] / (2 * h)
            print(f"      t = {tt}: r {cc[0]:.4f}, f {cc[1]:+.3f}, r+ {rp:+.3e}, r- {rm:+.3e}, 2 e^f r+ r- = {2 * math.exp(cc[1]) * rp * rm:+.3e}, 2M/r = {1 + 2 * math.exp(cc[1]) * rp * rm:.4f}")


def pulseE():
    c, rc, X = -1.5, 1.0, 50.0
    mem = ep.Member(c)
    for h in (0.02, 0.01):
        N = int(round(X / h))
        par = ed.params(mem)

        def run(A):
            d = ed.member_rays(mem, rc, N, h, pulse=dict(A=A, xc=3.0, w=0.5, ray="+"))
            res = ed.evolve(d, par, stop=True, tol=1e-2, rexp=2.0)
            ev = res["events"]
            return (1, ed.event_T(d, ev)[1]) if ev and "expand" in ev else ((-1, ed.event_T(d, ev)[1]) if ev else (0, None))
        lo, hi = -1e-3, 1e-3
        assert run(lo)[0] == 1 and run(hi)[0] == 0
        while hi - lo > 1e-14:
            mid = 0.5 * (lo + hi)
            if run(mid)[0] == 1:
                lo = mid
            else:
                hi = mid
        ps = 0.5 * (lo + hi)
        print(f"h = {h}: expansion event for A < A*_E = {ps:.6e}; time of the event for A = A*_E - 10^-e:")
        print("   " + " ".join(f"{e:.1f}:{run(ps - 10.0 ** (-e))[1]:.2f}" for e in np.arange(3.0, 11.6, 1.0)))
        # the mode amplitude's zero: Y at the diagonal at t = 12 is linear in A near A*_E
        st = ed.Static(mem, rc, N, h)
        Ys = []
        As = [ps - 1e-4, ps, ps + 1e-4]
        for A in As:
            d = ed.member_rays(mem, rc, N, h, pulse=dict(A=A, xc=3.0, w=0.5, ray="+"))
            old = np.ascontiguousarray(d.Sp.copy())
            new = np.empty_like(old)
            dold = np.zeros(N + 1, np.uint8)
            dnew = np.zeros(N + 1, np.uint8)
            j12 = int(12.0 / (SQ2 * h))
            for j in range(j12 + 1):
                new[:, 0] = d.Sm[:, j + 1]
                dnew[0] = 0
                ed.row_step(old, new, dold, dnew, h, par, d.eps, 0.05, 25.0, 1)
                c_, up, um = dn.strip_cells(old, new, h)
                if j == j12:
                    Ys.append((c_[0][j] - st.R[N]) - st.Q[N] * (c_[2][j] - st.P[N]))
                old, new = new, old
                dold, dnew = dnew, dold
        slope = (Ys[2] - Ys[0]) / (As[2] - As[0])
        A0 = ps - Ys[1] / slope
        print(f"   Y at the diagonal, t = 12, for A = A*_E - 1e-4, A*_E, + 1e-4: {Ys[0]:+.3e} {Ys[1]:+.3e} {Ys[2]:+.3e}; dY/dA = {slope:.2f}; the mode amplitude (the sign of Y) changes at A = {A0:.3e}")


def profile():
    h, X = 0.01, 22.0
    ls = [-8, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6]
    for c, rc in ((-1.0, 1.0), (-1.01, 1.0), (-1.5, 1.0)):
        mem = ep.Member(c)
        N = int(round(X / h))
        d = ed.member_rays(mem, rc, N, h)
        res = ed.evolve(d, ed.params(mem), stop=False, store=True)
        G = res["G"]
        ks = np.arange(-N, N + 1)
        rho, A, f = ed.tortoise_rho(mem, rc, ks * h / SQ2)
        i = np.arange(N + 1)[:, None]
        j = np.arange(N + 1)[None, :]
        kk = (i - j) + N
        rr = rho[kk]
        Y = (G[0] - mem.R(rho)[kk]) - (rr * mem.R(rr) / SQ2) * (G[2] - mem.phi(rho)[kk])
        t = (i + j) * h / SQ2 + 0 * kk
        l = (i - j) * h / SQ2 + 0 * kk

        def at(tc):
            out = []
            for l0 in ls:
                m = (np.abs(l - l0) < 0.6 * h / SQ2 + 1e-12) & (np.abs(t - tc) < 0.6 * h * SQ2)
                out.append(Y[m][0] if m.any() else np.nan)
            return np.array(out)
        P = {tc: at(tc) for tc in (10.0, 14.0, 18.0)}
        print(f"c = {c}: Y from the truncation seed, h = {h}")
        print("   l          " + " ".join(f"{x:+8.1f}" for x in ls))
        for tc in P:
            print(f"   Y(t = {tc:4.1f}) " + " ".join(f"{v:+.1e}"[:8] for v in P[tc]))
        print("   ratio 14/10  " + " ".join(f"{v:8.2f}" for v in P[14.0] / P[10.0]))
        if c == -1.0:
            Pm1 = P
        if c == -1.01:
            print("   profile at t = 14: Y(c = -1)/Y(c = -1.01) " + " ".join(f"{v:8.3f}" for v in Pm1[14.0] / P[14.0]))


def resp(c, rc):
    mem = ep.Member(c)
    X = 40.0
    out = {}
    for h in (0.02, 0.01):
        N = int(round(X / h))
        st = ed.Static(mem, rc, N, h)
        par = ed.params(mem)
        res = {}
        for Ap in (1e-3, -1e-3):
            d = ed.member_rays(mem, rc, N, h, pulse=dict(A=Ap, xc=3.0, w=0.5, ray="+"))
            old = np.ascontiguousarray(d.Sp.copy())
            new = np.empty_like(old)
            dold = np.zeros(N + 1, np.uint8)
            dnew = np.zeros(N + 1, np.uint8)
            Y = []
            for j in range(int(30.0 / (SQ2 * h))):
                new[:, 0] = d.Sm[:, j + 1]
                dnew[0] = 0
                ed.row_step(old, new, dold, dnew, h, par, d.eps, 0.05, 25.0, 1)
                c_, up, um = dn.strip_cells(old, new, h)
                Y.append((c_[0][j] - st.R[N]) - st.Q[N] * (c_[2][j] - st.P[N]))
                old, new = new, old
                dold, dnew = dnew, dold
            res[Ap] = np.array(Y)
        out[h] = ((res[1e-3] - res[-1e-3]) / 2, (res[1e-3] + res[-1e-3]) / 2)
    print(f"c = {c}, rho_c = {rc}: pulse of amplitude +-1e-3 on the + ray (centre 3, width 0.5); Y on the diagonal; odd part = (Y(+) - Y(-))/2 (the linear response), even part = (Y(+) + Y(-))/2 (the seed and second order)")
    tt = (4, 6, 8, 10, 12, 14, 16, 18, 20, 24)
    for h, (odd, even) in out.items():
        print(f"   h = {h}: t:      " + " ".join(f"{t:9d}" for t in tt))
        print(f"   h = {h}: odd:    " + " ".join(f"{odd[int(t / (SQ2 * h))]:+.2e}" for t in tt))
        print(f"   h = {h}: even:   " + " ".join(f"{even[int(t / (SQ2 * h))]:+.2e}" for t in tt))
    # local slopes of ln|odd| and ln|even| between t = 8 and 16 at h = 0.01
    for h in (0.02, 0.01):
        odd, even = out[h]
        for name, y in (("odd", odd), ("even", even)):
            m = np.array([k for k in range(len(y)) if 10.0 <= (2 * k + 1) * h / SQ2 <= 16.0])
            t = (2 * m + 1) * h / SQ2
            print(f"   h = {h}: slope of ln|{name}| over t = 10 to 16: {np.polyfit(t, np.log(np.abs(y[m]) + 1e-300), 1)[0]:+.3f}")


def resp2(c, rc, Aps=(1e-5, 1e-6)):
    mem = ep.Member(c)
    h, X = 0.01, 40.0
    N = int(round(X / h))
    st = ed.Static(mem, rc, N, h)
    par = ed.params(mem)

    def ser(Ap):
        d = ed.member_rays(mem, rc, N, h, pulse=dict(A=Ap, xc=3.0, w=0.5, ray="+"))
        old = np.ascontiguousarray(d.Sp.copy())
        new = np.empty_like(old)
        dold = np.zeros(N + 1, np.uint8)
        dnew = np.zeros(N + 1, np.uint8)
        Y = []
        for j in range(int(30.0 / (SQ2 * h))):
            new[:, 0] = d.Sm[:, j + 1]
            dnew[0] = 0
            ed.row_step(old, new, dold, dnew, h, par, d.eps, 0.05, 25.0, 1)
            c_, up, um = dn.strip_cells(old, new, h)
            Y.append((c_[0][j] - st.R[N]) - st.Q[N] * (c_[2][j] - st.P[N]))
            old, new = new, old
            dold, dnew = dnew, dold
        return np.array(Y)
    print(f"c = {c}, rho_c = {rc}, h = {h}: Y on the diagonal, odd part (Y(+Ap) - Y(-Ap))/2 per unit Ap (the linear response to the pulse); kappa = {mem.kappa:.4f}")
    tt = (4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28)
    print("   t:           " + " ".join(f"{t:9d}" for t in tt))
    for Ap in Aps:
        odd = (ser(Ap) - ser(-Ap)) / 2 / Ap
        print(f"   Ap = {Ap:.0e}: " + " ".join(f"{odd[int(t / (SQ2 * h))]:+.2e}" for t in tt))
        for ta, tb in ((10, 16), (16, 22), (22, 28)):
            m = np.array([k for k in range(len(odd)) if ta <= (2 * k + 1) * h / SQ2 <= tb])
            t = (2 * m + 1) * h / SQ2
            print(f"      slope of ln|odd| over t = {ta} to {tb}: {np.polyfit(t, np.log(np.abs(odd[m]) + 1e-300), 1)[0]:+.3f}")


def constr():
    h, X = 0.01, 24.0
    for c, rc in ((-1.0, 1.0), (-1.01, 1.0), (-1.5, 1.0), (-0.3, 0.0)):
        mem = ep.Member(c)
        N = int(round(X / h))
        d = ed.member_rays(mem, rc, N, h)
        res = ed.evolve(d, ed.params(mem), stop=False, store=True)
        G = res["G"]
        ks = np.arange(-N, N + 1)
        rho, A, f = ed.tortoise_rho(mem, rc, ks * h / SQ2)
        i = np.arange(N + 1)[:, None]
        j = np.arange(N + 1)[None, :]
        kk = (i - j) + N
        E = np.zeros_like(G)
        E[0], E[1], E[2] = mem.R(rho)[kk], f[kk], mem.phi(rho)[kk]

        def cons(G):
            r, ff, ph = G[0], G[1], G[2]
            rpp = (r[2:, :] - 2 * r[1:-1, :] + r[:-2, :]) / h ** 2
            fp = (ff[2:, :] - ff[:-2, :]) / (2 * h)
            rp = (r[2:, :] - r[:-2, :]) / (2 * h)
            pp = (ph[2:, :] - ph[:-2, :]) / (2 * h)
            Cpp = rpp + fp * rp + (r[1:-1, :] / 2) * (-1.0) * pp ** 2
            rmm = (r[:, 2:] - 2 * r[:, 1:-1] + r[:, :-2]) / h ** 2
            fm = (ff[:, 2:] - ff[:, :-2]) / (2 * h)
            rm = (r[:, 2:] - r[:, :-2]) / (2 * h)
            pm = (ph[:, 2:] - ph[:, :-2]) / (2 * h)
            Cmm = rmm + fm * rm + (r[:, 1:-1] / 2) * (-1.0) * pm ** 2
            return Cpp, Cmm
        Cp, Cm = cons(G)
        Cpe, Cme = cons(E)
        print(f"c = {c}: at the diagonal node (j, j); dC = C(numerical) - C(exact static member), same stencil; Y = dr - (rho R/sqrt2) dphi; df = f - f_exact")
        for tc in (6.0, 8.0, 10.0, 12.0, 14.0):
            jj = int(round(tc / (SQ2 * h)))
            if jj < 2 or jj >= N - 2:
                continue
            # Cpp is [i-1 index shifted]: Cp[i-1, j] is C++ at node (i, j); Cm[i, j-1] is C-- at node (i, j)
            dCp = Cp[jj - 1, jj] - Cpe[jj - 1, jj]
            dCm = Cm[jj, jj - 1] - Cme[jj, jj - 1]
            rr = rho[kk[jj, jj]]
            Yv = (G[0][jj, jj] - E[0][jj, jj]) - (rr * mem.R(rr) / SQ2) * (G[2][jj, jj] - E[2][jj, jj])
            dfv = G[1][jj, jj] - E[1][jj, jj]
            print(f"   t = {tc:4.1f}: Y {Yv:+.2e}, df {dfv:+.2e}, dC++ {dCp:+.2e}, dC-- {dCm:+.2e}; |dC|/|Y| = {abs(dCp) / abs(Yv):.2e}, {abs(dCm) / abs(Yv):.2e}")


def tuned(c, rc, h, X):
    mem = ep.Member(c)
    dR, dL, info = ed.extract_mode(mem, rc, h, X + 20.0)
    N = int(round(X / h))
    par, st = ed.params(mem), ed.Static(mem, rc, N, h)

    def sign(A):
        d = ed.mode_rays(mem, rc, N, h, A, dR, dL)
        dp = ed.evolve_dep(d, par, st, Y0=3e-3)["dep"]
        return 0 if dp is None else dp["sign"]
    lo, hi = -1e-4, 1e-4
    fl = sign(lo)
    while hi - lo > 1e-14:
        mid = 0.5 * (lo + hi)
        if sign(mid) == fl:
            lo = mid
        else:
            hi = mid
    ps = 0.5 * (lo + hi)
    W = {-0.3: 0.9700800261, -1.5: 0.4941711675}.get(c)
    print(f"c = {c}, rho_c = {rc}, h = {h}: A*_h = {ps:.6e}, kappa = {mem.kappa:.4f}, #212 W = {W}")
    for dA in (0.0, 1e-9, -1e-9, 1e-7, -1e-7, 1e-5, -1e-5):
        A = ps + dA
        d = ed.mode_rays(mem, rc, N, h, A, dR, dL)
        r = ed.evolve_dep(d, par, st, Y0=1e9, stop_dep=False)
        y = r["ymax"]
        tt = (2 * np.arange(len(y)) + 1.0) * h / SQ2
        out = []
        for ta, tb in ((4, 7), (7, 10), (10, 13), (13, 16), (16, 19), (19, 22)):
            m = (tt >= ta) & (tt <= tb) & (y > 0)
            out.append(f"[{ta},{tb}] {np.polyfit(tt[m], np.log(y[m]), 1)[0]:+.3f} (|Y| {y[m][0]:.1e} to {y[m][-1]:.1e})")
        print(f"   A - A* = {dA:+.0e}: local rate of ln|Y| on the diagonal: " + "; ".join(out))


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "near":
        near(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))
    elif cmd == "resp":
        resp(float(sys.argv[2]), float(sys.argv[3]))
    elif cmd == "resp2":
        resp2(float(sys.argv[2]), float(sys.argv[3]))
    elif cmd == "constr":
        constr()
    elif cmd == "tuned":
        tuned(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))
    else:
        globals()[cmd]()

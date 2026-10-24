"""
Regenerate every figure of paper/main.tex from data committed in this repo.

    python paper/figures/make_figures.py        (run from the repo root or anywhere)

Inputs (all in the repo):
  phase2/dg_spectrum/growth_rates.csv          Fig. 3
  phase2/radial/crossing_{p,m}.txt, scan_eps±1.json, strong_{p,m}.txt   Fig. 4
  phase5/spherical/runs/*.csv                   Fig. 5
Figs. 1 and 2 are computed here from the closed forms (Dotti–Gleiser chi polynomial, eq. 55 of arXiv:0809.3615).
Units |M| = 1 (M = -1) for Figs. 1–3; code units m = 1 for Figs. 4–5.
"""
import json, os, re
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = HERE

plt.rcParams.update({
    # Computer Modern (matplotlib's bundled cmr10 + mathtext "cm") so figure text matches the LaTeX paper body.
    "font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm",
    "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
    "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 7.5,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "lines.linewidth": 1.3,
    "figure.dpi": 150, "savefig.bbox": "tight", "axes.spines.top": False, "axes.spines.right": False,
})
COL = ["#1f4e79", "#c0504d", "#4f8a3c", "#8064a2", "#d18b2c"]
W1 = 3.4  # single-column width (in)

r = sp.Symbol("r", positive=True)
M = sp.Integer(-1)


def chi_poly(l):
    terms = []
    for n in range(1, l + 3):
        c = sp.Integer(n - 2) * ((n - 4) * l * (l + 1) + n - 1) * sp.gamma(l + n - 1) \
            / (2**n * sp.gamma(n)**2 * sp.gamma(l - n + 3) * (-M)**n)
        terms.append(c * r**n)
    return sp.expand(sum(terms))


def closed_forms(l):
    """Return numpy callables V̂(r), κ_c(R), and the numbers r_s, r_*, r_neg."""
    c = (l + 2) * (l - 1)
    chi = chi_poly(l)
    f = 1 - 2 * M / r
    lam = sp.Rational(c, 2)
    V = 2 * f * (lam**2 * r**2 * ((lam + 1) * r + 3 * M) + 9 * M**2 * (lam * r + M)) / (r**3 * (lam * r + 3 * M)**2)
    psi0 = chi / (6 * M + r * c)
    g = f * sp.diff(psi0, r) / psi0
    Vhat = sp.simplify(V - 2 * f * sp.diff(g, r))
    Ncrit = sp.expand(c * chi - (c * r - 6) * sp.diff(chi, r))
    kc = (1 + 2 / r) * Ncrit / ((c * r - 6) * chi)
    rs = 6.0 / c
    rstar = [float(z) for z in sp.Poly(Ncrit, r).real_roots() if z > rs][0]
    Vf = sp.lambdify(r, Vhat, "numpy")
    rr = np.geomspace(1e-3, 50, 200001)
    vv = Vf(rr)
    rneg = rr[vv < 0].max()
    return Vf, sp.lambdify(r, kc, "numpy"), rs, rstar, rneg


def fig1_vhat():
    fig, ax = plt.subplots(figsize=(W1, 2.4))
    for i, l in enumerate([2, 3, 4]):
        Vf, _, rs, rstar, rneg = closed_forms(l)
        c = (l + 2) * (l - 1)
        rho = np.linspace(0.6, 30, 3000)          # scaled radius rho = c r
        rr = rho / c
        k = c * (c + 2) / 6
        v = Vf(rr) / k**2                          # DG eigenvalue scale: E0 = -(k/2)^2
        ax.plot(rho, v, color=COL[i], label=fr"$\ell={l}$")
        neg = v < 0
        ax.fill_between(rho, v, 0, where=neg, color=COL[i], alpha=0.18, lw=0)
    ax.axhline(0, color="k", lw=0.6)
    ax.axvline(6, color="0.4", ls=":", lw=0.9)
    ax.text(6.2, 1.2, r"$r_s$", color="0.3", fontsize=8)
    ax.set_xlim(2, 22)
    ax.set_ylim(-1.6, 1.5)
    ax.set_xlabel(r"$\rho = (\ell+2)(\ell-1)\,r/|M|$")
    ax.set_ylabel(r"$\hat V\,|M|^2/k^2$")
    ax.legend(frameon=False, loc="upper right")
    fig.savefig(os.path.join(OUT, "fig1_vhat.pdf"))
    plt.close(fig)


def fig2_kappa():
    fig, ax = plt.subplots(figsize=(W1, 2.5))
    rows = []
    for i, l in enumerate([2, 3, 4, 5]):
        _, kcf, rs, rstar, rneg = closed_forms(l)
        rows.append((l, rs, rstar, rneg))
        for lo, hi in [(0.02, rs * 0.999), (rs * 1.001, 30)]:
            R = np.geomspace(lo, hi, 2000)
            y = kcf(R) * R
            ax.plot(R, y, color=COL[i], label=fr"$\ell={l}$" if lo == 0.02 else None)
        ax.plot([rstar], [0], "o", ms=3, color=COL[i])
    ax.axhline(0, color="k", lw=0.6)
    ax.fill_between([0.02, 30], 0, 10, color="0.85", alpha=0.5, lw=0)
    ax.text(0.03, 7.5, r"$\kappa\geq0$", fontsize=8, color="0.3")
    ax.set_xscale("log")
    ax.set_xlim(0.02, 30)
    ax.set_ylim(-10, 10)
    ax.set_xlabel(r"$R/|M|$")
    ax.set_ylabel(r"$\kappa_c(R)\,R$")
    ax.legend(frameon=False, loc="upper right", ncol=2)
    fig.savefig(os.path.join(OUT, "fig2_kappa_c.pdf"))
    plt.close(fig)
    with open(os.path.join(OUT, "thresholds.csv"), "w") as fh:
        fh.write("l,r_s,r_star,r_neg\n")
        for row in rows:
            fh.write("%d,%.6f,%.6f,%.6f\n" % row)


def fig3_growth():
    import csv
    data = list(csv.DictReader(open(os.path.join(ROOT, "phase2/dg_spectrum/growth_rates.csv"))))
    l = np.array([int(d["l"]) for d in data])
    series = [
        ("full_DG_k_over_2", r"point mass ($R\to0$)", "k", "-"),
        ("R=0.3_Neumann", r"$R=0.3|M|$, Neumann", COL[0], "-o"),
        ("R=1.0_Neumann", r"$R=|M|$, Neumann", COL[1], "-s"),
        ("R=10.0_kappa=-10", r"$R=10|M|$, $\kappa=-10/|M|$", COL[2], "-^"),
    ]
    fig, ax = plt.subplots(figsize=(W1, 2.4))
    for key, lab, c, st in series:
        y = np.array([float(d[key]) for d in data])
        m = y > 0
        ax.plot(l[m], y[m], st, color=c, ms=3, label=lab)
    ax.set_yscale("log")
    ax.set_xlabel(r"$\ell$")
    ax.set_ylabel(r"growth rate $\gamma\,|M|$")
    ax.legend(frameon=False, loc="upper left")
    fig.savefig(os.path.join(OUT, "fig3_growth.pdf"))
    plt.close(fig)


def _read_crossing(fn):
    rows = []
    for line in open(os.path.join(ROOT, "phase2/radial", fn)):
        p = line.split()
        if len(p) == 4 and re.match(r"^[\d.]+$", p[0]):
            rows.append([float(x) for x in p])
    a = np.array(rows)
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3]


def _read_strong(fn):
    out = []
    for line in open(os.path.join(ROOT, "phase2/radial", fn)):
        w = float(re.search(r"w=([\d.]+)", line).group(1))
        N = float(re.search(r"N=([\d.]+)", line).group(1))
        Mv = float(re.search(r"M=([+\-\d.]+)", line).group(1))
        s = re.search(r"sigma2=\['([^']+)'", line)
        out.append((w, N, Mv, float(s.group(1)) if s else np.nan))
    return np.array(out)


def fig4_radial():
    fig, axs = plt.subplots(1, 2, figsize=(2 * W1, 2.4))
    ax = axs[0]
    for eps, fn, js, c in [(+1, "crossing_p.txt", "scan_eps+1.json", COL[0]),
                           (-1, "crossing_m.txt", "scan_eps-1.json", COL[1])]:
        sc = json.load(open(os.path.join(ROOT, "phase2/radial", js)))
        w = np.array([d["w"] for d in sc]); N = np.array([d["N"] for d in sc]); Mv = np.array([d["M"] for d in sc])
        o = np.argsort(w)
        lab = r"$\varepsilon=+1$" if eps > 0 else r"$\varepsilon=-1$"
        ax.plot(w[o], N[o], "-", color=c, label=lab + r": $N$")
        ax.plot(w[o], np.abs(Mv[o]), "--", color=c, label=lab + r": $|M|$")
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel(r"$N,\ |M|$")
    ax.legend(frameon=False, fontsize=7)
    ax = axs[1]
    for eps, fn, c in [(+1, "crossing_p.txt", COL[0]), (-1, "crossing_m.txt", COL[1])]:
        w, Mv, N, s2 = _read_crossing(fn)
        ax.plot(w, s2 * 1e3, "o-", ms=2.5, color=c, label=r"$\varepsilon=%+d$" % eps)
        dN = np.gradient(N, w)
        w0 = np.interp(0, dN, w) if dN[0] < dN[-1] else np.interp(0, dN[::-1], w[::-1])
        ax.axvline(w0, color=c, ls=":", lw=0.9)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel(r"$10^3\,\sigma_0^2$")
    ax.legend(frameon=False)
    ax.text(0.03, 0.06, r"dotted: $dN/d\omega=0$", transform=ax.transAxes, fontsize=7, color="0.3")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig4_radial.pdf"))
    plt.close(fig)


def _run(name):
    fn = os.path.join(ROOT, "phase5/spherical/runs", name + ".csv")
    skip = sum(1 for line in open(fn) if line.startswith("#"))
    return np.genfromtxt(fn, delimiter=",", names=True, skip_header=skip)


def fig5_evolution():
    fig, axs = plt.subplots(1, 2, figsize=(2 * W1, 2.4))
    ax = axs[0]
    for name, lab, c, st in [("c750_h0.2", r"$\varepsilon=+1$, $h=0.2$", COL[0], ":"),
                             ("c750_h0.1", r"$\varepsilon=+1$, $h=0.1$", COL[0], "-"),
                             ("p7603_h0.2", r"$\varepsilon=-1$, $h=0.2$", COL[1], ":"),
                             ("p7603_h0.1", r"$\varepsilon=-1$, $h=0.1$", COL[1], "-")]:
        d = _run(name)
        rel = np.abs(d["m50"] / d["m50"][0] - 1)
        ax.semilogy(d["t"], np.maximum(rel, 1e-12), st, color=c, label=lab)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$|m(t,50)/m(0,50)-1|$")
    ax.set_ylim(1e-11, 1e-4)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax = axs[1]
    for star, sig, lab, c in [("c940", 0.07370, r"$\varepsilon=+1$", COL[0]), ("p940", 0.06668, r"$\varepsilon=-1$", COL[1])]:
        dp = _run(star + "_h0.05_pp"); dm = _run(star + "_h0.05_pm")
        n = min(len(dp), len(dm))
        odd = 0.5 * (dp["amp"][:n] - dm["amp"][:n])
        t = dp["t"][:n]
        ax.semilogy(t, np.abs(odd), color=c, label=lab + r" ($\omega=0.94$)")
        tt = np.linspace(20, 150, 50)
        i0 = np.argmin(np.abs(t - 80.0))          # anchor inside the linear-growth phase
        ax.semilogy(tt, np.abs(odd[i0]) * np.exp(sig * (tt - t[i0])), "--", color=c, lw=0.9)
    ax.set_xlabel(r"$t$")
    ax.set_xlim(0, 250)
    ax.set_ylabel(r"$|\delta\phi_c|$ (odd part)")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.text(0.97, 0.05, r"dashed: $e^{\sigma t}$, $\sigma$ from radial analysis", transform=ax.transAxes,
            fontsize=6.5, color="0.3", ha="right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig5_evolution.pdf"))
    plt.close(fig)


def fig6_wormhole():
    """Paper III. (a) W^2 a^2 - 2(D - 3) against 1/(D - 3): Torii-Shinkai's Table I (D <= 11) and the r-integration
    (D >= 20), with the two- and three-term large-D results. (b) The growing mode's profile delta r / delta r(a) at
    D = 10, 20, 100, 1000 against the large-D limit (r/a) e^{1 - r^2/a^2}. Data: notes/wormhole_largeD/growth.csv,
    profiles.csv (written by notes/wormhole_largeD/profiles.py)."""
    import csv
    d = os.path.join(ROOT, "notes", "wormhole_largeD")
    rows = list(csv.DictReader(open(os.path.join(d, "growth.csv"))))
    D = np.array([int(r["D"]) for r in rows], float)
    W2 = np.array([float(r["W2a2"]) for r in rows])
    ts = np.array([r["source"].startswith("TS13") for r in rows])
    x, y = 1 / (D - 3), W2 - 2 * (D - 3)
    c0 = -4 * np.log(2)
    c1 = np.pi**2 / 3 + 4 * np.log(2) ** 2
    BLUE, ORANGE, AQUA, INK2 = "#2a78d6", "#eb6834", "#1baf7a", "#52514e"
    fig, axs = plt.subplots(1, 2, figsize=(2 * W1, 2.5))
    ax = axs[0]
    xx = np.linspace(0, 0.62, 200)
    ax.axhline(c0, color=INK2, ls="--", lw=1.0)
    ax.plot(xx, c0 + c1 * xx, color=AQUA, lw=1.4)
    ax.plot(x[ts], y[ts], "o", ms=4.5, mfc="white", mec=ORANGE, mew=1.2, label="Torii & Shinkai (2013)")
    ax.plot(x[~ts], y[~ts], "o", ms=3.5, color=BLUE, label="this work")
    for Dl in (4, 5, 6, 8, 11, 20):
        i = int(np.where(D == Dl)[0][0])
        ax.annotate(f"$D={Dl}$", (x[i], y[i]), xytext=(4, -9 if Dl > 5 else 3), textcoords="offset points",
                    fontsize=6.5, color="0.25")
    ax.text(0.97, c0 - 0.09, r"$-4\ln 2$", transform=ax.get_yaxis_transform(), ha="right", va="top", fontsize=7,
            color="0.25")
    ax.text(0.45, -0.66, r"$-4\ln 2+\dfrac{\pi^2/3+4\ln^2 2}{D-3}$", ha="left", va="center", fontsize=7,
            color="0.25")
    ax.set_xlim(0, 1.05)
    ax.set_ylim(-3.0, -0.4)
    ax.set_xlabel(r"$1/(D-3)$")
    ax.set_ylabel(r"$W^2a^2-2(D-3)$")
    ax.legend(frameon=False, loc="center right", bbox_to_anchor=(1.0, 0.36), fontsize=7)
    ax.set_title("(a)", loc="left", fontsize=8)

    ax = axs[1]
    prow = list(csv.DictReader(open(os.path.join(d, "profiles.csv"))))
    r = np.array([float(q["r_over_a"]) for q in prow])
    ramp = {10: "#86b6ef", 20: "#3987e5", 100: "#1c5cab", 1000: "#0d366b"}
    for Dv, col in ramp.items():
        Y = np.array([float(q[f"Y_D{Dv}"]) for q in prow])
        ax.plot(r, Y, color=col, lw=1.3, label=f"$D={Dv}$")
        if Dv in (10, 20):                       # 100 and 1000 nearly coincide with the limit; the legend names them
            j = np.argmin(np.abs(r - 2.05))
            ax.annotate(f"$D={Dv}$", (r[j], Y[j]), xytext=(2, 3), textcoords="offset points", fontsize=6.5,
                        color="0.25")
    lim = np.array([float(q["Y_largeD_limit"]) for q in prow])
    ax.plot(r, lim, "--", color=ORANGE, lw=1.2, label=r"$D\to\infty$: $(r/a)\,e^{1-r^2/a^2}$")
    ax.set_xlim(1, 3)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel(r"$r/a$")
    ax.set_ylabel(r"$\delta r(r)\,/\,\delta r(a)$")
    ax.legend(frameon=False, loc="upper right", fontsize=7)
    ax.set_title("(b)", loc="left", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig6_wormhole_largeD.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    import sys
    figs = {"fig1": fig1_vhat, "fig2": fig2_kappa, "fig3": fig3_growth, "fig4": fig4_radial, "fig5": fig5_evolution,
            "fig6": fig6_wormhole}
    for name in (sys.argv[1:] or figs):          # all figures, or only those named (e.g. fig6)
        figs[name](); print(name)

#!/usr/bin/env python3
"""Figures of Paper IV, from the outputs of record in notes/black_universe (read through sources.py).

    python make_figures.py        writes fig_rates.pdf and fig_slopes.pdf beside this file

fig_rates.pdf   #212's growth rate W and surface gravity kappa against c (stability_output.txt), W/kappa, and W/kappa - 1 near
                c = -1; the c = -1 member is marked by the rate of the linear double-null solver (#232, linear_output.txt R3).
fig_slopes.pdf  #217's lifetime slope against the grid spacing h, relative to 1/W (c = -0.3, -1.5) and 1/kappa (c = -1), with the
                OLS standard error (evolve_output.txt, analyze).
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sources  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
BLUE, ORANGE, AQUA, INK, GRID = "#257de0", "#f47230", "#00b894", "#222222", "#d9d9d9"   # palette validated (dataviz, light)
plt.rcParams.update({"font.family": "serif", "font.size": 8, "mathtext.fontset": "cm", "axes.edgecolor": INK,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "legend.frameon": False, "pdf.fonttype": 42})


def style(ax):
    ax.grid(True, color=GRID, linewidth=0.4)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def fig_rates():
    st = sources.stability()
    lin = sources.linear()
    rows = [(float(r["c"]), float(r["kappa"]), None if r["W"] is None else float(r["W"])) for r in st["rows"]]
    rate1, kap1 = float(lin["R3"][0]), float(lin["R3"][1])            # c = -1: linear double-null rate and kappa
    mode = [(c, k, w) for c, k, w in rows if w is not None]
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.45), constrained_layout=True)

    ax = axs[0]
    ax.plot([-c for c, k, w in mode], [w for c, k, w in mode], "-o", color=BLUE, ms=3.2, lw=1.2, label="$W$")
    ax.plot([-c for c, k, w in rows], [k for c, k, w in rows], "--s", color=ORANGE, ms=2.6, lw=1.0, label=r"$\kappa$")
    ax.plot([1.0], [rate1], "*", color=AQUA, ms=8, mec=INK, mew=0.4, label=r"$c=-1$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$-c$")
    ax.set_ylabel(r"rate ($1/b$)")
    ax.set_ylim(0, 1.3)
    ax.legend(loc="upper right", fontsize=7, handlelength=1.8)
    ax.set_title("(a) growth rate and surface gravity", fontsize=8, loc="left")
    style(ax)

    ax = axs[1]
    ax.plot([-c for c, k, w in mode], [w / k for c, k, w in mode], "-o", color=BLUE, ms=3.2, lw=1.2, label=r"$W/\kappa$")
    ax.plot([1.0], [rate1 / kap1], "*", color=AQUA, ms=8, mec=INK, mew=0.4, label=r"$c=-1$")
    ax.axhline(1.0, color=INK, lw=0.6, ls=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_yticks([0.8, 1, 2, 4, 8])
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(r"$-c$")
    ax.set_ylabel(r"$W/\kappa$")
    ax.annotate(r"$c=-1$: rate $\kappa$" + "\n(horizon sector)", xy=(1.0, rate1 / kap1), xytext=(1.6, 2.6), fontsize=7,
                color=INK, arrowprops=dict(arrowstyle="-", color=INK, lw=0.5))
    ax.set_title(r"(b) $W/\kappa$ against $c$", fontsize=8, loc="left")
    style(ax)

    ax = axs[2]
    near = [(c, k, w) for c, k, w in mode if abs(c + 1) <= 0.11]
    ax.plot([c + 1 for c, k, w in near], [w / k - 1 for c, k, w in near], "o", color=BLUE, ms=3.2, label=r"$W/\kappa-1$")
    ax.plot([0.0], [rate1 / kap1 - 1], "*", color=AQUA, ms=8, mec=INK, mew=0.4)
    ax.axhline(0, color=INK, lw=0.6, ls=":")
    ax.axvline(0, color=INK, lw=0.6, ls=":")
    ax.set_xlabel(r"$c+1$")
    ax.set_ylabel(r"$W/\kappa-1$")
    ax.set_title(r"(c) near $c=-1$", fontsize=8, loc="left")
    style(ax)
    path = os.path.join(HERE, "fig_rates.pdf")
    fig.savefig(path)
    print("wrote", path)


def fig_slopes():
    ev = sources.evolve()
    fig, ax = plt.subplots(figsize=(4.6, 2.6), constrained_layout=True)
    marks = {-0.3: ("o", BLUE, r"$c=-0.3$, against $1/W$"), -1.5: ("s", ORANGE, r"$c=-1.5$, against $1/W$"),
             -1.0: ("D", AQUA, r"$c=-1$, against $1/\kappa$")}
    shift = {-0.3: 0.93, -1.5: 1.0, -1.0: 1.07}
    for c, (mk, col, lab) in marks.items():
        hs = sorted({h for (cc, h) in ev["runs"] if cc == c}, key=float)
        x = [float(h) * shift[c] for h in hs]
        ref = [float(ev["runs"][(c, h)]["ref"]) for h in hs]
        y = [float(ev["runs"][(c, h)]["slope"]) / r - 1 for h, r in zip(hs, ref)]
        e = [float(ev["runs"][(c, h)]["se"]) / r for h, r in zip(hs, ref)]
        ax.errorbar(x, y, yerr=e, fmt=mk + "-", color=col, ms=3.5, lw=1.0, elinewidth=0.8, capsize=2, label=lab)
    ax.axhline(0, color=INK, lw=0.6, ls=":")
    ax.set_xscale("log")
    ax.set_xticks([0.0025, 0.005, 0.01, 0.02])
    ax.set_xticklabels(["0.0025", "0.005", "0.01", "0.02"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(r"grid spacing $h$")
    ax.set_ylabel("slope / reference $-$ 1")
    ax.legend(loc="upper left", fontsize=7)
    style(ax)
    path = os.path.join(HERE, "fig_slopes.pdf")
    fig.savefig(path)
    print("wrote", path)


if __name__ == "__main__":
    fig_rates()
    fig_slopes()

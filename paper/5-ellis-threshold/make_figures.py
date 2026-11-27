#!/usr/bin/env python3
"""Figures of Paper V (paper/5-ellis-threshold), drawn from the runs of record in negative-mass-raw (read only).

    python3 make_figures.py [--raw=~/Source/negative-mass-raw]

Writes fig_domain.pdf, fig_lifetime.pdf, fig_locus.pdf, fig_mass.pdf and fig_expansion.pdf beside this script and
prints the plotted numbers as tables. Inputs:
  I40-dh610/runs/D*_S_h*/points.jsonl, summary.json     family S series of #40 (threshold.py)
  I225-dh610/runs/series/D*_G_h*/points_*.jsonl         family G series of #225 (throat_locus.py), diagonal and locus
  I41-rpi5/runs/h0.0025_d1e-08/series.npz, scan_h0.01   M_AH(x+) and the delta scan of S (#41, settled_mass.py)
  I253-rpi5/runs/h0.005_d1e-0{6,8}/series.npz, scan_G   the same for G (#253, settled_mass_g.py)
  I213-dh610/runs/D{4,5}_tuned_*.npz                     the throat diagonal of the expanding runs (#213, expansion.py)
  notes/wormhole_evolution/expansion_output.txt         tau_s of those runs
The slopes, masses and k quoted in main.tex come from the notes' output files through paper_numbers.py, not from here.
"""
import argparse
import glob
import json
import math
import os
import re

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, "..", "..", "notes", "wormhole_evolution")
ap = argparse.ArgumentParser()
ap.add_argument("--raw", default=os.path.expanduser("~/Source/negative-mass-raw"))
RAW = ap.parse_args().raw

W2 = {4: 1.39705243433350684, 5: 2.98495893035176025, 6: 4.68662054311946479}   # #195's crate, a = 1
W = {D: math.sqrt(v) for D, v in W2.items()}
FINE = {4: 0.005, 5: 0.00176777, 6: 0.00125}

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948", "#a6339a", "#7a7668"]
INK, INK_SOFT, GRID = "#0b0b0b", "#52514e", "#e7e6e2"
GRAY = SERIES[-1]
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"], "mathtext.fontset": "stix",
    "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.edgecolor": INK_SOFT, "axes.labelcolor": INK, "xtick.color": INK_SOFT, "ytick.color": INK_SOFT,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "legend.frameon": False, "pdf.fonttype": 42, "lines.linewidth": 1.4,
})


def save(fig, name):
    path = os.path.join(HERE, name)
    fig.savefig(path)
    plt.close(fig)
    print(f"wrote {path}")


def jsonl(path):
    return [json.loads(line) for line in open(path) if line.strip()]


# ---------------------------------------------------------------------------------------------------------------------
# 1. the domain of the double-null scheme
# ---------------------------------------------------------------------------------------------------------------------
def fig_domain():
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.9), gridspec_kw={"width_ratios": [1, 1.25]})
    for a in ax:
        a.set_aspect("equal")
        a.axis("off")
    a = ax[0]
    a.add_patch(plt.Rectangle((0, 0), 1, 1, fc="#f4f4f2", ec=INK_SOFT, lw=1))
    a.plot([0, 1], [0, 0], color=SERIES[0], lw=2.5)
    a.plot([0, 0], [0, 1], color=SERIES[0], lw=2.5)
    a.plot([0, 1], [0, 1], color=SERIES[7], lw=1.3, ls=(0, (4, 3)))
    a.plot(0, 0, "o", color=SERIES[7], ms=5)
    a.annotate("", xy=(0.62, 0.42), xytext=(0.62, 0.03), arrowprops=dict(arrowstyle="-|>", color=SERIES[3], lw=1.2))
    a.add_patch(matplotlib.patches.Ellipse((0.62, 0.0), 0.12, 0.05, fc=SERIES[3], ec="none", alpha=0.9))
    a.text(0.5, -0.09, r"$\Sigma^+$: $x^-=0$ (right universe)", ha="center", va="top", color=SERIES[0])
    a.text(-0.06, 0.5, r"$\Sigma^-$: $x^+=0$ (left universe)", ha="right", va="center", rotation=90, color=SERIES[0])
    a.text(0.47, 0.57, "throat $x^+=x^-$", rotation=45, ha="center", va="bottom", color=SERIES[7])
    a.text(0.66, 0.24, "pulse", color="#8a5a00", ha="left", va="center")
    a.text(-0.02, -0.03, "S", ha="right", va="top", color=SERIES[7])
    a.text(1.03, 0.5, r"edge $x^+=X$", rotation=90, ha="left", va="center", color=INK_SOFT)
    a.text(0.5, 1.03, r"edge $x^-=X$", ha="center", va="bottom", color=INK_SOFT)
    a.text(0.25, 0.80, "no outer boundary:\nthe domain is the\nfuture of S", ha="center", va="center", color=INK_SOFT,
           fontsize=8)
    a.set_xlim(-0.25, 1.15)
    a.set_ylim(-0.22, 1.12)
    a.set_title(r"(a) the rectangle in $(x^+, x^-)$", fontsize=9, loc="left")
    b = ax[1]
    s = 1 / math.sqrt(2)
    P = np.array([[0, 0], [s, s], [0, 2 * s], [-s, s], [0, 0]])
    b.add_patch(plt.Polygon(P[:4], fc="#f4f4f2", ec="none"))
    b.plot(P[:2, 0], P[:2, 1], color=SERIES[0], lw=2.5)
    b.plot([0, -s], [0, s], color=SERIES[0], lw=2.5)
    b.plot([s, 0, -s], [s, 2 * s, s], color=INK_SOFT, lw=1)
    b.plot([0, 0], [0, 2 * s], color=SERIES[7], lw=1.3, ls=(0, (4, 3)))
    b.plot(0, 0, "o", color=SERIES[7], ms=5)
    b.text(0.04, 0.62, "throat\nworldline,\n$t\\leq\\sqrt{2}\\,X$", color=SERIES[7], fontsize=8, va="center")
    b.text(0.45, 0.2, r"$\Sigma^+$, $l>0$", color=SERIES[0], ha="left")
    b.text(-0.45, 0.2, r"$\Sigma^-$, $l<0$", color=SERIES[0], ha="right")
    b.annotate("", xy=(0.95, 0.0), xytext=(-0.95, 0.0), arrowprops=dict(arrowstyle="-|>", color=INK_SOFT, lw=0.8))
    b.annotate("", xy=(-0.95, 1.5), xytext=(-0.95, 0.0), arrowprops=dict(arrowstyle="-|>", color=INK_SOFT, lw=0.8))
    b.text(0.97, 0.0, "$l$", va="center", color=INK_SOFT)
    b.text(-0.95, 1.53, "$t$", ha="center", va="bottom", color=INK_SOFT)
    b.set_xlim(-1.05, 1.05)
    b.set_ylim(-0.12, 1.62)
    b.set_title(r"(b) the same domain in $l=(x^+-x^-)/\sqrt{2}$, $t=(x^++x^-)/\sqrt{2}$", fontsize=9, loc="left")
    fig.subplots_adjust(left=0.02, right=0.99, top=0.9, bottom=0.03, wspace=0.08)
    save(fig, "fig_domain.pdf")


# ---------------------------------------------------------------------------------------------------------------------
# 2. the lifetime law, family S, finest h at each D
# ---------------------------------------------------------------------------------------------------------------------
def fig_lifetime():
    fig, ax = plt.subplots(1, 3, figsize=(6.5, 2.45), sharey=False)
    print("\nlifetime law, family S, finest h: per side, points with e <= 11, intercept at slope 1/W, rms about the line")
    for a, D in zip(ax, (4, 5, 6)):
        d = os.path.join(RAW, "I40-dh610", "runs", f"D{D}_S_h{FINE[D]:g}")
        sm = json.load(open(os.path.join(d, "summary.json")))
        pstar = 0.5 * (sm["lo"] + sm["hi"])
        pts = [p for p in jsonl(os.path.join(d, "points.jsonl")) if p.get("phase") == "series" and p["fate"] == p["side"]]
        for side, col, lab in ((-1, SERIES[0], "collapse"), (1, SERIES[1], "expansion")):
            fit = [p for p in pts if p["side"] == side and p["e"] <= 11]
            flo = [p for p in pts if p["side"] == side and p["e"] > 11]
            x = np.array([math.log(1 / abs(p["A"] - pstar)) for p in fit])
            T = np.array([p["T"] for p in fit])
            c = float(np.mean(T - x / W[D]))
            xs = np.array([x.min() - 1, x.max() + 1])
            a.plot(xs, c + xs / W[D], color=col, lw=1, alpha=0.8)
            a.plot(x, T, "o" if side < 0 else "s", ms=3.5, mfc=col, mec="white", mew=0.5, color=col, ls="none",
                   label=lab)
            if flo:
                xf = [math.log(1 / abs(p["A"] - pstar)) for p in flo]
                a.plot(xf, [p["T"] for p in flo], "o" if side < 0 else "s", ms=3.5, mfc="none", mec=GRAY, ls="none",
                       label="roundoff floor (not fitted)" if side < 0 else None)
            rms = float(np.sqrt(np.mean((T - c - x / W[D]) ** 2)))
            print(f"  D = {D} h = {FINE[D]:g} {lab:9s}: {len(x)} points, x in [{x.min():.2f}, {x.max():.2f}], "
                  f"c = {c:.4f}, rms {rms:.4f}")
        a.set_title(f"$D={D}$, $h={FINE[D]:g}$" if D != 5 else r"$D=5$, $h=0.0025/\sqrt{2}$", fontsize=9, loc="left")
        a.text(0.04, 0.93, f"lines: slope $1/W={1 / W[D]:.4f}\\,a$", transform=a.transAxes, fontsize=7.5,
               color=INK_SOFT, va="top")
        a.set_xlabel(r"$\ln(1/|p-p_{*h}|)$")
    ax[0].set_ylabel(r"lifetime $t$ of the event $(a)$")
    ax[0].legend(loc="lower right", handletextpad=0.3, borderaxespad=0.2, fontsize=7)
    fig.tight_layout(pad=0.4, w_pad=0.8)
    save(fig, "fig_lifetime.pdf")


# ---------------------------------------------------------------------------------------------------------------------
# 3. the one-sided family G: the diagonal against the throat's locus
# ---------------------------------------------------------------------------------------------------------------------
def fig_locus():
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.6))
    hh = {4: "0.005", 5: "0.00176777"}
    print("\nfamily G: lifetime on the diagonal (t, #40's events) and on the locus (tau at r_min = 0.5 or 2)")
    for a, D in zip(ax, (4, 5)):
        pts = []
        for f in sorted(glob.glob(os.path.join(RAW, "I225-dh610", "runs", "series", f"D{D}_G_h{hh[D]}", "points_*.jsonl"))):
            pts += [p for p in jsonl(f) if p.get("phase") == "series"]
        pts = [p for p in pts if p["e"] is not None and p["e"] <= 11 and p["fate40"] == p["side"] == p["fate_loc"]]
        for side, mk in ((-1, "o"), (1, "s")):
            q = sorted((p for p in pts if p["side"] == side), key=lambda p: p["e"])
            x = np.array([p["e"] * math.log(10) for p in q])
            Td = np.array([p["T40"] for p in q])
            key = "c:0.5" if side < 0 else "e:2.0"
            Tl = np.array([p["loc"][key][1] for p in q])
            for T, col, lab in ((Td, GRAY, "diagonal, $t$"), (Tl, SERIES[0] if side < 0 else SERIES[1], "locus, $\\tau$")):
                s, c = np.polyfit(x, T, 1)
                a.plot(x, T, mk, ms=3.5, mfc=col, mec="white", mew=0.5, ls="none",
                       label=f"{lab}, {'collapse' if side < 0 else 'expansion'}")
                print(f"  D = {D} side {side:+d} {lab:14s}: {len(x)} points, own-side least-squares slope {s:.4f}")
            xs = np.array([x.min(), x.max()])
            c0 = float(np.mean(Tl - x / W[D]))
            a.plot(xs, c0 + xs / W[D], color=INK, lw=0.8, ls=(0, (3, 2)))
        a.set_title(f"$D={D}$, $h={hh[D]}$" if D == 4 else r"$D=5$, $h=0.0025/\sqrt{2}$", fontsize=9, loc="left")
        a.set_xlabel(r"$\ln(1/|p-p_{*h}|)$")
        a.text(0.04, 0.93, f"dashed: slope $1/W={1 / W[D]:.4f}\\,a$", transform=a.transAxes, fontsize=7.5,
               color=INK_SOFT, va="top")
    ax[0].set_ylabel(r"lifetime $(a)$")
    ax[1].legend(loc="lower right", fontsize=7, handletextpad=0.3, borderaxespad=0.2)
    fig.tight_layout(pad=0.4, w_pad=1.0)
    save(fig, "fig_locus.pdf")


# ---------------------------------------------------------------------------------------------------------------------
# 4. the settled mass: M_AH(x+) for S and G, and the plateau against delta
# ---------------------------------------------------------------------------------------------------------------------
def fig_mass():
    def series(p):
        return dict(np.load(os.path.join(p, "series.npz")))

    S = series(os.path.join(RAW, "I41-rpi5", "runs", "h0.0025_d1e-08"))
    G8 = series(os.path.join(RAW, "I253-rpi5", "runs", "h0.005_d1e-08"))
    G6 = series(os.path.join(RAW, "I253-rpi5", "runs", "h0.005_d1e-06"))
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.7), gridspec_kw={"width_ratios": [1.35, 1]})
    a = ax[0]
    print("\nsettled mass: M_AH at x+ - x+_f = 0, 2.5, 5, end")
    for s, col, ls, lab in ((S, SERIES[1], "-", r"S, $h=0.0025$, $\delta=10^{-8}$"),
                            (G8, SERIES[0], "-", r"G, $h=0.005$, $\delta=10^{-8}$"),
                            (G6, SERIES[0], (0, (4, 2)), r"G, $h=0.005$, $\delta=10^{-6}$")):
        x = s["xp"] - s["xp"][0]
        a.plot(x, s["M"], color=col, ls=ls, label=lab)
        vals = [float(np.interp(v, x, s["M"])) for v in (0, 2.5, 5.0)] + [float(s["M"][-1])]
        print(f"  {lab:34s}: " + ", ".join(f"{v:.5f}" for v in vals) + f"; e^-f at the end {s['F'][-1]:.3g}")
    a.axhline(0.30, color=GRAY, lw=0.8, ls=(0, (1, 2)))
    a.text(0.15, 0.303, "0.30 [SH02]", ha="left", va="bottom", fontsize=7.5, color=INK_SOFT)
    a.axhline(0.225, color=GRAY, lw=0.8, ls=(0, (1, 2)))
    a.text(9.4, 0.221, r"$r_{\rm AH}/2=0.225$ [GGS09b]", ha="right", va="top", fontsize=7.5, color=INK_SOFT)
    a.set_xlabel(r"$x^+-x^+_f$ along the apparent horizon $(a)$")
    a.set_ylabel(r"$M_{\rm AH}=r_{\rm AH}/2$ $(a)$")
    a.set_xlim(0, 9.6)
    a.set_ylim(0.20, 0.51)
    a.legend(loc="upper right", borderaxespad=0.2, fontsize=7)
    ins = a.inset_axes([0.55, 0.30, 0.42, 0.30])
    for s, col, ls in ((S, SERIES[1], "-"), (G8, SERIES[0], "-"), (G6, SERIES[0], (0, (4, 2)))):
        x = s["xp"] - s["xp"][0]
        ins.plot(x, s["M"], color=col, ls=ls, lw=1.1)
    ins.set_xlim(4.5, 9.5)
    ins.set_ylim(0.2215, 0.2375)
    ins.tick_params(labelsize=6.5)
    ins.grid(True, lw=0.4)
    b = ax[1]
    s41 = [r for r in json.load(open(os.path.join(RAW, "I41-rpi5", "runs", "scan_h0.01.json"))) if "M_end" in r]
    s41 = sorted((r for r in s41 if r["delta"] <= 1e-3), key=lambda r: r["delta"])
    b.plot([r["delta"] for r in s41], [r["M_end"] for r in s41], "o-", color=SERIES[1], ms=3.5, mec="white", mew=0.5,
           label="S, $h=0.01$")
    print("\nplateau against delta (S: M_end, h = 0.01; G: plateau, h = 0.005)")
    print("  S: " + ", ".join(f"{r['delta']:.0e} {r['M_end']:.5f}" for r in s41))
    sg = [r for r in json.load(open(os.path.join(RAW, "I253-rpi5", "runs", "scan_G_h0.005.json"))) if "plat1" in r]
    sg = sorted((r for r in sg if r["delta"] <= 1e-3), key=lambda r: r["delta"])
    b.plot([r["delta"] for r in sg], [r["plat1"] for r in sg], "s-", color=SERIES[0], ms=3.5, mec="white", mew=0.5,
           label="G, $h=0.005$")
    print("  G: " + ", ".join(f"{r['delta']:.0e} {r['plat1']:.5f}" for r in sg))
    b.set_xscale("log")
    b.set_xlabel(r"distance from threshold $\delta=p_{*h}-p$")
    b.set_ylabel(r"settled $M_{\rm AH}$ $(a)$")
    b.legend(loc="center right", borderaxespad=0.2)
    fig.tight_layout(pad=0.4, w_pad=1.2)
    save(fig, "fig_mass.pdf")


# ---------------------------------------------------------------------------------------------------------------------
# 5. the expanding branch on the symmetric throat
# ---------------------------------------------------------------------------------------------------------------------
def fig_expansion():
    out = open(os.path.join(EV, "expansion_output.txt")).read()
    runs = ((4, "0.01000", "D4_tuned_h0.01.npz"), (5, "0.00354", "D5_tuned_h0.00353553.npz"))
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.6))
    print("\nexpanding branch, tuned members: H/W and k = H (tau_s - tau) at r/a = 10, 25, 50")
    for (D, htxt, fn), col in zip(runs, (SERIES[0], SERIES[1])):
        blk = out.split(f"--- D{D}, h = {htxt}, member tuned")[1]
        ts = float(re.search(r"1/H linear in tau: k = [0-9.]+, tau_s = ([0-9.]+)", blk).group(1))
        z = np.load(os.path.join(RAW, "I213-dh610", "runs", fn), allow_pickle=True)
        ok = z["alive_diag"] & (z["r"] >= 1.02)
        r, H, tau = z["r"][ok], z["H"][ok], z["tau"][ok]
        sel = np.unique(np.searchsorted(np.maximum.accumulate(r), np.geomspace(r[0], r[-1], 600)).clip(0, len(r) - 1))
        lab = f"$D={D}$, $\\tau_s={ts:.4f}\\,a$"
        ax[0].plot(r[sel], H[sel] / W[D], color=col, label=lab)
        ax[1].plot(r[sel], (H * (ts - tau))[sel], color=col, label=lab)
        vals = [(float(np.interp(v, r, H)) / W[D], float(np.interp(v, r, H * (ts - tau)))) for v in (10, 25, 50)]
        print(f"  D = {D}: tau_s {ts}; " + "; ".join(f"r/a = {v}: H/W {a_:.3f}, k {b_:.4f}" for v, (a_, b_) in
                                                     zip((10, 25, 50), vals)))
    a = ax[0]
    a.set_xscale("log")
    a.set_yscale("log")
    a.axhline(1, color=GRAY, lw=0.8, ls=(0, (1, 2)))
    a.text(1.1, 1.08, "$H=W$ (linear rate)", fontsize=7.5, color=INK_SOFT)
    a.set_xlabel(r"areal radius of the throat $r/a$")
    a.set_ylabel(r"$H/W$, $H=d\ln r/d\tau$")
    a.legend(loc="lower right", borderaxespad=0.2)
    b = ax[1]
    b.set_xscale("log")
    for kv in (2, 4):
        b.axhline(kv, color=GRAY, lw=0.8, ls=(0, (1, 2)))
    b.set_ylim(0, 5)
    b.set_xlabel(r"areal radius of the throat $r/a$")
    b.set_ylabel(r"$k=H\,(\tau_s-\tau)$")
    b.legend(loc="lower right", borderaxespad=0.2)
    fig.tight_layout(pad=0.4, w_pad=1.2)
    save(fig, "fig_expansion.pdf")


if __name__ == "__main__":
    fig_domain()
    fig_lifetime()
    fig_locus()
    fig_mass()
    fig_expansion()

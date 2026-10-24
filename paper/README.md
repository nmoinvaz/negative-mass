# Papers

The draft was split into two companion papers that share one figure script. Paper III, on ghost-scalar wormholes, is
separate and uses the same script (`python figures/make_figures.py fig6` makes only its figure).

| | Paper | Source | PDF |
|---|---|---|---|
| I | *The instability of negative-mass Schwarzschild is a point-mass effect: exact thresholds for finite stars* | `1-thresholds/main.tex` | [`1-thresholds/main.pdf`](1-thresholds/main.pdf) |
| II | *Stability of a negative-mass boson star* | `2-boson-star/main.tex` | [`2-boson-star/main.pdf`](2-boson-star/main.pdf) |
| III | *The instability of ghost-scalar wormholes at large dimension* | `3-wormhole-largeD/main.tex` | [`3-wormhole-largeD/main.pdf`](3-wormhole-largeD/main.pdf) |

Paper I is the mathematical result (Theorems 1–2, thresholds). Paper II applies it to Shirokov's phantom boson star and
cites Paper I for κ_c. The pre-split single draft is in git history (commit 2995273).

Build:
```
python figures/make_figures.py     # regenerates figures/*.pdf and figures/thresholds.csv from repo data
cd 1-thresholds && latexmk -pdf main.tex
cd ../2-boson-star && latexmk -pdf main.tex
cd ../3-wormhole-largeD && latexmk -pdf main.tex
```

Where each result comes from:

| paper | result | repo |
|---|---|---|
| I | Lemma 1, Theorems 1–2, Corollary 1, κ_c asymptotics, large-ℓ thresholds | `phase2/dg_proof/` |
| I | Spectral numerics, 1000-case truncated scan, growth rates (Fig. 3) | `phase2/dg_spectrum/` |
| I | Thresholds table, Figs. 1–2 | computed in `figures/make_figures.py` from the closed forms |
| II | Boson-star backgrounds, ADM vs. Shirokov's 5 %-radius mass | `phase2/radial/background.py`, `phase2/mass_check/` |
| II | Radial perturbation equations and spectra (Fig. 1) | `phase2/radial/` |
| II | Nonlinear 1+1 evolutions (Fig. 2) | `phase5/spherical/` |
| II | Non-radial (polar/axial) equations, growing-mode scans, κ_eff | `phase2/nonradial/` |
| both | Prior-work check | `notes/novelty/RESULTS.md` |
| III | Master equation, spectra, matched expansion, every quoted number, Fig. 6 data | `notes/wormhole_largeD/` (`paper_numbers.py`) |

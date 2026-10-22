# Paper

`main.tex` → `main.pdf`: *The instability of negative-mass Schwarzschild is a point-mass effect: exact thresholds for
finite stars and a radially stable negative-mass boson star.*

Build:
```
python figures/make_figures.py     # regenerates figures/*.pdf and figures/thresholds.csv from repo data
latexmk -pdf main.tex
```

Where each result comes from:

| paper | repo |
|---|---|
| Lemma 1, Theorems 1–2, Corollary 1, κ_c asymptotics, large-ℓ thresholds | `phase2/dg_proof/` |
| Spectral numerics, 1000-case truncated scan, growth rates (Fig. 3) | `phase2/dg_spectrum/` |
| Thresholds table, Figs. 1–2 | computed in `figures/make_figures.py` from the closed forms |
| Boson-star backgrounds, ADM vs. paper's 5 %-radius mass | `phase2/radial/background.py`, `phase2/mass_check/` |
| Radial perturbation equations and spectra (Fig. 4) | `phase2/radial/` |
| Nonlinear 1+1 evolutions (Fig. 5) | `phase5/spherical/` |
| Prior-work check | `notes/novelty/RESULTS.md` |

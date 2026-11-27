# A double-null evolution that passes the wormhole throat (#39)

Spherically symmetric Einstein gravity in D = n + 2 dimensions with any number of massless scalars of either sign, in double-null
coordinates, with both universes of a wormhole in the grid. Design and its rulings: [design.md](design.md). Equations, checked symbolically:
[derive.py](derive.py) (25 checks, 9 mutations; [derive_output.txt](derive_output.txt)). Code: [dn.py](dn.py). Checks:
[checks.py](checks.py). Raw output: negative-mass-raw `I39-dh610`, whose README holds the criteria (written before the runs of record), the
development runs and the runs of record. Status words as in AGENTS.md: *checked symbolically*, *numerical evidence*; "planted" means an exact
solution run through the same code.

## What the code does

- Metric ds² = −2e^{−f}dx⁺dx⁻ + r²dΩₙ², fields r, f, φₖ (εₖ = +1 normal, −1 ghost), 8πG = 1, a = 1. Evolved equations
  r₊₋ = −(n−1)(r₊r₋ + e^{−f}/2)/r, f₊₋ = (n/r)r₊₋ + Σεφ₊φ₋, φ₊₋ = −n(r₊φ₋ + r₋φ₊)/(2r). The constraints are monitored.
- The domain is the rectangle [x⁺₀, x⁺₀ + Nh] × [x⁻₀, x⁻₀ + Nh] with characteristic data on its two lower edges, no outer boundary.
  The corner of the Ellis data is the throat, so both universes are covered, and t = (x⁺ + x⁻)/√2 is the static time.
- One cell update (predictor + one corrector, second order), two kernels: `numba` (row sweep, production) and `numpy` (anti-diagonal sweep,
  reference; no numba needed: `DN_NO_NUMBA=1` emulates its absence). Cells with r < r_cut or e^{−f} > F_cut (0.05, 25) are excised, with their
  future.
- Data (`dn.py`): `ellis_rays` (static Ellis wormhole in any D, plus Gaussian pulses of any field on either ray, `kind='phi'` or the derivative
  form of [SH02] (3.1) `'dphi'`, optional mirror-symmetric pairs; the constraint is solved by RK4), `tangherlini_rays` (exact, Kruskal
  coordinates, horizon at x⁻ = −U₀ with r = (2M)^{1/(n−1)}), `roberts_rays` (the exact ghost solution of derive.py D5, D = 4).
- Diagnostics: trapped / anti-trapped classification and events (`cell_class`, `events_from_grid`; streaming when `stop=True`), the diagonal
  throat series r, f, H, τ, the minimum of r on the slices t = const (`throat_min`), zero contours of r₊ and r₋ with their slope (`contour`),
  the Misner–Sharp mass (`mass_nodes`), constraints C₊₊, C₋₋, the edge flux residual, T₊₊, T₋₋, T^θ_θ, ρ (`stress`).

```
python checks.py                      all checks, numba kernel (run from this directory; do not use python -I, it drops the script's directory)
DN_NO_NUMBA=1 python checks.py --kernel=numpy        the reference kernel with numba absent (K is skipped)
python checks.py --only=P1,P3 --mutate=theta_coef    break one input; the checks in MUTATIONS must fire
python probe_expansion.py             the expansion probe below (not a check)
```

## Checks and numbers

Run of record, `run01_checks` (numba, 4 min 30 s, one thread, load 18 of 20): **11 of 11 pass**. The numpy-kernel run, numba absent, passes 10 of 11
with K skipped and its log equals the numba log apart from timings and K (17 min 38 s). All 14 mutations fire exactly their declared checks.

| id | what | result (numerical evidence) |
|---|---|---|
| I1 | constraint solve returns the Ellis r(l), n = 2, 3, 4 | 2.9e-9, tolerance 1e-8 |
| K | numba against numpy on three grids (one collapsing, with excision) | difference 0 (bit-identical), same dead cells |
| P1 | static Ellis stays static to the truncation error, for t ≤ t_c (4.0 at D = 4, 2.7 at D = 5) | order 2.05, 1.98 (D = 4); 2.06, 2.00 (D = 5) |
| P2 | exact Tangherlini in Kruskal coordinates, D = 4, 5, 6: (r, f), M, the horizon position and r on it | (r, f) 1.97 to 2.01; horizon 1.998 to 2.000; M 1.77 to 1.95 (inside [1.75, 2.25], below 2) |
| P3 | exact Roberts ghost solution, D = 4, crossing r₊ = 0 and r₋ = 0 into the anti-trapped wedge | order 1.93, 1.97; contours 1.98 to 2.00; slopes 0.666668 (2/3), 1.499991 (3/2); 0 misclassified cells |
| C1, C2 | constraints and the mass-flux identity on a pulse run | order 1.99, 2.00 and 1.99, 2.00 |
| G1 | **growth rate**: tune p to p*ₕ, fit the departure of H | D = 4: 1.181988, 1.181934 against W = 1.181970 (+1.5e-5, −3.0e-5), h = 0.02, 0.01; D = 5: 1.727492, 1.727247 against 1.727703 (−1.2e-4, −2.6e-4), h = 0.00707, 0.00354. Bar 1 %, agreement 0.5 % |
| F1 | D = 4, one pulse of each sign | ghost with E₀ = −0.030 expands; ghost with E₀ = +0.025 collapses (first trapped cell at M = 0.4975); normal pulse, E₀ = +0.0022, collapses. [SH02]'s two fates and its E₀ sign rule |
| S | excision r_cut halved and F_cut doubled | same fate and event time (difference 0) |

W is √(W²a²) of #195 (`notes/wormhole_largeD/w2_output.txt`): 1/(Wa) = 0.846045 (D = 4) and 0.578803 (D = 5), the issue's controls.

**The seed.** Truncation error seeds the unstable mode: the tuned threshold p*ₕ scales as h² (ratio 3.999 at D = 4 and 4.002 at D = 5 per halving), and an
untuned run departs after (1/W) ln(1/(κh²)) with κ = 0.0042 (D = 4) and 0.084 (D = 5), the same at both resolutions (the departure time shifts by
1.181 against ln4/W = 1.173 at D = 4, and by 0.802 against 0.802 at D = 5; departure means |H| ≥ 2e-3 W on the diagonal).

**Fit.** ln|H| against τ over the points after the global minimum of |H| with |H|/W ∈ [1e-3, 1e-2], the two sides of p*ₕ averaged. Each side alone is
off by up to ±0.5 % (nonlinear bias, opposite signs), the average by 1.5e-5 to 2.6e-4. At D = 5 the coarser h = 0.0141 fails this fit (an O(h²) floor
in H of about 1e-3 hides the growth), which is why the pair is (0.00707, 0.00354). This was found in development, before the runs of record.

**Mutations** (`checks.py` MUTATIONS, each run of record fires exactly its list): theta_coef (P1 P2 P3 C1 C2 G1 F1), scalar_coef and ghost_sign
(P1 P3 C1 C2 G1 F1), no_corrector (P2 P3 C1 C2), numpy_Spert (K), data_constraint_sign and ellis_charge (I1 P1 C1 C2 G1 F1), constraint_sign (C1),
flux_sign (C2), mass_exponent (P2 C2), finder_sign (P3 G1 F1 S), finder_swap (P2 P3), time_unit (G1), cut_extreme (S). Silent by design: the kernel
mutations do not fire K (they act on both kernels); no_corrector does not fire P1 (the predictor-only update is second order on the static solution).

**What was built differently from the design.** P2 is the exact Kruskal Tangherlini that contains the horizon (the ruling asks for the horizon at
(2M)^{1/(n−1)}), not the exterior-only static case of design §8. F1 uses a Gaussian in ∂φ, not in φ (a φ-Gaussian has almost no linear overlap with the
static flux: E₀ = +0.0096A − 4.66A² at h = 0.02). G1's D = 5 pair is finer than D = 4's scaled by 1/√2. Details: raw README, "Changes made during development".

## Cost

One thread, numba kernel, `nice -n 10`, dh610-fedora, load average 18 of 20, 30 runs per bisection to |A − A*ₕ| ≤ 1e-10 (every run ends at its event).

| case | rectangle | one bisection (30 runs) | per run | numpy kernel, same bisection |
|---|---|---|---|---|
| D = 4, h = 0.02 | X = 20, N = 1000 | 9.8 s | 0.33 s | 51 s |
| D = 4, h = 0.01 | N = 2000 | 39.8 s | 1.3 s | 187 s |
| D = 5, h = 0.00707 | X = 15, N = 2121 | 43.4 s | 1.4 s | 218 s |
| D = 5, h = 0.00354 | N = 4243 | 156 s | 5.2 s | 581 s |

This is 5 to 6 times the design's kernel-only estimate (0.22 s per run at D = 4, h = 0.01): the Python per-row diagnostics and the machine's load
are the likely reasons (not separated). It does not change the conclusion of the design: a bisection is minutes, at 4 threads all of #40 and #150 is under an hour.
The numpy kernel is 3.7 to 5.2 times slower and needs the whole grid in memory (2.8 GB at N = 4243).

## What #40, #150 and #213 can do now

- **#40, #150.** `ellis_rays` with a pulse family, `evolve(d, stop=True, store=False)` and the bisection of `checks.bisect_A` give p*ₕ and the
  lifetime series; `fit_growth` and `check_G1` are the template for a fit over a declared window. To keep in mind:
  - Lifetimes depend on the event definition (the tolerance 1e-2 of the trapped test, expansion at r ≥ 2); #40 declares its own.
  - A slope in x⁻ is 1/(√2W); the code reports t and τ.
  - p*ₕ ∝ h², so report it at three h and use a window above the O(h²) floor of H (D = 5 needs h ≤ 0.0071 here).
  - The G1 family is symmetric two-sided and passes through the static wormhole (p* ≈ 0). A one-sided pulse also excites the asymmetric zero mode, so
    the critical static member may not be the symmetric one (#42). A generic family needs a second parameter.
  - The threshold mass is not measured here: F1 gives M at the first trapped cell for pulses far from threshold (0.4975, 0.4920), not a settled mass, and
    the grid stretches after a horizon forms (e^{−f} grows; excision begins at 25).
- **#213.** The edge mass `mass_nodes(res)[-1, :]` with `edge_flux_residual`, T₊₊, T₋₋, T^θ_θ and ρ (`stress`), the masks and contours, and both throat definitions
  (`res.diag` on x⁺ = x⁻, and `throat_min`) are in. What `probe_expansion.py` shows (`probe_expansion_output.txt`; D = 4 and 5, two h each, values agree
  between the h to 3 to 4 digits): the expanding case runs to t = 42 with no excision and max e^{−f} = 1.00, so this gauge has not reached a coordinate limit there.
  On the diagonal r grows about linearly in t (D = 4: 5.0, 13.5, 25.0 at t = 12, 24, 40) while its proper time saturates (7.43, 8.06, 8.32), and H
  rises (0.75, 1.21, 2.07, 2.88 at t = 8, 12, 24, 40) instead of tending to a constant. That is observed here, not interpreted. The pulse is one-sided, and F1 shows that the diagonal and the minimum-r slice differ at the end of the
  expansion run (r = 8.52 and 8.18). #213 should choose its observer and its definition of w before reading H(t).

**Update from #40:** the bisection and the lifetime law are in [threshold.md](threshold.md) (`threshold.py`, `threshold_output.txt`), with the
coordinator's second code for family S in `threshold_check.py` (`threshold_check_output.txt`). For S and O the slope is 1/W to 3e-6, 3e-5 and 5e-4
at D = 4, 5, 6 over [1e-9, 1e-4], converging at second order in h; the generic one-sided family G departs from 1/W (cause open).

**Update from #213:** the symmetric expanding branch, followed along the throat, is in [expansion.md](expansion.md)
(`expansion.py`, `expansion_output.txt`, figure `expansion_H.svg`). The throat's rate does not settle: H = k/(τ_s − τ) with k = 2.0000
(D = 4) and 3.99 (D = 5), not a de Sitter rate; the coordinator's second code reproduces k, H/W and τ_s.

**Update from #226:** the null-constraint maximum that stops falling at the finest h at D = 5 and 6 is roundoff carried along the rays,
not truncation, and it does not bear on p*ₕ or the slope ([constraint_plateau.md](constraint_plateau.md)).

**Not built** (design items left out of the minimal build): symmetric-reduction by reflection, the asymmetric static family and Tangherlini in the
tortoise gauge as data, a D = 6 growth control (W²a² = 4.68662054311946479 is known), the Bondi-type extrapolation of the edge mass, a regular centre
(#45, #46), long-run output beyond the diagonal series (the full grid is kept when `store=True`).

## Files

`dn.py`, `checks.py`, `probe_expansion.py` (+ `probe_expansion_output.txt`), `derive.py` (+ `derive_output.txt`), `cost_model.py` (+ `cost_model_output.txt`),
`bench/` (the design's kernel stand-ins, not the evolution), `design.md`, `domain.svg`; `threshold.py` (+ `threshold_output.txt`,
`threshold.md`, #40) and `threshold_check.py` (+ `threshold_check_output.txt`, the coordinator's review of #40); `expansion.py`
(+ `expansion_output.txt`, `expansion.md`, `expansion_H.svg`, #213); `constraint_plateau.py` (+ `constraint_plateau_output.txt`,
`constraint_plateau.md`, #226).

## Coordinator's review

Under review: `wip/i39-throat-evolution` at 3a233b6, raw `I39-dh610` at ae8bff6. Review folder
`negative-mass-raw/I39-mac-review`, whose README holds the criteria, committed before the first run of `ind39.py`. The
changed files are those the brief names (`notes/wormhole_evolution/` and the raw folder).

- **Rerun on the Mac:** `derive.py` 25 of 25. `checks.py` with numba 11 of 11. With numpy and numba absent, 10 of 11
  pass and K is skipped. These are the verdicts of #39's runs of record.
- **ind39 K1, the equations.** The review's own brute-force Ricci tensor of ds² = −2e^{−f}du dv + r²dΩ_n² (n = 2, 3), with
  R_ab = Σ ε ∂_aφ ∂_bφ and □φ = 0, gives the three evolution equations of `dn.py`'s docstring with sympy residual 0.
  `dn.F_np` agrees with them to 3e-16 at 40 random states, with a ghost and a normal field. The constraint comes out as
  r_uu + f_u r_u = −(r/n)Σ ε φ_u². The mutation `theta-coef` fails it.
- **ind39 K2, a case #39 did not run: the tuned growth rate at D = 6.** With `checks.py`'s own G1 protocol, at
  h = 0.005/√3 and 0.0025/√3, W_fit = 2.1604 and 2.1625 against W = 2.1649 (−0.2% and −0.1%). The two agree to 1e-3,
  and the A*_h ratio is 4.005. *Added by the coordinator:* the first pair the review chose, h = 0.01/√3 and 0.005/√3,
  failed (K2, not met as worded), because at h = 0.0058 the fit is meaningless. So D = 6 needs h ≤ about 0.003, and #40
  was told. K2b (amendment A1, same bars) passes.
- **Development departures** (reported by the runner, all before the runs of record): P2 in Kruskal form, F1's pulse in
  dφ, the D = 5 pair at h = 0.00707 and 0.00354, and a cost 5 to 6 times the kernel-only estimate. A `check()` that
  failed to count numpy bools was fixed before any run of record.
- *Added by the coordinator:* the expansion probe (one-sided, to t = 42) shows no coordinate limit. The diagonal's proper
  time saturates near 8.3 while H rises from 0.75 to 2.9. Like the runner, I give no reading of this; #213 chooses its
  observer first.
- **Follow-ups.** #40 runs on this code. The runner's recommended tasks go to the issues they serve:
  - task 1, the asymmetric static family and the one-sided threshold, to #42;
  - task 2, the expansion observer with a symmetric long run, to #213;
  - task 3, the settled black-hole mass, to #41.

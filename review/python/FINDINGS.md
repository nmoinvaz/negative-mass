# Review: Shirokov GRTeclyn — initial data and analysis

Repo: https://github.com/Nikchik-coder/GRTeclyn (branch `develop`, commit 09fbcd6). Paper: arXiv:2608.24577.
Read-only review; nothing from the repo was executed. Independent check scripts are in this directory
(data read with `python3 -I`). The paper PDF text layer was not readable, so "paper numbers" are taken from the
pack READMEs and `grteclyn-wrapper/.../bondi_dipole/article_figures.py` (which the pack calls the arbiter).

Masses: true ADM (our solvers) M+ = 0.0146265, |M−| = 0.0145803, M̄_true = 0.0146034;
paper (enclosed at φ = 5% φ_c) M̄ = 0.0143223 = 0.9808 M̄_true.

## Findings (ranked)

### F1 — HIGH: the a·d² ≈ M̄ agreement rests on the coarse 128³ grid (verified at d=10; suspected elsewhere)
`article_figures.py:107-111` (SCAN) and `tools/fit_mass_law.py:163-168` use only N128 cells for the separation
scan and the mass ladder. At d=10, a(N256)/a(N128) = 1.0697.

| d | a·d²/M̄ paper mass (N128) | true mass (N128) |
|---|---|---|
| 8 | 1.0309 | 1.0110 |
| 10 | 1.0111 | 0.9916 |
| 12 | 1.0030 | 0.9837 |
| 16 | 1.0014 | 0.9821 |
| 20 | 1.0023 | 0.9830 |

Converged d=10: N256 1.0816 (paper mass; the paper's "8% above M/d²"), **1.0608 with true mass**; Richardson
(p=2, N192/N256) 1.075 / 1.054. Unless the resolution error is confined to d=10, converged a·d²/M̄ at d ≥ 12 is
~1.06–1.08, not 1.000; "within 0.1% from d=12" and "converged d=10 is 8% high" cannot both describe converged
physics. The force-law exponent survives only if the grid error is d-independent — untested (no N256 at d≠10).

### F2 — MED-HIGH: on the finer grids the acceleration is not constant (verified, re-checked by coordinator)
Core-tracker midpoint acceleration in four quarters of t ∈ [5, 200] (×10⁻⁴):

| grid | Q1 | Q2 | Q3 | Q4 | core sep start→end |
|---|---|---|---|---|---|
| N128 | 1.468 | 1.445 | 1.460 | 1.469 | 10.000→10.107 |
| N192 | 1.478 | 1.524 | 1.594 | 1.662 | 10.000→10.029 |
| N256 | 1.480 | 1.522 | 1.580 | 1.641 | 10.000→10.013 |

~11% rise on both fine grids at fixed separation, so not 1/d². "Steady to 2%" and the four-window consistency
come from N128 (including the t=400 run). The headline 1.549×10⁻⁴ averages a rising acceleration.
Coordinator note (untested hypothesis): the paper reports the stars' field content growing 12–21% by t=200
(decreasing with resolution); a growing |M| would raise a ∝ M̄. Worth checking against the pack's mass/field diagnostics.

### F3 — MEDIUM: downstream effect of the 5%-radius mass (verified)
Bias is frequency-dependent (−1.9% at ω=0.75, −3.0% at ω=0.85). Mass ladder (d=20, `mass_law.csv`): slope
0.966±0.061 → **0.982±0.062** with true masses; mean a·d²/M 0.9993 → **0.980**; per rung 0.994/0.978/0.974/0.973.
Separation ratios all drop 1.9% (F1). Pull-ratio predictions with true masses: w0804 0.8026→0.8059, w088
0.5997→0.6079, heavy 0.7471→0.7516. Negligible: matched-mass ratio (0.9961→0.9968), sign-flip point (~0.995→0.996).

### F4 — MEDIUM (interpretive): the signed-momentum test measures rigidity, not conservation (arithmetic verified)
`extraction/sector_dynamics.py:135-150, 288`: both sectors' momenta use the same canonical-sign integral and
`px_total = px_c − px_p` (col 13 = col 11 − col 12 exactly). With identical KG dynamics and nearly identical stars,
p+ ≈ p− whenever v+ ≈ v−, so any rigid co-moving motion passes, artefacts included. p−/p+ tracks v−/v+
(0.9940 vs 0.9885 at N256, t=190). |P_tot/P+| for t ≥ 20: N128 max 6.65% / median 1.08%; N192 1.71% / 0.57%;
N256 0.73% / 0.39%; grows toward the end on every grid. "<1%" holds only at N256. P ≈ 1.11 M_true v.

### F5 — MEDIUM (interpretive): the mirror cell is an exact parity image (verified from params)
Phantom at +5, canonical at −5 on the same centred grid; drift ratio −1.0000218, pointwise mismatch 2.2×10⁻⁵.
It proves only that the code has no left–right bias; it cannot separate physics from artefact (contrary to the
README). The mass-ratio sign reversal is the discriminating test.

### F6 — LOW-MED: tracker mixing in quoted numbers (verified)
Pull ratios quoted (0.809/0.597/0.739; controls 0.973/1.052/0.725) reproduce only with the barycentre tracker,
but `article_figures.py:1055` uses the core tracker (0.810/0.618/0.756; controls 1.011/1.028/1.001). On core,
w088 is off by +3.1% (paper mass) / +1.7% (true mass), not 0.05%. "Four disjoint fits 1.473/1.442/1.456/1.469"
are core-tracker; the headline a = 1.4481×10⁻⁴ is barycentre, whose windows spread ±12.5%. Separation-scan
barycentre window half-spread per cell 11/13/14/24/40% (core 1.5–2.6%); the quoted ±0.01 slope error ignores it.
Slope: barycentre −2.028 (all), −2.011 (no d=8), −2.002 (d≥12); core −2.040±0.010, −2.025, −2.017.

### F7 — LOW (reproducibility): most `analysis/*.py` and CSVs belong to a retired campaign (verified)
`separation_scaling.py:33-36`, `newtonian_reference.py:33-41` (ω=0.55 masses 0.06395/−0.07696, `convA_*`),
`momentum_balance.py:37-42` (`pair_pm_v2`), `convergence_check.py:30-41`, `make_tables.py:33-38`, and
`wave_check.csv`/`constraint_check.csv` reference cells not in the pack; the README (301-305) says that campaign
was contaminated by an initial-data offset "the same size as the signal", yet README:189-204 presents its tables
as current. Current numbers come from `article_figures.py`/`fit_mass_law.py`, the latter recomputing masses with
the author's solver — not re-derivable from the pack alone.

### F8 — LOW: the two READMEs normalise differently
Campaign README divides by M+ = 0.014350; pack README by M̄. Campaign README lists the phantom as M = −0.014350;
the table header says −0.0142947.

### F9 — LOW: seam at r_m (verified, negligible)
`boson_star_ode.py:424-437`: lapse jumps by −4.3×10⁻⁵ (canonical) / +5.1×10⁻⁵ (phantom) with a 5% slope kink,
because the tail uses the 5%-radius M (0.01435) while the core lapse implies ≈0.01457. φ′ continuous to 0.3–0.4%;
painted tail within 0.5% to r=15; integrated matter beyond r_m within 0.1%. Lapse enters only via Π = −(ω/α)φ at
≤5×10⁻⁵; the table has no ψ column (GRTresna re-solves ψ), so evolved stars should carry ≈ true ADM mass and the
5% convention only mislabels them. The C++ table reader was not seen (patch hunks show linear interpolation).

## Checked and correct
ODE system (Lichnerowicz, maximal-slicing lapse, KG; ρ and ρ+S; ε only in metric sources; origin factors);
ω/α_inf and ψ_inf² r rescalings; ω iteration to 1e-5; painted profile vs ours ≤4×10⁻⁵ inside r_m; per-lump inputs
to GRTresna; every data column read in `article_figures.py` and `fit_mass_law.py`; the weighted fit. Reproduced
exactly: ladder drifts/accelerations (2.8815/3.0139/3.0016; 1.4481/1.5568/1.5491×10⁻⁴), last-two-thirds fits,
end gaps, slope −2.0278±0.0109 (−2.041 four-point), `mass_law.csv` and slope, mirror ratio, N256 momentum fractions,
mean-shift drag 1.07×10⁻⁵ → 4.0×10⁻⁴. Doubled box: core drift moves 0.2%; the README's 4% is barycentre sensitivity.

## Not checked
The paper text itself; the C++ GRTresna reader and evolution code (separate review); ψ₄ "r^−4.8" falloff (quick
fit gives −3.2 to −4.2 depending on time, R16/R40 ratio 26–55 vs paper 31 — convention not established);
quadrupole/luminosity budget; constraint, energy-condition and curvature numbers; merger times. Whether the N128
grid error depends on d — the crux of F1 — needs N256 runs at d = 16 or 20.

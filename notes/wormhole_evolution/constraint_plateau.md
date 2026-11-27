# Why the constraint maximum stops falling at the finest h (#226)

**Answer.** The plateau of #40's constraint maximum is roundoff, carried along the rays. Over a fixed slab the maximum falls by 4 per halving until the roundoff term, which grows as 2⁻⁵²/h², passes the truncation term; in binary128 every halving of the declared h is 4.0 to 4.2 in every set at D = 5 and 6, and a closed-form model of the roundoff matches the extra constraint to 0.71 to 1.14 (*numerical evidence*, not a proof). The ray data's RK4 step, a single cell and the corner are ruled out.

Code `constraint_plateau.py`, output `constraint_plateau_output.txt`, raw `negative-mass-raw/I226-rpi5` (README: development runs, criteria written before the runs of record, mutations). `dn.py`, `checks.py` and `threshold.py` are unchanged.

## What was measured

The fixed region is a slab of the (x⁺, x⁻) rectangle, the nodes with t = (x⁺ + x⁻)/√2 ≤ T_c, the same at every h. The streamed constraints are `dn.constraints`' formulas (control K0, difference 0.0), the evolution is `dn.row_step`, and the grid is never stored.

| set | what runs | T_c at D = 5, 6 | h at D = 5 | h at D = 6 |
|---|---|---|---|---|
| static | the Ellis wormhole, no pulse | 4, 4 | 0.00707, 0.00354, 0.00177 (+ 0.000884) | 0.005, 0.0025, 0.00125 (+ 0.000625) |
| S ∓ | family S at p*ₕ − 1e-8 and + 1e-8 | 10, 8 | same | same |
| bracket | #40's two amplitudes lo, hi (1e-14 apart) | 16, 12 | three | three |

The extension h is below #40's range and has its own p*ₕ, bisected here (D = 5: 4.50511e-6, D = 6: 1.734492e-5, brackets 1.8e-10 and 1.7e-10). T_c is below the first departure of the coarsest static run (the static
wormhole's unstable mode, seeded by truncation error, departs at t = 6.7, 7.5, 8.3 at D = 5 and 5.1, 5.7, 6.4 at D = 6) and below the first event of every run of the other sets.
A halving is second order if its ratio is in [3.364, 4.757], the band of #39's order test.

```
 x⁻
  ▲          throat x⁺ = x⁻ (r = 1)
  │        ╱
  │      ╱ ◉  max |C| of the slab: late edge, on the throat
  │    ╱   ▲
  │  ╱     │  C·rⁿ⁻¹ is constant up the ray; r falls from r₀ ≈ x⁺/√2 to 1
  │╱       │
  ●────────┴──────────────▶ x⁺      roundoff 2⁻⁵² r is made in every node of every row
```

## Results

Maximum of |C₊₊|, |C₋₋| over the slab (float64 unless said), ratios of the declared halvings in the last row.

| D = 5 | static | S −, float64 | S −, binary128 | S +, float64 | bracket lo, float64 | bracket lo, binary128 |
|---|---|---|---|---|---|---|
| h = 0.00707 | 4.578e-05 | 4.618e-05 | 4.618e-05 | 4.606e-05 | 4.606e-05 | 4.605e-05 |
| h = 0.00354 | 1.140e-05 | 1.165e-05 | 1.153e-05 | 1.159e-05 | 1.183e-05 | 1.146e-05 |
| h = 0.00177 | 2.896e-06 | 3.464e-06 | 2.882e-06 | 3.498e-06 | 6.251e-06 | 2.865e-06 |
| extension h = 0.000884 | 9.342e-07 | 5.883e-06 | 7.205e-07 | 4.258e-06 | not run | not run |
| ratios, declared | 4.015, 3.937 | 3.962, 3.364 | 4.004, 4.001 | 3.973, 3.314 | 3.893, 1.892 | 4.018, 4.001 |
| ratio, extension | 3.100 | 0.589 | 4.000 | 0.822 | | |

| D = 6 | static | S −, float64 | S −, binary128 | S +, float64 | bracket lo, float64 | bracket lo, binary128 |
|---|---|---|---|---|---|---|
| h = 0.005 | 7.841e-05 | 7.776e-05 | 7.775e-05 | 7.776e-05 | 7.776e-05 | 7.775e-05 |
| h = 0.0025 | 1.865e-05 | 1.909e-05 | 1.860e-05 | 1.905e-05 | 2.192e-05 | 1.860e-05 |
| h = 0.00125 | 4.963e-06 | 9.883e-06 | 4.610e-06 | 8.392e-06 | 2.055e-05 | 4.599e-06 |
| extension h = 0.000625 | 3.713e-06 | 2.561e-05 | 1.152e-06 | 2.662e-05 | not run | not run |
| ratios, declared | 4.204, 3.758 | 4.073, 1.932 | 4.181, 4.034 | 4.081, 2.270 | 3.547, 1.066 | 4.181, 4.044 |
| ratio, extension | 1.337 | 0.386 | 4.003 | 0.315 | | |

Bracket hi is the same as lo (D = 5: 3.900, 1.829; D = 6: 3.701, 1.017). The binary128 columns are the numbers in `constraint_plateau_output.txt`, C1, which has the other binary128 sets.

- **Static (R1, met).** The ratio is second order at the declared h in the departure-free slab (D = 5: 4.015, 3.937; D = 6: 4.204, 3.758). The plateau shows on the static case only at the extension h (3.100, 1.337) or with added noise.
  The maximum is on the diagonal x⁺ = x⁻ at the slab's late edge (declared h: t = 3.85 to 4.00, x = 2.72 to 2.84), where the truncation term has grown slowly.
- **S ∓ 1e-8 and the bracket (R2, R3, measurements).** At D = 6 the finest declared halving is 1.93 (S −), 2.27 (S +), 1.07 and 1.02 (bracket lo, hi); at D = 5 it is 3.364 (S −, on the band's edge), 3.314 (S +), 1.89 and 1.83. The ratio at the extension halving is 0.3 to 0.8: the maximum rises as h falls.
  The maximum is on the diagonal near the slab's late edge (D = 6, S −: x = 5.46, t = 7.72; bracket: t = 12.0; D = 5, S −: x = 6.85, t = 9.68; bracket: t = 15.8 to 16.0). #40's bracket events were at t of 13.7 or later (D = 6) and 17.3 or later (D = 5); the slabs end before them.
- **Binary128 (C1, met).** The same slabs in 113-bit precision, data included, fall by 4.0 to 4.2 at every declared halving of every set (static, S −, S +, bracket lo) and 4.000, 4.003 at the extension halving. Where the float64 ratio fails, the float64 to binary128 ratio of the maximum is
  1.2 (D = 5, h = 0.00177, S −) and 2.1 (D = 6, h = 0.00125), 8.2 and 22 at the extension h.
- **Time profile.** In binary128 the maximum per time bin is flat in t (D = 6, h = 0.00125: 2.94e-06 at t = 0.25, 4.6e-06 from t = 3.25 on). In float64 it grows with t: the float64 to binary128 ratio is 1.003, 1.06, 1.19, 1.29, 1.96 at t = 1.25, 3.25, 5.25, 6.25, 7.25, so the excess over binary128 (in quadrature)
  is 1.6e-06, 3.0e-06, 3.7e-06, 7.7e-06 at t = 3.25, 5.25, 6.25, 7.25. It grows as the amplification of a violation from the ray, (t/2)ⁿ⁻¹ (output, section P), but more slowly than that factor.
- **Added noise (C2, met).** Uniform relative noise a·2⁻⁵² added to every node as it is computed moves the finest declared h in proportion to a: D = 6, S −, a = 0, 1, 16, 256: 9.88e-06, 1.17e-05 (×1.19), 6.68e-05 (×6.8), 1.11e-03 (×113); D = 5: 3.46e-06, 3.91e-06 (×1.13), 1.50e-05 (×4.3), 1.97e-04 (×57).
  The coarsest h does not respond at a = 16 (×1.001, ×1.005) and does at a = 256 (×1.33, ×1.13). The issue's mutation, a = 1 (doubling the roundoff), is inside the 13 % scatter of a single maximum (D1), so a ≥ 16 is the resolved test.
- **Ray data (D1, met).** The data's RK4 constraint solve at h/4 and h/16 (its own step against the grid's) changes the finest-h maximum by ×1.08, ×1.03 (D = 5) and ×0.87, ×0.98 (D = 6); a planted bump of 1e-7 in r on the rays gives ×15.8 and ×27.5.
- **A single cell (L1, met).** On the slab's last diagonal at least 36 nodes lie within a factor 2 of the maximum in every run (S −: 179 and 120; bracket lo: 59 and 48 for C₊₊); a planted spike of 1e-10 in one node leaves 2.
- **The corner (L2, met).** The maximum is at x⁺, x⁻ ≥ 2.7 in every declared-h run; a spike at node (2, 2) puts it at x⁺ = 0.001.

## The cause: roundoff carried along the rays (numerical evidence)

1. Roundoff d (relative 2⁻⁵², about a·2⁻⁵² r with a ≈ 0.9) in a node of r is a step in the quadrant above it, because the update is N = E + W − S + h²F. C₊₊ is a second difference in x⁺, so the step becomes the pair +d/h², −d/h² in C₊₊ at the two columns beside the node on every later row.
2. On the evolution equations C₊₊ rⁿ⁻¹ is constant along x⁻ (∂₋C₊₊ = −(n − 1)(r₋/r)C₊₊, the coordinator's sympy check in #40). A pair made where r is large is amplified by (r₀/r)ⁿ⁻¹ when the ray falls to the throat: about 50 (D = 5) and 360 (D = 6) for r₀ ≈ 7, x⁺ = 10.
3. Independent noise in each row adds in quadrature:
   `rms C₊₊(i,j) = √(2/3) · a · 2⁻⁵² / h² · r(i,j)^−(n−1) · (Σ_{j'<j} r(i,j')^{2n})^{1/2}`, with the static r(l), l = (x⁺ − x⁻)/√2. This shares none of the evolution (M1).
4. Observed over the model, by |l| band at t = 6 (D = 5, 6; h₂, h₃; C₊₊ and C₋₋): natural roundoff (float64 minus binary128) 0.71 to 1.14, injected a = 16: 0.93 to 1.14, a = 256: 0.91 to 1.14 (48 values each). With exponent n for n − 1 the ratios are 2 to 4.5 for |l| ≥ 1 (mutation `model_exp`).
5. Truncation falls as h², roundoff rises as h^−5/2 at a node (ε/h² times √rows), so the maximum falls, crosses and rises: ×1.2 and ×2.1 above the truncation floor at the finest declared h, ×8 and ×22 at the extension h.

The effect is in the diagnostic, not the solution. Exploratory (not declared; `supplement/r_roundoff.py`): on the static diagonal t = 3.5 the roundoff error of r is 1.4e-11 (D = 6, h = 0.00125) and 4.2e-12 (D = 5, h = 0.00177) against a truncation-seeded difference
from the exact static r(l) of 4.1e-4 and 6.4e-5. That is why p*ₕ and the slope keep their second-order ratios while the constraint maximum does not. At t = 7 in family S the roundoff error of r has grown to 2.5e-8 (D = 6), a factor 1700 for 3.5 time units against e^{W·3.5} = 1950 (a rough comparison: the static case at t = 3.5 and family S at t = 7).
That growth is the roundoff seed of the unstable mode; it bears on #227's floor.

## What the notes should now say (proposals for the coordinator)

`threshold.md`, "Constraint diagnostic": replace "its plateau is not explained here. Candidates not tested: roundoff in second differences of a field that has grown from the 1e-13 seed, and a local maximum at one cell" with:

> At the finest h the maximum is roundoff. The second difference of r carries roundoff 2⁻⁵² r per node as ε/h², and C rⁿ⁻¹ is constant along x⁻, so a violation made at large r is amplified by (r₀/r)ⁿ⁻¹ at the throat (#226).
> In binary128 the same slabs fall by 4.0 to 4.2 per halving at D = 5 and 6, and noise of a·2⁻⁵² per node raises the finest-h maximum in proportion to a (×1.2 to ×2.1 above the truncation floor at the finest declared h, ×8 to ×22 at one more halving). The maximum is not a single cell and not the corner, and the ray data's RK4 step does not change it.
> The roundoff in r is 1e-11 against 4e-4, so the plateau does not bear on p*ₕ or the slope.

Add under the table: the finest-h entries at D = 5 and 6 are roundoff, not truncation. Recommended task 2 ("The constraint plateau at the finest h") is done by #226. Paper III's limits are unchanged.

## Disagreements and limits

- #40's bracket amplitudes do not reproduce #40's events on this machine (D = 5, h = 0.00707, lo: expansion at cell 2076 against #40's collapse at 1822), as expected from a floor of p*ₕ of 1e-12; the max |C| at the coarse h agrees (D = 6, h = 0.005: 7.7755e-05 against 7.7753e-05).
  The runs at p*ₕ ∓ 1e-8 do not depend on this. #40's event-past maxima at the finest h (D = 6: 5.28e-05) are above the bracket slabs here (2.06e-05, t ≤ 12) because the event's past reaches t = 13.7 and later, where the amplification is larger.
- R2 and R3 are not met at the finest declared halving, except D = 5, S − (3.364, on the edge). The declared prediction named D = 6 and the edge at D = 5; D = 5's bracket (1.89, 1.83) was not predicted and is not met either.
- The static region cannot be run past t ≈ 5 to 6 at D = 5 and 6, so the plateau does not appear on the static case at the declared h; it appears at the extension h and with a ≥ 16.
- The model is compared with the runs (within 0.71 to 1.14 here), not proved. Binary128 shares the scheme and the constraint formulas with the float64 runs and changes only the precision. A single maximum has a scatter of about 13 % (D1).
- Not done: D = 4, family G, the first-difference (flux) form of the constraint, other machines (the natural roundoff a ≈ 0.9 is aarch64 numba and numpy).

## Recommended tasks

1. *A constraint diagnostic without the 1/h² roundoff.* What: the mass-flux form of the constraint (`edge_flux_residual`, first differences only) on #40's S ∓ 1e-8 and bracket slabs at D = 5, 6. Why: the second-difference C is a poor diagnostic once ε/h² passes h² (about h = 0.0015).
   Acceptance: the maximum falls by 4 ± the order band at the three declared h and the extension h; noise a = 256 raises it less than it raises C (×113, ×57 here); a planted constraint violation of 1e-7 is seen. Depends on: this note.
2. *Roundoff seed of the unstable mode and the p*ₕ floor (#227).* What: measure |r_float64 − r_binary128| on the stall against t at D = 4, 5, 6 and compare its growth rate with W and its extrapolated size at t = 0 with the 1e-11 floor of #40. Why: 2.5e-8 at t = 7 is 1700 times the 1.4e-11 at t = 3.5, e^{W·3.5} = 1950.
   Acceptance: the growth rate agrees with W within 2 %; the floor in p* predicted from the seed agrees with #227's within a factor of 3. Depends on: #227.

## Files

`constraint_plateau.py` (controls, runs, analysis; `python3 constraint_plateau.py analyze --raw=<raw>/I226-rpi5/runs`, `controls`, `--mutate=NAME`), `constraint_plateau_output.txt` (from raw `I226-rpi5/make_output.sh`, run from this directory), this note.
Raw: `I226-rpi5/` (`README.md`, `jobs_wave1.txt`, `jobs_wave2.txt`, `launch_I226.sh`, `pstar_table.json`, `pstar_ext/`, `runs/`, `logs/`, `dev/`, `supplement/`).

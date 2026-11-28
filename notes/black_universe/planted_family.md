# The exact family of #269 as planted truth for #217's nonlinear double-null code (#277)

Code [planted_family.py](planted_family.py), output of record [planted_family_output.txt](planted_family_output.txt), assembled from negative-mass-raw `I277-mac/runs`; the criteria, written
before the runs of record, and the development logs are in `I277-mac/README.md` and `I277-mac/dev`. Builds on #269 ([exact_mode.md](exact_mode.md)) and #217 ([evolution.md](evolution.md),
[evolution_design.md](evolution_design.md)). Status words as in AGENTS.md: every error, ratio and v\* below is *numerical evidence*; D1b, D2 and D3 are identities checked at the stated precision.
No sentence here says a result is new (this note was not searched against the literature).

## Answer

**#217's kernel reproduces the exact nonlinear, time-dependent family of #269 from characteristic data on two rays.** At c = −1.5 and −1, β = ±0.05, the sup error of the evolved r against
λ(v)R(ρ) over the whole rectangle falls by 3.977 to 4.003 per halving of h = 0.02 … 0.00125 (16 of 16 ratios inside 4.0 ± 0.3), to 7.4e-8, 1.7e-6, 2.1e-6, 6.9e-7 at h = 0.00125.
For β < 0, v\* = ln(−L/β)/K is found from the evolved r² cos²(φ/√2) to 2.0e-4 at the coarsest h and 7.9e-7 at the finest at the worse member (bar 1e-3). The β sign mutation fails both. 12 of 12 mutations fire exactly
their declared checks; the family's residuals in #217's equations are 2.0e-14 at worst.

## 1. The chart map

The family is in the ingoing null chart (v, ρ); the code is double null with data r, f, φ on x⁻ = 0 and on x⁺ = 0. The map, from the family's F:

- **x⁺ = (v − v₀)/√2** everywhere (x⁺ and v are both constant on ingoing rays); the ray x⁻ = 0 is the outgoing ray through ρ_c = 1, parametrized by it.
- **x⁻ labels the outgoing rays**, dρ/dv = F(v, ρ)/2. On the ingoing ray v = v₀ its parameter is the tortoise-like one, dρ/dx⁻ = −F(v₀, ρ)/√2, that is e^{−f} = F(v₀, ρ) there
  (at β = 0 this is #217's e^{−f} = A).
- **The gauge on the other ray follows.** The metric is 2ρ_{x⁻} dv dx⁻, so e^{−f} = −√2 w with w = ∂ρ/∂x⁻ at fixed v, and w_v = F′w/2. Along an outgoing ray
  dρ/dv = F/2 and d ln e^{−f}/dv = F′/2, with F′ = ∂_ρF = A′ − 2λ̇/λ.
- **Read-out without the map.** Λ² = r² cos²(φ/√2) = r²/(1 + ρ²) = λ(v)² at every x⁻, from the evolved r and φ alone.

```
 x⁻ (label of the outgoing rays)
 16 ┬ ρ₀ → ρ_F(v₀), F(v₀, ρ_F) = 0:  e^{−f} → 0
    │   left ray v = v₀:  dρ/dx⁻ = −F(v₀, ρ)/√2,  e^{−f} = F(v₀, ρ),  r = λ(v₀)R,  φ = √2 arctan ρ
    │  ┌────────────────────────────────────────────────────────┐
    │  │ evolve.row_step, unchanged                             │  exact: ρ(v; ρ₀) from dρ/dv = F/2,
    │  │ v = √2 x⁺,  r = λ(v)R(ρ),  φ = √2 arctan ρ             │         d ln e^{−f}/dv = F′/2
  0 ┼──┴─ corner (v₀ = 0, ρ_c = 1) ─── right ray x⁻ = 0:  ρ(v), e^{−f} by the same flow ───┬─→ x⁺
    0                                                                                     2.16 (c = −1.5)
                                                                                          3.20 (c = −1)
```

The rectangles end where λ² = 0.076 (c = −1.5, v ≤ 3.0547) and 0.108 (c = −1, v ≤ 4.5255) at β = −0.05; v\* = 3.137124 and 4.705685. Both are chosen from development runs (README, dev_05).

## 2. The data path

| id | what | result (bar) |
|---|---|---|
| D1 | closed forms r, F, and V, V′ of `evolve_potential.pot` against exact_mode's sympy forms, 50 digits | 2.2e-15 (1e-13) |
| D1b | hand-written jets (13 of them, mixed and second derivatives) against the chain-rule jets of exact_mode's forms, 50 digits | 6.5e-15 (1e-12) |
| D2 | the family's residuals in #217's evolution equations: `F_np` and the numba `_F` on the jets, r₊₋, f₊₋, φ₊₋, 400 points per case | 2.0e-14 (1e-12); per case below |
| D3 | null constraints C₊₊, C₋₋ on the jets, 400 points per case | 8.2e-14 (1e-12) |
| D4 | C₊₊ on the stored right-ray data, central differences, relative to its largest term, at h = 0.00125 | 3.5e-6, 6.1e-5, 5.3e-7, 1.4e-5 (1e-3); falls by at least 3.1 per halving (2.5) |
| D5 | the exact map with tolerances tightened to 1e-13 | 5.1e-12 (1e-10) |
| D6 | the module's row sweep against `evolve.run_grid`, four cases | numba bit-identical; numpy 0 (1e-12) |
| PL1 | planted errors C hᵖ through the ratio band | p = 2 accepted (4.000); p = 1.8, 2.2, 1 rejected (3.48, 4.60, 2.00) |
| PL2 | planted Λ² = 1 + βe^{Kv} through the v\* finder | recovered to 2e-15; noise 1e-5, 40 seeds: 3.7e-6, 6.3e-6; the finder with the rate × 1.01 misses by 3.0e-3, 4.8e-3 |

D2 per case (max relative residual over 400 points): F_np r₊₋ / f₊₋ / φ₊₋, numba `_F`.

| c | β | F_np r₊₋ | f₊₋ | φ₊₋ | `_F` |
|---|---|---|---|---|---|
| −1.5 | +0.05 | 2.0e-14 | 1.3e-15 | 1.1e-15 | 2.0e-14 |
| −1.5 | −0.05 | 1.2e-14 | 2.5e-15 | 1.6e-15 | 1.2e-14 |
| −1 | +0.05 | 6.0e-15 | 8.1e-16 | 6.8e-16 | 6.0e-15 |
| −1 | −0.05 | 6.2e-15 | 1.7e-15 | 7.3e-16 | 6.2e-15 |

**What is shared.** The runs use one `Family` for the data and for the exact fields, so an error in F would not show in O. It is caught by D1 and D1b (against exact_mode's sympy forms, which share nothing
with this module), by D2 (F's jets must satisfy the kernel's equations) and by D3 and D4 (the constraints). The kernel and `F_np` are #217's; the coordinator's `ind217.py` ties `F_np` to a brute-force Ricci
tensor (1.1e-15) and `ind269.py` ties the family to one (9e-40). The potential is shared between the code and exact_mode only through its formula, written twice (D1).

## 3. The evolution: error against λ(v)R(ρ)

Sup |r_code − λR| over every node of the rectangle (the whole of x⁺ ∈ [0, 2.16] or [0, 3.20], x⁻ ∈ [0, 16]; f reaches 11 to 14), and its ratio per halving of h.

| c | β | h = 0.02 | 0.01 | 0.005 | 0.0025 | 0.00125 | ratios |
|---|---|---|---|---|---|---|---|
| −1.5 | +0.05 | 1.8986e-05 | 4.7454e-06 | 1.1862e-06 | 2.9653e-07 | 7.4132e-08 | 4.001, 4.000, 4.000, 4.000 |
| −1.5 | −0.05 | 4.1930e-04 | 1.0542e-04 | 2.6397e-05 | 6.6021e-06 | 1.6508e-06 | 3.977, 3.994, 3.998, 3.999 |
| −1 | +0.05 | 5.3727e-04 | 1.3438e-04 | 3.3599e-05 | 8.3999e-06 | 2.1000e-06 | 3.998, 4.000, 4.000, 4.000 |
| −1 | −0.05 | 1.7700e-04 | 4.4213e-05 | 1.1048e-05 | 2.7615e-06 | 6.9029e-07 | 4.003, 4.002, 4.001, 4.000 |

**O: met** (16 of 16 ratios in [3.7, 4.3]). The ratios of f and φ, reported and not criteria, are 3.977 to 4.006 as well. At h = 0.00125 the sup sits at the last row (the latest v) at three cases and at v = 4.04
(row 2286 of 2560) at c = −1, β = −0.05; it sits at x⁻ = 16 except at c = −1.5, β = +0.05 (x⁻ = 5.9). At β < 0 the error is concentrated near v\*: c = −1.5, h = 0.00125, the running sup over the first three
quarters of v is 1.45e-8, 3.3e-8, 4.2e-8 and over all of v 1.65e-6. The sup error relative to r is 5.5e-8 to 1.4e-6 at h = 0.00125. The error is the discretization
error of the kernel against an independently integrated map (D5: 5e-12); it is not a fit error.

## 4. v\* for β < 0

v\* is the zero of a least-squares fit a + b e^{Kv} to the evolved Λ² averaged over the interior (rows i ≥ 1, columns j ≥ 1), v\* = ln(−a/b)/K. It uses no data row and not β.

| c | v\* | h = 0.02 | 0.01 | 0.005 | 0.0025 | 0.00125 |
|---|---|---|---|---|---|---|
| −1.5 | 3.137124 | −3.58e-05 | −8.85e-06 | −2.20e-06 | −5.49e-07 | −1.37e-07 |
| −1 | 4.705685 | +2.00e-04 | +5.03e-05 | +1.26e-05 | +3.16e-06 | +7.91e-07 |

**V: met** (bar 1e-3 at every h). The errors fall by 4.0 per halving (4.04, 4.02, 4.01, 4.01 and 3.97, 3.98, 3.99, 4.00), and the fit's (a, b) at the finest h is (1.000000, −0.050000) at both members. The direct
estimate from the last row, v − ln(1 − Λ²)/K (needs L = 1), gives −5.5e-5 … −2.2e-7 and +6.8e-5 … +2.6e-7.

## 5. Mutations (12 of 12 fire exactly the declared sets)

| mutation | what | fires | effect at h = 0.00125 |
|---|---|---|---|
| beta_sign | data from −β, compared with +β | O, V | sup error of r 0.44 (c = −1.5, β = −0.05), ratios 1.000; the fit has b = +0.05, no root |
| kh | coefficient of 1/λ² − 1: 1 → 1.01 | D1, D1b, D2, D3, D4, O, V | 1.2e-2, ratios 1.00; v\* off by 6.5e-3 |
| kf | coefficient of 2ρλ′/λ: 2 → 2.02 | D1, D1b, D2, D3, D4, O, V | 4.7e-2, ratios 1.00; v\* off by 1.1e-2 |
| rate | K → 1.001K inside λ² only | D1, D1b, D3, D4 | O, D2, V pass: ratios 4.00; v\* off by 3.0e-4, 4.9e-4 from the mutated family's own v\*, inside the bar |
| gauge_frozen | right ray's e^{−f} = F(v₀, ρ(v)) | D4, O, V | 0.18, ratios 1.00; v\* off by 3.5e-2 |
| th_coef, vterm_sign, vterm_coef, Vp_sign, V_sign_f | #217's own kernel mutations | D2, O, V | 37, 4.6, 0.35, 1.2, 0.26 |
| eps_sign | the kernel run with ε = +1 | D2, O, V | 7.8 |
| vfit_K | the v\* finder's rate × 1.01 | PL2, V | v\* off by −3.0e-3, −4.9e-3 |

**S: met.** The β sign mutation exits 1 and fires O and V in the run of record.

The `rate` row is a result. With the rate scaled by 1.001 the family still satisfies the three evolution equations at every point (D2 passes) and fails the constraint C₊₊ (D3, D4); the code reproduces it at
second order (ratios 4.00). That matches exact_mode.md Sec. 2, where the rate comes from the (vv) equation and the other four hold for every λ(v). O and V cannot see a wrong rate; D1, D3 and D4 can.

## 6. Where the exact map ends (`planted_family.py rays`, reported, added after the runs of record)

Outgoing rays labelled by ρ₀ at v₀ = 0, integrated with the family's flow.

- **β < 0.** Every ray runs to ρ = +∞ as v → v\*, with (1 + Kρ)√(K(v\* − v)) tending to a ray-dependent constant (2.66274, 1.37948, 1.37587 at c = −1.5 for x⁻ = 0, 8, 16;
  2.54692, 1.18706, 1.18228 at c = −1). Along each ray r = λR tends to a positive constant (2.788, 1.444, 1.441 at c = −1.5; 4.001, 1.865, 1.857 at c = −1). So r → 0 at v\* holds at fixed ρ, as
  exact_mode.md says, and does not hold along the outgoing rays, which the rectangle consists of. This is why v\* is read from Λ² and not seen as an event of the grid. By hand, from F ~ (1 + Kρ)/λ²
  and φ′² = 2/R⁴, the Ricci scalar along such a ray tends to 0; not computed here.
- **β > 0.** The rays run to ρ = −1000 at v = 6.92, 4.47, 4.43 (c = −1.5; x⁻ = 0, 8, 16) and 11.93, 5.35, 5.28 (c = −1). The rectangles end at v = 3.05 and 4.53. At c = −1.5 with x⁺ up to 4 (v = 5.66) the
  flow stops (dev_02).

Both ends are limits of the chart (v, ρ) of the family, not of the code.

## 7. Proposed text (the files are not touched here)

`evolution.md`, in the checks table after E7, and in "What this does not show" (first bullet) after "(a family not through it would, as #40's G)": *Update from #277
([planted_family.md](planted_family.md)): the kernel reproduces the exact time-dependent family of #269, r = λ(v)R, F = A − 2ρλ′/λ + 1/λ² − 1, from data on two rays at c = −1.5 and −1, β = ±0.05: the
sup error of r over the rectangle falls by 3.977 to 4.003 per halving of h = 0.02 … 0.00125, and v\* of the β < 0 members is found from the evolved r² cos²(φ/√2) to 2.0e-4 (h = 0.02) and 7.9e-7
(h = 0.00125) (numerical evidence). The rectangles reach v = 3.05 and 4.53, against lifetimes of T ≈ 15 to 20 here; nothing is said about the growth of the code's own unstable mode on this background.*

Paper IV (`paper/4-black-universes/main.tex`), Sec. "Nonlinear evolutions", after "the static members stay static to second order in the step $h$":

```latex
The code also reproduces the exact time-dependent solution of Sec.~\ref{sec:cauchy} from characteristic data on two rays. With $r=\lambda(v)R(\rho)$,
$\lambda^2=1+\beta e^{Kv}$ and $\beta=\pm0.05$ at $c=-1.5$ and $-1$, the supremum of $|r-\lambda R|$ over the rectangle falls by $3.98$ to $4.00$ per halving of $h$
(five steps from $0.02$ to $0.00125$), to between $7\times10^{-8}$ and $2\times10^{-6}$, and for $\beta<0$ the zero $v_*$ of the evolved $r^2\cos^2(\phi/\sqrt2)$ is found to
$2\times10^{-4}$ at the coarsest step and $8\times10^{-7}$ at the finest. The rectangles end where $\lambda^2\approx0.1$ and reach $v\le3.05$ and $4.53$, far short of the lifetimes below.
```

and in Sec. "Computations", the "Evolutions" item, append: *The exact family of Sec.~\ref{sec:cauchy} as planted truth: the code's error against it falls by $4.0\pm0.3$ per halving in four cases, the sign
mutation of $\beta$ fires, and twelve mutations each fire their declared checks; a wrong rate $K$ in $\lambda(v)$ is seen by the constraints and not by the evolution.*

Paper IV's sentence "For $\beta<0$ the areal radius vanishes at the finite $v_*$" holds at fixed $\rho$; Sec. 6 above suggests adding "at fixed $\rho$": *along the outgoing null rays $\rho\to\infty$ as
$(v_*-v)^{-1/2}$ and $r$ tends to a positive constant (numerical evidence)*. The label `sec:cauchy` is that of the subsection holding the family's paragraph; it was not compiled here.

## 8. What this does not show

- **Long times.** The rectangles reach v = 3.05 and 4.53. The family leaves the chart at v\* (β < 0) or at v ≈ 4.4 to 5.3 on the boundary rays (β > 0). The lifetimes of #217 (T ≈ 15 to 20) and the growth of the code's own
  unstable mode (rate W) on a time-dependent background are not tested.
- **Through v\*.** The data do not exist beyond v\*; v\* is found by extrapolating Λ² to zero from rows at λ² ≥ 0.076. In development, with the rectangle ending at λ² = 0.007 at c = −1, the error ratios per halving from
  h = 0.02 to 0.0025 are 3.49, 3.77, 3.92 (dev_05; not a criterion).
- **One corner, one L, two β.** v₀ = 0, ρ_c = 1, L = 1, β = ±0.05, c = −1.5 and −1, D = 4, one potential. The left-ray gauge is the tortoise-like one; other gauges of the data are not tried.
- **The representation of the data.** The data are exact at the nodes; interpolation of data onto a grid is not tested.
- **Asymptotically flat data.** The family is not asymptotically flat for β ≠ 0 and linearizes to Y_ref (exact_mode.md Sec. 4). Nothing here bears on #250.
- **The wrong-rate blind spot is closed only by D1, D3, D4**, which test the family's formula, not the code.

## Recommended tasks

1. **The family at small |β|, at #217's lifetimes.** *What:* run `check`'s runs at β = ±1e-3 and ±1e-4 (v\* = 7.2 and 9.6 at c = −1.5; 10.9 and 14.5 at c = −1) on rectangles to v ≈ 0.9 v\* (β < 0) and to the
   escape of the boundary ray (β > 0), h = 0.01 … 0.00125. *Why:* the exact background then stays near the static member until v ≈ 8 to 14, where #217's truncation-seeded unstable mode (rate W = 0.494 at c = −1.5, κ = K at c = −1) grows; this
   tests its amplification on an exact nonlinear solution. *Acceptance:* the error ratio is 4.0 ± 0.3 over the rectangle, or the first v at which the sup error leaves the h² law is reported and compared with #217's
   floor (T ≈ 16 to 20); the β sign mutation fires. *Depends on:* #277.
2. **How close to v\* the code stays second order.** *What:* β = −0.05 at both members with λ²_end = 0.04, 0.02, 0.01, 0.005 and h down to 0.000625. *Why:* dev_05 shows the first ratio at 3.49 when the rectangle ends at
   λ² = 0.007 (c = −1), and the v\* finder improves as the rectangle approaches v\*. *Acceptance:* the ratio at the three finest h is 4.0 ± 0.3 for λ²_end ≥ 0.01 at both members, or the λ²_end
   where it leaves the band is reported; the finder's error keeps ≤ 1e-3. *Depends on:* #277.
3. **The linear companion on the family's first-order data.** *What:* the β → 0 linearization (δr = (β/2)R e^{Kv}, δφ = 0, δF = −β(1 + Kρ)e^{Kv}) as data for #232's linear double-null solver, compared with
   this code at β = ±1e-4 divided by β. *Why:* exact_mode.md Sec. 8 lists "the double-null solver of #232 reproduces the first-order mode" as not run. *Acceptance:* agreement to 1e-6 relative at two h with second-order
   convergence; a mutation of K by 1.001 fails it. *Depends on:* #232, #269, #277.
4. **Planted family in `evolve.py check`.** *What:* add D2, D3 and a run at β = −0.05, c = −1.5, h = 0.01 and 0.005 to #217's E-series as E12. *Why:* the `rate` mutation shows that a time-dependent solution can pass
   O and V with a wrong rate; D3 keeps the constraint in the regression. *Acceptance:* E12 passes in under 10 s; the `kf`, `kh` and `rate` mutations fire it. *Depends on:* #277.

## Reproduce

From this directory, `~/Source/negative-mass/.venv/bin/python`, `NUMBA_NUM_THREADS=1`, `nice -n 10`, one process:

```
planted_family.py check [--mutate=beta_sign|kh|kf|rate|gauge_frozen|th_coef|vterm_sign|vterm_coef|Vp_sign|V_sign_f|eps_sign|vfit_K] [--out=DIR]    # 47 s; exit 1 when a check fails
planted_family.py mutations DIR check_       # fired sets of DIR/check_<name>.log against the declared ones
planted_family.py rays                       # Sec. 6
planted_family.py case --c=-1.5 --beta=-0.05 --h=0.005   # one run, JSON
```

`planted_family_output.txt` is `check`, `mutations`, the mutations' effects and `rays` of the runs of record (`I277-mac/runs`, assembled by `assemble.sh`), in that order. The code of the runs of record is commit 4ee1c6b; `rays` was added in
c883a16 and the final script's `check` reproduces the run of record to the printed digits (`final_check_none.log`).

## Coordinator's review

Review of branch `wip/i277-planted-family` at 5f2b81c and raw `I277-mac` at 63f0f279, recorded in negative-mass-raw `I277-mac-review`.

- **Scope.** The branch adds the script, its output and this note.
- **Rerun.** `planted_family.py check` on the Mac gives 11 of 11, every verdict line equal.
- **What it checks.** This is itself the independent check of #217's code: an exact nonlinear solution (#269), confirmed by the
  coordinator's brute-force Ricci check (`I269-mac-review/ind269.py`), reproduced at second order in h.
- **Findings, as reported.**
  - The `rate` mutation (K × 1.001) still solves the evolution equations, so only the constraints (D1, D3, D4) see it.
  - For β < 0, r → 0 at v* holds at fixed ρ only; along the grid's outgoing rays r tends to a positive constant. Paper IV's sentence now
    says "at fixed ρ", and the point bears on #278.
- **Applied.** `evolution.md` and Paper IV (Sec. 5 and the computations list) carry the result, and `paper_numbers.py` checks its numbers
  (139 of 139).

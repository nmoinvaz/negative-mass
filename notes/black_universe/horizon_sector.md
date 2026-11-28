# The horizon sector of the c = −1 black universe in Cauchy data (#248)

Code [horizon_sector.py](horizon_sector.py) (formulation, Kruskal-chart solver, functional, checks), output of record
[horizon_sector_output.txt](horizon_sector_output.txt) assembled from negative-mass-raw `I248-dh610/runs`; the criteria, written before the runs of record, are in
`I248-dh610/README.md`. Builds on #232 ([linear.md](linear.md)), #212 ([README.md](README.md)) and #217 ([evolution.md](evolution.md)). Status words as in AGENTS.md:
the algebra of the D checks is *checked symbolically*, the solver results are *numerical evidence*. No sentence here says a result is new.

## Answer

**c₂ is free in regular Cauchy data on a slice through the horizon.** It is the number `∂_V δr` at the point where the slice crosses the future horizon, the same
number at every point of the horizon down to the bifurcation surface S; no constraint restricts it, 20 of 20 random regular data have c₂ ≠ 0, and the growing solution
`Y = c₂ V G(Π)` (*update from #262:* in closed form `Y = c₂ V R(ρ(Π))`, R the areal radius; this note's `G(Π)` there is not dρ/dΠ) is smooth across S and through the interior (ρ < 0). So #232's horizon sector belongs to physical initial data: at c = −1 the linear problem has a
growing solution of rate κ, regular on the whole slice, excited by regular data that change δr at the horizon crossing, and not excited by data supported away from it.
The count of #212 (none) stands for the class `u = Y/ρ → 0`, which has c₁ = c₂ = 0. Status: the horizon law, the constraint propagation, the series at S and the gauge
statements are *checked symbolically*; the solver results (functional against #232's solutions to 3e-10 after Richardson, planted case, control, regularity, random data) are *numerical evidence*.

```
          Kruskal chart (V, U),  Pi = -V U,  rho(Pi) analytic,  rho = 0 on both horizons

                        interior (T-region)   V, U > 0,  rho < 0
                        ends at Pi = -0.37  (rho -> -oo)
          V = 0, U > 0     ╲                  ╱    U = 0, V > 0
          left future       ╲                ╱     future horizon H+ :  Y_h = c1 + c2 V
          horizon            ╲              ╱
                              ╲            ╱
     left exterior             ╲          ╱             right exterior
     V < 0 < U, rho > 0         ╲        ╱              U < 0 < V, rho > 0
              ─────────────────── S ─────────────────   slice T = 0  (V + U = 0)
                                ╱  bifurcation surface:  rho = 0, area minimum
                               ╱        ╲
                              ╱  past interior ╲
                             ╱    V, U < 0      ╲
```

## What c₂ is, as a functional of slice data

- **Chart.** `x⁺ = V = e^{bp}`, `x⁻ = U = −Q(σ)` of #232's chart (b = √2κ, κ = 2/π), `Π = −VU = e^{2κx}`, metric `−2e^{−F} dV dU + r² dΩ²`, `e^{−F} = G/κ`,
  `G = dρ/dΠ = A/(2κΠ)`, `r = R(ρ)`, `φ = √2 arctan ρ`. ρ(Π) is analytic through Π = 0, so the chart contains both horizons, S and the interior. #232's chart
  starts at p = 0 (V = 1) and does not contain S. The member solves the three evolution equations and both constraints in this chart, for symbolic κ (D1: 15 of 15
  residuals identically zero). The interior ends at Π ≈ −0.37; the solver rectangle V ∈ [−0.05, 3], U ∈ [−2, 0.06] reaches Π = −0.18 (ρ ≈ −1.5).
- **Gauge.** `Y = δr − (ρR/√2) δφ` is invariant under every linear diffeomorphism, not only the residual null reparametrizations: `L_ξ r = ξ^μ∂_μρ R′`,
  `L_ξ φ = ξ^μ∂_μρ Φ′`, and `R′ − (ρR/√2)Φ′ = 0` (D3a). On the future horizon of the c = −1 member ρ = 0, so `g = ρR/√2 = 0` and `Y_h = δr_h` (D2). The horizon is
  the area minimum, `R′(0) = 0`, so `r_U = r_V = 0` there and no shift of the horizon's position changes δr_h: c₂ is not a gauge artifact. At c = −1.5, where `R′(ρ_h) = 0.3216`,
  `r_U = −R′(ρ_h)G(ρ_h)V` on the horizon and a shift of the horizon produces exactly a `c₂V` term in δr_h (D3b); this is the planted case below.
- **The law.** On U = 0 the linearized constraint is `δr_VV = 0` (D2: `r_V = φ_V = F_V = 0` there), so `Y_h(V) = c₁ + c₂V` holds from S outward, and `c₂ = ∂_V δr` at any
  point of H⁺. The law is for δr and, at c = −1, for Y. A scalar field has no such law: a Schwarzschild test field's horizon values of δφ are not linear in V (N6), and at c ≠ −1 `Y_h = δr_h − g δφ_h` contains the free δφ_h.
- **The functional.** For the spacelike slice `V + U = const` crossing U = 0 at V₀, with tangent `T = ∂_V − ∂_U` and normal direction `N = ∂_V + ∂_U`:

  ```
  c₂ = ∂_V δr (V₀, 0) = ½ (∂_T δr + ∂_N δr),      c₁ = δr(V₀, 0) − V₀ c₂
  ```

  `∂_T` is the derivative of the slice values, `∂_N` the normal derivative (the second Cauchy datum). The code uses 4th-order central stencils on the nodes
  `(j₀ + k, i₀ − k)` and `(j₀ + k, i₀ + k)`, k = ±1, ±2. Through S the slice is V + U = 0 (T = 0). Invariantly, `c₂ = ℓ^μ∇_μ Y` with ℓ = ∂_V the horizon generator.
- **What c₁ is.** `c₁ = Y(S) = δr(S)` is the first-order change of the area of the bifurcation surface, δA_S = 8π r_S c₁ (r_S = 1 at c = −1), and κ c₁ = δM_ADM (#259, [static_c1.md](static_c1.md)). It is the coefficient of the static mass mode, the exact solution r = √L R, F = A + 1/L − 1 of #269 (c₁ = δL/2, c₂ = 0, δM_ADM = δL/π), and it does not grow. The classes A0, A1R, A1L carry mass; G, L, R do not.

## The check on #232's solutions (acceptance 1)

#232's solver at c = −1 (σ* = 16, left-ray pulses at σ₀ = 2.5, 3.5, 4.5 of check L3a, P = 10), h = 0.005 and 0.0025, Richardson. c₂ is the coefficient of e^{bp} = V. The
functional is evaluated on the slices `p + σ = const` crossing the horizon at p = 1, 2, 4, 6 and compared with the least-squares fit of the horizon values and with #232's own
two-point c₂ (`c12_from`):

| σ₀ | c₂ | slices p = 1, 2, 4, 6 against the fit | against #232's two-point c₂ | raw at h = 0.005 |
|---|---|---|---|---|
| 2.5 | −1.987994719e-04 | 3.3e-10, 3.4e-10, 3.4e-10, 3.5e-10 | 3.4e-10 | 2.2e-06 |
| 3.5 | −1.040481039e-04 | 3.1e-10, 3.2e-10, 3.2e-10, 3.3e-10 | 3.2e-10 | 3.0e-06 |
| 4.5 | −4.723353344e-05 | 3.0e-10, 3.1e-10, 3.1e-10, 3.1e-10 | 3.0e-10 | 3.5e-06 |

The bar of 1e-6 is met after Richardson by three orders; the raw value at h = 0.005 alone (2.2e-6 to 3.5e-6) is not below 1e-6. The same c₂ at four crossings is the e^{bp} law
(N5). In the Kruskal solver through S the functional at V₀ = 0 (through S), 0.5 and 1 equals the horizon fit to ≤ 7e-8 (N4).

## Regularity across S and into the interior

| check | what | result |
|---|---|---|
| D6 | formal series of the Y equation at S, degree ≤ 7 (c = −1, M = 0): free coefficients exactly `y_mn` with m = 1 or n = 1 (13), consistent; `y_00` and `y_0n`, `y_m0` (n, m ≠ 1) vanish | 32 equations, 23 determined unknowns, 13 forced to zero |
| D7 | Frobenius series of Y at W = κ, exponent 1/2: no logarithm, `g₂ = π²/128 + 5/12`; so `Y_ref e^{−κx}` is analytic in ρ and in Π | solvable at every order ≤ 5 |
| D8 | no singular point on ρ ∈ [−200, 0) ∪ (0, 200]: A < 0 for ρ < 0, A > 0 for ρ > 0; the ρ^{1/2} solution integrates to ρ = ±40 | finite; max \|Y\| 6.3e12 at ρ = 40 (static), 56.6 (interior) |
| N2 | exact pure-gauge solutions reproduced through S and the interior, c = −1, −1.5, Schwarzschild | last ratio 3.98, 4.00, 4.00 |
| N8 | δr, δF, δφ and Y at S, two interior points (−Π = 0.04, 0.12), the static region and the left side, classes G and A0 | Cauchy ratios 4.00, 4.00 |

D6 says that the smooth solutions of the Y equation near S contain the two non-smooth-`u` branches `Y = y₁₀V` and `Y = y₀₁U` (u = Y/ρ with a simple pole at U = 0 or V = 0) and nothing else
outside the smooth-u class `Y = O(UV)`. `Y = V G(Π)` (= V R(ρ(Π)), #262) is the W = κ branch: smooth at S, equal to V on H⁺ exactly, continued analytically into Π < 0. In the exterior it grows as `e^{κt}` at fixed ρ
(`V = e^{κ(t+x)}`) and toward large ρ, so it is not an L² bound state in the tortoise coordinate: that is why #212's count misses it.

## The data are free, and physical

- **Free.** Propagation: `∂_U δC_VV = −(r_U/r) δC_VV` and `∂_V δC_UU = −(r_V/r) δC_UU` on solutions of the evolution equations, for any background solving the field equations (D4, exact rational
  residual 0 at generic jets). So the two constraints restrict second derivatives only (`a_VV`, `a_UU`), and the first-order jet `(δr, ∂_Vδr, ∂_Uδr, …)` at a point of H⁺ is unrestricted:
  ∂_Vδr = c₂ is free.
- **Constraint-satisfying.** The evolved solutions satisfy C_VV and C_UU at truncation level: ratios 4.00 and 4.01 (G), 3.95 and 4.14 (A0) per halving, below 7.5e-4 of the largest term (N10).
- **Random data.** 20 random smooth data (corner constants and pulses, amplitudes ≤ 1e-3): |c₂| between 3.4e-4 and 2.1e-2, median 7.7e-3; the slice functional equals the fit to 1.4e-8 (N11).
- **Classes** (c = −1, Kruskal chart through S, amplitudes 1e-3, Richardson of h = 0.0025, 0.00125; the fit and the slices V = 0 (through S), 0.5, 1 agree to the digits shown):

  | class | data | c₁ = Y(S) | c₂ |
  |---|---|---|---|
  | L | δφ pulse on the ray V = −0.05, centred on the horizon; zero on the ray U = −2 | −9.3e-13 | −9.610862e-04 |
  | G | L and a pulse on the ray U = −2 at V = 1 | −4.8e-12 | −9.610862e-04 |
  | A0 | corner δr | 1.658047e-03 | 1.227043e-02 |
  | A1R | corner δr_V | 5.236768e-05 | 1.047354e-03 |
  | A1L | corner δr_U | 2.094707e-03 | 1.007157e-02 |
  | R | pulse on the ray U = −2 only (exterior, away from the horizon crossing) | −3.9e-12 | 5e-13 (seed 1.8e-9 at h = 0.00125, falling by 4 per halving) |

  Class L carries no data on the exterior ray (δr = δφ = 0 on U = −2 for all V), and has c₂ = −9.6e-4. Class R gives c₂ = 0 up to the truncation seed: c₂ depends only on data in the
  past of the crossing, and the same number holds all along H⁺, so data supported away from the crossing and from S do not excite it (*update from #258:* this holds for one-sided data; on a full slice the constraints also tie c₂ to the data on the other side of S, and two-sided data supported away from S have c₂ ≠ 0 in general, −3.18e-4 at amplitude 1e-3, because the constraint solution is non-zero between the supports). This is the unperturbed-horizon class of #232.

## Planted case and control

- **Planted (N6).** Pure gauge `L_ξ(member)` at c = −1.5 has `δr_h = −η(0)R′(ρ_h)G(ρ_h)V` with G(ρ_h) from a 60-digit quadrature (known coefficient −1.12925016e-02 for η(0) = 0.03), and in Schwarzschild (r_s = 1, φ = 0, V = 0)
  `δr_h = −η(0)e^{−1}V`. The functional returns the known coefficient to 3.6e-9 (c = −1.5) and 1.2e-10 (Schwarzschild) on exact data, and 1.4e-8 and 4.8e-10 on the marched solution; the
  gauge-invariant version `∂_V Y` returns 0 (5e-9 and 2e-8 of the coefficient at c = −1.5). A Schwarzschild test field (a pulse on the right ray) leaves δr = 0 exactly, and its horizon values of δφ obey
  no linear law (fit residual 0.79 of the maximum at h = 0.01 and 0.005, not falling), against 4e-6 falling by 4 for Y_h at c = −1. So the extraction recovers a coefficient where no growing sector exists, and the linear law is specific to Y at c = −1.
- **Control (N7).** Exact pure gauge at c = −1 (r stationary on the horizon): |c₂|/max|δr| = 2.7e-6 at h = 0.005, falling by 16.1 from 0.01 (the stencil's h⁴). Class R: seeds 1.2e-7, 2.9e-8, 7.3e-9, 1.8e-9 (ratios 3.96, 3.99, 4.00),
  1.9e-6 of class L's c₂. Superposition `L + 0.917633 A1R` with λ from the horizon fits: c₂ = 1.5e-13 of class L's by the fit and 9.0e-10, 2.8e-9, 4.8e-8 by the slice functional at V = 0, 0.5, 1.

## Checks (D: `horizon_sector.py derive`, 8 of 8; N: `horizon_sector.py check`, 11 of 11)

| id | what | result |
|---|---|---|
| D1 | the Kruskal-chart member solves the 3 evolution equations and 2 constraints, symbolic κ, c = −1, −3/2, −3/10 | 0 of 15 residuals nonzero |
| D2 | on U = 0, C_VV is `a_VV = 0`; g = 0 there | exact |
| D3a, D3b | Y invariant under every diffeomorphism; r_U on H⁺ zero at c = −1, nonzero at −1.5 | exact; R′(ρ_h) = 0.3216 |
| D4 | constraint propagation, generic jets | residual exactly 0 |
| D6, D7, D8 | series at S, no logarithm at W = κ, no singular point | see above |
| N1 | Kruskal jets against `linear.bg` | 5.8e-16 |
| N2 | pure gauge through S, three backgrounds | last ratios 3.98, 4.00, 4.00 |
| N3 | horizon values linear in V on [−0.05, 3], classes L, G, A0, A1R, A1L; residual at h = 0.00125 | 4.7e-7, 3.4e-6, 4.3e-6, 4.1e-6, 2.2e-6; ratios 3.99 to 4.00 |
| N4 | slice functional against the fit, V₀ = 0, 0.5, 1 | ≤ 6.8e-8 |
| N5 | #232's solutions | 3.5e-10 after Richardson |
| N6, N7 | planted and control | above |
| N8 | convergence at S, interior, left | 4.00, 4.00 |
| N9 | Kruskal solver against `linear.py` on corner-constant data at (V, U) = (1, −1) | 1.3e-6, 1.8e-6, 4.6e-7 (bar 1e-5) |
| N10 | constraints of the evolved solutions | ratios 3.95 to 4.14 |
| N11 | 20 random data | slice against fit 1.4e-8; \|c₂\| ≥ 3.4e-4 |

Mutations: 7 for `derive` (d_emF, d_Gode, d_hz, d_prop, d_gauge, d_res, d_log) and 9 for `check` (k_emF, k_FV, k_cross, k_con, fn_aT, fn_dirU, fn_b, plant_sign, ctl_lambda), each with the checks it must fire declared
before the runs of record (README of the raw folder); in the runs of record the fired sets equal the declared ones (16 of 16; `horizon_sector_output.txt`). Two mutations of the solver (`k_cross`, `k_con`) do not fire N8 or N5: N8 tests convergence, not the equation,
and N5 uses `linear.py`'s solver.

## What this means for the c = −1 member, in the terms of #232's correction

#232 corrected #212's "no growing mode at c = −1" to "none for Y_h = 0" and left open whether Y_h ≠ 0 is a class of physical data. It is. `Y_h = 0` means `c₁ = c₂ = 0`: two numbers, `Y(S)` (static) and
`∂_V δr` at the horizon (growing). c₂ is a free number of regular Cauchy data, so:

- the linear problem about the c = −1 member has one growing direction of rate κ (#232's rate, numerical evidence there), a solution that is regular across S and through the interior, and is excited by regular data
  with c₂ ≠ 0, such as a δφ pulse crossing the horizon (class L, c₂ = −9.6e-4 with amplitude 1e-3 and no data on the exterior ray);
- the rest of the data space, c₂ = 0, is a hyperplane at the horizon crossing; it contains the data supported away from the horizon crossing (class R), which decay, as #212 and #232's class R;
- #212's L² count misses the growing direction because `Y ∝ V` is `u ∼ 1/U`, not normalizable in the tortoise coordinate; the horizon sector adds one growing direction, W = κ (the logarithm of the Frobenius recursion vanishes only there, #232 D3 and D7 here), so c = −1 then has one like the other 19 members.

Nothing is said here about the nonlinear fate: #217's lifetime law with slope 1/κ at c = −1 is the nonlinear observation consistent with this (read from #232), not tested here.

## What README.md and critical_ghost.md §1.1 should now say (proposed; the files are not touched here)

`README.md`, Result, first bullet, replace the sentence "Update from #232: …": *Update from #232 and #248: with Y regular and nonzero on the horizon the linear problem at c = −1 has a growing solution of rate W = κ
(`Y ∝ e^{κv}`, regular across the bifurcation surface and the interior). Its coefficient c₂ = ∂_V δr on the horizon is a free number of regular Cauchy data on a slice through the horizon (checked symbolically; numerical evidence); #212's count of zero is
the count for `u = Y/ρ → 0` (c₁ = c₂ = 0), which also holds for data supported away from the horizon crossing.* Sec. 8, knife edge, last sentence: *The rate is κ at c = −1 as well; the member decays only for data with c₂ = 0, a single condition on the data at the horizon crossing, which includes data supported away from the horizon (#248).*

`critical_ghost.md` §1.1, the black-universe row, replace #232's sentence: *At c = −1 the count of #212 is zero for perturbations with `u = Y/ρ → 0`; the linear problem has a growing solution of rate κ whose coefficient c₂ = ∂_V δr on the horizon is a free number of regular Cauchy data
(#248; checked symbolically and numerical evidence). The member decays for data with c₂ = 0, which include data supported away from the horizon.*

## What this does not show

- A full Cauchy slice to spatial infinity: *built in #258* ([full_slice.md](full_slice.md)). Compactly supported δφ and ∂_Nδφ on the slice V + U = 0 through S
  solve both constraints with a solution that is exactly zero outside the support (δM_ADM = 0), and c₂ from the slice equals c₂ of the evolved horizon fit to 2.7e-8 (numerical evidence).
- The left exterior beyond a sliver (V ≥ −0.05, U ≤ 0.06); the interior past −Π = 0.18 (the analytic solution is followed to ρ = −40, D8, but I⁺ was not examined).
- The late-time growth: the Kruskal chart cannot reach V ∼ e^{bP}. The rate κ is #232's; here only the horizon law and its coefficient are tested. The decaying partner `Y ∝ U` (the coefficient ∂_U Y(S) on the past horizon) is not computed separately.
- *Done in #259:* c₁ is the static mass mode, not a neighbouring member δc (the potential of the member changes with c); κ c₁ = δM_ADM to 7.4e-8 for A0, A1R, A1L.
- N2 at c = −1 has a first ratio of 5.13 (h = 0.02 to 0.01, pre-asymptotic); the bar is on the last ratio, 3.98. N5's raw values at h = 0.005 are 2.2e-6 to 3.5e-6; the 1e-6 bar is met by Richardson.
- Radial perturbations of D = 4, one potential; nothing about polar modes (#218) or the nonlinear fate.

## Recommended tasks

1. **A compactly supported datum on a full slice.** *What:* solve the two constraints along a spacelike slice from the left through S to large ρ, with δφ and ∂_N δφ compactly supported near the crossing and δM_ADM = 0, and evolve. *Why:* closes the infinity gap above.
   *Acceptance:* the constraint solution decays to pure gauge outside the support; c₂ from the slice equals c₂ of the evolved horizon fit to 1e-6; a mutation with the wrong asymptotic condition changes c₂. *Depends on:* #248.
2. **The count with the Kruskal-regular class.** *What:* #249's recount with the W = κ solution `Y ∝ V G(Π)` admitted at c = −1 and the c ≠ −1 members unchanged. *Why:* the horizon sector is now a datum, not an option. *Acceptance:* as #249's, plus c₂ from this note's functional for the eigenfunction. *Depends on:* #249.
3. *Done in #259.* **The static sector.** *What:* identify c₁ (data A0, A1L) with the static M and with a static perturbation of the member. *Why:* c₁ ≠ 0 is regular data; its role in #217's family is open. *Acceptance:* c₁ as a functional of the slice data, checked on #232's solutions to 1e-6. *Depends on:* #232, #248.
4. **#217's family along c₂.** *What:* #250's lifetime law along class R and a family with c₂ tuned, now that c₂ is a datum. *Acceptance:* #250's. *Depends on:* #250.

## Reproduce

From this directory, `~/Source/negative-mass-venv314/bin/python`, one thread, about 11 minutes in all:

```
horizon_sector.py derive                          # 8 checks, 1.5 min:  --mutate=d_emF|d_Gode|d_hz|d_prop|d_gauge|d_res|d_log
horizon_sector.py check [--only=N4,N5]            # 11 checks, 9 min:   --mutate=k_emF|k_FV|k_cross|k_con|fn_aT|fn_dirU|fn_b|plant_sign|ctl_lambda
horizon_sector.py table                           # the class table and the samples, 1 min
horizon_sector.py mutations DIR PREFIX            # fired sets of stored logs DIR/PREFIX<name>.log against the declared ones
```

`horizon_sector_output.txt` is `derive`, `check`, `table` and `mutations` of the runs of record (`I248-dh610/runs`, code commit 69de2ea), in that order. Memory 718 MB at most (N5, h = 0.0025).

## Coordinator's review

Review of branch `wip/i248-horizon-sector` at 0dd8c98 and raw `I248-dh610` at 88ac882, recorded in negative-mass-raw `I248-mac-review`
(criteria committed before the first run; results 120e0775). On the Mac, nice 10.

- **Scope.** The branch adds the script, its output and this note, nothing else.
- **Rerun.** `derive` gives 8 of 8 and `check` 11 of 11. `table` gives class L's c₂ = −9.610862e-04 at the fit and on the three slices.
  The verdict lines equal dh610's, apart from roundoff-level digits (N5's slice deviations of 3.3e-10 against 3.4e-10).
- **Independent check** (`ind248.py`). It computes c₂ of class L by transport along the data ray V = −0.05, which runs inside the horizon
  (Π from −0.1 to 0). It uses no 2D solver and none of this note's or #232's code.
  - The transverse derivatives (δr_V, δF_V, δφ_V) obey the linearized evolution equations along the ray, in the coordinator's forms
    verified in `ind217.py`, with the background in mpmath at 30 digits. c₂ = δr_V at U = 0.
  - Result: c₂ = −9.61086221626e-04, which agrees with the 2D Kruskal solver to 2.3e-8, the note's printed digits.
  - The same transport reproduces an exact pure-gauge solution along the whole ray to 3.3e-14, and gives c₂ = 0 for it, since r_UV = 0
    on U = 0 at c = −1.
  - Mutations fire as declared: the scalar's sign (K1, K2), the potential dropped (K1, K2), and the pulse centred at U = 0.01 (K1 alone).
  - Shared with the code under test: the member, the chart's normalization Π(ρ = 1) = 1, class L's definition and the gauge δF = 0 on the
    data rays.
- **What this settles, and what not.** That c₂ is set by data at and before the horizon crossing, and is nonzero for a generic pulse
  there, now rests on two methods. The free-data argument (constraint propagation, D4) and the regularity at S (D6, D7) are *checked
  symbolically* here and were not redone. The open gap is the note's own: no compactly supported datum on a full slice to spatial
  infinity, with δM_ADM = 0, was built. That is recommended task 1, filed below. The nonlinear fate is not addressed.
- **Applied.** The black-universe README (Result, §5, §8) and `critical_ghost.md` §1.1 carry the update. The row's lead no longer calls
  c = −1 an exception.
- **Follow-ups.**
  - Task 1 is #258.
  - Task 2 is added to #249.
  - Task 3 is #259.
  - Task 4 is added to #250.
- *Update from #262 (coordinator):* the Frobenius solution of D7 with g₁ = 0 is Y_ref − (π/8) Y_b. Y_ref = R e^{κx} itself has g₁ = π/8 and
  g₂ = 5/12 + 3π²/128. The Kruskal-regular solutions at W = κ form a two-dimensional space {Y_ref + βY_b} with the same c₂, and the condition at
  spatial infinity fixes β (`cminus1_mode.md`). #262's normalization has e^{2κx} = ρ(1 + …) at 0, and this note's has Π(ρ = 1) = 1, so they differ by
  Π(1) = 2.4725459594.

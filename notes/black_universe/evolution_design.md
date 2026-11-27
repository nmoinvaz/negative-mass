# Evolving a black universe along its growing mode: design (#217)

Design note, written before the build. Equations: [evolve_derive.py](evolve_derive.py), output [evolve_derive_output.txt](evolve_derive_output.txt). Potential and
member: [evolve_potential.py](evolve_potential.py). Prototype, not the build: [evolve_dev.py](evolve_dev.py); its runs [evolve_dev_runs.py](evolve_dev_runs.py) and
[evolve_dev_runs2.py](evolve_dev_runs2.py); the linear master equation in the time domain [evolve_master1d.py](evolve_master1d.py). Logs: negative-mass-raw `I217-mac/dev`
(`dev_NN`). Status words as in AGENTS.md. Every number from a development run is *numerical evidence*. No criterion is declared yet; no run of record exists.
**Built and run: see [evolution.md](evolution.md).** **After the coordinator's rulings (c89cb61):** Sec. 12 reports the diagnosis of rulings 4 and 5. It supersedes the c = −1.5 numbers and the departure threshold of Secs 5 to 7 (the old values stay in the logs `dev_40` to `dev_49`).

## Rulings requested (answered by the coordinator; the findings are in Sec. 12)

1. **Family.** Member + A × δ̂, with δ̂ the growing mode's shape in the code's own gauge, taken from the member's own truncation-seeded departure in the linear regime (Sec. 5).
   Not a transformation from #212's gauge. The one-sided pulse of #40 (family O) stays as a cross-check.
2. **Kernel.** A new kernel in `evolve.py` beside `dn.py`; `dn.py` stays unchanged (Sec. 3).
3. **Lifetime and finders.** Lifetime = the departure time of the gauge-invariant Y at the diagonal. Outcomes are named by the finders, with the raw test of `dn.cell_class`
   weighted by e^{f/2}, because the raw one fails in the tortoise gauge (Sec. 6).
4. **c = −1 control.** In the double-null prototype it does **not** behave as #212 says (a threshold with slope 1/κ). An independent time-domain run of #212's master
   equation does behave as #212 says (Sec. 8). I ask for a time-boxed diagnosis in a horizon-regular gauge before criteria are declared, and for the control to be reported as found.
5. **c = −1.5.** The lifetime law is clean only for |A − A\*| above about 3e-6, where T saturates near 15 whatever h is; A\*ₕ is not ∝ h²; and the slope has a 6 to 8 % fit error. The same diagnosis may remove all three (Sec. 7).
   Say whether to spend it before declaring, or to declare windows above the floor and report the member as found.

## 1. Equations

`R_μν − ε ∂_μφ ∂_νφ − V(φ) g_μν = 0`, `ε □φ = V′(φ)`, ε = −1, 8πG = 1, D = 4, `ds² = −2e^{−f}dx⁺dx⁻ + r²dΩ₂²` (notation of #39: p = x⁺, m = x⁻).
Derived from the brute-force Ricci tensor of the full four-dimensional metric and the Laplacian of φ, with a generic V (`evolve_derive.py`):

```
r_{+-}  = −(r₊r₋ + e^{−f}/2)/r + ½ V r e^{−f}
f_{+-}  = (2/r) r_{+-} + ε φ₊φ₋ − V e^{−f}          [ = −2(r₊r₋ + e^{−f}/2)/r² + ε φ₊φ₋ : the V terms cancel ]
φ_{+-}  = −(r₊φ₋ + r₋φ₊)/r − (ε/2) e^{−f} V′(φ)
C₊₊ = r₊₊ + f₊r₊ + (r/2) ε φ₊² = 0,   C₋₋ likewise           [ no V: g₊₊ = g₋₋ = 0 ]
```

| check | what | result |
|---|---|---|
| A1 | brute-force equations equal the three forms above, V general, ε = ±1 | pass |
| A2 | the two constraints, no V term | pass |
| A3 | V = 0 equals `dn.F_np` at 40 random states | 4.4e-16 |
| B1 | ∂₋C₊₊ = −(r₋/r) C₊₊ on solutions, V general (the potential drops out; this fixes the sign of the V′ term against the Einstein equations) | pass |
| C1 | on the static ansatz the evolution residuals are r_{+-} = (A/2R)·[angular eq.], C₊₊ = (A²R/4)·(5), φ_{+-} = −(εA/2R²)·(3), f_{+-} = −(A/2R²)·(6) + (A²/2)·(5), of #212 | pass |
| C2 | the family of [BF06], c = −1.5, −1, −0.3, 0: three evolution equations and both constraints, 40 digits | 1.8e-40 |

Mutations `V_sign_r`, `V_coef_r`, `V_sign_f`, `Vp_sign` (each fires A1, B1, C1, C2) and `eps_sign` (A1, A2, A3, C1, C2) fire exactly the declared checks
(`dev_01` to `dev_06`; the targets were computed in those runs before they were written into the script). n = 2 only: for general D the V terms carry 2/n, not derived here.

## 2. Static member in the tortoise null gauge

`x = (x⁺ − x⁻)/√2`, `t = (x⁺ + x⁻)/√2`, `dρ/dx = A(ρ)`: **e^{−f} = A, r = R, φ = √2 arctan ρ**, r_± = ±R′A/√2. At c = 0 this is #39's Ellis data: the rays agree with
`dn.ellis_rays` to 1.7e-9 in r (its RK4 constraint solve, #39's I1), 5e-13 in φ, 0 in f (`dev_31`). In this gauge r₊r₋ = −(R′A)²/2 ≤ 0 in the whole static region, so no
trapped or anti-trapped cell exists in the grid for the static member.

## 3. Where the potential term goes

| option | cost | risk |
|---|---|---|
| flag in `dn.py`: extra entries in `par`, a branch in `_F` and `F_np` | about 12 lines | touches the file #39, #40 and their reviews rest on; every #39 and #40 result has to be reproduced again |
| **new kernel in `evolve.py`**: `_F`, `row_step`, `evolve` copied with the V term, importing `dn.Data`, `cell_class`, `strip_cells` | about 150 duplicated lines | drift between the copies |

Choice: the new kernel. The drift is caught by a check that can fail: at c = 0 (V ≡ 0 exactly) with the potential flag off or on, the prototype's G array is **bit-identical** to
`dn.evolve`'s on a 300 × 300 grid (max |ΔG| = 0, same dead mask; `dev_07`). The build repeats it on three grids with both kernels, and adds a mutation that applies the V term with a wrong sign at c = 0.
The potential is in `evolve_potential.py`: `pot(phi, c, K) → (V, V′)` (numba-compatible) and `Member(c)`; **#222 can import it** with no dependence on `dn.py`.
`dV/dψ = −2c sin 2ψ − 2K(1 + cos 2ψ + ψ sin 2ψ)`, ψ = φ/√2; V(ψ = π/2) = 0 and V(−π/2) = −6c.

## 4. Data for each member

| member | ρ_h | κ | W | corner | rays |
|---|---|---|---|---|---|
| c = −0.3 | −0.926754 | 0.6013 | 0.97008 | the throat ρ = 0, on the diagonal l = 0 | right: ρ from 0 to ∞; left: ρ from 0 toward ρ_h |
| c = −1.5 | 0.339617 | 0.5517 | 0.49417 | ρ_c = 1, in the static region, outside the horizon | right: ρ from 1 to ∞; left: ρ from 1 toward ρ_h |
| c = −1 (control) | 0 | 0.6366 | none (#212) | ρ_c = 1 | as c = −1.5 |

- **Construction.** ρ(x) comes from `dρ/dx = A` (DOP853, rtol 1e-13). On the left ray the variable is ζ = ln(ρ − ρ_h), `dζ/dx = A/(ρ − ρ_h)`, with `A/z` from a Taylor series for z < 1e-3,
  so A and f = −ζ − ln(A/z) keep full relative accuracy. r, f, φ on each ray are exact closed forms of ρ. No constraint solve is needed for the member.
- **The horizon is not reached.** In this gauge it sits at l → −∞. At the end of the left ray (X = 40 and 50), z = ρ − ρ_h = 1.6e-15 (c = −0.3, l = −28.3) and 7e-18 (c = −1.5, l = −35.4), f ≈ 33 and 39.
  A perturbation that makes the trapped region grow outward shows as trapped cells at finite coordinates (Sec. 6); the finders never see the static member's own horizon.
- **Rectangle.** X = 40 for c = −0.3 and X = 50 for c = −1.5 and −1: the latest time is √2 X = 57 and 71, above the longest departure (about 25) and above the post-departure naming runs.
  No outer boundary. Memory is O(N) with `store=False` (the dev runs with `store=True` at N = 2400 peak at 0.9 GB).
- **Corner for c = −1.5.** ρ_c = 1 and 2 give slopes that agree within the fit error (Sec. 7); ρ_c = 0.6 failed its bracket (`dev_47`, kept).

## 5. The family

**Chosen: static member + A × δ̂.** Recipe (`evolve_dev.extract_mode`, prototype):
1. Run the untuned member. The truncation error seeds the unstable mode.
2. At the first row (after t = 6) where |Y| = 2e-4 at the diagonal (2e-3 in `dev_40` to `dev_49`; Sec. 12), take the data on the two rays through that diagonal corner and subtract the exact static member's.
   The static member's rays through any l = 0 corner are the same, since the member is static.
3. Divide by Y at the corner node, so that **A is Y at the corner**. Data = static + A × δ̂ on each ray.

`Y = δr − (R′/φ₀′) δφ = δr − (ρR/√2) δφ` is #212's δφ = 0 variable and is invariant under the residual relabelings of the double-null gauge to first order. In the extraction run the
mode is linear (Y = 2e-3 at the corner node). The extraction times are t_s = 8.1, 9.6, 11.0 at c = −0.3 and 13.4, 16.2, 19.0 at c = −1.5 for h = 0.02, 0.01, 0.005.
δ̂ is mostly gauge: per unit Y its largest |δr|, |δf|, |δφ| are 3.8, 7.3, 1.1 (right ray) and 1.0, 6.7, 3.1 (left ray) at c = −0.3, h = 0.02; and 9.3, 10.8, 2.3 and 5.2, 20.1, 13.2 at c = −1.5, h = 0.01.
Data of this form violate the nonlinear constraint at O(A²); the build measures it (a check).

*Why not the one-sided pulse (#40's O).* It has a position, a width and a transient. At c = −1.5 its own Y at the diagonal at t = 12 is 4.0e-3 at the expansion boundary (h = 0.02), above the departure threshold
of 3e-3; and the pulse family has two boundaries at O(h²) from each other, the expansion event at A\*_E = −1.9077e-4 and the sign change of Y at A = +6.0e-5 (h = 0.02; −5.742e-5 and +1.6e-5 at h = 0.01; `dev_52`).
The extracted family has the sign change only. Both pass through the static member, so both have A\* = O(h²). The pulse is kept as the cross-check that the slope does not depend on the family.

*What the extraction shares with the code:* all of it. Independent checks of the mode, all in the build: (i) W from the code against #212's table; (ii) the shape Y(l) of the extracted mode against #212's Y(ρ(l)),
shooting from the horizon and from infinity; (iii) the 1D master equation (Sec. 8). The transformation from #212's gauge to ray data (Y(ρ) on the whole static region with a complex-arc continuation
through the throat at c = −0.3, δν from the scalar equation, two gauge quadratures per ray) is not built; it is a possible task.

**A\*ₕ is the amplitude of δ̂ that cancels the truncation seed** (sign of the departure on the diagonal, bisection to 1e-14; linear-level settings of Sec. 12: extraction at |Y| = 2e-4, departure at |Y| = 2e-5; `dev_80` to `dev_88`):

| member | h = 0.02 | 0.01 | 0.005 | ratios |
|---|---|---|---|---|
| c = −0.3 | −7.89399e-07 | −1.92364e-07 | −4.75546e-08 | 4.104, 4.045 |
| c = −1.5, ρ_c = 1 | −2.85416e-06 | −7.04113e-07 | −1.77257e-07 | 4.053, 3.972 |
| c = −1 | −3.39873e-06 | −8.29464e-07 | −2.05469e-07 | 4.098, 4.037 |

Second order at all three members; Richardson (order 2) gives A\* = −1.6e-9 at c = −1.5 and +8e-10 at c = −0.3, zero to the extrapolation error. With the threshold of the first design (|Y| = 3e-3) the c = −1.5 values were
−6.64e-06, −3.60e-06, −4.28e-06 (`dev_43` to `dev_45`): not second order, because |Y| = 3e-3 is already nonlinear there (Sec. 12).

## 6. Outcomes, finders, lifetime

- **The raw finder is gauge dependent.** `dn.cell_class` tests r₊ and r₋ against an absolute 1e-2. In the tortoise gauge f reaches 3 to 20 near threshold, so r₊, r₋ fall below it on solutions that are expanding.
  Weighted test: the signs of r₊ and r₋ and `2 e^f r₊r₋ ≥ 0.02` (2M/r − 1; identical to the raw test at f ≈ 0). At c = −0.3, h = 0.02, A = A\* + 1e-8 the raw finder is silent through the whole run
  and the weighted one fires at t = 17.3 (`dev_50`).
- **c = −0.3** (mode family, `dev_50`, weighted finders, |A − A\*| from 1e-3 to 1e-8): A < A\* forms a trapped region (C, collapse; first cell at l = +0.06 to +0.17, t = 4.7 → 16.7).
  A > A\* forms an anti-trapped region (E, expansion; first cell at l = −0.28 to −0.14, t = 5.3 → 17.3). Both outcomes are found and named by the finders over 5 decades. On the E side the diagonal cell is not in the
  anti-trapped region later: at A\* + 1e-8 it has r = 15.9, r₋ = −8e-10 and 2M/r = 0.12 at t = 50.
- **c = −1.5** (`dev_51`, h = 0.01): A < A\* collapses at every |A − A\*| from 1e-3 to 1e-8 (trapped cells at l = −0.2 to +1.3, t = 8 → 19). For A > A\* and |A − A\*| ≥ 1e-5 there is no trapped cell, and the weighted
  finder fires late (t = 50 to 63, l = −2.2 … +4.2); the raw finder is silent. For |A − A\*| ≤ 1e-6 **both sides collapse at t ≈ 18 to 19**: that is the floor of Sec. 7.
- **Lifetime.** The departure of Y at the diagonal: the first row with |Y| ≥ Y₀, fate = sign of Y there, T = t of that cell. Y₀ = 3e-3 was the ruling; Sec. 12 finds it nonlinear at c = −1.5 and proposes Y₀ = 2e-5. It is a linear-regime measure, defined on both sides at every |A − A\*|, and independent
  of what the nonlinear stage does. At c = −0.3 the weighted finder event follows the departure by a constant 3.6 (collapse) and 4.1 to 4.3 (expansion) over the five decades (`dev_50`), so the slope is the same. The finder times are
  reported as a second series where they exist.
- **Naming runs** (build): for each member, both sides, |A − A\*| = 1e-4, 1e-5, 1e-6, run to the end of the grid: first trapped and anti-trapped cell (raw and weighted), r on the diagonal against t, the Misner–Sharp mass
  against the member's ADM mass m = −2c/(3π) (0.0637 and 0.318), the signs of r₊ and r₋. If a side stays without a finder event it is reported so.

## 7. Static check, growth seed, fit

**Planted static check** (`dev_08`, X = 10, h = 0.04, 0.02, 0.01; max over r, f, φ against the exact member, all four members):

| c | error at t ≤ 8, h = 0.04, 0.02, 0.01 | ratios |
|---|---|---|
| 0 | 1.6e-1, 3.7e-2, 9.2e-3 | 4.31, 4.05 |
| −0.3 | 2.1e-2, 5.4e-3, 1.4e-3 | 3.92, 3.97 |
| −1.5 | 2.0e-3, 5.1e-4, 1.3e-4 | 3.96, 3.98 |
| −1 | 7.6e-3, 1.9e-3, 4.8e-4 | 3.95, 3.98 |

The declared window will be t ≤ 8, orders 1.75 to 2.25. The truncation seed is A\*ₕ above. After the window the error grows at the member's rate: `dev_09`, the gauge-invariant Y from the truncation seed alone,
h = 0.01, fit over t = 8 to 14 at l = −2, −0.5, +1:

| c | Y rate | #212 W | κ |
|---|---|---|---|
| −0.3 | 0.963, 0.951, 0.948 | 0.97008 | 0.6013 |
| −1.5 | 0.494, 0.491, 0.490 | 0.49417 | 0.5517 |
| −1.01 | 0.633, 0.631, 0.629 | 0.63315 | 0.6351 |
| −0.99 | 0.640, 0.638, 0.636 | 0.64012 | 0.6381 |
| **−1** | **0.637, 0.634, 0.632** | none | 0.6366 |

The raw errors in r, f, φ grow faster (0.57 to 0.75 at c = −1.5 and −1): they carry the gauge part.

**Fit** (as #40): bisection of the sign of the departure to hi − lo ≤ 1e-14, then A = A\* ∓ 10^−e, e = 2.5 … 11.5 in steps of 0.5; T = c_side + s ln(1/|A − A\*|), common slope, one constant per side;
five windows × (joint, each side) = 15 fits; three error bars named as in #40: the OLS SE; the envelope of s ± 2 SE over the fits (fit spread); the difference of the two finest h (bar 0.5 %).
**Planted series** through the same `bisect`, `series`, `fit_set` as `threshold.py` (40 noise seeds, slope 1/W recovered, 1/(√2 W) rejected, as there), with the time-unit mutation.

**Development result** (linear-level settings of Sec. 12, `dev_80` to `dev_88`; slope against 1/W (1/κ at c = −1), window |A − A\*| ∈ [1e-7, 1e-5], joint fit):

| member | h = 0.02 | 0.01 | 0.005 |
|---|---|---|---|
| c = −0.3 (1/W = 1.03084) | 1.02200 ± 0.04103 (−8.6e-3) | 1.03060 ± 0.01029 (−2.3e-4) | 1.03091 ± 0.00272 (+6.5e-5) |
| c = −1.5 (1/W = 2.02359) | 2.08454 ± 0.06273 (+3.0e-2) | 2.03541 ± 0.02560 (+5.8e-3) | 2.02558 ± 0.01141 (+9.9e-4) |
| c = −1 (1/κ = 1.57080) | 1.58214 ± 0.05754 (+7.2e-3) | 1.57354 ± 0.01394 (+1.8e-3) | 1.57201 ± 0.00299 (+7.7e-4) |

The relative errors fall with h at every member (by 2 to 37 per halving, not smoothly), so the slopes converge to 1/W at c = −0.3 and −1.5 and to 1/κ at c = −1. The wider window [1e-9, 1e-5] is
less good (−1.3e-4 at c = −0.3, −9.97e-2 at c = −1.5, h = 0.005): at c = −1.5 the series bends and saturates below |A − A\*| ≈ 3e-8 (T = 15.0 and 16.6 at h = 0.01, `dev_61`; Sec. 12).
With the first design's threshold (3e-3) the c = −1.5 floor was near 3e-6 (`dev_44`, `dev_45`) and the fit errors 6 to 8 %. A window starts above the floor: [1e-7, 1e-5] at h ≤ 0.01.

## 8. The control: c = −1

#212: no growing mode (the knife edge). Expected for a family through this member: no divergent lifetime at A = 0; small perturbations decay or stay static to T_max.

Found in development:
- **The master equation agrees with #212.** `evolve_master1d.py` evolves u_tt = u_xx − V̂u (V̂ in closed form, compact data, leapfrog, x ∈ [−90, 90], dx = 0.02, t to 60), with no double-null code, sharing only V̂.
  Growth rate of ln|u(x = 0)| over t = 45 to 60 (`dev_30`):

  | c | −2 | −1.5 | −1.1 | −1.01 | −1.001 | **−1** |
  |---|---|---|---|---|---|---|
  | measured | 0.3985 | 0.4942 | 0.6032 | 0.6331 | 0.6363 | **−0.037** (decaying) |
  | #212 W | 0.39845 | 0.49417 | 0.60319 | 0.63315 | 0.63627 | none |

  This confirms #212's table at five members to four digits, and its zero count at c = −1, in the linear problem.
- **The double-null prototype does not agree at c = −1.** Y grows at 0.632 to 0.637 from the truncation seed (table in Sec. 7), with the same shape as c = −1.01 to within 5 % along l = −8 … 6 and a factor
  12.5 to 13.1 in 4 time units along l = −6 … 6 (`dev_53`). The mode family has A\*ₕ ∝ h² there (ratio 3.91) and a lifetime slope 1.571 ± 0.009 at h = 0.01, which is 1/κ to 4e-4 (`dev_48`, `dev_49`).
  **The control as worded would fail in this code.**
- It is not a failure of #212's count (the independent run agrees with it) and the code reproduces W at four other members. Both facts together say the growth belongs to the double-null evolution in the tortoise
  gauge or to its characteristic data, and not to the master equation. Candidate: it is the same A-independent growth that floors c = −1.5; its rate there is not measured.

Proposed diagnosis, time-boxed to 2 h, before criteria are declared (done; the results are in Sec. 12):
(a) a **horizon-regular gauge**, x⁻ → Kruskal U = −e^{−κ√2 x⁻}/(κ√2) with x⁺ unchanged, in which the grid reaches and crosses the horizon, and the growth either stays or goes;
(b) the A-independent growth's rate at c = −1.5 by the untuned run, against W and κ;
(c) c = −1 at h = 0.005.
If the growth persists in the regular gauge it is in the evolution equations as posed and needs a separate look at #212's count with these data; if it goes, the tortoise gauge is the cause, c = −1.5's floor probably goes
with it, and the hybrid gauge becomes the gauge of the black members. The control is reported as found either way.

## 9. Grids and cost (timed on this Mac, one thread, `nice 10`)

Full run without a stop, c = −1.5, X = 40 (`dev_11`):

| h | N | time |
|---|---|---|
| 0.02 | 2000 | 0.42 s |
| 0.01 | 4000 | 1.40 s |
| 0.005 | 8000 | 5.31 s |
| 0.0025 | 16000 | 21.4 s |

A job (extraction, one bisection, 38 series points: about 75 runs, each ending at its departure) took 6.3, 21.6, 73.5 s at h = 0.02, 0.01, 0.005 for c = −0.3 (X = 40) and 6.0, 22.8, 83.9 s for c = −1.5 (X = 50).
Proposed grids for the runs of record: h = 0.01, 0.005, 0.0025 for c = −0.3, −1.5 and −1. The h = 0.0025 jobs are about 5 and 6 minutes. All runs of record together are well under 2 h on two processes; the 2 h stop applies per job.

## 10. Build plan

Files: `evolve.py` (kernel, data, extraction, departure, finders, bisection, series, fits, planted series, controls, checks), `evolve_output.txt`, `evolution.md`, RAW `I217-mac/README.md` with the criteria
written before the runs of record.

| id | check | target |
|---|---|---|
| E0 | kernel with V off or on at c = 0 equals `dn.evolve`, numba and numpy, three grids | max abs 0 (dev: 0 on one grid) |
| E1 | kernel right-hand side equals the sympy equations with V at 40 random states | 1e-14 |
| E2 | numba against numpy kernel, V on, three grids | 1e-12 |
| E3 | c = 0 rays against `dn.ellis_rays` | r 2e-9 (#39's I1), φ 1e-12, f 0 |
| E4 | static members c = 0, −0.3, −1.5, −1: second order for t ≤ 8, three h | 1.75 to 2.25 |
| E5 | Y growth from the truncation seed against #212's W, c = −0.3 and −1.5, two h | within 3 % (dev: 0.3 to 2.3 %) |
| E6 | 1D master equation against #212's W, c = −2 … −1.001; c = −1 decays | 1e-3 |
| E7 | extracted mode: Y(l) against #212's Y(ρ(l)) | declared after the first run |
| E8 | A\*ₕ at three h: order | 1.75 to 2.25 (dev: met at c = −0.3, c = −1; not at c = −1.5) |
| E9 | lifetime slope against 1/W, joint fit, declared window, two finest h, 0.5 % agreement | as #40 |
| E10 | planted series (slope 1/W recovered, 1/(√2 W) rejected, 40 seeds), time-unit mutation | as #40 |
| E11 | the control at c = −1 | reported as found |
| E12 | naming runs for each member and side (Sec. 6) | reported |

Mutations (targets computed in the development runs before declaration): `vterm_sign`, `vterm_coef`, `Vp_sign`, `eps_sign` (kernel), `time_unit`, `pstar_off` (fits), `mode_norm` (A scaled by 2: A\*ₕ changes by 2, the
slope does not), `gauge_part` (δ̂ with its f row zeroed: E7 and the departure fire).

## 11. What this does not settle

- The cause of the c = −1 growth in the double-null code, and of the c = −1.5 floor and its A\*ₕ. Sec. 8 has a plan.
- What state the side without a finder event ends in (Sec. 6): the naming runs report the mass and the signs, not an identification.
- The mode is extracted from the code; the analytic link to #212's gauge is checked in profile (E7), not built.
- Only the radial mode of two members and a control; nothing about polar modes (#218) or the sourced black bounce (#223).
- No sentence here or in the note to follow calls the threshold by a type name; only outcomes and slopes are reported.

## 12. Diagnosis of the c = −1 growth and of the c = −1.5 floor (rulings 4 and 5)

Prototype `evolve_dev_reg.py` (two step sizes), `evolve_dev_runs2.py`, logs `dev_60` to `dev_88`.

**(a) Horizon-regular gauge.** The kernel takes separate steps in x⁺ (tortoise, hp = 0.01) and in the second null coordinate. The tortoise run (z = y, hm = hp) and the Kruskal run (z = −e^{−κ√2 y}/(κ√2),
e^{−f} = A e^{κ√2 y}, regular on z = 0) use the same code and the same observer, Y on the diagonal, untuned member, slope of ln|Y| over t = 8 to 14 (`dev_70` to `dev_75`):

| member | tortoise | Kruskal x⁻ | #212 W | κ |
|---|---|---|---|---|
| c = −1.01 | 0.6300 | 0.6322 | 0.63315 | 0.6351 |
| c = −1.5 | 0.4905 | 0.4908 | 0.49417 | 0.5517 |
| **c = −1** | **0.6334** | **0.6357** | none | 0.6366 |

The Kruskal grid has hm = 4e-6 (2.8e5 rows) at c = −1.01 and −1 and 1.5e-5 (8.5e4 rows) at c = −1.5, so the dynamic range e^{κu} up to u = 14 is resolved to 3 % per cell; it ends at z = −1.5e-4 and −5.4e-4
and does **not cross** the horizon (the interior data are not built). **The growth at c = −1 is the same in the regular gauge.** By ruling 4 I stop and report here.

**What else is known** (all development, none explains it):
- **The growing solution satisfies the constraints.** C₊₊ and C₋₋ (the stencil of `dn.constraints`, minus the exact static member's) stay at their truncation level while Y grows 160-fold at c = −1:
  dC₊₊ = −1.9e-8 → +5.8e-9, dC₋₋ = 2.35e-7 → 2.27e-7 for Y = 3.8e-5 → 6.0e-3 (t = 6 → 14; `dev_76`). The same at c = −1.01, −1.5, −0.3. It is not a constraint-violating mode of the free evolution.
- **The linear response to a compact pulse is linear in its amplitude** (`dev_78`: A_p = 1e-5 and 1e-6 give the same Y per unit A_p to 1 to 2 % up to Y ≈ 0.08), so the odd part of Y in ±A_p is the linear response with the
  truncation seed removed. At c = −1.5 it grows at 0.495 over t = 10 to 16 (W = 0.494), as expected. At c = −1 and −1.01 it grows by three to four orders of magnitude between t = 12 and 20 (0.19 → 1.2e3 and −0.96 → +1.0e3 per unit A_p) and then **saturates at Y ≈ 7.7e3 A_p** (7.70e3 and 7.62e3 at t = 28),
  which no mode of the master equation does. At c = −1.5 it changes sign between t = 24 and 26. So the double-null linear response departs from the master equation's e^{Wt} at late times at every member, and at c = −1 it grows where
  the master equation decays (`dev_30`, to t = 60). I have no explanation. A pulse of 1e-3 is already nonlinear by t = 12 (`dev_77`).
- **The 3e-3 threshold is nonlinear at c = −1.5.** δ̂ is mostly gauge (δf ≈ 10 to 20 × Y), so |Y| = 3e-3 is a 5 % perturbation in δf. On the negative branch Y reaches −3e-3 and turns positive near t = 16 (`dev_60`:
  the tuned run grows at 0.489, 0.492, 0.467 over t = 4 to 13, stalls at |Y| = 1.5e-3 and then grows at 2.7). The sign of Y at 3e-3 is then not the sign of the mode. This is where the floor and the erratic A\*ₕ of c = −1.5 came from.
- **With the threshold at the linear level** (extraction at |Y| = 2e-4, departure at |Y| = 2e-5) c = −1.5 has A\*ₕ at second order and a slope at 1/W to 1e-3 (Secs 5 and 7).

**(b)** The "A-independent growth" of ruling 5 is the same mode: at the tuned amplitude of the first design its rate is W (0.489 to 0.492), not κ (`dev_60`). At c = −1.5 a saturation remains below |A − A\*| ≈ 3e-8: T = 15.0 and 16.6 at h = 0.01 (`dev_61`, `dev_86`) and 16.8 and 18.5 at h = 0.005 (`dev_83`), later by 1.8 per halving
(ln 4/W = 2.8 would be a seed ∝ h², so it scales as about h^1.3). I have not explained it. The windows start above it.

**(c)** c = −1 at three h: A\*ₕ in second order and the slope at 1/κ to 7.7e-4 at h = 0.005 (Sec. 7 tables).

**What follows for the criteria.** The control as worded fails in this code. It fails in the regular gauge too, with the constraints satisfied, and the linear response departs from the master equation's growth at late times
at every member. I propose: (1) take Y₀ = 2e-5 (extraction 2e-4) for the departure, which is the change from the ruling's 3e-3, windows [1e-7, 1e-5] at h ≤ 0.01; (2) declare the control's expected behavior as #212's
(no divergent lifetime at A = 0), report that the code shows a lifetime law with slope 1/κ, and name the open question: the double-null linear problem and the master equation disagree at c = −1 and at late times;
(3) the interior of the horizon (crossing it) and a direct linear double-null solver are the two next diagnostics, outside this issue. I wait for your ruling before I build.

## Superseded development

`dev_12` to `dev_29` used an extraction from a run with a pulse, at t ≥ 14. `dev_40` to `dev_49` (extraction at |Y| = 2e-3, departure at 3e-3) are superseded for c = −1.5 by Sec. 12. At c = −0.3 that was nonlinear (the corner node differed from the static member by 0.7). I found it by printing the raw corner value
and replaced it with the extraction of Sec. 5 (`dev_40` to `dev_49`). The logs stay; their c = −0.3 slope (1/W to 4e-4 at h = 0.005) and c = −1.5 slope were not used here.

## References

- [BF06] K. A. Bronnikov, J. C. Fabris, PRL 96, 251101 (2006), arXiv:gr-qc/0511109. [BKZ12] K. A. Bronnikov, R. A. Konoplya, A. Zhidenko, PRD 86, 024028 (2012), arXiv:1205.2224. Both as in `README.md` of this directory.

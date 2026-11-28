# The expanding symmetric throat's late state in closed form, its exponent k(n), and tests at D = 6 and X = 800 (#230)

Driver: [throat_scaling.py](throat_scaling.py), output of record: [throat_scaling_output.txt](throat_scaling_output.txt)
(`python3 throat_scaling.py analyze --i213=<negative-mass-raw>/I213-dh610 --i230=<negative-mass-raw>/I230-dh610`, exit 1 because four checks are not met). Raw output, with the
declarations committed before the runs: negative-mass-raw `I230-dh610` (README with the declaration and the results, `dev/`, `runs/`, `scripts/`, `figures/throat_scaling.svg`).
The evolution is #213's `expansion.stream` (#39's `dn.py`), unchanged. Status words as in AGENTS.md: every number is *numerical evidence* unless marked *proved*, *checked symbolically* or
*conditional*. No sentence here says a result is new.

Commits: the checks, bars, members and mutations were committed at 9b6b1f8 and written into the raw README (9843f9a0) before the first run of record. 977e873 adds report-only rows (the
"post hoc" section below) after the runs; the verdict lines of `analyze` at 9b6b1f8 and at 977e873 are identical (`runs/analysis_9b6b1f8.log`, `runs/analysis.log`).

## Answer

The late state of the mirror-symmetric expanding throat is an exact solution of the field equations with the sphere term dropped, with one free constant A (*checked symbolically*, below):

- r^n = x⁺ + x⁻ in suitable null coordinates, φ = A arcsin((x⁺ − x⁻)/(x⁺ + x⁻)), and e^f = C (x⁺x⁻)^{A²} ((x⁺ + x⁻)/2)^{(n−1)/n − A²}.
- On the throat it gives **H = k/(τ_s − τ), k = 2/(nA² − n − 1), |ρ|/H² = c = n²A²/2**, and r ∝ (τ_s − τ)^{−k}. The throat relation c = n(1/k + (n+1)/2) holds identically.

A is the amplitude of the field's step across the throat. **If A = (2/π) φ_∞**, with φ_∞ = (π/2) √(n/(n−1)) the end value of the static wormhole's field (hypothesis **H1**, not derived), then
**k = 2(n−1) and c = n³/(2(n−1))**. So the conjecture of #213 is now a *conditional* statement with its hypothesis named.

| D | n | k = 2(n−1) | c = n³/(2(n−1)) |
|---|---|---|---|
| 4 | 2 | 2 | 4 |
| 5 | 3 | 4 | 27/4 |
| 6 | 4 | 6 | 32/3 |
| 7 | 5 | 8 | 125/8 |
| 8 | 6 | 10 | 108/5 |
| 9 | 7 | 12 | 343/12 |
| 10 | 8 | 14 | 256/7 |

What the runs show (details and failures below):

- **D = 6 (the prediction): met.** k_fit = 5.99908, 5.99981, 5.99846 for A = 0.05, 0.2, 0.01 at h = 0.003 (bar 1e-2), and |ρ|/H² at r = 50a is 10.66756, 10.66688, 10.66814 against 32/3 = 10.66667 (bar 1e-3).
- **D = 4 at X = 800 (#213's task 2): met.** τ_s = 17.890219 (the X = 200 run gives 17.890221), k = 2.000000, K(τ_s − τ)² = 4.00000 at the last cell, τ_s − τ = 0.155.
- **The closed form follows the data, not only the Ellis family.** Adding a step b to the field's end value (the jump test) moves k to 2/(nA² − n − 1) with A = (2/π)(φ_∞ + b): k from 0.64 to 12.6 at
  D = 4, 5, 6, agreeing to 2e-4 (relative, k_inf). On the stored columns r^n = P(x⁺) + P(x⁻) holds to 4.5e-4 and φ = A arcsin(...) to 1.3e-3, with no fitted parameter.
- **Four of the 25 declared checks are not met:** K5F (#213's fitted k at D = 5, 3.1e-3 against the bar 1e-3, expected), JTC (1.06e-3 against 1e-3 at one jump-test member, expected), and K6X and WS (the extrapolated exponent
  of one fine-grid run, D = 6, A = 0.2, is noise-limited). Each is explained below and none is changed.

## Derivation (*checked symbolically*; by hand)

Conventions of `dn.py`: ds² = −2e^{−f} dx⁺dx⁻ + r² dΩₙ², n = D − 2, R_ab = −∂_aφ ∂_bφ (ghost), 8πG = 1. The field equations are

- (E1) r₊₋ = −(n−1)(r₊r₋ + e^{−f}/2)/r,
- (E2) f₊₋ = (n/r) r₊₋ − φ₊φ₋,
- (E3) φ₊₋ = −n(r₊φ₋ + r₋φ₊)/(2r),
- (C±) r_{±±} + f_± r_± = (r/n) φ_±².

The mirror-symmetric data have φ(x⁺, x⁻) = −φ(x⁻, x⁺), and r, f symmetric. On the throat x⁺ = x⁻: ρ = −e^f φ₊², H = e^{f/2}(r₊ + r₋)/(√2 r), dτ = √2 e^{−f/2} dx.

1. **The sphere term.** Drop e^{−f}/2 in (E1). Then r r₊₋ + (n−1) r₊r₋ = 0, which is (rⁿ)₊₋ = 0, so rⁿ = P(x⁺) + P(x⁻) (the mirror symmetry makes the two functions equal). Reparametrize
   x^± → P(x^±), which keeps the form of the equations and the symmetry. Then **rⁿ = x⁺ + x⁻ = 2T**, with X = (x⁺ − x⁻)/2.
2. **The field.** (E3) becomes φ_TT + φ_T/T = φ_XX, the cylindrical wave equation (T the radius, X the time). Its odd solution of degree 0 is **φ = A arcsin(X/T)**: for φ = Φ(X/T) the equation is
   (ξ² − 1)Φ'' + ξΦ' = 0, so Φ' ∝ (1 − ξ²)^{−1/2}, and the other solution is a constant.
3. **The conformal factor.** (C±) give f_± = (n−1)/(n S) + S φ_±² with S = x⁺ + x⁻, a closed one-form by (E3). With the φ above,
   **e^f = C (x⁺x⁻)^{A²} ((x⁺ + x⁻)/2)^{(n−1)/n − A²}**. (E1) to (C±) are then satisfied exactly with the sphere term dropped (`throat_scaling.py derive`, check CF1: residuals identically zero for
   n = 2 to 6 with A free).
4. **The throat.** On x⁺ = x⁻ = x: e^f = C x^p with p = A² + (n−1)/n, H = e^{f/2}/(√2 n x), dτ = √2 e^{−f/2} dx. So 1/H = √2 n x e^{−f/2}, and τ_s − τ = ∫ₓ^∞ √2 e^{−f/2} dx' is finite iff p > 2, that is
   nA² > n + 1. Then 1/H = (nA² − n − 1)(τ_s − τ)/2, which is the line in the Answer (CF2). Also ρ/H² = −n²A²/2 and r ∝ x^{1/n} ∝ (τ_s − τ)^{−k}.
5. **Size of the dropped term.** On the diagonal e^{−f}/(2 r₊r₋) ∝ x^{−(A² − (n−1)/n)}, which tends to zero when nA² > n + 1. The data confirm it: the sphere term of the throat relation, (n−1)/(2r²H²), at the last cell is below 3.4e-7 (D = 4), 2.3e-5 (D = 5) and 1.0e-5 (D = 6) for the Ellis-family runs,
   and below 2.9e-4 for the jump-test runs (the largest, 1.3e-4 and 2.9e-4, are D = 5 and 6 at b = −0.03, where k is 5.5 and 12.6).

(The equation in step 2 is a textbook one; no claim is made about where this solution has or has not been written down.)

## H1: what fixes A, and what would be needed to derive it

The solution tends to φ = ±(π/2)A on the two null edges x⁻ → 0 and x⁺ → 0 of the quadrant, so A is twice the edge value over π. **H1:** the field on the null boundary of the anti-trapped region, as x⁺ → ∞, takes its value at
spatial infinity, φ_∞. For the Ellis family φ_∞ is fixed by the static solution: r'² = 1 − (a/r)^{2(n−1)}, φ' = Q/rⁿ with Q = √(n(n−1)) a^{n−1}, so φ_∞ = Q ∫ dl/rⁿ = (π/2) √(n/(n−1)) (CF3: `dn.ellis_dense` agrees
to 2e-11 for n = 2 to 6). Then A² = n/(n−1).

Evidence for H1, all *numerical evidence*:

- c at r = 50a is within 3.7e-4 (D = 4, 5) and 1.4e-4 (D = 6) of n³/(2(n−1)) for the nine Ellis-family runs of the gated sets, and the extrapolated exponents agree with 2(n−1) to 1e-5 at D = 4 and 5, and the throat-relation series to 2e-6 at D = 6.
- The spatial test: φ = √(n/(n−1)) arcsin((P⁺ − P⁻)/(P⁺ + P⁻)) on 22 runs to 1.3e-3 (absolute), with P from the run's own diagonal and no free parameter (SP2).
- The jump test is the sharpest: data whose field end value is shifted by b follow A = (2/π)(φ_∞ + b). A formula for k that is only a coincidence of the Ellis numbers would not move with b.

What a derivation would need (open; none of it is done here):

1. **A decay estimate in the exterior.** φ(x⁺ → ∞, x⁻) = φ_∞ for every x⁻ up to the contour r₋ = 0. In the exterior the field is φ_∞ − Q/((n−1) r^{n−1}) plus a radiative part; one needs the charge Q = rⁿφ_l to be
   conserved at the end and the radiative part to vanish at null infinity, for the nonlinear solution up to the contour.
2. **A matching argument between the static throat and the late state.** At the contour the sphere term e^{−f}/2 is of the size of r₊r₋, so the quadrant solution does not apply there. The layer has to be solved (or bounded) to show it
   does not change the edge value. The development data say it changes A by less than 1e-5 at D = 4 (the late values of |ρ|/H² are 4.00000).
3. **An attractor statement for the cylindrical wave equation.** The Goursat problem in the quadrant with edge data tending to ±φ_∞ has φ → (2φ_∞/π) arcsin(X/T) at fixed X/T. The central slope is a causal functional of the edge data (an
   Abel-type transform of the edge values); its correction terms need a bound. For an edge tail δ ∝ x̃^{−γ} the odd similarity response diverges at ξ → 1 like (1 − ξ)^{1/2−γ} when γ ≥ 1/2, so the approach to the limit is not a
   similarity solution there, and the observed exponent of the approach (g = 1.7 to 2.9 in r, below) is not derived.

## The runs

Family S of #40 at the same pulse (φ-Gaussian, centre 3, width 0.5, ghost), `expansion.stream` unchanged. k_fit is #213's fit (1/H linear in τ over the last e-fold of r), k_inf is the local exponent k_loc = H²/(dH/dτ)
extrapolated as k_inf − a (r/r_end)^{−g} over the same window, c is |ρ|/H² at r = 50a. Declared in the raw README at 9843f9a0.

### Checks (the output file has every line)

| id | what | verdict | number |
|---|---|---|---|
| CF1, CF2 | the closed form solves the flat-limit equations and constraints; k, c and the throat relation on its diagonal (sympy, n = 2 to 6) | pass | residuals 0 |
| CF3, CF4 | Ellis end value and charge; k = 2(n−1), c = n³/(2(n−1)) from A² = n/(n−1) | pass | 2.3e-11, 0 |
| PF, PX, PS | planted fit, extrapolation, spatial test | pass | 5e-5 and 1e-5; 2e-5; exact 4e-9, +2 % in A gives 2.3e-2 |
| K4F, K4X, C4 | D = 4, four stored members: k_fit, k_inf, c against 2 and 4 | pass | 2.0e-4, 4e-6, 3.7e-4 |
| **K5F** | D = 5, #213's fitted k against 4 (bar 1e-3) | **not met** | 3.13e-3 (3.98749, 3.98973) |
| K5X, C5 | D = 5: k_inf against 4, c against 27/4 | pass | 8e-6 (4.00003), 3.7e-4 |
| K6F, C6, RES6 | D = 6, h = 0.003, three members: k_fit within 1e-2 of 6, c within 1e-3 of 32/3; h = 0.006 against 0.003 | pass | 1.5e-3, 1.4e-4; H 3.9e-3, k 6e-6 |
| **K6X** | D = 6: k_inf within 1e-2 of 6 | **not met** | 2.8e-2 (A = 0.2: 6.02838; A = 0.05: 6.00163; A = 0.01: 6.00033) |
| **WS** | k_inf over the last 1.0 and 1.5 e-folds within 1e-3 (relative) | **not met** | 4.7e-3 (A = 0.2 at D = 6; every other run below 8e-4) |
| T800, K800, KT800 | D = 4 tuned, X = 800: τ_s within 1e-3 of 17.8902, k = 2.000 ± 1e-3, K (τ_s − τ)² = 4.000 | pass | 17.890219; 2.000000; 4.00000 |
| JTK | jump test, 11 runs: k_inf against 2/(nA² − n − 1) | pass | 2.2e-4 |
| **JTC** | jump test: c at r = 50a against n²A²/2 (bar 1e-3) | **not met** | 1.06e-3 (D = 4, b = 0.4; the next worst is 6.4e-4) |
| SP1, SP2 | r^n = P⁺ + P⁻ and φ = A arcsin(...) on the stored columns, 22 runs | pass | 4.5e-4 relative; 1.3e-3 absolute |

21 of 25 pass. The mutations: 18 declared, each fires exactly its declared set (`runs/mut/rec_*.log`).

### D = 4 and 5: the stored runs

The coordinator's ruling is that K5F is kept as declared: it compares #213's fitted k, whose window ends before τ_s, so it is a pre-asymptotic estimate, and it is **not met** (3.1e-3 and 2.6e-3 below 4 against 1e-3). The local exponent
and c show where it comes from. At r = 50a: k_loc = 3.988 and c = 6.7525 (3.7e-4 above 27/4). The throat relation gives δk/k = −n² δc/c for c = n³/(2(n−1)), so a 3.7e-4 offset in c is a 3.3e-3 offset in k at D = 5 and 1.6e-2 at D = 6.
The extrapolated exponent over the same window, 4.00003 (D = 5, both members) and 2.00000 (D = 4), is a **second estimator applied to the same runs**, not an independent measurement; it is reported beside K5F as the evidence that the closed form's
k = 4 holds, together with C5.

### D = 6: the prediction

| member | h | r_end | τ_s − τ_end | k_fit | k_inf | k from the throat relation, extrapolated (post hoc) | c at r = 50a | K (τ_s − τ)² (k² = 36) |
|---|---|---|---|---|---|---|---|---|
| A = 0.05 | 0.003 | 168.6 | 1.93 | 5.99908 | 6.00163 | 5.99999 | 10.66756 | 35.9889 |
| A = 0.2 | 0.003 | 186.5 | 0.59 | 5.99981 | 6.02838 | 6.00000 | 10.66688 | 35.9977 |
| A = 0.01 | 0.003 | 166.7 | 2.56 | 5.99846 | 6.00033 | 6.00000 | 10.66814 | 35.9814 |
| A = 0.05 | 0.006 | 168.6 | 1.93 | 5.99909 | 5.99996 | 5.99999 | 10.66755 | 35.9890 |

The two resolutions agree on c to 1e-6 and on k_fit to 1e-5; τ_s to 6e-3 (8.17309 and 8.16739), which is the fit's sensitivity at τ_s − τ = 1.9. The approach to k = 6 is a power of r: c/c_pred − 1 = 4.1e-3, 7.5e-4, 8.4e-5, 1.6e-5 at
r = 10, 20, 50, 100 (A = 0.05), g = 2.4 from the throat-relation series.

**Why K6X and WS fail.** k_loc = H²/(dH/dτ) differentiates H, which carries a cell-to-cell noise of about 3e-9 (relative, rms of the second difference) at h = 0.003. H changes by only 3e-6 to 9e-6 per cell, so k_loc scatters by 1.8e-3
(3e-4 at h = 0.006; in D = 4 at h = 0.01, 1e-5 at X = 200 and 3.5e-4 at X = 800). The mean of k_loc over the window equals the mean of k from the throat relation (no differentiation) to 3e-7, so there is no bias. For A = 0.2 the approach is below the scatter over most of the window
(about 6e-5 at r_end: c is within 3e-6 of 32/3 at r = 100a), the three-parameter extrapolation is ill-conditioned, and g runs to its lower bound 0.05 (k_inf = 6.028). A planted series with this scatter and a 7e-5 approach at r_end gives errors of rms 7e-3 and up to 1.7e-2, three of
20 seeds beyond 1e-2 (post hoc row of the output); the declared PX used scatter 1e-4, which is why it could not show this. So the failures are the estimator's noise limit at fine h, not a deviation of the run. The throat-relation series
extrapolates to 5.99999 to 6.00000 with a fit rms of 2e-7 (post hoc).

### D = 4 at X = 800 (#213's task 2)

N = 80000, h = 0.01, t = 1131, 88 min at one thread; the last cell is at τ_s − τ = 0.155 (0.315 at X = 200), r/a = 790. The fit gives τ_s = 17.890219 for windows of 0.5, 1 and 2 e-folds (the X = 200 run: 17.890220 to 17.890228), k = 2.000000,
K (τ_s − τ)² = 4.00000 and r (τ_s − τ)² = 18.9049 at τ_s − τ = 0.5, 0.3 and 0.155 (r ∝ (τ_s − τ)^{−2}, to the printed digits). The differenced k_loc at the last cell reads 1.99870 (scatter 3.5e-4 at this length, the noise discussed under D = 6); no bar uses it.

### The jump test

The ordinary S pulse (A = 0.05) plus a `dphi` pulse (centre 3, width 0.5) that raises the field's end value by b on each ray, with the constraint solved for the total. A = (2/π)(φ_∞ + b):

| D | b | A | k predicted | k_inf | c predicted | c at r = 50a | c extrapolated (post hoc) |
|---|---|---|---|---|---|---|---|
| 4 | 0.05 | 1.44604 | 1.6919 | 1.6919 | 4.1821 | 4.1817 | 4.1821 |
| 4 | 0.1 | 1.47788 | 1.4617 | 1.4617 | 4.3682 | 4.3675 | 4.3682 |
| 4 | 0.2 | 1.54154 | 1.1411 | 1.1411 | 4.7527 | 4.7510 | 4.7527 |
| 4 | 0.4 | 1.66886 | 0.7782 | 0.7782 | 5.5702 | **5.5643** | 5.5702 |
| 4 | −0.05 | 1.38238 | 2.4332 | 2.4332 | 3.8220 | 3.8222 | 3.8220 |
| 5 | 0.2 | 1.35207 | 1.3475 | 1.3475 | 8.2264 | 8.2288 | 8.2264 |
| 5 | 0.5 | 1.54305 | 0.6363 | 0.6364 | 10.7146 | 10.7179 | 10.7143 |
| 5 | −0.03 | 1.20565 | 5.5440 | 5.5452 | 6.5411 | 6.5453 | 6.5411 |
| 6 | 0.1 | 1.21836 | 2.1330 | 2.1330 | 11.8753 | 11.8754 | 11.8753 |
| 6 | 0.3 | 1.34569 | 0.8915 | 0.8915 | 14.4870 | 14.4870 | 14.4870 |
| 6 | −0.03 | 1.13560 | 12.6289 | 12.6290 | 10.3167 | 10.3204 | 10.3168 |

JTC as declared reads c at r = 50a and misses its bar once (the bold entry, 1.06e-3 against 1e-3): for A² = 2.78 the approach is slower, and at r = 117a the same run reads 5.5691. The extrapolated c column is the coordinator-allowed supplementary row,
labelled post hoc; it agrees to 2.4e-5 or better and gates nothing. Negative b of −0.1 and −0.2 at D = 4 collapses (trapped cells) in the development run and is not part of the test. For k < 1 (b = 0.4, 0.5, 0.3 above) r ∝ (τ_s − τ)^{−k} still has
1/r² = o(H²).

### The spatial test

On the stored columns x⁺ = X/4, X/2, X at x⁻ ≥ 20 (22 runs), r^n/(P⁺ + P⁻) − 1 is below 4.5e-4 and |φ − A arcsin(ξ)| below 1.3e-3. D = 4 is the sharpest (9.8e-6 and 1.5e-6 at x⁺ = 800 in the X = 800 run). The r^n deviation halves per doubling of x⁺ (D = 4, tuned: 1.6e-4, 8.6e-5, 4.1e-5 at x⁺ = 50, 100, 200).

### The approach to the limit (reported, not gated)

k_loc(r) and c(r) tend to the closed-form values as a power of r: fitted exponents g = 1.7 to 2.9 (D = 4, 5, 6 and the jump test). c/c_pred − 1 at r = 10, 20, 50a is
2.0e-3, 1.4e-4, 8e-7 (D = 4), 1.1e-2, 2.6e-3, 3.7e-4 (D = 5) and 4.1e-3, 7.5e-4, 8e-5 (D = 6, A = 0.05). The approach is slowest at D = 5. The exponent is not derived (item 3 of the H1 list above).

## Discrepancies and things that failed

1. **K5F not met** (3.13e-3 against 1e-3), as declared and expected from #213's table; the cause is the pre-asymptotic window (above). Not changed.
2. **JTC not met** at D = 4, b = 0.4 (1.06e-3 against 1e-3), expected from the development run. The coordinator ruled to keep it as declared. Not changed.
3. **K6X and WS not met**, both from the D = 6, A = 0.2 run at h = 0.003: the declared extrapolation is noise-limited. The declared PX control used noise 1e-4 and could not have shown it. Not changed; post hoc rows added.
4. **The declaration post on the issue said 29 checks.** The script has 25 checks plus BUILD in the mutation runs. The raw README lists them one by one.
5. **Appendix B of Paper V** says 1/r² = o(H²) for k > 1. With H = k/(τ_s − τ) and r ∝ (τ_s − τ)^{−k}, 1/(r²H²) ∝ (τ_s − τ)^{2k+2}, which tends to zero for every k > −1. The statement is true as written (sufficient), but the
   jump test has k < 1 members (0.64 to 0.89), for which it is silent.
6. **Figure.** One is in the raw folder (`figures/throat_scaling.svg`, `scripts/make_figure.py`): k_loc/k against r for D = 4, 5, 6, and k against A with the closed form's curves and the jump-test points. The graph skill's
   "better" arrow is omitted because neither axis has a better direction.
7. **Memory.** Sampled by `ps` at five minutes of run time: 0.21 to 0.26 GB per run (X = 800 included). Not recorded later in the runs.

## Proposed text (for the coordinator)

**`expansion.md`**, a note under "Remark" in the review and in "Recommended tasks": tasks 1 and 2 are done in #230. "k = 2(n−1) is *conditional* on H1 (`throat_scaling.md`): the closed form r^n = x⁺ + x⁻, φ = A arcsin((x⁺ − x⁻)/(x⁺ + x⁻)) is *checked
symbolically*, gives H = k/(τ_s − τ) with k = 2/(nA² − n − 1) and c = n²A²/2, and A = (2/π)φ_∞ of the static wormhole gives k = 2(n−1), c = n³/(2(n−1)). At D = 4, X = 800: τ_s = 17.890219, k = 2.000000, K(τ_s − τ)² = 4.00000 at τ_s − τ = 0.155."
Also: "the D = 5 fitted k of 3.9875 is a pre-asymptotic value: the extrapolated exponent is 4.00003 (same runs)."

**`critical_ghost.md` §1.3**, after the #225 update:

- **Update from #230** (`wormhole_evolution/throat_scaling.md`; *numerical evidence* except where marked): the late state of the symmetric expanding throat is a closed-form solution of the field equations with the sphere term dropped (*checked symbolically*): rⁿ = x⁺ + x⁻,
  φ = A arcsin((x⁺ − x⁻)/(x⁺ + x⁻)), H = k/(τ_s − τ) with k = 2/(nA² − n − 1) and |ρ|/H² = n²A²/2. The constant A is the step of the field across the throat; with A = (2/π)φ_∞ of the static wormhole (hypothesis H1, not derived) k = 2(n−1): 2, 4, 6 at D = 4, 5, 6.
- D = 6 (new run, h = 0.003, three amplitudes): k_fit = 5.998 to 6.000, |ρ|/H² = 10.667 (32/3). D = 4 to τ_s − τ = 0.155: τ_s = 17.890219, k = 2.000000.
- Adding a step b to the field's end value changes k as the formula says (k from 0.64 to 12.6 at D = 4, 5, 6 to 2e-4). So the exponent is a property of the field's total step, not a universal constant.
- Not found: a derivation of H1 (needs a decay estimate in the exterior and a matching argument across the layer where the sphere term is of the size of r₊r₋).

**Paper V**, `main.tex` (read only here):

- Replace the `\pending{\#230}{...}` after the Reading paragraph with: "\emph{The late state.} Near $\tau_s$ the sphere term of the radial equation is negligible. In null coordinates in which $r^n=x^++x^-$ the field obeys the cylindrical wave equation, and
  $\phi=A\arcsin[(x^+-x^-)/(x^++x^-)]$ with $e^f\propto(x^+x^-)^{A^2}[(x^++x^-)/2]^{(n-1)/n-A^2}$ solves all five field equations of Appendix~\ref{app:equations} (\emph{checked symbolically}). On the throat it gives Eq.~\eqref{eq:ft} with
  $k=2/(nA^2-n-1)$ and $c=n^2A^2/2$. Taking $A=(2/\pi)\phi_\infty=\sqrt{n/(n-1)}$ from the static wormhole's end value (hypothesis H1, not derived) gives $k=2(n-1)$ and $c=n^3/[2(n-1)]$: $2$, $4$, $6$ and $4$, $27/4$, $32/3$ in $D=4$, $5$, $6$.
  A run at $D=6$, $h=0.003$ gives $k_{\rm fit}=5.998$ to $6.000$ and $|\rho|/H^2=10.667$ at $r=50a$; $D=4$ to $\tau_s-\tau=0.155a$ gives $\tau_s=17.890219a$ and $k=2.000000$. A step $b$ added to the field's end value moves $k$ to $2/(nA^2-n-1)$ with $A=(2/\pi)(\phi_\infty+b)$, for
  $k$ from $0.64$ to $12.6$, to $2\times10^{-4}$."
- Table~\ref{tab:expand}: a $D=6$ row ($A=0.05$, $h=0.003$: $k_{\rm fit}=5.9991$, $|\rho|/H^2=10.668$ at $r=50a$). The tabulated $D=5$ fits ($3.9875$, $3.9897$) are pre-asymptotic: the extrapolated exponent over the same window is $4.0000$; say so beside the table.
- Reading paragraph and the Limits bullet: "$k=2(n-1)$ is a \emph{conjecture} (\#230 pending)" becomes "$k=2(n-1)$ is \emph{conditional} on H1 (the field's end value fixes $A$)".
- Appendix B: "so $1/r^2=o(H^2)$ for $k>1$" can read "for $k>-1$".

## Limits

- One observer, the mirror-symmetric diagonal. One-sided families (#225, #42) have two different end values; the closed form is odd and does not apply to them as written.
- The sphere term is dropped. Its size is shown small on the diagonal (below 2.3e-5 at the ends of the Ellis-family runs, 2.9e-4 for the jump-test runs) and the spatial test holds for x⁻ ≥ 20 only; near the contours r_± = 0 it is of the size of r₊r₋ and the closed form does not apply.
- H1 is a hypothesis, supported by the runs and the jump test (3 values of D, 11 steps) but not derived.
- τ_s is an extrapolation and the throat is not followed to it: the last cell is at τ_s − τ = 0.155 (D = 4), 0.59 to 2.6 (D = 6). Whether the throat reaches τ_s is not measured.
- The extrapolated exponent assumes a power-law approach; the planted control covers noise 1e-4 only, and at the fine-h scatter (1.8e-3) the declared estimator fails (K6X, WS). The throat-relation series is the better-conditioned estimator.
- The runs share the evolution code. What the closed form shares with it is the equations of `dn.py` (CF1 uses them as written); the checks that read the runs are independent of the sympy derivation, and the spatial test reads the run's columns, not the diagonal.
- The jump test steps are small perturbations of the data with the constraint solved for the total; D = 4 steps below −0.05 collapse.

## Recommended tasks

1. **Derive H1.** What: the three items under "H1" above. Why: it turns k = 2(n−1) from *conditional* into *proved*. Acceptance (can fail): a statement with its hypotheses, checked on the stored columns by the spatial test with A from the derivation, and the jump test's
   k and c reproduced from it. Depends on: none.
2. **The approach exponent.** What: derive or measure the first correction to the closed form (g = 1.7 to 2.9 in r) from the linearization about it. Why: it sets how far a run must go for k_fit to reach the bar (D = 5's K5F). Acceptance: a predicted g and amplitude for the D = 5 and 6 runs
   within 20 % of the extrapolation's. Depends on: none.
3. **Asymmetric ends.** What: one-sided data (#225, #42) with end values φ_R ≠ −φ_L; whether k follows the total step (φ_R − φ_L)/π through the same formula. Why: tests H1 off the mirror symmetry. Acceptance: k of the locus within 1e-3 of the prediction or the discrepancy named. Depends on: #225.
4. **A better-conditioned exponent estimator for fine grids.** What: declare the throat-relation series (k = H²/Ḣ from ρ, no numerical derivative) as the extrapolation input, with the planted control at the observed scatter. Why: K6X and WS fail only from the differenced series' noise. Acceptance: planted case with scatter 2e-3 recovered
   within 1e-3. Depends on: none.
5. **Copy the figure.** `figures/throat_scaling.svg` (raw folder) can be copied next to this note. Not written into the repository because the brief restricts the files.

## Files

`throat_scaling.py`, `throat_scaling_output.txt`, this note; raw `negative-mass-raw/I230-dh610` (README, `dev/`, `runs/`, `scripts/`, `figures/throat_scaling.svg`).

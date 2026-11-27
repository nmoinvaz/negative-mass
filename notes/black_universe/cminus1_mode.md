# The growing mode of the c = −1 black universe in closed form, and the slope of its rate (#262)

Code [cminus1_mode.py](cminus1_mode.py), output of record [cminus1_mode_output.txt](cminus1_mode_output.txt), assembled from negative-mass-raw `I262-mac/runs`; the criteria,
written before the runs of record, are in `I262-mac/README.md`. Builds on #212 ([README.md](README.md)), #232 ([linear.md](linear.md)) and #248 ([horizon_sector.md](horizon_sector.md)).
Status words as in AGENTS.md: the exact solution and the monodromy entries are *checked symbolically*, the numbers are *numerical evidence*. No sentence here says a result is new
(#260 reads the literature).

## Answer

**1. The growing mode has a closed form.** With x = ∫dρ/A the tortoise coordinate and R = √(1 + ρ²),

```
Y_ref = R e^{κx},        Ψ = Y_ref e^{−κx} = R,        Y = c₂ V R(ρ(Π))        (c = −1, W = κ = 2/π)
```

solves #212's equation, *checked symbolically* (S1). The residual of Y = R e^{Wx} is `2(πW + 2c) e^{Wx}/(ρR(…))`, so the same formula solves the equation for **every** member c at the
rate W = K = −2c/π = 3m (b = 1); at c = −1, K = κ. At c = 0 (K = 0, W = 0) it is `Y = R`, the variation of the throat radius b in the gauge δφ = 0 (S2). The solution is in the
Kruskal-regular class (Ψ = R is analytic in ρ and in Π), equals `c₂V` on the future horizon, and has the series `ρ^{1/2}(1 + (π/8)ρ + (5/12 + 3π²/128)ρ² + …)`. Against that
series it agrees to 5.9e-61 (bar 1e-30), and against the numerical integration of the master equation on ρ ∈ [−40, 40], static region and interior, to 4.3e-42 (bar 1e-20);
a coefficient changed by 1% (R² → 1 + 1.01ρ²) gives 5e-5 and 5e-3 against the bars 1e-30 and 1e-20.

**2. The slope.** The equation has the exact solution at W = K for every c, and the rate through c = −1 is W(c) = K(c) + (c + 1)/I₀ + O((c+1)²):

```
dW/dc|₋₁ = 1/I₀ − 2/π = 0.34869547355073775,        d(W/κ)/dc|₋₁ = −(π/2)β* = π/(2I₀) − π²/8 = 0.31402901888733605,
I₀ = ∫₀^∞ 2Kρ² e^{−2κx}/(A R²) dρ = 1.0149036099284089 ,        β* = π/4 − 1/I₀ = −0.19991708252087078
```

The derivation is first-order perturbation theory of the monodromy about the c = −1 mode space (below). The slope agrees with the fit of #212's rates to 2.9e-5 (second-order fit, s₁ = 0.31403797) and 2.7e-8 (third-order fit, 0.31402903) (bar 1e-3), and
dW/dc with the central difference to 2.7e-8 (central difference 0.3486954830). PSLQ at 50 digits finds no relation for 1/I₀ or for the slope in either declared basis (B1 = {1, π, π², π³, 1/π, 1/π²}, B2 = {1, π, π², 1/π, ln 2, ln π, G, γ}; the planted number is recovered and the control returns none), so neither number is identified: I₀ stays a quadrature of the closed form.

**3. What this means.** The "knife edge" is one smooth branch. At c = −1 the solution that is regular on the horizon and outgoing at infinity,
`Y_∞ = Y_ref − Y_b/I₀ = Y_ref ∫_ρ^∞ 2Kr²e^{−2κx}/(A R²) dr / I₀`, exists and has W = κ; for c ≠ −1 the one growing mode of #212 is its continuation, with W(c)/κ(c) smooth through
c = −1 (slope 0.3140; #212's table has 0.31). #212's "none" at c = −1 is the count for the single branch u = Y/ρ → 0, which is Y_b (Y ∼ ρ^{3/2}); at W = κ the whole
two-dimensional solution space is Kruskal-regular, and the condition at infinity picks Y_∞ in it. Y_ref itself grows as e^{κx} toward large ρ.

```
 c < −1 (ρ_h > 0)       c = −1 (ρ_h = 0)       c > −1 (ρ_h < 0, throat inside the static region)
 W < K                  W = κ = K              W > K          (#212's table: 8 members below, 11 above; `table`)
 W(c) − K(c) = (c + 1)/I₀ + O((c + 1)²),   1/I₀ = 0.98532
```

## 1. The equation and the closed form

#212 (README Sec. 2): Y = δR in the gauge δφ = 0, e^{Wt},

```
Y'' = P Y' + Q Y,   P = 2/ρ − A'/A,   Q = 1/R⁴ − A'/(AρR²) − 2/(AR²) + W²/A²,   A = 1 + cR² + K(ρ + R² arctan ρ),   K = −2c/π,   dx = dρ/A.
```

N0: the equation as coded here is #212's `coef_x` to 3.2e-15. With Y = R e^{Wx} (E = e^{Wx}, E′ = WE/A) the terms in W²/A² cancel and what is left is

```
(Y″ − P Y′ − Q Y)/e^{Wx} = 2 (πW + 2c) / ( ρ R (2cρ² arctan ρ − πcρ² + 2cρ + 2c arctan ρ − πc − π) ).
```

At W = K it is zero for all ρ and c (S1, sympy). Interpretation: in the gauge δφ = 0, Y_K = R e^{Kv} is the perturbation δr = ε r e^{3m v}; at c = 0 it is the static variation of b (S2). Why it exists for
every c is not understood here (recommended task 3).

**How it was found.** Not by a search over the class below: the formal series of Y at ρ = +∞, run as a convergence test (N3), had the coefficients 1, 0, 1/2, 0, −1/8, … of √(1 + ρ²) for the e^{+κx} solution,
and a ratio growing like n for the e^{−κx} one. The raw README records the order.

**Normalization and the other solutions.** e^{κx} = ρ^{1/2} exp(H(ρ)), H = ∫₀^ρ (K/A − 1/(2ρ′)) dρ′ (analytic at 0), so Y_ref = ρ^{1/2}(1 + g₁ρ + …) with g₁ = 1/(4K) = π/8. The
Frobenius recursion at W = κ has g₁ free (the exponents 1/2 and 3/2 differ by 1 and the logarithm vanishes, #248 D7): g₂ = (40K² + 24Kg₁ + 3)/(96K²). In terms of g₁:

| g₁ | solution | g₂ |
|---|---|---|
| 0 | #248 D7, D8's start | 5/12 + π²/128 |
| −1/(4K) = −π/8 | Y_nat, the continuation in W of the series, g₁(w) = −w/(4K) | 5/12 − π²/128 |
| 1/(4K) = π/8 | Y_ref = R e^{κx} | 5/12 + 3π²/128 |

The second solution, with Y ∼ ρ^{3/2}, is Y_b = Y_ref ∫₀^ρ 2Kr² e^{−2κx}/(A R²) dr (Abel: the Wronskian is 2Kρ²/A, and Y_ref² = R² e^{2κx}). Y_ref = Y_nat + (π/4)Y_b.

**Classification (declared class, README of the raw folder).** 𝒦 = {e^{kκx} A^s F, F ∈ ℂ(ρ, R, arctan ρ)}, which holds the issue's example. Every member has a convergent expansion
e^{kκx}ρ^pΣc_nρ^{−n} at +∞; the equation has two formal solutions there, e^{±κx}ρΣc_n^±ρ^{−n}, unique up to scale (N3). The e^{+κx} series is that of R (120 terms to 1.9e-122 at
300 digits), so Y_ref ∈ 𝒦. The e^{−κx} series diverges: |c_{n+1}/c_n|/n = 0.797 at n = 100 against 1/(2K) = π/4 = 0.7854 (*numerical evidence*), so no member of 𝒦 is a solution
of that type. Therefore 𝒦 ∩ solutions = ℂ·Y_ref. Y_b and Y_∞ are quadratures of Y_ref, not closed forms of this kind.

## 2. The slope by first-order monodromy

Let ε = c + 1. The Y equation has regular singular points at ρ_h(ε) ≈ −(π/4)ε (A = 0; exponents ±W/A′(ρ_h) = ±w′/2, w′ = W/κ_c) and at ρ = 0 (exponents 0, 3, both solutions
analytic, #212), which merge at c = −1 into the point with exponents 1 ∓ w/2. Take the loop |ρ| = r < 1, counterclockwise, away from both. Its monodromy M(ε, W) in a fixed basis
is analytic in (ε, W); ρ = 0 has none, so M is the monodromy about ρ_h alone: eigenvalues e^{±iπw′}, det M = 1. The ingoing solution (Y = e^{Wx}·analytic at ρ_h) is the eigenvector of
e^{+iπw′}.

At (ε, W) = (0, κ): M = −1 (no logarithm), basis Y_a = Y_nat, Y_b. Write M = −(1 + n), n = εJ + δN₂ + O(2), δ = W − κ, with

```
J = ∮ Φ₀⁻¹ B₁ Φ₀ dρ  (B₁ = ∂_ε of the first-order system at fixed W, Φ₀ = [[Y_a, Y_b], [Y_a′, Y_b′]]),        N₂ = (iπ/κ) diag(−1, +1)   (from e^{∓iπw}).
```

J is a residue at ρ = 0 of products of Frobenius series (S3, exact in π, series to order 8):

```
J = [[ iπ(π² − 8)/8,  −iπ² ], [ 0,  −iπ(π² − 8)/8 ]]       J_aa = iπ κ̇/κ,   κ̇ = dκ/dc = π/4 − 2/π,   dρ_h/dc = −π/4.
```

The entries J_aa = −J_bb = iπκ̇/κ and J_ba = 0 are what det M = 1 and the exponents require; they were not used to compute J. J_ab = −iπ² is the entry the slope uses. The eigenvector equation
(n − iπ(w′ − 1))v = 0, v = (v_a, v_b), gives v_a/v_b = εJ_abκ/(2iπ(δ − κ̇ε)). The eigenfunction is Y_∞ = Y_a + β*Y_b, so v_a/v_b = 1/β*, hence

```
w′ − 1 = (δ − κ̇ε)/κ = ε J_ab β*/(2πi) = −(π/2) β* ε,        d(W/κ)/dc = −(π/2) β*.
```

**An exact check that J did not use (S4).** At W = K(c), the exact solution Y_ref is the ingoing solution for every c. The formula puts it at v_b/v_a = 2πi(K̇ − κ̇)/(J_ab κ) =
(π/4), and Y_ref = Y_nat + (π/4)Y_b: equal. **β* from the closed form (N4).** Y_∞ ∝ Y_ref ∫_ρ^∞ 2Kr²e^{−2κx}/(A R²) dr, so Y_∞ = Y_ref − Y_b/I₀ up to scale and β* = π/4 − 1/I₀.
Shooting Y_nat and Y_b to ρ = 100 and taking −Y_nat/Y_b gives β* to 1.1e-41 of the quadrature value. dW/dc follows from W = κ w′: dW/dc = κ̇ + κ(−π/2)β* = 1/I₀ − 2/π (the π/4
cancels).

**Mutations that must fail.** The sign of the perturbing operator (B₁ → −B₁) gives −0.3140 (J flips sign; S3, S4, N5 fire). The slope in the g₁ = 0 basis gives 0.9309; dropping dκ/dc, 0.5477.

## 3. The numbers

50-digit values are in `cminus1_mode_output.txt` (`table`); the first 20 digits:

| quantity | value |
|---|---|
| I₀ = ∫₀^∞ 2Kρ²e^{−2κx}/(AR²)dρ | 1.0149036099284088921 |
| 1/I₀ = d(W − K)/dc | 0.98531524591831909048 |
| β* = π/4 − 1/I₀ (quadrature; shooting to ρ = 100 agrees to 1.1e-41) | −0.19991708252087078086 |
| dW/dc | 0.34869547355073774741 |
| dκ/dc = π/4 − 2/π, dK/dc = −2/π | 0.14877839102986696654, −0.63661977236758134308 |
| d(W/κ)/dc = −(π/2)β* | 0.31402901888733605393 |
| e^{2κx} at ρ = 1 (e^{κx} = ρ^{1/2}(1 + …) at 0) | 2.4725459594 |

#212's rates (R2, `stability.py` imported) and the fit W/κ − 1 = s₁ε + s₂ε², ε = c + 1:

| c | κ | W | W/κ − 1 |
|---|---|---|---|
| −1.01 | 0.635122950037 | 0.6331476235206 | −3.110148e-03 |
| −1.001 | 0.636470902342 | 0.6362712253621 | −3.137252e-04 |
| −1.0001 | 0.636604893611 | 0.6365849043053 | −3.139986e-05 |
| −0.9999 | 0.636634649289 | 0.6366546434019 | 3.140595e-05 |
| −0.999 | 0.636768458844 | 0.6369686164035 | 3.143333e-04 |
| −0.99 | 0.638098237758 | 0.6401216240773 | 3.170964e-03 |

| fit | s₁ | s₂ | s₃ | against −(π/2)β* |
|---|---|---|---|---|
| R2, second order | 0.31403797 | 0.30408 | | 2.9e-5 |
| R2, third order | 0.31402903 | 0.30408 | 0.2657 | 2.7e-8 |
| R1, second order | 0.31405838 | 0.30411 | | 9.3e-5 |

The central difference of W at c = −1 ± 1e-4 is 0.3486954830 (R2); R1 gives 0.3487345755 (1.1e-4 higher: R1's W has rtol 1e-9 and differs from R2's by up to 8e-9). W against K = −2c/π at #212's members: W < K at the eight members
c = −10, −3, −2, −1.5, −1.1, −1.01, −1.001, −1.0001, and W > K at the eleven members c = −0.9999 … −0.01 (`table`).

## 4. Checks (derive: S1–S4; check: N0–N8)

| id | what | result |
|---|---|---|
| S1 | Y = R^{pw} e^{Wx} for symbolic c: residual ∝ (πW + 2c), zero at W = K | 0 at W = K; −5.2 at W = 2K, c = −1, ρ = 1/2 |
| S2 | c = 0: Y_K = R is ∂_b of the static member in the gauge δφ = 0 | Y = √(ρ² + 1) |
| S3 | J from residues at ρ = 0, exact in π | J_aa = iπ(π² − 8)/8 = iπκ̇/κ, J_bb = −J_aa, J_ba = 0, J_ab = −iπ² |
| S4 | the ingoing line at W = K(c) predicted by J is Y_ref | v_b/v_a = π/4 = g₁(Y_ref) − g₁(Y_nat) |
| N0 | the equation as coded equals #212's `coef_x` | 3.2e-15 (bar 1e-12) |
| N1 | closed form against the Frobenius series at ρ = ±0.01, ±0.05, ±0.1 | 5.9e-61 (bar 1e-30) |
| N2 | closed form against the integration at W = κ on ρ = ±0.3 … ±40 | 4.3e-42 (bar 1e-20) |
| N3 | formal series at +∞: e^{+κx} is R's; e^{−κx} diverges | 1.9e-122 (bar 1e-40); ratio/n = 0.797 against 0.785 (bar 5%) |
| N4 | β*: shooting to ρ = 100 against π/4 − 1/I₀ | 1.1e-41 (bar 1e-30) |
| N5 | derived d(W/κ)/dc against the fit to #212's R2 rates | 2.9e-5 second order, 2.7e-8 third (bar 1e-3) |
| N6 | 1/I₀ − 2/π against the central difference of W at c = −1 ± 1e-4 | 2.7e-8 (bar 1e-4) |
| N7 | PSLQ planted recovered, control none (both bases) | relation 6t − 9 + 10π − 12π² = 0; none, none |
| N8 | planted synthetic rates through the same fit | 2.9e-5 second order, 5.8e-7 third (bar 1e-3) |

All 13 pass in the runs of record (derive 4 of 4, check 9 of 9; 51 s and 128 s). Mutations: 15 for the script, each with the checks it must fire declared before the runs of record (README of the raw folder; the `j_g1` set was corrected from {S3, S4, N5} to {S3} after its development run, before the README); in the runs of record the fired sets equal the declared ones (15 of 15; `cminus1_mode_output.txt`).

## 5. Proposed text (the files are not touched here)

`README.md` Sec. 5, after the first paragraph (Update from #248): *Update from #262: at c = −1 the growing mode is Y = c₂ V R(ρ(Π)), Y_ref = R e^{κx} (x the tortoise
coordinate), an exact solution for every member at W = K = −2c/π (checked symbolically); the solution that is regular on the horizon and decays at infinity is Y_∞ = Y_ref − Y_b/I₀
with I₀ = 1.01490361. The rate W(c) is smooth through c = −1, W(c) = K(c) + (c+1)/I₀ + O((c+1)²), dW/dc = 0.34869547 and d(W/κ)/dc = 0.31402902 at c = −1 (first-order monodromy
perturbation theory, numerical evidence; the fit of #212's rates gives 0.314029 at third order). The two-sided 0.31 of the table is this slope.*

`horizon_sector.md` D7, D8: *The Frobenius solution of D7 with g₁ = 0 is Y_ref − (π/8) Y_b; Y_ref = R e^{κx} itself has g₁ = π/8, g₂ = 5/12 + 3π²/128 and Ψ = Y_ref e^{−κx} = R(ρ). The
Kruskal-regular solutions at W = κ form a two-dimensional space {Y_ref + βY_b}; the horizon coefficient c₂ is the same for all of them, and the condition at spatial infinity fixes β.* In the
Answer replace "the growing solution `Y = c₂ V G(Π)`" by "`Y = c₂ V R(ρ(Π))`, with G(Π) = dρ/dΠ unrelated". Normalization: here e^{2κx} = ρ(1 + …) at 0; #248's chart has Π(ρ = 1) = 1, so Π_{#248} = Π/Π(1) with Π(1) = 2.4725459594.

## 6. What this does not show

- A closed form of x(ρ) = ∫dρ/A, of I₀ or of the slope: PSLQ in the declared bases found none (below); e^{κx} is a quadrature. Whether x is elementary was not examined.
- The origin of the exact solution Y_K = R e^{Kx} for every c, and whether it is the linearization of an exact nonlinear family (b varying with v?). Nothing is said about the nonlinear fate.
- A proof of the first-order monodromy formula beyond its two exact consequences (S3, S4) and the agreement with the fit; the second-order coefficient (s₂ = 0.304, fitted only).
- Which member of {Y_ref + βY_b} a given regular Cauchy datum excites: c₂ is the same for all; β is set by what happens at large ρ (#258's slice to infinity).
- The class 𝒦 is the declared class; the divergence of the e^{−κx} series is a ratio test at n = 100, not a proof.

## Recommended tasks

1. **The second-order coefficient.** *What:* second-order monodromy perturbation theory about c = −1 (J₂ by residues, W = κ + (dW/dc)ε + ½(d²W/dc²)ε²). *Why:* the fit gives s₂ = 0.304
   (W/κ − 1 = s₁ε + s₂ε²); a derived value tests the method past first order. *Acceptance:* s₂ against the third-order fit to 1e-2; a sign mutation fails. *Depends on:* #262.
2. **#249's recount with the two-dimensional regular space at W = κ.** *What:* the count at c = −1 with Y_ref + βY_b and the outgoing condition: one growing mode, Y_∞, for all 20 members.
   *Acceptance:* #249's, plus Y_∞ by quadrature of Y_ref against the shooting solution to 1e-20. *Depends on:* #249, #262.
3. **Why Y = R e^{Kx} is exact for every c.** *What:* find the symmetry or the exact nonlinear family behind it (at c = 0 it is ∂_b; for c ≠ 0 a time-dependent b?), by linearizing a candidate
   family and by the linear solver of #232. *Acceptance:* a candidate whose linearization equals Y_K in the double-null solver to 1e-6; a perturbed candidate fails. *Depends on:* #262.
4. **A bound from the factorization.** *What:* with w = K − A/(ρR²), H + K² = 𝔄†𝔄, 𝔄 = −∂_x + w: show W < K for c < −1 and find the throat term that reverses it at c > −1 (observed at
   all 19 members of #212's table). *Acceptance:* the inequality derived or a counterexample at a member not in the table. *Depends on:* #262.

## Reproduce

From this directory, `~/Source/negative-mass/.venv/bin/python`, one thread, about 6 minutes in all:

```
cminus1_mode.py derive [--mutate=sym_pw|sym_gauge|j_sign|j_g1]                 # S1 - S4, 1 min
cminus1_mode.py check [--only=N5,N6] [--mutate=eq_P|cf_R2|cf_k|cf_xsign|inf_P|beta_g1|s_g1|s_kappa|dwdc_K|pslq_plant|fit_eps|j_sign|j_g1]   # N0 - N8, 3 min
cminus1_mode.py table                                                           # the constants, the fits, PSLQ, 3 min
cminus1_mode.py mutations DIR PREFIX                                            # fired sets of stored logs DIR/PREFIX<name>.log against the declared ones
```

`cminus1_mode_output.txt` is `derive`, `check`, `table` and `mutations` of the runs of record (`I262-mac/runs`), in that order.

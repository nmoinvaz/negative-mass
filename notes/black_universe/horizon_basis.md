# The horizon basis (c₁, c₂, c₂′) of the c = −1 black universe and the decaying partner of #262's mode (#302)

Code [horizon_basis.py](horizon_basis.py), output of record [horizon_basis_output.txt](horizon_basis_output.txt), assembled from negative-mass-raw `I302-mac/runs`; the criteria, written before the runs of record, are in
`I302-mac/README.md`. Builds on #248 ([horizon_sector.md](horizon_sector.md)), #258 ([full_slice.md](full_slice.md)), #259 ([static_c1.md](static_c1.md)), #262 ([cminus1_mode.md](cminus1_mode.md)) and
#269 ([exact_mode.md](exact_mode.md)). Status words as in AGENTS.md: S1 and S2 are *checked symbolically*, the smoothness condition of Sec. 3 is *conditional* on #258's constraint M (derived there by hand, checked
on exact pure gauge and again here on evolved data), every number from the solver is *numerical evidence*. No sentence here says a result is new (the literature was not read for it).

## Answer

**1. The three numbers at S have three exact modes, but only the static one is a smooth solution through S.** The decaying partner is `Y = c₂′ U R(ρ(Π))`, the time reverse of #262's mode. Through #248's Kruskal
solver, from exact ray data (the static mode's closed form; the growing and decaying modes in a double-null gauge built from #269's family), the static mode gives (c₁, c₂, c₂′) = (1e-3, −4e-11, −1e-12), the growing mode
c₂ = 1e-3 and the decaying mode c₂′ = −1e-3 to 1.2e-7 relative, by the horizon fits and by the slice functionals. For 8 classes and 8 random data, c₂′ by the slice functional at S equals the fit of the horizon values on V = 0 to 8.8e-8, and the two horizons' `dr_h` are reproduced by `c₁ + c₂V` (U = 0) and
`c₁ + c₂′U` (V = 0), with (c₁, c₂, c₂′) read at S, to 3.8e-8. Swapping V and U (in the functional, in the law, in the fit) fails every time it is tried.

**2. The exact growing and decaying modes cannot be put through S, as #258's and #259's tasks worded it.** For every perturbation with C² fields at S the slice constraint M gives
`∂_s² Y_N (0) = −F_Π(0) Y_N(0)` (Y_N = ∂_N Y on V + U = 0, F_Π(0) = (π/2)G₀ = 3.88387). It holds on evolved smooth data to 3e-6 (seven classes; 3.5e-4 for the narrow-pulse class Rb) (L1). The mode has `Y_N = c₂R = c₂(1 + O(s⁴))`, so it fails by `F_Π c₂`, and its
`∂_N δφ = −(π/√2) c₂ ln|s| + O(1)` at S. **The log is the solution's, not the slice gauge's:** the condition is on the gauge-invariant `Y`, and no gauge (smooth or double-null) removes the log
(Sec. 3). #269's family in a double-null gauge has the same log with the same coefficient (L2). Sentences of the repository that say the mode is smooth across S are listed in Sec. 6; none is edited here.

```
        V = 0, U > 0                          U = 0, V > 0  =  H⁺ :  Y_h = c₁ + c₂ V        (mode Y = c₂ V R)
                     ╲    future interior    ╱
                      ╲                     ╱
      left exterior    ╲                   ╱    right exterior (V > 0 > U)
                        ╲                 ╱
        ───────────────────────  S  ───────────────────────   slice V + U = 0;  c₁ = Y(S): static mass mode
                        ╱                 ╲
                       ╱                   ╲
      U = 0, V < 0    ╱   past interior     ╲   V = 0, U < 0  =  H⁻ :  Y_h = c₁ + c₂′ U     (mode Y = c₂′ U R)
```

Y = c₂VR is smooth on U = 0 for V > 0 and has `ln V` in δφ_N as V → 0; Y = c₂′UR is its image under T, smooth on V = 0 for U < 0, with `ln|U|` as U → 0.

| acceptance | result |
|---|---|
| I1: c₂′ of the slice equals the mode's to 1e-6 | decaying mode, ray data on U ≤ −0.05: c₂′ = −1.000000085e-3 (fit on V = 0), −1.000000083e-3 and −1.000000123e-3 (functional d_U at U₀ = −0.5, −1); exact −1e-3 (P4). Slice at S against the fit on V = 0: 8.8e-8 over 16 data (H2). **Not met through S itself:** the exact mode has no smooth representation there (Sec. 3) |
| I2: a random datum's (c₁, c₂, c₂′) is reproduced by the three modes to 1e-6 | `dr_h(V, 0) − (c₁ + c₂V)` and `dr_h(0, U) − (c₁ + c₂′U)` over every node of both horizons: 3.8e-8 of max(max\|dr_h\|, 1e-3), 8 classes and 8 random data (H3); the modes' own horizon values: P2 (4e-8 of c₁), P3, P4 (1.2e-7) |
| I3: swapping V and U fails | `fn_swap` (d_U reads d_V): P4 reads 1.74e-3 and 4.27e-3 for −1e-3, H2 off by 2.7e2, H3 by 0.67; `law_x` (fit against −U): H2 off by 2.0, P4 fires; `s_swapfn` (#258's c₂, c₂′ swapped): S3 fires |
| C: the log, the solution's or the gauge's | the solution's (Sec. 3): S1, P1, L1, L2 |
| each check has a `--mutate=` with a computed effect | 14 of 14 fired sets equal the declared ones (Sec. 5) |

## 1. The decaying mode (S2, S3)

The Kruskal master equation of #232 (`full_slice.ModeData.M_at`) is `2κ²Π u_VU + ½V̂u = M`, u = Y/ρ. For `Y = V f(Π) ρ` and for `Y = U f(Π) ρ`, with f = R/ρ, `u_VU = −V (Πf)″` and `−U (Πf)″`: the same bracket.
With d/dΠ = (A/2κΠ) d/dρ it becomes `κA f′ + ½A(Af′)′ = ½V̂f`, with no Π left. S2 evaluates it symbolically (sympy: 0) and at 50 digits at six points (7.1e-49). So `Y = c₂′ U R(ρ(Π))` solves the equation iff `Y = c₂ V R(ρ(Π))` does.
It is #262's mode under `T: (V, U) → (−U, −V)`, `t → −t`: `c₂′ = −c₂`, and at fixed ρ it decays as e^{−κt}. On the slice V + U = 0, T leaves (δr, δφ, δr′) and flips (n, δφ_N); the slice code's functionals
`c₂ = ½(δr′ + n)`, `c₂′ = ½(n − δr′)` swap under it. S3: on datum G of #258 and three random compact data, `c₂(TD) = −c₂′(D)`, `c₂′(TD) = −c₂(D)`, `c₁(TD) = c₁(D)` to 0.0 (the data are scaled by −1).

## 2. The three exact modes through the solver (P2, P3, P4)

**Static** (#259's closed-form ray data, δL = 2e-3, rectangle of #248): c₁ = 9.999999965e-4 at S, 9.999999951e-4 (fit on U = 0), 9.999999966e-4 (fit on V = 0); c₂ = −4.0e-11 and −3.1e-11, c₂′ = −1.3e-12 and −9.6e-13 (slice functional, fit).

**Growing and decaying**, exact double-null data from #269's family. In the ingoing chart `ds² = −F dv² + 2 dv dρ`, the perturbation is δr = (β/2) R e^{κv}, δF = −β(1 + κρ)e^{κv}, δφ = 0 (c = −1, K = κ, e^{κv} = V). With `v = t + x`, `u = t − x`, `V = e^{κv}`,
`U = −e^{−κu}` the background is `−A dv du`, and the perturbation is `δg_vv = β V (1 + κρ)` only. The gauge vector ξ^v = 0, `∂_v ξ^u = δg_vv/A` removes it:

```
ξ^u(V, U) = (β/κ) ∫ (1 + κρ)/A  dV′  at fixed U (zero on a base ray);   ∂_u ξ^u = −κU ∂_U ξ^u,   ∂_U ξ^u = −(β/κ) ∫ V′ W_ρ G dV′,   W = (1 + κρ)/A,  G = dρ/dΠ
δr = ½ β R V − ξ^u Aρ/(2R),     δφ = −ξ^u A/(√2 R²),     δF = ½ A′ ξ^u − ∂_u ξ^u        (double-null gauge, Kruskal chart; ξ^u carries the factor β)
```

The data rays are the base ray (ξ = 0, closed form) and one line of constant U (quadrature). The decaying mode is the image under T, with the base ray V′ = 3.0 the image of the data ray U = −3.
The evolved Y is compared with `c₂ V R(ρ(Π))` and `c₂′ U R(ρ(Π))` on a sub-grid.

| mode | amplitude | Y error / max at h = 0.005, 0.0025, 0.00125 | last ratio | c by the fit | c by the functional at the crossings (Richardson of 0.0025, 0.00125) | exact |
|---|---|---|---|---|---|---|
| growing, V ≥ 0.05 | β = 2e-3 | 3.64e-3, 9.11e-4, 2.28e-4 | 4.00 | c₂ = 1.000000004e-3 | 1.000000009e-3 (V₀ = 0.5), 1.000000121e-3 (1), 1.0000025e-3 (2) | 1e-3 |
| decaying, U ≤ −0.05 | β = 2e-3 | 1.16e-2, 2.90e-3, 7.25e-4 | 4.00 | c₂′ = −1.000000085e-3 | −1.000000083e-3 (U₀ = −0.5), −1.000000123e-3 (−1), −1.000000264e-3 (−2) | −1e-3 |

c₁ of the fit is 1.8e-11 (growing) and −2.9e-10 (decaying). The functional's error grows with the lever arm (2.5e-6 at V₀ = 2, outside the bar, which names V₀ = 0.5, 1). The growing mode through the solver is smooth for V > 0,
and has the log of Sec. 3 at V → 0 (it was not evolved into V ≤ 0). `gauge_dphi` and `gauge_dF` (dropping the δφ or the ∂_uξ^u term) make the data a non-solution: P3 and P4 fire.

## 3. The log, and why it is the solution's (S1, P1, L1, L2)

**The condition.** On V = s, U = −s, Π = s², with n = ∂_N δr, z = ∂_N δφ, ε = −1, #258's M is `n′ + sF_Π n = −s r_Π δF_N − ε r s φ_Π z`. With `Y = δr − gδφ`, `g = ρR/√2` (invariant under every linear diffeomorphism, #248 D3a),
`n = Y_N + g z` and `Y_N = ∂_N Y`:

```
Y_N′ + s F_Π Y_N = −g z′ − (g′ + ε r s φ_Π + s F_Π g) z − s r_Π δF_N                       (M′)
g = a s² + …,  a = G₀/√2;   g′ + ε r s φ_Π = 2√2 s G ρ²/R = O(s⁵) (exact);   s F_Π g = O(s³);   r_Π = O(s²);   F_Π(0) = G₀/κ = (π/2)G₀
```

(S1, sympy: A′(0) = 2κ, A″(0) = −2, ρ = G₀Π + ½G₁Π² with G₁ = A″G₀²/2κ, F_Π(0) = −G₁/G₀, the identity above, ζ = −F_Π/a = −π/√2.) For fields C² at S the right side is O(s²), so the coefficient of s in the left side vanishes:

**`Y_N″(0) = −F_Π(0) Y_N(0)`**. L1 tests it on evolved smooth solutions, Y_N″(0) from 5-point stencils on the evolved array, Richardson of h = 0.00125, 0.000625:

| class | A0 | A1R | A1L | L | G | Lm | Rb | Gb |
|---|---|---|---|---|---|---|---|---|
| Y_N″(0) / (−F_Π Y_N(0)) | 1.000000 | 1.000003 | 1.000003 | 1.000000 | 1.000000 | 1.000000 | 1.000345 | 1.000032 |
| Y_N(0) | 2.98e-2 | 1.38e-3 | 1.68e-2 | −9.94e-4 | −9.94e-4 | −2.73e-3 | −2.36e-4 | −2.96e-3 |

(Rb has a pulse of width 0.03; its value at 0.00125 alone is 0.834.) The exact mode, `Y_N = c₂R(ρ(s²))`, has `R″(0) = 1.5e-9` against `−F_Π = −3.88387`: it violates the condition by `F_Π(0) c₂`. With δF_N = 0 the equation for z is first order,
`g z′ = −(Y_N′ + sF_Π Y_N) − (…)z`, so `z = −(F_Π/a) c₂ ln|s| + O(1) = −(π/√2) c₂ ln|s| + O(1)`: integrating it inward from s = 0.3 (the E9 data of #258) gives the slope ratio 1.0015386, 1.0000243, 1.0000003 at s = 10⁻², 10⁻³, 10⁻⁴ (P1).
(`Θ_gen` of the data is 1.5e-12 of its largest term: the mode carries no mass, which is what #258's M = 0 initial value says; with E9's initial value left at 0 it is 0.31.)

**The same log from #269's family.** In the double-null gauge of Sec. 2, ξ^u ≈ −β ln V/(2κ²G₀U) near V = 0, so `δφ ≈ −(β/(√2κ)) V ln V` and `∂_Vδφ = −(β/(√2κ)) ln V = −(π/√2) c₂ ln V` (c₂ = β/2, κ = 2/π).
L2 evaluates the exact fields at V = 10⁻³ … 10⁻⁶; the slope of ∂_Vδφ against ln V over −(π/√2)c₂ is

| U, base ray V_b | V = 1e-3 | 1e-4 | 1e-5 | 1e-6 |
|---|---|---|---|---|
| −0.5, 0.05 | 1.01530 | 1.00243 | 1.00033 | 1.00004 |
| −0.5, 1.0 | 1.03822 | 1.00473 | 1.00056 | 1.00006 |
| −1.0, 0.05 | 1.03155 | 1.00497 | 1.00068 | 1.00009 |
| −1.0, 1.0 | 1.10263 | 1.01212 | 1.00139 | 1.00016 |

Two formalisms (the slice constraint with the background at S; the ingoing family and a gauge vector) give one number. The agreement uses A″(0) = 2c = −2, i.e. it is specific to the member.

**Why it is not the gauge's.** (1) The condition is on `Y`, which every gauge shares. (2) A representation with smooth (δr, δφ, δF) in any gauge goes to a smooth double-null representation through a smooth ξ
(`∂_V ξ^U = h_VV/2g_VU`, `∂_U ξ^V = h_UU/2g_VU` are integrals of smooth functions), and `h′ = h + L_ξ g` is smooth only for smooth ξ. (3) Two double-null representations differ by `ξ^V(V), ξ^U(U)`;
these shift δφ by `−Φ′G (U ξ^V(V) + V ξ^U(U))`. The V ln V coefficient above does not depend on U (L2 at U = −0.5, −1), and a term `U ξ^V(V)` with ξ^V = V ln V has a coefficient proportional to U and puts `ln V` into
δF. So no smooth representation exists. The condition L1 verifies is the one that excludes it. What is not shown: the hand argument (2) assumes the representation is a tensor field smooth in the Kruskal chart; a gauge in which the fields are only continuous at S was not examined.

**What survives.** The horizon laws (δr_VV = 0 on U = 0, #248 D2; the time reverse on V = 0, H1), c₂ and c₂′ as free numbers of smooth data (N11, H2), the growth rate κ, and the numbers of #248, #258, #259 (all from smooth data; the
slice code of #259 holds smooth data, and the exact mode is not among them). What changes is the reading of `Y = c₂VR` as "the regular growing solution": the smooth data excite `c₁Y_static + c₂VR + c₂′UR + Y_rem`, where
`Y_rem` vanishes on both horizons (H3) and carries the `Y_N` condition above (D6's free coefficients y₁₂ − y₂₁). That `Y_rem` decays, so that the late growth is c₂VR, was not tested.

## 4. The basis on both horizons (H1, H2, H3)

Rectangle V ∈ [−0.1, 3], U ∈ [−2, 0.06]; the past strip is wide enough for a pulse on the ray U = −2 to reach V = 0 (Π = −0.2 at the corner). Slice readings at S and fits, Richardson of h = 0.00125, 0.000625:

| datum | c₁ (S) | c₂ (S) | c₂′ (S) | c₂′ (fit on V = 0) |
|---|---|---|---|---|
| A0 | 4.067796e-3 | 2.837763e-2 | 1.418882e-3 | 1.418882e-3 |
| A1R | 1.337668e-4 | 1.337668e-3 | 3.853927e-5 | 3.853927e-5 |
| A1L | 2.675336e-3 | 1.541571e-2 | 1.337668e-3 | 1.337668e-3 |
| Rb (pulse on the ray U = −2, V ≈ −0.05) | −5.7e-13 | −1.9e-11 | −2.359099e-4 | −2.359099e-4 |
| L (pulse on the ray V = V₀ at U = 0) | −2.1e-17 | −9.942330e-4 | −8.2e-15 | −2.3e-14 |
| rand0 … rand7 | −4.1770e-3 … 2.0663e-3 | −2.4005e-2 … 1.4914e-2 | −1.9589e-3 … 8.2351e-4 (\|c₂′\| ≥ 5.0e-5) | equal to 6 digits |

(Full table in the output file.) The pulse on the left ray (L, Lm) feeds c₂ and not c₂′; the pulse in the strip on the bottom ray (Rb) feeds c₂′ and not c₂. H1: the residual of the linear fit of `dr(0, U)` falls by 3.98 … 4.00 per halving
(A0 5.7e-6, 1.4e-6, 3.6e-7; A1R 1.7e-5, 4.3e-6, 1.1e-6; Rb 3.3e-5, 8.2e-6, 2.1e-6 at h = 0.005, 0.0025, 0.00125). The classes L, G and Lm have `max|dr| ≤ 1e-8` on V = 0: only the truncation seed, ratios not scored.
H2: the functional at S, U₀ = −0.5, −1 equals the fit to 8.8e-8 of max(|c₂′|, 1e-5) (worst of the 16; Rb and Gb give 8.2e-8 at U₀ = −1). H3: 3.8e-8.

## 5. Checks and mutations

| id | what | result |
|---|---|---|
| S1 | series at S, F_Π(0) = G₀/κ, ζ = −π/√2, identity of g′ + ε r s φ_Π | exact |
| S2 | Y = V R and Y = U R solve the master equation | 0 (sympy); 7.1e-49 |
| S3 | time reversal of the slice code, 4 data | 0.0 |
| P1 | E9 data: Θ_gen 1.5e-12, M 2.8e-11; z slope 1.0000003 at s = 1e-4 | pass |
| L1 | Y_N″(0) = −F_Π Y_N(0), 8 classes; exact mode R″(0) = 1.5e-9 | 4e-6 (3.5e-4 for Rb) |
| L2 | log coefficient of the exact double-null δφ | 1.00004 … 1.00016 at V = 1e-6 |
| P2 | static mode (c₁, c₂, c₂′) | c₂, c₂′ ≤ 4e-8 of c₁; c₁ to 4e-9 |
| P3, P4 | growing and decaying modes: Y error ratio; c₂, c₂′ | 4.00; 1.2e-7 |
| H1 | law on V = 0, 8 classes | residual ratios 3.98 … 4.00 |
| H2 | c₂′ slice against fit, 8 classes + 8 random | 8.8e-8 |
| H3 | both horizons' dr_h from (c₁, c₂, c₂′) at S | 3.8e-8 |

All 12 pass in the runs of record (`derive` 3 of 3, 2 s; `check` 9 of 9, 135 s, 715 MB; code d1f849b). 14 mutations, each with the checks it must fire declared before the runs of record (raw README; computed in
development runs `mut_d1`, `mut_d2`, `mut_d3`): in the runs of record the fired sets equal the declared ones (14 of 14; `horizon_basis_output.txt`).

| mutation | fires | effect in the run of record |
|---|---|---|
| s_zeta, s_kappa, s_noflip, s_swapfn | S1, S2, S3, S3 | F_Π doubled; κ·1.01; (n, φ_N) not flipped; #258's `fn_sign` |
| fn_swap | P2, P4, H2, H3 | P4 reads 1.74e-3, 4.27e-3 for −1e-3; H2 2.7e2; H3 0.67 |
| fn_aT | P2, P4, H2, H3 | P4 reads −2.74e-3, −5.27e-3; H2 2.7e2; H3 0.67 |
| law_c1 | H3 | 5.5e-2 (c₁ missing from the law on V = 0) |
| law_x | P4, H2 | H2 2.0 (the fit against −U) |
| law_off | P2, P4, H1, H2, H3 | the line V = 0.025: H2 1.0, H3 1.7e-2 |
| gauge_dF, gauge_dphi | P3, P4; L2, P3, P4 | the exact data are not a solution; L2's slope 0 |
| alpha_off | P1 | Θ_gen 0.31 of its largest term, M 13 of its scale before the fix |
| sc_off | L1 | F_Π(0) = 0 in the condition |
| zeta_F | P1 | slope ratio 2.0 |

Blind spots, kept: `law_c1` does not fire P4 or H2 (c₂′ does not depend on c₁); `fn_swap` does not fire P3; L1 and L2 are not fired by the reading mutations (they do not use the functionals).
**What is shared and what is not.** P3, P4 and L2 share the background ρ(Π) of the chart (through jets) and the closed form of A with the solver; the exact fields are algebra on #269's family (a different formalism from the
solver's equations: the solver evolves its own linearized system, and the exact data satisfy it to the rate 4.00). S1 and L2 give the log coefficient from M (slice, background at S) and from the family (gauge vector); S2 is sympy
and shares only A. H1 to H3 use #248's solver and #248's functionals, read at S and on V = 0, a line #248 did not read.

## 6. Proposed text and discrepancies (the files are not touched here)

Sentences that call the mode smooth or regular through S, with what they should say (Y is smooth; the perturbation with that Y is not; smooth data have Y = c₁Y_static + c₂VR + c₂′UR + Y_rem):

| file, line | now | change |
|---|---|---|
| `paper/4-black-universes/main.tex` 399–401 | "*Regular.* The growing solution is Y = c₂VR(ρ(Π)) … analytic in Π: smooth across S and into the interior. … the series at W = κ has no logarithm." | say *the gauge-invariant Y*; add that smooth data satisfy ∂_s²Y_N(0) = −F_Π(0)Y_N(0), which Y = c₂VR violates, so the perturbation with this Y has ∂_Nδφ ∼ −(π/√2)c₂ ln\|s\| at S in every gauge; the series statement is about Y |
| `horizon_sector.md` 12–13 | "is smooth across S and through the interior (ρ < 0). So #232's horizon sector belongs to physical initial data: … a growing solution of rate κ, regular on the whole slice, excited by regular data …" | "Y is a smooth function across S …; regular data have Y_h = c₁ + c₂V on H⁺ and c₁ + c₂′U on H⁻, c₂ free, and excite a growth of rate κ with coefficient c₂ (#302: the perturbation with Y = c₂VR is not smooth at S)" |
| `horizon_sector.md` 145, 155 | "a solution that is regular across S and through the interior", "regular across the bifurcation surface and the interior" | the same qualification (155 is the text proposed for README §5; the applied sentence is `README.md` line 13, next row) |
| `README.md` 13 | "the solution Y = c₂ V G(Π) is smooth across the bifurcation surface and into the interior" | "the function Y = c₂ V R(ρ(Π)) …, not the Y of a smooth perturbation (#302)" |
| `horizon_lemma.md` 81 | "the c = −1 mode Y = c₂VΨ(Π) is smooth across S and grows at …" | "is a smooth function across S (not the Y of a smooth perturbation, #302) and grows at …" |
| `horizon_literature.md` 41 | "… G(Π) is smooth across S." | the same |
| `cminus1_mode.md` 18, 35, 174; `horizon_sector.md` 175, 222 | "Kruskal-regular class", "Kruskal-regular solutions at W = κ" | add "(as solutions of the Y equation; #302)" |
| `full_slice.md` Sec. 8, bullet "A planted exact solution through S" | "the equation for δφ has g ∼ s² in front of its highest derivative, singular at S" | the singular equation is a symptom: no smooth δφ_N exists (Sec. 3); the exact mode is planted on [0.3, 1.6] only |
| `full_slice.md` Sec. 9 task 2, `static_c1.md` Sec. 9 task 2 | "put the exact mode through S", "put the time reverse … through the slice code and read c₂′ = ½(n − δr′)" | *Done in #302*: the decaying partner, the law and the functional on H⁻; the exact modes cannot pass through S |
| `horizon_sector.md` "What this does not show", bullet "The decaying partner Y ∝ U … is not computed separately" | | *Done in #302* (H1, H2, P4) |

`horizon_sector.md`, after "What c₁ is": *The horizon laws are Y_h(V) = c₁ + c₂V on H⁺ (U = 0) and Y_h(U) = c₁ + c₂′U on H⁻ (V = 0), with c₂′ = ∂_U δr(S) = ½(n − δr′). The three numbers have the exact modes c₁Y_static (smooth), c₂VR(ρ(Π))
and c₂′UR(ρ(Π)) (solutions of the master equation; singular at S as perturbations: ∂_Nδφ ∼ −(π/√2)c₂ ln\|s\|), whose horizon values are c₁ + c₂V on H⁺ and c₁ + c₂′U on H⁻ (#302).*

`static_c1.md`, recommended task 2: *Done in #302.* **Paper IV**, after the paragraph "The static coefficient (#259)":

```
\paragraph{The decaying partner (\#302).} The time reverse of the growing solution is $Y=c_2'UR(\rho(\Pi))$, and $c_2'=\partial_U\delta r(S)=\tfrac12(n-\delta r')$. The horizon laws are $Y_h=c_1+c_2V$ on $H^+$ and $Y_h=c_1+c_2'U$ on $H^-$; for
evolved smooth data both hold, with $c_1,c_2,c_2'$ read at $S$, to $4\times10^{-8}$ of the amplitude (8 classes, 8 random data). The exact modes are not smooth perturbations through $S$: smooth data satisfy
$\partial_s^2Y_N(0)=-F_\Pi(0)Y_N(0)$ (to $10^{-6}$ on evolved data), and $Y=c_2VR$ has $Y_N=c_2R$, so $\partial_N\delta\phi=-(\pi/\sqrt2)\,c_2\ln|s|+O(1)$ there, in every gauge.
```

## 7. What this does not show

- That `Y_rem` of smooth data decays. The growth of smooth data at rate κ with coefficient c₂ is #232's and #248's observation (N-checks to V = 3); the split into c₂VR and a remainder was not evolved to late times.
- A gauge-free proof. Sec. 3's argument (2) assumes tensor fields smooth in the Kruskal chart. The condition rests on #258's M (checked on pure gauge there and on evolved data here, 1e-6), not re-derived from the Einstein equations.
- The other members. A″(0) = 2c is used; at c ≠ −1 the horizon is not at the area minimum and ρ_h ≠ 0 (#301), so S1 does not carry over.
- The static mode through S as a *planted* case with the three readings at once: P2 does it (through the solver); the growing and decaying modes are tested away from S only.
- The nonlinear fate: #269's family for β ≠ 0 is singular at V = 0 (at linear order, ln V in δφ); whether the nonlinear solution is smooth there was not examined.
- Literature: nothing was read for the log or for the horizon laws on H⁻.

## 8. Recommended tasks

1. **The remainder of smooth data.** *What:* evolve classes A0, Rb and random data to V ≈ 10 (the chart of #232), subtract c₁Y_static + c₂VR + c₂′UR, and test that Y_rem decays. *Why:* Sec. 3 changes the reading of "the mode excited by smooth data". *Acceptance:* `|Y_rem| / |c₂VR|` falls by 4 ± 0.5 per halving of h and below 1e-2 at V = 10; a mutation dropping the F_Π condition fails. *Depends on:* #302, #232.
2. **D6's free coefficients from the condition.** *What:* the formal series of Y at S with the condition ∂_s²Y_N(0) = −F_ΠY_N(0) imposed: which combination of y₁₂, y₂₁ it fixes in terms of c₂ + c₂′. *Acceptance:* the coefficient against the fit of evolved Y to 1e-6; a wrong F_Π fails. *Depends on:* #248, #302.
3. **A Θ for slices that miss S** (#259 task 3, unchanged): c₁ from one point of a slice like #232's pulses'. *Depends on:* #259, #232.

## 9. Reproduce

From this directory, `~/Source/negative-mass/.venv/bin/python`, one thread:

```
horizon_basis.py derive [--mutate=s_zeta|s_kappa|s_noflip|s_swapfn]                                   # 3 checks, 2 s
horizon_basis.py check [--only=P1,L1,L2,P2,P3,P4,H1,H2,H3] [--cache=DIR] [--mutate=fn_swap|fn_aT|law_c1|law_x|law_off|gauge_dF|gauge_dphi|alpha_off|sc_off|zeta_F]   # 9 checks, 135 s cold, 15 s with the cache
horizon_basis.py table [--cache=DIR]                                                                # the numbers of this note
horizon_basis.py mutations DIR PREFIX                                                               # fired sets of the stored logs DIR/PREFIX<name>.log against the declared ones
```

`horizon_basis_output.txt` is `derive`, `check`, `table` and `mutations` of the runs of record (`I302-mac/runs`, code commit d1f849b), in that order. The cache holds readouts of evolved rectangles; no mutation changes an evolution.

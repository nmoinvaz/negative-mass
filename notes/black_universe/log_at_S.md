# The logarithm of the c = −1 growing mode at S, rederived from the field equations (#321)

Code [log_at_S.py](log_at_S.py), output of record [log_at_S_output.txt](log_at_S_output.txt), assembled from negative-mass-raw `I321-mac/runs`; the criteria, written before the runs of record, are in `I321-mac/README.md`. Checks #302
([horizon_basis.md](horizon_basis.md) Sec. 3) with code that does not import or reuse `horizon_basis.py`, the slice constraint M or the Kruskal solvers. Status words as in AGENTS.md: the series results are *checked symbolically*
(exact linear algebra over ℚ(κ) on truncated formal series, which is not a proof for all orders), the numbers of Route B are *numerical evidence*. No sentence here says a result is new (the literature was not read for it).

## Answer

**1. #302's coefficient is reproduced exactly by a route that shares only the background with it.** Route A writes my own Ricci tensor of the double-null metric, linearizes it about the member, and solves the formal series of the five
equations with an ln V factor allowed in every field. With Y = c₂VR(ρ(Π)) the coefficient of V ln V in δθ (θ = φ/√2) is q = −c₂/κ for every solution (gauge freedom included), so `∂_Nδφ = −(√2/κ) c₂ ln|s| = −(π/√2) c₂ ln|s|` at
κ = 2/π: #302's number, with no difference. Route B builds a double-null chart from the outgoing rays of #269's exact family and integrates them at 30 digits: the slope of `∂_Vδφ` against `ln V` over `−(π/√2)c₂` is
1.0000018, 1.0000071, 1.0000009, 1.0000018 at V = 1e-8 for four (V_b, U) (bar 1e-4), with c₂ = 1/2 per unit β read from the family.

**2. A smooth solution with δr = c₂V on H⁺ and δr = 0 on H⁻ exists, and its Y differs from c₂VR at third order.** For every smooth solution (formal series; double-null gauge, and the general smooth gauge with δg_VV and δg_UU) the Y
coefficients obey `y₂₁ = (F_Π/2) y₁₀`, `y₁₂ = (F_Π/2) y₀₁` and `y_i0 = y_0i = 0` for i ≥ 2 (Y = Σ y_ij Vⁱ Uʲ, F_Π(0) = 1/κ in my normalization Π = ρ + …). So with c₁ = c₂′ = 0,
`Y = c₂ (V + (F_Π/2) V²U + …)`, while `VR(ρ(Π)) = V + O(V³U²)` has y₂₁ = 0. #302's condition `Y_N″(0) = −F_Π Y_N(0)` is the sum of the two relations; the relations hold in each branch separately. Without logarithms there is no
solution with Y = c₂VR, and none with Y = c₂′UR; with t in `Y = V + t V²U + …` the log coefficient is `q = 2t − F_Π` and a smooth solution exists exactly at t = F_Π/2.

**3. The planted solution has no logarithm and the control follows the prediction.** The static mass mode (Y = c₁R) has q = 0 for every solution of the same system (S4), and in Route B, in its own Kruskal chart, a slope of 6e-8 and 7e-7
of the log scale at V = 1e-8 (B3). In the chart of the background (V = e^{κv}) it has a log proportional to U, the gauge `ξ^V = (δκ/κ) V ln V`, equal to a closed form to 5e-6 (B4). Members with A″(0) = −2/b² and κ independent
(b = 2, 3/2, 1/2; also G₀ = 3, 5/2) have `q = A″(0) c_h/(2κ)` exactly (S5).

| acceptance of the issue | result |
|---|---|
| I1: the coefficient of ln\|s\| (ln V) is −(π/√2)c₂ within 1 % | Route A: q = −c₂/κ exactly, i.e. −(π/√2)c₂ at κ = 2/π (S3, double-null and general smooth gauge); Route B: ratio 1 + 1.7e-6 … 7.1e-6 at V = 1e-8 (B2, X1); a smooth solution exists (S1b) and differs from c₂VR at V²U (S1, S6) |
| I2: a planted smooth solution, the static mass mode, has no logarithm through the same code | S4 (series: q = 0, no spread); B3 (rays: 6e-8, 7e-7 of the scale; bar 1e-6); B4 (the log it does show in the wrong chart is the closed-form gauge term) |
| I3: a control with A″(0) changed moves the coefficient as predicted | S5: q = A″(0)/(2κ) c_h for A″(0) = −2, −1/2, −8/9, −8, equal as rational functions of κ; S6: q = 2y₂₁ − F_Πy₁₀ |
| each check has a `--mutate=` with a computed effect | 18 of 18 fired sets equal the declared ones (Sec. 7) |

## 1. The equations and the series (derive: D1–D4)

My Ricci tensor (sympy, brute force) of `ds² = −2Ω dV dU + h_VV dV² + h_UU dU² + r² dΩ²` with `R_ab = −2θ_aθ_b + 𝒱 g_ab` (ε = −1, θ = φ/√2) and `−2□θ = 𝒱_θ`. The linearization is d/dη at
`r₀ + ηδr, Ω₀ + ηδΩ, θ₀ + ηδθ, ηh_VV, ηh_UU` (h_VV = h_UU = 0 is double-null gauge). The background is the member (b, κ) with the horizon at the throat, R² = b² + ρ², θ = arctan(ρ/b),
A = 1 − R²/b² + κ(ρ + R² arctan(ρ/b)/b), A′(0) = 2κ, A″(0) = −2/b², Π = e^{2κx} with ρ = G₀Π + …, Ω₀ = A/(2κ²Π) (from `−A du dv`); b = 1, κ = 2/π is the c = −1 black universe, and other κ are local solutions with the
potential of that κ. The potential is not assumed: D1 solves E_θθ for 𝒱 and checks E_tt, E_ρρ and the scalar equation of the static metric symbolically (b, κ symbols), and that 𝒱 = (1 − (ρA)′)/R². D2: the background series solves
all five double-null equations (no residual at degree ≤ 8; first nonzero at 9 or 10, the truncation). D3: pure-gauge perturbations (`ξ^V, ξ^U` arbitrary polynomials, δΩ = ξ·∂Ω₀ + Ω₀(∂_Vξ^V + ∂_Uξ^U), h_VV = −2Ω₀∂_Vξ^U,
h_UU = −2Ω₀∂_Uξ^V) solve the five linearized equations: the operator, the Kruskal normalization and the Lie derivative are tested against one another. D4: `Y = δr − g δθ`, g = r₀_Π/θ₀_Π, vanishes on them.

The solver takes fields as polynomials in V, U of total degree ≤ N (times 1 or ln V), expands the equations to degree N − 2 and solves exactly over ℚ(κ) (`DomainMatrix`). Questions put to it: the relations among the y_ij valid for
all solutions (left nullspace of the map solutions → Y coefficients); whether a solution with prescribed Y exists; the log coefficient; and whether that coefficient changes over the nullspace.

## 2. Smooth solutions (S1, S1b, S2)

Relations of the double-null system (N = 7, nullity 33), degree ≤ 4, as the output prints them:

| relation | reads |
|---|---|
| y₁₀ − 2κ y₂₁ = 0, y₀₁ − 2κ y₁₂ = 0 | y₂₁ = y₁₀/(2κ) = (F_Π/2) y₁₀, y₁₂ = (F_Π/2) y₀₁ |
| y₂₀ = y₃₀ = y₄₀ = 0, y₀₂ = y₀₃ = y₀₄ = 0 | δr = c₁ + c₂V on H⁺ and c₁ + c₂′U on H⁻ (#261, #248 D2) |
| y₀₀ + y₁₁/κ − 2y₂₂ = 0 | one more at degree 4 (the static mode R satisfies it) |

The same two relations and the horizon laws hold for the general smooth gauge, with δg_VV and δg_UU as unknowns (N = 5, nullity 67; S1g). With Y_N = (∂_V + ∂_U)Y on V = s, U = −s,
Y_N(0) = y₁₀ + y₀₁ and Y_N″(0) = 2(3y₃₀ − y₂₁ − y₁₂ + 3y₀₃), so the two relations give `Y_N″(0) = −F_Π(0) Y_N(0)`, #302's condition (S1: in the span; checked with F_Π read from my Ω₀). `y₂₁ = 0`, which Y = VR has, is not in the span.

S1b: with only y₀₀ = y₀₁ = 0, y₁₀ = 1 imposed, the system has a solution (nullity 30). Its Y is `V + (F_Π/2)V²U + …`; the first difference from `VR(ρ(Π)) = V + R₂V³U² + …` is the V²U term, with coefficient
F_Π(0)/2 = 1/(2κ) (my normalization; in general `y₂₁/y₁₀ = F_Π(0)/2`). S2: with the smooth ansatz and Y = VR (all degrees ≤ 7, or only degrees ≤ 3) there is no solution; likewise for Y = UR; with Y = R there is one.

The relations are per branch, which is more than #302 tested: its L1 reads the sum y₁₂ + y₂₁. In `horizon_basis.md` Sec. 3 ("D6's free coefficients y₁₂ − y₂₁") and its recommended task 2 ("which combination of y₁₂, y₂₁ the condition fixes"):
for smooth solutions of the full system both are fixed, y₂₁ by c₂ and y₁₂ by c₂′; the combination y₁₂ − y₂₁ is `(F_Π/2)(y₀₁ − y₁₀)`.

## 3. The logarithm (S3, S6)

With ln V allowed in every field, Y = c₂VR has solutions (nullity 31 in double-null gauge; 103 with h_VV, h_UU). The coefficient of V ln V in δθ is `−1/κ` per unit c₂ (c₂ = ∂_Vδr on H⁺) in every one of them: adding any
nullspace vector (pure gauge with or without logs, the constraint-free parts, the top-degree unknowns) leaves it unchanged. The other log terms of δθ of degree ≤ 1 are zero. Along the slice V = s, U = −s,
∂_Nδθ = q(ln s + 1) + O(s ln s), so the coefficient of ln s is q, and in φ = √2θ at κ = 2/π, `−(√2/κ)c₂ = −(π/√2)c₂`. A logarithm in δΩ (V ln V) is present in the particular solution but not invariant (it moves with a gauge ξ^V ∝ V ln V).

S6 ties the two results: prescribing Y = V + t V²U (degrees ≤ 3), t = 0, ½, 1, 2 times F_Π/2, gives q = −F_Π, −F_Π/2, 0, F_Π, i.e. `q = 2y₂₁ − F_Π y₁₀`, and a smooth solution only at q = 0.
The logarithm is the violation of the smooth relation, with coefficient proportional to it.

## 4. Planted solution and control (S4, S5, B3, B4)

**S4.** Y = c₁R (the static mass mode, δr = (δL/2)R in the static gauge) through the same log ansatz: a solution exists, q = 0 for every solution, and the particular solution has no log term in δθ.

**S5.** Members with other (b, G₀), κ a symbol; A″(0) read from my series:

| b | G₀ | A″(0) | q (series) | A″(0)/(2κ) |
|---|---|---|---|---|
| 1 | 1 | −2 | −1/κ | −1/κ |
| 2 | 1 | −1/2 | −1/(4κ) | −1/(4κ) |
| 3/2 | 1 | −8/9 | −4/(9κ) | −4/(9κ) |
| 1/2 | 1 | −8 | −4/κ | −4/κ |
| 1 | 3 | −2 | −1/κ | −1/κ |
| 2 | 5/2 | −1/2 | −1/(4κ) | −1/(4κ) |

(c_h = ∂_Vδr on H⁺ = 1 by construction, `Y = V R/R(0)`.) q follows A″(0) and κ, and not the normalization G₀ of Π. For the asymptotically flat member κ = 2/(πb), so `q = −π/(2b)`: the coefficient scales as 1/b, as dimensional analysis requires.

**B3, B4.** In #269's static mass mode `r = √L R`, `F = A + 1/L − 1` the horizon moves to ρ_h = δL/(2κ) and the surface gravity to κ_L = κ − δL/(2κ). In its own chart `V = e^{κ_L v}`, U measured from ρ_h, the
slope of `∂_Vδφ` against ln V is 6.3e-8 and 6.8e-7 of (π/√2)/2 at V = 1e-8 (U = −0.5, −1; falling like V). In the chart `V = e^{κv}` of the background the same mode has `∂_Vδφ ⊃ (√2U_K/(2κ²)) ln V`, `U_K = −Π(ρ_b)/V_b`, per unit
δL (the gauge ξ^V = (δκ/κ)V ln V): measured −1.33537 and −4.3139115, predicted −1.335368 and −4.3138906. That term is proportional to U and vanishes on U = 0, as the U-independent coefficient of Route A requires.

## 5. Route B: the exact family in a chart built from the rays (B1, B5, B2, X1)

#269's family (b = 1, c = −1, K = κ = 2/π): `ds² = −F dv² + 2dv dρ + r²dΩ²`, `r = λR`, `λ² = 1 + βe^{κv}`, `F = A − 2ρλ′/λ + 1/λ² − 1`. B1: it solves vv, vρ, ρρ, θθ and the scalar equation of my Ricci tensor in this chart to 3.7e-31
(2 values of β, 3 points; this repeats #269's check with my own tensor). The chart: V = e^{κv}; U = −ρ_b/V_b, the label of the outgoing ray `dρ/dv = F/2` by its ρ at the base ray V_b (zero gauge vector there; smooth in the Kruskal U).
`δφ(V, U) = ∂_βφ` at fixed (V, U) by a central difference in β (step 1e-9, and 1e-10 in X1); `∂_Vφ = Φ′F/(2κV)` along the ray, so no derivative in V is taken numerically. B5: on U = 0, δr = c₂V with c₂ = 0.5000000000 and c₁ = 0.

Slope of `∂_Vδφ` against ln V over `−(π/√2)c₂` between decades (the label is the last decade's V):

| V_b, U | 1e-4 | 1e-5 | 1e-6 | 1e-7 | 1e-8 |
|---|---|---|---|---|---|
| 1, −0.5 | 1.00881026 | 1.00109839 | 1.00013149 | 1.00001531 | 1.00000175 |
| 1, −1 | 1.04244313 | 1.00495788 | 1.00056586 | 1.00006358 | 1.00000706 |
| 0.05, −0.5 | 1.00288906 | 1.00043333 | 1.00005775 | 1.00000722 | 1.00000087 |
| 0.05, −1 | 1.00596401 | 1.00089155 | 1.00011859 | 1.0000148 | 1.00000177 |

The deviation falls like V ln V; it does not depend on U in the limit, nor on the base ray (X1: 1.00000175 for β-step 1e-9 and 1e-10). The mechanism is a resonance: on the horizon the source of the ray equation is `−(β/2)e^{κv}` against the
homogeneous solution `e^{κv}`, so `ρ = V(C − (β/2κ) ln V)`; `δφ = √2δρ`, `c₂ = β/2` give `−(π/√2)c₂`. A rate K ≠ κ removes it (`rate` fires B2).

## 6. What is shared with #302 and what is not

| | shared | not shared |
|---|---|---|
| background | A, R, θ = arctan(ρ/b) (D1 solves my equations with them), Kruskal V = e^{κv}, U = −e^{−κu}, Π = −VU, Ω₀ = A/(2κ²Π) (D2) | the series for ρ(Π), Ω₀, 𝒱 (my own recursion) |
| equations | the field equations R_ab = ε φ_aφ_b + 𝒱 g_ab, ε□φ = 𝒱′ | the Ricci tensor, the linearization, the Lie derivative (D3), the potential (from E_θθ) |
| the mode | Y = δr − (dr₀/dθ₀)δθ (D4: invariant), the target Y = c₂VR (S3: the log-extended system is consistent with it; with Y = V it is not, `rhs`) | the slice constraint M, #258's slice code, #248's solver and functionals, #302's gauge-vector quadrature |
| Route B | #269's family (B1 checks it) | the ray integration, the chart, the label, the static planted case with its own κ_L |

## 7. Checks and mutations

| id | what | result |
|---|---|---|
| D1 | static chart, own Ricci, potential from E_θθ | residuals 0 |
| D2 | background series solves the 5 equations | no residual at degree ≤ 8 |
| D3 | pure gauge solves the linearized equations | no residual at degree ≤ 8 |
| D4 | Y vanishes on pure gauge | zero to degree 10 |
| S1, S1g | relations y₂₁, y₁₂, horizon laws, Y_N″ condition; y₂₁ = 0 excluded | in the span; excluded |
| S1b | smooth solution with c₁ = c₂′ = 0 exists | consistent |
| S2 | no smooth Y = VR, UR; Y = R yes | inconsistent, inconsistent, inconsistent (degree ≤ 3); consistent |
| S3 | log coefficient | −1/κ, no spread, 2 gauges |
| S4 | planted Y = R | q = 0 |
| S5 | A″(0) control, 6 members | equal |
| S6 | q = 2t − F_Π; smooth iff t = F_Π/2 | equal |
| B1 | family against my Ricci | 3.7e-31 |
| B5 | c₁, c₂ | 0, 0.5 |
| B2 | slope ratio, 4 (V_b, U) | 1.0000018 … 1.0000071 at 1e-8 |
| X1 | two routes; step independence | 1.00000175, 1.00000175 |
| B3 | planted static mode, own chart | 6.3e-8, 6.8e-7 |
| B4 | static mode in the background's chart | 5e-6 |

All 18 pass in the runs of record (`derive` 4 of 4, 7 s, 84 MB; `series` 8 of 8, 80 s, 88 MB; `check` 6 of 6, 211 s, 77 MB; code b1ef08e). Eighteen mutations, each with the checks it must fire declared before the runs of record
(raw README; computed in development runs `mut_d1` to `mut_d6`):

| mutation | command | fires | effect in the run of record |
|---|---|---|---|
| eps | derive | D1, D2, D3 | canonical scalar: the static member no longer solves E_tt, E_ρρ |
| a2 | derive | D2, D3 | A″(0) × 1.01 in the background series: residual at low degree |
| om | derive | D2, D3 | Ω₀ = A/(κ²Π): residual at low degree |
| gauge | derive | D3 | sign of Ω₀∂ξ in the Lie derivative |
| gth | derive | D4 | g × 1.01: Y ≠ 0 on pure gauge |
| fpi | series | S1, S1g, S6 | F_Π doubled: the forms leave the span; S6's t shifts |
| nolog | series | S3, S4, S5, S6 | the log ansatz dropped: S3, S5, S6 inconsistent, S4 has no log column |
| pred | series | S5 | A″(0) × 1.01 in the prediction |
| rhs | series | S3, S5, S6 | Y = V (R → 1) is not a solution of the log system |
| hor | series | S1b | y₂₀ = 1 imposed against the horizon law: inconsistent |
| s2t | series | S2 | the V²U term smooth solutions need makes the degree ≤ 3 case consistent |
| rate | check | B1, B2, B5, X1 | K = 1.01κ: not a solution; the resonance is gone (B2 ratios 0.9985 at V = 1e-8); c₂ = 0.4256 |
| h | check | B1, B2, X1 | 1/λ² − 1 dropped |
| beta | check | B1, B2, X1 | sign of β in F: ratio −1 |
| phi | check | B1, B2, B4, X1 | dθ/dρ × 1.01 |
| chart | check | B3 | the static mode in the background's chart: slope / scale −1.2 and −3.9 (bar 1e-6) |
| pgauge | check | B4 | prediction with the sign of δκ/κ flipped |
| step | check | B5, X1 | β-step 0.3: c₂ = 0.50000007 (bar 1e-8), X1 differs between steps |

`log_at_S.py mutations`: 18 of 18 fired sets equal the declared ones (`log_at_S_output.txt`).

Blind spots, kept: S1b and S2 are fired only by the contrived `hor` and `s2t`; `a2`, `om` corrupt the background and were run under `derive` only; B1 is not fired by `chart`, `pgauge`, `step`; X1 fires with B2 and alone only through the step.

## 8. What this does not show

- Anything beyond formal series. The relations are statements about Taylor series of C^∞ perturbations at S (and the log solutions are formal); the series are truncated at N = 5 to 7 (in development q = −1/κ also at N = 6, 7, 8, probe `p9`, `dev/probes`), so a relation that appears only at higher degree was not looked for.
- Gauges in which the fields are not smooth in the Kruskal chart. S1g covers every smooth gauge (all five fields); a gauge with continuous but not smooth fields, or with a log in the gauge vector (which puts 1/V into δg_VV), was not examined.
- That Y = c₂VR is a solution in a stated function space. It is a solution of the log-extended formal system and the linearization of the exact family; its regularity class at S (δφ ∈ C⁰, ∂δφ ∼ ln V) is read off the series, not proved.
- Other members. κ and b are free in the series but the horizon is at the throat, A″(0) = 2c; at c ≠ −1 the horizon is not at the area minimum (#301). Route B is b = 1, κ = 2/π, first order in β (central difference), and not the nonlinear solution at finite β.
- The decaying partner. Y = c₂′UR has no smooth solution (S2) and by the V ↔ U symmetry of the equations a ln U with the same coefficient; that was not run separately.
- Decay of Y_rem (#322) and the time-dependent smooth data: nothing here evolves anything.

## 9. Text of the repository

The result agrees with #302; no sentence of Paper IV or of the notes changes. Two refinements for the coordinator to decide on, not applied:

| file | now | could say |
|---|---|---|
| `horizon_basis.md` Sec. 3, "What survives" | "carries the `Y_N` condition above (D6's free coefficients y₁₂ − y₂₁)" | the condition holds per branch: y₂₁ = (F_Π/2)c₂ and y₁₂ = (F_Π/2)c₂′ (#321, S1); both coefficients are fixed, not a combination |
| `horizon_basis.md` Sec. 8, task 2 | "which combination of y₁₂, y₂₁ it fixes in terms of c₂ + c₂′" | *Done in #321*: both, separately (y₂₁ ↔ c₂, y₁₂ ↔ c₂′) |
| Paper IV, "Regular data, singular mode" | "every perturbation smooth at S satisfies ∂_s²Y_N(0) = −F_Π(0)Y_N(0)" | true as written; it is the sum of the two per-branch relations, and the statement holds in every smooth gauge (S1g) |

## 10. Recommended tasks

1. **The finite-β family at S.** *What:* the nonlinear family (β = ±0.3) in the chart of Route B, ρ(V, U; β) near V = 0 from the same ray ODE. *Why:* #302 Sec. 7 left whether the nonlinear solution is smooth at S open; the ODE has the same resonance, so ρ ∋ −(β/2κ)V ln V at every β. *Acceptance:* the V ln V coefficient of ρ against `−β/(2κ)` to 1e-6 at β = ±0.3, ±0.05; a reparametrization V → V(1 + aβ ln V) shifts only U-proportional terms; a mutation K → 1.01κ removes it. *Depends on:* #269, #321.
2. **The mode's finite-energy regularity.** *What:* whether ∫(∂_Nδφ)² ds converges at S for the log mode (∂_Nδφ ∼ ln s is square integrable), and the energy flux of c₂VR through V = ε. *Why:* the mode is not smooth but may still be a finite-energy perturbation, which decides how Paper IV words it. *Acceptance:* the closed form of the integral against the series solution to 1e-6; a mutation scaling q by 1.01 fails. *Depends on:* #321.

(#322 is the remainder of smooth data and is not repeated.)

## Reproduce

From this directory, `~/Source/negative-mass/.venv/bin/python`, one thread:

```
log_at_S.py derive [--mutate=eps|a2|om|gauge|gth]                                 # D1–D4, 7 s
log_at_S.py series [--only=S1,..] [--mutate=fpi|nolog|pred|rhs|hor|s2t]           # S1–S6, 80 s
log_at_S.py check  [--only=B1,B2,..] [--mutate=rate|h|beta|phi|chart|pgauge|step] # B1–B5, X1, 211 s
log_at_S.py mutations DIR PREFIX                                                  # fired sets of the stored logs DIR/PREFIX<name>.log against the declared ones
```

`log_at_S_output.txt` is `derive`, `series`, `check` and `mutations` of the runs of record (`I321-mac/runs`, code commit b1ef08e), in that order.

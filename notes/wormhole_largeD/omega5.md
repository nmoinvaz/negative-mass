# #263: Ω₅ from the order-7 matching integral

Status: Ω₅ is **computed** from G₇ = 0 by method B on two grids at two precisions, agreeing to 382.7 digits, and is **identified** by PSLQ at 330 digits in the weight ≤ 6 basis Q6 (and in Q5, Q7, Q8 and, with weight-7 constants, Q9) as

    Ω₅ = (256/45)ln⁶2 + (64/9)π²ln⁴2 − (2560/3)ζ(3)ln³2 + (176/45)π⁴ln²2 − 3072ζ(3)ln²2 + 640π²ζ(3)ln 2 + (9728/45)π⁴ln 2 − 12288ζ(3)ln 2 + 41728ζ(5)ln 2
         − (38812/315)π⁶ + 3200ζ(3)² − (256/3)π²ζ(3) + 74240ζ(5) + (39424/45)π⁴ − 36864ζ(3)
       = 40860.455354525687312544064826033704278039372215373284079454102667838994220332959565447528072780104912968505385…,

and the same integral is reduced **by hand** to that closed form: each step is checked symbolically, or numerically at 30 to 100 digits (checks RD1 to RD8), and the closed form comes out of `reduce7` and is not typed in. The reduction uses E1, E2, E3 (#26, #208, re-derived in #209), E4
(#220) at one level higher, and two new steps, E2W (an exchange identity for L₄ with a weight) and E5 (for Y₄). **The closed form was in hand, and equal to a 148-digit development grid to 148.1 digits, before the real PSLQ ran** (raw README, amendments A1 and A2): PSLQ confirms integer coefficients
that were already known, to 330 digits; it was not a blind search. Every number is in `order7_output.txt` (the runs of record, raw folder `negative-mass-raw/I263-rpi5`, README first) unless it names another raw file.
No sentence here claims novelty; the prior equations used are README.md's order equations (matched.py, order3.py … order6.py), the closed forms and exchange identities of #26, #208, #209 and #220, Clausen's functions, the Beta-function moments of ln cos, Parseval's identity and classical Euler sums.

## Result

| evaluation | method | N, w, V, working digits | Ω₅ against the closed form |
|---|---|---|---|
| G1 | the hierarchy from its definitions, 192 panels | 285, 1/4, 48, 445 | 418.9 digits |
| G2 | same, 157 panels | 275, 3/10, 47, 410 | 382.7 digits |
| B3, B2, B1 (reported) | same | 170, 1/4, 38, 285; 100, 1/4, 30, 170; 60, 2/5, 24, 100 | 251.0, 148.1, 77.5 digits |
| A7 | closed forms of Y₀ … Y₃, L₀ … L₃, E1 … E5, quadrature at 40 digits | 40 + 150 guard | 39.1 digits |
| X3 | #195's 64 values of W² (qd, 60 digits), polynomial extrapolation in 1/s with Ω₀ … Ω₄ at their closed forms | windows K = 14 … 44 | 7.5e-28 … 3.3e-25 in absolute difference |

- **Ω₅** is the value above (G1; G2 agrees to 382.7 digits, |difference| 8.99e-379). The difference between evaluations is a spread over grid and precision, an estimator range, **not a rigorous bound**; nothing here is a rigorous enclosure. G₇(Ω₅ = 0) = 4πΩ₅ = 513467.625456446519397870844712… (RD6).
- **In 1/(D − 3).** Ω₅/32 = 1276.8892298289277285170020258135532586887…, with the closed form (`order7_output.txt`, run_07)
  Ω₅/32 = (8/45)ln⁶2 + (2/9)π²ln⁴2 − (80/3)ζ(3)ln³2 + (11/90)π⁴ln²2 − 96ζ(3)ln²2 + 20π²ζ(3)ln 2 + (304/45)π⁴ln 2 − 384ζ(3)ln 2 + 1304ζ(5)ln 2 − (9703/2520)π⁶ + 100ζ(3)² − (8/3)π²ζ(3) + 2320ζ(5) + (1232/45)π⁴ − 1152ζ(3). The series is

  ```math
  W^2a^2 = 2(D-3) - 4\ln 2 + \frac{\Omega_1/2}{D-3} + \frac{\Omega_2/4}{(D-3)^2} + \frac{\Omega_3/8}{(D-3)^3} + \frac{\Omega_4/16}{(D-3)^4} + \frac{\Omega_5/32}{(D-3)^5} + O((D-3)^{-6}),
  ```
  with Ω₁/2 = 5.2116801893692586, Ω₂/4 = 18.592346681983329, Ω₃/8 = −81.605013446257016, Ω₄/16 = −799.65624125144966, Ω₅/32 = 1276.8892298289277.
- **The planted orders by the same code.** `omega_chain` finds Ω₀ … Ω₅ in turn from G₂ … G₇ = 0 with no closed form of Y_k or L_k and no input value; Ω₀ … Ω₄ against their closed forms (PL0 … PL4): G1 419.3 to 420.3 digits, G2 383.0 to 384.5 (the issue's bar: 100).
  The slope of G_k in Ω_{k−2} is −4π to 420.5 (G1) and 384.1 (G2) digits at every k = 2 … 7 (SL), and |G_k| at the chain's root is below 1.1e-436 and 4.9e-402 times |slope| (RES).
- **Rounding and grids.** The grids follow digits ≈ N log₁₀(7.7/w) of #208 at every order including the seventh: N = 60, 100, 170, 275, 285 give 77.5, 148.1, 251.0, 382.7, 418.9 digits of Ω₅ against the closed form. The cached chain has 48 integration units
  for k = 2 … 7 against 120 uncached; **the working digits need only the expected digits plus about 20**, where #220 used 1.6 times as many (`dev_02` … `dev_09`).
  Run times (rpi5, three processes at nice 10): G1 1 h 34 m 46 s (panel setup 8 m 22 s), G2 1 h 4 m 52 s (6 m 52 s), B1 … B3 19 m 10 s; peak resident sets 1.02 GB, 0.79 GB, 0.44 GB. The issue estimated about 2 h each at N ≈ 200.
- **#195's W²** (the issue's independent reference): `w2_output.txt` has no fit of Ω₅. A refit of its 64 values (qd, 60 digits, D = 50 … 10⁴) with Ω₀ … Ω₄ fixed at the closed forms, by polynomial extrapolation in 1/s through the K largest-s points, K = 14, 18, 22, 26, 30, 36, 44,
  gives Ω₅ to 7.5e-28, 5.8e-27, 2.3e-26, 5.9e-26, 1.1e-25, 2.1e-25, 3.3e-25 (X3, bar 1e-22, set knowing a development extrapolation; an estimator range). Two methods that share nothing agree to 7.5e-28 or better, 31 digits relative.
- **Remainder** (rem = W²a² − [s − 4 ln 2 + Ω₁/s + … + Ω₄/s⁴], s = 2(D − 3), against #195's qd W², this Ω₅ and the six-term series of #220):

  | D | s | rem·s⁵ | (rem·s⁵ − Ω₅)·s |
  |---|---|---|---|
  | 1000 | 1994 | 42708.395 | 3684791.8 |
  | 2000 | 3994 | 41781.483 | 3678586.1 |
  | 5103 | 10200 | 41220.725 | 3674753.1 |
  | 10000 | 19994 | 41044.187 | 3673534.1 |

  rem·s⁵ approaches Ω₅ (the distance falls like 1/s, and the last column is rem₇·s⁶ with this Ω₅, which approaches Ω₆). The two-point extrapolation from D = 5103 and 10⁴ gives Ω₅ ≈ 40860.33 (R6: relative difference −3.05e-6 from the computed value) and **Ω₆ ≈ 3.672·10⁶** (reported, an estimate, not a bound):
  the seven-term series leaves Ω₆/s⁶. R6 and R7 see an error of Ω₅ of about 4 or 0.13 or more, not less; X3 is the sharp check.

- **PSLQ.** `mp.pslq` at d = 330 digits (tolerance 1e-318, maxsteps 10⁶, Ω₅ of G2; the two grids share 382.7 digits, so every relation is also checked at 350 digits), bases declared before the run (raw README), each with a planted relation recovered exactly and a perturbed basis
  (a constant + √2·10⁻³) that finds nothing. The weight-7 constants of Q9 and the polylogarithm products of Q2, Q7, Q8 are in the table's last column:

  | basis | constants besides Ω₅ | n | H | budget | result |
  |---|---|---|---|---|---|
  | Q1 pure weight 6 | 8 monomials | 8 | 10¹² | 111.7 | not found |
  | Q2 pure weight 6 with Li₄(½)ln²2, Li₄(½)π², Li₅(½)ln 2, Li₆(½) | Q1 and 4 | 12 | 10¹² | 160.9 | not found |
  | Q3 weights 6, 5 | Q1 and 6 | 14 | 10¹¹ | 170.5 | not found |
  | Q4 weights 6 … 4 | Q3 and 4 | 18 | 10¹⁰ | 196.7 | not found |
  | Q5 weights 6 … 3 | Q4 and 3 | 21 | 10¹⁰ | 227.6 | **found** |
  | Q6 weight ≤ 6 | Q5 and ln²2, π², ln 2, 1 | 25 | 10⁹ | 242.8 | **found**, the same relation |
  | Q7 weight ≤ 6 with Li₄(½){1, ln 2, ln²2, π²}, Li₅(½){1, ln 2}, Li₆(½) | Q6 and 7 | 32 | 10⁹ | 307.9 | **found**, the same relation, zero on the 7 Li products |
  | Q8 weights 6 … 3 with the same 7 | Q5 and 7 | 28 | 10⁹ | 270.7 | **found**, the same relation, zero on the 7 |
  | Q9 weight ≤ 7 | Q6 and the 11 monomials of weight 7 (ln⁷2, …, π²ζ(5), ζ(7)) | 36 | 10⁸ | 308.1 | **found**, the same relation, zero on the 11 |

  The relation (Q6's vector, residual below 1e-318, holds at 350 digits on G1 and G2 to 3.35e-344, bar 1e-338; the perturbed basis finds nothing) is
  315 Ω₅ = 1792 ln⁶2 + 2240π²ln⁴2 + 1232π⁴ln²2 − 38812π⁶ − 268800ζ(3)ln³2 + 201600π²ζ(3)ln 2 + 1008000ζ(3)² + 13144320ζ(5)ln 2 + 68096π⁴ln 2 − 967680ζ(3)ln²2 − 26880π²ζ(3) + 23385600ζ(5) + 275968π⁴ − 3870720ζ(3)ln 2 − 11612160ζ(3),
  with zero coefficients on ln⁵2, π²ln³2, ln⁴2, π²ln²2, ln³2, π²ln 2, every constant of weight ≤ 2, **every Li₄(½), Li₅(½), Li₆(½) product** and **every weight-7 constant, π²ζ(5) and ζ(7) among them**. It is 315 times the closed form of the reduction (RD8: the integer vectors are equal). Q1 … Q4 find nothing because the lower weights are there.
  **The task text names π²ζ(5) among the weight ≤ 6 constants; it has weight 7 (π²: 2, ζ(5): 5), so it is in Q9 and not in Q1 … Q8.** I believe the alternating Euler sums of weight 6 have one irreducible beyond Li₆(½) (*unverified*: dimensions 1, 1, 2, 3, 5, 8, 13 for weights 0 … 6, so 33 up to weight 6 against the 32 constants of Q7); it is in no basis, and a "not found" is for these bases only. The relation found does not need it.
- **Weights.** Ω₅ is not of uniform weight (as Ω₃ and Ω₄ were not): **the pattern of Ω₄ continues, one weight up**. Its weight-6 part is −78718.5817…, weight-5 part 88791.0309…, weight-4 part 75100.6318… and weight-3 part −44312.6257… (run_09; the parts are larger than Ω₅ = 40860.46 and cancel).
  All 8 monomials of weight 6 are present; of weight 5 four of six (ln⁵2 and π²ln³2 are absent); of weight 4 two of four (π⁴ and ζ(3)ln 2 only, as in Ω₄); of weight 3 only ζ(3) (as in Ω₃ and Ω₄); nothing of weight ≤ 2. No Li₄(½), Li₅(½) or Li₆(½) enters Ω₁ … Ω₅.
  The pure ζ(3) terms of Ω₃, Ω₄, Ω₅ are −384, 3072, −36864 (the same in F below, f₄, f₅, f₆), and 96(−4)^{k−3}(k − 3)! gives them for k = 4, 5, 6: three numbers fitted by a one-parameter formula, an observation and not a claim about the series.

## The order equations at order 7

README.md's equations (check E0: the ε-expansion of the master equation, sympy, residual 0, now through k = 7):

    dY_k/dl = F_{k−1} = Σ_{j<k} L_{k−1−j} l^j/j!,    (R − 1)L_k′ − L_k/2 = h_k = (1 − 1/R)L_{k−1} + Σ_{j<k} Y_j c_{k−j} + Y_k/2,
    L_k = −2Ω_{k−1} + (R − 1)^{1/2}C_k,    G_k = C_k(∞) = ∫₁^∞ (h_k − Ω_{k−1})(R − 1)^{−3/2}dR = 0,

with c₁ … c₆ as in `order6.py` and the new

    c₇ = Ω₀l⁶/720 + Ω₁l⁵/120 + Ω₂l⁴/24 + Ω₃l³/6 + Ω₄l²/2 + Ω₅l + Ω₆ + l⁷/3360 − l⁶/720 + l⁶/(720R),

the coefficient of ε⁷ of (1 + εΩ)e^{εl} − (½ + ε)(1 − 1/R)e^{−εl} − e^{−εl}/(2R). G_k is affine in Ω_{k−2} with slope −4π; Ω₆ drops out of G₇ (it is h₇(1), subtracted), so **G₇ = 0 gives Ω₅**.

## The reduction by hand

Notation as omega4.md: θ ∈ [0, π/2), R = sec²θ, S = sin²θ, u = tan θ, ℓ = ln cos θ = −l/2, C = Cl₂(π − 2θ), D₃ = Cl₃(π − 2θ) + ¾ζ(3), Lc2 = ∫₀^θℓ², Ω₀ = −4 ln 2.
The derivatives close: ℓ′ = −u, u′ = 1 + u², C′ = 2ℓ + 2 ln 2, D₃′ = 2C, Lc2′ = ℓ².

1. **E1, E2, E3 at k = 7, m = 6** (#209; they need G_j = 0 for j ≤ 6):
   G₇ = 8∫H₆ℓ/S + 2∫(F₆ − L₆) + 2∫(S′₇ − Ω₆)/S − 4∫(ℓ − Ω₀)F₅, H₆ = S L₅ + Σ_{j≤5}Y_jc_{6−j}, S′₇ = Σ_{j≤5}Y_jc_{7−j}.
   Y₆ and L₆ are gone. What is left is Y₀ … Y₅ and L₀ … L₅; Y₄, Y₅, L₄, L₅ have no closed form.
2. **L₅.** Its weight is 8ℓ (from H₆) + 2l (from F₆ − L₆) − 4(ℓ − Ω₀) (from F₅) = 4Ω₀, and E2 at m = 5 gives 4Ω₀∫L₅ = 8Ω₀∫h₅ℓ/S, h₅ = S L₄ + Σ_{j≤4}Y_jc_{5−j} + ½Y₅.
3. **E4 at one level up: Y₅.** Its terms are 8ℓc₁/S, 2c₂/S and 4Ω₀ℓ/S, so its weight is w₅ = (2Ω₁ + 8Ω₀ℓ − 22ℓ²)/S − 12ℓ: the weight of Y₄ at k = 6 (the pair c₁, c₂ and the collapse of the top L repeat at every k).
   With W₄ = cot θ (22ℓ² − 8Ω₀ℓ − 2Ω₁) + 16C, **∫Y₅w₅ = −2∫(22ℓ² − 8Ω₀ℓ − 2Ω₁ + 16uC)F₄**, F₄ = L₄ + L₃l + ½L₂l² + ⅙L₁l³ + (1/24)L₀l⁴ (check E4_7).
4. **L₄: E2 with a weight (E2W).** L₄ now has the weight g₄ = 12ℓ² (from F₆ − L₆, F₅ and h₅) − 2(22ℓ² − 8Ω₀ℓ − 2Ω₁ + 16uC) = −32ℓ² + 16Ω₀ℓ + 4Ω₁ − 32uC. From (L cot θ)′ = 2h/S (the equation for L_m in θ) and
   L cot θ → 0 at π/2 (G_m = 0): ∫L_m g dθ = −2∫P h_m/S with P′ = g tan θ, P(0) = 0 (E2 is g = 1). Here P = 32ℓ³/3 − 8Ω₀ℓ² − 4Ω₁ℓ − 32(uC + ℓ² + 2 ln 2 ℓ − D₃/2), because ∫u²C = uC + ℓ² + 2 ln 2 ℓ − D₃/2
   (checks E2Wa for the weight 12ℓ², E2Wb for g₄). L₄ is gone; ½Y₄ stays in h₄.
5. **E5: Y₄.** Collect the Y₄ terms: 8ℓc₂/S, 2c₃/S, 8Ω₀ℓc₁/S (from step 2) and −P/S (from step 4): w₄′ = [2Ω₂ + (8Ω₁ + 8Ω₀² + 64 ln 2)ℓ + (32 − 28Ω₀)ℓ² − (20/3)ℓ³ + 32uC − 16D₃]/S − 20ℓ² − 8Ω₀ℓ.
   Its primitive is in the generators of order 6: with Lc1 = (C − 2 ln 2 θ)/2, ∫ℓ^b/S = −ℓ^b cot θ − b Lc_{b−1} (Lc₀ = θ), ∫C tan θ = −Cℓ + 2Lc2 + 2 ln 2 Lc1, and **Z = ∫C cot θ cancels** between 32∫uC/S = 32(Z + ∫C tan θ) and −16∫D₃/S = 16D₃cot θ − 32Z:
   W = −2Ω₂cot θ − a₁(ℓ cot θ + θ) − a₂(ℓ²cot θ + 2Lc1) − a₃(ℓ³cot θ + 3Lc2) + 32(−Cℓ + 2Lc2 + 2 ln 2 Lc1) + 16D₃cot θ − 20Lc2 − 8Ω₀Lc1 + const,
   a₁ = 8Ω₁ + 8Ω₀² + 64 ln 2, a₂ = 32 − 28Ω₀, a₃ = −20/3, and the constant fixed by W(π/2) = 0 (W′ = w₄′: sympy residual 0, check XS). Since dY₄/dθ = 2uF₃,
   **∫Y₄w₄′ = −2∫W u F₃**, F₃ = L₃ + L₂l + ½L₁l² + ⅙L₀l³ (check E5_7). Boundary terms: [WY₄] is O(θ) at 0 (W ~ cot θ, Y₄ = O(θ²)) and O(x lⁿ) at π/2 (x = π/2 − θ).
6. **So G₇ = ∫₀^{π/2}Ĩ₇dθ with Ĩ₇ built from the closed forms of Y₀ … Y₃, L₀ … L₃** (#26, #208; L₃ through Lc2) and the weights above (check G7R on the hierarchy, RD1 and A7 on the closed forms).
   Its decomposition in C^a D₃^b Lc2^c (sympy, RD1) has the monomials 1, C, D₃, Lc2, C², CD₃, CLc2 (as at order 6) and **three new ones, D₃², D₃Lc2, Lc2²**.
7. **Finite parts and the generators.** The masters are free of ln 2 in lt = ln(2 cos θ) = ℓ + ln 2, C′ = 2lt, and K₂ = ∫₀^θlt², K₃ = ∫₀^θlt³ (K₂ = Lc2 + ln 2 C − ln²2 θ). The integrand is converted with ℓ = lt − ln 2 and Lc2 = K₂ − ln 2 C + ln²2 θ, and the symbolic Ω₁ … Ω₄ and J₂ are inserted
   (RD1: the converted integrand against the numeric integrand of closed forms, 45 digits). Several pieces diverge at π/2 and cancel only in the sum, so ⟨f⟩ is the finite part (the constant term of x⁰ln⁰x, x = π/2 − θ) of ∫f, as in omega4.md step 5.
   For a primitive F, ⟨F′⟩ = B[F] := FP F(π/2) − FP F(0), with the series of the generators at π/2 and 0 to order 10 (RD4 checks them at x = θ = 10⁻³ to 8.9e-36 and 8.9e-33).
   One step is I ← I − dF with const += B[F]; the integral of the final I is a sum of masters. The identity I_start − Σ dF − I_final = 0 is checked symbolically (RD2).
8. **The 37 integrations by parts** (the names of RD2), in the order they are applied:
   - **LL**: F = a·uK₂², a the coefficient of u²K₂²: ⟨u²K₂²⟩ = B[uK₂²] − 2⟨u lt²K₂⟩ − ⟨K₂²⟩.
   - **DL**: F = −a·lt D₃K₂ for a·uD₃K₂.
   - **DD**: the coefficient of D₃² is −512[(3lt + ln 2)(1 + u⁻²) + 3], with primitive A = 512(3lt + ln 2)cot θ (the θ terms cancel); ⟨cD₃²⟩ = B[AD₃²] − 4⟨A D₃C⟩. **The 1/u part of the coefficient of CD₃ cancels exactly** (RD3), as the cot θ parts did at order 6.
   - **CD, CL, CC**: F = −G·CD₃ (G = ∫g dlt), F = u g̃ CK₂ and F = u f C², with g̃ − g̃′ = g and f − f′ = p the polynomials of the u² parts of the coefficients (as M1, M3 of omega4.md).
   - **D_u, D**: the u θ D₃ part, then the primitive of the rest of the D₃ coefficient: ⟨c D₃⟩ = B[AD₃] − 2⟨AC⟩, A in θ, C/2, K₂, K₃ and lt^b cot θ; K₃ is new at this order (from ∫lt⁴/S). **The 1/u part of the C coefficient cancels against the cot θ part of −2AC** (RD3).
   - **K3C**: ⟨CK₃⟩ = B[D₃K₃/2] − ½⟨lt³D₃⟩.
   - **C_u_lt_b, C_u2th_lt_b**: ⟨u lt^bC⟩ = B[−lt^{b+1}C/(b + 1)] + 2⟨lt^{b+2}⟩/(b + 1); ⟨u²θ lt^bC⟩ = B[uθ lt^bC] − ⟨θ lt^bC⟩ − ⟨u lt^bC⟩ + b⟨u²θ lt^{b−1}C⟩ − ⟨uθ lt^b·2lt⟩.
   - **K_u_lt_b, K_u2th**: the same for K₂. **P_u2, P_ui_lt_b, P_u_th_j_lt_b**: the pure terms: ⟨u²θ²⟩, ⟨lt^b cot²θ⟩ = ⟨lt^b/S⟩ − J̃_b with ⟨lt^b/S⟩ = −b J̃_{b−1} + B[−lt^b cot θ], and ⟨uθ^j lt^b⟩.
   The 1/u parts and the u parts cancel or reduce, and **no u, 1/u, or new generator is left**: I_final is a polynomial in lt, θ, C, D₃, K₂ with 18 monomials.
9. **The 18 master integrals** (7 of them the moments J̃₀ … J̃₆; RD5: each against quadrature at 30 digits, 25.6 to 30.0 digits). J̃_n = ∫₀^{π/2}lnⁿ(2cos θ) dθ = n![aⁿ](π/2)Γ(a + 1)/Γ(a/2 + 1)² (Beta function, free of ln 2): J̃₀ = π/2, J̃₁ = 0, J̃₂ = π³/24, J̃₃ = −3πζ(3)/4, J̃₄ = 19π⁵/480, J̃₅ = −(5π/8)(18ζ(5) + π²ζ(3)),
   J̃₆ = (5π/2688)(6048ζ(3)² + 55π⁶): ζ(3)² enters first here. The others, with sums over n of the Fourier series lt = Σ(−1)^{k+1}cos 2kθ/k, C = Σ(−1)^{n+1}sin 2nθ/n², D₃ = ¾ζ(3) + Σ(−1)ⁿcos 2nθ/n³:
   - ⟨θ²⟩ = π³/24; ⟨C²⟩ = π⁵/360 (Parseval, #209); ⟨lt C²⟩ = 0 and ⟨θ lt C⟩ = −⟨C²⟩/4 (by (C²)′ = 4 lt C and C(π/2) = 0); ⟨θC⟩ = πζ(3)/4.
   - ⟨lt²D₃⟩ = (π/4)(5ζ(5) − 2ζ(2)ζ(3)) + π³ζ(3)/32, from #220's ∫ℓ²D₃ with lt = ℓ + ln 2 (the ln 2 terms cancel).
   - ⟨CK₂⟩ = B[D₃K₂/2] − ½⟨lt²D₃⟩ = (7/8)ζ(3)J̃₂ − ½⟨lt²D₃⟩; ⟨lt CK₂⟩ = ¼B[C²K₂] − ¼⟨lt²C²⟩ = −⟨lt²C²⟩/4.
   - ⟨θ²lt²⟩ = 11π⁵/1440 (**new**): lt² = π²/12 + Σα_n cos 2nθ, α_n = ((−1)ⁿ/n)(2H_n − 1/n); ∫θ²cos 2nθ = π(−1)ⁿ/(4n²); so ⟨θ²lt²⟩ = π⁵/288 + (π/4)Σ(2H_n − 1/n)/n³, with Euler's ΣH_n/n³ = π⁴/72. ⟨θK₂⟩ = (π²/8)J̃₂ − ½⟨θ²lt²⟩.
   - ⟨K₂²⟩ = π⁷/3456 − π⁷/2880 + (π/16)((61/6)ζ(6) − 6ζ(3)²) (**new**): K₂ = (π²/12)θ + Σα_n sin 2nθ/(2n) and Parseval, with ΣH_n/n⁵ = (7/4)ζ(6) − ½ζ(3)² and ΣH_n²/n⁴ = (97/24)ζ(6) − 2ζ(3)².
   - ⟨lt³D₃⟩ = π((9/16)ζ(3)² − (11/4)ζ(6)) (**new**): from ∫₀^{π/2}(2cos θ)^a cos 2nθ dθ = (π/2)Γ(a + 1)/(Γ(a/2 + n + 1)Γ(a/2 − n + 1)): with b = a/2, ⟨(2cos θ)^aD₃⟩ = (π/2)G(a)[¾ζ(3) − b Σ_n (1 − b)_{n−1}/((1 + b)_n n³)], G = Γ(a + 1)/Γ(a/2 + 1)²,
     and the sum is ζ(4) + … + b²(ζ(6) − 2ΣH_n/n⁵ + 2ΣH_n²/n⁴) = … + b²((67/12)ζ(6) − 3ζ(3)²); the coefficient of a³ gives the value. It uses the same two Euler sums as ⟨K₂²⟩.
   - **⟨lt²C²⟩ = π((3/4)ζ(3)² − (95/96)ζ(6)) is identified, not derived.** PSLQ at 55 digits in {π⁶, ζ(3)²} for ⟨lt²C²⟩/π gives the relation (−18144, −19, 13608), and the value holds at 100 digits (RD5b). By Parseval ⟨lt²C²⟩ = (π/4)Σ_m s_m², where lt·C = Σs_m sin 2mθ,
     s_m = (−1)^m[(H_m^{(2)} − ζ(2))/m + 2H_m/m² − 3/(2m³)] (derived by hand, not checked here): a sum of Euler sums of weight 6 that I did not evaluate.
10. **Assembly.** ⟨F′⟩ sums: G₇ = Σ_i B[F_i] + Σ (coefficient)·(master). Its slope in Ω₅ is −4π exactly (RD3), so Ω₅ = G₇(Ω₅ = 0)/(4π), which is the closed form of the Result. The number G₇(Ω₅ = 0) = 513467.625456446519397870844712… equals the quadrature of the integrand of closed forms to 39 digits (RD6); the closed form equals the grids to their digits (RD7).

## The structure of the series in ln 2

All the ln 2 in Ω₀ … Ω₅ is that of exp(−4 ln 2/s): with ε = 1/s,

    s (1 + Ω₀ε + Ω₁ε² + … ) = s e^{−4 ln 2 ε} F(ε),   F(ε) = 1 + f₂ε² + f₃ε³ + f₄ε⁴ + f₅ε⁵ + f₆ε⁶ + O(ε⁷),

with f₂ … f₆ free of ln 2 (check LN2, sympy: d/d ln 2 of each coefficient is 0 through ε⁶). So W²a² = 2(D − 3) 2^{−2/(D−3)}F(1/(2(D − 3))), or W² ã² = s F(1/s) with ã = a 2^{1/(D−3)}. In the integrals this is the variable lt = ln(2 cos θ) of the masters, in which no ln 2 appears.
This is a statement about the six terms; nothing is claimed for the series beyond them.

## What each status word covers

*Derived* (by hand, each step checked): steps 1 to 10, given E1, E2, E3 (#209), the closed forms of Y₀ … Y₃, L₀ … L₃ and Ω₀ … Ω₄ (#26, #208, #220), the Beta-function moments, Parseval's identity, and **seven classical Euler sums**
(ΣH_n/n³, ΣH_n/n⁵, ΣH_n²/n⁴, ΣH_n^{(2)}/n², ΣH_n^{(2)}/n⁴, Σ(H_n^{(2)})²/n², ΣH_nH_n^{(2)}/n³), which RD5c checks numerically at 40 digits and which are not proved here.
*Checked symbolically*: XS, RD2, RD3, the exact assembly in RD5c, LN2. *Numerical evidence*: RD1, RD4 to RD7, E1_7 … G7R, A7, and the agreement of G1, G2 and #195's fit (X3), which use none of the reduction.
Not derived here: E5's and E2W's boundary terms are argued as E3's in #209 and checked on the hierarchy; the finite-part calculus beyond its statement (step 7), whose series RD4 checks numerically; the Euler sums.
**The reduction is not independent of the closed-form evaluation (A7 shares E1 … E5 and the closed forms); G1 and G2, which use none of them, are the independent check**, and #195's refit shares no machinery with either.

## Defects and discrepancies

- **The task text and one of its constants.** It names π²ζ(5) among the weight ≤ 6 constants; the weight of π²ζ(5) is 7. It is in Q9 (weight ≤ 7), where PSLQ gives it a zero coefficient.
- **The README of the raw folder said d ≤ 340; the code has 330** (`PSLQ_DIGITS`, committed before the README; amendment A2). The runs use 330, every budget is below 318, so no criterion changes.
- **First versions that failed or were stopped** (development, all logs kept): `dev_12` had E123_7 relative to G₇, which is ≈ 0 at the root (the identity held; the check now divides by its largest term, `dev_14`);
  `dev_21_t1_wrong_Om0` … `wrong_Om4` died at their output with `KeyError: 'resid2'` and `dev_21_t1_stale_cache` with `KeyError: 2` (fixed, `dev_21b_*`); `dev_22_t2_*` on the control grids had PQ5 fail in the nominal run for a budget above d = 88 (the T2 controls use the grids of record);
  `dev_23_rd_nominal` was stopped after 28 minutes inside the Euler-sum check, which is now the part `euler7` (amendment A2). I used `pkill -f` once, with a pattern naming only my own `reduce7` job (its unit was already stopped); it matched nothing else.
- **Controls that were silent where I expected a firing** (declared results, raw README A1 to A3): `e5_const` (the constant of W is zero: W(π/2) = 0 identically), `wrong_c7` (E1 … E5 hold for any table), `rd_kit` and `rd_Dstep` for RD6 and RD7 (a wrong primitive is still an exact step I ← I − dF with const += B[F]; RD2 sees it), and `e2w_g` reaches only E2Wa.
- **The reduction's two new steps** (E2W, E5) and the finite-part chain rest on boundary terms argued as E3's in #209 and checked on the hierarchy and by quadrature; none has been re-derived independently (#264 does that for E4). The seven Euler sums are classical results checked here to 34 digits, not proved.
- **The closed form was in hand before the real PSLQ ran** (amendments A1, A2): PSLQ confirms known integer coefficients to 330 digits. The bases were declared before the reduction was coded.
- **The bars of R6, R7 and X3 were set knowing a development Ω₅** (a 148-digit grid): 1e-4, 2500 and 1e-22 (declared in the raw README). R6 and R7 do not see an error of Ω₅ below about 4 and 0.13; X3 does not see an error below 1e-22.
- **Shared machinery.** G1, G2 and B1 … B3 use the same `PanelB`, table and cached chain; the independent check of Ω₅ is #195's refit (X3, 7.5e-28) and the closed-form reduction (RD7, 418.9 and 382.7 digits), which share with the grids only the table of c_m (E0).
  The cached chain equals the uncached chain bit for bit (CACHE).

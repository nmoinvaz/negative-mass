# #220: Ω₄ from the order-6 matching integral

Status: Ω₄ is **computed** from G₆ = 0 by method B on two grids at two precisions, agreeing to 209.4 digits:
Ω₄ = −12794.4998600231946327920124203810859338447198413163781740569024444356749414644670831342162805609723015636824062925937864476888763120695778347255155553…,
it is **identified** by PSLQ at 200 digits in the weight ≤ 5 basis (a relation of height 469440, residual 6.4e-196) as

    Ω₄ = −10432 ζ(5) − 160π² ζ(3) + 640 ln²2 ζ(3) − (64/9)π² ln³2 − (88/45)π⁴ ln 2 − (128/15) ln⁵2 − (2432/45)π⁴ + 1536 ln 2 ζ(3) + 3072 ζ(3),

and the same integral is reduced **by hand** to that closed form, so the identification has a **derivation** (each step checked symbolically, or numerically at 30 to 45
digits, checks RD0 to RD8; the closed form comes out of `reduce6` and is not typed in). The reduction uses the exchange identities E1, E2, E3 of #26, #208 (re-derived in #209) and one new
step, E4, that removes Y₄ and L₄. Every number is in `order6_output.txt` (the runs of record, raw folder `negative-mass-raw/I220-rpi5`, README first) unless it names another raw file.
No sentence here claims novelty; the prior equations used are README.md's order equations (matched.py, order3.py … order5.py), #26's and #208's closed forms and exchange identities,
Clausen's functions, the Beta-function moments of ln cos, Parseval's identity and Euler's sum Σ H_n/n⁴ = 3ζ(5) − ζ(2)ζ(3).

## Result

| evaluation | method | working digits | Ω₄ against the closed form |
|---|---|---|---|
| G1 | the hierarchy from its definitions, N = 170, w = 1/4, V = 38 (152 panels) | 400 | 250.7 digits |
| G2 | same, N = 150, w = 3/10, V = 36 (120 panels) | 360 | 209.4 digits |
| B3, B2, B1 | N = 100, 60, 40 (reported) | 260, 200, 120 | 148.5, 77.8, 48.3 digits |
| A6 | closed forms of Y₀ … Y₃, L₀ … L₃, E1 … E4, quadrature at 40 digits | 40 + 80 guard | 35.9 digits |

- **Ω₄** is the value above (G1; G2 agrees to 209.4 digits, |difference| 5.25e-206). The difference between evaluations is a spread over grid and precision, an estimator range,
  **not a rigorous bound**; nothing here is a rigorous enclosure. G₆(Ω₄ = 0) = 4πΩ₄ = −160780.4270664180246996853460695….
- **In 1/(D − 3).** Ω₄/16 = −652ζ(5) − 10π²ζ(3) + 40 ln²2 ζ(3) − (4/9)π² ln³2 − (11/90)π⁴ ln 2 − (8/15) ln⁵2 − (152/45)π⁴ + 96 ln 2 ζ(3) + 192ζ(3) = −799.65624125144966454950077627…,
  against the fit's −12794.499860023/16 = −799.6562412514375 (14 digits). The series is

  ```math
  W^2a^2 = 2(D-3) - 4\ln 2 + \frac{\Omega_1/2}{D-3} + \frac{\Omega_2/4}{(D-3)^2} + \frac{\Omega_3/8}{(D-3)^3} + \frac{\Omega_4/16}{(D-3)^4} + O((D-3)^{-5}),
  ```
  with Ω₁/2 = 5.2116801893692586, Ω₂/4 = 18.592346681983329, Ω₃/8 = −81.605013446257016, Ω₄/16 = −799.65624125144966.
- **The planted orders by the same code.** `omega_chain` finds Ω₀, Ω₁, Ω₂, Ω₃, Ω₄ in turn from G₂ … G₆ = 0 with no closed form of Y_k or L_k and no input value; Ω₀ … Ω₃ against −4 ln 2, 2π²/3 + 8 ln²2,
  80ζ(3) − (8/3)π² ln 2 − (32/3) ln³2 and #208's closed form: G1 250.8, 250.3, 250.0, 250.0 digits and G2 209.7, 209.1, 209.1, 209.6 (the issue's bar: 100). The slope of G_k in Ω_{k−2} is −4π
  to 251.4 (G1) and 210.7 (G2) digits at every k = 2 … 6. With Ω₀ … Ω₃ shifted by 1/100 in the orders above them, the next order moves by 0.027626, 0.027726, 0.027726, 0.027726 (controls `wrong_Om0` … `wrong_Om3`, on the control grids).
- **Rounding and grids.** The grids follow digits ≈ N log₁₀(7.7/w) of #208 at every order including the sixth: N = 100, w = 1/4 gives 148.8 digits (B3), N = 150, w = 3/10 gives 209.7, N = 170, w = 1/4 gives 250.8.
  Run times (rpi5, three processes at nice 10): G1 35 min 41 s, G2 19 min 20 s, B1 … B3 7 min 32 s; peak resident sets 443, 318 and 210 MB.
- **#195's fit** (the issue's bar, "agrees to its digits"): `w2_output.txt` prints Ω₄ = −12794.499860023 (14 digits, the range over its six cells 1.5e-30). Ω₄ differs from it by 1.9463e-10, which is inside half a unit of its
  last printed digit (X1, bar 5e-10). A refit of #195's 64 values of W² (qd, 60 digits, D = 50 … 10⁴) with Ω₀ … Ω₃ fixed at the closed forms, by polynomial extrapolation in 1/s through the K largest-s points,
  K = 14, 18, 22, 26, 30, 36, 44, gives Ω₄ to 5.4e-32, 4.3e-31, 1.7e-30, 4.4e-30, 8.4e-30, 1.6e-29, 2.4e-29 (X2, bar 1e-25; an estimator range). Two methods that share nothing agree to 2.4e-29 or better, 32 to 35 digits relative.
- **Remainder** (rem = W²a² − [s − 4 ln 2 + Ω₁/s + Ω₂/s² + Ω₃/s³ + Ω₄/s⁴], s = 2(D − 3), against #195's qd W², this Ω₄):

  | D | s | rem | rem·s⁵ |
  |---|---|---|---|
  | 1000 | 1994 | 1.35484e-12 | 42708.395 |
  | 2000 | 3994 | 4.11096e-14 | 41781.483 |
  | 5103 | 10200 | 3.73349e-16 | 41220.725 |
  | 10000 | 19994 | 1.28456e-17 | 41044.187 |

  rem·s⁵ approaches Ω₅ (the distance falls like 1/s). The two-point extrapolation from D = 5103 and 10⁴ gives Ω₅ ≈ 40860 (R6: an estimate, between 4.0e4 and 4.2e4 as declared, not a bound). **The remainder is Ω₅/s⁵:
  the six-term series with this Ω₄ leaves the D⁻⁵ remainder.** R5 and R6 see an error of Ω₄ of about 0.05 or more, not less (a shift of 0.0277 passes both); X1 and X2 are the sharp checks.
- **PSLQ.** `mp.pslq` at d = 200 digits (tolerance 1e-188, maxsteps 10⁶, Ω₄ of G2), bases declared before the run (raw README), each with a planted relation recovered exactly and a perturbed basis
  (a constant + √2·10⁻³) that finds nothing:

  | basis | constants besides Ω₄ | n | H | budget | result |
  |---|---|---|---|---|---|
  | P1 pure weight 5 | π⁴ ln 2, π² ln³2, ln⁵2, π²ζ(3), ζ(3) ln²2, ζ(5) | 7 | 10⁹ | 66.1 | not found |
  | P2 weight ≤ 5 | P1 and π⁴, π² ln²2, ln⁴2, ζ(3) ln 2, ζ(3), π² ln 2, ln³2, π², ln²2, ln 2, 1 | 18 | 10⁸ | 150.4 | **found** |
  | P3 weight ≤ 5 with Li₄(½), Li₄(½) ln 2, Li₅(½) | P2 and these three | 21 | 10⁷ | 154.3 | **found**, the same relation |
  | P4 pure weight 5 with Li₄(½) ln 2, Li₅(½) | P1 and these two | 9 | 10⁹ | 84.7 | not found |

  P2 and P3 return (−45, −88, −320, −384, −7200, 28800, −469440, −2432, 0, 0, 69120, 138240, 0, …, 0): 45Ω₄ = −88π⁴ ln 2 − 320π² ln³2 − 384 ln⁵2 − 7200π²ζ(3) + 28800 ζ(3) ln²2 − 469440ζ(5) − 2432π⁴ + 69120 ζ(3) ln 2
  + 138240ζ(3), residual 6.42e-196, with zero coefficients on π² ln²2, ln⁴2, every constant of weight ≤ 2, Li₄(½), Li₄(½) ln 2 and Li₅(½). Its budget is below d, and it is the closed form above. **The relation was
  already derived by hand when PSLQ ran** (below; raw README A1): the PSLQ run confirms known integer coefficients to 196 digits, and was not a blind search. The pure weight 5 bases fail because the lower
  weights are there: the weight-4 part is −(2432/45)π⁴ + 1536 ln 2 ζ(3), the weight-3 part 3072ζ(3).
- **Weights.** Ω₄ is not of uniform weight (as Ω₃ was not). Its weight-5 part is −12502.597…, its weight-4 part −3984.621… and its weight-3 part 3692.719…. Nothing is claimed about Ω₅ or the series.

## The order equations at order 6

README.md's equations (check E0: the ε-expansion of the master equation, sympy, residual 0, now through k = 6):

    dY_k/dl = F_{k−1} = Σ_{j<k} L_{k−1−j} l^j/j!,    (R − 1)L_k′ − L_k/2 = h_k = (1 − 1/R)L_{k−1} + Σ_{j<k} Y_j c_{k−j} + Y_k/2,
    L_k = −2Ω_{k−1} + (R − 1)^{1/2}C_k,    G_k = C_k(∞) = ∫₁^∞ (h_k − Ω_{k−1})(R − 1)^{−3/2}dR = 0,

with c₁ … c₅ as in `order5.py` and the new

    c₆ = Ω₀l⁵/120 + Ω₁l⁴/24 + Ω₂l³/6 + Ω₃l²/2 + Ω₄l + Ω₅ + l⁶/1440 + l⁵/120 − l⁵/(120R),

the coefficient of ε⁶ of (1 + εΩ)e^{εl} − (½ + ε)(1 − 1/R)e^{−εl} − e^{−εl}/(2R). G_k is affine in Ω_{k−2} with slope −4π; Ω₅ drops out of G₆ (it is h₆(1), subtracted), so **G₆ = 0 gives Ω₄**. Method B, as in omega2.md and omega3.md:
Y_k from the left, the tail C_right from the right, L_k = −2h_k(1) − u C_right, G_k = C_right(0), on Chebyshev panels in v = √l (`PanelB` of order4.py, shared; `PanelK` replaces its two sums by `mp.fdot` and caches l^j/j! and
c_m on the nodes, which changes timing and the last rounding only). What the evaluations share: the table (E0), the analytic-at-throat condition, `PanelB` and mpmath; G1 and G2 differ in N, w, V and precision. #195's fit (the full equation
in a Rust crate, no matching) shares none of it and is the independent reference; the closed forms of Ω₀ … Ω₃ are the planted truth.

## The reduction by hand

Notation as omega2.md and omega3.md: θ ∈ [0, π/2), R = sec²θ, S = sin²θ, u = tan θ, ℓ = ln cos θ = −l/2, q = ln 2, C = Cl₂(π − 2θ), D₃ = Cl₃(π − 2θ) + ¾ζ(3), Lc2 = ∫₀^θ ℓ², J_n = ∫₀^{π/2}ℓⁿdθ, Ω₀ = −4q.
The derivatives close: ℓ′ = −u, u′ = 1 + u², C′ = 2ℓ + 2q, D₃′ = 2C, Lc2′ = ℓ² (RD0). Closed forms Y₀ … Y₃, L₀ … L₃ as in #26 and #208 (L₃ through Lc2).

1. **E1, E2, E3 at k = 6, m = 5** (#209; they need G_j = 0 for j ≤ 5, which the lower orders give):
   G₆ = 8∫H₅ℓ/S + 2∫(F₅ − L₅) + 2∫(S′₆ − Ω₅)/S − 4∫(ℓ − Ω₀)F₄, with H₅ = S L₄ + Σ_{j≤4}Y_jc_{5−j}, S′₆ = Σ_{j≤4}Y_jc_{6−j}. Y₅ and L₅ are gone (E3 removes Y₅). What is left is Y₀ … Y₄ and L₀ … L₄,
   and **Y₄ and L₄ have no closed form** (#26: Y₄ is only an integral). The reduction goes on.
2. **L₄.** It enters through ∫L₄[8ℓ + 2l − 4(ℓ − Ω₀)] = 4Ω₀∫L₄ (l = −2ℓ), and E2 at m = 4 (it needs G₄ = 0) gives 4Ω₀∫L₄ = 8Ω₀∫h₄ℓ/S, with h₄ = S L₃ + Σ_{j≤3}Y_jc_{4−j} + Y₄/2. L₄ is gone; ½Y₄ stays.
3. **E4: Y₄.** Collect the Y₄ terms: 8ℓc₁/S (in H₅), 2c₂/S (in S′₆), 4Ω₀ℓ/S (from step 2). With c₁ = Ω₀ − 3ℓ − S and c₂ = Ω₁ − 2Ω₀ℓ + ℓ² − 2Sℓ they are Y₄w₄, w₄ = (2Ω₁ + 8Ω₀ℓ − 22ℓ²)/S − 12ℓ.
   With Ω₀ = −4q, W₄ = cot θ (22ℓ² − 8Ω₀ℓ − 2Ω₁) + 16C satisfies W₄′ = w₄ (the θ terms cancel: 8Ω₀θ + 16(2qθ − C) = −16C). Since dY₄/dθ = 2uF₃ and W₄ u = 22ℓ² − 8Ω₀ℓ − 2Ω₁ + 16uC,

   ```math
   \int_0^{\pi/2} Y_4 w_4\,d\theta = -2\int_0^{\pi/2}\big(22\ell^2 - 8\Omega_0\ell - 2\Omega_1 + 16\,u\,C\big)F_3\,d\theta,\qquad F_3 = L_3 + L_2 l + \tfrac12 L_1 l^2 + \tfrac16 L_0 l^3 .
   ```
   Boundary terms [W₄Y₄]: at 0, W₄ = −2Ω₁/θ + O(θ) and Y₄ = O(θ²), so O(θ); at π/2, W₄ = O(x ln²x) and Y₄ is polynomial in l (#209, F3), so O(x lⁿ) → 0. Checked on the hierarchy, with no closed form: E4_6 to 44.2 digits, and the reduced
   form of G₆ (Y₀ … Y₃, L₀ … L₃ only) against G₆: G6R to 43.9 digits, at Ω₄ = its root, 0 and 1/2 (E-grid 36, 1/2, 16; E1_6, E2_5, E2_4, E3_6: 45.7, 43.7, 43.7, 43.7 digits). **So G₆ = ∫₀^{π/2}Ĩ₆dθ with Ĩ₆ built from the closed forms.**
4. **The decomposition** (sympy, RD1: 45 digits against the numeric integrand of closed forms). With E = 2ℓ³ + 72qℓ² + (64q² + 2Ω₁)ℓ + Ω₂ and the u-Laurent form (1/S = 1 + u⁻²):

   Ĩ₆ = a_CC C² + a_CD C D₃ + a_CL C Lc2 + c_C C + c_D D₃ + c_L Lc2 + c_R,

   a_CC = 1024u²(3ℓ + 2q + 1), a_CD = −512u, a_CL = −2048u²,
   c_C = 64E/u + (128/3)u(121ℓ³ + (420q + 33)ℓ² + (288q² + 120q − 6Ω₁)ℓ − 6qΩ₁ − 3Ω₁ + 3Ω₂) + 256Ω₁θu²,
   c_D = −32E/S − 32(6ℓ² + 48qℓ − 2Ω₁), c_L = −256u(8ℓ² + 16qℓ − Ω₁),

   and c_R, a polynomial in ℓ with u⁻², u⁰ and θu parts (`dev_decomposition_of_the_integrand.txt`), whose top terms are (4/45)·862ℓ⁶/S and 8Ω₄ℓ/S. The slope −4π is its Ω₄ term, 8Ω₄∫ℓ/S = −4πΩ₄.
5. **Finite parts.** Several pieces diverge at π/2 (u²C Lc2 ∼ ln x/x; u θ, u ℓⁿ, u² ℓⁿ), and cancel only in the sum. ⟨f⟩ is the constant term (x⁰ ln⁰x, x = π/2 − θ) of ∫₀^{π/2−x}f as x → 0; for an integrable f it
   is the integral, ⟨·⟩ is linear, and for a primitive F, ⟨F′⟩ = B[F] := FP F(π/2) − FP F(0). So every integration by parts below holds as an identity of finite parts, and the sum of the finite parts is G₆.
   The boundary terms are the constant terms of the series of the generators (order 9; RD4), e.g. u C Lc2 → 2(1 − ln x − q)J₂, so FP = 2J₂(1 − q).
6. **The three quadratic terms.** Each is the derivative of a product, M1: F = u f(ℓ) C², M2: F = −ℓ C D₃, M3: F = u C Lc2 (RD2: the derivatives and the coefficients, sympy residual 0):

   ⟨u²(f − f_ℓ)C²⟩ = B[u f C²] − ⟨f C²⟩ − 4⟨u f(ℓ + q)C⟩, f = 3072ℓ + 2048q + 4096 (so that a_CC = u²(f − f_ℓ)); B = 0.
   ⟨u C D₃⟩ = B[−ℓ C D₃] + ⟨ℓ(2ℓ + 2q)D₃⟩ + 2⟨ℓ C²⟩; B = 0.
   ⟨u²C Lc2⟩ = B[u C Lc2] − ⟨C Lc2⟩ − 2⟨u(ℓ + q)Lc2⟩ − ⟨u ℓ²C⟩; B = 2J₂(1 − q).

   They add to c_C: −4u f(ℓ + q) + 2048uℓ²; to c_D: −512ℓ(2ℓ + 2q); to c_L: 4096u(ℓ + q); and constants through ⟨C²⟩, ⟨ℓC²⟩, ⟨C Lc2⟩. No θ enters.
7. **The D₃ term.** c_D = c₀(ℓ) + c₁(ℓ)/S with c₀ of degree 2 and c₁ = −32E. Its primitive is A = Σa_bP_b − Σβ_bF_b with P₀ = θ, P₁ = (C − 2qθ)/2, P₂ = Lc2, F_b = ∫ℓ^b/S = −ℓ^b cot θ − bP_{b−1}
   (RD2: A′ = c_D, residual 0), so ⟨c_D D₃⟩ = B[A D₃] − 2⟨A C⟩. B = (7/4)ζ(3)A(π/2) = 0: the D₃ coefficient integrates to zero. **The cot θ terms of −2A C cancel the 1/u part of c_C exactly (64E/u·C: RD3)**, as the cot θ and θ terms did at order 5.
   What is left of −2⟨A C⟩ is ⟨θC⟩, ⟨C²⟩ and ⟨C Lc2⟩ with constant coefficients.
8. **The C terms.** c_C has u¹ℓ^b (b ≤ 3) and u²θ parts. ⟨uℓ^bC⟩ = (2J_{b+2} + 2qJ_{b+1})/(b + 1) (by parts, C′ = 2ℓ + 2q, boundary 0; #208). ⟨u²θC⟩ = B[uθC] − ⟨θC⟩ − ⟨uC⟩ − 2⟨uθ(ℓ + q)⟩ with
   B[uθC] = π(1 − q) and ⟨uθ(ℓ + q)⟩ = B[−½θ(ℓ + q)²] + ½⟨(ℓ + q)²⟩ = −πq²/4 + ½(J₂ + 2qJ₁ + q²J₀).
9. **The Lc2 and the rest.** ⟨u g(ℓ)Lc2⟩ = B[−G Lc2] + ⟨Gℓ²⟩, G = ∫g dℓ with G(0) = 0, B = 0. For c_R: ⟨ℓ^b/S⟩ = −bJ_{b−1}, ⟨u⁻²ℓ^b⟩ = −bJ_{b−1} − J_b, ⟨uθℓ^b⟩ = B[−θℓ^{b+1}/(b + 1)] + J_{b+1}/(b + 1) with B = 0.
   **J₆, and with it ζ(3)², cancels** (the ℓ⁶ terms add up to 862ℓ⁶/S, which gives −6J₅; RD3).
10. **The master integrals** (RD5: each against quadrature at 40 digits, 35.1 to 40.0 digits):
    - J_n from B(a) = ∫cos^aθ dθ = (√π/2)Γ((a + 1)/2)/Γ(a/2 + 1) = Σ J_n aⁿ/n!; J₅ = (π/96)(−19π⁴ ln 2 − 1080ζ(5) − 60π²ζ(3) − 720 ln²2 ζ(3) − 40π² ln³2 − 48 ln⁵2) is where ζ(5) enters first.
    - ∫C² = π⁵/360 (Parseval for sin 2nθ, #209); ∫ℓC² = −q∫C², since ℓ + q = C′/2 and ∫C′C² = [C³/3] = 0; ∫θC = Σ((−1)^{n+1}/n²)((−1)^{n+1}π/(4n)) = (π/4)ζ(3); ∫uC = π³/12 (#26).
    - **∫ℓ²D₃ dθ = ¾ζ(3)J₂ + (π/4)(2qζ(4) + 5ζ(5) − 2ζ(2)ζ(3)).** From ln cos θ = −q − Σ_{k≥1}(−1)^k cos 2kθ/k, the coefficient of cos 2nθ in ℓ² is a_n = ((−1)ⁿ/n)(2q + 2H_n − 1/n)
      (from 2q b_n + Σ_m b_m b_{m+n} + ½Σ_{m<n}b_m b_{n−m}, b_k = (−1)^k/k, Σ_m 1/(m(m + n)) = H_n/n, Σ_{m<n}1/(m(n − m)) = 2H_{n−1}/n), and D₃ = ¾ζ(3) + Σ(−1)ⁿcos 2nθ/n³, so by Parseval
      ∫ℓ²D₃ = ¾ζ(3)J₂ + (π/4)Σ_n(2q + 2H_n − 1/n)/n⁴, and Σ H_n/n⁴ = 3ζ(5) − ζ(2)ζ(3) (Euler). This is the second place where ζ(5) enters.
    - ∫C Lc2 = [½D₃ Lc2] − ½∫ℓ²D₃ = (7/8)ζ(3)J₂ − ½∫ℓ²D₃, with D₃(π/2) = (7/4)ζ(3), Lc2(π/2) = J₂, D₃(0) = Lc2(0) = 0.
11. **Assembly** (sympy, exact). With Ω₁, Ω₂, Ω₃ at their closed forms and Ω₄ = 0, the sum of the finite parts is G₆(0) = 4πΩ₄ᶜ with Ω₄ᶜ the closed form at the top; its slope in Ω₄ is −4π (RD3), so G₆ = 4π(Ω₄ᶜ − Ω₄). **Derived:** G₆(Ω₄ = 0) = −160780.42706641802469968534606….
    The closed form against the grids: G1 250.7 digits, G2 209.4 (RD7). The integer relation of PSLQ is 45 × this closed form's coefficients (RD8).

**What each status word covers.** *Derived* (by hand, each step checked): steps 1 to 11, given E1, E2, E3 (#209: derived, with their boundary terms), the closed forms of Y₀ … Y₃ and L₀ … L₃ and Ω₀ … Ω₃ (#26, #208),
Parseval, Euler's sum and the Beta-function moments. *Checked symbolically*: RD2, RD3 (the by-parts identities, the cancellations), RD8's coefficient comparison. *Numerical evidence*: RD0, RD1, RD4 to RD7, E1_6 … G6R, and the agreement of G1,
G2, A6 and #195's fit. Not derived here: E4's boundary terms at π/2 (argued as E3's in #209, checked on the hierarchy), the finite-part calculus beyond its statement (steps 5), and the series of the generators (RD4 checks them numerically).
**A6 is not independent of the derivation** (same E1 … E4 and closed forms); **G1 and G2, which use none of them, are the independent check**, and #195's fit and its refit (X1, X2) share no machinery with either.

## Controls and checks

`order6.py` has 30 checks in the runs of record: E0; PL0 … PL3, SL, M4, X1, X2, R5, R6 (`values`, `remainder`); E1_6, E2_5, E2_4, E3_6, E4_6, G6R (`exch6`); Q1 … Q4 (`pslq`); RD0 … RD8 (`reduce6`); PLW is the control's own statement
and reports only under `wrong_Om*`. **All pass** (exit statuses 0, 0, 0, 0, 0, 0). Each has a control (`--mutate=NAME`, 25 of them), whose targets were computed in development runs before being declared (raw README and its amendment A1).
**Controls of record: 27 mutated runs and 4 nominal runs, 31 of 31 on target** (`ctl_t1_*`: 8 grid-dependent mutations on the control grids; `ctl_t2_*`: 4 post-processing mutations on the saved Ω's; `ctl_ex_*`: 7 mutations of the hierarchy checks;
`ctl_red_*`: 8 mutations of the reduction, each family with its nominal run). On target means that the set of checks that fire is exactly the declared set, so every declared target fired and no declared-silent check fired
(`order6_output.txt`, last section). The declared silences that are results:
- the exchange checks E1_6 … E4_6 hold for any table: `wrong_c6` is silent on all six, and only E0, X1 and X2 see it;
- `wrong_Om3` and `wrong_c6` (an error of 0.028 and 0.015 in Ω₄) are silent on R5 and R6; `coarse_B` is silent on X1; `short_tail` and the `wrong_Om*` inputs, which both grids share, are silent on M4;
- `rd_kit` (the sign of ln 2 in C's expansion at π/2) fires only RD4: the closed form does not depend on it.

## Defects and discrepancies

- **R5 and R6 are coarse.** A wrong Ω₄ of 0.0277 (`wrong_Om3`) or 0.0152 (`wrong_c6`) passes both (sensitivity about 0.05); X1 and X2 fire. `wrong_Om4_series` needed a shift of 1 to fire them.
- **Controls that were silent where I expected a firing** (results, raw README A1): `rd_kit` (the sign of ln 2 in C's expansion at π/2) and a first `rd_kit2` (D₃(π/2) = 3ζ(3)/2) change nothing in the closed form. The first because the finite parts of u θ C and u C Lc2 move together
  and cancel (8J₂ = πΩ₁/2), the second because the D₃ coefficient integrates to zero. `coarse_B` is silent on X1 (the grid's error 1e-11 is below the 14 digits of the fit).
- **First versions that failed** (development, all logs kept): PSLQ P2 and P3 on the control grids at d = 113 with H = 10⁸ and 10⁷ have budgets above d, so the perturbed basis found a chance relation (`dev_04`); the control grids now use H = 10⁵ and 10⁴.
  `wrong_Om4_series` at +1/100 fired nothing (`dev_14`). RD8 first compared the relation's vector with the closed form's coefficients with the wrong sign (`dev_21`).
- **The closed form was in hand before the real PSLQ ran** (above); PSLQ is a confirmation of integer coefficients to 196 digits, not a discovery. The bases were declared before the reduction was coded.
- The bars of R5 and R6 were set knowing the order of magnitude of #208's table (declared in the raw README); X2's bar was set knowing a development extrapolation (5e-32 to 2e-29).
- A shared limit: G1, G2 and the control grids use the same `PanelB` and table; the independent check of Ω₄ is #195's fit and the derivation.

## What README.md's "Result" and "Limits", #195, #210 and Paper III should now say

- **README.md, Result.** Add the sixth term of the series: Ω₄/16 = −652ζ(5) − 10π²ζ(3) + 40 ln²2 ζ(3) − (4/9)π² ln³2 − (11/90)π⁴ ln 2 − (8/15) ln⁵2 − (152/45)π⁴ + 96 ln 2 ζ(3) + 192ζ(3) = −799.65624125144966…, **derived** from the order-6 matching integral
  (this note: G₆ = 0, E1, E2, E3 and a new step E4 that removes Y₄ and L₄, then finite-part integration by parts to Beta-function moments, Parseval and Euler's Σ H_n/n⁴), confirmed to 209 and 250 digits by two evaluations that use none of the reduction,
  to 14 digits by #195's fit and to 32 to 35 digits by a refit of its W². The remainder statement becomes: the six terms leave Ω₅/(2(D − 3))⁵ with Ω₅ ≈ 4.09·10⁴ (an estimate from rem·s⁵; the number is not computed).
  Ω₄ is not of uniform weight; it has no Li₄(½) or Li₅(½) term (PSLQ P3).
- **README.md, Limits.** Replace "Ω₄ is a fit (#195), not computed or identified" by: Ω₅ is an estimate (≈ 4.09·10⁴), not computed.
- **#195** (W² to 60 digits and the fit): its Ω₄ = −12794.499860023 is now equal to the derived value to the 14 digits it prints.
- **#210 (Paper III).** The sixth term can be quoted as derived. Its derivation needs E1 … E3, the new E4, the finite-part bookkeeping of step 5, Parseval, Euler's sum and the Beta-function moments up to J₅; ζ(5) enters through J₅ and ∫ℓ²D₃.
  `omega_graph.py` and `paper_numbers.py` carry Ω₄ as #195's fit; the derived value replaces it (the coordinator's files). The paper's "how this work was done" section should record that the closed forms and their checks were done by a language model with sympy and mpmath, as AGENTS.md asks.

## Files

- `order6.py` (parts `sym`, `evals`, `values`, `exch6`, `pslq`, `remainder`, `reduce6`; `grid` for development; `--mutate=NAME`, `--cgrids`, `--grids=`, `--save`, `--load`, `--w2raw`, `--no-real`), `order6_output.txt` (the runs of record concatenated, then the control summaries),
  this note. It imports `PanelB`, the closed forms and `Lc2_poly` from `order4.py`, which is unchanged.
- Raw: `negative-mass-raw/I220-rpi5` (README with the criteria, controls and amendment A1; `dev_*` the development runs; `run_01` … `run_06` the runs of record; `ctl_*` the controls of record; `G1.json`, `G2.json`, `B.json` the saved Ω's;
  `declared_*.json`; `scratch_prototypes/` the sympy and mpmath prototypes of the reduction and of X2; `assemble_output.py`).

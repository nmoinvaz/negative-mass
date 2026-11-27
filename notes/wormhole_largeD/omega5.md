# #263: Ω₅ from the order-7 matching integral

DRAFT, filled in as the runs of record finish.

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

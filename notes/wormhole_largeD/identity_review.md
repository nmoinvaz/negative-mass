# E4, E2W and E5 derived again, and the masters of Ω₄ and Ω₅ checked by quadrature (#264)

November 2026. Issue #264 (milestone G5, Paper III), the independent review of #220's identity E4 and, by the coordinator's widening, of #263's E2W and E5.
The derivations of §2 to §5 were written from the order equations and the definitions alone and committed (raw `I264-rpi5`, `derivations.md`, 38869df)
before any script of this review ran and before `omega4.md` and `omega5.md` were read beyond their statement lines (§1); §7 is the comparison
made afterwards. The script is [`identity_review.py`](identity_review.py), its output [`identity_review_output.txt`](identity_review_output.txt),
the raw logs `I264-rpi5` of negative-mass-raw.

Status words as in AGENTS.md: *proved* (by hand, each step with its boundary terms), *checked symbolically* (sympy residual 0), *numerical evidence*
(agreement of integrals), *identified* (a closed form found by PSLQ). A digit count is −log₁₀ of a relative difference. It is not an error bar.

## Result

E4 (at k = 6 on Y₄ and at k = 7 on Y₅), E2W (for L₄, weights 12ℓ² and g₄) and E5 (on Y₄) are **proved** with every boundary term and every hypothesis (§3, §4),
**checked symbolically** (§5) and **checked numerically** on a Taylor-jet hierarchy of my own that shares no code with `PanelK` (§6). No sign, constant or
result differs from `omega4.md` or `omega5.md` (§7).

| item | status | evidence (`identity_review_output.txt`) |
|---|---|---|
| E4, k = 6: ∫Y₄w₄ = −2∫(22ℓ² − 8Ω₀ℓ − 2Ω₁ + 16uC)F₃ | proved; symbolic SY_W4, SY_w4; numerical E4_6 | 63.5 digits; 91.7 with its boundary term; bar 30 |
| E4 one level up (Y₅, F₄) | proved; numerical E4_7 | 62.3; 91.9 |
| E2W: ∫g₄L₄ = −2∫P h₄/S, and for 12ℓ² | proved; SY_g4, SY_P, SY_Pa; E2Wb, E2Wa | 62.6, 62.6; 92.1, 92.1 |
| E5: ∫Y₄w₄′ = −2∫W u F₃, with W(π/2) = 0 | proved; SY_w4p, SY_W, SY_K; E5 | 62.7; 91.9 |
| E1, E2, E3 at the orders used (k = 6, 7) | the lemmas of §2 give them; numerical E1_6 … E3_7 | 60.7 to 68.3; 91.3 to 93.6 |
| reduced G₆ and G₇ (only Y₀ … Y₃, L₀ … L₃) equal G₆, G₇ of the definition | numerical G6R, G7R; slope −4π at both | 64.5, 62.8; 91.6, 92.1 |
| closed forms Ω₄, Ω₅ (and Ω₀ … Ω₃) are roots | numerical H2, H3 on the same hierarchy, no closed form of Y or L used | G₆ = −1.7e-59, G₇ = 6.9e-57 (tail-limited, §6) |
| ∫ℓ²D₃ = ¾ζ(3)J₂ + (π/4)(2 ln 2 ζ(4) + 5ζ(5) − 2ζ(2)ζ(3)) | numerical M1, no series | 72.2 digits (bar 40) |
| ⟨lt²C²⟩ = π(¾ζ(3)² − (95/96)ζ(6)), ⟨lt³D₃⟩ = π(9ζ(3)²/16 − (11/4)ζ(6)) | identified by PSLQ at 70 digits in my own bases, then checked at 125 digits; the notes' forms checked, NM2, NM3 | 125.8 and 119.3; 57.7 and 53.7 digits |
| the other 22 masters of both notes | numerical, quadrature of the definitions | 49.3 to 58.8 digits, or equal at 58 |
| the finite-part calculus (four statements) | proved (§8); two planted splits FP1, FP2 | 59.1, 59.2 |
| the 11 steps (k = 6) and 37 integrations by parts (k = 7) of the finite-part chain | **not re-derived** | covered end to end by the row on Ω₄, Ω₅ |
| controls | 18 of 18 on target | §9 |

Findings that are new, not in `omega4.md` or `omega5.md` (§7): W(π/2) = 0 in E5 holds only at the roots of Ω₀ and Ω₁ (the constant is
−4π(Ω₁ − 2π²/3 − 8 ln²2) at Ω₀ = −4 ln 2), so E5 depends on G₃ = 0 through that constant; E2W needs G₄ = 0: with Ω₂ moved by 1e-3 its two sides agree to 0.7 digits
while E4 and E5 hold to 63; each by-parts identity is an exact identity on [0, θ_e] with its boundary term, which
the hierarchy reproduces to 91 digits.

## 1. What I read, and when

Before the derivations: AGENTS.md, the issue and its widening comment, the handovers `i209`, `i220`, `i263` of my auto-memory, `README.md` of this
directory, `exchange_review.md` (#209, used as established for E1 to E3 and the Clausen integrals). From `omega4.md` lines 68 to 83 and 93 to 98,
from `omega5.md` lines 74 to 86 and 98 to 105: the order equations, the definition of c₆ and c₇, and the statements of E4, E2W and E5. Those
lines carry short sketches: for E2W the equation (L cot θ)′ = 2h/S and the sentence "∫L_m g = −2∫P h_m/S with P′ = g tan θ, P(0) = 0", for E5 the cancellation of
∫C cot θ, and a one-line boundary remark. I saw them, so the key equation of E2W and the shape of E5's primitive are not independent. What was
written from the order equations alone: the two lemmas with their conditions (§2), the growth and throat bounds (§1 of `derivations.md`), the weights w₄, g₄,
w₄′ by tracking every term through the chain, P and W, the constant of W, every boundary term (§3). The handovers `i220`, `i263` also state w₄, W₄ and name E2W, E5.

After the declared runs of the first set (`run_02`, raw e12aee3): `omega4.md` lines 84 to 204 and `omega5.md` lines 1 to 195. Not read at any point:
`order6.py`, `order7.py`, `PanelK`, the raw folders `I220-rpi5`, `I263-rpi5`, the coordinator's `I220-mac-review`, `I263-mac-review`.

## 2. Setting, the order equations and two lemmas

Variables as in #209: R = sec²θ, u = tan θ, l = ln R = −2ℓ, ℓ = ln cos θ, S = sin²θ = 1 − 1/R, x = π/2 − θ, lt = ln(2 cos θ) = ℓ + ln 2, C = Cl₂(π − 2θ)
(C′ = 2 lt), D₃ = 2∫₀^θ C = ¾ζ(3) + Cl₃(π − 2θ). From README's exact master equation, with d/dl = R d/dR, the order equations are

  dY_k/dl = F_{k−1} = Σ_{j<k} L_{k−1−j} l^j/j!,   S L_k′ − L_k/2 = h_k = S L_{k−1} + S_k + ½Y_k,   S_k = Σ_{j<k} Y_j c_{k−j},

with c_m = [ε^m] of K(ε) = (1+εΩ)e^{εl} − ½e^{−εl} − εS e^{−εl}. In ℓ:
c_m = Σ_{i<m} Ω_i(−2ℓ)^{m−1−i}/(m−1−i)! + [1 − ½(−1)^m](−2ℓ)^m/m! − S(2ℓ)^{m−1}/(m−1)! (SY_c, c₁ … c₇ against the expansion; SY_c67, c₆ and c₇ against the
definitions printed in the two notes). So c₁ = Ω₀ − 3ℓ − S, c₂ = Ω₁ − 2Ω₀ℓ + ℓ² − 2Sℓ, c₃ = Ω₂ − 2Ω₁ℓ + 2Ω₀ℓ² − 2ℓ³ − 2Sℓ², and c_m(0, R = 1) = Ω_{m−1}.
The solution analytic at the throat is L_k = −2Ω_{k−1} + uC_k, C_k = 2∫₀^θ(h_k − Ω_{k−1})/S, G_k = C_k(π/2) = 0. **H_k:** G_j = 0 for 2 ≤ j ≤ k − 1; Ω_{k−2} free.

**Throat (T).** h_k − Ω_{k−1} = (R − 1)a_k with a_k analytic, so C_k = O(θ), L_k = −2Ω_{k−1} + O(θ²), Y_k = O(l) = O(θ²).
**Growth (G).** Under H_k, L_j, h_j, Y_j, F_j are polynomials in l at most, and C_j = O(x l^q): since G_j = 0, C_j(θ) = −2∫_θ^{π/2}(h_j − Ω_{j−1})/S, and on [π/4, π/2]
1/S ≤ 2 and l is increasing, so |C_j| ≤ 4K x l^q, uC_j = O(l^q) (tan θ ≤ 1/x), Y_{j+1} = ∫F_j dl = O(l^{q+1}). An induction on j, valid to j = 6 under H₇.
Numerically Y_k/l^k at l = 256 is −1.0, −1.48, 0.53, 1.70, 0.66, −0.86, −1.13 for k = 1 … 7 (GROW, no bar), so the bound is not sharp.

**(L cot θ)′ = 2h/S** follows from the L equation (R − 1)dL/dR − L/2 = h, since (R − 1)dL/dR = ½ sin θ cos θ L_θ (SY_L).

**Lemma Y.** If W′ = w, then ∫₀^{π/2} w Y_m dθ = [W Y_m] − 2∫ u W F_{m−1} dθ (dY_m/dθ = 2uF_{m−1}). The boundary term vanishes if W = O(1/θ) at 0 (Y_m = O(θ²))
and W(π/2) = 0 with w polylogarithmic (W = −∫_θ^{π/2}w = O(x l^s), Y_m = O(l^p)). **W + K is not allowed**: K Y_m(π/2⁻) diverges, so the constant is fixed by W(π/2) = 0.

**Lemma L.** If P′ = g tan θ, P(0) = 0, then ∫₀^{π/2} g L_m dθ = [P L_m cot θ] − 2∫ P h_m/S dθ. At 0: P = O(θ²), L_m cot θ = O(1/θ). At π/2: L_m cot θ = −2Ω_{m−1}cot θ + C_m
= O(x l^q) **requires G_m = 0**; without it L_m ≈ G_m/x, and ∫gL_m diverges. E1, E2, E3 are the cases w = 1/S (W = −cot θ), g = 1 (P = −ℓ), w = (4ℓ + 2c₁)/S = 2P_E′ with P_E = cot θ(ℓ − Ω₀) (SY_E3).

## 3. The boundary terms of every step

| step | boundary term | at θ = 0 | at π/2 | uses |
|---|---|---|---|---|
| E1 at k | [−cot θ Y_k] | 2Ω_{k−2}θ | O(x l^p) | T, G |
| E2 at m | [−ℓ L_m cot θ] | −Ω_{m−1}θ | O(x l^{q+1}) | G_m = 0 |
| E3 at m | [2P_E Y_m] | 4Ω₀Ω_{m−2}θ | O(x l^{p+1}) | T, G |
| E4 on Y₄ (k = 6) and Y₅ (k = 7) | [W₄ Y_n] | 4Ω₁Ω_{n−2}θ | W₄ = O(x l²) + 16C, C = Cl₂(2x) = O(x ln x) | Ω₀, Ω₁ roots (θ term and W₄(π/2) = 0) |
| E2W on L₄ | [P L₄ cot θ] | −4Ω₁Ω₃θ | O(x l^{q+3}) | G₄ = 0 |
| E5 on Y₄ | [W Y₄] | 4Ω₂²θ | W(π/2) = 0, O(x l^s) | Ω₀, Ω₁ roots (K = 0) |

The leading coefficients at the throat are checked symbolically (SY_bnd0, SY_bnd1) and at the two nodes with the smallest l, θ = 4.4e-4 and 2.3e-3: the ratio to the prediction is 1 to 4.6e-6 for the
11 terms (BND_throat, bar 1e-4). At the outer end b/(x l^q) is bounded by 10 at l = 30, 60, 100 (BND_outer, bar 1e3; the table is in the output). All of them tend to 0 under H_k.

## 4. The chain at k = 6 and k = 7

**k = 6.** E1 at k = 6; E2 at m = 5 (G₅ = 0); E3 at m = 5. The L₄ weight is 8ℓ + 2l − 4(ℓ − Ω₀) = 4Ω₀, and E2 at m = 4 (G₄ = 0) turns 4Ω₀∫L₄ into 8Ω₀∫ℓh₄/S. The Y₄ terms are
4Ω₀ℓ/S, 8ℓc₁/S, 2c₂/S:

  **w₄ = (2Ω₁ + 8Ω₀ℓ − 22ℓ²)/S − 12ℓ,   W₄ = cot θ(22ℓ² − 8Ω₀ℓ − 2Ω₁) + 16C − (8Ω₀ + 32 ln 2)θ,**

and the last term vanishes at Ω₀ = −4 ln 2. Lemma Y gives **E4**: ∫Y₄w₄ = −2∫(22ℓ² − 8Ω₀ℓ − 2Ω₁ + 16uC)F₃. It needs Ω₀ and Ω₁ at their roots and not Ω₂. The reduced G₆ has only Y₀ … Y₃, L₀ … L₃
(`derivations.md` §4.9).

**k = 7.** E1 at k = 7; E2 at m = 6; E3 at m = 6; the L₅ weight is again 4Ω₀ and E2 at m = 5 follows. The Y₅ terms have the same weight w₄, so E4 holds one level up (needs L₄ bounded:
Ω₂ too). The L₄ weight collects 8Ω₀ℓ (E2, m = 5) + l² (E1) + 8ℓ(ℓ − Ω₀) (E3) − 2(22ℓ² − 8Ω₀ℓ − 2Ω₁ + 16uC) (E4):

  **g₄ = −32ℓ² + 16Ω₀ℓ + 4Ω₁ − 32uC** (12ℓ² from the first three), **P = (32/3)ℓ³ − 8Ω₀ℓ² − 4Ω₁ℓ − 32(uC + ℓ² + 2 ln 2 ℓ − D₃/2)**,

with ∫u²C = uC + ℓ² + 2 ln 2 ℓ − D₃/2. Lemma L gives **E2W**: ∫g₄L₄ = −2∫P h₄/S (G₄ = 0), and for 12ℓ², P = −4ℓ³. The Y₄ terms are 8Ω₀ℓc₁/S, 8ℓc₂/S, 2c₃/S and −P/S:

  **w₄′ = [2Ω₂ + a₁ℓ + a₂ℓ² + a₃ℓ³ + 32uC − 16D₃]/S − 20ℓ² − 8Ω₀ℓ,** a₁ = 8Ω₁ + 8Ω₀² + 64 ln 2, a₂ = 32 − 28Ω₀, a₃ = −20/3.

The primitive: ∫uC/S = Z + ∫C tan θ and ∫D₃/S = −cot θ D₃ + 2Z with Z = ∫C cot θ, so **Z cancels** in 32∫uC/S − 16∫D₃/S. W is the expression of the two notes (`derivations.md` §5.9) with
W(π/2) = K₀ + K, **K₀ = −4π(Ω₁ − 2π²/3 − 8 ln²2) at Ω₀ = −4 ln 2** (SY_K), so K = 0 only at the root of Ω₁. Lemma Y gives **E5**: ∫Y₄w₄′ = −2∫W u F₃. The reduced G₇ is in `derivations.md` §5.11.

## 5. Symbolic checks (sympy residual 0)

SY_c, SY_cS, SY_c67, SY_c0 (the c_m); SY_L ((L cot θ)′ = 2h/S); SY_w4; SY_E3; SY_W4 (W₄′ = w₄); SY_g4 (g₄ from the chain); SY_P (P′ = u g₄); SY_Pa; SY_w4p (the stated form of w₄′
from the c_m and P); SY_W (W′ = w₄′ with Lc₁′ = ℓ, Lc₂′ = ℓ², C′ = 2(ℓ + ln 2), D₃′ = 2C); SY_K (K₀); SY_bnd0, SY_bnd1 (throat expansions); SY_J (J₁, J₂ from B(a)). 17 of 17.

## 6. The hierarchy and its checks

Own code, in the variable l: the order equations of §2 with c_m expanded by sympy at run time; the Frobenius recursion at the throat (230 terms), then power series about l = 1.5, 2.5, … propagated by steps of 1
(160 terms), all in Python integers at 440 bits; orders Y₀ … Y₇, L₀ … L₆; 13 s on rpi5 to l = 345. The θ integrals are Gauss–Legendre panels in v = √l (dθ = v dv/u, analytic in v), 40 nodes on panels of 0.5 to v = 18.5
(l = 342.25), at 100 digits; C and D₃ from mpmath (M0 checks them against nested quadrature, 2.8e-76). It shares with the notes only README's master equation, the closed forms of Ω₀ … Ω₅ as inputs to be tested,
mpmath and sympy.

**Planted truth and convergence.** H0: Y₁, L₁, Y₂, L₂ against README's closed forms (K(u) by quadrature) at 20 nodes, l up to 50: 1.0e-81. H1: a second setting (500 bits, 190 terms, 44 nodes on panels of 0.4, l = 338.6)
changes each of 11 integrals (E4, E2W, E5, E1_7, E3_7 and the two reduced forms) by less than 9.7e-61 relatively. H3: G_k by quadrature of 2(h_k − Ω_{k−1})/S equals (L_k + 2Ω_{k−1})/u at the cutoff to 91 digits for k = 2 … 6.

**Each by-parts step is exact on [0, θ_e].** With the boundary term b of §3 at the cutoff θ_e, "lhs − rhs = b" holds to 91 to 94 digits for all 19 rows below (the right column); the raw difference is the boundary term, which is O(x l^q) ~ 1e-57 at l_e = 342.
All Ω's at their closed forms, so G₂ … G₇ ≈ 0; "scale" is the integral of the absolute values of the larger side.

| id | identity | lhs | raw digits | digits with b | b at the cutoff | scale |
|---|---|---|---|---|---|---|
| E1_6 | E1 at k = 6 (lhs = G₆) | −1.71616079519e-59 | 64.7 | 93.2 | 6.505e-60 | 3.38e+5 |
| E1_7 | E1 at k = 7 (lhs = G₇) | 6.88753635796e-57 | 63.2 | 93.6 | 2.978e-57 | 4.33e+6 |
| E2_2 … E2_6 | E2 at m = 2 … 6 | 29.27, −346.76, −2683.46, 37488.7, 736689 | 68.3, 66.6, 65.5, 63.7, 62.6 | 92.8, 91.4, 91.6, 91.6, 92.2 | 3.9e-67 … −2.9e-57 | 80 … 1.3e+6 |
| E3_3 … E3_7 | E3 at m = 3 … 7 | 118.59, 4251.11, −17646.1, −904683, −1430053 | 67.0, 65.1, 63.7, 62.6, 60.7 | 91.3, 91.8, 91.6, 93.0, 92.5 | −3.4e-65 … 1.0e-54 | 3e+2 … 5.4e+6 |
| **E4_6** | ∫Y₄w₄ = −2∫(…)F₃ | −232044.844217 | 63.5 | 91.7 | 7.241e-59 | 2.45e+5 |
| E4_7 | the same on Y₅, F₄ | 1927597.49904 | 62.3 | 91.9 | 9.84e-57 | 1.99e+6 |
| **E2Wa** | ∫12ℓ²L₄ = −2∫(−4ℓ³)h₄/S | −1010550.94245 | 62.6 | 92.1 | −2.848e-57 | 1.05e+6 |
| **E2Wb** | ∫g₄L₄ = −2∫P h₄/S | 2970980.70448 | 62.6 | 92.1 | 7.635e-57 | 3.03e+6 |
| **E5** | ∫Y₄w₄′ = −2∫W u F₃ | 1553061.99865 | 62.7 | 91.9 | −4.107e-57 | 1.89e+6 |
| G6R | G₆ (definition) = reduced form | −1.71616079519e-59 | 64.5 | 91.6 | −9.327e-60 | 3.18e+5 |
| G7R | G₇ (definition) = reduced form | 6.88753635796e-57 | 62.8 | 92.1 | 7.02e-57 | 4.05e+6 |

**Constants off their roots** (`shifts` part). Ω₄ + 1/10: G₆ = −1.25663706144 = −4π/10 (SLOPE_G6), and E1_6, E2_4, E2_5, E3_5, E3_6, E4_6, E4_7, E2Wa, E2Wb, E5, G6R hold to 62 to 65 digits with G₆ ≠ 0.
Ω₅ + 1/10: G₇ = −4π/10 (SLOPE_G7), all rows hold, G7R included (62.8 digits). So the reduced forms carry the slope −4π.
Ω₂ + 1e-3 (HYP): G₄ = −0.0125663706144 = −4π·1e-3. **E2Wb fails, by 0.73 digits**, and equals its boundary term at the cutoff to 58.9 digits (the growing solution G₄ tan θ); E2Wa also (0.70), E2_4, E2_5, E2_6, E3_5 … E3_7, E4_7, G6R, G7R fail;
**E4_6 holds to 63.5 digits and E5 to 62.7**, as §3 says (they use G₂ = G₃ = 0 only). With a table that has no root at all (`c_wrong`, §9) the same rows fail.

**G_k = 0 at the closed forms.** G₂ … G₇ by quadrature: 2.3e-69, 7.6e-67, −1.4e-64, −1.2e-61, −1.7e-59, 6.9e-57. At the second cutoff (l = 338.6, scratch `scratch_tail.py`) they are 6.2, 6.1, 6.1, 6.0, 5.9, 5.9 times larger,
against e^{(342.25 − 338.56)/2} = 6.3: they are the tails x l^q of the cutoff, not an error of Ω. So Ω₀ … Ω₅ are roots to better than |G_k|/(4π) ≈ 1e-60 … 6e-58 (an upper bound that the tail sets, from a hierarchy that uses
no closed form of Y or L and no part of the finite-part reduction).

## 7. Comparison with `omega4.md` and `omega5.md`

Made after the first set of runs was committed (§1). I read the reductions, status paragraphs and defect lists of both notes.

**No difference in a sign, a constant or a result.** The two notes and this one agree on: c₁ and c₂ in the S-form; the Y₄ weight w₄ and its primitive W₄ (up to the wording below); ∫Y₄w₄ = −2∫(22ℓ² − 8Ω₀ℓ − 2Ω₁ + 16uC)F₃; the L₄ collapse to 4Ω₀∫L₄;
g₄ = −32ℓ² + 16Ω₀ℓ + 4Ω₁ − 32uC and its split 12ℓ² − 2(22ℓ² − …); ∫u²C = uC + ℓ² + 2 ln 2 ℓ − D₃/2 and P; w₄′ term by term, a₁, a₂, a₃ and −20ℓ² − 8Ω₀ℓ; the primitive W, the cancellation of Z = ∫C cot θ and the vanishing constant;
E4 one level up with the same weight; the reduced forms of G₆ and G₇ (E4_6 … G7R pass with no change of coefficient); the masters ⟨lt²C²⟩ and ⟨lt³D₃⟩ (§8); the finite-part calculus (a, b, c of §8 are
the statements of `omega4.md` step 5).

**Differences in what is stated** (none changes a result):

1. **Boundary terms.** `omega4.md` ("Not derived here: E4's boundary terms at π/2") and `omega5.md` ("E5's and E2W's boundary terms are argued as E3's in #209", and the same in Defects) leave them argued.
   §3 derives all six, with the throat coefficients and the outer rate, and the hierarchy reproduces each as the difference of the two sides to 91 digits.
2. **The constant of W.** `omega5.md` (Defects, `e5_const`): "the constant of W is zero: W(π/2) = 0 identically". It is zero at the roots of Ω₀ and Ω₁; as a function of Ω₁ it is K₀ = −4π(Ω₁ − 2π²/3 − 8 ln²2) (SY_K, control `e5_K`).
   So E5 uses G₃ = 0 a second time, and `e5_const`'s silence is a statement about the roots.
3. **Hypotheses.** E4 needs only G₂ = G₃ = 0 (at k = 6, Ω₂ is free: HYP shows E4_6 at 63.5 digits with Ω₂ + 1e-3). E4 on Y₅ at k = 7 also needs G₄ = 0 (L₄ in F₄ must be bounded); E2W needs G₄ = 0 and not merely a bounded L₄ (G₄ ≠ 0 gives
   L₄ ≈ G₄/x, §2). `omega5.md` step 3 states no hypothesis for E4 one level up.
4. **The θ term of W₄.** `omega4.md`: "the θ terms cancel: 8Ω₀θ + 16(2qθ − C) = −16C". The primitive is cot θ(…) + 16C − (8Ω₀ + 32 ln 2)θ (§4); the bracket there is the negative of that θ and C part, and it vanishes by G₂ = 0.
   The result is the same.
5. **The L cot θ → 0 step.** `omega5.md` (step 4) justifies the boundary term at π/2 by "L cot θ → 0". The product with P also needs P = O(l³), which is true (ℓ³ and uC = O(l)).
6. **The last Defects line of `omega5.md`** ("none has been re-derived independently (#264 does that for E4)") is superseded by this note for E4, E2W and E5.

Nothing here disagrees with the repository's text or with another computation.

## 8. The masters and the finite parts

**Finite parts.** Let FP∫f be the constant term (x⁰ ln⁰x) of ∫₀^{π/2−x}f as x → 0, for f with an expansion in x^a ln^b x at π/2. (a) FP is linear and equals ∫ for integrable f. (b) A convergent integrand written as a sum of divergent pieces:
∫ = ΣFP∫ pieces, since at each cutoff the integral of the sum is the sum of the integrals and the constant term of a finite sum is the sum of constant terms. (c) By parts, exactly at the cutoff:
FP∫fg′ = FP[fg](π/2) − FP[fg](0) − FP∫f′g. (d) The generators ℓ, C = Cl₂(2x) = 2x(1 − ln 2x) + O(x³), D₃ = (7/4)ζ(3) + O(x² ln x), cot θ, tan θ have expansions in x with polynomial-in-ln x coefficients, so FP is defined.
These are the statements of `omega4.md` step 5 and `omega5.md` step 7 (*proved*). Planted splits: ∫u(π/2 − θ) = (π/2)ln 2 as (π/2)FP∫u − FP∫uθ (FP∫u = 0, FP∫uθ = J₁), and ∫u(J₁ − Lc₁) = −J₂: the rule and the
quadrature agree to 59.1 and 59.2 digits (FP1, FP2). **Not re-derived:** the 11 steps of the k = 6 chain and the 37 integrations by parts of the k = 7 chain. Their output, the closed forms, is tested in §6 (the Ω's are roots).

**Masters, by quadrature of the definitions** (mpmath, tanh–sinh, 58 digits; no series).

- **M1.** ∫₀^{π/2}ℓ²D₃ = 3.98960647090248246096726447258… equals ¾ζ(3)J₂ + (π/4)(2 ln 2 ζ(4) + 5ζ(5) − 2ζ(2)ζ(3)) to 72.2 digits (bar 40). PSLQ on the value finds
  M1/π = (5/4)ζ(5) − (5/96)π²ζ(3) + (π⁴ ln 2)/180 + (3/8)ζ(3)ln²2, the same number: ¾ζ(3)J₂/π gives (3/8)ζ(3)ln²2 + π²ζ(3)/32 and the rest gives −π²ζ(3)/12 + …, which sum to −5π²ζ(3)/96.
- **M2, M3.** ⟨lt²C²⟩ = 0.24177653206372068633365757411286… and ⟨lt³D₃⟩ = −6.2357919357623579752559660515… were found by PSLQ at 70 digits in my own bases
  (π⁶, ζ(3)², π²ζ(3), π⁴, ζ(5), ζ(3), π², 1) with the relations (−18144, −19, 13608) and (−15120, −44, 8505) on (M/π, π⁶, ζ(3)²), zero elsewhere:
  M2 = π(¾ζ(3)² − 19π⁶/18144), M3 = π(9ζ(3)²/16 − 11π⁶/3780). A quadrature at 125 digits meets them to 125.8 and 119.3 digits (M2, M3). With ζ(6) = π⁶/945 these are the forms of `omega5.md`,
  π(¾ζ(3)² − (95/96)ζ(6)) and π((9/16)ζ(3)² − (11/4)ζ(6)), which the quadrature meets to 57.7 and 53.7 digits (NM2, NM3). *Identified* here (no Fourier series, no Parseval); *derived* there through seven Euler sums.
- **The others** (`masters2`: 24 checks, of which NM2 and NM3 are the two averages against the notes' forms; the other 22 all pass): J₁ … J₅ against n![aⁿ]B(a) and J₅ against `omega4.md`'s closed form (51.0 digits); J̃₂ … J̃₆ (49.3 to 55.4); ⟨θ²⟩, ⟨C²⟩, ⟨ℓC²⟩, ⟨lt C²⟩ = −1e-65 (closed form 0), ⟨θ lt C⟩, ⟨θC⟩ (equal at 58 digits);
  ⟨uC⟩ (57.6), ⟨lt²D₃⟩ (55.3), ⟨C Lc₂⟩ (58.5), ⟨CK₂⟩ (58.8), ⟨lt C K₂⟩ (equal), ⟨CK₃⟩ (57.7), ⟨θ²lt²⟩ (55.3), ⟨θK₂⟩ (58.6), ⟨K₂²⟩ (58.5), with K₂ = ∫₀^θlt², K₃ = ∫₀^θlt³.

## 9. Controls

Every check has a `--mutate=` control (18). Their targets were computed in development runs (raw `dev_*`), declared in the raw README (A1, A2) and in `TARGETS`, then run (`ctl2_*`, three at a time). The summary
(`CT_*` in the output) reads the logs: **18 of 18 on target**, meaning the set of checks that fire equals the declared set and the exit code is 1. (A first set `ctl_*` of 15 on the earlier script commit is also 15 of 15.)

| control | breaks | fires (all others silent) |
|---|---|---|
| `e4_20`, `e4_C` | E4's 22ℓ² as 20ℓ² (the task's example); 16uC as 8uC | SY_W4, E4_6, E4_7, G6R, G7R |
| `e4_Om1` | −2Ω₁ as +2Ω₁ | the same and SY_bnd0 |
| `e2w_P` | P without D₃/2 | SY_P, SY_bnd0, E2Wb, G7R |
| `e2w_g` | g₄'s −32ℓ² as −30ℓ² | SY_P, E2Wb |
| `e5_noD3`, `e5_a1` | W without 16D₃ cot θ; a₁'s 64 ln 2 as 32 ln 2 | SY_W, E5, G7R |
| `e5_K` | W + 1/100 | SY_K, E5, G7R |
| `hyp_G4` | Ω₂ + 1e-3 with every check applied | E1_6, E1_7, E2_4 … E2_6, E3_5 … E3_7, E4_7, E2Wa, E2Wb, G6R, G7R, H2_G5, H2_G6 |
| `c_wrong` | c₃'s −2ℓ³ as −3ℓ³ | H2_G3 … H2_G6, E1_6, E1_7, E2_3 … E2_6, E3_4 … E3_7, E4_6, E4_7, E2Wa, E2Wb, E5, G6R, G7R, BND_outer, H1 (silent: H0, BND_throat, E2_2, E3_3) |
| `g6r_slope`, `g7r_slope` | the reduced form without its Ω₄ (Ω₅) term | G6R, G7R |
| `m1_basis`, `m1_J2` | 5ζ(5) as 4ζ(5); J₂ + 0.01 | M1 |
| `fp_const` | FP∫u = −ln 2 | FP1, FP2 |
| `m2_K2K2`, `m2_J5`, `m2_M3` | (61/6)ζ(6) as 10; 1080ζ(5) as 1008; (11/4)ζ(6) as (10/4) | MK2K2; MJ5; NM3 |

Silences that are results: `e2w_g` leaves G7R silent (the reduced form uses P, not g₄). `e5_K` leaves SY_W silent (a constant has no derivative); SY_K sees it. `c_wrong` leaves H0 silent (c₃ first enters h₃; the README's prediction of the control, H0, was wrong:
raw README A1). With `c_wrong`, E2_2 and E3_3 hold (they need only G₀, G₁), and E4_6 fails because G₃ ≠ 0 makes L₃ unbounded. `hyp_G4` shows the hypotheses are needed: E4_6 and E5 stay silent, as §3 predicts.

## 10. What Paper III may say

- E4, E2W and E5 are **proved**, with their boundary terms and hypotheses, in an independent derivation, and checked on an independent hierarchy to 62 to 63 digits (91 to 92 with the boundary term). The sentence "boundary terms argued as E3's"
  can be replaced by "derived in #264". The hypotheses G₂ = G₃ = 0 (E4, E5's constant) and G₄ = 0 (E2W, E4 on Y₅) can be stated.
- Ω₄ and Ω₅ are roots of G₆ and G₇ to better than 1e-58 on a hierarchy that uses neither the identities nor the closed forms of Y₀ … Y₃, L₀ … L₃ nor the finite-part reduction. The reduced integrands equal G₆ and G₇ to 62 digits.
  So the wording "derived" for Ω₄ and Ω₅ is supported at the level of the identities and of the end result. It does **not** rest, in this review, on a second derivation of the finite-part chain (11 and 37 steps) or of the seven Euler sums;
  those are supported by the 22 quadratures of §8, the planted FP splits and the end-to-end row, not re-derived.
- ∫ℓ²D₃, ⟨lt²C²⟩, ⟨lt³D₃⟩ are confirmed by a route that uses neither Fourier series nor Parseval (quadrature, then PSLQ and a 125-digit check). That is evidence for the values, not a second derivation of the Parseval and Euler-sum steps.
- Not to say: that E1 … E5 were machine-proved (the growth bound is a hand induction; the hierarchy checks it); that the master equation of README was re-derived (this review takes it as given, as the notes do).

## 11. Limits

1. The hierarchy uses README's master equation (c_m from its ε expansion), the same equation as the notes. The identities hold for any table that has its own roots (§9, `c_wrong`), so this review tests the identities and the closed forms,
   not the equation. #195's crate and the notes' E0 are its checks.
2. The boundary terms and the hypotheses are proved by hand (§2, §3). The growth bound is an induction on j under H_k, valid to j = 6 under H₇; GROW and BND_outer are numerical consistency checks of it, not a proof.
3. "Digits" are differences between two evaluations on one hierarchy, limited by the cutoff (e^{−171}·l^q ≈ 1e-57) for the raw rows and by the quadrature (about 1e-92) for the rows with the boundary term. They are not error bars.
4. The statement lines of §1 carried sketches of E2W's key equation and E5's primitive. Their independence is partial; the weights, the chain, the boundary terms, the constant of W and the hypotheses are mine.
5. M2 and M3 by PSLQ are identifications; their agreement with the notes' Euler-sum forms is evidence that those derivations are right, not a proof of either.

## Files

- `identity_review.py`, `identity_review_output.txt`: parts `sym exch shifts bnd masters masters2 fp ctlsum` and the control summary (151 checks and 18 controls; 70 s on rpi5 with the node caches, about 4 min without;
  peak 202 MB). With mpmath and sympy installed (the project's venv), run `python3 identity_review.py sym exch shifts bnd masters masters2 fp ctlsum --ctldir=../../../negative-mass-raw/I264-rpi5 --ctlprefix=ctl2` from this directory (the raw clone beside the checkout, its controls run by `ctl2_queue.sh`);
  `--cache=<dir>` keeps the derived nodes between runs; `--mutate=NAME` for the controls.
- Raw `I264-rpi5` (negative-mass-raw): `derivations.md` (written first), `README.md` (criteria, routes, amendments A1 and A2), `dev_*`, `run_01` … `run_03`, `ctl_*` (first set), `ctl2_*` (second set), `scratch_tail.py`.

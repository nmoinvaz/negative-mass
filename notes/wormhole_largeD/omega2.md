# #26: Ω₂ from the order-4 matching integral

Status: Ω₂ is **computed** from G₄ = 0 by four evaluations (two methods, two precisions, two grids) as
74.36938672793331759538467029584661726673414341502…, and the same integral is reduced **by hand** to the closed form Ω₂
= 80ζ(3) − (8/3)π² ln 2 − (32/3) ln³2, so the identification of #195 now has a **derivation** (each step checked
symbolically or in mpmath at 30 digits, checks H0 to H8; it takes Ω₀ = −4 ln 2 from README.md and uses the exchange
identities E1, E2, which are derived here and checked numerically). The same steps at order 3 give Ω₁ = 2π²/3 + 8 ln²2
with no Clausen function at all (a by-product for #25, which is not closed by this note). Every number is in
`order4_output.txt` (the five runs of record, raw folder `negative-mass-raw/I26-dh610`, README first) unless it names
another raw file. No sentence here claims novelty; the prior equations used are README.md's order equations (matched.py,
order3.py), Clausen's functions, and the Beta-function moments of ln cos.

## Result

| evaluation | method | precision | Ω₂ agrees with A70 | Ω₂ against 80ζ(3) − (8/3)π² ln 2 − (32/3) ln³2 |
|---|---|---|---|---|
| A50 | closed forms, one integral | 50 digits | 77.6 digits | 77.6 digits |
| A70 | same | 70 digits | (reference) | 80.0 digits, the cap of the 80-digit comparison |
| B1 | the hierarchy from its definitions, N = 40, w = ½, V = 18 | 120 digits | 48.2 digits | 48.2 digits |
| B2 | same, N = 50, w = 2/5, V = 20 | 140 digits | 64.4 digits | 64.4 digits |

- **Ω₂** = 74.36938672793331759538467029584661726673414341502 (A70). Against #195's fit 74.3693867279333175953847 the
  difference is 3.0e-23 (half a unit of its last digit); against its 38-digit median
  74.369386727933317595384670295846617267 it is 2.7e-37, which is the rounding of that string to 36 decimals (half a
  unit, 5e-37); to those digits the two are equal. The difference between evaluations is a spread over method, precision
  and grid, an estimator range and **not a rigorous bound** (nothing here is a rigorous enclosure). The issue's bar, 20
  digits between two grids, is met 48.2 digits (B1 against B2) and, between the two methods, 64.4.
- **The planted order 3 by the same code.** Ω₀ = −4 ln 2 and Ω₁ = 2π²/3 + 8 ln²2 are reproduced to 80 digits by A
  (order3.py reached 10) and to 48.7/48.2 (B1) and 65.2/64.7 (B2); B finds both from G₂ = 0, G₃ = 0 with no input. The
  slope of G_k in Ω_{k−2} is −4π to 79 to 80 digits in A, 49 (B1) and 65 (B2) digits. With Ω₀ + 1/100 in the orders ≥ 3,
  Ω₁ moves by 0.0276 in both methods (control `wrong_Om0`, as #25 asks; development run `dev_07_wrong_Om0`, P3).
- **PSLQ** (identification) on {Ω₂, ζ(3), π² ln 2, ln³2} at tolerance 1e-18, coefficients ≤ 1000, Ω₂ at 60 digits: the
  relation (3, −240, 8, 32), residual 2.9e-59; the perturbed basis (ζ(3) + 10⁻³) finds none. On the eight-element weight
  ≤ 3 basis {Ω₂, ζ(3), π² ln 2, ln³2, π², ln²2, ln 2, 1} at 1e-30: the same relation with zeros on the four lower terms,
  residual 2.9e-59; the perturbed basis finds none. The relation is (3Ω₂ − 240ζ(3) + 8π² ln 2 + 32 ln³ 2 = 0): the same
  as #195's, now found on the integral's own Ω₂ to 60 digits instead of #195's 45.
- **Remainder** against #195's W² (the qd values of `w2_output.txt`, 21 digits), rem = W²a² − [s − 4 ln 2 + Ω₁/s +
  Ω₂/s²], s = 2(D − 3), with this note's Ω₁ and Ω₂:

  | D | s | rem | rem·s³ | rem·s³ − Ω₃ |
  |---|---|---|---|---|
  | 200 | 394 | −1.11994e-5 | −684.9893 | −32.149 |
  | 500 | 994 | −6.77794e-7 | −665.6667 | −12.827 |
  | 1000 | 1994 | −8.31519e-8 | −659.2459 | −6.406 |
  | 2000 | 3994 | −1.02969e-8 | −656.0409 | −3.201 |

  rem·s³ approaches #195's Ω₃ = −652.84010757006, and the distance to it falls like 1/s (32.1, 12.8, 6.4, 3.2; its
  coefficient is Ω₄ = −12794.5). A two-point extrapolation from D = 1000 and 2000,
  (s₂x₂ − s₁x₁)/(s₂ − s₁), gives −652.8456, 0.0055 from Ω₃. The 60-digit qd values at D = 200 … 2000 (raw
  I195-dh610/run_02) differ from the 21-digit ones by 7e-12 to 2.8e-7 in rem·s³.
  Over the 64 values of D (50 to 10⁴) the number (rem·s³ − Ω₃)·s goes from −11930.7 at D = 50 to −12792.4 at D = 10⁴,
  tending to Ω₄ = −12794.5. **The remainder is Ω₃/s³ with Ω₃ = −652.84 as #195 fitted:
  the four-term series with this note's Ω₂ leaves the D⁻³ remainder the issue asked for.**

## The order equations and the closed forms

The order equations are README.md's (check E0: the ε-expansion of the master equation, sympy, residual 0, k = 0 … 4):

    dY_k/dl = F_{k−1} = Σ_{j<k} L_{k−1−j} l^j/j!,    (R − 1)L_k′ − L_k/2 = h_k = (1 − 1/R)L_{k−1} + Σ_{j<k} Y_j c_{k−j} + Y_k/2,
    L_k = −2Ω_{k−1} + (R − 1)^{1/2} C_k,    G_k = C_k(∞) = ∫₁^∞ (h_k − Ω_{k−1})(R − 1)^{−3/2}dR = 0,

with θ = arctan √(R − 1) (so R = sec²θ, l = −2ℓ, ℓ = ln cos θ, G_k = 2∫₀^{π/2}(h_k − Ω_{k−1})/sin²θ dθ). G_k is affine
in Ω_{k−2} with slope −4π. c_1 … c_4 are in the raw README. The closed forms, with u = tan θ, K = 2θ ln 2 − Cl₂(π − 2θ),
T = ∫₀^θ φ ln cos φ dφ = (θ/2)Cl₂(π − 2θ) − ¼Cl₃(π − 2θ) − (3/16)ζ(3) − (ln 2)θ²/2:

    Y₁ = −l, Y₂ = −2Ω₀ l − (3/2)l², L₀ = −1, L₁ = −2Ω₀ − 2l, L₂ = −2Ω₁ + 2Ω₀ l + 4l² − 8Ω₀ uθ − 16uK     (README.md, checked S2, N3g, N3)
    Y₃ = −2Ω₁ l + l³/2 + 8l² + 8Ω₀ l + 8Ω₀θ² − 16Ω₀ uθ − 32uK + 32Kθ + 64T                               (new; S1, N2)

Y₃ comes from ∫uθ dl = 2uθ − l − θ², ∫uK dl = 2uK − l²/2 − 2Kθ − 4T. At Ω₀ = −4 ln 2, with C = Cl₂(π − 2θ) and D₃ =
Cl₃(π − 2θ) + ¾ζ(3) (so that K = −(Ω₀/2)θ − C, T = θC/2 − D₃/4 − (ln 2/2)θ²), the θ's cancel and (H0, sympy residual 0)

    L₂ = −2Ω₁ − 4Ω₀ℓ + 16ℓ² + 16uC,      Y₃ = 4Ω₁ℓ − 4ℓ³ + 32ℓ² − 16Ω₀ℓ + 32uC − 16D₃.

**L₃.** L₃ = −2h₃ + 2uB(θ), B = ∫₀^θ h₃,φ cot φ dφ (S3, residual 0, any h₃). With the forms above B′ = −16(3 − 2c²)ℓ² −
32 ln 2 (5 − 2c²)ℓ + 32u(1 + c²)C − 64 ln²2 − 2Ω₁(1 + 2c²), c² = cos²θ, no D₃ and no 1/u, and (S4, sympy residual 0; N5,
45 digits at θ = 0.4, 1.0, 1.4)

    B = −32(ℓ + ln 2)C + 16 sin²θ C + 32 Lc2 + 16 sin θ cos θ ℓ(ℓ + 2 ln 2) − 4Ω₁θ − 2Ω₁ sin θ cos θ,    Lc2(θ) = ∫₀^θ ln² cos φ dφ.

With t = e^{2iθ} the log-cosine integral is Lc2 = (1/2i)[F(t) − F(1)], F(t) = ln t ln²(1 + t) − (2/3)ln³(1 + t) − 2 ln(1
+ t)Li₂(1/(1 + t)) − 2Li₃(1/(1 + t)) + ln²2 ln t + (ln t)³/12 + 2 ln 2 Li₂(−t) + ln t Li₂(−t) − Li₃(−t) + (ln 2/2)ln² t
(from ln cos φ = ln(1 + t) − ln 2 − ½ ln t and the standard antiderivatives; N6: 45 digits against quadrature, and L₃
with it against its integral 32.6 digits). So **Y₃ is closed in Cl₂ and Cl₃, and L₃ in Cl₂, Cl₃ and Li₂, Li₃ of 1/(1 +
e^{2iθ}) and −e^{2iθ}** (through Lc2). Y₄ = ∫(L₃ + L₂l + L₁l²/2 + L₀l³/6)dl is not given in closed form: G₄ needs only
L₃ and Y₃ (below), and method B computes Y₄ numerically. **How far the closed forms got:** Y₁ … Y₃ and L₀ … L₃ closed
(L₃ with Li₂, Li₃); Y₄ and h₄ only as integrals; G₄ closed (below).

## The two evaluations of G₄

- **Method A.** G_k as one integral of closed forms: with (E1) G_k = 4∫₀^{π/2}L_{k−1}dθ + 2∫Σ_{j=1}^{k−1}L_{k−1−j}l^j/j!
  dθ + 2∫(Σ_{j<k}Y_j c_{k−j} − Ω_{k−1})/sin²θ dθ (from sin²θ = 1 − 1/R and ∫Y_k/sin²θ dθ = 2∫F_{k−1}dθ) and (E2)
  ∫₀^{π/2}L_{k−1}dθ = 2∫h_{k−1} ln cos θ/sin²θ dθ (from L = −2h + 2 tan θ ∫₀^θ h_φ cot φ dφ and G_{k−1} = 0). E2 is
  checked at k − 1 = 2 to 44.5 digits (N4a) and at k − 1 = 3 by nested quadrature: N4bg, 30 guard digits, 27.5 digits
  (bar 25, met); N4b as written, no guard digits, 22.0 digits (not met). Ω₀, Ω₁, Ω₂ are solved in turn by the same code
  from G₂, G₃, G₄ (Gauss–Legendre on [0, 1], tanh-sinh in π/2 − θ beyond, 30 guard digits).
- **Method B.** The hierarchy from its definitions: Chebyshev panels in v = √l, cumulative integrals for Y_k, the tail
  integral for L_k (the analytic solution when G_k = 0, which avoids the u ~ e^{v²/2} amplification of a left-hand
  integral's error), Ω_k from G_{k+2} = 0. No Cl₂, Cl₃, E1, E2 and no input value of Ω₀ or Ω₁.
- **Shared and not shared.** Both use the table of c_m and h_k (one sympy expansion, E0) and the analytic-at-throat
  condition. A and B share no closed form, no integration by parts and no grid. #195's fit (the full equation in a Rust
  crate, no matching) shares none of the matched machinery. That the A evaluations agree with B to 64 digits is the
  check of E1 and E2 at order 4 that does not use them.

## The reduction of G₄ by hand

Notation: S = sin²θ, ℓ = ln cos θ, u = tan θ, C = Cl₂(π − 2θ), D₃ = Cl₃(π − 2θ) + ¾ζ(3), q = ln 2, Ω₀ = −4q, J_n =
∫₀^{π/2}ℓ^n dθ. Each line is a check of `order4.py reduce` (H0 … H8).

    G₄ = ∫₀^{π/2} I₄ dθ,    I₄ = A(ℓ)/S + B(ℓ) + γ C + δ D₃                                             (H1: against method A's integrand, 45 digits)
    A  = (398/3)ℓ⁴ − (392Ω₀/3 + 64)ℓ³ + (16Ω₀² + 96Ω₀ − 24Ω₁)ℓ² + (16Ω₀Ω₁ − 32Ω₀²)ℓ + 8Ω₂ℓ,   B = 104ℓ³ − (72Ω₀ + 64)ℓ² + (32Ω₀ − 16Ω₁)ℓ
    δ = −32(Ω₀ − S − ℓ)/S,    γ = 64u(Ω₀ − S − ℓ + Sℓ)/S

1. **The D₃ term, by parts.** With P = cot θ (ℓ − Ω₀), P′ = (Ω₀ − S − ℓ)/S, so δ = −32P′ (H2a, residual 0); D₃′ = 2C and
   D₃(0) = 0, P(π/2) = 0, so ∫δD₃ = 64∫C P dθ (H2b, 45 digits; the boundary terms vanish).
2. **The C terms combine.** γ + 64P = 64(Ω₀ − 1)u (H3, residual 0). So the Clausen content of G₄ is 64(Ω₀ − 1)∫₀^{π/2}
   tan θ C dθ.
3. **The Clausen integral.** C′ = 2ℓ + 2q (checked to 45 digits), ℓC → 0 at both ends, so ∫tan θ C dθ = ∫ℓ C′ dθ = 2J₂ +
   2qJ₁ = π³/12 (H4, 45 digits against the quadrature of its definition). The C-part of G₄ is −(16π³/3)(4 ln 2 + 1).
4. **The logarithmic part.** ∫ℓ^b/S dθ = −bJ_{b−1} (b ≥ 1, by parts, H5, 43 digits) and A(0) = 0, so ∫(A/S + B)dθ = −Σ b
   a_b J_{b−1} + Σ b_b J_b.
5. **The moments.** J_n from ∫₀^{π/2}cos^a = (√π/2)Γ((a + 1)/2)/Γ(a/2 + 1) (sympy series, H6, 43 digits): J₀ = π/2, J₁ =
   −(π/2)q, J₂ = (π/2)(q² + π²/12), J₃ = −(π/2)(q³ + (π²/4)q + (3/2)ζ(3)).
6. **Assembly** (H7, sympy, difference 0; H8, against method A at 44 digits), with Ω₁ = 2π²/3 + 8q² from order 3 below
   and Ω₂ = 0 in A:

   | piece | value |
   |---|---|
   | ∫A/S dθ | (π/3)(1194ζ(3) − 596q³ − 96q² + 19π²q + 24π²) = 1631.0111775700… |
   | ∫B dθ | (π/3)(−234ζ(3) + 468q³ + 96q² + 13π²q − 8π²) = −72.59694306292… |
   | Clausen part | −(64π³/3)q − 16π³/3 = −623.86095852130… |
   | G₄(Ω₂ = 0) | π(320ζ(3) − (128/3)q³ − (32/3)π²q) = 934.5532759858143274892985… |

   and G₄ = G₄(Ω₂ = 0) − 4πΩ₂, so **Ω₂ = 80ζ(3) − (8/3)π² ln 2 − (32/3) ln³ 2**. The π³ terms cancel, as does everything
   of weight 5 (there is no ℓ⁴ without 1/S).
7. **Order 3 by the same steps.** G₃ has no C or D₃ at all: L₂ is replaced by E2 in terms of h₂, which has only L₁, Y₁,
   Y₂. With A₃ = −28ℓ³ − 24Ω₀ℓ² + (8Ω₀² + 8Ω₁)ℓ and B₃ = −20ℓ² − 16Ω₀ℓ: G₃ = (4π/3)(−3Ω₀² − 12qΩ₀ − 3Ω₁ + 24q² + 2π²),
   which with Ω₀ = −4q is −4π(Ω₁ − 2π²/3 − 8q²): **Ω₁ = 2π²/3 + 8 ln²2** (H7). With Ω₀ → Ω₀ + 1/100 it moves (control).

**What each status word covers.** *Derived* (by hand, each step checked): steps 1 to 7 given E1, E2, the closed forms
(checked symbolically) and Ω₀ = −4 ln 2 (README.md, from G₂ = 0); the independent confirmation of the value is method B
at 64 digits and #195's fit at 38. *Numerical evidence*: the four evaluations and their agreement. *Identification*: the
PSLQ relations, which the derivation supersedes. What is not derived here: Ω₀ is README.md's derivation, E2's boundary
terms are checked numerically (N4a, N4bg), the Beta-function moments are standard.

## Controls and checks

`order4.py` has 35 checks in the runs of record (E0, S1 … S4, N1 … N6 with N3g, N4a, N4b, N4bg, H0 … H8 with H2a, H2b,
P1, P2, M1 … M3, X1, X2, Q1, Q2, R1, R2; P3 is the control `wrong_Om0`'s own statement and reports only under it). **33
are met, N3 and N4b as written are not.** Each has a control (`--mutate=NAME`) whose targets were computed in
development runs (raw README, amendments A2 to A4). Silent targets are listed there: `wrong_Om1` and the pair (A, B)
take the same wrong input and so do not move M1 to M3; `T_const` and `K_sign` do not enter the planted orders (K, T
enter at k = 4); `short_tail` is shared by both B grids so M2 stays; `coarse_B` fires Q2 only because PSLQ is given the
14 digits shared with B (A70 itself is untouched); `wrong_h`, `bad_exchange` and the others fire as listed.

## Defects and discrepancies

- **N3 and N4b as written are not met** and stay so (`run_01`, `run_02`; the same 25.8 and 22.0 digits as in the
  development runs): N3 (L₂ closed form against its integral) reads 25.8 digits against 30, N4b (E2 at order 3 by nested
  quadrature, dps 30) 22.0 against 25. The cause is the checks' own cancellation (the integrand near the throat, and u·K
  near π/2, with no guard digits): N3 with 30 guard digits reads 40.1 digits (N3g). N3g and N4bg were added by the
  coordinator's ruling (A1) with the same bars.
- A **first B** had a wrong integration matrix (T₀ taken as 1), giving slope −19197 instead of −4π; found in the first
  development run and fixed before any run of record. A PSLQ crash on disagreeing digits (`prec cannot be less than 53`)
  was found in the controls and fixed. The N4b control under `T_const` was stopped after 58 min. An earlier N6 ran its
  quadrature at 75 digits and read 2.4 digits (tanh-sinh nodes at 1e-70 with the integrand's cancellation); fixed by
  guarding the integrand only (raw README A4).
- **#25's text** asks the closed form to involve K(u) and reduce to polylogarithms; with E2 the order-3 integral has no
  K at all, only J₀ … J₂. That reduction is a finding for #25, not its acceptance.
- **growth.csv** is not used (#195: its large-D rows are low by 8e-10 to 6e-6).

## What README.md's "Result" and "Limits", #25 and Paper III should now say

- **README.md, Result.** Replace "Both remain identifications, and the derivation is #25's and #26's" by: Ω₂ = 80ζ(3) −
  (8/3)π² ln 2 − (32/3) ln³2 is **derived** from the order-4 matching integral (this note: G₄ = 0, one integral of
  closed forms, reduced to the Clausen integral ∫ tan θ Cl₂(π − 2θ) = π³/12 and the moments ∫ ln^n cos), and confirmed
  to 64 digits by an evaluation that uses none of the reduction and to 38 digits by #195's fit. The fifth term is Ω₃ =
  −652.84010757006, a fit (#195). The sentence "the three terms leave about 73/(2(D − 3))²" is superseded: with Ω₂ the
  remainder is Ω₃/(2(D − 3))³ (table above).
- **README.md, Limits.** "The third term's closed form is an identification, not a derivation" stays for Ω₁ until #25 is
  accepted (the reduction above is available there) and becomes false for Ω₂. Add: Ω₃ is a fit, not identified in
  weight-4 constants (#195).
- **README.md, Method.** Add order 4: the closed forms Y₃, L₃ (above), the exchange identities E1, E2, and G₄ = 0.
- **#25.** The same reduction at order 3 (no Clausen term) gives Ω₁ = 2π²/3 + 8 ln²2 in a few lines (H7), with each
  integral checked at 43 digits (H5, H6); the control with a wrong Ω₀ moves Ω₁ by 0.0276 (P3). Whether that meets #25's
  acceptance is the coordinator's to decide.
- **Paper III.** The fourth term can be quoted as derived with the same matching integral as the third; its derivation
  needs E1, E2 and one Clausen integral. The remainder statement is the one in the table above. The paper's "how this
  work was done" section should record that the closed forms and their checks were done by a language model with sympy
  and mpmath, as AGENTS.md asks.

## Files

- `order4.py` (parts `sym`, `closed`, `lc2`, `reduce`, `n4b`, `n4bg`, `evals`, `remainder`; `--mutate=NAME`),
  `order4_output.txt` (the five runs of record concatenated), this note.
- Raw: `negative-mass-raw/I26-dh610` (README with the criteria and amendments A1 to A3; `dev_*` the development runs;
  `run_01` … `run_05` the runs of record; `inputs_w2_growthD.txt`; the scratch derivations `dev_12_reduce1` …
  `dev_17_l3ansatz`).

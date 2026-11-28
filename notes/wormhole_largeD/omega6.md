# #282: Ω₆ from the order-8 matching integral, and the ln 2 factor of F at ε⁷

Status: Ω₆ is **computed** from G₈ = 0 by method B on two grids at two precisions, agreeing to 430.8 digits. **The ln 2 structure of #263 holds one order further, as an identification:** F's seventh coefficient f₇ is a rational combination of ln 2-free monomials of π², ζ(3), ζ(5) and ζ(7), found by PSLQ at 390 digits in the ln 2-free basis F6 and, as the same relation with zero coefficients on every ln 2 monomial, in the three bases X1, X2, X3 and in the full basis XF (36 constants, up to ln⁷2). It is *identified* (PSLQ); it is a *conjecture* until a by-hand reduction meets it, and **no by-hand reduction of G₈ is attempted here**. Ω₆'s closed form follows from f₇ by exact algebra, with Ω₀ … Ω₅ taken at their closed forms (Ω₅ is derived in #263; E2W and E5 are under review in #264).
Every number is in `order8_output.txt` (the runs of record, raw folder `negative-mass-raw/I282-rpi5`, README first) unless it names another raw file. No sentence here claims novelty; the prior equations used are README.md's order equations (matched.py, order3.py … order7.py), the closed forms of #26, #208, #220 and #263, and #195's values of W².

## Result

| evaluation | method | N, w, V, working digits | Ω₆ against G1 |
|---|---|---|---|
| G1 | the hierarchy from its definitions, 200 panels | 300, 1/4, 50, 465 | (reference) |
| G2 | same, 164 panels | 310, 3/10, 49, 455 | 430.8 digits |
| B3, B2, B1 (reported) | same | 170, 1/4, 40, 285; 100, 1/4, 32, 170; 60, 2/5, 26, 100 | 251.1, 149.0, 78.1 digits |
| closed form (below) | exact algebra from the F6 relation | evaluated at 400 digits | 400.0 digits on G1 and on G2 (the working precision of the check) |
| X4 | #195's 64 values of W² (qd, 60 digits), polynomial extrapolation in 1/s with Ω₀ … Ω₅ at their closed forms | windows K = 14 … 44 | 1.1e-23 … 4.4e-21 in absolute difference |

- **Ω₆** = 3672259.162594390489657859369620626570276236684090990352394926781482389956927483434947129858452968847371327045151799415381323739262… (G1; G2 agrees to 430.8 digits, |difference| 6.39e-425). The difference between evaluations is a spread over grid and precision, an estimator range, **not a rigorous bound**; nothing here is a rigorous enclosure.
- **In 1/(D − 3).** Ω₆/64 = 57379.0494155373514009040526503222901605661981889217242561707…, and the series is

  ```math
  W^2a^2 = 2(D-3) - 4\ln 2 + \frac{\Omega_1/2}{D-3} + \frac{\Omega_2/4}{(D-3)^2} + \frac{\Omega_3/8}{(D-3)^3} + \frac{\Omega_4/16}{(D-3)^4} + \frac{\Omega_5/32}{(D-3)^5} + \frac{\Omega_6/64}{(D-3)^6} + O((D-3)^{-7}),
  ```
  with Ω₁/2 = 5.2116801893692586, Ω₂/4 = 18.592346681983329, Ω₃/8 = −81.605013446257016, Ω₄/16 = −799.65624125144966, Ω₅/32 = 1276.8892298289277, Ω₆/64 = 57379.049415537351.
- **The planted orders by the same code.** `omega_chain` finds Ω₀ … Ω₆ in turn from G₂ … G₈ = 0 with no closed form of Y_k or L_k and no input value; Ω₀ … Ω₅ against their closed forms (PL0 … PL5): G1 439.6 to 440.2 digits, G2 429.7 to 430.8 (the issue's bar: 100).
  The slope of G_k in Ω_{k−2} is −4π to 441.4 (G1) and 431.2 (G2) digits at every k = 2 … 8 (SL), and |G_k| at the chain's root is below 6.75e-454 (G1) and 3.17e-443 (G2) times |slope| (RES).
- **Rounding and grids.** The digits law of #208 holds at order 8: N = 60, 100, 170 give 78.1, 149.0, 251.1 digits of Ω₆ against G1, and G2 (N = 310) 430.8 (expected 77.1, 148.9, 253.1, 436.9; Ω₀ … Ω₅ of G1 against their closed forms 439.6 … 440.2 digits for the expected 446.6: the grids lose about 6 digits, as at order 7). The working digits are the expected digits plus 18.
  The cached chain has 56 integration units for k = 2 … 8 against 154 uncached, and equals the uncached chain bit for bit (CACHE).
  Run times (rpi5, three processes at nice 10): G1 2 h 14 m 06 s (panel setup 598 s), G2 1 h 55 m 27 s (setup 644 s), B1 … B3 22 m 59 s; peak resident sets 1.21 GB, 1.04 GB, 0.49 GB. The development run of B3 took 1055 s (one process beside two others); the estimate before the grids of record was 2.5 h and 2.3 h.
- **#195's W²** (the independent reference): the refit of its 64 values (qd, 60 digits, D = 50 … 10⁴) with Ω₀ … Ω₅ fixed at the closed forms, by polynomial extrapolation in 1/s through the K largest-s points, K = 14, 18, 22, 26, 30, 36, 44, gives Ω₆ to 1.1e-23, 8.0e-23, 3.1e-22, 7.9e-22, 1.5e-21, 2.8e-21, 4.4e-21 (X4, bar 1e-18, set knowing a 148-digit development value; an estimator range).
  Two methods that share nothing agree to 1.1e-23 or better in absolute difference (3.0e-30 relative at K = 14, 1.2e-27 at K = 44). The remainder of the seven-term series (rem₇ = W²a² − [s − 4 ln 2 + Ω₁/s + … + Ω₅/s⁵], s = 2(D − 3), against #195's qd W²):

  | D | s | rem₇·s⁶ | (rem₇·s⁶ − Ω₆)·s |
  |---|---|---|---|
  | 1000 | 1994 | 3684791.7782838 | 24990035.6846 |
  | 2000 | 3994 | 3678586.1009805 | 25269791.9143 |
  | 5103 | 10200 | 3674753.0692428 | 25437847.8133 |
  | 10000 | 19994 | 3673534.0735349 | 25490569.3438 |

  rem₇·s⁶ approaches Ω₆ (R7: the two-point extrapolation from D = 5103, 10⁴ gives 3672264.5, relative difference 1.47e-6 from the computed value; an estimate), and the last column is rem₈·s⁷ with this Ω₆, which approaches Ω₇ ≈ 2.55 × 10⁷ (two-point extrapolation 25545476; **an estimate, not computed**).
  R8 (the eight-term series is stable between D = 5103 and 10⁴: |rem₈·s⁷(10⁴) − rem₈·s⁷(5103)| = 52722, bar 10⁵) sees an error of Ω₆ of about 5 or more; X4 is the sharp check.

## The closed form of Ω₆ (identified)

With f₇ the coefficient of ε⁷ of F(ε) = e^{4 ln 2 ε}(1 + Ω₀ε + Ω₁ε² + …) (so f₇ = Ω₆ + Σ_{j≤5}Ω_j(4 ln 2)^{6−j}/(6−j)! + (4 ln 2)⁷/7!, with Ω₀ … Ω₅ exact), PSLQ finds, at d = 390 digits (run_06a, run_06b, run_06d), the one relation (up to the sign of the vector)

    f₇ = (7656704/2835)π⁶ + (328672/45)π⁴ζ(3) − (74752/5)π⁴ + 30720π²ζ(3) − 2176π²ζ(5) + 36864ζ(3)² + 565248ζ(3) − 307200ζ(5) + 974336ζ(7)

(the integer vector of F5, F6, X1, X2, X3 and XF is (−2835, 20706336, −6168960, 2762242560, 7656704, 104509440, 87091200, −870912000, −42384384, 1602478080, 0, …, 0) on f₇, π⁴ζ(3), π²ζ(5), ζ(7), π⁶, ζ(3)², π²ζ(3), ζ(5), π⁴, ζ(3), π², 1 and the ln 2 monomials; height 2.76e9), and then

    Ω₆ = (155248/315)π⁶ln 2 + (7656704/2835)π⁶ − (704/135)π⁴ln³2 − (19456/45)π⁴ln²2 − (157696/45)π⁴ln 2 + (328672/45)π⁴ζ(3) − (74752/5)π⁴
         − (256/45)π²ln⁵2 − 1280π²ζ(3)ln²2 + (1024/3)π²ζ(3)ln 2 + 30720π²ζ(3) − 2176π²ζ(5)
         − (1024/315)ln⁷2 + (2560/3)ζ(3)ln⁴2 + 4096ζ(3)ln³2 + 24576ζ(3)ln²2 − 83456ζ(5)ln²2 − 12800ζ(3)²ln 2 + 147456ζ(3)ln 2 − 296960ζ(5)ln 2
         + 36864ζ(3)² + 565248ζ(3) − 307200ζ(5) + 974336ζ(7)

and in 1/(D − 3), Ω₆/64 = (9703/1260)π⁶ln 2 + (119636/2835)π⁶ − (11/135)π⁴ln³2 − (304/45)π⁴ln²2 − (2464/45)π⁴ln 2 + (10271/90)π⁴ζ(3) − (1168/5)π⁴ − (4/45)π²ln⁵2 − 20π²ζ(3)ln²2 + (16/3)π²ζ(3)ln 2 + 480π²ζ(3) − 34π²ζ(5) − (16/315)ln⁷2 + (40/3)ζ(3)ln⁴2 + 64ζ(3)ln³2 + 384ζ(3)ln²2 − 1304ζ(5)ln²2 − 200ζ(3)²ln 2 + 2304ζ(3)ln 2 − 4640ζ(5)ln 2 + 576ζ(3)² + 8832ζ(3) − 4800ζ(5) + 15224ζ(7).
The expression equals the grids' Ω₆ to 400 digits (ID6; its evaluation is exact algebra, so the check is the relation's, at the 410 digits of its verification: 4.0e-401 on G1 and on G2, bar 1e-398).

**What it is not.** It was not derived: the closed form rests on PSLQ and on Ω₀ … Ω₅. Ω₆ has 24 terms against Ω₅'s 15; **the weights of Ω₆ are 7, 6, 5, 4, 3** (the pattern of Ω₅ continues one weight up: weights 6, 5, 4, 3):
all 11 monomials of weight 7 are present (part 2082208.554…), 6 of the 8 of weight 6 (part 2420538.593…; ln⁶2 and π²ln⁴2 are absent), 4 of the 6 of weight 5 (part −176504.056…; ln⁵2 and π²ln³2 absent), 2 of the 4 of weight 4 (part −1333444.188…; π⁴ and ζ(3)ln 2 only, as in Ω₄ and Ω₅) and ζ(3) alone of weight 3 (part 679460.260…);
nothing of weight ≤ 2. The parts partly cancel; the sum of the five parts is Ω₆ (2082208.554 + 2420538.593 − 176504.056 − 1333444.188 + 679460.260 = 3672259.163). The relation has no Li₄(½), Li₅(½), Li₆(½) or Li₇(½) term and needs none; none of them was in a basis (see "What each basis covers").

## The structure of F, one order further

F(ε) = 1 + f₂ε² + f₃ε³ + f₄ε⁴ + f₅ε⁵ + f₆ε⁶ + f₇ε⁷ + O(ε⁸), with (f₂ … f₆ of #263, checked here symbolically from the closed forms of Ω₀ … Ω₅: LN2)

    f₂ = 2π²/3,  f₃ = 80ζ(3),  f₄ = 22π⁴/45 − 384ζ(3),  f₅ = −(2432/45)π⁴ − 160π²ζ(3) + 3072ζ(3) − 10432ζ(5),
    f₆ = −(38812/315)π⁶ + (39424/45)π⁴ − (256/3)π²ζ(3) + 3200ζ(3)² − 36864ζ(3) + 74240ζ(5),
    f₇ = (7656704/2835)π⁶ + (328672/45)π⁴ζ(3) − (74752/5)π⁴ + 30720π²ζ(3) − 2176π²ζ(5) + 36864ζ(3)² + 565248ζ(3) − 307200ζ(5) + 974336ζ(7).

So **W²a² = s·2^{−4/s}F(1/s) with F free of ln 2 through ε⁷** (the first six by symbolic checks, f₇ by PSLQ). f₇ has the weights 7 (π⁴ζ(3), π²ζ(5), ζ(7)), 6 (π⁶, ζ(3)²), 5 (π²ζ(3), ζ(5)), 4 (π⁴), 3 (ζ(3)) and no π² or 1: the weights of f_k are {k} … {3} as before, ζ(7) enters first at k = 7 (as ζ(5) did at k = 5, ζ(3) at k = 3).
**The pure ζ(3) coefficients** of f₄, f₅, f₆ are −384, 3072, −36864 and #263 noted that 96(−4)^{k−3}(k − 3)! gives them for k = 4, 5, 6; **f₇'s is 565248, and the formula gives 589824** (k = 7): the one-parameter pattern fitted to three numbers does not continue (the difference is 24576 = 3·2¹³). It was an observation, not a claim, and it is withdrawn.

| basis | constants besides f₇ | n | H | budget | result |
|---|---|---|---|---|---|
| F1 weight 7 | π⁴ζ(3), π²ζ(5), ζ(7) | 3 | 10¹⁶ | 66.2 | not found |
| F2 weights 7, 6 | F1, π⁶, ζ(3)² | 5 | 10¹⁶ | 98.8 | not found |
| F3 weights 7 … 5 | F2, π²ζ(3), ζ(5) | 7 | 10¹⁶ | 131.4 | not found |
| F4 weights 7 … 4 | F3, π⁴ | 8 | 10¹⁶ | 147.7 | not found |
| F5 weights 7 … 3 | F4, ζ(3) | 9 | 10¹⁶ | 164.0 | **found** |
| F6 weight ≤ 7 | F5, π², 1 | 11 | 10¹⁶ | 196.6 | **found**, the same relation, zero on π², 1 |
| X1 | F6 and ln 2 × 8 monomials of weight ≤ 6 | 19 | 10¹² | 247.0 | **found**, the same relation, zero on the 8 ln 2 monomials |
| X2 | X1 and ln²2 × 6 of weight ≤ 5 | 25 | 10¹¹ | 294.8 | **found**, the same, zero on the 14 |
| X3 | X2 and ln³2 × 4 of weight ≤ 4 | 29 | 10¹⁰ | 310.0 | **found**, the same, zero on the 18 |
| XF full, weight ≤ 7 (run_06c, code 0ac9b7e) | X3 and ln⁴2 × 3, ln⁵2 × 2, ln⁶2, ln⁷2 | 36 | 10⁸ | 308.1 | not found: the relation's height 2.76e9 is above H |
| **XF full, weight ≤ 7 (run_06d, the declared H)** | the same | 36 | 10¹⁰ | 382.1 | **found**, the same relation, zero on all 25 ln 2 monomials |

PSLQ at d = 390 digits (`mp.pslq`, tolerance 1e-378, maxsteps 10⁶, f₇ of G2; the two grids share 430.8 digits, so every relation is also verified at 410 digits on both grids), bases declared before the run (raw README), each with a planted relation recovered exactly and a perturbed basis (a constant + √2·10⁻³) that finds nothing; a found relation counts when its budget (n + 1)log₁₀(2H + 1) + 1 is below d and the perturbed basis finds nothing.
The relation is the same in every basis that found it, residual 1.77e-381 (bar 1e-378), holds at 410 digits to 4.0e-401 on G1 and G2 (bar 1e-398).
**The README's H for X1, X2, X3 and XF were typed from memory and differ from the code's (X1 10¹², X2 10¹¹, X3 10¹⁰, XF 10⁸ at first):** the runs of X1, X2, X3 use the code's (found; the height 2.76e9 is below all of them), XF ran twice, at the code's H (not found) and at the declared (found); see the raw README, A2 and A3.

**Verdict on the ln 2 factor at ε⁷.** *Identified* (PSLQ, 390 digits, ln 2-free basis F6 and the full basis XF with ln 2 to ln⁷2 returning the same relation): F is free of ln 2 at ε⁷. It is a *conjecture* for the series beyond, and it is not derived: the by-hand reduction of G₈ would have to meet it.

**What each basis covers.** The ln 2-free monomials of weight ≤ 7 span the non-alternating Euler sums of weight ≤ 7 (Zagier's dimensions 1, 0, 1, 1, 1, 2, 2, 3; independence is conjectural). Odd powers of π, Catalan's constant, Li₄(½) and its products with ln 2 and π², and the alternating sums of weight 6 and 7 beyond them are in no basis: **a "not found" in F1 … F4 and XF (H = 10⁸) is for these bases only**, and the "found" does not need them.

## The order equations at order 8

README.md's equations (check E0: the ε-expansion of the master equation, sympy, residual 0, now through k = 8):

    dY_k/dl = F_{k−1} = Σ_{j<k} L_{k−1−j} l^j/j!,    (R − 1)L_k′ − L_k/2 = h_k = (1 − 1/R)L_{k−1} + Σ_{j<k} Y_j c_{k−j} + Y_k/2,
    L_k = −2Ω_{k−1} + (R − 1)^{1/2}C_k,    G_k = C_k(∞) = ∫₁^∞ (h_k − Ω_{k−1})(R − 1)^{−3/2}dR = 0,

with c₁ … c₇ as in `order7.py` and the new

    c₈ = Ω₀l⁷/5040 + Ω₁l⁶/720 + Ω₂l⁵/120 + Ω₃l⁴/24 + Ω₄l³/6 + Ω₅l²/2 + Ω₆l + Ω₇ + l⁸/80640 + l⁷/5040 − l⁷/(5040R),

the coefficient of ε⁸ of (1 + εΩ)e^{εl} − (½ + ε)(1 − 1/R)e^{−εl} − e^{−εl}/(2R). G_k is affine in Ω_{k−2} with slope −4π; Ω₇ drops out of G₈ (it is h₈(1), subtracted), so **G₈ = 0 gives Ω₆**. Method B (order4.py's `PanelB`, shared) and the cached chain are as in omega5.md.
**The by-hand reduction of G₈ is not attempted.** It would need E1 … E3, E4 two levels up, E2W and E5 at the next level and new steps for L₅ and Y₅, and the finite-part integration by parts to masters of weight 7 (J̃₇, new Euler sums); none of that was done, and nothing here says how far it would get.

## What each status word covers

*Computed*: Ω₆ on two grids (G1, G2) that use none of a reduction and share only the code, the table of c_m (E0), the analytic-at-throat condition and `PanelB`; the independent check is #195's W² refit (X4).
*Numerical evidence*: the agreement of G1, G2, the closed form and #195's fit. *Checked symbolically*: E0 (k ≤ 8) and LN2 (f₂ … f₆ ln 2-free). *Identified* (PSLQ): f₇ and, with Ω₀ … Ω₅ at their closed forms, Ω₆. *Derived*: nothing at order 8.
The closed form of Ω₆ inherits what is under review for Ω₅ (E2W, E5: #264) only through the planted values (the chain computes Ω₅ itself, equal to the closed form to 429.7 … 440.2 digits).

## Defects and discrepancies

- **Three H of the raw README's PSLQ table were typed from memory** (X1, X2, X3: 10¹⁴, 10¹³, 10¹² against the code's 10¹², 10¹¹, 10¹⁰; XF 10¹⁰ against 10⁸ in the code of the first run). The consequences are stated above and in A2, A3; XF at the declared H is `run_06d`; every found relation has a budget below d and a perturbed-basis control. The H are not a result.
- **The bars of R7, R8 and X4 were set knowing a development Ω₆** (148 digits): 5e-5, 10⁵ and 1e-18 (declared in the raw README). R7 does not see an error of Ω₆ below about 180, R8 below about 5; X4 does not see an error below 1e-18 absolute.
- **First versions that failed or were stopped:** `dev_08` (the planted PSLQ tool checks at cap 200) was stopped by me after 11 minutes inside the full basis (`dev_08_pslq_planted.stopped`); `run_06c` is the XF run at 10⁸ (not found, by height).
- **Controls that were silent where a firing is natural** (declared results, raw README A1): `wrong_c8` and `wrong_Om5` fire only X4: the change in Ω₆ is below R7's and R8's sensitivity (as `wrong_c7` and `wrong_Om4` at order 7).
- **A wrong number in a post:** the first progress comment on the issue gave the resident sets as "about 1.0 and 0.9 GB" before they were measured (0.58 and 0.56 GB at the time); it was corrected by comment id.
- **Shared machinery.** G1 and G2 use the same `PanelK` and table; the independent check of Ω₆ is #195's refit (X4, 1.1e-23). The closed form is PSLQ's, not an independent evaluation.

## What README.md's "Result" and "Limits", #195, #210 and Paper III should now say

- **README.md, Result.** Add the eighth term of the series: Ω₆/64 = (9703/1260)π⁶ln 2 + (119636/2835)π⁶ − (11/135)π⁴ln³2 − (304/45)π⁴ln²2 − (2464/45)π⁴ln 2 + (10271/90)π⁴ζ(3) − (1168/5)π⁴ − (4/45)π²ln⁵2 − 20π²ζ(3)ln²2 + (16/3)π²ζ(3)ln 2 + 480π²ζ(3) − 34π²ζ(5) − (16/315)ln⁷2 + (40/3)ζ(3)ln⁴2 + 64ζ(3)ln³2 + 384ζ(3)ln²2 − 1304ζ(5)ln²2 − 200ζ(3)²ln 2 + 2304ζ(3)ln 2 − 4640ζ(5)ln 2 + 576ζ(3)² + 8832ζ(3) − 4800ζ(5) + 15224ζ(7) = 57379.0494155373514…,
  **computed** from the order-8 matching integral (G₈ = 0) on two grids that agree to 430.8 digits and **identified** by PSLQ at 390 digits (the closed form is a conjecture until a by-hand reduction meets it; none is attempted), confirmed to 27 … 30 digits relative by #195's W² values.
  The remainder statement becomes: the eight terms leave Ω₇/(2(D − 3))⁷ with Ω₇ ≈ 2.55 × 10⁷ (an estimate from rem₈·s⁷; not computed). Ω₆ has weights 7, 6, 5, 4, 3 (Ω₅: 6, 5, 4, 3); ζ(7) enters first here.
  **All the ln 2 of Ω₀ … Ω₆ is that of 2^{−4/s}:** W²a² = s 2^{−4/s}F(1/s) with F free of ln 2 through ε⁷ (six terms checked symbolically; the seventh identified by PSLQ in the full basis with ln 2 to ln⁷2, as the same relation as the ln 2-free one).
- **README.md, Limits.** Replace "Ω₆ is an estimate (≈ 3.7 × 10⁶), not computed" by: Ω₇ is an estimate (≈ 2.55 × 10⁷), not computed. The "pure ζ(3) tower" −384, 3072, −36864 is not continued by f₇ (565248, not 589824); if Paper III quotes the pattern it should not.
- **#195** (W² to 60 digits and the fit): its values give Ω₆ to 1.1e-23 with Ω₀ … Ω₅ fixed (X4); it printed no fit of Ω₆.
- **#210 (Paper III).** The eighth term can be quoted as computed and identified, not derived. The ln 2 factor can be quoted as identified through ε⁷ (PSLQ in the full basis weight ≤ 7). `omega_graph.py` and `paper_numbers.py` carry Ω₆ as an estimate; the identified value replaces it (the coordinator's files).

## Recommended tasks

1. **Ω₇ from order 9 and the ln 2 factor at ε⁸.** Add c₉, run the cached chain to k = 9 (one more step of 8 units after the 56: about +15 % per grid; the grids above took 2 h 14 m and 1 h 55 m, so N ≈ 330 for the same 100-digit margin and about 3 h); Ω₀ … Ω₆ planted to 100 digits with Ω₆'s identified form above, two grids agree on Ω₇ to 100 digits (Ω₇ ≈ 2.55 × 10⁷ is an estimate).
   PSLQ for f₈ in the ln 2-free basis of weight ≤ 8 and in bases with ln 2: **the weight-8 monomials π⁸, π²ζ(3)², ζ(3)ζ(5) are not a complete basis of the non-alternating Euler sums (Zagier's dimension 4 at weight 8), so the declaration must add the double zeta value ζ(5,3), or state that a "not found" is for the monomials only**. Why: the structure holds through ε⁷ as an identification. Acceptance: perturbed-basis control for every basis; a relation not found is reported as not found. Depends on #282.
2. **A by-hand reduction of G₈** to masters of weight 7. What: E1 … E3, E4 two levels up, E2W and E5 at the next level and the new steps for L₅ and Y₅, then the finite-part integration by parts (masters J̃₇ and new Euler sums). Why: it would turn f₇ and Ω₆ from identified to derived, and explain why the masters stay ln 2-free in lt = ln(2 cos θ).
   Acceptance: the closed form of f₇ above to 100 digits from the reduction alone; each step checked symbolically or at 30 digits with a control that fires. Depends on #264 (E2W, E5 reviewed).

## Files

- `order8.py` (parts `sym`, `struct`, `cache`, `evals`, `values`, `pslq`, `closed6`, `remainder`; `grid` for development; `--mutate=NAME`, `--cgrids`, `--grids=`, `--bases=`, `--pslq-cap=`, `--save=`, `--save-rel=`, `--rel=`, `--load=`, `--w2raw=`, `--no-real`), `order8_output.txt` (the runs of record concatenated, then the control summaries), this note.
  It imports `PanelB` from `order4.py`, which is unchanged. The commands run as written from `notes/wormhole_largeD` with mpmath and sympy installed, for example `python3 order8.py sym struct cache` (15 s), `python3 order8.py values remainder --load=…/G1.json,…/G2.json,…/B.json --w2raw=…/inputs_w2_growthD.txt` (a second) and
  `python3 order8.py pslq --bases=F1,F2,F3,F4,F5,F6,X1 --load=…/G1.json,…/G2.json --save-rel=rel_a.json` (4 minutes), `python3 order8.py closed6 --rel=rel_a.json,rel_b.json,rel_c.json,rel_d.json --load=…/G1.json,…/G2.json` (a second).
- Raw: `negative-mass-raw/I282-rpi5` (README with the criteria, grids, bases, controls and amendments A1 to A4; `dev_*` the development runs; `run_01` … `run_07` the runs of record; `ctl_*` the controls of record; `G1.json`, `G2.json`, `B.json` the saved Ω's; `rel_*.json` the PSLQ relations;
  `declared_*.json`; `dev_planted_f7.json`, the planted-f₇ end-to-end test; `assemble_output.py`, `queue.sh`, `ctl_summary.py`).

## Controls and checks

`order8.py` has 28 checks and 29 results in the runs of record (PXF ran at H = 10⁸ and at the declared 10¹⁰): E0, CACHE, LN2 (`run_01`); PL0 … PL5, SL, RES, M6, RNG, R7, R8, X4 (`run_05`); PF1 … PF6, PX1 (`run_06a`), PX2, PX3 (`run_06b`), PXF (`run_06c`, `run_06d`); ID6, ZL (`run_07`).
**All pass** (ten runs, all exit status 0). Each check has a control (`--mutate=NAME`, 21 mutations) except PLW, which is the controls' own statement; the targets were computed in development runs before they were declared (raw README, A1, A2, A4).
**Controls of record: 36 runs, 36 on target.** On target means that the set of checks that fire is exactly the declared set (`order8_output.txt`, last section). The groups: CHK 7 (E0, CACHE, LN2: `wrong_h`, `wrong_c6`, `wrong_c7`, `wrong_c8` fire E0, `stale_cache` CACHE, `ln2_struct` LN2),
T1 14 (the 13 grid-dependent mutations and the nominal run on the control grids), T2 5 and T2b 3 and T2bxf 3 (the PSLQ and series mutations on the saved Ω's of record, at caps 150 and 250), and ST 4 (the structure verdict on the real G1, G2: `ln2_leak` fires ID6 and ZL, `f7_rel` ID6, `ln2_zl` ZL).
- **Off target: none.** The three H of the README that differed from the code (X1, X2, X3, XF) were found before the controls of record and are in the defects above; they did not move a target.
- **Declared silences that are results:** `wrong_c8` and `wrong_Om5` are seen only by X4 (Ω₆ moves by less than R7 and R8 can see; as `wrong_c7`, `wrong_Om4` at order 7); `wrong_c7` and `wrong_Om4` do not fire PL4 (Ω₄ does not depend on them); M6 fires only under `coarse_B` (both grids share the other errors); RES fires only under `stale_cache`;
  PLW is silent in every `wrong_Om` run (the wrong input moved Ω_{i+1}); the planted-f₇ tests (`dev_13`, `dev_14`) are the stand-in for the real relation in the ST development runs.

## Coordinator's review

Merged at 57a1cac. Review folder `negative-mass-raw/I282-mac-review` at d951b02b, criteria committed (d725ea4f) before the first run of
`ind282.py`. The changed files are those the brief names.

- **ind282**, 3 of 3. It shares only the printed Ω₃ … Ω₆ with this note's code. From F = (1 + Σ Ω_{k−1}ε^k)·2^{4ε} it rebuilds f₂ and
  f₃ exactly (a planted check of the series) and f₇, which equals the identified closed form to 7.9e-115 relative; the 107 printed
  digits of Ω₃ set that limit. The mutation `zeta3` (the tower's 589824 in place of 565248) and `ln2` fire.
- Not rerun: the two grids (2 h each on rpi5). Their agreement to 430.8 digits and the refit of #195's independent W² to 1.1e-23 are
  the runner's evidence for Ω₆; the review checks the identification against the printed values.
- The README's H for X1 … XF differed from the code's, reported by the runner (A2, A3); XF ran at the declared H and found the
  relation.
- Applied: `README.md` ("Result", "Limits") and Paper III (Ω₆ computed and its closed form identified; F free of ln 2 at ε⁷
  identified; `paper_numbers.py` passes). `omega_graph.py` still marks Ω₆ as an estimate (a follow-up).

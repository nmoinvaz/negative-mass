# #208: Ω₃ from the order-5 matching integral

Status: Ω₃ is **computed** from G₅ = 0 by method B on two grids at two precisions, agreeing to 112.2 digits:
Ω₃ = −652.84010757005612543523930379964736108046174754557825736454474541984957…; it is **identified** by PSLQ at 110
digits in the weight ≤ 4 basis (a relation of height 17280, residual 2.7e-107) as

    Ω₃ = 22π⁴/45 + (16/3)π² ln²2 + (32/3)ln⁴2 − 320 ζ(3) ln 2 − 384 ζ(3),

and the same integral is reduced **by hand** to that closed form, so the identification has a **derivation** (each step
checked symbolically or in mpmath at 30+ digits, checks RD0 to RD8 and a second numerical evaluation, A5). The derivation
uses the exchange identities E1, E2 of #26 and one new step, E3, that removes Y₄. Every number is in `order5_output.txt`
(the four runs of record, raw folder `negative-mass-raw/I208-dh610`, README first) unless it names another raw file. No
sentence here claims novelty; the prior equations used are README.md's order equations (matched.py, order3.py, order4.py),
#26's closed forms and exchange identities, Clausen's functions, the Beta-function moments of ln cos, and Parseval's identity.

## Result

| evaluation | method | working digits | Ω₃ against the closed form |
|---|---|---|---|
| B3 | the hierarchy from its definitions, N = 80, w = 3/10, V = 28 (94 panels) | 180 | 112.2 digits |
| B4 | same, N = 90, w = 1/4, V = 30 (120 panels) | 240 | 133.4 digits |
| B1, B2 | #26's grids (N = 40, w = ½, V = 18; N = 50, w = 2/5, V = 20), reported | 120, 140 | 48.9, 64.4 digits |
| A5 | closed forms of #26 and one integral (E1, E2, E3), 40 digits | 40 + 80 guard | 37.0 digits (against B4 also 37.0) |

- **Ω₃** is the value above (B4; B3 agrees to 112.2 digits, |difference| 3.8e-110). The difference between evaluations is
  a spread over grid and precision, an estimator range, **not a rigorous bound**; nothing here is a rigorous enclosure.
  Against #195's fit it is 3.8e-38 from the 40-digit median −652.8401075700561254352393037996473610805 (the fit's range is
  2.74e-35). In 1/(D − 3): Ω₃/8 = (11/180)π⁴ + (2/3)π² ln²2 + (4/3)ln⁴2 − 40ζ(3) ln 2 − 48ζ(3) = −81.60501344625701568…,
  #195's −81.605013446257016.
- **The planted orders by the same code.** `omega_chain` finds Ω₀, Ω₁, Ω₂, Ω₃ in turn from G₂ … G₅ = 0 with no closed form
  of Y_k or L_k and no input value; Ω₀, Ω₁, Ω₂ against −4 ln 2, 2π²/3 + 8 ln²2 and 80ζ(3) − (8/3)π² ln 2 − (32/3) ln³2:
  B3 113.0 / 112.4 / 112.2 digits, B4 133.8 / 133.3 / 133.2 (the issue's bar: 40). The slope of G_k in Ω_{k−2} is −4π
  to 113.7 (B3) and 135.5 (B4) digits at every k = 2 … 5. With Ω₀ + 1/100 in the orders ≥ 3, Ω₁ moves by 0.027626; with
  Ω₁ + 1/100 in the orders ≥ 4, Ω₂ moves by 0.027726 (controls `wrong_Om0`, `wrong_Om1`; both grids take the same wrong
  input, so the agreement between grids does not see it, the planted checks do).
- **Rounding.** One grid (N = 60, w = 2/5, V = 24) at 200 and at 240 digits agrees to 198.5, 198.6, 200.5, 198.7 digits on
  Ω₀ … Ω₃: rounding costs about 1.5 digits at every order, the fifth included. The grid is the limit (empirically digits ≈ N log₁₀(7.7/w):
  48 at N = 40, 77 at 60, 95 at 70, 112 at 80, 133 at 90).
- **The issue's bar against #195's fit, as worded, is not met** (X1): the string −652.8401075700561254352393037996 has 28 decimals
  and is #195's median truncated, so any Ω₃ that equals the median differs from it by 4.7e-29, against the bar 1e-30. The coordinator
  ruled that X1 stays as worded and is reported not met; X1b, the same bar against the 40-digit median, is met by 3.8e-38.
- **Remainder** (rem = W²a² − [s − 4 ln 2 + Ω₁/s + Ω₂/s² + Ω₃/s³], s = 2(D − 3), against #195's qd W² (60 digits), this Ω₃ from B4):

  | D | s | rem | rem·s⁴ | rem·s⁴ − Ω₄ | (rem·s⁴ − Ω₄)·s |
  |---|---|---|---|---|---|
  | 1000 | 1994 | −8.0797e-10 | −12773.0814 | 21.4185 | 42708.4 |
  | 2000 | 3994 | −5.0238e-11 | −12784.0388 | 10.4611 | 41781.5 |
  | 5103 | 10200 | −1.1816e-12 | −12790.4586 | 4.0413 | 41220.7 |
  | 10000 | 19994 | −8.0049e-14 | −12792.4470 | 2.0528 | 41044.2 |

  rem·s⁴ approaches #195's fitted Ω₄ = −12794.499860023 (the distance falls like 1/s; 0.016 % of |Ω₄| at D = 10⁴); the
  two-point extrapolation from D = 5103 and 10⁴ gives −12794.5179, 0.018 from Ω₄. (rem·s⁴ − Ω₄)·s tends to Ω₅, about
  4.1·10⁴ over the 64 values of D (50 to 10⁴, `order5_output.txt`). **The remainder is Ω₄/s⁴ with Ω₄ = −12794.4999 as #195
  fitted: the five-term series with this Ω₃ leaves the D⁻⁴ remainder.** Ω₄ is #195's fit; it is not computed here.

## The order equations at order 5

README.md's equations (check E0: the ε-expansion of the master equation, sympy, residual 0, now through k = 5):

    dY_k/dl = F_{k−1} = Σ_{j<k} L_{k−1−j} l^j/j!,    (R − 1)L_k′ − L_k/2 = h_k = (1 − 1/R)L_{k−1} + Σ_{j<k} Y_j c_{k−j} + Y_k/2,
    L_k = −2Ω_{k−1} + (R − 1)^{1/2}C_k,    G_k = C_k(∞) = ∫₁^∞ (h_k − Ω_{k−1})(R − 1)^{−3/2}dR = 0,

with c₁ … c₄ as in omega2.md and the new

    c₅ = Ω₀l⁴/24 + Ω₁l³/6 + Ω₂l²/2 + Ω₃l + Ω₄ + l⁵/80 − l⁴/24 + l⁴/(24R),

the coefficient of ε⁵ of (1 + εΩ)e^{εl} − (½ + ε)(1 − 1/R)e^{−εl} − e^{−εl}/(2R). G_k is affine in Ω_{k−2} with slope −4π; Ω₄ drops
out of G₅ (it is h₅(1), subtracted), so **G₅ = 0 gives Ω₃**. Method B, as in omega2.md, at order k: Y_k from the left, the tail C_right from
the right, L_k = −2h_k(1) − u C_right, G_k = C_right(0), on Chebyshev panels in v = √l (`PanelB` of order4.py, shared). What the evaluations
share: the table (E0), the analytic-at-throat condition, `PanelB` and mpmath; B3 and B4 differ in N, w, V and precision. #195's fit (the full
equation in a Rust crate, no matching) shares none of it and is the independent reference; the closed forms of Ω₀ … Ω₂ (#26) are the planted truth.

## PSLQ (identification)

mpmath's `pslq`, Ω₃ of B4 at d = 110 digits (the grids share 112.2), tolerance 1e-98, `maxsteps` 10⁶; budget n·log₁₀(2H + 1) + 1.

| basis | constants besides Ω₃ | n | H | budget | result |
|---|---|---|---|---|---|
| P1 pure weight 4 | π⁴, π² ln²2, ln⁴2, ζ(3) ln 2, Li₄(½) | 6 | 10⁹ | 56.8 | not found |
| P2 weight ≤ 4 | P1 and ζ(3), π² ln 2, ln³2, π², ln²2, ln 2, 1 | 13 | 10⁶ | 82.9 | **found** (below) |
| P3 weight 4, extended | P1 and πζ(3), π³ ln 2, π ln³2, π²G, β(4), G ln²2, πG ln 2 | 13 | 10⁶ | 82.9 | not found |

P2: (45, −22, −240, −480, 14400, 0, 17280, 0, 0, 0, 0, 0, 0), i.e. 45Ω₃ − 22π⁴ − 240π² ln²2 − 480 ln⁴2 + 14400 ζ(3) ln 2 + 17280 ζ(3) = 0, residual 2.73e-107; its
own budget 13·log₁₀(34561) + 1 = 60 and the declared 82.9 are below 110; the basis with π⁴ + √2·10⁻³ finds nothing. Each basis recovers its planted relation
(a random integer combination of its constants, 3 to 4 digit coefficients, `random.Random(seed)`) exactly and finds nothing in the perturbed basis: the tool
sees relations of the size found. The relation holds to the digits of the grids, not only of PSLQ: B3 112.2, B4 133.4 digits (X3). Li₄(½) does not appear.
The weight-3 term ζ(3) is why the pure weight-4 basis P1 found nothing at #195's 35 digits (`dev_25_pslq_omega3`) or here, and why #195's 13-element try
returned a chance relation: it needed about 60 digits for coefficients up to 17280, and had 35. Clausen values at π/2 (G, β(4), π²G) do not occur
(P3's first, declared guess): the arguments of every special function here are π − 2θ and e^{2iθ} at θ = 0 and π/2 only.

## The reduction by hand

Notation as omega2.md: S = sin²θ, ℓ = ln cos θ, u = tan θ, C = Cl₂(π − 2θ), D₃ = Cl₃(π − 2θ) + ¾ζ(3), Lc2 = ∫₀^θ ln² cos, q = ln 2, J_n = ∫₀^{π/2}ℓⁿdθ,
Ω₀ = −4q. Each line is a check of `order5.py reduce5 a5` (RD0 … RD8, A5; raw README, amendment A3).

1. **E1, E2** (#26) give G₅ = ∫₀^{π/2}I₅dθ with I₅ = 8ℓh₄/S + 2Σ_{j=1}^{4}L_{4−j}l^j/j! + 2(Σ_{j<5}Y_jc_{5−j} − Ω₄)/S; h₄ holds Y₄, which has no closed form.
2. **E3: Y₄ by parts.** Y₄ enters only through (4ℓ + 2c₁)Y₄/S, c₁ = Ω₀ − 3ℓ − S, which is 2P′Y₄ with P = cot θ(ℓ − Ω₀), P′ = (Ω₀ − ℓ − S)/S (RD0, residual 0).
   Then ∫2P′Y₄dθ = −2∫PY₄′dθ = −4∫(ℓ − Ω₀)F₃dθ, since Y₄′ = 2 tan θ F₃ and P tan θ = ℓ − Ω₀. The boundary terms vanish: Y₄ = O(θ²), P ~ 1/θ at 0; Y₄ = O(ℓ⁴), P ~ x ln x at π/2.
   Checked at j = 3 on the closed Y₃ and F₂ (RD1, 43.5 digits); at j = 4 the identity is the same by parts and its use is checked by A5 against B. I₅ now needs only Y₀ … Y₃
   and L₀ … L₃, all closed in #26 (L₃ through Lc2).
3. **The decomposition.** At Ω₀ = −4q: I₅ = rest + cC·C + cD·D₃ + cL·Lc2 with cL = −1024 q u, cD = 32(6ℓ − 8q) + 32(11ℓ² + 8qℓ − Ω₁)/S, cC = 64u(3Sℓ² + 24Sℓq − 6Sℓ + 16Sq² + 8Sq − 11ℓ² − 8ℓq + Ω₁)/S,
   and rest = (A(ℓ) + S·B(ℓ))/S + 128 q Ω₁θu, A(0) = 0, affine in Ω₃ with coefficient 8ℓ/S (slope −4π) (RD2: against method A's numerical integrand, 45 digits).
4. **The D₃ term, by parts.** With K′ = −2ℓ (K = 2qθ − C) and F_n = ∫ℓⁿ/S (F₀ = −cot θ, F₁ = −ℓ cot θ − θ, F₂ = −ℓ² cot θ + K), cD has the antiderivative
   A_D = 32(−3K − 8qθ) + 32(11F₂ + 8qF₁ − Ω₁F₀) (A_D′ = cD, sympy residual 0), A_D(π/2) = 0 and D₃(0) = 0, so ∫cD D₃ = −2∫A_D C (RD3, 45 digits).
5. **The C terms combine.** With K = 2qθ − C: cC − 2A_D = u(64(16q² + 8q + Ω₁) + 128(8q − 3)ℓ − 512ℓ²) + 512 C (RD4, sympy residual 0): the cot θ and θ terms cancel exactly,
   as the C-terms did in #26.
6. **The Clausen integrals.** ∫uℓⁿC dθ = (2J_{n+2} + 2qJ_{n+1})/(n + 1) (by parts, C′ = 2ℓ + 2q, ∫uℓⁿ = −ℓ^{n+1}/(n + 1), boundary terms 0), and ∫₀^{π/2}C²dθ = π⁵/360 (Parseval: ½∫₀^π Cl₂² = (π/4)ζ(4)); RD5, 44.9 and 45 digits.
7. **The Lc2 and θ terms.** −1024q·u·Lc2 and 128qΩ₁·θu each diverge at π/2 but together u·[−1024q(Lc2 − J₂) + 128qΩ₁(θ − π/2)] does not: the constants cancel because −1024qJ₂ + 64πqΩ₁ = 0
   with Ω₁ = 2π²/3 + 8q² (sympy 0). By parts ∫u[…]dθ = −1024qJ₃ + 128qΩ₁J₁ (RD6, 45 digits).
8. **Moments and assembly.** ∫ℓᵇ/S = −bJ_{b−1} (b ≤ 5) and J_n (n ≤ 4) from the Beta function (RD7, 42 to 43 digits), so G₅ = 4π(Ω₃ᶜ − Ω₃) with Ω₃ᶜ the closed form above (RD8, sympy, difference 0, slope −4π).
   Derived: G₅(Ω₃ = 0) = −8203.830743643434750828152377131025800843.
9. **A5**: method A at 40 digits, G₅ by quadrature of the Y₄-free integrand with the closed forms of order4.py, gives Ω₃ = −652.84010757005612543523930379964736101 (37.0 digits against B4, the closed
   form and the derived G₅(0); slope −4π to 38.4 digits). At 20 and 60 digits (`dev_28` … `dev_33`): 3e-15 and 7.4e-54 from B4.

**What each status word covers.** *Derived* (by hand, each step checked): steps 1 to 8, given E1, E2 (#26: N4a, N4bg), E3 (RD1 and the boundary argument), the closed forms of Y₃ and L₃ (#26),
Ω₀, Ω₁, Ω₂ (#26's derivations), Parseval for ∫C² and the Beta-function moments. *Numerical evidence*: the evaluations B3, B4, A5 and their agreements. *Identification*: the PSLQ relation, which the
derivation supersedes. Not derived here: E2's boundary terms (checked numerically in #26), the boundary argument of E3 at j = 4 (analytic, not checked at j = 4 numerically), Y₄ and L₄ in closed form
(not needed). A5 is not independent of the derivation (same E1 to E3 and closed forms); **B3 and B4, which use none of them, are the independent check**, and #195's fit (shares no machinery) agrees to 3.8e-38.

## Controls and checks

`order5.py` has 24 checks in the runs of record (E0; PL0 … PL2, SL, M3, X1, X1b, Q1 … Q3, R3, R4; X3; RD0 … RD8, A5; PL3 is the control's own statement and reports only under `wrong_Om0`, `wrong_Om1`). **23 are
met; X1 as worded is not.** Each has a control (`--mutate=NAME`, 15 of them) whose targets were computed in development runs before being declared (raw README). Silent targets are listed there:
`wrong_c5` and `wrong_Om2` do not move the planted orders (both enter at k = 5 only) or M3 (both grids share the input); `short_tail` and `wrong_Om*` leave M3 silent for the same reason; `coarse_B` fires Q1 … Q3 only
because PSLQ is given the 15 digits the grids share; `wrong_AD` is seen by RD3 and RD4 but not RD8, which uses RD4's coefficients as printed.

## Defects and discrepancies

- **X1 as worded is not met** (above): a defect of the criterion's wording, not of Ω₃; X1b was added by the coordinator's ruling (raw README A1).
- A first version of the PSLQ perturbed-basis control shifted a constant by 10⁻³; P2 contains 1, so the true relation scaled by 1000 was found (coefficients ≤ 10⁶). Found by running the real-target branch on
  synthetic targets (`dev_23`), changed to √2·10⁻³ before the runs of record (`dev_25`).
- The first remainder code read D = 1000, 2000 from the 64-value list, which does not have them (`KeyError`); fixed by the 60-digit qd values of #195's `run_02`, as #26 did (`inputs_w2_growthD.txt`).
- #195's `w2_output.txt` gives Ω₃ to 14 digits only; the 40-digit median is in the raw log `I195-dh610/dev_25_pslq_omega3.log` (line 1), which is a development log of #195.
- The bars of R3 and R4 (falling and 0.05 %; 0.1) were set knowing the order of magnitude of #26's table (declared in the raw README).
- The identification was found after the criteria were written; X3 and the reduction's checks RD0 … A5 were declared (amendments A2, A3) before their runs, with bars below the development values.

## What README.md's "Result" and "Limits", #201, #25 and Paper III should now say

- **README.md, Result.** Add to the series the fourth term of the 1/(D − 3) expansion: Ω₃/8 = (11/180)π⁴ + (2/3)π² ln²2 + (4/3)ln⁴2 − 40ζ(3) ln 2 − 48ζ(3) = −81.6050134462570…, **derived** from the order-5
  matching integral (this note: G₅ = 0, one integral of closed forms after E1, E2, E3, reduced to Clausen integrals and log-cosine moments), confirmed to 112 digits by an evaluation that uses none of the
  reduction and to 38 digits by #195's fit. The remainder statement becomes: the five terms leave Ω₄/(2(D − 3))⁴ with Ω₄ = −12794.4999 (#195's fit), table above.
- **README.md, Limits.** "The third term's closed form is an identification, not a derivation" stays for Ω₁ until #25 is accepted. Add: Ω₄ is a fit (#195), not computed or identified.
- **#201** (identify Ω₃ in weight-4 constants, fit route at 60 digits): the target is met by this issue to 112 digits (B3 against B4; B4 alone is 133 digits from the closed form); Ω₃ is in weight ≤ 4, not pure weight 4 (the ζ(3) term). The fit route is no longer needed for Ω₃;
  its value now is Ω₄, for which Ω₃ is exact (a fit of the 64 values of W² with Ω₀ … Ω₃ fixed).
- **#25.** Unchanged by this note; E3 does not enter order 3.
- **Paper III.** The fifth term can be quoted as derived with the same matching integral as the third and fourth; its derivation needs E1, E2, E3, one Clausen integral (∫C² by Parseval) and the log-cosine
  moments. The paper's "how this work was done" section should record that the closed forms and their checks were done by a language model with sympy and mpmath, as AGENTS.md asks.

## Files

- `order5.py` (parts `sym`, `evals`, `pslq`, `remainder`, `closedform`, `reduce5`, `a5`; `grid` for development; `--mutate=NAME`, `--save`, `--load`, `--w2raw`), `order5_output.txt` (the four runs of record
  concatenated), this note. It imports `PanelB`, the closed forms and the moments from `order4.py`, which is unchanged.
- Raw: `negative-mass-raw/I208-dh610` (README with the criteria, controls and amendments A1 to A3; `dev_*` the development runs including the sympy scripts `dev_34`, `dev_35` and the numerical method A `dev_28`;
  `run_01` … `run_04` the runs of record; `inputs_w2_growthD.txt`; `assemble_output.py`).

## Coordinator's review

Under review: `wip/i208-omega3` at fd15081, raw `I208-dh610` at ace0b32. Review folder
`negative-mass-raw/I208-mac-review` (at e9cf98a), whose README holds the criteria, committed before the first run of
`ind208.py`. The changed files are those the brief names.

- **Rerun on the Mac** of the four runs of record, from their `.env` files. Every line equals dh610's except the
  timings, with the same exit codes. X1 fails as worded there too.
- **ind208**, run 1, 3 of 3. It is the review's ind26 carried to order 5, and shares with this note only mpmath, sympy
  and README.md's matching conditions.
  - K1: the order equations for k = 0 … 5, derived again from #195's ℓ-system, equal `order5.py`'s h_k and dY_k/dl (sympy
    residual 0). This includes the new c₅.
  - K2: a Taylor integration of that hierarchy in √l, with C_k from the throat outward, gives
    Ω₃ = −652.840107570056125435239303799647361080405. That is within 5.6e-38 of the closed form, and within 9.5e-38 of
    #195's median. Ω₀ … Ω₂ match their closed forms as in #26's review.
  - K3: the slopes are −4π at k = 2 … 5, and G₅(Ω₃ = 0)/(4π) matches the closed form to 8.6e-41 relative.
  - The mutations `k1-B`, `short-V` and `drop-z3` fail exactly their declared targets. Without its −384ζ(3) term the
    closed form misses the computed Ω₃ by 462.
- **Not met, as ruled:** X1 as worded. Its reference string was the coordinator's truncation of #195's median to 28
  decimals; X1b, against the 40-digit median, is met.
- *Added by the coordinator:* Ω₃ is not of uniform weight. Its −384ζ(3) term has weight 3, so README.md's observation
  about Ω₀ … Ω₂ stops at three terms; README.md says so.
- **What the review does not cover.** As in #26's review, ind208 checks the end point of the reduction, not the
  exchange identities. E3 is new here, and #209 now covers E1, E2 and E3.
- **Follow-ups.** #201's question, Ω₃ in weight-4 constants, is answered here, and #201 is closed. #210 (Paper III) now
  includes the fifth term. Ω₄ remains a fit; no task computes it from G₆ yet.

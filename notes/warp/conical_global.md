# Do conical shifts with det R₃ = 0 exist on all of x/ρ ∈ ℝ?

November 2026, issue [#167](https://github.com/nmoinvaz/negative-mass/issues/167). Script
[`conical_global.py`](conical_global.py), its stdout [`conical_global_output.txt`](conical_global_output.txt), the Rust crate
[`conical_global/`](conical_global/). It continues [`conical_limit.md`](conical_limit.md) (#159) and
[`det_classify.md`](det_classify.md) §6. Status words as in AGENTS.md: **proved**, **checked symbolically**, **numerical
evidence**, **conditional**, **conjecture**. Nothing here is certified: no enclosure is computed. Every number quoted below is in
the output file.

## Answer

**No.** I found no real solution of N = 0 on all of η ∈ ℝ, regular at η = 0, in the two-parameter family F = F_∞ − 2F_∞/η² +
(C + F_∞(F_∞² + 7) ln η)/η⁴ + … at both ends. For 187 pairs (F_∞, C) the inward integration from the series at η = 200 or beyond
never reaches η = 0 regularly: it ends at a singular η* > 0 in 172 cases, and in 15 cases it reaches η = 10⁻⁶ with F′ ∝ 1/η, which is
not regular either. The solutions that do cross η = 0 regularly (found from η = 0 outward) do not tend to a limit F_∞: for F(0) > 1
and F′(0) ≥ 1 they grow like ln ln η. The constants F ≡ F_∞ are global, regular and even, and are outside the family (a₂ = 0); they are
not counted, as in #159.

| question | answer | status |
|---|---|---|
| a global F in the (F_∞, C) family, same limit at both ends | none on the grid F_∞ ∈ {±0.25 … 4}, C ∈ [−3000, 3000], 187 points; see the next three rows | **numerical evidence** |
| first singular η of each grid point | F → 0, F′ → ∞ at 166 points, η* from 4.6e-5 to 109.4, F ~ (η − η*)^α with α = 2/(η*² + 4) (measured to 1.2e-3 relative); F″ ~ ln\|η − η*\| at the curve F² = η² + 1 at 3 points; an F″ pole at η* ≈ 1e-6 to 4e-6 at 3 points | **numerical evidence**; the exponent's formula by hand from the printed N |
| inward flows that reach η = 10⁻⁶ | 15 points (F_∞ = 0.25, 0.5, 0.75 and their mirror, large \|C\|): F′ ≈ p₀/η with p₀ = ηF′ of 9e-5 to 6e-2, falling as \|C\| grows; p₀ → 0 only as C → ∞, where F tends to the constant F_∞ | **numerical evidence** |
| solutions through η = 0 | with F′(0) ≠ 0 and F(0) ∉ {0, ±1}: an analytic two-parameter family, F″(0) = −F′(0)²/F(0). F′(0) = 0 and F(0) ∉ {0, ±1}: only F ≡ const. F(0) = 0: none but F ≡ 0. F(0) = ±1: a one-parameter family with F′(0) ≠ 0, and one even solution, F″(0) = ±4/5 | **checked symbolically**; the exclusions rest on **hypothesis H₀** (a leading term c η^p), so they are **conditional** |
| the even solution through F(0) = 1 | F = 1 + 2/5 η² − 22/125 η⁴ + 998/10625 η⁶ − …; walked from η = 1/20 to 200 with no singularity, and there F ~ ln ln η grows without a limit | exact series **checked symbolically**; the growth **numerical evidence** |
| forward map from η = 0 | 100 solutions (F(0), F′(0)) walked to ±200: singular, or log-log growth; none FAMILY on either side | **numerical evidence** |
| the regime F ~ ln ln η | an exact solution of the ν⁰ part A₀ = p²(p² − q)(p² + q) = 0, q = −p² (p = 1/(ln η + t₀)); unbounded, so outside #159's family and outside its assumption H (F − F_∞ → 0) | exact solution by hand; the approach of regular solutions to it **numerical evidence** |

Not covered (§5): the parts of (F_∞, C) space off the grid, solutions that are only C² or C¹ near the singular curves, and any
approach to η = 0 that is not a leading power. Closure on a global solution is moot: there is none to close.

## 1. The equation and how it is integrated

N is `det_classify.conical_N()`, cubic in F″: N = a₃F″³ + a₂F″² + a₁F″ + a₀ with

    a₃ = 2Fη²(F² − η² − 1)
    a₂ = F′η(4F³ + F²F′η³ + 6F²F′η − 6Fη² − 4F − F′η⁵ − 3F′η³ − 2F′η)
    a₁ = 2F′²(F³ + F²F′η³ + 6F²F′η + 2FF′²η⁴ + 4FF′²η² − 3Fη² − F − F′η⁵ − 3F′η³ − 2F′η)
    a₀ = F′³(F²F′η² + 2F²F′ + 4FF′²η³ + 8FF′²η − 2Fη + F′³η⁶ + 4F′³η⁴ + 4F′³η² − F′η⁴ − 3F′η² − 2F′)

The crate takes the state (F, F′, F″) and, at each step, the Taylor series about the step's start from the implicit equation: order
k ≥ 1 of N = 0 is solved for the coefficient F″ₖ, with one division by ∂N/∂F″ at the step's start (`poly.rs`). The order-0 term of
N is not imposed, so N is constant along the flow, and N = 0 comes from the end data. No hand-differentiated F‴ is typed. The step
is multishoot's rule (the last two terms below 1e-32, order 30, double-double). The end states at η = ±L are the series of
`conical_limit.formal_series` recomputed here through η⁻¹⁷ (`series_stages`: every stage solved, one free parameter a₄,₀, a₄,₁ =
F_∞(F_∞² + 7), 36 terms, and `src/tables.rs` equals the regenerated table). #159's series stopped at n ≤ 9; n ≤ 13 failed A2 as
first worded (3.2e-19 at (1, 0)), n ≤ 17 gives 4e-25.

**The matching.** With N = 0 at both ends, F″ at the matching node is redundant with (F, F′), so the left end carries a nuisance
unknown e (F″ → F″ + e at η = −L), which the matching returns as 0. The unknowns are C⁻, e, C⁺ and the global F_∞, pinned. Planted
truth through the same code: the exact solution α + (a + bη)/(1 + η²) of (1 + η²)F″ + 4ηF′ + 2(F − α) = 0, matched between its
two truncated expansions at ±40, is recovered to 6e-28 (C⁻), 3e-29 (C⁺) and 6e-34 (e), and the same with the cubic M·N_p
(M = 1 + F′² + F″²) to 8e-28. No grid point needed the matching, because none was regular (`connect`: not applicable).

## 2. Local analysis at η = 0 and at the curve F² = η² + 1

N(η = 0) = 2F′²(F² − 1)(FF″ + F′²) and ∂N/∂F″(η = 0) = 2F′²F(F² − 1) (`eta0_local`, **checked symbolically**). So η = 0 is a
singular point of the implicit equation exactly where F′(0) = 0 (and where F(0) ∈ {0, ±1}); the claim that F′(0) = 0 makes a solution
even does not follow from uniqueness, and the coordinator's hint that it does was withdrawn. What passes through η = 0:

| local form at η = 0 | result | status |
|---|---|---|
| F′(0) ≠ 0, F(0) ∉ {0, ±1} | F″(0) = −F′(0)²/F(0); each order solves linearly for the next Taylor coefficient (two rational points, 7 orders): an analytic two-parameter family (F(0), F′(0)) | **checked symbolically** |
| F − f ~ cη^p, p > 1, f ∉ {0, ±1} | leading coefficient 2f(f² − 1)c³p⁵(p − 1) ≠ 0: only c = 0, i.e. F ≡ f | **conditional** on H₀ |
| F(0) = 0 | N(0) = −2F′⁴ forces F′(0) = 0; then F ~ cη^p has leading coefficient −2c⁴p⁵(2p − 1), zero only at p = 1/2: F ~ c√η, F′ unbounded, not regular | **conditional** on H₀ |
| F(0) = ±1 | order 1 gives F″(0) = 1/6 − F′(0)² (for F = 1) for each F′(0) ≠ 0, one parameter; the even branch F = 1 + b₁η² + …, order 4: 64b₁³(5b₁ − 2) = 0, so b₁ = 2/5 (F″(0) = 4/5) or F ≡ 1; the whole series follows from b₁ = 2/5, with no free coefficient (exact rationals to η⁶⁰, N = 0 exactly through order 62) | order 1 and order 4 **checked symbolically**; the exclusion of other exponents **conditional** on H₀ |
| F′ = p₀/η (F ~ p₀ ln η) | not a solution by itself: N(F = f + p₀ ln η, p₀/η, −p₀/η²) = p₀⁴(p₀²(η² + 2)² − 4F²)/η⁴; the flows that reach η = 10⁻⁶ show it as a slowly varying p = ηF′ | symbolic identity; the flows **numerical evidence** |

**Hypothesis H₀.** F − f has a leading term cη^p (real c ≠ 0, p > 0) whose derivatives are the termwise ones. A logarithmic
factor or an oscillating term is not covered.

So a global solution must cross η = 0 with F′(0) ≠ 0, or be one of the constants, or be the special solution through F(0) = ±1. The
special solution has no limit (§4). A crossing with F′(0) ≠ 0 is not even unless F′ changes sign, so the two ends have
different F_∞ unless a coincidence holds; the forward scan (§4) found no solution with a limit on either side.

**The curve κ = F² − η² − 1 = 0, where a₃ vanishes.** At a point of it the cubic loses a root to infinity and the other two roots are
analytic branches (∂N/∂F″ ≠ 0 there in general), so a solution on a finite branch crosses the curve regularly. The branch the
inward flows follow is the escaping root: at 3 grid points ((1.5, 3000), (2, 3000) and (−2, −3000)) it ends with F″ ~ ln|η − η*|
(fitted exponent of F″: γ = 0.99, 0.96, 0.96), F and F′ finite, and \|F² − η² − 1\| < 1e-8 there. At seven other points the first pass
(stop at R ≤ 1e-6) also reported such an end, and the refined pass (R ≤ 1e-8) walked through it: the nearest singularity of the
Taylor series was a complex pair closer than 1e-6 to the real axis, and the solution crosses κ = 0 analytically. That is why the pass
of record is the refined one (README A1). A solution that is only C¹ or C² at the log points could continue, and is not covered.

## 3. The scan of (F_∞, C)

The grid: F_∞ ∈ {0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, −0.5, −1, −2} and C ∈ ±{0, 1, 3, 10, 30, 100, 300, 1000, 3000}, 187 points. Each
point starts at the smallest L ∈ {200, 800, 3200, 12800, 51200} where the series residual |N|/Σ|terms| is ≤ 1e-20 (115 points at 200,
58 at 800, 14 at 3200), is walked inward, and ends as listed in the output. The singular η is located by the Taylor
coefficients (Domb–Sykes extrapolation of the ratios cₖ/cₖ₋₁ in 1/k; the error bar is the spread over the last five steps, **not an
enclosure**), to a bar of at most 1.9e-12 at every singular point of record. The kinds:

| kind | count | what it is |
|---|---|---|
| `F->0` | 166 | F → 0 and F′ → ∞ at η*, F ~ c(η − η*)^α. The leading balance of N at small F is −η²(η² + 1)F″²[(η² + 2)F′² + 2FF″] = 0, which gives α(η*² + 4) = 2. Measured against 2/(η*² + 4): (1, 0): η* = 1.28193326312, 0.35417 against 0.35440; (4, 0): 3.03548392532, 0.15134 against 0.15135; (1, 100): 0.0118195947859, 0.50056 against 0.49998; worst over the 166 points 1.2e-3. η* is smooth in C and falls to 5e-5 at C ≈ 1000 for F_∞ = 1 |
| `a3->0` | 3 | the escaping root at F² = η² + 1, F″ ~ ln\|Δ\| (§2) |
| `other` | 3 | η* = 1.2e-6 (0.75, 300), 3.6e-6 (±1, ±3000): F″ has a simple pole, γ ≈ 0, with F = 0.13 to 0.25, F′ ≈ 1e5: the end of the F′ ≈ p₀/η flows just above the floor |
| FLOOR, `log` | 15 | the walk reached η = 10⁻⁶ with F′ ≈ p₀/η; p₀ = 9e-5 (0.25, 3000) … 6e-2 (0.5, 100) at 10⁻⁶ |

No point is THROUGH, so there is no GLOBAL point and nothing for the shooting solve to refine. The bisections between adjacent
grid values of different class end at η* ≈ 2e-5 to 5e-5 on both sides (a label boundary of the kinds, with η* continuous in C), or
between F->0 at η* = 0.56 and a3->0 at 1.08 at F_∞ = 1.5, C ≈ 1373.5, and F_∞ = 2, C ≈ 1724.1; they are not separatrices of solutions that
are regular at 0. The symmetry (F_∞, C) → (−F_∞, −C) holds for all 51 compared pairs.

The honest limit of this table: where the singularity sits at η* < 10⁻⁶ the walk stops at the floor and calls it FLOOR, and the F_∞ of
the 15 FLOOR points are below 1, where F → 0 is the only singularity the walk could meet (a₃ has no real zero off η = 0 for |F| < 1).

## 4. The forward map from η = 0

For F(0) ∈ {0.1, 0.25, 0.5, 0.75, 1.5, 2, 3, 5, −0.5, −2} and F′(0) ∈ ±{0.01, 0.1, 1, 10, 100} with F″(0) = −F′²/F (max |N(0)| 4e-24), the
walk goes to +200 and −200. By the parity of N the minus side of (f, f₁) equals the plus side of (f, −f₁) (checked: 0 differences,
η* to 1e-8). The 200 walks end as: `F->0` 50, `fold` 56 (F″ has a square-root branch point, ∂N/∂F″ → 0 on the branch), `a3->0` 40, `other` 6,
the κ locus with overflow 18 (η* known to 1e-6 only), and 30 in the log-log regime. None is FAMILY.

The pattern: 0 < F(0) < 1: F′(0) > 0 ends at a fold at η* from 8.5e-4 (F′(0) = 100) to 1.6 (F′(0) = 0.01), F′(0) < 0
ends at F → 0 (at η* from 4e-3 to 19; the mirror cases for F(0) = −0.5); F(0) > 1: F′(0) < 0 ends at a3 → 0, F′(0) ∈ {0.01, 0.1} at the κ locus, and F′(0) ≥ 1 enters the log-log regime. So the solutions
that cross η = 0 regularly and reach η = 200 are in the regime p = ηF′ ≈ 1/ln η, q = η²F″ + ηF′ ≈ −p² (q/p² within 0.05 of −1; sign-flipped for F < 0), where
F ≈ ln ln η + const. This is the ν⁰ part A₀ = p²(p² − q)(p² + q) = 0 with q = dp/dt, t = ln η: p′ = −p², p = 1/(t + t₀). It is
unbounded, so it is not a limit solution, and #159's assumption H (F − F_∞ → 0) did not include it. The exact even solution through
F(0) = 1 is in this regime (F = 3.0222 and p = 0.1657 at η = 200, started from η₀ = 1/20 and 1/10 with the same F and p to 1e-12).

## 5. Not decided and not covered

- Off the grid: |F_∞| < 0.25 or > 4 (the series is singular as F_∞ → 0: a₆,₀ has 1/F_∞), |C| > 3000, and the continuum between the 187 points. A
  regular crossing would sit on a curve in (F_∞, C), and the scan cannot exclude a narrow one; the bisections and the smoothness of η*(C) suggest
  there is none in the range.
- Solutions that are only C¹ or C² at the log points of the curve κ = 0 (F″ ~ ln|Δ|), and those that cross it on the finite branches.
- Solutions whose approach to η = 0 is not a leading power (H₀), and oscillating approaches at ±∞ (H of #159).
- The forward grid is 100 points; the log-log regime is described, not proved, and whether a solution of it can be continued to a
  limit at huge η (the family needs |p| ≈ 4F_∞/η²) was not followed beyond η = 10⁴.
- Everything is numerical evidence; the exclusions at η = 0 are conditional on H₀.

## 6. Planted truth, checks and controls

17 checks pass in the run of record; the lines `NOT MET (as first worded)` in the output are results and do not gate the exit code.

| check | what it tests | result |
|---|---|---|
| `n_table` | the crate's N at 40 rational points against `det_classify.conical_N()` and `conical_limit.build_N()` | 2.1e-32 (threshold 1e-26) |
| `series_stages`, `series_order` | the formal series through η⁻¹⁷; N of it vanishes through ν²³ and not at ν²⁴ at 12 exact rational pairs | pass |
| `series_residual` | the residual at ±L and the crate's state against mpmath at 50 digits | residual ≤ 1e-20 at the adaptive L; difference 2.3e-32. **As first worded (literal start 200): 2 of 6 pairs fail**, (4, 3000) at 2.9e-14 and (0.25, −3000) at 3.8e-5 |
| `plant_linear`, `plant_cubic`, `plant_blowup` | the matching and the locator through the same code | 6e-28, 8e-28; the locator on F″ = 2F³ and FF″ + F′² = 0 to 8e-17 and 7e-12 |
| `step_halving` | order 30, 1e-32 against order 36, 1e-38 at η = 3 (10 for (0.5, −30)) | 4.9e-27 at most |
| `scan`, `n_along_path`, `symmetry`, `bisect`, `connect` | the grid | §3. **As first worded: the start residual exceeds 1e-18 at 58 of 187 points, and two points have a bar 2.4e-8 at the declared stop** |
| `eta0_local`, `forward`, `even_pm1` | §2, §4 and the even solution | pass |
| `independent` | own mpmath integration (25 digits, `odefun`) of (F, F′) against the crate from η = 200 to 4 or 10 | 9.4e-18 at most; shares N and the series with the crate |

| control `--mutate=` | change | target check |
|---|---|---|
| `eq-coeff` | one coefficient of N | `n_table`, `series_residual` |
| `series-coeff` | a₄,₁ × 1.01 | `series_residual` |
| `a16-coeff` | a₁₆,₀ × 1.001 | `series_order` |
| `match-sign` | the mirror end takes +F′ | `plant_linear`, `plant_cubic` |
| `plant-locator` | the side test inverted | `plant_blowup` |
| `halving` | first plan order 14, tol 1e-12 | `step_halving` |
| `n-drift` | order 3 of N not imposed | `n_along_path` |
| `indep-eq` | −6 → −5.9999 on FF′²F″η² | `independent` |
| `sym-break` | (−F_∞, +C) | `symmetry` |
| `bisect-end` | the wrong class at the upper end | `bisect` |
| `eta0-sign` | N(0) expected with 3 in place of 2 | `eta0_local` |
| `fwd-class` | log-log tested at the wrong sign | `forward` |
| `even-series` | b₃ × 1.001 | `even_pm1` |

The controls are in the output with their exit codes. A control runs only the checks it names. From the repository root:

    python3 notes/warp/conical_global.py --fresh        # the checks and the 13 controls: 1334 s on the Mac, peak 205 MB (series 185 s of it)
    python3 notes/warp/conical_global.py --no-controls  # the 17 checks, 774 s with the series cached
    python3 notes/warp/conical_global.py --mutate=NAME --only=CHECK

The crate builds with `cargo build --release` (multishoot 4.2.0, multisolver 0.2.2, multidouble 5.1.0, Rust 1.94.0). The run of record
exited 0: 17 checks passed, 13 of 13 controls fired on their target only, and the three `NOT MET (as first worded)` lines are printed.

## 7. What the code proves, and what is assumed

Checked symbolically (sympy): the identities of §2 (`eta0_local`), the series stages and the order of N's first nonzero term
(`series_stages`, `series_order`), the exact even series (`even_pm1`, exact rationals). Numerical evidence in double-double: every
walk, every η*. By hand, not coded: the leading balance at small F that gives α = 2/(η*² + 4); the exact solution p = 1/(t + t₀) of A₀ = 0.
Assumed: N from #65's D as in #159 (not rederived here); H₀ for the local exclusions; sympy 1.14.0, mpmath 1.3.0.

## 8. What `conical_limit.md` §6 and `det_classify.md` §6 should now say

For `conical_limit.md` §6, items (a) and (b): "(a) A global solution on η ∈ ℝ with the same F_∞ at both ends in the family of §3: none
found on 187 (F_∞, C) with \|F_∞\| ∈ [0.25, 4], \|C\| ≤ 3000: each ends at a singular η* > 0 (F → 0 with F′ → ∞ in 166 cases) or at the floor
with F′ ∝ 1/η ([`conical_global.md`](conical_global.md), #167; **numerical evidence**). (b) At η = 0 the equation is singular where F′(0) = 0:
a solution through η = 0 with F(0) ∉ {0, ±1} has F′(0) ≠ 0 and F″(0) = −F′(0)²/F(0), or is constant (conditional on a leading-power form); the only nonconstant even
regular solution has F(0) = ±1, F″(0) = ±4/5, and grows like ln ln η. Solutions that cross η = 0 regularly and reach large η follow
q = −p² (F ~ ln ln η), which is outside the family and outside assumption H."

For `det_classify.md` §6, after the "Update from #159" paragraph: "**Update from #167:** the limits of #159 exist, but no solution in that family was
found on all of x/ρ ∈ ℝ on the scanned grid; the solutions regular at x = 0 found by integrating outward have no limit (F ~ ln ln(x/ρ)). Localization in x for conical
shifts is therefore not realized by the family of #159 on the scanned set, and is not excluded for F_∞ and C off it."

## 9. Lemmas for Lean (#135)

Polynomial identities proved by `ring` over ℚ, stated for a later formalization (no Lean written): (L1) N(F, F′, F″; 0) = 2F′²(F² − 1)(FF″ + F′²);
(L2) ∂N/∂F″(F, F′, F″; 0) = 2F′²F(F² − 1); (L3) with X = η^p, N(f + cX, pcX/η, p(p − 1)cX/η²; η) has the X³ coefficient
−2c³fp⁵((f² − 1)(1 − p) + pη²)/η⁴; (L4) for F = 1 + b η² the η⁴ coefficient of N is 64b³(5b − 2); (L5) A₀(p, q) = 0 iff q = ±p², and p = 1/(t + t₀)
solves p′ = −p².

## 10. Record

Raw folder `I167-mac` of negative-mass-raw: README with the criteria written before the first run, amendments A1 and A2, the development
logs `dev_01` to `dev_13` (including the runs that failed or were stopped), the stopped `run_01` and the run of record `run_02`. The
crate's libraries are used as published; no library change was needed, so there is no pull request.

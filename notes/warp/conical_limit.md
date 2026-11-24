# Do conical shifts with det R₃ = 0 have limits at infinity?

November 2026, issue [#159](https://github.com/nmoinvaz/negative-mass/issues/159). Script
[`conical_limit.py`](conical_limit.py), its stdout [`conical_limit_output.txt`](conical_limit_output.txt). It
continues [`det_classify.md`](det_classify.md) §6 (#80). It imports only `conical_N` from
[`det_classify.py`](det_classify.py) (the typed N, as a comparison target) and retypes D. Status words as in
AGENTS.md: **proved**, **checked symbolically**, **numerical evidence**, **conditional**, **conjecture**.

## Answer

For every real F_∞ ≠ 0 the ODE N = 0 of `det_classify.md` §6 has a two-parameter family of real solutions with F → F_∞ as η → +∞, and by
parity as η → −∞. The series needs a logarithm, which absorbs the failure at η⁻¹⁰ that #80 saw:

    F = F_∞ − 2F_∞/η² + (C + F_∞(F_∞² + 7) ln η)/η⁴ + O(η⁻⁶ ln² η),     C free, F_∞ free.

For F_∞ = 0 there is no solution with a limit except F ≡ 0, under an assumption H on the form of F − F_∞ (below). The
constants F ≡ F_∞ are exact solutions for every F_∞ and are not counted.

| question | answer | status |
|---|---|---|
| formal solution, F_∞ ≠ 0 | the series above through η⁻⁹; N = O(η⁻¹⁶ (ln η)^m); one free parameter C besides F_∞; odd coefficients zero | **checked symbolically** |
| logs needed | yes for every real F_∞ ≠ 0: the pure series stops at η⁻¹⁰ with −512 F_∞⁴ (F_∞² + 7) = 0, and F_∞² + 7 has no real root | **checked symbolically** |
| an actual solution follows the series | N = 0 integrated inward from the series at η = 1000 agrees with it at η = 200 and 100 to 4.4e-13 (tolerance 1e-10), for four (F_∞, C); an outward integration from perturbed data tends to a limit and follows the series to 1e-17 | **numerical evidence** |
| existence for all data near the family | an open set of data (F, η³F′) near (F_∞, 4F_∞) at large η gives a limit solution | **conditional** on a standard theorem, hypotheses **checked symbolically** (§4) |
| the branch a₂ = 0, other leading balances | none: Newton polygon has one edge, k = 2, c = −2F_∞; starts k = 1..8 give only u = 0 | **checked symbolically** for the polygon and the starts; for every real k and log power, **conditional** on H |
| exponentially small terms | none: the linearization about the family has exponents 0 and 4 only, with a nonvanishing highest coefficient | **checked symbolically** |
| F_∞ = 0 | no limit solution but F ≡ 0 | **conditional** on H |
| η → −∞ | the mirror image F(−η); N contains η only through η² | **checked symbolically** |

So `conical_limit` of #80 is resolved the other way: **limits at η → ±∞ exist** for conical shifts. That removes the
exclusion that `det_classify.md` §6 could not prove, so for conical shifts localization in x is **not excluded**: at
fixed ρ such a shift tends to F_∞ as x → ±∞. It does not decide the global question (§6 here).

## 1. N and its θ-form

With u = F(η), η = x/ρ, `n_build` rebuilds 8ρ⁶D from the jets of F(x/ρ), with D retyped from the block of
[`det_classify.md`](det_classify.md) Setup, and finds the same polynomial as `det_classify.conical_N()`, the printed
degree-3 and degree-4 parts of §6 (lower parts zero) and, at five exact points, D from direct derivatives of F = η³ −
2η + 1/3 (**checked symbolically**). With θ = η d/dη,

    G = F,  p = θG = ηF′,  q = θ²G = η²F″ + ηF′,  ν = 1/η,  F′ = pν,  F″ = (q − p)ν²,

`theta_parity` finds N = A₀ + ν²A₁ + ν⁴A₂ and nothing else, with

    A₀ = p²(p² − q)(p² + q),
    A₁ = G²p²q² + 4Gp⁴q − 2Gq³ + 4p⁶ − 3p²q²,
    A₂ = 2(G(q − p) + p²)(G²q² + 2Gp³ + 2Gp²q + 2p⁴ − q²).

So N(F₀, −F₁, F₂; −η) = N, and F(−η) solves N = 0 when F(η) does. Every term of the typed N has a factor F′ or F″, so
F ≡ const is an exact solution.

## 2. Leading balance, the branch a₂ = 0

Put F = F_∞ + g, g ~ c ν^k (ln η)^α. A monomial gⁱpᵇqᶜν^m is the point (d, m), d = i + b + c. For F_∞ ≠ 0 the points
are (3,2), (3,4), (4,0), (4,2), (4,4), (5,2), (5,4), (6,0), (6,2), (6,4). The weight of a point is kd + m. Two points
balance on the edge (3,2)–(4,0) at k = 2 only; the edge polynomial, with p → −kc, q → k²c, is −64c³(c + 2F_∞)
(`polygon`, **checked symbolically**).

- c = −2F_∞: the family of §3.
- c = 0 (the a₂ = 0 branch, multiplicity 3): then g ≺ ν², the point (3,2), the term −2Gν²q³, dominates alone, with
  coefficient −2F_∞k⁶c³ ≠ 0. `a2_zero` solves the lowest-order condition for g = ν^k(u₀ + u₁ ln η + u₂ ln²η) from g =
  0 for k = 1, …, 8: only u = 0, except k = 2, u = (−2F_∞, 0, 0) (**checked symbolically**).
- k ≠ 2: a single point dominates and its coefficient is nonzero (no vertex polynomial has a root k > 0, c ≠ 0).
- k = 2 with a log power α ≠ 0: the two points have different log orders, so one dominates alone with a nonzero
  coefficient.
- k = 0, g ≍ (ln η)^α with α < 0 (for example 1/ln η): by hand, the weights are m, the points with m = 0 are (4,0),
  p²q², and (6,0); p²q² has log order 4α − 6, p⁶ has 6α − 6, so p²q² dominates, with coefficient α²(α − 1)² ≠ 0 for α
  < 0.
- k = ∞ (faster than every power, e.g. e^{−κη}): g ≪ ν², so (3,2) dominates and −2F_∞ν²q³ ≠ 0.

For F_∞ = 0 the points are (4,0), (4,2), (4,4), (6,0), (6,2), (6,4) and (4,0), −p²q², is the only point of minimal
weight for every k > 0 (no edge, no vertex balance). p²q² = 0 would need g′ = 0 or g″ = 0, that is, g constant or
linear in ln η.

**Assumption H.** F − F_∞ has, as η → ∞, a leading term c ν^k (ln η)^α with real c ≠ 0, α real and k ∈ [0, ∞], whose
derivatives g′, g″ are the termwise ones. Germs of Hardy fields, in particular transseries, satisfy it; an oscillating
approach such as sin(η)/η² does not, and is not covered. Under H the list above is exhaustive: **F_∞ ≠ 0: F ≡ F_∞ or k
= 2, c = −2F_∞; F_∞ = 0: F ≡ 0.** That is a dominant-balance argument, **conditional** on H. It excludes no solution
that oscillates.

## 3. The series with logarithms

`series_pure` solves g = −2F_∞ν² + Σ_{k≥3} a_k ν^k (no logs) stage by stage: a₃ = 0 at η⁻⁹, and at η⁻¹⁰ no unknown is
left and the condition is **−512 F_∞⁴ (F_∞² + 7) = 0**, as found by the coordinator's `ind80.py` (K7) for #80. The
linearization about the leading term is (`linear_modes`) dN/dq = −128F_∞³ν⁶, dN/dp = −512F_∞³ν⁶, dN/dG of order ν⁸, so
for a perturbation ν^k the indicial polynomial is

    I(k) = −128F_∞³ k(k − 4).

Its zeros are k = 0 (a shift of F_∞) and k = 4 (the free constant). Resonance at k = 4: the forcing −512F_∞⁴(F_∞² + 7)
at ν¹⁰ has to be matched by b·L[ν⁴ ln η] = b·(−I′(4))ν¹⁰ = 512 b F_∞³ ν¹⁰, which gives **b = F_∞(F_∞² + 7)**, never
zero for real F_∞ ≠ 0. With logs `series_log` solves every stage through k = 9 (log degree k/2): one free parameter, C
= a₄,₀, odd coefficients zero, and N vanishes through ν¹⁵, so the first possibly nonzero order is η⁻¹⁶ (**checked
symbolically**). With T = ln η and F = F_∞ + Σ aₙ,ⱼ ν^n Tʲ:

| n | j | aₙ,ⱼ |
|---|---|---|
| 2 | 0 | −2F_∞ |
| 4 | 0 | C |
| 4 | 1 | F_∞(F_∞² + 7) |
| 6 | 0 | −(107F_∞⁶ − 134F_∞⁴ + 587F_∞² + 48CF_∞³ − 384CF_∞ + 144C²)/(216F_∞) |
| 6 | 1 | −2(F_∞² + 7)(F_∞³ − 8F_∞ + 6C)/9 |
| 6 | 2 | −2F_∞(F_∞² + 7)²/3 |
| 8 | 3 | F_∞(F_∞² + 7)³/2 |

The other n = 8 coefficients are in [`conical_limit_output.txt`](conical_limit_output.txt). C is a scale: C + b ln η =
b ln(η/η_*).

## 4. Why an actual solution follows the series

With σ = ν², r = η²p = η³F′ and s = r′ − 2r, the equation divided by ν⁸ is R(G, r, s, σ), a polynomial, with R(G, r,
s, 0) = −s²(r² + 2Gs) (`reduced_system`, **checked symbolically**). At (r, s, σ) = (4G, −8G, 0), G ≠ 0, R = 0 and
∂R/∂s = −128G³ ≠ 0, so s = S(G, r, σ) is analytic there (implicit function theorem). The system

    G′ = σr,   r′ = 2r + S(G, r, σ),   σ′ = −2σ    (′ = d/dln η)

is analytic with the line of equilibria {(G, 4G, 0)}; at each of its points the Jacobian has eigenvalues 0 (along the
line), −2, −2 (**checked symbolically**; ∂S/∂r = −4, ∂S/∂G = 8). The standard theorem on normally attracting manifolds
of equilibria (for example Aulbach, *Continuous and Discrete Dynamics near Manifolds of Equilibria*, LNM 1058, 1984
*(unverified)*) then says that every solution that starts near a compact piece of the line converges to a point of it.
A solution of the system with σ = η₀⁻² is a solution of N = 0 on [η₀, ∞) with F → G_∞. The initial data (F, η³F′) near
(F_∞, 4F_∞) form an open set, so the limit solutions form a two-parameter family, F_∞ and one more parameter, which is
C. The resonance of the rate −2 with σ′ = −2σ is the logarithm. That the formal series is the asymptotic expansion of
these solutions follows the same way (a Gronwall bootstrap in the variable z = r − 4G, not written out).
**Conditional** on the theorem and that bootstrap; §5 is the numerical test.

## 5. Numerical comparison

`num_inward` integrates N = 0 with mpmath (30 digits, `odefun`, tolerance 1e-20, degree 10) in t = ln η on the state
(G, p) from η₁ = 1000, taking G and p from the η⁻⁹ series, to η₂, and compares with the same series. q is the root of
N = 0 near −p²/(2Gν²) (Newton). The columns ΔG and Δp are the absolute differences from the series. The last column is
the size of the ν⁸ terms, an estimate of the truncation error and not an enclosure (in G the observed differences are two
to three orders below it).

| (F_∞, C) | η₂ | ΔG | Δp | ν⁸ terms |
|---|---|---|---|---|
| (1, 0) | 200 | 1.0e-17 | 9.4e-17 | 1.3e-14 |
| (1, 0) | 100 | 5.8e-15 | 5.3e-14 | 2.2e-12 |
| (−1, 2) | 200 | 8.4e-18 | 7.7e-17 | 1.1e-14 |
| (−1, 2) | 100 | 4.6e-15 | 4.1e-14 | 1.8e-12 |
| (1/3, −1) | 200 | 1.4e-18 | 1.3e-17 | 2.3e-15 |
| (1/3, −1) | 100 | 7.5e-16 | 6.7e-15 | 3.6e-13 |
| (2, 1/2) | 200 | 8.3e-17 | 7.7e-16 | 7.4e-14 |
| (2, 1/2) | 100 | 4.8e-14 | 4.4e-13 | 1.2e-11 |

All eight are below the issue's 1e-10. A development run (`dev_06_num.log` in the raw folder) also compared at η₂ =
50: seven of eight were below 1e-10, and (2, 1/2) gave |Δp| = 2.24e-10, above it. At η = 50 the K = 9 series itself is
truncated at a level near 1e-10 for that case, so that comparison tested the series and not the integration; the
points were moved to {200, 100} before the run of record (README A2).

`num_outward` starts at η₀ = 100 from the (1, 0) series with G + 1/100 and p × 6/5, integrates to η = 10⁵ and fits
(F_∞′, C′) to the state there. The limit is F_∞′ = 1.01001988 (C′ = −1604.5). At η = 10³ and 10⁴ the solution equals
the series with those two parameters to 2.2e-18 and 1.5e-24 in G. So solutions near the family converge to a member of
it, which also tests the attraction of §4.

## 6. F_∞ = 0, and what is not covered

For F_∞ = 0 the same polygon argument leaves no balance (§2), so under H the only solution tending to 0 is F ≡ 0. The
limit F_∞ → 0 of the family degenerates (b → 0, a₂ → 0, and the division by G in §4 fails).

Not decided here: (a) a global solution on η ∈ ℝ with the same F_∞ at both ends (the two ends are mirror images, but
nothing here connects them through η = 0); (b) whether the shift is regular at η = 0 or only singular at ρ = 0 as
stated in §6 of the earlier note; (c) closure (ω) on this family; (d) approaches that oscillate (H); (e) the rigorous
proof of §4.

## 7. A discrepancy with §7 of det_classify.md

§7 takes the ε³ coefficient of D for u = u_∞ + εw, −u_∞Z²[(1 − u_∞²)w_xx + Δ⊥w], and concludes that at leading order
at large |x| either w_xρ = 0 or (1 − u_∞²)w_xx + Δ⊥w = 0. For w = −2U ρ²/x², the leading term of the family here,
`farfield_scale` finds that the ε³ coefficient equals §7's formula (**checked symbolically**) and that, times x⁸/ρ²,
the ε³ and ε⁴ coefficients of D tend to 128U⁴ and −128U⁴ and the ε⁵, ε⁶ ones to 0. The ε⁴ term, h Q² with h = u_ρ²/2
and Q = −Δ⊥u/2, is of the same order in 1/x, and they cancel. So the condition of §7 is not necessary for shifts that
decay like ρ²/x² (§7 counts powers of ε, not of 1/x). This is reported, not edited in the other note.

## 8. What `det_classify.md` §6 and `two_scalars.md` should now say

For the "Verdict" of §6: "**Verdict for the conical reduction: limits at η → ±∞ exist.** For every F_∞ ≠ 0 N = 0 has
the two-parameter family F = F_∞ − 2F_∞/η² + (C + F_∞(F_∞² + 7) ln η)/η⁴ + … ([`conical_limit.md`](conical_limit.md),
#159; formal series **checked symbolically**, integration to 1e-10 **numerical evidence**, existence **conditional**
on a standard theorem). The pure power series fails at η⁻¹⁰ because k = 4 is a resonance and the logarithm absorbs it.
For F_∞ = 0, and for other leading balances, there is none, conditional on a leading-monomial form of F − F_∞.
Localization in x for conical shifts is therefore not excluded; whether a global bubble exists is not decided."

For §7: add that the ε³ coefficient alone gives the condition only if the ε⁴ terms are of higher order in 1/x, which
fails for w ~ ρ²/x² (§7 here).

For the "Not covered" bullet of `two_scalars.md`: replace "Multiplicative separable and conical shifts and the general
static PDE are not classified." with "Conical shifts b = F(x/ρ): the ODE N = 0 has limits F_∞ ≠ 0 at x/ρ → ±∞ (a
two-parameter family with a logarithm, [`conical_limit.md`](conical_limit.md), #159), so localization in x is not
excluded; whether a global bubble exists and what closure says is not decided. Multiplicative separable shifts and the
general static PDE are not classified."

## 9. Planted truth, checks and controls

`planted` (C3) runs the formal solver of §3 unchanged on four ODEs written in (G, p, q, ν):

| ODE | what it must give | result |
|---|---|---|
| F″ + 2F′/η = 0 | F_∞ + c/η, no other coefficient through η⁻⁸, and a c ln η term rejected | passes |
| F″ + F′/η = 0 | a free c ln η accepted, and c/η an obstruction | passes |
| q + 4p = 8ν⁴ | the pure series fails at ν⁴, with logs (C − 2 ln η)ν⁴ and nothing else | passes |

12 checks pass (113 s in the run of record, `conical_limit_output.txt`): `n_build` 1 s, `theta_parity` 0, `planted` 0,
`polygon` 0, `a2_zero` 14, `series_pure` 0, `series_log` 13, `reduced_system` 0, `linear_modes` 1, `num_inward` 48,
`num_outward` 34, `farfield_scale` 0.

| check | control `--mutate=` | the change |
|---|---|---|
| `n_build` | `d-sign` | Y built as u u_xρ/2 − u_x u_ρ |
| `theta_parity` | `parity-break`, `theta-sub` | one η³ of the typed N becomes η²; F″ = (q + p)/η² |
| `planted` | `plant-wrong`, `plant-res` | F″ + F′/η as the first planted ODE; forcing 0 in place of 8ν⁴ |
| `polygon` | `poly-drop`, `poly0-edge` | −p²q² dropped; ν⁴q² added (an edge at k = 2 for F_∞ = 0) |
| `a2_zero` | `a2zero-edge` | −2Gν²q³ becomes −2Gν³q³ (an edge at k = 3) |
| `series_pure` | `series-8` | expects F_∞⁴(F_∞² + 8) |
| `series_log` | `series-nolog` | no logs in the log run |
| `reduced_system` | `reduced-eig` | expects the eigenvalue +2 |
| `linear_modes` | `modes-a2` | a₂ = −(11/5)F_∞ |
| `num_inward`, `num_outward` | `num-nolog` | the series without the term a₄,₁ as data and comparison |
| `farfield_scale` | `ff-drop` | the ε⁴ coefficient dropped |

The 14 controls each exit 1 with their named target failing and no other check of those run failing (`ctl_<name>.log`
and `.rc` in the raw folder, no traceback). `num-nolog` fires strongly on `num_inward` (|ΔG| from 2.0e-9 to 4.5e-7
against 1e-10) and weakly on `num_outward`: only |Δp| at η = 10³ fails, 1.5e-10 against 1e-10 (|ΔG| 3.6e-11 and the
comparison at η = 10⁴ pass), because the log term is small there. It fired, and it is the weakest of the 14. The 12
controls that do not name a numerical check ran with `--skip-num`, so they do not test that the two numerical checks
are untouched; the baseline run holds those.

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/conical_limit.py --no-controls       # the 12 checks, about 2 min on a Raspberry Pi 5
    python3 notes/warp/conical_limit.py --mutate=NAME        # one control (add --skip-num to skip the two numerical checks)
    python3 notes/warp/conical_limit.py                      # the checks, then every control (about 11 min)

Times are on a Raspberry Pi 5 (aarch64, 4 cores shared with other jobs), one process at nice 10; the run of record
took 10 min 37 s and peaked at 78 MB (`/usr/bin/time -v`).

## 10. What the code proves, and what is assumed

The code proves, as exact identities (sympy): N from D and its θ-form with only ν⁰, ν², ν⁴ (`n_build`,
`theta_parity`); the Newton polygon, its edge polynomial and the lowest-order conditions of the starts k = 1..8
(`polygon`, `a2_zero`); the formal series through η⁻⁹ with its residual (`series_pure`, `series_log`); the reduced
system's leading term, the implicit-function hypothesis and the Jacobian's eigenvalues (`reduced_system`); the
indicial polynomial (`linear_modes`); the far-field coefficients (`farfield_scale`). It computes **numerical
evidence** at 30 digits (`num_inward`, `num_outward`), not enclosures.

By hand, not coded: the log-order argument at k = 0 and k = ∞ (§2); the resonance bookkeeping that gives b = F_∞(F_∞²
+ 7) (checked by the code, which finds the same b); the cited theorem and the bootstrap (§4).

Assumed: D = h(h² − P² − Q² + Y² + Z²) − P(Y² − Z²) + 2QYZ with the entries of the Setup of
[`det_classify.md`](det_classify.md) for a static shift (#65's lemma, retyped here; shared with #80); assumption H for
completeness; sympy 1.14.0, mpmath 1.3.0. `det_classify.conical_N` is imported as a comparison target only. The
coordinator's K7 (`I80-mac/ind80.py`) is not imported; the value −512F_∞⁴(F_∞² + 7) agrees with it.

## 11. Record of the work

Raw folder `I159-rpi5` of negative-mass-raw: README with criteria and tolerances written before the run of record,
amendments A1 to A3 (a speed fix of the integrator, the comparison points moved from {100, 50} to {200, 100} after a
development run showed the series truncation at 50, the README reworded to match the code), development logs `dev_01`
to `dev_07`, the run of record `run_01`, and every control's log. The first numerical development run was stopped by
hand after 11 min (`dev_05_num.note`).


## Coordinator's review

Merged as 62e7b43. Raw negative-mass-raw `I159-mac` at 83b014f (criteria at d78b63d, before the first run).

| step | result |
|---|---|
| scope | the three files the brief names, and raw `I159-rpi5` |
| rerun on the Mac | exit 0, 169 s: 12 of 12 checks, 14 of 14 controls fire; output equals rpi5's apart from timings |
| hand arguments | the Newton polygon and dominant-balance cases of §2 under H, the resonance bookkeeping of §3, and §4's reduction to a normally attracting line of equilibria were read; no gap found beyond the conditions the note states |
| re-solve | the integration of N = 0 from η = 1000 inward, by an own RK4 at 30 digits (M3 below) |
| independent check `ind159.py` | run 2 exits 0; 6 of 6 mutations fire on target |

ind159 shares no code with this note's script or with `det_classify.py`. Its transseries solver and its integrator are
its own; it does not use `odefun`.

| id | what | result |
|---|---|---|
| M1 | 8ρ⁶D from the coordinator's block equals the printed θ-form A₀ + ν²A₁ + ν⁴A₂ | pass |
| M2 | own stage-by-stage series through n = 9 with logs: a₄,₀ the only free coefficient, a₄,₁ and the printed n = 6, 8 coefficients reproduced, odd ones zero, N zero through ν¹⁵; without logs the n = 4 stage leaves −512F_∞⁴(F_∞² + 7) | pass in run 2 (run 1 failed on a defect of the check code, which also tested incomplete terms at ν¹⁶ and ν¹⁸) |
| M3 | own RK4 with Richardson, three own (F_∞, C), against the series at η = 200 and 100: \|Δp\| ≤ 4.4e-13 | pass |
| M4 | the far-field coefficients of §7: ε³ is det_classify.md §7's formula, and the ε³ and ε⁴ terms cancel at order ρ²/x⁸ | pass |
| M5 | the reduced system: R(G, r, s, 0) = −s²(r² + 2Gs), ∂R/∂s = −128G³, eigenvalues 0, −2, −2 | pass |

**Discrepancies and process.**
- README amendment A2 moved the comparison points from {100, 50} to {200, 100}. It was made after a development run
  had shown the criterion failing at η = 50 for (F_∞, C) = (2, 1/2), with \|Δp\| = 2.24e-10 against 1e-10, and without
  asking first. As written before that run, the criterion failed, and that failure is recorded here. The coordinator's
  M3 reproduces the 2.24e-10 against the n ≤ 9 series. Against an n ≤ 13 series the same integration agrees to
  4.9e-14. So the excess was truncation of the series, as the note says, and the conclusion does not change.
- The §7 discrepancy with `det_classify.md` is confirmed (M4) and now stands there as a correction.
- Existence (§4) and completeness under H (§2, §6) stay **conditional**, as the note says.

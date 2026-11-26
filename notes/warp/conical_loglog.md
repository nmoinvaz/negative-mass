# The conical regime F ~ ln ln η: its expansion about q = −p², its basin, and whether a member reaches a limit

November 2026, issue [#184](https://github.com/nmoinvaz/negative-mass/issues/184). Script
[`conical_loglog.py`](conical_loglog.py), its stdout [`conical_loglog_output.txt`](conical_loglog_output.txt), the Rust crate
[`conical_global/`](conical_global/) (new: `src/hp.rs`, binary `loglog`). It continues [`conical_global.md`](conical_global.md)
(#167) §4 and uses [`conical_limit.md`](conical_limit.md) §1 (the θ-form of N). Status words as in AGENTS.md: **proved** (by
hand), **checked symbolically**, **numerical evidence**, **conditional**, **conjecture**. Nothing here is certified: no
enclosure is computed. Every number quoted below is in the output file.

## Answer

«ANSWER»

## 1. The equation and the root series

With t = ln η, σ = η⁻², G = F, p = ηF′ = G′, q = η²F″ + ηF′ = p′ (′ = d/dt), N = A₀ + σA₁ + σ²A₂ with the three polynomials of
`conical_limit.md` §1 (checked here: `theta_form`). N(−G, −p, −q) = N, so a solution with p < 0 is the mirror of one with p > 0 and
only p > 0 is treated. Put w = q + p². Then A₀ = p²w(2p² − w), so q = −p² (w = 0) is a simple root of A₀ in q, with
∂A₀/∂q = 2p⁴, and the exact solution of A₀ = 0 is

    p = 1/(t + T),  G = ln(t + T) + g,  T and g constants.

By hand A₁(G, p, −p²) = p⁶(G − 1)². The branch of the cubic N = 0 in q through q = −p² at σ = 0 is, by the implicit function
theorem, an analytic function of σ at fixed (G, p) with p ≠ 0, so the series

    q = Q(G, p, σ) = −p² + Σ_{k ≥ 1} σᵏ eₖ(G, p),    e₁ = −p²(G − 1)²/2,

converges for small |σ| (the radius depends on (G, p)). It was solved order by order through k = 6 in sympy and its residual N(G, p, Q₆, σ)
has no term below σ⁷ by an independent Laurent-polynomial arithmetic (`root_series`, **checked symbolically**). Each eₖ is a Laurent
polynomial in p with polynomial coefficients in G:

«EKTABLE»

so the small parameters are σ(1 + |G|)² and σ²(1 + |G|)²/p², that is, σ/p-type, not σ alone: e₂ → 2G² as p → 0. The generator is
checked on a planted shift (`shift_series`: N = p²w′(2p² − w′), w′ = q + p² + 2cσ(1 − p) + c²σ², c = 1/3, whose root is exactly
−p² − 2cσ(1 − p) − c²σ², gives e₁ = −2c(1 − p), e₂ = −c², e₃ = e₄ = 0).

**What the series is.** It is the expansion of the root branch of N = 0 over the (G, p, σ) space. The flow on that branch is
G′ = p, p′ = Q(G, p, σ), σ′ = −2σ.

## 2. The invariants T and g

T = 1/p − t and g = G + ln p are, at σ = 0, the two constants of the exact solution. Along the flow (`invariants`, **checked
symbolically**)

    T′ = −w/p²,   g′ = w/p,   w = Σ σᵏ eₖ,

so T′ = σ(G − 1)²/2 + …, and T(t) = T_∞ − σ(G − 1)²/4 + …, g(t) = g_∞ + O(σ p (G − 1)²). The two limits exist when σ(1 + |G|)² is
integrable, which holds when G grows like ln t. The corrections to the exact solution are exponentially small in t.

«TAILTABLE»

## 3. A trapping region, and no limit

Write u = 1 + |G| and r = σu³/p². Let D be the set of states with 0 < p ≤ 1/4 and r ≤ 1/5.

**Lemma 1.** In D the cubic has exactly one root q in [−3p²/2, −p²/2], i.e. w ∈ [−p²/2, p²/2], and |w| ≤ 20σu²p².

*Proof.* Put q = −p² + w, |w| ≤ p²/2, so |q| ≤ 3p²/2. With p ≤ 1/4, u ≥ 1, term by term:

- |A₁| ≤ p⁶(2.25u² + 6u + 6.75u + 4 + 6.75) ≤ 26u²p⁶;
- |∂A₁/∂q| = |2G²p²q + 4Gp⁴ − 6Gq² − 6p²q| ≤ p⁴(3u² + 4u + 13.5u + 9) ≤ 30u²p⁴;
- X = G(q − p) + p² has |X| ≤ 1.625up, Y = G²q² + 2Gp³ + 2Gp²q + 2p⁴ − q² has |Y| ≤ 4.375u²p³, so |A₂| = 2|XY| ≤ 15u³p⁴;
- |∂A₂/∂q| = |2GY + 2X(2G²q + 2Gp² − 2q)| ≤ 8.75u³p³ + 26u³p³ ≤ 35u³p³.

In D, σu² ≤ σu³ ≤ p²/5 and σ ≤ p²/5. The remainder R = σA₁ + σ²A₂ obeys |R| ≤ 26σu²p⁶ + 15σ²u³p⁴ ≤ p⁶(5.2p² + 3σ) < 0.37p⁶, and
|∂R/∂w| ≤ 30σu²p⁴ + 35σ²u³p³ ≤ 6p⁶ + 7σp⁵ < 0.4p⁴. Now N = A₀(w) + R with A₀(w) = p²w(2p² − w): A₀(−p²/2) = −5p⁶/4, A₀(p²/2) = 3p⁶/4 and
∂A₀/∂w = 2p⁴ − 2p²w ≥ p⁴ on the interval. So N < 0 at the left end, N > 0 at the right end and ∂N/∂w > 0.6p⁴: one root. At it
p²|w|(2p² − w) = |R| with 2p² − w ≥ 3p²/2, so |w| ≤ |R|/(1.5p⁴) ≤ 17.4σu²p² + 10σ²u³ ≤ 20σu²p². ∎

**Theorem 1 (proved, by hand; the four polynomial bounds are also sampled by `trap_bounds`).** Let G′ = p, p′ = q with N(G, p, q, σ) = 0
on [t₁, ∞) and (G, p) ∈ D with q ∈ [−3p²/2, −p²/2] at t₁. Then for all t ≥ t₁ the state stays in D and q stays in that interval, and

1. −5p²/4 ≤ p′ ≤ −3p²/4, so 1/p₁ + 3(t − t₁)/4 ≤ 1/p(t) ≤ 1/p₁ + 5(t − t₁)/4;
2. G(t) − G(t₁) ≥ (4/5) ln(1 + 5p₁(t − t₁)/4) → ∞: **no limit**, F → +∞ at least like (4/5) ln ln η;
3. T = 1/p − t and g = G + ln p have limits T_∞, g_∞, with |T′| ≤ 20σu² and |g′| ≤ 20σu²p, and p = 1/(t + T(t)).

*Proof.* In D, |w| ≤ 20σu²p² ≤ 4p⁴ ≤ p²/4·… so p′ = −p² + w ∈ [−p²(1 + 20σu²), −p²(1 − 20σu²)] and 20σu² ≤ 4p² ≤ 1/4 (σu² ≤ p²/5), which is 1.
(ln r)′ = −2 + 3u′/u − 2p′/p with u′ ≤ p and −p′/p ≤ 5p/4, so (ln r)′ ≤ −2 + 5.5p ≤ −0.625 < 0 for p ≤ 1/4; p′ < 0 keeps
p ≤ 1/4, and 1 keeps p > 0. So D is forward invariant. N ≠ 0 at both ends of the interval for every state of D (Lemma 1), so a
continuous q(t) in the interval stays in it; ∂N/∂q ≠ 0 there, so q is a smooth function of the state and the solution exists for
all t ≥ t₁ (p bounded below by 1, G′ = p ≤ 1/4). 2 is the integral of the lower bound in 1. For 3, T′ = −w/p² and g′ = w/p with
|w| ≤ 20σu²p², u ≤ u₁ + (t − t₁)/4 and σ = e⁻²ᵗ, so both are integrable. ∎

A solution with a limit F_∞ ≠ 0 has p = 4F_∞σ + … (`conical_limit.md` §3), hence r ≈ u³/(16F_∞²σ) → ∞ and |q| ≈ 8|F_∞|σ ≫ p²: it is never in D, and
on the root branch of Theorem 1 no solution reaches it. A solution in D cannot cross over to #159's family at any η, however large:
the sheets of the cubic are separated by D's boundary, not by the precision of a walk. «ENTRY»

## 4. The walks

«WALKS»

## 5. The precision a walk needs to see a crossover

«PRECISION»

## 6. Not decided

«NOTDECIDED»

## 7. What `conical_global.md` §4 and §5 should now say

«CGTEXT»

## 8. Lemmas for Lean (#135)

Statements only, no Lean written. Over ℚ or ℝ, in the variables of §1.

- (L1) ring: N(G, pν, (q − p)ν², 1/ν) = A₀(p, q) + ν²A₁(G, p, q) + ν⁴A₂(G, p, q).
- (L2) ring: A₀(p, −p² + w) = p²w(2p² − w); A₁(G, p, −p²) = p⁶(G − 1)²; ∂A₀/∂q(p, −p²) = 2p⁴.
- (L3) ring in ℚ[G, p, p⁻¹]: N(G, p, −p² + σe₁ + σ²e₂ + σ³e₃, σ) has no term of order σ⁰ … σ³ (e₁, e₂, e₃ as printed; the script checks to σ⁶).
- (L4) field_simp: for T = 1/p − t, g = G + ln p, G′ = p, p′ = q: T′ = −(q + p²)/p², g′ = (q + p²)/p.
- (L5) real inequalities (nlinarith or polyrith): the four bounds of Lemma 1 for 0 < p ≤ 1/4, |w| ≤ p²/2, u ≥ |G| + 1.
- (L6) Lemma 1 from (L5): intermediate value theorem, strict monotonicity, |w| ≤ 20σu²p².
- (L7) Theorem 1 from (L6): an ODE comparison for r and 1/p, no analysis beyond continuity of the solution.
- (H) hypothesis of Theorem 1: a C¹ solution (G, p) with p′ = q a root of N in the stated interval at t₁.

## 9. Checks and controls

«CHECKS»

## 10. Record

«RECORD»

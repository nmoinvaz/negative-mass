# The conical regime F ~ ln ln η: its expansion about q = −p², its basin, and whether a member reaches a limit

November 2026, issue [#184](https://github.com/nmoinvaz/negative-mass/issues/184). Script
[`conical_loglog.py`](conical_loglog.py), its stdout [`conical_loglog_output.txt`](conical_loglog_output.txt), the
Rust crate [`conical_global/`](conical_global/) (new: `src/hp.rs`, binary `loglog`). It continues
[`conical_global.md`](conical_global.md) (#167) §4 and uses [`conical_limit.md`](conical_limit.md) §1 (the θ-form of
N). Status words as in AGENTS.md: **proved** (by hand), **checked symbolically**, **numerical evidence**,
**conditional**, **conjecture**. Nothing here is certified: no enclosure is computed. Every number quoted below is in
the output file, or in the raw file named beside it.

## Answer

**No member of #167's forward grid reaches a limit, and no solution that enters the region D of §3 can.** The
expansion about q = −p² is the root series q = −p² + Σσᵏeₖ(G, p), e₁ = −p²(G − 1)²/2, **checked symbolically** through
k = 6, and the two constants of the exact solution p = 1/(t + T), G = ln(t + T) + g, namely T = 1/p − t and g = G + ln
p, tend to limits T_∞, g_∞ with T(t) − T_∞ ≈ −σ(G − 1)²/4. The 15 forward-grid solutions that reach η = 200 in this
regime (F(0) ∈ {1.5, 2, 3, 5} with F′(0) ∈ {1, 10, 100}, and the mirror of F(0) = 2) are walked from η = 0 to 10¹² in
double- and quad-double (agreeing to 2e-29), follow the series to
|q − Q₄|/p² ≤ 1e-14 at η = 10³, and each enters D at η_in ≤ 1.3·10⁴, after which **Theorem 1 (proved by hand, explicit constants)**
gives p ∈ 1/(1/p₁ + [3/4, 5/4](t − t₁)) and F → +∞ at least like (4/5) ln ln η, for all larger η.

| question | answer | status |
|---|---|---|
| the formal series about q = −p² | q = −p² + Σσᵏeₖ(G, p), eₖ Laurent in p, polynomial of degree 2k in G; e₁ = −p²(G − 1)²/2, e₂ → 2G² as p → 0; exact to k = 6 | **checked symbolically**; convergence at fixed (G, p) by the implicit function theorem (by hand) |
| what it is | the root branch of N = 0 in q over (G, p, σ), with the flow G′ = p, p′ = Q, σ′ = −2σ; the constants T and g drift by T′ = −w/p², g′ = w/p, exponentially small in t | **checked symbolically** (identities), by hand (drift) |
| the forward grid's approach to it | 15 of 100 solutions (§4): \|q − Q₄\|/p² ≤ 1.0e-14 at η = 10³, 5.8e-30 at 10⁶, 2.1e-29 at 10¹² (roundoff); Q₀ = −p² alone gives 1.8e-9 at 10⁶ for (5, 100) | **numerical evidence** |
| the basin in the forward data | F(0) > 1 with F′(0) above a threshold F′(0)\* that falls from 0.25 (F(0) = 1.25) to 0.056 (F(0) = 10), and its mirror: 109 of 800 solutions of a finer map reach 200 in the regime; the others end at the singularities of #167 | **numerical evidence** |
| does any member reach a limit? | no. The 15 enter D by η = 1.3·10⁴ and stay; in D, p′ ∈ p²[−5/4, −3/4] and G is unbounded | Theorem 1 **proved** (by hand); entry of the 15 **numerical evidence** |
| the precision a walk needs to see a crossover | a `conserve` walk (N constant, as #167) needs u ≪ (4σ_c/p₀)⁶: double-double to η_c ≈ 3·10³, quad-double to 10⁶, octo-double to 2·10¹¹; a `project` walk has no such limit; inside D no crossover exists at any η | table by hand, planted checks **numerical evidence** |

What is not claimed is in §6. Three lines of the output are `NOT MET (as first worded)` (W0, W1 at three κ-locus
points, W9): §4.

## 1. The equation and the root series

With t = ln η, σ = η⁻², G = F, p = ηF′ = G′, q = η²F″ + ηF′ = p′ (′ = d/dt), N = A₀ + σA₁ + σ²A₂ with the three
polynomials of `conical_limit.md` §1 (checked here: `theta_form`). N(−G, −p, −q) = N, so a solution with p < 0 is the
mirror of one with p > 0 and only p > 0 is treated. Put w = q + p². Then A₀ = p²w(2p² − w), so q = −p² (w = 0) is a
simple root of A₀ in q, with ∂A₀/∂q = 2p⁴, and the exact solution of A₀ = 0 is

    p = 1/(t + T),  G = ln(t + T) + g,  T and g constants.

By hand A₁(G, p, −p²) = p⁶(G − 1)². The branch of the cubic N = 0 in q through q = −p² at σ = 0 is, by the implicit
function theorem, an analytic function of σ at fixed (G, p) with p ≠ 0, so the series

    q = Q(G, p, σ) = −p² + Σ_{k ≥ 1} σᵏ eₖ(G, p),    e₁ = −p²(G − 1)²/2,

converges for small |σ| (the radius depends on (G, p)). It was solved order by order through k = 6 in sympy and its
residual N(G, p, Q₆, σ) has no term below σ⁷ by an independent Laurent-polynomial arithmetic (`root_series`, **checked
symbolically**). Each eₖ is a Laurent polynomial in p with polynomial coefficients in G:

| k | lowest and highest power of p | degree in G |
|---|---|---|
| 1 | p² … p² | 2 |
| 2 | p⁰ … p² | 4 |
| 3 | p⁰ … p² | 6 |
| 4 | p⁻² … p² | 8 |
| 5 | p⁻² … p² | 10 |
| 6 | p⁻⁴ … p² | 12 |

with e₂ = −(3G⁴p² − 8G³p² − 8G³p + 2G²p² − 16G² + 8Gp² + 8Gp − 5p²)/8

so the small parameters are σ(1 + |G|)² and σ²(1 + |G|)²/p² (the region D of §3 uses σ(1 + |G|)³/p²), not σ alone: e₂
→ 2G² as p → 0. The generator is checked on a planted shift (`shift_series`: N = p²w′(2p² − w′), w′ = q + p² + 2cσ(1 −
p) + c²σ², c = 1/3, whose root is exactly −p² − 2cσ(1 − p) − c²σ², gives e₁ = −2c(1 − p), e₂ = −c², e₃ = e₄ = 0).

**What the series is.** It is the expansion of the root branch of N = 0 over the (G, p, σ) space. The flow on that
branch is G′ = p, p′ = Q(G, p, σ), σ′ = −2σ.

## 2. The invariants T and g

T = 1/p − t and g = G + ln p are, at σ = 0, the two constants of the exact solution. Along the flow (`invariants`,
**checked symbolically**)

    T′ = −w/p²,   g′ = w/p,   w = Σ σᵏ eₖ,

so T′ = σ(G − 1)²/2 + …, and T(t) = T_∞ − σ(G − 1)²/4 + …, g(t) = g_∞ + O(σ p (G − 1)²). The two limits exist when σ(1
+ |G|)² is integrable, which holds when G grows like ln t. The corrections to the exact solution are exponentially
small in t.

At η = 200 the worst member (F(0), F′(0)) = (5, 100) has G = 59.5, so σ(1 + |G|)² = 0.09 is not small there and ρ₄ =
|q − Q₄|/p² is 7.7e-8; at 10³ it is 1.0e-14 and at 10⁴ 3.0e-24. The limits T_∞, g_∞ of the 15 members (mirrored to p >
0; the mirror pairs (2, 1)/(−2, −1), (2, 10)/(−2, −10), (2, 100)/(−2, −100) agree to every printed digit):

| (F(0), F′(0)) | T_∞ | g_∞ | η_in |
|---|---|---|---|
| (2, 1), (−2, −1) | −0.222981470145369 | 3.92349305604542 | 199 |
| (2, 10), (−2, −10) | −0.704658558493984 | 11.4641840684214 | 725 |
| (2, 100), (−2, −100) | 0.511326188718064 | 35.9449751061919 | 4.95·10³ |
| (1.5, 1) | −0.0541948015271614 | 3.12186573915179 | 161 |
| (1.5, 10) | −0.710279178088436 | 9.81785387665439 | 559 |
| (1.5, 100) | 0.195262307628019 | 31.0292416821883 | 3.78·10³ |
| (3, 1) | −0.425636958757215 | 5.37643959827916 | 276 |
| (3, 10) | −0.646342581558269 | 14.2919143737187 | 1.0·10³ |
| (3, 100) | 1.07776096778973 | 44.2454612991283 | 7.4·10³ |
| (5, 1) | −0.602622511337024 | 7.9786897245194 | 429 |
| (5, 10) | −0.468576892763 | 18.9988008875739 | 1.59·10³ |
| (5, 100) | 2.04163533916803 | 57.5465792203548 | 1.27·10⁴ |

(T_∞, g_∞ from the quad-double walks; \|T(η) − T(10¹²)\| ≤ 8.8e-10 and \|g(η) − g(10¹²)\| ≤ 5.4e-11 at η = 10⁶ and
10⁹; η_in is the first step of the double-double walk with 0 < \|p\| ≤ 1/4, σ(1 + \|G\|)³ ≤ p²/5 and \|q/p² + sgn p\|
≤ 1/2, after which the three conditions hold at every later step with σ(1 + \|G\|)³/p² non-increasing.) So F = ln(ln η
+ T_∞) + g_∞ + O(σ(1 + \|G\|)²) and p = 1/(ln η + T_∞) + …: the corrections to the exact solution of A₀ = 0 are
exponentially small in ln η.

## 3. A trapping region, and no limit

Write u = 1 + |G| and r = σu³/p². Let D be the set of states with 0 < p ≤ 1/4 and r ≤ 1/5.

**Lemma 1.** In D the cubic has exactly one root q in [−3p²/2, −p²/2], i.e. w ∈ [−p²/2, p²/2], and |w| ≤ 20σu²p².

*Proof.* Put q = −p² + w, |w| ≤ p²/2, so |q| ≤ 3p²/2. With p ≤ 1/4, u ≥ 1, term by term:

- |A₁| ≤ p⁶(2.25u² + 6u + 6.75u + 4 + 6.75) ≤ 26u²p⁶;
- |∂A₁/∂q| = |2G²p²q + 4Gp⁴ − 6Gq² − 6p²q| ≤ p⁴(3u² + 4u + 13.5u + 9) ≤ 30u²p⁴;
- X = G(q − p) + p² has |X| ≤ 1.625up, Y = G²q² + 2Gp³ + 2Gp²q + 2p⁴ − q² has |Y| ≤ 4.375u²p³, so |A₂| = 2|XY| ≤
  15u³p⁴;
- |∂A₂/∂q| = |2GY + 2X(2G²q + 2Gp² − 2q)| ≤ 8.75u³p³ + 26u³p³ ≤ 35u³p³.

In D, σu² ≤ σu³ ≤ p²/5 and σ ≤ p²/5. The remainder R = σA₁ + σ²A₂ obeys |R| ≤ 26σu²p⁶ + 15σ²u³p⁴ ≤ p⁶(5.2p² + 3σ) <
0.37p⁶, and
|∂R/∂w| ≤ 30σu²p⁴ + 35σ²u³p³ ≤ 6p⁶ + 7σp⁵ < 0.4p⁴. Now N = A₀(w) + R with A₀(w) = p²w(2p² − w): A₀(−p²/2) = −5p⁶/4, A₀(p²/2) = 3p⁶/4 and
∂A₀/∂w = 2p⁴ − 2p²w ≥ p⁴ on the interval. So N < 0 at the left end, N > 0 at the right end and ∂N/∂w > 0.6p⁴: one
root. At it p²|w|(2p² − w) = |R| with 2p² − w ≥ 3p²/2, so |w| ≤ |R|/(1.5p⁴) ≤ 17.4σu²p² + 10σ²u³ ≤ 20σu²p². ∎

**Theorem 1 (proved, by hand; the four polynomial bounds are also sampled by `trap_bounds`).** Let G′ = p, p′ = q with
N(G, p, q, σ) = 0 on [t₁, ∞) and (G, p) ∈ D with q ∈ [−3p²/2, −p²/2] at t₁. Then for all t ≥ t₁ the state stays in D
and q stays in that interval, and

1. −5p²/4 ≤ p′ ≤ −3p²/4, so 1/p₁ + 3(t − t₁)/4 ≤ 1/p(t) ≤ 1/p₁ + 5(t − t₁)/4;
2. G(t) − G(t₁) ≥ (4/5) ln(1 + 5p₁(t − t₁)/4) → ∞: **no limit**, F → +∞ at least like (4/5) ln ln η;
3. T = 1/p − t and g = G + ln p have limits T_∞, g_∞, with |T′| ≤ 20σu² and |g′| ≤ 20σu²p, and p = 1/(t + T(t)).

*Proof.* In D, |w| ≤ 20σu²p² ≤ 4p⁴ ≤ p²/4, so p′ = −p² + w ∈ [−p²(1 + 20σu²), −p²(1 − 20σu²)] and 20σu² ≤ 4p² ≤ 1/4
(σu² ≤ p²/5); this is 1. (ln r)′ = −2 + 3u′/u − 2p′/p with u′ ≤ p and −p′/p ≤ 5p/4, so (ln r)′ ≤ −2 + 5.5p ≤ −0.625 <
0 for p ≤ 1/4; p′ < 0 keeps p ≤ 1/4, and 1 keeps p > 0. So D is forward invariant. N ≠ 0 at both ends of the interval
for every state of D (Lemma 1), so a continuous q(t) in the interval stays in it; ∂N/∂q ≠ 0 there, so q is a smooth
function of the state and the solution exists for all t ≥ t₁ (p is bounded below by item 1, and G′ = p ≤ 1/4). 2 is
the integral of the lower bound in 1. For 3, T′ = −w/p² and g′ = w/p with |w| ≤ 20σu²p², u ≤ u₁ + (t − t₁)/4 and σ = e⁻²ᵗ, so both are integrable. ∎

A solution with a limit F_∞ ≠ 0 has p = 4F_∞σ + … (`conical_limit.md` §3), hence r ≈ u³/(16F_∞²σ) → ∞ and |q| ≈
8|F_∞|σ ≫ p²: it is never in D, while a solution of Theorem 1 never leaves D. So a solution in D cannot cross over to
#159's family at any η, however large; that is a property of the equation, not of the precision of a walk.

## 4. The walks

**The walker.** `conical_global/src/hp.rs` is #167's implicit Taylor recursion (`poly.rs`: order k of N = 0 solved for
F″ₖ, one division per step) written once over multidouble's `Precision`, so that it runs in double-, quad- and
octo-double lanes. Every constant is formed in the lane: 1/k, the coefficients of N about the step's x₀ (from exact
rationals and the f64 x₀), the division by ∂N/∂F″; no double-double constant enters a quad-double lane (control
`recip-dd`). Three things were needed and are not in #167's walker. (i) The step variable is s = (x − x₀)/L with L a
power of two near the previous step: the coefficients about x₀ scale like x₀⁻ᵏ and underflow the f64 range at x₀ ≳
3·10⁵ (order 60) and 10⁹ (order 30), which gave NaN (`dev_00` in the raw folder); L changes only the constants (1/k →
L/k, eⱼ → eⱼLʲ) and is exact. (ii) The planner takes the last two terms relative to each field's own magnitude and has
no largest step: #167's absolute rule with a cap of 10 reaches η = 3.9·10⁴ in 4000 steps, and an absolute rule lets
the step pass the radius of convergence once |F′| falls below the tolerance. (iii) A mode `project` (one Newton step
on F″ at each step's start, so N = 0 holds to the lane's rounding of the *current* scale) beside `conserve` (#167's:
order 0 of N is not imposed and the rounding of every step stays in N). Planted truth through the same code (`P1` to
`P5`, below): #167's linear, cubic and blow-up cases unchanged; q + p² = 0 to 10¹² with errors 8.4e-32 (DD), 1.1e-64
(QD), 1.2e-126 (OD, `project`); a planted σ-shift to 9.0e-31; a planted linear equation with a limit, F to 1.7e-32 and
verdict LIMIT.

**The forward grid.** #167's 100 solutions (F(0) ∈ {0.1, 0.25, 0.5, 0.75, 1.5, 2, 3, 5, −0.5, −2}, F′(0) ∈ ±{0.01,
0.1, 1, 10, 100}, F″(0) = −F′²/F), plus sides, whose minus side is the plus side of (f, −f₁). A rerun of #167's crate
here gives its class and kind labels at 100 of 100 points and η* within 1e-8 of #167's raw records at 82 of 85 points
that have one. The LOGLOG set is #167's 15. The three differences are κ-locus points (1.5, 0.01), (2, 0.01), (−2,
−0.01), where η* is the x of an F″ overflow that #167 states to 1e-6 only: they differ from the Mac's raw records by
4.7e-7 to 7.9e-7 (W1 as first worded is printed `NOT MET`). Walks of the 15 to 10¹²:

- double-double `project`: 0.2 s per member and quad-double `project` 6 to 7 s (raw `dev_04`); they agree to 2.0e-29
  in G, p, q at 200, 10³, 10⁴, 10⁶, 10⁹, 10¹² (`dd_qd`); N is zero to 2.5e-32 (DD) and 3.3e-64 (QD) of the sum of its
  terms at the step ends (`n_floor`).
- ρ_K = |q − Q_K|/p² for the worst member (5, 100): ρ₀ 4.5e-2, 1.7e-3, 1.7e-5, 1.8e-9 and 1.8e-21 at η = 200, 10³,
  10⁴, 10⁶, 10¹²; ρ₄ 7.7e-8, 1.0e-14, 3.0e-24 at 200, 10³, 10⁴, then roundoff (3e-30, 6e-31); ρ₆ 1.1e-10 at 200. The
  15 have ρ₄ ≤ 1.0e-14 (10³), 5.8e-30 (10⁶), 2.1e-29 (10¹²) against the bar 1e-10. The whole table is in the output.
- the verdict at 10¹² is LOGLOG for all 15 (none LIMIT): G(10¹²) − G(10⁶) is 0.675 to 0.720, equal to ln((t₂ + T₁)/(t₁
  + T₁)) with T₁ = 1/p − t at 10⁶ to 2.4e-11.
- an independent integration (30 digits, `mp.odefun`, no Rust) of the reduced system G′ = p, p′ = Q₆(G, p, e⁻²ᵗ) from
  η = 10³ (the quad-double state) matches the walk at 10⁶, 10⁹, 10¹² to 4.0e-22 at worst, for (2, 1), (5, 100) and
  (−2, −10). It shares N and the start state, not the Taylor recursion or the root solve.
- #167's walker and the new double-double `conserve` walker differ at η = 200 by up to 5.2e-18 (W0, bar 1e-24: **NOT
  MET as first worded**). The cause is the conserved N: both walks carry the rounding of N from their first steps
  (n_rel up to 3.1e-17 over the 15 walks) while the scale of N falls like p⁶, so both are only good to about their own
  n_rel; the same member walked by `project` equals the quad-double one to 7e-32 (the dev run in raw `dev_08`), and a
  quad-double walk with F′(0) changed by 1e-30 relative moves the state at 200 by 1e-31, so the map is well
  conditioned. W0 is printed, does not gate, and `dd_qd` is the accuracy gate. The same floor sits in #167's forward
  walks for F′(0) = 100; it does not change their classes.
- the parity of N: the 100 singular sides agree with their mirror pairs in class, kind and η* (0 difference), and the
  walk of (f, −f₁) to −200 equals that of (f, f₁) to +200 for the 15 members (difference 0; the walk is bitwise
  mirror-symmetric). W9 as worded (compare the two walks of a *singular* pair at ±200) cannot be done and is printed
  `NOT MET`.

**The finer map.** 20 values of F(0) in [0.1, 5] (log-spaced) and ±20 values of F′(0) in [0.01, 100], DD `conserve`
walks to 200: 109 of 800 reach 200 in the band q = ∓p², all at F(0) > 1 and F′(0) > 0 (the mirror at F(0) < −1, F′(0)
< 0 is by the symmetry); none at F(0) ≤ 0.963. The smallest F′(0) that reaches 200, by bisection with both ends
checked, is 0.1971 at F(0) = 1.05, 0.2513 at 1.25, 0.2372 at 1.5, 0.1995 at 2, 0.1492 at 3, 0.1004 at 5 and 0.05607 at
10 (at F(0) = 30 the interval [0.01, 100] does not bracket it).

## 5. The precision a walk needs to see a crossover

A crossover to #159's family means p falls from about 1/ln η to about 4F_∞/η². In `conserve` mode N is a constant of
the flow: the data and the first steps leave N = ε ≈ u S₀ (u the unit roundoff, S₀ the sum of the terms of N at the
start), and N has no term below degree 6 in (p, q), so on the family's side the scale is S ≈ (p_c/p₀)⁶ S₀ and the walk
follows N = 0 there only if u ≪ (4σ_c/p₀)⁶, p_c ≈ 4F_∞σ_c. With p₀ = 0.05 and F_∞ = 1 (`precision_table`, printed):

| crossover at η_c | needed u ≪ | double-double (1e-32) | quad-double (1e-64) | octo-double (1e-127) |
|---|---|---|---|---|
| 10³ | 2.6e-25 | yes | yes | yes |
| 10⁵ | 2.6e-49 | no | yes | yes |
| 10⁶ | 2.6e-61 | no | yes | yes |
| 10⁹ | 2.6e-97 | no | no | yes |
| 10¹² | 2.6e-133 | no | no | no |

("yes" means u ≤ the entry/100; with that margin double-double reaches η_c = 2.8·10³, quad-double 1.3·10⁶, octo-double
2.3·10¹¹.) So a `conserve` walk cannot see a crossover at 10¹² in any lane precision, and the `project` mode, which
restores N = 0 to the rounding of the current scale at every step, needs only the relative accuracy of (G, p). The
mechanism is tested on a planted equation, N = p⁴(q + 3p), where p falls by 2.7·10⁷ and the scale of N by 7·10⁻³⁸ from
η = 10 to 3000: the relative error of F′ at 3000 is 11 for DD `conserve` (the prediction u(p₀/p)⁵ is 1e5: the error
saturates), 1.1e-28 for QD `conserve` (prediction 1e-27), 1.0e-31 for DD `project` and 2.3e-64 for QD `project`: the
four outcomes are as predicted (`plant_scale`). On the real N the same floor shows in the log-log regime: at η ≥ 10⁹
the DD `conserve` walk of (5, 100) returns (q + p²)/p² with the wrong sign (+9.8e-15 against the true −1.8e-15, raw
`dev_04`), the `project` walk −1.796e-21 at 10¹², which is the series' −σ(G − 1)²/2.

What this means for the question. Inside D no precision is needed, because Theorem 1 excludes a crossover at every η.
Before η_in (at most 1.3·10⁴) a `conserve` walk would need u ≪ 1.1e-38 (the table's formula at η_c = 1.3·10⁴), which
quad-double meets and double-double does not. The production walks are `project`, whose floor is the rounding of (G,
p), and the double- and quad-double `project` walks agree to 2.0e-29, so the floor hides no crossover. The octo-double
lanes run (planted case to 10¹², 1.2e-126) and `Relaxed` is not used: a probe (raw `relaxed_probe`, multidouble rev
ba0ad101ce76) forms 0 squares at DD order 30 and QD order 60 (bitwise `cs_range`, no faster: 0.387 against 0.345 ms at
QD order 60) and 14 squares at OD order 120 (1.5 times faster), 39 at order 240 (2.2 times).

## 6. Not decided

- Solutions that never enter D: the 85 other grid solutions end at singularities (#167's table), and the finer map is
  a map to η = 200 in double-double `conserve`, not a classification. Whether a solution outside the forward map (data
  off F″(0) = −F′²/F, or an inward family member) can be on the small-root branch inside D is not covered.
- Entry of the 15 into D is numerical (the walk's state at η_in); Theorem 1 is for every later η. Its four polynomial
  bounds are hand estimates sampled at 4000 random states (largest ratio to the bound 0.30), not enclosed, and the
  proof has had no independent review. Nothing is certified.
- The root series is checked to k = 6; its radius of convergence in σ and the behavior of eₖ for large k are not
  studied.
- p < 0 is by the symmetry N(−G, −p, −q) = N; mirrored members were walked directly and equal their mirrors.
- The limit solutions of #159 lie outside D (r → ∞) on a different sheet of the cubic; whether a solution that is in
  neither D nor the family exists at large η (for example with p ≲ σ and q on the third root) is not decided.
- The statements about `conserve` walks at crossovers are a by-hand scaling, tested on a planted equation only.
- W0 (§4), W1 at three κ-locus points and W9 as worded are not met.

## 7. What `conical_global.md` §4 and §5 should now say

§4, after "This is the ν⁰ part A₀ = …": "The regime is expanded in `conical_loglog.md` (#184): q = −p² + Σσᵏeₖ(G, p)
with e₁ = −p²(G − 1)²/2 (**checked symbolically** to k = 6); T = 1/p − t and g = G + ln p tend to limits T_∞, g_∞ (T′
= −w/p², g′ = w/p, w = q + p²), so p = 1/(ln η + T_∞) + O(σ(1 + |G|)²) and F = ln(ln η + T_∞) + g_∞ + …. All 15
plus-side LOGLOG solutions enter the region D = {0 < p ≤ 1/4, σ(1 + |G|)³ ≤ p²/5} by η = 1.3·10⁴ and follow the series
to 1e-10 from η = 10³ (**numerical evidence**; walked to 10¹² in double- and quad-double). Theorem 1 of that note
(**proved**, by hand) shows that a solution on the root branch in D has p′ ∈ p²[−5/4, −3/4] and F → +∞ at least like
(4/5) ln ln η; it has no limit and cannot cross over to #159's family at any η."

§5, the fourth bullet, "the log-log regime is described, not proved, and whether a solution of it can be continued to
a limit at huge η … was not followed beyond η = 10⁴": replace by "the log-log regime is expanded and, for solutions in
D, shown to have no limit (`conical_loglog.md` §3, **proved**); the 15 forward-grid members enter D by η = 1.3·10⁴.
Not covered: solutions that are not in the forward grid and never enter D."

§2's last sentence ("a solution that crosses with F′(0) ≠ 0 … no solution with a limit on either side") stands; the §8
text for `conical_limit.md` §6 "grow like ln ln η" may add "(no limit: `conical_loglog.md` Theorem 1)".

## 8. Lemmas for Lean (#135)

Statements only, no Lean written. Over ℚ or ℝ, in the variables of §1.

- (L1) ring: N(G, pν, (q − p)ν², 1/ν) = A₀(p, q) + ν²A₁(G, p, q) + ν⁴A₂(G, p, q).
- (L2) ring: A₀(p, −p² + w) = p²w(2p² − w); A₁(G, p, −p²) = p⁶(G − 1)²; ∂A₀/∂q(p, −p²) = 2p⁴.
- (L3) ring in ℚ[G, p, p⁻¹]: N(G, p, −p² + σe₁ + σ²e₂ + σ³e₃, σ) has no term of order σ⁰ … σ³ (e₁, e₂, e₃ as printed;
  the script checks to σ⁶).
- (L4) field_simp: for T = 1/p − t, g = G + ln p, G′ = p, p′ = q: T′ = −(q + p²)/p², g′ = (q + p²)/p.
- (L5) real inequalities (nlinarith or polyrith): the four bounds of Lemma 1 for 0 < p ≤ 1/4, |w| ≤ p²/2, u ≥ |G| + 1.
- (L6) Lemma 1 from (L5): intermediate value theorem, strict monotonicity, |w| ≤ 20σu²p².
- (L7) Theorem 1 from (L6): an ODE comparison for r and 1/p, no analysis beyond continuity of the solution.
- (H) hypothesis of Theorem 1: a C¹ solution (G, p) with p′ = q a root of N in the stated interval at t₁.

## 9. Checks and controls

20 gating checks pass in the run of record and 17 of 17 controls fire on their target only; the three `NOT MET (as
first worded)` lines (W0 `regress`, W1 at three κ-locus points, W9 as worded) are results and do not gate the exit
code. The names are those of the raw README (criteria and bars written before the first run; amendments A1 and A2
after the development runs, in particular the coordinator's rulings on W0, W1 and W9).

| check | what it tests | result |
|---|---|---|
| `theta_form` | `conical_N` in (G, p, q, ν) equals A₀ + ν²A₁ + ν⁴A₂ exactly; the crate's N at 40 rational points | 0; 2.1e-32 (bar 1e-26) |
| `root_series` | e₁ … e₆ solved in sympy, residual by independent Laurent arithmetic; A₁(G, p, −p²) = p⁶(G − 1)² | no term below σ⁷; e₁ as printed |
| `shift_series` | the same generator on a planted shift | e₁ = −2c(1 − p), e₂ = −c², e₃ = e₄ = 0 exactly |
| `invariants` | T′, g′, and the exact solution, in sympy | all True |
| `basin` | #167's forward grid rerun against its table (labels) and raw records (η*); the 15 reach 10¹² as LOGLOG | 100 of 100; 82 of 85 within 1e-8; 15 of 15 |
| `approach` | ρ₄ ≤ 1e-10 at η = 10³, 10⁶, 10¹² for the 15 | 1.0e-14, 5.8e-30, 2.1e-29 |
| `no_limit` | verdict at 10¹²; G growth against ln((t₂ + T₁)/(t₁ + T₁)) | 15 LOGLOG; 2.4e-11 |
| `trap_bounds` | the four polynomial bounds of Lemma 1 at 4000 random states; the root, \|w\| ≤ 20σu²p², p′ ∈ p²[−5/4, −3/4], (ln r)′ < −3/5 at 300 states of D | worst ratios 0.12, 0.30, 0.09, 0.25; p′/p² ∈ [−1.001, −1.0]; (ln r)′ ≤ −1.003 |
| `trap_entry` | each of the 15 enters D at η_in ≤ 10⁶ and stays with r non-increasing | η_in ≤ 1.27·10⁴ |
| `invariants_conv` | \|T(η) − T(10¹²)\|, \|g(η) − g(10¹²)\| at 10⁶, 10⁹ (QD) | 8.8e-10, 5.4e-11 (bar 1e-8) |
| `n_floor` | N/Σ\|terms\| at the step ends, `project` walks | 2.5e-32 (DD), 3.3e-64 (QD); `conserve` 3.1e-17 printed |
| `dd_qd` | DD `project` against QD `project` at six checkpoints | 2.0e-29 (bar 1e-26) |
| `independent` | mpmath integration of the reduced system from 10³ against the QD walk | 4.0e-22 (bar 1e-12) |
| `parity` (W9′) | mirror pairs: 100 singular sides; 15 regular sides at ±200 | 100 of 100, η* 0; difference 0 |
| `plant_167` | #167's planted matching and blow-up cases | 6.5e-28 (C⁻), 3.5e-29 (C⁺), 7.5e-12 (η*) |
| `plant_loglog` | q + p² = 0 to 10¹², DD, QD, OD, both modes | 8.4e-32, 2.8e-30; 1.1e-64, 2.5e-62; 1.2e-126, 7.8e-123; LOGLOG |
| `plant_shift` | a planted σ-shift to 10¹² (DD `project`) | 9.0e-31; (p − 1/(t + T))/(cσ) = 1.0 at 10⁶ |
| `plant_limit` | a planted equation with a limit to 10¹² | F to 1.7e-32; LIMIT |
| `plant_scale` | N = p⁴(q + 3p): DD/QD, `conserve`/`project` | 11, 1.1e-28, 1.0e-31, 2.3e-64: as predicted |
| `recip_check` | QD step-halving with lane-precision constants | 7.6e-65 (bar 1e-55); the control 2.0e-33 |

Printed and not gating: W0 `regress` (5.2e-18 against the bar 1e-24: NOT MET), the W1 κ-locus points (4.7e-7 to
7.9e-7: NOT MET), W9 as worded (NOT MET), the 20 × ±20 basin map and its boundary, the precision table, the η_in
table.

| control `--mutate=` | change | target check |
|---|---|---|
| `theta-coeff` | 4Gp⁴q → 5Gp⁴q in A₁ | `theta_form`, `root_series` |
| `root-coef` | one monomial of e₂ perturbed in the residual | `root_series` |
| `shift-series` | e₂ expected as −c²/2 | `shift_series` |
| `inv-sign` | T′ expected as +w/p² | `invariants` |
| `class-sign` | the LOGLOG test reads q = +p² at p > 0 | `basin`, `no_limit` |
| `series-k` | K = 1 in place of 4 | `approach` |
| `walk-eq` | one coefficient of N changed in the walker | `approach`, `independent` |
| `indep-q` | the independent integration uses Q₀ = −p² | `independent` |
| `inv-conv` | T = 1/p + t | `invariants_conv` |
| `parity-sign` | F′ summed with the wrong sign | `parity` |
| `trap-bound` | 2.6u²p⁶ in place of 26u²p⁶ | `trap_bounds` |
| `trap-eps` | D with ε = 1e-12 | `trap_entry` |
| `plant-exact` | the comparison solution has g + 1e-3 (the start state does not) | `plant_loglog` |
| `plant-c` | c = 1003/3000 in the planted equation | `plant_shift` |
| `limit-class` | the LIMIT test reads \|p\| ≥ 1e-6 | `plant_limit` |
| `no-project` | `project` runs as `conserve` | `plant_scale`, `n_floor`, `dd_qd` |
| `recip-dd` | QD lanes with double-double constants | `recip_check` |

Two controls were wrong in the development run `dev_09` and were reported, not adjusted silently: `plant-exact` was
silent (the mutation also shifted the start state, and F does not enter q + p² = 0) and is now a comparison-only
shift; `limit-class` was predicted to fire P2 as well and does not (G changes by 0.6, so the |p| test alone does not
make it LIMIT). `walk-eq` raised a `KeyError` in `independent` when its mutated walk ended early; both checks now
count an early end as a failure. A control that makes the DD walks `conserve` (`no-project`) fails `dd_qd` at 4.1e-14.

Run from the repository root, in the environment of the repository's setup (Python 3.14.10, sympy 1.14.0, mpmath
1.3.0, numpy; Rust 1.100.0, `multishoot =4.2.0`, `multisolver =0.2.2`, `multidouble =5.1.0`, built by the script):

    python3 notes/warp/conical_loglog.py                     # the 20 checks and the 17 controls: 1534 s on this machine (680 s of checks), peak 100 MB, at nice 10 with 2 threads
    python3 notes/warp/conical_loglog.py --no-controls       # the checks only
    python3 notes/warp/conical_loglog.py --mutate=NAME --only=CHECK

## 10. Record

Raw folder `I184-rpi5` of negative-mass-raw: the README with the criteria, bars, planted cases with their predicted
outcomes and control names written before the first run (commit 1792163), amendments A1 (the development runs; the
coordinator's rulings on W0, W1, W9; the checks added after) and A2 (the control targets and what each did in
development), the development logs `dev_00` … `dev_12` (including the first unlogged runs, the failed and silent
ones), `run_01` (a failed start: no `cargo` on the systemd unit's PATH), `run_02` (the first complete run, kept) and
`run_03`, the run of record, and `relaxed_probe/`. No library change was needed, so there is no pull request;
multidouble's series module was probed, not used (§5).

The crate additions are `src/hp.rs` and `src/bin/loglog.rs`; `src/main.rs`, `poly.rs`, `walk.rs` and the tables are
unchanged, and #167's planted checks and `forward` give the same results.

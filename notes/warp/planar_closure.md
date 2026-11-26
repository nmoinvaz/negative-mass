# Closure for planar shifts that cross |u| = 1 or do not travel

November 2026, issue [#161](https://github.com/nmoinvaz/negative-mass/issues/161). Script
[`planar_closure.py`](planar_closure.py), its stdout [`planar_closure_output.txt`](planar_closure_output.txt). It
builds on [`det_classify.md`](det_classify.md) §2 (#80) and [`two_scalars.md`](two_scalars.md) (#65). It imports
`ricci` from [`phantom_source.py`](phantom_source.py) and `closure`, `dform` from [`two_scalars.py`](two_scalars.py),
and retypes #80's formulas (it does not import `det_classify.py`, which reads `--mutate` at import).

Status words as in AGENTS.md: **proved** (by hand), **checked symbolically**, **numerical evidence**, **not decided**.

## Answer

**Q2, traveling shifts (static u(ξ)): no localized one passes closure, whether or not it crosses |u| = 1.** On every
component of {|u| ≠ 1} the map s′ = 2/(1 − u²) is a bijection onto ℝ, while Φ = K(1 − u²) tends to 0 at both ends, and
closure forces Φ = C e^{αs}, which cannot do that unless C = 0 (Theorem T′, **proved**; steps **checked symbolically**
on u = 3/(1 + ξ²), a superluminal bubble that crosses twice). So a superluminal planar bubble that travels at constant
speed is excluded, under u ∈ C² (K bounded at the horizons), C³ where K ≠ 0, and limits with |u_±| ≠ 1.

**Q1, shifts that do not travel: a localized closure-passing b(t, x) exists, with unequal limits at ±∞.** For σ =
H(p)J(q) with H′, J′ ≥ 0 supported in [0, 1], the normal geodesics of the line p + q = −2 give a global PG chart with
b = 0 for t ≤ 0, b constant in x outside a compact set at every t, b_− = 0 and b_+(t) = −L′(t) nonincreasing to
−2ΔHΔJ, and 0 ≤ −b ≤ 2ΔHΔJ (Theorem E, **proved**; the example h = j = 3/2 is **numerical evidence**). It has |b| > 1
at late times, so it is a closure-passing planar shift that does cross |b| = 1. It is not a bubble: b_+ ≠ b_−. For
**equal limits** at every t the question is **not decided**; what is proved is that a sign-definite K (H′J′ ≥ 0) is
impossible, that ΔHΔJ = 0 is needed if the chart covers the null plane, and that at linear order no compactly
supported data work (Paley–Wiener).

Two criteria fail as written and stay failed by the coordinator's rulings (§7): `crit_pg_iv` (C2 (iv): K = 0 at one
sample point) and `E_hol` (E4: a quadrature slack of 1e-9 against an excess of 4.8e-6 that was s-quadrature error).
The script exits 1 because of them.

## 1. The criterion

Planar b(t, x): R = K g₂ ⊕ 0 with K = −a_x, a = b_t − b b_x (#80 §2, `crit_pg` (i)). #80 proved that the normal +
phantom closure holds iff σ_pq = F(p)G(q) for g₂ = −2e^{2σ} dp dq, i.e. g₂ = −2 exp(2H(p)J(q)) dp dq. In covariant
form, with □ the 2D d'Alembertian and C[b] = □ ln|K| − 2K:

    d_p d_q ln|R_pq| = −(e^{2σ}/2) C[b],   □f = −2 e^{−2σ} f_pq,   K = 2σ_pq e^{−2σ},   R_pq = −2σ_pq,

so closure holds on {K ≠ 0} iff C[b] = 0 (`crit_null` (i), **checked symbolically**, for an undefined σ). The six
planted σ of #80 agree with the code `closure`: C = 0 for the three passes, C ≠ 0 for the three failures (`crit_null`
(ii)).

For a static shift u(ξ) = b(x), C[u] = Θ′ exactly, Θ = (1 − u²)K′/K − 2uu′, and the obstruction of `closure` is
dω_{tξ} = Θ′/2 (`crit_pg` (ii), (iii), undefined u; **checked symbolically**). In a PG chart Θ′ = ((1 − u²)K′/K)′ − 2K
has no division by 1 − u², so the criterion is regular across a horizon wherever K ≠ 0. The null coordinates are not:
p or q runs to infinity there.

A formal planted pass, with the same code: for the jets u = 1/3, u′ = 1/2, u″ = 1/5, u‴ = 2/7 at ξ = 0 the Taylor
coefficients c₄ … c₈ can be solved so that C[u] = O(ξ⁵); with c₆ + 1/7 the ξ² coefficient of C[u] is nonzero
(`crit_pg` (v), exact rationals).

## 2. Q2: traveling planar shifts, with or without crossings

**Theorem T′ (proved).** Let u ∈ C²(ℝ), u ∈ C³ on {K ≠ 0}, K = (uu′)′, with limits u_± as ξ → ±∞, |u_±| ≠ 1, and u′,
u″ tending to 0. If closure holds on {K ≠ 0}, then K ≡ 0, so u is constant.

*Proof.* Let I be a component of {|u| ≠ 1}, with ends e₋ < e₊ (finite, where |u| = 1, or infinite). On I, s′ = 2/(1 −
u²) has one sign, so s is strictly monotone. At a finite end 1 − u² = O(|ξ − e|) since u ∈ C¹, so s′ is not integrable
and |s| → ∞. At an infinite end s′ → 2/(1 − u_±²) ≠ 0, so |s| → ∞. A monotone function with |s| → ∞ at both ends tends
to +∞ at one and −∞ at the other: s maps I onto ℝ. Next Φ = K(1 − u²) → 0 at both ends: at a finite end K is bounded
and 1 − u² → 0, at an infinite end K = u′² + uu″ → 0. Let J be a maximal interval of {K ≠ 0} inside I. On J closure is
Θ′ = 0, Θ = 2α, which is Φ′/Φ = 2α/(1 − u²) = α s′, i.e. Φ = C e^{αs}, C ≠ 0. If an end of J were an interior point
ξ_b of I, then K(ξ_b) = 0 by maximality while Φ → C e^{αs(ξ_b)} ≠ 0 (s is finite there), a contradiction with
continuity. So J = I. Then Φ = C e^{αs} tends to 0 where s → +∞ only if α < 0 and where s → −∞ only if α > 0, and α =
0 gives Φ = C ≠ 0. All three are impossible, so K ≡ 0 on I. Hence K = 0 on {|u| ≠ 1}, and by continuity on its
closure; on the interior of {|u| = 1} u = ±1 is constant and K = 0. So (uu′)′ ≡ 0, uu′ = c₁, and the limits force c₁ =
0 and u constant. ∎

This extends Theorem T of #80, which assumed |u| ≠ 1 on all of ℝ (one component, no horizons). Nothing about the
behaviour at the horizons is used except that K is bounded there; the closure condition across a horizon is not
needed.

**What Θ′ = 0 requires at a crossing** (`horizon` (i), (ii), exact). Let |u(ξ₀)| = 1 and K(ξ₀) ≠ 0. Θ is continuous
through ξ₀, and (1 − u²)K′/K → 0 there, so Θ(ξ₀) = −2uu′(ξ₀): the constant is **α = −u u′(ξ₀)**, and it is the same on
both sides. Near a simple crossing u = 1 + u₁η + …, a solution of (1 − u²)K′ = (2α + 2uu′)K has K ~ |η|^γ with **γ =
−(α + u₁)/u₁** (exact leading order), bounded iff α/u₁ ≤ −1. For the superluminal bubble u_s = 3/(1 + ξ²) (u_s = 1 at
ξ = ±√2, u_s > 1 between, K(±√2) ≠ 0) the two horizons give α = 2√2/3 and α = −2√2/3. One connected component of {K ≠
0} has one α, so it cannot contain both horizons: closure through both is impossible already for this reason. T′ shows
more: no closure at all, with or without a horizon in the component.

**Checks on explicit profiles** (`horizon` (iii), exact): on each of the three components of {|u_s| ≠ 1} ((−∞, −√2),
(−√2, √2), (√2, ∞)) and on the one component of u₂ = 1/(2(1 + ξ²)) the sign of s′ is constant, the residue or limit of
s′ at each end gives divergence of s of opposite signs at the two ends, Φ → 0 at both ends (exact limits), and K ≠ 0
at an interior point. So Φ = C e^{αs} with C ≠ 0 is excluded on all four, which is the proof's mechanism.

**Regularity.** Assumed: u ∈ C² (so K is bounded at a horizon), C³ where K ≠ 0, limits with u′, u″ → 0. Not covered: a
shift that is only C¹ at a horizon with K unbounded there. *By hand, not checked:* for α = 0 the closure ODE is (u²)″
= 2C/(1 − u²), which integrates to ((u²)′)² = −4C ln|1 − u²| + 2E, so a solution that reaches a horizon has u′ → ∞
there and is not C¹. For other α the case is **not decided**.

**Planted truth.** Theorem T's case, through the same code: the subluminal traveling bubbles u₁ = ½(1 − 1/(1 + ξ²)²)
(#80's planted failure), u₂ = 1/(2(1 + ξ²)) and u₃ = 1/(3(1 + ξ²)²) have K ≠ 0 and C[u] ≠ 0. As worded in the README
(5 fixed points ξ = 1/3, 1, 3/2, 2, 5/2) this holds for u₁ and u₂ and fails for u₃ at ξ = 1/3, where K = 4(9ξ² −
1)/(9(1 + ξ²)⁶) = 0: **`crit_pg_iv` FAILS** (§7). The check that came after it, `crit_pg_ivb`, takes the first five
points of the fixed list 1/5, 1/3, 1/2, 2/3, 1, 3/2, 2, 5/2, 3, 4 at which K ≠ 0 is verified exactly and finds C ≠ 0
at all 15 points (u₁, u₂: 1/5, 1/3, 1/2, 2/3, 1; u₃: 1/5, 1/2, 2/3, 1, 3/2).

## 3. Q1: shifts that do not travel

**Localization used.** (L1) limits: b(t, x) → b_±(t) as x → ±∞ at every t, here with b_x of compact support in x at
each t, so the limits are attained; (L2) equal limits b_+ = b_−, i.e. b → 0 at both ends after x → x + c(t), which is
the Alcubierre-like bubble; (L3) bounded b. Time interval: all t ∈ ℝ.

**PG charts of a 2D metric** (**proved** by hand, steps **checked symbolically** in `jacobi`). In a chart −dt² + (dx +
b dt)², g^{tt} = −1, so n = −∇t = ∂_t − b∂_x is a unit geodesic field. Hence a PG chart is the Gaussian normal chart
of any slice S₀ = {t = t₀}: t is proper time along the normal geodesics γ_s, x = X(t, s) is arclength along the slices
plus a gauge c(t), and with A = X_s, b = −X_t,

    g₂ = −dt² + A² ds²,  R = (A_tt/A) g₂,  b_x = −A_t/A,  n(b_x) − b_x² = −A_tt/A,  a = −X_tt,  K = −a_x = A_tt/A.

So A_tt = K A (Jacobi), b_x = −θ is minus the expansion of the geodesics, and with L(t) = ∫A ds the length of a slice
between two comoving observers, L′ = ∫A_t ds = b_− − b_+ and L″ = ∫K A ds = ∫K dx = a_− − a_+ (`jacobi` (iii) checks
the last identity for b = γ + β(1 + tanh x)/2, where it equals −β′). Conversely any spacelike curve S₀ gives a chart.
So the question is which S₀ give b → const at both ends at every t.

**Theorem E (proved).** Let H, J ∈ C^∞(ℝ) with H = J = 0 on (−∞, 0], H′ ≥ 0, J′ ≥ 0, H′ = J′ = 0 on [1, ∞), and H₁ =
H(1), J₁ = J(1) positive. Let σ = HJ and g₂ = −2e^{2σ} dp dq on ℝ². Then:

- (a) closure holds (σ_pq = H′(p)J′(q) is a product); K = 2H′J′e^{−2HJ} ≥ 0, K ≢ 0, supp K ⊂ [0, 1]²;
- (b) the unit normal geodesics of S₀ = {p + q = −2} define a PG chart on all of ℝ² (t, x ∈ ℝ);
- (c) b = 0 for t ≤ 0, and for every t there is X(t) with b(t, x) = 0 for x ≤ −X(t) and b(t, x) = −L′(t) for x ≥ X(t)
  in the gauge b_− ≡ 0; L′ is nondecreasing, L′ ≥ 0;
- (d) 0 ≤ −b ≤ 2H₁J₁ everywhere, and L′(t) = 2H₁J₁ once the slice lies above the diamond.

*Proof.* The region W = {p ≤ 0} ∪ {q ≤ 0} has σ = 0: it is Minkowski, g = −2 dp dq, and S₀ ⊂ W is a straight spacelike
line. Geodesics. For g₂ the nonzero Christoffel symbols are Γ^p_pp = 2σ_p, Γ^q_qq = 2σ_q (`exterior` (i)), so p̈ =
−2σ_p ṗ² and (1/ṗ)˙ = 2σ_p = 2H′(p)J(q) ≥ 0: ṗ stays positive and finite for all τ (ṗ never increases); the same for
q. The unit normal geodesic through (p₀, q₀) = (−1 − s/√2, −1 + s/√2) has ṗ = q̇ = 1/√2 initially and satisfies
2e^{2σ}ṗq̇ = 1. For τ < 0 it stays in {p + q < −2} ⊂ W, a straight line. So E(τ, s) is defined on ℝ². Jacobi field.
The separation A = e^{σ}√(−2p_s q_s) equals p_s/ṗ, and A_ττ = K A follows from the geodesic equations, the
normalization and the orthogonality p_s q̇ + ṗ q_s = 0 alone (`exterior` (ii), **checked symbolically** on jets). A(0)
= 1, A_τ(0) = 0, K ≥ 0 give A_τ ≥ 0 and A ≥ 1 for τ ≥ 0, and A = 1 for τ ≤ 0 (K = 0 on W). The Jacobian of E is p_τ
q_s − p_s q_τ = −A e^{−2σ} ≠ 0, so E is a local diffeomorphism. Properness. If E(τ_n, s_n) stays in a compact set,
then p and q increase along each geodesic, so the starting point (p₀, q₀) = (−1 − s/√2, −1 + s/√2) satisfies p₀ ≤ P,
q₀ ≤ Q, which bounds |s_n|; and (p + q)˙ ≥ 2√(ṗq̇) = √2 e^{−σ} ≥ c > 0 (σ ≤ H₁J₁), so τ_n is bounded. For τ_n < 0 the
geodesics are straight lines of W, with s = (q − p)/√2 and τ = (p + q + 2)/√2, both bounded. A proper local
diffeomorphism ℝ² → ℝ² is a covering map, hence a diffeomorphism since ℝ² is simply connected. This is (b). (Each
slice E(τ, ·) is also embedded directly: it is spacelike, so p decreases and q increases along it.) Exterior. In W the
geodesic is the straight line (p₀ + u, q₀ + u), u = τ/√2. If s ≥ 1/√2 it reaches p = 0 at q = √2 s ≥ 1, and q only
increases, so J′(q) = 0 and K = 0 along the whole geodesic: A ≡ 1, A_t = 0. If s ≤ −1/√2 it reaches q = 0 at p = −√2 s
≥ 1 and H′ = 0 along it. The geodesics with s ≥ 1/√2 sweep a connected region, so b_x = −A_t/A = 0 and b is constant
there at each t; the same for s ≤ −1/√2. For t ≤ 0, A ≡ 1 everywhere and b ≡ 0 in the gauge b_− = 0. This is (c).
Bounds. A_τ ≥ 0 makes b(t, s) = −∫_{−∞}^{s} A_t ds′ nonincreasing in s, so b ∈ [−L′(t), 0]. Since K ≥ 0 and dA = dt dx
in a PG chart, L′(t) = ∫_{−∞}^{t}∫K dx dt′ = ∫K dA over the part of ℝ² below the slice t, which increases to ∫K dA =
2∫H′ dp ∫J′ dq = 2H₁J₁ (`exterior` (iii): K e^{2σ} = 2H′J′, dA = e^{2σ} dp dq, ΔS = 1 for the step of the script), and
equals it as soon as the slice lies above [0, 1]². This is (d). ∎

By the equivalence of [`two_scalars.md`](two_scalars.md) (Setup, stages 1 and 2), R = a a − b b with a = ½(c dH −
(4/c) dJ), b = ½(c dH + (4/c) dJ), any real c ≠ 0 (so k = a + b = c dH, m = a − b = −(4/c) dJ, ½κ₁κ₂ = −2H′J′ = R_pq):
the scalars are φ₁ ∝ cH − (4/c)J and φ₂ ∝ cH + (4/c)J up to normalization, λ = 0, and these metrics are exact
solutions with one normal and one phantom scalar, planar in b(t, x). *Update from #171:* with α = 1/(4√(2π)) the
scalars φ₁ = α(cH − (4/c)J) (normal) and φ₂ = α(cH + (4/c)J) (phantom), V = 0, satisfy G_ab − 8π(T₁ − T₂)_ab = 0
exactly for H, J undefined and every c ≠ 0 (*checked symbolically*, [`planar_residual.md`](planar_residual.md), #171).
The normalizations and V = 0 are forced wherever H′J′ ≠ 0.

**Numerical illustration** (E1 to E5, **numerical evidence**; H = (3/2)S(p), J = (3/2)S(q), S(x) = f(x)/(f(x) + f(1 −
x)), f = e^{−1/x}; RK4, dt = 0.002, 601 geodesics on s ∈ [−6, 6], t ∈ [0, 10]; thresholds in the README of the raw
folder, fixed before the first run):

| check | statement | value |
|---|---|---|
| `E_norm` | 2e^{2σ}ṗq̇ = 1 along every geodesic | max error 2.4e-11 (≤ 1e-8) |
| `E_focus` | A ≥ 1 for all s and t; the lump acts | min A = 1.000000000000; max_s A(10, s) = 109.6 (≥ 2) |
| `E_ext` | A = 1 for \|s\| = 3 … 6 at all stored t | max \|A − 1\| = 0 (≤ 1e-12); A′ ≠ 0 only for \|s\| ≤ 0.60 at t = 10 |
| `E_hol` | L′(10) = 2hj within 1e-5, and b ≥ −2hj − 1e-9 | L′ = 0, 0.98740, 4.49018, 4.5, 4.5, 4.5 at t = 0, 2, …, 10; **FAILS the second part**, §7 |
| `E_hol_b` | the same with Gauss–Legendre in s, error estimates | L′(10) − 4.5 = +2.1e-12, spreads 2.4e-12 and 3.6e-13 |
| `E_geo` | A from the ODE equals e^σ√(−2p_s q_s) from neighboring geodesics | max relative difference 3.0e-7 (≤ 1e-5) |

So b is 0 for x below the lump and drops to −L′(t) above it, with L′ = 0, 0.987, 4.490 and 4.5000 at t = 0, 2, 4 and t
≥ 6; once the slice lies above the diamond the far observers separate at a constant relative rate, which is also the
boost rapidity 2H₁J₁ between the two flat charts of the exterior (the holonomy around the lump). For t ≥ 6, |b| = 4.5
exceeds 1 on the right, so this closure-passing planar shift **does cross |b| = 1**; it is time-dependent, not
traveling, so T′ does not apply.

**Which closure-passing metrics have a localized PG chart?** Sufficient: Theorem E (with the sign condition H′, J′ ≥ 0
and supports in a bounded interval). Excluded: those that are traveling (T′, §2). Not decided: H′J′ of both signs,
where the Jacobi field A may reach 0 (focusing) and the covering argument is not available; supports that are not
bounded.

## 4. Equal limits (L2)

**Proposition L (proved).** If b ∈ C² has b_x → 0 and a → b_±′ at ±∞, then ∫K dx = a_− − a_+ = b_−′ − b_+′ at each t
(`jacobi` (iii) for the example). So b_+ = b_− on an interval I of times requires ∫K dx = 0 on every slice of I.
Consequences:

1. *If K ≥ 0 (H′J′ ≥ 0) then (L2) fails on every interval on which K ≢ 0* (proved). In Theorem E, b_+ − b_− = −L′ ≠ 0
   as soon as the lump is crossed.
2. If the chart covers the null plane, L′(∞) − L′(−∞) = ∫K dA = 2ΔHΔJ (`exterior` (iii), for the step example), so
   (L2) for all t ∈ ℝ needs ΔH·ΔJ = 0, where ΔH = H(+∞) − H(−∞) (by hand, a hypothesis on the chart).
3. At linear order in the amplitude no compactly supported lump has (L2) for all t. Write b = ε b₁ + O(ε²); then C[b]
   = (−∂_t² + ∂_x²) ln|K₁| + O(ε) with K₁ = −b₁,tx (`limits_necessary` (i), exact), so K₁ = F(t − x) G(t + x) at
   leading order, and ∫K₁ dx = (F∗G)(2t) (`limits_necessary` (ii), exact, checked on the two pieces for F = G = u(1 −
   u)(1 − 2u) on [0, 1], where F∗G is a nonzero polynomial on each piece, `limits_necessary` (iii)). If F, G ∈ L¹ have
   compact support, F̂ and Ĝ are entire (Paley–Wiener), so F∗G ≡ 0 gives F̂Ĝ ≡ 0, hence F̂ ≡ 0 or Ĝ ≡ 0, i.e. F ≡ 0 or
   G ≡ 0 (**proved**, by hand). Data that are not compactly supported (band-limited with disjoint spectra) escape
   this.

The full nonlinear question with K of both signs is **not decided**: a lump whose ∫K dx vanishes on every slice would
be a localized bubble with equal limits, and neither a construction nor an exclusion is available here. No argument at
the onset of the lump is made here.

## 5. What is decided

| question | answer | status |
|---|---|---|
| closure criterion for planar shifts | C[b] = □ ln\|K\| − 2K = 0 on {K ≠ 0}; for static u it is Θ′ = 0 | **proved**, **checked symbolically** |
| traveling shift, \|u\| ≠ 1 on ℝ (Theorem T) | fails unless flat | **proved** (#80) |
| traveling shift that crosses \|u\| = 1, any number of times | fails unless flat (Theorem T′), so a superluminal planar bubble of constant speed is excluded | **proved**, steps **checked symbolically** |
| what Θ′ = 0 needs at a simple crossing | α = −uu′(ξ₀), K ~ \|η\|^γ, γ = −(α + u₁)/u₁; two horizons of a bubble give α of opposite sign | **checked symbolically** |
| non-traveling shift, limits at ±∞ (L1), all t | exists: Theorem E, with b_− = 0, b_+ → −2ΔHΔJ | **proved**; the example is **numerical evidence** |
| equal limits (L2), K ≥ 0 | impossible on any interval where K ≢ 0 | **proved** |
| equal limits (L2), linear order, compact support | impossible | **proved** (Paley–Wiener) |
| equal limits (L2), K of both signs, nonlinear | — | **not decided** |
| Theorem E for H′, J′ of both signs, or unbounded supports | — | **not decided** |
| C¹ shifts with K unbounded at a horizon | — | **not decided** (α = 0: u′ → ∞, by hand) |

## 6. Checks and controls

Run of record: `run_final.run1` of script commit 20d1fa4, in the raw folder `I161-rpi5` of negative-mass-raw (stdout
copied to `planar_closure_output.txt`). Checks only: `python3 notes/warp/planar_closure.py --no-controls` (91 s on a
Raspberry Pi 5 at nice 10); one control: `--mutate=NAME`; everything: no argument (checks, then 19 controls, 31 min on
the unit; one process; peak resident 95 MB for the checks, measured by `/usr/bin/time -v`, raw
`run_nocontrols.run1.time`). Exit 1: `crit_pg_iv` and `E_hol` fail as worded (§7), 12 of 14 checks pass.

| check | what it tests | result | s | controls |
|---|---|---|---|---|
| `crit_null` | covariant criterion in null coordinates; 3 planted passes, 3 planted failures vs `closure` | pass | 3 | crit-factor, crit-deform |
| `crit_pg` | PG chart: R = K g, C[u] = Θ′, obstruction Θ′/2, formal planted pass and its perturbation | pass | 28 (one function) | crit-static-g, crit-plant-jet |
| `crit_pg_iv` | Theorem T's case at 5 fixed points, K ≠ 0 and C ≠ 0 (**FAILS**: u₃ at ξ = 1/3) | **FAIL** | (same) | none (fails at baseline) |
| `crit_pg_ivb` | the same at the first 5 points of a fixed list with K ≠ 0 (15 points) | pass | (same) | ivb-flat |
| `horizon` | exponent and α at a crossing, opposite α at the two horizons of u_s, s and Φ on 4 components | pass | 0 | hor-orient, hor-alpha |
| `jacobi` | Gaussian normal chart = PG chart, A_tt = K A, b_x = −A_t/A, ∫K dx = a_− − a_+ | pass | 3 | jacobi-sign, jacobi-theta |
| `exterior` | Christoffels, A = p_s/p_τ and A_ττ = K A from geodesics, step S, ∫K dA = 2hj | pass | 0 | ext-hol, ext-null, ext-geo |
| `limits_necessary` | linearization of C, ∫K₁ dx = (F∗G)(2t), F∗G ≠ 0 | pass | 1 | lin-order, conv-zero |
| `E_norm` | unit normalization along 601 geodesics (numerical) | pass | 6 | E-norm |
| `E_focus` | A ≥ 1 for all s, t; the lump acts (numerical) | pass | 0 | E-kneg |
| `E_ext` | A = 1 for |s| = 3 … 6 (numerical) | pass | 0 | E-ext |
| `E_hol` | L′(10) = 2hj to 1e-5 and b ≥ −2hj − 1e-9 (**FAILS** the second part, §7) | **FAIL** | 0 | none (fails at baseline) |
| `E_hol_b` | the same with Gauss–Legendre in s, RK4 dt = 0.001, printed spreads | pass | 47 | E-hol |
| `E_geo` | A from the ODE vs the geometric separation of neighboring geodesics | pass | 2 | E-geo |

| control `--mutate=` | change | target | outcome | failed checks (baseline failures `crit_pg_iv`, `E_hol` omitted) |
|---|---|---|---|---|
| `crit-factor` | expects d_p d_q ln\|R_pq\| = -e^(2 sigma) C | `crit_null` | fired, exit 1 | `crit_null` |
| `crit-deform` | adds p q^2/10 to sigma of the planted passes | `crit_null` | fired, exit 1 | `crit_null` |
| `crit-static-g` | g^xx = 1 + u^2 in box | `crit_pg` | fired, exit 1 | `crit_pg` |
| `crit-plant-jet` | the planted pass jets get c6 + 1/7 | `crit_pg` | fired, exit 1 | `crit_pg` |
| `ivb-flat` | u1 replaced by the constant 1/3 (K = 0 at every point of the list) | `crit_pg_ivb` | fired, exit 1 | `crit_pg_ivb` |
| `hor-orient` | s' = 2/(1 + u^2) | `horizon` | fired, exit 1 | `horizon` |
| `hor-alpha` | alpha = +u u'(xi0) | `horizon` | fired, exit 1 | `horizon` |
| `jacobi-sign` | K = -A_tt/A | `jacobi` | fired, exit 1 | `jacobi` |
| `jacobi-theta` | b_x = +A_t/A | `jacobi` | fired, exit 1 | `jacobi` |
| `ext-hol` | expects integral of K dA = h j | `exterior` | fired, exit 1 | `exterior` |
| `ext-null` | Gamma^p_pp = sigma_p | `exterior` | fired, exit 1 | `exterior` |
| `ext-geo` | p'' = -sigma_p p'^2 in the derivation of A_tt = K A | `exterior` | fired, exit 1 | `exterior` |
| `lin-order` | K1 = -b1_tt | `limits_necessary` | fired, exit 1 | `limits_necessary` |
| `conv-zero` | G = 0 in the convolution | `limits_necessary` | fired, exit 1 | `limits_necessary` |
| `E-norm` | initial velocity (1, 1) | `E_norm` | fired+extra, exit 1 | `E_norm`, `E_hol_b`, `E_geo` |
| `E-kneg` | A'' = -K A | `E_focus` | fired+extra, exit 1 | `E_focus`, `E_hol_b`, `E_geo` |
| `E-ext` | support of J' is [0, 2] | `E_ext` | **silent**, exit 1 | none besides the baseline failures |
| `E-hol` | expects L'(10) = h j | `E_hol_b` | fired, exit 1 | `E_hol_b` |
| `E-geo` | finite-difference spacing 0.5 | `E_geo` | fired, exit 1 | `E_geo` |
| `E-ext5` (supplementary run, `run_ctl_E-ext5.run1.log`, commit 0707d9c) | support of J' is [0, 5] | `E_ext` | fired+extra, exit 1 | `E_ext`, `E_hol_b` |

Outcomes of the 19 controls of the run of record: 16 fire on their target and on no other check; `E-norm` and `E-kneg`
fire on target with extras, `E_hol_b` and `E_geo` (a unit-speed normalization and A″ = KA are inputs of both; the
README named one target per control and said extras are reported, so none was predicted); **`E-ext` is silent**.
Cause, computed afterwards: I had guessed that the geodesics with |s| = 3 … 6 enter the support of K when the support
of J′ is [0, 2]. A normal geodesic of S₀ reaches p = 0 at q = √2 s, so it meets the support only if √2 s is below its
end, i.e. for s < 1.41, and `E_ext` tests |s| = 3 … 6 only. The control `E-ext5` (support [0, 5]; in a scratch copy it
gave max |A − 1| = 1.61e-2) was added in the README (A8), committed alone (0707d9c) and run alone: it fires on target,
with `E_hol_b` as an extra. The two baseline failures have no control, since their targets already fail.

## 7. Criteria that fail as written, and defects fixed

**`crit_pg_iv` (C2 (iv)) FAILS.** The README fixed five points (ξ = 1/3, 1, 3/2, 2, 5/2) and asked for K ≠ 0 and C ≠ 0
at each for the three profiles u₁, u₂, u₃ of §2. u₁ and u₂ pass at all five. u₃ = 1/(3(1 + ξ²)²) has K = 4(9ξ² −
1)/(9(1 + ξ²)⁶), which is 0 at ξ = 1/3, so the criterion fails there. The reason: the wording did not exclude zeros of
K, where Φ = 0 and Θ is undefined. By the coordinator's ruling the check stays as worded and stays failed, the profile
and the points were not moved, and one further check was allowed with a rule fixed before its run: `crit_pg_ivb`
(README A3) passes. (iv) has no control, since its target fails at baseline.

**`E_hol` (E4) FAILS.** The README asked |L′(10) − 2hj| ≤ 1e-5 (passes: +4.79e-6) and b ≥ −2hj − 1e-9 (fails: b_min =
−4.5000048). The diagnosis (scratch study, raw `diag_E_hol.run1.log`, not a check): the excess is quadrature error of
the trapezoid sum in s, not the time step and not a violation of the bound. L′(10) − 4.5 against the s spacing at dt =
0.002:

| s spacing | 0.02 | 0.01 | 0.005 | 0.0025 | 0.005 with dt = 0.001 |
|---|---|---|---|---|---|
| L′(10) − 4.5 | +4.793e-6 | +8.047e-9 | −5.214e-12 | −5.963e-12 | +4.334e-13 |

The excess falls faster than any power of the spacing, as for a C^∞ integrand of compact support; the 5e-12 left is
the RK4 time-step error (it changes with dt, not with the grid). Domain and ends are not the cause: [−3, 3] gives the
same values as [−6, 6], and A′ = 0 exactly for |s| > 0.64 at t = 10. The true value L′(10) = 2hj is the total
curvature, by a Gauss–Bonnet-type hand argument: L′(t) = ∫K dA over the part of the plane below the slice (§3), K
e^{2σ} = 2H′J′, so the total is 2ΔHΔJ = 2hj, and it is reached once the slice lies above the diamond. The coordinator
refused a check with a slack of 1e-5 (a threshold picked after seeing the excess) and accepted `E_hol_b` (README A6,
A7), with the README's own slack 1e-9, composite Gauss–Legendre in s (order 10, 48 panels; 96 for the estimate) and
RK4 at dt = 0.001: L′(10) − 2hj = +2.1e-12, with error estimates e_q = 2.4e-12 (panels 48 vs 96) and e_t = 3.6e-13 (dt
0.001 vs 0.0005). These estimates are spreads over two resolutions, not enclosures. `E_hol_b` passes. `E_hol` stays
failed as worded.

**Defects fixed, logs kept** (raw `dev_*`): a helper compared a matrix with 0 (`dev_crit_null.run1.log`, C1 failed,
rerun passes); a `limit` with a symbolic exponent hung C3 (`dev_pg_horizon.run1.log`, stopped by hand after 13 min,
10:47 CPU minutes, 670 MB, no result); the leading-order coefficient is now written by hand and C3 takes 1 s.

The script exits 1 because of `crit_pg_iv` and `E_hol`.

## 8. What the code proves, and what is assumed

Exact (sympy, undefined functions or exact rationals): R = K g₂, K = 2σ_pq e^{−2σ}, □f = −2e^{−2σ}f_pq and the
identity d_p d_q ln|R_pq| = −(e^{2σ}/2) C for an undefined σ; the six planted σ; K = −a_x for a PG metric with
undefined b(t, x); C[u] = Θ′ and dω_{tξ} = Θ′/2 for undefined static u; the 15 + 15 exact points of `crit_pg_iv` and
`crit_pg_ivb`; the formal planted pass and fail; the leading order at a simple crossing and α = −uu′; the residues and
limits on the four components; R = (A_tt/A) g for −dt² + A² ds², the chain-rule identities of the PG chart of a normal
chart, the limits of a with b → b_±; the Christoffel symbols of −2e^{2σ} dp dq, A = p_s/p_τ and A_ττ = K A from the
geodesic equations on jets, the step S and ∫K dA = 2hj; the linearization of C and the convolution identity for F = G
= u(1 − u)(1 − 2u).

Numerical evidence: E1 to E5 and `E_hol_b` (floating point, RK4, quadrature in s, printed spreads).

By hand: Theorem T′, Theorem E (geodesic completeness, properness and the covering argument, the exterior, the
bounds), Proposition L and its consequences, the Paley–Wiener step, the identification of the scalar fields. The check
names in §§2–4 mark which steps the code supports.

Imported and trusted: `ricci` (`crit_null`, `crit_pg`, `jacobi`, `exterior`) and `closure`, `dform` (`crit_null`,
`crit_pg`) from the code of #57 and #65, reviewed in #59 and #79; `closure` and the null-coordinate form of #80 are
not independent of #65. The planted passes of `crit_null` are the ones of #80.

Assumed: two real, minimally coupled scalars of opposite sign with λ = 0; planar b(t, x); u ∈ C² in Theorem T′ (C³
where K ≠ 0); H, J ∈ C^∞ with H′, J′ ≥ 0 and supports in [0, 1] in Theorem E; sympy 1.14.0, numpy 2.5.4, Python 3.14.

Not tested: Theorem E for other H, J (one example is integrated); K of both signs; the nonlinear (L2) case; any
traveling shift with K unbounded at a horizon.

## 9. Lemmas for Lean (#135)

Statements only; ℝ throughout.

1. **Closure criterion** (`field_simp; ring` after clearing e^{2σ}). For σ(p, q) with σ_pq ≠ 0, d_p d_q ln|σ_pq| =
   −(e^{2σ}/2)(□ ln|K| − 2K) with K = 2σ_pq e^{−2σ}, □f = −2e^{−2σ}f_pq.
2. **Static reduction** (`ring`). With K = (uu′)′, ((1 − u²)(ln|K|)′)′ − 2K = Θ′, Θ = (1 − u²)K′/K − 2uu′.
3. **Jacobi from geodesics** (`ring` on jets). With p_ττ = −2σ_p p_τ², 2e^{2σ}p_τq_τ = 1 and p_s q_τ + p_τ q_s = 0,
   the quantity A = p_s/p_τ satisfies A_ττ = 2σ_pq e^{−2σ} A.
4. **Mechanism of T′** (elementary analysis). If s : (a, b) → ℝ is continuous and strictly monotone with |s| → ∞ at
   both ends, Φ = C e^{αs}, C ≠ 0, then Φ does not tend to 0 at both ends.
5. **Exponent at a crossing** (`ring` on series). The leading order of (1 − u²)K′ = (2α + 2uu′)K at u = 1 + u₁η + …
   for K = κη^γ gives γ = −(α + u₁)/u₁.
6. **Covering step of Theorem E**. A proper local diffeomorphism ℝ² → ℝ² is a bijection (covering space theory).
7. **Total curvature.** ∫∫ 2H′(p)J′(q) dp dq = 2(H(∞) − H(−∞))(J(∞) − J(−∞)) for H′, J′ ∈ L¹ (Fubini).
8. **Paley–Wiener step.** For F, G ∈ L¹ with compact support, F∗G = 0 implies F = 0 or G = 0 (entire functions have no
   zero divisors; the Fourier transform is injective on L¹).

## 10. What `det_classify.md` and `two_scalars.md` should now say

*Not edited here; for the coordinator.*

`det_classify.md` §2, last paragraph ("The theorem covers every subluminal planar bubble … is not decided."), could
read: "Theorem T covers traveling planar shifts with |u| ≠ 1 everywhere. [`planar_closure.md`](planar_closure.md)
(#161, Theorem T′, *proved*) extends it: a traveling planar shift with limits and bounded K fails closure unless flat,
whether or not it crosses |u| = 1, so a superluminal planar bubble of constant speed is excluded. Non-traveling planar
shifts: closure-passing metrics −2 exp(2H(p)J(q)) dp dq with H′, J′ ≥ 0 of bounded support have a global PG chart with
b = 0 for t ≤ 0, b constant in x outside a compact set at every t, and unequal limits b_− = 0, b_+ → −2ΔHΔJ (Theorem
E, *proved*; the example is *numerical evidence*). With equal limits at every t the question is *not decided*; it is
impossible for K ≥ 0 and at linear order for compact data."

`det_classify.md` §12 ("Not covered"): replace "non-traveling planar shifts with closure (§2); traveling planar shifts
that cross |u| = 1 (§2)" by "planar shifts with closure and equal limits at every t, with K of both signs (#161)".

`two_scalars.md`, "Not covered", the planar sentence: "Traveling planar shifts fail closure unless flat, whether or
not they cross |u| = 1 (*proved*, [`planar_closure.md`](planar_closure.md)). A non-traveling planar shift b(t, x) from
a pair −2 exp(2H(p)J(q)) dp dq with H′, J′ ≥ 0 passes closure and is localized in x with unequal limits (*proved*);
with equal limits it is *not decided*."

Nothing in #158's Theorem P changes: it concerns constant K, which closure of the one-sign pair allows, and here K is
not constant.



## Coordinator's review

Merged as dccb4a0. Raw negative-mass-raw `I161-mac` at 773f663 (criteria at 29e51f2, before the first run).

| step | result |
|---|---|
| scope | the three files the brief names, and raw `I161-rpi5` |
| rerun on the Mac | exit 1 as ruled: `crit_pg_iv` and `E_hol` fail as worded, 12 of 14 checks pass, `E-ext` silent, `E-ext5` fires; output equals rpi5's apart from timings and the `E-ext5` line |
| hand proofs | Theorem T′, Theorem E (geodesic completeness, the Jacobi field, properness and covering, the exterior, the bounds), Proposition L and the Paley–Wiener step read and rederived; no gap found |
| re-solve | Theorem E's example, L′(10) = 2H₁J₁ = 4.5, by the coordinator's own geodesic and Jacobi integration |
| independent check `ind161.py` | run 1 exits 1 (Q4 fails on one clause, below); 6 of 6 mutations fire on target |

ind161 shares no code with this note's script, `two_scalars.py` or `phantom_source.py`.

| id | what | result |
|---|---|---|
| Q0–Q2 | own Ricci of −2e^{2σ}dp dq and of the PG metric: K = 2σ_pq e^{−2σ}, K = −a_x; the closure criterion identity; C[u] = Θ′ | pass |
| Q3 | u_s = 3/(1 + ξ²): horizons ±√2 with α = ∓2√2/3; on each component s′ one-signed, s divergent at both ends, Φ → 0; γ = −(α + u₁)/u₁ | pass |
| Q4 | own integration of Theorem E's example: L′(10) = 4.4999999998 (within 1e-9), A ≥ 1, A_τ ≥ 0 | **fails** on the normalization clause: 1.6e-8 against the coordinator's 1e-10, which is RK4 truncation at the fixed steps (error ∝ dt⁴, diagnostic in the raw folder); the L′ and sign clauses pass |
| Q5 | the linearization of C; ∫K₁ dx = (F ∗ G)(2t) at three exact t; F ∗ G ≢ 0 | pass |

**Discrepancies and process.**
- `crit_pg_iv` and `E_hol` failed as worded and stay failed, by the coordinator's rulings during the run.
  `crit_pg_ivb` and `E_hol_b` were added under rules fixed before their runs, and `E_hol_b` kept the README's 1e-9
  slack. `E-ext` was silent because its effect was guessed; `E-ext5` replaced it after its effect was computed.
- The coordinator's own Q4 failed in the same way as `E_hol`: a tolerance tighter than the integration step allows. It
  stays failed, with its diagnostic.
- Not computed: the residual stage G − 8πT = 0 with the explicit scalars of Theorem E (§3, by hand). *Update from #171:* computed, zero exactly.

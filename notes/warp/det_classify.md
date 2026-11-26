# Which shifts have det R₃ = 0 on an open set?

November 2026, issue [#80](https://github.com/nmoinvaz/negative-mass/issues/80). Script
[`det_classify.py`](det_classify.py), its stdout [`det_classify_output.txt`](det_classify_output.txt). It builds on
[`two_scalars.md`](two_scalars.md) (#65, reviewed in #79 by [`review_two_scalars.md`](review_two_scalars.md)), imports
`ricci`, `ricci_at_point` from [`phantom_source.py`](phantom_source.py) and `closure`, `two_scalar_form`, `dform` from
[`two_scalars.py`](two_scalars.py), and re-derives D itself.

Status words as in AGENTS.md: **proved** (by hand), **checked symbolically**, **numerical evidence**, **conjecture**,
**not decided**.

## Answer

D = det R₃ vanishes on an open set for many shifts, so "classify" is done inside reductions. What is decided, and what is not:

| reduction | solutions of D = 0 on an open set | localized in x? | closure (ω) | status |
|---|---|---|---|---|
| planar b(t, x) (b_ρ ≡ 0) | **all of them**: D ≡ 0 | yes, e.g. b = −v f(x − vt) | traveling and \|u\| ≠ 1 everywhere: fails unless flat; other cases: see §2 | **proved**, **checked symbolically** |
| b(ρ) of #65 and its linear-in-x deformation u = c x + g(ρ) | g = −ε ln(α J₀(√2 c ρ) + β Y₀(√2 c ρ)); c = 0 is #65 | no (linear or constant in x) | fails on every open set, every member (#162) | **proved**, **checked symbolically** |
| additive b = X(x) + R(ρ) | exactly planar, or u = c x + d + g(ρ) as above | only the planar ones | as above | **proved** (hand argument), **checked symbolically** |
| multiplicative b = X(ξ) R(ρ) | none with X X′ R R′ ≠ 0, for every v (#160); X′ ≡ 0 is #65's b(ρ), R′ ≡ 0 is planar | — | not needed | **proved** in [`mult_classify.md`](mult_classify.md), steps **checked symbolically** (§5) |
| conical b = F(x/ρ) | ODE of second order, cubic in F″ | **not decided**: "no limit at η → ±∞" is not established (§6) | not run | ODE **checked symbolically** |
| static or traveling b(x, ρ), general | one second-order PDE, many local solutions | only inside the rows above | not run | not classified |

So the answer to the issue's question is: **yes**, shifts localized in x with D = 0 on an open set exist, and the ones found are the planar
ones, which are constant in ρ and so are not bubbles. Within the additive-separable reduction no shift that depends on ρ is
localized in x (**proved**). For the other reductions whether a shift with u_ρ ≢ 0 is localized in x is **not decided**.
The traveling bubble of constant speed passes on no open set (planted truth, §8).

Two criteria of the README failed and are reported as results: `bessel_pert` (§3) and `conical_limit` (§6). A third, the
hand formula for the degree-3 part of the conical ODE, failed first and was replaced by a check against a hand derivation (§6).

## Setup

ds² = −dt² + (dx + b dt)² + dy² + dz², ρ² = y² + z², n = ∂_t − b ∂_x. #65's block of the frame Ricci tensor R₃ on
(n, ∂_x, e_ρ) is

    R₃ = [[P − h, Q, Z], [Q, −P − h, Y], [Z, Y, h]],
    P = n(b_x) − b_x², Q = −½(b_ρρ + b_ρ/ρ), Z = ½ b_xρ, Y = −½ n(b_ρ) + b_x b_ρ, h = ½ b_ρ²,
    D = det R₃ = h(h² − P² − Q² + Y² + Z²) − P(Y² − Z²) + 2QYZ.

With a = b_t − b b_x (the comoving acceleration) the same entries are P = a_x and Y = ½(b_x b_ρ − a_ρ) (`d_block`: **checked
symbolically** against the coordinate Ricci tensor of an undefined b(t, x, ρ), with R(e_φ, ·) = 0). By #65's lemma
(**proved**; reviewed in #79), R = λη + a a − b b with λ real holds iff λ = 0 and D = 0. So "D = 0 on an open set" is the algebraic stage of
the normal-plus-phantom scalar test on that set, and closure (ω, [`two_scalars.md`](two_scalars.md) stage 2) decides the rest.

**Static and traveling shifts are one class** (`travel`, **checked symbolically**, polynomial identities in the jets). For
b = B(x − vt, ρ) the jets are b_t = −v B_ξ and so on, and D(b) = D_static(u) with u = B + v. Also b′(t, x, ρ) = b(t, x + ct, ρ) + c has the same D
(Galilean covariance). A bubble at constant speed v is a static shift u that tends to v as |x| → ∞.

**Localization**, my definition for this issue: on a strip ℝ × I, I an interval in ρ > 0, u(·, ρ) has limits u_±(ρ) as x → ±∞ for
each ρ in I (compact support in x is u_± = 0), and u_x ≢ 0. A shift is *transversally nontrivial* if u_ρ ≢ 0. "Bubble-like" means
both. Smoothness: C² where D is evaluated, C³ in the planar theorem, C^∞ not needed anywhere.

## 1. What a local solution looks like

Solving D = 0 for Q gives a quadratic with the leading coefficient −h, so where h ≠ 0 the equation can be solved for
b_ρρ in terms of the other jets. Cauchy–Kovalevskaya then gives analytic local solutions for analytic data b and b_ρ on ρ = const, two free
functions of (t, x) (of x for a static shift), where h ≠ 0 and the discriminant is positive (the other jets in D have ρ-order at most 1). So the full set of local solutions is not a finite-dimensional
family, and no classification without a reduction is attempted. (This paragraph is not a theorem about global solutions; it is
why the reductions below are the only way to say anything exact.)

## 2. The planar family b(t, x): a counterexample family

For b_ρ ≡ 0: h = Q = Z = Y = 0, so R₃ = diag(P, −P, 0) and **D ≡ 0 for every b(t, x)** (**proved**; `planar`). It is the product of a 2D
metric g₂ = −dt² + (dx + b dt)² with ℝ², so R = K g₂ ⊕ 0 with K = −P = (u u′)′ in the static case. Where K ≠ 0, R has rank 2
and signature (1, 1): the algebraic stage passes with λ = 0 (`planar`: 8 exact points of 3 rational planar shifts through
`two_scalar_form`, status "form", λ = 0; the same shifts plus ρ²/10 fail it with D ≠ 0 at all 8 points).
Every shift b(t, x) with compact support or limits in x is *x-localized* and has D ≡ 0. They are not bubbles: u_ρ ≡ 0.

**Closure for planar shifts.** In null coordinates g₂ = −2e^{2σ} dp dq, R_pq = −2σ_pq (the others vanish). Then (**proved**, and
`planar_closure` (i) is a planted test with three passes and three failures) `closure` returns ok iff

    ∂_p ∂_q ln|σ_pq| = 0,  i.e.  σ_pq = F(p) G(q),  i.e.  σ = H(p) J(q) + f(p) + g(q),  H′ = F, J′ = G.

Reparametrizing p and q absorbs f and g, so the planar closure-passing metrics are exactly g₂ = −2 e^{2H(p)J(q)} dp dq, two
arbitrary functions of one variable. Proof of the equivalence: k and m with R = ½(k ⊗ m + m ⊗ k) must be the two null covector
fields, k = κ₁ dp, m = κ₂ dq, and e^χ k, e^{−χ} m are closed iff e^χ κ₁ = A(p) and e^{−χ} κ₂ = B(q), which needs κ₁κ₂ =
2R_pq = A(p)B(q). This is a local statement.

**Traveling planar shifts.** For a static planar shift u(ξ) (a traveling shift, §Setup) the null coordinates are p = t − A(ξ), q = t + B(ξ),
A′ = 1/(1 − u), B′ = 1/(1 + u), with −(1 − u²) dp dq equal to the metric (`planar_traveling` (a), **checked symbolically**), and R_ab = (uu′)′ g_ab,
so R_pq = −Φ/2 with Φ = (uu′)′(1 − u²) (`planar_traveling` (b)). The obstruction of the ω construction is

    dω_{tξ} = ½ Θ′,  Θ = (1 − u²) Φ′/Φ     (`planar_closure` (ii), exact for an undefined u),

so closure holds on an interval of {K ≠ 0, |u| ≠ 1} iff Θ is constant, Θ = 2α, i.e. Φ = c e^{αs} with s′ = 2/(1 − u²).

**Theorem T (planar, traveling; proved).** Let u ∈ C³(ℝ) with |u| ≠ 1 on ℝ, limits u_± as ξ → ±∞ with |u_±| ≠ 1, and u′, u″ → 0.
If closure holds on {K ≠ 0}, then K ≡ 0, i.e. u is constant.

*Proof.* Suppose K ≢ 0 and let J = (ξ₁, ξ₂) be a maximal interval of {K ≠ 0}. On J, Φ = c e^{αs} ≠ 0. If ξ₂ is finite, K(ξ₂) = 0 by
continuity and maximality, so Φ(ξ₂) = 0; but s is continuous at ξ₂ (|u| ≠ 1), so Φ → c e^{αs(ξ₂)} ≠ 0, a contradiction. The same at ξ₁. So J = ℝ. Then
Φ → 0 at both ends (K → 0, 1 − u² bounded). Since s′ → 2/(1 − u_±²), s ~ 2ξ/(1 − u_±²). If |u| < 1 on ℝ then s′ > 0, s → ±∞ at ±∞, and
c e^{αs} → 0 at both ends needs α < 0 and α > 0. If |u| > 1 on ℝ then s′ < 0 and the same contradiction. α = 0 gives Φ = c ≠ 0, not → 0. ∎

The theorem covers every subluminal planar bubble (|u| ≤ v < 1) and every traveling planar shift that stays on one side of |u| = 1. It does not cover
a shift that crosses |u| = 1 (a superluminal bubble crosses it twice, at its horizons), and it does not cover non-traveling planar shifts b(t, x).
For those, closure-passing local solutions exist (the family above), and **whether one is x-localized is not decided**. As a planted
failure through the same code, the 1+1 bubble u = ½(1 − 1/(1 + ξ²)²) has Θ′ ≠ 0 at 5 exact points (`planar_closure` (iii)); a jet with Θ′(0) = 0 passes.

**Update from #161** ([`planar_closure.md`](planar_closure.md)): Theorem T′ (*proved*) extends Theorem T. A
traveling planar shift with limits |u_±| ≠ 1 and K bounded fails closure unless flat, whether or not it crosses
|u| = 1, so a superluminal planar bubble of constant speed is excluded. Non-traveling planar shifts: the
closure-passing metrics −2 exp(2H(p)J(q)) dp dq with H′, J′ ≥ 0 of bounded support have a global PG chart with b = 0
for t ≤ 0, b constant in x outside a compact set at every t, and unequal limits b_− = 0, b_+ → −2ΔHΔJ (Theorem E,
*proved*; the example is *numerical evidence*). With equal limits at every t the question is *not decided*; it is
impossible for K ≥ 0, and at linear order for compactly supported data.

## 3. The family u = c x + g(ρ), which contains #65's b(ρ)

For u = c x + g(ρ): D = g′²(2c² + g′² − S)(2c² + g′² + S)/8 with S = g″ + g′/ρ (`bessel` (a), **checked symbolically**, exact), so where g′ ≠ 0, D = 0 iff

    g″ + g′/ρ = ε (g′² + 2c²),  ε = ±1.

With g = −ε ln ψ this is ψ″ + ψ′/ρ + 2c²ψ = 0 (**proved**: R″ + R′/ρ = −εψ″/ψ + ε(ψ′/ψ)² − εψ′/(ρψ), and substituting the Bessel equation gives
ε(g′² + 2c²)), so

    u = c x − ε ln(α J₀(√2 c ρ) + β Y₀(√2 c ρ)) + d,   and for c = 0   u = −ε ln(C₁ + C₂ ln ρ) + d,

which is #65's b(ρ) family (its b = −ε ln|C − ε ln ρ| is C₁ = C, C₂ = −ε). D = 0 at 36 points of 6 parameter sets (c ∈ {½, 3/2, 0}, ε = ±1, 6 points each, random rational α, β), with the largest |D| 8.9e-117 at 130 working digits (`bessel`).

R = a a − b b on the ODE (`bessel` (c), exact), with a = (g′ dρ + 2c θ¹)/√2, b = √(g′² + 2c²) ((1 + εu) dt + ε dx)/√2, θ¹ = dx + u dt, for ε = ±1.
So the algebraic stage passes with λ = 0 on an open set: the normal gradient a is spacelike and the phantom gradient b is null (ν = (1 + εu)dt + ε dx is null). **Closure fails** (`bessel` (d)):
at the exact jets (c, g′(ρ₀)) = (2, 1), (6, 17), (0, ½), at ρ₀ = 5/7, x₀ = 1/3, and both signs of ε, the two 2-form equations ω ∧ k = −dk and ω ∧ m = dm of the ω construction have
no common solution, with the ODE substituted for g″ and g‴ (**checked symbolically**, at a point). Any jet (g, g′) at ρ₀ lies on a member of the family, so closure fails at a point of a member for each of
these parameter values. For c = 0 this agrees with #65's Frobenius argument (**proved** there). A Frobenius argument for c ≠ 0, which would give failure on every member, was not made.

**Update from #162** ([`bessel_closure.md`](bessel_closure.md)): decided. For every member with c ≠ 0, every ε, α, β,
closure fails on every open set (*proved*, identities *checked symbolically*). The null covectors k = a + b, m = a − b
of R have k̃ ∧ dk̃ = O_k dt ∧ dx ∧ dρ and m̃ ∧ dm̃ = O_m dt ∧ dx ∧ dρ with O_k + O_m = 2g′(g′² + 4c²), so the pair is
Frobenius-integrable only where g′ = 0, and there g″ = 2εc² ≠ 0 makes the zeros isolated. For c = 0 this is #65's
ν ∧ dν = b′ dx ∧ dρ ∧ dt.
These shifts are linear in x for c ≠ 0 and independent of x for c = 0. Neither is localized in x.

**`bessel_pert`: a failed criterion, kept.** The README asked, as a second part of (b), that for ψ + ρ²/10 the value |D| exceed 1e-3 at the same points, as a sensitivity test.
It fails: the smallest |D| over the 36 perturbed points is 9.06e-6, at c = 0, ε = −1, ρ = 23/7. The threshold was not changed, by the coordinator's ruling, and the failure gates the exit code.
Observation, not a criterion: the unperturbed |D| is below 1e-50 at all 36 points (largest 8.9e-117). The small values do not come from small ρ. They are at c = 0 at every ρ tested (1.6e-5 to 8.2e-5 for ε = +1,
9e-6 to 3e-3 for ε = −1) and at the largest ρ for c = ½ (6.1e-5 at ρ = 23/7 for ε = +1), while the points at ρ = 3/7 give 2e-5 and 3e-3 (c = 0), 2e-3 and 8e-3 (c = ½) and 0.7 and 2.3 (c = 3/2). I expect, and did not test, that the value is small where
h = ½ g′² is small, since D carries the factor h, and g′ = −ε ψ′/ψ is small for c = 0 and at large ρ.

## 4. Additive separable shifts b = X(x) + R(ρ)

**Theorem A (proved; `additive` checks each algebraic step).** Let u = X(x) + R(ρ) be C² on a box I_x × I_ρ ⊂ ℝ × (0, ∞) with D = 0, and R′ ≠ 0 on I_ρ.
Then X″ ≡ 0, and R″ + R′/ρ = ε (R′² + 2c²) with c = X′, so u is the family of §3. If R′ ≡ 0 on I_ρ, u = X(x) is planar and D ≡ 0 for every X.
Consequently, on the open set {R′ ≠ 0}, no additive-separable shift is localized in x; the localized ones are the planar ones.

*Proof.* With h = R′²/2, D = h·E where E = (h + X′²)² − ((X + R)X″)² − Q² (`additive` (a): exact identity in the jets; P + X′² = −(X + R)X″, Y² = 2X′²h, Z = 0).
Since h ≠ 0, E = 0, and

    E = f₁g₁ + … + f₅g₅,  f = (1, X′², X″², X X″², X′⁴ − X²X″²),  g = (A² − B, 2A, −R², −2R, 1),  A = R′²/2, B = (R″ + R′/ρ)²/4   (`additive` (b)).

(1) g₃, g₄, g₅ are independent since R is not constant (their Wronskian is 2R′³ ≠ 0, `additive` Wronskian check), so the space L of constants (c₁, …, c₅) with Σ c_i g_i ≡ 0 on I_ρ has dimension ≤ 2.
(2) For each x, f(x) ∈ L, so W = span{f(x)} has dimension ≤ 2, and ≥ 1 because f₁ ≡ 1.
(3) If dim W = 1, f is constant, X′² and X″² are constant, so X″ = 0 (X′ is continuous).
(4) If dim W = 2, pick x₀, x₁ with f(x₁) ≠ f(x₀) and β = f(x₁) − f(x₀) ∈ L, β₁ = 0. Then f(x) = f(x₀) + λ(x)β. Σ β_i g_i = 0 reads β₂R′² = β₃R² + 2β₄R − β₅. If β₂ = 0 then
β₃ = β₄ = β₅ = 0 (R nonconstant) and β = 0, so β₂ ≠ 0 and λ = (X′² − α₂)/β₂ is affine in y = X′². Then X″² = γ₀ + β₃y, XX″² = δ₀ + β₄y, y² − X²X″² = ε₀ + β₅y. On a subinterval where X″ ≠ 0 this
gives X = (δ₀ + β₄y)/(γ₀ + β₃y) and the polynomial identity (y² − ε₀ − β₅y)(γ₀ + β₃y) = (δ₀ + β₄y)² on the interval of values of y. Its y³ coefficient gives β₃ = 0, its y² coefficient γ₀ = β₄². So X″² = γ₀ =: m² is constant, X″ = s m, and X is quadratic.
For X = s m x²/2 + c₁x + c₀, y = (c₁ + s m x)² and X = (s/(2m)) y + const, so XX″² = m² X has slope β₄ = s m/2 in y, so γ₀ = β₄² = m²/4 against γ₀ = m², hence m = 0, a contradiction with X″ ≠ 0. In the script this is a
Gröbner basis computation: the equations together with m t = 1 generate the unit ideal for both signs of s (`additive`, elimination step). So dim W = 2 is impossible.
(5) With X″ = 0, X′ = c, E = ((R′² + 2c²)² − S²)/4 (`additive`), and E = 0 is the equation of §3. ∎

## 5. Multiplicative separable shifts b = X(ξ) R(ρ), u = X R + v

The reduced equation is 8D/R′² = E(X, X′, X″, R, R′, S; v), a polynomial in the jets printed in the script (`multiplicative` (a), **checked symbolically**; v ≠ 0 is the traveling case, v = 0 the static one).
For v = 0 and the log derivatives p₁ = X′/X, p₂ = (ln X)″, c₁ = R′/R, s₂ = R″/R + R′/(ρR): E = R²X⁴ [X²R² A + B] with T = 5p₁² + c₁² − 2p₂, A = T(9p₁² + 2c₁² − T), B = p₁²T − (s₂ + 3p₁²)² (`multiplicative` (b)).
The subcases R′ ≡ 0 (planar) and X′ ≡ 0 (b(ρ) of #65) are recovered (`multiplicative` (c)). Eliminating the exponential factor e^{2ψ+2χ} = −B/A leaves two conditions, one that fixes ψ‴ and one that
fixes (s₂)′, and their mutual compatibility gives *one* polynomial relation of degree 17 in (p₁, p₂, c₁, s₂, ρ) (after removing the factors c₁, p₁ and the two factors of A), because the numerators of the two compatibility conditions are equal (`multiplicative` (d)). I did not finish the prolongation.
**The multiplicative case is not classified.** A numerical formal-series search (README C8 (e)) was written and then withdrawn before any run of record (README A5): its converged starts were all of tiny amplitude, so it could not tell a solution from a small
shift, also on the additive equation that has solutions. The scratch evidence is in the raw folder (`scratch_series_search.log`). A properly normalized search is a recommended task.

**Update from #160** ([`mult_classify.md`](mult_classify.md)): decided. For every real v there is no solution of
D = 0 with X X′ R R′ ≠ 0 on an open box (*proved*, steps *checked symbolically*), by separating E into a sum of 6, 10
or 9 products and bounding the dimensions of the two spans. The relation of degree 17 above, which #160 counts as
degree 16 in (p₁, p₂, c₁, s₂, 1/ρ), was not needed, and the prolongation ideal itself was not computed.

## 6. Conical shifts b = F(x/ρ)

D = N(F, F′, F″; η)/(8ρ⁶) with η = x/ρ (`conical`, exact). N is cubic in F″ and printed in the script. The dilation (t, x, ρ) → s(t, x, ρ) is the only symmetry of D = 0 besides translations and Galilean boosts, so F(x/ρ) is the natural similarity form of a static shift. *Added by the coordinator:* no symmetry computation backs
"the only symmetry"; read it as a conjecture.
At fixed ρ the shift tends to limits as x → ±∞ iff F has limits as η → ±∞. The shift is singular at ρ = 0, so these would be open-set solutions only.

**A failure on record.** In the README I wrote the degree-3 part of N in (F′, F″) as 2F w²(F²F″ − ηw), w = F′ + ηF″, and its lowest part with F, F′, F″ all small as −w²((η² + 1)(η² + 2)F′² + 2ηFw). Both are wrong. The formulas came from the script's printed N first, and my
hand simplification dropped a term. The check C9 (b), (c) as worded failed in `dev_mid.run1.log`; that wording stands as the record in the README and the failed log is kept. By the coordinator's ruling it is replaced by C9b, written after the failure.

**C9b (i), re-derived by hand.** With u = F(η), η = x/ρ: u_x = F′/ρ, u_ρ = −ηF′/ρ, u_xx = F″/ρ², u_xρ = −w/ρ², u_ρρ + u_ρ/ρ = η w/ρ² (so Q = −ηw/(2ρ²)). So, times ρ²,

    P̂ = −(F′² + F F″),  Q̂ = −ηw/2,  Ẑ = −w/2,  Ŷ = −(ηF′² + F w/2),  ĥ = η²F′²/2,   8ρ⁶D = 8 D(P̂, Q̂, Ẑ, Ŷ, ĥ).

(The script checks these five against sympy's jets of F(x/ρ), and 8D(P̂, …) against the typed N: `conical_parts`.) Degree in (F′, F″) with F fixed: P̂ = P̂₁ + P̂₂ with P̂₁ = −F F″ of degree 1 and P̂₂ = −F′² of degree 2; Q̂, Ẑ are of degree 1; Ŷ = Ŷ₁ + Ŷ₂ with Ŷ₁ = −F w/2.
The lowest degree is 3: D̂₃ = −P̂₁(Ŷ₁² − Ẑ²) + 2Q̂Ŷ₁Ẑ, where Ŷ₁² − Ẑ² = (F² − 1)w²/4 and 2Q̂Ŷ₁Ẑ = −ηF w³/4. So 8D̂₃ = 2F w²((F² − 1)F″ − ηw). **The dropped term is −2F F″w²: it is −P̂₁ times the Ẑ² in Ŷ² − Z², which the README formula left out.**
With F, F′, F″ all small: Ŷ = −(ηF′² + F w/2) of degree 2, ĥ of degree 2; the lowest degree is 4, D̂₄ = ĥ(Ẑ² − Q̂²) + P̂Ẑ² + 2Q̂ŶẐ with P̂ = −(F′² + FF″), which gives 8D̂₄ = −w²[(η² + 1)(η² + 2)F′² + 2F(ηw + F″)].

**C9b (ii): limits at η → ∞, tested directly.** Substituting F = F_∞ + c η^{−k} (k ∈ {½, 1, 3/2, 2, 5/2, 3, 4, 5, 6}) or F_∞ + c e^{−κη} (κ ∈ {½, 1, 2}) with F_∞ ∈ {0, ⅓, 1, −1, 2} into N and reading the lowest order, `conical_limit` finds that
no exponential is allowed and that no power is allowed except **k = 2 with c = −2F_∞ for F_∞ = ⅓, 1, −1, 2** (and none for F_∞ = 0). So leading order does **not** exclude F = F_∞(1 − 2/η²). The criterion "no case is open" fails, as README's C9b (ii) words it, and the failure gates the exit code.
Observation, not a criterion (printed by `conical_limit`): in F = F_∞ + a₂/η² + a₃/η³ + … the order η^{−8} of N gives a₂ = 0 or a₂ = −2F_∞, the next order gives a₃ = 0, and the third order after the leading one gives
F_∞⁴(F_∞² + 7) = 0 with no unknown left, so for F_∞ ≠ 0 there is no pure power series on the a₂ = −2F_∞ branch. Logarithmic terms are not tested. *Added by the
coordinator:* `I80-mac/ind80.py` (K7) finds a₃ = 0 at η⁻⁹ and the condition −512F_∞⁴(F_∞² + 7) = 0 at η⁻¹⁰, which is the
second order after the leading η⁻⁸, not the third.
The planted detections work (`conical_planted`): the same scan applied to F″ + 2F′/η finds the limit F_∞ + c/η at k = 1 for every F_∞ and not at k = 2, and applied to F″ + F′ it finds F_∞ + c e^{−η} at κ = 1 and not at κ = 2.

**Verdict for the conical reduction: "no limit at η → ±∞" is not established; localization in x for conical shifts is not decided.** Evidence: leading order allows F = F_∞(1 − 2/η²) at k = 2 (**checked symbolically**); a pure power series in 1/η with a₂ = −2F_∞ fails at third order (an observation, **checked symbolically**);
log terms are not tested. Anything I had written about localization in x that rested on the failed (b), (c) is open for this reduction. For F_∞ = 0, for every other k tested and for all the exponential approaches tested, the exclusion at leading order holds.

**Update from #159** ([`conical_limit.md`](conical_limit.md)): limits at η → ±∞ exist. For every F_∞ ≠ 0, N = 0 has the
two-parameter family F = F_∞ − 2F_∞/η² + (C + F_∞(F_∞² + 7) ln η)/η⁴ + … (formal series **checked symbolically**,
integration to 1e-10 **numerical evidence**, existence **conditional** on a standard theorem). The pure power series
fails at η⁻¹⁰ because k = 4 is a resonance, and the logarithm absorbs it. For F_∞ = 0, and for other leading
balances, there is none, **conditional** on a leading-monomial form of F − F_∞. Localization in x for conical shifts is
therefore not excluded; whether a global bubble exists is not decided. *Update from #168* ([`conical_attraction.md`](conical_attraction.md)):
existence near the family and the exclusion of every other approach, including F_∞ = 0 and oscillation, are now
**proved**, with certified constants.

## 7. The far field of a localized shift

For u = u_∞ + ε w(x, ρ) the ε⁰, ε¹, ε² coefficients of D vanish and the ε³ coefficient is −u_∞ Z²[(1 − u_∞²) w_xx + Δ⊥ w] with Z = ½ w_xρ (`farfield`, exact). So for u_∞ ≠ 0 a shift that tends to u_∞ with all its derivatives has, at leading order at large |x|, either w_xρ = 0 or
(1 − u_∞²) w_xx + Δ⊥ w = 0, which is elliptic for |u_∞| < 1 and hyperbolic for |u_∞| > 1 in (x, ρ). This is a necessary condition on the leading order only; it does not classify anything and is not used in the verdicts above.

**Correction from #159** ([`conical_limit.md`](conical_limit.md) §7): the ε³ coefficient alone gives this condition
only if the ε⁴ terms are of higher order in 1/x. That fails for w ~ ρ²/x², the decay of the conical family: there the
ε³ and ε⁴ coefficients are 128U⁴ρ²/x⁸ and −128U⁴ρ²/x⁸ at leading order, and they cancel. So the condition is not
necessary for such shifts (`I159-mac/ind159.py`, M4).

## 8. Planted truth

- #65's b(ρ) families are found with D = 0 (§3, c = 0 with ε = ±1, 6 points each) and fail closure (6 exact jets).
- The bubble at constant speed passes on no open set (`bubble`): u = w(r), 16 r³D is a cubic in c² = cos²θ whose c⁰ and c⁶ coefficients equal the printed ones of [`two_scalars.md`](two_scalars.md), its coefficients and (2w − r w′, (w′ + r w″)²) generate the same ideal after inverting r w′
  (Gröbner bases, both inclusions), and 2w = r w′ with w′ + r w″ = 0 on an interval force w′ = 0 (differentiating the first gives w′ = r w″). This reuses #65's route with a re-derived D (not imported).
- The planar closure criterion and the closure code are tested on three planted passes and three planted failures (§2).
- The series search of §5 was to be tested on the additive equation; it was withdrawn (README A5).

## 9. Checks and controls

15 checks, 13 pass and 2 fail (`bessel_pert`, `conical_limit`: criteria that failed on the real data and stay as written, §3 and §6); 19 controls, each fires on its target and on no other check
(`det_classify_output.txt`; the 2 baseline failures are not counted as extra failures of a control; each control's own stdout is in the raw folder, `final_mut/`). The script exits 1 because of the two failed checks. Per check, in seconds:
`d_block` 2, `travel` 0, `planar` 2, `planar_closure` 13, `planar_traveling` 0, `bessel` and `bessel_pert` 14, `additive` 0, `multiplicative` 3, `conical` 1, `conical_parts` 0, `conical_limit` 5, `conical_planted` 0, `farfield` 0, `bubble` 2; 43 s in all.

| check | what it tests | control `--mutate=` | the control's change |
|---|---|---|---|
| `d_block` | block of R₃ for a generic b, P = a_x, Y = (b_xb_ρ − a_ρ)/2, det = D | `d-sign` | Y built as (b_xb_ρ + a_ρ)/2 |
| `travel` | traveling = static with u = B + v; Galilean covariance | `travel-offset` | u = B, without + v |
| `planar` | R₃ = diag(P, −P, 0), D = 0; 8 + 8 exact points | `planar-leak`, `planar-point` | + ρ²/10 in the symbolic check; + ρ²/10 in the 'form' points |
| `planar_closure` | null-coordinate criterion, dω = Θ′/2, planted jet, bubble | `closure-deform`, `theta-factor`, `planted-bubble` | + pq²/10 in the planted passes; Θ′/4; the bubble replaced by the jet with Θ′(0) = 0 |
| `planar_traveling` | null coordinates, R_ab = (uu′)′ g_ab, d_pd_q Λ(q − p) = −Λ″ | `pg-null` | A′ = 1/(1 + u) |
| `bessel` | factorization of D, D = 0 at 36 points, R = aa − bb, no ω at 6 jets | `bessel-freq`, `bessel-closed` | J₀(cρ) in place of J₀(√2 cρ); the closed pair (dρ, dt + dx) |
| `bessel_pert` | perturbed |D| > 1e-3 (**fails**) | none (its target already fails at baseline) | |
| `additive` | D = hE, E = Σ f_ig_i, Wronskian, unit ideal of the elimination step, X″ = 0 reduction | `add-drop`, `add-quad` | X² for (X + R)²; the y² coefficient equation dropped |
| `multiplicative` | the polynomial E, the log form, the subcases, equal numerators of the compatibility conditions | `mult-drop-v` | u = X R without + v |
| `conical` | 8ρ⁶D = N | `conical-term` | coefficient 2 → 3 of F³F″³η² |
| `conical_parts` | hand entries, corrected parts of degree 3 and 4 | `conical-hand` | the original wrong degree-3 part |
| `conical_limit` | no limit at η → ∞ by a power or an exponential (**fails**) | `conical-limit` **dropped** (target already fails at baseline), README A3 | |
| `conical_planted` | the same scan finds the limits of two planted equations | `conical-plant` | F″ + F′/η in place of F″ + 2F′/η |
| `farfield` | the ε³ coefficient of D | `far-sign` | + u_∞ |
| `bubble` | cubic in c², printed c⁰ and c⁶, ideal equality, ODE step | `bubble-w2`, `bubble-ode` | uu_xx dropped from P; rw″ = +w′ |

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/det_classify.py --no-controls      # the 15 checks
    python3 notes/warp/det_classify.py --mutate=NAME       # one control
    python3 notes/warp/det_classify.py                     # the checks, then every control (about 15 min on a Raspberry Pi 5)

Times are on a Raspberry Pi 5 (aarch64, 4 cores shared with other jobs), one process at nice 10; memory was sampled, not measured by `/usr/bin/time`: 73 MB resident at the start of the final run (systemd `MemoryPeak` 64 MB then).
The withdrawn controls and checks: `mult_series` and `mult-series-add` (README C8 (e), A5) and `conical-limit` (A3). Failed and stopped runs are kept in the raw folder (`dev_*`, README amendments A1 to A5).

## 10. What the code proves, and what is assumed

The code proves, as exact identities (sympy, rational or symbolic jets): the block of R₃ for a generic b against the coordinate Ricci tensor (`d_block`); the reductions of §Setup (`travel`); R₃ = diag(P, −P, 0) for planar shifts (`planar`);
the null-coordinate form of the planar metric and R_ab = (uu′)′ g_ab (`planar_traveling`); the obstruction identity dω_{tξ} = ½Θ′ (`planar_closure`); the factorization of D for u = c x + g(ρ) and R = a a − b b on the ODE (`bessel`); the identity D = h·E and the expansion of E (`additive`),
and the unit ideal of the elimination step by Gröbner bases (`additive`); the polynomial E for u = X R + v and the equal numerators of the compatibility conditions (`multiplicative`); the ODE N of the conical shifts, the hand entries, the corrected parts of degree 3 and 4, and the leading-order scan (`conical`, `conical_parts`, `conical_limit`);
the far-field coefficient (`farfield`); the bubble's coefficients and the ideal equality (`bubble`). Exact outcomes at the sample points: `planar` (8 + 8 points, exact rationals), `planar_closure` (iii) (5 points), `bessel` (d) (6 jets, exact algebraic numbers).

The 36 values of |D| in `bessel` (b) are **numerical evidence** at 130 working digits (not exact), as are the perturbed values.

By hand: the lemma (#65, reviewed in #79); the closure criterion for planar shifts (uniqueness of the null factors k, m); Theorem T; Theorem A, steps (1) to (3) and (5) (linear algebra; step (4) uses the Gröbner result for the polynomial identity and a hand computation of the slope);
the Bessel reduction of §3; the hand derivation of the conical parts (§6), which the script checks.

Imported and trusted: `ricci` (used in `d_block`, `bessel` (c), `planar_closure`, `planar_traveling`), `ricci_at_point` and `two_scalar_form` (`planar`), `closure` and `dform` (`planar_closure`, `bessel` (d)), all from the code of #57 and #65, reviewed in #59 and #79. The checks that call them are not independent of #65.
`closure` at a point (`bessel` (d)) evaluates ω at the point, so only the failure of the two 2-form equations is a sound conclusion there; the script marks any other outcome as not decided, and none occurred.
Assumed: two real, minimally coupled scalars with a Lorentzian kinetic form and any potential (so λ = 0 by the lemma), G_ab = 8π T_ab, axisymmetric b of class C² (C³ in Theorem T), sympy 1.14.0.

Not tested at all: closure on the conical family; the hypothesis |u| ≠ 1 of Theorem T beyond the cases named in §2; any non-polynomial asymptotic form at η → ±∞ in §6.

## 11. Lemmas for Lean (#135)

Statements only; `ℝ` throughout; D as in §Setup.

1. **Planar** (`ring`). `D(P, 0, 0, 0, 0) = 0`, i.e. with h = Q = Z = Y = 0 in the formula for D.
2. **Static and traveling** (`ring`). With jets (b, b_x, b_ρ, b_t, b_tx, b_xx, b_xρ, b_tρ, b_ρρ) of B(x − vt, ρ) the value D equals D_static(B + v).
3. **Galilean covariance** (`ring`). D(b + c; b_t + c b_x, b_tx + c b_xx, b_tρ + c b_xρ) = D(b; b_t, …).
4. **Additive identity** (`ring`). With u = X + R, D = h·((h + X′²)² − ((X + R)X″)² − Q²).
5. **Theorem A step (4)**. For real β₃, β₄, γ₀, ε₀, β₅, δ₀ with (y² − ε₀ − β₅y)(γ₀ + β₃y) = (δ₀ + β₄y)² for infinitely many y, β₃ = 0 and γ₀ = β₄² (`Polynomial.eq_zero_of_infinite_isRoot`). Then γ₀ = m² and β₄ = ±m/2 force m = 0.
6. **Theorem A step (1)–(3)**: linear algebra (relations among 1, R, R² when R is nonconstant, `LinearIndependent`).
7. **Bessel reduction**. If ψ″ + ψ′/ρ + 2c²ψ = 0 and ψ > 0 then g = −ε ln ψ satisfies g″ + g′/ρ = ε(g′² + 2c²), and D(c x + g) = 0.
8. **Planar closure** (assumed for the differential geometry): ω construction as in #79 §10 (statement 12), and dω_{tξ} = ½Θ′ for Θ = (1 − u²)Φ′/Φ.
9. **Theorem T**: a continuity argument as in the proof of §2, with the limits of e^{αs}.
10. **Conical parts** (`ring`). With the five entries of §6, the parts of degree 3 and 4 of 8D̂ are the stated ones.

## 12. Not covered, and what `two_scalars.md`'s "Not covered" bullet should say

Not classified here: the general static or traveling PDE D = 0; time-dependent shifts beyond planar b(t, x) and the traveling reduction; multiplicative separable shifts (§5); the conical reduction's behavior at η → ±∞ (§6); planar shifts with closure and equal limits at every t, with K of both signs (§2, #161); closure on the conical family. Two scalars of one sign, three scalars and non-minimal coupling are as in [`two_scalars.md`](two_scalars.md).

**For #158 (one-sign pairs), by hand, not in the script.** For a planar shift R₃ − λη₃ = diag(P + λ, −P − λ, −λ), which has rank 1 at λ = −P = K (the φ direction supplies the second square), so planar shifts pass the algebraic stage of
[`same_sign_scalars.md`](same_sign_scalars.md) with λ = K and gradients along y and z. M = R − λg = K (dy² + dz²) is invariant along the (t, x) kernel only if K is constant, which is then dS₂ × ℝ² or AdS₂ × ℝ² (#78's P1, P2) in a
Painlevé–Gullstrand chart; by the formulas of §2 such a shift has Θ′ = −2K ≠ 0, so the (1, 1) pair of this note fails closure on it, consistently. Whether #158 needs more than this is for that issue.
*Added by the coordinator:* the sign is M = R − λg = −K(dy² + dz²) (`I80-mac/ind80.py`, K9), as #78's Lemma S needs for
(+,+) at λ = K < 0.

Suggested replacement for the bullet "Shifts b(t, x, ρ) other than the bubble and b(ρ)": "Shifts b(t, x, ρ) other than the bubble and b(ρ): [`det_classify.md`](det_classify.md) (#80). D = 0 on an open set for every planar shift b(t, x) (*proved*), which can be localized in x but is constant in ρ.
Traveling planar shifts with |u| ≠ 1 everywhere fail closure unless flat (*proved*). Within the additive-separable reduction every shift that depends on ρ and solves D = 0 is u = c x + d + g(ρ) with g from a Bessel equation, which is not localized in x and fails closure at the jets tested (*proved*, *checked symbolically*).
Multiplicative separable and conical shifts, and the general static PDE, are not classified."

## Coordinator's review

Raw `I80-mac` (negative-mass-raw). Mac, nice 10, at most 2 processes.

| step | result |
|---|---|
| scope | the three files the issue names, and raw `I80-rpi5` |
| rerun `det_classify.py` | exit 1 in 6 min 22 s, as ruled: `bessel_pert` and `conical_limit` fail, 13 of 15 checks pass, 19 of 19 controls fire on their targets only. The output equals rpi5's apart from timings |
| hand proofs | Theorem T (planar, traveling) and Theorem A (additive) were read and their steps checked; no gap found |
| re-solve | the series on the branch a₂ = −2F_∞ (§6), order by order with the coordinator's own code (K7) |
| independent check `ind80.py` | run 1 exits 0 in 3 s, and all 10 mutations fire on target |

`ind80.py` imports nothing from the repository. It builds D from #65's formula, which `I79-mac/ind79.py` checked against
the coordinator's own Ricci routine, and it uses that routine for the planar metric. It shares with this note the
reductions, the printed formulas and the claims it tests.

| id | what | run 1 |
|---|---|---|
| K0 | planted: Painlevé–Gullstrand Ricci-flat | pass |
| K1 | P = a_x, Y = ½(b_x b_ρ − a_ρ); planar shifts have D ≡ 0 | pass |
| K2 | traveling = static with u = B + v; Galilean covariance | pass |
| K3 | planar R = (uu′)′g₂; the null coordinates; R_pq = −Φ/2 | pass |
| K4 | the factorization of D for u = cx + g(ρ), the Bessel reduction, \|D\| ≤ 4.4e-60 at 8 own points | pass |
| K5 | D = hE and E = Σ f_i g_i for additive shifts | pass |
| K6 | the conical entries, the corrected parts of degree 3 and 4, and the k = 2 balance c = −2F_∞ | pass |
| K7 | the series condition F_∞⁴(F_∞² + 7) = 0, at η⁻¹⁰ | pass |
| K8 | the far-field ε³ coefficient | pass |
| K9 | R − Kg = −K(dy² + dz²) for a planar static shift | pass |

So the reductions, the Bessel family, the additive identity, the conical parts and the far field are confirmed by code
that shares nothing with `det_classify.py`.

**Discrepancies.**
- Two criteria failed on the real data and stay failed: `bessel_pert` (C6 (b): perturbed \|D\| = 9.06e-6, under its
  1e-3 threshold) and `conical_limit` (C9b (ii)). The script exits 1, and the conical verdict is "not decided".
- C9 (b), (c), the hand formulas for the conical parts, failed first and were replaced by C9b after the failure.
- The control `conical-limit` was dropped because its target fails at baseline. The series search `mult_series` was
  withdrawn before any run of record.
- In §6 the series condition is at the second order after the leading one, not the third (coordinator sentence there).
- In §12 the remark for #158 has the wrong sign of M (coordinator sentence there).
- In §6 "the only symmetry" has no symmetry computation behind it (coordinator sentence there).
- Many prose lines here are much longer than AGENTS.md's 120 characters.

**Recommended tasks.**
- 1 is #159.
- 2 is #160.
- 3 is #161.
- 4 is #162.
- 5, the Lean statements of §11, is row 24 of #135.

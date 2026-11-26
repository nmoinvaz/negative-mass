# Attraction to the conical family's line of equilibria

November 2026, issue [#168](https://github.com/nmoinvaz/negative-mass/issues/168). Script
[`conical_attraction.py`](conical_attraction.py), its stdout
[`conical_attraction_output.txt`](conical_attraction_output.txt). It continues [`conical_limit.md`](conical_limit.md)
(#159), §2, §4 and §6. Status words as in AGENTS.md: **proved** (every step written here), **certified**
(computer-assisted: an exact rational computation whose code is named), **checked symbolically**, **numerical
evidence**, **conditional**. A computer-assisted proof needs an independent review before it is relied on; that review
has not been done (§10).

## Answer

Route: a written proof; no theorem is cited and none is used. With ν = 1/η, σ = ν², G = F, y = η³F′/F, ζ = y − 4 and u =
σ(1 + G²), the system of #159 §4 reduces to one scalar equation, y′ = y(4 − y)/2 − σy² + ê/2 (′ = d/d ln η), where ê is
O(u) on the "huge root" of the cubic N = 0 in F″. A Gronwall bootstrap with explicit constants proves the attraction.

| question | answer | status |
|---|---|---|
| attraction (Theorem A) | if \|ζ₀\| ≤ 1/4, u₀ ≤ 1/1000 and \|ê₀\| ≤ 1 at η₀, the solution exists on [η₀, ∞), \|ζ\| ≤ 0.451 (η₀/η)^(7/4), \|ζ\| ≤ 1.14 (\|ζ₀\| + 50ū ln(η/η₀)) (η₀/η)², F → F_∞ ≠ 0 with \|ln(F/F_∞)\| ≤ (9/4)η⁻² | **proved**; constants **certified** |
| rate | (1 + ln η) η⁻² is attained: the resonance is a Jordan block, which is the logarithm | **proved** (Theorem B) |
| open set of data | every (F(η₀), F′(η₀)) with \|ζ₀\| ≤ 1/4, u₀ ≤ 1/1000 and the huge root gives a limit solution: a 2-parameter family | **proved** |
| formal series = asymptotic expansion (Theorem B) | for every m ≥ 1, \|F − F^(2m)\| ≤ M_m η^{−2m−2} (1 + ln η)^m with C = C(solution) | **proved** |
| oscillation (Theorem C) | F′ = 0 forces F″ = 0; F is eventually constant or strictly monotone; a branch analysis of the cubic leaves only the family | **proved**; constants **certified** |
| the branch a₂ = 0, other leading balances | none, for every real C² solution on a half line with a limit; assumption H of #159 is not needed | **proved** |
| F_∞ = 0 | F ≡ 0 near infinity: no other solution has the limit 0 | **proved** |
| η → −∞ | mirror image (N contains η only as η²) | **checked symbolically** (#159) |
| cited theorem (Aulbach, LNM 1058) | not used; the book's data confirmed, theorem number and hypotheses *(unverified)* (§6) | *(unverified)* |
| numerical comparison | the inequalities of Theorem A hold along 15 integrations; the direct integration of N = 0 agrees to 9e-12 | **numerical evidence** |

Not claimed: that every pair (F_∞, C) is realized by an actual solution; anything about solutions without a limit, or
defined only on a bounded interval; the global question on all of η ∈ ℝ (#167); closure on the family (#172). The facts
N = 8ρ⁶D and D's entries are those of #65, #80 and #159 (assumed, §10).

## 1. Variables and identities

η = x/ρ, t = ln η, ν = 1/η, σ = ν². For F ∈ C²: G = F, p = ηF′, q = η²F″ + ηF′ (θ = η d/dη = d/dt, so p = θG, q = θp), r
= η³F′ = p/σ, s = η⁴F″ + η³F′ = q/σ = r′ − 2r. N = A₀ + ν²A₁ + ν⁴A₂ of `conical_limit.md` §1. For G ≠ 0: y = r/G, ζ = y
− 4, a = Gσ, b = G²σ, u = σ + b = σ(1 + G²); for p ≠ 0: ω = 2Gσ/p = 2/y, w = q/p². Put

    κ = 1 − b + σ,  μ = 1 + 3σ − b + 2σ² − 6bσ,  c = 1 + 2σ,  M = μ + (b − σ)ω,
    P(w) = −2aκ w³ − M w² + 4ac w + c² − ω².

**Identities** (M1, **checked symbolically**; the script retypes D, builds N from it and compares). (a) N = A₀ + ν²A₁ +
ν⁴A₂. (b) N = c₃q³ + c₂q² + c₁q + c₀ with c₃ = −2Gσκ, c₂ = −p(2G(G² − 1)σ² + pμ), c₁ = 4Gσc p⁴, c₀ = p⁴(c²p² − 4G²σ²);
every coefficient but c₃ is divisible by p. (c) N = p⁶ P(q/p²) with ω = 2Gσ/p. (d) N(−G, −p, −q) = N, so F ↦ −F is a
symmetry and G > 0 is no loss; y, ζ, u and ω are invariant, a and w change sign. (e) With r = σ⁻¹p, s = σ⁻¹q: R := N/σ⁴
= σ²r⁶ P(s/(σr²)).

**The scalar equation** (M2). Put ê := (r² + 2Gs)/G² = y²(1 + 2aw) (since s = σr²w). From G′ = p = σr and r′ = 2r + s, s
= (G²ê − r²)/(2G):

    (ln|G|)′ = σ y,      y′ = y(4 − y)/2 − σ y² + ê/2,      σ′ = −2σ,      i.e.  ζ′ = −2ζ − ζ²/2 + ê/2 − σ(4 + ζ)².     (1)

(checked symbolically: (a) of M2). Conversely, a solution (G, y) of (1) with ê = y²(1 + 2aw) and P(w) = 0 gives F = G ∈
C² with p = σGy and q = p′ = σG(ê − y²)/2 = p²w, hence N = p⁶P(w) = 0.

**Equilibria, eigenvalues, the resonance** (M2, M3, **checked symbolically**). With ρ(ê; y, σ, g₂) := R/G⁴ at r = Gy, s
= G(ê − y²)/2 and g₂ = G² (a polynomial): ∂ρ/∂ê = −y⁴/4 at (ê, σ) = (0, 0), so ê = ê*(y, σ, g₂) is analytic near σ = 0
with ê*(y, 0, g₂) = 0, and at y = 4: ê* = −8(g₂ + 3)σ + O(σ²). Then f(ζ, σ, g₂) := −ζ²/2 + ê*/2 − σ(4 + ζ)² has f(0, 0,
g₂) = 0, f_ζ(0, 0, g₂) = 0 and **f_σ(0, 0, g₂) = −4(g₂ + 7)**. The system (G′ = σr, r′ = 2r + S, σ′ = −2σ) of #159 §4
has the Jacobian [[0, 0, 4G], [8, −2, S_σ], [0, 0, −2]] at (G, 4G, 0), S_σ = −4G(G² + 3): eigenvalues 0, −2, −2, and the
2×2 minor of J + 2I is −8G(G² + 7) ≠ 0 for every real G ≠ 0, so the eigenvalue −2 has a Jordan block. The same factor G² +
7 is the logarithm coefficient b = F_∞(F_∞² + 7) of #159 §3.

## 2. The cubic and its roots

**Lemma 1 (shape).** Let G > 0, σ > 0, u = σ + b ≤ 1/500, a = Gσ (so a² = σb ≤ u²/4, a ≤ 1/1000), and |ω| ≤ 100. Then
(i) κ ∈ [0.998, 1.002], μ ≥ 0.99, c ∈ [1, 1.01], M ≥ 0.79; (ii) P has exactly one real root ≤ −1/(4a), simple; (iii)
every other real root lies in (−2, 2), and there is such a root only if |ω| ≤ 1.05; (iv) for |ω| ≤ 1/2 there are exactly
two: w₋ ∈ (−2, −1/2) and w₊ ∈ (1/2, 2).

*Proof.* The numbers of (i) and the inequalities (A)–(F) below are the content of check M4: exact rational interval
arithmetic over σ, b ∈ [0, 1/500], a ∈ [0, 1/1000], ω ∈ [−100, 100] (**certified**): (A) 4M − 3κ > 0 (≥ 0.1859); (B)
32a²P(−1/(4a)) = κ − 2M + 32a²(c² − c − ω²) < 0 (≤ −0.594); (C) P(2) < 0, P(−2) < 0 (≤ −2.1) and P′(2) < 0; (D) c² +
16aκ + 8ac ≤ 1.0321 ≤ 1.05²; (E) for |ω| ≤ 1/2: P(0), P(1/2), P(−1/2) > 0 (each ≥ 0.496). Symbolically P′(w) = −6aκw² −
2Mw + 4ac is a concave parabola with discriminant 4M² + 96a²κc > 0, P′(0) = 4ac > 0 and P′(−1/(4a)) = (4M − 3κ)/(8a) +
4ac. By (A) and (C), P′ > 0 at −1/(4a) and at 0 and P′(2) < 0, so its roots satisfy w_c < −1/(4a) < 0 < w_M < 2: P
decreases on (−∞, w_c], increases on [w_c, w_M], decreases on [w_M, ∞), and P(−∞) = +∞. By (B), P(−1/(4a)) < 0 and P is
increasing at that point, so P(w_c) < 0: exactly one root in (−∞, w_c), simple (P′ < 0 there), none in [w_c, −1/(4a)].
On [−1/(4a), −2] ⊂ [w_c, w_M], P is increasing and P(−2) < 0 by (C): no root. On [2, ∞), P decreases from P(2) < 0: no
root. This is (ii) and the first half of (iii); P(±2) < 0 makes the number of roots in (−2, 2) even. If w ∈ [−2, 2] is a
root, c² − ω² = 2aκw³ + Mw² − 4acw ≥ −16aκ − 8ac (M > 0), so ω² ≤ c² + 16aκ + 8ac ≤ 1.05² by (D): the second half of
(iii). For (iv): P(0) > 0 and P(−2) < 0 < P(−1/2), with P increasing on [−2, w_M]: one root in (−2, −1/2); P(1/2) > 0 >
P(2) and P decreasing on [w_M, 2], P(w_M) ≥ P(1/2) if w_M ≥ 1/2: one root in (1/2, 2); the count is even, so these are
all. ∎

The root ≤ −1/(4a) is the **huge root** w = −x, x ≥ 1/(4a). Physically q = p²w ≈ −p²/(2Gσ): the branch s ≈ −r²/(2G) of #159.

**Lemma 2 (size of ê).** For the huge root, y ≥ 1/50 and the hypotheses of Lemma 1:
|ê| ≤ u(αy² + βy), α = 51/25, β = 181/50, where ê = y²(1 + 2aw) = y²(1 − 2ax).

*Proof.* P(−x) = 0 gives, exactly (M5, **checked symbolically**), κ(1 − 2ax) = (κ − M) − 4ac/x − (ω² − c²)/x², and κ − μ
= −(2σ + 2σ² − 6bσ). With 1/x ≤ 4a, ω = 2/y, 16a² ≤ 4u² and |ω² − c²| ≤ ω² + c²:
|ê| ≤ y²[(2σ + 2σ² + 6bσ) + |b − σ|·2/y + 16a²(c + c²) + 64a²/y²]/κ. Over the simplex σ = uf, b = u(1 − f), u ≤ 1/500 (interval
arithmetic, **certified**): the y² coefficient is ≤ 2.0362u ≤ αu, the y coefficient ≤ 2.0041u ≤ (201/100)u, the constant
≤ 0.03206u ≤ (321/10000)u, and 321/10000 ≤ (β − 201/100)·(1/50) = 0.0322 turns the constant into a multiple of y for y ≥
1/50. ∎

At y ∈ [7/2, 9/2]: |ê| ≤ c_e u with c_e = α(9/2)² + β(9/2) = 57.6 ≤ 58.

## 3. Theorem A: attraction

**Theorem A.** Let F ∈ C²([η₀, η₁)) solve N = 0 with G₀ = F(η₀) > 0 (G₀ < 0: apply to −F). Put σ₀ = η₀⁻², u₀ = σ₀(1 +
G₀²), y₀ = η₀³F′(η₀)/G₀, ζ₀ = y₀ − 4, ê₀ = (r₀² + 2G₀s₀)/G₀². Assume |ζ₀| ≤ 1/4, u₀ ≤ 1/1000, |ê₀| ≤ 1. With τ =
ln(η/η₀), ū = e^{9σ₀/2}u₀ ≤ 1.0046·10⁻³, the following hold on [η₀, η₁): (a) q is the huge root and y ∈ [7/2, 9/2]; (b) |ζ|
≤ (|ζ₀| + 200ū)e^{−7τ/4} ≤ 0.451 e^{−7τ/4}; (c) |ζ| ≤ 1.14(|ζ₀| + 50ū τ)e^{−2τ}; (d) |ê| ≤ 58u, u ≤ ūe^{−2τ}; (e) G →
G_∞ as η → ∞ (if η₁ = ∞), G₀e^{7σ₀/4} ≤ G_∞ ≤ G₀e^{9σ₀/4} and |ln(G/G_∞)| ≤ (9/4)σ. Moreover every such (G₀, y₀) at
every η₀ has a unique solution on [η₀, ∞) with these properties, and F ∈ C².

*Proof.* **Step 1 (the root stays huge).** At η₀: w₀ = (ê₀ − y₀²)/(2a₀y₀²) ≤ (1/y₀² − 1)/(2a₀) ≤ (1/14.06 − 1)/(2a₀) ≤
−0.464/a₀ ≤ −1/(4a₀) (by hand, y₀ ≥ 15/4). On any interval where G > 0, y ≥ 1/50 and u ≤ 1/500, w = q/p² is continuous
and, by identity (c) of §1, a root of P. By Lemma 1 the real roots lie in (−∞, −1/(4a)] ∪ (−2, 2), separated by the gap
[−1/(4a), −2], which contains none. So w stays ≤ −1/(4a): it is the huge root, and ê = y²(1 + 2aw) obeys Lemma 2. **Step
2 (bootstrap).** Let T be the supremum of τ with |ζ| ≤ 1/2 on [0, τ] (T > 0 by continuity). On [0, T): y ∈ [7/2, 9/2],
(ln G)′ = σy ≤ (9/2)σ gives G² ≤ G₀²e^{9σ₀/2}, so u(τ) = σ(1 + G²) ≤ e^{9σ₀/2}σ₀(1 + G₀²)e^{−2τ} = ūe^{−2τ} ≤ 1/500, and
Lemma 2 gives |ê| ≤ c_e u, c_e = 58. By (1), with φ = e^{2τ}ζ: φ′ = −(ζ/2)φ + e^{2τ}(ê/2 − σy²), and
|ê/2 − σy²| ≤ (c_e/2 + (9/2)²)u ≤ c_f u, c_f = 50 (29 + 20.25 = 49.25), so |e^{2τ}(ê/2 − σy²)| ≤ c_fū.
**Stage 1.** |ζ/2| ≤ 1/4, so |φ(τ)| ≤ |ζ₀| + ∫₀^τ [(1/4)|φ| + c_fū]dτ′; Gronwall: |φ(τ)| ≤ (|ζ₀| + 4c_fū)e^{τ/4} (the
constant 4 is 1/(1/4)). With |ζ₀| ≤ 1/4 and ū ≤ 1.0046·10⁻³: |ζ| ≤ B e^{−7τ/4}, B = 1/4 + 200ū ≤ 0.451 < 1/2. So the
bound
|ζ| ≤ 1/2 is never reached and T is the end of the interval. This is (a), (b), (d).
**Stage 2.** ∫₀^∞ |ζ|/2 dτ ≤ B/(2·7/4) =: H ≤ 0.129, e^H ≤ 1.14. From |φ|′ ≤ (|ζ|/2)|φ| + c_fū and the integrating
factor,
|φ(τ)| ≤ e^H(|ζ₀| + c_fūτ). This is (c).
**Step 3 (G).** (ln G)′ = σy ∈ [(7/2)σ, (9/2)σ], σ = σ₀e^{−2τ}: ln G increases, ∫₀^∞σy dτ ∈ [(7/4)σ₀, (9/4)σ₀], so G_∞
exists with the stated bounds and |ln(G/G_∞)| = ∫_τ^∞σy ≤ (9/4)σ(τ). This is (e). **Step 4 (existence).** On the open
set {G > 0, y > 1/50, u < 1/500} the huge root x is a simple root of P (P′ < 0 there), so by the implicit function
theorem it, and ê = y²(1 − 2ax), are analytic in (G, y, σ); (1) is an analytic system and has a unique local solution.
By §1 each solution of (1) is a C² solution of N = 0, and by Steps 1–3 it stays in the compact set {|ζ| ≤ 0.451, |ln G −
ln G₀| ≤ (9/4)σ₀, τ ≤ T} for every T, so it continues to [η₀, ∞). Uniqueness: Step 1 shows any C² solution with these
data is on the huge root, which is the solution of (1). ∎

The constants of the chain (c_e ≤ 58, c_f = 50, ū ≤ 1.0046·10⁻³, B ≤ 0.451, H ≤ 0.129, e^H ≤ 1.14, all in exact rational
arithmetic with a rigorous upper bound of the exponential) are check M6 (**certified**); the planted chain on ζ′ = −2ζ +
ū e^{−2τ} (exact (ζ₀ + ūτ)e^{−2τ}) is P2: the bound is ≥ the exact solution and equals it up to the factor e^H. Check M7
(**numerical evidence**): fifteen integrations (G₀ ∈ {10⁻², 0.3, 1, 3, 30}, ζ₀ ∈ {−1/4, 0, 1/4}, u₀ = 1/1000, to τ = 20,
h = 1/400) satisfy (b)–(e) at every step (smallest bound/actual ratios 1.81, 1.14, 1.89, 1.09 for (b), (c), the ê bound,
(e)); a direct integration of N = 0 (Newton in q on the θ-form, which uses neither the cubic P nor ê) agrees with (1) to 9.0e-12
at t = 6 (h = 1/800; first worded at h = 1/200 it failed, §11).

*Sharpness.* If G_∞² + 7 ≠ 0 (always) the Jordan block makes ζ ≈ −4(G_∞² + 7)σ ln η: the factor τ in (c) cannot be
dropped. The neighborhood is far smaller than the true basin: the numerics of #159 and M16 converge from y₀ = 0.5 to 30.

## 4. Theorem B: the series is the expansion

Let (G, ζ) be a solution of Theorem A, t = ln η, σ = e^{−2t}, G_∞ its limit. For large t, ê = ê*(4 + ζ, σ, G²) (the
analytic branch of §1: both solve ρ = 0, the root near 0 is unique, and |ê| ≤ 58u → 0), so ζ′ = −2ζ + f(ζ, σ, G²) with f
analytic.

**Theorem B.** There is a number C such that, for every m ≥ 1, |F(η) − F^(2m)(η; G_∞, C)| ≤ M_m η^{−2m−2}(1 + ln η)^m
for η ≥ η₀, where F^(2m) is the series of `conical_limit.md` §3, F_∞ + Σ_{n ≤ 2m} a_{n,j} ν^n (ln η)^j with a_{2,0} =
−2F_∞, a_{4,0} = C, a_{4,1} = F_∞(F_∞² + 7), and the other coefficients determined (table in #159 §3).

*Proof.* **Induction.** Claim H_n (n ≥ 0): ζ = Σ_{k=1}^n σ^kπ_k(t) + O(σ^{n+1}(1 + t)^{n+1}) with polynomials π_k, deg
π_k ≤ k, and the same for ζ′. H₀ is Theorem A(c). H_n ⟹ H_{n+1}: (i) ln(G/G_∞) = −∫_t^∞σ(4 + ζ) dτ, and ∫_t^∞
e^{−2kτ}τ^j dτ = e^{−2kt}× a polynomial of degree j, so ln(G/G_∞) = Σ_{k=1}^{n+1}σ^kλ_k(t) + O(σ^{n+2}(1 + t)^{n+1}), λ₁
= −2, deg λ_k ≤ k − 1; hence G² = G_∞²(1 + Σσ^kγ_k(t)) + the same remainder. (ii) f(ζ, σ, G²) is analytic at (0, 0,
G_∞²); substituting (i) and H_n into its Taylor series gives f = Σ_{k=1}^{n+1}σ^kφ_k(t) + O(σ^{n+2}(1 + t)^{n+2}), φ₁ =
f_σ(0, 0, G_∞²) = −4(G_∞² + 7), deg φ_k ≤ k. (iii) ζ = e^{−2t}[c + ∫^t e^{2τ}f dτ]: the term k = 1 integrates to σφ₁t;
for k ≥ 2, ∫^t e^{−2(k−1)τ}φ_k dτ = const − e^{−2(k−1)t}ψ_k(t), deg ψ_k = deg φ_k; all the constants, and the one from
the remainder, are absorbed in c, and the tail ∫_t^∞ e^{2τ}O(σ^{n+2}τ^{n+2}) = O(σ^{n+1}t^{n+2}) gives ζ = σ(c + φ₁t) +
Σ_{k≥2}σ^kπ_k + O(σ^{n+2}(1 + t)^{n+2}): this is H_{n+1}, with π₁ = c − 4(G_∞² + 7)t. ζ′ follows from the equation.
**Expansion of F.** From H_n: F = G_∞(1 + Σ_{k=1}^{n+1}σ^kΓ_k(t)) + O(σ^{n+2}(1 + t)^{n+1}), Γ₁ = −2. With n + 1 = m
this is the form of the claim; the polynomial coefficients are what must be identified. **Identification.** p = σG(4 +
ζ) and q = σG(ê* − y²)/2 have expansions of the same kind (analytic composition, derivatives from (1)). Truncating at
ν^{2n+2}, N(G, p, q, ν) has remainder O(σ^{n+5}(1 + t)^{n+1}) since ∂N/∂q = O(σ³) (check M9); N of the truncated series
is a finite sum Σc_{k,j}ν^kT^j, so c_{k,j} = 0 for k < 2n + 10. Hence the coefficients solve the stage equations of the
formal solution up to stage k ≤ 2n + 3. **Uniqueness of the formal solution.** At stage k ≥ 3 the unknowns a_{k,0..J}
enter the coefficient of ν^{k+6}T^j through ∂N/∂q = −128F³ν⁶ + …, ∂N/∂p = −512F³ν⁶ + … (∂N/∂G is O(ν⁸), the nonlinear
terms of higher order): the equation is triangular in j with diagonal I(k) = −128F³ k(k − 4) (check M9, **checked
symbolically**), and the subdiagonal j·(−512F³ + 256F³k) (θ(ν^kT^j) = ν^k(−kT^j + jT^{j−1})). So for k ≠ 4 the stage has
a unique solution; at k = 4 the diagonal vanishes, a_{4,0} is free and a_{4,j+1} is forced by the T^j equation since
512F³ ≠ 0, which gives a_{4,1} = F(F² + 7). The exact solver confirms this at 12 pairs (F_∞, C) through k = 9: full
column rank at every stage but k = 4, null space spanned by a_{4,0}, the table of #159, N = O(ν¹⁶) (check M8, **checked
symbolically**). The stage k = 2 value a_{2,0} = −2F_∞ is Γ₁ = −2 above. So the actual solution's coefficients are the
formal series with C := a_{4,0}. ∎

The link between the series and the induction (check M10, exact): for the series, ζ̃ = η³F̃′/F̃ − 4 has ν⁰ coefficient
0, σt coefficient −4(F_∞² + 7) = f_σ and σ coefficient F_∞² + 15 − 4C/F_∞ (affine in C with slope −4/F_∞ ≠ 0): the free
constant c of the induction is C. Check M11 (**numerical evidence**, mpmath, 30 digits): from G₀ = 1, ζ₀ = 1/5, η₀ = 45
the integration to η = 10⁷ gives F_∞ = 1.00101194680117415 and C = −123.0139833; e_m(η) = |F − F^(2m)|/(ν^{2m+2}(1 + ln
η)^m) at η = 100, 200, 400, 800 lies between 9.0 and 2.5·10⁴ and decreases with η for every m = 1..4 (threshold 10⁵).

## 5. Oscillation and classification

**Lemma 3 (local).** Let F ∈ C² solve N = 0 on an interval J on which F ≠ 0 and κ(F, ν) = 1 − (F² − 1)ν² ≠ 0. If F′
vanishes at one point of J then F is constant on J.

*Proof.* At p = 0: N(G, 0, q, ν) = c₃q³ (every other coefficient has a factor p, §1(b)), c₃ = −2Gσκ ≠ 0, so q = 0: F′ =
F″ = 0. Near a point η₁ ∈ J with F′(η₁) = 0, |p| ≤ 1 and G, ν, c₃ vary in a compact set with |c₃| bounded below. Write
c₂ = p c̃₂, c₁ = p⁴c̃₁, c₀ = p⁴c̃₀ (§1(b), c̃_i polynomials, hence bounded there): |c₂/c₃| ≤ L|p|, |c₁/c₃|^{1/2} ≤
L|p|², |c₀/c₃|^{1/3} ≤ L|p|^{4/3} ≤ L|p|. Every root q of N = 0 obeys |q| ≤ 3Θ, Θ = max(|c₂/c₃|, |c₁/c₃|^{1/2}, |c₀/c₃|^{1/3}):
if |q| > 3Θ then |c₂q²| < |c₃||q|³/3, |c₁q| < |c₃||q|³/9,
|c₀| < |c₃||q|³/27, the sum is < (13/27)|c₃||q|³, contradicting N = 0. So |q| ≤ K|p|, i.e. |p′| ≤ K|p| (p′ = q), and
p(t₁) = 0 gives p ≡ 0 near t₁ (Gronwall). The set {F′ = 0} ∩ J is closed in J, nonempty and open: it is J. ∎

Check M12 (**checked symbolically** for the divisibility and N(G, 0, q, ν) = c₃q³; the bound |q| ≤ 3Θ tested on 30000
roots at 10⁴ random parameter points: 0 violations). Planted P4: for N_osc = ν²(q − p) + G (F″ + F = 0) the divisibility
and N(G, 0, q, ν) = c₃q³ both fail, and F = cos η has p = 0, q = π² ≠ 0 at η = π: the lemma must not apply, and the code
says so.

**Corollary 4.** Let F ∈ C²([η₀, ∞)) solve N = 0 with F → F_∞. (i) F_∞ ≠ 0: either F ≡ F_∞ on some [η₁, ∞), or F′ ≠ 0 on
[η₁, ∞): F is strictly monotone and F − F_∞ has one sign. (ii) F_∞ = 0: either F ≡ 0 on some [η₁, ∞), or F ≠ 0 and F′ ≠
0 on [η₁, ∞).

*Proof.* (i) Choose η₁ with |F − F_∞| < |F_∞|/2 and (F² − 1)ν² < 1/2 on J = [η₁, ∞); apply Lemma 3. (ii) Choose η₁ with |F|
≤ 1 (κ ≥ 1). On a component K of {F ≠ 0} ∩ [η₁, ∞), Lemma 3 applies: F is constant on K or F′ ≠ 0 on K. A constant c ≠ 0
on K contradicts F → 0 if K is unbounded and continuity if K has an endpoint where F = 0, so F is strictly monotone on
K. A bounded component (a, b) with F(a) = F(b) = 0 is impossible; so is an unbounded component (a, ∞) with F(a) = 0.
What remains: F ≠ 0 on all of [η₁, ∞), or one first component [η₁, b) and F ≡ 0 on [b, ∞). ∎

So F − F_∞ cannot change sign infinitely often: sin η/η² and every other sign-changing approach is excluded. What is
left is the size of r = η³F′: it could still vary within one sign. Theorem C closes that.

**Theorem C.** Let F ∈ C²([η₀, ∞)) solve N = 0 with F → F_∞. Then either F is constant on some [η₁, ∞), or F_∞ ≠ 0 and
there is η₂ such that the data of F at η₂ satisfy the hypotheses of Theorem A (for F_∞ < 0 apply it to −F); in
particular F = F_∞ − 2F_∞η⁻² + (C + F_∞(F_∞² + 7) ln η)η⁻⁴ + O(η⁻⁶ ln²η) and η³F′/F → 4. If F_∞ = 0, F ≡ 0 near
infinity.

*Proof.* By symmetry F_∞ ≥ 0. If F is not eventually constant, Corollary 4 gives G > 0 and p ≠ 0 of one sign on [η₁, ∞)
(for F_∞ = 0 and F ≢ 0 near infinity, G > 0 after a sign change of F). G → F_∞, so G ≤ g₂ < ∞ and u = σ(1 + G²) ≤ 1/1000
for t ≥ t₂. **Regime B (|y| ≤ 1/50).** Check M13 (**certified**): for |y| ≤ 1/50, a ≤ 1/1000, u ≤ 1/500, every real root
satisfies |s| ≤ 0.037|r|: by the bound of Lemma 3, 3Θ ≤ 3·max(0.01208, 3·10⁻⁵, 0.0035)|p| ≤ 0.0363|p|. So (ln|y|)′ = 2 +
s/r − σy ≥ 1.9 − σ|y| > 0 at |y| = 1/50 and |y| increases exponentially while ≤ 1/50: in finite time |y| ≥ 1/50 for good
(the set {|y| ≥ 1/50} is forward invariant). **Regime A (|y| ≥ 1/50, |ω| ≤ 100).** Lemma 1 holds; w(t) = q/p² is
continuous and a root of P, on the huge side (≤ −1/(4a)) or in (−2, 2) for all later t (the gap contains no root).
*Small side.* If w ∈ (−2, 2) at some time then |ω| ≤ 1.05 (Lemma 1(iii)). Since (ln|ω|)′ = p/G − 2 − pw (M14) and, for |ω|
≥ 1/2,
|w| ≤ 2: |p|(1/G + 2) ≤ (2Gσ/|ω|)(1/G + 2) ≤ 4σ + 8a ≤ 0.016 (**certified**), (ln|ω|)′ ≤ −1.98: after at most ln 2.1 < 0.75 we have |ω| ≤ 1/2
forever, and by Lemma 1(iv) w stays in (−2, −1/2) or in (1/2, 2) (continuity). Then p′ = p²w gives (1/p)′ = −w: p > 0, w
< −1/2: 1/p ≤ 1/p₄ + 2(t − t₄), p ≥ 1/(1/p₄ + 2(t − t₄)), ∫p = ∞, so G → ∞: contradiction. p > 0, w > 1/2: 1/p decreases
at rate ≥ 1/2 and reaches 0: p unbounded in finite time: contradiction. p < 0, w < −1/2: 1/p < 0 increases at rate ≥ 1/2
and reaches 0: contradiction. p < 0, w > 1/2: |p| ≥ 1/(1/|p₄| + 2(t − t₄)), ∫|p| = ∞, G → −∞: contradiction with G > 0.
*Huge side.* p < 0: (1/p)′ = −w ≥ 1/(4a) ≥ 1/(4g₂σ(t₃)) > 0, so 1/p < 0 reaches 0 in finite time: contradiction. p > 0:
y > 0, and with
|ê| ≤ u(αy² + βy) (Lemma 2) and u, σ ≤ 1/500 (check M15, **certified**): y′ ≥ 1.49y on [1/50, 1]; y′ ≥ 1/10 on [1, 15/4]; y′ ≤ −0.104y on [17/4, ∞)
(the bracket is increasing in y). So y reaches [15/4, 17/4] in finite time, where |ζ| ≤ 1/4; u ≤ 1/1000 and |ê| ≤ 58u ≤
1 hold there: the hypotheses of Theorem A, at that η₂. Theorem A and B give the expansion. **Leaving regime B.** At the
exit |ω| = 100 > 1.05 no small root is real, so w is on the huge side, for either sign of p, and the huge-side argument
applies. **F_∞ = 0.** G > 0 is decreasing (p < 0): contradiction as above, unless F ≡ 0 near infinity; p > 0 would make
G increase from a positive value. ∎

Checks: M12 to M15 above and M16 (**numerical evidence**, RK4 with the root of P tracked by continuity, η₀ = 100): the
huge root with p > 0 converges to y = 4 (|y − 4| ≤ 3e-12 at t = 15, from y₀ ∈ {0.5, 8, 30}, G₀ ∈ {0.5, 2}); with p < 0
it blows up at t = 1.10 (y₀ = −0.5) and 0.35 (y₀ = −4); the two small roots keep w in (−2, −1/2) and (1/2, 2) on [3, 15]
for y₀ = ±3, ±6.

## 6. The cited theorem

#159 §4 cited the theorem on normally attracting manifolds of equilibria (Aulbach, *Continuous and Discrete Dynamics near Manifolds of
Equilibria*, Lecture Notes in Mathematics 1058, Springer 1984). It is **not used** here. What was looked up (raw
`citation_attempts.md`): the bibliographic data (author, title, series volume 1058, 1984, DOI 10.1007/BFb0071569, ISBN
978-3-540-13329-2, XII + 144 pages) through a search-result summary of the Springer record. The Springer page itself
redirects to a cookie handshake (not followed) and two catalog pages returned an access-denied page or HTTP 403. So the
theorem number and its exact hypotheses stay *(unverified)*. The hypotheses that such a theorem would need are checked
here by `hyp_line` (check M3, **checked symbolically**): the system is analytic; the points (G, 4G, 0) are equilibria;
zero is a simple eigenvalue with the line as eigenvector; the other eigenvalues are −2, −2 (a Jordan block). The planted
cases P1 (the linear field with eigenvalues 0, −2, −2 and a block: attracting, with |x₂,₃| ≤ (1 + t)e^{−2t}|x₀| exact;
diag(0, −2, 2): not attracting; diag(0, 0, −2): kernel larger than the line) show that `hyp_line` can say no.

## 7. What `conical_limit.md` should now say

For the Answer table: *existence for all data near the family* → "**proved** (conical_attraction.md Theorem A): data
with |ζ₀| ≤ 1/4, u₀ ≤ 1/1000 on the huge root"; *the branch a₂ = 0, other leading balances* → "**proved**, for every
real C² solution on a half line with a limit (Theorem C); assumption H is not needed"; *F_∞ = 0* → "**proved**: F ≡ 0
near infinity". For §4: replace the last paragraph (theorem cited, bootstrap not written) by Theorems A and B. For §6:
(d) is decided (no oscillating approach, Corollary 4 and Theorem C) and (e) is Theorem A. For §2: H is replaced by
Theorem C for solutions with a limit; the dominant-balance argument stays as an explanation of why only k = 2 occurs.
Still **conditional** in `conical_limit.md`: nothing from this issue; the independent review (§10) and the assumed D
remain.

## 8. Lemmas for Lean (#135, not written)

All over ℝ, with the constants named. (L1) Shape: for 0 < a ≤ 1/1000, σ, b ≥ 0, σ + b ≤ 1/500, |ω| ≤ 100, P has exactly
one root ≤ −1/(4a), simple, and the other real roots lie in (−2, 2), only if |ω| ≤ 1.05; two of them, in (−2, −1/2) and
(1/2, 2), if |ω| ≤ 1/2. (L2) |ê| ≤ u(αy² + βy) on the huge root for y ≥ 1/50. (L3) Gronwall chain: φ′ ≤ hφ + k, h, k ≥
0, k constant gives φ(τ) ≤ e^{∫h}(φ(0) + kτ) and with h ≤ 1/4: φ ≤ (φ(0) + 4k)e^{τ/4}. (L4) Local lemma: N(G, 0, q, ν) =
c₃q³ and |q| ≤ 3Θ. (L5) Corollary 4 (a continuous function on [η₁, ∞) with |p′| ≤ K|p| on a clopen set). (L6) The four
sign cases of (1/p)′ = −w. (L7) y′ bounds on [1/50, 1], [1, 15/4], [17/4, ∞). (L8) The induction step H_n ⟹ H_{n+1}
(variation of constants). (L9) Uniqueness of the formal solution (triangular stage equations, I(k) = −128F³k(k − 4)).

## 9. Planted truth, checks and controls

The planted cases run first through the same functions: P1 `planted_hyp` (`hyp_line`), P2 `planted_gronwall` (the
bootstrap chain), P3 `planted_shape` (the root classifier and hypotheses on −2a(w − h)(w² − 1): huge root −50 and ±1;
with h = +50 the hypotheses fail), P4 `planted_local`, P5 `planted_series` (the stage solver on q + p = 0: F_∞ + cν, no
log; on q + 4p = 8ν⁴: obstruction without a log, a_{4,1} = −2 with it).

Run of record 01 (`run_final.sh 01`, one process at nice 10 under systemd): exit 0, 8 min 4 s, 98 MB peak
(`/usr/bin/time -v`). 21 of 21 checks pass (70 s; M11 takes 55 s of it), then 27 of 27 controls fire, each on its named
check and on no other, with no traceback.

| id | check | result of the run of record |
|---|---|---|
| P1–P5 | planted | all pass: Jordan field attracting with |x₂,₃| ≤ (1 + t)e^{−2t}|x₀|, diag(0, −2, 2) not attracting, diag(0, 0, −2) kernel larger than the line; chain ≥ exact and = e^H·exact; huge root −50 and ±1; N_osc not covered; F_∞ + cν and a_{4,1} = −2 |
| M1 | `n_identity` | (a)–(e) exact |
| M2, M3 | `reduction`, `hyp_line_real` | f_σ = −4(g₂ + 7); eigenvalues {0, −2, −2}; minor −8G(G² + 7); attracting with a Jordan block at G = −3, −1/2, 1, 2, 7 |
| M4 | `shape` | κ ∈ [0.998, 1.002], M ≥ 0.798, 4M − 3κ ≥ 0.186, 32a²P(−1/(4a)) ≤ −0.594, ω² ≤ 1.0321 |
| M5 | `ehat` | y² coefficient ≤ 2.0362u ≤ 51/25 u, y coefficient ≤ 2.0041u, constant ≤ 0.03206u |
| M6 | `theorem_a` | c_e = 57.6 ≤ 58, ū = 1.00451·10⁻³, B = 0.45090, H = 0.12883, e^H ≤ 1.13750 |
| M7 | `numeric_a` | 15 cases to τ = 20: smallest bound/actual 1.81, 1.14, 1.89, 1.09; direct vs reduced 8.98e-12; step study 6.37e-12 |
| M8 | `series` | 12 pairs (F_∞, C): stages k = 3..9 full rank but k = 4 (null space a_{4,0}); table; N = O(ν¹⁶) |
| M9, M10 | `stage_symbol`, `expansion_link` | I(k) = −128F³k(k − 4); σt coefficient −4(F² + 7), σ coefficient F² + 15 − 4C/F |
| M11 | `numeric_b` | F_∞ = 1.00101194680117415, C = −123.013983339; e_m ≤ 2.47·10⁴ (m = 4, η = 100), ≤ 10⁵ |
| M12–M15 | `local_lemma`, `regime_b`, `small_branch`, `huge_branch` | 0 violations in 30000 roots; 3Θ ≤ 0.0363 |p|; (ln|ω|)′ ≤ −1.984; y′ ≥ 1.4923y, ≥ 0.1062, ≤ −0.1042y |
| M16 | `classify_numeric` | huge, p > 0: |y − 4| ≤ 2.6e-12; huge, p < 0: blow-up at t = 1.10 and 0.35; small roots in their bands |

The controls fire on their targets: `d-sign`, `P-drop` (M1); `ehat-sign`, `f-sigma` (M2); `eig-plus` (M2 and M3);
`shape-3`, `shape-om` (M4); `ehat-alpha` (M5); `A-u0`, `A-cf` (M6); `A-kappa`, `plant-gron` (P2); `num-rate` (M7);
`series-a2`, `series-8` (M8); `stage-I` (M9); `link-slope` (M10); `numB-coef` (M11); `local-c2` (M12); `regimeb-y`
(M13); `small-sigma` (M14); `huge-thr` (M15); `class-sign` (M16); `plant-hyp`, `plant-shape`, `plant-local`,
`plant-series` (P1, P3, P4, P5). The weakest margins are the mutations that only just cross a bound: `small-sigma` (true
−1.984, claimed −1.99), `huge-thr` (1.4923 against 1.5), `ehat-alpha` (2.0362 against 2). `A-kappa` as coded changes
only the planted chain, so its target is `planted_gronwall`, not `theorem_a`. `num-rate` makes the rate-bound ratio 6.8e-5
against ≥ 1. `numB-coef` gives e₄ up to 4·10⁹.

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/conical_attraction.py --no-controls     # the 21 checks, about 70 s on a Raspberry Pi 5
    python3 notes/warp/conical_attraction.py --mutate=NAME      # one control (add --skip-num to skip M7, M11, M16)
    python3 notes/warp/conical_attraction.py                    # the checks, then every control (8 min)

## 10. What the code proves, and what is assumed

The code proves, as exact rational or symbolic statements: the identities of §1 (M1, M2); the Jacobian, its Jordan block
and f_σ = −4(g₂ + 7) (M2, M3); the inequalities of Lemma 1 and Lemma 2 (M4, M5, interval arithmetic over exact
Fractions: Python integers, no floating point, **certified**); the rational chain of Theorem A (M6, with a rigorous
upper bound of the exponential); the series, the table of #159, N = O(ν¹⁶) at 12 pairs, the rank of every stage and the
indicial polynomial (M8, M9); the link (M10); the divisibility and the root bound of Lemma 3 (M12, the bound tested
numerically); the constants of Regime B, the small branch and the huge branch (M13 to M15). It computes **numerical
evidence** only in M7, M11 and M16.

By hand, written in this note, not coded: the shape argument from the inequalities (concave P′, monotone pieces); the
Gronwall steps of Theorem A; the induction of Theorem B and the identification with the formal series (the orders of the
remainders, the triangularity for all k); Lemma 3's Gronwall step; Corollary 4's topology; the case analysis of Theorem
C; the one-line arithmetic w₀ ≤ −0.464/a₀ of Theorem A Step 1.

Assumed: D = h(h² − P² − Q² + Y² + Z²) − P(Y² − Z²) + 2QYZ with the entries of `det_classify.md` Setup for a static
shift (#65's lemma; the script retypes it and compares N with A₀ + ν²A₁ + ν⁴A₂, shared source with #80 and #159); sympy 1.14.0,
mpmath 1.3.0, numpy. The independent review of this computer-assisted proof is not done: it should re-derive the
interval constants of M4, M5, M6, M13 to M15 with other code.

## 11. Record of the work

Raw folder `I168-rpi5` of negative-mass-raw: README with route, constants, criteria and planted cases written before the
first run of record (6fcc19c); amendments A1 to A8 (before the run of record); development logs `dev_01` to `dev_13`;
the exploration scripts in `explore/`; the run of record and every control's log; `citation_attempts.md`.

**M7 (A1).** M7 as first worded (h = 1/200, direct integration of N = 0 against the reduced system to 1e-9) **failed**
in the development run `dev_04`: y differed by 2.5e-9. The cause is RK4 truncation of the direct integrator, which is
fourth order: the deviation is 4.1e-8, 2.5e-9, 1.6e-10, 9.0e-12, 2.6e-12 for h = 1/100 … 1/1600 (`dev_05`, factor 16 per
halving). The coordinator's condition that h = 1/400 and h = 1/800 differ by less than 1e-10 **failed as worded**
(`dev_09`: direct y differs by 1.48e-10; reduced system 8e-15); the threshold was the coordinator's error. The run of
record uses h = 1/800 (tolerance 1e-9 unchanged, deviation 9.0e-12) and a step study 1/800 versus 1/1600 below 1e-10 for
both integrators (6.4e-12); the inequality trajectories use h = 1/400. These settings were chosen after `dev_05` and
`dev_09` and known to pass. M7 is numerical evidence. Other amendments: A2 (a README sentence with a wrong inequality;
the code checks the right one), A3 to A5 (float slack 1e-12, wording), A6 (the stop rule of M16 ended small-root runs
early; criteria unchanged), A7 (P′(2) < 0 added to M4), A8 (a secant crash in a control run).

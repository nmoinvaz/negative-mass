# The two Frank–Loss gaps of #68, closed

Issue #74. Script [fl_gaps.py](fl_gaps.py), output [fl_gaps_output.txt](fl_gaps_output.txt), run from the repository root as
`python notes/warp/fl_gaps.py` (about 14 min in all, one thread, Normal priority, 206 MB peak resident). It exits 1 unless every check
passes and every `--mutate=` control fires. `--no-controls` runs the checks only, `--only=a,b` runs sections, and `--mutate=NAME`
runs one control. The author is the `negative-mass-win2` agent, and nothing here has been independently reviewed.

Read for this note, and not edited, are [review_fl.md](review_fl.md) (#68), [frank_loss.md](frank_loss.md) (#63), Frank–Loss
[FL22], arXiv:2201.03610v1 (Section 1 with Remark 4, Section 2 with Proposition 7, and Section 3), and the regularity theorem
and proof of their earlier paper [FL21], arXiv:2012.13646 (Section 2). Both were read from the arXiv source on 2026-11-13, and
title, authors and dates were checked against the arXiv API. The source of [FL21] is its version 2 (2022-01-10).

## Answer

Both gaps close. Each is a short argument by hand, and the script checks the operator-level content of each step.

| gap | what it says | status |
|---|---|---|
| G1 | [FL22] Remark 4's inequality alone gives ψ ∈ Ḣ¹, and the coefficient is used through its adjoint in exactly one step | *proved* (by hand), given [FL21]'s regularity proof as read |
| G2 | [FL22] Proposition 7 holds for every ψ ∈ Ḣ¹, not only where ∇f ∈ L^d | *proved* (by hand), the continuity test is *numerical evidence* |

G2 is on #63's path and G1 is not. #63 applies Proposition 7 to the ψ and χ of its §2, which lie in Ḣ¹ ∩ L³ and are not
smooth in general. #63 gets Ḣ¹ by construction, so G1 only settles the remark as printed (#68's discrepancy 3), as #68 said.

The script runs 17 checks and 13 controls. All checks pass and every control fires.

## G1. Remark 4 needs no hypothesis beyond the inequality

[FL22] Remark 4 assumes 0 ≢ ψ ∈ L^{2d/(d−1)} and |γ·∇ψ| ≤ |A||ψ| with A ∈ L^d. Its Section 3 proof uses ψ ∈ Ḣ¹, and the remark
does not say where that comes from. In Theorem 1, Ḣ¹ comes from [FL21]'s regularity theorem and the equation. This section
writes out the same argument for the inequality, at d = 3, which is the case #63 uses.

### Reading R

"|γ·∇ψ| ≤ |A||ψ|" is taken to say that the distribution g := σ·(−i∇)ψ is a function and that |g| ≤ |A||ψ| a.e. Nothing weaker
makes the inequality meaningful.

### Lemma G1

Let A ∈ L³(ℝ³, ℝ³), 0 ≢ ψ ∈ L³(ℝ³, ℂ²) and g := σ·(−i∇)ψ a function with |g| ≤ |A||ψ|. Then ψ ∈ L^r for every 3/2 < r < ∞,
ψ ∈ Ḣ¹, and ‖∇ψ‖₂ = ‖g‖₂.

### Proof, step 0, the matrix

Put M := gψ†/|ψ|² where ψ ≠ 0 and M := 0 where ψ = 0. Where ψ = 0, |g| ≤ |A||ψ| forces g = 0, so Mψ = g a.e. Since M has rank one,
its operator norm is |g|/|ψ| ≤ |A|, so |M| ∈ L³. Thus σ·(−i∇)ψ = Mψ with |M| ≤ |A|. [FL21] has the same equation with M = σ·A,
which is Hermitian. This M is not, unless g is parallel to ψ (`m_matrix`, median |M − M†|/|M| = 1.4 for unrelated g and ψ).

### Proof, step 1, the integral equation

Let Γ(x) := (i/4π)σ·x/|x|³, so that σ·(−i∇)Γ = δ and |Γ(x)| = 1/(4π|x|²). By Hölder, |g| ≤ |A||ψ| ∈ L^{3/2}. By
Hardy–Littlewood–Sobolev (HLS, the kernel |x|⁻² in ℝ³ maps L^{3/2} to L³), Γ∗g ∈ L³. Also σ·(−i∇)(Γ∗g) = g, by approximating g in
L^{3/2} by C_c^∞. So σ·(−i∇)(ψ − Γ∗g) = 0, hence −Δ(ψ − Γ∗g) = 0 since (σ·(−i∇))² = −Δ. A harmonic distribution in L³ is smooth by
Weyl's lemma, its averages over large balls tend to 0, and the mean value property makes it vanish. Hence ψ = Γ∗g = Γ∗(Mψ).

The script checks this equation at a point, for a smooth ψ whose g is known. It does so for Frank–Loss's extremal (g = Λψ, so M is
Hermitian) and for a generic ψ whose M is visibly non-Hermitian, with relative errors 4.0×10⁻¹⁵ and 5.7×10⁻¹⁵
(`integral_equation`).

### Proof, step 2, the a priori bound

This is [FL21]'s duality argument with the adjoint written. Fix 3/2 < r < ∞ and set 1/s = 1/r + 1/3, so 1 < s < 3. For Λ > 0 let

    S_Λ := sup{ |∫⟨φ, ψ⟩| : ‖φ‖_{r′} ≤ 1, ‖φ‖_{3/2} ≤ Λ },

which is at most Λ‖ψ‖₃, hence finite. Since |M| ≤ |A| ∈ L³, dominated convergence lets us write M = F + G with
F := M 1{|x| < R, |M| < R}, bounded with compact support, and ‖|G|‖₃ ≤ η, with η fixed below. Then ψ = Γ∗(Fψ) + Γ∗(Gψ), and for
admissible φ

    ∫⟨φ, ψ⟩ = ∫⟨φ, Γ∗(Fψ)⟩ + ∫⟨G†(Γ∗φ), ψ⟩.

The second term is where the coefficient enters. Fubini (justified by bilinear HLS, 2/3 + 2/3 + 2/3 = 2) and Γ(x − y)† = Γ(y − x)
make Γ∗ self-adjoint, and ⟨v, Gw⟩ = ⟨G†v, w⟩ is the definition of the adjoint. [FL21] writes G(Γ∗φ), which is the same thing
because their G is Hermitian.

With χ := G†(Γ∗φ) we have |χ| ≤ |G||Γ∗φ| since |G†| = |G|. HLS and Hölder give ‖χ‖_{3/2} ≤ C₁η‖φ‖_{3/2} ≤ C₁ηΛ and
‖χ‖_{r′} ≤ C₁η‖φ‖_{r′} ≤ C₁η. So χ/(C₁η) is admissible, and |∫⟨χ, ψ⟩| ≤ C₁η S_Λ. The first term is at most
‖Γ∗(Fψ)‖_r ≤ C‖Fψ‖_s ≤ C‖F‖_{3s/(3−s)}‖ψ‖₃ by HLS and Hölder. Hence S_Λ ≤ C₂‖F‖‖ψ‖₃ + C₁ηS_Λ. Choose η = 1/(2C₁), so that
S_Λ ≤ 2C₂‖F‖‖ψ‖₃ for every Λ. Letting Λ → ∞ and using duality, ψ ∈ L^r. The exponent bookkeeping is exact in rational
arithmetic over 6754 admissible pairs (p, r) (`exponents`).

### Proof, step 3, Ḣ¹

Take r = 6. By Hölder with 1/2 = 1/3 + 1/6, g ∈ L², so g ∈ L^{3/2} ∩ L². Take g_n ∈ C_c^∞ with g_n → g in both norms and put
ψ_n := Γ∗g_n. Each ψ_n is smooth and decays like |x|⁻², and ‖∇ψ_n‖₂ = ‖σ·(−i∇)ψ_n‖₂ = ‖g_n‖₂, because |σ·k v| = |k||v|
(integrate by parts, and the cross terms cancel by the Clifford relations). So (ψ_n) is Cauchy in Ḣ¹, with limit ψ̃ ∈ Ḣ¹ ⊂ L⁶.
By HLS, ψ_n → ψ in L³, and ψ_n → ψ̃ in L⁶, so ψ = ψ̃. Then ψ ∈ Ḣ¹ in [FL22]'s sense, since ∇ψ ∈ L² and {|ψ| > t} has finite
measure for ψ ∈ L⁶, and ‖∇ψ‖₂ = ‖g‖₂. ∎

The script checks the identity ∫|σ·∇ψ|² = ∫|∇ψ|² on non-smooth and smooth fields in Ḣ¹ (`h1_step`, relative difference
below 3×10⁻¹⁶).

### Where Hermiticity is used

Every use, in the order of the proof. "Coefficient" means M, or σ·A in [FL21].

| step | uses | what for |
|---|---|---|
| 0, the matrix, and the split M = F + G | nothing | only Mψ = g and \|M\| ≤ \|A\| |
| 1, Γ and Weyl | σ_j† = σ_j, Clifford relations | \|Γ(x)\| = 1/(4π\|x\|²), (σ·(−i∇))² = −Δ |
| 2, Γ∗ self-adjoint | σ_j† = σ_j | Γ(x)† = Γ(−x) |
| 2, the duality identity (eq. `decomp` of [FL21]) | **the adjoint G†** | χ := G†(Γ∗φ), the only use of the coefficient |
| 2, the bounds on χ | \|G†\| = \|G\| | norms only, so \|G\| ≤ \|M\| ≤ \|A\| is all that enters |
| 3, Ḣ¹ | σ·k Hermitian | \|σ·k v\| = \|k\|\|v\| |

Only the fourth row involves the coefficient, and it needs the adjoint and nothing more. The script evaluates each
operator-level step for a general complex 2×2 matrix field. It then evaluates the formula that holds for Hermitian G on the same
general field (`audit`).

| step | general M, adjoint used | general M, Hermitian-case formula |
|---|---|---|
| duality identity | passes (6.4×10⁻¹⁶) | **fails** (1.2) |
| pointwise bound on χ | passes | passes |
| split M = F + G | passes | passes |
| kernel norm, Γ(x)† = Γ(−x), symbol norm | passes | passes |

`audit` fails if a step fails for the general field, or if the set of steps that need the adjoint is not exactly the duality
identity. The Hermiticity of the Pauli matrices is separate, and `kernel` checks its three uses (0, 6.7×10⁻¹⁶, 1.8×10⁻¹⁵).

### What the checks do not cover

HLS, Weyl's lemma, Fubini and the density statements are standard and were recalled and not re-read. The decomposition of step 2 and
the duality argument are [FL21]'s, read from the source, and the only change is the adjoint. The case d > 3 is not examined.
[FL22] cite [FL21] for it, and the argument is the same with 3 replaced by d.

### Status of G1

*Proved* (by hand), given [FL21]'s regularity proof as read, and not independently reviewed. The operator-level steps are
*numerical evidence*.

## G2. Proposition 7 holds on all of Ḣ¹

### The gap

Step 3 of the proof of [FL22] Proposition 7 first proves, for a fixed f that is bounded, compactly supported and has ∇f ∈ L^d, the
identity

    ∫ f ⟨∂_kψ, γ_kγ_j∂_jψ⟩ + ∫ f ⟨∂_jψ, γ_jγ_k∂_kψ⟩ = −∫ ∂_kf ⟨ψ, γ_kγ_j∂_jψ⟩ − ∫ ∂_jf ⟨ψ, γ_jγ_k∂_kψ⟩,    j ≠ k,

by density of C_c^∞ in Ḣ¹ in ψ, with f held fixed. It then takes f = χ|ψ|_ε^{−2/(d−1)} with χ ∈ C_c^∞, "which also implies
∇f ∈ L^d". That f depends on ψ, and ∇f contains ∇|ψ|_ε, which is in L² and in general not in L^d.

A concrete case, T4 below, is ψ = (1 + r²)^{−3/2}(a + Σ_j x̂_j b_j) with a, b_j ∈ ℂ². It is homogeneous of degree 0 at the origin,
so |∇ψ| ~ 1/r, which is in L² and not in L³ near 0. For d = 3, f = |ψ|_ε^{−1} and ε = 0.5, the integral of |∇f|³ over the shell
10^{−(i+1)} < r < 10^{−i} is 0.58528 and 0.58122 for i = 1, 2 and 0.58118 for every decade i ≥ 3. So ∫_{s<|x|<1}|∇f|³ ≈ 0.2524 ln(1/s)
diverges, while ‖∇ψ‖₂² = 51.9856 is finite (`gradf_Ld`). Its control evaluates the smooth ψ_δ with δ = 0.05 instead, where the decades
give 0.113, 0.0006 and 0, and the check fails.

So Step 3 as printed does not apply to this f. Proposition 7 is still true for every ψ ∈ Ḣ¹, by the argument that follows.

### The repair

For ψ ∈ Ḣ¹(ℝ^d, ℂ^N), ε > 0, χ ∈ C_c^∞ and s := |ψ|_ε, Step 3 states

    ∫ s^{−2/(d−1)}|∇ψ|²χ = (2(d−1)/(d−2)²) ∫ |∇s^{(d−2)/(d−1)}|²χ + ∫ s^{−2/(d−1)}|γ·∇ψ|²χ
        − (2/(d−1)) ∫ s^{−1−2/(d−1)} Re⟨γ·(∇s)ψ, γ·∇ψ⟩χ + Σ_{j≠k} ∫ s^{−2/(d−1)} Re⟨γ_jψ, γ_k∂_kψ⟩ ∂_jχ.        (S3)

#### Lemma G2

(S3) holds for every ψ ∈ Ḣ¹.

#### Proof

First, for ψ ∈ C_c^∞ the function f := χ s^{−2/(d−1)} is smooth, bounded and compactly supported, so ∇f ∈ L^d. The proof of (S3) in
[FL22], which is the preliminary identity with this f followed by the chain rule for ∂_jf, applies unchanged.

Second, both sides of (S3) are continuous on Ḣ¹. This is shown below.

Third, C_c^∞ is dense in Ḣ¹ ([FL22] say so, by cut-off and mollification), so (S3) extends from C_c^∞ to Ḣ¹. ∎

#### Continuity

Write the terms with ∇s = Re⟨ψ, ∇ψ⟩/s. The first four terms have the form ∫ K(ψ)[∇ψ, ∇ψ] χ, with K continuous on ℂ^N and

    |K(ψ)| ≤ C ε^{−2/(d−1)}.

The coefficients are s^{−2/(d−1)}, β² s^{2β−4}|ψ|² ≤ β² s^{−2/(d−1)} for β = (d−2)/(d−1) (because 2β − 2 = −2/(d−1)), and
s^{−2−2/(d−1)}|ψ|² ≤ s^{−2/(d−1)}. All are at most ε^{−2/(d−1)}, since s ≥ ε. The last term has the form ∫ B(ψ)[∇ψ] ∂_jχ with
|B(ψ)| ≤ s^{−2/(d−1)}|ψ| ≤ s^{(d−3)/(d−1)}, which is bounded for d = 3 and grows sublinearly for d > 3.

Let ψ_n → ψ in Ḣ¹. Then ∇ψ_n → ∇ψ in L² and, by Sobolev, ψ_n → ψ in L^{2d/(d−2)}. Pass to a subsequence that converges a.e.
For the first four terms, with Q_n := ∇ψ_n ⊗ ∇ψ_n → Q in L¹,

    |K(ψ_n)Q_n − K(ψ)Q| ≤ ‖K‖_∞|Q_n − Q| + |K(ψ_n) − K(ψ)||Q|.

The first term tends to 0 in L¹. The second tends to 0 by dominated convergence, since it is bounded by 2‖K‖_∞|Q|. For the last
term, supp ∂_jχ is compact and ψ_n → ψ in L² there. Since |B(ψ_n)|² ≤ C(1 + |ψ_n|²), the generalised dominated convergence
theorem gives B(ψ_n) → B(ψ) in L², and B(ψ_n)∇ψ_n → B(ψ)∇ψ in L¹ by Cauchy–Schwarz. Every subsequence has a further subsequence
with the same limit, so the whole sequence converges.

The rest of the proof is as in [FL22]. Step 4 lets χ_R → 1, using only ∇ψ ∈ L², s ≥ ε, Hardy's inequality and monotone
convergence, all valid on Ḣ¹. Step 5 uses the pointwise identities (comp1) and (comp2), which hold a.e. by the chain rule for
weakly differentiable functions. So Proposition 7 holds for every ψ ∈ Ḣ¹ and every ε > 0.

### Numerical test of the continuity

The identity is evaluated pointwise on analytic spinor fields, with analytic gradients, and integrated in polar coordinates. The
radial integral is Gauss–Legendre in t = ln r, so that the singularity at the origin is resolved, and the angular integral
is a product Gauss rule on the sphere. The fields are ψ = P(r)a + Q(r)Σ_j x_j b_j, in four cases.

| case | ψ | note |
|---|---|---|
| T1 | \|x\|^{−1/4}(1 + \|x\|²)^{−3/2} φ₀, φ₀ = (0.6, 0.8i), the field of the issue | in Ḣ¹ ∩ L³, unbounded at 0 |
| T2 | the same profile times (1 + iσ·x)φ₀, the spinor structure of FL's extremal | unbounded at 0 |
| T3 | the same profile times a + Σ_j x_j b_j with random complex a, b_j | unbounded at 0 |
| T4 | (1 + \|x\|²)^{−3/2}(a + Σ_j x̂_j b_j) with random complex a, b_j | bounded, ∇ψ ∉ L³ |

The smooth approximants ψ_n replace |x|^{−1/4} by (|x|² + δ²)^{−1/8} and x̂ by x/√(|x|² + δ²), and multiply by a smooth cut-off
θ(|x|/R) that is 1 below R and 0 above 2R, so that ψ_n ∈ C_c^∞. The pairs (δ_n, R_n) run from (10⁻², 8) to (10⁻¹², 256) in six steps. At each ψ_n and at ψ,
the script evaluates both sides of Proposition 7 (d = 3, N = 2, with the bracket 4 − 3|ψ|²/|ψ|_ε²), both sides of (S3) with the radial
cut-off χ(r) = θ(r), and both sides of (comp3), which is (S3) with χ = 1. The values of ε are 0.5 and 0.05, except 0.5 and 0.3 for T3,
whose angular structure at 0.05 needs a finer sphere rule than the others.

#### Planted truth first

Three checks come before any real data (`planted_norms`, `planted_extremal`, `quad_converged`). For T1, ∫|∇ψ|² = 7.728546877334
agrees with an mpmath integral to 1.1×10⁻¹⁶. For FL's extremal Ψ = (1 + r²)^{−3/2}(1 + iσ·x)φ₀ at ε = 0 the closed forms are
∫|∇s^{1/2}|² = 3π²/4, ∫|σ·∇Ψ|²/s = 9π²/4 and ∫|∇Ψ|² = 9π²/8, and the quadrature reproduces them to 8.4×10⁻¹⁶. Both sides of
Proposition 7 are then (2/3)(9π²/4) − 2(3π²/4)(4 − 3) = 0, and the code returns 2.4×10⁻³⁰ and −2.2×10⁻¹⁵. For the non-smooth ψ of
all four cases at both values of ε, doubling the radial panels changes every integral by at most 1.1×10⁻¹¹ relative, and a 48 × 64
sphere rule agrees with the 32 × 48 one to 2.0×10⁻¹² for T3 and T4.

#### Result

For each case and ε, the table gives Proposition 7's value at ψ, where the two sides agree, and the Ḣ¹ distance ‖∇(ψ_n − ψ)‖₂ and
the two sides at the first and the last approximant. The full tables, with all six approximants, are in the output file.

| case | ε | value at ψ (LHS = RHS) | ‖ψ₁ − ψ‖₂ | LHS(ψ₁) = RHS(ψ₁) | ‖ψ₆ − ψ‖₂ | LHS(ψ₆) = RHS(ψ₆) |
|---|---|---|---|---|---|---|
| T1 | 0.5 | 1.9769262240 | 3.71×10⁻¹ | 1.9685432671 | 1.17×10⁻³ | 1.9769262237 |
| T1 | 0.05 | 2.7359331181 | 3.71×10⁻¹ | 2.7266779632 | 1.17×10⁻³ | 2.7359331178 |
| T2 | 0.5 | 5.7830572485 | 3.77×10⁻¹ | 5.7825569356 | 1.18×10⁻³ | 5.7830572915 |
| T2 | 0.05 | 2.7877172896 | 3.77×10⁻¹ | 2.8548094195 | 1.18×10⁻³ | 2.7877177219 |
| T3 | 0.5 | 7.5246040727 | 3.19×10⁻¹ | 7.5264557676 | 9.95×10⁻⁴ | 7.5246041209 |
| T3 | 0.3 | 8.7763119484 | 3.19×10⁻¹ | 8.7828909088 | 9.95×10⁻⁴ | 8.7763120290 |
| T4 | 0.5 | 11.7726856566 | 4.98×10⁻¹ | 11.6146648183 | 6.37×10⁻⁶ | 11.7726856567 |
| T4 | 0.05 | 14.9366611641 | 4.98×10⁻¹ | 14.7821221341 | 6.37×10⁻⁶ | 14.9366611646 |

All five checks of this section pass. On every smooth ψ_n, the two sides of Proposition 7, of (S3) and of (comp3) agree, with worst
relative differences 1.1×10⁻¹² and 7.1×10⁻¹² (`prop7_smooth`, `step3_smooth`). On the non-smooth ψ itself they agree to the same
accuracy (`prop7_nonsmooth`, `step3_nonsmooth`). The sides at ψ_n converge to their values at ψ while ‖∇(ψ_n − ψ)‖₂ falls from
0.37 to 1.2×10⁻³ for T1 and from 0.50 to 6.4×10⁻⁶ for T4. The largest ratio of the gap in the sides to the distance is 0.064, so the
sides move no faster than the distance does (`prop7_limit`).

#### The continuity test can fail

`--mutate=weakonly` adds to ψ_n the term n⁻¹e^{−r²}cos(nr)φ₀ with n = 8 to 256, for T1 at ε = 0.5. It tends to 0 in the supremum
norm, in every L^p and weakly in Ḣ¹, and its gradient has norm 0.99 for every n. The check `prop7_smooth` still passes, because the
identity holds for each smooth ψ_n. `prop7_limit` fails, because the distance stays near 0.99 and the sides stay 15 % away from
their values at ψ.

| n | δ_n | ‖∇(ψ_n − ψ)‖₂ | LHS(ψ_n) = RHS(ψ_n) | max gap |
|---|---|---|---|---|
| 8 | 10⁻² | 1.0799 | 2.1710529584 | 0.126 |
| 32 | 10⁻⁶ | 0.9943 | 2.1729990449 | 0.149 |
| 256 | 10⁻¹² | 0.9922 | 2.1727908033 | 0.149 |

The value at ψ is 1.9769262240.

#### What this does and does not show

It tests four one-parameter families at d = 3 and N = 2. The theorem is the continuity lemma above, and the numbers are consistent
with it. They do not test d > 3, N = 4, or families other than these.

### Status of G2

*Proved* (by hand), not independently reviewed. The test above is *numerical evidence*.

## Controls

Each control runs the same code with one input broken, and counts as fired only if its target checks fail.

| `--mutate=` | what it breaks | checks that fail |
|---|---|---|
| `mtranspose` | M = gψᵀ/\|ψ\|² (no conjugate) in place of gψ†/\|ψ\|² | `m_matrix`, `integral_equation` |
| `adjoint_drop` | χ := G(Γ∗φ), G taken Hermitian, in place of G†(Γ∗φ) | `duality`, `audit` |
| `hermpart` | bounds \|G†v\| by the norm of the Hermitian part of G | `norm_bound` |
| `gamma_nonherm` | σ₁ replaced by σ₁ + 0.4iσ₃ | `kernel`, `duality`, `h1_step` |
| `hls_exp` | 1/s = 1/r + 1/2 in place of 1/r + 1/3 | `exponents` |
| `gamma_i` | Γ without its factor i | `integral_equation` |
| `prop7coef` | 4.1 in place of 4 in the bracket of Proposition 7 | `planted_extremal`, `prop7_smooth`, `prop7_nonsmooth` |
| `step3coef` | 1.05 in place of 1 on the cross term of Step 3 | `step3_smooth`, `step3_nonsmooth` |
| `dropt5` | drops the cut-off term Σ_{j≠k} … ∂_jχ of Step 3 | `step3_smooth`, `step3_nonsmooth` |
| `weakonly` | approximants that converge to ψ weakly in Ḣ¹ but not in norm | `prop7_limit` |
| `smoothf` | evaluates ∇f for the smooth ψ_δ in place of the singular ψ | `gradf_Ld` |
| `quadbug` | sphere-rule weights off by 1 % | `planted_norms` |
| `coarse` | 4 Gauss nodes per radial panel in place of 16 | `quad_converged`, `planted_norms`, `planted_extremal` |

All 13 fired. `quadbug` does not fail `planted_extremal`, because that check compares two sides that are both 0 for the extremal, and a
uniform scale error cannot move 0. `planted_norms` catches it.

## Discrepancies and reading notes

1. [FL22]'s sentence "which also implies ∇f ∈ L^d" in Step 3 is false for general ψ ∈ Ḣ¹ (T4, `gradf_Ld`). The proposition is true,
   by Lemma G2. This is a gap in the exposition and not an error in the result.
2. G1 is off #63's path, which agrees with #68. The issue's text says both gaps are where the extension rests on [FL22] as read.
   That holds for G2, and for G1 only as far as the remark as printed.
3. #68 describes [FL21]'s proof as using only |G_ε| and |F_ε|, with Hermiticity entering "only in writing the adjoint". The audit
   confirms it, with one qualification. The adjoint is G_ε† and not |G_ε|, so the proof is not a statement about the norms alone.
   It is a statement about the pair G, G†, which have the same norm.
4. [FL21] was read as version 2. #68 lists it as 2020, and version 1 dates from 2020-12-25. Version 1 was not compared.
5. The 32 × 48 sphere rule that T3 uses is not converged for T3 at ε = 0.05, where it differs from the 48 × 64 rule by up to
   1.5×10⁻⁶ relative (printed as `info:` in the output, and not used in any check). T3 therefore uses ε = 0.3 as its second value.
   This limits the test and not the identity.

## Corrections made to the harness while building it

Each was found by a check or a control, and none loosened a pass criterion in order to pass.

1. `prop7_smooth` failed in the first full run, by 1.8×10⁻⁵ relative and for T4 only. The T4 approximants had no cut-off, while
   the integration stopped at r = 2R, so the sides were off by the boundary term of the truncation (about 3×10⁻⁵ at R = 8). The
   cause was the field and not the identity. T4's approximants now carry the cut-off θ(r/R), and the Ḣ¹ distances are integrated out
   to r = e^14 for every approximant, since the non-smooth ψ has a tail beyond 2R.
2. `quad_converged` first used a tolerance of 10⁻¹¹ that was written before all eight cases were measured. The extended check read
   1.1×10⁻¹¹, from the column ∫|∇f|³ for T1 at ε = 0.05, which no identity check uses. The columns that the identities use change by
   at most 7.1×10⁻¹². The tolerance is now a tenth of the identity tolerance (10⁻¹⁰), which is its purpose, and the printed
   1.1×10⁻¹¹ is within it.
3. Two controls were silent in the first full run, because one target list named a check that does not exist
   (`comp3_nonsmooth`) and because `planted_extremal` compares two sides that are both 0, which a uniform 1 % weight error
   cannot move. Both target lists were corrected, and the checks themselves were not changed.

## Recommended tasks

### Task 1, an independent review of #74's two proofs

What. A reader who has not seen this note reads [FL21] Section 2, [FL22] Remark 4 and Section 2, then Lemma G1 and Lemma G2 here,
and writes a verdict on steps 0 to 3 of Lemma G1, on the three parts of the proof of Lemma G2, and on the continuity estimate.
Why. Both are proofs by hand that only their author has read, and #63's √3 extension relies on G2. Acceptance, each can fail. Each
step is marked confirmed or not, with the line of [FL21] or [FL22] it was checked against. The reviewer also runs
`python notes/warp/fl_gaps.py` and reports exit 0 with every control fired. It fails if any step needs a hypothesis that this note
does not list, or if any control is silent. Depends on #74.

### Task 2, move the status text once Task 1 confirms

What. Cite fl_gaps.md from frank_loss.md (§4 and discrepancy 3), mass_bounds.md item 5 and prior_work.md §4, at the places that
say Remark 4 needs Ḣ¹ "from the equation" or that Step 3 of Proposition 7 is read and not proved. Why. Those sentences are the
two gaps, and they are out of date once the proofs are reviewed. Acceptance, each can fail. `grep -n "does not say how" notes/warp/*.md`
returns nothing outside fl_gaps.md and review_fl.md, and each changed sentence carries the status of Task 1, with "independently
reviewed" only if Task 1 says so. Depends on Task 1 and #68.

## References

- [FL22] R. L. Frank, M. Loss, "A sharp criterion for zero modes of the Dirac equation", arXiv:2201.03610v1 (2022). Remark 4,
  Section 2 (Proposition 7 and its proof) and Section 3 read from the arXiv source, checked against the arXiv API. Journal version
  *(unverified)*.
- [FL21] R. L. Frank, M. Loss, "Which magnetic fields support a zero mode?", arXiv:2012.13646v2. Section 2 (the regularity theorem for
  zero modes and its proof) read from the arXiv source, checked against the arXiv API. Journal version *(unverified)*.
- Hardy–Littlewood–Sobolev, Weyl's lemma, Sobolev's and Hardy's inequalities and the generalised dominated convergence theorem are
  standard and were recalled and not re-read *(unverified)*.

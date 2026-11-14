# Independent review of #63's Frank–Loss extension

Issue #68. Script [review_fl_check.py](review_fl_check.py), output [review_fl_check_output.txt](review_fl_check_output.txt),
run from the repository root as `python notes/warp/review_fl_check.py` (about 6 min, one thread, below-normal priority,
781 MB peak resident; it exits 1 unless every check passes and every `--mutate=` control fires; `--only=a,b` runs
sections). The reviewer is the `negative-mass-win` agent, which had not seen #63's work.

Reviewed: [frank_loss.md](frank_loss.md) §2–§8 at origin/main 2b4a121, [mass_bounds.md](mass_bounds.md) item 5, and
Frank–Loss [FL22], arXiv:2201.03610v1, read from the arXiv source (Sections 1–3, Proposition 7, the Appendix). For one
point their earlier paper [FL21], arXiv:2012.13646, was read (its regularity theorem). Both were fetched on 2026-11-13,
and title, authors and dates were checked against the arXiv API. [FL22] has a single arXiv version.

## Verdict

No error found. Every step of §2–§8 is confirmed. Two gaps lie in the exposition of [FL22], none in #63, and neither
changes a conclusion of #63. My side of G1 is a reading and not a written proof.

| # | Step of frank_loss.md | Verdict | Where |
|---|---|---|---|
| R1 | §2 The reduction D̃ε = 0 ⇔ σ·∇ψ = ½τψ and σ·∇χ = −½τχ, with \|ψ\|² + \|χ\|² = 2\|ε\|², for real and complex τ | confirmed | by hand, `reduction`, `reduction_norms` |
| R2 | §2 These are the pointwise inequalities of [FL22] Remark 4 with \|A\| = ½\|τ\|, N = 2 | confirmed | read, `prop7` |
| R3 | §3 A zero mode in Ḣ¹ lies in L³ (five steps) | confirmed | by hand, `greens_kernel`, `kernel_norm`, `planted_greens` |
| R4 | §4 The inequality as used, Proposition 7's display, the Section 3 chain, S₃, attainment | confirmed | read against the source, `prop7`, `prop7_4x4`, `sobolev_const`, `fl_pair`, `fl_norm`, `fl_local_min` |
| R5 | §4 and discrepancy 3: does Remark 4 need ψ ∈ Ḣ¹ | yes, a gap in [FL22]'s remark (G1), no effect on #63 | read |
| R6 | §5 No kernel for κ_K < √3, also for complex v | confirmed | by hand, `planted_ode` |
| R7 | §6 Compactness, Fredholm index 0, existence and uniqueness | confirmed | by hand |
| R8 | §7 ρ(T₁) ≤ 1/(√3 v_c), E discrete, symmetric, contained in \|v\| ≥ √3 v_c | confirmed | by hand, `block_form`, `planted_ode` |
| R9 | §8 ∫τ\|ε\|² = 0 for a zero mode at real v | confirmed, and true sector by sector (remark 2) | `divxi`, `sector_identity` |
| R10 | §9 The planted case, ρ(T₁) = 1 for τ₁ = 6/(1 + r²) | confirmed, exactly | `planted_ode`, `planted_nocomplex`, `planted_greens` |
| R11 | Answer, the table of v_c and √3 v_c | confirmed | `v_c` |

The chain, each node confirmed.

```
 §2 reduction (R1, R2)  ──►  ψ, χ ∈ Ḣ¹ with |σ·∇ψ| = ½|τ||ψ|
                                        │
 §3 zero modes in L³ (R3) ──────────────┤  ψ, χ ∈ Ḣ¹ ∩ L³
                                        ▼
 §4 [FL22] Remark 4 (R4, R5)  ──►  §5 no kernel for κ_K < √3, also for complex v (R6)
                                        │
 §6 compact, index 0 (R7) ──────────────┤  existence and uniqueness for κ_K < √3
                                        ▼
 §7 ρ(T₁) ≤ 1/(√3 v_c), E discrete, E ⊂ {|v| ≥ √3 v_c} (R8)
 §8 ∫τ|ε|² = 0 at real v (R9)  ──►  E = ∅ when τ has one sign
```

## Step by step

### R1, R2: the reduction (§2)

With ε = (a, b), γ⁰ = i diag(1, 1, −1, −1) and γ^j = [[0, σ_j], [σ_j, 0]], D̃ε = (σ·∇b + (i/2)τa, σ·∇a − (i/2)τb) =: (u, w).
Then w + iu = σ·∇ψ − ½τψ and w − iu = σ·∇χ + ½τχ, with u = (R₁ − R₂)/2i and w = (R₁ + R₂)/2, so ε ↦ (ψ, χ) is an
invertible linear change of variables, and both equations hold for any τ, real or complex. I derived them by hand, and
the script confirms them on 400 random complex samples with a real and a complex τ alternating (worst 2×10⁻¹⁵). The
control `--mutate=g0block` (lower block of γ⁰ flipped) fails `reduction`.

The same data give a sharper statement (`block_form`). In the coordinates (w + iu, w − iu), D = diag(S, S) and
B₁ = diag(−M, +M), with S = σ·∇ and M = ½τ₁. So T₁ = D⁻¹B₁ is −S⁻¹M on the ψ sector and +S⁻¹M on the χ sector, and
the spectrum of T₁ is that of S⁻¹M together with its negative. That is the exact reason E is symmetric under v → −v (§7),
and it reduces the planted spectrum to a 2-spinor problem (R10).

Theorem 1 of [FL22] is about σ·(−i∇ − A)ψ = 0 with real A. The ψ equation here, σ·∇ψ = ½τψ, reads
σ·(−i∇)ψ = −(i/2)τψ, an imaginary scalar potential, so Theorem 1 does not apply and Remark 4 is the right tool. #63 uses
Remark 4. Its numbering was checked in the typeset PDF (Theorem 1, Remarks 2–4, Theorem 5, Remark 6, Proposition 7;
remarks share the theorem counter).

### R3: L³ decay (§3)

Each of the five steps holds.

1. g = Dε = −Bε has |g| = ½|τ||ε|. Hölder with 2/3 = 1/2 + 1/6 puts g in L^{3/2}, and with 1/2 = 1/3 + 1/6 in L².
2. D² = Δ because the γ^j anticommute and square to 1. E = D(−1/(4π|x|)) = γ·x/(4π|x|³) and DE = δ. The script checks that
   ∂_j(x_k/r³) is symmetric (so only div(x/r³) survives), that div(x/r³) = 0 off the origin, and that the flux of
   x/(4π r³) is 1 (`greens_kernel`, whose control with 2π fails it). |E| = 1/(4π|x|²) since γ·x is Hermitian with square
   |x|² (`kernel_norm`). The Fourier transform (iγ·k)⁻¹ is stated, which gives E ∗ h = D⁻¹h for h ∈ C_c^∞ directly.
3. HLS with kernel |x|⁻² in ℝ³ maps L^{3/2} to L³, since 1/q = 1/p + 2/3 − 1 = 1/3 for p = 3/2 > 1.
4. The approximation h_n → g in L² and L^{3/2}, D⁻¹ bounded from L² to Ḣ¹ ⊂ L⁶, E∗ bounded from L^{3/2} to L³, and a common
   a.e. subsequence give D⁻¹g = E ∗ g.
5. D is injective on Ḣ¹, so ε = D⁻¹g = E ∗ g ∈ L³.

Numerically, `planted_greens` evaluates (E ∗ B₁ε)(x₀) for the Frank–Loss zero mode by quadrature in polar coordinates about
x₀ = (0.3, −0.2, 0.5), where the |y|⁻² of E cancels the Jacobian. It finds D⁻¹B₁ε = λε with λ = +i to the ten printed
digits and a residual of 5×10⁻¹⁶. So the Green's-function representation, its 1/(4π), and the convolution form of T₁ are
right on a zero mode in Ḣ¹ ∩ L³ (the mode decays like r⁻², which is in L^p for every p > 3/2). The control with 1/(2π)
fails.

### R4, R5: the Frank–Loss inequality (§4)

The statement used (ψ ∈ Ḣ¹ ∩ L³, |σ·∇ψ| ≤ |A||ψ|, then ‖A‖₃² ≥ 3S₃) is [FL22] Remark 4 at d = 3, with the Ḣ¹ that G1
below says the remark needs. The display of Proposition 7 in §4 is [FL22]'s with d = 3, checked term by term against the
source (d/(d−1) = 3/2, 1/d = 1/3, (d−1)/d = 2/3, (d−1)/(d−2)² = 2, (d−2)/(d−1) = 1/2, 2(d−1) = 4). Numerically `prop7`
evaluates both sides on a smooth two-lobe spinor at ε = 1 on a grid, and they agree to 1×10⁻¹² (relative, N = 88). With the
4×4 gammas of mass_bounds.md, N = 4, they agree to 5×10⁻¹² (`prop7_4x4`), which supports #63's remark that the algebra does
not need the reduction. The control that changes 4 to 4.1 fails both.

The chain of Section 3 reads as #63 summarises it, and I redid it. Dropping the left side and using 4 − 3|ψ|²/|ψ|_ε² ≥ 1
gives 3∫|∇|ψ|_ε^{1/2}|² ≤ ∫|σ·∇ψ|²/|ψ|_ε ≤ ∫|A|²|ψ|. Sobolev on |ψ|_ε^{1/2} − ε^{1/2} and Hölder with |A|² ∈ L^{3/2} and
|ψ| ∈ L³, then ε → 0 by monotone convergence, give S‖ψ‖₃ ≤ ⅓‖A‖₃²‖ψ‖₃. Dividing by ‖ψ‖₃ ∈ (0, ∞) gives 3S ≤ ‖A‖₃². ψ ∈ L³ is
needed to divide and to make the right side finite.

The constants. S₃ = (3/4)|S³|^{2/3} = 3(π/2)^{4/3} = 5.477904089531, and the Aubin–Talenti ratio of (1 + r²)^{−1/2} is the
same number, with three other trial functions above it (`sobolev_const`). The Frank–Loss appendix pair
Ψ = (1 + r²)^{−3/2}(1 + isσ·x)φ₀ with Λ = 3s/(1 + r²) satisfies σ·(−i∇)Ψ = ΛΨ symbolically for s = ±1 and any φ₀
(`fl_pair`). ‖Λ‖₃² = 3S to 10⁻¹⁰, and κ_K(τ₁ = 6/(1 + r²), v = 1) = √3 to 10⁻¹⁰, so √3 v_c = 1 (`fl_norm`).

Near the extremal the inequality behaves as a sharp inequality should (`fl_local_min`). With a = |σ·∇ψ|/|ψ| and ψ = Ψ + δη
for a smooth localised η, ‖a‖₃² − 3S is 0 at δ = 0, 0.134 at δ = ±0.05 and 0.531, 0.534 at δ = ±0.1. It is non-negative,
equal on both sides (the first variation vanishes) and quadratic (ratio 3.97). This is *numerical evidence* for one
family and not a test of the theorem.

#### G1, discrepancy 3: Remark 4 needs ψ ∈ Ḣ¹

Proposition 7 is stated for ψ ∈ Ḣ¹, its Step 3 integrates by parts against ∇ψ ∈ L², and Step 4 uses Hardy's inequality and
∇ψ ∈ L². For Theorem 1, Section 3 obtains Ḣ¹ in its first lines, from [FL21]'s regularity theorem (ψ ∈ L^{2d/(d−2)}) and the
equation (γ·Aψ ∈ L², so γ·(−i∇)ψ ∈ L², so ψ ∈ Ḣ¹). Remark 4 lists only ψ ∈ L^{2d/(d−1)} and the inequality, and says
"inspection of the proof shows". So the remark as printed omits a hypothesis its proof uses, and discrepancy 3 in #63 is
stated correctly.

Does Ḣ¹ follow from the inequality alone? I believe so, by the same argument. [FL21]'s regularity proof writes
ψ = Γ ∗ (σ·Aψ), splits σ·A = F_ε + G_ε with F_ε bounded and compactly supported and ‖G_ε‖₃ small, and bootstraps with HLS
and a duality argument. Its estimates use only |G_ε| and |F_ε|, and Hermiticity enters only in writing the adjoint. In
the inequality case g := σ·(−i∇)ψ has |g| ≤ |A||ψ|, so g = Mψ with M := gψ†/|ψ|² (zero where ψ = 0) and |M| ≤ |A|, and
the proof runs with M in place of σ·A. This is my reading and not a written proof, so it is *conditional* on that reading.
It does not matter for #63, because there ε ∈ Ḣ¹ is the hypothesis, so ψ, χ ∈ Ḣ¹ by construction (§2), and they are in
L³ by §3.

#### G2, a gap in the proof of [FL22] Proposition 7, Step 3

Step 3 proves a preliminary identity for fixed bounded, compactly supported f with ∇f ∈ L^d, which by density in Ḣ¹ is
enough because (∇f)ψ ∈ L² then. It then takes f = χ|ψ|_ε^{−2/(d−1)} and states that this "also implies ∇f ∈ L^d". For a
general ψ ∈ Ḣ¹ it does not, because ∇f contains ∇|ψ|_ε ∈ L² only. The identity of Proposition 7 is nonetheless true for
every ψ ∈ Ḣ¹, by a routine density argument that [FL22] do not write. Every term of the identity is an integral of a
quadratic form in ∇ψ whose coefficients are bounded continuous functions of ψ (because |ψ|_ε ≥ ε), so both sides are
continuous on Ḣ¹ (an a.e. convergent subsequence and dominated convergence), and the identity extends from C_c^∞, where f
is smooth with compact support, to Ḣ¹. The numerical check of Proposition 7 above is on smooth spinors and does not test
this extension. It matters for #63 because ψ and χ need not be smooth. It is a gap in exposition and not in truth, and
#63 cites [FL22] rather than re-proving it.

### R6: no kernel below √3 (§5)

Correct as written. ε ∈ Ḣ¹ with D̃ε = 0 gives ψ, χ ∈ Ḣ¹ ∩ L³ with |σ·∇ψ| = ½|τ||ψ|, and Remark 4 with |A| = ½|τ|e₁ gives
¼‖τ‖₃² ≥ 3S for ψ ≢ 0, that is κ_K ≥ √3. The arithmetic ‖τ‖₃ ≥ 2√(3S) ⇔ κ_K ≥ √3 is right. The complex-amplitude
statement is also right, because the reduction and the L³ step use only |τ|. In the planted case a kernel does appear at
|v| = 1 = √3 v_c exactly (`fl_pair`, `planted_greens`). The argument-principle count finds no kernel at |ℓ| < 0.95 in
sector 1 and at |ℓ| < 1.4 in sector 2 (`planted_nocomplex`), and the real-root scan finds none below |ℓ| = 1 in sectors 1
to 6 (`planted_ode`). So in those sectors the bound is attained and not exceeded.

### R7: Fredholm index 0, existence (§6)

Correct. Multiplication by τ ∈ L³ is compact from Ḣ¹ to L² (bounded, compactly supported τ by Hölder and
Rellich–Kondrachov, then a norm limit with ‖M_τ − M_{τ_m}‖ ≤ ‖τ − τ_m‖₃/√S). D⁻¹: L² → Ḣ¹ is bounded, so T = D⁻¹B is compact on
Ḣ¹ and D̃ = D(I + T) is Fredholm of index 0. For ε = ε₀ + ψ, D̃ψ = −Bε₀ with Bε₀ ∈ L² because τ ∈ L². Existence and
uniqueness follow when the kernel is trivial. The planted spectrum is consistent with compactness, since the eigenvalues
±i/ℓ_n with ℓ_n = 1 + 2n/3 accumulate only at 0.

### R8: the spectrum of T₁ and the set E (§7)

Correct. λ ≠ 0 is an eigenvalue iff D̃(−1/λ) has a kernel, and by `block_form` that kernel is
ker(S − ½vτ₁) ⊕ ker(S + ½vτ₁). §5 for complex v gives |λ| ≤ 1/(√3 v_c), hence ρ(T₁) ≤ 1/(√3 v_c). The nonzero spectrum of
a compact operator is a set of eigenvalues accumulating only at 0, so E ⊂ ℝ has no finite accumulation point and lies in
|v| ≥ √3 v_c. E is symmetric under v → −v because the two summands swap. The Picard statement follows from Gelfand's
formula, since the linear part of the iteration is −vT₁.

### R9: the identity at real v (§8)

Correct. The pointwise identity 2 Re⟨ε, γ⁰D̃ε⟩ = ∂_jξ^j − τ|ε|² holds to 7.5×10⁻¹⁵ at 40 points for random smooth ε and a
real τ of either sign (`divxi`, which the `g0block` control breaks, because γ⁰γ^j is then not Hermitian). The cut-off
argument is sound. ξ ∈ L^{3/2} since |ξ| ≤ |ε|² and ε ∈ L³ (§3), τ|ε|² ∈ L¹, and the shell term is at most
(C/R)‖ξ‖_{L^{3/2}(|x|>R)} (vol of the shell)^{1/3} ≤ C′‖ξ‖_{L^{3/2}(|x|>R)} → 0. For one sign of τ, τ|ε|² = 0 a.e. gives
Bε = 0, then Dε = 0, then ε = 0.

## Remarks on #63 (not gaps)

1. The reduction is not needed for the inequality, as #63 says. Proposition 7 and Section 3 use only Hermitian
   anticommuting gammas (`prop7_4x4`), and |Dε| = ½|τ||ε| holds for the 4-spinor too. [FL22] state N = 2^{[d/2]}, so
   applying it to the reducible 4×4 representation is a reading of the proof, which the reduction avoids. What the
   reduction does earn is `block_form`.
2. §8's identity is true sector by sector, which is stronger and shorter. For real v and τ₁,
   2 Re⟨ψ, (σ·∇ − ½vτ₁)ψ⟩ = ∂_j(ψ†σ_jψ) − vτ₁|ψ|² (`sector_identity`), so the same cut-off argument gives ∫τ₁|ψ|² = 0 and
   ∫τ₁|χ|² = 0 separately for a zero mode at real v ≠ 0, where §8 has their sum. For one-signed τ₁ this says that the
   weighted problem σ·(−i∇)ψ = ½ℓτ₁ψ is self-adjoint with a positive weight, so ℓ is real and the amplitude v = iℓ is
   imaginary. That is what the planted case shows (R10), and why "FL's extremal reaches the bound only at imaginary v" (§8).

## The independent check

The script shares no code with `frank_loss.py`. What the two share is inputs and not machinery: the gammas and conventions
of mass_bounds.md, the sharp S, and the formulae of [FL22]. If those conventions were wrong, both would fail together. The
sign of ∇̃ is #59's and `conventions.py`'s subject and is not examined here. The planted spectrum comes from an ODE and not
from a box discretisation, with no FFT, Krylov method or ARPACK. The grids that remain (`prop7`, `fl_local_min`) are
quadrature of analytic integrands with analytic derivatives.

### The planted case, without a grid (`planted_ode`)

By `block_form`, T₁ is −S⁻¹M on ψ and +S⁻¹M on χ with M = 3/(1 + r²). In a sector of total angular momentum j = k − ½ write
ψ = f Ω_k + g Ω_{−k}, with Ω_{±k} the spin-angular functions, σ·r̂ Ω_k = −Ω_{−k} and σ·L Ω_k = −(1 + k)Ω_k. With μ = iℓ and
f = iF, σ·∇ψ = μWψ, W = 3/(1 + r²), becomes the real system

    g′ + (1 − k)g/r = ℓ W F,        F′ + (1 + k)F/r = −ℓ W g,

regular at 0 as (g, F) ~ (r^{k−1}, r^k) and decaying at ∞ as (g, F) ~ (r^{−k−2}, r^{−k−1}), the other branch growing as
r^{k−1}. The script integrates both solutions to r = 1 in t = ln r, where each integration runs in its stable direction, and
takes the matching determinant. Its zeros ℓ are the eigenvalues, with λ = i/ℓ in the ψ sector (−i/ℓ in χ), so |λ| = 1/|ℓ|.

| k | roots ℓ found in (0, 4.2] | (2k + 1 + 2m)/3 |
|---|---|---|
| 1 | 1, 5/3, 7/3, 3, 11/3 | 1, 5/3, 7/3, 3, 11/3 |
| 2 | 5/3, 7/3, 3, 11/3 | same |
| 3 | 7/3, 3, 11/3 | same |
| 4, 5, 6 | lowest only: 3, 11/3, 13/3 | (2k + 1)/3 |

All roots agree with the formula to 1.4×10⁻¹². So the eigenvalues of T₁ are λ = ±i/ℓ_n with ℓ_n = 1 + 2n/3, n = 0, 1, 2, …,
and the largest modulus is **ρ(T₁) = 1/ℓ₀ = 1.000000000000**, from the sector k = 1. This is Frank–Loss's bound attained,
ρ(T₁) = 1/(√3 v_c) with √3 v_c = 1. The values are the Dirac eigenvalues ±(3/2 + n) of the round S³ (stereographic
projection, ℓ = (2/3)(3/2 + n)), which I recall and have not verified *(unverified)*. The ODE computation does not use
that spectrum and reproduces it in sectors 1 to 6. Each sector k holds 2k states, so the multiplicity of one sign is
Σ_{k ≤ n+1} 2k = (n + 1)(n + 2), and ψ and χ together double it.

| n | ℓ_n | \|λ\| = 1/ℓ_n | multiplicity of +i/ℓ_n | of −i/ℓ_n |
|---|---|---|---|---|
| 0 | 1 | 1 | 4 | 4 |
| 1 | 5/3 | 0.6 | 12 | 12 |
| 2 | 7/3 | 0.4286 | 24 | 24 |

This matches #63 §9, which finds eight eigenvalues ±i (fourfold each) and the ninth at 3/5. #63's 0.9910 at L = 10, N = 64 is
therefore a 0.9 % discretisation error on an exactly known 1, as #63 says.

`planted_nocomplex` counts the zeros of the same determinant in the complex ℓ plane by the argument principle (360-point
circles). The k = 1 determinant has winding number 0 on |ℓ| = 0.95 and 2 on |ℓ| = 1.3, and the k = 2 determinant 0 on
|ℓ| = 1.4 and 2 on |ℓ| = 1.95. So no eigenvalue lies off the real ℓ axis in those discs, as the self-adjoint structure of
remark 2 predicts, without using it.

### Checks and controls

| check | step | result |
|---|---|---|
| `reduction`, `reduction_norms`, `block_form` | R1, R2 | worst 2×10⁻¹⁵ over 400 samples |
| `clifford`, `divxi`, `sector_identity` | R9 | residuals 7.5×10⁻¹⁵ and 1.8×10⁻¹⁵ |
| `greens_kernel`, `kernel_norm`, `planted_greens` | R3 | flux 1, λ = +i, residual 5.5×10⁻¹⁶ |
| `sobolev_const`, `fl_pair`, `fl_norm` | R4 | S = 5.477904089531, ‖Λ‖₃² = 3S, κ_K = √3 |
| `prop7`, `prop7_4x4` | R4 | 1.1×10⁻¹², 4.8×10⁻¹² |
| `fl_local_min` | R4 | excess 0, 0.134, 0.531 |
| `planted_ode`, `planted_nocomplex` | R10, R6, R8 | ρ(T₁) = 1.000000000000, no complex root in the discs |
| `v_c` | R11 | worst difference 4.6×10⁻⁶, below the printed 5×10⁻⁶ |

Each control runs the same code with one input broken, and counts as fired only if its target checks fail.

| `--mutate=` | what it breaks | checks that fail |
|---|---|---|
| `g0block` | lower block of γ⁰ flipped | `reduction`, `divxi` (and `block_form`, `clifford`) |
| `half` | B without its ½ | `planted_ode`, `planted_greens` |
| `ode` | (2 + k)/r in place of (1 + k)/r | `planted_ode` |
| `sobolev` | 21/20 of the sharp S | `sobolev_const`, `fl_norm` |
| `greens` | 1/(2π) in the kernel | `greens_kernel`, `planted_greens` |
| `prop7` | 4.1 in place of 4 in Proposition 7 | `prop7`, `prop7_4x4` |
| `vcthree` | 2π in place of π in the angular integral | `v_c` |
| `extremal` | drops ∂(1 + iσ·x) from the extremal's gradient | `fl_local_min` |

All eight fired.

## Not reviewed

§1 (the setting) and §9 beyond the planted case and the v_c table, which means the ARPACK, GMRES and box numerics for the
tanh profile and the claim that E is empty for it (a *conjecture* in #63). Also the correspondence with [FNS11] in §2, the
HLS, Rellich–Kondrachov and Riesz–Schauder statements (standard, exponents checked here, statements not re-read), and
mass_bounds.md items other than the last sentences of item 5. [FL22]'s proof was read and not re-proved, apart from the
checks above, and the repair in G2 is described and not written out.

## Suggested status text (not edited here)

- mass_bounds.md item 5, "(*proved*, by hand, not yet independently reviewed)", can become "*proved*, by hand, given [FL22],
  and confirmed by #68's review ([review_fl.md](review_fl.md))". The dependence on [FL22] stays, with G1 and G2.
- prior_work.md §4 ("A possible sharpening, not checked here", lines 105–112, "this note's reading, not checked") and
  discrepancy 2 ("pending an independent review") are out of date since #63 and now since #68.
- frank_loss.md's "not yet independently reviewed" (Answer, §5, §6, §7, §8 and the closing list) can be replaced the same way.

## Recommended tasks

### Task 1, close G1 and G2 by hand

Write out the L^{2d/(d−2)} regularity for |γ·∇ψ| ≤ |A||ψ| through M = gψ†/|ψ|², and the density step of Proposition 7 on
Ḣ¹. These are the two places where the chain rests on [FL22] read and not proved. Acceptance, each can fail. The
density step comes with a numerical continuity test, both sides of Proposition 7 evaluated along smooth ψ_n → ψ in Ḣ¹ for
the non-smooth ψ = |x|^{−1/4}(1 + |x|²)^{−3/2}φ₀ (in Ḣ¹ ∩ L³, unbounded at 0), which fails if the two sides do not
converge to the same limit. The regularity step lists every hypothesis of [FL21]'s argument and where Hermiticity is used,
and fails if any step needs more than the adjoint. Depends on #68.

### Task 2, an independent computation of the tanh-profile spectrum

τ₁ = −f′(r) cos θ is axisymmetric about x, so with `block_form` the problem S⁻¹M reduces to a 2D problem per eigenvalue of
J_x. Compute its leading eigenvalues and the real-v zero modes by a method that is not a periodic FFT box (finite
elements, or a spherical-harmonic expansion with the radial system of `planted_ode` coupled across l). The claims
"ρ(T₁)v_c = 0.40938 at σR = 2" and "E is empty for the tanh profile" rest on one discretisation. Acceptance: ρ(T₁)v_c
agrees with 0.4094 and 0.3674 to three digits at σR = 2 and 4, the same code returns ρ = 1 for τ₁ = 6/(1 + r²), and the
sector conditions ∫τ₁|ψ|² = ∫τ₁|χ|² = 0 are evaluated on every real-λ candidate. Depends on #63, #68.

### Task 3, read [FNS11] Eq. 2.1

Compare the 2-spinor Witten operator with the reduction in §2 line by line, which #63 says was not done. The sentence
"parallels the 2-spinor Witten operator of [FNS11]" is unchecked. Acceptance: it is either confirmed against the equation,
with its sign convention, or corrected in frank_loss.md. Depends on #63.

## References

- [FL22] R. L. Frank, M. Loss, "A sharp criterion for zero modes of the Dirac equation", arXiv:2201.03610v1 (2022), checked
  against the arXiv API. Theorem 1, Remark 4, Proposition 7, Section 3 and the Appendix read from the arXiv source.
  Journal version *(unverified)*, as in [prior_work.md](prior_work.md).
- [FL21] R. L. Frank, M. Loss, "Which magnetic fields support a zero mode?", arXiv:2012.13646 (2020), checked against the
  arXiv API. The regularity theorem and its proof read from the arXiv source. Journal version *(unverified)*.
- [FNS11] J. Frauendiener, J. M. Nester, L. B. Szabados, arXiv:1105.5008, as in [../literature.md](../literature.md), not re-read.
- The Dirac spectrum ±(3/2 + n) of the round S³ is recalled and not read *(unverified)*. The ODE computation does not use it.

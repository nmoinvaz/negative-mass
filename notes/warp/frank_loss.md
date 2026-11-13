# The Witten-spinor existence threshold, sharpened with Frank–Loss

Issue #63. Script [frank_loss.py](frank_loss.py), output [frank_loss_output.txt](frank_loss_output.txt). Run from the
repository root; the commands that wrote the output are at its top. `python notes/warp/frank_loss.py --controls` runs
every `--mutate=` control (`g0block`, `gamma`, `ly`, `sobolev`, `half`, `vscale`, `Jsign`) and exits 1 unless each fires.

## Answer

The √3 extension holds. On a flat slice with tr K ∈ L² ∩ L³, D̃ = D + ½(tr K)γ⁰ has no kernel in Ḣ¹ when
‖tr K‖_{L³} < 2√(3S), that is κ_K < √3, and the Witten spinor then exists and is unique in ε₀ + Ḣ¹. Status: *proved*,
by hand, from the inequality of Frank and Loss [FL22, Section 3], which is cited and not re-proved; not yet
independently reviewed. The 4-to-2 reduction it rests on is *checked symbolically*. Every v_c of
[mass_bounds.md](mass_bounds.md) item 5 grows by √3 for existence, while the contraction (and its uniqueness and
convergence estimate) stays at κ_K < 1.

| σR | v_c, contraction (#58) | √3 v_c, Frank–Loss | thin-wall formula × √3 |
|---|---|---|---|
| 2 | 3.74490 | 6.48636 | 6.82624 |
| 4 | 2.47097 | 4.27985 | 4.30026 |
| 8 | 1.56243 | 2.70621 | 2.70900 |
| 16 | 0.98503 | 1.70612 | 1.70656 |

Two further statements answer #58's question of what happens beyond v_c, both *proved* by hand and not reviewed.

- For every real v outside a discrete set E ⊂ {|v| ≥ √3 v_c}, the Witten spinor exists and is unique (§7). E is the set
  of speeds at which D̃ has a zero mode.
- A zero mode at real v satisfies ∫ tr K |ε|² dV = 0 (§8). So if tr K ≥ 0 everywhere, or ≤ 0, E is empty and the
  spinor exists for every amplitude. Alcubierre's tr K = −v f′ cos θ changes sign, so this does not apply to it.

For the tanh profile at σR = 2 and 4 the numerics (*numerical evidence*, §9) find no real zero mode. GMRES solves the
Witten equation to a residual below 10⁻¹⁰ at every v tried in (v_c, √3 v_c) and at about 3.4 v_c.

## 1. Setting

The slice is ℝ³ with the flat metric. τ := tr K is a real function in L² ∩ L³. Spinors are ℂ⁴-valued, with the gammas of
mass_bounds.md: γ⁰ = i diag(1, 1, −1, −1), γ^j = [[0, σ_j], [σ_j, 0]]. Ḣ¹ is the completion of C_c^∞ in ‖∇u‖₂; by
Sobolev's inequality it is {u ∈ L⁶ : ∇u ∈ L²}, which is the space [FL22] call Ḣ¹ (their Section 2). D = γ^j∂_j maps
Ḣ¹ isometrically onto L² (its symbol iγ·k has |iγ·k ĥ| = |k||ĥ|). B := ½τγ⁰ and D̃ := D + B.

By Hölder and Sobolev, ‖Bu‖₂ = ½‖τ|u|‖₂ ≤ ½‖τ‖₃‖u‖₆ ≤ κ_K‖∇u‖₂ with κ_K = ‖τ‖₃/(2√S), so D̃ : Ḣ¹ → L² is bounded.
For Alcubierre data τ = vτ₁ with τ₁ = −f′(r) cos θ, and κ_K = v/v_c.

## 2. The reduction (*checked symbolically*)

Write ε = (a, b) with a, b ∈ ℂ². Then D̃ε = (u, w) with u = σ·∇b + (i/2)τa and w = σ·∇a − (i/2)τb. With
ψ = a + ib and χ = a − ib,

    w + iu = σ·∇ψ − ½τψ,    w − iu = σ·∇χ + ½τχ,

and conversely u = (R₁ − R₂)/(2i), w = (R₁ + R₂)/2 for R₁, R₂ the right sides. So, for any ε,

    D̃ε = 0   ⟺   σ·∇ψ = ½τψ  and  σ·∇χ = −½τχ,

and |ψ|² + |χ|² = 2|ε|², so ε = 0 iff ψ = χ = 0. Both are identities in sympy for an arbitrary real τ(x, y, z) and
arbitrary complex a, b (`reduction`, `norms`). Flipping the sign of the lower block of γ⁰ (`--mutate=g0block`) makes
`reduction` fail: then w + iu = σ·∇ψ − ½τχ, and the two equations couple. The reduction is the one #60 wrote by hand
(prior_work.md §4). It is the 2-spinor form of the Witten operator, in which only tr K enters, as [FNS11] note (their
Eq. 2.1, recorded in [../literature.md](../literature.md), not re-read here).

If ε ∈ Ḣ¹, then ψ, χ ∈ Ḣ¹(ℝ³, ℂ²) and, a.e., |σ·∇ψ| = ½|τ||ψ| and |σ·∇χ| = ½|τ||χ|. These are the pointwise
inequalities of [FL22] Remark 4 with |A| = ½|τ|, for the Pauli matrices that [FL22] use in d = 3 (N = 2^{[3/2]} = 2).

## 3. Zero modes lie in L³ (*proved*, by hand)

Claim. If τ ∈ L² ∩ L³ and ε ∈ Ḣ¹ with D̃ε = 0, then ε ∈ L³, with ‖ε‖₃ ≤ C‖τ‖₂‖ε‖₆.

Proof.
1. g := Dε = −Bε. Since ε ∈ L⁶ and τ ∈ L², Hölder (2/3 = 1/2 + 1/6) gives g ∈ L^{3/2}. Since τ ∈ L³, g ∈ L² (§1).
2. D² = Δ, because the γ^j anticommute and square to 1. The fundamental solution of Δ is −1/(4π|x|), so
   E(x) := D(−1/(4π|x|)) = γ·x/(4π|x|³) satisfies DE = δ. Its Fourier transform is −iγ·k/|k|² = (iγ·k)⁻¹, the symbol of
   D⁻¹. As a matrix, γ·x̂ is Hermitian with square 1, so |E(x)| = 1/(4π|x|²).
3. The Hardy–Littlewood–Sobolev inequality in ℝ³ with kernel |x|⁻² maps L^{3/2} to L³ (1/3 = 2/3 − 1/3). So E ∗ g ∈ L³
   with ‖E ∗ g‖₃ ≤ C‖g‖_{3/2}.
4. E ∗ h = D⁻¹h for h ∈ C_c^∞ (step 2). Take h_n ∈ C_c^∞ with h_n → g in L² and in L^{3/2} (truncate, then mollify).
   Then D⁻¹h_n → D⁻¹g in Ḣ¹, hence in L⁶, and E ∗ h_n → E ∗ g in L³. A subsequence converges a.e. in both, so
   D⁻¹g = E ∗ g a.e.
5. D is injective on Ḣ¹, so ε = D⁻¹(Dε) = D⁻¹g = E ∗ g ∈ L³. ∎

So ψ, χ ∈ Ḣ¹ ∩ L³. The decay assumption is the one item 5 already makes (τ ∈ L² for the source Bε₀ ∈ L²). For the tanh
profile τ decays exponentially and lies in every L^p. Without some decay of τ beyond L³, step 1 fails: Ḣ¹ alone gives
only L⁶, and a function like (1 + r²)^{−1/2} is in Ḣ¹ but not in L³.

## 4. The Frank–Loss inequality (cited, read from arXiv:2201.03610)

What [FL22] prove, in the form used here (d = 3, N = 2). Let ψ ∈ Ḣ¹(ℝ³, ℂ²) ∩ L³, A ∈ L³(ℝ³, ℝ³), and
|σ·∇ψ| ≤ |A||ψ| a.e. If ψ ≢ 0, then ‖A‖₃² ≥ 3S₃, with S₃ = (3/4)|S³|^{2/3} = 3(π/2)^{4/3} the S of item 5 (`S_3`).

Where it is. Theorem 1 states it for zero modes of σ·(−i∇ − A), and Remark 4 extends it to the pointwise inequality
for ψ ∈ L^{2d/(d−1)} = L³. The proof is Section 3, from Proposition 7. For d = 3 and every ε > 0, Proposition 7 gives,
for ψ ∈ Ḣ¹ and |ψ|_ε := (|ψ|² + ε²)^{1/2},

    ∫ Σ_j |[−i∂_j − ⅓σ_j σ·(−i∇)] (ψ/|ψ|_ε^{3/2})|² |ψ|_ε² = ⅔ ∫ |σ·∇ψ|²/|ψ|_ε − 2 ∫ |∇|ψ|_ε^{1/2}|² (4 − 3|ψ|²/|ψ|_ε²).

Section 3 drops the left side (≥ 0), uses 4 − 3|ψ|²/|ψ|_ε² ≥ 1, bounds |σ·∇ψ|² ≤ |A|²|ψ|² and |ψ|²/|ψ|_ε ≤ |ψ|, then
applies Sobolev to |ψ|_ε^{1/2} − ε^{1/2} and Hölder, ∫|A|²|ψ| ≤ ‖A‖₃²‖ψ‖₃. Letting ε → 0 by monotone convergence,
S‖ψ‖₃ ≤ ⅓‖A‖₃²‖ψ‖₃, and ψ ≢ 0 with ‖ψ‖₃ < ∞ gives the claim.

What the proof needs, as read. ψ ∈ Ḣ¹ (for Proposition 7 and its integrations by parts), ψ ∈ L³ (to divide by ‖ψ‖₃),
and the pointwise inequality. In Theorem 1 the Ḣ¹ membership comes from the equation, through FL's earlier paper;
Remark 4 does not say how it is obtained under the inequality alone. Our ψ and χ are in Ḣ¹ by construction (§2), so
this point does not arise here. The algebra of Proposition 7 uses only the anticommutation relations of Hermitian
gammas: the Penrose formula (their first display in the proof) and the cancellation in their Step 3. The script checks
the Penrose formula for the 2×2 and the 4×4 gammas (`penrose`), so the reduction is not needed for the algebra; it is
used so that the statement applies exactly as written, with N = 2.

The constant is attained, so √3 cannot be improved by this route. The Loss–Yau pair (FL Section 1) solves
σ·(−i∇ − 𝒜)Ψ = 0 with ‖𝒜‖₃² = 3S (`loss_yau`, `loss_yau_norm`). FL's appendix pair
σ·(−i∇)Ψ̃ = ΛΨ̃, Λ = 3s/(1 + r²), is a zero mode of D̃ with τ = 6/(1 + r²) at the imaginary speed v = is, and that τ has
κ_K(1) = √3, so √3 v_c = 1 = |v| exactly (`fl_scalar_extremal`). At real v that configuration is not available (§8).

## 5. No kernel when κ_K < √3 (*proved*, by hand, given §4)

Let ε ∈ Ḣ¹ with D̃ε = 0. By §2, ψ = a + ib and χ = a − ib are in Ḣ¹(ℝ³, ℂ²) and satisfy |σ·∇ψ| = ½|τ||ψ|,
|σ·∇χ| = ½|τ||χ| a.e. By §3 both are in L³. Apply §4 with A = ½|τ| e₁ ∈ L³. If ψ ≢ 0, then ¼‖τ‖₃² ≥ 3S, that is
‖τ‖₃ ≥ 2√(3S), that is κ_K ≥ √3; likewise for χ. So κ_K < √3 forces ψ = χ = 0, and ε = 0 by §2. The same holds for
complex amplitudes: with τ = vτ₁ and v ∈ ℂ, the reduction is linear and |A| = ½|v||τ₁|, so D̃(v) has no kernel when
|v| < √3 v_c.

## 6. Fredholm index 0, and existence (*proved*, by hand)

1. Multiplication by τ ∈ L³ is compact from Ḣ¹ to L². For τ bounded with support in a ball B: a bounded sequence in Ḣ¹
   is bounded in H¹(B) (its L²(B) norm is controlled by L⁶), so by Rellich–Kondrachov a subsequence converges in L²(B),
   and τu_n converges in L². For general τ ∈ L³, the truncations τ_m = τ 1{|x| < m, |τ| < m} → τ in L³ by dominated
   convergence, and ‖M_τ − M_{τ_m}‖_{Ḣ¹→L²} ≤ ‖τ − τ_m‖₃/√S → 0 (§1). A norm limit of compact operators is compact.
2. So T := D⁻¹B is compact on Ḣ¹, and D̃ = D(I + T). By the Riesz–Schauder theory, I + T is Fredholm of index 0 on Ḣ¹,
   and D is an isomorphism Ḣ¹ → L², so D̃ : Ḣ¹ → L² is Fredholm of index 0: injective iff surjective.
3. The Witten equation for ε = ε₀ + ψ with constant ε₀ is D̃ψ = −Bε₀, and Bε₀ = ½τγ⁰ε₀ ∈ L² because τ ∈ L². When D̃ has
   no kernel, it is a bijection with bounded inverse, so ψ = −D̃⁻¹Bε₀ ∈ Ḣ¹ exists and is unique.

With §5: κ_K < √3 gives a unique Witten spinor in ε₀ + Ḣ¹. This is the sharpened item 5. The contraction of item 5
gives more where it applies (κ_K < 1): convergence of the plain iteration at the rate κ_K per step in Ḣ¹ and the bound
‖∇ψ‖₂ ≤ ‖Bε₀‖₂/(1 − κ_K).

## 7. Beyond √3 v_c: the spectrum of T₁ (*proved*, by hand)

With τ = vτ₁, T(v) = vT₁ and T₁ := D⁻¹B₁ (B₁ = ½τ₁γ⁰) is compact on Ḣ¹. Its nonzero spectrum consists of eigenvalues of
finite multiplicity that accumulate only at 0, and λ ≠ 0 is an eigenvalue iff D̃(−1/λ) has a kernel.

- **Spectral radius.** By §5 (complex v), every eigenvalue has |λ| ≤ 1/(√3 v_c). So ρ(T₁) ≤ 1/(√3 v_c), while
  ‖T₁‖_{Ḣ¹} ≤ 1/v_c is all the contraction gives, and #59 found that norm bound attained for some τ. The Picard iteration
  ψ ← −vT₁(ε₀ + ψ) of #58 converges whenever ρ(vT₁) < 1 (Gelfand's formula), so in particular for v < √3 v_c, though
  not monotonically.
- **Real speeds.** D̃(v) fails to be invertible exactly when −1/v is a real eigenvalue of T₁. These speeds form a set E
  with no finite accumulation point, contained in |v| ≥ √3 v_c and symmetric under v → −v (the ψ and χ sectors of §2
  swap). For v ∉ E the Witten spinor exists and is unique (§6). At v ∈ E uniqueness fails, and existence depends on
  whether Bε₀ lies in the range.

Whether E is empty for the tanh profile is not proved. The numerics of §9 find no real eigenvalue among the largest ones.

## 8. A zero mode at real v has ∫ tr K |ε|² = 0 (*proved*, by hand; the pointwise identity *checked symbolically*)

With ξ^j = ⟨ε, γ⁰γ^jε⟩ (the vector of item 6), for any smooth ε,

    2 Re⟨ε, γ⁰D̃ε⟩ = ∂_j ξ^j − τ|ε|²,

because γ⁰γ^j is Hermitian and (γ⁰)² = −1 (`divxi`; `--mutate=g0block` makes it fail). The product rule holds for
ε ∈ H¹_loc, so the identity holds in L¹_loc. For a zero mode ε ∈ Ḣ¹ ∩ L³ (§3), |ξ| ≤ |ε|² gives ξ ∈ L^{3/2}, and
∂_jξ^j = τ|ε|² ∈ L¹. Integrate against a cut-off φ(x/R) with ∇φ supported in 1 ≤ |x| ≤ 2. The boundary term is at most
(C/R)‖ξ‖_{L^{3/2}(|x|>R)} (vol of the shell)^{1/3} ≤ C′‖ξ‖_{L^{3/2}(|x|>R)} → 0, so ∫τ|ε|² dV = 0. If τ ≥ 0 a.e. (or ≤ 0) and τ ≢ 0, then τ|ε|² = 0 a.e.,
so Bε = 0 and Dε = 0, and ε = 0 because D is injective on Ḣ¹. Hence for a slice whose tr K has one sign, D̃(v) has no
kernel at any real v, and the Witten spinor exists and is unique for every amplitude.

The same identity explains why FL's extremal τ = 6/(1 + r²) > 0 reaches the bound only at imaginary v (§4). For real v
and sign-changing τ the true first zero-mode speed can lie above √3 v_c, and for the tanh profile it does (§9).

## 9. Numerics (*numerical evidence*)

TBD

## What this does and does not show

TBD

## Discrepancies and reading notes

1. mass_bounds.md item 5 says "Beyond v_c this says nothing either way". That stays true of the contraction, but
   existence and uniqueness now hold for κ_K < √3 (§5–§6) and for all real v outside the discrete set E (§7). The item
   is not edited here; #62 is editing that file.
2. prior_work.md, discrepancy 2 and §4, said [FL22] "may" raise every v_c by √3, with the reduction, the L³ decay and the
   Fredholm step unchecked. They are checked here (§2, §3, §6), and the √3 holds. Its sentence on compactness and index 0
   is correct as written.
3. [FL22] Remark 4 states its hypothesis as ψ ∈ L³ with |γ·∇ψ| ≤ |A||ψ|. The Section 3 proof also uses ψ ∈ Ḣ¹, which
   Theorem 1 gets from the equation. The remark does not say how Ḣ¹ follows from the inequality alone. Here ψ ∈ Ḣ¹ is
   given, so nothing depends on it. This is a reading note on [FL22], not an error found in it.
4. #58's column |D̃ε|/|Bε₀| (1.8×10⁻¹⁰ at σR = 2, v = 1, N = 128) is its Picard stopping tolerance, not a floor of the
   discretisation. GMRES reaches 2.6×10⁻¹³ there at N = 96 (§9).

## References

- [FL22] R. L. Frank, M. Loss, "A sharp criterion for zero modes of the Dirac equation", arXiv:2201.03610 (2022);
  Theorem 1, Remark 4, Proposition 7 and Section 3 read from the arXiv source on 2026-11-13. Journal version
  *(unverified)*, as in [prior_work.md](prior_work.md).
- [FNS11] J. Frauendiener, J. M. Nester, L. B. Szabados, arXiv:1105.5008, CQG 28, 185004 (2011), as recorded in
  [../literature.md](../literature.md).
- Hardy–Littlewood–Sobolev, Rellich–Kondrachov, Riesz–Schauder and Gelfand's formula are standard and were recalled, not
  re-read (for example E. H. Lieb, M. Loss, *Analysis*, 2nd ed., AMS 2001, Theorem 4.3 for HLS *(unverified)*).

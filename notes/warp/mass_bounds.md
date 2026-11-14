# Alcubierre initial data and the signed-source mass bounds

Issue #58. Script [alcubierre_bounds.py](alcubierre_bounds.py), output [alcubierre_bounds_output.txt](alcubierre_bounds_output.txt),
run from the repository root as `python notes/warp/alcubierre_bounds.py` (about 70 s, exits 1 if a check fails). Every
check has a `--mutate=` control that makes it fail: `beta2`, `adm_m`, `adm_p`, `split_old`, `quad`, `Jsign`, `norm`.

## Setup

The Alcubierre metric [Alc94] is ds² = −dt² + (dx − v f(r_s) dt)² + dy² + dz², r_s² = (x − vt)² + y² + z², with
Alcubierre's profile f(r) = [tanh σ(r + R) − tanh σ(r − R)] / (2 tanh σR) and R = 1. Its t = 0 slice is flat, g_ij = δ_ij,
with lapse 1, shift β_x = −v f and K_ij = ½(∂_iβ_j + ∂_jβ_i). It is the K ≠ 0 initial data set that phase 4 lacks
(phase4/RESULTS.md §1, #13). μ and J are its Eulerian energy and momentum densities, with J_i = −n^a T_ai.

## Results

1. **The sources** (*checked symbolically*, by the 3+1 constraints and by G_ab of the 4-metric, which agree exactly):
   μ = −(v²/32π)((y² + z²)/r_s²) f′² and J_i = (v/16π)(∂_x∂_i f − δ_ix ∇²f). On the x-axis inside the wall μ = 0 while
   J_x = −v f′/(8π r_s) ≠ 0, so neither T nor −T satisfies μ ≥ |J|. μ is Alcubierre's Eq. (19) [Alc94], and J is the
   Eulerian flux (1/16π)(∇×∇×w)_i, w = −β, of [SSV22] Eqs. 4.7–4.8.
2. **E = P = 0** (*proved*). E's surface integral vanishes identically, since g_ij = δ_ij. P's integrand is linear in
   K, which decays exponentially, so its limit is 0 ([review.md](review.md) gap 7). Both are [BBV26] Theorem IV.19
   (iii). The surface integrals agree (*numerical evidence*): P is 1.8×10⁻³ at radius 2, 8×10⁻¹⁰ at 4 and 0 at 8
   (σR = 4, v = 1). The same code returns m ψ(R)³ for Schwarzschild and P for Bowen–York data. So these data cannot
   violate any lower bound on M. They test the bounds' hypotheses only.
3. **The minimal split** (*proved*, by hand, and checked numerically). Over splits T = T₊ − T₋ with μ± ≥ |J±|,
   the pointwise minimum is μ₋ = max(0, −μ, (|J| − μ)/2). The lower bound follows from μ₊ = μ + μ₋ ≥ 0 and
   |J| ≤ |J₊| + |J₋| ≤ μ + 2μ₋. It is attained with J₋ = −μ₋ Ĵ when |J| ≥ −μ, and with J₋ = −J otherwise. Along it,
   μ₋ + |J₋| = max(0, |J| − μ). It agrees with a Nelder–Mead minimisation at 300 random points to 2×10⁻¹³. The
   formula first written in #58, max(0, (|J| − μ)/2), fails where −μ > |J| (55 of the 300 points).
4. **E₋ = W₋** (ψ = 1), *numerical evidence*, quadrature converged to 4×10⁻⁵ or better. Selected rows:

   | σR | v | Eulerian ∫μ dV | E₋ (minimal split) | E₋ / (v R² σ) |
   |---|---|---|---|---|
   | 8 | 0.1 | −2.233×10⁻³ | 8.371×10⁻² | 0.1046 |
   | 8 | 1 | −0.2233 | 0.9432 | 0.1179 |
   | 16 | 0.1 | −4.450×10⁻³ | 0.1617 | 0.1011 |
   | 16 | 1 | −0.4450 | 1.832 | 0.1145 |
   | 16 | 5 | −11.13 | 15.45 | 0.1931 |

   The Eulerian negative energy grows like v², approaching −v²R²σ/36 for thin walls (a hand estimate the table
   approaches). The Eulerian integral itself is [PF97] Eq. 26, −(v²/12)∫r²f′² dr. E₋ grows like v at small v, driven
   by the momentum density of the wall, approaching about π v R²σ/32 ≈ 0.098 v R²σ (a hand estimate from ½∫|J| with
   the f″ term only). At v = 0.1, σR = 8, the minimal negative sector is 37 times the Eulerian negative energy.

   The local violation at O(v) is prior work. The null energy along ± the direction of motion has an exact term
   linear in v ([LV04a] Eqs. 15–19), boosted observers see negative energy at O(v) ([LV04b] Eq. 8), the flux behind it
   is [SSV22]'s (item 1), and [Le26] proves |j| ∝ v exactly, with a null-witness deficit linear in |j|. [LV04a]
   integrate the signed O(v) term, find zero, and call the net violation O(v²). That is a signed volume integral, while
   E₋ integrates a pointwise minimum over splits, so the two do not conflict ([prior_work.md](prior_work.md) §3).
5. **The Witten spinor exists for κ_K < 1** (*proved*, by hand, confirmed by #59's review, [review.md](review.md)). On a flat slice
   ∇̃_i = ∂_i + ½K_ij γ^jγ⁰, with the sign fixed in "Conventions", and K_ij γ^iγ^j = tr K (*checked symbolically*),
   so the Witten operator is D̃ = D + ½(tr K)γ⁰ with D = γ^i∂_i. That only tr K enters is [FNS11] Eq. 2.1, in
   two-spinor form. D is an isometry from Ḣ¹ onto L² (Fourier). With ε = ε₀ + ψ, the equation is
   ψ = −D⁻¹(½ tr K γ⁰(ε₀ + ψ)). Hölder, Kato's inequality and the sharp Sobolev inequality ‖h‖²₆ ≤ ‖∇h‖²₂ / S,
   S = 3(π/2)^{4/3}, make the map a contraction on Ḣ¹ with constant κ_K = ‖tr K‖_{L³}/(2√S), given tr K ∈ L² ∩ L³.
   This is the standard smallness argument that [FL22] recall for zero modes of γ·(−i∇ − A), with the diamagnetic
   inequality in the place of Kato's. So κ_K < 1 gives a unique Witten spinor in the class ε − ε₀ ∈ Ḣ¹, and D̃ has
   no kernel in Ḣ¹. Only tr K enters. For tr K = 0 (a maximal slice, as in Natário's zero-expansion drive [Nat02])
   ε = ε₀ is exact, as [FNS11] (Secs. 1 and 2.1) state for maximal, intrinsically flat slices. For the tanh profile,
   ‖tr K‖_{L³} ∝ v, and κ_K = 1 at v_c = 3.745, 2.471, 1.562 and 0.985 for σR = 2, 4, 8 and 16 (thin walls:
   v_c ≈ 2√S (15/(2πR²σ²))^{1/3}). The contraction says nothing beyond v_c. #63 ([frank_loss.md](frank_loss.md))
   extends existence and uniqueness with [FL22] to κ_K < √3, and beyond it to every real v outside a discrete set,
   which is empty when tr K has one sign (*proved*, by hand, given [FL22], and confirmed step by step by #68's
   independent review, [review_fl.md](review_fl.md)). GMRES solves converge to 1e-13 at up to 3.47 v_c for σR = 2
   and 4 (*numerical evidence*).
6. **The identity at E = P = 0** (*numerical evidence*). In these conventions, with ξ^j = ⟨ε, γ⁰γ^jε⟩, it reads
   (1/4π)∫|∇̃ε|² dV + ∫(μ|ε|² − J·ξ) dV = E − P·ξ₀ = 0. Its boundary terms need more decay than ε − ε₀ ∈ Ḣ¹ gives in
   general. Here the source of ψ decays exponentially, so ψ = O(r⁻²) and they vanish ([review.md](review.md) gap 2).
   Picard iteration of the FFT solve converges in 8–10 steps. The identity holds to 4×10⁻¹⁴, 4×10⁻¹⁴ and 1×10⁻¹²
   (relative to the gradient term G) at N = 128 for (σR, v) = (2, 0.5), (2, 1) and (4, 0.5), with spectral
   convergence in N (10⁻⁶, 10⁻¹⁰, 10⁻¹⁴ at N = 64, 96, 128). With the minimal split, ∫T₋(ξ_ε, n) = G + ∫T₊(ξ_ε, n):
   0.1635 = 0.0161 + 0.1474 at (2, 0.5). |ε|² stays within 4 % of 1.

## What this does and does not show

- Alcubierre data have M = 0, so they cannot test a lower bound on M, only its hypotheses. The sector-wise DEC can
  be met (item 3) at a cost in E₋ that is linear in v. The spinor exists below v_c (item 5).
- The identity's lower bound on the negative sector, ∫T₋(ξ_ε, n) ≥ G, is O(v²) (G ≈ |∫μ dV|), while the minimal E₋ is
  O(v). The spinor bound is therefore weak here, by the same factor that separates the Eulerian energy from E₋.
- The split uses the slice-wise condition μ± ≥ |J±|, which is what the Witten argument needs. A split whose sectors
  satisfy the full DEC for all observers would cost at least as much.
- Flat slices have E = 0 for every K, so they cannot carry the bound M ≥ −E₋/(1 − κ) of #13 with K ≠ 0. The existence
  argument of item 5 avoids Kato's inequality for ∇̃, the obstruction #13 names, only because the metric is flat.
- Nothing here claims novelty. [prior_work.md](prior_work.md) (#60) compares each result with prior work, and the
  items above cite that work where they use it.

## Conventions

γ⁰ = i diag(1, 1, −1, −1) (anti-Hermitian, γ⁰² = −1) and γ^j = [[0, σ_j], [σ_j, 0]] (Hermitian, γ^j² = 1), so
{γ^a, γ^b} = 2η^ab with η = diag(−1, 1, 1, 1). The positive inner product is ⟨a, b⟩ = a†b. K_ij = −g(∇_{e_i} n, e_j),
which on this slice is +½(∂_iβ_j + ∂_jβ_i), and J_i = −n^a T_ai. With these, the 4D route gives the J of item 1.

**The sign of ∇̃.** ∇̃ takes the K term with the + sign, ∇̃_i = ∂_i + ½K_ij γ^jγ⁰. It is the spacetime spin connection
¼ω_ab γ^aγ^b, restricted to the slice, after γ⁰ ↦ −γ⁰. That is, it is written in the gammas (−γ⁰, γ^j), which obey
the same Clifford relations, with ω_ab(X) = g(e_a, ∇_X e_b), e_0 = n and e_j = ∂_j. In the gammas above, the same
connection restricts to ∂_i − ½K_ij γ^jγ⁰ ([review.md](review.md) §3). [conventions.py](conventions.py) checks both
statements for lapse 1 and any shift β(t, x, y, z) (*checked symbolically*, output
[conventions_output.txt](conventions_output.txt)), and its control `--mutate=sign`, the other sign in ∇̃, fails.
With the other sign throughout, D̃ = D − ½(tr K)γ⁰ and item 6 holds with +J·ξ ([review.md](review.md) §4).

**The causal vector.** With ξ^j = ⟨ε, γ⁰γ^jε⟩, X = |ε|² n + ξ^j e_j is future causal, |ξ| ≤ |ε|², because the
γ⁰γ^j are Hermitian and anticommute. Then T(X, n) = μ|ε|² − J·ξ, which is item 6's integrand. Its time component
|ε|² is not ⟨ε, γ⁰γ⁰ε⟩, which is −|ε|² in these gammas. Flipping the sign of the J·ξ term (`--mutate=Jsign`) breaks
the identity, at 1.8×10⁻³ against 10⁻¹². phase4/RESULTS.md §1 uses these conventions.

## References

- [Alc94] M. Alcubierre, "The warp drive: hyper-fast travel within general relativity", Class. Quantum Grav. 11, L73–L77
  (1994), doi:10.1088/0264-9381/11/5/001, arXiv:gr-qc/0009013 (verified on arXiv).

The entries below are those of [`../literature.md`](../literature.md), "Warp drives", whose journal data were checked
against INSPIRE and arXiv identifiers against arXiv for #60.

- [BBV26] Barzegar, Buchert, Vigneron, arXiv:2602.16495 (2026, preprint).
- [FL22] Frank, Loss, arXiv:2201.03610 (2022); JEMS 28, 305 (2026) *(unverified)*.
- [FNS11] Frauendiener, Nester, Szabados, arXiv:1105.5008, CQG 28, 185004 (2011).
- [Le26] Le, arXiv:2602.18023, CQG 43, 215013 (2026).
- [LV04a] Lobo, Visser, arXiv:gr-qc/0406083, CQG 21, 5871 (2004).
- [LV04b] Lobo, Visser, arXiv:gr-qc/0412065 (2004, proceedings).
- [Nat02] Natário, arXiv:gr-qc/0110086, CQG 19, 1157 (2002).
- [PF97] Pfenning, Ford, arXiv:gr-qc/9702026, CQG 14, 1743 (1997).
- [SSV22] Santiago, Schuster, Visser, arXiv:2105.03079, PRD 105, 064038 (2022).

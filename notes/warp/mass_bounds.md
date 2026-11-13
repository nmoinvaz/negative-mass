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
   J_x = −v f′/(8π r_s) ≠ 0, so neither T nor −T satisfies μ ≥ |J|.
2. **E = P = 0** (*numerical evidence*, surface integrals). E vanishes identically on the flat slice. P is 1.8×10⁻³ at
   radius 2, 8×10⁻¹⁰ at 4 and 0 at 8 (σR = 4, v = 1). The same code returns m ψ(R)³ for Schwarzschild and P for
   Bowen–York data. So these data cannot violate any lower bound on M. They test the bounds' hypotheses only.
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
   approaches). E₋ grows like v at small v, driven by the momentum density of the wall, approaching about
   π v R²σ/32 ≈ 0.098 v R²σ (a hand estimate from ½∫|J| with the f″ term only). At v = 0.1, σR = 8, the minimal
   negative sector is 37 times the Eulerian negative energy.
5. **The Witten spinor exists for κ_K < 1** (*proved*, by hand, not yet independently reviewed). On a flat slice
   ∇̃_i = ∂_i + ½K_ij γ^jγ⁰ and K_ij γ^iγ^j = tr K (*checked symbolically*), so the Witten operator is
   D̃ = D + ½(tr K)γ⁰ with D = γ^i∂_i. D is an isometry from Ḣ¹ onto L² (Fourier). With ε = ε₀ + ψ, the equation
   is ψ = −D⁻¹(½ tr K γ⁰(ε₀ + ψ)). Hölder, Kato's inequality and the sharp Sobolev inequality ‖h‖²₆ ≤ ‖∇h‖²₂ / S,
   S = 3(π/2)^{4/3}, make the map a contraction on Ḣ¹ with constant κ_K = ‖tr K‖_{L³}/(2√S), given tr K ∈ L² ∩ L³.
   So κ_K < 1 gives a unique Witten spinor, and D̃ has no decaying kernel. Only tr K enters. For tr K = 0
   (a maximal slice, as in Natário's warp drive) ε = ε₀ is exact. For the tanh profile, ‖tr K‖_{L³} ∝ v, and
   κ_K = 1 at v_c = 3.745, 2.471, 1.562 and 0.985 for σR = 2, 4, 8 and 16 (thin walls: v_c ≈ 2√S (15/(2πR²σ²))^{1/3}).
   Beyond v_c this says nothing either way.
6. **The identity at E = P = 0** (*numerical evidence*). In these conventions, with ξ^j = ⟨ε, γ⁰γ^jε⟩, it reads
   (1/4π)∫|∇̃ε|² dV + ∫(μ|ε|² − J·ξ) dV = E − P·ξ₀ = 0. Picard iteration of the FFT solve converges in 8–10 steps. The
   identity holds to 4×10⁻¹⁴, 4×10⁻¹⁴ and 1×10⁻¹² (relative to the gradient term G) at N = 128 for (σR, v) = (2, 0.5),
   (2, 1) and (4, 0.5), with spectral convergence in N (10⁻⁶, 10⁻¹⁰, 10⁻¹⁴ at N = 64, 96, 128). With the minimal
   split, ∫T₋(ξ_ε, n) = G + ∫T₊(ξ_ε, n): 0.1635 = 0.0161 + 0.1474 at (2, 0.5). |ε|² stays within 4 % of 1.

## What this does and does not show

- Alcubierre data have M = 0, so they cannot test a lower bound on M, only its hypotheses. The sector-wise DEC can
  be met (item 3) at a cost in E₋ that is linear in v. The spinor exists below v_c (item 5).
- The identity's lower bound on the negative sector, ∫T₋(ξ_ε, n) ≥ G, is O(v²) (G ≈ |∫μ dV|), while the minimal E₋ is
  O(v). The spinor bound is therefore weak here, by the same factor that separates the Eulerian energy from E₋.
- The split uses the slice-wise condition μ± ≥ |J±|, which is what the Witten argument needs. A split whose sectors
  satisfy the full DEC for all observers would cost at least as much.
- Flat slices have E = 0 for every K, so they cannot carry the bound M ≥ −E₋/(1 − κ) of #13 with K ≠ 0. The existence
  argument of item 5 avoids Kato's inequality for ∇̃, the obstruction #13 names, only because the metric is flat.
- Nothing here claims novelty. [prior_work.md](prior_work.md) (#60) compares each result with prior work. It finds the
  flat-slice form of the Witten operator and the Sobolev–Hölder argument in earlier papers, and a sharper constant
  (Frank–Loss) that may raise every v_c by √3, not yet checked.

## Conventions

γ⁰ = i diag(1, 1, −1, −1) (anti-Hermitian, γ⁰² = −1) and γ^j = [[0, σ_j], [σ_j, 0]] (Hermitian, γ^j² = 1). The
positive inner product is ⟨a, b⟩ = a†b. With K_ij = +½(∂_iβ_j + ∂_jβ_i), the 4D route gives the J of item 1. Flipping the
sign of the J·ξ term (`--mutate=Jsign`) breaks the identity, at 1.8×10⁻³ against 10⁻¹². phase4/RESULTS.md §1 writes
T(ξ_ε, n) = μ|ε|² + J·ξ, which is this convention with ξ → −ξ.

## References

- [Alc94] M. Alcubierre, "The warp drive: hyper-fast travel within general relativity", Class. Quantum Grav. 11, L73–L77
  (1994), doi:10.1088/0264-9381/11/5/001, arXiv:gr-qc/0009013 (verified on arXiv).

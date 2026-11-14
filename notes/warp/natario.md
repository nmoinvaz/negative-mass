# Natário's zero-expansion drive: the sources, E₋ against v, and the Witten identity

November 2026, issue [#66](https://github.com/nmoinvaz/negative-mass/issues/66), after #58
([mass_bounds.md](mass_bounds.md)). Script [`natario.py`](natario.py), its stdout
[`natario_output.txt`](natario_output.txt), run from the repository root as `python notes/warp/natario.py` (5 min 16 s,
with the mutation runs four at a time, peak 2.7 GiB resident in one process, exits 1 if a check fails). It loads
[`alcubierre_bounds.py`](alcubierre_bounds.py) unmodified, without running it, for its profile, radial nodes, ADM
surface integrals and Witten solver.

Status words: **proved** (by hand), **checked symbolically** (sympy, exact, with a control that fails), **numerical
evidence**. Hand estimates are marked as such and compared with the numbers.

## Answer

The linear growth of E₋ with v survives tr K = 0, so it does not come from the expansion. For the tanh profile with
n = f/2, the local exponent d ln E₋ / d ln v is 1.00 for v ≤ 0.01 at every wall width tried, σR = 2 to 32 (**numerical
evidence**). The O(v) term is ½∫|J| dV (item 3), and J = −∇²X/16π does not vanish when ∇·X = 0.

The coefficient is larger than Alcubierre's, by a factor that grows like σR. At v = 0.001 the ratio of E₋ to #58's is
2.0, 3.3, 6.2, 12.3 and 24.7 for σR = 2, 4, 8, 16 and 32. For thin walls E₋/v → πR³σ²/(24√3), against #58's πR²σ/32, a
hand estimate that the table approaches to 0.1 % at σR = 128. The two thin-wall estimates cross at v ≈ 14/(σR), above
which E₋ tends to the Eulerian negative energy, ∝ v², which is −v²R⁴σ³/180 for thin walls (also a hand estimate, to
0.1 %).

The Witten spinor is ε₀ exactly, for every v and wall width, so the existence limit v_c of #58 item 5 does not arise.
With ε = ε₀ the identity is the integrated Hamiltonian constraint plus ∫J dV = P, and it holds at rounding level
pointwise and to quadrature error integrated (**proved**, then **numerical evidence**).

Scope. One profile family (n = f/2 with #58's tanh f, R = 1), constant v, axisymmetric about the direction of motion.

## Setup

[Nat02] writes the drive as ds² = −dt² + Σ(dxⁱ − Xⁱ dt)² with a divergence-free X. Its Sec. 2 takes X = −v d(f r²
sin²θ dφ) in the frame comoving with the bubble, which is X = −2v f cosθ e_r + v(2f + r f′) sinθ e_θ, with f = 0
inside and 1/2 outside. In the frame at rest outside the bubble, X is larger by v e_x. With n = 1/2 − f it is

    X = v ∇×(n(r) e_x × r⃗) = v [2n e_x + n′ (r e_x − x r̂)],    X^r = 2v n cosθ,   X^θ = −v(2n + r n′) sinθ,

with n = 1/2 inside the bubble and 0 outside. The t = 0 slice is flat with lapse 1, shift β = −X and, in the sign of #58,
K_ij = −½(∂_iX_j + ∂_jX_i), which is minus [Nat02]'s K. So ρ is the same and J changes sign, and the 4D route fixes the
sign of J without reference to a convention. The comoving frame has the same K, μ and J, because only ∂X enters.

Here n = f/2 with f the tanh profile of #58, R = 1, so that (σR, v) match #58's table. Then X = v f e_x + (v f′/2)(r
e_x − x r̂). Alcubierre's X = v f e_x [Alc94] lacks the term (v f′/2)(r e_x − x r̂), which cancels its divergence.

## Results

1. **The sources** (*checked symbolically*, by the 3+1 constraints and by G_ab of the 4-metric, which agree exactly, and
   by [Nat02]'s ρ). ∇·X = 0 and tr K = 0 identically, so μ = −|K|²/16π ≤ 0 with equality only where K = 0 ([Nat02]
   Sec. 1). In spherical components about the x-axis, K_rr = −2K_θθ = −2K_φφ = 2v n′ cosθ and
   K_rθ = −(v/2)(2n′ + r n″) sinθ, which are [Nat02] Sec. 2's components with n′ = −f′_N. Then

       μ = −(v²/16π) [6 n′² cos²θ + ½ (2n′ + r n″)² sin²θ]
       J = (v / 16π r) [−2w cosθ e_r + (r w)′ sinθ e_θ],    w = 4n′ + r n″,   (r w)′ = 4n′ + 6 r n″ + r² n‴

   and J = (1/16π)∇×∇×X = −∇²X/16π, which is the flux of [SSV22] Eqs. 4.7–4.8 with w = −β. With n = f/2 these are μ =
   −(v²/32π)[3f′² cos²θ + (f′ + r f″/2)² sin²θ], w = 2f′ + r f″/2 and (r w)′ = 2f′ + 3 r f″ + r² f‴/2. [Nat02] gives ρ
   and no momentum density. The 3+1 route differentiates X, the 4D route computes the Einstein tensor of the metric with
   X(x − v t, y, z) and projects with n^a = (1, Xⁱ), and each reproduces the closed forms exactly in the jets n, …, n‴.
   The closed μ equals [Nat02]'s ρ. A doubled shift in the 3+1 route fails both checks (`--mutate=beta2`). On the axis, μ
   = −3v²f′²/32π ≠ 0, where Alcubierre's μ vanishes, and |J| = v|w|/(8π r).

2. **E = P = 0** (*proved*, then *numerical evidence*). E's integrand vanishes since g_ij = δ_ij. P is the surface
   integral of K, which decays exponentially. At σR = 4, v = 1 the same code that returns m ψ(R)³ for Schwarzschild
   and P for Bowen–York data gives P_x = −1.1×10⁻² at radius 2, −1.1×10⁻⁸ at 4 and 0 at 8. By the momentum constraint
   ∫J dV = P, and the quadrature gives |∫J_x dV| / ∫|J| dV = 3×10⁻¹⁴. So these data, like Alcubierre's, have M = 0 and
   test the bounds' hypotheses only.

3. **The minimal split** (the formula is #58 item 3, *proved*; rechecked here). μ₋ = max(0, −μ, (|J| − μ)/2) agrees with a
   restarted Nelder–Mead minimisation to 2×10⁻¹⁶ at 150 random (μ, J) and to 9×10⁻¹⁶ at 150 points of the Natário wall
   (σR = 4, 8, 16, v = 0.01 to 10, 8 of them with −μ > |J|), and the named split is admissible at every point. Since μ ≤ 0
   here, μ₋ = max(|μ|, (|J| + |μ|)/2), so (*proved*, by hand)

       E₋ = ½ (∫|J| dV + |E_Eul|) + ½ ∫(|μ| − |J|)₊ dV,    E_Eul = ∫μ dV = −(v²/8) ∫ r² [4 n′² + (2/3)(2n′ + r n″)²] dr.

   So E₋ ≥ |E_Eul| and E₋ ≥ ½∫|J|, and E₋ = ½∫|J| + O(v²) at small v. The Eulerian integral is [Nat02]'s ρ integrated.

4. **E₋ against v and wall width** (*numerical evidence*). Selected rows of the 35 in `natario_output.txt`, with #58's
   E₋ at the same (σR, v) from the same quadrature, which reproduces the 20 rows of `alcubierre_bounds_output.txt` to
   3×10⁻⁶:

   | σR | v | Eulerian ∫μ dV | E₋ | E₋ of #58 | ratio | E₋ / (v σ²) | E₋ / (−∫μ dV) |
   |---|---|---|---|---|---|---|---|
   | 2 | 0.001 | −1.752×10⁻⁷ | 6.094×10⁻⁴ | 3.107×10⁻⁴ | 1.96 | 0.1524 | 3479 |
   | 2 | 1 | −0.1752 | 0.6974 | 0.3434 | 2.03 | 0.1743 | 3.98 |
   | 8 | 0.1 | −3.288×10⁻² | 0.5302 | 8.371×10⁻² | 6.33 | 0.0828 | 16.1 |
   | 8 | 1 | −3.288 | 7.034 | 0.9432 | 7.46 | 0.1099 | 2.14 |
   | 16 | 0.1 | −0.2364 | 2.088 | 0.1617 | 12.9 | 0.0816 | 8.83 |
   | 16 | 1 | −23.64 | 35.34 | 1.832 | 19.3 | 0.1381 | 1.495 |
   | 16 | 5 | −590.99 | 614.81 | 15.45 | 39.8 | 0.4803 | 1.040 |
   | 32 | 0.001 | −1.838×10⁻⁴ | 7.783×10⁻² | 3.156×10⁻³ | 24.7 | 0.0760 | 423 |
   | 32 | 1 | −183.8 | 218.1 | 3.634 | 60.0 | 0.2130 | 1.187 |
   | 32 | 5 | −4595 | 4663 | 30.77 | 151.5 | 0.9106 | 1.015 |

   The local exponent p = d ln E₋ / d ln v between successive v is

   | σR | 0.001–0.01 | 0.01–0.1 | 0.1–0.5 | 0.5–1 | 1–2 | 2–5 |
   |---|---|---|---|---|---|---|
   | 2 | 1.001 | 1.006 | 1.034 | 1.094 | 1.178 | 1.414 |
   | 4 | 1.001 | 1.008 | 1.046 | 1.130 | 1.257 | 1.536 |
   | 8 | 1.001 | 1.012 | 1.078 | 1.227 | 1.425 | 1.710 |
   | 16 | 1.002 | 1.023 | 1.150 | 1.410 | 1.655 | 1.865 |
   | 32 | 1.005 | 1.047 | 1.287 | 1.650 | 1.843 | 1.948 |

   E₋ is linear in v for v ≲ 0.01 and tends to ∝ v² for large v and thin walls, where it also tends to |E_Eul| (1.015 at
   σR = 32, v = 5). The deficit ∫(|J| − μ)₊ = ∫|J| + |E_Eul| is in the output. Doubling the nodes, (80, 24, 96) to
   (160, 32, 192) in (panels, nodes per panel, angular nodes), changes E_Eul, E₋ and the deficit by at most 1.3×10⁻⁵.

5. **Thin-wall hand estimates** (compared with *numerical evidence*). For f ≈ (1 − tanh σ(r − R))/2 near the wall, the
   n″² term of E_Eul gives −(v²/12) ∫ r⁴ n″² dr = −v²R⁴σ³/180, using ∫f″² dr = 4σ³/15. At small v, |J| ≈ (v/16π) r
   |n‴| sinθ, so E₋ ≈ ½∫|J| dV = (π v/32) ∫ r³ |n‴| dr = π v R³σ²/(24√3), using ∫|f‴| dr = 8σ²/(3√3). At σR = 128 the
   numbers are E_Eul / (−v²σ³/180) = 1.0006 and E₋ / (π v σ²/(24√3)) = 1.0008 at v = 10⁻⁴. The first scales like the
   sharp bound in the abstract of [BVA26], v²R⁴/Δ³ with Δ = 1/σ (abstract only). The two estimates cross at v =
   180·π/(24√3)/(σR) ≈ 14/(σR), which is consistent with p ≈ 1.5 near v ≈ 0.5 to 1 at σR = 16 and 32. Against #58's
   πvR²σ/32 the E₋ ratio is 4σR/(3√3), which gives 24.6 at σR = 32 against 24.7 in the table.

6. **The Witten identity** (*proved*, then *numerical evidence*). For tr K = 0, D̃ = D, so ε = ε₀ solves the Witten
   equation, and it is the only solution with ε − ε₀ ∈ Ḣ¹ because D is an isometry from Ḣ¹ onto L²
   ([mass_bounds.md](mass_bounds.md) item 5). Then κ_K = 0 for every v and wall width. With ∇̃_iε₀ = ½K_ij γ^jγ⁰ε₀,
   Σ_i |K_ij γ^jγ⁰ε|² = |K|²|ε|² for every symmetric K and spinor ε (*checked symbolically*), so (1/4π)∫|∇̃ε₀|² dV =
   (1/16π)∫|K|² dV = −∫μ dV, and the identity of [mass_bounds.md](mass_bounds.md) item 6 becomes −ξ₀·∫J dV = −P·ξ₀
   with E = 0. It holds because ∫J dV = P, and it needs no solve.

   The solve is run anyway. The ported solver of #58, with K as an argument, reproduces #58's own solve on Alcubierre
   data (σR = 2, v = 0.5, N = 64) with relative difference 0 and 8 iterations in both. On Natário data it takes one
   iteration and returns ε₀, with |ε − ε₀| ≤ 9×10⁻¹⁹ |K| and |D̃ε₀| ≤ 1.4×10⁻¹⁷ |K|, both rms, at rounding level. Over
   (σR, v) = (2, 0.5), (2, 1), (4, 0.5), (4, 2) and (8, 0.1), with ε₀ = e₁ (ξ₀ = 0) and the null spinors with ξ₀ = ∓e_x:

   - Pointwise, |∇̃ε₀|²/4π + μ|ε₀|² = 0 at every node to 8×10⁻¹⁶ of max |μ|, with K from the 3+1 route and μ from the
     closed form, which share no derivation.
   - Integrated, |G + ∫T(ξ, n)| is at most 3.2×10⁻¹⁶ of G for e₁, and for e± it is ∓ξ₀·∫J on the grid, spectral in N. At
     (2, 0.5) it is 1.70×10⁻⁶, 1.65×10⁻¹⁰ and 3.15×10⁻¹¹ of the sum of its terms G + |ξ₀|∫|J| at N = 64, 96 and 128.
     The worst at the finest grids is 7.8×10⁻¹⁰.
   - With the minimal split, ∫T₋(ξ₀, n) = G + ∫T₊(ξ₀, n), T₊ ≥ 0 pointwise: 4.298 = 2.325 + 1.973 for e₁, and
     4.185 = 2.325 + 1.860 and 4.410 = 2.325 + 2.085 for e₊ and e₋, at (4, 2). At (2, 0.5) it is 0.3265 = 0.0438 + 0.2828.
   - The grid E₋ agrees with the (r, cosθ) quadrature to 3.6×10⁻⁴ at worst, from the kinks of the max, and the grid E_Eul
     to 1.9×10⁻⁸, which is also how well G = |E_Eul| holds at the finest grids.

   The spinor bound ∫T₋(ξ₀, n) ≥ G = |E_Eul| is O(v²), against E₋ = O(v), so it is weak at small v by E₋/|E_Eul|
   (35 at σR = 2, v = 0.1). It is nearly sharp where E₋ → |E_Eul|, with T₊ → 0, at large v and thin walls.

## What this does and does not show

- E₋ ∝ v does not need an expansion. That the O(v) term comes from J is shown only by the formula of item 3. Nothing
  computed here separates it from other effects of the profile.
- J = (1/16π)∇×∇×w, so J = 0 for every zero-vorticity flat-slice drive, and item 3 would then give E₋ = ∫(−μ)₊, which
  is O(v²). That is a hand consequence and is not computed here, so whether the O(v) cost needs vorticity is open.
- Natário data have M = 0 like Alcubierre's, so they test the bounds' hypotheses only. The difference from #58 item 5
  is the spinor, which exists for all v here, and for Alcubierre only below v_c = 3.7, 2.5, 1.6 and 0.99 at σR = 2, 4,
  8 and 16.
- At tr K = 0 the identity is pointwise algebra plus ∫J = P, so the solve is a control of the solver, not an
  independent test of the identity. The pointwise check, the symbolic gamma identity and the grid E₋ do the testing.
- The integrated identity does not detect every error. `--mutate=gamdup` adds the term (K²)_xz, which is odd in z and
  integrates to zero, so only the pointwise check sees it.
- Nothing here claims novelty. ρ is [Nat02]'s, the flux is [SSV22]'s, and [Nat02] Secs. 1 and 2 were read in the arXiv
  v3 text for this issue. I have not compared E₋ or the split for Natário drives with prior work beyond
  [prior_work.md](prior_work.md).

## Controls and mutations

Every check has a mutation, and `--part=mutations` runs each in a subprocess, requires its named checks to fail, requires
the unmutated quick runs to pass, and requires every check to be named by some mutation.

| `--mutate=` | what it does | check that fails |
|---|---|---|
| `div` | halves X^y and X^z, so ∇·X ≠ 0 | divergence, [Nat02] field, 3+1 and 4D μ and J, and in the solve ε₀, the pointwise identity and the Witten identity |
| `beta2` | doubled shift in the 3+1 route | 3+1 μ and J |
| `Jsign` | opposite sign of J in the 3+1 route | 3+1 J |
| `nusign` | spatial part of the normal flipped in the 4D route | 4D μ and J |
| `natX`, `natrho` | sign of the X^θ term, 3 → 2 in ρ, in the transcription of [Nat02] | field, ρ |
| `rw1`, `mu6` | drops r²n‴ in (r w)′, 6 → 5 in μ, in the numpy forms | numeric against symbolic, and for `mu6` grid E_Eul |
| `prof3` | f‴ → −f‴ | profile derivatives |
| `gamdup` | γ³ replaced by γ¹ | gamma algebra, pointwise identity |
| `adm_tail`, `adm_m`, `adm_p` | Bowen–York tail on K, Schwarzschild m × 1.1, Bowen–York P × 1.1 | P(R = 8), the two controls |
| `split_old` | #58's first split formula | the three split checks, T₊ ≥ 0 |
| `quad` | angular weights × (1 + 0.01u²) | Gaussian control, #58's table, ∫J = 0 |
| `coarse`, `thin` | reported nodes (4, 4, 8), π/(24√3) → π/24 | convergence, thin-wall estimates |
| `port`, `norm`, `musign` | sign of the tr K term, 4π → 2π, μ → −μ | solver control, pointwise and Witten identity |
| `Jgrid` | J × (1 + 0.1 r) on the grid | grid ∫J, grid E₋ |

## Reported after the first runs

1. **Convergence protocol.** #58 compares nodes (40, 16, 64) with the reported (80, 24, 96). Here that comparison differs
   by up to 1.7×10⁻⁴ (σR ≥ 16, from the radial panels), above the 1×10⁻⁴ criterion. The criterion was kept and the test
   changed to doubling the reported nodes. The comparison is still printed in the output.
2. **Identity normalisation.** Normalised by G, the integrated residual was 2.5×10⁻⁸ at (8, 0.1), N = 128, above the 10⁻⁸
   criterion. It is ξ₀·∫J on the grid, 8×10⁻¹⁰ of ∫|J|, and ∫|J| is about 30 G at small v. The criterion now normalises by
   G + |ξ₀|∫|J|, and the pointwise check was added. Both normalisations are in the output.
3. **Mutations.** A `lapse` mutation changed nothing, since a constant lapse N rescales t, and with ν = (1, X) both
   projections are unchanged, so it was replaced by `nusign`. `prof3` and `gamdup` fire only the checks listed above.
4. **Nelder–Mead** stalled on the kink of the max, near the 10⁻⁶ criterion on one set of random points, so it restarts
   from its best point.

## References

- [Alc94] M. Alcubierre, "The warp drive: hyper-fast travel within general relativity", Class. Quantum Grav. 11, L73–L77
  (1994), arXiv:gr-qc/0009013.
- [Nat02] J. Natário, "Warp drive with zero expansion", Class. Quantum Grav. 19, 1157–1166 (2002), arXiv:gr-qc/0110086
  (v3 text read, Secs. 1–2).

The entries below are those of [`../literature.md`](../literature.md), "Warp drives" and "Energy conditions", checked
for #60.

- [BVA26] Bolívar, Vasilev, Abellán, arXiv:2609.36211 (2026, preprint, abstract read).
- [SSV22] Santiago, Schuster, Visser, arXiv:2105.03079, PRD 105, 064038 (2022).

# The tanh profile's Frank–Loss spectrum, without a periodic box

Issue #75. Script [tanh_spectrum.py](tanh_spectrum.py), output [tanh_spectrum_output.txt](tanh_spectrum_output.txt), raw
output in [negative-mass-raw](https://github.com/nmoinvaz/negative-mass-raw) folder `I75-dh170` (commit 93a55e2). Run from
the repository root, the commands that wrote the output are at its top. `python notes/warp/tanh_spectrum.py --controls`
runs every `--mutate=` control (`g0block`, `ladder`, `dterm`, `sobolev`, `half`, `lscale`, `rot`, `trunc`) and exits 1
unless each fires. The reviewed claims are those of [frank_loss.md](frank_loss.md) §9, "ρ(T₁)v_c = 0.40938" and "E is
empty for the tanh profile", which rest on #58's FFT box.

## Answer

The spectral radius is confirmed without a box, to a relative difference of 3.8×10⁻⁴ at σR = 2 and 7.0×10⁻⁴ at σR = 4.
Status: *numerical evidence*.

| σR | ρ(T₁)v_c here | two discretisations differ by | #63, N = 128 | relative difference | v_c |
|---|---|---|---|---|---|
| 2 | 0.409556 | 5×10⁻¹² | 0.4094 (L = 8) | 3.8×10⁻⁴ | 3.744900 |
| 4 | 0.367657 | 5×10⁻¹⁰ | 0.3674 (L = 5) | 7.0×10⁻⁴ | 2.470971 |

The planted weight τ = 6/(1 + r²) returns ρ = 1 to 7×10⁻¹⁵ before any real data is computed, and its translated and
dilated image, which is axisymmetric and not spherical, returns the same spectrum to 1.6×10⁻¹¹. The leading eigenvalues
of T₁ are the same complex quadruple #63 found, λv_c = ±0.02730 ± 0.40865i at σR = 2 and ±0.01590 ± 0.36731i at
σR = 4, each fourfold in T₁, so sixteen in the cluster as in #63. The next quadruples, ±0.00483 ± 0.29896i at σR = 2 and
±0.00278 ± 0.28852i at σR = 4, are #63's "±0.0048 ± 0.2989i and ±0.0028 ± 0.2885i".

No real eigenvalue of S⁻¹M, that is no real zero-mode speed, lies in a converged window. Every eigenvalue with
|μ|v_c ≥ 0.0830 (σR = 2) and ≥ 0.0969 (σR = 4) is converged in every J_z sector computed, m = ½ to 21/2, and the leading
eigenvalue of m = 21/2 is below that, so a real zero mode of D̃ has |v| > 12.05 v_c at σR = 2 and |v| > 10.32 v_c at
σR = 4. #63 reached 4 v_c with the grid operator. That E is empty for the tanh profile stays a *conjecture*.

The identity ∫τ₁|ψ|² = ∫τ₁|χ|² = 0 has nothing to evaluate on, since there is no real candidate, and for this profile it
would hold without one (§5).

## 1. What is computed

τ₁ = −f′(r) cos θ is axisymmetric about the x axis. With the reduction of [frank_loss.md](frank_loss.md) §2, rederived
here symbolically (`reduction`), D = diag(S, S) and B₁ = diag(−M, +M) in the coordinates (w + iu, w − iu) of the range
and (ψ, χ) of the domain, with S = σ·∇ and M = τ₁/2. So T₁ = −S⁻¹M on the ψ sector and +S⁻¹M on the χ sector, and

    spec(T₁) = Σ ∪ −Σ,   Σ = spec(S⁻¹M),   real zero mode of D̃ at speed v  ⇔  real μ = 1/v in Σ.

The cyclic relabelling (x, y, z) → (z, x, y) of the Pauli matrices makes the profile's axis the z axis, with no change to
σ·∇ or to the handedness. Then J_z = −i∂_φ + σ_z/2 commutes with S and with M, and each half-integer m = J_z is a
separate problem. The spectrum of m < 0 is that of m > 0 up to a sign, by a reflection through a plane containing the
axis, so only m > 0 is solved.

```
 σ·∇ψ = ½vτ₁ψ on R³
        │  J_z = m, one problem per m = ½, 3/2, …
        ▼
 spin-angular channels Ω_κ, κ = ±(m+½), …, ±(m+½+dK)     (σ·L Ω_κ = −(1+κ) Ω_κ)
        │  S(f Ω_κ) = (f′ + (1+κ)f/r) σ·r̂Ω_κ,   cos θ couples κ to its l ± 1 neighbours
        ▼
 radial functions f = y^l (1−y)² q(y),  y = r/(r+Lmap),  q Jacobi       (no box, no outer cut)
        │  Galerkin in the Ḣ¹ inner product ⟨Sa, Sb⟩
        ▼
 K x = v N x,   K_ij = ⟨Sφ_i, Sφ_j⟩,  N_ij = ⟨Sφ_i, Mφ_j⟩,  both real
        │  parity under x → −x splits the channels in two blocks, each a dense generalised eigenproblem
        ▼
 μ = 1/v, all eigenvalues, with eigenvectors for the leading ones
```

### Why the Galerkin eigenvalues are the right ones

S⁻¹M is compact on Ḣ¹ ([frank_loss.md](frank_loss.md) §6) and ⟨Sa, Sb⟩ = ⟨∇a, ∇b⟩ is the Ḣ¹ inner product, so
K⁻¹N is the compression P_h S⁻¹M of the compact operator to a subspace V_h ⊂ Ḣ¹. If V_h exhausts Ḣ¹, then P_h → I
strongly and P_h S⁻¹M → S⁻¹M in norm, and the eigenvalues of the compression converge to those of the operator with
multiplicities and with no spurious eigenvalues. That V_h exhausts Ḣ¹ is not proved here. The planted checks and the
convergence in every parameter (§4) are what stand in for it.

Every basis function lies in Ḣ¹ on the whole of R³: regular at 0 (f ~ r^l) and f ~ r⁻² at infinity, which is the free
decaying solution, so L² as well. There is no box and no outer truncation. The radial integrals are exact polynomials in
y for K and smooth in y for the weights.

### What the matrices' structure gives

- The matrices are real, because the Ω_κ are real at φ = 0. So eigenvalues are real or in conjugate pairs, and a real
  zero mode is a real eigenvalue of a real matrix.
- Under x → −x the weight is odd, so ψ(x) → ψ(−x) maps solutions at v to solutions at v, and the channels split by the
  parity of l into two blocks that do not couple (`pat` in the script, connected components of the angular
  coupling). Within one m the spectrum is symmetric under μ → −μ (diagnostic column of the tanh table, §4), since
  S⁻¹M has adjoint −M S⁻¹ and the matrices are real.
- For an eigenfunction of one parity class |ψ|² is even and τ₁ is odd, so ∫τ₁|ψ|² = 0 identically (§5).

## 2. What is shared with frank_loss.py

`frank_loss.py` and [review_fl.md](review_fl.md) were read, as the task directs, and no numerics, matrix, basis or solver
were taken from them. The script was written from the formulation in [frank_loss.md](frank_loss.md) §2 and §9.

| shared | where it came from | what happens to it here |
|---|---|---|
| the problem, D̃ = D + ½τγ⁰ with the gammas of [mass_bounds.md](mass_bounds.md) | mass_bounds.md "Conventions" | the reduction is rederived (`reduction`) |
| the profile f(r) = [tanh σ(r+R) − tanh σ(r−R)]/(2 tanh σR), R = 1, τ₁ = −f′cos θ | mass_bounds.md item 1 | rewritten, `v_c` recomputed by quadrature and compared to the printed values |
| S = 3(π/2)^{4/3} and v_c = 2√S/‖τ₁‖₃ | mass_bounds.md item 5 | the same closed form, not rederived |
| the planted weight 6/(1 + r²) and its exact levels μ = ±i/ℓ_n, ℓ_n = 1 + 2n/3 | [review_fl.md](review_fl.md) `planted_ode` | compared against, not used in the solver |
| the reporting helpers `check`, `peak_gb` and the control runner, and the names `reduction`, `v_c`, `planted` | frank_loss.py's layout | a pattern, not numerics |
| the claim ρ(T₁)v_c = 0.4094 and 0.3674, and the next quadruples | [frank_loss.md](frank_loss.md) §9 | compared against |
| the formulation, T₁ = D⁻¹B₁ compact on Ḣ¹ with the Ḣ¹ inner product | [frank_loss.md](frank_loss.md) §6–§7 | same operator, a different discretisation |

Not shared: the discretisation (a spin-angular and rational-radial Galerkin space, not a periodic FFT grid), the
eigensolver (dense generalised eigenproblems, not ARPACK on matrix-free operators), the spinor representation (2-spinors
in J_z sectors, not 4-spinors on a cube), and every line of code. What the two share beyond the table is the
mathematics: both assume the Ḣ¹ framework of §6, and both would fail together if the reduction or the sign conventions
of mass_bounds.md were wrong.

## 3. Checks and controls

Each check ran on a case whose answer is known before real data was computed. Planted truth runs first, and the tanh
part does not run when `planted`, `shift`, `angular` or a radial check has failed.

| check | what it tests | result |
|---|---|---|
| `reduction` | W′DW⁻¹ = diag(S, S), W′B₁W⁻¹ = diag(−M, M), symbolic | passes |
| `angular` | Ω_κ orthonormal, \|⟨Ω_i\|σ·r̂\|Ω_j⟩\| = δ_{i,−j} | 4×10⁻¹⁴ |
| `angular_u` | ⟨Ω_i\|cos θ\|Ω_j⟩ against the Gaunt formula | 4×10⁻¹⁴ |
| `angular_dirac` | σ·∇(f Ω_κ) against finite differences in Cartesian coordinates | 1×10⁻⁹ |
| `radial_op` | the radial operator in y against finite differences in r | 1×10⁻¹⁰ |
| `radial_int` | the y-integrals K, N, G against quadrature in r | 9×10⁻⁸ |
| `v_c` | v_c by quadrature against mass_bounds.md | 3.7×10⁻⁶ |
| `planted` | 6/(1 + r²), m = ½, 3/2, 5/2, levels n ≤ 4 with multiplicities, ρ = 1 | 7×10⁻¹⁵ |
| `shift` | the same weight translated by 0.5 and dilated by 1.2, levels n ≤ 4 | 1.6×10⁻¹¹ |
| `converged`, `conv` | two discretisations, and one parameter at a time | §4 |
| `ref63` | ρ(T₁)v_c against frank_loss.md, relative 10⁻³ | §4 |
| `real_identity`, `no_real_in_window` | §5 | §5 |

`planted` checks the radial machinery and the Ω basis for all j, with eigenvalue pairing level by level and the
multiplicity n + 3/2 − m per sign at fixed m > 0. `shift` is the check of the cos θ coupling: the translated and
dilated weight τ(x) = 6s/(1 + s²|x − a e_z|²) is a Möbius image of the planted one, so by the conformal covariance of
σ·∇ it has the same eigenvalues, and as a weight it needs every Legendre order P_l(cos θ) up to l = 26 in the angular
matrices that the tanh run uses for l = 1. The tanh weight itself is checked by `radial_int` and `v_c`.

Controls. Each mutation is run beside the same command without it, and counts as fired only if its target checks pass
without it and fail with it.

| mutation | what it breaks | checks that fail |
|---|---|---|
| `g0block` | lower block of γ⁰ flipped | `reduction` |
| `ladder` | sign of the σ·L ladder term | `angular`, `angular_dirac` |
| `dterm` | (1 + κ)f/r dropped in the radial operator | `radial_op`, `radial_int` |
| `sobolev` | 21/20 of the sharp S | `v_c` |
| `half` | ½ dropped in M | `planted`, and `ref63` on the tanh weight |
| `lscale` | l = 1 Legendre term of the shifted weight × 11/10 | `shift` |
| `rot` | N × i, so the planted eigenvalues are real | `real_identity` (and `planted`) |
| `trunc` | 3 angular levels and 8 radial functions | `converged` (and `ref63`) |

All eight fired (`--controls`, quick sizes, 21 s, nine runs because `half` has two).

## 4. The spectrum

Sizes. Discretisation A is dK = 28 angular levels above each sign's lowest, p = 80 radial functions at σR = 2 and 100
at σR = 4, Lmap = 2, 160 quadrature panels of 16 nodes. B has dK + 10, p + 10, Lmap = 3 and 240 panels. Eleven sectors,
m = ½ to 21/2, were solved at each σR. The eigenvalues of the two discretisations differ by a quantity that is an
estimator range, not an enclosure: Galerkin eigenvalues of a non-normal compact operator are not bounds.

The leading eigenvalues of S⁻¹M times v_c, by J_z sector. T₁ has them and their negatives. For m = ½ each of the four
values ±a ± ib occurs once in Σ at m and once at −m, and again with the opposite sign in T₁, so each is fourfold in T₁
and the cluster has sixteen, as in #63. "Converged above" is the largest |μ|v_c at which an eigenvalue of A has no
partner in B within 10⁻⁵ relative, so every eigenvalue above it is stable to that.

σR = 2, v_c = 3.744900, ρ(T₁)v_c = 0.40955601, √3 ρ v_c = 0.7094 (FL: ≤ 1).

| m | leading μv_c (a quadruple ±a ± ib) | \|μ\|v_c | next \|μ\|v_c | converged above |
|---|---|---|---|---|
| 0.5 | 0.02730 + 0.40865i | 0.409556 | 0.257195 | 0.0777 |
| 1.5 | 0.00483 + 0.29895i | 0.298994 | 0.207406 | 0.0803 |
| 2.5 | 0.00106 + 0.23298i | 0.232982 | 0.173753 | 0.0751 |
| 3.5 | 0.00026 + 0.19037i | 0.190371 | 0.149201 | 0.0711 |
| 4.5 | 0.00007 + 0.16080i | 0.160798 | 0.130580 | 0.0667 |
| 5.5 | 0.00002 + 0.13913i | 0.139127 | 0.116022 | 0.0688 |
| 6.5 | 0.00001 + 0.12258i | 0.122581 | 0.104349 | 0.0649 |
| 7.5 | 0.00000 + 0.10954i | 0.109541 | 0.094791 | 0.0611 |
| 8.5 | 0.00000 + 0.09900i | 0.099002 | 0.086825 | 0.0578 |
| 9.5 | 0.00000 + 0.09031i | 0.090309 | 0.080088 | 0.0595 |
| 10.5 | 0.00000 + 0.08302i | 0.083018 | 0.074316 | 0.0563 |

σR = 4, v_c = 2.470971, ρ(T₁)v_c = 0.36765732, √3 ρ v_c = 0.6368.

| m | leading μv_c (a quadruple ±a ± ib) | \|μ\|v_c | next \|μ\|v_c | converged above |
|---|---|---|---|---|
| 0.5 | 0.01590 + 0.36731i | 0.367657 | 0.250083 | 0.0969 |
| 1.5 | 0.00278 + 0.28852i | 0.288531 | 0.211765 | 0.0910 |
| 2.5 | 0.00060 + 0.23538i | 0.235378 | 0.182557 | 0.0858 |
| 3.5 | 0.00015 + 0.19848i | 0.198480 | 0.160023 | 0.0806 |
| 4.5 | 0.00004 + 0.17151i | 0.171512 | 0.142278 | 0.0840 |
| 5.5 | 0.00001 + 0.15097i | 0.150969 | 0.128000 | 0.0790 |
| 6.5 | 0.00000 + 0.13481i | 0.134809 | 0.116287 | 0.0745 |
| 7.5 | 0.00000 + 0.12177i | 0.121766 | 0.106515 | 0.0705 |
| 8.5 | 0.00000 + 0.11102i | 0.111020 | 0.098245 | 0.0727 |
| 9.5 | 0.00000 + 0.10201i | 0.102014 | 0.091157 | 0.0688 |
| 10.5 | 0.00000 + 0.09436i | 0.094358 | 0.085018 | 0.0653 |

The leading eigenvalue falls slowly with m, roughly like 1/m at large m (m|μ|v_c = 0.87 at m = 21/2, σR = 2). Its real
part falls faster, from 2.7×10⁻² at m = ½ to 7×10⁻⁵ at m = 9/2 and 4×10⁻⁷ at m = 17/2 (σR = 2). In each sector the
spectrum is symmetric under μ → −μ, to at most 1.7×10⁻⁶ of |μ| over the converged eigenvalues (the last column of the
output's tables).

Convergence. One parameter at a time from the base dK = 20, p = 80 (σR = 2) or 100 (σR = 4), Lmap = 2 and 160 panels,
m = ½, as `--part=conv`, the change in |μ|v_c.

| change | σR = 2, leading | σR = 2, fifth | σR = 4, leading | σR = 4, fifth |
|---|---|---|---|---|
| dK = 12 | +3×10⁻¹⁴ | −7×10⁻⁹ | −1×10⁻¹² | −1×10⁻⁷ |
| dK = 28 | +1×10⁻¹⁶ | +7×10⁻¹⁶ | +1×10⁻¹⁵ | +3×10⁻¹⁶ |
| p − 20 | +4×10⁻¹⁰ | −2×10⁻¹¹ | −1×10⁻⁹ | −4×10⁻⁹ |
| p + 20 | −6×10⁻¹² | −7×10⁻¹² | +5×10⁻¹² | +2×10⁻¹¹ |
| Lmap = 1 | −6×10⁻¹² | −7×10⁻¹² | +6×10⁻¹⁰ | +4×10⁻¹⁰ |
| Lmap = 4 | −5×10⁻¹² | −7×10⁻¹² | −3×10⁻¹⁰ | −1×10⁻⁹ |
| 80 panels | −4×10⁻¹⁶ | +1×10⁻¹⁵ | +7×10⁻¹⁶ | +4×10⁻¹⁶ |
| 320 panels | −1×10⁻¹⁶ | +2×10⁻¹⁶ | +6×10⁻¹⁶ | +1×10⁻¹⁶ |

The leading eigenvalue moves by at most 7×10⁻⁹ relative in any variant, over the two sectors m = ½ and 3/2 and both
σR (`conv`). The angular depth dK is what sets the converged window: dK = 12 moves the fifth eigenvalue by 7×10⁻⁹ to
1.7×10⁻⁶ over the four cases, dK = 28 by less than 10⁻¹⁴. The L² residual ‖Sψ − vMψ‖/‖Sψ‖ of the leading
eigenfunctions is 1×10⁻⁶ to 5×10⁻⁵ at σR = 2 and 2×10⁻⁴ to 1.4×10⁻³ at σR = 4, rising with m. The eigenvalue
converges far faster than the eigenfunction: A and B differ by 5×10⁻¹⁰ in ρ at σR = 4, where the residual is 2×10⁻⁴
and above.

Cost. The full run took 36 min of wall time with four workers at `nice -n 10` and one thread each, at 1.7 GiB (σR = 2)
and 2.6 GiB (σR = 4) of peak resident memory per worker, and the machine's load stayed at four.

## 5. Real eigenvalues and the identity

Σ is searched in every sector for real eigenvalues, the real zero-mode speeds v = 1/μ. A real eigenvalue of a real
matrix is stable under refinement, and a spurious one is not, so each is classified by the discretisation that does not
reproduce it.

| σR | sectors | eigenvalues in the converged windows (A / B) | real in the windows | real with \|μ\|v_c > 10⁻⁸ in A / in B, per sector | the largest \|μ\|v_c of them | real in A that B reproduces to 10⁻⁴ |
|---|---|---|---|---|---|---|
| 2 | 11 | 652 / 676 | 0 | 34 to 50 / 30 to 44 | 1.3×10⁻³ | 0 |
| 4 | 11 | 396 / 424 | 0 | 55 to 66 / 43 to 55 | 1.2×10⁻³ | 0 |

The window of a sector is \|μ\|v_c above its "converged above" value, 0.056 to 0.080 at σR = 2 and 0.065 to 0.097 at
σR = 4. There is no real eigenvalue in any of them, in A or in B. Below the windows there are 30 to 66 real eigenvalues
per sector, all at \|μ\|v_c ≤ 1.3×10⁻³, and none of them is reproduced by the other discretisation. They are the
unconverged tail of a dense real matrix. Taking the windows together,

    σR = 2:  no real eigenvalue with |μ|v_c ≥ 0.0830, so no real zero mode with |v| ≤ 12.05 v_c = 45.1
    σR = 4:  no real eigenvalue with |μ|v_c ≥ 0.0969, so no real zero mode with |v| ≤ 10.32 v_c = 25.5

provided the leading eigenvalue keeps falling with m beyond m = 21/2 as it does from m = ½ to 21/2 (0.083 and 0.094 at
the last sector). That is *numerical evidence*, a statement about the Galerkin eigenvalues of the discretisation that
two discretisations agree on, and not a bound on the spectrum.

### The identity ∫τ₁|ψ|² = ∫τ₁|χ|² = 0

If Sψ = vMψ with ψ ∈ Ḣ¹ ∩ L² (the eigenfunctions decay like r⁻², and the surface term r²|ψ|² then vanishes), then
⟨ψ, Sψ⟩ is imaginary, because S is anti-Hermitian, and ⟨ψ, Mψ⟩ is real, so Re(v)⟨ψ, Mψ⟩ = 0. For a real v this is
∫τ₁|ψ|² = 0, and the same holds for χ, whose equation is that of ψ at −v. This is #68's remark 2. The ψ and χ
eigenfunctions of a real λ are the eigenfunctions of S⁻¹M at −λ and at λ, both of which are real eigenvalues of Σ, so
evaluating the integral on every real eigenvalue of Σ covers both.

The script evaluates it on 80 to 112 eigenfunctions per sector, the leading 24 of each parity block and every real
eigenvalue above 10⁻⁸, together with the Gram matrix of M on the span of every eigenvector whose eigenvalue agrees to
10⁻¹⁰ (the check `real_identity`, the quantity |Re μ|/|μ| times the larger of the normalised integral and that Gram
matrix). The worst value is 0 to the digits printed, in all 22 sectors, and there is no real candidate in a window to
evaluate it on. The zero is structural. τ₁ is odd under x → −x, the inversion ψ(x) → ψ(−x) maps solutions at v to
solutions at v, and the channels split by the parity of l into two blocks, so an eigenfunction in either block has |ψ|²
even and ∫τ₁|ψ|² = 0 identically. For this profile the identity cannot exclude a real zero mode, as #63 says of any
τ that changes sign. It would constrain only an eigenvalue shared by both blocks with Re v ≠ 0.

The check can fail. Under `--mutate=rot` the planted eigenvalues are real and the identity is violated, since there
τ₁ > 0 and ∫τ₁|ψ|² > 0, and `real_identity` fires (§3). On the tanh weight it has not been tested by a case where it
would fail, since no real candidate exists to give it one.

## 6. The box error of #63's σR = 4 value

#63's σR = 4 value, 0.36732 to 0.36736, is the box L = 5. At σR = 2 its own runs show the box error, ρ(T₁)v_c rising
with L from 0.40938 at L = 8 to 0.40949 at L = 12 and 0.40952 at L = 16 at one spacing, dx = 0.25. For σR = 4 the
unmodified `frank_loss.py` was run at L = 7.5, N = 96 and L = 10, N = 128, both at the spacing dx = 0.15625 of L = 5,
N = 64 (`--part=krylov --sigmas=4 --Ns=<N> --L=<L> --vs=0.5`, logs in the raw folder).

| σR | L | N | dx | ρ(T₁)v_c, box | here | here − box |
|---|---|---|---|---|---|---|
| 4 | 5 | 64 | 0.156 | 0.36732 | 0.367657 | +3.4×10⁻⁴ |
| 4 | 5 | 128 | 0.078 | 0.36736 | 0.367657 | +3.0×10⁻⁴ |
| 4 | 7.5 | 96 | 0.156 | 0.36753 | 0.367657 | +1.3×10⁻⁴ |
| 4 | 10 | 128 | 0.156 | 0.36758 | 0.367657 | +7.7×10⁻⁵ |
| 2 | 8 | 128 | 0.125 | 0.40940 | 0.409556 | +1.6×10⁻⁴ |
| 2 | 16 | 128 | 0.25 | 0.40952 | 0.409556 | +3.6×10⁻⁵ |

The box value rises toward this one as L grows, by 2.1×10⁻⁴ from L = 5 to 7.5 and 5×10⁻⁵ from 7.5 to 10, a ratio of
0.24. If the increments keep that ratio the limit at dx = 0.156 is 0.36760, and #63's refinement at L = 5, N = 64 to 128
adds 4×10⁻⁵, which gives 0.36764 against 0.367657 here. At σR = 2 the increments +1.1×10⁻⁴ and +3×10⁻⁵ (ratio 0.27)
give 0.40953, and the +2×10⁻⁵ of N = 64 to 128 at L = 8 gives 0.40955 against 0.409556. The remaining differences are
1.7×10⁻⁵ at σR = 4 and 5×10⁻⁶ at σR = 2. These extrapolations assume a geometric ratio and are guesses, *numerical
evidence* at most. What the runs do show is that the disagreement of 7×10⁻⁴ at σR = 4 is the L = 5 box, and that it
shrinks with L at fixed spacing.

The frank_loss.py check `identity` fails in both runs and `#58 Picard value of G` in the L = 10 run. They compare with
#58's values at L = 5 and with an identity accurate only on the finer grid, and the frank_loss output says so of its
L = 12 and 16 runs. They are not the quantity used here.

## What this does and does not show

- ρ(T₁)v_c = 0.409556 at σR = 2 and 0.367657 at σR = 4, from a method with no periodic box, no FFT and no outer
  truncation, *numerical evidence*. Two discretisations differ by 5×10⁻¹² and 5×10⁻¹⁰, and eight parameter changes in
  the two leading sectors by at most 7×10⁻⁹. These are estimator ranges, not enclosures.
- They agree with #63's 0.4094 and 0.3674 to a relative 3.8×10⁻⁴ and 7.0×10⁻⁴ (8.1×10⁻⁴ against its five-digit
  0.36736). That is three digits in the sense of a relative difference below 10⁻³, and not in the sense of rounding at
  σR = 4 (0.368 against 0.367). The box runs of §6 show #63's value moving toward this one as L grows.
- The leading eigenvalues of #63 and the next quadruples are reproduced: ±0.0273 ± 0.4086i and ±0.0048 ± 0.2990i at
  σR = 2, ±0.0159 ± 0.3673i and ±0.0028 ± 0.2885i at σR = 4, each fourfold in T₁. FL's bound ρ√3 v_c ≤ 1 holds at 0.709
  and 0.637.
- No real eigenvalue of S⁻¹M with |μ|v_c ≥ 0.083 (σR = 2) or ≥ 0.097 (σR = 4), in sectors m = ½ to 21/2, so no real
  zero-mode speed |v| up to 12.05 v_c and 10.32 v_c, given that the leading eigenvalue keeps falling with m. That is
  two and a half to three times #63's 4 v_c. It is *numerical evidence*. That E is empty for the tanh profile is still
  a *conjecture*, and the window does not reach large v.
- The identity ∫τ₁|ψ|² = ∫τ₁|χ|² = 0 was evaluated on 80 to 112 eigenfunctions per sector and is zero identically, by
  the inversion symmetry of the odd weight. It has no real candidate to test and cannot exclude one for this profile.
- Nothing here changes #63's theorem or the existence statements. They depend on §5 to §7 of frank_loss.md, which this
  note does not examine. What it examines is §9's numerics, by a method that shares none of their machinery.
- No novelty is claimed. The reduction, the sector identity and the planted spectrum are #63's, #68's and Frank–Loss's.
  The Galerkin approximation of compact operators and the Ω_κ spin-angular functions are standard and were recalled, not
  re-read *(unverified)*, for example J. E. Osborn, "Spectral approximation for compact operators", Math. Comp. 29
  (1975), for the first.

## Discrepancies and reading notes

1. The three-digit agreement the issue asks for holds at σR = 2 (0.4096 and 0.4094, relative 3.8×10⁻⁴) and fails by
   rounding at σR = 4 (0.3677 and 0.3674, 0.368 against 0.367, relative 7.0×10⁻⁴). #63's σR = 4 value is the L = 5 box,
   and §6 shows it rising with L at fixed spacing. The disagreement is therefore #63's box error at σR = 4, which its
   own table measured at σR = 2 only (its closing list says "the box error estimated from L = 12 and 16 at σR = 2
   only"). frank_loss.md is not edited here.
2. frank_loss.md §9 gives the leading eigenvalues at σR = 2 as ±0.0270 ± 0.4085i (N = 64, L = 8) and ±0.0271 ± 0.4085i
   (N = 96, 128). Here they are ±0.02730 ± 0.40865i. The real part differs by 1.1%, the imaginary part by 4×10⁻⁴. At
   σR = 4 the real part is 0.0159 here and 0.0154 there, a 3% difference. The real part is one to two orders of
   magnitude smaller than the imaginary part, so the same absolute box error is a larger relative error. It does not
   change any conclusion.
3. #63 reads "no real zero mode of the grid operator up to 4 v_c" from s_min, and none of the 24 largest eigenvalues real,
   which excludes real speeds below 3.35 v_c and 3.47 v_c. The present windows contain the whole spectrum above
   |μ|v_c = 0.083 and 0.097, so they contain those 24 and no real eigenvalue. The two agree, and the windows reach
   12.05 v_c and 10.32 v_c.
4. The residual of the leading eigenfunctions at σR = 4 is 2×10⁻⁴ to 1.4×10⁻³ with p = 100, while the eigenvalues have
   converged to 5×10⁻¹⁰. A diagnostic that uses the eigenfunction, such as the identity, is not converged to more than
   that residual at σR = 4. The structural zero of §5 does not depend on it.

### Corrections made to the harness while building it

- The first full run failed `real_identity` in six sectors (σR = 2, m = 17/2, 19/2, 21/2, and σR = 4, m = 15/2, 19/2,
  21/2), with a worst value between 2×10⁻⁷ and 2×10⁻⁶. The cluster of "equal" eigenvalues over which the Gram matrix of
  M is formed had a tolerance of 10⁻⁶|μ|, which paired two distinct eigenvalues of the two parity blocks (examined at
  σR = 2, m = 17/2, rank 20 and 21, |μ|v_c = 0.0606) that differ by 3×10⁻⁸ and sit within 3×10⁻⁷ of the imaginary
  axis, so a nonzero cross term met a nonzero Re μ. For exact
  eigenvalues the identity holds on each eigenvector separately and on a degenerate eigenspace, not across two
  distinct eigenvalues. The tolerance is now 10⁻¹⁰|μ|. The failing log is kept in the raw folder as
  `first-run/run-before-cluster-fix.log`, its leading eigenvalues are those of the final run to 10⁻⁹, and the final run
  passes.
- The first version of `radial_int` divided each error by the integral plus 10⁻¹⁴, so integrals that vanish
  exactly gave a relative error of 1. It now divides by the largest entry of the same integral.
- At σR = 4 the first full run used p = 80 and its eigenfunction residual was 2×10⁻³ to 8×10⁻³. The run quoted uses
  p = 100, with residuals of 2×10⁻⁴ to 1.4×10⁻³, and the leading eigenvalues did not change beyond 10⁻⁹.
- A scratch convergence experiment, not part of the script, ran for about two minutes without `OMP_NUM_THREADS`
  and used the default BLAS threads, above the four-thread budget. It was stopped and repeated single-threaded. All
  runs of the script set one thread per process.

## Recommended tasks

### Task 1, an independent reimplementation of the eigenvalues

Compute ρ(T₁)v_c and the next two quadruples at σR = 2 and 4 by a method that shares neither this note's basis nor
#63's box, for example finite elements in (r, θ) per J_z sector. Why: the Galerkin approximation and the Ω_κ basis are
this note's only route, and the exhaustion of Ḣ¹ by V_h is not proved. Acceptance, each can fail: the planted
6/(1 + r²) returns ρ = 1 to 10⁻⁶ first, and ρ(T₁)v_c agrees with 0.409556 and 0.367657 to 10⁻⁶, or the difference is
reported with its cause. Depends on #75.

### Task 2, a hand bound on the leading eigenvalue of each J_z sector

The claim of §5 that no real eigenvalue exists up to |v| = 12.05 v_c needs the leading eigenvalue to keep falling with
m beyond 21/2, which is observed to m = 21/2 and not proved. Bound ‖T₁‖ on the J_z = m sector, for instance through the
centrifugal term in ‖∇ψ‖², and show it falls with m. Acceptance, each can fail: the bound exceeds every computed
leading |μ| in the output's tables (it fails if it undercuts one), and it gives a value of m beyond which no sector
reaches |μ|v_c = 0.083 at σR = 2 and 0.097 at σR = 4. Depends on #75.

### Task 3, the thin walls σR = 8 and 16

The same code with the radial basis enlarged until `conv` shows 10⁻⁸ (a shell of width 1/σ = 0.06 at σR = 16). Why: the
windows and the ρ(T₁)v_c column of #63 stop at σR = 4 and the walls of interest are thinner, and a real eigenvalue
could appear when the leading eigenvalues approach the imaginary axis more slowly. Acceptance, each can fail: the
planted checks pass at those sizes, A and B agree to 10⁻⁶, ρ(T₁)v_c agrees with the box value of frank_loss.py at its
largest affordable L within the box error, and the windows and any real candidates are reported with the identity
evaluated on each. Depends on #75.

### Task 4, set #69's real zero-mode search against these windows

When #69 reports, compare its largest v with no zero mode, per σR, with 12.05 v_c and 10.32 v_c here. Acceptance:
the two statements are consistent, or the discrepancy is written up with the sector, the eigenvalue and the
discretisation that produced it, and neither is dismissed. Depends on #69 and #75.

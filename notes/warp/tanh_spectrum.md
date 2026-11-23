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
#94 recomputed every sector's leading eigenvalue with B-splines in (r, θ), a basis that shares nothing with this one,
and matched this note's tables to 6×10⁻⁶ ([tanh_spectrum_check.md](tanh_spectrum_check.md)). #111 and #118 refined
that basis in all 22 sectors until its windows start below this note's 0.0830 and 0.0969, and matched the second
modulus to 5.0×10⁻⁶. It finds no real eigenvalue in any window, so it confirms the absence of real zero modes to
12.05 v_c at σR = 2 and to 10.60 v_c at σR = 4. The second line is set by the leading eigenvalue of m = 21/2 rather
than by a window, so it is above this note's 10.32 v_c. Both rest on the same condition as this note's, that the
leading eigenvalue keeps falling beyond m = 21/2. #120 computed it for m = 23/2 to 41/2, where it falls strictly and stays
below both lines, so the condition now concerns m > 41/2 only. #125 recomputed m = 23/2, 29/2, 35/2 and 41/2 with this
note's code, within 3.0×10⁻⁶ of #120's ([tanh_spectrum_check.md](tanh_spectrum_check.md)).

The identity ∫τ₁|ψ|² = ∫τ₁|χ|² = 0 has nothing to evaluate on, since there is no real candidate, and for this profile it
would hold without one (§5).

Thin walls (#96, §7). ρ(T₁)v_c = 0.29951080 at σR = 8 and 0.22508330 at σR = 16. A and B agree to 1.4×10⁻¹¹ and
1.5×10⁻¹⁰, and #94's B-spline basis gives the same leading moduli to 4×10⁻⁷. No real eigenvalue lies in any window, so
a real zero mode needs |v| > 9.57 v_c at σR = 8 and > 8.87 v_c at σR = 16, given that the leading eigenvalue keeps
falling beyond m = 21/2, which is assumed there. At σR = 8 the box value at L = 3.5 is 3.9×10⁻⁴ below this one, and the
box error does not clearly account for that difference. At σR = 16 there is no box value. Status: *numerical evidence*.

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
the last sector). #120 computed that fall for m = 23/2 to 41/2 in the B-spline basis
([tanh_spectrum_check.md](tanh_spectrum_check.md)), so it is assumed for m > 41/2 only. That is *numerical evidence*, a statement about the Galerkin eigenvalues of the discretisation that
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
N = 64, as `python notes/warp/frank_loss.py --part=krylov --sigmas=4 --Ns=96 --L=7.5 --vs=0.5` and the same with
`--Ns=128 --L=10`, at `nice -n 10` and one thread. The logs are in the raw folder.

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

## 7. Thin walls, σR = 8 and 16 (#96)

Issue #96. Output [tanh_spectrum_thin_output.txt](tanh_spectrum_thin_output.txt): the two sector logs, σR = 8 and 16, each
under a header with its machine, code commit and command. Raw output in
[negative-mass-raw](https://github.com/nmoinvaz/negative-mass-raw): `I96-win` (the laptop: sizes, pilot, σR = 8 sectors,
box runs; commit 7940b1e), `I96-mac` (box probes, a477e61), `I96-dh610` (the σR = 16 sectors, 68abc50) and `I96-rpi5`
(numbers that no log prints, 921d262). Code: branch `wip/i96-thin` at c5fd40d, which adds `--p`, `--npan` and `--convtol`
to tanh_spectrum.py and guards frank_loss.py's `import resource` so it starts on Windows (peak memory prints as nan there,
so the peaks below are the memory sampler's). The scripts are otherwise those of #75 and #63. Status: *numerical
evidence*.

Every number below is taken from a log or file of the raw folders, named beside it, or from
`I96-rpi5/sector_numbers_output.txt`. That script reads the sector logs and the `eig/*.npz` and prints the sector and
window tables of the two σR (below). The σR = 8 tables were drafted from `I96-win/table.py`'s output, and they equal the
rows of `sector_numbers.py`'s.

### What ran

Commands run from the repository root with `python` the project's venv, one thread per process (`OMP_NUM_THREADS = 1`,
`--procs=1`), `DIR` the raw folder's `eig`.

| run | where | command | exit | wall |
|---|---|---|---|---|
| sizes, σR = 8 | win | `tanh_spectrum.py --part=conv --sigmas=8 --p=8:P --npan=N --convtol=1e-8 --procs=1` | 1, 0 | 2073 s, 4160 s |
| pilot, σR = 16, m = ½ | win | `tanh_spectrum.py --part=vc,planted,tanh --sigmas=16 --ms=0.5 --p=16:320 --npan=320 --procs=1 --raw=DIR` | 0 | 8615 s |
| sizes, σR = 16 | win | `tanh_spectrum.py --part=conv --sigmas=16 --p=16:P --npan=N --convtol=1e-8 --procs=1` | 1, 0 | 12,635 s, 14,259 s |
| sectors, σR = 8 | win | `tanh_spectrum.py --part=vc,planted,shift,tanh --sigmas=8 --p=8:200 --npan=160 --procs=1 --raw=DIR` | 0 | 24,270 s |
| sectors, σR = 16 | dh610 | the same with `--sigmas=16 --p=16:340 --npan=320` | 0 | 97,906 s (27.2 h) |
| box, σR = 8 | win, Mac | `frank_loss.py --part=krylov --sigmas=8 --Ns=N --L=L --vs=0.5` | see the box table | |
| recount of the sector logs | rpi5 | `I96-rpi5/sector_numbers.py` | 1 (R1 fails at σR = 8, Reading notes 6) | 1 s |

Logs: `I96-win/conv-s8-p160-n160.log`, `conv-s8-p200-n160.log`, `pilot-s16.log`, `conv-s16-p320-n320.log`,
`conv-s16-p340-n320.log`, `sectors-s8.log` (each with its `.exit`), `I96-dh610/sectors-s16-p340.log` (with its `.exit`),
and the box logs of the box table. The laptop (4 physical cores, 32 GB, shared with another project's solver) ran at most
two of these at once, at below-normal priority, detached. Its peak working sets are the largest `peak=` of
`memory-sampler.log` for the run's processes, sampled once a minute (every 10 s for the box run of L = 3.5). Nothing was
stopped over memory. dh610 (20 cores, 62 GB, shared with that project's solver) ran the sectors at nice 10 as a systemd
user unit, with a memory watcher (`memwatch.log`, once a minute) that only reports.

### Sizes

The `conv` part changes one parameter at a time from the base (dK = 20, p = P, Lmap = 2, N panels), in m = ½ and 3/2,
and the check takes the largest relative move of the leading |μ|v_c, not counting the deliberately coarse dK = 12. The
tolerance is 10⁻⁸ here, against #75's 10⁻⁷ (`--convtol`).

| σR | attempt | P | N | largest move | variant | wall | peak working set | exit |
|---|---|---|---|---|---|---|---|---|
| 8 | 1 | 160 | 160 | 5.2×10⁻⁸ | Lmap = 4 (m = 3/2); next p = 140 at 2.5×10⁻⁸ | 2073 s | 2.54 GiB | 1, `conv` fails |
| 8 | 2 | 200 | 160 | 1.1×10⁻⁹ | Lmap = 4 (m = ½) | 4160 s | 3.90 GiB | 0 |
| 16 | 1 | 320 | 320 | 2.4×10⁻⁸ | Lmap = 4 (m = 3/2); next Lmap = 4 (m = ½) at 1.9×10⁻⁸ | 12,635 s | 9.78 GiB | 1, `conv` fails |
| 16 | 2 | 340 | 320 | 7.2×10⁻⁹ | Lmap = 4 (m = 3/2) | 14,259 s | 11.03 GiB | 0 |

The largest moves are the `conv` check's own line in each log. The variant and the runner-up are read from the logs'
printed changes, which have two digits, so a printed ratio can differ from the check's value in the second digit.

At σR = 8 the failing variants (Lmap = 4, p = 140) point to the radial size, not to the panels: 320 panels change the
leading eigenvalue by 6.1×10⁻¹⁴ (absolute, m = ½) at P = 160. So P rose by 40 (the task's rule) with N unchanged. At
σR = 16 the same single variant limits: at P = 320 every other variant moves the eigenvalue by at most 2.4×10⁻⁹, and at
P = 340 the variants p = 320 and p = 360 move it by at most 5.5×10⁻¹⁰, while Lmap = 4 moves it by 7.2×10⁻⁹. The step
from P = 320 was +20, not the task's +40: the coordinator changed it, since sectors at P = 360 would peak near 31 GiB
(24.43 GiB of the pilot scaled by (360/320)², an estimate) and do not fit the laptop. The σR = 16 sectors went to dh610
at P = 340. The pass has 28% of the tolerance to spare and rests on one variant.

Pilot (step 3, `pilot-s16.log`): σR = 16, m = ½, P = 320, N = 320, alone, 8605 s for the sector, peak working set
24.43 GiB against 20 GiB estimated. ρ(T₁)v_c = 0.22508330 in A and in B, differing by 4.4×10⁻¹⁰, window 0.1128 with no
real eigenvalue in it. The sector at P = 340 gives 0.22508330 again: its A differs from the pilot's A by 2.1×10⁻¹¹ and its
B (P = 350) from the pilot's B (P = 330) by 3.2×10⁻¹⁰ (`I96-rpi5/sector_numbers_output.txt`, R3).

### The σR = 8 sectors

`sectors-s8.log`: discretisation A has dK = 28, p = 200, Lmap = 2 and 160 panels of 16 nodes, B has dK = 38, p = 210,
Lmap = 3 and 240 panels. All 11 sectors, m = ½ to 21/2, exit 0, "all checks passed", wall 24,270 s (6.74 h). A sector took
2001 to 2521 s, beside a σR = 16 `conv` run from 14:53 on (`queue.log`). The peak working set was 9.93 GiB.

The checks that run first passed. `vc` matches mass_bounds.md at all four σR to 3.7×10⁻⁶. `planted` returns ρ = 1 to
7.3×10⁻¹⁵ (`real_identity` 2.3×10⁻¹⁶) and `shift` to 1.6×10⁻¹¹ (max |Re μ|/|μ| 7.7×10⁻¹²). These two run at sizes of
their own that the script fixes (planted dK = 12, shift dK = 22, p = 60, 120 panels, Lmap = 2, m = ½, 3/2, 5/2), not at
`--p` and `--npan`, and the planted weight does not depend on σR. They show that the code and the Ω basis are right at
their sizes, not that P = 200 resolves the wall.

σR = 8, v_c = 1.562434, ρ(T₁)v_c = 0.29951080 (A and B, differing by 1.4×10⁻¹¹), √3 ρ v_c = 0.5188 (FL: ≤ 1).

| m | leading μv_c (a quadruple ±a ± ib) | \|μ\|v_c, A | \|μ\|v_c, B | next \|μ\|v_c | stable to 10⁻⁵ above | symmetric to |
|---|---|---|---|---|---|---|
| 0.5 | 0.00723 + 0.29942i | 0.299511 | 0.299511 | 0.224930 | 0.1045 | 2.7×10⁻⁹ |
| 1.5 | 0.00116 + 0.25251i | 0.252510 | 0.252510 | 0.199123 | 0.0988 | 4.4×10⁻¹⁰ |
| 2.5 | 0.00023 + 0.21708i | 0.217079 | 0.217079 | 0.177495 | 0.1034 | 5.6×10⁻¹⁰ |
| 3.5 | 0.00005 + 0.19041i | 0.190406 | 0.190406 | 0.159828 | 0.0975 | 7.6×10⁻¹⁰ |
| 4.5 | 0.00001 + 0.16964i | 0.169640 | 0.169640 | 0.145264 | 0.0922 | 9.9×10⁻¹⁰ |
| 5.5 | 0.00000 + 0.15300i | 0.153001 | 0.153001 | 0.133092 | 0.0873 | 1.2×10⁻⁹ |
| 6.5 | 0.00000 + 0.13936i | 0.139358 | 0.139358 | 0.122779 | 0.0903 | 2.5×10⁻¹¹ |
| 7.5 | 0.00000 + 0.12796i | 0.127962 | 0.127962 | 0.113936 | 0.0856 | 3.4×10⁻¹¹ |
| 8.5 | 0.00000 + 0.11830i | 0.118296 | 0.118296 | 0.106273 | 0.0813 | 4.6×10⁻¹¹ |
| 9.5 | 0.00000 + 0.10999i | 0.109993 | 0.109993 | 0.099569 | 0.0773 | 9.3×10⁻¹¹ |
| 10.5 | 0.00000 + 0.10278i | 0.102782 | 0.102782 | 0.093657 | 0.0738 | 1.2×10⁻¹⁰ |

The leading eigenvalue falls with m, from 0.299511 to 0.102782, with m|μ|v_c = 1.08 at m = 21/2 (0.87 at σR = 2, §4).
Its real part falls from 7.2×10⁻³ at m = ½ to 5.1×10⁻⁵ at m = 7/2 and 3×10⁻⁶ at m = 11/2 (the log's first four
eigenvalues of A). The spectrum is symmetric under μ → −μ to at most 2.7×10⁻⁹ of |μ| over the converged eigenvalues.
The L² residual ‖Sψ − vMψ‖/‖Sψ‖ of the leading eigenfunctions rises with m from 2.3×10⁻⁵ at m = ½ to 2.8×10⁻³ at
m = 21/2, while A and B agree on every leading eigenvalue to the digits printed. As in §4 the eigenvalue converges far
faster than the eigenfunction, and a diagnostic made from an eigenfunction is no better than that residual.

### Real eigenvalues and the window, σR = 8

| m | window \|μ\|v_c ≥ | eigenvalues in it (A / B) | real in it (A / B) | real with \|μ\|v_c > 10⁻⁸ (A / B) | largest | reproduced by B | identity, worst | eigenfunctions (real) |
|---|---|---|---|---|---|---|---|---|
| 0.5 | 0.1045 | 34 / 36 | 0 / 0 | 88 / 92 | 7.8×10⁻⁴ | 0 | 0 | 135 (87) |
| 1.5 | 0.0988 | 34 / 36 | 0 / 0 | 103 / 94 | 7.4×10⁻⁴ | 0 | 0 | 146 (98) |
| 2.5 | 0.1034 | 26 / 28 | 0 / 0 | 102 / 94 | 6.9×10⁻⁴ | 0 | 0 | 146 (98) |
| 3.5 | 0.0975 | 24 / 28 | 0 / 0 | 102 / 94 | 6.5×10⁻⁴ | 0 | 0 | 145 (97) |
| 4.5 | 0.0922 | 26 / 28 | 0 / 0 | 102 / 94 | 6.2×10⁻⁴ | 0 | 0 | 146 (98) |
| 5.5 | 0.0873 | 26 / 28 | 0 / 0 | 102 / 94 | 5.9×10⁻⁴ | 0 | 0 | 144 (96) |
| 6.5 | 0.0903 | 16 / 20 | 0 / 0 | 102 / 94 | 5.6×10⁻⁴ | 0 | 0 | 145 (97) |
| 7.5 | 0.0856 | 16 / 20 | 0 / 0 | 102 / 93 | 5.3×10⁻⁴ | 1 | 0 | 144 (96) |
| 8.5 | 0.0813 | 16 / 20 | 0 / 0 | 102 / 93 | 5.1×10⁻⁴ | 0 | 0 | 144 (96) |
| 9.5 | 0.0773 | 18 / 20 | 0 / 0 | 101 / 92 | 4.8×10⁻⁴ | 0 | 0 | 144 (96) |
| 10.5 | 0.0738 | 18 / 20 | 0 / 0 | 101 / 92 | 4.6×10⁻⁴ | 0 | 0 | 144 (96) |

The window of a sector is |μ|v_c above its "stable to 10⁻⁵" value, 0.0738 to 0.1045. There is no real eigenvalue in
any window, in A or in B, and the windows hold 254 (A) and 284 (B) eigenvalues together. Below them each sector has 88
to 103 real eigenvalues in A and 92 to 94 in B, all with |μ|v_c ≤ 7.8×10⁻⁴. B reproduces one of A's to 10⁻⁴: at
m = 15/2, A's −2.16459×10⁻⁸ against B's −2.16478×10⁻⁸, a relative difference of 8.6×10⁻⁵ and an absolute one of
1.9×10⁻¹², just above the 10⁻⁸ floor of the search (`reals_check.py 8`). None of the other 1106 of A's 1107 real
eigenvalues is reproduced (Reading notes, 1). The identity `real_identity` is zero to the digits printed on the 135 to
146 eigenfunctions evaluated per sector, the leading 24 of each parity block and the real ones. For this profile the zero
is structural, as in §5, and no real candidate stands in a window to evaluate it on.

The window line of the log:

    window: every eigenvalue with |mu| v_c >= 0.1045 is converged in every m listed, and the leading one of m = 10.5 is
    0.1028; so a real zero mode needs |v| > 9.57 v_c = 14.95

That is |v| > 9.57 v_c at σR = 8, against 12.05 v_c at σR = 2 and 10.32 v_c at σR = 4 (§5), and 14.95 in v against 45.1
and 25.5. The window of σR = 8 starts above the leading eigenvalue of m = 21/2, so, unlike the two narrower walls, the
line is set by the converged window of m = ½ and not by the last sector. It carries §5's condition, that the leading
eigenvalue keeps falling with m beyond 21/2, and #120's computation of that fall (m = 23/2 to 41/2) covers σR = 2 and
4 only. At σR = 8 the condition is assumed for every m > 21/2. This is *numerical evidence* about the Galerkin
eigenvalues, not a bound.

### The σR = 16 sectors

`I96-dh610/sectors-s16-p340.log`: discretisation A has dK = 28, p = 340, Lmap = 2 and 320 panels of 16 nodes, B has
dK = 38, p = 350, Lmap = 3 and 480 panels. All 11 sectors, m = ½ to 21/2, exit 0, "all checks passed", wall 97,906 s
(27.2 h). A sector took 5635 to 20,362 s, m = ½ the longest. The worker's peak resident memory printed with each sector is
28.39 to 28.40 GiB. The memory watcher's largest sum of resident sets was 28.67 GiB and its lowest MemAvailable 19.34 GiB
(`I96-dh610/memwatch.log`, `I96-rpi5/sector_numbers_output.txt`, R4).

The checks that run first passed. `vc` matches mass_bounds.md to 3.7×10⁻⁶, `planted` returns ρ = 1 to 7.3×10⁻¹⁵
(`real_identity` 1.0×10⁻¹⁶) and `shift` to 1.6×10⁻¹¹ (max |Re μ|/|μ| 7.7×10⁻¹²). As at σR = 8 they run at the sizes of
their own that the script fixes, not at `--p` and `--npan`.

σR = 16, v_c = 0.985030, ρ(T₁)v_c = 0.22508330 (A and B, differing by 1.5×10⁻¹⁰), √3 ρ v_c = 0.3899 (FL: ≤ 1).

| m | leading μv_c (a quadruple ±a ± ib) | \|μ\|v_c, A | \|μ\|v_c, B | next \|μ\|v_c | stable to 10⁻⁵ above | symmetric to |
|---|---|---|---|---|---|---|
| 0.5 | 0.00278 + 0.22507i | 0.225083 | 0.225083 | 0.184364 | 0.1128 | 5.3×10⁻¹⁰ |
| 1.5 | 0.00037 + 0.20079i | 0.200790 | 0.200790 | 0.169416 | 0.1075 | 2.0×10⁻¹¹ |
| 2.5 | 0.00006 + 0.18081i | 0.180812 | 0.180812 | 0.155971 | 0.1025 | 2.2×10⁻¹¹ |
| 3.5 | 0.00001 + 0.16471i | 0.164708 | 0.164708 | 0.144418 | 0.1063 | 2.2×10⁻¹¹ |
| 4.5 | 0.00000 + 0.15143i | 0.151426 | 0.151426 | 0.134462 | 0.1011 | 2.5×10⁻¹¹ |
| 5.5 | 0.00000 + 0.14024i | 0.140244 | 0.140244 | 0.125805 | 0.0964 | 2.6×10⁻¹¹ |
| 6.5 | 0.00000 + 0.13067i | 0.130674 | 0.130674 | 0.118209 | 0.0920 | 5.5×10⁻¹¹ |
| 7.5 | 0.00000 + 0.12237i | 0.122374 | 0.122374 | 0.111489 | 0.0947 | 2.9×10⁻¹¹ |
| 8.5 | 0.00000 + 0.11510i | 0.115098 | 0.115098 | 0.105500 | 0.0904 | 6.1×10⁻¹¹ |
| 9.5 | 0.00000 + 0.10866i | 0.108659 | 0.108659 | 0.100127 | 0.0865 | 8.2×10⁻¹¹ |
| 10.5 | 0.00000 + 0.10292i | 0.102918 | 0.102918 | 0.095278 | 0.0829 | 9.7×10⁻¹¹ |

The leading eigenvalue falls with m, from 0.225083 to 0.102918, with m|μ|v_c = 1.08 at m = 21/2 (1.08 at σR = 8, 0.87 at
σR = 2). Its real part falls from 2.8×10⁻³ at m = ½ to 3.7×10⁻⁴ at m = 3/2 and 1.1×10⁻⁵ at m = 7/2 (the log's first four
eigenvalues of A). The spectrum is symmetric under μ → −μ to at most 5.3×10⁻¹⁰ of |μ| over the converged eigenvalues.
A and B differ by 9.6×10⁻¹¹ to 5.5×10⁻¹⁰ on the leading |μ|v_c of each sector (`sector_numbers.py`, R2).

The L² residual ‖Sψ − vMψ‖/‖Sψ‖ of the leading eigenfunctions rises with m at both σR, and in every sector it is larger
at σR = 16 (P = 340) than at σR = 8 (P = 200), from the logs:

| m | σR = 8 | σR = 16 |
|---|---|---|
| 0.5 | 2.3×10⁻⁵ | 8.3×10⁻⁵ |
| 1.5 | 2.6×10⁻⁵ | 1.1×10⁻⁴ |
| 2.5 | 3.2×10⁻⁵ | 1.7×10⁻⁴ |
| 3.5 | 4.4×10⁻⁵ | 3.2×10⁻⁴ |
| 4.5 | 7.6×10⁻⁵ | 6.0×10⁻⁴ |
| 5.5 | 1.5×10⁻⁴ | 1.1×10⁻³ |
| 6.5 | 3.0×10⁻⁴ | 1.8×10⁻³ |
| 7.5 | 5.7×10⁻⁴ | 2.9×10⁻³ |
| 8.5 | 1.0×10⁻³ | 4.4×10⁻³ |
| 9.5 | 1.7×10⁻³ | 6.5×10⁻³ |
| 10.5 | 2.8×10⁻³ | 9.0×10⁻³ |

This is data. No cause is proposed: no run here tests one. In the same sectors A and B agree on the leading eigenvalue to
5.5×10⁻¹⁰ at most, so, as in §4, the eigenvalue converges far faster than the eigenfunction, and a diagnostic made from an
eigenfunction is no better than that residual.

### Real eigenvalues and the window, σR = 16

| m | window \|μ\|v_c ≥ | eigenvalues in it (A / B) | real in it (A / B) | real with \|μ\|v_c > 10⁻⁸ (A / B) | largest | reproduced by B | identity, worst | eigenfunctions (real) |
|---|---|---|---|---|---|---|---|---|
| 0.5 | 0.1128 | 20 / 24 | 0 / 0 | 76 / 74 | 5.1×10⁻⁴ | 0 | 0 | 124 (76) |
| 1.5 | 0.1075 | 22 / 24 | 0 / 0 | 76 / 74 | 4.5×10⁻⁴ | 0 | 0 | 124 (76) |
| 2.5 | 0.1025 | 22 / 24 | 0 / 0 | 76 / 74 | 4.3×10⁻⁴ | 0 | 0 | 124 (76) |
| 3.5 | 0.1063 | 16 / 20 | 0 / 0 | 76 / 73 | 4.0×10⁻⁴ | 0 | 0 | 124 (76) |
| 4.5 | 0.1011 | 18 / 20 | 0 / 0 | 76 / 73 | 3.8×10⁻⁴ | 0 | 0 | 124 (76) |
| 5.5 | 0.0964 | 16 / 20 | 0 / 0 | 76 / 72 | 3.6×10⁻⁴ | 0 | 0 | 124 (76) |
| 6.5 | 0.0920 | 18 / 20 | 0 / 0 | 75 / 72 | 3.4×10⁻⁴ | 0 | 0 | 123 (75) |
| 7.5 | 0.0947 | 12 / 16 | 0 / 0 | 75 / 72 | 3.2×10⁻⁴ | 0 | 0 | 123 (75) |
| 8.5 | 0.0904 | 14 / 16 | 0 / 0 | 74 / 72 | 3.1×10⁻⁴ | 0 | 0 | 122 (74) |
| 9.5 | 0.0865 | 14 / 16 | 0 / 0 | 74 / 72 | 3.0×10⁻⁴ | 0 | 0 | 122 (74) |
| 10.5 | 0.0829 | 14 / 16 | 0 / 0 | 74 / 72 | 2.8×10⁻⁴ | 0 | 0 | 122 (74) |

The window of a sector is |μ|v_c above its "stable to 10⁻⁵" value, 0.0829 to 0.1128. There is no real eigenvalue in any
window, in A or in B, and the windows hold 186 (A) and 216 (B) eigenvalues together. Below them each sector has 74 to 76
real eigenvalues in A and 72 to 74 in B, all with |μ|v_c ≤ 5.1×10⁻⁴, and B reproduces none of A's. The identity
`real_identity` is zero to the digits printed on the 122 to 124 eigenfunctions evaluated per sector, the leading 24 of each
parity block and the real ones. The real ones among them equal the table's real count in every sector: the script's
unscaled floor 10⁻⁸ lies below the table's, 10⁻⁸/v_c, since v_c < 1, and no real eigenvalue falls between
(`I96-rpi5/reals_check_s16_output.txt`; at σR = 8 some do, Reading notes 2). For this profile the zero is structural, as in
§5, and no real candidate stands in a window to evaluate it on.

The window line of the log:

    window: every eigenvalue with |mu| v_c >= 0.1128 is converged in every m listed, and the leading one of m = 10.5 is
    0.1029; so a real zero mode needs |v| > 8.87 v_c = 8.74

That is |v| > 8.87 v_c at σR = 16, against 9.57 v_c at σR = 8, 10.32 v_c at σR = 4 and 12.05 v_c at σR = 2, and 8.74 in v
against 14.95, 25.5 and 45.1. As at σR = 8 the window starts above the leading eigenvalue of m = 21/2, so the line is set
by the converged window of m = ½ and not by the last sector. It carries §5's condition, that the leading eigenvalue keeps
falling with m beyond 21/2, which #120 covers at σR = 2 and 4 only. At σR = 16 the condition is assumed for every m > 21/2.
This is *numerical evidence* about the Galerkin eigenvalues, not a bound.

### The thin walls against σR = 2 and 4

| σR | v_c | ρ(T₁)v_c | A and B differ by | √3 ρ v_c (FL: ≤ 1) | window, \|μ\|v_c ≥ | real zero mode needs \|v\| > |
|---|---|---|---|---|---|---|
| 2 | 3.744900 | 0.409556 | 5×10⁻¹² | 0.7094 | 0.0830 | 12.05 v_c = 45.1 |
| 4 | 2.470971 | 0.367657 | 5×10⁻¹⁰ | 0.6368 | 0.0969 | 10.32 v_c = 25.5 |
| 8 | 1.562434 | 0.29951080 | 1.4×10⁻¹¹ | 0.5188 | 0.1045 | 9.57 v_c = 14.95 |
| 16 | 0.985030 | 0.22508330 | 1.5×10⁻¹⁰ | 0.3899 | 0.1128 | 8.87 v_c = 8.74 |

The σR = 16 pilot (m = ½ only, P = 320) gives the same 0.22508330. The σR = 2 and 4 rows are those of the Answer's table
and the two lines of §5.

### Box comparison, σR = 8

The wall is 1/8 wide, and dx = 2L/N = 1/32 puts four grid points across it. §6's method: hold dx and raise L,
`python notes/warp/frank_loss.py --part=krylov --sigmas=8 --Ns=N --L=L --vs=0.5`, one thread, c5fd40d. Raw folders
`I96-mac` (L = 1, 2, 2.5) and `I96-win` (L = 1 and L = 3, 3.5).

| L | N | ρ(T₁)v_c, box | here − box | `gmres` residual (target 10⁻¹⁰) | of it the k = 0 part | `identity` | peak memory |
|---|---|---|---|---|---|---|---|
| 1 | 64 | 0.21589 | +8.4×10⁻² | 1.1×10⁻² fail | 1.1×10⁻² | 2.4×10⁻¹ fail | 0.75 GiB, Mac |
| 2 | 128 | 0.29732 | +2.19×10⁻³ | 1.6×10⁻⁹ fail | 1.6×10⁻⁹ | 2.3×10⁻⁸ | 4.75 GiB, Mac |
| 2.5 | 160 | 0.29843 | +1.08×10⁻³ | 5.1×10⁻¹³ | 4.9×10⁻¹³ | 5.2×10⁻¹² | 8.37 GiB, Mac |
| 3 | 192 | 0.29890 | +6.1×10⁻⁴ | 1.3×10⁻¹³ | 1.5×10⁻¹⁶ | 2.5×10⁻¹⁶ | 12.77 GiB, win |
| 3.5 | 224 | 0.29912 | +3.9×10⁻⁴ | 1.3×10⁻¹³ | 1.4×10⁻¹⁹ | 1.4×10⁻¹⁵ | 20.69 GiB, win |

"here − box" is 0.29951080 minus the box value, which the logs print to five digits, so each difference is good to about
5×10⁻⁶. Logs: `I96-win/box-s8-L1-N64.log` (the same numbers as `I96-mac/box-s8-L1-N64.log`),
`I96-mac/box-s8-L2-N128.log`, `box-s8-L2.5-N160.log`, `I96-win/box-s8-L3-N192.log` and `box-s8-L3.5-N224.log`. The Mac
peaks are the maximum resident set (the script's own figure, in GiB), the laptop's are the sampler's peak working set,
so they are different quantities. L = 3.5 is the largest L whose N fits the laptop's 24 GiB rule alone: it took 1361 s
and 20.69 GiB, 2% above the 20.3 GiB that N³ growth from L = 3 gives, and the machine's free memory fell to 3.7 GiB.
L = 4 (N = 256) would need about 31 GiB (20.69 GiB × (256/224)³, an estimate) and was not run. The L = 3 run's exit code
was not captured (its `.exit` file is empty); its log ends "all checks passed" and its `.err` is empty. L = 3.5 exited 0.

L = 1 and 2 fail `gmres`, and the failure is entirely in the k = 0 part, which falls from 1.1×10⁻² at L = 1 to 1.5×10⁻¹⁶
at L = 3 (the target is unchanged). So these two values belong to a box that is too small, and only L ≥ 2.5 is a valid
solve. The box value rises toward the sectors' as L grows at this spacing, by 1.11×10⁻³ from L = 2 to 2.5, 4.7×10⁻⁴
from 2.5 to 3 and 2.2×10⁻⁴ from 3 to 3.5, each good to about 10⁻⁵. The ratios of successive increments are 0.42 and 0.47
(the second ±0.03), so the ratio has risen with L. Extrapolations of the box value at dx = 1/32 (`box_extrap.py`), all
estimates and none a bound:

| model | fitted to | limit at dx = 1/32 | sectors − limit |
|---|---|---|---|
| geometric, §6's, the ratio of the last two increments (0.47) | L = 2 to 3.5 | 0.29931 (0.29928 to 0.29935 over the rounding of the prints) | 2.0×10⁻⁴ |
| geometric, the ratio of the first two (0.42) | L = 2 to 3.5 | 0.29928 | 2.3×10⁻⁴ |
| power law ρ∞ − c L⁻ᵃ, a = 3.5 | L = 2.5, 3, 3.5 | 0.29943 | 8×10⁻⁵ |
| power law, a = 3.2 | L = 2, 2.5, 3 | 0.29949 | 2×10⁻⁵ |

So the gap to the box's limit at this spacing is somewhere between 2×10⁻⁵ and 2.3×10⁻⁴ (8×10⁻⁴ relative at most), and
four values do not choose between the models. From L = 2, 2.5, 3 alone the geometric limit was 0.29924, 2.7×10⁻⁴ below
the sectors', so the fourth value moved that estimate by 0.7×10⁻⁴. For σR = 4 and 2 §6 has remaining differences of
1.7×10⁻⁵ and 5×10⁻⁶ after the same kind of extrapolation and a refinement of the spacing at fixed L. No spacing
refinement was run at σR = 8, so the box's spacing error here is not measured, and this note proposes no cause for the
gap. The difference of 3.9×10⁻⁴ at L = 3.5 is set against a box error left at that L of 1.9×10⁻⁴ (geometric) to 3.1×10⁻⁴
and 3.7×10⁻⁴ (power law). That leaves 2.0×10⁻⁴ of the difference unaccounted for by the box error under the
geometric model, and 8×10⁻⁵ or 2×10⁻⁵ under the power law.

### Box comparison, σR = 16: not run

The wall is 1/16 wide, so four points across it need dx ≤ 1/64. At that spacing L = 2 is N = 256, estimated at 31 GiB
(20.69 GiB × (256/224)³), and L = 2.5 is N = 320, estimated at 60 GiB (the same scaling), both above the laptop's 24 GiB
rule and not run. N = 192 (L = 1.5) would fit, at the 12.77 GiB measured for that N, but L = 1.5 is below the L = 2 that
fails `gmres` at σR = 8, where L ≥ 2.5 was needed; whether σR = 16 needs the same L is not measured. So no box value
exists at σR = 16.

The coordinator decided not to run L = 2, N = 256 on a machine with the memory. Step 6's rule applies: a run that does
not fit the 24 GiB rule is reported with its estimated memory, and not run. The estimate is 31 GiB (30.9). It is the
peak working set of the σR = 8 run at N = 224, 20.69 GiB (21,189 MiB, `I96-win/memory-sampler.log`, sampled every 10 s),
times (256/224)³, since the box's memory goes as N³. Scaled the same way from the Mac's maximum resident set at N = 160,
8.37 GiB (`I96-mac/box-s8-L2.5-N160.log`), it is 34 GiB. Neither is a measurement at N = 256.

### Reading notes for §7

1. §5 says no real eigenvalue below the windows is reproduced by the other discretisation. At σR = 8 one is: m = 15/2,
   |μ|v_c = 2.2×10⁻⁸, A's −2.16459×10⁻⁸ against B's −2.16478×10⁻⁸ (relative 8.6×10⁻⁵, absolute 1.9×10⁻¹²). It is below
   every window and at the floor of the search. Nothing here tests why it matches.
2. §5 evaluates the identity on "every real eigenvalue above 10⁻⁸". The script picks those eigenfunctions with |μ| > 10⁻⁸
   in its own unscaled μ, while the real-eigenvalue column counts |μ|v_c > 10⁻⁸. At σR = 8 (v_c = 1.56) the two counts
   differ by 1 to 6 per sector (`reals_check.py 8`: 88 and 87 at m = ½, 103 and 98 at m = 3/2, and so on), the real
   eigenvalues with |μ|v_c from 1.0×10⁻⁸ to 1.6×10⁻⁸ not being evaluated. All are below the windows. The unscaled counts
   equal the table's "(real)" column in all 11 sectors.
3. The pilot's peak, 24.43 GiB, was above the 20 GiB estimated by scaling, and the σR = 16 `conv` runs peaked at 9.78 and
   11.03 GiB against 10 and about 11 GiB estimated. The box run of L = 3.5 peaked at 20.69 GiB against 20.3 GiB. The
   σR = 16 sectors at P = 340 peaked at 28.40 GiB (the worker's resident set) against the 27.6 GiB estimated
   (`I96-dh610/README.md`).
4. The L = 2.5 peak of the Mac probe is 8.37 GiB (8,982,675,456 bytes), which the 2026-11-22 table of #96 gave as 8.98 GB,
   the same bytes divided by 10⁹.
5. The helper `convmax.py` prints 5.1×10⁻⁸ for the first σR = 8 attempt where the `conv` check's line says 5.2×10⁻⁸ and
   7.0×10⁻⁹ for the second σR = 16 attempt against 7.2×10⁻⁹, since it divides the logs' two-digit changes. The note
   quotes the check's lines.
6. A recount of the sector logs from the stored arrays (`I96-rpi5/sector_numbers.py`, R1) agrees with all 132 values of
   the σR = 16 log and with 129 of the σR = 8 log's. The three that differ are in m = 7/2 and 15/2: the log's window of A
   holds 24 and 16 eigenvalues, the recount 26 and 18, and the symmetry value at m = 15/2 is 3.4×10⁻¹¹ in the log and
   5.5×10⁻¹¹ in the recount. In both sectors the log leaves out the two eigenvalues that define the window's cut (one and
   its conjugate), and the recount gives the log's 24, 16 and 3.4×10⁻¹¹ exactly when it leaves them out too. The script
   takes the cut as `abs` of one eigenvalue and the window as `np.abs(array) >= cut`. On the machine that made the recount
   those two differ in the last bit for 309,576 of the 825,000 stored eigenvalues, so whether the cut's own eigenvalues
   are in the window is decided by rounding. That the laptop's rounding put them out is inferred from its log and was
   not run there. No real eigenvalue enters or leaves a window, the cuts agree, and the totals above (254 and 284) follow
   the log; with the two included, A's total is 258. The σR = 8 tables are not changed.
7. The planted check's `real_identity` is 1.0×10⁻¹⁶ in the dh610 log and 2.3×10⁻¹⁶ in the laptop's. The planted weight
   does not depend on σR, and `planted` (7.3×10⁻¹⁵) and `shift` (1.6×10⁻¹¹, max |Re μ|/|μ| 7.7×10⁻¹²) are the same in
   both. Both values are at rounding level and the difference is not investigated.
8. The issue asks that ρ(T₁)v_c be "compared with frank_loss.py's box value at the largest affordable L" and the
   difference "set against the box error measured as in §6". Task 3 below says it "agrees with the box value … within the
   box error". The σR = 8 result meets the issue's wording. It does not show Task 3's agreement: 2×10⁻⁵ to 2.0×10⁻⁴ of
   the difference is left by the extrapolation models, and the spacing error of the box is not measured.

### Acceptance criteria

| criterion | σR = 8 | σR = 16 |
|---|---|---|
| `conv` moves the leading eigenvalue by less than 10⁻⁸ | pass at P = 200, N = 160: 1.1×10⁻⁹ (P = 160 fails at 5.2×10⁻⁸) | pass at P = 340, N = 320: 7.2×10⁻⁹ (P = 320 fails at 2.4×10⁻⁸) |
| planted and shift checks pass at the sizes used | pass, 7.3×10⁻¹⁵ and 1.6×10⁻¹¹, at their own fixed sizes, not at P = 200 | pass, 7.3×10⁻¹⁵ and 1.6×10⁻¹¹ (`real_identity` of `planted` 1.0×10⁻¹⁶), at their own fixed sizes, not at P = 340 |
| A and B agree on ρ(T₁)v_c to 10⁻⁶ | pass, 1.4×10⁻¹¹ (at most 6.5×10⁻¹⁰ relative over the 11 leading eigenvalues) | pass, 1.5×10⁻¹⁰ (at most 5.3×10⁻⁹ relative over the 11 leading eigenvalues) |
| ρ(T₁)v_c against the box at the largest affordable L, the difference against the box error | pass as the issue words it: 0.29951080 against 0.29912 at L = 3.5, +3.9×10⁻⁴; box error left 1.9×10⁻⁴ to 3.7×10⁻⁴ by model; unaccounted 2×10⁻⁵ to 2.0×10⁻⁴, not decided at dx = 1/32 (Reading notes 8) | fail, not run: the coordinator decided not to run L = 2, N = 256 at dx = 1/64, estimated at 31 GiB (34 GiB from the Mac), over the 24 GiB rule. No box value |
| windows and real eigenvalues reported, the identity on each | pass: window 0.1045, \|v\| > 9.57 v_c = 14.95; no real eigenvalue in any window of 11 sectors; identity 0 on 135 to 146 eigenfunctions a sector | pass: window 0.1128, \|v\| > 8.87 v_c = 8.74; no real eigenvalue in any window of 11 sectors (186 A and 216 B eigenvalues in them); identity 0 on 122 to 124 eigenfunctions a sector, 74 to 76 of them real |

Nine of the ten cells pass. The box comparison at σR = 16 does not: no box value exists. The σR = 16 sectors exit 0 with
"all checks passed" (`I96-dh610/sectors-s16-p340.log.exit` holds 0), and the σR = 8 ones too.

### Coordinator's review (#96)

Raw `I96-mac2` (negative-mass-raw). Mac, nice 10, at most 2 processes, one thread each.

| step | result |
|---|---|
| rerun `tanh_spectrum.py --part=vc,planted,shift --sigmas=8,16` | exit 0, the logs' values (planted's worst level error 1.1×10⁻¹⁴ here) |
| one mutation | `--part=shift --quick --mutate=lscale`: exit 1, `shift` fails, on target |
| re-solve one sector, σR = 8, m = 7/2, P = 200 | exit 0, 349 s, 7.31 GiB. The sector line equals the laptop's log, and the arrays match the stored ones to 9×10⁻¹⁵ of the leading modulus |
| rerun `I96-rpi5/sector_numbers.py` | exit 1, like rpi5, but the other way round (below) |
| independent check `ind96.py` | exit 1, G1 fails (below); G2 and G3 pass; its 3 mutations fire on target |

`ind96.py` drives #94's B-spline solver (`tanh_spectrum_check.py`, ARPACK). It shares the formulation and the profile
with `tanh_spectrum.py`, not the basis or the solver. At σR = 8 and 16, m = ½ and 21/2:
- its two discretisations agree to 2×10⁻⁶;
- the leading |μ|v_c agrees with the tables above to 1.7×10⁻⁵ and the second to 4.6×10⁻⁵;
- at σR = 16 the whole difference is that checker's v_c, which is 9.7×10⁻⁶ low (0.9850202). mpmath at 30 digits gives
  0.985029856 (`I96-mac2/vc96.txt`), as `tanh_spectrum.py` has. So criterion G1, v_c to 10⁻⁶, fails, and the fault is in
  the checker's 400-point quadrature, not in this section.

Divided by v_c, the leading |μ| agree to 4×10⁻⁷.

Findings, not fixed in the text above:
- Reading note 6's cause is now tested. On the Mac the recount matches all 132 values of the σR = 8 log, and differs from
  the σR = 16 log at m = 15/2 (14 against 12). np.abs(array) − abs(scalar) at the cut there is −1×10⁻¹⁷ where rpi5 has 0,
  and the reverse. So "σR = 16 recounts exactly" holds on rpi5 only. Each window count of a sector can move by 2 with
  the machine. No real eigenvalue is involved.
- The acceptance cells "planted and shift checks pass at the sizes used" are marked pass above. The coordinator counts
  them as not met: the criterion meant the sector sizes, and the checks ran at p = 60. That makes 7 of the 10 cells,
  not 9.

What this section does not show: that the box error accounts for the 3.9×10⁻⁴ at σR = 8, any box value at σR = 16, a
cause of the rise of the eigenfunction residual with m, that `planted` and `shift` pass at P = 200 or 340 (they run at
sizes of their own), or anything beyond m = 21/2, where the leading eigenvalue is assumed to keep falling. The status of
every number is *numerical evidence*. No novelty is claimed.

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
  zero-mode speed |v| up to 12.05 v_c and 10.32 v_c, given that the leading eigenvalue keeps falling with m (computed
  through m = 41/2 in #120). That is
  two and a half to three times #63's 4 v_c. It is *numerical evidence*. That E is empty for the tanh profile is still
  a *conjecture*, and the window does not reach large v.
- At the thin walls σR = 8 and 16 (§7), ρ(T₁)v_c = 0.29951080 and 0.22508330, and no real eigenvalue lies in a window up to
  9.57 v_c and 8.87 v_c, under the same condition on m > 21/2, there assumed. The box comparison is not closed at σR = 8
  (3.9×10⁻⁴ at L = 3.5) and absent at σR = 16, and the planted and shift checks ran only at their own sizes.
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

Done in #96 (§7). Not met: agreement with the box within its error at σR = 8, any box value at σR = 16, and the
planted checks at the sizes used.


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

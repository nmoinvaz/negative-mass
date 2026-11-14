# A search for real zero modes of the Alcubierre Witten operator, up to 10 v_c

Issue #69. Script [zero_modes.py](zero_modes.py), output [zero_modes_output.txt](zero_modes_output.txt), raw logs and
jsonl in [negative-mass-raw](https://github.com/nmoinvaz/negative-mass-raw) `I69-dh610/` at commit c48ec66. Run from
the repository root; the commands that wrote the output are at the top of each section of it. `python
notes/warp/zero_modes.py --controls --raw=../negative-mass-raw/I69-dh610` runs every `--mutate=` control (`half`,
`conj`, `adj`, `bscale`, `rot`, `lip`, `lyscale`, `lowiter`, `lipbig`, `mem`, `gridperturb`) and exits 1 unless each
fires. Follows [frank_loss.md](frank_loss.md) (#63).

## Answer

No real zero mode of D̃ = D + ½(tr K)γ⁰ was found for the tanh profile at σR = 2 and 4 for 0 < v ≤ 10 v_c. Status:
*numerical evidence*. That E is empty for these data stays a *conjecture*, and nothing here reaches beyond 10 v_c.

The smallest singular value s_min(v) of D̃D⁻¹ falls smoothly from 1 at v = 0 to 0.236 (σR = 2) and 0.275 (σR = 4) at
10 v_c, at the finest grid of the larger box, with no dip. Every one of the 40 steps of 0.25 v_c is certified by a
Lipschitz bound (§1), the smallest lower bound being 0.128 and 0.176, so the grid operator has no real zero mode in the
range. No minimum fell below 10⁻³, so none was refined. Peak memory was 2.27 GiB for one process in the scans and
4.35 GiB at N = 160, under the 8 GB limit.

| σR | v/v_c | L = 8 or 5, N = 64 | L = 8 or 5, N = 96 | L = 16 or 10, N = 96 | L = 16 or 10, N = 128 | 1/L extrapolation |
|---|---|---|---|---|---|---|
| 2 | 1 | 0.691615 | 0.691615 | 0.673207 | 0.673208 | 0.654801 |
| 2 | 2 | 0.531532 | 0.531532 | 0.508119 | 0.508140 | 0.484747 |
| 2 | 4 | 0.389574 | 0.389576 | 0.366724 | 0.366742 | 0.343910 |
| 2 | 6 | 0.322856 | 0.322861 | 0.301986 | 0.302107 | 0.281358 |
| 2 | 8 | 0.282046 | 0.282056 | 0.263178 | 0.263077 | 0.244107 |
| 2 | 10 | 0.253711 | 0.253712 | 0.237040 | 0.236173 | 0.218635 |
| 4 | 1 | 0.729964 | 0.729965 | 0.704532 | 0.704553 | 0.679141 |
| 4 | 2 | 0.580202 | 0.580201 | 0.546105 | 0.546166 | 0.512129 |
| 4 | 4 | 0.440720 | 0.440730 | 0.405416 | 0.405385 | 0.370050 |
| 4 | 6 | 0.374230 | 0.374240 | 0.341144 | 0.341100 | 0.307971 |
| 4 | 8 | 0.333125 | 0.333129 | 0.302212 | 0.302227 | 0.271328 |
| 4 | 10 | 0.304150 | 0.304151 | 0.275082 | 0.275173 | 0.246195 |

The full table, at every 0.25 v_c, is in the output. v_c is the contraction speed of #58, 3.74490 at σR = 2 and 2.47097
at σR = 4. The boxes are #58's (L = 8, 5) and twice them. The last column is 2 s(2L, N₂) − s(L, N₁) at equal spacing
(σR = 2: dx = 0.25, σR = 4: dx = 0.156), which removes a 1/L term. It is an extrapolation from two boxes, not a bound.
It assumes s_min = s_∞ + c/L, the form #63 saw at three boxes ([frank_loss.md](frank_loss.md) §9).

Before any real data the same code recovered planted truth (§2): the zero modes of τ = 6/(1 + r²) at v = ±i, no real
zero mode for that τ, and a real zero mode of a Loss–Yau potential at v = 1.

## 1. What is scanned

On the periodic box [−L, L)³, with #58's FFT discretisation, A(v) := D̃D⁻¹ = I + vM and M := BD⁻¹, with D⁻¹ the symbol
−iγ·k/|k|² and k = 0 set to zero. In the ℓ² norm the singular values of A are those of D̃ : Ḣ¹ → L² (#63 §9). M is
compact in the continuum (#63 §6), so A(v) is Fredholm of index 0 and D̃(v) has a kernel iff s_min(v) := σ_min(A(v)) = 0,
iff −1/v is an eigenvalue of M. On the grid, s_min(v) = 0 exactly when −1/v is an eigenvalue of the matrix M. So a real
zero mode is a dip of s_min(v) to 0.

**The operator.** #63 §2 reduces D̃ to two 2-spinor operators, σ·∇ψ = ½vτ₁ψ and σ·∇χ = −½vτ₁χ (τ₁ = −f′ cos θ). D̃D⁻¹ is
unitarily the direct sum of A(v) and A(−v) for the 2-spinor operator A(v) = I − ½vτ₁(σ·∇)⁻¹. The two blocks are
unitarily equivalent by the π rotation about the y axis, which flips τ₁. The scan uses the one block, which halves the
cost. `sectors` checks that the 4-spinor s_min, the block at +v and the block at −v agree to 3×10⁻¹⁴ at σR = 2, L = 8,
N = 64, v = 2 and 4 v_c. The 4-spinor operator is also in the code (`witten4`), with #63's gammas.

**Full space, not mean-free.** #63's `smin_Dt` restricts to mean-free fields. A periodic box cannot satisfy the k = 0
mode of Bψ, which lifts s_min near a zero. Here s_min is taken on the full space, whose zeros are exactly the
eigenvectors of M. For the tanh profile the two agree to the 6 digits of #63's table (`regress63`, difference from #63's
printed values 4.9×10⁻⁷). They differ for a slowly decaying τ (§2).

**Lanczos.** s_min² is the smallest eigenvalue of AᴴA = (I + v̄Mᴴ)(I + vM), with Mᴴ = −D⁻¹Bᴴ. It is found by ARPACK
(`eigsh`, ncv = 24, relative tolerance 10⁻⁶) from a random starting vector, which has a component in every angular
sector of the rotation symmetry about the x axis. Three starting vectors agree to 4×10⁻¹⁴ at v = 3, 7, 10 v_c, and every
run in every scan converged (`lanczos`, `lanczos_converged`).

**Certificate.** A(v) − A(v′) = (v − v′)M, and singular values are 1-Lipschitz, so |s_min(v) − s_min(v′)| ≤ ‖M‖|v − v′|.
On [v_a, v_b] this gives s_min ≥ (s_a + s_b − ‖M‖(v_b − v_a))/2. ‖M‖ is estimated by Lanczos from below, and 1.02‖M‖ is
used. The scan evaluates s_min at v = 0.25, 0.5, …, 10 v_c, certifies each interval by this bound, bisects any interval
where it is not positive (to 1/16 of the step), and minimises s_min over what remains. A minimum below 10⁻³ is a
candidate, and a candidate is refined at each of the four (L, N) pairs and kept only if its speed settles and s_min
stays below 10⁻³ (`confirm_zero`). The certificate rests on the estimated norm and on the Lanczos value, so it is
*numerical evidence* about the grid operator, not a proof about the continuum. ‖M‖v_c is 0.81 and 0.86 at σR = 2
(L = 8, 16), 0.73 and 0.79 at σR = 4 (L = 5, 10), under the contraction bound 1.

## 2. Planted truth, run before the real data

All checks of this section passed in one run (`python notes/warp/zero_modes.py --confirm`, exit 0).

| check | what it tests | result |
|---|---|---|
| `adjoint` | Mᴴ = −D⁻¹Bᴴ for witten4, witten2 and the Loss–Yau operator, by dot test | 2.6×10⁻¹⁷ |
| `dense` | Lanczos s_min against a dense SVD from explicit DFT matrices, Pauli and γ matrices and a pseudo-inverse of D, N = 8, real and complex v, full and mean-free | 1.9×10⁻¹⁵ |
| `vc` | v_c equal to #58's printed values | 3.7×10⁻⁶ |
| `regress63` | the mean-free s_min at N = 64 against #63's scan table, six points | 4.9×10⁻⁷ |
| `sectors` | 4-spinor and 2-spinor blocks at ±v | 2.7×10⁻¹⁴ |
| `lanczos` | three starting vectors, convergence | spread 3.8×10⁻¹⁴, 0 of 39 unconverged |
| `planted_pm_i` | τ₁ = 6/(1 + r²), L = 10, N = 64 | see below |
| `planted_conv` | the zero-mode speed against the box | see below |
| `planted_real` | the real-v scan of the planted τ | see below |
| `ly_real`, `ly_converged` | a real zero mode, planted | see below |
| `memory` | peak under 8 GB | 0.77 GiB |

The dense check shares with the FFT code only the array τ₁ and the definition of the spectral derivative (k at the
Nyquist index included). Everything else, the DFT as a matrix, the gammas as explicit arrays, the inverse of D, the
SVD and the mean-free projector, is built separately.

**τ = 6/(1 + r²) has its zero modes at v = ±i** (the FL appendix pair, #63 §4). The grid operator at L = 10, N = 64 has
eigenvalues ±0.9910i of M, fourfold in the 2-spinor block, so the zero-mode speeds are v* = ±1.0091i. s_min(v*) is
1×10⁻⁸ and 8×10⁻⁹, and a minimisation along the imaginary axis finds s_min = 7×10⁻⁶ at |v| = 1.00912. At the exact
v = ±i, s_min is 0.0071. The 0.9 % shortfall is the box, as in #63: 1/|λ| = 1.0091, 1.0040, 1.0013 at L = 10, 14, 20
(N = 64), tending to 1. For this τ the mean-free s_min cannot see the zero, it is 0.0355 at v*, which is the k = 0 lift.

**Its real-v scan finds no zero mode**, as §8 of #63 requires for a τ of one sign. s_min stays in 0.84–0.94, all 40
intervals are certified and the smallest value is 0.8446 at 10 v_c. The same function, `scan`, runs the tanh data.

**A real zero mode is found.** The Loss–Yau vector potential A (FL Section 1, `loss_yau` of #63) gives an operator
D − ivσ·A with a zero mode at v = 1, anti-Hermitian in B like the Witten one. The real scan to 3 v_c flags the
interval, bisects it and refines s_min = 5×10⁻⁶ at v* = 1.00989 (L = 10, N = 32). The speed settles toward the exact 1
as the box grows, 1.0099 at L = 10 and 1.0043 at L = 14 (N = 32 to 64, spread 6×10⁻³), with s_min < 10⁻⁵ at each. A
second zero near v = 1.70 is also found.

## 3. The scan

Eight runs, `--part=scan --sigma=S --L=L --N=N` at (σR, L, N) in {2, 4} × {(8, 64), (8, 96), (16, 96), (16, 128)} for
σR = 2 and {(5, 64), (5, 96), (10, 96), (10, 128)} for σR = 4. Each has 40 points. In all eight, all 40 intervals are
certified, no bisection was needed, no candidate appeared, and no Lanczos run failed to converge.

| σR | L | N | dx | ‖M‖v_c | smallest s_min | smallest lower bound | Lanczos applications | wall | peak |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 8 | 64 | 0.25 | 0.8108 | 0.253711 | 0.1519 | 5746 | 622 s | 0.35 GiB |
| 2 | 8 | 96 | 0.167 | 0.8108 | 0.253712 | 0.1519 | 6202 | 2288 s | 1.02 GiB |
| 2 | 16 | 96 | 0.333 | 0.8565 | 0.237040 | 0.1292 | 6202 | 2009 s | 1.02 GiB |
| 2 | 16 | 128 | 0.25 | 0.8564 | 0.236173 | 0.1284 | 5986 | 4561 s | 2.27 GiB |
| 4 | 5 | 64 | 0.156 | 0.7276 | 0.304150 | 0.2130 | 6418 | 614 s | 0.36 GiB |
| 4 | 5 | 96 | 0.104 | 0.7276 | 0.304151 | 0.2130 | 6274 | 2300 s | 1.02 GiB |
| 4 | 10 | 96 | 0.208 | 0.7879 | 0.275082 | 0.1761 | 6490 | 2234 s | 1.02 GiB |
| 4 | 10 | 128 | 0.156 | 0.7884 | 0.275173 | 0.1761 | 6466 | 5112 s | 2.27 GiB |

In every run the smallest value is at 10 v_c, the end of the range, so s_min is still falling there. It falls roughly
as v^(−1/2): s_min √(v/v_c) is 0.67 at v_c and 0.75 at 10 v_c for σR = 2, and 0.70 and 0.87 for σR = 4 (last column of
the full table), which rises slowly, so the fall is a little slower than v^(−1/2). If s_min √(v/v_c) stayed at its value at
10 v_c, s_min would reach 0.1 at 56 v_c (σR = 2) and 76 v_c (σR = 4). This is arithmetic on a trend, not a claim about
v > 10 v_c.

**Peak memory.** One process peaked at 2.27 GiB at N = 128 in the scans and at 4.35 GiB at N = 160 (§4), by the
script's `ru_maxrss` (`time` reports 2.37 GB resident for the N = 128 runs). The script labels GiB as GB, as #63's did.
Up to four processes ran at once, so the sum of their peaks was at most about 11 GiB (two N = 160 runs and two control runs), on a machine with 62 GB.

## 4. Box and grid

**The box dominates.** At one spacing, s_min falls as the box doubles, by 6.9 % at σR = 2 and 9.5 % at σR = 4 at
10 v_c (L = 8 to 16 and 5 to 10), by 2.7 % and 3.5 % at v_c, which is the slow tail of the spinor of #63 §9. The 1/L
extrapolation is 0.2186 and 0.2462 at 10 v_c. It agrees with #63's three-box extrapolation (0.53 at v = 6.3 = 1.68 v_c,
σR = 2): the table gives 0.5166 at 1.75 v_c and 0.5546 at 1.5 v_c. The box error does not bring s_min near zero.

**The grid converges at dx ≤ 0.25 and not at dx = 1/3.** At L = 8 and 5, N = 64 to 96 changes s_min by at most 1.0×10⁻⁵
and 1.2×10⁻⁵. At L = 10 (dx 0.208 to 0.156) it changes by 9.1×10⁻⁵. At L = 16 the pair is dx = 1/3 and 0.25, and s_min
changes by 8.7×10⁻⁴ at 10 v_c (1.8×10⁻⁵ at 4 v_c), because the field varies faster at large v and the coarse grid does
not resolve it. `agree_N` (largest change between the two grids of a box under 2×10⁻⁴) therefore FAILS in the table
run, at σR = 2 and L = 16. The threshold was written before the data and is not changed.

The N = 160 points settle which column is converged. At the larger box of each σR, at v/v_c = 4, 7, 10, N = 160
(dx = 0.2 and 0.125) differs from N = 128 by at most 8.3×10⁻⁶ (`agree_N160`, passes), so the N = 128 columns are
converged to about 10⁻⁵ and the (16, 96) column is not. Peak memory at N = 160 was 4.35 GiB.

| σR | L | v/v_c | N = 128 | N = 160 | change |
|---|---|---|---|---|---|
| 2 | 16 | 4 | 0.366742 | 0.366744 | 1.7×10⁻⁶ |
| 2 | 16 | 7 | 0.280552 | 0.280558 | 6.1×10⁻⁶ |
| 2 | 16 | 10 | 0.236173 | 0.236170 | 2.3×10⁻⁶ |
| 4 | 10 | 4 | 0.405385 | 0.405394 | 8.3×10⁻⁶ |
| 4 | 10 | 7 | 0.319663 | 0.319669 | 6.1×10⁻⁶ |
| 4 | 10 | 10 | 0.275173 | 0.275172 | 5.6×10⁻⁷ |

## 5. Controls

Each mutation is run beside the same command without it, and counts as fired only if its target checks pass without
it and fail with it.

| mutation | what it breaks | checks that fail |
|---|---|---|
| `half` | ½ dropped from B in the FFT operator | `dense`, `regress63`, `planted_pm_i` |
| `conj` | Gram operator uses v, not v̄ | `dense` |
| `adj` | the anti-Hermitian B of the Loss–Yau operator treated as Hermitian in Mᴴ | `adjoint` |
| `bscale` | B of the 2-spinor block scaled by 11/10 | `sectors` |
| `rot` | the planted τ scanned along the imaginary axis | `planted_real` |
| `lip` | a tenth of ‖M‖ in the certificate | `ly_real` |
| `lyscale` | Loss–Yau A scaled by 11/10 in the operator | `ly_real` |
| `lowiter` | Lanczos stopped after two restarts | `lanczos` |
| `lipbig` | 50‖M‖ in the certificate | `scan_resolved` |
| `mem` | memory limit set to 0.05 GB | `memory` |
| `gridperturb` | the (2L, 128) column of the table scaled by 1.01 | `agree_N160` |

All eleven fired in the final runs (`controls-half`, `controls-rest`, `controls-gridperturb` sections of the output).
The first run left `conj` and `lipbig` unfired (Discrepancies, item 3).

`conj` is seen only by `dense`. At a planted zero A φ = 0, so (I + wT)A φ = 0 for any w, and a wrong conjugate leaves
the zero in place, so `planted_pm_i` cannot see it. `dense` at complex v can.

## What this does and does not show

- Within 0 < v ≤ 10 v_c, the grid operator of the tanh profile at σR = 2 and 4 has s_min ≥ 0.236 and 0.275 (finest
  grid, larger box), certified ≥ 0.128 and ≥ 0.176 between samples. So it has no real zero mode there (*numerical
  evidence*). With the 1/L extrapolation the smallest value is 0.219 and 0.246.
- That E is empty for the tanh profile remains a *conjecture*. The range is 10 v_c and the smallest s_min is at its
  end. The scan says nothing about v > 10 v_c.
- The result is for the periodic grid operator at two boxes and two grids. A continuum zero mode decays like r⁻²
  (#63 §3), so it fits a periodic box with a truncation error that falls as the box grows. A dip to 0 would show up as
  a dip at both boxes. This is an argument, not a proof.
- The certificate depends on a Lanczos estimate of ‖M‖ (2 % margin) and on Lanczos for s_min. Both were cross-checked
  against a dense SVD at N = 8 and against #63's table, not at N = 64 by an independent solver.
- The planted tests show the detector sees a real zero mode (Loss–Yau) and a complex one (6/(1 + r²)). The Witten
  operator with a real zero mode and a sign-changing τ was not available, so the real-axis detector was tested on a
  different operator of the same form, an anti-Hermitian B.
- The scan covers the ψ block at v > 0 and relies on the π-rotation symmetry for the χ block and for v < 0. `sectors`
  checks it at two speeds.
- No novelty is claimed.

## Discrepancies and reading notes

1. `agree_N` fails at σR = 2, L = 16 (§4), by 8.7×10⁻⁴ at 10 v_c. It is reported, not loosened, so the table run exits 1.
   The coarse (16, 96) grid has dx = 1/3 and is under-resolved at large v, and N = 160 shows N = 128 is converged
   (`agree_N160`). Its control is `gridperturb` on `agree_N160`, since `agree_N` already fails on the data. Whether the
   (16, 96) column should be dropped from the table is the coordinator's call.
2. #63 took s_min on mean-free fields, this script on the full space. They agree to the printed digits for the tanh
   profile and differ for 6/(1 + r²) at its zero (§2). Nothing in #63 changes.
3. The first control run (`controls-noraw.log` in the raw folder, script commit 7058016) did not fire for two
   mutations. `conj` was aimed at `planted_pm_i` as well as `dense`, but a wrong conjugate cannot move a zero of A. And
   `lipbig` passed because a bisected interval that stayed uncertified was counted as resolved when any candidate had been
   refined. The targets of `conj` and the `scan_resolved` condition were changed, and both controls then fired. The
   scans of σR = 2 at (16, 96), (16, 128) and σR = 4 at (10, 128) ran before the `scan_resolved` change and print its old
   wording, with the same verdict, since all of their intervals were certified.
4. #63 §9 wrote that no real zero-mode speed lies below 3.35 v_c, from the top 24 eigenvalues of T₁. The present range
   is larger and uses s_min, not eigenvalues, and is consistent with it.

## Recommended tasks

1. **Axisymmetric reduction of A(v) to larger boxes and speeds.** The operator conserves the total angular momentum
   about the x axis, so A(v) splits into sectors, each a two-dimensional problem. Compute s_min per sector for L up to
   64, a three-box 1/L extrapolation, and v up to 40 v_c. Why: the box term is the largest error here (§4), and the range
   ends at the smallest value. Acceptance: reproduces this table at L = 8 and 16 to 10⁻⁴, finds the Loss–Yau zero in
   its sector, and a scan to 40 v_c with every interval certified or a candidate refined. Depends on #69.
2. **The large-v trend of s_min.** Fit s_min(v) beyond 10 v_c from task 1 and test whether it decays like v^(−1/2),
   with the smallest value at each v, to say whether the decay ever reaches 10⁻³. Acceptance: a fit with its range and
   error bar stated as a spread over fit choices, planted on a synthetic decay with known exponent. Depends on task 1.

## References

- [FL22] R. L. Frank, M. Loss, arXiv:2201.03610, as in [frank_loss.md](frank_loss.md); the appendix pair
  (τ = 6/(1 + r²)) and the Loss–Yau pair (Section 1) were checked symbolically in #63 (`fl_scalar_extremal`,
  `loss_yau`) and used here as planted truth.
- #58 (the contraction speeds v_c, [alcubierre_bounds.py](alcubierre_bounds.py)), #63 ([frank_loss.md](frank_loss.md)).

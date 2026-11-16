# An independent B-spline check of #75's sector spectra of S⁻¹M

Issue #94. Script [tanh_spectrum_check.py](tanh_spectrum_check.py), output
[tanh_spectrum_check_output.txt](tanh_spectrum_check_output.txt), raw output in
[negative-mass-raw](https://github.com/nmoinvaz/negative-mass-raw) folder `I94-win` (commit f1c7d58). It rechecks
[tanh_spectrum.md](tanh_spectrum.md) (#75), whose sector tables and windows rest on one basis.

## Answer

The leading |μ|v_c of every sector, m = ½ to 21/2 at σR = 2 and 4, agrees with #75's table to 6.0×10⁻⁶ relative or
better, and no real eigenvalue lies in any window, in either discretisation. The report exits 1, because seven of its 90
checks fail: the two discretisations' leading eigenvalues differ by 1.2×10⁻⁶ to 1.9×10⁻⁶ at σR = 4, m = 9/2 to 21/2,
against the script's 10⁻⁶. This basis converges a narrower window than #75's, so its speed bound is weaker, 5.47 v_c at
σR = 2 and 5.19 v_c at σR = 4 (10⁻⁵ tolerance) against #75's 12.05 v_c and 10.32 v_c. It neither reproduces nor contradicts them. Status:
*numerical evidence*.

| check | sectors | verdict |
|---|---|---|
| v_c against #75's 3.7449 and 2.470971 | 2 | pass |
| planted τ = 6/(1 + r²), levels and multiplicities to 10⁻⁵ | 11 | pass, worst 4.1×10⁻⁶ (m = ½), 1.9×10⁻⁶ (m = 21/2) |
| `--mutate=kterm` (drops (m + ½)B/ρ), m = ½ | 1 | fires, exit 1, worst error 5.0×10⁻¹ |
| `--mutate=half` (doubles the weight ½τ), m = ½ | 1 | fires, exit 1, worst error 1.0 |
| leading \|μ\|v_c against #75 to 10⁻⁵ (`lead`) | 22 | pass, worst 6.0×10⁻⁶ (σR = 4, m = 13/2) |
| A and B leading agree to 10⁻⁶ (`AB`) | 22 | σR = 2 pass (worst 9.8×10⁻⁷), σR = 4 fails at m = 9/2 to 21/2, 1.2×10⁻⁶ to 1.9×10⁻⁶ |
| symmetry μ → −μ in the window to 10⁻⁴ (`sym`) | 22 | pass, worst 4.2×10⁻⁶ |
| no real eigenvalue in the window, in A or B (`noreal`) | 22 | pass |

The seven failures are results, and the threshold was not changed. They are of A against B, not of either against #75,
and #75 prints six digits, so its table cannot resolve a comparison below 5×10⁻⁷/|μ|v_c, 6×10⁻⁶ at the last sector of
σR = 2. The `lead` agreement is therefore at the resolution of #75's printed values. By the A–B differences this basis's
own discretisation error in the leading eigenvalue is up to 10⁻⁶ at σR = 2 and 2×10⁻⁶ at σR = 4. #75 reports its two
discretisations differing by 5×10⁻¹² and 5×10⁻¹⁰ in ρ(T₁)v_c, at m = ½, where here they differ by 1.1×10⁻⁷, so this
basis is two to four orders less accurate in the leading eigenvalue.

## What was run

Machine win (surface-laptop-4, Windows 11, 8 logical CPUs, 32 GB), commit 6f7d199, Python 3.14.10 with numpy 2.5.4 and
scipy 1.18.1 (the pinned requirements.txt), from the repository root. Every process had `OMP_NUM_THREADS=1` and
`OPENBLAS_NUM_THREADS=1`, ran at below-normal priority, and at most four ran at a time. `<raw>` is the folder `I94-win`
of negative-mass-raw, cloned beside this checkout.

    python notes/warp/tanh_spectrum_check.py --part=planted
    python notes/warp/tanh_spectrum_check.py --part=planted --ms=0.5 --mutate=kterm        (exits 1)
    python notes/warp/tanh_spectrum_check.py --part=planted --ms=0.5 --mutate=half         (exits 1)
    python notes/warp/tanh_spectrum_check.py --part=solve --sigma=S --disc=D --out=../negative-mass-raw/I94-win/eig
    python notes/warp/tanh_spectrum_check.py --part=report --out=../negative-mass-raw/I94-win/eig > notes/warp/tanh_spectrum_check_output.txt

The solve ran for (S, D) = (2, A), (2, B), (4, A), (4, B), after a pilot of the largest sector, σR = 4, disc B, m =
21/2, whose file the (4, B) solve reused. The report command, run from the repository root with the raw repository
beside it, reproduces the output file apart from its wall-time line (checked). The exits and times:

| run | exit | wall | peak memory |
|---|---|---|---|
| planted, 11 sectors | 0 | 29 min | 3.7 GB |
| `kterm`, `half` controls | 1, 1 | 54 s, 57 s | 0.8 GB |
| pilot, σR = 4, B, m = 21/2, n = 13674 | 0 | 13 min | 5.0 GB |
| solve σR = 2, A | 0 | 57 min | 3.0 GB |
| solve σR = 2, B | 0 | 101 min | 5.0 GB |
| solve σR = 4, A | 0 | 75 min | 3.0 GB |
| solve σR = 4, B (m = 21/2 from the pilot) | 0 | 114 min | 4.6 GB |
| report | 1 | 2 s | |

The peaks are sampled once a minute (`<raw>/memory-sampler.log`), and the machine's free memory stayed above 12.6 GB.
Discretisation A is h = 0.1 min(2/σR, 1/(1 + m/10)) on [0, 4], then q = 1.15 to Rout = 800, and 24(1 + (m − ½)/20) θ
intervals. B is h × 0.9, q = 1.13, Rout = 1200 and 28(1 + (m − ½)/20) intervals. Matrix sizes n are in the tables.

## The sector tables

From the report, `A` and `B` are the two discretisations. "A vs #75" and "A vs B" are relative differences of the
leading |μ|v_c. A window is the |μ|v_c above the largest one of A or B without a partner in the other within 10⁻⁵ (#75's
tolerance) or 10⁻⁴, and "#75 window" is #75's "converged above". "Eigenvalues in window" counts A's in the 10⁻⁴ window,
and "real" counts the real ones in A / B, a real eigenvalue being one with |Im μ| ≤ 10⁻⁹|μ|. A ✗ marks a failed `AB`
check.

σR = 2
| m | n, A / B | \|μ\|v_c, A | B | #75 | A vs #75 | A vs B | window 10⁻⁵ | window 10⁻⁴ | #75 window | eigenvalues in window | real |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.5 | 4888 / 6540 | 0.409556 | 0.409556 | 0.409556 | 2.0×10⁻⁷ | 1.1×10⁻⁷ | 0.1830 | 0.1427 | 0.0777 | 24 | 0 / 0 |
| 1.5 | 5148 / 6840 | 0.298994 | 0.298994 | 0.298994 | 1.1×10⁻⁶ | 2.8×10⁻⁷ | 0.1577 | 0.1271 | 0.0803 | 24 | 0 / 0 |
| 2.5 | 5616 / 7616 | 0.232982 | 0.232982 | 0.232982 | 1.5×10⁻⁸ | 5.9×10⁻⁷ | 0.1738 | 0.1380 | 0.0751 | 12 | 0 / 0 |
| 3.5 | 6264 / 8250 | 0.190371 | 0.190371 | 0.190371 | 5.3×10⁻⁷ | 6.3×10⁻⁷ | 0.1492 | 0.1034 | 0.0711 | 24 | 0 / 0 |
| 4.5 | 6780 / 9030 | 0.160798 | 0.160798 | 0.160798 | 1.6×10⁻⁶ | 8.0×10⁻⁷ | 0.1306 | 0.1097 | 0.0667 | 12 | 0 / 0 |
| 5.5 | 7254 / 9648 | 0.139127 | 0.139127 | 0.139127 | 8.4×10⁻⁷ | 8.6×10⁻⁷ | 0.1160 | 0.0993 | 0.0688 | 12 | 0 / 0 |
| 6.5 | 7808 / 10286 | 0.122581 | 0.122581 | 0.122581 | 3.2×10⁻⁶ | 8.9×10⁻⁷ | 0.1043 | 0.0907 | 0.0649 | 12 | 0 / 0 |
| 7.5 | 8316 / 11232 | 0.109540 | 0.109540 | 0.109541 | 5.8×10⁻⁶ | 9.8×10⁻⁷ | 0.0948 | 0.0834 | 0.0611 | 12 | 0 / 0 |
| 8.5 | 9100 / 11920 | 0.099002 | 0.099002 | 0.099002 | 3.6×10⁻⁶ | 7.2×10⁻⁷ | 0.0868 | 0.0695 | 0.0578 | 20 | 0 / 0 |
| 9.5 | 9720 / 12936 | 0.090309 | 0.090309 | 0.090309 | 1.1×10⁻⁶ | 7.9×10⁻⁷ | 0.0801 | 0.0651 | 0.0595 | 20 | 0 / 0 |
| 10.5 | 10286 / 13674 | 0.083018 | 0.083018 | 0.083018 | 5.5×10⁻⁶ | 7.8×10⁻⁷ | 0.0743 | 0.0613 | 0.0563 | 20 | 0 / 0 |


σR = 4
| m | n, A / B | \|μ\|v_c, A | B | #75 | A vs #75 | A vs B | window 10⁻⁵ | window 10⁻⁴ | #75 window | eigenvalues in window | real |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.5 | 7124 / 9360 | 0.367657 | 0.367657 | 0.367657 | 6.8×10⁻⁷ | 1.1×10⁻⁷ | 0.1926 | 0.1539 | 0.0969 | 16 | 0 / 0 |
| 1.5 | 7124 / 9360 | 0.288531 | 0.288531 | 0.288531 | 1.5×10⁻⁶ | 3.7×10⁻⁷ | 0.1678 | 0.1386 | 0.0910 | 16 | 0 / 0 |
| 2.5 | 7398 / 9984 | 0.235378 | 0.235378 | 0.235378 | 4.6×10⁻⁷ | 8.3×10⁻⁷ | 0.1826 | 0.1491 | 0.0858 | 8 | 0 / 0 |
| 3.5 | 7946 / 10296 | 0.198480 | 0.198480 | 0.198480 | 1.5×10⁻⁶ | 9.2×10⁻⁷ | 0.1600 | 0.1339 | 0.0806 | 8 | 0 / 0 |
| 4.5 | 8220 / 10920 | 0.171512 | 0.171512 | 0.171512 | 2.5×10⁻⁶ | 1.2×10⁻⁶ ✗ | 0.1423 | 0.1214 | 0.0840 | 8 | 0 / 0 |
| 5.5 | 8494 / 11232 | 0.150969 | 0.150969 | 0.150969 | 1.9×10⁻⁷ | 1.4×10⁻⁶ ✗ | 0.1280 | 0.1110 | 0.0790 | 8 | 0 / 0 |
| 6.5 | 8768 / 11544 | 0.134808 | 0.134808 | 0.134809 | 6.0×10⁻⁶ | 1.5×10⁻⁶ ✗ | 0.1163 | 0.1021 | 0.0745 | 8 | 0 / 0 |
| 7.5 | 9042 / 12168 | 0.121766 | 0.121766 | 0.121766 | 3.3×10⁻⁶ | 1.8×10⁻⁶ ✗ | 0.1065 | 0.0946 | 0.0705 | 8 | 0 / 0 |
| 8.5 | 9590 / 12480 | 0.111020 | 0.111020 | 0.111020 | 2.3×10⁻⁶ | 1.6×10⁻⁶ ✗ | 0.0982 | 0.0797 | 0.0727 | 16 | 0 / 0 |
| 9.5 | 9864 / 13104 | 0.102014 | 0.102014 | 0.102014 | 1.4×10⁻⁶ | 1.9×10⁻⁶ ✗ | 0.0912 | 0.0750 | 0.0688 | 16 | 0 / 0 |
| 10.5 | 10286 / 13674 | 0.094358 | 0.094358 | 0.094358 | 4.4×10⁻⁶ | 1.9×10⁻⁶ ✗ | 0.0850 | 0.0708 | 0.0653 | 16 | 0 / 0 |

No real eigenvalue lies in any window of either discretisation, so the listing the issue asks for, of each real
eigenvalue in a window with both values, is empty.

## The speed bound

The report's last lines, the windows taken together over the sectors computed:

    σR = 2, tolerance 1e-05: every eigenvalue with |mu| v_c >= 0.1830 is converged ... so no real zero mode with |v| <= 5.47 v_c
    σR = 2, tolerance 1e-04: every eigenvalue with |mu| v_c >= 0.1427 is converged ... so no real zero mode with |v| <= 7.01 v_c
    σR = 4, tolerance 1e-05: every eigenvalue with |mu| v_c >= 0.1926 is converged ... so no real zero mode with |v| <= 5.19 v_c
    σR = 4, tolerance 1e-04: every eigenvalue with |mu| v_c >= 0.1539 is converged ... so no real zero mode with |v| <= 6.50 v_c

| σR | this check, 10⁻⁵ | this check, 10⁻⁴ | #75 |
|---|---|---|---|
| 2 | 5.47 v_c | 7.01 v_c | 12.05 v_c |
| 4 | 5.19 v_c | 6.50 v_c | 10.32 v_c |

The bound is set by m = ½ at both σR and both tolerances, where the window starts at 0.183 and 0.193, against #75's
0.078 and 0.097. In all 22 sectors, at both tolerances, this basis's window starts above #75's, by a factor of 2.0 to
2.4 at m = ½ and 1.3 at m = 21/2 (10⁻⁵). Both lines are conditional on the leading eigenvalue continuing to fall with m
beyond 21/2, which is not computed, and are statements about the Galerkin eigenvalues, not bounds on the spectrum. At
10⁻⁵ only m = 19/2 and 21/2 have windows below 1/12.05 = 0.0830 (σR = 2) and 1/10.32 = 0.0969 (σR = 4), so this check
confirms the absence of real eigenvalues only above 0.183 and 0.193 in |μ|v_c, and does not test #75's statement below
them.

## What the script shares with the code under check

From its docstring. Shared with [tanh_spectrum.py](tanh_spectrum.py): the formulation, a Galerkin compression of S⁻¹M in
Ḣ¹, and the profile f and τ₁ = −f′cos θ. Not shared: the basis (tensor cubic B-splines in (r, θ) against spin-weighted
angular functions and mapped radial functions), the coordinates, the reduction to a sector, and the solver. Nothing is
shared with [frank_loss.py](frank_loss.py) beyond the profile. The script grew from the review script in #75's closing
comment. Also shared, and recorded here as such:

| shared | what happens to it |
|---|---|
| #75's tables, typed into the script as `REF` (v_c, leading \|μ\|v_c, "converged above") | compared against, not used in the solver |
| the planted weight 6/(1 + r²) and its levels μ = ±i/ℓ_n, ℓ_n = 1 + 2n/3 | the same as #75's `planted` check, which took them from review_fl.md |
| S = 3(π/2)^{4/3} and v_c = 2√S/‖τ₁‖₃ of mass_bounds.md item 5 | the same closed form, by a 400-point Gauss–Legendre quadrature on [0, 15] |
| the "converged above" definition, a partner in the other discretisation within 10⁻⁵ | #75's, with 10⁻⁴ beside it |

The planted weight is spherical, so it tests the sector operator, including the (m + ½)B/ρ term that `kterm` removes,
and not the cos θ dependence of τ₁ that couples the angular functions. That dependence is tested only by the agreement
with #75's tables, so an error in τ₁ common to both would not be seen. The two discretisations A and B are one family, the
same splines, the same Dirichlet condition at Rout and the same Galerkin form, so their agreement shows stability
against the refinements made, not that the family is right. Both this check and #75 assume the reduction of D̃ to σ·∇ψ =
(v/2)τψ, and a wrong reduction or sign convention in mass_bounds.md would not be seen by either.

## Where it differs from tanh_spectrum.md

| | tanh_spectrum.md | here |
|---|---|---|
| leading \|μ\|v_c, 22 sectors | its §4 tables | agree to 6.0×10⁻⁶, at the resolution of its six digits |
| leading eigenvalue, A against B, m = ½ | 5×10⁻¹² (σR = 2), 5×10⁻¹⁰ (σR = 4) | 1.1×10⁻⁷ at both |
| leading eigenvalue, A against B, any sector | not printed | up to 9.8×10⁻⁷ (σR = 2), 1.9×10⁻⁶ (σR = 4) |
| window, m = ½ | 0.0777, 0.0969 | 0.1830, 0.1926 at 10⁻⁵ |
| window, m = 21/2 | 0.0563, 0.0653 | 0.0743, 0.0850 at 10⁻⁵ |
| speed bound | 12.05 v_c, 10.32 v_c | 5.47 v_c, 5.19 v_c at 10⁻⁵ |
| real eigenvalues in a window | none, A or B | none, A or B |
| symmetry μ → −μ over the window | at most 1.7×10⁻⁶ of \|μ\| | at most 4.2×10⁻⁶, over a different set of eigenvalues |
| cost | 36 min, 4 workers, 1.7 and 2.6 GiB each | 114 min the longest solve, 3.0 to 5.0 GB, dense n up to 13674 |

Nothing in the table is a disagreement about the physics. The leading eigenvalues, the absence of real eigenvalues in
the windows and v_c agree. The narrower windows and the larger A–B differences are properties of this basis, as the
script's docstring expects for the windows, and the failed `AB` checks are the one place the script's own threshold was
not met.

The coordinator's review (#94's closing comment) found the cause of the `AB` failures with a third, finer discretisation
at σR = 4, m = 9/2, `--part=solve --sigma=4 --disc=A --ms=4.5 --hfac=0.75 --nth=34`. Its leading |μ|v_c is
0.171511885, with B 4.9×10⁻⁷ and A 1.7×10⁻⁶ below it, so the eigenvalue converges from below. The failures are A's own
discretisation error, about 2×10⁻⁶. The threshold stays as written.

## What this does and does not show

Status words of AGENTS.md. These are *numerical evidence*. A Galerkin eigenvalue of a non-normal operator is not a bound
on the spectrum, nor is a window, and the two discretisations agreeing is the evidence of its stability, not a proof of
it. Nothing here is *proved*, *certified* or *checked symbolically*, and "no real zero mode below |v| ≈ 5 v_c" is
*conditional* on the leading eigenvalue's decline with m. That E is empty for the tanh profile is not shown by this
check either, and stays a *conjecture*.

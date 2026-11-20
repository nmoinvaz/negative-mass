# An independent B-spline check of #75's sector spectra of S⁻¹M

Issue #94. Script [tanh_spectrum_check.py](tanh_spectrum_check.py), output
[tanh_spectrum_check_output.txt](tanh_spectrum_check_output.txt), raw output in
[negative-mass-raw](https://github.com/nmoinvaz/negative-mass-raw) folder `I94-win` (commit f1c7d58). It rechecks
[tanh_spectrum.md](tanh_spectrum.md) (#75), whose sector tables and windows rest on one basis.

Issue #111 widens the windows and reports the second eigenvalue, in [Wider windows (#111)](#wider-windows-111), with
output [tanh_spectrum_check_wide_output.txt](tanh_spectrum_check_wide_output.txt) and raw folder `I111-win`. The other
sections are #94's, from the script at 6f7d199. The regenerated output file has #111's added columns.

## Answer

The leading |μ|v_c of every sector, m = ½ to 21/2 at σR = 2 and 4, agrees with #75's table to 6.0×10⁻⁶ relative or
better, and no real eigenvalue lies in any window, in either discretisation. The report exits 1, because seven of its 90
checks fail: the two discretisations' leading eigenvalues differ by 1.2×10⁻⁶ to 1.9×10⁻⁶ at σR = 4, m = 9/2 to 21/2,
against the script's 10⁻⁶. This basis converges a narrower window than #75's, so its speed bound is weaker at #94's sizes, 5.47 v_c at
σR = 2 and 5.19 v_c at σR = 4 (10⁻⁵ tolerance, widened in #111 below) against #75's 12.05 v_c and 10.32 v_c. It neither reproduces nor contradicts them. Status:
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

#111 widened the windows with finer grids and an ARPACK solver. At m = ½ to 17/2 the 10⁻⁵ window now starts at 0.0552 to
0.0697 (σR = 2) and 0.0562 to 0.0712 (σR = 4), below #75's 0.0830 and 0.0969 in every sector, with m = 19/2 and 21/2 still
from #94's files (0.0801, 0.0743 and 0.0912, 0.0850). No real eigenvalue lies in any window. The second |μ|v_c agrees with
#75's to 4.0×10⁻⁶ (re-solved sectors) and 2.1×10⁻⁵ (#94's two). The speed-bound lines are 12.05 v_c (σR = 2) and
10.60 v_c (σR = 4), against #75's 12.05 v_c and 10.32 v_c. The wide report exits 1 on two `AB` failures, in #94's files.

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

## Wider windows (#111)

Issue #111. Script at 0e3170f (branch `wip/i111-wide-windows`), output
[tanh_spectrum_check_wide_output.txt](tanh_spectrum_check_wide_output.txt), raw output in negative-mass-raw folder
`I111-win` (commit ae8f274). The sectors m = ½ to 17/2 were re-solved at σR = 2 and 4. The sectors 19/2 and 21/2 are
read from #94's `I94-win` files, whose windows already start below the targets. The `files` column of the tables says
which. Status: *numerical evidence*. Galerkin eigenvalues of a non-normal operator are not bounds.

### What was run

Machine win (surface-laptop-4, Windows 11, 8 logical CPUs, 32 GB, shared with choptuik-large-d), Python 3.14.10 with
numpy 2.5.4 and scipy 1.18.1 (the pinned requirements.txt), from the repository root. Every process had
`OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`, ran at below-normal priority and was launched detached with `start`.
At most two solve processes ran at a time. The raw repository is cloned beside this checkout. The sizes are the
coordinator's pilot's, and were not changed to chase a target. The commands, as bash:

    RAW=../negative-mass-raw/I111-win
    REF=../negative-mass-raw/I94-win/eig
    WIDE="--solver=arpack --floor=0.07 --fnth=4 --fh=0.4"
    MS=--ms=0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5
    CHK="python notes/warp/tanh_spectrum_check.py"
    $CHK --part=solve --sigma=2 --disc=A --ms=0.5 --solver=arpack --out=$RAW/base
    $CHK --part=solve --sigma=4 --disc=B --ms=0.5 --solver=arpack --out=$RAW/base
    $CHK --part=compare --out=$RAW/base --ref=$REF
    $CHK --part=compare --out=$RAW/base --ref=$REF --mutate=compare                 (exits 1)
    $CHK --part=planted $MS $WIDE
    $CHK --part=planted --ms=0.5 $WIDE --mutate=kterm                               (exits 1)
    $CHK --part=solve --sigma=4 --disc=B --ms=8.5 $WIDE --out=$RAW/eig             (the pilot)
    $CHK --part=solve --sigma=4 --disc=B $MS $WIDE --out=$RAW/eig                  (and (4, A), (2, B), (2, A))
    $CHK --part=report --out=$RAW/eig,$REF > notes/warp/tanh_spectrum_check_wide_output.txt
    $CHK --part=report --out=$REF > notes/warp/tanh_spectrum_check_output.txt
    $CHK --part=report --out=$RAW/eig,$REF --mutate=next                            (exits 1)

The solve ran for (S, D) = (4, B) and (4, A), then (2, B) and (2, A). The second pair started as slots freed (2, B at
12:20 when (4, A) exited, (2, A) at 13:32 when (4, B) exited), so two processes ran throughout. (4, B) reused the
pilot's m = 17/2 file. All three reports exit 1. The `--mutate=next` one must. The other two fail #94's `AB` check
(below). The times and peak working sets (sampled once a minute, `$RAW/memory-sampler.log`, in GiB):

| run | exit | wall | peak memory |
|---|---|---|---|
| ARPACK at #94's sizes, σR = 2 A and σR = 4 B, m = ½ | 0, 0 | 6 s, 14 s | |
| `compare` against #94's dense spectra | 0 | 0 s | |
| `compare --mutate=compare` | 1 | 0 s | |
| planted, m = ½ to 17/2, n = 32144 to 66960 | 0 | 56 min | 5.7 GiB |
| planted `--mutate=kterm`, m = ½ | 1 | 33 s | 2.0 GiB |
| pilot, σR = 4, B, m = 17/2, n = 93852 | 0 | 14 min | 7.8 GiB |
| solve σR = 4, B (m = 17/2 from the pilot) | 0 | 172 min | 9.7 GiB |
| solve σR = 4, A | 0 | 99 min | 7.9 GiB |
| solve σR = 2, B | 0 | 119 min | 7.9 GiB |
| solve σR = 2, A | 0 | 72 min | 5.6 GiB |
| report, #94 and #111 files | 1 | 1 s | |
| report, #94 files (regenerated) | 1 | 4 s | |
| report, #94 and #111 files, `--mutate=next` | 1 | 1 s | |

The pilot was inside its stop conditions (3 h, 16 GB). The (4, B) solve peaked higher than the pilot, in a sector other
than m = 17/2. From 13:31 a choptuik-large-d build of a Rust extension held 17 GB of the machine's memory at 13:48, and
the machine's free memory fell to 0.4 GB at the lowest sample. No process was stopped and no solve failed.

### Changes to the script

From its docstring paragraph "Wider windows (#111)" and the two #111 commits, 0d27123 and 0e3170f:

- `--fnth` and `--fh` scale both discretisations' θ intervals and core spacing, so A and B stay distinct. The sizes
  here are `--fnth=4 --fh=0.4`, which at m = ½ give n = 32144 (σR = 2, A) to 67716 (σR = 4, B), and 93852 at the largest.
- `--solver=arpack` computes the largest-modulus eigenvalues of K⁻¹N from a sparse LU of K, with k doubled from 64 until
  the smallest returned |μ|v_c is below `--floor`, and stores each pair's residual. The planted check also requires the
  residual below 10⁻⁸. The dense solver is the default and is unchanged.
- The report treats a truncated spectrum as known only above its smallest modulus, so a window cannot start below the
  larger of A's and B's smallest modulus (with 10⁻⁴ room). `--out` may list several directories, read in order per
  sector. The new columns are the second modulus (the largest |μ|v_c below the leading quadruple's by more than 10⁻⁶
  relative) with #75's, the signs of A − B and the directory a sector came from.
- New checks and mutations: `next` (the second modulus against #75's to 10⁻⁴ relative, failed by `--mutate=next`, which
  shifts #75's by 3×10⁻⁴), `--part=compare` (ARPACK against dense spectra at the same sizes, 10⁻¹⁰ relative, failed by
  `--mutate=compare`, which scales the ARPACK eigenvalues by 1 + 10⁻⁹) and the bound lines' "a window from" text.

### The sector tables

From the wide report. "window set by" is `floor` where the 10⁻⁵ window equals the ARPACK truncation floor, so every
eigenvalue above it has a partner within 10⁻⁵ and this basis's own window starts lower, and `match` where an eigenvalue
without a partner sets it. "A − B" gives the signs of |μ|_A − |μ|_B for the leading and the second modulus. The other
columns are as above, and "A vs #75" is relative.

σR = 2
| m | n, A / B | window 10⁻⁵ | window 10⁻⁴ | #75 window | 10⁻⁵ window set by | real, A / B | \|μ\|v_c, A | A vs #75 | A vs B | second \|μ\|v_c, A | A vs #75 | A − B, lead and second | files |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.5 | 32144/42408 | 0.0680 | 0.0680 | 0.0777 | floor | 0 / 0 | 0.409556 | 1.6×10⁻⁸ | 2.5×10⁻¹⁰ | 0.257195 | 1.4×10⁻⁶ | − − | I111 |
| 1.5 | 35496/47124 | 0.0644 | 0.0644 | 0.0803 | floor | 0 / 0 | 0.298994 | 6.5×10⁻⁷ | 2.9×10⁻¹⁰ | 0.207406 | 3.1×10⁻⁷ | − − | I111 |
| 2.5 | 39590/52080 | 0.0608 | 0.0608 | 0.0751 | floor | 0 / 0 | 0.232982 | 9.2×10⁻⁷ | 3.2×10⁻¹⁰ | 0.173753 | 2.3×10⁻⁶ | − − | I111 |
| 3.5 | 43512/57720 | 0.0572 | 0.0572 | 0.0711 | floor | 0 / 0 | 0.190371 | 5.4×10⁻⁷ | 3.3×10⁻¹⁰ | 0.149201 | 3.2×10⁻⁶ | − − | I111 |
| 4.5 | 47792/63180 | 0.0697 | 0.0697 | 0.0667 | floor | 0 / 0 | 0.160798 | 2.9×10⁻⁶ | 3.1×10⁻¹⁰ | 0.130580 | 1.7×10⁻⁶ | − − | I111 |
| 5.5 | 52514/69090 | 0.0654 | 0.0654 | 0.0688 | floor | 0 / 0 | 0.139127 | 4.8×10⁻⁷ | 2.9×10⁻¹⁰ | 0.116022 | 7.1×10⁻⁷ | − − | I111 |
| 6.5 | 57204/75558 | 0.0616 | 0.0616 | 0.0649 | floor | 0 / 0 | 0.122581 | 1.8×10⁻⁶ | 2.6×10⁻¹⁰ | 0.104349 | 1.4×10⁻⁶ | − − | I111 |
| 7.5 | 62094/81472 | 0.0582 | 0.0582 | 0.0611 | floor | 0 / 0 | 0.109541 | 4.4×10⁻⁶ | 2.4×10⁻¹⁰ | 0.094791 | 1.0×10⁻⁶ | − − | I111 |
| 8.5 | 66960/88480 | 0.0552 | 0.0552 | 0.0578 | floor | 0 / 0 | 0.099002 | 2.5×10⁻⁶ | 2.2×10⁻¹⁰ | 0.086825 | 4.0×10⁻⁶ | − − | I111 |
| 9.5 | 9720/12936 | 0.0801 | 0.0651 | 0.0595 | - | 0 / 0 | 0.090309 | 1.1×10⁻⁶ | 7.9×10⁻⁷ | 0.080086 | 2.0×10⁻⁵ | − − | I94 |
| 10.5 | 10286/13674 | 0.0743 | 0.0613 | 0.0563 | - | 0 / 0 | 0.083018 | 5.5×10⁻⁶ | 7.8×10⁻⁷ | 0.074315 | 1.9×10⁻⁵ | − − | I94 |

σR = 4
| m | n, A / B | window 10⁻⁵ | window 10⁻⁴ | #75 window | 10⁻⁵ window set by | real, A / B | \|μ\|v_c, A | A vs #75 | A vs B | second \|μ\|v_c, A | A vs #75 | A − B, lead and second | files |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.5 | 51548/67716 | 0.0712 | 0.0618 | 0.0969 | match | 0 / 0 | 0.367657 | 8.7×10⁻⁷ | 1.6×10⁻¹⁰ | 0.250083 | 9.7×10⁻⁷ | − − | I111 |
| 1.5 | 53652/70686 | 0.0679 | 0.0593 | 0.0910 | match | 0 / 0 | 0.288531 | 8.9×10⁻⁷ | 2.7×10⁻¹⁰ | 0.211765 | 1.3×10⁻⁶ | − − | I111 |
| 2.5 | 56282/73656 | 0.0648 | 0.0568 | 0.0858 | match | 0 / 0 | 0.235378 | 1.7×10⁻⁶ | 4.2×10⁻¹⁰ | 0.182557 | 5.5×10⁻⁷ | − − | I111 |
| 3.5 | 58386/77220 | 0.0620 | 0.0543 | 0.0806 | match | 0 / 0 | 0.198480 | 9.1×10⁻⁸ | 6.0×10⁻¹⁰ | 0.160023 | 2.0×10⁻⁶ | − − | I111 |
| 4.5 | 61016/80190 | 0.0685 | 0.0685 | 0.0840 | floor | 0 / 0 | 0.171512 | 6.1×10⁻⁷ | 7.7×10⁻¹⁰ | 0.142278 | 9.1×10⁻⁷ | − − | I111 |
| 5.5 | 63646/83754 | 0.0649 | 0.0649 | 0.0790 | floor | 0 / 0 | 0.150969 | 2.0×10⁻⁶ | 9.7×10⁻¹⁰ | 0.128000 | 1.5×10⁻⁶ | − − | I111 |
| 6.5 | 66276/87318 | 0.0617 | 0.0617 | 0.0745 | floor | 0 / 0 | 0.134809 | 3.6×10⁻⁶ | 1.2×10⁻⁹ | 0.116287 | 2.1×10⁻⁶ | − − | I111 |
| 7.5 | 68906/90288 | 0.0588 | 0.0588 | 0.0705 | floor | 0 / 0 | 0.121766 | 5.3×10⁻⁷ | 1.5×10⁻⁹ | 0.106515 | 2.0×10⁻⁶ | − − | I111 |
| 8.5 | 71010/93852 | 0.0562 | 0.0562 | 0.0727 | floor | 0 / 0 | 0.111020 | 4.9×10⁻⁷ | 1.8×10⁻⁹ | 0.098245 | 1.5×10⁻⁶ | − − | I111 |
| 9.5 | 9864/13104 | 0.0912 | 0.0750 | 0.0688 | - | 0 / 0 | 0.102014 | 1.4×10⁻⁶ | 1.9×10⁻⁶ | 0.091156 | 1.6×10⁻⁵ | − − | I94 |
| 10.5 | 10286/13674 | 0.0850 | 0.0708 | 0.0653 | - | 0 / 0 | 0.094358 | 4.4×10⁻⁶ | 1.9×10⁻⁶ | 0.085016 | 2.1×10⁻⁵ | − − | I94 |

The 10⁻⁵ window starts at 0.0552 to 0.0697 (σR = 2) and 0.0562 to 0.0712 (σR = 4) in the nine re-solved sectors, against
0.1830 and 0.1926 at m = ½ in #94's files. It is below #75's 0.0830 and 0.0969 in every sector, including the two from
#94's files (0.0801 and 0.0912 at most). It is below #75's own window in 17 of the 18 re-solved sectors. The exception
is σR = 2, m = 9/2, 0.0697 against 0.0667, where the ARPACK stop (k = 64 returns down to 0.0697, just under the 0.07
floor) set the window, not convergence. At σR = 2 every re-solved window is a `floor`. At σR = 4 the windows of
m = ½ to 7/2 are `match`, 0.0077 to 0.0094 above their floors. No real eigenvalue lies in any window, of A or B. ARPACK
keeps the largest 64 to 128 eigenvalues of each sector, so the real eigenvalues below the windows that #75 lists (30 to
66 per sector, unconverged) are not computed here.

### Each check

| check | sectors | verdict |
|---|---|---|
| `compare`, ARPACK against #94's dense spectra at #94's sizes (σR = 2 A, σR = 4 B, m = ½) | 2 | pass, worst partner distance 1.2×10⁻¹⁴ and 1.6×10⁻¹⁴ (126 and 124 eigenvalues) |
| `compare --mutate=compare` | 2 | fires, exit 1, 1.0×10⁻⁹ |
| planted τ = 6/(1 + r²) with ARPACK at WIDE, levels and multiplicities to 10⁻⁵ | 9 | pass, worst 2.5×10⁻⁷ (m = ½), ARPACK residual 1×10⁻¹⁴ or less |
| `--mutate=kterm`, m = ½ | 1 | fires, exit 1, worst error 5.0×10⁻¹ |
| v_c against #75's 3.7449 and 2.470971 | 2 | pass |
| leading \|μ\|v_c against #75 to 10⁻⁵ (`lead`) | 22 | pass, worst 5.5×10⁻⁶ (σR = 2, m = 21/2), all within #75's printed digits |
| second \|μ\|v_c against #75 to 10⁻⁴ (`next`) | 22 | pass, worst 4.0×10⁻⁶ in the 18 re-solved sectors and 2.1×10⁻⁵ in #94's four |
| `--mutate=next` | 22 | fires in all 22, exit 1 (the two `AB` failures also print) |
| A and B leading agree to 10⁻⁶ (`AB`) | 22 | 20 pass, worst 1.8×10⁻⁹ in the re-solved sectors. Fails at σR = 4, m = 19/2 and 21/2, 1.9×10⁻⁶ each |
| symmetry μ → −μ in the window to 10⁻⁴ (`sym`) | 22 | pass, worst 4.2×10⁻⁶ |
| no real eigenvalue in the window, in A or B (`noreal`) | 22 | pass |

The wide report has 112 checks and exits 1 on the two `AB` failures, both from #94's files, where #94 reported the same
1.9×10⁻⁶ (its review traced it to A's own discretisation error). The threshold stays as written. At the re-solved sizes
the leading eigenvalues of A and B agree to 3.3×10⁻¹⁰ (σR = 2) and 1.8×10⁻⁹ (σR = 4) in every sector, against #94's up to
9.8×10⁻⁷ and 1.9×10⁻⁶, and #75's 5×10⁻¹² and 5×10⁻¹⁰ at m = ½ (here 2.5×10⁻¹⁰ and 1.6×10⁻¹⁰). The `lead` agreement
with #75 is at the resolution of its six printed digits in all 22 sectors (the largest difference is 4.9×10⁻⁷ absolute).

#94's regenerated output has 112 checks, 105 pass and 7 fail, the same seven `AB` failures. Its old checks, table
columns, bound lines (apart from the added "a window from" text) and every other line are unchanged (checked by diff
against the committed file, apart from the wall-time line). The added lines are the 22 `next` checks, all passing, and
the new columns.

### The second modulus and the signs of A − B

The second |μ|v_c of the 18 re-solved sectors differs from #75's by at most 4.0×10⁻⁶ relative (σR = 2, m = 17/2, 0.086825),
which is #75's print resolution in all 18 (absolute difference 5×10⁻⁷ or less to rounding, checked per sector). In #94's
four sectors, m = 19/2 and 21/2, it differs by 1.6×10⁻⁵ to 2.1×10⁻⁵, 1.4×10⁻⁶ to 1.8×10⁻⁶ absolute, beyond the print
resolution, with A below #75's in all four and B at or below it. That is #94's discretisation error in the second
modulus, and the issue's unscripted "within 3×10⁻⁵" is this. All are inside the `next` tolerance of 10⁻⁴.

The sign of |μ|_A − |μ|_B is negative for both the leading and the second modulus in all 22 sectors, at both σR. A, the
coarser discretisation, lies below B in every sector, whether the difference is 10⁻⁶ (#94's files) or 10⁻¹⁰ (the
re-solved ones). This is the direction of convergence from below found in #94's review for the leading eigenvalue, and it
holds for the second modulus too. The sign is that of the difference of the unrounded values, at ARPACK residuals of
10⁻¹⁴ or less.

### The speed bound

The report's last lines, the windows taken together over the sectors computed:

    σR = 2, tolerance 1e-05: every eigenvalue with |mu| v_c >= 0.0830 is converged ... so no real zero mode with |v| <= 12.05 v_c
    σR = 2, tolerance 1e-04: every eigenvalue with |mu| v_c >= 0.0830 is converged ... so no real zero mode with |v| <= 12.05 v_c
    σR = 4, tolerance 1e-05: every eigenvalue with |mu| v_c >= 0.0944 is converged ... so no real zero mode with |v| <= 10.60 v_c
    σR = 4, tolerance 1e-04: every eigenvalue with |mu| v_c >= 0.0944 is converged ... so no real zero mode with |v| <= 10.60 v_c

| σR | #94, 10⁻⁵ | #111, 10⁻⁵ | #111, 10⁻⁴ | #75 |
|---|---|---|---|---|
| 2 | 5.47 v_c | 12.05 v_c | 12.05 v_c | 12.05 v_c |
| 4 | 5.19 v_c | 10.60 v_c | 10.60 v_c | 10.32 v_c |

No window starts above the leading eigenvalue of m = 21/2 (0.0830 and 0.0944) in any sector, so both lines are that
leading eigenvalue's reciprocal, 1/0.083018 and 1/0.094358. The largest window starts are 0.0801 and 0.0912, from #94's
m = 19/2 files, and 0.0697 and 0.0712 in the re-solved sectors. At σR = 2 the line equals #75's by this same definition.
At σR = 4 it is 2.7% above #75's 10.32 v_c, because #75's line was its m = ½ window of 0.0969 and the window here is 0.0712.
Both lines are conditional on the leading eigenvalue continuing to fall with m beyond 21/2, which is not computed, and
are statements about the Galerkin eigenvalues. The widening reached its aim, so no memory or time limit set a smaller
window.

### Where it differs from tanh_spectrum.md

| | tanh_spectrum.md | here |
|---|---|---|
| leading and second \|μ\|v_c | its §4 tables | agree to the resolution of its six digits in the 18 re-solved sectors (second: 2.1×10⁻⁵ in #94's four) |
| windows, "converged above" | 0.056 to 0.080 (σR = 2), 0.065 to 0.097 (σR = 4) | 0.055 to 0.070 and 0.056 to 0.071 in the re-solved sectors, mostly the ARPACK floor |
| window, σR = 2, m = 9/2 | 0.0667 | 0.0697, the ARPACK floor |
| speed bound at σR = 4 | 10.32 v_c | 10.60 v_c, since #75's m = ½ window is 0.0969 |
| leading eigenvalue, A against B, m = ½ | 5×10⁻¹² (σR = 2), 5×10⁻¹⁰ (σR = 4) | 2.5×10⁻¹⁰, 1.6×10⁻¹⁰ |
| real eigenvalues below the windows | 30 to 66 per sector | not computed, ARPACK keeps the largest 64 to 128 |
| real eigenvalues in a window | none | none |
| cost | 36 min, 4 workers | 172 min the longest solve, 5.6 to 9.7 GiB, n up to 93852 |

Nothing in the table is a disagreement about the physics. The window and bound differences follow from the two bases
converging differently and from the truncation floor. The only failed check is `AB` in #94's two coarse sectors.

## What this does and does not show

Status words of AGENTS.md. These are *numerical evidence*. A Galerkin eigenvalue of a non-normal operator is not a bound
on the spectrum, nor is a window, and the two discretisations agreeing is the evidence of its stability, not a proof of
it. Nothing here is *proved*, *certified* or *checked symbolically*, and "no real zero mode below |v| ≈ 5 v_c" (#94's
sizes) and "below 12.05 v_c and 10.60 v_c" (#111's) are *conditional* on the leading eigenvalue's decline with m. The
#111 windows are also truncated at the ARPACK floor in 14 of 18 re-solved sectors, so they say where this basis is known
to agree, not where it stops agreeing. That E is empty for the tanh profile is not shown by this check either, and stays
a *conjecture*.

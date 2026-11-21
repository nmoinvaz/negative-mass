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
#118 re-solved those two sectors at the wide sizes, and the wide report then exits 0 (last subsection of "Wider windows
(#111)").

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
beyond 21/2, which #120 later computed through m = 41/2 (last section), and are statements about the Galerkin eigenvalues, not bounds on the spectrum. At
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
which. #118 re-solved those two at the wide sizes (last subsection below). Status: *numerical evidence*. Galerkin
eigenvalues of a non-normal operator are not bounds.

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
(below). With #118's files the wide report exits 0 (last subsection). The times and peak working sets (sampled once a minute, `$RAW/memory-sampler.log`, in GiB):

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
#94's files (0.0801 and 0.0912 at most, and 0.0499 to 0.0538 re-solved in #118). It is below #75's own window in 17 of
the 18 re-solved sectors. The exception
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
1.9×10⁻⁶ (its review traced it to A's own discretisation error). The threshold stays as written. #118 re-solves the two
sectors at the wide sizes, where `AB` passes (last subsection). At the re-solved sizes
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
modulus, and the issue's unscripted "within 3×10⁻⁵" is this. All are inside the `next` tolerance of 10⁻⁴. At #118's wide
sizes the four differ by 6.1×10⁻⁷ to 5.0×10⁻⁶ (last subsection).

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
m = 19/2 files (re-solved in #118, which moves neither line), and 0.0697 and 0.0712 in the re-solved sectors. At σR = 2 the line equals #75's by this same definition.
At σR = 4 it is 2.7% above #75's 10.32 v_c, because #75's line was its m = ½ window of 0.0969 and the window here is 0.0712.
Both lines are conditional on the leading eigenvalue continuing to fall with m beyond 21/2, computed through m = 41/2 in
#120 (last section), and are statements about the Galerkin eigenvalues. The widening reached its aim, so no memory or time limit set a smaller
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
converging differently and from the truncation floor. The only failed check is `AB` in #94's two coarse sectors, which
#118 re-solves at the wide sizes, where it passes.

### m = 19/2 and 21/2 at the wide sizes (#118)

Issue #118. Script unmodified, at 69f5a61 (branch `wip/i118-wide-m19-21`), regenerated
[tanh_spectrum_check_wide_output.txt](tanh_spectrum_check_wide_output.txt), raw output in negative-mass-raw folder
`I118-win` (commit 566373f). The sectors m = 19/2 and 21/2 were re-solved at σR = 2 and 4, discretisations A and B, at
#111's wide sizes, and the wide report reads them from `I118-win` before #111's `I111-win` files. All 22 sectors of the
wide report are now from the wide sizes (4 from `I118-win`, 18 from `I111-win`, none from `I94-win`). The text above
describes #111's state, in which these two sectors came from #94's files and the wide report exited 1. Its tables are
left as they were. Status: *numerical evidence*.

#### What was run

Machine win, Python, thread settings and below-normal detached launches as in #111 above. The sizes are #111's
(`--fnth=4 --fh=0.4`). One solve process ran at a time, and a second started only while choptuik-large-d's processes
together used under 4 cores (CPU-seconds over 10 s, 1.2 to 1.7 cores at the three launches of a second process). The
commands, as bash:

    RAW=../negative-mass-raw/I118-win
    OLD=../negative-mass-raw/I111-win/eig
    WIDE="--solver=arpack --floor=0.07 --fnth=4 --fh=0.4"
    CHK="python notes/warp/tanh_spectrum_check.py"
    $CHK --part=planted --ms=9.5,10.5 $WIDE
    $CHK --part=solve --sigma=4 --disc=B --ms=10.5 $WIDE --out=$RAW/eig                 (the pilot)
    $CHK --part=solve --sigma=4 --disc=B --ms=9.5,10.5 $WIDE --out=$RAW/eig            (and (4, A), (2, B), (2, A))
    $CHK --part=report --out=$RAW/eig,$OLD > notes/warp/tanh_spectrum_check_wide_output.txt
    $CHK --part=report --out=$RAW/eig,$OLD --mutate=next                               (exits 1)

The planted check ran first and passed before any tanh weight was solved. The pilot ran alone, inside its stop condition
(3 h). (4, B) and (4, A) then started at 16:21, (2, B) at 16:47 when (4, B) exited and (2, A) at 17:05 when (4, A)
exited, so two solves ran together from about 16:22 to 17:30. (4, B) reused the pilot's m = 21/2 file. Every n equals the
size planned from the script's knot rules, and the largest, 102414 (B, m = 21/2, both σR), is 9% above #111's largest
93852. Peak working sets were sampled once a minute (`I118-win/memory-sampler.log`, in GiB):

| run | exit | wall | peak memory |
|---|---|---|---|
| planted, m = 19/2 and 21/2, n = 72240 and 78010 | 0 | 47 min | 6.5 GiB |
| pilot, σR = 4, B, m = 21/2, n = 102414 | 0 | 27 min | 8.7 GiB |
| solve σR = 4, B (m = 21/2 from the pilot) | 0 | 25 min | 8.2 GiB |
| solve σR = 4, A | 0 | 43 min | 6.5 GiB |
| solve σR = 2, B | 0 | 43 min | 8.7 GiB |
| solve σR = 2, A | 0 | 38 min | 6.5 GiB |
| report, #118 and #111 files | 0 | 0 s | |
| report, #118 and #111 files, `--mutate=next` | 1 | 0 s | |

A solved sector took 14 min (σR = 2, B, m = 19/2) to 29 min (σR = 2, B, m = 21/2). The machine's free memory stayed
between 12.9 and 26.7 GiB. No process was stopped and none failed.

#### The four sectors

From the wide report. The columns are as in the tables above, with σR as the first. "window set by" is `floor` in all
four.

| σR | m | n, A / B | window 10⁻⁵ | window 10⁻⁴ | #75 window | 10⁻⁵ window set by | real, A / B | \|μ\|v_c, A | A vs #75 | A vs B | second \|μ\|v_c, A | A vs #75 | A − B, lead and second | files |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 9.5 | 72240/94866 | 0.0524 | 0.0524 | 0.0595 | floor | 0 / 0 | 0.090309 | 2.3×10⁻⁶ | 2.0×10⁻¹⁰ | 0.080088 | 4.6×10⁻⁶ | − − | I118 |
| 2 | 10.5 | 78010/102414 | 0.0499 | 0.0499 | 0.0563 | floor | 0 / 0 | 0.083018 | 4.3×10⁻⁶ | 1.8×10⁻¹⁰ | 0.074316 | 4.0×10⁻⁶ | − − | I118 |
| 4 | 9.5 | 73640/96822 | 0.0538 | 0.0538 | 0.0688 | floor | 0 / 0 | 0.102014 | 1.9×10⁻⁶ | 2.1×10⁻⁹ | 0.091157 | 5.0×10⁻⁶ | − − | I118 |
| 4 | 10.5 | 78010/102414 | 0.0516 | 0.0516 | 0.0653 | floor | 0 / 0 | 0.094358 | 1.2×10⁻⁶ | 2.2×10⁻⁹ | 0.085018 | 6.1×10⁻⁷ | − − | I118 |

Against the same sectors in #94's files, as #111's tables above give them:

| σR | m | n, A / B, #94 | window 10⁻⁵, #94 → here | A vs B, #94 → here | second \|μ\|v_c vs #75, #94 → here |
|---|---|---|---|---|---|
| 2 | 9.5 | 9720/12936 | 0.0801 → 0.0524 | 7.9×10⁻⁷ → 2.0×10⁻¹⁰ | 2.0×10⁻⁵ → 4.6×10⁻⁶ |
| 2 | 10.5 | 10286/13674 | 0.0743 → 0.0499 | 7.8×10⁻⁷ → 1.8×10⁻¹⁰ | 1.9×10⁻⁵ → 4.0×10⁻⁶ |
| 4 | 9.5 | 9864/13104 | 0.0912 → 0.0538 | 1.9×10⁻⁶ → 2.1×10⁻⁹ | 1.6×10⁻⁵ → 5.0×10⁻⁶ |
| 4 | 10.5 | 10286/13674 | 0.0850 → 0.0516 | 1.9×10⁻⁶ → 2.2×10⁻⁹ | 2.1×10⁻⁵ → 6.1×10⁻⁷ |

ARPACK kept 64 eigenvalues per sector, down to |μ|v_c = 0.0524 and 0.0499 (σR = 2) and 0.0538 and 0.0516 (σR = 4), at
residuals of 2×10⁻¹⁴ or less, and 60 of them lie in each window. All four windows are the truncation floor, so this
basis's own windows start lower, as for 14 of #111's 18 (now 18 of 22). The 10⁻⁵ windows of the 22 sectors start at
0.0499 to 0.0697 (σR = 2) and 0.0516 to 0.0712 (σR = 4), and below #75's own window in 21 of 22 (the exception is
σR = 2, m = 9/2, as above). The leading and second |μ|v_c differ from #75's printed values by 5.2×10⁻⁸ to 4.5×10⁻⁷
absolute, inside the half unit of its sixth digit. The sign of |μ|_A − |μ|_B is negative for both, as in the other 18
sectors, now at 10⁻¹⁰ to 10⁻⁹ instead of #94's 10⁻⁶. A and B have equal second moduli to the six printed digits in all
four.

#### Each criterion

| criterion of #118 | verdict |
|---|---|
| planted check at m = 19/2 and 21/2 at the wide sizes, before the tanh weight | pass. Exit 0, worst \|μ\| error 5.5×10⁻¹⁰ and 4.2×10⁻¹⁰, residual 10⁻¹⁴ or less, max \|Re μ\| 2.9×10⁻¹⁷ |
| wide report, #118's files then #111's, exits 0 | met. Exit 0, 112 checks, all pass |
| `AB` at 10⁻⁶ in all 22 sectors | pass. Worst 2.2×10⁻⁹ (σR = 4, m = 21/2). The four new ones are 1.8×10⁻¹⁰ to 2.2×10⁻⁹ |
| `next` within 10⁻⁴ of #75's | pass in all 22. Worst 5.0×10⁻⁶ (σR = 4, m = 19/2), the four new ones 6.1×10⁻⁷ to 5.0×10⁻⁶ |
| no real eigenvalue in any window | pass. 0 / 0 in all 22 sectors, A and B |
| 10⁻⁵ windows of m = 19/2 and 21/2 below 0.0830 (σR = 2) and 0.0969 (σR = 4) | met. 0.0524 and 0.0499, 0.0538 and 0.0516, all `floor` |
| speed-bound lines against #111's 12.05 v_c and 10.60 v_c | 12.05 v_c and 10.60 v_c, unchanged to the printed digits |
| `--mutate=next` fails the check | fires. Exit 1, 22 failed `next` checks and no other |

The other checks of the 22 sectors: `lead` passes with worst 4.4×10⁻⁶ (σR = 2, m = 15/2), and `sym` with worst
3.5×10⁻⁷ (σR = 2, m = ½), the four new ones at 3.1×10⁻⁹ or less.

#### The speed bound

The report's last lines are identical to #111's, the windows taken together over the sectors computed:

    σR = 2, tolerance 1e-05 and 1e-04: ... so no real zero mode with |v| <= 12.05 v_c
    σR = 4, tolerance 1e-05 and 1e-04: ... so no real zero mode with |v| <= 10.60 v_c

They are the reciprocals of the leading eigenvalues of m = 21/2 (1/0.083018 and 1/0.094358), because no window starts above
them. The largest window starts in the 22 sectors are now 0.0697 (σR = 2, m = 9/2) and 0.0712 (σR = 4, m = ½), where they had
been 0.0801 and 0.0912 in #94's m = 19/2 files, which were below the same leading eigenvalues. So re-solving the two sectors
did not move the lines, as expected. Both stay conditional on the leading eigenvalue continuing to fall with m beyond 21/2,
which is not computed here (#120 computes it through m = 41/2, next section), and are statements about the Galerkin
eigenvalues.

#### What differs from the text above and from the previous output

- The regenerated output differs from the one at 69f5a61 in 22 lines: the four table rows (n, windows, the two moduli, the
  differences and the directory column, `I94-win/eig` to `I118-win/eig`), the `next`, `AB`, `sym` and `noreal` check lines
  of those four sectors, the wall-time line and the closing line (`FAILED: AB sR=4 m=9.5, AB sR=4 m=10.5` to
  `all checks passed`). The `lead` check lines, every other row and the speed-bound lines are unchanged.
- #94's regenerated report, `tanh_spectrum_check_output.txt`, is untouched and was not rerun. It reads only #94's files,
  so it still exits 1 on its seven `AB` failures at base sizes, which are grid A's own discretisation error. They are not
  an `AB` failure of the wide sizes.
- The planted check took 47 min and a solve 14 to 29 min, against the issue's estimate of about 20 min for the planted
  check and 20 to 30 min per sector.

## Sectors m = 23/2 to 41/2 (#120)

Issue #120. Script unmodified, at ad62a25 (branch `wip/i120-high-m`), output
[tanh_spectrum_check_high_m_output.txt](tanh_spectrum_check_high_m_output.txt), raw output in negative-mass-raw folder
`I120-win` (commit 7cb59a7). The leading eigenvalue of the ten sectors m = 23/2 to 41/2, at σR = 2 and 4, discretisations
A and B, at #94's base sizes (no `--fnth` or `--fh`) with ARPACK. Status: *numerical evidence*. Galerkin eigenvalues of a
non-normal operator are not bounds.

In all 20 pairs of (σR, m) the leading |μ|v_c of A and B agree to 1.8×10⁻⁶ at most, lie below 0.0830 and 0.0944, and fall
strictly with m from #75's m = 21/2 value. The largest are 0.076814 and 0.087769, at m = 23/2. m|μ|v_c rises at every
step. `--part=highm` exits 0 with 42 checks passing, and `--mutate=bump` fails `fall`.

### What was run

Machine win (surface-laptop-4, Windows 11, 8 logical CPUs, 32 GB, shared with choptuik-large-d), Python 3.14.10 with
numpy 2.5.4 and scipy 1.18.1, from the repository root. Thread settings and below-normal detached launches as in #111
above. Two solve processes ran at a time, the second started only while choptuik-large-d's processes together used under
4 cores (1.19 to 1.51 cores at the three readings). The sizes are #94's, with ARPACK at the default `--floor=0.07`. The
commands, as bash:

    RAW=../negative-mass-raw/I120-win
    MS=--ms=11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5
    CHK="python notes/warp/tanh_spectrum_check.py"
    $CHK --part=planted --ms=11.5,20.5 --solver=arpack
    $CHK --part=solve --sigma=4 --disc=B --ms=20.5 --solver=arpack --out=$RAW/eig            (the pilot)
    $CHK --part=solve --sigma=4 --disc=B $MS --solver=arpack --out=$RAW/eig                (and (4, A), (2, B), (2, A))
    $CHK --part=highm --out=$RAW/eig $MS > notes/warp/tanh_spectrum_check_high_m_output.txt
    $CHK --part=highm --out=$RAW/eig $MS --mutate=bump                                     (exits 1)
    python $RAW/trend.py > $RAW/trend.log                                                  (derived columns, no check)

The planted check ran first and passed before any tanh weight was solved. The pilot ran alone, inside its 1 h stop
condition. (4, B) and (4, A) then started at 18:34 and 18:35, (2, B) and (2, A) at 18:39 as they exited. (4, B) reused the
pilot's m = 41/2 file. Every n equals the size planned from the script's knot rules, and the largest, 23484 (B, m = 41/2),
is under a quarter of #118's 102414. The peak working sets were sampled once a minute (`I120-win/memory-sampler.log`):

| run | exit | wall | peak memory |
|---|---|---|---|
| planted, m = 23/2 and 41/2, n = 10868 and 17836 | 0 | 46 s | not sampled (under 100 MB or short) |
| pilot, σR = 4, B, m = 41/2, n = 23484 | 0 | 76 s | 1.5 GiB |
| solve σR = 4, B (m = 41/2 from the pilot) | 0 | 273 s | 1.4 GiB |
| solve σR = 4, A | 0 | 208 s | 1.1 GiB |
| solve σR = 2, B | 0 | 365 s | 1.5 GiB |
| solve σR = 2, A | 0 | 217 s | 1.1 GiB |
| `highm` | 0 | 0 s | |
| `highm --mutate=bump` | 1 | 0 s | |

A sector took 11 to 79 s, and the planted check, pilot and four solves took 14 min of wall time against the issue's
estimate of 30 min. The machine's free memory stayed between 23.7 and 25.9 GiB. No process was stopped and none died.
ARPACK kept 64 eigenvalues per sector, down to |μ|v_c = 0.0336 to 0.0495, at residuals of 2×10⁻¹⁴ or less. All of them
are below the 0.07 floor, so k was never doubled.

### The sectors

From the report. "A vs B" is the leading values' difference relative to A's, "leading / line" the larger of the two
over 0.0830 or 0.0944 (`I120-win/trend.log`). The coordinator's pilot (m = 23/2 and 41/2) is reproduced to every printed
digit.

| σR | m | n, A / B | leading \|μ\|v_c, A | B | A vs B | m\|μ\|v_c | leading / line |
|---|---|---|---|---|---|---|---|
| 2 | 23/2 | 10868 / 14344 | 0.076814 | 0.076814 | 7.6×10⁻⁷ | 0.8834 | 0.9255 |
| 2 | 25/2 | 11544 / 15548 | 0.071472 | 0.071472 | 8.0×10⁻⁷ | 0.8934 | 0.8611 |
| 2 | 27/2 | 12464 / 16262 | 0.066824 | 0.066824 | 5.9×10⁻⁷ | 0.9021 | 0.8051 |
| 2 | 29/2 | 13104 / 17444 | 0.062744 | 0.062744 | 6.1×10⁻⁷ | 0.9098 | 0.7559 |
| 2 | 31/2 | 13846 / 18300 | 0.059132 | 0.059132 | 5.9×10⁻⁷ | 0.9165 | 0.7124 |
| 2 | 33/2 | 14520 / 19074 | 0.055914 | 0.055914 | 5.7×10⁻⁷ | 0.9226 | 0.6737 |
| 2 | 35/2 | 15210 / 20352 | 0.053027 | 0.053027 | 5.8×10⁻⁷ | 0.9280 | 0.6389 |
| 2 | 37/2 | 16262 / 21276 | 0.050424 | 0.050424 | 4.3×10⁻⁷ | 0.9328 | 0.6075 |
| 2 | 39/2 | 17088 / 22624 | 0.048064 | 0.048064 | 4.4×10⁻⁷ | 0.9373 | 0.5791 |
| 2 | 41/2 | 17836 / 23484 | 0.045915 | 0.045915 | 4.2×10⁻⁷ | 0.9413 | 0.5532 |
| 4 | 23/2 | 10868 / 14344 | 0.087769 | 0.087769 | 1.8×10⁻⁶ | 1.0093 | 0.9298 |
| 4 | 25/2 | 11544 / 15548 | 0.082040 | 0.082040 | 1.7×10⁻⁶ | 1.0255 | 0.8691 |
| 4 | 27/2 | 12464 / 16262 | 0.077012 | 0.077012 | 1.3×10⁻⁶ | 1.0397 | 0.8158 |
| 4 | 29/2 | 13104 / 17444 | 0.072564 | 0.072564 | 1.2×10⁻⁶ | 1.0522 | 0.7687 |
| 4 | 31/2 | 13846 / 18300 | 0.068601 | 0.068602 | 1.1×10⁻⁶ | 1.0633 | 0.7267 |
| 4 | 33/2 | 14520 / 19074 | 0.065049 | 0.065049 | 1.0×10⁻⁶ | 1.0733 | 0.6891 |
| 4 | 35/2 | 15210 / 20352 | 0.061846 | 0.061846 | 9.9×10⁻⁷ | 1.0823 | 0.6552 |
| 4 | 37/2 | 16262 / 21276 | 0.058944 | 0.058944 | 7.9×10⁻⁷ | 1.0905 | 0.6244 |
| 4 | 39/2 | 17088 / 22624 | 0.056302 | 0.056302 | 7.6×10⁻⁷ | 1.0979 | 0.5964 |
| 4 | 41/2 | 17836 / 23484 | 0.053886 | 0.053886 | 7.0×10⁻⁷ | 1.1047 | 0.5708 |

The sign of |μ|_A − |μ|_B is negative in all 20 sectors, as in the sectors above, and 1.8×10⁻⁶ at most. #75's m = 21/2
values, from which `fall` starts, are 0.083018 (m|μ|v_c = 0.8717) and 0.094358 (0.9908).

### Each criterion

| criterion of #120 | verdict |
|---|---|
| planted check at m = 23/2 and 41/2 at the sizes used | pass. Exit 0 in 46 s, worst \|μ\| error 1.8×10⁻⁶ and 6.6×10⁻⁷, residual 2×10⁻¹⁴ and 1×10⁻¹⁴, max \|Re μ\| 3.1×10⁻¹⁷ and 1.9×10⁻¹⁷ |
| `ABhigh`: A and B leading agree to 10⁻⁵ in every sector | pass in all 20. Worst 1.8×10⁻⁶ (σR = 4, m = 23/2), 5.6 times inside the threshold. σR = 2: worst 8.0×10⁻⁷ (m = 25/2) |
| `below`: leading under 0.0830 (σR = 2) and 0.0944 (σR = 4) | pass in all 20. Largest 0.076814 and 0.087769 (m = 23/2), 0.9255 and 0.9298 of the lines |
| `fall`: leading falls strictly with m from #75's m = 21/2 value | pass at both σR. The largest step, the one nearest zero, is −2.15×10⁻³ (σR = 2) and −2.42×10⁻³ (σR = 4), both 39/2 to 41/2 |
| `--mutate=bump` fails `fall` | fires. Exit 1, `fall` fails at both σR (largest step +5.9×10⁻⁴ and +6.9×10⁻⁴). `ABhigh` also fails at m = 33/2 at both σR (6.4×10⁻² and 6.1×10⁻²), since the mutation changes A's value only. Four failed checks, no other |

`below` has no `--mutate` of its own, and `bump` leaves the bumped value under the line. As a control, `--part=highm` on
#118's m = 19/2 and 21/2 files (`I120-win/highm-control-below.log`, read-only) exits 1 and fails `below` for 0.090309
and 0.083018 (σR = 2) and 0.102014 (σR = 4), as it should against 0.0830 and 0.0944.

### The trend of m|μ|v_c

"Step" is the change of m|μ|v_c from the sector before, "exponent" the local d ln|μ| / d ln m between the two. The first
row starts from #75's m = 21/2 value (`I120-win/trend.log`).

| m | step, σR = 2 | exponent, σR = 2 | step, σR = 4 | exponent, σR = 4 |
|---|---|---|---|---|
| 23/2 | +0.01167 | −0.8538 | +0.01858 | −0.7957 |
| 25/2 | +0.01004 | −0.8644 | +0.01615 | −0.8096 |
| 27/2 | +0.00873 | −0.8737 | +0.01416 | −0.8218 |
| 29/2 | +0.00765 | −0.8818 | +0.01252 | −0.8325 |
| 31/2 | +0.00677 | −0.8889 | +0.01115 | −0.8420 |
| 33/2 | +0.00603 | −0.8952 | +0.00999 | −0.8505 |
| 35/2 | +0.00540 | −0.9008 | +0.00900 | −0.8581 |
| 37/2 | +0.00487 | −0.9058 | +0.00815 | −0.8650 |
| 39/2 | +0.00441 | −0.9104 | +0.00742 | −0.8712 |
| 41/2 | +0.00401 | −0.9145 | +0.00678 | −0.8769 |

m|μ|v_c rises at every step, from 0.8717 to 0.9413 at σR = 2 and from 0.9908 to 1.1047 at σR = 4, so the leading value
falls more slowly than 1/m throughout, as in the coordinator's pilot. The steps shrink at every sector and the exponent
moves toward −1 from above, but at m = 41/2 neither has reached it: the steps are still 0.0040 and 0.0068. Nothing is
extrapolated here, and the trend is a description of ten computed sectors.

### The speed bound

The lines are unchanged: 12.05 v_c (σR = 2) and 10.60 v_c (σR = 4) are the reciprocals of the leading values of m = 21/2.
Of the Galerkin eigenvalues of m = 23/2 to 41/2 that ARPACK returned (the 64 of largest modulus per sector), none has
modulus above 0.076814 (σR = 2) or 0.087769 (σR = 4), that is none with |v| below 13.02 v_c or 11.39 v_c, above the lines
at 12.05 and 10.60. That holds for the whole sector if ARPACK missed none, which the planted check supports (all 42
expected eigenvalues returned at m = 23/2 and 41/2). #111's comparison with #94's dense spectra at base sizes (1.6×10⁻¹⁴)
covered m = ½ only. This moves the range over which the
fall of the leading eigenvalue is computed from m ≤ 21/2 to m ≤ 41/2. For m > 41/2 it is still an assumption, which #95's
bound would remove. The text above describes the state before #120, and is left as it was.

### What it does not show, and what differs

- A and B share the script, the formulation, the solver and the profile. Their agreement to 10⁻⁶ tests the two grids
  against each other, not the formulation. #75's spectral basis stops at m = 21/2, so no second code covers m > 21/2, and
  the report has no `lead` check against #75 for these sectors. At m = 21/2 the base sizes agreed with #75 to 6.0×10⁻⁶
  (#94), which is far below the 7% margin to the nearest line.
- The report header prints the σR = 4 line as "1/0.0944 = 10.59 v_c", because the script's line is the rounded 0.0944. The
  speed-bound line is 1/0.094358 = 10.60 v_c. The checks do not change, since the nearest sector is 7.0% below its line.
  The same rounding puts #75's own σR = 2 value at m = 21/2 (0.083018) above the script's 0.0830, which the control above
  shows. Neither affects the 20 sectors here.
- The planted check, pilot and solves took 14 min against the estimate of 30 min, and each solve peaked at 1.1 to 1.5 GiB,
  against 5.6 to 9.7 GiB at the wide sizes of #111 and #118.

## Lower floor (#119)

Issue #119. Script unmodified, at a123433 (branch `wip/i119-floor`), regenerated
[tanh_spectrum_check_wide_output.txt](tanh_spectrum_check_wide_output.txt), raw output in negative-mass-raw folder
`I119-win` (commit a89e38d). The 18 sectors whose 10⁻⁵ window #111 or #118 found equal to the ARPACK floor of 0.07
(σR = 2, m = ½ to 21/2, and σR = 4, m = 9/2 to 21/2, both discretisations) were re-solved at `--floor=0.04`, with the
other sizes unchanged. Status: *numerical evidence*. Galerkin eigenvalues of a non-normal operator are not bounds.

Sixteen of the 18 windows are now `match`: set by an eigenvalue without a partner, no longer by the truncation. Two are
still `floor`: σR = 2, m = 19/2 and 21/2, at 0.0341 and 0.0330, the smallest moduli ARPACK returned there, so this
basis's own windows start lower. The wide report exits 1 on those two `match` checks and passes the other 154 of its 156.
All 22 windows are below #75's own (`own75`), including σR = 2, m = 9/2 (0.0475 against 0.0667), no real eigenvalue lies
in any window, and the speed lines are unchanged.

### What was run

Machine win (surface-laptop-4, Windows 11, 8 logical CPUs, 32 GB, shared with choptuik-large-d), Python 3.14.10 with
numpy 2.5.4 and scipy 1.18.1, from the repository root. Thread settings and below-normal detached launches as in #111
above. At most two solve processes ran at a time, the second started only while choptuik-large-d's processes together
used under 4 cores (0.04 to 3.62 cores at the 18 launches made while another solve ran, and four readings of 4.0 to
5.3 cores at 01:03 to 01:05 held the queue for two minutes). The sizes are #111's and #118's. The commands, as bash:

    RAW=../negative-mass-raw/I119-win
    OLD=../negative-mass-raw/I118-win/eig,../negative-mass-raw/I111-win/eig
    NEW="--solver=arpack --floor=0.04 --fnth=4 --fh=0.4"
    CHK="python notes/warp/tanh_spectrum_check.py"
    $CHK --part=solve --sigma=4 --disc=B --ms=10.5 $NEW --out=$RAW/eig                      (the pilot)
    $CHK --part=solve --sigma=2 --disc=A --ms=4.5 $NEW --out=$RAW/eig                       (stage 1, and B)
    $CHK --part=solve --sigma=4 --disc=B --ms=4.5,...,10.5 $NEW --out=$RAW/eig              (and (2, B), m = 0.5,...,10.5)
    $CHK --part=solve --sigma=2 --disc=A --ms=<one sector> $NEW --out=$RAW/eig              (17 single-sector A jobs)
    $CHK --part=report --out=$RAW/eig,$OLD > notes/warp/tanh_spectrum_check_wide_output.txt  (exits 1)
    $CHK --part=report --out=$RAW/eig,$OLD --mutate=floorup                                 (exits 1)
    $CHK --part=compare --out=$RAW/eig --ref=../negative-mass-raw/I111-win/eig              (and I118-win, both exit 0)
    $CHK --part=compare --out=$RAW/eig --ref=../negative-mass-raw/I111-win/eig --mutate=compare   (exits 1)

The pilot ran alone, inside its stop condition (4 h), from 19:08. Stage 1 (σR = 2, m = 9/2, A and B) followed, then (4, B)
and (2, B) with their full lists as slots freed, (4, B) reusing the pilot's m = 21/2 file and (2, B) stage 1's m = 9/2
file. The 17 A sectors not yet started then ran as single-sector jobs, longest expected first, so that the two slots
finished together and not with (2, A) alone for its last hours. The sector files are the same either way. The last solve
exited at 06:17, 11 h 9 min after the pilot started. Peak working sets were sampled once a minute
(`I119-win/memory-sampler.log`):

| run | exit | wall | peak memory |
|---|---|---|---|
| pilot, σR = 4, B, m = 21/2, n = 102414 | 0 | 3095 s | 9.06 GiB |
| stage 1, σR = 2, m = 9/2, A and B | 0, 0 | 1022 s, 2594 s | 2.98 GiB, 5.50 GiB |
| solve σR = 4, B (6 sectors, m = 21/2 from the pilot) | 0 | 285 min | 9.68 GiB |
| solve σR = 2, B (10 sectors, m = 9/2 from stage 1) | 0 | 351 min | 9.07 GiB |
| 17 single-sector A jobs, n = 32144 to 78010 | 0 in all 17 | 664 to 3127 s | 2.00 to 7.95 GiB |
| report on #118's and #111's files (the "before"), report, report `--mutate=floorup` | 1, 1, 1 | 0 s each | |
| `compare` against #111's and #118's files, `compare --mutate=compare` | 0, 0, 1 | 0 s each | |

All 22 solve jobs exited 0 and none died. The 36 solves took 20.9 h of solve time against 9.0 h for the same 36 at floor
0.07 in #111's and #118's logs (ratio 2.33, `I119-win/cost.log`), and the whole run 11 h 9 min of wall time against the
issue's estimate of 15 to 20 h. The ratio per sector is 1.65 to 3.01 where ARPACK stopped at k = 256 (26 files) and 2.87
to 4.67 at k = 512 (10 files). Free memory stayed between 10.3 and 25.7 GiB. No process was stopped.

### The 18 sectors

"Before" is the same report on #118's and #111's files (floor 0.07, `I119-win/report-old-i118-i111.log`), "after" the
wide report. The windows are the 10⁻⁵ ones, "set by" is `floor` or `match`, "real" the real eigenvalues in the window of A
and B, k the eigenvalues ARPACK returned for A and B, down to the smallest |μ|v_c listed (`I119-win/table.log`). Residuals
are 4.8×10⁻¹⁴ or less in all 36 files, the largest at σR = 4, A, m = 17/2 (the run's report said 4×10⁻¹⁴; corrected in
review).

| σR | m | window 10⁻⁵ before | set by | window 10⁻⁵ after | set by | #75's own | real, A / B | k, A / B | smallest \|μ\|v_c, A / B | n, A / B |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 0.5 | 0.0680 | floor | 0.0632 | match | 0.0777 | 0 / 0 | 512 / 512 | 0.0343 / 0.0343 | 32144/42408 |
| 2 | 1.5 | 0.0644 | floor | 0.0586 | match | 0.0803 | 0 / 0 | 512 / 512 | 0.0332 / 0.0332 | 35496/47124 |
| 2 | 2.5 | 0.0608 | floor | 0.0518 | match | 0.0751 | 0 / 0 | 512 / 512 | 0.0323 / 0.0323 | 39590/52080 |
| 2 | 3.5 | 0.0572 | floor | 0.0496 | match | 0.0711 | 0 / 0 | 512 / 512 | 0.0313 / 0.0313 | 43512/57720 |
| 2 | 4.5 | 0.0697 | floor | 0.0475 | match | 0.0667 | 0 / 0 | 512 / 512 | 0.0304 / 0.0304 | 47792/63180 |
| 2 | 5.5 | 0.0654 | floor | 0.0426 | match | 0.0688 | 0 / 0 | 256 / 256 | 0.0394 / 0.0394 | 52514/69090 |
| 2 | 6.5 | 0.0616 | floor | 0.0410 | match | 0.0649 | 0 / 0 | 256 / 256 | 0.0379 / 0.0379 | 57204/75558 |
| 2 | 7.5 | 0.0582 | floor | 0.0395 | match | 0.0611 | 0 / 0 | 256 / 256 | 0.0366 / 0.0366 | 62094/81472 |
| 2 | 8.5 | 0.0552 | floor | 0.0361 | match | 0.0578 | 0 / 0 | 256 / 256 | 0.0353 / 0.0353 | 66960/88480 |
| 2 | 9.5 | 0.0524 | floor | 0.0341 | floor | 0.0595 | 0 / 0 | 256 / 256 | 0.0341 / 0.0341 | 72240/94866 |
| 2 | 10.5 | 0.0499 | floor | 0.0330 | floor | 0.0563 | 0 / 0 | 256 / 256 | 0.0330 / 0.0330 | 78010/102414 |
| 4 | 4.5 | 0.0685 | floor | 0.0553 | match | 0.0840 | 0 / 0 | 256 / 256 | 0.0384 / 0.0384 | 61016/80190 |
| 4 | 5.5 | 0.0649 | floor | 0.0532 | match | 0.0790 | 0 / 0 | 256 / 256 | 0.0373 / 0.0373 | 63646/83754 |
| 4 | 6.5 | 0.0617 | floor | 0.0512 | match | 0.0745 | 0 / 0 | 256 / 256 | 0.0363 / 0.0363 | 66276/87318 |
| 4 | 7.5 | 0.0588 | floor | 0.0465 | match | 0.0705 | 0 / 0 | 256 / 256 | 0.0353 / 0.0353 | 68906/90288 |
| 4 | 8.5 | 0.0562 | floor | 0.0449 | match | 0.0727 | 0 / 0 | 256 / 256 | 0.0344 / 0.0344 | 71010/93852 |
| 4 | 9.5 | 0.0538 | floor | 0.0435 | match | 0.0688 | 0 / 0 | 256 / 256 | 0.0335 / 0.0335 | 73640/96822 |
| 4 | 10.5 | 0.0516 | floor | 0.0400 | match | 0.0653 | 0 / 0 | 256 / 256 | 0.0327 / 0.0327 | 78010/102414 |

The windows of the 18 are 0.0330 to 0.0632 (σR = 2) and 0.0400 to 0.0553 (σR = 4), against 0.0499 to 0.0697 and 0.0516 to
0.0685 before. A `match` window lies above its floor by 0.0008 at σR = 2, m = 17/2 (0.0361 against 0.0353), by 0.0029 to
0.0032 at σR = 2, m = 11/2 to 15/2, and by 0.0073 or more in the other 13 (to the printed digits), so the thinnest is a
margin of 2%. At σR = 2 the windows fall with m (0.0475, 0.0426, 0.0410, 0.0395, 0.0361 for m = 9/2 to 17/2), and the two
`floor` sectors are the two with the largest m, where k = 256 returned down to 0.0341 and 0.0330 without passing below
the basis's own window. At σR = 4 the window at m = 21/2, 0.0400, is 0.0073 above its floor, 0.0327.

ARPACK returned no real eigenvalue (|Im μ| ≤ 10⁻⁹|μ|, the report's tolerance) in any of the 36 files, down to
|μ|v_c = 0.0304 to 0.0394 (`I119-win/reals.log`). So #75's unconverged real eigenvalues below its windows (30 to 66 per
sector) are not among the ones returned here, and `noreal` passes with none to detect.

### Each criterion

| criterion of #119 | verdict |
|---|---|
| the report labels each window `floor` or `match`, from the script | met. 20 `match` and 2 `floor` of 22: 16 of the 18 re-solved are `match`, the 4 not re-solved (σR = 4, m = ½ to 7/2) are #111's `match` |
| `--mutate=floorup` makes every window `floor` and fails a check | fires. Exit 1, `match` fails in all 22 sectors. `own75` also fails in all 22 (the floor is above #75's windows). 112 of 156 checks pass |
| `compare` against #111's files passes above #111's floor | pass. 28 sector files (124 and 124 eigenvalues above #111's floor in 8, 60 and 60 in 20), worst partner distance 2.4×10⁻¹⁴. Against #118's files: 8 files, 2.1×10⁻¹⁴ |
| `compare --mutate=compare` fails | fires. Exit 1, all 28 files fail at 1.0×10⁻⁹ |
| each re-solved window is `match` and below #75's own, or is reported as it is | 16 are `match` and below #75's own. 2 are `floor` (σR = 2, m = 19/2 and 21/2, 0.0341 and 0.0330), below #75's own 0.0595 and 0.0563, reported as they are. The report exits 1 on their `match` |
| no real eigenvalue lies in any window | pass. `noreal` in all 22 sectors, 0 / 0 |
| `own75`, the 10⁻⁵ window below #75's own | pass in all 22, including σR = 2, m = 9/2 (0.0475 against 0.0667), which failed at #111 (0.0697) |

The other checks of the 22 sectors are the same as before: `lead` worst 4.4×10⁻⁶ (σR = 2, m = 15/2), `next` worst
5.0×10⁻⁶ (σR = 4, m = 19/2), `AB` worst 2.2×10⁻⁹ (σR = 4, m = 21/2), and `sym` worst 3.3×10⁻⁶ (σR = 2, m = ½, against
3.5×10⁻⁷ before, in a wider window now).

### The speed bound

The report's last lines are those of #111 and #118:

    σR = 2, tolerance 1e-05 and 1e-04: ... so no real zero mode with |v| <= 12.05 v_c
    σR = 4, tolerance 1e-05 and 1e-04: ... so no real zero mode with |v| <= 10.60 v_c

They are the reciprocals of the leading eigenvalues of m = 21/2 (1/0.083018 and 1/0.094358), because no window starts
above them. The largest window starts are 0.0632 (σR = 2, m = ½) and 0.0712 (σR = 4, m = ½, #111's file). Lowering the
floor moved neither line, as #120 found: they are set by the leading eigenvalue of m = 21/2, not by a window. Both stay
conditional on the leading eigenvalue continuing to fall with m beyond 21/2, computed through m = 41/2 in #120 (previous
section), and are statements about the Galerkin eigenvalues.

### What differs from the previous output

The regenerated output differs from the one at 41fee23 (main) in 105 added and 61 removed lines. They are the two table
headers (new "set by" column), the 22 table rows (the new column, and in the 18 re-solved rows the windows, the count of
eigenvalues in the window and the directory column), 18 `sym` and 18 `noreal` lines (the new window), 44 new check lines
(22 `match`, 22 `own75`) and the closing line (`all checks passed` to `FAILED: match sR=2 m=9.5, match sR=2 m=10.5`).
The 66 `lead`, `next` and `AB` lines and the speed-bound lines are identical, and in all 22 rows the |μ|v_c columns, the
second moduli, the signs of A − B and the differences from #75 are unchanged to the printed digits. One column changes
that carries no information: the sign of the imaginary part of the leading μ differs in 10 rows, ARPACK's pick within a
conjugate pair. #94's regenerated output is untouched.

### What it does not show, and what differs

- Two windows are still the truncation floor, so for σR = 2, m = 19/2 and 21/2 the basis's own window is unknown and below
  0.0341 and 0.0330. The doubling stops at the first k whose smallest modulus is under `--floor`, so a floor of 0.04 left
  k = 256 there and a lower `--floor` is needed. Nothing is extrapolated from the other sectors.
- The truncation remains in the other 16: the spectrum is known down to 0.0304 to 0.0394, and a `match` window says only
  that an eigenvalue above its floor has no partner within 10⁻⁵ in the other discretisation. The thinnest margin is
  0.0008 (σR = 2, m = 17/2).
- `compare` is the same code at two values of k, so it shares the script, formulation, grids and solver with what it
  checks, and tests that the Krylov dimension changes no eigenvalue above the older floor. It does not test the
  formulation, which the `lead` and `next` checks against #75's spectral basis do.
- `--mutate=floorup` fails `own75` as well as `match`, since the floor lies above #75's windows. `noreal` has no
  mutation of its own that fires it, and with no real eigenvalue returned anywhere it has not seen one here.
- The estimate was 15 to 20 h of wall time, and the run took 11 h 9 min. The pilot's ratio (1.95, k = 256) was below the
  coordinator's Mac pilot's 3.4 (m = ½, k = 512), which lies in the 2.87 to 4.67 found at k = 512 here.

## k = 512 at σR = 2, m = 19/2 and 21/2 (#127)

Issue #127. Script at 4b9a8b9 (branch `wip/i127-k512`, unmodified by this task), regenerated
[tanh_spectrum_check_wide_output.txt](tanh_spectrum_check_wide_output.txt), raw output in negative-mass-raw folder
`I127-win` (commit 3e1ccda). The two windows #119 left at the ARPACK floor (σR = 2, m = 19/2 and 21/2, both
discretisations) were re-solved with `--floor=0.025 --k0=512 --kmax=512`, at #119's sizes. Status: *numerical
evidence*. Galerkin eigenvalues of a non-normal operator are not bounds.

Both windows are now `match`: 0.0332 (m = 19/2) and 0.0322 (m = 21/2) at the 10⁻⁵ tolerance, below #119's floors of
0.0341 and 0.0330 and below #75's own 0.0595 and 0.0563. The wide report exits 0 with 156 of 156 checks. ARPACK stopped
at k = 512 above the requested floor (0.0266 and 0.0259, not under 0.025), and a window set by an eigenvalue above that
floor is what `match` asks for. No real eigenvalue lies in either window or is among the 512 returned, and the speed
lines are unchanged at 12.05 v_c and 10.60 v_c.

### What was run

Machine win (surface-laptop-4, Windows 11, 8 logical CPUs, 32 GB, shared with choptuik-large-d), Python 3.14.10 with
numpy 2.5.4 and scipy 1.18.1, from the repository root. Thread settings and below-normal detached launches as in #111 and
#119 above. At most two solve processes ran at a time, the second started once, while choptuik-large-d's processes
together used 2.46 cores (the threshold is 4). The sizes are #118's. The commands, as bash:

    RAW=../negative-mass-raw/I127-win
    OLD=../negative-mass-raw/I119-win/eig,../negative-mass-raw/I118-win/eig,../negative-mass-raw/I111-win/eig
    NEW="--solver=arpack --floor=0.025 --k0=512 --kmax=512 --fnth=4 --fh=0.4"
    CHK="python notes/warp/tanh_spectrum_check.py"
    $CHK --part=solve --sigma=2 --disc=B --ms=10.5 $NEW --out=$RAW/eig                      (the pilot, alone)
    $CHK --part=solve --sigma=2 --disc=B --ms=9.5 $NEW --out=$RAW/eig
    $CHK --part=solve --sigma=2 --disc=A --ms=9.5,10.5 $NEW --out=$RAW/eig
    $CHK --part=report --out=$RAW/eig,$OLD > notes/warp/tanh_spectrum_check_wide_output.txt  (exits 0)
    $CHK --part=report --out=$RAW/eig,$OLD --mutate=floorup                                 (exits 1)
    $CHK --part=report --out=$RAW/eig,$OLD --mutate=real                                    (exits 1)
    $CHK --part=compare --out=$RAW/eig --ref=../negative-mass-raw/I119-win/eig              (exits 0)
    $CHK --part=compare --out=$RAW/eig --ref=../negative-mass-raw/I119-win/eig --mutate=compare   (exits 1)

The pilot (the largest matrix, n = 102414) ran alone from 06:38 to 07:41, inside its stop condition (5 h). The other two
jobs then ran together from 07:42, B at m = 19/2 until 08:28 and A (m = 19/2, then 21/2) until 09:15. Peak working sets
were sampled once a minute (`I127-win/memory-sampler.log`, `runs.log`):

| run | exit | wall | peak memory |
|---|---|---|---|
| pilot, σR = 2, B, m = 21/2, n = 102414 | 0 | 3817 s | 9.87 GiB |
| σR = 2, B, m = 19/2, n = 94866 | 0 | 2752 s | 6.91 GiB |
| σR = 2, A, m = 19/2 and 21/2, n = 72240 and 78010 | 0 | 2950 s and 2620 s | 7.33 GiB (the job's) |
| report, `floorup`, `real`, `compare`, `compare --mutate=compare` | 0, 1, 1, 0, 1 | 0 s each | |

All three jobs exited 0 and none died. The four solves took 12139 s (3.4 h) of solve time and the whole run 2 h 37 min of
wall time, against the issue's 6 to 9 h and 3 to 5 h. Free memory stayed between 12.1 and 24.9 GiB. No process was
stopped.

### The two sectors

"#119" is the wide report of #119 (k = 256, floor 0.04, `I119-win`), "now" the wide report of this task
(`I127-win/table.log`). Windows are the 10⁻⁵ ones, "set by" is `floor` or `match`, "real" the real eigenvalues in the
window of A and B, k the eigenvalues ARPACK returned for A and B, down to the smallest |μ|v_c listed. Residuals are
3.5×10⁻¹⁴ to 5.0×10⁻¹⁴ in the four files (largest at B, m = 19/2).

| m | window 10⁻⁵, #119 | set by | window 10⁻⁵, now | set by | #75's own | real, A / B | k, A / B | smallest \|μ\|v_c, A / B | n, A / B | time, A / B | peak, A / B |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 19/2 | 0.0341 | floor | 0.0332 | match | 0.0595 | 0 / 0 | 512 / 512 | 0.0266 / 0.0266 | 72240 / 94866 | 2950 s / 2752 s | 7.33 / 6.91 GiB |
| 21/2 | 0.0330 | floor | 0.0322 | match | 0.0563 | 0 / 0 | 512 / 512 | 0.0259 / 0.0259 | 78010 / 102414 | 2620 s / 3817 s | 7.33 / 9.87 GiB |

The 10⁻⁴ windows are 0.0277 and 0.0260 (0.0341 and 0.0330 in #119). The eigenvalues in the window rise from 252 and
252 to 456 and 504. A `match` window lies above its floor by 0.0066 (m = 19/2) and 0.0063 (m = 21/2), 24 to 25 %,
against 2 % at σR = 2, m = 17/2 (0.0361 against 0.0353, still the thinnest). At σR = 2 the windows now fall with m
through all seven sectors from 9/2: 0.0475, 0.0426, 0.0410, 0.0395, 0.0361, 0.0332, 0.0322.

ARPACK returned no real eigenvalue (|Im μ| ≤ 10⁻⁹|μ|, the report's tolerance) in any of the four files, down to
|μ|v_c = 0.0259 and 0.0266 (`I127-win/reals.log`), as in #119's 36.

### Each criterion

| criterion of #127 | verdict |
|---|---|
| each new file holds k = 512 and a residual under 10⁻¹² | met. 4 files, 512 eigenvalues each, residual 3.5×10⁻¹⁴ to 5.0×10⁻¹⁴ |
| both windows `match`, the wide report exits 0 with 156 checks | met. 0.0332 and 0.0322, both `match`, exit 0, 156 of 156 pass. Not set by the floor: ARPACK's smallest moduli 0.0266 and 0.0259 lie below the windows |
| `compare` against #119's files passes, `--mutate=compare` fails | pass: exit 0, 4 files, 252 and 252 eigenvalues each above #119's floor, worst partner distance 2.8×10⁻¹⁴. Fires: exit 1, 4 of 4 fail at 1.0×10⁻⁹ |
| `--mutate=floorup` still fails `match` in 22 of 22 | fires. Exit 1, `match` fails in 22 of 22 (`own75` too, as in #119), 112 of 156 checks pass |
| `--mutate=real` exits 1, `noreal` and `sym` fail at σR = 2, m = 11/2 and nowhere else | fires. Exit 1, those two checks fail and 154 of 156 pass. The unmutated report fails nothing |
| the other 20 rows unchanged, the speed lines against 12.05 v_c and 10.60 v_c | met. 20 rows and their checks are identical, the speed lines are identical (12.05 v_c and 10.60 v_c) |

The other checks are the same as before: `lead` worst 4.4×10⁻⁶ (σR = 2, m = 15/2), `next` worst 5.0×10⁻⁶ (σR = 4,
m = 19/2), `AB` worst 2.2×10⁻⁹ (σR = 4, m = 21/2), `sym` worst 3.3×10⁻⁶ (σR = 2, m = ½). `sym` at the two new sectors is
3.5×10⁻⁷ and 2.7×10⁻⁷, against 4.4×10⁻⁸ and 3.3×10⁻⁸ in the narrower windows of #119.

### The speed bound

The report's last lines are those of #119, #111 and #118:

    σR = 2, tolerance 1e-05 and 1e-04: ... so no real zero mode with |v| <= 12.05 v_c
    σR = 4, tolerance 1e-05 and 1e-04: ... so no real zero mode with |v| <= 10.60 v_c

They are the reciprocals of the leading eigenvalues of m = 21/2 (1/0.083018 and 1/0.094358), because no window starts
above them. Lowering the windows of the last two sectors moved neither line, as #120 found: they are set by the leading
eigenvalue of m = 21/2, not by a window. Both stay conditional on the leading eigenvalue continuing to fall with m beyond
21/2 (computed through m = 41/2 in #120), and are statements about the Galerkin eigenvalues.

### What differs from the previous output

The regenerated output differs from the one at fa3f7b3 (main, #119's) in 11 added and 11 removed lines: the 2 table rows
of σR = 2, m = 19/2 and 21/2 (the 10⁻⁵ and 10⁻⁴ windows, `floor` to `match`, the eigenvalues in the window, and the
directory column `I119-win/eig` to `I127-win/eig`), 8 check lines (`sym`, `noreal`, `match` and `own75` of those two
sectors, `match` going from FAIL to PASS) and the closing line (`FAILED: match sR=2 m=9.5, match sR=2 m=10.5` to `all
checks passed`). The |μ|v_c columns of those two rows, their `lead`, `next` and `AB` lines, the other 20 rows and the
speed-bound lines are unchanged to the printed digits (`I127-win/diff-wide.log`). #94's regenerated output is untouched.

### What it does not show, and what differs

- ARPACK's smallest modulus at k = 512 is 0.0266 and 0.0259, not under the requested 0.025. From k = 256 the smallest
  modulus fell by a factor of 0.78 per doubling (0.0341 to 0.0266 and 0.0330 to 0.0259), where #127 expected about 0.71
  (k = 128 to 512 at σR = 2, m = ½ in #119), and so about 0.024 and 0.023. The sectors are `match` all the same, so no
  further k was needed; had the window equalled the smallest modulus, the report would have said `floor` there.
- The truncation remains: the spectrum is known down to 0.0259 and 0.0266, and a `match` window says only that an
  eigenvalue above its floor has no partner within 10⁻⁵ in the other discretisation. The margins of 0.0066 and 0.0063 are
  wider than the 0.0008 of m = 17/2.
- `compare` is the same code at two values of k, so it shares the script, formulation, grids and solver with what it
  checks, and tests that the Krylov dimension changes no eigenvalue above #119's floor. It does not test the formulation,
  which `lead` and `next` against #75's spectral basis do.
- `--mutate=real` shows that `noreal` fires on a real eigenvalue, here the leading one of σR = 2, m = 11/2 made real with
  its modulus kept (0.1391, far above the window's edge). It tests the check's detection and the pairing of A and B, not
  its sensitivity at the window's lower edge. No real eigenvalue is returned in any file: the four new ones, or the other
  32 current files of #119 (the four it superseded also had none).
- The estimate of 6 to 9 h of solve time was 3.4 h. One run at k = 512 took 0.93 to 1.12 times #119's k = 64 to 256 chain
  of the same sector (2950 s against 2635 s for A, m = 19/2; 2620 s against 2497 s for A, m = 21/2; 2752 s against 2607 s
  for B, m = 19/2; 3817 s against 4094 s for B, m = 21/2, which ran beside another solve there and alone here), not 2 to 3
  times a k = 256 run. The peak of 9.87 GiB (B, m = 21/2) is 0.8 GiB above the 9.07 GiB of #119's σR = 2 B job and 0.2 GiB
  above its largest peak, 9.68 GiB (σR = 4, B), not the 1 to 2 GiB more the issue expected.

## #75's code at m = 23/2 to 41/2 (#125)

Issue #125. Scripts unmodified, at 5ed8160 (branch `wip/i125-highm75`), output
[tanh_spectrum_highm_output.txt](tanh_spectrum_highm_output.txt), raw output in negative-mass-raw folder `I125-win`
(commit 6417f1c). #75's spectral code, `tanh_spectrum.py`, computed the leading |μ|v_c of m = 23/2, 29/2, 35/2 and 41/2
at σR = 2 and 4, and of m = 21/2 as a reproduction, at #75's sizes. `tanh_spectrum_highm.py` compares them with #120's
B-spline values. Status: *numerical evidence*. Galerkin eigenvalues of a non-normal operator are not bounds.

In all eight (σR, m) #75's leading value agrees with #120's, in A and in B, to 3.0×10⁻⁶ at most (the threshold is 10⁻⁴),
and lies above both in all 16 comparisons. #75's own two resolutions agree to 1.6×10⁻⁷ at most in all ten sectors. The
largest values are 0.07681409 and 0.08776927, at m = 23/2. `tanh_spectrum_highm.py` exits 0 with 28 checks passing, and
`--mutate=plant` and `--mutate=conv` each exit 1 with only their target failing.

### What was run

Machine win (surface-laptop-4, Windows 11, 8 logical CPUs, 32 GB, shared with choptuik-large-d), Python 3.14.10 with
numpy 2.5.4 and scipy 1.18.1, from the repository root. Thread settings and below-normal detached launches as in #111
above, with `--procs=1`, so a job is a main process and one worker. Two jobs ran at a time, the second started only
while choptuik-large-d's processes together used under 4 cores (2.02 cores at the one reading). The sizes are #75's, the
script's defaults without `--quick`: A has dK = 28, p = 80 (σR = 2) or 100 (σR = 4), Lmap = 2, npan = 160, and B has dK
= 38, p = 90 or 110, Lmap = 3, npan = 240. The commands, as bash:

    RAW=../negative-mass-raw/I125-win
    TS="python notes/warp/tanh_spectrum.py --part=vc,planted,tanh --procs=1 --raw=$RAW/eig"
    HM="python -X utf8 notes/warp/tanh_spectrum_highm.py --raw=$RAW/eig --ref=../negative-mass-raw/I120-win/eig"
    $TS --sigmas=4 --ms=20.5                                    (the pilot, alone)
    $TS --sigmas=2 --ms=10.5,11.5,14.5,17.5,20.5
    $TS --sigmas=4 --ms=10.5,11.5,14.5,17.5
    $HM > notes/warp/tanh_spectrum_highm_output.txt             (exits 0)
    $HM --mutate=plant                                           (exits 1)
    $HM --mutate=conv                                            (exits 1)

The pilot ran alone from 09:38 to 09:44, inside its 30 min stop condition. The two solves started at 09:45 and ended at
10:04 and 10:11, so the three solving runs took 33 min of wall time against the issue's 1 to 2 h. The peak working sets
were sampled once a minute during the pilot and every 10 s after (`I125-win/memory-sampler.log`, `runs.log`), the
worker's being the job's peak:

| run | exit | wall | peak memory |
|---|---|---|---|
| pilot, σR = 4, m = 41/2 | 0 | 342 s | 2.57 GiB |
| σR = 2, m = 21/2 to 41/2 | 0 | 1126 s | 1.69 GiB |
| σR = 4, m = 21/2 to 35/2 | 0 | 1539 s | 2.57 GiB |
| `tanh_spectrum_highm.py` | 0 | 0.2 s | |
| `--mutate=plant` | 1 | 0.2 s | |
| `--mutate=conv` | 1 | 0.2 s | |

The pilot's last 42 s, B's eigensolve, went unsampled. The σR = 4 solve has the same sizes and peaked at 2.57 GiB (2636
MiB against the pilot's 2627 MiB), so the pilot's value is not short. The machine's free memory stayed between 19.3 and
25.1 GiB. No process was stopped and none died. In each solving run the planted check passed first (worst error
7.3×10⁻¹⁵), `v_c` equals mass_bounds.md's to 3.7×10⁻⁶, and the script's `converged`, `no_real_in_window` and
`real_identity` checks pass. `ref63` is skipped with a printed line, since m = ½ is not among the sectors, and the
"peak resident" lines print nan on Windows.

### The sectors

From `I125-win/table.log`. "A vs B" is #75's two resolutions' difference relative to A's value, "#75 A vs #120" the
difference of #75's A from #120's A and B, relative to each. The time is the sector's, the peak its process's.
m = 21/2 is the reproduction, which #120 did not compute.

| σR | m | #75 A | #75 B | A vs B | #120 A | #120 B | #75 A vs #120 A / B | time, peak of its process |
|---|---|---|---|---|---|---|---|---|
| 2 | 21/2 | 0.08301764 | 0.08301764 | 3.1e-12 | — | — | — | 198 s, 1.69 GiB |
| 2 | 23/2 | 0.07681409 | 0.07681409 | 5.8e-12 | 0.07681400 | 0.07681406 | 1.2e-06 / 3.9e-07 | 201 s, 1.69 GiB |
| 2 | 29/2 | 0.06274360 | 0.06274360 | 2.2e-11 | 0.06274354 | 0.06274358 | 9.1e-07 / 2.9e-07 | 207 s, 1.69 GiB |
| 2 | 35/2 | 0.05302721 | 0.05302721 | 1.5e-10 | 0.05302717 | 0.05302720 | 8.3e-07 / 2.5e-07 | 227 s, 1.69 GiB |
| 2 | 41/2 | 0.04591550 | 0.04591550 | 4.8e-09 | 0.04591547 | 0.04591549 | 6.3e-07 / 2.0e-07 | 287 s, 1.69 GiB |
| 4 | 21/2 | 0.09435789 | 0.09435789 | 1.3e-09 | — | — | — | 360 s, 2.57 GiB |
| 4 | 23/2 | 0.08776927 | 0.08776927 | 2.7e-09 | 0.08776901 | 0.08776917 | 3.0e-06 / 1.2e-06 | 358 s, 2.57 GiB |
| 4 | 29/2 | 0.07256402 | 0.07256402 | 1.5e-08 | 0.07256388 | 0.07256397 | 2.0e-06 / 7.7e-07 | 448 s, 2.57 GiB |
| 4 | 35/2 | 0.06184636 | 0.06184635 | 5.5e-08 | 0.06184627 | 0.06184633 | 1.5e-06 / 4.9e-07 | 365 s, 2.57 GiB |
| 4 | 41/2 | 0.05388590 | 0.05388589 | 1.6e-07 | 0.05388585 | 0.05388589 | 8.9e-07 / 1.9e-07 | 335 s, 2.57 GiB |

The pilot's value at σR = 4, m = 41/2 (0.05388590 and 0.05388589) is the Mac test's, to every printed digit. #75's
values fall strictly with m through the five sectors computed, from 0.08301764 to 0.04591550 (σR = 2) and from
0.09435789 to 0.05388590 (σR = 4). #75's A is above #120's A and B in all 16 comparisons, and #120's B is above its A,
so the B-spline values approach #75's from below. #120's A vs B gap is then the difference of their distances from
#75's: at σR = 4, m = 23/2, 3.0×10⁻⁶ − 1.2×10⁻⁶ = 1.8×10⁻⁶. The distances from #75's A are 6.3×10⁻⁷ to 3.0×10⁻⁶ for
#120's A and 1.9×10⁻⁷ to 1.2×10⁻⁶ for its B, and fall with m. They are errors of the B-splines only where #75's own A vs
B (up to 1.6×10⁻⁷) is smaller than them, which it is not for B at σR = 4, m = 41/2 (1.9×10⁻⁷).

### Each criterion

| criterion of #125 | verdict |
|---|---|
| #75's code reproduces its m = 21/2 values (0.083018 and 0.094358) to its printed digits | pass. `repro`: 0.08301764 and 0.09435789 (A), which round to 0.083018 and 0.094358 |
| at each of the eight (σR, m) the leading \|μ\|v_c agrees with #120's to 10⁻⁴ relative, a planted factor 1.001 on one value fails | pass in 8 of 8 against #120's A and B, worst 3.0×10⁻⁶ (σR = 4, m = 23/2, against A), 33 times inside the threshold. `vc` agrees to 10⁻⁶ in all 8. `--mutate=plant` fires: exit 1, `vs120` fails at σR = 4, m = 29/2 (1.0×10⁻³) and nothing else fails |
| a convergence statement per sector: the leading value at two resolutions, and their difference | pass. A and B differ by 3.1×10⁻¹² to 4.8×10⁻⁹ (σR = 2) and 1.3×10⁻⁹ to 1.6×10⁻⁷ (σR = 4), under #75's `converged` threshold of 10⁻⁶ in all ten, 6.3 times inside at worst. `--mutate=conv` fires: exit 1, `conv` fails at σR = 2, m = 41/2 (1.0×10⁻⁵) and nothing else fails |

Of the 28 checks, 2 are `repro`, 10 `conv`, 8 `vc` and 8 `vs120`.

### The speed bound

The lines are unchanged: 12.05 v_c and 10.60 v_c, the reciprocals of #75's m = 21/2 values. #75's leading values at
m = 23/2 to 41/2 are all below them, the largest 0.07681409 (0.9253 of 0.083018) and 0.08776927 (0.9302 of 0.094358),
both at m = 23/2. The dense solver returns every eigenvalue of the Galerkin matrix, so for these ten sectors the leading
value is the largest of all of them, not of the 64 that ARPACK returned in #120. That #120's leading values match #75's
in eight sectors also shows that ARPACK did not miss the leading eigenvalue in them. This covers four of #120's ten
m per σR (23/2, 29/2, 35/2, 41/2) and no others. The range over which the fall is computed is unchanged (m ≤ 41/2), and
for m > 41/2 the decline is still an assumption.

### What the two codes share, and what it does not show

- What they share is the operator S⁻¹M, the tanh profile and v_c's definition. #75's basis (Jacobi polynomials in y =
  r/(r + Lmap) times spin-angular functions, dense eigensolver) shares none of the code of #120's (B-splines, ARPACK).
  Each script computes v_c by its own quadrature of one formula (#75's adaptive `quad`, #120's 400-point Gauss–Legendre
  on [0, 15]), and `vc` compares them to 10⁻⁶. The agreement therefore tests the discretisation, the solver and their
  sizes. It does not test the formulation: an error in S⁻¹M or in the profile would be in both codes, and only the
  planted checks of each (the exact spectrum of the Frank–Loss weight, 7.3×10⁻¹⁵ here) test the operator, at a different
  weight from tanh.
- The margins of the two thresholds narrow with m at #75's fixed sizes. A vs B grows from 1.3×10⁻⁹ (m = 21/2) to
  1.6×10⁻⁷ (m = 41/2) at σR = 4 and from 3.1×10⁻¹² to 4.8×10⁻⁹ at σR = 2, against the threshold 10⁻⁶. Nothing is
  extrapolated here, and sectors beyond m = 41/2 were not run.
- The script prints the L² residual ‖Sψ − vMψ‖/‖Sψ‖ of the leading eigenfunctions, with no threshold. At σR = 4 it is
  1.4×10⁻³ at m = 21/2, the top of the 2×10⁻⁴ to 1.4×10⁻³ that tanh_spectrum.md quotes for #75's run, and rises to
  8.6×10⁻³ at m = 41/2 (σR = 2: 5.0×10⁻⁵ to 6.8×10⁻³). tanh_spectrum.md notes that the eigenvalue converges far faster
  than the eigenfunction. The leading eigenvalues here agree with #120's to 3.0×10⁻⁶ at residuals up to 8.6×10⁻³.
- The script's `real_identity` check, which is not a criterion of #125, is 8.0×10⁻⁹ at σR = 2, m = 41/2, against its
  threshold of 10⁻⁸, and 0 in the other nine sectors. It passes. The quantity is the maximum, over the evaluated
  eigenfunctions, of max(|q_rel|, cluster Gram) × |Re μ|/|μ|.
- Twelve of #120's 20 sectors (m = 25/2, 27/2, 31/2, 33/2, 37/2 and 39/2 at both σR) were not re-solved with #75's code.
- On Windows, `tanh_spectrum_highm.py` as the task states the command fails before any check runs: `printed75` opens
  tanh_spectrum.md without an encoding, and Python 3.14 reads it as cp1252 (`UnicodeDecodeError`, exit 1,
  `I125-win/highm-cp1252-error.log`). All three `highm` runs above are `python -X utf8`. The script is unmodified.
- #94 and #111 found #75 and the B-spline basis to agree to 6.0×10⁻⁶ at m = ½ to 21/2. Here, at #120's sizes (#94's),
  they agree to 3.0×10⁻⁶ at m = 23/2 to 41/2.

## What this does and does not show

Status words of AGENTS.md. These are *numerical evidence*. A Galerkin eigenvalue of a non-normal operator is not a bound
on the spectrum, nor is a window, and the two discretisations agreeing is the evidence of its stability, not a proof of
it. Nothing here is *proved*, *certified* or *checked symbolically*, and "no real zero mode below |v| ≈ 5 v_c" (#94's
sizes) and "below 12.05 v_c and 10.60 v_c" (#111's) are *conditional* on the leading eigenvalue's decline with m. A
window set by the ARPACK floor says where this basis is known to agree, not where it stops agreeing: 14 of #111's 18
re-solved windows were, and after #119 two are, σR = 2, m = 19/2 and 21/2 (#127). That E is empty for the tanh profile
is not shown by this check either, and stays a *conjecture*.

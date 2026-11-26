# Step 3: the ghost sector counted at D = 10, 12 and 20.04, and the other sector's zero on Im μ = ω

Issue [#30](https://github.com/nmoinvaz/negative-mass/issues/30). Run on 2026-11-26 on the coordinator's Mac, read-only on
choptuik-large-d (CLD) as in [`wormhole_next.md`](../wormhole_next.md) §1.3 and §6, at 2 threads per job and two jobs at
once. All 15 jobs exited 0 in 16 min. Status of every result below: **numerical evidence**.

- **Launcher.** `wormhole_next/step3.sh` (queues A and B; `R2` was not needed and not run).
- **Table and checks.** `wormhole_next/step3_table.py`, output in `step3_output.txt` (15 checks PASS, exit 0; the
  self-test makes each of its 8 mutations fail exactly its target). Run it as `cd notes/wormhole_next && python3
  step3_table.py` with the raw clone beside the checkout (`--raw=DIR` otherwise).
- **Raw.** [`negative-mass-raw/I30-mac`](https://github.com/nmoinvaz/negative-mass-raw/tree/main/I30-mac) at b338967: its
  README holds the criteria, edges and planted pairs, written before the first run, and the run log. Counts, caches and
  logs are there.
- **Environment.** CLD at `5fcadda2`, `hp/` clean, kernel `modecount_cx.cpython-313-darwin.so` of 2026-11-08. That is not
  the 2026-10-25 build of §1.3 and §6, and `modecount_common.py` now carries #441's series gate. All three states pass the
  gate. The same environment line is in all 15 job logs. Against step 2's build, the real-axis zeros agree to 1.8×10⁻¹² or
  better (E1) and the two on-line controls to the last digit (C5).

## 1. Winding counts

Each count is CLD's `modecount_test.py` with the boxes B0, U0 and full (aLL < Re μ < aR, |Im μ| < 1.5ω), the left edge aLL
just right of the real-axis zero by step 2's rule (`pick_edges.py` on the stored zero), and a planted pair inside B0. Every
plain box is resolved (every phase step below π/4, the alias rule met, no verification mismatch) with winding of g equal
to 0 (C1). Every planted run counts the pair as +2 in B0 and leaves U0 and full unchanged (C2).

| D | state (M) | counted here | zeros besides μ = 0 | planted pair | determinants (plain, with plant) |
|---|---|---|---|---|---|
| 10 | up_D10 (48) | −0.596 < Re μ < 6.70 | 0 | −0.45 ± 0.6i in B0: +2 | 84, 108 |
| 12 | up_D12 (48) | −0.715 < Re μ < 6.70 | 0 | −0.55 ± 0.7i in B0: +2 | 92, 117 |
| 20.04 | up2_D20.0437 (32) | −0.794 < Re μ < 6.31 | 0 | −0.65 ± 0.8i in B0: +2 | 92, 112 |

The real-axis zeros that fix aLL are §1.3's, and were re-found in this environment before the counts (E1):

| D | μ stored | μ re-found | difference |
|---|---|---|---|
| 10 | −0.676124203429 | −0.676124203430 | 1.4×10⁻¹³ |
| 12 | −0.734531644203 | −0.734531644204 | 5.9×10⁻¹³ |
| 20.04 | −0.874313093543 | −0.874313093541 | 1.8×10⁻¹² |

## 2. The other parity sector, refined on the line

§6.2's scan (`--lo=-1.6 --hi=-0.02 --h=0.05`, `axis_s1.jsonl`) found one phase jump at D = 10 and at D = 12, in
(−1.05, −1.00), below the real-axis zero (C3). Its refinement ended at (−1.2936263734, Im/ω = 0.999583) and
(−1.3247896060, 0.999684), outside the interval and off the line.

**Those two points were not zeros.** Interpolated from the scan's own samples at the same Re, the line's ln|g| differs
from ln|g| at them by 0.002 and 0.004 (`step3_output.txt`). A zero sits tens of units below the line's trend: at each
refined zero below, ln|g| is 30 below the cell's end samples.

**Why the refinement left the line.** `modecount_test_axis.py` refines with `mpmath.findroot(solver="secant", tol=1e-24,
maxsteps=40, verify=False)` from the pair (cell midpoint, cell midpoint + 0.25h + 0.001i). The iteration is complex, so
nothing holds it on the line, and with `verify=False` its 40th iterate is returned unchecked. At h = 0.05 the start is
0.023 (D = 10) and 0.022 (D = 12) from the zero found below, farther than every §6.2 case that converged (up to
0.0199). The margin is thin, so the evidence is mostly the contrast of the two rungs below.

| run | grid | D = 10 | D = 12 |
|---|---|---|---|
| rung 0 (§6.2's start, cell (−1.05, −1.00)) | h = 0.05 | Re −1.2928, Im/ω − 1 = −4.1×10⁻⁴ | Re −1.3250, −3.2×10⁻⁴ |
| rung 1 (the record) | h = 0.005 | Re **−1.0016613614**, Im/ω − 1 = 0 | Re **−1.0033895309**, Im/ω − 1 = 0 |

Rung 0 fails the on-line and in-cell tests at both D (X1: they fire), so the old failure repeats. Rung 1 passes C4: one
jump, in the cell (−1.005, −1.000); |Im/ω − 1| < 10⁻⁵ (here below 10⁻¹⁵); the zero's Re inside the cell; ln|g| 30.3 and
29.6 below the cell's lower end sample. **Yes, the refinement stays on the line**, and the controls at D = 16.76 and
20.04 return §6.2's zeros to the last digit (−1.0033192792, −1.0026635066), also on the line (C5).

Two limits on the explanation:
- Start distance and first step both scale with h (0.25h), so the runs do not separate them. What they show is that the
  start pair of §6.2 fails and a start within 0.0025 works.
- Rung 0 reproduces §6.2's points only to 8×10⁻⁴ and 2×10⁻⁴ in Re, not bit for bit. The script multiplies out only the
  poles within 2 of the window, so a narrower window changes g and with it the secant's path.

The zeros stay on the line to 10⁻¹⁶ in every converged case. det M(conj μ) = conj det M(μ), with the folded period 2ω,
makes the zero set symmetric about the line, so a zero is on it or part of a pair straddling it. The phase jumps by π
along the line (−2.0797 − 1.0633 = −3.143 at D = 10, raw `axis_s1.jsonl`), not 2π, which says one zero, on the line. This is **argued** here,
not shown.

## 3. What these results do not show

- **One M and layout per D.** The counts and zeros use M = 48 (D = 10, 12) and M = 32 (D = 20.04). The resolution error of
  μ is ≤ 2×10⁻¹² at D = 5 and 6 (CLD #632), and nothing here repeats that at the new D.
- **The planted pair is in B0 only.** U0 and full have no planted control at these D (D = 4's full count has one). Their
  winding 0 rests on the resolved flag and on CLD's own mutation controls for the sampler.
- **Shared machinery.** The counts, the line scans and the real-axis zeros all use CLD's `test_det`, the states and the
  sampler. The planted pair tests the winding code and the box geometry, not the determinant. No independent determinant
  exists (CLD #646 is open).
- **Between the counted D.** Only 14 values of D from 4 to 67 are counted. CLD holds states at a spacing of about 1 in D
  from 8.25 to 60 (`hp/data/states`), most of which no count has used. Each would have to pass the series gate.

## 4. What wormhole_next.md should now say

For the coordinator; `wormhole_next.md` is not edited here.

**§6.1.** Add three rows to the table:

| D | state (M) | counted here | with CLD's | zeros besides μ = 0 | planted pair | evaluations |
|---|---|---|---|---|---|---|
| 10 | up_D10 (48) | −0.596 < Re μ < 6.70 | | 0 | −0.45 ± 0.6i: +2 | 84 |
| 12 | up_D12 (48) | −0.715 < Re μ < 6.70 | | 0 | −0.55 ± 0.7i: +2 | 92 |
| 20.04 | up2_D20.0437 (32) | −0.794 < Re μ < 6.31 | | 0 | −0.65 ± 0.8i: +2 | 92 |

and read "At these 14 values of D from 4 to 67" for "11" (also in the summary and in §6.5). The values of D between those
counted are still not counted, but D = 10, 12 and 20.04, the three with real-axis rates and no count, now are.

**§6.2.** Replace the row "10, 12 | (−1.05, −1.00) | (the refinement left the line)" with

| D | phase jump in | refined on the line |
|---|---|---|
| 10 | (−1.05, −1.00) | −1.001661 |
| 12 | (−1.05, −1.00) | −1.003390 |

and replace "From D = 10 on it lies just below −1 and approaches it, the departure falling faster than 1/n² at large D" by:

> It crosses −1 between D = 8.25 (−0.996904) and D = 10 (−1.001661). The departure below −1 is 0.00166 at D = 10, 0.00339
> at D = 12, 0.00332 at D = 16.76 and 0.00266 at D = 20.04, then falls: its maximum lies between D = 10 and 17 (D = 12 and
> 16.76 differ by 7×10⁻⁵) and is not located. The refinements at D = 10 and 12 use a 0.005 grid, because the secant of the 0.05 grid left the line
> (step3.md §2).

The claim that the departure falls faster than 1/n² at large D is unchanged and was not tested here.

**§6.5.** "The other parity sector's slowest mode is more damped at every D scanned" stands, now with a refined zero at
every D scanned from 5 to 67 (C3 holds at D = 10 and 12). Under "Remaining", "Values of D between the counted ones" stays,
with the states in §3 above. Add: the real-axis zeros at D = 10, 12 and 20.04 and the on-line zeros at D = 16.76 and 20.04
were re-found on a newer kernel build (2026-11-08 against 2026-10-25) and agree with the stored ones to 2×10⁻¹² or better.

## Coordinator's review

Reviewed 2026-11-26 at ad00567 (raw `I30-mac` b338967); the review's files are in `negative-mass-raw/I30-mac-review`
(raw 4a0dca0).

| item | result |
|---|---|
| scope | the four files the task named, and the raw folder; nothing written in choptuik-large-d |
| rerun | `step3_table.py` on the raw records: exit 0, every line equal to `step3_output.txt` |
| re-solve | CLD's `modecount_test_axis.py` at D = 10 with the cell h = 0.0025 (the record used 0.005): the zero −1.0016613614 + 2.4922604598i with Im/ω = 1.000000, the record's value to every printed digit |
| independent recount `ind30.py` | **fails as written, by the coordinator's design.** The records store arg det M, whose phase moves by up to 3.1 rad between samples. The count winds g = det M · exp(lnqk − p_c), so a winding of det M cannot be recovered from these samples, and an independent recount would need CLD's own normalization. The planted pairs (C2, +2 in every B0) remain the counts' control |

Findings: the answer stands as numerical evidence. The text for `wormhole_next.md` §6.1, §6.2 and §6.5, and the count in
PLAN.md, are applied. §6.2's old sentence ("the departure falling faster than 1/n²" from D = 10 on) was wrong about
D = 10–17, as this note shows; the large-D claim is untested.

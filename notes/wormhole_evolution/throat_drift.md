# The throat radius drift of family G at D = 5, its dependence on the normal pulse, and the expanding member at X = 160 (#255, #42)

Family G of #40 (a ghost pulse A tuned to p*ₕ plus a normal pulse B on the same ray) at D = 5, with the normal amplitude B = 0.04, 0.1 and 0.2 (and, post hoc, 0.02 and 0.07), and the expanding member A = 0.05 on a
rectangle of side X = 160. #39's evolution (`dn.row_step`, unchanged) is run row by row, the throat locus of #225 is tracked, nothing of the grid is stored. Driver: [throat_drift.py](throat_drift.py) (imports `dn.py`, `checks.py`,
`threshold.py`, `expansion.py` and `throat_locus.py` unchanged). Output of record: [throat_drift_output.txt](throat_drift_output.txt)
(`python3 throat_drift.py analyze --raw=<negative-mass-raw>/I255-dh610/runs`, exit 1: 11 of 48 checks fail, listed below). Raw output, with the prediction and the criteria written before the runs: negative-mass-raw `I255-dh610`.
Status words as in AGENTS.md: every number below is *numerical evidence*; no sentence says a result is new.

## Answer

**The cause of the drift is not established. The prediction that the drift is second order in B is half right: the threshold is (PS met), the throat radius and its speed are not (PR and PV not met).** At B = 0.04 and 0.1 the tuned amplitude
scales as p*ₕ ∝ B^2.04 and B^2.11 (h = 0.00354, 0.00177), so the normal pulse does enter the threshold at second order. The throat radius does not: r_th − 1 at t = 8, 10, 12 scales as B^2.6 to B^2.9 between 0.04 and 0.1, and the speed
ṙ∞ that dr_th/dτ approaches scales as B^3.16 (declared band 1.8 to 2.4 for both). The extra scan (post hoc) puts the local exponent of r_th − 1 at 2.2 to 2.4 (B = 0.02 to 0.04), 2.5 to 2.8 (0.04 to 0.07), 2.8 to 3.0 (0.07 to 0.1) and 2.2
to 2.4 (0.1 to 0.2), and ṙ∞ ≈ 0.110 E_net^(3/2) to 0.3 % at B = 0.02, 0.07 and 0.1 (0.093, −15 %, at B = 0.04), with E_net = (B² − A*²) I₂ the energy the two pulses leave. E_net^(3/2) ∝ |B|³: an odd power of B, which an
expansion in B² excludes, so the speed of the radius is not a perturbative response to the pulse (not established beyond that). What the runs do show:

- The tuned throat is quasi-static: on the locus the ghost's energy density is that of a static throat of the same radius, |ρ_g/ρ_eq − 1| ≤ 3e-4 (B = 0.04) and 3.6e-3 (B = 0.1) over t = 8 to 12 (QS met), 1.0e-2 at B = 0.2 (not met).
- It is not driven by the normal pulse (ρ_n/|ρ_g| = 2e-7 at B = 0.04 and 1e-5 at B = 0.1 at t = 10) and not by the unstable mode (the pair's odd part is at most 7e-5 at B = 0.1 and 3.4e-4 at B = 0.04 for t ≤ 12).
- It converges at second order in h: r_th(10) − 1 at the three h has observed order 2.02 for B = 0.04 and 0.1 (the two finest differ by 0.78 % and 0.095 %).
- The fitted static member with the data's ghost charge holds only where there is no drift to explain: Q/q₀ = 0.997 to 1.001 at B = 0.04 (t = 6 to 12) and 1.02 to 1.23 at B = 0.1.
- At D = 4 the radius relaxes (B = 0.04: r_th − 1 = +0.0009 at t = 8, −0.0009 at t = 12, ṙ∞ = 0.0000 ± 4e-6; B = 0.1: ṙ∞ = +0.0013); at D = 5, B = 0.1 it is +0.026, 20 times larger.

**The expanding member is class FT** (#213's classes, the locus's proper time): D = 5, G, A = 0.05, X = 160, k = 3.9957 at both h = 0.00707 and 0.00354, τ_s/τ_end = 1.1759 and 1.1760 (bound 1.2; #225's X = 100 gave 1.2025 and
class GROW). KX2 and RES2 are met.

| acceptance (issue #255) | verdict | numbers |
|---|---|---|
| r_th(t) explained against the data's energies, or its B dependence measured at two h, with the second-order convergence in h stated | B dependence measured at two h: **met**; explained against the energies: **not established** | B = 0.04, 0.1 at three h, B = 0.2 at two (the third was stopped at 2 h); orders 2.02, 2.02; PS met, PR and PV not met |
| the expanding member's class and k at two h | **met** | FT, k = 3.9957, 3.9957; τ_s 9.1290, 9.1226 |
| a planted case through the same locus code recovers its known value | **met** | static throat of radius 1.2: 4.6e-5 (bar 2e-4), second order; moving throat v = 0.13, −0.3: r_th = 1 to 1.2e-4 |

## What was run

| item | declaration (README of the raw folder, written before the runs) |
|---|---|
| family | G of #40, normal amplitude B = 0.04, 0.1, 0.2 (0.02, 0.07 added post hoc by the coordinator's ruling); Gaussian pulses in φ, centre 3, width 0.5, on the ray x⁻ = 0 |
| D = 5, h | 0.00707, 0.00354, 0.00177 (exact `threshold.CRIT['H'][3]`); gates use the two finest |
| p*ₕ | B = 0.1: #40's. Others: bisection of #40's events to 1e-12 (37 to 41 runs, 90 s to 57 min); bracket (0.3, 3) × (B/0.1)² p*ₕ(0.1, h) |
| the pair | A = p*ₕ ∓ 1e-10 kept whole; the tuned series is its mean; odd = (r₊ − r₋)/2 is the unstable mode's part |
| fit | dr_th/dτ = ṙ∞ + c e^(−κ(τ−τ₀)) over t in [6.5, 12] with |odd| ≤ 3e-4 |
| expanding member | D = 5, G (B = 0.1), A = 0.05, X = 160 (t to 226), h = 0.00707, 0.00354; 4 and 31 min; streamed (two grid rows and the locus's per-slice arrays: 306 MB at most over all jobs) |
| cost | 4 lanes, nice 10, load 12 to 23 of 20; the B = 0.2, h = 0.00177 bisection was stopped by its 2 h timeout after 33 runs (near threshold the + side event is at t = 38 to 56) and not resumed (coordinator's ruling) |

## The prediction and its verdicts

Written before any run with B ≠ 0.1 (raw README, Stage 1): the normal field is a test field at first order, so the tuned amplitude is O(B²) and every effect of the tuned data on the throat is O(B²) with the same time dependence.

| check | claim | result |
|---|---|---|
| PS | exponent of p*ₕ between B = 0.04 and 0.1 in [1.8, 2.4] | **met**: 2.044 (h = 0.00354), 2.107 (0.00177); 1.83 at h = 0.00707 (not gated) |
| PR | exponent of r_th − 1 at t = 8, 10, 12 in [1.8, 2.4] | **not met**: 2.610, 2.769, 2.857 (h = 0.00354); 2.602, 2.761, 2.851 (0.00177) |
| PV | exponent of ṙ∞ in [1.8, 2.4]; κ(0.04)/κ(0.1) in [0.75, 1.25] | **not met**: 3.163, 3.162; κ ratios 1.026, 1.033 are inside |
| CV | r_th(10) − 1 second order over three h | met for B = 0.04, 0.1 (order 2.02); B = 0.2 not evaluable (two h) |

The numbers behind PR and PV (h = 0.00177; r_th − 1 = R):

| B | p*ₕ | R(8) | R(10) | R(12) | ṙ∞ | κ |
|---|---|---|---|---|---|---|
| 0.04 | 1.554974e-3 | 0.00812 | 0.01031 | 0.01292 | +0.00144 ± 6e-6 | 0.497 |
| 0.1 | 1.071568e-2 | 0.08805 | 0.12944 | 0.17601 | +0.02602 ± 4e-5 | 0.481 |
| 0.2 (h = 0.00354) | 6.285376e-2 | 0.47036 | 0.64361 | 0.81343 | the fit does not converge | |

Predicted from B = 0.1 by B²: R(10) = 0.0207 and ṙ∞ = 0.0042 at B = 0.04; measured 0.0103 and 0.00144.

## The supplementary scan (post hoc, gates nothing)

D = 5, h = 0.00354 (h = 0.00707 gives exponents within 0.2 of these for R and 0.1 for ṙ∞; both are in the output):

| B | p*ₕ | R(8), R(10), R(12) | ṙ∞ | E_net | ṙ∞/E_net^(3/2) | d ln R/dτ at t = 10 | linear, exponential rms of R(τ) over t in [8, 12] |
|---|---|---|---|---|---|---|---|
| 0.02 | 5.179e-4 | 0.00172, 0.00206, 0.00245 | 0.00021 | 0.01548 | 0.1097 | 0.092 | 3.9e-3, 1.0e-3 |
| 0.04 | 1.6625e-3 | 0.00805, 0.01023, 0.01283 | 0.00143 | 0.06186 | 0.0932 | 0.119 | 5.9e-3, 2.3e-3 |
| 0.07 | 5.006e-3 | 0.03245, 0.04545, 0.06090 | 0.00898 | 0.18881 | 0.1095 | 0.162 | 7.8e-3, 7.2e-3 |
| 0.1 | 1.0819e-2 | 0.08794, 0.12932, 0.17587 | 0.02601 | 0.38277 | 0.1098 | 0.181 | 5.6e-3, 1.2e-2 |
| 0.2 | 6.285e-2 | 0.47036, 0.64361, 0.81343 | no fit | | | 0.163 | 2.4e-3, 1.4e-2 |

Local exponents (0.02 to 0.04, 0.04 to 0.07, 0.07 to 0.1, 0.1 to 0.2): p*ₕ 1.68, 1.97, 2.16, 2.54; R(10) 2.31, 2.67, 2.93, 2.32; ṙ∞ 2.76, 3.28, 2.98. At the small B the p*ₕ exponent is below 2, presumably because p*ₕ carries an
h² offset that weighs more at small p* (the coarse h gives 1.14 and 1.70); not separated here. The growth of r_th − 1 is closer to exponential at B = 0.02 and 0.04 and closer to linear (a constant speed) at B = 0.1 and 0.2. The transient rate κ of the fit is 0.65, 0.50, 0.38, 0.48 at
B = 0.02, 0.04, 0.07, 0.1: its B independence at 0.04 and 0.1 (PV's κ part) does not extend to the scan, so it is not established.

## What is ruled out, what is not

| candidate | verdict | evidence |
|---|---|---|
| flux of the normal pulse through the locus | ruled out after t ≈ 8 | ρ_n/\|ρ_g\| at t = 10 is 5e-8, 2e-7, 1e-6, 1e-5 at B = 0.02, 0.04, 0.07, 0.1 |
| the unstable mode | ruled out | odd part ≤ 6.8e-5 (B = 0.1), 3.5e-4 (B = 0.04, NO not met) for t ≤ 12, second order in the mean |
| resolution | ruled out for B = 0.04, 0.1 | order 2.02 over three h |
| the observer (the minimum of r on a slice) | ruled out | planted moving throat v = 0.13, −0.3: r_th = 1 to 1.2e-4, rate ≤ 2.2e-4 (PL2); static throat of radius 1.2 to 4.6e-5, second order (PL1) |
| a first-order (in B) response | ruled out | p*ₕ ∝ B^2.04 to B^2.11 (PS) |
| a second-order (in B) response of the radius | ruled out | PR, PV; the exponents rise from 2.3 to 3 with B |
| a static member with the data's ghost charge | holds at B = 0.04, fails at 0.1 | Q/q₀ = 0.997, 0.995, 0.997, 1.001 at t = 6, 8, 10, 12 (B = 0.04); 1.023, 1.061, 1.125, 1.230 (B = 0.1); 1.45 to 2.67 (B = 0.2) |
| the constraint plateau (#226) | not tested here | roundoff of 1e-11 against ṙ∞ = 0.026 |
| a generalized zero mode (a Jordan block) of the static family; a slow tail of the ghost field; a drift of the lapse relative to the far region | Jordan block ruled out; a slow tail is present in the linear problem (#272) | no bounded secular partner exists (the zero modes' partners t·Y are exact and Y^(1) grows as l³); the resolvent has a double pole with a logarithm at W = 0 at D = 5 (throat speed ∼ t^−0.2) and a simple pole at D = 4 (shift saturates, speed ∼ t^−3); no pole near κ ([ellis_spectrum.md](ellis_spectrum.md)); at D = 4, B = 0.02 the nonlinear speed is a power law with no drift dominating before t = 15, and its local rate falls from 0.21 at t = 10 to 0.14 at t = 17 ([throat_tail.md](throat_tail.md), #291) |

The tuned throat radius follows a quasi-static sequence of throats whose growth speed approaches 0.026 at B = 0.1 (D = 5); what fixes that speed (∝ |B|³ in the scan, ∝ E_net^(3/2)) is the open question. At D = 4 the same measurement gives
ṙ∞ ≤ 0.0013 and a relaxing radius (B = 0.04: r_th − 1 = +0.0009, −0.0002, −0.0009 at t = 8, 10, 12; B = 0.1: +0.0095, +0.0050, +0.0032; B = 0.2: the radius runs away, 2.1 at t = 12).

## The expanding member

| member | h | class | k | τ_s | τ_end | τ_s/τ_end | final window r/a | H/W at r/a = 10, 25, 50 |
|---|---|---|---|---|---|---|---|---|
| G, A = 0.05, X = 160 | 0.00707 | FT | 3.9957 | 9.1290 | 7.7634 | 1.1759 | 50.5 to 137.2 | 0.871, 1.105, 1.315 |
| G, A = 0.05, X = 160 | 0.00354 | FT | 3.9957 | 9.1226 | 7.7572 | 1.1760 | 50.5 to 137.2 | 0.871, 1.105, 1.316 |
| G, A = 0.05, X = 100 (#225) | 0.00707, 0.00354 | GROW | 3.9883 | 9.1260, 9.1196 | 7.5896, 7.5834 | 1.2024, 1.2026 | 31.3 to 85.0 | the same |

1/H is linear in τ to 7.3e-5 (X = 100: 2e-4); H max/min over the window 1.284, |ρ| max/min 1.649, w̄ = −0.5000 exactly; the two h agree on H at fixed r/a to 5e-5 and on k to all digits (RES2). The locus ends at l = +3.23 off the
diagonal (v = +0.0144). #211's slope test is not met (H/W = 0.871, 1.105, 1.315 against bar 10 %), as at X = 100 and for S. k = 3.9957 against #213's 3.9897 (S, A = 0.05, a window ending at r/a = 85) and 4 for (n, k) = (3, 4)
(`expansion.md`'s |ρ|/H² = 27/4): k moves toward 4 as the window moves out (3.9883 at r/a = 31 to 85, 3.9957 at 50 to 137). The X = 160 runs reproduce #225's X = 100 locus series bit for bit on every slice inside both rectangles (CT3).

## Checks and mutations

48 checks, 37 pass. The 11 that fail are results, not repaired: PR and PV (above); NO at B = 0.04 at all three h (max |odd| 3.4e-4 to 3.5e-4 against 3e-4, at t = 12; the fit masks the points above 3e-4, which removes the last 0.1 in t; the mean
differs from the tuned member by the odd part squared, 1e-7); QS at B = 0.2 at two h (1.0e-2, 9.9e-3 against 5e-3); JOB at B = 0.2 at all three h (the + side run does not reach #40's expansion event inside the triangle t ≤ 34 at h = 0.00707,
0.00354, and the locus's r_min = 2 level is crossed on both sides because the throat radius exceeds 2 during the stall, so the locus fates are (+1, +1); at h = 0.00177 the bisection was stopped at 2 h); CV at B = 0.2 (three h needed).
PL1, PL2, PL3, PL5, PX (exponents 2 pass; 1 and 3 fail all three gates), CT1 to CT4 (bit for bit: #40's 8 points, #225's pair at B = 0.1, #225's X = 100 locus, #40's p*ₕ at D = 4 to 1.4e-13), DL (1e-9 against 1e-10: ≤ 8.4e-8), RV (ṙ∞ = 0.02602 at
B = 0.1), KX2, RES2 pass.

13 mutations, each target computed in development before it was declared (raw README Stage 2, `dev/mut`). On the stored record every target fires except where the check already fails (PR and PV, so `b_label`, `r_not_r1` and
`rdot_swap` change PS, PX1, PX2 only) and `ts_bound` (the relaxed bound leaves KX2 passing); `runs/mut` has the logs.

## Discrepancies and things that failed

1. PR and PV not met (above): the brief's "second order" for the radius is wrong beyond the threshold.
2. B = 0.2 at h = 0.00177 not finished (2 h stop). Near threshold at B = 0.2 each run reaches t = 38 to 56 (the expansion side) and takes 10 to 30 min; 14 bisection steps were left. B = 0.2 is reported at two h.
3. The triangle cap of the pair runs (t ≤ 34, inherited from #225) is too short for the + side event at B = 0.2; the series in t ≤ 12 is not affected.
4. NO fails at B = 0.04 by 13 % (the bar was set from B = 0.1's odd part; the unstable mode grows faster at the smaller radius).
5. During development: the first control bisection stopped at a bracket of 4.7e-12 (triangle cap too short at D = 4); the cap was made to grow with the bracket (`dev_09`, `dev_10`). The first `analyze` over the finished jobs stopped on a report-only bug in the mass table (`runs/analysis_partial1.log`).
6. Report-only code after the declaration (6f060f4): 6603ac1 (the stopped job as not met, the post-hoc tables, the mass-table fix) and a44d066 (the energy table); no check, bar, window or mutation changed.

## What `throat_locus.md`, `critical_ghost.md` §1.3 and #42 should now say (proposals for the coordinator)

`throat_locus.md`, the table row "the expanding side on the locus" and the paragraph after it: D = 5, A = 0.05 at X = 160 is class FT, k = 3.9957, τ_s/τ_end = 1.176 at both h (#255): KX is met once the locus reaches τ_s/τ_end ≤ 1.2. Recommended task 5 is done,
task 3 is partly done: the drift's B dependence is measured (`throat_drift.md`), its cause is not established.

`critical_ghost.md` §1.3, the bullet on the generic one-sided family: at D = 5 the tuned state's throat radius grows at a speed that approaches ṙ∞ = 0.026 a per unit proper time at B = 0.1 (D = 4: ≤ 0.0013); the speed scales as |B|³ over B = 0.02 to 0.1 (post hoc),
the threshold p*ₕ as B²; the throat stays quasi-static (ρ_g within 4e-3 of a static throat's); the cause is not established (*numerical evidence*).

**#42:** the member G approaches at D = 5 is not a single member: the throat radius is 1.05 (t = 6) to 1.20 (t = 13) at B = 0.1 and 1.007 to 1.014 at B = 0.04 (the fitted member keeps the data's ghost charge to 0.5 % at B = 0.04 and exceeds it by 23 % at B = 0.1 by t = 12). A one-sided
generic family approaches the static wormhole's neighbourhood with a throat that drifts at a speed ∝ |B|³.

## Limits

- B = 0.2 has two h; its speed fit does not converge (the transient is not an exponential there), and its locus fates are ill-defined (r_th > 2 during the stall).
- ṙ∞ is an extrapolated constant of a fit over t in [6.5, 12] (about 9 proper time units); the runs end by t ≈ 18 to 22 (the pair at 1e-10), so whether the speed stays constant is not measured. The B = 0.04 point of the E_net^(3/2) law is 15 % low and unexplained.
- Only p*ₕ to 1e-12 (the roundoff floor of #40) and δ = 1e-10 are used; the pair's mean removes the unstable mode to second order in δ.
- The member fit (s_th, Q/q₀) rests on M(r) over r = 1.25 to 4 r_th, where its rms is 1e-5 to 3e-2 at D = 5 (#225's limit); Q/q₀ is a diagnostic, not a check.
- D = 4 values are reported from two h; D = 6 was not run (#228).

## Recommended tasks

1. *Done in #272* ([ellis_spectrum.md](ellis_spectrum.md)): no pole of the linear resolvent near κ at D = 4 or 5; the transient's "rate" is the local exponent of a power-law tail. **The linear spectrum of the symmetric Ellis wormhole at D = 5 beyond the unstable mode.** What: the even and odd spherically symmetric perturbations of the static wormhole (zero modes, decaying modes, any slowly growing mode), the generalized
   kernel of the scale mode. Why: ṙ∞ approaches a constant and the transient has rate about 0.5 here and 0.24 at D = 4; a Jordan block of the zero modes predicts a constant speed. Acceptance (can fail): the linear mode rates against κ = 0.495 ± 0.01 and 0.244; the scale mode's
   secular term against ṙ∞/E_net at B = 0.02. Depends on: this note, #195.
2. **Why ṙ∞ ∝ |B|³.** What: B = 0.01, 0.014 and 0.03 at h = 0.00354 (about 25 min each) and the same B at D = 4, with p*ₕ by bisection to 1e-12; fit ṙ∞ with the window extended to t = 16 (pairs at 1e-12). Why: an odd power of B cannot come from an expansion in B². Acceptance:
   the exponent between B = 0.01 and 0.04 with its fit error; the B = 0.04 point explained or reproduced at a third h. Depends on: this note.
3. **B = 0.2 at h = 0.00177** with a resumable bisection: needs about 2.5 h more and a waiver of the 2 h stop; not needed for the acceptance. Acceptance: the three-h order of r_th(10) − 1 at B = 0.2. Depends on: this note.
4. **G at D = 6 on the locus** (#228's third h): the drift speed and its B dependence at D = 6. Depends on: #228, this note.

## Files

`throat_drift.py`, `throat_drift_output.txt`, this note; raw `negative-mass-raw/I255-dh610` (README: prediction, criteria, rulings, development runs and results; `dev/`, `runs/` with the bisections, pairs, expansions, controls, extras and the mutation runs; `scripts/`).

## Coordinator's review

Review of branch `wip/i255-g-drift` at 572cde0 (declared code 6f060f4) and raw `I255-dh610` at a1cb65b2, recorded in negative-mass-raw
`I255-mac-review`. On the Mac, nice 10.

- **Scope.** The branch adds the script, its output and this note, nothing else.
- **Rerun.** `analyze` on the stored runs gives the same 48 verdicts (37 pass, 11 fail, as listed). Every verdict line equals the output of
  record.
- **Rulings during the task.** There was no waiver for B = 0.2 at h = 0.00177, so B = 0.2 is reported at two h. The supplementary B scan
  (0.02, 0.07) was allowed as post hoc and gates nothing. Both are recorded on the issue.
- **No independent code.** The coordinator's own march is the symmetric half domain and cannot run the one-sided family G. The result rests
  on the locus code of #225 (its planted moving throats) and the second-order convergence in h.
- **What holds.**
  - The threshold scales as B² (PS).
  - The throat radius and its speed do not: PR and PV are not met, and r_th − 1 grows as B^2.6 to B^2.9.
  - The post-hoc fit ṙ∞ ≈ 0.110 E_net^(3/2) ∝ |B|³ is an odd power, which no expansion in B² gives. Its cause is not established, and
    the B = 0.04 point is 15 % off it.
  - The expanding member at X = 160 is class FT with k = 3.9957, so #225's KX is met once the locus reaches τ_s/τ_end ≤ 1.2.
- **Applied.** `throat_locus.md`, `critical_ghost.md` §1.3, #42 and Paper V's #255 placeholder carry the update.
- **Follow-ups.**
  - Task 1 (the linear spectrum beyond the unstable mode, the scale mode's secular term) is #272.
  - Task 2 (why |B|³) is #273.
  - Task 3 is not filed: it is not needed for the acceptance.
  - Task 4 is added to #228.

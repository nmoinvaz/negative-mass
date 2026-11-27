# The one-sided threshold family's lifetime along the throat locus, and the static member at the stall (#225, #42)

Family G of #40 (a ghost pulse A tuned to p*ₕ plus a fixed normal pulse B = 0.1 on the same ray) at D = 4 and 5, with the lifetime read in the proper time of the minimum-r locus instead of on
the diagonal. #39's evolution (`dn.row_step`, unchanged) is run row by row on a triangle of the grid, with the locus of every slice t = const tracked on the way. Driver:
[throat_locus.py](throat_locus.py). Output of record: [throat_locus_output.txt](throat_locus_output.txt) (`python3 throat_locus.py analyze --raw=<negative-mass-raw>/I225-dh610/runs`,
exit 1: four checks fail, listed below). Raw output, with the criteria written before the runs of record: negative-mass-raw `I225-dh610`. Status words as in AGENTS.md: every slope,
member and class below is *numerical evidence*; no sentence says a result is new, and none says what the end state of an expansion is beyond the declared classes.

## Answer

**At D = 4 the +6 % of #40 is the observer, not the critical solution.** In the locus's proper time G's lifetime slope is 0.84954, 0.84954, 0.84949 at h = 0.02, 0.01, 0.005 (+4.1e-3 above
1/W = 0.846045; the two finest h agree to 5.5e-5), while #40's diagonal slope is reproduced bit for bit (0.898291, 0.898043, 0.899089). The near-critical state is a static ghost wormhole of
the two-parameter family with throat radius r_th = 1.0028, asymmetry γ₁ = 0.048 (end masses m₊ = +0.045, m₋ = −0.052) and no normal charge (|q_n²/q_g²| ≤ 1e-5); its e-folding time
T(γ₁) r_th = 0.8472 (#151's T) agrees with the growth rate of the unstable mode measured on the locus, 1/λ = 0.8478 at t = 12 (0.04 %), and with the series slope (0.2 %). **At D = 5 there is no
plateau**: the locus's throat radius drifts from 1.07 to 1.20 between t = 7 and 13, the growth rate of the mode follows the symmetric member's W/r_th(t) to 0.5 % along it, and the slope of the
series (0.7022, +21 % over 1/W; #40's diagonal +38 %) is an average of r_th(τ)/W over the departures, so it does not equal the e-fold of the member at the stall time (CMP2 fails, +14.6 %).
On the expanding side the locus gives class FT with k = 1.9993 to 1.9999 at D = 4 (S: 2.0000) for three members and two h, and k = 3.9834 (D = 5, tuned; class FT) and 3.9883 (A = 0.05; class
GROW only because τ_s/τ_end = 1.2025 exceeds the declared 1.2); a drift of 0.014 to 0.15 does not change k.

| acceptance (issue #225, #213's addition) | verdict | numbers |
|---|---|---|
| slope of T in the locus's τ agrees between two h within 0.5 % | met | D = 4: 5.5e-5 (h = 0.01, 0.005); D = 5: 6.7e-5 (h = 0.00354, 0.00177) |
| planted moving throat recovers its velocity and the unboosted e-fold | met | static throat, drift v = 0.1, 0.5, −0.3 and another gauge: l, r, f, τ, v to ≤ 3.5e-4, second order; family S with v = 0, 0.1, 0.5 (D = 4), 0.13 (D = 5): locus τ slope within 1.9e-5, 8.4e-6, 3.1e-4, 4.4e-4 of 1/W; diagonal +0.07 %, +7.8 %, +20 %, +11.5 % |
| control: #40's diagonal slope reproduced | met | all 40 points of each of the six jobs bit for bit; 0.898291, 0.898043, 0.899089 (D = 4), 0.798948, 0.800044, 0.800277 (D = 5) |
| the member at the stall: r_th, asymmetry, normal charge | measured | D = 4: r_th 1.0028, γ₁ 0.048, m₊ +0.045, m₋ −0.052, q_n²/q_g² ≤ 1e-5; D = 5: r_th 1.06 to 1.21 (drifting), s_th 0.12 to 0.17 |
| its e-fold compared with the locus slope, reported whether or not they agree | D = 4 agree (0.2 %); **D = 5 do not agree at the stall time (+14.6 %), they agree along the drift (0.5 %)** | CMP1–3, DR |
| the expanding side on the locus: classes, H/W at r/a = 10, 25, 50, k | measured | D = 4 FT, k 1.9993 to 1.9999; D = 5 FT (tuned) and GROW (A = 0.05), k 3.9834, 3.9883; H/W below |

## An error of the brief

The brief reads #40's stall table (`threshold.md`) as "dl/dt = +0.58 at D = 4, +0.25 at D = 5". The table's columns are the locus's position l = (x⁺ − x⁻)/√2 and then its drift dl/dt: +0.584 and
+0.0735 (D = 4), +0.254 and +0.127 (D = 5); the raw `I40-dh610/runs/D4_G_h0.01/stall.json` has `l` = 0.5842, `v` = 0.0735. The +0.58 and +0.25 are positions. This was found in the
development runs and reported to the coordinator, who confirmed it and kept the planted drifts chosen here (v = 0.1, 0.5, 0.13). The measured drift on the locus at the events is 0.060 to 0.095
(D = 4) and 0.072 to 0.141 (D = 5).

## What was measured

- **The locus.** On each slice t = (x⁺ + x⁻)/√2 = const (the anti-diagonals of the grid nodes) the minimum of r with `dn.throat_min`'s parabolic refinement; position l(t), drift v = dl/dt (a linear fit of l over ±0.5 in t),
  lapse f at the minimum, **proper time** dτ = e^{−f/2} √(1 − v²) dt. On the locus r₊ = r₋ (r_l = 0), so H = d ln r_min/dτ needs no field data. The streamed sweep reproduces `dn.throat_min`,
  `dn.evolve` and #40's events to roundoff (r_min 2e-16, l and events identical; bit-identical to the rectangle inside the triangle).
- **Events.** #40's own (first row with a trapped cell; first anti-trapped diagonal cell with r ≥ 2) run in the same sweep and give the fate and the control. The locus's events are the levels
  r_min = 0.5 (collapse) and 2 (expansion), by interpolation; alternatives 0.75 and 3. The lifetime of a point is τ at its own side's primary level.
- **The series.** p*ₕ ∓ 10^{−e}, e = 3, 3.5, …, 12, with p*ₕ from #40's raw (the same floats, so the same runs are repeated), `threshold.fit_set` unchanged (window [1e-11, 1e-3], 15 fits, spread = s ± 2 SE).
  Triangle t ≤ 34 (bracket ends 44). Resolutions: #40's three h at each D.
- **The member at the stall.** Two runs p*ₕ ∓ 1e-9 kept whole, with the profile of every slice t = 2, 3, …. On a slice: the Misner–Sharp mass M = rⁿ⁻¹(1 + 2e^f r₊r₋)/2 on the two sides of the locus at
  r = 1.25 … 4 r_th; the static ghost wormholes of the two-parameter family (r_th, s_th = r_th N′/N at the throat; `static_member` solves rr″ = −(n−1)(r′² − 1) − rr′N′/N, N′ = K/rⁿ) give M₊(r) and M₋(r), and s_th is the least-squares fit
  of M₊ − M₋. At D = 4, s_th = γ₁/√(1 + γ₁²) (GGS09a's family; the solver reproduces its closed form, m₊ to 2e-8). The growth rate of the unstable mode is λ(τ) = d ln|r_min(p*ₕ + 1e-9) − r_min(p*ₕ − 1e-9)|/dτ
  on the slices both runs have. The normal field's charge enters through (∂φ_k)² = −2e^f φ₊φ₋ (static: q_k²/(N rⁿ)²); the linear radial sector of ghost + normal static states is the pure-ghost problem with the net charge
  (the orthogonal combination of the two fields is a test field), so #151's T(γ₁) at the member's metric applies even if q_n ≠ 0; here q_n = 0.

## Why the diagonal reads high (planted, family S)

Family S in coordinates x⁺ → λx⁺, x⁻ → x⁻/λ (a gauge map: the metric −2e^{−f}dx⁺dx⁻ is invariant, the static throat drifts with v = (1 − λ²)/(1 + λ²)), through bisection to 1e-13, the series and the same fits (h = 0.02, D = 4):

| drift v | locus τ slope | rel. to 1/W | #40's diagonal slope (t) | rel. to 1/W |
|---|---|---|---|---|
| 0 | 0.846061 ± 0.000026 | +1.9e-5 | 0.846642 | +7.1e-4 |
| 0.1 | 0.846052 ± 0.000027 | +8.4e-6 | 0.912185 | **+7.8 %** |
| 0.5 | 0.845785 ± 0.000045 | −3.1e-4 | 1.017801 | +20 % |
| 0.13 (D = 5, h = 0.00707) | 0.578546 ± 0.000088 | −4.4e-4 | 0.645295 | +11.5 % |

A throat that drifts by 0.1 in the grid's coordinates is read 7.8 % high on #40's diagonal events, the size of G's +6.2 % at D = 4 (drift 0.06 to 0.09 at the events). On G the same runs give: locus τ slope 0.84954,
locus t slope 0.87088 (+2.9 %: the lapse e^{f/2} = 1.025 at the locus), diagonal 0.89829. Static solutions in other double-null gauges (PG: x̃⁺ = F(x⁺), f ≠ 0 and a variable drift) are recovered the same way.

## The slope in the locus's proper time

| D | h | locus τ slope ± OLS SE | collapse, expansion | fit spread | rel. to 1/W | locus t slope | diagonal t (#40) |
|---|---|---|---|---|---|---|---|
| 4 | 0.02 | 0.849542 ± 0.000249 | 0.849923, 0.849160 | [0.848178, 0.852161] | +4.13e-3 | 0.870881 | 0.898291 |
| 4 | 0.01 | 0.849536 ± 0.000247 | 0.849925, 0.849146 | [0.848179, 0.852151] | +4.13e-3 | 0.870901 | 0.898043 |
| 4 | 0.005 | 0.849489 ± 0.000231 | 0.849719, 0.849259 | [0.848180, 0.851639] | +4.07e-3 | 0.870859 | 0.899089 |
| 5 | 0.00707 | 0.701666 ± 0.005193 | 0.715935, 0.687397 | [0.665555, 0.753408] | +2.12e-1 | 0.714683 | 0.798948 |
| 5 | 0.00354 | 0.702232 ± 0.005234 | 0.716598, 0.687867 | [0.665883, 0.754458] | +2.13e-1 | 0.715293 | 0.800044 |
| 5 | 0.00177 | 0.702185 ± 0.005149 | 0.715880, 0.688490 | [0.665997, 0.752506] | +2.13e-1 | 0.715250 | 0.800277 |

Three different error bars: the OLS SE (statistical), the fit spread (window and side choice) and the resolution difference (5.5e-5 and 6.7e-5 between the two finest h). The alternative levels
(r_min = 0.75, 3) give 0.84992, 0.84991, 0.84987 (D = 4) and 0.7140, 0.7146, 0.7146 (D = 5), inside the spreads (ED).

At D = 4, 1/W lies 4e-3 below the primary slope (outside the spread by 2e-4 at the finest h); the same fits give +1.9e-5 for the planted S at v = 0, so this is the member's: r_th = 1.0028 and γ₁ = 0.048 give 0.8472 and
the lifetimes are read at the departure time (0.8465, 0.8481, 0.8510, 0.8531 over e = 3 to 5, 5 to 7, 7 to 9, 9 to 11 on the collapse side at h = 0.005).
**At D = 5 the log law is not straight:** the local slope per e-fold rises from 0.644 to 0.788 over e = 3 to 11 (collapse side, h = 0.00177; expansion side 0.628 to 0.759), the collapse and expansion slopes differ by 4 %,
and the fits average over it. The pair rate says why (next section).

## The member at the stall

Profile slices of the pair p*ₕ ∓ 1e-9 (h = 0.005, D = 4; `throat_locus_output.txt` has every slice, h and D):

| t | τ | r_th | l | v | f | s_th (γ₁) | M₊, M₋ at 2 r_th (symmetric 1/(4r_th)) | q_n²/q_g² | 1/λ | T(γ₁) r_th |
|---|---|---|---|---|---|---|---|---|---|---|
| 11 | 10.77 | 1.0038 | +0.475 | +0.0763 | +0.0501 | 0.0485 | 0.2669, 0.2334 (0.2491) | +7.9e-6 | 0.8484 | 0.8479 |
| 12 | 11.74 | 1.0031 | +0.552 | +0.0780 | +0.0485 | 0.0484 | 0.2667, 0.2333 (0.2492) | −8.0e-6 | 0.8478 | 0.8474 |
| 13 | 12.71 | 1.0028 | +0.631 | +0.0794 | +0.0470 | 0.0483 | 0.2665, 0.2333 (0.2493) | −8.0e-6 | 0.8475 | 0.8472 |
| 14 | 13.69 | 1.0028 | +0.711 | +0.0806 | +0.0457 | 0.0482 | 0.2665, 0.2333 (0.2493) | −1.2e-6 | 0.8476 | 0.8471 |

- **D = 4 (a plateau).** r_th relaxes from 1.024 (t = 5) to 1.0027 (t = 15), then leaves; the fit of s_th has an rms of 5e-5 or less over t = 11 to 15, and γ₁ = 0.0485, 0.0484, 0.0483 at t = 12, 13, 14 at each of the three h. The member's asymptotic masses are m₊ = +0.0449 and m₋ = −0.0523 (the fitted
  member's M(r) extrapolated in 1/r), γ₁ = 0.0483; T(0.0483) = 0.84478 (#151's method A at the member's γ₁), times r_th 1.0028 = 0.8471. The pair rate is 1.1799 against W = 1.18197 (the symmetric unit-radius value; the member's own
  1/(0.8471) = 1.1805). CMP1: 1/λ against the member's e-fold at t = 12, +0.000 (bar 2 %); CMP2: the series slope 0.8495 against it, +0.002; CMP3: over t = 7 to 15 the two agree within −0.0034 to +0.0012.
  DR: the throat radius changes by −0.23 % from t = 10 to 15. The normal field has no charge: |q_n²/q_g²| ≤ 8.1e-6 for t = 11 to 15 (1e-3 to 4e-3 at t = 5 to 6 while the pulse leaves), against the brief's remark that G's end state may carry one.
- **D = 5 (a drifting scale).** At h = 0.00177: r_th = 1.059 (t = 6), 1.088 (8), 1.129 (10), 1.175 (12), 1.197 (13), 1.213 (14); s_th = 0.12 to 0.17 (the member is moderately asymmetric; T(γ₁) is not available at D = 5, so the member's e-fold is the symmetric r_th/W);
  m₊ = +0.077, m₋ = −0.097 at t = 8 (fitted member). 1/λ against r_th/W: 0.981, 0.992, 0.995, 0.995, 0.996, 0.997, 0.998, 1.003 at t = 6, …, 13 (CMP3: −0.0051 to +0.0033 over t = 7 to 13). DR: +11.7 % from t = 7 to 13.
  CMP1 at t = 6: 1/λ = 0.6013 against 0.6130, **−1.9 %** (bar 2 %; the same at the three h). CMP2: the series slope 0.7022 against the member's 0.6130 at the stall time, **+14.6 %, not met** (all three h); the member's e-fold over t = 6 to 14 runs from 0.613 to 0.702,
  so the slope sits at its t = 14 end. The pair rate is not exactly the symmetric member's (the asymmetry would lower it by a percent or two if D = 5 behaved like D = 4: T(0.15)/T(0) ≈ 0.985), and that is not separated here.

## Expanding side on the locus (#213's observer and classes, the locus's proper time)

Rectangle runs to X = 200 (D = 4, t up to 283) and 100 (D = 5), G with A = p*ₕ + 1e-6 ("tuned"), 0.05 and 0.2 (D = 4), tuned and 0.05 (D = 5); the control is S (A = 0.05) through the same code.

| member | D | h | class | k | τ_s | H/W at r/a = 10, 25, 50 (locus; diagonal) | drift v on the locus; l at the end |
|---|---|---|---|---|---|---|---|
| S, A = 0.05 (control) | 4 | 0.02 | FT | 2.00004 | 8.43401 | 1.280, 2.023, 2.861 (same) = #213's | 0; 0 |
| G tuned | 4 | 0.02, 0.01 | FT | 1.9998 | 18.5102, 18.5103 | 1.226, 1.942, 2.746 (1.201, 1.930, 2.740) | +0.0404; +11.4 |
| G 0.05 | 4 | 0.02, 0.01 | FT | 1.9999 | 9.5141, 9.5125 | 1.230, 1.946, 2.753 (same) | +0.0244; +6.8 |
| G 0.2 | 4 | 0.02, 0.01 | FT | 1.9993 | 7.3742, 7.3740 | 1.424, 2.254, 3.188 (1.421, 2.249, 3.184) | −0.149; −36.5 |
| G tuned | 5 | 0.00707, 0.00354 | FT | 3.9834 | 17.0011, 17.0092 | 0.681, 0.865, 1.030 (0.657, 0.858, 1.027) | +0.0234; +4.2 |
| G 0.05 | 5 | 0.00707, 0.00354 | GROW | 3.9883 | 9.1260, 9.1196 | 0.871, 1.105, 1.315 (same) | +0.0144; +2.0 |

- Every member has w̄ = (1 − n)/(1 + n) (−1/3, −1/2) to 4 digits and |ρ|/H² → 4.00 (D = 4); #213's identities are the locus's too, with the velocity-dependent frame (PO checks the frame algebra to 2e-16).
  k and |ρ|/H² agree across the members (D = 4: k 1.9993 to 1.9999); H/W at fixed r/a depends on A (D = 4, r/a = 10: 1.226 tuned, 1.230 for A = 0.05, 1.424 for A = 0.2), and the tuned and A = 0.05 members are within 0.4 % of
  #213's S tuned member (1.230, 1.946, 2.752). A = 0.2 puts the locus 36 proper lengths from the diagonal at t = 283 (v = −0.149), and k changes by 7e-4.
- The slope test of #211 (H/W within 10 % of 1 at r/a = 10, 25, 50) is not met by any run, as for S. 1/H is linear in τ to 2e-5 (D = 4) and 3e-4 (D = 5).
- **KX is not met as declared**: D = 5 A = 0.05 is class GROW, because the series ends (the locus reaches the edge x⁺ = X) at τ_end = 7.59 with τ_s = 9.13, a ratio of 1.2025 against the declared bound 1.2 of #213's FT test; its 1/H is linear to 2e-4 and k = 3.9883 (S: 3.9897 at the same D). Not changed.
  The D = 5 resolution check RES passes (H at fixed r/a within 8e-4, k 3.9883/3.9883 and 3.9834/3.9834).

## Checks and mutations

59 checks, 55 pass; the four that fail are reported as results: **CMP2 at D = 5 (three h)** and **KX**. PB, PG, PM, PO, PL, PS, JOB, CTL, FA, ED, RA, CMP1 (D = 5 at −1.9 %, inside 2 %), CMP3, DR, RES and XC pass.
Each of the 20 mutations fires its declared set on the stored record, plus the two checks that already fail without a mutation (CMP2, KX), and for the 15 that act on a run each fires exactly its declared set on the reduced suite
(`runs/mut/`, `runs/mut_mini/`); the targets were computed in development before they were declared (MUTATIONS in `throat_locus.py`; `finder_sign` is silent on a stored record by design). XC: S through the locus code against #213's stored
diagonal run: k 2.00004 / 2.00004, τ_s 8.43401 / 8.43404, H/W to 1.5e-6.

Cost (dh610, load 12 to 24 of 20, one thread each, `nice -n 10`, four lanes): 31 jobs in 64 min wall (09:07 to 10:11), the longest 44 min (D = 5, h = 0.00177, one of four parts), peak memory 274 MB, none near the 2 h stop;
the analysis takes 2 min. The sweep stores no grid.

## Discrepancies and things that failed

1. The brief's "dl/dt = +0.58, +0.25" are positions (above).
2. CMP2 at D = 5 and KX (above). The criteria were not moved.
3. The brief's remark that G's static end state may carry a normal charge: not seen at D = 4 (≤ 1e-5), and at D = 5 the ratio is 2.6e-2 at t = 5, 7e-3 at t = 6, −6e-3 at t = 7 (the pulse's remnant) and below 1e-3 from t = 8. The pure-ghost T(γ₁) applies at the member's metric in any case (above).
4. During development: a first D = 5 series used h rounded to six figures (p*ₕ ∝ h², 1e-9 off; the log was overwritten by the rerun, the rule's slip); `points40` read #40's `ed:` points into its dictionary and the control CTL failed on 8 of 40 points at D = 4,
   h = 0.01 until fixed; the development units were stopped once to add two checks at the coordinator's ruling. All kept in the raw README.
5. The brief's planted velocity (the stall table's 0.58) would have been ≈ 8× the measured drift; v = 0.1 was kept (v = 0.5 for the lever arm).

## What `threshold.md`, `critical_ghost.md` §1.3 and #42 should now say (proposals for the coordinator)

`threshold.md`, "What the near-critical run looks like" and "Limits": the stall table's `l` and `dl/dt` columns are a position and a drift (this note's table of members gives them at the pair of runs). Replace "the lifetime of G is that of an observer that is not at the throat … not separated here"
by: G's lifetime in the proper time of the minimum-r locus has slope 0.8495 at D = 4 (+0.4 % over 1/W, two finest h agree to 5.5e-5), which is the e-fold of the static member at the stall (γ₁ = 0.048, r_th = 1.0028, no normal charge,
T(γ₁) r_th = 0.847); #40's +6 % is the diagonal observer (a planted drift of 0.1 reads +7.8 %). At D = 5 the locus slope is +21 % over 1/W and the throat radius drifts (1.07 to 1.20 over t = 7 to 13) with the mode's growth rate following W/r_th(t); the
slope is not a single e-fold. D = 6 was not run.

`critical_ghost.md` §1.3, the bullet "A generic one-sided family (G…) gives slopes 6 %, 38 % and 115 % above 1/W … Its throat drifts off the line where the lifetime is read, and the cause is open (#225, #42)":
at D = 4 the cause is the observer: with the lifetime in the locus's proper time the slope is 0.8495 and the near-critical state is the member γ₁ = 0.048, r_th = 1.003 (m₊ = +0.045, m₋ = −0.052 in the code's units, a = 1); at D = 5 the state's throat radius grows
with time (a drifting scale) and the series' slope averages over it (*numerical evidence*, `throat_locus.md`). The last paragraph of "A difference from textbook type I" ("which member generic data select … has not been studied") can say: for the family G at D = 4 the member has γ₁ = 0.048, r_th = 1.0028 a.

**#42:** its first acceptance bullet asks for the throat radius and asymmetry of the member two families approach. S (and O) approach γ₁ = 0, r_th = 1.0000 a (the mirror symmetry, *proved* in #213, and #40's stall table); G at D = 4 approaches γ₁ = 0.048 ± 0.001, r_th = 1.0028 a (the fit's rms ≤ 5e-5 over
t = 11 to 15; the h = 0.02, 0.01, 0.005 pairs agree to 3 digits) with the e-fold of #151's T(γ₁), and at D = 5 a member of growing radius (s_th 0.12 to 0.17). The second bullet (the second collapse window) is not touched.

## Limits

- G at D = 6 was not run (#228). D = 5's member e-fold is the symmetric r_th/W: T(γ₁) is known at D = 4 only (GGS09a's family, #151/#163); the asymmetric linear mode at D = 5 is not computed (task 2).
- The pair rate is read at p*ₕ ∓ 1e-9 and is a linear-regime quantity; the comparison with the series' slope (CMP2) is a statement about an average over a drifting member at D = 5.
- The asymmetry is fitted from M(r) on two sides over r = 1.25 to 4 r_th, in the region the pulse has left (the fit's rms is 1e-5 at D = 4 and 3e-3 at D = 5); the end masses quoted are the fitted member's, not a measured asymptote (the evolution does not reach r → ∞ in the member's frame).
- The locus is slicing dependent beyond linear order: it is the throat worldline of a static member (planted checks), and on the expanding runs it is the minimum-r curve of the grid's slices, which differs from #213's geodesic diagonal (H/W agrees to 0.3 % for the tuned member).
- The planted velocities are gauge maps, not physical boosts; the planted lifetime series test the code path (events, τ, fits), not the physics of drifting members. The asymmetric static member as *data* (a planted e-fold of T(γ₁) r_th) is not built (task 2).
- p*ₕ is #40's (D = 5, h = 0.00707: its bracket was 2.9e-12; the series is intact to e = 11 and the fits ignore the floor).
- #41's settled-mass chain for G was not run here (the chain is on main, not in this worktree, and needs its own p*ₕ and controls); recommended task 3.

## Recommended tasks

1. **The settled mass of G (#41's note on #225).** What: #41's chain (`settled_mass.py`) for family G at D = 4, h = 0.02, 0.01, 0.005 (and 0.0025 if p*ₕ is bisected there), with its planted cases unchanged. Why: the near-critical G state is the member γ₁ = 0.048,
   r_th = 1.0028 (this note), so a settled mass different from S's 0.2231 a would not come from the static member. Acceptance (can fail): M_AH(t) second order at three h; PL1–PL3 unchanged; the value compared with 0.2231 a. Depends on: this note, #41.
2. **The asymmetric static member as data, with its e-fold.** What: build the member (r_th, s_th) as characteristic data (the lapse f = −2 ln N on the rays), evolve it at D = 4 (γ₁ = 0.5 and the member found here) and D = 5, tune with a one-sided pulse and measure λ on the locus.
   Why: it plants #151's T(γ₁) r_th in the evolution (#39's recommended task 1, #42), and supplies T(γ₁) at D = 5, where CMP1 sits at −1.9 %. Acceptance: static at second order; λ within 1 % of T(γ₁) r_th (D = 4); at D = 5 λ against W/r_th and the asymmetric correction reported. Depends on: this note.
3. **The scale drift at D = 5.** What: find what sets dr_th/dt ≈ +0.02 per unit t at D = 5 (and the relaxation at D = 4) from the pulse's normal and ghost energy, and whether it persists as h falls and as B changes. Why: it makes the lifetime a non-log law at D = 5 and 6.
   Acceptance: r_th(t) of the member against a prediction from the data's energies, or the dependence on B (0.04, 0.1, 0.2) measured at D = 5. Depends on: this note.
4. **G at D = 6 on the locus** (with #228's third h): the locus slope and the member's drift at D = 6. Acceptance: the pair rate against W/r_th(t); RA at the two finest h. Depends on: #228.
5. **The class bound of the expanding side.** What: rerun D = 5 A = 0.05 with X large enough for the locus to reach τ_s/τ_end ≤ 1.2 (X ≈ 160), to see whether KX's FT class holds. Acceptance: class and k at two h. Depends on: this note.

## Files

`throat_locus.py`, `throat_locus_output.txt`, this note; raw `negative-mass-raw/I225-dh610` (README with the criteria and the development runs, `dev/`, `runs/` with the series, stall pairs, expanding runs and planted series, `scripts/`).

## Coordinator's review

Review of branch `wip/i225-throat-locus` at ed26e65 and raw `I225-dh610` at 982783a, recorded in negative-mass-raw
`I225-mac-review`. On the Mac, nice 10.

- **Scope.** The branch adds `throat_locus.py`, its output and this note, nothing else.
- **Rerun.** `analyze` over the runner's stored runs gives 55 of 59: CMP2 at D = 5 three times, and KX. Every line of the
  output is reproduced apart from the header.
- **Independent check** (`ind225.py`). T(γ₁) is even in γ₁. An even quartic through three published values gives
  T(0.048) = 0.84510 and e-fold 0.84747 at r_th = 1.0028: 2.7e-4 from this note's 0.8472, and 3.4e-5 from the measured
  1/λ. The three values are T(0) = 1/W (#195), T(1/2) = 0.758357 and T(1) = 0.675, and none of them is #151's code at
  γ₁ = 0.048.
- **The brief's error.** The brief quoted the stall's positions as drift speeds. This note's section "An error of the brief"
  is right; the correction is now in `threshold.md`.
- **Applied.**
  - `threshold.md` carries an update from #225 under the stall table.
  - `critical_ghost.md` §1.3 has the D = 4 cause and the D = 5 drift.
  - #42 has the member G approaches at D = 4.
- **Follow-ups.**
  - #253: the settled mass of G (task 1).
  - #254: the asymmetric member as data, with its e-fold at D = 4 and 5 (task 2).
  - #255: the scale drift at D = 5, and the class bound of the expanding side (tasks 3 and 5).
  - Task 4 is added to #228.

# Family G at D = 6 at its third resolution, on the throat locus, and its drift (#228, #225, #255)

#40's family G (a ghost pulse A tuned to p*ₕ plus a normal pulse B = 0.1 on the same ray) at D = 6, run at h = 0.00125 with #40's job, read on the
throat locus with #225's code, and its drift measured with #255's. Driver: [g_d6.py](g_d6.py), a thin driver that imports `threshold.py`, `throat_locus.py`
and `throat_drift.py` unchanged and adds the D = 6 analysis. Outputs of record: [threshold_output.txt](threshold_output.txt) (regenerated, `python3 g_d6.py
threshold-analyze --raw=<negative-mass-raw>/I228-dh610/runs40`) and [g_d6_output.txt](g_d6_output.txt) (`python3 g_d6.py analyze --base=<raw>/I228-dh610/runs --hks=1,2`,
exit 1: 12 of 41 checks fail, listed below). Raw output, with the criteria written before the runs and the settings before the analyses: negative-mass-raw `I228-dh610`.
Status words as in AGENTS.md: every number below is *numerical evidence*. No sentence says a result is new.

## Answer

The third h runs in 1.9 h, not 4 h, and closes the three-h order: p*ₕ = 1.265768e-2, 1.108798e-2, 1.069693e-2 at h = 0.005, 0.0025, 0.00125 has order 2.005
(band [1.75, 2.25]) and Richardson p* = 1.05672e-2, and equals #246's bisection of the same h bit for bit. G's diagonal slope is 0.990046, 0.994756, 0.996251
(the two finest differ by 1.5e-3, bar 5e-3), +116 % over 1/W = 0.461924, so the departure is resolution independent at D = 6 as at D = 4 and 5.

**On the locus, G at D = 6 has no lifetime law at #225's levels.** The locus's r_min reaches 0.5 on the collapse side only for e ≤ 4 (3 of 19 points), and for
e ≥ 10 both sides cross r_min = 2 at the same τ = 15.75 ± 0.03: the tuned throat's own radius grows through 2. FA and ED fail, and RA (7.1e-4, expansion side alone 1.3e-3)
compares numbers that saturate. What the locus does give is the mechanism: the unstable mode's growth rate follows r_th(t)/W to −0.8 … +1.5 % over t = 7 to 11 (CMP3, both h)
while r_th grows by 26 % (1.267 → 1.593) at a speed ṙ = dr_th/dτ ≈ 0.09. The collapse-side lifetime is the time at which the mode reaches 10 % of r_th plus −0.42 ± 0.02
(e ≥ 4.5), and its slope per decade rises from 0.71 to 1.05 over e = 3 to 11 as r_th/W does: G's fitted 0.99 is that average, not one e-fold.

**The drift** is not a constant speed at B = 0.1: ṙ falls 0.0938, 0.0872, 0.0790 at t = 8, 10, 12 (two h agree to 1e-5, window part 5e-4) and #255's fit does not converge (κ → 0).
Where it converges (B ≤ 0.07) ṙ∞ ∝ B^2.7 (0.02 to 0.04) and B^2.4 (0.04 to 0.07), against B^3.16 at D = 5; the local ṙ(8) ∝ B^2.7 (0.04 to 0.07), B^2.2 (0.07 to 0.1) and B^1.1 (0.1 to 0.2).
ṙ∞ / E_net^(3/2), constant at D = 5 (0.110), is 0.136, 0.112, 0.083 at B = 0.02, 0.04, 0.07 (post hoc).

| acceptance (issue #228 and the task) | verdict | numbers |
|---|---|---|
| the job finishes, bisection to 1e-14, no run static through X (BR, FM) | met | 44 runs, hi − lo 5.68e-15, 6854 s; 0 of 34 series points off their side |
| p*ₕ of the job against #246's (P246) | met | `lo`, `hi`, p*ₕ equal bit for bit at h = 0.00125 (and at 0.0025, 0.005) |
| p*ₕ's order over the three h in [1.75, 2.25] (PS) | met | 2.005; p* = 1.0567197e-2 |
| slope at the two finest h within 0.5 % (RA) | met | 0.994756, 0.996251; 1.50e-3 |
| 1/W inside the fit spread (SL, reported) | not met, as at the two coarser h | spread [0.798695, 1.209563], +1.16 |
| RA of the locus τ slopes (LRA) | met in form | 0.550627, 0.550239: 7.06e-4; expansion side alone 1.29e-3 (LRAx); FA and ED fail, so the slope is not a lifetime law |
| the pair rate against r_th(t)/W (CMP) | CMP3 met, CMP1 not met | t = 7 … 11: −0.8 … +1.5 %; t = 6: −3.1 % (bar 2 %) |
| the drift speed and its B exponent, each with its error | measured | tables below |

## What ran

dh610-fedora, one thread per run, `nice -n 10`, at most 4 threads, systemd units `nm-i228-*`, 2026-11-27 23:39 to 2026-11-28 03:40.

| run | h | cost (wall, one thread) |
|---|---|---|
| job: bisection, series (`threshold.py job`) | 0.00125 | 6854 s (bisection 3511 s, 44 runs) |
| tail: `constraint_runs`, `stall_runs` (`g_d6.py tail`) | 0.00125 | 903 s, 32.1 GB peak |
| locus series, 4 parts each, and the stall pair at p*ₕ ∓ 1e-9 | 0.0025; 0.00125 | 8 to 16 min per part and 2 min; 58 to 70 min per part and 17 min |
| drift bisections B = 0.04, 0.07 (1e-12) | 0.00125 | 3228 s (37 runs), 3937 s (38 runs) |
| drift bisections B = 0.02, 0.04, 0.07 | 0.0025 | 1028 s, 965 s, 954 s |
| drift bisection B = 0.2 | 0.0025 | stopped by its 2 h limit after 36 runs (rc 124): bracket 7.0e-12 |
| drift pairs p*ₕ ∓ 1e-10 | B = 0.02, 0.04, 0.07, 0.1, 0.2 at 0.0025; 0.04, 0.07, 0.1 at 0.00125 | 2 to 4 min; 10 to 23 min |
| development at h = 0.005 | 0.005 | not a run of record (`dev/`) |

The tail is #40's `constraint_runs` and `stall_runs` on the job's own bracket, run after the job so that a failure in the 32 GB step could not lose the bisection. On #40's own
D = 6, G, h = 0.005 bracket it writes `constraints.json` and `stall.json` byte for byte equal to #40's. `threshold.py analyze` skips this job by its declaration `not_run={(4, "G"): 2}`; the driver clears
that entry and calls `threshold.main()`: on #40's raw it reproduces `threshold_output.txt` line for line (header and timing aside), and on the merged record it differs in the D = 6 G rows, the tally and the timing only.

## p*ₕ, its order, and the diagonal slope (#40's quantities)

| h | p*ₕ | primary slope ± OLS SE | collapse, expansion | fit spread | rel to 1/W | M at first trapped cell (r) | max \|C\| |
|---|---|---|---|---|---|---|---|
| 0.005 | 1.2657679e-2 | 0.990046 ± 0.026391 | 0.854538, 1.125554 | [0.795196, 1.190386] | +1.14 | 2.57667 (1.727) | 7.371e-5 |
| 0.0025 | 1.1087978e-2 | 0.994756 ± 0.026714 | 0.857876, 1.131636 | [0.798040, 1.200660] | +1.15 | 2.60758 (1.734) | 1.828e-5 |
| 0.00125 | 1.0696934e-2 | 0.996251 ± 0.027127 | 0.857151, 1.135351 | [0.798695, 1.209563] | +1.16 | 2.61484 (1.736) | 4.294e-5 |

- Order 2.005, p* = 1.056720e-2 (D = 4: 1.4855e-2, D = 5: 1.0681e-2). The mass at the first trapped cell converges at second order too (differences 0.03091, 0.00726: ratio 4.26).
- The constraint maximum falls by 4.03 over the first halving and then rises (ratio 0.426): the roundoff plateau of #226 at the finest h, not truncation.
- The series is intact to e = 11 (0 of 34 off, 0 of 4 off below that). The collapse- and expansion-side slopes differ by 32 % at every h, as at the two coarser h: the two events are read at different phases of the departure (below).
- Stall of the two runs at p*ₕ ∓ 1e-8 (h = 0.00125; the same to 3 digits at the other h): collapse t = 11.92, r = 1.656, l = +1.308, v = +0.111; expansion t = 6.24, r = 1.209, l = +0.431, v = +0.2005.

## The locus (#225's code, h = 0.0025 and 0.00125)

**The controls pass.** CTL: all 40 points of each series (38 and the two bracket ends) reproduce #40's event, cell and time bit for bit and its diagonal slope 0.994756, 0.996251; CTLN: all 40 amplitudes are points of
#40's record (so CTL is not vacuous). PL1_6 and PL2_6 (planted static throat of radius 1.2 and moving throats, at n = 4, h = 0.005 … 0.00125): 2.2e-5 and 6.7e-5 against the bars 2e-4 and 5e-4, order 2.00.

**The levels.** At the primary levels r_min = 0.5 (collapse) and 2 (expansion):

| e | collapse side: locus event | τ(r_min = 2), collapse side | τ(r_min = 2), expansion side | T (diagonal) collapse, expansion |
|---|---|---|---|---|
| 3 | r_min = 0.5 at τ = 6.12 | none | 6.47 | 4.99, 7.52 |
| 5 | none | none | 9.11 | 8.27, 11.94 |
| 7 | none | none | 12.07 | 11.90, 16.95 |
| 9 | none | none | 14.80 | 16.11, 22.32 |
| 10 | none | 16.19 | 15.53 | 18.41, 25.21 |
| 11 | none | 15.77 | 15.72 | 20.88, 28.42 |
| 12 (floor) | none | 15.74 | 15.73 | 23.49, 29.65 |

The collapse side leaves the slice before the locus reaches 0.5 (#40's trapped cell forms at r ≈ 1.7), so the sweep runs to the end of the triangle (5 min per point at h = 0.00125). For e ≥ 10 the locus
crosses 2 at τ = 15.75 ± 0.03 on both sides: it is the time at which the tuned throat reaches 2, not a function of p − p*ₕ. Hence FA fails (14 of 34 points off at h = 0.00125, 15 at 0.0025) and ED fails
(alternative levels 0.75, 3 give 0.8146, against the primary spread [0.4529, 0.6201]). The primary locus τ slope is 0.550627 ± 0.018349 and 0.550239 ± 0.017929 (+19.2 %, +19.1 % over 1/W), and 0.550556 and 0.549844
on the expansion side alone. RA and LRAx pass, and CMP2 passes (series slope 0.5498 against the member's r_th/W = 0.5507 at t = 6), but those agree because the values saturate, not because the slope is the e-fold.

**The pair rate follows r_th(t)/W on a drifting throat** (pair p*ₕ ∓ 1e-9, h = 0.00125; the member is the symmetric one, e-fold r_th/W; h = 0.0025 agrees to 3 digits):

| t | τ | r_th | v | f | 1/λ | r_th/W | 1/λ ÷ (r_th/W) |
|---|---|---|---|---|---|---|---|
| 6 | 5.62 | 1.1921 | +0.202 | +0.295 | 0.5337 | 0.5507 | 0.969 |
| 7 | 6.48 | 1.2671 | +0.191 | +0.235 | 0.5805 | 0.5853 | 0.992 |
| 8 | 7.37 | 1.3501 | +0.173 | +0.194 | 0.6233 | 0.6236 | 1.000 |
| 9 | 8.27 | 1.4335 | +0.154 | +0.171 | 0.6661 | 0.6622 | 1.006 |
| 10 | 9.18 | 1.5149 | +0.137 | +0.157 | 0.7075 | 0.6998 | 1.011 |
| 11 | 10.10 | 1.5930 | +0.122 | +0.146 | 0.7468 | 0.7358 | 1.015 |
| 12 | 11.02 | 1.6674 | +0.110 | +0.138 | 0.7843 | 0.7702 | 1.018 |
| 14 | 12.89 | 1.8015 | +0.091 | +0.114 | 0.8576 | 0.8321 | 1.031 |
| 16 | 14.81 | 1.8766 | +0.078 | +0.006 | 0.9546 | 0.8669 | 1.101 |

CMP3 (t = 7 to 11, bar 2 %): −0.8 % to +1.5 %; DR: r_th grows by 25.7 %. CMP1 at t = 6 fails (−3.1 %, bar 2 %; −3.3 % at h = 0.005): the pulse has just left (the normal field's charge ratio is −1.6e-2 there, 8e-5 at t = 7).
The ratio drifts upward beyond t = 12, where the pair at 1e-9 stops being linear (|odd| = 7e-4 at t = 12). The normal field carries no charge from t = 7 (|q_n²/q_g²| ≤ 1.6e-3).

**What sets the slope of T** (post hoc, gates nothing; `PH` lines of the output). By linearity the deviation of the run at p*ₕ + δ from the tuned throat is δ/δ₀ times the pair's odd part, so its departure time t_dep (the
deviation reaches ε = 0.1 of r_th) is where |odd| = ε δ₀/δ. Against #40's collapse-side lifetimes of the same job:

| e | t_dep | τ_dep | r_th(τ_dep) | r_th/W | T(collapse) − t_dep | T(expansion) − t_dep |
|---|---|---|---|---|---|---|
| 3 | 5.60 | 5.29 | 1.165 | 0.538 | −0.62 | +1.92 |
| 5 | 8.68 | 7.98 | 1.407 | 0.650 | −0.41 | +3.25 |
| 7 | 12.32 | 11.32 | 1.691 | 0.781 | −0.41 | +4.63 |
| 9 | 16.55 | 15.28 | 1.959 | 0.905 | −0.43 | +5.78 |

- For e ≥ 4.5 the collapse-side lifetime is t_dep − 0.42 ± 0.02 (ε = 0.1; ε = 0.03 and 0.3 shift the constant to +0.4 … +0.8 and −1.5 … −1.2 with a spread of 0.4 and 0.3, so the constant is a choice of ε and its constancy is ε-dependent).
- The slope per decade of T(collapse) over e = 3–5, 4–6, 5–7, 6–8, 7–9 is 0.713, 0.740, 0.788, 0.850, 0.914, and that of t_dep is 0.669, 0.724, 0.789, 0.856, 0.919; that of τ_dep is 0.584, 0.653, 0.724, 0.795, 0.861,
  against r_th/W at the mid 0.594, 0.652, 0.716, 0.781, 0.843 (the pair's own rate, within 2.1 %).
- The expansion event on the diagonal lags the departure by a time that grows from 1.9 to 5.8: the throat has moved off the diagonal (l ≈ +1.3) by then.

So G's fitted slope is an average over a rising local slope, r_th(τ_dep)/W × dt/dτ, with r_th(τ_dep) growing 1.2 → 2.0 over e = 3 → 9. The cause of the growth of r_th is not touched here (#255, #273).

## The drift (#255's code, B = 0.04 … 0.2)

Pairs A = p*ₕ ∓ 1e-10 kept whole; the pair's mean is the tuned member. ṙ is dr_th/dτ by a linear fit over ±0.5 in τ (±0.25 and ±1 give the window part). Finest h at each B (0.0025 at B = 0.02 and 0.2):

| B | h | r_th − 1 at t = 8, 10, 12 | ṙ at t = 8, 10, 12 | fit ṙ∞ ± SE, κ | max \|odd\| (NO6) |
|---|---|---|---|---|---|
| 0.02 | 0.0025 | 0.00635, 0.00960, 0.01181 | 0.00154, 0.00170, (−0.0067) | 0.00173 ± 7e-6, 0.72 | 2.4e-2 |
| 0.04 | 0.00125 | 0.03611, 0.05590, 0.07745 | 0.00926, 0.01063, (0.00968) | 0.01144 ± 6e-5, 0.47 | 1.1e-2 |
| 0.07 | 0.00125 | 0.14970, 0.23362, 0.31877 | 0.04250, 0.04470, 0.04393 | 0.04487 ± 3e-5, 1.10 | 1.1e-3 |
| 0.1 | 0.00125 | 0.35006, 0.51492, 0.66809 | 0.09382, 0.08723, 0.07901 | not converged (κ → 0) | 6.9e-5 |
| 0.2 | 0.0025 | 0.88638, 1.18175, 1.45025 | 0.20289, 0.17524, 0.15140 | not converged (κ → 0) | 9.3e-7 |

Values in brackets are beyond the pair's linear window (NO6 fails for B ≤ 0.07, the unstable mode grows faster at a small radius); the fits mask the points above 3e-4 (t to 9.9, 10.2, 11.2). Errors:

- window part of ṙ (±0.25 and ±1 against ±0.5): ≤ 5e-4 at B = 0.1, ≤ 1.5e-4 at B = 0.07, ≤ 2e-5 at B = 0.04 (t = 8, 10);
- resolution (h = 0.0025 against 0.00125): ≤ 1e-5 at B = 0.1, 4e-5 at B = 0.07, 1.0e-4 at B = 0.04;
- the fit's SE: as in the table (a statistical error of an exponential on a window of 3 to 5 proper-time units; 0.7 % of ṙ∞ at B = 0.04).

**B exponents** (a ∝ B^α between neighbouring B, same h; the error of ṙ∞'s exponent is the two fits' SE propagated; the window and resolution parts are those of ṙ above):

| B | h | p*ₕ (the extrapolated p* in brackets) | r_th − 1 at t = 8, 10, 12 | local ṙ at t = 8, 10, 12 | ṙ∞ |
|---|---|---|---|---|---|
| 0.02 → 0.04 | 0.0025 | 1.098 | 2.47, 2.51, 2.69 | 2.58, 2.63, n/a | 2.725 ± 0.010 |
| 0.04 → 0.07 | 0.00125 | 2.035 (2.157) | 2.54, 2.56, 2.53 | 2.72, 2.57, (2.70) | 2.443 ± 0.010 |
| 0.07 → 0.1 | 0.00125 | 2.388 (2.438) | 2.38, 2.22, 2.08 | 2.22, 1.88, 1.65 | not converged |
| 0.1 → 0.2 | 0.0025 | 2.792 | 1.34, 1.20, 1.12 | 1.11, 1.01, 0.94 | not converged |

The resolution part of the local-ṙ exponents is −0.016, −0.013, +0.013 (0.04 → 0.07) and −0.003, −0.002, −0.002 (0.07 → 0.1); the window part is ≤ 0.004 at t = 8 and 10.
p*ₕ at B ≠ 0.1 moves with h (order 2.008, 2.007 at B = 0.04, 0.07 over h = 0.005, 0.0025, 0.00125; extrapolated p* = 1.32443e-3, 4.42972e-3), so the exponent of the raw p*ₕ carries an h² offset and the extrapolated one is the better estimate.

- **Against D = 5** (#255: p*ₕ ∝ B^2.04–2.11, r_th − 1 ∝ B^2.6–2.9, ṙ∞ ∝ B^3.16 over B = 0.04 to 0.1): at D = 6 the threshold exponent is 2.27 (extrapolated, 0.04 to 0.1), r_th − 1 goes as B^2.5 (0.04 to 0.07) and B^2.4 to 2.1 (0.07 to 0.1), and ṙ∞ as B^2.4–2.7 below B = 0.07.
  The local exponent of ṙ falls with B (2.7, 2.2 at t = 8 between 0.04–0.07 and 0.07–0.1; 1.1 between 0.1 and 0.2): the speed saturates at large B.
- ṙ∞ / E_net^(3/2) with E_net = (B² − p*ₕ²) I₂, I₂ = 136.5455 (n = 4; post hoc): 0.1357 (B = 0.02), 0.1122 (0.04), 0.0825 (0.07). D = 5's 0.110 constant (0.3 %) does not carry to D = 6.
- The speed is not the asymptote of a transient here: at B = 0.1 and 0.2 ṙ(τ) is a slowly falling function (0.094 → 0.079 over t = 8 to 12), and `fit_rdot` degenerates (κ → 0, ṙ∞ = −124 ± 5e4). PL3_6 (planted ṙ∞ 0.026, 0.05 with κ = 0.15, 0
  with κ = 0.3) recovers the local slope to 1e-5 and ṙ∞ within 1.8 SE, so the estimator is sound where the model holds; the model does not hold at B ≥ 0.1.

## Checks and mutations

`g_d6_output.txt`: 41 checks, 29 pass, 12 fail; the failures are results and no criterion was moved:

| check | verdict | why |
|---|---|---|
| FA, ED (#225's, h = 0.0025, 0.00125) | fail | the locus levels are not a lifetime at D = 6 (above) |
| CMP1 (both h) | fail | −3.1 % at t = 6, bar 2 % |
| NO6 at B = 0.02 (0.0025), 0.04 (both h), 0.07 (both h) | fail | max \|odd\| 2.4e-2, 1.1e-2, 1.1e-3 against 3e-4: the pair leaves its linear window before t = 12 (as #255's NO at B = 0.04 by 13 %) |
| JOB6 at B = 0.2 (h = 0.0025) | fail | the + side's #40 event is beyond the triangle t ≤ 34 (fates (−1, 0)); the bisection was stopped at 7.0e-12 |

Passing: P246 (2), JOB (2), CTL (2), CTLN (2), RA, LRAx, CMP2 (2), CMP3 (2), DR (2), JOB6 (7), NO6 (B = 0.1 at both h, 0.2), PL1_6, PL2_6, PL3_6. With #40's analysis: 115 checks, 103 pass, 12 fail (SL at G for D = 4, 5, 6, ED D = 4, BR D = 5, SL D = 6 S).

Mutations (effect computed in the development runs, then repeated on the frozen record; `runs/mut/`): `p246_off` → P246; `ctl_shift` → CTLN; `odd_scale6` → NO6; `ra_scale6` → LRAx; `tau_unit` → PL3_6; `rth_diag` → PL2_6; `pstar_off` → CMP2, CTL;
`locus_diag` → CMP2, CMP3; `drop_bracket` → JOB; `ra_scale` → RA. Each fires exactly its declared set on the record. One declared target was wrong: `rth_diag`, declared JOB6, NO6, PL2_6 from a development run in which the data
grew between the baseline and the mutated analysis; only PL2_6 changes on the frozen record (JOB6 and NO6 do not read the radius). #40's own mutations through the wrapped analysis of the merged record are not interpretable
(the record already fails SL, BR, ED) and are not claimed.

## Discrepancies and things that failed

1. `threshold.py analyze` cannot take this job by its own declaration (`not_run`); the driver clears it. The script is read-only here.
2. The job ran without `--cons --stall`; the tail ran after it (32.1 GB peak, 15 min), byte-identical to #40's on #40's own bracket. #40's flags in one job would have put the 32 GB step after 1.9 h of results.
3. The job took 1.9 h against the issue's estimate of about 4 h; the 5 h waiver was not used.
4. FA, ED, CMP1, NO6 and JOB6 failures above, and #255's fit not converging at B ≥ 0.1: the D = 5 declarations do not carry to D = 6 at B ≥ 0.07 or on the collapse side.
5. B = 0.2 at h = 0.0025: the bisection hit its 2 h limit (36 runs, bracket 7.0e-12 against the declared 1e-12); its pair used p* = 0.07681202535369897 from the stopped points (`pstar_stopped.json`), which straddles at δ = 1e-10. B = 0.2 at h = 0.00125 was not run (each run would take about 25 min).
6. Two pair runs (B = 0.07, h = 0.0025 and 0.00125) failed at once (rc 1): the queue's WAIT line is per lane, so another lane took the RUN line before its dependency existed. They were rerun as `_run2` (the failed logs are kept) with the dependency inside the run.
7. During development the planted evolutions PL1 and PL2 with #255's h = 0.02 … 0.005 missed their bars at n = 4 (3.6e-4 against 2e-4, second order); they are run at D = 6's own h with the bars unchanged. README part 2 records it.
8. README part 2 gave its writing time as 02:00; it was written at about 01:40 (corrected in the README; pushed before any analysis of record).
9. B = 0.02 and 0.2 have one h only in the record; the exponents that use them have no resolution part.

## What the other notes should now say (proposals; not applied)

**`threshold.md`**
- "Resolutions": replace "G at D = 6 was run at 0.005 and 0.0025 only (the third h is about 4 h, over the 2 h stop). Per the coordinator's ruling, D = 6 conclusions are drawn from h ≤ 0.0025" by: G at D = 6 was run at the third h in #228 (1.9 h in one thread);
  the D = 6 conclusions use the two finest h.
- The p*ₕ table's last row: `| 6 | G | 1.26577e-02, 1.10880e-02, 1.06969e-02 | 2.005 | 1.05672e-02 |`.
- The slope table: add `| 6 | G | 0.00125 | 0.996251 ± 0.027127 | 0.857151 | 1.135351 | [0.798695, 1.209563] | +1.16e+00 | **no** |`; the G resolution agreement is now 1.5e-3 (h = 0.0025, 0.00125) and 4.7e-3 (0.005, 0.0025).
- The threshold-mass row for G at D = 6: `2.57667, 2.60758, 2.61484 (r 1.727, 1.734, 1.736)`, second order (differences 0.03091, 0.00726, ratio 4.26); the constraint table row: `7.37e-05, 1.83e-05, 4.29e-05`, ratios 4.03 and 0.426 (the #226 plateau).
- The stall table's D = 6 G row is unchanged to 3 digits at the third h. "Limits": delete "G at D = 6 has two h; the third would be 4 h". "Recommended tasks" 4: done by #228.
- An update from #228 under the stall table: at D = 6 the locus is not a lifetime (`g_d6.md`): the throat radius grows 1.27 → 1.59 over t = 7 to 11, the pair rate follows r_th/W within 1.5 %, the collapse-side lifetime is the mode's departure time − 0.42,
  and G's slope is r_th(τ_dep)/W × dt/dτ rising from 0.71 to 1.05 over e = 3 to 11.

**`throat_locus.md`**
- "Limits", first bullet: "G at D = 6 was not run (#228)" → run in #228; at D = 6 the locus levels 0.5 and 2 do not mark a lifetime: the collapse side crosses 0.5 only for e ≤ 4 and both sides cross 2 at τ = 15.75 for e ≥ 10 (FA, ED fail; RA 7.1e-4 on saturating values).
- Recommended task 4 is done (`g_d6.md`): the pair rate against r_th(t)/W holds to −0.8 … +1.5 % over t = 7 to 11; RA at the two finest h is 7.1e-4.
- The statement "slope in the locus's proper time" for D = 6: not defined at #225's levels; see the departure table of `g_d6.md`.

**`throat_drift.md`**
- "Limits", last bullet: "D = 6 was not run (#228)" → measured in #228: the throat radius at B = 0.1 grows at ṙ = 0.094, 0.087, 0.079 (t = 8, 10, 12; two h agree to 1e-5), not toward a constant; ṙ∞ converges for B ≤ 0.07 (0.00173, 0.01144, 0.04487 at B = 0.02, 0.04, 0.07;
  ∝ B^2.7, B^2.4), the local ṙ(8) ∝ B^2.5 (0.04 to 0.1) and B^1.1 (0.1 to 0.2); p*ₕ ∝ B^2.27 (extrapolated). The |B|³ law and the E_net^(3/2) scaling of D = 5 do not carry (ṙ∞/E_net^(3/2) = 0.136, 0.112, 0.083).
- Recommended task 4 is done; task 2's D = 6 half is open (below).

**Paper V** (`paper/5-ellis-threshold/main.tex`)
- l. 210–211: "G in D = 6 was run at the first two only" → G in D = 6 is run at all three.
- Table `tab:pstar`: the last row becomes `6 & G & $1.2658\times10^{-2}$ & $1.1088\times10^{-2}$ & $1.0697\times10^{-2}$ & 2.005 & $1.0567\times10^{-2}$`; the caption's last sentence ("G in $D=6$ has two spacings, whose ratio is $1.1416$") goes; the text after it: G converges at second order to $1.0567\times10^{-2}$ in $D=6$ as well.
- Table `tab:life`: caption "(G in $D=6$: $h=0.0025$)" → $h=0.00125$; the row `6 & 0.461924 & G, line & $0.996251\pm0.027127$ & $[0.798695, 1.209563]$`; the resolution sentence for G: $1.5\times10^{-3}$ between the two finest $h$.
- l. 586: "G in $D=6$ was not read on the locus" → G in $D=6$ was read on the locus (\#228): the levels that give the lifetime at $D=4$ do not at $D=6$ (the collapse side never reaches $r_{\min}=0.5$ for $e\ge4.5$; both sides cross $2$ at $\tau=15.75$), the pair's growth rate follows
  $r_{\rm th}(t)/W$ within $1.5\%$ while the radius grows by $26\%$, and the slope $0.996$ is the average of a local slope that rises from $0.71$ to $1.05$.

## Recommended tasks

1. **A departure-time lifetime on the locus.** What: define G's lifetime on the locus as the τ at which the run's deviation from the tuned throat reaches ε of r_th (from a pair at two δ), at D = 4, 5, 6, and fit its slope against r_th(τ_dep)/W.
   Why: #225's levels fail at D = 6 (FA, ED) and the diagonal slope is read by an observer that drifts off. Acceptance (can fail): at D = 6 dτ_dep/dln(1/δ) = r_th(τ_dep)/W to 2 %; the slope at D = 4 reproduces #225's 0.8495; a planted drifting throat recovers its rate; ε = 0.03, 0.1, 0.3 agree on the slope to 2 %. Depends on: this note.
2. **The D = 6 half of #273 (why ṙ∞ has the exponent it has).** What: B = 0.01, 0.014, 0.03 and 0.15 at h = 0.0025 with p*ₕ to 1e-12, local ṙ at t = 8, 10 and the fit where it converges. Why: the exponent falls from 2.7 to 1.1 between B = 0.02 and 0.15 at D = 6
   (3.16 at D = 5) and ṙ∞/E_net^(3/2) is not constant. Acceptance: the local exponent against B with its window and resolution parts; the model of #272's resolvent against ṙ∞ at B = 0.02. Depends on: this note, #272, #273.
3. **B = 0.2 at D = 6 with a longer triangle.** What: p*ₕ to 1e-12 at h = 0.0025 with a 3 h limit and the triangle of the pair run to t = 70. Why: the + side event lies beyond t = 34 and the bisection hit its 2 h limit. Acceptance: JOB6 met at B = 0.2; a pair at two h. Needs a waiver of the 2 h stop. Depends on: this note.
4. **`threshold.py`'s `not_run` entry.** What: delete `not_run` from `CRIT` (one line) so `threshold.py analyze` takes the third h. Why: the driver clears it today. Acceptance: `threshold.py analyze --raw=…/runs40` reproduces `threshold_output.txt`. Depends on: this note.

## Files

`g_d6.py` (the driver: `tail`, `threshold-analyze`, `locus`, `drift`, `analyze`, `mutate`, `stopped-pstar`), `g_d6_output.txt`, `threshold_output.txt` (regenerated), this note; raw `negative-mass-raw/I228-dh610` (README with the criteria, settings and mutations,
`runs/` the job, tail, series, stalls, pairs and bisections, `runs/analysis/`, `runs/mut/`, `runs40/` the merged record, `dev/`, `scripts/`).

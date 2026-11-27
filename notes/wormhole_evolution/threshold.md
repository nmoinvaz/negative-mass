# The lifetime law at the Ellis threshold, D = 4, 5, 6 (#40, #150)

Bisection to the discrete threshold p*ₕ in three pulse families, the lifetime T against ln(1/|p − p*ₕ|), and its slope against 1/W
(0.846045, 0.578803, 0.461924 in a, from #195's W²a² = 1.39705243433350684, 2.98495893035176025, 4.68662054311946479). The evolution is #39's
double-null code ([README.md](README.md), `dn.py`, `checks.py`, unchanged). Driver: [threshold.py](threshold.py). Output of record:
[threshold_output.txt](threshold_output.txt) (`python3 threshold.py analyze --raw=<negative-mass-raw>/I40-dh610/runs`). Raw output, with the criteria
written before the first run of record: negative-mass-raw `I40-dh610`. Status words as in AGENTS.md: every slope and p*ₕ below is *numerical evidence*.
The prediction tested is Paper III's (notes/wormhole_largeD/README.md, "A prediction for threshold evolutions"); [SH02] and [ST17] made the earlier
runs (critical_ghost.md §1.3). No sentence here says a result is new.

## Answer

In the two families whose threshold data are the static wormhole (S, symmetric two-sided pulse; O, one-sided ghost pulse) the lifetime slope equals 1/W
to 6e-4 relative at D = 4, 4e-4 at D = 5 and 3e-3 to 5e-3 at D = 6 (h ≤ 0.0025), with 1/W inside the fit spread, over 8 decades of |p − p*ₕ| from 1e-11 to 1e-3.
p*ₕ converges at second order. The family that does not start at the static wormhole (G, a one-sided ghost pulse plus a fixed normal pulse, added by the
coordinator's ruling) has slopes 6 % (D = 4), 38 % (D = 5) and 115 % (D = 6) above 1/W, with the first trapped cell off the diagonal and a drifting throat.
That is a difference from 1/W, reported with B untuned; this note does not say what causes it.

## What was measured

- **Families** (a = 1, 8πG = 1, Gaussian pulses in φ, centre 3, width 0.5, on the ray x⁻ = 0): S ghost amplitude p mirrored onto both rays; O ghost amplitude p
  on one ray; G ghost amplitude p plus a normal field of fixed amplitude B = 0.1 on the same ray. O's threshold is the unperturbed wormhole to O(h²)
  (p*ₕ(O) = 2 p*ₕ(S) to 2e-3 at the coarsest h and 1e-4 at the finest, all D), so O is not a generic family; the brief's description of it as the generic case
  does not hold, and G is what supplies a threshold away from the static wormhole (fate scans at D = 4: dev_02, dev_03 in the raw folder).
- **Events and T.** Collapse: first row of the x⁻ sweep with a cell where r₊ and r₋ are both below −1e-2. Expansion: first row whose diagonal cell has both above
  1e-2 and r ≥ 2. T = t = (x⁺ + x⁻)/√2 at the event cell's centre, start at the corner S (t = 0). In x⁻ the slope is 1/(√2 W).
- **Bisection** to hi − lo ≤ 1e-14 (the brief asks 1e-10), then T at p = p*ₕ ∓ 10^−e, e = 3, 3.5, …, 12, both sides. Above e = 11 the series is not
  a smooth function of p (roundoff, below); those points are run and not fitted.
- **Fit.** T = c_side + s ln(1/|p − p*ₕ|), common slope, one constant per side. Primary window [1e-11, 1e-3] (8 decades); four more windows;
  each window joint, collapse side alone, expansion side alone: 15 fits.
- **Error bars, three different things.** *OLS SE* of the primary fit: statistical, assumes independent residuals (event-row quantization √2h/√12 and
  roundoff). *Fit spread*: the envelope of s ± 2 SE over the 15 fits, so window and side choice plus noise; "agrees within its fit error" is 1/W inside it.
  *Resolution difference*: the relative difference of the primary slopes at the two finest h, bar 0.5 %.
- Resolutions: D = 4 h = 0.02, 0.01, 0.005; D = 5 h = 0.01/√2, 0.005/√2, 0.0025/√2; D = 6 h = 0.005, 0.0025, 0.00125. G at D = 6 was run at 0.005 and 0.0025 only
  (the third h is about 4 h, over the 2 h stop). Per the coordinator's ruling, D = 6 conclusions are drawn from h ≤ 0.0025; h = 0.005 is the coarse member.

## p*ₕ and its order

| D | family | p*ₕ at h₁, h₂, h₃ | order | Richardson p* |
|---|---|---|---|---|
| 4 | S | 4.7528e-05, 1.1884e-05, 2.9712e-06 | 2.000 (and 2.000, 2.000 from the ratios, p* = 0) | −7.4e-10 |
| 4 | O | 9.5052e-05, 2.3768e-05, 5.9423e-06 | 2.000 | −1.7e-09 |
| 4 | G | 1.49458e-02, 1.48774e-02, 1.48603e-02 | 2.000 | 1.48546e-02 |
| 5 | S | 2.8848e-04, 7.2091e-05, 1.8021e-05 | 2.001 (2.001, 2.000) | 9.6e-09 |
| 5 | O | 5.7699e-04, 1.4418e-04, 3.6042e-05 | 2.001 | 2.1e-08 |
| 5 | G | 1.12330e-02, 1.08191e-02, 1.07157e-02 | 2.001 | 1.06812e-02 |
| 6 | S | 1.1155e-03, 2.7784e-04, 6.9396e-05 | 2.007 (2.005, 2.001) | 3.5e-07 |
| 6 | O | 2.2329e-03, 5.5579e-04, 1.3880e-04 | 2.008 | 8.1e-07 |
| 6 | G | 1.26577e-02, 1.10880e-02 (two h) | not computed; ratio 1.1416 | — |

The continuum p* of S and O is 0 (the unperturbed wormhole), and the Richardson values are 0 to the extrapolation error (≤ 8e-7). G's is 0.01485 (D = 4) and
0.01068 (D = 5). p*ₕ of S falls by 3.999 to 4.015 per halving in every D (G1 of #39 found 3.999 and 4.002 at D = 4, 5).

## The slope against 1/W

Primary joint slope ± OLS SE, collapse side alone and expansion side alone (same window), fit spread, relative deviation from 1/W, and whether 1/W is inside the spread.

| D | fam | h | slope ± OLS SE | collapse | expansion | fit spread | rel dev | inside |
|---|---|---|---|---|---|---|---|---|
| 4 | S | 0.02 | 0.846642 ± 0.000380 | 0.847335 | 0.845950 | [0.844435, 0.849102] | +7.05e-4 | yes |
| 4 | S | 0.01 | 0.846522 ± 0.000323 | 0.846853 | 0.846191 | [0.845217, 0.848691] | +5.63e-4 | yes |
| 4 | S | 0.005 | 0.846522 ± 0.000261 | 0.846853 | 0.846191 | [0.845723, 0.848441] | +5.63e-4 | yes |
| 4 | O | 0.02 | 0.846040 ± 0.000302 | 0.845890 | 0.846191 | [0.844393, 0.847632] | −6.3e-6 | yes |
| 4 | O | 0.01 | 0.845859 ± 0.000158 | 0.845829 | 0.845890 | [0.844957, 0.846741] | −2.20e-4 | yes |
| 4 | O | 0.005 | 0.845867 ± 0.000168 | 0.845423 | 0.846311 | [0.844639, 0.846716] | −2.11e-4 | yes |
| 4 | G | 0.02 | 0.898291 ± 0.005730 | 0.894814 | 0.901768 | [0.864041, 0.958733] | +6.18e-2 | **no** |
| 4 | G | 0.01 | 0.898043 ± 0.005462 | 0.895250 | 0.900835 | [0.863276, 0.957708] | +6.15e-2 | **no** |
| 4 | G | 0.005 | 0.899089 ± 0.005454 | 0.896620 | 0.901558 | [0.864070, 0.958277] | +6.27e-2 | **no** |
| 5 | S | 0.00707 | 0.578346 ± 0.000127 | 0.578080 | 0.578612 | [0.577478, 0.579151] | −7.89e-4 | yes |
| 5 | S | 0.00354 | 0.579011 ± 0.000220 | 0.579283 | 0.578740 | [0.577925, 0.580133] | +3.60e-4 | yes |
| 5 | S | 0.00177 | 0.579027 ± 0.000229 | 0.579368 | 0.578687 | [0.578002, 0.580219] | +3.88e-4 | yes |
| 5 | O | 0.00707 | 0.579049 ± 0.000236 | 0.579251 | 0.578846 | [0.578164, 0.581068] | +4.24e-4 | yes |
| 5 | O | 0.00354 | 0.578857 ± 0.000128 | 0.578740 | 0.578974 | [0.577989, 0.579507] | +9.3e-5 | yes |
| 5 | O | 0.00177 | 0.578966 ± 0.000285 | 0.578155 | 0.579778 | [0.577095, 0.581630] | +2.82e-4 | yes |
| 5 | G | 0.00707 | 0.798948 ± 0.015509 | 0.723239 | 0.874656 | [0.707438, 0.935869] | +3.80e-1 | **no** |
| 5 | G | 0.00354 | 0.800044 ± 0.015453 | 0.724761 | 0.875327 | [0.708978, 0.936973] | +3.82e-1 | **no** |
| 5 | G | 0.00177 | 0.800277 ± 0.015565 | 0.724332 | 0.876221 | [0.709379, 0.938761] | +3.83e-1 | **no** |
| 6 | S | 0.005 | 0.457989 ± 0.000776 | 0.456439 | 0.459540 | [0.449488, 0.460292] | −8.52e-3 | **no** |
| 6 | S | 0.0025 | 0.460247 ± 0.000666 | 0.458336 | 0.462159 | [0.454915, 0.463275] | −3.63e-3 | yes |
| 6 | S | 0.00125 | 0.459628 ± 0.000868 | 0.458904 | 0.460353 | [0.455143, 0.462600] | −4.97e-3 | yes |
| 6 | O | 0.005 | 0.458629 ± 0.001063 | 0.457959 | 0.459299 | [0.453675, 0.468485] | −7.13e-3 | yes |
| 6 | O | 0.0025 | 0.460609 ± 0.000497 | 0.460255 | 0.460963 | [0.457147, 0.464379] | −2.85e-3 | yes |
| 6 | O | 0.00125 | 0.460325 ± 0.000791 | 0.461599 | 0.459051 | [0.452896, 0.465270] | −3.46e-3 | yes |
| 6 | G | 0.005 | 0.990046 ± 0.026391 | 0.854538 | 1.125554 | [0.795196, 1.190386] | +1.14 | **no** |
| 6 | G | 0.0025 | 0.994756 ± 0.026714 | 0.857876 | 1.131636 | [0.798040, 1.200660] | +1.15 | **no** |

- **Resolution agreement** (two finest h run, bar 0.5 %): S 1.5e-12, 2.8e-5, 1.35e-3 at D = 4, 5, 6; O 8.9e-6, 1.9e-4, 6.2e-4; G 1.2e-3, 2.9e-4, 4.7e-3. All pass.
  The D = 4 S figure of 1.5e-12 is a lattice effect, not a convergence statement: event rows agree between h to within one coarse row, so the slope can take
  only values on a lattice of spacing about 3e-5, and two coincide (raw dev_11: perturbing T by 0.003 moves the slope in the fifth digit).
- **D = 6 does not tighten toward 1/W with h.** S: −8.5e-3 (h = 0.005), −3.6e-3, −5.0e-3; O: −7.1e-3, −2.9e-3, −3.5e-3. The primary slope is 2.0σ to 2.7σ (OLS) below
  1/W at h = 0.0025 and 0.00125, inside the fit spread (±0.8 %). Two finest h agree with each other, so the offset is not a resolution effect at these h. The collapse-side
  slope is 0.4583, 0.4589 (S) and the expansion-side 0.4622, 0.4604; development runs (dev_04, dev_05) showed the collapse side is not linear above |p − p*ₕ| = 1e-5 at D = 6,
  and the primary window reaches 1e-3. That is a candidate cause, not tested here (**conjecture**). The coarse member h = 0.005 of S misses 1/W
  (−0.85 %, outside its spread [0.4495, 0.4603]); its p*ₕ does not break the three-h convergence (order 2.007).
- **At D = 4 and 5 the finest-h slope is within 2σ of 1/W** (OLS): S +1.8σ and O −1.1σ at D = 4, S +1.0σ and O +0.6σ at D = 5.

## #40 and #150: acceptance

| criterion | verdict | numbers |
|---|---|---|
| slope agrees with 1/W within its fit error at each D (#40) | S and O: met at all h, except S at D = 6, h = 0.005; G: not met, reported | table above |
| two resolutions agree on the slope (#40) | met, all families, D | resolution agreement above |
| threshold mass in units of a, with its dependence on the family (#40) | reported as a formation value, not a settled mass | below |
| planted type I series recovers its slope within the fit spread (#150) | met | 40 of 40 seeds at D = 4, 5, 6; p* recovered to 4.3e-16 in 40 of 40 |
| √2 time-unit control makes the comparison with 1/W fail (#150) | met | planted slope 1/(√2 W) rejected in 40 of 40 seeds; `--mutate=time_unit` fires PL, PX, SL, G6 |
| ≥ 8 decades of ln|p − p*ₕ|, bisection to ≤ 1e-10 (#150) | met | window [1e-11, 1e-3]; bisection to 1e-14, valid to 1e-11 |
| a note states what is measured (numerical evidence) and what is assumed (#150) | this note | assumptions: the diagonal event, the window, 1/W from #195 |
| D = 6 growth control before any D = 6 bisection | met | W_fit 2.160530, 2.162551 against 2.164860 (−2.0e-3, −1.07e-3), h = 0.00289, 0.00144, X = 12 |

#150 asked for D = 4 only; the same measurement was made at D = 5 and 6.

## Threshold mass (formation values, not settled masses; #41 stays open)

M at the first trapped cell for p*ₕ − 1e-8, from the cell-centred Misner–Sharp mass r^{D−3}(1 + 2e^f r₊r₋)/2 (the code's normalization, 8πG = 1, a = 1; M/a^{D−3}),
and r at that cell. h-dependence: three h.

| D | S (h₁, h₂, h₃) | O | G |
|---|---|---|---|
| 4 | 0.49420, 0.49414, 0.49416 (r 0.9882, 0.9881, 0.9881) | 0.49415, 0.49418, 0.49416 | 0.49515, 0.49519, 0.49544 (r 0.9901, 0.9901, 0.9907) |
| 5 | 0.49163, 0.49194, 0.49201 (r 0.9915, 0.9918, 0.9919) | 0.49121, 0.49181, 0.49199 | 0.68602, 0.68680, 0.68719 (r 1.171, 1.172, 1.172) |
| 6 | 0.48254, 0.48853, 0.48999 (r 0.9882, 0.9922, 0.9932) | 0.47877, 0.48789, 0.48986 | 2.57667, 2.60758 (two h; r 1.727, 1.734) |

In the first trapped cell S and O have r ≈ 0.99 a, the throat radius of the static wormhole, and
M ≈ 0.49 a^{D−3}; the two agree to 2e-4 at the finest h, and G differs by 1.3e-3 at D = 4 and is larger at D = 5, 6 (M 0.687 and 2.6 with r 1.17 and 1.73). The h-dependence at D = 6 is second order
(differences 0.00599, 0.00146: ratio 4.1). These are the mass where the trapped test first holds, in a gauge whose grid stretches after a horizon forms; they are not
[SH02]'s 0.30 a or [GGS09b]'s 0.22 a, whose definitions and units are #41's question, and not a mass after the horizon has settled (#39's recommended task 3).

## What the near-critical run looks like (reported, bears on #42)

At p*ₕ ∓ 1e-8, the minimum of r on the slices t = const during the stall, and the event cell relative to the diagonal.

| D | family | r at the throat locus (stall) | position l = (x⁺ − x⁻)/√2, drift dl/dt | event cell offset from the diagonal, l (collapse side) |
|---|---|---|---|---|
| 4, 5, 6 | S, O | 1.0000 at D = 4; 0.9996 to 1.0000 at D = 5; 0.9989 to 0.9998 at D = 6, h ≤ 0.0025 (0.9951 to 0.9957 at h = 0.005) | 0, 0 | 0 (O D = 6 h = 0.005: −0.0035) |
| 4 | G | 1.0028 (collapse), 1.0033 (expansion) | +0.584, +0.0735 | +0.863 |
| 5 | G | 1.0618 | +0.254, +0.127 | +1.167 |
| 6 | G | 1.65 (collapse), 1.20 (expansion); stall span under 1 time unit, no plateau | +1.31 and +0.43, +0.11 and +0.20 | +1.51 |

For S and O the stall is the static wormhole of radius a, on the diagonal, for 9 to 16 time units. For G the throat locus sits off the diagonal and moves; at D = 6 there is
no plateau by the 2 % criterion. The expansion event is read on the diagonal and the collapse event is the first trapped cell anywhere in the row, so the lifetime of G
is that of an observer that is not at the throat. That bears on the slope (it adds the drift to T) and is the first thing to separate before reading G's departure from 1/W
as a property of the critical solution; it is not separated here. The same G slope at the three h (+6.2 %, +6.2 %, +6.3 % at D = 4) says it is not a resolution effect.
ED for G fails at D = 4: with the alternative events (trapped tolerance 0.1, and rexp = 3) the slopes are 0.9193 ± 0.0145 and 0.8614 ± 0.0075, 6 % apart, and the second lies
just below the primary spread [0.8633, 0.9577]; for S and O the alternatives lie inside the spread at every D, and for G at D = 5 and 6.

**Update from #225** ([throat_locus.md](throat_locus.md)):
- In the stall table above, `l` and `dl/dt` are the locus's position and its drift. The drift speeds are +0.0735 and
  +0.127 at D = 4 and 5, not the +0.58 and +0.25 of the position column.
- In the locus's proper time, G's slope at D = 4 is 0.8495, +0.4 % over 1/W, with the two finest h agreeing to 5.5e-5.
  That is the e-fold of the stall's static member: γ₁ = 0.048, r_th = 1.0028, T(γ₁) r_th = 0.847.
- So the +6 % above is the diagonal observer: a planted drift of 0.1 reads +7.8 %.
- At D = 5 the member's radius drifts (1.07 to 1.20). The rate follows W/r_th(t), and the slope is not a single
  e-fold.

## Constraint diagnostic (reported; no bar)

Max |C₊₊|, |C₋₋| over the nodes of the event cell's past, along the two runs that bracket p*ₕ (last collapse, last expansion), and the same divided by the largest sum of the
absolute values of the terms ("rel"); `repro` is the difference between `dn.constraints` and the sum assembled from its terms (≤ 2.2e-16). `checks.check_C` takes the unscaled
maximum; the coordinator's message says "scaled as C1 does", which it does not, so both are in the output.

| D | family | max \|C\| at h₁, h₂, h₃ | ratio h₁→h₂, h₂→h₃ (second order: 4) |
|---|---|---|---|
| 4 | S | 5.92e-05, 1.48e-05, 3.76e-06 | 3.998, 3.941 |
| 4 | O | 5.92e-05, 1.48e-05, 3.79e-06 | 3.999, 3.905 |
| 4 | G | 1.05e-04, 2.61e-05, 6.54e-06 | 3.997, 4.000 |
| 5 | S | 4.69e-05, 1.23e-05, 8.64e-06 | 3.807, **1.425** |
| 5 | O | 4.70e-05, 1.26e-05, 7.83e-06 | 3.731, **1.607** |
| 5 | G | 4.70e-05, 1.18e-05, 7.67e-06 | 3.995, **1.533** |
| 6 | S | 7.78e-05, 2.29e-05, 5.28e-05 | 3.401, **0.433** |
| 6 | O | 7.99e-05, 3.25e-05, 5.10e-05 | 2.459, **0.637** |
| 6 | G | 7.37e-05, 1.83e-05 (two h) | 4.032 |

The violation is truncation error where it converges: ratio 4 per halving at D = 4 in every family, and at the first halving at D = 5 and 6 (2.5 to 4.0, G at D = 6 4.03).
**At the finest h at D = 5 and 6 it stops falling** (ratios 1.4 to 1.6 at D = 5, 0.4 to 0.6 at D = 6, i.e. larger than at h/2). Over the same halving p*ₕ keeps its second-order
ratio (4.00) and the slope changes by 1.6e-5 to 2.3e-4 (D = 5) and −6.2e-4, −2.8e-4 (D = 6, S and O). So the constraint maximum at the finest h is not the quantity that
sets p*ₕ or the slope.

**Update from #226** ([constraint_plateau.md](constraint_plateau.md)): at the finest h the maximum is roundoff (*numerical evidence*). The second
difference of r carries the roundoff 2⁻⁵² r of each node as ε/h², and C rⁿ⁻¹ is constant along x⁻, so a violation made at large r is amplified
by (r₀/r)ⁿ⁻¹ at the throat. In binary128 the same slabs fall by 4.0 to 4.2 per halving at D = 5 and 6. Noise of a·2⁻⁵² per node raises the
finest-h maximum in proportion to a. The maximum is not a single cell and not the corner, and the ray data's RK4 step does not change it. The
roundoff in r is 1e-11 against a truncation error of 4e-4, so the plateau does not bear on p*ₕ or the slope. The finest-h entries of the table
at D = 5 and 6 are roundoff, not truncation. The coordinator's own march reproduces the static ratios (4.197, 3.833, 1.356 against 4.204, 3.758,
1.337) and the noise model to 1 to 7 % (`I226-mac-review`).

## Roundoff floor

T(p) is a log law for |p − p*ₕ| ≥ 1e-11 at D = 4, 5, 6 (1e-12 at D = 4 in the first development run). Below it T jitters by about 0.1 and, at e ≥ 12 to 13, points land
on the wrong side of p*ₕ (raw dev_05, bracket 5.7e-16, series to e = 14.5). The series of record stops at e = 12: the 4 floor points per job (e = 11.5, 12) all landed on their own side in
the 26 jobs, and the fate-monotone check (FM) covers e ≤ 11. The bisection to 1e-14 is therefore below what
p*ₕ means: p*ₕ is determined to about 1e-12 by the roundoff of the evolution. Compare with the brief's bar of 1e-10.

## Controls and mutations

- **Planted series** through the same `bisect`, `series`, `fit_set`: noise σ = 0.01 plus the event-row quantization, p* = 4.7e-5, T₀ = (−0.79, +4.05). 40 of 40 seeds recover p* to 4.3e-16 and have 1/W inside the spread
  at D = 4, 5, 6 (mean primary relative error ≈ −5e-5 to −1e-4, rms 4e-4 to 6e-4); a series planted with slope 1/(√2 W) is rejected in 40 of 40.
- **Mutations** (each target computed in development before it was declared; the runs of record repeat them at D = 4, S, three h; see the end of threshold_output.txt and `runs/mut` in the raw folder).

  | mutation | what it breaks | fired in the run of record (equal to the declared target) |
  |---|---|---|
  | `time_unit` | T scaled by √2 in the fits and in the growth-rate fit | PL, PX, SL on the stored runs; G6 in the mutated growth control (W_fit 1.5277, 1.5292: −29.4 % = 1/√2 − 1) |
  | `pstar_off` | p*ₕ shifted by 1e-9 in the fits | FM, PL, PX, RA, SL |
  | `bisect_loose` | bisection stopped at 1e-9 | BR, FM, RA |
  | `ed_unit` | the alternative events timed in x⁻ units | ED |
  | `pstar_h1` | p*ₕ replaced by p*ₕ h₁/h (first-order seed) | PS |
  | `no_corrector` | #39: predictor-only update | SL (at h = 0.02, 0.01) |
  | `theta_coef`, `ghost_sign`, `finder_sign`, `ellis_charge`, `scalar_coef`, `data_constraint_sign` | #39's input mutations | JOB (the bracket does not straddle) |

  The logs are `runs/mut/analysis_<name>.log` (`time_unit` in two: `_runs` and `_g6`); each ends with the expected and the fired set. The constraint, formation and stall tables are reports and have no mutation.
- **One check that did not fire in development was not declared**: `ed_event` (alternative expansion event at r ≥ 0.5) leaves the slope unchanged by design.

## What Paper III's threshold section and critical_ghost.md §1.3 should now say (proposals for the coordinator)

Paper III, §"Consequence for the collapse–expansion threshold" (Table 3 and the last paragraph of the section):

- Add a row to Table 3 for this work, in t and in x⁻ (÷√2): D = 4, 5, 6 predicted 0.846, 0.579, 0.462 (0.598, 0.409, 0.327 in x⁻), measured S 0.8465, 0.5790, 0.4596 (finest h), with the fit spread of about ±0.2 % (D = 4, 5) and ±0.8 % (D = 6), over 8 decades of |p − p*ₕ|. Say these are *numerical evidence* from a double-null evolution
  with p* determined by bisection, for a symmetric and a one-sided ghost pulse whose threshold data are the static wormhole.
- Replace "A bisection to the threshold in D = 5 and 6 … would test Eq. (lifetime) properly" by what was done, and state its limit: a family that passes through the static wormhole tests the linear
  growth rate in a nonlinear code, not that generic near-threshold data approach the wormhole with that e-fold.
- Add the one-sided generic family result: with a fixed normal pulse added (G), the slope is 6 %, 38 % and 115 % above 1/W at D = 4, 5, 6, with the first trapped cell displaced from the diagonal and a throat that
  drifts; the cause is not established here. Do not describe the generic case as agreeing.
- The "[ST17] agree within their uncertainties" sentence stays as the comparison with the two-run slopes.

critical_ghost.md §1.3 ("Has the boundary been treated as a critical phenomenon?"):

- "Nobody did it with the standard protocol" can become: a bisection to p* with the lifetime against ln|p − p*ₕ| over 8 decades, in a symmetric and a one-sided ghost family, gives the
  slope 1/W of the linear mode to 6e-4 (D = 4), 4e-4 (D = 5) and 5e-3 (D = 6) relative (this note); a generic one-sided family with a fixed normal pulse does not (+6 % to +115 %). Keep [SH02]'s two-point
  slopes 0.626 and 0.585 against 0.598 as they are.
- "Is the threshold black-hole mass finite?": the mass where the trapped test first holds at p*ₕ − 1e-8 is 0.494, 0.492, 0.490 a^{D−3} (D = 4, 5, 6; code's normalization), the same in S and O,
  with r ≈ 0.99 a at the cell; it is a formation value in a grid that stretches, not [SH02]'s 0.30 a or [GGS09b]'s 0.22 a, and does not settle #41.
- "A difference from textbook type I": the stall of S and O is the static wormhole of radius a to 1e-4 (r = 1.0000); G's stall is displaced and moving.

## Limits

- The diagonal observer. G's lifetime is read on the diagonal (expansion) or at the first trapped cell anywhere in the row (collapse); a lifetime along the throat locus is not measured.
- G at D = 6 has two h; the third would be 4 h. G at D = 5, h = 0.00707 (BR fails): X = 22 was too short for the coarsest h, the bisection ended at 2.9e-12 with "static through X"; the series
  is intact (FM passes) and the dev run with X = 30 reached 5.7e-15.
- 1/W is #195's; the planted check shares the fit code with the real runs, by design, and shares none of the evolution.
- The event definition, the window and the roundoff floor are the assumptions in the answer; the p-independence of the slope under the alternative events (tolerance 0.1, rexp = 3) holds for S and O at every D
  (ED), and fails for G at D = 4.
- Not measured: the settled black-hole mass (#41, #39's task 3); the asymmetric static members that the one-sided families may approach (#42); the second collapse window of [GGS09b].

## Recommended tasks

1. *Lifetime along the throat locus for one-sided families (#42, #40).* What: measure G's lifetime as proper time along the minimum-r locus, and the drift velocity, at D = 4 and 5. Why: G's slope
   is read by an observer that is off the throat. Acceptance: the slope in the locus's proper time agrees between two h; a planted boosted static case recovers its known boost; the diagonal slope is
   reproduced as a control. Depends on: this note.
2. *The constraint plateau at the finest h (#39).* Done by #226. What: find why max |C| stops falling at h = 0.00177 (D = 5) and 0.00125 (D = 6). Acceptance: a run at h and h/2 on a smooth static case shows ratio 4 or the floor
   is attributed to a named term; a mutation that doubles the roundoff changes it. Depends on: #39.
3. *The roundoff floor (#40).* What: bisect with extended precision (or perturb the data at 1e-16) to see whether the 1e-11 floor is the evolution's. Acceptance: the series is smooth to the declared level in the
   extended run. Depends on: none.
4. *G at D = 6 at the third h.* Cost about 4 h in one unit; needs a waiver of the 2 h stop. Acceptance: the three-h order of p*ₕ at D = 6 for G.

## Files

`threshold.py` (this driver), `threshold_output.txt` (its analysis of the runs of record), this note; raw `negative-mass-raw/I40-dh610` (README with the criteria, `dev/`, `runs/`, `scripts/`).
`python3 threshold.py job|g6|planted|analyze` as in its docstring, from this directory with numba available; `--kernel` is not an option: the production kernel is numba, as in #39.

## Coordinator's review

Review of branch `wip/i40-ellis-threshold` at 2e97643 (code 13a1301) and raw `I40-dh610` at 99423f6, recorded in negative-mass-raw `I40-mac-review`
(criteria committed at e7e61a8 before the first run; runs at 2376dd2). On the Mac (Apple M5), Python 3.14.7, numba 0.68.1, nice 10, two processes.

- **Scope.** The branch adds this note, `threshold.py` and `threshold_output.txt`, nothing else. `dn.py` and `checks.py` are unchanged.
- **Rerun.**
  - The analysis of record over the runner's 26 jobs prints every line of `threshold_output.txt` again. The only differences are one relative
    difference of two equal slopes (1.52e-12 against 1.51e-12) and the timing. The tally is the same: 100 of 111 checks, the same 11 failures.
  - 13 of the 26 jobs were rerun: D = 4 S, O and G at three h; D = 5 S at the two coarse h; D = 6 S at the two coarse h.
  - In each job p*ₕ equals the runner's within 2.9e-12, and every series point with e ≤ 9 has the runner's fate and event cell.
  - Points with e ≥ 9.5 to 10 differ between the two machines in 10 of the 13 jobs. The rerun's slopes over [1e-11, 1e-3] move by up to 1.4e-4 through
    them (D = 4 S, h = 0.005).
  - The review's first command rounded h at D = 5 to six figures. That moved p*ₕ by 1.8e-10, which is p*ₕ ∝ h² at the rounded h. The error is the
    review's; rerun 2 used the runner's h and agrees (amendment A1 there).
- **Independent check** (`threshold_check.py`, a copy of the review's `ind40.py`; output `threshold_check_output.txt`).
  - The coordinator wrote a second double-null code for family S. It integrates the background and the ray constraint itself. It uses its own diamond
    scheme on the half x⁺ ≥ x⁻, which the symmetry allows. Its events are r = 0.5 and r = 2 on the throat (the diagonal), interpolated in t.
  - It shares no code with `dn.py` or `threshold.py`. What it shares: the field equations (derived again in the review of #39), the family's
    definition, and W from #195.
  - 5 of 5 criteria pass:

    | D | slope against 1/W at h₁, h₂, h₃ ([1e-9, 1e-4]) | #40's S points at its finest h, refit over [1e-9, 1e-4] | p*ₕ order |
    |---|---|---|---|
    | 4 | −3.3e-5, −6.9e-6, −2.6e-6 | 0.846040 (−5.9e-6) | 2.000 |
    | 5 | −4.9e-4, −1.35e-4, −2.5e-5 | 0.578786 (−2.9e-5) | 2.001 |
    | 6 | −6.3e-3, −1.56e-3, −4.7e-4 | 0.461770 (−3.3e-4) | 2.007 |

  - The two codes agree to 4e-6 at D = 4 and 5 and to 1.4e-4 at D = 6, on the same window. Their p*ₕ agree to 2e-4 relative, so the two schemes'
    truncation seeds are close.
  - The slope in the throat's proper time τ equals the slope in t to 1e-5.
  - Mutations: `oracle-slope`, `fit-pstar` and `time-unit` fire exactly their targets. `time-unit` gives 0.598243 = 1/(√2W). `sphere-coef` stops at
    the bracket before any criterion is evaluated, as #40's six JOB mutations do.
- **What this changes in the answer.**
  - For S and O the slope is 1/W, converging at second order in h at D = 5 and 6 once the window stops at 1e-4. At the finest h it is within 3e-6,
    3e-5 and 5e-4 of 1/W at D = 4, 5, 6 (*numerical evidence*, two codes).
  - **Update from the review:** #40's D = 6 offset (−3.6e-3 and −5.0e-3, "does not tighten") comes from the decade 1e-4 to 1e-3 of its primary
    window. Its own points over [1e-9, 1e-4] give −3.3e-4. This confirms the runner's conjecture that the collapse side is not linear above 1e-5
    there.
  - #40's +5.6e-4 at D = 4 is the event-row quantization this note describes. With interpolated events and the same data window it is −5.9e-6.
  - G's departure from 1/W stands as reported, with its cause open.
- **Follow-ups.** The four recommended tasks are filed as #225 to #228 (see the issues). Task 1 carries the review's remark that the asymmetric
  static members have *shorter* e-folding times per throat radius (critical_ghost.md §1.1), so a different static member alone does not explain
  G's larger slope.

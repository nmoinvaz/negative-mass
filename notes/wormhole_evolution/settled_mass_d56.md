# The settled black-hole mass at D = 5 and 6, families S and G, and a planted Tangherlini (#246)

Code: [settled_mass_d56.py](settled_mass_d56.py), a thin driver that imports #41's `settled_mass.py`, #253's `settled_mass_g.py`, #40's `threshold.py` and #39's `dn.py` unchanged. Output of record:
[settled_mass_d56_output.txt](settled_mass_d56_output.txt) (`python3 settled_mass_d56.py analyze --raw=<negative-mass-raw>/I246-dh610`, exit 1: the four RC8a fail, as declared). Raw output, with the criteria
written before the runs of record (README at 726241f9, code at 0c52443; runs, analysis, mutated analyses and figure at 085c3667): negative-mass-raw `I246-dh610` (`README.md`, `runs/`, `mut/`, `dev/`, `posts/`, `figures/settled_mass_d56.svg`). Status words as in AGENTS.md: every mass
below is *numerical evidence*; what is not shown is marked *conjecture*. No sentence here says a result is new. Units: a = 1 is the throat radius of the initial static wormhole, 8πG = 1, D = n + 2, M = r^(n−1)/2 on the apparent
horizon (r₊ = 0), so M_throat = 1/2 and E/E_throat = 2M = r_AH^(n−1).

## Answer

**Family S settles at M_AH = 0.14587 a² (r_AH = 0.5401 a) at D = 5 and 0.10792 a³ (r_AH = 0.5996 a) at D = 6**, against 0.22308 a at D = 4 (#41): E/E_throat = 0.2917 and 0.2158 against 0.4462. These are second-order
Richardson values at three h, δ = 1e-8, to e^{−f} = 4.0e8 on the horizon (the finest-h plateaus are 0.145859 and 0.107793); the same for δ = 1e-6 to 1.4e-6 (D = 5) and 2.4e-5 (D = 6). **Family G does not settle at S's mass and
depends strongly on δ**: 0.26395 a² (D = 5) and 0.86286 a³ (D = 6) at δ = 1e-8, 0.23500 and 0.56793 at δ = 1e-6; G − S is +0.118 and +0.754, against +0.0125 at D = 4 (#253). [ST17]'s E_f is the static throat's own energy in
each D, and the mass on the horizon at its first trapped cell is 0.984 (D = 5) and 0.980 (D = 6) of it: their value is a formation value, as #41 found at D = 4.

| acceptance (issue #246) | verdict | numbers |
|---|---|---|
| planted Tangherlini recovered at second order | **met** (PT1, four cases) | max \|M − M_exact\| at the finest h: 3.2e-6, 3.3e-6 (D = 5, S-like, G-like), 3.3e-6, 5.8e-6 (D = 6); orders 1.88 to 2.08 over four h (band 1.75 to 2.25); bars 3e-4 M_exact |
| M_AH(x⁺) second order at three h up to a declared limit | **met** (RC6) | observed orders 1.92 to 2.13, three windows, both δ, all four cases, at e^{−f} = 4.0e8 and at #41's 1.57e6 (RC6a, same orders) |
| settled at the limit (RC8) | met at 4.0e8; **not met at 1.57e6** (RC8a, all four cases) | flux over the last 0.5 of x⁺ at 4.0e8: 5e-7 to 1.1e-6; at 1.57e6: 1.2e-5 (D = 5 G) to 7.2e-5 (D = 6 G), bar 1e-5; declared as expected (below) |
| settled values against [ST17]'s rows in one unit | **done** | E = c_D M, c_D = 12.566 (D = 5), 19.739 (D = 6); E/E_throat = 2M; table below; CONV passes |
| G against S per D, with the δ dependence | **done**; G differs | \|G − S\| = 0.118 > bar 0.029 (D = 5), 0.754 > 0.297 (D = 6); δ dependence +2.89e-2 and +2.95e-1, the same at all three h |

Controls: 69 counted checks, 65 pass, the 4 failing are RC8a of the four cases. PL1 to PL3 (#41's, through the driver) pass and are printed by #253's helper outside the count (see "Disagreements"). 26 of 26 mutations fire their
declared sets in the runs of record (`mut/`, exit 0 each). The roundoff floor of #226 is **not reached** by M_AH's convergence at D = 5 and 6 (ratios 0.003 to 0.12, below).

## What was run

- #40's p*ₕ of S and G at D = 5, 6 from its raw (the h of each job.log; hi − lo 9.1e-15 (S), 5.7e-15 (G), except D = 5, G, h = 0.00707: 2.9e-12). One input was not #40's: #40 has no run of G at D = 6, h = 0.00125, so p*ₕ = 0.01069693394962485
  was bisected here with `threshold.bisect` unchanged (hi − lo 5.7e-15, 44 runs, 2569 s; the three p*ₕ of G at D = 6 have differences 1.5699e-3 and 3.911e-4, ratio 4.014). #228 can reuse it (`runs/pstar_D6_G_k2.json`).
- #41's chain with two changes (raw README, "Changes made during development"): K = 28 (e^{−f} on the horizon to 4.03e8) and dxp = 14. Reason: the x⁺ span after formation shrinks with D (D = 4: 9.4; D = 5: 6.1 (S), 8.1 (G); D = 6: 4.6 (S), 8.5 (G)
  at #41's limit), so the last window at 1.57e6 still holds the decay tail. The K = 28 series contain the K = 20 and K = 24 series as bit-identical prefixes (PFX; M_end at K = 20, 24, 28: 0.145768, 0.145758, 0.145764 at D = 5 S;
  0.105828, 0.105841, 0.105833 at D = 6 S). Every limit check is made at both limits; the quoted values are the K = 28 ones (ruling of the coordinator).
- 24 runs (D = 5, 6; S, G; three h; δ = 1e-6, 1e-8), 8 prefix runs, 48 ulp-noise runs; 10 to 470 s each at one thread, median 80 s (the finest D = 6 runs 4 to 8 min and up to 9.2 GB; four runs show 3200 to 3700 s of wall time because they spanned a
  53-minute gap, about 20:16 to 21:09 local on 2026-11-27, in which dh610 did not answer; their results are unaffected). Level 0 of every chain equals `dn.evolve` bit for bit (L0ID).
  The first law on the horizon holds to 3.5e-5 of the mass change (RC4, general n: dM/dx⁺ = (r^n/n) e^f r₋ (φ₀₊² − φ₁₊²)). M_AH at the first trapped cell is within 9.4e-4 of #40's M_hor (RC3); #40's first trapped cells are reproduced (RC2).

## The planted Tangherlini (PT0, PT1)

The exact Tangherlini solution (r_h = λ, M = λ^(n−1)/2) in closed form for n = 3, 4 (the tortoise coordinate in elementary functions), in the gauge x⁻ = x_h + λU, V = κ e^{kp(x⁺−xs)}/(λ kp), in which the horizon is
x⁻ = x_h and e^{−f} on it is e^{kp(x⁺−xs)}: the stretching of the real runs. PT0: it agrees with `dn.kruskal_exact` (a different code) to 1.4e-13 in r and f at 400 points, and the horizon is where it should be to 0. PT1: the same chain (K = 28 levels) from the rays x⁻ = 0 and x⁺ = 0.
(λ, kp) of each case are the settled r_AH and the late stretching rate d ln e^{−f}/dx⁺ of the real S and G runs, rounded.

| case | M_exact | max \|M − M_exact\| at h = 2h₁, h₁, h₂, h₃ | orders | bar at h₃ |
|---|---|---|---|---|
| D = 5, S-like (λ 0.54, kp 2.7), h = 0.01414 … 0.00177 | 0.14580 | 2.11e-4, 5.33e-5, 1.27e-5, 3.19e-6 | 1.99, 2.07, 1.99 | 4.4e-5 |
| D = 5, G-like (λ 0.73, kp 1.9) | 0.26645 | 2.20e-4, 5.22e-5, 1.34e-5, 3.35e-6 | 2.08, 1.96, 2.00 | 8.0e-5 |
| D = 6, S-like (λ 0.60, kp 3.6), h = 0.01 … 0.00125 | 0.10800 | 1.99e-4, 4.76e-5, 1.30e-5, 3.26e-6 | 2.07, 1.88, 1.99 | 3.2e-5 |
| D = 6, G-like (λ 1.19, kp 1.75) | 0.84258 | 3.62e-4, 9.59e-5, 2.33e-5, 5.82e-6 | 1.92, 2.04, 2.00 | 2.5e-4 |

Every run reaches e^{−f} = 4.0e8 (status `limit`). The tracker's error at the finest real h (3e-6 to 6e-6) is below the h differences of the real plateaus (2e-5 to 2e-3, table "Error bars"), so those are the resolution of the
physical series (the rapid fall of M after formation), not of the horizon finder. The planted runs are exact data: the floor of #226 is not tested by them.

## M_AH(x⁺) after formation

Finest h (D = 5: 0.00177, D = 6: 0.00125), δ = 1e-8, from `settled_mass_d56_output.txt`; x⁺_f is the first trapped cell (D = 5: S 7.5157, G 9.9817; D = 6: S 6.2581, G 10.9206); the rows are matched in x⁺ − x⁺_f and not in t. Both δ and
the time t = (x⁺ + x⁻)/√2 on the horizon are in the output.

D = 5:

| x⁺ − x⁺_f | r_AH (S) | M_AH (S) | e^{−f} (S) | r_AH (G) | M_AH (G) | e^{−f} (G) |
|---|---|---|---|---|---|---|
| 0 | 0.99202 | 0.49205 | 1.022 | 1.17241 | 0.68728 | 1.006 |
| 0.5 | 0.97325 | 0.47361 | 1.079 | 1.15186 | 0.66340 | 1.108 |
| 1 | 0.91276 | 0.41657 | 1.306 | 1.09061 | 0.59472 | 1.392 |
| 1.5 | 0.75205 | 0.28279 | 2.653 | 0.96739 | 0.46793 | 2.379 |
| 2 | 0.57772 | 0.16688 | 15.32 | 0.82187 | 0.33774 | 6.639 |
| 2.5 | 0.54666 | 0.14942 | 82.63 | 0.74762 | 0.27947 | 22.91 |
| 3 | 0.54582 | 0.14896 | 344.7 | 0.73269 | 0.26842 | 69.21 |
| 4 | 0.54066 | 0.14616 | 5557 | 0.73091 | 0.26712 | 510.9 |
| 5 | 0.54014 | 0.14588 | 8.353e+04 | 0.72739 | 0.26454 | 3628 |
| 6 | 0.54012 | 0.14586 | 1.217e+06 | 0.72658 | 0.26396 | 2.526e+04 |
| 8 | 0.54011 | 0.14586 | 2.456e+08 | 0.72654 | 0.26393 | 1.169e+06 |

D = 6:

| x⁺ − x⁺_f | r_AH (S) | M_AH (S) | e^{−f} (S) | r_AH (G) | M_AH (G) | e^{−f} (G) |
|---|---|---|---|---|---|---|
| 0 | 0.99331 | 0.49003 | 1.022 | 1.73583 | 2.61513 | 1.011 |
| 0.5 | 0.97038 | 0.45688 | 1.109 | 1.70659 | 2.48516 | 1.184 |
| 1 | 0.86927 | 0.32842 | 1.712 | 1.61542 | 2.10777 | 1.647 |
| 1.5 | 0.64303 | 0.13294 | 13.6 | 1.44870 | 1.52023 | 3.287 |
| 2 | 0.61043 | 0.11373 | 122.5 | 1.28482 | 1.06046 | 10.13 |
| 2.5 | 0.60456 | 0.11048 | 785.8 | 1.22064 | 0.90935 | 31.68 |
| 3 | 0.60012 | 0.10806 | 5014 | 1.21225 | 0.89073 | 83.19 |
| 4 | 0.59971 | 0.10784 | 1.839e+05 | 1.20770 | 0.88073 | 497.5 |
| 5 | 0.59962 | 0.10779 | 6.507e+06 | 1.20047 | 0.86502 | 2964 |
| 6 | 0.59961 | 0.10779 | 2.271e+08 | 1.19945 | 0.86282 | 1.72e+04 |
| 8 | – | – | – | 1.19919 | 0.86225 | 5.592e+05 |

M_AH is monotone decreasing in each series (the ghost flux, r₋ < 0). The fall after formation is steeper than at D = 4 (S, between x⁺ − x⁺_f = 1 and 2: 0.471 to 0.378 at D = 4, 0.417 to 0.167 at D = 5, 0.328 to 0.114 at D = 6); M_AH is within 0.3 % of its settled value from
x⁺ − x⁺_f = 4 (D = 5, S), 3 (D = 6, S), 5 (D = 5, G) and 5.5 (D = 6, G). G lags S as at D = 4 (at D = 6 it reaches e^{−f} = 100 at x⁺ − x⁺_f ≈ 3.1 against 1.95) and forms a larger horizon: r_AH at the first trapped cell is 1.172 a (D = 5) and 1.736 a (D = 6), against 0.992 and 0.993
for S. The late stretching rates d ln e^{−f}/dx⁺ are 2.7 (S) and 1.9 (G) at D = 5, 3.6 and 1.8 at D = 6.

## The settled values and what their error bars are

Quantity: plateau mean of M_AH over the last 1.0 of x⁺ at e^{−f} = 4.0e8 (`plat1`), δ = 1e-8, finest h. The kinds are kept separate; none is a rigorous enclosure.

| case | value (finest h) | Richardson (second order, three h) | h: finest pair | h: spread over the three | δ: \|plat(1e-6) − plat(1e-8)\| | r_AH | E/E_throat |
|---|---|---|---|---|---|---|---|
| D = 5, S | 0.145859 | 0.145866 | 2.15e-5 | 1.09e-4 | 1.39e-6 | 0.54011 | 0.2917 |
| D = 5, G | 0.263924 | 0.263951 | 8.08e-5 | 4.17e-4 | 2.89e-2 (physical) | 0.72653 | 0.5278 |
| D = 6, S | 0.107793 | 0.107920 | 3.81e-4 | 1.96e-3 | 2.43e-5 | 0.59962 | 0.2156 |
| D = 6, G | 0.862180 | 0.862862 | 2.05e-3 | 1.10e-2 | 2.95e-1 (physical) | 1.19916 | 1.7244 |

- *h* terms are differences between resolutions of this code (the series is second order, RC6; Richardson is the better value, and at D = 6 it differs from the finest-h plateau by 1.3e-4 (S) and 6.8e-4 (G)).
- *δ* terms for G are a physical dependence on the distance from threshold and not an error; for S they are the (small) difference between two δ and are second order in h (D = 6: 3.85e-4, 9.33e-5, 2.43e-5 at the three h, observed orders
  2.05, 1.94; D = 5: 1.59e-5, 3.99e-6, 1.39e-6).
- *limit*: the plateau at #41's limit (1.57e6, the same series cut) is 4.5e-5 (D = 6, S) and 6.6e-5 (D = 6, G) higher than at 4.0e8 and 1.0e-5 (D = 5, S and G) higher, because its last window still holds the decay tail; the end values M_end differ by 4e-6
  or less between the limits (0.107796 and 0.107792, D = 6 S, h = 0.00125). RC8a fails there (plateau change 3.9e-4 (D = 6 S), 3.1e-4 (D = 6 G), 3.9e-5 (D = 5 S), 1.4e-5 (D = 5 G); flux 6.3e-5, 7.1e-5, 1.5e-5, 1.3e-5), and is reported and not loosened.
- *machine to machine*: #41's 2e-5 at D = 4 is not remeasured here.

## G against S per D, and the dependence on δ

| | D = 4 (#41, #253) | D = 5 | D = 6 |
|---|---|---|---|
| S settled M | 0.22308 | 0.145859 | 0.107793 |
| G settled M, δ = 1e-8 | 0.23555 | 0.263924 | 0.862180 |
| G settled M, δ = 1e-6 | 0.23498 | 0.234978 | 0.567448 |
| G − S (δ = 1e-8) | +0.0125 | +0.1181 | +0.7544 |
| G/S | 1.056 | 1.809 | 7.999 |
| declared bar (h + δ + 2e-5) | 6.0e-4 | 2.90e-2 | 2.97e-1 |
| verdict | differs | **differs** | **differs** |
| G, δ = 1e-8 against 1e-6 | +5.7e-4 (0.24 %) | +2.89e-2 (12 %), at all three h | +2.95e-1 (52 %), at all three h |
| r_AH at the first trapped cell, G (δ = 1e-8; δ = 1e-6) | 0.991 | 1.172 ; 1.098 | 1.736 ; 1.451 |

The difference from S grows steeply with D: 0.0125 a, 0.118 a², 0.754 a³. G's mass also depends on δ much more than at D = 4, and the same at the three h (DEP: 2.89e-2 at each h at D = 5; 0.291, 0.294, 0.295 at D = 6), so it is not truncation. G's
formation radius depends on δ too (1.172 against 1.098 at D = 5; 1.736 against 1.451 at D = 6): its near-critical state is not the static throat of radius a, and the stall member drifts (#255: at D = 5 the throat radius of G's stall member
exceeds a and grows with B; the drift has a longer time to act at smaller δ). The ratio of the settled r_AH to the formation r_AH is 0.620 and 0.624 at D = 5 for G (δ = 1e-8, 1e-6) and 0.691 and 0.719 at D = 6, against 0.544 and 0.604 for S at either δ (arithmetic on the printed first rows
and plateaus, finest h). *Conjecture, not tested here:* the settled mass follows the stall member's size at departure, as #253 and #267 test at D = 4; D = 5, 6 are the case where it would be visible, because the δ dependence is 50 times (D = 5) and 500 times (D = 6) larger than at D = 4.
G at D = 6, h = 0.00125 uses the p*ₕ bisected here; at the two other h the settled value is the same function of h (second order, RC6).

## The comparison with [ST17]

[ST17] (arXiv:1706.02070, read from the TeX source; hashes equal #41's `papers.sha256`) Table I, GR rows: E_f/a₀ = 3.14 (n = 4), 6.28 (n = 5), 9.87 (n = 6) for both pulses (c₁ = 0.25, 0.5), E evaluated at x⁺ = 5, from their (2.47)
E_n = (n−2)A_{n−2}/(2κ²Ω){k + 2e^f ϑ₊ϑ₋/(n−2)²}, Ω = 1/r, κ² = 1. **Conversion:** the Misner–Sharp mass here is M = r^{D−3}(1 − (∇r)²)/2 with 8πG = κ² = 1, so E = c_D M with c_D = (D − 2) V_{D−2}, where V_2 = π, V_3 = 4π/3,
V_4 = π²/2 are the volumes of the unit balls, which is what the printed A₂, A₃, A₅ are. The static throat (M = 1/2) then gives E = π, 2π, π² = 3.14, 6.28, 9.87 (CONV passes): [ST17]'s E_f equals E_throat in each D. The printed A₄ = π²/Ω⁴ is
twice V_4 and gives E_throat = 2π² = 19.74 at n = 6, against the table's 9.87 (A4P: ratio 1.9999); the table needs V_4 = π²/2. This is a disagreement inside the paper, not fixed here.

| D | quantity | M | E = c_D M | E/E_throat | [ST17] E_f | [ST17] E_f/E_throat |
|---|---|---|---|---|---|---|
| 5 | S, formation (first trapped cell) | 0.49205 | 6.183 | 0.984 | 6.28 | 0.9995 |
| 5 | S, settled | 0.14586 | 1.833 | 0.292 | 6.28 | 0.9995 |
| 5 | G, formation | 0.68728 | 8.637 | 1.375 | 6.28 | 0.9995 |
| 5 | G, settled (δ = 1e-8) | 0.26392 | 3.317 | 0.528 | 6.28 | 0.9995 |
| 6 | S, formation | 0.49003 | 9.673 | 0.980 | 9.87 | 1.0000 |
| 6 | S, settled | 0.10779 | 2.128 | 0.216 | 9.87 | 1.0000 |
| 6 | G, formation | 2.61513 | 51.62 | 5.230 | 9.87 | 1.0000 |
| 6 | G, settled (δ = 1e-8) | 0.86218 | 17.02 | 1.724 | 9.87 | 1.0000 |

[ST17]'s E_f/E_throat = 1.00 matches the formation value of S (0.984, 0.980; r_AH = 0.992 a, 0.993 a at the first trapped cell) to 2 %, not the settled one (0.292, 0.216, which is 3.4 and 4.6 times smaller), as at D = 4 (0.988 and 0.446). They read E at the
end of a grid x⁺ = 5, 2 to 3 a₀ after their x_H⁻; M_AH is still falling then (D = 5: 0.167 at x⁺ − x⁺_f = 2 and 0.149 at 2.5, against 0.14586 settled). Their pulse (c₁, c₂, c₃ = 0.25 or 0.5, 16, 0.7) is not the threshold family, so the comparison is of
the quantity read and when, not of data. Not compared: their x_H⁻ column.

## Is M_AH's convergence at the roundoff floor of #226?

No. Noise of 2⁻⁵² relative on every node of both rays (three seeds), h = the two finest per D, δ = 1e-8; ratio of the scatter max |M_noisy − M| to the h difference E between that h and its coarser neighbour (the declared reading: ≥ 1 reached,
0.1 to 1 approached, < 0.1 not reached; windows x⁺ − x⁺_f in [0, 1.5], [1.5, 3], [3, end]):

| case | h | scatter | E | ratio, three windows |
|---|---|---|---|---|
| D = 5, S | 0.00354 | 1.5e-5, 1.5e-5, 4.8e-6 | 1.1e-3, 1.1e-3, 1.2e-4 | 0.013, 0.013, 0.040 |
| D = 5, S | 0.00177 | 1.6e-5, 1.6e-5, 1.2e-6 | 2.8e-4, 2.8e-4, 3.0e-5 | 0.056, 0.057, 0.041 |
| D = 5, G | 0.00354 | 2.7e-5, 2.7e-5, 5.6e-6 | 2.6e-3, 2.6e-3, 3.9e-4 | 0.010, 0.010, 0.014 |
| D = 5, G | 0.00177 | 3.3e-5, 3.3e-5, 2.7e-6 | 6.2e-4, 6.3e-4, 9.5e-5 | 0.052, 0.052, 0.028 |
| D = 6, S | 0.0025 | 1.3e-4, 4.5e-5, 4.7e-6 | 1.6e-2, 5.8e-3, 1.6e-3 | 0.008, 0.008, 0.003 |
| D = 6, S | 0.00125 | 2.0e-4, 6.7e-5, 1.3e-6 | 3.9e-3, 1.5e-3, 4.0e-4 | 0.051, 0.044, 0.003 |
| D = 6, G | 0.0025 | 1.6e-3, 1.5e-3, 1.4e-4 | 7.1e-2, 6.6e-2, 1.0e-2 | 0.022, 0.023, 0.013 |
| D = 6, G | 0.00125 | 2.0e-3, 1.9e-3, 1.7e-4 | 1.6e-2, 1.5e-2, 2.4e-3 | 0.124, 0.125, 0.072 |

- Not reached everywhere except D = 6, G, h = 0.00125, first two windows, where it is approached (0.12). The second-order orders of RC6 are 1.92 to 2.13 at these h, so the floor does not show in the convergence yet.
- The scatter does not fall with h while E falls by 4 per halving, so the ratio rises 4 to 6 times per halving (D = 5 S: 0.013 to 0.056; D = 6 S: 0.008 to 0.051; D = 6 G: 0.022 to 0.124). If that continues, the ratio reaches 1 about 1.5 (D = 6 G) to 2.2
  halvings (D = 5, 6 S and D = 5 G, first window: 2.1 to 2.2) below the finest h: an extrapolation of two numbers, not tested (an extension h needs its own p*ₕ).
- The scatter sits where M falls fastest: its maximum is at x⁺ − x⁺_f = 1.2 to 1.5 (a scratch check on the stored runs, not in the output file), so it is a timing jitter of the formation, δM ≈ (dM/dx⁺) δx⁺ with δx⁺ of order 2e-3 at D = 6 G.
- What the scatter is: the difference between roundoff realizations of the run, not a calibrated response (#253): it does not grow from 2⁻⁵² to 64 × 2⁻⁵² (D = 5, S, h = 0.00177: 1.6e-5 and 1.2e-5; D = 6, G, h = 0.00125: 2.0e-3 and 2.9e-4). The estimator's control, RNR (the scatter rises at least 10 times at 2²⁰ × 2⁻⁵²), passes in all four
  cases (34 to 1420 times), unlike #253's at 64×. The zero-noise run reproduces the run of record bit for bit (RN0). A constraint diagnostic at the floor (#226) and M_AH's scatter are different quantities; this says only that the second-order differences of M_AH are 8 to 330 times the realization scatter at the finest two h.

## Disagreements and checks that did not pass

1. **RC8a** (settled at #41's limit 1.57e6) fails in all four cases (numbers above). Declared as expected. At the K = 28 limit all RC8 pass. The settled values at the two limits differ by 4e-6 or less at the end and by 1e-5 to 6.6e-5 in the plateau mean.
2. **Count.** README declared 72 checks; the analysis counts 69. PL1, PL2, PL3 are produced by #253's `planted_through_driver`, which records them in its own list: they print `[pass]` but are neither counted nor able to change the exit code. The declaration's
   count of 72 included them. They pass here. A fix (add `smg.results` to this driver's list) changes the analysis code only and is not made after the runs of record; it is recommended below.
3. **[ST17]'s printed A₄** (above): inside the paper.
4. **The comparison bar for G against S** adds a physical δ dependence to the error terms, as #253 did, and the verdict holds by a factor 4 (D = 5) and 2.5 (D = 6) against it; against the h terms alone it holds by factors 1500 and 370. The machine term 2e-5 is #41's at D = 4.
5. **p*ₕ at D = 5, G, h = 0.00707** has hi − lo = 2.9e-12 (#40's bisection ended on "static through X"); the bar for RC1 (1e-11) was declared for it; every δ is 1e3 times larger.
6. Interruptions (reboot at 19:59 and suspend at 21:23 on 2026-11-27) stopped the lanes; every interrupted run was kept (`*.run1_suspend`, `dev/mut/*.run1_reboot`) and rerun from the same queue and code (raw README, posts 01 to 03). The runs of record that finished are not rerun.

## Stale or conflicting text found (none edited)

- `paper/5-ellis-threshold/main.tex` line 554: "Masses were measured in D = 4 only (#246)" and the `\pending{\#246}` placeholder (Sec. 3.4).
- `settled_mass.md` Limits: "D = 4 and family S (plus SH02's data) only", and recommended task 4 ("D = 5, 6": done here).
- `critical_ghost.md` §1.3: "[ST17]: … their normalization of the mass is not reconciled here with SH's" is reconciled by the conversion above (E = (D − 2) V_{D−2} M); "Added by the coordinator (from #40): 0.494, 0.492 and 0.490 a^{D−3} at D = 4, 5, 6" agrees with the formation values here (0.492, 0.490).

## Proposed text (for the coordinator; none applied)

**`settled_mass.md`**, Limits, first bullet: after "…family S (plus SH02's data) only." add "D = 5, 6 are measured in #246 (`settled_mass_d56.md`): S settles at 0.14587 a² and 0.10792 a³, G at 0.2640 a² and 0.8629 a³ (δ = 1e-8)." and mark
recommended task 4 "done in #246".

**`critical_ghost.md` §1.3**, after the "Update from #253" item, add:

> *Update from #246* (`wormhole_evolution/settled_mass_d56.md`; *numerical evidence*). The same chain at D = 5 and 6 (M = r^{n−1}/2 on the apparent horizon; second order at three h; a planted Tangherlini recovered to 3e-6 to 6e-6; to e^{−f} = 4.0e8 on the horizon)
> gives, for the symmetric family S, a settled mass of 0.14587 a² (r_AH = 0.5401 a) at D = 5 and 0.10792 a³ (r_AH = 0.5996 a) at D = 6, E/E_throat = 0.292 and 0.216 against 0.446 at D = 4. [ST17]'s E_f is the static throat's energy in each D (E = (D − 2)V_{D−2} M with the unit-ball
> volumes; their printed A₄ = π² is twice the 4-ball's π²/2), and the horizon's mass at its first trapped cell is 0.984 and 0.980 of it: a formation value, as at D = 4. The one-sided family G settles higher (0.2640 a² and 0.8629 a³ at δ = 1e-8, 0.118 a² and 0.754 a³ above S's) and depends strongly on δ
> (0.2350 and 0.5679 at δ = 1e-6, the same at all three h), with a formation horizon larger than the throat (r_AH = 1.17 a and 1.74 a): the settled mass is not the same for every family at any D, and the dependence grows with D. That it follows the stall member's size at departure is a *conjecture*.

and in §4 (the bullet "The threshold mass is 0.30 a [SH02] or 0.22 a [GGS09b]") append: "*Update from #246:* S settles at 0.146 a² (D = 5) and 0.108 a³ (D = 6); [ST17]'s E_f is the formation value at D = 5, 6 too."

**Paper V**, Sec. 3.4, replacing `\pending{\#246}{…}` and the sentence "Masses were measured in $D=4$ only (\#246)" (line 554):

```latex
In $D=5$ and $6$ the same chain, with $M=r^{n-1}/2$ on the horizon, gives for the symmetric family a settled mass of
$0.1459\,a^{2}$ ($r_\AH=0.5401a$) and $0.1079\,a^{3}$ ($r_\AH=0.5996a$) at $\delta=10^{-8}$, that is $E/E_{\rm throat}=2M/a^{n-1}=0.292$
and $0.216$ against $0.446$ in $D=4$; the series are second order in $h$ at three resolutions to $e^{-f}=4\times10^{8}$ on the
horizon, and a planted Tangherlini is recovered to $3\times10^{-6}$ to $6\times10^{-6}$ in $M$ at second order
(\texttt{notes/wormhole\_evolution/settled\_mass\_d56.md}). Shinkai and Torii's final energy is the static throat's own, $E_f/E_{\rm throat}=1.00$
in each $D$ (their $E=(D-2)V_{D-2}M$ with the unit-ball volumes, which the table needs; the printed $A_4$ is twice that), and the mass on the
horizon at its first trapped cell is $0.984$ and $0.980$ of it: their value is the formation value in $D=5$ and $6$ as in $D=4$. The one-sided
family G settles at $0.264\,a^{2}$ and $0.863\,a^{3}$, $0.118\,a^{2}$ and $0.754\,a^{3}$ above the symmetric family, and its mass depends strongly
on $\delta$ ($0.235$ and $0.568$ at $\delta=10^{-6}$).
```

and in the "Reading" paragraph: "…and the dependence on the distance from threshold grows with the dimension (G: 0.0125a at $D=4$, 0.118a² at $D=5$, 0.754a³ at $D=6$ above S)." `paper_numbers.py` would check 0.14587, 0.10792, 0.26395, 0.86286, 0.984, 0.980 from `settled_mass_d56_output.txt`.

## Limits

- Three h per D (G at D = 6: two of #40's and one bisected here), two δ (1e-6, 1e-8), B = 0.1 for G, the throat-pulse families of #40 only; no event horizon, no I⁺. The settled value is to e^{−f} = 4.0e8 (the flux over the last 0.5 is 5e-7 to 1.1e-6); #41's flux that returned once at D = 4 was not seen here, and no extrapolation is made.
- The δ dependence of G rests on two δ (1e-6, 1e-8): no law and no limit as δ → 0 is claimed (#267 at D = 4).
- The planted cases are two (λ, kp) per D, not a scan; the planted data are exact, so they do not test the floor.
- The Richardson value assumes second order, which RC6 supports (1.92 to 2.13); its difference from the finest-h plateau is up to 7e-4 (D = 6, G).
- The machine term 2e-5 of the G–S bar is #41's at D = 4. All runs are on one machine (dh610); p*ₕ of #40 are from the same machine.
- Mutations run on D = 5 (S, or G for the normal-field ones) and the planted D = 5 case; PL1 to PL3 are #41's and are not exercised by them; `col_shift` is silent for RC6 at D = 5 and `normal_drop`/`normal_sign` for RC4 (as in #253).

## Recommended tasks

1. *G's settled mass against δ and the stall member at D = 5, 6 (#255, #267).* What: M_settled(δ) for δ = 1e-12 … 1e-4 at the two finest h, with the stall member's radius at departure (`throat_drift.py`'s measure) at each δ. Why: the dependence is 50 to 500 times larger than at D = 4.
   Acceptance: M_settled is one function of the member's radius at departure for six δ at both D, or the failure is reported. Depends on: #255, #267, this note.
2. *[ST17]'s own data through the chain at n = 5, 6.* What: their pulse π₊ = c₁ exp(−c₂(x⁺ − c₃)²), c₂ = 16, c₃ = 0.7, c₁ = 0.25, 0.5, on the static Ellis wormhole of radius a₀, run through this chain; E at x⁺ = 5 and at formation, and x_H⁻. Why: their Table I columns E_i, x_H⁻.
   Acceptance: their x_H⁻ (2.26, 1.66 at n = 5; 1.92, 1.46 at n = 6) within 5 %, their E_f/E_throat = 1.00 against ours at their reading; a planted Tangherlini recovered. Depends on: this note.
3. *The roundoff floor at the extension h.* What: p*ₕ at h/2 and h/4 of the finest h of S at D = 5, 6 and the noise runs there. Why: the ratio is 0.05 to 0.12 at the finest h and rises 4 to 6 per halving. Acceptance: the ratio ≥ 1 measured at some h and the order of RC6 reported there. Depends on: #226, this note.
4. *Count PL1 to PL3 in the exit code.* What: add `smg.results` to this driver's results (or record them with `check`) and rerun the analysis of record. Why: the declared 72 and the counted 69 differ by exactly these three. Acceptance: 72 counted, PL1 to PL3 fire under `mass_exponent`. Depends on: this note.
5. *Settled mass beyond K = 28 and the event horizon (#247).* Same as #41's tasks 2 and 5 at D = 5, 6. Depends on: this note, #247.

## Files

`settled_mass_d56.py` (commands `pstar`, `job`, `planted`, `analyze`, `tang_test`; `--mutate=NAME`), `settled_mass_d56_output.txt`, this note; raw `negative-mass-raw/I246-dh610` (`README.md`: the declaration and its changes; `runs/`: p*ₕ, the bisection, the 24 runs, the 8 prefix runs, the 48 noise runs, the
analysis log; `mut/`: the 26 mutated analyses of record; `dev/`: development runs, queues and the 26 development mutations; `posts/`: the issue posts; `scripts/`: lane scripts, `bisect_pstar.py`, `make_figure.py`, `data_graphs.py`; `figures/settled_mass_d56.svg` and its table).

## Coordinator's review

Merged at 1fca2ef. Review folder `negative-mass-raw/I246-mac-review` at 6f436e4c, criteria committed (9a59b906) before the first run of
`ind246.py`. The changed files are those the brief names.

- **Rerun on the Mac** from 5b3517b: `analyze` gives the same 72 verdict lines as the output of record (exit 1, RC8a). Every line
  read from the stored runs is identical. The planted Tangherlini runs, integrated afresh on the analysing machine, differ in
  the third or fourth printed digit between dh610 (x86-64) and the Mac (aarch64), each within its bar.
- **ind246**, 4 of 4. It shares only the series with this note's code. It recomputes the finest-h plateaus (to 4.2e-7), the
  Richardson values (to 3.6e-7; orders 1.99 to 2.13 in all eight cases), E/E_throat (0.2917, 0.2158), and the bit-identity of
  the K = 20 and 24 prefixes (8 of 8). The mutation `order1` fires. `last-half` and `ball` do not: the plateau is flat to 9e-7,
  and the unit-ball volume cancels in E/E_throat, so that check tests only M/(½a^{n−1}). The absolute throat energies π, 2π, π²
  were recomputed separately.
- The interruptions (a reboot at 03:59 UTC and a login-screen suspend at 05:23 UTC on dh610) are documented in the raw folder's
  posts. Every partial log is kept, and all 80 runs of record exit 0. dh610's sleep targets are now masked.
- Applied: `settled_mass.md`, `critical_ghost.md` §1.3 and §4, and Paper V's #246 paragraph and limits bullet
  (`paper_numbers.py` 251 of 251). Paper V quotes the Richardson values, so G − S is 0.755 a³ there against the finest-h 0.754.
